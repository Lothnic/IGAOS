#include "pdhg.hpp"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#define CK(x)                                                                   \
    do {                                                                        \
        auto e_ = (x);                                                          \
        if (e_ != cudaSuccess && e_ != CUBLAS_STATUS_SUCCESS &&                 \
            e_ != CUSPARSE_STATUS_SUCCESS) {                                    \
            std::fprintf(stderr, "[pdhg] error %d at %d\n", (int)e_, __LINE__);\
            sol.message = "CUDA error in PDHG engine";                          \
            return sol;                                                         \
        }                                                                       \
    } while (0)

namespace igaos::pdhg {
namespace {

// Fused PDHG iteration kernels: one launch per update instead of a chain of
// elementwise launches + blases (kernel-launch latency dominated small
// instances: ~8 launches/iter -> 3).

// xnew = P_[cl,cu](x - tau*(c + A'y));  xbar = 2*xnew - x
__global__ void kern_primal(int n, double tau, const double* x,
                            const double* c, const double* ata,
                            const double* lo, const double* hi, double* xnew,
                            double* xbar) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        double g = c[i] + ata[i];
        double xn = fmin(fmax(x[i] - tau * g, lo[i]), hi[i]);
        xnew[i] = xn;
        xbar[i] = 2.0 * xn - x[i];
    }
}

// t = y + sigma*Axbar;  ynew = t - sigma*P_[rmin,rmax](t/sigma)
__global__ void kern_dual(int m, double s, const double* y, const double* axb,
                          const double* rmin, const double* rmax,
                          double* ynew) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) {
        double t = y[i] + s * axb[i];
        ynew[i] = t - s * fmin(fmax(t / s, rmin[i]), rmax[i]);
    }
}

// Halpern averaging of both x and y in one launch.
__global__ void kern_mix2(int n, int m, double eta, const double* x,
                          double* xavg, const double* y, double* yavg) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) xavg[i] += eta * (x[i] - xavg[i]);
    else if (i < n + m) {
        int j = i - n;
        yavg[j] += eta * (y[j] - yavg[j]);
    }
}

struct Spmv {
    cusparseHandle_t h;
    cusparseSpMatDescr_t mat;
    cusparseDnVecDescr_t vin, vout;
    size_t buf = 0;
    void* dbuf = nullptr;

    void mv(const double* x, double* out) {
        cusparseDnVecSetValues(vin, (void*)x);
        cusparseDnVecSetValues(vout, (void*)out);
        double alpha = 1.0, beta = 0.0;
        cusparseSpMV(h, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat, vin,
                     &beta, vout, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, dbuf);
    }
};

static constexpr double INF_NAN = std::numeric_limits<double>::infinity();

struct Metrics {
    double pinf, dinf, gap, op, od;
    bool finite_res() const {
        return std::isfinite(pinf) && std::isfinite(dinf);
    }
    double err() const {
        double g = std::isfinite(gap) ? gap : 0.0;
        return std::max(pinf, std::max(dinf, g));
    }
};

Metrics evaluate(int m, int n, const std::vector<double>& rmin,
                 const std::vector<double>& rmax,
                 const std::vector<double>& cl, const std::vector<double>& cu,
                 const std::vector<double>& c, const std::vector<double>& x,
                 const std::vector<double>& y, const std::vector<double>& Ax,
                 const std::vector<double>& d) {
    const double INF = std::numeric_limits<double>::infinity();
    Metrics mt{};
    double pv = 0.0, pa = 0.0;
    for (int i = 0; i < m; ++i) {
        pv = std::max(pv, std::max(rmin[i] - Ax[i], Ax[i] - rmax[i]));
        pa = std::max(pa,
                      std::fabs(std::min(std::max(Ax[i], rmin[i]), rmax[i])));
    }
    mt.pinf = std::max(pv, 0.0) / (1.0 + pa);
    double dv = 0.0;
    for (int j = 0; j < n; ++j) {
        double bnd =
            std::fabs(std::isfinite(cu[j]) ? cu[j] : cl[j]);
        double tol = 1e-9 * (1.0 + bnd);
        double v;
        if (x[j] <= cl[j] + tol && d[j] > 0) v = 0.0;
        else if (x[j] >= cu[j] - tol && d[j] < 0) v = 0.0;
        else v = std::fabs(d[j]);
        dv = std::max(dv, v);
    }
    mt.dinf = dv / (1.0 + [&]{ double mx=0; for(double cj:c) mx=std::max(mx,std::fabs(cj)); return mx; }());
    mt.op = 0.0;
    for (int j = 0; j < n; ++j) mt.op += c[j] * x[j];
    // Lagrangian dual bound at y (code convention: d = c + A'y, prox maps y
    // into the row-sign-feasible cone):
    //   D(y) = min_{x in [cl,cu]} d'x  -  max_{z in [rmin,rmax]} y'z
    //       = sum_j [d_j>=0 ? d_j*cl_j : d_j*cu_j]
    //         + sum_i [y_i>0 ? -y_i*rmax_i : -y_i*rmin_i]
    // D(y) <= p* for every primal-feasible point, so (op - D(y))/(1+|op|)
    // is a valid relative gap. A wrong-sign d_j/y_i at an infinite bound
    // makes D(y) = -inf (gap not certifiable -> no termination on gap).
    double od = 0.0;
    bool bounded = true;
    double cmax = 0.0, ymax = 0.0;
    for (double cj : c) cmax = std::max(cmax, std::fabs(cj));
    for (double yi : y) ymax = std::max(ymax, std::fabs(yi));
    const double dz = 1e-8 * (1.0 + cmax);  // noise slack at infinite bounds
    const double yz = 1e-8 * (1.0 + ymax);
    for (int j = 0; j < n; ++j) {
        if (d[j] >= 0) {
            if (std::isfinite(cl[j])) od += d[j] * cl[j];
            else if (d[j] > dz) { bounded = false; break; }
        } else {
            if (std::isfinite(cu[j])) od += d[j] * cu[j];
            else if (-d[j] > dz) { bounded = false; break; }
        }
    }
    if (bounded) {
        for (int i = 0; i < m; ++i) {
            if (y[i] > yz) {
                if (std::isfinite(rmax[i])) od -= y[i] * rmax[i];
                else { bounded = false; break; }
            } else if (y[i] < -yz) {
                if (std::isfinite(rmin[i])) od -= y[i] * rmin[i];
                else { bounded = false; break; }
            }
        }
    }
    mt.od = bounded ? od : -INF;
    mt.gap = bounded ? std::max(0.0, mt.op - od) / (1.0 + std::fabs(mt.op))
                     : INF;
    return mt;
}

}  // namespace

Solution solve(const io::Model& lp, const Options& opt) {
    Solution sol;
    cublasHandle_t cb;
    cusparseHandle_t cs;
    if (cublasCreate(&cb) != CUBLAS_STATUS_SUCCESS ||
        cusparseCreate(&cs) != CUSPARSE_STATUS_SUCCESS) {
        sol.message = "failed to create CUDA library handles";
        return sol;
    }

    int m = lp.m, n = lp.n, nnz = lp.nnz();

    // Ruiz equilibration (rows dr, cols dc) + cost scaling gamma, PDLP-style.
    // Scaled problem: A~ = dr^-1 A dc, r~ = dr^-1 r, c~ = dc c / gamma,
    // c~l = dc^-1 cl, c~u = dc^-1 cu. Iterates map back x = dc x~, y = y~/dr.
    std::vector<double> dr(m, 1.0), dc(n, 1.0);
    for (int sweep = 0; sweep < 10; ++sweep) {
        for (int i = 0; i < m; ++i) {
            double mx = 0.0;
            for (int p = lp.ap[i]; p < lp.ap[i + 1]; ++p)
                mx = std::max(mx,
                              std::fabs(lp.ax[p]) * dc[lp.ai[p]] / dr[i]);
            if (!(mx > 0) || !std::isfinite(mx)) continue;
            dr[i] *= std::sqrt(mx);
        }
        for (int j = 0; j < n; ++j) {
            double mx = 0.0;
            for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                mx = std::max(mx,
                              std::fabs(lp.acx[p]) / dr[lp.ci[p]] * dc[j]);
            if (!(mx > 0) || !std::isfinite(mx)) continue;
            dc[j] /= std::sqrt(mx);  // values grow with dc: divide, not mult
        }
    }
    for (double& v : dr) if (!(v > 0) || !std::isfinite(v)) v = 1.0;
    for (double& v : dc) if (!(v > 0) || !std::isfinite(v)) v = 1.0;
    // Cost scale gamma: shrinks the dual travel distance y~* = dr*y*/gamma.
    // Mapping back (derived from d~ = Dc(c/gamma + A'y)):
    //   x = dc*x~,  y = gamma*y~/dr,  Ax = dr*(A~x~),  A'y = gamma*(A~'y~)/dc
    // Cost scale gamma stays 1: gamma=||dc*c||inf shrank the dual signal and
    // stalled convergence (kb2 -1735 vs -1749.9, afiro dinf plateau 2.8e-2).
    const double gamma = 1.0;
    std::vector<double> sax(nnz), sacx(nnz), sc(n), scl(n), scu(n), srmin(m),
        srmax(m);
    for (int i = 0; i < m; ++i) {
        srmin[i] = lp.rmin[i] / dr[i];
        srmax[i] = lp.rmax[i] / dr[i];
        for (int p = lp.ap[i]; p < lp.ap[i + 1]; ++p)
            sax[p] = lp.ax[p] * dc[lp.ai[p]] / dr[i];
    }
    for (int j = 0; j < n; ++j) {
        sc[j] = lp.c[j] * dc[j] / gamma;
        scl[j] = lp.cl[j] / dc[j];
        scu[j] = lp.cu[j] / dc[j];
        for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
            sacx[p] = lp.acx[p] * dc[j] / dr[lp.ci[p]];
    }

    int *d_ap, *d_ai, *d_cp, *d_ci;
    double *d_ax, *d_acx, *d_c, *d_cl, *d_cu, *d_rmin, *d_rmax;
    CK(cudaMalloc(&d_ap, (m + 1) * sizeof(int)));
    CK(cudaMalloc(&d_ai, nnz * sizeof(int)));
    CK(cudaMalloc(&d_ax, nnz * sizeof(double)));
    CK(cudaMalloc(&d_cp, (n + 1) * sizeof(int)));
    CK(cudaMalloc(&d_ci, nnz * sizeof(int)));
    CK(cudaMalloc(&d_acx, nnz * sizeof(double)));
    CK(cudaMalloc(&d_c, n * sizeof(double)));
    CK(cudaMalloc(&d_cl, n * sizeof(double)));
    CK(cudaMalloc(&d_cu, n * sizeof(double)));
    CK(cudaMalloc(&d_rmin, m * sizeof(double)));
    CK(cudaMalloc(&d_rmax, m * sizeof(double)));
    CK(cudaMemcpy(d_ap, lp.ap.data(), (m + 1) * sizeof(int),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_ai, lp.ai.data(), nnz * sizeof(int),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_ax, sax.data(), nnz * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cp, lp.cp.data(), (n + 1) * sizeof(int),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_ci, lp.ci.data(), nnz * sizeof(int),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_acx, sacx.data(), nnz * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_c, sc.data(), n * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cl, scl.data(), n * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cu, scu.data(), n * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_rmin, srmin.data(), m * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_rmax, srmax.data(), m * sizeof(double),
                  cudaMemcpyHostToDevice));

    double *d_x, *d_xnew, *d_xbar, *d_xavg, *d_y, *d_ynew, *d_yavg, *d_ata,
        *d_tm, *d_tmpn;
    CK(cudaMalloc(&d_x, n * sizeof(double)));
    CK(cudaMalloc(&d_xnew, n * sizeof(double)));
    CK(cudaMalloc(&d_xbar, n * sizeof(double)));
    CK(cudaMalloc(&d_xavg, n * sizeof(double)));
    CK(cudaMalloc(&d_y, m * sizeof(double)));
    CK(cudaMalloc(&d_ynew, m * sizeof(double)));
    CK(cudaMalloc(&d_yavg, m * sizeof(double)));
    CK(cudaMalloc(&d_ata, n * sizeof(double)));
    CK(cudaMalloc(&d_tm, m * sizeof(double)));
    CK(cudaMalloc(&d_tmpn, n * sizeof(double)));
    CK(cudaMemset(d_x, 0, n * sizeof(double)));
    CK(cudaMemset(d_y, 0, m * sizeof(double)));
    CK(cudaMemset(d_xavg, 0, n * sizeof(double)));
    CK(cudaMemset(d_yavg, 0, m * sizeof(double)));

    Spmv mvA{}, mvAT{};
    mvA.h = cs; mvAT.h = cs;
    CK(cusparseCreateCsr(&mvA.mat, m, n, nnz, d_ap, d_ai, d_ax,
                         CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                         CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
    CK(cusparseCreateCsr(&mvAT.mat, n, m, nnz, d_cp, d_ci, d_acx,
                         CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                         CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
    CK(cusparseCreateDnVec(&mvA.vin, n, d_tmpn, CUDA_R_64F));
    CK(cusparseCreateDnVec(&mvA.vout, m, d_tm, CUDA_R_64F));
    CK(cusparseCreateDnVec(&mvAT.vin, m, d_y, CUDA_R_64F));
    CK(cusparseCreateDnVec(&mvAT.vout, n, d_ata, CUDA_R_64F));
    size_t bA = 0, bAT = 0;
    double one = 1.0, zero = 0.0;
    CK(cusparseSpMV_bufferSize(cs, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                               mvA.mat, mvA.vin, &zero, mvA.vout, CUDA_R_64F,
                               CUSPARSE_SPMV_ALG_DEFAULT, &bA));
    CK(cusparseSpMV_bufferSize(cs, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                               mvAT.mat, mvAT.vin, &zero, mvAT.vout,
                               CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &bAT));
    CK(cudaMalloc(&mvA.dbuf, std::max(bA, (size_t)1)));
    CK(cudaMalloc(&mvAT.dbuf, std::max(bAT, (size_t)1)));

    double lam2 = 1.0;
    {
        std::vector<double> hv(n, 1.0 / std::sqrt((double)n));
        double *dv, *dw;
        CK(cudaMalloc(&dv, n * sizeof(double)));
        CK(cudaMalloc(&dw, n * sizeof(double)));
        CK(cudaMemcpy(dv, hv.data(), n * sizeof(double),
                      cudaMemcpyHostToDevice));
        for (int it = 0; it < 30; ++it) {
            mvA.mv(dv, d_tm);
            mvAT.mv(d_tm, dw);
            double nw = 0.0;
            CK(cublasDnrm2(cb, n, dw, 1, &nw));
            if (nw < 1e-300) break;
            double inv = 1.0 / nw;
            CK(cublasDscal(cb, n, &inv, dw, 1));
            CK(cudaMemcpy(dv, dw, n * sizeof(double),
                          cudaMemcpyDeviceToDevice));
        }
        mvA.mv(dv, d_tm);
        double na = 0.0;
        CK(cublasDnrm2(cb, m, d_tm, 1, &na));
        lam2 = na * na;
        if (!(lam2 > 0) || !std::isfinite(lam2)) lam2 = 1.0;
        cudaFree(dv);
        cudaFree(dw);
    }
    double ts_max = 1.2 / lam2, ts = 1.2 * ts_max;
    double tau = std::sqrt(ts), sigma = tau;

    const int CHECK = 25;
    // Under an explicit wall-clock budget the time limit governs; the 500k
    // default iteration cap would otherwise cut runs short with time left
    // (share2b was still improving at k=500k inside its 16s budget).
    long MAXIT = (opt.time_limit_s > 0 && opt.time_limit_s < 1e6)
                     ? 20000000L
                     : opt.max_iterations;
    const double TOL = opt.tolerance;
    long hn = 0;
    long since_restart = 0;
    double epoch_r0 = -1.0;
    long iters_run = 0;
    auto t0 = std::chrono::steady_clock::now();

    std::vector<double> hx(n), hy(m), hax(m), hd(n);
    // metrics_of evaluates KKT in ORIGINAL problem units: copy the scaled
    // iterates, unscale (x = dc*x~, y = y~/dr), and post-scale the SpMVs
    // (Ax = dr*(A~x~), A'y = (A~'y~)/dc).
    auto metrics_of = [&](const double* dx, const double* dy, Metrics& mt) {
        cudaError_t e1 = cudaMemcpy(hx.data(), dx, n * sizeof(double),
                                    cudaMemcpyDeviceToHost);
        cudaError_t e2 = cudaMemcpy(hy.data(), dy, m * sizeof(double),
                                    cudaMemcpyDeviceToHost);
        if (e1 != cudaSuccess || e2 != cudaSuccess) {
            std::fprintf(stderr, "[pdhg] sync failed %d %d\n", (int)e1,
                         (int)e2);
            mt.pinf = mt.dinf = INF_NAN;
            mt.gap = INF_NAN;
            mt.op = 0;
            mt.od = 0;
            return;
        }
        for (int j = 0; j < n; ++j) hx[j] *= dc[j];
        for (int i = 0; i < m; ++i) hy[i] *= gamma / dr[i];
        mvA.mv(dx, d_tm);
        cudaError_t e3 = cudaMemcpy(hax.data(), d_tm, m * sizeof(double),
                                    cudaMemcpyDeviceToHost);
        for (int i = 0; i < m; ++i) hax[i] *= dr[i];
        mvAT.mv(dy, d_ata);
        cudaError_t e4 = cudaMemcpy(hd.data(), d_ata, n * sizeof(double),
                                    cudaMemcpyDeviceToHost);
        for (int j = 0; j < n; ++j) hd[j] = lp.c[j] + gamma * hd[j] / dc[j];
        if (e3 != cudaSuccess || e4 != cudaSuccess) {
            std::fprintf(stderr, "[pdhg] sync2 failed %d %d\n", (int)e3,
                         (int)e4);
            mt.pinf = mt.dinf = INF_NAN;
            mt.gap = INF_NAN;
            mt.op = 0;
            mt.od = 0;
            return;
        }
        mt = evaluate(m, n, lp.rmin, lp.rmax, lp.cl, lp.cu, lp.c, hx, hy, hax,
                      hd);
    };

    Metrics last{};
    bool last_is_avg = false;
    long k = 1;
    for (; k <= MAXIT; ++k) {
        iters_run = k;
        mvAT.mv(d_y, d_ata);
        kern_primal<<<(n + 255) / 256, 256>>>(n, tau, d_x, d_c, d_ata, d_cl,
                                              d_cu, d_xnew, d_xbar);
        mvA.mv(d_xbar, d_tm);
        kern_dual<<<(m + 255) / 256, 256>>>(m, sigma, d_y, d_tm, d_rmin,
                                            d_rmax, d_ynew);
        std::swap(d_x, d_xnew);
        std::swap(d_y, d_ynew);
        double weta = 2.0 / ((double)hn + 2.0);
        ++hn;
        kern_mix2<<<(n + m + 255) / 256, 256>>>(n, m, weta, d_x, d_xavg, d_y,
                                                d_yavg);

        if (k % CHECK == 0) {
            Metrics mc, ma;
            metrics_of(d_x, d_y, mc);
            metrics_of(d_xavg, d_yavg, ma);
            double ec = mc.finite_res() ? mc.err() : INF_NAN;
            double ea = ma.finite_res() ? ma.err() : INF_NAN;
            last = (ec <= ea) ? mc : ma;
            last_is_avg = (ea < ec);
            if (opt.verbosity > 0 && k % (CHECK * 40) == 0)
                std::fprintf(stderr,
                             "[pdhg] k=%ld pinf=%.3e dinf=%.3e gap=%.3e\n",
                             k, last.pinf, last.dinf, last.gap);
            // Honest termination: residuals AND a certified duality gap,
            // on either the current or the averaged iterate. Residual-only
            // acceptance returns wrong objectives (sc205, adlittle):
            // dinf/pinf do not bound the objective error.
            auto passes = [&](const Metrics& t) {
                return t.pinf <= TOL && t.dinf <= TOL &&
                       std::isfinite(t.gap) && t.gap <= TOL;
            };
            if (passes(mc)) { last = mc; last_is_avg = false; break; }
            if (passes(ma)) { last = ma; last_is_avg = true; break; }

            // Adaptive step ratio (PDLP-style): every 500 iters (from k=1000)
            // nudge tau/sigma toward the measured pinf/dinf imbalance, product
            // tau*sigma pinned, move clamped to a factor of 1.4 and the ratio
            // banded to [1e-2, 1e2]. Without the band a primal-feasible
            // trapped iterate (pinf~0) multiplies tau by 0.7 every 500 iters
            // until it underflows and the primal freezes (kb2 denormal trap).
            if (k >= 1000 && k % 500 == 0 &&
                std::isfinite(last.pinf) && last.pinf > 0 &&
                std::isfinite(last.dinf) && last.dinf > 0) {
                double tgt = std::clamp(last.pinf / last.dinf, 0.1, 10.0);
                double cur = tau / sigma;
                double nxt = cur * std::clamp(tgt / cur, 0.7, 1.4);
                nxt = std::clamp(nxt, 0.1, 10.0);
                double prod = tau * sigma;
                tau = std::sqrt(prod * nxt);
                sigma = prod / tau;
            }
            double r_now = std::isfinite(ma.err()) ? ma.err()
                                                   : INF_NAN;
            if (epoch_r0 == INF_NAN) epoch_r0 = r_now;
            bool kkt_restart =
                std::isfinite(r_now) && r_now <= 0.3 * epoch_r0;
            if (kkt_restart || ++since_restart >= 2000) {
                hn = 0;
                CK(cudaMemcpy(d_xavg, d_x, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_yavg, d_y, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                since_restart = 0;
                if (kkt_restart && opt.verbosity > 1)
                    std::fprintf(stderr,
                                 "[pdhg] KKT restart @it=%ld r=%.3e\n",
                                 iters_run, r_now);
                epoch_r0 = INF_NAN;
            }

            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration<double>(now - t0).count() >
                opt.time_limit_s) {
                sol.status = Status::TimeLimit;
                sol.message = "wall-clock budget exhausted";
                break;
            }
        }
    }
    double secs = std::chrono::duration<double>(
                      std::chrono::steady_clock::now() - t0)
                      .count();

    auto certified = [&](const Metrics& t) {
        return t.pinf <= TOL && t.dinf <= TOL && std::isfinite(t.gap) &&
               t.gap <= TOL;
    };

    if (sol.status != Status::TimeLimit) {
        if (certified(last)) {
            sol.status = Status::NearOptimal;
            sol.message = "converged: residuals + duality gap @tol";
        } else {
            sol.status = Status::IterationLimit;
            sol.message = last.pinf <= TOL && last.dinf <= TOL
                              ? "iteration cap: residuals @tol but gap not "
                                "certified (dual bound -inf or > tol)"
                                : "iteration cap reached before tolerance";
        }
    }

    const double* fx = last_is_avg ? d_xavg : d_x;
    const double* fy = last_is_avg ? d_yavg : d_y;
    sol.objective = last.op + lp.obj_const;
    sol.pinf = last.pinf;
    sol.dinf = last.dinf;
    sol.rel_gap = last.gap;
    sol.iterations = iters_run;
    sol.solve_time_ms = secs * 1000.0;
    sol.x.resize(n);
    CK(cudaMemcpy(hx.data(), fx, n * sizeof(double),
                  cudaMemcpyDeviceToHost));
    for (int j = 0; j < n; ++j) sol.x[j] = hx[j] * dc[j];
    sol.y.resize(m);
    CK(cudaMemcpy(hy.data(), fy, m * sizeof(double),
                  cudaMemcpyDeviceToHost));
    for (int i = 0; i < m; ++i) sol.y[i] = gamma * hy[i] / dr[i];

    cudaFree(d_x); cudaFree(d_xnew); cudaFree(d_xbar); cudaFree(d_xavg);
    cudaFree(d_y); cudaFree(d_ynew); cudaFree(d_yavg);
    cudaFree(d_ata); cudaFree(d_tm); cudaFree(d_tmpn);
    cudaFree(d_ap); cudaFree(d_ai); cudaFree(d_ax);
    cudaFree(d_cp); cudaFree(d_ci); cudaFree(d_acx);
    cudaFree(d_c); cudaFree(d_cl); cudaFree(d_cu);
    cudaFree(d_rmin); cudaFree(d_rmax);
    cublasDestroy(cb);
    cusparseDestroy(cs);
    return sol;
}

}  // namespace igaos::pdhg
