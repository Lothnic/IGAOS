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

__global__ void kern_clamp(int n, const double* v, const double* lo,
                           const double* hi, double* out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = fmin(fmax(v[i], lo[i]), hi[i]);
}

__global__ void kern_two_minus(int n, const double* xnew, const double* xold,
                               double* xbar) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) xbar[i] = 2.0 * xnew[i] - xold[i];
}

__global__ void kern_add(int n, const double* a, const double* b, double* z) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) z[i] = a[i] + b[i];
}

__global__ void kern_dual(int m, const double* t, const double* rmin,
                          const double* rmax, double s, double* y) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) y[i] = t[i] - s * fmin(fmax(t[i] / s, rmin[i]), rmax[i]);
}

__global__ void kern_mix(int n, double eta, const double* cur, double* avg) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) avg[i] += eta * (cur[i] - avg[i]);
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
    double od = 0.0;
    bool bounded = true;
    double cmax = 0.0;
    for (double cj : c) cmax = std::max(cmax, std::fabs(cj));
    const double dz = 1e-10 * (1.0 + cmax);
    for (int j = 0; j < n; ++j) {
        if (d[j] >= 0) {
            if (!std::isfinite(cl[j])) {
                if (std::fabs(d[j]) > dz) { bounded = false; break; }
                continue;
            }
            od += d[j] * cl[j];
        } else {
            if (!std::isfinite(cu[j])) {
                if (std::fabs(d[j]) > dz) { bounded = false; break; }
                continue;
            }
            od += d[j] * cu[j];
        }
    }
    mt.od = bounded ? od : -INF;
    mt.gap = bounded ? std::fabs(mt.op - od) / (1.0 + std::fabs(mt.op)) : INF;
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
    CK(cudaMemcpy(d_ax, lp.ax.data(), nnz * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cp, lp.cp.data(), (n + 1) * sizeof(int),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_ci, lp.ci.data(), nnz * sizeof(int),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_acx, lp.acx.data(), nnz * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_c, lp.c.data(), n * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cl, lp.cl.data(), n * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cu, lp.cu.data(), n * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_rmin, lp.rmin.data(), m * sizeof(double),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_rmax, lp.rmax.data(), m * sizeof(double),
                  cudaMemcpyHostToDevice));

    double *d_x, *d_xnew, *d_xbar, *d_xavg, *d_y, *d_ynew, *d_t, *d_yavg,
        *d_ata, *d_grad, *d_tm, *d_tmpn;
    CK(cudaMalloc(&d_x, n * sizeof(double)));
    CK(cudaMalloc(&d_xnew, n * sizeof(double)));
    CK(cudaMalloc(&d_xbar, n * sizeof(double)));
    CK(cudaMalloc(&d_xavg, n * sizeof(double)));
    CK(cudaMalloc(&d_y, m * sizeof(double)));
    CK(cudaMalloc(&d_ynew, m * sizeof(double)));
    CK(cudaMalloc(&d_t, m * sizeof(double)));
    CK(cudaMalloc(&d_yavg, m * sizeof(double)));
    CK(cudaMalloc(&d_ata, n * sizeof(double)));
    CK(cudaMalloc(&d_grad, n * sizeof(double)));
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
    double ts_max = 0.9 / lam2, ts = 0.5 * ts_max;
    double tau = std::sqrt(ts), sigma = tau;

    const int CHECK = 25;
    long MAXIT = opt.max_iterations;
    const double TOL = opt.tolerance;
    double wsum = 0.0;
    long iters_run = 0;
    auto t0 = std::chrono::steady_clock::now();

    std::vector<double> hx(n), hy(m), hax(m), hd(n);
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
        mvA.mv(dx, d_tm);
        cudaError_t e3 = cudaMemcpy(hax.data(), d_tm, m * sizeof(double),
                                    cudaMemcpyDeviceToHost);
        mvAT.mv(dy, d_ata);
        kern_add<<<(n + 255) / 256, 256>>>(n, d_c, d_ata, d_grad);
        cudaError_t e4 = cudaMemcpy(hd.data(), d_grad, n * sizeof(double),
                                    cudaMemcpyDeviceToHost);
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
    double best_err = INF_NAN;
    int stag = 0;
    double *d_xbest, *d_ybest;
    CK(cudaMalloc(&d_xbest, n * sizeof(double)));
    CK(cudaMalloc(&d_ybest, m * sizeof(double)));
    CK(cudaMemcpy(d_xbest, d_x, n * sizeof(double),
                  cudaMemcpyDeviceToDevice));
    CK(cudaMemcpy(d_ybest, d_y, m * sizeof(double),
                  cudaMemcpyDeviceToDevice));
    for (; k <= MAXIT; ++k) {
        iters_run = k;
        mvAT.mv(d_y, d_ata);
        kern_add<<<(n + 255) / 256, 256>>>(n, d_c, d_ata, d_grad);
        const double neg_tau = -tau;
        CK(cublasDcopy(cb, n, d_x, 1, d_tmpn, 1));
        CK(cublasDaxpy(cb, n, &neg_tau, d_grad, 1, d_tmpn, 1));
        kern_clamp<<<(n + 255) / 256, 256>>>(n, d_tmpn, d_cl, d_cu, d_xnew);
        kern_two_minus<<<(n + 255) / 256, 256>>>(n, d_xnew, d_x, d_xbar);
        mvA.mv(d_xbar, d_tm);
        CK(cublasDaxpy(cb, m, &sigma, d_tm, 1, d_y, 1));
        kern_dual<<<(m + 255) / 256, 256>>>(m, d_y, d_rmin, d_rmax, sigma,
                                            d_ynew);
        std::swap(d_x, d_xnew);
        std::swap(d_y, d_ynew);
        double wk = std::min((double)k, 500.0);
        double weta = wk / (wsum + wk);
        wsum += wk;
        kern_mix<<<(n + 255) / 256, 256>>>(n, weta, d_x, d_xavg);
        kern_mix<<<(m + 255) / 256, 256>>>(m, weta, d_y, d_yavg);

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
            bool res_ok = last.pinf <= TOL && last.dinf <= TOL;
            bool gap_ok = std::isfinite(last.gap) && last.gap <= TOL;
            if (res_ok && (!opt.restarts || gap_ok)) break;

            if (opt.restarts) {
            if (mc.pinf > 1e20 || !std::isfinite(mc.pinf)) {
                CK(cudaMemcpy(d_x, d_xbest, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_y, d_ybest, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                tau = sigma = std::sqrt(ts);
                stag = 0;
            } else if (k >= 400 && ea < ec && ea <= 0.36 * best_err) {
                CK(cudaMemcpy(d_x, d_xavg, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_y, d_yavg, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_xbest, d_x, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_ybest, d_y, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                best_err = ea;
                tau = sigma = std::sqrt(ts);
                wsum = 0.0;
                stag = 0;
            } else if (ec < best_err) {
                best_err = ec;
                CK(cudaMemcpy(d_xbest, d_x, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_ybest, d_y, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                stag = 0;
            } else if (ec > 0.8 * best_err) {
                if (++stag >= 40) {
                    CK(cudaMemcpy(d_x, d_xbest, n * sizeof(double),
                                  cudaMemcpyDeviceToDevice));
                    CK(cudaMemcpy(d_y, d_ybest, m * sizeof(double),
                                  cudaMemcpyDeviceToDevice));
                    tau = sigma = std::sqrt(ts);
                    wsum = 0.0;
                    stag = 0;
                }
            }
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
    CK(cudaFree(d_xbest));
    CK(cudaFree(d_ybest));
    double secs = std::chrono::duration<double>(
                      std::chrono::steady_clock::now() - t0)
                      .count();

    if (sol.status != Status::TimeLimit) {
        bool res_ok = last.pinf <= TOL && last.dinf <= TOL;
        bool gap_ok = std::isfinite(last.gap) && last.gap <= TOL;
        if (res_ok && gap_ok) {
            sol.status = Status::NearOptimal;
            sol.message = last_is_avg
                              ? "converged (avg iterate): residuals + gap @tol"
                              : "converged: residuals + gap @tol";
        } else if (res_ok) {
            sol.status = Status::Feasible;
            sol.message = "residuals @tol; duality gap not certified";
        } else {
            sol.status = Status::IterationLimit;
            sol.message = "iteration cap reached before tolerance";
        }
    }

    const double* fx = last_is_avg ? d_xavg : d_x;
    const double* fy = last_is_avg ? d_yavg : d_y;
    sol.objective = last.op;
    sol.pinf = last.pinf;
    sol.dinf = last.dinf;
    sol.rel_gap = last.gap;
    sol.iterations = iters_run;
    sol.solve_time_ms = secs * 1000.0;
    sol.x.resize(n);
    CK(cudaMemcpy(sol.x.data(), fx, n * sizeof(double),
                  cudaMemcpyDeviceToHost));
    sol.y.resize(m);
    CK(cudaMemcpy(sol.y.data(), fy, m * sizeof(double),
                  cudaMemcpyDeviceToHost));

    cudaFree(d_x); cudaFree(d_xnew); cudaFree(d_xbar); cudaFree(d_xavg);
    cudaFree(d_y); cudaFree(d_ynew); cudaFree(d_t); cudaFree(d_yavg);
    cudaFree(d_ata); cudaFree(d_grad); cudaFree(d_tm); cudaFree(d_tmpn);
    cudaFree(d_ap); cudaFree(d_ai); cudaFree(d_ax);
    cudaFree(d_cp); cudaFree(d_ci); cudaFree(d_acx);
    cudaFree(d_c); cudaFree(d_cl); cudaFree(d_cu);
    cudaFree(d_rmin); cudaFree(d_rmax);
    cublasDestroy(cb);
    cusparseDestroy(cs);
    return sol;
}

}  // namespace igaos::pdhg
