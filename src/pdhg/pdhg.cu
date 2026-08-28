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
    double prow;  // max per-row violation / (1 + row's own scale)
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
    double pv = 0.0, pa = 0.0, prow = 0.0;
    for (int i = 0; i < m; ++i) {
        pv = std::max(pv, std::max(rmin[i] - Ax[i], Ax[i] - rmax[i]));
        pa = std::max(pa,
                      std::fabs(std::min(std::max(Ax[i], rmin[i]), rmax[i])));
        // per-row relative violation: pinf's (1+pa) normalizer hides big
        // absolute violations on large-scale rows (adlittle: 0.1 abs
        // violation at pinf 4e-5 bought 383 of objective). Certification
        // requires this per-row measure instead.
        double viol = std::max(
            0.0, std::max(rmin[i] - Ax[i], Ax[i] - rmax[i]));
        double sc = std::fabs(Ax[i]);
        if (std::isfinite(rmin[i]))
            sc = std::max(sc, std::fabs(rmin[i]));
        if (std::isfinite(rmax[i]))
            sc = std::max(sc, std::fabs(rmax[i]));
        prow = std::max(prow, viol / (1.0 + sc));
    }
    mt.pinf = std::max(pv, 0.0) / (1.0 + pa);
    mt.prow = prow;
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
// Cholesky factorization of a symmetric PD matrix in place (lower
// triangle, row-major, stride nr). Returns false if not PD.
bool chol_factor(std::vector<double>& G, int nr) {
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s = G[i * nr + j];
            for (int t = 0; t < j; ++t) s -= G[i * nr + t] * G[j * nr + t];
            if (i == j) {
                if (!(s > 0)) return false;
                G[i * nr + i] = std::sqrt(s);
            } else {
                G[i * nr + j] = s / G[j * nr + j];
            }
        }
    }
    return true;
}

// Solve with a chol_factor'd matrix; rhs in place.
void chol_apply(const std::vector<double>& G, int nr,
                std::vector<double>& rhs) {
    for (int i = 0; i < nr; ++i) {
        double s = rhs[i];
        for (int t = 0; t < i; ++t) s -= G[i * nr + t] * rhs[t];
        rhs[i] = s / G[i * nr + i];
    }
    for (int i = nr - 1; i >= 0; --i) {
        double s = rhs[i];
        for (int t = i + 1; t < nr; ++t) s -= G[t * nr + i] * rhs[t];
        rhs[i] = s / G[i * nr + i];
    }
}

bool dual_repair(const io::Model& lp, const std::vector<double>& x,
                 const std::vector<double>& Ax,
                 const std::vector<double>& y_it, double col_tol,
                 std::vector<double>& y) {
    const int m = lp.m, n = lp.n;
    const double NL = -std::numeric_limits<double>::infinity();
    const double PL = std::numeric_limits<double>::infinity();
    std::vector<double> lo(m), hi(m);
    for (int i = 0; i < m; ++i) {
        bool fin_lo = std::isfinite(lp.rmin[i]);
        bool fin_hi = std::isfinite(lp.rmax[i]);
        if (fin_lo && fin_hi) { lo[i] = NL; hi[i] = PL; }
        else if (fin_hi)      { lo[i] = 0.0; hi[i] = PL; }
        else if (fin_lo)      { lo[i] = NL; hi[i] = 0.0; }
        else                  { lo[i] = hi[i] = 0.0; }
    }
    // column constraint types: 0: d>=0 (cu=inf), 1: d<=0 (cl=-inf),
    // 2: d==0 (both inf), 3: free (both finite)
    std::vector<char> ct(n, 3);
    std::vector<double> cn2(n, 0.0);
    std::vector<int> ji;
    for (int j = 0; j < n; ++j) {
        bool up_inf = !std::isfinite(lp.cu[j]);
        bool lo_inf = !std::isfinite(lp.cl[j]);
        ct[j] = (up_inf && lo_inf) ? 2 : up_inf ? 0 : lo_inf ? 1 : 3;
        for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
            cn2[j] += lp.acx[p] * lp.acx[p];
        // interior: x strictly off both bounds -> d_j must be 0
        auto btol = [&](double tol, double b) {
            return tol * (1.0 + (std::isfinite(b) ? std::fabs(b) : 0.0));
        };
        if (x[j] - lp.cl[j] > btol(col_tol, lp.cl[j]) &&
            lp.cu[j] - x[j] > btol(col_tol, lp.cu[j]))
            ji.push_back(j);
    }
    const int k = (int)ji.size();
    // G = A_JI A_JI' + ridge I, factored once
    std::vector<double> G;
    if (k > 0) {
        G.assign(m * m, 0.0);
        for (int jj = 0; jj < k; ++jj) {
            const int j = ji[jj];
            for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p) {
                const int i1 = lp.ci[p];
                const double v1 = lp.acx[p];
                for (int q = lp.cp[j]; q < lp.cp[j + 1]; ++q)
                    G[i1 * m + lp.ci[q]] += v1 * lp.acx[q];
            }
        }
        double tr = 0.0;
        for (int i = 0; i < m; ++i) tr += G[i * m + i];
        double ridge = tr > 0 ? 1e-10 * tr / m : 1e-14;
        bool ok = false;
        for (int att = 0; att < 6 && !ok; ++att) {
            std::vector<double> Gc = G;
            for (int i = 0; i < m; ++i) Gc[i * m + i] += ridge;
            ok = chol_factor(Gc, m);
            if (ok) G = Gc;
            else ridge = ridge > 0 ? ridge * 100.0 : 1e-14;
        }
        if (!ok) return false;
    }
    double best_gap = std::numeric_limits<double>::infinity();
    bool have = false;
    std::vector<double> d(n), rhs(m), yc(m), ybest;
    double best_viol = std::numeric_limits<double>::infinity();
    std::vector<double> yviol;
    for (int anch = 0; anch < 2; ++anch) {
        yc = anch == 0 ? y_it : std::vector<double>(m, 0.0);
        for (int i = 0; i < m; ++i)
            yc[i] = std::min(std::max(yc[i], lo[i]), hi[i]);
        // decay-rate abort: hopeless projections (wrong active set, far
        // anchor) plateau; without this each failed attempt burns
        // ~0.3s of host time and march-dominated runs (sc205) lose ~40%
        // of their iteration budget
        double vref = -1.0;
        for (int it = 0; it < 300; ++it) {
            // affine projection: e = G^-1 (-A_JI d_JI)
            if (k > 0) {
                std::fill(rhs.begin(), rhs.end(), 0.0);
                double res2 = 0.0;
                for (int jj = 0; jj < k; ++jj) {
                    const int j = ji[jj];
                    double dj = lp.c[j];
                    for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                        dj += lp.acx[p] * yc[lp.ci[p]];
                    res2 += dj * dj;
                    if (dj != 0.0)
                        for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                            rhs[lp.ci[p]] -= lp.acx[p] * dj;
                }
                if (res2 > 0.0) {
                    chol_apply(G, m, rhs);
                    for (int i = 0; i < m; ++i) yc[i] += rhs[i];
                }
            }
            // Hildreth halfspace fixes for sign-constrained columns
            for (int j = 0; j < n; ++j) {
                if (ct[j] == 3 || cn2[j] == 0.0) continue;
                double dj = lp.c[j];
                for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                    dj += lp.acx[p] * yc[lp.ci[p]];
                double v = 0.0;
                if (ct[j] == 0 && dj < 0.0) v = -dj;
                else if (ct[j] == 1 && dj > 0.0) v = -dj;
                else if (ct[j] == 2 && dj != 0.0) v = -dj;
                if (v != 0.0) {
                    double st = v / cn2[j];
                    for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                        yc[lp.ci[p]] += st * lp.acx[p];
                }
            }
            // row cone
            for (int i = 0; i < m; ++i)
                yc[i] = std::min(std::max(yc[i], lo[i]), hi[i]);
            for (int j = 0; j < n; ++j) {
                double dj = lp.c[j];
                for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                    dj += lp.acx[p] * yc[lp.ci[p]];
                d[j] = dj;
            }
            if (it % 25 == 0) {
                double viol = 0.0;
                for (int j = 0; j < n; ++j) {
                    if (ct[j] == 3) continue;
                    double v = 0.0;
                    if (ct[j] == 0 && d[j] < 0) v = -d[j];
                    else if (ct[j] == 1 && d[j] > 0) v = d[j];
                    else if (ct[j] == 2) v = std::fabs(d[j]);
                    viol = std::max(viol, v);
                }
                for (int i = 0; i < m; ++i) {
                    double v = 0.0;
                    if (lo[i] == 0.0 && hi[i] == 0.0)
                        v = std::fabs(yc[i]);
                    else if (hi[i] == 0.0 && yc[i] > 0) v = yc[i];
                    else if (lo[i] == 0.0 && yc[i] < 0) v = -yc[i];
                    viol = std::max(viol, v);
                }
                if (it == 25) vref = viol;
                if (viol < best_viol) {
                    best_viol = viol;
                    yviol = yc;
                }
                if (it == 100 && vref > 0 && viol > 0.25 * vref) break;
            }
            Metrics t = evaluate(m, n, lp.rmin, lp.rmax, lp.cl, lp.cu, lp.c,
                                 x, yc, Ax, d);
            if (std::isfinite(t.gap) && t.gap < best_gap) {
                best_gap = t.gap;
                ybest = yc;
                have = true;
            }
        }
    }
    // no certifying dual found; still expose the best-by-violation POCS
    // point so the caller can use it as a dual-restart anchor (grow22
    // basin escape)
    y = have ? ybest : yviol;
    return have;
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
    double *d_xb, *d_yb;  // incumbent snapshot
    CK(cudaMalloc(&d_xb, n * sizeof(double)));
    CK(cudaMalloc(&d_yb, m * sizeof(double)));
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
    const double prod0 = tau * sigma;  // step product stays pinned
    const double rt0 = std::sqrt(prod0);

    const int CHECK = 50;
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
    Metrics inc{};
    bool have_inc = false;
    double inc_op = 0.0;
    long last_rollback = 0;
    double prev_pinf = std::numeric_limits<double>::infinity();
    double prev_win_op = 0.0;
    double boom_cap = std::numeric_limits<double>::infinity();
    std::vector<double> yprev(m, 0.0), xprev(n, 0.0);
    bool have_prev = false;
    long next_repair = 5000;   // dual-repair attempt schedule (backoff)
    bool repaired = false;
    long repair_spacing = 2000;
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
            if (opt.verbosity > 1 && k % 1000 == 0) {
                double dy = 0.0, dx = 0.0;
                if (have_prev) {
                    for (int i = 0; i < m; ++i)
                        dy = std::max(dy, std::fabs(hy[i] - yprev[i]));
                    for (int j = 0; j < n; ++j)
                        dx = std::max(dx, std::fabs(hx[j] - xprev[j]));
                }
                std::fprintf(stderr,
                             "[dbg] k=%ld tau=%.3e sig=%.3e dy=%.3e dx=%.3e "
                             "op=%.6g mc.pinf=%.2e ma.pinf=%.2e mc.dinf=%.2e "
                             "ma.dinf=%.2e\n",
                             k, tau, sigma, dy, dx, last.op, mc.pinf, ma.pinf,
                             mc.dinf, ma.dinf);
                std::copy(hy.begin(), hy.end(), yprev.begin());
                std::copy(hx.begin(), hx.end(), xprev.begin());
                have_prev = true;
            }
            // Honest termination: residuals AND a certified duality gap,
            // on either the current or the averaged iterate. Residual-only
            // acceptance returns wrong objectives (sc205, adlittle):
            // dinf/pinf do not bound the objective error.
            auto passes = [&](const Metrics& t) {
                return t.pinf <= TOL && t.dinf <= TOL && t.prow <= TOL &&
                       std::isfinite(t.gap) && t.gap <= TOL;
            };
            if (passes(mc)) { last = mc; last_is_avg = false; break; }
            if (passes(ma)) { last = ma; last_is_avg = true; break; }

            // Dual repair when only certification lags: a candidate's
            // residuals are at tolerance but its gap is not certifiable
            // (the iterate's y is ~1e-6 short of exact dual feasibility at
            // the infinite bounds, so the Lagrangian bound is -inf).
            // Exactly project the dual onto the feasible region (strictly
            // convex QP, dual active-set) and re-evaluate the bound. Weak
            // duality holds for ANY y and evaluate() re-verifies every
            // sign before a bound is accepted, so a failed projection can
            // only fail to certify — never certify a false gap.
            // ponytail: dense working-set solves on host; the dinf gate
            // keeps march-dominated runs (sc205, grow22) from paying for
            // hopeless attempts.
            if (k >= 5000 && k >= next_repair) {
                next_repair = k + repair_spacing;
                repair_spacing = std::min(repair_spacing * 3, 30000L);
                for (int c = 1; c <= 2; ++c) {
                    const Metrics& t = c == 1 ? ma : inc;
                    if (c == 2 && !have_inc) break;
                    if (!(t.pinf <= TOL && t.prow <= TOL &&
                          t.dinf <= 10.0 * TOL &&
                          !(std::isfinite(t.gap) && t.gap <= TOL)))
                        continue;
                    std::vector<double> rx(n), ry_it(m), ry(m);
                    const double* dxs = c == 1 ? d_xavg : d_xb;
                    double* dys = c == 1 ? d_yavg : d_y;
                    CK(cudaMemcpy(rx.data(), dxs, n * sizeof(double),
                                  cudaMemcpyDeviceToHost));
                    for (int j = 0; j < n; ++j) rx[j] *= dc[j];
                    CK(cudaMemcpy(ry_it.data(), dys, m * sizeof(double),
                                  cudaMemcpyDeviceToHost));
                    for (int i = 0; i < m; ++i)
                        ry_it[i] = gamma * ry_it[i] / dr[i];
                    std::vector<double> rAx(m, 0.0);
                    for (int j = 0; j < n; ++j)
                        for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                            rAx[lp.ci[p]] += lp.acx[p] * rx[j];
                    // col_tol ladder: the averaged point keeps small-x
                    // nonbasic columns interior; counting them makes the
                    // equality system inconsistent. A looser tier drops
                    // them and leaves the consistent basis system.
                    static const double ct[2] = {1e-6, 1e-1};
                    bool ok_rep = false;
                    Metrics mr{};
                    for (int tier = 0; tier < 2 && !ok_rep; ++tier) {
                        if (!dual_repair(lp, rx, rAx, ry_it, ct[tier], ry))
                            continue;
                        std::vector<double> rd(n);
                        for (int j = 0; j < n; ++j) {
                            rd[j] = lp.c[j];
                            for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                                rd[j] += lp.acx[p] * ry[lp.ci[p]];
                        }
                        Metrics t = evaluate(m, n, lp.rmin, lp.rmax, lp.cl,
                                             lp.cu, lp.c, rx, ry, rAx, rd);
                        // Acceptance for a REPAIRED pair is pinf + prow +
                        // gap: the gap is a valid weak-duality bound for
                        // any y, so together with x's feasibility it IS
                        // the certificate — the repaired dual's own dinf
                        // is redundant (it is not the iterate's residual).
                        if (t.pinf <= TOL && t.prow <= TOL &&
                            std::isfinite(t.gap) && t.gap <= TOL) {
                            ok_rep = true;
                            mr = t;
                        }
                    }
                    if (!ok_rep) continue;
                    // adopt: write the repaired dual back so the standard
                    // output path reports it verbatim
                    {
                        std::vector<double> sy(m);
                        for (int i = 0; i < m; ++i)
                            sy[i] = ry[i] * dr[i] / gamma;
                        CK(cudaMemcpy(dys, sy.data(), m * sizeof(double),
                                      cudaMemcpyHostToDevice));
                    }
                    if (c == 2)
                        CK(cudaMemcpy(d_x, d_xb, n * sizeof(double),
                                      cudaMemcpyDeviceToDevice));
                    if (opt.verbosity > 1)
                        std::fprintf(stderr,
                                     "[pdhg] dual repair c=%d @k=%ld "
                                     "gap=%.3e\n",
                                     c, k, mr.gap);
                    last = mr;
                    last_is_avg = (c == 1);
                    repaired = true;
                    break;
                }
                if (repaired) break;
            }

            // Step-size controller. In the march regime the product
            // tau*sigma only ever shrinks below prod0 (growth re-pins it);
            // shrinking tau while raising sigma with the product pinned is
            // the divergent corner — the huge dual step amplifies the xbar
            // extrapolation oscillation (grow22).
            //
            // March regime: the objective is still improving, the gap is
            // far from certifiable (pinf << dinf), and the iterate is
            // crawling along a feasible face — the march speed depends
            // strongly and non-monotonically on tau (grow22: tau~1.7*rt0
            // marches ~170x faster than rt0, tau~7*rt0 is slow again —
            // resonance with the box-clipping pattern). Ride the edge:
            // grow tau while the averaged iterate stays deep-feasible,
            // back off when it degrades. The gate is objective progress,
            // not pinf level, so it persists through the current iterate's
            // ~1e-3 oscillation and exits by itself at convergence. The
            // old relative-imbalance rule suppressed tau in exactly this
            // regime and slowed sc205/grow22 >10x.
            //
            // Outside the march regime: the proven PDLP-style imbalance
            // rule, plus an immediate shrink when pinf rises fast.
            if (k >= 1000 && std::isfinite(last.pinf) && last.pinf > 0 &&
                std::isfinite(last.dinf) && last.dinf > 0) {
                if (k % 500 == 0) {
                    // best-objective iterate: monotone through the march,
                    // unlike the err-selected `last` whose op oscillates and
                    // spuriously exits march mode mid-run (grow22 k=125k)
                    double best_op = std::isfinite(mc.op) ? mc.op : ma.op;
                    if (std::isfinite(ma.op) && ma.op < best_op)
                        best_op = ma.op;
                    double op_gain = prev_win_op - best_op;
                    bool marching = last.pinf < 0.1 * last.dinf &&
                                    op_gain > 1e-10 *
                                                   (1.0 + std::fabs(
                                                              last.op));
                    if (marching) {
                        if (ma.pinf > 1e-3) {
                            // feasibility degrading: pull tau down WITHOUT
                            // raising sigma (product shrinks = strictly
                            // safer); a pinned product here would explode
                            // sigma and the dual amplifies the xbar
                            // extrapolation into divergence (grow22)
                            tau = std::max(tau * 0.8, 0.3 * rt0);
                        } else if (op_gain <
                                   1e-6 * (1.0 + std::fabs(best_op))) {
                            // arrival: the march has slowed near its target;
                            // tighten the orbit (sigma unchanged — see
                            // above) so the endgame settles instead of
                            // oscillating (kb2/sc50a certification)
                            tau = std::max(tau * 0.7, 0.3 * rt0);
                        } else if (ma.pinf < 1e-4) {
                            tau = std::min(tau * 1.5, 8.0 * rt0);
                            sigma = prod0 / tau;
                        }
                        // boom clamp: never re-enter the step size that
                        // blew the state up (grow22 divergence loop);
                        // relax slowly so the controller can re-probe
                        tau = std::min(tau, boom_cap);
                        boom_cap *= 1.02;
                    } else {
                        double tgt = std::clamp(last.pinf / last.dinf, 0.1,
                                                10.0);
                        double cur = tau / sigma;
                        double nxt =
                            cur * std::clamp(tgt / cur, 0.7, 1.4);
                        nxt = std::clamp(nxt, 0.1, 10.0);
                        tau = std::sqrt(prod0 * nxt);
                    }
                    sigma = prod0 / tau;
                    // ponytail: last.op (not best_op) — the two variants
                    // measure slightly different window gains and sc205's
                    // chaotic endgame is sensitive to which one gates
                    // march mode; this one is the measured-good choice
                    prev_win_op = last.op;
                }
                if (ma.pinf >= 1e-2 &&
                    std::isfinite(prev_pinf) && ma.pinf > 4.0 * prev_pinf) {
                    tau = std::max(tau * 0.5, 0.3 * rt0);
                    sigma = prod0 / tau;
                }
                prev_pinf = ma.pinf;
            }

            // Incumbent: best (lowest) objective among deep-feasible
            // iterates, current or averaged. Feasible points cannot beat
            // the true optimum, so this is an honest bound to report on
            // time-limit. Also the rollback anchor when the state blows up.
            const double* inc_x_src = nullptr;
            const double* inc_y_src = nullptr;
            for (int c = 0; c < 2; ++c) {
                const Metrics& t = c ? ma : mc;
                // 1e-4: tight enough that infeasibility cannot buy more than
                // ~1e-4-ish of "superoptimal" objective, loose enough to
                // catch the feasible dips the endgame oscillation produces.
                if (t.pinf <= 1e-4 && std::isfinite(t.op) &&
                    (!have_inc || t.op < inc_op)) {
                    have_inc = true;
                    inc_op = t.op;
                    inc = t;
                    inc_x_src = c ? d_xavg : d_x;
                    inc_y_src = c ? d_yavg : d_y;
                }
            }
            if (inc_x_src) {
                CK(cudaMemcpy(d_xb, inc_x_src, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_yb, inc_y_src, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
            }
            // Self-heal: a catastrophically blown-up state (pinf >= 0.5)
            // never recovers the march regime within the budget — restore
            // the incumbent and re-anchor averaging there with safer steps,
            // and clamp tau to half the size that blew up so the restored
            // state does not instantly re-diverge (grow22 divergence loop).
            // Milder excursions are left alone: the shrunk tau lets them
            // re-settle on their own and they keep their objective gains.
            if (mc.pinf >= 0.5 && have_inc && k - last_rollback > 2000) {
                CK(cudaMemcpy(d_x, d_xb, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_y, d_yb, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_xavg, d_xb, n * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                CK(cudaMemcpy(d_yavg, d_yb, m * sizeof(double),
                              cudaMemcpyDeviceToDevice));
                boom_cap = std::min(boom_cap, 0.5 * tau);
                tau = 0.3 * rt0;  // balanced-small: the original stable corner
                sigma = prod0 / tau;
                hn = 0;
                since_restart = 0;
                last_rollback = k;
                if (opt.verbosity > 1)
                    std::fprintf(stderr, "[pdhg] rollback @k=%ld to inc\n", k);
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

    // Final status: pinf + prow + a certified duality gap. The gap is a
    // weak-duality bound, so together with x's feasibility it bounds the
    // objective error — dinf is redundant for the guarantee (the loop's
    // passes() still keeps it as a stricter early-exit gate for iterate
    // pairs, but a repaired dual's dinf measures a different y).
    auto certified = [&](const Metrics& t) {
        return t.pinf <= TOL && t.prow <= TOL && std::isfinite(t.gap) &&
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
    // Uncertified end-of-run: report the incumbent if its feasible
    // objective beats the last iterate, or if the last iterate is not even
    // feasible at tolerance (e.g. a diverged tail — grow22). Honest — the
    // incumbent's own residuals/gap are reported alongside, status stays
    // time-limit.
    if (!certified(last) && have_inc && inc.pinf <= TOL &&
        std::isfinite(inc.op) &&
        (inc.op < last.op || last.pinf > TOL)) {
        last = inc;
        last_is_avg = false;
        fx = d_xb;
        fy = d_yb;
    }
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
    cudaFree(d_xb); cudaFree(d_yb);
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
