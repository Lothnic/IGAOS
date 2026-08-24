#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include "mps.hpp"

#define CK(x)                                                                       \
    do {                                                                            \
        auto e_ = (x);                                                              \
        if (e_ != cudaSuccess && e_ != CUBLAS_STATUS_SUCCESS &&                     \
            e_ != CUSPARSE_STATUS_SUCCESS) {                                        \
            std::fprintf(stderr, "CUDA error %d at %s:%d\n", (int)e_, __FILE__,     \
                         __LINE__);                                                 \
            std::exit(1);                                                           \
        }                                                                           \
    } while (0)

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

__global__ void kern_sub(int n, const double* a, const double* b, double* z) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) z[i] = a[i] - b[i];
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
        CK(cusparseSpMV(h, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat, vin,
                        &beta, vout, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, dbuf));
    }
};

static Spmv make_spmv(cusparseHandle_t h, int rows, int cols, int nnz,
                      int* rp, int* ri, double* rv) {
    Spmv s;
    s.h = h;
    CK(cusparseCreateCsr(&s.mat, rows, cols, nnz, rp, ri, rv,
                         CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                         CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
    return s;
}

struct Metrics {
    double pinf, dinf, gap, op, od;
    double err() const {
        double g = std::isfinite(gap) ? gap : 0.0;
        return std::max(pinf, std::max(dinf, g));
    }
};

template <class F>
static Metrics evaluate(int m, int n, const std::vector<double>& rmin,
                        const std::vector<double>& rmax,
                        const std::vector<double>& cl,
                        const std::vector<double>& cu,
                        const std::vector<double>& c, const std::vector<double>& x,
                        const std::vector<double>& y, const std::vector<double>& Ax,
                        const std::vector<double>& d, F&& finf_is_bad) {
    const double INF = std::numeric_limits<double>::infinity();
    Metrics mt{};
    double pv = 0.0, pa = 0.0;
    for (int i = 0; i < m; ++i) {
        double lo = std::max(rmin[i] - Ax[i], 0.0);
        double hi = std::max(Ax[i] - rmax[i], 0.0);
        pv = std::max(pv, std::max(lo, hi));
        pa = std::max(pa, std::fabs(std::min(std::max(Ax[i], rmin[i]), rmax[i])));
    }
    mt.pinf = pv / (1.0 + pa);
    double dv = 0.0;
    for (int j = 0; j < n; ++j) {
        double tol = 1e-9 * (1.0 + std::min(std::fabs(cl[j]), std::isfinite(cu[j])
                                                                ? std::fabs(cu[j])
                                                                : std::fabs(cl[j])));
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
    bool dual_bounded = true;
    for (int j = 0; j < n; ++j) {
        double term = (d[j] >= 0) ? d[j] * cl[j] : d[j] * cu[j];
        if (!std::isfinite(term)) { dual_bounded = false; break; }
        od += term;
    }
    mt.od = dual_bounded ? od : -INF;
    mt.gap = dual_bounded ? std::fabs(mt.op - od) / (1.0 + std::fabs(mt.op))
                          : INF;
    return mt;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s model.mps\n", argv[0]); return 2; }
    LP lp;
    parse_mps(argv[1], lp);
    int m = lp.m, n = lp.n, nnz = (int)lp.ax.size();
    std::string name = argv[1];
    auto slash = name.find_last_of('/');
    if (slash != std::string::npos) name = name.substr(slash + 1);
    auto dot = name.find_last_of('.');
    if (dot != std::string::npos) name = name.substr(0, dot);

    cublasHandle_t cb;
    CK(cublasCreate(&cb));
    cusparseHandle_t cs;
    CK(cusparseCreate(&cs));

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
    CK(cudaMemcpy(d_ap, lp.ap.data(), (m + 1) * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_ai, lp.ai.data(), nnz * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_ax, lp.ax.data(), nnz * sizeof(double), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cp, lp.cp.data(), (n + 1) * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_ci, lp.ci.data(), nnz * sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_acx, lp.acx.data(), nnz * sizeof(double), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_c, lp.c.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cl, lp.cl.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cu, lp.cu.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_rmin, lp.rmin.data(), m * sizeof(double), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_rmax, lp.rmax.data(), m * sizeof(double), cudaMemcpyHostToDevice));

    double *d_x, *d_xnew, *d_xbar, *d_xavg, *d_y, *d_ynew, *d_t, *d_yavg, *d_ata,
        *d_grad, *d_axb, *d_tmpn, *d_tm;
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
    CK(cudaMalloc(&d_axb, m * sizeof(double)));
    CK(cudaMalloc(&d_tmpn, n * sizeof(double)));
    CK(cudaMalloc(&d_tm, m * sizeof(double)));
    CK(cudaMemset(d_x, 0, n * sizeof(double)));
    CK(cudaMemset(d_y, 0, m * sizeof(double)));
    CK(cudaMemset(d_xavg, 0, n * sizeof(double)));
    CK(cudaMemset(d_yavg, 0, m * sizeof(double)));

    Spmv mvA = make_spmv(cs, m, n, nnz, d_ap, d_ai, d_ax);
    Spmv mvAT = make_spmv(cs, n, m, nnz, d_cp, d_ci, d_acx);
    cusparseDnVecDescr_t vn_in, vn_out, vm_in, vm_out;
    CK(cusparseCreateDnVec(&vn_in, n, d_tmpn, CUDA_R_64F));
    CK(cusparseCreateDnVec(&vm_out, m, d_tm, CUDA_R_64F));
    CK(cusparseCreateDnVec(&vm_in, m, d_y, CUDA_R_64F));
    CK(cusparseCreateDnVec(&vn_out, n, d_ata, CUDA_R_64F));
    size_t bA = 0, bAT = 0;
    double one = 1.0, zero = 0.0;
    CK(cusparseSpMV_bufferSize(cs, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                               mvA.mat, vn_in, &zero, vm_out, CUDA_R_64F,
                               CUSPARSE_SPMV_ALG_DEFAULT, &bA));
    CK(cusparseSpMV_bufferSize(cs, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                               mvAT.mat, vm_in, &zero, vn_out, CUDA_R_64F,
                               CUSPARSE_SPMV_ALG_DEFAULT, &bAT));
    mvA.buf = std::max(bA, (size_t)1);
    mvAT.buf = std::max(bAT, (size_t)1);
    CK(cudaMalloc(&mvA.dbuf, mvA.buf));
    CK(cudaMalloc(&mvAT.dbuf, mvAT.buf));
    mvA.vin = vn_in; mvA.vout = vm_out;
    mvAT.vin = vm_in; mvAT.vout = vn_out;

    double lam2 = 1.0;
    {
        std::vector<double> hv(n, 1.0);
        std::vector<double> hout(std::max(m, n), 0.0);
        for (int j = 0; j < n; ++j) hv[j] = (double)((j * 37) % 17) - 8.0;
        double* dv; double* dout;
        CK(cudaMalloc(&dv, n * sizeof(double)));
        CK(cudaMalloc(&dout, std::max(m, n) * sizeof(double)));
        CK(cudaMemcpy(dv, hv.data(), n * sizeof(double), cudaMemcpyHostToDevice));
        mvA.mv(dv, dout);
        CK(cudaMemcpy(hout.data(), dout, m * sizeof(double),
                      cudaMemcpyDeviceToHost));
        double err1 = 0.0;
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int p = lp.ap[i]; p < lp.ap[i + 1]; ++p)
                s += lp.ax[p] * hv[lp.ai[p]];
            err1 = std::max(err1, std::fabs(s - hout[i]));
        }
        mvAT.mv(dv, dout);
        CK(cudaMemcpy(hout.data(), dout, n * sizeof(double),
                      cudaMemcpyDeviceToHost));
        double err2 = 0.0;
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int p = lp.cp[j]; p < lp.cp[j + 1]; ++p)
                s += lp.acx[p] * hv[lp.ci[p]];
            err2 = std::max(err2, std::fabs(s - hout[j]));
        }
        std::fprintf(stderr, "selftest spmvA=%.3e spmvAT=%.3e\n", err1, err2);
        cudaFree(dv); cudaFree(dout);
    }
    {
        std::vector<double> hv(n, 1.0 / std::sqrt((double)n));
        double *dv, *dw;
        CK(cudaMalloc(&dv, n * sizeof(double)));
        CK(cudaMalloc(&dw, n * sizeof(double)));
        CK(cudaMemcpy(dv, hv.data(), n * sizeof(double), cudaMemcpyHostToDevice));
        for (int it = 0; it < 30; ++it) {
            mvA.mv(dv, d_tm);
            mvAT.mv(d_tm, dw);
            double nw = 0.0;
            CK(cublasDnrm2(cb, n, dw, 1, &nw));
            if (nw < 1e-300) break;
            double inv = 1.0 / nw;
            CK(cublasDscal(cb, n, &inv, dw, 1));
            CK(cudaMemcpy(dv, dw, n * sizeof(double), cudaMemcpyDeviceToDevice));
        }
        mvA.mv(dv, d_tm);
        double na = 0.0;
        CK(cublasDnrm2(cb, m, d_tm, 1, &na));
        lam2 = na * na;
        if (!(lam2 > 0) || !std::isfinite(lam2)) lam2 = 1.0;
        cudaFree(dv); cudaFree(dw);
    }
    double ts_max = 0.9 / lam2, ts = 0.5 * ts_max;
    double tau = std::sqrt(ts), sigma = tau;
    std::fprintf(stderr, "lam2=%.6g ts_max=%.6g\n", lam2, ts_max);

    std::vector<double> hx(n), hy(m), hax(m), hd(n);
    auto snapshot_metrics = [&](const double* dx, const double* dy,
                                Metrics& mt) {
        CK(cudaMemcpy(hx.data(), dx, n * sizeof(double), cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(hy.data(), dy, m * sizeof(double), cudaMemcpyDeviceToHost));
        mvA.mv(dx, d_tm);
        CK(cudaMemcpy(hax.data(), d_tm, m * sizeof(double), cudaMemcpyDeviceToHost));
        mvAT.mv(dy, d_ata);
        kern_add<<<(n + 255) / 256, 256>>>(n, d_c, d_ata, d_grad);
        CK(cudaMemcpy(hd.data(), d_grad, n * sizeof(double),
                      cudaMemcpyDeviceToHost));
        mt = evaluate(m, n, lp.rmin, lp.rmax, lp.cl, lp.cu, lp.c, hx, hy, hax,
                      hd, [](bool){ return 1e30; });
    };

    double *d_xbest, *d_ybest;
    CK(cudaMalloc(&d_xbest, n * sizeof(double)));
    CK(cudaMalloc(&d_ybest, m * sizeof(double)));

    const int CHECK = 25;
    const int MAXIT = 500000;
    const double TOL = 1e-4, TOL_TIGHT = 1e-6;
    double best_err = std::numeric_limits<double>::infinity();
    double achieved = std::numeric_limits<double>::infinity();
    Metrics last{};
    int stag = 0, iters = 0, tight_hit = -1;
    double wsum = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    const double budget_s = 120.0;
    std::string status = "MAXITER";

    for (int k = 1; k <= MAXIT; ++k) {
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
        iters = k;

        if (k % CHECK == 0) {
            Metrics mc, ma;
            snapshot_metrics(d_x, d_y, mc);
            snapshot_metrics(d_xavg, d_yavg, ma);
            double ec = (std::isfinite(mc.pinf) && std::isfinite(mc.dinf))
                            ? mc.err() : 1e30;
            double ea = (std::isfinite(ma.pinf) && std::isfinite(ma.dinf))
                            ? ma.err() : 1e30;
            last = (ec <= ea) ? mc : ma;
            achieved = std::min(achieved, last.err());
            bool cur_conv = ec <= TOL;
            bool avg_conv = ea <= TOL;
            if ((cur_conv || avg_conv) && tight_hit < 0) tight_hit = k;

            if (last.err() <= TOL) {
                status = (ec <= ea) ? "OPT@1e-4" : "OPT@1e-4(avg)";
                break;
            }
            auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration<double>(now - t0).count() > budget_s) {
                status = "TIMEOUT";
                break;
            }
        }
    }

    double secs = std::chrono::duration<double>(
                      std::chrono::steady_clock::now() - t0)
                      .count();
    if (status == "MAXITER" && tight_hit > 0) status = "OPT@1e-6(late)";
    printf("%s,%d,%d,%d,%d,%.1f,%.10g,%.10g,%.3e,%.3e,%.3e,%s\n", name.c_str(),
           m, n, nnz, iters, secs * 1000.0, last.op, last.od, last.gap,
           last.pinf, last.dinf, status.c_str());
    return 0;
}
