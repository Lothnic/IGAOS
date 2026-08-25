#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "dense_lu.hpp"
#include "scaling.hpp"

using namespace igaos::simplex;

int main() {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> ud(-1.0, 1.0);

    int n = 40;
    std::vector<double> A(n * n, 0.0), M(n * n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                M[i * n + j] += ud(rng) * ud(rng);
    for (int i = 0; i < n; ++i) M[i * n + i] += n;

    DenseLU lu;
    lu.factor(M);
    assert(lu.ok);

    std::vector<double> xr(n), xb(n, 0.0);
    for (int i = 0; i < n; ++i) xr[i] = ud(rng);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) xb[i] += M[i * n + j] * xr[j];
    std::vector<double> xs;
    lu.solve(xb, xs);
    double err = 0.0;
    for (int i = 0; i < n; ++i) err = std::max(err, std::fabs(xs[i] - xr[i]));
    printf("solve residual: %.3e\n", err);
    assert(err < 1e-8);

    std::vector<double> bt(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) bt[j] += M[i * n + j] * xr[i];
    std::vector<double> xt;
    lu.solve_transpose(bt, xt);
    double errt = 0.0;
    for (int i = 0; i < n; ++i) errt = std::max(errt, std::fabs(xt[i] - xr[i]));
    printf("transpose residual: %.3e\n", errt);
    assert(errt < 1e-8);

    int m = 12, nc = 10;
    std::vector<int> ap(m + 1, 0), cp(nc + 1, 0);
    std::vector<int> ai, ci;
    std::vector<double> ax, acx;
    auto push = [&](int r, int c, double v) {
        ai.push_back(c);
        ax.push_back(v);
        ++ap[r + 1];
        ci.push_back(r);
        acx.push_back(v);
        ++cp[c + 1];
    };
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < nc; ++j) {
            double v = ud(rng) * std::pow(10.0, (i % 5) - 2);
            if (std::fabs(v) > 1e-4) {
                push(i, j, v);
            }
        }
    for (int i = 0; i < m; ++i) ap[i + 1] += ap[i];
    for (int j = 0; j < nc; ++j) cp[j + 1] += cp[j];

    Scaling sc;
    geometric_mean_scaling(ap, ai, ax, cp, ci, acx, m, nc, 12, sc);
    power_of_two_snap(sc.row);
    power_of_two_snap(sc.col);

    auto entry = [&](int i, int j) -> double {
        for (int p = ap[i]; p < ap[i + 1]; ++p)
            if (ai[p] == j) return ax[p];
        return 0.0;
    };
    auto geo_mean = [&](bool cols) {
        double worst = 1.0;
        int outer = cols ? nc : m;
        for (int o = 0; o < outer; ++o) {
            double prod = 0.0;
            int cnt = 0;
            if (cols) {
                for (int p = cp[o]; p < cp[o + 1]; ++p) {
                    double v = std::fabs(acx[p]) * sc.row[ci[p]] * sc.col[o];
                    if (v > 0) { prod += std::log(v); ++cnt; }
                }
            } else {
                for (int p = ap[o]; p < ap[o + 1]; ++p) {
                    double v = std::fabs(ax[p]) * sc.col[ai[p]] * sc.row[o];
                    if (v > 0) { prod += std::log(v); ++cnt; }
                }
            }
            if (cnt > 1) {
                double gm = std::exp(prod / cnt);
                double dev = std::max(gm, 1.0 / gm);
                worst = std::max(worst, dev);
            }
        }
        return worst;
    };
    double gm_rows = geo_mean(false);
    int total = 0, inband = 0;
    for (int i = 0; i < m; ++i)
        for (int p = ap[i]; p < ap[i + 1]; ++p) {
            double v = std::fabs(ax[p]) * sc.col[ai[p]] * sc.row[i];
            if (v <= 0) continue;
            ++total;
            if (v >= 1e-3 && v <= 1e3) ++inband;
        }
    double frac = total ? (double)inband / total : 1.0;
    fprintf(stderr, "rows gm-dev %.3f | entries in [1e-3,1e3]: %.1f%%\n",
            gm_rows, 100.0 * frac);
    printf("entries within band: %.1f%% (rows gm-dev %.3f)\n", 100.0 * frac,
           gm_rows);
    assert(gm_rows < 4.0);
    assert(frac >= 0.9);

    printf("linalg smoke OK\n");
    return 0;
}
