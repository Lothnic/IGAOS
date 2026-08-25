#pragma once

#include <cmath>
#include <vector>

namespace igaos::simplex {

struct Scaling {
    std::vector<double> row, col;
};

inline void geometric_mean_scaling(const std::vector<int>& ap,
                                   const std::vector<int>& ai,
                                   const std::vector<double>& ax,
                                   const std::vector<int>& cp,
                                   const std::vector<int>& ci,
                                   const std::vector<double>& acx, int m,
                                   int n, int rounds, Scaling& sc) {
    sc.row.assign(m, 1.0);
    sc.col.assign(n, 1.0);

    auto quality = [&]() {
        double q = 0.0;
        for (int j = 0; j < n; ++j) {
            double prod = 0.0;
            int c = 0;
            for (int p = cp[j]; p < cp[j + 1]; ++p) {
                double v = std::fabs(acx[p]) * sc.row[ci[p]] * sc.col[j];
                if (v > 0) { prod += std::log(v); ++c; }
            }
            if (c > 1) q += std::fabs(prod / c);
        }
        for (int i = 0; i < m; ++i) {
            double prod = 0.0;
            int c = 0;
            for (int p = ap[i]; p < ap[i + 1]; ++p) {
                double v = std::fabs(ax[p]) * sc.col[ai[p]] * sc.row[i];
                if (v > 0) { prod += std::log(v); ++c; }
            }
            if (c > 1) q += std::fabs(prod / c);
        }
        return q;
    };
    double prev_q = quality();
    std::vector<double> prow = sc.row, pcol = sc.col;

    for (int r = 0; r < rounds; ++r) {
        for (int j = 0; j < n; ++j) {
            double prod = 0.0;
            int cnt = 0;
            for (int p = cp[j]; p < cp[j + 1]; ++p) {
                double v = std::fabs(acx[p]) * sc.row[ci[p]];
                if (v > 0) {
                    prod += std::log(v);
                    ++cnt;
                }
            }
            if (cnt > 1) {
                double d = std::exp(-prod / cnt);
                if (d > 1e-8 && d < 1e8 && std::isfinite(d)) sc.col[j] *= d;
            }
        }
        bool done = true;
        for (int i = 0; i < m; ++i) {
            double prod = 0.0;
            int cnt2 = 0;
            double mx = 0.0, mn = 1e300;
            for (int p = ap[i]; p < ap[i + 1]; ++p) {
                double v =
                    std::fabs(ax[p]) * sc.col[ai[p]] * sc.row[i];
                if (v > 0) {
                    prod += std::log(v);
                    ++cnt2;
                    mx = std::max(mx, v);
                    mn = std::min(mn, v);
                }
            }
            if (cnt2 > 1) {
                double d = std::exp(-prod / cnt2);
                if (d > 1e-8 && d < 1e8 && std::isfinite(d)) sc.row[i] *= d;
            }
            if (cnt2 > 1 && mn > 0 && mx / mn > 100.0) done = false;
        }
        double q = quality();
        if (q >= prev_q) {
            sc.row = prow;
            sc.col = pcol;
            break;
        }
        prev_q = q;
        prow = sc.row;
        pcol = sc.col;
        if (done) break;
    }
}

inline void power_of_two_snap(std::vector<double>& v) {
    for (double& x : v) {
        if (!(x > 0) || !std::isfinite(x)) continue;
        int e;
        double mm = std::frexp(x, &e);
        x = (mm >= 0.7071067811865476) ? std::ldexp(1.0, e)
                                       : std::ldexp(0.5, e);
    }
}

}  // namespace igaos::simplex
