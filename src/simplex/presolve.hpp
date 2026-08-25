#pragma once

#include "model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace igaos::simplex {

struct PresolveInfo {
    int singletons_tightened = 0;
    int empty_rows = 0;
    double max_bound_shift = 0.0;
};

inline int strengthen_bounds(io::Model& M, PresolveInfo& info,
                             int rounds = 3) {
    const double INF = std::numeric_limits<double>::infinity();
    int changes = 0;
    std::vector<int> row_nz(M.m, 0);
    std::vector<int> row_single_col(M.m, -1);
    std::vector<double> row_single_a(M.m, 0.0);

    auto recount = [&]() {
        std::fill(row_nz.begin(), row_nz.end(), 0);
        std::fill(row_single_col.begin(), row_single_col.end(), -1);
        for (int j = 0; j < M.n; ++j)
            for (int p = M.cp[j]; p < M.cp[j + 1]; ++p) {
                int r = M.ci[p];
                ++row_nz[r];
                if (row_nz[r] == 1) {
                    row_single_col[r] = j;
                    row_single_a[r] = M.acx[p];
                } else if (row_nz[r] == 2) {
                    row_single_a[r] = 0.0;
                }
            }
    };
    recount();

    for (int rd = 0; rd < rounds; ++rd) {
        bool any = false;

        for (int j = 0; j < M.n; ++j) {
            bool is_int =
                !M.integ.empty() && M.integ[j] != 0;
            double nlo = M.cl[j], nup = M.cu[j];
            if (is_int && std::isfinite(nlo)) {
                double c = std::ceil(nlo - 1e-9);
                if (c > nlo) nlo = c;
            }
            if (is_int && std::isfinite(nup)) {
                double f2 = std::floor(nup + 1e-9);
                if (f2 < nup) nup = f2;
            }
            if (nlo > M.cl[j] + 1e-9 || nup < M.cu[j] - 1e-9) {
                info.max_bound_shift =
                    std::max(info.max_bound_shift,
                             std::max(nlo - M.cl[j], M.cu[j] - nup));
                M.cl[j] = nlo;
                M.cu[j] = nup;
                ++changes;
                any = true;
            }
        }

        for (int i = 0; i < M.m; ++i) {
            if (row_nz[i] != 1) continue;
            int j = row_single_col[i];
            double a = row_single_a[i];
            if (j < 0 || std::fabs(a) < 1e-12) continue;
            double d_lo = std::isfinite(M.rmin[i]) ? M.rmin[i] / a : -INF;
            double d_up = std::isfinite(M.rmax[i]) ? M.rmax[i] / a : INF;
            if (a < 0) std::swap(d_lo, d_up);
            double nlo = std::max(M.cl[j], d_lo);
            double nup = std::min(M.cu[j], d_up);
            if (nlo > M.cu[j] + 1e-6 || nup < M.cl[j] - 1e-6)
                return changes;
            if (nlo > M.cl[j] + 1e-9 || nup < M.cu[j] - 1e-9) {
                info.max_bound_shift =
                    std::max(info.max_bound_shift,
                             std::max(std::fabs(nlo - M.cl[j]),
                                      std::fabs(nup - M.cu[j])));
                M.cl[j] = nlo;
                M.cu[j] = nup;
                ++changes;
                any = true;
            }
        }

        if (!any) break;
        if (rd + 1 < rounds) recount();
    }
    info.singletons_tightened = changes;
    return changes;
}

}  // namespace igaos::simplex
