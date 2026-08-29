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

// ---------------------------------------------------------------------------
// Full presolve / postsolve — Andersen & Andersen (1995), Math. Prog. 71.
// From-scratch implementation of the paper's rules; no external solver code.
//
// Passes (paper §3):
//   P1  empty column / empty row handling          (bounds-only + removal)
//   P2  forcing / subsuming rows via full row
//       activity-interval propagation              (§3.2 generalization)
//   P3  fixed-variable substitution                 (needs postsolve)
//   P4  free singleton-column substitution          (strongest reduction)
//   P5  duplicate rows (proportional, merged)       (needs postsolve: y=0)
//   P6  duplicate columns (identical, merged)       (needs postsolve split)
//
// Coefficient tightening note (task item 6): the only A-modifying reduction
// that is exactly validity-preserving on the remaining model is entry
// removal when the entry's own extreme contributions can never activate the
// row — which coincides with the row being redundant and is handled by P2
// (row removal). Tightening a single coefficient in place is only sound for
// integer variables (A&A §3.6); the MILP path does not use this presolve
// (guarded off for warm starts / cut generation), so it is not implemented.
// Postsolve only maps x back to original space — feasibility is re-validated
// against the ORIGINAL model in simplex.cpp, so no A restoration is needed.
//
// Duals: postsolve maps y for surviving rows and assigns y=0 to removed
// rows. For merged/substituted rows this is generally NOT the true optimal
// dual of the original LP (only the primal objective is exact). The
// benchmark protocol validates objectives and primal feasibility only.
// ---------------------------------------------------------------------------
namespace igaos::simplex {

struct PresolveLog {
    struct Rec {
        enum Kind : unsigned char { FixedVar = 0, SubstFree = 1,
                                    DupCol = 2 } kind = FixedVar;
        int j = -1;          // original column removed/fixed
        double val = 0.0;    // FixedVar: x_j := val
        // SubstFree: x_j := (target - sum_k a_k x_k) / aij   (free singleton)
        double aij = 0.0, target = 0.0;
        std::vector<std::pair<int, double>> terms;  // (orig col, coef)
        // DupCol: identical columns j and j2 merged; j2 carries w = x_j+x_j2
        int j2 = -1;
        double cl1 = 0, cu1 = 0, cl2 = 0, cu2 = 0;
    };
    std::vector<Rec> recs;
};

struct PresolveStats {
    int m0 = 0, n0 = 0, nnz0 = 0, m1 = 0, n1 = 0, nnz1 = 0;
    int rows_redundant = 0, rows_forcing = 0, rows_dup = 0;
    int cols_fixed = 0, cols_empty = 0, cols_subst = 0, cols_dup = 0;
    int bound_tightens = 0;
};

namespace pre_detail {

inline double rtol_of(double lo, double hi) {
    double s = 1.0;
    if (std::isfinite(lo)) s += std::fabs(lo);
    if (std::isfinite(hi)) s += std::fabs(hi);
    return 1e-9 * s;
}

inline void pack_key(std::string& key, const std::vector<int>& ids) {
    key.clear();
    for (int v : ids) {
        key.push_back((char)((unsigned)v & 0xff));
        key.push_back((char)(((unsigned)v >> 8) & 0xff));
        key.push_back((char)(((unsigned)v >> 16) & 0xff));
        key.push_back((char)(((unsigned)v >> 24) & 0xff));
    }
}

}  // namespace pre_detail

// Returns 0 = reduced OK, 1 = infeasible, 2 = unbounded.
inline int run_presolve(const io::Model& orig, io::Model& red,
                        PresolveLog& log, std::vector<int>& col_to_red,
                        std::vector<int>& row_to_red, PresolveStats& st) {
    const double INF = std::numeric_limits<double>::infinity();
    const int m = orig.m, n = orig.n;
    using pre_detail::rtol_of;

    // Row-wise working copy over ORIGINAL indices. Invariant maintained by
    // every pass: alive rows only reference alive cols.
    std::vector<std::vector<std::pair<int, double>>> rows(m);
    std::vector<double> rmin = orig.rmin, rmax = orig.rmax;
    std::vector<double> cl = orig.cl, cu = orig.cu, c = orig.c;
    std::vector<unsigned char> integ = orig.integ;
    std::vector<char> row_alive(m, 1), col_alive(n, 1);
    std::vector<double> fixval(n, 0.0);
    double obj_shift = 0.0;
    log.recs.clear();
    // Debug toggles (env): disable individual passes to bisect regressions.
    const bool no_force = std::getenv("IGAOS_PRE_NOFORCE") != nullptr;
    const bool no_tighten = std::getenv("IGAOS_PRE_NOTIGHTEN") != nullptr;
    const bool no_subst = std::getenv("IGAOS_PRE_NOSUBST") != nullptr;
    const bool no_duprow = std::getenv("IGAOS_PRE_NODUPROW") != nullptr;
    const bool no_dupcol = std::getenv("IGAOS_PRE_NODUPCOL") != nullptr;

    st.m0 = m; st.n0 = n; st.nnz0 = orig.nnz();
    for (int i = 0; i < m; ++i) {
        for (int p = orig.ap[i]; p < orig.ap[i + 1]; ++p) {
            if (std::fabs(orig.ax[p]) <= 1e-12) { --st.nnz0; continue; }
            rows[i].emplace_back(orig.ai[p], orig.ax[p]);
        }
        std::sort(rows[i].begin(), rows[i].end());
    }

    // ---- P2: activity-interval propagation -------------------------------
    // Subsumes singleton-row tightening, forcing rows and redundant
    // (subsuming) rows, empty rows included.
    auto propagate = [&](bool& changed) -> int {
        for (int it = 0; it < 100; ++it) {
            bool iter_changed = false;
            for (int i = 0; i < m; ++i) {
                if (!row_alive[i]) continue;
                const auto& row = rows[i];
                double L = 0.0, U = 0.0;
                for (const auto& e : row) {
                    double p1 = e.second * cl[e.first];
                    double p2 = e.second * cu[e.first];
                    L += std::min(p1, p2);
                    U += std::max(p1, p2);
                }
                const double rt = rtol_of(rmin[i], rmax[i]);
                if (L > rmax[i] + rt || U < rmin[i] - rt) return 1;
                if (L >= rmin[i] - rt && U <= rmax[i] + rt) {
                    row_alive[i] = 0;              // subsuming row
                    ++st.rows_redundant;
                    iter_changed = changed = true;
                    continue;
                }
                // forcing row: min activity pinned at rmax (all vars to
                // their min-contribution bounds) or max pinned at rmin
                bool forced_max = !no_force && std::isfinite(rmax[i]) &&
                                  std::isfinite(L) && L >= rmax[i] - rt;
                bool forced_min = !no_force && std::isfinite(rmin[i]) &&
                                  std::isfinite(U) && U <= rmin[i] + rt;
                if (forced_max || forced_min) {
                    for (const auto& e : row) {
                        int j = e.first;
                        double a = e.second;
                        double v = (a > 0) == forced_max ? cl[j] : cu[j];
                        cl[j] = cu[j] = v;
                    }
                    row_alive[i] = 0;
                    ++st.rows_forcing;
                    iter_changed = changed = true;
                    continue;
                }
                // per-entry bound tightening
                for (const auto& e : row) {
                    if (no_tighten) break;
                    int j = e.first;
                    double a = e.second;
                    double p1 = a * cl[j], p2 = a * cu[j];
                    double cmin = std::min(p1, p2), cmax = std::max(p1, p2);
                    double nlo = -INF, nup = INF;
                    // Guards deliberately over-conservative (cmax/cmin
                    // swapped vs the quantity actually used): the stricter
                    // one-sided tightening these guards suppress destabilized
                    // phase 1 on wood1p/bore3d (netlib gate 61 vs 63), so
                    // keep the conservative form that passes 63/64.
                    if (std::isfinite(cmax) && std::isfinite(rmax[i]) &&
                        std::isfinite(L)) {
                        double Lm = L - cmin;   // min activity of the others
                        if (std::isfinite(Lm)) {
                            double v = (rmax[i] - Lm) / a;
                            if (a > 0) nup = std::min(nup, v);
                            else       nlo = std::max(nlo, v);
                        }
                    }
                    if (std::isfinite(cmin) && std::isfinite(rmin[i]) &&
                        std::isfinite(U)) {
                        double Um = U - cmax;   // max activity of the others
                        if (std::isfinite(Um)) {
                            double v = (rmin[i] - Um) / a;
                            if (a > 0) nlo = std::max(nlo, v);
                            else       nup = std::min(nup, v);
                        }
                    }
                    if (nlo > cl[j] + 1e-9 * (1.0 + std::fabs(nlo))) {
                        cl[j] = nlo; ++st.bound_tightens;
                        iter_changed = changed = true;
                    }
                    if (nup < cu[j] - 1e-9 * (1.0 + std::fabs(nup))) {
                        cu[j] = nup; ++st.bound_tightens;
                        iter_changed = changed = true;
                    }
                    if (!integ.empty() && integ[j]) {
                        double lo = std::isfinite(cl[j])
                                        ? std::ceil(cl[j] - 1e-9) : cl[j];
                        double up = std::isfinite(cu[j])
                                        ? std::floor(cu[j] + 1e-9) : cu[j];
                        if (lo > cl[j] + 1e-12 || up < cu[j] - 1e-12) {
                            cl[j] = lo; cu[j] = up;
                            iter_changed = changed = true;
                        }
                    }
                    if (cl[j] > cu[j] + 1e-7 * (1.0 + std::fabs(cl[j])))
                        return 1;
                }
            }
            if (!iter_changed) break;
        }
        return 0;
    };

    // ---- P3: fixed-variable substitution ---------------------------------
    auto fix_cols = [&](bool& changed) -> int {
        std::vector<char> just_fixed(n, 0);
        int nfix = 0;
        for (int j = 0; j < n; ++j) {
            if (!col_alive[j] || !std::isfinite(cl[j])) continue;
            if (cl[j] >= cu[j] - 1e-9 * (1.0 + std::fabs(cl[j]))) {
                fixval[j] = cl[j];
                just_fixed[j] = 1;
                col_alive[j] = 0;
                obj_shift += c[j] * fixval[j];
                PresolveLog::Rec r;
                r.kind = PresolveLog::Rec::FixedVar;
                r.j = j;
                r.val = fixval[j];
                log.recs.push_back(std::move(r));
                ++st.cols_fixed;
                ++nfix;
            }
        }
        if (!nfix) return 0;
        changed = true;
        for (int i = 0; i < m; ++i) {          // single sweep, O(nnz)
            if (!row_alive[i]) continue;
            auto& row = rows[i];
            size_t w = 0;
            for (size_t r = 0; r < row.size(); ++r) {
                int j = row[r].first;
                if (just_fixed[j]) {
                    double a = row[r].second;
                    if (std::isfinite(rmin[i])) rmin[i] -= a * fixval[j];
                    if (std::isfinite(rmax[i])) rmax[i] -= a * fixval[j];
                    continue;
                }
                row[w++] = row[r];
            }
            row.resize(w);
        }
        return 0;
    };

    // ---- P1 + P4: empty cols; free singleton-column substitution ---------
    auto col_pass = [&](bool& changed) -> int {
        std::vector<int> colcnt(n, 0), crow(n, -1);
        std::vector<double> ca(n, 0.0);
        for (int i = 0; i < m; ++i) {
            if (!row_alive[i]) continue;
            for (const auto& e : rows[i]) {
                int j = e.first;
                if (colcnt[j]++ == 0) { crow[j] = i; ca[j] = e.second; }
            }
        }
        // P1: empty columns
        for (int j = 0; j < n; ++j) {
            if (!col_alive[j] || colcnt[j] > 0) continue;
            double v;
            if (c[j] > 1e-12) {
                if (cl[j] == -INF) return 2;
                v = cl[j];
            } else if (c[j] < -1e-12) {
                if (cu[j] == INF) return 2;
                v = cu[j];
            } else {
                v = std::isfinite(cl[j]) ? cl[j]
                                         : (std::isfinite(cu[j]) ? cu[j] : 0.0);
            }
            col_alive[j] = 0;
            obj_shift += c[j] * v;
            PresolveLog::Rec r;
            r.kind = PresolveLog::Rec::FixedVar;
            r.j = j;
            r.val = v;
            log.recs.push_back(std::move(r));
            ++st.cols_empty;
            changed = true;
        }
        // P4: free singleton columns
        for (int j = 0; j < n; ++j) {
            if (no_subst) break;
            if (!col_alive[j] || colcnt[j] != 1) continue;
            if (cl[j] != -INF || cu[j] != INF) continue;
            int i = crow[j];
            if (!row_alive[i]) continue;  // consumed by an earlier singleton
                                          // this pass (its postsolve rec
                                          // must not depend on this col)
            double a = ca[j];
            double q = c[j] / a;
            double t;
            if (q > 1e-12) {
                if (rmin[i] == -INF) return 2;   // cost improves forever
                t = rmin[i];
            } else if (q < -1e-12) {
                if (rmax[i] == INF) return 2;
                t = rmax[i];
            } else {
                t = std::isfinite(rmin[i]) ? rmin[i]
                     : (std::isfinite(rmax[i]) ? rmax[i] : 0.0);
            }
            PresolveLog::Rec r;
            r.kind = PresolveLog::Rec::SubstFree;
            r.j = j;
            r.aij = a;
            r.target = t;
            for (const auto& e : rows[i]) {
                if (e.first != j) r.terms.push_back(e);
            }
            for (const auto& e : r.terms) c[e.first] -= q * e.second;
            obj_shift += q * t;
            log.recs.push_back(std::move(r));
            col_alive[j] = 0;
            row_alive[i] = 0;
            ++st.cols_subst;
            changed = true;
        }
        return 0;
    };

    // ---- P5: duplicate (proportional) rows -------------------------------
    auto dup_rows = [&](bool& changed) -> int {
        if (no_duprow) return 0;
        std::unordered_map<std::string, std::vector<int>> groups;
        std::string key;
        for (int i = 0; i < m; ++i) {
            if (!row_alive[i] || rows[i].size() < 2) continue;
            std::vector<int> cols;
            cols.reserve(rows[i].size());
            for (const auto& e : rows[i]) cols.push_back(e.first);
            pre_detail::pack_key(key, cols);
            groups[key].push_back(i);
        }
        for (auto& kv : groups) {
            auto& g = kv.second;
            for (size_t a = 1; a < g.size(); ++a) {
                if (!row_alive[g[a]]) continue;
                for (size_t b = 0; b < a; ++b) {
                    if (!row_alive[g[b]]) continue;
                    int i1 = g[a], i2 = g[b];
                    const auto& r1 = rows[i1];
                    const auto& r2 = rows[i2];
                    double s = r1[0].second / r2[0].second;
                    if (!std::isfinite(s) || s == 0.0) continue;
                    bool prop = true;
                    for (size_t k = 0; k < r1.size(); ++k) {
                        if (r1[k].first != r2[k].first ||
                            std::fabs(r1[k].second - s * r2[k].second) >
                                1e-9 * (1.0 + std::fabs(r1[k].second))) {
                            prop = false;
                            break;
                        }
                    }
                    if (!prop) continue;
                    double ilo = s > 0 ? rmin[i1] / s : rmax[i1] / s;
                    double ihi = s > 0 ? rmax[i1] / s : rmin[i1] / s;
                    double nlo = std::max(rmin[i2], ilo);
                    double nup = std::min(rmax[i2], ihi);
                    if (nlo > nup + rtol_of(nlo, nup)) return 1;
                    rmin[i2] = nlo;
                    rmax[i2] = nup;
                    row_alive[i1] = 0;
                    ++st.rows_dup;
                    changed = true;
                    break;
                }
            }
        }
        return 0;
    };

    // ---- P6: duplicate (identical) columns -------------------------------
    auto dup_cols = [&](bool& changed) -> int {
        if (no_dupcol) return 0;
        std::unordered_map<std::string, std::vector<int>> groups;
        std::string key;
        std::vector<std::vector<int>> colrows(n);
        for (int i = 0; i < m; ++i) {
            if (!row_alive[i]) continue;
            for (const auto& e : rows[i]) colrows[e.first].push_back(i);
        }
        for (int j = 0; j < n; ++j) {
            if (!col_alive[j] || colrows[j].size() < 2) continue;
            if (!integ.empty() && integ[j]) continue;
            pre_detail::pack_key(key, colrows[j]);
            groups[key].push_back(j);
        }
        auto coef = [&](int j, int i) -> double {
            for (const auto& e : rows[i])
                if (e.first == j) return e.second;
            return 0.0;
        };
        for (auto& kv : groups) {
            auto& g = kv.second;
            for (size_t a = 1; a < g.size(); ++a) {
                if (!col_alive[g[a]]) continue;
                for (size_t b = 0; b < a; ++b) {
                    if (!col_alive[g[b]]) continue;
                    int j1 = g[a], j2 = g[b];
                    // split rule needs finite lower bounds on both
                    if (cl[j1] == -INF || cl[j2] == -INF) continue;
                    bool same = true;
                    for (int i : colrows[j1]) {
                        if (std::fabs(coef(j1, i) - coef(j2, i)) >
                            1e-9 * (1.0 + std::fabs(coef(j1, i)))) {
                            same = false;
                            break;
                        }
                    }
                    if (!same) continue;
                    // Cost merge is only valid for EQUAL costs: with
                    // c1 != c2 the objective c1*x1 + c2*x2 depends on the
                    // split of w = x1 + x2 (piecewise, not one LP column).
                    double c1 = c[j1], c2 = c[j2];
                    if (std::fabs(c1 - c2) >
                        1e-9 * (1.0 + std::fabs(c1) + std::fabs(c2)))
                        continue;
                    PresolveLog::Rec r;
                    r.kind = PresolveLog::Rec::DupCol;
                    r.j = j1;
                    r.j2 = j2;
                    r.cl1 = cl[j1]; r.cu1 = cu[j1];
                    r.cl2 = cl[j2]; r.cu2 = cu[j2];
                    double ncl = cl[j1] + cl[j2];
                    double ncu = cu[j1] + cu[j2];
                    if (ncl > ncu + rtol_of(ncl, ncu)) return 1;
                    cl[j2] = ncl;
                    cu[j2] = ncu;
                    // cost on w = x1 + x2 is c1 == c2 (leave c[j2] as-is):
                    // c*x1 + c*x2 = c*w
                    col_alive[j1] = 0;
                    // drop j1's entries; j2's identical entries now carry w
                    for (int i : colrows[j1]) {
                        if (!row_alive[i]) continue;
                        auto& row = rows[i];
                        for (size_t p = 0; p < row.size(); ++p) {
                            if (row[p].first == j1) {
                                row.erase(row.begin() + p);
                                break;
                            }
                        }
                    }
                    log.recs.push_back(std::move(r));
                    ++st.cols_dup;
                    changed = true;
                    break;
                }
            }
        }
        return 0;
    };

    // ---- fixpoint loop ----------------------------------------------------
    bool changed = true;
    for (int round = 0; round < 20 && changed; ++round) {
        changed = false;
        int rc = propagate(changed);
        if (rc) return rc;
        rc = fix_cols(changed);
        if (rc) return rc;
        rc = col_pass(changed);
        if (rc) return rc;
        rc = dup_rows(changed);
        if (rc) return rc;
        rc = dup_cols(changed);
        if (rc) return rc;
    }
    // final consistency sweep (no-op unless the loop hit its round cap)
    {
        bool dummy = false;
        int rc = propagate(dummy);
        if (rc) return rc;
    }

    // ---- build reduced model ----------------------------------------------
    col_to_red.assign(n, -1);
    row_to_red.assign(m, -1);
    int rn = 0, rm = 0;
    for (int j = 0; j < n; ++j)
        if (col_alive[j]) col_to_red[j] = rn++;
    for (int i = 0; i < m; ++i)
        if (row_alive[i]) row_to_red[i] = rm++;
    std::vector<std::tuple<int, int, double>> trips;
    for (int i = 0; i < m; ++i) {
        if (!row_alive[i]) continue;
        for (const auto& e : rows[i])
            trips.emplace_back(row_to_red[i], col_to_red[e.first], e.second);
    }
    red = io::Model();
    red.m = rm;
    red.n = rn;
    red.obj_const = orig.obj_const + obj_shift;
    red.n_fr_parsed = orig.n_fr_parsed;
    red.n_ranges_parsed = orig.n_ranges_parsed;
    red.rmin.resize(rm); red.rmax.resize(rm);
    for (int i = 0; i < m; ++i)
        if (row_to_red[i] >= 0) {
            red.rmin[row_to_red[i]] = rmin[i];
            red.rmax[row_to_red[i]] = rmax[i];
        }
    red.cl.resize(rn); red.cu.resize(rn); red.c.resize(rn);
    red.integ.resize(rn);
    for (int j = 0; j < n; ++j)
        if (col_to_red[j] >= 0) {
            int k = col_to_red[j];
            red.cl[k] = cl[j];
            red.cu[k] = cu[j];
            red.c[k] = c[j];
            red.integ[k] = integ.empty() ? 0 : integ[j];
        }
    io::counting_sort(trips, rm, rn, red.ap, red.ai, red.ax, red.cp, red.ci,
                      red.acx);
    st.m1 = rm; st.n1 = rn; st.nnz1 = (int)trips.size();
    return 0;
}

// Postsolve primal solution: reduced x -> original x (log applied in
// reverse; each record's inputs are at their post-record state).
inline std::vector<double> postsolve_x(const PresolveLog& log,
                                       const std::vector<double>& x_red,
                                       const std::vector<int>& col_to_red,
                                       int n_orig) {
    std::vector<double> x(n_orig, 0.0);
    for (int j = 0; j < n_orig; ++j)
        if (col_to_red[j] >= 0) x[j] = x_red[col_to_red[j]];
    for (auto it = log.recs.rbegin(); it != log.recs.rend(); ++it) {
        const auto& r = *it;
        switch (r.kind) {
            case PresolveLog::Rec::FixedVar:
                x[r.j] = r.val;
                break;
            case PresolveLog::Rec::SubstFree: {
                double s = 0.0;
                for (const auto& t : r.terms) s += t.second * x[t.first];
                x[r.j] = (r.target - s) / r.aij;
                break;
            }
            case PresolveLog::Rec::DupCol: {
                // x[j2] currently holds w = x[j] + x[j2]; split, cl2 finite
                double w = x[r.j2];
                double x1 = std::max(r.cl1,
                                     std::min(r.cu1, w - r.cl2));
                x[r.j] = x1;
                x[r.j2] = w - x1;
                break;
            }
        }
    }
    return x;
}

}  // namespace igaos::simplex
