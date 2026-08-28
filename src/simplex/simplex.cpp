#include "simplex.hpp"

#include "dense_lu.hpp"
#include "presolve.hpp"
#include "scaling.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <limits>
#include <random>
#include <vector>

namespace igaos::simplex {
namespace {

using std::vector;

constexpr double INF = std::numeric_limits<double>::infinity();

enum VarKind : unsigned char { KIND_X = 0, KIND_SLACK = 1, KIND_ART = 2 };
enum NBStatus : unsigned char { AT_LOWER = 0, AT_UPPER = 1, BASIC = 2,
                                FREE_NB = 3 };

struct Engine {
    int m = 0;
    int nx = 0;
    long total = 0;
    vector<double> lo, up, cost;
    vector<unsigned char> kind;
    vector<int> art_row;
    vector<double> art_sigma;
    vector<int> cp, ci, ap, ai;
    vector<double> cv, av;
    Scaling sc;
    vector<double> xs_cost;
    double cmax_x = 0.0;

    void col_dense(int j, vector<double>& out) const {
        out.assign(m, 0.0);
        switch (kind[j]) {
            case KIND_X:
                for (int p = cp[j]; p < cp[j + 1]; ++p) out[ci[p]] = cv[p];
                break;
            case KIND_SLACK:
                out[j - nx] = -1.0;
                break;
            case KIND_ART:
                out[art_row[j - nx - m]] = art_sigma[j - nx - m];
                break;
        }
    }

    double reduced_cost_part(int j, const vector<double>& y) const {
        switch (kind[j]) {
            case KIND_X: {
                double d = cost[j];
                for (int p = cp[j]; p < cp[j + 1]; ++p)
                    d -= y[ci[p]] * cv[p];
                return d;
            }
            case KIND_SLACK:
                return cost[j] + y[j - nx];
            case KIND_ART:
                return cost[j] -
                       y[art_row[j - nx - m]] * art_sigma[j - nx - m];
        }
        return 0.0;
    }

    void accumulate_row_rhs(int j, double vj, vector<double>& rhs) const {
        if (vj == 0.0) return;
        switch (kind[j]) {
            case KIND_X:
                for (int p = cp[j]; p < cp[j + 1]; ++p)
                    rhs[ci[p]] -= vj * cv[p];
                break;
            case KIND_SLACK:
                rhs[j - nx] += vj;
                break;
            case KIND_ART:
                rhs[art_row[j - nx - m]] += vj * art_sigma[j - nx - m];
                break;
        }
    }

    void basis_col(int bj, vector<double>& out) const {
        out.assign(m, 0.0);
        switch (kind[bj]) {
            case KIND_X:
                for (int p = cp[bj]; p < cp[bj + 1]; ++p)
                    out[ci[p]] = cv[p];
                break;
            case KIND_SLACK:
                out[bj - nx] = -1.0;
                break;
            case KIND_ART:
                out[art_row[bj - nx - m]] = art_sigma[bj - nx - m];
                break;
        }
    }
};

}  // namespace

Solution solve(const io::Model& model, const Options& opt,
               const WarmStart* warm, WarmStart* warm_out) {
    Solution sol;
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             t0)
            .count();
    };

    io::Model work = model;
    PresolveInfo pinfo;
    {
        int nch = strengthen_bounds(work, pinfo, 3);
        if (opt.verbosity > 0)
            std::fprintf(stderr,
                         "[simplex] presolve: %d bound tightenings\n", nch);
    }

    Engine E;
    E.m = work.m;
    E.nx = work.n;

    geometric_mean_scaling(work.ap, work.ai, work.ax, work.cp, work.ci,
                           work.acx, work.m, work.n, 10, E.sc);
    power_of_two_snap(E.sc.row);
    power_of_two_snap(E.sc.col);

    E.cp = work.cp;
    E.ci = work.ci;
    E.cv.resize(work.nnz());
    for (int j = 0; j < work.n; ++j)
        for (int p = work.cp[j]; p < work.cp[j + 1]; ++p)
            E.cv[p] = work.acx[p] * E.sc.row[work.ci[p]] * E.sc.col[j];
    E.ap = work.ap;
    E.ai = work.ai;
    E.av.resize(work.nnz());
    for (int i = 0; i < work.m; ++i)
        for (int p = work.ap[i]; p < work.ap[i + 1]; ++p)
            E.av[p] = work.ax[p] * E.sc.col[work.ai[p]] * E.sc.row[i];

    E.kind.assign((size_t)work.n + work.m, KIND_SLACK);
    for (int j = 0; j < work.n; ++j) E.kind[j] = KIND_X;
    E.lo.resize((size_t)work.n + work.m);
    E.up.resize((size_t)work.n + work.m);
    E.cost.assign((size_t)work.n + work.m, 0.0);
    E.xs_cost.resize(work.n);
    for (int j = 0; j < work.n; ++j) {
        E.lo[j] = std::isfinite(work.cl[j])
                      ? work.cl[j] / E.sc.col[j]
                      : -INF;
        E.up[j] = std::isfinite(work.cu[j])
                      ? work.cu[j] / E.sc.col[j]
                      : INF;
        E.xs_cost[j] = work.c[j] * E.sc.col[j];
        E.cmax_x = std::max(E.cmax_x, std::fabs(E.xs_cost[j]));
    }
    for (int i = 0; i < work.m; ++i) {
        E.lo[(size_t)work.n + i] =
            std::isfinite(work.rmin[i]) ? work.rmin[i] * E.sc.row[i] : -INF;
        E.up[(size_t)work.n + i] =
            std::isfinite(work.rmax[i]) ? work.rmax[i] * E.sc.row[i] : INF;
        E.cost[(size_t)work.n + i] = 0.0;
    }

    vector<int> basis(work.m, -1);
    vector<unsigned char> st((size_t)work.n + work.m, AT_LOWER);
    vector<double> z((size_t)work.n + work.m, 0.0);

    // Warm-start restore: reuse a parent solve's basis/states for a model
    // with the same matrix (B&B child, bounds-only difference). Basics
    // violating the child bounds are snapped to the violated bound and
    // their slots become artificials — the existing elastic phase 1 then
    // restores feasibility from a nearly-good basis.
    bool warm_ok = warm != nullptr &&
                   (int)warm->basis.size() == work.m &&
                   (int)warm->nb_state.size() == work.n + work.m;
    if (warm_ok) {
        for (int i = 0; i < work.m; ++i) {
            long v = warm->basis[i];
            basis[i] = (v >= 0 && v < (long)(work.n + work.m)) ? (int)v : -1;
        }
        // nonbasic states from snapshot; z recomputed from bounds so the
        // child's tighter bounds are honored
        for (long j = 0; j < (long)(work.n + work.m); ++j) {
            bool in_basis = false;
            for (int i = 0; i < work.m && !in_basis; ++i)
                in_basis = basis[i] == (int)j;
            if (in_basis) { st[j] = BASIC; continue; }
            unsigned char s = warm->nb_state[j];
            st[j] = s;
            if (s == AT_LOWER && std::isfinite(E.lo[j])) z[j] = E.lo[j];
            else if (s == AT_UPPER && std::isfinite(E.up[j])) z[j] = E.up[j];
            else if (!std::isfinite(E.lo[j]) && std::isfinite(E.up[j]))
                { st[j] = AT_UPPER; z[j] = E.up[j]; }
            else if (std::isfinite(E.lo[j])) { st[j] = AT_LOWER; z[j] = E.lo[j]; }
            else { st[j] = FREE_NB; z[j] = 0.0; }
        }
        // evaluate basic values via a factorization of the restored basis
        DenseLU wlu;
        {
            vector<double> Bm((size_t)E.m * E.m, 0.0), col;
            for (int c = 0; c < E.m; ++c) {
                if (basis[c] < 0) { Bm[(size_t)c * E.m + c] = 1.0; continue; }
                E.basis_col(basis[c], col);
                for (int r = 0; r < E.m; ++r)
                    if (col[r] != 0.0) Bm[(size_t)r * E.m + c] = col[r];
            }
            wlu.factor(std::move(Bm));
        }
        if (!wlu.ok) {
            warm_ok = false;  // singular restored basis — fall back cold
        } else {
            vector<double> wzb(work.m, 0.0);
            {
                vector<double> rhs(work.m, 0.0);
                for (long j = 0; j < (long)(work.n + work.m); ++j) {
                    if (st[j] == BASIC || z[j] == 0.0) continue;
                    E.accumulate_row_rhs((int)j, z[j], rhs);
                }
                wlu.solve(rhs, wzb);
            }
            for (int i = 0; i < work.m; ++i) {
                int v = basis[i];
                if (v < 0) continue;
                double val = wzb[i];
                if (val < E.lo[v]) {
                    st[v] = AT_LOWER; z[v] = E.lo[v]; basis[i] = -1;
                } else if (val > E.up[v]) {
                    st[v] = AT_UPPER; z[v] = E.up[v]; basis[i] = -1;
                } else {
                    z[v] = val;
                }
            }
            // re-validate the augmented basis: snapped slots will become
            // art columns (+e_i), and replacing a basis column by e_i can
            // be singular even when the original was not
            {
                vector<double> Bm((size_t)E.m * E.m, 0.0), col;
                for (int c = 0; c < E.m; ++c) {
                    if (basis[c] < 0) { Bm[(size_t)c * E.m + c] = 1.0; continue; }
                    E.basis_col(basis[c], col);
                    for (int r = 0; r < E.m; ++r)
                        if (col[r] != 0.0) Bm[(size_t)r * E.m + c] = col[r];
                }
                DenseLU chk;
                chk.factor(std::move(Bm));
                if (!chk.ok) warm_ok = false;
            }
        }
    }
    if (!warm_ok) {
        for (int j = 0; j < work.n; ++j) {
            if (std::isfinite(E.lo[j])) { st[j] = AT_LOWER; z[j] = E.lo[j]; }
            else if (std::isfinite(E.up[j])) { st[j] = AT_UPPER; z[j] = E.up[j]; }
            else { st[j] = FREE_NB; z[j] = 0.0; }
        }
        vector<double> ax_val(work.m, 0.0);
        for (int j = 0; j < work.n; ++j)
            if (z[j] != 0.0)
                for (int p = E.cp[j]; p < E.cp[j + 1]; ++p)
                    ax_val[E.ci[p]] += E.cv[p] * z[j];

        for (int i = 0; i < work.m; ++i) {
            int sv = work.n + i;
            if (ax_val[i] < E.lo[sv]) { st[sv] = AT_LOWER; z[sv] = E.lo[sv]; }
            else if (ax_val[i] > E.up[sv]) { st[sv] = AT_UPPER; z[sv] = E.up[sv]; }
            else { st[sv] = BASIC; basis[i] = sv; z[sv] = ax_val[i]; }
        }
    }
    vector<double> ax_val(work.m, 0.0);
    if (!warm_ok) {
        for (int j = 0; j < work.n; ++j)
            if (z[j] != 0.0)
                for (int p = E.cp[j]; p < E.cp[j + 1]; ++p)
                    ax_val[E.ci[p]] += E.cv[p] * z[j];
    }
    for (int i = 0; i < work.m; ++i) {
        if (basis[i] >= 0) continue;
        // warm path: art value is a placeholder, recomputed by the
        // refresh right after factorization
        double rho = warm_ok ? 0.0 : ax_val[i] - z[work.n + i];
        E.art_row.push_back(i);
        E.art_sigma.push_back(rho > 0 ? -1.0 : 1.0);
        int av_idx = (int)E.kind.size();
        E.kind.push_back(KIND_ART);
        E.lo.push_back(0.0);
        E.up.push_back(INF);
        E.cost.push_back(1.0);
        st.push_back(BASIC);
        z.push_back(std::fabs(rho));
        basis[i] = av_idx;
    }
    E.total = (long)E.kind.size();

    auto bpos_of = [&](long j) -> int {
        for (int i = 0; i < E.m; ++i)
            if (basis[i] == (int)j) return i;
        return -1;
    };

    EtaFile lu;
    bool basis_dirty = true;  // eta file must be rebuilt from scratch
    vector<double> zb(E.m, 0.0);
    auto build_and_factor = [&]() {
        vector<double> Bm((size_t)E.m * E.m, 0.0), col;
        for (int c = 0; c < E.m; ++c) {
            E.basis_col(basis[c], col);
            for (int r = 0; r < E.m; ++r)
                if (col[r] != 0.0) Bm[(size_t)r * E.m + c] = col[r];
        }
        lu.factor(std::move(Bm));
        basis_dirty = false;  // eta file now matches `basis`
        return lu.ok;
    };
    auto refresh_values = [&]() {
        vector<double> rhs(E.m, 0.0);
        for (long j = 0; j < E.total; ++j) {
            if (st[j] == BASIC || z[j] == 0.0) continue;
            E.accumulate_row_rhs((int)j, z[j], rhs);
        }
        lu.solve(rhs, zb);
        for (int i = 0; i < E.m; ++i) z[basis[i]] = zb[i];
    };

    bool inv_logged = false;

    std::mt19937 rng(opt.seed ? (unsigned)opt.seed : 20260825u);
    long iter = 0;
    int degen_streak = 0;
    bool bland = false;
    bool perturbed = false;
    // Columns blocked this basis: their only ratio-test blockers sit below
    // the Harris stability threshold, so pivoting on them is numerically
    // unsafe (observed: singular basis on perold). Skipped in pricing,
    // cleared after every pivot/flip; if pricing exhausts anyway, the
    // unsafe pivot is accepted as a last resort.
    std::vector<char> skip_col(E.total, 0);
    bool accept_unsafe = false;
    const int REFACTOR_EVERY = 64;

    auto run_simplex = [&](bool phase1) -> int {
        std::fill(skip_col.begin(), skip_col.end(), 0);
        accept_unsafe = false;
        basis_dirty = true;  // force a factorization at phase entry
        auto verify_pivot_invariant = [&](const char* tag, int enter,
                                          int lpos, long lvar, double th,
                                          double dr) {
            if (inv_logged) return;
            vector<int> pos(E.total, -1);
            bool ghost = false;
            for (int i = 0; i < E.m; ++i) {
                if (basis[i] < 0 || basis[i] >= (int)E.total) {
                    ghost = true;
                    continue;
                }
                if (pos[basis[i]] >= 0) ghost = true;
                pos[basis[i]] = i;
            }
            vector<double> act(E.m, 0.0);
            double mag = 1.0;
            for (long j = 0; j < E.total; ++j) {
                long pj = (st[j] == BASIC) ? pos[j] : -1;
                double vj = (pj >= 0) ? zb[pj] : z[j];
                if (st[j] == BASIC && pj < 0) ghost = true;
                if (vj == 0.0) continue;
                switch (E.kind[j]) {
                    case KIND_X:
                        for (int p = E.cp[j]; p < E.cp[j + 1]; ++p) {
                            double t = E.cv[p] * vj;
                            act[E.ci[p]] += t;
                            mag = std::max(mag, std::fabs(t));
                        }
                        break;
                    case KIND_SLACK: {
                        double t = -vj;
                        act[j - E.nx] += t;
                        mag = std::max(mag, std::fabs(t));
                        break;
                    }
                    case KIND_ART: {
                        int k = (int)(j - E.nx - E.m);
                        double t = E.art_sigma[k] * vj;
                        act[E.art_row[k]] += t;
                        mag = std::max(mag, std::fabs(t));
                        break;
                    }
                }
            }
            double resid = ghost ? INF : 0.0;
            for (int r = 0; r < E.m; ++r)
                resid = std::max(resid, std::fabs(act[r]));
            if (resid <= 1e-6 * mag) return;
            inv_logged = true;
            std::fprintf(stderr,
                         "[simplex] PIVOT INVARIANT BREAK (%s) it=%ld "
                         "ph1=%d |Mz|=%.6g enter=%ld leave=%ld@%d "
                         "theta=%.10g dir=%.6g zb[lp]=%.10g "
                         "z[enter]=%.10g\n",
                         tag, iter, (int)phase1, resid, (long)enter, lvar,
                         lpos, th, dr, lpos >= 0 ? zb[lpos] : 0.0,
                         z[enter]);
        };
        bool retry_done = false;
        // pivot-recovery snapshot: last (basis, st) before the most recent
        // basis-changing pivot, plus the entering column — used to roll
        // back a pivot whose basis turns out singular at refactorization
        vector<int> rec_basis;
        vector<unsigned char> rec_st;
        int rec_enter = -1;
        for (; iter <= opt.max_iterations; ++iter) {
            {
                for (long j = 0; j < E.total; ++j) {
                    if (st[j] == BASIC) continue;
                    if (st[j] == AT_LOWER && std::isfinite(E.lo[j]))
                        z[j] = E.lo[j];
                    else if (st[j] == AT_UPPER && std::isfinite(E.up[j]))
                        z[j] = E.up[j];
                    else if (std::isfinite(E.lo[j])) { st[j] = AT_LOWER; z[j] = E.lo[j]; }
                    else if (std::isfinite(E.up[j])) { st[j] = AT_UPPER; z[j] = E.up[j]; }
                    else st[j] = FREE_NB;
                }

                if (basis_dirty || lu.needs_refactor()) {
                    if (!build_and_factor()) {
                        if (rec_enter >= 0) {
                            // the last pivot singularized the basis —
                            // roll it back and reject that entering column
                            basis = rec_basis;
                            st = rec_st;
                            skip_col[rec_enter] = 1;
                            rec_enter = -1;
                            basis_dirty = true;  // etas are stale
                            continue;
                        }
                        sol.status = Status::Error;
                        sol.message = "basis factorization failed";
                        return 5;
                    }
                    basis_dirty = false;
                }
                refresh_values();
                refresh_values();
            }
            if (opt.verbosity > 2) {
                vector<int> pos(E.total, -1);
                for (int i = 0; i < E.m; ++i)
                    if (basis[i] < (int)E.total && basis[i] >= 0)
                        pos[basis[i]] = i;
                vector<double> act(E.m, 0.0);
                vector<double> col;
                double zbad = 0.0;
                long zbad_j = -1;
                for (long j = 0; j < E.total; ++j) {
                    double vj;
                    if (st[j] == BASIC) {
                        if (pos[j] < 0) {
                            zbad = std::numeric_limits<double>::infinity();
                            zbad_j = j;
                            continue;
                        }
                        vj = zb[pos[j]];
                    } else {
                        vj = z[j];
                        if (st[j] == AT_LOWER &&
                            std::fabs(vj - E.lo[j]) > 1e-6 &&
                            std::isfinite(E.lo[j]))
                            { zbad = std::fabs(vj - E.lo[j]); zbad_j = j; }
                        if (st[j] == AT_UPPER &&
                            std::fabs(vj - E.up[j]) > 1e-6 &&
                            std::isfinite(E.up[j]))
                            { zbad = std::fabs(vj - E.up[j]); zbad_j = j; }
                    }
                    if (vj != 0.0) {
                        E.col_dense((int)j, col);
                        for (int r = 0; r < E.m; ++r) act[r] += col[r] * vj;
                    }
                }
                double res = 0.0;
                int res_row = -1;
                for (int r = 0; r < E.m; ++r)
                    if (std::fabs(act[r]) > res) {
                        res = std::fabs(act[r]);
                        res_row = r;
                    }
                if (res > 1e-7 || zbad > 1e-6 || zbad_j >= 0)
                    std::fprintf(stderr,
                                 "[simplex] PROBE it=%ld ph1=%d "
                                 "res=%.3g@row%d zbad=%.3g@var%ld\n",
                                 iter, (int)phase1, res, res_row, zbad,
                                 zbad_j);
            }

            if (phase1) {
                double asum = 0.0;
                double bviol = 0.0;
                for (int i = 0; i < E.m; ++i) {
                    if (E.kind[basis[i]] == KIND_ART)
                        asum += std::fabs(zb[i]);
                    bviol = std::max(
                        bviol,
                        std::max(E.lo[basis[i]] - zb[i],
                                 zb[i] - E.up[basis[i]]));
                }
                if (asum <= 1e-8 * (1.0 + E.cmax_x) &&
                    bviol <= 1e-7 * (1.0 + E.cmax_x))
                    return 0;
            }

            vector<double> cb(E.m, 0.0);
            for (int i = 0; i < E.m; ++i) cb[i] = E.cost[basis[i]];
            vector<double> y;
            lu.solve_transpose(cb, y);

            int enter0 = -1;
            double best_score = opt.tolerance * 1e-2;
            double dir0 = 1.0;
            for (long j = 0; j < E.total; ++j) {
                if (st[j] == BASIC) continue;
                if (skip_col[j]) continue;
                if (std::isfinite(E.lo[j]) && std::isfinite(E.up[j]) &&
                    E.up[j] - E.lo[j] <= 1e-12)
                    continue;
                double dj = E.reduced_cost_part((int)j, y);
                if (st[j] == FREE_NB) {
                    double score = std::fabs(dj);
                    if (score > best_score) {
                        best_score = score;
                        enter0 = (int)j;
                        dir0 = dj > 0 ? -1.0 : 1.0;
                    }
                } else if (st[j] == AT_LOWER && dj < -best_score) {
                    best_score = -dj;
                    enter0 = (int)j;
                    dir0 = 1.0;
                } else if (st[j] == AT_UPPER && dj > best_score) {
                    best_score = dj;
                    enter0 = (int)j;
                    dir0 = -1.0;
                }
            }
            if (enter0 < 0 && bland) {
                for (long j = 0; j < E.total; ++j) {
                    if (st[j] == BASIC) continue;
                    if (skip_col[j]) continue;
                    if (std::isfinite(E.lo[j]) && std::isfinite(E.up[j]) &&
                        E.up[j] - E.lo[j] <= 1e-12)
                        continue;
                    double dj = E.reduced_cost_part((int)j, y);
                    if (st[j] == FREE_NB && std::fabs(dj) > 1e-9) {
                        enter0 = (int)j;
                        dir0 = dj > 0 ? -1.0 : 1.0;
                        break;
                    }
                    if (st[j] == AT_LOWER && dj < -1e-9) {
                        enter0 = (int)j;
                        dir0 = 1.0;
                        break;
                    }
                    if (st[j] == AT_UPPER && dj > 1e-9) {
                        enter0 = (int)j;
                        dir0 = -1.0;
                        break;
                    }
                }
            }
            if (enter0 < 0 && !skip_col.empty() && !accept_unsafe) {
                bool any = false;
                for (long j = 0; j < E.total && !any; ++j) any = skip_col[j];
                if (any) {
                    // every improving column is blocked below stability —
                    // accept an unsafe pivot rather than stall
                    accept_unsafe = true;
                    std::fill(skip_col.begin(), skip_col.end(), 0);
                    continue;
                }
            }
            if (enter0 < 0) return 0;

            vector<double> ae;
            E.col_dense(enter0, ae);
            vector<double> alpha;
            lu.solve(ae, alpha);
            double amax = 0.0;
            for (double a : alpha) amax = std::max(amax, std::fabs(a));

            double theta = INF;
            double theta_rel = INF;
            double theta_strict_min = INF;
            int leave_pos = -1;
            double leave_bound = 0.0;
            bool leave_to_upper = false;
            const double ALPHA_FLOOR = 1e-9;
            const double DELTA = 1e-8;
            const double STABILITY = 1e-7 * amax;
            if (amax > ALPHA_FLOOR) {
                for (int i = 0; i < E.m; ++i) {
                    if (std::fabs(alpha[i]) <= ALPHA_FLOOR) continue;
                    double rate = -alpha[i] * dir0;
                    bool to_upper = rate > 0;
                    double dist = to_upper ? (E.up[basis[i]] - zb[i])
                                           : (zb[i] - E.lo[basis[i]]);
                    if (!std::isfinite(dist)) continue;
                    if (dist < 0) dist = 0.0;
                    double bmag = std::fabs(to_upper ? E.up[basis[i]]
                                                     : E.lo[basis[i]]);
                    double relax = DELTA * (1.0 + bmag);
                    double t = (dist + relax) /
                               (std::fabs(rate) * (1.0 + DELTA));
                    if (t < theta_rel) theta_rel = t;
                    double ts = dist / std::fabs(rate);
                    if (ts < theta_strict_min) theta_strict_min = ts;
                }
                double best_abs = STABILITY;
                for (int i = 0; i < E.m; ++i) {
                    if (std::fabs(alpha[i]) < best_abs ||
                        std::fabs(alpha[i]) <= ALPHA_FLOOR)
                        continue;
                    double rate = -alpha[i] * dir0;
                    bool to_upper = rate > 0;
                    double dist = to_upper ? (E.up[basis[i]] - zb[i])
                                           : (zb[i] - E.lo[basis[i]]);
                    if (!std::isfinite(dist)) continue;
                    if (dist < 0) dist = 0.0;
                    double bmag = std::fabs(to_upper ? E.up[basis[i]]
                                                     : E.lo[basis[i]]);
                    double relax = DELTA * (1.0 + bmag);
                    double t_relaxed = (dist + relax) /
                                       (std::fabs(rate) * (1.0 + DELTA));
                    if (t_relaxed <= theta_rel * (1.0 + 1e-9) +
                                         1e-9) {
                        best_abs = std::fabs(alpha[i]);
                        leave_pos = i;
                        leave_to_upper = to_upper;
                        leave_bound = to_upper ? E.up[basis[i]]
                                               : E.lo[basis[i]];
                    }
                }
                if (bland) {
                    double tmin = INF;
                    for (int i = 0; i < E.m; ++i) {
                        if (std::fabs(alpha[i]) <= ALPHA_FLOOR) continue;
                        double rate2 = -alpha[i] * dir0;
                        double dist2 =
                            rate2 > 0 ? (E.up[basis[i]] - zb[i])
                                      : (zb[i] - E.lo[basis[i]]);
                        if (!std::isfinite(dist2)) continue;
                        if (dist2 < 0) dist2 = 0.0;
                        double ts = dist2 / std::fabs(rate2);
                        if (ts < tmin) tmin = ts;
                    }
                    if (std::isfinite(tmin)) {
                        double lim = tmin * (1.0 + 1e-12) + 1e-12;
                        long bestvar = std::numeric_limits<long>::max();
                        int bpos = -1;
                        bool bu = false;
                        double bb = 0.0;
                        for (int i = 0; i < E.m; ++i) {
                            if (std::fabs(alpha[i]) <= ALPHA_FLOOR)
                                continue;
                            double rate2 = -alpha[i] * dir0;
                            double dist2 =
                                rate2 > 0 ? (E.up[basis[i]] - zb[i])
                                          : (zb[i] - E.lo[basis[i]]);
                            if (!std::isfinite(dist2)) continue;
                            if (dist2 < 0) dist2 = 0.0;
                            double ts = dist2 / std::fabs(rate2);
                            if (ts <= lim && basis[i] < bestvar) {
                                bestvar = basis[i];
                                bpos = i;
                                bu = rate2 > 0;
                                bb = bu ? E.up[basis[i]] : E.lo[basis[i]];
                            }
                        }
                        if (bpos >= 0) {
                            theta = tmin;
                            leave_pos = bpos;
                            leave_to_upper = bu;
                            leave_bound = bb;
                        }
                    }
                }
                // Last-resort Harris fallback: the stability threshold
                // (1e-7 * amax) filtered out every blocker, but pass 1 saw
                // a finite ratio. Only safe to pivot here when the column
                // cannot be skipped (pricing exhausted) — otherwise the
                // column is rejected and pricing moves on.
                if (leave_pos < 0 && std::isfinite(theta_rel) &&
                    accept_unsafe) {
                    double fb_abs = ALPHA_FLOOR;
                    for (int i = 0; i < E.m; ++i) {
                        if (std::fabs(alpha[i]) <= ALPHA_FLOOR) continue;
                        double rate = -alpha[i] * dir0;
                        bool to_upper = rate > 0;
                        double dist = to_upper
                                          ? (E.up[basis[i]] - zb[i])
                                          : (zb[i] - E.lo[basis[i]]);
                        if (!std::isfinite(dist)) continue;
                        if (dist < 0) dist = 0.0;
                        double bmag = std::fabs(
                            to_upper ? E.up[basis[i]] : E.lo[basis[i]]);
                        double relax = DELTA * (1.0 + bmag);
                        double t_relaxed =
                            (dist + relax) /
                            (std::fabs(rate) * (1.0 + DELTA));
                        if (t_relaxed <= theta_rel * (1.0 + 1e-9) + 1e-9 &&
                            std::fabs(alpha[i]) >= fb_abs) {
                            fb_abs = std::fabs(alpha[i]);
                            leave_pos = i;
                            leave_to_upper = to_upper;
                            leave_bound = to_upper ? E.up[basis[i]]
                                                   : E.lo[basis[i]];
                        }
                    }
                }

                if (leave_pos >= 0) {
                    double rate =
                        -alpha[leave_pos] * dir0;
                    double dist = leave_to_upper
                                      ? (E.up[basis[leave_pos]] -
                                         zb[leave_pos])
                                      : (zb[leave_pos] -
                                         E.lo[basis[leave_pos]]);
                    if (dist < 0) dist = 0.0;
                    theta = dist / std::fabs(rate);
                } else {
                    theta = theta_rel;
                }
            }

            double range = (st[enter0] == FREE_NB)
                               ? INF
                               : (E.up[enter0] - E.lo[enter0]);

            if ((leave_pos < 0 || theta > range) && !std::isfinite(range)) {
                if (std::isfinite(theta_rel) && !accept_unsafe) {
                    // blocked (finite ratio exists) but no stable pivot —
                    // reject this column, reprice
                    skip_col[enter0] = 1;
                    continue;
                }
                if (!retry_done) {
                    retry_done = true;
                    continue;
                }
                if (opt.verbosity > 0) {
                    std::fprintf(stderr,
                                 "[RAY] enter=%d kind=%d lo=%.4g up=%.4g "
                                 "amax=%.3g\n",
                                 enter0, (int)E.kind[enter0], E.lo[enter0],
                                 E.up[enter0], amax);
                    for (int i = 0; i < E.m; ++i)
                        if (std::fabs(alpha[i]) > 1e-12)
                            std::fprintf(stderr,
                                         " p%d v%ld k%d a%.4g z%.4g "
                                         "l%.3g u%.3g\n",
                                         i, basis[i],
                                         (int)E.kind[basis[i]], alpha[i],
                                         zb[i], E.lo[basis[i]],
                                         E.up[basis[i]]);
                }
                sol.status = Status::Unbounded;
                sol.message = phase1 ? "phase 1 unbounded (anomaly)"
                                     : "unbounded: no blocking variable";
                return 2;
            }

            bool flip = std::isfinite(range) &&
                        (leave_pos < 0 ||
                         theta_strict_min >= range);
            if (flip) theta = range;

            for (int i = 0; i < E.m; ++i) zb[i] -= alpha[i] * theta * dir0;
            z[enter0] += theta * dir0;

            if (opt.verbosity > 2) {
                int bmin = INT_MAX, bmax = INT_MIN;
                for (int i = 0; i < E.m; ++i) {
                    bmin = std::min(bmin, basis[i]);
                    bmax = std::max(bmax, basis[i]);
                }
                std::fprintf(stderr,
                             "[simplex] it=%ld basis range [%d,%d] "
                             "total=%ld\n",
                             iter, bmin, bmax, E.total);
                const vector<double>& chk = zb;
                double worst_v = 0.0;
                int worst_i = -1;
                for (int i = 0; i < E.m; ++i) {
                    double v =
                        std::max(E.lo[basis[i]] - chk[i],
                                 chk[i] - E.up[basis[i]]);
                    v = std::max(v, 0.0);
                    if (v > worst_v) { worst_v = v; worst_i = i; }
                }
                bool bad_nb = false;
                int bad_j = -1;
                for (long j = 0; j < E.total && !bad_nb; ++j) {
                    if (st[j] == BASIC) continue;
                    if (st[j] == AT_LOWER &&
                        std::fabs(z[j] - E.lo[j]) > 1e-6)
                        bad_nb = true, bad_j = (int)j;
                    if (st[j] == AT_UPPER &&
                        std::fabs(z[j] - E.up[j]) > 1e-6)
                        bad_nb = true, bad_j = (int)j;
                }
                if ((worst_v > 1e-7 || bad_nb) && opt.verbosity > 2) {
                    std::fprintf(stderr,
                                 "[simplex] INVARIANT BREAK it=%ld "
                                 "worst_basic=%.3g@%d badnb=%d\n",
                                 iter, worst_v, worst_i, bad_j);
                    if (bad_j >= 0)
                        std::fprintf(stderr,
                                     "   var=%d kind=%d st=%u z=%.6g "
                                     "lo=%.6g up=%.6g\n",
                                     bad_j, (int)E.kind[bad_j],
                                     (unsigned)st[bad_j], z[bad_j],
                                     E.lo[bad_j], E.up[bad_j]);
                }
            }

            if (theta > 1e-9) {
                degen_streak = 0;
                bland = false;
            } else {
                ++degen_streak;
                if (degen_streak >= 20 && !bland) {
                    bland = true;
                    degen_streak = 0;
                }
                // deep stall: even Bland is not escaping — perturb costs
                // once (HiGHS-recipe flavor, robustness-strategy rung 4);
                // the post-solve cleanup pass restores them
                if (degen_streak >= 500 && !perturbed && !phase1) {
                    perturbed = true;
                    degen_streak = 0;
                    for (int j = 0; j < E.nx; ++j) {
                        if (E.cost[j] == 0.0) continue;
                        double mag = std::fabs(E.cost[j]);
                        double r = (double)rng() / (double)rng.max() - 0.5;
                        E.cost[j] += r * 1e-5 * mag;
                    }
                }
            }
            if (!skip_col.empty()) std::fill(skip_col.begin(),
                                             skip_col.end(), 0);
            accept_unsafe = false;

            if (flip) {
                st[enter0] = (dir0 > 0) ? AT_UPPER : AT_LOWER;
                z[enter0] = (dir0 > 0) ? E.up[enter0] : E.lo[enter0];
                verify_pivot_invariant("flip", enter0, -1, -1, theta, dir0);
                continue;
            }

            int leaving = basis[leave_pos];
            // snapshot for pivot rollback (see factor-failure recovery)
            rec_basis = basis;
            rec_st = st;
            rec_enter = enter0;
            z[leaving] = leave_bound;
            st[leaving] = leave_to_upper ? AT_UPPER : AT_LOWER;
            basis[leave_pos] = enter0;
            st[enter0] = BASIC;
            zb[leave_pos] = z[enter0];
            // product-form update: B_new = B_old · E(leave_pos, alpha)
            lu.update(leave_pos, alpha);
            verify_pivot_invariant("pivot", enter0, leave_pos, leaving,
                                   theta, dir0);

            if (elapsed() > opt.time_limit_s) {
                sol.status = Status::TimeLimit;
                sol.message = "wall-clock budget exhausted";
                return 3;
            }
        }
        sol.status = Status::IterationLimit;
        sol.message = "iteration cap reached";
        return 4;
    };

    if (getenv("IGAOS_DUMP")) {
        {
            FILE* g = fopen("/tmp/opencode/scaled_final.txt", "w");
            std::fprintf(g, "SCOLS\n");
            for (int j = 0; j < work.n; ++j)
                std::fprintf(g, "%d %.17g\n", j, E.sc.col[j]);
            std::fprintf(g, "SROWS\n");
            for (int i = 0; i < work.m; ++i)
                std::fprintf(g, "%d %.17g\n", i, E.sc.row[i]);
            std::fprintf(g, "FINAL\n");
            for (long j = 0; j < E.total; ++j)
                std::fprintf(g, "%ld %.17g %d\n", j,
                             st[j] == BASIC ? zb[bpos_of(j)] : z[j],
                             (int)st[j]);
            fclose(g);
        }
        FILE* f = fopen("/tmp/opencode/hav_model.txt", "w");
        std::fprintf(f, "m=%d n=%d nnz=%d\n", work.m, work.n,
                     work.nnz());
        for (int i = 0; i < work.m; ++i)
            std::fprintf(f, "ROW %d lo=%.17g up=%.17g\n", i, work.rmin[i],
                         work.rmax[i]);
        for (int j = 0; j < work.n; ++j)
            std::fprintf(f, "COL %d lo=%.17g up=%.17g c=%.17g\n", j,
                         work.cl[j], work.cu[j], work.c[j]);
        for (int j = 0; j < work.n; ++j)
            for (int p = work.cp[j]; p < work.cp[j + 1]; ++p)
                std::fprintf(f, "E %d %d %.17g\n", work.ci[p], j,
                             work.acx[p]);
        fclose(f);
    }


    auto snap_and_resolve = [&]() {
        for (long j = 0; j < E.total; ++j) {
            if (st[j] == BASIC) continue;
            if (st[j] == AT_LOWER && std::isfinite(E.lo[j]))
                z[j] = E.lo[j];
            else if (st[j] == AT_UPPER && std::isfinite(E.up[j]))
                z[j] = E.up[j];
        }
        build_and_factor();
        refresh_values();
    };

    int rc1 = run_simplex(true);
    if (rc1 != 0) {
        sol.iterations = iter;
        sol.solve_time_ms = elapsed() * 1000.0;
        return sol;
    }
    snap_and_resolve();
    {
        // Phase-1 optimum with positive artificials = infeasibility
        // certificate (elastic phase 1 cannot drive row violations to
        // zero). Without this check the art drive-out below silently
        // pins infeasible arts to zero and the final original-model
        // guard reports a bogus Error.
        double asum = 0.0;
        for (int i = 0; i < E.m; ++i)
            if (E.kind[basis[i]] == KIND_ART)
                asum += std::fabs(zb[i]);
        if (asum > 1e-6 * (1.0 + E.cmax_x)) {
            sol.status = Status::Infeasible;
            sol.message = "infeasible: phase-1 optimum retains " +
                          std::to_string(asum) + " artificial magnitude";
            sol.iterations = iter;
            sol.solve_time_ms = elapsed() * 1000.0;
            return sol;
        }
    }
    {
        double bviol = 0.0;
        for (int i = 0; i < E.m; ++i)
            bviol = std::max(bviol,
                             std::max(E.lo[basis[i]] - zb[i],
                                      zb[i] - E.up[basis[i]]));
        if (bviol > 1e-6 * (1.0 + E.cmax_x)) {
            sol.status = Status::Error;
            sol.message = "phase 1 ended with basic bound violation " +
                          std::to_string(bviol);
            sol.iterations = iter;
            sol.solve_time_ms = elapsed() * 1000.0;
            return sol;
        }
    }

    if (opt.verbosity > 0) {
        int nb = 0;
        for (int i = 0; i < E.m; ++i)
            if (E.kind[basis[i]] == KIND_ART) ++nb;
        std::fprintf(stderr, "[simplex] phase1 done iter=%ld arts_basic=%d\n",
                     iter, nb);
        for (long j = 0; j < E.total; ++j) {
            if (st[j] == BASIC) continue;
            bool bad = (st[j] == AT_LOWER &&
                        std::fabs(z[j] - E.lo[j]) > 1e-6 &&
                        std::isfinite(E.lo[j])) ||
                       (st[j] == AT_UPPER &&
                        std::fabs(z[j] - E.up[j]) > 1e-6 &&
                        std::isfinite(E.up[j]));
            if (bad)
                std::fprintf(stderr,
                             "  [pre-driveout DESYNC] var=%ld kind=%d "
                             "st=%u z=%.6g lo=%.6g up=%.6g\n",
                             j, (int)E.kind[j], (unsigned)st[j], z[j],
                             E.lo[j], E.up[j]);
        }
    }

    for (int i = 0; i < E.m; ++i) {
        if (E.kind[basis[i]] != KIND_ART) continue;
        vector<double> unit(E.m, 0.0);
        unit[i] = 1.0;
        vector<double> erow;
        lu.solve_transpose(unit, erow);
        int pick = -1;
        double best_abs = 1e-9;
        for (long j = 0; j < E.total; ++j) {
            if (st[j] == BASIC || E.kind[j] == KIND_ART) continue;
            double aij = 0.0;
            switch (E.kind[j]) {
                case KIND_X:
                    for (int p = E.cp[j]; p < E.cp[j + 1]; ++p)
                        aij += erow[E.ci[p]] * E.cv[p];
                    break;
                case KIND_SLACK:
                    aij = -erow[j - E.nx];
                    break;
                default:
                    break;
            }
            if (std::fabs(aij) > best_abs) {
                best_abs = std::fabs(aij);
                pick = (int)j;
            }
        }
        if (pick >= 0) {
            int old_art = basis[i];
            basis[i] = pick;
            st[pick] = BASIC;
            z[old_art] = 0.0;
            st[old_art] = AT_LOWER;
            build_and_factor();
            refresh_values();
        } else {
            E.lo[basis[i]] = 0.0;
            E.up[basis[i]] = 0.0;
            z[basis[i]] = std::min(std::max(zb[i], E.lo[basis[i]]),
                                   E.up[basis[i]]);
        }
    }

    for (long j = (long)work.n + work.m; j < E.total; ++j) {
        E.lo[j] = 0.0;
        E.up[j] = 0.0;
        E.cost[j] = 0.0;
        if (st[j] != BASIC) z[j] = 0.0;
    }
    for (int j = 0; j < work.n; ++j) E.cost[j] = E.xs_cost[j];

    snap_and_resolve();
    bool skip_phase2 = getenv("IGAOS_NO_PHASE2") != nullptr;
    if (skip_phase2) {
        build_and_factor();
        refresh_values();
    } else {
        int rc2 = run_simplex(false);
        if (rc2 != 0) {
            sol.iterations = iter;
            sol.solve_time_ms = elapsed() * 1000.0;
            return sol;
        }
        if (perturbed) {
            // cleanup pass: restore true costs and re-optimize from the
            // perturbed optimum (bounded damage — the perturbed point is
            // near-optimal, so few pivots)
            for (int j = 0; j < E.nx; ++j) E.cost[j] = E.xs_cost[j];
            perturbed = false;
            int rc3 = run_simplex(false);
            if (rc3 != 0) {
                sol.iterations = iter;
                sol.solve_time_ms = elapsed() * 1000.0;
                return sol;
            }
        }
    }

    snap_and_resolve();
    sol.x.assign(work.n, 0.0);
    for (int i = 0; i < E.m; ++i)
        if (E.kind[basis[i]] == KIND_X)
            sol.x[basis[i]] = zb[i] * E.sc.col[basis[i]];
    for (long j = 0; j < (long)work.n; ++j)
        if (st[j] != BASIC) sol.x[j] = z[j] * E.sc.col[j];
    if (getenv("IGAOS_DUMP")) {
        FILE* f = fopen("/tmp/opencode/simplex_final.txt", "w");
        std::fprintf(f, "SCALED SOLUTION (zb):\n");
        for (int i = 0; i < E.m; ++i)
            std::fprintf(f, "BAS pos=%d var=%ld kind=%d val=%.17g "
                            "lo=%.17g up=%.17g\n",
                         i, basis[i], (int)E.kind[basis[i]], zb[i],
                         E.lo[basis[i]], E.up[basis[i]]);
        for (long j = 0; j < E.total; ++j)
            if (st[j] != BASIC)
                std::fprintf(f, "NB var=%ld kind=%d st=%d val=%.17g "
                                "lo=%.17g up=%.17g\n",
                             j, (int)E.kind[j], (int)st[j], z[j], E.lo[j],
                             E.up[j]);
        fclose(f);
    }
    vector<double> cb(E.m, 0.0);
    for (int i = 0; i < E.m; ++i) cb[i] = E.cost[basis[i]];
    vector<double> ysc;
    lu.solve_transpose(cb, ysc);
    sol.y.resize(work.m);
    sol.row_activity.assign(work.m, 0.0);
    double obj = 0.0;
    for (int j = 0; j < work.n; ++j) obj += work.c[j] * sol.x[j];
    obj += work.obj_const;
    for (int i = 0; i < work.m; ++i) {
        sol.y[i] = E.sc.row[i] * ysc[i];
        for (int p = work.ap[i]; p < work.ap[i + 1]; ++p)
            sol.row_activity[i] += work.ax[p] * sol.x[work.ai[p]];
    }
    double internal_resid = 0.0;
    {
        vector<double> act(E.m, 0.0), col;
        vector<int> bpos(E.total, -1);
        for (int i = 0; i < E.m; ++i)
            if (basis[i] >= 0 && basis[i] < (int)E.total)
                bpos[basis[i]] = i;
        for (long j = 0; j < E.total; ++j) {
            double vj = (st[j] == BASIC) ? zb[bpos[j]] : z[j];
            if (vj == 0.0) continue;
            E.col_dense((int)j, col);
            for (int r = 0; r < E.m; ++r) act[r] += col[r] * vj;
        }
        for (int r = 0; r < E.m; ++r)
            internal_resid = std::max(internal_resid,
                                      std::fabs(act[r]));
    }
    if (opt.verbosity > 0 || internal_resid > 1e-6)
        std::fprintf(stderr, "[simplex] internal |Mz|inf=%.6g\n",
                     internal_resid);
    double orig_viol = 0.0;
    for (int j = 0; j < model.n; ++j)
        orig_viol = std::max(orig_viol,
                             std::max(model.cl[j] - sol.x[j],
                                      sol.x[j] - model.cu[j]));
    for (int i = 0; i < model.m; ++i)
        orig_viol = std::max(orig_viol,
                             std::max(model.rmin[i] - sol.row_activity[i],
                                      sol.row_activity[i] - model.rmax[i]));
    if (opt.verbosity > 0)
        std::fprintf(stderr, "[simplex] final orig_viol=%.6g\n", orig_viol);
    if (orig_viol > 1e-6 * (1.0 + std::fabs(obj))) {
        if (opt.verbosity > 0) {
            int shown = 0;
            for (int i = 0; i < model.m && shown < 5; ++i)
                if (model.rmin[i] - sol.row_activity[i] > 1e-6 ||
                    sol.row_activity[i] - model.rmax[i] > 1e-6) {
                    std::fprintf(stderr,
                                 "[simplex] VIOL row %d act=%.4f "
                                 "[%.4f,%.4f]\n",
                                 i, sol.row_activity[i], model.rmin[i],
                                 model.rmax[i]);
                    ++shown;
                }
        }
        sol.status = Status::Error;
        sol.message = "final solution violates original model by " +
                      std::to_string(orig_viol);
    }
    sol.objective = obj;
    sol.pinf = 0.0;
    sol.dinf = 0.0;
    sol.rel_gap = 0.0;
    if (orig_viol <= 1e-6 * (1.0 + std::fabs(obj))) {
        sol.status = Status::Optimal;
    } else {
        sol.status = Status::Error;
        sol.message = "final solution violates original model by " +
                      std::to_string(orig_viol);
    }
    sol.message = perturbed ? "optimal (cost perturbation applied)"
                            : "optimal";
    sol.iterations = iter;
    sol.solve_time_ms = elapsed() * 1000.0;
    if (warm_out != nullptr && sol.status == Status::Optimal) {
        // snapshot: structural+slack basis slots (arts recorded as -1 —
        // restore treats them as empty slots) and nonbasic states
        warm_out->basis.assign(work.m, -1);
        for (int i = 0; i < work.m; ++i)
            if (basis[i] < work.n + work.m) warm_out->basis[i] = basis[i];
        warm_out->nb_state.assign(st.begin(),
                                  st.begin() + (work.n + work.m));
    }
    return sol;
}

}  // namespace igaos::simplex
