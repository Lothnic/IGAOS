#include "simplex.hpp"

#include "dense_lu.hpp"
#include "sparse_lu.hpp"
#include "presolve.hpp"
#include "scaling.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <limits>
#include <map>
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

    // append basis column bj to a CSC assembly (ri/rv) — never dense
    void basis_col_sparse(int bj, vector<int>& ri, vector<double>& rv) const {
        switch (kind[bj]) {
            case KIND_X:
                for (int p = cp[bj]; p < cp[bj + 1]; ++p) {
                    if (cv[p] == 0.0) continue;
                    ri.push_back(ci[p]);
                    rv.push_back(cv[p]);
                }
                break;
            case KIND_SLACK:
                ri.push_back(bj - nx);
                rv.push_back(-1.0);
                break;
            case KIND_ART:
                ri.push_back(art_row[bj - nx - m]);
                rv.push_back(art_sigma[bj - nx - m]);
                break;
        }
    }
};

}  // namespace

    // factor the basis in CSC form — sparse path, never materializes the
    // dense m x m matrix (the O(m^2) memory wall on hanoi5-class instances)
    bool factor_basis_csc(const Engine& E, const vector<int>& basis,
                          SparseLU& lu) {
        vector<int> cp(E.m + 1, 0), ri;
        vector<double> rv;
        for (int c = 0; c < E.m; ++c) {
            if (basis[c] < 0) {
                ri.push_back(c);  // empty slot: unit column e_c
                rv.push_back(1.0);
            } else {
                E.basis_col_sparse(basis[c], ri, rv);
            }
            cp[c + 1] = (int)ri.size();
        }
        lu.factor_csc(E.m, cp, ri, rv);
        return lu.ok;
    }

static Solution solve_impl(const io::Model& model, const Options& opt,
                           const WarmStart* warm, WarmStart* warm_out,
                           std::vector<CutRow>* cuts_out,
                           const std::vector<unsigned char>* true_integ) {
    Solution sol;
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             t0)
            .count();
    };

    io::Model work = model;
    PresolveInfo pinfo;
    PresolveLog plog;
    std::vector<int> col_to_red, row_to_red;
    // Full A&A presolve only on the cold, matrix-stable path: warm starts
    // and Gomory cut generation index the original matrix, and the
    // postsolve wiring below belongs to the plain-solve tail.
    bool full_pre = opt.presolve && warm == nullptr && warm_out == nullptr &&
                    cuts_out == nullptr && !model.has_quadratic();
    if (full_pre) {
        io::Model red;
        PresolveStats pst;
        int prc = run_presolve(model, red, plog, col_to_red, row_to_red, pst);
        if (prc == 1) {
            sol.status = Status::Infeasible;
            sol.message = "infeasible: detected in presolve";
            sol.solve_time_ms = elapsed() * 1000.0;
            return sol;
        }
        if (prc == 2) {
            sol.status = Status::Unbounded;
            sol.message = "unbounded: detected in presolve";
            sol.solve_time_ms = elapsed() * 1000.0;
            return sol;
        }
        if (red.n == 0 && red.m > 0) {
            full_pre = false;  // unexpected shape; play safe
        } else {
            work = std::move(red);
            if (opt.verbosity > 0 || std::getenv("IGAOS_PRESOLVE_REPORT"))
                std::fprintf(stderr,
                             "[simplex] presolve: %dx%d/%d -> %dx%d/%d "
                             "(rows: %d red %d forced %d dup; cols: %d fix "
                             "%d empty %d subst %d dup; %d tightens)\n",
                             pst.m0, pst.n0, pst.nnz0, pst.m1, pst.n1,
                             pst.nnz1, pst.rows_redundant, pst.rows_forcing,
                             pst.rows_dup, pst.cols_fixed, pst.cols_empty,
                             pst.cols_subst, pst.cols_dup, pst.bound_tightens);
        }
    }
    if (!full_pre) {
        int nch = strengthen_bounds(work, pinfo, 3);
        if (opt.verbosity > 0)
            std::fprintf(stderr,
                         "[simplex] presolve: %d bound tightenings\n", nch);
    }

    if (full_pre && work.m == 0) {
        // Presolve removed every row: reduced problem is a box LP (or
        // empty). Solve directly, then postsolve and validate against the
        // original model — same acceptance rule as the main path.
        bool unb = false;
        double obj = work.obj_const;
        std::vector<double> xr(work.n, 0.0);
        for (int j = 0; j < work.n && !unb; ++j) {
            double v = 0.0;
            if (work.c[j] > 1e-12) {
                if (work.cl[j] == -INF) unb = true;
                else v = work.cl[j];
            } else if (work.c[j] < -1e-12) {
                if (work.cu[j] == INF) unb = true;
                else v = work.cu[j];
            } else {
                v = std::isfinite(work.cl[j])
                        ? work.cl[j]
                        : (std::isfinite(work.cu[j]) ? work.cu[j] : 0.0);
            }
            if (!unb) { xr[j] = v; obj += work.c[j] * v; }
        }
        if (unb) {
            sol.status = Status::Unbounded;
            sol.message = "unbounded: detected in presolve";
        } else {
            sol.x = postsolve_x(plog, xr, col_to_red, model.n);
            sol.y.assign(model.m, 0.0);
            sol.row_activity.assign(model.m, 0.0);
            for (int i = 0; i < model.m; ++i)
                for (int p = model.ap[i]; p < model.ap[i + 1]; ++p)
                    sol.row_activity[i] +=
                        model.ax[p] * sol.x[model.ai[p]];
            double viol = 0.0;
            for (int j = 0; j < model.n; ++j)
                viol = std::max(viol, std::max(model.cl[j] - sol.x[j],
                                               sol.x[j] - model.cu[j]));
            for (int i = 0; i < model.m; ++i)
                viol = std::max(viol,
                                std::max(model.rmin[i] - sol.row_activity[i],
                                         sol.row_activity[i] - model.rmax[i]));
            sol.objective = obj;
            sol.pinf = sol.dinf = sol.rel_gap = 0.0;
            if (viol <= 1e-6 * (1.0 + std::fabs(obj))) {
                sol.status = Status::Optimal;
                sol.message = "optimal (solved in presolve)";
            } else {
                sol.status = Status::Error;
                sol.message =
                    "presolve solution violates original model by " +
                    std::to_string(viol);
            }
        }
        sol.iterations = 0;
        sol.solve_time_ms = elapsed() * 1000.0;
        return sol;
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
    // with the same matrix (B&B child, bounds-only difference). Two paths:
    //  - dual_warm: snapshot basis is complete (no arts) — the parent
    //    optimum is dual-feasible, so the DUAL simplex repairs the child's
    //    primal bound violations directly (no phase 1, no art machinery).
    //  - otherwise: snap out-of-bound basics to bounds, empty slots become
    //    artificials, elastic phase 1 restores feasibility.
    bool warm_ok = warm != nullptr &&
                   (int)warm->basis.size() == work.m &&
                   (int)warm->nb_state.size() == work.n + work.m;
    bool dual_warm = false;
    if (warm_ok) {
        for (int i = 0; i < work.m; ++i) {
            long v = warm->basis[i];
            basis[i] = (v >= 0 && v < (long)(work.n + work.m)) ? (int)v : -1;
        }
        bool complete = true;
        for (int i = 0; i < work.m; ++i)
            if (basis[i] < 0) complete = false;
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
        SparseLU wlu;
        factor_basis_csc(E, basis, wlu);
        if (!wlu.ok) {
            warm_ok = false;  // singular restored basis — fall back cold
        } else if (complete) {
            dual_warm = true;  // dual simplex takes it from here
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
                SparseLU chk;
                if (!factor_basis_csc(E, basis, chk)) warm_ok = false;
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

    // Sparse base + eta file (dense fallback inside SparseLU, plus the
    // IGAOS_DENSE_LU=1 escape hatch that forces the old dense factor)
    const bool force_dense_lu = getenv("IGAOS_DENSE_LU") != nullptr;
    EtaFileT<SparseLU> lu;
    long iter = 0;  // declared before the lambdas that trace/capture it
    bool basis_dirty = true;  // eta file must be rebuilt from scratch
    vector<double> zb(E.m, 0.0);
    auto build_and_factor = [&]() {
        lu.reset();
        lu.m = E.m;
        if (force_dense_lu) {
            vector<double> Bm((size_t)E.m * E.m, 0.0), col;
            for (int c = 0; c < E.m; ++c) {
                E.basis_col(basis[c], col);
                for (int r = 0; r < E.m; ++r)
                    if (col[r] != 0.0) Bm[(size_t)r * E.m + c] = col[r];
            }
            lu.base.factor_dense(std::move(Bm));
        } else {
            factor_basis_csc(E, basis, lu.base);
        }
        lu.ok = lu.base.ok;
        if (getenv("IGAOS_TRACE_PIVOTS")) {
            double dmin = INF, dmax = 0.0;
            for (int i = 0; i < E.m; ++i) {
                double d = lu.base.pivot_mag(i);
                dmin = std::min(dmin, d);
                dmax = std::max(dmax, d);
            }
            std::fprintf(stderr, "FACT it=%ld ok=%d dmin=%.3e dmax=%.3e "
                                 "ratio=%.3e\n",
                         iter, (int)lu.ok, dmin, dmax,
                         dmax > 0 ? dmin / dmax : 0.0);
        }
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
    int degen_streak = 0;
    int crawl_steps = 0;  // cumulative eps-steps since last perturbation
    double expand_tol = 0.0;  // EXPAND ratio-test tolerance (degenerate runs)
    bool bland = false;
    bool perturbed = false;
    // bound perturbation state (deep-stall anti-degeneracy, see below)
    vector<double> lo0, up0;
    bool bounds_perturbed = false;
    int bound_perturb_count = 0;
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
        expand_tol = 0.0;
        basis_dirty = true;  // force a factorization at phase entry
        // Devex pricing (Forrest–Goldfarb 1992): reference weights w_j >= 1,
        // entering column maximizes dj^2/w_j. Fresh framework (all 1) at each
        // phase entry; reset when weights blow up or on deep degeneracy.
        // IGAOS_DANZIG=1 in the environment falls back to plain Dantzig.
        const bool use_devex = (getenv("IGAOS_DANZIG") == nullptr);
        bool devex_on = use_devex;
        vector<double> dw((size_t)E.total, 1.0);
        const double DEVEX_RESET = getenv("IGAOS_DEVEX_RESET")
                                       ? atof(getenv("IGAOS_DEVEX_RESET"))
                                       : 1e6;
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
        bool rolled_back = false;
        // checkpoint: the basis/st at the last SUCCESSFUL factorization.
        // When a refactorization fails and single-pivot rollback doesn't
        // help (the corruption entered >1 pivot ago), restore this — it is
        // guaranteed to refactor (it did) and z/zb are recomputed from
        // scratch, clearing any accumulated eta-file desync.
        vector<int> chk_basis;
        vector<unsigned char> chk_st;
        bool chk_valid = false;
        // entering columns committed since the last verified factorization —
        // when a checkpoint restore fires, these are blocked: one of them
        // poisoned the basis, and without blocking them pricing re-attempts
        // the identical pivot sequence forever (observed on wood1p)
        vector<int> entered_since_chk;
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
                            // roll it back and reject that entering column.
                            // Block level 2 (persistent): clearing it on the
                            // next pivot attempt would re-admit the column
                            // from the SAME restored basis and loop forever
                            // (observed: alternating singular pivots on
                            // scsd1 cols 66/88). Cleared only after a
                            // refactorization-valid pivot changes the basis.
                            basis = rec_basis;
                            st = rec_st;
                            skip_col[rec_enter] = 2;
                            rec_enter = -1;
                            rolled_back = true;
                            basis_dirty = true;  // etas are stale
                            continue;
                        }
                        if (chk_valid) {
                            // the bad pivot predates the last snapshot:
                            // fall back to the last verified basis and block
                            // every column committed since — one of them
                            // poisoned it
                            basis = chk_basis;
                            st = chk_st;
                            for (int j : entered_since_chk)
                                if (j >= 0 && j < (int)skip_col.size())
                                    skip_col[j] = 2;
                            entered_since_chk.clear();
                            basis_dirty = true;
                            continue;
                        }
                        sol.status = Status::Error;
                        sol.message = "basis factorization failed";
                        return 5;
                    }
                    basis_dirty = false;
                    chk_basis = basis;
                    chk_st = st;
                    chk_valid = true;
                    entered_since_chk.clear();
                    if (rolled_back) {
                        // this factorization validates the RESTORED basis —
                        // keep level-2 blocks; the next refactorization that
                        // follows committed pivots clears them
                        rolled_back = false;
                    } else {
                        for (auto& s : skip_col)
                            if (s == 2) s = 0;
                    }
                }
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
            double best_devex = 0.0;
            double dir0 = 1.0;
            for (long j = 0; j < E.total; ++j) {
                if (st[j] == BASIC) continue;
                if (skip_col[j]) continue;
                if (std::isfinite(E.lo[j]) && std::isfinite(E.up[j]) &&
                    E.up[j] - E.lo[j] <= 1e-12)
                    continue;
                double dj = E.reduced_cost_part((int)j, y);
                if (!devex_on) {
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
                    continue;
                }
                // Devex: same eligibility as Dantzig (|dj| above tolerance),
                // rank by dj^2/w_j (FREE_NB uses |dj|^2/w_j).
                double mag;
                if (st[j] == FREE_NB) {
                    mag = std::fabs(dj);
                } else if (st[j] == AT_LOWER) {
                    mag = -dj;
                } else {
                    mag = dj;
                }
                if (mag <= best_score) continue;
                double score = (dj * dj) / dw[j];
                if (score > best_devex) {
                    best_devex = score;
                    enter0 = (int)j;
                    dir0 = (st[j] == FREE_NB) ? (dj > 0 ? -1.0 : 1.0)
                              : (st[j] == AT_LOWER ? 1.0 : -1.0);
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
            if (enter0 < 0) {
                // Terminal honesty guard: reduced costs from a long eta
                // file carry ~1e-6 noise (observed: an exactly-zero dj on a
                // paired column reading as -1e-6 and manufacturing a false
                // unbounded ray). Never declare optimality off stale etas —
                // refactor once and reprice. Terminates: the refactor
                // empties the eta file, so the second visit passes.
                if (!lu.etas.empty()) {
                    basis_dirty = true;
                    continue;
                }
                return 0;
            }

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
            // Absolute pivot floor: pivoting on |alpha| below this wrecks
            // the basis (observed: a 1.9e-8 pivot producing B^{-1} entries
            // of 1e8 and garbage duals on scsd1). The relative Harris
            // stability check (1e-7*amax) does not catch these when amax is
            // itself blown up or when the Bland/accept_unsafe fallbacks
            // bypass it — every pivot selection path enforces this floor.
            const double PIVOT_ABS = 1e-7;
            const double DELTA_BASE = 1e-8;
            // EXPAND anti-degeneracy: expand_tol (updated post-pivot) grows
            // each degenerate iteration so zero-distance blockers become
            // eligible pivots one by one — the basis sequence cannot repeat
            const double DELTA = DELTA_BASE + expand_tol;
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
                // Harris pass 2: largest |alpha| among candidates within the
                // relaxed minimum ratio (best pivot).
                double best_abs = STABILITY;
                for (int i = 0; i < E.m; ++i) {
                    if (std::fabs(alpha[i]) < best_abs ||
                        std::fabs(alpha[i]) <= PIVOT_ABS)
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
                        if (std::fabs(alpha[i]) <= PIVOT_ABS) continue;
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
                            if (std::fabs(alpha[i]) <= PIVOT_ABS)
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
                    double fb_abs = PIVOT_ABS;
                    for (int i = 0; i < E.m; ++i) {
                        if (std::fabs(alpha[i]) <= PIVOT_ABS) continue;
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
                if (phase1) {
                    // phase-1 unboundedness is impossible (art objective
                    // bounded below by 0) — always a numerical artifact;
                    // skip the column and reprice
                    skip_col[enter0] = 1;
                    continue;
                }
                if (std::isfinite(theta_rel) && !accept_unsafe) {
                    // blocked (finite ratio exists) but no stable pivot —
                    // reject this column, reprice
                    skip_col[enter0] = 1;
                    continue;
                }
                if (!retry_done) {
                    retry_done = true;
                    // re-examine the column from a FRESH factorization —
                    // unboundedness must never be declared off stale etas
                    // (see the terminal honesty guard above)
                    basis_dirty = true;
                    continue;
                }
                if (opt.verbosity > 0) {
                    double djchk = E.reduced_cost_part(enter0, y);
                    double ymax = 0.0;
                    for (double v : y) ymax = std::max(ymax, std::fabs(v));
                    std::fprintf(stderr,
                                 "[RAY] enter=%d kind=%d lo=%.4g up=%.4g "
                                 "amax=%.3g dj=%.6g ymax=%.3g etas=%zu "
                                 "bland=%d\n",
                                 enter0, (int)E.kind[enter0], E.lo[enter0],
                                 E.up[enter0], amax, djchk, ymax,
                                 lu.etas.size(), (int)bland);
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

            // Epsilon-crawl counts as degenerate: an eps-scale step
            // (theta ~ 1e-7) is zero progress for practical purposes —
            // d6cube phase 2 takes 100k+ such steps, each resetting the
            // streak, so the anti-degeneracy machinery never fired. Any
            // step below 1e-6 of the entering variable's own scale feeds
            // the streak; only real progress resets it.
            double eps_step = 1e-6 * (1.0 + std::fabs(z[enter0]));
            if (theta > eps_step) {
                degen_streak = 0;
                bland = false;
                expand_tol = 0.0;
            } else {
                ++degen_streak;
                ++crawl_steps;
                expand_tol = std::min(expand_tol + DELTA_BASE,
                                      1e-4 * (1.0 + E.cmax_x));
                if (degen_streak >= 20 && !bland) {
                    bland = true;
                    degen_streak = 0;
                    // degenerate stall: the reference framework has drifted —
                    // restart Devex weights from a fresh one
                    if (devex_on) std::fill(dw.begin(), dw.end(), 1.0);
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
                // deeper stall: primal degeneracy crawl (eps steps on >70%
                // of pivots, e.g. forplan/d6cube). Expand every finite
                // bound of structural+slack variables by a random epsilon
                // ABOVE the crawl scale (crawl steps run ~1e-7; a 1e-7
                // perturbation dissolves nothing): blocking basics sitting
                // on a bound gain a nonzero distance, so the ratio-test
                // tie structure that produces zero-progress pivots
                // dissolves. True bounds are restored after each phase and
                // a cleanup pass re-verifies optimality (bounds do not
                // affect reduced costs, so the perturbed optimum is
                // eps-optimal for the true problem).
                // CUMULATIVE crawl counter (not consecutive streak): real
                // steps interleave the crawl and reset a streak long
                // before 1000, so the consecutive trigger fired once per
                // 80k iterations on d6cube. Re-perturbed (fresh randomness
                // off the TRUE bounds) when the crawl resumes — a new draw
                // escapes a different tie structure.
                if (crawl_steps >= 2000 && bound_perturb_count < 16) {
                    // deep degeneracy crawl that Devex weights are not
                    // helping (observed: modszk1 4k -> 31k iterations) —
                    // drop the weight heuristic for the rest of the phase
                    devex_on = false;
                    if (!bounds_perturbed) {
                        bounds_perturbed = true;
                        lo0 = E.lo;
                        up0 = E.up;
                    } else {
                        E.lo = lo0;
                        E.up = up0;
                    }
                    ++bound_perturb_count;
                    if (opt.verbosity > 0)
                        std::fprintf(stderr,
                                     "[simplex] bound perturbation #%d "
                                     "applied it=%ld ph1=%d\n",
                                     bound_perturb_count, iter, (int)phase1);
                    for (long j = 0; j < E.total; ++j) {
                        if (E.kind[j] == KIND_ART) continue;
                        bool fin_lo = std::isfinite(E.lo[j]);
                        bool fin_up = std::isfinite(E.up[j]);
                        if (fin_lo && fin_up &&
                            E.up[j] - E.lo[j] <= 1e-12)
                            continue;  // fixed vars stay fixed
                        if (fin_lo) {
                            double r = (double)rng() / (double)rng.max();
                            E.lo[j] -= (1e-7 + 9e-7 * r) *
                                       (1.0 + std::fabs(E.lo[j]));
                        }
                        if (fin_up) {
                            double r = (double)rng() / (double)rng.max();
                            E.up[j] += (1e-7 + 9e-7 * r) *
                                       (1.0 + std::fabs(E.up[j]));
                        }
                    }
                    degen_streak = 0;
                    crawl_steps = 0;
                }
            }
            // stability-skips (level 1) are retried on every new pivot;
            // rollback-blocks (level 2) persist until a verified basis
            // change (see the factorization block above)
            for (auto& s : skip_col)
                if (s == 1) s = 0;
            accept_unsafe = false;

            if (flip) {
                st[enter0] = (dir0 > 0) ? AT_UPPER : AT_LOWER;
                z[enter0] = (dir0 > 0) ? E.up[enter0] : E.lo[enter0];
                verify_pivot_invariant("flip", enter0, -1, -1, theta, dir0);
                continue;
            }

            // Devex weight update (Forrest–Goldfarb 1992): with pivot row r,
            // entering q, pivot alpha_rq:
            //   w_j = max(w_j, (alpha_rj/alpha_rq)^2 * w_q)  (nonbasic j)
            //   w_leaving = max(w_q/alpha_rq^2, 1)
            // alpha_rj = (e_r^T B^{-1}) A_j; rho = B^{-T} e_r gives the row
            // via the row-major matrix copy (slack/art columns handled
            // directly). Framework resets when any weight exceeds 1e6.
            if (devex_on && leave_pos >= 0) {
                vector<double> unit(E.m, 0.0);
                unit[leave_pos] = 1.0;
                vector<double> rho;
                lu.solve_transpose(unit, rho);
                const double arq = alpha[leave_pos];
                const double wq = dw[enter0];
                const double wq_arq2 = wq / (arq * arq);
                vector<double> prow((size_t)E.total, 0.0);
                for (int i = 0; i < E.m; ++i) {
                    if (rho[i] == 0.0) continue;
                    const double ri = rho[i];
                    for (int p = E.ap[i]; p < E.ap[i + 1]; ++p)
                        prow[E.ai[p]] += ri * E.av[p];
                }
                double wmax = 1.0;
                for (long j = 0; j < E.total; ++j) {
                    if (st[j] == BASIC) continue;
                    double arj = 0.0;
                    if (j < E.nx) {
                        arj = prow[j];
                    } else if (E.kind[j] == KIND_SLACK) {
                        arj = -rho[j - E.nx];
                    } else {  // KIND_ART
                        long k = j - E.nx - E.m;
                        arj = rho[E.art_row[k]] * E.art_sigma[k];
                    }
                    if (arj != 0.0)
                        dw[j] = std::max(dw[j], arj * arj * wq_arq2);
                    wmax = std::max(wmax, dw[j]);
                }
                int leaving_var = basis[leave_pos];
                dw[leaving_var] = std::max(wq_arq2, 1.0);
                wmax = std::max(wmax, dw[leaving_var]);
                if (wmax > DEVEX_RESET)
                    std::fill(dw.begin(), dw.end(), 1.0);
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
            entered_since_chk.push_back(enter0);
            zb[leave_pos] = z[enter0];
            // product-form update: B_new = B_old · E(leave_pos, alpha)
            lu.update(leave_pos, alpha);
            if (getenv("IGAOS_TRACE_PIVOTS"))
                std::fprintf(stderr,
                             "PIV it=%ld ph1=%d enter=%d leave=%ld@%d "
                             "theta=%.3e ap=%.3e amax=%.3e\n",
                             iter, (int)phase1, enter0, leaving, leave_pos,
                             theta, std::fabs(alpha[leave_pos]),
                             [&] {
                                 double a = 0.0;
                                 for (double v : alpha)
                                     a = std::max(a, std::fabs(v));
                                 return a;
                             }());
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

    // Bounded-variable dual simplex (from mathematical foundations):
    // start from a dual-feasible basis (parent optimum), repair primal
    // bound violations. Leaving row = largest scaled violation; entering
    // column = dual ratio-test winner among eligible nonbasics. Dual
    // unboundedness (no eligible column) certifies primal infeasibility.
    auto run_dual_simplex = [&]() -> int {
        std::fill(skip_col.begin(), skip_col.end(), 0);
        accept_unsafe = false;
        expand_tol = 0.0;
        basis_dirty = true;
        vector<double> rho(E.m), cb(E.m), y, unit(E.m, 0.0),
            ae(E.m), acol(E.m);
        for (; iter <= opt.max_iterations; ++iter) {
            if (elapsed() > opt.time_limit_s) {
                sol.status = Status::TimeLimit;
                sol.message = "dual simplex: wall-clock budget exhausted";
                return 3;
            }
            for (long j = 0; j < E.total; ++j) {
                if (st[j] == BASIC) continue;
                if (st[j] == AT_LOWER && std::isfinite(E.lo[j]))
                    z[j] = E.lo[j];
                else if (st[j] == AT_UPPER && std::isfinite(E.up[j]))
                    z[j] = E.up[j];
            }
            if (basis_dirty || lu.needs_refactor()) {
                if (!build_and_factor()) {
                    sol.status = Status::Error;
                    sol.message = "dual simplex: basis factorization failed";
                    return 2;
                }
            }
            refresh_values();

            // leaving row: largest scaled bound violation
            int r = -1;
            double best_v = 1e-9;  // violation ladder floor
            for (int i = 0; i < E.m; ++i) {
                long v = basis[i];
                double viol = std::max(E.lo[v] - zb[i], zb[i] - E.up[v]);
                if (viol <= 0.0) continue;
                double bmag = 1.0 + std::max(
                    std::isfinite(E.lo[v]) ? std::fabs(E.lo[v]) : 0.0,
                    std::isfinite(E.up[v]) ? std::fabs(E.up[v]) : 0.0);
                double sv = viol / bmag;
                if (sv > best_v) { best_v = sv; r = i; }
            }
            if (r < 0) return 0;  // primal feasible + dual feasible

            // dual ratio test on row r
            std::fill(unit.begin(), unit.end(), 0.0);
            unit[r] = 1.0;
            lu.solve_transpose(unit, rho);
            for (int i = 0; i < E.m; ++i) cb[i] = E.cost[basis[i]];
            lu.solve_transpose(cb, y);

            bool below = zb[r] < E.lo[basis[r]];
            int enter = -1;
            double best_ratio = INF, best_arj = 0.0;
            for (long j = 0; j < E.total; ++j) {
                if (st[j] == BASIC) continue;
                if (skip_col[j]) continue;
                if (std::isfinite(E.lo[j]) && std::isfinite(E.up[j]) &&
                    E.up[j] - E.lo[j] <= 1e-12)
                    continue;
                // alpha_rj = rho . A_j
                double arj = 0.0;
                switch (E.kind[j]) {
                    case KIND_X:
                        for (int p = E.cp[j]; p < E.cp[j + 1]; ++p)
                            arj += rho[E.ci[p]] * E.cv[p];
                        break;
                    case KIND_SLACK:
                        arj = -rho[j - E.nx];
                        break;
                    case KIND_ART:
                        arj = rho[E.art_row[j - E.nx - E.m]] *
                              E.art_sigma[j - E.nx - E.m];
                        break;
                }
                if (std::fabs(arj) <= 1e-9) continue;
                // eligibility: d zb_r = -arj * t_j must be positive (below)
                // resp. negative (above); t_j > 0 from AT_LOWER, < 0 from
                // AT_UPPER, either for FREE_NB
                bool elig;
                if (st[j] == FREE_NB) {
                    elig = true;
                } else if (below) {
                    elig = (st[j] == AT_LOWER) ? (arj < 0) : (arj > 0);
                } else {
                    elig = (st[j] == AT_LOWER) ? (arj > 0) : (arj < 0);
                }
                if (!elig) continue;
                double dj = E.reduced_cost_part((int)j, y);
                double ratio = std::fabs(dj) / std::fabs(arj);
                // stability: among near-ties prefer the larger |arj|
                if (ratio < best_ratio * (1.0 - 1e-9) ||
                    (ratio < best_ratio * (1.0 + 1e-9) &&
                     std::fabs(arj) > std::fabs(best_arj))) {
                    best_ratio = ratio;
                    best_arj = arj;
                    enter = (int)j;
                }
            }
            if (enter < 0) {
                // dual unbounded: no column can repair row r
                sol.status = Status::Infeasible;
                sol.message =
                    "infeasible: dual simplex found no repair column for "
                    "violated row";
                return 1;
            }

            // pivot: entering takes slot r; leaving exits to its violated
            // bound
            int leaving = basis[r];
            z[leaving] = below ? E.lo[leaving] : E.up[leaving];
            st[leaving] = below ? AT_LOWER : AT_UPPER;
            basis[r] = enter;
            st[enter] = BASIC;
            // eta: full column B^{-1} a_enter
            E.col_dense(enter, ae);
            lu.solve(ae, acol);
            lu.update(r, acol);
        }
        sol.status = Status::IterationLimit;
        sol.message = "dual simplex: iteration cap reached";
        return 4;
    };

    bool dual_done = false;
    if (dual_warm) {
        for (int j = 0; j < E.nx; ++j) E.cost[j] = E.xs_cost[j];
        int rcd = run_dual_simplex();
        if (rcd == 0) {
            dual_done = true;
        } else {
            // dual failed (time/error/iteration): honest error — the B&B
            // caller cold-solves on Error
            sol.iterations = iter;
            sol.solve_time_ms = elapsed() * 1000.0;
            if (rcd != 1) return sol;  // infeasible already carries status
            return sol;
        }
    }

    int rc1 = -1;
    if (!dual_done) rc1 = run_simplex(true);
    if (!dual_done && rc1 != 0) {
        sol.iterations = iter;
        sol.solve_time_ms = elapsed() * 1000.0;
        return sol;
    }
    if (bounds_perturbed) {
        // Phase 1 ran on perturbed bounds: the phase-1 optimum is
        // eps-infeasible for the TRUE bounds (observed: 2e-6 basic
        // violation on d6cube failing the strict post-phase check) and
        // phase 2 would inherit eps-out-of-bounds basics — the ratio test
        // then sees negative distances clamped to zero, an instant
        // degenerate crawl. Restore true bounds and clean up with a short
        // phase-1 re-run from the near-feasible point.
        E.lo = lo0;
        E.up = up0;
        bounds_perturbed = false;
        bound_perturb_count = 0;
        basis_dirty = true;
        int rcb = run_simplex(true);
        if (rcb != 0) {
            sol.iterations = iter;
            sol.solve_time_ms = elapsed() * 1000.0;
            return sol;
        }
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
        if (bounds_perturbed) {
            // cleanup pass: restore true bounds. The perturbed optimum is
            // dual-feasible with true costs (reduced costs are
            // bounds-independent) but eps-primal-infeasible — exactly the
            // starting state the bounded-variable DUAL simplex repairs.
            E.lo = lo0;
            E.up = up0;
            bounds_perturbed = false;
            basis_dirty = true;
            int rc4 = run_dual_simplex();
            // rc4==1 ("dual unbounded") after a bound restore is NOT a
            // trustworthy infeasibility certificate (eps-scale violations,
            // numerics): fall through and let the final original-model
            // violation check decide honestly. Other failures (time/iter/
            // factorization) return as-is.
            if (rc4 != 0 && rc4 != 1) {
                sol.iterations = iter;
                sol.solve_time_ms = elapsed() * 1000.0;
                return sol;
            }
        }
    }

    snap_and_resolve();
    vector<double> xw(work.n, 0.0);  // reduced-space primal solution
    for (int i = 0; i < E.m; ++i)
        if (E.kind[basis[i]] == KIND_X)
            xw[basis[i]] = zb[i] * E.sc.col[basis[i]];
    for (long j = 0; j < (long)work.n; ++j)
        if (st[j] != BASIC) xw[j] = z[j] * E.sc.col[j];
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
    // Objective of the reduced solve: presolve folded the removed columns'
    // constant contributions into work.obj_const, so this is already the
    // original-space objective.
    double obj = work.obj_const;
    for (int j = 0; j < work.n; ++j) obj += work.c[j] * xw[j];
    if (full_pre) {
        sol.x = postsolve_x(plog, xw, col_to_red, model.n);
        sol.y.assign(model.m, 0.0);
        for (int i = 0; i < model.m; ++i)
            if (row_to_red[i] >= 0)
                sol.y[i] = E.sc.row[row_to_red[i]] * ysc[row_to_red[i]];
        // ponytail: duals of removed rows stay 0 — objective/feasibility
        // are exact, true original duals are not reconstructed (documented
        // limitation; benchmarks validate objective + feasibility only)
    } else {
        sol.x = xw;
        sol.y.resize(work.m);
        for (int i = 0; i < work.m; ++i)
            sol.y[i] = E.sc.row[i] * ysc[i];
    }
    // row activity always recomputed in original space
    sol.row_activity.assign(model.m, 0.0);
    for (int i = 0; i < model.m; ++i)
        for (int p = model.ap[i]; p < model.ap[i + 1]; ++p)
            sol.row_activity[i] += model.ax[p] * sol.x[model.ai[p]];
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
    const std::vector<unsigned char>& cut_integ =
        true_integ != nullptr ? *true_integ : model.integ;
    if (cuts_out != nullptr && sol.status == Status::Optimal) {
        // Gomory mixed-integer cuts from fractional basic integer
        // variables. Tableau rows fold nonbasic slacks (s_i = A_i x)
        // back into structural space; nonbasic artificials are fixed at
        // 0 and skipped. Scaling is power-of-two so unscaling is exact.
        const int MAX_CUTS = 20;
        // instrumentation: why rows are rejected (IGAOS_DEBUG_CUTS)
        long dbg_rows = 0, dbg_nonstruct = 0, dbg_nonint = 0,
             dbg_f0 = 0, dbg_art = 0, dbg_empty = 0, dbg_dyn = 0,
             dbg_dense = 0, dbg_emit = 0, dbg_artcols = 0;
        for (int i = 0; i < E.m; ++i)
            if (E.kind[basis[i]] == KIND_ART) ++dbg_artcols;
        vector<double> rho(E.m), unit(E.m);
        for (int r = 0; r < E.m && (int)cuts_out->size() < MAX_CUTS; ++r) {
            ++dbg_rows;
            long k = basis[r];
            if (k < 0 || k >= E.nx) { ++dbg_nonstruct; continue; }
            if (!cut_integ[k]) { ++dbg_nonint; continue; }
            double beta = zb[r] * E.sc.col[k];              // original value
            double f0 = beta - std::floor(beta);
            if (f0 < 0.01 || f0 > 0.99) { ++dbg_f0; continue; }
            std::fill(unit.begin(), unit.end(), 0.0);
            unit[r] = 1.0;
            lu.solve_transpose(unit, rho);
            // tableau row over nonbasic structurals and slacks (slack i is
            // the row activity s_i = A_i x, so its coefficient folds into
            // structural space). Nonbasic ARTIFICIALS are SKIPPED, not
            // disqualifying: after the drive-out every art is fixed at
            // [0,0] with value 0, so its tableau coefficient multiplies a
            // constant 0 — it contributes nothing to the cut and any
            // integer-feasible point of the original model extends to the
            // augmented system with arts at 0, so dropping the term is
            // exact. (This was the silent-generator bug: rows whose
            // B^-1 e_r touched an art row were all rejected.)
            bool usable = true;
            vector<std::pair<int, double>> row_terms;  // (j, a_j^orig)
            vector<std::pair<int, double>> slack_terms;  // (row i, a^orig)
            for (long j = 0; j < E.total && usable; ++j) {
                if (st[j] == BASIC) continue;
                if (E.kind[j] == KIND_ART) {
                    // sanity: art must be fixed at 0 to be droppable
                    if (!(E.lo[j] == 0.0 && E.up[j] == 0.0 &&
                          z[j] == 0.0))
                        usable = false;
                    continue;
                }
                double aj = 0.0;
                switch (E.kind[j]) {
                    case KIND_X:
                        for (int p = E.cp[j]; p < E.cp[j + 1]; ++p)
                            aj += rho[E.ci[p]] * E.cv[p];
                        break;
                    case KIND_SLACK:
                        aj = -rho[j - E.nx];
                        break;
                    case KIND_ART:
                        break;
                }
                if (std::fabs(aj) <= 1e-9) continue;
                if (st[j] == FREE_NB) {
                    // a free nonbasic is unconstrained in u-space — a
                    // nonzero coefficient would make the cut invalid
                    usable = false;
                    continue;
                }
                double sigma = (st[j] == AT_UPPER) ? -1.0 : 1.0;
                if (j < E.nx) {
                    // unscale into original space (powers of two: exact)
                    double a_orig = aj * sigma * E.sc.col[k] / E.sc.col[j];
                    row_terms.emplace_back((int)j, a_orig);
                } else {
                    int i = (int)(j - E.nx);
                    // scaled slack value = sc.row[i] * (A_i x), so the
                    // original-space tableau coefficient is
                    // alpha * sc.row[i] * sc.col[k] (was divided by
                    // sc.row[i] — silently invalid cuts on rows whose
                    // power-of-two row scale != 1)
                    double a_orig =
                        aj * sigma * E.sc.col[k] * E.sc.row[i];
                    slack_terms.emplace_back(i, a_orig);
                }
            }
            if (!usable || (row_terms.empty() && slack_terms.empty())) {
                if (!usable) ++dbg_art;  // unfixed art (should never fire)
                else ++dbg_empty;
                continue;
            }
            if (const char* dp = getenv("IGAOS_DUMP_CUTS")) {
                // raw tableau row for offline identity checks:
                // x_k(x) = beta - sum a_j sigma_j (x_j - x_j^cur)
                //          - sum a_si sigma_i (s_i(x) - s_i^cur)
                FILE* f = fopen(dp, "a");
                if (f) {
                    std::fprintf(f, "ROW k=%ld beta=%.17g f0=%.17g "
                                        "xk=%.17g\n",
                                 k, zb[r] * E.sc.col[k], f0, sol.x[k]);
                    for (auto& [j, a] : row_terms)
                        std::fprintf(f, "T %d %.17g %.17g %d\n", j, a,
                                     z[j] * E.sc.col[j], (int)st[j]);
                    for (auto& [i, a] : slack_terms) {
                        double act = 0.0;
                        for (int p = work.ap[i]; p < work.ap[i + 1]; ++p)
                            act += work.ax[p] * sol.x[work.ai[p]];
                        std::fprintf(f, "S %d %.17g %.17g %d\n", i, a, act,
                                     (int)st[E.nx + i]);
                    }
                    fclose(f);
                }
            }
            // GMI coefficients over structurals (integer formula) and
            // slacks (continuous formula). Cut: sum pi_j u_j >= f0 with
            // u_j = sigma_j (x_j - x_j^cur) >= 0; slacks substitute
            // s_i = A_i x into structural space. The stored tableau
            // coefficient a IS the u-space coefficient (alpha_j *
            // sigma_j, sigma folded at row-build time) — the formula is
            // applied to it directly and the x-space coefficient is
            // pi * sigma_j. The integer formula is only valid when u_j
            // is integer, i.e. the nonbasic integer sits at an INTEGRAL
            // value; at a fractional bound it must use the continuous
            // formula.
            //
            // Continuous coefficients (row convention x_k = beta -
            // sum a_j u_j): pi = a for a > 0, |a|*f0/(1-f0) for a < 0.
            // The previous form (a*f0/(1-f0) / |a|*(1-f0)/f0) is
            // INVALID for f0 > 0.5 with negative coefficients — it cut
            // off the true optimum on flugpl; brute-force validated on
            // random rows (see /tmp harness in the commit message).
            CutRow cut;
            double lhs = f0;
            std::map<int, double> coef;  // column -> accumulated coefficient
            for (auto& [j, a] : row_terms) {
                double sigma = (st[j] == AT_UPPER) ? -1.0 : 1.0;
                double z_orig = z[j] * E.sc.col[j];
                bool int_form =
                    cut_integ[j] &&
                    std::fabs(z_orig - std::floor(z_orig + 0.5)) < 1e-6;
                double pi;
                if (int_form) {
                    double fj = a - std::floor(a);
                    pi = (fj <= f0) ? fj : f0 * (1.0 - fj) / (1.0 - f0);
                } else {
                    pi = (a > 0) ? a : -a * f0 / (1.0 - f0);
                }
                if (std::fabs(pi) < 1e-9) continue;
                coef[j] += pi * sigma;
                lhs += pi * sigma * z_orig;
            }
            for (auto& [i, a] : slack_terms) {
                // slack i is continuous with original value = row activity
                long sv = E.nx + i;
                double sigma = (st[sv] == AT_UPPER) ? -1.0 : 1.0;
                double pi = (a > 0) ? a : -a * f0 / (1.0 - f0);
                if (std::fabs(pi) < 1e-9) continue;
                double act = 0.0;
                for (int p = work.ap[i]; p < work.ap[i + 1]; ++p) {
                    int j = work.ai[p];
                    coef[j] += pi * sigma * work.ax[p];
                    act += work.ax[p] * sol.x[j];
                }
                lhs += pi * sigma * act;
            }
            for (auto& [j, c] : coef)
                if (std::fabs(c) > 1e-9) cut.coeffs.emplace_back(j, c);
            cut.lhs = lhs;
            // dynamism guard: reject wild cuts
            double amax = 0.0;
            for (auto& [j, c] : cut.coeffs)
                amax = std::max(amax, std::fabs(c));
            if (amax > 1e6 || cut.coeffs.empty()) { ++dbg_dyn; continue; }
            // density cap: a dense GMI row (common on set-covering models
            // where every tableau alpha is O(1)) multiplies the LP's nnz
            // and every BTRAN/FTRAN crawls — measured on mod010: 100
            // fully-dense rows made node LPs ~35x slower than the bound
            // closure was worth. Sparse cuts only.
            if ((int)cut.coeffs.size() > 100) { ++dbg_dense; continue; }
            ++dbg_emit;
            if (const char* dp = getenv("IGAOS_DUMP_CUTS")) {
                FILE* f = fopen(dp, "a");
                if (f) {
                    std::fprintf(f, "CUT nnz=%zu lhs=%.17g\n",
                                 cut.coeffs.size(), cut.lhs);
                    for (auto& [j, c] : cut.coeffs)
                        std::fprintf(f, "%d %.17g\n", j, c);
                    fclose(f);
                }
            }
            cuts_out->push_back(std::move(cut));
        }
        if (getenv("IGAOS_DEBUG_CUTS"))
            std::fprintf(stderr,
                         "[cuts] rows=%ld nonstruct=%ld nonint=%ld f0=%ld "
                         "art_unfixed=%ld empty=%ld dyn=%ld dense=%ld "
                         "emitted=%ld arts_basic=%ld\n",
                         dbg_rows, dbg_nonstruct, dbg_nonint, dbg_f0,
                         dbg_art, dbg_empty, dbg_dyn, dbg_dense, dbg_emit,
                         dbg_artcols);
    }
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

Solution solve(const io::Model& model, const Options& opt,
               const WarmStart* warm, WarmStart* warm_out,
               std::vector<CutRow>* cuts_out,
               const std::vector<unsigned char>* true_integ) {
    Solution sol = solve_impl(model, opt, warm, warm_out, cuts_out,
                              true_integ);
    // Presolve-reduced models occasionally walk the simplex into a
    // near-singular-basis cascade the recovery machinery can't escape
    // (observed: wood1p reduced 244x2594 -> 172x1802, sparse-LU pivot
    // floor rejects, rollback blocks exhaust -> "basis factorization
    // failed"). Honest fallback: an ERROR on a presolved cold solve is
    // retried on the unreduced model — a retry can only spend time,
    // never mask a wrong answer (Error is loud by construction).
    if (sol.status == Status::Error && opt.presolve && warm == nullptr &&
        warm_out == nullptr && cuts_out == nullptr && !model.has_quadratic()) {
        Options o2 = opt;
        o2.presolve = false;
        sol = solve_impl(model, o2, warm, warm_out, cuts_out, true_integ);
        if (sol.status != Status::Error)
            sol.message += " (recovered: presolve-off retry)";
    }
    return sol;
}

}  // namespace igaos::simplex
