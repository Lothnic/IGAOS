#include "milp.hpp"
#include <limits>

#include "simplex.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <tuple>
#include <vector>

namespace igaos::milp {

using std::numeric_limits;
constexpr double INF = std::numeric_limits<double>::infinity();

Solution solve(const io::Model& model, const Options& opt) {
    Solution sol;
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             t0)
            .count();
    };

    std::vector<int> ivars;
    for (int j = 0; j < model.n; ++j)
        if (model.integ[j]) ivars.push_back(j);
    if (ivars.empty()) {
        sol.status = Status::Error;
        sol.message = "no integer variables in model";
        return sol;
    }

    double best_obj = INF;
    std::vector<double> best_x;
    bool has_incumbent = false;
    // heuristic incumbents prune the tree but do NOT flip the pop rule to
    // best-bound: DFS diving is what actually finds good incumbents on
    // weak-bound instances (gt2/ran13x13). Only a tree-proven incumbent
    // (integral leaf) switches to best-bound for the proof.
    bool tree_incumbent = false;
    long nodes = 0, solves = 0, lp_failures = 0, lp_infeasible = 0;

    // best-bound + dive (#19): before any incumbent exists we dive
    // (deepest node) to reach integer feasibility fast; afterwards the
    // best-bound node is popped to close the proof. Nodes live in a pool;
    // the dive stack and bound multiset hold indices, dead entries are
    // skipped lazily.
    struct Node {
        std::vector<double> cl, cu;
        double bound;
        int depth;
        bool alive = true;
        simplex::WarmStart warm;      // parent's final basis (empty = cold)
        bool warm_valid = false;
        int bvar = -1;                // var this node branched on (parent's)
        bool bdown = false;           // which side of the parent branch
        double bfrac = 0.0;           // parent's fractionality of bvar
    };
    std::vector<Node> pool;
    std::vector<int> dive;  // LIFO stack of pool indices
    // ponytail: index multiset bound-major — fine to 1e5 nodes; a pairing
    // heap if node counts grow 100x
    auto bound_worse = [&](int a, int b) {
        if (pool[a].bound != pool[b].bound) return pool[a].bound > pool[b].bound;
        return pool[a].depth < pool[b].depth;
    };
    std::multiset<int, decltype(bound_worse)> by_bound(bound_worse);

    auto push_node = [&](Node&& nd) {
        pool.push_back(std::move(nd));
        int idx = (int)pool.size() - 1;
        dive.push_back(idx);
        by_bound.insert(idx);
    };

    Node root;
    root.cl = model.cl;
    root.cu = model.cu;
    root.bound = -INF;
    root.depth = 0;
    push_node(std::move(root));

    const double INT_TOL = 1e-6;
    const double GAP_TOL = 1e-9;  // bound pruning epsilon
    const double FEAS_TOL = 1e-6;

    // Pseudocost branching (#20): accumulate per-var unit objective gains
    // from the tree's own child LP solves (free — no extra LPs). Selection
    // score is the classic product form (eps + f*U)(eps + (1-f)*D);
    // uninitialized vars default to 1.0 which degrades gracefully to
    // most-fractional. ponytail: no strong branching — add shallow-depth
    // probing LPs if the free variant stalls short of the target.
    struct Pseudo {
        double dsum = 0, usum = 0;
        long dn = 0, un = 0;
    };
    std::vector<Pseudo> pc(model.n);
    auto pc_down = [&](int j) {
        return pc[j].dn ? pc[j].dsum / pc[j].dn : 1.0;
    };
    auto pc_up = [&](int j) {
        return pc[j].un ? pc[j].usum / pc[j].un : 1.0;
    };

    // Incumbent acceptance gate (#20): a heuristic candidate only counts if
    // it satisfies EVERY original row/bound/integrality of the model as
    // parsed — never trust the heuristic LP's own feasibility claim.
    auto validate_vs_model = [&](const std::vector<double>& x) {
        for (int j = 0; j < model.n; ++j) {
            if (x[j] < model.cl[j] - FEAS_TOL || x[j] > model.cu[j] + FEAS_TOL)
                return false;
            if (model.integ[j] &&
                std::fabs(x[j] - std::floor(x[j] + 0.5)) > INT_TOL)
                return false;
        }
        for (int i = 0; i < model.m; ++i) {
            double lhs = 0.0;
            for (int p = model.ap[i]; p < model.ap[i + 1]; ++p)
                lhs += model.ax[p] * x[model.ai[p]];
            if (lhs < model.rmin[i] - FEAS_TOL * (1.0 + std::fabs(model.rmin[i])))
                return false;
            if (lhs > model.rmax[i] + FEAS_TOL * (1.0 + std::fabs(model.rmax[i])))
                return false;
        }
        return true;
    };

    auto try_incumbent = [&](const std::vector<double>& x) {
        if (!validate_vs_model(x)) return -1;
        double obj = 0.0;
        for (int j = 0; j < model.n; ++j) obj += model.c[j] * x[j];
        if (has_incumbent && obj >= best_obj - GAP_TOL) return 0;
        best_obj = obj;
        best_x = x;
        has_incumbent = true;
        return 1;
    };

    // Root cut loop (#19 Rung-2): generate Gomory MI cuts at the root LP
    // optimum, add them as rows, re-solve warm — repeat while the bound
    // improves. Cut rows are appended to root_lp and every child builds
    // on root_lp (bounds-only deltas), so cuts propagate tree-wide.
    io::Model root_lp = model;
    for (int j : ivars) root_lp.integ[j] = 0;
    simplex::WarmStart root_ws;  // last root-round basis (seeds root node)
    if (std::getenv("IGAOS_NO_CUTS") == nullptr) {
        simplex::WarmStart& ws = root_ws;
        std::vector<simplex::CutRow> cuts;
        const int ROOT_ROUNDS = 8;
        double prev_bound = -INF;  // min: bound must INCREASE to continue
        double bound0 = INF;  // pre-cut root bound (cut-worthiness gate)
        for (int round = 0; round < ROOT_ROUNDS; ++round) {
            // root-cut budget: rounds after the mandatory root LP get at
            // most a quarter of the wall clock, and no round may outrun
            // the global budget (each simplex::solve runs its own clock)
            if (round > 0 && elapsed() > 0.25 * opt.time_limit_s) break;
            simplex::WarmStart ws_next;
            Options ropt = opt;
            ropt.time_limit_s =
                std::max(1.0, opt.time_limit_s - elapsed());
            Solution r = simplex::solve(root_lp, ropt,
                                        ws.basis.empty() ? nullptr : &ws,
                                        &ws_next, &cuts, &model.integ);
            if (r.status != Status::Optimal) {
                if (std::getenv("IGAOS_DEBUG_CUTS"))
                    std::fprintf(stderr, "[cuts] round=%d status=%d msg=%s\n",
                                 round, (int)r.status, r.message.c_str());
                break;
            }
            if (std::getenv("IGAOS_DEBUG_CUTS")) {
                std::fprintf(stderr, "[cuts] round=%d bound=%.6f cuts=%zu\n",
                             round, r.objective, cuts.size());
                for (const auto& c : cuts) {
                    double lhs_eval = 0.0;
                    for (auto& [j, co] : c.coeffs)
                        lhs_eval += co * r.x[j];
                    std::fprintf(stderr, "  cut: nnz=%zu lhs=%.6g eval=%.6g "
                                         "violation=%.3g\n",
                                 c.coeffs.size(), c.lhs, lhs_eval,
                                 lhs_eval - c.lhs);
                }
            }
            if (round == 0) bound0 = r.objective;
            // minimization: a bound INCREASE is the improvement — the
            // loop continues only while the cuts keep lifting the bound
            if (r.objective <= prev_bound + 1e-9 || cuts.empty()) break;
            prev_bound = r.objective;
            // append cut rows AND materialize them into the matrix —
            // the next round solves the augmented model, so the rebuild
            // must happen inside the loop
            int oldm = root_lp.m;
            std::vector<std::tuple<int, int, double>> all;
            for (int j = 0; j < root_lp.n; ++j)
                for (int p = root_lp.cp[j]; p < root_lp.cp[j + 1]; ++p)
                    all.emplace_back(root_lp.ci[p], j, root_lp.acx[p]);
            for (const auto& cut : cuts) {
                int newm = root_lp.m++;
                root_lp.rmin.push_back(cut.lhs);
                root_lp.rmax.push_back(INF);
                for (auto& [j, c] : cut.coeffs)
                    all.emplace_back(newm, j, c);
            }
            io::counting_sort(all, root_lp.m, root_lp.n, root_lp.ap,
                              root_lp.ai, root_lp.ax, root_lp.cp,
                              root_lp.ci, root_lp.acx);
            cuts.clear();
            // cut-row budget: cut rows at most double the model — on
            // small LPs (flugpl m=18) unbounded rounds piled 76 extra
            // rows on an 18-row LP and the tree drowned
            if (root_lp.m > 2 * model.m) break;
            // SOUND warm start across the appended rows: extend the
            // snapshot with the new rows' slack variables basic in their
            // own slots. The augmented basis matrix is block triangular
            // (old B + -I), the old optimum's reduced costs are
            // unchanged (cut slacks cost 0), so the restored basis is
            // DUAL feasible and primal infeasible exactly where the cuts
            // bite — the dual simplex repairs it through the standard
            // validated warm path (factor check, final orig-model guard).
            if (ws_next.basis.size() == (size_t)oldm) {
                for (int i = oldm; i < root_lp.m; ++i) {
                    ws_next.basis.push_back(root_lp.n + i);
                    ws_next.nb_state.push_back(2);  // engine BASIC
                }
                ws = std::move(ws_next);
            } else {
                ws = simplex::WarmStart();  // malformed: cold next round
            }
        }
        // cut-worthiness gate: keep cuts only when they closed the root
        // bound by >= 5% — otherwise the larger child LPs and the changed
        // branching trajectory cost more than the closure buys.
        // Measured on the MIPLIB slice: flugpl (0.8% closure, was PROVEN
        // optimal without cuts, finds no incumbent with them) and
        // ran13x13 (2.3%, slightly worse incumbent) must be dropped;
        // gt2 (53% closure, incumbent 77383 -> 32084) is the profile
        // cuts are for. ponytail: heuristic; per-instance cut selection
        // (aging, orthogonality) if the bar demands it
        bool cuts_worthwhile =
            std::isfinite(bound0) && std::isfinite(prev_bound) &&
            prev_bound > bound0 + 0.05 * (1.0 + std::fabs(bound0));
        if (!cuts_worthwhile && root_lp.m > model.m) {
            root_lp = model;               // drop the cut rows entirely
            for (int j : ivars) root_lp.integ[j] = 0;
        }
    }
    // root was moved-from by the first push_node — refill its bounds
    // (identical to model's unless worthwhile cuts are in root_lp)
    root.cl = root_lp.cl;
    root.cu = root_lp.cu;
    // seed the root node with the last round's basis when it matches the
    // final matrix (mismatch = cuts were dropped / rounds stopped early —
    // sizes fail the warm check and the node cold-solves)
    root.warm = std::move(root_ws);
    root.warm_valid =
        (int)root.warm.basis.size() == root_lp.m &&
        (int)root.warm.nb_state.size() == root_lp.n + root_lp.m;
    pool.clear();
    dive.clear();
    by_bound.clear();
    push_node(std::move(root));

    auto pop_node = [&]() -> int {
        auto try_dive = [&]() -> int {
            while (!dive.empty()) {
                int idx = dive.back();
                dive.pop_back();
                if (pool[idx].alive) return idx;
            }
            return -1;
        };
        auto try_bound = [&]() -> int {
            while (!by_bound.empty()) {
                int idx = *by_bound.begin();
                by_bound.erase(by_bound.begin());
                if (pool[idx].alive) return idx;
            }
            return -1;
        };
        // primary rule per #19 (dive before incumbent, best-bound after),
        // but never abandon live nodes in the other container. Mode switch
        // keys on tree_incumbent, not has_incumbent — see its comment.
        if (!tree_incumbent) {
            int i = try_dive();
            return i >= 0 ? i : try_bound();
        }
        int i = try_bound();
        return i >= 0 ? i : try_dive();
    };

    while (true) {
        if (elapsed() > opt.time_limit_s) {
            sol.status = Status::TimeLimit;
            sol.message = "wall-clock budget exhausted";
            break;
        }
        if ((long)(nodes + 1) > opt.max_iterations) {
            sol.status = Status::IterationLimit;
            sol.message = "node cap reached";
            break;
        }
        int idx = pop_node();
        if (idx < 0) break;
        Node nd = pool[idx];  // by value: push_node may reallocate the pool
        pool[idx].alive = false;
        by_bound.erase(idx);
        ++nodes;

        // prune on inherited bound before paying for an LP solve
        if (has_incumbent && nd.bound > best_obj - GAP_TOL) continue;

        io::Model lp = root_lp;  // includes root Gomory cuts
        lp.cl = nd.cl;
        lp.cu = nd.cu;
        for (int j : ivars) lp.integ[j] = 0;

        ++solves;
        simplex::WarmStart child_warm;
        bool warm_here = nd.warm_valid &&
                         std::getenv("IGAOS_NO_WARM") == nullptr;
        double t_node0 = elapsed();
        Solution r = simplex::solve(lp, opt, warm_here ? &nd.warm : nullptr,
                                    &child_warm);
        if (warm_here && (r.status == Status::Error ||
                          r.status == Status::Infeasible)) {
            // Warm starts land primal-infeasible (child bounds); the
            // art-based phase 1 cannot always repair basic bound
            // violations. Failures are loud (Error) or untrustworthy
            // (Infeasible from a bogus phase-1 optimum) — fall back to a
            // cold solve; warm still pays off on the nodes it completes.
            r = simplex::solve(lp, opt, nullptr, &child_warm);
        }
        if (std::getenv("IGAOS_DEBUG_NODES")) {
            double mv = 0.0;
            if ((int)r.x.size() == model.n) {
                for (int i = 0; i < model.m; ++i) {
                    double lhs = 0.0;
                    for (int p = model.ap[i]; p < model.ap[i + 1]; ++p)
                        lhs += model.ax[p] * r.x[model.ai[p]];
                    mv = std::max(mv, std::max(model.rmin[i] - lhs,
                                               lhs - model.rmax[i]));
                }
            }
            std::fprintf(stderr, "[node] %ld warm=%d st=%d obj=%.4f %.3fs "
                                "maxviol=%.3g\n",
                         nodes, (int)warm_here, (int)r.status, r.objective,
                         elapsed() - t_node0, mv);
        }
        bool child_warm_valid =
            r.status == Status::Optimal && !child_warm.basis.empty();
        // Root verdicts: the root LP is solved cold (the warm fallback
        // above re-solves cold on Error/Infeasible, so a surviving
        // Infeasible is a trustworthy phase-1 certificate; cuts in
        // root_lp are valid, so augmented-infeasible implies the MIP is
        // infeasible). Unbounded root LP can NOT certify an unbounded
        // MIP without a feasible integer point — report Error, honestly.
        if (nd.depth == 0 && nodes == 1) {
            if (r.status == Status::Infeasible) {
                sol.status = Status::Infeasible;
                sol.message = "root LP infeasible";
                sol.iterations = nodes;
                sol.solve_time_ms = elapsed() * 1000.0;
                return sol;
            }
            if (r.status == Status::Unbounded) {
                sol.status = Status::Error;
                sol.message = "root LP unbounded: MIP is unbounded or "
                              "infeasible (not certified)";
                sol.iterations = nodes;
                sol.solve_time_ms = elapsed() * 1000.0;
                return sol;
            }
        }
        if (r.status != Status::Optimal && r.status != Status::NearOptimal &&
            r.status != Status::Feasible) {
            ++lp_failures;
            // Infeasible survived the cold fallback above, so it is a
            // trustworthy phase-1 certificate for this subtree — count it
            // separately from engine Errors so the final status can tell
            // an infeasible MIP from a broken engine.
            if (r.status == Status::Infeasible) ++lp_infeasible;
            continue;
        }

        // LP bound from this node's own relaxation
        if (has_incumbent && r.objective > best_obj - GAP_TOL) continue;

        // free pseudocost update: this node's LP gain over its parent's
        // bound, per unit branch distance (bfrac is direction-relative:
        // f for a down child, 1-f for an up child)
        if (nd.bvar >= 0 && r.objective > nd.bound + GAP_TOL) {
            double unit =
                (r.objective - nd.bound) / std::max(nd.bfrac, 1e-6);
            Pseudo& p = pc[nd.bvar];
            if (nd.bdown) { p.dsum += unit; ++p.dn; }
            else          { p.usum += unit; ++p.un; }
        }

        int frac_var = -1;
        double frac_dist = INT_TOL;
        double best_score = -1.0;
        for (int j : ivars) {
            double x = r.x[j];
            double f = x - std::floor(x);  // up-distance to floor+1 is 1-f
            double fm = std::min(f, 1.0 - f);
            if (fm <= INT_TOL) continue;
            if (fm > frac_dist) frac_dist = fm;  // fallback: most fractional
            // product score: est. down gain = f*D, est. up gain = (1-f)*U
            double score = (1e-6 + f * pc_down(j)) *
                           (1e-6 + (1.0 - f) * pc_up(j));
            // fractionality factor so uninitialized vars (D=U=1, equal
            // scores) degrade to most-fractional
            score *= (0.01 + fm);
            if (score > best_score) {
                best_score = score;
                frac_var = j;
            }
        }

        if (frac_var >= 0) {
            // --- Diving heuristic (#20) -----------------------------------
            // Two cheap shots at a first incumbent before branching:
            //   (a) free rounding of the LP solution (no LP solves)
            //   (b) rounding dive: fix the most-fractional integer var to
            //       its nearest integer, warm-resolve, repeat until integral
            //       or infeasible (with one direction flip before giving
            //       up). Candidates are gated by validate_vs_model — only
            //       points feasible for the ORIGINAL rows/bounds count.
            // Runs only while NO incumbent exists: the root / early-DFS
            // dives are what buy the first feasible point cheaply. Once
            // any incumbent is in, every LP goes to the proof tree —
            // measured on flugpl, post-incumbent dives cost 3x wall time
            // for ~zero accepted improvements (bound-aborted dives at
            // every-16th node included).
            if (!has_incumbent) {
                // (a) pure rounding, three biases (nearest / floor /
                // ceil) — no LPs, validation-gated; covering models tend
                // to round up, packing models down
                {
                    const double bias[3] = {0.5, 1.0 - 1e-9, 1e-9};
                    for (int b = 0; b < 3; ++b) {
                        if (b > 0 && nodes > 100) break;  // nearest-only later on
                        std::vector<double> xr = r.x;
                        for (int j : ivars) {
                            double v = std::floor(xr[j] + bias[b]);
                            xr[j] = std::min(std::max(v, model.cl[j]),
                                            model.cu[j]);
                        }
                        if (try_incumbent(xr) == 1) break;
                    }
                }
                // (b) rounding dive — cheap warm LPs off the node basis.
                // At the root both directions run (nearest first: blend2's
                // incumbent 8.10 vs 40.1 guided; guided second: ran13x13
                // 3566 -> 3487); deeper nodes run nearest only
                if (child_warm_valid &&
                    elapsed() < opt.time_limit_s - 0.5) {
                  for (int dpass = 0; dpass < (nodes <= 1 ? 2 : 1); ++dpass) {
                    static long dbg_dives = 0, dbg_fail = 0, dbg_int = 0,
                                dbg_rej = 0, dbg_notbetter = 0, dbg_acc = 0;
                    ++dbg_dives;
                    std::vector<double> xd = r.x;
                    simplex::WarmStart ws = child_warm, ws_next;
                    const int DIVE_MAX = (int)ivars.size() + 4;
                    for (int step = 0; step < DIVE_MAX; ++step) {
                        if (elapsed() > opt.time_limit_s - 0.25) break;
                        int fv = -1;
                        double fd = INT_TOL;
                        for (int j : ivars) {
                            double f = std::fabs(xd[j] - std::floor(xd[j] + 0.5));
                            if (f > fd) { fd = f; fv = j; }
                        }
                        if (fv < 0) {
                            ++dbg_int;
                            {
                                int tr = try_incumbent(xd);
                                if (tr < 0) ++dbg_rej; else if (tr == 0) ++dbg_notbetter; else ++dbg_acc;
                            }
                            break;
                        }
                        double fl = std::floor(xd[fv]);
                        double fr = xd[fv] - fl;
                        // pass 1 (root only): pseudocost-guided — dive
                        // toward the smaller estimated objective gain
                        bool down = dpass == 0
                                        ? fr < 0.5
                                        : fr * pc_down(fv) <=
                                              (1.0 - fr) * pc_up(fv);
                        // last-feasible snapshot for the final-rounding
                        // fallback when every fix direction dies
                        std::vector<double> xlast = xd;
                        double old_cl = lp.cl[fv], old_cu = lp.cu[fv];
                        bool solved = false;
                        for (int dir = 0; dir < 2 && !solved; ++dir) {
                            bool go_down = dir == 0 ? down : !down;
                            lp.cl[fv] = old_cl;
                            lp.cu[fv] = old_cu;
                            if (go_down) lp.cu[fv] = fl;
                            else         lp.cl[fv] = fl + 1.0;
                            if (lp.cl[fv] > lp.cu[fv]) continue;
                            ++solves;
                            Solution dr = simplex::solve(lp, opt, &ws,
                                                         &ws_next);
                            if (dr.status == Status::Error ||
                                dr.status == Status::Infeasible) {
                                // warm start landed primal-infeasible (same
                                // phase-1 limitation as child nodes) — one
                                // cold retry before declaring this
                                // direction dead
                                ++solves;
                                dr = simplex::solve(lp, opt, nullptr,
                                                    &ws_next);
                            }
                            if (dr.status != Status::Optimal &&
                                dr.status != Status::NearOptimal &&
                                dr.status != Status::Feasible)
                                continue;
                            ws = ws_next;
                            xd = dr.x;
                            solved = true;
                            // bound-abort: a dive whose LP bound can no
                            // longer beat the incumbent is dead
                            if (has_incumbent &&
                                dr.objective > best_obj - GAP_TOL)
                                step = DIVE_MAX;
                        }
                        if (!solved) {
                            lp.cl[fv] = old_cl;
                            lp.cu[fv] = old_cu;
                            ++dbg_fail;
                            if (std::getenv("IGAOS_DEBUG_DIVE"))
                                std::fprintf(stderr,
                                             "[dive] abort step=%d var=%d "
                                             "val=%.4f\n",
                                             step, fv, xlast[fv]);
                            // final rounding: last feasible LP point with
                            // the remaining fractional vars rounded —
                            // validation gate rejects it if rows break
                            for (int j : ivars) {
                                double v = std::floor(xlast[j] + 0.5);
                                xlast[j] = std::min(std::max(v, model.cl[j]),
                                                    model.cu[j]);
                            }
                            ++dbg_int;
                            {
                                int tr = try_incumbent(xlast);
                                if (tr < 0) ++dbg_rej; else if (tr == 0) ++dbg_notbetter; else ++dbg_acc;
                            }
                            break;
                        }
                    }
                    if (std::getenv("IGAOS_DEBUG_DIVE") && dbg_dives % 200 == 1)
                        std::fprintf(stderr,
                                     "[dive] dives=%ld fail=%ld int=%ld "
                                     "rej=%ld notbetter=%ld acc=%ld\n",
                                     dbg_dives, dbg_fail, dbg_int, dbg_rej,
                                     dbg_notbetter, dbg_acc);
                  }
                }
            }
            // ------------------------------------------------------------
            double fl = std::floor(r.x[frac_var]);
            double f = r.x[frac_var] - fl;  // direction-relative distances
            Node down = nd, upn = nd;
            down.cu[frac_var] = std::min(nd.cu[frac_var], fl);
            down.bound = r.objective;
            upn.cl[frac_var] = std::max(nd.cl[frac_var], fl + 1.0);
            upn.bound = r.objective;
            down.depth = upn.depth = nd.depth + 1;
            down.warm = upn.warm = child_warm;
            down.warm_valid = upn.warm_valid = child_warm_valid;
            down.bvar = upn.bvar = frac_var;
            down.bdown = true;  upn.bdown = false;
            down.bfrac = f;     upn.bfrac = 1.0 - f;
            if (down.cl[frac_var] <= down.cu[frac_var])
                push_node(std::move(down));
            if (upn.cl[frac_var] <= upn.cu[frac_var])
                push_node(std::move(upn));
            continue;
        }

        if (!has_incumbent || r.objective < best_obj) {
            best_obj = r.objective;
            best_x = r.x;
            has_incumbent = true;
            tree_incumbent = true;
        }
    }

    char b[160];
    std::snprintf(b, sizeof(b), "%ld nodes / %ld LPs explored", nodes, solves);
    if (sol.status == Status::TimeLimit ||
        sol.status == Status::IterationLimit) {
        if (has_incumbent) {
            sol.status = Status::Feasible;
            sol.objective = best_obj + model.obj_const;
            sol.x = best_x;
            sol.message = b;
        }
    } else if (has_incumbent) {
        sol.status = Status::Optimal;
        sol.objective = best_obj + model.obj_const;
        sol.x = best_x;
        sol.message = b;
    } else if (lp_failures > 0) {
        if (lp_infeasible == lp_failures) {
            // every open subtree was LP-infeasible (cold certificates) —
            // the MIP has no integer-feasible point
            sol.status = Status::Infeasible;
            sol.message = std::string(b) + "; all open subtrees "
                          "LP-infeasible";
        } else {
            // LP engine failed on open nodes — not a proof of
            // infeasibility
            sol.status = Status::Error;
            sol.message = std::string(b) + "; " +
                          std::to_string(lp_failures) +
                          " node LPs unsolved";
        }
    } else {
        sol.status = Status::Infeasible;
        sol.message = "no integer-feasible solution found";
    }
    sol.iterations = nodes;
    sol.solve_time_ms = elapsed() * 1000.0;
    return sol;
}

}  // namespace igaos::milp
