#include "milp.hpp"
#include <limits>

#include "simplex.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
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
    long nodes = 0, solves = 0, lp_failures = 0;

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
        // but never abandon live nodes in the other container
        if (!has_incumbent) {
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

        io::Model lp = model;
        lp.cl = nd.cl;
        lp.cu = nd.cu;
        for (int j : ivars) lp.integ[j] = 0;

        ++solves;
        simplex::WarmStart child_warm;
        bool warm_here = nd.warm_valid &&
                         std::getenv("IGAOS_NO_WARM") == nullptr;
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
        bool child_warm_valid =
            r.status == Status::Optimal && !child_warm.basis.empty();
        if (r.status != Status::Optimal && r.status != Status::NearOptimal &&
            r.status != Status::Feasible) {
            ++lp_failures;
            continue;
        }

        // LP bound from this node's own relaxation
        if (has_incumbent && r.objective > best_obj - GAP_TOL) continue;

        int frac_var = -1;
        double frac_dist = INT_TOL;
        for (int j : ivars) {
            double f = std::fabs(r.x[j] - std::floor(r.x[j] + 0.5));
            if (f > frac_dist && f > INT_TOL) {
                frac_dist = f;
                frac_var = j;
            }
        }

        if (frac_var >= 0) {
            double fl = std::floor(r.x[frac_var]);
            Node down = nd, upn = nd;
            down.cu[frac_var] = std::min(nd.cu[frac_var], fl);
            down.bound = r.objective;
            upn.cl[frac_var] = std::max(nd.cl[frac_var], fl + 1.0);
            upn.bound = r.objective;
            down.depth = upn.depth = nd.depth + 1;
            down.warm = upn.warm = child_warm;
            down.warm_valid = upn.warm_valid = child_warm_valid;
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
        // LP engine failed on open nodes — not a proof of infeasibility
        sol.status = Status::Error;
        sol.message = std::string(b) + "; " + std::to_string(lp_failures) +
                      " node LPs unsolved";
    } else {
        sol.status = Status::Infeasible;
        sol.message = "no integer-feasible solution found";
    }
    sol.iterations = nodes;
    sol.solve_time_ms = elapsed() * 1000.0;
    return sol;
}

}  // namespace igaos::milp
