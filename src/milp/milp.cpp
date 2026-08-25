#include "milp.hpp"
#include <limits>

#include "simplex.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
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
    long nodes = 0;

    struct Node {
        std::vector<double> cl, cu;
    };
    std::vector<Node> stack;
    Node root;
    root.cl = model.cl;
    root.cu = model.cu;
    stack.push_back(root);

    const double INT_TOL = 1e-6;

    while (!stack.empty()) {
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
        Node nd = stack.back();
        stack.pop_back();
        ++nodes;

        io::Model lp = model;
        lp.cl = nd.cl;
        lp.cu = nd.cu;
        for (int j : ivars) lp.integ[j] = 0;

        Solution r = simplex::solve(lp, opt);
        if (r.status != Status::Optimal && r.status != Status::NearOptimal &&
            r.status != Status::Feasible)
            continue;

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
            if (has_incumbent && r.objective > best_obj - 1e-9) continue;
            double fl = std::floor(r.x[frac_var]);
            Node down = nd, upn = nd;
            down.cu[frac_var] = std::min(nd.cu[frac_var], fl);
            upn.cl[frac_var] = std::max(nd.cl[frac_var], fl + 1.0);
            if (down.cl[frac_var] <= down.cu[frac_var])
                stack.push_back(down);
            if (upn.cl[frac_var] <= upn.cu[frac_var])
                stack.push_back(upn);
            continue;
        }

        if (!has_incumbent || r.objective < best_obj) {
            best_obj = r.objective;
            best_x = r.x;
            has_incumbent = true;
        }
    }

    if (sol.status == Status::TimeLimit ||
        sol.status == Status::IterationLimit) {
        if (has_incumbent) {
            sol.status = Status::Feasible;
            sol.objective = best_obj;
            sol.x = best_x;
            char b[128];
            std::snprintf(b, sizeof(b), "%ld nodes explored", nodes);
            sol.message = b;
        }
    } else if (has_incumbent) {
        sol.status = Status::Optimal;
        sol.objective = best_obj;
        sol.x = best_x;
        char b[128];
        std::snprintf(b, sizeof(b), "%ld nodes explored", nodes);
        sol.message = b;
    } else if (sol.status != Status::TimeLimit &&
               sol.status != Status::IterationLimit) {
        sol.status = Status::Infeasible;
        sol.message = "no integer-feasible solution found";
    }
    sol.iterations = nodes;
    sol.solve_time_ms = elapsed() * 1000.0;
    return sol;
}

}  // namespace igaos::milp
