#pragma once

// Shared engine dispatch — used by both the CLI and the pybind11
// bindings so the two surfaces cannot drift.

#include <string>

#include "igaos/options.h"
#include "igaos/solution.h"
#include "igaos/status.h"
#include "model.hpp"
#include "milp.hpp"
#include "simplex.hpp"
#ifdef IGAOS_HAS_PDHG
#include "pdhg.hpp"
#endif

namespace igaos {

// engine: "auto" | "simplex" | "pdhg" | "milp" | "qp"
inline Solution solve_with_engine(const io::Model& model,
                                  const Options& opts,
                                  const std::string& engine) {
    Solution sol;
    if (engine == "qp") {
        sol.status = Status::Error;
        sol.message = "engine not yet wired: qp";
    } else if (engine == "milp") {
        sol = milp::solve(model, opts);
    } else if (engine == "simplex") {
        sol = simplex::solve(model, opts);
    } else if (engine == "pdhg") {
#ifdef IGAOS_HAS_PDHG
        sol = pdhg::solve(model, opts);
#else
        sol.status = Status::Error;
        sol.message = "built without CUDA: PDHG engine unavailable";
#endif
    } else {
        sol = simplex::solve(model, opts);
        if (sol.status == Status::Error)
            sol.message += "; auto fallback failed";
    }
    return sol;
}

}  // namespace igaos
