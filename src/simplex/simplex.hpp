#pragma once

#include "igaos/options.h"
#include "igaos/solution.h"
#include "model.hpp"

namespace igaos::simplex {

// Warm start: final basis + nonbasic states from a previous solve of a
// model with the SAME matrix (B&B child nodes differ only in bounds).
// Variable indexing is the engine's: structural [0, n), slacks [n, n+m).
struct WarmStart {
    std::vector<int> basis;                 // size m: basic var per slot
    std::vector<unsigned char> nb_state;    // size n+m: nonbasic states
};

Solution solve(const io::Model& model, const Options& options,
               const WarmStart* warm = nullptr,
               WarmStart* warm_out = nullptr);

}  // namespace igaos::simplex
