#pragma once

#include "igaos/options.h"
#include "igaos/solution.h"
#include "model.hpp"

namespace igaos::qp {

Solution solve(const io::Model& model, const Options& options);

}
