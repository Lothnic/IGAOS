#pragma once

#include "igaos/options.h"
#include "igaos/solution.h"
#include "model.hpp"

namespace igaos::pdhg {

Solution solve(const io::Model& model, const Options& options);

}
