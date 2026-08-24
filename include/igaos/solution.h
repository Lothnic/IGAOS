#pragma once

#include <string>
#include <vector>

#include "igaos/status.h"

namespace igaos {

struct Solution {
    Status status = Status::Error;
    double objective = 0.0;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> row_activity;
    double pinf = 0.0;
    double dinf = 0.0;
    double rel_gap = 0.0;
    long iterations = 0;
    double solve_time_ms = 0.0;
    std::string message;
};

}  // namespace igaos
