#pragma once

namespace igaos {

struct Options {
    double time_limit_s = 120.0;
    double tolerance = 1e-4;
    long max_iterations = 500000;
    int seed = 0;
    bool presolve = true;
    int verbosity = 0;
};

}  // namespace igaos
