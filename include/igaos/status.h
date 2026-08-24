#pragma once

namespace igaos {

enum class Status {
    Optimal,
    NearOptimal,
    Feasible,
    Infeasible,
    Unbounded,
    IterationLimit,
    TimeLimit,
    Error,
};

inline const char* status_name(Status s) {
    switch (s) {
        case Status::Optimal: return "optimal";
        case Status::NearOptimal: return "near-optimal";
        case Status::Feasible: return "feasible";
        case Status::Infeasible: return "infeasible";
        case Status::Unbounded: return "unbounded";
        case Status::IterationLimit: return "iteration-limit";
        case Status::TimeLimit: return "time-limit";
        case Status::Error: return "error";
    }
    return "unknown";
}

}  // namespace igaos
