// pybind11 surface — the shape documented in python/README.md:
//   igaos.solve(path, ...) -> Solution
//   igaos.read_mps(path)   -> Model handle
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "engine.hpp"
#include "model.hpp"

namespace py = pybind11;
using igaos::io::Model;
using igaos::Options;
using igaos::Solution;

namespace {

Solution py_solve(const std::string& path, double time_limit,
                  double tolerance, long max_iterations, bool presolve,
                  long seed, const std::string& engine) {
    Model mdl;
    igaos::io::parse_mps(path, mdl);
    Options opts;
    opts.time_limit_s = time_limit;
    opts.tolerance = tolerance;
    opts.max_iterations = max_iterations;
    opts.presolve = presolve;
    opts.seed = static_cast<unsigned>(seed);
    return igaos::solve_with_engine(mdl, opts, engine);
}

}  // namespace

PYBIND11_MODULE(igaos, m) {
    m.doc() = "IGAOS — sovereign LP/MILP solver core (SIH26119)";

    py::class_<Model>(m, "Model")
        .def_readonly("m", &Model::m)
        .def_readonly("n", &Model::n)
        .def_readonly("obj_const", &Model::obj_const)
        .def_readonly("n_fr_parsed", &Model::n_fr_parsed)
        .def_readonly("n_ranges_parsed", &Model::n_ranges_parsed)
        .def("nnz", &Model::nnz);

    py::class_<Solution>(m, "Solution")
        .def_readonly("status", &Solution::status)
        .def_readonly("objective", &Solution::objective)
        .def_readonly("x", &Solution::x)
        .def_readonly("y", &Solution::y)
        .def_readonly("row_activity", &Solution::row_activity)
        .def_readonly("pinf", &Solution::pinf)
        .def_readonly("dinf", &Solution::dinf)
        .def_readonly("rel_gap", &Solution::rel_gap)
        .def_readonly("iterations", &Solution::iterations)
        .def_readonly("solve_time_ms", &Solution::solve_time_ms)
        .def_readonly("message", &Solution::message);

    py::enum_<igaos::Status>(m, "Status")
        .value("Optimal", igaos::Status::Optimal)
        .value("NearOptimal", igaos::Status::NearOptimal)
        .value("Feasible", igaos::Status::Feasible)
        .value("Infeasible", igaos::Status::Infeasible)
        .value("Unbounded", igaos::Status::Unbounded)
        .value("IterationLimit", igaos::Status::IterationLimit)
        .value("TimeLimit", igaos::Status::TimeLimit)
        .value("Error", igaos::Status::Error);

    m.def("solve", &py_solve, py::arg("path"),
          py::arg("time_limit") = 120.0, py::arg("tolerance") = 1e-4,
          py::arg("max_iterations") = 500000L, py::arg("presolve") = true,
          py::arg("seed") = 0L, py::arg("engine") = "auto");

    m.def("read_mps",
          [](const std::string& path) {
              auto mdl = std::make_unique<Model>();
              igaos::io::parse_mps(path, *mdl);
              return mdl;
          },
          py::arg("path"), py::return_value_policy::move);
}
