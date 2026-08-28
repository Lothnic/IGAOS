// pybind11 surface — the shape documented in python/README.md:
//   igaos.solve(path, ...) -> Solution
//   igaos.read_mps(path)   -> Model handle
// CLI parity: status is the same string `igaos solve --out` emits
// (igaos::status_name), and Model.stats mirrors the `igaos info` shape.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "engine.hpp"
#include "model.hpp"

namespace py = pybind11;
using igaos::io::Model;
using igaos::Options;
using igaos::Solution;

namespace {

const double INF = std::numeric_limits<double>::infinity();

Solution py_solve(const std::string& path, double time_limit,
                  double tolerance, long max_iterations, bool presolve,
                  long seed, const std::string& engine) {
    static const char* kEngines[] = {"auto",    "simplex", "pdhg",
                                     "milp",    "qp"};
    bool valid = false;
    for (const char* e : kEngines) valid |= engine == e;
    if (!valid)
        throw std::invalid_argument(
            "engine must be one of auto|simplex|pdhg|milp|qp, got '" +
            engine + "'");
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

// Same counting rules as cmd_info() in igaos_cli.cpp so the two surfaces
// cannot drift. rows_l excludes ranged rows (info prints L = nL - ranged).
py::dict model_stats(const Model& m) {
    int e = 0, l = 0, g = 0, ranged = 0;
    for (int i = 0; i < m.m; ++i) {
        if (m.rmin[i] == m.rmax[i]) ++e;
        else {
            bool fl = std::isfinite(m.rmin[i]), fu = std::isfinite(m.rmax[i]);
            if (fl) ++g;
            if (fu) ++l;
            if (fl && fu) ++ranged;
        }
    }
    int nfree = 0, nfixed = 0, nboxed = 0, nonesided = 0;
    for (int j = 0; j < m.n; ++j) {
        bool fl = std::isfinite(m.cl[j]), fu = std::isfinite(m.cu[j]);
        if (!fl && !fu) ++nfree;
        else if (fl && fu && m.cl[j] == m.cu[j]) ++nfixed;
        else if (fl && fu) ++nboxed;
        else ++nonesided;
    }
    int n_int = 0;
    for (unsigned char v : m.integ) n_int += v != 0;

    py::dict d;
    d["rows_e"] = e;
    d["rows_l"] = l - ranged > 0 ? l - ranged : 0;
    d["rows_g"] = g;
    d["rows_ranged"] = ranged;
    d["cols_free"] = nfree;
    d["cols_fixed"] = nfixed;
    d["cols_boxed"] = nboxed;
    d["cols_one_sided"] = nonesided;
    d["n_int"] = n_int;
    d["n_fr_parsed"] = m.n_fr_parsed;
    d["n_ranges_parsed"] = m.n_ranges_parsed;
    d["obj_const"] = m.obj_const;
    return d;
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
        .def("nnz", &Model::nnz)
        .def_property_readonly("stats", &model_stats)
        .def("__repr__", [](const Model& m) {
            return "<igaos.Model m=" + std::to_string(m.m) + " n=" +
                   std::to_string(m.n) + " nnz=" +
                   std::to_string(m.nnz()) + ">";
        });

    py::class_<Solution>(m, "Solution")
        .def_property_readonly(
            "status", [](const Solution& s) { return igaos::status_name(s.status); })
        .def_property_readonly(
            "status_enum", [](const Solution& s) { return s.status; })
        .def_readonly("objective", &Solution::objective)
        .def_readonly("x", &Solution::x)
        .def_readonly("y", &Solution::y)
        .def_readonly("row_activity", &Solution::row_activity)
        .def_readonly("pinf", &Solution::pinf)
        .def_readonly("dinf", &Solution::dinf)
        .def_readonly("rel_gap", &Solution::rel_gap)
        .def_readonly("iterations", &Solution::iterations)
        .def_readonly("solve_time_ms", &Solution::solve_time_ms)
        .def_readonly("message", &Solution::message)
        .def("__repr__", [](const Solution& s) {
            return "<igaos.Solution status='" +
                   std::string(igaos::status_name(s.status)) +
                   "' objective=" + std::to_string(s.objective) + ">";
        });

    py::enum_<igaos::Status>(m, "Status")
        .value("Optimal", igaos::Status::Optimal)
        .value("NearOptimal", igaos::Status::NearOptimal)
        .value("Feasible", igaos::Status::Feasible)
        .value("Infeasible", igaos::Status::Infeasible)
        .value("Unbounded", igaos::Status::Unbounded)
        .value("IterationLimit", igaos::Status::IterationLimit)
        .value("TimeLimit", igaos::Status::TimeLimit)
        .value("Error", igaos::Status::Error);

#ifdef IGAOS_HAS_PDHG
    m.attr("has_pdhg") = true;
#else
    m.attr("has_pdhg") = false;
#endif

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
