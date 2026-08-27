// Reader + engine smoke tests: synthetic MPS fixtures covering the edge
// cases the robustness suite demands (negative UP, MI, FR, RANGES,
// objective constant, QUADOBJ) plus tiny known-optimum LP/MILP/QP solves.
// Assert-based, no framework — same pattern as linalg_smoke.cpp.
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>

#include "model.hpp"
#include "milp.hpp"
#include "qp.hpp"
#include "simplex.hpp"

using namespace igaos;
using io::Model;
static const double INF = std::numeric_limits<double>::infinity();

static const char* TMP = "/tmp/igaos_test.mps";

static Model from_string(const std::string& mps) {
    { std::ofstream f(TMP); f << mps; }
    Model m;
    io::parse_mps(TMP, m);
    return m;
}

static void test_negative_up() {
    // UP with negative value on a default-lo-0 column: lb becomes -inf
    // (MPS convention). x is free below, min x with x >= -5 via row.
    Model m = from_string(
        "NAME T1\nROWS\n N obj\n G r1\nCOLUMNS\n x obj 1.0 r1 1.0\n"
        "RHS\n RHS1 r1 -5.0\nBOUNDS\n UP BND x -3.0\nENDATA\n");
    assert(m.cl[0] == -INF);
    assert(m.cu[0] == -3.0);
    Solution s = simplex::solve(m, {});
    assert(s.status == Status::Optimal);
    assert(std::abs(s.objective - (-5.0)) < 1e-8);
    std::puts("negative-UP: ok");
}

static void test_mi_bound() {
    // MI: lb = -inf. min x s.t. x >= -2 -> -2
    Model m = from_string(
        "NAME T2\nROWS\n N obj\n G r1\nCOLUMNS\n x obj 1.0 r1 1.0\n"
        "RHS\n RHS1 r1 -2.0\nBOUNDS\n MI BND x\nENDATA\n");
    assert(m.cl[0] == -INF && m.cu[0] == INF);
    Solution s = simplex::solve(m, {});
    assert(s.status == Status::Optimal);
    assert(std::abs(s.objective - (-2.0)) < 1e-8);
    std::puts("MI-bound: ok");
}

static void test_ranges_and_objconst() {
    // RANGES on G row: [b, b+|r|]; objective constant -7.113 via RHS on N
    Model m = from_string(
        "NAME T3\nROWS\n N obj\n G r1\nCOLUMNS\n x obj 1.0 r1 1.0\n"
        "RHS\n RHS1 obj -7.113 r1 2.0\nRANGES\n RNG r1 3.0\nENDATA\n");
    assert(m.n_ranges_parsed == 1);
    assert(std::abs(m.obj_const - 7.113) < 1e-12);  // negated per MPS
    assert(m.rmin[0] == 2.0 && m.rmax[0] == 5.0);
    Solution s = simplex::solve(m, {});
    assert(s.status == Status::Optimal);
    assert(std::abs(s.objective - (2.0 + 7.113)) < 1e-8);
    std::puts("RANGES + obj-const: ok");
}

static void test_free_column() {
    // FR parsed; y's freedom bounded by a row (otherwise genuinely
    // unbounded): min x + 2y, x + y = 4, y >= 1, x >= 0 -> y=1, x=3, obj 5
    Model m = from_string(
        "NAME T4\nROWS\n N obj\n E r1\n G r2\nCOLUMNS\n x obj 1.0 r1 1.0\n"
        " x r2 0.0\n y obj 2.0 r1 1.0\n y r2 1.0\nRHS\n RHS1 r1 4.0 r2 1.0\n"
        "BOUNDS\n FR BND y\nENDATA\n");
    assert(m.n_fr_parsed == 1);
    Solution s = simplex::solve(m, {});
    assert(s.status == Status::Optimal);
    assert(std::abs(s.objective - 5.0) < 1e-8);
    std::puts("FR-column: ok");
}

static void test_infeasible() {
    Model m = from_string(
        "NAME T5\nROWS\n N obj\n G r1\n L r2\nCOLUMNS\n x obj 1.0 r1 1.0\n"
        " x r2 1.0\nRHS\n RHS1 r1 10.0 r2 5.0\nENDATA\n");
    Solution s = simplex::solve(m, {});
    assert(s.status == Status::Infeasible);
    std::puts("infeasible: ok");
}

static void test_milp_toy() {
    // max x (min -x) s.t. 2x <= 5, x integer -> x = 2
    Model m = from_string(
        "NAME T6\nROWS\n N obj\n L r1\nCOLUMNS\n MARK0000 'MARKER' 'INTORG'\n"
        " x obj -1.0 r1 2.0\n MARK0001 'MARKER' 'INTEND'\nRHS\n RHS1 r1 5.0\n"
        "ENDATA\n");
    assert(m.integ[0] == 1);
    Solution s = milp::solve(m, {});
    assert(s.status == Status::Optimal);
    assert(std::abs(s.objective - (-2.0)) < 1e-8);
    std::puts("milp-toy: ok");
}

static void test_qp() {
    // min x^2 + x s.t. x >= 1: analytic optimum x = 1, obj = 2
    Model m = from_string(
        "NAME T7\nROWS\n N obj\n G r1\nCOLUMNS\n x obj 1.0 r1 1.0\n"
        "RHS\n RHS1 r1 1.0\nQUADOBJ\n x x 2.0\nENDATA\n");
    assert(m.q_v.size() == 1);
    Options o;
    o.tolerance = 1e-6;
    Solution s = qp::solve(m, o);
    assert(s.status == Status::Optimal);
    assert(std::abs(s.objective - 2.0) < 1e-3);
    std::puts("qp-toy: ok");
}

int main() {
    test_negative_up();
    test_mi_bound();
    test_ranges_and_objconst();
    test_free_column();
    test_infeasible();
    test_milp_toy();
    test_qp();
    std::puts("engine_smoke: ALL PASS");
    return 0;
}
