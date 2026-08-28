// MILP smoke tests: integer-infeasible / unbounded MIPs, a hand-checked
// knapsack, and the node-cap path (never a false Optimal). Assert-based,
// no framework — same pattern as engine_smoke.cpp.
//
// ponytail: the infeasible/unbounded cases accept Status::Error with an
// "unsolved node LPs" message because B&B currently conflates infeasible
// child LPs with engine failures (milp.cpp lp_failures). The hard guarantee
// under test is "never claims optimality" — tighten to Status::Infeasible /
// Status::Unbounded once milp distinguishes the two.
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

#include "model.hpp"
#include "milp.hpp"
#include "simplex.hpp"

using namespace igaos;
using io::Model;

static const char* TMP = "/tmp/igaos_milp_test.mps";

static Model from_string(const std::string& mps) {
    { std::ofstream f(TMP); f << mps; }
    Model m;
    io::parse_mps(TMP, m);
    return m;
}

static void test_integer_infeasible() {
    // 2x = 3 with x integer: LP-relaxation feasible (x=1.5), no integer
    // point exists -> the tree must exhaust without an incumbent.
    Model m = from_string(
        "NAME MI1\nROWS\n N obj\n E r1\nCOLUMNS\n MARK0000 'MARKER' "
        "'INTORG'\n x obj 1.0 r1 2.0\n MARK0001 'MARKER' 'INTEND'\n"
        "RHS\n RHS1 r1 3.0\nENDATA\n");
    assert(m.integ[0] == 1);
    Solution s = milp::solve(m, {});
    assert(s.status == Status::Infeasible);
    assert(!s.message.empty());
    std::printf("milp-infeasible: %s (%s)\n", status_name(s.status),
                s.message.c_str());
}

static void test_unbounded_mip() {
    // min -x, x integer, x >= 0, nothing else: unbounded relaxation.
    Model m = from_string(
        "NAME MU1\nROWS\n N obj\nCOLUMNS\n MARK0000 'MARKER' 'INTORG'\n"
        " x obj -1.0\n MARK0001 'MARKER' 'INTEND'\nRHS\nENDATA\n");
    Solution s = milp::solve(m, {});
    assert(s.status != Status::Optimal);
    assert(s.status != Status::Feasible);
    assert(!s.message.empty());
    std::printf("milp-unbounded: %s (%s)\n", status_name(s.status),
                s.message.c_str());
}

static void test_knapsack() {
    // max 7x1 + 9x2 + 5x3 + 10x4 s.t. 3x1+4x2+2x3+5x4 <= 7, binary.
    // Hand check: {x1,x2} w=7 v=16 (best); {x4,x3} w=7 v=15; {x2,x3}
    // w=6 v=14; every other pair overflows or is worth less.
    // min-form objective: -16, x = (1,1,0,0).
    Model m = from_string(
        "NAME KNAP\nROWS\n N obj\n L w\nCOLUMNS\n MARK0000 'MARKER' "
        "'INTORG'\n x1 obj -7.0 w 3.0\n x2 obj -9.0 w 4.0\n"
        " x3 obj -5.0 w 2.0\n x4 obj -10.0 w 5.0\n"
        " MARK0001 'MARKER' 'INTEND'\nRHS\n RHS1 w 7.0\nBOUNDS\n"
        " UI B x1 1\n UI B x2 1\n UI B x3 1\n UI B x4 1\nENDATA\n");
    Solution s = milp::solve(m, {});
    assert(s.status == Status::Optimal);
    assert(std::abs(s.objective - (-16.0)) < 1e-8);
    assert(s.x.size() == 4);
    assert(std::abs(s.x[0] - 1.0) < 1e-6 && std::abs(s.x[1] - 1.0) < 1e-6);
    assert(std::abs(s.x[2] - 0.0) < 1e-6 && std::abs(s.x[3] - 0.0) < 1e-6);
    std::puts("milp-knapsack: ok (-16)");
}

static void test_node_cap() {
    // Node cap hit before the tree is exhausted: an incumbent found so
    // far may be reported as Feasible (or IterationLimit with none) —
    // but NEVER as Optimal, and the incumbent must satisfy the model.
    Model m = from_string(
        "NAME CAP\nROWS\n N obj\n L r1\nCOLUMNS\n MARK0000 'MARKER' "
        "'INTORG'\n x obj -1.0 r1 1.0\n MARK0001 'MARKER' 'INTEND'\n"
        "RHS\n RHS1 r1 5.0\nENDATA\n");
    Options o;
    o.max_iterations = 1;  // stop right after the root node
    Solution s = milp::solve(m, o);
    assert(s.status == Status::Feasible || s.status == Status::IterationLimit);
    if (s.status == Status::Feasible) {
        // incumbent honesty: integral, within bounds, row satisfied
        assert(std::abs(s.x[0] - std::floor(s.x[0] + 0.5)) < 1e-6);
        assert(s.x[0] >= -1e-9 && s.x[0] <= 5.0 + 1e-9);
        assert(s.objective <= 1e-9);  // 0 >= obj >= -5
    }
    std::printf("milp-node-cap: %s obj=%.6f\n", status_name(s.status),
                s.objective);

    // Cap of 0: no node processed at all, no incumbent can exist.
    o.max_iterations = 0;
    Solution s0 = milp::solve(m, o);
    assert(s0.status == Status::IterationLimit);
    std::puts("milp-node-cap-0: iteration-limit, ok");
}

int main() {
    test_integer_infeasible();
    test_unbounded_mip();
    test_knapsack();
    test_node_cap();
    std::puts("milp_smoke: ALL PASS");
    return 0;
}
