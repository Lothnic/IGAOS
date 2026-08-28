// Warm-start smoke tests: solve, snapshot the final basis via
// simplex::solve's warm_out, tighten a bound, re-solve warm — the warm
// objective must match a cold solve of the SAME tightened model. Guards
// the dual/basis warm machinery (zb / st-z desync class of bugs).
// Assert-based, no framework — same pattern as engine_smoke.cpp.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "model.hpp"
#include "simplex.hpp"

using namespace igaos;
using io::Model;

static const char* TMP = "/tmp/igaos_warm_test.mps";

static Model from_string(const std::string& mps) {
    { std::ofstream f(TMP); f << mps; }
    Model m;
    io::parse_mps(TMP, m);
    return m;
}

// Hand-checked LP: min -x1 - 2x2 s.t. x1 + x2 <= 4, 0 <= x <= 3.
// Optimum x = (1, 3), obj = -7. Tightening cu[x1] = 0 moves the optimum
// to x = (0, 3), obj = -6 — a real optimum change for the warm basis to
// absorb (x1 leaves the basis, x2 stays basic).
static void test_warm_bound_tightening() {
    Model m = from_string(
        "NAME W1\nROWS\n N obj\n L r1\nCOLUMNS\n x1 obj -1.0 r1 1.0\n"
        " x2 obj -2.0 r1 1.0\nRHS\n RHS1 r1 4.0\nBOUNDS\n"
        " UI B x1 3\n UI B x2 3\nENDATA\n");
    Solution cold1 = simplex::solve(m, {});
    assert(cold1.status == Status::Optimal);
    assert(std::abs(cold1.objective - (-7.0)) < 1e-8);

    simplex::WarmStart ws;
    Solution s1 = simplex::solve(m, {}, nullptr, &ws);
    assert(s1.status == Status::Optimal);
    assert(ws.basis.size() == (size_t)m.m);
    assert(!ws.nb_state.empty());

    m.cu[0] = 0.0;  // tighten x1 away from its optimal value
    Solution warm2 = simplex::solve(m, {}, &ws);
    Solution cold2 = simplex::solve(m, {});
    assert(warm2.status == Status::Optimal);
    assert(cold2.status == Status::Optimal);
    assert(std::abs(warm2.objective - cold2.objective) < 1e-8);
    assert(std::abs(warm2.objective - (-6.0)) < 1e-8);  // hand-checked
    std::printf("warm-tighten: cold=%g warm=%g ok\n", cold2.objective,
                warm2.objective);

    // Second hop from the SAME snapshot: tighten x2 too. Optimum is now
    // x = (0, 0), obj = 0 — both basics are gone, warm must survive it.
    m.cu[1] = 0.0;
    Solution warm3 = simplex::solve(m, {}, &ws);
    Solution cold3 = simplex::solve(m, {});
    assert(warm3.status == Status::Optimal);
    assert(std::abs(warm3.objective - cold3.objective) < 1e-8);
    std::printf("warm-tighten-2: cold=%g warm=%g ok\n", cold3.objective,
                warm3.objective);
}

// Real-basis scale: afiro. Fix a column that is basic at the optimum to
// its optimal value (the old optimum stays feasible, so the objective is
// unchanged) — warm and cold must agree with the original optimum.
static void test_warm_real_model() {
#ifdef IGAOS_SRC_DIR
    Model m;
    io::parse_mps(IGAOS_SRC_DIR "/demo/models/netlib_full/afiro.mps", m);
    Solution cold = simplex::solve(m, {});
    assert(cold.status == Status::Optimal);
    assert(std::abs(cold.objective - (-464.7531429)) < 1e-4);

    simplex::WarmStart ws;
    Solution s = simplex::solve(m, {}, nullptr, &ws);
    assert(s.status == Status::Optimal);

    // pick a column strictly between its bounds (basic at optimum) and
    // fix it there; a second round fixes a second such column
    for (int round = 0; round < 2; ++round) {
        int pick = -1;
        for (int j = 0; j < m.n && pick < 0; ++j)
            if (m.cl[j] < m.cu[j] && cold.x[j] > m.cl[j] + 1e-6 &&
                cold.x[j] < m.cu[j] - 1e-6)
                pick = j;
        if (pick < 0) break;  // no interior column left to fix
        double v = cold.x[pick];
        m.cl[pick] = m.cu[pick] = v;
        Solution warm = simplex::solve(m, {}, &ws);
        Solution cold2 = simplex::solve(m, {});
        assert(warm.status == Status::Optimal);
        assert(cold2.status == Status::Optimal);
        assert(std::abs(warm.objective - cold2.objective) <
               1e-6 * (1.0 + std::abs(cold2.objective)));
        assert(std::abs(warm.objective - cold.objective) <
               1e-6 * (1.0 + std::abs(cold.objective)));
        std::printf("warm-afiro round %d: fix col %d = %.6g -> obj %.9g "
                    "ok\n", round, pick, v, warm.objective);
    }
#else
    std::puts("warm-real-model: skipped (IGAOS_SRC_DIR not defined)");
#endif
}

int main() {
    test_warm_bound_tightening();
    test_warm_real_model();
    std::puts("warm_smoke: ALL PASS");
    return 0;
}
