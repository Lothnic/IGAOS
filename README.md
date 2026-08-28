# IGAOS — Indigenous GPU-Accelerated Optimization Solver

A sovereign LP / MILP / QP solver core built from mathematical foundations for
[SIH 2026 problem statement SIH26119](https://github.com/Lothnic/IGAOS/issues/1)
(MRPL): revised simplex (primal + dual), first-order GPU methods (PDHG),
branch-and-bound with Gomory cuts, and ADMM QP — no existing
optimization-solver library as a base.

- Wayfinder map: [issue #1](https://github.com/Lothnic/IGAOS/issues/1)
- Problem-statement analysis: [`docs/SIH26119-RESEARCH-REPORT.md`](docs/SIH26119-RESEARCH-REPORT.md)
- Research sheets: `docs/research/` (PDHG algorithm · simplex design · benchmark protocol · refinery cases)
- Dependency policy: [`docs/DEPENDENCIES.md`](docs/DEPENDENCIES.md) · vocabulary: [`CONTEXT.md`](CONTEXT.md)

## Layout

```
src/
  common/   shared types, numerics utilities
  linalg/   sparse/dense linear algebra; swappable CPU/GPU backends
  simplex/  revised simplex (primal + dual, eta updates, warm starts)
  pdhg/     first-order GPU LP engine
  milp/     branch-and-bound + Gomory cuts
  qp/       OSQP-style ADMM QP engine
  io/       MPS reader (LP/MILP/QP), solution writers
  api/      CLI + pybind11 surface
python/     Python bindings: igaos.solve() / igaos.read_mps()
benchmarks/ harness per docs/research/benchmark-protocol.md
tests/      assert-based engine smoke tests
```

## Build

```sh
cmake -S . -B build
cmake --build build
```

Builds CPU-only automatically when no CUDA toolchain is present (the PDHG
engine requires CUDA).

## Solve

```sh
$ ./build/src/api/igaos solve model.mps --engine auto --time-limit 60
{
  "instance": "model.mps",
  "status": "optimal",
  "objective": -464.7531429,
  ...
}
```

Engines: `auto | simplex | pdhg | milp | qp`. All four engine classes are
live and verified against pinned baselines — current scores: Netlib
52/64 exact vs HiGHS, MIPLIB starters 6/20 @1e-4, robustness suite
10/15 per-class gates, Haverly QP three-way verified. Details and
honest failure records: `docs/research/`.

Python:

```python
import igaos
sol = igaos.solve("model.mps", time_limit=60, engine="milp")
sol.status, sol.objective, sol.x
```

