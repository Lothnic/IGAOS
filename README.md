# IGAOS — Indigenous GPU-Accelerated Optimization Solver

A sovereign LP / MILP / QP solver core built from mathematical foundations for
[SIH 2026 problem statement SIH26119](https://github.com/Lothnic/IGAOS/issues/1)
(MRPL): revised simplex, first-order GPU methods (PDHG), and branch-and-bound —
no existing optimization-solver library as a base.

- Wayfinder map: [issue #1](https://github.com/Lothnic/IGAOS/issues/1)
- Problem-statement analysis: [`docs/SIH26119-RESEARCH-REPORT.md`](docs/SIH26119-RESEARCH-REPORT.md)
- Research sheets: `docs/research/` (PDHG algorithm · simplex design · benchmark protocol · refinery cases)
- Dependency policy: [`docs/DEPENDENCIES.md`](docs/DEPENDENCIES.md) · vocabulary: [`CONTEXT.md`](CONTEXT.md)

## Layout

```
src/
  common/   shared types, numerics utilities
  linalg/   sparse/dense linear algebra; swappable CPU/GPU backends
  simplex/  revised simplex engine
  pdhg/     first-order GPU LP engine
  milp/     branch-and-bound + cuts + heuristics
  io/       MPS/LP readers, solution writers
  api/      CLI + pybind11 surface
python/     Python bindings (lands with the first solver target)
benchmarks/ harness per docs/research/benchmark-protocol.md
tests/
```

Component libraries are CMake interface targets wired along the chain above until real
sources land.

## Build

```sh
cmake -S . -B build
cmake --build build
```

Builds CPU-only automatically when no CUDA toolchain is present (the PDHG
engine requires CUDA).

## Solve

```sh
$ ./build/src/api/igaos solve model.mps --time-limit 60
{
  "instance": "model.mps",
  "status": "near-optimal",
  "objective": -464.772,
  "pinf": 1.76e-05,
  ...
}
```

First wired engine: PDHG (first-order GPU LP). Simplex, MILP and QP follow
their locked specs (#5/#19/#20 resolutions).
