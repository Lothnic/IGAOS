# Idea-Submission Deck — Content & Claims Ledger

Structure locked in #22 (10 slides, SIH rubric weights). Every claim below
maps to a repo artifact — the ledger is the traceability requirement.
MIPLIB table cells marked **[MIPLIB]** pending the full 20-instance run.

## Slide 1 — Title

- PS SIH26119 (MRPL). Pitch: *a sovereign LP/MILP/QP solver core, built
  from mathematical foundations — no existing solver library as a base.*
- Team: 6 (CUDA + linear algebra background), this repo.

## Slide 2 — The dependency problem

- Refinery planning (crude blending, pooling) runs on foreign solvers
  (Gurobi/CPLEX): license cost, opacity, no sovereign control.
- PS explicitly forbids building on existing open-source solvers → the
  capability must be built from scratch.
- Artifact: `docs/SIH26119-RESEARCH-REPORT.md`

## Slide 3 — Why existing options fall short

- Commercial: black box, per-seat licensing, export/control risk.
- Open source (HiGHS/CBC/GLPK): mature but CPU-first; GPU first-order
  methods (PDLP class) are where the field is moving.
- Artifact: `docs/research/pdlp-algorithm-sheet.md`

## Slide 4 — Architecture

io (MPS reader) → presolve-lite → {revised simplex · GPU PDHG · B&B ·
ADMM-QP (P2)} behind one `Status/Solution/Options` interface; CLI +
pybind11 surface.
- Artifacts: `src/`, `include/igaos/`, `src/api/engine.hpp`, `python/`

## Slide 5 — Proof of life (verified numbers)

| Evidence | Result | Artifact |
|---|---|---|
| Netlib LP ladder (64 inst ≤700 rows) | **52/64 exact vs HiGHS** (81%) | `docs/research/sweep_results.csv`, `demo/regress_netlib.py` |
| False-unbounded bug class | fixed (pilot4 exact; perold honest error) | commit `4c49d09` |
| Infeasibility detection | klein1 + refinery **proven** infeasible | `docs/research/robustness_results.csv` |
| MPS edge cases | RANGES/FR/objective-constant correct (e226 exact) | `src/io/model.hpp` |
| Python parity | bindings == CLI JSON on 4 cases | `demo/test_bindings.py` |
| GPU PDHG spike set | **7/8 obj-match ≤1e-3 @16s** (afiro/kb2/adlittle gap-certified; grow22 documented basin boundary) | `spike/pdhg-spike/` |

## Slide 6 — Innovation

- First-order GPU LP (restarted adaptive PDHG: Halpern averaging,
  KKT-error restarts, residual-balanced steps) in a from-scratch
  sovereign core — every engine implemented from mathematical
  foundations: primal+dual revised simplex, product-form eta updates,
  Gomory MI cuts, OSQP-style ADMM QP.
- Honest crossover: GPU wins above ~1M-nnz; small LPs stay CPU.

## Slide 7 — Benchmark honesty

- Pinned HiGHS v1.15.1 baseline; Netlib + MIPLIB 2017 (v36 solufile)
  protocol; `solved@1e-4 / tight@1e-6` ladder.
- Robustness suite: 10/15 per-class PASS (degeneracy, ill-conditioning,
  RANGES, free columns, infeasibility).
- MIPLIB starter subset: **6/20 @1e-4** (flugpl 3.4s, air03 5.9s,
  khb05250 32.2s proven optimal; 10/20 honest feasible incumbents) —
  `docs/research/miplib_results.txt`. Target ≥12/20.
- Artifacts: `docs/research/benchmark-protocol.md`,
  `demo/bench_robustness.py`, `benchmarks/suites/robustness.yaml`

## Slide 8 — MRPL fit

- Williams refinery LP solved to published optimum (−211365.13).
- Haverly pooling at three honesty levels: L0 LP relaxation AND L1
  convex QP both solved; L1 verified three ways (KKT, HiGHS cross-check,
  L0-bound consistency: −496.00 / −2196.00 vs HiGHS identical).
- MRPL-domain infeasible LP (`refinery`) correctly proven infeasible.
- Artifacts: `demo/models/`, `docs/research/refinery-case-studies.md`

## Slide 9 — Roadmap with kill criteria

- Sep 20: submission (this deck). Oct 31: QP/ADMM kill gate (freeze at
  L0 if unverified). Nov 15: component kill ladder (QP → GPU story →
  simplex+B&B floor). Dec: finale.
- Artifact: issue #7 (MVP scope), #16 (go/no-go).

## Slide 10 — Team + ask

- Ask: mentorship on industrial LP/MILP instances; MRPL pilot data.
- Team roles: `docs/research/` research sheets + this repo's history.

## Demo video plan (2 min)

1. `igaos info` + `igaos solve` on a Netlib instance (JSON output).
2. Robustness suite scroll (`bench_robustness.py` output).
3. Refinery model solve + infeasible-refinery proof.
Script: `demo/run_demo.sh` covers 1 and 3 today.
