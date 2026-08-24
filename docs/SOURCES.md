# Sources — SIH26119 Research

*Compiled 24 August 2026. All URLs verified accessible at compile time.*

---

## 1. Problem Statement (primary source)

- **SIH26119 full PS text** (official mirror, scraped from sih.gov.in 2026-08-21, CC-BY-4.0):
  https://github.com/vedantchalke36/sih-2026-problem-statements/blob/main/ps_2026/SIH26119.md
- Official SIH portal: https://www.sih.gov.in/
- SIH 2026 PS page: https://www.sih.gov.in/sih2026PS
- Full dataset mirror (all 226 PS, JSON/CSV): https://github.com/vedantchalke36/sih-2026-problem-statements

## 2. SIH 2026 Schedule, Themes & Evaluation Context

- Times Now — "Smart India Hackathon 2026: Over 220 Problem Statements Released" (launch, schedule table, 17 themes):
  https://www.timesnownews.com/education/smart-india-hackathon-2026-over-220-problem-statements-released-article-155949432
- Ministry of Education Innovation Cell — official launch recording of SIH 2026 (YouTube):
  https://www.youtube.com/watch?v=2RaZ7moTzoQ
- CodeHunters Academy — SIH 2026 PS Explorer (evaluator-lens scoring, verdicts; analysis by a 3-time official evaluator):
  https://www.codehuntersacademy.com/sih-2026-ps
- Neurosignal — "Smart India Hackathon 2026: Complete Guide to Win" (evaluation dynamics, finale format ~1000 teams, 48–72h on-site builds):
  https://neurosignal.tech/smart-india-hackathon-2026-complete-guide-win/
- Reskilll — "SIH 2026: 20 Project Ideas That Actually Win" (rubric weights: understanding 30%, technical 25%, innovation 20%, feasibility 15%, presentation 10%):
  https://blogs.reskilll.com/smart-india-hackathon-2026-20-project-ideas-that-actually-win/
- Reskilll — "SIH 2026 Problem Statements Released" early analysis:
  https://reskilll.com/blogs/sih-2026-problem-statements-released-categories-themes-early-analysis/

## 3. GPU-Accelerated Optimization — State of the Art (technical evidence)

- NVIDIA Technical Blog — "How cuOpt Accelerates Mixed Integer Optimization using Primal Heuristics" (GPU feasibility pump + domain propagation; new feasible solutions on previously unsolved MIPLIB instances; paper ref arXiv 2510.20499):
  https://developer.nvidia.com/blog/learn-how-nvidia-cuopt-accelerates-mixed-integer-optimization-using-primal-heuristics/
- arXiv 2508.16094 — "GPU Implementation of Second-Order Linear and Nonlinear Programming Solvers" (IPM + cuDSS, >10× speedups at medium precision; KKT formulation survey):
  https://arxiv.org/html/2508.16094
- CHAP: A Hybrid GPU-CPU Heuristic for MIP (arXiv 2605.05086) — beats Gurobi (47 vs 44 found) and cuOpt (43) under 5-min limit on 50-instance Land–Doig 2026 competition benchmark:
  https://arxiv.org/html/2605.05086v1
- PSR Energy × NVIDIA — "GPU-based algorithms for large-scale optimization" (PDLP vs HiGHS IPM on 757k-row dispatch LP: 464s vs 1400s ≈ 3×; small LPs still favor CPU):
  https://www.psr-inc.com/en/analytics-report/post/gpu-based-algorithms-for-large-scale-optimization/
- CuClarabel (ACM TOMPECS) — GPU interior-point conic solver with cuDSS, mixed-precision results:
  https://doi.org/10.1145/3815420
- HiOp porting to accelerators (arXiv 2605.13736) — sparse IPM realized on GPU-dominant systems via dense compression:
  https://arxiv.org/html/2605.13736v1

## 4. Baseline Solvers & Benchmarks (comparison targets)

- HiGHS documentation — solvers overview (simplex/IPM/PDLP incl. GPU PDLP via cuPDLP-C; MIP branch-and-cut; multithreaded MIP prototype Feb 2026):
  https://ergo-code.github.io/HiGHS/stable/solvers/
- HiGHS GPU acceleration guide:
  https://ergo-code.github.io/HiGHS/stable/guide/gpu/
- Benchmark libraries referenced in the PS:
  - MIPLIB: https://miplib.zib.de/
  - Netlib LP: https://netlib.org/lp/
  - Mittelmann benchmarks: http://plato.asu.edu/bench.html
  - QPLIB: http://qplib.zib.de/

## 5. Prior Art: Student GPU Simplex Implementations

- alaintis/gpu-accelerated-simplex-solver (C++/CUDA revised simplex, GPU residency, Sherman–Morrison basis updates, CUDA Graphs, Harris ratio test, validated on Netlib vs HiGHS) — proof that the target architecture is student-buildable:
  https://github.com/alaintis/gpu-accelerated-simplex-solver

## 6. Solver Ecosystem Reference (context only — NOT to be used as base per PS constraint)

- COIN-OR CBC, GLPK, SCIP — named in the PS as existing open-source alternatives that still lag commercial solvers.
- CPMpy (solver-agnostic modeling layer listing commercial + open-source solvers):
  https://github.com/cpmpy/cpmpy

---

### Note on the "from scratch" clause

The PS forbids building **upon any existing open source solver library**. CUDA numerical primitives (cuSPARSE/cuBLAS/cuDSS) are math libraries, not solvers — using them is analogous to using BLAS. **Action item:** obtain written confirmation from MRPL/the SPOC before relying on them.
