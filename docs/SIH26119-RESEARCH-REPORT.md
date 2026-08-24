# SIH26119 Research Report: Indigenous GPU-Accelerated Optimization Solver

**Smart India Hackathon 2026 — Problem Statement S.No. 119**
*Report date: 24 August 2026*

---

## 1. Verified Problem Statement (from official mirror of sih.gov.in)

| Field | Value |
|---|---|
| **PS ID** | SIH26119 (S.No. 119) |
| **Title** | Indigenous GPU-Accelerated Optimization Solver (Sovereign Alternative to Xpress/CPLEX)* |
| **Organization** | MRPL (Mangalore Refinery and Petrochemicals Ltd) |
| **Department** | Mangalore Refinery and Petrochemicals Limited (MRPL) |
| **Category / Theme** | Software / Miscellaneous |
| **Idea deadline** | **20 September 2026** (ideas: 0/500 as of scrape) |
| **Schedule** | Ideas Aug 21–Sep 20 · Evaluation Sep 10–Oct 30 · Results Nov wk 1 · Mentoring Nov 10–30 · Grand Finale Dec 2026 (tentative) |

*\*The official title contains typos ("Express / CEPLEX") — it means FICO Xpress and IBM CPLEX.*

### What MRPL actually demands (full PS text summary)

- A **solver core**, not a modeling environment. Initial scope: **LP + MILP + QP**, with modular architecture extensible later to MIQP, NLP, MINLP.
- Algorithms explicitly named: revised simplex, interior-point methods, branch-and-bound, branch-and-cut, cutting planes, presolve, heuristics, advanced node selection.
- Sparse matrix techniques, efficient numerical linear algebra, multi-core parallelization, with **GPU acceleration where it provides measurable benefits**.
- **Hard constraint:** "It shall not be built upon any existing open source solver library but shall be built from scratch from mathematical foundation."
- Benchmarks: standard sets — **MIPLIB, Netlib, Mittelmann instances, QPLIB (QP)** — with solution quality and runtime compared against **at least one established commercial or open-source solver**; plus representative refinery scheduling, crude blending, production planning and supply-chain case studies from open literature.
- Target scale: industrial problems with **thousands to millions of variables/constraints**, including highly degenerate models, ill-conditioned matrices, weak LP relaxations, hard mixed-integer formulations.
- Deliverable: robust engine with **basic API or CLI**; polished GUI explicitly not required.
- Domain emphasis (MRPL's world): refinery scheduling, crude blending, process optimization, production planning, logistics, power system dispatch, transportation, supply chain management.

---

## 2. Buildability Verdict

**As literally stated (match Gurobi/CPLEX at millions of variables): NO — not for anyone.**
CPLEX/Gurobi each represent hundreds of person-years of algorithmic tuning. HiGHS, the strongest open-source solver, took academic teams decades of incremental work. No hackathon team closes that gap in one cycle — and evaluators know it.

**As realistically evaluable: YES — conditionally.**
The PS itself leaves the opening: *"optimal **or near-optimal** solutions"* and comparison against *at least one established solver*. A strong 6-person team can credibly deliver:

| Component | Feasible scope | Effort |
|---|---|---|
| LP — revised simplex (sparse, Harris ratio test, Devex pricing, scaling) | Optimal on most of Netlib | Hard but doable |
| LP — GPU first-order method (PDLP-style primal-dual hybrid gradient, from scratch in CUDA) | Large sparse LPs; genuine GPU speedup story at scale | Doable — sparse mat-vecs only, maps naturally to GPU |
| QP — ADMM/operator-splitting (implementing the published OSQP *algorithm*, not its code) | Medium QPs | Moderate |
| MILP — branch-and-bound + presolve + Gomory cuts + rounding/diving heuristics + **GPU feasibility-pump heuristic** | MIPLIB easy/medium subset, near-optimal gaps | The long pole |
| Benchmark harness + honest comparison tables vs HiGHS | Fully controllable | Easy |

### Critical insight from current literature (pitch ammunition)

- Sequential simplex does **not** parallelize well on GPUs — do not claim it will.
- What genuinely wins on GPU:
  - **First-order PDHG methods**: NVIDIA's cuPDLP beats CPU interior-point ~3× on a 757k-row power-dispatch LP (PSR/NVIDIA benchmark).
  - **Interior-point with GPU sparse direct solvers (cuDSS)**: >10× speedups at medium precision reported (arXiv 2508.16094).
  - **GPU primal heuristics for MIP**: the CHAP paper (2026 Land–Doig MIP Competition) shows a GPU-CPU hybrid heuristic framework finding feasible solutions on **47/50 instances vs Gurobi's 44** under a 5-minute limit, also beating cuOpt heuristics-only mode (43).
- The ownable student-team niche: **fast feasible/near-optimal solutions on huge instances within time limits** — not exact solves.

---

## 3. Architecture (all from scratch)

```
┌─ CLI + Python API (pybind11), MPS/LP file reader ─────────────┐
│ Presolve/postsolve (own implementation)                        │
│ ┌──────────────┬───────────────────┬────────────────────────┐ │
│ │ Simplex (CPU)│ PDLP first-order  │ ADMM-QP                │ │
│ │ sparse       │ (CUDA, from       │ (CPU/GPU)              │ │
│ │ revised      │ scratch)          │                        │ │
│ └──────┬───────┴─────────┬─────────┴───────────┬────────────┘ │
│        └──── Branch & Bound ← cuts, diving, GPU feas.-pump ──┘ │
│ Numerics layer: scaling, perturbation, bound flipping,         │
│ iterative refinement, mixed precision w/ double fallback       │
│ Benchmark harness: Netlib / MIPLIB / Mittelmann vs HiGHS       │
└────────────────────────────────────────────────────────────────┘
```

- **Language:** C++20 + CUDA C++; Python bindings for usability.
- **Allowed gray zone — resolve early:** CUDA math libraries (cuSPARSE/cuBLAS/cuDSS) are numerical primitives, not solver libraries. Get written clarification from MRPL/SPOC that using them does not violate "not built upon any existing open source solver library" (it shouldn't — otherwise you'd have to write your own BLAS too).
- **Numerical robustness kit** (this is what the PS really tests):
  - Geometric-mean scaling of rows/columns
  - Harris two-pass ratio test for stability
  - Cost shifting + perturbation against degeneracy
  - Regularized KKT systems for ill-conditioning
  - Tolerance ladder with exact/double-precision fallback
  - Iterative refinement of solutions

---

## 4. Plan Mapped to the SIH Calendar

| Phase | Window | Deliverables |
|---|---|---|
| **0. Lock-in** | Now → Sep 20 | Team of 6 (2 OR/algorithms, 2 CUDA/HPC, 1 systems/benchmarking, 1 pitch/math-writing); faculty OR mentor; working spike: dense simplex + toy B&B; idea PPT with benchmark strategy; submit idea |
| **1. Core engine** | Sep–Oct (internal hackathon + evaluation window) | Sparse MPS reader, revised simplex solving Netlib end-to-end, presolve, scaling; PDLP-CUDA prototype; B&B skeleton; benchmark harness v1 |
| **2. Differentiate** | Nov 10–30 (mentoring window) | Tuned GPU PDLP (restarts, adaptive steps), cuts + heuristics incl. GPU feasibility pump, robustness suite (degenerate/ill-conditioned instances), **crude-blending refinery demo model**, results tables |
| **3. Finale** | Dec 2026 (36–48 h) | Live benchmark runs on stage machine, gap-table dashboard, rehearsed pitch |

### SIH 2026 official schedule (for reference)

| Activity | Timeline |
|---|---|
| Launch, PS announcement, idea/solution submission | Aug 21 – Sep 20, 2026 |
| Evaluation | Sep 10 – Oct 30, 2026 |
| Result announcement | First week of November 2026 |
| Training/mentoring for finale teams | Nov 10 – Nov 30, 2026 |
| Grand Finale | December 2026 (tentative) |

---

## 5. Competition Landscape & Winning Strategy

- **Pickup will be low** — this PS filters out most teams (requires CUDA + numerical optimization literacy). Quality-per-team will be high, and MRPL evaluators wrote the PS knowing weak implementations fail; they will probe numerics hard.
- **Most competitors will fail in one of two ways:**
  1. Wrapping HiGHS/CBC behind a UI → direct violation of the from-scratch clause.
  2. Toy dense-simplex demos that collapse on real MIPLIB instances.
- **Your winning wedge:** algorithmic depth + measured honesty. Show the crossover chart:
  - Your CPU simplex matches HiGHS on small/medium instances.
  - Your from-scratch GPU first-order method pulls ahead on large instances.
  - Your GPU heuristics find good MILP incumbents fast under time limits.
  - Cite the Land–Doig 2026 competition result as proof this approach beats commercial solvers in time-boxed feasibility regimes.
- **SIH rubric fit:** Problem understanding 30% (refinery framing, degeneracy discussion) · Technical implementation 25% (real benchmarks) · Innovation 20% (GPU-first-order + GPU heuristics from scratch) · Feasibility/scalability 15% · Presentation 10%.

---

## 6. Top Risks

1. **Numerics rabbit holes** — budget them; ship tolerance ladders/fallbacks before polish.
2. **GPU claims scrutiny** — never show an unexplained GPU loss; state upfront that small problems lose to CPU and show where the crossover happens.
3. **"From scratch" ambiguity** — get written clarification on cuDSS/cuSPARSE/cuBLAS usage before building on them.
4. **Team skill gap** — without someone who has implemented simplex/B&B before, immediately cut QP and GPU-MILP-heuristic scope.

---

## 7. Bottom Line

Buildable as a scoped sovereign-solver **demonstration**: exact LP via your own simplex, scalable LP via your own GPU first-order method, respectable MILP via your own branch-and-bound + GPU heuristics, all benchmarked honestly against HiGHS. Not buildable as a true Gurobi competitor — and saying so explicitly in the pitch earns evaluator trust rather than costing points.
