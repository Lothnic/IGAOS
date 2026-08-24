# PDLP / PDHG Algorithm Sheet — First-Order GPU LP Solver

**Ticket:** wayfinder research #4 (part of map #1) · **Repo:** Lothnic/IGAOS · **Date:** August 2026
**Purpose:** specification precise enough to implement a from-scratch first-order GPU LP solver against. All claims cite primary sources (papers + public source code).

---

## 1. Problem formulation

Solve the primal-dual LP pair ([Applegate et al. 2022, §2](https://arxiv.org/abs/2106.04756); [Lu et al. 2023, §2](https://arxiv.org/abs/2312.14832)):

```
min c'x                    max q'y + l'λ⁺ − u'λ⁻
s.t. Gx ≥ h                s.t. c − K'y = λ
     Ax = b                     y_{1:m₁} ≥ 0        (equality-block duals free)
     l ≤ x ≤ u                  λ ∈ Λ
```

with `K' = (G', A')`, `q = (h; b)`, and the bound-dual sign set

```
Λ_i = {0}   if l_i = −∞, u_i = +∞     ℝ⁻ if l_i = −∞ only
      ℝ⁺    if u_i = +∞ only          ℝ   otherwise (two-sided bound)
```

Equivalently the saddle-point problem, which is what PDHG actually solves:

```
min_{x∈X} max_{y∈Y}  L(x,y) := c'x − y'Kx + q'y
X = {x : l ≤ x ≤ u},   Y = ℝ^{m₁} × ℝ₊^{m₂}   (m₂ = inequality rows)
```

A saddle point of (2) recovers an optimal primal-dual pair by LP duality. Any general-form LP (including Netlib/MPS `RANGES`, free rows, etc.) converts into this form in presolve.

---

## 2. Baseline PDHG (Chambolle–Pock) iteration

PDHG applied to the saddle point ([Chambolle & Pock 2011](https://doi.org/10.1007/s10851-010-0251-1); specialized to LP in [Applegate et al. 2022, eq. (3)](https://arxiv.org/abs/2106.04756)):

```
x^{k+1} = proj_X( x^k − τ (c − K'y^k) )              # primal: gradient step + box projection
y^{k+1} = proj_Y( y^k + σ (q − K(2x^{k+1} − x^k)) )  # dual: over-relaxed step + sign projection
```

- `proj_X` = clamp to `[l, u]` elementwise. `proj_Y` = clamp equality-row duals to free, inequality duals to ≥ 0.
- The `2x^{k+1} − x^k` extrapolation is what makes it converge; plain Arrow–Hurwicz does not.
- Convergence guarantee: `τσ‖K‖²₂ ≤ 1` ⟺ `η ≤ 1/‖K‖₂` under the reparameterization below.
- Reparameterization used everywhere downstream: **`τ = η/ω` and `σ = ωη`**, where `η` is the step size and `ω > 0` the *primal weight* balancing primal vs dual progress. The weighted norm ‖z‖_ω = √(ω‖x‖² + ‖y‖²/ω) appears throughout the heuristics.

Per iteration this costs exactly **one SpMV with K and one with K′** plus O(n+m) vector ops. No factorization ever. Vanilla PDHG alone is *not* competitive on LP — baseline PDHG solved only 50/383 MIPLIB-relaxation instances at 1e−8 vs 283/383 for enhanced PDLP ([Applegate et al. 2022, §4.2](https://arxiv.org/abs/2106.04756)). The enhancements below are not optional; they are the solver.

---

## 3. Full algorithm (cuPDLP-style restarted adaptive PDHG)

This is the GPU variant: identical to PDLP except that restarts use a **KKT-error metric instead of the normalized duality gap** (the trust-region subproblem for the latter is sequential and GPU-hostile) — [Lu & Yang 2023, §3.4](https://arxiv.org/abs/2311.12180).

### Pseudocode

```
precondition: K̃ = D₁KD₂ via Ruiz equilibration (10 sweeps) then Pock–Chambolle α=1
presolve (optional but recommended): remove empty/fixed/duplicate rows & cols

initialize:
  x⁰ = y⁰ = 0
  η̂ = 1/‖K‖∞                       # initial step size (power-iteration-free bound)
  ω  = ‖c‖₂/‖q‖₂  if both norms > ε_zero else 1
  n ← 0 (outer/restart epoch), k ← 0 (total iterations)

repeat  ───────────────── outer loop (epochs) ─────────────────
  repeat  ── inner loop (until restart/termination), checked every ~40 iters ──
    1. PRIMAL STEP:  x^{t+1} = proj_X(x^t − (η/ω)(c − K'y^t))
    2. DUAL STEP:    y^{t+1} = proj_Y(y^t + ηω(q − K(2x^{t+1} − x^t)))
    3. ADAPTIVE η:   η̄ = ‖Δz‖²_ω / (2 Δy'KΔx),  Δz = z^{t+1} − z^t
         if η ≤ η̄: accept step
         else:      η ← min((1−(k+1)^{−0.3})·η̄, (1+(k+1)^{−0.6})·η) and retry step
       (guarantees the convergence inequality η ≤ ‖Δz‖²_ω/(2Δy'KΔx);
        provably η ≥ (1−o(1))/‖K‖₂ as k→∞ — near-maximal steps without estimating ‖K‖₂)
    4. AVERAGE:      S += η_k·z^{t+1};  W += η_k ;   z̄ = S/W
    5. every CHECK_INTERVAL iterations (40 in cuPDLP-C):
       - compute termination residuals for BOTH last iterate and average
       - restart candidate:  z_c = argmin_KKTω( z^{t+1}, z̄ )
       - RESTART to epoch start z_c if any of:
           (i)   KKT(z_c) ≤ 0.2 · KKT(z^{n,0})                      [sufficient decay]
           (ii)  KKT(z_c) ≤ 0.8 · KKT(z^{n,0}) and KKT increased
                 since previous check                                [stall safeguard]
           (iii) t ≥ 0.36·k                                          [artificial/long loop]
       - check INFEASIBILITY certificates from iterate differences
  until restart or termination
  ω ← exp( θ·log(Δy/Δx) + (1−θ)·log ω ),  θ = 0.5
      where Δx, Δy = ‖x^{n,0} − x^{n−1,0}‖, ‖y^{n,0} − y^{n−1,0}‖  (skip if either < ε_zero)
until termination or limits
return better of {last iterate, average iterate}
```

Sources: [Applegate et al. 2022, Algorithms 1–3, §3](https://arxiv.org/abs/2106.04756); GPU restart variant [Lu & Yang 2023, Alg. 1–2, §3.4](https://arxiv.org/abs/2311.12180); cuPDLP-C restatement [Lu et al. 2023, §2](https://arxiv.org/abs/2312.14832).

### Termination criteria (exact formulas)

Terminate when all three relative criteria hold, evaluated on the **original (unscaled) instance**:

```
primal feasibility:   ‖ Ax − b ; [h − Gx]⁺ ‖₂ ≤ ε(1 + ‖q‖₂)
dual feasibility:     ‖ c − K'y − λ ‖₂      ≤ ε(1 + ‖c‖₂)
relative duality gap: |q'y + l'λ⁺ − u'λ⁻ − c'x| ≤ ε(1 + |c'x| + |q'y + l'λ⁺ − u'λ⁻|)
bound duals:          λ = proj_Λ(c − K'y)   (never carried as iterates)
```

([cuPDLP-C paper §2](https://arxiv.org/abs/2312.14832); [cuPDLP.jl paper §4.1](https://arxiv.org/abs/2311.12180); identical in [cuPDLP-C README](https://github.com/COPT-Public/cuPDLP-C).)

### Weighted-average iterates — implementation detail

Verified in source ([cuPDLP-C `cupdlp_step.c`](https://github.com/COPT-Public/cuPDLP-C/blob/main/cupdlp/cupdlp_step.c), [cuPDLP.jl `saddle_point_gpu.jl`](https://github.com/jinwen-yang/cuPDLP.jl/blob/master/src/saddle_point_gpu.jl)):

- Maintain incrementally: `xSum += w_k·x^{k}`, `ySum += w_k·y^{k}`; per-epoch weight sums `W`. In cuPDLP-C the weight is the geometric mean of the current primal/dual steps, `w_k = √(τ_k σ_k)` (= η). Average = sum/W.
- Also maintain `K·x̄` and `K'·ȳ` as the same weighted sums of per-iterate products (`sum_primal_product`, `sum_dual_product`) so residual checks never re-multiply the average — one SpMV saved per check.
- Drift control: there is **no separate drift-correction step in the public papers/code**; drift is bounded by (a) resetting sums at every restart, (b) falling back avg ← current iterate whenever a numerical error is flagged or the count is zero (cuPDLP.jl main loop), and (c) testing both last and average iterates at termination so a degraded average can never be returned over a good last iterate. Implement all three.

### Primal weight update

At each restart only (frequent small updates destabilize PDHG): `ω ← exp(0.5·log(Δy/Δx) + 0.5·log ω_prev)` — exponential smoothing in log space balances weighted primal/dual distance traveled. Initialization `ω₀ = ‖c‖/‖q‖` gives scale invariance ([Applegate et al. 2022, §3.3](https://arxiv.org/abs/2106.04756)).

### Preconditioning (do this before anything else)

Diagonal scaling `K̃ = D₁KD₂` with transformed data (`c̃ = D₂c`, `(b̃,h̃) = D₁(b,h)`, bounds divided by `D₂`): default is **10 Ruiz equilibration sweeps (‖row‖∞ → 1) followed by Pock–Chambolle α=1 scaling** ((D₁)_jj = √‖K_j,·‖₁-ish row norm, (D₂)_ii = √‖K_·,i‖₁ col norm) ([Applegate et al. 2022, §3.5](https://arxiv.org/abs/2106.04756)). cuPDLP-C defaults confirm: `ifRuizScaling=true` (10×), `ifPcScaling=true`, `ifL2Scaling=false` ([README parameter table](https://github.com/COPT-Public/cuPDLP-C)). Presolve contributed an additional ~2–5× on MIPLIB relaxations in the cuPDLP-C tables (SGM10 10.28 → 5.43 with COPT presolve).

---

## 4. Parameter defaults table

Defaults from the [cuPDLP-C README](https://github.com/COPT-Public/cuPDLP-C) and papers; these are the numbers to ship v1 with.

| Parameter | Default | Notes |
|---|---|---|
| Step-size rule `eLineSearchMethod` | 2 = adaptive | alternatives: 0 fixed, 1 Malitsky–Pock |
| Initial step size η̂ | `1/‖K‖∞` | cheap upper bound; adaptive rule takes over |
| Restart rule `eRestartMethod` | 1 = KKT-based | 0 = none (never ship 0) |
| β_sufficient / β_necessary / β_artificial | 0.2 / 0.8 / 0.36 | GPU restart thresholds ([Lu & Yang 2023 §3.4](https://arxiv.org/abs/2311.12180)); CPU-PDLP normalized-gap variants are 0.9/0.1/0.5 |
| Primal-weight smoothing θ | 0.5 | update only at restarts |
| Ruiz iterations | 10 | then Pock–Chambolle α=1 |
| Check interval (residuals + restart) | every 40 iters (+ first 10 iters always) | `CUPDLP_RELEASE_INTERVAL=40` in [`cupdlp_defs.h`](https://github.com/COPT-Public/cuPDLP-C/blob/main/cupdlp/cupdlp_defs.h) |
| `dPrimalTol` / `dDualTol` / `dGapTol` | 1e−4 / 1e−4 / 1e−4 | papers use 1e−4 "moderate", 1e−8 "high quality" |
| Infeasibility tolerance `dFeasTol` | 1e−8 | |
| Iteration limit | INT_MAX (default); time limit 3600 s | always also cap wall time |
| Init point | zeros | all references use x=y=0 |

---

## 5. Why this maps to GPU (and where it wins/loses)

**The kernel inventory is exactly SpMV + vector ops.** Per iteration: two sparse mat-vecs (K·x, K'·y), elementwise axpy/clamp/dot products, and two norm reductions. No LU/Cholesky, no triangular solves, no pivoting — the three things that made every prior attempt at GPU simplex/IPM fail ([Lu & Yang 2023, §1](https://arxiv.org/abs/2311.12180): factorizations are memory-hungry and inherently sequential). Implementation facts from cuPDLP.jl ([§3](https://arxiv.org/abs/2311.12180)):

- Matrix stored CSR; SpMV via `cusparseSpMV()` with `CUSPARSE_SPMV_CSR_ALG2` (deterministic across runs); custom kernels assign **one thread per coordinate**; 1-D thread config suffices.
- Keep everything device-resident: exactly **two host↔device transfers per solve** (upload scaled instance, download solution).
- Restarts/residual reductions need grid-wide sync → CUDA cooperative groups grid sync ([NVIDIA cuOpt blog 2024](https://developer.nvidia.com/blog/accelerate-large-linear-programming-problems-with-nvidia-cuopt/)).
- The algorithm is memory-bandwidth-bound: performance scales directly with GPU HBM bandwidth (H100 FP64: 26 TFLOPS, 2 TB/s vs 16-core CPU 256 GFLOPS, 137 GB/s — a ~15× bandwidth ratio that shows up directly in solve time).

**Convergence behavior.** Restarted PDHG converges linearly on LP, matching the worst-case lower bound with restarts ([Applegate et al., Math. Program. 2023](https://doi.org/10.1007/s10107-022-01901-9)); Lu–Yang proved plain PDHG is linear on LP with a two-stage behavior (finite identification of active variables, then linear solve of a homogeneous system) which explains slow tails on near-degenerate instances ([2307.03664](https://arxiv.org/abs/2307.03664)). Practical consequences: fast early progress, tail-off at tight accuracy; degenerate/ill-conditioned instances are its weak spot.

**When it wins:** large sparse instances (≥ ~1M nonzeros keeps a GPU fed; crossover vs CPU IPM observed above roughly the 4-hour-horizon unit-commitment scale) and anything too big to factorize. **When it loses:** small/medium instances (kernel-launch overhead dominates; CPU IPM/simplex win outright), tight-accuracy demands (see §6), and numerically nasty instances (8/49 Mittelmann public instances time out even for NVIDIA's tuned implementation at 1e−4).

---

## 6. Accuracy: what PDHG delivers vs what MRPL's PS demands

- The SIH26119 PS explicitly accepts "**optimal or near-optimal** solutions" benchmarked for solution quality against at least one established solver — near-optimal is in-scope by the problem text itself ([SIH26119 report §1](../../SIH26119-RESEARCH-REPORT.md)).
- What PDHG delivers cheaply: relative tolerances **1e−4 (moderate) to 1e−6**, per NVIDIA's own guidance ("cuOpt can very quickly solve LPs to low accuracy, e.g., relative tolerance 1e−4 to 1e−6"; many applications needing 1e−8 should use barrier) ([NVIDIA cuOpt docs/blog](https://developer.nvidia.com/blog/solve-linear-programs-using-the-gpu-accelerated-barrier-method-in-nvidia-cuopt/)). Artelys/FICO note cuOpt PDLP ships a 4-digit default vs Xpress barrier's 6-digit ([Artelys news 2025](https://www.artelys.com/news/artelys-nvidia-fico-collaboration-power-system-optimisation/)).
- Cost of tightness (cuPDLP-C, MIPLIB-383 SGM10 seconds, COPT presolve): **5.43 @1e−4 → 18.53 @1e−8 (~3.4×)**, and instances solved drop 379 → 369 ([Table 1, Lu et al. 2023](https://arxiv.org/abs/2312.14832)). Plan a tolerance ladder: 1e−4 default, 1e−6 stretch, 1e−8 only if time remains.
- Reliability caveat to engineer for: moderate feasibility violations can be materially wrong — SCS at 1e−4 returned an objective off by **13%** on `irish-electricity`; this is precisely why PDLP's authors added high-accuracy machinery and why we must always post-check residuals/objective bounds before reporting ([PDLP journal paper, §1](https://arxiv.org/abs/2501.07018)).
- Bridge to exactness when needed: PDHG solutions can be polished/crossover'd to basic feasible solutions (OR-Tools and cuOpt both offer crossover from PDLP output); alternatively hand the PDHG solution to our own CPU simplex as a warm start. This covers MRPL workflows that need vertex solutions or exact optima on medium instances.

---

## 7. Practical numbers to calibrate against

| Setting | Result | Source |
|---|---|---|
| MIPLIB-2017 LP relaxations (383, nnz > 100k) | COPT(IPM) 383 solved @1e−8; cuPDLP-C 379 @1e−4 / 369 @1e−8; cuPDLP-C ≈ 2–4× slower than COPT @1e−4 | [Lu et al. 2023, Tables 1,4](https://arxiv.org/abs/2312.14832) |
| Mittelmann LP benchmark (49 public) | COPT 48 @1e−8; cuPDLP-C 46 @1e−4; cuOpt PDLP times out on 8/49 @1e−4 (1 h limit) | [Lu et al. 2023 Table 1](https://arxiv.org/abs/2312.14832); [NVIDIA blog 2024](https://developer.nvidia.com/blog/accelerate-large-linear-programming-problems-with-nvidia-cuopt/) |
| Huge instance `zib03`: 19.7M × 29.1M, 104M nnz | cuPDLP-C **916 s** on H100 vs COPT **16.5 h** (CPLEX needed 139 days in 2009) | [Lu et al. 2023, Table 2](https://arxiv.org/abs/2312.14832) |
| PageRank LP, 10M nodes (80M nnz) | 44 s on H100; IPM/simplex cannot run at this size | [Lu et al. 2023, Table 2](https://arxiv.org/abs/2312.14832) |
| Supply chain inv-60: 16.2M × 14.5M, 60M nnz | 11,102 s @ gap 1e−5 | [Lu et al. 2023, Table 2](https://arxiv.org/abs/2312.14832) |
| QAP LP relaxations (Adams–Johnson, n=50) | cuPDLP-C solves all 5 in ≤ 437 s; COPT simplex+barrier fail entirely | [Lu et al. 2023, §3.3](https://arxiv.org/abs/2312.14832) |
| 11 giant instances, 125M–**6.3B** nnz (PDLP journal) | CPU PDLP solves 8/11 to 1% gap in ≤ 6 days single machine; Gurobi barrier exceeds 1 TB RAM on 8; simplex solves 3 | [Applegate et al. 2025 abstract](https://arxiv.org/abs/2501.07018) |
| Power-system dispatch LP, 757k rows / 1.42M cols | cuOpt GPU PDLP **464 s** vs HiGHS CPU IPM 1400 s (**3.0×**); objectives agree within ~1e−3 relative; CPU wins below ~4-hour horizon | [PSR/NVIDIA benchmark 2026](https://www.psr-inc.com/en/analytics-report/post/gpu-based-algorithms-for-large-scale-optimization/) |
| Unit commitment hybrid Xpress-presolve + cuOpt PDLP (B200) | up to **20×** vs CPU-only barrier pipeline (>6 h → <20 min) | [Artelys/FICO/NVIDIA 2025](https://www.artelys.com/news/artelys-nvidia-fico-collaboration-power-system-optimisation/) |
| One MCF instance, GPU vs best CPU solver @1e−4 | up to 5000× speedup reported (outlier class); 10–300× typical on MCF | [NVIDIA cuOpt blog 2024](https://developer.nvidia.com/blog/accelerate-large-linear-programming-problems-with-nvidia-cuopt/) |

Iteration counts to expect: tens of thousands of KKT passes for hard Mittelmann instances at 1e−8 within the 100k-pass budget used in the PDLP ablations; PageRank-scale instances converge in thousands of iterations. Report progress in **KKT passes** (= one K + one K′ product), the standard hardware-noise-free metric ([Applegate et al. 2022, §4.1](https://arxiv.org/abs/2106.04756)).

---

## 8. Reference implementations (read, don't copy — license MIT for cuPDLP-C)

- **cuPDLP-C**: <https://github.com/COPT-Public/cuPDLP-C> — C + CUDA, MIT. Files to study: `cupdlp_step.c` (adaptive step + average), `cupdlp_restart.c` (restart logic), `cupdlp_solver.c` (main loop, check intervals), `cupdlp_scaling.c` (Ruiz + PC).
- **cuPDLP.jl**: <https://github.com/jinwen-yang/cuPDLP.jl> — Julia/CUDA reference prototype.
- **cuPDLPx**: <https://github.com/MIT-Lu-Lab/cuPDLPx> — successor using restarted **Halpern** PDHG + PID-controlled primal weight + constant step-size rule; **2.5–5× faster than cuPDLP on MIPLIB, 3–6.8× on Mittelmann** ([Lu, Peng, Yang 2025](https://arxiv.org/abs/2507.14051)). Worth a v2 upgrade path after v1 matches cuPDLP-C behavior.
- **Google OR-Tools PDLP** (CPU C++, production): `ortools/pdlp/primal_dual_hybrid_gradient.cc` in <https://github.com/google/or-tools> — best-commented reference for termination/infeasibility bookkeeping; math doc: <https://developers.google.com/optimization/lp/pdlp_math>.

## Sources

1. Chambolle & Pock, *A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging*, J. Math. Imaging Vision 40 (2011). https://doi.org/10.1007/s10851-010-0251-1
2. Applegate, Díaz, Hinder, Lu, Lubin, O'Donoghue, Schudy, *Practical Large-Scale Linear Programming using Primal-Dual Hybrid Gradient*, NeurIPS 2021 / arXiv:2106.04756. https://arxiv.org/abs/2106.04756
3. Applegate, Hinder, Lu, Lubin, *Faster First-Order Primal-Dual Methods for Linear Programming Using Restarts and Sharpness*, Math. Program. 201 (2023). https://arxiv.org/abs/2105.12715
4. Applegate, Díaz, Lu, Lubin, *Infeasibility Detection with Primal-Dual Hybrid Gradient for Large-Scale Linear Programming*, SIAM J. Optim. / arXiv:2102.04592. https://arxiv.org/abs/2102.04592
5. Lu & Yang, *cuPDLP.jl: A GPU Implementation of Restarted Primal-Dual Hybrid Gradient for Linear Programming in Julia*, arXiv:2311.12180. https://arxiv.org/abs/2311.12180
6. Lu, Yang, Hu, Huangfu, Liu, Liu, Ye, Zhang, Ge, *cuPDLP-C: A Strengthened Implementation of cuPDLP for Linear Programming by C language*, arXiv:2312.14832. https://arxiv.org/abs/2312.14832 · code: https://github.com/COPT-Public/cuPDLP-C
7. Lu & Yang, *Restarted Halpern PDHG for Linear Programming*, arXiv:2407.16144. https://arxiv.org/abs/2407.16144
8. Lu, Peng, Yang, *cuPDLPx: A Further Enhanced GPU-Based First-Order Solver for Linear Programming*, arXiv:2507.14051. https://arxiv.org/abs/2507.14051 · code: https://github.com/MIT-Lu-Lab/cuPDLPx
9. Applegate et al., *PDLP: A Practical First-Order Method for Large-Scale Linear Programming* (journal version), arXiv:2501.07018. https://arxiv.org/abs/2501.07018
10. Google Research blog, *Scaling up linear programming with PDLP* (Beale–Orchard-Hays Prize 2024). https://research.google/blog/scaling-up-linear-programming-with-pdlp/
11. OR-Tools PDLP math background. https://developers.google.com/optimization/lp/pdlp_math
12. NVIDIA cuOpt blogs/docs: PDLP acceleration (2024) https://developer.nvidia.com/blog/accelerate-large-linear-programming-problems-with-nvidia-cuopt/ · barrier-vs-PDLP accuracy guidance (2025) https://developer.nvidia.com/blog/solve-linear-programs-using-the-gpu-accelerated-barrier-method-in-nvidia-cuopt/ · LP/QP features https://docs.nvidia.com/cuopt/user-guide/latest/lp-qp-features.html
13. PSR Inc., *GPU-based algorithms for large-scale optimization* (Brazilian dispatch LP benchmark, 2026). https://www.psr-inc.com/en/analytics-report/post/gpu-based-algorithms-for-large-scale-optimization/
14. Artelys, *How Artelys, FICO and NVIDIA cuOpt Join Efforts to Scale Power System Optimisation* (2025). https://www.artelys.com/news/artelys-nvidia-fico-collaboration-power-system-optimisation/
