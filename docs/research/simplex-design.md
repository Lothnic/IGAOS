# Design Notes: From-Scratch CPU Revised-Simplex Solver

**Wayfinder ticket:** Lothnic/IGAOS#5 (parent map: #1) · **Date:** August 2026
**Scope:** LP core for SIH26119 (indigenous GPU-accelerated optimization solver). Audience: a student team with linear-algebra experience and **zero prior simplex implementation experience**. This document pins implementable decisions; it deliberately favors "what to build, in what order, with which tolerances" over survey breadth.

All claims are cited inline; full reference list with URLs in §10.

---

## 0. Pinned design (TL;DR)

| # | Decision | Pinned choice |
|---|---|---|
| D1 | Algorithm | **Bounded-variable dual revised simplex**, primal simplex kept only as cleanup/fallback |
| D2 | Basis representation | **Sparse LU (Markowitz + threshold pivoting) + product-form (PF) eta-file updates**, refactorize every 50–100 pivots |
| D3 | Pricing (dual) | **Devex reference-framework pricing**; exact dual steepest edge (Forrest–Goldfarb) as phase-3 upgrade |
| D4 | Ratio test | **Harris two-pass ratio test with δ=10⁻⁷ relaxation**, then bound-flipping ratio test (BFRT) on top |
| D5 | Degeneracy/cycling | **Cost perturbation à la HiGHS (base ≈ 5e-7·max\|c\|) + cost shifting at cleanup**; random tie-breaking in pricing |
| D6 | LU factorization | **Own sparse LU**: Markowitz count-minimization + threshold pivoting (τ=0.1), Gilbert–Peierls triangular solves; **dense LU fallback below m≈2000** |
| D7 | Scaling | **Iterative geometric-mean row/column scaling in powers of two before anything else runs** |

Rationale per item below.

---

## 1. Solver architecture

Implement the bounded-variable revised simplex in computational form

> min cᵀx s.t. Ax = 0, l ≤ x ≤ u,

with slack/logical columns folded into `A` so every variable carries bounds — this is the form used by HiGHS/hsol and makes the dual simplex uniform ([HH18], §2.1). The dual simplex is preferred as the main algorithm because dual steepest-edge-style pricing and the BFRT (both 1990s inventions) made it dramatically faster than the primal variant on practical instances, especially warm starts inside branch-and-bound later in the project ([HH18], §1).

Per iteration ([HH18], §2.2):

1. **chuzr** — choose leaving variable p among primal-infeasible basics (weighted by pricing rule).
2. **btran + spmv** — compute pivotal row αₚᵀ = eₚᵀB⁻¹A.
3. **chuzc** — Harris two-pass + BFRT choose entering variable q (§4 below).
4. **ftran** — compute B⁻¹a_q; update primal values and reduced costs.
5. **update-factor** — update the B⁻¹ representation; periodically refactorize.

Measured component costs from hsol on a real test set ([HH18], Table 1): ftran-dse 26.4%, spmv 18.4%, invert 13.3%, ftran 10.8%, btran 8.7%. Lesson: **the LU solves dominate runtime**, which is why D2/D6 get the most engineering care; the ratio-test bookkeeping is comparatively cheap (chuzc1+2 ≈ 8.8%).

---

## 2. Representing B⁻¹ (ticket question #1)

Options considered:

| Scheme | Idea | Sparsity of solves | Implement difficulty | Stability | Verdict |
|---|---|---|---|---|---|
| Explicit inverse / tableau | Store B⁻¹ or full tableau | Dense always | Trivial | Poor at scale | Reject beyond toy sizes; O(m²) memory |
| Product form (PF) eta file [Dantzig–Orchard-Hays 1954] | Keep initial LU₀; append one rank-1 eta E = I+(â_q−e_p)e_pᵀ per pivot; B_k⁻¹ = E_k⁻¹…E₁⁻¹B₀⁻¹ | Degrades linearly with eta count | Easy — reuses ftran result â_q | OK if refactorized regularly ("very rare for a good PF-based solver to fail", [HH14] §3.1) | **Pinned (D2)** |
| Forrest–Tomlin (FT) [1972] | Modify U directly (delete spiked row/col via row-eta R), push spike down with elimination; augment L | Best fill-in control; hyper-sparsity exploitable | Hard — dynamic data structures for deletion/insertion in U [HH14] | Good | Phase-4 stretch goal |
| Suhl–Suhl update [1993] | Rank-1 update of both L and U | Comparable/better than FT (Koberstein found FT "generally inferior … since it produces more fill-in" [Kob05] §8.2.4) | Medium-hard | Good | Optional alternative if FT attempted |
| Middle product form (MPF) [HH14] | PF variant storing T = I+u vᵀ etas | Approaches FT performance | Easy-medium | Same caveats as PF | Upgrade candidate after plain PF works |

Justification for pinning PF + periodic refactorization:

- Huangfu & Hall state plainly that MPF/PF-style updates are "**very much easier to implement than the FT update** … an attractive update procedure when developing a simple, relatively efficient implementation of the revised simplex method," with solve-time density approaching FT on many models ([HH14], §5 conclusions + Table 1 densities).
- The cost of PF degradation is capped by refactorizing every ~50–100 pivots; reinversion is only 13.3% of hsol's time even at that frequency ([HH18], Table 1). Refactorization is *also* your numerical reset button and your singularity-recovery point.
- The team's first solver must be correct before fast; PF has no dynamic data structures to debug, and the same ftran result needed anyway supplies the eta vector.
- Migration path is real: HiGHS itself uses Markowitz-LU + FT updates ([COINKB]); when (and only when) profiling shows eta-file solves dominating, swap in FT/Suhl behind the same FTRAN/BTRAN interface.

**Interface contract to code against from day 1:** `ftran(a) -> B⁻¹a`, `btran(e_p) -> e_pᵀB⁻¹`, `update(p, q)`, `invert()`. All four schemes above fit behind these four calls.

---

## 3. Pricing rules (ticket question #2)

For the dual simplex, "pricing" = choosing the leaving row p by weighted infeasibility Δx_i/w_i:

- **Dantzig:** w_i ≡ 1 (largest raw infeasibility). Cheapest per iteration; well known to need many more iterations on realistic problems because steps are short and degenerate pivots frequent ([Lingo17] discussion; standard result).
- **Devex** [Harris 1973]: approximate steepest-edge weights maintained cheaply from a reference framework, reset when weights drift >~2× reality. Cost per iteration barely above Dantzig; iteration counts close to exact steepest edge on most problems.
- **Exact dual steepest edge (DSE)** [Forrest–Goldfarb 1992]: w_i = ‖e_iᵀB⁻¹‖² updated each iteration via τ = B⁻¹ê_p (an extra ftran — that is hsol's single largest component, 26.4% of time, [HH18] Table 1). Fewest iterations, highest per-iteration cost; wins on large/hard instances.

What production codes do:

- **HiGHS default** (`simplex_dual_edge_weight_strategy=-1` "choose") resolves to exact DSE for the dual simplex with Forrest–Goldfarb weight updates, falling back to Devex when accumulated weight error exceeds a threshold (`dual_steepest_edge_weight_log_error_threshold` default 10); the primal simplex defaults to Devex ([AMPL26 options]; [COINKB]; verified in `HEkkDual.cpp`: unit weights for logical basis, `computeDualSteepestEdgeWeights` otherwise, `switchToDevex()` on error growth).
- CPLEX likewise auto-selects, offering steepest-edge variants and Devex as alternates ([IBM docs]).

**Pinned (D3): Devex first.** It delivers most of the iteration-count benefit at a fraction of the implementation risk, needs no extra ftran, and its reference-framework reset logic is forgiving of beginner bugs. Port the Forrest–Goldfarb DSE update formulas (given verbatim in [HH18] §2.2.3) once Devex passes Netlib smoke tests — the code path shares the same weight array. Do not ship Dantzig except as a debugging baseline.

---

## 4. Ratio tests (ticket question #3)

### 4.1 Harris two-pass ratio test [Harris 1973]

Why: a textbook ratio test picks the exact minimum θ, which routinely selects tiny-magnitude pivot elements → unstable bases → garbage iterates. Harris' fix: trade a *tiny, controlled* amount of dual feasibility for a much better pivot.

Mechanics (dual version):

1. Compute θ_max = max allowed dual step keeping all reduced costs within their tolerance band.
2. **Pass 1:** compute relaxed bound θ̃ = θ_max·(1+δ) (HiGHS-class implementations use δ≈1e-7 … 1e-9; Koberstein details parameterization [Kob05] ch. 5–6). Collect all candidates j with |α_pj| > 0 whose ratios fall within θ̃.
3. **Pass 2:** among those candidates, pick q maximizing |α_pq| (best pivot magnitude) rather than minimizing ratio. The resulting small dual infeasibilities are then repaired by the shifting/perturbation machinery of §5.

This is exactly how HiGHS structures chuzc1/chuzc2 ([HH18] §2.2.2), and MOSEK's rewritten simplex lists "bound-flipping ratio test with Harris tolerance" as a core ingredient [Fri25].

**Pinned (D4a): Harris two-pass with δ = 1e-7, minimum pivot magnitude floor ~1e-7 (cf. HiGHS `dual_simplex_pivot_growth_tolerance` 1e-9, `factor_pivot_threshold` 0.1 for LU).**

### 4.2 Bound flipping ratio test (BFRT) [Fourer 1994; Kob05 §3.2]

Briefly: when a nonbasic boxed variable would become blocking *before* the chosen breakpoint, flip it to its opposite bound instead of stopping there; its contribution to the pivotal row RHS is folded into an aggregate column a_F, requiring one extra ftran (ftran-bfrt) per iteration where flips occur. Dual step length grows substantially → fewer iterations, at negligible risk. Both Huangfu–Hall and Koberstein credit BFRT + DSE as the reasons the dual simplex dominates today ([HH18] §1, §2.2.2).

**Pinned (D4b): implement BFRT after plain Harris works.** It touches only chuzc and one extra ftran call — safe increment.

---

## 5. Degeneracy and cycling (ticket question #4)

Practical anti-cycling kit, in the order implemented by modern solvers:

1. **Cost perturbation (primary workhorse).** Perturb costs slightly so ties/degenerate steps vanish statistically. HiGHS' concrete recipe (verified in `HEkk.cpp ::initialiseCost`, perturb=true path):
   - base = multiplier × 5e-7 × max|c_j|; if max|c_j| > 100 use √√(max|c_j|) instead; if <1% of variables are boxed cap base at 1.0;
   - for each nonbasic column: Δc_j = ±(1 + rand) × (|c_j|+1) × base, sign toward its active bound (+ at lower, − at upper, sign(c_j) if boxed); skip free and fixed variables; row costs additionally jittered at ~1e-12 scale;
   - option `dual_simplex_cost_perturbation_multiplier` default 1 (=on).
   Huangfu–Hall's paper documents this design lineage from hsol through Xpress [HH18].
2. **Cleanup without silent lies.** At termination, remove the perturbation, recompute duals, and repair any residual dual infeasibilities by **cost shifting** (add a shift to c_q of the offending nonbasic column so its reduced cost becomes feasible; shifts are reported, and HiGHS falls back to a primal-simplex "cleanup" solve if shifting alone is insufficient — see `HEkkDual::cleanup()` and the `kSolvePhaseOptimalCleanup` path in source). Never declare optimality under perturbed costs.
3. **Random tie-breaking** in pricing reduces stalling; Koberstein reports simple randomized tie resolution "very useful to reduce stalling and diminish the risk of cycling" [Kob05] §6.1.
4. **Theory backstop:** Bland's rule guarantees finite termination but is far too slow to run by default; keep it as a debug flag only. DSE/Devex + Harris already greatly reduce degenerate iterations [Kob05] §6.1.

**Pinned (D5): deterministic-seed random cost perturbation per HiGHS recipe + shift-at-cleanup + random tie-breaks.** Bound perturbation for the primal side can wait until the primal cleanup path exists.

---

## 6. Own sparse LU factorization, no external libs (ticket question #5)

The hardest component. Required pieces, all classical:

1. **Markowitz pivoting** [Markowitz 1957]: at each elimination step choose the pivot entry minimizing (row_count−1)×(col_count−1) among candidates. Exact Markowitz is expensive; the standard approximation processes columns in ascending nonzero count and considers only a few cheapest rows/columns — this is what Suhl & Suhl describe for LP bases [SS90].
2. **Threshold pivoting** for stability: accept a Markowitz candidate only if |pivot| ≥ τ × max|entries in its column|; HiGHS ships `factor_pivot_threshold` = **0.1** (range [0.0008, 0.5]) [GAMS-HIGHS]. Lower τ ⇒ sparser but less stable.
3. **Triangular solves:** Gilbert–Peierls DFS-based solve costs O(nnz(L)+nnz(U)), not O(m) [GP88]; add hyper-sparse handling (stop when the solve result stays sparse; exploit it in btran/spmv) — this is what makes large sparse LPs tractable [HH18] §3.4, Hall–McKinnon hyper-sparsity.
4. **Singular-basis recovery:** keep the last good basis; on refactorization failure repair via the standard basis-repair pass (Koberstein ch. 5; HiGHS `getNonsingularInverse`).

Realistic effort estimate for students (first sparse-LU ever, C++20, no external libs):

| Component | Estimate |
|---|---|
| Dense LU + partial pivoting, triangular solves | 1 week (do it first, it is also the fallback) |
| Sparse data structures (CSC/CSC-with-row-swap) | 1 week |
| Markowitz + threshold factorization | 3–5 weeks including tests |
| GP-style sparse FTRAN/BTRAN + hyper-sparse mode | 2–3 weeks |
| Update integration (PF eta file) + refactor cadence | 1–2 weeks |
| **Total** | **≈ 8–12 person-weeks** — schedule as the critical path |

Calibration: MOSEK's professional rewrite of a simplex engine took years even with experienced staff (>100 kLOC legacy replaced starting 2020) [Fri25]; a student scope of Netlib-class robustness (not Mittelmann-scale dominance) at 8–12 weeks is consistent with what prior student projects achieved (§7).

**Dense fallback threshold (D6):** below **m ≈ 2000** rows, just use dense LU (or even dense revised simplex with explicit B⁻¹): dense factorization is O(m³) ≈ 8e9 flops at m=2000 (~seconds, amortized over hundreds of iterations), memory m² = 4M doubles = 32 MB, and correctness is easy to guarantee. Above m ≈ 2000 the m² storage and cubic refactor cost start to dominate iteration time, and sparsity pays for its own complexity. Concretely: build the dense path first (it unlocks afiro/sc205/etamacro-class Netlib testing immediately), gate the sparse path behind `m > 2000 || nnz(A)/m < 50`.

---

## 7. Feasibility signals from prior student work (ticket question #7)

Reference repo: https://github.com/alaintis/gpu-accelerated-simplex-solver (C++/CUDA, validated on Netlib against cuOpt/HiGHS as golden reference). Verified from README/repo structure:

- **What they built:** GPU-native *revised* simplex with full GPU residency, **Sherman–Morrison basis updates** (PF-equivalent rank-1 update), CUDA Graphs to strip kernel-launch overhead, **parallel Harris ratio test**, sparse MPS pipeline, CPU/GPU swappable backends, plus a correctness runner and NSight-based profiling workflow [alaintis]. A related published effort similarly reports effective GPU revised simplex with new memory management and cycle avoidance [GPU18].
- **What worked / transfers to our CPU project:**
  - Scope discipline: they solved *regularly shaped* LPs (Ax ≤ b, x ≥ 0) first, added MPS + Netlib harness, then tuned. Mirror that ladder.
  - Sherman–Morrison/PF updates were enough to reach optimal solutions on Netlib instances — corroborates D2 for a student team.
  - Harris ratio test was implementable in parallel/simple form — corroborates D4a as early milestone.
  - Golden-reference testing vs HiGHS/cuOpt per instance (objective delta = 0 checks) — adopt identical methodology: `highs` binary as oracle, compare objective values and basis status on every Netlib run.
- **Cautionary signal:** sequential simplex on GPU does not beat good CPU implementations generally; their win condition was launch-overhead elimination at small scale [HH18] §1 notes parallel/tableau simplex loses to good sequential revised simplex on sparse LPs. Our project's GPU story lives in PDLP (separate ticket), not here.

Conclusion: the pinned scope (§0) is demonstrably achievable by a team of this profile within a semester-length effort.

---

## 8. Scaling (ticket question #6)

Run scaling **before anything else**; unscaled industrial matrices wreck both LU stability and ratio-test tolerances.

- **Pinned (D7): iterative geometric-mean scaling.** Repeat a few sweeps of r_i ← r_i/geom-mean(row i), c_j ← c_j/geom-mean(col j) using **powers of two only** (exact, no roundoff). Tomlin's classic comparative study concluded geometric-mean scaling, optionally followed by equilibration or Curtis–Reid, is the best combination [Tom75]; HiGHS' default `simplex_scale_strategy` = 2 implements this family (equilibration/choose) with power-of-two factors bounded by `allowed_matrix_scale_factor`=20 exponents [GAMS-HIGHS]; Xpress historically defaulted to power-of-two equilibration with Curtis–Reid as numerics-rescue option [BH21].
- Stop criterion: matrix coefficient range (max/min |a_ij|) stops improving or < ~10 sweeps.
- Scale costs too (power-of-two objective factor; HiGHS `cost_scale_factor`) so dual tolerances behave uniformly.
- Solve scaled; unscale solution; optionally verify residuals unscaled (HiGHS `simplex_unscaled_solution_strategy`).
- Curtis–Reid least-squares scaling [CR72] is the documented rescue for pathological instances — optional phase-3 item, not day-one scope.

---

## 9. Phased build order (with exit criteria)

| Phase | Build | Exit criterion |
|---|---|---|
| 0 | Dense revised simplex, Dantzig, big-M-free bounded form, MPS reader | Optimal on afiro + 20 small Netlib, matches HiGHS objective |
| 1 | Scaling (D7) + Harris two-pass (D4a) + Devex (D3) | Iteration counts within ~2× of HiGHS on sc50a-class set |
| 2 | PF updates + refactor-every-100 (D2), dense→sparse gate at m<2000 | Solves 50a–105a range; refactor count logged and sane |
| 3 | Cost perturbation + shifting cleanup (D5), BFRT (D4b) | No cycling on degenerate set (e.g., degme class); BFRT cuts dual iterations measurably |
| 4 | Sparse LU (Markowitz+threshold, GP solves, hyper-sparse) (D6) | Netlib medium set (m up to ~15k) completes |
| 5 (stretch) | Exact DSE weights (FG update), then FT/Suhl update behind same interface | Iteration/time parity trend vs HiGHS |

Effort totals: phases 0–3 ≈ 6–8 person-weeks; phase 4 is the 8–12 person-week critical path from §6; phase 5 optional.

---

## 10. References (with URLs)

Primary algorithmic sources:

1. **Huangfu & Hall, "Parallelizing the dual revised simplex method," Math. Prog. Comp. 10(1):119–142, 2018.** arXiv: https://arxiv.org/abs/1503.01889 (HTML: https://arxiv.org/html/1503.01889v1) — component timing table, DSE weight updates, Harris+BFRT structure, hyper-sparsity. **Start here.**
2. **Huangfu & Hall, "Novel update techniques for the revised simplex method," COAP 60:587–608, 2014** (ERGO-13-001). PDF: https://optimization-online.org/wp-content/uploads/2013/02/3774.pdf ; https://link.springer.com/article/10.1007/s10589-014-9689-1 — PF vs FT vs MPF difficulty/perf; the justification for D2.
3. **Koberstein, "The dual simplex method, techniques for a fast and stable implementation," PhD thesis, Paderborn, 2005.** https://digital.ub.uni-paderborn.de/hs/download/pdf/3885?originalFilename=true (mirror: https://d-nb.info/978580478/34) — ch. 3 dual simplex + BFRT; ch. 5 LU/update; ch. 6 anti-degeneracy/stability. **Read second, selectively.**
4. **Maros, "Computational Techniques of the Simplex Method," Kluwer, 2003** (ISORMS vol. 61) — background reading for revised simplex mechanics, PFI, LU updating, degeneracy chapters; use as the theory companion while coding phases 0–2. Entry: https://link.springer.com/book/10.1007/978-1-4615-0257-9
5. **Fourer, "Notes on the dual simplex method," 1994 (unpublished).** http://users.iems.northwestern.edu/~4er/WRITINGS/dual.pdf — BFRT origin.
6. **Harris, "Pivot selection methods of the Devex LP code," Math. Programming 5:1–28, 1973.** https://link.springer.com/article/10.1007/BF01580108 — Devex + Harris two-pass ratio test.
7. **Forrest & Goldfarb, "Steepest-edge simplex algorithms for linear programming," Math. Programming 57:341–374, 1992.** https://link.springer.com/article/10.1007/BF01581089 — exact DSE updates (phase 5).

Sparse LU:

8. **Markowitz, "The elimination form of the inverse and its applications to linear programming," Management Science 3:255–269, 1957** — pivot criterion.
9. **Suhl & Suhl, "Computing sparse LU factorizations for large-scale LP bases," INFORMS J. Computing 2(4):325–335, 1990** — practical Markowitz+threshold recipe. Companion update paper: Annals of OR 43:33–47, 1993.
10. **Gilbert & Peierls, "Sparse partial pivoting in time proportional to arithmetic operations," SIAM J. Sci. Stat. Comput. 9:862–874, 1988** — O(nnz) triangular solves.
11. Elble & Sahinidis, "A review of the LU update in the simplex algorithm," IJMOR 4(4):366–399, 2012 — survey backing §2 table.

Scaling:

12. Tomlin, "On scaling linear programming problems," Mathematical Programming Study 4:146–166, 1975 — comparative study favoring geometric mean (+equilibration/Curtis-Reid). Context summary: https://www.researchgate.net/publication/304744500_The_impact_of_scaling_on_simplex_type_algorithms
13. Curtis & Reid, "On the automatic scaling of matrices for Gaussian elimination," IMA J. Appl. Math. 10:118–124, 1972.
14. Berthold & Hendel, "Learning to scale mixed-integer programs," AAAI 2021. https://ojs.aaai.org/index.php/AAAI/article/view/16482 — modern practice context (Xpress standard vs CR scaling).

Solver ground truth (defaults & behavior):

15. HiGHS options (defaults incl. `simplex_dual_edge_weight_strategy=-1`, `simplex_scale_strategy=2`, `factor_pivot_threshold=0.1`, `dual_simplex_cost_perturbation_multiplier=1`): https://www.gams.com/latest/docs/S_HIGHS.html ; https://dev.ampl.com/solvers/highs/options.html
16. HiGHS source (cost perturbation recipe, DSE↔Devex switching, cleanup/shifting): https://github.com/ERGO-Code/HiGHS/blob/master/highs/simplex/HEkk.cpp , https://github.com/ERGO-Code/HiGHS/blob/master/highs/simplex/HEkkDual.cpp
17. COIN-OR Knowledge Base on HiGHS internals (dual DSE w/ FG weights, primal Devex, Markowitz LU + FT updates): https://monistowl.github.io/coin-or-kb/libraries/highs/
18. Friberg (MOSEK), "On reimplementing the simplex optimizers in MOSEK," SIAM OP26 slides, 2026. https://docs.mosek.com/slides/2026/siopt/talk-simplex.pdf — ingredient checklist + realistic timeline calibration.

Prior student work:

19. alaintis, gpu-accelerated-simplex-solver (README + structure). https://github.com/alaintis/gpu-accelerated-simplex-solver

CPLEX comparison points (for benchmark framing only — never wrapped): IBM CPLEX dual-pricing parameter docs https://www.ibm.com/docs/en/icos/22.1.2?topic=parameters-dual-simplex-pricing-algorithm
