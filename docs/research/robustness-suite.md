# Robustness Test Suite — degeneracy / ill-conditioning / edge-case LPs

**Ticket:** Lothnic/IGAOS#21 · **Date:** 2026-08-25 · **Status:** DRAFT (awaiting human lock)
**Scope:** a 15-instance Netlib(+Chinneck) suite that exercises the adopted robustness strategy — geometric-mean scaling → Harris two-pass ratio test (δ=1e−7) → bound-flipping ratio test → cost shifting + perturbation vs degeneracy → iterative refinement on residual failure; acceptance ladder `solved@1e−4` / `tight@1e−6`.
**Relationship to benchmark protocol:** this suite *reuses* the conventions of [benchmark-protocol.md](benchmark-protocol.md) (tiers, TLs, canonical table, honest-losses rules) but applies **robustness-specific PASS criteria** (§3). Instances overlapping the main ladder are deliberate: same file, different gate.

Everything below verified against primary sources August 2026 (sources §6).

---

## 1. Instance table

Sizes and optimal values from the Netlib index (`netlib.org/lp/data/readme`; rows/cols/nnz exclude slacks/RHS per Netlib convention). Infeasible-family sizes from Chinneck's index (`netlib.org/lp/infeas/readme`). Degeneracy percentages = MINOS 5.3 "steps degenerate", as tabulated in the Netlib index for the Tomlin-supplied problems. All files confirmed present in the SkyLiu0/netlib mirror (`github.com/SkyLiu0/netlib`, tree probed via GitHub API 2026-08-25): feasible under `feasible/<name>.mps`, infeasible under `infeasible/<name>.mps`.

| # | Instance | Class | m | n | nnz | Stress mechanism | Optimal value (source) | Availability |
|---|---|---|---:|---:|---:|---|---|---|
| 1 | sc50a | B dual-degen | 51 | 48 | 131 | Set-covering LP; many alternative optima expected; tiny ⇒ cheap CI canary for alt-optima detection | −6.4575077059E+01 (netlib) | ✓ mirror |
| 2 | sc50b | B dual-degen | 51 | 48 | 119 | Twin of sc50a, different data — isolates structure-driven vs data-driven ties | −7.0000000000E+01 (netlib) | ✓ mirror |
| 3 | kb2 | A primal-degen (+empty RHS) | 44 | 41 | 291 | ~35% degenerate steps (MINOS); empty RHS section — reader + phase-1 edge | −1.7499001299E+03 (netlib) | ✓ mirror |
| 4 | klein1 | F infeasible | 55 | 54 | 696 | Tiny infeasible (Chinneck/Klotz); fast smoke test for infeasibility proof path | infeasible (Chinneck idx) | ✓ mirror |
| 5 | capri | E free/FX bounds | 272 | 353 | 1786 | 14 FR columns + FX bounds (probe-verified); Stanford SOL model | 2.6900129138E+03 (netlib) | ✓ mirror |
| 6 | israel | C ill-conditioned | 175 | 142 | 2358 | Wide aspect ratio, huge objective range (−8.97e5), dense columns; classic scaling stress | −8.9664482186E+05 (netlib) | ✓ mirror |
| 7 | modszk1 | A primal-degen | 688 | 1620 | 4158 | Maros: "very degenerate"; dual simplex may need up to 10× fewer iters than primal (his quote, netlib idx); 2 FR cols | 3.2061972906E+02 (netlib) | ✓ mirror |
| 8 | boeing1 | D RANGES-heavy | 351 | 384 | 3865 | 45 RANGES entries (probe-verified); BR-flagged; originally had INTORG markers (removed) | −3.3521356751E+02 (netlib) | ✓ mirror |
| 9 | degen2 | A primal-degen | 445 | 534 | 4449 | ~57% degenerate steps (MINOS) — the named stress test #1; cost negated by netlib (minimization) | −1.4351780000E+03 (netlib) | ✓ mirror |
| 10 | pilot4 | E free/PL bounds | 411 | 1000 | 5145 | 88 FR columns + FX + PL types (probe-verified) — widest bound-type coverage in lp/data | −2.5811392641E+03 (netlib) | ✓ mirror |
| 11 | refinery | F infeasible | 324 | 464 | 1694 | Doctored petrochemical plant (volume/quality), volatile IIS; MRPL-domain demo tie-in ([refinery-case-studies.md](refinery-case-studies.md)) | infeasible (Chinneck idx) | ✓ mirror |
| 12 | pilot87 | C ill-conditioned | 2031 | 4883 | 73804 | Lustig (quoting Stone, netlib idx): "harder than PILOT because of the bad scaling in the numerics" | 3.0171072827E+02 (netlib) | ✓ mirror |
| 13 | degen3 | A primal-degen | 1504 | 1818 | 26230 | ~53% degenerate steps at scale — stress test #2; largest dense-ish degenerate core | −9.8729400000E+02 (netlib) | ✓ mirror |
| 14 | nesm | D RANGES-heavy | 663 | 2923 | 13988 | 44 RANGES entries (probe-verified); explicit zeros present (emps caveat); naval logistics model | 1.4076073035E+07 (netlib) | ✓ mirror |
| 15 | klein3 | F infeasible | 995 | 88 | 12107 | Tall-thin infeasible (Klotz); medium-scale phase-1 / conflict-detection stress | infeasible (Chinneck idx) | ✓ mirror |

Class legend: **A** heavy primal degeneracy · **B** dual degeneracy / multiple optima · **C** ill-conditioned matrices · **D** RANGES stress · **E** free/bounded-column edge cases · **F** infeasibility detection.

### Class rationale & known folklore (per class)

- **A — primal degeneracy.** Primary source is the Netlib index's own table (MINOS 5.3 on Tomlin's problems): DEGEN2 56.74%, DEGEN3 52.51%, KB2 35.37%; MODSZK1 has no percentage but Maros' submitted note calls it "very degenerate" with dual-simplex iteration ratios up to 10× better — making it our cleanest Harris/BFRT-vs-perturbation discriminator. Expected behavior: no cycling, iteration counts within the harness guard (§3); cost shifting + perturbation should visibly cut degenerate steps on degen2/degen3.
- **B — dual degeneracy / multiple optima.** No Netlib instance carries an official "multiple optima" flag; following Gamrath–Berthold–Salvagnin (2020), dual degeneracy = multiple optimal bases, detected operationally as nonbasic variables with zero reduced cost at optimum. The sc* set-covering family (sc50a/b here) structurally invites ties; we *measure* rather than assume (criterion B-PASS below). sc50a/sc50b also give the cheapest possible reproduction of the classic "same basis, different data" trap.
- **C — ill-conditioning.** ISRAEL: wide (142 cols vs 175 rows) with objective coefficients spanning orders of magnitude; documented solver-spread values do not exist for it, so it gates at 1e−6 rel. PILOT87: bad scaling explicitly blamed by its contributor (primary-source quote above). DFL001 (Bixby's κ∞(scaled basis)=213737 note, tolerance-dependent optima) stays on the main ladder stretch tier; add it here as optional extension if the suite proves cheap to run.
- **D — RANGES stress.** The Netlib index flags exactly five feasible problems with ranges (BR column): BOEING1, BOEING2, FORPLAN, SEBA, NESM. We take the two heaviest by entry count (45 and 44, probe-verified against the mirror MPS files). These catch range-conversion bugs (GE→E flips, range slack bookkeeping) that no other class exposes.
- **E — free/bounded columns.** CAPRI (14 FR + FX) and PILOT4 (88 FR + FX + **PL**, unique in lp/data) probe reader and basis-handling paths for free columns. **Negative-UP-bounds (UP<0 ⇒ lb=−inf) and MI bound types could NOT be verified in any Netlib feasible instance**: the index's own BOUND-TYPE TABLE contains no MI anywhere, and a grep probe of 11 candidate MPS files found zero negative UP entries. This edge case therefore gets a synthetic unit fixture in the reader tests, not a suite slot.
- **F — infeasibility detection.** Chinneck collection (sizes from his summary table): KLEIN1 tiny smoke, KLEIN3 tall/thin medium, REFINERY domain-relevant volatile-IIS case. Gate is *proven* detection, never a timeout or fake optimum.

---

## 2. Suite definition (ordered list, as consumed by the harness)

Order = tiny-first interleaving so CI fails cheap and early; tiers/TLs inherited from benchmark-protocol §6.2 (T0/T1 60 s, T2 300 s).

```yaml
# benchmarks/suites/robustness.yaml
name: robustness
acceptance_ladder: {solved: 1e-4, tight: 1e-6}   # CONTEXT.md "near-optimal" vocabulary
instances:
  - {id: sc50a,    class: B, tier: T0, mps: netlib/mps/sc50a.mps,      expect: optimal,   z_ref: -6.4575077059e+01}
  - {id: sc50b,    class: B, tier: T0, mps: netlib/mps/sc50b.mps,      expect: optimal,   z_ref: -7.0e+01}
  - {id: kb2,      class: A, tier: T0, mps: netlib/mps/kb2.mps,        expect: optimal,   z_ref: -1.7499001299e+03}
  - {id: klein1,   class: F, tier: T0, mps: netlib/infeas/klein1.mps,  expect: infeasible}
  - {id: capri,    class: E, tier: T1, mps: netlib/mps/capri.mps,      expect: optimal,   z_ref: 2.6900129138e+03, fr_cols: 14}
  - {id: israel,   class: C, tier: T1, mps: netlib/mps/israel.mps,     expect: optimal,   z_ref: -8.9664482186e+05}
  - {id: modszk1,  class: A, tier: T2, mps: netlib/mps/modszk1.mps,    expect: optimal,   z_ref: 3.2061972906e+02}
  - {id: boeing1,  class: D, tier: T1, mps: netlib/mps/boeing1.mps,    expect: optimal,   z_ref: -3.3521356751e+02, ranges: 45}
  - {id: degen2,   class: A, tier: T2, mps: netlib/mps/degen2.mps,     expect: optimal,   z_ref: -1.4351780000e+03}
  - {id: pilot4,   class: E, tier: T2, mps: netlib/mps/pilot4.mps,     expect: optimal,   z_ref: -2.5811392641e+03, fr_cols: 88}
  - {id: refinery, class: F, tier: T2, mps: netlib/infeas/refinery.mps,expect: infeasible}
  - {id: pilot87,  class: C, tier: T3, mps: netlib/mps/pilot87.mps,    expect: optimal,   z_ref: 3.0171072827e+02}
  - {id: degen3,   class: A, tier: T3, mps: netlib/mps/degen3.mps,     expect: optimal,   z_ref: -9.8729400000e+02}
  - {id: nesm,     class: D, tier: T2, mps: netlib/mps/nesm.mps,       expect: optimal,   z_ref: 1.4076073035e+07, ranges: 44}
  - {id: klein3,   class: F, tier: T2, mps: netlib/infeas/klein3.mps,  expect: infeasible}
```

Notes:
- `mps:` paths resolve under `benchmarks/` (gitignored; fetched via the benchmark-protocol §3.1 `emps` procedure for feasible files, and from `netlib.org/lp/infeas` for class F — both mirrored at SkyLiu0/netlib as fallback).
- degen2 value guard: use `-1.4351780000e+03` (the netlib index value); the `-…E+00` figure in benchmark-protocol §3.2 row is a typo already flagged there.
- Optional extension slots (not part of the locked 15): `dfl001` (class C), `cycle`/`tuff`/`woodw` (class A, 47%/46%/39% degenerate, empty-RHS variant via tuff), `klein2`, `boeing2`.

---

## 3. PASS criteria (harness-measurable, per class)

Global gates for every row: status matches `expect`; time within tier TL; run recorded in canonical table format (§4). Ladder: acceptance at `solved@1e−4`, record `tight@1e−6`.

| Class | PASS criterion (all must hold) |
|---|---|
| A primal-degen | Reaches @1e−4 within TL **and** \|obj − z_ref\| ≤ 1e−9·max(1,\|z_ref\|) **and** iteration guard: ≤ 5× HiGHS v1.15.1 iterations on same machine (guard, not publication metric) **and** no stall detector fired (no cycling: N consecutive degenerate pivots without objective/rhs change ⇒ abort+report) |
| B dual-degen | Same as A **plus** perturbation-stability check: re-solve with costs perturbed c′=c+δ, δ~U(±1e−7·‖c‖∞), seed=1; objective must match z_ref within 1e−9 rel (basis may differ — difference is reported, not penalized) **and** harness reports zero-rc count: #nonbasic with \|rc\| ≤ 1e−9 at optimum |
| C ill-conditioned | Reaches @1e−4 within TL **and** \|obj − z_ref\| ≤ 1e−6·max(1,\|z_ref\|) (relaxed: netlib itself documents cross-solver spread at the 8th digit for this family) **and** appended KKT residual columns (primal/dual, unscaled) reported **and** iterative-refinement fallback engaged if residuals > 1e−7 after solve — engagement logged, not gated |
| D RANGES | Reader assertion: #RANGES entries consumed == manifest `ranges` value **and** same objective gate as A (1e−9 rel) |
| E free/bounds | Reader assertions: #FR columns parsed == manifest `fr_cols`; PL accepted as ub=+inf; no spurious lb=0 applied to FR columns **and** same objective gate as A |
| F infeasible | Proven infeasible status (presolve conflict **or** phase-1 dual ray certificate), not timeout/cycle-out **and** wall time ≤ TL with iteration-cap guard **and** (diagnostic only, unpunished) IIS size if isolation implemented |

Failure taxonomy mirrors protocol §6.1 statuses (`error` for reader/assertion failures, `timeout-no-incumbent`, etc.) so robustness regressions are distinguishable from perf regressions.

## 4. Rendering in the canonical table

Ten fixed columns exactly as benchmark-protocol §7 (`instance,m,n,nnz,our_time_s,highs_time_s,our_objective,highs_objective,rel_gap,status`); class-F rows carry `status=infeasible` with blank objectives/gap. Appended diagnostic columns (right side only, per protocol rule):

```
...,class,zero_rc_nonbasic,perturbed_obj_rel_diff,ranges_parsed,fr_cols_parsed,kkt_primal_res,kkt_dual_res,ir_fallback_used,iis_rows
```

Suite aggregate line (per published table footer): #PASS/class, shifted geometric mean runtime (shift 10 s), plus HiGHS version, our git SHA, machine-sheet ID, date, thread count — identical footnoting discipline as the main ladder.

## 5. Integration notes

- Runs as `igaos bench --suite robustness` intent: loads `benchmarks/suites/robustness.yaml`, resolves MPS files under `benchmarks/netlib/`, executes ordered list cold-start, `random_seed=0`, single-thread primary track (protocol §2 flags), emits one canonical-table CSV per configuration under `benchmarks/logs/<date>-<machine-id>/`.
- HiGHS comparison column runs the identical pinned command line; no per-instance tuning anywhere (protocol §9).
- The perturbation re-solve (class B) and reader assertions (D/E) execute inside the harness runner, not inside the solver; the solver binary sees ordinary files.
- Suite is additive to the Netlib ladder: overlapping instances (kb2, modszk1, degen2, degen3, israel, pilot87) share downloaded artifacts; no second download pipeline.

## 6. Sources

- Netlib LP index (sizes, optima, Tomlin degeneracy table, Maros MODSZK1 note, Lustig PILOT87 quote, Bixby DFL001 κ note, BR flags, BOUND-TYPE TABLE): http://www.netlib.org/lp/data/readme
- Chinneck infeasible-LP collection (descriptions + summary table): http://www.netlib.org/lp/infeas/readme
- SkyLiu0/netlib mirror availability probe (tree API, 2026-08-25): https://github.com/SkyLiu0/netlib (114 feasible + 29 infeasible MPS)
- Gamrath, Berthold, Salvagnin: "An exploratory computational analysis of dual degeneracy in mixed-integer programming", EURO J. Comput. Optim. 8 (2020), doi:10.1007/s13675-020-00130-z (dual-degeneracy definition, zero-rc metric)
- Lustig, "An Analysis of an Available Set of Linear Programming Test Problems", SOL 87-11 / Comput. Opns. Res. 16 (1989) (referenced from netlib index; MINOS scaling/pricing behavior)
- Repo: docs/research/benchmark-protocol.md (tiers, TLs, canonical format, §9 fair-comparison), docs/SIH26119-RESEARCH-REPORT.md §2–3 (PS robustness emphasis, numerics kit), docs/research/refinery-case-studies.md (refinery instance dossier)
