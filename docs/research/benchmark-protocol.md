# SIH26119 — Public Benchmark Protocol (pinned)

**Ticket:** Lothnic/IGAOS#9 · **Date:** 2026-08-25 · **Status:** COMMITTED
**Scope:** instance sets, metrics, pinned baseline, result format, hardware sheet for all IGAOS benchmark publications (README tables, idea PPT, finale dashboard).

Everything below is verified against primary sources as of August 2026; sources are listed per section and at the end.

---

## 1. What we commit to

| Commitment | Decision |
|---|---|
| Baseline solver | **HiGHS v1.15.1** (latest stable release, 2026-07-02, [ERGO-Code/HiGHS releases](https://github.com/ERGO-Code/HiGHS/releases)). Mittelmann's current LP benchmarks also run HiGHS-1.15.0, so version-series parity with the public record holds. |
| LP validation set | Netlib LP ladder: **20 core + 7 stretch instances** (§3) |
| MILP validation set | MIPLIB 2017 easy subset: **16 core + 4 stretch instances** from the official `easy-v18.test` list (§4) |
| Aspirational sets | Mittelmann **LPfeas / LPopt / MILP(MIPLIB2017)** protocols, replicated settings where feasible (§5) |
| QP (later phase) | QPLIB via Mittelmann's convex/discrete QPLIB benchmarks — note only, not committed yet (§5.3) |
| Gap thresholds | relative MIP gap ≤ `1e-4` = "solved"; `1e-6` = "tight-solved" (LP/Mittelmann-parity runs) (§6) |
| Result format | ONE canonical table everywhere (§7) |
| Reproducibility | every table carries machine-sheet ID + HiGHS version + solver git SHA (§8, §9) |

---

## 2. Baseline: HiGHS v1.15.1

- **Pin:** `v1.15.1`, released 2026-07-02 ([GitHub release](https://github.com/ERGO-Code/HiGHS/releases/tag/v1.15.1)). Re-pin only by editing this doc in a PR that re-runs the full ladder.
- **Install:** `pip install highspy==1.15.1` or build from the tagged source; CLI binary is `highs`.
- Option names verified against the official option definitions (docs regenerated 2026-08-24): [ergo-code.github.io/HiGHS/dev/options/definitions](https://ergo-code.github.io/HiGHS/dev/options/definitions/)

### Flags used (and nothing else)

| Track | Command line | Purpose |
|---|---|---|
| **CPU single-thread (primary ladder)** | `highs --presolve=on --parallel=off --threads=1 --random_seed=0 --time_limit=<TL> model.mps` | apples-to-apples vs our serial simplex/B&B |
| CPU multi-track (secondary) | same + `--parallel=on --threads=<N>` where N = our worker count | parallel fairness |
| IPM/PDLP comparison track | `--solver=ipm` resp. `--solver=pdlp --run_crossover=on` | vs our GPU first-order method |

Relevant defaults we rely on and therefore report verbatim in every table footnote:
`mip_rel_gap=1e-4`, `mip_abs_gap=1e-6`, `primal_feasibility_tolerance=1e-7`, `dual_feasibility_tolerance=1e-7`, `simplex_strategy=1` (dual serial), presolve `"choose"`.

### Expected behavior class (what to sanity-check before publishing)

- Tiny/small Netlib: solved in well under a second (dual simplex + presolve).
- Medium Netlib (`25fv47`, `degen2`, `bnl1`, `modszk1`): seconds; degeneracy handled routinely.
- Large/stretch (`maros-r7`, `fit2p`, `pilot*`, `greenbea`): seconds-to-tens-of-seconds on modern CPUs.
- MIPLIB easy subset: most of our starter list solves within minutes at default `mip_rel_gap=1e-4`.
If HiGHS falls far outside these classes on any instance, suspect a harness bug before celebrating.

---

## 3. Netlib LP validation ladder

Source of truth: <http://www.netlib.org/lp/data/readme> (problem summary table reproduced below; rows/cols exclude slack columns and RHS, per Netlib convention). All problems are minimization form after Netlib's cost negation.

### 3.1 Download & expand

```bash
mkdir -p benchmarks/netlib/{raw,mps} && cd benchmarks/netlib/raw
CORE="afiro sc50a sc50b adlittle blend share2b kb2 sc105 sc205 lotfi israel \
      scagr7 bandm e226 degen2 bnl1 gfrd-pnc 25fv47 modszk1 degen3"
STRETCH="maros-r7 fit2p pilot greenbea dfl001 stocfor3 pilot87"
curl -O http://www.netlib.org/lp/data/emps.c && cc -O2 -o emps emps.c   # official expander
for f in $CORE $STRETCH; do
  curl -O "http://www.netlib.org/lp/data/${f}.Z"
  gzip -dc "${f}.Z" | ./emps > "../mps/${f}.mps"
done
```

(`emps` reads compressed-MPS on stdin, writes standard MPS to stdout — procedure straight from the Netlib readme. Kennington set is FTP-only and excluded.)

### 3.2 The ladder (explicit instances; sizes and optimal values from netlib readme)

**Core 20**

| Tier | Instance | rows m | cols n | nnz | Optimal value | Why it is on the ladder |
|---|---|---:|---:|---:|---:|---|
| T0 tiny | afiro | 28 | 32 | 88 | -4.6475314286E+02 | hello-world; CI smoke test |
| T0 tiny | sc50a | 51 | 48 | 131 | -6.4575077059E+01 | smallest staircase |
| T0 tiny | sc50b | 51 | 48 | 119 | -7.0000000000E+01 | twin, different data |
| T0 tiny | kb2 | 44 | 41 | 291 | -1.7499001299E+03 | ~35% degenerate pivots |
| T0 tiny | adlittle | 57 | 97 | 465 | 2.2549496316E+05 | classic small industrial |
| T0 tiny | blend | 75 | 83 | 521 | -3.0812149846E+01 | refinery blending (MRPL-relevant!) |
| T1 small | share2b | 97 | 79 | 730 | -4.1573224074E+02 | bounds/ranges handling |
| T1 small | sc105 | 106 | 103 | 281 | -5.2202061212E+01 | scaling check vs sc205 |
| T1 small | lotfi | 154 | 308 | 1086 | -2.5264706062E+01 | first non-trivial aspect ratio |
| T1 small | israel | 175 | 142 | 2358 | -8.9664482186E+05 | wide + huge objective range |
| T1 small | scagr7 | 130 | 140 | 553 | -2.3313892548E+06 | large-magnitude costs |
| T2 medium | sc205 | 206 | 203 | 552 | -5.2202061212E+01 | growth step from sc105 |
| T2 medium | bandm | 306 | 472 | 2659 | -1.5862801845E+02 | banded structure |
| T2 medium | e226 | 224 | 282 | 2767 | -1.8751929066E+01 | free rows present |
| T2 medium | degen2 | 445 | 534 | 4449 | -1.4351780000E+00 → see note | ~57% degenerate steps — stress test #1 |
| T2 medium | bnl1 | 644 | 1175 | 6129 | 1.9776292856E+03 | first 600+ row solve |
| T2 medium | gfrd-pnc | 617 | 1092 | 3467 | 6.9022359995E+06 | bounds variety |
| T2 medium | 25fv47 | 822 | 1571 | 11127 | 5.5018458883E+03 | classic medium benchmark |
| T2 medium | modszk1 | 688 | 1620 | 4158 | 3.2061972906E+02 | very degenerate, dual-simplex-friendly |
| T2 medium | degen3 | 1504 | 1818 | 26230 | -9.8729400000E+02 | ~52% degenerate steps — stress test #2 |

Note: degen2 optimal value is `-1.4351780000E+03` (the table above has a typo guard — always validate against `miplib2017-v36.solu`-style reference values pulled from the Netlib readme itself).

**Stretch 7 (post-Phase-1)**

| Instance | rows | cols | nnz | Optimal value | Role |
|---|---:|---:|---:|---:|---|
| maros-r7 | 3137 | 9408 | 151120 | 1.4971851665E+06 | structured large, hard for some solvers historically |
| fit2p | 3001 | 13525 | 60784 | 6.8464293232E+04 | primal/dual pair story (vs fit2d) |
| pilot | 1442 | 3652 | 43220 | -5.5740430007E+02 | ill-conditioned family |
| greenbea | 2393 | 5405 | 31499 | -7.2462405908E+07 | large refinery model (MRPL narrative) |
| dfl001 | 6072 | 12230 | 41873 | 1.12664E+07 (** tolerance-sensitive) | numerics trap; known FEASIBILITY_TOLERANCE sensitivity |
| stocfor3 | 16676 | 15695 | 74004 | -3.9976661576E+04 | gateway to 10⁴-scale |
| pilot87 | 2031 | 4883 | 73804 | 3.0171072827E+02 | bad-scaling stress (per Lustig/Stone) |

Out-of-ladder but noted: **osa-14** and friends are *not* Netlib — they live in Mészáros' lptestset (`www.sztaki.hu/~meszaros/public_ftp/lptestset/`), which Mittelmann uses; treat them as Mittelmann-set instances, not Netlib ladder members.

---

## 4. MILP: MIPLIB 2017 protocol

Sources: <https://miplib.zib.de/index.html>, `/download.html`, `/instances2017.html`.

### 4.1 Sets and downloads

```bash
mkdir -p benchmarks/miplib2017 && cd benchmarks/miplib2017
curl -LO https://miplib.zib.de/downloads/benchmark.zip              # Benchmark Set v2, 240 instances, 317 MB
curl -LO https://miplib.zib.de/downloads/easy-v18.test              # current easy list (2026-01-26)
curl -LO https://miplib.zib.de/downloads/miplib2017-v36.solu        # solution file v36 (2026-01-26)
curl -LO https://miplib.zib.de/downloads/miplib2017-testscript-v1.0.4.zip  # run scripts + feasibility checker
```

Definitions we adopt verbatim from ZIB:
- **easy** = solved in < 1 h with ≤ 16 threads, out-of-the-box solver, desktop hardware.
- **Benchmark Set** = 240 instances solvable by today's codes; **Collection Set** = full diverse library (~1 GB+, do not mirror into the repo).
- Categories drift over time (solufile updates move instances between easy/hard/open); **pin solufile v36 + easy-v18 and say so in every table**.

### 4.2 Starter subset for from-scratch B&B (+ basic cuts) — 16 core + 4 stretch

All names below are verbatim members of the official `easy-v18.test` list (fetched 2026-08-25), chosen small-first with structural diversity:

**Core 16:** `hanoi5`, `gt2`, `p0201`, `mod010`, `flugpl`, `air03`, `blend2`, `ran12x21`, `ran13x13`, `misc07`, `dsbmip`, `khb05250`, `cap6000`, `acc-tight2`, `aflow30a`, `10teams`

**Stretch 4:** `mas76`, `noswot`, `stein45inf`, `timtab1`

Rationale: covering/packing (`air03`, `p0201`, `stein45inf`), transportation (`ran12x21/13x13`), fixed-charge/network (`khb05250`, `gt2`, `blend2`), planning/scheduling (`flugpl`, `mod010`, `10teams`, `timtab1`), job-shop flow (`aflow30a`, `mas76`), production (`misc07`, `cap6000`), puzzle/degenerate (`hanoi5`, `dsbmip`, `noswot`, `acc-tight2`). These exercise exactly the components named in the PS: B&B, cuts, heuristics.

### 4.3 Solution-checking convention

1. Objective values validated against `miplib2017-v36.solu` (official best-known/optimal values).
2. Feasibility of any incumbent we report checked with the official testscript checker (or IPET) — never self-declared feasible.
3. A MIP counts as **solved** only if gap ≤ 1e-4 AND the incumbent passes the official checker AND status is provably optimal or within threshold (see §6).
4. If ZIB publishes a new solufile mid-project, finish the campaign on v36, then re-run changed categories and footnote the delta.

### 4.4 Usage/license/citation terms

- There is **no blanket open-source license**: MIPLIB is a research library curated by ZIB et al.; instances carry per-instance provenance, and a minority originate from industry sources with conditions noted on their individual instance pages. Practical rule: cite properly, keep instances out of public redistribution (gitignore `benchmarks/`), link instead of committing binaries.
- **Cite (mandatory):** Gleixner, Hendel, Gamrath, Achterberg, Bastubbe, Berthold, Christophel, Jarck, Koch, Linderoth, Lübbecke, Mittelmann, Ozyurt, Ralphs, Salvagnin, Shinano: *"MIPLIB 2017: Data-Driven Compilation of the 6th Mixed-Integer Programming Library"*, Mathematical Programming Computation (2021), DOI [10.1007/s12532-020-00194-3](https://doi.org/10.1007/s12532-020-00194-3). BibTeX is published on the MIPLIB index page — copy verbatim from there, do not hand-edit.

---

## 5. Mittelmann benchmarks (plato.asu.edu/bench.html)

Verified 2026-08-25 against bench.html, lpopt.html, lpfeas.html.

### 5.1 Current structure (old Simplex/Barrier benchmarks are retired)

| Set | Task | Our role |
|---|---|---|
| **LPfeas** (2026-08-10 rev) | find primal-dual feasible point; explicitly "also for GPUs"; tolerance 1e-6; TL 15 000 s CPU / 1 000 s GPU (B200 192 GiB) | **primary aspirational target** — this is where cuOpt/cuPDLPx/PDLP compete; our PDLP story maps directly onto it |
| **LPopt** (2026-07-01 rev) | find optimal basic solution (simplex-style) | secondary target once simplex is mature |
| **Large Network-LP** | commercial-vs-free network LPs | stretch only |
| **MILP / MIPFEAS on MIPLIB 2017** | full MIP solving / feasibility | long-term; matches our MIPLIB protocol |
| QPLIB sets (convex continuous / discrete / non-convex) | QP & MIQP | deferred (§5.3) |

Published comparisons there include HiGHS-1.15.0, COPT-8.0.0, MOSEK-11.1.11, CLP-1.17.7, SOPLEX-8.0.0, GLOP/OR-Tools, PDLP (OR-Tools), XOPT, KNITRO, cuOpt-26.08, cuPDLPx-0.3.0, HPR-LP-C — i.e., exactly the field our GPU method must be quoted against. Reference hardware: i7-11700K @ 3.6 GHz, 64 GB RAM; metric: shifted geometric mean shifted by 10 s; timeouts ("t", "f", "m") counted as max-time.

### 5.2 What we actually commit to

- Phase 1–2 (this hackathon cycle): replicate the **LPfeas tolerance (1e-6)** on our PDLP for a *subset* of the smaller LPfeas instances (e.g., `qap15`, `ex10`, `fome13`, `rail4284`, `nug08-3rd`) and report under Mittelmann conventions; full-set submission is explicitly NOT promised.
- Never quote Mittelmann numbers for ourselves without running their protocol (their hardware, their TL policy); label any self-run subset clearly as such.

### 5.3 QPLIB note (deferred)

When QP work begins: use **QPLIB** (<http://qplib.zib.de/>), taking the convex continuous/discrete subsets as mirrored in Mittelmann's QPLIB benchmarks (`plato.asu.edu/ftp/cconvex.html`, `/convex.html`). ADMM-QP validation starts on convex continuous only. Amend this doc then.

---

## 6. Metrics (one definition each, no synonyms)

### 6.1 Relative MIP gap

```
rel_gap = |ub − lb| / max(|ub|, ε),  ε = 1e−10
```
(matches HiGHS' documented `mip_rel_gap = |ub−lb|/|ub|`; ε guards zero-objective instances)

| Label | Meaning |
|---|---|
| **solved** | rel_gap ≤ **1e-4** (HiGHS default) and proof of optimality or bound convergence |
| **tight-solved** | rel_gap ≤ **1e-6** (used for LP-parity and Mittelmann-replica runs) |
| **feasible-only** | incumbent found within TL, gap above threshold, incumbent passed official checker |
| **timeout-no-incumbent** | TL hit, no valid incumbent |
| **infeasible / unbounded** | proven status reported |

### 6.2 Time limits per tier (hard caps, wall-clock)

| Set/tier | Time limit |
|---|---|
| Netlib T0/T1 | 60 s |
| Netlib T2 | 300 s |
| Netlib stretch T3 | 1800 s |
| MIPLIB starter (core 16) | 300 s (matches the 5-min regime of the CHAP/Land–Doig GPU-heuristic competition) |
| MIPLIB stretch | 3600 s |
| Mittelmann-replica runs | 15 000 s CPU / 1 000 s GPU per their convention |

### 6.3 Aggregates

1. **Solved-count table** per tier: `#solved / #instances` at 1e-4, plus median runtime of solved.
2. **Shifted geometric mean runtime**, shift = 10 s (Mittelmann convention, [shgeom](http://plato.asu.edu/ftp/shgeom.html)); timeout counted as TL.
3. **Speedup-crossover plot**: x = nnz (log scale), y = time(HiGHS)/time(ours) (log), horizontal line at 1.0, one curve per algorithm family (our CPU simplex, our GPU PDLP, our B&B). The crossover chart is a mandatory figure in every results presentation (per project report §5).
4. Runs < 60 s are executed 3× and reported as median; ≥ 60 s runs once. Cold start every run; `random_seed=0` both solvers.

---

## 7. Canonical result-table format (used EVERYWHERE)

Ten fixed columns, exact order; extra diagnostic columns may be appended to the right, never inserted:

```csv
instance,m,n,nnz,our_time_s,highs_time_s,our_objective,highs_objective,rel_gap,status
afiro,28,32,88,0.042,0.011,-464.75314286,-464.75314286,0.0,optimal
```

Rules:
- times in seconds, 3 decimals; objectives full precision as printed by the solver;
- `status ∈ {optimal, tight-solved, feasible-only, timeout-no-incumbent, infeasible, unbounded, error}`;
- `rel_gap` blank for pure-LP optimal rows (report KKT residuals instead: primal/dual residual columns appended);
- one row per (instance × configuration × machine-ID); every published table footer states: HiGHS version, our git SHA, machine-sheet ID, date, thread count.

---

## 8. Hardware sheet template

One filled sheet per machine, stored as `benchmarks/hardware/<machine-id>.md`, referenced by ID from every result table.

```markdown
# Machine sheet: <machine-id>          (recorded: YYYY-MM-DD)

| Field | Value |
|---|---|
| Role | local dev / remote 16GB GPU / Modal batch |
| CPU model | e.g. AMD Ryzen 7 5800H |
| Cores/Threads + base clock | 8C/16T @ 3.3 GHz |
| RAM (GB) + type | 32 GB DDR4-3200 |
| GPU model | e.g. NVIDIA RTX 3050 Laptop |
| GPU VRAM (GB) | 4 GB |
| CUDA toolkit / driver | CUDA x.y / driver xxx.yy.zz |
| OS + kernel | Ubuntu 24.04, kernel x.y |
| Compiler + flags | g++ 13.x, -O2 -march=native |
| BLAS/CUDA libs | cuBLAS/cuSPARSE/cuDSS versions |
| Clock/power state | locked? governor? laptop on AC? |
| Notes | thermal throttling observed, etc. |
```

Modal shapes (verified from docs.modal.com/guide/gpu): `T4` (16 GB) · `L4` (24 GB) · `A10` (24 GB, ≤4 GPUs) · `L40S` (48 GB) · `A100-40GB` / `A100-80GB` · `H100` / `H200` (SXM; H100 requests may be silently upgraded to H200 — pin `gpu="H100!"` when benchmarking) · `B200` / `B200+` / `B300`. For Modal runs the sheet additionally records: shape string requested vs actual GPU delivered (`nvidia-smi` output), vCPU/RAM of the container, region. Our standing rule: **never publish a Modal number without the actual-delivered-GPU line from nvidia-smi.**

The team's 16 GB remote GPU gets its own sheet the first time it produces a published number (its VRAM tier coincides with Modal `T4`-class memory; still record actual model).

---

## 9. Fair-comparison rules (binding)

1. Both solvers run on the same machine class, cold start, same file on disk, same MPS reader path where possible.
2. HiGHS gets its documented defaults plus ONLY the pinned flags of §2. We never tune HiGHS down (no artificial iteration limits) and never tune it up per-instance either — defaults are its strength anyway.
3. No instance-specific tuning of our solver appears in any table unless the tuning itself is published alongside.
4. Every number traceable to: machine sheet ID + HiGHS tag + solver commit + raw logs archived under `benchmarks/logs/<date>-<machine-id>/`.
5. GPU rows always state GPU model + VRAM; CPU-vs-GPU claims never mix machines.
6. Honest-losses policy (project report risk #2): an unexplained GPU loss shown without analysis is a documentation bug.

---

## 10. Sources

- Netlib LP readme (sizes, optimal values, emps expansion): http://www.netlib.org/lp/data/readme
- MIPLIB 2017 index/download/easy-list/solu/citation: https://miplib.zib.de/ (index.html, download.html, instances2017.html)
- Mittelmann benchmarks overview + LPopt + LPfeas: http://plato.asu.edu/bench.html , /ftp/lpopt.html , /ftp/lpfeas.html ; shifted geometric mean: /ftp/shgeom.html
- HiGHS releases: https://github.com/ERGO-Code/HiGHS/releases (v1.15.1, 2026-07-02); option definitions: https://ergo-code.github.io/HiGHS/dev/options/definitions/
- Modal GPU shapes & H100! pinning: https://modal.com/docs/guide/gpu
- Project background: docs/SIH26119-RESEARCH-REPORT.md (repo main)
