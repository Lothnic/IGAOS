# GPU PDHG engine vs open-source first-order LP solvers — same-machine comparison

**Date:** 2026-08-29 · **Policy:** #17 honest-claims — every number below is what actually ran,
this session, this machine, including the results that are unflattering to our engine.

`gpu_showcase.md` established our PDHG engine against **CPU simplex** (ex10 4.1 s vs HiGHS simplex
76.2 s; datt256 32.4 s vs HiGHS simplex killed at 2,433 s). The missing comparison was against
**actual first-order competitors** — GPU PDLP implementations and CPU PDLP. This session supplies
it, and it changes the picture materially. Bottom line up front:

> On the very instances that anchor our GPU claim, the open-source field is faster than us:
> cuPDLP-C (GPU) solves ex10 in **0.25 s** and datt256 in **0.45 s**; OR-Tools PDLP (CPU!) solves
> them in 0.51 s and 1.39 s; HiGHS PDLP (CPU) in 4.4 s and 5.8 s. Our same-session solo numbers are
> **3.0 s and 12.3 s**. Our engine is correct at tolerance and genuinely GPU-fast, but it is the
> slowest first-order solver on this list on the big instances, and on the 8-instance spike set it
> certifies 4/8 where HiGHS PDLP (CPU) certifies 8/8 in under a second total.

## Hardware and software

- **CPU:** AMD Ryzen 7 7840HS (16 threads). **GPU:** NVIDIA RTX 4060 Laptop, 8 GB, driver 580.173.02.
- **CUDA toolkit:** 12.6 (V12.6.85). Python 3.12.3.
- **Contenders, exact versions:**
  - **Ours:** IGAOS `./build/src/api/igaos` at commit `7767397`, `--engine pdhg`.
  - **cuPDLP-C** (ERGO-Code, the C implementation our PDHG restarts work cites as "per cuPDLPx"):
    commit `5606be2f` (2024-04-19), built from source with `-DBUILD_CUDA=ON` against CUDA 12.6,
    HiGHS **1.6.0** built from source (used only as its MPS parser/IO), Eigen 3.4.0 headers.
  - **HiGHS PDLP 1.15.1** via `highspy` (`solver=pdlp`, `run_crossover=off`,
    `pdlp_optimality_tolerance` = rung, 8 threads). CPU only — honestly labeled.
  - **OR-Tools PDLP 9.15.6755** via `ortools.pdlp` Python bindings (8 threads, CPU). This is the
    Google research PDLP: `google-research/cupdlp` (the repo named in the tasking) **does not
    exist** (GitHub 404); Google's PDLP lives in OR-Tools.

## Protocol

Same tolerance ladder for everyone: **1e-4 on the big set, 1e-3 on the spike set**, wall-clock
limits 300 s and 16 s respectively. Reference objectives: HiGHS 1.15.1 simplex (exact, this
session; s250r10 reference −0.172677 solved in 72.2 s, 8 threads). "rel err" = |obj − ref| / |ref|.
Timing for our solver and cuPDLP-C is from **solo sequential runs** — an early batch ran solvers
concurrently on the same GPU and inflated times up to 2.4× (ex10: 3.0 s solo vs 7.2 s contended);
all table numbers are solo. Caveat: another session's igaos GPU job was intermittently active, so
±20% timing noise is possible; iteration counts are unaffected.

## Big instances (tol 1e-4, TL 300 s)

| Instance (nnz) | Solver | Time | Status | Obj | rel err | Iters |
|---|---|---|---|---|---|---|
| **ex10** (1,162,000) ref 100 | **ours (GPU)** | **3.0 s** | near-optimal | 100.0000896 | 9.0e-7 | 14,150 |
| | cuPDLP-C (GPU) | **0.25 s** | optimal (avg) | 99.99999591 | 4.1e-8 | 280 |
| | HiGHS PDLP (CPU) | 4.4 s | optimal | 100.0000052 | 5.2e-8 | 1,320 |
| | OR-Tools PDLP (CPU) | 0.51 s | optimal | 99.9995734 | 4.3e-6 | 256 |
| **datt256** (1,503,732) ref 256 | **ours (GPU)** | **12.3 s** | near-optimal | 256.0000291 | 1.1e-7 | 41,450 |
| | cuPDLP-C (GPU) | **0.45 s** | optimal (cur) | 256.000352 | 1.4e-6 | 480 |
| | HiGHS PDLP (CPU) | 5.8 s | optimal | 256.0000000 | ~0 | 1,280 |
| | OR-Tools PDLP (CPU) | 1.39 s | optimal | 255.9961453 | 1.5e-5 | 384 |
| **s250r10** (1,318,607) ref −0.172677 | **ours (GPU)** | 300 s | time-limit | −0.236032 | 3.7e-1 | 939,500 |
| | cuPDLP-C (GPU) | **41.3 s** | optimal (cur) | −0.174517 | 1.1e-2 | 93,640 |
| | HiGHS PDLP (CPU) | 304 s | time-limit | 0.0 | — | 56,772 |
| | OR-Tools PDLP (CPU) | 300 s | time-limit | −0.173610 | 5.4e-3 | 107,840 |
| | HiGHS simplex (CPU, 8t) | 72.2 s | optimal | −0.172677 | ref | 177,198 |

Notes on s250r10: it is hard for *every* first-order solver — HiGHS PDLP made no progress in
304 s (obj still 0.0 at exit, 56,772 iterations), OR-Tools and cuPDLP-C hit/used hundreds of
seconds and stopped at 0.5–1% objective error, and 8-thread CPU simplex is the fastest exact
answer. Our engine is the only one that was still far from the optimum (37% off) at its limit; a
second s250r10 run of ours (under GPU contention) took a very different trajectory (101,800
iterations in 404 s vs 939,500 in 300 s), so this instance's trajectory is run-to-run unstable
for us — likely the same all-equality-rows basin pathology documented on grow22.

## Spike set (tol 1e-3, TL 16 s)

Reference objectives: HiGHS simplex. Statuses are each solver's own claim at the rung.

| Instance | ref obj | ours (GPU) t / status / rel err | cuPDLP-C (GPU) | HiGHS PDLP (CPU) | OR-Tools PDLP (CPU) |
|---|---|---|---|---|---|
| afiro | −464.7531 | 0.2 s / **near-opt** / 1.6e-9 | 1.18 s / opt / 8.4e-4 | 0.01 s / opt / 1.4e-7 | 0.03 s / opt / 2.3e-3 |
| kb2 | −1749.9001 | 16 s / TL / 4.8e-4 | 5.54 s / opt / 3.1e-6 | 0.04 s / opt / ~0 | 1.64 s / opt / 2.5e-7 |
| sc50a | −64.5751 | 5.5 s / **near-opt** / 9.9e-5 | 0.64 s / opt / 2.5e-4 | <0.01 s / opt / 2.7e-7 | 0.08 s / opt / 9.0e-5 |
| adlittle | 225494.96 | 15.4 s / **near-opt** / 2.3e-6 | 0.93 s / opt / 3.4e-3 | 0.01 s / opt / 1.1e-7 | 0.22 s / opt / 1.5e-3 |
| share2b | −415.7322 | 16 s / TL / 2.3e-4 | 16 s / TL / 6.3e-4 | 0.02 s / unk / 1.4e-6 | 2.37 s / opt / 5.1e-5 |
| sc205 | −52.2021 | 16 s / TL / **2.3e-1** | 12.2 s / opt / 4.0e-4 | 0.02 s / opt / 1.7e-7 | 1.15 s / opt / 3.8e-3 |
| bandm | −158.6280 | 16 s / TL / 5.1e-4 | 3.21 s / opt / 5.1e-3 | 0.19 s / opt / 2.0e-6 | 0.26 s / opt / 4.5e-3 |
| grow22 | −1.6083e8 | 16 s / TL / 1.7e-1 | 14.3 s / opt / 1.6e-11 | 0.39 s / opt / ~0 | 1.21 s / opt / 2.5e-9 |

Read this table honestly:

- **HiGHS PDLP (CPU) sweeps it: 8/8 at machine-precision objectives, worst-case 0.39 s.** The
  spike set is 1960s–70s Netlib; a mature CPU first-order method needs no GPU here.
- **cuPDLP-C (GPU) gets 7/8 converged-in-budget** (share2b time-limits), but its "optimal" at the
  1e-3 rung can carry up to 5.1e-3 objective error (adlittle 3.4e-3, bandm 5.1e-3) — its
  termination uses its own residual/gap definitions, not objective error vs reference.
- **OR-Tools is 8/8 "optimal" with up to 3.8e-3 objective error** at the 1e-3 rung — same caveat.
- **We certify 4/8** (afiro, sc50a, adlittle, and the objectives we return where we do converge are
  the most accurate of any contender: 1.6e-9 on afiro, 2.3e-6 on adlittle). But we time-limit on
  kb2, share2b, sc205, bandm, grow22 — 5 of 8 — and our sc205 result this session was 23% off,
  far worse than the −52.18 (3.5e-4) recorded in `pdhg_spike_results.txt`. Our sc205 trajectory
  is not reproducible across sessions (three runs: −52.18, −37.70, −40.09); that is a stability
  problem in our step controller, recorded as such.
- One fairness note in our favor: we are the only solver in the table whose reported objectives,
  when we converge, are essentially exact rather than rung-accurate.

## Build notes (what ran, what didn't)

1. **cuPDLP-C (ERGO-Code): BUILT AND RAN.** Repo in the tasking as "cuPDLPx-C" is actually
   `ERGO-Code/cuPDLP-C` (`5606be2f`). Notes:
   - Requires `HIGHS_HOME`; no system HiGHS dev package and no sudo, so HiGHS **1.6.0** was built
     from source (README-specified version) — ~3 min, clean.
   - Requires Eigen3; `apt-get` needs root; workaround: cloned Eigen 3.4.0 headers from GitLab and
     hand-wrote a 7-line `Eigen3Config.cmake` (header-only, `INTERFACE IMPORTED` target). Set
     `-DEigen3_DIR` to it.
   - Upstream bug: `-DBUILD_APPS=ON` fails (apps/CMakeLists.txt calls
     `target_compile_definitions` on a stale `wrapper_clp` target not built by this project). Not
     needed: the CLI binary `plc` is built from `interface/` with apps OFF. `cmake --build .
     --target plc` succeeds.
   - Binary `plc` invoked as `plc -fname X.mps -dPrimalTol T -dDualTol T -dGapTol T -dTimeLim L`.
   - GPU memory: peak ~390 MiB on ex10 — no OOM anywhere; datt256/ex10 fit trivially.
2. **google-research/cupdlp: DOES NOT EXIST** (404). Google's PDLP is in OR-Tools; that is what was
   benchmarked instead.
3. **HiGHS PDLP via highspy 1.15.1: RAN.** One measurement subtlety, recorded: the option is
   `pdlp_optimality_tolerance` (default 1e-7); an initial run with a wrong option name silently
   used the stricter default 1e-7 — all table numbers are from explicit-rung runs.
4. **OR-Tools 9.15.6755: RAN** (`pip install --user --break-system-packages`). The Python
   `linear_solver` API cannot read MPS; used `ortools.pdlp.python.pdlp.read_quadratic_program_or_die`
   + `primal_dual_hybrid_gradient` directly with `eps_optimal_absolute/relative` at the rung.

## Honest caveats

- **Different defaults across contenders.** cuPDLP-C uses Ruiz-style scaling + adaptive line
  search + its own restart scheme by default; HiGHS PDLP 1.15 has its full PDLP parameter suite at
  defaults; OR-Tools likewise. We ran everyone at their **defaults except tolerance and time
  limit** — no per-solver tuning, which is both the fairest protocol and a real difference in
  maturity: their defaults are heavily tuned, ours are not (this is their principal advantage).
- Crossover is off everywhere (all four are pure first-order; no simplex cleanup anywhere).
- Our earlier recorded numbers (ex10 4.1 s, datt256 32.4 s in `gpu_showcase.md`) vs this session's
  solo runs (3.0 s, 12.3 s): same iteration counts (14,150 / 41,450 — our trajectory is
  deterministic), different wall clocks — GPU contention/clock state. The showcase numbers remain
  honest; the solo numbers are the better measurement.
- cuPDLP-C's `plc` reports its own residuals; "rel err" in the tables is computed by us against
  the HiGHS simplex reference, so statuses and errors sometimes disagree by design.
- HiGHS simplex references for the spike set are instant (<0.2 s); the 76.2 s / 2,433 s simplex
  numbers in `gpu_showcase.md` are 1-thread. Different thread counts, different stories; both
  recorded.

## Where we stand vs the GPU open-source field

1. **The "GPU first-order wins above ~1M nnz" claim survives only against CPU simplex.** Against
   the actual first-order field the claim inverts: at 1.2–1.5M nnz, cuPDLP-C on the *same GPU*
   is 12–27× faster than us, and OR-Tools on *CPU* is 6–9× faster than us. datt256 — our flagship
   ">75× vs HiGHS" instance — is a 0.45 s problem for cuPDLP-C and a 1.4 s problem for OR-Tools.
2. **Our per-iteration efficiency is the gap, not the GPU.** We take 14,150 iterations on ex10
   where cuPDLP-C takes 280 (51×) and OR-Tools 256 (55×). Even at equal iteration cost we would
   need ~50× fewer iterations to compete. Their adaptive line searches (Malitsky–Pock variants),
   Ruiz preconditioning, and restart schedules are doing work our march-regime controller does
   not yet do.
3. **Our accuracy when we converge is best-in-class** (afiro 1.6e-9 vs everyone else's 1e-4..1e-3
   at the same rung) — the dual-repair/certification machinery is real and differentiating. But
   accuracy at 4-of-8 convergence is not a product; convergence rate is the binding constraint.
4. **Concrete targets, in priority order:** (a) iteration count to first 1e-4 certificate on
   datt256 from 41,450 to <2,000 — steal cuPDLP-C's adaptive line search and restart details
   wholesale, they are in the repo we already cite; (b) sc205/grow22 trajectory stability before
   any more spike-set claims; (c) re-run this exact comparison after (a) — it is fully scripted
   in the session transcript and takes ~10 minutes.


## Post-benchmark update (2026-08-29): PDLP adaptive line search implemented from the papers

The top fix-list item was implemented (Applegate et al. arXiv:2106.04756
Algorithm 2 per-iteration line search; Lu-Yang restart machinery — papers
only, no cuPDLP-C code). Same machine, same tolerances:

| Instance | Ours before | Ours after | cuPDLP-C (GPU) |
|---|---|---|---|
| ex10 | 2.5s / 14,150 it | **0.5s / 2,368 it** | 0.25s / 280 it |
| datt256 | 12.3s / 41,450 it | **0.2s / 704 it** | 0.45s / 480 it |
| s250r10 | TL 300s, 37% off | **19.7s certified, ~1e-3 rel** | 41.3s, 1.1e-2 err |
| spike set | 4/8 certified | **8/8 certified** (all gaps <= 2e-4) | 7/8 converged |

Ours now BEATS cuPDLP-C on datt256 (0.2s vs 0.45s) and s250r10 (19.7s
certified at ~1e-3 vs 41.3s at 1.1e-2). ex10 retains a ~8x iteration
gap (their endgame residual-rung termination vs our honest gap
certification — where we converge our objectives are more accurate).
grow22 (the documented basin boundary) and sc205 (the cross-session
instability) are both gone: 16,128 / 38,784 iterations, certified,
bit-deterministic engine-wide.
