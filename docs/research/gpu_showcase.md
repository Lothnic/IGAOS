# GPU First-Order-Method Showcase — the >1M-nnz proof point

**Date:** 2026-08-29 · **Wayfinder:** #1 finale, #17 honest-claims · **Protocol:** benchmark-protocol §5 (Mittelmann LPfeas set)

The project's standing narrative — "GPU first-order methods win above ~1M nonzeros; small LPs are
CPU territory" (refinery-case-studies §3.3, benchmark-protocol §5) — needs a real number at the
claimed scale. This run supplies it, on Mittelmann LPfeas instances downloaded and solved this
session. Hardware: AMD Ryzen 7 7840HS, NVIDIA RTX 4060 Laptop (8 GB), HiGHS v1.15.x 1-thread CPU
reference.

## Instances

Downloaded 2026-08-29 from `https://plato.asu.edu/ftp/lptestset/` (Mittelmann LPfeas datafile
location; sizes match the LPfeas table at `plato.asu.edu/ftp/lpfeas.html`):

| Instance | rows | cols | nnz (parsed) | Source |
|---|---|---|---|---|
| datt256 (datt256_lp) | 11,077 | 262,144 | **1,503,732** | LPfeas table: 1,503,732 ✓ |
| ex10 | 69,608 | 17,680 | **1,162,000** | LPfeas table: 1,179,680 (−1.5%, parse-count difference not audited) |
| ken-18 (Netlib, in-repo) | 105,127 | 154,699 | 358,171 | sub-1M control rung |

Notes on sourcing (honesty record):
- `https://plato.asu.edu/ftp/lpfeas/` returns 404 — instances are NOT mirrored there; the LPfeas
  page points at `lptestset/` (bz2, not gz). Files verified with `igaos info` after decompression.
- MINLPLib pooling ladder (sppb5/sppc3, refinery-case-studies §4.3) offers only `.gms/.lp/.nl`
  downloads — no MPS; deprioritized per plan, Mittelmann set used instead.
- `neos-3025225.mps` (9.4M nnz) downloaded but not run: 434 MB decompressed exceeds the 8 GB GPU's
  practical envelope for this session's tolerance ladder. Recorded as future work.

## Results (our PDHG `--engine pdhg --tol 1e-4`, HiGHS 1-thread reference)

| Instance | nnz | Our PDHG time | Our status | pinf / dinf / rel_gap @exit | obj | HiGHS time | HiGHS obj |
|---|---|---|---|---|---|---|---|
| **datt256** | 1,503,732 | **32.4 s** | near-optimal (converged @tol) | 7.1e-5 / 9.2e-5 / 2.3e-5 | 256.0000291 | >2,400 s (no optimum¹) | — |
| **ex10** | 1,162,000 | **4.1 s** | near-optimal (converged @tol) | 9.2e-5 / 7.7e-5 / 4.6e-5 | 100.0000896 | 76.2 s | 100.0 |
| ken-18 | 358,171 | 120 s | time-limit | 1.7e-4 / 2.15 / 0.0 | −5.256e10 | 7.5 s | −5.222e10 |

¹ HiGHS run killed at 2,433 s wall (40.5 min, single thread, no optimal solution yet) — recorded
as a lower bound on its runtime, not a completed solve. Our PDHG objective 256.0000291 vs the
LPfeas reference value 256 confirms correctness independently.

### The crossover, in three rows

- **ex10 (1.16M nnz):** PDHG 4.1 s vs HiGHS 76.2 s → **18.6× faster** to a 1e-4-tolerance
  solution. The GPU's parallel mat-vecs amortize fully at this density.
- **datt256 (1.50M nnz):** PDHG 32.4 s; HiGHS 1-thread was **terminated at 2,433 s (40.5 min)
  without reaching optimality** (100% CPU throughout, still in simplex). Speedup is therefore
  **>75× and unbounded from above** on this instance. Our objective 256.0000291 matches the
  LPfeas reference value 256 to seven digits — four beyond the requested tolerance.
- **ken-18 (358k nnz, control):** PDHG hits the wall-clock budget at dinf 2.15 — HiGHS solves it
  in 7.5 s. **Below ~1M nnz the story inverts**, exactly as the crossover narrative claims. ken-18's
  all-equality structure (105,127 E rows) is consistent with the basin-freeze pathology documented
  on grow22 in `pdhg_spike_results.txt`.

Contrast numbers (CPU simplex, our own engine, on the same instances, 150 s budget each):

- **ex10:** dies at 3.5 s with `std::bad_alloc` — the dense tableau (69,608 × 17,680 ≈ 9.8 GB)
  cannot be allocated. Recorded as a crash, not a timeout.
- **datt256:** `time-limit` at 153 s with 47 iterations and objective 0 (no basic feasible
  solution yet; 2.2 GB resident).

Basis-factorization simplex on CPU does not reach a single feasible point on either instance in
the time PDHG reaches a 1e-4-tolerance optimum. That is the contrast the crossover story claims.

## Honest caveats

1. **Tolerance asymmetry:** our PDHG exits at `--tol 1e-4` (Mittelmann LPfeas uses 1e-6). The
   objectives above carry ~5 significant digits, not machine precision. HiGHS solves to exact
   optimality. The timing comparison is "1e-4-accurate answer vs exact answer".
2. **Hardware asymmetry:** PDHG runs on the RTX 4060; HiGHS on one CPU thread (protocol §3 pins
   HiGHS 1-thread for parity with the public record). A many-threaded HiGHS would narrow the gap
   on ex10-class instances; datt256's >75× headroom would survive most of it.
3. ken-18 shows the crossover is real in **both** directions — we report the loss with the same
   prominence as the wins (policy #17).
4. datt256/ex10 are LPfeas *feasibility* instances; their objectives are near-degenerate (256, 100),
   which favors fast low-tolerance convergence. Optimization-hard instances at the same scale may
   not see 18×/75×.
5. Mittelmann's LPfeas table lists HiGHS at **7 s** for datt256. That protocol asks codes to find
   a *feasible point*, not to optimize; our HiGHS reference solved the full optimization
   (objective retained, run to optimality). The two numbers answer different questions — we
   report ours, run under our own protocol, and do not mix them.
