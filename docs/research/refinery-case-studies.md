# Refinery Demo Case Studies — Selection Report

**Wayfinder ticket:** Lothnic/IGAOS#10 (parent map #1) · **Project:** SIH26119 Indigenous GPU-Accelerated Optimization Solver · **Sponsor:** MRPL
**Date:** 25 August 2026 · **Scope:** research only — select one crude-blending/refinery-scheduling **LP** and one pooling/blending **QP** demo case from open literature, with fully available numeric data.

---

## 1. Selection criteria

A case study is adoptable only if it passes all five gates:

| # | Gate | Why |
|---|------|-----|
| 1 | **Complete numeric data freely downloadable** (no book-only or paywalled data) | We must encode and solve the exact published instance; evaluators must be able to reproduce it |
| 2 | **Reads as refinery-real to MRPL** (crude selection, distillation yields, octane/sulfur specs, blending) | The PS explicitly names "refinery scheduling, crude blending, production planning" as target domains |
| 3 | **Known reference optimum** | Enables honest "our solver vs HiGHS/Gurobi" tables required by the PS benchmark clause |
| 4 | **Ladder of scales** (tiny hand-checkable → 100s–1000s variables) | PS demands thousands-to-millions variable capability; a single fixed-size toy cannot demo scale |
| 5 | **Fits solver roadmap** LP (simplex/PDLP) and QP (ADMM) modules from `docs/SIH26119-RESEARCH-REPORT.md` §3 | Demos must exercise the engines we are actually building |

---

## 2. Candidate comparison

Sizes for pooling instances below are exact counts pulled from MINLPLib's machine-readable instance database (`minlplib.org/instancedata.csv`, PQ/TP reformulation rows); optima are MINLPLib primal bounds.

### 2a. Crude-blending / refinery-scheduling LP candidates

| Candidate | Type / size | Data availability | Refinery-realism story | Verdict |
|---|---|---|---|---|
| **Williams "Refinery Optimization"**, *Model Building in Mathematical Programming*, 5th ed., example #6 (pp. 258, 306–310), open implementation: [Gurobi/modeling-examples/refinery](https://github.com/Gurobi/modeling-examples/tree/master/refinery) (`refinery.ipynb`) | LP, ~29 constraints + bounds, ~30–40 vars; optimal profit **$211,365.13** stated in notebook | **Complete**: all parameters inline in Python in an open GitHub notebook (crude purchase limits, distillation splitting coefficients, reforming/cracking yields, octane numbers, vapor-pressure constants, product profits) | Full refinery flow: buy crudes → distillation splits → naphtha reforming → oil cracking → lube oil (0.5×residuum) → blend premium/regular/jet/fuel-oil/lube under octane ≥94/≥84, vapor pressure, premium ≥ 40% of regular | ✅ **FINAL LP PICK** — canonical "how to run an oil refinery", zero-friction reproduction |
| DOcplex `oil_blending` ([IBM docs](https://ibmdecisionoptimization.github.io/docplex-doc/mp/oil_blending.html)) | Tiny LP, 9 blend vars + 3 advert vars; full numpy data arrays printed on page | Complete inline | Buy 3 crudes, blend into super/regular/diesel under octane & lead specs, maximize revenue − costs | 🥈 runner-up — ideal 10-line "hello world" warm-up before the main demo |
| Food Manufacture I/II ([Gurobi/modeling-examples/food_manufacturing](https://github.com/Gurobi/modeling-examples/tree/master/food_manufacturing)) | MILP, multi-period blending with storage/hardening | Complete open notebooks | Edible-oil refining/blending over months — structurally identical to multi-period crude blending | 🥉 runner-up — natural MILP extension once branch-and-bound lands |
| Netlib `lp/infeas/refinery` ([Netlib](https://netlib.org/lp/infeas/), mirror: [SuiteSparse LPnetlib/lpi_refinery](https://sparse.tamu.edu/LPnetlib/lpi_refinery)) | LP **323 rows × 464 cols, 1,626 nonzeros** (T. Baker, Chesapeake Decision Sciences, 1993) | Complete MPS file on Netlib | Medium petrochemical-plant model, deliberately doctored infeasible (quality/volume) | ⭐ bonus pick — perfect presolve/infeasibility-diagnostics (IIS) demo; PS stresses numerical robustness |
| Alabi & Castro integrated refinery planning (IRP) ([Optimization Online 2008](https://optimization-online.org/2008/09/2080/)) | LP up to **1M constraints × 2.5M variables × 59M nonzeros** (2–300 day horizon) | Paper/report public; **instance files not downloadable** → fails gate 1 as a primary instance | Industrial-scale horizontal integration, crude purchase → product distribution | ❌ as instance, ✅ as citation proving multi-period refinery LPs reach millions of vars |
| EMRPS refinery-planning benchmark ([GitHub EMRPS/refinery-planning-benchmark](https://github.com/EMRPS/refinery-planning-benchmark), paper [arXiv:2503.22057](https://arxiv.org/html/2503.22057v1)) | NLP/MINLP; Case 1: 2 CDUs × 25 crudes × 25 units × 24 products; Case 3: multi-period | **Complete**: all sets/parameters as `.gdx`/`.txt`/`.gms` in repo | Real-world refinery-petrochemical complex, perturbed industrial data | 🏅 stretch goal — adopt after NLP/MINLP extensions (beyond current LP+QP scope) |

### 2b. Pooling / blending QP-NLP candidates

All instances downloadable from [MINLPLib](https://minlplib.org/) in `gms/lp/mod/nl/py` formats; original source collection: Alfaki & Haugland, ["Standard Pooling Problem Instances"](http://www.ii.uib.no/~mohammeda/spooling/) (UiB), referenced in *Strong formulations for the pooling problem*, J. Global Optimization 56:897–916, 2013.

| Instance family | Vars / cons (PQ form) | Bilinear terms | Known optimum (min form) | Story | Verdict |
|---|---|---|---|---|---|
| **Haverly 1/2/3** (Haverly 1978, ACM SIGMAP Bull. 45:9–28) | 10 / 13 each | 4 | **−400 / −600 / −750** | The canonical pooling problem: one pool blends two sulfur-graded feeds into two spec-bound products | ✅ **FINAL QP PICK** (canonical tier) |
| Ben-Tal 4 / 5 (Ben-Tal, Eiger, Gershovitz, Math. Program. 63:193–212, 1994) | 13 / 16 and 92 / 86 | 6 / 60 | −450 / −3500 | Q-formulation origin family; bt5 already "hundreds-scale" | ✅ included in ladder |
| Adhya 1–4 (Adhya, Tawarmalani, Sahinidis, I&EC Res. 38:1956–1972, 1999) | 33–58 / 49–77 | 20–40 | −549.80 … −877.65 | Aviation-gasoline blending instances — literally fuel blending | ✅ included in ladder |
| Foulds 2–4 (Foulds, Haugland, Jörnsten, Optimization 24:165–180, 1992) | 36→**672 / 571** | 16→512 | −1100 / −8 | First genuinely large classical family | ✅ scale rung (~700 vars) |
| SPP random families sppa0→sppc3 (Alfaki–Haugland generator) | 500→**10,567 / 10,876** | 329→9,116 | −35,812→−122,383 | Modern large-scale standard pooling problems | ✅ scale ceiling (10⁴ vars) demos GPU first-order methods |

Pooling problem status note: the pooling problem is **NP-hard** (bilinear quality-balances; Alfaki & Haugland 2013; also stated in Ceccon & Misener, [arXiv:2105.01687](https://arxiv.org/abs/2105.01687)). It is therefore presented honestly at three levels rather than mislabeled as a pure QP:

1. **LP level** — Haverly with pool-quality bounds relaxed/fixed = plain LP relaxation (weak but linear; fine for simplex/PDLP demos);
2. **Bilinear NLP level** — the true Haverly form (4 bilinears) — used later for NLP module demos;
3. **Convex-QP level (our defensible QP)** — see §4.

---

## 3. FINAL PICK — LP: Williams Refinery Optimization (+ time-expanded scaling)

### 3.1 Formulation sketch (as implemented in the open notebook)

Sets: crudes `c ∈ C`; distillation products `d` (light/medium/heavy naphtha, light/heavy oil, residuum); final products `p ∈ P` (premium fuel, regular fuel, jet fuel, lube oil, fuel oil).

Data (all numeric values in [`Gurobi/modeling-examples/refinery/refinery.ipynb`](https://github.com/Gurobi/modeling-examples/blob/master/refinery/refinery.ipynb)): purchase limit per crude, distillation capacity; splitting coefficients `distillation_splitting_coefficients[d,c]` (crude → fractions); reforming capacity + naphtha→reformed-gasoline yield coefficients; cracking capacity + (light,heavy)-oil→(cracked oil, cracked gasoline) yields; lube-oil factor 0.5 (per unit residuum); octane numbers per component and minima 94 (premium) / 84 (regular); vapor-pressure constants for jet fuel; profit vector per end product; lube-oil min/max production; premium ≥ 0.40·regular.

Decision variables: crude purchases `crudes[c]`, intermediate flows (reform usage, cracking usage, blending allocations), end-product quantities.

```
maximize   Σ_p profit[p] · end_products[p]
subject to Σ_c crudes[c] ≤ distillation_capacity                 (1)
           Σ_d reform_usage[d] ≤ reform_capacity                 (2)
           cracking_usage[light]+cracking_usage[heavy] ≤ crack_capacity   (3)
           Σ_c split[d,c]·crudes[c] == distill_products[d]       (4.1–4.6)
           reform_usage · reform_yields == reformed_gasoline     (4.7)
           cracking_usage · crack_yields == crack_products       (4.8–4.9)
           0.5 · residuum_to_lube == lube_oil                    (4.10)
           blending mass balances for motor fuels / jet / fuel-oil fixed-proportion (5.x)
           Σ_naphthas flows == available naphthas                (mass conservation)
           premium_fuel ≥ 0.40 · regular_fuel                    (7)
           Σ_i octane[i]·blend[i,fuel] ≥ octane_min[fuel]·fuel   (8.1–8.2)
           vapor-pressure bound on jet-fuel blend                (9)
           all vars ≥ 0; lube oil within [lbo_min, lbo_max]
```

Reference solution: optimal profit **$211,365.13** (printed in the notebook).

### 3.2 Objective story (one line for the pitch)

> *"Buy the right crude slate, run the right process units, and blend products that meet octane and vapor-pressure specs — maximize revenue minus crude cost."*

This is verbatim MRPL's monthly planning question (Mangalore refinery crude slate selection across ~dozens of imported + domestic crudes).

### 3.3 Scale ladder for the demo (100s–1000s variables)

Base model ≈ 30–40 vars. Standard textbook time-expansion (multi-period pattern, cf. AMPL book's `steelT` multi-period model, [dev.ampl.com/ampl/books/ampl/examples](https://dev.ampl.com/ampl/books/ampl/examples/index.html)) replicates it over `T` periods with inventory carryover:

- vars ≈ `T × (2 crudes + 11 intermediates + 5 products + inventories)`; constraints grow identically.
- `T=26` weeks → **≈ 1,000 vars / ~800 constraints**; `T=52` → **≈ 2,000+ vars**.
- Crude-slate widening `C: 2→12` crudes (mirroring real slate choice; cf. EMRPS Case 1 using 25 crudes) multiplies further → **5k–10k vars** while every number stays derived from the published base data.
- Narrative support: Alabi & Castro (2008) show real integrated refinery-planning LPs reach **10⁶ constraints × 2.5×10⁶ variables** — exactly where our GPU PDLP module claims its crossover advantage (report §2).
- Bonus robustness demo on the same domain: Netlib `lp/infeas/refinery` (323×464, doctored infeasible petrochemical plant) exercises presolve/IIS diagnostics — a PS-stated requirement ("highly degenerate models, ill-conditioned matrices").

### 3.4 Why MRPL cares

MRPL's business is buying crude (slate optimization under API/sulfur/yield trade-offs), running CDU/reforming/cracking/lube units, and blending BS-VI fuels to tight specs. Every constraint class in this model maps to a named MRPL department function; the shadow prices/duals make a compelling stage moment ("this dual = marginal value of one more barrel of distillation capacity").

---

## 4. FINAL PICK — QP: Haverly pooling family, presented at three levels

### 4.1 Canonical data (Haverly 1978; verified against GAMS Model Library model `haverly`, mirrored at [minlplib.org/haverly.html](https://www.minlplib.org/haverly.html))

Network: feeds A and B → single pool → two products; plus a bypass stream sold directly.

| Item | Value |
|---|---|
| Feed A cost | $6/unit, sulfur content 3 |
| Feed B cost | $16/unit, sulfur content 1 |
| Purchased stream cost | $10/unit, sulfur content 2 |
| Product 1 price | $9/unit, max demand 100, sulfur spec ≤ 2.5 |
| Product 2 price | $15/unit, max demand 200, sulfur spec ≤ 1.5 |
| Pool sulfur `p` | variable in [1,3], set by mix of A/B |

P-formulation (canonical statement):

```
maximize  9·x + 15·y − 6·xA − 16·xB − 10·z            x,y,z ≥ 0
s.t.      p·(y_pool + ȳ_pool) == 3·xA + 1·xB          (pool quality balance — BILINEAR)
          p·y_pool + 2·z_1 ≤ 2.5·(y_pool + z_1)       (product-1 sulfur spec — BILINEAR)
          p·ȳ_pool + 2·z_2 ≤ 1.5·(ȳ_pool + z_2)      (product-2 sulfur spec — BILINEAR)
          mass balances on pool and purchased stream
          y + z_1 ≤ 100 ;  ȳ + z_2 ≤ 200              (demand caps)
```

Instance deltas (verified by diffing [`pooling_haverly{1,2,3}pq.gms`](https://minlplib.org/pooling_haverly1pq.html)):
- **Haverly1**: demand cap 100/200 → optimum profit **400**
- **Haverly2**: caps relaxed to 600/800 → optimum **600**
- **Haverly3**: base caps kept, objective coefficients shifted (price structure change) → optimum **750**

### 4.2 The three-level presentation (defensible QP strategy)

Because pooling is bilinear/NP-hard, we present it as a difficulty ladder on **identical data**:

| Level | What we solve | Solver engine exercised |
|---|---|---|
| L0 | **LP relaxation**: replace bilinear terms by McCormick envelope bounds on `p ∈ [1,3]` (standard convex relaxation; cf. GALINI paper §2) | Simplex / GPU PDLP |
| L1 | **Convex QP (our headline QP)**: same network + linear constraints, objective augmented with quadratic penalty `λ·Σ_j w_j·(quality_j − spec_target_j)²` — variance-minimization blending QP penalizing deviation of blended quality from spec midpoint | ADMM-QP module |
| L2 | **True bilinear NLP** Haverly1/2/3 (4 bilinear terms) | future NLP module; today: show HiGHS/Gurobi agreement vs our L0/L1 bounds |

Justification for L1: (i) it is strictly convex QP — squarely inside the PS-mandated QP scope; (ii) it is industrially meaningful — refiners penalize spec variance because give-away quality (octane giveaway) is lost margin; (iii) it reuses the exact published Haverly numbers, so no invented data; (iv) it gives our ADMM-QP a stage-demo with a visible optimum and KKT residuals, paired against the same instance solved as LP at L0.

### 4.3 Scale ladder (same family, same download location)

| Rung | Instance | Vars / Cons | Optimum |
|---|---|---|---|
| Hand-check | haverly1/2/3 | 10 / 13 | −400 / −600 / −750 |
| Small | bental4 | 13 / 16 | −450 |
| Medium | adhya4 | 58 / 77 | −877.65 |
| Large | foulds3 | **672 / 571** (512 bilinears) | −8 |
| X-large | sppb5 | 8,991 / 8,751 | −60,635.7 |
| XX-large | sppc3 | **10,567 / 10,876** | −122,383.0 |

All rows: MINLPLib instance pages (e.g. [pooling_sppc3pq](https://minlplib.org/instances.html)) with `.lp/.gms/.nl/.py` downloads — the L0 relaxations of the large rungs are ready-made large sparse **LPs** for GPU-PDLP benchmarking, and their McCormick envelopes give large sparse **convex QCQPs** adjacent to our QP scope.

### 4.4 Why MRPL cares

Gasoline/diesel blending is pooling in the wild: MRPL blends dozens of intermediate components into BS-VI grades where sulfur, octane (RON/MON), density and vapor pressure pool nonlinearly. The famous Haverly pathology — LP recursion converging to wrong operating points because pool quality depends on the blend itself (the very reason Haverly wrote the 1978 paper) — is a story every refinery optimizer knows; presenting it shows we understand why commercial solvers (and our sovereign alternative) matter here.

---

## 5. Recommended adoption summary

| Slot | Pick | Primary data link | First action |
|---|---|---|---|
| **LP demo** | Williams Refinery Optimization (base) → time-expanded ×26 weeks (scale) | [refinery.ipynb](https://github.com/Gurobi/modeling-examples/blob/master/refinery/refinery.ipynb) | Encode in own MPS writer + solve with simplex; verify $211,365.13 |
| **QP demo** | Haverly1→3 ladder; L0 LP-relaxation + L1 variance-min QP | [minlplib.org/pooling_haverly1pq.html](https://www.minlplib.org/pooling_haverly1pq.html) | Verify optima 400/600/750 via L0 vs known values |
| Warm-up | DOcplex oil_blending | [docplex-doc](https://ibmdecisionoptimization.github.io/docplex-doc/mp/oil_blending.html) | CLI smoke test |
| Robustness | Netlib lp/infeas `refinery` | [netlib.org/lp/infeas](https://netlib.org/lp/infeas/) | Presolve/IIS demo |
| Stretch | EMRPS benchmark (NLP/MINLP era) | [EMRPS repo](https://github.com/EMRPS/refinery-planning-benchmark) | Post-MVP |

## 6. References

1. Haverly, C.A., *Studies of the Behavior of Recursion for the Pooling Problem*, ACM SIGMAP Bulletin 45:9–28, 1978. doi:10.1145/1111237.1111238
2. Williams, H.P., *Model Building in Mathematical Programming*, 5th ed., Wiley, 2013 — Example #6 "Refinery Optimization", pp. 258, 306–310.
3. Gurobi Modeling Examples, "Refinery". https://github.com/Gurobi/modeling-examples/tree/master/refinery · https://gurobi.github.io/modeling-examples/refinery/
4. MINLPLib — pooling instances + machine-readable metadata. https://minlplib.org/instances.html · https://minlplib.org/instancedata.csv
5. Alfaki, M., Haugland, D., *Strong formulations for the pooling problem*, J. Global Optim. 56:897–916, 2013 (source collection: http://www.ii.uib.no/~mohammeda/spooling/)
6. Ben-Tal, A., Eiger, G., Gershovitz, V., Math. Program. 63:193–212, 1994. Adhya et al., I&EC Res. 38:1956–1972, 1999. Foulds et al., Optimization 24:165–180, 1992.
7. Ceccon, F., Misener, R., *Solving the pooling problem at scale with extensible solver GALINI*, Comput. Chem. Eng. 2022. https://arxiv.org/abs/2105.01687
8. Misener, R., Floudas, C.A., *Advances for the pooling problem: modeling, global optimization, and computational studies*, 2009 (applications survey: petroleum refining, wastewater, supply chain).
9. IBM DOcplex, *Maximizing the profit of an oil company*. https://ibmdecisionoptimization.github.io/docplex-doc/mp/oil_blending.html
10. Chinneck, J.W. (ed.), Netlib lp/infeas collection (`refinery`, T. Baker, 1993). https://netlib.org/lp/infeas/ · https://sparse.tamu.edu/LPnetlib/lpi_refinery
11. Alabi, A., Castro, J., *Dantzig-Wolfe and block coordinate-descent decomposition in large-scale integrated refinery-planning*, UPC DR 2008/01. https://optimization-online.org/2008/09/2080/
12. Du, W. et al., *A production planning benchmark for real-world refinery-petrochemical complexes*. https://arxiv.org/abs/2503.22057 · https://github.com/EMRPS/refinery-planning-benchmark
13. AMPL book example files (multi-period pattern). https://dev.ampl.com/ampl/books/ampl/examples/index.html · https://www.netlib.org/ampl/models/
