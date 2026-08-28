#!/usr/bin/env bash
# IGAOS finale demo — 2-minute tour (issue #22 video plan, issue #1 finale sequence)
# Sections: Williams duals -> Haverly honesty ladder -> infeasible refinery proof
#           -> MIPLIB MILP -> GPU crossover showcase.
# Every number is computed live or read from docs/research/gpu_showcase.md (recorded runs).
set -euo pipefail
cd "$(dirname "$0")/.."

if [ -x ./build-fin/src/api/igaos ]; then
    BIN=./build-fin/src/api/igaos
else
    BIN=./build/src/api/igaos
fi
PY=/usr/bin/python3
BIG=benchmarks/lpfeas/ex10.mps   # 1.16M-nnz Mittelmann LPfeas instance (gpu_showcase.md)
SHOWCASE_DOC=docs/research/gpu_showcase.md

banner() { echo ""; echo "--- $1 ---------------------------------------------------------"; }

echo "==================================================================="
echo " IGAOS — sovereign optimization solver · finale demo"
echo "==================================================================="

if [ ! -x "$BIN" ]; then
    echo "[demo] building..."
    cmake -S . -B build-fin >/dev/null
    cmake --build build-fin -j >/dev/null
fi
if [ ! -f demo/models/williams_refinery.mps ]; then
    $PY demo/generate_models.py >/dev/null
fi

# ---------------------------------------------------------------------------
banner "1 · Williams refinery LP — optimal plan + the duals that price it"
$BIN info demo/models/williams_refinery.mps
if $PY -c "import sys; sys.path.insert(0,'python'); import igaos" 2>/dev/null; then
$PY - <<'EOF'
import sys; sys.path.insert(0, "python")
import igaos
r = igaos.solve("demo/models/williams_refinery.mps", engine="simplex", time_limit=10)
names = []
with open("demo/models/williams_refinery.mps") as f:
    in_rows = False
    for line in f:
        if line.startswith("ROWS"): in_rows = True; continue
        if in_rows and line[0] not in " \t": break
        if in_rows:
            p = line.split()
            if len(p) == 2 and p[0] != "N": names.append(p[1])
print(f"status      : {r.status}")
print(f"profit      : $ {-r.objective:,.2f}   (reference: $211,365.13, Gurobi notebook)")
print(f"simplex iters: {r.iterations}   time: {r.solve_time_ms:.1f} ms")
print("top duals   :")
for name, y in sorted(zip(names, r.y), key=lambda t: -abs(t[1]))[:4]:
    print(f"  {name:<12s} {y:+10.2f}   ($/unit marginal value of this constraint)")
EOF
else
    # ponytail: pybind11 module not built — CLI-only, skip the duals display
    $BIN solve demo/models/williams_refinery.mps --engine simplex --time-limit 10 | \
        $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"status: {d['status']}  profit: $ {-d['objective']:,.2f}  ({d['iterations']} iters, {d['solve_time_ms']:.1f} ms)\"); print('[duals display needs the python module: build with pybind11 present]')"
fi
echo "TAKEAWAY: refinery-wide optimum in 1 simplex pass; duals price each unit's capacity."

# ---------------------------------------------------------------------------
banner "2 · Haverly pooling honesty ladder — LP -500 vs QP -496"
$BIN solve demo/models/haverly1_l0.mps --engine simplex --time-limit 10 | \
    $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"L0 LP relax : obj {d['objective']:.1f}  ({d['iterations']} iters, {d['solve_time_ms']:.1f} ms, {d['engine']})\")"
$BIN solve demo/models/haverly1_l1.mps --engine qp --time-limit 10 | \
    $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"L1 QP (bilinear penalty): obj {d['objective']:.1f}  ({d['iterations']} iters, {d['solve_time_ms']:.1f} ms, {d['engine']})\")"
echo "  (min form; literature optimum of the true bilinear NLP: -400 — pooling is NP-hard,"
echo "   we present it honestly at the levels our engines actually solve)"
echo ""
echo "2b · infeasibility proof — Netlib lp/infeas/refinery.mps (doctored petrochemical plant)"
$BIN solve benchmarks/netlib/infeas/refinery.mps --engine simplex --time-limit 20 | \
    $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"status: {d['status']}  ({d['message']})\")"
echo "TAKEAWAY: three honesty levels on one pooling family, plus a certified infeasible verdict."

# ---------------------------------------------------------------------------
banner "3 · Real MILP — MIPLIB 2017 khb05250, branch-and-bound to proven optimality"
$BIN solve benchmarks/miplib2017/khb05250.mps --engine milp --time-limit 60 | \
    $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"status: {d['status']}  obj {d['objective']:.0f}  ({d['message']}, {d['solve_time_ms']/1000:.1f} s)\")"
echo "TAKEAWAY: integer optima with a proof certificate, not just a heuristic answer."

# ---------------------------------------------------------------------------
banner "4 · GPU crossover — PDHG vs simplex/HiGHS"
echo "afiro (Netlib, 27 rows — CPU territory, engine contrast):"
$BIN solve spike/pdhg-spike/data/afiro.mps --engine simplex --time-limit 20 | \
    $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"  simplex : {d['status']}  obj {d['objective']:.4f}  {d['solve_time_ms']:.1f} ms\")"
$BIN solve spike/pdhg-spike/data/afiro.mps --engine pdhg --time-limit 20 | \
    $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"  PDHG/GPU: {d['status']}  obj {d['objective']:.4f}  {d['solve_time_ms']:.1f} ms  pinf {d['pinf']:.1e}\")"
echo ""
echo "ex10 (Mittelmann LPfeas, 1,162,000 nnz — GPU territory):"
if [ -f "$BIG" ]; then
    $BIN solve "$BIG" --engine pdhg --time-limit 60 --tol 1e-4 | \
        $PY -c "import json,sys; d=json.load(sys.stdin); print(f\"  PDHG/GPU: {d['status']}  obj {d['objective']:.4f}  pinf {d['pinf']:.1e}  gap {d['rel_gap']:.1e}  {d['solve_time_ms']/1000:.1f} s  ({d['iterations']} iters)\")"
    echo "  HiGHS 1-thread reference (recorded run): 76.2 s to exact optimum (obj 100.0)"
else
    echo "  [instance not downloaded — recorded numbers from $SHOWCASE_DOC]"
    echo "  PDHG/GPU: near-optimal  obj 100.0000  pinf 9.2e-05  gap 4.6e-05  4.1 s  (14150 iters)"
    echo "  HiGHS 1-thread reference (recorded run): 76.2 s to exact optimum (obj 100.0)"
fi
echo "TAKEAWAY: below ~1M nnz the CPU simplex wins (afiro: sub-ms vs ~0.1-0.2 s); at 1.16M nnz"
echo "          GPU PDHG is ~19x faster than HiGHS — the crossover is measured, not claimed."
echo ""
echo "Full table + honest caveats: $SHOWCASE_DOC"
