#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=./build/src/api/igaos
PY=/usr/bin/python3
MODELS="williams_refinery haverly1_l0 haverly2_l0"
NETLIB="afiro sc50a adlittle"

echo "==================================================================="
echo " IGAOS internal-hackathon demo — sovereign optimization solver"
echo "==================================================================="

if [ ! -x "$BIN" ]; then
    echo "[demo] building..."
    cmake -S . -B build >/dev/null
    cmake --build build >/dev/null
fi

if [ ! -f demo/models/williams_refinery.mps ]; then
    echo "[demo] generating refinery models..."
    $PY demo/generate_models.py >/dev/null
fi

for m in $MODELS; do
    f="demo/models/$m.mps"
    echo ""
    echo "--- $m ---------------------------------------------------------"
    "$BIN" info "$f"
    echo ""
    "$BIN" solve "$f" --engine auto --time-limit 20
done

echo ""
echo "--- reference optima (pinned baseline: HiGHS v1.15.1, 1 thread) --"
printf "%-22s %14s %10s\n" instance highs_obj highs_ms
$PY spike/pdhg-spike/highs_baseline.py \
    $(for m in $MODELS; do echo "demo/models/$m.mps"; done) | \
    while IFS=, read -r n t o s; do
        printf "%-22s %14.4f %10.1f\n" "$n" "$o" "$t"
    done

echo ""
echo "--- netlib ladder (GPU PDHG vs HiGHS) ----------------------------"
$PY spike/pdhg-spike/highs_baseline.py \
    $(for m in $NETLIB; do echo "spike/pdhg-spike/data/$m.mps"; done) > /tmp/opencode/demo_highs.csv
printf "%-9s %12s %12s %9s %s\n" instance ours_ms highs_ms ratio status
for m in $NETLIB; do
    out=$("$BIN" solve "spike/pdhg-spike/data/$m.mps" --engine auto --time-limit 20)
    echo "$out" | {
        $PY -c "
import json, sys
d = json.load(sys.stdin)
ref = {}
for line in open('/tmp/opencode/demo_highs.csv'):
    n, t, o, s = line.strip().split(',')
    ref[n] = float(t)
n = d['instance'].split('/')[-1][:-4]
ratio = d['solve_time_ms']/ref[n] if ref[n] > 0 else float('inf')
print(f\"{n:<9} {d['solve_time_ms']:>12.0f} {ref[n]:>12.1f} {ratio:>8.0f}x {d['status']}\")
"
    }
done

echo ""
echo "Honesty note: GPU wins live above ~1M-nnz instances; small LPs are"
echo "CPU territory by design (crossover documented in docs/research/)."
