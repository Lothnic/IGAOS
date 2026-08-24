#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

ARCH=${CUDA_ARCH:-native}
nvcc -O3 -std=c++17 -arch="$ARCH" -lcublas -lcusparse pdhg.cu -o pdhg_spike

INSTANCES="afiro kb2 adlittle share2b sc50a sc205 bandm grow22"

echo "instance,m,n,nnz,iters,time_ms,obj_p,obj_d,rel_gap,pinf,dinf,status" > spike_results.csv
for i in $INSTANCES; do ./pdhg_spike "data/$i.mps" >> spike_results.csv; done

/usr/bin/python3 highs_baseline.py $(for i in $INSTANCES; do echo "data/$i.mps"; done) > highs_results.csv

/usr/bin/python3 - <<'PYEOF'
import csv

ours = {r["instance"]: r for r in csv.DictReader(open("spike_results.csv"))}
highs = {}
for line in open("highs_results.csv"):
    name, tms, obj, st = line.strip().split(",")
    highs[name] = (float(tms), float(obj), st)

hdr = f"{'instance':<9} {'m':>5} {'n':>5} {'iters':>7} {'pdhg_ms':>9} {'highs_ms':>9} {'ratio':>7} {'rel_gap':>9} {'pinf':>8} {'dinf':>8} {'status':<14} {'obj_match':>10}"
print(hdr)
print("-" * len(hdr))
for name, r in ours.items():
    ratio = float(r["time_ms"]) / highs[name][0] if highs[name][0] > 0 else float("inf")
    match = "n/a"
    if highs[name][1] != 0:
        match = f"{abs(float(r['obj_p']) - highs[name][1]) / abs(highs[name][1]):.2e}"
    print(f"{name:<9} {r['m']:>5} {r['n']:>5} {r['iters']:>7} {float(r['time_ms']):>9.1f} "
          f"{highs[name][0]:>9.1f} {ratio:>6.0f}x {float(r['rel_gap']):>9.2e} "
          f"{float(r['pinf']):>8.1e} {float(r['dinf']):>8.1e} {r['status']:<14} {match:>10}")
PYEOF
