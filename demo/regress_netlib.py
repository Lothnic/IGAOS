#!/usr/bin/env python3
"""Netlib regression sweep vs recorded HiGHS references.

Usage: regress_netlib.py            (all 64)
       regress_netlib.py perold e226 ...
Reads docs/research/sweep_results.csv for HiGHS reference objectives and
re-solves every instance with the current simplex. No CSV writes — a
read-only gate for changes to the solver core.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FULL = ROOT / "demo/models/netlib_full"
BIN = Path(__import__("os").environ.get(
    "IGAOS_BIN", ROOT / "build/src/api/igaos"))
REF = ROOT / "docs/research/sweep_results.csv"
TL = 30  # generous uniform cap; earlier runs used 20

ref = {}
for r in csv.DictReader(open(REF)):
    if r["solver"] == "HiGHS":
        ref[r["instance"]] = float(r["objective"])

names = sys.argv[1:] or sorted(ref)
n = ok = mism = other = 0
for name in names:
    n += 1
    try:
        out = subprocess.run(
            [str(BIN), "solve", str(FULL / f"{name}.mps"),
             "--engine", "simplex", "--time-limit", str(TL)],
            capture_output=True, text=True, timeout=TL + 15)
        d = json.loads(out.stdout)
        st, obj = d["status"], d["objective"]
    except Exception:
        st, obj = "Timeout", float("nan")
    rv = ref[name]
    rel = abs(obj - rv) / max(1.0, abs(rv)) if obj == obj else float("inf")
    if st == "optimal" and rel <= 1e-6:
        ok += 1
    elif st == "optimal":
        mism += 1
        print(f"  MISMATCH {name}: {obj} vs {rv}")
    else:
        other += 1
        print(f"  {st:<12} {name}")
print(f"\n{ok}/{n} exact (recorded best 49), {mism} optimal-mismatch, "
      f"{other} non-optimal")
sys.exit(0 if ok >= 49 else 1)
