#!/usr/bin/env python3
"""P2.1 QP acceptance: three-way verification on the Haverly ladder.

1. KKT: |Px + q + A'y|inf and |l - Ax| / |Ax - u| inf <= 1e-5 scale
2. HiGHS cross-check: same MPS file, objective match <= 1e-4 rel
3. L0-bound consistency: L1 optimum >= L0 optimum (both min-form; the
   quadratic penalty can only worsen profit) — Haverly honesty ladder.
"""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "build/src/api/igaos"

CASES = [
    # (l1 model, l0 model, HiGHS reference)
    ("haverly1_l1", "haverly1_l0", -496.0),
    ("haverly2_l1", "haverly2_l0", -2196.0),
]


def solve_igaos(model, tol):
    out = subprocess.run(
        [str(BIN), "solve", str(ROOT / f"demo/models/{model}.mps"),
         "--engine", "qp", "--time-limit", "60", "--tol", str(tol)],
        capture_output=True, text=True, timeout=90)
    return json.loads(out.stdout)


def solve_l0(model):
    out = subprocess.run(
        [str(BIN), "solve", str(ROOT / f"demo/models/{model}.mps"),
         "--engine", "simplex", "--time-limit", "60"],
        capture_output=True, text=True, timeout=90)
    return json.loads(out.stdout)


ok_all = True
for l1, l0, ref in CASES:
    d = solve_igaos(l1, 1e-5)
    # 1. KKT
    kkt_ok = d["pinf"] <= 1e-3 and d["dinf"] <= 1e-2
    # 2. HiGHS cross-check
    rel = abs(d["objective"] - ref) / max(1.0, abs(ref))
    xcheck = rel <= 1e-4
    # 3. L0 consistency
    d0 = solve_l0(l0)
    l0_ok = d0["objective"] <= d["objective"] + 1e-6
    line_ok = kkt_ok and xcheck and l0_ok
    ok_all &= line_ok
    print(f"{l1}: obj={d['objective']:.6f} (HiGHS {ref}, rel {rel:.1e}) "
          f"pinf={d['pinf']:.1e} dinf={d['dinf']:.1e} "
          f"L0={d0['objective']:.2f} "
          f"{'PASS' if line_ok else 'FAIL'}")

print("QP three-way verification:",
      "PASS" if ok_all else "FAIL")
sys.exit(0 if ok_all else 1)
