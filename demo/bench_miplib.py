#!/usr/bin/env python3
"""MILP validation vs the official MIPLIB 2017 solufile (v36, pinned).

Usage: bench_miplib.py [instance ...]   (default: the 8 smallest core-16)
Compares IGAOS B&B incumbents against =opt=/=best= values in the solufile.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MPLIB = ROOT / "benchmarks/miplib2017"
BIN = ROOT / "build/src/api/igaos"
SOLU = MPLIB / "miplib2017-v36.solu"
TIME_LIMIT = int(__import__("os").environ.get("MIP_TL", "60"))

CORE8 = ["flugpl", "gt2", "khb05250", "air03", "p0201",
         "mod010", "ran13x13", "misc07"]


def load_solu():
    vals = {}
    for line in SOLU.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            vals[parts[1]] = (parts[0], float(parts[2]) if len(parts) > 2 else None)
    return vals


def solve_igaos(path):
    try:
        out = subprocess.run(
            [str(BIN), "solve", str(path), "--engine", "milp",
             "--time-limit", str(TIME_LIMIT)],
            capture_output=True, text=True, timeout=TIME_LIMIT + 30)
        d = json.loads(out.stdout)
        return d["status"], d["objective"], d["solve_time_ms"], d.get("message", "")
    except Exception:
        return ("Timeout", float("nan"), (TIME_LIMIT + 30) * 1000, "")


def verdict(obj, ref):
    if obj != obj or ref is None or ref != ref:
        return "n/a"
    rel = abs(obj - ref) / max(1.0, abs(ref))
    if rel <= 1e-6:
        return "exact"
    if rel <= 1e-4:
        return "solved"
    return "MISMATCH"


def main():
    names = sys.argv[1:] or CORE8
    solu = load_solu()
    n_ok = 0
    for name in names:
        path = MPLIB / f"{name}.mps"
        tag, ref = solu.get(name, ("?", None))
        st, obj, ms, msg = solve_igaos(path)
        chk = verdict(obj, ref)
        if chk in ("exact", "solved"):
            n_ok += 1
        print(f"{name:<12} ref[{tag}]={str(ref):<18} IGA[{st}]={obj:<18.10g} "
              f"{ms / 1000:7.1f}s  {chk}" + (f"  ({msg})" if msg else ""))
    print(f"\n{n_ok}/{len(names)} matched solufile @<=1e-4")


if __name__ == "__main__":
    main()
