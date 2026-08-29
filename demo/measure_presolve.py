#!/usr/bin/env python3
"""Presolve reduction + timing sweep: netlib_full (64) + 3 MIPLIB instances.

Runs each instance with the full presolve (on) and the legacy bound-only
path (off), collecting the reduction report line from stderr and
solve_time_ms from stdout JSON.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "build-pre/src/api/igaos"
NETLIB = ROOT / "demo/models/netlib_full"
MIPLIB = ROOT / "benchmarks/miplib2017"
TL = 30

names = sorted(p.stem for p in NETLIB.glob("*.mps"))
mip = ["khb05250", "p0201", "mod010"]

REP = re.compile(
    r"presolve: (\d+)x(\d+)/(\d+) -> (\d+)x(\d+)/(\d+) "
    r"\(rows: (\d+) red (\d+) forced (\d+) dup; cols: (\d+) fix (\d+) empty "
    r"(\d+) subst (\d+) dup; (\d+) tightens\)")

rows = []
for name in names + mip:
    d = NETLIB if name not in mip else MIPLIB
    path = d / f"{name}.mps"
    # presolve on
    env = dict(os.environ, IGAOS_PRESOLVE_REPORT="1")
    out = subprocess.run(
        [str(BIN), "solve", str(path), "--engine", "simplex",
         "--time-limit", str(TL)],
        capture_output=True, text=True, timeout=TL + 15, env=env)
    try:
        j = json.loads(out.stdout)
        st, obj, t_on = j["status"], j["objective"], j["solve_time_ms"]
    except Exception:
        st, obj, t_on = "crash", float("nan"), float("nan")
    rep = REP.search(out.stderr)
    # presolve off (legacy bound-only path)
    out2 = subprocess.run(
        [str(BIN), "solve", str(path), "--engine", "simplex",
         "--time-limit", str(TL), "--presolve", "off"],
        capture_output=True, text=True, timeout=TL + 15)
    try:
        j2 = json.loads(out2.stdout)
        t_off, st2 = j2["solve_time_ms"], j2["status"]
    except Exception:
        t_off, st2 = float("nan"), "crash"
    if rep:
        g = [int(x) for x in rep.groups()]
        rows.append((name, st, obj, t_on, t_off, st2, g))
    else:
        rows.append((name, st, obj, t_on, t_off, st2, None))

print(f"{'instance':<12} {'status':<12} {'m x n / nnz':>22} "
      f"{'-> reduced':>18} {'t_off':>8} {'t_on':>8}")
tot = {"m0": 0, "m1": 0, "n0": 0, "n1": 0, "nnz0": 0, "nnz1": 0}
passes = [0] * 8
t_off_sum = t_on_sum = 0.0
for name, st, obj, t_on, t_off, st2, g in rows:
    if g is None:
        print(f"{name:<12} {st:<12} (no report line)")
        continue
    m0, n0, nnz0, m1, n1, nnz1, red, forced, rdup, fix, empty, subst, cdup, tight = g
    tot["m0"] += m0; tot["m1"] += m1
    tot["n0"] += n0; tot["n1"] += n1
    tot["nnz0"] += nnz0; tot["nnz1"] += nnz1
    for i, v in enumerate((red, forced, rdup, fix, empty, subst, cdup, tight)):
        passes[i] += v
    t_off_sum += t_off if t_off == t_off else 0
    t_on_sum += t_on if t_on == t_on else 0
    print(f"{name:<12} {st:<12} {m0:>7}x{n0:<6} {nnz0:>7} "
          f"{m1:>6}x{n1:<6} {nnz1:>6} {t_off:>7.1f} {t_on:>7.1f}")

print(f"\nTOTAL rows {tot['m0']} -> {tot['m1']} "
      f"({100.0*(1-tot['m1']/tot['m0']):.1f}% cut), "
      f"cols {tot['n0']} -> {tot['n1']} "
      f"({100.0*(1-tot['n1']/tot['n0']):.1f}% cut), "
      f"nnz {tot['nnz0']} -> {tot['nnz1']} "
      f"({100.0*(1-tot['nnz1']/tot['nnz0']):.1f}% cut)")
print("pass totals: rows_redundant=%d rows_forcing=%d rows_dup=%d "
      "cols_fixed=%d cols_empty=%d cols_subst=%d cols_dup=%d tightens=%d"
      % tuple(passes))
print(f"time sum: off={t_off_sum:.0f}ms on={t_on_sum:.0f}ms")
