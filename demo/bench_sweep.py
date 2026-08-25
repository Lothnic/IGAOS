#!/usr/bin/env python3
"""Batch sweep over the full Netlib feasible mirror.

Usage: bench_sweep.py <start_idx> <end_idx>   (0-based, into sorted list)
Appends rows to docs/research/sweep_results.csv and prints summary.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import highspy

ROOT = Path(__file__).resolve().parent.parent
FULL = ROOT / "demo/models/netlib_full"
BIN = ROOT / "build/src/api/igaos"
CSV = ROOT / "docs/research/sweep_results.csv"
TIME_LIMIT = int(os.environ.get("SWEEP_TL", "20"))
MAX_ROWS = int(os.environ.get("SWEEP_MAX_ROWS", "600"))


def solve_highs(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    st = h.readModel(str(path))
    lp = h.getLp()
    dims = (lp.num_row_, lp.num_col_)
    t0 = time.perf_counter()
    h.run()
    dt = (time.perf_counter() - t0) * 1000
    return (h.modelStatusToString(h.getModelStatus()),
            h.getInfo().objective_function_value, dt, dims)


def solve_gurobi(path):
    import gurobipy as gp
    gp.setParam("OutputFlag", 0)
    gp.setParam("Threads", 1)
    gp.setParam("TimeLimit", TIME_LIMIT)
    m = gp.read(str(path))
    t0 = time.perf_counter()
    m.optimize()
    dt = (time.perf_counter() - t0) * 1000
    st = {2: "Optimal", 3: "Infeasible", 4: "Unbounded",
          9: "TimeLimit"}.get(m.Status, f"S{m.Status}")
    obj = m.ObjVal if m.SolCount > 0 else float("nan")
    return st, obj, dt


def solve_igaos(path, engine="simplex"):
    try:
        out = subprocess.run(
            [str(BIN), "solve", str(path), "--engine", engine,
             "--time-limit", str(TIME_LIMIT)],
            capture_output=True, text=True, timeout=TIME_LIMIT + 20)
        d = json.loads(out.stdout)
        return d["status"], d["objective"], d["solve_time_ms"]
    except Exception:
        return ("Timeout", float("nan"), float(TIME_LIMIT + 20) * 1000)


def verdict(obj, ref):
    if obj != obj or ref != ref:
        return "n/a"
    rel = abs(obj - ref) / max(1.0, abs(ref))
    if rel <= 1e-6:
        return "exact"
    if rel <= 1e-4:
        return "ok"
    if rel <= 1e-3:
        return "near"
    return "MISMATCH"


def main():
    a, b = int(sys.argv[1]), int(sys.argv[2])
    allf = sorted(FULL.glob("*.mps"))
    files = []
    skipped = []
    import highspy as _hs
    for p in allf[a:b]:
        try:
            hh = _hs.Highs(); hh.setOptionValue('output_flag', False)
            hh.readModel(str(p))
            if hh.getLp().num_row_ <= MAX_ROWS:
                files.append(p)
            else:
                skipped.append((p.stem, hh.getLp().num_row_))
        except Exception:
            pass
    if skipped:
        print(f"skipped {len(skipped)} instances > {MAX_ROWS} rows: "
              + ', '.join(f'{s}[{r}]' for s, r in skipped[:10]))
    print(f"sweeping {len(files)} of {len(allf[a:b])}")
    newfile = not CSV.exists() or CSV.stat().st_size == 0
    fcsv = open(CSV, "a")
    if newfile:
        fcsv.write("instance,m,n,solver,status,objective,ms,check\n")
    for p in files:
        try:
            hst, hobj, hms, dims = solve_highs(p)
        except Exception as e:
            print(f"{p.stem}: highs failed {e}")
            continue
        ist, iobj, ims = solve_igaos(p)
        chk = verdict(iobj, hobj)
        print(f"{p.stem:<14} {str(dims):<12} ref={hobj:<14.8g} "
              f"IGA={iobj:<14.8g} {chk}")
        w = fcsv.write
        w(f"{p.stem},{dims[0]},{dims[1]},HiGHS,{hst},{hobj:.10g},{hms:.1f},ref\n")
        w(f"{p.stem},{dims[0]},{dims[1]},IGAOS-simplex,{ist},{iobj:.10g},{ims:.1f},{chk}\n")
        fcsv.flush()
    fcsv.close()


if __name__ == "__main__":
    main()
