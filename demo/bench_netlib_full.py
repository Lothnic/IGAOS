#!/usr/bin/env python3
"""Full Netlib mirror run (all ~114 instances, no row cap).

Unlike regress_netlib.py (the 64-instance ≤700-row gate against recorded
HiGHS references), this runs the COMPLETE mirror: HiGHS 1-thread reference
fresh on every instance, then IGAOS simplex, same verdict ladder. Writes
docs/research/netlib_full_results.csv.

Usage: bench_netlib_full.py [instance ...]   (default: all)
"""
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FULL = ROOT / "demo/models/netlib_full"
BIN = Path(__import__("os").environ.get(
    "IGAOS_BIN", ROOT / "build/src/api/igaos"))
OUT = ROOT / "docs/research/netlib_full_results.csv"
TL = int(__import__("os").environ.get("NL_TL", "60"))
HI_TL = int(__import__("os").environ.get("NL_HI_TL", "300"))


def solve_highs(path):
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("time_limit", float(HI_TL))
    h.readModel(str(path))
    t0 = time.perf_counter()
    h.run()
    dt = (time.perf_counter() - t0) * 1000
    lp = h.getLp()
    return (h.modelStatusToString(h.getModelStatus()),
            h.getInfo().objective_function_value, dt,
            (lp.num_row_, lp.num_col_))


def solve_igaos(path):
    try:
        out = subprocess.run(
            [str(BIN), "solve", str(path), "--engine", "simplex",
             "--time-limit", str(TL)],
            capture_output=True, text=True, timeout=TL + 30)
        d = json.loads(out.stdout)
        return d
    except Exception:
        return {"status": "Timeout", "objective": float("nan"),
                "solve_time_ms": (TL + 30) * 1000}


def verdict(obj, ref):
    if obj != obj or ref != ref:
        return "n/a"
    rel = abs(obj - ref) / max(1.0, abs(ref))
    if rel <= 1e-6:
        return "exact"
    if rel <= 1e-4:
        return "ok"
    return "MISMATCH"


def main():
    names = sys.argv[1:] or sorted(
        p.stem for p in FULL.glob("*.mps"))
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines()[1:]:
            if line.strip():
                done.add(line.split(",")[0])
    names = [n for n in names if n not in done]
    print(f"resuming: {len(done)} done, {len(names)} to go")
    fcsv = open(OUT, "a")
    if not done:
        fcsv.write("instance,m,n,our_time_s,highs_time_s,our_objective,"
                   "highs_objective,rel_gap,status,check\n")
    n_ok = n_mism = n_other = n_ref_fail = 0
    for name in names:
        path = FULL / f"{name}.mps"
        try:
            hst, hobj, hms, dims = solve_highs(path)
        except Exception as e:
            print(f"{name}: HiGHS failed {e}")
            n_ref_fail += 1
            continue
        d = solve_igaos(path)
        chk = verdict(d["objective"], hobj)
        if chk in ("exact", "ok"):
            n_ok += 1
        elif chk == "n/a":
            n_other += 1
        else:
            n_mism += 1
        print(f"{name:<14} {str(dims):<14} ref={hobj:<16.8g} "
              f"IGA[{d['status']}]={d['objective']:<16.8g} "
              f"{chk:<9} {d['solve_time_ms']/1000:7.1f}s "
              f"(HiGHS {hms/1000:6.1f}s)")
        w = fcsv.write
        w(f"{name},{dims[0]},{dims[1]},{d['solve_time_ms']/1000:.2f},"
          f"{hms/1000:.2f},{d['objective']:.10g},{hobj:.10g},"
          f"{chk},{d['status']},{chk}\n")
        fcsv.flush()
    fcsv.close()
    n = len(names) - n_ref_fail
    print(f"\n{n_ok}/{n} matched @<=1e-4 ({n_mism} mismatch, "
          f"{n_other} no-verdict, {n_ref_fail} ref-failed)")
    sys.exit(0 if n_mism == 0 else 1)


if __name__ == "__main__":
    main()
