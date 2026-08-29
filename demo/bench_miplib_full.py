#!/usr/bin/env python3
"""Full MIPLIB 2017 Benchmark Set run (all 240 instances).

Solves every instance with --engine milp, verdicts against the pinned
v36 solufile (=opt=/=best=/=inf=). Resume-capable; writes
docs/research/miplib_full_results.csv.

Usage: bench_miplib_full.py [instance ...]   (default: all 240)
"""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MPLIB = ROOT / "benchmarks/miplib2017"
ALL = MPLIB / "all240"
BIN = Path(__import__("os").environ.get(
    "IGAOS_BIN", ROOT / "build/src/api/igaos"))
SOLU = MPLIB / "miplib2017-v36.solu"
OUT = ROOT / "docs/research/miplib_full_results.csv"
TL = int(__import__("os").environ.get("MIP_TL", "60"))


def load_solu():
    vals = {}
    for line in SOLU.read_text().splitlines():
        p = line.split()
        if len(p) >= 2:
            vals[p[1]] = (p[0], float(p[2]) if len(p) > 2 else None)
    return vals


def solve_igaos(path):
    try:
        out = subprocess.run(
            [str(BIN), "solve", str(path), "--engine", "milp",
             "--time-limit", str(TL)],
            capture_output=True, text=True, timeout=TL + 90)
        return json.loads(out.stdout)
    except Exception:
        return {"status": "Timeout", "objective": float("nan"),
                "solve_time_ms": (TL + 90) * 1000, "message": "harness timeout"}


def main():
    solu = load_solu()
    names = sys.argv[1:] or sorted(
        p.name[:-7] for p in ALL.glob("*.mps.gz"))
    # materialize mps on demand
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines()[1:]:
            if line.strip():
                done.add(line.split(",")[0])
    todo = [n for n in names if n not in done]
    print(f"resuming: {len(done)} done, {len(todo)} to go")
    fcsv = open(OUT, "a")
    if not done:
        fcsv.write("instance,status,our_objective,our_time_s,ref_tag,"
                   "ref_objective,verdict\n")
    n_match = n_inf = n_mism = n_other = 0
    for name in todo:
        gz = ALL / f"{name}.mps.gz"
        mps = MPLIB / f"{name}.mps"
        if not mps.exists():
            subprocess.run(["gunzip", "-kc", str(gz)], stdout=open(mps, "w"),
                           check=True)
        tag, ref = solu.get(name, ("?", None))
        d = solve_igaos(mps)
        obj = d["objective"]
        mps.unlink(missing_ok=True)  # keep disk use flat
        if d["status"] == "infeasible" and tag == "=inf=":
            verdict, n_inf = "inf-proven", n_inf + 1
        elif ref is not None and obj == obj and d["status"] in (
                "optimal", "feasible", "near-optimal"):
            rel = abs(obj - ref) / max(1.0, abs(ref))
            if rel <= 1e-4:
                verdict = "matched"
                n_match += 1
            elif rel <= 1e-2:
                verdict = "near"
                n_other += 1
            else:
                verdict = "MISMATCH"
                n_mism += 1
        else:
            verdict = "no-answer"
            n_other += 1
        print(f"{name:<24} {d['status']:<14} obj={obj:<18.8g} "
              f"{verdict:<11} {d['solve_time_ms']/1000:6.1f}s "
              f"(ref {tag} {ref})")
        fcsv.write(
            f"{name},{d['status']},{obj:.10g},{d['solve_time_ms']/1000:.2f},"
            f"{tag},{ref},{verdict}\n")
        fcsv.flush()
    fcsv.close()
    n = n_match + n_inf + n_mism + n_other
    print(f"\nTOTAL {n}: matched {n_match}, infeasible-proven {n_inf}, "
          f"near {n_other}, mismatch {n_mism}")


if __name__ == "__main__":
    main()
