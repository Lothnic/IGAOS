#!/usr/bin/env python3
"""Robustness suite runner (docs/research/robustness-suite.md).

Usage: bench_robustness.py [--tl T0_TL] [instance-id ...]
Emits the canonical table with appended diagnostic columns and a
per-class PASS/FAIL summary. PASS criteria per robustness-suite.md §3.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
SUITE = ROOT / "benchmarks/suites/robustness.yaml"
BIN = Path(__import__("os").environ.get(
    "IGAOS_BIN", ROOT / "build/src/api/igaos"))
OUT = ROOT / "docs/research/robustness_results.csv"

TIER_TL = {"T0": 60, "T1": 60, "T2": 300, "T3": 300}
# ponytail: global cap keeps interactive runs sane; override for full runs
GLOBAL_TL = int(__import__("os").environ.get("ROB_TL", "300"))

# objective tolerance per class (rel): A/B/D/E 1e-9, C 1e-6 (documented
# cross-solver spread), F n/a
OBJ_TOL = {"A": 1e-9, "B": 1e-9, "C": 1e-6, "D": 1e-9, "E": 1e-9, "F": None}


def load():
    return yaml.safe_load(SUITE.read_text())["instances"]


def igaos_info(path):
    out = subprocess.run([str(BIN), "info", str(path)],
                         capture_output=True, text=True, timeout=60)
    fr = ranges = 0
    for line in out.stdout.splitlines():
        if line.startswith("parsed"):
            parts = line.split()
            fr = int(parts[1].split("=")[1])
            ranges = int(parts[2].split("=")[1])
    return fr, ranges


def igaos_solve(path, tl):
    try:
        out = subprocess.run(
            [str(BIN), "solve", str(path), "--engine", "simplex",
             "--time-limit", str(tl)],
            capture_output=True, text=True, timeout=tl + 30)
        d = json.loads(out.stdout)
        return d
    except Exception:
        return {"status": "Timeout", "objective": float("nan"),
                "solve_time_ms": (tl + 30) * 1000, "iterations": 0,
                "message": "harness timeout"}


def highs_solve(path):
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.readModel(str(path))
    t0 = time.perf_counter()
    h.run()
    dt = (time.perf_counter() - t0) * 1000
    return (h.modelStatusToString(h.getModelStatus()),
            h.getInfo().objective_function_value, dt)


def judge(inst, d, hi, fr, ranges, tl):
    """Return (pass: bool, why: str)."""
    cls, expect = inst["class"], inst["expect"]
    st, obj, ms = d["status"], d["objective"], d["solve_time_ms"]
    if ms / 1000 > tl:
        return False, f"time {ms/1000:.0f}s > TL {tl}s"
    if expect == "infeasible":
        if st == "infeasible":
            return True, "proven infeasible"
        return False, f"expected infeasible, got {st}"
    if st != "optimal":
        return False, f"expected optimal, got {st}"
    zref = inst["z_ref"]
    rel = abs(obj - zref) / max(1.0, abs(zref))
    if rel > OBJ_TOL[cls]:
        return False, f"obj rel err {rel:.2e} > {OBJ_TOL[cls]:.0e}"
    if cls == "D" and inst.get("ranges") is not None:
        if ranges != inst["ranges"]:
            return False, f"RANGES parsed {ranges} != {inst['ranges']}"
    if cls == "E" and inst.get("fr_cols") is not None:
        if fr != inst["fr_cols"]:
            return False, f"FR parsed {fr} != {inst['fr_cols']}"
    return True, "ok"


def main():
    ids = [a for a in sys.argv[1:] if not a.startswith("--")]
    tl_cap = GLOBAL_TL
    for a in sys.argv[1:]:
        if a.startswith("--tl="):
            tl_cap = int(a.split("=")[1])
    insts = load()
    if ids:
        insts = [i for i in insts if i["id"] in ids]
    rows = []
    n_pass = {}
    for inst in insts:
        path = ROOT / "benchmarks" / inst["mps"]
        if not path.exists():
            path = ROOT / "demo/models" / inst["mps"]
        tl = min(TIER_TL[inst["tier"]], tl_cap)
        fr, ranges = igaos_info(path)
        d = igaos_solve(path, tl)
        try:
            hi = highs_solve(path)
        except Exception:
            hi = ("n/a", float("nan"), float("nan"))
        ok, why = judge(inst, d, hi, fr, ranges, tl)
        n_pass.setdefault(inst["class"], [0, 0])
        n_pass[inst["class"]][1] += 1
        if ok:
            n_pass[inst["class"]][0] += 1
        print(f"{inst['id']:<10} {inst['class']} {inst['tier']} "
              f"{'PASS' if ok else 'FAIL':<5} {why:<44} "
              f"obj={d['objective']:<14.8g} {d['solve_time_ms']/1000:7.1f}s "
              f"(HiGHS {hi[2]/1000:6.1f}s)")
        rows.append((inst, d, hi, fr, ranges, ok, why))

    with open(OUT, "w") as f:
        f.write("instance,m,n,our_time_s,highs_time_s,our_objective,"
                "highs_objective,rel_gap,status,class,fr_parsed,"
                "ranges_parsed,verdict\n")
        for inst, d, hi, fr, ranges, ok, why in rows:
            zref = inst.get("z_ref", "")
            rel = ""
            if zref != "" and d["objective"] == d["objective"]:
                rel = abs(d["objective"] - zref) / max(1.0, abs(zref))
                rel = f"{rel:.2e}"
            f.write(f"{inst['id']},,,{d['solve_time_ms']/1000:.2f},"
                    f"{hi[2]/1000:.2f},{d['objective']:.10g},"
                    f"{hi[1]:.10g},{rel},{d['status']},{inst['class']},"
                    f"{fr},{ranges},{'PASS' if ok else 'FAIL'}\n")
    print(f"\nresults -> {OUT.relative_to(ROOT)}")
    print("per class: " + "  ".join(
        f"{c}:{p}/{t}" for c, (p, t) in sorted(n_pass.items())))
    total_p = sum(p for p, t in n_pass.values())
    total_n = sum(t for p, t in n_pass.values())
    print(f"TOTAL: {total_p}/{total_n} PASS")


if __name__ == "__main__":
    main()
