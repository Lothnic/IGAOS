#!/usr/bin/env python3
"""Industry comparison bench: IGAOS engines vs Gurobi vs HiGHS.

Run from repo root:  /usr/bin/python3 demo/bench_industry.py [instance.mps ...]

Reference objective = pinned HiGHS v1.15.1 (benchmark protocol).
Validity marks: 'exact' within 1e-6 rel, 'ok' within 1e-4,
'~near' within 1e-2, else MISMATCH.
"""

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "build/src/api/igaos"

DEFAULT_SET = [
    "spike/pdhg-spike/data/afiro.mps",
    "spike/pdhg-spike/data/sc50a.mps",
    "spike/pdhg-spike/data/kb2.mps",
    "spike/pdhg-spike/data/adlittle.mps",
    "spike/pdhg-spike/data/share2b.mps",
    "spike/pdhg-spike/data/sc205.mps",
    "spike/pdhg-spike/data/bandm.mps",
    "spike/pdhg-spike/data/grow22.mps",
    "demo/models/williams_refinery.mps",
    "demo/models/haverly1_l0.mps",
    "demo/models/haverly2_l0.mps",
]


def solve_highs(path):
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


def solve_gurobi(path):
    import gurobipy as gp
    gp.setParam("OutputFlag", 0)
    gp.setParam("Threads", 1)
    gp.setParam("TimeLimit", 60)
    m = gp.read(str(path))
    t0 = time.perf_counter()
    m.optimize()
    dt = (time.perf_counter() - t0) * 1000
    st = {2: "Optimal", 3: "Infeasible", 4: "Unbounded",
          9: "TimeLimit", 11: "Interrupt"}.get(m.Status, f"Status{m.Status}")
    obj = m.ObjVal if m.SolCount > 0 else float("nan")
    return st, obj, dt


def solve_igaos(path, engine):
    out = subprocess.run(
        [str(BIN), "solve", str(path), "--engine", engine, "--time-limit", "60"],
        capture_output=True, text=True, timeout=180)
    import json
    try:
        d = json.loads(out.stdout)
        return d["status"], d["objective"], d["solve_time_ms"]
    except Exception:
        return ("Error", float("nan"), float("nan"))


def verdict(obj, ref):
    if obj != obj or ref != ref:
        return "n/a"
    rel = abs(obj - ref) / max(1.0, abs(ref))
    if rel <= 1e-6:
        return "exact"
    if rel <= 1e-4:
        return "ok"
    if rel <= 1e-2:
        return "~near"
    return "MISMATCH"


def main():
    paths = [Path(p) for p in (sys.argv[1:] or DEFAULT_SET)]
    print(f"{'instance':<20} {'solver':<14} {'status':<16} "
          f"{'objective':>16}  {'ms':>9}  check")
    print("-" * 92)

    md = ["# Industry Comparison Run", "",
          "| instance | solver | status | objective | ms | check |",
          "|---|---|---|---|---|---|"]

    for p in paths:
        rows = []
        rows.append(("HiGHS 1.15.1", solve_highs(p)))
        try:
            rows.append(("Gurobi 13.0.3", solve_gurobi(p)))
        except Exception as e:
            rows.append(("Gurobi 13.0.3", ("license/error", float("nan"), 0)))
        rows.append(("IGAOS simplex", solve_igaos(p, "simplex")))
        rows.append(("IGAOS pdhg", solve_igaos(p, "pdhg")))

        ref = rows[0][1][1]
        for name, (st, obj, ms) in rows:
            chk = verdict(obj, ref)
            print(f"{p.stem:<20} {name:<14} {st:<16} {obj:>16.8g}  "
                  f"{ms:>9.1f}  {chk}")
            md.append(f"| {p.stem} | {name} | {st} | {obj:.8g} | {ms:.1f} | "
                      f"{chk} |")
        print("-" * 92)
        md.append("|---|")

    out = ROOT / "docs/research/industry-comparison.md"
    header = ("# Industry Comparison Run\n\n"
              "IGAOS engines vs Gurobi 13.0.3 (restricted-size evaluation "
              "license) vs pinned HiGHS v1.15.1 baseline.\n"
              "Reference column = HiGHS. `exact`<=1e-6, `ok`<=1e-4, "
              "`~near`<=1e-2 relative objective difference.\n")
    body = "\n".join(md[:-1])
    note = ("\n\n## License note\n\n"
            "Gurobi runs under its free restricted-size evaluation license "
            "(2000 vars / 2000 cons), which covers this whole ladder. A "
            "full commercial-license comparison at finale scale is a "
            "reasonable ask for the MRPL partnership discussion.\n")
    out.write_text(header + body + note + "\n")
    print(f"\nsaved: {out}")


if __name__ == "__main__":
    main()
