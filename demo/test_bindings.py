#!/usr/bin/env python3
"""Python bindings vs CLI parity: all engines, status strings, read_mps
stats (verified against `igaos info`), and honest pdhg degradation."""
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
import igaos  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BIN = Path(os.environ.get("IGAOS_BIN", ROOT / "build/src/api/igaos"))

# (model, kwargs) — status string, objective, and x length must agree
# with the CLI JSON for every engine.
CASES = [
    ("demo/models/netlib_full/afiro.mps", {"engine": "simplex"}),
    ("demo/models/knap25.mps", {"engine": "milp"}),
    ("demo/models/haverly1_l1.mps", {"engine": "qp"}),
    ("demo/models/williams_refinery.mps", {}),      # auto -> simplex
    ("demo/models/netlib_full/e226.mps", {}),       # obj_const path
]


def cli_solve(path, kw):
    flags = sum(([f"--{k.replace('_', '-')}", str(v)]
                 for k, v in kw.items()), [])
    out = subprocess.run(
        [str(BIN), "solve", str(ROOT / path), "--time-limit", "30"] + flags,
        capture_output=True, text=True, timeout=120)
    return json.loads(out.stdout)


for path, kw in CASES:
    cli = cli_solve(path, kw)
    py = igaos.solve(str(ROOT / path), time_limit=30, **kw)
    # status is a STRING matching the CLI's status strings exactly
    assert isinstance(py.status, str), (path, type(py.status))
    assert py.status == cli["status"], (path, py.status, cli["status"])
    assert isinstance(py.status_enum, igaos.Status), (path, py.status_enum)
    assert abs(py.objective - cli["objective"]) <= 1e-6 * max(
        1, abs(cli["objective"])), (path, py.objective, cli["objective"])
    assert len(py.x) == cli["n"], (path, len(py.x), cli["n"])
    print(f"par: {Path(path).name:<26} engine={kw.get('engine', 'auto'):<8} "
          f"{cli['status']:<10} obj={cli['objective']:.10g} OK")

# --- read_mps stats: williams_refinery vs `igaos info` -----------------
WILLIAMS = "demo/models/williams_refinery.mps"
mdl = igaos.read_mps(str(ROOT / WILLIAMS))
assert (mdl.m, mdl.n, mdl.nnz()) == (29, 36, 106), (mdl.m, mdl.n, mdl.nnz())
info = subprocess.run([str(BIN), "info", str(ROOT / WILLIAMS)],
                      capture_output=True, text=True, timeout=60).stdout
parsed = {}
for line in info.splitlines():
    key, _, rest = line.partition("  ")
    parsed[key] = rest
rows = dict(kv.split("=") for kv in parsed["row types"].split())
cols = dict(kv.split("=") for kv in parsed["col bounds"].split())
pr = dict(kv.split("=") for kv in parsed["parsed"].split())
s = mdl.stats
assert (s["rows_e"], s["rows_l"], s["rows_g"], s["rows_ranged"]) == (
    int(rows["E"]), int(rows["L"]), int(rows["G"]), int(rows["ranged"])), s
assert (s["cols_free"], s["cols_fixed"], s["cols_boxed"],
        s["cols_one_sided"]) == (
    int(cols["free"]), int(cols["fixed"]), int(cols["boxed"]),
    int(cols["one-sided"])), s
assert (s["n_fr_parsed"], s["n_ranges_parsed"]) == (
    int(pr["FR"]), int(pr["RANGES"])), s
assert s["obj_const"] == mdl.obj_const == float(pr["obj_const"]), s
knap = igaos.read_mps(str(ROOT / "demo/models/knap25.mps"))
assert knap.stats["n_int"] == 25, knap.stats
print(f"read_mps stats: williams m={mdl.m} n={mdl.n} nnz={mdl.nnz()} "
      f"E/L/G/ranged={s['rows_e']}/{s['rows_l']}/{s['rows_g']}/"
      f"{s['rows_ranged']} OK; knap25 n_int=25 OK")

# --- pdhg: real engine when built with CUDA, honest error otherwise ----
pd = igaos.solve(str(ROOT / "demo/models/netlib_full/afiro.mps"),
                 engine="pdhg", time_limit=30)
if igaos.has_pdhg:
    assert pd.status != "error", pd.message
    print(f"pdhg: built with CUDA, {pd.status} obj={pd.objective:.6f} OK")
else:
    assert pd.status == "error" and "CUDA" in pd.message, pd.message
    print("pdhg: no CUDA build, honest error OK")

# --- invalid engine: ValueError, not a silent fallback -----------------
try:
    igaos.solve(str(ROOT / "demo/models/knap25.mps"), engine="bogus")
except ValueError as e:
    print(f"invalid engine: ValueError OK ({e})")
else:
    raise AssertionError("engine='bogus' did not raise")

print("bindings parity: PASS")
