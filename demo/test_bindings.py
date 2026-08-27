#!/usr/bin/env python3
"""P1.2 acceptance: python bindings vs CLI JSON parity on afiro + williams."""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
import igaos  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CASES = [
    ("demo/models/williams_refinery.mps", {}),
    ("demo/models/netlib_full/afiro.mps", {}),
    ("demo/models/netlib_full/e226.mps", {}),  # obj_const path
    ("demo/models/knap25.mps", {"engine": "milp"}),
]

for path, kw in CASES:
    cli = json.loads(subprocess.run(
        [str(ROOT / "build/src/api/igaos"), "solve", str(ROOT / path),
         "--time-limit", "30"] +
        sum(([f"--{k.replace('_', '-')}", str(v)]
             for k, v in kw.items()), []),
        capture_output=True, text=True, timeout=60).stdout)
    py = igaos.solve(str(ROOT / path), time_limit=30, **kw)
    assert str(py.status).split(".")[-1].lower() == cli["status"], \
        (path, py.status, cli["status"])
    assert abs(py.objective - cli["objective"]) <= 1e-6 * max(
        1, abs(cli["objective"])), (path, py.objective, cli["objective"])
    assert len(py.x) == cli["n"], (path, len(py.x), cli["n"])
    print(f"par: {Path(path).name:<24} {cli['status']:<10} "
          f"obj={cli['objective']:.10g} OK")

print("bindings parity: PASS")
