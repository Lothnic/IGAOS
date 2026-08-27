#!/usr/bin/env python3
"""Generate refinery demo LPs as free-format MPS (minimization convention).

Sources (locked in docs/research/refinery-case-studies.md):
  - williams_refinery.mps : Williams, Model Building in Mathematical
    Programming 5th ed., example 6, encoded exactly as implemented in
    Gurobi/modeling-examples/refinery. Published optimum: maximize
    $211365.13 -> min-form optimum -211365.13.
  - haverly{1,2}_l0.mps : Haverly pooling, level-L0 McCormick LP relaxation
    of the canonical data (doc section 4). Known min-form optima -400/-600.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "models")


class MPS:
    def __init__(self, name):
        self.name = name
        self.rows = []
        self.cols = {}
        self.rhs = {}
        self.bounds = []
        self.quad = []  # (colname, colname, value) upper-triangular QUADOBJ

    def row(self, typ, name, rhs=0.0):
        self.rows.append((typ, name))
        if rhs != 0.0:
            self.rhs[name] = rhs

    def _var(self, v):
        return self.cols.setdefault(v, [])

    def obj(self, v, c):
        e = self._var(v)
        for i, (r, cur) in enumerate(e):
            if r == "COST":
                e[i] = ("COST", cur + c)
                return
        e.append(("COST", c))

    def con(self, rowname, v, c):
        assert any(r[1] == rowname for r in self.rows), f"unknown row {rowname}"
        e = self._var(v)
        for i, (r, cur) in enumerate(e):
            if r == rowname:
                e[i] = (rowname, cur + c)
                return
        e.append((rowname, c))

    def bound(self, kind, v, val=None):
        self.bounds.append((kind, v, val))

    def write(self, path):
        lines = [f"NAME {self.name}", "ROWS", " N COST"]
        for typ, nm in self.rows:
            lines.append(f" {typ} {nm}")
        lines.append("COLUMNS")
        for v, entries in self.cols.items():
            emitted = False
            for r, c in entries:
                if c == 0.0:
                    continue
                lines.append(f" {v} {r} {c:.14g}")
                emitted = True
            if not emitted and self.rows:
                lines.append(f" {v} {self.rows[0][1]} 0")
        lines.append("RHS")
        for rw, val in self.rhs.items():
            lines.append(f" RHS1 {rw} {val:.14g}")
        if self.bounds:
            lines.append("BOUNDS")
            for kind, v, val in self.bounds:
                if val is None:
                    lines.append(f" {kind} BND {v}")
                else:
                    lines.append(f" {kind} BND {v} {val:.14g}")
        if self.quad:
            lines.append("QUADOBJ")
            for i, j, v in self.quad:
                lines.append(f" {i} {j} {v:.14g}")
        lines.append("ENDATA")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")


def williams():
    m = MPS("WILLIAMS_REFINERY")

    m.row("L", "DistillCap", 45000)
    m.con("DistillCap", "CR1", 1.0)
    m.con("DistillCap", "CR2", 1.0)
    m.row("L", "ReformCap", 10000)
    m.con("ReformCap", "RU_LN", 1.0)
    m.con("ReformCap", "RU_MN", 1.0)
    m.con("ReformCap", "RU_HN", 1.0)
    m.row("L", "CrackCap", 8000)
    m.con("CrackCap", "CLO", 1.0)
    m.con("CrackCap", "CHO", 1.0)

    splits = {
        "LN": (0.10, 0.15), "MN": (0.20, 0.25), "HN": (0.20, 0.18),
        "LO": (0.12, 0.08), "HO": (0.20, 0.19), "RES": (0.13, 0.12),
    }
    for dpn, (a, b) in splits.items():
        m.row("E", f"Yield_{dpn}")
        m.con(f"Yield_{dpn}", "CR1", a)
        m.con(f"Yield_{dpn}", "CR2", b)
        m.con(f"Yield_{dpn}", dpn, -1.0)

    m.row("E", "Yield_RG")
    for v, c in (("RU_LN", 0.60), ("RU_MN", 0.52), ("RU_HN", 0.45)):
        m.con("Yield_RG", v, c)
    m.con("Yield_RG", "RG", -1.0)

    m.row("E", "Yield_CG")
    m.con("Yield_CG", "CLO", 0.28)
    m.con("Yield_CG", "CHO", 0.20)
    m.con("Yield_CG", "CG", -1.0)
    m.row("E", "Yield_CO")
    m.con("Yield_CO", "CLO", 0.68)
    m.con("Yield_CO", "CHO", 0.75)
    m.con("Yield_CO", "CO", -1.0)

    m.row("E", "Yield_LUBE")
    m.con("Yield_LUBE", "LUBIN", 0.5)
    m.con("Yield_LUBE", "LUBE", -1.0)

    m.row("E", "Sum_PF")
    for v in ("PLN", "PMN", "PHN", "PRG", "PCG"):
        m.con("Sum_PF", v, 1.0)
    m.con("Sum_PF", "PF", -1.0)
    m.row("E", "Sum_RF")
    for v in ("RLN", "RMN", "RHN", "RRG", "RCG"):
        m.con("Sum_RF", v, 1.0)
    m.con("Sum_RF", "RF", -1.0)
    m.row("E", "Sum_JF")
    for v in ("JLO", "JHO", "JRES", "JCO"):
        m.con("Sum_JF", v, 1.0)
    m.con("Sum_JF", "JF", -1.0)

    for nm, terms in (
        ("Mass_LN", [("RU_LN", 1.0), ("RLN", 1.0), ("PLN", 1.0),
                     ("LN", -1.0)]),
        ("Mass_MN", [("RU_MN", 1.0), ("RMN", 1.0), ("PMN", 1.0),
                     ("MN", -1.0)]),
        ("Mass_HN", [("RU_HN", 1.0), ("RHN", 1.0), ("PHN", 1.0),
                     ("HN", -1.0)]),
        ("Mass_LO", [("CLO", 1.0), ("JLO", 1.0), ("FO", 0.55),
                     ("LO", -1.0)]),
        ("Mass_HO", [("CHO", 1.0), ("JHO", 1.0), ("FO", 0.17),
                     ("HO", -1.0)]),
        ("Mass_CO", [("JCO", 1.0), ("FO", 0.22), ("CO", -1.0)]),
        ("Mass_RES", [("LUBIN", 1.0), ("JRES", 1.0), ("FO", 0.055),
                      ("RES", -1.0)]),
        ("Mass_CG", [("RCG", 1.0), ("PCG", 1.0), ("CG", -1.0)]),
        ("Mass_RG", [("RRG", 1.0), ("PRG", 1.0), ("RG", -1.0)]),
    ):
        m.row("E", nm)
        for v, c in terms:
            m.con(nm, v, c)

    m.row("G", "Prem2Reg")
    m.con("Prem2Reg", "PF", 1.0)
    m.con("Prem2Reg", "RF", -0.40)

    octane = (("RLN", 90), ("RMN", 80), ("RHN", 70), ("RRG", 115),
              ("RCG", 105))
    octane_p = (("PLN", 90), ("PMN", 80), ("PHN", 70), ("PRG", 115),
                ("PCG", 105))
    m.row("G", "Octane_RF")
    for v, c in octane:
        m.con("Octane_RF", v, float(c))
    m.con("Octane_RF", "RF", -84.0)
    m.row("G", "Octane_PF")
    for v, c in octane_p:
        m.con("Octane_PF", v, float(c))
    m.con("Octane_PF", "PF", -94.0)

    m.row("L", "VaporPress")
    m.con("VaporPress", "JLO", 1.0)
    m.con("VaporPress", "JHO", 0.60)
    m.con("VaporPress", "JCO", 1.50)
    m.con("VaporPress", "JRES", 0.05)
    m.con("VaporPress", "JF", -1.0)

    profit = {"PF": 7.0, "RF": 6.0, "JF": 4.0, "FO": 3.5, "LUBE": 1.5}
    for v, p in profit.items():
        m.obj(v, -p)

    m.bound("UP", "CR1", 20000)
    m.bound("UP", "CR2", 30000)
    m.bound("LO", "LUBE", 500)
    m.bound("UP", "LUBE", 1000)

    m.write(os.path.join(MODELS, "williams_refinery.mps"))


def haverly(demand1, demand2, tag, quad=None, name_suffix="_l0"):
    m = MPS(f"HAVERLY{tag}{name_suffix.upper()}")

    m.row("E", "PoolMass")
    for v, c in (("Y", 1.0), ("W", 1.0), ("XA", -1.0), ("XB", -1.0)):
        m.con("PoolMass", v, c)

    m.row("E", "QualBal")
    for v, c in (("B1", -1.0), ("XA", 3.0), ("XB", 1.0)):
        m.con("QualBal", v, c)

    m.row("E", "Def_T")
    m.con("Def_T", "T", 1.0)
    m.con("Def_T", "Y", -1.0)
    m.con("Def_T", "W", -1.0)

    m.row("L", "SpecP1")
    m.con("SpecP1", "B2", 1.0)
    m.con("SpecP1", "Z1", -0.5)
    m.con("SpecP1", "Y", -2.5)

    m.row("L", "SpecP2")
    m.con("SpecP2", "B3", 1.0)
    m.con("SpecP2", "Z2", 0.5)
    m.con("SpecP2", "W", -1.5)

    m.row("L", "DemP1", float(demand1))
    m.con("DemP1", "Y", 1.0)
    m.con("DemP1", "Z1", 1.0)
    m.row("L", "DemP2", float(demand2))
    m.con("DemP2", "W", 1.0)
    m.con("DemP2", "Z2", 1.0)

    def mc_cormick(b, q, lo_q, up_q):
        lo_p, up_p = 1.0, 3.0
        m.row("G", f"{b}_mc1",
              rhs=-lo_p * lo_q)
        m.con(f"{b}_mc1", b, 1.0)
        m.con(f"{b}_mc1", q, -lo_p)
        m.con(f"{b}_mc1", "P", -lo_q)
        m.row("G", f"{b}_mc2", rhs=-up_p * up_q)
        m.con(f"{b}_mc2", b, 1.0)
        m.con(f"{b}_mc2", q, -up_p)
        m.con(f"{b}_mc2", "P", -up_q)
        m.row("L", f"{b}_mc3", rhs=-up_p * lo_q)
        m.con(f"{b}_mc3", b, 1.0)
        m.con(f"{b}_mc3", q, -up_p)
        m.con(f"{b}_mc3", "P", -lo_q)
        m.row("L", f"{b}_mc4", rhs=-lo_p * up_q)
        m.con(f"{b}_mc4", b, 1.0)
        m.con(f"{b}_mc4", q, -lo_p)
        m.con(f"{b}_mc4", "P", -up_q)

    mc_cormick("B1", "T", 0.0, float(demand1 + demand2))
    mc_cormick("B2", "Y", 0.0, float(demand1))
    mc_cormick("B3", "W", 0.0, float(demand2))

    m.row("E", "QualLink")
    m.con("QualLink", "B1", 1.0)
    m.con("QualLink", "B2", -1.0)
    m.con("QualLink", "B3", -1.0)

    m.obj("Y", -9.0)
    m.obj("Z1", -9.0)
    m.obj("W", -15.0)
    m.obj("Z2", -15.0)
    m.obj("XA", 6.0)
    m.obj("XB", 16.0)
    m.obj("Z1", 10.0)
    m.obj("Z2", 10.0)

    m.bound("LO", "P", 1.0)
    m.bound("UP", "P", 3.0)
    m.bound("UP", "Y", float(demand1))
    m.bound("UP", "W", float(demand2))

    if quad:
        lam = quad
        m.quad.append(("P", "P", 2.0 * lam))
    m.write(os.path.join(MODELS, f"haverly{tag}{name_suffix}.mps"))


if __name__ == "__main__":
    williams()
    haverly(100, 200, "1")
    haverly(600, 800, "2")
    haverly(100, 200, "1", quad=1.0, name_suffix="_l1")
    haverly(600, 800, "2", quad=1.0, name_suffix="_l1")
    for f in sorted(os.listdir(MODELS)):
        print(f, os.path.getsize(os.path.join(MODELS, f)))
