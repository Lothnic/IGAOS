# IGAOS

Language for the SIH26119 solver effort: one team, one engine, one vocabulary.

## Language

**Solver core**:
The from-scratch engine that turns a model into a solution (simplex, PDHG, branch-and-bound). Deliberately excludes any modeling environment or GUI.
_Avoid_: framework, platform, stack

**From-scratch boundary**:
The rule that every algorithm inside the solver core is implemented by us from its mathematical foundations; no existing optimization-solver code may be a base. Operationalized in `docs/DEPENDENCIES.md`.
_Avoid_: "built on", wrapper approach

**Primitive**:
A generic numerical linear-algebra routine (SpMV, triangular solves, factorization kernels) with no optimization logic of its own; usable per the dependency policy while the MRPL ruling on cuSPARSE/cuBLAS/cuDSS is pending.
_Avoid_: library (unqualified), helper

**Near-optimal**:
A solution accepted within the published tolerance ladder of the true optimum (solved@1e-4, tight@1e-6 per the benchmark protocol); the honest answer class reported instead of blanket optimality claims.
_Avoid_: good-enough, approximate

**Presolve / postsolve**:
Reductions applied to a model before solving, and their reversal on the returned solution so results refer to the original model.
_Avoid_: preprocessing (ambiguous with ML usage)
