# Dependency Policy

The problem statement (SIH26119) forbids building upon any existing open-source solver
library. This page records where IGAOS draws the line.

## Allowed

| Class | Examples | Status |
|---|---|---|
| BLAS-class numerical primitives | cuSPARSE, cuBLAS, cuDSS; dense header linear algebra of the same class | Pending written MRPL ruling — tracked in [issue #2](https://github.com/Lothnic/IGAOS/issues/2); do not hard-depend until it lands |
| Dev tools (never ship in the solver) | pybind11, Catch2/doctest, fmt, CMake | Free |

## Forbidden

Anything containing optimization strategy or solver algorithms: OSQP, SCS, OR-Tools,
HiGHS, CBC, GLPK, SCIP, Gurobi/Xpress/CPLEX bindings — as linked libraries, vendored
code, or translated source. Published *algorithms* (papers) are fine to implement from
their mathematical foundations; published *code* is not a base.

## Rule of thumb

If the artifact decides how an optimization is solved, we write it.
If it only computes standard numerical operations we call, it is a primitive.
