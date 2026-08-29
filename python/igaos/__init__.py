"""IGAOS — sovereign LP/MILP/QP solver core.

A from-scratch optimization suite: revised primal+dual simplex with
sparse LU basis factorization, Devex pricing, Andersen-Andersen presolve,
branch-and-bound with Gomory mixed-integer cuts, OSQP-style ADMM for
convex QP, and (when built with CUDA) a GPU first-order PDHG engine.

Quick start:
    >>> import igaos
    >>> sol = igaos.solve("model.mps", engine="auto", time_limit=60)
    >>> sol.status, sol.objective
"""

from ._core import *  # noqa: F401,F403
from ._core import (  # noqa: F401  (explicit re-exports for IDEs)
    Model,
    Status,
    has_pdhg,
    read_mps,
    solve,
)

__version__ = "0.1.0"
