# Python surface (target — lands with pybind11 bindings)

This is the interface the `igaos` Python module will expose once the pybind11
binding target exists. React to this shape; internals churn freely underneath.

```python
import igaos

sol = igaos.solve("model.mps",
                  time_limit=60.0,     # seconds
                  tolerance=1e-4,      # termination ladder rung
                  max_iterations=500_000,
                  presolve=True,
                  seed=0)

sol.status        # "optimal" | "near-optimal" | "feasible" | ... (mirrors Status enum)
sol.objective     # float
sol.x             # list[float], primal solution (len = n)
sol.y             # list[float], duals (len = m)
sol.row_activity  # list[float] (len = m)
sol.pinf          # relative primal infeasibility
sol.dinf          # relative dual infeasibility
sol.rel_gap       # relative duality gap
sol.iterations    # int
sol.solve_time_ms # float
sol.message       # str, engine notes/errors

mdl = igaos.read_mps("model.mps")   # -> Model handle (m, n, nnz, stats) for inspection
```

CLI parity: everything above is also available via the `igaos` executable
(`igaos solve model.mps --time-limit 60 --out sol.json`), same defaults,
same JSON fields.
