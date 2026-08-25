# Industry Comparison Run

IGAOS engines vs Gurobi 13.0.3 (restricted-size evaluation license) vs pinned HiGHS v1.15.1 baseline.
Reference column = HiGHS. `exact`<=1e-6, `ok`<=1e-4, `~near`<=1e-2 relative objective difference.
# Industry Comparison Run

| instance | solver | status | objective | ms | check |
|---|---|---|---|---|---|
| afiro | HiGHS 1.15.1 | Optimal | -464.75314 | 3.4 | exact |
| afiro | Gurobi 13.0.3 | Optimal | -464.75314 | 0.4 | exact |
| afiro | IGAOS simplex | optimal | -464.75314 | 0.5 | exact |
| afiro | IGAOS pdhg | feasible | -464.77152 | 678.2 | ok |
|---|
| sc50a | HiGHS 1.15.1 | Optimal | -64.575077 | 0.5 | exact |
| sc50a | Gurobi 13.0.3 | Optimal | -64.575077 | 0.4 | exact |
| sc50a | IGAOS simplex | optimal | -64.575077 | 3.5 | exact |
| sc50a | IGAOS pdhg | feasible | -64.572191 | 2292.1 | ok |
|---|
| kb2 | HiGHS 1.15.1 | Optimal | -1749.9001 | 0.6 | exact |
| kb2 | Gurobi 13.0.3 | Optimal | -1749.9001 | 0.7 | exact |
| kb2 | IGAOS simplex | optimal | -1749.9001 | 6.9 | exact |
| kb2 | IGAOS pdhg | iteration-limit | -32.443361 | 16665.3 | MISMATCH |
|---|
| adlittle | HiGHS 1.15.1 | Optimal | 225494.96 | 1.0 | exact |
| adlittle | Gurobi 13.0.3 | Optimal | 225494.96 | 0.8 | exact |
| adlittle | IGAOS simplex | optimal | 225494.96 | 13.1 | exact |
| adlittle | IGAOS pdhg | feasible | 225241.5 | 6105.6 | ~near |
|---|
| share2b | HiGHS 1.15.1 | Optimal | -415.73224 | 1.3 | exact |
| share2b | Gurobi 13.0.3 | Optimal | -415.73224 | 1.2 | exact |
| share2b | IGAOS simplex | optimal | -415.73224 | 45.4 | exact |
| share2b | IGAOS pdhg | iteration-limit | -470.64544 | 16348.5 | MISMATCH |
|---|
| sc205 | HiGHS 1.15.1 | Optimal | -52.202061 | 1.5 | exact |
| sc205 | Gurobi 13.0.3 | Optimal | -52.202061 | 1.0 | exact |
| sc205 | IGAOS simplex | optimal | -52.202061 | 304.2 | exact |
| sc205 | IGAOS pdhg | iteration-limit | -3.1196611 | 17376.6 | MISMATCH |
|---|
| bandm | HiGHS 1.15.1 | Optimal | -158.62802 | 6.3 | exact |
| bandm | Gurobi 13.0.3 | Optimal | -158.62802 | 4.1 | exact |
| bandm | IGAOS simplex | optimal | -158.62802 | 9666.0 | exact |
| bandm | IGAOS pdhg | iteration-limit | -145.18672 | 18119.2 | MISMATCH |
|---|
| grow22 | HiGHS 1.15.1 | Optimal | -1.6083434e+08 | 58.4 | exact |
| grow22 | Gurobi 13.0.3 | Optimal | -1.6083434e+08 | 41.7 | exact |
| grow22 | IGAOS simplex | optimal | -1.965e+08 | 3792.1 | MISMATCH |
| grow22 | IGAOS pdhg | iteration-limit | -80371761 | 18903.5 | MISMATCH |
|---|
| williams_refinery | HiGHS 1.15.1 | Optimal | -211365.13 | 0.5 | exact |
| williams_refinery | Gurobi 13.0.3 | Optimal | -211365.13 | 0.6 | exact |
| williams_refinery | IGAOS simplex | optimal | -215338.29 | 1.6 | MISMATCH |
| williams_refinery | IGAOS pdhg | iteration-limit | -98042.371 | 17161.1 | MISMATCH |
|---|
| haverly1_l0 | HiGHS 1.15.1 | Optimal | -500 | 0.4 | exact |
| haverly1_l0 | Gurobi 13.0.3 | Optimal | -500 | 0.2 | exact |
| haverly1_l0 | IGAOS simplex | optimal | 1082.4219 | 0.4 | MISMATCH |
| haverly1_l0 | IGAOS pdhg | iteration-limit | -393.57355 | 16931.3 | MISMATCH |
|---|
| haverly2_l0 | HiGHS 1.15.1 | Optimal | -2200 | 0.5 | exact |
| haverly2_l0 | Gurobi 13.0.3 | Optimal | -2200 | 0.2 | exact |
| haverly2_l0 | IGAOS simplex | optimal | 9567.1875 | 0.3 | MISMATCH |
| haverly2_l0 | IGAOS pdhg | iteration-limit | -209.6902 | 17247.7 | MISMATCH |

## License note

Gurobi runs under its free restricted-size evaluation license (2000 vars / 2000 cons), which covers this whole ladder. A full commercial-license comparison at finale scale is a reasonable ask for the MRPL partnership discussion.

