# Industry Comparison Run

IGAOS engines vs Gurobi 13.0.3 (restricted-size evaluation license) vs pinned HiGHS v1.15.1 baseline.
Reference column = HiGHS. `exact`<=1e-6, `ok`<=1e-4, `~near`<=1e-2 relative objective difference.
# Industry Comparison Run

| instance | solver | status | objective | ms | check |
|---|---|---|---|---|---|
| afiro | HiGHS 1.15.1 | Optimal | -464.75314 | 0.4 | exact |
| afiro | Gurobi 13.0.3 | Optimal | -464.75314 | 0.6 | exact |
| afiro | IGAOS simplex | optimal | -464.75314 | 0.5 | exact |
| afiro | IGAOS pdhg | feasible | -464.77152 | 687.4 | ok |
|---|
| sc50a | HiGHS 1.15.1 | Optimal | -64.575077 | 0.6 | exact |
| sc50a | Gurobi 13.0.3 | Optimal | -64.575077 | 0.4 | exact |
| sc50a | IGAOS simplex | optimal | -64.575077 | 3.4 | exact |
| sc50a | IGAOS pdhg | feasible | -64.572191 | 2331.8 | ok |
|---|
| kb2 | HiGHS 1.15.1 | Optimal | -1749.9001 | 0.6 | exact |
| kb2 | Gurobi 13.0.3 | Optimal | -1749.9001 | 0.7 | exact |
| kb2 | IGAOS simplex | optimal | -1749.9001 | 4.7 | exact |
| kb2 | IGAOS pdhg | iteration-limit | -32.443361 | 16617.9 | MISMATCH |
|---|
| adlittle | HiGHS 1.15.1 | Optimal | 225494.96 | 1.0 | exact |
| adlittle | Gurobi 13.0.3 | Optimal | 225494.96 | 0.7 | exact |
| adlittle | IGAOS simplex | optimal | 225494.96 | 13.0 | exact |
| adlittle | IGAOS pdhg | feasible | 225241.5 | 5741.6 | ~near |
|---|
| share2b | HiGHS 1.15.1 | Optimal | -415.73224 | 1.1 | exact |
| share2b | Gurobi 13.0.3 | Optimal | -415.73224 | 0.9 | exact |
| share2b | IGAOS simplex | optimal | -415.73224 | 31.0 | exact |
| share2b | IGAOS pdhg | iteration-limit | -470.64544 | 16709.5 | MISMATCH |
|---|
| sc205 | HiGHS 1.15.1 | Optimal | -52.202061 | 1.4 | exact |
| sc205 | Gurobi 13.0.3 | Optimal | -52.202061 | 0.9 | exact |
| sc205 | IGAOS simplex | optimal | -52.202061 | 266.1 | exact |
| sc205 | IGAOS pdhg | iteration-limit | -3.1196611 | 17489.2 | MISMATCH |
|---|
| bandm | HiGHS 1.15.1 | Optimal | -158.62802 | 5.7 | exact |
| bandm | Gurobi 13.0.3 | Optimal | -158.62802 | 3.8 | exact |
| bandm | IGAOS simplex | optimal | -158.62802 | 6802.3 | exact |
| bandm | IGAOS pdhg | iteration-limit | -145.18672 | 18076.2 | MISMATCH |
|---|
| grow22 | HiGHS 1.15.1 | Optimal | -1.6083434e+08 | 53.8 | exact |
| grow22 | Gurobi 13.0.3 | Optimal | -1.6083434e+08 | 40.8 | exact |
| grow22 | IGAOS simplex | optimal | -1.6083434e+08 | 4356.9 | exact |
| grow22 | IGAOS pdhg | iteration-limit | -80371761 | 18573.7 | MISMATCH |
|---|
| williams_refinery | HiGHS 1.15.1 | Optimal | -211365.13 | 0.4 | exact |
| williams_refinery | Gurobi 13.0.3 | Optimal | -211365.13 | 0.5 | exact |
| williams_refinery | IGAOS simplex | optimal | -211365.13 | 1.3 | exact |
| williams_refinery | IGAOS pdhg | iteration-limit | -98042.371 | 16443.0 | MISMATCH |
|---|
| haverly1_l0 | HiGHS 1.15.1 | Optimal | -500 | 0.4 | exact |
| haverly1_l0 | Gurobi 13.0.3 | Optimal | -500 | 0.2 | exact |
| haverly1_l0 | IGAOS simplex | optimal | -500 | 0.4 | exact |
| haverly1_l0 | IGAOS pdhg | iteration-limit | -393.57355 | 16610.0 | MISMATCH |
|---|
| haverly2_l0 | HiGHS 1.15.1 | Optimal | -2200 | 0.4 | exact |
| haverly2_l0 | Gurobi 13.0.3 | Optimal | -2200 | 0.2 | exact |
| haverly2_l0 | IGAOS simplex | optimal | -2200 | 0.3 | exact |
| haverly2_l0 | IGAOS pdhg | iteration-limit | -209.6902 | 16462.8 | MISMATCH |

## License note

Gurobi runs under its free restricted-size evaluation license (2000 vars / 2000 cons), which covers this whole ladder. A full commercial-license comparison at finale scale is a reasonable ask for the MRPL partnership discussion.

