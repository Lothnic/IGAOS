NAME HAVERLY2_L1
ROWS
 N COST
 E PoolMass
 E QualBal
 E Def_T
 L SpecP1
 L SpecP2
 L DemP1
 L DemP2
 G B1_mc1
 G B1_mc2
 L B1_mc3
 L B1_mc4
 G B2_mc1
 G B2_mc2
 L B2_mc3
 L B2_mc4
 G B3_mc1
 G B3_mc2
 L B3_mc3
 L B3_mc4
 E QualLink
COLUMNS
 Y PoolMass 1
 Y Def_T -1
 Y SpecP1 -2.5
 Y DemP1 1
 Y B2_mc1 -1
 Y B2_mc2 -3
 Y B2_mc3 -3
 Y B2_mc4 -1
 Y COST -9
 W PoolMass 1
 W Def_T -1
 W SpecP2 -1.5
 W DemP2 1
 W B3_mc1 -1
 W B3_mc2 -3
 W B3_mc3 -3
 W B3_mc4 -1
 W COST -15
 XA PoolMass -1
 XA QualBal 3
 XA COST 6
 XB PoolMass -1
 XB QualBal 1
 XB COST 16
 B1 QualBal -1
 B1 B1_mc1 1
 B1 B1_mc2 1
 B1 B1_mc3 1
 B1 B1_mc4 1
 B1 QualLink 1
 T Def_T 1
 T B1_mc1 -1
 T B1_mc2 -3
 T B1_mc3 -3
 T B1_mc4 -1
 B2 SpecP1 1
 B2 B2_mc1 1
 B2 B2_mc2 1
 B2 B2_mc3 1
 B2 B2_mc4 1
 B2 QualLink -1
 Z1 SpecP1 -0.5
 Z1 DemP1 1
 Z1 COST 1
 B3 SpecP2 1
 B3 B3_mc1 1
 B3 B3_mc2 1
 B3 B3_mc3 1
 B3 B3_mc4 1
 B3 QualLink -1
 Z2 SpecP2 0.5
 Z2 DemP2 1
 Z2 COST -5
 P B1_mc2 -1400
 P B1_mc4 -1400
 P B2_mc2 -600
 P B2_mc4 -600
 P B3_mc2 -800
 P B3_mc4 -800
RHS
 RHS1 DemP1 600
 RHS1 DemP2 800
 RHS1 B1_mc2 -4200
 RHS1 B1_mc4 -1400
 RHS1 B2_mc2 -1800
 RHS1 B2_mc4 -600
 RHS1 B3_mc2 -2400
 RHS1 B3_mc4 -800
BOUNDS
 LO BND P 1
 UP BND P 3
 UP BND Y 600
 UP BND W 800
QUADOBJ
 P P 2
ENDATA
