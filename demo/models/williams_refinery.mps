NAME WILLIAMS_REFINERY
ROWS
 N COST
 L DistillCap
 L ReformCap
 L CrackCap
 E Yield_LN
 E Yield_MN
 E Yield_HN
 E Yield_LO
 E Yield_HO
 E Yield_RES
 E Yield_RG
 E Yield_CG
 E Yield_CO
 E Yield_LUBE
 E Sum_PF
 E Sum_RF
 E Sum_JF
 E Mass_LN
 E Mass_MN
 E Mass_HN
 E Mass_LO
 E Mass_HO
 E Mass_CO
 E Mass_RES
 E Mass_CG
 E Mass_RG
 G Prem2Reg
 G Octane_RF
 G Octane_PF
 L VaporPress
COLUMNS
 CR1 DistillCap 1
 CR1 Yield_LN 0.1
 CR1 Yield_MN 0.2
 CR1 Yield_HN 0.2
 CR1 Yield_LO 0.12
 CR1 Yield_HO 0.2
 CR1 Yield_RES 0.13
 CR2 DistillCap 1
 CR2 Yield_LN 0.15
 CR2 Yield_MN 0.25
 CR2 Yield_HN 0.18
 CR2 Yield_LO 0.08
 CR2 Yield_HO 0.19
 CR2 Yield_RES 0.12
 RU_LN ReformCap 1
 RU_LN Yield_RG 0.6
 RU_LN Mass_LN 1
 RU_MN ReformCap 1
 RU_MN Yield_RG 0.52
 RU_MN Mass_MN 1
 RU_HN ReformCap 1
 RU_HN Yield_RG 0.45
 RU_HN Mass_HN 1
 CLO CrackCap 1
 CLO Yield_CG 0.28
 CLO Yield_CO 0.68
 CLO Mass_LO 1
 CHO CrackCap 1
 CHO Yield_CG 0.2
 CHO Yield_CO 0.75
 CHO Mass_HO 1
 LN Yield_LN -1
 LN Mass_LN -1
 MN Yield_MN -1
 MN Mass_MN -1
 HN Yield_HN -1
 HN Mass_HN -1
 LO Yield_LO -1
 LO Mass_LO -1
 HO Yield_HO -1
 HO Mass_HO -1
 RES Yield_RES -1
 RES Mass_RES -1
 RG Yield_RG -1
 RG Mass_RG -1
 CG Yield_CG -1
 CG Mass_CG -1
 CO Yield_CO -1
 CO Mass_CO -1
 LUBIN Yield_LUBE 0.5
 LUBIN Mass_RES 1
 LUBE Yield_LUBE -1
 LUBE COST -1.5
 PLN Sum_PF 1
 PLN Mass_LN 1
 PLN Octane_PF 90
 PMN Sum_PF 1
 PMN Mass_MN 1
 PMN Octane_PF 80
 PHN Sum_PF 1
 PHN Mass_HN 1
 PHN Octane_PF 70
 PRG Sum_PF 1
 PRG Mass_RG 1
 PRG Octane_PF 115
 PCG Sum_PF 1
 PCG Mass_CG 1
 PCG Octane_PF 105
 PF Sum_PF -1
 PF Prem2Reg 1
 PF Octane_PF -94
 PF COST -7
 RLN Sum_RF 1
 RLN Mass_LN 1
 RLN Octane_RF 90
 RMN Sum_RF 1
 RMN Mass_MN 1
 RMN Octane_RF 80
 RHN Sum_RF 1
 RHN Mass_HN 1
 RHN Octane_RF 70
 RRG Sum_RF 1
 RRG Mass_RG 1
 RRG Octane_RF 115
 RCG Sum_RF 1
 RCG Mass_CG 1
 RCG Octane_RF 105
 RF Sum_RF -1
 RF Prem2Reg -0.4
 RF Octane_RF -84
 RF COST -6
 JLO Sum_JF 1
 JLO Mass_LO 1
 JLO VaporPress 1
 JHO Sum_JF 1
 JHO Mass_HO 1
 JHO VaporPress 0.6
 JRES Sum_JF 1
 JRES Mass_RES 1
 JRES VaporPress 0.05
 JCO Sum_JF 1
 JCO Mass_CO 1
 JCO VaporPress 1.5
 JF Sum_JF -1
 JF VaporPress -1
 JF COST -4
 FO Mass_LO 0.55
 FO Mass_HO 0.17
 FO Mass_CO 0.22
 FO Mass_RES 0.055
 FO COST -3.5
RHS
 RHS1 DistillCap 45000
 RHS1 ReformCap 10000
 RHS1 CrackCap 8000
BOUNDS
 UP BND CR1 20000
 UP BND CR2 30000
 LO BND LUBE 500
 UP BND LUBE 1000
ENDATA
