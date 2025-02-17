# model: losartan_body
Autogenerated ODE System from SBML with [sbmlutils](https://github.com/matthiaskoenig/sbmlutils).
```
time: [min]
substance: [mmol]
extent: [mmol]
volume: [l]
area: [m^2]
length: [m]
```

## Parameters `p`
```
BW = 75.0  # [kg] body weight [kg]  
COBW = 1.548  # [ml/s/kg] cardiac output per bodyweight [ml/s/kg]  
COHRI = 150.0  # [ml] increase of cardiac output per heartbeat [ml/min*min]  
FQgu = 0.18  # [-] gut fractional tissue blood flow  
FQh = 0.215  # [-] hepatic (venous side) fractional tissue blood flow  
FQki = 0.19  # [-] kidney fractional tissue blood flow  
FQlu = 1.0  # [-] lung fractional tissue blood flow  
FVar = 0.0257  # [l/kg] arterial fractional tissue volume  
FVgu = 0.0171  # [l/kg] gut fractional tissue volume  
FVhv = 0.001  # [l/kg] hepatic venous fractional tissue volume  
FVki = 0.0044  # [l/kg] kidney fractional tissue volume  
FVli = 0.021  # [l/kg] liver fractional tissue volume  
FVlu = 0.0076  # [l/kg] lung fractional tissue volume  
FVpo = 0.001  # [l/kg] portal fractional tissue volume  
FVve = 0.0514  # [l/kg] venous fractional tissue volume  
Fblood = 0.02  # [-] blood fraction of organ volume  
GU__F_los_abs = 0.6  # [-] fraction absorbed losartan  
GU__Ka_dis_los = 2.0  # [1/hr] Ka_dis [1/hr] dissolution losartan  
GU__LOSABS_k = 0.01  # [1/min] rate of losartan absorption  
GU__LOSEFL_Km_los = 0.01  # [mmol/l] Km for losartan efflux (PG)  
GU__LOSEFL_Vmax = 1.0  # [mmol/min/l] Vmax for losartan efflux (PG)  
GU__LOSEX_Km_los = 0.01  # [mmol/l] Km for losartan export in plasma  
GU__LOSEX_Vmax = 1.0  # [mmol/min/l] Vmax for losartan export in plasma  
GU__METEXC_k = 0.01  # [1/min] rate of metabolite feces excretion  
GU__Mr_los = 422.911  # [g/mol] Molecular weight losartan [g/mole]  
GU__Vapical = nan  # [m^2] apical membrane (intestinal membrane enterocytes)  
GU__Vbaso = nan  # [m^2] basolateral membrane (intestinal membrane enterocytes)  
GU__Ventero = 1.0  # [l] intestinal lining (enterocytes)  
GU__Vlumen = 1.15425  # [l] intestinal lumen (inner part of intestine)  
GU__Vstomach = 1.0  # [l] stomach  
GU__f_OATP2B1 = 1.0  # [-] parameter for OATP2B1 activity  
GU__f_abcb1 = 1.0  # [-] parameter for P-glycoprotein activity  
HCT = 0.51  # [-] hematocrit  
HEIGHT = 170.0  # [cm] height [cm]  
HR = 70.0  # [1/min] heart rate [1/min]  
HRrest = 70.0  # [1/min] heart rate [1/min]  
KI__E3174EX_k = 0.2  # [1/min] rate urinary excretion of e3174  
KI__L158EX_k = 0.5  # [1/min] rate urinary excretion of l158  
KI__LOSEX_k = 0.1  # [1/min] rate urinary excretion of losartan  
KI__Vmem = nan  # [m^2] plasma membrane  
KI__f_renal_function = 1.0  # [-] parameter for renal function  
LI__E3174EX_Km_e3174 = 0.01  # [mmol/l] Km E3174 export  
LI__E3174EX_Vmax = 1.0  # [mmol/min/l] Vmax E3174 export  
LI__L158EX_Km_l158 = 0.01  # [mmol/l] Km L158 export  
LI__L158EX_Vmax = 1.0  # [mmol/min/l] Vmax L158 export  
LI__LOS2E3174_Km_los = 0.0011  # [mmol/l] Km losartan conversion  
LI__LOS2E3174_Vmax = 0.01  # [mmol/min/l] Vmax losartan conversion  
LI__LOS2L158_Km_los = 0.1  # [mmol/l] Km losartan conversion  
LI__LOS2L158_Vmax = 0.02  # [mmol/min/l] Vmax losartan conversion  
LI__LOSIM_Km_los = 0.1  # [mmol/l] Km losartan import  
LI__LOSIM_Vmax = 0.01  # [mmol/min/l] Vmax losartan import  
LI__MBIEX_k = 0.0001  # [1/min] rate for export in bile  
LI__Vapical = nan  # [m^2] apical membrane  
LI__Vbi = 1.0  # [l] bile  
LI__Vlumen = 1.15425  # [l] intestinal lumen (inner part of intestine)  
LI__Vmem = nan  # [m^2] plasma membrane  
LI__f_cyp2c9 = 1.0  # [-] activity of CYP2C9  
MAP = 100.0  # [133.32239 N/m^2] mean arterial pressure [mmHg]  
Mr_e3174 = 436.894  # [g/mol] Molecular weight e3174 [g/mole]  
Mr_l158 = 436.894  # [g/mol] Molecular weight l158 [g/mole]  
Mr_los = 422.911  # [g/mol] Molecular weight los [g/mole]  
Ri_e3174 = 0.0  # [mg/min] Ri [mg/min] rate of infusion e3174  
Ri_los = 0.0  # [mg/min] Ri [mg/min] rate of infusion los  
Vfeces = 1.0  # [l] feces  
Vstomach = 1.0  # [l] stomach  
Vurine = 1.0  # [l] urine  
conversion_min_per_day = 1440.0  # [min/day] Conversion factor min to hours  
f_cardiac_function = 1.0  # [-] heart function  
f_cirrhosis = 0.0  # [-] severity of cirrhosis [0, 0.95]  
mr_los_e317s_feces = 0.0  # [-]   
mr_los_e317s_urine = 0.0  # [-]   
ti_e3174 = 10.0  # [s] injection time e3174 [s]  
ti_los = 10.0  # [s] injection time los [s]  
```

## Initial conditions `x0`
```
Afeces_e3174 = 0.0  # [mmol] E3174 (feces) in Vfeces  
Afeces_l158 = 0.0  # [mmol] L158 (feces) in Vfeces  
Afeces_los = 0.0  # [mmol] losartan (feces) in Vfeces  
Aurine_e3174 = 0.0  # [mmol] E3174 (urine) in Vurine  
Aurine_l158 = 0.0  # [mmol] L158 (urine) in Vurine  
Aurine_los = 0.0  # [mmol] losartan (urine) in Vurine  
Car_e3174 = 0.0  # [mmol/l] E3174 (arterial blood plasma) in Var  
Car_l158 = 0.0  # [mmol/l] L158 (arterial blood plasma) in Var  
Car_los = 0.0  # [mmol/l] losartan (arterial blood plasma) in Var  
Cgu_e3174 = 0.0  # [mmol/l] E3174 (gut) in Vgu  
Cgu_l158 = 0.0  # [mmol/l] L158 (gut) in Vgu  
Cgu_los = 0.0  # [mmol/l] losartan (gut) in Vgu  
Cgu_plasma_e3174 = 0.0  # [mmol/l] E3174 (gut plasma) in Vgu_plasma  
Cgu_plasma_l158 = 0.0  # [mmol/l] L158 (gut plasma) in Vgu_plasma  
Cgu_plasma_los = 0.0  # [mmol/l] losartan (gut plasma) in Vgu_plasma  
Chv_e3174 = 0.0  # [mmol/l] E3174 (hepatic vein plasma) in Vhv  
Chv_l158 = 0.0  # [mmol/l] L158 (hepatic vein plasma) in Vhv  
Chv_los = 0.0  # [mmol/l] losartan (hepatic vein plasma) in Vhv  
Cki_plasma_e3174 = 0.0  # [mmol/l] E3174 (kidney plasma) in Vki_plasma  
Cki_plasma_l158 = 0.0  # [mmol/l] L158 (kidney plasma) in Vki_plasma  
Cki_plasma_los = 0.0  # [mmol/l] losartan (kidney plasma) in Vki_plasma  
Cli_plasma_e3174 = 0.0  # [mmol/l] E3174 (liver plasma) in Vli_plasma  
Cli_plasma_l158 = 0.0  # [mmol/l] L158 (liver plasma) in Vli_plasma  
Cli_plasma_los = 0.0  # [mmol/l] losartan (liver plasma) in Vli_plasma  
Clu_plasma_e3174 = 0.0  # [mmol/l] E3174 (lung plasma) in Vlu_plasma  
Clu_plasma_l158 = 0.0  # [mmol/l] L158 (lung plasma) in Vlu_plasma  
Clu_plasma_los = 0.0  # [mmol/l] losartan (lung plasma) in Vlu_plasma  
Cpo_e3174 = 0.0  # [mmol/l] E3174 (portal vein plasma) in Vpo  
Cpo_l158 = 0.0  # [mmol/l] L158 (portal vein plasma) in Vpo  
Cpo_los = 0.0  # [mmol/l] losartan (portal vein plasma) in Vpo  
Cre_plasma_e3174 = 0.0  # [mmol/l] E3174 (rest plasma) in Vre_plasma  
Cre_plasma_l158 = 0.0  # [mmol/l] L158 (rest plasma) in Vre_plasma  
Cre_plasma_los = 0.0  # [mmol/l] losartan (rest plasma) in Vre_plasma  
Cve_e3174 = 0.0  # [mmol/l] E3174 (venous blood plasma) in Vve  
Cve_l158 = 0.0  # [mmol/l] L158 (venous blood plasma) in Vve  
Cve_los = 0.0  # [mmol/l] losartan (venous blood plasma) in Vve  
GU__los = 0.0  # [mmol/l] losartan (enterocytes) in GU__Ventero  
GU__los_stomach = 0.0  # [mmol] losartan (stomach) in GU__Vstomach  
IVDOSE_e3174 = 0.0  # [mg] IV bolus dose e3174 [mg]  
IVDOSE_los = 0.0  # [mg] IV bolus dose los [mg]  
LI__e3174 = 0.0  # [mmol/l] E3174 (liver) in Vli_tissue  
LI__e3174_bi = 0.0  # [mmol] E3174 (bile) in LI__Vbi  
LI__l158 = 0.0  # [mmol/l] L158 (liver) in Vli_tissue  
LI__l158_bi = 0.0  # [mmol] L158 (bile) in LI__Vbi  
LI__los = 0.0  # [mmol/l] losartan (liver) in Vli_tissue  
LI__los_bi = 0.0  # [mmol] losartan (bile) in LI__Vbi  
PODOSE_los = 0.0  # [mg] oral dose los [mg]  
cum_dose_e3174 = 0.0  # [mg] Cumulative dose due to infusion e3174  
cum_dose_los = 0.0  # [mg] Cumulative dose due to infusion los  
```

## ODE system
```
# y
BSA = 0.024265 * (BW / 1)**0.5378 * (HEIGHT / 1)**0.3964  # [m^2] body surface area [m^2]  
CO = f_cardiac_function * BW * COBW + (HR - HRrest) * COHRI / 60  # [ml/s] cardiac output [ml/s]  
FQre = 1 - (FQki + FQh)  # [-] rest of body fractional tissue blood flow  
FVre = 1 - (FVgu + FVki + FVli + FVlu + FVve + FVar)  # [l/kg] rest of body fractional tissue volume  
GU__dissolution_los = (GU__Ka_dis_los / 60) * PODOSE_los / GU__Mr_los  # [mmol/min] dissolution losartan  
Ki_e3174 = (0.693 / ti_e3174) * 60  # [1/min] injection rate IV  
Ki_los = (0.693 / ti_los) * 60  # [1/min] injection rate IV  
Var = BW * FVar - (FVar / (FVar + FVve)) * BW * Fblood * (1 - FVve - FVar)  # [l] arterial blood  
Vgu = BW * FVgu  # [l] gut  
Vhv = (1 - HCT) * (BW * FVhv - (FVhv / (FVar + FVve + FVpo + FVhv)) * BW * Fblood * (1 - (FVar + FVve + FVpo + FVhv)))  # [l] hepatic venous plasma  
Vki = BW * FVki  # [l] kidney  
Vli = BW * FVli  # [l] liver  
Vlu = BW * FVlu  # [l] lung  
Vpo = (1 - HCT) * (BW * FVpo - (FVpo / (FVar + FVve + FVpo + FVhv)) * BW * Fblood * (1 - (FVar + FVve + FVpo + FVhv)))  # [l] portal plasma  
Vve = BW * FVve - (FVve / (FVar + FVve)) * BW * Fblood * (1 - FVve - FVar)  # [l] venous blood  
Xfeces_e3174 = Afeces_e3174 * Mr_e3174  # [mg] E3174 amount (feces)  
Xfeces_l158 = Afeces_l158 * Mr_l158  # [mg] L158 amount (feces)  
Xfeces_los = Afeces_los * Mr_los  # [mg] losartan amount (feces)  
Xurine_e3174 = Aurine_e3174 * Mr_e3174  # [mg] E3174 amount (urine) [mg]  
Xurine_l158 = Aurine_l158 * Mr_l158  # [mg] L158 amount (urine) [mg]  
Xurine_los = Aurine_los * Mr_los  # [mg] losartan amount (urine) [mg]  
f_shunts = f_cirrhosis  # [-] fraction of portal venous blood shunted by the liver  
f_tissue_loss = f_cirrhosis  # [-] fraction of lost parenchymal liver volume  
mr_e3174_los_feces = Afeces_e3174 / (Afeces_los + 1e-12)  # [-]   
mr_e3174_los_plasma = Cve_e3174 / (Cve_los + 1e-12)  # [-]   
mr_e3174_los_urine = Aurine_e3174 / (Aurine_los + 1e-12)  # [-]   
mr_los_e3174_feces = Afeces_los / (Afeces_e3174 + 1e-12)  # [-]   
mr_los_e3174_plasma = Cve_los / (Cve_e3174 + 1e-12)  # [-]   
mr_los_e3174_urine = Aurine_los / (Aurine_e3174 + 1e-12)  # [-]   
Aar_e3174 = Car_e3174 * Var  # [mmol] E3174 amount (arterial blood) [mmole]  
Aar_l158 = Car_l158 * Var  # [mmol] L158 amount (arterial blood) [mmole]  
Aar_los = Car_los * Var  # [mmol] losartan amount (arterial blood) [mmole]  
Ahv_e3174 = Chv_e3174 * Vhv  # [mmol] E3174 amount (hepatic vein) [mmole]  
Ahv_l158 = Chv_l158 * Vhv  # [mmol] L158 amount (hepatic vein) [mmole]  
Ahv_los = Chv_los * Vhv  # [mmol] losartan amount (hepatic vein) [mmole]  
Apo_e3174 = Cpo_e3174 * Vpo  # [mmol] E3174 amount (portal vein) [mmole]  
Apo_l158 = Cpo_l158 * Vpo  # [mmol] L158 amount (portal vein) [mmole]  
Apo_los = Cpo_los * Vpo  # [mmol] losartan amount (portal vein) [mmole]  
Ave_e3174 = Cve_e3174 * Vve  # [mmol] E3174 amount (venous blood) [mmole]  
Ave_l158 = Cve_l158 * Vve  # [mmol] L158 amount (venous blood) [mmole]  
Ave_los = Cve_los * Vve  # [mmol] losartan amount (venous blood) [mmole]  
GU__E3174EXC = GU__METEXC_k * Vgu * Cgu_e3174  # [mmol/min] excretion e3174 (feces)  
GU__L158EXC = GU__METEXC_k * Vgu * Cgu_l158  # [mmol/min] excretion l158 (feces)  
GU__LOSEFL = GU__f_abcb1 * GU__LOSEFL_Vmax * Vgu * GU__los / (GU__los + GU__LOSEFL_Km_los)  # [mmol/min] efflux losartan (PG)  
GU__LOSEX = GU__LOSEX_Vmax * Vgu * GU__los / (GU__los + GU__LOSEX_Km_los)  # [mmol/min] losartan export plasma  
GU__absorption = GU__LOSABS_k * Vgu * Cgu_los  # [mmol/min] absorption losartan  
QC = (CO / 1000) * 60  # [l/min] cardiac output [L/hr]  
Vgu_plasma = Vgu * Fblood * (1 - HCT)  # [l] plasma volume of gut  
Vgu_tissue = Vgu * (1 - Fblood)  # [l] tissue volume of gut  
Vki_plasma = Vki * Fblood * (1 - HCT)  # [l] plasma volume of kidney  
Vki_tissue = Vki * (1 - Fblood)  # [l] tissue volume of kidney  
Vli_plasma = Vli * Fblood * (1 - HCT)  # [l] plasma volume of liver  
Vli_tissue = Vli * (1 - f_tissue_loss) * (1 - Fblood)  # [l] tissue volume of liver  
Vlu_plasma = Vlu * Fblood * (1 - HCT)  # [l] plasma volume of lung  
Vlu_tissue = Vlu * (1 - Fblood)  # [l] tissue volume of lung  
Vre = BW * FVre  # [l] rest of body  
iv_e3174 = Ki_e3174 * IVDOSE_e3174 / Mr_e3174  # [mmol/min] iv E3174  
iv_los = Ki_los * IVDOSE_los / Mr_los  # [mmol/min] iv losartan  
Agu_plasma_e3174 = Cgu_plasma_e3174 * Vgu_plasma  # [mmol] E3174 amount (gut) [mmole]  
Agu_plasma_l158 = Cgu_plasma_l158 * Vgu_plasma  # [mmol] L158 amount (gut) [mmole]  
Agu_plasma_los = Cgu_plasma_los * Vgu_plasma  # [mmol] losartan amount (gut) [mmole]  
Aki_plasma_e3174 = Cki_plasma_e3174 * Vki_plasma  # [mmol] E3174 amount (kidney) [mmole]  
Aki_plasma_l158 = Cki_plasma_l158 * Vki_plasma  # [mmol] L158 amount (kidney) [mmole]  
Aki_plasma_los = Cki_plasma_los * Vki_plasma  # [mmol] losartan amount (kidney) [mmole]  
Ali_plasma_e3174 = Cli_plasma_e3174 * Vli_plasma  # [mmol] E3174 amount (liver) [mmole]  
Ali_plasma_l158 = Cli_plasma_l158 * Vli_plasma  # [mmol] L158 amount (liver) [mmole]  
Ali_plasma_los = Cli_plasma_los * Vli_plasma  # [mmol] losartan amount (liver) [mmole]  
Alu_plasma_e3174 = Clu_plasma_e3174 * Vlu_plasma  # [mmol] E3174 amount (lung) [mmole]  
Alu_plasma_l158 = Clu_plasma_l158 * Vlu_plasma  # [mmol] L158 amount (lung) [mmole]  
Alu_plasma_los = Clu_plasma_los * Vlu_plasma  # [mmol] losartan amount (lung) [mmole]  
GU__LOSABS = GU__f_OATP2B1 * GU__F_los_abs * GU__absorption  # [mmol/min] absorption losartan  
GU__LOSEXC = (1 - GU__F_los_abs) * GU__absorption  # [mmol/min] excretion losartan (feces)  
KI__E3174EX = KI__f_renal_function * KI__E3174EX_k * Vki_tissue * Cki_plasma_e3174  # [mmol/min] e3174 excretion (E3174EX)  
KI__L158EX = KI__f_renal_function * KI__L158EX_k * Vki_tissue * Cki_plasma_l158  # [mmol/min] l158 excretion (L158EX)  
KI__LOSEX = KI__f_renal_function * KI__LOSEX_k * Vki_tissue * Cki_plasma_los  # [mmol/min] losartan excretion (LOSEX)  
LI__E3174BIEX = LI__MBIEX_k * Vli_tissue * LI__e3174  # [mmol/min] E3174 bile export  
LI__E3174EX = (LI__E3174EX_Vmax / LI__E3174EX_Km_e3174) * Vli_tissue * (LI__e3174 - Cli_plasma_e3174) / (1 + Cli_plasma_e3174 / LI__E3174EX_Km_e3174 + LI__e3174 / LI__E3174EX_Km_e3174)  # [mmol/min] E3174 export (E3174EX)  
LI__L158BIEX = LI__MBIEX_k * Vli_tissue * LI__l158  # [mmol/min] L158 bile export  
LI__L158EX = (LI__L158EX_Vmax / LI__L158EX_Km_l158) * Vli_tissue * (LI__l158 - Cli_plasma_l158) / (1 + Cli_plasma_l158 / LI__L158EX_Km_l158 + LI__l158 / LI__L158EX_Km_l158)  # [mmol/min] L158 export (L158EX)  
LI__LOS2E3174 = LI__f_cyp2c9 * LI__LOS2E3174_Vmax * Vli_tissue * LI__los / (LI__los + LI__LOS2E3174_Km_los)  # [mmol/min] E3174 formation (CYP2C9/3A4)  
LI__LOS2L158 = LI__LOS2L158_Vmax * Vli_tissue * LI__los / (LI__los + LI__LOS2L158_Km_los)  # [mmol/min] L158 formation (UGT)  
LI__LOSBIEX = LI__MBIEX_k * Vli_tissue * LI__los  # [mmol/min] Losartan bile export  
LI__LOSIM = (LI__LOSIM_Vmax / LI__LOSIM_Km_los) * Vli_tissue * (Cli_plasma_los - LI__los) / (1 + Cli_plasma_los / LI__LOSIM_Km_los + LI__los / LI__LOSIM_Km_los)  # [mmol/min] losartan import (LOSIM)  
Mar_e3174 = (Aar_e3174 / Var) * Mr_e3174  # [mg/l] E3174 concentration (arterial blood) [mg/l]  
Mar_l158 = (Aar_l158 / Var) * Mr_l158  # [mg/l] L158 concentration (arterial blood) [mg/l]  
Mar_los = (Aar_los / Var) * Mr_los  # [mg/l] losartan concentration (arterial blood) [mg/l]  
Mhv_e3174 = (Ahv_e3174 / Vhv) * Mr_e3174  # [mg/l] E3174 concentration (hepatic vein) [mg/l]  
Mhv_l158 = (Ahv_l158 / Vhv) * Mr_l158  # [mg/l] L158 concentration (hepatic vein) [mg/l]  
Mhv_los = (Ahv_los / Vhv) * Mr_los  # [mg/l] losartan concentration (hepatic vein) [mg/l]  
Mpo_e3174 = (Apo_e3174 / Vpo) * Mr_e3174  # [mg/l] E3174 concentration (portal vein) [mg/l]  
Mpo_l158 = (Apo_l158 / Vpo) * Mr_l158  # [mg/l] L158 concentration (portal vein) [mg/l]  
Mpo_los = (Apo_los / Vpo) * Mr_los  # [mg/l] losartan concentration (portal vein) [mg/l]  
Mve_e3174 = (Ave_e3174 / Vve) * Mr_e3174  # [mg/l] E3174 concentration (venous blood) [mg/l]  
Mve_l158 = (Ave_l158 / Vve) * Mr_l158  # [mg/l] L158 concentration (venous blood) [mg/l]  
Mve_los = (Ave_los / Vve) * Mr_los  # [mg/l] losartan concentration (venous blood) [mg/l]  
Qgu = QC * FQgu  # [l/min] gut blood flow  
Qh = QC * FQh  # [l/min] hepatic (venous side) blood flow  
Qki = QC * FQki  # [l/min] kidney blood flow  
Qlu = QC * FQlu  # [l/min] lung blood flow  
Qre = QC * FQre  # [l/min] rest of body blood flow  
Vre_plasma = Vre * Fblood * (1 - HCT)  # [l] plasma volume of rest  
Vre_tissue = Vre * (1 - Fblood)  # [l] tissue volume of rest  
Xar_e3174 = Aar_e3174 * Mr_e3174  # [mg] E3174 amount (arterial blood) [mg]  
Xar_l158 = Aar_l158 * Mr_l158  # [mg] L158 amount (arterial blood) [mg]  
Xar_los = Aar_los * Mr_los  # [mg] losartan amount (arterial blood) [mg]  
Xhv_e3174 = Ahv_e3174 * Mr_e3174  # [mg] E3174 amount (hepatic vein) [mg]  
Xhv_l158 = Ahv_l158 * Mr_l158  # [mg] L158 amount (hepatic vein) [mg]  
Xhv_los = Ahv_los * Mr_los  # [mg] losartan amount (hepatic vein) [mg]  
Xpo_e3174 = Apo_e3174 * Mr_e3174  # [mg] E3174 amount (portal vein) [mg]  
Xpo_l158 = Apo_l158 * Mr_l158  # [mg] L158 amount (portal vein) [mg]  
Xpo_los = Apo_los * Mr_los  # [mg] losartan amount (portal vein) [mg]  
Xve_e3174 = Ave_e3174 * Mr_e3174  # [mg] E3174 amount (venous blood) [mg]  
Xve_l158 = Ave_l158 * Mr_l158  # [mg] L158 amount (venous blood) [mg]  
Xve_los = Ave_los * Mr_los  # [mg] losartan amount (venous blood) [mg]  
Are_plasma_e3174 = Cre_plasma_e3174 * Vre_plasma  # [mmol] E3174 amount (rest) [mmole]  
Are_plasma_l158 = Cre_plasma_l158 * Vre_plasma  # [mmol] L158 amount (rest) [mmole]  
Are_plasma_los = Cre_plasma_los * Vre_plasma  # [mmol] losartan amount (rest) [mmole]  
Flow_ar_gu_e3174 = Qgu * Car_e3174  # [mmol/min] inflow gut E3174  
Flow_ar_gu_l158 = Qgu * Car_l158  # [mmol/min] inflow gut L158  
Flow_ar_gu_los = Qgu * Car_los  # [mmol/min] inflow gut losartan  
Flow_ar_ki_e3174 = Qki * Car_e3174  # [mmol/min] inflow kidney E3174  
Flow_ar_ki_l158 = Qki * Car_l158  # [mmol/min] inflow kidney L158  
Flow_ar_ki_los = Qki * Car_los  # [mmol/min] inflow kidney losartan  
Flow_ar_re_e3174 = Qre * Car_e3174  # [mmol/min] inflow rest E3174  
Flow_ar_re_l158 = Qre * Car_l158  # [mmol/min] inflow rest L158  
Flow_ar_re_los = Qre * Car_los  # [mmol/min] inflow rest losartan  
Flow_gu_po_e3174 = Qgu * Cgu_plasma_e3174  # [mmol/min] outflow gut E3174  
Flow_gu_po_l158 = Qgu * Cgu_plasma_l158  # [mmol/min] outflow gut L158  
Flow_gu_po_los = Qgu * Cgu_plasma_los  # [mmol/min] outflow gut losartan  
Flow_hv_ve_e3174 = Qh * Chv_e3174  # [mmol/min] outflow hepatic vein E3174  
Flow_hv_ve_l158 = Qh * Chv_l158  # [mmol/min] outflow hepatic vein L158  
Flow_hv_ve_los = Qh * Chv_los  # [mmol/min] outflow hepatic vein losartan  
Flow_ki_ve_e3174 = Qki * Cki_plasma_e3174  # [mmol/min] outflow kidney E3174  
Flow_ki_ve_l158 = Qki * Cki_plasma_l158  # [mmol/min] outflow kidney L158  
Flow_ki_ve_los = Qki * Cki_plasma_los  # [mmol/min] outflow kidney losartan  
Flow_lu_ar_e3174 = Qlu * Clu_plasma_e3174  # [mmol/min] outflow lung E3174  
Flow_lu_ar_l158 = Qlu * Clu_plasma_l158  # [mmol/min] outflow lung L158  
Flow_lu_ar_los = Qlu * Clu_plasma_los  # [mmol/min] outflow lung losartan  
Flow_re_ve_e3174 = Qre * Cre_plasma_e3174  # [mmol/min] outflow rest E3174  
Flow_re_ve_l158 = Qre * Cre_plasma_l158  # [mmol/min] outflow rest L158  
Flow_re_ve_los = Qre * Cre_plasma_los  # [mmol/min] outflow rest losartan  
Flow_ve_lu_e3174 = Qlu * Cve_e3174  # [mmol/min] inflow lung E3174  
Flow_ve_lu_l158 = Qlu * Cve_l158  # [mmol/min] inflow lung L158  
Flow_ve_lu_los = Qlu * Cve_los  # [mmol/min] inflow lung losartan  
LI__E3174EHC = LI__E3174BIEX  # [mmol/min] E3174 enterohepatic circulation  
LI__L158EHC = LI__L158BIEX  # [mmol/min] L158 enterohepatic circulation  
LI__LOSEHC = LI__LOSBIEX  # [mmol/min] Losartan enterohepatic circulation  
Mgu_plasma_e3174 = (Agu_plasma_e3174 / Vgu_plasma) * Mr_e3174  # [mg/l] E3174 concentration (gut) [mg/l]  
Mgu_plasma_l158 = (Agu_plasma_l158 / Vgu_plasma) * Mr_l158  # [mg/l] L158 concentration (gut) [mg/l]  
Mgu_plasma_los = (Agu_plasma_los / Vgu_plasma) * Mr_los  # [mg/l] losartan concentration (gut) [mg/l]  
Mki_plasma_e3174 = (Aki_plasma_e3174 / Vki_plasma) * Mr_e3174  # [mg/l] E3174 concentration (kidney) [mg/l]  
Mki_plasma_l158 = (Aki_plasma_l158 / Vki_plasma) * Mr_l158  # [mg/l] L158 concentration (kidney) [mg/l]  
Mki_plasma_los = (Aki_plasma_los / Vki_plasma) * Mr_los  # [mg/l] losartan concentration (kidney) [mg/l]  
Mli_plasma_e3174 = (Ali_plasma_e3174 / Vli_plasma) * Mr_e3174  # [mg/l] E3174 concentration (liver) [mg/l]  
Mli_plasma_l158 = (Ali_plasma_l158 / Vli_plasma) * Mr_l158  # [mg/l] L158 concentration (liver) [mg/l]  
Mli_plasma_los = (Ali_plasma_los / Vli_plasma) * Mr_los  # [mg/l] losartan concentration (liver) [mg/l]  
Mlu_plasma_e3174 = (Alu_plasma_e3174 / Vlu_plasma) * Mr_e3174  # [mg/l] E3174 concentration (lung) [mg/l]  
Mlu_plasma_l158 = (Alu_plasma_l158 / Vlu_plasma) * Mr_l158  # [mg/l] L158 concentration (lung) [mg/l]  
Mlu_plasma_los = (Alu_plasma_los / Vlu_plasma) * Mr_los  # [mg/l] losartan concentration (lung) [mg/l]  
Qha = Qh - Qgu  # [l/min] hepatic artery blood flow  
Qpo = Qgu  # [l/min] portal blood flow  
Xgu_plasma_e3174 = Agu_plasma_e3174 * Mr_e3174  # [mg] E3174 amount (gut) [mg]  
Xgu_plasma_l158 = Agu_plasma_l158 * Mr_l158  # [mg] L158 amount (gut) [mg]  
Xgu_plasma_los = Agu_plasma_los * Mr_los  # [mg] losartan amount (gut) [mg]  
Xki_plasma_e3174 = Aki_plasma_e3174 * Mr_e3174  # [mg] E3174 amount (kidney) [mg]  
Xki_plasma_l158 = Aki_plasma_l158 * Mr_l158  # [mg] L158 amount (kidney) [mg]  
Xki_plasma_los = Aki_plasma_los * Mr_los  # [mg] losartan amount (kidney) [mg]  
Xli_plasma_e3174 = Ali_plasma_e3174 * Mr_e3174  # [mg] E3174 amount (liver) [mg]  
Xli_plasma_l158 = Ali_plasma_l158 * Mr_l158  # [mg] L158 amount (liver) [mg]  
Xli_plasma_los = Ali_plasma_los * Mr_los  # [mg] losartan amount (liver) [mg]  
Xlu_plasma_e3174 = Alu_plasma_e3174 * Mr_e3174  # [mg] E3174 amount (lung) [mg]  
Xlu_plasma_l158 = Alu_plasma_l158 * Mr_l158  # [mg] L158 amount (lung) [mg]  
Xlu_plasma_los = Alu_plasma_los * Mr_los  # [mg] losartan amount (lung) [mg]  
Flow_arli_hv_e3174 = f_shunts * Qha * Car_e3174  # [mmol/min] flow arterial shunts  
Flow_arli_hv_l158 = f_shunts * Qha * Car_l158  # [mmol/min] flow arterial shunts  
Flow_arli_hv_los = f_shunts * Qha * Car_los  # [mmol/min] flow arterial shunts  
Flow_arli_li_e3174 = (1 - f_shunts) * Qha * Car_e3174  # [mmol/min] arterial inflow liver E3174  
Flow_arli_li_l158 = (1 - f_shunts) * Qha * Car_l158  # [mmol/min] arterial inflow liver L158  
Flow_arli_li_los = (1 - f_shunts) * Qha * Car_los  # [mmol/min] arterial inflow liver losartan  
Flow_li_hv_e3174 = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_e3174  # [mmol/min] outflow liver E3174  
Flow_li_hv_l158 = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_l158  # [mmol/min] outflow liver L158  
Flow_li_hv_los = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_los  # [mmol/min] outflow liver losartan  
Flow_po_hv_e3174 = f_shunts * Qpo * Cpo_e3174  # [mmol/min] flow portal shunts  
Flow_po_hv_l158 = f_shunts * Qpo * Cpo_l158  # [mmol/min] flow portal shunts  
Flow_po_hv_los = f_shunts * Qpo * Cpo_los  # [mmol/min] flow portal shunts  
Flow_po_li_e3174 = (1 - f_shunts) * Qpo * Cpo_e3174  # [mmol/min] outflow po E3174  
Flow_po_li_l158 = (1 - f_shunts) * Qpo * Cpo_l158  # [mmol/min] outflow po L158  
Flow_po_li_los = (1 - f_shunts) * Qpo * Cpo_los  # [mmol/min] outflow po losartan  
Mre_plasma_e3174 = (Are_plasma_e3174 / Vre_plasma) * Mr_e3174  # [mg/l] E3174 concentration (rest) [mg/l]  
Mre_plasma_l158 = (Are_plasma_l158 / Vre_plasma) * Mr_l158  # [mg/l] L158 concentration (rest) [mg/l]  
Mre_plasma_los = (Are_plasma_los / Vre_plasma) * Mr_los  # [mg/l] losartan concentration (rest) [mg/l]  
Xre_plasma_e3174 = Are_plasma_e3174 * Mr_e3174  # [mg] E3174 amount (rest) [mg]  
Xre_plasma_l158 = Are_plasma_l158 * Mr_l158  # [mg] L158 amount (rest) [mg]  
Xre_plasma_los = Are_plasma_los * Mr_los  # [mg] losartan amount (rest) [mg]  

# odes
d Afeces_e3174/dt = GU__E3174EXC  # [mmol/min] E3174 (feces)  
d Afeces_l158/dt = GU__L158EXC  # [mmol/min] L158 (feces)  
d Afeces_los/dt = GU__LOSEXC  # [mmol/min] losartan (feces)  
d Aurine_e3174/dt = KI__E3174EX  # [mmol/min] E3174 (urine)  
d Aurine_l158/dt = KI__L158EX  # [mmol/min] L158 (urine)  
d Aurine_los/dt = KI__LOSEX  # [mmol/min] losartan (urine)  
d Car_e3174/dt = (-Flow_ar_ki_e3174 / Var - Flow_arli_li_e3174 / Var - Flow_arli_hv_e3174 / Var) + Flow_lu_ar_e3174 / Var - Flow_ar_gu_e3174 / Var - Flow_ar_re_e3174 / Var  # [mmol/l/min] E3174 (arterial blood plasma)  
d Car_l158/dt = (-Flow_ar_ki_l158 / Var - Flow_arli_li_l158 / Var - Flow_arli_hv_l158 / Var) + Flow_lu_ar_l158 / Var - Flow_ar_gu_l158 / Var - Flow_ar_re_l158 / Var  # [mmol/l/min] L158 (arterial blood plasma)  
d Car_los/dt = (-Flow_ar_ki_los / Var - Flow_arli_li_los / Var - Flow_arli_hv_los / Var) + Flow_lu_ar_los / Var - Flow_ar_gu_los / Var - Flow_ar_re_los / Var  # [mmol/l/min] losartan (arterial blood plasma)  
d Cgu_e3174/dt = LI__E3174EHC / Vgu - GU__E3174EXC / Vgu  # [mmol/l/min] E3174 (gut)  
d Cgu_l158/dt = LI__L158EHC / Vgu - GU__L158EXC / Vgu  # [mmol/l/min] L158 (gut)  
d Cgu_los/dt = ((LI__LOSEHC / Vgu - GU__LOSABS / Vgu) + GU__LOSEFL / Vgu - GU__LOSEXC / Vgu) + GU__dissolution_los / Vgu  # [mmol/l/min] losartan (gut)  
d Cgu_plasma_e3174/dt = Flow_ar_gu_e3174 / Vgu_plasma - Flow_gu_po_e3174 / Vgu_plasma  # [mmol/l/min] E3174 (gut plasma)  
d Cgu_plasma_l158/dt = Flow_ar_gu_l158 / Vgu_plasma - Flow_gu_po_l158 / Vgu_plasma  # [mmol/l/min] L158 (gut plasma)  
d Cgu_plasma_los/dt = (Flow_ar_gu_los / Vgu_plasma - Flow_gu_po_los / Vgu_plasma) + GU__LOSEX / Vgu_plasma  # [mmol/l/min] losartan (gut plasma)  
d Chv_e3174/dt = Flow_arli_hv_e3174 / Vhv + Flow_po_hv_e3174 / Vhv + Flow_li_hv_e3174 / Vhv - Flow_hv_ve_e3174 / Vhv  # [mmol/l/min] E3174 (hepatic vein plasma)  
d Chv_l158/dt = Flow_arli_hv_l158 / Vhv + Flow_po_hv_l158 / Vhv + Flow_li_hv_l158 / Vhv - Flow_hv_ve_l158 / Vhv  # [mmol/l/min] L158 (hepatic vein plasma)  
d Chv_los/dt = Flow_arli_hv_los / Vhv + Flow_po_hv_los / Vhv + Flow_li_hv_los / Vhv - Flow_hv_ve_los / Vhv  # [mmol/l/min] losartan (hepatic vein plasma)  
d Cki_plasma_e3174/dt = Flow_ar_ki_e3174 / Vki_plasma - Flow_ki_ve_e3174 / Vki_plasma - KI__E3174EX / Vki_plasma  # [mmol/l/min] E3174 (kidney plasma)  
d Cki_plasma_l158/dt = Flow_ar_ki_l158 / Vki_plasma - Flow_ki_ve_l158 / Vki_plasma - KI__L158EX / Vki_plasma  # [mmol/l/min] L158 (kidney plasma)  
d Cki_plasma_los/dt = Flow_ar_ki_los / Vki_plasma - Flow_ki_ve_los / Vki_plasma - KI__LOSEX / Vki_plasma  # [mmol/l/min] losartan (kidney plasma)  
d Cli_plasma_e3174/dt = (Flow_arli_li_e3174 / Vli_plasma + Flow_po_li_e3174 / Vli_plasma - Flow_li_hv_e3174 / Vli_plasma) + LI__E3174EX / Vli_plasma  # [mmol/l/min] E3174 (liver plasma)  
d Cli_plasma_l158/dt = (Flow_arli_li_l158 / Vli_plasma + Flow_po_li_l158 / Vli_plasma - Flow_li_hv_l158 / Vli_plasma) + LI__L158EX / Vli_plasma  # [mmol/l/min] L158 (liver plasma)  
d Cli_plasma_los/dt = Flow_arli_li_los / Vli_plasma + Flow_po_li_los / Vli_plasma - Flow_li_hv_los / Vli_plasma - LI__LOSIM / Vli_plasma  # [mmol/l/min] losartan (liver plasma)  
d Clu_plasma_e3174/dt = Flow_ve_lu_e3174 / Vlu_plasma - Flow_lu_ar_e3174 / Vlu_plasma  # [mmol/l/min] E3174 (lung plasma)  
d Clu_plasma_l158/dt = Flow_ve_lu_l158 / Vlu_plasma - Flow_lu_ar_l158 / Vlu_plasma  # [mmol/l/min] L158 (lung plasma)  
d Clu_plasma_los/dt = Flow_ve_lu_los / Vlu_plasma - Flow_lu_ar_los / Vlu_plasma  # [mmol/l/min] losartan (lung plasma)  
d Cpo_e3174/dt = (-Flow_po_li_e3174 / Vpo - Flow_po_hv_e3174 / Vpo) + Flow_gu_po_e3174 / Vpo  # [mmol/l/min] E3174 (portal vein plasma)  
d Cpo_l158/dt = (-Flow_po_li_l158 / Vpo - Flow_po_hv_l158 / Vpo) + Flow_gu_po_l158 / Vpo  # [mmol/l/min] L158 (portal vein plasma)  
d Cpo_los/dt = (-Flow_po_li_los / Vpo - Flow_po_hv_los / Vpo) + Flow_gu_po_los / Vpo  # [mmol/l/min] losartan (portal vein plasma)  
d Cre_plasma_e3174/dt = Flow_ar_re_e3174 / Vre_plasma - Flow_re_ve_e3174 / Vre_plasma  # [mmol/l/min] E3174 (rest plasma)  
d Cre_plasma_l158/dt = Flow_ar_re_l158 / Vre_plasma - Flow_re_ve_l158 / Vre_plasma  # [mmol/l/min] L158 (rest plasma)  
d Cre_plasma_los/dt = Flow_ar_re_los / Vre_plasma - Flow_re_ve_los / Vre_plasma  # [mmol/l/min] losartan (rest plasma)  
d Cve_e3174/dt = (iv_e3174 / Vve + Flow_ki_ve_e3174 / Vve + Flow_hv_ve_e3174 / Vve - Flow_ve_lu_e3174 / Vve) + Flow_re_ve_e3174 / Vve  # [mmol/l/min] E3174 (venous blood plasma)  
d Cve_l158/dt = (Flow_ki_ve_l158 / Vve + Flow_hv_ve_l158 / Vve - Flow_ve_lu_l158 / Vve) + Flow_re_ve_l158 / Vve  # [mmol/l/min] L158 (venous blood plasma)  
d Cve_los/dt = (iv_los / Vve + Flow_ki_ve_los / Vve + Flow_hv_ve_los / Vve - Flow_ve_lu_los / Vve) + Flow_re_ve_los / Vve  # [mmol/l/min] losartan (venous blood plasma)  
d GU__los/dt = GU__LOSABS / GU__Ventero - GU__LOSEFL / GU__Ventero - GU__LOSEX / GU__Ventero  # [mmol/l/min] losartan (enterocytes)  
d GU__los_stomach/dt = 0  # [mmol/min] losartan (stomach)  
d IVDOSE_e3174/dt = -iv_e3174 * Mr_e3174 + Ri_e3174  # [mg/min] IV bolus dose e3174 [mg]  
d IVDOSE_los/dt = -iv_los * Mr_los + Ri_los  # [mg/min] IV bolus dose los [mg]  
d LI__e3174/dt = LI__LOS2E3174 / Vli_tissue - LI__E3174EX / Vli_tissue - LI__E3174BIEX / Vli_tissue  # [mmol/l/min] E3174 (liver)  
d LI__e3174_bi/dt = LI__E3174BIEX - LI__E3174EHC  # [mmol/min] E3174 (bile)  
d LI__l158/dt = LI__LOS2L158 / Vli_tissue - LI__L158EX / Vli_tissue - LI__L158BIEX / Vli_tissue  # [mmol/l/min] L158 (liver)  
d LI__l158_bi/dt = LI__L158BIEX - LI__L158EHC  # [mmol/min] L158 (bile)  
d LI__los/dt = LI__LOSIM / Vli_tissue - LI__LOS2E3174 / Vli_tissue - LI__LOS2L158 / Vli_tissue - LI__LOSBIEX / Vli_tissue  # [mmol/l/min] losartan (liver)  
d LI__los_bi/dt = LI__LOSBIEX - LI__LOSEHC  # [mmol/min] losartan (bile)  
d PODOSE_los/dt = -GU__dissolution_los * GU__Mr_los  # [mg/min] oral dose los [mg]  
d cum_dose_e3174/dt = Ri_e3174  # [mg/min] Cumulative dose due to infusion e3174  
d cum_dose_los/dt = Ri_los  # [mg/min] Cumulative dose due to infusion los  
```