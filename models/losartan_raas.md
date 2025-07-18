# model: losartan_raas
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
ALDSEC_k = 1.0105027952686e-06  # [mmol/min] rate aldosterone secretion  
ANGGEN2ANG1_k = 0.100300143544081  # [l/min] rate angen to ang1 conversion  
BP_ald_fe = 0.311999262694736  # [-] effect of aldosterone on blood pressure [-]  
DBP_ref = 80.0  # [133.32239 N/m^2] reference diastolic blood pressure [mmHg]  
E50_e3174 = 0.000291072740111984  # [mmol/l] half-maximum effect concentration exp3174  
RENSEC_fa = 5.0  # [-] activation renin secretion via exp3174  
RENSEC_k = 0.1  # [mmol/min] rate renin secretion  
SBP_ref = 120.0  # [133.32239 N/m^2] reference systolic blood pressure [mmHg]  
Vplasma = 5.0  # [l] plasma  
ald_ref = 2.5e-07  # [mmol/l] reference concentration of aldosterone  
ang1_ref = 1e-08  # [mmol/l] reference concentration of angiotensin I  
ang2_ref = 8e-09  # [mmol/l] reference concentration of angiotensin II  
anggen_ref = 1e-08  # [mmol/l] reference concentration of angiotensinogen  
ren_ref = 1e-09  # [mmol/l] reference concentration of renin  
```

## Initial conditions `x0`
```
ald = <libsbml.ASTNode; proxy of <Swig Object of type 'ASTNode *' at 0x737bbffbd2f0> >  # [mmol/l] aldosterone in Vplasma  
ang1 = <libsbml.ASTNode; proxy of <Swig Object of type 'ASTNode *' at 0x737bbffbd200> >  # [mmol/l] angiotensin I in Vplasma  
ang2 = <libsbml.ASTNode; proxy of <Swig Object of type 'ASTNode *' at 0x737bbffbc510> >  # [mmol/l] angiotensin II in Vplasma  
anggen = <libsbml.ASTNode; proxy of <Swig Object of type 'ASTNode *' at 0x737bbff7a8b0> >  # [mmol/l] angiotensinogen in Vplasma  
e3174 = 0.0  # [mmol/l] e3174 in Vplasma  
ren = <libsbml.ASTNode; proxy of <Swig Object of type 'ASTNode *' at 0x737bbff97960> >  # [mmol/l] renin in Vplasma  
```

## ODE system
```
# y
ALDDEG_k = ALDSEC_k / ald_ref  # [l/min] rate aldosterone degradation  
ANG1ANG2_k = ANGGEN2ANG1_k * anggen_ref / ang1_ref  # [l/min] rate ang1 to ang2 conversion  
ANG2DEG_k = ANGGEN2ANG1_k * anggen_ref / ang2_ref  # [l/min] rate aldosterone degradation  
ANGGEN2ANG1 = ANGGEN2ANG1_k * anggen * ren / ren_ref  # [mmol/min] angiotensinogen to angiotensin I (renin)  
DBP = DBP_ref + BP_ald_fe * DBP_ref * (ald - ald_ref) / ald_ref  # [133.32239 N/m^2] diastolic blood pressure [mmHg]  
RENDEG_k = RENSEC_k / ren_ref  # [l/min] rate renin degradation  
SBP = SBP_ref + BP_ald_fe * SBP_ref * (ald - ald_ref) / ald_ref  # [133.32239 N/m^2] systolic blood pressure [mmHg]  
ald_change = ald - ald_ref  # [mmol/l] aldosterone change  
ald_ratio = ald / ald_ref  # [-] aldosterone ratio  
ang1_change = ang1 - ang1_ref  # [mmol/l] angiotensin I change  
ang1_ratio = ang1 / ang1_ref  # [-] angiotensin I ratio  
ang2_change = ang2 - ang2_ref  # [mmol/l] angiotensin II change  
ang2_ratio = ang2 / ang2_ref  # [-] angiotensin II ratio  
fe_e3174 = e3174 / (E50_e3174 + e3174)  # [-] effect via exp3174 in [0, 1]  
ren_change = ren - ren_ref  # [mmol/l] renin change  
ren_ratio = ren / ren_ref  # [-] renin ratio  
ALDDEG = ALDDEG_k * ald  # [mmol/min] aldosterone degradation (ALDDEG)  
ALDSEC = ALDSEC_k * (1 - fe_e3174)  # [mmol/min] aldosterone secretion (ALDSEC)  
ANG1ANG2 = ANG1ANG2_k * ang1  # [mmol/min] angiotensin I to angiotensin II (ACE)  
ANG2DEG = ANG2DEG_k * ang2  # [mmol/min] angiotensin II degradation (ANG2DEG)  
DBP_change = DBP - DBP_ref  # [133.32239 N/m^2] e3174 change  
DBP_ratio = DBP / DBP_ref  # [-] e3174 ratio  
MAP = DBP + (SBP - DBP) / 3  # [133.32239 N/m^2] mean arterial pressure [mmHg]  
RENDEG = RENDEG_k * ren  # [mmol/min] renin degradation (RENDEG)  
RENSEC = RENSEC_k * (1 + RENSEC_fa * fe_e3174)  # [mmol/min] renin secretion (RENSEC)  
SBP_change = SBP - SBP_ref  # [133.32239 N/m^2] e3174 change  
SBP_ratio = SBP / SBP_ref  # [-] e3174 ratio  

# odes
d ald/dt = ALDSEC / Vplasma - ALDDEG / Vplasma  # [mmol/l/min] aldosterone  
d ang1/dt = ANGGEN2ANG1 / Vplasma - ANG1ANG2 / Vplasma  # [mmol/l/min] angiotensin I  
d ang2/dt = ANG1ANG2 / Vplasma - ANG2DEG / Vplasma  # [mmol/l/min] angiotensin II  
d anggen/dt = 0  # [mmol/l/min] angiotensinogen  
d e3174/dt = 0  # [mmol/l/min] e3174  
d ren/dt = RENSEC / Vplasma - RENDEG / Vplasma  # [mmol/l/min] renin  
```