# model: losartan_liver
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
E3174EX_k = 0.0112001182272967  # [1/min] rate E3174 export  
E3174L158_k = 0.00113400909936956  # [1/min] rate losartan conversion  
L158EX_k = 1.0  # [1/min] rate L158 export  
LOS2E3174_Km_los = 0.0011  # [mmol/l] Km losartan conversion  
LOS2E3174_Vmax = 0.00072514631731673  # [mmol/min/l] Vmax losartan conversion  
LOSIM_k = 1.0  # [1/min] rate losartan import  
MBIEX_k = 0.066315021777011  # [1/min] rate for export in bile  
Vapical = nan  # [m^2] apical membrane  
Vbi = 1.0  # [l] bile  
Vext = 1.5  # [l] plasma  
Vli = 1.5  # [l] liver  
Vlumen = 1.15425  # [l] intestinal lumen (inner part of intestine)  
Vmem = nan  # [m^2] plasma membrane  
f_cyp2c9 = 1.0  # [-] activity of CYP2C9  
```

## Initial conditions `x0`
```
e3174 = 0.0  # [mmol/l] E3174 (liver) in Vli  
e3174_bi = 0.0  # [mmol] E3174 (bile) in Vbi  
e3174_ext = 0.0  # [mmol/l] E3174 (plasma) in Vext  
e3174_lumen = 0.0  # [mmol/l] E3174 (lumen) in Vlumen  
l158 = 0.0  # [mmol/l] L158 (liver) in Vli  
l158_bi = 0.0  # [mmol] L158 (bile) in Vbi  
l158_ext = 0.0  # [mmol/l] L158 (plasma) in Vext  
l158_lumen = 0.0  # [mmol/l] L158 (lumen) in Vlumen  
los = 0.0  # [mmol/l] losartan (liver) in Vli  
los_bi = 0.0  # [mmol] losartan (bile) in Vbi  
los_ext = 0.0  # [mmol/l] losartan (plasma) in Vext  
los_lumen = 0.0  # [mmol/l] losartan (lumen) in Vlumen  
```

## ODE system
```
# y
E3174BIEX = MBIEX_k * Vli * e3174  # [mmol/min] E3174 bile export  
E3174EX = E3174EX_k * Vli * (e3174 - e3174_ext)  # [mmol/min] E3174 export (E3174EX)  
E3174L158 = E3174L158_k * Vli * e3174  # [mmol/min] L158 formation (UGT)  
L158BIEX = MBIEX_k * Vli * l158  # [mmol/min] L158 bile export  
L158EX = L158EX_k * Vli * (l158 - l158_ext)  # [mmol/min] L158 export (L158EX)  
LOS2E3174 = f_cyp2c9 * LOS2E3174_Vmax * Vli * los / (los + LOS2E3174_Km_los)  # [mmol/min] E3174 formation (CYP2C9/3A4)  
LOSBIEX = MBIEX_k * Vli * los  # [mmol/min] Losartan bile export  
LOSIM = LOSIM_k * Vli * (los_ext - los)  # [mmol/min] losartan import (LOSIM)  
E3174EHC = E3174BIEX  # [mmol/min] E3174 enterohepatic circulation  
L158EHC = L158BIEX  # [mmol/min] L158 enterohepatic circulation  
LOSEHC = LOSBIEX  # [mmol/min] Losartan enterohepatic circulation  

# odes
d e3174/dt = LOS2E3174 / Vli - E3174L158 / Vli - E3174EX / Vli - E3174BIEX / Vli  # [mmol/l/min] E3174 (liver)  
d e3174_bi/dt = E3174BIEX - E3174EHC  # [mmol/min] E3174 (bile)  
d e3174_ext/dt = E3174EX / Vext  # [mmol/l/min] E3174 (plasma)  
d e3174_lumen/dt = E3174EHC / Vlumen  # [mmol/l/min] E3174 (lumen)  
d l158/dt = E3174L158 / Vli - L158EX / Vli - L158BIEX / Vli  # [mmol/l/min] L158 (liver)  
d l158_bi/dt = L158BIEX - L158EHC  # [mmol/min] L158 (bile)  
d l158_ext/dt = L158EX / Vext  # [mmol/l/min] L158 (plasma)  
d l158_lumen/dt = L158EHC / Vlumen  # [mmol/l/min] L158 (lumen)  
d los/dt = LOSIM / Vli - LOS2E3174 / Vli - LOSBIEX / Vli  # [mmol/l/min] losartan (liver)  
d los_bi/dt = LOSBIEX - LOSEHC  # [mmol/min] losartan (bile)  
d los_ext/dt = -LOSIM / Vext  # [mmol/l/min] losartan (plasma)  
d los_lumen/dt = LOSEHC / Vlumen  # [mmol/l/min] losartan (lumen)  
```