# Model for hepatic losartan metabolism.

Model `losartan_liver`, SBML Level 3 Version 2, read from losartan\_liver.xml, written by sbmlode 0.3.0.

Model for the hepatic metabolism of losartan.

losartan exchanges fast between the blood of the liver and the hepatocytes (LOSIM)

losartan -&gt; E3174 (CYP2C9, LOS2E3174)

losartan -&gt; other metabolites (LOS2L158, CYP2C9 for the fraction f\_l158\_cyp2c9, CYP3A4 and UGTs for the rest)

export of E3174 and the other metabolites into the blood and the bile

The concentrations of the liver model are blood concentrations (the species of the whole-body model), the rate constants are intrinsic clearances of blood per liver volume. CYP2C9 forms about 55-60 % of the hepatic clearance: the loss of CYP2C9 (\*3/\*3, bucolome) increases the AUC of losartan 1.7-2.4 fold and decreases the AUC of E3174 to 0.1-0.2 \[Bae2011, Kobayashi2008\].

The content of this model has been carefully created in a manual research effort.

Terms of use

Copyright © 2026 Matthias König.

This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.

Redistribution and use of any part of this model, with or without modification, are permitted provided that the following conditions are met:

Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in a different form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Units

| Time | Substance | Extent | Volume | Area | Length |
| --- | --- | --- | --- | --- | --- |
| min | mmol | mmol | l | m<sup>2</sup> | m |

## Compartments

| Symbol | Id | Name | Size | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{Vext}$ | `Vext` | blood | $1.5$ | l | ✓ |
| $\mathrm{Vli}$ | `Vli` | liver | $1.5$ | l | ✓ |
| $\mathrm{Vmem}$ | `Vmem` | plasma membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vapical}$ | `Vapical` | apical membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vbi}$ | `Vbi` | bile | $1$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{los}_{\mathrm{ext}}$ | `los_ext` | losartan (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{los}$ | `los` | losartan (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{e3174}_{\mathrm{ext}}$ | `e3174_ext` | E3174 (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $e_{3174}$ | `e3174` | E3174 (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{e3174}_{\mathrm{bi}}$ | `e3174_bi` | E3174 (bile) | $\mathrm{Vbi}$ |  | mmol | amount |
| $\mathrm{e3174}_{\mathrm{lumen}}$ | `e3174_lumen` | E3174 (gut lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{l158}_{\mathrm{ext}}$ | `l158_ext` | other losartan metabolites (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $l_{158}$ | `l158` | other losartan metabolites (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{l158}_{\mathrm{bi}}$ | `l158_bi` | other losartan metabolites (bile) | $\mathrm{Vbi}$ |  | mmol | amount |
| $\mathrm{l158}_{\mathrm{lumen}}$ | `l158_lumen` | other losartan metabolites (gut lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{cyp2c9}}$ | `f_cyp2c9` | activity of CYP2C9 | $1$ | \- | ✓ |
| $f_{\mathrm{l158\_cyp2c9}}$ | `f_l158_cyp2c9` | fraction of the formation of the other metabolites by CYP2C9 | $0.5$ | \- | ✓ |
| $\mathrm{LOSIM}_{\mathrm{k}}$ | `LOSIM_k` | rate losartan import | $50$ | 1/min | ✓ |
| $\mathrm{LOS2E3174}_{\mathrm{k}}$ | `LOS2E3174_k` | rate losartan → E3174 | $0.412801474839129$ | 1/min | ✓ |
| $\mathrm{LOS2L158}_{\mathrm{k}}$ | `LOS2L158_k` | rate losartan → other metabolites | $1.78636406676044$ | 1/min | ✓ |
| $\mathrm{E3174EX}_{\mathrm{k}}$ | `E3174EX_k` | rate E3174 export | $0.00411279765463081$ | 1/min | ✓ |
| $\mathrm{E3174IM}_{\mathrm{k}}$ | `E3174IM_k` | rate E3174 import | $0.0614100473485477$ | 1/min | ✓ |
| $\mathrm{L158EX}_{\mathrm{k}}$ | `L158EX_k` | rate export of the other metabolites | $1$ | 1/min | ✓ |
| $\mathrm{E3174BIEX}_{\mathrm{k}}$ | `E3174BIEX_k` | rate E3174 bile export | $0.00274339281823762$ | 1/min | ✓ |
| $\mathrm{L158BIEX}_{\mathrm{k}}$ | `L158BIEX_k` | rate bile export of the other metabolites | $0.257303343770967$ | 1/min | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{e3174}_{\mathrm{bi}} &= 0 \cdot \mathrm{Vbi} \\[1ex]
\mathrm{l158}_{\mathrm{bi}} &= 0 \cdot \mathrm{Vbi}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{LOSIM}}$ | `LOSIM` | losartan import (LOSIM) | $\mathrm{los}_{\mathrm{ext}} \rightleftharpoons \mathrm{los}$ |
| $v_{\mathrm{LOS2E3174}}$ | `LOS2E3174` | E3174 formation (CYP2C9) | $\mathrm{los} \longrightarrow e_{3174}$ |
| $v_{\mathrm{LOS2L158}}$ | `LOS2L158` | other metabolites formation (CYP2C9, CYP3A4, UGT) | $\mathrm{los} \longrightarrow l_{158}$ |
| $v_{\mathrm{E3174EX}}$ | `E3174EX` | E3174 export (E3174EX) | $e_{3174} \rightleftharpoons \mathrm{e3174}_{\mathrm{ext}}$ |
| $v_{\mathrm{L158EX}}$ | `L158EX` | other metabolites export (L158EX) | $l_{158} \rightleftharpoons \mathrm{l158}_{\mathrm{ext}}$ |
| $v_{\mathrm{E3174BIEX}}$ | `E3174BIEX` | E3174 bile export | $e_{3174} \longrightarrow \mathrm{e3174}_{\mathrm{bi}}$ |
| $v_{\mathrm{L158BIEX}}$ | `L158BIEX` | other metabolites bile export | $l_{158} \longrightarrow \mathrm{l158}_{\mathrm{bi}}$ |
| $v_{\mathrm{E3174EHC}}$ | `E3174EHC` | E3174 bile into the gut lumen | $\mathrm{e3174}_{\mathrm{bi}} \longrightarrow \mathrm{e3174}_{\mathrm{lumen}}$ |
| $v_{\mathrm{L158EHC}}$ | `L158EHC` | other losartan metabolites bile into the gut lumen | $\mathrm{l158}_{\mathrm{bi}} \longrightarrow \mathrm{l158}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{LOSIM}} &= \mathrm{LOSIM}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathopen{}\left(\mathrm{los}_{\mathrm{ext}} - \mathrm{los}\right) \\[1ex]
v_{\mathrm{LOS2E3174}} &= f_{\mathrm{cyp2c9}} \cdot \mathrm{LOS2E3174}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{los} \\[1ex]
v_{\mathrm{LOS2L158}} &= \mathopen{}\left(f_{\mathrm{l158\_cyp2c9}} \cdot f_{\mathrm{cyp2c9}} + 1 - f_{\mathrm{l158\_cyp2c9}}\right) \cdot \mathrm{LOS2L158}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{los} \\[1ex]
v_{\mathrm{E3174EX}} &= \mathrm{Vli} \cdot \mathopen{}\left(\mathrm{E3174EX}_{\mathrm{k}} \cdot e_{3174} - \mathrm{E3174IM}_{\mathrm{k}} \cdot \mathrm{e3174}_{\mathrm{ext}}\right) \\[1ex]
v_{\mathrm{L158EX}} &= \mathrm{L158EX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathopen{}\left(l_{158} - \mathrm{l158}_{\mathrm{ext}}\right) \\[1ex]
v_{\mathrm{E3174BIEX}} &= \mathrm{E3174BIEX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot e_{3174} \\[1ex]
v_{\mathrm{L158BIEX}} &= \mathrm{L158BIEX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot l_{158} \\[1ex]
v_{\mathrm{E3174EHC}} &= v_{\mathrm{E3174BIEX}} \\[1ex]
v_{\mathrm{L158EHC}} &= v_{\mathrm{L158BIEX}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{los}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{LOSIM}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{los}}{\mathrm{d} t} &= \frac{v_{\mathrm{LOSIM}} - v_{\mathrm{LOS2E3174}} - v_{\mathrm{LOS2L158}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{e3174}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{E3174EX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} e_{3174}}{\mathrm{d} t} &= \frac{v_{\mathrm{LOS2E3174}} - v_{\mathrm{E3174EX}} - v_{\mathrm{E3174BIEX}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{e3174}_{\mathrm{bi}}}{\mathrm{d} t} &= v_{\mathrm{E3174BIEX}} - v_{\mathrm{E3174EHC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{e3174}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{v_{\mathrm{E3174EHC}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{l158}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{L158EX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} l_{158}}{\mathrm{d} t} &= \frac{v_{\mathrm{LOS2L158}} - v_{\mathrm{L158EX}} - v_{\mathrm{L158BIEX}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{l158}_{\mathrm{bi}}}{\mathrm{d} t} &= v_{\mathrm{L158BIEX}} - v_{\mathrm{L158EHC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{l158}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{v_{\mathrm{L158EHC}}}{\mathrm{Vlumen}}
\end{aligned}
$$
