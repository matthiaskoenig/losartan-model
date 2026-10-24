# Model for losartan absorption in the small intestine

Model `losartan_intestine`, SBML Level 3 Version 2, read from losartan\_intestine.xml, written by sbmlode 0.3.0.

Model for losartan absorption

dissolution of the tablet and transit into the gut lumen

absorption of losartan (LOSABS, OATP2B1) into the enterocytes, efflux back into the lumen (LOSEFL, P-glycoprotein/ABCB1) and transport into the portal blood (LOSEX); the absorption is nearly complete, the bioavailability of 33 % is the hepatic first pass \[Lo1995\]

the fecal excretion of the content of the gut lumen with the transit through the colon (METEXC\_k): unabsorbed losartan and the E3174 and other metabolites of the bile

The species of the gut lumen (\*\_lumen) are concentrations in the lumen volume Vlumen.

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
| $\mathrm{Vext}$ | `Vext` | blood | $1$ | l | ✓ |
| $\mathrm{Vgu}$ | `Vgu` | intestine | $1.2825$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Ventero}$ | `Ventero` | intestinal lining (enterocytes) | $1$ | l |  |
| $\mathrm{Vapical}$ | `Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vbaso}$ | `Vbaso` | basolateral membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vstomach}$ | `Vstomach` | stomach | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{los}_{\mathrm{stomach}}$ | `los_stomach` | losartan (stomach) | $\mathrm{Vstomach}$ |  | mmol | amount, boundary |
| $\mathrm{los}_{\mathrm{lumen}}$ | `los_lumen` | losartan (gut lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{los}$ | `los` | losartan (enterocytes) | $\mathrm{Ventero}$ | $0$ | mmol/l | concentration |
| $\mathrm{los}_{\mathrm{ext}}$ | `los_ext` | losartan (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{los}_{\mathrm{feces}}$ | `los_feces` | losartan (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{e3174}_{\mathrm{lumen}}$ | `e3174_lumen` | E3174 (gut lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{e3174}_{\mathrm{feces}}$ | `e3174_feces` | E3174 (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{l158}_{\mathrm{lumen}}$ | `l158_lumen` | other losartan metabolites (gut lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{l158}_{\mathrm{feces}}$ | `l158_feces` | other losartan metabolites (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{LOSABS}_{\mathrm{k}}$ | `LOSABS_k` | rate of losartan absorption | $0.0476385224387545$ | 1/min | ✓ |
| $\mathrm{METEXC}_{\mathrm{k}}$ | `METEXC_k` | rate of feces excretion | $0.0006$ | 1/min | ✓ |
| $\mathrm{PODOSE}_{\mathrm{los}}$ | `PODOSE_los` | oral dose losartan \[mg\] | $0$ | mg |  |
| $\mathrm{Ka}_{\mathrm{dis\_los}}$ | `Ka_dis_los` | rate dissolution losartan | $2$ | 1/hr | ✓ |
| $\mathrm{Mr}_{\mathrm{los}}$ | `Mr_los` | Molecular weight losartan \[g/mole\] | $422.911$ | g/mol | ✓ |
| $f_{\mathrm{OATP2B1}}$ | `f_OATP2B1` | absorption activity | $1$ | \- | ✓ |
| $f_{\mathrm{abcb1}}$ | `f_abcb1` | PG activity | $1$ | \- | ✓ |
| $f_{\mathrm{LOSEFL\_k}}$ | `f_LOSEFL_k` | fractional rate losartan efflux (PG) | $1$ | \- | ✓ |
| $\mathrm{LOSEFL}_{\mathrm{k}}$ | `LOSEFL_k` |  |  | 1/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{los}_{\mathrm{stomach}} &= 0 \cdot \mathrm{Vstomach} \\[1ex]
\mathrm{los}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{e3174}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{l158}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{LOSEFL}_{\mathrm{k}} &= f_{\mathrm{LOSEFL\_k}} \cdot \mathrm{LOSABS}_{\mathrm{k}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{LOSABS}}$ | `LOSABS` | absorption losartan | $\mathrm{los}_{\mathrm{lumen}} \longrightarrow \mathrm{los}$ |
| $v_{\mathrm{LOSEFL}}$ | `LOSEFL` | efflux losartan (PG) | $\mathrm{los} \longrightarrow \mathrm{los}_{\mathrm{lumen}}$ |
| $v_{\mathrm{LOSEX}}$ | `LOSEX` | losartan export plasma | $\mathrm{los} \longrightarrow \mathrm{los}_{\mathrm{ext}}$ |
| $v_{\mathrm{LOSEXC}}$ | `LOSEXC` | excretion losartan (feces) | $\mathrm{los}_{\mathrm{lumen}} \longrightarrow \mathrm{los}_{\mathrm{feces}}$ |
| $v_{\mathrm{E3174EXC}}$ | `E3174EXC` | excretion E3174 (feces) | $\mathrm{e3174}_{\mathrm{lumen}} \longrightarrow \mathrm{e3174}_{\mathrm{feces}}$ |
| $v_{\mathrm{L158EXC}}$ | `L158EXC` | excretion other metabolites (feces) | $\mathrm{l158}_{\mathrm{lumen}} \longrightarrow \mathrm{l158}_{\mathrm{feces}}$ |
| $v_{\mathrm{dissolution\_los}}$ | `dissolution_los` | dissolution losartan | $\mathrm{los}_{\mathrm{stomach}} \longrightarrow \mathrm{los}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{LOSABS}} &= f_{\mathrm{OATP2B1}} \cdot \mathrm{LOSABS}_{\mathrm{k}} \cdot \mathrm{Vgu} \cdot \mathrm{los}_{\mathrm{lumen}} \\[1ex]
v_{\mathrm{LOSEFL}} &= f_{\mathrm{abcb1}} \cdot \mathrm{LOSEFL}_{\mathrm{k}} \cdot \mathrm{Vgu} \cdot \mathrm{los} \\[1ex]
v_{\mathrm{LOSEX}} &= \mathrm{LOSABS}_{\mathrm{k}} \cdot \mathrm{Vgu} \cdot \mathrm{los} \\[1ex]
v_{\mathrm{LOSEXC}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vgu} \cdot \mathrm{los}_{\mathrm{lumen}} \\[1ex]
v_{\mathrm{E3174EXC}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vgu} \cdot \mathrm{e3174}_{\mathrm{lumen}} \\[1ex]
v_{\mathrm{L158EXC}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vgu} \cdot \mathrm{l158}_{\mathrm{lumen}} \\[1ex]
v_{\mathrm{dissolution\_los}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{dis\_los}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{los}}}{\mathrm{Mr}_{\mathrm{los}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{los}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{LOSABS}} + v_{\mathrm{LOSEFL}} - v_{\mathrm{LOSEXC}} + v_{\mathrm{dissolution\_los}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{los}}{\mathrm{d} t} &= \frac{v_{\mathrm{LOSABS}} - v_{\mathrm{LOSEFL}} - v_{\mathrm{LOSEX}}}{\mathrm{Ventero}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{los}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LOSEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{los}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{LOSEXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{e3174}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{E3174EXC}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{e3174}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{E3174EXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{l158}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{L158EXC}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{l158}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{L158EXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{los}}}{\mathrm{d} t} &= -v_{\mathrm{dissolution\_los}} \cdot \mathrm{Mr}_{\mathrm{los}} \qquad \text{(rate rule)}
\end{aligned}
$$
