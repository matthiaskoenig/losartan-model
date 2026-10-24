# Model for renal excretion of losartan, E3174, and the other metabolites.

Model `losartan_kidney`, SBML Level 3 Version 2, read from losartan\_kidney.xml, written by sbmlode 0.3.0.

Model for renal excretion of losartan, E3174, and the other metabolites.

version 2

The renal clearances are given relative to the glomerular filtration rate (GFR) and scale with the renal function (f\_renal\_function): the renal clearance of plasma of losartan and E3174 is &lt;X&gt;EX\_f times the filtration of their unbound fraction in plasma (fu\_&lt;x&gt; \* GFR), the one of the other metabolites L158EX\_f times the GFR. The concentrations of the kidney model are blood concentrations, the plasma concentration is the blood concentration over the blood to plasma ratio BP\_&lt;x&gt;.

Changelog

version 2

renal clearances relative to the GFR, plasma by the blood to plasma ratio

version 1

initial model

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
| $\mathrm{Vki}$ | `Vki` | kidney | $0.3$ | l | ✓ |
| $\mathrm{Vmem}$ | `Vmem` | plasma membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{los}_{\mathrm{ext}}$ | `los_ext` | losartan (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{los}_{\mathrm{urine}}$ | `los_urine` | losartan (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{e3174}_{\mathrm{ext}}$ | `e3174_ext` | E3174 (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{e3174}_{\mathrm{urine}}$ | `e3174_urine` | E3174 (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{l158}_{\mathrm{ext}}$ | `l158_ext` | other losartan metabolites (blood) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{l158}_{\mathrm{urine}}$ | `l158_urine` | other losartan metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{renal\_function}}$ | `f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{GFR}_{\mathrm{base}}$ | `GFR_base` | glomerular filtration rate (healthy) | $101$ | ml/min | ✓ |
| $\mathrm{fu}_{\mathrm{los}}$ | `fu_los` | fraction unbound in plasma losartan | $0.013$ | \- | ✓ |
| $\mathrm{fu}_{\mathrm{e3174}}$ | `fu_e3174` | fraction unbound in plasma E3174 | $0.002$ | \- | ✓ |
| $\mathrm{BP}_{\mathrm{los}}$ | `BP_los` | blood to plasma ratio losartan | $0.53$ | \- | ✓ |
| $\mathrm{BP}_{\mathrm{e3174}}$ | `BP_e3174` | blood to plasma ratio E3174 | $0.58$ | \- | ✓ |
| $\mathrm{BP}_{\mathrm{l158}}$ | `BP_l158` | blood to plasma ratio other metabolites | $1$ | \- | ✓ |
| $\mathrm{LOSEX}_{\mathrm{f}}$ | `LOSEX_f` | renal clearance of unbound losartan relative to the GFR | $59.6320404780412$ | \- | ✓ |
| $\mathrm{E3174EX}_{\mathrm{f}}$ | `E3174EX_f` | renal clearance of unbound E3174 relative to the GFR | $124.916497816236$ | \- | ✓ |
| $\mathrm{L158EX}_{\mathrm{f}}$ | `L158EX_f` | renal clearance of the other metabolites relative to the GFR | $5.18711908659398$ | \- | ✓ |
| $\mathrm{GFR}$ | `GFR` | glomerular filtration rate |  | ml/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{los}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{e3174}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{l158}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{GFR} &= f_{\mathrm{renal\_function}} \cdot \mathrm{GFR}_{\mathrm{base}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{LOSEX}}$ | `LOSEX` | losartan excretion (LOSEX) | $\mathrm{los}_{\mathrm{ext}} \longrightarrow \mathrm{los}_{\mathrm{urine}}$ |
| $v_{\mathrm{E3174EX}}$ | `E3174EX` | E3174 excretion (E3174EX) | $\mathrm{e3174}_{\mathrm{ext}} \longrightarrow \mathrm{e3174}_{\mathrm{urine}}$ |
| $v_{\mathrm{L158EX}}$ | `L158EX` | other metabolites excretion (L158EX) | $\mathrm{l158}_{\mathrm{ext}} \longrightarrow \mathrm{l158}_{\mathrm{urine}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{LOSEX}} &= \frac{\mathrm{GFR} \cdot 0.001 \cdot \mathrm{fu}_{\mathrm{los}} \cdot \mathrm{LOSEX}_{\mathrm{f}} \cdot \mathrm{los}_{\mathrm{ext}}}{\mathrm{BP}_{\mathrm{los}}} \\[1ex]
v_{\mathrm{E3174EX}} &= \frac{\mathrm{GFR} \cdot 0.001 \cdot \mathrm{fu}_{\mathrm{e3174}} \cdot \mathrm{E3174EX}_{\mathrm{f}} \cdot \mathrm{e3174}_{\mathrm{ext}}}{\mathrm{BP}_{\mathrm{e3174}}} \\[1ex]
v_{\mathrm{L158EX}} &= \frac{\mathrm{GFR} \cdot 0.001 \cdot \mathrm{L158EX}_{\mathrm{f}} \cdot \mathrm{l158}_{\mathrm{ext}}}{\mathrm{BP}_{\mathrm{l158}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{los}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{LOSEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{los}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{LOSEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{e3174}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{E3174EX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{e3174}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{E3174EX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{l158}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{L158EX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{l158}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{L158EX}}
\end{aligned}
$$
