# Model for RAAS system of blood pressure regulation.

Model `losartan_raas`, SBML Level 3 Version 2, read from losartan\_raas.xml, written by sbmlode 0.3.0.

Model for RAAS system of blood pressure regulation.

version 3

Renin cleaves angiotensinogen to angiotensin I, ACE converts angiotensin I to angiotensin II. Angiotensin II acts on the AT1 receptors, which E3174 blocks (fe\_e3174, the fraction of blocked receptors): the AT1 signal f\_at1 is the angiotensin II relative to its baseline times the fraction of free receptors. The AT1 signal drives the secretion of aldosterone, inhibits the secretion of renin (short loop feedback, RENSEC\_n) and raises the blood pressure (BP\_emax). The turnover of renin, aldosterone and the angiotensins is given by their half-lives in plasma (renin 15 min, aldosterone 20 min, angiotensin I 0.5 min and angiotensin II 0.4 min); the concentrations do not depend on the volume of the compartment. The baselines (\*\_ref) are steady states, the initial concentrations are the baselines (initial assignments), which the experiments set per study.

Changelog

version 3

turnover from the half-lives, angiotensinogen 1 µM

AT1 signal of angiotensin II and E3174: aldosterone secretion, renin feedback and blood pressure

version 2

better handling of baselines and scaling

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
| $\mathrm{Vplasma}$ | `Vplasma` | plasma | $3$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{ald}$ | `ald` | aldosterone | $\mathrm{Vplasma}$ | $2.5 \times 10^{-7}$ | mmol/l | concentration |
| $\mathrm{ren}$ | `ren` | renin | $\mathrm{Vplasma}$ | $10^{-9}$ | mmol/l | concentration |
| $\mathrm{anggen}$ | `anggen` | angiotensinogen | $\mathrm{Vplasma}$ | $0.001$ | mmol/l | concentration, boundary |
| $\mathrm{ang}_{1}$ | `ang1` | angiotensin I | $\mathrm{Vplasma}$ | $10^{-8}$ | mmol/l | concentration |
| $\mathrm{ang}_{2}$ | `ang2` | angiotensin II | $\mathrm{Vplasma}$ | $8 \times 10^{-9}$ | mmol/l | concentration |
| $e_{3174}$ | `e3174` | e3174 | $\mathrm{Vplasma}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{ald}_{\mathrm{change}}$ | `ald_change` | aldosterone change |  | mmol/l |  |
| $\mathrm{ald}_{\mathrm{ratio}}$ | `ald_ratio` | aldosterone ratio |  | \- |  |
| $\mathrm{ren}_{\mathrm{change}}$ | `ren_change` | renin change |  | mmol/l |  |
| $\mathrm{ren}_{\mathrm{ratio}}$ | `ren_ratio` | renin ratio |  | \- |  |
| $\mathrm{ang1}_{\mathrm{change}}$ | `ang1_change` | angiotensin I change |  | mmol/l |  |
| $\mathrm{ang1}_{\mathrm{ratio}}$ | `ang1_ratio` | angiotensin I ratio |  | \- |  |
| $\mathrm{ang2}_{\mathrm{change}}$ | `ang2_change` | angiotensin II change |  | mmol/l |  |
| $\mathrm{ang2}_{\mathrm{ratio}}$ | `ang2_ratio` | angiotensin II ratio |  | \- |  |
| $\mathrm{anggen}_{\mathrm{ref}}$ | `anggen_ref` | reference concentration of angiotensinogen | $0.001$ | mmol/l | ✓ |
| $\mathrm{ang1}_{\mathrm{ref}}$ | `ang1_ref` | reference concentration of angiotensin I | $10^{-8}$ | mmol/l | ✓ |
| $\mathrm{ang2}_{\mathrm{ref}}$ | `ang2_ref` | reference concentration of angiotensin II | $8 \times 10^{-9}$ | mmol/l | ✓ |
| $\mathrm{ren}_{\mathrm{ref}}$ | `ren_ref` | reference concentration of renin | $10^{-9}$ | mmol/l | ✓ |
| $\mathrm{ald}_{\mathrm{ref}}$ | `ald_ref` | reference concentration of aldosterone | $2.5 \times 10^{-7}$ | mmol/l | ✓ |
| $\mathrm{E50}_{\mathrm{e3174}}$ | `E50_e3174` | half-maximum effect E3174 | $0.000206885369335131$ | mmol/l | ✓ |
| $\mathrm{ren}_{\mathrm{t12}}$ | `ren_t12` | half-life of renin | $15$ | min | ✓ |
| $\mathrm{ald}_{\mathrm{t12}}$ | `ald_t12` | half-life of aldosterone | $20$ | min | ✓ |
| $\mathrm{ang1}_{\mathrm{t12}}$ | `ang1_t12` | half-life of angiotensin I | $0.5$ | min | ✓ |
| $\mathrm{RENSEC}_{\mathrm{n}}$ | `RENSEC_n` | feedback of the AT1 signal on the renin secretion | $2.96806352767082$ | \- | ✓ |
| $\mathrm{SBP}_{\mathrm{ref}}$ | `SBP_ref` | reference systolic blood pressure | $120$ | (133.322 N)/m<sup>2</sup> | ✓ |
| $\mathrm{DBP}_{\mathrm{ref}}$ | `DBP_ref` | reference diastolic blood pressure | $80$ | (133.322 N)/m<sup>2</sup> | ✓ |
| $\mathrm{SBP}$ | `SBP` | systolic blood pressure |  | (133.322 N)/m<sup>2</sup> |  |
| $\mathrm{DBP}$ | `DBP` | diastolic blood pressure |  | (133.322 N)/m<sup>2</sup> |  |
| $\mathrm{MAP}$ | `MAP` | mean arterial pressure |  | (133.322 N)/m<sup>2</sup> |  |
| $\mathrm{BP}_{\mathrm{emax}}$ | `BP_emax` | maximal fractional fall of the blood pressure | $0.244346661432525$ | \- | ✓ |
| $\mathrm{SBP}_{\mathrm{change}}$ | `SBP_change` | systolic blood pressure change |  | (133.322 N)/m<sup>2</sup> |  |
| $\mathrm{SBP}_{\mathrm{ratio}}$ | `SBP_ratio` | systolic blood pressure ratio |  | \- |  |
| $\mathrm{DBP}_{\mathrm{change}}$ | `DBP_change` | diastolic blood pressure change |  | (133.322 N)/m<sup>2</sup> |  |
| $\mathrm{DBP}_{\mathrm{ratio}}$ | `DBP_ratio` | diastolic blood pressure ratio |  | \- |  |
| $\mathrm{RENSEC}_{\mathrm{k}}$ | `RENSEC_k` | rate renin secretion |  | mmol/min |  |
| $\mathrm{RENDEG}_{\mathrm{k}}$ | `RENDEG_k` | rate renin degradation |  | l/min |  |
| $\mathrm{ALDSEC}_{\mathrm{k}}$ | `ALDSEC_k` | rate aldosterone secretion |  | mmol/min |  |
| $\mathrm{ALDDEG}_{\mathrm{k}}$ | `ALDDEG_k` | rate aldosterone degradation |  | l/min |  |
| $\mathrm{ANGGEN2ANG1}_{\mathrm{k}}$ | `ANGGEN2ANG1_k` | rate angen → ang1 |  | l/min |  |
| $\mathrm{ANG1ANG2}_{\mathrm{k}}$ | `ANG1ANG2_k` | rate ang1 → ang2 |  | l/min |  |
| $\mathrm{ANG2DEG}_{\mathrm{k}}$ | `ANG2DEG_k` | rate angiotensin II degradation |  | l/min |  |
| $\mathrm{fe}_{\mathrm{e3174}}$ | `fe_e3174` | AT1 receptor blockade by E3174 in \[0, 1\] |  | \- |  |
| $f_{\mathrm{at1}}$ | `f_at1` | AT1 signal relative to the baseline |  | \- |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{ald} &= \mathrm{ald}_{\mathrm{ref}} \\[1ex]
\mathrm{ren} &= \mathrm{ren}_{\mathrm{ref}} \\[1ex]
\mathrm{anggen} &= \mathrm{anggen}_{\mathrm{ref}} \\[1ex]
\mathrm{ang}_{1} &= \mathrm{ang1}_{\mathrm{ref}} \\[1ex]
\mathrm{ang}_{2} &= \mathrm{ang2}_{\mathrm{ref}}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{RENDEG}_{\mathrm{k}} &= \frac{\ln\mathopen{}\left(2\right)}{\mathrm{ren}_{\mathrm{t12}}} \cdot \mathrm{Vplasma} \\[1ex]
\mathrm{RENSEC}_{\mathrm{k}} &= \mathrm{RENDEG}_{\mathrm{k}} \cdot \mathrm{ren}_{\mathrm{ref}} \\[1ex]
\mathrm{ALDDEG}_{\mathrm{k}} &= \frac{\ln\mathopen{}\left(2\right)}{\mathrm{ald}_{\mathrm{t12}}} \cdot \mathrm{Vplasma} \\[1ex]
\mathrm{ALDSEC}_{\mathrm{k}} &= \mathrm{ALDDEG}_{\mathrm{k}} \cdot \mathrm{ald}_{\mathrm{ref}} \\[1ex]
\mathrm{ANG1ANG2}_{\mathrm{k}} &= \frac{\ln\mathopen{}\left(2\right)}{\mathrm{ang1}_{\mathrm{t12}}} \cdot \mathrm{Vplasma} \\[1ex]
\mathrm{ANGGEN2ANG1}_{\mathrm{k}} &= \frac{\mathrm{ANG1ANG2}_{\mathrm{k}} \cdot \mathrm{ang1}_{\mathrm{ref}}}{\mathrm{anggen}_{\mathrm{ref}}} \\[1ex]
\mathrm{ANG2DEG}_{\mathrm{k}} &= \frac{\mathrm{ANG1ANG2}_{\mathrm{k}} \cdot \mathrm{ang1}_{\mathrm{ref}}}{\mathrm{ang2}_{\mathrm{ref}}} \\[1ex]
\mathrm{ald}_{\mathrm{change}} &= \mathrm{ald} - \mathrm{ald}_{\mathrm{ref}} \\[1ex]
\mathrm{ald}_{\mathrm{ratio}} &= \frac{\mathrm{ald}}{\mathrm{ald}_{\mathrm{ref}}} \\[1ex]
\mathrm{ren}_{\mathrm{change}} &= \mathrm{ren} - \mathrm{ren}_{\mathrm{ref}} \\[1ex]
\mathrm{ren}_{\mathrm{ratio}} &= \frac{\mathrm{ren}}{\mathrm{ren}_{\mathrm{ref}}} \\[1ex]
\mathrm{ang1}_{\mathrm{change}} &= \mathrm{ang}_{1} - \mathrm{ang1}_{\mathrm{ref}} \\[1ex]
\mathrm{ang1}_{\mathrm{ratio}} &= \frac{\mathrm{ang}_{1}}{\mathrm{ang1}_{\mathrm{ref}}} \\[1ex]
\mathrm{ang2}_{\mathrm{change}} &= \mathrm{ang}_{2} - \mathrm{ang2}_{\mathrm{ref}} \\[1ex]
\mathrm{ang2}_{\mathrm{ratio}} &= \frac{\mathrm{ang}_{2}}{\mathrm{ang2}_{\mathrm{ref}}} \\[1ex]
\mathrm{fe}_{\mathrm{e3174}} &= \frac{e_{3174}}{\mathrm{E50}_{\mathrm{e3174}} + e_{3174}} \\[1ex]
f_{\mathrm{at1}} &= \frac{\mathrm{ang}_{2}}{\mathrm{ang2}_{\mathrm{ref}}} \cdot \mathopen{}\left(1 - \mathrm{fe}_{\mathrm{e3174}}\right) \\[1ex]
\mathrm{SBP} &= \mathrm{SBP}_{\mathrm{ref}} \cdot \mathopen{}\left(1 - \mathrm{BP}_{\mathrm{emax}} \cdot \mathopen{}\left(1 - f_{\mathrm{at1}}\right)\right) \\[1ex]
\mathrm{DBP} &= \mathrm{DBP}_{\mathrm{ref}} \cdot \mathopen{}\left(1 - \mathrm{BP}_{\mathrm{emax}} \cdot \mathopen{}\left(1 - f_{\mathrm{at1}}\right)\right) \\[1ex]
\mathrm{MAP} &= \mathrm{DBP} + \frac{\mathrm{SBP} - \mathrm{DBP}}{3} \\[1ex]
\mathrm{SBP}_{\mathrm{change}} &= \mathrm{SBP} - \mathrm{SBP}_{\mathrm{ref}} \\[1ex]
\mathrm{SBP}_{\mathrm{ratio}} &= \frac{\mathrm{SBP}}{\mathrm{SBP}_{\mathrm{ref}}} \\[1ex]
\mathrm{DBP}_{\mathrm{change}} &= \mathrm{DBP} - \mathrm{DBP}_{\mathrm{ref}} \\[1ex]
\mathrm{DBP}_{\mathrm{ratio}} &= \frac{\mathrm{DBP}}{\mathrm{DBP}_{\mathrm{ref}}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation | Modifiers |
| --- | --- | --- | --- | --- |
| $v_{\mathrm{RENSEC}}$ | `RENSEC` | renin secretion (RENSEC) | $\varnothing \longrightarrow \mathrm{ren}$ | $\mathrm{ang}_{2}, e_{3174}$ |
| $v_{\mathrm{RENDEG}}$ | `RENDEG` | renin degradation (RENDEG) | $\mathrm{ren} \longrightarrow \varnothing$ |  |
| $v_{\mathrm{ALDSEC}}$ | `ALDSEC` | aldosterone secretion (ALDSEC) | $\varnothing \longrightarrow \mathrm{ald}$ | $\mathrm{ang}_{2}, e_{3174}$ |
| $v_{\mathrm{ALDDEG}}$ | `ALDDEG` | aldosterone degradation (ALDDEG) | $\mathrm{ald} \longrightarrow \varnothing$ |  |
| $v_{\mathrm{ANGGEN2ANG1}}$ | `ANGGEN2ANG1` | angiotensinogen to angiotensin I (renin) | $\mathrm{anggen} \longrightarrow \mathrm{ang}_{1}$ | $\mathrm{ren}$ |
| $v_{\mathrm{ANG1ANG2}}$ | `ANG1ANG2` | angiotensin I to angiotensin II (ACE) | $\mathrm{ang}_{1} \longrightarrow \mathrm{ang}_{2}$ |  |
| $v_{\mathrm{ANG2DEG}}$ | `ANG2DEG` | angiotensin II degradation (ANG2DEG) | $\mathrm{ang}_{2} \longrightarrow \varnothing$ |  |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RENSEC}} &= \mathrm{RENSEC}_{\mathrm{k}} \cdot f_{\mathrm{at1}}^{-\mathrm{RENSEC}_{\mathrm{n}}} \\[1ex]
v_{\mathrm{RENDEG}} &= \mathrm{RENDEG}_{\mathrm{k}} \cdot \mathrm{ren} \\[1ex]
v_{\mathrm{ALDSEC}} &= \mathrm{ALDSEC}_{\mathrm{k}} \cdot f_{\mathrm{at1}} \\[1ex]
v_{\mathrm{ALDDEG}} &= \mathrm{ALDDEG}_{\mathrm{k}} \cdot \mathrm{ald} \\[1ex]
v_{\mathrm{ANGGEN2ANG1}} &= \frac{\mathrm{ANGGEN2ANG1}_{\mathrm{k}} \cdot \mathrm{anggen} \cdot \mathrm{ren}}{\mathrm{ren}_{\mathrm{ref}}} \\[1ex]
v_{\mathrm{ANG1ANG2}} &= \mathrm{ANG1ANG2}_{\mathrm{k}} \cdot \mathrm{ang}_{1} \\[1ex]
v_{\mathrm{ANG2DEG}} &= \mathrm{ANG2DEG}_{\mathrm{k}} \cdot \mathrm{ang}_{2}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{ald}}{\mathrm{d} t} &= \frac{v_{\mathrm{ALDSEC}} - v_{\mathrm{ALDDEG}}}{\mathrm{Vplasma}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{ren}}{\mathrm{d} t} &= \frac{v_{\mathrm{RENSEC}} - v_{\mathrm{RENDEG}}}{\mathrm{Vplasma}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{ang}_{1}}{\mathrm{d} t} &= \frac{v_{\mathrm{ANGGEN2ANG1}} - v_{\mathrm{ANG1ANG2}}}{\mathrm{Vplasma}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{ang}_{2}}{\mathrm{d} t} &= \frac{v_{\mathrm{ANG1ANG2}} - v_{\mathrm{ANG2DEG}}}{\mathrm{Vplasma}}
\end{aligned}
$$
