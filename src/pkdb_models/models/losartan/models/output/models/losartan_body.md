# Losartan body model

Model `losartan_body`, SBML Level 3 Version 2, read from losartan\_body.xml, written by sbmlode 0.3.0.

Losartan body model

Whole-body model of losartan, its active metabolite E3174 and its other (inactive) metabolites with the blood circulation, the tissue distribution, the submodels of the intestine, the liver and the kidneys, and the renin-angiotensin-aldosterone system (RAAS) on which E3174 acts.

The circulating species are whole blood concentrations: the blood flows transport the blood concentrations, the plasma concentration is the blood concentration over the blood to plasma ratio BP\_&lt;sid&gt; (losartan 0.53, E3174 0.58 \[Lo1995\], the other metabolites 1). The tissues exchange with the blood by the flow ftissue\_&lt;sid&gt; towards the tissue/blood partition coefficient Kp\_&lt;sid&gt;.

Literature values of the model \[Lo1995, FDA1995S60\]: plasma clearance of losartan 610 ml/min (renal 70 ml/min, 12 %), Vss 34 l, terminal half-life 2.1 hr, bioavailability 33 % by hepatic first pass with nearly complete absorption, 14 % of a dose is converted to E3174; E3174 plasma clearance 47 ml/min (renal 26 ml/min, 55 %), Vss 10 l, half-life 6-9 hr. After 14C losartan 43 % (IV) and 35 % (PO) of the radioactivity is recovered in the urine (losartan 11.4 % and 4.5 %, E3174 6.4 %), 50 % and 58 % in the feces.

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
| $\mathrm{Vre}$ | `Vre` | rest of body |  | l |  |
| $\mathrm{Vgu}$ | `Vgu` | gut |  | l |  |
| $\mathrm{Vki}$ | `Vki` | kidney |  | l |  |
| $\mathrm{Vli}$ | `Vli` | liver |  | l |  |
| $\mathrm{Vlu}$ | `Vlu` | lung |  | l |  |
| $\mathrm{Vve}$ | `Vve` | venous blood |  | l |  |
| $\mathrm{Var}$ | `Var` | arterial blood |  | l |  |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vstomach}$ | `Vstomach` | stomach | $1$ | l | ✓ |
| $\mathrm{Vpo}$ | `Vpo` | portal blood |  | l |  |
| $\mathrm{Vhv}$ | `Vhv` | hepatic venous blood |  | l |  |
| $\mathrm{Vki}_{\mathrm{tissue}}$ | `Vki_tissue` | ki tissue |  | l |  |
| $\mathrm{Vki}_{\mathrm{blood}}$ | `Vki_blood` | ki blood |  | l |  |
| $\mathrm{Vli}_{\mathrm{tissue}}$ | `Vli_tissue` | li tissue |  | l |  |
| $\mathrm{Vli}_{\mathrm{blood}}$ | `Vli_blood` | li blood |  | l |  |
| $\mathrm{Vlu}_{\mathrm{tissue}}$ | `Vlu_tissue` | lu tissue |  | l |  |
| $\mathrm{Vlu}_{\mathrm{blood}}$ | `Vlu_blood` | lu blood |  | l |  |
| $\mathrm{Vgu}_{\mathrm{tissue}}$ | `Vgu_tissue` | gu tissue |  | l |  |
| $\mathrm{Vgu}_{\mathrm{blood}}$ | `Vgu_blood` | gu blood |  | l |  |
| $\mathrm{Vre}_{\mathrm{tissue}}$ | `Vre_tissue` | re tissue |  | l |  |
| $\mathrm{Vre}_{\mathrm{blood}}$ | `Vre_blood` | re blood |  | l |  |
| $\mathrm{KI}_{\mathrm{\_Vmem}}$ | `KI__Vmem` | plasma membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{LI}_{\mathrm{\_Vmem}}$ | `LI__Vmem` | plasma membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{LI}_{\mathrm{\_Vapical}}$ | `LI__Vapical` | apical membrane | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{LI}_{\mathrm{\_Vbi}}$ | `LI__Vbi` | bile | $1$ | l | ✓ |
| $\mathrm{LI}_{\mathrm{\_Vlumen}}$ | `LI__Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |
| $\mathrm{GU}_{\mathrm{\_Vlumen}}$ | `GU__Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |
| $\mathrm{GU}_{\mathrm{\_Ventero}}$ | `GU__Ventero` | intestinal lining (enterocytes) | $1$ | l |  |
| $\mathrm{GU}_{\mathrm{\_Vapical}}$ | `GU__Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vbaso}}$ | `GU__Vbaso` | basolateral membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{ald}$ | `ald` | aldosterone | $\mathrm{Vplasma}$ | $2.5 \times 10^{-7}$ | mmol/l | concentration |
| $\mathrm{ren}$ | `ren` | renin | $\mathrm{Vplasma}$ | $10^{-9}$ | mmol/l | concentration |
| $\mathrm{anggen}$ | `anggen` | angiotensinogen | $\mathrm{Vplasma}$ | $0.001$ | mmol/l | concentration, boundary |
| $\mathrm{ang}_{1}$ | `ang1` | angiotensin I | $\mathrm{Vplasma}$ | $10^{-8}$ | mmol/l | concentration |
| $\mathrm{ang}_{2}$ | `ang2` | angiotensin II | $\mathrm{Vplasma}$ | $8 \times 10^{-9}$ | mmol/l | concentration |
| $e_{3174}$ | `e3174` | e3174 | $\mathrm{Vplasma}$ |  | mmol/l | concentration |
| $\mathrm{Cki}_{\mathrm{blood\_los}}$ | `Cki_blood_los` | losartan (kidney blood) | $\mathrm{Vki}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{blood\_los}}$ | `Cli_blood_los` | losartan (liver blood) | $\mathrm{Vli}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{blood\_los}}$ | `Clu_blood_los` | losartan (lung blood) | $\mathrm{Vlu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{blood\_los}}$ | `Cgu_blood_los` | losartan (gut blood) | $\mathrm{Vgu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{blood\_los}}$ | `Cre_blood_los` | losartan (rest blood) | $\mathrm{Vre}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{los}}$ | `Car_los` | losartan (arterial blood blood) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{los}}$ | `Cve_los` | losartan (venous blood blood) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{los}}$ | `Cpo_los` | losartan (portal vein blood) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{los}}$ | `Chv_los` | losartan (hepatic vein blood) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{los}}$ | `Clu_los` | losartan (lung) | $\mathrm{Vlu}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{los}}$ | `Cre_los` | losartan (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{los}}$ | `Aurine_los` | losartan (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{los}}$ | `Afeces_los` | losartan (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cki}_{\mathrm{blood\_e3174}}$ | `Cki_blood_e3174` | E3174 (kidney blood) | $\mathrm{Vki}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{blood\_e3174}}$ | `Cli_blood_e3174` | E3174 (liver blood) | $\mathrm{Vli}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{blood\_e3174}}$ | `Clu_blood_e3174` | E3174 (lung blood) | $\mathrm{Vlu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{blood\_e3174}}$ | `Cgu_blood_e3174` | E3174 (gut blood) | $\mathrm{Vgu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{blood\_e3174}}$ | `Cre_blood_e3174` | E3174 (rest blood) | $\mathrm{Vre}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{e3174}}$ | `Car_e3174` | E3174 (arterial blood blood) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{e3174}}$ | `Cve_e3174` | E3174 (venous blood blood) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{e3174}}$ | `Cpo_e3174` | E3174 (portal vein blood) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{e3174}}$ | `Chv_e3174` | E3174 (hepatic vein blood) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{e3174}}$ | `Clu_e3174` | E3174 (lung) | $\mathrm{Vlu}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{e3174}}$ | `Cre_e3174` | E3174 (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{e3174}}$ | `Aurine_e3174` | E3174 (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{e3174}}$ | `Afeces_e3174` | E3174 (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cki}_{\mathrm{blood\_l158}}$ | `Cki_blood_l158` | other losartan metabolites (kidney blood) | $\mathrm{Vki}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{blood\_l158}}$ | `Cli_blood_l158` | other losartan metabolites (liver blood) | $\mathrm{Vli}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{blood\_l158}}$ | `Clu_blood_l158` | other losartan metabolites (lung blood) | $\mathrm{Vlu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{blood\_l158}}$ | `Cgu_blood_l158` | other losartan metabolites (gut blood) | $\mathrm{Vgu}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{blood\_l158}}$ | `Cre_blood_l158` | other losartan metabolites (rest blood) | $\mathrm{Vre}_{\mathrm{blood}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{l158}}$ | `Car_l158` | other losartan metabolites (arterial blood blood) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{l158}}$ | `Cve_l158` | other losartan metabolites (venous blood blood) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{l158}}$ | `Cpo_l158` | other losartan metabolites (portal vein blood) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{l158}}$ | `Chv_l158` | other losartan metabolites (hepatic vein blood) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{l158}}$ | `Aurine_l158` | other losartan metabolites (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{l158}}$ | `Afeces_l158` | other losartan metabolites (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cgu}_{\mathrm{los}}$ | `Cgu_los` | losartan (gut lumen) | $\mathrm{Vgu}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{e3174}}$ | `Cgu_e3174` | E3174 (gut lumen) | $\mathrm{Vgu}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{l158}}$ | `Cgu_l158` | other losartan metabolites (gut lumen) | $\mathrm{Vgu}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cveplasma}_{\mathrm{los}}$ | `Cveplasma_los` | losartan (venous plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Cveplasma}_{\mathrm{e3174}}$ | `Cveplasma_e3174` | E3174 (venous plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Cveplasma}_{\mathrm{l158}}$ | `Cveplasma_l158` | other losartan metabolites (venous plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Cveplasma}_{\mathrm{los\_e3174}}$ | `Cveplasma_los_e3174` | losartan and E3174 (venous plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Cveplasma}_{\mathrm{total}}$ | `Cveplasma_total` | losartan and its metabolites (venous plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_los}}$ | `LI__los` | losartan (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_e3174}}$ | `LI__e3174` | E3174 (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_e3174\_bi}}$ | `LI__e3174_bi` | E3174 (bile) | $\mathrm{LI}_{\mathrm{\_Vbi}}$ |  | mmol | amount |
| $\mathrm{LI}_{\mathrm{\_l158}}$ | `LI__l158` | other losartan metabolites (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_l158\_bi}}$ | `LI__l158_bi` | other losartan metabolites (bile) | $\mathrm{LI}_{\mathrm{\_Vbi}}$ |  | mmol | amount |
| $\mathrm{GU}_{\mathrm{\_los\_stomach}}$ | `GU__los_stomach` | losartan (stomach) | $\mathrm{Vstomach}$ |  | mmol | amount, boundary |
| $\mathrm{GU}_{\mathrm{\_los}}$ | `GU__los` | losartan (enterocytes) | $\mathrm{GU}_{\mathrm{\_Ventero}}$ | $0$ | mmol/l | concentration |

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
| $\mathrm{BW}$ | `BW` | body weight | $75$ | kg | ✓ |
| $\mathrm{BW}_{\mathrm{ref}}$ | `BW_ref` | reference body weight | $75$ | kg | ✓ |
| $\mathrm{COBW}$ | `COBW` | cardiac output per bodyweight | $1.548$ | ml/s/kg | ✓ |
| $f_{\mathrm{cardiac\_function}}$ | `f_cardiac_function` | cardiac function | $1$ | \- |  |
| $\mathrm{CO}$ | `CO` | cardiac output |  | ml/s |  |
| $\mathrm{QC}$ | `QC` | cardiac output |  | l/min |  |
| $\mathrm{Fblood}$ | `Fblood` | blood fraction of organ volume | $0.02$ | \- |  |
| $f_{\mathrm{GFR\_Qki}}$ | `f_GFR_Qki` | relationship GFR and Qki | $76.3105006271061$ | ml/l | ✓ |
| $\mathrm{GFR}_{\mathrm{base}}$ | `GFR_base` | glomerular filtration rate (healthy) |  | ml/min |  |
| $\mathrm{FVgu}$ | `FVgu` | gut fractional tissue volume | $0.0171$ | l/kg | ✓ |
| $\mathrm{FVki}$ | `FVki` | kidney fractional tissue volume | $0.0044$ | l/kg | ✓ |
| $\mathrm{FVli}$ | `FVli` | liver fractional tissue volume | $0.021$ | l/kg | ✓ |
| $\mathrm{FVlu}$ | `FVlu` | lung fractional tissue volume | $0.0076$ | l/kg | ✓ |
| $\mathrm{FVre}$ | `FVre` | rest of body fractional tissue volume |  | l/kg |  |
| $\mathrm{FVve}$ | `FVve` | venous fractional tissue volume | $0.0514$ | l/kg | ✓ |
| $\mathrm{FVar}$ | `FVar` | arterial fractional tissue volume | $0.0257$ | l/kg | ✓ |
| $\mathrm{FVpo}$ | `FVpo` | portal fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FVhv}$ | `FVhv` | hepatic venous fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FQgu}$ | `FQgu` | gut fractional tissue blood flow | $0.18$ | \- | ✓ |
| $\mathrm{FQki}$ | `FQki` | kidney fractional tissue blood flow | $0.19$ | \- | ✓ |
| $\mathrm{FQh}$ | `FQh` | hepatic fractional tissue blood flow | $0.215$ | \- | ✓ |
| $\mathrm{FQlu}$ | `FQlu` | lung fractional tissue blood flow | $1$ | \- | ✓ |
| $\mathrm{FQre}$ | `FQre` | rest of body fractional tissue blood flow |  | \- |  |
| $\mathrm{conversion}_{\mathrm{min\_per\_day}}$ | `conversion_min_per_day` | Conversion factor day to min | $1440$ | min/day | ✓ |
| $f_{\mathrm{cirrhosis}}$ | `f_cirrhosis` | severity of cirrhosis \[0, 0.95\] | $0$ | \- | ✓ |
| $f_{\mathrm{shunts}}$ | `f_shunts` | fraction of portal venous blood shunted by the liver |  | \- |  |
| $f_{\mathrm{tissue\_loss}}$ | `f_tissue_loss` | fraction of lost parenchymal liver volume |  | \- |  |
| $\mathrm{PODOSE}_{\mathrm{los}}$ | `PODOSE_los` | oral dose los \[mg\] | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{los}}$ | `Mr_los` | Molecular weight los | $422.911$ | g/mol | ✓ |
| $\mathrm{BP}_{\mathrm{los}}$ | `BP_los` | blood to plasma ratio los | $0.53$ | \- | ✓ |
| $\mathrm{ftissue}_{\mathrm{los}}$ | `ftissue_los` | tissue distribution LOS | $1.03525478749496$ | l/min | ✓ |
| $\mathrm{Kp}_{\mathrm{los}}$ | `Kp_los` | tissue/blood partition coefficient LOS | $1.16735657203939$ | \- | ✓ |
| $\mathrm{IVDOSE}_{\mathrm{los}}$ | `IVDOSE_los` | IV bolus dose los \[mg\] | $0$ | mg |  |
| $\mathrm{ti}_{\mathrm{los}}$ | `ti_los` | injection time los | $10$ | s | ✓ |
| $\mathrm{Ki}_{\mathrm{los}}$ | `Ki_los` | Ki \[1/min\] injection los |  | 1/min |  |
| $\mathrm{Ri}_{\mathrm{los}}$ | `Ri_los` | Ri \[mg/min\] rate of infusion los | $0$ | mg/min | ✓ |
| $\mathrm{cum}_{\mathrm{dose\_los}}$ | `cum_dose_los` | Cumulative dose due to infusion los | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{e3174}}$ | `Mr_e3174` | Molecular weight e3174 | $436.894$ | g/mol | ✓ |
| $\mathrm{BP}_{\mathrm{e3174}}$ | `BP_e3174` | blood to plasma ratio e3174 | $0.58$ | \- | ✓ |
| $\mathrm{ftissue}_{\mathrm{e3174}}$ | `ftissue_e3174` | tissue distribution E3174 | $0.00217361854744583$ | l/min | ✓ |
| $\mathrm{Kp}_{\mathrm{e3174}}$ | `Kp_e3174` | tissue/blood partition coefficient E3174 | $0.185665446932989$ | \- | ✓ |
| $\mathrm{IVDOSE}_{\mathrm{e3174}}$ | `IVDOSE_e3174` | IV bolus dose e3174 \[mg\] | $0$ | mg |  |
| $\mathrm{ti}_{\mathrm{e3174}}$ | `ti_e3174` | injection time e3174 | $10$ | s | ✓ |
| $\mathrm{Ki}_{\mathrm{e3174}}$ | `Ki_e3174` | Ki \[1/min\] injection e3174 |  | 1/min |  |
| $\mathrm{Ri}_{\mathrm{e3174}}$ | `Ri_e3174` | Ri \[mg/min\] rate of infusion e3174 | $0$ | mg/min | ✓ |
| $\mathrm{cum}_{\mathrm{dose\_e3174}}$ | `cum_dose_e3174` | Cumulative dose due to infusion e3174 | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{l158}}$ | `Mr_l158` | Molecular weight l158 | $438.9$ | g/mol | ✓ |
| $\mathrm{BP}_{\mathrm{l158}}$ | `BP_l158` | blood to plasma ratio l158 | $1$ | \- | ✓ |
| $\mathrm{mr}_{\mathrm{e3174\_los\_plasma}}$ | `mr_e3174_los_plasma` |  |  | \- |  |
| $\mathrm{mr}_{\mathrm{e3174\_los\_urine}}$ | `mr_e3174_los_urine` |  |  | \- |  |
| $\mathrm{mr}_{\mathrm{e3174\_los\_feces}}$ | `mr_e3174_los_feces` |  |  | \- |  |
| $\mathrm{mr}_{\mathrm{los\_e3174\_plasma}}$ | `mr_los_e3174_plasma` |  |  | \- |  |
| $\mathrm{mr}_{\mathrm{los\_e3174\_urine}}$ | `mr_los_e3174_urine` |  |  | \- |  |
| $\mathrm{mr}_{\mathrm{los\_e3174\_feces}}$ | `mr_los_e3174_feces` |  |  | \- |  |
| $\mathrm{Aurine}_{\mathrm{los\_e3174}}$ | `Aurine_los_e3174` |  |  | mmol |  |
| $\mathrm{Aurine}_{\mathrm{total}}$ | `Aurine_total` |  |  | mmol |  |
| $\mathrm{Afeces}_{\mathrm{total}}$ | `Afeces_total` |  |  | mmol |  |
| $\mathrm{RENSEC}_{\mathrm{k}}$ | `RENSEC_k` | rate renin secretion |  | mmol/min |  |
| $\mathrm{RENDEG}_{\mathrm{k}}$ | `RENDEG_k` | rate renin degradation |  | l/min |  |
| $\mathrm{ALDSEC}_{\mathrm{k}}$ | `ALDSEC_k` | rate aldosterone secretion |  | mmol/min |  |
| $\mathrm{ALDDEG}_{\mathrm{k}}$ | `ALDDEG_k` | rate aldosterone degradation |  | l/min |  |
| $\mathrm{ANGGEN2ANG1}_{\mathrm{k}}$ | `ANGGEN2ANG1_k` | rate angen → ang1 |  | l/min |  |
| $\mathrm{ANG1ANG2}_{\mathrm{k}}$ | `ANG1ANG2_k` | rate ang1 → ang2 |  | l/min |  |
| $\mathrm{ANG2DEG}_{\mathrm{k}}$ | `ANG2DEG_k` | rate angiotensin II degradation |  | l/min |  |
| $\mathrm{fe}_{\mathrm{e3174}}$ | `fe_e3174` | AT1 receptor blockade by E3174 in \[0, 1\] |  | \- |  |
| $f_{\mathrm{at1}}$ | `f_at1` | AT1 signal relative to the baseline |  | \- |  |
| $\mathrm{Qgu}$ | `Qgu` | gut blood flow |  | l/min |  |
| $\mathrm{Qki}$ | `Qki` | kidney blood flow |  | l/min |  |
| $\mathrm{Qh}$ | `Qh` | hepatic (venous side) blood flow |  | l/min |  |
| $\mathrm{Qha}$ | `Qha` | hepatic artery blood flow |  | l/min |  |
| $\mathrm{Qlu}$ | `Qlu` | lung blood flow |  | l/min |  |
| $\mathrm{Qre}$ | `Qre` | rest of body blood flow |  | l/min |  |
| $\mathrm{Qpo}$ | `Qpo` | portal blood flow |  | l/min |  |
| $\mathrm{KI}_{\mathrm{\_f\_renal\_function}}$ | `KI__f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_fu\_los}}$ | `KI__fu_los` | fraction unbound in plasma losartan | $0.013$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_fu\_e3174}}$ | `KI__fu_e3174` | fraction unbound in plasma E3174 | $0.002$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_LOSEX\_f}}$ | `KI__LOSEX_f` | renal clearance of unbound losartan relative to the GFR | $59.6320404780412$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_E3174EX\_f}}$ | `KI__E3174EX_f` | renal clearance of unbound E3174 relative to the GFR | $124.916497816236$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_L158EX\_f}}$ | `KI__L158EX_f` | renal clearance of the other metabolites relative to the GFR | $5.18711908659398$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_GFR}}$ | `KI__GFR` | glomerular filtration rate |  | ml/min |  |
| $\mathrm{LI}_{\mathrm{\_f\_cyp2c9}}$ | `LI__f_cyp2c9` | activity of CYP2C9 | $1$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_f\_l158\_cyp2c9}}$ | `LI__f_l158_cyp2c9` | fraction of the formation of the other metabolites by CYP2C9 | $0.5$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_LOSIM\_k}}$ | `LI__LOSIM_k` | rate losartan import | $50$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_LOS2E3174\_k}}$ | `LI__LOS2E3174_k` | rate losartan → E3174 | $0.412801474839129$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_LOS2L158\_k}}$ | `LI__LOS2L158_k` | rate losartan → other metabolites | $1.78636406676044$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_E3174EX\_k}}$ | `LI__E3174EX_k` | rate E3174 export | $0.00411279765463081$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_E3174IM\_k}}$ | `LI__E3174IM_k` | rate E3174 import | $0.0614100473485477$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_L158EX\_k}}$ | `LI__L158EX_k` | rate export of the other metabolites | $1$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_E3174BIEX\_k}}$ | `LI__E3174BIEX_k` | rate E3174 bile export | $0.00274339281823762$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_L158BIEX\_k}}$ | `LI__L158BIEX_k` | rate bile export of the other metabolites | $0.257303343770967$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_LOSABS\_k}}$ | `GU__LOSABS_k` | rate of losartan absorption | $0.0476385224387545$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_METEXC\_k}}$ | `GU__METEXC_k` | rate of feces excretion | $0.0006$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_Ka\_dis\_los}}$ | `GU__Ka_dis_los` | rate dissolution losartan | $2$ | 1/hr | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_OATP2B1}}$ | `GU__f_OATP2B1` | absorption activity | $1$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_abcb1}}$ | `GU__f_abcb1` | PG activity | $1$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_LOSEFL\_k}}$ | `GU__f_LOSEFL_k` | fractional rate losartan efflux (PG) | $1$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_LOSEFL\_k}}$ | `GU__LOSEFL_k` |  |  | 1/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{ald} &= \mathrm{ald}_{\mathrm{ref}} \\[1ex]
\mathrm{ren} &= \mathrm{ren}_{\mathrm{ref}} \\[1ex]
\mathrm{anggen} &= \mathrm{anggen}_{\mathrm{ref}} \\[1ex]
\mathrm{ang}_{1} &= \mathrm{ang1}_{\mathrm{ref}} \\[1ex]
\mathrm{ang}_{2} &= \mathrm{ang2}_{\mathrm{ref}} \\[1ex]
\mathrm{Aurine}_{\mathrm{los}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{los}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{Aurine}_{\mathrm{e3174}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{e3174}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{Aurine}_{\mathrm{l158}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{l158}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{LI}_{\mathrm{\_e3174\_bi}} &= 0 \cdot \mathrm{LI}_{\mathrm{\_Vbi}} \\[1ex]
\mathrm{LI}_{\mathrm{\_l158\_bi}} &= 0 \cdot \mathrm{LI}_{\mathrm{\_Vbi}} \\[1ex]
\mathrm{GU}_{\mathrm{\_los\_stomach}} &= 0 \cdot \mathrm{Vstomach}
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
f_{\mathrm{shunts}} &= f_{\mathrm{cirrhosis}} \\[1ex]
f_{\mathrm{tissue\_loss}} &= f_{\mathrm{cirrhosis}} \\[1ex]
\mathrm{FVre} &= 1 - \mathrm{FVgu} - \mathrm{FVki} - \mathrm{FVli} \\
&\quad {} - \mathrm{FVlu} - \mathrm{FVve} - \mathrm{FVar} - \mathrm{FVpo} \\
&\quad {} - \mathrm{FVhv} \\[1ex]
\mathrm{FQre} &= 1 - \mathrm{FQki} - \mathrm{FQh} \\[1ex]
\mathrm{CO} &= f_{\mathrm{cardiac\_function}} \cdot \mathrm{BW} \cdot \mathrm{COBW} \\[1ex]
\mathrm{QC} &= \frac{\mathrm{CO}}{1000} \cdot 60 \\[1ex]
\mathrm{Vgu} &= \mathrm{BW} \cdot \mathrm{FVgu} \\[1ex]
\mathrm{Vki} &= \mathrm{BW} \cdot \mathrm{FVki} \\[1ex]
\mathrm{Vli} &= \mathrm{BW} \cdot \mathrm{FVli} \\[1ex]
\mathrm{Vlu} &= \mathrm{BW} \cdot \mathrm{FVlu} \\[1ex]
\mathrm{Vre} &= \mathrm{BW} \cdot \mathrm{FVre} \\[1ex]
\mathrm{Vve} &= \mathrm{BW} \cdot \mathrm{FVve} - \frac{\mathrm{FVve}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Var} &= \mathrm{BW} \cdot \mathrm{FVar} - \frac{\mathrm{FVar}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Vpo} &= \mathrm{BW} \cdot \mathrm{FVpo} - \frac{\mathrm{FVpo}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right) \\[1ex]
\mathrm{Vhv} &= \mathrm{BW} \cdot \mathrm{FVhv} - \frac{\mathrm{FVhv}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right) \\[1ex]
\mathrm{Qgu} &= \mathrm{QC} \cdot \mathrm{FQgu} \\[1ex]
\mathrm{Qki} &= \mathrm{QC} \cdot \mathrm{FQki} \\[1ex]
\mathrm{Qh} &= \mathrm{QC} \cdot \mathrm{FQh} \\[1ex]
\mathrm{Qha} &= \mathrm{Qh} - \mathrm{Qgu} \\[1ex]
\mathrm{Qlu} &= \mathrm{QC} \cdot \mathrm{FQlu} \\[1ex]
\mathrm{Qre} &= \mathrm{QC} \cdot \mathrm{FQre} \\[1ex]
\mathrm{Qpo} &= \mathrm{Qgu} \\[1ex]
\mathrm{GFR}_{\mathrm{base}} &= f_{\mathrm{GFR\_Qki}} \cdot \mathrm{Qki} \\[1ex]
\mathrm{Vki}_{\mathrm{blood}} &= \mathrm{Vki} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vki}_{\mathrm{tissue}} &= \mathrm{Vki} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{blood}} &= \mathrm{Vli} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vli}_{\mathrm{tissue}} &= \mathrm{Vli} \cdot \mathopen{}\left(1 - f_{\mathrm{tissue\_loss}}\right) \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{blood}} &= \mathrm{Vlu} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vlu}_{\mathrm{tissue}} &= \mathrm{Vlu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{blood}} &= \mathrm{Vgu} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vgu}_{\mathrm{tissue}} &= \mathrm{Vgu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{blood}} &= \mathrm{Vre} \cdot \mathrm{Fblood} \\[1ex]
\mathrm{Vre}_{\mathrm{tissue}} &= \mathrm{Vre} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Ki}_{\mathrm{los}} &= \frac{0.693}{\mathrm{ti}_{\mathrm{los}}} \cdot 60 \\[1ex]
\mathrm{Ki}_{\mathrm{e3174}} &= \frac{0.693}{\mathrm{ti}_{\mathrm{e3174}}} \cdot 60 \\[1ex]
\mathrm{mr}_{\mathrm{e3174\_los\_urine}} &= \frac{\mathrm{Aurine}_{\mathrm{e3174}}}{\mathrm{Aurine}_{\mathrm{los}} + 10^{-12}} \\[1ex]
\mathrm{mr}_{\mathrm{e3174\_los\_feces}} &= \frac{\mathrm{Afeces}_{\mathrm{e3174}}}{\mathrm{Afeces}_{\mathrm{los}} + 10^{-12}} \\[1ex]
\mathrm{mr}_{\mathrm{los\_e3174\_urine}} &= \frac{\mathrm{Aurine}_{\mathrm{los}}}{\mathrm{Aurine}_{\mathrm{e3174}} + 10^{-12}} \\[1ex]
\mathrm{mr}_{\mathrm{los\_e3174\_feces}} &= \frac{\mathrm{Afeces}_{\mathrm{los}}}{\mathrm{Afeces}_{\mathrm{e3174}} + 10^{-12}} \\[1ex]
\mathrm{Aurine}_{\mathrm{los\_e3174}} &= \mathrm{Aurine}_{\mathrm{los}} + \mathrm{Aurine}_{\mathrm{e3174}} \\[1ex]
\mathrm{Cveplasma}_{\mathrm{los}} &= \frac{\mathrm{Cve}_{\mathrm{los}}}{\mathrm{BP}_{\mathrm{los}}} \\[1ex]
\mathrm{Cveplasma}_{\mathrm{e3174}} &= \frac{\mathrm{Cve}_{\mathrm{e3174}}}{\mathrm{BP}_{\mathrm{e3174}}} \\[1ex]
\mathrm{mr}_{\mathrm{e3174\_los\_plasma}} &= \frac{\mathrm{Cveplasma}_{\mathrm{e3174}}}{\mathrm{Cveplasma}_{\mathrm{los}} + 10^{-12}} \\[1ex]
\mathrm{mr}_{\mathrm{los\_e3174\_plasma}} &= \frac{\mathrm{Cveplasma}_{\mathrm{los}}}{\mathrm{Cveplasma}_{\mathrm{e3174}} + 10^{-12}} \\[1ex]
\mathrm{Cveplasma}_{\mathrm{l158}} &= \frac{\mathrm{Cve}_{\mathrm{l158}}}{\mathrm{BP}_{\mathrm{l158}}} \\[1ex]
\mathrm{Cveplasma}_{\mathrm{los\_e3174}} &= \frac{\mathrm{Cve}_{\mathrm{los}}}{\mathrm{BP}_{\mathrm{los}}} + \frac{\mathrm{Cve}_{\mathrm{e3174}}}{\mathrm{BP}_{\mathrm{e3174}}} \\[1ex]
\mathrm{Cveplasma}_{\mathrm{total}} &= \frac{\mathrm{Cve}_{\mathrm{los}}}{\mathrm{BP}_{\mathrm{los}}} + \frac{\mathrm{Cve}_{\mathrm{e3174}}}{\mathrm{BP}_{\mathrm{e3174}}} + \frac{\mathrm{Cve}_{\mathrm{l158}}}{\mathrm{BP}_{\mathrm{l158}}} \\[1ex]
\mathrm{Aurine}_{\mathrm{total}} &= \mathrm{Aurine}_{\mathrm{los}} + \mathrm{Aurine}_{\mathrm{e3174}} + \mathrm{Aurine}_{\mathrm{l158}} \\[1ex]
\mathrm{Afeces}_{\mathrm{total}} &= \mathrm{Afeces}_{\mathrm{los}} + \mathrm{Afeces}_{\mathrm{e3174}} + \mathrm{Afeces}_{\mathrm{l158}} \\[1ex]
e_{3174} &= \mathrm{Cveplasma}_{\mathrm{e3174}} \\[1ex]
\mathrm{fe}_{\mathrm{e3174}} &= \frac{e_{3174}}{\mathrm{E50}_{\mathrm{e3174}} + e_{3174}} \\[1ex]
f_{\mathrm{at1}} &= \frac{\mathrm{ang}_{2}}{\mathrm{ang2}_{\mathrm{ref}}} \cdot \mathopen{}\left(1 - \mathrm{fe}_{\mathrm{e3174}}\right) \\[1ex]
\mathrm{SBP} &= \mathrm{SBP}_{\mathrm{ref}} \cdot \mathopen{}\left(1 - \mathrm{BP}_{\mathrm{emax}} \cdot \mathopen{}\left(1 - f_{\mathrm{at1}}\right)\right) \\[1ex]
\mathrm{DBP} &= \mathrm{DBP}_{\mathrm{ref}} \cdot \mathopen{}\left(1 - \mathrm{BP}_{\mathrm{emax}} \cdot \mathopen{}\left(1 - f_{\mathrm{at1}}\right)\right) \\[1ex]
\mathrm{MAP} &= \mathrm{DBP} + \frac{\mathrm{SBP} - \mathrm{DBP}}{3} \\[1ex]
\mathrm{SBP}_{\mathrm{change}} &= \mathrm{SBP} - \mathrm{SBP}_{\mathrm{ref}} \\[1ex]
\mathrm{SBP}_{\mathrm{ratio}} &= \frac{\mathrm{SBP}}{\mathrm{SBP}_{\mathrm{ref}}} \\[1ex]
\mathrm{DBP}_{\mathrm{change}} &= \mathrm{DBP} - \mathrm{DBP}_{\mathrm{ref}} \\[1ex]
\mathrm{DBP}_{\mathrm{ratio}} &= \frac{\mathrm{DBP}}{\mathrm{DBP}_{\mathrm{ref}}} \\[1ex]
\mathrm{KI}_{\mathrm{\_GFR}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot \mathrm{GFR}_{\mathrm{base}} \\[1ex]
\mathrm{GU}_{\mathrm{\_LOSEFL\_k}} &= \mathrm{GU}_{\mathrm{\_f\_LOSEFL\_k}} \cdot \mathrm{GU}_{\mathrm{\_LOSABS\_k}}
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
| $v_{\mathrm{transport\_lu\_los}}$ | `transport_lu_los` | transport losartan | $\mathrm{Clu}_{\mathrm{blood\_los}} \rightleftharpoons \mathrm{Clu}_{\mathrm{los}}$ |  |
| $v_{\mathrm{transport\_re\_los}}$ | `transport_re_los` | transport losartan | $\mathrm{Cre}_{\mathrm{blood\_los}} \rightleftharpoons \mathrm{Cre}_{\mathrm{los}}$ |  |
| $v_{\mathrm{iv\_los}}$ | `iv_los` | iv losartan | $\varnothing \longrightarrow \mathrm{Cve}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_los}}$ | `Flow_ar_ki_los` | inflow kidney losartan | $\mathrm{Car}_{\mathrm{los}} \longrightarrow \mathrm{Cki}_{\mathrm{blood\_los}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_los}}$ | `Flow_ki_ve_los` | outflow kidney losartan | $\mathrm{Cki}_{\mathrm{blood\_los}} \longrightarrow \mathrm{Cve}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_los}}$ | `Flow_arli_li_los` | arterial inflow liver losartan | $\mathrm{Car}_{\mathrm{los}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_los}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_los}}$ | `Flow_arli_hv_los` | flow arterial shunts | $\mathrm{Car}_{\mathrm{los}} \longrightarrow \mathrm{Chv}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_los}}$ | `Flow_po_li_los` | outflow po losartan | $\mathrm{Cpo}_{\mathrm{los}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_los}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_los}}$ | `Flow_po_hv_los` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{los}} \longrightarrow \mathrm{Chv}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_los}}$ | `Flow_li_hv_los` | outflow liver losartan | $\mathrm{Cli}_{\mathrm{blood\_los}} \longrightarrow \mathrm{Chv}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_los}}$ | `Flow_hv_ve_los` | outflow hepatic vein losartan | $\mathrm{Chv}_{\mathrm{los}} \longrightarrow \mathrm{Cve}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_los}}$ | `Flow_ve_lu_los` | inflow lung losartan | $\mathrm{Cve}_{\mathrm{los}} \longrightarrow \mathrm{Clu}_{\mathrm{blood\_los}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_los}}$ | `Flow_lu_ar_los` | outflow lung losartan | $\mathrm{Clu}_{\mathrm{blood\_los}} \longrightarrow \mathrm{Car}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_los}}$ | `Flow_ar_gu_los` | inflow gut losartan | $\mathrm{Car}_{\mathrm{los}} \longrightarrow \mathrm{Cgu}_{\mathrm{blood\_los}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_los}}$ | `Flow_gu_po_los` | outflow gut losartan | $\mathrm{Cgu}_{\mathrm{blood\_los}} \longrightarrow \mathrm{Cpo}_{\mathrm{los}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_los}}$ | `Flow_ar_re_los` | inflow rest losartan | $\mathrm{Car}_{\mathrm{los}} \longrightarrow \mathrm{Cre}_{\mathrm{blood\_los}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_los}}$ | `Flow_re_ve_los` | outflow rest losartan | $\mathrm{Cre}_{\mathrm{blood\_los}} \longrightarrow \mathrm{Cve}_{\mathrm{los}}$ |  |
| $v_{\mathrm{transport\_lu\_e3174}}$ | `transport_lu_e3174` | transport E3174 | $\mathrm{Clu}_{\mathrm{blood\_e3174}} \rightleftharpoons \mathrm{Clu}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{transport\_re\_e3174}}$ | `transport_re_e3174` | transport E3174 | $\mathrm{Cre}_{\mathrm{blood\_e3174}} \rightleftharpoons \mathrm{Cre}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{iv\_e3174}}$ | `iv_e3174` | iv E3174 | $\varnothing \longrightarrow \mathrm{Cve}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_e3174}}$ | `Flow_ar_ki_e3174` | inflow kidney E3174 | $\mathrm{Car}_{\mathrm{e3174}} \longrightarrow \mathrm{Cki}_{\mathrm{blood\_e3174}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_e3174}}$ | `Flow_ki_ve_e3174` | outflow kidney E3174 | $\mathrm{Cki}_{\mathrm{blood\_e3174}} \longrightarrow \mathrm{Cve}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_e3174}}$ | `Flow_arli_li_e3174` | arterial inflow liver E3174 | $\mathrm{Car}_{\mathrm{e3174}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_e3174}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_e3174}}$ | `Flow_arli_hv_e3174` | flow arterial shunts | $\mathrm{Car}_{\mathrm{e3174}} \longrightarrow \mathrm{Chv}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_e3174}}$ | `Flow_po_li_e3174` | outflow po E3174 | $\mathrm{Cpo}_{\mathrm{e3174}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_e3174}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_e3174}}$ | `Flow_po_hv_e3174` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{e3174}} \longrightarrow \mathrm{Chv}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_e3174}}$ | `Flow_li_hv_e3174` | outflow liver E3174 | $\mathrm{Cli}_{\mathrm{blood\_e3174}} \longrightarrow \mathrm{Chv}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_e3174}}$ | `Flow_hv_ve_e3174` | outflow hepatic vein E3174 | $\mathrm{Chv}_{\mathrm{e3174}} \longrightarrow \mathrm{Cve}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_e3174}}$ | `Flow_ve_lu_e3174` | inflow lung E3174 | $\mathrm{Cve}_{\mathrm{e3174}} \longrightarrow \mathrm{Clu}_{\mathrm{blood\_e3174}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_e3174}}$ | `Flow_lu_ar_e3174` | outflow lung E3174 | $\mathrm{Clu}_{\mathrm{blood\_e3174}} \longrightarrow \mathrm{Car}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_e3174}}$ | `Flow_ar_gu_e3174` | inflow gut E3174 | $\mathrm{Car}_{\mathrm{e3174}} \longrightarrow \mathrm{Cgu}_{\mathrm{blood\_e3174}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_e3174}}$ | `Flow_gu_po_e3174` | outflow gut E3174 | $\mathrm{Cgu}_{\mathrm{blood\_e3174}} \longrightarrow \mathrm{Cpo}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_e3174}}$ | `Flow_ar_re_e3174` | inflow rest E3174 | $\mathrm{Car}_{\mathrm{e3174}} \longrightarrow \mathrm{Cre}_{\mathrm{blood\_e3174}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_e3174}}$ | `Flow_re_ve_e3174` | outflow rest E3174 | $\mathrm{Cre}_{\mathrm{blood\_e3174}} \longrightarrow \mathrm{Cve}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_l158}}$ | `Flow_ar_ki_l158` | inflow kidney other losartan metabolites | $\mathrm{Car}_{\mathrm{l158}} \longrightarrow \mathrm{Cki}_{\mathrm{blood\_l158}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_l158}}$ | `Flow_ki_ve_l158` | outflow kidney other losartan metabolites | $\mathrm{Cki}_{\mathrm{blood\_l158}} \longrightarrow \mathrm{Cve}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_l158}}$ | `Flow_arli_li_l158` | arterial inflow liver other losartan metabolites | $\mathrm{Car}_{\mathrm{l158}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_l158}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_l158}}$ | `Flow_arli_hv_l158` | flow arterial shunts | $\mathrm{Car}_{\mathrm{l158}} \longrightarrow \mathrm{Chv}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_l158}}$ | `Flow_po_li_l158` | outflow po other losartan metabolites | $\mathrm{Cpo}_{\mathrm{l158}} \longrightarrow \mathrm{Cli}_{\mathrm{blood\_l158}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_l158}}$ | `Flow_po_hv_l158` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{l158}} \longrightarrow \mathrm{Chv}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_l158}}$ | `Flow_li_hv_l158` | outflow liver other losartan metabolites | $\mathrm{Cli}_{\mathrm{blood\_l158}} \longrightarrow \mathrm{Chv}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_l158}}$ | `Flow_hv_ve_l158` | outflow hepatic vein other losartan metabolites | $\mathrm{Chv}_{\mathrm{l158}} \longrightarrow \mathrm{Cve}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_l158}}$ | `Flow_ve_lu_l158` | inflow lung other losartan metabolites | $\mathrm{Cve}_{\mathrm{l158}} \longrightarrow \mathrm{Clu}_{\mathrm{blood\_l158}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_l158}}$ | `Flow_lu_ar_l158` | outflow lung other losartan metabolites | $\mathrm{Clu}_{\mathrm{blood\_l158}} \longrightarrow \mathrm{Car}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_l158}}$ | `Flow_ar_gu_l158` | inflow gut other losartan metabolites | $\mathrm{Car}_{\mathrm{l158}} \longrightarrow \mathrm{Cgu}_{\mathrm{blood\_l158}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_l158}}$ | `Flow_gu_po_l158` | outflow gut other losartan metabolites | $\mathrm{Cgu}_{\mathrm{blood\_l158}} \longrightarrow \mathrm{Cpo}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_l158}}$ | `Flow_ar_re_l158` | inflow rest other losartan metabolites | $\mathrm{Car}_{\mathrm{l158}} \longrightarrow \mathrm{Cre}_{\mathrm{blood\_l158}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_l158}}$ | `Flow_re_ve_l158` | outflow rest other losartan metabolites | $\mathrm{Cre}_{\mathrm{blood\_l158}} \longrightarrow \mathrm{Cve}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{KI\_\_LOSEX}}$ | `KI__LOSEX` | losartan excretion (LOSEX) | $\mathrm{Cki}_{\mathrm{blood\_los}} \longrightarrow \mathrm{Aurine}_{\mathrm{los}}$ |  |
| $v_{\mathrm{KI\_\_E3174EX}}$ | `KI__E3174EX` | E3174 excretion (E3174EX) | $\mathrm{Cki}_{\mathrm{blood\_e3174}} \longrightarrow \mathrm{Aurine}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{KI\_\_L158EX}}$ | `KI__L158EX` | other metabolites excretion (L158EX) | $\mathrm{Cki}_{\mathrm{blood\_l158}} \longrightarrow \mathrm{Aurine}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{LI\_\_LOSIM}}$ | `LI__LOSIM` | losartan import (LOSIM) | $\mathrm{Cli}_{\mathrm{blood\_los}} \rightleftharpoons \mathrm{LI}_{\mathrm{\_los}}$ |  |
| $v_{\mathrm{LI\_\_LOS2E3174}}$ | `LI__LOS2E3174` | E3174 formation (CYP2C9) | $\mathrm{LI}_{\mathrm{\_los}} \longrightarrow \mathrm{LI}_{\mathrm{\_e3174}}$ |  |
| $v_{\mathrm{LI\_\_LOS2L158}}$ | `LI__LOS2L158` | other metabolites formation (CYP2C9, CYP3A4, UGT) | $\mathrm{LI}_{\mathrm{\_los}} \longrightarrow \mathrm{LI}_{\mathrm{\_l158}}$ |  |
| $v_{\mathrm{LI\_\_E3174EX}}$ | `LI__E3174EX` | E3174 export (E3174EX) | $\mathrm{LI}_{\mathrm{\_e3174}} \rightleftharpoons \mathrm{Cli}_{\mathrm{blood\_e3174}}$ |  |
| $v_{\mathrm{LI\_\_L158EX}}$ | `LI__L158EX` | other metabolites export (L158EX) | $\mathrm{LI}_{\mathrm{\_l158}} \rightleftharpoons \mathrm{Cli}_{\mathrm{blood\_l158}}$ |  |
| $v_{\mathrm{LI\_\_E3174BIEX}}$ | `LI__E3174BIEX` | E3174 bile export | $\mathrm{LI}_{\mathrm{\_e3174}} \longrightarrow \mathrm{LI}_{\mathrm{\_e3174\_bi}}$ |  |
| $v_{\mathrm{LI\_\_L158BIEX}}$ | `LI__L158BIEX` | other metabolites bile export | $\mathrm{LI}_{\mathrm{\_l158}} \longrightarrow \mathrm{LI}_{\mathrm{\_l158\_bi}}$ |  |
| $v_{\mathrm{LI\_\_E3174EHC}}$ | `LI__E3174EHC` | E3174 bile into the gut lumen | $\mathrm{LI}_{\mathrm{\_e3174\_bi}} \longrightarrow \mathrm{Cgu}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{LI\_\_L158EHC}}$ | `LI__L158EHC` | other losartan metabolites bile into the gut lumen | $\mathrm{LI}_{\mathrm{\_l158\_bi}} \longrightarrow \mathrm{Cgu}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{GU\_\_LOSABS}}$ | `GU__LOSABS` | absorption losartan | $\mathrm{Cgu}_{\mathrm{los}} \longrightarrow \mathrm{GU}_{\mathrm{\_los}}$ |  |
| $v_{\mathrm{GU\_\_LOSEFL}}$ | `GU__LOSEFL` | efflux losartan (PG) | $\mathrm{GU}_{\mathrm{\_los}} \longrightarrow \mathrm{Cgu}_{\mathrm{los}}$ |  |
| $v_{\mathrm{GU\_\_LOSEX}}$ | `GU__LOSEX` | losartan export plasma | $\mathrm{GU}_{\mathrm{\_los}} \longrightarrow \mathrm{Cgu}_{\mathrm{blood\_los}}$ |  |
| $v_{\mathrm{GU\_\_LOSEXC}}$ | `GU__LOSEXC` | excretion losartan (feces) | $\mathrm{Cgu}_{\mathrm{los}} \longrightarrow \mathrm{Afeces}_{\mathrm{los}}$ |  |
| $v_{\mathrm{GU\_\_E3174EXC}}$ | `GU__E3174EXC` | excretion E3174 (feces) | $\mathrm{Cgu}_{\mathrm{e3174}} \longrightarrow \mathrm{Afeces}_{\mathrm{e3174}}$ |  |
| $v_{\mathrm{GU\_\_L158EXC}}$ | `GU__L158EXC` | excretion other metabolites (feces) | $\mathrm{Cgu}_{\mathrm{l158}} \longrightarrow \mathrm{Afeces}_{\mathrm{l158}}$ |  |
| $v_{\mathrm{GU\_\_dissolution\_los}}$ | `GU__dissolution_los` | dissolution losartan | $\mathrm{GU}_{\mathrm{\_los\_stomach}} \longrightarrow \mathrm{Cgu}_{\mathrm{los}}$ |  |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{RENSEC}} &= \mathrm{RENSEC}_{\mathrm{k}} \cdot f_{\mathrm{at1}}^{-\mathrm{RENSEC}_{\mathrm{n}}} \\[1ex]
v_{\mathrm{RENDEG}} &= \mathrm{RENDEG}_{\mathrm{k}} \cdot \mathrm{ren} \\[1ex]
v_{\mathrm{ALDSEC}} &= \mathrm{ALDSEC}_{\mathrm{k}} \cdot f_{\mathrm{at1}} \\[1ex]
v_{\mathrm{ALDDEG}} &= \mathrm{ALDDEG}_{\mathrm{k}} \cdot \mathrm{ald} \\[1ex]
v_{\mathrm{ANGGEN2ANG1}} &= \frac{\mathrm{ANGGEN2ANG1}_{\mathrm{k}} \cdot \mathrm{anggen} \cdot \mathrm{ren}}{\mathrm{ren}_{\mathrm{ref}}} \\[1ex]
v_{\mathrm{ANG1ANG2}} &= \mathrm{ANG1ANG2}_{\mathrm{k}} \cdot \mathrm{ang}_{1} \\[1ex]
v_{\mathrm{ANG2DEG}} &= \mathrm{ANG2DEG}_{\mathrm{k}} \cdot \mathrm{ang}_{2} \\[1ex]
v_{\mathrm{transport\_lu\_los}} &= \frac{\mathrm{ftissue}_{\mathrm{los}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Clu}_{\mathrm{blood\_los}} \cdot \mathrm{Kp}_{\mathrm{los}} - \mathrm{Clu}_{\mathrm{los}}\right) \\[1ex]
v_{\mathrm{transport\_re\_los}} &= \frac{\mathrm{ftissue}_{\mathrm{los}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Cre}_{\mathrm{blood\_los}} \cdot \mathrm{Kp}_{\mathrm{los}} - \mathrm{Cre}_{\mathrm{los}}\right) \\[1ex]
v_{\mathrm{iv\_los}} &= \frac{\mathrm{Ki}_{\mathrm{los}} \cdot \mathrm{IVDOSE}_{\mathrm{los}}}{\mathrm{Mr}_{\mathrm{los}}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_los}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_los}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{blood\_los}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_los}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_los}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_los}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_los}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_los}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{blood\_los}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_los}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_los}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_los}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{blood\_los}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_los}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_los}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{blood\_los}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_los}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{los}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_los}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{blood\_los}} \\[1ex]
v_{\mathrm{transport\_lu\_e3174}} &= \frac{\mathrm{ftissue}_{\mathrm{e3174}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Clu}_{\mathrm{blood\_e3174}} \cdot \mathrm{Kp}_{\mathrm{e3174}} - \mathrm{Clu}_{\mathrm{e3174}}\right) \\[1ex]
v_{\mathrm{transport\_re\_e3174}} &= \frac{\mathrm{ftissue}_{\mathrm{e3174}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Cre}_{\mathrm{blood\_e3174}} \cdot \mathrm{Kp}_{\mathrm{e3174}} - \mathrm{Cre}_{\mathrm{e3174}}\right) \\[1ex]
v_{\mathrm{iv\_e3174}} &= \frac{\mathrm{Ki}_{\mathrm{e3174}} \cdot \mathrm{IVDOSE}_{\mathrm{e3174}}}{\mathrm{Mr}_{\mathrm{e3174}}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_e3174}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_e3174}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{blood\_e3174}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_e3174}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_e3174}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_e3174}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_e3174}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_e3174}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{blood\_e3174}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_e3174}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_e3174}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_e3174}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{blood\_e3174}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_e3174}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_e3174}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{blood\_e3174}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_e3174}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_e3174}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{blood\_e3174}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_l158}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_l158}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{blood\_l158}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_l158}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_l158}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_l158}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_l158}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_l158}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{blood\_l158}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_l158}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_l158}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_l158}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{blood\_l158}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_l158}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_l158}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{blood\_l158}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_l158}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_l158}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{blood\_l158}} \\[1ex]
v_{\mathrm{KI\_\_LOSEX}} &= \frac{\mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_fu\_los}} \cdot \mathrm{KI}_{\mathrm{\_LOSEX\_f}} \cdot \mathrm{Cki}_{\mathrm{blood\_los}}}{\mathrm{BP}_{\mathrm{los}}} \\[1ex]
v_{\mathrm{KI\_\_E3174EX}} &= \frac{\mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_fu\_e3174}} \cdot \mathrm{KI}_{\mathrm{\_E3174EX\_f}} \cdot \mathrm{Cki}_{\mathrm{blood\_e3174}}}{\mathrm{BP}_{\mathrm{e3174}}} \\[1ex]
v_{\mathrm{KI\_\_L158EX}} &= \frac{\mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_L158EX\_f}} \cdot \mathrm{Cki}_{\mathrm{blood\_l158}}}{\mathrm{BP}_{\mathrm{l158}}} \\[1ex]
v_{\mathrm{LI\_\_LOSIM}} &= \mathrm{LI}_{\mathrm{\_LOSIM\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{Cli}_{\mathrm{blood\_los}} - \mathrm{LI}_{\mathrm{\_los}}\right) \\[1ex]
v_{\mathrm{LI\_\_LOS2E3174}} &= \mathrm{LI}_{\mathrm{\_f\_cyp2c9}} \cdot \mathrm{LI}_{\mathrm{\_LOS2E3174\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_los}} \\[1ex]
v_{\mathrm{LI\_\_LOS2L158}} &= \mathopen{}\left(\mathrm{LI}_{\mathrm{\_f\_l158\_cyp2c9}} \cdot \mathrm{LI}_{\mathrm{\_f\_cyp2c9}} + 1 - \mathrm{LI}_{\mathrm{\_f\_l158\_cyp2c9}}\right) \cdot \mathrm{LI}_{\mathrm{\_LOS2L158\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_los}} \\[1ex]
v_{\mathrm{LI\_\_E3174EX}} &= \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{LI}_{\mathrm{\_E3174EX\_k}} \cdot \mathrm{LI}_{\mathrm{\_e3174}} - \mathrm{LI}_{\mathrm{\_E3174IM\_k}} \cdot \mathrm{Cli}_{\mathrm{blood\_e3174}}\right) \\[1ex]
v_{\mathrm{LI\_\_L158EX}} &= \mathrm{LI}_{\mathrm{\_L158EX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathopen{}\left(\mathrm{LI}_{\mathrm{\_l158}} - \mathrm{Cli}_{\mathrm{blood\_l158}}\right) \\[1ex]
v_{\mathrm{LI\_\_E3174BIEX}} &= \mathrm{LI}_{\mathrm{\_E3174BIEX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_e3174}} \\[1ex]
v_{\mathrm{LI\_\_L158BIEX}} &= \mathrm{LI}_{\mathrm{\_L158BIEX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_l158}} \\[1ex]
v_{\mathrm{LI\_\_E3174EHC}} &= v_{\mathrm{LI\_\_E3174BIEX}} \\[1ex]
v_{\mathrm{LI\_\_L158EHC}} &= v_{\mathrm{LI\_\_L158BIEX}} \\[1ex]
v_{\mathrm{GU\_\_LOSABS}} &= \mathrm{GU}_{\mathrm{\_f\_OATP2B1}} \cdot \mathrm{GU}_{\mathrm{\_LOSABS\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{Cgu}_{\mathrm{los}} \\[1ex]
v_{\mathrm{GU\_\_LOSEFL}} &= \mathrm{GU}_{\mathrm{\_f\_abcb1}} \cdot \mathrm{GU}_{\mathrm{\_LOSEFL\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{GU}_{\mathrm{\_los}} \\[1ex]
v_{\mathrm{GU\_\_LOSEX}} &= \mathrm{GU}_{\mathrm{\_LOSABS\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{GU}_{\mathrm{\_los}} \\[1ex]
v_{\mathrm{GU\_\_LOSEXC}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{Cgu}_{\mathrm{los}} \\[1ex]
v_{\mathrm{GU\_\_E3174EXC}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{Cgu}_{\mathrm{e3174}} \\[1ex]
v_{\mathrm{GU\_\_L158EXC}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{Vgu} \cdot \mathrm{Cgu}_{\mathrm{l158}} \\[1ex]
v_{\mathrm{GU\_\_dissolution\_los}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_dis\_los}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{los}}}{\mathrm{Mr}_{\mathrm{los}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{ald}}{\mathrm{d} t} &= \frac{v_{\mathrm{ALDSEC}} - v_{\mathrm{ALDDEG}}}{\mathrm{Vplasma}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{ren}}{\mathrm{d} t} &= \frac{v_{\mathrm{RENSEC}} - v_{\mathrm{RENDEG}}}{\mathrm{Vplasma}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{ang}_{1}}{\mathrm{d} t} &= \frac{v_{\mathrm{ANGGEN2ANG1}} - v_{\mathrm{ANG1ANG2}}}{\mathrm{Vplasma}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{ang}_{2}}{\mathrm{d} t} &= \frac{v_{\mathrm{ANG1ANG2}} - v_{\mathrm{ANG2DEG}}}{\mathrm{Vplasma}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{blood\_los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_los}} - v_{\mathrm{Flow\_ki\_ve\_los}} - v_{\mathrm{KI\_\_LOSEX}}}{\mathrm{Vki}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{blood\_los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_los}} + v_{\mathrm{Flow\_po\_li\_los}} - v_{\mathrm{Flow\_li\_hv\_los}} - v_{\mathrm{LI\_\_LOSIM}}}{\mathrm{Vli}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{blood\_los}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_lu\_los}} + v_{\mathrm{Flow\_ve\_lu\_los}} - v_{\mathrm{Flow\_lu\_ar\_los}}}{\mathrm{Vlu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{blood\_los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_los}} - v_{\mathrm{Flow\_gu\_po\_los}} + v_{\mathrm{GU\_\_LOSEX}}}{\mathrm{Vgu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{blood\_los}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_re\_los}} + v_{\mathrm{Flow\_ar\_re\_los}} - v_{\mathrm{Flow\_re\_ve\_los}}}{\mathrm{Vre}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{los}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_los}} - v_{\mathrm{Flow\_arli\_li\_los}} - v_{\mathrm{Flow\_arli\_hv\_los}} + v_{\mathrm{Flow\_lu\_ar\_los}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_los}} - v_{\mathrm{Flow\_ar\_re\_los}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{los}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vve}} \Bigl( v_{\mathrm{iv\_los}} + v_{\mathrm{Flow\_ki\_ve\_los}} + v_{\mathrm{Flow\_hv\_ve\_los}} - v_{\mathrm{Flow\_ve\_lu\_los}} \\
&\quad {} + v_{\mathrm{Flow\_re\_ve\_los}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{los}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_los}} - v_{\mathrm{Flow\_po\_hv\_los}} + v_{\mathrm{Flow\_gu\_po\_los}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_los}} + v_{\mathrm{Flow\_po\_hv\_los}} + v_{\mathrm{Flow\_li\_hv\_los}} - v_{\mathrm{Flow\_hv\_ve\_los}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_lu\_los}}}{\mathrm{Vlu}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_re\_los}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{los}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_LOSEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{los}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_LOSEXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{blood\_e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_e3174}} - v_{\mathrm{Flow\_ki\_ve\_e3174}} - v_{\mathrm{KI\_\_E3174EX}}}{\mathrm{Vki}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{blood\_e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_e3174}} + v_{\mathrm{Flow\_po\_li\_e3174}} - v_{\mathrm{Flow\_li\_hv\_e3174}} + v_{\mathrm{LI\_\_E3174EX}}}{\mathrm{Vli}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{blood\_e3174}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_lu\_e3174}} + v_{\mathrm{Flow\_ve\_lu\_e3174}} - v_{\mathrm{Flow\_lu\_ar\_e3174}}}{\mathrm{Vlu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{blood\_e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_e3174}} - v_{\mathrm{Flow\_gu\_po\_e3174}}}{\mathrm{Vgu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{blood\_e3174}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_re\_e3174}} + v_{\mathrm{Flow\_ar\_re\_e3174}} - v_{\mathrm{Flow\_re\_ve\_e3174}}}{\mathrm{Vre}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{e3174}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_e3174}} - v_{\mathrm{Flow\_arli\_li\_e3174}} - v_{\mathrm{Flow\_arli\_hv\_e3174}} + v_{\mathrm{Flow\_lu\_ar\_e3174}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_e3174}} - v_{\mathrm{Flow\_ar\_re\_e3174}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{e3174}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vve}} \Bigl( v_{\mathrm{iv\_e3174}} + v_{\mathrm{Flow\_ki\_ve\_e3174}} + v_{\mathrm{Flow\_hv\_ve\_e3174}} - v_{\mathrm{Flow\_ve\_lu\_e3174}} \\
&\quad {} + v_{\mathrm{Flow\_re\_ve\_e3174}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{e3174}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_e3174}} - v_{\mathrm{Flow\_po\_hv\_e3174}} + v_{\mathrm{Flow\_gu\_po\_e3174}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_e3174}} + v_{\mathrm{Flow\_po\_hv\_e3174}} + v_{\mathrm{Flow\_li\_hv\_e3174}} - v_{\mathrm{Flow\_hv\_ve\_e3174}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_lu\_e3174}}}{\mathrm{Vlu}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_re\_e3174}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{e3174}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_E3174EX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{e3174}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_E3174EXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{blood\_l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_l158}} - v_{\mathrm{Flow\_ki\_ve\_l158}} - v_{\mathrm{KI\_\_L158EX}}}{\mathrm{Vki}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{blood\_l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_l158}} + v_{\mathrm{Flow\_po\_li\_l158}} - v_{\mathrm{Flow\_li\_hv\_l158}} + v_{\mathrm{LI\_\_L158EX}}}{\mathrm{Vli}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{blood\_l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ve\_lu\_l158}} - v_{\mathrm{Flow\_lu\_ar\_l158}}}{\mathrm{Vlu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{blood\_l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_l158}} - v_{\mathrm{Flow\_gu\_po\_l158}}}{\mathrm{Vgu}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{blood\_l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_re\_l158}} - v_{\mathrm{Flow\_re\_ve\_l158}}}{\mathrm{Vre}_{\mathrm{blood}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{l158}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_l158}} - v_{\mathrm{Flow\_arli\_li\_l158}} - v_{\mathrm{Flow\_arli\_hv\_l158}} + v_{\mathrm{Flow\_lu\_ar\_l158}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_l158}} - v_{\mathrm{Flow\_ar\_re\_l158}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ki\_ve\_l158}} + v_{\mathrm{Flow\_hv\_ve\_l158}} - v_{\mathrm{Flow\_ve\_lu\_l158}} + v_{\mathrm{Flow\_re\_ve\_l158}}}{\mathrm{Vve}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{l158}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_l158}} - v_{\mathrm{Flow\_po\_hv\_l158}} + v_{\mathrm{Flow\_gu\_po\_l158}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_l158}} + v_{\mathrm{Flow\_po\_hv\_l158}} + v_{\mathrm{Flow\_li\_hv\_l158}} - v_{\mathrm{Flow\_hv\_ve\_l158}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{l158}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_L158EX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{l158}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_L158EXC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{los}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{GU\_\_LOSABS}} + v_{\mathrm{GU\_\_LOSEFL}} - v_{\mathrm{GU\_\_LOSEXC}} + v_{\mathrm{GU\_\_dissolution\_los}}}{\mathrm{Vgu}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_E3174EHC}} - v_{\mathrm{GU\_\_E3174EXC}}}{\mathrm{Vgu}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_L158EHC}} - v_{\mathrm{GU\_\_L158EXC}}}{\mathrm{Vgu}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_LOSIM}} - v_{\mathrm{LI\_\_LOS2E3174}} - v_{\mathrm{LI\_\_LOS2L158}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_e3174}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_LOS2E3174}} - v_{\mathrm{LI\_\_E3174EX}} - v_{\mathrm{LI\_\_E3174BIEX}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_e3174\_bi}}}{\mathrm{d} t} &= v_{\mathrm{LI\_\_E3174BIEX}} - v_{\mathrm{LI\_\_E3174EHC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_l158}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_LOS2L158}} - v_{\mathrm{LI\_\_L158EX}} - v_{\mathrm{LI\_\_L158BIEX}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_l158\_bi}}}{\mathrm{d} t} &= v_{\mathrm{LI\_\_L158BIEX}} - v_{\mathrm{LI\_\_L158EHC}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_los}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_LOSABS}} - v_{\mathrm{GU\_\_LOSEFL}} - v_{\mathrm{GU\_\_LOSEX}}}{\mathrm{GU}_{\mathrm{\_Ventero}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{los}}}{\mathrm{d} t} &= -v_{\mathrm{GU\_\_dissolution\_los}} \cdot \mathrm{Mr}_{\mathrm{los}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{IVDOSE}_{\mathrm{los}}}{\mathrm{d} t} &= -v_{\mathrm{iv\_los}} \cdot \mathrm{Mr}_{\mathrm{los}} + \mathrm{Ri}_{\mathrm{los}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{cum}_{\mathrm{dose\_los}}}{\mathrm{d} t} &= \mathrm{Ri}_{\mathrm{los}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{IVDOSE}_{\mathrm{e3174}}}{\mathrm{d} t} &= -v_{\mathrm{iv\_e3174}} \cdot \mathrm{Mr}_{\mathrm{e3174}} + \mathrm{Ri}_{\mathrm{e3174}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{cum}_{\mathrm{dose\_e3174}}}{\mathrm{d} t} &= \mathrm{Ri}_{\mathrm{e3174}} \qquad \text{(rate rule)}
\end{aligned}
$$
