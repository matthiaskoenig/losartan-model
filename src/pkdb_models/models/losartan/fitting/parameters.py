"""Parameters to optimize in the losartan fits.

The start values were set by hand from the literature
(Lo1995, FDA1995S60): plasma clearance of losartan 610 ml/min (renal 70),
Vss 34 l, half-life 2.1 hr, bioavailability 33 %, conversion to E3174 14 %;
E3174 clearance 47 ml/min (renal 26), Vss 10 l, half-life 6-9 hr. The bounds
contain these values.
"""

from sbmlsim.fit import FitParameter

#: the tissue distribution of losartan and E3174 (Vss 34 l and 10 l [Lo1995]).
#: The E3174 of the tissues is a slow compartment: its terminal half-life is
#: 6.3 hr after the infusion of E3174 and 6-9 hr after oral losartan
#: [Lo1995]. Without the bounds of a slow exchange (`ftissue_e3174` up to
#: 0.1 l/min) and of a tissue store (`Kp_e3174` from 0.05), the first fit put
#: the E3174 of the tissues to nothing (Kp 0.014, ftissue not identifiable):
#: Vss 7.8 l, a terminal half-life of 3.8 hr and the E3174 at 24-30 hr after
#: the infusion 2-3 fold below the data, which weigh little in the cost
parameters_pk_distribution: list[FitParameter] = [
    FitParameter(
        pid="ftissue_los",
        lower_bound=0.05,
        start_value=0.6,
        upper_bound=5,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_los",
        lower_bound=0.1,
        start_value=0.9,
        upper_bound=5,
        unit="dimensionless",
    ),
    FitParameter(
        pid="ftissue_e3174",
        lower_bound=0.001,
        start_value=0.03,
        upper_bound=0.1,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_e3174",
        lower_bound=0.05,
        start_value=0.2,
        upper_bound=2,
        unit="dimensionless",
    ),
]

#: the absorption of losartan (tmax about 1 hr). The fraction absorbed is
#: nearly complete (absorption against the colonic transit `GU__METEXC_k`,
#: which is not fitted), the efflux (`GU__f_LOSEFL_k`) is not identifiable
#: separately from the absorption and is not fitted
parameters_pk_absorption: list[FitParameter] = [
    FitParameter(
        pid="GU__LOSABS_k",
        lower_bound=1e-3,
        start_value=0.05,
        upper_bound=0.5,
        unit="1/min",
    ),
]

#: the hepatic metabolism of losartan, the hepatic transport and the biliary
#: excretion of E3174 and the other metabolites
parameters_pk_liver: list[FitParameter] = [
    FitParameter(
        pid="LI__LOS2E3174_k",
        lower_bound=0.01,
        start_value=0.45,
        upper_bound=10,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__LOS2L158_k",
        lower_bound=0.1,
        start_value=1.6,
        upper_bound=20,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__E3174EX_k",
        lower_bound=1e-4,
        start_value=0.006,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__E3174IM_k",
        lower_bound=1e-3,
        start_value=0.08,
        upper_bound=10,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__E3174BIEX_k",
        lower_bound=1e-5,
        start_value=0.0025,
        upper_bound=0.1,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__L158BIEX_k",
        lower_bound=1e-3,
        start_value=0.045,
        upper_bound=1,
        unit="1/min",
    ),
]

#: the renal excretion relative to the GFR (renal clearance of plasma of
#: losartan 70 ml/min, E3174 26 ml/min [Lo1995])
parameters_pk_kidney: list[FitParameter] = [
    FitParameter(
        pid="KI__LOSEX_f",
        lower_bound=5,
        start_value=59,
        upper_bound=500,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__E3174EX_f",
        lower_bound=10,
        start_value=129,
        upper_bound=1000,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__L158EX_f",
        lower_bound=0.05,
        start_value=0.5,
        upper_bound=10,
        unit="dimensionless",
    ),
]

#: the parameters of the pharmacokinetics
parameters_pk: list[FitParameter] = [
    *parameters_pk_distribution,
    *parameters_pk_absorption,
    *parameters_pk_liver,
    *parameters_pk_kidney,
]

#: the parameters of the pharmacodynamics, fitted on the adopted
#: pharmacokinetics: the AT1 blockade by E3174, the feedback on the renin
#: secretion and the effect on the blood pressure. The turnover of the RAAS
#: (`ren_t12`, `ald_t12`, `ang1_t12`) is the literature and is not fitted.
parameters_pd: list[FitParameter] = [
    FitParameter(
        pid="E50_e3174",
        lower_bound=1e-5,
        start_value=2e-4,
        upper_bound=1e-2,
        unit="mM",
    ),
    FitParameter(
        pid="RENSEC_n",
        lower_bound=0.2,
        start_value=2.0,
        upper_bound=10,
        unit="dimensionless",
    ),
    FitParameter(
        pid="BP_emax",
        lower_bound=0.01,
        start_value=0.25,
        upper_bound=0.6,
        unit="dimensionless",
    ),
]
