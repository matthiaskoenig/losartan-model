"""FitParameters for losartan."""
import copy

from sbmlsim.fit import FitParameter


parameters_pk = [
    # tissue distribution
    FitParameter(
        pid="ftissue_los",
        lower_bound=0.01,
        start_value=0.1,
        upper_bound=10,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_los",
        lower_bound=1,
        start_value=20,
        upper_bound=200,
        unit="dimensionless",
    ),
    # FitParameter(
    #     pid="ftissue_e3174",
    #     lower_bound=0.01,
    #     start_value=0.1,
    #     upper_bound=50,
    #     unit="l/min",
    # ),
    # FitParameter(
    #     pid="Kp_e3174",
    #     lower_bound=1,
    #     start_value=20,
    #     upper_bound=200,
    #     unit="dimensionless",
    # ),

    # absorption
    # FitParameter(
    #     pid="GU__F_los_abs",
    #     lower_bound=0.5,
    #     start_value=0.6,
    #     upper_bound=0.7,
    #     unit="dimensionless",
    # ),
    FitParameter(
        pid="GU__LOSABS_k",
        lower_bound=1E-4,
        start_value=0.02,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter( # efflux relative to absorption
        pid="GU__f_LOSEFL_k",
        lower_bound=0.1,
        start_value=1,
        upper_bound=10,
        unit="dimensionless",
    ),
    FitParameter(  # excretion feces
        pid="GU__METEXC_k",
        lower_bound=1E-5,
        start_value=0.01,
        upper_bound=1E-1,
        unit="1/min",
    ),

    # hepatic metabolism


    FitParameter( # e3174 transport
        pid="LI__E3174EX_k",
        lower_bound=1E-3,
        start_value=1.0,
        upper_bound=10,
        unit="1/min",
    ),

    FitParameter( # cyp2c9 activity
        pid="LI__LOS2E3174_Vmax",
        lower_bound=1E-5,
        start_value=1.0,
        upper_bound=100,
        unit="mmol/min/l",
    ),
    FitParameter(
        pid="LI__E3174L158_k",
        lower_bound=1E-5,
        start_value=1.0,
        upper_bound=10,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__MBIEX_k",
        lower_bound=1E-5,
        start_value=1E-4,
        upper_bound=1,
        unit="1/min",
    ),

    # kidney removal
    FitParameter(
        pid="KI__LOSEX_k",
        lower_bound=1E-4,
        start_value=1E-1,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="KI__E3174EX_k",
        lower_bound=1E-4,
        start_value=1E-1,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="KI__L158EX_k",
        lower_bound=1E-4,
        start_value=1E-1,
        upper_bound=1,
        unit="1/min",
    ),
]

#  --- Pharmacodynamics ---
parameters_pd = [
    # FitParameter(
    #     pid="RENSEC_k",
    #     lower_bound=1E-6,
    #     start_value=0.1,
    #     upper_bound=1,
    #     unit="mmole/min",
    # ),
    FitParameter(
        pid="ANGGEN2ANG1_k",
        lower_bound=1E-2,
        start_value=100,
        upper_bound=1E5,
        unit="l/min",
    ),
    FitParameter(
        pid="E50_e3174",
        lower_bound=5E-7,
        start_value=1E-3,
        upper_bound=5E-2,
        unit="mM",
    ),
    FitParameter(
        pid="ALDSEC_k",
        lower_bound=1E-6,
        start_value=0.1,
        upper_bound=1,
        unit="mmole/min",
    ),
    FitParameter(
        pid="BP_ald_fe",
        lower_bound=0.1,
        start_value=0.2,
        upper_bound=0.6,
        unit="dimensionless",
    ),

    # Parameter(
    #     sid="Emax_e3174",
    #     name="maximum effect of exp3174",
    #     value=0.91,  # 0.91 [0.66 - 1.00] {Munafo1992}
    #     unit=U.dimensionless,
    #     sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    # ),

    # Parameter(
    #     sid="Egamma_e3174",
    #     name="hill coefficient exp3174",
    #     value=0.9,  # 0.9 [0.13 - 3] dimensionless {Munafo1992}
    #     unit=U.dimensionless,
    #     sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    # )

]

parameters_control = parameters_pk + parameters_pd
parameters_all = copy.deepcopy(parameters_control)
