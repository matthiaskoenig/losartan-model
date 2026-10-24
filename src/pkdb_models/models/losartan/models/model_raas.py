"""Model of RAAS blood pressure regulation."""

from dataclasses import dataclass

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    FactoryResult,
    InitialAssignment,
    Model,
    Parameter,
    Reaction,
    Species,
    UnitDefinition,
    ValidationOptions,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.losartan.models import annotations, templates
from pkdb_models.models.losartan.models.fits import FITS


class U(templates.U):
    """Units of the RAAS model."""

    mmHg = UnitDefinition("mmHg", "133.32239 N/m^2")


mid = "losartan_raas"
version = 3

_m = Model(
    sid=mid,
    name="Model for RAAS system of blood pressure regulation.",
    notes=f"""
    Model for RAAS system of blood pressure regulation.

    **version** {version}

    Renin cleaves angiotensinogen to angiotensin I, ACE converts angiotensin I
    to angiotensin II. Angiotensin II acts on the AT1 receptors, which E3174
    blocks (`fe_e3174`, the fraction of blocked receptors): the AT1 signal
    `f_at1` is the angiotensin II relative to its baseline times the
    fraction of free receptors. The AT1 signal drives the secretion of
    aldosterone, inhibits the secretion of renin (short loop feedback,
    `RENSEC_n`) and raises the blood pressure (`BP_emax`). The turnover of
    renin, aldosterone and the angiotensins is given by their half-lives in
    plasma (renin 15 min, aldosterone 20 min, angiotensin I 0.5 min and
    angiotensin II 0.4 min); the concentrations do not depend on the volume
    of the compartment. The baselines (`*_ref`) are steady states, the
    initial concentrations are the baselines (initial assignments), which the
    experiments set per study.

    ## Changelog
    **version 3**

    - turnover from the half-lives, angiotensinogen 1 µM
    - AT1 signal of angiotensin II and E3174: aldosterone secretion, renin
      feedback and blood pressure

    **version 2**

    - better handling of baselines and scaling

    **version 1**

    - initial model

    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        (BQB.HAS_PROPERTY, "NCIT:C15720"),  # pharmacodynamics
        (BQB.HAS_PROPERTY, "NCIT:C91434"),  # Renin-Angiotensin Pathway
        *annotations.model,
    ],
)

_m.compartments = [
    Compartment(
        "Vplasma",
        value=3.0,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True,
        notes="""
        The plasma of the RAAS, about 3 l. The rates scale with the volume,
        the concentrations do not depend on it.
        """,
    ),
]


@dataclass
class Substance:
    """Definition of a substance of the RAAS model."""

    sid: str
    name: str
    init: float  # initial concentration [mM]
    unit: str
    boundary: bool = False
    port: bool = False


# ---------------------------------------------------------------------------------------------------------------------
# Pharmacodynamics
# ---------------------------------------------------------------------------------------------------------------------
ald_ref = 250e-9  # [mM]
ren_ref = 1e-9  # [mM]
anggen_ref = 1e-3  # [mM], 1 µM (25-60 µg/ml, Mr 52670)
ang1_ref = 10e-9  # [mM]
ang2_ref = 8e-9  # [mM]

substances: list[Substance] = [
    # Mr: 360.444 [g/mole]
    # 100 pg/ml {ramipril/Bussien1985} = 100/360.444*1000 [pmol/l] = 277.4 [pmol/l]
    # 80-90 pg/ml {losartan/Christen1991a} = 221.9 - 249.7 [pmol/l]
    # 80-100 pg/ml {losartan/Ohtawa1993} = 221.9 - 277.4 [pmol/l]
    # 754 ± 351 pg/ml {losartan/Doig1993} (salt depletion) = 2092 [pmol/l]
    Substance("ald", "aldosterone", init=ald_ref, unit="mmole"),
    # Mr: 45057 [g/mole]
    # ~ 50-65 pg/ml {losartan/Azizi1999} = 50-65/45.057 [pmol/l] = 1.11 - 1.44 [pmol/l]
    # ~ 45-75 pg/ml {losartan/Doig1993}
    # the plasma renin activity (PRA, ng angiotensin I/ml/hr: Goldberg1995,
    # Ohtawa1993, Christen1991a, Doig1993) is not a renin concentration, it
    # is compared as the ratio to its baseline (`ren_ratio`)
    Substance("ren", "renin", init=ren_ref, unit="mmole"),
    # Mr: 52670 [g/mole]
    Substance(
        "anggen", "angiotensinogen", init=anggen_ref, unit="mmole", boundary=True
    ),
    # Mr: 1296.5 [g/mole]
    # ~ 11-13 [pg/ml] {losartan/Azizi1999} = 11-13 / 1.2965 [pmol/l] = 8.5 - 10.0 [pmol/l]
    # 10 [pmol/l] {ramipril/Manhem1985}
    Substance("ang1", "angiotensin I", init=ang1_ref, unit="mmole"),
    # Mr: 1046.2 [g/mole]
    # ~ 6-10 [pg/ml] {losartan/Azizi1999} = 6-10 / 1.0462 [pmol/l] = 5.7 - 9.6 [pmol/l]
    # 10 [pmol/l] {ramipril/Manhem1985}
    # 5-6 [fmol/ml] {losartan/Christen1991a} = 5-6 [pmol/l]
    # 4-8 [pg/ml] {losartan/Ohtawa1993} = 4-8 / 1.0462 [pmol/l] = 3.8 - 9.2 [pmol/l]
    Substance("ang2", "angiotensin II", init=ang2_ref, unit="mmole"),
    # ACE is not modelled; Mr: 149715 [g/mole]
    # 10 pmol/l: assumption, same order as ang1 and ang and renin
    Substance("e3174", "e3174", init=0, unit="mmole", port=True),
]

for s in substances:
    _m.species.append(
        Species(
            s.sid,
            name=s.name,
            initialConcentration=s.init,
            compartment="Vplasma",
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,
            boundaryCondition=s.boundary,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[s.sid],
            port=s.port,
        )
    )
    if s.sid not in {"e3174"}:
        _m.assignments.append(
            # the initial concentration is the baseline of the study
            InitialAssignment(s.sid, f"{s.sid}_ref", unit=U.mM)
        )

    # changes in variables
    if s.sid in ["ang1", "ang2", "ren", "ald"]:
        # absolute change
        _m.parameters.append(
            Parameter(
                f"{s.sid}_change",
                name=f"{s.name} change",
                value=np.nan,
                unit=U.mM,
                annotations=annotations.species[s.sid],
                constant=False,
                notes=f"Absolute change to baseline {s.name}",
            )
        )
        _m.rules.append(
            AssignmentRule(f"{s.sid}_change", f"{s.sid}-{s.sid}_ref", unit=U.mM)
        )
        # ratio to baseline
        _m.parameters.append(
            Parameter(
                f"{s.sid}_ratio",
                name=f"{s.name} ratio",
                value=np.nan,
                unit=U.dimensionless,
                annotations=annotations.species[s.sid],
                constant=False,
                notes=f"Ratio relative to baseline {s.name}",
            )
        )
        _m.rules.append(
            AssignmentRule(
                f"{s.sid}_ratio", f"{s.sid}/{s.sid}_ref", unit=U.dimensionless
            )
        )

_m.parameters.extend(
    [
        Parameter(
            "anggen_ref",
            anggen_ref,
            U.mM,
            name="reference concentration of angiotensinogen",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "ang1_ref",
            ang1_ref,
            U.mM,
            name="reference concentration of angiotensin I",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "ang2_ref",
            ang2_ref,
            U.mM,
            name="reference concentration of angiotensin II",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "ren_ref",
            ren_ref,
            U.mM,
            name="reference concentration of renin",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "ald_ref",
            ald_ref,
            U.mM,
            name="reference concentration of aldosterone",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        # pharmacodynamics of exp3174
        Parameter(
            sid="E50_e3174",
            name="half-maximum effect E3174",  # Mr = 436.894 g/mol
            # ~ 32 [1-99] ng/ml {Munafo1992} = 32/436.894*1000 [nmol/l] = 73.24 [2.29 - 226.6] [nmol/l]
            **FITS.parameter("E50_e3174", unit=U.mM),
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The E3174 of plasma at which half of the AT1 receptors are blocked.
            Fitted to the renin, angiotensin II and aldosterone of the PD
            data: a blockade of 0.83 at the peak and 0.09 at 24 hr after
            100 mg losartan. The inhibition of the pressor response to
            angiotensin II is about 0.85 at the peak and 0.25-0.40 at 24 hr
            [FDA label]; the E3174 of the model falls faster than measured
            (terminal half-life 3.4 hr against 6.3 hr [Lo1995]).
            """,
        ),
        Parameter(
            "ren_t12",
            15,
            U.min,
            name="half-life of renin",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="10-15 min (active renin in plasma).",
        ),
        Parameter(
            "ald_t12",
            20,
            U.min,
            name="half-life of aldosterone",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="15-20 min.",
        ),
        Parameter(
            "ang1_t12",
            0.5,
            U.min,
            name="half-life of angiotensin I",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Seconds to a minute; it sets the flux through the angiotensins,
            the half-life of angiotensin II is ang1_t12 * ang2_ref/ang1_ref.
            """,
        ),
        Parameter(
            "RENSEC_n",
            **FITS.parameter("RENSEC_n", unit=U.dimensionless),
            name="feedback of the AT1 signal on the renin secretion",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The renin secretion is the AT1 signal to the power -RENSEC_n (the
            short loop feedback of angiotensin II on the juxtaglomerular
            cells). In the steady state renin and angiotensin II rise by
            (1 - fe)^(-n/(1+n)) and the AT1 signal and aldosterone fall to
            (1 - fe)^(1/(1+n)): PRA x2-3 and aldosterone -20 to -40 % at the
            peak of 50-100 mg losartan.
            """,
        ),
    ]
)
_m.rules.extend(
    [
        AssignmentRule(
            "fe_e3174",
            "e3174/(E50_e3174 + e3174)",
            unit=U.dimensionless,
            name="AT1 receptor blockade by E3174 in [0, 1]",
        ),
        AssignmentRule(
            "f_at1",
            "ang2/ang2_ref * (1 dimensionless - fe_e3174)",
            unit=U.dimensionless,
            name="AT1 signal relative to the baseline",
            notes="""
            Angiotensin II relative to its baseline times the fraction of the
            AT1 receptors which E3174 does not block; 1 at the baseline.
            """,
        ),
    ]
)


_m.reactions = [
    # renin turnover
    Reaction(
        sid="RENSEC",
        name="renin secretion (RENSEC)",
        equation="-> ren [ang2, e3174]",
        compartment="Vplasma",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RENSEC_k",
                "RENDEG_k * ren_ref",
                U.mmole_per_min,
                name="rate renin secretion",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
                notes="The steady state of the baseline renin.",
            ),
        ],
        formula=("RENSEC_k * power(f_at1, -RENSEC_n)", U.mmole_per_min),
        notes="""
        Angiotensin II inhibits the secretion of renin via AT1 (short loop
        feedback); the blockade of AT1 by E3174 removes the inhibition.
        """,
    ),
    Reaction(
        sid="RENDEG",
        name="renin degradation (RENDEG)",
        equation="ren ->",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "RENDEG_k",
                "ln(2 dimensionless)/ren_t12 * Vplasma",
                U.l_per_min,
                name="rate renin degradation",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            )
        ],
        formula=("RENDEG_k * ren", U.mmole_per_min),
    ),
    # aldosterone turnover
    Reaction(
        sid="ALDSEC",
        name="aldosterone secretion (ALDSEC)",
        equation="-> ald [ang2, e3174]",
        compartment="Vplasma",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "ALDSEC_k",
                "ALDDEG_k * ald_ref",
                U.mmole_per_min,
                name="rate aldosterone secretion",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
                notes="The steady state of the baseline aldosterone.",
            )
        ],
        formula=("ALDSEC_k * f_at1", U.mmole_per_min),
        notes="""
        Angiotensin II drives the secretion of aldosterone via AT1 receptors
        of the adrenal cortex; E3174 blocks the receptors.
        """,
    ),
    Reaction(
        sid="ALDDEG",
        name="aldosterone degradation (ALDDEG)",
        equation="ald ->",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "ALDDEG_k",
                "ln(2 dimensionless)/ald_t12 * Vplasma",
                U.l_per_min,
                name="rate aldosterone degradation",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            )
        ],
        formula=("ALDDEG_k * ald", U.mmole_per_min),
    ),
    # anggen -> ang1 -> ang2 ->
    Reaction(
        sid="ANGGEN2ANG1",
        name="angiotensinogen to angiotensin I (renin)",
        equation="anggen -> ang1 [ren]",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "ANGGEN2ANG1_k",
                "ANG1ANG2_k * ang1_ref/anggen_ref",
                U.l_per_min,
                name="rate angen → ang1",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            ),
        ],
        formula=("ANGGEN2ANG1_k * anggen * ren/ren_ref", U.mmole_per_min),
    ),
    Reaction(
        sid="ANG1ANG2",
        name="angiotensin I to angiotensin II (ACE)",
        equation="ang1 -> ang2",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "ANG1ANG2_k",
                "ln(2 dimensionless)/ang1_t12 * Vplasma",
                U.l_per_min,
                name="rate ang1 → ang2",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            ),
        ],
        formula=("ANG1ANG2_k * ang1", U.mmole_per_min),
    ),
    Reaction(
        sid="ANG2DEG",
        name="angiotensin II degradation (ANG2DEG)",
        equation="ang2 ->",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "ANG2DEG_k",
                "ANG1ANG2_k * ang1_ref/ang2_ref",
                U.l_per_min,
                name="rate angiotensin II degradation",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            )
        ],
        formula=("ANG2DEG_k * ang2", U.mmole_per_min),
    ),
]

# blood pressure model
_m.parameters.extend(
    [
        Parameter(
            "SBP_ref",
            120,
            U.mmHg,
            name="reference systolic blood pressure",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "DBP_ref",
            80,
            U.mmHg,
            name="reference diastolic blood pressure",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "SBP",
            np.nan,
            U.mmHg,
            name="systolic blood pressure",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            constant=False,
        ),
        Parameter(
            "DBP",
            np.nan,
            U.mmHg,
            name="diastolic blood pressure",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            constant=False,
        ),
        Parameter(
            "MAP",
            np.nan,
            U.mmHg,
            constant=False,
            name="mean arterial pressure",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "BP_emax",
            **FITS.parameter("BP_emax", unit=U.dimensionless),
            name="maximal fractional fall of the blood pressure",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The fall of the blood pressure relative to its baseline at the
            complete loss of the AT1 signal (vasoconstriction, sodium
            retention). Fitted, but not identified by the blood pressure of
            the data: -7 mmHg of the diastolic blood pressure at the peak of
            100 mg losartan (-6 to -10 mmHg in hypertension [FDA label]).
            """,
        ),
    ]
)
_m.rules.extend(
    [
        AssignmentRule(
            "SBP",
            "SBP_ref * (1 dimensionless - BP_emax * (1 dimensionless - f_at1))",
            unit=U.mmHg,
        ),
        AssignmentRule(
            "DBP",
            "DBP_ref * (1 dimensionless - BP_emax * (1 dimensionless - f_at1))",
            unit=U.mmHg,
        ),
        AssignmentRule("MAP", "DBP + (SBP - DBP)/3 dimensionless", unit=U.mmHg),
    ]
)

# blood pressure changes
for bp_sid, bp_name in [
    ("SBP", "systolic blood pressure"),
    ("DBP", "diastolic blood pressure"),
]:
    # absolute change
    _m.parameters.append(
        Parameter(
            f"{bp_sid}_change",
            name=f"{bp_name} change",
            value=np.nan,
            unit=U.mmHg,
            constant=False,
            notes=f"Absolute change to baseline {bp_name}",
        )
    )
    _m.rules.append(
        AssignmentRule(f"{bp_sid}_change", f"{bp_sid}-{bp_sid}_ref", unit=U.mmHg)
    )
    # relative change
    _m.parameters.append(
        Parameter(
            f"{bp_sid}_ratio",
            name=f"{bp_name} ratio",
            value=np.nan,
            unit=U.dimensionless,
            constant=False,
            notes=f"Ratio relative to baseline {bp_name}",
        )
    )
    _m.rules.append(
        AssignmentRule(
            f"{bp_sid}_ratio", f"{bp_sid}/{bp_sid}_ref", unit=U.dimensionless
        )
    )


model_raas = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.losartan import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_raas,
        filepath=MODEL_BASE_PATH / f"{model_raas.sid}.xml",
        sbml_level=3,
        sbml_version=2,
        validation_options=ValidationOptions(units_consistency=True),
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
