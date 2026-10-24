"""Kidney model for losartan.

The renal clearances are relative to the glomerular filtration rate (GFR):

- losartan: 70-72 ml/min of plasma, 12 % of its clearance, fu 1.3 % (net
  tubular secretion) [Lo1995]; 4-6 % of an oral dose in the urine.
- E3174: 26 ml/min of plasma, 55 % of its clearance, fu 0.2 % [Lo1995]; 6 % of
  an oral dose of losartan in the urine.
- the other metabolites: about 18 % of an intravenous dose in the urine
  [FDA1995S60].
"""

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    FactoryResult,
    Model,
    Parameter,
    Reaction,
    Species,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.losartan.models import annotations, templates
from pkdb_models.models.losartan.models.fits import FITS


class U(templates.U):
    """Units of the kidney model."""

    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    l_per_ml = UnitDefinition("l_per_ml", "l/ml")


mid = "losartan_kidney"
version = 2

_m = Model(
    sid=mid,
    name="Model for renal excretion of losartan, E3174, and the other metabolites.",
    notes=f"""
    Model for renal excretion of losartan, E3174, and the other metabolites.

    **version** {version}

    The renal clearances are given relative to the glomerular filtration rate
    (GFR) and scale with the renal function (`f_renal_function`): the renal
    clearance of plasma of losartan and E3174 is `<X>EX_f` times the
    filtration of their unbound fraction in plasma (`fu_<x> * GFR`), the one of
    the other metabolites `L158EX_f` times the GFR. The concentrations of the
    kidney model are blood concentrations, the plasma concentration is the
    blood concentration over the blood to plasma ratio `BP_<x>`.

    ## Changelog

    **version 2**

    - renal clearances relative to the GFR, plasma by the blood to plasma ratio

    **version 1**

    - initial model

    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),
        (BQB.OCCURS_IN, "ncit/C12415"),
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
        *annotations.model,
    ],
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="blood",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["blood"],
        port=True,
    ),
    Compartment(
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True,
    ),
    Compartment(
        "Vmem",
        value=np.nan,
        unit=U.m2,
        name="plasma membrane",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma membrane"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vurine",
        1.0,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["urine"],
    ),
]

# ---------------------------------------------------------------------------------------------------------------------
# Pharmacokinetics
# ---------------------------------------------------------------------------------------------------------------------
_m.species = [
    Species(
        "los_ext",
        name="losartan (blood)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        port=True,
    ),
    Species(
        "los_urine",
        name="losartan (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,  # this is in amounts
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        port=True,
    ),
    Species(
        "e3174_ext",
        name="E3174 (blood)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["e3174"],
        port=True,
    ),
    Species(
        "e3174_urine",
        name="E3174 (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,  # this is in amounts
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["e3174"],
        port=True,
    ),
    Species(
        "l158_ext",
        name="other losartan metabolites (blood)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["l158"],
        port=True,
    ),
    Species(
        "l158_urine",
        name="other losartan metabolites (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,  # this is in amounts
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["l158"],
        port=True,
    ),
]

_m.parameters.extend(
    [
        Parameter(
            "f_renal_function",
            name="renal function",
            value=1.0,
            unit=U.dimensionless,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The GFR relative to the healthy GFR `GFR_base` (101 ml/min at
            75 kg): a measured creatinine clearance of X ml/min is the renal
            function X/101. 1.0 is the normal renal function, below 1.0 an
            impaired one.
            """,
        ),
        Parameter(
            "GFR_base",
            101,
            U.ml_per_min,
            name="glomerular filtration rate (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
        ),
        Parameter(
            "fu_los",
            0.013,
            U.dimensionless,
            name="fraction unbound in plasma losartan",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="1.3 % [Lo1995, FDA label]",
        ),
        Parameter(
            "fu_e3174",
            0.002,
            U.dimensionless,
            name="fraction unbound in plasma E3174",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="0.2 % [Lo1995, FDA label]",
        ),
        Parameter(
            "BP_los",
            0.53,
            U.dimensionless,
            name="blood to plasma ratio losartan",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
        ),
        Parameter(
            "BP_e3174",
            0.58,
            U.dimensionless,
            name="blood to plasma ratio E3174",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
        ),
        Parameter(
            "BP_l158",
            1.0,
            U.dimensionless,
            name="blood to plasma ratio other metabolites",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
        ),
        Parameter(
            "LOSEX_f",
            **FITS.parameter("KI__LOSEX_f", unit=U.dimensionless),
            name="renal clearance of unbound losartan relative to the GFR",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of losartan is LOSEX_f times the filtration of
            its unbound fraction, fu_los * GFR: 1 is filtration alone, above 1
            tubular secretion. 70 ml/min of plasma [Lo1995] is about 59 (the
            kidney extracts from its blood).
            """,
        ),
        Parameter(
            "E3174EX_f",
            **FITS.parameter("KI__E3174EX_f", unit=U.dimensionless),
            name="renal clearance of unbound E3174 relative to the GFR",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of E3174 is E3174EX_f times the filtration of
            its unbound fraction, fu_e3174 * GFR. 26 ml/min of plasma [Lo1995]
            is about 129.
            """,
        ),
        Parameter(
            "L158EX_f",
            **FITS.parameter("KI__L158EX_f", unit=U.dimensionless),
            name="renal clearance of the other metabolites relative to the GFR",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The other metabolites are not measured in plasma, only the
            fraction of the urine against the feces is known (about 1:1
            [FDA1995S60]); their binding is unknown (fu taken as 1).
            """,
        ),
    ]
)
_m.rules.append(
    AssignmentRule(
        "GFR",
        "f_renal_function * GFR_base",
        unit=U.ml_per_min,
        name="glomerular filtration rate",
    ),
)

_m.reactions = [
    Reaction(
        sid="LOSEX",
        name="losartan excretion (LOSEX)",
        equation="los_ext -> los_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        # [ml/min] * [l/ml] * [-] * [mmole/l]
        formula=(
            "GFR * 0.001 l_per_ml * fu_los * LOSEX_f * los_ext / BP_los",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="E3174EX",
        name="E3174 excretion (E3174EX)",
        equation="e3174_ext -> e3174_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "GFR * 0.001 l_per_ml * fu_e3174 * E3174EX_f * e3174_ext / BP_e3174",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="L158EX",
        name="other metabolites excretion (L158EX)",
        equation="l158_ext -> l158_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "GFR * 0.001 l_per_ml * L158EX_f * l158_ext / BP_l158",
            U.mmole_per_min,
        ),
    ),
]

model_kidney = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.losartan import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
