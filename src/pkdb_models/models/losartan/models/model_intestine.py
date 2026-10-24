"""Losartan intestine model."""

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    FactoryResult,
    Model,
    Parameter,
    RateRule,
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
    """Units of the intestine model."""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "losartan_intestine",
    name="Model for losartan absorption in the small intestine",
    notes="""
    # Model for losartan absorption

    - dissolution of the tablet and transit into the gut lumen
    - absorption of losartan (`LOSABS`, OATP2B1) into the enterocytes, efflux
      back into the lumen (`LOSEFL`, P-glycoprotein/ABCB1) and transport into
      the portal blood (`LOSEX`); the absorption is nearly complete, the
      bioavailability of 33 % is the hepatic first pass [Lo1995]
    - the fecal excretion of the content of the gut lumen with the transit
      through the colon (`METEXC_k`): unabsorbed losartan and the E3174 and
      other metabolites of the bile

    The species of the gut lumen (`*_lumen`) are concentrations in the lumen
    volume `Vlumen`.
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "ncit/C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine
        (BQB.HAS_PROPERTY, "ncit/C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
        *annotations.model,
    ],
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["blood"],
    ),
    Compartment(
        "Vgu",
        1.2825,  # 0.0171 [l/kg] * 75 kg
        name="intestine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["gu"],
    ),
    Compartment(
        "Vlumen",
        1.2825 * 0.9,  # 0.0171 [l/kg] * 75 kg * 0.9,
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Ventero",
        1.0,
        name="intestinal lining (enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["apical"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vbaso",
        np.nan,
        name="basolateral membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["basolateral"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["stomach"],
    ),
]


_m.species = [
    Species(
        "los_stomach",
        metaId="meta_los_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name="losartan (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        boundaryCondition=True,
    ),
    Species(
        "los_lumen",
        initialConcentration=0.0,
        name="losartan (gut lumen)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        port=True,
    ),
    Species(
        "los",
        initialConcentration=0.0,
        name="losartan (enterocytes)",
        compartment="Ventero",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
    ),
    Species(
        "los_ext",
        initialConcentration=0.0,
        name="losartan (blood)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        port=True,
    ),
    Species(
        "los_feces",
        initialConcentration=0.0,
        name="losartan (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        port=True,
    ),
    Species(
        "e3174_lumen",
        initialConcentration=0.0,
        name="E3174 (gut lumen)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["e3174"],
        port=True,
    ),
    Species(
        "e3174_feces",
        initialConcentration=0.0,
        name="E3174 (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["e3174"],
        port=True,
    ),
    Species(
        "l158_lumen",
        initialConcentration=0.0,
        name="other losartan metabolites (gut lumen)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["l158"],
        port=True,
    ),
    Species(
        "l158_feces",
        initialConcentration=0.0,
        name="other losartan metabolites (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["l158"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        "LOSABS_k",
        **FITS.parameter("GU__LOSABS_k", unit=U.per_min),
        name="rate of losartan absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="Fitted: tmax of losartan about 1 hr.",
    ),
    Parameter(
        "METEXC_k",
        6e-4,
        unit=U.per_min,
        name="rate of feces excretion",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        The transit of the content of the gut lumen into the feces, the
        colonic transit time of 24-36 hr (mean residence 1/METEXC_k = 28 hr).
        Not fitted: the data of the feces are amounts at the end of the
        collection (96-168 hr), which do not identify the rate.
        """,
    ),
]

_m.reactions = [
    Reaction(
        "LOSABS",
        name="absorption losartan",
        equation="los_lumen -> los",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
            Parameter(
                sid="f_OATP2B1",
                name="absorption activity",
                value=1.0,
                unit=U.dimensionless,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
                parameter controls activity (1.0 reference activity,
                < 1.0 reduced activity, > 1.0 increased activity)
                """,
            ),
        ],
        formula=("f_OATP2B1 * LOSABS_k * Vgu * los_lumen", U.mmole_per_min),
        notes="""
        Import via OATP2B1 at the apical membrane of the enterocytes.
        `f_OATP2B1` is an interaction factor, set by the experiments (no
        experiment uses it at the moment).
        """,
    ),
    # efflux back into the lumen
    Reaction(
        sid="LOSEFL",
        name="efflux losartan (PG)",
        compartment="Vapical",
        equation="los -> los_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                sid="f_abcb1",
                name="PG activity",
                value=1.0,
                unit=U.dimensionless,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
               parameter controls activity (1.0 reference activity,
               < 1.0 reduced activity, > 1.0 increased activity)
               """,
            ),
            Parameter(
                "f_LOSEFL_k",
                1.0,
                unit=U.dimensionless,
                name="fractional rate losartan efflux (PG)",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                The efflux relative to the absorption, not fitted: the
                absorption is nearly complete, so the efflux delays the
                absorption (the ABCB1 genotypes change Cmax and the early AUC,
                not the AUC [Shin2020]) and is not identifiable separately
                from `LOSABS_k`.
                """,
            ),
        ],
        rules=[
            AssignmentRule("LOSEFL_k", "f_LOSEFL_k * LOSABS_k", unit=U.per_min),
        ],
        formula=(
            "f_abcb1 * LOSEFL_k * Vgu * los",
            U.mmole_per_min,
        ),
        notes="""
        Efflux is catalyzed by P-glycoprotein.
        This is calibrated relative to the absorption rate.
        """,
    ),
    Reaction(
        "LOSEX",
        name="losartan export plasma",
        equation="los -> los_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vbaso",
        pars=[],
        formula=(
            "LOSABS_k * Vgu * los",
            U.mmole_per_min,
        ),
        notes="""
        Assumption: the transport over the basolateral membrane of the
        enterocytes has the rate constant of the absorption.
        """,
    ),
    Reaction(
        sid="LOSEXC",
        name="excretion losartan (feces)",
        compartment="Vlumen",
        equation="los_lumen -> los_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            "METEXC_k * Vgu * los_lumen",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="E3174EXC",
        name="excretion E3174 (feces)",
        compartment="Vlumen",
        equation="e3174_lumen -> e3174_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            "METEXC_k * Vgu * e3174_lumen",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="L158EXC",
        name="excretion other metabolites (feces)",
        compartment="Vlumen",
        equation="l158_lumen -> l158_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            "METEXC_k * Vgu * l158_lumen",
            U.mmole_per_min,
        ),
    ),
]

_m.parameters.extend(
    [
        Parameter(
            "PODOSE_los",
            0,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="oral dose losartan [mg]",
            port=True,
        ),
        Parameter(
            "Ka_dis_los",
            2.0,
            U.per_hr,
            constant=True,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="rate dissolution losartan",
            port=True,
        ),
        Parameter(
            "Mr_los",
            422.911,
            U.g_per_mole,
            constant=True,
            name="Molecular weight losartan [g/mole]",
            sboTerm=SBO.MOLECULAR_MASS,
            port=True,
        ),
    ]
)

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid="dissolution_los",
            name="dissolution losartan",
            formula=(
                "Ka_dis_los/60 min_per_hr * PODOSE_los/Mr_los",
                U.mmole_per_min,
            ),
            equation="los_stomach -> los_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """,
        ),
    ]
)
_m.rate_rules.append(
    RateRule("PODOSE_los", "-dissolution_los * Mr_los", U.mg_per_min),
)

model_intestine = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.losartan import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
