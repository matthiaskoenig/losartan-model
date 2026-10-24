"""Liver model for losartan.

Losartan is cleared by the liver: 14 % of a dose is oxidized to the active
carboxylic acid E3174 (CYP2C9, minor CYP3A4), most of the rest forms inactive
metabolites (the hydroxylated L-158,795/796, L-158,783 and glucuronides, by
CYP2C9, CYP3A4 and UGTs) [Lo1995, FDA1995S60]. E3174 and the other metabolites
are excreted into the urine and the bile.
"""

import numpy as np
from sbmlutils.factory import (
    Compartment,
    FactoryResult,
    Model,
    Parameter,
    Reaction,
    Species,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.losartan.models import annotations, templates
from pkdb_models.models.losartan.models.fits import FITS


class U(templates.U):
    """Units of the liver model."""


mid = "losartan_liver"
version = 2

_m = Model(
    sid=mid,
    name="Model for hepatic losartan metabolism.",
    notes="""
    Model for the hepatic metabolism of losartan.

    - losartan exchanges fast between the blood of the liver and the
      hepatocytes (`LOSIM`)
    - losartan -> E3174 (CYP2C9, `LOS2E3174`)
    - losartan -> other metabolites (`LOS2L158`, CYP2C9 for the fraction
      `f_l158_cyp2c9`, CYP3A4 and UGTs for the rest)
    - export of E3174 and the other metabolites into the blood and the bile

    The concentrations of the liver model are blood concentrations (the
    species of the whole-body model), the rate constants are intrinsic
    clearances of blood per liver volume. CYP2C9 forms about 55-60 % of the
    hepatic clearance: the loss of CYP2C9 (*3/*3, bucolome) increases the
    AUC of losartan 1.7-2.4 fold and decreases the AUC of E3174 to 0.1-0.2
    [Bae2011, Kobayashi2008].
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "ncit/C12392"),
        (BQB.HAS_PROPERTY, "ncit/C79371"),  # Pharmacokinetics: Metabolism
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
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
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
        "Vapical",
        np.nan,
        name="apical membrane",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["apical"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vbi",
        1.0,
        name="bile",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["bi"],
        port=True,
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
]

info = {
    "los": "losartan",
    "e3174": "E3174",
    "l158": "other losartan metabolites",
}
#: the substances which are excreted into the bile
BILE: tuple[str, ...] = ("e3174", "l158")

for sid, name in info.items():
    _m.species.extend(
        [
            Species(
                f"{sid}_ext",
                name=f"{name} (blood)",
                initialConcentration=0.0,
                compartment="Vext",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,  # this is a concentration
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
                port=True,
            ),
            Species(
                sid,
                name=f"{name} (liver)",
                initialConcentration=0.0,
                compartment="Vli",
                substanceUnit=U.mmole,
                hasOnlySubstanceUnits=False,  # this is a concentration
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            ),
        ]
    )
    if sid in BILE:
        _m.species.extend(
            [
                Species(
                    f"{sid}_bi",
                    initialConcentration=0.0,
                    name=f"{name} (bile)",
                    compartment="Vbi",
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=True,
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                    annotations=annotations.species[sid],
                    notes=f"Bile {name} amount.",
                ),
                Species(
                    f"{sid}_lumen",
                    initialConcentration=0.0,
                    name=f"{name} (gut lumen)",
                    compartment="Vlumen",
                    substanceUnit=U.mmole,
                    hasOnlySubstanceUnits=False,
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                    annotations=annotations.species[sid],
                    port=True,
                ),
            ]
        )

_m.parameters = [
    Parameter(
        "f_cyp2c9",
        1.0,
        unit=U.dimensionless,
        name="activity of CYP2C9",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Relative activity of CYP2C9: 1.0 normal activity (wildtype), <1.0
        reduced activity (genotypes, inhibition), >1.0 increased activity
        (induction). Set by the experiments (genotype, interactions), not
        fitted.
        """,
    ),
    Parameter(
        "f_l158_cyp2c9",
        0.5,
        unit=U.dimensionless,
        name="fraction of the formation of the other metabolites by CYP2C9",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        The oxidative metabolites of losartan (L-158,795/796) are formed by
        CYP2C9 and CYP3A4, the glucuronides by UGTs. 0.5 gives together with
        the formation of E3174 a CYP2C9 share of about 0.55 of the hepatic
        clearance: losartan AUC x2.4 and E3174 AUC x0.12 in CYP2C9 *3/*3
        [Bae2011], losartan AUC x1.67 and E3174 AUC x0.165 with bucolome
        [Kobayashi2008].
        """,
    ),
]

_m.reactions = [
    Reaction(
        sid="LOSIM",
        name="losartan import (LOSIM)",
        equation="los_ext <-> los",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "LOSIM_k",
                50,
                U.per_min,
                name="rate losartan import",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            )
        ],
        formula="LOSIM_k * Vli * (los_ext - los)",
        notes="""
        Exchange of losartan between the blood and the hepatocytes, partly by
        OATP transporters (OATP2B1). Assumption: fast exchange, i.e., the
        concentration of the hepatocytes is the one of the blood and the
        hepatic clearance is limited by the metabolism and the blood flow
        (extraction 0.68 [Lo1995]).
        """,
    ),
    Reaction(
        sid="LOS2E3174",
        name="E3174 formation (CYP2C9)",
        equation="los -> e3174",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "LOS2E3174_k",
                **FITS.parameter("LI__LOS2E3174_k", unit=U.per_min),
                name="rate losartan → E3174",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                First order: the Km of CYP2C9 for losartan (about 1.1 µM
                unbound [Maekawa2009, Wang2014]) is far above the unbound
                concentrations (fu 1.3 %), the exposure is dose proportional
                from 25 to 200 mg [Ohtawa1993]. 14 % of a dose is converted
                [Lo1995].
                """,
            ),
        ],
        formula="f_cyp2c9 * LOS2E3174_k * Vli * los",
        notes="""
        Oxidation of losartan via the aldehyde E3179 to the carboxylic acid
        E3174, mainly by CYP2C9 (minor CYP3A4).
        """,
    ),
    Reaction(
        sid="LOS2L158",
        name="other metabolites formation (CYP2C9, CYP3A4, UGT)",
        equation="los -> l158",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "LOS2L158_k",
                **FITS.parameter("LI__LOS2L158_k", unit=U.per_min),
                name="rate losartan → other metabolites",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                The hepatic clearance which does not form E3174, about 85 % of
                the non-renal clearance [Lo1995, FDA1995S60].
                """,
            ),
        ],
        formula=(
            "(f_l158_cyp2c9 * f_cyp2c9 + (1 dimensionless - f_l158_cyp2c9))"
            " * LOS2L158_k * Vli * los"
        ),
        notes="""
        The inactive metabolites of losartan: the hydroxylated L-158,795 and
        L-158,796 and L-158,783 (CYP2C9, CYP3A4) and the glucuronides (UGT).
        They are lumped into one species.
        """,
    ),
    Reaction(
        sid="E3174EX",
        name="E3174 export (E3174EX)",
        equation="e3174 <-> e3174_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "E3174EX_k",
                **FITS.parameter("LI__E3174EX_k", unit=U.per_min),
                name="rate E3174 export",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
                notes="""
                Export of E3174 from the hepatocytes into the blood, fitted. It
                is slow and, with the bile export, sets the retention of the
                E3174 formed in the hepatocytes, i.e., the delayed appearance
                of E3174 in plasma (tmax 3-4 hr after intravenous and oral
                losartan [Lo1995]).
                """,
            ),
            Parameter(
                "E3174IM_k",
                **FITS.parameter("LI__E3174IM_k", unit=U.per_min),
                name="rate E3174 import",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
                notes="""
                Uptake of E3174 from the blood into the hepatocytes (OATP),
                concentrative. With the bile export it sets the non-renal
                clearance of E3174 (21 ml/min of plasma, 45 % of its clearance
                [Lo1995]).
                """,
            ),
        ],
        formula="Vli * (E3174EX_k * e3174 - E3174IM_k * e3174_ext)",
    ),
    Reaction(
        sid="L158EX",
        name="other metabolites export (L158EX)",
        equation="l158 <-> l158_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "L158EX_k",
                1.0,
                U.per_min,
                name="rate export of the other metabolites",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            ),
        ],
        formula="L158EX_k * Vli * (l158 - l158_ext)",
        notes="Assumption: fast transport.",
    ),
    Reaction(
        sid="E3174BIEX",
        name="E3174 bile export",
        equation="e3174 -> e3174_bi",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
            Parameter(
                "E3174BIEX_k",
                **FITS.parameter("LI__E3174BIEX_k", unit=U.per_min),
                name="rate E3174 bile export",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                The non-renal clearance of E3174 is 45 % of its plasma
                clearance (47 ml/min, renal 26 ml/min [Lo1995]), about 10 % of
                a dose of losartan is found as E3174 in the feces
                [FDA1995S60].
                """,
            ),
        ],
        formula=("E3174BIEX_k * Vli * e3174", U.mmole_per_min),
    ),
    Reaction(
        sid="L158BIEX",
        name="other metabolites bile export",
        equation="l158 -> l158_bi",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
            Parameter(
                "L158BIEX_k",
                **FITS.parameter("LI__L158BIEX_k", unit=U.per_min),
                name="rate bile export of the other metabolites",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                The metabolites of losartan which are not E3174 are recovered
                in about equal parts in the urine and the feces (17.9 and
                17.7 % of an intravenous dose [FDA1995S60]).
                """,
            ),
        ],
        formula=("L158BIEX_k * Vli * l158", U.mmole_per_min),
    ),
]

# the bile flows into the gut lumen
for sid in BILE:
    _m.reactions.append(
        Reaction(
            f"{sid.upper()}EHC",
            name=f"{info[sid]} bile into the gut lumen",
            equation=f"{sid}_bi -> {sid}_lumen",
            sboTerm=SBO.TRANSPORT_REACTION,
            compartment="Vlumen",
            formula=(f"{sid.upper()}BIEX", U.mmole_per_min),
        ),
    )


model_liver = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.losartan import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
