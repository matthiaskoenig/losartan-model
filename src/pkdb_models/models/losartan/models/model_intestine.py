"""Losartan intestine model."""
import numpy as np
import pandas as pd

from sbmlutils import cytoscape as cyviz
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.losartan.models import annotations
from pkdb_models.models.losartan.models import templates


class U(templates.U):
    """UnitDefinitions"""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "losartan_intestine",
    name="Model for losartan absorption in the small intestine",
    notes="""
    # Model for losartan absorption

    - absorption losartan (los), ~60% fraction absorbed (bioavailability ~33%, due to liver metabolism)
    - enterohepatic circulation of los, e3174 and l158
    - ABCB1 (P-glycoprotein) influences pharmacokinetics, possible role in intestinal uptake
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
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
        f"los_stomach",
        metaId=f"meta_los_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name=f"losartan (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        boundaryCondition=True,
    ),
    Species(
        "los_lumen",
        initialConcentration=0.0,
        name="losartan (lumen)",
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
        name="losartan (plasma)",
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
        name="E3174 (lumen)",
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
        name="L158 (lumen)",
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
        name="L158 (feces)",
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
        f"F_los_abs",
        0.6,
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"fraction absorbed losartan",
        notes="""
        Fraction absorbed, i.e., only a fraction of the losartan in the intestinal lumen
        is absorbed. This parameter determines how much of the losartan is excreted.
        
        `F_los_abs` of dose is absorbed. `(1-F_los_abs)` is excreted in feces.
        """,
    ),
    Parameter(
        "LOSABS_k",
        0.03617722689911161,
        unit=U.per_min,
        name="rate of losartan absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "METEXC_k",
        2.6964160819935065e-05,
        unit=U.per_min,
        name="rate of feces excretion",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
]

_m.rules.append(
    AssignmentRule(
        "absorption",
        value="LOSABS_k * Vgu * los_lumen",
        unit=U.mmole_per_min,
        name="absorption losartan",
    ),
)

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
        formula=("f_OATP2B1 * F_los_abs * absorption", U.mmole_per_min),
        notes="""
        Import via OATP2B1 is irreversible at apical membrane of enterocytes.
        """,
        #formula=("F_los_abs * absorption", U.mmole_per_min),
    ),
    #fraction excreted (not available for absorption)
    Reaction(
        sid="LOSEFL",
        name=f"efflux losartan (PG)",
        compartment="Vapical",
        equation=f"los -> los_lumen",
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
                1.020035540752619,
                unit=U.dimensionless,
                name="fractional rate losartan efflux (PG)",
                sboTerm=SBO.MAXIMAL_VELOCITY,
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
        Assumption: transport over basolateral membrane is similar to apical memmrane
        in the enterocytes; reusing rate constant for absorption.
        """
    ),

    Reaction(
        sid="LOSEXC",
        name=f"excretion losartan (feces)",
        compartment="Vlumen",
        equation=f"los_lumen -> los_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"(1 dimensionless - F_los_abs) * absorption",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="E3174EXC",
        name=f"excretion e3174 (feces)",
        compartment="Vlumen",
        equation=f"e3174_lumen -> e3174_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"METEXC_k * Vgu * e3174_lumen",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="L158EXC",
        name=f"excretion l158 (feces)",
        compartment="Vlumen",
        equation=f"l158_lumen -> l158_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[],
        formula=(
            f"METEXC_k * Vgu * l158_lumen",
            U.mmole_per_min,
        ),
    ),
]

_m.parameters.extend([
    Parameter(
        f"PODOSE_los",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose losartan [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_dis_los",
        2.0,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"rate dissolution losartan",
        port=True
    ),
    Parameter(
        f"Mr_los",
        422.911,
        U.g_per_mole,
        constant=True,
        name=f"Molecular weight losartan [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
        port=True,
    ),
])

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"dissolution_los",
            name=f"dissolution losartan",
            formula=(
                f"Ka_dis_los/60 min_per_hr * PODOSE_los/Mr_los",
                U.mmole_per_min,
            ),
            equation=f"los_stomach -> los_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """
        ),
    ]
)
_m.rate_rules.append(
    RateRule(f"PODOSE_los", f"-dissolution_los * Mr_los", U.mg_per_min),
)

model_intestine = _m

def intestine_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition."""

    positions = [
        # sid, x, y

        ["los_ext", 0, 0],

        ["LOSEX", 0, 0.5 * dy],

        ["los", 0, dy],

        ["LOSABS", -0.5*dx, 1.5*dy],
        ["LOSEFL", 0.5*dx, 1.5*dy],

        ["los_lumen", 0, 2 * dy],
        ["e3174_lumen", 1.5 * dx, 2 * dy],
        ["l158_lumen", 2.5 * dx, 2 * dy],

        ["dissolution_los", -0.5 * dx, 2.5 * dy],
        ["LOSEXC", 0.5 * dx, 2.5 * dy],
        ["E3174EXC", 1.5 * dx, 2.5 * dy],
        ["L158EXC", 2.5 * dx, 2.5 * dy],

        ["los_stomach", -0.5 * dx, 3 *dy],
        ["los_feces", 0.5 * dx, 3 * dy],
        ["e3174_feces", 1.5 * dx, 3 * dy],
        ["l158_feces", 2.5 * dx, 3 *dy],

    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df

def intestine_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'plasma','enterocytes', intestinal 'lumen' and 'feces'."""

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }
    annotations = [

        #plasma
        cyviz.AnnotationShape(
            x_pos= -dx,
            y_pos= -0.25 * dy,
            width= 2 * dx,
            height=0.75 * dy,
            fill_color="#FF0000",  # pale pink
            **kwargs
        ),

        # enterocytes
        cyviz.AnnotationShape(
            x_pos= -dx,
            y_pos= 0.5 * dy,
            width= 2 * dx,
            height=1 * dy,
            fill_color="#FFC0CB", # pastel pink
            **kwargs
        ),
        # lumen
        cyviz.AnnotationShape(
            x_pos= -dx,
            y_pos= 1.5 * dy,
            width= 4 * dx,
            height= 1 * dy,
            fill_color="#0000FF", # blue
            **kwargs
        ),
        # stomach
        cyviz.AnnotationShape(
            x_pos= -dx,
            y_pos= 2.5 * dy,
            width= 1 * dx,
            height= 0.75 * dy,
            fill_color="#FFA500", # orange
            **kwargs
        ),
        # feces
        cyviz.AnnotationShape(
            x_pos=0,
            y_pos= 2.5 * dy,
            width= 3 * dx,
            height= 0.75 * dy,
            fill_color="#E0FFD1", # light green
            **kwargs
        ),
    ]
    return annotations

if __name__ == "__main__":
    from pkdb_models.models.losartan import MODEL_BASE_PATH
    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3, sbml_version=2
    )

    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=intestine_layout())
    cyviz.add_annotations(annotations=intestine_annotations())

