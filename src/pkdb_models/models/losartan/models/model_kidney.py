"""Kidney model for losartan."""
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

    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")


mid = "losartan_kidney"
version = 1

_m = Model(
    sid=mid,
    name="Model for renal excretion of losartan, E3174, and L158.",
    notes=f"""
    Model for renal excretion of losartan, E3174, and L158.
    
    **version** {version}
    
    ## Changelog
    
    **version 1**
    
    - initial model
        
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
    Compartment(
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True
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
        name="losartan (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["los"],
        port=True
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
        port=True
    ),
    Species(
        "e3174_ext",
        name="e3174 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["e3174"],
        port=True
    ),
    Species(
        "e3174_urine",
        name="e3174 (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,  # this is in amounts
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["e3174"],
        port=True
    ),
    Species(
        "l158_ext",
        name="l158 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["l158"],
        port=True
    ),
    Species(
        "l158_urine",
        name="l158 (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,  # this is in amounts
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["l158"],
        port=True
    ),
]

_m.parameters.extend([
    Parameter(
        "f_renal_function",
        name="renal function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""scaling factor for renal function. 1.0: normal renal function; 
        <1.0: impaired renal function; >1.0 increased renal function
        """
    ),
])

_m.reactions = [
    Reaction(
        sid="LOSEX",
        name="losartan excretion (LOSEX)",
        equation="los_ext -> los_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "LOSEX_k",
                0.0773896174258953,
                U.per_min,
                name="rate urinary excretion losartan",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * LOSEX_k * Vki * los_ext"  # [mmole/min]
        )
    ),
    Reaction(
        sid="E3174EX",
        name="e3174 excretion (E3174EX)",
        equation="e3174_ext -> e3174_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "E3174EX_k",
                0.028927173263331843,
                U.per_min,
                name="rate urinary excretion E3174",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * E3174EX_k * Vki * e3174_ext"  # [mmole/min]
        )
    ),
    Reaction(
        sid="L158EX",
        name="l158 excretion (L158EX)",
        equation="l158_ext -> l158_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "L158EX_k",
                0.28890834506773044,
                U.per_min,
                name="rate urinary excretion L158",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * L158EX_k * Vki * l158_ext"  # [mmole/min]
        )
    ),

]

model_kidney = _m

def kidney_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition."""

    positions = [
        # sid, x, y

        ["e3174_ext", -dx, 0],
        ["los_ext", 0, 0],
        ["l158_ext", dx, 0],

        ["E3174EX", -dx, 0.5*dy],
        ["LOSEX", 0, 0.5*dy],
        ["L158EX", dx, 0.5*dy],

        ["e3174_urine", -dx, dy],
        ["los_urine", 0, dy],
        ["l158_urine", dx, dy],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df

def kidney_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'plasma' and 'urine'."""

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }
    annotations = [

        #plasma
        cyviz.AnnotationShape(
            x_pos=-1.5 * dx,
            y_pos=-0.25 * dy,
            width=3 * dx,
            height=0.75 * dy,
            fill_color="#FF0000",  # pale pink
            **kwargs
        ),

        # urine
        cyviz.AnnotationShape(
            x_pos=-1.5 * dx,
            y_pos=0.5 * dy,
            width=3 * dx,
            height=0.75 * dy,
            fill_color="#FFFFE0", # pale yellow
            **kwargs
        ),
    ]
    return annotations

if __name__ == "__main__":
    from pkdb_models.models.losartan import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=kidney_layout())
    cyviz.add_annotations(annotations=kidney_annotations())
