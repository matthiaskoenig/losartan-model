"""Liver model for losartan.
"""

import numpy as np
import pandas as pd


from sbmlutils import cytoscape as cyviz
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *
from pkdb_models.models.templates import terms_of_use

from pkdb_models.models.losartan.models import annotations
from pkdb_models.models.losartan.models import templates


class U(templates.U):
    """UnitDefinitions"""

    pass


mid = "losartan_liver"
version = 1

_m = Model(
    sid=mid,
    name="Model for hepatic losartan metabolism.",
    notes=f"""
    Model for losartan metabolism.
    
    - metabolism of losartan -> exp3174 (liver CYP2C9)
    """ + terms_of_use,
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
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
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
    "l158": "L158",
}

for sid, name in info.items():
    _m.species.extend([
        Species(
            f"{sid}_ext",
            name=f"{name} (plasma)",
            initialConcentration=0.0,
            compartment="Vext",
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,  # this is a concentration
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
            port=True
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
        Species(
            f"{sid}_bi",
            initialConcentration=0.0,
            name=f"{name} (bile)",
            compartment="Vbi",
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=True,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
            notes=f"""
            Bile {name} amount.
            """,
        ),
        Species(
            f"{sid}_lumen",
            initialConcentration=0.0,
            name=f"{name} (lumen)",
            compartment="Vlumen",
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
            port=True,
        ),
    ])

_m.parameters = [
    Parameter(
        "f_cyp2c9",
        1.0,
        unit=U.dimensionless,
        name="activity of CYP2C9",
        sboTerm=SBO.MAXIMAL_VELOCITY,
        notes="""Relative activity of the enzyme: 1.0 - normal activity (wildtype); 
    <1.0 reduced activity; >1.0 increased activity."""
    ),
    Parameter(
        "MBIEX_k",
        0.066315021777011,
        U.per_min,
        name="rate for export in bile",
        sboTerm=SBO.KINETIC_CONSTANT,
    )
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
                1,
                U.per_min,
                name="rate losartan import",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            )
        ],
        formula=(
            "LOSIM_k * Vli * (los_ext - los)"
        ),
        notes="""
        Organic Anion Transporting Polypeptide 2B1 (OATP2B1): Recent findings suggest that losartan is also a substrate for OATP2B1, a transporter that plays a crucial role in the uptake of various drugs into hepatocytes. This transporter facilitates the entry of losartan into liver cells, where it can be metabolized by cytochrome P450 enzymes. [perplexity]
        Assumption: fast transport.
        """

    ),
    Reaction(
        sid="LOS2E3174",
        name="E3174 formation (CYP2C9/3A4)",
        equation="los -> e3174",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "LOS2E3174_Vmax",
                0.0007251463173167296,
                U.mmole_per_min_l,
                name="Vmax losartan conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "LOS2E3174_Km_los",
                0.0011,  # 1.1 µM
                U.mM,
                name="Km losartan conversion",
                sboTerm=SBO.MICHAELIS_CONSTANT,
                notes="""FIXME: check in vitro data
                
                1.12 +- 0.13 µM [Maekawa2009]
                1.134 µM [Thu2017]
                1.05 +- 0.017 µM [Wang2014]
                """
            )
        ],
        formula=(
            "f_cyp2c9 * LOS2E3174_Vmax * Vli * los/(los + LOS2E3174_Km_los)"
        ),
        notes="""
        Likely catalyzed CYP2C9 and CYP3A4.
        Currently only modified via CYP2C9 factor.
        """
    ),
    Reaction(
        sid="E3174L158",
        name="L158 formation (UGT)",
        equation="e3174 -> l158",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "E3174L158_k",
                0.0011340090993695589,
                U.per_min,
                name="rate losartan conversion",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            ),
        ],
        formula=(
            "E3174L158_k * Vli * e3174"
        ),
        notes="""
        Formation of losartan glucuronides, likely via UGT enzymes.
        Assumption: e3174 metabolite is converted; not much information available on pathway.
        """
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
                0.011200118227296746,
                U.per_min,
                name="rate E3174 export",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            ),
        ],
        formula=(
            "E3174EX_k * Vli * (e3174 - e3174_ext)"
        ),
        notes="""assumption fast transport"""
    ),
    Reaction(
        sid="L158EX",
        name="L158 export (L158EX)",
        equation="l158 <-> l158_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "L158EX_k",
                1.0,
                U.per_min,
                name="rate L158 export",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            ),
        ],
        formula=(
            "L158EX_k * Vli * (l158 - l158_ext)"
        ),
        notes="""assumption fast transport"""
    ),
]

# Enterohepatic circulation (EHC), bile excretion
for sid, name in info.items():
    _m.reactions.extend([
        Reaction(
            sid=f"{sid.upper()}BIEX",
            name=f"{name.title()} bile export",
            equation=f"{sid} -> {sid}_bi",
            sboTerm=SBO.TRANSPORT_REACTION,
            compartment="Vapical",
            formula=(
                f"MBIEX_k * Vli * {sid}",
                U.mmole_per_min,
            ),
            notes="""
            Assumption: identical bile transport of los, e3174 and l158.
            """
        ),
        Reaction(
            f"{sid.upper()}EHC",
            name=f"{name.title()} enterohepatic circulation",
            equation=f"{sid}_bi -> {sid}_lumen",
            sboTerm=SBO.TRANSPORT_REACTION,
            compartment="Vlumen",
            formula=(f"{sid.upper()}BIEX", U.mmole_per_min),
        ),
    ])


model_liver = _m

def liver_layout(dx=200, dy=200) -> pd.DataFrame:
    """Layout definition."""

    positions = [
        # sid, x, y

        ["los_ext", -dx, 0],
        ["e3174_ext", 0, 0],
        ["l158_ext", dx, 0],

        ["LOSIM", -dx, 0.5*dy],
        ["E3174EX", 0, 0.5*dy],
        ["L158EX", dx, 0.5*dy],

        ["los", -dx, dy],
        ["e3174", 0, dy],
        ["l158", dx, dy],

        ["LOS2E3174", -0.5*dx, 1.25*dy],
        ["E3174L158", 0.5*dx, 1.25*dy],

        ["LOSBIEX", -dx, 1.5 * dy],
        ["E3174BIEX", 0, 1.5 * dy],
        ["L158BIEX", dx, 1.5 * dy],

        ["los_bi", -dx, 2*dy],
        ["e3174_bi", 0, 2*dy],
        ["l158_bi", dx, 2*dy],

        ["LOSEHC", -dx, 2.5 * dy],
        ["E3174EHC", 0, 2.5 * dy],
        ["L158EHC", dx, 2.5 * dy],

        ["los_lumen", -dx, 3 * dy],
        ["e3174_lumen", 0, 3 * dy],
        ["l158_lumen", dx, 3 * dy],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df

def liver_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'plasma','liver', 'bile' and liver 'lumen'."""

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

        # liver
        cyviz.AnnotationShape(
            x_pos=-1.5 * dx,
            y_pos=0.5 * dy,
            width=3 * dx,
            height=1 * dy,
            fill_color="#FFFFFF",
            **kwargs
        ),
        # bile
        cyviz.AnnotationShape(
            x_pos=-1.5 * dx,
            y_pos=1.5 * dy,
            width=3 * dx,
            height=1 * dy,
            fill_color="#00FF00",
            **kwargs
        ),
        # lumen
        cyviz.AnnotationShape(
            x_pos=-1.5 * dx,
            y_pos=2.5 * dy,
            width=3 * dx,
            height=0.75 * dy,
            fill_color="#0000FF",
            **kwargs
        ),
    ]
    return annotations

if __name__ == "__main__":
    from pkdb_models.models.losartan import MODEL_BASE_PATH
    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3, sbml_version=2,
    )
    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=liver_layout())
    cyviz.add_annotations(annotations=liver_annotations())
