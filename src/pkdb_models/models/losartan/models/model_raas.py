"""Model of RAAS blood pressure regulation."""
from dataclasses import dataclass

import numpy as np
import pandas as pd


from sbmlutils import cytoscape as cyviz
from sbmlutils.converters import odefac
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.losartan.models import annotations
from pkdb_models.models.losartan.models import templates


class U(templates.U):
    """UnitDefinitions"""
    mmHg = UnitDefinition("mmHg", "133.32239 N/m^2")


mid = "losartan_raas"
version = 2

_m = Model(
    sid=mid,
    name="Model for RAAS system of blood pressure regulation.",
    notes=f"""
    Model for RAAS system of blood pressure regulation.
    
    **version** {version}
    
    ## Changelog
    **version 2**
    
    - better handling of baselines and scaling
    
    **version 1**
    
    - initial model
        
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
)

_m.compartments = [
    Compartment(
        "Vplasma",
        value=5.0,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
]

@dataclass
class Substance:
    sid: str
    name: str
    init: float
    unit: str
    boundary: bool = False
    port: bool = False


# ---------------------------------------------------------------------------------------------------------------------
# Pharmacodynamics
# ---------------------------------------------------------------------------------------------------------------------
ald_ref = 250E-9  # [mM]
ren_ref = 1E-9  # [mM]
anggen_ref = 10E-9  # [mM]
ang1_ref = 10E-9  # [mM]
ang2_ref = 8E-9  # [mM]

substances: list[Substance] = [


    # Mr: 360.444 [g/mole]
    # 100 pg/ml {ramipril/Bussien1985} = 100/360.444*1000 [pmol/l] = 277.4 [pmol/l]
    # 80-90 pg/ml {losartan/Christen1991a} = 221.9 - 249.7 [pmol/l]
    # 80-100 pg/ml {losartan/Ohtawa1993} = 221.9 - 277.4 [pmol/l]
    # 600-800 pg/ml {losartan/Doig1993} FIXME: check unit = 1942.0 - 2219.5 [pmol/l]  #
    Substance("ald", "aldosterone", init=ald_ref, unit="mmole"),

    # Mr: 45057 [g/mole]
    # ~ 50-65 pg/ml {losartan/Azizi1999} = 50-65/45.057 [pmol/l] = 1.11 - 1.44 [pmol/l]
    # ~ 45-75 pg/ml {losartan/Doig1993}
    # 1 ng/mL/h PRA ≈ 7.6 pg/ml (5.5 - 9.7)?
    # ~2 ng/ml/hr {losartan/Goldberg1995} = 15.2 pg/ml = 0.337 [pmol/l]
    # ~2 ng/ml/hr {losartan/Tsuruoka2005} = 0.337 [pmol/l]
    # 1-2 ng/ml/hr {losartan/Ohtawa1993} = 0.169 - 0.337 [pmol/l]
    # 1-2 ng/ml/hr renin activity {losartan/Christen1991a} = 0.169 - 0.337 [pmol/l]
    # 5-6 ng/ml/hr renin activity {losartan/Doig1993} FIXME: check unit = 0.843 - 1.01 [pmol/l]
    Substance("ren", "renin", init=ren_ref, unit="mmole"),

    # Mr: 52670 [g/mole]
    Substance("anggen", "angiotensinogen", init=anggen_ref, unit="mmole", boundary=True),

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

    # Mr: 149715 [g/mole]
    # 10 pmol/l: assumption, same order as ang1 and ang and renin
    # Substance("ace", "ACE", init=1E-9, unit="mmole"),
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
            # FIXME: better handling of initial concentrations
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
            AssignmentRule(
                f"{s.sid}_change", f"{s.sid}-{s.sid}_ref", unit=U.mM
            )
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

_m.parameters.extend([
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
        name="half-maximum effect concentration exp3174",  # Mr = 436.894 g/mol
        # value=73.24E-6,  # ~ 32 [1-99] ng/ml {Munafo1992} = 32/436.894*1000 [nmol/l] = 73.24 [2.29 - 226.6] [nmol/l]
        value = 0.00029107274011198355,
        unit=U.mM,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
])
_m.rules.extend([
    AssignmentRule(
        "fe_e3174",
        "e3174/(E50_e3174 + e3174)",
        unit=U.dimensionless,
        name="effect via exp3174 in [0, 1]"
    ),
])


_m.reactions = [
    # renin turnover
    Reaction(
        sid="RENSEC",
        name="renin secretion (RENSEC)",
        equation="-> ren [e3174]",  # FIXME: ang2 should affect the renin secretion !?
        compartment="Vplasma",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RENSEC_k",
                0.1,
                U.mmole_per_min,
                name="rate renin secretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
            Parameter(
                "RENSEC_fa",
                5,
                U.dimensionless,
                name="activation renin secretion via exp3174",
                sboTerm=SBO.KINETIC_CONSTANT,
            )
        ],
        formula=("RENSEC_k * (1 dimensionless + RENSEC_fa * fe_e3174)", U.mmole_per_min)
    ),
    Reaction(
        sid="RENDEG",
        name="renin degradation (RENDEG)",
        equation="ren ->",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            # in steady state: kd = ks/ren0
            Parameter(
                "RENDEG_k",
                f"RENSEC_k/ren_ref",
                U.l_per_min,
                name="rate renin degradation",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            )
        ],
        formula=("RENDEG_k * ren", U.mmole_per_min)
    ),

    # aldosterone turnover
    Reaction(
        sid="ALDSEC",
        name="aldosterone secretion (ALDSEC)",
        equation="-> ald [e3174]",
        compartment="Vplasma",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "ALDSEC_k",
                1.0105027952685962e-06,
                U.mmole_per_min,
                name="rate aldosterone secretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            )
        ],
        # formula=("ALDSEC_k * (1 dimensionless - fe_e3174) * piecewise(ang2/ang2_ref, e3174 > 0.01 dimensionless * E50_e3174, 1 dimensionless)", U.mmole_per_min),
        formula=("ALDSEC_k * (1 dimensionless - fe_e3174)", U.mmole_per_min),
        notes="""
        - exp3174 inhibits aldosterone secretion; probably via the inhibition of AT1 receptors, resulting in 
          reduced aldosterone secretion. This is modelled via the fi dependency.
        - angiotensin II activates aldosterone secretion by binding to AT1 receptors
        """
    ),
    Reaction(
        sid="ALDDEG",
        name="aldosterone degradation (ALDDEG)",
        equation="ald ->",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        # in steady state: kd = ks/ren0
        pars=[
            Parameter(
                "ALDDEG_k",
                f"ALDSEC_k/ald_ref",
                U.l_per_min,
                name="rate aldosterone degradation",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            )
        ],
        formula=("ALDDEG_k * ald", U.mmole_per_min)
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
                "ANGGEN2ANG1_k",  # k1
                0.10030014354408065,
                U.l_per_min,
                name="rate angen to ang1 conversion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("ANGGEN2ANG1_k * anggen * ren/ren_ref", U.mmole_per_min)
    ),
    Reaction(
        sid="ANG1ANG2",
        name="angiotensin I to angiotensin II (ACE)",
        equation="ang1 -> ang2",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "ANG1ANG2_k",  # k2
                "ANGGEN2ANG1_k * anggen_ref/ang1_ref",
                U.l_per_min,
                name="rate ang1 to ang2 conversion",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            ),
        ],
        formula=("ANG1ANG2_k * ang1", U.mmole_per_min)
    ),
    Reaction(
        sid="ANG2DEG",
        name="angiotensin II degradation (ANG2DEG)",
        equation="ang2 ->",
        compartment="Vplasma",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "ANG2DEG_k",  # k3
                "ANGGEN2ANG1_k * anggen_ref/ang2_ref",
                U.l_per_min,
                name="rate aldosterone degradation",
                sboTerm=SBO.KINETIC_CONSTANT,
                constant=False,
            )
        ],
        formula=("ANG2DEG_k * ang2", U.mmole_per_min)
    ),
]

# blood pressure model
_m.parameters.extend([
    Parameter(
        "SBP_ref",
        120,
        U.mmHg,
        name="reference systolic blood pressure [mmHg]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        "DBP_ref",
        80,
        U.mmHg,
        name="reference diastolic blood pressure [mmHg]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),

    Parameter(
        "SBP",
        np.nan,
        U.mmHg,
        name="systolic blood pressure [mmHg]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        constant=False,
    ),
    Parameter(
        "DBP",
        np.nan,
        U.mmHg,
        name="diastolic blood pressure [mmHg]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        constant=False,
    ),
    Parameter(
        "MAP",
        np.nan,
        U.mmHg,
        constant=False,
        name="mean arterial pressure [mmHg]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        "BP_ald_fe",
        0.3119992626947359,
        U.dimensionless,
        name="effect of aldosterone on blood pressure [-]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
])
_m.rules.extend([
    AssignmentRule("SBP", "SBP_ref + BP_ald_fe * SBP_ref * (ald-ald_ref)/ald_ref", unit=U.mmHg),
    AssignmentRule("DBP", "DBP_ref + BP_ald_fe * DBP_ref * (ald-ald_ref)/ald_ref", unit=U.mmHg),
    AssignmentRule("MAP", "DBP + (SBP - DBP)/3 dimensionless", unit=U.mmHg),
])

# blood pressure changes
for s.sid in ["SBP", "DBP"]:
    # absolute change
    _m.parameters.append(
        Parameter(
            f"{s.sid}_change",
            name=f"{s.name} change",
            value=np.nan,
            unit=U.mmHg,
            constant=False,
            notes=f"Absolute change to baseline {s.name}",
        )
    )
    _m.rules.append(
        AssignmentRule(
            f"{s.sid}_change", f"{s.sid}-{s.sid}_ref", unit=U.mmHg
        )
    )
    # relative change
    _m.parameters.append(
        Parameter(
            f"{s.sid}_ratio",
            name=f"{s.name} ratio",
            value=np.nan,
            unit=U.dimensionless,
            constant=False,
            notes=f"Ratio relative to baseline {s.name}",
        )
    )
    _m.rules.append(
        AssignmentRule(
            f"{s.sid}_ratio", f"{s.sid}/{s.sid}_ref", unit=U.dimensionless
        )
    )



model_raas = _m

# def raas_layout(dx=200, dy=200) -> pd.DataFrame:
#     """Layout definition."""
#
#     positions = [
#         # sid, x, y
#
#         ["e3174_ext", -dx, 0],
#         ["los_ext", 0, 0],
#         ["l158_ext", dx, 0],
#
#         ["E3174EX", -dx, 0.5*dy],
#         ["LOSEX", 0, 0.5*dy],
#         ["L158EX", dx, 0.5*dy],
#
#         ["e3174_urine", -dx, dy],
#         ["los_urine", 0, dy],
#         ["l158_urine", dx, dy],
#     ]
#
#     df = pd.DataFrame(positions, columns=["id", "x", "y"])
#     df.set_index("id", inplace=True)
#     return df
#
# def raas_annotations(dx=200, dy=200) -> list:
#     """Bounding boxes for 'plasma' and 'urine'."""
#
#     kwargs = {
#         "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
#         "opacity": 20,
#         "border_color": "#000000",
#         "border_thickness": 2,
#     }
#     annotations = [
#
#         #plasma
#         cyviz.AnnotationShape(
#             x_pos=-1.5 * dx,
#             y_pos=-0.25 * dy,
#             width=3 * dx,
#             height=0.75 * dy,
#             fill_color="#FF0000",  # pale pink
#             **kwargs
#         ),
#
#         # urine
#         cyviz.AnnotationShape(
#             x_pos=-1.5 * dx,
#             y_pos=0.5 * dy,
#             width=3 * dx,
#             height=0.75 * dy,
#             fill_color="#FFFFE0", # pale yellow
#             **kwargs
#         ),
#     ]
#     return annotations



if __name__ == "__main__":
    from pkdb_models.models.losartan import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_raas,
        filepath=MODEL_BASE_PATH / f"{model_raas.sid}.xml",
        sbml_level=3, sbml_version=2,
        validation_options=ValidationOptions(units_consistency=True)
    )
    # create differential equations
    md_path = MODEL_BASE_PATH / f"{model_raas.sid}.md"
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=md_path)

    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    # cyviz.apply_layout(layout=kidney_layout())
    # cyviz.add_annotations(annotations=kidney_annotations())
