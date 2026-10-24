"""Losartan body model."""

from typing import NotRequired, TypedDict

from sbmlutils.factory import (
    PORT_SUFFIX,
    AssignmentRule,
    Compartment,
    ExternalModelDefinition,
    Model,
    Parameter,
    RateRule,
    Reaction,
    ReplacedElement,
    Species,
    Submodel,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.losartan.models import annotations, templates
from pkdb_models.models.losartan.models.fits import FITS


# -----------------------------------------------------------------------------
# Whole Body Metabolism
# -----------------------------------------------------------------------------
class U(templates.U):
    """Units of the whole-body model."""

    cm = UnitDefinition("cm")
    hr = UnitDefinition("hr")
    kg = UnitDefinition("kg")

    l_per_hr = UnitDefinition("l_per_hr", "l/hr")
    l_per_kg = UnitDefinition("l_per_kg", "l/kg")
    l_per_min_kg = UnitDefinition("l_per_min_kg", "l/min/kg")
    l_per_ml = UnitDefinition("l_per_ml", "l/ml")

    m2_per_kg = UnitDefinition("m2_per_kg", "meter^2/kg")

    mg_per_day = UnitDefinition("mg_per_day", "mg/day")
    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    mg_per_hr = UnitDefinition("mg_per_hr", "mg/hr")
    mg_per_l = UnitDefinition("mg_per_l", "mg/l")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")

    min_per_day = UnitDefinition("min_per_day", "min/day")
    min_per_hr = UnitDefinition("min_per_hr", "min/hr")

    ml = UnitDefinition("ml")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    ml_per_s = UnitDefinition("ml_per_s", "ml/s")
    ml_per_s_kg = UnitDefinition("ml_per_s_kg", "ml/s/kg")

    mmole_per_hr = UnitDefinition("mmole_per_hr", "mmole/hr")
    mmole_per_l = UnitDefinition("mmole_per_l", "mmole/l")
    mmole_per_min_kg = UnitDefinition("mmole_per_min_kg", "mmole/min/kg")
    mmole_per_hr_ml = UnitDefinition("mmole_per_hr_ml", "mmole/hr/ml")

    mul_per_g = UnitDefinition("mulitre_per_g", "microlitre/g")
    mul_per_min_mg = UnitDefinition("mul_per_min_mg", "microliter/min/mg")

    per_hr = UnitDefinition("per_hr", "1/hr")
    per_min_kg = UnitDefinition("per_min_kg", "1/min/kg")

    s_per_min = UnitDefinition("s_per_min", "s/min")
    mmHg = UnitDefinition("mmHg", "133.32239 N/m^2")


_m = Model(
    "losartan_body",
    name="Losartan body model",
    notes="""
    # Losartan body model

    Whole-body model of losartan, its active metabolite E3174 and its other
    (inactive) metabolites with the blood circulation, the tissue
    distribution, the submodels of the intestine, the liver and the kidneys,
    and the renin-angiotensin-aldosterone system (RAAS) on which E3174 acts.

    The circulating species are whole blood concentrations: the blood flows
    transport the blood concentrations, the plasma concentration is the blood
    concentration over the blood to plasma ratio `BP_<sid>` (losartan 0.53,
    E3174 0.58 [Lo1995], the other metabolites 1). The tissues exchange with
    the blood by the flow `ftissue_<sid>` towards the tissue/blood partition
    coefficient `Kp_<sid>`.

    Literature values of the model [Lo1995, FDA1995S60]: plasma clearance of
    losartan 610 ml/min (renal 70 ml/min, 12 %), Vss 34 l, terminal half-life
    2.1 hr, bioavailability 33 % by hepatic first pass with nearly complete
    absorption, 14 % of a dose is converted to E3174; E3174 plasma clearance
    47 ml/min (renal 26 ml/min, 55 %), Vss 10 l, half-life 6-9 hr. After 14C
    losartan 43 % (IV) and 35 % (PO) of the radioactivity is recovered in the
    urine (losartan 11.4 % and 4.5 %, E3174 6.4 %), 50 % and 58 % in the feces.
    """
    + templates.terms_of_use,
    units=U,
    model_units=templates.model_units,
    creators=templates.creators,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "BTO:0001489"),  # whole body
        (BQB.OCCURS_IN, "FMA:256135"),  # body
        (BQB.HAS_PROPERTY, "ncit/C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "ncit/C79371"),  # Pharmacokinetics: Metabolism
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
        (BQB.HAS_PROPERTY, "ncit/C79370"),  # Pharmacokinetics: Distribution
        *annotations.model,
    ],
)

_m.ports = []
_m.deletions = []
_m.replaced_elements = []


class SubstanceBody(TypedDict):
    """Definition of a substance transported in the whole-body model."""

    name: str
    unit: UnitDefinition
    cinit: float  # initial concentration [mmole/l]
    Mr: float  # molecular weight [g/mole]
    BP: float  # blood to plasma ratio [-]
    IVDOSE: NotRequired[float]  # intravenous dose [mg]
    PODOSE: NotRequired[float]  # oral dose [mg]
    tissue_distribution: NotRequired[bool]  # ftissue and Kp, fitted


#: the compartments without a tissue store: the blood and the organs of their
#: own submodels
NO_TISSUE: frozenset[str] = frozenset({"ve", "ar", "li", "ki", "gu", "po", "hv"})

SUBSTANCES_BODY: dict[str, SubstanceBody] = {
    "los": {
        "name": "losartan",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # Molecular weight (free acid; the tablets are losartan potassium,
        # 461.0 g/mole)
        "Mr": 422.911,  # [g/mole]
        # blood to plasma ratio 0.53 [Lo1995]
        "BP": 0.53,
        # iv doses [mg]
        "IVDOSE": 0,  # venous blood
        # oral absorption
        "PODOSE": 0,  # dose
        # the exchange flow (ftissue) and the tissue/blood partition
        # coefficient (Kp) are fitted, see `fits.FITS`; Vss 34 l
        # (0.45 l/kg) [Lo1995]; fu 1.3 %, albumin bound
        "tissue_distribution": True,
    },
    "e3174": {
        "name": "E3174",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # Molecular weight
        "Mr": 436.894,  # [g/mole]
        # blood to plasma ratio 0.58 [Lo1995]
        "BP": 0.58,
        # iv doses [mg]
        "IVDOSE": 0,  # venous blood
        # the exchange flow (ftissue) and the tissue/blood partition
        # coefficient (Kp) are fitted, see `fits.FITS`; Vss 10 l
        # (0.13 l/kg) [Lo1995]; fu 0.2 %, albumin bound
        "tissue_distribution": True,
    },
    "l158": {
        "name": "other losartan metabolites",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # Molecular weight of the hydroxylated losartan (L-158,795/796,
        # C22H23ClN6O2); the amounts of the data are molar (% of the dose)
        "Mr": 438.9,  # [g/mole]
        # not measured in blood, assumed to distribute as plasma water
        "BP": 1.0,
    },
}

# -----------------------------------------------------------------------------
# Submodels
# -----------------------------------------------------------------------------
COMPARTMENTS_BODY = {
    "ki": "kidney",
    "li": "liver",
    "lu": "lung",
    "gu": "gut",
    "re": "rest",
    "ar": "arterial blood",
    "ve": "venous blood",
    "po": "portal vein",
    "hv": "hepatic vein",
}

SUBMODEL_SID_DICT = {  # tissue to submodel mapping
    "ki": "KI",  # kidney
    "li": "LI",  # liver
    "gu": "GU",  # gut/intestine
}

kidney_id = "losartan_kidney"
liver_id = "losartan_liver"
gut_id = "losartan_intestine"

_m.external_model_definitions = [
    ExternalModelDefinition(
        sid="kidney",
        name="kidney model definition",
        source=f"{kidney_id}.xml",
        modelRef=kidney_id,
    ),
    ExternalModelDefinition(
        sid="liver",
        name="liver model definition",
        source=f"{liver_id}.xml",
        modelRef=liver_id,
    ),
    ExternalModelDefinition(
        sid="gut", name="gut model definition", source=f"{gut_id}.xml", modelRef=gut_id
    ),
]

_m.submodels = [
    Submodel(sid=SUBMODEL_SID_DICT["ki"], modelRef="kidney", name="kidney submodel"),
    Submodel(sid=SUBMODEL_SID_DICT["li"], modelRef="liver", name="liver submodel"),
    Submodel(sid=SUBMODEL_SID_DICT["gu"], modelRef="gut", name="gut submodel"),
]

# -------------------------------------------------------------------------------------------------
# Compartments
# -------------------------------------------------------------------------------------------------
_m.compartments = [
    Compartment(
        "Vre",
        metaId="meta_Vre",
        value=1,
        unit=U.liter,
        constant=False,
        name="rest of body",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["re"],
    ),
    Compartment(
        "Vgu",
        metaId="meta_Vgu",
        value=1,
        unit=U.liter,
        constant=False,
        name="gut",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["gu"],
    ),
    Compartment(
        "Vki",
        metaId="meta_Vki",
        value=1,
        unit=U.liter,
        constant=False,
        name="kidney",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
    ),
    Compartment(
        "Vli",
        metaId="meta_Vli",
        value=1,
        unit=U.liter,
        constant=False,
        name="liver",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
    ),
    Compartment(
        "Vlu",
        metaId="meta_Vlu",
        value=1,
        unit=U.liter,
        constant=False,
        name="lung",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["lu"],
    ),
    Compartment(
        "Vve",
        metaId="meta_Vve",
        value=1,
        unit=U.liter,
        constant=False,
        name="venous blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ve"],
    ),
    Compartment(
        "Var",
        metaId="meta_Var",
        value=1,
        unit=U.liter,
        constant=False,
        name="arterial blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ar"],
    ),
    Compartment(
        "Vurine",
        metaId="meta_Vurine",
        value=1,
        unit=U.liter,
        constant=True,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["urine"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["stomach"],
    ),
    Compartment(
        "Vpo",
        metaId="meta_Vpo",
        value=1,
        unit=U.liter,
        constant=False,
        name="portal blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["po"],
    ),
    Compartment(
        "Vhv",
        metaId="meta_Vhv",
        value=1,
        unit=U.liter,
        constant=False,
        name="hepatic venous blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["hv"],
    ),
]

# create blood and tissue compartments of the organs
for cid in COMPARTMENTS_BODY:
    if cid not in ["ar", "ve", "po", "hv"]:
        _m.compartments.extend(
            [
                Compartment(
                    f"V{cid}_tissue",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} tissue",
                    metaId=f"meta_V{cid}_tissue",
                    port=True,
                    annotations=annotations.compartments["parenchyma"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                ),
                Compartment(
                    f"V{cid}_blood",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} blood",
                    metaId=f"meta_V{cid}_blood",
                    port=True,
                    annotations=annotations.compartments["blood"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                ),
            ]
        )

# replace volumes of submodels
_m.replaced_elements.extend(
    [
        # kidney
        ReplacedElement(
            sid="Vki_tissue_RE",
            metaId="Vki_tissue_RE",
            elementRef="Vki_tissue",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vki{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vki_blood_RE",
            metaId="Vki_blood_RE",
            elementRef="Vki_blood",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vurine_RE",
            metaId="Vurine_RE",
            elementRef="Vurine",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vurine{PORT_SUFFIX}",
        ),
        # liver
        ReplacedElement(
            sid="Vli_tissue_RE",
            metaId="Vli_tissue_RE",
            elementRef="Vli_tissue",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"Vli{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vli_blood_RE",
            metaId="Vli_blood_RE",
            elementRef="Vli_blood",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        # intestine
        ReplacedElement(
            sid="Vgu_blood_RE",
            metaId="Vgu_blood_RE",
            elementRef="Vgu_blood",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vgu_RE",
            metaId="Vgu_RE",
            elementRef="Vgu",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vgu{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vfeces_RE",
            metaId="Vfeces_RE",
            elementRef="Vfeces",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vfeces{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vstomach_RE",
            metaId="Vstomach_RE",
            elementRef="Vstomach",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vstomach{PORT_SUFFIX}",
        ),
    ]
)

# -------------------------------------------------------------------------------------------------
# Species
# -------------------------------------------------------------------------------------------------
_m.species = []
for sid, sdict in SUBSTANCES_BODY.items():
    for cid, cname in COMPARTMENTS_BODY.items():
        if cid in ["ve", "ar", "po", "hv"]:
            sid_ex = f"C{cid}_{sid}"
            cid_ex = f"V{cid}"
        else:
            # blood of the organ
            sid_ex = f"C{cid}_blood_{sid}"
            cid_ex = f"V{cid}_blood"

        _m.species.append(
            Species(
                sid_ex,
                metaId=f"meta_{sid_ex}",
                initialConcentration=sdict["cinit"],
                compartment=cid_ex,
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} ({cname} blood)",
                hasOnlySubstanceUnits=False,
                port=True,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )

    # tissue species
    for cid, cname in COMPARTMENTS_BODY.items():
        if sdict.get("tissue_distribution") and cid not in NO_TISSUE:
            sid_ex = f"C{cid}_{sid}"
            cid_ex = f"V{cid}_tissue"

            _m.species.extend(
                [
                    Species(
                        sid_ex,
                        metaId=f"meta_{sid_ex}",
                        initialConcentration=sdict["cinit"],
                        compartment=cid_ex,
                        substanceUnit=sdict["unit"],
                        name=f"{sdict['name']} ({cname})",
                        hasOnlySubstanceUnits=False,
                        port=True,
                        sboTerm=SBO.SIMPLE_CHEMICAL,
                        annotations=annotations.species[sid],
                    )
                ]
            )

    # urine metabolites dummy species
    if sid in ["los", "e3174", "l158"]:
        _m.species.extend(
            [
                Species(
                    f"Aurine_{sid}",
                    metaId=f"meta_Aurine_{sid}",
                    initialConcentration=0,
                    compartment="Vurine",
                    substanceUnit=sdict["unit"],
                    name=f"{sdict['name']} (urine)",
                    hasOnlySubstanceUnits=True,
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                    annotations=annotations.species[sid],
                ),
            ]
        )

    # feces metabolites
    if sid in ["los", "e3174", "l158"]:
        _m.species.append(
            Species(
                f"Afeces_{sid}",
                metaId=f"meta_Afeces_{sid}",
                initialConcentration=0,
                compartment="Vfeces",
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} (feces)",
                hasOnlySubstanceUnits=True,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )

for sid in ["los", "e3174", "l158"]:
    sdict = SUBSTANCES_BODY[sid]
    _m.species.append(
        Species(
            f"Cgu_{sid}",
            metaId=f"meta_Cgu_{sid}",
            initialConcentration=0,
            compartment="Vgu",
            substanceUnit=sdict["unit"],
            name=f"{sdict['name']} (gut lumen)",
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
        )
    )

# plasma concentrations: the blood concentration over the blood to plasma ratio
PLASMA_SPECIES: dict[str, tuple[str, str]] = {
    "los": ("losartan (venous plasma)", "Cve_los/BP_los"),
    "e3174": ("E3174 (venous plasma)", "Cve_e3174/BP_e3174"),
    "l158": ("other losartan metabolites (venous plasma)", "Cve_l158/BP_l158"),
    "los_e3174": (
        "losartan and E3174 (venous plasma)",
        "Cve_los/BP_los + Cve_e3174/BP_e3174",
    ),
    "total": (
        "losartan and its metabolites (venous plasma)",
        "Cve_los/BP_los + Cve_e3174/BP_e3174 + Cve_l158/BP_l158",
    ),
}
for key, (name, _) in PLASMA_SPECIES.items():
    _m.species.append(
        Species(
            f"Cveplasma_{key}",
            metaId=f"meta_Cveplasma_{key}",
            initialConcentration=0,
            compartment="Vve",
            substanceUnit=U.mmole,
            name=name,
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species.get(key, []),
        )
    )

# replace species
replaced_species = {
    "li": ["los", "e3174", "l158"],
    "ki": ["los", "e3174", "l158"],
    "gu": ["los"],
}

# blood to external species
for tkey, skey_list in replaced_species.items():
    for skey in skey_list:
        _m.replaced_elements.extend(
            [
                ReplacedElement(
                    sid=f"C{tkey}_blood_{skey}_RE",
                    metaId=f"C{tkey}_blood_{skey}_RE",
                    elementRef=f"C{tkey}_blood_{skey}",
                    submodelRef=SUBMODEL_SID_DICT[tkey],
                    portRef=f"{skey}_ext{PORT_SUFFIX}",
                ),
            ]
        )

# urine species replacements
for skey in ["los", "e3174", "l158"]:
    _m.replaced_elements.append(
        ReplacedElement(
            sid=f"Aurine_{skey}_RE",
            metaId=f"Aurine_{skey}_RE",
            elementRef=f"Aurine_{skey}",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"{skey}_urine{PORT_SUFFIX}",
        )
    )

# feces species replacements
for skey in ["los", "e3174", "l158"]:
    _m.replaced_elements.append(
        ReplacedElement(
            sid=f"Afeces_{skey}_RE",
            metaId=f"Afeces_{skey}_RE",
            elementRef=f"Afeces_{skey}",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"{skey}_feces{PORT_SUFFIX}",
        )
    )

# intestine species replacements
for skey in ["los", "e3174", "l158"]:
    _m.replaced_elements.extend(
        [
            ReplacedElement(
                sid=f"Cgu_{skey}_RE",
                metaId=f"Cgu_{skey}_RE",
                elementRef=f"Cgu_{skey}",
                submodelRef=SUBMODEL_SID_DICT["gu"],
                portRef=f"{skey}_lumen{PORT_SUFFIX}",
            ),
        ]
    )

# liver species replacement: the bile of the metabolites flows into the gut lumen
for skey in ["e3174", "l158"]:
    _m.replaced_elements.extend(
        [
            ReplacedElement(
                sid=f"Cgu_{skey}_RE2",
                metaId=f"Cgu_{skey}_RE2",
                elementRef=f"Cgu_{skey}",
                submodelRef=SUBMODEL_SID_DICT["li"],
                portRef=f"{skey}_lumen{PORT_SUFFIX}",
            ),
        ]
    )

# -------------------------------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------------------------------
_m.parameters.extend(
    [
        # whole body data
        Parameter(
            "BW",
            75,
            U.kg,
            constant=True,
            name="body weight",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "BW_ref",
            75,
            U.kg,
            constant=True,
            name="reference body weight",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=(
                "Body weight at which the exchange flows of the tissue "
                "distribution (`ftissue_<sid>`) are given; they scale with the "
                "body weight as the volumes, the blood flows and the renal "
                "clearance do."
            ),
        ),
        Parameter(
            "COBW",
            1.548,
            U.ml_per_s_kg,
            constant=True,
            name="cardiac output per bodyweight",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_cardiac_function",
            1.0,
            U.dimensionless,
            constant=False,
            name="cardiac function",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            1.0: normal function;
            <1.0: reduced function as reduced cardiac output;
            """,
        ),
        Parameter(
            "CO",
            108.33,
            U.ml_per_s,
            constant=False,
            name="cardiac output",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "QC",
            108.33 / 1000 * 60,
            U.l_per_min,
            constant=False,
            name="cardiac output",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue volumes
        Parameter(
            "Fblood",
            0.02,
            U.dimensionless,
            constant=False,
            name="blood fraction of organ volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_GFR_Qki",
            101 / 1.32354,
            U.ml_per_l,
            constant=True,
            name="relationship GFR and Qki",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            GFR per renal blood flow (76.3 ml/l). Calibrated so that the healthy
            GFR is 101 ml/min at the renal blood flow of 1.32354 l/min of 75 kg,
            the normal GFR of the classes of renal impairment
            (`pkpdsim.impairment`).
            """,
        ),
        Parameter(
            "GFR_base",
            101,
            U.ml_per_min,
            constant=False,
            name="glomerular filtration rate (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The GFR of healthy kidneys, which scales with the renal blood flow
            and so with the body weight; the GFR is `KI__f_renal_function`
            times `GFR_base`.
            """,
        ),
        # fractional tissue volumes
        Parameter(
            "FVgu",
            0.0171,
            U.l_per_kg,
            constant=True,
            name="gut fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVki",
            0.0044,
            U.l_per_kg,
            constant=True,
            name="kidney fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVli",
            0.0210,
            U.l_per_kg,
            constant=True,
            name="liver fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVlu",
            0.0076,
            U.l_per_kg,
            constant=True,
            name="lung fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVre",
            0,
            U.l_per_kg,
            constant=False,
            name="rest of body fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="Set by its assignment rule: the rest of the body weight.",
        ),
        # calculated based on other tissues
        Parameter(
            "FVve",
            0.0514,
            U.l_per_kg,
            constant=True,
            name="venous fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVar",
            0.0257,
            U.l_per_kg,
            constant=True,
            name="arterial fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVpo",
            0.001,
            U.l_per_kg,
            constant=True,
            name="portal fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVhv",
            0.001,
            U.l_per_kg,
            constant=True,
            name="hepatic venous fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue blood flows
        Parameter(
            "FQgu",
            0.18,
            U.dimensionless,
            constant=True,
            name="gut fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
        Assumed FQgu is the sum of FQgu + FQsp + FQpa for correct portal blood flow.
        with FQgu=0.146, FQsp=0.017, FQpa=0.017
        """,
        ),
        Parameter(
            "FQki",
            0.190,
            U.dimensionless,
            constant=True,
            name="kidney fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQh",
            0.215,
            U.dimensionless,
            constant=True,
            name="hepatic fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQlu",
            1,
            U.dimensionless,
            constant=True,
            name="lung fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQre",
            0,
            U.dimensionless,
            constant=False,
            name="rest of body fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="Set by its assignment rule: the rest of the cardiac output.",
        ),
        Parameter(
            "conversion_min_per_day",
            1440,
            U.min_per_day,
            constant=True,
            name="Conversion factor day to min",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_cirrhosis",
            0,
            U.dimensionless,
            constant=True,
            name="severity of cirrhosis [0, 0.95]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            annotations=[
                (BQB.IS, "NCIT:C2951"),  # Cirrhosis
                (BQB.IS, "efo/0001422"),  # cirrhosis of liver
            ],
            notes="""
            healthy = 0.0; measures decrease of cirrhosis often via CPT score.
            "Control": 0,
            "Mild cirrhosis": 0.399,  CPT A
            "Moderate cirrhosis": 0.698,  CPT B
            "Severe cirrhosis": 0.813,  CPT C
            """,
        ),
        Parameter(
            "f_shunts",
            0,
            U.dimensionless,
            constant=False,
            name="fraction of portal venous blood shunted by the liver",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            annotations=[
                (
                    BQB.IS,
                    "hp/HP:0002629",
                ),  # Gastrointestinal arteriovenous malformation
            ],
        ),
        Parameter(
            "f_tissue_loss",
            0,
            U.dimensionless,
            constant=False,
            name="fraction of lost parenchymal liver volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    ]
)


for sid, sdict in SUBSTANCES_BODY.items():
    if "PODOSE" in sdict:
        # Parameter for intestine model
        _m.parameters.extend(
            [
                Parameter(
                    f"PODOSE_{sid}",
                    sdict["PODOSE"],
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"oral dose {sid} [mg]",
                ),
            ]
        )
        _m.replaced_elements.extend(
            [
                ReplacedElement(
                    sid=f"PODOSE_{sid}_RE",
                    metaId=f"PODOSE_{sid}_RE",
                    elementRef=f"PODOSE_{sid}",
                    submodelRef=SUBMODEL_SID_DICT["gu"],
                    portRef=f"PODOSE_{sid}{PORT_SUFFIX}",
                ),
            ]
        )

_m.rules = [
    *_m.rules,
    AssignmentRule("f_shunts", "f_cirrhosis", unit=U.dimensionless),
    AssignmentRule("f_tissue_loss", "f_cirrhosis", unit=U.dimensionless),
]

replaced_parameters = {
    "li": [],
    "gu": ["Mr_los"],
    "ki": ["GFR_base", "BP_los", "BP_e3174", "BP_l158"],
}

for ckey, pids in replaced_parameters.items():
    for pid in pids:
        _m.replaced_elements.append(
            ReplacedElement(
                sid=f"{pid}_{ckey}_RE",
                metaId=f"{pid}_{ckey}_RE",
                elementRef=f"{pid}",
                submodelRef=SUBMODEL_SID_DICT[ckey],
                portRef=f"{pid}{PORT_SUFFIX}",
            )
        )

# species specific parameters
for sid, sdict in SUBSTANCES_BODY.items():
    _m.parameters.extend(
        [
            # molecular weights
            Parameter(
                f"Mr_{sid}",
                sdict["Mr"],
                U.g_per_mole,
                constant=True,
                name=f"Molecular weight {sid}",
                sboTerm=SBO.MOLECULAR_MASS,
            ),
            # blood to plasma ratio, for the plasma concentrations
            Parameter(
                f"BP_{sid}",
                sdict["BP"],
                U.dimensionless,
                constant=True,
                name=f"blood to plasma ratio {sid}",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            ),
        ]
    )

    # tissue distribution
    if sdict.get("tissue_distribution"):
        _m.parameters.extend(
            [
                Parameter(
                    f"ftissue_{sid}",
                    **FITS.parameter(f"ftissue_{sid}", unit=U.l_per_min),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"tissue distribution {sid.upper()}",
                    notes="Exchange flow at the reference body weight `BW_ref`.",
                ),
                Parameter(
                    f"Kp_{sid}",
                    **FITS.parameter(f"Kp_{sid}", unit=U.dimensionless),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"tissue/blood partition coefficient {sid.upper()}",
                ),
            ]
        )

    # dosing
    if "IVDOSE" in sdict:
        _m.parameters.extend(
            [
                Parameter(
                    f"IVDOSE_{sid}",
                    sdict["IVDOSE"],
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"IV bolus dose {sid} [mg]",
                ),
                # iv kinetics after application
                Parameter(
                    f"ti_{sid}",
                    10,
                    U.second,
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"injection time {sid}",
                ),
                Parameter(
                    f"Ki_{sid}",
                    0.02,
                    U.per_min,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Ki [1/min] injection {sid}",
                ),
                # continuous infusion
                Parameter(
                    f"Ri_{sid}",
                    0,
                    U.mg_per_min,
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Ri [mg/min] rate of infusion {sid}",
                ),
                Parameter(
                    f"cum_dose_{sid}",
                    0,
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Cumulative dose due to infusion {sid}",
                ),
            ]
        )

# -------------------------------------------------------------------------------------------------
# AssignmentRules
# -------------------------------------------------------------------------------------------------
_m.rules = [
    *_m.rules,
    # Rest body volume
    AssignmentRule(
        "FVre",
        "1.0 l_per_kg - (FVgu + FVki + FVli + FVlu + FVve + FVar + FVpo + FVhv)",
        U.l_per_kg,
    ),
    # Rest body perfusion
    AssignmentRule("FQre", "1.0 dimensionless - (FQki + FQh)", U.dimensionless),
    # cardiac output (depending on heart rate and bodyweight)
    AssignmentRule("CO", "f_cardiac_function * BW * COBW", U.ml_per_s),
    # cardiac output (depending on bodyweight)
    AssignmentRule("QC", "CO/1000 ml_per_l * 60 s_per_min", U.l_per_min),
    # volumes
    AssignmentRule("Vgu", "BW*FVgu", U.liter),
    AssignmentRule("Vki", "BW*FVki", U.liter),
    AssignmentRule("Vli", "BW*FVli", U.liter),
    AssignmentRule("Vlu", "BW*FVlu", U.liter),
    AssignmentRule("Vre", "BW*FVre", U.liter),
    # venous and arterial blood volume (corrected for tissue blood volumes)
    AssignmentRule(
        "Vve",
        "BW*FVve - FVve/(FVar+FVve) * BW * Fblood * (1 l_per_kg - FVve - FVar)",
        U.liter,
    ),
    AssignmentRule(
        "Var",
        "BW*FVar - FVar/(FVar+FVve) * BW * Fblood * (1 l_per_kg - FVve - FVar)",
        U.liter,
    ),
    AssignmentRule(
        "Vpo",
        "BW*FVpo - FVpo/(FVar+FVve+FVpo+FVhv) "
        "* BW * Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv))",
        U.liter,
    ),
    AssignmentRule(
        "Vhv",
        "BW*FVhv - FVhv/(FVar+FVve+FVpo+FVhv) "
        "* BW * Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv))",
        U.liter,
    ),
    # blood flows
    AssignmentRule("Qgu", "QC*FQgu", U.l_per_min, name="gut blood flow"),
    AssignmentRule("Qki", "QC*FQki", U.l_per_min, name="kidney blood flow"),
    AssignmentRule(
        "Qh", "QC*FQh", U.l_per_min, name="hepatic (venous side) blood flow"
    ),
    AssignmentRule("Qha", "Qh - Qgu", U.l_per_min, name="hepatic artery blood flow"),
    AssignmentRule("Qlu", "QC*FQlu", U.l_per_min, name="lung blood flow"),
    AssignmentRule("Qre", "QC*FQre", U.l_per_min, name="rest of body blood flow"),
    AssignmentRule("Qpo", "Qgu", U.l_per_min, name="portal blood flow"),
    # renal function
    AssignmentRule(
        "GFR_base",
        "f_GFR_Qki * Qki",
        U.ml_per_min,
        name="glomerular filtration rate (healthy)",
    ),
]

# Volumes for explicit tissue models
for cid, cname in COMPARTMENTS_BODY.items():
    if cid not in ["ve", "ar", "po", "hv"]:
        _m.rules.append(
            # blood volume of the organ
            AssignmentRule(
                f"V{cid}_blood",
                value=f"V{cid} * Fblood",
                unit=U.liter,
                name=f"blood volume of {cname}",
            ),
        )

        if cid == "li":
            _m.rules.append(
                # Cirrhosis: Adjustment of fractional liver tissue volume
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - f_tissue_loss) * (1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                ),
            )
        else:
            _m.rules.append(
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                ),
            )

for sid, sdict in SUBSTANCES_BODY.items():
    # injection
    if "IVDOSE" in sdict:
        _m.rules.extend(
            [
                AssignmentRule(
                    f"Ki_{sid}",
                    f"0.693 dimensionless/ti_{sid} * 60 s_per_min",
                    U.per_min,
                    name="injection rate IV",
                    notes="""
                log(2) ~ 0.693 used for conversion of half-life injection time to
                injection rate.
                """,
                ),
            ]
        )

# --------------------------------------------------------------------------------------------------
# Reactions
# --------------------------------------------------------------------------------------------------
_m.reactions = []

for sid, sdict in SUBSTANCES_BODY.items():
    sname = sdict["name"]

    # --------------------
    # tissue distribution
    # --------------------
    if sdict.get("tissue_distribution"):
        for cid in COMPARTMENTS_BODY:
            if cid not in NO_TISSUE:
                _m.reactions.append(
                    Reaction(
                        sid=f"transport_{cid}_{sid}",
                        name=f"transport {sname}",
                        formula=(
                            f"ftissue_{sid} * BW/BW_ref * (C{cid}_blood_{sid}*Kp_{sid} - C{cid}_{sid})",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} <-> C{cid}_{sid}",
                        compartment=f"V{cid}_tissue",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                )

    if "IVDOSE" in sdict:
        _m.reactions.extend(
            [
                # --------------------
                # iv dose application
                # --------------------
                Reaction(
                    sid=f"iv_{sid}",
                    name=f"iv {sname}",
                    formula=(f"Ki_{sid}*IVDOSE_{sid}/Mr_{sid}", U.mmole_per_min),
                    equation=f"-> Cve_{sid}",
                    sboTerm=SBO.TRANSPORT_REACTION,
                    annotations=[(BQB.IS, "NCIT:C38276")],
                    compartment="Vve",
                ),
            ]
        )

    for cid, cname in COMPARTMENTS_BODY.items():
        if cid in ["ve", "ar", "po", "hv"]:
            continue
        # --------------------
        # ve -> lung -> ar
        # --------------------
        if cid == "lu":
            rid_in = f"Flow_ve_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"

            rid_out = f"Flow_{cid}_ar_{sid}"
            name_out = f"outflow {cname} {sname}"

            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Cve_{sid}", U.mmole_per_min),
                        equation=f"Cve_{sid} -> C{cid}_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} -> Car_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> organ -> ve
        # --------------------
        if cid in ["ki", "re"]:
            rid_in = f"Flow_ar_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"
            rid_out = f"Flow_{cid}_ve_{sid}"
            name_out = f"outflow {cname} {sname}"

            # only distribution in the blood of the organ
            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> organ -> po
        # --------------------
        if cid in ["gu"]:
            rid_in = f"Flow_ar_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"
            rid_out = f"Flow_{cid}_po_{sid}"
            name_out = f"outflow {cname} {sname}"

            _m.reactions.extend(
                [
                    # only distribution in the blood of the organ
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} -> Cpo_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> li
        # po -> li
        # li -> hv -> ve
        # --------------------
        if cid == "li":
            # liver
            _m.reactions.extend(
                [
                    Reaction(
                        sid=f"Flow_arli_li_{sid}",
                        name=f"arterial inflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qha*Car_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Car_{sid} -> Cli_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # shunted arterial flow
                    Reaction(
                        sid=f"Flow_arli_hv_{sid}",
                        name="flow arterial shunts",
                        formula=(f"f_shunts*Qha*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # (unshunted) portal flow
                    Reaction(
                        sid=f"Flow_po_li_{sid}",
                        name=f"outflow po {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qpo*Cpo_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cpo_{sid} -> Cli_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # portal shunts
                    Reaction(
                        sid=f"Flow_po_hv_{sid}",
                        name="flow portal shunts",
                        formula=(f"f_shunts*Qpo*Cpo_{sid}", U.mmole_per_min),
                        equation=f"Cpo_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_li_hv_{sid}",
                        name=f"outflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*(Qpo+Qha)*Cli_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cli_blood_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_hv_ve_{sid}",
                        name=f"outflow hepatic vein {sname}",
                        formula=(f"Qh*Chv_{sid}", U.mmole_per_min),
                        equation=f"Chv_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )

# --------------------------------------------------------------------------------------------------
# RateRules
# --------------------------------------------------------------------------------------------------
_m.rate_rules = []

for sid, sdict in SUBSTANCES_BODY.items():
    if "IVDOSE" in sdict:
        _m.rate_rules.extend(
            [
                # injection of dose
                RateRule(
                    f"IVDOSE_{sid}", f"-iv_{sid}*Mr_{sid} + Ri_{sid}", U.mg_per_min
                ),
                # cumulative infusion dose
                RateRule(f"cum_dose_{sid}", f"Ri_{sid}", U.mg_per_min),
            ]
        )

# --------------------------------------------------------------------------------------------------
# Metabolic ratios
# --------------------------------------------------------------------------------------------------
_m.parameters.extend(
    [
        Parameter(
            sid="mr_e3174_los_plasma",
            value=0,
            unit=U.dimensionless,
            notes="""Metabolic ratio E3174/losartan""",
            constant=False,
        ),
        Parameter(
            sid="mr_e3174_los_urine",
            value=0,
            unit=U.dimensionless,
            notes="""Metabolic ratio E3174/losartan""",
            constant=False,
        ),
        Parameter(
            sid="mr_e3174_los_feces",
            value=0,
            unit=U.dimensionless,
            notes="""Metabolic ratio E3174/losartan""",
            constant=False,
        ),
        Parameter(
            sid="mr_los_e3174_plasma",
            value=0,
            unit=U.dimensionless,
            notes="""Metabolic ratio losartan/E3174""",
            constant=False,
        ),
        Parameter(
            sid="mr_los_e3174_urine",
            value=0,
            unit=U.dimensionless,
            notes="""Metabolic ratio losartan/E3174""",
            constant=False,
        ),
        Parameter(
            sid="mr_los_e3174_feces",
            value=0,
            unit=U.dimensionless,
            notes="""Metabolic ratio losartan/E3174""",
            constant=False,
        ),
        Parameter(
            sid="Aurine_los_e3174",
            value=0,
            unit=U.mmole,
            notes="""Sum losartan and E3174 urine""",
            constant=False,
        ),
        Parameter(
            sid="Aurine_total",
            value=0,
            unit=U.mmole,
            notes="""Sum losartan and E3174 and L158 urine""",
            constant=False,
        ),
        Parameter(
            sid="Afeces_total",
            value=0,
            unit=U.mmole,
            notes="""Sum losartan and E3174 and L158 feces""",
            constant=False,
        ),
    ]
)
_m.rules.extend(
    [
        AssignmentRule(
            "mr_e3174_los_plasma",
            "Cveplasma_e3174/(Cveplasma_los + 1E-12 mM)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "mr_e3174_los_urine",
            "Aurine_e3174/(Aurine_los +  1E-12 mmole)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "mr_e3174_los_feces",
            "Afeces_e3174/(Afeces_los +  1E-12 mmole)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "mr_los_e3174_plasma",
            "Cveplasma_los/(Cveplasma_e3174 + 1E-12 mM)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "mr_los_e3174_urine",
            "Aurine_los/(Aurine_e3174 +  1E-12 mmole)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "mr_los_e3174_feces",
            "Afeces_los/(Afeces_e3174 +  1E-12 mmole)",
            unit=U.dimensionless,
        ),
        AssignmentRule("Aurine_los_e3174", "Aurine_los + Aurine_e3174", unit=U.mmole),
        *(
            AssignmentRule(
                f"Cveplasma_{key}",
                formula,
                unit=U.mM,
                name=name,
                notes="The blood concentration over the blood to plasma ratio (BP).",
            )
            for key, (name, formula) in PLASMA_SPECIES.items()
        ),
        AssignmentRule(
            "Aurine_total", "Aurine_los + Aurine_e3174 + Aurine_l158", unit=U.mmole
        ),
        AssignmentRule(
            "Afeces_total", "Afeces_los + Afeces_e3174 + Afeces_l158", unit=U.mmole
        ),
    ]
)

# coupling of the PK and the RAAS model: the E3174 of the RAAS model, which is
# merged into the whole-body model, is the E3174 of the venous plasma, the
# concentration at the AT1 receptors of the vessels, the adrenal gland and the
# kidney
_m.rules.append(AssignmentRule("e3174", "Cveplasma_e3174", unit=U.mM))

model_body: Model = _m


if __name__ == "__main__":
    from pkdb_models.models.losartan import MODEL_BASE_PATH
    from pkdb_models.models.losartan.models.model_raas import model_raas

    result = create_model(
        filepath=MODEL_BASE_PATH / f"{model_body.sid}.xml",
        model=[model_raas, model_body],
        sbml_level=3,
        sbml_version=2,
    )
