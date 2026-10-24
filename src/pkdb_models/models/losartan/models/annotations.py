"""Annotations of the losartan models: the model and its species.

The compartments are the shared ones of `pkpdsim.annotations`.
"""

from pkpdsim.annotations import compartments
from sbmlutils.metadata import BQB

__all__ = ["compartments", "model", "species"]

model = [
    # taxonomy
    (BQB.HAS_TAXON, "taxonomy/9606"),  # human
    (BQB.HAS_TAXON, "snomedct/337915000"),  # human
    # modelling approach
    (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model
    # biological process explained by model (GO/NCIT)
    (BQB.HAS_PROPERTY, "NCIT:C66869"),  # losartan
    (BQB.HAS_PROPERTY, "NCIT:C29165"),  # Losartan Potassium
    (BQB.HAS_PROPERTY, "NCIT:C66930"),  # Angiotensin II Receptor Antagonist
    # model relevance to a particular area
    (BQB.HAS_PROPERTY, "NCIT:C3117"),  # Hypertension
    (BQB.HAS_PROPERTY, "NCIT:C2985"),  # Diabetes Mellitus
    (BQB.HAS_PROPERTY, "NCIT:C99532"),  # Diabetes Therapy
    # model reference
    (BQB.IS_DESCRIBED_BY, "doi/10.3390/pharmaceutics18020262"),
    (BQB.IS_DESCRIBED_BY, "doi/10.36903/physiome.31359823"),
    (BQB.IS_DESCRIBED_BY, "doi/10.5281/zenodo.14761114"),
]

species = {
    "los": [  # losartan
        (BQB.IS, "pubchem.compound/11751549"),
        (BQB.IS, "chebi/CHEBI:6541"),
        (BQB.IS, "snomedct/108582002"),
        (BQB.IS, "NCIT:C66869"),
    ],
    "e3174": [  # E3174 (EXP3174), active carboxylic acid metabolite
        (BQB.IS, "pubchem.compound/108185"),
        (BQB.IS, "chebi/CHEBI:74125"),
    ],
    "l158": [],
    "anggen": [
        (BQB.IS, "chebi/CHEBI:2720"),
        (BQB.IS, "NCIT:C136723"),
        (BQB.IS, "uniprot/P01019"),
    ],
    "ang1": [
        (BQB.IS, "chebi/CHEBI:147350"),
        (BQB.IS, "NCIT:C248"),
    ],
    "ang2": [
        (BQB.IS, "chebi/CHEBI:58506"),
        (BQB.IS, "NCIT:C107562"),
    ],
    "ald": [
        (BQB.IS, "chebi/CHEBI:27584"),
        (BQB.IS, "NCIT:C219"),
    ],
    "ren": [
        (BQB.IS, "uniprot/P00797"),
        (BQB.IS, "NCIT:C113591"),
    ],
    "ace": [
        (BQB.IS, "uniprot/P12821"),
    ],
}
