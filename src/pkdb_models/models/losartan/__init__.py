"""Losartan: PBPK/PD model of losartan, its metabolite E3174 and the RAAS.

The paths of the model and its model object `MODEL`, which the commands
`run_simulation`, `run_fit`, `run_report` and `run_identifiability` of pkpdsim
run by its package `pkdb_models.models.losartan`.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim.release.settings import Author, ReleaseSettings

from pkdb_models.models.losartan.tolerances import ABSOLUTE_TOLERANCE

LOSARTAN_PATH = Path(__file__).parent

MODEL_BASE_PATH: Path = LOSARTAN_PATH / "models" / "output" / "models"
MODEL_PATH: Path = MODEL_BASE_PATH / "losartan_body_flat.xml"

DATA_PATH_BASE = LOSARTAN_PATH / "data"
DATA_PATH_LOSARTAN = DATA_PATH_BASE / "losartan"
DATA_PATH_DEXTROMETHORPHAN = DATA_PATH_BASE / "dextromethorphan"
DATA_PATH_CAFFEINE = DATA_PATH_BASE / "caffeine"
DATA_PATH_CHLORZOXAZONE = DATA_PATH_BASE / "chlorzoxazone"
DATA_PATHS = [
    DATA_PATH_LOSARTAN,
    DATA_PATH_DEXTROMETHORPHAN,
    DATA_PATH_CAFFEINE,
    DATA_PATH_CHLORZOXAZONE,
]

MODEL = PKPDModel(
    name="losartan",
    package=__name__,
    title="Losartan PBPK/PD model",
    path=LOSARTAN_PATH,
    model_path=MODEL_PATH,
    data_paths=tuple(DATA_PATHS),
    archive_path=MODEL_BASE_PATH.parent / "losartan_model.omex",
    factory=".models.factory:create_losartan_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    # the adopted fits, the factory takes their values (`models/fits.py`)
    fits=("PK_20261019_113633__d9465", "PD_20261019_134029__a3a50"),
    simulator=SimulatorSettings(
        absolute_tolerance=ABSOLUTE_TOLERANCE, relative_tolerance=1e-10
    ),
    # the resolution of the PNG figures; at 600 dpi the large figures need up to
    # 12 GB of memory to render
    figures=FigureSettings(formats=("svg",), legend_fontsize=10, fig_dpi=300),
    # the release of the model repository
    release=ReleaseSettings(
        title=(
            "Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD)"
            " model of losartan"
        ),
        repository="matthiaskoenig/losartan-model",
        doi="10.5281/zenodo.14761114",
        authors=(
            # the ORCIDs: Tensil from the Zenodo record of losartan-model,
            # which gives Myshkina the same ORCID; Myshkina from the Zenodo
            # records of the other models (apixaban, lenvatinib)
            Author("Tensil, Ennie", "0009-0001-1433-4293", "Charité Berlin"),
            Author(
                "Myshkina, Mariia",
                "0009-0000-7769-2154",
                (
                    "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver"
                ),
            ),
            Author(
                "König, Matthias",
                "0000-0003-1725-179X",
                (
                    "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering; University of Southern Denmark, Clinical Pharmacology, Pharmacy and Environmental Medicine, Department of Public Health; University Hospital Schleswig-Holstein, Campus Lübeck, First Department of Medicine, Lübeck, Germany"
                ),
            ),
        ),
        keywords=(
            "model",
            "losartan",
            "E3174",
            "angiotensin II receptor blocker",
            "renin-angiotensin-aldosterone system",
            "pharmacokinetics",
            "pharmacodynamics",
            "PBPK/PD",
            "physiologically based pharmacokinetic/pharmacodynamic model",
            "SBML",
        ),
        docker_image="matthiaskoenig/losartan",
        funding=(
            "Mariia Myshkina was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within Priority Programme SPP 2311 by grant number 465194077 (Subproject SimLivA). "
            "Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B)."
        ),
        description=(
            "This repository provides the losartan physiologically based"
            " pharmacokinetic/pharmacodynamic (PBPK/PD) model: the oral"
            " absorption of losartan, its distribution, its hepatic metabolism"
            " to the active metabolite E3174 (CYP2C9) and to inactive"
            " metabolites, the renal and biliary excretion of losartan and its"
            " metabolites, and the blockade of the AT1 receptor by E3174 with"
            " its effects on renin, angiotensin, aldosterone and the blood"
            " pressure."
        ),
        first_year=2024,
    ),
)
"""The losartan model: its commands, see `pkpdsim.model`."""
