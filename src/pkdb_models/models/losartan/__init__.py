from pathlib import Path

LOSARTAN_PATH = Path(__file__).parent

MODEL_BASE_PATH: Path = LOSARTAN_PATH / "models" / "results" / "models"
MODEL_PATH: Path = MODEL_BASE_PATH / "losartan_body_flat.xml"

RESULTS_PATH = LOSARTAN_PATH / "results"
RESULTS_PATH_SIMULATION = RESULTS_PATH / "simulation"
RESULTS_PATH_FIT = RESULTS_PATH / "fit"

# DATA_PATH_BASE = LOSARTAN_PATH.parents[3] / "pkdb_data" / "studies"

DATA_PATH_BASE = LOSARTAN_PATH / "data"

DATA_PATH_LOSARTAN = DATA_PATH_BASE / "losartan"
DATA_PATHS = [
     DATA_PATH_LOSARTAN,
     DATA_PATH_BASE / "dextromethorphan",
     DATA_PATH_BASE / "caffeine",
     DATA_PATH_BASE / "chlorzoxazone",
]
