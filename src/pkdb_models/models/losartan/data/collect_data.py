"""Collect the data of the losartan studies from the PK-DB study folders.

The studies are copied from a `pkdb_data` checkout (`PKDB_DATA`, else
`../pkdb_data`, see `pkpdsim.data.studies_dir`) into the data folder of the
model, see `pkpdsim.data.collect_tsv_files`.
Besides the losartan studies, cocktail studies of other substances have
losartan groups: Kim2016 (dextromethorphan), Oh2012 and Tanaka2014 (caffeine)
and Puris2019 (chlorzoxazone).
"""

from pathlib import Path

from pkpdsim.data import collect_tsv_files
from pkpdsim.data import studies_dir as pkdb_studies_dir

from pkdb_models.models.losartan import (
    DATA_PATH_CAFFEINE,
    DATA_PATH_CHLORZOXAZONE,
    DATA_PATH_DEXTROMETHORPHAN,
    DATA_PATH_LOSARTAN,
)

#: the studies of other substances which the losartan experiments use
DEXTROMETHORPHAN_STUDIES = {"Kim2016"}
CAFFEINE_STUDIES = {"Oh2012", "Tanaka2014"}
CHLORZOXAZONE_STUDIES = {"Puris2019"}


def collect_data(studies_dir: Path | None = None) -> None:
    """Copy the data of the losartan studies into the data folder of the model.

    Args:
        studies_dir: the `studies` folder of pkdb_data, by default the one of
            the checkout `pkpdsim.data.studies_dir` finds.
    """
    if studies_dir is None:
        studies_dir = pkdb_studies_dir()
    collect_tsv_files(
        source_dir=studies_dir / "losartan", target_dir=DATA_PATH_LOSARTAN
    )
    collect_tsv_files(
        source_dir=studies_dir / "dextromethorphan",
        target_dir=DATA_PATH_DEXTROMETHORPHAN,
        filter_study=lambda study: study in DEXTROMETHORPHAN_STUDIES,
    )
    collect_tsv_files(
        source_dir=studies_dir / "caffeine",
        target_dir=DATA_PATH_CAFFEINE,
        filter_study=lambda study: study in CAFFEINE_STUDIES,
    )
    collect_tsv_files(
        source_dir=studies_dir / "chlorzoxazone",
        target_dir=DATA_PATH_CHLORZOXAZONE,
        filter_study=lambda study: study in CHLORZOXAZONE_STUDIES,
    )


if __name__ == "__main__":
    collect_data()
