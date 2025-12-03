from pathlib import Path

from pkdb_models.models.data import collect_tsv_files

def collect_losartan_data():
    common_parent: Path = Path(__file__).parents[5]
    source_dir = common_parent / "pkdb_data" / "studies" / "losartan"
    target_dir = Path(__file__).parent / "losartan"

    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)

    # collect  (dextromethorphan)
    def is_dextromethorphan(study_name) -> bool:
        return study_name in ["Kim2016"]

    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "dextromethorphan",
        target_dir=Path(__file__).parent / "dextromethorphan",
        filter_study=is_dextromethorphan,
    )

    # collect  (caffeine)
    def is_caffeine(study_name) -> bool:
        return study_name in ["Tanaka2014", "Oh2012"]

    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "caffeine",
        target_dir=Path(__file__).parent / "caffeine",
        filter_study=is_caffeine,
    )

    # collect  (chlorzoxazone)
    def is_chlorzoxazone(study_name) -> bool:
        return study_name in ["Puris2019", "Kim2016"]

    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "chlorzoxazone",
        target_dir=Path(__file__).parent / "chlorzoxazone",
        filter_study=is_chlorzoxazone,
    )

if __name__ == "__main__":
    collect_losartan_data()

