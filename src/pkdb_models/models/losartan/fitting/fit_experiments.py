"""Parameter fit problems for losartan."""
from typing import Dict, List
from sbmlsim.fit.helpers import f_fitexp, filter_empty
from sbmlutils.console import console
from sbmlutils.log import get_logger

from sbmlsim.fit import FitExperiment, FitMapping

from pkdb_models.models.losartan import LOSARTAN_PATH, DATA_PATHS
from pkdb_models.models.losartan.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, LosartanMappingMetaData, Coadministration, Genotype
)
from pkdb_models.models.losartan.experiments.studies import *


logger = get_logger(__name__)


# --- Filters ---
def filter_control(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Return control experiments/mappings."""

    metadata: LosartanMappingMetaData = fit_mapping.metadata

    # only PO and IV (no SL, MU, RE)
    if metadata.route not in {Route.PO, Route.IV}:
        return False

    # filter coadminstration
    if metadata.coadministration != Coadministration.NONE:
        return False

    # filter health (no renal, cardiac impairment, ...)
    if metadata.health not in {Health.HEALTHY, Health.HYPERTENSION}:
        return False

    # filter multiple dosing (only single dosing)
    if metadata.dosing == Dosing.MULTIPLE:
        return False

    # only fasted subjects
    # FIXME: possible influence of food, ignoring for now: see Lv2014, Shah2009
    # if metadata.fasting not in {Fasting.FASTED, Fasting.NR}:
    #     return False

    # filter genotypes
    if metadata.genotype not in {Genotype.CYP2C9_1_1, Genotype.ABCB1_GG_CC, Genotype.NR}:
        return False

    # remove outliers
    if metadata.outlier is True:
        return False

    return True


# --- Fit experiments ---
f_fitexp_kwargs = dict(
    # PK fitted with: FDA1995S60, FDA1995S67, Goldberg1995a, Munafo1992, Ohtawa1993

    experiment_classes=[
        # Azizi1999,
        # Doig1993,
        # Bae2011,
        # # Christen1991a, # FIXME: ang1 and ang2 protocols
        # Donzelli2014,
        # Fischer2002,
        # FDA1995S60,
        # FDA1995S67,
        # Goldberg1995,
        # Goldberg1995a,
        # Han2009a,
        # Huang2021,
        # Kim2016,
        # Kobayashi2008,
        # Lee2003b,
        # Li2009,
        # Lo1995,
        # Munafo1992,
        # Oh2012,
        # Ohtawa1993,
        # Puris2019,
        # Sekino2003,
        # Shin2020,
        # Sica1995,
        # Tanaka2014,
        # Yasar2002a,

        # PD data
        Azizi1999,
        Doig1993,
        Goldberg1995,
        # Goldberg1995a,
        Ohtawa1993,
        Munafo1992,
        Sekino2003,
    ],
    base_path=LOSARTAN_PATH,
    data_path=DATA_PATHS,
)


def filter_pharmacodynamics(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only pharmacodynamics data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if "ren" in yid:
        return True
    if "ang1" in yid:
        return True
    if "ang2" in yid:
        return True
    if "ald" in yid:
        return True
    if "SBP" in yid:
        return True
    if "DBP" in yid:
        return True
    if "MAP" in yid:
        return True

    return False

def filter_pharmacokinetics(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only pharmacokinetics data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if "los" in yid:
        return True
    if "e3174" in yid:
        return True

    return False


def f_fitexp_all():
    """All data."""
    return f_fitexp(metadata_filters=filter_empty, **f_fitexp_kwargs)


def f_fitexp_control() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(metadata_filters=filter_control, **f_fitexp_kwargs)

def f_fitexp_pk() -> Dict[str, List[FitExperiment]]:
    """Pharmacokinetic control data."""
    return f_fitexp(metadata_filters=[filter_control, filter_pharmacokinetics], **f_fitexp_kwargs)

def f_fitexp_pd() -> Dict[str, List[FitExperiment]]:
    """Pharmacodynamic control data."""
    return f_fitexp(metadata_filters=[filter_control, filter_pharmacodynamics], **f_fitexp_kwargs)


if __name__ == "__main__":
    """Test construction of FitExperiments."""

    for f in [
        f_fitexp_all,
        # f_fitexp_control,
        # f_fitexp_pk,
        f_fitexp_pd,
    ]:
        console.rule(style="white")
        console.print(f"{f.__name__}")
        fitexp = f()
