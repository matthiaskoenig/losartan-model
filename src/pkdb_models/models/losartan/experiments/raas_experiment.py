"""Base class of the simulation experiments of the RAAS model alone."""

from collections.abc import Mapping
from pathlib import Path
from typing import ClassVar, override

from pkpdsim.experiment import PKPDSimulationExperiment
from sbmlsim.data import Data
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task
from sbmlsim.units import Quantity

from pkdb_models.models.losartan import MODEL, MODEL_BASE_PATH


class RaasSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the simulation experiments of the RAAS model alone.

    The RAAS model is the part of the whole-body model which the E3174 in
    plasma acts on: renin, angiotensinogen, angiotensin I and II, aldosterone
    and the blood pressure.
    """

    pkpd_model = MODEL

    font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 22}
    scan_font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 13
    suptitle_font_size = 40

    units: ClassVar[Mapping[str, str]] = {
        "time": "hr",
        "[e3174]": "nM",  # e3174 in Vplasma [0-1000 nM]
        "[anggen]": "pM",  # angiotensinogen in Vplasma
        "[ang1]": "pM",  # angiotensin I in Vplasma
        "[ang2]": "pM",  # angiotensin II in Vplasma
        "[ren]": "pM",  # renin in Vplasma
        "ren_ratio": "dimensionless",
        "[ald]": "pM",  # aldosterone in Vplasma
        "ald_ratio": "dimensionless",
        "fe_e3174": "dimensionless",  # [-] effect of E3174
        "ANGGEN2ANG1": "pmol/min",  # angiotensinogen to angiotensin I (renin)
        "ANG1ANG2": "pmol/min",  # angiotensin I to angiotensin II (ACE)
        "ANG2DEG": "pmol/min",  # angiotensin II degradation (ANG2DEG)
        "RENSEC": "pmol/min",  # renin secretion (RENSEC)
        "RENDEG": "pmol/min",  # renin degradation (RENDEG)
        "ALDSEC": "pmol/min",  # aldosterone secretion (ALDSEC)
        "ALDDEG": "pmol/min",  # aldosterone degradation (ALDDEG)
        "SBP": "mmHg",
        "DBP": "mmHg",
        "MAP": "mmHg",
        "SBP_change": "mmHg",
        "DBP_change": "mmHg",
        "SBP_ratio": "dimensionless",
        "DBP_ratio": "dimensionless",
    }
    labels: ClassVar[Mapping[str, str]] = {
        "time": "hr",
        "[e3174]": "e3174",  # e3174 in Vplasma
        "[anggen]": "angiotensinogen",  # angiotensinogen in Vplasma
        "[ang1]": "angiotensin I",  # angiotensin I in Vplasma
        "[ang2]": "angiotensin II",  # angiotensin II in Vplasma
        "[ren]": "renin",  # renin in Vplasma
        "ren_ratio": "renin ratio",
        "[ald]": "aldosterone",  # aldosterone in Vplasma
        "ald_ratio": "aldosterone ratio",
        "fe_e3174": "E3174 effect",  # [-] effect of E3174
        "ANGGEN2ANG1": "ANGGEN2ANG1",  # angiotensinogen to angiotensin I (renin)
        "ANG1ANG2": "ANG1ANG2",  # angiotensin I to angiotensin II (ANG1ANG2)
        "ANG2DEG": "ANG2DEG",  # angiotensin II degradation (ANG2DEG)
        "RENSEC": "renin secretion\n(RENSEC)",  # renin secretion (RENSEC)
        "RENDEG": "renin degradation\n(RENDEG)",  # renin degradation (RENDEG)
        "ALDSEC": "aldosterone secretion\n(ALDSEC)",  # aldosterone secretion (ALDSEC)
        "ALDDEG": "aldosterone degradation\n(ALDDEG)",  # aldosterone degradation (ALDDEG)
        "SBP": "systolic blood pressure",
        "DBP": "diastolic blood pressure",
        "MAP": "mean arterial blood pressure",
        "SBP_change": "systolic blood pressure change",
        "DBP_change": "diastolic blood pressure change",
        "SBP_ratio": "systolic blood pressure ratio",
        "DBP_ratio": "diastolic blood pressure ratio",
    }

    @override
    def models(self) -> dict[str, AbstractModel | Path]:
        return {
            "model": AbstractModel(
                source=MODEL_BASE_PATH / "losartan_raas.xml",
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    def default_changes(self) -> dict[str, Quantity | float]:
        """Changes of every simulation of the experiments, none at the moment.

        The parameters are part of the model and are not changed in the
        simulations: a parameter which is fitted must not be set by a simulation,
        it would override the value of the fit.
        """
        return {}

    @override
    def tasks(self) -> dict[str, Task]:
        return {
            f"task_{key}": Task(model="model", simulation=key)
            for key in self.simulations()
        }

    @override
    def data(self) -> dict[str, Data]:
        self.add_selections_data(
            selections=[
                "time",
                "[e3174]",
                "[anggen]",
                "[ang1]",
                "[ang2]",
                "[ren]",
                "ren_ratio",
                "[ald]",
                "ald_ratio",
                "fe_e3174",
                "ANGGEN2ANG1",
                "ANG1ANG2",
                "ANG2DEG",
                "RENSEC",
                "RENDEG",
                "ALDSEC",
                "ALDDEG",
                "SBP",
                "DBP",
                "MAP",
                "SBP_change",
                "DBP_change",
                "SBP_ratio",
                "DBP_ratio",
            ]
        )
        return {}
