"""Reusable functionality for multiple simulation experiments."""

from typing import Dict

from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

from pkdb_models.models.losartan import MODEL_BASE_PATH


class RaasSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 13
    suptitle_font_size = 40

    units: Dict[str, str] = {
        "time": "hr",
        "[e3174]": "nM",  # e3174 in Vplasma [0-1000 nM]

        "[anggen]": "pM",  # angiotensinogen in Vplasma
        "[ang1]": "pM",  # angiotensin I in Vplasma
        "[ang2]": "pM",  # angiotensin II in Vplasma

        "[ren]": "pM",  # renin in Vplasma
        "[ald]": "pM",  # aldosterone in Vplasma

        "fa_e3174": "dimensionless",  # [-] activation via exp3174
        "fi_e3174": "dimensionless",  # [-] inhibition via exp3174

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
        # "MAP_change": "mmHg",

        "SBP_ratio": "dimensionless",
        "DBP_ratio": "dimensionless",
        # "MAP_ratio": "dimensionless",
    }
    labels: Dict[str, str] = {
        "time": "hr",
        "[e3174]": "e3174",  # e3174 in Vplasma

        "[anggen]": "angiotensinogen",  # angiotensinogen in Vplasma
        "[ang1]": "angiotensin I",  # angiotensin I in Vplasma
        "[ang2]": "angiotensin II",  # angiotensin II in Vplasma

        "[ren]": "renin",  # renin in Vplasma
        "[ald]": "aldosterone",  # aldosterone in Vplasma

        "fa_e3174": "E3174 activation",  # [-] activation via exp3174
        "fi_e3174": "E3174 inhibition",  # [-] inhibition via exp3174

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
        # "MAP_change": "mean arterial blood pressure change",

        "SBP_ratio": "systolic blood pressure ratio",
        "DBP_ratio": "diastolic blood pressure ratio",
        # "MAP_ratio": "mean arterial blood pressure ratio",
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_BASE_PATH / "losartan_raas.xml",
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {}

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return RaasSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",
                "[e3174]",

                "[anggen]",
                "[ang1]",
                "[ang2]",

                "[ren]",
                "[ald]",

                "fa_e3174",
                "fi_e3174",

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
                # "MAP_change",

                "SBP_ratio",
                "DBP_ratio",
                # "MAP_ratio",
            ]
        )
        return {}
