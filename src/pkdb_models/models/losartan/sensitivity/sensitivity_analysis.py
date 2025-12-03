"""Sensitivity analysis of the losartan model."""

import numpy as np
from roadrunner._roadrunner import NamedArray
from sbmlsim.sensitivity.global_sensitivity import (
    SensitivitySimulation, SensitivityAnalysis, LocalSensitivityAnalysis, GlobalSobolSensitivityAnalysis,
    SamplingSensitivityAnalysis
)
from sbmlsim.sensitivity.parameters import (
    parameters_for_sensitivity_analysis
)
from sbmlutils.console import console
from pkdb_analysis.pk.pharmacokinetics import TimecoursePK
from pint import UnitRegistry


class LosartanSensitivitySimulation(SensitivitySimulation):
    """Simulation for sensitivity calculation."""
    tend = 5 * 24 * 60  # [min]

    def simulate(self, changes: dict[str, float]) -> dict[str, float]:
        # apply changes and simulate
        all_changes = {
            **self.changes_simulation,
            **changes
        }
        self.apply_changes(all_changes, reset_all=True)
        s: NamedArray = self.rr.simulate(start=0, end=self.tend)
        self._plot(s)

        # pharmacokinetic parameters
        y: dict[str, float] = {}

        # pharmacokinetics
        ureg = UnitRegistry()
        Q_ = ureg.Quantity

        # losartan
        time = Q_(s["time"], "min")
        tcpk = TimecoursePK(
            time=time,
            concentration=Q_(s["[Cve_los]"], "mM"),
            substance="losartan",
            ureg=ureg,
            dose=Q_(10, "mg"),
        )
        pk_dict = tcpk.pk.to_dict()
        # console.print(pk_dict)
        for pk_key in [
            "aucinf",
            "cmax",
            "thalf",
            "kel",
        ]:
            y[f"[Cve_los]_{pk_key}"] = pk_dict[pk_key]

        # E3174, L158
        for sid in ["[Cve_e3174]", "[Cve_l158]"]:
            tcpk = TimecoursePK(
                time=time,
                concentration=Q_(s[sid], "mM"),
                substance="losartan",
                ureg=ureg,
                dose=None,
            )
            pk_dict = tcpk.pk.to_dict()
            # console.print(pk_dict)
            for pk_key in [
                "aucinf",
                "cmax",
                "thalf",
                "kel",
            ]:
                y[f"{sid}_{pk_key}"] = pk_dict[pk_key]

        # pharmacodynamics
        for sid, f in [
            ("[ang1]", "max"),
            ("[ang2]", "max"),
            ("[ren]", "max"),
            ("[ald]", "min"),
            ("SBP", "min"),
            ("DBP", "min"),
            ("MAP", "min"),
        ]:
            # minimal and maximal value of readout
            if f == "max":
                y[f"{sid}_max"] = np.min(s[sid])
            elif f == "min":
                y[f"{sid}_min"] = np.max(s[sid])

        return y

    def _plot(self, s: NamedArray) -> None:

        # plotting
        from matplotlib import pyplot as plt
        plt.plot(
            s["time"], s["[Cve_los]"],
            # df["time"], df["[Cve_los]"],
            marker="o",
            markeredgecolor="black",
            color="tab:blue",
        )
        plt.show()

if __name__ == "__main__":
    from pkdb_models.models.losartan import MODEL_PATH, MODEL_BASE_PATH

    sensitivity_simulation = LosartanSensitivitySimulation(
        model_path=MODEL_PATH,
        selections=[
            "time",
            "[Cve_los]",
            "[Cve_e3174]",
            "[Cve_l158]",
            "[ang1]",
            "[ang2]",
            "[ren]",
            "[ald]",
            "SBP",
            "DBP",
            "MAP",
        ],
        changes_simulation = {
            "PODOSE_los": 10.0,  # [mg]
        }
    )
    console.print(sensitivity_simulation.outputs)

    # FIXME: exclude zeros

    # parameters for sensitivity analysis
    parameters = parameters_for_sensitivity_analysis(
        sbml_path=MODEL_PATH,
        exclude_ids={
            "conversion_min_per_day",  # constant conversion factor
            "Mr_los",  # molecular weight
            "Mr_e3174",  # molecular weight
            "Mr_l158",  # molecular weight
            "PODOSE_los", # dose
        }
    )
    console.print(parameters)

    #
    # SamplingSensitivityAnalysis
    # sa = GlobalSobolSensitivityAnalysis(
    #     sensitivity_simulation=sensitivity_simulation,
    #     parameters=list(parameters.keys()),
    # )

    sa = LocalSensitivityAnalysis(
        sensitivity_simulation=sensitivity_simulation,
        parameters=list(parameters.keys()),
        difference=0.1,
    )
    sa.create_samples()
    sa.simulate_samples()
    sa.calculate_sensitivity()
    sa.plot_sensitivity()