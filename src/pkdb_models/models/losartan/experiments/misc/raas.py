"""The RAAS model alone with constant E3174 concentrations."""

from typing import ClassVar, override

import numpy as np
from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.raas_experiment import (
    RaasSimulationExperiment,
)


class RaasExperiment(RaasSimulationExperiment):
    """Effect of E3174 on the RAAS and the blood pressure.

    The RAAS model alone at its reference state for 5 hr, then a constant
    plasma concentration of E3174 (0 to 500 nM) for 24 hr and no E3174 for
    another 24 hr.
    """

    e3174_values: ClassVar[np.ndarray] = np.linspace(0.0, 500.0, num=6)  # [nM]
    # without E3174 black, then the blues of ColorBrewer's Blues (7 classes)
    # from the third on, so that the lowest concentration is visible on white
    colors: ClassVar[list[str]] = [
        "black",
        "#9ecae1",
        "#6baed6",
        "#4292c6",
        "#2171b5",
        "#084594",
    ]

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        for k, e3174 in enumerate(self.e3174_values):
            # E3174 from 0 hr to 24 hr after 5 hr without
            tcsims[f"raas_{k}"] = Simulation(
                time_unit="min",
                end=(5 + 24 + 24) * 60,
                preinit_changes={
                    **self.default_changes(),
                    "[e3174]": Q(0, "nM"),
                },
                changes=[
                    Change(5 * 60, {"[e3174]": Q(e3174, "nM")}),
                    Change(29 * 60, {"[e3174]": Q(0, "nM")}),
                ],
                time_shift=-5 * 60,
            )

        return tcsims

    @override
    def figures(self) -> dict[str, Figure]:

        info = [
            ("[e3174]", 0),  # e3174 in Vplasma [0-1000 nM]
            ("fe_e3174", 1),  # [-] effect of E3174
            # the renin secretion is activated, the aldosterone secretion inhibited
            ("ren_ratio", 2),
            ("ald_ratio", 2),
            ("[anggen]", 3),  # angiotensinogen in Vplasma
            ("[ang1]", 4),  # angiotensin I in Vplasma
            ("[ang2]", 5),  # angiotensin II in Vplasma
            ("ANGGEN2ANG1", 6),  # angiotensinogen to angiotensin I (renin)
            ("ANG1ANG2", 7),  # angiotensin I to angiotensin II (ACE)
            ("ANG2DEG", 8),  # angiotensin II degradation (ANG2DEG)
            ("[ren]", 9),  # renin in Vplasma
            ("RENSEC", 10),  # renin secretion (RENSEC)
            ("RENDEG", 11),  # renin degradation (RENDEG)
            ("[ald]", 12),  # aldosterone in Vplasma
            ("ALDSEC", 13),  # aldosterone secretion (ALDSEC)
            ("ALDDEG", 14),  # aldosterone degradation (ALDDEG)
            ("SBP", 15),
            ("MAP", 15),
            ("DBP", 15),
            ("SBP_change", 16),
            ("DBP_change", 16),
            ("SBP_ratio", 17),
            ("DBP_ratio", 17),
        ]

        fig_sid = "Fig_RAAS"
        fig = Figure(
            experiment=self,
            sid=fig_sid,
            num_rows=6,
            num_cols=3,
            name="RAAS",
        )
        # every panel draws the same concentrations of E3174, which the legend of
        # the first panel names; a legend in every panel would cover the curves
        plots = fig.create_plots(xaxis=Axis("Time", unit="hr"), legend=False)
        plots[0].legend = True

        # the axis labels of the panels with several variables
        panel_labels = {
            2: "renin|aldosterone ratio",
            15: "SBP|MAP|DBP",
            16: "SBP|DBP change",
            17: "SBP|DBP ratio",
        }
        for sid, ksid in info:
            legend = ksid == 0
            label = panel_labels.get(ksid, self.labels[sid])
            plots[ksid].set_yaxis(label=label, unit=self.units[sid])
            yaxis = plots[ksid].yaxis
            if sid == "fe_e3174" and yaxis is not None:  # set by `set_yaxis`
                yaxis.min = -0.02
                yaxis.max = 1.02
            for k, e3174 in enumerate(self.e3174_values):
                plots[ksid].add_data(
                    task=f"task_raas_{k}",
                    xid="time",
                    yid=sid,
                    label=f"{e3174:.1f} nM" if legend else None,
                    color=self.colors[k],
                )

        return {fig_sid: fig}


if __name__ == "__main__":
    MODEL.run_experiments(RaasExperiment, output_dir=RaasExperiment.__name__)
