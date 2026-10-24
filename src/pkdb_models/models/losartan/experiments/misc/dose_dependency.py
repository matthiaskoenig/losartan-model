"""Dose dependency of the pharmacokinetics and pharmacodynamics of losartan."""

from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)


class DoseDependencyExperiment(LosartanSimulationExperiment):
    """Single doses of 0 to 200 mg losartan (oral, intravenous) and E3174 over 24 hr."""

    # the routes by substance
    routes: ClassVar[dict[str, list[str]]] = {
        "los": ["PO", "IV"],
        "e3174": ["IV"],
    }
    doses: ClassVar[list[int]] = [0, 25, 50, 100, 200]  # [mg]
    colors: ClassVar[list[str]] = [
        "black",
        "tab:orange",
        "tab:blue",
        "tab:red",
        "tab:green",
    ]

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        for substance, routes in self.routes.items():
            for route in routes:
                for dose in self.doses:
                    tcsims[f"los_{substance}_{route}_{dose}"] = Simulation(
                        time_unit="min",
                        end=24 * 60,
                        preinit_changes={
                            **self.default_changes(),
                            f"{route}DOSE_{substance}": Q(dose, "mg"),
                        },
                    )

        return tcsims

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> dict[str, Figure]:
        """Figures of the pharmacokinetics in plasma, urine and feces."""
        figures = []
        for substance, routes in self.routes.items():
            for route in routes:
                fig = Figure(
                    experiment=self,
                    sid=f"Fig_dose_dependency_pk_{substance}_{route}",
                    num_rows=3,
                    num_cols=4,
                    name=f"Dose dependency PK ({substance.upper()} {route})",
                    # the two-line y labels keep off the panel to their left
                    settings={"legend_scope": "row", "fig_subplots_wspace": 0.4},
                )
                plots = fig.create_plots(
                    xaxis=Axis(self.label_time, unit="hr"), legend=True
                )
                sids = [
                    # plasma
                    "[Cveplasma_los]",
                    "[Cveplasma_e3174]",
                    "[Cveplasma_l158]",
                    "mr_e3174_los_plasma",
                    # urine,
                    "Aurine_los",
                    "Aurine_e3174",
                    "Aurine_l158",
                    "mr_e3174_los_urine",
                    # feces
                    "Afeces_los",
                    "Afeces_e3174",
                    "Afeces_l158",
                    "mr_e3174_los_feces",
                ]
                for ksid, sid in enumerate(sids):
                    plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for ksid, sid in enumerate(sids):
                    for kval, dose in enumerate(self.doses):
                        # no ratio of the metabolite without a dose
                        if dose == 0 and sid.startswith("mr_"):
                            continue
                        plots[ksid].add_data(
                            task=f"task_los_{substance}_{route}_{dose}",
                            xid="time",
                            yid=sid,
                            label=f"{substance.upper()} {route} {dose}",
                            color=self.colors[kval],
                        )

                figures.append(fig)

        return self.figures_by_sid(*figures)

    def figure_pd(self) -> dict[str, Figure]:
        """Figures of the pharmacodynamics, the RAAS and the blood pressure."""
        figures = []
        for substance, routes in self.routes.items():
            for route in routes:
                fig = Figure(
                    experiment=self,
                    sid=f"Fig_dose_dependency_{substance}_{route}",
                    num_rows=6,
                    num_cols=3,
                    name=f"Dose dependency ({substance.upper()} {route})",
                    settings={"legend_scope": "row"},
                )
                plots = fig.create_plots(
                    xaxis=Axis(self.label_time, unit="hr"), legend=True
                )
                sids = [
                    "[e3174]",  # e3174 in Vplasma [0-1000 nM]
                    "fe_e3174",  # [-] effect via exp3174
                    "ren_ratio",  # with ald_ratio, see below
                    "[anggen]",  # angiotensinogen in Vplasma
                    "[ang1]",  # angiotensin I in Vplasma
                    "[ang2]",  # angiotensin II in Vplasma
                    "ANGGEN2ANG1",  # angiotensinogen to angiotensin I (renin)
                    "ANG1ANG2",  # angiotensin I to angiotensin II (ACE)
                    "ANG2DEG",  # angiotensin II degradation (ANG2DEG)
                    "[ren]",  # renin in Vplasma
                    "RENSEC",  # renin secretion (RENSEC)
                    "RENDEG",  # renin degradation (RENDEG)
                    "[ald]",  # aldosterone in Vplasma
                    "ALDSEC",  # aldosterone secretion (ALDSEC)
                    "ALDDEG",  # aldosterone degradation (ALDDEG)
                    "SBP",
                    "DBP",
                    "MAP",
                ]
                for ksid, sid in enumerate(sids):
                    plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for ksid, sid in enumerate(sids):
                    for kval, dose in enumerate(self.doses):
                        plots[ksid].add_data(
                            task=f"task_los_{substance}_{route}_{dose}",
                            xid="time",
                            yid=sid,
                            label=f"{substance.upper()} {route} {dose}",
                            color=self.colors[kval],
                        )

                # the renin secretion is activated, the aldosterone secretion
                # inhibited: their ratios to the baseline share a panel
                kratio = sids.index("ren_ratio")
                plots[kratio].set_yaxis(
                    label="renin|aldosterone ratio", unit=self.units["ren_ratio"]
                )
                for kval, dose in enumerate(self.doses):
                    plots[kratio].add_data(
                        task=f"task_los_{substance}_{route}_{dose}",
                        xid="time",
                        yid="ald_ratio",
                        label=None,
                        color=self.colors[kval],
                    )

                figures.append(fig)

        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(
        DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__
    )
