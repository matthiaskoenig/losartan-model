"""Simulation experiment of the study Munafo1992."""

from typing import ClassVar, override

import numpy as np
from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)


class Munafo1992(LosartanSimulationExperiment):
    """Simulation experiment of Munafo1992.

    Single oral doses of placebo and 40, 80 and 120 mg losartan potassium
    dissolved in water in six young healthy male volunteers (66.5 kg,
    crossover at intervals of one week, fasting not reported): losartan and
    E3174 in plasma (Fig 1) and the change of plasma aldosterone from the
    baseline (Fig 6, mean and SEM for placebo and 120 mg, the copy stores the
    SD). Aldosterone falls after placebo as well (circadian, about -25 pg/ml),
    which the model does not describe: the losartan arms are mapped as the
    change corrected by the placebo (the difference of the means, the SD of
    the difference of independent means where both have one), the placebo arm
    is mapped as it is.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "ald_change": "aldosterone",
    }
    interventions: ClassVar[dict[str, int]] = {
        "Placebo": 0,
        "LOS40": 40,
        "LOS80": 80,
        "LOS120": 120,
    }
    colors: ClassVar[dict[int, str]] = {
        0: "black",
        40: "tab:blue",
        80: "tab:red",
        120: "tab:green",
    }
    bodyweight = 66.5  # kg

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig6"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.startswith("losartan_"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174_"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                elif label.startswith("aldosterone_"):
                    dset.unit_conversion("mean", 1 / self.Mr.ald)
                dsets[label] = dset

        # aldosterone change corrected by the placebo
        df = self.load_dataframe("Fig6")
        placebo = df[df.label == "aldosterone_Placebo"].set_index("time")
        for intervention in self.interventions:
            if intervention == "Placebo":
                continue
            df_los = df[df.label == f"aldosterone_{intervention}"].copy()
            pla = placebo.loc[df_los.time]
            df_los["mean"] = df_los["mean"].to_numpy() - pla["mean"].to_numpy()
            df_los["sd"] = np.sqrt(
                df_los["sd"].to_numpy() ** 2 + pla["sd"].to_numpy() ** 2
            )
            df_los["se"] = np.nan
            df_los["mean_pm"] = np.nan
            df_los["label"] = f"aldosterone_corrected_{intervention}"
            dset = DataSet.from_df(df_los, self.ureg)
            dset.unit_conversion("mean", 1 / self.Mr.ald)
            dsets[f"aldosterone_corrected_{intervention}"] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        for intervention, dose in self.interventions.items():
            tcsims[f"po_{intervention}"] = Simulation(
                time_unit="min",
                end=36 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "PODOSE_los": self.losartan_dose(Q(dose, "mg")),
                },
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}

        for sid, name in self.info.items():
            for intervention in self.interventions:
                if name in {"losartan", "exp3174"} and intervention == "Placebo":
                    # no placebo data
                    continue
                dataset = f"{name}_{intervention}"
                if name == "aldosterone" and intervention != "Placebo":
                    dataset = f"aldosterone_corrected_{intervention}"

                mappings[f"fm_po_{dataset}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=dataset,
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1_fig2(),
            **self.fig6(),
        }

    def fig1_fig2(self) -> dict[str, Figure]:
        """Losartan and E3174 in plasma after the single doses."""
        fig = Figure(
            experiment=self,
            sid="Fig1_Fig2",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(
            self.labels["[Cveplasma_los]"], unit=self.units["[Cveplasma_los]"]
        )
        plots[1].set_yaxis(
            self.labels["[Cveplasma_e3174]"], unit=self.units["[Cveplasma_e3174]"]
        )

        for k, (sid, name) in enumerate(self.info.items()):
            if k == 2:
                continue

            for intervention, dose in self.interventions.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {dose} mg",
                    color=self.colors[dose],
                )
                # data
                if name in {"losartan", "exp3174"} and intervention == "Placebo":
                    # no placebo data
                    continue
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    color=self.colors[dose],
                )

        return self.figures_by_sid(fig)

    def fig6(self) -> dict[str, Figure]:
        """Change of plasma aldosterone, as measured and corrected by placebo."""
        fig = Figure(
            experiment=self,
            sid="Fig6",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["ald_change"], unit=self.units["ald_change"])
        plots[1].set_yaxis(
            "aldosterone change\ncorrected by placebo", unit=self.units["ald_change"]
        )

        for intervention, dose in self.interventions.items():
            for k in [0, 1]:
                if k == 1 and intervention == "Placebo":
                    continue
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid="ald_change",
                    label=f"Sim {dose} mg",
                    color=self.colors[dose],
                )
                # data
                dataset = (
                    f"aldosterone_{intervention}"
                    if k == 0
                    else f"aldosterone_corrected_{intervention}"
                )
                plots[k].add_data(
                    dataset=dataset,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    color=self.colors[dose],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Munafo1992, output_dir=Munafo1992.__name__)
