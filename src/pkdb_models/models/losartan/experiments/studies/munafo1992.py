from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models import losartan
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)
from pkdb_models.models.losartan.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, LosartanMappingMetaData, Coadministration, Genotype,
)

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.losartan.helpers import run_experiments


class Munafo1992(LosartanSimulationExperiment):
    """Simulation experiment of Munafo1992."""

    info = {
        "[Cve_los]": "losartan",
        "[Cve_e3174]": "exp3174",
        # "SBP_ratio": "SBP_ratio",
        "ald_change": "aldosterone",
    }
    interventions = {
        "Placebo": 0,
        "LOS40": 40,
        "LOS80": 80,
        "LOS120": 120
    }
    colors = {0: "black", 40: "tab:blue", 80: "tab:red", 120: "tab:green"}
    bodyweight = 66.5  # kg

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig6"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.startswith("losartan_"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174_"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                elif label.startswith("aldosterone_"):
                    dset.unit_conversion("mean", 1 / self.Mr.ald)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for intervention, dose in self.interventions.items():
            tcsims[f"po_{intervention}"] = TimecourseSim(
            Timecourse(
                start=0,
                end=36 * 60,  # [min]
                steps=1000,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "PODOSE_los": Q_(dose, "mg"),
                    },
                ),
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention, dose in self.interventions.items():

                if name in {"losartan", "exp3174"} and intervention == "Placebo":
                    # no placebo data
                    continue


                mappings[f"task_po_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_{intervention}", xid="time", yid=sid,
                    ),
                    metadata=LosartanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.fig1_fig2(),
            **self.fig6(),
        }

    def fig1_fig2(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_Fig2",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["[Cve_los]"], unit=self.units["[Cve_los]"])
        plots[1].set_yaxis(self.labels["[Cve_e3174]"], unit=self.units["[Cve_e3174]"])

        for k, sid in enumerate(self.info):
            if k == 2:
                continue

            name = self.info[sid]
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

        return {
            fig.sid: fig,
        }


    def fig6(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig6",
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["ald_change"], unit=self.units["ald_change"])

        for k, sid in enumerate(self.info):
            if k < 2:
                continue
            name = self.info[sid]

            for intervention, dose in self.interventions.items():
                # simulation
                plots[0].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {dose} mg",
                    color=self.colors[dose],
                )
                # data
                plots[0].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    color=self.colors[dose],
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = losartan.RESULTS_PATH_SIMULATION / Munafo1992.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Munafo1992, output_dir=Munafo1992.__name__)
