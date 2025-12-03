from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models import losartan
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)
from pkdb_models.models.losartan.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, LosartanMappingMetaData, Coadministration, Genotype,
)

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.losartan.helpers import run_experiments


class Christen1991a(LosartanSimulationExperiment):
    """Simulation experiment of Christen1991a.

    FIXME: implement the angiotensin I and angiotensin II protocol with pressure response
    """

    info = {
        "[ren]": "renin",
        "[ang2]": "ang2",
        "[ald]": "aldosterone",
        "SBP": "sbp",
        "DBP": "dbp",
        "SBP_ratio": "sbp_ratio",
    }
    interventions_single = {
        "placebo": 0,
        "LOS10": 10,
        "LOS20": 20,
        "LOS40": 40
    }
    interventions_multi = {
        "placeboM": 0,
        "LOS5M": 5,
        "LOS10M": 10,
        "LOS20M": 20,
        "LOS40M": 40
    }
    colors = {
        0: "black",
        5: "tab:red",
        10: "tab:green",
        20: "tab:blue",
        40: "tab:brown"
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig4", "Fig5", "Fig6"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.endswith("_renin"):
                    dset.unit_conversion("mean", 1 / self.Mr.ren)
                # elif label.endswith("_ang2"):
                #     dset.unit_conversion("mean", 1 / self.Mr.ang2)
                elif label.endswith("_aldosterone"):
                    dset.unit_conversion("mean", 1 / self.Mr.ald)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        # single dose
        for intervention, dose in self.interventions_single.items():
            tcsims[f"po_{intervention}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=35 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(73, "kg"),
                        "PODOSE_los": Q_(dose, "mg"),
                        "ang2_ref": Q_(5.2, "fmole/ml"), #placebo
                        "[ang2]": Q_(5.2, "fmole/ml"),  #placebo
                        "ald_ref": Q_(78.4, "pg/ml") /self.Mr.ald, #placebo
                        "[ald]": Q_(78.4, "pg/ml") /self.Mr.ald, #placebo
                    },
                ),
            )
        for intervention, dose in self.interventions_multi.items():
            # multiple dose
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(70, "kg"),
                    "PODOSE_los": Q_(dose, "mg"),
                    "ang2_ref": Q_(4.06, "fmole/ml"), #placebo
                    "[ang2]": Q_(4.06, "fmole/ml"), #placebo
                    "ren_ref": Q_(5.9, "pg/ml") / self.Mr.ren,  # placebo
                    "[ren]": Q_(5.9, "pg/ml") / self.Mr.ren,  # placebo
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_los": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=32 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_los": Q_(dose, "mg"),
                },
            )

            tcsims[f"po_{intervention}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(6)] + [tc2],
                # time_offset=-6*24*60,
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for k, sid in enumerate(["SBP_ratio", "[ang2]", "[ald]"]):
            for intervention, dose in self.interventions_single.items():
                name = self.info[sid]
                mappings[f"task_po_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{intervention}_{name}",
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
                        outlier=True,  # FIXME: ang1 and ang2 infusion protocol
                    ),
                )

        for k, sid in enumerate(["SBP_ratio", "[ren]", "[ang2]"]):
            for intervention, dose in self.interventions_multi.items():
                name = self.info[sid]
                mappings[f"task_po_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{intervention}_{name}",
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
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                        outlier=True,  # FIXME: ang1 and ang2 infusion protocol
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.fig1_fig5(),
            **self.fig4_fig6(),
        }

    def fig1_fig5(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_Fig5",
            num_rows=3,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["SBP_ratio"], unit=self.units["SBP_ratio"])
        plots[1].set_yaxis(self.labels["[ang2]"], unit=self.units["[ang2]"])
        plots[2].set_yaxis(self.labels["[ald]"], unit=self.units["[ald]"])

        for k, sid in enumerate(["SBP_ratio", "[ang2]", "[ald]"]):
            for intervention, dose in self.interventions_single.items():
                name = self.info[sid]

                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {dose} mg",
                    color=self.colors[dose],
                )
                # data
                plots[k].add_data(
                    dataset=f"{intervention}_{name}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=f"{dose} mg",
                    color=self.colors[dose],
                )

        return {
            fig.sid: fig,
        }

    def fig4_fig6(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig4_Fig6",
            num_rows=3,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["SBP_ratio"], unit=self.units["SBP_ratio"])
        plots[1].set_yaxis(self.labels["[ren]"], unit=self.units["[ren]"])
        plots[2].set_yaxis(self.labels["[ang2]"], unit=self.units["[ang2]"])

        for k, sid in enumerate(["SBP_ratio", "[ren]", "[ang2]"]):
            for intervention, dose in self.interventions_multi.items():
                name = self.info[sid]
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {dose} mg",
                    color=self.colors[dose],
                )
                # data
                plots[k].add_data(
                    dataset=f"{intervention}_{name}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    color=self.colors[dose],
                    linestyle="",
                )

        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    out = losartan.RESULTS_PATH_SIMULATION / Christen1991a.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Christen1991a, output_dir=Christen1991a.__name__)
