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


class Azizi1999(LosartanSimulationExperiment):
    """Simulation experiment of Azizi1999."""

    info = {
        "[Cve_los]": "losartan",
        "[Cve_e3174]": "exp3174",
        "[ren]": "renin",
        "[ang1]": "ang1",
        "[ang2]": "ang2",
        "MAP": "map",
    }
    colors = {
        "placebo": "black",
        "LOS50": "tab:blue",
    }
    losp_doses = {
        "placebo": 0,  # [mg]
        "LOS50": 50,  # [mg]
    }
    interventions = list(losp_doses.keys())


    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.endswith("_renin"):
                    dset.unit_conversion("mean", 1 / self.Mr.ren)
                elif label.endswith("_ang1"):
                    dset.unit_conversion("mean", 1 / self.Mr.ang1)
                elif label.endswith("_ang2"):
                    dset.unit_conversion("mean", 1 / self.Mr.ang2)
                elif label.endswith("_losartan"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.endswith("_exp3174"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        map_value = 87  # mmHg
        for intervention in self.interventions:
            tcsims[f"po_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=25 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "PODOSE_los": Q_(self.losp_doses[intervention], "mg") * self.Mr.los/self.Mr.losp,
                        "ren_ref" : Q_(58.5, "pg/ml") / self.Mr.ren, # placebo
                        "[ren]": Q_(58.5, "pg/ml") / self.Mr.ren,  # placebo
                        "ang1_ref" : Q_(11.8, "pg/ml") / self.Mr.ang1,  # placebo
                        "[ang1]": Q_(11.8, "pg/ml") / self.Mr.ang1,  # placebo
                        "ang2_ref" : Q_(7.2, "pg/ml") / self.Mr.ang2,  # placebo
                        "[ang2]": Q_(7.2, "pg/ml") / self.Mr.ang2,  # placebo
                        # MAP 87  MAP = DBP + (SBP - DBP)/3  =>
                        "SBP_ref" : Q_(120, "mmHg"),
                        "DBP_ref" : Q_((3 * map_value - 120) / 2, "mmHg"),
                    },
                )]
            )
        # console.print(tcsims)
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention in self.interventions:
                if intervention == "placebo" and sid in {"[Cve_los]", "[Cve_e3174]"}:
                    continue

                mappings[f"fm_po_{intervention}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{intervention}_{name}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
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
                        fasting=Fasting.FED,
                        coadministration=Coadministration.NONE,
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.fig1_2_3(),
        }

    def fig1_2_3(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_2_3",
            num_rows=3,
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.labels["[ren]"], unit=self.units["[ren]"], min=0)
        plots[3].set_yaxis(self.labels["[ang1]"], unit=self.units["[ang1]"], min=0)
        plots[4].set_yaxis(self.labels["[ang2]"], unit=self.units["[ang2]"], min=0)
        plots[5].set_yaxis(self.labels["MAP"], unit=self.units["MAP"], min=60)

        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for intervention in self.interventions:
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {intervention}",
                    color=self.colors[intervention],
                )
                # data
                if intervention == "placebo" and sid in {"[Cve_los]", "[Cve_e3174]"}:
                    continue
                plots[k].add_data(
                    dataset=f"{intervention}_{name}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=f"{intervention}",
                    color=self.colors[intervention],
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = losartan.RESULTS_PATH_SIMULATION / Azizi1999.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Azizi1999, output_dir=Azizi1999.__name__)
