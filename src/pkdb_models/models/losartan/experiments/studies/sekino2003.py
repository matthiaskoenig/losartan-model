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


class Sekino2003(LosartanSimulationExperiment):
    """Simulation experiment of Sekino2003."""

    info = {
        "mr_e3174_los_plasma": "exp3174_losartan_plasma",
        "mr_e3174_los_urine": "exp3174_losartan_urine",
        "SBP_change": "sbp_change",
        "DBP_change": "dbp_change",
    }

    genotypes = {
        "1": "*1/*1",
        "3": "*1/*3",
    }

    bodyweights = {
        "1": 65.7,
        "3": 61.7,
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1","Tab1A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                dsets[label] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for number, genotype in self.genotypes.items():
            tcsims[f"po_los25_{number}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=15 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweights[number], "kg"),
                        "LI__f_cyp2c9": Q_(self.cyp2c9_activity[genotype], "dimensionless"),
                        "PODOSE_los": Q_(25, "mg") * self.Mr.los/self.Mr.losp,
                    },
                )]
            )
        # console.print(tcsims)
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for number, genotype in self.genotypes.items():
                mappings[f"fm_po_los25_{name}_{number}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{number}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_los25_{number}", xid="time", yid=sid,
                    ),
                    metadata=LosartanMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                        genotype=Genotype(genotype),
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.fig1(),
            **self.tab1a(),
        }

    def tab1a(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Tab1A",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_mr_los_e3174, unit=self.unit_mr)
        plots[1].set_yaxis(self.label_mr_los_e3174_urine, unit=self.unit_mr)

        for k, sid in enumerate(["mr_e3174_los_plasma", "mr_e3174_los_urine"]):
            name = self.info[sid]

            for number, genotype in self.genotypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los25_{number}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{number}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )

        return {
            fig.sid: fig,
        }

    def fig1(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["SBP_change"], unit=self.units["SBP_change"])
        plots[1].set_yaxis(self.labels["DBP_change"], unit=self.units["DBP_change"])


        for k, sid in enumerate(["SBP_change", "DBP_change"]):
            name = self.info[sid]

            for number, genotype in self.genotypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los25_{number}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{number}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = losartan.RESULTS_PATH_SIMULATION / Sekino2003.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Sekino2003, output_dir=Sekino2003.__name__)
