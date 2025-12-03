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


class Huang2021(LosartanSimulationExperiment):
    """Simulation experiment of Huang2021."""

    interventions = ["placebo", "breviscapine"]
    colors = {
        "placebo": "black",
        "breviscapine": "tab:blue",
    }
    markers = {
        "placebo": "s",
        "breviscapine": "o",
    }

    info = {
        "los": "losartan",
        "e3174": "exp3174",
    }

    genotypes = {
        "1": "*1/*1",
        "3": "*1/*3",
    }
    bodyweights = {
        "all": 54,
        "1": 54,
        "3": 52,
    }  # [kg]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.startswith("losartan_"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174_"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        # Fig2
        tcsims[f"po_los50"] = TimecourseSim(
            [Timecourse(
                start=0,
                end=40 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweights["all"], "kg"),
                    "PODOSE_los": Q_(50, "mg") * self.Mr.los / self.Mr.losp,
                },
            )]
        )

        # Fig3
        for number, genotype in self.genotypes.items():
            tcsims[f"po_los50_{number}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=40* 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweights[number], "kg"),
                        "LI__f_cyp2c9": Q_(self.cyp2c9_activity[genotype], "dimensionless"),
                        "PODOSE_los": Q_(50, "mg") * self.Mr.los/self.Mr.losp,
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
                coadministration = Coadministration.NONE
                if "breviscapine" in intervention:
                    coadministration = Coadministration.BREVISCAPINE

                # Fig2
                mappings[f"fm_po_los50_{sid}_{intervention}"] = FitMapping(
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
                        self, task=f"task_po_los50", xid="time", yid=f"[Cve_{sid}]",
                    ),
                    metadata=LosartanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=coadministration,
                        genotype=Genotype.CYP2C9_1_1_AND_1_3,
                    ),
                )

                # Fig3
                for number, genotype in self.genotypes.items():
                    mappings[f"fm_po_los50_{sid}_{intervention}_{number}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{name}_{intervention}_{number}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self, task=f"task_po_los50_{number}", xid="time", yid=f"[Cve_{sid}]",
                        ),
                        metadata=LosartanMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.FASTED,
                            coadministration=coadministration,
                            genotype=Genotype(genotype),
                        ),
                    )

        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.fig2(),
            **self.fig3(),
        }

    def fig2(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)

        for k, sid in enumerate(self.info):
            name = self.info[sid]
            # simulation
            plots[k].add_data(
                task=f"task_po_los50",
                xid="time",
                yid=f"[Cve_{sid}]",
                label=f"Sim 50 mg",
                color="black",
            )
            for intervention in self.interventions:

                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{intervention} 50 mg",
                    color=self.colors[intervention],
                    marker=self.markers[intervention],
                )

        return {
            fig.sid: fig,
        }

    def fig3(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)

        for k, sid in enumerate(self.info):
            name = self.info[sid]
            for number, genotype in self.genotypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los50_{number}",
                    xid="time",
                    yid=f"[Cve_{sid}]",
                    label=f"Sim {genotype} 50 mg",
                    color=self.cyp2c9_colors[genotype],
                )
                for intervention in self.interventions:
                    # data
                    plots[k].add_data(
                        dataset=f"{name}_{intervention}_{number}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=f"{genotype} {intervention} 50 mg",
                        color=self.cyp2c9_colors[genotype],
                        marker=self.markers[intervention],
                    )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = losartan.RESULTS_PATH_SIMULATION / Huang2021.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Huang2021, output_dir=Huang2021.__name__)
