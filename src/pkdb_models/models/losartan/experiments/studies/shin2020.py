from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

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


class Shin2020(LosartanSimulationExperiment):
    """Simulation experiment of Shin2020."""

    info = {
        "[Cve_los]": "losartan",
        "[Cve_e3174]": "exp3174",
        "[Cve_los_e3174]": "losartan_exp3174",
        "Aurine_los_e3174": "losartan_exp3174_urine"
    }

    genotypes = {
        "GG_CC": "GG/CC",
        "GT_CT": "GT/CT",
        "TT_TT": "TT/TT"
    }

    bodyweight = 67.4  # kg

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab4"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for number, genotype in self.genotypes.items():
            tcsims[f"po_los50_{number}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=20 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "GU__f_abcb1": Q_(self.abcb1_activity[genotype], "dimensionless"),
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
            for number, genotype in self.genotypes.items():
                mappings[f"fm_po_los50_{name}_{number}"] = FitMapping(
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
                        self, task=f"task_po_los50_{number}", xid="time", yid=sid,
                    ),
                    metadata=LosartanMappingMetaData(
                        tissue=Tissue.URINE if sid == "Aurine_los_e3174" else Tissue.PLASMA,
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
            **self.fig1_tab4(),
        }

    def fig1_tab4(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_Tab4",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (ABCB1)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_los_e3174, unit=self.unit_los_e3174)
        plots[3].set_yaxis(self.label_los_e3174_urine, unit=self.unit_los_e3174_urine)

        for k, sid in enumerate(self.info):
            name = self.info[sid]

            for number, genotype in self.genotypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los50_{number}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {genotype} 50 mg",
                    color=self.abcb1_colors[genotype],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{number}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype} 50 mg",
                    color=self.abcb1_colors[genotype],
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Shin2020, output_dir=Shin2020.__name__)
