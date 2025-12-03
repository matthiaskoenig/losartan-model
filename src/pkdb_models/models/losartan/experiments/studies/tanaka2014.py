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


class Tanaka2014(LosartanSimulationExperiment):
    """Simulation experiment of Tanaka2014."""

    info = {
        "[Cve_los]": "losartan",
        "[Cve_e3174]": "exp3174",
        "mr_e3174_los_plasma": "exp3174_losartan",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3B"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label == "losartan":
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label == "exp3174":
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                elif label == "exp3174_losartan":
                    # remove unit of ratio
                    dset.unit_conversion("mean", 1 / self.Q_(1, "ng/ml"))

                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        # single dose
        tcsims[f"po_los50"] = TimecourseSim(
            Timecourse(
                start=0,
                end=12 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "PODOSE_los": Q_(50, "mg") * self.Mr.los / self.Mr.losp,
                    },
                ),
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            mappings[f"task_po_{name}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{name}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_los50", xid="time", yid=sid,
                ),
                metadata=LosartanMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.COCKTAIL,
                ),
            )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.fig3(),
        }

    def fig3(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=3,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_mr_e3174_los, unit=self.unit_mr)

        for k, sid in enumerate(self.info):
            name = self.info[sid]

            # simulation
            plots[k].add_data(
                task=f"task_po_los50",
                xid="time",
                yid=sid,
                label=f"Sim 50 mg",
                color="black",
            )
            # data
            plots[k].add_data(
                dataset=f"{name}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"50 mg",
                color="black",
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = losartan.RESULTS_PATH_SIMULATION / Tanaka2014.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Tanaka2014, output_dir=Tanaka2014.__name__)
