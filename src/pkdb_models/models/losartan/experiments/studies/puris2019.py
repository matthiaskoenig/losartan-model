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


class Puris2019(LosartanSimulationExperiment):
    """Simulation experiment of Puris2019.

    FIXME: still healthy if obese?
    """

    info = {
        "[Cve_los]": "losartan",
        "[Cve_e3174]": "exp3174",
    }

    groups = [
        "obese",
        "healthy",
        "one_year_after"
    ]

    colors = {
        "healthy": "black",
        "obese": "tab:blue",
        "one_year_after": "tab:green",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.startswith("losartan"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[f"{label}"] = dset

        # console.print(dsets)
        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        tcsims[f"po_los12.5"] = TimecourseSim(
            Timecourse(
                start=0,
                end=15 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    # "BW": Q_(self.bodyweight, "kg"), ' FIXME: no bodyweight only obese bmi
                    "PODOSE_los": Q_(12.5, "mg") * self.Mr.los / self.Mr.losp,
                    },
                ),
        )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for group in self.groups:
                mappings[f"task_po_{name}_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_los12.5", xid="time", yid=sid,
                    ),
                    metadata=LosartanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,  # FIXME: Check
                        fasting=Fasting.NR,
                        coadministration=Coadministration.COCKTAIL,
                    ),
                )
        # console.print(mappings)
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
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
                task=f"task_po_los12.5",
                xid="time",
                yid=sid,
                label=f"Sim 12.5 mg",
                color="black",
            )
            for group in self.groups:
                # data
                plots[k].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{group}",
                    color=self.colors[group],
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = losartan.RESULTS_PATH_SIMULATION / Puris2019.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Puris2019, output_dir=Puris2019.__name__)
