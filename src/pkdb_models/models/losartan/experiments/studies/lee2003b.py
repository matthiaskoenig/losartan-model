"""Simulation experiment of the study Lee2003b."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Genotype,
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


class Lee2003b(LosartanSimulationExperiment):
    """Simulation experiment of Lee2003b.

    Single oral dose of losartan potassium 50 mg as tablet after an overnight
    fast in 15 healthy men and women (79 kg) with the CYP2C9 genotypes *1/*1,
    *1/*2 and *1/*3 (five each): losartan and E3174 in plasma (Fig 1, mean and
    SD). The points on the floor of the log axis (0.1 ng/ml, below the limit of
    quantification) are not data and are removed from the copy.
    """

    info: ClassVar[dict[str, str]] = {
        "los": "losartan",
        "e3174": "exp3174",
    }
    genotypes: ClassVar[dict[str, str]] = {
        "1": "*1/*1",
        "2": "*1/*2",
        "3": "*1/*3",
    }
    bodyweight = 79.0  # kg

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.startswith("losartan_"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174_"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for number, genotype in self.genotypes.items():
            simulations[f"po_los50_{number}"] = Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "LI__f_cyp2c9": Q(self.cyp2c9_activity[genotype], "dimensionless"),
                    "BW": Q(self.bodyweight, "kg"),
                    "PODOSE_los": self.losartan_dose(Q(50, "mg")),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for number, genotype in self.genotypes.items():
                mappings[f"fm_po_los50_{sid}_{number}"] = FitMapping(
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
                        self,
                        task=f"task_po_los50_{number}",
                        xid="time",
                        yid=f"[Cveplasma_{sid}]",
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                        genotype=Genotype(f"CYP2C9 {genotype}"),
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
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

        for k, (sid, name) in enumerate(self.info.items()):
            for number, genotype in self.genotypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los50_{number}",
                    xid="time",
                    yid=f"[Cveplasma_{sid}]",
                    label=f"Sim {genotype} 50 mg",
                    color=self.cyp2c9_colors[genotype],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{number}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype} 50 mg",
                    color=self.cyp2c9_colors[genotype],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Lee2003b, output_dir=Lee2003b.__name__)
