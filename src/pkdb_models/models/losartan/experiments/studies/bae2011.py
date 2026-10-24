"""Simulation experiment of the study Bae2011."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    DataType,
    Dosing,
    Fasting,
    Genotype,
    Health,
    PKPDData,
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


class Bae2011(LosartanSimulationExperiment):
    """Simulation experiment of Bae2011.

    Single oral dose of losartan potassium 50 mg (tablet) in 13 fasted healthy
    Korean male subjects with the CYP2C9 genotypes *1/*1 (12 subjects) and
    *3/*3 (one subject): losartan and E3174 in plasma.
    """

    info: ClassVar[dict[str, str]] = {
        "los": "losartan",
        "e3174": "exp3174",
    }

    genotypes: ClassVar[dict[str, str]] = {
        "1": "*1/*1",
        "3": "*3/*3",
    }
    #: the standard deviation of the curves, none for the single *3/*3 subject
    yid_sds: ClassVar[dict[str, str | None]] = {
        "1": "mean_sd",
        "3": None,
    }

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
        tcsims = {}
        for number, genotype in self.genotypes.items():
            tcsims[f"po_los50_{number}"] = Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "LI__f_cyp2c9": Q(self.cyp2c9_activity[genotype], "dimensionless"),
                    "PODOSE_los": self.losartan_dose(Q(50, "mg")),
                },
            )
        return tcsims

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
                        yid_sd=self.yid_sds[number],
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
                        data=PKPDData.PK,
                        data_type=DataType.MEAN
                        if number == "1"
                        else DataType.INDIVIDUAL,
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
                    yid_sd=self.yid_sds[number],
                    count="count",
                    label=f"{genotype} 50 mg",
                    color=self.cyp2c9_colors[genotype],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Bae2011, output_dir=Bae2011.__name__)
