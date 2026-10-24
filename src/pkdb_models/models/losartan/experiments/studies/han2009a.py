"""Simulation experiment of the study Han2009a."""

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


class Han2009a(LosartanSimulationExperiment):
    """Simulation experiment of Han2009a.

    Single oral dose of losartan potassium 50 mg as tablet in twelve healthy
    male subjects with the CYP2C9 genotypes *1/*1 and *1/*3, after placebo or
    silymarin 140 mg three times daily for 14 days: losartan and E3174 in
    plasma (mean +/- SE in the figures, stored as SD). Silymarin (silybin, a
    mechanism-based CYP2C9 inhibitor) lowers the activity of CYP2C9
    (`cyp2c9_interactions`, metabolic ratio x0.45 in *1/*1; x0.76 in *1/*3),
    applied on top of the activity of the genotype; its arms are validation
    data.
    """

    genotypes: ClassVar[dict[str, str]] = {
        "1": "*1/*1",
        "3": "*1/*3",
    }
    interventions: ClassVar[list[str]] = ["", "_SIL140"]
    colors: ClassVar[dict[str, str]] = {
        "1": LosartanSimulationExperiment.cyp2c9_colors["*1/*1"],
        "1_SIL140": "tab:blue",
        "3": LosartanSimulationExperiment.cyp2c9_colors["*1/*3"],
        "3_SIL140": "tab:blue",
    }

    info: ClassVar[dict[str, str]] = {
        "los": "losartan",
        "e3174": "exp3174",
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
                if label.startswith("losartan"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for number, genotype in self.genotypes.items():
            for intervention in self.interventions:
                f_cyp2c9 = self.cyp2c9_activity[genotype]
                if "SIL" in intervention:
                    f_cyp2c9 *= self.cyp2c9_interactions["silymarin"]
                simulations[f"po_los50_{number}{intervention.lower()}"] = Simulation(
                    time_unit="min",
                    end=25 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "LI__f_cyp2c9": Q(f_cyp2c9, "dimensionless"),
                        "PODOSE_los": self.losartan_dose(Q(50, "mg")),
                    },
                )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for number, genotype in self.genotypes.items():
                for intervention in self.interventions:
                    mappings[f"fm_po_los50_{name}_{number}{intervention}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{name}_{number}{intervention}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=f"task_po_los50_{number}{intervention.lower()}",
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
                            coadministration=Coadministration.SILYMARIN
                            if "SIL" in intervention
                            else Coadministration.NONE,
                            genotype=Genotype(f"CYP2C9 {genotype}"),
                            data=PKPDData.PK,
                            data_type=DataType.MEAN,
                        ),
                    )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1(),
        }

    def fig1(self) -> dict[str, Figure]:
        """Figure of losartan and E3174 in plasma by genotype (Fig1)."""
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )

        for k, (sid, name) in enumerate(self.info.items()):
            for col, (number, genotype) in enumerate(self.genotypes.items()):
                plot = plots[2 * k + col]
                plot.set_yaxis(
                    self.label_los if sid == "los" else self.label_e3174,
                    unit=self.unit_los if sid == "los" else self.unit_e3174,
                )
                if k == 0:
                    plot.set_title(f"CYP2C9 {genotype}")
                for intervention in self.interventions:
                    label = genotype + (
                        " placebo" if not intervention else " silymarin"
                    )
                    # simulation
                    plot.add_data(
                        task=f"task_po_los50_{number}{intervention.lower()}",
                        xid="time",
                        yid=f"[Cveplasma_{sid}]",
                        label=f"Sim {label}",
                        color=self.colors[f"{number}{intervention}"],
                    )
                    # data
                    plot.add_data(
                        dataset=f"{name}_{number}{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=label,
                        color=self.colors[f"{number}{intervention}"],
                    )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Han2009a, output_dir=Han2009a.__name__)
