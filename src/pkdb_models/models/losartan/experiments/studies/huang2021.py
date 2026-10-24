"""Simulation experiment of the study Huang2021."""

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


class Huang2021(LosartanSimulationExperiment):
    """Simulation experiment of Huang2021.

    Single oral dose of losartan potassium 50 mg as tablet after an overnight
    fast in eleven healthy Chinese subjects (6 CYP2C9 *1/*1, 5 *1/*3), after
    placebo or breviscapine 40 mg three times daily for 14 days: losartan and
    E3174 in plasma of all subjects (Fig2) and by genotype (Fig3). The assay
    is calibrated with losartan potassium, losartan is converted with the
    salt. The group of all subjects is simulated with the mean CYP2C9 activity
    of its subjects. Breviscapine has no effect on the pharmacokinetics
    (all ratios not significant): its data are compared to the simulations
    without breviscapine.
    """

    interventions: ClassVar[list[str]] = ["placebo", "breviscapine"]
    colors: ClassVar[dict[str, str]] = {
        "placebo": "black",
        "breviscapine": "tab:blue",
    }
    markers: ClassVar[dict[str, str]] = {
        "placebo": "s",
        "breviscapine": "o",
    }

    info: ClassVar[dict[str, str]] = {
        "los": "losartan",
        "e3174": "exp3174",
    }

    genotypes: ClassVar[dict[str, str]] = {
        "1": "*1/*1",
        "3": "*1/*3",
    }
    bodyweights: ClassVar[dict[str, int]] = {
        "all": 54,
        "1": 54,
        "3": 52,
    }  # [kg]
    #: the mean CYP2C9 activity of all subjects, 6 *1/*1 and 5 *1/*3
    cyp2c9_all = (
        6 * LosartanSimulationExperiment.cyp2c9_activity["*1/*1"]
        + 5 * LosartanSimulationExperiment.cyp2c9_activity["*1/*3"]
    ) / 11

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.startswith("losartan_"):
                    dset.unit_conversion("mean", 1 / self.Mr.losp)
                elif label.startswith("exp3174_"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        # Fig2
        simulations["po_los50"] = Simulation(
            time_unit="min",
            end=40 * 60,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweights["all"], "kg"),
                "LI__f_cyp2c9": Q(self.cyp2c9_all, "dimensionless"),
                "PODOSE_los": self.losartan_dose(Q(50, "mg")),
            },
        )

        # Fig3
        for number, genotype in self.genotypes.items():
            simulations[f"po_los50_{number}"] = Simulation(
                time_unit="min",
                end=40 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[number], "kg"),
                    "LI__f_cyp2c9": Q(self.cyp2c9_activity[genotype], "dimensionless"),
                    "PODOSE_los": self.losartan_dose(Q(50, "mg")),
                },
            )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
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
                        self,
                        task="task_po_los50",
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
                        coadministration=coadministration,
                        genotype=Genotype.CYP2C9_1_1_AND_1_3,
                        data=PKPDData.PK,
                        data_type=DataType.MEAN,
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
                            coadministration=coadministration,
                            genotype=Genotype(f"CYP2C9 {genotype}"),
                            data=PKPDData.PK,
                            data_type=DataType.MEAN,
                        ),
                    )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig2(),
            **self.fig3(),
        }

    def fig2(self) -> dict[str, Figure]:
        """Figure of losartan and E3174 in plasma of all subjects (Fig2)."""
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

        for k, (sid, name) in enumerate(self.info.items()):
            # simulation
            plots[k].add_data(
                task="task_po_los50",
                xid="time",
                yid=f"[Cveplasma_{sid}]",
                label="Sim 50 mg",
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

        return self.figures_by_sid(fig)

    def fig3(self) -> dict[str, Figure]:
        """Figure of losartan and E3174 in plasma by genotype (Fig3)."""
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

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Huang2021, output_dir=Huang2021.__name__)
