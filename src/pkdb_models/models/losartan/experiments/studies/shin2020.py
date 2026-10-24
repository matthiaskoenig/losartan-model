"""Simulation experiment of the study Shin2020."""

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


class Shin2020(LosartanSimulationExperiment):
    """Simulation experiment of Shin2020.

    Single oral dose of 50 mg losartan potassium after an overnight fast in 38
    healthy Korean subjects grouped by the ABCB1 c.2677G>T/c.3435C>T diplotype
    (GG/CC n=13, GT/CT n=12, TT/TT n=13): losartan, E3174 and their sum in
    plasma (Fig 1) and the urinary excretion of losartan plus E3174 (Tab 4).

    Table 4 gives the excretion in the intervals 0-4 and 4-8 h, the copy of
    the data holds the cumulative amounts (8 h: the sum of both intervals) for
    the cumulative observable `Aurine_los_e3174`.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "[Cveplasma_los_e3174]": "losartan_exp3174",
        "Aurine_los_e3174": "losartan_exp3174_urine",
    }

    genotypes: ClassVar[dict[str, str]] = {
        "GG_CC": "GG/CC",
        "GT_CT": "GT/CT",
        "TT_TT": "TT/TT",
    }

    #: the genotypes of the metadata: the diplotypes of c.2677G>T and c.3435C>T
    metadata_genotypes: ClassVar[dict[str, Genotype]] = {
        "GG_CC": Genotype.ABCB1_2677GG_3435CC,
        "GT_CT": Genotype.ABCB1_2677GT_3435CT,
        "TT_TT": Genotype.ABCB1_2677TT_3435TT,
    }

    #: body weight of the groups [kg]
    bodyweights: ClassVar[dict[str, float]] = {
        "GG_CC": 68.5,
        "GT_CT": 66.5,
        "TT_TT": 67.1,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab4"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for number, genotype in self.genotypes.items():
            tcsims[f"po_los50_{number}"] = Simulation(
                time_unit="min",
                end=12 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[number], "kg"),
                    "GU__f_abcb1": Q(self.abcb1_activity[genotype], "dimensionless"),
                    "PODOSE_los": self.losartan_dose(Q(50, "mg")),
                },
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for number, genotype in self.metadata_genotypes.items():
                mappings[self.to_sid(f"fm_po_los50_{name}_{number}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{number}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_los50_{number}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE
                        if sid == "Aurine_los_e3174"
                        else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                        genotype=genotype,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1_tab4(),
        }

    def fig1_tab4(self) -> dict[str, Figure]:
        """Losartan and E3174 in plasma and urine by ABCB1 diplotype."""
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

        for k, (sid, name) in enumerate(self.info.items()):
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
                    dataset=self.to_sid(f"{name}_{number}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype} 50 mg",
                    color=self.abcb1_colors[genotype],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Shin2020, output_dir=Shin2020.__name__)
