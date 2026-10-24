"""Simulation experiment of the study Kim2016."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    DataType,
    Dosing,
    Fasting,
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


class Kim2016(LosartanSimulationExperiment):
    """Simulation experiment of Kim2016.

    Single oral dose of losartan potassium 50 mg (Cozaar) as tablet in a
    fasting cocktail of probe drugs (caffeine, losartan, omeprazole,
    dextromethorphan, midazolam; fexofenadine on the next day) in 15 healthy
    Korean men (71.9 kg, no CYP2C9 poor metabolizers), before (control) and
    after two weeks of fermented red ginseng: losartan and E3174 in plasma
    (mean +/- SE in the figure, stored as SD). Both arms are a cocktail, the
    control arm is not a control of the fits. Ginseng weakly inhibits CYP2C9
    (`cyp2c9_interactions`, metabolic ratio x0.77), its arm is validation data.
    """

    info: ClassVar[dict[str, str]] = {
        "los": "losartan",
        "e3174": "exp3174",
    }
    interventions: ClassVar[list[str]] = ["control", "ginseng"]
    colors: ClassVar[dict[str, str]] = {
        "control": "black",
        "ginseng": "tab:blue",
    }
    bodyweight = 71.9  # [kg]

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
        for intervention in self.interventions:
            f_cyp2c9 = 1.0
            if intervention == "ginseng":
                f_cyp2c9 = self.cyp2c9_interactions["ginseng"]
            simulations[f"po_los50_{intervention}"] = Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "LI__f_cyp2c9": Q(f_cyp2c9, "dimensionless"),
                    "PODOSE_los": self.losartan_dose(Q(50, "mg")),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                coadministration = Coadministration.COCKTAIL
                if "ginseng" in intervention:
                    coadministration = Coadministration.GINSENG

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
                        task=f"task_po_los50_{intervention}",
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
                        data=PKPDData.PK,
                        data_type=DataType.MEAN,
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
            for intervention in self.interventions:
                # simulation
                plots[k].add_data(
                    task=f"task_po_los50_{intervention}",
                    xid="time",
                    yid=f"[Cveplasma_{sid}]",
                    label=f"Sim {intervention} 50 mg",
                    color=self.colors[intervention],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{intervention} 50 mg",
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kim2016, output_dir=Kim2016.__name__)
