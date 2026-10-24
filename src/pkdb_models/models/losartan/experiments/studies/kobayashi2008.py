"""Simulation experiment of the study Kobayashi2008."""

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


class Kobayashi2008(LosartanSimulationExperiment):
    """Simulation experiment of Kobayashi2008.

    Single oral dose of losartan potassium 25 mg as tablet after an overnight
    fast in six healthy Japanese men (62.3 kg), alone and after bucolome 300 mg
    once daily for 7 days: losartan and E3174 in serum (samples 0-24 hr).
    Bucolome is a potent CYP2C9 inhibitor (`cyp2c9_interactions`, losartan AUC
    x1.67, E3174 AUC x0.165), its arm is validation data.
    """

    bodyweight = 62.3  # [kg]
    info: ClassVar[dict[str, str]] = {
        "los": "losartan",
        "e3174": "exp3174",
    }
    interventions: ClassVar[list[str]] = ["LOS", "LOS_BUC"]
    #: the legend names of the arms
    legends: ClassVar[dict[str, str]] = {
        "LOS": "losartan 25 mg",
        "LOS_BUC": "losartan 25 mg + bucolome 300 mg",
    }
    colors: ClassVar[dict[str, str]] = {
        "LOS": "black",
        "LOS_BUC": "tab:blue",
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
        simulations = {}
        for intervention in self.interventions:
            f_cyp2c9 = 1.0
            if "BUC" in intervention:
                f_cyp2c9 = self.cyp2c9_interactions["bucolome"]
            simulations[f"po_los25_{intervention.lower()}"] = Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "LI__f_cyp2c9": Q(f_cyp2c9, "dimensionless"),
                    "PODOSE_los": self.losartan_dose(Q(25, "mg")),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                coadministration = Coadministration.NONE
                if "BUC" in intervention:
                    coadministration = Coadministration.BUCOLOME

                mappings[f"fm_po_los25_{sid}_{intervention}"] = FitMapping(
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
                        task=f"task_po_los25_{intervention.lower()}",
                        xid="time",
                        yid=f"[Cveplasma_{sid}]",
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.SERUM,
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
                    task=f"task_po_los25_{intervention.lower()}",
                    xid="time",
                    yid=f"[Cveplasma_{sid}]",
                    label=f"Sim {self.legends[intervention]}",
                    color=self.colors[intervention],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.legends[intervention],
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kobayashi2008, output_dir=Kobayashi2008.__name__)
