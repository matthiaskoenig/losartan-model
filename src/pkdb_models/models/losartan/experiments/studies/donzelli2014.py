"""Simulation experiment of the study Donzelli2014."""

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


class Donzelli2014(LosartanSimulationExperiment):
    """Simulation experiment of Donzelli2014.

    Single oral dose of losartan potassium 12.5 mg (Cosaar tablet) in 16 fasted
    healthy male subjects as part of the Basel phenotyping cocktail, given
    together with omeprazole 10 mg and midazolam 2 mg (A) or with the whole
    cocktail (E, in addition caffeine 100 mg, metoprolol 12.5 mg and efavirenz
    50 mg): losartan and E3174 in plasma (LC-MS/MS). The cocktail does not
    change the losartan AUC (ratio E/A 0.99, 0.90-1.10), so both arms are one
    simulation. The error bars are curated as SD (not stated in the paper).

    The CYP2C9 genotypes of the group are 10 *1/*1, 2 *1/*2 and 4 *1/*3,
    simulated with the mean activity of the subjects (`f_cyp2c9`).
    """

    #: the CYP2C9 genotypes of the subjects: the number of subjects
    cyp2c9_genotypes: ClassVar[dict[str, int]] = {
        "*1/*1": 10,
        "*1/*2": 2,
        "*1/*3": 4,
    }

    interventions: ClassVar[list[str]] = ["los_ome_mid", "basel_cocktail"]
    colors: ClassVar[dict[str, str]] = {
        "los_ome_mid": "tab:blue",
        "basel_cocktail": "tab:orange",
    }
    legend_labels: ClassVar[dict[str, str]] = {
        "los_ome_mid": "12.5 mg, omeprazole, midazolam",
        "basel_cocktail": "12.5 mg, Basel cocktail",
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
                if label.startswith("losartan_"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174_"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        f_cyp2c9 = sum(
            self.cyp2c9_activity[genotype] * n
            for genotype, n in self.cyp2c9_genotypes.items()
        ) / sum(self.cyp2c9_genotypes.values())
        tcsims = {}
        tcsims["po_los12p5"] = Simulation(
            time_unit="min",
            end=25 * 60,
            preinit_changes={
                **self.default_changes(),
                "PODOSE_los": self.losartan_dose(Q(12.5, "mg")),
                "LI__f_cyp2c9": Q(f_cyp2c9, "dimensionless"),
            },
        )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                mappings[self.to_sid(f"fm_po_los12.5_{sid}_{intervention}")] = (
                    FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=self.to_sid(f"{name}_{intervention}"),
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task="task_po_los12p5",
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
                            coadministration=Coadministration.COCKTAIL,
                            data=PKPDData.PK,
                            data_type=DataType.MEAN,
                        ),
                    )
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
            # simulation
            plots[k].add_data(
                task="task_po_los12p5",
                xid="time",
                yid=f"[Cveplasma_{sid}]",
                label="Sim 12.5 mg",
                color="black",
            )

            for intervention in self.interventions:
                # data
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.legend_labels[intervention],
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Donzelli2014, output_dir=Donzelli2014.__name__)
