"""Simulation experiment of the study Tanaka2014."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
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


class Tanaka2014(LosartanSimulationExperiment):
    """Simulation experiment of Tanaka2014.

    Single oral dose of 50 mg losartan (potassium, the Inje cocktail) with
    caffeine, omeprazole, dextromethorphan and midazolam in 4 healthy Japanese
    men (34-41 years, no body weight reported): losartan, E3174 and the
    E3174/losartan ratio in plasma (Fig 3B) and the amounts of losartan and
    E3174 in the 0-8 h urine (Tab 5), mean and SD.

    The ratio of the figure is the ratio of the mass concentrations, the
    observable is the molar ratio: the data are multiplied by
    Mr(losartan)/Mr(E3174). The ratio is derived from the two concentration
    curves (the same subjects twice), its points at 0.5 and 1 h are ratios of
    concentrations near the limit of quantification and are dropped.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "mr_e3174_los_plasma": "exp3174_losartan",
        "Aurine_los": "losartan_urine",
        "Aurine_e3174": "exp3174_urine",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3B", "Tab5"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label == "exp3174_losartan":
                    # ratios of concentrations near the limit of quantification
                    df_label = df_label[df_label.time > 1.0]
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l and mole
                if label.startswith("losartan"):
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label.startswith("exp3174_losartan"):
                    # the mass ratio, curated in ng/ml, as molar ratio
                    dset.unit_conversion(
                        "mean", self.Mr.los / self.Mr.e3174 / Q(1, "ng/ml")
                    )
                elif label.startswith("exp3174"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)

                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        # single dose
        tcsims["po_los50"] = Simulation(
            time_unit="min",
            end=12 * 60,
            preinit_changes={
                **self.default_changes(),
                "PODOSE_los": self.losartan_dose(Q(50, "mg")),
            },
        )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            mappings[self.to_sid(f"fm_po_los50_{name}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(name),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_los50",
                    xid="time",
                    yid=sid,
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.URINE if sid.startswith("Aurine") else Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.COCKTAIL,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig3_tab5(),
        }

    def fig3_tab5(self) -> dict[str, Figure]:
        """Losartan, E3174 and their ratio in plasma, the amounts in urine."""
        fig = Figure(
            experiment=self,
            sid="Fig3_Tab5",
            num_rows=2,
            num_cols=3,
            name=f"{self.__class__.__name__} (healthy, cocktail)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_mr_e3174_los, unit=self.unit_mr)
        plots[3].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
        plots[4].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)
        plots[5].xaxis = None  # empty panel

        for k, (sid, name) in enumerate(self.info.items()):
            # simulation
            plots[k].add_data(
                task="task_po_los50",
                xid="time",
                yid=sid,
                label="Sim 50 mg",
                color="black",
            )
            # data
            plots[k].add_data(
                dataset=self.to_sid(name),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label="50 mg",
                color="black",
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Tanaka2014, output_dir=Tanaka2014.__name__)
