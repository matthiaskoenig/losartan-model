"""Simulation experiment of the study Oh2012."""

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


class Oh2012(LosartanSimulationExperiment):
    """Simulation experiment of Oh2012.

    Single oral microdose of 2 mg losartan (potassium) in a capsule with the
    cocktail of caffeine, omeprazole, dextromethorphan and midazolam in 13
    healthy volunteers (body weight not reported, fasting not reported):
    losartan and E3174 in plasma (Fig 3, the error type is not stated). The
    E3174/losartan ratio is not reported by the study: the curation derived it
    from the two mean curves (mass ratio, converted here to the molar ratio),
    it counts the same subjects again.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "mr_e3174_los_plasma": "exp3174_losartan",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label not in self.info.values():
                    # the other drugs of the cocktail
                    continue
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label == "losartan":
                    dset.unit_conversion("mean", 1 / self.Mr.los)
                elif label == "exp3174":
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                elif label == "exp3174_losartan":
                    # the ratio of the mass concentrations is curated in ng/ml,
                    # the molar ratio is mass ratio * Mr(losartan)/Mr(E3174)
                    dset.unit_conversion(
                        "mean", self.Mr.los / self.Mr.e3174 / Q(1, "ng/ml")
                    )

                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        # single dose
        tcsims["po_los2"] = Simulation(
            time_unit="min",
            end=15 * 60,
            preinit_changes={
                **self.default_changes(),
                "PODOSE_los": Q(2, "mg") * self.Mr.los / self.Mr.losp,
            },
        )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            mappings[f"fm_po_los2_{name}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=name,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_po_los2",
                    xid="time",
                    yid=sid,
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.CAPSULE,
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
            **self.fig3(),
        }

    def fig3(self) -> dict[str, Figure]:
        """Losartan, E3174 and their ratio in plasma."""
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=3,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_mr_e3174_los, unit=self.unit_mr)

        for k, (sid, name) in enumerate(self.info.items()):
            # simulation
            plots[k].add_data(
                task="task_po_los2",
                xid="time",
                yid=sid,
                label="Sim 2 mg",
                color="black",
            )
            # data
            plots[k].add_data(
                dataset=name,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label="2 mg",
                color="black",
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Oh2012, output_dir=Oh2012.__name__)
