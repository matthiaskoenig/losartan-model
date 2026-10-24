"""Simulation experiment of the study Sekino2003."""

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


class Sekino2003(LosartanSimulationExperiment):
    """Simulation experiment of Sekino2003.

    Single oral dose of 25 mg losartan (potassium) after an overnight fast in
    seven healthy normotensive Japanese subjects with the CYP2C9 genotypes
    *1/*1 (n=4) and *1/*3 (n=3): change of the blood pressure (Fig 1) and the
    E3174/losartan ratio in plasma at 4 and 6 h and in the urine collected
    0-4, 0-8 and 0-12 h (Tab 1).

    The urinary ratio is the cumulative E3174 over the cumulative losartan
    excreted in the interval (Methods), the observable `mr_e3174_los_urine` is
    the ratio of the amounts excreted since the dose: the ratios 0-4, 0-8 and
    0-12 h are mapped, the ratios of the intervals 4-8 and 8-12 h are not (the
    copy keeps them as `exp3174_losartan_urine_interval_*`). The ratios at
    time 0 (0/0) carry no information and are dropped. The paper has no plasma
    concentration curves; losartan AUC0-12 149 (*1/*1) and 233 ng*h/ml
    (*1/*3), E3174 AUC0-12 1314 and 1035 ng*h/ml are given in the text only.
    """

    info: ClassVar[dict[str, str]] = {
        "mr_e3174_los_plasma": "exp3174_losartan_plasma",
        "mr_e3174_los_urine": "exp3174_losartan_urine",
        "SBP_change": "sbp_change",
        "DBP_change": "dbp_change",
    }

    genotypes: ClassVar[dict[str, str]] = {
        "1": "*1/*1",
        "3": "*1/*3",
    }

    bodyweights: ClassVar[dict[str, float]] = {
        "1": 65.7,
        "3": 61.7,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab1A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if "interval" in label:
                    # urinary ratios of single intervals, not an observable
                    continue
                if label.startswith("exp3174_losartan"):
                    # the ratio at time 0 is 0/0
                    df_label = df_label[df_label.time > 0]
                dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for number, genotype in self.genotypes.items():
            tcsims[f"po_los25_{number}"] = Simulation(
                time_unit="min",
                end=15 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[number], "kg"),
                    "LI__f_cyp2c9": Q(self.cyp2c9_activity[genotype], "dimensionless"),
                    "PODOSE_los": self.losartan_dose(Q(25, "mg")),
                },
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for number, genotype in self.genotypes.items():
                if name.startswith("exp3174_losartan"):
                    tissue = Tissue.URINE if "urine" in name else Tissue.PLASMA
                else:
                    tissue = Tissue.NONE  # blood pressure
                mappings[self.to_sid(f"fm_po_los25_{name}_{number}")] = FitMapping(
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
                        task=f"task_po_los25_{number}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=tissue,
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
        return {
            **self.fig1(),
            **self.tab1a(),
        }

    def tab1a(self) -> dict[str, Figure]:
        """E3174/losartan ratio in plasma and urine (cumulative)."""
        fig = Figure(
            experiment=self,
            sid="Tab1A",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_mr_e3174_los, unit=self.unit_mr)
        plots[1].set_yaxis(self.label_mr_e3174_los_urine, unit=self.unit_mr)

        for k, sid in enumerate(["mr_e3174_los_plasma", "mr_e3174_los_urine"]):
            name = self.info[sid]

            for number, genotype in self.genotypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los25_{number}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )
                # data
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{number}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )

        return self.figures_by_sid(fig)

    def fig1(self) -> dict[str, Figure]:
        """Change of the systolic and diastolic blood pressure."""
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["SBP_change"], unit=self.units["SBP_change"])
        plots[1].set_yaxis(self.labels["DBP_change"], unit=self.units["DBP_change"])

        for k, sid in enumerate(["SBP_change", "DBP_change"]):
            name = self.info[sid]

            for number, genotype in self.genotypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los25_{number}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )
                # data
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{number}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype} 25 mg",
                    color=self.cyp2c9_colors[genotype],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Sekino2003, output_dir=Sekino2003.__name__)
