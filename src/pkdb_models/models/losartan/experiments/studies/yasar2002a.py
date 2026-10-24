"""Simulation experiment of the study Yasar2002a."""

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


class Yasar2002a(LosartanSimulationExperiment):
    """Simulation experiment of Yasar2002a.

    Single oral dose of 50 mg losartan potassium (Cozaar) after an overnight
    fast (a standardized lunch 4 h after the dose) in 22 healthy Swedish
    subjects with the CYP2C9 genotypes *1/*1 (n=6), *1/*2 (3), *1/*3 (5),
    *2/*2 (3), *2/*3 (4) and *3/*3 (1): losartan and E3174 in plasma (Fig 1,
    *1/*1, *2/*2, *3/*3), their 24 h urinary recovery and the urinary
    losartan/E3174 ratio 0-8 and 0-24 h (Tab 2). Single oral dose of 25 mg
    losartan potassium after 3-4 h of fasting in 17 healthy Spanish subjects
    (*1/*1 n=8, *1/*2 2, *1/*3 5, *3/*3 2): the losartan/E3174 ratio in the
    0-8 h urine (Tab 3).

    The ratios are mapped as reported, losartan/E3174 (`mr_los_e3174_urine`),
    the reciprocal ratios of the copy are not used; the ratios at time 0 are
    0/0 and dropped. The urinary amounts are the recoveries of the molar dose
    of losartan potassium. The pooled ratios of Table 4 count the same
    subjects again and are not used. Body weights are reported as ranges only
    (52-91 kg, 45-90 kg), the simulations have the reference weight.
    """

    selections = (
        *LosartanSimulationExperiment.selections,
        "mr_los_e3174_urine",
    )

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "e3174",
        "Aurine_los": "losartan_urine",
        "Aurine_e3174": "e3174_urine",
        "mr_los_e3174_urine": "losartan_exp3174_urine",
    }

    genotypes: ClassVar[dict[str, str]] = {
        "1_1": "*1/*1",
        "1_2": "*1/*2",
        "1_3": "*1/*3",
        "2_2": "*2/*2",
        "2_3": "*2/*3",
        "3_3": "*3/*3",
    }
    genotypes_fig1: ClassVar[dict[str, str]] = {
        "1_1": "*1/*1",
        "2_2": "*2/*2",
        "3_3": "*3/*3",
    }

    genotypes_tab3: ClassVar[dict[str, str]] = {
        "1_1": "*1/*1",
        "1_2": "*1/*2",
        "1_3": "*1/*3",
        "3_3": "*3/*3",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Tab2", "Tab2A", "Tab3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label.startswith("exp3174_losartan"):
                    # the reciprocal of the reported ratio
                    continue
                if label.startswith("losartan_exp3174"):
                    # the ratio at time 0 is 0/0
                    df_label = df_label[df_label.time > 0]
                dsets[self.to_sid(label)] = DataSet.from_df(df_label, self.ureg)

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        # Swedish 50 mg, Spanish 25 mg
        for dose, genotypes in [(50, self.genotypes), (25, self.genotypes_tab3)]:
            for number, genotype in genotypes.items():
                tcsims[f"po_los{dose}_{number}"] = Simulation(
                    time_unit="min",
                    end=25 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "LI__f_cyp2c9": Q(
                            self.cyp2c9_activity[genotype], "dimensionless"
                        ),
                        "PODOSE_los": self.losartan_dose(Q(dose, "mg")),
                    },
                )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        # Swedish, 50 mg losartan
        for sid, name in self.info.items():
            gtypes = self.genotypes if "urine" in sid else self.genotypes_fig1

            for number, genotype in gtypes.items():
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
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                        genotype=Genotype(f"CYP2C9 {genotype}"),
                    ),
                )

        # Spanish
        sid = "mr_los_e3174_urine"
        name = "losartan_exp3174_urine"
        for number, genotype in self.genotypes_tab3.items():
            mappings[self.to_sid(f"fm_po_los25_{name}_{number}_spanish")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"{name}_{number}_spanish"),
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
                    tissue=Tissue.URINE,
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
            **self.fig1_tab2(),
            **self.tab3(),
        }

    def fig1_tab2(self) -> dict[str, Figure]:
        """Swedish subjects: losartan and E3174 in plasma and urine."""
        fig = Figure(
            experiment=self,
            sid="Fig1_Tab2",
            num_rows=3,
            num_cols=2,
            name=f"{self.__class__.__name__} (Swedish, healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
        plots[3].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)
        # the ratio of the cumulative amounts starts at 0/0, log scale
        plots[4].set_yaxis(
            self.label_mr_los_e3174_urine,
            unit=self.unit_mr,
            scale="log",
            min=0.1,
            max=200,
        )
        plots[5].xaxis = None  # empty panel

        for k, (sid, name) in enumerate(self.info.items()):
            if "urine" in sid:
                gtypes = self.genotypes
                linestyle = "dashdot"
            else:
                gtypes = self.genotypes_fig1
                linestyle = "--"

            for number, genotype in gtypes.items():
                # simulation
                plots[k].add_data(
                    task=f"task_po_los50_{number}",
                    xid="time",
                    yid=sid,
                    label=None if "urine" in sid else f"Sim {genotype}",
                    color=self.cyp2c9_colors[genotype],
                )
                # data
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{number}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{genotype}",
                    color=self.cyp2c9_colors[genotype],
                    linestyle=linestyle,
                )

        return self.figures_by_sid(fig)

    def tab3(self) -> dict[str, Figure]:
        """Spanish subjects: urinary E3174/losartan ratio."""
        fig = Figure(
            experiment=self,
            sid="Tab3",
            num_rows=1,
            name=f"{self.__class__.__name__} (Spanish, healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        # the ratio of the cumulative amounts starts at 0/0, log scale
        plots[0].set_yaxis(
            self.label_mr_los_e3174_urine,
            unit=self.unit_mr,
            scale="log",
            min=0.1,
            max=200,
        )

        sid = "mr_los_e3174_urine"
        name = "losartan_exp3174_urine"

        for number, genotype in self.genotypes_tab3.items():
            # simulation
            plots[0].add_data(
                task=f"task_po_los25_{number}",
                xid="time",
                yid=sid,
                label=None,
                color=self.cyp2c9_colors[genotype],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"{name}_{number}_spanish"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{genotype} 25 mg",
                color=self.cyp2c9_colors[genotype],
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Yasar2002a, output_dir=Yasar2002a.__name__)
