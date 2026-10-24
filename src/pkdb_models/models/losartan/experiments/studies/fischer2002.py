"""Simulation experiment of the study Fischer2002."""

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
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)


class Fischer2002(LosartanSimulationExperiment):
    """Simulation experiment of Fischer2002.

    Losartan potassium 50 mg (tablet) once daily for 10 days, alone or with
    phenytoin 4 mg/kg/day at bedtime, in 16 healthy subjects (14 CYP2C9 *1/*1,
    2 *1/*2), the last dose after an overnight fast: losartan and E3174 in
    plasma (mean + SD, log axis in the paper) and urine (0-24 hr) after the
    10th dose. The ten doses are simulated, the time is relative to the 10th
    dose. Phenytoin inhibits CYP2C9 (`cyp2c9_interactions`), its arm is
    validation data.
    """

    interventions: ClassVar[list[str]] = ["LOS50", "LOS50_PHE4"]
    colors: ClassVar[dict[str, str]] = {
        "LOS50": "black",
        "LOS50_PHE4": "tab:blue",
    }
    intervention_labels: ClassVar[dict[str, str]] = {
        "LOS50": "50 mg",
        "LOS50_PHE4": "50 mg + phenytoin",
    }
    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "Aurine_los": "losartan_urine",
        "Aurine_e3174": "exp3174_urine",
    }

    bodyweight = 72  # [kg]
    n_doses = 10
    #: the mean CYP2C9 activity of the subjects, 14 *1/*1 and 2 *1/*2
    cyp2c9_group = (
        14 * LosartanSimulationExperiment.cyp2c9_activity["*1/*1"]
        + 2 * LosartanSimulationExperiment.cyp2c9_activity["*1/*2"]
    ) / 16

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Tab2A"]:
            df = self.load_dataframe(fig_id)
            # the predose losartan of LOS50 (0.1 ng/ml, no SD) is the bottom of
            # the log axis, i.e. below the LOQ of 1 ng/ml
            df = df[~((df.label == "losartan_LOS50") & (df.time == 0))]
            for label, df_label in df.groupby("label"):
                label = str(label)
                if "recovery" in label:
                    # the same amounts in percent of the dose
                    continue
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if fig_id in {"Fig1", "Fig2"}:
                    if label.startswith("losartan"):
                        dset.unit_conversion("mean", 1 / self.Mr.los)
                    elif label.startswith("exp3174"):
                        dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        dose = self.losartan_dose(Q(50, "mg"))
        end = (self.n_doses - 1) * 24 + 25  # [hr]
        tcsims = {}
        for intervention in self.interventions:
            f_cyp2c9 = self.cyp2c9_group
            if "PHE" in intervention:
                f_cyp2c9 *= self.cyp2c9_interactions["phenytoin"]
            # a dose every 24 hr, the urine is collected after the 10th dose
            tcsims[f"po_{intervention.lower()}"] = Simulation(
                time_unit="hr",
                end=end,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "LI__f_cyp2c9": Q(f_cyp2c9, "dimensionless"),
                    "PODOSE_los": dose,
                },
                changes=[
                    self.add_doses(
                        [24 * k for k in range(1, self.n_doses)], {"PODOSE_los": dose}
                    ),
                    Change(
                        [24 * k for k in range(1, self.n_doses)],
                        {
                            "Aurine_los": Q(0, "mmole"),  # reset
                            "Aurine_e3174": Q(0, "mmole"),  # reset
                        },
                    ),
                ],
                time_shift=-(self.n_doses - 1) * 24,
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                coadministration = Coadministration.NONE
                if "PHE" in intervention:
                    coadministration = Coadministration.PHENYTOIN

                mappings[f"fm_po_los50_{name}_{intervention}"] = FitMapping(
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
                        task=f"task_po_{intervention.lower()}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=coadministration,
                        genotype=Genotype.NR,
                        data=PKPDData.PK,
                        data_type=DataType.MEAN,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1_2(),
        }

    def fig1_2(self) -> dict[str, Figure]:
        """Figures 1 and 2 and table 2: losartan and E3174 in plasma and urine."""
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy, 10th dose)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-1, max=25),
            legend=True,
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los, scale="log", min=0.1)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174, scale="log", min=1)
        plots[2].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
        plots[3].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)

        for k, (sid, name) in enumerate(self.info.items()):
            for intervention in self.interventions:
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention.lower()}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {self.intervention_labels[intervention]}",
                    color=self.colors[intervention],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.intervention_labels[intervention],
                    color=self.colors[intervention],
                    linestyle="" if "urine" in sid else "--",
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Fischer2002, output_dir=Fischer2002.__name__)
