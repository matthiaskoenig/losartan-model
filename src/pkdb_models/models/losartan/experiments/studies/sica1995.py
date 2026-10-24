"""Simulation experiment of the study Sica1995."""

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
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)


class Sica1995(LosartanSimulationExperiment):
    """Simulation experiment of Sica1995.

    Oral doses of 100 mg losartan (potassium) once daily for 7 days in
    hypertensive patients with normal renal function (group I, creatinine
    clearance >= 75 ml/min), mild to moderate (group II, 30-74 ml/min) and
    severe renal insufficiency (group III, 10-29 ml/min), six per group:
    losartan (Fig 1) and E3174 (Fig 2) in plasma and their urinary excretion
    0-24 h (Tab 2) after the last dose.

    The patients continued an ongoing diuretic therapy (Study design), the
    dose was given at 08:00 without a statement on food. The plasma
    concentrations are expressed as losartan potassium and E3174 free acid
    (Sample analysis). The renal function is the creatinine clearance over
    101 ml/min, the GFR of the model. The digitized curves have no error
    bars; the first three E3174 points of group III (0-20 min) are a rough
    guess of the curation and not used. The AUC0-24 (Fig 3) and the renal
    clearances (Fig 4) are summaries of the same curves and urines and not
    used.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "Aurine_los": "losartan_urine",
        "Aurine_e3174": "exp3174_urine",
    }
    groups: ClassVar[dict[str, str]] = {
        "group1": "Normal renal function",
        "group2": "Mild renal impairment",
        "group3": "Severe renal impairment",
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "group1": 84.6,
        "group2": 75.7,
        "group3": 75.4,
    }
    #: creatinine clearance [ml/min] (Table 1) relative to the GFR of 101 ml/min
    renal_functions: ClassVar[dict[str, float]] = {
        "group1": 112.2 / 101,
        "group2": 58.8 / 101,
        "group3": 16.6 / 101,
    }
    health: ClassVar[dict[str, Health]] = {
        "group1": Health.HYPERTENSION,
        "group2": Health.RENAL_IMPAIRMENT,  # and hypertension
        "group3": Health.RENAL_IMPAIRMENT,  # and hypertension
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Tab2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                # points of the curation which are a rough guess
                df_label = df_label[df_label.comments != "very rough guess"]
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l, losartan as potassium salt
                if fig_id == "Fig1" and label.startswith("losartan_"):
                    dset.unit_conversion("mean", 1 / self.Mr.losp)
                elif fig_id == "Fig2" and label.startswith("exp3174_"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for group in self.groups:
            dose = self.losartan_dose(Q(100, "mg"))
            # a dose every 24 hr for 7 days, the urine is collected from dose
            # to dose, the time is relative to the last dose
            tcsims[f"po_los100_{group}"] = Simulation(
                time_unit="hr",
                end=6 * 24 + 65,
                preinit_changes={
                    **self.default_changes(),
                    "KI__f_renal_function": Q(
                        self.renal_functions[group], "dimensionless"
                    ),
                    "BW": Q(self.bodyweights[group], "kg"),
                    "PODOSE_los": dose,
                },
                changes=[
                    self.add_doses([24 * k for k in range(1, 7)], {"PODOSE_los": dose}),
                    Change(
                        [24 * k for k in range(1, 7)],
                        {
                            "Aurine_los": Q(0, "mmole"),  # reset
                            "Aurine_e3174": Q(0, "mmole"),  # reset
                        },
                    ),
                ],
                time_shift=-6 * 24,
            )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            urine = sid.startswith("Aurine")
            for group in self.groups:
                mappings[f"fm_po_los100_{name}_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        # the digitized plasma curves have no error bars
                        yid_sd="mean_sd" if urine else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_los100_{group}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if urine else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=self.health[group],
                        fasting=Fasting.NR,
                        # ongoing diuretic therapy
                        coadministration=Coadministration.COMEDICATION,
                        genotype=Genotype.NR,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1_2(),
            **self.tab2(),
        }

    def fig1_2(self) -> dict[str, Figure]:
        """Losartan and E3174 in plasma after the last dose."""
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=2,
            name=f"{self.__class__.__name__} (renal impairment)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)

        for k, sid in enumerate(["[Cveplasma_los]", "[Cveplasma_e3174]"]):
            name = self.info[sid]
            for group, group_name in self.groups.items():
                color = self.renal_colors[group_name]

                # simulation
                plots[k].add_data(
                    task=f"task_po_los100_{group}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {group}",
                    color=color,
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=f"{group}",
                    color=color,
                )

        return self.figures_by_sid(fig)

    def tab2(self) -> dict[str, Figure]:
        """Urinary excretion of losartan and E3174 after the last dose."""
        fig = Figure(
            experiment=self,
            sid="Tab2",
            num_rows=2,
            name=f"{self.__class__.__name__} (renal impairment)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24), legend=True
        )
        plots[0].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
        plots[1].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)

        for k, sid in enumerate(["Aurine_los", "Aurine_e3174"]):
            name = self.info[sid]
            for group, group_name in self.groups.items():
                color = self.renal_colors[group_name]

                # simulation
                plots[k].add_data(
                    task=f"task_po_los100_{group}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {group}",
                    color=color,
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{group}",
                    color=color,
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Sica1995, output_dir=Sica1995.__name__)
