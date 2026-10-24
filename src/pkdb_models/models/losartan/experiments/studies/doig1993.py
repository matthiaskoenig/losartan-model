"""Simulation experiment of the study Doig1993."""

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
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)


class Doig1993(LosartanSimulationExperiment):
    """Simulation experiment of Doig1993.

    Dose ranging study in 12 healthy male subjects (76.5 kg), sodium depleted
    before each study day (40 mmol sodium per day and furosemide 40 mg twice
    daily for 3 days), fasted and supine: single oral doses of losartan
    potassium (DuP 753/MK-954) 5, 10, 25, 50 and 100 mg or placebo with 200 ml
    of clear fluid (the formulation is not reported), meals at 2, 5 and 9 hr.
    Plasma renin activity (PRA) and aldosterone (Fig2) and the supine systolic
    and diastolic blood pressure of placebo, 50 and 100 mg (Fig1). Mean and SD
    (the SD only for placebo and 100 mg).

    The PRA (ng angiotensin I/ml/hr) is the renin activity, not the renin
    concentration: it is mapped relative to the value of its arm before the
    dose (`ren_ratio`); the renin rows of the curation (PRA x 7.6) are not used.
    The sodium depletion is described only by the baselines of each arm. The
    placebo arm changes with the meals (blood pressure -3 mmHg, PRA rises at
    5-10 hr), which the model does not describe.
    """

    info: ClassVar[dict[str, str]] = {
        "pra": "ren_ratio",
        "aldosterone": "[ald]",
        "supine_sbp": "SBP",
        "supine_dbp": "DBP",
    }

    colors: ClassVar[dict[str, str]] = {
        "placebo": "black",
        "LOS5": "tab:brown",
        "LOS10": "tab:orange",
        "LOS25": "tab:red",
        "LOS50": "tab:green",
        "LOS100": "tab:blue",
    }

    losp_doses: ClassVar[dict[str, int]] = {
        "placebo": 0,  # [mg] losartan potassium
        "LOS5": 5,  # [mg] losartan potassium
        "LOS10": 10,  # [mg] losartan potassium
        "LOS25": 25,  # [mg] losartan potassium
        "LOS50": 50,  # [mg] losartan potassium
        "LOS100": 100,  # [mg] losartan potassium
    }
    #: the arms with the blood pressure (Fig1)
    bp_interventions: ClassVar[list[str]] = ["placebo", "LOS50", "LOS100"]

    #: the baselines of the arms before the dose: aldosterone [pg/ml] at 0 hr
    #: (Fig2) and the supine SBP and DBP [mmHg], the mean of -0.5 and 0 hr
    #: (Fig1); the arms without blood pressure take the values of placebo
    baselines: ClassVar[dict[str, dict[str, float]]] = {
        "placebo": {"ald": 774.1, "SBP": 112.4, "DBP": 62.9},
        "LOS5": {"ald": 771.9, "SBP": 112.4, "DBP": 62.9},
        "LOS10": {"ald": 874.0, "SBP": 112.4, "DBP": 62.9},
        "LOS25": {"ald": 699.8, "SBP": 112.4, "DBP": 62.9},
        "LOS50": {"ald": 636.5, "SBP": 113.0, "DBP": 65.4},
        "LOS100": {"ald": 716.4, "SBP": 112.1, "DBP": 61.9},
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label.endswith("_renin"):
                    # PRA x 7.6, not a measurement
                    continue
                if label.endswith("_pra"):
                    # relative to the value of the arm before the dose
                    df_ratio = df_label.copy()
                    baseline = df_ratio.loc[df_ratio.time == 0, "mean"].to_numpy()[0]
                    for key in ["mean", "sd", "se"]:
                        df_ratio[key] = df_ratio[key] / baseline
                    df_ratio["unit"] = "dimensionless"
                    dsets[label] = DataSet.from_df(df_ratio, self.ureg)
                    continue

                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.endswith(("_aldosterone", "_aldosterone_change_abs")):
                    dset.unit_conversion("mean", 1 / self.Mr.ald)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        for intervention, dose in self.losp_doses.items():
            baseline = self.baselines[intervention]
            ald = Q(baseline["ald"], "pg/ml") / self.Mr.ald
            # the baseline for 2 hr, the dose at 0 hr
            tcsims[f"po_{intervention}"] = Simulation(
                time_unit="min",
                end=27 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(76.5, "kg"),
                    "ald_ref": ald,
                    "[ald]": ald,
                    "SBP_ref": Q(baseline["SBP"], "mmHg"),
                    "DBP_ref": Q(baseline["DBP"], "mmHg"),
                },
                changes=[
                    Change(
                        2 * 60,
                        {"PODOSE_los": self.losartan_dose(Q(dose, "mg"))},
                    )
                ],
                time_shift=-2 * 60,
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for name, sid in self.info.items():
            for intervention in self.losp_doses:
                if sid in {"SBP", "DBP"} and intervention not in self.bp_interventions:
                    continue
                mappings[f"fm_po_{intervention}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{intervention}_{name}",
                        xid="time",
                        yid="mean",
                        # the SD only for placebo and 100 mg: all unweighted
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.NONE if sid in {"SBP", "DBP"} else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.NR,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        # furosemide 40 mg twice daily on the 3 days before
                        coadministration=Coadministration.FUROSEMIDE,
                        data=PKPDData.PD,
                        data_type=DataType.MEAN,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1_2(),
            **self.fig2_aldosterone(),
        }

    def fig1_2(self) -> dict[str, Figure]:
        """Figures 1 and 2: PRA, aldosterone and blood pressure."""
        fig = Figure(
            experiment=self,
            sid="Fig1_Fig2",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy, sodium depleted)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis("PRA ratio\n", unit=self.units["ren_ratio"])
        plots[1].set_yaxis(self.labels["[ald]"], unit=self.units["[ald]"])
        plots[2].set_yaxis(self.labels["SBP"], unit=self.units["SBP"])
        plots[3].set_yaxis(self.labels["DBP"], unit=self.units["DBP"])

        for k, (name, sid) in enumerate(self.info.items()):
            for intervention, dose in self.losp_doses.items():
                if sid in {"SBP", "DBP"} and intervention not in self.bp_interventions:
                    continue
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {dose} mg",
                    color=self.colors[intervention],
                )
                # data
                plots[k].add_data(
                    dataset=f"{intervention}_{name}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)

    def fig2_aldosterone(self) -> dict[str, Figure]:
        """Figure 2: aldosterone, its absolute and its relative change."""
        fig = Figure(
            experiment=self,
            sid="Fig2_aldosterone",
            num_rows=3,
            num_cols=1,
            name=f"{self.__class__.__name__} (healthy, sodium depleted)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.labels["[ald]"], unit=self.units["[ald]"])
        plots[1].set_yaxis(self.labels["ald_change"], unit=self.units["ald_change"])
        plots[2].set_yaxis(self.labels["ald_ratio"], unit=self.units["ald_ratio"])

        datasets = {
            "[ald]": "aldosterone",
            "ald_change": "aldosterone_change_abs",
            "ald_ratio": "aldosterone_change_rel",
        }
        for intervention, dose in self.losp_doses.items():
            for k, (sid, name) in enumerate(datasets.items()):
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {dose} mg",
                    color=self.colors[intervention],
                )
                # data
                plots[k].add_data(
                    dataset=f"{intervention}_{name}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{dose} mg",
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Doig1993, output_dir=Doig1993.__name__)
