"""Simulation experiment of the study Azizi1999."""

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


class Azizi1999(LosartanSimulationExperiment):
    """Simulation experiment of Azizi1999.

    Single oral dose of losartan potassium 50 mg (tablet) or placebo in 16
    healthy male subjects after a light breakfast (the candesartan cilexetil
    arms of the study are not simulated): losartan and E3174 in plasma, plasma
    renin (active renin concentration), angiotensin I and II and the mean
    arterial pressure (change from the baseline, curated as 87 mmHg plus the
    change).

    The subjects were mildly sodium depleted (furosemide 40 mg at 9 pm on the
    evening before and a diet of 30 mmol sodium for 36 hr), which the model
    describes only by the baselines of the RAAS of each arm. Losartan in plasma
    is reported as the potassium salt (HPLC of Merck). The RAAS of the placebo
    arm falls at 12 hr with the meals, which the model does not describe.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "[ren]": "renin",
        "[ang1]": "ang1",
        "[ang2]": "ang2",
        "MAP": "map",
    }
    pk_sids: ClassVar[set[str]] = {"[Cveplasma_los]", "[Cveplasma_e3174]"}
    colors: ClassVar[dict[str, str]] = {
        "placebo": "black",
        "LOS50": "tab:blue",
    }
    losp_doses: ClassVar[dict[str, int]] = {
        "placebo": 0,  # [mg] losartan potassium
        "LOS50": 50,  # [mg] losartan potassium
    }
    interventions: ClassVar[list[str]] = list(losp_doses.keys())
    legend_labels: ClassVar[dict[str, str]] = {
        "placebo": "placebo",
        "LOS50": "50 mg",
    }

    #: the baselines of the arms, the values before the dose (Fig1, Fig2):
    #: renin, angiotensin I and II [pg/ml] and the mean arterial pressure
    #: [mmHg] (the change from the baseline is plotted, 87 mmHg is the baseline
    #: of the curation)
    baselines: ClassVar[dict[str, dict[str, float]]] = {
        "placebo": {"ren": 58.5, "ang1": 11.8, "ang2": 7.2, "MAP": 87.0},
        "LOS50": {"ren": 70.9, "ang1": 14.6, "ang2": 11.1, "MAP": 87.0},
    }
    #: the systolic blood pressure is not reported, an assumption which splits
    #: the mean arterial pressure into SBP and DBP (only MAP is mapped)
    sbp_ref: ClassVar[float] = 120.0  # [mmHg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l
                if label.endswith("_renin"):
                    dset.unit_conversion("mean", 1 / self.Mr.ren)
                elif label.endswith("_ang1"):
                    dset.unit_conversion("mean", 1 / self.Mr.ang1)
                elif label.endswith("_ang2"):
                    dset.unit_conversion("mean", 1 / self.Mr.ang2)
                elif label.endswith("_losartan"):
                    # reported as losartan potassium
                    dset.unit_conversion("mean", 1 / self.Mr.losp)
                elif label.endswith("_exp3174"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for intervention in self.interventions:
            baseline = self.baselines[intervention]
            ren = Q(baseline["ren"], "pg/ml") / self.Mr.ren
            ang1 = Q(baseline["ang1"], "pg/ml") / self.Mr.ang1
            ang2 = Q(baseline["ang2"], "pg/ml") / self.Mr.ang2
            # MAP = DBP + (SBP - DBP)/3
            dbp = (3 * baseline["MAP"] - self.sbp_ref) / 2
            tcsims[f"po_{intervention}"] = Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "PODOSE_los": self.losartan_dose(
                        Q(self.losp_doses[intervention], "mg")
                    ),
                    "ren_ref": ren,
                    "[ren]": ren,
                    "ang1_ref": ang1,
                    "[ang1]": ang1,
                    "ang2_ref": ang2,
                    "[ang2]": ang2,
                    "SBP_ref": Q(self.sbp_ref, "mmHg"),
                    "DBP_ref": Q(dbp, "mmHg"),
                },
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for intervention in self.interventions:
                if intervention == "placebo" and sid in self.pk_sids:
                    continue

                mappings[f"fm_po_{intervention}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{intervention}_{name}",
                        xid="time",
                        yid="mean",
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
                        tissue=Tissue.NONE if sid == "MAP" else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        # a light breakfast 2 hr before the dose
                        fasting=Fasting.LIGHT_BREAKFAST,
                        # furosemide 40 mg on the evening before the dose
                        coadministration=Coadministration.FUROSEMIDE,
                        data=PKPDData.PK if sid in self.pk_sids else PKPDData.PD,
                        data_type=DataType.MEAN,
                    ),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1_2_3(),
        }

    def fig1_2_3(self) -> dict[str, Figure]:
        """Figures 1, 2 and 3: pharmacokinetics, RAAS and blood pressure."""
        fig = Figure(
            experiment=self,
            sid="Fig1_2_3",
            num_rows=3,
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy, sodium depleted)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.labels["[ren]"], unit=self.units["[ren]"], min=0)
        plots[3].set_yaxis(self.labels["[ang1]"], unit=self.units["[ang1]"], min=0)
        plots[4].set_yaxis(self.labels["[ang2]"], unit=self.units["[ang2]"], min=0)
        plots[5].set_yaxis(self.labels["MAP"], unit=self.units["MAP"], min=60)

        for k, (sid, name) in enumerate(self.info.items()):
            for intervention in self.interventions:
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {self.legend_labels[intervention]}",
                    color=self.colors[intervention],
                )
                # data
                if intervention == "placebo" and sid in self.pk_sids:
                    continue
                plots[k].add_data(
                    dataset=f"{intervention}_{name}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=self.legend_labels[intervention],
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Azizi1999, output_dir=Azizi1999.__name__)
