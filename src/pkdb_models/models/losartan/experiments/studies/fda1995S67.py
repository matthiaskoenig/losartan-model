"""Simulation experiment of the study FDA1995S67."""

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
from sbmlsim.units import Quantity

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)


class FDA1995S67(LosartanSimulationExperiment):
    """Simulation experiment of FDA1995S67.

    Study 067 of the FDA review of NDA 020386: single doses of losartan
    potassium 50 mg oral (tablet) and 10 mg as a 20 min intravenous infusion,
    and of E3174 10 mg as a 20 min intravenous infusion in 9 fasted patients
    with alcoholic liver disease (82.3 kg, Child-Pugh scores 5 to 9: 3 class A,
    6 class B). Losartan and E3174 in plasma (Fig1) and the fraction of the
    dose in urine after 48 hr (Tab1), simulated with cirrhosis and without for
    comparison.

    The intravenous doses are the actual doses (FDA review p. 647): 9.2 mg
    losartan potassium and 8.9 mg E3174 (free acid); the plasma concentrations
    are not normalized to the nominal doses (E3174 at 20 min 2093 ng/ml in Fig1,
    2406 ng/ml normalized to 10 mg in Tab1). Losartan is reported as the
    potassium salt (HPLC of Merck). The amounts in urine are the fractions of
    the dose times the actual dose in moles (the curated amounts use the
    nominal doses and the molecular weight of losartan).
    """

    #: the actual doses [mg]: losartan potassium, E3174 free acid
    doses: ClassVar[dict[str, float]] = {"LOS50PO": 50.0, "LOS10IV": 9.2, "E10IV": 8.9}
    routes: ClassVar[dict[str, str]] = {
        "LOS50PO": "PO",
        "LOS10IV": "IV",
        "E10IV": "IV",
    }
    substances: ClassVar[dict[str, str]] = {
        "LOS50PO": "los",
        "LOS10IV": "los",
        "E10IV": "e3174",
    }
    #: the duration of the intravenous infusions
    infusion_time: ClassVar[float] = 20.0  # [min]
    colors: ClassVar[dict[str, str]] = {
        "LOS50PO": "tab:blue",
        "LOS10IV": "tab:orange",
        "E10IV": "tab:green",
    }
    #: the number of patients by Child-Pugh class (scores 5, 5, 5 and 7, 7, 7,
    #: 7, 8, 9, TabIndividuals)
    child_pugh: ClassVar[dict[str, int]] = {
        "Mild cirrhosis": 3,
        "Moderate cirrhosis": 6,
    }

    interventions: ClassVar[list[str]] = list(doses.keys())

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "Aurine_los": "losartan_urine",
        "[Cveplasma_e3174]": "exp3174",
        "Aurine_e3174": "exp3174_urine",
    }

    def f_cirrhosis(self) -> float:
        """Get the degree of cirrhosis: the mean of the patients' classes."""
        return sum(
            self.cirrhosis_map[key] * n for key, n in self.child_pugh.items()
        ) / sum(self.child_pugh.values())

    def dose_amount(self, intervention: str) -> Quantity:
        """Get the actual dose of an intervention in moles."""
        mr = self.Mr.e3174 if self.substances[intervention] == "e3174" else self.Mr.losp
        return (Q(self.doses[intervention], "mg") / mr).to("µmole")

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}

        # plasma: losartan as the potassium salt, E3174 as the free acid
        df = self.load_dataframe("Fig1")
        for label, df_label in df.groupby("label"):
            label = str(label)
            dset = DataSet.from_df(df_label, self.ureg)
            if label.startswith("losartan_"):
                dset.unit_conversion("mean", 1 / self.Mr.losp)
            elif label.startswith("exp3174_"):
                dset.unit_conversion("mean", 1 / self.Mr.e3174)
            dsets[label] = dset

        # urine: the fraction of the dose times the actual dose in moles
        df = self.load_dataframe("Tab1B")
        for label, df_label in df.groupby("label"):
            label = str(label)
            if "_recovery_" not in label:
                # the curated amounts, see the class
                continue
            intervention = label.split("_")[-1]
            amount = self.dose_amount(intervention)
            df_amount = df_label.copy()
            for key in ["mean", "sd", "se"]:
                df_amount[key] = df_amount[key] / 100 * amount.magnitude
            df_amount["unit"] = "µmole"
            dsets[label.replace("_recovery", "")] = DataSet.from_df(
                df_amount, self.ureg
            )

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        f_cirrhosis = self.f_cirrhosis()
        for intervention in self.interventions:
            route = self.routes[intervention]
            substance = self.substances[intervention]
            dose = Q(self.doses[intervention], "mg")
            if substance == "los":
                dose = self.losartan_dose(dose)

            changes = []
            if route == "PO":
                dosing = {"PODOSE_los": dose}
            else:
                rate = (dose / Q(self.infusion_time, "min")).to("mg/min")
                dosing = {f"Ri_{substance}": rate}
                changes = [
                    Change(self.infusion_time, {f"Ri_{substance}": Q(0.0, "mg/min")})
                ]

            for state, f in [("cirrhosis", f_cirrhosis), ("normal", 0.0)]:
                tcsims[f"{intervention}_{state}"] = Simulation(
                    time_unit="min",
                    end=50 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "BW": Q(82.3, "kg"),
                        "f_cirrhosis": Q(f, "dimensionless"),
                        **dosing,
                    },
                    changes=changes,
                )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            route = self.routes[intervention]
            for k, (sid, name) in enumerate(self.info.items()):
                if intervention == "E10IV" and k < 2:
                    continue

                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if "urine" in name else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{intervention}_cirrhosis",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO if route == "PO" else Route.IV,
                        application_form=ApplicationForm.TABLET
                        if route == "PO"
                        else ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.CIRRHOSIS,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                        genotype=Genotype.NR,
                        data=PKPDData.PK,
                        data_type=DataType.MEAN,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig1_tab1(),
        }

    def fig1_tab1(self) -> dict[str, Figure]:
        """Figure 1 and table 1B: plasma and urine, a figure per intervention."""
        figures = []
        for intervention in self.interventions:
            if intervention in ["LOS50PO", "LOS10IV"]:
                fig = Figure(
                    experiment=self,
                    sid=f"Fig1_Tab1_{intervention}",
                    num_rows=2,
                    num_cols=2,
                    name=f"{self.__class__.__name__} (cirrhosis)",
                )
                plots = fig.create_plots(
                    xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
                )
                plots[0].set_yaxis(self.label_los, unit=self.unit_los)
                plots[1].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
                plots[2].set_yaxis(self.label_e3174, unit=self.unit_e3174)
                plots[3].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)

            else:  # E10IV
                fig = Figure(
                    experiment=self,
                    sid=f"Fig1_Tab1_{intervention}",
                    num_rows=1,
                    num_cols=2,
                    name=f"{self.__class__.__name__} (cirrhosis)",
                )
                plots = fig.create_plots(
                    xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
                )
                plots[0].set_yaxis(self.label_e3174, unit=self.unit_e3174)
                plots[1].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)

            for k, (sid, name) in enumerate(self.info.items()):
                if intervention == "E10IV":
                    if k < 2:
                        continue
                    kp = k - 2

                else:
                    kp = k

                # simulation: black without cirrhosis, colored with cirrhosis
                for state in ["normal", "cirrhosis"]:
                    plots[kp].add_data(
                        task=f"task_{intervention}_{state}",
                        xid="time",
                        yid=sid,
                        label=f"Sim {intervention} ({state})",
                        color="black"
                        if state == "normal"
                        else self.colors[intervention],
                    )

                # data
                plots[kp].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd" if "urine" in name else None,
                    count="count",
                    label=f"{intervention} (cirrhosis)",
                    color=self.colors[intervention],
                )

            figures.append(fig)

        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(FDA1995S67, output_dir=FDA1995S67.__name__)
