"""Simulation experiment of the study Lo1995."""

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


class Lo1995(LosartanSimulationExperiment):
    """Simulation experiment of Lo1995.

    Two studies in healthy men after an overnight fast. Study 1 (n=12, 75.6
    kg): losartan potassium 50 mg as oral tablet and 20 mg as intravenous
    infusion over 20 min, the assayed doses 49.1 and 19.7 mg (Table I). Study
    2 (n=6, 78.6 kg): 14C-losartan 100 mg as oral solution and 30 mg as
    intravenous infusion over 20 min, the actual doses 90.3 and 26.4 mg
    losartan potassium (Table II), and E3174 20 mg (as hydrochloride) as
    intravenous infusion over 20 min, the actual dose 16.6 mg free acid
    (Table III). Losartan (as potassium salt) and E3174 (as free acid) in
    plasma (Fig 1, 2) and their amounts in urine (0-48 hr in study 1 and
    after E3174, 0-96 hr after losartan in study 2), the amounts in the copy
    are the recoveries (% of dose, molar) times the actual doses.
    """

    info_dict: ClassVar[dict[str, dict[str, str]]] = {
        "LOS50PO": {
            "[Cveplasma_los]": "losartan",
            "[Cveplasma_e3174]": "exp3174",
            "Aurine_los": "losartan_urine",
            "Aurine_e3174": "exp3174_urine",
        },
        "LOS20IV": {
            "[Cveplasma_los]": "losartan",
            "[Cveplasma_e3174]": "exp3174",
            "Aurine_los": "losartan_urine",
            "Aurine_e3174": "exp3174_urine",
        },
        "LOS30IV": {
            "[Cveplasma_los]": "losartan",
            "[Cveplasma_e3174]": "exp3174",
            "Aurine_los": "losartan_urine",
            "Aurine_e3174": "exp3174_urine",
        },
        "E20IV": {
            "[Cveplasma_e3174]": "exp3174",
            "Aurine_e3174": "exp3174_urine",
        },
        "LOS100PO": {
            "Aurine_los": "losartan_urine",
            "Aurine_e3174": "exp3174_urine",
        },
    }
    interventions: ClassVar[list[str]] = list(info_dict.keys())

    #: actual doses [mg], losartan potassium and E3174 free acid
    doses: ClassVar[dict[str, float]] = {
        "LOS50PO": 49.1,
        "LOS20IV": 19.7,
        "LOS30IV": 26.4,
        "E20IV": 16.6,
        "LOS100PO": 90.3,
    }
    forms: ClassVar[dict[str, ApplicationForm]] = {
        "LOS50PO": ApplicationForm.TABLET,
        "LOS20IV": ApplicationForm.SOLUTION,
        "LOS30IV": ApplicationForm.SOLUTION,
        "E20IV": ApplicationForm.SOLUTION,
        "LOS100PO": ApplicationForm.SOLUTION,
    }
    routes: ClassVar[dict[str, str]] = {
        "LOS50PO": "PO",
        "LOS20IV": "IV",
        "LOS30IV": "IV",
        "E20IV": "IV",
        "LOS100PO": "PO",
    }
    colors: ClassVar[dict[str, str]] = {
        "LOS50PO": "tab:blue",
        "LOS20IV": "tab:orange",
        "LOS30IV": "tab:red",
        "E20IV": "tab:green",
        "LOS100PO": "tab:brown",
    }

    studies: ClassVar[dict[str, str]] = {
        "LOS50PO": "study1",
        "LOS20IV": "study1",
        "LOS30IV": "study2",
        "E20IV": "study2",
        "LOS100PO": "study2",
    }

    bodyweights: ClassVar[dict[str, float]] = {
        "study1": 75.6,
        "study2": 78.6,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Tab1B", "Tab2B", "Tab3B"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)

                # unit conversion to mole/l, losartan is measured as potassium
                # salt, E3174 as free acid
                if fig_id in {"Fig1", "Fig2"}:
                    if label.startswith("losartan_"):
                        dset.unit_conversion("mean", 1 / self.Mr.losp)
                    elif label.startswith("exp3174_"):
                        dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            route = self.routes[intervention]
            substance = "los" if intervention.startswith("LOS") else "e3174"
            study = self.studies[intervention]
            dose = Q(self.doses[intervention], "mg")
            if substance == "los":
                dose = self.losartan_dose(dose)

            if route == "PO":
                simulations[intervention] = Simulation(
                    time_unit="min",
                    end=100 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "BW": Q(self.bodyweights[study], "kg"),
                        f"PODOSE_{substance}": dose,
                    },
                )
            else:
                # infusion over 20 min
                simulations[intervention] = Simulation(
                    time_unit="min",
                    end=100 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "BW": Q(self.bodyweights[study], "kg"),
                        f"Ri_{substance}": dose / Q(20, "min"),
                    },
                    changes=[Change(20, {f"Ri_{substance}": Q(0, "mg/min")})],
                )

        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            route = self.routes[intervention]

            for sid, name in self.info_dict[intervention].items():
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
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO if route == "PO" else Route.IV,
                        application_form=self.forms[intervention],
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                        genotype=Genotype.NR,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.fig_LOSIV(),
            **self.fig_LOSPO(),
            **self.fig_E3174IV(),
        }

    def fig_LOSIV(self) -> dict[str, Figure]:
        """Figure of the intravenous infusions of losartan."""
        fig = Figure(
            experiment=self,
            sid="Fig_LOSIV",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
        plots[3].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)
        for k in [0, 1]:
            xaxis = plots[k].xaxis
            if xaxis is not None:  # set by create_plots
                xaxis.min = -2
                xaxis.max = 35

        for intervention in ["LOS20IV", "LOS30IV"]:
            for k, (sid, name) in enumerate(self.info_dict[intervention].items()):
                # simulation
                plots[k].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {intervention}",
                    color=self.colors[intervention],
                )

                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd" if "urine" in name else None,
                    count="count",
                    label=f"{intervention}",
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)

    def fig_LOSPO(self) -> dict[str, Figure]:
        """Figure of the oral doses of losartan."""
        fig = Figure(
            experiment=self,
            sid="Fig_LOSPO",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
        plots[3].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)

        for k in [0, 1]:
            xaxis = plots[k].xaxis
            if xaxis is not None:  # set by create_plots
                xaxis.min = -2
                xaxis.max = 35

        for intervention in ["LOS50PO", "LOS100PO"]:
            # LOS100PO has only the urine, which are the panels 2 and 3
            start = 2 if intervention == "LOS100PO" else 0
            for k, (sid, name) in enumerate(
                self.info_dict[intervention].items(), start=start
            ):
                # simulation
                plots[k].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {intervention}",
                    color=self.colors[intervention],
                )

                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd" if "urine" in name else None,
                    count="count",
                    label=f"{intervention}",
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)

    def fig_E3174IV(self) -> dict[str, Figure]:
        """Figure of the intravenous infusion of E3174."""
        fig = Figure(
            experiment=self,
            sid="Fig_E3174IV",
            num_rows=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[1].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)
        xaxis = plots[0].xaxis
        if xaxis is not None:  # set by create_plots
            xaxis.min = -2
            xaxis.max = 35

        intervention = "E20IV"
        for k, (sid, name) in enumerate(self.info_dict[intervention].items()):
            # simulation
            plots[k].add_data(
                task=f"task_{intervention}",
                xid="time",
                yid=sid,
                label=f"Sim {intervention}",
                color=self.colors[intervention],
            )

            # data
            plots[k].add_data(
                dataset=f"{name}_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd" if "urine" in name else None,
                count="count",
                label=f"{intervention}",
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Lo1995, output_dir=Lo1995.__name__)
