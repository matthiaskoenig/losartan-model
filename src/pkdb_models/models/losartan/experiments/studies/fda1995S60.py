"""Simulation experiment of the study FDA1995S60."""

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


class FDA1995S60(LosartanSimulationExperiment):
    """Simulation experiment of FDA1995S60.

    Study 060 of the FDA review of NDA 020386 (the study 2 of Lo1995, the same
    6 fasted healthy male subjects, 78.6 kg): 14C-losartan potassium 100 mg as
    an oral solution (30 mg labelled and 70 mg unlabelled) and 30 mg as a
    20 min intravenous infusion, and E3174 20 mg as a 20 min intravenous
    infusion. Losartan and E3174 in plasma (Fig1), the radioactivity in urine
    (96 hr) and feces (168 hr) and its losartan, E3174 and L-158 fractions
    (Tab60_2).

    The doses are the actual doses of Lo1995 (Tab2, Tab4): 90.3 mg (oral) and
    26.4 mg (intravenous) losartan potassium and 16.6 mg E3174 (free acid); the
    plasma concentrations are not normalized to the nominal doses. Losartan is
    reported as the potassium salt (HPLC of Merck, Lo1995). The amounts in
    urine and feces are the fractions of the dose (radioactivity, Tab60_2)
    times the actual dose in moles (the curated amounts use 30 mg for the oral
    dose, the nominal doses and the molecular weight of losartan).

    The L-158 fraction (L-158,795, L-158,796 and L-158,783) is a part of the
    other metabolites of the model (`l158`, all metabolites but E3174); the
    losartan and E3174 in feces after the intravenous dose (2 subjects) come
    from a biliary excretion of losartan, which the model does not have.
    """

    #: the actual doses [mg]: losartan potassium, E3174 free acid
    doses: ClassVar[dict[str, float]] = {
        "LOS100PO": 90.3,
        "LOS30IV": 26.4,
        "E20IV": 16.6,
    }
    routes: ClassVar[dict[str, str]] = {
        "LOS100PO": "PO",
        "LOS30IV": "IV",
        "E20IV": "IV",
    }
    substances: ClassVar[dict[str, str]] = {
        "LOS100PO": "los",  # 14C losartan
        "LOS30IV": "los",  # 14C losartan
        "E20IV": "e3174",
    }
    #: the duration of the intravenous infusions
    infusion_time: ClassVar[float] = 20.0  # [min]
    colors: ClassVar[dict[str, str]] = {
        "LOS100PO": "tab:blue",
        "LOS30IV": "tab:red",
        "E20IV": "tab:green",
    }
    interventions: ClassVar[list[str]] = list(doses.keys())

    infos: ClassVar[dict[str, dict[str, str | None]]] = {
        "LOS100PO": {
            "[Cveplasma_los]": "losartan",
            "Aurine_los": "losartan_urine",
            "Afeces_los": None,
            "[Cveplasma_e3174]": "exp3174",
            "Aurine_e3174": "exp3174_urine",
            "Afeces_e3174": None,
            "[Cveplasma_l158]": None,
            "Aurine_l158": "l158_urine",
            "Afeces_l158": None,
            "[Cveplasma_total]": None,
            "Aurine_total": "losartan_total_urine",
            "Afeces_total": "losartan_total_feces",
        },
        "LOS30IV": {
            "[Cveplasma_los]": "losartan",
            "Aurine_los": "losartan_urine",
            "Afeces_los": "losartan_feces",
            "[Cveplasma_e3174]": "exp3174",
            "Aurine_e3174": "exp3174_urine",
            "Afeces_e3174": "exp3174_feces",
            "[Cveplasma_l158]": None,
            "Aurine_l158": "l158_urine",
            "Afeces_l158": "l158_feces",
            "[Cveplasma_total]": None,
            "Aurine_total": "losartan_total_urine",
            "Afeces_total": "losartan_total_feces",
        },
        "E20IV": {
            "[Cveplasma_e3174]": "exp3174",
        },
    }

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

        # urine and feces: the fraction of the dose (radioactivity) times the
        # actual dose of losartan potassium in moles
        df = self.load_dataframe("Tab60_2B")
        for label, df_label in df.groupby("label"):
            label = str(label)
            if "_recovery_" not in label:
                # the curated amounts, see the class
                continue
            intervention = label.split("_")[-1]
            amount = (Q(self.doses[intervention], "mg") / self.Mr.losp).to("µmole")
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
        for intervention in self.interventions:
            route = self.routes[intervention]
            substance = self.substances[intervention]
            dose = Q(self.doses[intervention], "mg")
            if substance == "los":
                dose = self.losartan_dose(dose)

            changes = []
            if route == "PO":
                preinit_changes = {"PODOSE_los": dose}
            else:
                rate = (dose / Q(self.infusion_time, "min")).to("mg/min")
                preinit_changes = {f"Ri_{substance}": rate}
                changes = [
                    Change(self.infusion_time, {f"Ri_{substance}": Q(0.0, "mg/min")})
                ]

            tcsims[f"{intervention}"] = Simulation(
                time_unit="min",
                end=170 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(78.6, "kg"),  # value from Lo1995
                    **preinit_changes,
                },
                changes=changes,
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}

        for intervention in self.interventions:
            route = self.routes[intervention]
            info = self.infos[intervention]
            for sid, name in info.items():
                if not name:
                    continue
                tissue = Tissue.PLASMA
                if "urine" in name:
                    tissue = Tissue.URINE
                elif "feces" in name:
                    tissue = Tissue.FECES

                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=tissue,
                        route=Route.PO if route == "PO" else Route.IV,
                        # an oral solution, an intravenous infusion
                        application_form=ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
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
            **self.fig_los(),
            **self.fig_e3174(),
        }

    def fig_los(self) -> dict[str, Figure]:
        """Figure 1 and table 2B: oral and iv losartan."""
        figures = []
        for intervention in ["LOS100PO", "LOS30IV"]:
            info = self.infos[intervention]
            fig = Figure(
                experiment=self,
                sid=f"Fig1_{intervention}",
                num_cols=3,
                num_rows=4,
                name=f"{self.__class__.__name__} (healthy)",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            plots[0].set_yaxis(self.label_los, unit=self.unit_los)
            plots[1].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
            plots[2].set_yaxis(self.label_los_feces, unit=self.unit_los_feces)

            plots[3].set_yaxis(self.label_e3174, unit=self.unit_e3174)
            plots[4].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)
            plots[5].set_yaxis(self.label_e3174_feces, unit=self.unit_e3174_feces)

            plots[6].set_yaxis(self.label_l158, unit=self.unit_l158)
            plots[7].set_yaxis(self.label_l158_urine, unit=self.unit_l158_urine)
            plots[8].set_yaxis(self.label_l158_feces, unit=self.unit_l158_feces)

            plots[9].set_yaxis(self.label_total, unit=self.unit_total)
            plots[10].set_yaxis(self.label_total_urine, unit=self.unit_total_urine)
            plots[11].set_yaxis(self.label_total_feces, unit=self.unit_total_feces)

            for k in [0, 3, 6, 9]:
                xaxis = plots[k].xaxis
                if xaxis is not None:  # set by `create_plots`
                    xaxis.max = 32

            for k, (sid, name) in enumerate(info.items()):
                # simulation
                plots[k].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {intervention}",
                    color=self.colors[intervention],
                )
                # data
                if name:
                    plots[k].add_data(
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                        label=f"{intervention}",
                        color=self.colors[intervention],
                    )
            figures.append(fig)

        return self.figures_by_sid(*figures)

    def fig_e3174(self) -> dict[str, Figure]:
        """Figure 1: iv E3174."""
        intervention = "E20IV"
        info = self.infos[intervention]
        fig = Figure(
            experiment=self,
            sid=f"Fig1_{intervention}",
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, max=32), legend=True
        )
        plots[0].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        for k, (sid, name) in enumerate(info.items()):
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
                yid_sd=None,
                count="count",
                label=f"{intervention}",
                color=self.colors[intervention],
            )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(FDA1995S60, output_dir=FDA1995S60.__name__)
