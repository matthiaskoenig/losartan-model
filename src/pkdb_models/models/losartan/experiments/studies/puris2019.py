"""Simulation experiment of the study Puris2019."""

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


class Puris2019(LosartanSimulationExperiment):
    """Simulation experiment of Puris2019.

    Single oral dose of 12.5 mg losartan potassium with a cocktail of probe
    drugs (melatonin, nicotine, bupropion, repaglinide, omeprazole,
    chlorzoxazone, midazolam, dextromethorphan), given as capsules compounded
    from crushed tablets, in 8 obese patients (BMI 41.0 kg/m2) before and one
    year after Roux-en-Y gastric bypass surgery (BMI 34.2 kg/m2) and in 6
    non-obese healthy controls: losartan and E3174 in serum at 1, 2, 4 and 6 h
    (Fig 1, mean and SEM) and the 6 h urinary recovery of the patients (Tab S4,
    geometric mean).

    The groups differ in the model only by the body weight: the study reports
    the BMI, not the weight, the weight of the patients is the BMI times the
    square of a height of 1.685 m (6 women of 1.65 m and 2 men of 1.79 m, the
    mean heights of Finnish adults); the controls have the reference weight.
    The 6 h urinary losartan/E3174 ratio of Table 3 is a summary of the same
    urine samples and not used.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
    }
    info_urine: ClassVar[dict[str, str]] = {
        "Aurine_los": "losartan_urine",
        "Aurine_e3174": "exp3174_urine",
    }

    groups: ClassVar[list[str]] = ["healthy", "obese", "one_year_after"]
    #: the groups with urine data (Tab S4)
    groups_urine: ClassVar[list[str]] = ["obese", "one_year_after"]

    health: ClassVar[dict[str, Health]] = {
        "healthy": Health.HEALTHY,
        "obese": Health.OBESE,
        "one_year_after": Health.OBESE,  # BMI 34.2 kg/m2 after the surgery
    }

    #: body weight [kg] from the BMI and a height of 1.685 m, the controls
    #: (BMI 18-25 kg/m2, sex not reported) have the reference weight
    bodyweights: ClassVar[dict[str, float | None]] = {
        "healthy": None,
        "obese": 41.0 * 1.685**2,  # 116.4 kg
        "one_year_after": 34.2 * 1.685**2,  # 97.1 kg
    }

    colors: ClassVar[dict[str, str]] = {
        "healthy": "black",
        "obese": "tab:blue",
        "one_year_after": "tab:green",
    }
    #: the labels of the groups in the figures
    group_labels: ClassVar[dict[str, str]] = {
        "healthy": "healthy",
        "obese": "obese",
        "one_year_after": "one year after surgery",
    }

    dose = Q(12.5, "mg")  # losartan potassium

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig1")
        for label, df_label in df.groupby("label"):
            label = str(label)
            # the serum concentrations of the cocktail drugs, mean and SEM
            if not label.startswith(("losartan", "exp3174")):
                continue
            dset = DataSet.from_df(df_label, self.ureg)
            if label.startswith("losartan"):
                dset.unit_conversion("mean", 1 / self.Mr.los)
            elif label.startswith("exp3174"):
                dset.unit_conversion("mean", 1 / self.Mr.e3174)
            dsets[self.to_sid(label)] = dset

        # urinary recovery in 6 h [% of the dose], geometric mean (range), as
        # amount of the molar dose of losartan potassium
        dose_umole = (self.dose / self.Mr.losp).to("µmole").magnitude
        df = self.load_dataframe("TabS4")
        for substance in ["losartan", "exp3174"]:
            for group in self.groups_urine:
                df_group = df[(df.substance == substance) & (df.group == group)].copy()
                for key in ["geomean", "min", "max"]:
                    df_group[key] = df_group[key] / 100 * dose_umole
                df_group["mean"] = df_group["geomean"]
                df_group["unit"] = "µmole"
                dsets[self.to_sid(f"{substance}_urine_{group}")] = DataSet.from_df(
                    df_group, self.ureg
                )

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for group in self.groups:
            bodyweight = self.bodyweights[group]
            tcsims[f"po_los12p5_{group}"] = Simulation(
                time_unit="min",
                end=12 * 60,
                preinit_changes={
                    **self.default_changes(),
                    **({"BW": Q(bodyweight, "kg")} if bodyweight else {}),
                    "PODOSE_los": self.losartan_dose(self.dose),
                },
            )
        return tcsims

    def _metadata(self, group: str, tissue: Tissue) -> PKPDMappingMetaData:
        """Metadata of the mappings of a group."""
        return PKPDMappingMetaData(
            tissue=tissue,
            route=Route.PO,
            application_form=ApplicationForm.CAPSULE,
            dosing=Dosing.SINGLE,
            health=self.health[group],
            fasting=Fasting.NR,
            coadministration=Coadministration.COCKTAIL,
        )

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for group in self.groups:
                mappings[self.to_sid(f"fm_po_los12p5_{name}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_los12p5_{group}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=self._metadata(group, Tissue.SERUM),
                )
        for sid, name in self.info_urine.items():
            for group in self.groups_urine:
                mappings[self.to_sid(f"fm_po_los12p5_{name}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_po_los12p5_{group}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=self._metadata(group, Tissue.URINE),
                )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_TabS4",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__} (obesity, cocktail)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis(self.label_los_urine, unit=self.unit_los_urine)
        plots[3].set_yaxis(self.label_e3174_urine, unit=self.unit_e3174_urine)

        for k, (sid, name) in enumerate({**self.info, **self.info_urine}.items()):
            urine = sid in self.info_urine
            for group in self.groups:
                # simulation
                plots[k].add_data(
                    task=f"task_po_los12p5_{group}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {self.group_labels[group]}",
                    color=self.colors[group],
                )
                if urine and group not in self.groups_urine:
                    continue
                # data
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{group}"),
                    xid="time",
                    yid="mean",
                    yid_sd=None if urine else "mean_sd",
                    count="count",
                    label=self.group_labels[group],
                    color=self.colors[group],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Puris2019, output_dir=Puris2019.__name__)
