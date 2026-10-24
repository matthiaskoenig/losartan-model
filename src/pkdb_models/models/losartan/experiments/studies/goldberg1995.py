"""Simulation experiment of the study Goldberg1995."""

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


class Goldberg1995(LosartanSimulationExperiment):
    """Simulation experiment of Goldberg1995.

    Single oral dose of losartan potassium 100 mg (MK-954) as tablet after an
    overnight fast in eight healthy male subjects, alone and on day 6 of
    cimetidine 400 mg four times daily, and a control day without losartan:
    losartan and E3174 in plasma (Merck HPLC assay, losartan as the potassium
    salt) and the plasma renin activity (PRA). The PRA is mapped as the ratio
    to the predose value of its day (`ren_ratio`); the PRA is proportional to
    renin at a constant angiotensinogen. Cimetidine (CYP3A4, losartan AUC x1.18,
    E3174 AUC x1.00) is not part of the model: its arm is compared to the
    simulation without cimetidine. The control day is flat in the model.
    """

    info: ClassVar[dict[str, str]] = {
        "[Cveplasma_los]": "losartan",
        "[Cveplasma_e3174]": "exp3174",
        "ren_ratio": "pra_ratio",
    }
    colors: ClassVar[dict[str, str]] = {
        "control": "black",
        "LOS100": "tab:blue",
        "LOS100_CIM400": "tab:orange",
    }
    #: the legend names of the arms
    legends: ClassVar[dict[str, str]] = {
        "control": "control",
        "LOS100": "losartan 100 mg",
        "LOS100_CIM400": "losartan 100 mg + cimetidine 400 mg",
    }
    losp_doses: ClassVar[dict[str, int]] = {
        "control": 0,  # [mg]
        "LOS100": 100,  # [mg]
    }
    #: the simulation of every arm of the data
    arms: ClassVar[dict[str, str]] = {
        "control": "control",
        "LOS100": "LOS100",
        "LOS100_CIM400": "LOS100",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label.endswith("_renin"):
                    # PRA x 7.6, not a measurement of renin
                    continue
                if label.endswith("_pra"):
                    # the PRA relative to the predose value of the day
                    df_ratio = df_label.copy()
                    baseline = df_label[df_label.time == 0]["mean"].to_numpy()[0]
                    df_ratio["mean"] = df_label["mean"] / baseline
                    df_ratio["unit"] = "dimensionless"
                    df_ratio["label"] = f"{label}_ratio"
                    dsets[f"{label}_ratio"] = DataSet.from_df(df_ratio, self.ureg)
                    continue

                dset = DataSet.from_df(df_label, self.ureg)
                # unit conversion to mole/l
                if label.endswith("_losartan"):
                    dset.unit_conversion("mean", 1 / self.Mr.losp)
                elif label.endswith("_exp3174"):
                    dset.unit_conversion("mean", 1 / self.Mr.e3174)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention, dose in self.losp_doses.items():
            simulations[f"po_{intervention}"] = Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "PODOSE_los": self.losartan_dose(Q(dose, "mg")),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.info.items():
            for arm, intervention in self.arms.items():
                if arm == "control" and sid != "ren_ratio":
                    continue
                if arm == "LOS100_CIM400" and sid == "ren_ratio":
                    # cimetidine is not part of the model
                    continue

                mappings[f"fm_po_{arm}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{arm}_{name}",
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
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.CIMETIDINE
                        if "CIM" in arm
                        else Coadministration.NONE,
                        data=PKPDData.PD if sid == "ren_ratio" else PKPDData.PK,
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
        """Figure of losartan and E3174 in plasma and of the PRA (Fig1, Fig2)."""
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_rows=1,
            num_cols=3,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_los, unit=self.unit_los)
        plots[1].set_yaxis(self.label_e3174, unit=self.unit_e3174)
        plots[2].set_yaxis("PRA ratio\n", unit=self.units["ren_ratio"])

        for k, (sid, name) in enumerate(self.info.items()):
            for intervention in self.losp_doses:
                if intervention == "control" and sid != "ren_ratio":
                    continue
                # simulation
                plots[k].add_data(
                    task=f"task_po_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"Sim {self.legends[intervention]}",
                    color=self.colors[intervention],
                )
            for arm in self.arms:
                if arm == "control" and sid != "ren_ratio":
                    continue
                if arm == "LOS100_CIM400" and sid == "ren_ratio":
                    continue
                # data
                plots[k].add_data(
                    dataset=f"{arm}_{name}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=self.legends[arm],
                    color=self.colors[arm],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Goldberg1995, output_dir=Goldberg1995.__name__)
