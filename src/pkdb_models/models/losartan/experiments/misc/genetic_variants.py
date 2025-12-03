from copy import deepcopy
from typing import Dict

from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)
from pkdb_models.models.losartan.helpers import run_experiments


class GeneticVariantsExperiment(LosartanSimulationExperiment):
    """Tests po application."""

    variant_definitions = {
        "CYP2C9": {
            "parameter": "LI__f_cyp2c9",
            "activities": LosartanSimulationExperiment.cyp2c9_activity,
            "colors": LosartanSimulationExperiment.cyp2c9_colors,
        },
        "ABCB1": {
            "parameter": "GU__f_abcb1",
            "activities": LosartanSimulationExperiment.abcb1_activity,
            "colors": LosartanSimulationExperiment.abcb1_colors,
        },
    }

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for key, info in self.variant_definitions.items():
            for variant, activity in info["activities"].items():
                tcsims[f"los_{variant}"] = TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,  # [min]
                        steps=300,
                        changes={
                            **self.default_changes(),
                            f"PODOSE_los": Q_(50, "mg"),
                            info["parameter"]: Q_(activity, "dimensionless"),
                        },
                    )
                )

        return tcsims

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> Dict[str, Figure]:

        figures = {}
        for key, info in self.variant_definitions.items():

            fig = Figure(
                experiment=self,
                sid=f"Fig_genetic_variants_pk_{key}",
                num_rows=8,
                num_cols=4,
                name=f"Genetic variants ({key})"
            )
            plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)
            for k in [30, 31]:
                plots[k].xaxis = None
            sids = [
                # plasma
                "[Cve_los]",
                "[Cve_e3174]",
                "[Cve_l158]",
                "mr_e3174_los_plasma",
                # "mr_los_e3174_plasma",

                # urine,
                "Aurine_los",
                "Aurine_e3174",
                "Aurine_l158",
                "mr_e3174_los_urine",
                # "mr_los_e3174_urine",

                # feces
                "Afeces_los",
                "Afeces_e3174",
                "Afeces_l158",
                "mr_e3174_los_feces",
                # "mr_los_e3174_feces",

                "[e3174]",  # e3174 in Vplasma [0-1000 nM]
                "fe_e3174",  # [-] effect via exp3174
                None,

                "[anggen]",  # angiotensinogen in Vplasma
                "[ang1]",  # angiotensin I in Vplasma
                "[ang2]",  # angiotensin II in Vplasma

                "ANGGEN2ANG1",  # angiotensinogen to angiotensin I renin)
                "ANG1ANG2",  # angiotensin I to angiotensin II ACE)
                "ANG2DEG",  # angiotensin II degradation ANG2DEG)

                "[ren]",  # renin in Vplasma
                "RENSEC",  # renin secretion RENSEC)
                "RENDEG",  # renin degradation RENDEG)

                "[ald]",  # aldosterone in Vplasma
                "ALDSEC",  # aldosterone secretion ALDSEC)
                "ALDDEG",  # aldosterone degradation (ALDDEG),

                "SBP",
                "DBP",
                "MAP",

            ]
            for ksid, sid in enumerate(sids):
                if sid:
                    plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            for ksid, sid in enumerate(sids):
                if sid:
                    for variant in info["activities"]:
                        plots[ksid].add_data(
                            task=f"task_los_{variant}",
                            xid="time",
                            yid=sid,
                            label=variant,
                            color=info["colors"][variant],
                        )

            figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    run_experiments(GeneticVariantsExperiment, output_dir=GeneticVariantsExperiment.__name__)
