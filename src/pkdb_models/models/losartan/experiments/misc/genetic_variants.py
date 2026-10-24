"""Effect of the genetic variants of CYP2C9 and ABCB1 on losartan."""

from typing import ClassVar, NamedTuple, override

from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.losartan import MODEL
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)


class Variants(NamedTuple):
    """The variants of a gene: the parameter of the activity and its values."""

    parameter: str
    activities: dict[str, float]
    colors: dict[str, str]


class GeneticVariantsExperiment(LosartanSimulationExperiment):
    """Single oral dose of 50 mg losartan with the genotypes of CYP2C9 and ABCB1."""

    variant_definitions: ClassVar[dict[str, Variants]] = {
        "CYP2C9": Variants(
            parameter="LI__f_cyp2c9",
            activities=LosartanSimulationExperiment.cyp2c9_activity,
            colors=LosartanSimulationExperiment.cyp2c9_colors,
        ),
        "ABCB1": Variants(
            parameter="GU__f_abcb1",
            activities=LosartanSimulationExperiment.abcb1_activity,
            colors=LosartanSimulationExperiment.abcb1_colors,
        ),
    }

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        for info in self.variant_definitions.values():
            for variant, activity in info.activities.items():
                tcsims[self.to_sid(f"los_{variant}")] = Simulation(
                    time_unit="min",
                    end=24 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "PODOSE_los": Q(50, "mg"),
                        info.parameter: Q(activity, "dimensionless"),
                    },
                )

        return tcsims

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> dict[str, Figure]:
        """Figures of the pharmacokinetics and pharmacodynamics by gene."""
        figures = []
        for key, info in self.variant_definitions.items():
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_genetic_variants_pk_{key}"),
                num_rows=8,
                num_cols=4,
                name=f"Genetic variants ({key})",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit="hr"), legend=True
            )
            for k in [30, 31]:
                plots[k].xaxis = None
            sids = [
                # plasma
                "[Cveplasma_los]",
                "[Cveplasma_e3174]",
                "[Cveplasma_l158]",
                "mr_e3174_los_plasma",
                # urine,
                "Aurine_los",
                "Aurine_e3174",
                "Aurine_l158",
                "mr_e3174_los_urine",
                # feces
                "Afeces_los",
                "Afeces_e3174",
                "Afeces_l158",
                "mr_e3174_los_feces",
                "[e3174]",  # e3174 in Vplasma [0-1000 nM]
                "fe_e3174",  # [-] effect via exp3174
                "ren_ratio",  # with ald_ratio, see below
                "[anggen]",  # angiotensinogen in Vplasma
                "[ang1]",  # angiotensin I in Vplasma
                "[ang2]",  # angiotensin II in Vplasma
                "ANGGEN2ANG1",  # angiotensinogen to angiotensin I (renin)
                "ANG1ANG2",  # angiotensin I to angiotensin II (ACE)
                "ANG2DEG",  # angiotensin II degradation (ANG2DEG)
                "[ren]",  # renin in Vplasma
                "RENSEC",  # renin secretion (RENSEC)
                "RENDEG",  # renin degradation (RENDEG)
                "[ald]",  # aldosterone in Vplasma
                "ALDSEC",  # aldosterone secretion (ALDSEC)
                "ALDDEG",  # aldosterone degradation (ALDDEG)
                "SBP",
                "DBP",
                "MAP",
            ]
            for ksid, sid in enumerate(sids):
                plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            # a variant with the activity of an earlier one (CYP2C9 *1/*13 and
            # *1/*3) is dashed, its curves would hide those of the earlier one
            linestyles: dict[str, str] = {}
            for variant, activity in info.activities.items():
                earlier = {info.activities[v] for v in linestyles}
                linestyles[variant] = "--" if activity in earlier else "-"

            for ksid, sid in enumerate(sids):
                for variant in info.activities:
                    plots[ksid].add_data(
                        task=self.to_sid(f"task_los_{variant}"),
                        xid="time",
                        yid=sid,
                        label=variant,
                        color=info.colors[variant],
                        linestyle=linestyles[variant],
                    )

            # the renin secretion is activated, the aldosterone secretion
            # inhibited: their ratios to the baseline share a panel
            kratio = sids.index("ren_ratio")
            plots[kratio].set_yaxis(
                label="renin|aldosterone ratio", unit=self.units["ren_ratio"]
            )
            for variant in info.activities:
                plots[kratio].add_data(
                    task=self.to_sid(f"task_los_{variant}"),
                    xid="time",
                    yid="ald_ratio",
                    label=None,
                    color=info.colors[variant],
                    linestyle=linestyles[variant],
                )

            figures.append(fig)

        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(
        GeneticVariantsExperiment, output_dir=GeneticVariantsExperiment.__name__
    )
