from typing import Dict

import numpy as np
from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.losartan.experiments.raas_experiment import RaasSimulationExperiment
from pkdb_models.models.losartan.helpers import run_experiments


class RaasExperiment(RaasSimulationExperiment):
    """Effect of exp3174."""

    e3174_values = np.linspace(0.0, 500.0, num=6)  # [nM]
    colors = ['black', '#eff3ff','#bdd7e7','#6baed6','#3182bd','#08519c']

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        ren_ref = Q_(5, "pM")
        ang2_ref = Q_(20, "pM")
        ald_ref = Q_(400, "pM")

        baseline_changes = {
            # "ren_ref": ren_ref,
            # "[ren]": ren_ref,
            # "ang2_ref": ang2_ref,
            # "[ang2]": ang2_ref,
            # "ald_ref": ald_ref,
            # "[ald]": ald_ref,
            #
            # "SBP_ref": Q_(140, "mmHg"),
            # "DBP_ref": Q_(100, "mmHg"),
        }

        for k, e3174 in enumerate(self.e3174_values):
            tcsims[f"raas_{k}"] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=5 * 60,  # [min] # simulate 1 day
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "[e3174]": Q_(0, "nM"),
                        **baseline_changes,
                    },
                ),
                Timecourse(
                    start=0,
                    end=24 * 60,  # [min] # simulate 1 day
                    steps=2000,
                    changes={
                        "[e3174]": Q_(e3174, "nM")
                    },
                ),
                Timecourse(
                    start=0,
                    end=24 * 60,  # [min] # simulate 1 day
                    steps=2000,
                    changes={
                        "[e3174]": Q_(0, "nM")
                    },
                ),
                ],
                time_offset=-5*60
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:

        info = [
            ("[e3174]", 0),  # e3174 in Vplasma [0-1000 nM]
            ("fa_e3174", 1),  # [-] activation via exp3174
            ("fi_e3174", 2),  # [-] inhibition via exp3174

            ("[anggen]", 3),  # angiotensinogen in Vplasma
            ("[ang1]", 4),  # angiotensin I in Vplasma
            ("[ang2]", 5),  # angiotensin II in Vplasma

            ("ANGGEN2ANG1", 6),  # angiotensinogen to angiotensin I (renin)
            ("ANG1ANG2", 7),  # angiotensin I to angiotensin II (ACE)
            ("ANG2DEG", 8),  # angiotensin II degradation (ANG2DEG)

            ("[ren]", 9),  # renin in Vplasma
            ("RENSEC", 10),  # renin secretion (RENSEC)
            ("RENDEG", 11),  # renin degradation (RENDEG)

            ("[ald]", 12),  # aldosterone in Vplasma
            ("ALDSEC", 13),  # aldosterone secretion (ALDSEC)
            ("ALDDEG", 14),  # aldosterone degradation (ALDDEG)

            ("SBP", 15),
            ("MAP", 15),
            ("DBP", 15),

            ("SBP_change", 16),
            # ("MAP_change", 16),
            ("DBP_change", 16),

            ("SBP_ratio", 17),
            # ("MAP_ratio", 17),
            ("DBP_ratio", 17),
        ]

        fig = Figure(
            experiment=self,
            sid=f"Fig_RAAS",
            num_rows=6,
            num_cols=3,
            name=f"RAAS",
        )
        plots = fig.create_plots(xaxis=Axis("time", unit="hour"), legend=True)

        for sid, ksid in info:
            plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])
            if sid in ["fa", "fi"]:
                yaxis = plots[ksid].yaxis
                yaxis.min = -0.02
                yaxis.max = 1.02
            for k, e3174 in enumerate(self.e3174_values):
                plots[ksid].add_data(
                    task=f"task_raas_{k}",
                    xid="time",
                    yid=sid,
                    label=f"{e3174:.1f} nM" if sid not in {"MAP", "DBP"} else None,
                    color=self.colors[k],
                )

        return {fig.sid: fig}


if __name__ == "__main__":
    run_experiments(
       RaasExperiment, output_dir=RaasExperiment.__name__
    )
