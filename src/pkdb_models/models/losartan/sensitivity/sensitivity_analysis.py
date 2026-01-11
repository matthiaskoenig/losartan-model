"""Sensitivity analysis."""
from __future__ import annotations

import dill
import numpy as np
import pandas as pd
import roadrunner
from roadrunner._roadrunner import NamedArray

from sbmlsim.sensitivity.analysis import (
    SensitivitySimulation,
    SensitivityOutput,
    LocalSensitivityAnalysis,
    SobolSensitivityAnalysis,
    SamplingSensitivityAnalysis, AnalysisGroup,
)

from sbmlsim.sensitivity.parameters import (
    SensitivityParameter,
    parameters_for_sensitivity_analysis, ParameterType,
)
from sbmlutils.console import console
from pkdb_analysis.pk.pharmacokinetics import TimecoursePK
from pint import UnitRegistry


class LosartanSensitivitySimulation(SensitivitySimulation):
    """Simulation for sensitivity calculation."""


    def simulate(self, r: roadrunner.RoadRunner, changes: dict[str, float]) -> dict[str, float]:
        tend = 20 * 24 * 60  # [min]
        steps = 3000

        # apply changes and simulate
        all_changes = {
            **self.changes_simulation,
            **changes
        }
        self.apply_changes(r, all_changes, reset_all=True)
        s: NamedArray = r.simulate(start=0, end=tend, steps=steps)

        # pharmacokinetic parameters
        y: dict[str, float] = {}

        # pharmacokinetics
        ureg = UnitRegistry()
        Q_ = ureg.Quantity

        # losartan
        Mr_los = Q_(461.0, "g/mole")
        time = Q_(s["time"], "min")
        tcpk = TimecoursePK(
            time=time,
            concentration=Q_(s["[Cve_los]"], "mM"),
            substance="losartan",
            ureg=ureg,
            dose=Q_(10, "mg")/Mr_los,
        )
        pk_dict = tcpk.pk.to_dict()
        # console.print(pk_dict)
        for pk_key in [
            "aucinf",
            "cmax",
            "thalf",
            "vd",
            "cl",
            "kel",
        ]:
            y[f"[Cve_los]_{pk_key}"] = pk_dict[pk_key]

        # E3174, L158
        for sid in [
            "[Cve_e3174]",
            "[Cve_l158]"
        ]:
            tcpk = TimecoursePK(
                time=time,
                concentration=Q_(s[sid], "mM"),
                substance="losartan",
                ureg=ureg,
                dose=None,
            )
            pk_dict = tcpk.pk.to_dict()
            # console.print(pk_dict)
            for pk_key in [
                "aucinf",
                "cmax",
                "thalf",
                # "kel",
            ]:
                y[f"{sid}_{pk_key}"] = pk_dict[pk_key]

        # pharmacodynamics
        for sid, f in [
            ("[ang1]", "max"),
            ("[ang2]", "max"),
            ("[ren]", "max"),
            ("[ald]", "min"),
            ("SBP", "min"),
            ("DBP", "min"),
            ("MAP", "min"),
        ]:
            # minimal and maximal value of readout
            if f == "max":
                y[f"{sid}_max"] = np.max(s[sid])
            elif f == "min":
                y[f"{sid}_min"] = np.min(s[sid])

        return y

    def _plot(self, s: NamedArray) -> None:

        # plotting
        from matplotlib import pyplot as plt
        plt.plot(
            s["time"], s["[Cve_los]"],
            # df["time"], df["[Cve_los]"],
            marker="o",
            markeredgecolor="black",
            color="tab:blue",
        )
        plt.show()

    @staticmethod
    def sensitivity_simulation() -> LosartanSensitivitySimulation:
        sensitivity_simulation = LosartanSensitivitySimulation(
            model_path=MODEL_PATH,
            selections=[
                "time",
                "[Cve_los]",
                "[Cve_e3174]",
                "[Cve_l158]",
                "[ang1]",
                "[ang2]",
                "[ren]",
                "[ald]",
                "SBP",
                "DBP",
                "MAP",
            ],
            changes_simulation = {
                "PODOSE_los": 10,  # [mg]
                "f_cirrhosis": 0,  # [-]
                "KI__f_renal_function": 1.0,  # [-]
            },
            outputs=[
                SensitivityOutput(uid='[Cve_los]_aucinf', name='LOS AUC∞', unit="mM*min"),
                SensitivityOutput(uid='[Cve_los]_cmax', name='LOS Cmax', unit="mM"),
                SensitivityOutput(uid='[Cve_los]_thalf', name='LOS Half-life', unit="min"),
                SensitivityOutput(uid='[Cve_los]_vd', name='LOS Vd', unit="l"),
                SensitivityOutput(uid='[Cve_los]_cl', name='LOS CL', unit="mole/min/mM"),
                SensitivityOutput(uid='[Cve_los]_kel', name='LOS kel', unit="1/min"),

                SensitivityOutput(uid='[Cve_e3174]_aucinf', name='E3174 AUC∞', unit="mM*min"),
                SensitivityOutput(uid='[Cve_e3174]_cmax', name='E3174 Cmax', unit="mmole/l"),
                SensitivityOutput(uid='[Cve_e3174]_thalf', name='E3174 Half-life', unit="min"),
                SensitivityOutput(uid='[Cve_l158]_aucinf', name='L158 AUC∞', unit="mM*min"),
                SensitivityOutput(uid='[Cve_l158]_cmax', name='L158 Cmax', unit="mmole/l"),
                SensitivityOutput(uid='[Cve_l158]_thalf', name='L158 Half-life', unit="min"),

                SensitivityOutput(uid='[ang1]_max', name='angiotensin 1 max', unit="mmole/l"),
                SensitivityOutput(uid='[ang2]_max', name='angiotensin 2 max', unit="mmole/l"),
                SensitivityOutput(uid='[ren]_max', name='renin max', unit="mmole/l"),
                SensitivityOutput(uid='[ald]_min', name='aldosterone min', unit="mmole/l"),
                SensitivityOutput(uid='SBP_min', name='SBP min', unit="mmHg"),
                SensitivityOutput(uid='DBP_min', name='DBP min', unit="mmHg"),
                SensitivityOutput(uid='MAP_min', name='MAP min', unit="mmHg"),
            ]
        )
        console.rule("Outputs", style="white")
        console.print(sensitivity_simulation.outputs)

        return sensitivity_simulation

    def sensitivity_parameters(self) -> list[SensitivityParameter]:
        """Definition of parameters and bounds for sensitivity analysis."""
        console.rule("Parameters", style="white")
        # parameters for sensitivity analysis
        parameters: list[SensitivityParameter] = parameters_for_sensitivity_analysis(
            sbml_path=self.model_path,
            exclude_ids={
                # conversion factors
                "conversion_min_per_day",

                # molecular weights
                "Mr_los",
                "Mr_e3174",
                "Mr_l158",

                # unchangable values
                "FQlu",
                "FVhv",
                "FVpo",

                # fast transports
                "LI__L158EX_k",
                "LI__LOSIM_k",

                # dosing parameters
                "PODOSE_los",
                "ti_los",
                "ti_e3174",

                # unused volumes
                "Vurine",
                "Vfeces",
                "Vplasma",
                "Vstomach",
                "LI__Vbi",

                # reference values
                "SBP_ref",
                "DBP_ref",
                "ald_ref",
                "ren_ref",
                "ang1_ref",
                "ang2_ref",
            },
            exclude_na=True,
            exclude_zero=True,
        )

        # bounds from fitted parameters
        fit_bounds = {
            'ftissue_los': [0.01, 10],
            'Kp_los': [1, 200],
            'GU__LOSABS_k': [0.0001, 1],
            'GU__f_LOSEFL_k': [0.1, 10],
            'GU__METEXC_k': [1e-05, 0.1],
            'LI__E3174EX_k': [0.001, 10],
            'LI__LOS2E3174_Vmax': [1e-05, 100],
            'LI__E3174L158_k': [1e-05, 10],
            'LI__MBIEX_k': [1e-05, 1],
            'KI__LOSEX_k': [0.0001, 1],
            'KI__E3174EX_k': [0.0001, 1],
            'KI__L158EX_k': [0.0001, 1],
            'ANGGEN2ANG1_k': [0.01, 100000.0],
            'E50_e3174': [5e-07, 0.05],
            'ALDSEC_k': [1e-06, 1],
            'BP_ald_fe': [0.1, 0.6],
        }
        fit_bounds = [
            (k, v[0], v[1], ParameterType.FIT) for k, v in fit_bounds.items()
        ]
        SensitivityParameter.parameters_set_bounds(parameters, bounds=fit_bounds)

        # bounds from scaled parameters
        bounds_fraction = 0.15  # fraction of bounds relative to value
        uids_scaling = [
            "GU__f_OATP2B1",
            "GU__f_abcb1",
            "LI__f_cyp2c9",
            "KI__f_renal_function",
        ]
        scaling_bounds = [
            (uid, 1 - bounds_fraction, 1 + bounds_fraction, ParameterType.SCALING) for uid in uids_scaling
        ]
        SensitivityParameter.parameters_set_bounds(parameters, bounds=scaling_bounds)

        # references for values
        reference_data={
            "HCT": r"\cite{Mondal2025, Fiseha2023}",
            "BW": r"\cite{Ogden2004, Jones2013, Thompson2009, Brown1997}",
            "FQgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
            "FQh": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
            "FQki": r"\cite{Jones2013, Thompson2009, Brown1997}",
            "FVar": r"\cite{Jones2013, Thompson2009, Brown1997}",
            "FVgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
            "FVki": r"\cite{Jones2013, Thompson2009, Brown1997}",
            "FVli": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
            "FVlu": r"\cite{Jones2013, Thompson2009, Brown1997}",
            "FVve": r"\cite{Jones2013, Thompson2009, Brown1997}",
            "LI__LOS2E3174_Km_los": r"\cite{Maekawa2009, Thu2017, Wang2014}",
            "COBW": r"\cite{Cattermole2017, Patel2021, Collis2001}",
            "HR": r"\cite{Cattermole2017, Patel2021, Collis2001, Thompson2009}",
            "GU__F_los_abs": r"\cite{Lo1995}", # r"60\% fraction absorbed (oral bioavailability 35\%, due to liver metabolism)~\cite{Lo1995}",
            "GU__f_abcb1": r"\cite{Siegmund2002, Hoffmeyer2000}",
            "KI__f_renal_function": r"\cite{Stevens2024}",
            "LI__f_cyp2c9": r"\cite{Kusama2009, Wang2014, Maekawa2009}",
        }
        p_dict = {p.uid: p for p in parameters}
        for pid, reference in reference_data.items():
            p = p_dict[pid]
            p.reference = reference
            if p.type == ParameterType.NA:
                p.type = ParameterType.DATA

        # setting missing bounds;
        for p in parameters:
            if np.isnan(p.lower_bound) and np.isnan(p.upper_bound):
                p.lower_bound = p.value * (1 - bounds_fraction)
                p.upper_bound = p.value * (1 + bounds_fraction)

        # print parameters
        pd.options.display.float_format = "{:.5g}".format
        df_parameters = SensitivityParameter.parameters_to_df(parameters)
        console.print(df_parameters)

        return parameters

    def sensitivity_groups(self) -> list[AnalysisGroup]:
        groups = [
            AnalysisGroup(
                uid="control",
                name="Control",
                changes={},
                color="dimgrey",
            ),
            AnalysisGroup(
                uid="mildRI",
                name="Mild renal impairment",
                changes={"KI__f_renal_function": 0.69},
                color="#66c2a4",
            ),
            AnalysisGroup(
                uid="modRI",
                name="Moderate renal impairment",
                changes={"KI__f_renal_function": 0.32},
                color="#2ca25f",
            ),
            AnalysisGroup(
                uid="sevRI",
                name="Severe renal impairment",
                changes={"KI__f_renal_function": 0.19},
                color="#006d2c",
            ),
            AnalysisGroup(
                uid="CPT A",
                name="Mild cirrhosis (CPT A)",
                changes={"f_cirrhosis": 0.399},
                color="#74a9cf",
            ),
            AnalysisGroup(
                uid="CPT B",
                name="Moderate cirrhosis (CPT B)",
                changes={"f_cirrhosis": 0.698},
                color="#2b8cbe",
            ),
            AnalysisGroup(
                uid="CPT C",
                name="Severe cirrhosis (CPT C)",
                changes={"f_cirrhosis": 0.813},
                color="#045a8d",
            )
        ]
        return groups


def local_sensitivity_analysis():
    """Local sensitivity analysis"""
    console.rule("LOSARTAN LOCAL SENSITIVITY ANALYSIS", style="blue bold", align="center")

    sensitivity_simulation = LosartanSensitivitySimulation.sensitivity_simulation()
    parameters = sensitivity_simulation.sensitivity_parameters()
    groups = sensitivity_simulation.sensitivity_groups()
    sa = LocalSensitivityAnalysis(
        sensitivity_simulation=sensitivity_simulation,
        parameters=parameters,
        groups=groups,
        results_path=RESULTS_PATH / "sensitivity",
        seed=1234,
        difference=0.01,
    )

    console.rule("Samples", style="white")
    sa.create_samples()

    console.rule("Results", style="white")
    sa.simulate_samples()
    console.print(sa.results)

    console.rule("Sensitivity", style="white")
    sa.calculate_sensitivity()
    console.print(sa.sensitivity)

    console.rule("Plotting", style="white")
    for kg, group in enumerate(sa.groups):
        sa.plot_sensitivity(
            group_id=group.uid,
            sensitivity_key="normalized",
            # title=f"{group.name}",
            cutoff=0.05,
            cluster_rows = False,
            cmap = "seismic",
            vcenter=0.0,
            vmin=-2.0,
            vmax=2.0,
            fig_path=sa.results_path / f"local_sensitivity_{kg:>02}_{group.uid}_{sa.difference}.png",
        )


def global_sensitivity_analysis():
    """Global sensitivity analysis"""

    console.rule("LOSARTAN GLOBAL SENSITIVITY ANALYSIS", style="blue bold", align="center")
    sensitivity_simulation = LosartanSensitivitySimulation.sensitivity_simulation()
    parameters = sensitivity_simulation.sensitivity_parameters()
    groups = sensitivity_simulation.sensitivity_groups()
    # only analysis on control group
    groups = [g for g in groups if g.uid == "control"]

    sa = SobolSensitivityAnalysis(
        sensitivity_simulation=sensitivity_simulation,
        parameters=parameters,
        groups=groups,
        results_path=RESULTS_PATH / "sensitivity",
        N=4096,
        seed=1234,
    )

    console.rule("Samples", style="white")
    sa.create_samples()
    console.print(sa.samples)

    console.rule("Results", style="white")
    cache_results: bool = True
    results_path = sa.results_path / f"sobol_sensitivity_N{sa.N}_results.pkl"
    if not cache_results or (cache_results and not results_path.exists()):
        sa.simulate_samples()
        with open(results_path, 'wb') as f:
            dill.dump(sa.results, f)
    else:
        with open(results_path, 'rb') as f:
            sa.results = dill.load(f)

    console.print(sa.results)

    console.rule("Sensitivity", style="white")
    cache_sensitivity: bool = True
    sensitivity_path = sa.results_path / f"sobol_sensitivity_N{sa.N}.pkl"
    if not cache_sensitivity or (cache_sensitivity and not sensitivity_path.exists()):
        sa.calculate_sensitivity()
        with open(sensitivity_path, 'wb') as f:
            dill.dump(sa.sensitivity, f)
    else:
        with open(sensitivity_path, 'rb') as f:
            sa.sensitivity = dill.load(f)

    console.print(sa.sensitivity)

    console.rule("Plotting", style="white")
    # Heatmaps
    for kg, group in enumerate(sa.groups):
        for key in ["ST", "S1"]:
            sa.plot_sensitivity(
                group_id=group.uid,
                sensitivity_key=key,
                # title=f"{key} {group.name}",
                cutoff=0.05,
                cluster_rows=False,
                cmap="viridis",
                vcenter=0.5,
                vmin=0.0,
                vmax=1.0,
                fig_path=sa.results_path / f"sobol_sensitivity_N{sa.N}_{kg:>02}_{group.uid}_{key}.png"
            )

        # Barplots
        sa.plot_sobol_indices(
            fig_path=sa.results_path / f"sobol_sensitivity_N{sa.N}_{kg:>02}_{group.uid}.png",
        )

def sampling_sensitivity_analysis():
    """Sampling sensitivity/uncertainty analysis"""

    console.rule("LOSARTAN SAMPLING SENSITIVITY ANALYSIS", style="blue bold", align="center")
    sensitivity_simulation = LosartanSensitivitySimulation.sensitivity_simulation()
    parameters = sensitivity_simulation.sensitivity_parameters()
    groups = sensitivity_simulation.sensitivity_groups()

    sa = SamplingSensitivityAnalysis(
        sensitivity_simulation=sensitivity_simulation,
        parameters=parameters,
        results_path=RESULTS_PATH / "sensitivity",
        N=1000,
        seed=1234,
        groups=groups,
    )

    console.rule("Samples", style="white")
    sa.create_samples()
    console.print(sa.samples)

    console.rule("Results", style="white")
    cache_results: bool = True
    results_path = sa.results_path / f"sampling_sensitivity_N{sa.N}_results.pkl"
    if not cache_results or (cache_results and not results_path.exists()):
        sa.simulate_samples()
        with open(results_path, 'wb') as f:
            dill.dump(sa.results, f)
    else:
        with open(results_path, 'rb') as f:
            sa.results = dill.load(f)

    console.print(sa.results)

    console.rule("Sensitivity", style="white")
    sa.calculate_sensitivity()
    console.print(sa.sensitivity)

    sa.df_sampling_sensitivity(
        df_path=sa.results_path / f"sampling_sensitivity_N{sa.N}_statistics.tsv"
    )

    console.rule("Plotting", style="white")
    sa.plot_sampling_sensitivity(
        fig_path=sa.results_path / f"sampling_sensitivity_N{sa.N}.png",
    )

def parameter_table():
    console.rule("LOSARTAN PARAMETER TABLE", style="blue bold", align="center")
    sensitivity_simulation = LosartanSensitivitySimulation.sensitivity_simulation()
    parameters = sensitivity_simulation.sensitivity_parameters()
    tex_path = RESULTS_PATH / "sensitivity" / "parameter_table.tex"
    df = SensitivityParameter.parameters_to_df(parameters)
    tex_str = df.to_latex(
        None, index=False, float_format="{:.3g}".format
    )
    tex_str = tex_str.replace("_", r"\_")

    with open(tex_path, 'w') as f:
        f.write(tex_str)


if __name__ == "__main__":
    from pkdb_models.models.losartan import MODEL_PATH, RESULTS_PATH
    from matplotlib import pyplot as plt

    # parameters
    parameter_table()

    # sensitivity analysis
    # local_sensitivity_analysis()
    # sampling_sensitivity_analysis()
    global_sensitivity_analysis()



