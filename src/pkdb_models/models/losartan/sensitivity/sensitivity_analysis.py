"""Sensitivity analysis of the losartan model.

The outputs are the pharmacokinetic parameters of a single oral dose of 10 mg
(the non-compartmental analysis of `pkpdutils`) of losartan and its metabolites
E3174 and the other metabolites and the extrema of the AT1 blockade, the RAAS
and the blood pressure. The parameters are the parameters of the model with
bounds of +-15 % around their value, the fitted parameters (the adopted fits,
`models/fits.py`) among them. The analyses are run as a script, the local
analysis by default:

    python -m pkdb_models.models.losartan.sensitivity.sensitivity_analysis
    python -m pkdb_models.models.losartan.sensitivity.sensitivity_analysis local sampling

The analyses run on `pkpdsim.sensitivity.SensitivityRunner`, which writes the
table of the parameters into `results/sensitivity` of the model and the tables
and figures of an analysis into `results/sensitivity/<analysis>`.
An analysis does not vary a parameter which its simulation or one of its groups
sets (`pkpdsim.sensitivity.free_parameters`): every analysis leaves out
`KI__f_renal_function`, which the simulation sets (`f_cirrhosis`, 0 in the
model, is no parameter of the analyses). The groups are the classes of renal
impairment and cirrhosis of the simulation experiments (`renal_map`,
`cirrhosis_map` of `LosartanSimulationExperiment`). The environment variable
`PKPDSIM_SENSITIVITY_N` replaces the numbers of samples of the global analyses
for a quick run, see `pkpdsim.sensitivity.sample_size`.
"""

from collections.abc import Mapping
from typing import Any, NamedTuple

import numpy as np
from pkpdsim.sensitivity import (
    LocalAnalysis,
    ParameterType,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityParameter,
    SensitivityRunner,
    SobolAnalysis,
    parameters_from_sbml,
    set_bounds,
    set_relative_bounds,
)
from sbmlsim import Q
from sbmlsim.simulation import PK, Custom, Formula, ObservableKind, Simulation

from pkdb_models.models.losartan import MODEL, MODEL_PATH
from pkdb_models.models.losartan.experiments.base_experiment import (
    LosartanSimulationExperiment,
)
from pkdb_models.models.losartan.fitting.parameters import (
    parameters_pd,
    parameters_pk,
)

#: oral dose of losartan in the analysis [mg]
DOSE_LOS = 10
#: fraction of the bounds of a parameter relative to its value
BOUNDS_FRACTION = 0.15

#: outputs of the non-compartmental analysis of pkpdutils -> unit of the output
PK_OUTPUTS: dict[str, str] = {
    "auc_inf_obs": "mM*min",
    "cmax": "mM",
    "thalf": "min",
    "vz_f": "l",
    "cl_f": "ml/min",
    "lambda_z": "1/min",
}
#: the outputs of the metabolites, which have no dose of their own
PK_OUTPUTS_METABOLITES: list[str] = ["auc_inf_obs", "cmax", "thalf"]
#: limit of quantification relative to the maximal concentration: the
#: concentrations of the 20 days which are smaller are the noise of the
#: integrator, they would be the terminal phase
LLOQ_FRACTION = 1e-8
#: the outputs whose unit has the amount of the dose: the dose is a mass
#: (`PODOSE_los` in mg), so they are divided by the molecular weight (its value
#: at the start, a value per simulation like the parameter of the analysis)
PK_PER_DOSE_MASS: frozenset[str] = frozenset({"vz_f", "cl_f"})
#: the labels of the pharmacokinetic parameters in the names of the outputs
PK_LABELS: dict[str, str] = {
    "auc_inf_obs": "AUC∞",
    "cmax": "Cmax",
    "thalf": "Half-life",
    "vz_f": "Vz/F",
    "cl_f": "CL/F",
    "lambda_z": "λz",
}
#: the pharmacodynamic readouts and their extremum: the blood pressure as its
#: change from the baseline, whose fall is small against its value
PD_OUTPUTS: dict[str, str] = {
    "fe_e3174": "max",
    "[ang1]": "max",
    "[ang2]": "max",
    "[ren]": "max",
    "[ald]": "min",
    "SBP_change": "min",
    "DBP_change": "min",
}
#: the names of the pharmacodynamic readouts
PD_NAMES: dict[str, str] = {
    "fe_e3174": "AT1 blockade",
    "[ang1]": "angiotensin I",
    "[ang2]": "angiotensin II",
    "[ren]": "renin",
    "[ald]": "aldosterone",
    "SBP_change": "change SBP",
    "DBP_change": "change DBP",
}
#: the units of the pharmacodynamic readouts
PD_UNITS: dict[str, str] = {
    "fe_e3174": "dimensionless",
    "[ang1]": "mM",
    "[ang2]": "mM",
    "[ren]": "mM",
    "[ald]": "mM",
    "SBP_change": "mmHg",
    "DBP_change": "mmHg",
}


def quantified(time: np.ndarray, values: dict[str, Any]) -> np.ndarray:
    """Get a concentration above its limit of quantification, missing below it.

    The limit is `LLOQ_FRACTION` of the maximum of the timecourse, as the
    limit of quantification of the old analysis (`Timecourse.lloq`), whose
    values below it the non-compartmental analysis leaves out.

    Args:
        time: the time points of the simulation.
        values: the concentration, the one symbol of the observable.

    Returns:
        The concentration, `NaN` below the limit.
    """
    (concentration,) = values.values()
    c = np.asarray(concentration, dtype=float)
    return np.where(c < LLOQ_FRACTION * np.max(c), np.nan, c)


class Substance(NamedTuple):
    """A substance of the pharmacokinetic outputs, its concentration in plasma."""

    #: the concentration, which the non-compartmental analysis reads
    selection: str
    #: the name of the substance in the names of the outputs
    name: str
    #: the id of the non-compartmental analysis, no id of the model
    pk_id: str
    #: the outputs of the analysis -> unit
    outputs: dict[str, str]
    #: the dose, of losartan only
    dose: str | None

    @property
    def quantified_id(self) -> str:
        """Get the id of the concentration above the limit of quantification."""
        return f"{self.selection.strip('[]')}_quantified"


#: the substances of the pharmacokinetic outputs
SUBSTANCES: list[Substance] = [
    Substance("[Cveplasma_los]", "LOS", "pk_los", PK_OUTPUTS, "PODOSE_los"),
    *(
        Substance(
            f"[Cveplasma_{sid}]",
            name,
            f"pk_{sid}",
            {key: PK_OUTPUTS[key] for key in PK_OUTPUTS_METABOLITES},
            None,
        )
        for sid, name in [("e3174", "E3174"), ("l158", "other metabolites")]
    ),
]

#: the outputs: the id of a pharmacokinetic output names the concentration with
#: its brackets (`[Cveplasma_los]_cmax`, `[ang1]_max`), as the tables of sbmlsim
#: 0.8.5, its observable has the id without them (`Cveplasma_los_cmax`, `ang1_max`)
OUTPUTS: list[SensitivityOutput] = [
    *(
        SensitivityOutput(
            f"{s.selection}_{key}",
            f"{s.name} {PK_LABELS[key]}",
            f"{s.selection.strip('[]')}_{key}",
        )
        for s in SUBSTANCES
        for key in s.outputs
    ),
    *(
        SensitivityOutput(
            f"{sid}_{extremum}",
            f"{PD_NAMES[sid]} {extremum}",
            f"{sid.strip('[]')}_{extremum}",
        )
        for sid, extremum in PD_OUTPUTS.items()
    ),
]


def _pk_formula(substance: Substance, key: str) -> str:
    """Get the formula of a pharmacokinetic output of a substance."""
    if key in PK_PER_DOSE_MASS:
        return f"{substance.pk_id}.{key} / at(Mr_los, 0)"
    return f"{substance.pk_id}.{key}"


#: the observables of a simulation, the outputs with their units
OBSERVABLES = [
    *(
        Custom(
            s.quantified_id,
            quantified,
            "mM",
            symbols=[s.selection],
            kind=ObservableKind.TIMECOURSE,
        )
        for s in SUBSTANCES
    ),
    *(
        PK(
            s.pk_id,
            s.quantified_id,
            dose=s.dose,
            route="oral" if s.dose else None,
            parameters=list(s.outputs),
        )
        for s in SUBSTANCES
    ),
    *(
        Formula(f"{s.selection.strip('[]')}_{key}", _pk_formula(s, key), unit=unit)
        for s in SUBSTANCES
        for key, unit in s.outputs.items()
    ),
    *(
        Formula(
            f"{sid.strip('[]')}_{extremum}", f"{extremum}({sid})", unit=PD_UNITS[sid]
        )
        for sid, extremum in PD_OUTPUTS.items()
    ),
]


def sensitivity_simulation(changes: Mapping[str, float]) -> Simulation:
    """Get the simulation of a group: a single oral dose of losartan.

    The simulation sets the healthy liver and kidney, which a group changes.
    """
    return Simulation(
        end=20 * 24 * 60,
        time_unit="min",
        preinit_changes={
            "PODOSE_los": Q(DOSE_LOS, "mg"),
            "f_cirrhosis": Q(0, "dimensionless"),
            "KI__f_renal_function": Q(1.0, "dimensionless"),
            **changes,
        },
    )


#: the groups the analyses are run for: the control and the stages of renal
#: impairment and cirrhosis
SENSITIVITY_GROUPS: list[SensitivityGroup] = [
    SensitivityGroup("control", "Control", sensitivity_simulation({}), "dimgrey"),
    *(
        SensitivityGroup(
            uid,
            name,
            sensitivity_simulation(
                {"KI__f_renal_function": LosartanSimulationExperiment.renal_map[name]}
            ),
            LosartanSimulationExperiment.renal_colors[name],
        )
        for uid, name in [
            ("mildRI", "Mild renal impairment"),
            ("modRI", "Moderate renal impairment"),
            ("sevRI", "Severe renal impairment"),
        ]
    ),
    *(
        SensitivityGroup(
            uid,
            f"{name} ({uid})",
            sensitivity_simulation(
                {"f_cirrhosis": LosartanSimulationExperiment.cirrhosis_map[name]}
            ),
            LosartanSimulationExperiment.cirrhosis_colors[name],
        )
        for uid, name in [
            ("CPT A", "Mild cirrhosis"),
            ("CPT B", "Moderate cirrhosis"),
            ("CPT C", "Severe cirrhosis"),
        ]
    ),
]


def sensitivity_parameters() -> list[SensitivityParameter]:
    """Get the parameters of the analyses with their bounds and references."""
    parameters = parameters_from_sbml(
        MODEL_PATH,
        exclude_ids={
            # conversion factors
            "conversion_min_per_day",
            # molecular weights
            "Mr_los",
            "Mr_e3174",
            "Mr_l158",
            # unchangeable values
            "FQlu",
            "FVhv",
            "FVpo",
            # fast transports
            "LI__L158EX_k",
            "LI__LOSIM_k",
            # the RAAS turnover (fast against the pharmacokinetics)
            "ren_t12",
            "ald_t12",
            "ang1_t12",
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
            # the GFR is calibrated to the renal blood flow
            "f_GFR_Qki",
            # reference values
            "BW_ref",
            "anggen_ref",
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

    # the fitted parameters, their bounds are set relative to the value below
    set_bounds(
        parameters,
        [
            (fp.pid, np.nan, np.nan, ParameterType.FIT)
            for fp in [*parameters_pk, *parameters_pd]
        ],
    )
    # the activities of the transporters and enzymes and the renal function
    set_bounds(
        parameters,
        [
            (uid, 1 - BOUNDS_FRACTION, 1 + BOUNDS_FRACTION, ParameterType.SCALING)
            for uid in [
                "GU__f_OATP2B1",
                "GU__f_abcb1",
                "LI__f_cyp2c9",
                "KI__f_renal_function",
            ]
        ],
    )

    # the references of the values from the literature
    references = {
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
        "COBW": r"\cite{Cattermole2017, Patel2021, Collis2001}",
        "BP_los": r"\cite{Lo1995}",
        "BP_e3174": r"\cite{Lo1995}",
        "KI__fu_los": r"\cite{Lo1995}",
        "KI__fu_e3174": r"\cite{Lo1995}",
        "GU__f_abcb1": r"\cite{Siegmund2002, Hoffmeyer2000}",
        "KI__f_renal_function": r"\cite{Stevens2024}",
        "LI__f_cyp2c9": r"\cite{Kusama2009, Wang2014, Maekawa2009}",
    }
    by_uid = {p.uid: p for p in parameters}
    for uid, reference in references.items():
        p = by_uid[uid]
        p.reference = reference
        if p.type == ParameterType.NA:
            p.type = ParameterType.DATA

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, BOUNDS_FRACTION)
    return parameters


#: the analyses of the script, with their groups and numbers of samples
SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=SENSITIVITY_GROUPS, delta=0.01),
        "sampling": SamplingAnalysis(groups=SENSITIVITY_GROUPS, n=1000),
        "sobol": SobolAnalysis(groups=SENSITIVITY_GROUPS[:1], n=4096),
    },
    default=("local",),
)


if __name__ == "__main__":
    SENSITIVITY.main()
