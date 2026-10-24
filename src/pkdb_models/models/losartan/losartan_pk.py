"""Pharmacokinetic and pharmacodynamic parameters of the losartan scans.

The pharmacokinetic parameters are the non-compartmental analysis of
`pkpdutils`, named like the variables of its `NCAResult` (e.g. `auc_last`,
`auc_inf_obs`, `cmax`, `thalf`, `lambda_z`, `cl_f`, `vz_f`), with the unit of
every parameter in `<name>_unit`, for losartan and its metabolites E3174 and
the other metabolites in plasma (`[Cveplasma_<sid>]`). Losartan is analysed
with its dose (oral or intravenous), the metabolites without a dose (`cl_f`
and `vz_f` are not a number). The concentrations below `LLOQ_FRACTION` of the
maximum are below the limit of quantification, so that the noise of the
integrator is not the terminal phase. The renal and fecal clearances are the
amounts excreted until the end of the simulation over `auc_last`. The
pharmacodynamic parameters are the minimum and the maximum of the RAAS and the
blood pressure.

A scan has one dimension, the scanned parameter; a run over the draws of the
fitted parameters (`run_simulation pkdb_models.models.losartan --uncertainty`)
adds the dimension of the draws. The parameters are calculated for every point
of every dimension, a row per point with its position along each dimension
(the columns `point` and, with draws, `draw`), see `summarize`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pkpdsim.pk import nca_frame, scan_nca
from pkpdsim.scan import end_values, scan_dimensions, start_values
from pkpdutils import Route, ureg
from sbmlsim import Q

if TYPE_CHECKING:
    from sbmlsim.result import ScanResult

    from pkdb_models.models.losartan.experiments.base_experiment import (
        LosartanSimulationExperiment,
    )

#: the doses of losartan of the scans by their route
DOSES: dict[str, Route] = {
    "PODOSE_los": Route.ORAL,
    "IVDOSE_los": Route.IV_BOLUS,
}
#: the substances by id of the concentration in plasma (`[Cveplasma_<id>]`)
SUBSTANCES: dict[str, str] = {
    "los": "losartan",
    "e3174": "E3174",
    "l158": "other metabolites",
}
#: the readouts whose minimum and maximum are the pharmacodynamic parameters
EFFECTS: list[str] = [
    "[ang1]",
    "[ang2]",
    "[ren]",
    "[ald]",
    "SBP",
    "DBP",
    "MAP",
]
#: the limit of quantification relative to the maximal concentration
LLOQ_FRACTION = 1e-4


def _names(res: ScanResult, across: str | None) -> dict[str, str]:
    """Get the columns of the dimensions of a scan, `point` and `draw`.

    The scan has one dimension besides the draws `across`.

    Raises:
        ValueError: if the scan has not one dimension besides the draws.
    """
    dims = scan_dimensions(res, across)
    if len(dims) != 1:
        raise ValueError(
            f"A scan of one dimension is required, not of {dims} (draws: {across})."
        )
    return {dims[0]: "point", **({across: "draw"} if across is not None else {})}


def _doses(
    experiment: LosartanSimulationExperiment, res: ScanResult
) -> tuple[np.ndarray, np.ndarray]:
    """Get the dose of losartan of every simulation of a scan and its route.

    The dose is the first of `DOSES` at the start of the simulation, in the
    substance unit of the concentration, which `cl_f` and `vz_f` require.

    Returns:
        The amounts [mmole] and the routes over the dimensions of the scan,
        `NaN` and `None` for a simulation without a dose.
    """
    shape = tuple(res.ds.sizes[dim] for dim in res.dims)
    amounts = np.full(shape, np.nan)
    routes = np.full(shape, None, dtype=object)
    # in reverse, so that the first dose of `DOSES` wins
    for key, route in reversed(DOSES.items()):
        mass = Q(start_values(res, key).values, res.units[key])
        given = mass.magnitude > 0
        amounts = np.where(
            given, (mass / experiment.Mr.los).to("mmole").magnitude, amounts
        )
        routes[given] = route
    return amounts, routes


def pk_parameters(
    experiment: LosartanSimulationExperiment,
    res: ScanResult,
    sid: str,
) -> pd.DataFrame:
    """Calculate the pharmacokinetic parameters of a substance in every simulation of a scan.

    The non-compartmental analysis of all simulations with a dose of one route
    at once (`scan_nca`), and of all simulations without a dose (the
    metabolites) at once. The renal and fecal clearances are the amounts
    excreted in the urine and the feces until the end of the simulation over
    `auc_last`, the total clearance by excretion is their sum.

    Args:
        experiment: experiment of the scan, for the molecular weight.
        res: result of the scan.
        sid: id of the substance, a key of `SUBSTANCES`.

    Returns:
        One row per simulation of the scan (in the order of the dimensions of
        the scan), with its position along every dimension (a column named
        after the dimension), `substance` (its id), the parameters and their
        units (`<name>_unit`).
    """
    concentration = f"[Cveplasma_{sid}]"
    shape = tuple(res.ds.sizes[dim] for dim in res.dims)
    amounts, routes = (
        _doses(experiment, res)
        if sid == "los"
        else (np.full(shape, np.nan), np.full(shape, None, dtype=object))
    )
    excretions = {"cl_renal": f"Aurine_{sid}", "cl_fecal": f"Afeces_{sid}"}
    # an analysis per route (in the order of the first simulation of a route),
    # of which the simulations of the route are taken
    labels = routes.ravel()
    frames = []
    for route in dict.fromkeys(labels):
        result = scan_nca(
            res,
            concentration,
            dose=None if route is None else Q(amounts, "mmole"),
            route=route,
            substance=SUBSTANCES[sid],
            lloq_fraction=LLOQ_FRACTION,
        )
        df = nca_frame(result).drop(columns="nca_flags")[labels == route]
        position = tuple(df[dim].to_numpy() for dim in res.dims)
        auc_last = ureg.Quantity(df["auc_last"].to_numpy(), df["auc_last_unit"].iloc[0])
        excreted = {
            amount: ureg.Quantity(
                end_values(res, amount).values[position], str(res.units[amount])
            )
            for amount in excretions.values()
        }
        clearances = {
            cl_key: (excreted[amount] / auc_last).to("ml/min")
            for cl_key, amount in excretions.items()
        }
        clearances["cl_total"] = clearances["cl_renal"] + clearances["cl_fecal"]
        for key, value in {**excreted, **clearances}.items():
            df[key] = value.magnitude
            df[f"{key}_unit"] = str(value.units)
        frames.append(df)
    df = pd.concat(frames).sort_index()
    df.insert(len(res.dims), "substance", sid)
    return df


def calculate_losartan_pk(
    experiment: LosartanSimulationExperiment,
    res: ScanResult,
    across: str | None = None,
) -> pd.DataFrame:
    """Calculate the pharmacokinetic parameters of every simulation of a scan.

    Args:
        experiment: experiment of the scan.
        res: result of the scan.
        across: the dimension of the draws of the fitted parameters, `None`
            for a scan without draws.

    Returns:
        One row per substance and simulation of the scan, in the order of the
        simulations for every substance, with its position (`point`, and
        `draw` with draws).
    """
    names = _names(res, across)
    return pd.concat(
        [
            pk_parameters(experiment, res, sid).rename(columns=names)
            for sid in SUBSTANCES
        ],
        ignore_index=True,
    )


def calculate_losartan_pd(
    res: ScanResult, across: str | None = None
) -> dict[str, pd.DataFrame]:
    """Calculate the pharmacodynamic parameters of every simulation of a scan.

    The minimum and the maximum of the timecourse of every simulation (the
    state after a change at its time, the padding of a ragged result left
    out), all simulations at once.

    Args:
        res: result of the scan.
        across: the dimension of the draws of the fitted parameters, `None`
            for a scan without draws.

    Returns:
        By readout, one row per simulation with the minimum and the maximum
        and its position (`point`, and `draw` with draws).
    """
    names = _names(res, across)
    shape = tuple(res.ds.sizes[dim] for dim in res.dims)
    position = {
        names[dim]: index.ravel()
        for dim, index in zip(res.dims, np.indices(shape), strict=True)
    }
    frames: dict[str, pd.DataFrame] = {}
    for sid in EFFECTS:
        timecourses = res.to_timecourses(sid)
        values = timecourses.values.reshape(-1, timecourses.n_time)
        frames[sid] = pd.DataFrame(
            {
                **position,
                "sid": sid,
                "min": np.nanmin(values, axis=-1),
                "max": np.nanmax(values, axis=-1),
                "unit": str(res.units[sid]),
            }
        )
    return frames


def summarize(
    df: pd.DataFrame, key: str, quantiles: tuple[float, float]
) -> pd.DataFrame:
    """Get the median and the quantiles of a parameter over the draws of every point.

    Args:
        df: the parameters with the columns `point` and `draw`.
        key: the parameter.
        quantiles: the lower and the upper quantile.

    Returns:
        One row per point (the index), the columns `lower`, `median`, `upper`.
    """
    low, high = quantiles
    grouped = df.groupby("point")[key]
    return pd.DataFrame(
        {
            "lower": grouped.quantile(low),
            "median": grouped.median(),
            "upper": grouped.quantile(high),
        }
    ).astype(np.float64)
