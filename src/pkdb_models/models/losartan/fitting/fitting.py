"""Losartan parameter fitting.

The fit problems of the losartan model: which fit mapping collections enter a
fit and which parameters are adjusted, see `FIT_DEFINITIONS`. The commands
`run_fit`, `run_report` and `run_identifiability` of pkpdsim run them with the
package of the model:

    run_fit pkdb_models.models.losartan --subset=PK --runs=3 --cores=3 --seed=1234 --method=DE
    run_fit pkdb_models.models.losartan --subset=PD --runs=3 --cores=3 --seed=1234 --method=DE

The pharmacokinetics are fitted first, the pharmacodynamics on the adopted
values of the pharmacokinetics (`MODEL.fits`).

The results, figures and the HTML report are written into the fit folder of the
model (`MODEL.fit_dir`), unless `--output_dir` is given. The identifiability of
parameters which were fitted earlier is analysed without fitting again with the
`run_identifiability` command, a report of stored parameters is created with
the `run_report` command:

    run_identifiability pkdb_models.models.losartan --subset=PK results/fit/PK_<fit id>/parameters.json
    run_report pkdb_models.models.losartan --subset=PK results/fit/PK_<fit id>/parameters.json
"""

from collections.abc import Callable

from pkpdsim.fit import default_fit_settings
from sbmlsim.fit import FitMappingCollection, FitParameter
from sbmlsim.fit.cli import FitDefinition

from pkdb_models.models.losartan import (
    DATA_PATHS,
    LOSARTAN_PATH,
)
from pkdb_models.models.losartan.fitting.mapping_collections import (
    f_collections_pd,
    f_collections_pk,
)
from pkdb_models.models.losartan.fitting.parameters import (
    parameters_pd,
    parameters_pk,
)
from pkdb_models.models.losartan.tolerances import absolute_tolerance

#: settings of the fits, stored with the result and read back by the report
#: the species and the RAAS get tolerances of their own scale, see `tolerances`
FIT_SETTINGS = default_fit_settings(
    absolute_tolerance=absolute_tolerance(1e-6, species=1e-10, raas=1e-14)
)


def _definition(
    mapping_collections: Callable[[], dict[str, list[FitMappingCollection]]],
    parameters: list[FitParameter],
) -> FitDefinition:
    """Get a fit problem of the losartan model with the settings of the model."""
    return FitDefinition(
        mapping_collections=mapping_collections,
        parameters=parameters,
        base_path=LOSARTAN_PATH,
        data_path=DATA_PATHS,
        settings=FIT_SETTINGS,
    )


#: the fit problems of the losartan model, the keys are the `--subset` argument
FIT_DEFINITIONS: dict[str, FitDefinition] = {
    # the pharmacokinetics of losartan, E3174 and the other metabolites
    "PK": _definition(f_collections_pk, parameters_pk),
    # the pharmacodynamics (RAAS, blood pressure) on the adopted
    # pharmacokinetics
    "PD": _definition(f_collections_pd, parameters_pd),
}
