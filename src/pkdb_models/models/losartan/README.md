# Losartan

Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of the angiotensin II receptor blocker losartan and its active metabolite E3174 (EXP3174): the oral absorption of losartan (dissolution, uptake by OATP2B1 and efflux by ABCB1) against the colonic transit, its distribution in blood (blood concentrations with the blood to plasma ratios of Lo1995, the plasma is a readout) and in the tissues, the hepatic first pass and metabolism of losartan to E3174 (CYP2C9) and to the other, inactive metabolites (half of them by CYP2C9), the hepatic uptake, export and biliary excretion of E3174 and of the other metabolites, and the renal excretion of losartan and its metabolites relative to the glomerular filtration rate. E3174 blocks the AT1 receptor: the AT1 signal (angiotensin II on the unblocked receptors) drives the aldosterone secretion, inhibits the renin secretion and sets the blood pressure, so that losartan raises renin and angiotensin and lowers aldosterone and the blood pressure. The model describes the CYP2C9 genotypes (*2, *3, *13), the ABCB1 2677/3435 genotypes, cirrhosis, renal impairment and the interactions with drugs which inhibit CYP2C9 (bucolome, phenytoin, silymarin, red ginseng).

The published version of the model is [losartan-model](https://github.com/matthiaskoenig/losartan-model) ([doi:10.5281/zenodo.14761114](https://doi.org/10.5281/zenodo.14761114)).

## Structure

- `__init__.py`: the paths of the model and its model object `MODEL` with its release metadata and its adopted fits, which the commands below run (pkpdsim).
- `models/`: the tissue models (`model_intestine.py`, `model_liver.py`, `model_kidney.py`), the renin-angiotensin-aldosterone system (`model_raas.py`) and the whole-body model (`model_body.py`). `models/factory.py` creates the SBML models, the flattened whole-body model which is simulated, the ODEs and the COMBINE archive in `models/output/`; the fitted values come from the adopted fits (`models/fits.py`, `models/output/fits/`).
- `experiments/`: the simulation experiments of the 26 clinical studies (`studies/`), the dose dependency, the genetic variants and the RAAS (`misc/`) and the parameter scans of cirrhosis, renal function, CYP2C9 and ABCB1 activity and dose (`scans/`). The fit mappings are described with the vocabulary of `pkpdsim.metadata`.
- `fitting/`: the fit problems of the pharmacokinetics (`PK`) and of the pharmacodynamics on the adopted pharmacokinetics (`PD`), their parameters and the selection of the training, validation and outlier data (`mapping_collections.py`).
- `sensitivity/`: the sensitivity analysis of the pharmacokinetics of losartan and its metabolites, the AT1 blockade, the RAAS and the blood pressure.
- `losartan_pk.py`: the pharmacokinetic and pharmacodynamic parameters of the scans, the non-compartmental analysis of `pkpdutils`.
- `tolerances.py`: the absolute tolerances of the integrator per state.
- `simulations.py`: the simulation experiments and their groups, `MODEL` runs them.
- `data/`: the data of the studies, a subset of the PK-DB studies (`data/collect_data.py`).
- `latex/`: the tables of the studies and of the models of the publication.
- `TODO.md`: the open issues of the model.

## Commands

```bash
# models, simulations
run_simulation pkdb_models.models.losartan --action factory
run_simulation pkdb_models.models.losartan --action list_experiments
run_simulation pkdb_models.models.losartan --action simulate --experiments all
run_simulation pkdb_models.models.losartan --action simulate --experiments all --uncertainty

# fits, adoption of a fit, reports of stored parameters, identifiability
run_fit pkdb_models.models.losartan --subset=PK --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.losartan --subset=PD --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.losartan --adopt PK_<fit id>
run_report pkdb_models.models.losartan --subset=PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.losartan --subset=PK results/fit/PK_<fit id>/parameters.json

# sensitivity analysis (local by default)
python -m pkdb_models.models.losartan.sensitivity.sensitivity_analysis
```

The results are written into `results/` of the model.

## Adopted fits and uncertainty

The values of the fitted parameters come from the adopted fits of the model,
not from numbers in the code: the pharmacokinetics from
`PK_20261019_113633__d9465`, the pharmacodynamics (fitted on the adopted
pharmacokinetics) from `PD_20261019_134029__a3a50`.
`run_fit pkdb_models.models.losartan --adopt <fit id>` copies the record of a
finished fit (its parameters, problem, profile likelihood and Fisher
information) into `models/output/fits/<fit id>/`, which is tracked; the ids of
the adopted fits are listed in `fits` of the model object (`__init__.py`). The
factory takes the fitted values and their uncertainty from the records
(`models/fits.py`) and writes the uncertainty into the SBML
(`distrib:uncertainty`). `run_simulation ... --uncertainty` runs the
experiments over draws of these distributions and draws the simulations as the
median with a band of quantiles into `results/simulation_uncertainty`. A
parameter which the profile likelihood does not identify is drawn uniformly
between the bounds of its fit (the tissue distribution of E3174 and the renal
excretion of the other metabolites), which widens the bands of E3174. A refit
is: fit, adopt, replace the id in `fits`, rebuild the model. See the
documentation of the repository on fitting (`docs/fitting.md`).
