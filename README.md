[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14761114.svg)](https://doi.org/10.5281/zenodo.14761114)
[![Release](https://img.shields.io/github/v/release/matthiaskoenig/losartan-model)](https://github.com/matthiaskoenig/losartan-model/releases)
[![Documentation](https://github.com/matthiaskoenig/losartan-model/actions/workflows/docs.yml/badge.svg)](https://matthiaskoenig.github.io/losartan-model/)
[![CI](https://github.com/matthiaskoenig/losartan-model/actions/workflows/ci.yml/badge.svg)](https://github.com/matthiaskoenig/losartan-model/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/docker/v/matthiaskoenig/losartan?sort=semver&label=docker)](https://hub.docker.com/r/matthiaskoenig/losartan)
[![Python](https://img.shields.io/badge/python-3.14-blue)](https://www.python.org)
[![Code license: MIT](https://img.shields.io/badge/code%20license-MIT-blue)](https://opensource.org/license/MIT)
[![Model license: CC BY-SA 4.0](https://img.shields.io/badge/model%20license-CC%20BY--SA%204.0-lightgrey)](https://creativecommons.org/licenses/by-sa/4.0/)

# Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of losartan

This repository provides the losartan physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model: the oral absorption of losartan, its distribution, its hepatic metabolism to the active metabolite E3174 (CYP2C9) and to inactive metabolites, the renal and biliary excretion of losartan and its metabolites, and the blockade of the AT1 receptor by E3174 with its effects on renin, angiotensin, aldosterone and the blood pressure.

Documentation: <https://matthiaskoenig.github.io/losartan-model/>

## Quick start

```bash
git clone https://github.com/matthiaskoenig/losartan-model.git && cd losartan-model
uv sync
uv run run_simulation pkdb_models.models.losartan --action all -r results
```

## How to cite

> Tensil, E., Myshkina, M. & König, M. (2026). *Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of losartan*. Zenodo. https://doi.org/10.5281/zenodo.14761114

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

&copy; 2024-2026 Ennie Tensil, Mariia Myshkina and Matthias König
