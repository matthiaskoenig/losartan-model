[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14761114.svg)](https://doi.org/10.5281/zenodo.14761114)

# losartan model
This repository provides the losartan physiologically based pharmacokinetics (PBPK) model.

The model is distributed as [SBML](http://sbml.org) available from [`losartan_body_flat.xml`](./models/losartan_body_flat.xml) with 
corresponding SBML4humans model report at [https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_body_flat.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_body_flat.xml) and equations from [`losartan_body_flat.md`](./models/losartan_body_flat.md).

The COMBINE archive is available from [`losartan_model.omex`](./losartan_model.omex).

![PK/PD model overview](./figures/losartan_model.png)


### Comp submodels
The liver submodel is available from [`losartan_liver.xml`](./models/losartan_liver.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_liver.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_liver.xml) and equations from [`losartan_liver.md`](./models/losartan_liver.md).

The kidney submodel is available from [`losartan_kidney.xml`](./models/losartan_kidney.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_kidney.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_kidney.xml) and equations from [`losartan_kidney.md`](./models/losartan_kidney.md).

The intestine submodel is available from [`losartan_intestine.xml`](./models/losartan_intestine.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_intestine.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_intestine.xml) and equations from [`losartan_intestine.md`](./models/losartan_intestine.md).

The whole-body submodel is available from [`losartan_body.xml`](./models/losartan_body.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_body.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/losartan-model/main/models/losartan_body.xml) and equations from [`losartan_body.md`](./models/losartan_body.md).

## How to cite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14761114.svg)](https://doi.org/10.5281/zenodo.14761114)

> Tensil, E., & König, M. (2025).
> *Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of losartan.*   
> Zenodo. [https://doi.org/10.5281/zenodo.14761114](https://doi.org/10.5281/zenodo.14761114)

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

## Run simulations
### python
Clone the repository 
```bash
git clone https://github.com/matthiaskoenig/losartan-model.git
cd losartan-model
```

#### uv
Setup environment with uv (https://docs.astral.sh/uv/getting-started/installation/)
```bash
uv venv
uv sync
```
Run the complete analysis:
```bash
uv run run_losartan -a all -r results
```

#### pip
If you use pip install the package via
```bash
pip install -e .
```
Run the complete analysis in the environment via:
```bash
run run_losartan -a all -r results
```

### docker
Simulations can also be run within a docker container:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/losartan:latest /bin/bash
```

Run the complete analysis:
```bash
uv run run_losartan -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"

Funding
=======
Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within LiSyM by grant number 031L0054 and ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B). 

© 2024-2025 Ennie Tensil & Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)
