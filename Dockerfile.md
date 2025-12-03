# Dockerfile for Glimepiride Model

## Build image
To build the latest image use:
```bash
docker build -f Dockerfile -t matthiaskoenig/losartan:0.9.6 -t matthiaskoenig/losartan:latest .
```

## Push images
The image is pushed to dockerhub: [Docker Hub – Dapagliflozin](https://hub.docker.com/repository/docker/matthiaskoenig/losartan/general)

```
docker login
docker push --all-tags matthiaskoenig/losartan
```

## Run container
To use the latest container version interactively use:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/losartan:latest /bin/bash
```

To use a specific container version provide the version tag:
```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/losartan:0.9.6 /bin/bash
```

## Run simulations
Run the complete analysis:
```bash
uv run run_losartan -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```
## Funding

Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within LiSyM by grant number 031L0054 and ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B). 
Mariia Myshkina was supported by the Federal Ministry of Education and Research (BMBF, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA).

© 2024-2025 Ennie Tensil, Mariia Myshkina, Michelle Elias & Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)