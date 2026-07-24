# Installation

## Requirements

- Linux with [Apptainer/Singularity](https://apptainer.org/) (for the tools container)
- [conda/mamba](https://github.com/conda-forge/miniforge)
- [Snakemake](https://snakemake.readthedocs.io/) >= 7
- GNU coreutils, `awk`, `bash`

The heavy binaries (art_modern, samtools, rsem-tbam2gbam, seqkit) run inside a
Singularity image so you do not have to build them into your host environment.
The R statistical steps run from a conda environment.

## 1. R environment

```bash
mamba env create -f workflow/envs/r.yaml
conda activate scARTist-r
# scDesign2 is not on bioconda:
Rscript -e 'remotes::install_github("JSB-UCLA/scDesign2")'
```

Alternatively let Snakemake manage it per-rule with `--use-conda`.

## 2. Tools container

```bash
apptainer build resources/art.sif containers/art_modern.def
```

Pin the art_modern commit/tag in `containers/art_modern.def` before release.
If you already have a working `art.sif`, just point `container.sif` in
`config/config.yaml` at it and skip this step.

## 3. Snakemake

```bash
pip install snakemake
snakemake --version
```

## Verify

```bash
snakemake -n            # dry run: prints the planned jobs
```
