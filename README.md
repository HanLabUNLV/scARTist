# scARTist

**sc** + **ART**  + **ist** = single-cell ART isoform simulation tool.

Single-cell isoform-resolved read simulation for benchmarking differential
splicing / differential transcript usage (DTU) methods.

`scARTist` re-simulates a **real single-cell isoform count matrix** with
[scDesign2](https://github.com/JSB-UCLA/scDesign2) (copula model, learns
gene-gene correlation, overdispersion and dropout), injects **ground-truth
differential expression and differential transcript usage**, renders the
per-cell counts into realistic paired-end short reads with
[art_modern](https://github.com/YU-Zhejian/art_modern), and maps the
transcriptome alignments to genome coordinates with
[RSEM](https://github.com/deweylab/RSEM). The output is a set of per-cell
genome-space BAMs (and optional FastQs) with a known isoform-usage ground truth.

> Status: pre-release scaffold ported from research scripts. The workflow
> encodes the full DAG but has not yet been executed end-to-end; smoke-test on
> the bundled example before a full run.

## Pipeline

```
real B-cell isoform count matrix
   |  scDesign2 fit + inject DE/DTU ground truth        [Stage 1: once]
   v
per-cell isoform counts, swept over depth x cell-count  [Stage 2]
   |  attach transcript sequences (PBSIM3 table)        [Stage 3]
   v
per cell: art_modern reads -> samtools -> rsem-tbam2gbam [Stage 4]
   v
per-cell genome BAMs  +  ground-truth tables
```

## Quick start

```bash
# 1. environment
mamba env create -f workflow/envs/r.yaml
conda activate scARTist-r
Rscript -e 'remotes::install_github("JSB-UCLA/scDesign2")'
pip install snakemake

# 2. build the tools container (art_modern + samtools + rsem + seqkit)
apptainer build resources/art.sif containers/art_modern.def

# 3. put your reference inputs in resources/ and edit config/config.yaml
#    (see docs/inputs.md)

# 4. dry-run, then run
snakemake -n
snakemake --cores 40 --use-conda
```

Stop after count simulation (skip reads):

```bash
snakemake counts_only --cores 8 --use-conda
```

## Documentation

- [docs/installation.md](docs/installation.md) - dependencies and container build
- [docs/inputs.md](docs/inputs.md) - reference inputs you must supply
- [docs/usage.md](docs/usage.md) - configuration and running the sweep
- [docs/outputs.md](docs/outputs.md) - output files and ground-truth format

## Scope and related tools

`scARTist` simulates **full-length (Smart-seq2-like) per-cell short-read
libraries**; it does **not** model UMIs or cell barcodes and is therefore not a
Smart-seq3 or droplet-protocol simulator. 

## Citing

If you use `scARTist`, please cite this pipeline and its core dependencies:
scDesign2 (Sun et al., Genome Biol 2021), art_modern (Yu et al., 2026), and RSEM
(Li & Dewey, BMC Bioinformatics 2011).

## License

See [LICENSE](LICENSE) - to be finalized before public release.
