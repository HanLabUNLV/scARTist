# Usage

## Configure the sweep

Edit `config/config.yaml`:

```yaml
sweep:
  read_depths: [10000000, 41000000, 83000000]
  n_cells: [80, 150]
  arms: ["null", "deds"]
```

The workflow runs the full cross-product: `read_depths x n_cells x arms`, and
within each config fans out one job **per cell**. With the defaults that is
3 x 2 = 6 count configs and, for the read stage,
`(80 + 150) x 3 depths x 2 arms = 1380` per-cell read jobs.

## Run

```bash
snakemake -n                      # dry run - inspect the plan
snakemake --cores 40 --use-conda  # full run
```

Useful targets and flags:

```bash
snakemake counts_only --cores 8 --use-conda        # stop after Stage 2
snakemake results/reads/rd_10000000_cells_80/deds/Bcell_genome.1.bam --cores 4
snakemake --dag | dot -Tsvg > dag.svg              # visualize the DAG
```

On a cluster, add a Snakemake profile (`--profile <profile>`); each per-cell job
requests `reads.art_threads` cores and the model steps request
`model.fit_ncores`.

## Output location (`--directory`)

All outputs are written under `results/` **relative to the working directory**,
including the shared `results/groundtruth/` (one fit per study). To put a study's
outputs somewhere else -- or to run two independent studies without them sharing
(and clobbering) the same `results/groundtruth/` -- give each its own working
directory with Snakemake's built-in `--directory`:

```bash
snakemake --snakefile /path/to/scARTist/workflow/Snakefile \
          --directory /path/to/studyA_out \
          --configfile /path/to/config.yaml \
          --cores 40 --use-singularity
```

Everything (`groundtruth/`, `sim_counts/`, `reads/`) then lands under
`/path/to/studyA_out/results/...`. Package resources (`r.sif`, `art.sif`, the
read-simulation script) are resolved relative to the workflow directory, not the
working directory, so `--directory` relocates outputs without breaking the run.
Point `inputs:` at absolute paths (or paths valid from the working directory),
and make sure `container.bind` (and any `--singularity-args --bind`) covers the
chosen output directory so the container can write there.

## Stages

| Stage | Rule | Runs | Output |
|---|---|---|---|
| 1 fit | `fit_model` | once | `results/groundtruth/` |
| 2 counts | `simulate_counts` | per (depth, cells) | `results/sim_counts/rd_*_cells_*/simulation_{null,deds}.txt` |
| 3 pbsim3 | `build_pbsim3_table` | per (depth, cells, arm) | `results/pbsim3/.../transcripts.pbsim3_table.tsv` |
| 4 reads | `simulate_cell_reads` | per cell | `results/reads/.../Bcell_genome.{cell}.bam` |

## Reproducibility

Set `model.seed` for deterministic fits/simulations. After the first run, set
`model.use_cached_biomart: true` so the gene->chromosome annotation no longer
depends on a live network call.
