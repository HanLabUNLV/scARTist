# Outputs

All outputs are written under `results/` (git-ignored).

```
results/
  groundtruth/
    estimated.copula_result.RData        fitted scDesign2 model (Stage 1)
    isoform_with_gene_fold_change.txt     per-isoform table incl. DE fold changes + truth flags
    ds_genes.txt                          genes chosen for differential transcript usage
    original_marginal_param*.txt          fitted NB marginals (provenance)
  sim_counts/rd_<rd>_cells_<nc>/
    simulation_deds.txt                   per-cell isoform counts, DE + DTU arm
    simulation_null.txt                   per-cell isoform counts, null arm
  marginal_params/rd_<rd>_cells_<nc>/
    {null,deds}_marginal_param*.txt       refit marginals for each simulated arm
  pbsim3/rd_<rd>_cells_<nc>/<arm>/
    transcripts.pbsim3_table.tsv          transcriptID + per-cell counts + sequence
  reads/rd_<rd>_cells_<nc>/<arm>/
    Bcells_trans_pbsim3_pe.<cell>.fastq   simulated paired-end reads
    Bcell_genome.<cell>.bam               genome-coordinate BAM (terminal deliverable)
```

## Ground-truth columns

In `simulation_deds.txt` (and the groundtruth table), the truth labels are:

- `gene_de_status` - 1 if the gene was given a differential-expression fold change.
- `fold_change` - the DE fold change applied to the gene.
- `gene_ds_status` - 1 if the gene has differential transcript usage.
- `transcript_ds_status` - 1 for the isoforms whose proportions were swapped.
- `diff_IsoPct` - gap between the two most-abundant isoform proportions.
- `IsoPct`, `TPM`, `FPKM` - per-isoform relative/absolute abundance.
- Per-cell count columns named `<cell_type>.1 ... <cell_type>.N`.

The `null` arm carries no injected DE/DTU signal and is the negative control for
benchmarking false-positive rates.

## What downstream (separate project) consumes

The per-cell `Bcell_genome.<cell>.bam` files are the input to the differential
splicing benchmarking project (DEXSeq / rMATS / GrASE / glmmTMB), which is not
part of this pipeline.
