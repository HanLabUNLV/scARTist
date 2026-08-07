#!/usr/bin/env bash
# Single-cell smoke run for scARTist, fully containerized.
#   - R rules (fit_model, simulate_counts, build_pbsim3_table) run inside r.sif
#     via Snakemake's container: directive + --use-singularity.
#   - the read rule calls `singularity exec art.sif ...` itself (no container:
#     directive, so Snakemake does not double-wrap it).
# Snakemake itself comes from the beers2 env (7.24.0; working attrs). The host
# only needs snakemake + apptainer/singularity.
set -euo pipefail
cd /mnt/data1/home/mirahan/scRNAsim/scARTist

SMK=/mnt/data1/home/mirahan/miniconda3/envs/beers2/bin/snakemake
TARGET=results/reads/rd_10000000_cells_80/null/Bcell_genome.1.bam
# The data dir is a SIBLING of the launch dir, so Snakemake's default cwd bind is
# not enough for containerized R rules -- bind the project root explicitly.
PROJ_ROOT=/mnt/data1/home/mirahan/scRNAsim

echo "[$(date '+%F %T')] runner snakemake: $(${SMK} --version)"
echo "[$(date '+%F %T')] singularity: $(command -v singularity)"
echo "[$(date '+%F %T')] starting smoke run -> ${TARGET}"
"${SMK}" --cores 20 --use-singularity \
    --singularity-args "--bind ${PROJ_ROOT}" \
    -p "${TARGET}"
echo "[$(date '+%F %T')] smoke run finished"
ls -la "${TARGET}"
