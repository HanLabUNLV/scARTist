#!/usr/bin/env bash
# Timed COLD full-pipeline run: 80-cell null-arm config, S1..S4 end to end.
# Wipes this config's results, then builds all 80 genome BAMs in one invocation
# with cells fanned across cores (20 concurrent x 4 threads = 80 cores). The
# overall /usr/bin/time gives cold wall-clock + peak RSS; the Snakemake log's
# per-rule timestamps decompose it into S1 (fit) / S2 (counts) / S3+S4 (reads).
set -euo pipefail
cd /mnt/data1/home/mirahan/scRNAsim/scARTist

SMK=/mnt/data1/home/mirahan/miniconda3/envs/beers2/bin/snakemake
PROJ_ROOT=/mnt/data1/home/mirahan/scRNAsim
CORES=80

# fresh cold start: remove all derived results for this config
rm -rf results/groundtruth results/sim_counts/rd_10000000_cells_80 \
       results/pbsim3/rd_10000000_cells_80 results/reads/rd_10000000_cells_80

TARGETS=""
for i in $(seq 1 80); do
  TARGETS="$TARGETS results/reads/rd_10000000_cells_80/null/Bcell_genome.${i}.bam"
done

echo "[$(date '+%F %T')] COLD full-pipeline timing: 80 cells, ${CORES} cores"
/usr/bin/time -v "${SMK}" --cores "${CORES}" --use-singularity \
    --singularity-args "--bind ${PROJ_ROOT}" \
    ${TARGETS}
echo "[$(date '+%F %T')] done"
echo "=== genome BAMs produced ==="
ls results/reads/rd_10000000_cells_80/null/*.bam | wc -l
