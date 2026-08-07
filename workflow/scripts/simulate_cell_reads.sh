#!/bin/bash
# Stage 4 worker: simulate reads for a single cell.
# Ported from simsc.sh. Changes vs original:
#   - the sequence column is the LAST column (NF), computed at run time, instead
#     of the stale hardcoded column 502
#   - the cell's count column is passed in (--count-col) instead of `cut -f<n>`
#     on a fixed layout
#   - all paths and art/rsem parameters are arguments, no hardcoded --bind path
set -euo pipefail

# ---- parse args -----------------------------------------------------------
# Default: keep ONLY the genome BAM. The per-read FASTQ (~7x the genome BAM size:
# 3.7G vs 500M/cell) is never read downstream -- rMATS/DEXSeq/GrASE use the genome
# BAM (built from --o-sam), not the reads. Pass --keep-fastq to retain it.
KEEP_FASTQ=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep-fastq)   KEEP_FASTQ=1; shift 1;;
        --table)        TABLE="$2"; shift 2;;
        --cell)         CELL="$2"; shift 2;;
        --count-col)    COUNT_COL="$2"; shift 2;;
        --outdir)       OUTDIR="$2"; shift 2;;
        --sif)          SIF="$2"; shift 2;;
        --bind)         BIND="$2"; shift 2;;
        --rsem-prefix)  RSEM_PREFIX="$2"; shift 2;;
        --qual1)        QUAL1="$2"; shift 2;;
        --qual2)        QUAL2="$2"; shift 2;;
        --read-len)     READ_LEN="$2"; shift 2;;
        --fcov)         FCOV="$2"; shift 2;;
        --frag-mean)    FRAG_MEAN="$2"; shift 2;;
        --frag-sd)      FRAG_SD="$2"; shift 2;;
        --art-threads)  ART_THREADS="$2"; shift 2;;
        --rsem-threads) RSEM_THREADS="$2"; shift 2;;
        *) echo "Unknown argument: $1" >&2; exit 1;;
    esac
done

mkdir -p "${OUTDIR}"
PBSIM3="${OUTDIR}/cell.${CELL}.pbsim3"
TBAM="${OUTDIR}/Bcells_trans_pbsim3_pe.${CELL}.bam"
SORTED="${OUTDIR}/Bcells_trans_pbsim3_pe.sorted.${CELL}.bam"
GBAM="${OUTDIR}/Bcell_genome.${CELL}.bam"
FASTQ="${OUTDIR}/Bcells_trans_pbsim3_pe.${CELL}.fastq"

# ---- build this cell's 4-column PBSIM3 table ------------------------------
# Columns of TABLE: 1=transcriptID, 2..N+1=per-cell counts, NF=sequence.
# PBSIM3 transcripts format: ID <tab> COV_POS <tab> COV_NEG <tab> SEQ
awk -v c="${COUNT_COL}" 'BEGIN{OFS="\t"} {print $1, $c, 0, $NF}' "${TABLE}" > "${PBSIM3}"

# singularity exec wrapper (add --bind only if provided)
BIND_ARG=()
[[ -n "${BIND}" ]] && BIND_ARG=(--bind "${BIND}")
sing() { singularity exec "${BIND_ARG[@]}" "${SIF}" "$@"; }

# ---- art_modern: PBSIM3 transcripts -> paired-end reads + transcriptome BAM
sing art_modern \
    --mode trans \
    --lc pe \
    --i-file "${PBSIM3}" \
    --o-fastq "${FASTQ}" \
    --o-sam-write_bam \
    --o-sam "${TBAM}" \
    --i-type pbsim3_transcripts \
    --qual_file_1 "${QUAL1}" \
    --qual_file_2 "${QUAL2}" \
    --read_len "${READ_LEN}" \
    --parallel "${ART_THREADS}" \
    --i-fcov "${FCOV}" \
    --pe_frag_dist_mean "${FRAG_MEAN}" \
    --pe_frag_dist_std_dev "${FRAG_SD}"

# ---- name-sort, then map transcriptome BAM to genome coordinates ----------
sing samtools sort -n -@ "${ART_THREADS}" "${TBAM}" -o "${SORTED}"
sing rsem-tbam2gbam "${RSEM_PREFIX}" "${SORTED}" "${GBAM}" -p "${RSEM_THREADS}"

# ---- tidy per-cell intermediates ------------------------------------------
# Always drop the pbsim table + transcriptome BAMs (never needed again). Drop the
# FASTQ too unless --keep-fastq -- it is the largest artifact (3.7G vs 500M genome
# BAM) and nothing downstream reads it.
rm -f "${PBSIM3}" "${TBAM}" "${SORTED}"
[[ "${KEEP_FASTQ}" -eq 0 ]] && rm -f "${FASTQ}"
echo "Done: ${GBAM}$([[ ${KEEP_FASTQ} -eq 1 ]] && echo ' (+fastq kept)')"
