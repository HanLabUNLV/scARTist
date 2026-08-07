# Stage 4: per-cell read simulation.
# For cell {cell}: extract that cell's coverage + sequence from the combined
# table, run art_modern (PBSIM3 transcripts mode), name-sort, and map the
# transcriptome BAM to genome space with rsem-tbam2gbam.

rule simulate_cell_reads:
    input:
        table="results/pbsim3/rd_{rd}_cells_{nc}/{arm}/transcripts.pbsim3_table.tsv",
        sif=ARTSIF,
        qual1=config["inputs"]["qual_profile_r1"],
        qual2=config["inputs"]["qual_profile_r2"],
    output:
        genome_bam="results/reads/rd_{rd}_cells_{nc}/{arm}/Bcell_genome.{cell}.bam",
    params:
        workdir="results/reads/rd_{rd}_cells_{nc}/{arm}",
        # cell {cell} is 1-indexed; its count column in the table is 1 + cell.
        count_col=lambda w: 1 + int(w.cell),
        rsem_prefix=config["inputs"]["rsem_reference_prefix"],
        bind=config["container"]["bind"],
        read_len=config["reads"]["read_len"],
        fcov=config["reads"]["fcov"],
        frag_mean=config["reads"]["pe_frag_dist_mean"],
        frag_sd=config["reads"]["pe_frag_dist_std_dev"],
        art_threads=config["reads"]["art_threads"],
        rsem_threads=config["reads"]["rsem_threads"],
        # 1 keeps the per-read FASTQ; 0 (default) deletes it -- it is ~7x the genome BAM
        # (3.7G vs 500M/cell) and nothing downstream reads it (genome BAM comes from --o-sam).
        keep_fastq=(1 if config["reads"].get("keep_fastq", False) else 0),
    threads: config["reads"]["art_threads"]
    log:
        "logs/simulate_reads/rd_{rd}_cells_{nc}_{arm}_cell{cell}.log",
    # Inlined per-cell read sim (formerly scripts/simulate_cell_reads.sh): build this cell's
    # 4-col PBSIM3 table -> art_modern (transcripts PE) -> name-sort -> rsem-tbam2gbam into
    # genome space. Kept as ONE job per cell (3 container exec calls: art, sort, tbam2gbam).
    # Intermediates deleted inline; the FASTQ is dropped unless keep_fastq.
    shell:
        r"""
        exec > {log} 2>&1
        set -euo pipefail
        OUTDIR={params.workdir}
        mkdir -p "$OUTDIR"
        PBSIM3="$OUTDIR/cell.{wildcards.cell}.pbsim3"
        TBAM="$OUTDIR/Bcells_trans_pbsim3_pe.{wildcards.cell}.bam"
        SORTED="$OUTDIR/Bcells_trans_pbsim3_pe.sorted.{wildcards.cell}.bam"
        FASTQ="$OUTDIR/Bcells_trans_pbsim3_pe.{wildcards.cell}.fastq"
        GBAM={output.genome_bam}

        # this cell's PBSIM3 table: ID <tab> COV_POS <tab> COV_NEG(0) <tab> SEQ (seq = last col)
        awk -v c={params.count_col} 'BEGIN{{OFS="\t"}} {{print $1, $c, 0, $NF}}' {input.table} > "$PBSIM3"

        BIND_ARG=()
        [ -n "{params.bind}" ] && BIND_ARG=(--bind "{params.bind}")
        sing() {{ singularity exec "${{BIND_ARG[@]}}" {input.sif} "$@"; }}

        # art_modern: PBSIM3 transcripts -> paired-end reads (FASTQ) + transcriptome BAM
        sing art_modern \
            --mode trans --lc pe \
            --i-file "$PBSIM3" \
            --o-fastq "$FASTQ" --o-sam-write_bam --o-sam "$TBAM" \
            --i-type pbsim3_transcripts \
            --qual_file_1 {input.qual1} --qual_file_2 {input.qual2} \
            --read_len {params.read_len} --parallel {params.art_threads} --i-fcov {params.fcov} \
            --pe_frag_dist_mean {params.frag_mean} --pe_frag_dist_std_dev {params.frag_sd}

        # name-sort, then map transcriptome BAM to genome coordinates
        sing samtools sort -n -@ {params.art_threads} "$TBAM" -o "$SORTED"
        sing rsem-tbam2gbam {params.rsem_prefix} "$SORTED" "$GBAM" -p {params.rsem_threads}

        # tidy intermediates; drop FASTQ (largest, unused downstream) unless keep_fastq
        rm -f "$PBSIM3" "$TBAM" "$SORTED"
        if [ {params.keep_fastq} -eq 0 ]; then rm -f "$FASTQ"; fi
        echo "Done: $GBAM"
        """
