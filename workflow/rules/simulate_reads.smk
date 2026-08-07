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
        keep_fastq_flag=("--keep-fastq" if config["reads"].get("keep_fastq", False) else ""),
        script=READS_SCRIPT,
    threads: config["reads"]["art_threads"]
    log:
        "logs/simulate_reads/rd_{rd}_cells_{nc}_{arm}_cell{cell}.log",
    shell:
        r"""
        bash {params.script} \
            --table {input.table} \
            --cell {wildcards.cell} \
            --count-col {params.count_col} \
            --outdir {params.workdir} \
            --sif {input.sif} \
            --bind "{params.bind}" \
            --rsem-prefix {params.rsem_prefix} \
            --qual1 {input.qual1} \
            --qual2 {input.qual2} \
            --read-len {params.read_len} \
            --fcov {params.fcov} \
            --frag-mean {params.frag_mean} \
            --frag-sd {params.frag_sd} \
            --art-threads {params.art_threads} \
            --rsem-threads {params.rsem_threads} \
            {params.keep_fastq_flag} \
            > {log} 2>&1
        """
