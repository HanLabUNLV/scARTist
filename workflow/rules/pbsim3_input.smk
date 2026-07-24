# Stage 3: attach transcript sequences to the simulated counts, producing a
# combined PBSIM3-ready table: col1 = transcriptID, cols 2..N+1 = per-cell
# counts, last col = transcript nucleotide sequence.
#
# This replaces join_cnts_seq.R. The original hardcoded 80 cells via the column
# slice [,c(2,22:101,103)]; here count columns are selected by name (^Bcell)
# so any cell number works.

rule build_pbsim3_table:
    input:
        counts="results/sim_counts/rd_{rd}_cells_{nc}/simulation_{arm}.txt",
        gencode_tab=config["inputs"]["gencode_transcripts_tab"],
    output:
        table="results/pbsim3/rd_{rd}_cells_{nc}/{arm}/transcripts.pbsim3_table.tsv",
    conda:
        "../envs/r.yaml"
    log:
        "logs/build_pbsim3/rd_{rd}_cells_{nc}_{arm}.log",
    script:
        "../scripts/build_pbsim3_table.R"
