# Stage 2: simulate per-cell isoform counts for one (read_depth, n_cells) config.
# One run writes BOTH arms (null and deds) plus the refit marginal parameters.

rule simulate_counts:
    input:
        copula="results/groundtruth/estimated.copula_result.RData",
        fold_change="results/groundtruth/isoform_with_gene_fold_change.txt",
        ds_genes="results/groundtruth/ds_genes.txt",
    output:
        deds="results/sim_counts/rd_{rd}_cells_{nc}/simulation_deds.txt",
        null="results/sim_counts/rd_{rd}_cells_{nc}/simulation_null.txt",
    params:
        groundtruth_dir="results/groundtruth",
        out_counts_dir="results/sim_counts/rd_{rd}_cells_{nc}",
        out_marginal_dir="results/marginal_params/rd_{rd}_cells_{nc}",
        read_depth=lambda w: int(w.rd),
        n_cells=lambda w: int(w.nc),
        seed=config["model"]["seed"],
        cell_type=config["model"]["cell_type"],
        marginal=config["model"]["marginal"],
        ncores=config["model"]["fit_ncores"],
    threads: config["model"]["fit_ncores"]
    conda:
        "../envs/r.yaml"
    log:
        "logs/simulate_counts/rd_{rd}_cells_{nc}.log",
    script:
        "../scripts/simulate_counts.R"
