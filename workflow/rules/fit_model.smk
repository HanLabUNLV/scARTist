# Stage 1: fit the scDesign2 model on the real matrix (runs ONCE).
# Produces the copula model + ground-truth DE/DS tables consumed by all configs.

rule fit_model:
    input:
        matrix=config["inputs"]["real_count_matrix"],
        mapped_events=config["inputs"]["mapped_events"],
    output:
        copula="results/groundtruth/estimated.copula_result.RData",
        fold_change="results/groundtruth/isoform_with_gene_fold_change.txt",
        ds_genes="results/groundtruth/ds_genes.txt",
    params:
        outdir="results/groundtruth",
        seed=config["model"]["seed"],
        cell_type=config["model"]["cell_type"],
        marginal=config["model"]["marginal"],
        read_len=config["model"]["read_len"],
        nbr_diff_expr=config["model"]["nbr_diff_expr"],
        nbr_diff_spliced=config["model"]["nbr_diff_spliced"],
        ncores=config["model"]["fit_ncores"],
    threads: config["model"]["fit_ncores"]
    container:
        RSIF
    log:
        "logs/fit_model.log",
    script:
        "../scripts/fit_scdesign2.R"
