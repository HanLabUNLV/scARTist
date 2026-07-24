## Stage 2: simulate per-cell isoform counts for one (read_depth, n_cells) config.
## Ported from diff_params_4/simulate_sc_isoforms.han.R.
## Changes vs original:
##   - read_depth and n_cells come from Snakemake wildcards (not argv)
##   - all input/output directories come from Snakemake params
##   - output dirs are created before writing

suppressPackageStartupMessages({
  library(gtools)
  library(dplyr)
  library(copula)
  library(scDesign2)
  library(tictoc)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")
tic.clearlog()

p <- snakemake@params
gt_dir       <- p[["groundtruth_dir"]]
counts_dir   <- p[["out_counts_dir"]]
marginal_dir <- p[["out_marginal_dir"]]
rd_num       <- as.numeric(p[["read_depth"]])
n_cell_new   <- as.integer(p[["n_cells"]])
seed         <- as.integer(p[["seed"]])
cell_type    <- p[["cell_type"]]
marginal     <- p[["marginal"]]
ncores       <- as.integer(p[["ncores"]])

dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(marginal_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)
isoform_with_gene_fold_change <- read.table(
  file.path(gt_dir, "isoform_with_gene_fold_change.txt"), header = TRUE, sep = "\t")
ds_genes <- read.table(file.path(gt_dir, "ds_genes.txt"), header = TRUE, sep = "\t")[, 1]
load(file.path(gt_dir, "estimated.copula_result.RData"))  # -> copula_result

write_marginals <- function(fit, prefix) {
  ct <- fit[[cell_type]]
  write.table(ct$marginal_param1, file.path(marginal_dir, paste0(prefix, "_marginal_param1.txt")),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(ct$marginal_param2, file.path(marginal_dir, paste0(prefix, "_marginal_param2.txt")),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  for (i in 1:3) {
    write.table(ct[[paste0("gene_sel", i)]],
                file.path(marginal_dir, sprintf("%s_gene_sel%d.txt", prefix, i)),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
}

## ---- null arm -------------------------------------------------------------
message("Simulating null isoform counts...")
sim_count_null <- simulate_count_scDesign2(copula_result, n_cell_new,
  total_count_new = rd_num, cell_type_prop = 1, sim_method = "copula")
null_fit <- fit_model_scDesign2(sim_count_null, cell_type, sim_method = "copula",
  marginal = marginal, ncores = ncores)
write_marginals(null_fit, "null")

## ---- DE arm ---------------------------------------------------------------
message("Simulating isoform counts with DE...")
sim_count_de <- simulate_count_scDesign2(copula_result, n_cell_new,
  total_count_new = rd_num, cell_type_prop = 1, sim_method = "copula",
  fold_change = isoform_with_gene_fold_change$fold_change)
deds_fit <- fit_model_scDesign2(sim_count_de, cell_type, sim_method = "copula",
  marginal = marginal, ncores = ncores)
write_marginals(deds_fit, "deds")

## ---- summarize both arms to TPM/IsoPct ------------------------------------
summarize_arm <- function(counts) {
  df <- as.data.frame(isoform_with_gene_fold_change)
  df$expected_count <- rowMeans(counts)
  df$frac <- df$expected_count / sum(df$expected_count)
  df$FPKM <- df$frac * 1e9 / df$effective_length
  df$FPKM[df$effective_length == 0] <- 0
  df$TPM <- df$FPKM / sum(df$FPKM) * 1e6
  df <- df %>% group_by(geneID) %>% mutate(IsoPct = TPM * 100 / sum(TPM)) %>% ungroup()
  df$IsoPct[is.na(df$IsoPct)] <- 0
  df
}
simCountsNull <- summarize_arm(sim_count_null)
simCountsDE   <- summarize_arm(sim_count_de)

final_null <- cbind.data.frame(simCountsNull, sim_count_null)
final_deds <- cbind.data.frame(simCountsDE, sim_count_de)

## ---- inject differential splicing into the deds arm -----------------------
## Swap the two most-abundant isoform proportions for the chosen ds_genes.
message("Introducing differential splicing...")
mutate_IsoPct <- function(w, ref) {
  o <- order(ref, decreasing = TRUE)[1:2]
  w[rev(o)] <- w[o]
  w
}
mutated_IsoPct <- function(ref) {
  o <- order(ref, decreasing = TRUE)[1:2]
  w <- rep(0, length(ref)); w[o] <- 1; w
}

isoform_summary_nonds <- final_deds[!(final_deds$geneID %in% ds_genes), ]
isoform_summary_ds    <- final_deds[final_deds$geneID %in% ds_genes, ]
cn <- paste(cell_type, seq_len(n_cell_new), sep = ".")
colnames(isoform_summary_ds)[colnames(isoform_summary_ds) == cell_type] <- cn
colnames(isoform_summary_nonds)[colnames(isoform_summary_nonds) == cell_type] <- cn

isoform_summary_ds <- isoform_summary_ds %>% group_by(geneID) %>%
  mutate(across(contains(cell_type), ~ mutate_IsoPct(.x, IsoPct)),
         across(c("expected_count", "frac", "FPKM", "TPM"), ~ mutate_IsoPct(.x, IsoPct)),
         diff_IsoPct = -diff(sort(IsoPct, decreasing = TRUE))[1] / 100,
         IsoPct = mutate_IsoPct(IsoPct, IsoPct),
         gene_ds_status = 1,
         transcript_ds_status = mutated_IsoPct(IsoPct))
isoform_summary_nonds <- isoform_summary_nonds %>% group_by(geneID) %>%
  mutate(diff_IsoPct = -diff(sort(IsoPct, decreasing = TRUE))[1] / 100,
         gene_ds_status = 0, transcript_ds_status = 0)
isoform_summary_nonds$diff_IsoPct[is.na(isoform_summary_nonds$diff_IsoPct)] <- 0

final_deds <- rbind(as.data.frame(isoform_summary_nonds), as.data.frame(isoform_summary_ds))

## ---- write ----------------------------------------------------------------
message("Writing result files...")
write.table(final_deds, file.path(counts_dir, "simulation_deds.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(final_null, file.path(counts_dir, "simulation_null.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
message("Stage 2 complete.")
