## Stage 1: fit scDesign2 model + build ground-truth DE/DS tables (run once).
## Ported from diff_params_4/estimate_sc_isoforms.han.R.
## Changes vs original:
##   - all paths/params come from Snakemake (no hardcoded /mnt/storage/... paths)
##   - the output directory is created BEFORE anything is written to it
##     (original called dir.create("groundtruth") AFTER writing into it)
##   - dropped the live biomaRt gene->chromosome lookup: it fetched
##     chromosome_name into gene_summary but that column was never used or
##     written (dead in the original too). Removing it makes Stage 1 offline and
##     leaves every output byte-identical.

suppressPackageStartupMessages({
  library(gtools)
  library(dplyr)
  library(fastglm)
  library(copula)
  library(scDesign2)
  library(tictoc)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")
tic.clearlog()
print(sessionInfo())

## ---- parameters -----------------------------------------------------------
p <- snakemake@params
isoform_scCounts_file <- snakemake@input[["matrix"]]
mapped_events_file    <- snakemake@input[["mapped_events"]]
outdir                <- p[["outdir"]]

seed             <- as.integer(p[["seed"]])
cell_type        <- p[["cell_type"]]
marginal         <- p[["marginal"]]
read_len         <- as.numeric(p[["read_len"]])
nbr_diff_expr    <- as.integer(p[["nbr_diff_expr"]])
nbr_diff_spliced <- as.integer(p[["nbr_diff_spliced"]])
ncores           <- as.integer(p[["ncores"]])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---- read real matrix + parse GENCODE headers -----------------------------
data_mat <- read.table(isoform_scCounts_file, sep = "\t", header = TRUE)
geneinfo <- data_mat[, 1:3]
hdr <- strsplit(geneinfo[, 1], "[|]")
geneinfo$transcriptID <- vapply(hdr, `[`, character(1), 1)
geneinfo$geneID       <- vapply(hdr, `[`, character(1), 2)
geneinfo$length       <- as.numeric(vapply(hdr, `[`, character(1), 7))
geneinfo$effective_length <- geneinfo$length - read_len + 1
data_mat <- round(data.matrix(data_mat[, 4:ncol(data_mat)]))

expected_count <- rowMeans(data_mat)
isoform.scCounts <- data.frame(expected_count)
rownames(isoform.scCounts) <- geneinfo$transcriptID
isoform.scCounts <- cbind.data.frame(
  isoform.scCounts,
  geneinfo[c("transcriptID", "geneID", "length", "effective_length")])

## reconstruct TPM / isoform percentage
isoform.scCounts$frac <- isoform.scCounts$expected_count / sum(isoform.scCounts$expected_count)
isoform.scCounts$FPKM <- isoform.scCounts$frac * 1e9 / isoform.scCounts$effective_length
isoform.scCounts$FPKM[isoform.scCounts$effective_length == 0] <- 0
isoform.scCounts$TPM  <- isoform.scCounts$FPKM / sum(isoform.scCounts$FPKM) * 1e6
isoform.scCounts <- isoform.scCounts %>% group_by(geneID) %>%
  mutate(IsoPct = TPM * 100 / sum(TPM)) %>% ungroup() %>% as.data.frame()
isoform.scCounts$IsoPct[is.na(isoform.scCounts$IsoPct)] <- 0

message("Creating gene summary...")
gene_summary <- isoform.scCounts %>% group_by(geneID) %>%
  summarise(expected_gene_count_gr1 = sum(expected_count),
            effective_gene_length = sum(effective_length * IsoPct / 100),
            nbr_isoforms = length(IsoPct),
            nbr_expr_isoforms = length(which(IsoPct > 0)),
            nbr_expr_isoforms10 = length(which(IsoPct > 10)))

## ---- introduce differential expression ------------------------------------
set.seed(seed)
fold_changes <- (2 + rexp(nbr_diff_expr, rate = 1))^(c(-1, 1)[round(runif(nbr_diff_expr)) + 1])
gene_summary$expected_gene_count_gr2 <- gene_summary$expected_gene_count_gr1
gene_summary$gene_de_status <- 0
fold_changes_arr <- rep(1, nrow(gene_summary))
if (nbr_diff_expr > 0) {
  diff_expr_genes <- sample(seq_len(nrow(gene_summary)), nbr_diff_expr, replace = FALSE)
  gene_summary$expected_gene_count_gr2[diff_expr_genes] <-
    gene_summary$expected_gene_count_gr1[diff_expr_genes] * fold_changes
  fold_changes_arr[diff_expr_genes] <- fold_changes_arr[diff_expr_genes] * fold_changes
  gene_summary$gene_de_status[diff_expr_genes] <- 1
}
gene_summary$fold_change <- fold_changes_arr
isoform_with_gene_fold_change <- isoform.scCounts %>% left_join(gene_summary, by = "geneID")

## ---- rMATS detectability (from mapped events) -----------------------------
data <- read.table(mapped_events_file, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
tlist <- strsplit(data$transcripts, ",[ ]*", perl = TRUE)
gene2transcript <- setNames(tlist, data$GeneID)
map <- unique(stack(gene2transcript))
transcript2gene <- with(map, split(as.character(ind), values))

gene_summary$rMATSfound_gene <-
  !sapply(gene2transcript[gene_summary$geneID], is.null)
isoform_with_gene_fold_change$rMATSfound_transcript <-
  !sapply(transcript2gene[isoform_with_gene_fold_change$transcriptID], is.null)

isoform_with_gene_fold_change <- isoform_with_gene_fold_change %>% group_by(geneID) %>%
  mutate(diff_IsoPct = -diff(sort(IsoPct, decreasing = TRUE))[1] / 100,
         gene_ds_status = 0, transcript_ds_status = 0)

## ---- choose differential-splicing genes -----------------------------------
message("Determining genes with differential splicing...")
sampling_set <- intersect(which(isoform_with_gene_fold_change$nbr_expr_isoforms10 >= 2),
                          which(isoform_with_gene_fold_change$expected_gene_count_gr1 > 50))
sampling_set <- intersect(sampling_set, which(isoform_with_gene_fold_change$rMATSfound_transcript))
if (nbr_diff_spliced > length(sampling_set)) nbr_diff_spliced <- length(sampling_set)
ds_genes <- isoform_with_gene_fold_change$geneID[sample(sampling_set, nbr_diff_spliced, replace = FALSE)]

write.table(isoform_with_gene_fold_change, file.path(outdir, "isoform_with_gene_fold_change.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(ds_genes, file.path(outdir, "ds_genes.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## ---- fit scDesign2 copula model -------------------------------------------
message("Estimating scDesign2 isoform-count model...")
colnames(data_mat) <- rep(cell_type, ncol(data_mat))
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
tic()
copula_result <- fit_model_scDesign2(data_mat, cell_type, sim_method = "copula",
                                     marginal = marginal, ncores = ncores)
toc(log = TRUE, quiet = TRUE)

## persist the fitted marginal parameters (provenance)
ct <- copula_result[[cell_type]]
write.table(ct$marginal_param1, file.path(outdir, "original_marginal_param1.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(ct$marginal_param2, file.path(outdir, "original_marginal_param2.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
for (i in 1:3) {
  write.table(ct[[paste0("gene_sel", i)]],
              file.path(outdir, sprintf("original_gene_sel%d.txt", i)),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

save(copula_result, file = file.path(outdir, "estimated.copula_result.RData"))
message("Stage 1 complete.")
