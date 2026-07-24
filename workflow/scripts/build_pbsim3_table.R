## Stage 3: attach transcript sequences to simulated counts -> combined table.
## Ported from join_cnts_seq.R.
## Original bug fixed: it selected columns by hardcoded index [,c(2,22:101,103)],
## which only worked for 80 cells. Here per-cell count columns are selected by
## name (they are named "<cell_type>" or "<cell_type>.<n>"), so any cell count works.
##
## Output layout (no header): transcriptID <tab> count_cell_1 ... count_cell_N <tab> sequence

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

counts_file <- snakemake@input[["counts"]]
gencode_tab <- snakemake@input[["gencode_tab"]]
out_file    <- snakemake@output[["table"]]
cell_type   <- snakemake@config[["model"]][["cell_type"]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

cnts <- read.table(counts_file, header = TRUE, sep = "\t", check.names = FALSE)

gtab <- read.table(gencode_tab, header = FALSE, sep = "\t",
                   quote = "", comment.char = "", stringsAsFactors = FALSE)
gtab$transcript_id <- str_extract(gtab$V1, "^ENST\\d+\\.\\d+(_PAR_Y)?")
gtab <- gtab[, c("transcript_id", "V2")]

merged <- cnts %>% left_join(gtab, by = c("transcriptID" = "transcript_id"))

## select the per-cell count columns by name prefix (robust to cell number)
count_cols <- grep(paste0("^", cell_type), colnames(cnts), value = TRUE)
if (length(count_cols) == 0)
  stop("No per-cell count columns matched prefix '", cell_type, "'")
message(sprintf("Found %d per-cell count columns.", length(count_cols)))

out <- merged[, c("transcriptID", count_cols, "V2")]

## drop transcripts with no sequence (cannot be simulated)
missing_seq <- is.na(out$V2) | out$V2 == ""
if (any(missing_seq))
  message(sprintf("Dropping %d transcripts with no sequence.", sum(missing_seq)))
out <- out[!missing_seq, ]

write.table(out, out_file, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
message("Stage 3 complete: ", out_file)
