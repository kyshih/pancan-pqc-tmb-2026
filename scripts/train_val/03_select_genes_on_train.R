#!/usr/bin/env Rscript

# 03_select_genes_on_train.R
# Apply BH FDR on train-fit p-values and select posBeta axis genes
# (estimate > 0, fdr_bh < threshold). Matches full-data axis selection logic.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
})

option_list <- list(
  make_option("--input_genewise", type = "character",
              help = "Path to genewise regression output CSV"),
  make_option("--output_with_fdr", type = "character",
              help = "Output CSV with fdr_bh column appended"),
  make_option("--output_axis_genes", type = "character",
              help = "Output TSV of selected axis genes"),
  make_option("--fdr_threshold", type = "double", default = 0.05,
              help = "BH FDR threshold for axis gene selection")
)
opt <- parse_args(OptionParser(option_list = option_list))

df <- read_csv(opt$input_genewise, show_col_types = FALSE) %>%
  filter(term == "log_mut")

df$fdr_bh <- p.adjust(df$p.value, method = "BH")

write_csv(df, opt$output_with_fdr)

selected <- df %>%
  filter(fdr_bh < opt$fdr_threshold, estimate > 0) %>%
  arrange(gene) %>%
  select(gene, estimate, std.error, t.value, p.value, fdr_bh, conf.low, conf.high)

write_tsv(selected, opt$output_axis_genes)

message(sprintf(
  "BH FDR selection: %d / %d genes pass (estimate > 0, fdr_bh < %.2f)",
  nrow(selected), nrow(df), opt$fdr_threshold
))
message(sprintf("Saved FDR results → %s", opt$output_with_fdr))
message(sprintf("Saved axis genes   → %s", opt$output_axis_genes))
