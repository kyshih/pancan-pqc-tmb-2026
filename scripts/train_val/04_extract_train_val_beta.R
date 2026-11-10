#!/usr/bin/env Rscript

# 04_extract_train_val_beta.R
# Extract log_mut coefficients from train-fit genewise GLMM results.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
})

option_list <- list(
  make_option("--input_genewise", type = "character",
              help = "Path to genewise regression output CSV"),
  make_option("--output_beta", type = "character",
              help = "Output TSV with gene and estimate columns"),
  make_option("--gene_list", type = "character", default = NULL,
              help = "Optional TSV/CSV with gene column to restrict output")
)
opt <- parse_args(OptionParser(option_list = option_list))

beta_df <- read_csv(opt$input_genewise, show_col_types = FALSE) %>%
  filter(term == "log_mut")

if (!is.null(opt$gene_list)) {
  if (grepl("\\.csv$", opt$gene_list, ignore.case = TRUE)) {
    gene_list <- read_csv(opt$gene_list, show_col_types = FALSE)
  } else {
    gene_list <- read_tsv(opt$gene_list, show_col_types = FALSE)
  }
  gene_list <- gene_list %>% pull(gene) %>% unique()
  beta_df <- beta_df %>% filter(gene %in% gene_list)
  message(sprintf("Restricted to %d genes from %s", length(gene_list), opt$gene_list))
}

beta_df <- beta_df %>%
  select(gene, estimate, std.error, t.value, p.value, any_of("fdr_bh"), conf.low, conf.high) %>%
  arrange(gene)

write_tsv(beta_df, opt$output_beta)
message(sprintf("Saved %d gene betas → %s", nrow(beta_df), opt$output_beta))
