#!/usr/bin/env Rscript

# ===============================================================
# 08_run_permutation_regression.R
# ---------------------------------------------------------------
# This script performs a within-cancer-type permutation of log_mut
# (mutation load) for each gene and re-runs the mixed-effects
# regression to generate a null distribution of slopes.
# ===============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(arrow)
})

# ---- 1. Parse arguments ----
option_list <- list(
  make_option("--input_data_long", type = "character", help = "Input long-format data (parquet)"),
  make_option("--output_data_long_perm", type = "character", help = "Output permuted long-format data (parquet)"),
  make_option("--output_genewise_perm", type = "character", help = "Output file for permuted regression results (csv)"),
  make_option("--genes", type = "character", default = "all", help = "Comma-separated list of genes, or 'all'")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- 2. Read input data ----
message("📥 Reading input data...")
df_long <- read_parquet(opt$input_data_long)
message(sprintf("✅ Loaded %d rows, %d columns", nrow(df_long), ncol(df_long)))

# ---- 3. Source regression functions ----
source("/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected/scripts/r/regression_functions.R")
source("/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected/scripts/r/append_skipped_genes.R")

# ---- 4. Determine genes to run ----
if (opt$genes == "all") {
  genes_to_run <- "all"
} else {
  genes_to_run <- strsplit(opt$genes, ",")[[1]] |> trimws()
}

# ---- 5. Permute log_mut within gene × cancer_type ----
set.seed(2025)
message("🔀 Permuting log_mut within each gene × cancer_type...")
df_perm <- df_long %>%
  group_by(gene, cancer_type) %>%
  mutate(log_mut = sample(log_mut)) %>%
  ungroup()

write_parquet(df_perm, opt$output_data_long_perm, compression = "zstd")

# ---- 6. Run gene-wise regression on permuted data ----
message("🏃 Running permuted regressions...")
genewise_res_perm <- run_genewise_regression(df_perm, genes_to_run, opt$output_genewise_perm)

message("🎉 Done! Permuted results written to:")
message(opt$output_genewise_perm)
