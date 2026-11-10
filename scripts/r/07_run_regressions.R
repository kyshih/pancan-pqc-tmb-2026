#!/usr/bin/env Rscript

# 07_run_regressions.R
# Runs both global and per-gene mixed-effects regressions on TCGA expression data.

# ---- 1. Load libraries ----
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(arrow)
})

# ---- 2. Parse arguments ----
option_list <- list(
  make_option("--input", type = "character", help = "Path to input CSV"),
  make_option("--output_data_long", type = "character", help = "Output melted df directory"),
  make_option("--output_global", type = "character", help = "Output global regression result"),
  make_option("--output_genewise", type = "character", help = "Output per-gene regression result"),
  make_option("--genes", type = "character",
            default = "all",
            help = "Comma-separated list of genes to analyze, or 'all' for all genes")

)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- 3. Load wide data ----
message("📂 Reading input file: ", opt$input)
#df <- fread(opt$input)

# ---- 4. Melt to long format ----
metadata_cols <- c("barcode", "cancer_type", "purity", "log_mut")

# df_long <- melt(
#  df,
#  id.vars = metadata_cols,       # columns to keep fixed
#  variable.name = "gene",        # column name for former headers
#  value.name   = "expr_z"        # column name for values
# )

# ---- 5. (optional) Save melted version ----
# write_parquet(df_long, opt$output_data_long, compression = "zstd")
# Load long data
df_long <- read_parquet(opt$output_data_long)
head(df_long)
dim(df_long)

# ---- 6. Source regression functions ----
source("/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected/scripts/r/regression_functions.R")
source("/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected/scripts/r/append_skipped_genes.R")

# ---- 7. Get genes of interest ----
if (opt$genes == "all") {
  genes_to_run <- "all"
} else {
  genes_to_run <- strsplit(opt$genes, ",")[[1]] |> trimws()
}

# ---- 4. Run regressions ----
# global_res <- run_global_regression(df_long, opt$output_global)
#print(opt$output_global)
genewise_res <- run_genewise_regression(df_long, genes_to_run, opt$output_genewise)
#print(genes_to_run)
#print(opt$output_genewise)