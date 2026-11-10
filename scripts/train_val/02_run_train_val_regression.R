#!/usr/bin/env Rscript

# 02_run_train_val_regression.R
# Refit per-gene pancan GLMM on training samples only.

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
  library(readr)
  library(optparse)
})

option_list <- list(
  make_option("--input_expr_long", type = "character",
              help = "Path to long-format expression parquet"),
  make_option("--input_split", type = "character",
              help = "Path to train/val split CSV"),
  make_option("--input_beta", type = "character", default = NULL,
              help = "Optional path to gene list (TSV/CSV with gene column). Default: all genes."),
  make_option("--reg_helper_dir", type = "character",
              help = "Directory containing regression helper R scripts"),
  make_option("--output_genewise", type = "character",
              help = "Output path for per-gene regression results")
)
opt <- parse_args(OptionParser(option_list = option_list))

message("Reading expression: ", opt$input_expr_long)
df <- read_parquet(opt$input_expr_long)

message("Reading split: ", opt$input_split)
split_df <- read_csv(opt$input_split, show_col_types = FALSE)

train_barcodes <- split_df %>%
  filter(split == "train") %>%
  pull(barcode)

df_train <- df %>%
  filter(barcode %in% train_barcodes)

message(sprintf(
  "Training samples: %d rows (%d unique barcodes, %d cancer types)",
  nrow(df_train),
  n_distinct(df_train$barcode),
  n_distinct(df_train$cancer_type)
))

source(file.path(opt$reg_helper_dir, "regression_functions.R"))
source(file.path(opt$reg_helper_dir, "append_skipped_genes.R"))

genes_to_run <- "all"
if (!is.null(opt$input_beta)) {
  beta_df <- read_csv(opt$input_beta, show_col_types = FALSE)
  genes_to_run <- unique(beta_df$gene)
  message(sprintf("Restricting to %d genes from %s", length(genes_to_run), opt$input_beta))
}

run_genewise_regression(
  df = df_train,
  genes = genes_to_run,
  output_file = opt$output_genewise
)
