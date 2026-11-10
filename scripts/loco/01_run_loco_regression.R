#!/usr/bin/env Rscript

# 05_run_loco_regressions.R
#

# ---- 1. Load libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
  library(readr)
  library(optparse)
})

# ---- 2. Parse arguments ----
option_list <- list(
  make_option("--input_expr_long", type = "character", help = "Path to input expression file long format"),
  make_option("--input_beta", type = "character", help = "Path to input GLMM beta"),
  make_option("--reg_helper_dir", type = "character", help = "Dir path to regression functions"),
  make_option("--output_dir", type = "character", help = "Output loco directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- 3. Load input data ----
message("📂 Reading input file: ", opt$input_expr_long)

df <- read_parquet(opt$input_expr_long)
pos_beta_df <- read_csv(opt$input_beta, show_col_types = FALSE)
reg_helper_dir <- opt$reg_helper_dir
output_dir <- opt$output_dir

pos_beta_genes = unique(pos_beta_df$gene)

# ---- 4. Source regression functions ----
source(file.path(reg_helper_dir, "regression_functions.R"))
source(file.path(reg_helper_dir, "append_skipped_genes.R"))

# ---- 5. Run loco regressions ----
cancers <- unique(df$cancer_type)

for (ct in cancers){

  message(paste("Running LOCO — leaving out:", ct))

  df_train <- df %>%
    filter(cancer_type != ct)

  output_file <- paste0(output_dir, "/", "loco_", ct, ".csv")

  run_genewise_regression(
    df = df_train,
    genes = pos_beta_genes,
    output_file = output_file
  )
}