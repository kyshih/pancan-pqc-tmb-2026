# ---- 1. Load libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(optparse)
})

# directory containing LOCO outputs
# ---- 2. Parse arguments ----
option_list <- list(
  make_option("--loco_dir", type = "character", help = "Path to directory containing LOCO outputs"),
  make_option("--output_combined", type = "character", help = "output combined file name"),
  make_option(
    "--input_beta",
    type = "character",
    default = NULL,
    help = "Optional posBeta gene list; only these genes are kept in loco_combined"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

loco_dir <- opt$loco_dir
outfname  <- opt$output_combined

allowed_genes <- NULL
if (!is.null(opt$input_beta)) {
  beta_df <- read_csv(opt$input_beta, show_col_types = FALSE)
  allowed_genes <- unique(beta_df$gene)
  message(sprintf(
    "Filtering combined output to %d genes from %s",
    length(allowed_genes), basename(opt$input_beta)
  ))
}

# list files
files <- list.files(
  loco_dir,
  pattern = "^loco_[A-Z]+\\.csv$",
  full.names = TRUE
)

print(files)

# read and combine
loco_combined <- map_dfr(files, function(f){

  df <- read_csv(f)

  held_out <- basename(f) |>
    str_remove("loco_") |>
    str_remove("\\.csv")

    out <- df %>%
      filter(term == "log_mut") %>%
      select(gene, estimate) %>%
      mutate(held_out_cancer = held_out)

    if (!is.null(allowed_genes)) {
      out <- out %>% filter(gene %in% allowed_genes)
    }

    out
})

# save combined table
write.csv(
  loco_combined,
  outfname,
  row.names = FALSE
)