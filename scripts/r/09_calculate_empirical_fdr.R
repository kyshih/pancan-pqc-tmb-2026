#!/usr/bin/env Rscript

# ============================================================
# 09_calculate_empirical_fdr.R
# Empirical FDR estimation from permuted β(TMB) regression
# ------------------------------------------------------------
# Inputs:
#   --observed   Path to observed regression results (.csv or .parquet)
#   --permuted   Path to permuted regression results (.csv or .parquet)
#   --output     Path to save output with FDR column
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(arrow)
  library(ggplot2)
})

# ---- CLI ----
option_list <- list(
  make_option("--observed", type="character", help="Observed regression results (CSV or Parquet)"),
  make_option("--permuted", type="character", help="Permuted regression results (CSV or Parquet)"),
  make_option("--output",   type="character", help="Output CSV with empirical FDR appended")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- 1. Load data ----
read_input <- function(path) {
  if (grepl("\\.parquet$", path)) {
    read_parquet(path)
  } else {
    read.csv(path)
  }
}

df_obs  <- read_input(opt$observed)
df_perm <- read_input(opt$permuted)

# ---- 1b. Filter to log10(TMB) term ----
df_obs  <- df_obs  %>% filter(term == "log_mut")
df_perm <- df_perm %>% filter(term == "log_mut")

message("✅ Filtered to term == 'log_mut' (", nrow(df_obs), " observed rows, ",
        nrow(df_perm), " permuted rows)")

# Expect columns: gene, beta (log10(TMB) coefficient)
stopifnot(all(c("gene", "estimate") %in% names(df_obs)))
stopifnot(all(c("gene", "estimate") %in% names(df_perm)))

# ---- 2. Compute empirical FDR curve ----
abs_obs  <- abs(df_obs$estimate)
abs_perm <- abs(df_perm$estimate)

n_obs_total  <- length(abs_obs)
n_perm_total <- length(abs_perm)

t_seq <- seq(0, max(abs_obs, abs_perm), length.out = 500)

fdr_table <- data.frame(
  threshold = t_seq,
  n_obs  = sapply(t_seq, function(t) sum(abs_obs  >= t)),
  n_perm = sapply(t_seq, function(t) sum(abs_perm >= t))
) %>%
  mutate(
    prop_obs  = n_obs  / n_obs_total,
    prop_perm = n_perm / n_perm_total,
    FDR_raw   = prop_perm / pmax(prop_obs, 1e-6)
  ) %>%
  mutate(FDR = cummin(pmin(FDR_raw, 1)))  # monotonic in correct direction

# ---- 2b. Save FDR curve for diagnostics ----
fdr_table_path <- file.path(dirname(opt$output), "fdr_curve.csv")
write.csv(fdr_table, fdr_table_path, row.names = FALSE)
message("📈 Saved FDR threshold table to: ", fdr_table_path)
                  
# ---- 3. Assign per-gene FDR ----
get_fdr_for_beta <- function(b) {
  idx <- which.min(abs(fdr_table$threshold - abs(b)))
  fdr_table$FDR[idx]
}
df_obs$FDR <- sapply(df_obs$estimate, get_fdr_for_beta)

# ---- 4. Diagnostic plots ----
plot_dir <- file.path(dirname(opt$output), "plots")
dir.create(plot_dir, showWarnings = FALSE)

# === 4a. Density: observed vs permuted β(TMB) ===
plot_path_density <- file.path(plot_dir, basename(sub("\\.csv$", "_density.pdf", opt$output)))
message("📊 Saving density plot to: ", plot_path_density)

pdf(plot_path_density, width = 6, height = 4)
ggplot() +
  geom_density(aes(x = df_perm$estimate), fill = "gray70", alpha = 0.6) +
  geom_density(aes(x = df_obs$estimate), color = "red", size = 1) +
  theme_bw() +
  labs(
    title = expression("Observed vs Permuted " * beta * "(TMB) Distribution"),
    x = expression(beta~"(TMB)"),
    y = "Density"
  )
dev.off()

# === 4b. FDR vs |β(TMB)| threshold ===
plot_path_fdr <- file.path(plot_dir, basename(sub("\\.csv$", "_fdr_curve.pdf", opt$output)))
message("📈 Saving empirical FDR curve to: ", plot_path_fdr)

pdf(plot_path_fdr, width = 6, height = 4)
ggplot(fdr_table, aes(x = threshold, y = FDR)) +
  geom_line(color = "blue", size = 1) +
  theme_bw() +
  labs(
    title = expression("Empirical FDR vs |" * beta * "(TMB)| Threshold"),
    x = expression("|" * beta * "(TMB)| threshold"),
    y = "Empirical FDR"
  )
dev.off()

# === 4c. Cumulative proportion (observed vs permuted) ===
plot_path_prop <- file.path(plot_dir, basename(sub("\\.csv$", "_prop_curves.pdf", opt$output)))
message("📉 Saving cumulative proportion plot to: ", plot_path_prop)

library(tidyr)
fdr_plot_long <- fdr_table %>%
  select(threshold, prop_obs, prop_perm) %>%
  pivot_longer(cols = c(prop_obs, prop_perm),
               names_to = "type", values_to = "proportion")

library(Cairo)
CairoPDF(plot_path_prop, width = 6, height = 4)
ggplot(fdr_plot_long, aes(x = threshold, y = proportion, color = type)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(
      title = "Cumulative Proportion ≥ |β(TMB)| Threshold",
      x = "|β(TMB)| threshold",
      y = "Proportion ≥ threshold",
      color = "Distribution"
)
dev.off()

# ---- 5. Write output ----
write.csv(df_obs, opt$output, row.names = FALSE)
message("✅ Saved empirical FDR results to: ", opt$output)
