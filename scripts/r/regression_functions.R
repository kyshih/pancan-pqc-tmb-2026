# ============================================================
# Regression functions for TCGA Pan-Cancer analyses
#
# 1. run_global_regression():  mixed-effects model across all genes
# 2. run_genewise_regression(): mixed-effects model per gene
#
# Both output coefficient tables with confidence intervals
# ============================================================

library(lme4)
library(dplyr)
library(tibble)
library(lmerTest)
# library(broom.mixed)

# ------------------------------------------------------------
# 1️⃣ Global regression: overall effect of TMB on expression
# ------------------------------------------------------------
# run_global_regression <- function(df, output_file) {
#   message("Running global regression of log_mut on z-scored expression...")

#   # Ensure expected columns
#   stopifnot(all(c("expr_z", "log_mut", "purity", "cancer_type", "gene") %in% colnames(df)))

#   # Recompute global z-score across all genes
#   # df <- df %>%
#   #   mutate(expr_z_global = scale(expr_z, center = TRUE, scale = TRUE)[, 1])

#   # Fit model: global average effect
#   # model <- lmer(expr_z_global ~ log_mut + purity + (1 | cancer_type) + (1 | gene), data = df)
#   model <- lmer(expr_z ~ log_mut + purity + (1 | cancer_type) + (1 | gene), data = df)

#   # Extract summary + CIs
#   coefs <- as.data.frame(coef(summary(model)))
#   ci <- confint(model, method = "Wald")[c("(Intercept)", "log_mut", "purity"), ] %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term")

#   results <- coefs %>%
#     tibble::rownames_to_column("term") %>%
#     left_join(ci, by = "term") %>%
#     mutate(gene = "Global")

#   write.table(results, file = output_file, sep = ",",
#               row.names = FALSE, col.names = !file.exists(output_file), append = TRUE)

#   message("✅ Global regression complete.")
#   return(results)
# }

run_global_regression <- function(df, output_file, residuals_file = NULL, residual_genes = NULL) {
  message("Running global regression of log_mut on per-gene z-scored expression...")

  # Ensure expected columns
  stopifnot(all(c("expr_z", "log_mut", "purity", "cancer_type", "gene") %in% colnames(df)))

  # ---- Fit model: global mixed effects ----
  # Use per-gene z-scored expression (expr_z)
  # Random intercepts for both gene and cancer_type
  model <- lmerTest::lmer(expr_z ~ log_mut + purity + (1 | cancer_type) + (1 | gene), data = df)

  # ---- Extract fixed-effect coefficients (with p-values and SEs) ----
  coefs <- as.data.frame(summary(model)$coefficients) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::rename(
      estimate  = Estimate,
      std.error = `Std. Error`,
      t.value   = `t value`,
      p.value   = `Pr(>|t|)`
    )

  # ---- Add 95% confidence intervals ----
  ci <- confint(model, method = "Wald")[c("(Intercept)", "log_mut", "purity"), , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::rename(conf.low = `2.5 %`, conf.high = `97.5 %`)

  # ---- Merge results ----
  results <- dplyr::left_join(coefs, ci, by = "term") %>%
    dplyr::mutate(gene = "Global") %>%
    dplyr::select(gene, term, estimate, std.error, t.value, p.value, conf.low, conf.high)

  # ---- Write results ----
  write.table(
    results,
    file = output_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(output_file),
    append = TRUE
  )

  # ---- Optional: residual summaries for selected genes ----
  if (!is.null(residuals_file) && !is.null(residual_genes)) {
    message("Computing residual summaries for selected genes...")

    resid_df <- data.frame(
      barcode = df$barcode,
      cancer_type = df$cancer_type,
      gene = df$gene,
      expr_z = df$expr_z,
      fitted = fitted(model),
      resid = resid(model)
    )

    # Filter to selected genes
    resid_sel <- resid_df %>%
      dplyr::filter(gene %in% residual_genes)

    # Summarize residuals per cancer type
    resid_summary <- resid_sel %>%
      dplyr::group_by(gene, cancer_type) %>%
      dplyr::summarise(
        mean_resid = mean(resid, na.rm = TRUE),
        mean_abs_resid = mean(abs(resid), na.rm = TRUE),
        sd_resid = sd(resid, na.rm = TRUE),
        .groups = "drop"
      )

    # Save summary
    data.table::fwrite(resid_summary, residuals_file)
    message(paste("Residual summaries written to:", residuals_file))
  }

  message("✅ Global regression complete.")
  return(results)
}


# ------------------------------------------------------------
# 2️⃣ Per-gene regression: independent model for each gene
# ------------------------------------------------------------
# run_genewise_regression <- function(df, genes = "all", output_file) {
#   message("Running per-gene regressions...")
#   stopifnot(all(c("expr_z", "log_mut", "purity", "cancer_type", "gene") %in% colnames(df)))

#   if (identical(genes, "all")) genes <- unique(df$gene) # skipp if provide genes

#   results_all <- list()

#   for (g in genes) {
#     message(paste0(" → ", g))
#     df_gene <- df %>% filter(gene == g)

#     if (nrow(df_gene) < 30) {
#       warning(paste("Skipping", g, "- too few samples (", nrow(df_gene), ")"))
#       next
#     }

#     tryCatch({
#       model <- lmer(expr_z ~ log_mut + purity + (1 | cancer_type), data = df_gene)

#       coefs <- as.data.frame(coef(summary(model)))
#       ci <- confint(model, method = "Wald")[c("(Intercept)", "log_mut", "purity"), ] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("term")

#       results <- coefs %>%
#         tibble::rownames_to_column("term") %>%
#         left_join(ci, by = "term") %>%
#         mutate(gene = g)

#       results_all[[g]] <- results

#       write.table(results, file = output_file, sep = ",",
#                   row.names = FALSE, col.names = !file.exists(output_file), append = TRUE)
#     },
#     error = function(e) message(paste("❌ Skipping", g, ":", e$message)))
#   }

#   message("✅ Finished per-gene regressions.")
#   return(bind_rows(results_all))
# }

run_genewise_regression <- function(df, genes = "all", output_file) {
  message("Running per-gene regressions (lmerTest + tidy output)...")

  stopifnot(all(c("expr_z", "log_mut", "purity", "cancer_type", "gene") %in% colnames(df)))

  # ---- Select which genes to run ----
  all_genes <- unique(df$gene)
  if (!identical(genes, "all")) {
    missing <- setdiff(genes, all_genes)
    if (length(missing) > 0)
      warning(sprintf("These genes not found and will be skipped: %s", paste(missing, collapse = ", ")))
    genes <- intersect(genes, all_genes)
  } else genes <- all_genes

  # ---- Identify already completed genes ----
  if (file.exists(output_file)) {
    done_genes <- tryCatch({
      unique(read.csv(output_file, colClasses = c(gene = "character"))$gene)
    }, error = function(e) character(0))
    genes <- setdiff(genes, done_genes)
    message(sprintf("Skipping %d completed genes already in %s",
                    length(done_genes), basename(output_file)))
  }

  message(sprintf("Running regression for %d remaining genes...", length(genes)))

  results_all <- list()
  skipped_file <- sub("\\.csv$", "_skipped.csv", output_file)

  for (g in genes) {
    message(sprintf(" → %s", g))
    df_gene <- dplyr::filter(df, gene == g)

    # ---- Skip small genes ----
    if (nrow(df_gene) < 30) {
      reason <- sprintf("too few samples (%d)", nrow(df_gene))
      warning(sprintf("Skipping %s — %s", g, reason))

      # Write skipped gene immediately
      append_skipped_genes(df,
        data.frame(gene = g, reason = reason),
        skipped_file
      )
      next
    }

    # ---- Model fitting ----
    tryCatch({
      model <- lmerTest::lmer(expr_z ~ log_mut + purity + (1 | cancer_type), data = df_gene)

      coefs <- as.data.frame(summary(model)$coefficients) %>%
        tibble::rownames_to_column("term") %>%
        dplyr::rename(
          estimate = Estimate, std.error = `Std. Error`,
          t.value = `t value`, p.value = `Pr(>|t|)`
        )

      ci <- confint(model, method = "Wald")[c("(Intercept)", "log_mut", "purity"), , drop = FALSE] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("term") %>%
        dplyr::rename(conf.low = `2.5 %`, conf.high = `97.5 %`)

      results <- dplyr::left_join(coefs, ci, by = "term") %>%
        dplyr::mutate(gene = g) %>%
        dplyr::select(gene, term, estimate, std.error, t.value, p.value, conf.low, conf.high)

      # Write regression result immediately
      write.table(results, file = output_file, sep = ",",
                  row.names = FALSE, col.names = !file.exists(output_file), append = TRUE)

      results_all[[g]] <- results
      message(sprintf("✅ Finished %s", g))
    },
    error = function(e) {
      reason <- paste("model error:", e$message)
      message(sprintf("❌ Skipping %s: %s", g, e$message))

      # Write skipped gene immediately
      append_skipped_genes(df,
        data.frame(gene = g, reason = reason),
        skipped_file
      )
    })
  }

  # ---- Final skipped summary (prints once) ----
if (file.exists(skipped_file)) {
  append_skipped_genes(df,
    read.csv(skipped_file, colClasses = c(gene = "character")),
    skipped_file,
    quiet = FALSE
  )
} else {
  message("✅ No genes were skipped.")
}
                           
  message("🏁 All per-gene regressions finished.")
  return(dplyr::bind_rows(results_all))
}
