# append_skipped_genes <- function(df, skipped_genes, skipped_file) {
#   if (nrow(skipped_genes) == 0) {
#     message("✅ No skipped genes to log.")
#     return(invisible(NULL))
#   }
    
#   # Ensure gene is always a character
#   skipped_genes$gene <- as.character(skipped_genes$gene)

#   # ---- Add n_samples and n_cancer_types (per skipped gene) ----
#   skipped_genes <- skipped_genes %>%
#     dplyr::mutate(gene = as.character(gene)) %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(
#       n_samples = sum(df$gene == gene),
#       n_cancer_types = dplyr::n_distinct(df$cancer_type[df$gene == gene])
#     ) %>%
#     dplyr::ungroup() %>%
#     dplyr::distinct()

#   # ---- Append only new genes ----
#   if (file.exists(skipped_file)) {
#     existing <- tryCatch(read.csv(skipped_file, colClasses = c(gene = "character")), error = function(e) data.frame()
# )
#     new <- dplyr::anti_join(skipped_genes, existing, by = "gene")

#     if (nrow(new) > 0) {
#       write.table(new, skipped_file, sep = ",", row.names = FALSE,
#                   col.names = FALSE, append = TRUE)
#       message(sprintf("⚠️ Appended %d new skipped genes to %s", nrow(new), skipped_file))
#     } else {
#       message("✅ No new skipped genes to append.")
#     }
#   } else {
#     write.table(skipped_genes, skipped_file, sep = ",", row.names = FALSE, col.names = TRUE)
#     message(sprintf("🆕 Created skipped gene log: %s", skipped_file))
#   }
# }

# append_skipped_genes <- function(df_long, skipped_genes, skipped_file) {
#   if (nrow(skipped_genes) == 0) {
#     message("✅ No skipped genes to log.")
#     return(invisible(NULL))
#   }

#   # Ensure character consistency
#   skipped_genes$gene <- as.character(skipped_genes$gene)
#   df_long$gene <- as.character(df_long$gene)

#   # ---- Compute compact per-gene stats ----
#   gene_stats <- df_long %>%
#     dplyr::group_by(gene) %>%
#     dplyr::summarise(
#       expr_mean = mean(expr_z, na.rm = TRUE),
#       expr_sd = sd(expr_z, na.rm = TRUE),
#       prop_na = mean(is.na(expr_z)),
#       prop_nonzero = mean(expr_z > 0, na.rm = TRUE),
#       # Correlation only if predictor exists
#       corr_with_log_mut = if ("log_mut" %in% names(df_long))
#         suppressWarnings(cor(expr_z, log_mut, use = "pairwise.complete.obs"))
#       else NA_real_,
#       .groups = "drop"
#     )

#   # ---- Merge with skipped genes ----
#   skipped_genes <- skipped_genes %>%
#     dplyr::left_join(gene_stats, by = "gene")

#   # ---- Append only new genes ----
#   if (file.exists(skipped_file)) {
#     existing <- tryCatch(
#       read.csv(skipped_file, colClasses = c(gene = "character")),
#       error = function(e) data.frame()
#     )

#     new <- dplyr::anti_join(skipped_genes, existing, by = "gene")

#     if (nrow(new) > 0) {
#       write.table(new, skipped_file, sep = ",", row.names = FALSE,
#                   col.names = FALSE, append = TRUE)
#       message(sprintf("⚠️ Appended %d new skipped genes to %s", nrow(new), skipped_file))
#     } else {
#       message("✅ No new skipped genes to append.")
#     }
#   } else {
#     write.table(skipped_genes, skipped_file, sep = ",", row.names = FALSE, col.names = TRUE)
#     message(sprintf("🆕 Created skipped gene log: %s", skipped_file))
#   }

#   # ---- Quick diagnostic summary ----
#   n_zero_sd <- sum(skipped_genes$expr_sd == 0, na.rm = TRUE)
#   n_all_na  <- sum(skipped_genes$prop_na == 1, na.rm = TRUE)
#   frac_zero_sd <- n_zero_sd / nrow(skipped_genes)
#   frac_all_na  <- n_all_na / nrow(skipped_genes)

#   message(sprintf(
#     "📊 Summary of skipped genes: %d total | %.1f%% zero variance | %.1f%% all NA",
#     nrow(skipped_genes), frac_zero_sd * 100, frac_all_na * 100
#   ))

#   if (frac_zero_sd > 0.1) {
#     message("⚠️ More than 10% of skipped genes had zero variance — likely biological, not code-related.")
#   }

#   invisible(skipped_genes)
# }

append_skipped_genes <- function(df_long, skipped_genes, skipped_file, quiet = TRUE) {
  if (nrow(skipped_genes) == 0) {
    if (!quiet) message("✅ No skipped genes to log.")
    return(invisible(NULL))
  }

  # Ensure character consistency
  skipped_genes$gene <- as.character(skipped_genes$gene)
  df_long$gene <- as.character(df_long$gene)

  # ---- Compute compact per-gene stats ----
  gene_stats <- df_long %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      expr_mean = mean(expr_z, na.rm = TRUE),
      expr_sd = sd(expr_z, na.rm = TRUE),
      prop_na = mean(is.na(expr_z)),
      prop_nonzero = mean(expr_z > 0, na.rm = TRUE),
      corr_with_log_mut = if ("log_mut" %in% names(df_long))
        suppressWarnings(cor(expr_z, log_mut, use = "pairwise.complete.obs"))
      else NA_real_,
      .groups = "drop"
    )

  # ---- Merge with skipped genes ----
  skipped_genes <- skipped_genes %>%
    dplyr::left_join(gene_stats, by = "gene")

  # ---- Append only new genes ----
  if (file.exists(skipped_file)) {
    existing <- tryCatch(
      read.csv(skipped_file, colClasses = c(gene = "character")),
      error = function(e) data.frame()
    )

    new <- dplyr::anti_join(skipped_genes, existing, by = "gene")

    if (nrow(new) > 0) {
      write.table(new, skipped_file, sep = ",", row.names = FALSE,
                  col.names = FALSE, append = TRUE)
      message(sprintf("⚠️ Appended %d new skipped genes to %s", nrow(new), skipped_file))
    } else {
      message("✅ No new skipped genes to append.")
    }
  } else {
    write.table(skipped_genes, skipped_file, sep = ",", row.names = FALSE, col.names = TRUE)
    message(sprintf("🆕 Created skipped gene log: %s", skipped_file))
  }

  # ---- Optional: Quick diagnostic summary (only if quiet=FALSE) ----
  if (!quiet) {
    skipped_full <- tryCatch(read.csv(skipped_file), error = function(e) skipped_genes)

    n_zero_sd <- sum(skipped_full$expr_sd == 0, na.rm = TRUE)
    n_all_na  <- sum(skipped_full$prop_na == 1, na.rm = TRUE)
    frac_zero_sd <- n_zero_sd / nrow(skipped_full)
    frac_all_na  <- n_all_na / nrow(skipped_full)

    message(sprintf(
      "📊 Summary of skipped genes: %d total | %.1f%% zero variance | %.1f%% all NA",
      nrow(skipped_full), frac_zero_sd * 100, frac_all_na * 100
    ))

    if (frac_zero_sd > 0.1) {
      message("⚠️ More than 10% of skipped genes had zero variance — likely biological, not code-related.")
    }
  }

  invisible(skipped_genes)
}

