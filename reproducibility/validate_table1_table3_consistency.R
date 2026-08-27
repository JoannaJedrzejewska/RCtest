# =============================================================================
# Automated check for Table 1 and Table 3
# Output:
#   reproducibility/output/table1_table3_consistency_check.csv
#   Console warning() if the check fails.
# =============================================================================
output_dir <- file.path("reproducibility", "output")

table1_path <- file.path(output_dir, "table1_AR1_raw.csv")
table3_path <- file.path(output_dir, "table3_AR1_raw.csv")

table3_mae_path <- file.path(output_dir, "table3_AR1_MAE_raw.csv")
table3_mase_path <- file.path(output_dir, "table3_AR1_MASE_raw.csv")

if (!file.exists(table1_path) || !file.exists(table3_path)) {
  stop(
    "Run table1_AR1.R and table3_AR1.R first -- missing: ",
    paste(c(table1_path, table3_path)[!file.exists(c(table1_path, table3_path))],
          collapse = ", ")
  )
}

table1_raw <- read.csv(table1_path, stringsAsFactors = FALSE)
table3_mse_raw <- read.csv(table3_path, stringsAsFactors = FALSE)

alpha <- 0.05
MAJORITY_THRESHOLD <- 0.5
WRC_SPA_SUSPICIOUS_THRESHOLD <- 0.99

# -----------------------------------------------------------------------------
# Helper: summarize a Table 3 per-model data frame.
# -----------------------------------------------------------------------------

summarize_table3 <- function(table3_df) {
  p_value_col <- grep("^P_Value", names(table3_df), value = TRUE)[1]
  if (is.na(p_value_col)) {
    stop("Could not locate a P_Value column in the supplied Table 3 data frame.")
  }
  
  frac_positive <- mean(table3_df$Mean_Loss_Diff > 0, na.rm = TRUE)
  frac_significant <- mean(table3_df[[p_value_col]] <= alpha, na.rm = TRUE)
  frac_significant_and_positive <- mean(
    table3_df$Mean_Loss_Diff > 0 & table3_df[[p_value_col]] <= alpha,
    na.rm = TRUE
  )
  
  list(
    frac_positive = frac_positive,
    frac_significant = frac_significant,
    frac_significant_and_positive = frac_significant_and_positive
  )
}

# -----------------------------------------------------------------------------
# Resolve which Table 3 source to use for each loss measure, and note
# whether it's a genuine measure-specific comparison or a fallback.
# -----------------------------------------------------------------------------

loss_measures <- c("MSE", "MAE", "MASE")

table3_sources <- list(
  MSE = list(data = table3_mse_raw, is_fallback = FALSE),
  MAE = if (file.exists(table3_mae_path)) {
    list(data = read.csv(table3_mae_path, stringsAsFactors = FALSE), is_fallback = FALSE)
  } else {
    list(data = table3_mse_raw, is_fallback = TRUE)
  },
  MASE = if (file.exists(table3_mase_path)) {
    list(data = read.csv(table3_mase_path, stringsAsFactors = FALSE), is_fallback = FALSE)
  } else {
    list(data = table3_mse_raw, is_fallback = TRUE)
  }
)

# -----------------------------------------------------------------------------
# Step through each loss measure and check Table 1 vs. the corresponding
# Table 3 summary.
# -----------------------------------------------------------------------------

results_list <- lapply(loss_measures, function(measure) {
  wrc_row <- table1_raw[table1_raw$Test == paste("WRC", measure), ]
  spa_row <- table1_raw[table1_raw$Test == paste("SPA", measure, "consistent"), ]
  cpa_row <- table1_raw[table1_raw$Test == paste("CPA", measure), ]
  
  if (nrow(wrc_row) == 0 || nrow(spa_row) == 0) {
    warning(sprintf(
      "Could not locate 'WRC %s' / 'SPA %s consistent' rows in table1_AR1_raw.csv -- skipping.",
      measure, measure
    ))
    return(NULL)
  }
  
  wrc_p <- wrc_row$P_value[1]
  spa_p <- spa_row$P_value[1]
  cpa_p <- if (nrow(cpa_row) > 0) cpa_row$P_value[1] else NA_real_
  
  table3_info <- table3_sources[[measure]]
  summary_stats <- summarize_table3(table3_info$data)
  
  contradiction_detected <- (summary_stats$frac_significant_and_positive > MAJORITY_THRESHOLD) &&
    (wrc_p >= WRC_SPA_SUSPICIOUS_THRESHOLD || spa_p >= WRC_SPA_SUSPICIOUS_THRESHOLD)
  
  data.frame(
    Loss_Measure = measure,
    Table3_Source = if (table3_info$is_fallback) {
      "FALLBACK: MSE-based table3_AR1_raw.csv (measure-specific file not found)"
    } else {
      "Measure-specific Table 3 file"
    },
    Frac_Models_Positive_And_Significant = round(summary_stats$frac_significant_and_positive, 4),
    WRC_P_Value = round(wrc_p, 4),
    SPA_P_Value_Consistent = round(spa_p, 4),
    CPA_P_Value = round(cpa_p, 4),
    Contradiction_Detected = contradiction_detected,
    stringsAsFactors = FALSE
  )
})

consistency_summary <- do.call(rbind, results_list[!sapply(results_list, is.null)])

cat("\n=== TABLE 1 / TABLE 3 CONSISTENCY CHECK (ALL LOSS MEASURES) ===\n")
cat(sprintf(
  "Rule: FAIL if >%.0f%%%% of models are individually significant AND superior,\n",
  MAJORITY_THRESHOLD * 100
))
cat(sprintf("      while the corresponding WRC/SPA p-value >= %.2f\n\n", WRC_SPA_SUSPICIOUS_THRESHOLD))
print(consistency_summary, row.names = FALSE)

any_fallback <- any(grepl("FALLBACK", consistency_summary$Table3_Source))
if (any_fallback) {
  cat("\nNOTE: one or more loss measures used the MSE-based Table 3 as a\n")
  cat("fallback because a measure-specific per-model DM file was not found.\n")
  cat("This means the MAE/MASE row(s) above check Table 1's MAE/MASE joint\n")
  cat("p-values against MSE-based per-model significance, NOT true MAE/MASE\n")
  cat("per-model results. Generate table3_AR1_MAE_raw.csv and\n")
  cat("table3_AR1_MASE_raw.csv (by rerunning compute_per_model_statistics()\n")
  cat("on loss_diff_mae / loss_diff_mase, analogous to table3_AR1.R) for a\n")
  cat("fully independent check on each measure.\n")
}

write.csv(
  consistency_summary,
  file.path(output_dir, "table1_table3_consistency_check.csv"),
  row.names = FALSE
)

n_contradictions <- sum(consistency_summary$Contradiction_Detected)

if (n_contradictions > 0L) {
  warning(
    n_contradictions, " of ", nrow(consistency_summary),
    " loss measure(s) still show the Major Comment 1 contradiction pattern: ",
    "a majority of models individually significant and superior, while the ",
    "joint WRC/SPA p-value sits at or above ", WRC_SPA_SUSPICIOUS_THRESHOLD, ". ",
    "See the Loss_Measure column in table1_table3_consistency_check.csv."
  )
} else {
  message(
    "No contradiction detected for any of MSE, MAE, or MASE. ",
    if (any_fallback) "(Some measures used a fallback Table 3 source -- see note above.)" else ""
  )
}
