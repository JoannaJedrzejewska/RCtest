# =============================================================================
# VaR backtesting example using RCtest's compute_kupiec(), with a direct
# numerical cross-check against:
#   (a) a hand-coded closed-form Kupiec LR statistic, and
#   (b) ExactVaRTest::backtest_lr(x, alpha, type = "uc"), which returns the
#       EXACT finite-sample p-value via dynamic-programming enumeration
#       rather than the chi-square asymptotic approximation.
#
# Outputs:
#   reproducibility/output/var_kupiec_backtest_results.csv
#   reproducibility/output/var_kupiec_validation_against_ExactVaRTest.csv
#   reproducibility/output/var_kupiec_validation_against_ExactVaRTest_manuscript.csv
#   reproducibility/output/var_kupiec_duplicate_statistic_diagnostic.csv
# =============================================================================

library(RCtest)

if (!requireNamespace("ExactVaRTest", quietly = TRUE)) {
  stop("Package 'ExactVaRTest' is required for this validation script.")
}

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260824)

forecast_model_cols <- seq_len(14L)
realized_col <- match("HA", colnames(metals))
benchmark_name <- "AR_1"
benchmark_col_full <- match(benchmark_name, colnames(metals))

if (is.na(realized_col)) {
  stop("HA was not found as the realized series.")
}

if (is.na(benchmark_col_full) || !benchmark_col_full %in% forecast_model_cols) {
  stop("AR_1 must be one of the first 14 forecast-model columns.")
}

forecast_matrix <- metals[, forecast_model_cols, drop = FALSE]
realized <- metals[, realized_col]
benchmark_col <- match(benchmark_name, colnames(forecast_matrix))
competitor_names <- setdiff(colnames(forecast_matrix), benchmark_name)

if (length(competitor_names) != 13L) {
  stop("Expected 13 competing forecast models after excluding AR_1.")
}

alpha_var <- 0.05
window_size <- 20L

forecast_variance <- estimate_forecast_variance(
  forecast_matrix = metals,
  realized = realized,
  benchmark_col = benchmark_col_full,
  window_size = window_size
)

sd_model_cols <- match(competitor_names, colnames(metals))
forecast_sd_models <- sqrt(forecast_variance[, sd_model_cols, drop = FALSE])
colnames(forecast_sd_models) <- competitor_names

if (ncol(forecast_sd_models) != ncol(forecast_matrix) - 1L) {
  stop("forecast_sd_models must contain exactly one column per competing model.")
}

rctest_kupiec <- compute_kupiec(
  forecast_matrix = forecast_matrix,
  forecast_sd_models = forecast_sd_models,
  realized = realized,
  benchmark_col = benchmark_col,
  alpha = alpha_var
)

rctest_kupiec_table <- do.call(
  rbind,
  lapply(names(rctest_kupiec), function(model_name) {
    result <- rctest_kupiec[[model_name]]
    
    data.frame(
      Model = model_name,
      RCtest_LR_Statistic = unname(result$statistic),
      RCtest_P_Value = result$p.value,
      RCtest_Actual_Exceedances = result$actual_exceedances,
      RCtest_Expected_Exceedances = result$expected,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  rctest_kupiec_table,
  file.path(output_dir, "var_kupiec_backtest_results.csv"),
  row.names = FALSE
)

hand_coded_kupiec_lr <- function(violations, alpha) {
  n <- length(violations)
  x <- sum(violations)
  p_hat <- x / n
  
  if (p_hat == 0 || p_hat == 1) {
    return(list(
      LR = NA_real_,
      p_value = NA_real_,
      x = x,
      n = n,
      p_hat = p_hat
    ))
  }
  
  log_lik_null <- (n - x) * log(1 - alpha) + x * log(alpha)
  log_lik_alt <- (n - x) * log(1 - p_hat) + x * log(p_hat)
  LR <- -2 * (log_lik_null - log_lik_alt)
  
  list(
    LR = LR,
    p_value = 1 - pchisq(LR, df = 1),
    x = x,
    n = n,
    p_hat = p_hat
  )
}

z_alpha <- qnorm(alpha_var)

cross_check_results <- lapply(competitor_names, function(model_name) {
  model_col <- match(model_name, colnames(forecast_matrix))
  model_sd <- forecast_sd_models[, model_name]
  
  valid_idx <- which(!is.na(model_sd) & is.finite(model_sd) & model_sd > 0)
  if (length(valid_idx) < 30L) {
    return(NULL)
  }
  
  forecast_values <- forecast_matrix[valid_idx, model_col]
  realized_values <- realized[valid_idx]
  sd_values <- model_sd[valid_idx]
  
  var_estimate <- forecast_values + z_alpha * sd_values
  violations <- as.integer(realized_values < var_estimate)
  
  hand_result <- hand_coded_kupiec_lr(violations, alpha_var)
  
  exact_result <- ExactVaRTest::backtest_lr(
    violations,
    alpha = alpha_var,
    type = "uc"
  )
  
  rctest_row <- rctest_kupiec_table[
    rctest_kupiec_table$Model == model_name,
    ,
    drop = FALSE
  ]
  
  data.frame(
    Model = model_name,
    N_Obs = hand_result$n,
    N_Violations = hand_result$x,
    Empirical_Violation_Rate = hand_result$p_hat,
    Nominal_Rate = alpha_var,
    RCtest_LR = rctest_row$RCtest_LR_Statistic,
    RCtest_P_Value_Asymptotic = rctest_row$RCtest_P_Value,
    HandCoded_LR = hand_result$LR,
    HandCoded_P_Value_Asymptotic = hand_result$p_value,
    ExactVaRTest_LR = unname(exact_result$stat),
    ExactVaRTest_P_Value_Exact = exact_result$pval,
    stringsAsFactors = FALSE
  )
})

cross_check_table <- do.call(
  rbind,
  cross_check_results[!vapply(cross_check_results, is.null, logical(1))]
)

write.csv(
  cross_check_table,
  file.path(output_dir, "var_kupiec_validation_against_ExactVaRTest.csv"),
  row.names = FALSE
)

lr_diff_handcoded <- abs(
  cross_check_table$RCtest_LR - cross_check_table$HandCoded_LR
)

lr_diff_exact <- abs(
  cross_check_table$HandCoded_LR - cross_check_table$ExactVaRTest_LR
)

max_diff_handcoded <- max(lr_diff_handcoded, na.rm = TRUE)
max_diff_exact <- max(lr_diff_exact, na.rm = TRUE)

cat("\n=== KUPIEC CROSS-CHECK ===\n")
print(cross_check_table, row.names = FALSE, digits = 6)

cat("\n=== AGREEMENT SUMMARY ===\n")
cat(sprintf("Max |RCtest LR - hand-coded LR|:      %.12g\n", max_diff_handcoded))
cat(sprintf("Max |hand-coded LR - ExactVaRTest LR|: %.12g\n", max_diff_exact))

numerical_tolerance <- 1e-10

if (max_diff_handcoded > numerical_tolerance) {
  stop(
    "Validation failed: RCtest and hand-coded Kupiec LR statistics differ by ",
    format(max_diff_handcoded, scientific = TRUE),
    ", exceeding tolerance ", format(numerical_tolerance, scientific = TRUE), "."
  )
}

if (max_diff_exact > numerical_tolerance) {
  stop(
    "Validation failed: hand-coded and ExactVaRTest Kupiec LR statistics differ by ",
    format(max_diff_exact, scientific = TRUE),
    ", exceeding tolerance ", format(numerical_tolerance, scientific = TRUE), "."
  )
}

message(
  "Validation passed: RCtest, the hand-coded formula, and ExactVaRTest agree ",
  "on the LR-UC statistic to numerical precision."
)

rounded_stat <- round(cross_check_table$RCtest_LR, 12)
duplicate_groups <- split(cross_check_table$Model, rounded_stat)
duplicate_groups <- duplicate_groups[vapply(duplicate_groups, length, integer(1)) > 1L]

duplicate_diagnostic <- if (length(duplicate_groups) > 0L) {
  do.call(rbind, lapply(names(duplicate_groups), function(stat_key) {
    models <- duplicate_groups[[stat_key]]
    rows <- cross_check_table[match(models, cross_check_table$Model), , drop = FALSE]
    
    data.frame(
      Shared_LR_Statistic = as.numeric(stat_key),
      Models = paste(models, collapse = "; "),
      N_Models = length(models),
      Shared_Violation_Count = length(unique(rows$N_Violations)) == 1L,
      Violation_Counts = paste(rows$N_Violations, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(
    Shared_LR_Statistic = numeric(0),
    Models = character(0),
    N_Models = integer(0),
    Shared_Violation_Count = logical(0),
    Violation_Counts = character(0)
  )
}

write.csv(
  duplicate_diagnostic,
  file.path(output_dir, "var_kupiec_duplicate_statistic_diagnostic.csv"),
  row.names = FALSE
)

if (nrow(duplicate_diagnostic) > 0L && any(!duplicate_diagnostic$Shared_Violation_Count)) {
  warning(
    "Some models share an LR statistic but not a violation count. Review ",
    "var_kupiec_duplicate_statistic_diagnostic.csv."
  )
}

# -----------------------------------------------------------------------------
# Manuscript-format output: values are formatted ONLY after numerical
# validation, preserving full precision in the raw validation CSV.
# -----------------------------------------------------------------------------

cross_check_manuscript <- cross_check_table

integer_columns <- c("N_Obs", "N_Violations")
for (column_name in integer_columns) {
  cross_check_manuscript[[column_name]] <- as.integer(
    cross_check_manuscript[[column_name]]
  )
}

numeric_columns <- c(
  "Empirical_Violation_Rate",
  "Nominal_Rate",
  "RCtest_LR",
  "RCtest_P_Value_Asymptotic",
  "HandCoded_LR",
  "HandCoded_P_Value_Asymptotic",
  "ExactVaRTest_LR",
  "ExactVaRTest_P_Value_Exact"
)

for (column_name in numeric_columns) {
  cross_check_manuscript[[column_name]] <- sprintf(
    "%.4f",
    cross_check_manuscript[[column_name]]
  )
}

write.csv(
  cross_check_manuscript,
  file.path(
    output_dir,
    "var_kupiec_validation_against_ExactVaRTest_manuscript.csv"
  ),
  row.names = FALSE
)