# =============================================================================
# block_length_sensitivity_AR1.R
#
# Sensitivity of WRC, SPA, and CPA results to MBB block length.
# AR_1 is the forecast-comparison benchmark; HA is the realised series.
#
# Outputs:
#   reproducibility/output/block_length_sensitivity_AR1_results_raw.csv
#   reproducibility/output/block_length_sensitivity_AR1_results.csv
# =============================================================================

library(RCtest)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

forecast_model_cols <- seq_len(14L)

realized_col <- match("HA", colnames(metals))

benchmark_name <- "AR_1"
benchmark_idx <- match(benchmark_name, colnames(metals))

if (is.na(realized_col)) {
  stop("HA was not found as the realised series.")
}

if (is.na(benchmark_idx) || !benchmark_idx %in% forecast_model_cols) {
  stop("AR_1 must be one of the first 14 forecast-model columns.")
}

competitor_idx <- setdiff(forecast_model_cols, benchmark_idx)
competitor_names <- colnames(metals)[competitor_idx]

if (length(competitor_idx) != 13L) {
  stop("Expected 13 competing forecast models after excluding AR_1.")
}

realized <- metals[, realized_col]

benchmark_error <- metals[, benchmark_idx] - realized

competitor_errors <- sweep(
  metals[, competitor_idx, drop = FALSE],
  1,
  realized,
  FUN = "-"
)

# Positive values indicate that the competing model has lower MSE than AR_1.
loss_diff_mse <- benchmark_error^2 - competitor_errors^2
colnames(loss_diff_mse) <- competitor_names

# CPA conditioning variable: realised magnitude.
weighting_vector <- abs(realized)

block_lengths <- c(3L, 5L, 8L, 12L)
n_boot <- 999L
alpha <- 0.05

set.seed(20260822)

sensitivity_results <- lapply(block_lengths, function(block_length) {
  set.seed(20260822 + block_length)
  
  wrc <- white_reality_check(
    loss_differences = loss_diff_mse,
    n_simulations = n_boot,
    block_length = block_length,
    alpha = alpha
  )
  
  spa <- superior_predictive_ability_test(
    loss_differences = loss_diff_mse,
    block_length = block_length,
    num_bootstrap_replications = n_boot,
    alpha = alpha
  )
  
  cpa <- white_reality_check_conditional(
    loss_differences = loss_diff_mse,
    weighting_vector = weighting_vector,
    block_length = block_length,
    num_bootstrap_replications = n_boot,
    alpha = alpha
  )
  
  data.frame(
    Block_Length = block_length,
    WRC_Statistic = unname(wrc$statistic),
    WRC_p_value = wrc$p.value,
    SPA_Statistic = unname(spa$statistic),
    SPA_Consistent_p_value = spa$p_consistent,
    SPA_Conservative_p_value = spa$p_conservative,
    CPA_Statistic = unname(cpa$statistic),
    CPA_p_value = cpa$p.value,
    stringsAsFactors = FALSE
  )
})

sensitivity_results <- do.call(rbind, sensitivity_results)

# Save raw numeric results for auditability and later recalculation.
write.csv(
  sensitivity_results,
  file.path(
    output_dir,
    "block_length_sensitivity_AR1_results_raw.csv"
  ),
  row.names = FALSE
)

format_bootstrap_p <- function(p, n_boot) {
  lower_bound <- 1 / (n_boot + 1)
  upper_bound <- 1 - lower_bound
  
  if (is.na(p)) {
    return(NA_character_)
  }
  
  if (p <= lower_bound) {
    return(sprintf("< %.3f", lower_bound))
  }
  
  if (p >= upper_bound) {
    return(sprintf(">= %.3f", upper_bound))
  }
  
  sprintf("%.3f", p)
}

# Manuscript-ready table: p-values are replaced directly with formatted values.
manuscript_results <- sensitivity_results

p_columns <- grep(
  "_p_value$",
  names(manuscript_results),
  value = TRUE
)

for (column_name in p_columns) {
  manuscript_results[[column_name]] <- vapply(
    manuscript_results[[column_name]],
    format_bootstrap_p,
    character(1),
    n_boot = n_boot
  )
}

statistic_columns <- grep(
  "_Statistic$",
  names(manuscript_results),
  value = TRUE
)

for (column_name in statistic_columns) {
  manuscript_results[[column_name]] <- sprintf(
    "%.4f",
    manuscript_results[[column_name]]
  )
}

print(manuscript_results, row.names = FALSE)

write.csv(
  manuscript_results,
  file.path(
    output_dir,
    "block_length_sensitivity_AR1_results.csv"
  ),
  row.names = FALSE
)

message(
  "Saved:\n",
  "  reproducibility/output/block_length_sensitivity_AR1_results_raw.csv\n",
  "  reproducibility/output/block_length_sensitivity_AR1_results.csv\n"
)