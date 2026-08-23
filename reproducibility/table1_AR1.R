# =============================================================================
# Table 1: Joint forecast-comparison tests using AR_1 as benchmark.
#
# HA is the realised outcome series.
# The comparison includes the 13 forecast methods other than AR_1.
#
# Outputs:
#   reproducibility/output/table1_AR1_raw.csv
#   reproducibility/output/table1_AR1.csv
# =============================================================================

library(RCtest)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260823)

forecast_model_cols <- seq_len(14L)
realized_col <- match("HA", colnames(metals))

benchmark_name <- "AR_1"
benchmark_col <- match(benchmark_name, colnames(metals))

if (is.na(realized_col)) {
  stop("HA was not found as the realised outcome series.")
}

if (is.na(benchmark_col) || !benchmark_col %in% forecast_model_cols) {
  stop("AR_1 must be one of the first 14 forecast-model columns.")
}

competitor_cols <- setdiff(forecast_model_cols, benchmark_col)
competitor_names <- colnames(metals)[competitor_cols]

if (length(competitor_cols) != 13L) {
  stop("Expected 13 competing models after excluding AR_1.")
}

realized <- metals[, realized_col]

benchmark_error <- metals[, benchmark_col] - realized

competitor_errors <- sweep(
  metals[, competitor_cols, drop = FALSE],
  1,
  realized,
  FUN = "-"
)

# Positive loss differentials indicate that the competitor beat AR_1.
loss_diff_mse <- benchmark_error^2 - competitor_errors^2

loss_diff_mae <- abs(benchmark_error) - abs(competitor_errors)

mase_scale <- mean(abs(diff(realized)), na.rm = TRUE)

loss_diff_mase <- loss_diff_mae / mase_scale

colnames(loss_diff_mse) <- competitor_names
colnames(loss_diff_mae) <- competitor_names
colnames(loss_diff_mase) <- competitor_names

block_length <- 5L
n_boot <- 999L
alpha <- 0.05

# -----------------------------------------------------------------------------
# Predictive-density inputs
#
# estimate_forecast_variance() returns a P x 15 matrix. The SD matrix passed
# into compute_klic() and compute_zp() must contain every column except AR_1,
# in exactly that order, because AR_1 is their designated benchmark column.
# -----------------------------------------------------------------------------

forecast_variance <- estimate_forecast_variance(
  forecast_matrix = metals,
  realized = realized,
  benchmark_col = benchmark_col,
  window_size = 20
)

sd_model_cols <- setdiff(seq_len(ncol(metals)), benchmark_col)

forecast_sd_models <- sqrt(
  forecast_variance[, sd_model_cols, drop = FALSE]
)

if (ncol(forecast_sd_models) != ncol(metals) - 1L) {
  stop("forecast_sd_models does not contain one SD column per non-benchmark series.")
}

threshold_value <- quantile(realized, probs = 0.05, na.rm = TRUE)

klic_loss <- compute_klic(
  forecast_matrix = metals,
  forecast_sd_models = forecast_sd_models,
  benchmark_col = benchmark_col
)

zp_loss <- compute_zp(
  forecast_matrix = metals,
  forecast_sd_models = forecast_sd_models,
  threshold = threshold_value,
  benchmark_col = benchmark_col
)

# Keep only the 13 forecast competitors; don't include HA as a competitor.
klic_diff <- klic_loss[, benchmark_col] -
  klic_loss[, competitor_cols, drop = FALSE]

zp_diff <- zp_loss[, benchmark_col] -
  zp_loss[, competitor_cols, drop = FALSE]

colnames(klic_diff) <- competitor_names
colnames(zp_diff) <- competitor_names

run_loss_tests <- function(loss_differences, label) {
  wrc <- white_reality_check(
    loss_differences = loss_differences,
    n_simulations = n_boot,
    block_length = block_length,
    alpha = alpha
  )
  
  spa <- superior_predictive_ability_test(
    loss_differences = loss_differences,
    block_length = block_length,
    num_bootstrap_replications = n_boot,
    alpha = alpha
  )
  
  cpa <- white_reality_check_conditional(
    loss_differences = loss_differences,
    weighting_vector = abs(realized),
    block_length = block_length,
    num_bootstrap_replications = n_boot,
    alpha = alpha
  )
  
  data.frame(
    Test = c(
      paste("WRC", label),
      paste("SPA", label, "consistent"),
      paste("SPA", label, "conservative"),
      paste("CPA", label)
    ),
    Statistic = c(
      unname(wrc$statistic),
      unname(spa$statistic),
      unname(spa$statistic),
      unname(cpa$statistic)
    ),
    P_value = c(
      wrc$p.value,
      spa$p_consistent,
      spa$p_conservative,
      cpa$p.value
    ),
    stringsAsFactors = FALSE
  )
}

results_mse <- run_loss_tests(loss_diff_mse, "MSE")
results_mae <- run_loss_tests(loss_diff_mae, "MAE")
results_mase <- run_loss_tests(loss_diff_mase, "MASE")

# -----------------------------------------------------------------------------
# Distributional joint tests.
# -----------------------------------------------------------------------------

zp_test <- reality_check_zp_test(
  zp_loss_differences = zp_diff,
  block_length = block_length,
  num_bootstrap_replications = n_boot,
  alpha = alpha
)

klic_test <- kullback_leibler_test(
  log_likelihood_differences = klic_diff,
  block_length = block_length,
  num_bootstrap_replications = n_boot,
  alpha = alpha
)

cdf_rc_test <- white_reality_check_cdf_approx(
  loss_differences = loss_diff_mse,
  block_length = block_length,
  num_bootstrap_replications = n_boot,
  alpha = alpha
)

results_distributional <- data.frame(
  Test = c("ZP", "KLIC", "CDF-RC"),
  Statistic = c(
    unname(zp_test$statistic),
    unname(klic_test$statistic),
    unname(cdf_rc_test$statistic)
  ),
  P_value = c(
    zp_test$p_consistent,
    klic_test$p.value,
    cdf_rc_test$p.value
  ),
  stringsAsFactors = FALSE
)

table1_raw <- rbind(
  results_mse,
  results_mae,
  results_mase,
  results_distributional
)

table1_raw$Decision <- ifelse(
  table1_raw$P_value <= alpha,
  "Reject H0",
  "Fail to reject H0"
)

write.csv(
  table1_raw,
  file.path(output_dir, "table1_AR1_raw.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Manuscript formatting.
# -----------------------------------------------------------------------------

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

table1_manuscript <- table1_raw

table1_manuscript$Statistic <- sprintf(
  "%.4f",
  table1_manuscript$Statistic
)

table1_manuscript$P_value <- vapply(
  table1_manuscript$P_value,
  format_bootstrap_p,
  character(1),
  n_boot = n_boot
)

print(table1_manuscript, row.names = FALSE)

write.csv(
  table1_manuscript,
  file.path(output_dir, "table1_AR1.csv"),
  row.names = FALSE
)
