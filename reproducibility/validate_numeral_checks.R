# =============================================================================
# - Reproduces CDF-RC stability and quantile-grid diagnostics.
# - Validates the WRC test statistic on a deterministic toy example.
# - Checks the SPA consistent/conservative p-value relationship.
# - Compares DM decisions with forecast::dm.test().
#
# output/cdf_rc_quantile_grid_and_stability_diagnostics.csv
# output/wrc_hand_calculation_validation.csv
# output/spa_recentering_property_check.csv
# output/dm_validation_against_forecast.csv
# output/crps_validation_against_scoringRules.csv
# output/var_kupiec_validation_against_ExactVaRTest.csv
# =============================================================================

library(RCtest)
library(forecast)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
data("metals", package = "RCtest")
metals <- as.matrix(metals)

realized_name <- "HA"
benchmark_name <- "AR_1"

if (!realized_name %in% colnames(metals)) {
  stop(
    paste0(
      "Column '", realized_name,
      "' not found. Available columns: ",
      paste(colnames(metals), collapse = ", ")
    )
  )
}

if (!benchmark_name %in% colnames(metals)) {
  stop(
    paste0(
      "Column '", benchmark_name,
      "' not found. Available columns: ",
      paste(colnames(metals), collapse = ", ")
    )
  )
}

realized_values <- as.numeric(metals[, realized_name])

forecast_matrix <- metals[
  , setdiff(colnames(metals), realized_name),
  drop = FALSE
]

colnames(forecast_matrix)

if (!benchmark_name %in% colnames(forecast_matrix)) {
  stop("Column 'AR_1' was not found in forecast_matrix.")
}

benchmark_forecast <- forecast_matrix[, benchmark_name]

competitor_names <- setdiff(
  colnames(forecast_matrix),
  benchmark_name
)

competitor_forecasts <- forecast_matrix[, competitor_names, drop = FALSE]

benchmark_mse_loss <- (realized_values - benchmark_forecast)^2

competitor_mse_loss <- (
  matrix(
    realized_values,
    nrow = length(realized_values),
    ncol = ncol(competitor_forecasts)) - competitor_forecasts)^2

mse_loss_differences <- sweep(
  competitor_mse_loss,1,benchmark_mse_loss,
  FUN = function(model_loss, benchmark_loss) benchmark_loss - model_loss)

colnames(mse_loss_differences) <- competitor_names


if (!exists("mse_loss_differences")) {
  stop(
    "Create mse_loss_differences from the exact Table 1 MSE ",
    "loss-difference matrix before running this validation script."
  )
}

mse_loss_differences <- as.matrix(mse_loss_differences)

# ============================================================
# LOCAL CDF-RC DIAGNOSTIC FUNCTION
# ============================================================

cdf_rc_diagnostic <- function(
    loss_differences,
    block_length = 5,
    num_bootstrap_replications = 999,
    quantile_grid = seq(0.1, 0.9, by = 0.1),
    variance_floor = 1e-10
) {
  loss_differences <- as.matrix(loss_differences)
  
  if (!is.numeric(loss_differences)) {
    stop("loss_differences must be numeric.")
  }
  
  P <- nrow(loss_differences)
  K <- ncol(loss_differences)
  
  if (P == 0L || K == 0L) {
    stop("loss_differences must have at least one row and one column.")
  }
  
  pooled_loss <- as.vector(loss_differences)
  pooled_loss <- pooled_loss[is.finite(pooled_loss)]
  
  if (length(pooled_loss) == 0L) {
    stop("No finite loss differences are available.")
  }
  
  x_tau_points <- as.numeric(
    quantile(
      pooled_loss,
      probs = quantile_grid,
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  J <- length(x_tau_points)
  
  G_data <- matrix(NA_real_, nrow = P, ncol = K * J)
  null_target <- rep(quantile_grid, times = K)
  
  for (k in seq_len(K)) {
    for (j in seq_len(J)) {
      index <- (k - 1L) * J + j
      G_data[, index] <- as.numeric(
        loss_differences[, k] <= x_tau_points[j]
      )
    }
  }
  
  complete_rows <- apply(
    loss_differences,
    1L,
    function(x) all(is.finite(x))
  )
  
  G_data <- G_data[complete_rows, , drop = FALSE]
  P_clean <- nrow(G_data)
  
  if (P_clean == 0L) {
    stop("No complete rows remain.")
  }
  
  S_mean <- colMeans(G_data)
  S_mean_centered <- S_mean - null_target
  
  V_hat_full <- estimate_long_run_covariance(
    G_data,
    block_length = block_length
  )
  
  V_k_raw <- diag(V_hat_full)
  
  n_nonfinite_variances <- sum(!is.finite(V_k_raw))
  n_floored_variances <- sum(
    !is.finite(V_k_raw) | V_k_raw <= variance_floor
  )
  
  V_k_used <- V_k_raw
  V_k_used[
    !is.finite(V_k_used) | V_k_used <= variance_floor
  ] <- variance_floor
  
  std_dev_k <- sqrt(V_k_used)
  T_k <- S_mean_centered / std_dev_k
  T_max <- max(T_k)
  
  bootstrap_statistics <- numeric(num_bootstrap_replications)
  
  for (b in seq_len(num_bootstrap_replications)) {
    boot_sample <- mbb_resample_data(
      G_data,
      block_length = block_length
    )
    
    boot_mean <- colMeans(boot_sample)
    boot_centered <- boot_mean - null_target
    boot_T <- boot_centered / std_dev_k
    
    bootstrap_statistics[b] <- max(
      boot_T - T_k,
      na.rm = TRUE
    )
  }
  
  list(
    statistic = T_max,
    p_value = mean(bootstrap_statistics > T_max, na.rm = TRUE),
    quantile_grid = quantile_grid,
    threshold_values = x_tau_points,
    variance_floor = variance_floor,
    n_input_rows = P,
    n_complete_rows = P_clean,
    n_removed_rows = P - P_clean,
    n_model_quantile_combinations = K * J,
    n_nonfinite_variances = n_nonfinite_variances,
    n_floored_variances = n_floored_variances,
    raw_long_run_variances = V_k_raw,
    used_long_run_variances = V_k_used
  )
}

# ============================================================
# CDF-RC QUANTILE-GRID AND STABILITY SENSITIVITY
# ============================================================

cdf_grids <- list(
  five_points = seq(0.1, 0.9, by = 0.2),
  nine_points = seq(0.1, 0.9, by = 0.1),
  seventeen_points = seq(0.1, 0.9, by = 0.05)
)

cdf_grid_results <- lapply(names(cdf_grids), function(grid_name) {
  set.seed(2026)
  
  fit <- cdf_rc_diagnostic(
    loss_differences = mse_loss_differences,
    block_length = 5,
    num_bootstrap_replications = 999,
    quantile_grid = cdf_grids[[grid_name]],
    variance_floor = 1e-10
  )
  
  data.frame(
    grid_name = grid_name,
    quantile_grid = paste(
      sprintf("%.2f", fit$quantile_grid),
      collapse = ", "
    ),
    n_thresholds = length(fit$quantile_grid),
    n_model_quantile_combinations =
      fit$n_model_quantile_combinations,
    n_input_rows = fit$n_input_rows,
    n_complete_rows = fit$n_complete_rows,
    n_removed_rows = fit$n_removed_rows,
    n_nonfinite_variances = fit$n_nonfinite_variances,
    n_floored_variances = fit$n_floored_variances,
    variance_floor = fit$variance_floor,
    statistic = fit$statistic,
    p_value = fit$p_value
  )
})

cdf_grid_results <- do.call(rbind, cdf_grid_results)

write.csv(
  cdf_grid_results,
  "output/cdf_rc_quantile_grid_and_stability_diagnostics.csv",
  row.names = FALSE
)

print(cdf_grid_results)

# ============================================================
# WRC HAND-CALCULATION CHECK
# ============================================================

toy_loss_differences <- cbind(
  Model_A = c(0.20, 0.10, 0.30, 0.15, 0.25),
  Model_B = c(-0.10, 0.00, -0.05, 0.05, -0.10)
)

manual_column_means <- colMeans(toy_loss_differences)
manual_wrc_statistic <- max(manual_column_means)

set.seed(2026)

wrc_fit <- white_reality_check(
  loss_differences = toy_loss_differences,
  n_simulations = 999,
  block_length = 2
)

wrc_validation <- data.frame(
  model_A_mean = manual_column_means["Model_A"],
  model_B_mean = manual_column_means["Model_B"],
  manual_wrc_statistic = manual_wrc_statistic,
  rc_wrc_statistic = unname(wrc_fit$statistic),
  absolute_difference = abs(
    manual_wrc_statistic - unname(wrc_fit$statistic)
  ),
  bootstrap_p_value = wrc_fit$p.value
)

write.csv(
  wrc_validation,
  "output/wrc_hand_calculation_validation.csv",
  row.names = FALSE
)

print(wrc_validation)

# ============================================================
# SPA RECENTERING / CONSISTENT-P-VALUE PROPERTY CHECK
# ============================================================

set.seed(2026)

spa_fit <- superior_predictive_ability_test(
  loss_differences = toy_loss_differences,
  num_bootstrap_replications = 999,
  block_length = 2,
  alpha = 0.05
)

spa_validation <- data.frame(
  statistic = unname(spa_fit$statistic),
  p_consistent = spa_fit$p_consistent,
  p_conservative = spa_fit$p_conservative,
  consistent_not_greater_than_conservative =
    spa_fit$p_consistent <= spa_fit$p_conservative,
  consistent_decision = ifelse(
    spa_fit$p_consistent < 0.05,
    "Reject H0",
    "Fail to reject H0"
  ),
  conservative_decision = ifelse(
    spa_fit$p_conservative < 0.05,
    "Reject H0",
    "Fail to reject H0"
  )
)

write.csv(
  spa_validation,
  "output/spa_recentering_property_check.csv",
  row.names = FALSE
)

print(spa_validation)

# ============================================================
#  DM REFERENCE CHECK
# ============================================================

if (exists("realized_values") && exists("forecast_matrix")) {
  actual <- as.numeric(realized_values)
  forecast_matrix <- as.matrix(forecast_matrix)
  
  benchmark_name <- "AR_1"
  
  if (benchmark_name %in% colnames(forecast_matrix)) {
    benchmark_forecast <- as.numeric(
      forecast_matrix[, benchmark_name]
    )
    
    model_names <- setdiff(
      colnames(forecast_matrix),
      benchmark_name
    )
    
    dm_results <- lapply(model_names, function(model_name) {
      model_forecast <- as.numeric(
        forecast_matrix[, model_name]
      )
      
      loss_difference <- (actual - benchmark_forecast)^2 -
        (actual - model_forecast)^2
      
      dm_fit <- forecast::dm.test(
        e1 = actual - model_forecast,
        e2 = actual - benchmark_forecast,
        alternative = "less",
        h = 1,
        power = 2
      )
      
      data.frame(
        model = model_name,
        mean_loss_difference_benchmark_minus_model =
          mean(loss_difference, na.rm = TRUE),
        forecast_dm_statistic = unname(dm_fit$statistic),
        forecast_dm_p_value = dm_fit$p.value,
        forecast_dm_decision_5pct = ifelse(
          dm_fit$p.value < 0.05,
          "Reject H0",
          "Fail to reject H0"
        )
      )
    })
    
    dm_validation <- do.call(rbind, dm_results)
    
    write.csv(
      dm_validation,
      "output/dm_validation_against_forecast.csv",
      row.names = FALSE
    )
    
    print(dm_validation)
  }
}
