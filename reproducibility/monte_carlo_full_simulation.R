# =============================================================================
# monte_carlo_full_simulation.R
#
# Monte Carlo size and power assessment for RCtest routines.
#
# Joint tests:
#   WRC, SPA, CPA, KLIC, ZP, CDF-RC
# report the fraction of Monte Carlo replications rejecting at alpha.
#
# Per-model tests:
#   DM_mean_rejection_rate and Kupiec_mean_rejection_rate report the average
#   fraction of individual model tests rejecting across MC replications.
# They are not family-wise "any rejection" probabilities.
#
# The loss-differential DGP has null and alternative cases. The Kupiec DGP
# separately has calibrated and deliberately miscalibrated VaR scenarios.
# =============================================================================

library(RCtest)

set.seed(20260819)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_mc <- 500L
n_boot <- 499L
alpha <- 0.05
K_models <- 5L

# -----------------------------------------------------------------------------
# Generic AR(1) loss-differential DGP
# -----------------------------------------------------------------------------

simulate_loss_differentials <- function(
    P,
    phi,
    correlated,
    heteroskedastic,
    alternative_means = rep(0, K_models),
    K = K_models
) {
  if (length(alternative_means) != K) {
    stop("alternative_means must have length K.")
  }
  
  if (correlated) {
    rho <- 0.5
    sigma_matrix <- matrix(rho, nrow = K, ncol = K)
    diag(sigma_matrix) <- 1
    cholesky_factor <- chol(sigma_matrix)
  } else {
    cholesky_factor <- diag(K)
  }
  
  shocks <- matrix(
    rnorm(P * K),
    nrow = P,
    ncol = K
  ) %*% cholesky_factor
  
  if (heteroskedastic) {
    innovation_scale <- 1 + 0.7 * sin(seq_len(P) / (P / 4))
    shocks <- shocks * innovation_scale
  }
  
  loss_differences <- matrix(0, nrow = P, ncol = K)
  
  for (model_index in seq_len(K)) {
    innovation <- shocks[, model_index]
    series <- numeric(P)
    series[1] <- innovation[1]
    
    for (time_index in 2:P) {
      series[time_index] <- (
        phi * series[time_index - 1] +
          innovation[time_index]
      )
    }
    
    loss_differences[, model_index] <- (
      series + alternative_means[model_index]
    )
  }
  
  colnames(loss_differences) <- paste0("M", seq_len(K))
  
  loss_differences
}

# -----------------------------------------------------------------------------
# VaR / Kupiec DGP
#
# Under kupiec_miscalibration = 1 and heteroskedastic = FALSE, the supplied
# forecast SD equals the AR(1) error process's unconditional SD. Therefore,
# each model has correctly calibrated Gaussian one-step VaR under the null.
#
# Under kupiec_miscalibration < 1, forecast SD is deliberately too small:
# VaR bands are too narrow and exceedances occur too often. This is a
# controlled alternative-power configuration for the Kupiec UC test.
# -----------------------------------------------------------------------------

simulate_forecast_levels_for_kupiec <- function(
    P,
    phi,
    correlated,
    heteroskedastic,
    kupiec_miscalibration = 1,
    K = K_models
) {
  if (kupiec_miscalibration <= 0) {
    stop("kupiec_miscalibration must be positive.")
  }
  
  realized <- cumsum(rnorm(P, sd = 1))
  
  forecast_errors <- simulate_loss_differentials(
    P = P,
    phi = phi,
    correlated = correlated,
    heteroskedastic = heteroskedastic,
    alternative_means = rep(0, K + 1L),
    K = K + 1L
  )
  
  forecast_matrix <- matrix(
    realized,
    nrow = P,
    ncol = K + 1L
  ) + forecast_errors
  
  colnames(forecast_matrix) <- c(
    paste0("M", seq_len(K)),
    "Benchmark"
  )
  
  benchmark_col <- K + 1L
  
  # For homoskedastic innovations with variance 1, the stationary AR(1)
  # forecast-error SD is sqrt(1 / (1 - phi^2)).
  unconditional_error_sd <- sqrt(1 / (1 - phi^2))
  
  if (heteroskedastic) {
    innovation_scale <- 1 + 0.7 * sin(seq_len(P) / (P / 4))
    
    sd_path <- unconditional_error_sd * innovation_scale
  } else {
    sd_path <- rep(unconditional_error_sd, P)
  }
  
  forecast_sd_models <- matrix(
    rep(sd_path * kupiec_miscalibration, K),
    nrow = P,
    ncol = K,
    dimnames = list(NULL, paste0("M", seq_len(K)))
  )
  
  list(
    forecast_matrix = forecast_matrix,
    forecast_sd_models = forecast_sd_models,
    realized = realized,
    benchmark_col = benchmark_col
  )
}

# -----------------------------------------------------------------------------
# Test battery
# -----------------------------------------------------------------------------

run_full_battery <- function(
    loss_diff,
    P,
    phi,
    correlated,
    heteroskedastic,
    block_length,
    n_boot,
    kupiec_miscalibration
) {
  K <- ncol(loss_diff)
  
  results <- list(
    WRC = NA_real_,
    SPA = NA_real_,
    CPA = NA_real_,
    KLIC = NA_real_,
    ZP = NA_real_,
    CDF_RC = NA_real_,
    DM_mean_rejection_rate = NA_real_,
    Kupiec_mean_rejection_rate = NA_real_
  )
  
  wrc <- tryCatch(
    white_reality_check(
      loss_differences = loss_diff,
      n_simulations = n_boot,
      block_length = block_length,
      alpha = alpha
    ),
    error = function(error) NULL
  )
  
  if (!is.null(wrc)) {
    results$WRC <- wrc$p.value
  }
  
  spa <- tryCatch(
    superior_predictive_ability_test(
      loss_differences = loss_diff,
      block_length = block_length,
      num_bootstrap_replications = n_boot,
      alpha = alpha
    ),
    error = function(error) NULL
  )
  
  if (!is.null(spa)) {
    results$SPA <- spa$p_consistent
  }
  
  conditioning_variable <- abs(rnorm(P))
  
  cpa <- tryCatch(
    white_reality_check_conditional(
      loss_differences = loss_diff,
      weighting_vector = conditioning_variable,
      block_length = block_length,
      num_bootstrap_replications = n_boot,
      alpha = alpha
    ),
    error = function(error) NULL
  )
  
  if (!is.null(cpa)) {
    results$CPA <- cpa$p.value
  }
  
  klic_input <- loss_diff / 2 +
    matrix(rnorm(P * K, sd = 0.1), nrow = P, ncol = K)
  
  klic <- tryCatch(
    kullback_leibler_test(
      log_likelihood_differences = klic_input,
      block_length = block_length,
      num_bootstrap_replications = n_boot,
      alpha = alpha
    ),
    error = function(error) NULL
  )
  
  if (!is.null(klic)) {
    results$KLIC <- klic$p.value
  }
  
  zp <- tryCatch(
    reality_check_zp_test(
      zp_loss_differences = loss_diff,
      block_length = block_length,
      num_bootstrap_replications = n_boot,
      alpha = alpha
    ),
    error = function(error) NULL
  )
  
  if (!is.null(zp)) {
    results$ZP <- zp$p_consistent
  }
  
  cdf_rc <- tryCatch(
    white_reality_check_cdf_approx(
      loss_differences = loss_diff,
      block_length = block_length,
      num_bootstrap_replications = n_boot,
      alpha = alpha
    ),
    error = function(error) NULL
  )
  
  if (!is.null(cdf_rc)) {
    results$CDF_RC <- cdf_rc$p.value
  }
  
  dm <- tryCatch(
    compute_per_model_statistics(
      loss_differences = loss_diff,
      model_names = colnames(loss_diff),
      n_boot = n_boot,
      block_length = block_length,
      alpha = alpha,
      H1 = "same"
    ),
    error = function(error) NULL
  )
  
  if (!is.null(dm)) {
    results$DM_mean_rejection_rate <- mean(
      dm$Significant,
      na.rm = TRUE
    )
  }
  
  kupiec_simulation <- simulate_forecast_levels_for_kupiec(
    P = P,
    phi = phi,
    correlated = correlated,
    heteroskedastic = heteroskedastic,
    kupiec_miscalibration = kupiec_miscalibration,
    K = K
  )
  
  kupiec_results <- tryCatch(
    compute_kupiec(
      forecast_matrix = kupiec_simulation$forecast_matrix,
      forecast_sd_models = kupiec_simulation$forecast_sd_models,
      realized = kupiec_simulation$realized,
      benchmark_col = kupiec_simulation$benchmark_col,
      alpha = alpha
    ),
    error = function(error) NULL
  )
  
  if (!is.null(kupiec_results) && length(kupiec_results) > 0L) {
    kupiec_p_values <- vapply(
      kupiec_results,
      function(result) result$p.value,
      numeric(1)
    )
    
    results$Kupiec_mean_rejection_rate <- mean(
      kupiec_p_values <= alpha,
      na.rm = TRUE
    )
  }
  
  results
}

# -----------------------------------------------------------------------------
# Simulation scenarios
# -----------------------------------------------------------------------------

null_baseline <- list(
  Scenario = "Null",
  P = 150L,
  phi = 0.3,
  correlated = TRUE,
  heteroskedastic = FALSE,
  block_length = 5L,
  alternative_means = rep(0, K_models),
  kupiec_miscalibration = 1
)

alternative_baseline <- modifyList(
  null_baseline,
  list(
    Scenario = "Alternative",
    alternative_means = c(0.30, 0.30, 0.15, 0, 0),
    kupiec_miscalibration = 0.60
  )
)

make_config_grid <- function(base_config) {
  configs <- list()
  
  for (P_value in c(100L, 150L, 250L)) {
    configs[[length(configs) + 1L]] <- modifyList(
      base_config,
      list(P = P_value)
    )
  }
  
  for (phi_value in c(0.0, 0.6)) {
    configs[[length(configs) + 1L]] <- modifyList(
      base_config,
      list(phi = phi_value)
    )
  }
  
  configs[[length(configs) + 1L]] <- modifyList(
    base_config,
    list(correlated = FALSE)
  )
  
  configs[[length(configs) + 1L]] <- modifyList(
    base_config,
    list(heteroskedastic = TRUE)
  )
  
  for (block_length_value in c(3L, 8L)) {
    configs[[length(configs) + 1L]] <- modifyList(
      base_config,
      list(block_length = block_length_value)
    )
  }
  
  configs[[length(configs) + 1L]] <- base_config
  
  config_keys <- vapply(
    configs,
    function(config) {
      paste(
        config$Scenario,
        config$P,
        config$phi,
        config$correlated,
        config$heteroskedastic,
        config$block_length,
        paste(config$alternative_means, collapse = "_"),
        config$kupiec_miscalibration,
        sep = "_"
      )
    },
    character(1)
  )
  
  configs[!duplicated(config_keys)]
}

configs <- c(
  make_config_grid(null_baseline),
  make_config_grid(alternative_baseline)
)

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

results_file <- file.path(
  output_dir,
  "monte_carlo_results.csv"
)

if (file.exists(results_file)) {
  file.remove(results_file)
}

test_names <- c(
  "WRC",
  "SPA",
  "CPA",
  "KLIC",
  "ZP",
  "CDF_RC",
  "DM_mean_rejection_rate",
  "Kupiec_mean_rejection_rate"
)

joint_test_names <- c(
  "WRC",
  "SPA",
  "CPA",
  "KLIC",
  "ZP",
  "CDF_RC"
)

per_model_test_names <- c(
  "DM_mean_rejection_rate",
  "Kupiec_mean_rejection_rate"
)

for (config_index in seq_along(configs)) {
  config <- configs[[config_index]]
  
  cat(
    sprintf(
      paste0(
        "\n[%d/%d] Scenario=%s, P=%d, phi=%.1f, correlated=%s, ",
        "heteroskedastic=%s, block_length=%d, means=(%s), VaR-SD multiplier=%.2f\n"
      ),
      config_index,
      length(configs),
      config$Scenario,
      config$P,
      config$phi,
      config$correlated,
      config$heteroskedastic,
      config$block_length,
      paste(config$alternative_means, collapse = ", "),
      config$kupiec_miscalibration
    )
  )
  
  battery_values <- matrix(
    NA_real_,
    nrow = n_mc,
    ncol = length(test_names),
    dimnames = list(NULL, test_names)
  )
  
  for (simulation_index in seq_len(n_mc)) {
    loss_diff <- simulate_loss_differentials(
      P = config$P,
      phi = config$phi,
      correlated = config$correlated,
      heteroskedastic = config$heteroskedastic,
      alternative_means = config$alternative_means,
      K = K_models
    )
    
    battery <- run_full_battery(
      loss_diff = loss_diff,
      P = config$P,
      phi = config$phi,
      correlated = config$correlated,
      heteroskedastic = config$heteroskedastic,
      block_length = config$block_length,
      n_boot = n_boot,
      kupiec_miscalibration = config$kupiec_miscalibration
    )
    
    for (test_name in test_names) {
      battery_values[simulation_index, test_name] <- battery[[test_name]]
    }
    
    if (simulation_index %% 50L == 0L) {
      cat(sprintf("  ...%d/%d\n", simulation_index, n_mc))
    }
  }
  
  rejection_rate <- rep(NA_real_, length(test_names))
  n_valid_draws <- integer(length(test_names))
  names(rejection_rate) <- test_names
  names(n_valid_draws) <- test_names
  
  for (test_name in joint_test_names) {
    values <- battery_values[, test_name]
    valid <- !is.na(values)
    n_valid_draws[test_name] <- sum(valid)
    
    if (n_valid_draws[test_name] > 0L) {
      rejection_rate[test_name] <- mean(values[valid] <= alpha)
    }
  }
  
  for (test_name in per_model_test_names) {
    values <- battery_values[, test_name]
    valid <- !is.na(values)
    n_valid_draws[test_name] <- sum(valid)
    
    if (n_valid_draws[test_name] > 0L) {
      rejection_rate[test_name] <- mean(values[valid])
    }
  }
  
  monte_carlo_se <- sqrt(
    rejection_rate * (1 - rejection_rate) / n_valid_draws
  )
  
  config_results <- data.frame(
    Scenario = config$Scenario,
    P = config$P,
    phi = config$phi,
    correlated = config$correlated,
    heteroskedastic = config$heteroskedastic,
    block_length = config$block_length,
    Alternative_Means = paste(config$alternative_means, collapse = ";"),
    Kupiec_SD_Multiplier = config$kupiec_miscalibration,
    Test = test_names,
    Rejection_Rate = rejection_rate,
    MC_SE = monte_carlo_se,
    N_Valid_Draws = n_valid_draws,
    stringsAsFactors = FALSE
  )
  
  write.table(
    config_results,
    results_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(results_file),
    append = file.exists(results_file)
  )
  
  cat(
    sprintf(
      "  Done: %s\n",
      paste(
        sprintf(
          "%s=%.3f(SE=%.3f)",
          test_names,
          rejection_rate,
          monte_carlo_se
        ),
        collapse = ", "
      )
    )
  )
}

cat("\n==== Complete. Results saved to ====\n")
cat(normalizePath(results_file), "\n")

print(read.csv(results_file), row.names = FALSE)