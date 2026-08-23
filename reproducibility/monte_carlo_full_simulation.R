library(RCtest)

set.seed(20260819)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_mc <- 200L
n_boot <- 999L
alpha <- 0.05
K_models <- 5L

# -----------------------------------------------------------------------------
# Data-generating processes
# -----------------------------------------------------------------------------

simulate_loss_differentials <- function(
    P,
    phi,
    correlated,
    heteroskedastic,
    K = K_models
) {
  if (correlated) {
    rho <- 0.5
    sigma_matrix <- matrix(rho, nrow = K, ncol = K)
    diag(sigma_matrix) <- 1
    cholesky_factor <- chol(sigma_matrix)
  } else {
    cholesky_factor <- diag(K)
  }
  
  shocks <- matrix(rnorm(P * K), nrow = P, ncol = K) %*% cholesky_factor
  
  if (heteroskedastic) {
    scale <- 1 + 0.7 * sin(seq_len(P) / (P / 4))
    shocks <- shocks * scale
  }
  
  loss_differences <- matrix(0, nrow = P, ncol = K)
  
  for (model_index in seq_len(K)) {
    innovation <- shocks[, model_index]
    series <- numeric(P)
    series[1] <- innovation[1]
    
    for (time_index in 2:P) {
      series[time_index] <- phi * series[time_index - 1] + innovation[time_index]
    }
    
    loss_differences[, model_index] <- series
  }
  
  colnames(loss_differences) <- paste0("M", seq_len(K))
  
  loss_differences
}

simulate_forecast_levels_for_kupiec <- function(
    P,
    phi,
    correlated,
    heteroskedastic,
    K = K_models
) {
  realized <- cumsum(rnorm(P, sd = 1))
  
  forecast_errors <- simulate_loss_differentials(
    P = P,
    phi = phi,
    correlated = correlated,
    heteroskedastic = heteroskedastic,
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
  
  # A positive time-varying predictive SD for each non-benchmark model.
  sd_base <- sd(diff(realized))
  forecast_sd_models <- matrix(
    sd_base * (1 + 0.1 * abs(rnorm(P * K))),
    nrow = P,
    ncol = K,
    dimnames = list(NULL, paste0("M", seq_len(K)))
  )
  
  list(
    forecast_matrix = forecast_matrix,
    forecast_sd_models = forecast_sd_models,
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
    n_boot
) {
  K <- ncol(loss_diff)
  
  results <- list(
    WRC = NA_real_,
    SPA = NA_real_,
    CPA = NA_real_,
    KLIC = NA_real_,
    ZP = NA_real_,
    CDF_RC = NA_real_,
    DM_any_reject = NA_real_,
    Kupiec_any_reject = NA_real_
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
  
  # Simulated NLS-difference input for the KLIC Reality Check.
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
    results$DM_any_reject <- as.numeric(
      any(dm$Significant, na.rm = TRUE)
    )
  }
  
  kupiec_simulation <- simulate_forecast_levels_for_kupiec(
    P = P,
    phi = phi,
    correlated = correlated,
    heteroskedastic = heteroskedastic,
    K = K
  )
  
  kupiec_results <- tryCatch(
    compute_kupiec(
      forecast_matrix = kupiec_simulation$forecast_matrix,
      forecast_sd_models = kupiec_simulation$forecast_sd_models,
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
    
    results$Kupiec_any_reject <- as.numeric(
      any(kupiec_p_values <= alpha, na.rm = TRUE)
    )
  }
  
  results
}

# -----------------------------------------------------------------------------
# Simulation configurations
# -----------------------------------------------------------------------------

baseline <- list(
  P = 150L,
  phi = 0.3,
  correlated = TRUE,
  heteroskedastic = FALSE,
  block_length = 5L
)

configs <- list()

for (P_value in c(100L, 150L, 250L)) {
  configs[[length(configs) + 1L]] <- modifyList(
    baseline,
    list(P = P_value)
  )
}

for (phi_value in c(0.0, 0.6)) {
  configs[[length(configs) + 1L]] <- modifyList(
    baseline,
    list(phi = phi_value)
  )
}

configs[[length(configs) + 1L]] <- modifyList(
  baseline,
  list(correlated = FALSE)
)

configs[[length(configs) + 1L]] <- modifyList(
  baseline,
  list(heteroskedastic = TRUE)
)

for (block_length_value in c(3L, 8L)) {
  configs[[length(configs) + 1L]] <- modifyList(
    baseline,
    list(block_length = block_length_value)
  )
}

configs[[length(configs) + 1L]] <- baseline

config_keys <- vapply(
  configs,
  function(config) {
    paste(
      config$P,
      config$phi,
      config$correlated,
      config$heteroskedastic,
      config$block_length,
      sep = "_"
    )
  },
  character(1)
)

configs <- configs[!duplicated(config_keys)]

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
  "DM_any_reject",
  "Kupiec_any_reject"
)

for (config_index in seq_along(configs)) {
  config <- configs[[config_index]]
  
  cat(
    sprintf(
      "\n[%d/%d] P=%d, phi=%.1f, correlated=%s, heteroskedastic=%s, block_length=%d\n",
      config_index,
      length(configs),
      config$P,
      config$phi,
      config$correlated,
      config$heteroskedastic,
      config$block_length
    )
  )
  
  p_values <- matrix(
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
      heteroskedastic = config$heteroskedastic
    )
    
    battery <- run_full_battery(
      loss_diff = loss_diff,
      P = config$P,
      phi = config$phi,
      correlated = config$correlated,
      heteroskedastic = config$heteroskedastic,
      block_length = config$block_length,
      n_boot = n_boot
    )
    
    for (test_name in test_names) {
      p_values[simulation_index, test_name] <- battery[[test_name]]
    }
    
    if (simulation_index %% 50L == 0L) {
      cat(sprintf("  ...%d/%d\n", simulation_index, n_mc))
    }
  }
  
  rejection_rate <- colMeans(
    p_values <= alpha,
    na.rm = TRUE
  )
  
  n_valid_draws <- colSums(!is.na(p_values))
  
  monte_carlo_se <- sqrt(
    rejection_rate * (1 - rejection_rate) / n_valid_draws
  )
  
  config_results <- data.frame(
    P = config$P,
    phi = config$phi,
    correlated = config$correlated,
    heteroskedastic = config$heteroskedastic,
    block_length = config$block_length,
    Test = test_names,
    Rejection_Rate = round(rejection_rate, 4),
    MC_SE = round(monte_carlo_se, 4),
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