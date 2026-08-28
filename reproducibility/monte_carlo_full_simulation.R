#
# Monte Carlo size and power assessment for RCtest forecast-evaluation routines.
#
# Point-forecast procedures use simulated AR(1) loss-difference series.
# KLIC and CRPS/CDF-RC use a comparative Gaussian density-forecast DGP.
# ZP uses a separate tail-probability DGP.
# Kupiec uses a separate variance-miscalibration DGP.
#
# For all joint procedures, a positive loss differential means that a
# competitor has lower loss than the benchmark.
# =============================================================================

library(RCtest)

set.seed(20260819)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_mc <- 500L
n_boot <- 499L
alpha <- 0.05
K_models <- 5L
n_crps_samples <- 500L

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
    Sigma <- matrix(rho, nrow = K, ncol = K)
    diag(Sigma) <- 1
    L <- chol(Sigma)
  } else {
    L <- diag(K)
  }
  
  shocks <- matrix(rnorm(P * K), nrow = P, ncol = K) %*% L
  
  if (heteroskedastic) {
    scale_t <- 1 + 0.7 * sin(seq_len(P) / (P / 4))
    shocks <- shocks * scale_t
  }
  
  loss_diff <- matrix(0, nrow = P, ncol = K)
  loss_diff[1L, ] <- shocks[1L, ]
  
  if (P > 1L) {
    for (t in 2:P) {
      loss_diff[t, ] <- phi * loss_diff[t - 1L, ] + shocks[t, ]
    }
  }
  
  loss_diff <- sweep(loss_diff, 2, alternative_means, FUN = "+")
  colnames(loss_diff) <- paste0("M", seq_len(K))
  loss_diff
}

simulate_comparative_density_dgp <- function(
    P,
    phi,
    correlated,
    heteroskedastic,
    density_advantage,
    competitor_sd_multiplier = 1,
    benchmark_bias = 0,
    benchmark_sd_multiplier = 1,
    K = K_models
) {
  if (length(density_advantage) != K) {
    stop("density_advantage must have length K.")
  }
  
  if (competitor_sd_multiplier <= 0 || benchmark_sd_multiplier <= 0) {
    stop("Standard-deviation multipliers must be positive.")
  }
  
  if (heteroskedastic) {
    scale_t <- 1 + 0.7 * sin(seq_len(P) / (P / 4))
  } else {
    scale_t <- rep(1, P)
  }
  
  innovations <- rnorm(P) * scale_t
  actual <- numeric(P)
  actual[1L] <- innovations[1L]
  
  if (P > 1L) {
    for (t in 2:P) {
      actual[t] <- phi * actual[t - 1L] + innovations[t]
    }
  }
  
  unconditional_sd <- sqrt(1 / (1 - phi^2))
  baseline_sd <- unconditional_sd * scale_t
  
  if (correlated) {
    rho <- 0.5
    Sigma <- matrix(rho, nrow = K + 1L, ncol = K + 1L)
    diag(Sigma) <- 1
    L <- chol(Sigma)
  } else {
    L <- diag(K + 1L)
  }
  
  forecast_shocks <- matrix(
    rnorm(P * (K + 1L)),
    nrow = P,
    ncol = K + 1L
  ) %*% L
  
  forecast_shocks <- forecast_shocks * scale_t
  
  forecast_errors <- matrix(0, nrow = P, ncol = K + 1L)
  forecast_errors[1L, ] <- forecast_shocks[1L, ]
  
  if (P > 1L) {
    for (t in 2:P) {
      forecast_errors[t, ] <- phi * forecast_errors[t - 1L, ] +
        forecast_shocks[t, ]
    }
  }
  
  benchmark_mean <- actual + forecast_errors[, K + 1L] + benchmark_bias
  
  competitor_means <- forecast_errors[, seq_len(K), drop = FALSE]
  competitor_means <- sweep(competitor_means, 1, actual, FUN = "+")
  competitor_means <- sweep(competitor_means, 2, density_advantage, FUN = "-")
  colnames(competitor_means) <- paste0("M", seq_len(K))
  
  benchmark_sd <- baseline_sd * benchmark_sd_multiplier
  
  competitor_sds <- matrix(
    rep(baseline_sd * competitor_sd_multiplier, K),
    nrow = P,
    ncol = K,
    dimnames = list(NULL, colnames(competitor_means))
  )
  
  list(
    actual = actual,
    competitor_means = competitor_means,
    competitor_sds = competitor_sds,
    benchmark_mean = benchmark_mean,
    benchmark_sd = benchmark_sd
  )
}

compute_klic_comparative <- function(
    actual,
    competitor_means,
    competitor_sds,
    benchmark_mean,
    benchmark_sd
) {
  P <- length(actual)
  K <- ncol(competitor_means)
  
  if (nrow(competitor_means) != P ||
      nrow(competitor_sds) != P ||
      ncol(competitor_sds) != K ||
      length(benchmark_mean) != P ||
      length(benchmark_sd) != P) {
    stop("Incompatible dimensions in comparative KLIC inputs.")
  }
  
  benchmark_nls <- -dnorm(
    actual,
    mean = benchmark_mean,
    sd = benchmark_sd,
    log = TRUE
  )
  
  competitor_nls <- matrix(
    NA_real_,
    nrow = P,
    ncol = K,
    dimnames = list(NULL, colnames(competitor_means))
  )
  
  for (k in seq_len(K)) {
    competitor_nls[, k] <- -dnorm(
      actual,
      mean = competitor_means[, k],
      sd = competitor_sds[, k],
      log = TRUE
    )
  }
  
  sweep(
    competitor_nls,
    1,
    benchmark_nls,
    FUN = function(model_loss, benchmark_loss) benchmark_loss - model_loss
  )
}

compute_zp_comparative <- function(
    actual,
    threshold,
    competitor_means,
    competitor_sds,
    benchmark_mean,
    benchmark_sd
) {
  P <- length(actual)
  K <- ncol(competitor_means)
  
  if (length(threshold) == 1L) {
    threshold <- rep(threshold, P)
  }
  
  if (length(threshold) != P ||
      nrow(competitor_means) != P ||
      nrow(competitor_sds) != P ||
      ncol(competitor_sds) != K ||
      length(benchmark_mean) != P ||
      length(benchmark_sd) != P) {
    stop("Incompatible dimensions in comparative ZP inputs.")
  }
  
  event <- as.numeric(actual <= threshold)
  
  benchmark_probability <- pnorm(
    threshold,
    mean = benchmark_mean,
    sd = benchmark_sd
  )
  
  benchmark_loss <- (event - benchmark_probability)^2
  
  competitor_loss <- matrix(
    NA_real_,
    nrow = P,
    ncol = K,
    dimnames = list(NULL, colnames(competitor_means))
  )
  
  for (k in seq_len(K)) {
    competitor_probability <- pnorm(
      threshold,
      mean = competitor_means[, k],
      sd = competitor_sds[, k]
    )
    
    competitor_loss[, k] <- (event - competitor_probability)^2
  }
  
  sweep(
    competitor_loss,
    1,
    benchmark_loss,
    FUN = function(model_loss, benchmark_loss) benchmark_loss - model_loss
  )
}

compute_crps_loss_differences <- function(
    competitor_means,
    competitor_sds,
    benchmark_mean,
    benchmark_sd,
    actual,
    n_samples
) {
  P <- length(actual)
  K <- ncol(competitor_means)
  
  competitor_crps <- matrix(
    NA_real_,
    nrow = P,
    ncol = K,
    dimnames = list(NULL, colnames(competitor_means))
  )
  
  benchmark_crps <- numeric(P)
  
  for (t in seq_len(P)) {
    benchmark_samples <- rnorm(
      n_samples,
      mean = benchmark_mean[t],
      sd = benchmark_sd[t]
    )
    
    benchmark_crps[t] <- compute_crps(benchmark_samples, actual[t])
    
    for (k in seq_len(K)) {
      model_samples <- rnorm(
        n_samples,
        mean = competitor_means[t, k],
        sd = competitor_sds[t, k]
      )
      
      competitor_crps[t, k] <- compute_crps(model_samples, actual[t])
    }
  }
  
  sweep(
    competitor_crps,
    1,
    benchmark_crps,
    FUN = function(model_score, benchmark_score) benchmark_score - model_score
  )
}

simulate_zp_tail_dgp <- function(
    P,
    phi,
    correlated,
    heteroskedastic,
    alternative,
    K = K_models
) {
  if (heteroskedastic) {
    scale_t <- 1 + 0.7 * sin(seq_len(P) / (P / 4))
  } else {
    scale_t <- rep(1, P)
  }
  
  innovations <- rnorm(P) * scale_t
  actual <- numeric(P)
  actual[1L] <- innovations[1L]
  
  if (P > 1L) {
    for (t in 2:P) {
      actual[t] <- phi * actual[t - 1L] + innovations[t]
    }
  }
  
  unconditional_sd <- sqrt(1 / (1 - phi^2))
  threshold <- qnorm(alpha) * unconditional_sd
  
  if (correlated) {
    rho <- 0.5
    Sigma <- matrix(rho, nrow = K + 1L, ncol = K + 1L)
    diag(Sigma) <- 1
    L <- chol(Sigma)
  } else {
    L <- diag(K + 1L)
  }
  
  forecast_noise <- matrix(
    rnorm(P * (K + 1L)),
    nrow = P,
    ncol = K + 1L
  ) %*% L
  
  forecast_noise <- 0.10 * forecast_noise
  
  competitor_means <- forecast_noise[, seq_len(K), drop = FALSE]
  colnames(competitor_means) <- paste0("M", seq_len(K))
  benchmark_mean <- forecast_noise[, K + 1L]
  
  competitor_sds <- matrix(
    unconditional_sd,
    nrow = P,
    ncol = K,
    dimnames = list(NULL, colnames(competitor_means))
  )
  
  benchmark_sd <- rep(unconditional_sd, P)
  
  if (alternative) {
    benchmark_sd <- rep(0.45 * unconditional_sd, P)
  }
  
  list(
    actual = actual,
    threshold = threshold,
    competitor_means = competitor_means,
    competitor_sds = competitor_sds,
    benchmark_mean = benchmark_mean,
    benchmark_sd = benchmark_sd
  )
}

run_full_battery <- function(
    loss_diff,
    P,
    phi,
    correlated,
    heteroskedastic,
    block_length,
    n_boot,
    alternative_means,
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
    error = function(e) NULL
  )
  
  if (!is.null(wrc)) results$WRC <- wrc$p.value
  
  spa <- tryCatch(
    superior_predictive_ability_test(
      loss_differences = loss_diff,
      block_length = block_length,
      num_bootstrap_replications = n_boot,
      alpha = alpha
    ),
    error = function(e) NULL
  )
  
  if (!is.null(spa)) results$SPA <- spa$p_consistent
  
  conditioning_variable <- abs(rnorm(P))
  
  cpa <- tryCatch(
    white_reality_check_conditional(
      loss_differences = loss_diff,
      weighting_vector = conditioning_variable,
      block_length = block_length,
      num_bootstrap_replications = n_boot,
      alpha = alpha
    ),
    error = function(e) NULL
  )
  
  if (!is.null(cpa)) results$CPA <- cpa$p.value
  
  dm <- tryCatch(
    compute_per_model_statistics(
      loss_differences = loss_diff,
      model_names = colnames(loss_diff),
      n_boot = n_boot,
      block_length = block_length,
      alpha = alpha,
      H1 = "same"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(dm)) {
    results$DM_mean_rejection_rate <- mean(dm$Significant, na.rm = TRUE)
  }
  
  is_alternative <- any(alternative_means != 0)
  
  benchmark_bias <- if (is_alternative) 0.40 else 0
  benchmark_sd_multiplier <- if (is_alternative) 1.75 else 1
  
  density_dgp <- simulate_comparative_density_dgp(
    P = P,
    phi = phi,
    correlated = correlated,
    heteroskedastic = heteroskedastic,
    density_advantage = alternative_means,
    competitor_sd_multiplier = 1,
    benchmark_bias = benchmark_bias,
    benchmark_sd_multiplier = benchmark_sd_multiplier,
    K = K
  )
  
  klic_diff <- tryCatch(
    compute_klic_comparative(
      actual = density_dgp$actual,
      competitor_means = density_dgp$competitor_means,
      competitor_sds = density_dgp$competitor_sds,
      benchmark_mean = density_dgp$benchmark_mean,
      benchmark_sd = density_dgp$benchmark_sd
    ),
    error = function(e) NULL
  )
  
  if (!is.null(klic_diff)) {
    klic <- tryCatch(
      kullback_leibler_test(
        log_likelihood_differences = klic_diff,
        block_length = block_length,
        num_bootstrap_replications = n_boot,
        alpha = alpha
      ),
      error = function(e) NULL
    )
    
    if (!is.null(klic)) results$KLIC <- klic$p.value
  }
  
  zp_dgp <- simulate_zp_tail_dgp(
    P = P,
    phi = phi,
    correlated = correlated,
    heteroskedastic = heteroskedastic,
    alternative = is_alternative,
    K = K
  )
  
  zp_diff <- tryCatch(
    compute_zp_comparative(
      actual = zp_dgp$actual,
      threshold = zp_dgp$threshold,
      competitor_means = zp_dgp$competitor_means,
      competitor_sds = zp_dgp$competitor_sds,
      benchmark_mean = zp_dgp$benchmark_mean,
      benchmark_sd = zp_dgp$benchmark_sd
    ),
    error = function(e) NULL
  )
  
  if (!is.null(zp_diff)) {
    zp <- tryCatch(
      reality_check_zp_test(
        zp_loss_differences = zp_diff,
        block_length = block_length,
        num_bootstrap_replications = n_boot,
        alpha = alpha
      ),
      error = function(e) NULL
    )
    
    if (!is.null(zp)) results$ZP <- zp$p_consistent
  }
  
  crps_diff <- tryCatch(
    compute_crps_loss_differences(
      competitor_means = density_dgp$competitor_means,
      competitor_sds = density_dgp$competitor_sds,
      benchmark_mean = density_dgp$benchmark_mean,
      benchmark_sd = density_dgp$benchmark_sd,
      actual = density_dgp$actual,
      n_samples = n_crps_samples
    ),
    error = function(e) NULL
  )
  
  if (!is.null(crps_diff)) {
    cdf_rc <- tryCatch(
      white_reality_check_cdf_approx(
        loss_differences = crps_diff,
        block_length = block_length,
        num_bootstrap_replications = n_boot,
        alpha = alpha
      ),
      error = function(e) NULL
    )
    
    if (!is.null(cdf_rc)) results$CDF_RC <- cdf_rc$p.value
  }
  
  kupiec_dgp <- simulate_comparative_density_dgp(
    P = P,
    phi = phi,
    correlated = correlated,
    heteroskedastic = heteroskedastic,
    density_advantage = rep(0, K),
    competitor_sd_multiplier = if (is_alternative) 0.60 else 1,
    benchmark_bias = 0,
    benchmark_sd_multiplier = 1,
    K = K
  )
  
  kupiec_matrix <- cbind(
    kupiec_dgp$competitor_means,
    Benchmark = kupiec_dgp$benchmark_mean
  )
  
  kupiec_results <- tryCatch(
    compute_kupiec(
      forecast_matrix = kupiec_matrix,
      forecast_sd_models = kupiec_dgp$competitor_sds,
      realized = kupiec_dgp$actual,
      benchmark_col = ncol(kupiec_matrix),
      alpha = alpha
    ),
    error = function(e) NULL
  )
  
  if (!is.null(kupiec_results) && length(kupiec_results) > 0L) {
    kupiec_p_values <- vapply(
      kupiec_results,
      function(x) x$p.value,
      numeric(1)
    )
    
    results$Kupiec_mean_rejection_rate <- mean(
      kupiec_p_values <= alpha,
      na.rm = TRUE
    )
  }
  
  return(results)
}

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
    configs[[length(configs) + 1L]] <- modifyList(base_config, list(P = P_value))
  }
  
  for (phi_value in c(0.0, 0.6)) {
    configs[[length(configs) + 1L]] <- modifyList(base_config, list(phi = phi_value))
  }
  
  configs[[length(configs) + 1L]] <- modifyList(base_config, list(correlated = FALSE))
  configs[[length(configs) + 1L]] <- modifyList(base_config, list(heteroskedastic = TRUE))
  
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

results_file <- file.path(output_dir, "monte_carlo_results.csv")

if (file.exists(results_file)) file.remove(results_file)

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

joint_test_names <- c("WRC", "SPA", "CPA", "KLIC", "ZP", "CDF_RC")
per_model_test_names <- c("DM_mean_rejection_rate", "Kupiec_mean_rejection_rate")

for (config_index in seq_along(configs)) {
  config <- configs[[config_index]]
  
  cat(sprintf(
    "\n[%d/%d] Scenario=%s, P=%d, phi=%.1f, correlated=%s, heteroskedastic=%s, block_length=%d, means=(%s), Kupiec-SD multiplier=%.2f\n",
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
  ))
  
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
      alternative_means = config$alternative_means,
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
    valid <- is.finite(values)
    n_valid_draws[test_name] <- sum(valid)
    
    if (n_valid_draws[test_name] > 0L) {
      rejection_rate[test_name] <- mean(values[valid] <= alpha)
    }
  }
  
  for (test_name in per_model_test_names) {
    values <- battery_values[, test_name]
    valid <- is.finite(values)
    n_valid_draws[test_name] <- sum(valid)
    
    if (n_valid_draws[test_name] > 0L) {
      rejection_rate[test_name] <- mean(values[valid])
    }
  }
  
  monte_carlo_se <- sqrt(rejection_rate * (1 - rejection_rate) / n_valid_draws)
  
  config_results <- data.frame(
    Scenario = config$Scenario,
    P = config$P,
    phi = config$phi,
    correlated = config$correlated,
    heteroskedastic = config$heteroskedastic,
    block_length = config$block_length,
    Alternative_Means = paste(config$alternative_means, collapse = ";"),
    Kupiec_SD_Multiplier = config$kupiec_miscalibration,
    CRPS_Samples = n_crps_samples,
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
  
  cat(sprintf(
    "  Done: %s\n",
    paste(
      sprintf("%s=%.3f(SE=%.3f)", test_names, rejection_rate, monte_carlo_se),
      collapse = ", "
    )
  ))
}

cat("\n==== Complete. Results saved to ====\n")
cat(normalizePath(results_file), "\n")
print(read.csv(results_file), row.names = FALSE)