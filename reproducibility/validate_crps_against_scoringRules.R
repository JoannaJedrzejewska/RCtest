# =============================================================================
# Direct numerical validation of RCtest::compute_crps() against the analytic
# Gaussian CRPS from scoringRules::crps_norm().
#
# This script separates two distinct issues:
#
#   1. Score arithmetic: compute_crps() receives a genuine predictive sample
#      from a known Gaussian distribution and is compared with its analytic
#      Gaussian CRPS. This validates the sample/empirical CRPS calculation.
#
#   2. Predictive-distribution interpretation: a cross-sectional vector of
#      point forecasts from DIFFERENT models is not automatically a
#      model-specific predictive distribution. The second output therefore
#      documents the properties of an ensemble-disagreement proxy; it does
#      NOT claim to validate it as conventional, model-specific CRPS.
#
# Output:
#   reproducibility/output/crps_validation_against_scoringRules.csv
#   reproducibility/output/crps_validation_against_scoringRules_summary.csv
#   reproducibility/output/crps_ensemble_disagreement_proxy_demo.csv
# =============================================================================

library(RCtest)

if (!requireNamespace("scoringRules", quietly = TRUE)) {
  stop("Package 'scoringRules' is required for this validation script.")
}

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260824)

closed_form_gaussian_crps <- function(y, mu, sigma) {
  scoringRules::crps_norm(y, mean = mu, sd = sigma)
}

validation_configs <- data.frame(
  Config = seq_len(4L),
  mu = c(0, 0, 5, -2),
  sigma = c(1, 3, 1, 0.5),
  y = c(0.5, 2, 6, -3)
)

n_sample <- 20000L

arith_results <- lapply(seq_len(nrow(validation_configs)), function(i) {
  mu <- validation_configs$mu[i]
  sigma <- validation_configs$sigma[i]
  y <- validation_configs$y[i]
  
  predictive_sample <- rnorm(n_sample, mean = mu, sd = sigma)
  
  rctest_crps <- compute_crps(
    forecast_density = predictive_sample,
    target_realization = y
  )
  
  reference_crps <- closed_form_gaussian_crps(y, mu, sigma)
  abs_diff <- abs(rctest_crps - reference_crps)
  
  data.frame(
    Config = validation_configs$Config[i],
    mu = mu,
    sigma = sigma,
    y = y,
    N_Sample = n_sample,
    RCtest_CRPS = rctest_crps,
    scoringRules_Closed_Form_CRPS = reference_crps,
    Abs_Diff = abs_diff,
    Rel_Diff_Pct = 100 * abs_diff / reference_crps,
    stringsAsFactors = FALSE
  )
})

arith_validation_table <- do.call(rbind, arith_results)

max_abs_diff <- max(arith_validation_table$Abs_Diff)
max_rel_diff <- max(arith_validation_table$Rel_Diff_Pct)
mean_rel_diff <- mean(arith_validation_table$Rel_Diff_Pct)

validation_summary <- data.frame(
  N_Configurations = nrow(arith_validation_table),
  Predictive_Sample_Size = n_sample,
  Maximum_Absolute_Difference = max_abs_diff,
  Maximum_Relative_Difference_Pct = max_rel_diff,
  Mean_Relative_Difference_Pct = mean_rel_diff,
  Interpretation = paste(
    "Monte Carlo approximation of analytic Gaussian CRPS;",
    "not a model-specific-density validation."
  ),
  stringsAsFactors = FALSE
)

cat("\n=== EMPIRICAL compute_crps() VS ANALYTIC scoringRules::crps_norm() ===\n")
cat("Predictive draws are generated from the known Gaussian used by the analytic\n")
cat("reference. This validates CRPS arithmetic for genuine predictive samples.\n\n")
print(arith_validation_table, row.names = FALSE, digits = 6)
cat("\n=== VALIDATION SUMMARY ===\n")
print(validation_summary, row.names = FALSE, digits = 6)

write.csv(
  arith_validation_table,
  file.path(output_dir, "crps_validation_against_scoringRules.csv"),
  row.names = FALSE
)

write.csv(
  validation_summary,
  file.path(output_dir, "crps_validation_against_scoringRules_summary.csv"),
  row.names = FALSE
)

realized_value <- 10

cross_model_point_forecasts <- c(
  7.0, 8.2, 9.4, 9.8,
  10.3, 10.9, 12.1, 13.5
)

genuine_predictive_sample <- rnorm(
  n_sample,
  mean = 10,
  sd = 1.5
)

proxy_demo <- data.frame(
  Input_Type = c(
    "Genuine predictive sample from one Gaussian model",
    "Cross-model point-forecast ensemble (disagreement proxy)"
  ),
  N_Values = c(length(genuine_predictive_sample), length(cross_model_point_forecasts)),
  Mean_Input = c(mean(genuine_predictive_sample), mean(cross_model_point_forecasts)),
  Input_SD = c(sd(genuine_predictive_sample), sd(cross_model_point_forecasts)),
  Mean_Absolute_Distance_To_Realized = c(
    mean(abs(genuine_predictive_sample - realized_value)),
    mean(abs(cross_model_point_forecasts - realized_value))
  ),
  RCtest_CRPS = c(
    compute_crps(genuine_predictive_sample, realized_value),
    compute_crps(cross_model_point_forecasts, realized_value)
  ),
  Interpretation = c(
    "Conventional sample-based CRPS for one specified predictive distribution.",
    "Ensemble-disagreement proxy; not a model-specific predictive distribution."
  ),
  stringsAsFactors = FALSE
)

cat("\n=== INPUT-INTERPRETATION DEMONSTRATION ===\n")
cat("The cross-model row is reported as an ensemble-disagreement proxy. Its CRPS\n")
cat("reflects both the point-forecast ensemble's location relative to the realized\n")
cat("value and its cross-sectional spread; it is not conventional model-specific\n")
cat("CRPS unless each model supplies its own predictive sample/distribution.\n\n")
print(proxy_demo, row.names = FALSE, digits = 6)

write.csv(
  proxy_demo,
  file.path(output_dir, "crps_ensemble_disagreement_proxy_demo.csv"),
  row.names = FALSE
)

relative_error_threshold_pct <- 5

if (max_rel_diff > relative_error_threshold_pct) {
  warning(
    "Empirical compute_crps() deviates from scoringRules::crps_norm() by ",
    round(max_rel_diff, 2), "% at most, exceeding the pre-specified ",
    relative_error_threshold_pct, "% Monte Carlo validation threshold."
  )
} else {
  message(sprintf(
    "CRPS arithmetic validation passed: maximum relative difference from ",
    "the analytic scoringRules reference is %.2f%% across %d configurations.",
    max_rel_diff,
    nrow(arith_validation_table)
  ))
}