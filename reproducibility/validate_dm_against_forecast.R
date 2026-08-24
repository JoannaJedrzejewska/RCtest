# =============================================================================
# Direct numerical validation of compute_per_model_statistics() against the
# independently implemented forecast::dm.test().
#
# For each competing model, run RCtest's per-model DM statistic and
# forecast::dm.test() on IDENTICAL forecast-error series, under an explicit,
# documented sign convention, and report both the raw agreement and a
# pass/fail flag per model.
#
# Output:
#   reproducibility/output/dm_validation_against_forecast.csv
# =============================================================================

library(RCtest)

if (!requireNamespace("forecast", quietly = TRUE)) {
  stop("Package 'forecast' is required for this validation script.")
}

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

forecast_model_cols <- seq_len(14L)
realized_col <- match("HA", colnames(metals))
benchmark_name <- "AR_1"
benchmark_col <- match(benchmark_name, colnames(metals))

if (is.na(realized_col)) stop("HA was not found as the realised series.")
if (is.na(benchmark_col) || !benchmark_col %in% forecast_model_cols) {
  stop("AR_1 must be one of the first 14 forecast-model columns.")
}

competitor_cols <- setdiff(forecast_model_cols, benchmark_col)
competitor_names <- colnames(metals)[competitor_cols]

realized <- metals[, realized_col]

# -----------------------------------------------------------------------------
# Sign-convention contract (must match the definition used in table1_AR1.R /
# table2_AR1.R): loss_difference = benchmark_loss - competitor_loss, so that a
# POSITIVE value means the competitor is superior (lower loss than AR_1).
#
# forecast::dm.test(e1, e2, alternative, h, power, varestimator):
#   e1, e2 are FORECAST ERRORS (forecast - realized), not losses.
#   alternative = "less"    -> H1: method 2 (e2) is LESS accurate than method 1 (e1)
#   alternative = "greater" -> H1: method 2 (e2) is MORE accurate than method 1 (e1)
#
# To test "competitor beats benchmark" (our RCtest convention's implicit H1),
# set e1 = benchmark errors, e2 = competitor errors, alternative = "greater".
# A significant result under alternative = "greater" means the competitor (e2)
# is more accurate than the benchmark (e1) -- i.e. loss_difference > 0 side.
# -----------------------------------------------------------------------------

block_length <- 5L
n_boot <- 999L
alpha <- 0.05

benchmark_error <- metals[, benchmark_col] - realized

results <- lapply(competitor_names, function(model_name) {
  model_col <- match(model_name, colnames(metals))
  competitor_error <- metals[, model_col] - realized
  
  loss_diff <- matrix(
    benchmark_error^2 - competitor_error^2,
    ncol = 1,
    dimnames = list(NULL, model_name)
  )
  
  set.seed(20260824)
  rctest_result <- compute_per_model_statistics(
    loss_differences = loss_diff,
    model_names = model_name,
    n_boot = n_boot,
    block_length = block_length,
    alpha = alpha,
    h = 1,
    H1 = "same"
  )
  
  dm_two_sided <- forecast::dm.test(
    e1 = benchmark_error,
    e2 = competitor_error,
    alternative = "two.sided",
    h = 1,
    power = 2,
    varestimator = "acf"
  )
  
  dm_greater <- forecast::dm.test(
    e1 = benchmark_error,
    e2 = competitor_error,
    alternative = "greater",
    h = 1,
    power = 2,
    varestimator = "acf"
  )
  
  rctest_mean_diff <- rctest_result$Mean_Loss_Diff[1]
  rctest_t_stat <- rctest_result$T_Stat[1]
  rctest_p <- rctest_result[[grep("^P_Value", names(rctest_result))[1]]][1]
  
  sign_agrees <- sign(rctest_mean_diff) == sign(dm_greater$statistic) ||
    (rctest_mean_diff > 0 && dm_greater$statistic > 0) ||
    (rctest_mean_diff < 0 && dm_greater$statistic < 0)
  
  data.frame(
    Model = model_name,
    RCtest_Mean_Loss_Diff = rctest_mean_diff,
    RCtest_T_Stat = rctest_t_stat,
    RCtest_P_Value = rctest_p,
    DM_forecast_Statistic_TwoSided = unname(dm_two_sided$statistic),
    DM_forecast_P_Value_TwoSided = dm_two_sided$p.value,
    DM_forecast_Statistic_Greater = unname(dm_greater$statistic),
    DM_forecast_P_Value_Greater = dm_greater$p.value,
    Sign_Agreement = sign_agrees,
    T_Stat_Abs_Diff = abs(abs(rctest_t_stat) - abs(unname(dm_two_sided$statistic))),
    stringsAsFactors = FALSE
  )
})

validation_table <- do.call(rbind, results)

n_disagree <- sum(!validation_table$Sign_Agreement)

cat("\n=== DM VALIDATION AGAINST forecast::dm.test() ===\n")
print(validation_table, row.names = FALSE, digits = 4)
cat(sprintf(
  "\n%d of %d models show sign disagreement between RCtest and forecast::dm.test().\n",
  n_disagree, nrow(validation_table)
))

if (n_disagree > 0L) {
  warning(
    "Sign disagreement detected between compute_per_model_statistics() and ",
    "forecast::dm.test() for ", n_disagree, " model(s). Investigate the sign ",
    "convention in compute_per_model_statistics() before treating Table 2 as valid."
  )
} else {
  message("All models agree in sign between RCtest and forecast::dm.test().")
}

write.csv(
  validation_table,
  file.path(output_dir, "dm_validation_against_forecast.csv"),
  row.names = FALSE
)