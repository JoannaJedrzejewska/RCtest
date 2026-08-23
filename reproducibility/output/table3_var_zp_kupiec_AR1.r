# =============================================================================
# table3_var_zp_kupiec_AR1.R
#
# Recreates Table 3: VaR coverage, Kupiec backtesting, and mean ZP loss,
# using AR_1 as the benchmark. This replaces the old manuscript Table 3,
# which incorrectly listed AR_1 as a competing model while AR_1 is now the
# benchmark in revised Tables 1 and 2.
#
# Outputs (written to working directory):
# - table3_var_zp_kupiec_AR1.csv clean combined manuscript table
# - table3_kupiec_raw_AR1.rds raw output for audit/reproducibility
# - table3_zp_loss_raw_AR1.csv period-level ZP loss matrix
# =============================================================================

library(RCtest)
data(metals)

P <- nrow(metals)
K_total <- ncol(metals)
all_model_names <- colnames(metals)[1:14]
benchmark_name <- "AR_1"
benchmark_idx <- which(all_model_names == benchmark_name)

if (length(benchmark_idx) != 1L) {
  stop("Could not uniquely identify AR_1 among the first 14 model columns.")
}

competitor_idx <- setdiff(seq_len(14), benchmark_idx)
competitor_names <- all_model_names[competitor_idx]
realized <- metals[, K_total]
window_size <- 20
alpha <- 0.05
zp_quantile <- 0.05

forecast_variance <- estimate_forecast_variance(
  metals,
  realized = realized,
  benchmark_col = benchmark_idx,
  window_size = window_size
)

all_nonbenchmark_cols <- setdiff(seq_len(K_total), benchmark_idx)
forecast_sd_all <- sqrt(forecast_variance[, all_nonbenchmark_cols, drop = FALSE])

forecast_matrix <- metals
threshold <- quantile(realized, zp_quantile, na.rm = TRUE)


kupiec_raw <- compute_kupiec(
  forecast_matrix,
  forecast_sd_all,
  benchmark_col = benchmark_idx,
  alpha = alpha
)
saveRDS(kupiec_raw, "table3_kupiec_raw_AR1.rds")

cat("=== RAW KUPIEC OUTPUT ===\n")
print(kupiec_raw)
cat("\n=== RAW KUPIEC STRUCTURE ===\n")
str(kupiec_raw)

zp_loss_all <- compute_zp(
  forecast_matrix,
  forecast_sd_all,
  threshold = threshold,
  benchmark_col = benchmark_idx
)
write.csv(zp_loss_all, "table3_zp_loss_raw_AR1.csv", row.names = FALSE)

mean_zp_loss <- colMeans(zp_loss_all[, competitor_idx, drop = FALSE], na.rm = TRUE)

sd_col_for_forecast_col <- match(competitor_idx, all_nonbenchmark_cols)
model_sd <- forecast_sd_all[, sd_col_for_forecast_col, drop = FALSE]
model_forecasts <- forecast_matrix[, competitor_idx, drop = FALSE]
model_var <- sweep(model_forecasts, 2, qnorm(alpha), FUN = function(mu, z) mu + z * 0) # placeholder

model_var <- model_forecasts + qnorm(alpha) * model_sd
violations <- colSums(sweep(model_var, 1, realized, FUN = function(var_value, y) y < var_value), na.rm = TRUE)
expected_violations <- P * alpha

table3 <- data.frame(
  Model = competitor_names,
  Violations = as.integer(violations),
  Expected = rep(expected_violations, length(competitor_names)),
  Kupiec_p_value = NA_real_,
  Mean_ZP_Loss = as.numeric(mean_zp_loss),
  stringsAsFactors = FALSE
)

if (is.data.frame(kupiec_raw)) {
  possible_p_columns <- grep("p.?value|p_value", names(kupiec_raw),
                             ignore.case = TRUE, value = TRUE)
  possible_model_columns <- grep("model", names(kupiec_raw),
                                 ignore.case = TRUE, value = TRUE)
  
  if (length(possible_p_columns) >= 1L && length(possible_model_columns) >= 1L) {
    kupiec_df <- kupiec_raw
    names(kupiec_df)[names(kupiec_df) == possible_model_columns[1]] <- "Model"
    names(kupiec_df)[names(kupiec_df) == possible_p_columns[1]] <- "Kupiec_p_value"
    table3 <- merge(
      table3[, c("Model", "Violations", "Expected", "Mean_ZP_Loss")],
      kupiec_df[, c("Model", "Kupiec_p_value")],
      by = "Model", all.x = TRUE, sort = FALSE
    )
    table3 <- table3[, c("Model", "Violations", "Expected", "Kupiec_p_value", "Mean_ZP_Loss")]
  }
}

table3$Model <- factor(table3$Model, levels = competitor_names)
table3 <- table3[order(table3$Model), ]
table3$Model <- as.character(table3$Model)

cat("\n=== TABLE 3: AR_1 BENCHMARK ===\n")
print(table3, row.names = FALSE)
write.csv(table3, "table3_var_zp_kupiec_AR1.csv", row.names = FALSE)

cat("\nIMPORTANT:\n")
cat("1. Inspect the RAW KUPIEC OUTPUT printed above.\n")
cat("2. If Kupiec_p_value remains NA in the clean table, paste the raw output\n")
cat(" here before putting Table 3 in the manuscript; do not guess the field.\n")
cat("3. AR_1 must NOT appear as a table row: it is the benchmark.\n")