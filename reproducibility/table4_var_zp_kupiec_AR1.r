# =============================================================================
# Table 4: VaR coverage, Kupiec UC backtests, and mean ZP loss.
#
# AR_1 is the forecast-comparison benchmark.
# HA is the realised outcome series.
#
# Outputs:
#   reproducibility/output/table4_var_zp_kupiec_AR1_raw.csv
#   reproducibility/output/table4_var_zp_kupiec_AR1.csv
#   reproducibility/output/table4_kupiec_raw_AR1.rds
#   reproducibility/output/table4_zp_loss_raw_AR1.csv
# =============================================================================
library(RCtest)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

P <- nrow(metals)

forecast_model_cols <- seq_len(14L)

benchmark_name <- "AR_1"
realized_name <- "HA"

benchmark_col <- match(benchmark_name, colnames(metals))
realized_col <- match(realized_name, colnames(metals))

if (is.na(benchmark_col) || !benchmark_col %in% forecast_model_cols) {
  stop("AR_1 must be one of the first 14 forecast-model columns.")
}

if (is.na(realized_col)) {
  stop("HA was not found as the realised outcome column.")
}

competitor_cols <- setdiff(forecast_model_cols, benchmark_col)
competitor_names <- colnames(metals)[competitor_cols]

if (length(competitor_cols) != 14L) {
  stop("Expected 13 competing forecasts after excluding AR_1.")
}

realized <- metals[, realized_col]

window_size <- 20L
alpha <- 0.05
zp_quantile <- 0.05

forecast_variance <- estimate_forecast_variance(
  forecast_matrix = metals,
  realized = realized,
  benchmark_col = benchmark_col,
  window_size = window_size
)

nonbenchmark_cols <- setdiff(seq_len(ncol(metals)), benchmark_col)

forecast_sd_models <- sqrt(
  forecast_variance[, nonbenchmark_cols, drop = FALSE]
)

if (ncol(forecast_sd_models) != ncol(metals) - 1L) {
  stop("Expected one predictive-SD column for each non-AR_1 series.")
}

threshold_value <- quantile(
  realized,
  probs = zp_quantile,
  na.rm = TRUE
)

kupiec_raw <- compute_kupiec(
  forecast_matrix = metals,
  forecast_sd_models = forecast_sd_models,
  realized = realized,
  benchmark_col = benchmark_col,
  alpha = alpha
)

saveRDS(
  kupiec_raw,
  file.path(output_dir, "table4_kupiec_raw_AR1.rds")
)

cat("=== RAW KUPIEC OUTPUT ===\n")
print(kupiec_raw)

kupiec_table <- do.call(
  rbind,
  lapply(names(kupiec_raw), function(model_name) {
    model_result <- kupiec_raw[[model_name]]
    
    data.frame(
      Model = model_name,
      Violations = unname(model_result$actual_exceedances),
      Expected = unname(model_result$expected),
      Kupiec_LR_UC = unname(model_result$statistic),
      Kupiec_p_value = unname(model_result$p.value),
      stringsAsFactors = FALSE
    )
  })
)

row.names(kupiec_table) <- NULL

kupiec_table <- kupiec_table[
  kupiec_table$Model %in% competitor_names,
  ,
  drop = FALSE
]

zp_loss_all <- compute_zp(
  forecast_matrix = metals,
  forecast_sd_models = forecast_sd_models,
  threshold = threshold_value,
  benchmark_col = benchmark_col
)

write.csv(
  zp_loss_all,
  file.path(output_dir, "table4_zp_loss_raw_AR1.csv"),
  row.names = FALSE
)

if (
  !is.null(colnames(zp_loss_all)) &&
  all(competitor_names %in% colnames(zp_loss_all))
) {
  mean_zp_loss <- colMeans(
    zp_loss_all[, competitor_names, drop = FALSE],
    na.rm = TRUE
  )
} else {
  mean_zp_loss <- colMeans(
    zp_loss_all[, competitor_cols, drop = FALSE],
    na.rm = TRUE
  )
  
  names(mean_zp_loss) <- competitor_names
}

table4_raw <- merge(
  kupiec_table,
  data.frame(
    Model = names(mean_zp_loss),
    Mean_ZP_Loss = as.numeric(mean_zp_loss),
    stringsAsFactors = FALSE
  ),
  by = "Model",
  all.x = TRUE,
  sort = FALSE
)

table4_raw$Model <- factor(
  table4_raw$Model,
  levels = competitor_names
)

table4_raw <- table4_raw[
  order(table4_raw$Model),
  ,
  drop = FALSE
]

table4_raw$Model <- as.character(table4_raw$Model)

write.csv(
  table4_raw,
  file.path(output_dir, "table4_var_zp_kupiec_AR1_raw.csv"),
  row.names = FALSE
)

format_p_value <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }
  
  if (p < 0.001) {
    return("< 0.001")
  }
  
  sprintf("%.3f", p)
}

table4_manuscript <- table4_raw

table4_manuscript$Kupiec_LR_UC <- sprintf(
  "%.4f",
  table4_manuscript$Kupiec_LR_UC
)

table4_manuscript$Kupiec_p_value <- vapply(
  table4_manuscript$Kupiec_p_value,
  format_p_value,
  character(1)
)

table4_manuscript$Mean_ZP_Loss <- sprintf(
  "%.6f",
  table4_manuscript$Mean_ZP_Loss
)

cat("\n=== TABLE 4: AR_1 BENCHMARK ===\n")
print(table4_manuscript, row.names = FALSE)

write.csv(
  table4_manuscript,
  file.path(output_dir, "table4_var_zp_kupiec_AR1.csv"),
  row.names = FALSE
)
