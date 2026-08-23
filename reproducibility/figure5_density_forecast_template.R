# =============================================================================
# Figure 5: Predictive density forecast.
# =============================================================================

library(RCtest)
library(ggplot2)

set.seed(20260822)
data(metals)
output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Choose the forecast model and the evaluation-period forecast displayed.
selected_model <- "BDMM-K"
selected_period <- nrow(metals)
benchmark_name <- "AR_1"
n_draws <- 10000

if (!selected_model %in% colnames(metals)) stop("selected_model not found.")
if (!benchmark_name %in% colnames(metals)) stop("AR_1 benchmark not found.")

# Forecast residuals are defined relative to the HA reference series, which is
# column 15 of metals. The empirical residual SD supplies predictive spread.
reference_series <- metals[, "HA"]
residuals <- metals[, selected_model] - reference_series
residual_sd <- sd(residuals, na.rm = TRUE)
point_forecast <- metals[selected_period, selected_model]

# Documented Gaussian predictive distribution: N(point forecast, residual SD^2).
full_distribution <- rnorm(n_draws, mean = point_forecast, sd = residual_sd)

p <- plot_density_forecast(
  full_distribution = full_distribution,
  point_forecast = point_forecast,
  title = paste("Predictive Density Forecast:", selected_model),
  ci_level = 0.90
)
print(p)
ggsave(file.path(output_dir, "figure5_density_forecast.png"),
       plot = p, width = 9, height = 5.8, dpi = 300)

write.csv(
  data.frame(Value = full_distribution),
  file.path(output_dir, "figure5_density_forecast_draws.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    Model = selected_model,
    Period = selected_period,
    Point_Forecast = point_forecast,
    Reference_Series = "HA",
    Residual_SD = residual_sd,
    Distribution = "Gaussian N(point forecast, residual SD^2)",
    Seed = 20260822,
    Draws = n_draws,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "figure5_density_forecast_metadata.csv"),
  row.names = FALSE
)
