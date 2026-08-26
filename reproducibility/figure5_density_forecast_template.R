# =============================================================================
# Figure 5 — Cross-sectional predictive density forecast
#
# The density is constructed from the 14 forecast models at one selected period.
# The cross-sectional spread is recentered to the cross-sectional mean, matching
# the empirical-density convention used for compute_crps().
# =============================================================================

library(RCtest)
library(ggplot2)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

forecast_model_cols <- seq_len(14L)
model_names <- colnames(metals)[forecast_model_cols]

selected_period <- 100L

if (selected_period < 1L || selected_period > nrow(metals)) {
  stop("selected_period must be between 1 and ", nrow(metals), ".")
}

ci_level <- 0.90

forecast_values <- as.numeric(
  metals[selected_period, forecast_model_cols]
)

cross_sectional_mean <- mean(
  forecast_values,
  na.rm = TRUE
)

full_distribution <- forecast_values -
  cross_sectional_mean +
  cross_sectional_mean

point_forecast <- cross_sectional_mean

p <- plot_density_forecast(
  full_distribution = full_distribution,
  point_forecast = point_forecast,
  title = NULL,
  ci_level = ci_level
) +
  labs(
    title = NULL,
    x = "Forecast Value",
    y = "Density"
  ) +
  theme(
    plot.title = element_blank()
  )

print(p)

figure_file <- file.path(
  output_dir,
  "figure5_cross_sectional_predictive_density.png"
)

ggsave(
  filename = figure_file,
  plot = p,
  width = 9,
  height = 5.8,
  units = "in",
  dpi = 300
)

density_data <- data.frame(
  Period = selected_period,
  Model = model_names,
  Forecast = forecast_values,
  Cross_Sectional_Mean = cross_sectional_mean,
  Recentered_Draw = full_distribution,
  stringsAsFactors = FALSE
)

write.csv(
  density_data,
  file.path(
    output_dir,
    "figure5_cross_sectional_predictive_density_data.csv"
  ),
  row.names = FALSE
)

ci_bounds <- quantile(
  full_distribution,
  probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
  na.rm = TRUE
)

metadata <- data.frame(
  Period = selected_period,
  Models = length(model_names),
  Point_Forecast = point_forecast,
  Reference = "Cross-sectional mean of 14 model forecasts",
  Confidence_Interval_Level = ci_level,
  Lower_CI = unname(ci_bounds[1]),
  Upper_CI = unname(ci_bounds[2]),
  Distribution_Construction = paste(
    "Empirical cross-sectional distribution of 14 forecasts,",
    "recentered to the cross-sectional mean"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    output_dir,
    "figure5_cross_sectional_predictive_density_metadata.csv"
  ),
  row.names = FALSE
)

cat(
  "Figure 5 saved to: ",
  normalizePath(figure_file),
  "\n",
  sep = ""
)