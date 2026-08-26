# =============================================================================
# Figure 5 — Model specific Gaussian predictive distribution
#
# Illustrates a  predictive distribution for BDMM-SS at one evaluation
# period. The distribution is a Gaussian approximation with:
#   mean = the model's point forecast;
#   SD   = its rolling error SD estimated from past realized errors.
# =============================================================================

library(RCtest)
library(ggplot2)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

forecast_model_cols <- seq_len(14L)
benchmark_name <- "AR_1"
model_name <- "BDMM-SS"
realized_name <- "HA"

selected_period <- 100L
window_size <- 20L
interval_level <- 0.90

benchmark_col <- match(benchmark_name, colnames(metals))
model_col <- match(model_name, colnames(metals))
realized_col <- match(realized_name, colnames(metals))

if (anyNA(c(benchmark_col, model_col, realized_col))) {
  stop("Could not locate model, benchmark, or realized series in metals.")
}

if (selected_period < 1L || selected_period > nrow(metals)) {
  stop("selected_period must be between 1 and ", nrow(metals), ".")
}

realized <- metals[, realized_col]

forecast_variance <- estimate_forecast_variance(
  forecast_matrix = metals,
  realized = realized,
  benchmark_col = benchmark_col,
  window_size = window_size
)

model_sd <- sqrt(forecast_variance[selected_period, model_col])

if (!is.finite(model_sd) || model_sd <= 0) {
  stop("Model-specific predictive SD must be finite and positive.")
}

point_forecast <- metals[selected_period, model_col]
realized_value <- realized[selected_period]

set.seed(20260824)

n_predictive_draws <- 10000L

predictive_draws <- rnorm(
  n_predictive_draws,
  mean = point_forecast,
  sd = model_sd
)

p <- plot_density_forecast(
  full_distribution = predictive_draws,
  point_forecast = point_forecast,
  title = NULL,
  ci_level = interval_level
) +
  geom_vline(
    xintercept = realized_value,
    linetype = "solid",
    linewidth = 0.7,
    colour = "black"
  ) +
  labs(
    title = NULL,
    x = "Forecast value",
    y = "Predictive density"
  ) +
  theme(
    plot.title = element_blank()
  )

print(p)

figure_file <- file.path(
  output_dir,
  "figure5_bdmm_ss_predictive_density.png"
)

ggsave(
  filename = figure_file,
  plot = p,
  width = 9,
  height = 5.8,
  units = "in",
  dpi = 300
)

alpha_lower <- (1 - interval_level) / 2
interval_bounds <- qnorm(
  c(alpha_lower, 1 - alpha_lower),
  mean = point_forecast,
  sd = model_sd
)

metadata <- data.frame(
  Period = selected_period,
  Model = model_name,
  Benchmark = benchmark_name,
  Realized_Series = realized_name,
  Point_Forecast = point_forecast,
  Realized_Value = realized_value,
  Predictive_SD = model_sd,
  Interval_Level = interval_level,
  Lower_Predictive_Interval = interval_bounds[1],
  Upper_Predictive_Interval = interval_bounds[2],
  Distribution = "Gaussian approximation using model-specific rolling forecast-error SD",
  Window_Size = window_size,
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    output_dir,
    "figure5_bdmm_ss_predictive_density_metadata.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Predictive_Draw = predictive_draws
  ),
  file.path(
    output_dir,
    "figure5_bdmm_ss_predictive_density_draws.csv"
  ),
  row.names = FALSE
)

cat(
  "Figure 5 saved to: ",
  normalizePath(figure_file),
  "\n",
  sep = ""
)