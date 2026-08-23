# =============================================================================
# Figure 1: Forecasts from all competing models and the benchmark
# =============================================================================

library(RCtest)
library(ggplot2)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

forecast_models <- colnames(metals)[1:14]
benchmark_name <- "AR_1"
periods <- seq_len(nrow(metals))

plot_data <- data.frame(
  Period = rep(periods, length(forecast_models)),
  Forecast = as.vector(metals[, forecast_models, drop = FALSE]),
  Model = rep(forecast_models, each = length(periods))
)

p <- ggplot() +
  geom_line(
    data = subset(plot_data, Model != benchmark_name),
    aes(x = Period, y = Forecast, group = Model),
    colour = "#A6CEE3",
    linewidth = 0.35,
    alpha = 0.70
  ) +
  geom_line(
    data = subset(plot_data, Model == benchmark_name),
    aes(x = Period, y = Forecast),
    colour = "#B63A2B",
    linewidth = 0.95
  ) +
  labs(
    x = "Period",
    y = "Price Index"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p <- p +
  annotate("segment", x = 60, xend = 66, y = 145.8, yend = 145.8,
           colour = "#B63A2B", linewidth = 1.1) +
  annotate("text", x = 67, y = 145.8, label = "Benchmark (AR_1)",
           hjust = 0, size = 4) +
  annotate("segment", x = 98, xend = 104, y = 145.8, yend = 145.8,
           colour = "#A6CEE3", linewidth = 0.8) +
  annotate("text", x = 105, y = 145.8, label = "Competing",
           hjust = 0, size = 4)

print(p)

ggsave(
  file.path(output_dir, "figure1_forecast_comparison_AR1.png"),
  p,
  width = 11,
  height = 7,
  dpi = 300
)