# =============================================================================
# Figure 3 -Cumulative MSE loss difference relative to AR(1) benchmark
# =============================================================================


library(RCtest)
library(ggplot2)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

P <- nrow(metals)
K_total <- ncol(metals)

forecast_model_cols <- seq_len(14L)
all_model_names <- colnames(metals)[forecast_model_cols]

benchmark_name <- "AR_1"

benchmark_idx_in_models <- match(
  benchmark_name,
  all_model_names
)

if (is.na(benchmark_idx_in_models)) {
  stop(
    "Could not find AR_1 among the first 14 forecast columns. ",
    "Actual names: ",
    paste(all_model_names, collapse = ", ")
  )
}

realized <- metals[, "HA"]

errors <- sweep(
  as.matrix(metals[, forecast_model_cols, drop = FALSE]),
  1,
  realized,
  FUN = "-"
)

colnames(errors) <- all_model_names

fig3_ar1 <- plot_cumulative_loss(
  data_matrix = errors,
  benchmark_col = benchmark_idx_in_models
)

# Remove only the title supplied by plot_cumulative_loss().
fig3_ar1 <- fig3_ar1 +
  labs(title = NULL) +
  theme(
    plot.title = element_blank()
  )

print(fig3_ar1)

figure_file <- file.path(
  output_dir,
  "figure3_cumulative_loss_AR1_benchmark.png"
)

ggsave(
  filename = figure_file,
  plot = fig3_ar1,
  width = 6.30,
  height = 4.92,
  units = "in",
  dpi = 300
)

cat("Figure 3 saved to:", normalizePath(figure_file), "\n")
cat(
  "Benchmark:",
  benchmark_name,
  "(forecast-column position",
  benchmark_idx_in_models,
  ")\n"
)