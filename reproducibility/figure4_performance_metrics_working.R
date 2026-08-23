# =============================================================================
# Figure 4: Performance metrics and ERC weights
# =============================================================================

library(RCtest)
library(grid)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

K_models <- ncol(metals) - 1L
custom_weights <- (K_models:1) / sum(K_models:1)

p <- plot_performance_metrics(
  forecast_matrix = metals,
  weights = custom_weights,
  benchmark_col = 15
)
grid.newpage()
grid.draw(p)

png(
  filename = file.path(output_dir, "figure4_performance_metrics.png"),
  width = 3000,
  height = 2200,
  res = 300
)
grid.newpage()
grid.draw(p)
dev.off()

# Save the values used in the four panels.
realized <- metals[, 15]
models <- metals[, 1:14, drop = FALSE]
errors <- sweep(models, 1, realized, "-")
rmse_vals <- sqrt(colMeans(errors^2, na.rm = TRUE))
mae_vals <- colMeans(abs(errors), na.rm = TRUE)
nrmse_vals <- rmse_vals / mean(abs(realized), na.rm = TRUE)
naive_mae <- mean(abs(diff(realized)), na.rm = TRUE)
mase_vals <- mae_vals / naive_mae

figure4_data <- data.frame(
  Model = colnames(models),
  RMSE = as.numeric(rmse_vals),
  NRMSE = as.numeric(nrmse_vals),
  MAE = as.numeric(mae_vals),
  MASE = as.numeric(mase_vals),
  ERC_Weight = as.numeric(custom_weights),
  Risk_Contribution = as.numeric(rmse_vals * custom_weights),
  stringsAsFactors = FALSE
)

write.csv(
  figure4_data,
  file.path(output_dir, "figure4_performance_metrics_data.csv"),
  row.names = FALSE
)
