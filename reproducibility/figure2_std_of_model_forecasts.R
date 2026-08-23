# =============================================================================
# Figure 2: Standard deviation of model forecasts over time
# =============================================================================

library(RCtest)
library(ggplot2)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

forecast_models <- colnames(metals)[1:14]

figure2_data <- data.frame(
  Period = seq_len(nrow(metals)),
  Cross_Sectional_SD = apply(
    metals[, forecast_models, drop = FALSE],
    1,
    sd,
    na.rm = TRUE
  )
)

p <- ggplot(figure2_data, aes(x = Period, y = Cross_Sectional_SD)) +
  geom_area(fill = "#CFE8F3", alpha = 0.75) +
  geom_line(colour = "#2C7FB8", linewidth = 0.70) +
  labs(
    x = "Period",
    y = "Cross-Sectional SD"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  file.path(output_dir, "figure2_cross_sectional_disagreement.png"),
  p,
  width = 10,
  height = 7,
  dpi = 300
)

write.csv(
  figure2_data,
  file.path(output_dir, "figure2_cross_sectional_disagreement_data.csv"),
  row.names = FALSE
)