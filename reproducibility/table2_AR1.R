# =============================================================================
# Table 2: Per-model forecast-comparison statistics using AR_1 as benchmark.
#
# HA is the realised outcome series.
# The 13 forecast methods other than AR_1 are evaluated against AR_1.
#
# Outputs:
#   reproducibility/output/table2_AR1_raw.csv
#   reproducibility/output/table2_AR1.csv
# =============================================================================

library(RCtest)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260823)

forecast_model_cols <- seq_len(14L)
realized_col <- match("HA", colnames(metals))

benchmark_name <- "AR_1"
benchmark_col <- match(benchmark_name, colnames(metals))

if (is.na(realized_col)) {
  stop("HA was not found as the realised outcome series.")
}

if (is.na(benchmark_col) || !benchmark_col %in% forecast_model_cols) {
  stop("AR_1 must be one of the first 14 forecast-model columns.")
}

competitor_cols <- setdiff(forecast_model_cols, benchmark_col)
competitor_names <- colnames(metals)[competitor_cols]

if (length(competitor_cols) != 13L) {
  stop("Expected 13 competing forecast models after excluding AR_1.")
}

realized <- metals[, realized_col]

benchmark_loss <- (metals[, benchmark_col] - realized)^2

competitor_loss <- sweep(
  metals[, competitor_cols, drop = FALSE],
  1,
  realized,
  FUN = function(forecast, outcome) {
    (forecast - outcome)^2
  }
)

# Positive values mean the competing model has lower squared error than AR_1.
loss_differences <- sweep(
  competitor_loss,
  1,
  benchmark_loss,
  FUN = function(model_loss, benchmark_loss_value) {
    benchmark_loss_value - model_loss
  }
)

colnames(loss_differences) <- competitor_names

block_length <- 5L
n_boot <- 999L
alpha <- 0.05

# H1 = "same" gives the manuscript two-sided per-model comparison.
table2_raw <- compute_per_model_statistics(
  loss_differences = loss_differences,
  model_names = competitor_names,
  n_boot = n_boot,
  block_length = block_length,
  alpha = alpha,
  h = 1,
  H1 = "same"
)


if ("Frac_Better" %in% names(table2_raw)) {
  names(table2_raw)[
    names(table2_raw) == "Frac_Better"
  ] <- "Frac_Better_Than_AR_1"
}

table2_raw <- table2_raw[
  order(table2_raw$Mean_Loss_Diff, decreasing = TRUE),
  ,
  drop = FALSE
]

write.csv(
  table2_raw,
  file.path(output_dir, "table2_AR1_raw.csv"),
  row.names = FALSE
)

format_bootstrap_p <- function(p, n_boot) {
  lower_bound <- 1 / (n_boot + 1)
  upper_bound <- 1 - lower_bound
  
  if (is.na(p)) {
    return(NA_character_)
  }
  
  if (p <= lower_bound) {
    return(sprintf("< %.3f", lower_bound))
  }
  
  if (p >= upper_bound) {
    return(sprintf(">= %.3f", upper_bound))
  }
  
  sprintf("%.3f", p)
}

table2_manuscript <- table2_raw

numeric_six_columns <- intersect(
  c("Mean_Loss_Diff"),
  names(table2_manuscript)
)

for (column_name in numeric_six_columns) {
  table2_manuscript[[column_name]] <- sprintf(
    "%.6f",
    table2_manuscript[[column_name]]
  )
}

numeric_four_columns <- intersect(
  c("T_Stat"),
  names(table2_manuscript)
)

for (column_name in numeric_four_columns) {
  table2_manuscript[[column_name]] <- sprintf(
    "%.4f",
    table2_manuscript[[column_name]]
  )
}

fraction_columns <- grep(
  "^Frac_Better",
  names(table2_manuscript),
  value = TRUE
)

for (column_name in fraction_columns) {
  table2_manuscript[[column_name]] <- sprintf(
    "%.3f",
    table2_manuscript[[column_name]]
  )
}

p_columns <- grep(
  "^P_Value",
  names(table2_manuscript),
  value = TRUE
)

for (column_name in p_columns) {
  table2_manuscript[[column_name]] <- vapply(
    table2_raw[[column_name]],
    format_bootstrap_p,
    character(1),
    n_boot = n_boot
  )
}

print(table2_manuscript, row.names = FALSE)

write.csv(
  table2_manuscript,
  file.path(output_dir, "table2_AR1.csv"),
  row.names = FALSE
)
