# =============================================================================
# runtime_benchmarks.R
#
# Reproduces:
#   Table 4: Runtime scaling by model count and bootstrap replications.
#   Table 5: Runtime by core RCtest function.
#   Memory-allocation profile for one full WRC/SPA/CPA battery.
# =============================================================================


library(RCtest)
library(microbenchmark)

data(metals)

output_dir <- file.path("reproducibility", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== SESSION / VERSION INFO ===\n")
cat("R version:", R.version.string, "\n")

packages_to_report <- c(
  "RCtest",
  "microbenchmark",
  "ggplot2",
  "gridExtra",
  "ggrepel",
  "rlang"
)

installed_info <- installed.packages()

available_packages <- packages_to_report[
  packages_to_report %in% rownames(installed_info)
]

print(installed_info[available_packages, "Version"])

cat("Random seed: 20260822\n")
cat("Forecast-comparison benchmark: AR_1\n")
cat("Realised outcome series: HA\n\n")

write.csv(
  data.frame(
    Package = available_packages,
    Version = installed_info[available_packages, "Version"],
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "runtime_package_versions.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste("R version:", R.version.string),
    "Parallel computing: NONE (single-threaded)",
    "Random seed: 20260822",
    "Forecast-comparison benchmark: AR_1",
    "Realised outcome series: HA",
    paste("Date:", Sys.Date())
  ),
  con = file.path(output_dir, "runtime_session_information.txt")
)

set.seed(20260822)

P <- nrow(metals)

forecast_model_cols <- seq_len(14L)

realized_name <- "HA"
benchmark_name <- "AR_1"

realized_col <- match(realized_name, colnames(metals))
benchmark_col <- match(benchmark_name, colnames(metals))

if (is.na(realized_col)) {
  stop("HA was not found in the metals data.")
}

if (is.na(benchmark_col) || !benchmark_col %in% forecast_model_cols) {
  stop("AR_1 must be one of the first 14 forecast-model columns.")
}

competitor_cols <- setdiff(forecast_model_cols, benchmark_col)
competitor_names <- colnames(metals)[competitor_cols]

if (length(competitor_cols) != 13L) {
  stop("Expected 13 forecast competitors after excluding AR_1.")
}

realized <- metals[, realized_col]

# Benchmark and competitor losses are all calculated against realised HA.
benchmark_loss <- (metals[, benchmark_col] - realized)^2

competitor_loss <- sweep(
  metals[, competitor_cols, drop = FALSE],
  1,
  realized,
  FUN = function(forecast, outcome) {
    (forecast - outcome)^2
  }
)

# Positive values mean a competing model has lower squared loss than AR_1.
loss_diff <- sweep(
  competitor_loss,
  1,
  benchmark_loss,
  FUN = function(model_loss, benchmark_loss_value) {
    benchmark_loss_value - model_loss
  }
)

colnames(loss_diff) <- competitor_names

forecast_variance <- estimate_forecast_variance(
  forecast_matrix = metals,
  realized = realized,
  benchmark_col = benchmark_col,
  window_size = 20
)

nonbenchmark_cols <- setdiff(seq_len(ncol(metals)), benchmark_col)

forecast_sd_models <- sqrt(
  forecast_variance[, nonbenchmark_cols, drop = FALSE]
)

if (ncol(forecast_sd_models) != ncol(metals) - 1L) {
  stop("Expected one SD column for each non-AR_1 data column.")
}

threshold_val <- quantile(
  realized,
  probs = 0.05,
  na.rm = TRUE
)

weighting_vector <- abs(realized)

klic_loss <- compute_klic(
  forecast_matrix = metals,
  forecast_sd_models = forecast_sd_models,
  benchmark_col = benchmark_col
)

log_lik_diff <- klic_loss[, benchmark_col] -
  klic_loss[, competitor_cols, drop = FALSE]

colnames(log_lik_diff) <- competitor_names

zp_loss <- compute_zp(
  forecast_matrix = metals,
  forecast_sd_models = forecast_sd_models,
  threshold = threshold_val,
  benchmark_col = benchmark_col
)

zp_diff <- zp_loss[, benchmark_col] -
  zp_loss[, competitor_cols, drop = FALSE]

colnames(zp_diff) <- competitor_names

cat("Data preparation is OUTSIDE all timed blocks.\n")
cat("Bootstrap RNG inside each test call is included in runtime.\n")
cat("Table 4/5 empirical benchmark: AR_1 vs 13 competing forecasts.\n\n")


N_REPS <- 20L
alpha <- 0.05
block_length <- 5L
table5_boot_reps <- 999L

# =============================================================================
# TABLE 4 — Runtime scaling by model count and bootstrap replications
# =============================================================================

# The AR_1 workflow has 13 observed competitor columns. For K = 30, append
# 17 synthetic loss-difference series. Construction is done once before
# timing and is explicitly excluded from runtime measurements.

set.seed(20260822)

max_models <- 30L
observed_model_count <- ncol(loss_diff)

if (observed_model_count != 13L) {
  stop("Expected 13 observed AR_1 competitor loss-difference columns.")
}

loss_diff_scaling <- loss_diff

if (max_models > observed_model_count) {
  extra_model_count <- max_models - observed_model_count
  
  source_columns <- sample(
    seq_len(observed_model_count),
    size = extra_model_count,
    replace = TRUE
  )
  
  additional_loss_diff <- sapply(
    source_columns,
    function(source_column) {
      base_series <- loss_diff[, source_column]
      
      innovation_sd <- 0.05 * sd(base_series, na.rm = TRUE)
      
      innovation <- rnorm(
        P,
        mean = 0,
        sd = innovation_sd
      )
      
      perturbation <- numeric(P)
      perturbation[1] <- innovation[1]
      
      for (time_index in 2:P) {
        perturbation[time_index] <-
          0.30 * perturbation[time_index - 1] +
          innovation[time_index]
      }
      
      base_series + perturbation
    }
  )
  
  if (is.null(dim(additional_loss_diff))) {
    additional_loss_diff <- matrix(
      additional_loss_diff,
      ncol = extra_model_count
    )
  }
  
  colnames(additional_loss_diff) <- paste0(
    "Synthetic_Model_",
    seq_len(extra_model_count)
  )
  
  loss_diff_scaling <- cbind(
    loss_diff,
    additional_loss_diff
  )
}

if (ncol(loss_diff_scaling) != max_models) {
  stop("The scaling matrix does not contain 30 competitor columns.")
}

model_subsets <- list(
  `5` = seq_len(5L),
  `13` = seq_len(13L),
  `30` = seq_len(30L)
)

run_full_battery_once <- function(K, n_boot) {
  selected_columns <- model_subsets[[as.character(K)]]
  
  loss_diff_current <- loss_diff_scaling[
    ,
    selected_columns,
    drop = FALSE
  ]
  
  white_reality_check(
    loss_differences = loss_diff_current,
    n_simulations = n_boot,
    block_length = block_length,
    alpha = alpha
  )
  
  superior_predictive_ability_test(
    loss_differences = loss_diff_current,
    num_bootstrap_replications = n_boot,
    block_length = block_length,
    alpha = alpha
  )
  
  white_reality_check_conditional(
    loss_differences = loss_diff_current,
    weighting_vector = weighting_vector,
    block_length = block_length,
    num_bootstrap_replications = n_boot,
    alpha = alpha
  )
  
  invisible(NULL)
}

model_counts <- c(5L, 13L, 30L)
bootstrap_replications <- c(199L, 999L, 4999L)

table4_results <- expand.grid(
  Models = model_counts,
  Bootstrap_Replications = bootstrap_replications,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

table4_results$Minimum_Runtime_Sec <- NA_real_
table4_results$Median_Runtime_Sec <- NA_real_
table4_results$Mean_Runtime_Sec <- NA_real_
table4_results$Maximum_Runtime_Sec <- NA_real_
table4_results$Runs <- N_REPS

set.seed(20260822)

for (row_index in seq_len(nrow(table4_results))) {
  K <- table4_results$Models[row_index]
  
  n_boot_current <-
    table4_results$Bootstrap_Replications[row_index]
  
  timings <- vapply(
    seq_len(N_REPS),
    function(replication_index) {
      system.time(
        run_full_battery_once(
          K = K,
          n_boot = n_boot_current
        )
      )[["elapsed"]]
    },
    numeric(1)
  )
  
  table4_results$Minimum_Runtime_Sec[row_index] <- min(timings)
  table4_results$Median_Runtime_Sec[row_index] <- median(timings)
  table4_results$Mean_Runtime_Sec[row_index] <- mean(timings)
  table4_results$Maximum_Runtime_Sec[row_index] <- max(timings)
}

table4_results$Benchmark <- benchmark_name
table4_results$Realized_Series <- realized_name

table4_results$Data_Source <- ifelse(
  table4_results$Models <= 13L,
  "Observed AR_1-versus-competitor loss differentials",
  "Observed AR_1 loss differentials plus 17 synthetic AR(1)-perturbed series"
)

table4_results$Protocol_Note <- paste0(
  "Median/minimum/mean/maximum of ",
  N_REPS,
  " full WRC+SPA+CPA battery runs; setup excluded from timing."
)

cat("\n=== TABLE 4 SUMMARY ===\n")
print(table4_results, row.names = FALSE)

write.csv(
  table4_results,
  file.path(output_dir, "runtime_table4_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# TABLE 5 — Runtime by core RCtest function
# =============================================================================

set.seed(20260822)

table5_benchmark <- microbenchmark(
  white_reality_check = white_reality_check(
    loss_differences = loss_diff,
    n_simulations = table5_boot_reps,
    block_length = block_length,
    alpha = alpha
  ),
  
  superior_predictive_ability_test = superior_predictive_ability_test(
    loss_differences = loss_diff,
    num_bootstrap_replications = table5_boot_reps,
    block_length = block_length,
    alpha = alpha
  ),
  
  white_reality_check_conditional = white_reality_check_conditional(
    loss_differences = loss_diff,
    weighting_vector = weighting_vector,
    block_length = block_length,
    num_bootstrap_replications = table5_boot_reps,
    alpha = alpha
  ),
  
  kullback_leibler_test = kullback_leibler_test(
    log_likelihood_differences = log_lik_diff,
    block_length = block_length,
    num_bootstrap_replications = table5_boot_reps,
    alpha = alpha
  ),
  
  reality_check_zp_test = reality_check_zp_test(
    zp_loss_differences = zp_diff,
    block_length = block_length,
    num_bootstrap_replications = table5_boot_reps,
    alpha = alpha
  ),
  
  compute_crps = compute_crps(
    forecast_density = as.numeric(metals[1, forecast_model_cols]),
    target_realization = metals[1, realized_col]
  ),
  
  compute_kupiec = compute_kupiec(
    forecast_matrix = metals,
    forecast_sd_models = forecast_sd_models,
    benchmark_col = benchmark_col,
    alpha = alpha
  ),
  
  times = N_REPS
)

cat("\n=== TABLE 5 RAW BENCHMARK RESULTS ===\n")
print(table5_benchmark)

table5_results <- summary(table5_benchmark, unit = "s")

table5_results$Median_Runtime_Sec <- table5_results$median
table5_results$Minimum_Runtime_Sec <- table5_results$min
table5_results$Maximum_Runtime_Sec <- table5_results$max
table5_results$Mean_Runtime_Sec <- table5_results$mean
table5_results$Runs <- table5_results$neval

table5_results$Bootstrap_Replications <- ifelse(
  table5_results$expr %in% c(
    "white_reality_check",
    "superior_predictive_ability_test",
    "white_reality_check_conditional",
    "kullback_leibler_test",
    "reality_check_zp_test"
  ),
  table5_boot_reps,
  NA_integer_
)

table5_results <- table5_results[, c(
  "expr",
  "Runs",
  "Bootstrap_Replications",
  "Minimum_Runtime_Sec",
  "Median_Runtime_Sec",
  "Mean_Runtime_Sec",
  "Maximum_Runtime_Sec"
)]

names(table5_results)[names(table5_results) == "expr"] <- "Function"

table5_results$Benchmark <- benchmark_name
table5_results$Realized_Series <- realized_name
table5_results$Competing_Models <- 13L

cat("\n=== TABLE 5 SUMMARY ===\n")
print(table5_results, row.names = FALSE)

write.csv(
  table5_results,
  file.path(output_dir, "runtime_table5_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# MEMORY PROFILING
# =============================================================================

cat("\n=== MEMORY PROFILING ===\n")

memory_profile_file <- file.path(
  output_dir,
  "runtime_memory_profile.out"
)

memory_profile_available <- TRUE

tryCatch(
  {
    Rprofmem(memory_profile_file, threshold = 0)
    
    run_full_battery_once(
      K = 13L,
      n_boot = 999L
    )
    
    Rprofmem(NULL)
  },
  error = function(error) {
    cat(
      "Rprofmem() unavailable on this R build: ",
      conditionMessage(error),
      "\n",
      sep = ""
    )
    
    memory_profile_available <<- FALSE
  }
)

if (memory_profile_available && file.exists(memory_profile_file)) {
  memory_lines <- readLines(memory_profile_file)
  
  byte_values <- suppressWarnings(
    as.numeric(
      vapply(
        strsplit(memory_lines, " "),
        function(parts) parts[1],
        character(1)
      )
    )
  )
  
  total_allocated_bytes <- sum(
    byte_values,
    na.rm = TRUE
  )
  
  memory_summary <- data.frame(
    Benchmark = benchmark_name,
    Realized_Series = realized_name,
    Competing_Models = 13L,
    Bootstrap_Replications = 999L,
    Tests = "WRC + SPA + CPA",
    Allocation_Metric = "Total bytes allocated during one complete battery call",
    Total_Bytes = total_allocated_bytes,
    Approx_MB = round(total_allocated_bytes / 1e6, 2),
    stringsAsFactors = FALSE
  )
  
  cat(
    "Total allocated bytes during one full battery call: ",
    total_allocated_bytes,
    "\n",
    sep = ""
  )
  
  cat(
    "Approximate allocated MB: ",
    memory_summary$Approx_MB,
    "\n",
    sep = ""
  )
  
  write.csv(
    memory_summary,
    file.path(output_dir, "runtime_memory_profile_summary.csv"),
    row.names = FALSE
  )
} else {
  message(
    "Memory profiling could not be completed. ",
    "Report this limitation explicitly."
  )
}

cat("\n=== COMPLETED ===\n")
cat(
  "Saved benchmark outputs in:\n",
  normalizePath(output_dir),
  "\n",
  sep = ""
)