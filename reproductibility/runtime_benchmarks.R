# Reproduces Table 4 (runtime by model count / bootstrap replications),
# Table 5 (runtime by core test function), and the memory-usage figure
# reported in the manuscript's Software functionalities / Illustrative
# examples sections from the article on SoftwareX.
#
# This script explicitly documents its own benchmarking protocol, addressing
# every item requested by the reviewer:
#   1. Repeated timing runs: 20 per configuration (N_REPS)
#   2. Compilation/data prep/RNG: data prep and RNG happen ONCE, OUTSIDE the
#      timed block, before microbenchmark() is called -- only the function
#      call itself is timed. R has no separate compilation step (interpreted).
#   3. Mean/median/min: median reported (robust to scheduling noise)
#   4. R version and package versions: printed explicitly below
#   5. Parallel computing: NONE used (single-threaded)
#   6. Random seed: fixed and disclosed (20260822)
#   7. Benchmark code: THIS SCRIPT
#   8. Memory use: ACTUALLY MEASURED via Rprofmem().
#
# Requires: RCtest, microbenchmark 

library(RCtest)
library(microbenchmark)
data(metals)

cat("=== SESSION / VERSION INFO (disclosed per reviewer request) ===\n")
cat("R version:", R.version.string, "\n")
print(installed.packages()[c("RCtest", "ggplot2", "gridExtra", "ggrepel", "rlang"),
                           "Version"])
cat("Parallel computing used: NONE (single-threaded)\n")
cat("Random seed: 20260822\n\n")

# -----------------------------------------------------------------------------
# DATA PREPARATION AND RNG 
# -----------------------------------------------------------------------------
set.seed(20260822)

P <- nrow(metals)
K_total <- ncol(metals)
realized <- metals[, K_total]

bench_loss <- (metals[, K_total] - realized)^2
model_loss <- sweep(metals[, 1:14], 1, realized, FUN = function(x, y) (x - y)^2)
loss_diff <- bench_loss - model_loss

forecast_variance <- estimate_forecast_variance(metals, realized = realized,
                                                benchmark_col = K_total,
                                                window_size = 20)
forecast_sd_models <- sqrt(forecast_variance[, 1:14])
threshold_val <- quantile(realized, 0.05)
weighting_vector <- abs(realized)

klic_loss <- compute_klic(metals, forecast_sd_models, benchmark_col = K_total)
log_lik_diff <- klic_loss[, K_total] - klic_loss[, 1:14]
zp_loss <- compute_zp(metals, forecast_sd_models, threshold = threshold_val,
                      benchmark_col = K_total)
zp_diff <- zp_loss[, K_total] - zp_loss[, 1:14]

cat("Data preparation and RNG draws above are OUTSIDE the timed block.\n")
cat("Only the function calls themselves are timed below.\n\n")

# -----------------------------------------------------------------------------
# TABLE 5 -- TIMED BLOCK. N_REPS repeated timings, median reported.
# RNG used INSIDE bootstrap-based functions (WRC/SPA/CPA/KLIC/ZP resampling)
# IS included in each timed call, since that RNG cost is an intrinsic part
# of running the function, not a setup step -- this is disclosed explicitly.
# -----------------------------------------------------------------------------
N_REPS <- 20

set.seed(20260822)
bm <- microbenchmark(
  white_reality_check = white_reality_check(loss_diff, n_simulations = 999, block_length = 5),
  superior_predictive_ability_test = superior_predictive_ability_test(
    loss_diff, num_bootstrap_replications = 999, block_length = 5, alpha = 0.05),
  white_reality_check_conditional = white_reality_check_conditional(
    loss_diff, weighting_vector = weighting_vector, block_length = 5,
    num_bootstrap_replications = 999, alpha = 0.05),
  kullback_leibler_test = kullback_leibler_test(
    log_lik_diff, block_length = 5, num_bootstrap_replications = 999, alpha = 0.05),
  reality_check_zp_test = reality_check_zp_test(
    zp_diff, block_length = 5, num_bootstrap_replications = 999, alpha = 0.05),
  compute_crps = compute_crps(as.numeric(metals[1, 1:14]), metals[1, K_total]),
  compute_kupiec = compute_kupiec(metals, forecast_sd_models,
                                  benchmark_col = K_total, alpha = 0.05),
  times = N_REPS
)

print(bm)
bm_summary <- summary(bm, unit = "s")
write.csv(bm_summary, "runtime_table5_summary.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# TABLE 4 -- scaling by model count and bootstrap reps, same disclosed protocol
# -----------------------------------------------------------------------------
run_full_battery_once <- function(K, n_boot) {
  cols <- sample(seq_len(14), K)
  ld <- loss_diff[, cols, drop = FALSE]
  white_reality_check(ld, n_simulations = n_boot, block_length = 5)
  superior_predictive_ability_test(ld, num_bootstrap_replications = n_boot,
                                   block_length = 5, alpha = 0.05)
  white_reality_check_conditional(ld, weighting_vector = weighting_vector,
                                  block_length = 5, num_bootstrap_replications = n_boot,
                                  alpha = 0.05)
}

model_counts <- c(5, 14)
boot_reps <- c(199, 999, 4999)
table4_results <- expand.grid(Models = model_counts, Reps = boot_reps)
table4_results$Median_Runtime_Sec <- NA_real_

set.seed(20260822)
for (i in seq_len(nrow(table4_results))) {
  K <- table4_results$Models[i]; n_boot <- table4_results$Reps[i]
  timings <- sapply(seq_len(N_REPS), function(x)
    system.time(run_full_battery_once(K, n_boot))["elapsed"])
  table4_results$Median_Runtime_Sec[i] <- median(timings)
}
print(table4_results)
write.csv(table4_results, "runtime_table4_summary.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# MEMORY PROFILING
# Uses utils::Rprofmem() to log every memory allocation during one full-battery
# call at 14 models / 999 reps, then parses total bytes allocated.
# NOTE: Rprofmem requires R built with --enable-memory-profiling (most CRAN
# binary builds have this enabled by default; if this errors, memory profiling
# is unavailable on this build and should be reported as such rather than
# silently omitted).
# -----------------------------------------------------------------------------
cat("\n=== MEMORY PROFILING (14 models, 999 bootstrap reps) ===\n")
memprof_file <- "runtime_memory_profile.out"
memprof_available <- TRUE
tryCatch({
  Rprofmem(memprof_file, threshold = 0)
  run_full_battery_once(14, 999)
  Rprofmem(NULL)
}, error = function(e) {
  cat("Rprofmem() unavailable on this R build:", conditionMessage(e), "\n")
  memprof_available <<- FALSE
})

if (memprof_available && file.exists(memprof_file)) {
  memprof_lines <- readLines(memprof_file)
  byte_values <- suppressWarnings(as.numeric(sapply(strsplit(memprof_lines, " "), `[`, 1)))
  total_bytes <- sum(byte_values, na.rm = TRUE)
  cat("Total bytes allocated during one full battery call:", total_bytes, "\n")
  cat("Approx. MB allocated:", round(total_bytes / 1e6, 2), "\n")
  write.csv(data.frame(Total_Bytes = total_bytes, Approx_MB = round(total_bytes / 1e6, 2)),
            "runtime_memory_profile_summary.csv", row.names = FALSE)
} else {
  cat("Memory profiling could not be completed on this build.\n")
  cat("Report this limitation explicitly rather than restating the original\n")
  cat("qualitative claim -- see manuscript text guidance.\n")
}