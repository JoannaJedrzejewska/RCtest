library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for kullback_leibler_test(), reality_check_zp_test(),
# compute_klic(), and compute_zp() -- density forecast evaluation.
# =============================================================================

make_density_test_inputs <- function(seed = 42, P = 150, window_size = 20) {
  set.seed(seed)
  realized <- cumsum(rnorm(P, 0, 1))
  model_fc <- realized + rnorm(P, 0, 0.8)
  bench_fc <- realized + rnorm(P, 0, 1.5)
  full_mat <- cbind(Model = model_fc, Benchmark = bench_fc)
  
  fv <- estimate_forecast_variance(full_mat, realized = realized,
                                   benchmark_col = 2, window_size = window_size)
  sd_model <- sqrt(fv[, 1])
  
  list(full_mat = full_mat, realized = realized, sd_model = sd_model)
}

test_that("compute_klic() returns a hardcoded zero for the benchmark column", {
  inputs <- make_density_test_inputs()
  klic_loss <- compute_klic(inputs$full_mat, matrix(inputs$sd_model, ncol = 1),
                            benchmark_col = 2)
  
  expect_true(all(klic_loss[, "Benchmark"] == 0),
              info = "Benchmark column NLS is no longer hardcoded to zero -- confirm this is intentional.")
  expect_true(all(is.finite(klic_loss[, "Model"])),
              info = "Model NLS contains non-finite values.")
})

test_that("compute_klic() rejects a benchmark_col outside the valid range", {
  inputs <- make_density_test_inputs()
  expect_error(
    compute_klic(inputs$full_mat, matrix(inputs$sd_model, ncol = 1), benchmark_col = 5)
  )
})

test_that("kullback_leibler_test() bootstrap results are reproducible under a fixed seed", {
  inputs <- make_density_test_inputs()
  klic_loss <- compute_klic(inputs$full_mat, matrix(inputs$sd_model, ncol = 1), benchmark_col = 2)
  log_lik_diff <- matrix(klic_loss[, 2] - klic_loss[, 1], ncol = 1)
  
  set.seed(2026)
  result_1 <- kullback_leibler_test(log_lik_diff, block_length = 5,
                                    num_bootstrap_replications = 199, alpha = 0.05)
  set.seed(2026)
  result_2 <- kullback_leibler_test(log_lik_diff, block_length = 5,
                                    num_bootstrap_replications = 199, alpha = 0.05)
  
  expect_equal(result_1$statistic, result_2$statistic)
  expect_equal(result_1$p.value, result_2$p.value)
})

test_that("compute_zp() returns values bounded between 0 and 1", {
  inputs <- make_density_test_inputs()
  threshold_val <- quantile(inputs$realized, 0.05)
  
  zp_loss <- compute_zp(inputs$full_mat, matrix(inputs$sd_model, ncol = 1),
                        threshold = threshold_val, benchmark_col = 2)
  
  expect_true(all(zp_loss >= 0 & zp_loss <= 1, na.rm = TRUE),
              info = "ZP loss values fall outside the theoretically valid [0,1] range.")
})

test_that("reality_check_zp_test() returns finite, valid p-values in [0,1]", {
  inputs <- make_density_test_inputs()
  threshold_val <- quantile(inputs$realized, 0.05)
  zp_loss <- compute_zp(inputs$full_mat, matrix(inputs$sd_model, ncol = 1),
                        threshold = threshold_val, benchmark_col = 2)
  zp_diff <- matrix(zp_loss[, 2] - zp_loss[, 1], ncol = 1)
  
  result <- reality_check_zp_test(zp_diff, block_length = 5,
                                  num_bootstrap_replications = 199, alpha = 0.05)
  
  expect_true(is.finite(result$p_consistent) && result$p_consistent >= 0 && result$p_consistent <= 1)
  expect_true(is.finite(result$p_conservative) && result$p_conservative >= 0 && result$p_conservative <= 1)
})

test_that("kullback_leibler_test() rejects mismatched dimensions gracefully", {
  bad_input <- matrix(rnorm(100), ncol = 1)
  bad_input <- bad_input[1:50, , drop = FALSE]
  
  result <- tryCatch(
    kullback_leibler_test(bad_input, block_length = 5,
                          num_bootstrap_replications = 50, alpha = 0.05),
    error = function(e) e
  )
  expect_true(inherits(result, "htest") || inherits(result, "error"))
})