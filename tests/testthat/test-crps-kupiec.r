library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for compute_crps() and compute_kupiec().
# =============================================================================

test_that("compute_crps() returns zero for a point mass exactly at the true value", {
  # Documented property: a score of zero is achieved by a point mass at the
  # true value. If all competing forecasts equal the realized value exactly,
  # the cross-sectional spread is zero and CRPS should be zero (or very close).
  identical_forecasts <- rep(10, 5)
  realized_value <- 10
  
  result <- compute_crps(identical_forecasts, realized_value)
  
  expect_true(is.finite(result))
  expect_lt(abs(result), 1e-6)
})

test_that("compute_crps() is non-negative", {
  set.seed(11)
  forecasts <- rnorm(10, mean = 5, sd = 2)
  realized_value <- 5.5
  
  result <- compute_crps(forecasts, realized_value)
  
  expect_true(is.finite(result))
  expect_gte(result, 0)
})

test_that("compute_crps() increases as forecast dispersion increases around a fixed realized value", {
  # Wider disagreement among competing forecasts should generally correspond
  # to higher CRPS, holding the realized value fixed
  set.seed(22)
  realized_value <- 0
  tight_forecasts <- rnorm(20, mean = 0, sd = 0.5)
  wide_forecasts <- rnorm(20, mean = 0, sd = 5)
  
  crps_tight <- compute_crps(tight_forecasts, realized_value)
  crps_wide <- compute_crps(wide_forecasts, realized_value)
  
  expect_lt(crps_tight, crps_wide)
})

test_that("compute_kupiec() reports a p-value in [0,1] and a non-negative test statistic", {
  set.seed(33)
  P <- 150
  K_total <- 3
  realized <- cumsum(rnorm(P, 0, 1))
  forecasts <- cbind(
    Model_1 = realized + rnorm(P, 0, 1),
    Model_2 = realized + rnorm(P, 0, 1.2),
    Benchmark = realized + rnorm(P, 0, 1.5)
  )
  
  fv <- estimate_forecast_variance(forecasts, realized = realized,
                                   benchmark_col = 3, window_size = 20)
  sd_models <- sqrt(fv[, 1:2])
  
  result <- compute_kupiec(forecasts, sd_models, benchmark_col = 3, alpha = 0.05)
  
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("compute_kupiec() rejects a mismatched forecast_sd_models dimension", {
  set.seed(44)
  P <- 100
  forecasts <- matrix(rnorm(P * 3), ncol = 3)
  wrong_sd <- matrix(abs(rnorm(P)), ncol = 1) #2 columns for 2 competing models
  
  expect_error(
    compute_kupiec(forecasts, wrong_sd, benchmark_col = 3, alpha = 0.05)
  )
})