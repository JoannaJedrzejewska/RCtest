library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for compute_crps() and compute_kupiec().
# =============================================================================

test_that("compute_crps() returns zero for a point mass exactly at the true value", {
  identical_forecasts <- rep(10, 5)
  realized_value <- 10
  
  result <- compute_crps(
    forecast_density = identical_forecasts,
    target_realization = realized_value
  )
  
  expect_true(is.finite(result))
  expect_lt(abs(result), 1e-6)
})

test_that("compute_crps() is non-negative", {
  set.seed(11)
  
  forecasts <- rnorm(10, mean = 5, sd = 2)
  realized_value <- 5.5
  
  result <- compute_crps(
    forecast_density = forecasts,
    target_realization = realized_value
  )
  
  expect_true(is.finite(result))
  expect_gte(result, 0)
})

test_that("compute_crps() increases as forecast dispersion increases around a fixed realized value", {
  set.seed(22)
  
  realized_value <- 0
  
  tight_forecasts <- rnorm(
    20,
    mean = 0,
    sd = 0.5
  )
  
  wide_forecasts <- rnorm(
    20,
    mean = 0,
    sd = 5
  )
  
  crps_tight <- compute_crps(
    forecast_density = tight_forecasts,
    target_realization = realized_value
  )
  
  crps_wide <- compute_crps(
    forecast_density = wide_forecasts,
    target_realization = realized_value
  )
  
  expect_lt(crps_tight, crps_wide)
})

test_that("compute_kupiec() returns valid htest results for all competing models", {
  set.seed(33)
  
  P <- 150L
  
  realized <- cumsum(
    rnorm(P, mean = 0, sd = 1)
  )
  
  forecasts <- cbind(
    Model_1 = realized + rnorm(P, mean = 0, sd = 1.0),
    Model_2 = realized + rnorm(P, mean = 0, sd = 1.2),
    Benchmark = realized + rnorm(P, mean = 0, sd = 1.5)
  )
  
  forecast_variance <- estimate_forecast_variance(
    forecast_matrix = forecasts,
    realized = realized,
    benchmark_col = 3,
    window_size = 20
  )
  
  sd_models <- sqrt(
    forecast_variance[, 1:2, drop = FALSE]
  )
  
  result <- compute_kupiec(
    forecast_matrix = forecasts,
    forecast_sd_models = sd_models,
    realized = realized,
    benchmark_col = 3,
    alpha = 0.05
  )
  
  expect_type(result, "list")
  expect_named(result, c("Model_1", "Model_2"))
  expect_length(result, 2)
  
  for (model_name in names(result)) {
    model_result <- result[[model_name]]
    
    expect_s3_class(model_result, "htest")
    
    expect_true(
      is.finite(unname(model_result$statistic))
    )
    
    expect_gte(
      unname(model_result$statistic),
      0
    )
    
    expect_true(
      is.finite(model_result$p.value)
    )
    
    expect_gte(model_result$p.value, 0)
    expect_lte(model_result$p.value, 1)
    
    expect_true(
      is.numeric(model_result$actual_exceedances)
    )
    
    expect_gte(model_result$actual_exceedances, 0)
    expect_lte(model_result$actual_exceedances, P)
    
    expect_equal(
      model_result$expected,
      P * 0.05
    )
    
    expect_equal(
      model_result$null.value,
      c("expected violation rate" = 0.05)
    )
  }
})

test_that("compute_kupiec() rejects a mismatched forecast_sd_models dimension", {
  set.seed(44)
  
  P <- 100L
  
  realized <- cumsum(
    rnorm(P, mean = 0, sd = 1)
  )
  
  forecasts <- cbind(
    Model_1 = realized + rnorm(P, mean = 0, sd = 1.0),
    Model_2 = realized + rnorm(P, mean = 0, sd = 1.2),
    Benchmark = realized + rnorm(P, mean = 0, sd = 1.5)
  )
  
  wrong_sd <- matrix(
    abs(rnorm(P)),
    ncol = 1
  )
  
  expect_error(
    compute_kupiec(
      forecast_matrix = forecasts,
      forecast_sd_models = wrong_sd,
      realized = realized,
      benchmark_col = 3,
      alpha = 0.05
    ),
    "forecast_sd_models"
  )
})