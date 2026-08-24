library(testthat)
library(RCtest)

# =============================================================================
# Tests for invalid inputs and missing values
# =============================================================================

test_that("white_reality_check() rejects non-numeric input", {
  bad_input <- data.frame(a = letters[1:10], b = letters[11:20])
  expect_error(white_reality_check(bad_input, n_simulations = 10, block_length = 5))
})

test_that("white_reality_check() rejects a plain list", {
  bad_input <- list(1, 2, 3)
  expect_error(white_reality_check(bad_input, n_simulations = 10, block_length = 5))
})

test_that("functions handle NA values without crashing silently or returning misleading output", {
  set.seed(1)
  loss_diff_with_na <- matrix(rnorm(100), ncol = 1)
  loss_diff_with_na[5] <- NA
  loss_diff_with_na[50] <- NA

  result <- tryCatch(
    white_reality_check(loss_diff_with_na, n_simulations = 50, block_length = 5),
    error = function(e) e,
    warning = function(w) w
  )

  if (inherits(result, "error") || inherits(result, "warning")) {
    succeed("NA input produces an informative error/warning (acceptable).")
  } else {
    expect_true(is.finite(result$p.value),
                info = "NA input silently produced a non-finite p-value with no warning.")
  }
})

test_that("estimate_forecast_variance() rejects mismatched realized length", {
  forecast_matrix <- matrix(rnorm(200), ncol = 2)
  wrong_length_realized <- rnorm(50)

  expect_error(
    estimate_forecast_variance(forecast_matrix, realized = wrong_length_realized,
                                benchmark_col = 2, window_size = 20)
  )
})

test_that("compute_klic() rejects a forecast_sd_models matrix with the wrong number of columns", {
  forecast_matrix <- matrix(rnorm(300), ncol = 3)
  wrong_sd_matrix <- matrix(abs(rnorm(100)), ncol = 1)

  expect_error(
    compute_klic(forecast_matrix, wrong_sd_matrix, benchmark_col = 3)
  )
})

test_that("white_reality_check() does not silently return a finite result on zero competing models", {
  single_col <- matrix(rnorm(50), ncol = 1)
  zero_col_input <- single_col[, integer(0), drop = FALSE]

  expect_warning(
    result <- white_reality_check(zero_col_input, n_simulations = 10, block_length = 5)
  )
  expect_false(is.finite(unname(result$statistic)),
               info = "Expected a non-finite (-Inf) statistic on zero competing models, given the current warning-based behavior.")
})

test_that("functions reject a block_length that is not a positive integer", {
  loss_diff <- matrix(rnorm(100), ncol = 1)
  expect_error(white_reality_check(loss_diff, n_simulations = 10, block_length = 0))
  expect_error(white_reality_check(loss_diff, n_simulations = 10, block_length = -3))
})
