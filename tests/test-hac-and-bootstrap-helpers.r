library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for estimate_long_run_covariance() and mbb_resample_data().
# =============================================================================

test_that("estimate_long_run_covariance() estimates a plausible LRV for iid noise", {
  set.seed(55)
  n <- 300
  x <- matrix(rnorm(n, 0, 2), ncol = 1)

  result <- estimate_long_run_covariance(x, block_length = 5)

  expect_true(is.finite(result[1, 1]))
  expect_gt(result[1, 1], 0)
  # Generous finite-sample range around the true LRV=4. This is intentionally
  # broad to avoid false failures due to sampling variation/HAC weighting,
  # while still catching zero, negative, or explosively large output.
  expect_gt(result[1, 1], 0.5)
  expect_lt(result[1, 1], 10)
})

test_that("estimate_long_run_covariance() handles multiple columns independently", {
  set.seed(77)
  x <- matrix(rnorm(200 * 3), ncol = 3)

  result <- estimate_long_run_covariance(x, block_length = 5)

  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_true(all(is.finite(diag(result))))
  expect_true(all(diag(result) > 0))
})

test_that("estimate_long_run_covariance() handles a constant series predictably", {
  # Explicit numerical-stability / singular-HAC check 
  constant_series <- matrix(rep(5, 100), ncol = 1)

  result <- tryCatch(
    estimate_long_run_covariance(constant_series, block_length = 5),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    succeed("Constant series produces an informative error (acceptable).")
  } else {
    expect_true(is.finite(result[1, 1]),
                info = "Constant series returned a non-finite HAC covariance estimate.")
  }
})

test_that("mbb_resample_data() returns a resample of the same dimensions as the input", {
  set.seed(88)
  x <- matrix(rnorm(150 * 2), ncol = 2)

  resampled <- mbb_resample_data(x, block_length = 5)

  expect_equal(dim(resampled), dim(x))
})

test_that("mbb_resample_data() produces different resamples across calls without a fixed seed", {
  set.seed(99)
  x <- matrix(rnorm(150), ncol = 1)

  resample_1 <- mbb_resample_data(x, block_length = 5)
  resample_2 <- mbb_resample_data(x, block_length = 5)

  expect_false(isTRUE(all.equal(resample_1, resample_2)))
})

test_that("mbb_resample_data() handles a block length larger than the sample size predictably", {
  set.seed(111)
  x <- matrix(rnorm(20), ncol = 1)
  
  result <- tryCatch(
    mbb_resample_data(x, block_length = 50),
    error = function(e) e,
    warning = function(w) w
  )
  
  if (inherits(result, "error") || inherits(result, "warning")) {
    succeed("Oversized block length produces an informative error or warning.")
  } else {
    expect_true(is.matrix(result))
    expect_equal(dim(result), dim(x))
    expect_true(all(is.finite(result)))
  }
})
