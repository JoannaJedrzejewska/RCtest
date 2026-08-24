library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for white_reality_check_conditional() (CPA test).
# =============================================================================

make_loss_differences <- function(seed = 123, n = 200L) {
  set.seed(seed)
  benchmark_loss <- rnorm(n, mean = 4, sd = 1)^2
  cbind(
    Better = benchmark_loss - rnorm(n, mean = 2, sd = 1)^2,
    Worse  = benchmark_loss - rnorm(n, mean = 6, sd = 1)^2
  )
}

test_that("CPA with a constant weighting vector reduces to a WRC-like result", {
  loss_differences <- make_loss_differences()
  n <- nrow(loss_differences)
  constant_weights <- rep(1, n)

  result <- white_reality_check_conditional(
    loss_differences, weighting_vector = constant_weights,
    block_length = 5, num_bootstrap_replications = 199, alpha = 0.05
  )

  expect_s3_class(result, "htest")
  expect_true(is.finite(unname(result$statistic)))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("CPA rejects a weighting vector of the wrong length", {
  loss_differences <- make_loss_differences()
  wrong_length_weights <- rep(1, 50)  # loss_differences has 200 rows

  expect_error(
    white_reality_check_conditional(
      loss_differences, weighting_vector = wrong_length_weights,
      block_length = 5, num_bootstrap_replications = 99, alpha = 0.05
    )
  )
})

test_that("CPA bootstrap results are reproducible under a fixed seed", {
  loss_differences <- make_loss_differences()
  weights <- abs(rnorm(nrow(loss_differences)))

  set.seed(2026)
  result_1 <- white_reality_check_conditional(
    loss_differences, weighting_vector = weights,
    block_length = 5, num_bootstrap_replications = 199, alpha = 0.05
  )
  set.seed(2026)
  result_2 <- white_reality_check_conditional(
    loss_differences, weighting_vector = weights,
    block_length = 5, num_bootstrap_replications = 199, alpha = 0.05
  )

  expect_equal(result_1$statistic, result_2$statistic)
  expect_equal(result_1$p.value, result_2$p.value)
})

test_that("CPA with an absolute-realization instrument produces finite output", {
  set.seed(99)
  loss_differences <- make_loss_differences(seed = 99)
  realized_proxy <- abs(rnorm(nrow(loss_differences), mean = 5, sd = 2))

  result <- white_reality_check_conditional(
    loss_differences, weighting_vector = realized_proxy,
    block_length = 5, num_bootstrap_replications = 99, alpha = 0.05
  )

  expect_true(is.finite(unname(result$statistic)))
  expect_true(is.finite(result$p.value))
})
