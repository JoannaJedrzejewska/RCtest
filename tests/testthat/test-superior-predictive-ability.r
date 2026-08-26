library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for superior_predictive_ability_test(), following the same
# pattern established for white_reality_check().
# =============================================================================

make_loss_differences <- function(seed = 123, n = 200L) {
  set.seed(seed)
  benchmark_loss <- rnorm(n, mean = 4, sd = 1)^2
  cbind(
    Better = benchmark_loss - rnorm(n, mean = 2, sd = 1)^2,
    Worse  = benchmark_loss - rnorm(n, mean = 6, sd = 1)^2
  )
}

test_that("SPA loss-differential sign convention is respected", {
  loss_differences <- make_loss_differences()

  result <- superior_predictive_ability_test(
    loss_differences, num_bootstrap_replications = 199, block_length = 5, alpha = 0.05
  )

  expect_gt(mean(loss_differences[, "Better"]), 0)
  expect_lt(mean(loss_differences[, "Worse"]), 0)
  expect_true(is.finite(unname(result$statistic)))
  expect_true(result$p_consistent >= 0 && result$p_consistent <= 1)
  expect_true(result$p_conservative >= 0 && result$p_conservative <= 1)
})

test_that("SPA consistent p-value is never greater than the conservative p-value", {
  # Documented property from the manuscript: the consistent p-value is
  # always no greater than the conservative counterpart.
  loss_differences <- make_loss_differences(seed = 321)

  result <- superior_predictive_ability_test(
    loss_differences, num_bootstrap_replications = 199, block_length = 5, alpha = 0.05
  )

  expect_lte(result$p_consistent, result$p_conservative)
})

test_that("SPA bootstrap results are reproducible under a fixed seed", {
  loss_differences <- make_loss_differences()

  set.seed(2026)
  result_1 <- superior_predictive_ability_test(
    loss_differences, num_bootstrap_replications = 199, block_length = 5, alpha = 0.05
  )
  set.seed(2026)
  result_2 <- superior_predictive_ability_test(
    loss_differences, num_bootstrap_replications = 199, block_length = 5, alpha = 0.05
  )

  expect_equal(result_1$statistic, result_2$statistic)
  expect_equal(result_1$p_consistent, result_2$p_consistent)
  expect_equal(result_1$p_conservative, result_2$p_conservative)
})

test_that("SPA handles constant loss differentials without producing non-finite output", {
  loss_differences <- cbind(
    Equal = rep(0, 100),
    Inferior = rep(-0.25, 100)
  )

  result <- superior_predictive_ability_test(
    loss_differences, num_bootstrap_replications = 99, block_length = 5, alpha = 0.05
  )

  expect_s3_class(result, "htest")
  expect_true(is.finite(unname(result$statistic)))
  expect_true(is.finite(result$p_consistent))
})

test_that("SPA handles a single competing model", {
  set.seed(7)
  loss_differences <- matrix(
    rnorm(100, mean = 0.1, sd = 1), ncol = 1, dimnames = list(NULL, "Model_1")
  )

  result <- superior_predictive_ability_test(
    loss_differences, num_bootstrap_replications = 99, block_length = 5, alpha = 0.05
  )

  expect_s3_class(result, "htest")
  expect_length(result$statistic, 1)
  expect_true(is.finite(unname(result$statistic)))
})
