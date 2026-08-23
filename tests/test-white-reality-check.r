library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for white_reality_check().
# Covers sign convention, statistic direction, bootstrap reproducibility,
# constant/degenerate loss differentials, and a single competing model.
# =============================================================================

make_loss_differences <- function(seed = 123, n = 200L) {
  set.seed(seed)

  benchmark_loss <- rnorm(n, mean = 4, sd = 1)^2

  cbind(
    Better = benchmark_loss - rnorm(n, mean = 2, sd = 1)^2,
    Worse  = benchmark_loss - rnorm(n, mean = 6, sd = 1)^2
  )
}

test_that("loss-differential sign convention is respected", {
  loss_differences <- make_loss_differences()

  result <- white_reality_check(
    loss_differences,
    n_simulations = 199,
    block_length = 5
  )

  expect_gt(mean(loss_differences[, "Better"]), 0)
  expect_lt(mean(loss_differences[, "Worse"]), 0)
  expect_true(is.finite(unname(result$statistic)))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("test statistic equals the maximum sample mean loss differential", {
  loss_differences <- make_loss_differences()

  result <- white_reality_check(
    loss_differences,
    n_simulations = 99,
    block_length = 5
  )

  expect_equal(
    unname(result$statistic),
    max(colMeans(loss_differences)),
    tolerance = 1e-12
  )
})

test_that("reversing the loss-difference sign reverses the statistic as expected", {
  loss_differences <- make_loss_differences()

  forward <- white_reality_check(
    loss_differences,
    n_simulations = 99,
    block_length = 5
  )

  reverse <- white_reality_check(
    -loss_differences,
    n_simulations = 99,
    block_length = 5
  )

  expect_equal(
    unname(forward$statistic),
    max(colMeans(loss_differences)),
    tolerance = 1e-12
  )

  expect_equal(
    unname(reverse$statistic),
    -min(colMeans(loss_differences)),
    tolerance = 1e-12
  )
})

test_that("bootstrap results are reproducible under a fixed seed", {
  loss_differences <- make_loss_differences()

  set.seed(2026)
  result_1 <- white_reality_check(
    loss_differences,
    n_simulations = 199,
    block_length = 5
  )

  set.seed(2026)
  result_2 <- white_reality_check(
    loss_differences,
    n_simulations = 199,
    block_length = 5
  )

  expect_equal(result_1$statistic, result_2$statistic)
  expect_equal(result_1$p.value, result_2$p.value)
  expect_equal(result_1$reject_null, result_2$reject_null)
})

test_that("constant loss differentials produce valid htest output", {
  loss_differences <- cbind(
    Equal = rep(0, 100),
    Inferior = rep(-0.25, 100)
  )

  result <- white_reality_check(
    loss_differences,
    n_simulations = 99,
    block_length = 5
  )

  expect_s3_class(result, "htest")
  expect_true(is.finite(unname(result$statistic)))
  expect_true(is.finite(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("one competing model is accepted", {
  set.seed(7)

  loss_differences <- matrix(
    rnorm(100, mean = 0.1, sd = 1),
    ncol = 1,
    dimnames = list(NULL, "Model_1")
  )

  result <- white_reality_check(
    loss_differences,
    n_simulations = 99,
    block_length = 5
  )

  expect_s3_class(result, "htest")
  expect_length(result$statistic, 1)
  expect_true(is.finite(unname(result$statistic)))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})