library(testthat)
library(RCtest)

# =============================================================================
# Unit tests for white_reality_check_cdf_approx() and
# run_comprehensive_erc_analysis().
# =============================================================================

make_loss_differences <- function(seed = 123, n = 200L) {
  set.seed(seed)
  
  benchmark_loss <- rnorm(n, mean = 4, sd = 1)^2
  
  cbind(
    Better = benchmark_loss - rnorm(n, mean = 2, sd = 1)^2,
    Worse  = benchmark_loss - rnorm(n, mean = 6, sd = 1)^2
  )
}

test_that(
  "white_reality_check_cdf_approx() returns finite statistic and valid p-value",
  {
    loss_differences <- make_loss_differences()
    
    result <- white_reality_check_cdf_approx(
      loss_differences,
      block_length = 5,
      num_bootstrap_replications = 199,
      alpha = 0.05
    )
    
    expect_s3_class(result, "htest")
    expect_true(is.finite(unname(result$statistic)))
    expect_true(is.finite(result$p.value))
    expect_true(result$p.value >= 0 && result$p.value <= 1)
  }
)

test_that(
  "white_reality_check_cdf_approx() is reproducible under a fixed seed",
  {
    loss_differences <- make_loss_differences()
    
    set.seed(2026)
    result_1 <- white_reality_check_cdf_approx(
      loss_differences,
      block_length = 5,
      num_bootstrap_replications = 199,
      alpha = 0.05
    )
    
    set.seed(2026)
    result_2 <- white_reality_check_cdf_approx(
      loss_differences,
      block_length = 5,
      num_bootstrap_replications = 199,
      alpha = 0.05
    )
    
    expect_equal(result_1$statistic, result_2$statistic)
    expect_equal(result_1$p.value, result_2$p.value)
  }
)

test_that(
  "white_reality_check_cdf_approx() handles a single competing model",
  {
    set.seed(7)
    
    loss_differences <- matrix(
      rnorm(100, mean = 0.1, sd = 1),
      ncol = 1,
      dimnames = list(NULL, "Model_1")
    )
    
    result <- white_reality_check_cdf_approx(
      loss_differences,
      block_length = 5,
      num_bootstrap_replications = 99,
      alpha = 0.05
    )
    
    expect_s3_class(result, "htest")
    expect_true(is.finite(unname(result$statistic)))
  }
)

test_that(
  "run_comprehensive_erc_analysis() completes end-to-end on a minimal valid dataset",
  {
    set.seed(101)
    
    P <- 60L
    dataset_name <- "synthetic"
    R_start <- 0L
    
    realized_values <- cumsum(rnorm(P, mean = 0, sd = 1))
    
    forecasts <- cbind(
      Model_A = realized_values + rnorm(P, mean = 0, sd = 1.0),
      Model_B = realized_values + rnorm(P, mean = 0, sd = 1.2),
      Benchmark = realized_values + rnorm(P, mean = 0, sd = 1.5)
    )
    
    data_list_prepared <- list(
      synthetic = list(R_start = R_start)
    )
    
    y_hat_all <- list(
      synthetic = list(NULL, NULL, forecasts)
    )
    
    y_raw <- list(
      synthetic = realized_values
    )
    
    # Minimal model-specification matrix aligned to forecasts.
    # Each row describes the corresponding model in `forecasts`.
    mods_matrix <- matrix(
      c(
        1, 0,
        0, 1,
        1, 1
      ),
      nrow = 3,
      ncol = 2,
      byrow = TRUE,
      dimnames = list(
        colnames(forecasts),
        c("Feature_1", "Feature_2")
      )
    )
    
    # Small values keep this an integration/smoke test rather than a
    # computational manuscript-results run.
    set.seed(2026)
    
    result <- run_comprehensive_erc_analysis(
      data_list_prepared = data_list_prepared,
      mods_matrix = mods_matrix,
      alpha_grid = 0.05,
      window_size = 10,
      y_hat_all = y_hat_all,
      y_raw = y_raw,
      block_length = 5,
      n_boot = 19,
      zp_quantile = 0.05,
      n_crps_samples = 5,
      benchmark_col = 3
    )
    
    expect_true(is.list(result))
    
    expect_named(
      result,
      c("aggregate_results", "per_model_results")
    )
    
    expect_true(dataset_name %in% names(result$aggregate_results))
    expect_true(dataset_name %in% names(result$per_model_results))
    
    test_results <- result$aggregate_results[[dataset_name]]
    
    expect_true(
      any(grepl("White_Reality_Check", names(test_results)))
    )
    
    expect_true(
      any(grepl("Superior_Predictive_Ability", names(test_results)))
    )
    
    expect_true(
      any(grepl("Conditional_Predictive_Ability", names(test_results)))
    )
    
    expect_true("ZP" %in% names(test_results))
    expect_true("Kullback_Leibler" %in% names(test_results))
    expect_true("VaR_Backtests" %in% names(test_results))
  }
)