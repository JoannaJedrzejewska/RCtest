#' @importFrom stats rnorm quantile
#' @importFrom utils globalVariables write.csv write.table
utils::globalVariables(c("Test", "Dataset", "Reject_H0", "."))

#' @title Forecasts of Base Metals Prices
#'
#' @description A dataset of 165 observations of base metals price indices. It includes
#' point forecasts from 15 predictive models (Bayesian and shrinkage-based). This dataset
#' is based on the research methodology and benchmarks presented in Drachal and
#' Jedrzejewska (2025).
#'
#' @format A matrix with 165 rows and 15 columns:
#' \itemize{
#'   \item \code{metals[,1:14]}: Point forecasts from 14 competing models (including
#'     variations of Bayesian Dynamic Mixture Models (BDMM), Dynamic Model Averaging (DMA),
#'     Time-Varying Parameters (TVP), as well as Bayesian LASSO, Bayesian RIDGE,
#'     Autoregressive (AR1), and Linear Regression).
#'   \item \code{metals[,15]}: \code{HA} - Historical Average (Benchmark model).
#' }
#'
#' @details
#' Monthly World Bank base metals price index (2010 = 100, USD), including aluminium,
#' copper, lead, nickel, tin and zinc, forecasts for the period from 03/2011 to 11/2024.
#'
#' @source
#' Drachal, K. (2025). Base metals price index forecasts. \emph{figshare}.
#' \doi{10.6084/m9.figshare.28382480.v1}
#'
#' @references
#' Drachal, K., & Jedrzejewska, J. (2025). Forecasting Base Metals Prices: A Comparison
#' of Various Bayesian-Based Methods. \emph{Sinteza 2025}, 175--183.
#' \doi{10.15308/Sinteza-2025-175-183}
#'
#' World Bank. (2025). Commodity markets.
#' \url{https://www.worldbank.org/en/research/commodity-markets}
#'
#' @examples
#' data(metals)
#' head(metals)
#' plot(metals[, 1], type = "l", main = "Competing Model vs Benchmark")
#' lines(metals[, ncol(metals)], col = "red", lty = 2)
"metals"

#' @title Run Comprehensive Forecast Evaluation Analysis
#'
#' @description Runs a full suite of reality check and density forecast evaluation tests
#' on one or more datasets. Tests included: White's Reality Check (WRC), Superior
#' Predictive Ability (SPA), Conditional Predictive Ability (CPA), ZP Quantile Loss
#' test, Kullback-Leibler Information Criterion (KLIC) test, Kupiec Unconditional
#' Coverage (UC) test, CRPS-based CDF comparison, and per-model Diebold-Mariano
#' statistics across four error metrics (MSE, MAE, MASE).
#'
#' \strong{Technical Abbreviations:}
#' \itemize{
#'   \item \strong{WRC:} White's Reality Check (White, 2000). Tests whether any competing
#'     forecast has lower expected loss than the benchmark; controls family-wise error rate.
#'   \item \strong{SPA:} Superior Predictive Ability test (Hansen, 2005). A studentized
#'     extension of WRC with improved power that corrects for irrelevant forecasts.
#'   \item \strong{CPA:} Conditional Predictive Ability test (Giacomini & White, 2006).
#'     Tests whether loss differentials are predictable by a conditioning variable.
#'   \item \strong{ZP:} Quantile Loss test (Corradi & Swanson, 2006). Evaluates whether
#'     any competing forecast better calibrates the probability of a left-tail event
#'     defined by the \code{zp_quantile} threshold.
#'   \item \strong{KLIC:} Kullback-Leibler Information Criterion based density test
#'     (Corradi & Swanson, 2006). Selects the forecast whose predictive density is closest
#'     to the true density in terms of KLIC distance, evaluated via Negative
#'     Log-Likelihood Scores (NLS) under a Gaussian predictive density assumption.
#'   \item \strong{CRPS:} Continuous Ranked Probability Score (Gneiting & Raftery, 2007).
#'     Jointly rewards calibration and sharpness of the predictive distribution.
#'   \item \strong{UC:} Kupiec Unconditional Coverage test (Kupiec, 1995).
#'   \item \strong{MSE:} Mean Squared Error.
#'   \item \strong{MAE:} Mean Absolute Error.
#'   \item \strong{MASE:} Mean Absolute Scaled Error.
#' }
#'
#' @param data_list_prepared Named list of prepared data frames, one per dataset. Each
#'   element must be a list containing at least a field \code{R_start}: a non-negative
#'   integer specifying how many observations to skip from the start of \code{y_raw}
#'   before aligning with the forecast matrix. Set \code{R_start = 0} to use all
#'   available observations. This is used as a warm-up offset when the raw series is
#'   longer than the forecast evaluation window.
#' @param mods_matrix A legacy placeholder matrix retained for interface compatibility
#'   with external pipelines in which forecasts were previously defined as a parameter
#'   matrix. Its contents are not read or used anywhere in this function — the actual
#'   forecast structure is derived entirely from the forecast matrices supplied in
#'   \code{y_hat_all}. Pass \code{matrix(0)} when calling the function directly.
#' @param alpha_grid Numeric scalar or vector of significance levels. Only the first
#'   element (\code{alpha_grid[1]}) is used as the significance level for all tests.
#' @param window_size Integer window size for rolling variance estimation passed to
#'   \code{\link{estimate_forecast_variance}}.
#' @param y_hat_all Named list of forecast results, one per dataset. Each element must
#'   be a list of length at least 3, where the third element (\code{[[3]]}) is a numeric
#'   matrix of dimension \code{P x K_total}: columns \code{1:(K_total-1)} are the
#'   competing model forecasts and column \code{K_total} (or the column indicated by
#'   \code{benchmark_col}) is the benchmark forecast.
#' @param y_raw Named list of raw realized value vectors, one per dataset. Each element
#'   is a numeric vector whose length must be at least \code{R_start + P}.
#' @param block_length Integer block length for Moving Block Bootstrap (MBB) used in
#'   WRC, SPA, CPA, ZP, and KLIC tests. Default is 5. A commonly used rule of thumb is
#'   \eqn{T^{1/3}} (Politis & Romano, 1994); for \code{P = 165} this gives approximately
#'   5--6.
#' @param n_boot \code{\link[base]{integer}} number of MBB bootstrap replications.
#'   Default \code{999}; see Davidson & MacKinnon (2000).
#' @param zp_quantile Numeric quantile level used to define the left-tail threshold
#'   \eqn{\tau} for the ZP test, computed as
#'   \code{quantile(realizations, zp_quantile)}. Default is 0.05 (5th percentile).
#' @param n_crps_samples Integer number of Monte Carlo samples drawn from the
#'   Gaussian predictive distribution \eqn{N(\hat{y}_{k,t}, \hat{\sigma}_{k,t}^2)}
#'   to approximate the Continuous Ranked Probability Score (CRPS) for each forecast and 
#'   time period.Default is 10 (fast, suitable for examples only). For reliable results use
#'   at least 500; for publication-quality estimates use 1000 or more. Higher
#'   values reduce Monte Carlo variance but increase computation time linearly.
#' @param benchmark_col Index or name of the benchmark column in the forecast matrix.
#'   Defaults to the last column (\code{NULL}).
#'
#' @return \code{\link[base]{list}} containing:
#'  \itemize{
#'   \item \code{aggregate_results}: Named list of test results per dataset. Each
#'     dataset element contains named \code{htest} objects for each test and metric
#'     combination, plus \code{VaR_Backtests} (a list of per-model Kupiec \code{htest}
#'     objects).
#'   \item \code{per_model_results}: Named list of per-model Diebold-Mariano statistics
#'     per dataset, returned as \code{data.frame} objects from
#'     \code{\link{compute_per_model_statistics}}.
#' }
#'
#' @seealso
#' \code{\link{create_unified_summary}},
#' \code{\link{generate_comprehensive_report}},
#' \code{\link{white_reality_check}},
#' \code{\link{superior_predictive_ability_test}},
#' \code{\link{white_reality_check_conditional}},
#' \code{\link{reality_check_zp_test}},
#' \code{\link{kullback_leibler_test}},
#' \code{\link{compute_kupiec}}
#'
#' @references
#' Davidson, R., & MacKinnon, J. G. (2000).
#' Bootstrap tests: How many bootstraps?
#' \emph{Econometric Reviews}, 19(1), 55--68.
#' \doi{10.1080/07474930008800459}
#' 
#' White, H. (2000). A reality check for data snooping. \emph{Econometrica}, 68(5),
#' 1097--1126. \doi{10.1111/1468-0262.00152}
#'
#' Hansen, P. R. (2005). A Test for Superior Predictive Ability.
#' \emph{Journal of Business & Economic Statistics}, 23(4), 365--380.
#' \doi{10.1198/073500105000000063}
#'
#' Giacomini, R., & White, H. (2006). Tests of Conditional Predictive Ability.
#' \emph{Econometrica}, 74(6), 1545--1578. \doi{10.1111/j.1468-0262.2006.00718.x}
#'
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional confidence
#' interval accuracy tests. \emph{Journal of Econometrics}, 135(1--2), 187--228.
#' \doi{10.1016/j.jeconom.2005.07.026}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its recent
#' extensions. In \emph{Festschrift in honor of Halbert L. White}.
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules, Prediction,
#' and Estimation. \emph{Journal of the American Statistical Association}, 102(477),
#' 359--378. \doi{10.1198/016214506000001437}
#'
#' Kupiec, P. H. (1995). Techniques for Verifying the Accuracy of Risk Measurement
#' Models. \emph{The Journal of Derivatives}, 3(2), 173--184.
#' \doi{10.3905/jod.1995.407942}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' Diebold, F. X., & Mariano, R. S. (1995). Comparing Predictive Accuracy.
#' \emph{Journal of Business & Economic Statistics}, 13(3), 253--263.
#' \doi{10.1080/07350015.1995.10524599}
#'
#' @examples
#' \donttest{
#' data(metals)
#' ds_name           <- "Dataset1"
#' realizations_list <- list(Dataset1 = metals[, ncol(metals)])
#' prep_list         <- list(Dataset1 = list(R_start = 0))
#' forecasts_list    <- list(list(NULL, NULL, metals))
#' names(forecasts_list) <- ds_name
#'
#' res <- run_comprehensive_erc_analysis(
#'   data_list_prepared = prep_list,
#'   mods_matrix        = matrix(0),
#'   alpha_grid         = c(0.05),
#'   window_size        = 20,
#'   y_hat_all          = forecasts_list,
#'   y_raw              = realizations_list,
#'   n_boot             = 99,
#'   n_crps_samples     = 5
#' )
#' print(res$aggregate_results$Dataset1$White_Reality_Check_MSE)
#'}
#' @export
run_comprehensive_erc_analysis <- function(data_list_prepared, mods_matrix, alpha_grid,
                                           window_size, y_hat_all, y_raw,
                                           block_length   = 5,
                                           n_boot         = 999,
                                           zp_quantile    = 0.05,
                                           n_crps_samples = 10,
                                           benchmark_col  = NULL) {
  all_results           <- list()
  all_per_model_results <- list()
  alpha                 <- alpha_grid[1]
  
  for (dataset_name in names(data_list_prepared)) {
    results_for_dataset <- list()
    per_model_dataset   <- list()
    
    current_data    <- data_list_prepared[[dataset_name]]
    R_start_current <- current_data$R_start
    forecasts_list  <- y_hat_all[[dataset_name]]
    
    if (is.null(forecasts_list) || length(forecasts_list) < 3) {
      message(paste("No forecasts for:", dataset_name))
      next
    }
    
    forecasts <- forecasts_list[[3]]
    K_total   <- ncol(forecasts)
    
    if (is.null(benchmark_col)) {
      bench_idx <- K_total
    } else {
      if (is.character(benchmark_col)) {
        bench_idx <- which(colnames(forecasts) == benchmark_col)
      } else {
        bench_idx <- benchmark_col
      }
    }
    
    comp_cols <- setdiff(seq_len(K_total), bench_idx)
    K_models  <- length(comp_cols)
    
    realizations_raw <- as.numeric(
      y_raw[[dataset_name]]
    )[(R_start_current + 1):length(y_raw[[dataset_name]])]
    
    if (is.null(forecasts) || nrow(forecasts) != length(realizations_raw)) next
    
    forecast_var <- estimate_forecast_variance(forecasts,
                                               realized      = realizations_raw,
                                               benchmark_col = bench_idx,
                                               window_size   = window_size)
    forecast_sd_models <- sqrt(forecast_var[, comp_cols, drop = FALSE])
    forecast_sd_benchmark <- sqrt(forecast_var[, bench_idx])
    benchmark_series   <- forecasts[, bench_idx]
    
    metrics <- c("MSE", "MAE", "MASE")
    for (metric in metrics) {
      loss_matrix <- matrix(NA, nrow(forecasts), K_models)
      
      for (i in seq_len(K_models)) {
        k   <- comp_cols[i]
        err <- realizations_raw - forecasts[, k]
        
        if      (metric == "MSE")  loss_matrix[, i] <- err^2
        else if (metric == "MAE")  loss_matrix[, i] <- abs(err)
        else if (metric == "MASE") loss_matrix[, i] <- abs(err) / mean(abs(diff(realizations_raw)), na.rm = TRUE)
      }
      
      bench_err <- realizations_raw - benchmark_series
      if      (metric == "MSE")  bench_loss <- bench_err^2
      else if (metric == "MAE")  bench_loss <- abs(bench_err)
      else if (metric == "MASE") bench_loss <- abs(bench_err) / mean(abs(diff(realizations_raw)), na.rm = TRUE)
      
      loss_diffs <- bench_loss - loss_matrix
      colnames(loss_diffs) <- colnames(forecasts)[comp_cols]
      
      results_for_dataset[[paste0("White_Reality_Check_",            metric)]] <- white_reality_check(loss_diffs, block_length, n_simulations = n_boot)
      results_for_dataset[[paste0("Superior_Predictive_Ability_",    metric)]] <- superior_predictive_ability_test(loss_diffs, block_length, n_boot, alpha)
      results_for_dataset[[paste0("Conditional_Predictive_Ability_", metric)]] <- white_reality_check_conditional(loss_diffs, abs(realizations_raw), block_length, n_boot, alpha)
      
      per_model_dataset[[metric]] <- compute_per_model_statistics(loss_diffs, colnames(forecasts)[comp_cols])
      
      if (metric == "MSE") {
        ld_crps <- matrix(NA, nrow(forecasts), K_models)
        for (t in seq_len(nrow(forecasts))) {
          for (i in seq_len(K_models)) {
            ld_crps[t, i] <- compute_crps(
              rnorm(n_crps_samples, forecasts[t, comp_cols[i]], forecast_sd_models[t, i]),
              realizations_raw[t]
            )
          }
        }
        
        bench_crps <- numeric(nrow(forecasts))
        for (t in seq_len(nrow(forecasts))) {
          bench_crps[t] <- compute_crps(
            rnorm(n_crps_samples, benchmark_series[t], forecast_sd_benchmark[t]),
            realizations_raw[t]
          )
        }
        ld_crps_diff <- bench_crps - ld_crps
        results_for_dataset$Cumulative_Distribution_Function <-
          white_reality_check_cdf_approx(ld_crps_diff, block_length, n_boot)
      }
    }
    
    threshold     <- quantile(realizations_raw, zp_quantile)
    zp_matrix     <- cbind(forecasts[, comp_cols], forecasts[, bench_idx])
    new_bench_idx <- ncol(zp_matrix)
    
    lm_zp   <- compute_zp(zp_matrix, forecast_sd_models, threshold, benchmark_col = new_bench_idx)
    zp_diff <- lm_zp[, new_bench_idx] - lm_zp[, 1:K_models]
    results_for_dataset$ZP <- reality_check_zp_test(zp_diff, block_length, n_boot)
    
    lm_klic   <- compute_klic(zp_matrix, forecast_sd_models, benchmark_col = new_bench_idx)
    klic_diff <- lm_klic[, new_bench_idx] - lm_klic[, 1:K_models]
    results_for_dataset$Kullback_Leibler <- kullback_leibler_test(klic_diff, block_length, n_boot)
    
    results_for_dataset$VaR_Backtests <- compute_kupiec(zp_matrix, forecast_sd_models,
                                                        realized = realizations_raw,
                                                        benchmark_col = new_bench_idx)
    
    all_results[[dataset_name]]           <- results_for_dataset
    all_per_model_results[[dataset_name]] <- per_model_dataset
  }
  
  return(list(aggregate_results  = all_results,
              per_model_results   = all_per_model_results))
}
#' @title Flatten Results for Export
#'
#' @description Prepares the comprehensive results list for export (e.g., to Microsoft
#' Excel). Converts all \code{htest} objects and per-model results into flat
#' \code{data.frame} objects, one per dataset.
#'
#' @param comprehensive_results \code{\link[base]{list}} output from
#'   \code{\link{run_comprehensive_erc_analysis}}.
#' @param alpha Numeric significance level used to determine \code{Reject_H0}.
#'   Default is 0.05.
#'
#' @return \code{\link[base]{list}} of data frames, one per dataset, each with columns:
#'   \code{Model}, \code{Test_Type}, \code{P_Value}, \code{Statistic},
#'   \code{Actual_Violations}, \code{Reject_H0}.
#' @seealso \code{\link{run_comprehensive_erc_analysis}}, \code{\link{create_unified_summary}},
#'   \code{\link{generate_comprehensive_report}}
#' @examples
#' \donttest{
#' data(metals)
#' realizations <- list(M = metals[, ncol(metals)])
#' prep_list    <- list(M = list(R_start = 0))
#' f_hat        <- list(list(NULL, NULL, metals))
#' names(f_hat) <- "M"
#' res <- run_comprehensive_erc_analysis(
#'   data_list_prepared = prep_list,
#'   mods_matrix        = matrix(0),
#'   alpha_grid         = 0.05,
#'   window_size        = 20,
#'   y_hat_all          = f_hat,
#'   y_raw              = realizations,
#'   block_length       = 5,
#'   n_boot             = 10,
#'   zp_quantile        = 0.05
#' )
#' excel_data <- extract_and_flatten_results_aggregated(res$aggregate_results, alpha = 0.05)
#' head(excel_data[[1]])
#' }
#' @export
extract_and_flatten_results_aggregated <- function(comprehensive_results, alpha = 0.05) {
  all_data_for_excel <- list()
  
  safe_v <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA)
    return(x)
  }
  
  safe_round <- function(x, digits = 4) {
    val <- safe_v(x)
    if (is.na(val) || !is.numeric(val)) return(NA)
    return(round(as.numeric(val), digits))
  }
  
  for (dataset_name in names(comprehensive_results)) {
    dataset_results     <- comprehensive_results[[dataset_name]]
    dataset_metric_list <- list()
    
    for (metric_name in names(dataset_results)) {
      res    <- dataset_results[[metric_name]]
      if (is.null(res)) next
      df_row <- NULL
      
      if (startsWith(metric_name, "Superior_Predictive_Ability_")) {
        stat_val <- safe_round(res$statistic, 2)
        df_row   <- data.frame(
          Model             = rep("Benchmark", 2),
          Test_Type         = c(paste0("SPA_Consistent_",   metric_name),
                                paste0("SPA_Conservative_", metric_name)),
          P_Value           = c(safe_round(res$p.value), safe_round(res$p_conservative)),
          Statistic         = rep(stat_val, 2),
          Actual_Violations = NA,
          Reject_H0         = c(safe_v(res$p.value        <= alpha),
                                safe_v(res$p_conservative <= alpha)),
          stringsAsFactors  = FALSE
        )
      } else if (metric_name == "VaR_Backtests") {
        df_list <- lapply(names(res), function(nm) {
          item <- res[[nm]]
          data.frame(
            Model             = nm,
            Test_Type         = "Kupiec_UC",
            P_Value           = safe_round(item$p.value),
            Statistic         = safe_round(item$statistic, 2),
            Actual_Violations = safe_v(item$actual_exceedances),
            Reject_H0         = safe_v(!is.na(item$p.value) && item$p.value <= alpha),
            stringsAsFactors  = FALSE
          )
        })
        df_row <- do.call(rbind, df_list)
      } else if (inherits(res, "htest")) {
        df_row <- data.frame(
          Model             = "Benchmark",
          Test_Type         = res$method,
          P_Value           = safe_round(res$p.value),
          Statistic         = safe_round(res$statistic, 2),
          Actual_Violations = NA,
          Reject_H0         = safe_v(!is.na(res$p.value) && res$p.value <= alpha),
          stringsAsFactors  = FALSE
        )
      }
      
      if (!is.null(df_row)) dataset_metric_list[[metric_name]] <- df_row
    }
    
    if (length(dataset_metric_list) > 0) {
      aggregated_df <- do.call(rbind, dataset_metric_list)
      rownames(aggregated_df) <- NULL
      sheet_name <- substr(gsub(" ", "_", dataset_name), 1, 31)
      all_data_for_excel[[sheet_name]] <- aggregated_df
    }
  }
  
  return(all_data_for_excel)
}

#' @title Create Unified Summary
#'
#' @description Creates a summary data frame consolidating p-values and conclusions for
#' all statistical tests across all datasets and error metrics. Covers White's Reality
#' Check (WRC), Superior Predictive Ability (SPA), Conditional Predictive Ability (CPA),
#' ZP Quantile Loss test, and Kullback-Leibler (KLIC) test.
#'
#' @param comprehensive_results \code{\link[base]{list}} output from
#'   \code{\link{run_comprehensive_erc_analysis}}.
#' @param alpha Numeric significance level for determining conclusions. Default is 0.05.
#'
#' @return \code{\link[base]{list}} containing a data frame named \code{summary} with
#'   columns: \code{Dataset}, \code{Test}, \code{P_Value}, \code{Statistic},
#'   \code{Conclusion} ("H0 rejected" or "H0 accepted").
#'
#' @seealso \code{\link{run_comprehensive_erc_analysis}},
#'   \code{\link{generate_comprehensive_report}},
#'   \code{\link{compute_per_model_statistics}}
#'   
#' @references
#' White, H. (2000). A reality check for data snooping. \emph{Econometrica}, 68(5),
#' 1097--1126. \doi{10.1111/1468-0262.00152}
#'
#' Hansen, P. R. (2005). A Test for Superior Predictive Ability.
#' \emph{Journal of Business & Economic Statistics}, 23(4), 365--380.
#' \doi{10.1198/073500105000000063}
#'
#' Giacomini, R., & White, H. (2006). Tests of Conditional Predictive Ability.
#' \emph{Econometrica}, 74(6), 1545--1578. \doi{10.1111/j.1468-0262.2006.00718.x}
#'
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional confidence
#' interval accuracy tests. \emph{Journal of Econometrics}, 135(1--2), 187--228.
#' \doi{10.1016/j.jeconom.2005.07.026}
#'
#' @examples
#' \donttest{
#' data(metals)
#' realizations <- list(M = metals[, ncol(metals)])
#' prep_list    <- list(M = list(R_start = 0))
#' f_hat        <- list(list(NULL, NULL, metals))
#' names(f_hat) <- "M"
#' res <- run_comprehensive_erc_analysis(
#'   data_list_prepared = prep_list,
#'   mods_matrix        = matrix(0),
#'   alpha_grid         = 0.05,
#'   window_size        = 20,
#'   y_hat_all          = f_hat,
#'   y_raw              = realizations,
#'   block_length       = 5,
#'   n_boot             = 10,
#'   zp_quantile        = 0.05
#' )
#' summary_table <- create_unified_summary(res$aggregate_results)
#' print(summary_table$summary)
#' }
#' @export
create_unified_summary <- function(comprehensive_results, alpha = 0.05) {
  
  safe_extract <- function(obj, field_name) {
    if (is.null(obj)) return(NA)
    if (field_name == "p_value" && !is.null(obj$p.value)) {
      val <- obj$p.value
    } else if (field_name == "test_statistic" && !is.null(obj$statistic)) {
      val <- obj$statistic
    } else {
      return(NA)
    }
    if (is.na(val) || !is.numeric(val)) return(NA)
    
    if (field_name == "p_value") {
      return(as.numeric(val))
    }
    
    return(round(as.numeric(val), 4))
  }
  
  all_summaries <- list()
  
  for (dataset_name in names(comprehensive_results)) {
    dataset_results <- comprehensive_results[[dataset_name]]
    metrics         <- c("MSE", "MAE", "MASE")
    
    for (metric in metrics) {
      test_names <- c(
        paste0("White_Reality_Check_",            metric),
        paste0("Superior_Predictive_Ability_",    metric),
        paste0("Conditional_Predictive_Ability_", metric)
      )
      
      for (test_name in test_names) {
        if (test_name %in% names(dataset_results)) {
          res_obj    <- dataset_results[[test_name]]
          p_val      <- safe_extract(res_obj, "p_value")
          stat_val   <- safe_extract(res_obj, "test_statistic")
          conclusion <- if (!is.na(p_val)) {
            if (p_val <= alpha) "H0 rejected" else "Fail to reject H0"
          } else "Insufficient data"
          
          all_summaries[[paste0(dataset_name, "_", test_name)]] <- data.frame(
            Dataset = dataset_name, Test = test_name,
            P_Value = p_val, Statistic = stat_val,
            Conclusion = conclusion, stringsAsFactors = FALSE
          )
          
          if (startsWith(test_name, "Superior_Predictive_Ability_") &&
              !is.null(res_obj$p_conservative)) {
            
            p_cons <- as.numeric(res_obj$p_conservative)
            
            conc_cons <- if (!is.na(p_cons)) {
              if (p_cons <= alpha) "H0 rejected" else "Fail to reject H0"
            } else "Insufficient data"
            
            test_cons <- paste0(test_name, "_Conservative")
            
            all_summaries[[paste0(dataset_name, "_", test_cons)]] <- data.frame(
              Dataset = dataset_name, Test = test_cons,
              P_Value = p_cons, Statistic = stat_val,
              Conclusion = conc_cons, stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    dist_tests <- c("ZP", "Kullback_Leibler", "Cumulative_Distribution_Function")
    
    for (test_name in dist_tests) {
      if (test_name %in% names(dataset_results)) {
        res_obj <- dataset_results[[test_name]]
        p_val   <- safe_extract(res_obj, "p_value")
        
        all_summaries[[paste0(dataset_name, "_", test_name)]] <- data.frame(
          Dataset          = dataset_name,
          Test             = test_name,
          P_Value          = p_val,
          Statistic        = safe_extract(res_obj, "test_statistic"),
          Conclusion       = if (!is.na(p_val)) (if (p_val <= alpha) "H0 rejected" else "Fail to reject H0") else "N/A",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(all_summaries) == 0) return(list(summary = data.frame()))
  
  final_summary_df           <- do.call(rbind, all_summaries)
  rownames(final_summary_df) <- NULL
  
  return(list(summary = final_summary_df))
}
#' @title Generate Comprehensive Markdown Report
#'
#' @description Generates an automatic summary report in Markdown format covering all
#' error metrics (MSE, MAE, MASE) and distributional tests (ZP, KLIC, Kupiec).
#' For the ZP and KLIC sections, the report lists superior competing forecasts (i.e.
#' those whose forecasts are found to be more accurate than the benchmark forecast)
#'  or states that no superior forecasts were found. For the Kupiec section, forecasts 
#'  with correct VaR coverage (\code{Reject_H0 == FALSE}) are listed, or a message
#'  is printed if none passed.
#'
#'\strong{Technical Abbreviations:}
#' \itemize{
#'   \item \strong{WRC:} White's Reality Check (White, 2000). Tests whether any competing
#'     forecast has lower expected loss than the benchmark forecast; controls family-wise error rate.
#'   \item \strong{SPA:} Superior Predictive Ability test (Hansen, 2005). A studentized
#'     extension of WRC with improved power that corrects for irrelevant forecasts.
#'   \item \strong{CPA:} Conditional Predictive Ability test (Giacomini & White, 2006).
#'     Tests whether loss differentials are predictable by a conditioning variable.
#'   \item \strong{ZP:} Quantile Loss test (Corradi & Swanson, 2006). Evaluates whether
#'     any competing forecast better calibrates the probability of a left-tail event
#'     defined by the \code{zp_quantile} threshold.
#'   \item \strong{KLIC:} Kullback-Leibler Information Criterion based density test
#'     (Corradi & Swanson, 2006). Selects the forecast whose predictive density is closest
#'     to the true density in terms of KLIC distance, evaluated via Negative
#'     Log-Likelihood Scores (NLS) under a Gaussian predictive density assumption.
#'   \item \strong{CRPS:} Continuous Ranked Probability Score (Gneiting & Raftery, 2007).
#'     Jointly rewards calibration and sharpness of the predictive distribution.
#'   \item \strong{UC:} Kupiec Unconditional Coverage test (Kupiec, 1995).
#'   \item \strong{MSE:} Mean Squared Error.
#'   \item \strong{MAE:} Mean Absolute Error.
#'   \item \strong{MASE:} Mean Absolute Scaled Error.
#' }
#' @references
#' White, H. (2000). A reality check for data snooping. \emph{Econometrica}, 68(5),
#' 1097--1126. \doi{10.1111/1468-0262.00152}
#'
#' Hansen, P. R. (2005). A Test for Superior Predictive Ability.
#' \emph{Journal of Business & Economic Statistics}, 23(4), 365--380.
#' \doi{10.1198/073500105000000063}
#'
#' Giacomini, R., & White, H. (2006). Tests of Conditional Predictive Ability.
#' \emph{Econometrica}, 74(6), 1545--1578. \doi{10.1111/j.1468-0262.2006.00718.x}
#'
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional confidence
#' interval accuracy tests. \emph{Journal of Econometrics}, 135(1--2), 187--228.
#' \doi{10.1016/j.jeconom.2005.07.026}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its recent
#' extensions. In \emph{Festschrift in honor of Halbert L. White}.
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules, Prediction,
#' and Estimation. \emph{Journal of the American Statistical Association}, 102(477),
#' 359--378. \doi{10.1198/016214506000001437}
#'
#' Kupiec, P. H. (1995). Techniques for Verifying the Accuracy of Risk Measurement
#' Models. \emph{The Journal of Derivatives}, 3(2), 173--184.
#' \doi{10.3905/jod.1995.407942}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' Diebold, F. X., & Mariano, R. S. (1995). Comparing Predictive Accuracy.
#' \emph{Journal of Business & Economic Statistics}, 13(3), 253--263.
#' \doi{10.1080/07350015.1995.10524599}
#'
#' @param summary_df \code{\link[base]{data.frame}} unified summary from
#'   \code{\link{create_unified_summary}}.
#' @param zp_models_df \code{\link[base]{data.frame}} with columns \code{Model},
#'   \code{P_Value}, and \code{Dataset}. Models with \code{P_Value <= alpha} are
#'   considered superior and listed in the report; if none are found, a message
#'   \"No superior models found\" is printed. Pass all competing models to ensure
#'   complete and unbiased reporting.
#' @param klic_models_df \code{\link[base]{data.frame}} with columns \code{Model},
#'   \code{P_Value}, and \code{Dataset}. Models with \code{P_Value <= alpha} are
#'   considered superior and listed in the report; if none are found, a message
#'   \"No superior models found\" is printed. Pass all competing models to ensure
#'   complete and unbiased reporting.
#' @param kupiec_models_df \code{\link[base]{data.frame}} with columns \code{Model},
#'   \code{Reject_H0}, and \code{Dataset}. Contains all models tested under the Kupiec
#'   UC test. Models with \code{Reject_H0 == FALSE} are considered to have correct VaR
#'   coverage and are listed in the report; models with \code{Reject_H0 == TRUE} are
#'   excluded. Pass all competing models to ensure complete and unbiased reporting.
#' @param dataset_name \code{\link[base]{character}} the name of the dataset to be used
#'   in the report header.
#' @param alpha \code{\link[base]{numeric}} significance level (default 0.05).
#'
#' @return \code{\link[base]{character}} string in Markdown format.
#'
#' @examples
#' \donttest{
#' data(metals)
#' ds_name       <- "Dataset1"
#' prep_list     <- list(Dataset1 = list(R_start = 0))
#' realizations  <- list(Dataset1 = metals[, ncol(metals)])
#' f_hat         <- list(list(NULL, NULL, metals))
#' names(f_hat)  <- ds_name
#' res <- run_comprehensive_erc_analysis(
#'   data_list_prepared = prep_list,
#'   mods_matrix        = matrix(0),
#'   alpha_grid         = 0.05,
#'   window_size        = 20,
#'   y_hat_all          = f_hat,
#'   y_raw              = realizations,
#'   block_length       = 5,
#'   n_boot             = 10,
#'   zp_quantile        = 0.05
#' )
#' unified_summ <- create_unified_summary(res$aggregate_results)
#'
#' zp_table     <- data.frame(Model   = colnames(metals)[1:(ncol(metals) - 1)],
#'                             P_Value = 0.1, Dataset = ds_name)
#' klic_table   <- zp_table
#' kupiec_table <- data.frame(Model     = colnames(metals)[1:(ncol(metals) - 1)],
#'                             Reject_H0 = rep(FALSE, ncol(metals) - 1),
#'                             Dataset   = ds_name)
#'
#' report <- generate_comprehensive_report(
#'   summary_df       = unified_summ$summary,
#'   zp_models_df     = zp_table,
#'   klic_models_df   = klic_table,
#'   kupiec_models_df = kupiec_table,
#'   dataset_name     = ds_name
#' )
#' cat(report)
#' }
#' @export
generate_comprehensive_report <- function(summary_df, zp_models_df, klic_models_df,
                                          kupiec_models_df, dataset_name, alpha = 0.05) {
  if (is.list(summary_df) && !is.data.frame(summary_df)) summary_df <- summary_df$summary
  summary_df       <- as.data.frame(summary_df)
  zp_models_df     <- as.data.frame(zp_models_df)
  klic_models_df   <- as.data.frame(klic_models_df)
  kupiec_models_df <- as.data.frame(kupiec_models_df)
  
  safe_format <- function(val) {
    if (length(val) == 0 || is.na(val[1])) return("N/A")
    
    p <- as.numeric(val[1])
    
    if (p < 0.0001) {
      return("p < 0.0001")
    }
    
    return(paste0("p = ", formatC(p, format = "f", digits = 4)))
  }
  
  get_conc <- function(p) {
    if (length(p) == 0 || is.na(p[1]) || p[1] == 1) return("Accepted")
    if (as.numeric(p[1]) <= alpha) "**Rejected**" else "Accepted"
  }
  
  summary_data <- if (NROW(summary_df) > 0) summary_df[summary_df$Dataset == dataset_name, ] else data.frame()
  
  report <- paste0("## Forecasting Performance Summary: ", dataset_name, "\n\n")
  report <- paste0(report, "### 1. Reality Check Results (WRC, SPA & CPA)\n")
  report <- paste0(report, "H0: No competing forecast is more accurate than the benchmark forecast.\n\n")
  
  for (m in c("MSE", "MAE", "MASE")) {
    report <- paste0(report, "#### **Metric: ", m, "**\n")
    
    wrc_p <- summary_data$P_Value[summary_data$Test == paste0("White_Reality_Check_",            m)]
    spa_p <- summary_data$P_Value[summary_data$Test == paste0("Superior_Predictive_Ability_",    m)]
    cpa_p <- summary_data$P_Value[summary_data$Test == paste0("Conditional_Predictive_Ability_", m)]
    
    report <- paste0(report, "* **White's Reality Check (WRC):** H0 ",          get_conc(wrc_p), " (", safe_format(wrc_p), ")\n")
    report <- paste0(report, "* **Superior Predictive Ability (SPA):** H0 ",     get_conc(spa_p), " (", safe_format(spa_p), ")\n")
    report <- paste0(report, "* **Conditional Predictive Ability (CPA):** H0 ", get_conc(cpa_p), " (", safe_format(cpa_p), ")\n\n")
  }
  
  report <- paste0(report, "---\n### 2. Pairwise & Distributional Evaluation\n\n")
  
  zp_p   <- summary_data$P_Value[summary_data$Test == "ZP"]
  klic_p <- summary_data$P_Value[summary_data$Test == "Kullback_Leibler"]
  
  report <- paste0(report, "#### **ZP Quantile Test** (Corradi & Swanson, 2006)\n")
  if (length(zp_p) > 0 && !is.na(zp_p[1])) {
    report <- paste0(report, "* H0: ", get_conc(zp_p), " (", safe_format(zp_p), ")\n")
  } else {
    report <- paste0(report, "* P-Value: N/A\n")
  }
  if (nrow(zp_models_df) > 0 && "Model" %in% colnames(zp_models_df) &&
      "P_Value" %in% colnames(zp_models_df)) {
    zp_superior <- zp_models_df$Model[zp_models_df$P_Value <= alpha]
    if (length(zp_superior) > 0) {
      report <- paste0(report, "* Superior models: ", paste(zp_superior, collapse = ", "), "\n\n")
    } else {
      report <- paste0(report, "* Superior models: No superior models found.\n\n")
    }
  } else {
    report <- paste0(report, "\n")
  }
  
  report <- paste0(report, "#### **Kullback-Leibler Information Criterion (KLIC)** (Corradi & Swanson, 2006)\n")
  if (length(klic_p) > 0 && !is.na(klic_p[1])) {
    report <- paste0(report, "* H0: ", get_conc(klic_p), " (", safe_format(klic_p), ")\n")
  } else {
    report <- paste0(report, "* P-Value: N/A\n")
  }
  if (nrow(klic_models_df) > 0 && "Model" %in% colnames(klic_models_df) &&
      "P_Value" %in% colnames(klic_models_df)) {
    klic_superior <- klic_models_df$Model[klic_models_df$P_Value <= alpha]
    if (length(klic_superior) > 0) {
      report <- paste0(report, "* Superior models: ", paste(klic_superior, collapse = ", "), "\n\n")
    } else {
      report <- paste0(report, "* Superior models: No superior models found.\n\n")
    }
  } else {
    report <- paste0(report, "\n")
  }
  
  report <- paste0(report, "---\n### 3. Coverage Analysis (Kupiec Unconditional Test)\n\n")
  passed_models <- kupiec_models_df$Model[kupiec_models_df$Reject_H0 == FALSE]
  failed_models <- kupiec_models_df$Model[kupiec_models_df$Reject_H0 == TRUE]
  
  if (length(passed_models) > 0) {
    report <- paste0(report, "* Forecasts with correct coverage (H0 not rejected): ",
                     paste(passed_models, collapse = ", "), "\n")
  } else {
    report <- paste0(report, "* No forecasts passed the coverage test.\n")
  }
  
  if (length(failed_models) > 0) {
    report <- paste0(report, "* Forecasts with incorrect coverage (H0 rejected): ",
                     paste(failed_models, collapse = ", "), "\n")
  }
  
  return(report)
}