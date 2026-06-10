#' @importFrom stats qnorm pchisq var dnorm pnorm pt
NULL

#' @title Estimate Forecast Variance via Rolling Window
#'
#' @description Estimates forecast variance from historical forecast errors relative to
#' the benchmark using a rolling window. Used to approximate the time-varying predictive
#' standard deviation for each competing model forecast, required by \code{\link{compute_klic}},
#' \code{\link{compute_zp}}, and \code{\link{compute_kupiec}}.
#'
#' @param forecast_matrix \code{\link[base]{matrix}} of dimension \code{P x K_total},
#'   where \code{P} is the number of forecast periods and \code{K_total} is the total
#'   number of columns (competing forecasts plus the benchmark).
#' @param benchmark_col Index or name of the benchmark column. Defaults to the last
#'   column.
#' @param window_size \code{\link[base]{integer}} rolling window size. For the first
#'   \code{window_size} periods, the full available history is used instead (expanding
#'   window). From period \code{window_size + 1} onwards, a rolling window of exactly
#'   \code{window_size} observations is used.
#'
#' @details
#' For each competing forecast \code{k} and period \code{t}, the forecast error is defined
#' as \eqn{e_{t,k} = \text{benchmark}_t - \hat{y}_{t,k}}. The variance of these errors
#' is estimated over a rolling window:
#' \itemize{
#'   \item If \code{t < window_size}: variance is computed over observations
#'     \code{1:t} (expanding window).
#'   \item If \code{t >= window_size}: variance is computed over observations
#'     \code{max(1, t - window_size):t} (rolling window of size \code{window_size}).
#' }
#' For \code{t = 1} the variance of a single observation is undefined (\code{NA}).
#' Estimated variances that are \code{NA}, zero, or negative are replaced by
#' \code{1e-6} to ensure numerical stability in downstream computations.
#' The benchmark column in the returned matrix is set to zero throughout.
#'
#' @return \code{\link[base]{matrix}} of dimension \code{P x K_total} containing
#'   estimated variances. Columns correspond to the same forecasts as
#'   \code{forecast_matrix}; the benchmark column contains zeros.
#'
#' @seealso \code{\link{compute_klic}}, \code{\link{compute_zp}},
#'   \code{\link{compute_kupiec}}
#'
#' @examples
#' data(metals)
#' forecast_variance <- estimate_forecast_variance(metals, benchmark_col = 15,
#'                                                 window_size = 20)
#' head(forecast_variance)
#' @export
estimate_forecast_variance <- function(forecast_matrix,
                                       benchmark_col = ncol(forecast_matrix),
                                       window_size   = 20) {
  if (is.character(benchmark_col))
    benchmark_col <- which(colnames(forecast_matrix) == benchmark_col)
  if (length(benchmark_col) == 0)
    stop("The specified benchmark_col was not found in the matrix.")
  
  P          <- nrow(forecast_matrix)
  K_total    <- ncol(forecast_matrix)
  model_cols <- setdiff(seq_len(K_total), benchmark_col)
  K_models   <- length(model_cols)
  
  benchmark         <- forecast_matrix[, benchmark_col]
  forecast_variance <- matrix(0, P, K_total)
  colnames(forecast_variance) <- colnames(forecast_matrix)
  
  for (i in seq_len(K_models)) {
    k      <- model_cols[i]
    errors <- benchmark - forecast_matrix[, k]
    for (t in seq_len(P)) {
      if (t < window_size) {
        forecast_variance[t, k] <- var(errors[1:t], na.rm = TRUE)
      } else {
        start_idx               <- max(1L, t - window_size)
        forecast_variance[t, k] <- var(errors[start_idx:t], na.rm = TRUE)
      }
    }
  }
  
  for (i in seq_len(K_models)) {
    k <- model_cols[i]
    forecast_variance[
      is.na(forecast_variance[, k]) | forecast_variance[, k] <= 0, k
    ] <- 1e-6
  }
  
  return(forecast_variance)
}

#' @title Value-at-Risk (VaR) Unconditional Coverage Test (Kupiec)
#'
#' @description Performs Kupiec's (1995) Unconditional Coverage (UC) test for evaluating
#' Value-at-Risk (VaR) forecasts from competing forecast against realized values.
#'
#' \strong{Hypotheses:}
#' \itemize{
#'   \item \strong{H0:} The forecast correctly captures VaR — violations occur with the
#'     expected frequency \code{alpha}.
#'   \item \strong{H1:} The forecast fails to correctly capture VaR — the observed
#'     frequency of violations differs significantly from \code{alpha}.
#' }
#'
#' @param forecast_matrix \code{\link[base]{matrix}} of dimension \code{P x K_total}.
#'   Columns contain point forecasts for each model; the benchmark column supplies
#'   the realized values.
#' @param forecast_sd_models \code{\link[base]{matrix}} of dimension \code{P x K},
#'   where \code{K = K_total - 1}. Contains time-varying forecast standard deviations,
#'   typically from \code{\link{estimate_forecast_variance}}.
#' @param benchmark_col Index or name of the benchmark column. Defaults to the last
#'   column.
#' @param alpha \code{\link[base]{numeric}} VaR significance level (e.g.,
#'   \code{0.05} for 95\% VaR). A violation occurs when the realized value falls below
#'   the estimated VaR.
#'
#' @details
#' For each competing forecast \code{k}, the VaR at level \code{alpha} is:
#' \deqn{VaR_{t,k} = \hat{y}_{t,k} + \Phi^{-1}(\alpha) \cdot \hat{\sigma}_{t,k}}
#' where \eqn{\Phi^{-1}} is the standard normal quantile function. A violation occurs
#' when the realized value falls below \eqn{VaR_{t,k}}. The likelihood-ratio statistic
#' \eqn{LR_{UC}} follows a \eqn{\chi^2(1)} distribution under H0 (Kupiec, 1995).
#' \strong{Failing to reject H0} (large p-value) indicates correctly calibrated VaR;
#' \strong{rejecting H0} (small p-value) indicates the forecast under- or over-estimates
#' tail risk.
#'
#' @return A named list (one element per competing forecast) of \code{htest} objects,
#'   each containing:
#' \describe{
#'   \item{\code{statistic}}{The LR-UC test statistic (\eqn{\chi^2}-distributed under
#'     H0).}
#'   \item{\code{p.value}}{P-value from the \eqn{\chi^2(1)} distribution. A large
#'     p-value indicates correctly calibrated VaR coverage.}
#'   \item{\code{actual_exceedances}}{Observed number of VaR violations.}
#'   \item{\code{expected}}{Expected number of violations
#'     (\code{P * alpha}).}
#' }
#'
#' @seealso \code{\link{estimate_forecast_variance}},
#'   \code{\link{run_comprehensive_erc_analysis}}
#'
#' @references
#' Kupiec, P. H. (1995). Techniques for Verifying the Accuracy of Risk Measurement
#' Models. \emph{The Journal of Derivatives}, 3(2), 173--184.
#' \doi{10.3905/jod.1995.407942}
#'
#' @examples
#' data(metals)
#' benchmark_col      <- 15
#' K_total            <- ncol(metals)
#' comp_cols          <- setdiff(seq_len(K_total), benchmark_col)
#' forecast_variance  <- estimate_forecast_variance(metals,
#'                         benchmark_col = benchmark_col, window_size = 20)
#' forecast_sd_models <- sqrt(forecast_variance[, comp_cols])
#' coverage_results   <- compute_kupiec(metals, forecast_sd_models,
#'                         benchmark_col = benchmark_col, alpha = 0.05)
#' print(coverage_results[[1]])
#' @export
compute_kupiec <- function(forecast_matrix, forecast_sd_models,
                           benchmark_col = ncol(forecast_matrix),
                           alpha   = 0.05) {
  P          <- nrow(forecast_matrix)
  K_total    <- ncol(forecast_matrix)
  model_cols <- setdiff(seq_len(K_total), benchmark_col)
  K_models   <- length(model_cols)
  
  if (ncol(forecast_sd_models) != K_models)
    stop("forecast_sd_models must have K = K_total - 1 columns (one per competing model).")
  
  benchmark  <- forecast_matrix[, benchmark_col]
  VaR_matrix <- forecast_matrix
  
  for (i in seq_len(K_models)) {
    k <- model_cols[i]
    VaR_matrix[, k] <- forecast_matrix[, k] +
      qnorm(alpha) * forecast_sd_models[, i]
  }
  
  kupiec_results <- list()
  
  for (i in seq_len(K_models)) {
    k          <- model_cols[i]
    violations <- sum(benchmark < VaR_matrix[, k], na.rm = TRUE)
    expected   <- P * alpha
    
    p_hat <- violations / P
    if (p_hat == 0) p_hat <- 1e-10
    if (p_hat == 1) p_hat <- 1 - 1e-10
    
    LR_UC_stat <- -2 * log(
      (alpha^violations * (1 - alpha)^(P - violations)) /
        (p_hat^violations    * (1 - p_hat)^(P - violations))
    )
    
    p_value_UC <- 1 - pchisq(LR_UC_stat, df = 1)
    model_name <- colnames(forecast_matrix)[k]
    
    res <- list(
      statistic          = c("LR-UC" = LR_UC_stat),
      p.value            = p_value_UC,
      method             = paste("Kupiec UC Test for forecast:", model_name),
      data.name          = paste0("Benchmark: column ", benchmark_col,
                                  if (!is.null(colnames(forecast_matrix)))
                                    paste0(" (", colnames(forecast_matrix)[benchmark_col],
                                           ")")
                                  else ""),
      null.value         = c("expected violation rate" = alpha),
      alternative        = "violation rate differs from expected",
      actual_exceedances = violations,
      expected           = expected
    )
    class(res) <- "htest"
    kupiec_results[[model_name]] <- res
  }
  
  return(kupiec_results)
}

#' @title Compute Kullback-Leibler Information Criterion (KLIC) Negative Log-Likelihood Scores
#'
#' @description Computes the per-period Negative Log-Likelihood Score (NLS) loss matrix
#' under a Gaussian predictive density assumption. The NLS is the loss function
#' corresponding to minimisation of the Kullback-Leibler Information Criterion (KLIC)
#' distance from the true density (Corradi & Swanson, 2006).
#'
#' @param forecast_matrix \code{\link[base]{matrix}} of dimension \code{P x K_total}.
#'   The benchmark column supplies the realized values \eqn{y_t}.
#' @param forecast_sd_models \code{\link[base]{matrix}} of dimension
#'   \code{P x (K_total - 1)}, containing time-varying forecast standard deviations,
#'   typically from \code{\link{estimate_forecast_variance}}.
#' @param benchmark_col Index or name of the benchmark column. Defaults to the last
#'   column.
#'
#' @details
#' For each competing forecast \code{k} and period \code{t}:
#' \deqn{NLS_{t,k} = -\log \phi(y_t \mid \hat{y}_{t,k},\, \hat{\sigma}_{t,k})}
#' where \eqn{\phi} denotes the Gaussian density, \eqn{y_t} is the realized value,
#' \eqn{\hat{y}_{t,k}} is the point forecast, and \eqn{\hat{\sigma}_{t,k}} is the
#' forecast standard deviation. Minimising the average NLS is equivalent to minimising
#' the KLIC distance between the forecast's predictive density and the true density
#' (Corradi & Swanson, 2006). \strong{Lower NLS values are better.} The benchmark
#' column in the returned matrix is set to zero.
#'
#' @return \code{\link[base]{matrix}} of dimension \code{P x K_total} containing NLS
#'   values. Lower values indicate better density forecast accuracy. The benchmark
#'   column is set to zero.
#'
#' @seealso \code{\link{kullback_leibler_test}},
#'   \code{\link{estimate_forecast_variance}}
#'
#' @references
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional confidence
#' interval accuracy tests. \emph{Journal of Econometrics}, 135(1--2), 187--228.
#' \doi{10.1016/j.jeconom.2005.07.026}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its
#' recent extensions. In \emph{Festschrift in honor of Halbert L. White}.
#'
#' @examples
#' data(metals)
#' benchmark_col      <- 15
#' K_total            <- ncol(metals)
#' comp_cols          <- setdiff(seq_len(K_total), benchmark_col)
#' forecast_variance  <- estimate_forecast_variance(metals,
#'                         benchmark_col = benchmark_col)
#' forecast_sd_models <- sqrt(forecast_variance[, comp_cols])
#' klic_loss <- compute_klic(metals, forecast_sd_models,
#'                           benchmark_col = benchmark_col)
#' head(klic_loss)
#' @export
compute_klic <- function(forecast_matrix, forecast_sd_models,
                         benchmark_col = ncol(forecast_matrix)) {
  P          <- nrow(forecast_matrix)
  K_total    <- ncol(forecast_matrix)
  model_cols <- setdiff(seq_len(K_total), benchmark_col)
  K_models   <- length(model_cols)
  
  if (ncol(forecast_sd_models) != K_models)
    stop("forecast_sd_models must have K = K_total - 1 columns (one per competing model forecast).")
  
  benchmark <- forecast_matrix[, benchmark_col]
  nls_loss  <- matrix(NA, P, K_total)
  colnames(nls_loss) <- colnames(forecast_matrix)
  
  for (i in seq_len(K_models)) {
    k <- model_cols[i]
    nls_loss[, k] <- -dnorm(benchmark,
                            mean = forecast_matrix[, k],
                            sd   = forecast_sd_models[, i],
                            log  = TRUE)
  }
  nls_loss[, benchmark_col] <- 0
  
  return(nls_loss)
}

#' @title Compute ZP Quantile Loss
#'
#' @description Computes the per-period ZP quantile loss matrix based on the squared
#' difference between the indicator of a tail event and the forecast's predicted
#' probability of that event (Corradi & Swanson, 2006, eq. 7).
#'
#' @param forecast_matrix \code{\link[base]{matrix}} of dimension \code{P x K_total}.
#'   The benchmark column supplies the realized values \eqn{y_t}.
#' @param forecast_sd_models \code{\link[base]{matrix}} of dimension \code{P x K},
#'   where \code{K = K_total - 1}. Contains time-varying forecast standard deviations,
#'   typically from \code{\link{estimate_forecast_variance}}.
#' @param threshold \code{\link[base]{numeric}} tail threshold \eqn{\tau}. The ZP loss
#'   measures how well each model predicts the probability of \eqn{y_t \leq \tau}.
#'   Typically set to a low quantile of the realized series, e.g.,
#'   \code{quantile(realized, 0.05)} for the 5th-percentile left tail.
#'   In \code{\link{run_comprehensive_erc_analysis}} this is computed automatically
#'   as \code{quantile(realizations, zp_quantile)}.
#' @param benchmark_col Index or name of the benchmark column. Defaults to the last
#'   column.
#'
#' @details
#' For each competing forecast \code{k} and period \code{t}:
#' \deqn{ZP_{t,k} = \left(\mathbf{1}(y_t \leq \tau) -
#'   \Phi\!\left(\frac{\tau - \hat{y}_{t,k}}{\hat{\sigma}_{t,k}}\right)\right)^2}
#' where \eqn{y_t} is the realized value, \eqn{\tau} is the threshold,
#' \eqn{\hat{y}_{t,k}} is the point forecast, and \eqn{\hat{\sigma}_{t,k}} is the
#' forecast standard deviation. \strong{Lower ZP values are better.}
#'
#' Choosing \eqn{\tau} at the 5th percentile focuses evaluation on whether forecasts
#' correctly predict the risk of falling into the worst 5\% of outcomes. The benchmark
#' column is assigned a point-mass predictive distribution (\eqn{\hat{\sigma} = 10^{-6}}),
#' which approximates the Brier score for the tail indicator and serves as a conservative
#' reference. When the benchmark is the Historical Average (HA), the ZP test thus
#' evaluates whether any competing forecast's calibrated tail probability outperforms the
#' HA's point prediction of the tail event.
#'
#' @return \code{\link[base]{matrix}} of dimension \code{P x K_total} containing ZP
#'   loss values. Lower values indicate better left-tail probability calibration.
#'   The benchmark column uses \eqn{\hat{\sigma} = 10^{-6}}.
#'
#' @seealso \code{\link{reality_check_zp_test}},
#'   \code{\link{estimate_forecast_variance}}
#'
#' @references
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional confidence
#' interval accuracy tests. \emph{Journal of Econometrics}, 135(1--2), 187--228.
#' \doi{10.1016/j.jeconom.2005.07.026}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its
#' recent extensions. In \emph{Festschrift in honor of Halbert L. White}.
#'
#' @examples
#' data(metals)
#' benchmark_col      <- 15
#' K_total            <- ncol(metals)
#' comp_cols          <- setdiff(seq_len(K_total), benchmark_col)
#' forecast_variance  <- estimate_forecast_variance(metals,
#'                         benchmark_col = benchmark_col)
#' forecast_sd_models <- sqrt(forecast_variance[, comp_cols])
#' threshold_val      <- quantile(metals[, benchmark_col], 0.10)
#' zp_loss <- compute_zp(metals, forecast_sd_models,
#'                       threshold     = threshold_val,
#'                       benchmark_col = benchmark_col)
#' head(zp_loss)
#' @export
compute_zp <- function(forecast_matrix, forecast_sd_models, threshold,
                       benchmark_col = ncol(forecast_matrix)) {
  P          <- nrow(forecast_matrix)
  K_total    <- ncol(forecast_matrix)
  model_cols <- setdiff(seq_len(K_total), benchmark_col)
  K_models   <- length(model_cols)
  
  if (ncol(forecast_sd_models) != K_models)
    stop("forecast_sd_models must have K = K_total - 1 columns (one per competing model).")
  
  benchmark <- forecast_matrix[, benchmark_col]
  zp_loss   <- matrix(NA, P, K_total)
  colnames(zp_loss) <- colnames(forecast_matrix)
  
  for (t in seq_len(P)) {
    actual_violation <- as.numeric(benchmark[t] <= threshold)
    for (i in seq_len(K_models)) {
      k      <- model_cols[i]
      prob_k <- pnorm(threshold,
                      mean = forecast_matrix[t, k],
                      sd   = forecast_sd_models[t, i])
      zp_loss[t, k] <- (actual_violation - prob_k)^2
    }
    prob_benchmark            <- pnorm(threshold,
                                       mean = forecast_matrix[t, benchmark_col],
                                       sd   = 1e-6)
    zp_loss[t, benchmark_col] <- (actual_violation - prob_benchmark)^2
  }
  
  return(zp_loss)
}

#' @title Per-Model Diebold-Mariano Test (HAC + MBB Bootstrap)
#'
#' @description Performs a Diebold-Mariano (1995) test for each competing forecast
#' separately, testing whether the predictive accuracy of that forecast differs
#' significantly from the benchmark forecast. The test statistic is the sample
#' mean loss differential standardised by a Newey-West HAC standard error;
#' p-values are provided both analytically and via Moving Block Bootstrap (MBB).
#' The direction of the alternative hypothesis is controlled by the
#' \code{H1} argument: \code{"same"} for a two-sided test (\eqn{H_1: \bar{d}_k \neq 0}),
#' \code{"more"} for the one-sided test that forecast \eqn{k} is more accurate than
#' the benchmark (\eqn{H_1: \bar{d}_k > 0}), and \code{"less"} for the one-sided
#' test that forecast \eqn{k} is less accurate than the benchmark
#' (\eqn{H_1: \bar{d}_k < 0}).
#' @details
#' For each forecast \eqn{k}, the loss differential series is:
#' \deqn{d_{t,k} = g(e_{0,t}) - g(e_{k,t})}
#' where \eqn{g(\cdot)} is the loss function used to construct \code{loss_differences}.
#' In the standard workflow of this package
#' (\code{\link{run_comprehensive_erc_analysis}}), \eqn{g} is the squared error loss:
#' \deqn{d_{t,k} = (y_t - \hat{y}_{t,0})^2 - (y_t - \hat{y}_{t,k})^2}
#' The Diebold-Mariano test statistic is:
#' \deqn{DM_k = \frac{\bar{d}_k}{\hat{SE}_{HAC,k}}}
#' For multi-step forecasts (\eqn{h > 1}), the Harvey, Leybourne & Newbold (1997)
#' small-sample correction is applied:
#' \deqn{DM^*_k = DM_k \times \sqrt{\frac{T + 1 - 2h + \frac{1}{T}h(h-1)}{T}}}
#' For \eqn{h = 1} the correction reduces to \eqn{\sqrt{(T-1)/T}}, which approaches
#' 1 as \eqn{T \to \infty}. For \eqn{h > 1} the correction inflates the statistic,
#' improving finite-sample size control. The corrected statistic \eqn{DM^*_k} is
#' compared to a \eqn{t(T-1)} distribution.
#' The analytic p-value (\code{P_Value}) uses the t-distribution with \eqn{P - 1}
#' degrees of freedom. The bootstrap p-value (\code{P_Value_Boot}) uses MBB
#' resampling (Kunsch, 1989) with recentering at the sample mean \eqn{\bar{d}_k},
#' placing the bootstrap distribution under H0. The Harvey correction is applied
#' consistently to both the analytic and bootstrap statistics. The p-values are
#' computed according to the alternative hypothesis specified by \code{H1}:
#' \subsection{Alternative Hypothesis and P-values (\code{H1})}{
#'   \describe{
#'     \item{\code{"same"} -- two-sided}{\eqn{p = 2 \cdot P(t_{T-1} < -|DM^*_k|)}}
#'     \item{\code{"more"} -- one-sided right}{\eqn{p = P(t_{T-1} > DM^*_k)}}
#'     \item{\code{"less"} -- one-sided left}{\eqn{p = P(t_{T-1} < DM^*_k)}}
#'   }
#'   The bootstrap analogue replaces the t-distribution probability with the
#'   empirical proportion of bootstrap statistics falling in the appropriate tail.
#' }
#' where \eqn{T} follows a t-distribution with \eqn{P - 1} degrees of freedom.
#' The bootstrap analogue replaces the t-distribution tail probability with the
#' empirical proportion of bootstrap statistics falling in the appropriate tail.
#'
#' This function performs \eqn{K} individual tests and does \emph{not} control for
#' multiple comparisons. For a joint test controlling the family-wise error rate,
#' use \code{\link{white_reality_check}} or
#' \code{\link{superior_predictive_ability_test}}.
#'
#' \strong{Note on MASE:} When \code{loss_differences} are constructed from
#' Mean Absolute Scaled Errors, the scaling (division by the naive benchmark
#' MAE, i.e. \code{mean(abs(diff(realizations)))}) must be applied
#' \emph{before} passing the loss differentials to this function.
#' \code{compute_per_model_statistics} receives pre-computed loss differentials
#' and applies no internal rescaling — the caller is responsible for ensuring
#' that MASE-based \code{loss_differences} already contain scaled errors.
#' In \code{\link{run_comprehensive_erc_analysis}} this is handled automatically.
#' @param loss_differences \code{\link[base]{matrix}} of dimension \code{P x K},
#'   where \code{P} is the number of forecast periods and \code{K} is the number of
#'   competing forecasts. Each column \eqn{k} contains the loss differential series
#'   \eqn{d_{t,k} = g(e_{0,t}) - g(e_{k,t})} for a generic loss function \eqn{g}.
#'   In the standard workflow of this package
#'   (\code{\link{run_comprehensive_erc_analysis}}), \eqn{g} is the squared error
#'   loss, so
#'   \eqn{d_{t,k} = (y_t - \hat{y}_{t,0})^2 - (y_t - \hat{y}_{t,k})^2}.
#'   A positive value of \eqn{d_{t,k}} means the forecast from model \eqn{k}
#'   is more accurate than the benchmark forecast in period \eqn{t}.
#' @param model_names \code{\link[base]{character}} vector of length \code{K} with
#'   names of the competing model forecasts.
#' @param block_length \code{\link[base]{integer}} block length for HAC and MBB.
#'   Rule of thumb: \eqn{T^{1/3}}; for \code{P = 165} approximately 5--6.
#'   Default is \code{5}.
#' @param n_boot \code{\link[base]{integer}} number of MBB replications for
#'   \code{P_Value_Boot}. Default \code{999}; see Davidson & MacKinnon (2000).
#' @param alpha \code{\link[base]{numeric}} significance level. Default \code{0.05}.
#' @param h \code{\link[base]{integer}} forecast horizon (number of steps ahead).
#'   Default is \code{1} (one-step-ahead). Passed to the Harvey, Leybourne &
#'   Newbold (1997) small-sample correction; see \emph{Details}.
#' @param H1 \code{\link[base]{character}} alternative hypothesis: \code{"same"}
#'   (two-sided, default), \code{"more"} (one-sided, forecast \eqn{k} better),
#'   or \code{"less"} (one-sided, forecast \eqn{k} worse). See \emph{Details}.
#' @return \code{\link[base]{data.frame}} with one row per competing model forecast:
#' \tabular{ll}{
#'   \code{Model}                      \tab Model name. \cr
#'   \code{Mean_Loss_Diff}             \tab Sample mean of \eqn{d_{t,k}}. \cr
#'   \code{Frac_Better_Than_Benchmark} \tab Fraction of periods where \eqn{d_{t,k} > 0}. \cr
#'   \code{T_Stat}                     \tab Harvey-corrected DM statistic \eqn{DM^*_k}. \cr
#'   \code{P_Value}                    \tab Analytic p-value (t-distribution, \eqn{T-1} df). \cr
#'   \code{P_Value_Boot}               \tab MBB bootstrap p-value. \cr
#'   \code{Significant}                \tab \code{TRUE} if \code{P_Value <= alpha}. \cr
#'   \code{Significant_Boot}           \tab \code{TRUE} if \code{P_Value_Boot <= alpha}. \cr
#' }
#' @seealso \code{\link{white_reality_check}},
#'   \code{\link{superior_predictive_ability_test}},
#'   \code{\link{estimate_long_run_covariance}}
#'
#' @references
#' Davidson, R., & MacKinnon, J. G. (2000).
#'   Bootstrap tests: How many bootstraps?
#'   \emph{Econometric Reviews}, 19(1), 55--68.
#'   \doi{10.1080/07474930008800459}
#'
#' Diebold, F. X., & Mariano, R. S. (1995). Comparing Predictive Accuracy.
#' \emph{Journal of Business & Economic Statistics}, 13(3), 253--263.
#' \doi{10.1080/07350015.1995.10524599}
#'
#' Harvey, D., Leybourne, S., & Newbold, P. (1997). Testing the equality of
#' prediction mean squared errors. \emph{International Journal of Forecasting},
#' \emph{13}(2), 281--291. \doi{10.1016/S0169-2070(96)00719-4}
#'
#' Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary
#' observations. \emph{The Annals of Statistics}, 17(3), 1217--1241.
#' \doi{10.1214/aos/1176347265}
#'
#' Newey, W. K., & West, K. D. (1987). A Simple, Positive Semi-Definite,
#' Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.
#' \emph{Econometrica}, 55(3), 703--708. \doi{10.2307/1913610}
#'
#' @examples
#' data(metals)
#' P       <- nrow(metals)
#' K_total <- ncol(metals)
#' K       <- K_total - 1
#' # A small offset (+0.5) avoids degenerate zero loss differences (illustration only).
#' realized         <- c(metals[-1, K_total], metals[P, K_total]) + 0.5
#' bench_loss       <- (metals[, K_total] - realized)^2
#' model_loss       <- (metals[, 1:K]     - realized)^2
#' loss_differences <- bench_loss - model_loss
#' model_names      <- colnames(metals)[1:K]
#' # Two-sided test (default)
#' result_df <- compute_per_model_statistics(loss_differences, model_names,
#'                                           n_boot = 10)
#' print(result_df)
#' # One-sided test: H1 = forecast is more accurate than benchmark
#' result_more <- compute_per_model_statistics(loss_differences, model_names,
#'                                             n_boot = 10, H1 = "more")
#' print(result_more)
#' @export
compute_per_model_statistics <- function(loss_differences, model_names,
                                         n_boot       = 999,
                                         block_length = 5,
                                         alpha        = 0.05,
                                         h            = 1,
                                         H1           = "same") {
  H1 <- match.arg(H1, choices = c("same", "more", "less"))
  P  <- nrow(loss_differences)
  K  <- ncol(loss_differences)
  
  per_model_results <- list()
  
  for (k in seq_len(K)) {
    model_name  <- model_names[k]
    loss_k      <- loss_differences[, k]
    mean_loss   <- mean(loss_k, na.rm = TRUE)
    frac_better <- mean(loss_k > 0, na.rm = TRUE)
    
    lrc     <- estimate_long_run_covariance(matrix(loss_k, ncol = 1),
                                            block_length = block_length)
    se_loss <- sqrt(as.numeric(lrc) / sum(!is.na(loss_k)))

    if (!is.na(se_loss) && se_loss > 1e-9) {
      T_eff  <- sum(!is.na(loss_k))
      DM     <- mean_loss / se_loss * sqrt((T_eff + 1 - 2 * h + (1 / T_eff) * h * (h - 1)) / T_eff)
      t_stat <- DM
      if (H1 == "same") { p_value <- 2 * min(pt(q = DM, df = T_eff - 1, lower.tail = FALSE),
                                             1 - pt(q = DM, df = T_eff - 1, lower.tail = FALSE)) }
      if (H1 == "less") { p_value <- pt(q = DM, df = T_eff - 1, lower.tail = FALSE) }
      if (H1 == "more") { p_value <- 1 - pt(q = DM, df = T_eff - 1, lower.tail = FALSE) }
    } else {
      t_stat  <- NA
      p_value <- NA
    }
    
    loss_k_mat   <- matrix(loss_k, ncol = 1)
    boot_t_stats <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      boot_sample   <- mbb_resample_data(loss_k_mat, block_length)[, 1]
      boot_centered <- boot_sample - mean_loss
      boot_lrc      <- estimate_long_run_covariance(
        matrix(boot_centered, ncol = 1),
        block_length = block_length)
      boot_se <- sqrt(as.numeric(boot_lrc) / sum(!is.na(boot_centered)))
      boot_t_stats[b] <- if (!is.na(boot_se) && boot_se > 1e-9) {
        T_boot        <- sum(!is.na(boot_centered))
        harvey_boot   <- sqrt((T_boot + 1 - 2 * h + (1 / T_boot) * h * (h - 1)) / T_boot)
        ((mean(boot_sample, na.rm = TRUE) - mean_loss) / boot_se) * harvey_boot
      } else {
        NA_real_
      }
    }
    
    p_value_boot <- if (!is.na(t_stat)) {
      if (H1 == "same") {
        mean(abs(boot_t_stats) >= abs(t_stat), na.rm = TRUE)
      } else if (H1 == "more") {
        mean(boot_t_stats >= t_stat, na.rm = TRUE)
      } else {
        mean(boot_t_stats <= t_stat, na.rm = TRUE)
      }
    } else {
      NA_real_
    }
    
    per_model_results[[model_name]] <- data.frame(
      Model                      = model_name,
      Mean_Loss_Diff             = round(mean_loss,   6),
      Frac_Better_Than_Benchmark = round(frac_better, 3),
      T_Stat                     = round(ifelse(is.na(t_stat),       NA, t_stat),       4),
      P_Value                    = round(ifelse(is.na(p_value),      NA, p_value),      4),
      P_Value_Boot               = round(ifelse(is.na(p_value_boot), NA, p_value_boot), 4),
      Significant                = ifelse(is.na(p_value),      NA, p_value      <= alpha),
      Significant_Boot           = ifelse(is.na(p_value_boot), NA, p_value_boot <= alpha),
      stringsAsFactors           = FALSE,
      row.names                  = NULL
    )
  }
  
  result_df           <- do.call(rbind, per_model_results)
  rownames(result_df) <- NULL
  return(result_df)
}