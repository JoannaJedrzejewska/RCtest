#' @importFrom stats dnorm pnorm pchisq pt quantile var sd
#' @importFrom utils head

#' @title Compute Continuous Ranked Probability Score (CRPS)
#'
#' @description Calculates the Continuous Ranked Probability Score (CRPS) using the
#' energy score (Monte Carlo) approximation for a single forecast period.
#'
#' @param forecast_density \code{\link[base]{numeric}} vector of simulated forecasts
#'   (density samples) representing the predictive distribution for a single time period.
#' @param target_realization \code{\link[base]{numeric}} scalar representing the realized
#'   value against which the forecast density is evaluated.
#'
#' @details
#' The CRPS is a strictly proper scoring rule that jointly rewards calibration and sharpness
#' of a probabilistic forecast. It is computed via the energy score identity:
#' \deqn{CRPS = E|X - y| - \frac{1}{2} E|X - X'|}
#' where \eqn{X, X'} are independent draws from the forecast distribution and \eqn{y} is
#' the realization. \strong{Lower values are better}: a CRPS of 0 indicates a perfect
#' point-mass forecast at the true realization.
#'
#' @return \code{\link[base]{numeric}} scalar representing the CRPS loss, or \code{NA} if
#'   input is invalid. Lower values indicate better probabilistic forecast accuracy.
#'
#' @references
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules, Prediction, and
#' Estimation. \emph{Journal of the American Statistical Association}, 102(477), 359--378.
#' \doi{10.1198/016214506000001437}
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#'
#' # CRPS for forecast 1, period 1:
#' # Use the cross-sectional spread of all competing forecasts at period t=1 as the density
#' density_samples <- as.numeric(metals[1, 1:14])
#' realized_value <- metals[1, 15]
#' compute_crps(density_samples, realized_value)
#'
#' # In practice, iterate over all forecasts and periods.
#' # For forecast k and period t, the predictive density is approximated by shifting the
#' # cross-sectional spread of all K competing forecasts so that it is centred at the
#' # cross-sectional mean of forecasts at period t. Specifically, for each forecast k:
#' #   density_samples_tk = (forecasts of all K forecast at t) - forecast_k(t) +
#'                        mean # (all K forecasts at t)
#' # This preserves the spread (diversity) across forecasts while recentring around the
#' # cross-sectional mean rather than around forecast k's own point forecast. It is an
#' # empirical approximation to the predictive distribution when no parametric density
#' # is available.
#' P <- nrow(metals)
#' K <- ncol(metals) - 1L # 14 competing forecasts
#' crps_matrix <- matrix(NA_real_, nrow = P, ncol = K,
#'                       dimnames = list(NULL, colnames(metals)[1:K]))
#' for (t in seq_len(P)) {
#'   for (k in seq_len(K)) {
#'     density_samples_tk <- as.numeric(metals[t, 1:K]) - metals[t, k] + mean(metals[t, 1:K])
#'     crps_matrix[t, k]  <- compute_crps(as.numeric(density_samples_tk),
#'                                        target_realization = metals[t, ncol(metals)])
#'   }
#' }
#' head(crps_matrix)
#' @export
compute_crps <- function(forecast_density, target_realization) {
  if (is.na(target_realization) || length(forecast_density) == 0 || all(is.na(forecast_density))) {
    return(NA)
  }
  clean_density <- forecast_density[is.finite(forecast_density) & !is.na(forecast_density)]
  n_sim <- length(clean_density)
  if (n_sim < 1) return(NA)
  f_sorted <- sort(clean_density)
  crps_val <- mean(abs(f_sorted - target_realization)) - 0.5 * mean(abs(outer(f_sorted, f_sorted, `-`)))
  return(crps_val)
}

#' @title Long-Run Covariance Estimator via Bartlett Kernel (HAC)
#'
#' @description Estimates the long-run covariance matrix using the Newey-West (1987) approach
#' with a Bartlett kernel. Provides Heteroskedasticity and Autocorrelation Consistent (HAC)
#' variance estimates used for studentizing Reality Check test statistics.
#'
#' @param loss_differences A \code{\link[base]{numeric}} matrix (\code{P x K}) of loss
#'   differences (benchmark loss minus forecast loss), where \code{P} is the number of
#'   forecast periods and \code{K} is the number of competing forecasts.
#' @param block_length \code{\link[base]{integer}}. The truncation lag \eqn{l} for the
#'   Bartlett kernel, numerically set equal to the MBB block length used elsewhere in
#'   this package for consistency. In HAC estimation this controls how many autocovariance
#'   lags are included; in MBB it controls block size -- both capture the same dependence
#'   horizon. A commonly used rule of thumb is \eqn{l \approx T^{1/3}}
#'   (Politis & Romano, 1994). For \code{P = 165} this gives approximately 5--6.
#'
#' @details
#' Implements the Newey-West (1987) HAC covariance matrix estimator with Bartlett kernel
#' weights \eqn{w_j = 1 - j / (l + 1)} for lags \eqn{j = 1, \ldots, l}, where \eqn{l}
#' denotes the truncation lag (following the notation of Newey & West, 1987, and
#' Politis & Romano, 1994), here set equal to \code{block_length}. This is essential
#' for accounting for serial dependence in time-series forecast evaluations.
#'
#' @return A symmetric positive semi-definite \code{\link[base]{matrix}} of dimensions
#'   \code{K x K} representing the estimated long-run covariance.
#'
#' @references
#' Newey, W. K., & West, K. D. (1987). A Simple Positive Semi-Definite Heteroskedasticity
#' and Autocorrelation Consistent Covariance Matrix. \emph{Econometrica}, 55(3), 703--708.
#' \doi{10.2307/1913610}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#' # A small offset (+0.5) is added to the lagged benchmark to avoid degenerate zero
#' # loss differences when forecasts equal the realized value exactly (illustration only).
#' P <- nrow(metals)
#' K_total <- ncol(metals)
#' K <- K_total - 1 # 14 competing forecasts
#' realized <- c(metals[-1, K_total], metals[P, K_total]) + 0.5
#' benchmark_loss <- (metals[, K_total] - realized)^2
#' model_loss     <- (metals[, 1:K] - realized)^2
#' loss_diff      <- benchmark_loss - model_loss
#' lrc_result <- estimate_long_run_covariance(loss_diff, block_length = 5)
#' print(round(lrc_result[1:3, 1:3], 6))
#' @export
estimate_long_run_covariance <- function(loss_differences, block_length) {
  T_obs    <- nrow(loss_differences)
  K_models <- ncol(loss_differences)
  data_centered <- loss_differences - matrix(colMeans(loss_differences, na.rm = TRUE),
                                             nrow = T_obs, ncol = K_models, byrow = TRUE)
  V_hat <- matrix(0, K_models, K_models)
  V_hat <- V_hat + crossprod(data_centered, data_centered) / T_obs
  
  for (j in 1:(T_obs - 1)) {
    if (j > block_length) break
    w       <- 1 - j / (block_length + 1)
    Gamma_j <- crossprod(data_centered[(j+1):T_obs, , drop = FALSE],
                         data_centered[1:(T_obs-j), , drop = FALSE]) / T_obs
    V_hat   <- V_hat + w * (Gamma_j + t(Gamma_j))
  }
  return(V_hat)
}

#' @title Moving Block Bootstrap (MBB) Resampler
#'
#' @description Generates a bootstrap resample of a time series matrix using the Moving
#' Block Bootstrap (MBB) method of Kunsch (1989).
#'
#' @param data_series \code{\link[base]{matrix}} where rows are observations (\code{P}) and
#'   columns are variables (\code{K}).
#' @param block_length \code{\link[base]{integer}} block length for the resampler. Overlapping
#'   blocks of this length are sampled with replacement. A commonly used rule of thumb is
#'   \code{block_length} \eqn{\approx T^{1/3}} (Politis & Romano, 1994). For \code{P = 165},
#'   this gives approximately 5--6.
#'
#' @details
#' Resamples overlapping blocks of \code{block_length} consecutive rows with replacement,
#' then concatenates them to produce a bootstrap sample of the same length \code{P} as
#' the original series. This preserves the short-run autocorrelation structure of the data,
#' which is required for valid inference in the Reality Check and SPA-type tests.
#'
#' @return \code{\link[base]{matrix}} of bootstrap resampled data with the same dimensions
#'   as \code{data_series}.
#'
#' @references
#' Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations.
#' \emph{The Annals of Statistics}, 17(3), 1217--1241. \doi{10.1214/aos/1176347265}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its recent
#' extensions. In \emph{Festschrift in honor of Halbert L. White}.
#' 
#' Liu, R. Y., & Singh, K. (1992). Moving blocks jackknife and bootstrap capture weak
#' dependence. In R. LePage & L. Billard (Eds.),
#' \emph{Exploring the Limits of Bootstrap} (pp. 225--248). Wiley.
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-3 are the first three competing forecasts
#' mbb_resample_data(metals[, 1:3], block_length = 5)
#' @export
mbb_resample_data <- function(data_series, block_length) {
  P          <- nrow(data_series)
  num_blocks <- ceiling(P / block_length)
  start_indices <- sample(1:(P - block_length + 1), num_blocks, replace = TRUE)
  boot_indices  <- unlist(sapply(start_indices, function(i) i:(i + block_length - 1)))
  boot_indices  <- boot_indices[1:P]
  return(data_series[boot_indices, , drop = FALSE])
}

#' @title White's Reality Check (WRC)
#'
#' @description Implements White's (2000) Reality Check (WRC) for comparing forecast
#' accuracy of multiple competing forecasts against a benchmark forecast based on mean 
#' loss differences. The test controls the family-wise error rate across all forecast 
#' comparisons simultaneously, avoiding data-snooping bias.
#'
#' \strong{Hypotheses:}
#' \itemize{
#'   \item \strong{H0:} \eqn{\max_{k} E[g(u_{0,t}) - g(u_{k,t})] \leq 0} -- no competing
#'     forecast produces a strictly lower expected loss than the benchmark forecast.
#'   \item \strong{H1:} At least one competing forecast has strictly lower expected loss
#'     than the benchmark forecast.
#' }
#'
#' @param loss_differences A \code{\link[base]{numeric}} matrix (\code{P x K}) of loss
#'   differences (benchmark loss minus forecast loss), where \code{P} is the number of
#'   forecast periods and \code{K} is the number of competing forecasts.
#'   A positive entry means the competing forecast outperforms the benchmark
#'   forecast in that period.
#' @param n_simulations \code{\link[base]{integer}}. The number of MBB bootstrap
#'   replications. Default \code{999}; see Davidson & MacKinnon (2000).
#' @param block_length \code{\link[base]{integer}}. The block length for the Moving Block
#'   Bootstrap (MBB). A commonly used rule of thumb is \eqn{T^{1/3}} (Politis & Romano,
#'   1994). For \code{P = 165}, this gives approximately 5--6.
#' @param alpha \code{\link[base]{numeric}}. The significance level
#'   (default \code{0.05}).
#'
#' @details
#' The test statistic is \eqn{\hat{S}_P = \max_k \overline{d}_k}, where
#' \eqn{\overline{d}_k} is the sample mean of the loss differential series for forecast
#' \eqn{k} (White, 2000, eq. 2). Bootstrap p-values are obtained via the MBB of
#' Kunsch (1989) by recentring each bootstrap statistic at the sample mean, following
#' the procedure in Corradi & Swanson (2011). This is an \emph{unstudentized} test; for
#' a studentized version with improved power against irrelevant forecasts, see
#' \code{\link{superior_predictive_ability_test}}.
#'
#' @return An object of class \code{"htest"}. The printed output shows the
#'   test statistic (maximum mean loss differential), the bootstrap p-value, and the
#'   test name. A small p-value (below \code{alpha}) leads to rejection of H0,
#'   indicating that at least one competing forecast is significantly more accurate
#'   than the benchmark forecast.
#'
#' @references
#' White, H. (2000). A reality check for data snooping. \emph{Econometrica}, 68(5),
#' 1097--1126. \doi{10.1111/1468-0262.00152}
#'
#' Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary
#' observations. \emph{The Annals of Statistics}, 17(3), 1217--1241.
#' \doi{10.1214/aos/1176347265}
#' 
#' Davidson, R., & MacKinnon, J. G. (2000). Bootstrap tests: How many bootstraps?
#'   \emph{Econometric Reviews}, 19(1), 55--68.
#'   \doi{10.1080/07474930008800459}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its recent
#' extensions. In \emph{Festschrift in honor of Halbert L. White}.
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#' # A small offset (+0.5) is added to the lagged benchmark to avoid degenerate zero
#' # loss differences when forecasts equal the realized value exactly (illustration only).
#' P <- nrow(metals)
#' K_total <- ncol(metals)
#' K <- K_total - 1 # 14 competing forecasts
#' realized       <- c(metals[-1, K_total], metals[P, K_total]) + 0.5
#' benchmark_loss <- (metals[, K_total] - realized)^2
#' model_loss     <- (metals[, 1:K] - realized)^2
#' loss_diff      <- benchmark_loss - model_loss
#' res <- white_reality_check(loss_diff, block_length = 5, n_simulations = 50)
#' print(res)
#' @export
white_reality_check <- function(loss_differences,n_simulations = 999,
                                block_length = 5, alpha = 0.05) {
  SP_k           <- colMeans(loss_differences, na.rm = TRUE)
  test_statistic <- max(SP_k, na.rm = TRUE)
  
  bootstrap_stats <- numeric(n_simulations)
  for (b in 1:n_simulations) {
    boot_sample        <- mbb_resample_data(loss_differences, block_length)
    boot_SP_k          <- colMeans(boot_sample, na.rm = TRUE)
    boot_SP_k_centered <- boot_SP_k - SP_k
    bootstrap_stats[b] <- max(boot_SP_k_centered, na.rm = TRUE)
  }
  p_value <- mean(bootstrap_stats > test_statistic, na.rm = TRUE)
  
  res <- list(
    statistic   = c("max mean loss diff" = test_statistic),
    p.value     = p_value,
    method      = "White's Reality Check (WRC)",
    data.name   = deparse(substitute(loss_differences)),
    null.value  = c("max mean loss differential" = 0),
    alternative = "at least one competing forecast outperforms the benchmark forecast",
    reject_null = p_value <= alpha
  )
  class(res) <- "htest"
  return(res)
}

#' @title Superior Predictive Ability (SPA) Test
#'
#' @description Implements the Hansen (2005) Superior Predictive Ability (SPA) test, a
#' studentized extension of White's (2000) Reality Check that corrects for the
#' inclusion of irrelevant (poor) forecasts to reduce conservatism.
#'
#' \strong{Hypotheses:}
#' \itemize{
#'   \item \strong{H0:} \eqn{\max_{k} E[g(u_{0,t}) - g(u_{k,t})] \leq 0} -- no competing
#'     forecast produces strictly lower expected loss than the benchmark forecast.
#'   \item \strong{H1:} At least one competing forecast has strictly lower expected loss
#'     than the benchmark forecast.
#' }
#'
#' @param loss_differences A \code{\link[base]{numeric}} matrix (\code{P x K}) of loss
#'   differences (benchmark loss minus forecast loss).
#' @param block_length \code{\link[base]{integer}}. The block length for MBB and HAC
#'   estimation. A commonly used rule of thumb is \eqn{T^{1/3}} (Politis & Romano,
#'   1994). For \code{P = 165}, this gives approximately 5--6.
#' @param num_bootstrap_replications \code{\link[base]{integer}} number of MBB
#'   bootstrap replications. Default \code{999}; use at least \code{999} for
#'   reliable inference (Davidson & MacKinnon, 2000).
#' @param alpha \code{\link[base]{numeric}}. The significance level
#'   (default \code{0.05}).
#' @details
#' The SPA statistic studentizes each mean loss differential by its HAC standard
#' deviation (estimated via \code{\link{estimate_long_run_covariance}}), then takes
#' the maximum across forecasts. Two p-values are returned, corresponding to two choices
#' of the null distribution (Hansen, 2005, Section 3):
#' \itemize{
#'   \item \code{p_consistent}: uses the sample-dependent null estimator
#'     \eqn{\hat{\mu}^c}, which recentres the bootstrap statistic at the sample mean
#'     \eqn{\bar{d}_k} for each forecast. This is the \strong{recommended} p-value.
#'   \item \code{p_conservative}: uses the Least Favourable Configuration (LFC)
#'     \eqn{\hat{\mu}^u = 0} for all forecasts -- equivalent to White's (2000) Reality
#'     Check bootstrap, where no recentring is applied. This provides an upper bound
#'     on the true p-value and is always \eqn{\geq} \code{p_consistent}.
#' }
#' @return An object of class \code{"htest"}. Additionally contains
#'   \code{p_consistent} and \code{p_conservative} for the two SPA bootstrap variants.
#' @references
#' Hansen, P. R. (2005). A Test for Superior Predictive Ability.
#' \emph{Journal of Business & Economic Statistics}, 23(4), 365--380.
#' \doi{10.1198/073500105000000063}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its recent
#' extensions. In \emph{Festschrift in honor of Halbert L. White}.
#' 
#' Davidson, R., & MacKinnon, J. G. (2000).
#' Bootstrap tests: How many bootstraps?
#' \emph{Econometric Reviews}, 19(1), 55--68.
#' \doi{10.1080/07474930008800459}
#' 
#' Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary
#' observations. \emph{The Annals of Statistics}, 17(3), 1217--1241.
#' \doi{10.1214/aos/1176347265}
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#' # A small offset (+0.5) is added to the lagged benchmark to avoid degenerate zero
#' # loss differences when forecasts equal the realized value exactly (illustration only).
#' P <- nrow(metals)
#' K_total <- ncol(metals)
#' K <- K_total - 1 # 14 competing forecasts
#' realized       <- c(metals[-1, K_total], metals[P, K_total]) + 0.5
#' benchmark_loss <- (metals[, K_total] - realized)^2
#' model_loss     <- (metals[, 1:K] - realized)^2
#' loss_diff      <- benchmark_loss - model_loss
#' res <- superior_predictive_ability_test(loss_diff, block_length = 5,
#'                                         num_bootstrap_replications = 50,
#'                                         alpha = 0.05)
#' print(res)
#' @export
superior_predictive_ability_test <- function(loss_differences, block_length,
                                             num_bootstrap_replications, alpha) {
  P    <- nrow(loss_differences)
  K    <- ncol(loss_differences)
  SP_k <- colMeans(loss_differences, na.rm = TRUE)
  V_hat_full <- estimate_long_run_covariance(loss_differences, block_length)
  V_k        <- diag(V_hat_full)
  V_k[V_k <= 1e-10] <- 1e-10
  std_dev_k  <- sqrt(V_k)
  
  T_k_unscaled     <- SP_k / std_dev_k
  T_SPA_Consistent <- max(T_k_unscaled)
  
  bootstrap_consistent_all_k   <- matrix(NA, nrow = num_bootstrap_replications, ncol = K)
  bootstrap_conservative_all_k <- matrix(NA, nrow = num_bootstrap_replications, ncol = K)
  
  for (b in 1:num_bootstrap_replications) {
    boot_sample <- mbb_resample_data(loss_differences, block_length)
    boot_SP_k   <- colMeans(boot_sample, na.rm = TRUE)
    boot_T_k    <- boot_SP_k / std_dev_k
    bootstrap_consistent_all_k[b, ]   <- boot_T_k - T_k_unscaled
    bootstrap_conservative_all_k[b, ] <- boot_T_k
  }
  
  T_C_boot_max         <- apply(bootstrap_consistent_all_k,   1, max, na.rm = TRUE)
  T_R_boot_max         <- apply(bootstrap_conservative_all_k, 1, max, na.rm = TRUE)
  p_value_consistent   <- mean(T_C_boot_max > T_SPA_Consistent, na.rm = TRUE)
  p_value_conservative <- mean(T_R_boot_max > T_SPA_Consistent, na.rm = TRUE)

  res <- list(
    statistic      = c("T-SPA" = T_SPA_Consistent),
    p.value        = p_value_consistent,
    method         = "Superior Predictive Ability (SPA) Test",
    data.name      = deparse(substitute(loss_differences)),
    alternative    = "at least one competing forecast outperforms the benchmark forecast",
    p_consistent   = p_value_consistent,
    p_conservative = p_value_conservative,
    reject_null    = p_value_consistent <= alpha
  )
  class(res) <- "htest"
  return(res)
}

#' @title Conditional Predictive Ability (CPA) Reality Check Test
#'
#' @description Implements the Conditional Predictive Ability (CPA) test of Giacomini &
#' White (2006), extended to a multiple-forecast setting via a studentized Reality Check
#' statistic. Tests whether any competing forecast's predictive advantage over the benchmark
#' is state-dependent, i.e., predictable from a conditioning variable \eqn{h_t} known
#' at the time the forecast is made.
#'
#' \strong{Hypotheses:}
#' \itemize{
#'   \item \strong{H0:} \eqn{E[h_t \cdot (g(u_{0,t}) - g(u_{k,t}))] = 0} for all
#'     \eqn{k = 1,\ldots,K} -- no competing forecast's loss differential with the benchmark
#'     is predictable using the conditioning information \eqn{h_t}.
#'   \item \strong{H1:} At least one forecast's loss differential
#'     \eqn{d_{k,t} = g(u_{0,t}) - g(u_{k,t})} is predictable by \eqn{h_t}, i.e.,
#'     \eqn{E[h_t \cdot d_{k,t}] \neq 0} for some \eqn{k}.
#' }
#'
#' @param loss_differences A \code{\link[base]{numeric}} matrix (\code{P x K}) of loss
#'   differences (benchmark loss minus forecast loss), where \code{P} is the number of
#'   forecast periods and \code{K} is the number of competing forecasts.
#'   A positive entry means the competing forecast outperforms the benchmark
#'   forecast in that period.
#' @param weighting_vector \code{\link[base]{numeric}} vector of length \code{P}
#'   serving as the conditioning instrument \eqn{h_t} in the CPA test. At each
#'   period \eqn{t}, the test checks whether the loss differential \eqn{d_{k,t}}
#'   covaries with \eqn{h_t}, i.e., whether \eqn{E[h_t \cdot d_{k,t}] \neq 0}.
#'   See the \emph{Conditioning Instrument} section in Details for interpretation,
#'   requirements, and recommended choices.
#' @param block_length \code{\link[base]{integer}}. The block length for MBB and HAC
#'   estimation. A commonly used rule of thumb is \eqn{T^{1/3}} (Politis & Romano,
#'   1994). For \code{P = 165}, this gives approximately 5--6.
#' @param num_bootstrap_replications \code{\link[base]{integer}} number of MBB
#'   bootstrap replications. Default \code{999}; see Davidson & MacKinnon (2000).
#' @param alpha \code{\link[base]{numeric}}. The significance level
#'   (default \code{0.05}).
#' @details
#' The test multiplies each column of \code{loss_differences} element-wise by
#' \code{weighting_vector} to form the weighted loss differential series
#' \eqn{h_t \cdot d_{k,t}}. The unconditional mean of this product,
#' \eqn{E[h_t \cdot d_{k,t}]}, equals zero under H0 by the law of iterated
#' expectations when \eqn{h_t} is a valid instrument. The test statistic is the
#' maximum studentized mean across all \eqn{K} forecasts:
#' \deqn{\hat{T}_{CPA} = \max_{k} \frac{\frac{1}{P}\sum_t h_t d_{k,t}}
#'   {\hat{\sigma}_{k,h}}}
#' where \eqn{\hat{\sigma}_{k,h}} is the HAC standard deviation of \eqn{h_t d_{k,t}}
#' estimated via \code{\link{estimate_long_run_covariance}}. Bootstrap p-values are
#' obtained via the MBB of Kunsch (1989) with recentring, following the SPA-type
#' procedure of Hansen (2005) applied to the weighted series.
#' \subsection{Conditioning Instrument (\code{weighting_vector})}{
#'
#'   A significant result means that knowing \eqn{h_t} allows one to predict which
#'   forecast will perform better in period \eqn{t} -- the benchmark's advantage (or
#'   disadvantage) is state-dependent and potentially exploitable. This is a strictly
#'   stronger statement than the unconditional WRC: a forecast can fail the WRC (no
#'   unconditional improvement) yet pass the CPA test (conditional improvement in
#'   specific states).
#'
#'   \eqn{h_t} must be measurable with respect to the information set available at
#'   time \eqn{t} (Giacomini & White, 2006, Assumption 1) -- it must not use
#'   information from period \eqn{t+1} or later. The scale of \eqn{h_t} does not
#'   affect the test result because the statistic is studentized by its own HAC
#'   standard deviation.
#'
#'   Recommended choices:
#'   \describe{
#'     \item{\code{abs(realized)} -- absolute realised values}{
#'       Tests whether forecast performance depends on outcome magnitude -- a natural
#'       proxy for market volatility or economic uncertainty. Default in
#'       \code{\link{run_comprehensive_erc_analysis}}.}
#'     \item{\code{c(realized[1], realized[-length(realized)])} -- lagged realised values}{
#'       Tests whether the previous period's outcome predicts which forecast wins next
#'       period. Relevant when forecast errors are autocorrelated.}
#'     \item{\code{rep(1, P)} -- constant vector}{
#'       The product \eqn{h_t \cdot d_{k,t}} reduces to \eqn{d_{k,t}}, making the
#'       CPA test equivalent to the unconditional WRC. Use as a sanity check: results
#'       should be consistent with \code{\link{white_reality_check}}.}
#'     \item{External economic indicator}{
#'       E.g., a recession dummy, VIX level, lagged interest rate spread, or monetary
#'       policy stance dummy. Tests whether one forecast systematically outperforms
#'       during specific regimes. Must be lagged one period to ensure \eqn{h_t} is
#'       in the information set at the time of the forecast.}
#'   }
#' }
#' The scale of \code{weighting_vector} has no effect on inference because both the
#' test statistic and its bootstrap distribution are studentized by the same
#' \eqn{\hat{\sigma}_{k,h}}.
#'
#' @return An object of class \code{"htest"} with the following components:
#' \tabular{ll}{
#'   \code{statistic}   \tab Maximum studentized weighted mean loss differential
#'                           across all \eqn{K} forecasts, labelled \code{"T-CPA"}. \cr
#'   \code{p.value}     \tab Bootstrap p-value from the MBB procedure. \cr
#'   \code{method}      \tab \code{"Conditional Predictive Ability (CPA) Test"}. \cr
#'   \code{null.value}  \tab Named scalar \code{"max studentized weighted mean
#'                           loss differential"} = 0. \cr
#'   \code{alternative} \tab Direction of the alternative hypothesis. \cr
#'   \code{reject_null} \tab Logical: \code{TRUE} if \code{p.value <= alpha}. \cr
#' }
#' A small p-value indicates that at least one forecast's loss differential is
#' predictable from the conditioning variable \eqn{h_t}. Failure to reject H0
#' means no evidence of state-dependent predictive ability for the chosen instrument.
#'
#' @references
#' Giacomini, R., & White, H. (2006). Tests of Conditional Predictive Ability.
#' \emph{Econometrica}, 74(6), 1545--1578. \doi{10.1111/j.1468-0262.2006.00718.x}
#'
#' Hansen, P. R. (2005). A Test for Superior Predictive Ability.
#' \emph{Journal of Business & Economic Statistics}, 23(4), 365--380.
#' \doi{10.1198/073500105000000063}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its
#' recent extensions. In \emph{Festschrift in honor of Halbert L. White}.
#' 
#' Davidson, R., & MacKinnon, J. G. (2000).
#'   Bootstrap tests: How many bootstraps?
#'   \emph{Econometric Reviews}, 19(1), 55--68.
#'   \doi{10.1080/07474930008800459}
#'
#' Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary
#' observations. \emph{The Annals of Statistics}, 17(3), 1217--1241.
#' \doi{10.1214/aos/1176347265}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' @seealso
#' \code{\link{white_reality_check}} for the unconditional WRC test (equivalent to
#' CPA with a constant \code{weighting_vector});
#' \code{\link{superior_predictive_ability_test}} for the studentized unconditional test.
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#' # A small offset (+0.5) is added to the lagged benchmark to avoid degenerate zero
#' # loss differences when forecasts equal the realized value exactly (illustration only).
#' P       <- nrow(metals)
#' K_total <- ncol(metals)
#' K       <- K_total - 1L  # 14 competing forecasts
#' realized       <- c(metals[-1, K_total], metals[P, K_total]) + 0.5
#' benchmark_loss <- (metals[, K_total] - realized)^2
#' model_loss     <- (metals[, 1:K]     - realized)^2
#' loss_diff      <- benchmark_loss - model_loss
#'
#' # Example 1: absolute realised values as conditioning variable (volatility proxy)
#' res1 <- white_reality_check_conditional(
#'   loss_differences           = loss_diff,
#'   weighting_vector           = abs(realized),
#'   block_length               = 5,
#'   num_bootstrap_replications = 50,
#'   alpha                      = 0.05
#' )
#' print(res1)
#'
#' # Example 2: constant vector - should give results consistent with white_reality_check()
#' res2 <- white_reality_check_conditional(
#'   loss_differences           = loss_diff,
#'   weighting_vector           = rep(1, P),
#'   block_length               = 5,
#'   num_bootstrap_replications = 50,
#'   alpha                      = 0.05
#' )
#' print(res2)
#' @export
white_reality_check_conditional <- function(loss_differences, weighting_vector,
                                            block_length, num_bootstrap_replications,
                                            alpha) {
  P <- nrow(loss_differences)
  if (length(weighting_vector) != P)
    stop("Weighting vector length must match the number of forecast periods P.")
  
  weighted_diff_series <- loss_differences * weighting_vector
  V_hat_full_h <- estimate_long_run_covariance(weighted_diff_series, block_length)
  V_k_h        <- diag(V_hat_full_h)
  V_k_h[V_k_h <= 1e-10] <- 1e-10
  std_dev_k_h  <- sqrt(V_k_h)
  
  SP_k_h          <- colMeans(weighted_diff_series, na.rm = TRUE)
  T_k_unscaled_h  <- SP_k_h / std_dev_k_h
  T_CPA_statistic <- max(T_k_unscaled_h, na.rm = TRUE)
  
  bootstrap_stats_all_k <- matrix(NA, nrow = num_bootstrap_replications,
                                  ncol = ncol(loss_differences))
  for (b in 1:num_bootstrap_replications) {
    boot_weighted_diff    <- mbb_resample_data(weighted_diff_series, block_length)
    boot_SP_k_h           <- colMeans(boot_weighted_diff, na.rm = TRUE)
    boot_T_k_unscaled_h   <- boot_SP_k_h / std_dev_k_h
    boot_T_k_h_centered   <- boot_T_k_unscaled_h - T_k_unscaled_h
    bootstrap_stats_all_k[b, ] <- boot_T_k_h_centered
  }
  
  T_CPA_Consistent_boot_max <- apply(bootstrap_stats_all_k, 1, max, na.rm = TRUE)
  p_value_consistent <- mean(T_CPA_Consistent_boot_max > T_CPA_statistic, na.rm = TRUE)
  
  res <- list(
    statistic   = c("T-CPA" = T_CPA_statistic),
    p.value     = p_value_consistent,
    method      = "Conditional Predictive Ability (CPA) Test",
    data.name   = deparse(substitute(loss_differences)),
    null.value  = c("max studentized weighted mean loss differential" = 0),
    alternative = "at least one forecast's loss differential is predictable by the conditioning variable",
    reject_null = p_value_consistent <= alpha
  )
  class(res) <- "htest"
  return(res)
}

#' @title White's Reality Check via Expected Loss CDF Comparison (CDF-RC)
#'
#' @description Implements a studentized Reality Check test that compares competing
#' \strong{forecasts} against a benchmark across the entire distribution of loss
#' differences, not only the mean. For each forecast \eqn{k} and each quantile
#' threshold \eqn{x_\tau} (derived from the pooled empirical distribution of loss
#' differences), the test evaluates whether the empirical CDF of \strong{forecast
#' \eqn{k}'s} loss differences, evaluated at \eqn{x_{\tau_j}}, exceeds the nominal
#' quantile level \eqn{\tau_j} itself -- i.e., whether the competing forecast falls
#' below that threshold \emph{more often than expected} under a correctly calibrated
#' loss distribution, indicating stochastic dominance of the competing forecast's
#' loss distribution over the benchmark's.
#'
#' \strong{Hypotheses:}
#' \itemize{
#'   \item \strong{H0:} \eqn{\max_{k,j}\big(E[\mathbf{1}(d_{k,t} \leq x_{\tau_j})] - \tau_j\big) \leq 0}
#'     for all \eqn{k = 1,\ldots,K} and all quantile thresholds
#'     \eqn{x_{\tau_j},\ j = 1,\ldots,J} -- no competing forecast's empirical CDF of
#'     loss differences exceeds its nominal quantile level \eqn{\tau_j} at any
#'     evaluation point, i.e., no competing forecast is stochastically dominant
#'     over the benchmark anywhere in the loss distribution.
#'   \item \strong{H1:} At least one competing forecast's empirical CDF of loss
#'     differences significantly exceeds its nominal quantile level \eqn{\tau_j} at
#'     some threshold \eqn{x_{\tau_j}}, i.e., the benchmark is stochastically
#'     dominated in terms of loss differences at that point.
#' }
#' Note that because \eqn{x_{\tau_j}} is itself defined as the \eqn{\tau_j}-quantile
#' of the \emph{pooled} loss-difference distribution (across all \eqn{K} forecasts),
#' each column's empirical CDF value \eqn{\bar G_{k,j}} is centred at \eqn{\tau_j}
#' under a null of no systematic difference across forecasts -- \strong{not at
#' zero}. The test statistic and its bootstrap null distribution are therefore
#' constructed relative to \eqn{\tau_j}, not zero; see \dQuote{Details}.
#'
#' @param loss_differences A \code{\link[base]{numeric}} matrix (\code{P x K}) of loss
#'   differences (benchmark loss minus forecast forecast loss), where \code{P} is the number of
#'   forecast periods and \code{K} is the number of competing forecasts. A positive entry
#'   means the competing forecast is more accurate than the benchmark forecast in that period.
#' @param block_length \code{\link[base]{integer}}. The block length for the Moving Block
#'   Bootstrap (MBB) and for HAC variance estimation via
#'   \code{\link{estimate_long_run_covariance}}. A commonly used rule of thumb is
#'   \eqn{T^{1/3}} (Politis & Romano, 1994). For \code{P = 165}, this gives
#'   approximately 5--6.
#' @param num_bootstrap_replications \code{\link[base]{integer}} number of MBB
#'   bootstrap replications. Default \code{999}; see Davidson & MacKinnon (2000).
#' @param alpha \code{\link[base]{numeric}}. The significance level (default \code{0.05}).
#'
#' @details
#' The test proceeds in four steps.
#'
#' \strong{Step 1 -- Quantile grid.} A grid of \eqn{J = 9} evaluation points
#' \eqn{x_{\tau_1}, \ldots, x_{\tau_9}} is constructed as the
#' \eqn{\tau_j \in \{0.1, 0.2, \ldots, 0.9\}} quantiles of the \emph{pooled}
#' empirical distribution of all loss differences (across all forecasts and all periods).
#' Using quantiles of the data rather than a fixed grid ensures that the evaluation
#' points are always in the support of the observed loss differences.
#'
#' \strong{Step 2 -- Indicator matrix.} For each forecast \eqn{k} and each threshold
#' \eqn{x_{\tau_j}}, the binary indicator
#' \deqn{G_{k,j,t} = \mathbf{1}(d_{k,t} \leq x_{\tau_j})}
#' is formed, where \eqn{d_{k,t}} is the loss difference for forecast \eqn{k} at
#' period \eqn{t}. This yields a \eqn{P \times (K \cdot J)} matrix \code{Gdata}
#' with \eqn{K \times J = 14 \times 9 = 126} columns (for the \code{metals} dataset).
#' The column mean \eqn{\bar{G}_{k,j} = \frac{1}{P}\sum_t G_{k,j,t}} estimates the
#' empirical CDF of forecast \eqn{k}'s loss differences evaluated at \eqn{x_{\tau_j}}.
#'
#' \strong{Step 3 -- Null-centring (\emph{critical for correct inference}).} Because
#' \eqn{x_{\tau_j}} is defined as the \eqn{\tau_j}-quantile of the \emph{pooled}
#' loss-difference sample, \eqn{\bar{G}_{k,j}} is mechanically close to \eqn{\tau_j}
#' for every forecast \eqn{k} whenever forecasts are exchangeable under the null --
#' it is \strong{not} centred at zero. The relevant test quantity is therefore the
#' \emph{excess} empirical CDF over its nominal level,
#' \deqn{\bar{G}_{k,j} - \tau_j,}
#' which is genuinely centred at zero when no forecast is stochastically dominant.
#' Studentizing the raw \eqn{\bar{G}_{k,j}} (without subtracting \eqn{\tau_j}) before
#' comparing to a bootstrap distribution of \emph{deviations} from the observed
#' statistic produces a test statistic and a null distribution on incompatible
#' scales, and causes the test to reject in virtually all samples regardless of the
#' data -- this was corrected in the current package version (see \dQuote{Note}).
#'
#' \strong{Step 4 -- Studentized test statistic.} Each null-centred column mean is
#' studentized by its HAC standard deviation (from
#' \code{\link{estimate_long_run_covariance}}), and the test statistic is the
#' maximum studentized excess CDF value across all \eqn{K \times J} columns:
#' \deqn{\hat{T} = \max_{k,j} \frac{\bar{G}_{k,j} - \tau_j}{\hat{\sigma}_{k,j}}}
#' Bootstrap p-values are obtained via the MBB of Kunsch (1989) with recentring,
#' following the WRC procedure of White (2000) and Corradi & Swanson (2011); the
#' same \eqn{\tau_j} offset is subtracted from each bootstrap replicate's column
#' mean before studentizing, so that both the observed statistic and its bootstrap
#' distribution are computed on the same, correctly null-centred scale.
#'
#' \strong{Relationship to Corradi & Swanson (2006).} This test is a loss-difference
#' analogue of the predictive CDF comparison in Corradi & Swanson (2006, Section 4).
#' Rather than comparing forecast CDFs against the true conditional distribution
#' (as in the ZP test), it compares empirical CDFs of \emph{loss differences} against
#' their own nominal quantile levels, assessing stochastic dominance of a competing
#' forecast over the benchmark in terms of loss. It complements
#' \code{\link{white_reality_check}} (which tests only the mean) by detecting cases
#' where one forecast is better in the tails but not on average.
#'
#' \strong{Lower p-values are more informative:} rejection of H0 indicates that at
#' least one forecast stochastically dominates the benchmark at some point of the loss
#' distribution. Failure to reject does not preclude dominance at specific quantiles
#' -- it means no single \eqn{(k,j)} combination is significant after controlling
#' for multiple comparisons.
#'
#' @note
#' \strong{Package versions prior to the current release contained two
#' implementation errors in this function, both now corrected:}
#' \enumerate{
#'   \item \strong{Redundant variance division.}
#'     \code{\link{estimate_long_run_covariance}} already divides its output by the
#'     sample size internally; this function additionally divided
#'     \code{diag(V_hat_full)} by \code{P_clean} a second time before studentizing,
#'     inflating every studentized statistic by a factor of approximately
#'     \eqn{\sqrt{P}}. This has been removed.
#'   \item \strong{Missing null-centring.} As described in Step 3 above, the test
#'     statistic and bootstrap null distribution were computed on incompatible
#'     scales (raw empirical CDF values compared against a distribution of
#'     deviations from those same values), causing the test to reject in
#'     essentially 100\% of samples regardless of the data -- verified via an
#'     isolated Monte Carlo check under a pure i.i.d. null with a single forecast
#'     and a single quantile threshold. This function now subtracts the nominal
#'     quantile level \eqn{\tau_j} from each column's empirical CDF before
#'     studentizing, restoring rejection rates near the nominal significance level
#'     under simulated null data.
#' }
#' Results obtained from this function in prior package versions should be treated
#' as unreliable and are superseded by output from the current version.
#'
#' @return An object of class \code{"htest"} with the following components:
#' \tabular{ll}{
#'   \code{statistic}   \tab Maximum studentized excess CDF value (empirical CDF
#'                           minus nominal quantile level) across all
#'                           \eqn{K \times J} forecast-quantile combinations,
#'                           labelled \code{"KS-type"}. \cr
#'   \code{p.value}     \tab Bootstrap p-value from the MBB procedure. \cr
#'   \code{method}      \tab \code{"Expected Loss CDF Comparison Test"}. \cr
#'   \code{null.value}  \tab Named scalar
#'                           \code{"max studentized excess CDF value (empirical
#'                           CDF minus nominal quantile level)"} = 0. \cr
#'   \code{alternative} \tab Description of the alternative hypothesis. \cr
#' }
#' A small p-value (below \code{alpha}) leads to rejection of H0, indicating that
#' at least one competing forecast's empirical CDF of loss differences significantly
#' exceeds its nominal quantile level at some threshold -- i.e., the competing
#' forecast more frequently produces smaller losses than the benchmark forecast in
#' some region of the loss distribution than would be expected by chance.
#'
#' @references
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional
#' confidence interval accuracy tests. \emph{Journal of Econometrics},
#' 135(1--2), 187--228. \doi{10.1016/j.jeconom.2005.07.026}
#'
#' Davidson, R., & MacKinnon, J. G. (2000).
#' Bootstrap tests: How many bootstraps?
#' \emph{Econometric Reviews}, 19(1), 55--68.
#' \doi{10.1080/07474930008800459}
#'
#' White, H. (2000). A reality check for data snooping. \emph{Econometrica},
#' 68(5), 1097--1126. \doi{10.1111/1468-0262.00152}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its
#' recent extensions. In \emph{Festschrift in honor of Halbert L. White}.
#'
#' Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary
#' observations. \emph{The Annals of Statistics}, 17(3), 1217--1241.
#' \doi{10.1214/aos/1176347265}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' @seealso
#' \code{\link{white_reality_check}} for the mean-based WRC test;
#' \code{\link{reality_check_zp_test}} for a distributional test based on the
#' true conditional CDF; \code{\link{superior_predictive_ability_test}} for the
#' studentized mean-based SPA test.
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#' # A small offset (+0.5) is added to the lagged benchmark to avoid degenerate zero
#' # loss differences when forecasts equal the realized value exactly (illustration only).
#' P       <- nrow(metals)
#' K_total <- ncol(metals)
#' K       <- K_total - 1L  # 14 competing forecasts
#' realized       <- c(metals[-1, K_total], metals[P, K_total]) + 0.5
#' benchmark_loss <- (metals[, K_total] - realized)^2
#' model_loss     <- (metals[, 1:K]     - realized)^2
#' loss_diff      <- benchmark_loss - model_loss
#' res <- white_reality_check_cdf_approx(loss_diff,
#'                                       block_length               = 5,
#'                                       num_bootstrap_replications = 50)
#' print(res)
#' @export
white_reality_check_cdf_approx <- function(loss_differences, block_length,
                                           num_bootstrap_replications, alpha = 0.05) {
  P <- nrow(loss_differences)
  K <- ncol(loss_differences)
  
  quantile_grid <- seq(0.1, 0.9, by = 0.1)
  pooled_loss   <- as.vector(loss_differences)
  pooled_loss   <- pooled_loss[!is.na(pooled_loss)]
  x_tau_points  <- quantile(pooled_loss, probs = quantile_grid, na.rm = TRUE)
  J             <- length(x_tau_points)
  
  G_data      <- matrix(NA, nrow = P, ncol = K * J)
  null_target <- numeric(K * J)
  for (k in 1:K) {
    for (j in 1:J) {
      G_index             <- (k - 1) * J + j
      G_data[, G_index]   <- as.numeric(loss_differences[, k] <= x_tau_points[j])
      null_target[G_index] <- quantile_grid[j]
    }
  }
  
  na_rows <- apply(loss_differences, 1, anyNA)
  G_data  <- G_data[!na_rows, , drop = FALSE]
  P_clean <- nrow(G_data)
  if (P_clean == 0) return(NULL)
  
  S_mean          <- colMeans(G_data, na.rm = TRUE)
  S_mean_centered <- S_mean - null_target
  
  V_hat_full <- estimate_long_run_covariance(G_data, block_length)
  V_k_new    <- diag(V_hat_full)
  V_k_new[V_k_new <= 1e-10] <- 1e-10
  std_dev_k_new <- sqrt(V_k_new)
  
  T_k_unscaled <- S_mean_centered / std_dev_k_new
  T_max_ks     <- max(T_k_unscaled)
  
  bootstrap_stats <- numeric(num_bootstrap_replications)
  for (b in 1:num_bootstrap_replications) {
    boot_sample          <- mbb_resample_data(G_data, block_length)
    boot_S_mean          <- colMeans(boot_sample, na.rm = TRUE)
    boot_S_mean_centered <- boot_S_mean - null_target
    boot_T_unscaled      <- boot_S_mean_centered / std_dev_k_new
    boot_T_centered       <- boot_T_unscaled - T_k_unscaled
    bootstrap_stats[b]    <- max(boot_T_centered, na.rm = TRUE)
  }
  
  p_value <- mean(bootstrap_stats > T_max_ks, na.rm = TRUE)
  
  res <- list(
    statistic   = c("KS-type" = T_max_ks),
    p.value     = p_value,
    method      = "Expected Loss CDF Comparison Test",
    data.name   = deparse(substitute(loss_differences)),
    null.value  = c("max studentized excess CDF value (empirical CDF minus nominal quantile level)" = 0),
    alternative = "at least one competing forecast's empirical CDF of loss differences exceeds its nominal quantile level at some threshold"
  )
  class(res) <- "htest"
  return(res)
}

#' @title Kullback-Leibler Information Criterion (KLIC) Test
#'
#' @description Implements the Reality Check using Negative Log-Likelihood Scores (NLS) to
#' evaluate predictive densities in terms of their Kullback-Leibler divergence from the
#' true density. Based on Corradi & Swanson (2006).
#'\itemize{
#' \item \strong{H0:} \eqn{\max_k E[\log f_1(y_t) - \log f_k(y_t)] \leq 0} -- no
#'   competing forecast achieves a higher average log-likelihood (lower KLIC distance)
#'   than the benchmark density \eqn{f_1}.
#' \item \strong{H1:} At least one competing forecast achieves strictly higher average
#'   log-likelihood than the benchmark.
#'}
#' @param log_likelihood_differences A \code{\link[base]{numeric}} matrix (\code{P x K})
#'   of Negative Log-Likelihood Score (NLS) differences: benchmark NLS minus forecast NLS.
#'   A positive entry means the forecast's density assigns higher probability to the
#'   observed outcome than the benchmark density does.
#' @param block_length \code{\link[base]{integer}}. The block length for MBB and HAC
#'   estimation. A commonly used rule of thumb is
#'   \code{block_length} \eqn{\approx T^{1/3}}, where \eqn{T} is the number of
#'   observations (Politis & Romano, 1994). For \code{P = 165}, this gives
#'   approximately 5--6.
#' @param num_bootstrap_replications \code{\link[base]{integer}} number of MBB
#'   bootstrap replications. Default \code{999}; see Davidson & MacKinnon (2000).
#' @param alpha \code{\link[base]{numeric}}. The significance level (default \code{0.05}).
#'
#' @details
#' The KLIC between the true density \eqn{f_0} and a forecast density \eqn{f_k} is
#' \eqn{E[\log f_0(y) - \log f_k(y)]}. Minimising KLIC is equivalent to maximising
#' expected log-likelihood. This test therefore selects the forecast with the smallest
#' KLIC distance from the true density. \strong{Lower NLS values are better}: a forecast
#' with lower NLS assigns higher average probability to events that actually occurred
#' (Corradi & Swanson, 2006).
#'
#' The NLS loss matrix is constructed via \code{\link{compute_klic}} assuming normal
#' predictive densities parameterised by a point forecast and a rolling-window standard
#' deviation. The test statistic is the maximum studentized mean NLS differential,
#' with p-values obtained via MBB following the SPA bootstrap of Hansen (2005).
#'
#' @return An object of class \code{"htest"}. A small p-value (below \code{alpha})
#'   leads to rejection of H0, indicating that at least one competing forecast has a
#'   lower Kullback-Leibler distance from the true density than the benchmark.
#'   Failure to reject H0 suggests no forecast provides significantly better density
#'   fit than the benchmark.
#' @references
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional confidence
#' interval accuracy tests. \emph{Journal of Econometrics}, 135(1--2), 187--228.
#' \doi{10.1016/j.jeconom.2005.07.026}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its recent
#' extensions. In \emph{Festschrift in honor of Halbert L. White}.
#' 
#' Davidson, R., & MacKinnon, J. G. (2000). Bootstrap tests: How many bootstraps?
#' \emph{Econometric Reviews}, 19(1), 55--68.
#' \doi{10.1080/07474930008800459}
#' 
#' Hansen, P. R. (2005). A Test for Superior Predictive Ability.
#' \emph{Journal of Business & Economic Statistics}, 23(4), 365--380.
#' \doi{10.1198/073500105000000063}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#' benchmark_col      <- 15
#' P                  <- nrow(metals)
#' K_total            <- ncol(metals)
#' K                  <- K_total - 1 # 14 competing forecasts
#' realized           <- metals[, benchmark_col]
#' forecast_variance  <- estimate_forecast_variance(metals, realized = realized,
#'                         benchmark_col = benchmark_col, window_size = 20)
#' comp_cols          <- setdiff(seq_len(K_total), benchmark_col)
#' forecast_sd_models <- sqrt(forecast_variance[, comp_cols])
#' nls_loss  <- compute_klic(metals, forecast_sd_models, benchmark_col = K_total)
#' nls_diff  <- nls_loss[, K_total] - nls_loss[, comp_cols]
#' kullback_leibler_test(nls_diff, block_length = 5, num_bootstrap_replications = 50)
#' @export
kullback_leibler_test <- function(log_likelihood_differences, block_length,
                                  num_bootstrap_replications, alpha = 0.05) {
  P          <- nrow(log_likelihood_differences)
  K          <- ncol(log_likelihood_differences)
  Mean_L_k   <- colMeans(log_likelihood_differences, na.rm = TRUE)
  V_hat_full <- estimate_long_run_covariance(log_likelihood_differences, block_length)
  V_k        <- diag(V_hat_full)
  V_k[V_k <= 1e-10] <- 1e-10
  std_dev_k  <- sqrt(V_k)
  
  T_k_unscaled     <- Mean_L_k / std_dev_k
  T_KLIC_statistic <- max(T_k_unscaled, na.rm = TRUE)
  
  bootstrap_stats_all_k <- matrix(NA, nrow = num_bootstrap_replications, ncol = K)
  for (b in 1:num_bootstrap_replications) {
    boot_sample       <- mbb_resample_data(log_likelihood_differences, block_length)
    boot_Mean_L_k     <- colMeans(boot_sample, na.rm = TRUE)
    boot_T_k_unscaled <- boot_Mean_L_k / std_dev_k
    boot_T_k_centered <- boot_T_k_unscaled - T_k_unscaled
    bootstrap_stats_all_k[b, ] <- boot_T_k_centered
  }
  
  T_C_boot_max       <- apply(bootstrap_stats_all_k, 1, max, na.rm = TRUE)
  p_value_consistent <- mean(T_C_boot_max > T_KLIC_statistic, na.rm = TRUE)
  
  res <- list(
    statistic   = c("T-KLIC" = T_KLIC_statistic),
    p.value     = p_value_consistent,
    method      = "Kullback-Leibler Information Criterion (KLIC) Test",
    data.name   = deparse(substitute(log_likelihood_differences)),
    null.value  = c("max studentized mean NLS differential" = 0),
    alternative = "at least one forecast has lower Kullback-Leibler distance from the true density than the benchmark"
  )
  class(res) <- "htest"
  return(res)
}

#' @title ZP Quantile Loss Reality Check Test
#'
#' @description Implements the studentized Reality Check test for comparing predictive
#' densities based on the ZP quantile loss function of Corradi & Swanson (2006). The test
#' evaluates whether any competing forecast more accurately predicts the probability of the
#' outcome falling below a specified threshold than the benchmark forecast.
#' \itemize{
#'   \item \strong{H0:} \eqn{\max_k E[\mu^2_1(u) - \mu^2_k(u)] \leq 0} -- no competing
#'     forecast has a lower expected squared probability forecast error at threshold \eqn{u}
#'     than the benchmark (Corradi & Swanson, 2006, eq. 7).
#'   \item \strong{H1:} At least one competing forecast has lower expected ZP-loss than the
#'     benchmark.
#' }
#' @param zp_loss_differences A \code{\link[base]{numeric}} matrix (\code{P x K}) of
#'   ZP-loss differences (benchmark ZP-loss minus forecast ZP-loss), where a positive
#'   entry means the competing forecast outperforms the benchmark forecast at that period.
#' @param block_length \code{\link[base]{integer}}. The block length for MBB and HAC
#'   estimation. A commonly used rule of thumb is
#'   \code{block_length} \eqn{\approx T^{1/3}}, where \eqn{T} is the number of
#'   observations (Politis & Romano, 1994). For \code{P = 165}, this gives
#'   approximately 5--6.
#' @param num_bootstrap_replications \code{\link[base]{integer}} number of MBB
#'   bootstrap replications. Default \code{999}; see Davidson & MacKinnon (2000).
#' @param alpha \code{\link[base]{numeric}}. The significance level (default \code{0.05}).
#'
#' @details
#' The ZP loss for forecast \eqn{k} at period \eqn{t} is
#' \deqn{ZP_{t,k} = \left(\mathbf{1}(y_t \leq \tau) - F_k(\tau \mid \hat{y}_{t,k},
#' \hat{\sigma}_{t,k})\right)^2}
#' where \eqn{\tau} is a threshold, \eqn{F_k(\cdot)} is the forecast's predictive CDF at
#' period \eqn{t}, and \eqn{\mathbf{1}(y_t \leq \tau)} is the indicator for a tail event.
#' Interpretively, the threshold \eqn{\tau} defines a tail event of interest (e.g., the
#' 5th percentile of realizations). The ZP loss penalises the squared difference between
#' the predicted probability of this event and whether it actually occurred. A forecast with
#' \strong{lower} ZP loss more accurately calibrates the left-tail probability. The
#' threshold is typically set to a quantile of the realized series (e.g.,
#' \code{quantile(realized, 0.05)}); a lower threshold focuses the test more sharply on
#' extreme left-tail events.
#'
#' The benchmark is treated as a degenerate (point-mass) predictive distribution
#' with \eqn{\hat{\sigma} = 10^{-6}}, which is a conservative choice ensuring the
#' benchmark's ZP loss approximates the Brier score for the tail indicator.
#'
#' Two p-values are returned: \code{p_consistent} and \code{p_conservative}, analogous
#' to those in \code{\link{superior_predictive_ability_test}}.
#'
#' @return An object of class \code{"htest"}. Lower p-values indicate that
#'   at least one competing forecast is significantly better calibrated in the
#'   left-tail than the benchmark forecast.
#'   Also contains \code{p_consistent} and \code{p_conservative}.
#'   Failure to reject H0 suggests no forecast provides significantly better density
#'   fit than the benchmark forecast.
#' @references
#' Corradi, V., & Swanson, N. R. (2006). Predictive density and conditional confidence
#' interval accuracy tests. \emph{Journal of Econometrics}, 135(1--2), 187--228.
#' \doi{10.1016/j.jeconom.2005.07.026}
#'
#' Corradi, V., & Swanson, N. R. (2011). The White Reality Check and some of its recent
#' extensions. In \emph{Festschrift in honor of Halbert L. White}.
#'
#' Davidson, R., & MacKinnon, J. G. (2000). Bootstrap tests: How many bootstraps?
#' \emph{Econometric Reviews}, 19(1), 55--68.
#' \doi{10.1080/07474930008800459}
#' 
#' Hansen, P. R. (2005). A Test for Superior Predictive Ability.
#' \emph{Journal of Business & Economic Statistics}, 23(4), 365--380.
#' \doi{10.1198/073500105000000063}
#'
#' Politis, D. N., & Romano, J. P. (1994). The stationary bootstrap.
#' \emph{Journal of the American Statistical Association}, 89(428), 1303--1313.
#' \doi{10.1080/01621459.1994.10476870}
#'
#' @examples
#' data(metals)
#' # metals: 165 x 15; columns 1-14 are competing forecasts, column 15 is the benchmark
#' benchmark_col      <- 15
#' P                  <- nrow(metals)
#' K_total            <- ncol(metals)
#' K                  <- K_total - 1 # 14 competing forecasts
#' realized           <- metals[, benchmark_col]
#' forecast_variance  <- estimate_forecast_variance(metals, realized = realized,
#'                         benchmark_col = benchmark_col, window_size = 20)
#' comp_cols          <- setdiff(seq_len(K_total), benchmark_col)
#' forecast_sd_models <- sqrt(forecast_variance[, comp_cols])
#' threshold_val      <- quantile(metals[, K_total], 0.05)
#' zp_loss  <- compute_zp(metals, forecast_sd_models,
#'                        threshold = threshold_val, benchmark_col = K_total)
#' zp_diff  <- zp_loss[, K_total] - zp_loss[, comp_cols]
#' reality_check_zp_test(zp_diff, block_length = 5, num_bootstrap_replications = 50)
#' @export
reality_check_zp_test <- function(zp_loss_differences, block_length,
                                  num_bootstrap_replications, alpha = 0.05) {
  P          <- nrow(zp_loss_differences)
  K          <- ncol(zp_loss_differences)
  SP_k       <- colMeans(zp_loss_differences, na.rm = TRUE)
  V_hat_full <- estimate_long_run_covariance(zp_loss_differences, block_length)
  V_k        <- diag(V_hat_full)
  V_k[V_k <= 1e-10] <- 1e-10
  std_dev_k  <- sqrt(V_k)
  
  T_k_unscaled    <- SP_k / std_dev_k
  T_ZP_Consistent <- max(T_k_unscaled, na.rm = TRUE)
  
  bootstrap_consistent_all_k   <- matrix(NA, nrow = num_bootstrap_replications, ncol = K)
  bootstrap_conservative_all_k <- matrix(NA, nrow = num_bootstrap_replications, ncol = K)
  
  for (b in 1:num_bootstrap_replications) {
    boot_sample       <- mbb_resample_data(zp_loss_differences, block_length)
    boot_SP_k         <- colMeans(boot_sample, na.rm = TRUE)
    boot_T_k          <- boot_SP_k / std_dev_k
    bootstrap_consistent_all_k[b, ]   <- boot_T_k - T_k_unscaled
    bootstrap_conservative_all_k[b, ] <- boot_T_k
  }
  
  T_C_boot_max       <- apply(bootstrap_consistent_all_k,   1, max, na.rm = TRUE)
  T_R_boot_max       <- apply(bootstrap_conservative_all_k, 1, max, na.rm = TRUE)
  p_value_consistent   <- mean(T_C_boot_max > T_ZP_Consistent, na.rm = TRUE)
  p_value_conservative <- mean(T_R_boot_max > T_ZP_Consistent, na.rm = TRUE)
  
  res <- list(
    statistic            = c("T-ZP" = T_ZP_Consistent),
    p.value              = p_value_consistent,
    method               = "ZP Quantile Loss Test (Corradi & Swanson, 2006)",
    data.name            = deparse(substitute(zp_loss_differences)),
    null.value           = c("max studentized mean ZP-loss differential" = 0),
    alternative          = "at least one forecast better calibrates the left-tail probability than the benchmark",
    p_consistent         = p_value_consistent,
    p_conservative       = p_value_conservative,
    reject_consistent_zp = p_value_consistent <= alpha
  )
  class(res) <- "htest"
  return(res)
}
