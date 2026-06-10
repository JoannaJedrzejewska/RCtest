#' @importFrom stats aggregate
#' @importFrom rlang .data
#' @importFrom utils globalVariables
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom ggrepel geom_text_repel
NULL

utils::globalVariables(c("Significant", "value", "Value", "Model",
                         "Time", "ERC_Weight", "Risk_Contribution"))

# @noRd
get_top_bottom_models <- function(data, sort_var, n_top = 2, n_bottom = 2) {
  top_models    <- head(data[order(-data[[sort_var]]), ], n_top)
  bottom_models <- head(data[order( data[[sort_var]]), ], n_bottom)
  unique(rbind(top_models, bottom_models))
}

#' @title Plot Performance Metrics Comparison
#'
#' @description Generates a grid of scatter plots for Root Mean Squared Error
#' (RMSE), Normalized Root Mean Squared Error (N-RMSE), Mean Absolute Error
#' (MAE), and Mean Absolute Scaled Error (MASE) plotted against the Equally
#' Weighted Risk Contribution (ERC) Weight of each competing forecast. These
#' metrics are used for visualization only and are computed directly from
#' forecast errors; they are distinct from the loss functions used in the
#' WRC/SPA/CPA hypothesis tests (see \code{\link{run_comprehensive_erc_analysis}}).
#'
#' Each panel shows:
#' \itemize{
#'   \item \strong{Points}: one per model. The radius (size) of each circle is
#'     proportional to the model's \emph{Risk Contribution}, defined as
#'     \eqn{\text{RMSE}_k \times \text{ERC Weight}_k}. Larger circles indicate
#'     forecasts that contribute more to the portfolio-level risk.
#'   \item \strong{Dashed line}: an OLS regression line of the error metric
#'     (x-axis) on the ERC Weight (y-axis), showing whether higher-weighted
#'     forecasts tend to have systematically lower or higher error.
#' }
#'
#' \strong{Note on MASE scaling:} The MASE denominator used here is
#' \code{mean(abs(diff(benchmark)))}, where \code{benchmark} is the
#' benchmark forecast column extracted from \code{forecast_matrix}.
#' This differs from the scaling used in \code{\link{run_comprehensive_erc_analysis}},
#' where the denominator is \code{mean(abs(diff(realizations_raw)))} computed
#' from the actual realised values. The two denominators coincide when the
#' benchmark is a random-walk or historical-average forecast whose one-step-ahead
#' forecasts equal the lagged realised value, but will diverge otherwise.
#' The \code{plot_performance_metrics} function does not accept a separate
#' realisation vector, so the benchmark forecast series is used as a
#' proxy for the naive forecast scale. MASE values produced by this function
#' are therefore intended for \emph{visual comparison across forecasts only}
#' and should not be directly compared to MASE values from the hypothesis
#' tests in \code{\link{run_comprehensive_erc_analysis}}.
#' 
#' @section Equally Weighted Risk Contribution (ERC) Weight:
#' ERC weights are portfolio weights assigned so that every competing forecast
#' contributes an equal share to the total portfolio risk (measured here by
#' forecast error dispersion). Formally, weights \eqn{w_k} are chosen so that
#' \eqn{w_k \cdot \sigma_k = c} for all \eqn{k},
#' where \eqn{\sigma_k} is a measure of forecast \eqn{k}'s risk and
#' \eqn{c = \frac{1}{K}\sum_{k=1}^{K} w_k \sigma_k} is the common
#' per-forecast risk budget determined endogenously by the equal-contribution
#' constraint (Maillard et al., 2010).
#' When \code{weights} are supplied by the user, they are treated as
#' pre-computed ERC weights and normalised to sum to one. When
#' \code{weights = NULL}, equal weights \eqn{1/K} are used as a baseline.
#'
#' @param forecast_matrix Matrix or data frame of dimension \code{P x K_total},
#'   where \code{P} is the number of forecast periods, columns \code{1} to
#'   \code{K_total - 1} are competing model forecasts, and the last column
#'   (or \code{benchmark_col}) is the benchmark.
#' @param weights Optional numeric vector of length \code{K_total - 1} giving
#'   ERC weights for each competing forecast. If \code{NULL} (default), equal
#'   weights \eqn{1/K} are used.
#' @param benchmark_col Index or name of the benchmark column. Defaults to the
#'   last column.
#'
#' @return A \code{\link[gtable]{gtable}} object produced by
#'   \code{\link[gridExtra]{grid.arrange}} containing four panels arranged 
#'   in a 2x2 grid: RMSE (top-left), N-RMSE (top-right), MAE (bottom-left),
#'   MASE (bottom-right). Each panel is a
#'   \code{ggplot} object and can be extracted individually if needed.
#'
#' @seealso \code{\link{run_comprehensive_erc_analysis}},
#'   \code{\link{plot_cumulative_loss}}
#'
#' @references
#' Maillard, S., Roncalli, T., & Teïletche, J. (2010). The Properties of
#' Equally Weighted Risk Contributions Portfolios. \emph{The Journal of
#' Portfolio Management}, 36(4), 60--70. \doi{10.3905/jpm.2010.36.4.060}
#'
#' @examples
#' data(metals)
#' K_models       <- ncol(metals) - 1
#' custom_weights <- (K_models:1) / sum(K_models:1)
#' plot_performance_metrics(metals, weights = custom_weights, benchmark_col = 15)
#' @export
plot_performance_metrics <- function(forecast_matrix, weights = NULL,
                                     benchmark_col = ncol(forecast_matrix)) {
  if (is.matrix(forecast_matrix)) forecast_matrix <- as.data.frame(forecast_matrix)
  
  if (is.character(benchmark_col))
    benchmark_col <- which(colnames(forecast_matrix) == benchmark_col)
  
  if (length(benchmark_col) == 0) stop("Benchmark column not found in forecast_matrix.")
  
  K_total       <- ncol(forecast_matrix)
  benchmark     <- as.numeric(forecast_matrix[, benchmark_col])
  model_indices <- setdiff(seq_len(K_total), benchmark_col)
  models_only   <- forecast_matrix[, model_indices, drop = FALSE]
  K_models      <- length(model_indices)
  
  errors     <- sweep(as.matrix(models_only), 1, benchmark, "-")
  rmse_vals  <- sqrt(colMeans(errors^2, na.rm = TRUE))
  mae_vals   <- colMeans(abs(errors), na.rm = TRUE)
  nrmse_vals <- rmse_vals / mean(abs(benchmark), na.rm = TRUE)
  naive_mae  <- mean(abs(diff(benchmark)), na.rm = TRUE)
  mase_vals  <- mae_vals / naive_mae
  
  if (is.null(weights)) {
    weights <- rep(1 / K_models, K_models)
  } else {
    if (length(weights) != K_models)
      stop(paste("'weights' must have length equal to the number of competing forecasts:",
                 K_models))
    weights <- weights / sum(weights)
  }
  
  df <- data.frame(
    Model             = colnames(models_only),
    RMSE              = rmse_vals,
    NRMSE             = nrmse_vals,
    MAE               = mae_vals,
    MASE              = mase_vals,
    ERC_Weight        = weights,
    Risk_Contribution = rmse_vals * weights,
    stringsAsFactors  = FALSE
  )
  
  pad   <- 0.05
  y_min <- floor(  min(df$ERC_Weight) / pad) * pad - pad
  y_max <- ceiling(max(df$ERC_Weight) / pad) * pad + pad
  
  rmse_labeled  <- get_top_bottom_models(df, "RMSE")
  nrmse_labeled <- get_top_bottom_models(df, "NRMSE")
  mae_labeled   <- get_top_bottom_models(df, "MAE")
  mase_labeled  <- get_top_bottom_models(df, "MASE")
  
  create_p <- function(data, x_var, title, x_lab, labeled_data, color) {
    ggplot(data, aes(x = .data[[x_var]], y = ERC_Weight)) +
      geom_point(aes(size = Risk_Contribution), color = color, alpha = 0.7) +
      geom_smooth(method = "lm", color = "#566573", linetype = "dashed",
                  se = FALSE) +
      geom_text_repel(data          = labeled_data,
                      aes(label     = Model),
                      size          = 2.8, color = "black",
                      box.padding   = 0.4, point.padding = 0.2,
                      segment.color = "grey50", segment.size = 0.4,
                      max.overlaps  = Inf, force = 2, direction = "both") +
      labs(title = title, x = x_lab, y = "ERC Weight") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title      = element_text(size = 11, hjust = 0.5),
            axis.title      = element_text(size = 10)) +
      scale_size_continuous(range = c(2, 8)) +
      coord_cartesian(ylim = c(y_min, y_max))
  }
  
  colors <- c("#2874A6", "#8E44AD", "#148F77", "#D35400")
  p1 <- create_p(df, "RMSE",  "Root Mean Squared Error (RMSE)",    "RMSE",
                 rmse_labeled,  colors[1])
  p2 <- create_p(df, "NRMSE", "Normalized RMSE (N-RMSE)",          "N-RMSE",
                 nrmse_labeled, colors[2])
  p3 <- create_p(df, "MAE",   "Mean Absolute Error (MAE)",         "MAE",
                 mae_labeled,   colors[3])
  p4 <- create_p(df, "MASE",  "Mean Absolute Scaled Error (MASE)", "MASE",
                 mase_labeled,  colors[4])
  
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
}

#' @title Plot Density Forecast
#'
#' @description Generates a kernel density plot of a predictive distribution,
#' highlighting the point forecast and a symmetric credible interval at a
#' specified level. The predictive distribution can be constructed by treating
#' the cross-sectional spread of model forecasts at a given period as a proxy
#' for forecast uncertainty (see \code{\link{compute_crps}} for the centring
#' convention).
#'
#' @param full_distribution Numeric vector of forecast samples drawn from the
#'   predictive distribution (e.g., forecasts from all models at one time
#'   period, optionally centred and shifted as in \code{\link{compute_crps}}).
#' @param point_forecast Single numeric value: the point forecast to highlight
#'   (e.g., the mean or median of \code{full_distribution}, or one model's
#'   forecast).
#' @param title Character string for the plot title. Default is
#'   \code{"Density Forecast"}.
#' @param ci_level Numeric confidence interval level strictly between 0 and 1.
#'   Default is 0.90, producing lower and upper quantiles at
#'   \code{(1 - ci_level) / 2} and \code{1 - (1 - ci_level) / 2}.
#'
#' @details
#' To build a pseudo-predictive distribution from the \code{metals} dataset,
#' centre the cross-sectional model forecasts at period \code{t} around their
#' mean, following the same convention used in \code{\link{compute_crps}}:
#' \code{dist_t <- forecasts[t, ] - forecasts[t, k] + mean(forecasts[t, ])}.
#' This preserves cross-sectional spread while recentring on the cross-sectional
#' mean rather than on model \code{k}'s own point forecast.
#'
#' @return A \code{ggplot} object. The plot shows a kernel density curve
#'   (blue fill), a red dashed vertical line at \code{point_forecast}, and
#'   orange dotted vertical lines at the lower and upper quantiles of the
#'   credible interval. The subtitle reports the numeric bounds of the interval.
#'
#' @seealso \code{\link{compute_crps}}, \code{\link{run_comprehensive_erc_analysis}}'
#' @examples
#' data(metals)
#' t       <- 100
#' K       <- ncol(metals) - 1
#' dist_t  <- as.numeric(metals[t, 1:K]) - metals[t, 7] +
#'             mean(as.numeric(metals[t, 1:K]))
#' pt_fcst <- mean(as.numeric(metals[t, 1:K]))
#' plot_density_forecast(dist_t, pt_fcst,
#'                       title    = "Predictive Density at t = 100",
#'                       ci_level = 0.90)
#' @export
plot_density_forecast <- function(full_distribution, point_forecast,
                                  title    = "Density Forecast",
                                  ci_level = 0.90) {
  df     <- data.frame(value = full_distribution)
  lower  <- (1 - ci_level) / 2
  upper  <- 1 - lower
  q_low  <- quantile(full_distribution, lower, na.rm = TRUE)
  q_high <- quantile(full_distribution, upper, na.rm = TRUE)
  
  ci_label <- paste0(round(ci_level * 100), "% Confidence Interval (CI):")
  
  ggplot(df, aes(x = value)) +
    geom_density(fill = "lightblue", alpha = 0.7, color = "darkblue") +
    geom_vline(xintercept = point_forecast,   color = "red",    linetype = "dashed") +
    geom_vline(xintercept = c(q_low, q_high), color = "orange", linetype = "dotted") +
    labs(title    = title,
         subtitle = paste(ci_label, round(q_low, 2), "-", round(q_high, 2)),
         x        = "Forecast Value",
         y        = "Density") +
    theme_minimal()
}

#' @title Plot Cumulative Loss Differences
#'
#' @description Generates a time-series plot of cumulative squared error (MSE)
#' loss differences between each competing forecast and the benchmark. A positive
#' value at time \eqn{t} means the competing forecast has accumulated lower squared
#' errors than the benchmark up to that point (i.e., the forecast is outperforming
#' the benchmark cumulatively).
#'
#' @details
#' The cumulative loss difference for forecast \eqn{k} at time \eqn{t} is:
#' \deqn{\text{CLD}_{k,t} =
#'   \sum_{s=1}^{t} \left( e_{\text{bench},s}^2 - e_{k,s}^2 \right)}
#' where \eqn{e_{k,s} = y_s - \hat{y}_{k,s}} is the forecast error of forecast
#' \eqn{k} at period \eqn{s}. The function squares the values in
#' \code{data_matrix} directly, so \strong{pre-computed forecast errors}
#' (not raw forecasts) must be passed for the result to represent cumulative
#' MSE differences.
#'
#' A positive \eqn{\text{CLD}_{k,t}} means forecast \eqn{k} has accumulated lower
#' squared errors than the benchmark up to period \eqn{t}. A negative value
#' means the benchmark has been more accurate up to that point.
#'
#' @param data_matrix Matrix or data frame of dimension \code{P x K_total},
#'   where rows are time periods and columns are \strong{pre-computed forecast
#'   errors} (\eqn{y_t - \hat{y}_{k,t}}) for each forecast including the
#'   benchmark. \strong{Do not pass raw forecasts}: the function squares column
#'   values directly, so passing raw forecasts produces cumulative sums of
#'   squared forecast levels, not cumulative MSE differences.
#' @param benchmark_col Index or name of the benchmark column. Defaults to the
#'   last column.
#'
#' @return A \code{ggplot} object. The y-axis shows the cumulative MSE
#'   difference; forecasts above zero at the right edge have outperformed the
#'   benchmark over the full evaluation window. The top 2 and bottom 2 forecasts
#'   (by final cumulative loss difference) are labelled directly on the plot.
#'
#' @seealso \code{\link{white_reality_check}},
#'   \code{\link{plot_performance_metrics}}
#'
#' @examples
#' data(metals)
#' P        <- nrow(metals)
#' K_total  <- ncol(metals)
#' realized <- c(metals[-1, K_total], metals[P, K_total])
#' errors   <- sweep(as.matrix(metals), 1, realized, "-")
#' p <- plot_cumulative_loss(errors, benchmark_col = K_total)
#' print(p)
#' @export
plot_cumulative_loss <- function(data_matrix,
                                 benchmark_col = ncol(data_matrix)) {
  if (is.matrix(data_matrix)) data_matrix <- as.data.frame(data_matrix)
  
  if (is.character(benchmark_col)) {
    bench_idx <- which(colnames(data_matrix) == benchmark_col)
  } else {
    bench_idx <- benchmark_col
  }
  
  if (length(bench_idx) == 0 || bench_idx > ncol(data_matrix))
    stop("Benchmark column not found.")
  
  bench_name <- colnames(data_matrix)[bench_idx]
  P          <- nrow(data_matrix)
  n_cols     <- ncol(data_matrix)
  model_cols <- setdiff(seq_len(n_cols), bench_idx)
  K_models   <- length(model_cols)
  
  bench_sq  <- data_matrix[, bench_idx]^2
  model_mat <- as.matrix(data_matrix[, model_cols, drop = FALSE])
  cum_diff  <- apply(model_mat, 2, function(col) cumsum(bench_sq - col^2))
  
  df_long <- data.frame(
    Value = as.vector(cum_diff),
    Model = rep(colnames(data_matrix)[model_cols], each = P),
    Time  = rep(seq_len(P), length(model_cols))
  )
  
  df_long <- df_long[df_long$Model != bench_name, ]
  
  df_final     <- df_long[df_long$Time == max(df_long$Time), ]
  top_2        <- head(df_final[order(-df_final$Value), ], 2)
  bot_2        <- head(df_final[order( df_final$Value), ], 2)
  labeled_data <- unique(rbind(top_2, bot_2))
  
  colors <- c(
    "#2874A6", "#8E44AD", "#148F77", "#D35400", "#566573", "#95A5A6",
    "#C0392B", "#27AE60", "#F39C12", "#2980B9", "#884EA0", "#16A085",
    "#E67E22", "#7F8C8D", "#34495E", "#1ABC9C", "#9B59B6", "#F1C40F",
    "#E74C3C"
  )
  
  unique_models      <- unique(df_long$Model)
  colors_used        <- rep(colors, length.out = length(unique_models))
  names(colors_used) <- unique_models
  
  p <- ggplot(df_long, aes(x = Time, y = Value, color = Model, group = Model)) +
    geom_line(alpha = 0.9) +
    geom_text_repel(data          = labeled_data,
                    aes(label     = Model),
                    size          = 3, color = "black",
                    box.padding   = 0.5, point.padding = 0.3,
                    nudge_x       = 5,  direction = "y",
                    segment.color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = colors_used) +
    labs(title = paste("Cumulative Loss vs", bench_name),
         x     = "Periods",
         y     = "Cumulative MSE Difference",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title    = element_text(size = 9),
          legend.text     = element_text(size = 8))
  
  return(p)
}