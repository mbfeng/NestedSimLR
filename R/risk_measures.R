################################################################################
# risk_measures.R
#
# The nested statistics zeta(mu(theta)) of Section 7, computed from a vector
# of (estimated or true) conditional means, and the error summaries (bias,
# variance, MSE) reported in Tables 1 and 3.
#
# Paper notation -> code:
#   mu(theta_1), ..., mu(theta_M)      -> vector x
#   VaR level alpha                    -> alpha
#   threshold u in I(mu > u), (mu - u)^+ and ((mu - u)^+)^2 -> threshold
################################################################################

# The four risk measures of Section 7.1 for one vector x of conditional means:
#   VaR          the alpha-quantile of x (R's default "type 7" quantile,
#                as in the original code)
#   indicator    mean of I(x > threshold)
#   hockeystick  mean of (x - threshold)^+
#   square       mean of ((x - threshold)^+)^2
risk_measures <- function(x, alpha, threshold) {
  excess <- pmax(x - threshold, 0)
  c(
    VaR = unname(quantile(x, alpha)),
    indicator = mean(x > threshold),
    hockeystick = mean(excess),
    square = mean(excess^2)
  )
}

# Risk measures for many macro runs at once.
# estimates: matrix with one column per macro run (rows = outer scenarios).
# Returns a data frame with one row per macro run.
risk_measures_by_macro <- function(estimates, alpha, threshold) {
  n_macro <- ncol(estimates)
  result <- matrix(NA, n_macro, 4,
                   dimnames = list(NULL, c("VaR", "indicator", "hockeystick",
                                           "square")))
  for (k in seq_len(n_macro)) {
    result[k, ] <- risk_measures(estimates[, k], alpha, threshold)
  }
  as.data.frame(result)
}

# Bias, variance and MSE of risk-measure estimates over K macro runs, with
# the Monte Carlo standard error of the MSE.
#
# estimates: data frame or matrix, one row per macro run, one column per
#   risk measure; truth: named vector of the true values.
# Returns a data frame with one row per risk measure.
#
# The variance uses 1/K (not 1/(K - 1)), as in the original code, so that
# MSE = bias^2 + variance exactly.
#
# top1_share is the fraction of the MSE that comes from the worst 1% of the
# macro runs. Above 0.5 the squared errors are so heavy-tailed that the MSE
# and its standard error are dominated by a few runs and are unreliable
# (see README, Verification).
error_summary <- function(estimates, truth) {
  estimates <- as.matrix(estimates)
  measures <- colnames(estimates)
  n_macro <- nrow(estimates)
  n_worst <- ceiling(n_macro / 100)
  result <- data.frame(measure = measures, bias = NA, var = NA, mse = NA,
                       mse_se = NA, top1_share = NA)
  for (i in seq_along(measures)) {
    x <- estimates[, measures[i]]
    error <- x - truth[[measures[i]]]
    squared_error <- error^2
    worst <- sort(squared_error, decreasing = TRUE)[seq_len(n_worst)]
    result$bias[i] <- mean(error)
    result$var[i] <- mean((x - mean(x))^2)
    result$mse[i] <- mean(squared_error)
    result$mse_se[i] <- standard_error(squared_error)
    result$top1_share[i] <- sum(worst) / sum(squared_error)
  }
  result
}
