################################################################################
# risk_measures.R                                                              #
# The nested statistics zeta(mu(theta)) of Section 7, computed from a vector   #
# of (estimated or true) conditional means, and the error summaries (bias,     #
# variance, MSE) reported in Tables 1 and 3                                    #
# Paper notation -> code:                                                      #
# -- mu(theta_1), ..., mu(theta_M) -> vector x                                 #
# -- VaR level alpha -> alpha                                                  #
# -- threshold u in I(mu > u), (mu - u)^+ and ((mu - u)^+)^2 -> threshold      #
################################################################################

risk_measures <- function(x, alpha, threshold) {
  # The four risk measures of Section 7.1 for one vector of conditional means
  # Inputs:
  # -- x: vector of conditional means
  # -- alpha: VaR level
  # -- threshold: loss threshold u
  # Output: named vector of 4 values
  # -- VaR: the alpha-quantile of x (R's default type-7 quantile, as in the
  #    original code)
  # -- indicator: mean of I(x > u)
  # -- hockeystick: mean of (x - u)^+
  # -- square: mean of ((x - u)^+)^2

  excess <- pmax(x - threshold, 0) # (x - u)^+

  return(c(
    VaR = unname(quantile(x, alpha)),
    indicator = mean(x > threshold),
    hockeystick = mean(excess),
    square = mean(excess^2)
  ))
}

risk_measures_by_macro <- function(estimates, alpha, threshold) {
  # Risk measures for many macro runs at once
  # Inputs:
  # -- estimates: matrix, one column per macro run (rows = outer scenarios)
  # -- alpha, threshold: as in risk_measures()
  # Output:
  # -- data frame, one row per macro run

  n_macro <- ncol(estimates) # number of macro runs
  result <- matrix(NA, n_macro, 4,
                   dimnames = list(NULL, c("VaR", "indicator", "hockeystick",
                                           "square"))) # placeholder
  for (k in 1:n_macro) {
    result[k, ] <- risk_measures(estimates[, k], alpha, threshold)
  }

  return(as.data.frame(result))
}

error_summary <- function(estimates, truth) {
  # Bias, variance and MSE of risk-measure estimates over K macro runs, with
  # the Monte Carlo standard error of the MSE
  # Inputs:
  # -- estimates: data frame or matrix, one row per macro run, one column
  #    per risk measure
  # -- truth: named vector of the true values
  # Output:
  # -- data frame, one row per risk measure: bias, var, mse, mse_se,
  #    top1_share
  # The variance uses 1/K (not 1/(K - 1)), as in the original code, so that
  # MSE = bias^2 + variance exactly.
  # top1_share is the fraction of the MSE that comes from the worst 1% of the
  # macro runs. Above 0.5 the squared errors are so heavy-tailed that the MSE
  # and its standard error are dominated by a few runs and are unreliable
  # (see README, Verification).

  estimates <- as.matrix(estimates)
  measures <- colnames(estimates) # names of the risk measures
  n_macro <- nrow(estimates) # number of macro runs (K)
  n_worst <- ceiling(n_macro / 100) # the worst 1% of the macro runs
  result <- data.frame(measure = measures, bias = NA, var = NA, mse = NA,
                       mse_se = NA, top1_share = NA) # placeholder
  for (i in seq_along(measures)) {
    x <- estimates[, measures[i]]
    error <- x - truth[[measures[i]]]
    squared_error <- error^2
    worst <- sort(squared_error, decreasing = TRUE)[1:n_worst]
    result$bias[i] <- mean(error)
    result$var[i] <- mean((x - mean(x))^2)
    result$mse[i] <- mean(squared_error)
    result$mse_se[i] <- standard_error(squared_error)
    result$top1_share[i] <- sum(worst) / sum(squared_error)
  }

  return(result)
}
