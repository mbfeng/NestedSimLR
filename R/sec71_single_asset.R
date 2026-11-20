################################################################################
# sec71_single_asset.R
#
# Section 7.1: single-asset enterprise risk management (ERM) example.
#
# The portfolio is a straddle (a call plus a put, both with strike K = 110)
# on one stock in the Black-Scholes model.
#   outer scenario theta  = stock price S_tau at the risk horizon tau
#   inner simulation      = terminal price S_T given S_tau = theta (risk
#                           neutral); output g = discounted straddle payoff
#   conditional mean      = mu(theta) = straddle price at tau (closed form)
#   LR input              = S_T (only its distribution depends on theta)
#
# Designs compared (names as in the paper):
#   optimal     the optimal LR design (LP (11) + pooled LR estimator)
#   SNS         standard nested simulation with the optimal design's budget
#               Gamma: ceiling(Gamma^(2/3)) outer x ceiling(Gamma^(1/3)) inner
#   SNS+        standard nested simulation with N inner replications at all
#               M scenarios
#   regression  one inner replication at each of ceiling(Gamma) scenarios,
#               least squares on weighted Laguerre polynomials of degree 0-2
################################################################################

# Model parameters of Section 7.1.
single_asset_params <- function() {
  params <- list(
    s0 = 100,      # initial stock price
    mu = 0.05,     # real-world drift (outer simulation)
    rf = 0.02,     # risk-free rate (inner simulation)
    div = 0,       # dividend yield
    vol = 0.3,     # volatility
    k = 110,       # strike of the call and the put
    big_t = 2,     # maturity T
    tau = 1 / 4    # risk horizon
  )
  params$t_in <- params$big_t - params$tau # time from tau to maturity
  # Outer simulation (real world): log(S_tau / s0) ~ Normal(meanlog_out, sdlog_out^2)
  params$meanlog_out <- (params$mu - params$div - params$vol^2 / 2) * params$tau
  params$sdlog_out <- params$vol * sqrt(params$tau)
  # Inner simulation (risk neutral): log(S_T / theta) ~ Normal(meanlog_in, sdlog_in^2)
  params$meanlog_in <- (params$rf - params$div - params$vol^2 / 2) * params$t_in
  params$sdlog_in <- params$vol * sqrt(params$t_in)
  params
}

# The M outer scenarios: the 1/(M+1), 2/(M+1), ..., M/(M+1) quantiles of S_tau
# (fixed, not random).
single_asset_scenarios <- function(n_outer, params) {
  probabilities <- seq_len(n_outer) / (n_outer + 1)
  qlnorm(probabilities, log(params$s0) + params$meanlog_out, params$sdlog_out)
}

# True conditional means mu(theta): the straddle price at time tau.
single_asset_true_value <- function(theta, params) {
  bs_price_straddle(theta, params$k, params$rf, params$vol, params$t_in,
                    params$div)
}

# n inner replications at ONE scenario theta_j. Returns the outputs y
# (discounted straddle payoffs) and the LR input log(S_T).
single_asset_inner <- function(theta_j, n, params) {
  log_terminal <- log(theta_j) + rnorm(n, params$meanlog_in, params$sdlog_in)
  payoff <- abs(exp(log_terminal) - params$k) # straddle: |S_T - K|
  list(y = exp(-params$rf * params$t_in) * payoff, log_x = log_terminal)
}

# The optimal design for the scenarios theta (deterministic, no simulation).
# The LR input log(S_T) has log-mean log(theta) + meanlog_in under scenario
# theta and log-sd sdlog_in under every scenario.
single_asset_design <- function(theta, target_ess, params) {
  log_means <- log(theta) + params$meanlog_in
  psi <- psi_lognormal(log_means, params$sdlog_in)
  optimal_design(psi, target_ess, section = "7.1/7.3")
}

# One macro run of the optimal design: estimates of mu(theta_i) for all M
# scenarios, from inner replications at the sampled scenarios only.
single_asset_optimal <- function(theta, design, params) {
  log_means <- log(theta) + params$meanlog_in
  inner_samples <- list()
  for (b in seq_along(design$sampled)) {
    j <- design$sampled[b]
    sim <- single_asset_inner(theta[j], design$n_inner[b], params)
    log_lr <- log_lr_lognormal(sim$log_x, log_means, log_means[j],
                               params$sdlog_in)
    inner_samples[[b]] <- list(y = sim$y, log_lr = log_lr)
  }
  pooled_estimate(inner_samples, design$weights)
}

# Standard nested simulation: the average of n_inner replications at every
# scenario. All replications are drawn at once; column i = scenario i.
single_asset_standard <- function(theta, n_inner, params) {
  n_outer <- length(theta)
  log_returns <- matrix(rnorm(n_inner * n_outer, params$meanlog_in,
                              params$sdlog_in), nrow = n_inner)
  terminal <- exp(log_returns) * rep(theta, each = n_inner) # S_T = theta e^Z
  payoff <- abs(terminal - params$k)
  exp(-params$rf * params$t_in) * colMeans(payoff)
}

# Weighted Laguerre polynomials of degree 0, ..., degree (at most 3),
# evaluated at x: L_n(x) = exp(-x/2) * (Laguerre polynomial of degree n).
# Returns a data frame with columns L0, L1, ...
laguerre_basis <- function(x, degree = 2) {
  weight <- exp(-x / 2)
  basis <- data.frame(L0 = weight)
  if (degree >= 1) basis$L1 <- weight * (1 - x)
  if (degree >= 2) basis$L2 <- weight * (1 - 2 * x + x^2 / 2)
  if (degree >= 3) basis$L3 <- weight * (1 - 3 * x + 3 * x^2 / 2 - x^3 / 6)
  basis
}

# Regression design (Longstaff-Schwartz type): one inner replication at each
# design point, least-squares fit of the outputs on the Laguerre basis of
# theta / s0, then prediction at the scenarios theta.
single_asset_regression <- function(theta, design_points, params) {
  y <- single_asset_standard(design_points, 1, params)
  training <- laguerre_basis(design_points / params$s0)
  training$y <- y
  fit <- lm(y ~ L0 + L1 + L2, data = training)
  prediction_basis <- laguerre_basis(theta / params$s0)
  unname(predict(fit, newdata = prediction_basis))
}

# Sizes of the comparison designs that use the optimal design's budget.
single_asset_matched_sizes <- function(budget) {
  list(sns_outer = ceiling(budget^(2 / 3)),
       sns_inner = ceiling(budget^(1 / 3)),
       regression_points = ceiling(budget))
}

# True risk measures, from n_true equally spaced quantiles of S_tau (10^8 in
# the paper). The conditional means are computed in chunks to limit memory.
single_asset_true_risks <- function(alpha, threshold, params, n_true = 1e8,
                                    chunk = 1e7) {
  mu_true <- numeric(n_true)
  for (start in seq(1, n_true, by = chunk)) {
    i <- start:min(start + chunk - 1, n_true)
    theta <- qlnorm(i / (n_true + 1), log(params$s0) + params$meanlog_out,
                    params$sdlog_out)
    mu_true[i] <- single_asset_true_value(theta, params)
  }
  risk_measures(mu_true, alpha, threshold)
}
