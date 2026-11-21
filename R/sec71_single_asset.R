################################################################################
# sec71_single_asset.R                                                         #
# Section 7.1: single-asset enterprise risk management (ERM) example           #
# A straddle (a call plus a put, both with strike K = 110) on one stock in the #
# Black-Scholes model                                                          #
# -- outer scenario theta = stock price S_tau at the risk horizon tau          #
# -- inner simulation = terminal price S_T given S_tau = theta (risk neutral); #
#    output g = discounted straddle payoff                                     #
# -- conditional mean mu(theta) = straddle price at tau (closed form)          #
# -- LR input = S_T (only its distribution depends on theta)                   #
# Designs compared (names as in the paper):                                    #
# -- optimal: the optimal LR design (LP (11) + pooled LR estimator)            #
# -- SNS: standard nested simulation with the optimal design's budget Gamma,   #
#    ceiling(Gamma^(2/3)) outer x ceiling(Gamma^(1/3)) inner                   #
# -- SNS+: standard nested simulation, N inner replications at all M scenarios #
# -- regression: one inner replication at each of ceiling(Gamma) scenarios,    #
#    least squares on weighted Laguerre polynomials of degree 0-2              #
################################################################################

# ===== Model
single_asset_params <- function() {
  # Model parameters of Section 7.1
  # Output:
  # -- list of parameters, including the lognormal parameters of the outer
  #    and inner simulations

  params <- list(
    s0 = 100, # initial stock price
    mu = 0.05, # real-world drift (outer simulation)
    rf = 0.02, # risk-free rate (inner simulation)
    div = 0, # dividend yield
    vol = 0.3, # volatility
    k = 110, # strike of the call and the put
    big_t = 2, # maturity T
    tau = 1 / 4 # risk horizon
  )
  params$t_in <- params$big_t - params$tau # time from tau to maturity
  #-- outer simulation (real world): log(S_tau/s0) ~ N(meanlog_out, sdlog_out^2)
  params$meanlog_out <- (params$mu - params$div - params$vol^2 / 2) * params$tau
  params$sdlog_out <- params$vol * sqrt(params$tau)
  #-- inner simulation (risk neutral): log(S_T/theta) ~ N(meanlog_in, sdlog_in^2)
  params$meanlog_in <- (params$rf - params$div - params$vol^2 / 2) * params$t_in
  params$sdlog_in <- params$vol * sqrt(params$t_in)

  return(params)
}

single_asset_scenarios <- function(n_outer, params) {
  # The M outer scenarios (fixed, not random): the 1/(M+1), ..., M/(M+1)
  # quantiles of S_tau
  # Inputs:
  # -- n_outer: number of scenarios M
  # -- params: single_asset_params()
  # Output:
  # -- vector of M stock prices

  probabilities <- (1:n_outer) / (n_outer + 1)

  return(qlnorm(probabilities, log(params$s0) + params$meanlog_out,
                params$sdlog_out))
}

single_asset_true_value <- function(theta, params) {
  # True conditional means mu(theta): the straddle price at time tau
  # Inputs:
  # -- theta: vector of stock prices at tau
  # -- params: single_asset_params()
  # Output:
  # -- vector of straddle prices

  return(bs_price_straddle(theta, params$k, params$rf, params$vol,
                           params$t_in, params$div))
}

single_asset_inner <- function(theta_j, n, params) {
  # n inner replications at ONE scenario theta_j
  # Inputs:
  # -- theta_j: stock price at tau
  # -- n: number of replications
  # -- params: single_asset_params()
  # Output: a list of 2 items
  # -- y: discounted straddle payoffs
  # -- log_x: simulated log(S_T), the LR input

  log_terminal <- log(theta_j) + rnorm(n, params$meanlog_in, params$sdlog_in)
  payoff <- abs(exp(log_terminal) - params$k) # straddle: |S_T - K|

  return(list(y = exp(-params$rf * params$t_in) * payoff, log_x = log_terminal))
}

# ===== Designs
single_asset_design <- function(theta, target_ess, params) {
  # The optimal design for the scenarios theta (deterministic)
  # Inputs:
  # -- theta: vector of M stock prices at tau
  # -- target_ess: N
  # -- params: single_asset_params()
  # Output:
  # -- design, see optimal_design()
  # The LR input log(S_T) has log-mean log(theta) + meanlog_in under scenario
  # theta and log-sd sdlog_in under every scenario.

  log_means <- log(theta) + params$meanlog_in
  psi <- psi_lognormal(log_means, params$sdlog_in)

  return(optimal_design(psi, target_ess, section = "7.1/7.3"))
}

single_asset_optimal <- function(theta, design, params) {
  # One macro run of the optimal design
  # Inputs:
  # -- theta: vector of M stock prices at tau
  # -- design: single_asset_design()
  # -- params: single_asset_params()
  # Output:
  # -- vector of M estimates of mu(theta_i), from inner replications at the
  #    sampled scenarios only

  log_means <- log(theta) + params$meanlog_in # log-means of the LR input
  inner_samples <- list() # placeholder, one element per sampled scenario
  for (b in seq_along(design$sampled)) {
    j <- design$sampled[b] # index of the b-th sampled scenario
    sim <- single_asset_inner(theta[j], design$n_inner[b], params)
    log_lr <- log_lr_lognormal(sim$log_x, log_means, log_means[j],
                               params$sdlog_in)
    inner_samples[[b]] <- list(y = sim$y, log_lr = log_lr)
  }

  return(pooled_estimate(inner_samples, design$weights))
}

single_asset_standard <- function(theta, n_inner, params) {
  # Standard nested simulation: the average of n_inner replications at every
  # scenario
  # Inputs:
  # -- theta: vector of stock prices at tau
  # -- n_inner: number of inner replications per scenario
  # -- params: single_asset_params()
  # Output:
  # -- vector of estimates of mu(theta_i)
  # All replications are drawn at once; column i = scenario i.

  n_outer <- length(theta) # number of scenarios
  log_returns <- matrix(rnorm(n_inner * n_outer, params$meanlog_in,
                              params$sdlog_in), nrow = n_inner)
  terminal <- exp(log_returns) * rep(theta, each = n_inner) # S_T = theta e^Z
  payoff <- abs(terminal - params$k) # straddle payoffs

  return(exp(-params$rf * params$t_in) * colMeans(payoff))
}

laguerre_basis <- function(x, degree = 2) {
  # Weighted Laguerre polynomials of degree 0, ..., degree (at most 3)
  # Inputs:
  # -- x: vector of points
  # -- degree: highest degree
  # Output:
  # -- data frame with columns L0, L1, ...
  # formula: L_n(x) = exp(-x/2) * (Laguerre polynomial of degree n)

  weight <- exp(-x / 2)
  basis <- data.frame(L0 = weight)
  if (degree >= 1) basis$L1 <- weight * (1 - x)
  if (degree >= 2) basis$L2 <- weight * (1 - 2 * x + x^2 / 2)
  if (degree >= 3) basis$L3 <- weight * (1 - 3 * x + 3 * x^2 / 2 - x^3 / 6)

  return(basis)
}

single_asset_regression <- function(theta, design_points, params) {
  # Regression design (Longstaff-Schwartz type): one inner replication at each
  # design point, least-squares fit of the outputs on the Laguerre basis of
  # theta / s0, prediction at the scenarios theta
  # Inputs:
  # -- theta: vector of M stock prices at tau
  # -- design_points: stock prices at which the replications are run
  # -- params: single_asset_params()
  # Output:
  # -- vector of M estimates of mu(theta_i)

  y <- single_asset_standard(design_points, 1, params) # one rep per point
  training <- laguerre_basis(design_points / params$s0)
  training$y <- y
  fit <- lm(y ~ L0 + L1 + L2, data = training)
  prediction_basis <- laguerre_basis(theta / params$s0)

  return(unname(predict(fit, newdata = prediction_basis)))
}

single_asset_matched_sizes <- function(budget) {
  # Sizes of the comparison designs that use the optimal design's budget
  # Input:
  # -- budget: Gamma, the optimal design's budget
  # Output: a list of 3 items
  # -- sns_outer, sns_inner: SNS outer scenarios and inner replications
  # -- regression_points: number of regression design points

  return(list(sns_outer = ceiling(budget^(2 / 3)),
              sns_inner = ceiling(budget^(1 / 3)),
              regression_points = ceiling(budget)))
}

# ===== True risk measures
single_asset_true_risks <- function(alpha, threshold, params, n_true = 1e8,
                                    chunk = 1e7) {
  # True risk measures from n_true equally spaced quantiles of S_tau
  # Inputs:
  # -- alpha, threshold: as in risk_measures()
  # -- params: single_asset_params()
  # -- n_true: number of quantiles (10^8 in the paper)
  # -- chunk: number of conditional means computed at a time (memory)
  # Output:
  # -- named vector of the four true risk measures

  mu_true <- numeric(n_true) # placeholder for the true conditional means
  for (start in seq(1, n_true, by = chunk)) {
    i <- start:min(start + chunk - 1, n_true)
    theta <- qlnorm(i / (n_true + 1), log(params$s0) + params$meanlog_out,
                    params$sdlog_out)
    mu_true[i] <- single_asset_true_value(theta, params)
  }

  return(risk_measures(mu_true, alpha, threshold))
}
