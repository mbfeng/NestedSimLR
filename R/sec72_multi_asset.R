################################################################################
# sec72_multi_asset.R                                                          #
# Section 7.2: multi-asset ERM example                                         #
# 10 correlated stocks in the multivariate Black-Scholes model; on EACH stock  #
# the portfolio holds 30 options                                               #
# -- 10 European options: 5 puts and 5 calls                                   #
# -- 10 geometric average Asian options: 5 puts and 5 calls                    #
# -- 10 down-and-out put barrier options                                       #
# Weekly time steps (h = 1/52), risk horizon tau = 4 weeks, maturity T = 2     #
# years (100 weekly steps after tau)                                           #
# -- outer scenario theta = the 10 stock prices at tau (one row of theta)      #
# -- inner simulation = 100 weekly steps of every stock after tau              #
# -- conditional mean mu(theta) = portfolio value at tau (closed form)         #
# The optimal design is applied to each stock separately: the portfolio value  #
# is a sum over stocks, and the inner simulation of a stock depends only on    #
# that stock's price at tau. The LR input is the stock's price after the FIRST #
# inner step; all later steps have the same distribution in every scenario     #
# and cancel in the LR.                                                        #
################################################################################

# ===== Model
multi_asset_params <- function() {
  # Model, portfolio and experiment parameters of Section 7.2
  # Output:
  # -- list of parameters
  # The correlation matrix is drawn with randcorr() right after
  # set.seed(2023) with R's default (Mersenne-Twister) generator, exactly as in
  # the original code, so it is the matrix used in the paper and the same in
  # all macro runs. The generator type in use before the call is restored.

  params <- list(
    n_stocks = 10, # number of stocks
    s0 = rep(100, 10), # initial prices
    rf = 0.02, # risk-free rate
    big_t = 2, # maturity T (years)
    tau = 4 / 52, # risk horizon (4 weeks)
    h = 1 / 52 # time step (1 week)
  )
  params$mu <- seq(0.05, 0.14, length.out = params$n_stocks) # drifts
  params$vol <- seq(0.30, 0.48, length.out = params$n_stocks) # volatilities

  #-- random correlation matrix, as in the original code
  old_generator <- RNGkind() # generator in use before this call
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(2023)
  correlation <- randcorr::randcorr(params$n_stocks)
  RNGkind(old_generator[1], old_generator[2], old_generator[3])
  params$sigma <- diag(params$vol) %*% correlation %*% diag(params$vol)

  #-- time steps and per-step log-return distributions
  params$t_in <- params$big_t - params$tau
  params$n_steps_out <- round(params$tau / params$h) # 4 steps up to tau
  params$n_steps_in <- round(params$t_in / params$h) # 100 steps after tau
  params$meanlog_out <- (params$mu - params$vol^2 / 2) * params$h # real world
  params$meanlog_in <- (params$rf - params$vol^2 / 2) * params$h # risk neutral
  params$sdlog_step <- params$vol * sqrt(params$h)

  #-- options held on EACH stock (all long positions)
  params$k_euro <- seq(85, 130, by = 5) # European strikes
  params$type_euro <- rep(c("p", "c"), each = 5) # 5 puts, then 5 calls
  params$k_asian <- seq(85, 130, by = 5) # geometric Asian strikes
  params$type_asian <- rep(c("p", "c"), each = 5)
  params$k_barrier <- seq(85, 130, by = 5) # down-and-out put strikes
  params$h_barrier <- seq(75, 120, by = 5) # and their barriers

  return(params)
}

multi_asset_scenarios <- function(n, params) {
  # n outer scenarios
  # Inputs:
  # -- n: number of scenarios
  # -- params: multi_asset_params()
  # Output:
  # -- n-by-10 matrix of stock prices at tau

  log_returns <- mvtnorm::rmvnorm(n, mean = params$meanlog_out * params$n_steps_out,
                                  sigma = params$sigma * params$tau)
  initial_prices <- matrix(params$s0, n, params$n_stocks, byrow = TRUE)

  return(exp(log_returns) * initial_prices)
}

multi_asset_true_value_stock <- function(theta_stock, stock, params) {
  # True value at tau of the 30 options on one stock
  # Inputs:
  # -- theta_stock: vector of the stock's prices at tau
  # -- stock: index of the stock
  # -- params: multi_asset_params()
  # Output:
  # -- vector of option values, one per price

  vol <- params$vol[stock] # the stock's volatility
  value <- 0 # running sum over the options
  for (k in seq_along(params$k_euro)) {
    value <- value + bs_price(theta_stock, params$k_euro[k], params$rf, vol,
                              params$t_in, params$type_euro[k])
  }
  for (k in seq_along(params$k_asian)) {
    value <- value + geo_asian_price(theta_stock, params$k_asian[k], params$rf,
                                     vol, params$t_in, params$n_steps_in,
                                     params$type_asian[k])
  }
  for (k in seq_along(params$k_barrier)) {
    value <- value + barrier_dop_price(theta_stock, params$k_barrier[k],
                                       params$h_barrier[k], params$rf, 0, vol,
                                       params$h, params$t_in)
  }

  return(value)
}

multi_asset_true_value <- function(theta, params) {
  # True portfolio value mu(theta) for every scenario (row of theta)
  # Inputs:
  # -- theta: n-by-10 matrix of stock prices at tau
  # -- params: multi_asset_params()
  # Output:
  # -- vector of n portfolio values

  value_by_stock <- matrix(0, nrow(theta), params$n_stocks) # placeholder
  for (stock in 1:params$n_stocks) {
    value_by_stock[, stock] <- multi_asset_true_value_stock(theta[, stock],
                                                            stock, params)
  }

  return(rowSums(value_by_stock))
}

option_payoff <- function(x, k, type) {
  # Payoff of a call or put
  # Inputs:
  # -- x: price(s) of the underlying
  # -- k: strike
  # -- type: "c" = call, "p" = put
  # Output:
  # -- payoff(s)

  if (type == "c") {
    return(pmax(x - k, 0))
  } else {
    return(pmax(k - x, 0))
  }
}

# ===== Simulation
multi_asset_inner <- function(start_prices, stock, params) {
  # Inner simulation of ONE stock: one path per element of start_prices (the
  # stock's price at tau on that path), 100 weekly steps each
  # Inputs:
  # -- start_prices: vector of prices at tau, one per path
  # -- stock: index of the stock
  # -- params: multi_asset_params()
  # Output: a list of 2 items
  # -- y: discounted payoff of the stock's 30 options on every path
  # -- log_x: log price after the first step, the LR input

  n_paths <- length(start_prices) # number of paths
  n_steps <- params$n_steps_in # number of weekly steps
  log_returns <- matrix(rnorm(n_steps * n_paths, params$meanlog_in[stock],
                              params$sdlog_step[stock]), nrow = n_steps)
  #-- log_paths[step, path] = log price after `step` steps
  log_start <- matrix(log(start_prices), n_steps, n_paths, byrow = TRUE)
  log_paths <- matrixStats::colCumsums(log_returns) + log_start
  terminal <- exp(log_paths[n_steps, ]) # S_T
  geometric_average <- exp(colMeans(log_paths)) # over the 100 weekly dates

  #-- European and Asian options
  payoff <- numeric(n_paths) # running sum over the options
  for (k in seq_along(params$k_euro)) {
    payoff <- payoff + option_payoff(terminal, params$k_euro[k],
                                     params$type_euro[k])
  }
  for (k in seq_along(params$k_asian)) {
    payoff <- payoff + option_payoff(geometric_average, params$k_asian[k],
                                     params$type_asian[k])
  }

  #-- down-and-out puts: pay (K - S_T)^+ if the path never fell below the
  #-- barrier; survival is only checked where (K - S_T)^+ > 0
  n_barrier <- length(params$k_barrier) # number of barrier options
  put_payoff <- matrix(0, n_paths, n_barrier) # placeholder
  for (k in 1:n_barrier) {
    put_payoff[, k] <- pmax(params$k_barrier[k] - terminal, 0)
  }
  survived <- barrier_survive(log_paths, params$h_barrier, params$vol[stock],
                              params$h, needed = put_payoff > 0)
  survived[is.na(survived)] <- FALSE # not evaluated; the payoff is 0 there
  payoff <- payoff + rowSums(put_payoff * survived)

  discount <- exp(-params$rf * params$t_in) # discount factor from T to tau

  return(list(y = discount * payoff, log_x = log_paths[1, ]))
}

multi_asset_standard <- function(theta, n_inner, params, batch_size = 10) {
  # Standard nested simulation: the average of n_inner portfolio replications
  # at every scenario (row of theta)
  # Inputs:
  # -- theta: n-by-10 matrix of stock prices at tau
  # -- n_inner: number of replications per scenario
  # -- params: multi_asset_params()
  # -- batch_size: scenarios simulated at a time (memory)
  # Output:
  # -- vector of n estimates of mu(theta_i)

  n_outer <- nrow(theta) # number of scenarios
  estimate <- numeric(n_outer) # placeholder, running sum over the stocks
  for (first in seq(1, n_outer, by = batch_size)) {
    rows <- first:min(first + batch_size - 1, n_outer) # scenarios in batch
    for (stock in 1:params$n_stocks) {
      start_prices <- rep(theta[rows, stock], each = n_inner)
      y <- multi_asset_inner(start_prices, stock, params)$y
      y_by_scenario <- matrix(y, nrow = n_inner) # column = scenario
      estimate[rows] <- estimate[rows] + colMeans(y_by_scenario)
    }
  }

  return(estimate)
}

# ===== Designs
multi_asset_design <- function(theta_stock, stock, target_ess, params) {
  # The optimal design for one stock (Section 7.2 LP set-up)
  # Inputs:
  # -- theta_stock: vector of the stock's prices at tau
  # -- stock: index of the stock
  # -- target_ess: N
  # -- params: multi_asset_params()
  # Output:
  # -- design, see optimal_design()
  # The LR input (log price after the first inner step) has log-mean
  # log(theta) + meanlog_in and log-sd sdlog_step.

  log_means <- log(theta_stock) + params$meanlog_in[stock]
  psi <- psi_lognormal(log_means, params$sdlog_step[stock])

  return(optimal_design(psi, target_ess, section = "7.2"))
}

multi_asset_optimal <- function(theta, target_ess, params) {
  # One macro run of the optimal design, stock by stock
  # Inputs:
  # -- theta: n-by-10 matrix of stock prices at tau
  # -- target_ess: N
  # -- params: multi_asset_params()
  # Output: a list of 5 items
  # -- estimate: vector of estimates of mu(theta_i)
  # -- budget, n_sampled: per-stock budgets and numbers of sampled scenarios
  # -- seconds_lp, seconds_simulation: time spent solving LPs and simulating

  n_outer <- nrow(theta) # number of scenarios
  estimate <- numeric(n_outer) # placeholder, running sum over the stocks
  budget <- numeric(params$n_stocks) # placeholder
  n_sampled <- numeric(params$n_stocks) # placeholder
  seconds_lp <- 0
  seconds_simulation <- 0
  for (stock in 1:params$n_stocks) {
    #-- the design for this stock
    start <- Sys.time()
    design <- multi_asset_design(theta[, stock], stock, target_ess, params)
    seconds_lp <- seconds_lp + seconds_since(start)
    budget[stock] <- design$budget
    n_sampled[stock] <- length(design$sampled)

    #-- inner replications at the sampled scenarios, pooled with LRs
    start <- Sys.time()
    log_means <- log(theta[, stock]) + params$meanlog_in[stock]
    inner_samples <- list() # placeholder, one element per sampled scenario
    for (b in seq_along(design$sampled)) {
      j <- design$sampled[b] # index of the b-th sampled scenario
      start_prices <- rep(theta[j, stock], design$n_inner[b])
      sim <- multi_asset_inner(start_prices, stock, params)
      log_lr <- log_lr_lognormal(sim$log_x, log_means, log_means[j],
                                 params$sdlog_step[stock])
      inner_samples[[b]] <- list(y = sim$y, log_lr = log_lr)
    }
    estimate <- estimate + pooled_estimate(inner_samples, design$weights)
    seconds_simulation <- seconds_simulation + seconds_since(start)
  }

  return(list(estimate = estimate, budget = budget, n_sampled = n_sampled,
              seconds_lp = seconds_lp, seconds_simulation = seconds_simulation))
}

multi_asset_laguerre <- function(theta, params, degree = 2) {
  # Weighted Laguerre basis for all stocks, of the prices relative to s0
  # Inputs:
  # -- theta: n-by-10 matrix of stock prices
  # -- params: multi_asset_params()
  # -- degree: highest degree
  # Output:
  # -- data frame with columns L0_s1, L1_s1, L2_s1, L0_s2, ...

  basis <- NULL # placeholder
  for (stock in 1:params$n_stocks) {
    stock_basis <- laguerre_basis(theta[, stock] / params$s0[stock], degree)
    names(stock_basis) <- paste0(names(stock_basis), "_s", stock)
    if (is.null(basis)) {
      basis <- stock_basis
    } else {
      basis <- cbind(basis, stock_basis)
    }
  }

  return(basis)
}

multi_asset_regression <- function(theta, n_design, params, degrees = 2) {
  # Regression design: one inner replication of the whole portfolio at each
  # of n_design new outer scenarios, one least-squares regression on the
  # Laguerre basis of all stocks (with degree 2: 3 x 10 + intercept = 31
  # coefficients), predicted at the scenarios theta
  # Inputs:
  # -- theta: n-by-10 matrix of stock prices at tau
  # -- n_design: number of design points
  # -- params: multi_asset_params()
  # -- degrees: degree(s) of the basis; several degrees are fitted to the SAME
  #    simulated data (a sensitivity check of the regression basis)
  # Output: a list of 3 items
  # -- estimates: list, one vector per degree, named "deg2", "deg3", ...
  # -- seconds_simulation: simulation time
  # -- seconds_fit: fitting time per degree

  #-- one replication at each design point
  start <- Sys.time()
  design_points <- multi_asset_scenarios(n_design, params)
  y <- numeric(n_design) # placeholder, running sum over the stocks
  for (stock in 1:params$n_stocks) {
    y <- y + multi_asset_inner(design_points[, stock], stock, params)$y
  }
  seconds_simulation <- seconds_since(start)

  #-- least squares for each degree
  estimates <- list() # placeholder
  seconds_fit <- numeric(0) # placeholder
  for (degree in degrees) {
    start <- Sys.time()
    training <- multi_asset_laguerre(design_points, params, degree)
    training$y <- y
    fit <- lm(y ~ ., data = training) # y on all basis columns + intercept
    prediction_basis <- multi_asset_laguerre(theta, params, degree)
    name <- paste0("deg", degree)
    estimates[[name]] <- unname(predict(fit, newdata = prediction_basis))
    seconds_fit[name] <- seconds_since(start)
  }

  return(list(estimates = estimates, seconds_simulation = seconds_simulation,
              seconds_fit = seconds_fit))
}

multi_asset_sns <- function(budget, params) {
  # SNS with the optimal design's budget (average per-stock budget)
  # Inputs:
  # -- budget: Gamma
  # -- params: multi_asset_params()
  # Output: a list of 5 items
  # -- estimate, truth: estimates and true values at its own scenarios
  # -- n_outer, n_inner: ceiling(Gamma^(2/3)) NEW scenarios x
  #    ceiling(Gamma^(1/3)) inner replications
  # -- seconds: simulation time

  n_outer <- ceiling(budget^(2 / 3))
  n_inner <- ceiling(budget^(1 / 3))
  theta <- multi_asset_scenarios(n_outer, params)
  truth <- multi_asset_true_value(theta, params)
  start <- Sys.time()
  estimate <- multi_asset_standard(theta, n_inner, params)

  return(list(estimate = estimate, truth = truth, n_outer = n_outer,
              n_inner = n_inner, seconds = seconds_since(start)))
}

multi_asset_stats <- function(x, alpha = 0.995, threshold = 2800) {
  # The statistics of Table 3: the mean of mu(theta) and the four risk
  # measures (VaR level 0.995 and threshold 2800, as in the original code)
  # Inputs:
  # -- x: vector of conditional means
  # -- alpha, threshold: as in risk_measures()
  # Output:
  # -- named vector: mean, VaR, indicator, hockeystick, square

  return(c(mean = mean(x), risk_measures(x, alpha, threshold)))
}
