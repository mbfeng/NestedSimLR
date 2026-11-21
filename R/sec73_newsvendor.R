################################################################################
# sec73_newsvendor.R                                                           #
# Section 7.3: input uncertainty quantification for a 10-product newsvendor    #
# The daily demand X_l of product l is Poisson(lambda_l). The means lambda_l   #
# are unknown; with a Gamma(0.001, 0.001) prior and n_l observed demands, the  #
# posterior of lambda_l is Gamma(shape = 0.001 + sum of the observations,      #
# rate = 0.001 + n_l).                                                         #
# -- outer scenario theta = lambda = (lambda_1, ..., lambda_10), a posterior   #
#    draw (one row of the matrix lambda)                                       #
# -- inner simulation = demands X_1, ..., X_10 given lambda                    #
# -- output = one-day profit                                                   #
#    g(X) = sum_l price_l min(X_l, k_l) - cost sum_l k_l                       #
# -- conditional mean = mu(lambda), the expected profit (closed form)          #
# -- nested statistic = credible interval (CrI) of mu(lambda)                  #
################################################################################

# ===== Model
newsvendor_setup <- function() {
  # Problem parameters and the posterior of lambda
  # Output:
  # -- list of parameters, including posterior_rate and posterior_shape
  # The observed demands are simulated after set.seed(2020) with R's default
  # (Mersenne-Twister) generator, exactly as in the original code, so the
  # posterior is the one used in the paper. The generator is deliberately
  # left in its state after these draws: the original scripts continued
  # drawing from the same stream (the 10^6 posterior draws for Table 4 and the
  # fixed scenarios of Figure 3(a)), and the scripts do the same.

  product <- 1:10
  params <- list(
    n_products = 10, # number of products
    mean_demand = 5 + product, # true demand means (only used to simulate data)
    k = 9 + product, # stocking levels
    price = 7 + 3 * product, # unit sale prices
    cost = 2, # unit cost (all products)
    n_obs = 50 + 5 * product # number of observed demands
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(2020)
  prior_rate <- 0.001 # Gamma(0.001, 0.001) prior
  prior_shape <- 0.001
  params$posterior_rate <- numeric(params$n_products) # placeholder
  params$posterior_shape <- numeric(params$n_products) # placeholder
  for (l in product) {
    observed <- rpois(params$n_obs[l], params$mean_demand[l]) # observed data
    params$posterior_rate[l] <- prior_rate + params$n_obs[l]
    params$posterior_shape[l] <- prior_shape + sum(observed)
  }

  return(params)
}

newsvendor_draw_lambda <- function(n, params) {
  # n posterior draws of lambda, product by product in the same order as the
  # original code
  # Inputs:
  # -- n: number of draws
  # -- params: newsvendor_setup()
  # Output:
  # -- n-by-10 matrix of Poisson means

  lambda <- matrix(0, n, params$n_products) # placeholder
  for (l in 1:params$n_products) {
    lambda[, l] <- rgamma(n, rate = params$posterior_rate[l],
                          shape = params$posterior_shape[l])
  }

  return(lambda)
}

poisson_pmf_matrix <- function(lambda, values) {
  # Poisson probabilities P(X = values[v]) for every mean lambda[i]
  # Inputs:
  # -- lambda: vector of Poisson means
  # -- values: vector of values
  # Output:
  # -- matrix, [i, v] = P(X = values[v]) when X ~ Poisson(lambda[i])

  pmf <- matrix(0, length(lambda), length(values)) # placeholder
  for (v in seq_along(values)) {
    pmf[, v] <- dpois(values[v], lambda)
  }

  return(pmf)
}

newsvendor_true_profit <- function(lambda, params) {
  # Exact conditional mean mu(lambda), the expected profit, for every row of
  # lambda
  # Inputs:
  # -- lambda: n-by-10 matrix of Poisson means
  # -- params: newsvendor_setup()
  # Output:
  # -- vector of n expected profits
  # formula: E[min(X, k)] = sum_{j=1}^{k} j P(X = j) + k P(X > k)

  total <- -params$cost * sum(params$k) # running sum over the products
  for (l in 1:params$n_products) {
    values <- 1:params$k[l]
    pmf <- poisson_pmf_matrix(lambda[, l], values)
    tail_prob <- ppois(params$k[l], lambda[, l], lower.tail = FALSE) # P(X > k)
    expected_sold <- as.vector(pmf %*% values) + params$k[l] * tail_prob
    total <- total + params$price[l] * expected_sold
  }

  return(total)
}

newsvendor_true_variance <- function(lambda, params) {
  # Exact conditional variance of the profit for every row of lambda (used for
  # the exact variance of SNS in Figure 3(a))
  # Inputs:
  # -- lambda: n-by-10 matrix of Poisson means
  # -- params: newsvendor_setup()
  # Output:
  # -- vector of n variances

  total <- 0 # running sum over the (independent) products
  for (l in 1:params$n_products) {
    values <- 0:params$k[l]
    pmf <- poisson_pmf_matrix(lambda[, l], values)
    tail_prob <- 1 - rowSums(pmf) # P(X > k)
    first_moment <- as.vector(pmf %*% values) + params$k[l] * tail_prob
    second_moment <- as.vector(pmf %*% values^2) + params$k[l]^2 * tail_prob
    total <- total + params$price[l]^2 * (second_moment - first_moment^2)
  }

  return(total)
}

# ===== Simulation
newsvendor_profit <- function(demand, params) {
  # One-day profit for every row of a demand matrix
  # Inputs:
  # -- demand: matrix, rows = replications, columns = products
  # -- params: newsvendor_setup()
  # Output:
  # -- vector of profits, one per replication
  # formula: sum_l price_l min(X_l, k_l) - cost sum_l k_l

  sold <- demand # placeholder for the units sold
  for (l in 1:params$n_products) {
    sold[, l] <- pmin(demand[, l], params$k[l])
  }

  return(as.vector(sold %*% params$price) - params$cost * sum(params$k))
}

newsvendor_demand <- function(lambda_rows, params) {
  # Poisson demands for many replications at once
  # Inputs:
  # -- lambda_rows: matrix, one row of means per replication
  # -- params: newsvendor_setup()
  # Output:
  # -- matrix, one row per replication, one column per product

  n <- nrow(lambda_rows) # number of replications
  demand <- matrix(0, n, params$n_products) # placeholder
  for (l in 1:params$n_products) {
    demand[, l] <- rpois(n, lambda_rows[, l]) # rpois accepts a vector of means
  }

  return(demand)
}

newsvendor_inner <- function(lambda_j, n, params) {
  # n inner replications at ONE scenario lambda_j
  # Inputs:
  # -- lambda_j: vector of the 10 Poisson means
  # -- n: number of replications
  # -- params: newsvendor_setup()
  # Output: a list of 2 items
  # -- y: profits
  # -- x: n-by-10 matrix of demands, the LR input

  lambda_rows <- matrix(lambda_j, n, params$n_products, byrow = TRUE)
  demand <- newsvendor_demand(lambda_rows, params)

  return(list(y = newsvendor_profit(demand, params), x = demand))
}

newsvendor_standard <- function(lambda, n_inner, params) {
  # Standard nested simulation: average of n_inner replications at every
  # scenario (row of lambda)
  # Inputs:
  # -- lambda: n-by-10 matrix of Poisson means
  # -- n_inner: number of replications per scenario
  # -- params: newsvendor_setup()
  # Output:
  # -- vector of n estimates of mu(lambda_i)
  # The loop is over scenarios; the n_inner replications of a scenario are
  # simulated as one vector. Simulating all scenarios at once was tried and
  # is not faster (the time goes into drawing the Poisson numbers), so the
  # loop is kept: it draws the same random numbers as the 2024 code.

  estimate <- numeric(nrow(lambda)) # placeholder
  for (i in 1:nrow(lambda)) {
    estimate[i] <- mean(newsvendor_inner(lambda[i, ], n_inner, params)$y)
  }

  return(estimate)
}

# ===== Designs
newsvendor_design <- function(lambda, target_ess) {
  # The optimal design for the scenarios lambda (deterministic)
  # Inputs:
  # -- lambda: M-by-10 matrix of Poisson means
  # -- target_ess: N
  # Output:
  # -- design, see optimal_design()

  return(optimal_design(psi_poisson(lambda), target_ess, section = "7.1/7.3"))
}

newsvendor_optimal <- function(lambda, design, params) {
  # One macro run of the optimal design
  # Inputs:
  # -- lambda: M-by-10 matrix of Poisson means
  # -- design: newsvendor_design()
  # -- params: newsvendor_setup()
  # Output:
  # -- vector of M estimates of mu(lambda_i), from inner replications at the
  #    sampled scenarios only

  inner_samples <- list() # placeholder, one element per sampled scenario
  for (b in seq_along(design$sampled)) {
    j <- design$sampled[b] # index of the b-th sampled scenario
    sim <- newsvendor_inner(lambda[j, ], design$n_inner[b], params)
    log_lr <- log_lr_poisson(sim$x, lambda, lambda[j, ])
    inner_samples[[b]] <- list(y = sim$y, log_lr = log_lr)
  }

  return(pooled_estimate(inner_samples, design$weights))
}

newsvendor_regression_basis <- function(lambda) {
  # Regression basis: lambda_l and lambda_l^2 for each product (quadratic
  # polynomial without cross terms)
  # Input:
  # -- lambda: n-by-10 matrix of Poisson means
  # Output:
  # -- data frame with columns x1, ..., x10, x1_sq, ..., x10_sq

  basis <- data.frame(lambda, lambda^2)
  names(basis) <- c(paste0("x", 1:10), paste0("x", 1:10, "_sq"))

  return(basis)
}

newsvendor_regression <- function(lambda, n_design, params) {
  # Regression design: one replication at each of n_design NEW posterior
  # draws, least-squares fit on the quadratic basis, prediction at lambda
  # Inputs:
  # -- lambda: M-by-10 matrix of Poisson means (the scenarios)
  # -- n_design: number of design points
  # -- params: newsvendor_setup()
  # Output:
  # -- vector of M estimates of mu(lambda_i)
  # Vectorized for speed: the replications at all design points are simulated
  # at once, about 100 times faster than one design point at a time as in the
  # 2024 code (one replication per call is dominated by call overhead). The
  # random numbers are drawn in a different order than in 2024, with the
  # same distribution. Changed November 2026; see README change log.

  design_points <- newsvendor_draw_lambda(n_design, params)
  demand <- newsvendor_demand(design_points, params) # one row per design point
  y <- newsvendor_profit(demand, params)
  training <- newsvendor_regression_basis(design_points)
  training$y <- y
  fit <- lm(y ~ ., data = training) # y on all basis columns + intercept

  return(unname(predict(fit, newdata = newsvendor_regression_basis(lambda))))
}

# ===== Credible intervals
newsvendor_cri <- function(estimates, true_sorted, levels = c(0.90, 0.95, 0.99)) {
  # Credible intervals from the estimated conditional means, and their
  # coverage
  # Inputs:
  # -- estimates: vector of M estimated conditional means
  # -- true_sorted: large sorted sample of TRUE mu(lambda) values (10^6
  #    posterior draws in the paper)
  # -- levels: CrI levels 1 - alpha
  # Output:
  # -- data frame: level, lower, upper, width, coverage
  # For each level, the CrI is the [alpha/2, 1 - alpha/2] sample quantile of
  # the M estimates (type-7 quantiles, as in the original code); its coverage
  # is the probability that the true mu(lambda) of a posterior draw falls
  # inside it.

  n_levels <- length(levels) # number of CrI levels
  probabilities <- c((1 - levels) / 2, 1 - (1 - levels) / 2)
  q <- quantile(estimates, probabilities, names = FALSE)
  lower <- q[1:n_levels]
  upper <- q[n_levels + (1:n_levels)]
  # number of true values in [lower, upper]: count of values <= upper minus
  # count of values < lower
  n_inside <- findInterval(upper, true_sorted) -
    findInterval(lower, true_sorted, left.open = TRUE)

  return(data.frame(level = levels, lower = lower, upper = upper,
                    width = upper - lower,
                    coverage = n_inside / length(true_sorted)))
}
