################################################################################
# sec73_newsvendor.R
#
# Section 7.3: input uncertainty quantification for a 10-product newsvendor.
#
# The daily demand X_l of product l is Poisson(lambda_l). The means lambda_l
# are unknown; with a Gamma(0.001, 0.001) prior and n_l observed demands, the
# posterior of lambda_l is Gamma(shape = 0.001 + sum of the observations,
# rate = 0.001 + n_l).
#   outer scenario theta  = lambda = (lambda_1, ..., lambda_10), a posterior
#                           draw (one row of the matrix `lambda`)
#   inner simulation      = demands X_1, ..., X_10 given lambda
#   output                = one-day profit
#                           g(X) = sum_l price_l min(X_l, k_l) - cost sum_l k_l
#   conditional mean      = mu(lambda), the expected profit (closed form)
#   nested statistic      = credible interval (CrI) of mu(lambda)
################################################################################

# Problem parameters and the posterior of lambda.
#
# The observed demands are simulated after set.seed(2020) with R's default
# (Mersenne-Twister) generator, exactly as in the original code, so the
# posterior is the one used in the paper. The generator is deliberately left
# in its state after these draws: the original scripts continued drawing
# from the same stream (the 10^6 posterior draws for Table 4 and the fixed
# scenarios of Figure 3(a)), and the scripts do the same.
newsvendor_setup <- function() {
  product <- 1:10
  params <- list(
    n_products = 10,
    mean_demand = 5 + product, # true demand means (only used to simulate data)
    k = 9 + product,           # stocking levels
    price = 7 + 3 * product,   # unit sale prices
    cost = 2,                  # unit cost (all products)
    n_obs = 50 + 5 * product   # number of observed demands
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(2020)
  prior_rate <- 0.001
  prior_shape <- 0.001
  params$posterior_rate <- numeric(params$n_products)
  params$posterior_shape <- numeric(params$n_products)
  for (l in product) {
    observed <- rpois(params$n_obs[l], params$mean_demand[l])
    params$posterior_rate[l] <- prior_rate + params$n_obs[l]
    params$posterior_shape[l] <- prior_shape + sum(observed)
  }
  params
}

# n posterior draws of lambda: an n x 10 matrix. Drawn product by product,
# in the same order as the original code.
newsvendor_draw_lambda <- function(n, params) {
  lambda <- matrix(0, n, params$n_products)
  for (l in seq_len(params$n_products)) {
    lambda[, l] <- rgamma(n, rate = params$posterior_rate[l],
                          shape = params$posterior_shape[l])
  }
  lambda
}

# The probabilities P(X = values[v]) of a Poisson(lambda[i]) variable, for
# every lambda (rows) and value (columns).
poisson_pmf_matrix <- function(lambda, values) {
  pmf <- matrix(0, length(lambda), length(values))
  for (v in seq_along(values)) {
    pmf[, v] <- dpois(values[v], lambda)
  }
  pmf
}

# Exact conditional mean mu(lambda) (expected profit) for every row of lambda:
# E[min(X, k)] = sum_{j=1}^{k} j P(X = j) + k P(X > k).
newsvendor_true_profit <- function(lambda, params) {
  total <- -params$cost * sum(params$k)
  for (l in seq_len(params$n_products)) {
    values <- seq_len(params$k[l]) # 1, ..., k
    pmf <- poisson_pmf_matrix(lambda[, l], values)
    tail_prob <- ppois(params$k[l], lambda[, l], lower.tail = FALSE)
    expected_sold <- as.vector(pmf %*% values) + params$k[l] * tail_prob
    total <- total + params$price[l] * expected_sold
  }
  total
}

# Exact conditional variance of the profit for every row of lambda (used for
# the exact variance of SNS in Figure 3(a)).
newsvendor_true_variance <- function(lambda, params) {
  total <- 0
  for (l in seq_len(params$n_products)) {
    values <- 0:params$k[l]
    pmf <- poisson_pmf_matrix(lambda[, l], values)
    tail_prob <- 1 - rowSums(pmf) # P(X > k)
    first_moment <- as.vector(pmf %*% values) + params$k[l] * tail_prob
    second_moment <- as.vector(pmf %*% values^2) + params$k[l]^2 * tail_prob
    total <- total + params$price[l]^2 * (second_moment - first_moment^2)
  }
  total
}

# One-day profit for every row of a demand matrix (rows = replications,
# columns = products): sum_l price_l min(X_l, k_l) - cost sum_l k_l.
newsvendor_profit <- function(demand, params) {
  sold <- demand
  for (l in seq_len(params$n_products)) {
    sold[, l] <- pmin(demand[, l], params$k[l])
  }
  as.vector(sold %*% params$price) - params$cost * sum(params$k)
}

# Poisson demands for many replications at once: replication r has means
# lambda_rows[r, ] (one row per replication). Returns a matrix with one row
# per replication and one column per product.
newsvendor_demand <- function(lambda_rows, params) {
  n <- nrow(lambda_rows)
  demand <- matrix(0, n, params$n_products)
  for (l in seq_len(params$n_products)) {
    demand[, l] <- rpois(n, lambda_rows[, l]) # rpois accepts a vector of means
  }
  demand
}

# n inner replications at ONE scenario lambda_j (vector of 10 means).
# Returns the profits y and the demands x (n x 10, the LR input).
newsvendor_inner <- function(lambda_j, n, params) {
  lambda_rows <- matrix(lambda_j, n, params$n_products, byrow = TRUE)
  demand <- newsvendor_demand(lambda_rows, params)
  list(y = newsvendor_profit(demand, params), x = demand)
}

# Standard nested simulation: average of n_inner replications at every
# scenario (row of lambda).
#
# The loop is over scenarios; the n_inner replications of a scenario are
# simulated as one vector. Simulating all scenarios at once was tried and is
# not faster (the time goes into drawing the Poisson numbers), so the loop is
# kept: it draws the same random numbers as the 2024 code.
newsvendor_standard <- function(lambda, n_inner, params) {
  estimate <- numeric(nrow(lambda))
  for (i in seq_len(nrow(lambda))) {
    estimate[i] <- mean(newsvendor_inner(lambda[i, ], n_inner, params)$y)
  }
  estimate
}

# The optimal design for the scenarios lambda (deterministic).
newsvendor_design <- function(lambda, target_ess) {
  optimal_design(psi_poisson(lambda), target_ess, section = "7.1/7.3")
}

# One macro run of the optimal design: estimates of mu(lambda_i) for all
# scenarios, from inner replications at the sampled scenarios only.
newsvendor_optimal <- function(lambda, design, params) {
  inner_samples <- list()
  for (b in seq_along(design$sampled)) {
    j <- design$sampled[b]
    sim <- newsvendor_inner(lambda[j, ], design$n_inner[b], params)
    log_lr <- log_lr_poisson(sim$x, lambda, lambda[j, ])
    inner_samples[[b]] <- list(y = sim$y, log_lr = log_lr)
  }
  pooled_estimate(inner_samples, design$weights)
}

# Regression basis: lambda_l and lambda_l^2 for each product (quadratic
# polynomial without cross terms).
newsvendor_regression_basis <- function(lambda) {
  basis <- data.frame(lambda, lambda^2)
  names(basis) <- c(paste0("x", 1:10), paste0("x", 1:10, "_sq"))
  basis
}

# Regression design: one replication at each of n_design NEW posterior draws,
# least-squares fit on the quadratic basis, prediction at the scenarios lambda.
#
# Vectorized for speed: the replications at all design points are simulated
# at once, about 100 times faster than one design point at a time as in the
# 2024 code (one replication per call is dominated by call overhead). The
# random numbers are drawn in a different order than in 2024, with the same
# distribution. Changed November 2026; see README change log.
newsvendor_regression <- function(lambda, n_design, params) {
  design_points <- newsvendor_draw_lambda(n_design, params)
  demand <- newsvendor_demand(design_points, params) # one row per design point
  y <- newsvendor_profit(demand, params)
  training <- newsvendor_regression_basis(design_points)
  training$y <- y
  fit <- lm(y ~ ., data = training) # y on all basis columns + intercept
  unname(predict(fit, newdata = newsvendor_regression_basis(lambda)))
}

# Credible intervals from the estimated conditional means, and their coverage.
#
# For each level 1 - alpha, the CrI is [alpha/2, 1 - alpha/2] sample quantile
# of the M estimates (type-7 quantiles, as in the original code). Its
# coverage is the probability that the TRUE mu(lambda) of a posterior draw
# falls inside it, estimated with a large sorted sample `true_sorted` of true
# values (10^6 posterior draws in the paper).
#
# Returns a data frame: level, lower, upper, width, coverage.
newsvendor_cri <- function(estimates, true_sorted, levels = c(0.90, 0.95, 0.99)) {
  n_levels <- length(levels)
  probabilities <- c((1 - levels) / 2, 1 - (1 - levels) / 2)
  q <- quantile(estimates, probabilities, names = FALSE)
  lower <- q[seq_len(n_levels)]
  upper <- q[n_levels + seq_len(n_levels)]
  # number of true values in [lower, upper]: count of values <= upper minus
  # count of values < lower
  n_inside <- findInterval(upper, true_sorted) -
    findInterval(lower, true_sorted, left.open = TRUE)
  data.frame(level = levels, lower = lower, upper = upper,
             width = upper - lower, coverage = n_inside / length(true_sorted))
}
