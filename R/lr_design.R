################################################################################
# lr_design.R
#
# The optimal likelihood-ratio (LR) design of the paper (Section 4) and the
# pooled self-normalized LR estimator (Section 3). Used by all three examples.
#
# Paper notation -> code:
#   outer scenarios theta_1, ..., theta_M    -> rows/elements of `theta`
#   target effective sample size N           -> target_ess
#   Psi[i, j] = 1 / E_j[W_ij^2]              -> psi[i, j]   (i = target,
#                                                            j = sampling)
#   LP (11):  minimize   sum_j c_j
#             subject to sum_j Psi[i, j] c_j >= 1 for every i, c_j >= 0
#   inner sample size at scenario j          -> n_inner = ceiling(N c_j)
#   pooling weights gamma_ij (Eq. (10)),
#     proportional to n_j Psi[i, j]          -> design$weights[i, b]
#   self-normalized LR estimate of mu(theta_i)
#     from the samples at theta_j            -> self_normalized_estimate()
#   pooled estimate sum_j gamma_ij mu~_ij    -> pooled_estimate()
#
# "Sampled scenarios" are the scenarios j with c_j > 0: only these are
# simulated. In the code they are numbered b = 1, ..., J, and
# design$sampled[b] is the index j of the b-th sampled scenario.
################################################################################

#--- Psi matrices ---------------------------------------------------------------

# Psi for lognormal inputs with a common log standard deviation.
#
# If the LR input is lognormal with log-mean a_i under scenario i and log-sd
# s under every scenario, then E_j[W_ij^2] = exp((a_i - a_j)^2 / s^2), so
# Psi[i, j] = exp(-(a_i - a_j)^2 / s^2).
#
# log_means: vector (a_1, ..., a_M); log_sd: s.
psi_lognormal <- function(log_means, log_sd) {
  differences <- outer(log_means, log_means, "-") # [i, j] = a_i - a_j
  exp(-differences^2 / log_sd^2)
}

# Psi for independent Poisson inputs (the newsvendor example).
#
# For one Poisson input with mean l_t under the target scenario and l_s under
# the sampling scenario, E_s[W^2] = exp(l_t^2 / l_s + l_s - 2 l_t). With
# independent inputs the second moments multiply, so their logs add up.
#
# lambda_target: (number of targets) x (number of inputs) matrix of means.
# lambda_sampling: same for the sampling scenarios (default: the targets).
psi_poisson <- function(lambda_target, lambda_sampling = lambda_target) {
  n_target <- nrow(lambda_target)
  n_sampling <- nrow(lambda_sampling)
  log_second_moment <- matrix(0, n_target, n_sampling)
  for (input in seq_len(ncol(lambda_target))) {
    # [i, j] entries: l_t = target mean of scenario i, l_s = sampling mean of
    # scenario j (whole matrices at once, for speed with M = 10,000)
    l_t <- matrix(lambda_target[, input], n_target, n_sampling)
    l_s <- matrix(lambda_sampling[, input], n_target, n_sampling, byrow = TRUE)
    log_second_moment <- log_second_moment + (l_t^2 / l_s + (-2 * l_t + l_s))
  }
  exp(-log_second_moment)
}

#--- The optimal design: LP (11) ----------------------------------------------------

# Solve LP (11) with lpSolve (used in Sections 7.1 and 7.3).
# Returns the vector c = (c_1, ..., c_M).
solve_lp_lpsolve <- function(psi) {
  m <- nrow(psi)
  solution <- lpSolve::lp(
    direction = "min",
    objective.in = rep(1, m),       # minimize sum_j c_j
    const.mat = psi,                # constraint i: sum_j psi[i, j] c_j >= 1
    const.dir = rep(">=", m),
    const.rhs = rep(1, m)
  )
  if (solution$status != 0) {
    stop("lpSolve failed with status ", solution$status)
  }
  solution$solution
}

# Solve LP (11) with lpSolveAPI and the bounds 0 <= c_j <= 1 (used in
# Section 7.2). With tuned = TRUE the solver settings of the original code
# are used (dual simplex, devex pricing, Curtis-Reid scaling).
# Returns list(status, c); status 0 means success.
#
# Note: the original code created the model with make.lp(m, m), i.e., with m
# empty rows before adding the m constraints. This is kept on purpose: the
# LP can have several optimal solutions, and keeping the same set-up makes
# the solver find the same one as in the paper's experiments.
solve_lp_lpsolveapi <- function(psi, tuned = TRUE) {
  m <- nrow(psi)
  model <- lpSolveAPI::make.lp(m, m)
  lpSolveAPI::set.objfn(model, rep(1, m))
  for (i in seq_len(m)) {
    lpSolveAPI::add.constraint(model, psi[i, ], ">=", 1)
  }
  lpSolveAPI::set.bounds(model, lower = rep(0, m), upper = rep(1, m))
  if (tuned) {
    lpSolveAPI::lp.control(model, pivoting = "devex", simplextype = "dual",
                           scaling = "curtisreid")
  }
  status <- solve(model)
  c_values <- NULL
  if (status == 0) {
    c_values <- lpSolveAPI::get.variables(model)
  }
  list(status = status, c = c_values)
}

# Section 7.2 version: lpSolveAPI with the original settings, with two
# fall-backs. The tuned settings occasionally fail numerically (status 5):
# when two outer scenarios almost coincide, Psi has two nearly identical
# columns. (The original code then stopped the whole batch of macro runs.)
# The fall-backs solve the same LP: first lpSolveAPI with its default
# settings, then lpSolve with geometric scaling and c_j <= 1 as extra rows.
solve_lp_bounded <- function(psi) {
  result <- solve_lp_lpsolveapi(psi, tuned = TRUE)
  if (result$status == 0) {
    return(result$c)
  }
  message("LP: lpSolveAPI status ", result$status,
          " with the tuned settings; retrying with the default settings")
  result <- solve_lp_lpsolveapi(psi, tuned = FALSE)
  if (result$status == 0) {
    return(result$c)
  }
  message("LP: retrying with lpSolve")
  m <- nrow(psi)
  solution <- lpSolve::lp(
    direction = "min",
    objective.in = rep(1, m),
    const.mat = rbind(psi, diag(m)),             # Psi c >= 1 and c <= 1
    const.dir = c(rep(">=", m), rep("<=", m)),
    const.rhs = rep(1, 2 * m),
    scale = 4                                    # geometric scaling
  )
  if (solution$status != 0) {
    stop("all LP solvers failed")
  }
  solution$solution
}

# The optimal design: solve LP (11), choose the sampled scenarios, their
# inner sample sizes, and the pooling weights.
#
# Arguments:
#   psi         M x M matrix Psi
#   target_ess  N, the target effective sample size
#   section     "7.1/7.3" (plain LP, lpSolve) or "7.2" (the Section 7.2
#               set-up of the original code, see below)
#
# Section 7.2 set-up: Psi entries below 1e-6 are set to 0 (numerical
# stability), bounds 0 <= c_j <= 1, every scenario with c_j > 0 is sampled,
# and every sampled scenario gets at least 10 inner replications.
# Sections 7.1 and 7.3: no bounds, scenarios with c_j > 1e-8 are sampled.
#
# Returns a list:
#   sampled   indices j of the sampled scenarios (length J)
#   n_inner   their inner sample sizes n_j (length J)
#   budget    total number of inner replications, sum(n_inner)
#   weights   M x J matrix of pooling weights gamma_ij (each row sums to 1)
optimal_design <- function(psi, target_ess, section = c("7.1/7.3", "7.2")) {
  section <- match.arg(section)
  if (section == "7.2") {
    psi[psi < 1e-6] <- 0
    c_values <- solve_lp_bounded(psi)
    sampled <- which(c_values > 0)
    min_reps <- 10
  } else {
    c_values <- solve_lp_lpsolve(psi)
    sampled <- which(c_values > 1e-8)
    min_reps <- 0
  }
  n_inner <- pmax(ceiling(c_values[sampled] * target_ess), min_reps)

  # gamma_ij proportional to n_j Psi[i, j], normalized so that each row sums
  # to 1 (Eq. (10), with the rounded-up n_j)
  weights <- psi[, sampled, drop = FALSE]
  for (b in seq_along(sampled)) {
    weights[, b] <- weights[, b] * n_inner[b]
  }
  weights <- weights / rowSums(weights)

  list(sampled = sampled, n_inner = n_inner, budget = sum(n_inner),
       weights = weights)
}

#--- The pooled self-normalized LR estimator -------------------------------------

# Self-normalized LR estimates of mu(theta_i) for every target i, from the
# n inner replications of ONE sampled scenario j:
#   mu~_ij = sum_k W_ijk y_k / sum_k W_ijk.
#
# y: the n outputs; log_lr: M x n matrix of log W_ijk.
#
# Each row of log_lr is shifted by its maximum before exponentiating. The
# estimate does not change (the factor cancels in the ratio), and this
# prevents overflow and underflow.
self_normalized_estimate <- function(y, log_lr) {
  log_lr <- log_lr - matrixStats::rowMaxs(log_lr)
  lr <- exp(log_lr)
  as.vector(lr %*% y) / rowSums(lr)
}

# Pooled estimate mu^(theta_i) = sum_j gamma_ij mu~_ij (Eq. (8)).
#
# inner_samples: list with one element per sampled scenario (in the order of
#   design$sampled); each element is list(y = outputs, log_lr = M x n matrix).
# weights: design$weights.
pooled_estimate <- function(inner_samples, weights) {
  estimates <- matrix(0, nrow(weights), ncol(weights)) # [i, b] = mu~_ij
  for (b in seq_along(inner_samples)) {
    estimates[, b] <- self_normalized_estimate(inner_samples[[b]]$y,
                                               inner_samples[[b]]$log_lr)
  }
  rowSums(estimates * weights)
}

#--- Log likelihood ratios ----------------------------------------------------------

# log W_ij(x) for lognormal inputs with a common log-sd: the log of the
# lognormal density with log-mean a_i minus that with log-mean a_j.
#
# log_x: the n log inputs simulated at the sampling scenario j.
# target_log_means: (a_1, ..., a_M); sampling_log_mean: a_j.
# Returns the M x n matrix of log W_ijk.
#
# Uses (x - a_j)^2 - (x - a_i)^2 = (a_i - a_j)(2x - a_i - a_j).
log_lr_lognormal <- function(log_x, target_log_means, sampling_log_mean,
                             log_sd) {
  d <- target_log_means - sampling_log_mean # a_i - a_j, length M
  (outer(d, 2 * log_x) - d * (target_log_means + sampling_log_mean)) /
    (2 * log_sd^2)
}

# log W_ij(x) for independent Poisson inputs:
#   log W_ijk = sum_l [ x_kl log(lambda_il / lambda_jl) - (lambda_il - lambda_jl) ].
#
# x: n x d matrix of inputs simulated at the sampling scenario j.
# lambda_target: M x d matrix; lambda_sampling: length-d vector (scenario j).
# Returns the M x n matrix of log W_ijk.
log_lr_poisson <- function(x, lambda_target, lambda_sampling) {
  n_target <- nrow(lambda_target)
  log_sampling <- matrix(log(lambda_sampling), n_target, length(lambda_sampling),
                         byrow = TRUE)
  log_ratio <- log(lambda_target) - log_sampling            # M x d
  shift <- rowSums(lambda_target) - sum(lambda_sampling)    # length M
  log_ratio %*% t(x) - shift
}
