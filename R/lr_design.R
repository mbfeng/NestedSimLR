################################################################################
# lr_design.R                                                                  #
# The optimal likelihood-ratio (LR) design of the paper (Section 4) and the    #
# pooled self-normalized LR estimator (Section 3), used by all three examples  #
# Paper notation -> code:                                                      #
# -- outer scenarios theta_1, ..., theta_M -> rows/elements of theta           #
# -- target effective sample size N -> target_ess                              #
# -- Psi[i, j] = 1 / E_j[W_ij^2] (i = target, j = sampling) -> psi[i, j]       #
# -- LP (11): minimize sum_j c_j s.t. sum_j Psi[i, j] c_j >= 1, c_j >= 0       #
# -- inner sample size at scenario j -> n_inner = ceiling(N c_j)               #
# -- pooling weights gamma_ij (Eq. (10)) -> design$weights[i, b]               #
# -- self-normalized LR estimate mu~_ij -> self_normalized_estimate()          #
# -- pooled estimate sum_j gamma_ij mu~_ij -> pooled_estimate()                #
# "Sampled scenarios" are the scenarios j with c_j > 0; only these are         #
# simulated. They are numbered b = 1, ..., J, and design$sampled[b] is the     #
# index j of the b-th sampled scenario.                                        #
################################################################################

# ===== Psi matrices
psi_lognormal <- function(log_means, log_sd) {
  # Psi for lognormal inputs with a common log standard deviation
  # Inputs:
  # -- log_means: vector (a_1, ..., a_M), log-means of the LR input
  # -- log_sd: common log-sd s of the LR input
  # Output:
  # -- M-by-M matrix Psi
  # formula: E_j[W_ij^2] = exp((a_i - a_j)^2 / s^2), so
  #          Psi[i, j] = exp(-(a_i - a_j)^2 / s^2)

  differences <- outer(log_means, log_means, "-") # [i, j] = a_i - a_j

  return(exp(-differences^2 / log_sd^2))
}

psi_poisson <- function(lambda_target, lambda_sampling = lambda_target) {
  # Psi for independent Poisson inputs (the newsvendor example)
  # Inputs:
  # -- lambda_target: (number of targets)-by-(number of inputs) matrix of means
  # -- lambda_sampling: same for the sampling scenarios (default: the targets)
  # Output:
  # -- matrix Psi, [i, j] for target i and sampling scenario j
  # formula: for one input, E_s[W^2] = exp(l_t^2 / l_s + l_s - 2 l_t); with
  #          independent inputs the second moments multiply (logs add up)

  n_target <- nrow(lambda_target) # number of target scenarios
  n_sampling <- nrow(lambda_sampling) # number of sampling scenarios
  log_second_moment <- matrix(0, n_target, n_sampling) # placeholder
  for (input in 1:ncol(lambda_target)) {
    # [i, j] entries: l_t = target mean of scenario i, l_s = sampling mean of
    # scenario j (whole matrices at once, for speed with M = 10,000)
    l_t <- matrix(lambda_target[, input], n_target, n_sampling)
    l_s <- matrix(lambda_sampling[, input], n_target, n_sampling, byrow = TRUE)
    log_second_moment <- log_second_moment + (l_t^2 / l_s + (-2 * l_t + l_s))
  }

  return(exp(-log_second_moment))
}

# ===== The optimal design: LP (11)
solve_lp_lpsolve <- function(psi) {
  # Solve LP (11) with lpSolve (used in Sections 7.1 and 7.3)
  # Input:
  # -- psi: M-by-M matrix Psi
  # Output:
  # -- vector c = (c_1, ..., c_M)

  m <- nrow(psi) # number of scenarios
  solution <- lpSolve::lp(
    direction = "min",
    objective.in = rep(1, m), # minimize sum_j c_j
    const.mat = psi, # constraint i: sum_j psi[i, j] c_j >= 1
    const.dir = rep(">=", m),
    const.rhs = rep(1, m)
  )
  if (solution$status != 0) {
    stop("lpSolve failed with status ", solution$status)
  }

  return(solution$solution)
}

solve_lp_lpsolveapi <- function(psi, tuned = TRUE) {
  # Solve LP (11) with lpSolveAPI and the bounds 0 <= c_j <= 1 (Section 7.2)
  # Inputs:
  # -- psi: M-by-M matrix Psi
  # -- tuned: TRUE = solver settings of the original code (dual simplex,
  #    devex pricing, Curtis-Reid scaling); FALSE = lpSolveAPI defaults
  # Output: a list of 2 items
  # -- status: solver status, 0 = success
  # -- c: vector c = (c_1, ..., c_M), NULL if the solver failed
  # Note: the original code created the model with make.lp(m, m), i.e., with
  # m empty rows before adding the m constraints. This is kept on purpose:
  # the LP can have several optimal solutions, and the same set-up makes the
  # solver find the same one as in the paper's experiments.

  m <- nrow(psi) # number of scenarios
  model <- lpSolveAPI::make.lp(m, m)
  lpSolveAPI::set.objfn(model, rep(1, m)) # minimize sum_j c_j
  for (i in 1:m) {
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

  return(list(status = status, c = c_values))
}

solve_lp_bounded <- function(psi) {
  # Section 7.2 version: lpSolveAPI with the original settings and two
  # fall-backs for the same LP
  # Input:
  # -- psi: M-by-M matrix Psi
  # Output:
  # -- vector c = (c_1, ..., c_M)
  # The tuned settings occasionally fail numerically (status 5) when two
  # outer scenarios almost coincide (Psi has two nearly identical columns);
  # the original code then stopped the whole batch of macro runs. Fall-backs:
  # lpSolveAPI with its default settings, then lpSolve with geometric scaling
  # and c_j <= 1 as extra rows.

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
  m <- nrow(psi) # number of scenarios
  solution <- lpSolve::lp(
    direction = "min",
    objective.in = rep(1, m),
    const.mat = rbind(psi, diag(m)), # Psi c >= 1 and c <= 1
    const.dir = c(rep(">=", m), rep("<=", m)),
    const.rhs = rep(1, 2 * m),
    scale = 4 # geometric scaling
  )
  if (solution$status != 0) {
    stop("all LP solvers failed")
  }

  return(solution$solution)
}

optimal_design <- function(psi, target_ess, section = c("7.1/7.3", "7.2")) {
  # The optimal design: solve LP (11), choose the sampled scenarios, their
  # inner sample sizes and the pooling weights
  # Inputs:
  # -- psi: M-by-M matrix Psi
  # -- target_ess: N, the target effective sample size
  # -- section: "7.1/7.3" (plain LP with lpSolve; scenarios with c_j > 1e-8
  #    are sampled) or "7.2" (set-up of the original Section 7.2 code: Psi
  #    entries below 1e-6 set to 0, bounds 0 <= c_j <= 1, every scenario with
  #    c_j > 0 sampled, at least 10 inner replications each)
  # Output: a list of 4 items
  # -- sampled: indices j of the sampled scenarios (length J)
  # -- n_inner: their inner sample sizes n_j (length J)
  # -- budget: total number of inner replications, sum(n_inner)
  # -- weights: M-by-J matrix of pooling weights gamma_ij (rows sum to 1)

  section <- match.arg(section)
  if (section == "7.2") {
    psi[psi < 1e-6] <- 0 # for numerical stability
    c_values <- solve_lp_bounded(psi)
    sampled <- which(c_values > 0)
    min_reps <- 10
  } else {
    c_values <- solve_lp_lpsolve(psi)
    sampled <- which(c_values > 1e-8)
    min_reps <- 0
  }
  n_inner <- pmax(ceiling(c_values[sampled] * target_ess), min_reps)

  #-- pooling weights gamma_ij proportional to n_j Psi[i, j], each row
  #-- normalized to sum to 1 (Eq. (10), with the rounded-up n_j)
  weights <- psi[, sampled, drop = FALSE]
  for (b in seq_along(sampled)) {
    weights[, b] <- weights[, b] * n_inner[b]
  }
  weights <- weights / rowSums(weights)

  return(list(sampled = sampled, n_inner = n_inner, budget = sum(n_inner),
              weights = weights))
}

# ===== The pooled self-normalized LR estimator
self_normalized_estimate <- function(y, log_lr) {
  # Self-normalized LR estimates of mu(theta_i) for every target i, from the
  # n inner replications of ONE sampled scenario j
  # Inputs:
  # -- y: the n outputs
  # -- log_lr: M-by-n matrix of log W_ijk
  # Output:
  # -- vector of length M, mu~_ij for every target i
  # formula: mu~_ij = sum_k W_ijk y_k / sum_k W_ijk
  # Each row of log_lr is shifted by its maximum before exponentiating: the
  # estimate does not change (the factor cancels in the ratio), and this
  # prevents overflow and underflow.

  log_lr <- log_lr - matrixStats::rowMaxs(log_lr) # shift each row
  lr <- exp(log_lr)

  return(as.vector(lr %*% y) / rowSums(lr))
}

pooled_estimate <- function(inner_samples, weights) {
  # Pooled estimate mu^(theta_i) = sum_j gamma_ij mu~_ij (Eq. (8))
  # Inputs:
  # -- inner_samples: list, one element per sampled scenario (in the order of
  #    design$sampled), each list(y = outputs, log_lr = M-by-n matrix)
  # -- weights: design$weights
  # Output:
  # -- vector of length M, the pooled estimates

  estimates <- matrix(0, nrow(weights), ncol(weights)) # [i, b] = mu~_ij
  for (b in seq_along(inner_samples)) {
    estimates[, b] <- self_normalized_estimate(inner_samples[[b]]$y,
                                               inner_samples[[b]]$log_lr)
  }

  return(rowSums(estimates * weights))
}

# ===== Log likelihood ratios
log_lr_lognormal <- function(log_x, target_log_means, sampling_log_mean,
                             log_sd) {
  # log W_ij(x) for lognormal inputs with a common log-sd: the log density
  # with log-mean a_i minus that with log-mean a_j
  # Inputs:
  # -- log_x: the n log inputs simulated at the sampling scenario j
  # -- target_log_means: (a_1, ..., a_M)
  # -- sampling_log_mean: a_j
  # -- log_sd: common log-sd
  # Output:
  # -- M-by-n matrix of log W_ijk
  # formula: (x - a_j)^2 - (x - a_i)^2 = (a_i - a_j)(2x - a_i - a_j)

  d <- target_log_means - sampling_log_mean # a_i - a_j, length M
  log_lr <- (outer(d, 2 * log_x) - d * (target_log_means + sampling_log_mean)) /
    (2 * log_sd^2)

  return(log_lr)
}

log_lr_poisson <- function(x, lambda_target, lambda_sampling) {
  # log W_ij(x) for independent Poisson inputs
  # Inputs:
  # -- x: n-by-d matrix of inputs simulated at the sampling scenario j
  # -- lambda_target: M-by-d matrix of target means
  # -- lambda_sampling: length-d vector of the means of scenario j
  # Output:
  # -- M-by-n matrix of log W_ijk
  # formula: log W_ijk = sum_l [x_kl log(lambda_il / lambda_jl)
  #                             - (lambda_il - lambda_jl)]

  n_target <- nrow(lambda_target) # number of target scenarios (M)
  log_sampling <- matrix(log(lambda_sampling), n_target, length(lambda_sampling),
                         byrow = TRUE) # same row for every target
  log_ratio <- log(lambda_target) - log_sampling # M-by-d
  shift <- rowSums(lambda_target) - sum(lambda_sampling) # length M

  return(log_ratio %*% t(x) - shift)
}
