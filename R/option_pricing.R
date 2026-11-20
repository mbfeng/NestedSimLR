################################################################################
# option_pricing.R
#
# Option prices in the Black-Scholes model, used as the true conditional
# means mu(theta) of the ERM examples (Sections 7.1 and 7.2), and the
# barrier-survival step of the inner simulation of Section 7.2.
#
# All prices are vectorized over the stock price `s`, so the options are
# priced at all outer scenarios in one call.
################################################################################

#--- Closed-form prices --------------------------------------------------------

# European call (type = "c") or put (type = "p").
#   s: stock price(s); k: strike; r: risk-free rate; vol: volatility;
#   t: time to maturity; div: continuous dividend yield.
bs_price <- function(s, k, r, vol, t, type = c("c", "p"), div = 0) {
  type <- match.arg(type)
  discounted_stock <- s * exp(-div * t)
  discounted_strike <- k * exp(-r * t)
  d1 <- (log(discounted_stock / discounted_strike) + vol^2 * t / 2) /
    (vol * sqrt(t))
  d2 <- d1 - vol * sqrt(t)
  if (type == "c") {
    discounted_stock * pnorm(d1) - discounted_strike * pnorm(d2)
  } else {
    discounted_strike * pnorm(-d2) - discounted_stock * pnorm(-d1)
  }
}

# Straddle = call + put with the same strike (the Section 7.1 portfolio).
bs_price_straddle <- function(s, k, r, vol, t, div = 0) {
  call_price <- bs_price(s, k, r, vol, t, "c", div)
  put_price <- bs_price(s, k, r, vol, t, "p", div)
  call_price + put_price
}

# Fixed-strike geometric Asian option, averaging over n equally spaced
# monitoring dates (the current price is not included in the average).
# Closed form, because the geometric average of lognormal prices is
# lognormal.
geo_asian_price <- function(s, k, r, vol, t, n, type = c("c", "p")) {
  type <- match.arg(type)
  vol_avg <- (vol / n) * sqrt((n + 1) * (2 * n + 1) / 6) # vol of the average
  drift_avg <- (r - vol^2 / 2) * (n + 1) / (2 * n) + vol_avg^2 / 2
  d1 <- (log(s / k) + (drift_avg + vol_avg^2 / 2) * t) / (vol_avg * sqrt(t))
  d2 <- d1 - vol_avg * sqrt(t)
  if (type == "c") {
    s * exp((drift_avg - r) * t) * pnorm(d1) - k * exp(-r * t) * pnorm(d2)
  } else {
    k * exp(-r * t) * pnorm(-d2) - s * exp((drift_avg - r) * t) * pnorm(-d1)
  }
}

# Bivariate standard normal CDF: P(X < x, Y < y) with corr(X, Y) = rho.
#
# Vectorized over x and y (rho is a single number). This is a port of
# A. Genz's routine BVNU (Drezner-Wesolowsky method with Gauss-Legendre
# quadrature, accurate to about 1e-15). Its branch for |rho| >= 0.925 is not
# needed here; for such rho the function calls mvtnorm::pmvnorm() instead.
pbvnorm <- function(x, y, rho) {
  n <- max(length(x), length(y))
  x <- rep_len(x, n)
  y <- rep_len(y, n)

  if (abs(rho) >= 0.925) {
    correlation <- matrix(c(1, rho, rho, 1), 2)
    result <- numeric(n)
    for (i in seq_len(n)) {
      result[i] <- mvtnorm::pmvnorm(upper = c(x[i], y[i]), corr = correlation)
    }
    return(result)
  }

  # P(X < x, Y < y) = P(X > -x, Y > -y); Genz's code computes the latter
  h <- -x
  k <- -y

  # Gauss-Legendre weights and nodes (6, 12 or 20 points, by |rho|)
  if (abs(rho) < 0.3) {
    weight <- c(0.1713244923791705, 0.3607615730481384, 0.4679139345726904)
    node <- c(0.9324695142031522, 0.6612093864662647, 0.2386191860831970)
  } else if (abs(rho) < 0.75) {
    weight <- c(0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                0.2031674267230659, 0.2334925365383547, 0.2491470458134029)
    node <- c(0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
              0.5873179542866171, 0.3678314989981802, 0.1252334085114692)
  } else {
    weight <- c(0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
                0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
                0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
                0.1527533871307259)
    node <- c(0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
              0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
              0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
              0.07652652113349733)
  }
  weight <- c(weight, weight)
  node <- c(1 - node, 1 + node)

  hk <- h * k
  hs <- (h * h + k * k) / 2
  half_asin <- asin(rho) / 2
  sn <- sin(half_asin * node) # quadrature points
  # integrand at every (h, k) pair (rows) and quadrature point (columns):
  # exp((sn hk - hs) / (1 - sn^2))
  exponent <- (outer(hk, sn) - hs) /
    matrix(1 - sn^2, n, length(sn), byrow = TRUE)
  integral <- as.vector(exp(exponent) %*% weight)
  probability <- integral * half_asin / (2 * pi) + pnorm(-h) * pnorm(-k)
  pmin(pmax(probability, 0), 1)
}

# The "d" terms of the Black-Scholes-type formulas below:
#   d_minus = (log(ratio) + (r - q) time) / (vol sqrt(time)) - vol sqrt(time) / 2
#   d_plus  = d_minus + vol sqrt(time)
d_minus <- function(ratio, r, q, vol, time) {
  (log(ratio) + (r - q) * time) / (vol * sqrt(time)) - vol * sqrt(time) / 2
}

# Down-and-out put with strike k and barrier h_bar, maturity t, monitored
# continuously from time t1 on (t1 = one time step in the paper), dividend
# yield q. Same formula as BarrierOptionPrice_DOP() in the original code.
barrier_dop_price <- function(s, k, h_bar, r, q, vol, t1, t) {
  rho <- sqrt(t1 / t) # correlation of the Brownian motion at t1 and at t
  power <- 2 * (r - q) / vol^2
  sqrt_t1 <- vol * sqrt(t1)
  sqrt_t <- vol * sqrt(t)

  # d_minus terms (name: dm_<ratio>_<time>); d_plus = d_minus + vol sqrt(time)
  dm_s_h_t1 <- d_minus(s / h_bar, r, q, vol, t1)
  dm_s_k_t <- d_minus(s / k, r, q, vol, t)
  dm_s_h_t <- d_minus(s / h_bar, r, q, vol, t)
  dm_h_s_t1 <- d_minus(h_bar / s, r, q, vol, t1)
  dm_hh_sk_t <- d_minus(h_bar^2 / (s * k), r, q, vol, t)
  dm_h_s_t <- d_minus(h_bar / s, r, q, vol, t)

  term1 <- k * exp(-r * t) *
    (pbvnorm(dm_s_h_t1, -dm_s_k_t, -rho) -
       pbvnorm(dm_s_h_t1, -dm_s_h_t, -rho))
  term2 <- s * exp(-q * t) *
    (pbvnorm(dm_s_h_t1 + sqrt_t1, -(dm_s_k_t + sqrt_t), -rho) -
       pbvnorm(dm_s_h_t1 + sqrt_t1, -(dm_s_h_t + sqrt_t), -rho))
  term3 <- (h_bar / s)^(power - 1) * k * exp(-r * t) *
    (pbvnorm(-dm_h_s_t1, -dm_hh_sk_t, rho) -
       pbvnorm(-dm_h_s_t1, -dm_h_s_t, rho))
  term4 <- (h_bar / s)^(power + 1) * s * exp(-q * t) *
    (pbvnorm(-(dm_h_s_t1 + sqrt_t1), -(dm_hh_sk_t + sqrt_t), rho) -
       pbvnorm(-(dm_h_s_t1 + sqrt_t1), -(dm_h_s_t + sqrt_t), rho))
  term1 - term2 - term3 + term4
}

#--- Barrier survival in the inner simulation ------------------------------------

# Did each simulated path stay above each down-and-out barrier?
# (Continuous monitoring between the grid points of the simulated paths.)
#
# Arguments:
#   log_paths  n_steps x n_paths matrix of simulated log prices
#   barriers   vector of barrier levels H
#   vol, dt    volatility and time step
#   needed     n_paths x length(barriers) logical matrix: only these
#              path/barrier pairs are evaluated (where the put pays off);
#              the others are returned as NA and must not be used
# Returns an n_paths x length(barriers) logical matrix.
#
# Method. Between two grid points with log prices b1, b2 above log(H), a
# Brownian path stays above log(H) with probability
#   1 - exp(-2 (b1 - log H)(b2 - log H) / (vol^2 dt)),
# independently over the steps. A path survives with the product of these
# probabilities, so one uniform random number per path decides. (The
# original code drew the minimum of the Brownian bridge on every step, one
# random number per step; the result has the same distribution, see
# tests/test_barrier_survival.R.) Further:
#   * a path with a grid point at or below log(H) is knocked out;
#   * a path whose lowest grid point is more than sqrt(23 vol^2 dt) above
#     log(H) is treated as surviving: its knock-out probability per step is
#     below exp(-46), about 1e-20.
# Each barrier uses its own random numbers, as in the original code.
barrier_survive <- function(log_paths, barriers, vol, dt, needed) {
  n_paths <- ncol(log_paths)
  n_steps <- nrow(log_paths)
  lowest <- matrixStats::colMins(log_paths)
  safe_distance <- sqrt(23 * vol^2 * dt)
  step_end <- log_paths[-1, , drop = FALSE]          # b1 on every step
  step_start <- log_paths[-n_steps, , drop = FALSE]  # b2 on every step

  survived <- matrix(NA, n_paths, length(barriers))
  for (barrier in seq_along(barriers)) {
    log_h <- log(barriers[barrier])
    paths <- which(needed[, barrier])
    alive <- lowest[paths] > log_h               # never at or below the barrier
    near <- alive & (lowest[paths] - log_h < safe_distance)
    paths_near <- paths[near]
    if (length(paths_near) > 0) {
      above_end <- step_end[, paths_near, drop = FALSE] - log_h
      above_start <- step_start[, paths_near, drop = FALSE] - log_h
      log_survival <- matrixStats::colSums2(
        log1p(-exp(-2 * above_end * above_start / (vol^2 * dt)))
      )
      alive[near] <- runif(length(paths_near)) < exp(log_survival)
    }
    survived[paths, barrier] <- alive
  }
  survived
}
