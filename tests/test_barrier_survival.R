################################################################################
# test_barrier_survival.R
#
# Checks the fast barrier-survival method of barrier_survive() (one random
# number per path) against
#   (1) the method of the original 2024 code (the Brownian-bridge minimum
#       drawn on every time step, bridge_min() below), and
#   (2) the closed-form price of the down-and-out put, barrier_dop_price().
#
#   Rscript tests/test_barrier_survival.R      (about a minute)
################################################################################

source(file.path("R", "utilities.R"))
source_library("R")

# The original code's method: on every step between grid points b1 and b2,
# draw the minimum of the Brownian bridge,
#   min = (b1 + b2 - sqrt((b1 - b2)^2 - 2 vol^2 dt log(U))) / 2,
# and return the path minimum (one uniform U per step and path).
bridge_min <- function(log_paths, vol, dt) {
  step_start <- log_paths[-nrow(log_paths), , drop = FALSE]
  step_end <- log_paths[-1, , drop = FALSE]
  u <- matrix(runif(length(step_end)), nrow = nrow(step_end))
  step_min <- ((step_end + step_start) -
                 sqrt((step_end - step_start)^2 - 2 * vol^2 * dt * log(u))) / 2
  exp(matrixStats::colMins(step_min))
}

set.seed(2026)
s0 <- 100
strike <- 110
barrier <- 90
vol <- 0.35
rf <- 0.02
dt <- 1 / 52
n_steps <- 100
n_paths <- 2e5

# Paths start one step after time 0, as in the inner simulation of Section 7.2
log_returns <- matrix(rnorm(n_steps * n_paths, (rf - vol^2 / 2) * dt,
                            vol * sqrt(dt)), nrow = n_steps)
log_paths <- log(s0) + matrixStats::colCumsums(log_returns)
terminal <- exp(log_paths[n_steps, ])
put_payoff <- pmax(strike - terminal, 0)
discount <- exp(-rf * n_steps * dt)

survived_new <- barrier_survive(log_paths, barrier, vol, dt,
                                needed = matrix(put_payoff > 0, ncol = 1))
survived_new[is.na(survived_new)] <- FALSE
survived_old <- bridge_min(log_paths, vol, dt) > barrier

closed_form <- barrier_dop_price(s0, strike, barrier, rf, 0, vol, dt,
                                 n_steps * dt)
for (method in c("new", "old")) {
  survived <- if (method == "new") survived_new[, 1] else survived_old
  payoff <- discount * put_payoff * survived
  estimate <- mean(payoff)
  se <- standard_error(payoff)
  cat(sprintf("%s method: %.4f (SE %.4f), closed form %.4f, z = %.2f\n",
              method, estimate, se, closed_form, (estimate - closed_form) / se))
}
