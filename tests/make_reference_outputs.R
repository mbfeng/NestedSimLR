################################################################################
# make_reference_outputs.R
#
# Runs every library function on small inputs with fixed random seeds and
# saves the outputs. tests/run_tests.R compares the current code's outputs
# with the saved reference (tests/reference_outputs.rds), so that a change
# to the code that changes any result is noticed.
#
#   Rscript tests/make_reference_outputs.R [output file]
#     (default tests/reference_outputs.rds; only remake the reference when a
#     change of results is intended)
################################################################################

make_reference_outputs <- function() {
  out <- list()
  seed <- function() {
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  }

  #--- lr_design.R
  seed()
  a <- sort(rnorm(40, 4.6, 0.3))
  out$psi_lognormal <- psi_lognormal(a, 0.2)
  seed()
  lambda <- matrix(rgamma(30, 5, 1), 15, 2)
  out$psi_poisson <- psi_poisson(lambda)
  out$design_lpsolve <- optimal_design(psi_lognormal(a, 0.2), 50, "7.1/7.3")
  out$design_api <- optimal_design(psi_lognormal(a, 0.2), 50, "7.2")
  seed()
  x <- rnorm(7, 4.6, 0.2)
  out$log_lr_lognormal <- log_lr_lognormal(x, a, a[3], 0.2)
  seed()
  x_poisson <- matrix(rpois(14, 5), 7, 2)
  out$log_lr_poisson <- log_lr_poisson(x_poisson, lambda, lambda[2, ])

  #--- option_pricing.R
  s <- c(60, 90, 100, 120, 150)
  out$bs <- list(bs_price(s, 100, 0.02, 0.3, 1.5, "c"),
                 bs_price(s, 100, 0.02, 0.3, 1.5, "p"),
                 bs_price_straddle(s, 110, 0.02, 0.3, 1.75))
  out$asian <- list(geo_asian_price(s, 100, 0.02, 0.3, 1.9, 100, "c"),
                    geo_asian_price(s, 100, 0.02, 0.3, 1.9, 100, "p"))
  out$pbvnorm <- list(pbvnorm(c(-1, 0, 1.5), c(0.3, -0.2, 2), 0.1),
                      pbvnorm(c(-1, 0, 1.5), c(0.3, -0.2, 2), -0.5),
                      pbvnorm(c(-1, 0, 1.5), c(0.3, -0.2, 2), 0.9),
                      pbvnorm(c(-1, 0, 1.5), c(0.3, -0.2, 2), 0.95))
  out$barrier <- barrier_dop_price(s, 110, 90, 0.02, 0, 0.35, 1 / 52, 1.9)
  seed()
  steps <- matrix(rnorm(20 * 300, 0, 0.05), 20)
  paths <- log(100) + matrixStats::colCumsums(steps)
  seed()
  out$barrier_survive <- barrier_survive(paths, c(85, 95, 105), 0.35, 1 / 52,
                                         needed = matrix(TRUE, 300, 3))

  #--- risk_measures.R
  seed()
  values <- matrix(rnorm(300, 48, 3), 100, 3)
  out$risk <- risk_measures(values[, 1], 0.99, 49)
  out$risk_by_macro <- risk_measures_by_macro(values, 0.99, 49)
  out$error_summary <- error_summary(
    out$risk_by_macro, c(VaR = 50, indicator = 0.4, hockeystick = 1, square = 3))

  #--- Section 7.1
  params <- single_asset_params()
  theta <- single_asset_scenarios(128, params)
  out$sa_true <- single_asset_true_value(theta, params)
  design <- single_asset_design(theta, 128, params)
  out$sa_design <- design
  seed()
  out$sa_inner <- single_asset_inner(theta[5], 9, params)
  seed()
  out$sa_opt <- single_asset_optimal(theta, design, params)
  seed()
  out$sa_std <- single_asset_standard(theta, 10, params)
  sizes <- single_asset_matched_sizes(design$budget)
  out$sa_sizes <- sizes
  seed()
  regression_points <- single_asset_scenarios(sizes$regression_points, params)
  out$sa_reg <- single_asset_regression(theta, regression_points, params)
  out$laguerre <- laguerre_basis(c(0.5, 1, 1.5), 3)
  out$sa_true_risks <- single_asset_true_risks(0.99, 49, params, n_true = 1e5,
                                               chunk = 3e4)

  #--- Section 7.2
  params <- multi_asset_params()
  out$ma_sigma <- params$sigma
  seed()
  theta <- multi_asset_scenarios(60, params)
  out$ma_theta <- theta
  out$ma_true <- multi_asset_true_value(theta, params)
  seed()
  out$ma_inner <- multi_asset_inner(rep(theta[1, 3], 50), 3, params)
  seed()
  out$ma_std <- multi_asset_standard(theta[1:12, ], 7, params, batch_size = 5)
  seed()
  optimal <- multi_asset_optimal(theta, 60, params)
  out$ma_opt <- optimal[c("estimate", "budget", "n_sampled")]
  seed()
  regression <- multi_asset_regression(theta, 80, params, degrees = c(2, 3))
  out$ma_reg <- regression$estimates
  out$ma_laguerre <- multi_asset_laguerre(theta[1:4, ], params, 2)
  seed()
  sns <- multi_asset_sns(40, params)
  out$ma_sns <- sns[c("estimate", "truth", "n_outer", "n_inner")]
  out$ma_stats <- multi_asset_stats(out$ma_true)

  #--- Section 7.3
  params <- newsvendor_setup()
  out$nv_posterior <- cbind(params$posterior_rate, params$posterior_shape)
  lambda <- newsvendor_draw_lambda(50, params) # continues the set.seed(2020) stream
  out$nv_lambda <- lambda
  out$nv_true <- newsvendor_true_profit(lambda, params)
  out$nv_var <- newsvendor_true_variance(lambda, params)
  seed()
  out$nv_inner <- newsvendor_inner(lambda[1, ], 6, params)
  seed()
  out$nv_inner1 <- newsvendor_inner(lambda[1, ], 1, params)
  seed()
  out$nv_std <- newsvendor_standard(lambda, 5, params)
  design <- newsvendor_design(lambda, 50)
  out$nv_design <- design
  seed()
  out$nv_opt <- newsvendor_optimal(lambda, design, params)
  seed()
  out$nv_reg <- newsvendor_regression(lambda, 60, params)
  out$nv_cri <- newsvendor_cri(out$nv_opt, sort(out$nv_true))

  #--- tables.R, from the published run
  published <- file.path("results", "published_run")
  out$table1 <- table1(published)
  out$table2 <- table2(published)
  out$table3 <- table3(published)
  out$table4 <- table4(published)
  out$text73 <- text73(published)

  #--- utilities.R: the same results on 1 and 2 workers
  one_macro_run <- function(k) {
    rnorm(2)
  }
  setup_parallel(1)
  out$macros1 <- run_macros(4, one_macro_run, seed = 9)
  setup_parallel(2)
  out$macros2 <- run_macros(4, one_macro_run, seed = 9)
  setup_parallel(1)
  out
}

if (sys.nframe() == 0) { # run as a script, not sourced
  args <- commandArgs(trailingOnly = TRUE)
  output_file <- file.path("tests", "reference_outputs.rds")
  if (length(args) >= 1) {
    output_file <- args[1]
  }
  source(file.path("R", "utilities.R"))
  source_library("R")
  saveRDS(make_reference_outputs(), output_file)
  cat("reference outputs written to", output_file, "\n")
}
