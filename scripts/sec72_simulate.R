################################################################################
# sec72_simulate.R                                                             #
# Simulations of Section 7.2 (multi-asset ERM) of the paper: Tables 2 and 3    #
# Writes summary files (same format as results/published_run/) to out_dir:     #
# -- sec72_per_macro.csv: per macro run and design, the mean squared error     #
#    of mu^(theta_i) over the scenarios (for the AMSE) and the nested          #
#    statistics                                                                #
# -- sec72_budgets_per_macro.csv: optimal budgets, sampled scenarios, SNS      #
#    sizes                                                                     #
# -- sec72_cpu_times_per_macro.csv: wall-clock time per design                 #
# -- sec72_truth.csv: true nested statistics, from the true conditional        #
#    means of all macro runs together                                          #
# Usage (from the Code/ folder):                                               #
#   Rscript scripts/sec72_simulate.R [n_macro] [n_cores] [out_dir] [mode]      #
# -- n_macro: number of macro runs (default 1000, as in the paper)             #
# -- n_cores: number of parallel workers (default 1)                           #
# -- out_dir: output folder (default results/rerun)                            #
# -- mode: full = all four designs (default); regression = only the optimal    #
#    budgets and the regression design, fitted with Laguerre degrees 2 AND 3   #
#    to the same simulated data (a sensitivity check of the regression         #
#    basis); output files end in _degree_check                                 #
# Run time (full mode) is dominated by SNS+ (10^6 paths of 100 steps per       #
# macro run): about 200 s per macro run per core (Xeon Gold 6126); 1,000       #
# macro runs took 2.4 hours on 24 cores (README, Run times).                   #
################################################################################

args <- commandArgs(trailingOnly = TRUE)
n_macro <- if (length(args) >= 1) as.integer(args[1]) else 1000
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1
out_dir <- if (length(args) >= 3) args[3] else file.path("results", "rerun")
mode <- if (length(args) >= 4) args[4] else "full"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#-- load the library and set up the parallel workers
source(file.path("R", "utilities.R"))
source_library("R")
setup_parallel(n_cores)
params <- multi_asset_params()

#-- experiment parameters
N_OUTER <- 1000 # outer scenarios per macro run (M)
TARGET_ESS <- 1000 # N: target ESS of the optimal design, inner reps of SNS+
if (mode == "full") {
  SEED <- 72001
  file_suffix <- ""
  regression_degrees <- 2
} else {
  SEED <- 72002
  file_suffix <- "_degree_check"
  regression_degrees <- c(2, 3)
}

design_statistics <- function(estimate, true_values) {
  # The statistics stored for one design in one macro run
  # Inputs:
  # -- estimate: estimates of mu(theta_i)
  # -- true_values: true mu(theta_i)
  # Output:
  # -- named vector: mean squared error of the estimates (for the AMSE) and
  #    the nested statistics

  return(c(sq_err_mean = mean((estimate - true_values)^2),
           multi_asset_stats(estimate)))
}

one_macro_run <- function(k) {
  # One macro replication: new outer scenarios and all designs
  # Input:
  # -- k: index of the macro replication
  # Output: a list of 4 items
  # -- statistics: list, one vector of statistics per design
  # -- truth: true mu(theta_i) of the scenarios
  # -- seconds: time per design
  # -- budgets: optimal budget, number of sampled scenarios, SNS sizes

  theta <- multi_asset_scenarios(N_OUTER, params) # outer scenarios
  truth <- multi_asset_true_value(theta, params) # true conditional means
  statistics <- list() # placeholder, one element per design
  seconds <- c() # placeholder, time per design

  #-- SNS+ and the optimal design

  if (mode == "full") {
    start <- Sys.time()
    estimate <- multi_asset_standard(theta, TARGET_ESS, params)
    seconds["SNS+"] <- seconds_since(start)
    statistics[["SNS+"]] <- design_statistics(estimate, truth)

    optimal <- multi_asset_optimal(theta, TARGET_ESS, params)
    statistics[["optimal"]] <- design_statistics(optimal$estimate, truth)
    seconds["optimal"] <- optimal$seconds_lp + optimal$seconds_simulation
    seconds["optimal_lp"] <- optimal$seconds_lp
    budget <- optimal$budget
    n_sampled <- optimal$n_sampled
  } else {
    # budgets only: solve the LP for each stock, no inner simulation
    budget <- numeric(params$n_stocks) # placeholder
    n_sampled <- numeric(params$n_stocks) # placeholder
    for (stock in 1:params$n_stocks) {
      design <- multi_asset_design(theta[, stock], stock, TARGET_ESS, params)
      budget[stock] <- design$budget
      n_sampled[stock] <- length(design$sampled)
    }
  }
  average_budget <- mean(budget) # average per-stock optimal budget

  #-- regression: one replication at ceiling(average_budget) design points
  regression <- multi_asset_regression(theta, ceiling(average_budget), params,
                                       regression_degrees)
  for (degree_name in names(regression$estimates)) {
    if (mode == "full") {
      design_name <- "regression"
    } else {
      design_name <- paste0("regression_", degree_name) # e.g. regression_deg3
    }
    statistics[[design_name]] <- design_statistics(
      regression$estimates[[degree_name]], truth)
  }
  seconds["regression"] <- regression$seconds_simulation +
    regression$seconds_fit[1]

  #-- SNS with the optimal design's budget
  sns_outer <- NA
  sns_inner <- NA
  if (mode == "full") {
    sns <- multi_asset_sns(average_budget, params)
    statistics[["SNS"]] <- design_statistics(sns$estimate, sns$truth)
    seconds["SNS"] <- sns$seconds
    sns_outer <- sns$n_outer
    sns_inner <- sns$n_inner
  }

  return(list(statistics = statistics, truth = truth, seconds = seconds,
              budgets = c(budget_optimal = average_budget,
                          n_sampled_optimal = mean(n_sampled),
                          outer_SNS = sns_outer, inner_SNS = sns_inner)))
}

#-- run the macro replications
start <- Sys.time()
results <- run_macros(n_macro, one_macro_run, seed = SEED)
message(sprintf("Section 7.2 (%s): %d macro runs in %.0f s", mode, n_macro,
                seconds_since(start)))

#-- collect the results of all macro runs into tables
statistics_rows <- list() # placeholder, one row per macro run and design
budget_rows <- list() # placeholder, one row per macro run
time_rows <- list() # placeholder, one block per macro run
true_values <- list() # placeholder, true values of every macro run
for (k in 1:n_macro) {
  result <- results[[k]]
  for (design_name in names(result$statistics)) {
    statistics_rows[[length(statistics_rows) + 1]] <- data.frame(
      macro = k, design = design_name, t(result$statistics[[design_name]]),
      check.names = FALSE)
  }
  budget_rows[[k]] <- data.frame(macro = k, t(result$budgets))
  time_rows[[k]] <- data.frame(macro = k, design = names(result$seconds),
                               seconds = unname(result$seconds))
  true_values[[k]] <- result$truth
}
budgets <- stack_rows(budget_rows)
budgets$budget_SNS <- budgets$outer_SNS * budgets$inner_SNS
truth <- multi_asset_stats(unlist(true_values)) # all macro runs together

#-- write the summary files
output_file <- function(name) {
  # Path of an output file (with the suffix of the mode)

  return(file.path(out_dir, paste0(name, file_suffix, ".csv")))
}
write.csv(stack_rows(statistics_rows), output_file("sec72_per_macro"),
          row.names = FALSE)
write.csv(budgets, output_file("sec72_budgets_per_macro"), row.names = FALSE)
write.csv(stack_rows(time_rows), output_file("sec72_cpu_times_per_macro"),
          row.names = FALSE)
write.csv(data.frame(measure = names(truth), value = unname(truth)),
          output_file("sec72_truth"), row.names = FALSE)
