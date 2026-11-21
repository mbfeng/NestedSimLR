################################################################################
# sec73_simulate.R                                                             #
# Simulations of Section 7.3 (newsvendor IUQ) of the paper                     #
# Writes summary files (same format as results/published_run/) to out_dir:     #
# -- sec73_table4_per_macro.csv.gz: Table 4, CrI coverage and width per        #
#    macro run                                                                 #
# -- sec73_table4_budgets.csv: optimal budget per macro run                    #
# -- sec73_fig3a_var_ratio.csv: Figure 3(a), V(SNS) / V(optimal)               #
# -- sec73_fig3b_budgets.csv: Figure 3(b), optimal budget vs M                 #
# Usage (from the Code/ folder):                                               #
#   Rscript scripts/sec73_simulate.R [part] [n_cores] [out_dir]                #
# -- part: table4 | fig3a | fig3b | all (default all)                          #
# The posterior of lambda, the 10^6 posterior draws used to evaluate CrI       #
# coverage, and the fixed scenarios of Figure 3(a) are drawn after             #
# set.seed(2020) exactly as in the original code, so they are identical to     #
# the paper's. Macro runs use their own random-number streams (run_macros()).  #
################################################################################

args <- commandArgs(trailingOnly = TRUE)
part <- if (length(args) >= 1) args[1] else "all"
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1
out_dir <- if (length(args) >= 3) args[3] else file.path("results", "rerun")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "utilities.R"))
source_library("R")
setup_parallel(n_cores)

N_OUTER <- 1000     # M
TARGET_ESS <- 1000  # N
N_MACRO <- 1000
CRI_LEVELS <- c(0.90, 0.95, 0.99)

#-- Table 4
if (part %in% c("all", "table4")) {
  params <- newsvendor_setup()                        # set.seed(2020), data
  lambda_true <- newsvendor_draw_lambda(1e6, params)  # next draws, same stream
  true_sorted <- sort(newsvendor_true_profit(lambda_true, params))
  rm(lambda_true)

  one_macro_run <- function(k) {
    # One macro run: new posterior draws, all designs, their CrIs
    # Input:
    # -- k: index of the macro replication
    # Output: a list of 2 items
    # -- rows: data frame of CrIs, one row per design and level
    # -- budget: the optimal design's budget

    lambda <- newsvendor_draw_lambda(N_OUTER, params) # outer scenarios
    design <- newsvendor_design(lambda, TARGET_ESS)
    budget <- design$budget
    estimates <- list()
    estimates$oracle <- newsvendor_true_profit(lambda, params)
    estimates$optimal <- newsvendor_optimal(lambda, design, params)
    estimates$`SNS+` <- newsvendor_standard(lambda, TARGET_ESS, params)
    # SNS with the optimal budget: budget^(2/3) new scenarios x budget^(1/3) reps
    lambda_sns <- newsvendor_draw_lambda(ceiling(budget^(2 / 3)), params)
    estimates$SNS <- newsvendor_standard(lambda_sns, ceiling(budget^(1 / 3)),
                                         params)
    estimates$regression <- newsvendor_regression(lambda, budget, params)

    rows <- list()
    for (design_name in names(estimates)) {
      cri <- newsvendor_cri(estimates[[design_name]], true_sorted, CRI_LEVELS)
      rows[[design_name]] <- data.frame(design = design_name, macro = k, cri)
    }
    return(list(rows = stack_rows(rows), budget = budget))
  }
  start <- Sys.time()
  results <- run_macros(N_MACRO, one_macro_run, seed = 73001)
  message(sprintf("Table 4: %d macro runs in %.0f s", N_MACRO,
                  seconds_since(start)))

  rows <- list()
  budgets <- numeric(N_MACRO)
  for (k in seq_len(N_MACRO)) {
    rows[[k]] <- results[[k]]$rows
    budgets[k] <- results[[k]]$budget
  }
  per_macro <- stack_rows(rows)
  write.csv(per_macro[, c("design", "level", "macro", "coverage", "width")],
            gzfile(file.path(out_dir, "sec73_table4_per_macro.csv.gz")),
            row.names = FALSE)
  write.csv(data.frame(macro = seq_len(N_MACRO), budget_optimal = budgets),
            file.path(out_dir, "sec73_table4_budgets.csv"), row.names = FALSE)
}

#-- Figure 3(a): variance ratio at a fixed set of scenarios
if (part %in% c("all", "fig3a")) {
  params <- newsvendor_setup()
  lambda <- newsvendor_draw_lambda(N_OUTER, params) # the paper's scenarios
  design <- newsvendor_design(lambda, TARGET_ESS)
  one_macro_run <- function(k) {
    # One macro run of the optimal design at the fixed scenarios

    return(newsvendor_optimal(lambda, design, params))
  }
  start <- Sys.time()
  results <- run_macros(N_MACRO, one_macro_run, seed = 73002)
  message(sprintf("Figure 3(a): %d macro runs in %.0f s", N_MACRO,
                  seconds_since(start)))

  estimates <- matrix(0, N_OUTER, N_MACRO) # [i, k] = estimate in macro run k
  for (k in seq_len(N_MACRO)) {
    estimates[, k] <- results[[k]]
  }
  # SNS with N inner replications: its variance is known exactly
  var_sns <- newsvendor_true_variance(lambda, params) / TARGET_ESS
  var_optimal <- matrixStats::rowVars(estimates)
  write.csv(data.frame(scenario = seq_len(N_OUTER), var_sns = var_sns,
                       var_optimal = var_optimal,
                       ratio = var_sns / var_optimal),
            file.path(out_dir, "sec73_fig3a_var_ratio.csv"), row.names = FALSE)
}

#-- Figure 3(b): growth of the optimal budget with M = N
if (part %in% c("all", "fig3b")) {
  params <- newsvendor_setup()
  m_values <- (1:10) * 1000
  n_repetitions <- 100
  # The LP takes long for large M (M = 10,000: about 13 minutes on 24 cores),
  # so the results are saved after each M, and values of M that are already
  # in the output file are skipped: an interrupted run continues where it
  # stopped.
  output_file <- file.path(out_dir, "sec73_fig3b_budgets.csv")
  done <- NULL
  if (file.exists(output_file)) {
    done <- read.csv(output_file)
  }
  # Each worker holds dense M x M matrices (Psi and lpSolve's copy, about
  # 32 M^2 bytes). The number of workers is capped so that they fit in the
  # memory budget: environment variable MEM_BUDGET_GB (default 24).
  memory_budget <- as.numeric(Sys.getenv("MEM_BUDGET_GB", "24")) * 1e9

  for (case in seq_along(m_values)) {
    n_outer <- m_values[case]
    if (!is.null(done) && n_outer %in% done$M) {
      next
    }
    workers <- max(1, min(n_cores, floor(memory_budget / (32 * n_outer^2))))
    setup_parallel(workers)
    one_macro_run <- function(k) {
      # The optimal budget for one set of n_outer posterior draws

      lambda <- newsvendor_draw_lambda(n_outer, params)
      return(newsvendor_design(lambda, target_ess = n_outer)$budget)
    }
    start <- Sys.time()
    results <- run_macros(n_repetitions, one_macro_run, seed = 73100 + case)
    message(sprintf("Figure 3(b), M = %d: %.0f s (%d workers)", n_outer,
                    seconds_since(start), workers))
    new_rows <- data.frame(M = n_outer, macro = seq_len(n_repetitions),
                           budget_optimal = unlist(results))
    done <- rbind(done, new_rows)
    write.csv(done, output_file, row.names = FALSE)
  }
  setup_parallel(n_cores)
}
