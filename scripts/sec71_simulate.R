################################################################################
# sec71_simulate.R                                                             #
# Simulations of Section 7.1 (single-asset ERM) of the paper                   #
# Writes summary files (same format as results/published_run/) to out_dir:     #
# -- sec71_fig2_bands.csv: Figure 2, 2.5%/97.5% bands and means                #
# -- sec71_fig2_design.txt: optimal budget for M = N = 1000                    #
# -- sec71_table1_per_macro.csv.gz: Table 1, risk measures per macro run       #
# -- sec71_table1_truth.csv: true risk measures (10^8 scenarios)               #
# -- sec71_table1_budgets.csv: optimal budget for every M                      #
# Usage (from the Code/ folder):                                               #
#   Rscript scripts/sec71_simulate.R [part] [n_macro] [n_cores] [out_dir] [M]  #
# -- part: fig2 | table1 | truth | all (default all)                           #
# -- n_macro: number of macro runs (default 10000, as in the paper)            #
# -- n_cores: number of parallel workers (default 1)                           #
# -- out_dir: output folder (default results/rerun)                            #
# -- M: table1 only, comma-separated values of M (default 128,256,...,4096)    #
################################################################################

args <- commandArgs(trailingOnly = TRUE)
part <- if (length(args) >= 1) args[1] else "all"
n_macro <- if (length(args) >= 2) as.integer(args[2]) else 10000
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 1
out_dir <- if (length(args) >= 4) args[4] else file.path("results", "rerun")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "utilities.R"))
source_library("R")
setup_parallel(n_cores)
params <- single_asset_params()

ALPHA <- 0.99       # VaR level
THRESHOLD <- 49     # threshold of the indicator, hockey-stick and square losses
SEED_FIGURE2 <- 71001
SEED_TABLE1 <- 71100 # + index of the value of M

#-- Figure 2: M = N = 1000
if (part %in% c("all", "fig2")) {
  n_outer <- 1000
  theta <- single_asset_scenarios(n_outer, params)
  design <- single_asset_design(theta, target_ess = 1000, params)
  sizes <- single_asset_matched_sizes(design$budget)
  regression_points <- single_asset_scenarios(sizes$regression_points, params)

  one_macro_run <- function(k) {
    # One macro run: estimates of mu(theta_i) by three designs
    # Input:
    # -- k: index of the macro replication
    # Output:
    # -- list of three vectors of M estimates (optimal, SNS+, regression)

    return(list(
      optimal = single_asset_optimal(theta, design, params),
      `SNS+` = single_asset_standard(theta, 1000, params),
      regression = single_asset_regression(theta, regression_points, params)
    ))
  }
  start <- Sys.time()
  results <- run_macros(n_macro, one_macro_run, seed = SEED_FIGURE2)
  message(sprintf("Figure 2: %d macro runs in %.0f s", n_macro,
                  seconds_since(start)))

  true_values <- single_asset_true_value(theta, params)
  bands <- list()
  for (design_name in c("optimal", "SNS+", "regression")) {
    # estimates[i, k] = estimate of mu(theta_i) in macro run k
    estimates <- matrix(0, n_outer, n_macro)
    for (k in seq_len(n_macro)) {
      estimates[, k] <- results[[k]][[design_name]]
    }
    bands[[design_name]] <- data.frame(
      design = design_name, theta = theta, true = true_values,
      mean = rowMeans(estimates),
      lower = matrixStats::rowQuantiles(estimates, probs = 0.025),
      upper = matrixStats::rowQuantiles(estimates, probs = 0.975)
    )
  }
  write.csv(stack_rows(bands), file.path(out_dir, "sec71_fig2_bands.csv"),
            row.names = FALSE)
  writeLines(sprintf("optimal design budget = %d, non-zero scenarios = %d",
                     design$budget, length(design$sampled)),
             file.path(out_dir, "sec71_fig2_design.txt"))
}

#-- Table 1: M = N for several values of M
if (part %in% c("all", "table1")) {
  # 128, 256, ..., 4096 as in the original code; the paper's Table 1 shows all
  # of them except M = 256.
  m_values <- 2^(7:12)
  if (length(args) >= 5) {
    m_values <- as.numeric(strsplit(args[5], ",")[[1]])
  }
  designs <- c("optimal", "SNS", "SNS+", "regression")
  rows <- list()
  budgets <- numeric(length(m_values))
  for (case in seq_along(m_values)) {
    n_outer <- m_values[case]
    theta <- single_asset_scenarios(n_outer, params)
    design <- single_asset_design(theta, target_ess = n_outer, params)
    budgets[case] <- design$budget
    sizes <- single_asset_matched_sizes(design$budget)
    theta_sns <- single_asset_scenarios(sizes$sns_outer, params)
    regression_points <- single_asset_scenarios(sizes$regression_points, params)

    one_macro_run <- function(k) {
      # One macro run: the four risk measures of each design
      # Input:
      # -- k: index of the macro replication
      # Output:
      # -- 4-by-4 matrix, rows = designs in the order of `designs`,
      #    columns = risk measures

      optimal <- single_asset_optimal(theta, design, params)
      sns <- single_asset_standard(theta_sns, sizes$sns_inner, params)
      sns_plus <- single_asset_standard(theta, n_outer, params)
      regression <- single_asset_regression(theta, regression_points, params)
      return(rbind(optimal = risk_measures(optimal, ALPHA, THRESHOLD),
                   SNS = risk_measures(sns, ALPHA, THRESHOLD),
                   `SNS+` = risk_measures(sns_plus, ALPHA, THRESHOLD),
                   regression = risk_measures(regression, ALPHA, THRESHOLD)))
    }
    start <- Sys.time()
    results <- run_macros(n_macro, one_macro_run, seed = SEED_TABLE1 + case)
    message(sprintf("Table 1, M = %d (budget %d): %d macro runs in %.0f s",
                    n_outer, design$budget, n_macro, seconds_since(start)))

    for (design_name in designs) {
      measures <- matrix(0, n_macro, 4,
                         dimnames = list(NULL, colnames(results[[1]])))
      for (k in seq_len(n_macro)) {
        measures[k, ] <- results[[k]][design_name, ]
      }
      rows[[length(rows) + 1]] <- data.frame(M = n_outer, design = design_name,
                                             macro = seq_len(n_macro), measures)
    }
  }
  write.csv(stack_rows(rows),
            gzfile(file.path(out_dir, "sec71_table1_per_macro.csv.gz")),
            row.names = FALSE)
  write.csv(data.frame(M = m_values, budget_optimal = budgets),
            file.path(out_dir, "sec71_table1_budgets.csv"), row.names = FALSE)
}

#-- True risk measures from 10^8 scenarios (deterministic)
if (part %in% c("all", "truth")) {
  start <- Sys.time()
  truth <- single_asset_true_risks(ALPHA, THRESHOLD, params, n_true = 1e8)
  message(sprintf("True risk measures: %.0f s", seconds_since(start)))
  write.csv(data.frame(measure = names(truth), value = unname(truth)),
            file.path(out_dir, "sec71_table1_truth.csv"), row.names = FALSE)
}
