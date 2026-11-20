################################################################################
# tables.R
#
# Tables 1-4 of the paper (and the numbers quoted in the text of Section 7),
# computed from the summary files in a results folder (results/published_run/
# or a rerun folder). Every entry comes with its Monte Carlo standard error,
# which verify_against_paper.R uses.
################################################################################

# Read one summary file (CSV, possibly gzipped) from a results folder.
read_summary <- function(results_dir, file_name) {
  read.csv(file.path(results_dir, file_name), check.names = FALSE,
           stringsAsFactors = FALSE)
}

# The truth file (columns measure, value) as a named vector.
read_truth <- function(results_dir, file_name) {
  truth_table <- read_summary(results_dir, file_name)
  truth <- truth_table$value
  names(truth) <- truth_table$measure
  truth
}

RISK_MEASURES <- c("VaR", "indicator", "hockeystick", "square")

#--- Table 1 -------------------------------------------------------------------

# Table 1: MSE (and its standard error) of the four risk measures, for every
# number of scenarios M and every design. One row per (M, design, measure).
table1 <- function(results_dir) {
  per_macro <- read_summary(results_dir, "sec71_table1_per_macro.csv.gz")
  truth <- read_truth(results_dir, "sec71_table1_truth.csv")
  rows <- list()
  for (m in sort(unique(per_macro$M))) {
    for (design in unique(per_macro$design)) {
      runs <- per_macro[per_macro$M == m & per_macro$design == design, ]
      summary <- error_summary(runs[, RISK_MEASURES], truth)
      rows[[length(rows) + 1]] <- data.frame(M = m, design = design,
                                             n_macro = nrow(runs), summary)
    }
  }
  result <- stack_rows(rows)
  result[order(result$measure, result$design, result$M), ]
}

#--- Tables 2 and 3 --------------------------------------------------------------

# Table 2: simulation budget, number of outer scenarios at which inner
# replications are run, and time per macro run, for each design.
table2 <- function(results_dir) {
  budgets <- read_summary(results_dir, "sec72_budgets_per_macro.csv")
  # regression: one inner replication at ceiling(budget) design points
  regression_budget <- ceiling(budgets$budget_optimal)
  result <- data.frame(
    design = c("SNS+", "optimal", "SNS", "regression"),
    budget = c(1e6, mean(budgets$budget_optimal), mean(budgets$budget_SNS),
               mean(regression_budget)),
    budget_se = c(0, standard_error(budgets$budget_optimal),
                  standard_error(budgets$budget_SNS),
                  standard_error(regression_budget)),
    n_outer = c(1000, mean(budgets$n_sampled_optimal),
                mean(budgets$outer_SNS), mean(regression_budget)),
    n_outer_se = c(0, standard_error(budgets$n_sampled_optimal),
                   standard_error(budgets$outer_SNS),
                   standard_error(regression_budget))
  )

  # Time per macro run (wall-clock seconds). Reruns store
  # the time of every macro run; the published run stored totals per seed.
  result$cpu_seconds <- NA
  if (file.exists(file.path(results_dir, "sec72_cpu_times_per_macro.csv"))) {
    times <- read_summary(results_dir, "sec72_cpu_times_per_macro.csv")
    for (i in seq_len(nrow(result))) {
      design_times <- times$seconds[times$design == result$design[i]]
      result$cpu_seconds[i] <- mean(design_times)
    }
  } else {
    times <- read_summary(results_dir, "sec72_cpu_times_per_seed.csv")
    for (i in seq_len(nrow(result))) {
      is_design <- times$design == result$design[i]
      result$cpu_seconds[i] <- sum(times$sum_time[is_design]) /
        sum(times$n_macro[is_design])
    }
  }
  result$n_macro <- nrow(budgets)
  result
}

# Table 3: AMSE of the conditional means (average over the scenarios of the
# squared error of mu^(theta_i)) and the MSEs of the nested statistics.
table3 <- function(results_dir) {
  per_macro <- read_summary(results_dir, "sec72_per_macro.csv")
  truth <- read_truth(results_dir, "sec72_truth.csv")
  rows <- list()
  for (design in sort(unique(per_macro$design))) {
    runs <- per_macro[per_macro$design == design, ]
    amse <- data.frame(measure = "AMSE", bias = NA, var = NA,
                       mse = mean(runs$sq_err_mean),
                       mse_se = standard_error(runs$sq_err_mean),
                       top1_share = NA)
    summary <- rbind(amse, error_summary(runs[, RISK_MEASURES], truth))
    rows[[length(rows) + 1]] <- data.frame(design = design,
                                           n_macro = nrow(runs), summary)
  }
  stack_rows(rows)
}

#--- Table 4 and the numbers of Section 7.3 -------------------------------------

# Table 4: average CrI coverage and width over the macro runs, with their
# standard errors sd / sqrt(K) (K = number of macro runs).
table4 <- function(results_dir) {
  per_macro <- read_summary(results_dir, "sec73_table4_per_macro.csv.gz")
  rows <- list()
  for (design in sort(unique(per_macro$design))) {
    for (level in sort(unique(per_macro$level))) {
      runs <- per_macro[per_macro$design == design & per_macro$level == level, ]
      rows[[length(rows) + 1]] <- data.frame(
        design = design, level = level, n_macro = nrow(runs),
        coverage = mean(runs$coverage),
        coverage_se = standard_error(runs$coverage),
        width = mean(runs$width),
        width_se = standard_error(runs$width)
      )
    }
  }
  stack_rows(rows)
}

# Numbers quoted in the text of Section 7.3: the optimal design's average
# budget and its ratio to SNS+'s 10^6, and the variance ratios of Figure 3(a).
text73 <- function(results_dir) {
  budget <- read_summary(results_dir, "sec73_table4_budgets.csv")$budget_optimal
  ratio <- read_summary(results_dir, "sec73_fig3a_var_ratio.csv")$ratio
  list(budget_mean = mean(budget), budget_se = standard_error(budget),
       budget_ratio = 1e6 / mean(budget),
       ratio_mean = mean(ratio), ratio_max = max(ratio), ratio_min = min(ratio))
}

#--- Printing ---------------------------------------------------------------------

# Rearrange a long table into the paper's layout: one row per value of
# `row_column`, one column per value of `col_column`, entries from
# `value_column` (rounded to `digits` significant digits).
format_wide <- function(table, value_column, row_column, col_column,
                        digits = 3) {
  row_values <- unique(table[[row_column]])
  col_values <- unique(table[[col_column]])
  wide <- data.frame(row_values)
  names(wide) <- row_column
  for (col_value in col_values) {
    column <- numeric(length(row_values))
    for (r in seq_along(row_values)) {
      is_entry <- table[[row_column]] == row_values[r] &
        table[[col_column]] == col_value
      column[r] <- signif(table[[value_column]][is_entry], digits)
    }
    wide[[col_value]] <- column
  }
  wide
}
