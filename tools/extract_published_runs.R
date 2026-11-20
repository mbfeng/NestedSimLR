################################################################################
# extract_published_runs.R
#
# One-off provenance script (run once, November 2026). It reads the large
# .RData workspaces saved by the ORIGINAL 2024 code (kept offline, see
# README) and writes the compact per-macro summaries in
# results/published_run/. Those summaries are all that is needed to redraw
# every table and figure of the paper, so the multi-GB workspaces are no
# longer needed by the archive.
#
# Usage: Rscript tools/extract_published_runs.R <dir with original results> [part]
#   <dir> is the folder with the original code's subfolders listed below;
#   part is all (default), fig2, table1, tables23 or sec73.
#
# The original files and what they are:
#   Section_7_1_SingleAssetERM/SingleAssetERM_1_1000Inner1000Outer10000Macro.RData
#       Figure 2 (10,000 macro runs, M = N = 1000)
#   Section_7_1_SingleAssetERM/SingleAssetERM_2_6Cases10000Macros.RData
#       Table 1 (6 values of M, 10,000 macro runs)
#   Section_7_2_MultiAssetERM/MultiAssetERM_10Stocks1000Outer40Macro60xxSeed.RData
#       Tables 2-3, seeds 6001-6025 x 40 macro runs = 1,000 macro runs
#   Section_7_3_Newsvendor/Newsvendor_1000Inner1000Outer1000Macro.RData
#       Table 4
#   Section_7_3_Newsvendor/Newsvendor_VarRatio_1000Inner1000Outer1000Macro.RData
#       Figure 3(a)
#   Section_7_3_Newsvendor/Newsvendor_SimBudgetGrowth.RData
#       Figure 3(b)
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript tools/extract_published_runs.R <dir with original results>")
}
src <- args[1]
which_part <- if (length(args) >= 2) args[2] else "all"
out <- file.path("results", "published_run")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "risk_measures.R"))

load_env <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e
}

#===============================================================================
# Section 7.1, Figure 2: 95% bands of the conditional-mean estimates
#===============================================================================
if (which_part %in% c("all", "fig2")) {
  e <- load_env(file.path(
    src, "Section_7_1_SingleAssetERM",
    "SingleAssetERM_1_1000Inner1000Outer10000Macro.RData"
  ))
  band <- function(val, design) {
    data.frame(
      design = design, theta = e$s_tau, true = e$val_true,
      mean = rowMeans(val),
      lower = apply(val, 1, quantile, probs = 0.025),
      upper = apply(val, 1, quantile, probs = 0.975)
    )
  }
  fig2 <- rbind(
    band(e$val_opt_nested, "optimal"),
    band(e$val_std_plus, "SNS+"),
    band(e$val_reg_nested, "regression")
  )
  write.csv(fig2, file.path(out, "sec71_fig2_bands.csv"), row.names = FALSE)
  writeLines(
    sprintf("optimal design budget = %d, non-zero scenarios = %d",
            e$sim_budget_opt_nested, e$n_nz),
    file.path(out, "sec71_fig2_design.txt")
  )
  rm(e)
  gc()
  message("Figure 2 summary written.")
}

#===============================================================================
# Section 7.1, Table 1: per-macro risk-measure estimates
#===============================================================================
if (which_part %in% c("all", "table1")) {
  e <- load_env(file.path(
    src, "Section_7_1_SingleAssetERM",
    "SingleAssetERM_2_6Cases10000Macros.RData"
  ))
  designs <- list(
    optimal = "val_opt_nested_list", SNS = "val_std_nested_list",
    `SNS+` = "val_std_plus_list", regression = "val_reg_nested_list"
  )
  rows <- list()
  for (case in seq_along(e$n_outer_set)) {
    for (d in names(designs)) {
      r <- risk_measures_by_macro(e[[designs[[d]]]][[case]],
                                  e$alpha, e$threshold)
      rows[[length(rows) + 1]] <- data.frame(
        M = e$n_outer_set[case], design = d, macro = seq_len(nrow(r)), r
      )
    }
  }
  tab1 <- do.call(rbind, rows)
  write.csv(tab1, gzfile(file.path(out, "sec71_table1_per_macro.csv.gz")),
            row.names = FALSE)
  truth <- data.frame(
    measure = c("VaR", "indicator", "hockeystick", "square"),
    value = as.numeric(e$risks_true)
  )
  write.csv(truth, file.path(out, "sec71_table1_truth.csv"), row.names = FALSE)
  write.csv(data.frame(M = e$n_outer_set, budget_optimal = e$sim_budget),
            file.path(out, "sec71_table1_budgets.csv"), row.names = FALSE)
  rm(e)
  gc()
  message("Table 1 summary written.")
}

#===============================================================================
# Section 7.2, Tables 2-3: per-macro summaries, seeds 6001-6025
#===============================================================================
if (which_part %in% c("all", "tables23")) {
  dir72 <- file.path(src, "Section_7_2_MultiAssetERM")
  files <- sort(list.files(dir72,
    pattern = "^MultiAssetERM_10Stocks1000Outer40Macro60[0-9]{2}Seed\\.RData$",
    full.names = TRUE
  ))
  stopifnot(length(files) == 25)
  per_macro <- list()
  true_all <- list()
  for (f in files) {
    e <- load_env(f)
    seed <- as.integer(sub(".*Macro(\\d+)Seed.*", "\\1", f))
    for (m in seq_len(e$num_macro)) {
      true_m <- e$val_true[, m]
      true_all[[length(true_all) + 1]] <- true_m
      ests <- list(
        `SNS+` = e$val_nested_plus[, m], optimal = e$val_opt_nested[, m],
        SNS = e$val_std_nested[[m]], regression = e$val_reg_nested[, m]
      )
      truths <- list(
        `SNS+` = true_m, optimal = true_m,
        SNS = e$val_true_nested[[m]], regression = true_m
      )
      for (d in names(ests)) {
        per_macro[[length(per_macro) + 1]] <- data.frame(
          seed = seed, macro = m, design = d,
          sq_err_mean = mean((ests[[d]] - truths[[d]])^2), # for AMSE
          mean = mean(ests[[d]]),
          VaR = unname(quantile(ests[[d]], 0.995)),
          indicator = mean(ests[[d]] > 2800),
          hockeystick = mean(pmax(ests[[d]] - 2800, 0)),
          square = mean(pmax(ests[[d]] - 2800, 0)^2)
        )
      }
    }
    budgets <- data.frame(
      seed = seed, macro = seq_len(e$num_macro),
      budget_optimal = colMeans(e$sim_budgets),
      n_sampled_optimal = colMeans(e$total_num_samp),
      outer_SNS = e$sim_budget_std_nested[1, ],
      inner_SNS = e$sim_budget_std_nested[2, ],
      budget_SNS = e$sim_budget_std_nested[3, ]
    )
    times <- data.frame(
      seed = seed,
      design = c("SNS+", "optimal", "SNS", "regression"),
      sum_time = c(e$time_nested_plus[1], e$time_opt_nested_tot[1],
                   e$time_std_nested[1], e$time_reg_nested_tot[1]),
      n_macro = e$num_macro
    )
    write.table(budgets, file.path(out, "sec72_budgets_per_macro.csv"),
      sep = ",", row.names = FALSE, append = f != files[1],
      col.names = f == files[1]
    )
    write.table(times, file.path(out, "sec72_cpu_times_per_seed.csv"),
      sep = ",", row.names = FALSE, append = f != files[1],
      col.names = f == files[1]
    )
    rm(e)
    gc()
    message("  ", basename(f))
  }
  write.csv(do.call(rbind, per_macro),
            file.path(out, "sec72_per_macro.csv"), row.names = FALSE)
  v <- unlist(true_all)
  truth <- data.frame(
    measure = c("mean", "VaR", "indicator", "hockeystick", "square"),
    value = c(mean(v), unname(quantile(v, 0.995)), mean(v > 2800),
              mean(pmax(v - 2800, 0)), mean(pmax(v - 2800, 0)^2))
  )
  write.csv(truth, file.path(out, "sec72_truth.csv"), row.names = FALSE)
  message("Tables 2-3 summary written.")
}

#===============================================================================
# Section 7.3, Table 4 and Figure 3
#===============================================================================
if (which_part %in% c("all", "sec73")) {
  dir73 <- file.path(src, "Section_7_3_Newsvendor")
  e <- load_env(file.path(dir73, "Newsvendor_1000Inner1000Outer1000Macro.RData"))
  designs <- c(oracle = "oracle", optimal = "opt_nested", `SNS+` = "std_plus",
               SNS = "std_nested", regression = "reg_nested")
  rows <- list()
  for (d in names(designs)) {
    cov <- e[[paste0("ci_coverage_", designs[[d]])]]
    wid <- e[[paste0("ci_width_", designs[[d]])]]
    for (i in seq_along(e$pctile)) {
      rows[[length(rows) + 1]] <- data.frame(
        design = d, level = e$pctile[i], macro = seq_len(ncol(cov)),
        coverage = cov[i, ], width = wid[i, ]
      )
    }
  }
  write.csv(do.call(rbind, rows),
            gzfile(file.path(out, "sec73_table4_per_macro.csv.gz")),
            row.names = FALSE)
  write.csv(data.frame(macro = seq_along(e$sim_budget_opt_nested),
                       budget_optimal = e$sim_budget_opt_nested),
            file.path(out, "sec73_table4_budgets.csv"), row.names = FALSE)
  rm(e)

  e <- load_env(file.path(dir73,
    "Newsvendor_VarRatio_1000Inner1000Outer1000Macro.RData"))
  write.csv(data.frame(scenario = seq_along(e$var_ratio),
                       var_sns = e$var_oracle, var_optimal = e$var_opt_nested,
                       ratio = e$var_ratio),
            file.path(out, "sec73_fig3a_var_ratio.csv"), row.names = FALSE)
  rm(e)

  e <- load_env(file.path(dir73, "Newsvendor_SimBudgetGrowth.RData"))
  write.csv(data.frame(M = rep(e$n_outer, each = e$n_macro),
                       macro = rep(seq_len(e$n_macro), length(e$n_outer)),
                       budget_optimal = as.vector(e$sim_budget)),
            file.path(out, "sec73_fig3b_budgets.csv"), row.names = FALSE)
  message("Section 7.3 summaries written.")
}
