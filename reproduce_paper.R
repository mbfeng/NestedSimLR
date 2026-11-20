################################################################################
# reproduce_paper.R
#
# Reproduces all tables and figures of Section 7 of
#   Feng, B. M. and Song, E. (2025). Efficient nested simulation experiment
#   design via the likelihood ratio method. INFORMS Journal on Computing
#   37(3):723-742. https://doi.org/10.1287/ijoc.2022.0392
#
# Run from this folder:   Rscript reproduce_paper.R
#
# By default (RERUN_SIMULATIONS = FALSE) the tables and figures are redrawn
# from the summaries of the simulation runs reported in the paper
# (results/published_run/, about 13 MB) in under a minute, and then checked
# against every published number (verify_against_paper.R).
#
# With RERUN_SIMULATIONS = TRUE all simulations are run again from scratch
# (about 4.5 hours on 24 cores; see README, Run times), the summaries are
# written to results/rerun/, and the tables and figures are drawn from those.
# Reruns use different random numbers from the 2024 runs, so they agree with
# the paper within Monte Carlo error, not digit for digit.
################################################################################

RERUN_SIMULATIONS <- FALSE
N_CORES <- max(1, parallel::detectCores() - 2) # used only for reruns

source(file.path("R", "utilities.R"))
source_library("R")

rscript <- file.path(R.home("bin"), "Rscript")

# Run one simulation script with the given command-line arguments.
run_script <- function(script, script_args) {
  command <- paste(shQuote(rscript), script, paste(script_args, collapse = " "))
  message("Running: ", command)
  if (system(command) != 0) {
    stop(script, " failed")
  }
}

if (RERUN_SIMULATIONS) {
  results_dir <- file.path("results", "rerun")
  run_script("scripts/sec71_simulate.R", c("all", 10000, N_CORES, results_dir))
  run_script("scripts/sec72_simulate.R", c(1000, N_CORES, results_dir, "full"))
  run_script("scripts/sec73_simulate.R", c("all", N_CORES, results_dir))
} else {
  results_dir <- file.path("results", "published_run")
}

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

#--- Figures ---------------------------------------------------------------------
figure1("figures")
figure2(results_dir, "figures")
figure3(results_dir, "figures")

#--- Tables ----------------------------------------------------------------------
t1 <- table1(results_dir)
write.csv(t1, "tables/Table1_SingleAssetERM_MSE.csv", row.names = FALSE)
t2 <- table2(results_dir)
write.csv(t2, "tables/Table2_MultiAssetERM_Budget_CPU.csv", row.names = FALSE)
t3 <- table3(results_dir)
write.csv(t3, "tables/Table3_MultiAssetERM_MSE.csv", row.names = FALSE)
t4 <- table4(results_dir)
write.csv(t4, "tables/Table4_Newsvendor_CrI.csv", row.names = FALSE)

cat("\n==== Table 1 (MSE; rows = M, columns = designs) ====\n")
for (measure in c("VaR", "indicator", "hockeystick", "square")) {
  cat("\n--", measure, "\n")
  shown <- t1[t1$measure == measure & t1$M != 256, ] # the paper omits M = 256
  print(format_wide(shown, "mse", "M", "design"), row.names = FALSE)
}
cat("\n==== Table 2 ====\n")
print(t2[, c("design", "budget", "cpu_seconds", "n_outer")], digits = 5,
      row.names = FALSE)
cat("\n==== Table 3 (MSE) ====\n")
print(format_wide(t3, "mse", "measure", "design", digits = 4),
      row.names = FALSE)
cat("\n==== Table 4 ====\n")
print(t4[, c("design", "level", "coverage", "coverage_se", "width",
             "width_se")], digits = 4, row.names = FALSE)
numbers73 <- text73(results_dir)
cat(sprintf(paste0("\nSection 7.3: optimal budget %.1f (SE %.1f); ",
                   "Figure 3(a) ratio mean %.2f, max %.2f\n"),
            numbers73$budget_mean, numbers73$budget_se,
            numbers73$ratio_mean, numbers73$ratio_max))

#--- Check against the paper -------------------------------------------------------
system(paste(shQuote(rscript), "verify_against_paper.R", results_dir))
