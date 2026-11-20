################################################################################
# assemble_rerun.R -- collect the November 2026 reruns into results/rerun_full
#
# Table 1 at 10,000 macro runs was run in two parts (M = 128-1024 in
# results/rerun_10000_smallM, M = 2048-4096 in results/rerun_10000_largeM).
# This script stacks them into results/rerun_full, which already holds the
# Figure 2, Section 7.2 and Section 7.3 reruns, and draws the rerun figures
# in figures/rerun_full/.
#
#   Rscript tools/assemble_rerun.R
################################################################################

source(file.path("R", "utilities.R"))
source_library("R")
out_dir <- file.path("results", "rerun_full")
parts <- file.path("results", c("rerun_10000_smallM", "rerun_10000_largeM"))

per_macro <- list()
budgets <- list()
for (i in seq_along(parts)) {
  per_macro[[i]] <- read_summary(parts[i], "sec71_table1_per_macro.csv.gz")
  budgets[[i]] <- read_summary(parts[i], "sec71_table1_budgets.csv")
}
per_macro <- stack_rows(per_macro)
budgets <- stack_rows(budgets)
stopifnot(!anyDuplicated(budgets$M)) # each M in exactly one part

write.csv(per_macro, gzfile(file.path(out_dir, "sec71_table1_per_macro.csv.gz")),
          row.names = FALSE)
write.csv(budgets[order(budgets$M), ],
          file.path(out_dir, "sec71_table1_budgets.csv"), row.names = FALSE)
file.copy(file.path(parts[1], "sec71_table1_truth.csv"), out_dir,
          overwrite = TRUE)
n_per_m <- table(per_macro$M[per_macro$design == "optimal"])
cat("Table 1: M =", names(n_per_m), "; macro runs per M:", unique(n_per_m), "\n")

figure_dir <- file.path("figures", "rerun_full")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure1(figure_dir)
figure2(out_dir, figure_dir)
figure3(out_dir, figure_dir)
cat("figures written to", figure_dir, "\n")
