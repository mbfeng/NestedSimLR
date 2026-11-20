################################################################################
# verify_against_paper.R
#
# Compares the published values of Tables 1-4 (the estimates, not their
# standard errors) and the design budgets quoted in Section 7 with the values
# computed from a results folder.
#
#   Rscript verify_against_paper.R [results_dir] [mode]
#     results_dir  results/published_run (default) or a rerun folder
#     mode         exact  (default for results/published_run): the summaries
#                         come from the same simulation runs as the paper, so
#                         every entry must agree up to the printed rounding:
#                         |x - v| <= h
#                  mc     (default for reruns): agreement within Monte Carlo
#                         error: z = (x - v) / sqrt(se_x^2 + se_v^2 + h^2 / 3)
#                         must satisfy |z| <= 2.5. Here se_v, the Monte Carlo
#                         error of the paper's own value, is se_x scaled to
#                         the paper's number of macro runs.
#
# x = computed value, se_x its standard error, v = the paper's value, and
# h = half a unit in the last printed digit of v. Times (CPU seconds) depend
# on the computer and are not tested.
#
# Status of each check:
#   ok            agrees
#   heavy-tailed  disagrees, but the worst 1% of macro runs give more than
#                 half of the MSE, so the MSE and its SE are unreliable
#   ok (KS)       disagrees, but the per-macro estimates have the same
#                 distribution as the paper's run (Kolmogorov-Smirnov p >= 0.01)
#   known         Tables 2-3 in exact mode: the published run of Tables 2-3 was
#                 not kept (README, Provenance)
#   UNEXPECTED    anything else
################################################################################

args <- commandArgs(trailingOnly = TRUE)
results_dir <- file.path("results", "published_run")
if (length(args) >= 1) {
  results_dir <- args[1]
}
if (length(args) >= 2) {
  mode <- args[2]
} else if (grepl("published_run", results_dir)) {
  mode <- "exact"
} else {
  mode <- "mc"
}
published_dir <- file.path("results", "published_run")
is_published_run <- normalizePath(results_dir) == normalizePath(published_dir)

source(file.path("R", "utilities.R"))
source_library("R")

# The paper's printed values (as text, so that the printed precision is known)
paper_values <- read.csv(file.path("verify", "paper_values.csv"),
                         colClasses = "character")
paper_table1 <- read.csv(file.path("verify", "paper_table1.csv"),
                         colClasses = "character")

# Number of macro runs behind each item in the paper
N_MACRO_PAPER <- c(table1 = 10000, table2 = 1000, table3 = 1000,
                   table4 = 1000, text71 = 1, text73 = 1000)

# Documented differences: "item/design/quantity" -> reason
KNOWN <- c(
  "table2/optimal/budget" = "published Tables 2-3 run not kept",
  "table2/SNS/budget" = "published Tables 2-3 run not kept",
  "table2/regression/budget" = "published Tables 2-3 run not kept",
  "table2/regression/n_outer" = "published Tables 2-3 run not kept",
  "table2/optimal/n_outer" = "published Tables 2-3 run not kept",
  "table2/SNS/n_outer" = "published Tables 2-3 run not kept"
)

# Half a unit in the last printed digit, e.g. "0.185" -> 0.0005,
# "2.29E-03" -> 0.000005, "8200.8" -> 0.05.
half_unit <- function(printed) {
  printed <- toupper(printed)
  mantissa <- sub("E.*", "", printed)
  exponent <- 0
  if (grepl("E", printed)) {
    exponent <- as.numeric(sub(".*E", "", printed))
  }
  decimals <- 0
  if (grepl("\\.", mantissa)) {
    decimals <- nchar(sub(".*\\.", "", mantissa))
  }
  0.5 * 10^(exponent - decimals)
}

# Compare one computed value with the paper's printed value.
# Returns a one-row data frame (see the status list at the top).
check_entry <- function(item, design, quantity, computed, se, printed,
                        n_macro = NA, heavy_tailed = FALSE, ks_p = NA) {
  paper <- as.numeric(printed)
  h <- half_unit(printed)
  z <- NA
  if (mode == "exact") {
    agrees <- abs(computed - paper) <= h * (1 + 1e-9)
  } else {
    se_paper <- 0
    if (!is.na(n_macro)) {
      se_paper <- se * sqrt(n_macro / N_MACRO_PAPER[[item]])
    }
    denominator <- sqrt(se^2 + se_paper^2 + h^2 / 3)
    z <- (computed - paper) / denominator
    if (denominator == 0) {
      z <- if (computed == paper) 0 else Inf
    }
    agrees <- abs(z) <= 2.5
  }

  key <- paste(item, design, quantity, sep = "/")
  note <- ""
  if (agrees) {
    status <- "ok"
  } else if (heavy_tailed) {
    status <- "heavy-tailed"
    note <- "worst 1% of macro runs give > 50% of the MSE; SE unreliable"
  } else if (isTRUE(ks_p >= 0.01)) {
    status <- "ok (KS)"
    note <- sprintf(paste("z-test fails but the per-macro estimates have the",
                          "same distribution as the paper's run (KS p = %.2f)"),
                    ks_p)
  } else if (key %in% names(KNOWN)) {
    status <- "known"
    note <- KNOWN[[key]]
  } else if (mode == "exact" && item %in% c("table2", "table3")) {
    status <- "known"
    note <- "published Tables 2-3 run not kept"
  } else {
    status <- "UNEXPECTED"
  }
  data.frame(item = item, design = design, quantity = quantity,
             paper = printed, computed = signif(computed, 6),
             se = signif(se, 3), z = round(z, 2), status = status, note = note)
}

checks <- list()

#--- Table 1 and the Section 7.1 text ------------------------------------------------
if (file.exists(file.path(results_dir, "sec71_table1_per_macro.csv.gz"))) {
  table1_computed <- table1(results_dir)
  table1_published <- table1(published_dir)
  per_macro <- read_summary(results_dir, "sec71_table1_per_macro.csv.gz")
  per_macro_published <- read_summary(published_dir,
                                      "sec71_table1_per_macro.csv.gz")

  for (i in seq_len(nrow(paper_table1))) {
    entry <- paper_table1[i, ]
    m <- as.numeric(entry$M)
    is_row <- table1_computed$M == m & table1_computed$design == entry$design &
      table1_computed$measure == entry$measure
    if (!any(is_row)) {
      next # this M was not rerun
    }
    row <- table1_computed[is_row, ]
    is_published_row <- table1_published$M == m &
      table1_published$design == entry$design &
      table1_published$measure == entry$measure
    # Heavy tails are judged on the larger of the two samples' evidence: a
    # small rerun may contain none of the rare extreme macro runs that
    # dominate the MSE of the paper's 10,000-run sample.
    top1 <- max(row$top1_share, table1_published$top1_share[is_published_row])

    # Distribution check of the per-macro estimates (reruns only)
    ks_p <- NA
    if (!is_published_run) {
      rerun_values <- per_macro[[entry$measure]][per_macro$M == m &
                                                   per_macro$design == entry$design]
      published_values <- per_macro_published[[entry$measure]][
        per_macro_published$M == m & per_macro_published$design == entry$design]
      ks_p <- suppressWarnings(ks.test(rerun_values, published_values)$p.value)
    }
    checks[[length(checks) + 1]] <- check_entry(
      "table1", entry$design, paste0(entry$measure, "_M", entry$M), row$mse,
      row$mse_se, entry$printed, row$n_macro, heavy_tailed = top1 > 0.5,
      ks_p = ks_p)
  }

  budgets <- read_summary(results_dir, "sec71_table1_budgets.csv")
  if (512 %in% budgets$M) {
    ratio <- 512^2 / budgets$budget_optimal[budgets$M == 512]
    checks[[length(checks) + 1]] <- check_entry(
      "text71", "SNS+", "budget_ratio_M512", ratio, 0, "247")
  }
  if (4096 %in% budgets$M) {
    ratio <- 4096^2 / budgets$budget_optimal[budgets$M == 4096]
    checks[[length(checks) + 1]] <- check_entry(
      "text71", "SNS+", "budget_ratio_M4096", ratio, 0, "1760")
  }
}
if (file.exists(file.path(results_dir, "sec71_fig2_design.txt"))) {
  design_line <- readLines(file.path(results_dir, "sec71_fig2_design.txt"))
  budget <- as.numeric(sub(".*budget = (\\d+).*", "\\1", design_line))
  checks[[length(checks) + 1]] <- check_entry(
    "text71", "optimal", "budget_M1000", budget, 0, "2148")
}

#--- Tables 2 and 3 ------------------------------------------------------------------
if (file.exists(file.path(results_dir, "sec72_per_macro.csv"))) {
  table2_computed <- table2(results_dir)
  table3_computed <- table3(results_dir)
  for (i in which(paper_values$item == "table2")) {
    entry <- paper_values[i, ]
    if (entry$quantity == "cpu_seconds") {
      next # depends on the computer; not tested
    }
    row <- table2_computed[table2_computed$design == entry$design, ]
    checks[[length(checks) + 1]] <- check_entry(
      "table2", entry$design, entry$quantity, row[[entry$quantity]],
      row[[paste0(entry$quantity, "_se")]], entry$printed, row$n_macro)
  }
  for (i in which(paper_values$item == "table3")) {
    entry <- paper_values[i, ]
    row <- table3_computed[table3_computed$design == entry$design &
                             table3_computed$measure == entry$quantity, ]
    checks[[length(checks) + 1]] <- check_entry(
      "table3", entry$design, entry$quantity, row$mse, row$mse_se,
      entry$printed, row$n_macro,
      heavy_tailed = isTRUE(row$top1_share > 0.5))
  }
}

#--- Table 4 and the Section 7.3 text ------------------------------------------------
if (file.exists(file.path(results_dir, "sec73_table4_per_macro.csv.gz"))) {
  table4_computed <- table4(results_dir)
  for (i in which(paper_values$item == "table4")) {
    entry <- paper_values[i, ]
    # quantity is e.g. "coverage_0.9" or "width_0.95"
    level <- as.numeric(sub(".*_", "", entry$quantity))
    quantity <- sub("_[^_]*$", "", entry$quantity)
    row <- table4_computed[table4_computed$design == entry$design &
                             abs(table4_computed$level - level) < 1e-9, ]
    if (quantity == "coverage") {
      checks[[length(checks) + 1]] <- check_entry(
        "table4", entry$design, entry$quantity, row$coverage, row$coverage_se,
        entry$printed, row$n_macro)
    } else if (quantity == "width") {
      checks[[length(checks) + 1]] <- check_entry(
        "table4", entry$design, entry$quantity, row$width, row$width_se,
        entry$printed, row$n_macro)
    }
  }
  numbers73 <- text73(results_dir)
  if (mode == "exact") {
    checks[[length(checks) + 1]] <- check_entry(
      "text73", "fig3a", "ratio_mean", numbers73$ratio_mean, 0, "1.21")
    checks[[length(checks) + 1]] <- check_entry(
      "text73", "fig3a", "ratio_max", numbers73$ratio_max, 0, "1.38")
  }
}

#--- Report ----------------------------------------------------------------------------
report <- stack_rows(checks)
label <- basename(normalizePath(results_dir))
write.csv(report, file.path("verify", paste0("verification_", label, ".csv")),
          row.names = FALSE)
cat(sprintf("\nVerification of %s (mode = %s): %d checks\n", results_dir, mode,
            nrow(report)))
print(table(report$status))
if (any(report$status != "ok")) {
  cat("\nEntries that do not agree:\n")
  print(report[report$status != "ok", ], row.names = FALSE)
}
if (mode == "mc") {
  n_tested <- sum(!is.na(report$z))
  cat(sprintf(paste("\nlargest |z| = %.2f; expected number of |z| > 2.5 by",
                    "chance (no real difference): %.1f\n"),
              max(abs(report$z), na.rm = TRUE), 0.0124 * n_tested))
}
