################################################################################
# utilities.R
#
# Small helpers used by all scripts:
#   * loading the library (all files in R/),
#   * running macro replications reproducibly, serially or in parallel,
#   * timing, and stacking a list of data frames into one.
#
# Reproducibility: macro run k always uses the k-th random-number stream
# derived from the experiment's seed (L'Ecuyer-CMRG streams, created by the
# future.apply package). So the results are the same whether the macro runs
# are executed one after another or in parallel, on any number of cores.
################################################################################

# Source every .R file in the folder R/ (the library).
source_library <- function(dir = "R") {
  files <- sort(list.files(dir, pattern = "\\.R$", full.names = TRUE))
  for (file in files) {
    source(file)
  }
}

# Choose how many R processes run the macro replications.
# n_cores = 1 runs everything in the current R session.
setup_parallel <- function(n_cores = 1) {
  if (n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }
}

# Run one_macro_run(1), one_macro_run(2), ..., one_macro_run(n_macro) and
# return their results as a list (element k = result of macro run k).
#
# one_macro_run must be a function of the macro-run index k only. Each call
# gets its own random-number stream, derived from `seed`, so the results do
# not depend on the number of cores. The work is sent to the workers in
# about four chunks per worker.
run_macros <- function(n_macro, one_macro_run, seed) {
  n_workers <- future::nbrOfWorkers()
  chunk_size <- max(1, ceiling(n_macro / (4 * n_workers)))
  future.apply::future_lapply(
    seq_len(n_macro), one_macro_run,
    future.seed = as.integer(seed),
    future.chunk.size = chunk_size
  )
}

# Wall-clock time: start <- Sys.time(); ...; seconds_since(start)
seconds_since <- function(start) {
  as.numeric(difftime(Sys.time(), start, units = "secs"))
}

# Stack a list of data frames (with the same columns) into one data frame.
stack_rows <- function(list_of_data_frames) {
  result <- do.call(rbind, list_of_data_frames)
  rownames(result) <- NULL
  result
}

# Standard error of the mean of x.
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}
