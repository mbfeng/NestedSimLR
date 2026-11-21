################################################################################
# utilities.R                                                                  #
# Small helpers used by all scripts                                            #
# -- loading the library (all files in R/)                                     #
# -- running macro replications reproducibly, serially or in parallel          #
# -- timing, stacking data frames, standard errors                             #
# Reproducibility: macro run k always uses the k-th random-number stream       #
# derived from the experiment's seed (L'Ecuyer-CMRG streams of the             #
# future.apply package), so the results are the same on any number of cores    #
################################################################################

source_library <- function(dir = "R") {
  # Source every .R file in the folder dir (the library)
  # Input:
  # -- dir: folder of the library
  # Output: none (the functions are defined in the calling session)

  files <- sort(list.files(dir, pattern = "\\.R$", full.names = TRUE))
  for (file in files) {
    source(file)
  }
}

setup_parallel <- function(n_cores = 1) {
  # Choose how many R processes run the macro replications
  # Input:
  # -- n_cores: number of worker processes, 1 = the current R session
  # Output: none

  if (n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }
}

run_macros <- function(n_macro, one_macro_run, seed) {
  # Run one_macro_run(1), ..., one_macro_run(n_macro), serially or in
  # parallel, each with its own random-number stream derived from seed
  # Inputs:
  # -- n_macro: number of macro replications
  # -- one_macro_run: function of the macro index k only
  # -- seed: integer seed of the experiment
  # Output:
  # -- list, element k = result of macro run k

  n_workers <- future::nbrOfWorkers() # number of parallel workers
  chunk_size <- max(1, ceiling(n_macro / (4 * n_workers))) # ~4 chunks/worker
  results <- future.apply::future_lapply(
    1:n_macro, one_macro_run,
    future.seed = as.integer(seed),
    future.chunk.size = chunk_size
  )

  return(results)
}

seconds_since <- function(start) {
  # Wall-clock seconds since start (start <- Sys.time())

  return(as.numeric(difftime(Sys.time(), start, units = "secs")))
}

stack_rows <- function(list_of_data_frames) {
  # Stack a list of data frames with the same columns into one data frame

  result <- do.call(rbind, list_of_data_frames)
  rownames(result) <- NULL

  return(result)
}

standard_error <- function(x) {
  # Standard error of the mean of x

  return(sd(x) / sqrt(length(x)))
}
