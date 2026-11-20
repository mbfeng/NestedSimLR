################################################################################
# run_tests.R -- check that the library still gives the reference results
#
#   Rscript tests/run_tests.R
#
# Runs every library function on small inputs with fixed seeds (see
# tests/make_reference_outputs.R) and compares the outputs with
# tests/reference_outputs.rds. Run it after changing any file in R/. A
# difference means the change alters results; if that is intended, remake
# the reference with  Rscript tests/make_reference_outputs.R  and say why in
# the README's change log.
################################################################################

source(file.path("R", "utilities.R"))
source_library("R")
source(file.path("tests", "make_reference_outputs.R"))

reference <- readRDS(file.path("tests", "reference_outputs.rds"))
current <- make_reference_outputs()

n_failed <- 0
for (name in names(reference)) {
  check <- all.equal(reference[[name]], current[[name]], tolerance = 1e-12)
  if (isTRUE(check)) {
    cat(sprintf("ok      %s\n", name))
  } else {
    n_failed <- n_failed + 1
    cat(sprintf("FAILED  %s: %s\n", name, paste(check, collapse = "; ")))
  }
}
cat(sprintf("\n%d of %d checks passed\n", length(reference) - n_failed,
            length(reference)))
if (n_failed > 0) {
  quit(status = 1)
}
