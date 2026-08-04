# ----------------------------------------------------------------------------
# TISCA v2 R test helper.
# Sources every package module into the global environment so tests can run
# without a full R CMD INSTALL (works from a checkout). Also defines the
# module-loading used by the parity harness.
# ----------------------------------------------------------------------------

.tisca_R_dir <- function() {
  # Called from tisca/ (getwd()/R) or tisca/tests/ (getwd()/../R).
  here <- normalizePath(getwd())
  candidates <- c(file.path(here, "R"), file.path(dirname(here), "R"))
  for (cand in candidates) {
    if (dir.exists(cand)) return(cand)
  }
  # fall back to package install path if NOT running from a checkout
  pkg_path <- system.file("R", package = "tisca")
  if (nzchar(pkg_path)) return(pkg_path)
  stop("cannot locate the tisca/R module directory (looked under ", here, ")")
}

#' Source all TISCA v2 R modules. Idempotent within a session.
#' @param verbose whether to report each loaded module
#' @return invisible character vector of sourced files
source_tisca_R <- function(verbose = FALSE) {
  mods <- c("tisca.R", "estimands.R", "inference.R", "mcs.R",
            "multiplicity.R", "planning.R", "procedure.R", "validate.R")
  dir <- .tisca_R_dir()
  for (m in mods) {
    f <- file.path(dir, m)
    if (!file.exists(f)) stop("missing module file: ", f)
    source(f, local = FALSE)
    if (verbose) message("sourced ", m)
  }
  invisible(mods)
}

#' Assert a value is within an inclusive interval.
#' @param object value to check
#' @param lo lower bound
#' @param hi upper bound
expect_between <- function(object, lo, hi) {
  testthat::expect_true(object >= lo & object <= hi,
                        info = paste0("expected ", object, " in [", lo, ", ", hi, "]"))
}
