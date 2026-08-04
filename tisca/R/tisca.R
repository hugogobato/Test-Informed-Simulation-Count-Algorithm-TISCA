# ----------------------------------------------------------------------------
# TISCA v2 -- tisca.R
# Package entry / shared utilities for the R port (P2-T2).
#
# This is the R twin of tisca/python/tisca/__init__.py. It provides:
#   * the shared `%||%` null-coalescing helper (defines it once here; the
#     per-module files use it without redefining),
#   * the version constant,
#   * a top-level `tisca_power()` convenience and a small self-test hook.
#
# Module load order (alphabetical, each defines its own functions):
#   estimands, inference, mcs, multiplicity, planning, procedure, validate
# ----------------------------------------------------------------------------

#' TISCA v2 version.
tisca_version <- "2.0.0-dev"

#' Null-coalescing operator: returns b when a is NULL.
#' @param a value
#' @param b fallback
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Convenience wrapper matching the Python `power()` dispatch, for interactive use
#' and for the R/Python parity tests.
#' @param mode one of "M1".."M5"
#' @param J replication count (can be a vector)
#' @param delta planning alternative
#' @param sigma_D residual std dev of the contrast
#' @param alpha level
#' @param Delta margin (M3/M4/M5)
#' @param ... forwarded to power()
#' @return numeric power
power_ <- function(mode, J, delta, sigma_D, alpha, Delta = NULL, ...) {
  power(mode, J, delta, sigma_D, alpha, Delta = Delta, ...)
}

#' Report which TISCA v2 modules are loadable without their optional deps.
#' @return invisible named logical vector
tisca_ready <- function() {
  mods <- c("estimands","inference","mcs","multiplicity","planning","procedure","validate")
  ok <- vapply(mods, function(m) {
    tryCatch({ source(system.file("R", paste0(m, ".R"), package = "tisca", mustWork = FALSE),
                      local = TRUE); TRUE }, error = function(e) FALSE)
  }, logical(1))
  invisible(ok)
}
