# ----------------------------------------------------------------------------
# TISCA v2 -- validate.R
# Input-rejection rules (spec Section 8, IJDA #11e). The implementation must
# REFUSE metric structures it cannot validly test. All validators return a list
# of booleans + messages so the caller can aggregate and emit a hard stop.
# ----------------------------------------------------------------------------

## ---- rejection rules ----------------------------------------------------------

#' Validate a contrast request against the spec's rejection rules.
#'
#' @param metric_a name (or column) of method A's metric
#' @param metric_b name (or column) of method B's metric
#' @param metric_type kind of metric being contrasted: one of
#'        "loss" (lower better, scalar) | "binary" (0/1, raw) |
#'        "bounded_prop" (row 5, [0,1], raw) | "calibration" (row 8) |
#'        "rmse_ate" (row 3b, derived) | "runtime" (row 11)
#' @param lower_is_better logical: how the caller intends to treat the metric
#' @param mode planned/test hypothesis mode ("M1".."M5")
#' @param data optional data frame to check column presence and values
#' @return list(pass = logical, problems = character vector, infos = character)
validate_contrast <- function(metric_a, metric_b, metric_type,
                              lower_is_better = TRUE,
                              mode = "M1",
                              data = NULL) {
  problems <- character(0)
  infos <- character(0)
  pass <- TRUE

  # Rule 1: across-replication aggregate passed as per-replication value.
  # Detect common aggregate names; the caller must hand the per-replication column.
  agg_markers <- c("_mean$", "_avg$", "mean_", "avg_", "^mean$", "^avg$",
                   "^rmse_ate$", "replicate_means")
  agg_pat <- paste(agg_markers, collapse = "|")
  is_agg <- function(nm) grepl(agg_pat, nm, ignore.case = TRUE)
  if (is_agg(metric_a) || is_agg(metric_b)) {
    problems <- c(problems,
      sprintf("Rule 1: %s looks like an across-replication aggregate; require the per-replication column",
              ifelse(is_agg(metric_a), metric_a, metric_b)))
  }

  # Rule 2: per-replication "RMSE_ATE" is not a valid per-replication scalar.
  # Route through the squared-ATE-error column (row 3a) and re-root (row 3b).
  if (metric_type == "rmse_ate") {
    problems <- c(problems,
      "Rule 2: a per-replication 'RMSE_ATE' scalar is not well-defined; contrast on the squared-ATE-error column (row 3a) and re-root (row 3b).")
  }

  # Rule 3: a metric with a nominal target must not be run as "lower is better".
  if (metric_type %in% c("binary","bounded_prop","calibration") && isTRUE(lower_is_better)) {
    problems <- c(problems,
      "Rule 3: coverage / binary / bounded metrics have a nominal target; they must not be tested lower-is-better. Use the calibration deviation (row 8) or the interval score (row 9).")
  }

  # Rule 4: mode/level mismatch — the planned mode must equal the test mode (C4).
  # (The matched-mode check lives in procedure(); here we record it as the caller
  #  declares a single `mode` used for both planning and the final test.)

  # Rule 6: M5 with |delta| >= Delta is infeasible (signalled at planning time in
  # solve_power_J; here we only note the constraint dependency).
  if (mode == "M5") {
    infos <- c(infos,
      "M5 requires a margin Delta; genuine infeasibility is |delta| >= Delta (reported by the planner, not inferred from a small-J zero-power eval).")
  }

  # Data-level checks
  if (!is.null(data)) {
    for (nm in c(metric_a, metric_b)) {
      if (!nm %in% colnames(data)) {
        problems <- c(problems, sprintf("Column '%s' not found in data", nm))
      }
    }
    present <- c(metric_a, metric_b) %in% colnames(data)
    if (all(present)) {
      va <- data[[metric_a]]; vb <- data[[metric_b]]
      # per-replication finite length
      if (metric_type == "binary") {
        if (!all(va %in% c(0,1,NA)) || !all(vb %in% c(0,1,NA)))
          problems <- c(problems, "binary metric columns must be 0/1")
      } else if (metric_type == "bounded_prop") {
        if (!all(is.na(va) | (va >= 0 & va <= 1)) || !all(is.na(vb) | (vb >= 0 & vb <= 1)))
          problems <- c(problems, "bounded proportion columns must be in [0,1]")
      }
    }
  }

  list(pass = length(problems) == 0, problems = problems, infos = infos)
}

#' Validate a family of contrasts (all primary claims) before any planning.
#'
#' @param contrasts list; each element a list(metric_a, metric_b, metric_type,
#'        lower_is_better, mode, delta, Delta = NULL)
#' @param data optional data frame
#' @param K_expected expected family size
#' @return list(pass, problems, n_contrasts)
validate_family <- function(contrasts, data = NULL, K_expected = NULL) {
  problems <- character(0)
  for (i in seq_along(contrasts)) {
    ctr <- contrasts[[i]]
    res <- validate_contrast(ctr$metric_a, ctr$metric_b,
                             ctr$metric_type %||% "loss",
                             ctr$lower_is_better %||% TRUE,
                             ctr$mode %||% "M1",
                             data = data)
    if (!res$pass) {
      problems <- c(problems, sprintf("contrast %d: %s", i, paste(res$problems, collapse = "; ")))
    }
  }
  if (!is.null(K_expected) && length(contrasts) != K_expected) {
    problems <- c(problems, sprintf("family size %d does not match expected K = %d",
                                    length(contrasts), K_expected))
  }
  list(pass = length(problems) == 0, problems = problems,
       n_contrasts = length(contrasts))
}
# `%||%` is defined once in tisca.R and available package-wide.
