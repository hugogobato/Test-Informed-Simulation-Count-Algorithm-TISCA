# ----------------------------------------------------------------------------
# TISCA v2 -- estimands.R
# Per-replication estimands from the reference estimand table (docs/estimand_table.md,
# P1-T1 rev. 2) and the E3 per-replication metric columns. This module turns raw
# per-unit or per-replication data into the replicate-level loss L_j on which
# paired contrasts are built.
#
# Rows implemented here map to the estimand table as follows (nominal level
# parameterised per Section 3.2):
#   1  PEHE (rooted):            pehe(tau_hat, tau)
#   2  CATE MSE (unrooted):      cate_mse(tau_hat, tau)
#   3a ATE squared error:        ate_sq_err(ate_hat, ate_true)
#   3b RMSE_ATE estimand + MCSE: rmse_ate(mcse = ...)  (derived; see row 3b)
#   4  ATE coverage indicator:   ate_cov(ate_hat, ate_true, ci_lo, ci_hi)
#   5  CATE coverage proportion: cate_cov(tau_hat, tau, ci_lo_i, ci_hi_i)
#   6  CIL (width):              cil(ci_lo, ci_hi)
#   7  Bias (CATE):              cate_bias(tau_hat, tau)
#   8  Calibration deviation flow is handled at analysis time (apply abs to the
#      estimate, NEVER per replication); see Section 3.1 of the estimand table.
#   9  Interval score:           interval_score(x, lo, hi, c)
#   10 CRPS:                     crps(...)  -- supplied externally (scoringRules);
#      this module only carries the column as a loss L_j.
#   11 Runtime:                  passed through as a replicate-level column.
#   12 Convergence/failure flag: passed through as a binary column.
# ----------------------------------------------------------------------------

## ---- unit-level / replicate-level estimators --------------------------------

#' Rooted PEHE: PEHE_j = sqrt(mean_i (tau_hat_{j,i} - tau_{j,i})^2).
#' @param tau_hat vector of per-unit CATE estimates on one replication
#' @param tau true per-unit CATE on the same replication
#' @return scalar rooted PEHE for the replication
pehe <- function(tau_hat, tau) {
  sqrt(mean((tau_hat - tau)^2, na.rm = FALSE))
}

#' Unrooted CATE MSE: Q_j = mean_i (tau_hat_{j,i} - tau_{j,i})^2.
#' @inheritParams pehe
#' @return scalar unrooted MSE for the replication
cate_mse <- function(tau_hat, tau) {
  mean((tau_hat - tau)^2, na.rm = FALSE)
}

#' ATE squared error: Q_{ATE,j} = (ate_hat - ate_true)^2.
#' @param ate_hat estimated ATE for the replication
#' @param ate_true true ATE for the replication
#' @return scalar squared error
ate_sq_err <- function(ate_hat, ate_true) {
  (ate_hat - ate_true)^2
}

#' RMSE_ATE estimand is sqrt(E[Q_ATE]); this computes the estimator and its
#' delta-method MCSE from the raw per-replication squared errors (row 3b).
#'
#' The MCSE on the rooted scale is s_Q / (2*sqrt(J)*sqrt(mean(Q))). This is the
#' corrected formula; dropping the inner root inflates/understates by 2x-5x in the
#' PEHE-scale regime (estimand table Section 3, row 3b, verified numerically).
#' The recommended contrast path is to pair-test on the squared scale (row 3a)
#' and re-root the estimate, so this routine is provided for the plain report.
#'
#' @param q_ate vector of per-replication squared ATE errors
#' @return list(rmse = sqrt(mean(q)), mcse = delta-method MCSE)
rmse_ate <- function(q_ate) {
  mq <- mean(q_ate)
  if (length(q_ate) < 2) stop("need >= 2 replications for an MCSE")
  s <- stats::sd(q_ate)
  list(
    rmse  = sqrt(mq),
    mcse  = if (mq > 0) s / (2 * sqrt(length(q_ate)) * sqrt(mq)) else 0
  )
}

#' ATE coverage indicator (row 4): 1(ate_true in [lo, hi]).
#' @param ate_true true ATE
#' @param ci_lo lower credible/CI bound on ATE
#' @param ci_hi upper credible/CI bound on ATE
#' @return 0/1 indicator
ate_cov <- function(ate_true, ci_lo, ci_hi) {
  as.integer(ate_true >= ci_lo & ate_true <= ci_hi)
}

#' CATE coverage proportion (row 5): mean_i 1(tau_i in [lo_i, hi_i]).
#' The inputs are parallel per-unit vectors; a correctly calibrated midpoint of
#' the interval is not assumed, the bounds are used as given.
#' @param tau_i true per-unit CATE vector
#' @param ci_lo_i per-unit lower bound
#' @param ci_hi_i per-unit upper bound
#' @return scalar proportion covered in [0,1]
cate_cov <- function(tau_i, ci_lo_i, ci_hi_i) {
  n <- length(tau_i)
  if (n == 0) return(NA_real_)
  mean(tau_i >= ci_lo_i & tau_i <= ci_hi_i, na.rm = FALSE)
}

#' Mean credible/CI width (row 6): mean_i (hi_i - lo_i) [CATE], or (hi-lo) [ATE].
#' @param ci_lo numeric vector (per-unit or scalar)
#' @param ci_hi numeric vector (per-unit or scalar)
#' @return scalar mean width
cil <- function(ci_lo, ci_hi) {
  mean(ci_hi - ci_lo, na.rm = FALSE)
}

#' CATE bias (row 7): mean_i (tau_hat_i - tau_i).
#' @param tau_hat per-unit CATE estimates
#' @param tau per-unit true CATEs
#' @return scalar bias
cate_bias <- function(tau_hat, tau) {
  mean(tau_hat - tau, na.rm = FALSE)
}

#' Interval score (row 9), Winkler form at nominal level 1-c.
#' IS = (u-l) + (2/c)(l - x)1(x<l) + (2/c)(x - u)1(x>u).
#' Vectorised over units; returns per-unit vector. Lower is better.
#' @param x observed value (may be a vector)
#' @param lo lower bound (vector)
#' @param hi upper bound (vector)
#' @param c nominal miscoverage of the interval (0.05 for the 95% level, 0.50 for the 50% level)
#' @return numeric vector of per-unit interval scores
interval_score <- function(x, lo, hi, c = 0.05) {
  term <- (hi - lo) +
          (2 / c) * (lo - x) * (x < lo) +
          (2 / c) * (x - hi) * (x > hi)
  term
}

#' Mean interval score over units (the scalar D-loss the tables use).
#' @param x observed value (vector)
#' @param lo lower bound (vector)
#' @param hi upper bound (vector)
#' @param c nominal miscoverage
#' @return scalar mean interval score for the replication
mean_interval_score <- function(x, lo, hi, c = 0.05) {
  mean(interval_score(x, lo, hi, c = c), na.rm = FALSE)
}

## ---- validation helpers (used by validate.R) --------------------------------

#' Check a candidate L_j vector is a valid per-replication loss for a contrast.
#' @param vec candidate replicate-level column
#' @return TRUE if length >= 2 and all finite
is_perrep_loss <- function(vec) {
  length(vec) >= 2 && all(is.finite(vec))
}

#' Check a column is a valid binary (0/1) replicate-level indicator (rows 4, 12).
#' @param vec candidate indicator column
#' @return TRUE if all values in {0,1} and length >= 1
is_binary_indicator <- function(vec) {
  length(vec) >= 1 && all(vec %in% c(0, 1)) && all(is.finite(vec))
}

#' Check a bounded proportion column (row 5), 0 <= L_j <= 1.
#' @param vec candidate proportion column
#' @return TRUE if all values in [0,1] and length >= 2
is_bounded_proportion <- function(vec) {
  length(vec) >= 2 && all(is.finite(vec)) &&
    all(vec >= 0 & vec <= 1)
}
