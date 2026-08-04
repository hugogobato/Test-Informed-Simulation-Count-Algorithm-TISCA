# ----------------------------------------------------------------------------
# TISCA v2 -- planning.R
# Design/planning layer: hypothesis modes M1-M5, precision and power targets,
# variance-uncertainty propagation, and the J solver.
#
# Reference spec: docs/tisca_v2_spec.md (P1-T2 rev. 2). All formulas here
# inherit that document's notation verbatim. This is the R port and must stay
# in numeric parity with the Python reference (tisca/python/tisca/planning.py).
#
# Conventions:
#   D_j = L_{A,j} - L_{B,j}   (lower is better), theta = E[D_j], sigma_D^2 = Var(D_j)
#   delta = planning alternative, Delta = margin, alpha = (possibly adjusted) level
#   ncp = sqrt(J) * (theta_shifted) / sigma_D    (noncentral t with J-1 df)
# ----------------------------------------------------------------------------

## ---- helpers ---------------------------------------------------------------

#' Noncentral-t CDF power relation.
#'
#' In R `pt(q, df, ncp)` is already the correct noncentral-t CDF (unlike the v1
#' SciPy misuse flagged in the plan's Section 1.6), so we keep `pt` directly and
#' do NOT reimplement it. This wrapper exists only to name the primitive and to
#' be the single call site the Python parity harness compares against.
#'
#' @inheritParams stats::pt
#' @return vector of probabilities
pt_dncp <- function(q, df, ncp = 0, lower.tail = TRUE) {
  stats::pt(q = q, df = df, ncp = ncp, lower.tail = lower.tail)
}

## ---- hypothesis modes --------------------------------------------------------

#' M1 -- two-sided equality power (H0: theta = 0 vs H1: theta != 0).
#' @param J replication count
#' @param delta planning alternative (value of theta at which power is evaluated)
#' @param sigma_D residual std dev of the contrast (after any inflation)
#' @param alpha level (may be adjusted)
#' @return power in [0,1]
power_m1 <- function(J, delta, sigma_D, alpha) {
  if (length(J) == 0) return(numeric(0))
  crit <- stats::qt(1 - alpha / 2, J - 1)
  ncp  <- sqrt(J) * delta / sigma_D
  1 - stats::pt(crit, df = J - 1, ncp = ncp) +
      stats::pt(-crit, df = J - 1, ncp = ncp)
}

#' M2 -- directional superiority, lower is better (H0: theta >= 0 vs H1: theta < 0).
#' Reject iff T < t_{alpha, J-1}.
#' @param J replication count
#' @param delta planning alternative (delta < 0)
#' @param sigma_D residual std dev of the contrast
#' @param alpha level
#' @return power in [0,1]
power_m2 <- function(J, delta, sigma_D, alpha) {
  crit <- stats::qt(alpha, J - 1)
  ncp  <- sqrt(J) * delta / sigma_D
  stats::pt(crit, df = J - 1, ncp = ncp)
}

#' M3 -- minimum-effect superiority (H0: theta >= -Delta vs H1: theta < -Delta).
#' Reject iff T < t_{alpha, J-1} with T centered on the inner boundary -Delta.
#' @param J replication count
#' @param delta planning alternative, delta < -Delta
#' @param sigma_D residual std dev of the contrast
#' @param alpha level
#' @param Delta margin (must be > 0)
#' @return power in [0,1]
power_m3 <- function(J, delta, sigma_D, alpha, Delta) {
  crit <- stats::qt(alpha, J - 1)
  ncp  <- sqrt(J) * (delta + Delta) / sigma_D
  stats::pt(crit, df = J - 1, ncp = ncp)
}

#' M4 -- non-inferiority (H0: theta >= Delta vs H1: theta < Delta).
#' Reject iff T' < t_{alpha, J-1}, T' = sqrt(J)*(Dbar - Delta)/sigma.
#' @param J replication count
#' @param delta planning alternative, delta < Delta
#' @param sigma_D residual std dev of the contrast
#' @param alpha level
#' @param Delta margin (must be > 0)
#' @return power in [0,1]
power_m4 <- function(J, delta, sigma_D, alpha, Delta) {
  crit <- stats::qt(alpha, J - 1)
  ncp  <- sqrt(J) * (delta - Delta) / sigma_D
  stats::pt(crit, df = J - 1, ncp = ncp)
}

#' M5 -- equivalence via two one-sided tests (TOST), known-sigma planning form.
#'
#' Power of the combined TOST rejection at planning alternative delta with
#' |delta| < Delta, using the known-sigma approximation of the acceptance
#' boundary m = Delta - t_{1-alpha,J-1} * sigma/sqrt(J).
#'
#' For `exact = FALSE` the closed-form central-t difference is returned (spec
#' Section 1 M5). For `exact = TRUE` the quadrature form
#'   E_{s_D}[ Phi((m(s)-delta)sqrt(J)/sigma) - Phi((-m(s)-delta)sqrt(J)/sigma) ],
#'   m(s) = max(Delta - t_{1-alpha,J-1}*s/sqrt(J), 0), s_D ~ sigma*sqrt(chi2/(J-1))
#' is used. The approximation is accurate to <= 0.009 in power over the planning
#' range and returns the same planned J at the 0.80 target; exact is for reporting.
#'
#' @param J replication count
#' @param delta planning alternative, |delta| < Delta
#' @param sigma_D residual std dev of the contrast
#' @param alpha level
#' @param Delta margin (must be > 0)
#' @param exact if TRUE use the quadrature (exact) TOST power
#' @param n_quad number of quadrature nodes (exact form only)
#' @return power in [0,1]
power_m5 <- function(J, delta, sigma_D, alpha, Delta, exact = FALSE, n_quad = 2001L) {
  alpha <- alpha[[1]]; delta <- delta[[1]]; sigma_D <- sigma_D[[1]]
  Delta <- Delta[[1]]
  if (!(is.finite(delta) && is.finite(Delta)) || abs(delta) >= Delta) {
    return(0)
  }
  if (!exact) {
    t <- stats::qt(1 - alpha, J - 1)
    b_lo <- (-Delta + t * sigma_D / sqrt(J) - delta) * sqrt(J) / sigma_D
    b_hi <- ( Delta - t * sigma_D / sqrt(J) - delta) * sqrt(J) / sigma_D
    return(stats::pt(b_hi, df = J - 1) - stats::pt(b_lo, df = J - 1))
  }
  # exact quadrature over s_D ~ sigma_D * sqrt(chi2_{J-1}/(J-1))
  warning("TOST exact power is not in parity scope; use for reporting only.")
  t  <- stats::qt(1 - alpha, J - 1)
  df <- J - 1
  qs <- stats::qchisq((seq_len(n_quad) - 0.5) / n_quad, df = df)
  ss <- sigma_D * sqrt(qs / df)
  m  <- pmax(Delta - t * ss / sqrt(J), 0)
  mean(stats::pnorm((m - delta) * sqrt(J) / sigma_D) -
       stats::pnorm((-m - delta) * sqrt(J) / sigma_D))
}

#' Dispatch power() over the five modes.
#' @param mode one of "M1","M2","M3","M4","M5" (or 1L..5L)
#' @param J replication count
#' @param delta planning alternative
#' @param sigma_D residual std dev of the contrast (post-inflation)
#' @param alpha level (may be adjusted)
#' @param Delta margin, required for M3/M4/M5
#' @param ... passed to power_m5 (e.g. exact)
#' @return power in [0,1]
power <- function(mode, J, delta, sigma_D, alpha, Delta = NULL, ...) {
  mode <- mode[[1]]
  if (is.character(mode)) mode <- match.arg(mode, c("M1","M2","M3","M4","M5"))
  switch(as.character(mode),
    "M1" = power_m1(J, delta, sigma_D, alpha),
    "M2" = power_m2(J, delta, sigma_D, alpha),
    "M3" = {
      if (is.null(Delta)) stop("M3 requires a margin Delta")
      power_m3(J, delta, sigma_D, alpha, Delta)
    },
    "M4" = {
      if (is.null(Delta)) stop("M4 requires a margin Delta")
      power_m4(J, delta, sigma_D, alpha, Delta)
    },
    "M5" = {
      if (is.null(Delta)) stop("M5 requires a margin Delta")
      power_m5(J, delta, sigma_D, alpha, Delta, ...)
    },
    stop("unknown mode ", as.character(mode))
  )
}

## ---- variance-uncertainty propagation (Section 3 of the spec) ----------------

#' One-sided upper confidence bound on sigma_D from a noisy pilot (IJDA #8b).
#'
#' sigma_hat_UB = s_D * sqrt((J0-1)/chi2_{gamma, J0-1}), gamma = 0.20 default.
#' This is an assurance argument (raises Pr(achieved power >= target) to ~1-gamma),
#' NOT a bias correction on E[J]. The pilot row count J0 is the number of
#' complete (listwise-valid) replications used to estimate s_D.
#'
#' @param sd_D pilot sample sd of the contrast
#' @param J0 pilot size used to estimate sd_D
#' @param gamma assurance level (0.20 default: 80% confidence the true sigma is at most the bound)
#' @return inflated sd
sigma_ub <- function(sd_D, J0, gamma = 0.20) {
  chi <- stats::qchisq(gamma, df = J0 - 1)
  if (chi <= 0) return(Inf)
  sd_D * sqrt((J0 - 1) / chi)
}

## ---- solvers ----------------------------------------------------------------

#' Smallest J such that the Monte Carlo SE criterion holds.
#' MCSE target: s_D / sqrt(J) <= m.
#' @param sigma_D std dev of the contrast (post-inflation)
#' @param m absolute MCSE target
#' @param J_min starting search point (>= 2)
#' @param J_max budget cap
#' @return list with J and achieved MCSE
solve_mcse_J <- function(sigma_D, m, J_min = 2L, J_max = 1e6L) {
  if (sigma_D <= 0) {
    # degenerate: no replications needed, target achieved; report the cell (P3-T7 / spec Section 8.5)
    return(list(J = 2L, achieved = 0, degenerate = TRUE))
  }
  if (m <= 0) stop("MCSE target m must be > 0")
  j <- max(J_min, 2L)
  while (j <= J_max && sigma_D / sqrt(j) > m) j <- j + 1L
  if (j > J_max) {
    return(list(J = J_max, achieved = sigma_D / sqrt(J_max), capped = TRUE))
  }
  list(J = j, achieved = sigma_D / sqrt(j), degenerate = FALSE)
}

#' Smallest J such that the CI half-width target holds.
#' t_{1-alpha/2, J-1} * s_D / sqrt(J) <= h. Monotone in J for fixed s_D.
#' The default tolerance is machine-scale so we never accept a value that
#' strictly exceeds the target (parity: returns the same smallest J as Python).
#' @param sigma_D std dev of the contrast (post-inflation)
#' @param h target half-width
#' @param alpha level (for the central 1-alpha interval)
#' @param J_min starting search point (>= 2)
#' @param J_max budget cap
#' @param tolerance stopping tolerance on the half-width (machine-scale default)
#' @return list with J, achieved half-width
solve_halfwidth_J <- function(sigma_D, h, alpha = 0.05, J_min = 2L,
                              J_max = 1e6L, tolerance = 1e-12) {
  if (sigma_D <= 0) {
    return(list(J = 2L, achieved = 0, degenerate = TRUE))
  }
  if (h <= 0) stop("half-width target h must be > 0")
  j <- max(J_min, 2L)
  hw <- function(j) stats::qt(1 - alpha / 2, j - 1) * sigma_D / sqrt(j)
  while (j <= J_max && hw(j) > h + tolerance) j <- j + 1L
  if (j > J_max) {
    return(list(J = J_max, achieved = hw(J_max), capped = TRUE))
  }
  list(J = j, achieved = hw(j), degenerate = FALSE)
}

#' Smallest J achieving a target power at a given mode and adjusted alpha.
#'
#' Power uses the noncentral t with J-1 df and the adjusted level alpha_adj.
#' Increments J until power >= target. Handles the degenerate s_D = 0 cell by
#' reporting target as satisfied (no further replications), and, for M5, keeps
#' increasing J past a zero-power small-J artifact (spec Section 1 M5 feasibility).
#'
#' @param mode one of "M1".."M5"
#' @param delta planning alternative
#' @param sigma_D std dev of the contrast (post-inflation)
#' @param alpha_adj adjusted level
#' @param target_power desired power (1-beta)
#' @param Delta margin (M3/M4/M5)
#' @param J_min starting search (>= 2)
#' @param J_max budget cap
#' @param ... passed through to power() (e.g. exact for M5)
#' @return list(J, achieved_power, capped, degenerate)
solve_power_J <- function(mode, delta, sigma_D, alpha_adj, target_power = 0.80,
                          Delta = NULL, J_min = 2L, J_max = 1e6L, ...) {
  sigma_D <- sigma_D[[1]]
  if (!is.finite(sigma_D) || sigma_D <= 0) {
    # degenerate s_D: deterministic contrast, all power targets achieved (spec Section 8.5)
    return(list(J = 2L, achieved_power = 1, capped = FALSE, degenerate = TRUE))
  }
  alpha_adj <- alpha_adj[[1]]; delta <- delta[[1]]; target_power <- target_power[[1]]

  # M5 genuine infeasibility: |delta| >= Delta (spec Section 1 M5, Section 8.6)
  if (mode == "M5" || identical(as.character(mode), "M5")) {
    if (is.null(Delta) || is.finite(Delta) == FALSE) stop("M5 requires a margin Delta")
    if (abs(delta) >= Delta) {
      return(list(J = NA_integer_, achieved_power = NA_real_,
                  capped = FALSE, degenerate = FALSE, infeasible = TRUE,
                  message = "M5 infeasible: |delta| >= Delta, planning alternative is not inside the equivalence margin"))
    }
  }

  j <- max(J_min, 2L)
  p <- suppressWarnings(power(mode, j, delta, sigma_D, alpha_adj, Delta, ...))
  while (j <= J_max && length(p) >= 1L && p < target_power) {
    j <- j + 1L
    p <- suppressWarnings(power(mode, j, delta, sigma_D, alpha_adj, Delta, ...))
  }
  if (j > J_max) {
    return(list(J = J_max, achieved_power = p, capped = TRUE, degenerate = FALSE,
                message = "target not reached within budget J_max"))
  }
  list(J = j, achieved_power = p, capped = FALSE, degenerate = FALSE)
}

#' Combine active design targets and comparisons into a single J_final.
#'
#' J_final = min( max over active targets and comparisons of per-detail J, J_max ).
#' Per the spec Section 2.3, if the closed-form J exceeds J_max we stop at J_max
#' and report the reduced achieved precision/power via the `capped` flag.
#'
#' @param J_details numeric vector of per-comparison/per-target Js
#' @param J_max budget cap
#' @return list(J_final, capped)
combine_J <- function(J_details, J_max) {
  J_final <- max(c(J_details, 0), na.rm = TRUE)
  capped  <- J_final >= J_max
  list(J_final = as.integer(min(J_final, J_max)), capped = capped)
}

# R-visible exports are managed by the package NAMESPACE (roxygen); these
# source files can also be loaded directly with source() for ad-hoc use.
