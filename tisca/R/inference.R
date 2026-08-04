# ----------------------------------------------------------------------------
# TISCA v2 -- inference.R
# Paired inference on the contrast D_j = L_{A,j} - L_{B,j}.
#
# C1 (plan Section 1.2 / spec Section 7.1): all inference on paired contrasts
# from common replications with common random numbers. NAs are dropped LISTWISE
# across the pair, with the count reported and a sensitivity note available.
# Per-column NA dropping is forbidden (v1 defect, plan Section 1.6).
#
# Routines provided:
#   contrast_from_columns(A, B)        -> list(D, dropped, n_kept)  listwise deletion
#   paired_t(D, alpha)                 -> one-sample t on D, estimate-first output
#   studentized_paired_bootstrap(D, B, alpha, mcse) -> bootstrap CI/rejection
#   paired_bootstrap_bca(D, B, alpha)  -> BCa CIs (auxiliary)
#   mcnemar_exact(x_A, x_B, alpha)     -> exact/binomial McNemar on paired 0/1 indicators
# ----------------------------------------------------------------------------

## ---- listwise deletion ------------------------------------------------------

#' Build the paired contrast with LISTWISE NA deletion across the pair.
#'
#' A replication j is kept for the contrast (A,B) only if BOTH L_{A,j} and L_{B,j}
#' are present and finite. The number dropped is reported. Per-column dropping is
#' not performed anywhere in this module.
#'
#' @param A numeric vector of replication-level losses for method A
#' @param B numeric vector of replication-level losses for method B
#' @return list(D = D_j, dropped = number dropped, n_kept = length(D),
#'              na_A = count of NAs in A, na_B = count of NAs in B,
#'              paired = TRUE)
contrast_from_columns <- function(A, B) {
  if (length(A) != length(B)) {
    stop("contrast_from_columns: A and B must be parallel (same replication index)")
  }
  ok   <- stats::complete.cases(A, B)
  D    <- A[ok] - B[ok]
  list(
    D       = D,
    dropped = sum(!ok),
    n_kept  = length(D),
    na_A    = sum(is.na(A) | is.nan(A)),
    na_B    = sum(is.na(B) | is.nan(B)),
    paired  = TRUE
  )
}

## ---- paired one-sample t ----------------------------------------------------

#' Paired one-sample t-test on D_j (df = J-1), the default fast inference.
#'
#' Returns estimate and Monte Carlo CI first and the p-value second, per IJDA #1e
#' (estimate-first reporting). Correctly removes the `paired=FALSE` defect of v1.
#'
#' @param D numeric vector of paired contrasts
#' @param alpha level (default 0.05)
#' @param mu null hypothesis value of E[D] (default 0)
#' @param alternative "two.sided","less","greater" (must match the planned mode)
#' @return list(estimate, mcse, ci_lower, ci_upper, t_stat, p_value, df, n)
paired_t <- function(D, alpha = 0.05, mu = 0, alternative = c("two.sided","less","greater")) {
  alternative <- match.arg(alternative)
  n <- length(D)
  if (n < 2) stop("paired_t: need at least 2 replications")
  dbar <- mean(D)
  s    <- stats::sd(D)
  mcse <- s / sqrt(n)
  tstat <- (dbar - mu) / mcse
  df <- n - 1
  p <- switch(alternative,
    two.sided = 2 * stats::pt(-abs(tstat), df),
    less      = stats::pt(tstat, df),
    greater   = stats::pt(tstat, df, lower.tail = FALSE)
  )
  crit <- stats::qt(1 - alpha / 2, df)
  half <- crit * mcse
  list(
    estimate   = dbar,
    mcse       = mcse,
    ci_lower   = dbar - half,
    ci_upper   = dbar + half,
    t_stat     = tstat,
    p_value    = p,
    df         = df,
    n          = n,
    method     = "paired t",
    dropped    = NA_integer_
  )
}

## ---- studentized paired bootstrap ------------------------------------------

#' Studentized paired bootstrap on D_j (recommended check for skewed D).
#'
#' The test statistic is T = (Dbar - mu)/se(D). The bootstrap re-estimates
#' Dbar_b and se_b for each resample and accumulates the studentized pivot
#' (Dbar_b - Dbar)/se_b. The reported CI and p-value are based on that pivot,
#' which is second-order accurate and handles the skew/long tails that break
#' the CLT-based t (P3-T4). Report the bootstrap MCSE of the p-value/CI bound.
#'
#' @param D numeric vector of paired contrasts
#' @param B bootstrap resamples (default 4999)
#' @param alpha level
#' @param mu null value
#' @param alternative sidedness
#' @param seed optional RNG seed for reproducibility
#' @param margin margin pass-through for the combined statistic (unused here, kept for API parity)
#' @return list(estimate, ci_lower, ci_upper, p_value, pivot_mcse, B, seed)
studentized_paired_bootstrap <- function(D, B = 4999L, alpha = 0.05, mu = 0,
                                         alternative = c("two.sided","less","greater"),
                                         seed = NULL, margin = NULL) {
  alternative <- match.arg(alternative)
  n <- length(D)
  if (n < 2) stop("studentized_paired_bootstrap: need >= 2 replications")
  if (!is.null(seed)) set.seed(seed)
  dbar <- mean(D)
  se   <- stats::sd(D) / sqrt(n)
  if (se <= 0) {
    # degenerate deterministic contrast
    p <- if (alternative == "two.sided") (dbar != mu) else (if (alternative == "less") dbar < mu else dbar > mu)
    return(list(estimate = dbar, ci_lower = dbar, ci_upper = dbar,
                p_value = as.numeric(p), pivot_mcse = 0, B = B, seed = seed,
                degenerate = TRUE))
  }
  # fast bootstrap resampling
  idx <- matrix(sample.int(n, size = n * B, replace = TRUE), nrow = B, ncol = n)
  db  <- rowMeans(matrix(D[idx], nrow = B, ncol = n))
  vb  <- apply(matrix(D[idx], nrow = B, ncol = n), 1,
               function(x) stats::sd(x) / sqrt(length(x)))
  piv <- (db - dbar) / vb
  piv <- piv[is.finite(piv)]
  B_eff <- length(piv)
  if (B_eff == 0) stop("studentized_paired_bootstrap: no finite bootstrap pivots")

  lv <- stats::quantile(piv, alpha / 2, type = 6)
  uv <- stats::quantile(piv, 1 - alpha / 2, type = 6)
  ci_lower <- dbar - uv * se
  ci_upper <- dbar - lv * se

  # p-value: fraction of pivots on the far side of the observed pivot (0 in pivot units)
  obs <- (mu - dbar) / se  # pivot centred on mu for the null
  pv <- switch(alternative,
    two.sided = mean(abs(piv) >= abs(obs)),
    less      = mean(piv <= obs),
    greater   = mean(piv >= obs)
  )
  # bootstrap MCSE of the p-value (rate estimate)
  p_mcse <- sqrt(pv * (1 - pv) / B_eff)
  list(estimate = dbar, ci_lower = ci_lower, ci_upper = ci_upper,
       p_value = pv, pivot_mcse = p_mcse, B = B, seed = seed,
       degenerate = FALSE, n_kept = n)
}

#' BCa bootstrap CIs on the mean (auxiliary; not the primary studentized path).
#' @param D numeric vector of paired contrasts
#' @param B resamples
#' @param alpha level
#' @param seed optional seed
#' @return list(estimate, ci_lower, ci_upper)
paired_bootstrap_bca <- function(D, B = 4999L, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  boot_mean <- function(data) {
    res <- replicate(B, mean(sample(data, replace = TRUE)))
    # BCa
    z0 <- stats::qnorm(mean(res < mean(data)))
    jack <- numeric(length(data))
    for (i in seq_along(data)) jack[i] <- mean(data[-i])
    a <- sum((mean(data) - jack)^3) / (6 * (sum((mean(data) - jack)^2))^1.5)
    z1 <- z0 + z0 + stats::qnorm(1 - alpha / 2)
    za <- z0 + (z0 + stats::qnorm(1 - alpha / 2)) / (1 - a * (z0 + stats::qnorm(1 - alpha / 2)))
    p_low <- stats::pnorm(z0 + (z0 - stats::qnorm(1 - alpha / 2)) / (1 - a * (z0 - stats::qnorm(1 - alpha / 2))))
    list(lo = stats::quantile(res, p_low, type = 6),
         hi = stats::quantile(res, 1 - p_low, type = 6),
         est = mean(data))
  }
  r <- boot_mean(D)
  list(estimate = r$est, ci_lower = r$lo[[1]], ci_upper = r$hi[[1]])
}

## ---- McNemar for paired binary rows (4, 12) ----------------------------------

#' Exact/binomial McNemar test on paired binary indicators.
#'
#' For two binary columns x_A, x_B (e.g. ATE coverage indicators), only the
#' discordant pairs (1,0) and (0,1) carry information. The test is the sign test
#' on the discordant pairs: b = count of (A=1,B=0), c = count of (A=0,B=1),
#' H0: b/(b+c) = 0.5, with an exact binomial two-sided p-value (R's binom.test).
#' A mid-p correction is optionally applied (small-sample exact calibration).
#'
#' @param x_A 0/1 replicate-level indicator for A
#' @param x_B 0/1 replicate-level indicator for B
#' @param alpha level
#' @param alternative sidedness of the difference in coverage
#' @param midp apply the mid-p correction
#' @return list(b, c, n_discordant, p_value, estimate, ci_lower, ci_upper, method)
mcnemar_exact <- function(x_A, x_B, alpha = 0.05,
                          alternative = c("two.sided","less","greater"),
                          midp = FALSE) {
  alternative <- match.arg(alternative)
  ok  <- stats::complete.cases(x_A, x_B)
  x_A <- x_A[ok]; x_B <- x_B[ok]
  b <- sum(x_A == 1 & x_B == 0)
  c <- sum(x_A == 0 & x_B == 1)
  n_disc <- b + c
  est <- if (n_disc > 0) b / n_disc else NA_real_
  if (n_disc == 0) {
    return(list(b = b, c = c, n_discordant = 0, p_value = 1,
                estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                method = "McNemar exact (no discordant pairs)"))
  }
  # p-value under H0: X ~ Bin(n_disc, 0.5). Two-sided p via doubling the smaller tail.
  # (Also handles mid-p correction by subtracting half the point mass at b.)
  pv <- switch(alternative,
    two.sided = 2 * pmin(stats::pbinom(b, n_disc, 0.5),
                         stats::pbinom(b - 1, n_disc, 0.5, lower.tail = FALSE)),
    less      = stats::pbinom(b, n_disc, 0.5),
    greater   = stats::pbinom(b - 1, n_disc, 0.5, lower.tail = FALSE)
  )
  if (alternative == "two.sided") pv <- pmin(pv, 1)  # clamp the doubled tail
  if (midp) {
    pmf <- stats::dbinom(b, n_disc, 0.5)
    pv  <- pmin(pmax(pv - 0.5 * pmf, 0), 1)
  }
  ci <- stats::binom.test(b, n_disc, 0.5, conf.level = 1 - alpha)$conf.int
  list(b = b, c = c, n_discordant = n_disc, p_value = pv,
       estimate = est, ci_lower = ci[[1]], ci_upper = ci[[2]],
       dropped = sum(!ok), method = "McNemar exact (sign test on discordant pairs)")
}
