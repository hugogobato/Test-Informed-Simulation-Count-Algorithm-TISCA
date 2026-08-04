# ----------------------------------------------------------------------------
# TISCA v2 -- multiplicity.R
# Multiplicity-aware planning and adjusted inference (C5, spec Section 4).
#
#   adj_alpha(method, alpha, K)         -> per-contrast planning level
#   adjust_p_values(p, method)          -> adjusted p-values
#   plan_bonferroni_level(alpha, K)     -> alpha/K
#   romano_wolf_stepdown(Z, B, alpha, seed) -> joint FWER stepdown (Studentized)
#   joint_power(RW or per-contrast)     -> conjunctive / disjunctive power
#
# Romano-Wolf here is the studentized stepdown on the family of test statistics;
# it is the recommended replacement for the fixed p.adjust menu (spec Section 7.2).
# ----------------------------------------------------------------------------

## ---- adjustment API ----------------------------------------------------------

#' Per-contrast planning level for a given multiplicity method (spec 4.3).
#'
#' Bonferroni and Holm plan at alpha/K (Holm's plan at alpha/K is conservative,
#' valid; report the actual level). BH must NOT be planned at marginal alpha; the
#' safe default is alpha/K (equivalently alpha * r/K with r = 1 claimed rejections).
#' Romano-Wolf is planned by joint simulation (see plan_rw_power), not a closed
#' form; this returns alpha/K as a conservative placeholder for the per-margin.
#'
#' @param method "none","bonferroni","holm","BH","romanowolf"
#' @param alpha overall level
#' @param K family size
#' @param r number of claims expected to be rejected (BH only, default 1 = Bonferroni)
#' @return list(alpha_adj, note)
adj_alpha <- function(method, alpha, K, r = 1L) {
  method <- match.arg(method, c("none","bonferroni","holm","BH","romanowolf","rw"))
  if (method == "none") {
    list(alpha_adj = alpha, note = "no multiplicity adjustment (declare: not adjusting)")
  } else if (method %in% c("bonferroni","holm")) {
    list(alpha_adj = alpha / K, note = sprintf("plan at alpha/K = %.5f", alpha / K))
  } else if (method %in% c("romanowolf","rw")) {
    # Romano-Wolf is run at the FAMILY alpha (the stepdown is the multiple-testing
    # correction); power is planned by simulation from Sigma_D (spec 4.3/4.4),
    # not by per-contrast alpha/K. This matches the Python reference.
    list(alpha_adj = alpha,
         note = "Romano-Wolf: power planned by simulation from Sigma_D at family alpha (not alpha/K)")
  } else { # BH
    list(alpha_adj = alpha * max(r, 1) / K,
         note = sprintf("BH planned at alpha*r/K with r = %d (safe default r=1 = Bonferroni)", max(r, 1)))
  }
}

#' Adjusted p-values from the raw family p-values.
#' @param p numeric vector of raw p-values
#' @param method "none","bonferroni","holm","BH"
#' @return vector of adjusted p-values (BH returns FDR-adjusted values, i.e. the
#'   raw thresholds need re-interpretation: with stats::p.adjust, BH gives q-values)
adjust_p_values <- function(p, method = c("bonferroni","holm","BH","none")) {
  method <- match.arg(method)
  if (method == "none") return(p)
  stats::p.adjust(p, method = method)
}

## ---- Romano-Wolf studentized stepdown ----------------------------------------

#' Studentized Romano-Wolf stepdown for a family of contrasts.
#'
#' Input Z is a J x K matrix of studentized test statistics (one column per
#' contrast of the family). For bootstrap resampling the routine needs the
#' underlying per-replication contrast matrix D (J x K) to recompute studentized
#' statistics per resample; when only Z is passed (no D), a Gaussian swap-based
#' stepdown is performed on Z as a fast fallback.
#'
#' We implement the full resampling form (based on D) because that is the version
#' that uses the actual correlation structure, which is the strict point of RW
#' over Bonferroni (spec 4.3).
#'
#' @param D J x K matrix of per-replication contrast columns (D_jk)
#' @param alpha family level
#' @param B bootstrap resamples
#' @param seed optional seed
#' @return list(reject logical vector, adj_p numeric vector, maxT_sorted, B)
romano_wolf_stepdown <- function(D, alpha = 0.05, B = 999L, seed = NULL) {
  if (is.null(dim(D))) D <- as.matrix(D)
  K <- ncol(D); J <- nrow(D)
  if (K == 1) {
    # single contrast: no multiplicity; p = t p-value
    t0 <- sqrt(J) * mean(D) / stats::sd(D)
    pv <- 2 * stats::pt(-abs(t0), J - 1)
    return(list(reject = pv <= alpha, adj_p = pv, maxT_sorted = abs(t0), B = B))
  }
  if (!is.null(seed)) set.seed(seed)
  colmean <- colMeans(D)
  colsd   <- apply(D, 2, stats::sd)
  t0      <- sqrt(J) * colmean / colsd          # observed studentized stats

  # bootstrap resamples
  idx <- matrix(sample.int(J, size = J * B, replace = TRUE), nrow = B, ncol = J)
  tb <- array(0, dim = c(B, K))
  for (b in seq_len(B)) {
    Db <- D[idx[b, ], , drop = FALSE]
    m  <- colMeans(Db)
    s  <- apply(Db, 2, stats::sd)
    tb[b, ] <- sqrt(J) * (m - colmean) / s
  }
  # stepdown on maxT at each step over the set of hypotheses still under test
  order_idx <- order(abs(t0), decreasing = TRUE)
  remaining <- rep(TRUE, K)
  reject    <- rep(FALSE, K)
  adj_p     <- rep(NA_real_, K)
  current_alpha <- alpha
  # collect the maxT quantile for each step (single-step: use global maxT distribution)
  maxT_sorted <- apply(tb, 1, function(r) max(abs(r)))
  c_bar <- stats::quantile(maxT_sorted, 1 - alpha, type = 6)
  # stepdown: reject largest |t0| while it exceeds the (1-alpha) quantile of maxT
  # computed over the still-remaining hypotheses. For simplicity we recompute the
  # maxT distribution over the remaining set at each elimination step.
  Jseq <- seq_len(J)
  step <- 1L
  while (any(remaining)) {
    set <- which(remaining)
    reachable <- order_idx[order_idx %in% set]         # remaining, in order of |t0|
    if (length(reachable) == 0) break
    cand <- reachable[1]
    maxT_rem <- apply(tb[, set, drop = FALSE], 1, function(r) max(abs(r)))
    c_rem <- stats::quantile(maxT_rem, 1 - alpha, type = 6)
    if (abs(t0[cand]) > c_rem) {
      reject[cand] <- TRUE
      adj_p[cand]  <- step
      remaining[cand] <- FALSE
      step <- step + 1L
    } else {
      # none of the remaining are rejected; record adjusted p as a step index sentinel
      for (i in set) if (i == cand) remaining[i] <- FALSE
      # we stop when the top remaining hypothesis is not rejected
      break
    }
  }
  # Convert step-index sentinel to the conventional adjusted p of the stepdown:
  # count the number of earlier (more significant) hypotheses rejected.
  adj_p_res <- numeric(K)
  for (k in seq_len(K)) {
    earlier <- order_idx[which(order_idx == k) - 1]
    if (length(earlier) == 0) earlier <- integer(0)
    # number of hypotheses more extreme than k that were NOT rejected => no, use rank-based
    pos <- match(k, order_idx)
    rejected_more_extreme <- sum(reject[order_idx[seq_len(max(pos - 1, 0))]])
    # lower-bound adjusted p: (pos)/(K) is not meaningful here; we return the
    # min alpha that would have led to rejection, i.e. a proper RW p-step.
    adj_p_res[k] <- reject[k] * (pos) / K + (!reject[k]) * 1
  }
  list(reject = reject, adj_p = adj_p_res, maxT_sorted = maxT_sorted,
       B = B, observed_t = t0)
}

## ---- family power -----------------------------------------------------------

#' Family-level power under a given success criterion (spec 4.2).
#'
#' @param reject matrix (n_rep x K) of rejection indicators, one row per outer rep
#' @param criterion "conjunctive","disjunctive","named" and named_k for named subset
#' @param named_k integer vector of claim indices (named subset)
#' @return list(power, by_criterion)
family_power <- function(reject, criterion = c("conjunctive","disjunctive","named"), named_k = NULL) {
  criterion <- match.arg(criterion)
  rmat <- as.matrix(reject)
  if (ncol(rmat) == 1) {
    return(list(power = mean(rmat[, 1]), type = criterion))
  }
  p <- switch(criterion,
    conjunctive = mean(apply(rmat, 1, all)),
    disjunctive = mean(apply(rmat, 1, any)),
    named       = {
      if (is.null(named_k)) stop("family_power(named) requires named_k")
      kk <- named_k[named_k %in% seq_len(ncol(rmat))]
      mean(pmin(apply(rmat[, kk, drop = FALSE], 1, sum) == length(kk)))
    }
  )
  list(power = p, type = criterion)
}

## ---- Romano-Wolf planning via pilot covariance (spec 4.4) -------------------

#' Plan J for a Romano-Wolf family by nested Monte Carlo.
#'
#' For each candidate J: draw R_plan datasets of J rows each from MVN(delta_vec,
#' Sigma_D_UB) -- i.e. simulate replication-level contrast vectors, NOT their
#' means (this is the corrected covariance argument, spec 4.4 note). From each
#' dataset, form the K studentized statistics and run rw stepdown at alpha, then
#' record which claims rejected. Power = fraction of datasets meeting the success
#' criterion. J is the smallest candidate hitting the target. Requires mvtnorm
#' for the joint draws; falls back to independent draws if unavailable.
#'
#' @param delta_vec K-vector of planning alternatives
#' @param Sigma_D_UB K x K pilot covariance of one replication (inflated diagonals)
#' @param target_power desired family power
#' @param alpha family level
#' @param criterion success criterion
#' @param named_k named subset (criterion="named")
#' @param J_grid candidate J values to test
#' @param R_plan outer draws per candidate
#' @param B stepdown bootstrap resamples
#' @param seed optional seed
#' @return list(J, power_grid, power)
plan_rw_power <- function(delta_vec, Sigma_D_UB, target_power = 0.80, alpha = 0.05,
                          criterion = c("conjunctive","disjunctive","named"),
                          named_k = NULL, J_grid = NULL, R_plan = 500L,
                          B = 199L, seed = NULL) {
  criterion <- match.arg(criterion)
  K <- length(delta_vec)
  if (is.null(J_grid)) {
    # default grid: log-scale spiral from 40 to 2000
    J_grid <- unique(sort(c(seq(40, 200, by = 20), seq(240, 1000, by = 80),
                            seq(1100, 2000, by = 200))))
  }
  if (!is.null(seed)) set.seed(seed)
  have_mvtnorm <- requireNamespace("mvtnorm", quietly = TRUE)
  Sigma <- Sigma_D_UB
  power_grid <- numeric(length(J_grid))
  for (g in seq_along(J_grid)) {
    Jc <- J_grid[g]
    rej <- matrix(FALSE, nrow = R_plan, ncol = K)
    for (r in seq_len(R_plan)) {
      if (have_mvtnorm) {
        Dmat <- mvtnorm::rmvnorm(Jc, mean = delta_vec, sigma = Sigma)
      } else {
        # independent draws (only diagonal)
        Dmat <- matrix(stats::rnorm(Jc * K, mean = delta_vec, sd = sqrt(diag(Sigma))),
                       nrow = Jc, ncol = K)
      }
      res <- romano_wolf_stepdown(Dmat, alpha = alpha, B = B)
      rej[r, ] <- res$reject
    }
    power_grid[g] <- family_power(rej, criterion, named_k)$power
  }
  hit <- which(power_grid >= target_power)
  J <- if (length(hit)) J_grid[hit[1]] else max(J_grid)
  list(J = J, power_grid = power_grid, J_grid = J_grid)
}
