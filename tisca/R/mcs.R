# ----------------------------------------------------------------------------
# TISCA v2 -- mcs.R
# Bootstrap family inference for comparing more than two models (spec Section 7.2,
# plan Section 2.3): Hansen's SPA, White's Reality Check, and the Model Confidence
# Set (MCS) with T_R and T_max statistics. Romano-Wolf is in multiplicity.R.
#
# All routines operate on a J x M matrix L of per-replication scalar losses
# (LOWER is better): L[j, m] = loss of model m on replication j.
#
# Because replications are i.i.d. here (block bootstrap length 1), the time-series
# block bootstrap of the original MCS collapses to an ordinary resample of rows,
# which is the correct and simpler form for simulation designs (plan 2.3).
#
# Cross-check oracle: the CRAN package `MCS` (Bernardi & Catania). This is the
# primary implementation; the CRAN package is used only in unit tests.
# ----------------------------------------------------------------------------

## ---- common helpers ----------------------------------------------------------

#' Bootstrap resample of loss differentials about their sample mean (null set).
#' Resamples whole rows (i.i.d.), retaining the joint (dependence) structure.
#' @param d J x R differential matrix
#' @param B bootstrap resamples
#' @param seed optional seed
#' @return a J*B x R matrix (or a 3D array structure); for memory-safety, returns
#'   a function that draws a single resample each call instead of B at once.
make_resampler <- function(d, B, seed = NULL) {
  J <- nrow(d)
  idx <- if (is.null(seed)) NULL else {
    set.seed(seed)
  }
  function() {
    i <- sample.int(J, size = J, replace = TRUE)
    d[i, , drop = FALSE]
  }
}

## ---- White's Reality Check ----------------------------------------------------

#' White's Reality Check on a nominated champion against all benchmarks.
#'
#' H0: the champion does not beat every benchmark. Statistics on the advantage
#' vector d_m = L_champ - L_bench_m. T_max = max_m \bar d_m / se(\bar d_m).
#' Bootstrap null is centered (no advantage), so we resample centred differentials
#' about 0. Reject H0 => champion beats ALL benchmarks (with the selection
#' correction built in).
#'
#' @param L J x M loss matrix
#' @param champ index or name of the champion model
#' @param B bootstrap resamples
#' @param seed optional seed
#' @return list(p_value, T_max, indiv_stats, B, rep_est)
reality_check <- function(L, champ, B = 999L, seed = NULL) {
  if (is.null(dim(L))) L <- matrix(L, ncol = 1)
  if (!is.null(seed)) set.seed(seed)
  champ_idx <- if (is.character(champ)) match(champ, colnames(L)) else champ
  champ_vec <- L[, champ_idx]
  jbench <- setdiff(seq_len(ncol(L)), champ_idx)
  if (length(jbench) == 0) {
    return(list(p_value = NA_real_, T_max = 0, indiv_stats = numeric(0),
                B = B, rep_est = numeric(0)))
  }
  # d_m = L_bench_m - L_champ (advantage of the champion; positive = better)
  d <- sweep(L[, jbench, drop = FALSE], 1, champ_vec, FUN = "-")
  J <- nrow(d)
  m_bar <- colMeans(d)
  se <- apply(d, 2, stats::sd) / sqrt(J)
  Tmax_obs <- if (length(m_bar)) max(m_bar / se) else 0
  d_cent <- sweep(d, 2, m_bar)                 # centre on sample mean (null)
  idx <- matrix(sample.int(J, size = J * B, replace = TRUE), nrow = B, ncol = J)
  Tmax_boot <- numeric(B)
  for (b in seq_len(B)) {
    db <- d_cent[idx[b, ], , drop = FALSE]
    mb <- colMeans(db)
    seb <- apply(db, 2, stats::sd) / sqrt(J)
    Tmax_boot[b] <- if (length(mb)) max(mb / seb) else 0
  }
  list(p_value = mean(Tmax_boot >= Tmax_obs),
       T_max = Tmax_obs, indiv_stats = m_bar / se, B = B, rep_est = m_bar)
}

## ---- Hansen's SPA --------------------------------------------------------------

#' Hansen's SPA test that the proposed model beats all benchmarks.
#' Uses the centred bootstrap for the null and Hansen's sample-variance estimate.
#' Kept dependency-free so it works identically in R and Python; the CRAN `MCS`
#' package is only a unit-test oracle.
#'
#' @param L J x M loss matrix, or a J x k advantage matrix when `precomputed_d = TRUE`
#' @param champ index or name of the champion model (required when precomputed_d = FALSE)
#' @param B bootstrap resamples
#' @param seed optional seed
#' @param precomputed_d if TRUE, `L` is already a J x (M-1) advantage matrix
#' @return list(p_value, T_max, B)
spa_test <- function(L, champ = NULL, B = 999L, seed = NULL, precomputed_d = FALSE) {
  if (is.null(dim(L))) L <- matrix(L, ncol = 1)
  if (!is.null(seed)) set.seed(seed)
  if (!precomputed_d) {
    if (is.null(champ)) stop("spa_test requires `champ` when precomputed_d = FALSE")
    champ_idx <- if (is.character(champ)) match(champ, colnames(L)) else champ
    champ_vec <- L[, champ_idx]
    jbench <- setdiff(seq_len(ncol(L)), champ_idx)
    d <- sweep(L[, jbench, drop = FALSE], 1, champ_vec, FUN = "-")
  } else {
    d <- L[, , drop = FALSE]
  }
  J <- nrow(d)
  if (J < 2) stop("spa_test needs >= 2 replications")
  m_bar <- colMeans(d)
  var_h <- colMeans((d - rep(m_bar, each = J))^2)  # Hansen's variance estimator
  se <- sqrt(var_h / J)
  Tmax_obs <- max(m_bar / se)
  d_cent <- d - rep(m_bar, each = J)               # centre on sample mean (null)
  idx <- matrix(sample.int(J, size = J * B, replace = TRUE), nrow = B, ncol = J)
  Tb <- numeric(B)
  for (b in seq_len(B)) {
    db <- d_cent[idx[b, ], , drop = FALSE]
    mb <- colMeans(db)
    vh <- colMeans((db - rep(mb, each = J))^2)
    Tb[b] <- max(mb / sqrt(vh / J))
  }
  list(p_value = mean(Tb >= Tmax_obs), T_max = Tmax_obs, B = B)
}

## ---- Model Confidence Set ------------------------------------------------------

#' Model Confidence Set (Hansen, Lunde & Nason 2011), T_R and T_max variants.
#'
#' Sequentially eliminates the worst model (largest mean loss) while a statistic
#' rejects the hypothesis that all remaining models are equal in expected loss,
#' controlling FWER over the elimination path. Returns the set of models
#' indistinguishable from the best at level alpha and the elimination path.
#'
#' Statistics on the loss matrix:
#'   T_{jk} = \bar d_{.jk} / sqrt(\hat var(d_{.jk})/J),  d_{ijk} = L_{ij} - L_{ik}
#'   T_R (range) = max_{j,k} |T_{jk}|;  T_max = max over j of the max mean.
#' Elimination by the largest mean loss.
#'
#' @param L J x M loss matrix (lower is better)
#' @param B bootstrap resamples (B = 4999 for the final analysis MCSE)
#' @param alpha confidence level (set is 1-alpha "confidence"); complement size alpha
#' @param statistic "T_R" (default) or "T_max"
#' @param seed optional seed
#' @return object of class "tisca_mcs" with models_kept, eliminated, p_values, B,
#'   statistic, mcse
mcs <- function(L, B = 999L, alpha = 0.05,
                statistic = c("T_R","T_max"), seed = NULL) {
  statistic <- match.arg(statistic)
  if (is.null(dim(L))) L <- matrix(L, ncol = 1)
  colnames(L) <- if (is.null(colnames(L))) paste0("M", seq_len(ncol(L))) else colnames(L)
  if (is.null(seed)) seed <- 1L
  J <- nrow(L)
  models <- seq_len(ncol(L))
  eliminated <- character(0)
  p_elim <- numeric(0)
  while (length(models) > 1) {
    subL <- L[, models, drop = FALSE]
    stat_obs <- mcs_stat(subL, statistic)
    stat_boot <- mcs_stat_boot(subL, statistic, B, seed)
    p <- mean(stat_boot >= stat_obs)
    p_elim <- c(p_elim, p)
    if (p <= alpha) {
      worst <- which.max(colMeans(subL))
      eliminated <- c(eliminated, colnames(subL)[worst])
      models <- models[-worst]
      models <- models[order(models)]
    } else {
      break
    }
  }
  mcse <- if (length(p_elim)) sqrt(p_elim[length(p_elim)] * (1 - p_elim[length(p_elim)]) / B) else NA_real_
  structure(list(
    models_kept = colnames(L)[models],
    eliminated  = eliminated,
    p_values    = p_elim,
    B           = B,
    statistic   = statistic,
    mcse        = mcse
  ), class = "tisca_mcs")
}

print.tisca_mcs <- function(x, ...) {
  cat("Model Confidence Set (", x$statistic, "), B = ", x$B, "\n", sep = "")
  cat("kept (indistinguishable from best): ", paste(x$models_kept, collapse = ", "), "\n", sep = "")
  if (length(x$eliminated)) {
    cat("eliminated (path): ", paste(paste0(x$eliminated, " (p=", round(x$p_values, 3), ")"),
                                       collapse = ", "), "\n", sep = "")
  }
  cat("bootstrap-p MCSE: ", ifelse(is.na(x$mcse), "NA", round(x$mcse, 4)), "\n")
  invisible(x)
}

#' Range / max studentized statistic on a loss matrix.
#' @param subL J x k loss matrix
#' @param statistic "T_R" or "T_max"
#' @return scalar statistic
mcs_stat <- function(subL, statistic) {
  k <- ncol(subL); J <- nrow(subL)
  if (k == 1) return(0)
  if (statistic == "T_R") {
    # T_{jk} over all ordered pairs
    tjk <- matrix(0, k, k)
    for (p in seq_len(k)) for (q in seq_len(k)) if (p != q) {
      d <- subL[, p] - subL[, q]
      dbar <- mean(d)
      se <- stats::sd(d) / sqrt(J)
      tjk[p, q] <- if (se > 0) dbar / se else 0
    }
    max(abs(tjk))
  } else { # T_max
    # max of the model-mean studentized statistics (about the grand mean)
    means <- colMeans(subL)
    d_cent <- sweep(subL, 1, rowMeans(subL))
    var_m <- colMeans(d_cent^2)
    se <- sqrt(var_m / J)
    max(means / se)
  }
}

#' Bootstrap null distribution of the MCS statistic (i.i.d. row resample).
#' @param subL J x k loss matrix
#' @param statistic "T_R" or "T_max"
#' @param B resamples
#' @param seed optional seed
#' @return numeric vector of length B of bootstrap statistics
mcs_stat_boot <- function(subL, statistic, B, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  J <- nrow(subL); k <- ncol(subL)
  keep <- which(!is.na(colSums(subL)))
  subL <- subL[, keep, drop = FALSE]
  k <- ncol(subL)
  stats <- numeric(B)
  if (k == 1) return(stats)
  means <- colMeans(subL)
  d0 <- subL - rep(means, each = J)   # centre on the null (all means equal)
  for (b in seq_len(B)) {
    i <- sample.int(J, size = J, replace = TRUE)
    db <- d0[i, , drop = FALSE]
    cb <- colMeans(db)
    # centre bootstrap losses on the bootstrap grand mean (enforce null)
    dbc <- db - rep(mean(cb), each = J)
    if (statistic == "T_R") {
      tjk <- matrix(0, k, k)
      for (p in seq_len(k)) for (q in seq_len(k)) if (p != q) {
        dd <- dbc[, p] - dbc[, q]
        se <- stats::sd(dd) / sqrt(J)
        tjk[p, q] <- if (se > 0) mean(dd) / se else 0
      }
      stats[b] <- max(abs(tjk))
    } else {
      var_m <- colMeans(dbc^2)
      se <- sqrt(var_m / J)
      stats[b] <- max(sapply(seq_len(k), function(m) if (se[m] > 0) cb[m] / se[m] else 0))
    }
  }
  stats
}

## ---- MCS with CRAN oracle (unit-test path only) --------------------------------

#' Cross-check MCS implementation against the CRAN `MCS` package (oracle).
#' @param L J x M loss matrix
#' @param B resamples
#' @param alpha level
#' @param statistic one of "T_R","T_max" as accepted by MCS::MCSprocedure
#' @return the CRAN `MCS` result object
mcs_cran_oracle <- function(L, B = 1000L, alpha = 0.05,
                            statistic = c("T_R","T_max")) {
  if (!requireNamespace("MCS", quietly = TRUE)) {
    stop("CRAN package `MCS` is not installed; cannot run the oracle.")
  }
  statistic <- match.arg(statistic)
  MCS::MCSprocedure(Loss = L, B = B, alpha = alpha, statistic = statistic)
}
