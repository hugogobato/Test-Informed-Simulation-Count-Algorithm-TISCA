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
#' Statistics on the advantage vector d_m = L_bench_m - L_champ (positive =
#' champion better); T_max = max_m \bar d_m / se(\bar d_m). The bootstrap null
#' is centred on the sample mean.
#'
#' READ THE NULL CAREFULLY. H0 is max_m E[L_bench_m - L_champ] <= 0, i.e. the
#' champion is (weakly) the WORST of the set. Rejecting it licenses only
#' "the champion beats AT LEAST ONE benchmark" -- this is a max-type test and it
#' does not, and cannot, establish dominance over every competitor.
#'
#' The positive claim "the champion beats every benchmark" is the
#' intersection-union statement min_m E[L_bench_m - L_champ] > 0, established by
#' rejecting ALL K paired contrasts under FWER control
#' (`romano_wolf_stepdown`), not by this test. Section 3.6 of the manuscript
#' must state it that way.
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

#' Hansen's SPA: the studentized max-type test on paired loss differentials.
#' Same null and same interpretation caveat as `reality_check` above -- it is
#' NOT a test that the champion beats all benchmarks.
#' Uses full recentring (Hansen's liberal SPA_l, equivalently the studentized
#' Reality Check) and Hansen's sample-variance estimate. SPA_l is an upper bound
#' on the SPA_c p-value, so it is conservative for rejection; report the variant.
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

#' Pairwise loss differences and model-average differences.
#'
#' This is the CRAN MCS package's `GetD` calculation.  The diagonal is zero,
#' `mD_ij_bar[i, j] = mean(L_i - L_j)`, and `vD_i_bar` is the average over the
#' other models for each model i.
GetD <- function(mL) {
  mL <- as.matrix(mL)
  iM <- ncol(mL)
  mD_ij_bar <- matrix(0, nrow = iM, ncol = iM)
  if (iM > 1L) {
    for (i in seq_len(iM)) for (j in i:iM) if (i != j) {
      v <- mean(mL[, i] - mL[, j])
      mD_ij_bar[i, j] <- v
      mD_ij_bar[j, i] <- -v
    }
  }
  vD_i_bar <- if (iM > 1L) rowSums(mD_ij_bar) / (iM - 1L) else numeric(iM)
  list(mD_ij_bar = mD_ij_bar, vD_i_bar = vD_i_bar)
}

#' Python/CRAN MCS statistic calculation for one active model set.
#'
#' `bootstrap_indices` is a J x B matrix of 1-based row indices.  The variance
#' estimates are the bootstrap second moments about the observed statistic,
#' exactly as in the Python reference implementation.
.mcs_stats <- function(mL, bootstrap_indices) {
  mL <- as.matrix(mL)
  iT <- nrow(mL); iM <- ncol(mL); B <- ncol(bootstrap_indices)
  d <- GetD(mL)
  mD_ij_bar <- d$mD_ij_bar
  vD_i_bar <- d$vD_i_bar

  mean_res <- matrix(0, nrow = B, ncol = iM)
  for (b in seq_len(B)) mean_res[b, ] <- colMeans(mL[bootstrap_indices[, b], , drop = FALSE])

  aD_ij_bar_res <- array(0, dim = c(B, iM, iM))
  mD_i_bar_res <- matrix(0, nrow = B, ncol = iM)
  for (b in seq_len(B)) {
    aD_ij_bar_res[b, , ] <- outer(mean_res[b, ], mean_res[b, ], "-")
    mD_i_bar_res[b, ] <- rowSums(aD_ij_bar_res[b, , ]) / (iM - 1L)
  }

  mD_ij_bar_var <- matrix(0, nrow = iM, ncol = iM)
  for (i in seq_len(iM)) for (j in seq_len(iM)) {
    mD_ij_bar_var[i, j] <- mean((aD_ij_bar_res[, i, j] - mD_ij_bar[i, j])^2)
  }
  vD_i_bar_var <- vapply(seq_len(iM), function(i)
    mean((mD_i_bar_res[, i] - vD_i_bar[i])^2), numeric(1))

  tij <- matrix(0, nrow = iM, ncol = iM)
  ok <- mD_ij_bar_var > 0
  tij[ok] <- mD_ij_bar[ok] / sqrt(mD_ij_bar_var[ok])
  ti <- numeric(iM)
  ok_i <- vD_i_bar_var > 0
  ti[ok_i] <- vD_i_bar[ok_i] / sqrt(vD_i_bar_var[ok_i])

  tij_res <- array(0, dim = c(iM, iM, B))
  ti_res <- matrix(0, nrow = iM, ncol = B)
  for (b in seq_len(B)) {
    tij_res[, , b][ok] <- (aD_ij_bar_res[b, , ][ok] - mD_ij_bar[ok]) /
      sqrt(mD_ij_bar_var[ok])
    if (any(ok_i)) ti_res[ok_i, b] <-
      (mD_i_bar_res[b, ok_i] - vD_i_bar[ok_i]) / sqrt(vD_i_bar_var[ok_i])
  }
  list(mD_ij_bar = mD_ij_bar, vD_i_bar = vD_i_bar,
       mD_ij_bar_var = mD_ij_bar_var, vD_i_bar_var = vD_i_bar_var,
       tij = tij, ti = ti, tij_res = tij_res, ti_res = ti_res)
}

#' Model Confidence Set (Hansen, Lunde & Nason 2011), port of Python `mcs`.
#'
#' The elimination path always runs down to one model.  `p_H0` records the
#' bootstrap p-value at the step where a model is removed, while `p_mcs` is the
#' running maximum of those p-values.  Hansen's confidence set is selected from
#' `p_mcs`, not from the individual step p-values.
#'
#' @param L J x M loss matrix (lower is better)
#' @param B number of bootstrap resamples
#' @param alpha complement of the confidence level
#' @param statistic "Tmax", "TR", "T_max", or "T_R"
#' @param seed optional seed used once to draw all bootstrap indices
#' @param bootstrap_indices optional J x B matrix of 1-based row indices
#' @return object of class "tisca_mcs"
mcs <- function(L, B = 999L, alpha = 0.05,
                statistic = "T_R", seed = NULL,
                bootstrap_indices = NULL) {
  if (is.null(dim(L))) L <- matrix(L, ncol = 1)
  L <- as.matrix(L)
  if (!is.numeric(L) || any(!is.finite(L))) stop("L must be a finite numeric matrix")
  if (nrow(L) < 1L || ncol(L) < 1L) stop("L must have at least one row and column")
  if (length(B) != 1L || B <= 0) stop("B must be positive")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be in [0, 1]")
  }
  aliases <- c(Tmax = "Tmax", T_max = "Tmax", TR = "TR", T_R = "TR")
  canonical <- unname(aliases[statistic])
  if (length(statistic) != 1L || length(canonical) != 1L || is.na(canonical)) {
    stop("statistic must be one of 'Tmax', 'TR', 'T_max', or 'T_R'")
  }
  statistic <- canonical
  B <- as.integer(B)
  J <- nrow(L); iM <- ncol(L)
  if (is.null(colnames(L))) colnames(L) <- paste0("M", seq_len(iM))
  model_names <- colnames(L)

  if (is.null(bootstrap_indices)) {
    if (!is.null(seed)) set.seed(seed)
    bootstrap_indices <- matrix(sample.int(J, size = J * B, replace = TRUE),
                                nrow = J, ncol = B)
  } else {
    bootstrap_indices <- as.matrix(bootstrap_indices)
    if (!is.numeric(bootstrap_indices) ||
        !identical(dim(bootstrap_indices), c(J, B)) ||
        any(!is.finite(bootstrap_indices)) ||
        any(bootstrap_indices != floor(bootstrap_indices)) ||
        any(bootstrap_indices < 1 | bootstrap_indices > J)) {
      stop("bootstrap_indices must be a J x B matrix of 1-based row indices")
    }
    bootstrap_indices <- matrix(as.integer(bootstrap_indices), nrow = J, ncol = B)
  }

  p_H0 <- rep(NA_real_, iM)
  p_mcs <- rep(NA_real_, iM)
  avg_loss <- colMeans(L)
  elimination_order <- character(0)
  elimination_pvalues <- numeric(0)
  working <- L
  working_names <- model_names

  while (ncol(working) > 1L) {
    st <- .mcs_stats(working, bootstrap_indices)
    if (statistic == "Tmax") {
      obs <- max(st$ti)
      boot <- apply(st$ti_res, 2, max)
      e <- which.max(st$ti)
    } else {
      obs <- max(abs(st$tij))
      boot <- apply(abs(st$tij_res), 3, max)
      e <- which.max(apply(st$tij, 1, max))
    }
    p <- mean(boot > obs)
    eliminated <- working_names[e]
    pos <- match(eliminated, model_names)
    p_H0[pos] <- p
    p_mcs[pos] <- max(p_H0, na.rm = TRUE)
    elimination_order <- c(elimination_order, eliminated)
    elimination_pvalues <- c(elimination_pvalues, p)
    working <- working[, -e, drop = FALSE]
    working_names <- working_names[-e]
  }

  p_H0[is.na(p_H0)] <- 1.0
  p_mcs[is.na(p_mcs)] <- 1.0
  included <- model_names[p_mcs > alpha]
  excluded <- model_names[p_mcs <= alpha]
  mcse <- if (length(elimination_pvalues)) {
    p <- elimination_pvalues[length(elimination_pvalues)]
    sqrt(p * (1 - p) / B)
  } else NA_real_
  tab_order <- order(p_mcs)
  table <- cbind(avg_loss = avg_loss[tab_order], p_H0 = p_H0[tab_order],
                 p_mcs = p_mcs[tab_order])
  rownames(table) <- model_names[tab_order]

  structure(list(
    models_kept = included,
    eliminated = elimination_order,
    p_values = elimination_pvalues,
    included = included,
    excluded = excluded,
    p_H0 = setNames(p_H0, model_names),
    p_mcs = setNames(p_mcs, model_names),
    elimination_order = elimination_order,
    elimination_pvalues = elimination_pvalues,
    avg_loss = setNames(avg_loss, model_names),
    table = table,
    B = B, alpha = alpha, statistic = statistic, seed = seed, mcse = mcse
  ), class = "tisca_mcs")
}

print.tisca_mcs <- function(x, ...) {
  cat("Model Confidence Set (", x$statistic, "), B = ", x$B, "\n", sep = "")
  cat("confidence set: ", paste(x$included, collapse = ", "), "\n", sep = "")
  cat("excluded: ", paste(x$excluded, collapse = ", "), "\n", sep = "")
  if (length(x$elimination_order)) {
    cat("elimination path: ", paste(
      paste0(x$elimination_order, " (p_H0=", round(x$elimination_pvalues, 3), ")"),
      collapse = ", "), "\n", sep = "")
  }
  cat("bootstrap-p MCSE: ", ifelse(is.na(x$mcse), "NA", round(x$mcse, 4)), "\n")
  invisible(x)
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
