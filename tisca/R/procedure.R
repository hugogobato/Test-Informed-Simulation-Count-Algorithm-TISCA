# ----------------------------------------------------------------------------
# TISCA v2 -- procedure.R
# The simulation-design procedure: two-stage default (Algorithm 1) and the
# optional adaptive mode (Algorithm 2), per spec Sections 5-6 (C7, C8).
#
# Algorithm 1 (default): independent-seed pilot (seed block [J0_pilot_start,
# J0_pilot_start + J0 - 1]) estimates sigma_D / Sigma_D_UB; solve closed-form J;
# discard pilot rows; run J confirmatory replications (seed block [1, J]); final
# test once at alpha_adj. Because the pilot is discarded, unconditional Type I
# error is exactly alpha_adj by iterated expectation (spec Section 5).
#
# The caller supplies sim_func(seed) -> data.frame of per-replication metric
# columns (model A and B). This module computes the contrasts, plans J, runs the
# confirmatory phase, and runs the final paired test.
# ----------------------------------------------------------------------------

## ---- Algorithm 1: two-stage default -------------------------------------------

#' Build paired contrast(s) from a per-replication metric data frame and report
#' listwise drop counts (reviewer suggestion 2b: never drop per-column).
#'
#' @param frame data.frame with one row per replication; columns = metrics
#' @param primary list of list(metric_a, metric_b, mode, delta, Delta, target_type,
#'        target_m, target_h, target_power, metric_type, lower_is_better)
#' @return list(contrasts = list of list(D, dropped, n_kept), families)
make_contrasts <- function(frame, primary) {
  out <- vector("list", length(primary))
  for (i in seq_along(primary)) {
    ctr <- primary[[i]]
    A <- frame[[ctr$metric_a]]
    B <- frame[[ctr$metric_b]]
    out[[i]] <- contrast_from_columns(A, B)
  }
  out
}

#' One pass of the design solve for a single contrast under a given inflation.
#' @inheritParams solve_power_J
#' @param ctr a single primary-contrast spec (list with mode, delta, Delta, target_type, ...)
#' @param sd_D pilot sd of the contrast
#' @param alpha_adj adjusted level
#' @return list(J, target, achieved, kind)
solve_contrast_J <- function(ctr, sd_D, alpha_adj) {
  sd_ub <- sigma_ub(sd_D, ctr$J0 %||% 50, gamma = ctr$gamma %||% 0.20)
  kind <- ctr$target_type %||% "power"
  if (kind == "power") {
    if (is.null(ctr$target_power)) stop("power target requires target_power")
    r <- solve_power_J(ctr$mode, ctr$delta, sd_ub, alpha_adj,
                       target_power = ctr$target_power,
                       Delta = ctr$Delta, J_max = ctr$J_max %||% 1e6)
    list(J = r$J, kind = "power", achieved = r$achieved_power,
         capped = isTRUE(r$capped), degenerate = isTRUE(r$degenerate),
         infeasible = isTRUE(r$infeasible), message = r$message %||% "")
  } else if (kind == "mcse") {
    r <- solve_mcse_J(sd_ub, ctr$target_m, J_max = ctr$J_max %||% 1e6)
    list(J = r$J, kind = "mcse", achieved = r$achieved,
         capped = isTRUE(r$capped), degenerate = isTRUE(r$degenerate),
         infeasible = FALSE, message = "")
  } else if (kind == "halfwidth") {
    r <- solve_halfwidth_J(sd_ub, ctr$target_h, alpha = ctr$alpha %||% alpha_adj,
                           J_max = ctr$J_max %||% 1e6)
    list(J = r$J, kind = "halfwidth", achieved = r$achieved,
         capped = isTRUE(r$capped), degenerate = isTRUE(r$degenerate),
         infeasible = FALSE, message = "")
  } else {
    stop("unknown target_type: ", kind)
  }
}

#' TISCA v2 -- two-stage default (Algorithm 1).
#'
#' @param sim_func function(seed) -> data.frame of per-replication columns
#' @param primary list of primary-contrast specs (fields: metric_a, metric_b, mode,
#'        delta, Delta, target_type, target_m/target_h/target_power, J0, gamma, J_max,
#'        metric_type, lower_is_better)
#' @param J0 pilot size (>= 2)
#' @param pilot_seed_start first pilot seed (default 1000001 per seed protocol)
#' @param confirm_seed_start first confirmatory seed (default 1)
#' @param J_max hard budget cap
#' @param alpha family level
#' @param correction multiplicity method ("none","bonferroni","holm","BH","romanowolf")
#' @param gamma assurance level (0.20 default)
#' @param success_criterion "conjunctive","disjunctive","named"
#' @param rm_nas logical, always TRUE internally (listwise); kept for API parity
#' @param verbose print progress
#' @param B bootstrap resamples for the optional final inference
#' @param seed optional master RNG seed for reproducibility of the bootstrap only
#' @return list with the full audit trail: pilot_sizes, J_detail, J_final, capped,
#'         per-contrast final inference, and the pilot discard note.
run_two_stage <- function(sim_func, primary, J0 = 50L,
                          pilot_seed_start = 1000001L,
                          confirm_seed_start = 1L,
                          J_max = 1e6L, alpha = 0.05,
                          correction = "bonferroni",
                          gamma = 0.20, success_criterion = "conjunctive",
                          rm_nas = TRUE, verbose = TRUE, B = 4999L,
                          seed = NULL) {
  stopifnot(is.function(sim_func))
  K <- length(primary)

  # ----- validation ------------------------------------------------
  vfamily <- validate_family(primary, K_expected = K)
  if (!vfamily$pass) stop(paste("TISCA v2 family validation failed:\n",
                                paste(vfamily$problems, collapse = "\n")))

  # ----- multiplicity-adjusted level -------------------------------
  adj <- adj_alpha(correction, alpha, K)
  alpha_adj <- adj$alpha_adj

  # ----- Stage 1: independent-seed pilot ---------------------------
  pilot_rows <- vector("list", J0)
  for (j in seq_len(J0)) pilot_rows[[j]] <- sim_func(pilot_seed_start + j - 1L)
  pilot_frame <- do.call(rbind, pilot_rows)

  # build contrasts (listwise) and estimate pilot s_D + sigma_D_UB
  pilot_contrasts <- make_contrasts(pilot_frame, primary)
  p_sd <- sapply(pilot_contrasts, function(cc) stats::sd(cc$D))
  p_dbar <- sapply(pilot_contrasts, function(cc) mean(cc$D))
  p_dropped <- sapply(pilot_contrasts, function(cc) cc$dropped)

  # ----- Stage 2: solve J per contrast, combine ---------------------
  details <- lapply(seq_len(K), function(i) {
    solve_contrast_J(primary[[i]], p_sd[i], alpha_adj)
  })
  J_detail <- sapply(details, function(d) d$J)
  combined <- combine_J(J_detail, J_max)
  J_final <- combined$J_final

  if (verbose) {
    cat("TISCA v2 two-stage | pilot J0 =", J0, "| alpha_adj =", round(alpha_adj, 5), "\n")
    cat("pilot sd_D:", paste(round(p_sd, 4), collapse = ", "), "\n")
    cat("J_detail:", paste(J_detail, collapse = ", "), "-> J_final =", J_final,
        if (combined$capped) "(CAPPED at J_max)" else "", "\n")
  }

  # ----- Confirmatory run (seeds [confirm_seed_start, .. + J_final)) -------
  confirm_rows <- vector("list", J_final)
  for (j in seq_len(J_final)) confirm_rows[[j]] <- sim_func(confirm_seed_start + j - 1L)
  confirm_frame <- do.call(rbind, confirm_rows)

  # ----- Final inference (ONCE at alpha_adj), estimate-first ---------------
  final <- vector("list", K)
  for (i in seq_len(K)) {
    ctr <- primary[[i]]
    cc <- contrast_from_columns(confirm_frame[[ctr$metric_a]],
                                confirm_frame[[ctr$metric_b]])
    D <- cc$D
    mode <- ctr$mode %||% "M1"
    alternative <- if (mode %in% c("M2","M3","M4")) "less" else
                   if (mode == "M5") "two.sided" else "two.sided"
    bin_type <- identical(ctr$metric_type %||% "loss", "binary")
    if (bin_type) {
      res <- mcnemar_exact(confirm_frame[[ctr$metric_a]],
                           confirm_frame[[ctr$metric_b]],
                           alpha = alpha_adj, alternative = alternative)
      res$estimate <- res$estimate  # share of discordant favouring A
    } else {
      res <- paired_t(D, alpha = alpha_adj, alternative = alternative)
      # recommend studentized bootstrap when skew is non-trivial
      g1 <- sum(mean((D - mean(D))^3) / stats::sd(D)^3 + 0)  # skewness
      res$skew <- g1
    }
    res$n_kept <- cc$n_kept
    res$dropped <- cc$dropped
    final[[i]] <- res
  }

  list(
    J_final       = J_final,
    capped        = combined$capped,
    J_detail      = J_detail,
    alpha_adj     = alpha_adj,
    correction    = correction,
    gamma         = gamma,
    p_sd_pilot    = p_sd,
    p_dbar_pilot  = p_dbar,
    p_dropped     = p_dropped,
    pilot_frame   = pilot_frame,
    confirm_frame = confirm_frame,
    final         = final,
    pilot_discarded = TRUE,
    stage_description = "two-stage (Algorithm 1); pilot discarded from final inference"
  )
}

## ---- Algorithm 2: adaptive mode (optional) ------------------------------------

#' TISCA v2 -- adaptive / internal-pilot re-estimation (Algorithm 2, optional).
#'
#' Reuses the pilot rows in the final test (that is what distinguishes it from
#' Algorithm 1 and is the source of its error-rate risk). Kieser & Friede (2000)
#' and Proschan (2005) are the cited references for the non-nominal Type I error.
#' The algorithm is OFF by default; the caller must explicitly request it and
#' accept the validation caveat.
#'
#' @inheritParams run_two_stage
#' @param nmax_looks max re-estimation look times
#' @param look_schedule explicit vector of cumulative J at which to re-estimate
#'        (overrides nmax_looks for timings)
#' @param accept_caveat must be TRUE to confirm the caller accepts the risk
#' @return list like run_two_stage, flagged as adaptive with pilot rows REUSED.
run_adaptive <- function(sim_func, primary, J0 = 50L,
                         pilot_seed_start = 1000001L,
                         confirm_seed_start = 1L,
                         J_max = 1e6L, alpha = 0.05,
                         correction = "bonferroni", gamma = 0.20,
                         success_criterion = "conjunctive",
                         nmax_looks = 3L, look_schedule = NULL,
                         accept_caveat = FALSE, verbose = TRUE,
                         B = 4999L, seed = NULL) {
  if (!isTRUE(accept_caveat)) {
    stop("Algorithm 2 (adaptive) reuses pilot rows in the final test and its ",
         "unconditional error rates must be validated (P3-T2/P3-T3) before use. ",
         "Set accept_caveat = TRUE to proceed; the two-stage Algorithm 1 is the default.")
  }
  K <- length(primary)
  adj <- adj_alpha(correction, alpha, K)
  alpha_adj <- adj$alpha_adj

  # run pilot as "block 1" replications; adaptively accumulate
  block_size <- J0
  acc_rows <- vector("list", block_size)
  for (j in seq_len(block_size)) acc_rows[[j]] <- sim_func(pilot_seed_start + j - 1L)
  frame <- do.call(rbind, acc_rows)
  J <- block_size

  if (is.null(look_schedule)) {
    look_schedule <- unique(sort(block_size * 2^(seq_len(nmax_looks) - 1)))
  }
  progressed <- TRUE
  while (progressed && J <= J_max) {
    contrast <- make_contrasts(frame, primary)
    sds  <- sapply(contrast, function(cc) stats::sd(cc$D))
    dets <- lapply(seq_len(K), function(i) solve_contrast_J(primary[[i]], sds[i], alpha_adj))
    req_J <- max(sapply(dets, function(d) if (is.na(d$J)) J_max else d$J))
    req_J <- min(req_J, J_max)
    next_look <- min(look_schedule[look_schedule > J], req_J)
    if (length(next_look) == 0 || is.na(next_look)) next_look <- req_J
    if (req_J <= J) {
      progressed <- FALSE
    } else if (next_look > J_max) {
      progressed <- FALSE
    } else {
      add <- max(1L, min(next_look - J, J_max - J))
      for (jj in seq_len(add)) {
        idx <- length(acc_rows) + jj
        acc_rows[[idx]] <- sim_func(confirm_seed_start + idx - 1L)
      }
      frame <- do.call(rbind, acc_rows)
      J <- length(acc_rows)
    }
  }

  J_final <- min(J, J_max)
  capped <- J_final >= J_max

  # final test ONCE at alpha_adj on the FULL accumulated set INCLUDING pilot
  final <- vector("list", K)
  for (i in seq_len(K)) {
    ctr <- primary[[i]]
    cc <- contrast_from_columns(frame[[ctr$metric_a]], frame[[ctr$metric_b]])
    alternative <- if (ctr$mode %in% c("M2","M3","M4")) "less" else "two.sided"
    res <- paired_t(cc$D, alpha = alpha_adj, alternative = alternative)
    res$n_kept <- cc$n_kept; res$dropped <- cc$dropped
    final[[i]] <- res
  }

  list(
    J_final = J_final, capped = capped, alpha_adj = alpha_adj,
    correction = correction, gamma = gamma,
    pilot_frame = frame,
    final = final,
    pilot_discarded = FALSE,
    stage_description = paste0("adaptive (Algorithm 2); pilot rows REUSED in the ",
                               "final test; error rates are NOT nominal by default ",
                               "(Kieser-Friede 2000, Proschan 2005)")
  )
}

# `%||%` is defined once in tisca.R and available package-wide.
