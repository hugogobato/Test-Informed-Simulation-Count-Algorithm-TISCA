#!/usr/bin/env Rscript
#
# run_cell.R -- one shard of the MVBCF case-study re-run (plan P3-T5[a]).
#
# Adapted from the original authors' simulation driver, GitHub_DGP1.R, located
# upstream at https://github.com/Nathan-McJames/MVBCF_Paper (the authoritative
# spec of the DGP, the model settings, and the recorded metrics). That file is
# BROUGHT IN and not copied into this repository; the TISCA repo commits only
# this derived driver. See §0.1 of REVISION_PLAN.md.
#
# Intentional deviations from the upstream script (REVISION_PLAN.md §1.3, §1.7,
# P3-T5[a], and docs/seed_rng_protocol.md):
#   1. Fixed integer seeds passed as CLI args; as.numeric(Sys.time()) is
#      deleted. n is a CLI argument, not sample(c(100,500,1000),1).
#   2. One BCF benchmark, calibrated: stochtree::bcf with num_gfr = 50 AND
#      McJames's priors (50 prognostic / 20 treatment trees, outcome-specific
#      sigma2_leaf_init). The v1 cold-start and warm-start-default-prior
#      variants are both deleted. Gated by P3-T5(e).
#   3. Separate L'Ecuyer-CMRG RNG streams for data generation and model fits,
#      interleaved as substreams (2j, 2j+1) from one cell master (§3.3).
#   4. MVBCF fitted via fast_bart() from the upstream MVBCF_Code.cpp through
#      sourceCpp (as the original), not mvbcf::run_mvbcf() (§1.7c, P3-T5[a].5).
#   5. The full original metric set is kept and four superset columns added:
#      cate_mse, ate_sq_err, fit_seconds, converged_flag/error_message (§4.5).
#   6. Stochtree/dbarts fits run single-threaded (nthread = 1) so mc.cores is
#      the only source of parallelism (reproducibility, §3.5).
#
# CLI:
#   Rscript run_cell.R <dgp> <n> <seed_start> <seed_end> --out <csv> \
#       [--cores <k>] [--mode confirmatory|pilot]
#
# Each completed replication is appended to --out under mutex immediately (one
# row per seed), so a killed session loses at most one replication. Row order is
# not guaranteed under mclapply; the seed column is authoritative.
suppressPackageStartupMessages({
  library(parallel)
  library(mvtnorm)
  library(bartCause)
  library(skewBART)
  library(dbarts)
  library(stochtree)
  library(scoringRules)
  library(matrixStats)
  library(digest)
  library(Rcpp)
})

# Compile the upstream MVBCF C++ INTO THIS PROCESS. fast_bart() is a plain
# Rcpp[[export]] from MVBCF_Code.cpp, which is BROUGHT IN (never copied to this
# repo) and must be present when this script runs. The Colab notebook sets
# TISCA_MVBCF_CPP to the on-session copy it downloaded from the upstream raw URL.
cpp_path <- Sys.getenv("TISCA_MVBCF_CPP", NA_character_)
if (!is.na(cpp_path) && nzchar(cpp_path) && file.exists(cpp_path)) {
  Rcpp::sourceCpp(cpp_path)
} else if (!exists("fast_bart", mode = "function")) {
  stop("MVBCF_Code.cpp not found. Set TISCA_MVBCF_CPP to the (downloaded) path to ",
       "the upstream file, or ensure fast_bart() is already compiled into the session.")
}

# -----------------------------------------------------------------------------
# Argument parsing (base R only)
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
dgp_arg <- as.integer(args[[1]])
n_arg   <- as.integer(args[[2]])
s_start <- as.integer(args[[3]])
s_end   <- as.integer(args[[4]])
a <- args[-(1:4)]
get_opt <- function(key, default = NULL) {
  for (i in seq_along(a)) if (a[i] == key) return(a[i + 1L])
  default
}
out_path <- get_opt("--out", stop("--out required"))
cores    <- as.integer(get_opt("--cores", "1"))
mode     <- get_opt("--mode", "confirmatory")

stopifnot(!is.na(dgp_arg) && dgp_arg %in% 1:3)
stopifnot(!is.na(n_arg) && n_arg %in% c(100, 500))
stopifnot(!is.na(s_start), !is.na(s_end), s_end >= s_start, s_start >= 0)
stopifnot(cores >= 1L)
mkdir_out <- dirname(out_path)
if (!dir.exists(mkdir_out)) dir.create(mkdir_out, recursive = TRUE)

# Cell master and seed accounting (seed protocol §5): a cell uses the SAME
# L'Ecuyer master for both phases (confirmatory master 1, pilot master 2) and
# the disjoint index blocks do the work. The `seed` column carries the
# human-auditable labels: confirmatory 0..J-1, pilot 1_000_001..1_000_000+J0,
# while stream construction uses the small per-master indices. as.numeric(
# Sys.time()) never appears in this driver.
cell_master <- if (mode == "pilot") 2L else 1L

# -----------------------------------------------------------------------------
# RNG protocol: L'Ecuyer-CMRG interleaved data/fit substreams. Never derive a
# stream by arithmetic on .Random.seed (§3.1).
# -----------------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")

make_streams <- function(master, j_from, j_to) {
  set.seed(master)
  s <- .Random.seed
  n <- j_to - j_from + 1L
  data_s <- fit_s <- vector("list", n)
  idx <- 0L
  for (j in 0:j_to) {
    d <- s; s <- nextRNGStream(s)      # substream 2j   -> data
    f <- s; s <- nextRNGStream(s)      # substream 2j+1 -> fit
    if (j >= j_from) { idx <- idx + 1L; data_s[[idx]] <- d; fit_s[[idx]] <- f }
  }
  list(data = data_s, fit = fit_s)
}

set_rand <- function(x) assign(".Random.seed", x, envir = .GlobalEnv)

# Per-model substream of a replication's fit stream. Resetting the RNG to the
# SAME state before every model (which is what an unqualified
# `set_rand(f_stream)` does) makes each model draw the identical seed and share
# one state, so `model_seed_*` is constant across models and the two BCF
# outcome fits run on the same stream. Advancing the fit stream by `i`
# L'Ecuyer substreams gives every model its own independent, reproducible
# stream while keeping each fit independent of what ran before it.
model_stream <- function(f_stream, i) {
  s <- f_stream
  if (i > 0L) for (k in seq_len(i)) s <- nextRNGStream(s)
  s
}
# Fixed substream index per model, so a re-run of any single model is exact.
MODEL_SUBSTREAM <- c(propensity = 0L, mvbcf = 1L, bcf1 = 2L, bcf2 = 3L,
                     bart1 = 4L, bart2 = 5L, mvbart = 6L)
use_model_stream <- function(f_stream, name) {
  set_rand(model_stream(f_stream, MODEL_SUBSTREAM[[name]]))
  sample.int(.Machine$integer.max, 1L)
}

# -----------------------------------------------------------------------------
# DGP generation. Formulas VERBATIM from GitHub_DGP1.R (44-98), GitHub_DGP2.R
# (43-98), GitHub_DGP3.R (43-98).
# -----------------------------------------------------------------------------
generate_data <- function(dgp, n) {
  X1 <- runif(n); X2 <- runif(n); X3 <- runif(n); X4 <- runif(n); X5 <- runif(n)
  X6 <- rbinom(n, 1, 0.5); X7 <- rbinom(n, 1, 0.5); X8 <- rbinom(n, 1, 0.5)
  X9 <- sample(c(0, 1, 2, 3, 4), n, replace = TRUE)
  X10 <- sample(c(0, 1, 2, 3, 4), n, replace = TRUE)
  if (dgp == 1L) {
    X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
    Mu1 <- (11*sin(pi*X1*X2)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
    Mu2 <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+14*X4+8*X6+X9)*10+300
    Tau1 <- (2*X4+2*X5)*10
    Tau2 <- (1*X4+3*X5)*10
  } else if (dgp == 2L) {
    X <- cbind(X1, X2, X3, X5, X6, X7, X8, X9, X10)
    Mu1 <- (11*sin(pi*X1*X2)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
    Mu2 <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+0*X4+8*X6+X9)*10+300
    Tau1 <- (2*X4+2*X5)*10
    Tau2 <- (1*X4+3*X5)*10
  } else {
    X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
    Mu1 <- (11*sin(pi*X4*X5)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
    Mu2 <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+14*X4+8*X6+X9)*10+300
    Tau1 <- (2*X4+2*X2)*10
    Tau2 <- (1*X3+3*X5)*10
  }
  true_propensity <- X4
  Z <- rbinom(n, 1, true_propensity)
  Y <- cbind(Mu1 + Z*Tau1, Mu2 + Z*Tau2) +
    mvtnorm::rmvnorm(n, c(0, 0), matrix(c(50^2, 0, 0, 50^2), nrow = 2, byrow = TRUE))
  list(X = X, Y = Y, Z = Z, Tau1 = Tau1, Tau2 = Tau2)
}

generate_test <- function(dgp, n_test = 1000L) {
  X1 <- runif(n_test); X2 <- runif(n_test); X3 <- runif(n_test); X4 <- runif(n_test); X5 <- runif(n_test)
  X6 <- rbinom(n_test, 1, 0.5); X7 <- rbinom(n_test, 1, 0.5); X8 <- rbinom(n_test, 1, 0.5)
  X9 <- sample(c(0, 1, 2, 3, 4), n_test, replace = TRUE)
  X10 <- sample(c(0, 1, 2, 3, 4), n_test, replace = TRUE)
  if (dgp == 1L) {
    X_test <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
    Mu1_test <- (11*sin(pi*X1*X2)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
    Mu2_test <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+14*X4+8*X6+X9)*10+300
    Tau1_test <- (2*X4+2*X5)*10
    Tau2_test <- (1*X4+3*X5)*10
  } else if (dgp == 2L) {
    X_test <- cbind(X1, X2, X3, X5, X6, X7, X8, X9, X10)
    Mu1_test <- (11*sin(pi*X1*X2)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
    Mu2_test <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+0*X4+8*X6+X9)*10+300
    Tau1_test <- (2*X4+2*X5)*10
    Tau2_test <- (1*X4+3*X5)*10
  } else {
    X_test <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
    Mu1_test <- (11*sin(pi*X4*X5)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
    Mu2_test <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+14*X4+8*X6+X9)*10+300
    Tau1_test <- (2*X4+2*X2)*10
    Tau2_test <- (1*X3+3*X5)*10
  }
  true_propensity_test <- X4
  Z_test <- rbinom(n_test, 1, true_propensity_test)
  Y_test <- cbind(Mu1_test + Z_test*Tau1_test, Mu2_test + Z_test*Tau2_test) +
    mvtnorm::rmvnorm(n_test, c(0, 0), matrix(c(50^2, 0, 0, 50^2), nrow = 2, byrow = TRUE))
  list(X_test = X_test, Y_test = Y_test, Z_test = Z_test,
       Tau1_test = Tau1_test, Tau2_test = Tau2_test)
}

# -----------------------------------------------------------------------------
# Metric helpers (identical to upstream)
# -----------------------------------------------------------------------------
in_cred <- function(samples, value, interval) {
  q1 <- quantile(samples, 0 + (1 - interval)/2, names = FALSE)
  q2 <- quantile(samples, 1 - (1 - interval)/2, names = FALSE)
  ifelse(value >= q1 & value <= q2, TRUE, FALSE)
}
cred_width <- function(samples, interval) {
  diff(quantile(samples, c(0 + (1 - interval)/2, 1 - (1 - interval)/2), names = FALSE))
}

# -----------------------------------------------------------------------------
# Model constants (fixed to match the original / plan §1.7)
# -----------------------------------------------------------------------------
n_test <- 1000L
n_tree_mu <- 50L; n_tree_tau <- 20L; n_trees_total <- n_tree_mu + n_tree_tau
n_iter <- 1000L; n_burn <- 500L; n_mcmc <- n_iter - n_burn
mu_val <- 1; tau_val <- 0.375; v_val <- 1; wish_val <- 1; min_val <- 1
nthread_global <- 1L

# -----------------------------------------------------------------------------
# Append a data.frame to a CSV under a file lock (only the parent collects
# results, so concurrency is not between workers but is kept safe anyway).
# -----------------------------------------------------------------------------
append_csv <- function(df, path) {
  write.table(df, path, row.names = FALSE, sep = ",",
              col.names = !file.exists(path), append = file.exists(path))
}
append_csv_locked <- function(df, path) {
  lk <- paste0(path, ".lock")
  wait <- 0
  while (file.exists(lk) && wait < 300) { Sys.sleep(0.2); wait <- wait + 1 }
  if (!file.exists(lk)) file.create(lk)
  on.exit(unlink(lk), add = TRUE)
  append_csv(df, path)
}

# -----------------------------------------------------------------------------
# Build one output row from a completed replication.
# -----------------------------------------------------------------------------
make_stub <- function() {
  # A named LIST (not `c(...)`): c() would coerce everything up to character
  # because of the NA_character_ fields, silently stringifying every numeric
  # column. A list preserves integer/double/character type per field.
  stub <- list(
    seed = NA_integer_, n = NA_integer_, dgp = NA_integer_,
    seed_cell_master = NA_integer_, rng_kind = "L'Ecuyer-CMRG",
    seed_data_hash = NA_character_, seed_fit_hash = NA_character_,
    model_seed_propensity = NA_integer_, model_seed_mvbcf = NA_integer_,
    model_seed_bcf1 = NA_integer_, model_seed_bcf2 = NA_integer_,
    model_seed_bart1 = NA_integer_, model_seed_bart2 = NA_integer_,
    model_seed_mvbart = NA_integer_,
    fit_seconds_mvbcf = NA_real_, fit_seconds_bcf1 = NA_real_,
    fit_seconds_bcf2 = NA_real_, fit_seconds_bart1 = NA_real_,
    fit_seconds_bart2 = NA_real_, fit_seconds_mvbart = NA_real_,
    converged_flag = NA_integer_, error_message = NA_character_,
    hostname = NA_character_, git_sha = NA_character_,
    session_hash = NA_character_, seq_phase = NA_character_,
    replication_seconds = NA_real_
  )
  models <- c("mvbcf", "bcf", "bart", "mvbart")
  meas <- c("pehe1","pehe2","bias1","bias2","ate1","ate2","crmse1","crmse2",
            "atesq1","atesq2","crps1","crps2",
            "cov501","cov951","cov502","cov952","wid501","wid951","wid502","wid952")
  for (m in models) for (mm in meas) stub[[paste0(m, "_", mm)]] <- NA_real_
  ate <- c("ate_l951","ate_u951","ate_l952","ate_u952","ate_l501","ate_u501",
           "ate_l502","ate_u502","ate_wid951","ate_wid952","ate_wid501","ate_wid502",
           "ate_cov951","ate_cov952","ate_cov501","ate_cov502")
  for (m in models) for (mm in ate) stub[[paste0(m, "_", mm)]] <- NA_real_
  stub
}
STUB <- make_stub()
COL_NAMES <- names(STUB)

# -----------------------------------------------------------------------------
# One replication.
#   d_stream / f_stream : L'Ecuyer substreams (already indexed to this seed)
#   out_path            : shard CSV to append this row to
# -----------------------------------------------------------------------------
run_one <- function(seed_emit, d_stream, f_stream, out_path) {
  replication_start <- proc.time()[["elapsed"]]
  row <- STUB
  row[["seed"]] <- seed_emit
  row[["n"]] <- n_arg
  row[["dgp"]] <- dgp_arg
  row[["seed_cell_master"]] <- cell_master
  row[["seq_phase"]] <- mode
  row[["converged_flag"]] <- 1L
  err_msg <- NA_character_

  # ---- data generation on its own stream ----
  set_rand(d_stream)
  row[["seed_data_hash"]] <- paste(.Random.seed, collapse = "_")
  d <- generate_data(dgp_arg, n_arg)
  t <- generate_test(dgp_arg, n_test)
  X <- d$X; Y <- d$Y; Z <- d$Z
  Tau1 <- d$Tau1; Tau2 <- d$Tau2
  X_test <- t$X_test; Z_test <- t$Z_test
  Tau1_test <- t$Tau1_test; Tau2_test <- t$Tau2_test
  X2 <- X; X2_test <- X_test

  # ---- propensity score (from here on, the fit stream) ----
  set_rand(f_stream)
  row[["seed_fit_hash"]] <- paste(.Random.seed, collapse = "_")
  f_seed <- use_model_stream(f_stream, "propensity")
  row[["model_seed_propensity"]] <- f_seed
  # GitHub_DGP1.R:101 calls bart(x.train, y.train, x.test, k = 3) with the
  # dbarts DEFAULT sampler settings (ndpost = 1000, nskip = 100). The
  # propensity feeds every downstream model, so the calibration gate against
  # the published Table 2 is only meaningful if it is estimated the same way:
  # do not substitute the 500/500 schedule used for the outcome models.
  p_mod <- tryCatch(
    dbarts::bart(x.train = X, y.train = Z, x.test = X_test, k = 3,
                 n.threads = nthread_global, seed = f_seed),
    error = function(e) { err_msg <<- conditionMessage(e); NULL })
  if (is.null(p_mod)) {
    row[["converged_flag"]] <- 0L
    row[["error_message"]] <- paste0("propensity: ", ifelse(is.na(err_msg), "failed", err_msg))
    row[["replication_seconds"]] <- as.numeric(proc.time()[["elapsed"]] - replication_start)
    append_csv_locked(as.data.frame(row, stringsAsFactors = FALSE)[COL_NAMES], out_path)
    return(invisible(NULL))
  }
  p <- colMeans(pnorm(p_mod$yhat.train)); p_test <- colMeans(pnorm(p_mod$yhat.test))
  # Augment X with the propensity exactly as the original does (X <- cbind(X, p)),
  # keeping clean copies X2/X2_test (no propensity) for the tau/B CF components.
  X <- cbind(X, p); X_test <- cbind(X_test, p_test)

  # ---- MVBCF via fast_bart() ----
  f_seed <- use_model_stream(f_stream, "mvbcf")
  row[["model_seed_mvbcf"]] <- f_seed
  t_mv <- system.time({
    mvbcf_mod <- tryCatch(
      # Signature (MVBCF_Code.cpp:659): X_con, y, Z, X_mod, X_con_test,
      # X_mod_test, alpha, beta, alpha_tau, beta_tau, sigma_mu, sigma_tau,
      # v_0, sigma_0, n_iter, n_tree, n_tree_tau, min_nodesize.
      # X_mod_test (the tau-part test matrix) is REQUIRED and is the clean
      # covariate matrix, exactly as GitHub_DGP1.R passes it.
      fast_bart(X, Y, cbind(Z, Z), X2, X_test, X2_test,
                0.95, 2, 0.25, 3,
                diag((mu_val)^2/n_tree_mu, 2), diag((tau_val)^2/n_tree_tau, 2),
                v_val, diag(wish_val, 2), n_iter, n_tree_mu, n_tree_tau, min_val),
      error = function(e) { err_msg <<- conditionMessage(e); NULL })
  })
  row[["fit_seconds_mvbcf"]] <- as.numeric(t_mv[["elapsed"]])
  if (!is.null(mvbcf_mod)) {
    raw1 <- mvbcf_mod$predictions_tau_test[, 1, -c(1:n_burn)]
    stopifnot(nrow(raw1) == n_test, ncol(raw1) == n_mcmc)
    raw2 <- mvbcf_mod$predictions_tau_test[, 2, -c(1:n_burn)]
    row[c("mvbcf_pehe1","mvbcf_pehe2")] <- list(sqrt(mean((Tau1_test - rowMeans(raw1))^2)),
                                                sqrt(mean((Tau2_test - rowMeans(raw2))^2)))
    row[c("mvbcf_ate1","mvbcf_ate2")] <- list(mean(rowMeans(raw1)), mean(rowMeans(raw2)))
    row[["mvbcf_bias1"]] <- row[["mvbcf_ate1"]] - mean(Tau1_test)
    row[["mvbcf_bias2"]] <- row[["mvbcf_ate2"]] - mean(Tau2_test)
    row[["mvbcf_crmse1"]] <- mean((Tau1_test - rowMeans(raw1))^2)
    row[["mvbcf_crmse2"]] <- mean((Tau2_test - rowMeans(raw2))^2)
    row[["mvbcf_atesq1"]] <- (row[["mvbcf_ate1"]] - mean(Tau1_test))^2
    row[["mvbcf_atesq2"]] <- (row[["mvbcf_ate2"]] - mean(Tau2_test))^2
    row[["mvbcf_crps1"]] <- mean(crps_norm(Tau1_test, mean = rowMeans(raw1), sd = rowSds(raw1)))
    row[["mvbcf_crps2"]] <- mean(crps_norm(Tau2_test, mean = rowMeans(raw2), sd = rowSds(raw2)))
    row[c("mvbcf_cov501","mvbcf_cov951","mvbcf_cov502","mvbcf_cov952")] <-
      list(mean(diag(apply(raw1, 1, in_cred, Tau1_test, 0.5))),
           mean(diag(apply(raw1, 1, in_cred, Tau1_test, 0.95))),
           mean(diag(apply(raw2, 1, in_cred, Tau2_test, 0.5))),
           mean(diag(apply(raw2, 1, in_cred, Tau2_test, 0.95))))
    row[c("mvbcf_wid501","mvbcf_wid951","mvbcf_wid502","mvbcf_wid952")] <-
      list(mean(apply(raw1, 1, cred_width, 0.5)), mean(apply(raw1, 1, cred_width, 0.95)),
           mean(apply(raw2, 1, cred_width, 0.5)), mean(apply(raw2, 1, cred_width, 0.95)))
    a1 <- colMeans(raw1); a2 <- colMeans(raw2)
    row[c("mvbcf_ate_l951","mvbcf_ate_u951","mvbcf_ate_l952","mvbcf_ate_u952")] <-
      list(quantile(a1, 0.025, names=FALSE), quantile(a1, 0.975, names=FALSE),
           quantile(a2, 0.025, names=FALSE), quantile(a2, 0.975, names=FALSE))
    row[c("mvbcf_ate_l501","mvbcf_ate_u501","mvbcf_ate_l502","mvbcf_ate_u502")] <-
      list(quantile(a1, 0.25, names=FALSE), quantile(a1, 0.75, names=FALSE),
           quantile(a2, 0.25, names=FALSE), quantile(a2, 0.75, names=FALSE))
    row[c("mvbcf_ate_wid951","mvbcf_ate_wid952","mvbcf_ate_wid501","mvbcf_ate_wid502")] <-
      list(row[["mvbcf_ate_u951"]]-row[["mvbcf_ate_l951"]], row[["mvbcf_ate_u952"]]-row[["mvbcf_ate_l952"]],
           row[["mvbcf_ate_u501"]]-row[["mvbcf_ate_l501"]], row[["mvbcf_ate_u502"]]-row[["mvbcf_ate_l502"]])
    row[c("mvbcf_ate_cov951","mvbcf_ate_cov952","mvbcf_ate_cov501","mvbcf_ate_cov502")] <-
      list(as.integer(mean(Tau1_test)>=row[["mvbcf_ate_l951"]] & mean(Tau1_test)<=row[["mvbcf_ate_u951"]]),
           as.integer(mean(Tau2_test)>=row[["mvbcf_ate_l952"]] & mean(Tau2_test)<=row[["mvbcf_ate_u952"]]),
           as.integer(mean(Tau1_test)>=row[["mvbcf_ate_l501"]] & mean(Tau1_test)<=row[["mvbcf_ate_u501"]]),
           as.integer(mean(Tau2_test)>=row[["mvbcf_ate_l502"]] & mean(Tau2_test)<=row[["mvbcf_ate_u502"]]))
  } else {
    row[["converged_flag"]] <- 0L
    row[["error_message"]] <- paste0("mvbcf: ", ifelse(is.na(err_msg), "failed", err_msg))
  }

  # ---- BCF (calibrated stochtree; outcome-specific prior scale) ----
  for (k in 1:2) {
    f_seed <- use_model_stream(f_stream, paste0("bcf", k))
    row[[paste0("model_seed_bcf", k)]] <- f_seed
    t_bc <- system.time({
      bc <- tryCatch(
        stochtree::bcf(X2, Z, Y[, k], X_test = X2_test,
                       prognostic_forest_params = list(num_trees = n_tree_mu, sigma2_leaf_init = sd(Y[, k])^2),
                       treatment_effect_forest_params = list(num_trees = n_tree_tau, sigma2_leaf_init = (0.375*sd(Y[, k]))^2),
                       Z_test = Z_test, propensity_train = p, propensity_test = p_test,
                       num_burnin = n_burn, num_mcmc = n_mcmc, num_gfr = 50,
                       nthread = nthread_global, random_seed = f_seed),
        error = function(e) { err_msg <<- conditionMessage(e); NULL })
    })
    row[[paste0("fit_seconds_bcf", k)]] <- as.numeric(t_bc[["elapsed"]])
    ss <- if (k == 1L) 1L else 2L
    tt <- if (ss == 1L) Tau1_test else Tau2_test
    if (!is.null(bc)) {
      raw <- bc$tau_hat_test                 # rows = test obs, cols = MCMC draws
      stopifnot(nrow(raw) == n_test, ncol(raw) == n_mcmc)
      tau <- rowMeans(raw)                   # per-observation posterior mean
      ate_draw <- colMeans(raw)              # ATE posterior draw vector (average over test)
      row[[paste0("bcf_pehe", ss)]] <- sqrt(mean((tt - tau)^2))
      row[[paste0("bcf_ate", ss)]] <- mean(ate_draw)
      row[[paste0("bcf_bias", ss)]] <- mean(ate_draw) - mean(tt)
      row[[paste0("bcf_crmse", ss)]] <- mean((tt - tau)^2)
      row[[paste0("bcf_atesq", ss)]] <- (mean(ate_draw) - mean(tt))^2
      row[[paste0("bcf_crps", ss)]] <- mean(crps_norm(tt, mean = tau, sd = rowSds(raw)))
      row[c(paste0("bcf_cov50", ss), paste0("bcf_cov95", ss))] <-
        list(mean(diag(apply(raw, 1, in_cred, tt, 0.5))), mean(diag(apply(raw, 1, in_cred, tt, 0.95))))
      row[c(paste0("bcf_wid50", ss), paste0("bcf_wid95", ss))] <-
        list(mean(apply(raw, 1, cred_width, 0.5)), mean(apply(raw, 1, cred_width, 0.95)))
      row[paste0("bcf_ate_l95", ss)] <- quantile(ate_draw, 0.025, names = FALSE)
      row[paste0("bcf_ate_u95", ss)] <- quantile(ate_draw, 0.975, names = FALSE)
      row[paste0("bcf_ate_l50", ss)] <- quantile(ate_draw, 0.25, names = FALSE)
      row[paste0("bcf_ate_u50", ss)] <- quantile(ate_draw, 0.75, names = FALSE)
      row[paste0("bcf_ate_wid95", ss)] <- row[[paste0("bcf_ate_u95", ss)]] - row[[paste0("bcf_ate_l95", ss)]]
      row[paste0("bcf_ate_wid50", ss)] <- row[[paste0("bcf_ate_u50", ss)]] - row[[paste0("bcf_ate_l50", ss)]]
      row[paste0("bcf_ate_cov95", ss)] <- as.integer(mean(tt) >= row[[paste0("bcf_ate_l95", ss)]] & mean(tt) <= row[[paste0("bcf_ate_u95", ss)]])
      row[paste0("bcf_ate_cov50", ss)] <- as.integer(mean(tt) >= row[[paste0("bcf_ate_l50", ss)]] & mean(tt) <= row[[paste0("bcf_ate_u50", ss)]])
    } else {
      row[["converged_flag"]] <- 0L
      row[["error_message"]] <- paste0("bcf", k, ": ", ifelse(is.na(err_msg), "failed", err_msg))
    }
  }

  # ---- BART (univariate bartCause, one model per outcome) ----
  # GitHub_DGP1.R:154-157 names the design columns V1..Vp before bartc(), so
  # that predict(bart_mod, X_test) can match train and test columns by name.
  # Without this bartCause matches positionally and warns (or errors) on an
  # unnamed matrix. Set them at the same point in the sequence as the original.
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  colnames(X_test) <- paste0("V", seq_len(ncol(X_test)))

  # bartCause returns `icate` as (posterior draws x test observations) -- the
  # TRANSPOSE of the stochtree / fast_bart convention. GitHub_DGP1.R:158-166,
  # 261, 266 therefore uses colMeans() for the per-observation tau, rowMeans()
  # for the ATE posterior, apply(margin = 2) for coverage/width and colSds()
  # for the CRPS scale. Keep those margins here; swapping them silently
  # computes every BART metric on the wrong axis.
  for (k in 1:2) {
    f_seed <- use_model_stream(f_stream, paste0("bart", k))
    row[[paste0("model_seed_bart", k)]] <- f_seed
    t_bt <- system.time({
      bm <- tryCatch(
        bartc(Y[, k], Z, X, p.scoreAsCovariate = FALSE, n.chains = 1, n.threads = 1,
              keepTrees = TRUE, n.trees = n_trees_total, seed = f_seed),
        error = function(e) { err_msg <<- conditionMessage(e); NULL })
    })
    row[[paste0("fit_seconds_bart", k)]] <- as.numeric(t_bt[["elapsed"]])
    if (!is.null(bm)) {
      ss <- if (k == 1L) 1L else 2L
      ic <- predict(bm, X_test, type = "icate")   # rows = draws, cols = test obs
      # Fail loudly rather than silently averaging over the wrong axis. Every
      # BART metric below depends on this orientation, and n_test (1000) and
      # n_mcmc (500) are both plausible dimensions, so a transposed matrix would
      # not error -- it would recycle and write meaningless numbers into every
      # shard. Costs microseconds; protects the whole campaign.
      stopifnot(ncol(ic) == n_test, nrow(ic) == n_mcmc)
      tt <- if (ss == 1L) Tau1_test else Tau2_test
      tau <- colMeans(ic)                         # per-observation posterior mean
      ate_draw <- rowMeans(ic)                    # per-draw ATE (average over test)
      row[[paste0("bart_pehe", ss)]] <- sqrt(mean((tt - tau)^2))
      row[[paste0("bart_ate", ss)]] <- mean(ate_draw)
      row[[paste0("bart_bias", ss)]] <- mean(ate_draw) - mean(tt)
      row[[paste0("bart_crmse", ss)]] <- mean((tt - tau)^2)
      row[[paste0("bart_atesq", ss)]] <- (mean(ate_draw) - mean(tt))^2
      row[[paste0("bart_crps", ss)]] <- mean(crps_norm(tt, mean = tau, sd = colSds(ic)))
      row[c(paste0("bart_cov50", ss), paste0("bart_cov95", ss))] <-
        list(mean(diag(apply(ic, 2, in_cred, tt, 0.5))), mean(diag(apply(ic, 2, in_cred, tt, 0.95))))
      row[c(paste0("bart_wid50", ss), paste0("bart_wid95", ss))] <-
        list(mean(apply(ic, 2, cred_width, 0.5)), mean(apply(ic, 2, cred_width, 0.95)))
      row[c(paste0("bart_ate_l95", ss), paste0("bart_ate_u95", ss))] <-
        list(quantile(ate_draw, 0.025, names = FALSE), quantile(ate_draw, 0.975, names = FALSE))
      row[c(paste0("bart_ate_l50", ss), paste0("bart_ate_u50", ss))] <-
        list(quantile(ate_draw, 0.25, names = FALSE), quantile(ate_draw, 0.75, names = FALSE))
      row[paste0("bart_ate_wid95", ss)] <- quantile(ate_draw, 0.975, names = FALSE) - quantile(ate_draw, 0.025, names = FALSE)
      row[paste0("bart_ate_wid50", ss)] <- quantile(ate_draw, 0.75, names = FALSE) - quantile(ate_draw, 0.25, names = FALSE)
      row[paste0("bart_ate_cov95", ss)] <- as.integer(mean(tt) >= quantile(ate_draw, 0.025, names = FALSE) & mean(tt) <= quantile(ate_draw, 0.975, names = FALSE))
      row[paste0("bart_ate_cov50", ss)] <- as.integer(mean(tt) >= quantile(ate_draw, 0.25, names = FALSE) & mean(tt) <= quantile(ate_draw, 0.75, names = FALSE))
    } else {
      row[["converged_flag"]] <- 0L
      row[["error_message"]] <- paste0("bart", k, ": ", ifelse(is.na(err_msg), "failed", err_msg))
    }
  }

  # ---- MVBART (skewBART) ----
  f_seed <- use_model_stream(f_stream, "mvbart")
  row[["model_seed_mvbart"]] <- f_seed
  t_mb <- system.time({
    mb <- tryCatch({
      hypers <- Hypers(X = cbind(X, Z), Y = Y, num_tree = n_trees_total)
      opts <- Opts(num_burn = n_burn, num_save = n_mcmc)
      X_test_z <- rbind(cbind(X_test, 0*Z_test), cbind(X_test, 0*Z_test + 1))
      MultiskewBART(X = cbind(X, Z), Y = Y, test_X = X_test_z,
                    hypers = hypers, opts = opts, do_skew = FALSE)
    }, error = function(e) { err_msg <<- conditionMessage(e); NULL })
  })
  row[["fit_seconds_mvbart"]] <- as.numeric(t_mb[["elapsed"]])
  if (!is.null(mb)) {
    z0 <- mb$y_hat_test[1:n_test, , ]; z1 <- mb$y_hat_test[-c(1:n_test), , ]
    z10 <- z1 - z0
    stopifnot(dim(z10)[1] == n_test, dim(z10)[2] == 2L, dim(z10)[3] == n_mcmc)
    ic1 <- rowMeans(z10[, 1, ]); ic2 <- rowMeans(z10[, 2, ])
    row[c("mvbart_pehe1","mvbart_pehe2")] <-
      list(sqrt(mean((Tau1_test - ic1)^2)), sqrt(mean((Tau2_test - ic2)^2)))
    row[c("mvbart_ate1","mvbart_ate2")] <- list(mean(ic1), mean(ic2))
    row[["mvbart_bias1"]] <- mean(ic1) - mean(Tau1_test)
    row[["mvbart_bias2"]] <- mean(ic2) - mean(Tau2_test)
    row[["mvbart_crmse1"]] <- mean((Tau1_test - ic1)^2)
    row[["mvbart_crmse2"]] <- mean((Tau2_test - ic2)^2)
    row[["mvbart_atesq1"]] <- (mean(ic1) - mean(Tau1_test))^2
    row[["mvbart_atesq2"]] <- (mean(ic2) - mean(Tau2_test))^2
    # DEVIATION FROM UPSTREAM (deliberate, documented). GitHub_DGP1.R:262 and
    # :293 use colSds(z_1_0_preds[,k,]) here, but that slice is
    # (test obs x draws): colSds returns one value PER DRAW (length n_mcmc),
    # while crps_norm needs one sd PER TEST OBSERVATION (length n_test). With
    # n_test = 1000 and n_mcmc = 500 the upstream call silently recycles and
    # the mvbart CRPS column is meaningless. Every other model in the upstream
    # script does supply a per-observation sd. CRPS is promoted to the MCS loss
    # for the uncertainty comparison (plan 1.7b, 2.3), so this one must be
    # right: use rowSds, matching the mvbcf branch above.
    row[["mvbart_crps1"]] <- mean(crps_norm(Tau1_test, mean = ic1, sd = rowSds(z10[, 1, ])))
    row[["mvbart_crps2"]] <- mean(crps_norm(Tau2_test, mean = ic2, sd = rowSds(z10[, 2, ])))
    row[c("mvbart_cov501","mvbart_cov951","mvbart_cov502","mvbart_cov952")] <-
      list(mean(diag(apply(z10[, 1, ], 1, in_cred, Tau1_test, 0.5))),
           mean(diag(apply(z10[, 1, ], 1, in_cred, Tau1_test, 0.95))),
           mean(diag(apply(z10[, 2, ], 1, in_cred, Tau2_test, 0.5))),
           mean(diag(apply(z10[, 2, ], 1, in_cred, Tau2_test, 0.95))))
    row[c("mvbart_wid501","mvbart_wid951","mvbart_wid502","mvbart_wid952")] <-
      list(mean(apply(z10[, 1, ], 1, cred_width, 0.5)), mean(apply(z10[, 1, ], 1, cred_width, 0.95)),
           mean(apply(z10[, 2, ], 1, cred_width, 0.5)), mean(apply(z10[, 2, ], 1, cred_width, 0.95)))
    a1 <- colMeans(z10[, 1, ]); a2 <- colMeans(z10[, 2, ])
    row[c("mvbart_ate_l951","mvbart_ate_u951","mvbart_ate_l952","mvbart_ate_u952")] <-
      list(quantile(a1, 0.025, names=FALSE), quantile(a1, 0.975, names=FALSE),
           quantile(a2, 0.025, names=FALSE), quantile(a2, 0.975, names=FALSE))
    row[c("mvbart_ate_l501","mvbart_ate_u501","mvbart_ate_l502","mvbart_ate_u502")] <-
      list(quantile(a1, 0.25, names=FALSE), quantile(a1, 0.75, names=FALSE),
           quantile(a2, 0.25, names=FALSE), quantile(a2, 0.75, names=FALSE))
    row[c("mvbart_ate_wid951","mvbart_ate_wid952","mvbart_ate_wid501","mvbart_ate_wid502")] <-
      list(row[["mvbart_ate_u951"]]-row[["mvbart_ate_l951"]], row[["mvbart_ate_u952"]]-row[["mvbart_ate_l952"]],
           row[["mvbart_ate_u501"]]-row[["mvbart_ate_l501"]], row[["mvbart_ate_u502"]]-row[["mvbart_ate_l502"]])
    row[c("mvbart_ate_cov951","mvbart_ate_cov952","mvbart_ate_cov501","mvbart_ate_cov502")] <-
      list(as.integer(mean(Tau1_test)>=row[["mvbart_ate_l951"]] & mean(Tau1_test)<=row[["mvbart_ate_u951"]]),
           as.integer(mean(Tau2_test)>=row[["mvbart_ate_l952"]] & mean(Tau2_test)<=row[["mvbart_ate_u952"]]),
           as.integer(mean(Tau1_test)>=row[["mvbart_ate_l501"]] & mean(Tau1_test)<=row[["mvbart_ate_u501"]]),
           as.integer(mean(Tau2_test)>=row[["mvbart_ate_l502"]] & mean(Tau2_test)<=row[["mvbart_ate_u502"]]))
  } else {
    row[["converged_flag"]] <- 0L
    row[["error_message"]] <- paste0("mvbart: ", ifelse(is.na(err_msg), "failed", err_msg))
  }

  row[["hostname"]] <- Sys.info()[["nodename"]]
  row[["git_sha"]] <- git_sha
  row[["session_hash"]] <- session_hash
  row[["replication_seconds"]] <- as.numeric(proc.time()[["elapsed"]] - replication_start)
  df <- as.data.frame(row, stringsAsFactors = FALSE)[COL_NAMES]
  append_csv_locked(df, out_path)
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Static run metadata captured once per process
# -----------------------------------------------------------------------------
git_sha <- Sys.getenv("TISCA_GIT_SHA", "not-a-git-build")
session_hash <- digest::digest(sessionInfo())
cat("run_cell dgp=", dgp_arg, " n=", n_arg, " seeds ", s_start, "..", s_end,
    " out=", out_path, " cores=", cores, " mode=", mode,
    " session_hash=", session_hash, "\n", sep = "")

# -----------------------------------------------------------------------------
# Pre-build the (cheap) stream list for the whole shard. The CLI seed range is
# the small per-master index range (confirmatory 0.., pilot 0..J0-1); the emit
# labels used in the `seed` column are shifted for pilots so blocks are visually
# disjoint (confirmatory 0.., pilot 1_000_001..1_000_000+J0). seed_cell_master
# records the master actually used, so a shard can always rebuild its streams.
# -----------------------------------------------------------------------------
stream <- make_streams(cell_master, s_start, s_end)
# NB the offset is applied to the WHOLE index range. `s_start + 1000000L + 1L`
# collapses the pilot to a length-1 vector, so a 50-seed pilot shard would run
# exactly one replication and the calibration gate would be sized at J0 = 1.
seeds_emit <- if (mode == "pilot") (s_start:s_end) + 1000000L + 1L else s_start:s_end
seeds_emit <- as.integer(seeds_emit)
stopifnot(length(seeds_emit) == s_end - s_start + 1L,
          length(seeds_emit) == length(stream$data),
          length(seeds_emit) == length(stream$fit),
          !anyDuplicated(seeds_emit))

cat("running seeds", s_start, "..", s_end, "(", length(seeds_emit), " reps)\n", sep = "")
idx <- seq_along(seeds_emit)
res <- mclapply(idx, function(i)
  run_one(seeds_emit[i], stream$data[[i]], stream$fit[[i]], out_path),
  mc.cores = cores, mc.set.seed = FALSE)
done <- file.exists(out_path)
cat("done; shard rows appended to", out_path, "(exists=", done, ")\n")
