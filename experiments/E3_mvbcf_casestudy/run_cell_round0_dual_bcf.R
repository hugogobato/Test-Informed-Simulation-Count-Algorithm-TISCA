#!/usr/bin/env Rscript
#
# run_cell_round0_dual_bcf.R -- paired Round 0 BCF benchmark.
#
# Each replication generates one DGP/test data set and one BART propensity
# estimate, then fits both the closest paper-specification stochtree::bcf and
# the original bcf package used by GitHub_DGP1.R. The two methods are kept in a
# separate driver and CSV schema so existing E3 checkpoints are unchanged.
#
# CLI:
#   Rscript run_cell_round0_dual_bcf.R <dgp> <n> <seed_start> <seed_end> \
#       --out <csv> [--cores <k>] [--mode pilot|confirmatory]
#
suppressPackageStartupMessages({
  library(parallel)
  library(mvtnorm)
  library(dbarts)
  library(stochtree)
  library(bcf)
  library(matrixStats)
  library(digest)
})

# -----------------------------------------------------------------------------
# Argument parsing and constants
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("expected <dgp> <n> <seed_start> <seed_end>")
dgp_arg <- as.integer(args[[1]])
n_arg <- as.integer(args[[2]])
s_start <- as.integer(args[[3]])
s_end <- as.integer(args[[4]])
a <- args[-(1:4)]
get_opt <- function(key, default = NULL) {
  for (i in seq_along(a)) if (a[i] == key) return(a[i + 1L])
  default
}
out_path <- get_opt("--out", stop("--out required"))
cores <- as.integer(get_opt("--cores", "1"))
mode <- get_opt("--mode", "confirmatory")

stopifnot(!is.na(dgp_arg), dgp_arg %in% 1:3,
          !is.na(n_arg), n_arg %in% c(100L, 500L),
          !is.na(s_start), !is.na(s_end), s_start >= 0L, s_end >= s_start,
          !is.na(cores), cores >= 1L,
          mode %in% c("pilot", "confirmatory"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

n_test <- 1000L
n_tree_mu <- 50L
n_tree_tau <- 20L
n_burn <- 500L
n_mcmc <- 500L
nthread_global <- 1L
cell_master <- if (mode == "pilot") 2L else 1L

# -----------------------------------------------------------------------------
# Reproducible L'Ecuyer-CMRG streams. The data stream is shared by both
# methods. Fit streams are disjoint by model.
# -----------------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
make_streams <- function(master, j_from, j_to) {
  set.seed(master)
  s <- .Random.seed
  n <- j_to - j_from + 1L
  data_s <- fit_s <- vector("list", n)
  idx <- 0L
  for (j in 0:j_to) {
    d <- s; s <- nextRNGStream(s)
    f <- s; s <- nextRNGStream(s)
    if (j >= j_from) {
      idx <- idx + 1L
      data_s[[idx]] <- d
      fit_s[[idx]] <- f
    }
  }
  list(data = data_s, fit = fit_s)
}
set_rand <- function(x) assign(".Random.seed", x, envir = .GlobalEnv)
model_stream <- function(f_stream, i) {
  s <- f_stream
  if (i > 0L) for (k in seq_len(i)) s <- nextRNGStream(s)
  s
}
MODEL_SUBSTREAM <- c(propensity = 0L, stochtree1 = 1L, stochtree2 = 2L,
                     paper1 = 3L, paper2 = 4L)
use_model_stream <- function(f_stream, name) {
  set_rand(model_stream(f_stream, MODEL_SUBSTREAM[[name]]))
  sample.int(.Machine$integer.max, 1L)
}

# -----------------------------------------------------------------------------
# DGPs and test-set formulas from the authors' GitHub_DGP1.R--GitHub_DGP3.R.
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
  Z_test <- rbinom(n_test, 1, X4)
  list(X_test = X_test, Z_test = Z_test,
       Tau1_test = Tau1_test, Tau2_test = Tau2_test)
}

# -----------------------------------------------------------------------------
# Metrics. Stochtree is observation x draw; bcf is draw x observation.
# -----------------------------------------------------------------------------
in_cred <- function(samples, value, interval) {
  q <- quantile(samples, c((1 - interval)/2, 1 - (1 - interval)/2), names = FALSE)
  value >= q[[1]] && value <= q[[2]]
}
cred_width <- function(samples, interval) {
  diff(quantile(samples, c((1 - interval)/2, 1 - (1 - interval)/2), names = FALSE))
}
crps_normal <- function(observed, mean, sd) {
  # Algebraically identical to scoringRules::crps_norm, kept local so this
  # diagnostic driver can also smoke-test in a minimal R installation.
  sd <- pmax(as.numeric(sd), .Machine$double.eps)
  z <- (observed - mean) / sd
  sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
}
set_metrics <- function(row, prefix, raw, truth, orientation) {
  if (orientation == "observation_draw") {
    stopifnot(nrow(raw) == n_test, ncol(raw) == n_mcmc)
    tau <- rowMeans(raw); ate_draw <- colMeans(raw); posterior_sd <- rowSds(raw)
    coverage <- function(level) mean(vapply(seq_len(n_test), function(i)
      in_cred(raw[i, ], truth[i], level), logical(1)))
    width <- function(level) mean(vapply(seq_len(n_test), function(i)
      cred_width(raw[i, ], level), numeric(1)))
  } else {
    stopifnot(nrow(raw) == n_mcmc, ncol(raw) == n_test)
    tau <- colMeans(raw); ate_draw <- rowMeans(raw); posterior_sd <- colSds(raw)
    coverage <- function(level) mean(vapply(seq_len(n_test), function(i)
      in_cred(raw[, i], truth[i], level), logical(1)))
    width <- function(level) mean(vapply(seq_len(n_test), function(i)
      cred_width(raw[, i], level), numeric(1)))
  }
  row[[paste0(prefix, "_pehe")]] <- sqrt(mean((truth - tau)^2))
  row[[paste0(prefix, "_ate")]] <- mean(ate_draw)
  row[[paste0(prefix, "_bias")]] <- mean(ate_draw) - mean(truth)
  row[[paste0(prefix, "_crmse")]] <- mean((truth - tau)^2)
  row[[paste0(prefix, "_atesq")]] <- (mean(ate_draw) - mean(truth))^2
  row[[paste0(prefix, "_crps")]] <- mean(crps_normal(
    truth, mean = tau, sd = posterior_sd))
  row[[paste0(prefix, "_cov50")]] <- coverage(0.50)
  row[[paste0(prefix, "_cov95")]] <- coverage(0.95)
  row[[paste0(prefix, "_wid50")]] <- width(0.50)
  row[[paste0(prefix, "_wid95")]] <- width(0.95)
  row[[paste0(prefix, "_ate_l95")]] <- quantile(ate_draw, 0.025, names = FALSE)
  row[[paste0(prefix, "_ate_u95")]] <- quantile(ate_draw, 0.975, names = FALSE)
  row[[paste0(prefix, "_ate_cov95")]] <- as.integer(
    mean(truth) >= row[[paste0(prefix, "_ate_l95")]] &&
    mean(truth) <= row[[paste0(prefix, "_ate_u95")]])
  row
}

# bcf computes this lambda when lambda=NULL. Match that residual IG prior in
# stochtree's standardized outcome units as closely as its API permits.
paper_lambda <- function(y, z, x_control, pihat, nu = 3, sigq = 0.9) {
  ystd <- (y - mean(y)) / sd(y)
  x_c <- cbind(x_control, pihat)
  sighat <- summary(lm(ystd ~ z + as.matrix(x_c)))$sigma
  (sighat^2 * qchisq(1 - sigq, nu)) / nu
}

# The paper package serializes trees before its prediction method is called.
paper_bcf_predict <- function(y, z, x_control, x_moderate, pihat,
                              x_predict_control, x_predict_moderate,
                              pi_pred, z_pred, seed) {
  tree_dir <- tempfile("tisca_paper_bcf_")
  dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tree_dir, recursive = TRUE, force = TRUE), add = TRUE)
  fit <- bcf::bcf(
    y = y, z = z, x_control = x_control, x_moderate = x_moderate,
    pihat = pihat, random_seed = seed, n_chains = 1, n_threads = 1,
    nburn = n_burn, nsim = n_mcmc,
    ntree_control = n_tree_mu, ntree_moderate = n_tree_tau,
    sd_control = sd(y), sd_moderate = 0.375 * sd(y),
    include_pi = "control", use_muscale = FALSE, use_tauscale = FALSE,
    no_output = FALSE, save_tree_directory = tree_dir,
    log_file = file.path(tree_dir, "fit.log"), verbose = FALSE)
  # Calling the namespace function is equivalent to predict.bcf() and is
  # robust to installations that do not expose the S3 method in methods().
  pred <- bcf:::predict.bcf(
    object = fit, x_predict_control = x_predict_control,
    x_predict_moderate = x_predict_moderate, pi_pred = pi_pred,
    z_pred = z_pred, save_tree_directory = tree_dir,
    log_file = file.path(tree_dir, "predict.log"), n_cores = 1,
    n_threads = 1)
  pred$tau
}

# -----------------------------------------------------------------------------
# Checkpointing and output schema
# -----------------------------------------------------------------------------
append_csv <- function(df, path) {
  has_contents <- file.exists(path) && file.info(path)$size > 0L
  if (has_contents) {
    header <- names(read.csv(path, nrows = 0L, check.names = FALSE))
    if (!identical(header, names(df))) stop("dual BCF checkpoint schema mismatch")
  }
  write.table(df, path, row.names = FALSE, sep = ",",
              col.names = !has_contents, append = has_contents)
}
append_csv_locked <- function(df, path) {
  lock <- paste0(path, ".lock")
  started <- proc.time()[["elapsed"]]
  repeat {
    if (dir.create(lock, showWarnings = FALSE)) break
    if (proc.time()[["elapsed"]] - started > 300) stop("checkpoint lock timeout")
    Sys.sleep(0.05)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  append_csv(df, path)
}

make_stub <- function() {
  row <- list(
    seed = NA_integer_, n = NA_integer_, dgp = NA_integer_,
    seed_cell_master = NA_integer_, rng_kind = "L'Ecuyer-CMRG",
    seed_data_hash = NA_character_, seed_fit_hash = NA_character_,
    model_seed_propensity = NA_integer_,
    model_seed_stochtree_bcf1 = NA_integer_, model_seed_stochtree_bcf2 = NA_integer_,
    model_seed_paper_bcf1 = NA_integer_, model_seed_paper_bcf2 = NA_integer_,
    paper_lambda1 = NA_real_, paper_lambda2 = NA_real_,
    fit_seconds_stochtree_bcf1 = NA_real_, fit_seconds_stochtree_bcf2 = NA_real_,
    fit_seconds_paper_bcf1 = NA_real_, fit_seconds_paper_bcf2 = NA_real_,
    converged_flag = NA_integer_, error_message = NA_character_,
    hostname = NA_character_, git_sha = NA_character_, session_hash = NA_character_,
    seq_phase = NA_character_, replication_seconds = NA_real_
  )
  measures <- c("pehe", "ate", "bias", "crmse", "atesq", "crps",
                "cov50", "cov95", "wid50", "wid95", "ate_l95", "ate_u95", "ate_cov95")
  for (prefix in c("stochtree_bcf1", "stochtree_bcf2", "paper_bcf1", "paper_bcf2")) {
    for (measure in measures) row[[paste0(prefix, "_", measure)]] <- NA_real_
  }
  row
}
STUB <- make_stub()
COL_NAMES <- names(STUB)

run_one <- function(seed_emit, d_stream, f_stream, out_path) {
  started <- proc.time()[["elapsed"]]
  row <- STUB
  row[c("seed", "n", "dgp", "seed_cell_master", "seq_phase", "hostname", "git_sha", "session_hash")] <-
    list(seed_emit, n_arg, dgp_arg, cell_master, mode, Sys.info()[["nodename"]], git_sha, session_hash)
  row[["converged_flag"]] <- 1L
  errors <- character()
  fail <- function(tag, e) {
    row[["converged_flag"]] <<- 0L
    errors <<- c(errors, paste0(tag, ": ", conditionMessage(e)))
    row[["error_message"]] <<- paste(errors, collapse = " | ")
  }

  set_rand(d_stream)
  row[["seed_data_hash"]] <- paste(.Random.seed, collapse = "_")
  d <- generate_data(dgp_arg, n_arg)
  t <- generate_test(dgp_arg, n_test)
  X2 <- d$X; Y <- d$Y; Z <- d$Z
  X2_test <- t$X_test; Z_test <- t$Z_test
  Tau_test <- list(t$Tau1_test, t$Tau2_test)

  # Both BCF implementations use this shared propensity estimate.
  set_rand(f_stream)
  row[["seed_fit_hash"]] <- paste(.Random.seed, collapse = "_")
  p_seed <- use_model_stream(f_stream, "propensity")
  row[["model_seed_propensity"]] <- p_seed
  p_mod <- tryCatch(
    dbarts::bart(x.train = X2, y.train = Z, x.test = X2_test, k = 3,
                 nthread = nthread_global, seed = p_seed),
    error = function(e) { fail("propensity", e); NULL })
  if (is.null(p_mod)) {
    row[["replication_seconds"]] <- proc.time()[["elapsed"]] - started
    append_csv_locked(as.data.frame(row, stringsAsFactors = FALSE)[COL_NAMES], out_path)
    return(invisible(NULL))
  }
  p <- colMeans(pnorm(p_mod$yhat.train))
  p_test <- colMeans(pnorm(p_mod$yhat.test))

  # The original paper script passes X_test after cbind(X_test, p_test), even
  # though predict.bcf then appends pi_pred for include_pi="control". Preserve
  # this literal paper-code behavior in the paper benchmark. Stochtree gets
  # clean X2_test and propensity_test separately.
  X_test_paper <- cbind(X2_test, p_test)
  paper_seed <- use_model_stream(f_stream, "paper1")
  row[["model_seed_paper_bcf1"]] <- paper_seed
  row[["model_seed_paper_bcf2"]] <- paper_seed

  for (k in 1:2) {
    truth <- Tau_test[[k]]
    y <- Y[, k]
    lambda <- tryCatch(paper_lambda(y, Z, X2, p),
                       error = function(e) { fail(paste0("lambda", k), e); NA_real_ })
    row[[paste0("paper_lambda", k)]] <- lambda

    sto_seed <- use_model_stream(f_stream, paste0("stochtree", k))
    row[[paste0("model_seed_stochtree_bcf", k)]] <- sto_seed
    elapsed <- system.time({
      sto_fit <- tryCatch(
        stochtree::bcf(
          X_train = X2, Z_train = Z, y_train = y,
          propensity_train = p, X_test = X2_test, Z_test = Z_test,
          propensity_test = p_test, num_gfr = 0, num_burnin = n_burn,
          num_mcmc = n_mcmc,
          general_params = list(
            standardize = TRUE, sample_sigma2_global = TRUE,
            sigma2_global_shape = 3 / 2,
            sigma2_global_scale = if (is.finite(lambda)) 3 * lambda / 2 else 1,
            propensity_covariate = "prognostic", adaptive_coding = FALSE,
            num_chains = 1, num_threads = nthread_global, random_seed = sto_seed),
          prognostic_forest_params = list(
            num_trees = n_tree_mu, alpha = 0.95, beta = 2,
            min_samples_leaf = 1, max_depth = -1,
            sample_sigma2_leaf = FALSE, sigma2_leaf_init = 1^2 / n_tree_mu),
          treatment_effect_forest_params = list(
            num_trees = n_tree_tau, alpha = 0.25, beta = 3,
            min_samples_leaf = 1, max_depth = -1,
            sample_sigma2_leaf = FALSE,
            sigma2_leaf_init = 0.375^2 / n_tree_tau)),
        error = function(e) { fail(paste0("stochtree", k), e); NULL })
    })
    row[[paste0("fit_seconds_stochtree_bcf", k)]] <- elapsed[["elapsed"]]
    if (!is.null(sto_fit)) {
      row <- set_metrics(row, paste0("stochtree_bcf", k), sto_fit$tau_hat_test,
                         truth, "observation_draw")
    }

    elapsed <- system.time({
      paper_tau <- tryCatch(
        paper_bcf_predict(
          y = y, z = Z, x_control = X2, x_moderate = X2, pihat = p,
          x_predict_control = X_test_paper, x_predict_moderate = X2_test,
          pi_pred = p_test, z_pred = Z_test, seed = paper_seed),
        error = function(e) { fail(paste0("paper_bcf", k), e); NULL })
    })
    row[[paste0("fit_seconds_paper_bcf", k)]] <- elapsed[["elapsed"]]
    if (!is.null(paper_tau)) {
      row <- set_metrics(row, paste0("paper_bcf", k), paper_tau,
                         truth, "draw_observation")
    }
  }

  row[["replication_seconds"]] <- proc.time()[["elapsed"]] - started
  if (identical(row[["converged_flag"]], 1L)) {
    character_columns <- c("rng_kind", "seed_data_hash", "seed_fit_hash",
                           "error_message", "hostname", "git_sha", "session_hash", "seq_phase")
    numeric_columns <- setdiff(COL_NAMES, character_columns)
    nonfinite <- numeric_columns[!vapply(row[numeric_columns], function(x)
      length(x) == 1L && is.numeric(x) && is.finite(x), logical(1))]
    if (length(nonfinite)) {
      row[["converged_flag"]] <- 0L
      row[["error_message"]] <- paste0("non-finite outputs: ", paste(nonfinite, collapse = ", "))
    }
  }
  append_csv_locked(as.data.frame(row, stringsAsFactors = FALSE)[COL_NAMES], out_path)
  invisible(NULL)
}

git_sha <- Sys.getenv("TISCA_GIT_SHA", "not-a-git-build")
session_hash <- digest::digest(sessionInfo())
cat("dual BCF dgp=", dgp_arg, " n=", n_arg, " seeds ", s_start, "..", s_end,
    " out=", out_path, " cores=", cores, " mode=", mode,
    " session_hash=", session_hash, "\n", sep = "")
stream <- make_streams(cell_master, s_start, s_end)
seeds_emit <- if (mode == "pilot") (s_start:s_end) + 1000001L else s_start:s_end
seeds_emit <- as.integer(seeds_emit)
stopifnot(length(seeds_emit) == length(stream$data), length(stream$data) == length(stream$fit),
          !anyDuplicated(seeds_emit))
idx <- seq_along(seeds_emit)
res <- mclapply(idx, function(i)
  run_one(seeds_emit[[i]], stream$data[[i]], stream$fit[[i]], out_path),
  mc.cores = cores, mc.set.seed = FALSE)
worker_errors <- which(vapply(res, function(x) inherits(x, "try-error"), logical(1)))
if (length(worker_errors)) stop("worker errors: ", paste(seeds_emit[worker_errors], collapse = ", "))
cat("done; dual BCF rows appended to ", out_path, "\n", sep = "")
