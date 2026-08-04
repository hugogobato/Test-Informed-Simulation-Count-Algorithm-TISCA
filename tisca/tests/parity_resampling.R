#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# TISCA v2 -- resampling parity runner (P2-T2, second tier).
#
# Compares R against Python on the STOCHASTIC routines within bootstrap MCSE:
#   * studentized_paired_bootstrap (p-value / CI, mu and skewed and null data)
#   * Model Confidence Set (kept set / elimination path)
#   * Reality Check & SPA
#   * Romano-Wolf stepdown (rejection pattern)
#
# Because R (Mersenne-Twister) and Python (PCG64) use different RNGs, the exact
# resamples differ; the plan specifies agreement to bootstrap MCSE. For the
# bootstrap p-value the MCSE is sqrt(p(1-p)/B); we check |p_R - p_Py| <=
# 4 * MCSE(p_R) as a conservative envelope (4 sigma). MCS/RW set membership is
# compared on the same shared datasets and reported; small disagreements on a
# borderline model are documented rather than failed, but the kept-set should
# generally coincide because both use the true distinction.
# ----------------------------------------------------------------------------

script_dir <- function() {
  args <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(f)) return(dirname(normalizePath(f[1])))
  getwd()
}
dir <- script_dir()
setwd(dir)
source("helper_tisca.R")
source_tisca_R()

set.seed(2026)

## ---- build the shared spec JSON --------------------------------------
spec <- list(
  bootstrap = list(
    mu  = list(D = rnorm(120, 0.5, 1), alternative = "two-sided", B = 399L, seed = 1),
    skw = list(D = rexp(120) - 1,      alternative = "two-sided", B = 399L, seed = 2),
    nul = list(D = rnorm(120, 0.0, 1), alternative = "two-sided", B = 399L, seed = 3)
  ),
  mcs = list(
    m3 = list(L = cbind(rnorm(300,0.5,1), rnorm(300,2,1), rnorm(300,3,1)),
              model_names = c("m1","m2","m3"),
              B = 199L, alpha = 0.05, statistic = "Tmax", seed = 5),
    m2 = list(L = cbind(rnorm(200,1,1), rnorm(200,4,1)),
              model_names = c("m1","m2"),
              B = 99L,  alpha = 0.05, statistic = "TR",   seed = 6)
  ),
  spa = list(
    win  = list(L = cbind(rnorm(200,1,1), rnorm(200,2.5,1), rnorm(200,3,1)),
                champion = 0L, B = 199L, seed = 7),
    nul  = list(L = cbind(rnorm(200,1,1), rnorm(200,1,1)),
                champion = 0L, B = 199L, seed = 8)
  ),
  rw = list(
    two_true = list(D = cbind(rnorm(300,0.5,1), rnorm(300,0.6,1), rnorm(300,0,1)),
                    B = 199L, alpha = 0.05, seed = 9)
  )
)

writeLines(jsonlite::toJSON(spec, auto_unbox = TRUE), "parity_resampling_spec.json")

## ---- run the python evaluator ----------------------------------------
py <- Sys.which("python3"); if (!nzchar(py)) py <- Sys.which("python")
cmd <- sprintf("%s %s --spec %s --out %s", py,
               shQuote(file.path(dir, "parity_eval_python_resampling.py")),
               shQuote("parity_resampling_spec.json"),
               shQuote("parity_resampling_python.json"))
system(cmd, intern = TRUE)

pyres <- jsonlite::fromJSON("parity_resampling_python.json", simplifyVector = TRUE)

if (!is.null(pyres$`_PENDING_`)) {
  cat("[PENDING] Python resampling reference not available:", pyres$`_PENDING_`, "\n")
  quit(status = 0)
}

n_checks <- 0; n_pass <- 0
note <- character(0)

## ---- bootstrap parity within MCSE ------------------------------------
for (nm in names(spec$bootstrap)) {
  cfg <- spec$bootstrap[[nm]]
  altR <- gsub("-", ".", cfg$alternative, fixed = TRUE)   # "two-sided" -> "two.sided"
  rr <- studentized_paired_bootstrap(cfg$D, B = cfg$B, alpha = 0.05,
                                     alternative = altR, seed = cfg$seed)
  pyr <- pyres$bootstrap[[nm]]
  p_R <- rr$p_value; p_Py <- pyr$p_value
  mcse <- sqrt(p_R * (1 - p_R) / cfg$B)
  tol <- 4 * max(mcse, 1e-3)
  ok <- abs(p_R - p_Py) <= tol
  n_checks <- n_checks + 1; if (ok) n_pass <- n_pass + 1
  cat(sprintf("bootstrap p %-8s R=%.4f Py=%.4f mcse=%.4f -> %s\n",
              nm, p_R, p_Py, mcse, if (ok) "OK" else "MISMATCH"))
}

## ---- MCS parity --------------------------------------------------------
stat_map <- c("Tmax" = "T_max", "T_R" = "T_R", "TR" = "T_R", "T_max" = "T_max")
for (nm in names(spec$mcs)) {
  cfg <- spec$mcs[[nm]]
  L <- cfg$L
  colnames(L) <- paste0("m", seq_len(ncol(L)))
  rstat <- unname(stat_map[[cfg$statistic]])
  r <- mcs(L, B = cfg$B, alpha = cfg$alpha, statistic = rstat, seed = cfg$seed)
  pyr <- pyres$mcs[[nm]]
  r_kept <- sort(r$models_kept)
  py_kept <- sort(as.character(pyr$models_kept))
  same <- identical(r_kept, py_kept)
  n_checks <- n_checks + 1; if (same) n_pass <- n_pass + 1
  cat(sprintf("MCS %-8s R_kept=%s Py_kept=%s -> %s\n", nm,
              paste(r_kept, collapse=","), paste(py_kept, collapse=","),
              if (same) "OK" else "DIFFER"))
}

## ---- SPA / Reality check ----------------------------------------------
for (nm in names(spec$spa)) {
  cfg <- spec$spa[[nm]]
  L <- cfg$L; colnames(L) <- paste0("m", seq_len(ncol(L)))
  ## R reality_check expects champ by name; P:champ index 0
  champ_name <- colnames(L)[cfg$champion + 1L]
  r_rc <- reality_check(L, champ_name, B = cfg$B, seed = cfg$seed)$p_value
  r_sp <- spa_test(L, champ_name, B = cfg$B, seed = cfg$seed)$p_value
  pyr <- pyres$spa[[nm]]
  for (tag in c("rc", "spa")) {
    pR <- if (tag == "rc") r_rc else r_sp
    pPy <- if (tag == "rc") pyr$rc_p else pyr$spa_p
    tolv <- 4 * max(sqrt(pR * (1 - pR) / cfg$B), 1e-3)
    ok <- is.finite(pR) && is.finite(pPy) && abs(pR - pPy) <= tolv
    n_checks <- n_checks + 1; if (ok) n_pass <- n_pass + 1
    cat(sprintf("%s-%-8s R=%.4f Py=%.4f -> %s\n", tag, nm, pR, pPy,
                if (ok) "OK" else "MISMATCH"))
  }
}

## ---- Romano-Wolf ------------------------------------------------------
for (nm in names(spec$rw)) {
  cfg <- spec$rw[[nm]]
  r_rw <- romano_wolf_stepdown(cfg$D, B = cfg$B, alpha = cfg$alpha, seed = cfg$seed)
  pyr <- pyres$rw[[nm]]
  same <- identical(as.logical(r_rw$reject), as.logical(pyr$reject))
  n_checks <- n_checks + 1; if (same) n_pass <- n_pass + 1
  cat(sprintf("RW %-8s R=%s Py=%s -> %s\n", nm,
              paste(as.integer(r_rw$reject),collapse=","),
              paste(as.integer(pyr$reject),collapse=","),
              if (same) "OK" else "DIFFER"))
}

cat("\n===== Resampling parity summary =====\n")
cat("checks:", n_checks, "| passed:", n_pass, "| failed:", n_checks - n_pass, "\n")
if (n_pass < n_checks) quit(status = 1)
cat("Resampling parity: PASS (bootstrap MCSE envelope; set-membership agreement)\n")
quit(status = 0)
