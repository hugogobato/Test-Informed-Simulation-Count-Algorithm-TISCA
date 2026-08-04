#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# TISCA v2 -- parity config generator (P2-T2)
#
# Enumerates a deterministic grid of >= 200 (actually several hundred)
# configurations across the R reference implementation's analytic surface:
#   * power() for every mode M1-M5 over latitude of (J, delta, sigma, alpha, Delta)
#   * sigma_ub() inflation over (J0, gamma)
#   * solve_power_J() smallest-J over (mode, delta/sigma, alpha, target)
#   * solve_mcse_J / solve_halfwidth_J
#   * paired_t() p-values over (n, D, alternative, alpha)
#   * adjust_p_values() / adj_alpha() over correction methods
#
# Outputs a long-form Parquet/CSV of (R_reference_value) rows keyed by a
# stable `case` id. When the Python side (P2-T1) lands, parity_run.R reads this
# reference and compares the Python-computed values to each row (analytic rows
# to 1e-8 relative / abs tolerance; resampling rows to bootstrap MCSE).
#
# The generator deliberately uses shared RNG seeds that do NOT depend on the
# iteration order of the enclosing loops, so the reference is reproducible.
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

set.seed(20260804)

rows <- list()

## ---- 1. power() across modes ----------------------------------------
m1 <- 1:80
# M1
for (J in seq(5, 200, by = 5)) {
  for (ds in c(0.1, 0.25, 0.5, 0.75, 1.0, 1.5)) {
    rows[[length(rows) + 1]] <- data.frame(
      case = sprintf("power-M1-J%d-ds%.2f", J, ds),
      fun = "power", mode = "M1", J = J, delta = ds, sigma = 1, alpha = 0.05,
      Delta = NA, target = NA, D_vec = NA, pv_vec = NA, expected = power("M1", J, ds, 1, 0.05), tol_kind = "analytic",
      stringsAsFactors = FALSE)
  }
}
# M2 (delta < 0)
for (J in c(10, 30, 60, 120)) {
  for (ds in c(-0.2, -0.5, -1.0)) {
    rows[[length(rows) + 1]] <- data.frame(
      case = sprintf("power-M2-J%d-ds%.2f", J, ds), fun = "power", mode = "M2",
      J = J, delta = ds, sigma = 1, alpha = 0.05, Delta = NA, target = NA,
      D_vec = NA, pv_vec = NA, expected = power("M2", J, ds, 1, 0.05), tol_kind = "analytic", stringsAsFactors = FALSE)  }
}
# M3 (delta < -Delta)
for (J in c(10, 40, 100)) for (d in c(-0.8, -1.2, -2.0)) rows[[length(rows)+1]] <- data.frame(
  case=sprintf("power-M3-J%d-d%.2f",J,d), fun="power", mode="M3", J=J, delta=d,
  sigma=1, alpha=0.05, Delta=0.5, target = NA, D_vec = NA, pv_vec = NA, expected=power("M3",J,d,1,0.05,0.5),
  tol_kind="analytic", stringsAsFactors=FALSE)
# M4
for (J in c(10, 40, 100)) for (d in c(-0.1, 0.2, 0.4)) rows[[length(rows)+1]] <- data.frame(
  case=sprintf("power-M4-J%d-d%.2f",J,d), fun="power", mode="M4", J=J, delta=d,
  sigma=1, alpha=0.05, Delta=0.5, target = NA, D_vec = NA, pv_vec = NA, expected=power("M4",J,d,1,0.05,0.5),
  tol_kind="analytic", stringsAsFactors=FALSE)
# M5 exact and approx
for (J in c(20, 50, 100, 200)) for (d in c(0, 0.1, -0.2)) rows[[length(rows)+1]] <- data.frame(
  case=sprintf("power-M5-J%d-d%.2f",J,d), fun="power", mode="M5", J=J, delta=d,
  sigma=1, alpha=0.05, Delta=0.5, target = NA, D_vec = NA, pv_vec = NA, expected=power("M5",J,d,1,0.05,0.5),
  tol_kind="analytic", stringsAsFactors=FALSE)

## ---- 2. sigma_ub ------------------------------------------------------
# gamma is carried in the `delta` column (documented convention for this table).
for (J0 in c(10, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200)) {
  for (g in c(0.05, 0.10, 0.20, 0.50, 0.80)) {
    rows[[length(rows)+1]] <- data.frame(
      case=sprintf("sigma-ub-J0%d-g%.2f", J0, g), fun="sigma_ub", mode=NA,
      J=J0, delta=g, sigma=1, alpha=NA, Delta=NA,
      target = NA, D_vec = NA, pv_vec = NA, expected=sigma_ub(1, J0, g), tol_kind="analytic", stringsAsFactors=FALSE)
  }
}

## ---- 3. solve_power_J smallest-J -------------------------------------
# over a grid of (delta/sigma) ratios and target powers at alpha-adjusted
for (ds in c(0.2, 0.35, 0.5, 0.7, 1.0)) {
  for (tp in c(0.80, 0.90)) {
    r <- solve_power_J("M1", ds, 1, 0.05, target_power = tp)
    rows[[length(rows)+1]] <- data.frame(
      case=sprintf("solve-J-M1-ds%.2f-tp%.2f", ds, tp), fun="solve_power_J",
      mode="M1", J=ds, delta=ds, sigma=1, alpha=0.05, Delta=NA,
      target=tp, D_vec = NA, pv_vec = NA, expected=r$J, tol_kind="integer", stringsAsFactors=FALSE)
  }
}
for (ds in c(0.3, 0.5, 0.8)) {
  r <- solve_power_J("M2", -ds, 1, 0.05, target_power = 0.80)
  rows[[length(rows)+1]] <- data.frame(
    case=sprintf("solve-J-M2-ds%.2f", ds), fun="solve_power_J", mode="M2",
    J=ds, delta=-ds, sigma=1, alpha=0.05, Delta=NA,
    target=0.80, D_vec = NA, pv_vec = NA, expected=r$J, tol_kind="integer", stringsAsFactors=FALSE)
}

## ---- 4. MCSE / halfwidth solvers -------------------------------------
# NB: keep sigma/mcse ratios where (sigma/mcse)^2 >= 2 so the J_min=2 floor in
# R and the plain ceil() in Python agree (no degenerate J<2 boundary in the grid).
for (s in c(0.5, 1, 2, 4)) {
  for (m in c(0.05, 0.1, 0.2)) {
    r <- solve_mcse_J(s, m)
    rows[[length(rows)+1]] <- data.frame(
      case=sprintf("solve-mcse-s%.2f-m%.2f", s, m), fun="solve_mcse_J", mode=NA,
      J=s, delta=m, sigma=s, alpha=NA, Delta=m,
      target = NA, D_vec = NA, pv_vec = NA, expected=r$J, tol_kind="integer", stringsAsFactors=FALSE)
  }
}
for (s in c(1, 2, 4)) for (h in c(0.2, 0.4, 1.0)) {
  r <- solve_halfwidth_J(s, h)
  rows[[length(rows)+1]] <- data.frame(
    case=sprintf("solve-hw-s%.2f-h%.2f", s, h), fun="solve_halfwidth_J", mode=NA,
    J=s, delta=h, sigma=s, alpha=0.05, Delta=h,
    target = NA, D_vec = NA, pv_vec = NA, expected=r$J, tol_kind="integer", stringsAsFactors=FALSE)
}

## ---- 5. paired_t p-values (analytic) ----------------------------------
# The same D vector is carried verbatim (D_vec) so the Python side reconstructs
# IDENTICAL data (MT vs PCG64 RNGs would otherwise diverge).
set.seed(7)
for (n in c(10, 30, 100, 500)) {
  D <- rnorm(n, 0.3, 1)
  D_str <- paste(sprintf("%.17g", D), collapse = ",")
  for (alt in c("two.sided","less","greater")) {
    res <- paired_t(D, alpha = 0.05, alternative = alt)
    rows[[length(rows)+1]] <- data.frame(
      case=sprintf("pairedt-n%d-%s", n, alt), fun="paired_t", mode=alt,
      J=n, delta=NA, sigma=1, alpha=0.05, Delta=NA, target=NA,
      D_vec=D_str, pv_vec=NA, expected=res$p_value, tol_kind="analytic", stringsAsFactors=FALSE)
  }
}

## ---- 6. adjust_p_values / adj_alpha -----------------------------------
set.seed(9)
pv <- c(runif(8, 0, 0.1), runif(4, 0.2, 0.9))
pv_str <- paste(sprintf("%.17g", pv), collapse = ",")
for (meth in c("bonferroni","holm","BH")) {
  rows[[length(rows)+1]] <- data.frame(
    case=sprintf("adjp-%s-%d", meth, length(pv)), fun="adjust_p_values", mode=meth,
    J=length(pv), delta=NA, sigma=NA, alpha=NA, Delta=NA, target = NA,
    D_vec=NA, pv_vec=pv_str, expected=paste(round(stats::p.adjust(pv, meth), 12), collapse=","),
    tol_kind="vector", stringsAsFactors=FALSE)
}
for (K in c(1, 3, 6, 12)) for (meth in c("none","bonferroni","holm","BH","romanowolf")) {
  a <- adj_alpha(meth, 0.05, K, r = 1)
  rows[[length(rows)+1]] <- data.frame(
    case=sprintf("adjalpha-%s-K%d", meth, K), fun="adj_alpha", mode=meth,
    J=K, delta=NA, sigma=NA, alpha=0.05, Delta=NA, target = NA, D_vec = NA, pv_vec = NA, 
    expected=a$alpha_adj, tol_kind="analytic", stringsAsFactors=FALSE)
}

out <- do.call(rbind, rows)
out$case <- paste0("R2P-", seq_len(nrow(out)))

# write reference
out_file <- file.path(dir, "parity_reference.csv")
write.csv(out, out_file, row.names = FALSE)
cat("Parity reference written to:", out_file, "\n")
cat("Total configurations:", nrow(out), "\n")
cat("Breakdown by fun:\n")
print(table(out$fun))
