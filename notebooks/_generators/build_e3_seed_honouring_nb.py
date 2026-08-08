#!/usr/bin/env python3
"""Generate ``E3_seed_honouring.ipynb`` (CALIBRATION.md §3.6, the last open box).

``E3_seed_verification.ipynb`` closed P3-T5 and P0-T2 §3.5 but explicitly did
**not** close §3.6, and could not: ``run_cell.R`` positions a fixed L'Ecuyer
substream in the *global* RNG state immediately before every fit, so re-running a
seed reproduces a model's columns whether the model honoured its explicit seed
argument or ignored it and drew from the global stream. Both hypotheses predict
identical output, so that comparison has no power to separate them.

This notebook runs the discriminating comparison instead: **same model, same data,
same explicit seed, two different global RNG states.**

Two things narrow the job a long way, and both were settled by reading source
rather than by spending a session:

* ``fast_bart()`` takes no seed argument. Its signature at ``MVBCF_Code.cpp:659``
  has 18 parameters, none of them a seed, and every draw inside the file goes
  through R's own generator (``R::runif``, plus ``rmvnorm``/``riwish``/``sample``
  from RcppDist and Rcpp sugar), with the ``Rcpp::export`` ``RNGScope`` reading and
  writing ``.Random.seed``. There is no independent generator in the file: no
  ``mt19937``, no ``rand()``, no ``arma::randu``/``randn``. ``MultiskewBART`` is
  likewise called without a seed. For those two, §3.6 is inapplicable rather than
  open, and §3.3 stream positioning (already validated by the §3.5 four-way
  identity test) is the whole mechanism.
* That leaves two libraries that do take a seed: ``stochtree::bcf``
  (``general_params$random_seed``) and ``dbarts``, used directly for the propensity
  fit and underneath ``bartCause::bartc``.

**No E3 replication is re-run here.** The notebook reads none of the 4,000 stored
rows; it generates one DGP1 ``n = 100`` dataset and issues nine fits. Nothing about
the shipped data is contingent on the outcome, because exact reproduction is already
an observed fact across ten hosts, two worker counts and two shard offsets. What the
outcome decides is only what the Data Availability statement may attribute that
reproducibility *to*.

**Read the verdict carefully; it is asymmetric.** A PASS (output invariant to the
global state at a fixed explicit seed) is conclusive: the seed argument alone
determines the fit. A FAIL does **not** mean the seed argument is ignored. It means
the fit also consumes the global stream somewhere, which is a normal thing for an R
model wrapper to do and is exactly the case §3.3 positioning exists to cover. The
positive control in each block guards against the opposite error: if two *different*
explicit seeds also produce identical output, the fit is seed-insensitive and the
whole comparison is vacuous, so that is reported as INCONCLUSIVE rather than PASS.

Runtime is about 20 minutes, most of it environment setup: one Colab session.

Regenerate with::

    python notebooks/_generators/build_e3_seed_honouring_nb.py \\
        --bundle-folder-url '<drive folder>' --bundle-sha256 '<sha>'
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

from _nbcommon import SETUP, md, code, write

_SPEC = importlib.util.spec_from_file_location(
    "_e3gen", Path(__file__).with_name("build_e3_notebooks.py"))
_E3 = importlib.util.module_from_spec(_SPEC)
sys.modules["_e3gen"] = _E3
_SPEC.loader.exec_module(_E3)


# The probe driver. Deliberately standalone: it does not import, source or modify
# run_cell.R, so nothing here can perturb the driver that produced the confirmatory
# shards. The DGP block and every model call are copied verbatim from run_cell.R
# (generate_data/generate_test at :139-199, the propensity fit at :341-344, the BCF
# specification at :425-449, the bartc call at :502-505) so the answer is about this
# study's fits, not about the libraries in the abstract.
PROBE_R = r'''
PROBE_R_SOURCE = r"""
suppressWarnings(suppressMessages({
  library(dbarts); library(bartCause); library(stochtree); library(mvtnorm)
}))

RNGkind("L'Ecuyer-CMRG")

# ---- settings, verbatim from run_cell.R:217-221 ---------------------------- #
n_train <- 100L
n_test  <- 1000L
n_tree_mu <- 50L; n_tree_tau <- 20L; n_trees_total <- n_tree_mu + n_tree_tau
n_iter <- 1000L; n_burn <- 500L; n_mcmc <- n_iter - n_burn
nthread_global <- 1L

# ---- DGP1, verbatim from run_cell.R:139-199 -------------------------------- #
generate_data <- function(n) {
  X1 <- runif(n); X2 <- runif(n); X3 <- runif(n); X4 <- runif(n); X5 <- runif(n)
  X6 <- rbinom(n, 1, 0.5); X7 <- rbinom(n, 1, 0.5); X8 <- rbinom(n, 1, 0.5)
  X9 <- sample(c(0, 1, 2, 3, 4), n, replace = TRUE)
  X10 <- sample(c(0, 1, 2, 3, 4), n, replace = TRUE)
  X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  Mu1 <- (11*sin(pi*X1*X2)+18*(X3-0.5)^2+10*X4+12*X6+X9)*10+300
  Mu2 <- (9*sin(pi*X1*X2)+22*(X3-0.5)^2+14*X4+8*X6+X9)*10+300
  Tau1 <- (2*X4+2*X5)*10
  Tau2 <- (1*X4+3*X5)*10
  Z <- rbinom(n, 1, X4)
  Y <- cbind(Mu1 + Z*Tau1, Mu2 + Z*Tau2) +
    mvtnorm::rmvnorm(n, c(0, 0), matrix(c(50^2, 0, 0, 50^2), nrow = 2, byrow = TRUE))
  list(X = X, Y = Y, Z = Z)
}

paper_lambda <- function(y, z, x_control, pihat, nu = 3, sigq = 0.9) {
  ystd <- (y - mean(y)) / sd(y)
  x_c <- cbind(x_control, pihat)
  sighat <- summary(lm(ystd ~ z + as.matrix(x_c)))$sigma
  (sighat^2 * qchisq(1 - sigq, nu)) / nu
}

set_rand <- function(x) assign(".Random.seed", x, envir = .GlobalEnv)

# Two global RNG states that are certainly different: independent L'Ecuyer
# substreams, which is what run_cell.R hands to consecutive replications. If a fit
# is invariant across these two, it is invariant across the whole design.
set.seed(11L); STATE_A <- .Random.seed
STATE_B <- parallel::nextRNGStream(parallel::nextRNGStream(STATE_A))
stopifnot(!identical(STATE_A, STATE_B))

SEED_1 <- 424242L   # the explicit seed held fixed across STATE_A / STATE_B
SEED_2 <- 989898L   # a different explicit seed, for the positive control

# ---- data: generated ONCE, from a state neither probe state uses ----------- #
set.seed(20260808L)
d <- generate_data(n_train)
X <- d$X; Y <- d$Y; Z <- d$Z
td <- generate_data(n_test)
X_test <- td$X; Z_test <- td$Z

# ---- fits, each a closure of (explicit seed, global state) ------------------ #
# Each returns a numeric vector: the object run_cell.R actually turns into the
# recorded metrics, flattened. Comparing the full posterior object rather than a
# scalar summary means a difference cannot hide behind an averaging step.

fit_dbarts <- function(seed_val, state) {
  set_rand(state)
  m <- dbarts::bart(x.train = X, y.train = Z, x.test = X_test, k = 3,
                    nthread = nthread_global, seed = seed_val, verbose = FALSE)
  as.numeric(m$yhat.test)
}

fit_bcf <- function(seed_val, state) {
  set_rand(state)
  # The propensity is fixed input here, not part of what is under test: compute it
  # once, deterministically, so a difference cannot be inherited from upstream.
  p <- P_FIXED; p_test <- P_TEST_FIXED
  lambda <- paper_lambda(Y[, 1], Z, X, p)
  m <- stochtree::bcf(
    X_train = X, Z_train = Z, y_train = Y[, 1],
    propensity_train = p, X_test = X_test, Z_test = Z_test,
    propensity_test = p_test, num_gfr = 0,
    num_burnin = n_burn, num_mcmc = n_mcmc,
    general_params = list(
      standardize = TRUE, sample_sigma2_global = TRUE,
      sigma2_global_shape = 3 / 2,
      sigma2_global_scale = if (is.finite(lambda)) 3 * lambda / 2 else 1,
      propensity_covariate = "prognostic", adaptive_coding = FALSE,
      num_chains = 1, num_threads = nthread_global, random_seed = seed_val),
    prognostic_forest_params = list(
      num_trees = n_tree_mu, alpha = 0.95, beta = 2,
      min_samples_leaf = 1, max_depth = -1,
      sample_sigma2_leaf = FALSE,
      sigma2_leaf_init = 1^2 / n_tree_mu),
    treatment_effect_forest_params = list(
      num_trees = n_tree_tau, alpha = 0.25, beta = 3,
      min_samples_leaf = 1, max_depth = -1,
      sample_sigma2_leaf = FALSE,
      sigma2_leaf_init = 0.375^2 / n_tree_tau))
  as.numeric(m$tau_hat_test)
}

fit_bartc <- function(seed_val, state) {
  set_rand(state)
  Xa <- cbind(X, P_FIXED)
  colnames(Xa) <- paste0("V", seq_len(ncol(Xa)))
  Xt <- cbind(X_test, P_TEST_FIXED)
  colnames(Xt) <- colnames(Xa)
  m <- bartCause::bartc(Y[, 1], Z, Xa, p.scoreAsCovariate = FALSE, n.chains = 1,
                        n.threads = 1, keepTrees = TRUE, n.trees = n_trees_total,
                        seed = seed_val, verbose = FALSE)
  as.numeric(predict(m, Xt, type = "icate"))
}

# Propensity, computed once from a third state so every later fit sees the same
# numbers. This is input to bcf and bartc, so it must not vary between arms.
set_rand(STATE_A)
p_mod <- dbarts::bart(x.train = X, y.train = Z, x.test = X_test, k = 3,
                      nthread = nthread_global, seed = 7777L, verbose = FALSE)
P_FIXED <- colMeans(pnorm(p_mod$yhat.train))
P_TEST_FIXED <- colMeans(pnorm(p_mod$yhat.test))

MODELS <- list(
  "dbarts::bart (propensity fit)" = fit_dbarts,
  "stochtree::bcf"                = fit_bcf,
  "bartCause::bartc (dbarts)"     = fit_bartc)

rows <- list()
for (nm in names(MODELS)) {
  f <- MODELS[[nm]]
  cat("fitting:", nm, "\n"); flush.console()
  t0 <- proc.time()[["elapsed"]]
  a1 <- f(SEED_1, STATE_A)          # reference
  a2 <- f(SEED_1, STATE_B)          # same explicit seed, DIFFERENT global state
  ctrl <- f(SEED_2, STATE_A)        # different explicit seed, same global state
  stopifnot(length(a1) == length(a2), length(a1) == length(ctrl))
  d_state <- max(abs(a1 - a2))
  d_seed  <- max(abs(a1 - ctrl))
  rows[[length(rows) + 1L]] <- data.frame(
    model = nm,
    n_values = length(a1),
    max_abs_diff_same_seed_other_state = d_state,
    max_abs_diff_other_seed_same_state = d_seed,
    invariant_to_global_state = identical(a1, a2),
    sensitive_to_explicit_seed = !identical(a1, ctrl),
    fit_seconds = round(proc.time()[["elapsed"]] - t0, 1),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
out$verdict <- ifelse(!out$sensitive_to_explicit_seed, "INCONCLUSIVE",
               ifelse(out$invariant_to_global_state, "HONOURS_SEED",
                      "ALSO_USES_GLOBAL_STREAM"))
write.csv(out, OUT_PATH, row.names = FALSE)
print(out)
cat("PROBE_OK\n")
"""
'''

WRITE_AND_RUN = r'''
import os
import subprocess
import time

import pandas as pd

os.makedirs("/content/e3", exist_ok=True)
OUT_CSV = "/content/e3/seed_honouring.csv"
script = ('OUT_PATH <- "%s"\n' % OUT_CSV) + PROBE_R_SOURCE
with open("/content/e3/seed_probe.R", "w") as f:
    f.write(script)

env = dict(os.environ)
env["R_LIBS"] = LIBDIR + ":" + env.get("R_LIBS", "")

t0 = time.time()
p = subprocess.run(["Rscript", "/content/e3/seed_probe.R"],
                   capture_output=True, text=True,
                   encoding="utf-8", errors="replace", env=env)
print(p.stdout[-6000:])
print(f"probe wall-clock: {time.time() - t0:.0f}s")
if p.returncode != 0 or "PROBE_OK" not in p.stdout:
    print(p.stderr[-6000:])
    raise RuntimeError("seed probe failed; do not tick the 3.6 box")

probe = pd.read_csv(OUT_CSV)
print()
print(probe.to_string(index=False))
'''

VERDICT = r'''
# --------------------------------------------------------------------------- #
# Verdict                                                                       #
# --------------------------------------------------------------------------- #
# Three outcomes, and only one of them is a defect-shaped finding:
#
#   HONOURS_SEED             the fit is bit-identical at a fixed explicit seed across
#                            two different global RNG states, and differs when the
#                            seed changes. The seed argument alone determines the fit.
#   ALSO_USES_GLOBAL_STREAM  the fit responds to its seed BUT also consumes the global
#                            stream. Normal for an R model wrapper, and precisely the
#                            case docs/seed_rng_protocol.md 3.3 positioning covers.
#                            run_cell.R sets that stream deterministically per fit, so
#                            reproducibility is unaffected -- the ATTRIBUTION changes,
#                            not the guarantee.
#   INCONCLUSIVE             two different explicit seeds gave identical output, so
#                            the fit is not seed-sensitive at all and this comparison
#                            has no power. Investigate before reporting anything.
#
# fast_bart() (MVBCF) and MultiskewBART are absent by construction: neither exposes a
# seed argument, so there is nothing for them to honour. They are governed entirely by
# global stream positioning, which the 3.5 four-way identity test already validated.

probe["attribution"] = probe["verdict"].map({
    "HONOURS_SEED": "explicit seed argument",
    "ALSO_USES_GLOBAL_STREAM": "global stream positioning (run_cell.R 3.3)",
    "INCONCLUSIVE": "UNDETERMINED -- comparison had no power",
})
inconclusive = probe[probe["verdict"] == "INCONCLUSIVE"]

os.makedirs(os.path.join(RESULTS, "E3"), exist_ok=True)
dest = os.path.join(RESULTS, "E3", "seed_honouring.csv")
probe.to_csv(dest, index=False)
download(dest)

print(probe[["model", "verdict", "attribution",
             "max_abs_diff_same_seed_other_state",
             "max_abs_diff_other_seed_same_state"]].to_string(index=False))
print()
if len(inconclusive):
    print("*** INCONCLUSIVE rows above: do NOT tick the 3.6 box for them ***")
else:
    print("Every model under test responded to its explicit seed, so the "
          "comparison had power and the verdicts stand.")

honours = probe.loc[probe["verdict"] == "HONOURS_SEED", "model"].tolist()
global_too = probe.loc[probe["verdict"] == "ALSO_USES_GLOBAL_STREAM", "model"].tolist()

print()
print("--- paste into experiments/E3_mvbcf_casestudy/CALIBRATION.md ---")
print(f"""
### P0-T2 3.6, per-model seed honouring, {pd.Timestamp.today().date()}

Run by `notebooks/E3_seed_honouring.ipynb`; table in
`results/E3/seed_honouring.csv`. One DGP1 n = 100 dataset, nine fits, no stored
replication read or re-run.

- [{'x' if not len(inconclusive) else ' '}] Discriminating comparison performed:
      same model, same data, same explicit seed, two different global RNG states,
      with a different-seed positive control confirming the comparison had power.
- Determined by the explicit seed argument alone: {honours or 'none'}.
- Responds to the seed but also consumes the global stream: {global_too or 'none'}.
      Not a defect: `run_cell.R` positions a fixed L'Ecuyer substream before every
      fit, which is what 3.3 requires and what the 3.5 four-way identity test
      already validated on the real driver.
- Out of scope by construction: `fast_bart()` (MVBCF) and `MultiskewBART` expose no
      seed argument, so `model_seed_mvbcf` and `model_seed_mvbart` are recorded
      labels rather than arguments passed. Their reproducibility rests entirely on
      3.3 positioning.

Data Availability wording this licenses: reproducibility of the shipped rows is
established empirically (three seeds, ten hosts, 159 columns, zero differences) and
is attributed to deterministic global-stream positioning{', with the seed argument sufficient on its own for ' + ', '.join(honours) if honours else ''}.
""")
'''


def build(bundle_source, bundle_sha):
    cells = []
    md(cells, """
        # E3 per-model seed honouring (CALIBRATION.md §3.6)

        `E3_seed_verification.ipynb` closed P3-T5 and P0-T2 §3.5. It could not close
        §3.6, and said so: `run_cell.R` positions a fixed L'Ecuyer substream in the
        **global** RNG state immediately before every fit, so re-running a seed
        reproduces a model's columns whether that model honoured its explicit seed
        argument or ignored it and drew from the global stream. Both hypotheses
        predict identical output, so the comparison has no power to separate them.

        This notebook runs the comparison that does: **same model, same data, same
        explicit seed, two different global RNG states.**

        **No E3 replication is re-run.** Not one of the 4,000 stored rows is read.
        The notebook generates one DGP1 `n = 100` dataset and issues nine fits, and
        nothing about the shipped data depends on the answer: exact reproduction is
        already an observed fact across ten hosts, two worker counts and two shard
        offsets. What the answer decides is what the Data Availability statement may
        attribute that reproducibility *to*.

        Two of the four models named in §3.6 are out of scope by construction, and
        source inspection settles that without spending a session. `fast_bart()`
        takes no seed argument: its signature at `MVBCF_Code.cpp:659` has 18
        parameters, none of them a seed, and every draw in the file goes through R's
        own generator (`R::runif`, plus `rmvnorm`, `riwish` and `sample` from
        RcppDist and Rcpp sugar), with the `Rcpp::export` `RNGScope` reading and
        writing `.Random.seed`. There is no independent generator anywhere in it: no
        `mt19937`, no `rand()`, no `arma::randu`/`randn`. `MultiskewBART` is likewise
        called without a seed. Both are governed entirely by §3.3 stream positioning,
        which the §3.5 four-way identity test already validated.

        That leaves `stochtree::bcf` and `dbarts` (used directly for the propensity
        fit and underneath `bartCause::bartc`).

        **The verdict is asymmetric, so read it carefully.** A PASS is conclusive:
        the explicit seed alone determines the fit. A FAIL does *not* mean the seed
        argument is ignored, it means the fit also consumes the global stream, which
        is ordinary behaviour for an R model wrapper and exactly the case §3.3
        positioning exists to cover. Each block also carries a positive control at a
        *different* explicit seed: if that comes back identical too, the fit is not
        seed-sensitive at all and the result is reported as INCONCLUSIVE rather than
        as a pass.

        Runtime is about 20 minutes, most of it environment setup.
        """)
    # SETUP first (REPO_ROOT, RESULTS, FIGURES, download()); then the E3 shared cells
    # for hardware report, R install and bundle restore. The fast_bart compile cell is
    # kept even though MVBCF is out of scope, because it is what verifies the restored
    # bundle is the same one the confirmatory shards used.
    code(cells, SETUP)
    row = {"round": "seed-honouring", "dgp": 1, "n": 100, "mode": "confirmatory",
           "seed_start": 0, "seed_end": 0}
    for c in _E3._shared_cells(bundle_source, bundle_sha, row, None)[1:]:
        cells.append(c)
    md(cells, """
        ## The probe driver

        Standalone on purpose: it does not import, source or modify `run_cell.R`, so
        nothing here can perturb the driver that produced the confirmatory shards.
        The DGP block and every model call are copied verbatim from it, so the answer
        is about this study's fits rather than about the libraries in the abstract.
        """)
    code(cells, PROBE_R)
    code(cells, WRITE_AND_RUN)
    md(cells, "## Verdict, and the block to paste into CALIBRATION.md")
    code(cells, VERDICT)
    write("E3_seed_honouring.ipynb", cells)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    _E3._bundle_args(ap)
    args = ap.parse_args()
    source, sha = _E3._validated_bundle(args)
    build(source, sha)


if __name__ == "__main__":
    main()
