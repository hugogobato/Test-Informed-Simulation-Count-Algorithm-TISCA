#!/usr/bin/env Rscript
# TISCA v2 — R dependency installer (P0-T1 environment specification).
#
# This is the source of truth for the R environment until the P0-T4 library
# bundle has been built. At that point `renv::snapshot()` is run over the built
# library and the resulting `renv.lock`, with fully resolved versions and
# hashes, is committed alongside this script (see notebooks/P0T4_build_rlib_bundle.ipynb).
#
# A lockfile whose package records carry no `Version` is not a valid renv
# lockfile and fails `renv::restore()`, which is why one is not committed
# before the versions are actually resolved.
#
# Usage:
#   Rscript env/install_R_dependencies.R [lib_path]
# Default lib_path is the first entry of .libPaths().

args <- commandArgs(trailingOnly = TRUE)
LIB <- if (length(args) >= 1) args[[1]] else .libPaths()[1]
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"),
        Ncpus = max(1L, parallel::detectCores()))

cat("R:", R.version.string, "\n")
cat("target library:", LIB, "\n\n")

CRAN_PKGS <- c(
  "stochtree",     # the calibrated BCF benchmark (E3)
  "dbarts",        # BART / MVBART
  "bartCause",     # treatment-effect BART
  "mvtnorm",       # DGP sampling
  "scoringRules",  # CRPS (the original authors' metric)
  "matrixStats",   # per-replication row/column statistics
  "progress",      # driver progress reporting
  "digest",        # RNG-stream digests required by docs/seed_rng_protocol.md
  "MCS"            # CRAN cross-check oracle for the v2 MCS layer
)

GITHUB_PKGS <- c(
  skewBART = "Seungha-Um/skewBART",
  mvbcf    = "Nathan-McJames/mvbcf"
)

for (p in CRAN_PKGS) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat("have    ", p, as.character(utils::packageVersion(p)), "\n")
  } else {
    cat("install ", p, "\n")
    install.packages(p, lib = LIB)
  }
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = LIB)

for (nm in names(GITHUB_PKGS)) {
  if (requireNamespace(nm, quietly = TRUE)) {
    cat("have    ", nm, as.character(utils::packageVersion(nm)), "\n")
  } else {
    cat("install ", GITHUB_PKGS[[nm]], "\n")
    remotes::install_github(GITHUB_PKGS[[nm]], lib = LIB,
                            upgrade = "never", build_vignettes = FALSE)
  }
}

# --- manifest: exactly what the manuscript has to report (IJDA #14a) ---------
ip <- utils::installed.packages(lib.loc = LIB)
keep <- intersect(rownames(ip), c(CRAN_PKGS, names(GITHUB_PKGS)))
man <- data.frame(package = keep,
                  version = unname(ip[keep, "Version"]),
                  remote_sha = vapply(keep, function(p) {
                    d <- utils::packageDescription(p, lib.loc = LIB)
                    if (!is.null(d$RemoteSha)) d$RemoteSha else NA_character_
                  }, character(1)),
                  stringsAsFactors = FALSE)
utils::write.csv(man, file.path(dirname(LIB), "R_DEPENDENCIES.csv"), row.names = FALSE)
print(man, row.names = FALSE)
cat("\nR version :", R.version.string, "\n")
cat("platform  :", R.version$platform, "\n")
cat("OS        :", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")
cat("\nWrote R_DEPENDENCIES.csv next to the library.\n")
cat("Once the P0-T4 bundle exists, run renv::snapshot() over it and commit renv.lock.\n")
