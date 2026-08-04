#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# TISCA v2 R test runner (P2-T2).
# Usage (from repo root or from tisca/):
#   Rscript tisca/tests/run_tests.R
# or with testthat:
#   Rscript -e 'testthat::test_local("tisca", reporter = "summary")'
# ----------------------------------------------------------------------------

argv <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1])
if (!length(this_file)) this_file <- "tisca/tests/run_tests.R"

repo_root <- normalizePath(dirname(this_file))
while (!file.exists(file.path(repo_root, "tisca", "R")) && nzchar(basename(repo_root))) {
  repo_root <- dirname(repo_root)
}
test_dir <- file.path(repo_root, "tisca", "tests")
setwd(test_dir)

suppressPackageStartupMessages(library(testthat))

test_files <- list.files(test_dir, pattern = "^test_.*\\.R$", full.names = TRUE)
if (!length(test_files)) stop("no test files found under ", test_dir)

cat("TISCA v2 R unit tests in:", test_dir, "\n")
res <- testthat::test_dir(test_dir, reporter = "summary", stop_on_failure = FALSE)
if (length(res) && any(as.logical(res$failed))) {
  quit(status = 1)
}
cat("All TISCA v2 R unit tests passed.\n")
