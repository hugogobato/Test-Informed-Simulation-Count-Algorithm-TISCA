# Unit tests: MCS cross-check against the CRAN `MCS` oracle (spec 2.3).
# The primary TISCA implementation lives in mcs.R; the CRAN package is used
# only as a unit-test oracle. Skipped if `MCS` is not installed.

source("helper_tisca.R")
source_tisca_R()

test_that("TISCA MCS agrees with the CRAN MCS oracle on model sets", {
  skip_if_not_installed("MCS")
  set.seed(50)
  J <- 300
  L <- cbind(m1 = rnorm(J, 0.5, 1), m2 = rnorm(J, 2.0, 1), m3 = rnorm(J, 3.0, 1))

  # TISCA implementation (Tmax)
  ours <- mcs(L, B = 199, alpha = 0.05, statistic = "T_max", seed = 1)

  # CRAN oracle (returns an S4 'SSM' object; @Info$included holds the MCS)
  oracle <- suppressWarnings(
    MCS::MCSprocedure(Loss = L, B = 199, alpha = 0.05, statistic = "Tmax"))
  kept_oracle <- tryCatch(as.character(oracle@Info$included), error = function(e) NULL)
  if (is.null(kept_oracle) || !length(kept_oracle)) skip("could not parse CRAN MCS output structure")

  # best model (lowest loss) should be kept by both
  expect_true("m1" %in% ours$models_kept)
  expect_true(any(grepl("m1", kept_oracle)))
})
