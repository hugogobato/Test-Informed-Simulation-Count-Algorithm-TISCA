# Unit tests: R MCS against the installed CRAN `MCS` oracle.
#
# The same `MCS:::GetIndices` matrix drives both procedures.  This tests the
# statistic and elimination implementation itself, rather than comparing two
# independent Monte Carlo draws.

source("helper_tisca.R")
source_tisca_R()

test_that("TISCA MCS matches CRAN on three matrices and both statistics", {
  skip_if_not_installed("MCS")
  set.seed(50)
  cases <- list(
    cbind(m1 = rnorm(120, 0.0, 1), m2 = rnorm(120, 0.4, 1),
          m3 = rnorm(120, 0.8, 1), m4 = rnorm(120, 1.2, 1)),
    cbind(m1 = rexp(140, 1), m2 = rexp(140, 1) + 0.15,
          m3 = rexp(140, 1) + 0.30),
    cbind(m1 = rgamma(160, 2), m2 = rgamma(160, 2) + 0.25,
          m3 = rgamma(160, 2) + 0.50, m4 = rgamma(160, 2) + 0.75)
  )
  B <- 399L
  alpha <- 0.05

  for (statistic in c("Tmax", "TR")) {
    for (L in cases) {
      set.seed(100 + ncol(L) + nrow(L))
      idx <- MCS:::GetIndices(nrow(L), 0L, B)
      ours <- mcs(L, B = B, alpha = alpha, statistic = statistic,
                  bootstrap_indices = idx)
      oracle <- suppressWarnings(MCS::MCSprocedure(
        Loss = L, B = B, alpha = alpha, statistic = statistic,
        k = 0L, verbose = FALSE, seed = 100 + ncol(L) + nrow(L)))

      tab <- oracle@show
      p0 <- as.numeric(tab[, "p-Value for H_{0,M_k}"])
      names(p0) <- rownames(tab)
      pm <- as.numeric(tab[, "MCS p-Value"])
      names(pm) <- rownames(tab)
      mcse <- sqrt(p0 * (1 - p0) / B)
      expect_equal(unname(ours$p_H0[names(p0)]), unname(p0),
                   tolerance = max(4 * mcse, 1e-12))
      expect_equal(unname(ours$p_mcs[names(pm)]), unname(pm),
                   tolerance = max(4 * sqrt(pm * (1 - pm) / B), 1e-12))
      expect_setequal(ours$included, as.character(oracle@Info$included))
      expect_setequal(ours$excluded, as.character(oracle@Info$excluded))
      expect_equal(ours$elimination_order, ours$eliminated)
      expect_length(ours$elimination_order, ncol(L) - 1L)
    }
  }
})

test_that("MCS accepts all documented statistic aliases", {
  L <- cbind(m1 = 1:40, m2 = (1:40) + 0.1)
  idx <- matrix(rep(seq_len(nrow(L)), 30), nrow = nrow(L), ncol = 30)
  for (statistic in c("Tmax", "T_max", "TR", "T_R")) {
    out <- mcs(L, B = 30, statistic = statistic, bootstrap_indices = idx)
    expect_s3_class(out, "tisca_mcs")
    expect_length(out$p_H0, 2)
  }
})
