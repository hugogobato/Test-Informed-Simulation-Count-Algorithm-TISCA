# Unit tests: mcs.R (Reality Check, SPA, Model Confidence Set).

source("helper_tisca.R")
source_tisca_R()

set.seed(103)

test_that("reality_check: a clear winner is significant, no-difference is not", {
  J <- 200
  L_win <- cbind(champ = rnorm(J, 1, 1), b1 = rnorm(J, 2, 1), b2 = rnorm(J, 3, 1))
  r_win <- reality_check(L_win, "champ", B = 399, seed = 1)
  expect_lt(r_win$p_value, 0.05)

  L_null <- cbind(champ = rnorm(J, 1, 1), b1 = rnorm(J, 1, 1), b2 = rnorm(J, 1, 1))
  r_null <- reality_check(L_null, "champ", B = 399, seed = 1)
  expect_gt(r_null$p_value, 0.05)
})

test_that("spa_test: same directionality as the reality check", {
  J <- 200
  L_win <- cbind(champ = rnorm(J, 1, 1), b = rnorm(J, 2, 1))
  s_win <- spa_test(L_win, "champ", B = 399, seed = 1)
  expect_lt(s_win$p_value, 0.05)
  L_null <- cbind(champ = rnorm(J, 1, 1), b = rnorm(J, 1, 1))
  s_null <- spa_test(L_null, "champ", B = 399, seed = 1)
  expect_gt(s_null$p_value, 0.05)
  # precomputed_d and (L,champ) agree
  d <- L_win[, "b"] - L_win[, "champ"]      # bench - champ = advantage
  dmat <- matrix(d, ncol = 1)
  s2 <- spa_test(dmat, precomputed_d = TRUE, B = 399, seed = 1)
  expect_equal(s2$T_max, s_win$T_max, tolerance = 1e-12)
})

test_that("MCS recovers the true best model and reports the set", {
  J <- 300
  L <- cbind(M1 = rnorm(J, 0.5, 1), M2 = rnorm(J, 2.0, 1), M3 = rnorm(J, 3.0, 1))
  m <- mcs(L, B = 199, alpha = 0.05, statistic = "T_R", seed = 1)
  expect_true("M1" %in% m$models_kept)
  expect_gte(length(m$models_kept), 1)
  # MCSE is a valid non-negative rate MCSE (or NA when nothing eliminated)
  expect_true(is.na(m$mcse) || m$mcse >= 0)
})

test_that("MCS respects the lower-is-better direction", {
  J <- 200
  # worst model has the highest mean loss -> eliminated
  L <- cbind(good = rnorm(J, 1, 1), bad = rnorm(J, 5, 1))
  m <- mcs(L, B = 99, alpha = 0.05, statistic = "T_R", seed = 1)
  expect_true("good" %in% m$models_kept)
})

test_that("mcs and reality_check never error on a single model", {
  L <- matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "solo"))
  m <- mcs(L, B = 50, alpha = 0.05, statistic = "T_R", seed = 1)
  expect_equal(m$models_kept, "solo")
  rc <- reality_check(L, "solo", B = 50, seed = 1)
  expect_true(is.finite(rc$p_value) || is.na(rc$p_value))
})
