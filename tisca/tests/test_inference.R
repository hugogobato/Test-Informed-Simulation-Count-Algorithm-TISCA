# Unit tests: inferface.R (listwise contrasts, paired t, paired bootstrap, McNemar).

source("helper_tisca.R")
source_tisca_R()

set.seed(100)

test_that("contrast_from_columns drops NAs LISTWISE, never per-column (C1/2b)", {
  A <- c(1, 2, NA, 4, 5)
  B <- c(1, NA, 3, 4, 5)
  cc <- contrast_from_columns(A, B)
  expect_equal(cc$n_kept, 3)      # rep 2 and 3 dropped (A3 or B2 missing)
  expect_equal(cc$dropped, 2)
  expect_equal(cc$na_A, 1)
  expect_equal(cc$na_B, 1)
  expect_equal(cc$D, c(0, 0, 0))   # only rows with BOTH present
  expect_true(cc$paired)
})

test_that("contrast_from_columns enforces parallel lengths", {
  expect_error(contrast_from_columns(c(1, 2), c(1, 2, 3)))
})

test_that("paired_t returns estimate before p-value and correct t", {
  D <- rnorm(60, 0.5, 1)
  res <- paired_t(D, alpha = 0.05, alternative = "two.sided")
  expect_equal(res$estimate, mean(D))
  expect_equal(res$mcse, stats::sd(D) / sqrt(60))
  # manual p
  t_manual <- mean(D) / (stats::sd(D) / sqrt(60))
  p_manual <- 2 * stats::pt(-abs(t_manual), 59)
  expect_equal(res$p_value, p_manual)
  expect_equal(names(res)[1], "estimate")   # estimate-first report
})

test_that("paired_t matches t.test on D (paired semantics)", {
  A <- rnorm(40, 3, 1); B <- rnorm(40, 3.3, 1)
  res <- paired_t(A - B, alternative = "two.sided")
  tt <- stats::t.test(A, B, paired = TRUE)$p.value
  expect_equal(res$p_value, tt, tolerance = 1e-12)
})

test_that("studentized paired bootstrap is finite and centred correctly", {
  D <- rnorm(80, 0.4, 1)
  res <- studentized_paired_bootstrap(D, B = 299, alpha = 0.05, seed = 1)
  expect_true(is.finite(res$p_value))
  expect_equal(res$estimate, mean(D))
  expect_true(res$ci_lower < res$estimate & res$ci_upper > res$estimate)
})

test_that("studentized bootstrap rejects a clear non-null", {
  D <- rnorm(100, 2, 1)
  res <- studentized_paired_bootstrap(D, B = 299, alpha = 0.05, seed = 2)
  expect_lt(res$p_value, 0.01)
})

test_that("Mcnemar exact on paired binary indicators", {
  # no discordance -> p = 1
  r1 <- mcnemar_exact(rep(1, 20), rep(1, 20))
  expect_equal(r1$p_value, 1)
  # strong one-directional discordance
  xA <- c(rep(1, 15), rep(0, 5)); xB <- rep(0, 20)
  r2 <- mcnemar_exact(xA, xB, alternative = "two.sided")
  expect_true(r2$b > 0)
  expect_lt(r2$p_value, 0.05)
  # symmetric discordance -> high p
  xA2 <- c(rep(1, 8), rep(0, 12)); xB2 <- c(rep(0, 12), rep(1, 8))
  r3 <- mcnemar_exact(xA2, xB2)
  expect_gt(r3$p_value, 0.05)
})

test_that("interval score penalises non-coverage and width", {
  # point far outside -> big score
  s_out <- mean(interval_score(c(5), c(0), c(1), c = 0.05))
  # inside -> just the width
  s_in <- mean(interval_score(c(0.5), c(0), c(1), c = 0.05))
  expect_gt(s_out, s_in)
  expect_equal(s_in, 1)   # width = 1, no penalty
})
