# Unit tests: estimands.R (per-replication estimands and derived RMSE).

source("helper_tisca.R")
source_tisca_R()

set.seed(101)

test_that("pehe and cate_mse differ (E[sqrt] != sqrt(E))", {
  x <- rnorm(100)      # errors
  q <- mean(x^2)
  expect_equal(pehe(0 * x, x), sqrt(q))
  expect_equal(cate_mse(0 * x, x), q)
  expect_false(isTRUE(all.equal(pehe(0 * x, x), cate_mse(0 * x, x))))
  # By Jensen, sqrt() is concave so E[sqrt(Q)] <= sqrt(E[Q]); equivalently the
  # rooted PEHE (RMS) is >= the arithmetic-mean magnitude.
  # mean(|x|) <= RMS = sqrt(mean(x^2)) is the definite AM <= RMS bound.
  expect_gte(sqrt(mean(x^2)), mean(abs(x)))
})

test_that("rmse_ate delta-method MCSE is on the rooted scale", {
  Q <- (rnorm(200, 0, 15))^2    # PEHE-scale regime, E[Q] ~ 225 -> sqrt ~ 15
  r <- rmse_ate(Q)
  expect_equal(r$rmse, sqrt(mean(Q)))
  # s_Q / (2*sqrt(J)*sqrt(mean Q))
  expect_equal(r$mcse, stats::sd(Q) / (2 * sqrt(200) * sqrt(mean(Q))),
               tolerance = 1e-12)
})

test_that("ate_cov and cate_cov give proportions in {0,1} / [0,1]", {
  expect_equal(ate_cov(2, c(0), c(3)), 1)
  expect_equal(ate_cov(5, c(0), c(3)), 0)
  tau <- seq(-1, 1, length.out = 20)
  lo  <- tau - 1; hi <- tau + 1
  expect_equal(cate_cov(tau, lo, hi), 1)
  # shift bounds to cover only half
  lo2 <- tau - 1; hi2 <- tau - 0.5
  p <- cate_cov(tau, lo2, hi2)
  expect_between(p, 0, 1)
})

test_that("cil computes the mean width", {
  expect_equal(cil(c(0, 2), c(5, 5)), 4)
})

test_that("cate_bias sign convention", {
  tau <- rnorm(50)
  expect_equal(cate_bias(tau + 0.5, tau), 0.5, tolerance = 1e-12)
})

test_that("is_perrep_loss / binary / bounded validators", {
  expect_true(is_perrep_loss(c(1, 2, 3)))
  expect_false(is_perrep_loss(c(1)))               # need >= 2
  expect_false(is_perrep_loss(c(1, NA, 3)))
  expect_true(is_binary_indicator(c(0, 1, 1, 0)))
  expect_false(is_binary_indicator(c(0, 2, 1)))
  expect_true(is_bounded_proportion(c(0.9, 0.95, 0.99)))
  expect_false(is_bounded_proportion(c(0.9, 1.2)))
})
