# Unit tests: planning.R (hypothesis modes M1-M5, variance inflation, solvers).
# All numeric expectations here come from docs/tisca_v2_spec.md rev. 2 and
# notebooks/P1_math_verification.ipynb. Running: from repo root,
#   Rscript -e 'testthat::test_local("tisca")'
# or with the runner in tests/run_tests.R.

source("helper_tisca.R")
source_tisca_R()

test_that("M1 two-sided power matches the spec's J=34 anchor", {
  # delta/sigma_D = 0.5 at 80% power requires J = 34 (spec Section 1 M1)
  p <- power("M1", 34, 0.5, 1, 0.05)
  expect_true(p > 0.80)
  p33 <- power("M1", 33, 0.5, 1, 0.05)
  expect_true(p33 <= 0.80)
})

test_that("M1 uses the noncentral t WITH lower.tail (pt semantics)", {
  # monotonic increasing in both J and |delta|
  Js <- 2:80
  ps <- vapply(Js, function(j) power("M1", j, 0.5, 1, 0.05), numeric(1))
  expect_true(all(diff(ps) >= -1e-10))
  expect_true(power("M1", 40, 0.8, 1, 0.05) > power("M1", 40, 0.5, 1, 0.05))
})

test_that("M2 directional superiority (lower is better)", {
  # H0: theta >= 0, reject iff T < t_{alpha}. Power grows as delta -> -inf.
  p0 <- power("M2", 50, -0.5, 1, 0.05)
  p1 <- power("M2", 50, -1.0, 1, 0.05)
  expect_gt(p1, p0)
  expect_between(p1, 0, 1)
})

test_that("M3 minimum-effect requires the inner boundary shift", {
  # H0: theta >= -Delta, power at delta < -Delta with ncp=sqrt(J)(delta+Delta)/s
  delta <- -1.0; Delta <- 0.5
  p_large <- power("M3", 100, delta, 1, 0.05, Delta)
  expect_between(p_large, 0, 1)
  expect_true(is.finite(p_large))
  # when delta == -Delta (boundary of H0/H1), power ~ alpha (unbounded by margin)
  p_bound <- power("M3", 200, -0.5, 1, 0.05, 0.5)
  expect_true(abs(p_bound - 0.05) < 0.02)
})

test_that("M4 non-inferiority power depends on the margin", {
  delta <- 0.2; Delta <- 0.5
  p <- power("M4", 60, delta, 1, 0.05, Delta)
  expect_between(p, 0, 1)
  # nil-power sanity when delta == Delta (null boundary)
  expect_lt(power("M4", 200, 0.5, 1, 0.05, 0.5), 0.10)
})

test_that("M5 TOST power is positive inside the margin and infeasible outside", {
  # genuine infeasibility |delta| >= Delta must be caught by the solver
  r1 <- solve_power_J("M5", 0.5, 1, 0.05, Delta = 0.5)
  expect_true(isTRUE(r1$infeasible))
  r2 <- solve_power_J("M5", 0.2, 1, 0.05, Delta = 0.5)
  expect_false(isTRUE(r2$infeasible))
  expect_true(is.finite(r2$J))
})

test_that("M5 TOST power increases with J past the small-J zero trap", {
  # Verified in spec: at Delta=0.5, sig=1, alpha=0.05, J=10 -> power 0 (m<0),
  # J=40 -> power 0.86. The plan must keep increasing J (solvers guarantee this).
  p10 <- power("M5", 10, 0.0, 1, 0.05, 0.5)
  p40 <- power("M5", 40, 0.0, 1, 0.05, 0.5)
  expect_lte(p10, 0.001)
  expect_gt(p40, 0.80)
})

test_that("sigma_ub inflation factors match the spec", {
  expect_equal(round(sigma_ub(1, 25, 0.20), 3), 1.153, tolerance = 0.002)
  expect_equal(round(sigma_ub(1, 50, 0.20), 3), 1.099, tolerance = 0.002)
  expect_equal(round(sigma_ub(1, 100, 0.20), 3), 1.067, tolerance = 0.002)
})

test_that("variance inflation is monotone (more pilot J0 -> less inflation)", {
  expect_gt(sigma_ub(1, 25), sigma_ub(1, 50))
  expect_gt(sigma_ub(1, 50), sigma_ub(1, 100))
  # higher gamma (weaker assurance) -> smaller bound
  expect_gt(sigma_ub(1, 50, 0.10), sigma_ub(1, 50, 0.50))
})

test_that("solve_power_J uses adjusted alpha and returns the smallest J", {
  # At alpha/6 and inflated sigma_D (J0=50) the spec 4.5 PEHE-Y1 cell is ~382.
  s_raw <- 2.5433      # calibrated so raw J ~= 205 at alpha = 0.05
  s_inf <- sigma_ub(s_raw, 50, 0.2)
  j_adj <- solve_power_J("M1", 0.5, s_raw, 0.05, 0.8)$J
  j_inf <- solve_power_J("M1", 0.5, s_inf, 0.05 / 6, 0.8)$J
  # the raw J at alpha = 0.05 should anchor ~205 (boundary-sensitive +-1)
  expect_equal(j_adj, 205, tolerance = 1)
  # inflated at alpha/6 should land ~382 (spec 4.5), off-by-one allowed
  expect_equal(j_inf, 382, tolerance = 2)
})

test_that("solve_power_J is monotone in sigma_D and target power", {
  j1 <- solve_power_J("M1", 0.5, 1.0, 0.05, 0.80)$J
  j2 <- solve_power_J("M1", 0.5, 1.5, 0.05, 0.80)$J
  j3 <- solve_power_J("M1", 0.5, 1.0, 0.05, 0.95)$J
  expect_gte(j2, j1)
  expect_gte(j3, j1)
})

test_that("degenerate s_D = 0 terminates without an infinite loop (P3-T7)", {
  r <- solve_power_J("M1", 0, 0, 0.05, 0.80)
  expect_true(isTRUE(r$degenerate))
  expect_true(is.finite(r$J))
})

test_that("MCSE and half-width solvers respect their monotone forms", {
  m1 <- solve_mcse_J(2, 0.1)$achieved
  expect_lte(m1, 0.1)
  h1 <- solve_halfwidth_J(2, 0.4, alpha = 0.05)$achieved
  expect_lte(h1, 0.4 + 1e-6)
  # smaller target -> larger J
  expect_gte(solve_mcse_J(2, 0.05)$J, solve_mcse_J(2, 0.1)$J)
})

test_that("combine_J applies the J_max cap and flags it", {
  c1 <- combine_J(c(100, 50), 200)
  expect_equal(c1$J_final, 100); expect_false(c1$capped)
  c2 <- combine_J(c(300, 50), 200)
  expect_equal(c2$J_final, 200); expect_true(c2$capped)
})
