# Unit tests: multiplicity.R (alpha planning levels, p.adjust, Romano-Wolf, family power).

source("helper_tisca.R")
source_tisca_R()

set.seed(102)

test_that("adj_alpha maps methods to per-contrast levels", {
  expect_equal(adj_alpha("none", 0.05, 6)$alpha_adj, 0.05)
  expect_equal(adj_alpha("bonferroni", 0.05, 6)$alpha_adj, 0.05 / 6)
  expect_equal(adj_alpha("holm", 0.05, 6)$alpha_adj, 0.05 / 6)
  # Romano-Wolf is run at the family alpha (planned by simulation from Sigma_D)
  expect_equal(adj_alpha("romanowolf", 0.05, 6)$alpha_adj, 0.05)
  # BH with r=1 -> Bonferroni; with r=2 -> alpha*2/K
  expect_equal(adj_alpha("BH", 0.05, 6, r = 1)$alpha_adj, 0.05 / 6)
  expect_equal(adj_alpha("BH", 0.05, 6, r = 3)$alpha_adj, 0.05 * 3 / 6)
})

test_that("adjust_p_values matches stats::p.adjust", {
  p <- c(0.01, 0.04, 0.2, 0.05)
  expect_equal(adjust_p_values(p, "bonferroni"), stats::p.adjust(p, "bonferroni"))
  expect_equal(adjust_p_values(p, "holm"), stats::p.adjust(p, "holm"))
  expect_equal(adjust_p_values(p, "BH"), stats::p.adjust(p, "BH"))
  expect_equal(adjust_p_values(p, "none"), p)
})

test_that("romano_wolf_stepdown rejects true positives and is FWE-bounded under null", {
  # Two true contrasts (mean .5), one null (mean 0)
  D <- cbind(p1 = rnorm(300, 0.5, 1), p2 = rnorm(300, 0.6, 1), null = rnorm(300, 0, 1))
  rw <- romano_wolf_stepdown(D, alpha = 0.05, B = 299, seed = 1)
  expect_true(rw$reject[1] && rw$reject[2])
  # Under the global null, FWE <= alpha with high probability
  nsim <- 40
  fwe <- mean(replicate(nsim, {
    D0 <- matrix(rnorm(300 * 3), 300, 3)
    any(romano_wolf_stepdown(D0, alpha = 0.05, B = 99, seed = NULL)$reject)
  }))
  expect_lte(fwe, 0.10)   # allowance for MC
})

test_that("family_power handles conjunctive / disjunctive / named", {
  # rows: T,T,F / T,F,T / F,T,F / F,F,F
  rej <- cbind(c(T, T, F, F), c(T, F, T, F), c(F, F, F, F))
  expect_equal(family_power(rej, "conjunctive")$power, 0)
  expect_equal(family_power(rej, "disjunctive")$power, 0.75)  # rows 1,2,3 have >=1 TRUE
  expect_equal(family_power(rej, "named", named_k = 1:2)$power, 0.25)
})

test_that("plan_rw_power returns a monotone power grid and a J", {
  # K=3 marginally identical contrasts with delta=0.5, off-diag corr 0.3
  Delta_vec <- c(0.5, 0.5, 0.5)
  S <- matrix(0.3, 3, 3); diag(S) <- 1
  res <- plan_rw_power(Delta_vec, S, target_power = 0.80, alpha = 0.05,
                       criterion = "disjunctive", R_plan = 40, B = 19,
                       J_grid = c(50, 100, 150), seed = 1)
  expect_true(all(diff(res$power_grid) >= -1e-12))
  expect_true(is.integer(res$J) || is.numeric(res$J))
})
