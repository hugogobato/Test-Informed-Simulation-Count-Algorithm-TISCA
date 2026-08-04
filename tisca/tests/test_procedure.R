# Unit tests: procedure.R (two-stage default, adaptive mode) and validate.R.

source("helper_tisca.R")
source_tisca_R()

set.seed(104)

# deterministic controllable simulation: A is better than B by `effect`
make_sim <- function(effect = 0.5, sd = 1) {
  function(seed) {
    set.seed(seed)
    data.frame(pehe_A = rnorm(1, 3 + effect, sd), pehe_B = rnorm(1, 3, sd))
  }
}

test_that("two-stage returns a valid plan and correct inference", {
  sim <- make_sim(effect = -0.5)
  primary <- list(list(metric_a = "pehe_A", metric_b = "pehe_B", mode = "M1",
                       delta = -0.5, target_type = "power", target_power = 0.80,
                       J0 = 40, gamma = 0.20, metric_type = "loss"))
  res <- run_two_stage(sim, primary, J0 = 40, correction = "none", alpha = 0.05,
                       seed = 3, verbose = FALSE)
  expect_true(res$J_final >= 2)
  expect_false(res$capped)
  expect_true(isTRUE(res$pilot_discarded))
  # confirmatory rows must equal J_final
  expect_equal(nrow(res$confirm_frame), res$J_final)
  # the estimate should reflect A-better-than-B (negative theta in our loss framing)
  expect_lt(res$final[[1]]$estimate, 0)
})

test_that("two-stage: pilot count and confirmatory count are auditable", {
  sim <- make_sim()
  primary <- list(list(metric_a = "pehe_A", metric_b = "pehe_B", mode = "M2",
                       delta = -0.5, target_type = "power", target_power = 0.80,
                       J0 = 30, gamma = 0.20, metric_type = "loss"))
  res <- run_two_stage(sim, primary, J0 = 30, correction = "none", alpha = 0.05)
  expect_equal(nrow(res$pilot_frame), 30)
  expect_equal(res$p_dropped, 0)             # no NAs here
})

test_that("two-stage with a J_max cap flags capped=True", {
  sim <- make_sim(effect = -0.5)
  primary <- list(list(metric_a = "pehe_A", metric_b = "pehe_B", mode = "M1",
                       delta = -0.5, target_type = "power", target_power = 0.80,
                       J0 = 30, gamma = 0.20, metric_type = "loss", J_max = 15))
  res <- run_two_stage(sim, primary, J0 = 30, J_max = 15, correction = "none")
  expect_equal(res$J_final, 15)
  expect_true(res$capped)
})

test_that("adaptive mode refuses to run without accepting the caveat", {
  sim <- make_sim()
  primary <- list(list(metric_a = "pehe_A", metric_b = "pehe_B", mode = "M1",
                       delta = -0.5, target_type = "power", target_power = 0.80,
                       J0 = 30, gamma = 0.20, metric_type = "loss"))
  expect_error(run_adaptive(sim, primary, J0 = 30, accept_caveat = FALSE))
  # with caveat accepted, runs to completion
  res <- run_adaptive(sim, primary, J0 = 30, accept_caveat = TRUE,
                      nmax_looks = 2, J_max = 200, correction = "none")
  expect_false(isTRUE(res$pilot_discarded))     # pilot reused
  expect_match(res$stage_description, "REUSED")
})

test_that("validate_contrast rejection rules fire (spec Section 8)", {
  # Rule 1: aggregate name
  r1 <- validate_contrast("pehe", "pehe_mean", "loss", lower_is_better = TRUE)
  expect_false(r1$pass)
  expect_match(paste(r1$problems, collapse = " "), "aggregate")
  # Rule 2: per-replication RMSE_ATE
  r2 <- validate_contrast("a", "b", "rmse_ate")
  expect_false(r2$pass)
  expect_match(paste(r2$problems, collapse = " "), "squared-ATE-error")
  # Rule 3: coverage as lower-is-better
  r3 <- validate_contrast("covA", "covB", "bounded_prop", lower_is_better = TRUE)
  expect_false(r3$pass)
  expect_match(paste(r3$problems, collapse = " "), "nominal target")
  # valid request passes
  rv <- validate_contrast("pehe_A", "pehe_B", "loss", lower_is_better = TRUE)
  expect_true(rv$pass)
})

test_that("validate_family enforces family size", {
  p1 <- list(list(metric_a = "pehe_A", metric_b = "pehe_B", mode = "M1"))
  v <- validate_family(p1, K_expected = 2)
  expect_false(v$pass)     # family size mismatch
  expect_equal(v$n_contrasts, 1)
})
