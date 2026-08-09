# Numerical validation of the v2 interval-score helpers, extracted verbatim
# from run_cell_v2.R, against (i) a naive scalar loop and (ii) the v1 helpers.
src <- readLines("experiments/E3_mvbcf_casestudy/run_cell_v2.R")
lo <- grep("^in_cred <- function", src)
hi <- grep('^E3_INTERVAL_LEVELS', src)
eval(parse(text = paste(src[lo:hi], collapse = "\n")))

naive_is <- function(draws_by_obs, truth, level) {
  a <- 1 - level
  s <- numeric(length(truth))
  for (i in seq_along(truth)) {
    q <- quantile(draws_by_obs[[i]], c(a/2, 1 - a/2), names = FALSE)
    pen <- 0
    if (truth[i] < q[1]) pen <- (2/a) * (q[1] - truth[i])
    if (truth[i] > q[2]) pen <- (2/a) * (truth[i] - q[2])
    s[i] <- (q[2] - q[1]) + pen
  }
  mean(s)
}

set.seed(11)
n_obs <- 137; n_draw <- 500
fails <- 0
report <- function(label, ok, detail = "") {
  cat(sprintf("%-58s %s %s\n", label, if (ok) "PASS" else "FAIL", detail))
  if (!ok) fails <<- fails + 1
}

for (level in c(0.5, 0.95)) {
  # Rows = observations (stochtree / fast_bart / MultiskewBART convention).
  M <- matrix(rnorm(n_obs * n_draw), nrow = n_obs)
  truth <- rnorm(n_obs, sd = 1.4)              # deliberately overdispersed:
  truth[1:20] <- truth[1:20] + 6               # forces genuine misses
  rows_list <- lapply(seq_len(n_obs), function(i) M[i, ])
  got <- interval_score_mat(M, 1L, truth, level)
  want <- naive_is(rows_list, truth, level)
  report(sprintf("margin 1 score vs naive loop (level %.2f)", level),
         abs(got$score - want) < 1e-12, sprintf("|d|=%.3g", abs(got$score - want)))

  # Columns = observations (bartCause icate convention). Same numbers.
  gotT <- interval_score_mat(t(M), 2L, truth, level)
  report(sprintf("margin 2 on transpose reproduces margin 1 (level %.2f)", level),
         abs(gotT$score - got$score) < 1e-12)

  # Width component must equal the v1 cred_width() column exactly.
  v1_width <- mean(apply(M, 1, cred_width, level))
  report(sprintf("width component == v1 cred_width mean (level %.2f)", level),
         abs(got$width - v1_width) < 1e-12, sprintf("|d|=%.3g", abs(got$width - v1_width)))

  # Decomposition identity.
  report(sprintf("score == width + penalty (level %.2f)", level),
         abs(got$score - (got$width + got$penalty)) < 1e-12)

  # Penalty must be zero exactly when every truth is inside its interval, and
  # then the score collapses to the width. Uses v1 in_cred() as the oracle.
  inside <- mean(diag(apply(M, 1, in_cred, truth, level)))
  clean <- interval_score_mat(M, 1L, apply(M, 1, median), level)
  report(sprintf("penalty == 0 when all truths interior (level %.2f)", level),
         clean$penalty == 0 && abs(clean$score - clean$width) < 1e-12,
         sprintf("(v1 coverage on the miss case: %.3f)", inside))

  # Mean miss distance is recoverable from the stored penalty column.
  a <- 1 - level
  q <- apply(M, 1, quantile, probs = c(a/2, 1 - a/2), names = FALSE)
  miss <- mean(pmax(q[1, ] - truth, 0) + pmax(truth - q[2, ], 0))
  report(sprintf("penalty/(2/alpha) recovers mean miss distance (%.2f)", level),
         abs(got$penalty / (2/a) - miss) < 1e-12)

  # Properness sanity: at the true predictive distribution the score should
  # beat both an inflated and a deflated interval, on average.
  y <- rnorm(4000)
  score_of <- function(scale) {
    D <- matrix(rnorm(4000 * n_draw, sd = scale), nrow = 4000)
    interval_score_mat(D, 1L, y, level)$score
  }
  s_true <- score_of(1.0); s_wide <- score_of(2.5); s_narrow <- score_of(0.35)
  report(sprintf("calibrated beats over/under-dispersed (level %.2f)", level),
         s_true < s_wide && s_true < s_narrow,
         sprintf("true=%.3f wide=%.3f narrow=%.3f", s_true, s_wide, s_narrow))

  # ATE-level scalar version.
  d <- rnorm(n_draw); tr <- 2.2
  qq <- quantile(d, c(a/2, 1 - a/2), names = FALSE)
  want_vec <- (qq[2] - qq[1]) + (2/a) * max(0, tr - qq[2], qq[1] - tr)
  report(sprintf("interval_score_vec matches closed form (level %.2f)", level),
         abs(interval_score_vec(d, tr, level) - want_vec) < 1e-12)
  report(sprintf("interval_score_vec is width when covered (level %.2f)", level),
         abs(interval_score_vec(d, median(d), level) - (qq[2] - qq[1])) < 1e-12)
}

# Orientation guard: a transposed matrix with mismatched dims must ERROR, not
# recycle. This is the failure mode that silently corrupted the upstream
# mvbart CRPS column.
bad <- tryCatch({ interval_score_mat(matrix(rnorm(500 * 1000), nrow = 500),
                                     1L, rnorm(1000), 0.95); FALSE },
                error = function(e) TRUE)
report("wrong obs_margin errors instead of recycling", bad)

# Known-answer check, computed by hand. Draws 1..100, level 0.95 -> type-7
# quantiles at 2.5% / 97.5% are 3.475 and 97.525; truth 200 is above.
d <- as.numeric(1:100)
q <- quantile(d, c(0.025, 0.975), names = FALSE)
hand <- (q[2] - q[1]) + 40 * (200 - q[2])
# Relative tolerance: the hand value uses the exact multiplier 40, while the
# implementation uses 2/(1 - 0.95) = 39.999999999999964. The gap is 1e-15
# relative and is a property of binary floating point, not of the score.
report("hand-computed known answer",
       abs(interval_score_vec(d, 200, 0.95) - hand) < 1e-10 * hand,
       sprintf("value=%.4f", hand))


# -----------------------------------------------------------------------------
# Schema integration. The helpers can be correct while fill_interval_scores()
# writes to a column name the stub does not define -- R would create the column
# silently, the row would then have the wrong length, and append_csv() would
# reject the shard only after a six-hour fit. Exercise the four model branches
# with their real orientations against the real schema.
# -----------------------------------------------------------------------------
cat("\n-- schema integration --\n")
block <- function(start_pattern, end_pattern) {
  s <- grep(start_pattern, src)[1]
  stopifnot(!is.na(s))
  e <- s - 1 + grep(end_pattern, src[s:length(src)])[1]
  paste(src[s:e], collapse = "\n")
}
# make_stub() through the close of fill_interval_scores(), which also brings in
# CHARACTER_COLUMNS. Slicing by pattern rather than by line number so this test
# breaks loudly if the driver is restructured, instead of testing stale code.
eval(parse(text = block("^E3_SCHEMA_VERSION", "^E3_INTERVAL_LEVELS")))
eval(parse(text = block("^make_stub <- function", "^\\}")))
eval(parse(text = block("^CHARACTER_COLUMNS <- ", "\"interval_levels\"\\)")))
eval(parse(text = block("^fill_interval_scores <- function", "^\\}")))

STUB <- make_stub()
COL_NAMES <- names(STUB)
n_obs <- 200; n_draw <- 300
truth <- rnorm(n_obs, sd = 3)
by_obs <- matrix(rnorm(n_obs * n_draw), nrow = n_obs)   # obs x draws
by_draw <- t(by_obs)                                     # draws x obs
ate_draw <- colMeans(by_obs)

row <- STUB
row[["width_audit_max_abs_dev"]] <- 0
# Populate the v1 width columns the audit compares against, exactly as each
# model branch does.
for (model in c("mvbcf", "bcf", "bart", "mvbart")) {
  for (ss in 1:2) for (lv in c(0.5, 0.95)) {
    tag <- if (lv == 0.5) "50" else "95"
    m <- if (model == "bart") by_draw else by_obs
    row[[paste0(model, "_wid", tag, ss)]] <-
      mean(apply(m, if (model == "bart") 2 else 1, cred_width, lv))
  }
}
# bart is the transposed branch; every other model is (obs x draws).
row <- fill_interval_scores(row, "mvbcf",  1L, by_obs,  1L, truth, ate_draw)
row <- fill_interval_scores(row, "mvbcf",  2L, by_obs,  1L, truth, ate_draw)
row <- fill_interval_scores(row, "bcf",    1L, by_obs,  1L, truth, ate_draw)
row <- fill_interval_scores(row, "bcf",    2L, by_obs,  1L, truth, ate_draw)
row <- fill_interval_scores(row, "bart",   1L, by_draw, 2L, truth, ate_draw)
row <- fill_interval_scores(row, "bart",   2L, by_draw, 2L, truth, ate_draw)
row <- fill_interval_scores(row, "mvbart", 1L, by_obs,  1L, truth, ate_draw)
row <- fill_interval_scores(row, "mvbart", 2L, by_obs,  1L, truth, ate_draw)

report("fill_interval_scores creates no columns outside the schema",
       identical(names(row), COL_NAMES),
       paste(setdiff(names(row), COL_NAMES), collapse = ","))

filled <- grep("_(is|pen)(50|95)[12]$|_ate_is(50|95)[12]$", COL_NAMES, value = TRUE)
report("every interval-score column in the schema was written",
       all(vapply(row[filled], function(v) length(v) == 1L && is.finite(v),
                  logical(1))),
       sprintf("(%d columns)", length(filled)))

report("width audit is at floating-point noise",
       row[["width_audit_max_abs_dev"]] < 1e-9,
       sprintf("max dev = %.3g", row[["width_audit_max_abs_dev"]]))

identity_ok <- TRUE
for (model in c("mvbcf", "bcf", "bart", "mvbart")) {
  for (ss in 1:2) for (tag in c("50", "95")) {
    lhs <- row[[paste0(model, "_is", tag, ss)]]
    rhs <- row[[paste0(model, "_wid", tag, ss)]] + row[[paste0(model, "_pen", tag, ss)]]
    identity_ok <- identity_ok && abs(lhs - rhs) < 1e-10 * max(1, abs(lhs))
  }
}
report("is == wid + pen holds for all 16 model/level/outcome cells", identity_ok)

report("the transposed bart branch matches the untransposed branches",
       abs(row[["bart_is951"]] - row[["mvbcf_is951"]]) < 1e-12)

report("character columns are declared as such",
       all(c("schema_version", "interval_levels") %in% CHARACTER_COLUMNS))

cat(if (fails == 0) "\nALL CHECKS PASSED\n" else sprintf("\n%d CHECK(S) FAILED\n", fails))
quit(status = if (fails == 0) 0 else 1)
