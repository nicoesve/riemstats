# Statistical correctness of riem_anova's p-value construction.
#
# These tests were added on 2026-08-27 alongside the fix for two defects in
# riem_anova() found by operation package-test-hardening (T6, 2026-08-20):
#
#   * the p-value omitted the observed labelling from the reference set,
#     returning #{v < t_obs} / nperm rather than
#     (1 + #{at least as extreme}) / (nperm + 1); and
#   * the tail was hard-coded to the left, which is correct for
#     log_wilks_lambda and backwards for pillais_trace, so Pillai returned
#     1 - p silently.
#
# The tail is now fixed unconditionally. The p-value construction is exposed as
# an argument, p_method, defaulting to "plugin" -- the historical construction
# -- so that results published from earlier versions reproduce by the default
# call (Control's ruling, 2026-08-27). p_method = "exact" is the
# level-preserving estimator. See ?riem_anova.
#
# Skipped on CRAN (time budget); run locally and in CI.

library(riemtan)
data("airm")

gen_conn <- function(p, sd = 0.3, shift = 1) {
  Q <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  d <- shift * exp(stats::rnorm(p, 0, sd))
  M <- Q %*% diag(d, p, p) %*% t(Q)
  Matrix::pack(Matrix::nearPD((M + t(M)) / 2)$mat)
}
gen_group <- function(p, n, shift = 1) {
  lapply(seq_len(n), function(i) gen_conn(p, shift = shift)) |>
    CSample$new(metric_obj = airm)
}

test_that("riem_anova's default reproduces its historical p-value construction", {
  skip_on_cran()
  # p_method = "plugin" is the default precisely so that results published from
  # earlier versions reproduce by the default call. This pins that: the default
  # must equal mean(stat_obs > perm_vals) computed independently here, which is
  # exactly what riem_anova() returned before p_method existed. If the default
  # is ever flipped, this is the test that must be changed deliberately rather
  # than silently.
  set.seed(11)
  ss <- list(gen_group(2, 6, shift = 1), gen_group(2, 6, shift = 2.5)) |>
    CSuperSample$new()

  # Warm the cached Within/Total so both routes below consume the RNG alike.
  invisible(suppressMessages(log_wilks_lambda(ss)))

  nperm <- 19L

  set.seed(11)
  p_default <- suppressMessages(riem_anova(ss, log_wilks_lambda, nperm = nperm))

  set.seed(11)
  stat_obs <- suppressMessages(log_wilks_lambda(ss))
  perm_vals <- vapply(
    seq_len(nperm),
    \(i) suppressMessages(one_permutation(ss, log_wilks_lambda)),
    numeric(1)
  )

  expect_equal(p_default, mean(stat_obs > perm_vals))
})

test_that("riem_anova's exact p_method includes the observed labelling", {
  skip_on_cran()
  # A correctly constructed permutation p-value can never be 0 -- the observed
  # labelling is always at least as extreme as itself -- so its minimum is
  # 1/(nperm + 1). That floor is the observable signature of the
  # (1 + k)/(B + 1) construction and is what p_method = "exact" buys.
  set.seed(4)
  ss <- list(
    gen_group(2, 6, shift = 1),
    gen_group(2, 6, shift = 5) # separated hard enough to sit in the tail
  ) |> CSuperSample$new()
  invisible(suppressMessages(log_wilks_lambda(ss)))

  nperm <- 19L

  # Same seed for both, so the two estimators see the same permutations and the
  # comparison below is about the construction, not about Monte-Carlo noise.
  set.seed(7)
  p_exact <- suppressMessages(
    riem_anova(ss, log_wilks_lambda, nperm = nperm, p_method = "exact")
  )
  set.seed(7)
  p_plugin <- suppressMessages(
    riem_anova(ss, log_wilks_lambda, nperm = nperm, p_method = "plugin")
  )

  expect_gte(p_exact, 1 / (nperm + 1))
  # (1 + #{<=})/(B + 1) >= #{<}/B for every attainable count, so the valid
  # estimator is never below the unbiased one on the same permutations.
  expect_gte(p_exact, p_plugin)
})

test_that("riem_anova orients its p-value to the statistic it is given", {
  skip_on_cran()
  # riem_anova's p-value used to be an unconditional LEFT-tail count. That is
  # correct for log_wilks_lambda (Lambda = |W|/|T| is SMALL when the groups
  # differ) and backwards for pillais_trace (tr((T-W)T^-1) is LARGE when they
  # differ), which is exported and documented as a stat_fun for this function.
  # Used that way it returned 1 - p, and it failed silently: a real difference
  # produced a p-value near 1 rather than an error.
  set.seed(2)
  ss_alt <- list(gen_group(2, 6, shift = 1), gen_group(2, 6, shift = 3)) |>
    CSuperSample$new()

  # Both statistics detect the same difference, so both p-values must be small.
  expect_lt(suppressMessages(riem_anova(ss_alt, log_wilks_lambda, nperm = 19)), 0.2)
  expect_lt(suppressMessages(riem_anova(ss_alt, pillais_trace, nperm = 19)), 0.2)
})

test_that("the shipped statistics declare their tail, and tail= overrides it", {
  skip_on_cran()
  # The orientation lives on the statistic, not in a lookup inside riem_anova,
  # so a user-supplied statistic can declare its own.
  expect_identical(attr(log_wilks_lambda, "riemstats_tail"), "left")
  expect_identical(attr(pillais_trace, "riemstats_tail"), "right")

  set.seed(2)
  ss_alt <- list(gen_group(2, 6, shift = 1), gen_group(2, 6, shift = 3)) |>
    CSuperSample$new()
  invisible(suppressMessages(pillais_trace(ss_alt)))

  nperm <- 19L
  set.seed(5)
  p_auto <- suppressMessages(riem_anova(ss_alt, pillais_trace, nperm = nperm))
  set.seed(5)
  p_forced <- suppressMessages(
    riem_anova(ss_alt, pillais_trace, nperm = nperm, tail = "left")
  )

  expect_lt(p_auto, 0.2)
  # Forcing the wrong tail recovers the old inverted behaviour. The two plugin
  # counts are exact complements only when no permuted statistic ties the
  # observed one -- ties are measure-zero for these statistics on Wishart draws
  # but are common for discrete ones -- so assert the inequality that always
  # holds rather than the equality that merely usually does.
  expect_gt(p_forced, 0.8)
  expect_lte(p_auto + p_forced, 1)
})
