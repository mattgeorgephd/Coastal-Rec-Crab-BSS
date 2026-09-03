###############################################################################
# ZERO-INFLATED NEGATIVE BINOMIAL POSTERIOR-PREDICTIVE ARITHMETIC
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   The catch likelihood became a MIXTURE on 2026-09-02 (`zi_catch = 1`):
#
#       P(Y = 0) = theta + (1 - theta) * NB2(0 | mu, r)
#       P(Y = y) =         (1 - theta) * NB2(y | mu, r)      for y > 0
#
#   Stan's `log_lik` was written against that mixture from the start, so `elpd_loo`,
#   the Pareto k diagnostics and every loo-derived number have always been correct.
#   The R-side posterior-predictive checks were NOT: `pit_block()` in
#   save_run_diagnostics.R and `calib()` in model_diagnostics.R both computed the
#   PIT, the coverage flags and `p_zero` from a plain NB2, i.e. from a likelihood the
#   model no longer used.
#
# WHAT THAT COST (2026-09-03 batch, stage Z1, shore all-gear catch stream)
#   ppc_byobs reported 676 observed zeros against 419.0 expected, z = +15.4, and the
#   batch recorded a REVIEW verdict reading "the zero bin got WORSE" against the
#   NB2 baseline's 605.4 / z = +3.8. The mixture expectation is 637.9, z = +2.0: the
#   zero bin narrows by about half, which is what the parameter was added to do. The
#   shore pot-closure replicate went from z = +2.9 to z = +0.8 the same way. The
#   feature was very nearly rejected on a diagnostic computed under the wrong model.
#   `pit_mean` (0.4325 against 0.4984), `pit_sd`, `coverage_50` (0.405 against 0.457)
#   and the `flag_pit_bias = TRUE` that Z1 raised on the shore all-gear fit are the
#   same artefact: with theta absorbing the excess zeros, lambda_C rises, and scoring
#   the observed zeros against the risen NB2 mean alone drags the PIT down.
#
# THE RULE THIS FILE ENFORCES
#   Every R-side predictive statistic is computed under the SAME likelihood Stan
#   sampled. `theta = NULL` reproduces the pre-2026-09-04 NB2 arithmetic bit for bit,
#   so every stream other than catch, and every run with `estimate_catch_zi = FALSE`,
#   is untouched. Only the catch stream of a zi_catch = 1 fit takes the mixture path.
#
# THE PIT IS NOT THE OBVIOUS FORMULA
#   The randomized PIT is F(y-1) + U * f(y), U ~ Uniform(0, 1), and its expectation
#   over U is F(y-1) + 0.5 * f(y). Under the mixture F(-1) = 0 while F(y) for y >= 0
#   carries the whole point mass theta, so the y = 0 case does NOT reduce to
#   theta + (1 - theta) * pit_nb2:
#
#       y = 0 :  0.5 * (theta + (1 - theta) * NB2_pmf(0))
#       y > 0 :  theta + (1 - theta) * (NB2_cdf(y-1) + 0.5 * NB2_pmf(y))
#
#   Folding the two lines into one expression is the natural mistake and it puts the
#   zero observations at roughly theta too high, which flatters the calibration of
#   exactly the observations the parameter was introduced to explain.
###############################################################################

# Randomized-PIT expectation for one observation, averaged over posterior draws.
#   y     scalar observed count
#   mu    numeric vector of draws of the NB2 mean
#   size  numeric vector of draws of the NB2 dispersion (same length as mu)
#   theta NULL (plain NB2) or a numeric vector of draws of the mixing weight
bss_zi_pit <- function(y, mu, size, theta = NULL) {
  if (is.null(theta)) {
    return(mean(stats::pnbinom(y - 1, size = size, mu = mu) +
                0.5 * stats::dnbinom(y, size = size, mu = mu)))
  }
  if (y == 0) {
    mean(0.5 * (theta + (1 - theta) * stats::dnbinom(0, size = size, mu = mu)))
  } else {
    mean(theta + (1 - theta) * (stats::pnbinom(y - 1, size = size, mu = mu) +
                                0.5 * stats::dnbinom(y, size = size, mu = mu)))
  }
}

# Model-implied P(Y = 0) for one observation, averaged over posterior draws. This is
# summed over ALL observations to form the expected zero count; see the note in
# save_run_diagnostics.R about why it must not be read off the observed zeros only.
bss_zi_p_zero <- function(mu, size, theta = NULL) {
  nb0 <- stats::dnbinom(0, size = size, mu = mu)
  if (is.null(theta)) mean(nb0) else mean(theta + (1 - theta) * nb0)
}

# Pull theta_C draws for the CATCH stream, or NULL when the fit is plain NB2.
#   Returns a numeric vector on the `use` draw subset, or NULL. NULL is returned for
#   any failure, so a fit that predates the parameter, or an extraction that throws,
#   degrades to the NB2 arithmetic rather than aborting the diagnostics block.
bss_zi_theta_draws <- function(fit, stan_data, use) {
  if (!identical(as.integer(stan_data$zi_catch %||% 0L), 1L)) return(NULL)
  th <- try(rstan::extract(fit, pars = "theta_C")$theta_C, silent = TRUE)
  if (inherits(th, "try-error") || is.null(th)) return(NULL)
  th <- as.numeric(th)
  if (length(th) < max(use)) return(NULL)
  th[use]
}
