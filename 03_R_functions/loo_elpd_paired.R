###############################################################################
# PAIRED elpd COMPARISON BETWEEN TWO FITS ON THE SAME OBSERVATIONS
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   On 2026-09-03 the ZINB prototype (stage Z1) was scored against the NB2 baseline
#   (stage Z0) by comparing the elpd DIFFERENCE, +14.8 nats, against `se_elpd_loo`
#   from each fit's own loo_summary, about 46 nats. The batch recorded NOT WORTH IT.
#   That comparison is wrong, and the error is not small: it rejected a real effect.
#
#   `se_elpd_loo` is the standard error of ONE model's total elpd. It is dominated by
#   variation ACROSS observations, i.e. by the fact that some interviews are easier to
#   predict than others. That variation is COMMON to both models, and it cancels in
#   the difference. The quantity a model comparison needs is the standard error of the
#   PAIRED difference:
#
#       d_i = elpd_loo_i(B) - elpd_loo_i(A)      (same observation i in both fits)
#       SE  = sd(d) * sqrt(n)
#
#   This is the standard `loo_compare()` statistic (Vehtari, Gelman & Gabry 2017,
#   "Practical Bayesian model evaluation using leave-one-out cross-validation and
#   WAIC", Statistics and Computing 27:1413-1432, section 3.3; and the loo package's
#   own `loo_compare` documentation, which reports elpd_diff alongside se_diff for
#   exactly this reason). On the Z1/Z0 shore all-gear catch stream the paired SE is
#   5.5 nats against the naive 46.7, so +14.8 nats is 2.69 SE, not 0.32 SE.
#
# WHAT IT STILL WILL NOT TELL YOU
#   Vehtari et al. caution that the normal approximation behind se_diff is unreliable
#   when n is small (they suggest roughly n < 100) or when the models are very close,
#   and that the SE is itself estimated. Treat 2 SE as evidence, not proof. The far
#   more informative output is the DECOMPOSITION: this function splits the difference
#   by whether the observation is zero, because a zero-inflation parameter that earns
#   its elpd purely at the zeros while degrading the positive counts is doing
#   something different from one that improves the fit everywhere, and the aggregate
#   number cannot distinguish them. On Z1/Z0 the split was +25.2 at the 676 zeros and
#   -10.5 at the 973 positives: net positive, but bought.
#
# INPUTS are the loo_pointwise_<stream>_<label>.csv files both runs already write.
###############################################################################

loo_elpd_paired <- function(path_a, path_b, label = NULL) {
  if (!file.exists(path_a) || !file.exists(path_b)) return(NULL)
  a <- utils::read.csv(path_a, stringsAsFactors = FALSE)
  b <- utils::read.csv(path_b, stringsAsFactors = FALSE)
  if (!all(c("obs_index", "elpd_loo", "observed") %in% names(a))) return(NULL)
  # Align on obs_index rather than assuming row order. A silently misaligned join would
  # produce a difference made of noise with a plausible-looking SE, which is the one
  # failure mode of this statistic that does not announce itself.
  m <- merge(a[, c("obs_index", "observed", "elpd_loo", "pareto_k")],
             b[, c("obs_index", "observed", "elpd_loo", "pareto_k")],
             by = "obs_index", suffixes = c("_a", "_b"))
  if (!nrow(m) || !identical(m$observed_a, m$observed_b)) return(NULL)
  d <- m$elpd_loo_b - m$elpd_loo_a
  n <- length(d)
  se <- stats::sd(d) * sqrt(n)
  iz <- m$observed_a == 0
  sub <- function(k) {
    if (sum(k) < 2) return(c(diff = NA_real_, se = NA_real_, n = sum(k)))
    c(diff = sum(d[k]), se = stats::sd(d[k]) * sqrt(sum(k)), n = sum(k))
  }
  list(label = label, n = n,
       elpd_a = sum(m$elpd_loo_a), elpd_b = sum(m$elpd_loo_b),
       elpd_diff = sum(d), se_diff = se,
       ratio = if (is.finite(se) && se > 0) sum(d) / se else NA_real_,
       # The naive statistic, reported ONLY so a reader comparing this against an older
       # review can see why the two disagree.
       se_naive_a = stats::sd(m$elpd_loo_a) * sqrt(n),
       se_naive_b = stats::sd(m$elpd_loo_b) * sqrt(n),
       zeros = sub(iz), positives = sub(!iz),
       # 2026-09-05: the decomposition that actually explains the sign of `positives`. On the
       # ZINB comparison "positives -10.5" was y=1 at -42.0 and every bin from 3 up POSITIVE,
       # with the 3+ bins carrying 87% of the catch. Read this table, not the two-way split.
       by_count = local({
         bins <- cut(m$observed_a, c(-1, 0, 1, 2, 4, 8, 16, Inf),
                     labels = c("0", "1", "2", "3-4", "5-8", "9-16", "17+"))
         t(vapply(levels(bins), function(l) {
           k <- bins == l
           c(n = sum(k), diff = if (any(k)) sum(d[k]) else NA_real_,
             se = if (sum(k) >= 2) stats::sd(d[k]) * sqrt(sum(k)) else NA_real_,
             catch_share = sum(m$observed_a[k]) / max(1, sum(m$observed_a)))
         }, numeric(4)))
       }),
       n_pareto_bad_a = sum(m$pareto_k_a > 0.7, na.rm = TRUE),
       n_pareto_bad_b = sum(m$pareto_k_b > 0.7, na.rm = TRUE))
}

# One-line rendering for a batch verdict's `observed` field.
loo_elpd_paired_str <- function(x) {
  if (is.null(x)) return("loo pointwise files not comparable")
  sprintf(paste0("elpd %.1f -> %.1f, diff %+.1f, paired SE %.1f (%.2f SE); ",
                 "naive per-fit SE would have read %.1f. Split: zeros n=%d %+.1f (%.1f SE), ",
                 "positives n=%d %+.1f (%.1f SE). Pareto k>0.7: %d -> %d"),
          x$elpd_a, x$elpd_b, x$elpd_diff, x$se_diff, x$ratio, x$se_naive_b,
          x$zeros[["n"]], x$zeros[["diff"]], x$zeros[["diff"]] / x$zeros[["se"]],
          x$positives[["n"]], x$positives[["diff"]], x$positives[["diff"]] / x$positives[["se"]],
          x$n_pareto_bad_a, x$n_pareto_bad_b)
}

# The by-count table as one line, for a verdict's `observed` field.
loo_elpd_by_count_str <- function(x) {
  if (is.null(x) || is.null(x$by_count)) return("")
  b <- x$by_count
  paste(sprintf("y=%s n=%d %+.1f (%.1f SE, %.0f%% of catch)", rownames(b), as.integer(b[, "n"]),
                b[, "diff"], b[, "diff"] / b[, "se"], 100 * b[, "catch_share"]), collapse = " | ")
}
