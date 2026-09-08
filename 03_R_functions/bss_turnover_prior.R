# -----------------------------------------------------------------------------
# Part of Coastal-Rec-Crab-BSS: recreational Dungeness crab creel estimation
# for Grays Harbor / Westport (WDFW).
# Copyright (C) 2024-2026 Washington Department of Fish and Wildlife.
#
# Adapted from CreelEstimates, the WDFW freshwater creel estimation framework:
#   https://github.com/dfw-wa/CreelEstimates   (licensed GPL-3.0).
# Substantial portions of the methodology, structure, and R/Stan code originate
# in CreelEstimates and remain (C) their authors under GPL-3.0; changes for
# recreational crab are by WDFW.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License, version 3, as published by the Free
# Software Foundation. It is distributed WITHOUT ANY WARRANTY; without even the
# implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program (see the LICENSE file);
# if not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
###############################################################################
# bss_turnover_prior.R  (2026-09-08, review item 3 + item 5)
#
# Resolve the BOAT deployment-turnover prior centre (tau_boat_prior_mu) from the
# season's own OSP/trailer overlap calibration, and hand the same number to the
# Point Estimator, so the prior, the PE and the BSS all expand the trailer count
# by one quantity.
#
# WHY THIS EXISTS
#   Until 2026-09-08 tau_boat_prior_mu was 1.2, set from two WBL ingress/egress
#   days. The shared-turnover model (tau_bar, adopted 2026-09-01) is identified by
#   the 61 OSP/trailer overlap days and lands at 2.60 [2.06, 3.25]; inverting the
#   posterior on the log scale, the likelihood alone puts it at about 2.97, so the
#   1.2 prior was shrinking the boat component by roughly 13%, toward a centre the
#   same season's data contradict. Meanwhile the PE still expanded on 1.2, so the
#   boat PE-vs-BSS comparison (3x apart) was mostly the two turnovers disagreeing,
#   not the two estimators.
#
# THE RULE
#   tau_boat_prior_mu = "calibration"  -> the implied turnover from
#       diagnose_osp_trailer_overlap() for THIS window, metric
#       params$tau_boat_calibration_metric (default "trailer_mean_per_visit").
#       The mean-per-visit metric is the consistent one: the Stan likelihood
#       treats EVERY trailer count as an observation of lambda_E / R_G_boat, so
#       the model's trailer level is the per-visit mean, not the daily maximum,
#       and the PE expands mean_count per day the same way.
#   tau_boat_prior_mu = <number>       -> used as is (the historical behaviour).
#   No overlap (OSP absent, or < 3 paired days) -> params$tau_boat_prior_mu_fallback
#       (SEASON-DERIVED in run_config.R; 2.7 = the 2024-25 max-per-day calibration).
#
# ON THE DOUBLE USE OF THE OVERLAP DAYS. The same 61 days sit in the likelihood, so
#   centring the prior on them is a mild empirical-Bayes step. It is harmless only
#   because the prior is WIDE: with tau_boat_prior_sigma = 0.5 the prior precision is
#   4 against a likelihood precision near 64 (about 6% weight), and the posterior SD
#   shrinks by roughly 3%. Do not tighten this prior below 0.3 without removing the
#   overlap days from one side or the other.
#
# The resolved value REPLACES params$tau_boat_prior_mu (numeric) so every existing
# consumer (bss_effort_spec, run_pe_pooled / run_pe_gear, pe_monthly_effort_share,
# diagnose_tau_boat_sensitivity, the gear driver's daily-combined series) reads one
# number. The provenance is kept in params$tau_boat_prior_source and
# params$tau_boat_prior_calibration_table.
###############################################################################

bss_resolve_tau_boat_prior <- function(params, osp_overlap = NULL, quiet = FALSE) {
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  raw <- params$tau_boat_prior_mu %||% 1.2
  metric   <- params$tau_boat_calibration_metric %||% "trailer_mean_per_visit"
  fallback <- suppressWarnings(as.numeric(params$tau_boat_prior_mu_fallback %||% 2.7))
  min_pairs <- as.integer(params$tau_boat_calibration_min_pairs %||% 3L)

  # Numeric: historical behaviour, nothing to resolve.
  if (is.numeric(raw) && length(raw) == 1L && is.finite(raw) && raw > 0) {
    params$tau_boat_prior_source <- "config (numeric)"
    return(params)
  }
  if (!(is.character(raw) && identical(tolower(raw), "calibration")))
    stop("params$tau_boat_prior_mu must be a positive number or the string \"calibration\" (got ",
         deparse(raw), ").", call. = FALSE)

  cal <- osp_overlap$calibration
  if (is.null(cal) || !is.data.frame(cal) || !nrow(cal) ||
      !all(c("trailer_metric", "n", "implied_turnover") %in% names(cal))) {
    .say(sprintf(paste0("  tau_boat prior: no OSP/trailer overlap calibration available; using the ",
                        "fallback centre %.2f (tau_boat_prior_mu_fallback).\n"), fallback))
    params$tau_boat_prior_mu     <- fallback
    params$tau_boat_prior_source <- "fallback (no overlap calibration)"
    return(params)
  }
  row <- cal[cal$trailer_metric == metric, , drop = FALSE]
  if (!nrow(row))
    stop("tau_boat_calibration_metric '", metric, "' is not in the overlap calibration; choose one of: ",
         paste(cal$trailer_metric, collapse = ", "), call. = FALSE)
  val <- suppressWarnings(as.numeric(row$implied_turnover[1]))
  n   <- suppressWarnings(as.integer(row$n[1]))
  if (!is.finite(val) || val <= 0 || !is.finite(n) || n < min_pairs) {
    .say(sprintf(paste0("  tau_boat prior: overlap calibration unusable (metric %s, n = %s, value %s); ",
                        "using the fallback centre %.2f.\n"), metric, n, val, fallback))
    params$tau_boat_prior_mu     <- fallback
    params$tau_boat_prior_source <- sprintf("fallback (calibration n = %s below %d)", n, min_pairs)
    return(params)
  }
  .say(sprintf(paste0("  tau_boat prior centre RESOLVED from the OSP/trailer overlap: %.3f ",
                      "(%s, n = %d paired days; prior log-SD %.2f). The PE expands on the same value.\n"),
               val, metric, n, as.numeric(params$tau_boat_prior_sigma %||% 0.5)))
  params$tau_boat_prior_mu                <- val
  params$tau_boat_prior_source            <- sprintf("calibration (%s, n = %d)", metric, n)
  params$tau_boat_prior_calibration_table <- cal
  params
}
