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
# crab_fraction.R  (Phase 2 + Phase 3 + improvement 8: directed-crabbing fraction f)
#
# f = the share of private boats at the port that are crabbing (vs. targeting other
# fisheries). The effort series (trailer counts, OSP boat totals) count ALL private
# boats, so f converts all-boat effort to crab effort. Without f the model implicitly
# sets f = 1 (every boat crabbing), which biases the boat catch high.
#
# ============================================================================
# IMPROVEMENT 8 (2026-08-25): TWO INFORMATION SOURCES, ONE f
# ============================================================================
#
# There are now two possible observations of the crab share, and they measure DIFFERENT
# things. Getting that difference right is the whole design.
#
#   (A) WPT/WBL EGRESS CLASSIFICATION  (ingress_egress.xlsx: boats_crabbing / boats_total)
#       Our own samplers classify boats as crabbing or not. A boat that crabs AND fishes
#       something else counts as CRABBING. This is a direct, unbiased observation of f:
#             n_crab[k] ~ Binomial(n_total[k], f[k])
#
#   (B) OSP CRAB-ONLY COUNTS  (WBL_boat_counts.xlsx: the crab-only column)
#       OSP records a fishery label per boat and does NOT record combo trips: a boat
#       crabbing AND fishing another fishery is labelled by the OTHER fishery. So the OSP
#       crab-only share is a LOWER BOUND on f, not an estimate of it. Treating it as an
#       estimate would bias the boat harvest DOWN by exactly the combo-trip rate, which
#       on a halibut or tuna day is the failure mode that matters most.
#
# The parameterization makes the bound structural rather than advisory:
#
#       f_lower[k] ~ Beta(1,1),   osp_crab_only[k] ~ Binomial(osp_total[k], f_lower[k])
#       f[k]       = f_lower[k] + (1 - f_lower[k]) * theta[k]
#
# theta[k] in [0,1] is the share of the NOT-crab-labelled boats that were also crabbing,
# i.e. the combo trips OSP cannot see. f can never fall below what OSP directly observed.
#
# IDENTIFIABILITY, STATED PLAINLY. OSP alone identifies f_lower, never theta: the data
# cannot distinguish "few crab boats" from "many combo trips". theta rests on its prior
# (crab_fraction_combo_share, a placeholder with the same standing f = 0.3 had before the
# egress pilot) until the egress classification covers the same stratum, at which point
# the Binomial in (A) pins f and therefore theta. That is the strongest operational
# argument for scheduling egress classification days on days OSP is also in port.
#
# HOW IT DEGRADES (the "defaults back to interview data" requirement):
#   OSP + egress   -> f_lower from OSP, theta pinned by egress. Both identified.
#   OSP only       -> f bounded below by data; the level rests on the combo prior.
#   egress only    -> f_lower = 0, f = theta, prior = the ordinary f prior. EXACTLY the
#                     Phase 2/3 behaviour.
#   neither        -> f = theta ~ Beta(set*kappa, (1-set)*kappa). EXACTLY today.
# The degradation is exact, not approximate: when a stratum has no OSP classification the
# R side hands Stan the ordinary f prior for theta instead of the combo prior, so the
# posterior for f is bit-identical to the pre-improvement-8 model.
#
# WHAT DOES NOT CHANGE. f still enters the Stan GENERATED QUANTITIES only. Both new
# likelihood terms are Binomials on OBSERVED boat counts, never on the latent effort, so
# the boat total stays exactly linear in f and the model CPUE stays invariant. That was
# the validated Phase 2/3 property and it is deliberately preserved.
#
# ============================================================================
# PHASE 3: f IS PER STRATUM
# ============================================================================
# crab_fraction_strata selects the stratification: "none" (one scalar = Phase 2),
# "month", "day_type", "month_day_type", and (improvement 4/8) "opener" /
# "month_opener", which key off whether a competing fishery was open that day. The
# opener strata exist because an opener EFFORT covariate and a constant f fight each
# other: fitting the halibut-day boat surge more faithfully and then multiplying it by a
# flat f makes the boat catch more biased on those days, not less.
#
# CONFIG (all optional; defaults shown)
#   use_crab_fraction         FALSE
#   crab_fraction_strata      "none"
#   crab_fraction_set         0.3     (set value = theta's prior mean when no OSP, and the fallback)
#   crab_fraction_prior_kappa 20      (Beta concentration for that prior)
#   crab_fraction_fixed       NA      (a number PINS f exactly, no uncertainty)
#   crab_fraction_min_obs     20      (min egress-classified boats before its Binomial binds)
#   use_osp_crab_lower        FALSE   (improvement 8; inert without a crab-only column)
#   crab_fraction_osp_min_obs 20      (min OSP-classified boats before the bound binds)
#   crab_fraction_combo_share 0.15    (theta's prior mean in an OSP-informed stratum)
#   crab_fraction_combo_kappa 8
###############################################################################

# Per-date stratum labels, derived from the DATE + config (so the day set, the egress
# classification rows and the OSP rows are labelled by one consistent rule).
# "opener" / "month_opener" additionally need params$opener_f_dates, a Date vector of
# days on which the chosen competing fishery was open (the driver builds it from
# params$opener_f_flag); absent -> every day is "shut", which reduces "opener" to a
# single stratum and "month_opener" to "month".
crab_fraction_strata_labels <- function(dates, params) {
  mode  <- params$crab_fraction_strata %||% "none"
  dates <- as.Date(dates)
  if (identical(mode, "none")) return(rep("all", length(dates)))
  mo <- format(dates, "%m")
  wknd_days <- params$days_wkend %||% c("Saturday", "Sunday")
  dt <- ifelse(weekdays(dates) %in% wknd_days, "wknd", "wkdy")
  hol <- params$crabbing_holiday_dates
  if (!is.null(hol)) dt[dates %in% as.Date(hol)] <- "wknd"   # holidays typed as weekend
  op_dates <- params$opener_f_dates %||% as.Date(character(0))
  op <- ifelse(dates %in% as.Date(op_dates), "open", "shut")
  switch(mode,
    "month"          = mo,
    "day_type"       = dt,
    "month_day_type" = paste(mo, dt, sep = "_"),
    "opener"         = op,
    "month_opener"   = paste(mo, op, sep = "_"),
    stop("params$crab_fraction_strata must be none|month|day_type|month_day_type|opener|month_opener (got '",
         mode, "')", call. = FALSE))
}


# Aggregate a per-row classification table (event_date + numerator + denominator) into
# per-stratum integer counts. `restrict_to` optionally limits the rows to a fit's own day
# set, which is what the OSP stream wants (each sub-season's f should reflect its own
# window). Returns list(n_total, n_crab), both length K.
.cf_aggregate <- function(rows, strata, params, num_col, den_col, restrict_to = NULL) {
  K <- length(strata)
  n_total <- rep(0L, K); n_crab <- rep(0L, K)
  if (is.null(rows) || !is.data.frame(rows) || nrow(rows) == 0) return(list(n_total = n_total, n_crab = n_crab))
  if (!all(c(num_col, den_col, "event_date") %in% names(rows))) return(list(n_total = n_total, n_crab = n_crab))
  if (!is.null(restrict_to)) rows <- rows[as.Date(rows$event_date) %in% as.Date(restrict_to), , drop = FALSE]
  if (nrow(rows) == 0) return(list(n_total = n_total, n_crab = n_crab))
  rlab <- crab_fraction_strata_labels(rows$event_date, params)
  for (k in seq_len(K)) {
    sel <- rlab == strata[k]
    if (any(sel)) {
      n_total[k] <- as.integer(round(sum(rows[[den_col]][sel], na.rm = TRUE)))
      n_crab[k]  <- as.integer(round(sum(rows[[num_col]][sel], na.rm = TRUE)))
      if (n_crab[k] > n_total[k]) n_crab[k] <- n_total[k]
    }
  }
  list(n_total = n_total, n_crab = n_crab)
}


# Build the Stan data for the (possibly per-stratum) crabbing fraction f.
#   is_shore : shore fits get apply_crab_fraction = 0 (f pinned to 1)
#   days     : the fit's day set (needs event_date); drives n_f_strata + f_stratum
#   params   : run_config, plus params$crab_fraction_rows (egress classification) and
#              params$osp_crab_rows (OSP total + crab-only per day), both lifted by the
#              driver from the readers' attributes.
# n_f_strata = 1 with no OSP rows reproduces the Phase 2 scalar exactly.
crab_fraction_stan_data <- function(is_shore, days, params, quiet = FALSE) {
  D <- nrow(days)
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  apply_cf <- as.integer(!is_shore && isTRUE(params$use_crab_fraction))

  labels <- crab_fraction_strata_labels(days$event_date, params)
  strata <- sort(unique(labels))
  K <- length(strata)
  f_stratum <- match(labels, strata)   # 1..K

  neutral <- list(
    apply_crab_fraction = 0L, crab_fraction_estimate = 0L,
    n_f_strata = 1L, f_stratum = as.array(rep(1L, D)),
    crab_fraction_value = as.array(1.0), crab_fraction_alpha0 = as.array(1.0),
    crab_fraction_beta0 = as.array(1.0), crab_fraction_n_total = as.array(0L),
    crab_fraction_n_crab = as.array(0L),
    osp_crab_lower = 0L, osp_f_n_total = as.array(0L), osp_f_n_crab = as.array(0L),
    OSPF_n = 0L, osp_f_stratum = integer(0), osp_f_total = integer(0),
    osp_f_crab = integer(0),
    osp_f_kappa_prior_mu = params$crab_fraction_osp_kappa_prior_mu %||% 20)
  if (apply_cf == 0L) return(neutral)

  .clamp <- function(x) pmin(pmax(as.numeric(x), 1e-4), 1 - 1e-4)
  set_val <- .clamp(params$crab_fraction_set %||% 0.3)
  kappa   <- params$crab_fraction_prior_kappa %||% 20
  min_obs <- params$crab_fraction_min_obs     %||% 20L
  fixed   <- params$crab_fraction_fixed

  base <- list(
    apply_crab_fraction = 1L, n_f_strata = as.integer(K),
    f_stratum = as.array(as.integer(f_stratum)))

  # Hard set value: pin f per stratum, no uncertainty (sensitivity lever). The OSP bound
  # is meaningless against a pinned f, so it is switched off here.
  if (!is.null(fixed) && length(fixed) == 1L && is.finite(suppressWarnings(as.numeric(fixed)))) {
    fv <- .clamp(fixed)
    .say(sprintf("  Crab fraction f PINNED at %.3f across %d stratum/strata (%s); boat catch x %.3f.\n",
                 fv, K, params$crab_fraction_strata %||% "none", fv))
    return(c(base, list(
      crab_fraction_estimate = 0L,
      crab_fraction_value  = as.array(rep(fv, K)),
      crab_fraction_alpha0 = as.array(rep(1.0, K)), crab_fraction_beta0 = as.array(rep(1.0, K)),
      crab_fraction_n_total = as.array(rep(0L, K)), crab_fraction_n_crab = as.array(rep(0L, K)),
      osp_crab_lower = 0L,
      osp_f_n_total = as.array(rep(0L, K)), osp_f_n_crab = as.array(rep(0L, K)),
      OSPF_n = 0L, osp_f_stratum = integer(0), osp_f_total = integer(0),
      osp_f_crab = integer(0),
      osp_f_kappa_prior_mu = params$crab_fraction_osp_kappa_prior_mu %||% 20)))
  }

  # --- (A) egress classification: crab-vs-total boats, combo trips INCLUDED -----------
  # Both f streams are aggregated over the SAME date window (2026-08-25). Previously
  # the egress rows were pooled season-wide (indeed across seasons) while the OSP rows
  # were restricted to the fit's day set, so under crab_fraction_strata = "none" the two
  # sub-season fits saw the same egress counts but different OSP counts -- the stratum
  # theta was pinned on and the stratum f_lower was bounded by described different date
  # sets. crab_fraction_restrict_to_fit = FALSE restores the old pooling. This is
  # behaviour-neutral today (the egress classification columns are still blank).
  .eg_window <- if (isTRUE(params$crab_fraction_restrict_to_fit %||% TRUE)) days$event_date else NULL
  eg  <- .cf_aggregate(params$crab_fraction_rows, strata, params,
                       num_col = "boats_crabbing", den_col = "boats_total",
                       restrict_to = .eg_window)
  use_eg <- eg$n_total >= min_obs
  eg_nt  <- ifelse(use_eg, eg$n_total, 0L)
  eg_nc  <- ifelse(use_eg, pmin(eg$n_crab, eg$n_total), 0L)

  # --- (B) OSP crab-only counts: the LOWER bound (improvement 8) ----------------------
  osp_on      <- isTRUE(params$use_osp_crab_lower) && isTRUE(params$use_osp_boat_counts)
  osp_min_obs <- params$crab_fraction_osp_min_obs %||% 20L
  osp <- if (osp_on) .cf_aggregate(params$osp_crab_rows, strata, params,
                                   num_col = "osp_crab_only", den_col = "osp_total",
                                   restrict_to = days$event_date)
         else list(n_total = rep(0L, K), n_crab = rep(0L, K))
  # pmax(osp_min_obs, 1) keeps this condition identical to Stan's `osp_f_n_total[k] > 0`
  # test in transformed parameters. Without it, setting crab_fraction_osp_min_obs = 0
  # would hand a stratum with ZERO OSP boats the combo prior (mean 0.15) while Stan
  # pinned its f_lower to 0, collapsing f onto the combo placeholder instead of the
  # ordinary f prior.
  use_osp <- osp$n_total >= pmax(as.numeric(osp_min_obs), 1)
  osp_nt  <- ifelse(use_osp, osp$n_total, 0L)
  osp_nc  <- ifelse(use_osp, pmin(osp$n_crab, osp$n_total), 0L)
  osp_flag <- as.integer(osp_on && any(use_osp))

  # PER-DAY rows for the beta-binomial likelihood. The per-stratum sums above are the
  # has-data flag and the reporting figure ONLY: a Binomial on them would treat every
  # boat-day as an independent trial and pin f_lower far harder than ~150 correlated
  # daily observations support (see the Stan data block). Only days in strata that
  # cleared the minimum are sent.
  ospf <- list(stratum = integer(0), total = integer(0), crab = integer(0))
  if (osp_flag == 1L) {
    rows <- params$osp_crab_rows
    rows <- rows[as.Date(rows$event_date) %in% as.Date(days$event_date), , drop = FALSE]
    if (nrow(rows) > 0) {
      rk <- match(crab_fraction_strata_labels(rows$event_date, params), strata)
      keep <- !is.na(rk) & use_osp[rk] &
              is.finite(rows$osp_total) & rows$osp_total > 0 &
              is.finite(rows$osp_crab_only) & rows$osp_crab_only >= 0
      if (any(keep)) {
        ospf$stratum <- as.integer(rk[keep])
        ospf$total   <- as.integer(round(rows$osp_total[keep]))
        ospf$crab    <- pmin(as.integer(round(rows$osp_crab_only[keep])), ospf$total)
      }
    }
  }

  # --- theta's prior, per stratum -----------------------------------------------------
  # OSP-informed stratum: theta is the COMBO share (of the not-crab-labelled boats), so
  # it takes the combo prior. Otherwise f = theta and theta takes the ordinary f prior,
  # which is what makes the no-OSP path bit-identical to Phase 2/3.
  combo_mu    <- .clamp(params$crab_fraction_combo_share %||% 0.15)
  combo_kappa <- params$crab_fraction_combo_kappa %||% 8
  a0 <- ifelse(use_osp & osp_flag == 1L, combo_mu * combo_kappa,       set_val * kappa)
  b0 <- ifelse(use_osp & osp_flag == 1L, (1 - combo_mu) * combo_kappa, (1 - set_val) * kappa)

  .say(sprintf(paste0("  Crab fraction f: %d stratum/strata (%s); %d informed by egress classification",
                      " (>= %d boats), %d bounded below by OSP crab-only (>= %d boats), rest on set value %.2f.\n"),
               K, params$crab_fraction_strata %||% "none", sum(use_eg), min_obs,
               sum(use_osp & osp_flag == 1L), osp_min_obs, set_val))
  if (osp_flag == 1L) {
    lb <- ifelse(osp_nt > 0, osp_nc / pmax(osp_nt, 1), NA_real_)
    .say(sprintf("    OSP crab-only lower bound by stratum: %s\n",
                 paste(sprintf("%s=%s", strata,
                               ifelse(is.na(lb), "-", sprintf("%.2f", lb))), collapse = " ")))
    .say(sprintf("    OSP daily observations feeding the bound: %d (beta-binomial, kappa prior centred %.0f)\n",
                 length(ospf$stratum), params$crab_fraction_osp_kappa_prior_mu %||% 20))
    if (sum(use_eg & use_osp) == 0)
      .say(paste0("    NOTE: no stratum has BOTH streams, so the combo-trip share theta rests entirely",
                  " on its prior (crab_fraction_combo_share). The OSP bound constrains f from below;",
                  " its level does not become data-driven until egress classification covers an",
                  " OSP-covered stratum.\n"))
  }

  c(base, list(
    crab_fraction_estimate = 1L,
    crab_fraction_value  = as.array(rep(set_val, K)),
    crab_fraction_alpha0 = as.array(as.numeric(a0)), crab_fraction_beta0 = as.array(as.numeric(b0)),
    crab_fraction_n_total = as.array(as.integer(eg_nt)), crab_fraction_n_crab = as.array(as.integer(eg_nc)),
    osp_crab_lower = osp_flag,
    osp_f_n_total = as.array(as.integer(osp_nt)), osp_f_n_crab = as.array(as.integer(osp_nc)),
    OSPF_n        = length(ospf$stratum),
    osp_f_stratum = as.integer(ospf$stratum),
    osp_f_total   = as.integer(ospf$total),
    osp_f_crab    = as.integer(ospf$crab),
    osp_f_kappa_prior_mu = params$crab_fraction_osp_kappa_prior_mu %||% 20))
}


# Per-DAY point f for the design-based PE: the central value of the BSS f for each day's
# stratum, so the boat PE and boat BSS sit on the same crab-directed basis. Returns
# rep(1, nrow(days)) for shore / when f is off.
#
# With the improvement-8 lower bound, E[f] is no longer a conjugate Beta mean, so it is
# computed by mirroring the Stan posterior on a 1-D grid over theta:
#     f_lower_hat = (1 + osp_crab) / (2 + osp_total)          Beta(1,1) posterior mean
#     p(theta) proportional to Beta(theta; a0, b0) x Binom(n_crab | n_total, f(theta))
#     E[f]        = f_lower_hat + (1 - f_lower_hat) * E[theta]
# Holding f_lower at its posterior mean is a deliberate approximation and a good one:
# the OSP denominators run to hundreds of boats per stratum, so f_lower's own sampling
# error is small next to theta's prior width. With no OSP rows this reduces EXACTLY to
# the previous Beta posterior mean.
crab_fraction_point_day <- function(is_boat, days, params) {
  D <- nrow(days)
  if (!isTRUE(is_boat) || !isTRUE(params$use_crab_fraction)) return(rep(1.0, D))
  cf <- crab_fraction_stan_data(is_shore = FALSE, days = days, params = params, quiet = TRUE)
  if (cf$crab_fraction_estimate == 0L) {
    fk <- as.numeric(cf$crab_fraction_value)                    # pinned per-stratum value
  } else {
    a0 <- as.numeric(cf$crab_fraction_alpha0); b0 <- as.numeric(cf$crab_fraction_beta0)
    nt <- as.numeric(cf$crab_fraction_n_total); nc <- as.numeric(cf$crab_fraction_n_crab)
    ot <- as.numeric(cf$osp_f_n_total);         oc <- as.numeric(cf$osp_f_n_crab)
    lower_on <- identical(as.integer(cf$osp_crab_lower), 1L)
    grid <- seq(1e-4, 1 - 1e-4, length.out = 2001)
    fk <- vapply(seq_along(a0), function(k) {
      fl <- if (lower_on && ot[k] > 0) (1 + oc[k]) / (2 + ot[k]) else 0
      if (nt[k] > 0) {
        f_grid <- fl + (1 - fl) * grid
        lw <- stats::dbeta(grid, a0[k], b0[k], log = TRUE) +
              stats::dbinom(nc[k], nt[k], f_grid, log = TRUE)
        w  <- exp(lw - max(lw)); w <- w / sum(w)
        e_theta <- sum(w * grid)
      } else {
        e_theta <- a0[k] / (a0[k] + b0[k])
      }
      fl + (1 - fl) * e_theta
    }, numeric(1))
  }
  fk[as.integer(cf$f_stratum)]
}
