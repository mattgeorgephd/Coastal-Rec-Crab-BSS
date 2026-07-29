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
# crab_fraction.R  (Phase 2 + Phase 3: directed-crabbing fraction f)
#
# f = the share of private boats at the launch that are crabbing (vs. targeting
# other fisheries). The effort series (trailer counts, OSP boat totals) count ALL
# private boats, so f converts all-boat effort to crab effort. Without f the model
# implicitly sets f = 1 (every boat crabbing), which biases the boat catch high.
#
# HYBRID DESIGN (matches R_G / tau_boat / L_effective): f is a Beta-prior parameter
# centered on a configurable SET VALUE, optionally updated by the crab-creel I/E
# crab-vs-total boat classification, applied ONLY in the Stan generated quantities
# (scaling boat E and C), so it never perturbs effort/CPUE sampling.
#
# PHASE 3: f is PER STRATUM. crab_fraction_strata in run_config selects the
# stratification ("none" = one scalar f = the Phase 2 behavior; "month";
# "day_type"; "month_day_type"). Each stratum gets its own Beta prior (centered on
# the set value) and its own I/E classification counts, so a season with a lower
# summer crab share and a higher shoulder share is captured instead of averaged.
# n_f_strata = 1 reproduces the Phase 2 scalar exactly.
#
# The I/E classification reaches here via params$crab_fraction_rows (the driver
# lifts attr(ie_data, "crab_fraction_rows"), emitted by fetch_ie_data, into params
# right after fetch_ie_data), a per-WBL-row tibble (event_date, boats_crabbing,
# boats_total). Empty (the pilot is in progress) -> every stratum uses the set value.
#
# CONFIG (all optional; defaults shown)
#   use_crab_fraction        FALSE
#   crab_fraction_strata     "none"  ("none" | "month" | "day_type" | "month_day_type")
#   crab_fraction_set        0.3     (set value = Beta prior mean, and the fallback)
#   crab_fraction_prior_kappa 20     (Beta concentration)
#   crab_fraction_fixed      NA      (a number PINS f exactly, no uncertainty)
#   crab_fraction_min_obs    20      (min classified boats in a stratum before its Binomial binds)
###############################################################################

# Per-date stratum labels, derived purely from the DATE + config (so the day set
# and the I/E classification rows are labeled by one consistent rule).
crab_fraction_strata_labels <- function(dates, params) {
  mode  <- params$crab_fraction_strata %||% "none"
  dates <- as.Date(dates)
  if (identical(mode, "none")) return(rep("all", length(dates)))
  mo <- format(dates, "%m")
  wknd_days <- params$days_wkend %||% c("Friday", "Saturday", "Sunday")
  dt <- ifelse(weekdays(dates) %in% wknd_days, "wknd", "wkdy")
  hol <- params$crabbing_holiday_dates
  if (!is.null(hol)) dt[dates %in% as.Date(hol)] <- "wknd"   # holidays typed as weekend
  switch(mode,
    "month"          = mo,
    "day_type"       = dt,
    "month_day_type" = paste(mo, dt, sep = "_"),
    stop("params$crab_fraction_strata must be none|month|day_type|month_day_type (got '",
         mode, "')", call. = FALSE))
}

# Build the Stan data for the (possibly per-stratum) crabbing fraction f.
#   is_shore : shore fits get apply_crab_fraction = 0 (f pinned to 1)
#   days     : the fit's day set (needs event_date); drives n_f_strata + f_stratum
#   params   : run_config (+ params$crab_fraction_rows lifted by the driver)
# Returns the Stan data entries; n_f_strata = 1 reproduces the Phase 2 scalar.
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
    crab_fraction_n_crab = as.array(0L))
  if (apply_cf == 0L) return(neutral)

  .clamp <- function(x) pmin(pmax(as.numeric(x), 1e-4), 1 - 1e-4)
  set_val <- .clamp(params$crab_fraction_set %||% 0.3)
  kappa   <- params$crab_fraction_prior_kappa %||% 20
  min_obs <- params$crab_fraction_min_obs     %||% 20L
  fixed   <- params$crab_fraction_fixed

  base <- list(
    apply_crab_fraction = 1L, n_f_strata = as.integer(K),
    f_stratum = as.array(as.integer(f_stratum)))

  # Hard set value: pin f per stratum, no uncertainty (sensitivity lever).
  if (!is.null(fixed) && length(fixed) == 1L && is.finite(suppressWarnings(as.numeric(fixed)))) {
    fv <- .clamp(fixed)
    .say(sprintf("  Crab fraction f PINNED at %.3f across %d stratum/strata (%s); boat catch x %.3f.\n",
                 fv, K, params$crab_fraction_strata %||% "none", fv))
    return(c(base, list(
      crab_fraction_estimate = 0L,
      crab_fraction_value  = as.array(rep(fv, K)),
      crab_fraction_alpha0 = as.array(rep(1.0, K)), crab_fraction_beta0 = as.array(rep(1.0, K)),
      crab_fraction_n_total = as.array(rep(0L, K)), crab_fraction_n_crab = as.array(rep(0L, K)))))
  }

  # Per-stratum I/E classification counts from the WBL classification rows.
  n_total <- rep(0L, K); n_crab <- rep(0L, K)
  rows <- params$crab_fraction_rows
  if (!is.null(rows) && is.data.frame(rows) && nrow(rows) > 0 &&
      all(c("boats_crabbing", "boats_total") %in% names(rows))) {
    rlab <- crab_fraction_strata_labels(rows$event_date, params)
    for (k in seq_len(K)) {
      sel <- rlab == strata[k]
      if (any(sel)) {
        n_total[k] <- as.integer(round(sum(rows$boats_total[sel],    na.rm = TRUE)))
        n_crab[k]  <- as.integer(round(sum(rows$boats_crabbing[sel], na.rm = TRUE)))
        if (n_crab[k] > n_total[k]) n_crab[k] <- n_total[k]
      }
    }
  }

  a0 <- rep(set_val * kappa, K); b0 <- rep((1 - set_val) * kappa, K)
  use <- n_total >= min_obs
  nt  <- ifelse(use, n_total, 0L)
  nc  <- ifelse(use, pmin(n_crab, n_total), 0L)
  .say(sprintf("  Crab fraction f: %d stratum/strata (%s); %d informed by I/E (>= %d classified boats), rest on set value %.2f.\n",
               K, params$crab_fraction_strata %||% "none", sum(use), min_obs, set_val))

  c(base, list(
    crab_fraction_estimate = 1L,
    crab_fraction_value  = as.array(rep(set_val, K)),
    crab_fraction_alpha0 = as.array(a0), crab_fraction_beta0 = as.array(b0),
    crab_fraction_n_total = as.array(as.integer(nt)), crab_fraction_n_crab = as.array(as.integer(nc))))
}

# Per-DAY point f for the design-based PE: the central value (Beta posterior mean)
# of the BSS f for each day's stratum, so the boat PE and boat BSS sit on the same
# crab-directed basis. Returns rep(1, nrow(days)) for shore / when f is off.
crab_fraction_point_day <- function(is_boat, days, params) {
  D <- nrow(days)
  if (!isTRUE(is_boat) || !isTRUE(params$use_crab_fraction)) return(rep(1.0, D))
  cf <- crab_fraction_stan_data(is_shore = FALSE, days = days, params = params, quiet = TRUE)
  if (cf$crab_fraction_estimate == 0L) {
    fk <- as.numeric(cf$crab_fraction_value)                    # pinned per-stratum value
  } else {
    a0 <- as.numeric(cf$crab_fraction_alpha0); b0 <- as.numeric(cf$crab_fraction_beta0)
    nt <- as.numeric(cf$crab_fraction_n_total); nc <- as.numeric(cf$crab_fraction_n_crab)
    fk <- (a0 + nc) / (a0 + b0 + nt)                            # Beta posterior mean = BSS central f
  }
  fk[as.integer(cf$f_stratum)]
}
