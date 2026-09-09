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
# crab_fraction.R  (Phase 2 + Phase 3 + improvement 8 + review item 1B: directed-crabbing fraction f)
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
#
# ============================================================================
# REVIEW ITEM 1B (2026-09-08): THE DYNAMIC f
# ============================================================================
# crab_fraction_dynamic = TRUE (the shipped run_config value; the FUNCTION default is FALSE
# so every caller that omits the key gets the legacy construction bit for bit) replaces the
# per-stratum Beta/Binomial above with a random walk on the log-odds of f across the strata
# in chronological order, observed per DAY:
#
#     logit f[k] = eta[k]
#     eta[k]     = f_level_mu + f_level_sd * z[k]                       anchored stratum
#     eta[k]     = eta[prev(k)] + sigma_f * sqrt(gap(k)) * z[k]         otherwise
#     z ~ N(0,1) or Student-t(crab_fraction_walk_df),  sigma_f ~ half-normal(walk_sd_prior)
#
#     sampler contacts, per day i:   crab_i ~ BetaBinomial(total_i, f[k(i)], kappa_I)
#     typed contacts, per day i:     combo_i ~ BetaBinomial(crabbing_i, c[k(i)], kappa_C)   (2026-09-09)
#     OSP crabbing-only, per day j:  osp_crab_j ~ BetaBinomial(osp_total_j, f[k(j)] (1 - c[k(j)]), kappa_O)
#     logit c[k]: the same random walk as f (own step SD sigma_c), the combo-trip share
#                 among crabbing boats, per stratum
#
# For the undergraduate reader: a random walk on the log-odds says "this month's crabbing
# share is probably close to last month's", which lets a month with five contacts borrow
# from its neighbours instead of falling back to a guess; sigma_f is learned from how far
# the well-sampled months actually move. Per-day beta-binomial rather than a monthly
# Binomial matters because the boats contacted on one day are not independent trials
# (weather and which other fisheries are open move all of them together). The two streams
# measure different things: a combo trip is CRABBING to the sampler (crab gear seen) and
# NOT crabbing to OSP (labelled by the other fishery), so the contacts observe f and OSP
# observes f(1 - c). 2026-09-09: the contacts now carry each boat's TRIP TYPE, so c is
# observed directly (combo of the typed crabbing boats, per day) and gets its own walk;
# OSP's column, when it arrives, becomes a check of the shift-time contacts against an
# all-day count rather than the only thing that could identify anything.
#
# THE WALK ORDER. Strata are sorted labels; "month" labels are YYYY-MM (chronological).
# f_walk_prev[k] is the index of the stratum k steps from (0 = anchored on the level
# prior) and f_walk_gap[k] the number of months between them: month -> k-1;
# month_day_type / month_opener -> the same day type's (opener class's) previous month,
# one chain per class; day_type / opener / none -> every stratum anchored, no walk.
#
# WHAT DOES NOT CHANGE. f still enters generated quantities only; both likelihoods are on
# OBSERVED boat counts. crab_fraction_min_obs is NOT applied to the contact stream under
# the dynamic model (every day enters; the walk does the pooling); it remains the PE's
# interpolation trigger in crab_fraction_point_day().
#
# CONFIG (dynamic; defaults shown)
#   crab_fraction_dynamic             FALSE in the function, TRUE in run_config.R
#   crab_fraction_level_sd            1.5   logit SD of the anchored-stratum prior (weak:
#                                           the level is learned; a tight prior at 0.30
#                                           would drag a 0.95 December down by a third)
#   crab_fraction_walk_sd_prior       1.5   half-normal scale of sigma_f (2024-25 monthly
#                                           logit steps have RMS ~1.5, some of it noise)
#   crab_fraction_walk_df             4     Student-t df of the steps (<= 0: Gaussian)
#   crab_fraction_contact_kappa_prior_mu 20 lognormal centre of kappa_I (log-SD 0.75)
#   crab_fraction_combo_level         0.3   centre of the level prior on an anchored stratum's c
#   crab_fraction_combo_level_sd      1.5   its logit SD (weak)
#   crab_fraction_combo_walk_sd_prior 1.5   half-normal scale of sigma_c
#   crab_fraction_combo_kappa_prior_mu 20   lognormal centre of kappa_C (log-SD 0.75)
#   crab_fraction_contact_areas       NULL  = boat_launch_areas; "all" for every private-boat
#                                           interview (the reader applies it, fetch_crab_data.R)
#   crab_fraction_pe_prior_kappa      1     Beta concentration of the PE's per-stratum
#                                           shrinkage (one pseudo-contact at the set value)
###############################################################################

# Months between two YYYY-MM labels (b - a).
.cf_months_between <- function(a, b) {
  ya <- as.integer(substr(a, 1, 4)); ma <- as.integer(substr(a, 6, 7))
  yb <- as.integer(substr(b, 1, 4)); mb <- as.integer(substr(b, 6, 7))
  (yb - ya) * 12L + (mb - ma)
}

# The walk order for a sorted stratum vector (review item 1B). Returns list(prev, gap,
# chain, pos): prev[k] is the predecessor index (0 = anchored), gap[k] the months to it
# (1 when anchored, unused), chain[k] an integer chain id and pos[k] the cumulative month
# position within the chain (both used by the PE's logit interpolation).
crab_fraction_walk_structure <- function(strata, mode = "none") {
  K <- length(strata)
  prev <- rep(0L, K); gap <- rep(1, K); chain <- seq_len(K); pos <- rep(0, K)
  if (K == 0) return(list(prev = prev, gap = gap, chain = chain, pos = pos))
  monthly <- mode %in% c("month", "month_day_type", "month_opener")
  if (!monthly) return(list(prev = prev, gap = gap, chain = chain, pos = pos))
  mo  <- substr(strata, 1, 7)                                   # YYYY-MM
  cls <- if (identical(mode, "month")) rep("all", K) else sub("^.{7}_", "", strata)
  ids <- match(cls, unique(cls))
  for (k in seq_len(K)) {
    earlier <- which(ids == ids[k] & seq_len(K) < k)
    chain[k] <- ids[k]
    if (length(earlier)) {
      p <- max(earlier)                                          # sorted, so the latest earlier one
      prev[k] <- p
      gap[k]  <- max(.cf_months_between(mo[p], mo[k]), 1L)
      pos[k]  <- pos[p] + gap[k]
    }
  }
  list(prev = as.integer(prev), gap = as.numeric(gap), chain = as.integer(chain), pos = as.numeric(pos))
}

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
  # 2026-09-08: year-qualified month labels, so strata sort CHRONOLOGICALLY (the dynamic f
  # of review item 1 walks across strata in order; "%m" alone put December after September
  # and merged the same month of two seasons on a span).
  mo <- format(dates, "%Y-%m")
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
#   params   : run_config, plus params$crab_fraction_rows (the classification rows:
#              sampler contacts and/or egress classification, see
#              crab_fraction_source_rows()) and params$osp_crab_rows (OSP total +
#              crab-only per day), both lifted by the driver from the readers.
# n_f_strata = 1 with no OSP rows reproduces the Phase 2 scalar exactly.
#
# review item 1B: the returned list ALWAYS carries the dynamic-f fields (the Stan data
# block declares them unconditionally); they are inert unless crab_fraction_dynamic = 1
# AND f is estimated. attr(., "f_strata") is a per-stratum audit table (label, walk
# link, contact and OSP sums) that the driver writes beside the posterior f.
crab_fraction_stan_data <- function(is_shore, days, params, quiet = FALSE) {
  D <- nrow(days)
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  apply_cf <- as.integer(!is_shore && isTRUE(params$use_crab_fraction))
  mode <- params$crab_fraction_strata %||% "none"

  labels <- crab_fraction_strata_labels(days$event_date, params)
  strata <- sort(unique(labels))
  K <- length(strata)
  f_stratum <- match(labels, strata)   # 1..K

  .clamp <- function(x) pmin(pmax(as.numeric(x), 1e-4), 1 - 1e-4)
  set_val <- .clamp(params$crab_fraction_set %||% 0.3)
  kappa   <- params$crab_fraction_prior_kappa %||% 20
  min_obs <- params$crab_fraction_min_obs     %||% 20L
  fixed   <- params$crab_fraction_fixed
  dyn_on  <- isTRUE(params$crab_fraction_dynamic %||% FALSE)

  # The dynamic-f fields, inert-valued. Every return path carries them (Stan declares
  # them unconditionally, and bss_assert_stan_data() refuses a list that lacks one).
  c_level <- .clamp(params$crab_fraction_combo_level %||% 0.3)
  .dyn_inert <- function(K) list(
    crab_fraction_dynamic = 0L,
    f_walk_prev = as.array(rep(0L, K)), f_walk_gap = as.array(rep(1, K)),
    f_level_mu = stats::qlogis(set_val),
    f_level_sd = as.numeric(params$crab_fraction_level_sd %||% 1.5),
    f_walk_sd_prior = as.numeric(params$crab_fraction_walk_sd_prior %||% 1.5),
    f_walk_df = as.numeric(params$crab_fraction_walk_df %||% 4),
    CFI_n = 0L, cfi_stratum = integer(0), cfi_total = integer(0), cfi_crab = integer(0),
    cfi_kappa_prior_mu = as.numeric(params$crab_fraction_contact_kappa_prior_mu %||% 20),
    # 2026-09-09: the combo-trip share c, a per-stratum walk observed from the trip types
    combo_dynamic = 0L,
    c_level_mu = stats::qlogis(c_level),
    c_level_sd = as.numeric(params$crab_fraction_combo_level_sd %||% 1.5),
    c_walk_sd_prior = as.numeric(params$crab_fraction_combo_walk_sd_prior %||% 1.5),
    CFC_n = 0L, cfc_stratum = integer(0), cfc_crab = integer(0), cfc_combo = integer(0),
    cfc_kappa_prior_mu = as.numeric(params$crab_fraction_combo_kappa_prior_mu %||% 20))

  neutral <- c(list(
    apply_crab_fraction = 0L, crab_fraction_estimate = 0L,
    n_f_strata = 1L, f_stratum = as.array(rep(1L, D)),
    crab_fraction_value = as.array(1.0), crab_fraction_alpha0 = as.array(1.0),
    crab_fraction_beta0 = as.array(1.0), crab_fraction_n_total = as.array(0L),
    crab_fraction_n_crab = as.array(0L),
    osp_crab_lower = 0L, osp_f_n_total = as.array(0L), osp_f_n_crab = as.array(0L),
    OSPF_n = 0L, osp_f_stratum = integer(0), osp_f_total = integer(0),
    osp_f_crab = integer(0),
    osp_f_kappa_prior_mu = params$crab_fraction_osp_kappa_prior_mu %||% 20),
    .dyn_inert(1L))
  if (apply_cf == 0L) {
    attr(neutral, "f_strata") <- tibble(stratum = 1L, label = "all", walk_prev = 0L, walk_gap = 1,
                                        contact_days = 0L, contacts = 0L, contacts_crabbing = 0L,
                                        typed_days = 0L, typed_crabbing = 0L, typed_combo = 0L,
                                        osp_days = 0L, osp_total = 0L, osp_crab_only = 0L)
    return(neutral)
  }

  base <- list(
    apply_crab_fraction = 1L, n_f_strata = as.integer(K),
    f_stratum = as.array(as.integer(f_stratum)))
  walk <- crab_fraction_walk_structure(strata, mode)

  # Hard set value: pin f per stratum, no uncertainty (sensitivity lever). The OSP bound
  # is meaningless against a pinned f, so it is switched off here; so is the walk.
  if (!is.null(fixed) && length(fixed) == 1L && is.finite(suppressWarnings(as.numeric(fixed)))) {
    fv <- .clamp(fixed)
    .say(sprintf("  Crab fraction f PINNED at %.3f across %d stratum/strata (%s); boat catch x %.3f.\n",
                 fv, K, mode, fv))
    out <- c(base, list(
      crab_fraction_estimate = 0L,
      crab_fraction_value  = as.array(rep(fv, K)),
      crab_fraction_alpha0 = as.array(rep(1.0, K)), crab_fraction_beta0 = as.array(rep(1.0, K)),
      crab_fraction_n_total = as.array(rep(0L, K)), crab_fraction_n_crab = as.array(rep(0L, K)),
      osp_crab_lower = 0L,
      osp_f_n_total = as.array(rep(0L, K)), osp_f_n_crab = as.array(rep(0L, K)),
      OSPF_n = 0L, osp_f_stratum = integer(0), osp_f_total = integer(0),
      osp_f_crab = integer(0),
      osp_f_kappa_prior_mu = params$crab_fraction_osp_kappa_prior_mu %||% 20),
      .dyn_inert(K))
    attr(out, "f_strata") <- tibble(stratum = seq_len(K), label = strata, walk_prev = 0L, walk_gap = 1,
                                    contact_days = 0L, contacts = 0L, contacts_crabbing = 0L,
                                    typed_days = 0L, typed_crabbing = 0L, typed_combo = 0L,
                                    osp_days = 0L, osp_total = 0L, osp_crab_only = 0L)
    return(out)
  }

  # --- (A) the classification rows: crab-vs-total boats, combo trips INCLUDED ---------
  # Both f streams are aggregated over the SAME date window (2026-08-25). Previously
  # the egress rows were pooled season-wide (indeed across seasons) while the OSP rows
  # were restricted to the fit's day set, so under crab_fraction_strata = "none" the two
  # sub-season fits saw the same egress counts but different OSP counts -- the stratum
  # theta was pinned on and the stratum f_lower was bounded by described different date
  # sets. crab_fraction_restrict_to_fit = FALSE restores the old pooling.
  .eg_window <- if (isTRUE(params$crab_fraction_restrict_to_fit %||% TRUE)) days$event_date else NULL
  eg  <- .cf_aggregate(params$crab_fraction_rows, strata, params,
                       num_col = "boats_crabbing", den_col = "boats_total",
                       restrict_to = .eg_window)
  use_eg <- eg$n_total >= min_obs
  # Legacy: strata under the minimum contribute no Binomial. Dynamic: the per-stratum
  # sums are reporting only (every day enters the per-day likelihood), so they are
  # carried ungated.
  eg_nt  <- if (dyn_on) eg$n_total else ifelse(use_eg, eg$n_total, 0L)
  eg_nc  <- if (dyn_on) pmin(eg$n_crab, eg$n_total) else ifelse(use_eg, pmin(eg$n_crab, eg$n_total), 0L)

  # PER-DAY contact rows for the dynamic model's beta-binomial (review item 1B). A row is
  # one source on one day (the interview contacts and the egress classification stay
  # separate rows when both exist: different protocols, different boats).
  cfi <- list(stratum = integer(0), total = integer(0), crab = integer(0))
  cfi_days <- rep(0L, K)
  # 2026-09-09: PER-DAY combo rows for the combo-trip share c: the typed crabbing boats of
  # a day and how many of them were combos (trip_type_class from the interview workbook).
  # A day with typed contacts but no crabbing boat carries no information about c.
  cfc <- list(stratum = integer(0), crab = integer(0), combo = integer(0))
  cfc_days <- rep(0L, K); cfc_crab_k <- rep(0L, K); cfc_combo_k <- rep(0L, K)
  if (dyn_on) {
    rows <- params$crab_fraction_rows
    if (!is.null(rows) && is.data.frame(rows) && nrow(rows) &&
        all(c("event_date", "boats_crabbing", "boats_total") %in% names(rows))) {
      if (!is.null(.eg_window)) rows <- rows[as.Date(rows$event_date) %in% as.Date(.eg_window), , drop = FALSE]
      if (nrow(rows)) {
        rk <- match(crab_fraction_strata_labels(rows$event_date, params), strata)
        tot <- suppressWarnings(as.numeric(rows$boats_total)); cr <- suppressWarnings(as.numeric(rows$boats_crabbing))
        keep <- !is.na(rk) & is.finite(tot) & tot > 0 & is.finite(cr) & cr >= 0
        if (any(keep)) {
          cfi$stratum <- as.integer(rk[keep])
          cfi$total   <- as.integer(round(tot[keep]))
          cfi$crab    <- pmin(as.integer(round(cr[keep])), cfi$total)
          cfi_days    <- as.integer(tabulate(cfi$stratum, nbins = K))
        }
        if (all(c("boats_crab_only", "boats_combo") %in% names(rows))) {
          co  <- suppressWarnings(as.numeric(rows$boats_crab_only)); cb <- suppressWarnings(as.numeric(rows$boats_combo))
          tc  <- co + cb
          keepc <- !is.na(rk) & is.finite(tc) & tc > 0 & is.finite(cb) & cb >= 0
          if (any(keepc)) {
            cfc$stratum <- as.integer(rk[keepc])
            cfc$crab    <- as.integer(round(tc[keepc]))
            cfc$combo   <- pmin(as.integer(round(cb[keepc])), cfc$crab)
            cfc_days    <- as.integer(tabulate(cfc$stratum, nbins = K))
            cfc_crab_k  <- as.integer(vapply(seq_len(K), function(k) sum(cfc$crab[cfc$stratum == k]), numeric(1)))
            cfc_combo_k <- as.integer(vapply(seq_len(K), function(k) sum(cfc$combo[cfc$stratum == k]), numeric(1)))
          }
        }
      }
    }
  }

  # --- (B) OSP crab-only counts: the LOWER bound (improvement 8) / the f(1-c) stream ---
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
  osp_days <- as.integer(tabulate(ospf$stratum, nbins = K))

  # --- theta's prior, per stratum (legacy construction) --------------------------------
  # OSP-informed stratum: theta is the COMBO share (of the not-crab-labelled boats), so
  # it takes the combo prior. Otherwise f = theta and theta takes the ordinary f prior,
  # which is what makes the no-OSP path bit-identical to Phase 2/3.
  combo_mu    <- .clamp(params$crab_fraction_combo_share %||% 0.15)
  combo_kappa <- params$crab_fraction_combo_kappa %||% 8
  a0 <- ifelse(use_osp & osp_flag == 1L, combo_mu * combo_kappa,       set_val * kappa)
  b0 <- ifelse(use_osp & osp_flag == 1L, (1 - combo_mu) * combo_kappa, (1 - set_val) * kappa)

  dyn <- .dyn_inert(K)
  if (dyn_on) {
    dyn$crab_fraction_dynamic <- 1L
    dyn$f_walk_prev <- as.array(walk$prev)
    dyn$f_walk_gap  <- as.array(walk$gap)
    dyn$CFI_n       <- length(cfi$stratum)
    dyn$cfi_stratum <- as.integer(cfi$stratum)
    dyn$cfi_total   <- as.integer(cfi$total)
    dyn$cfi_crab    <- as.integer(cfi$crab)
    # the combo-share walk is live whenever something observes c or f(1 - c): typed
    # contacts, or the OSP crabbing-only stream (where c is then the soft bound's width)
    dyn$combo_dynamic <- as.integer(length(cfc$stratum) > 0 || osp_flag == 1L)
    dyn$CFC_n       <- length(cfc$stratum)
    dyn$cfc_stratum <- as.integer(cfc$stratum)
    dyn$cfc_crab    <- as.integer(cfc$crab)
    dyn$cfc_combo   <- as.integer(cfc$combo)
    n_chains <- length(unique(walk$chain))
    .say(sprintf(paste0("  Crab fraction f (DYNAMIC, review item 1B): %d stratum/strata (%s) on %d walk chain(s);",
                        " %d contact days (%d boats, %d crabbing) enter per day; %d OSP crabbing-only days;",
                        " level prior logit(%.2f) +/- %.2f, sigma_f ~ half-N(%.2f), steps %s.\n"),
                 K, mode, n_chains, dyn$CFI_n, sum(cfi$total), sum(cfi$crab), length(ospf$stratum),
                 set_val, dyn$f_level_sd, dyn$f_walk_sd_prior,
                 if (dyn$f_walk_df > 0) sprintf("Student-t(%g)", dyn$f_walk_df) else "Gaussian"))
    if (K > 1) {
      sh <- ifelse(eg$n_total > 0, sprintf("%.2f (n=%d)", eg$n_crab / pmax(eg$n_total, 1), eg$n_total), "-")
      .say(sprintf("    contact share by stratum: %s\n", paste(sprintf("%s=%s", strata, sh), collapse = " ")))
    }
    if (dyn$CFC_n > 0) {
      cs <- ifelse(cfc_crab_k > 0, sprintf("%.2f (n=%d)", cfc_combo_k / pmax(cfc_crab_k, 1), cfc_crab_k), "-")
      .say(sprintf(paste0("    combo-trip share c (2026-09-09): %d typed contact days, %d crabbing boats, %d combos;",
                          " per-stratum walk, level prior logit(%.2f) +/- %.2f, sigma_c ~ half-N(%.2f); by stratum: %s\n"),
                   dyn$CFC_n, sum(cfc$crab), sum(cfc$combo), c_level, dyn$c_level_sd, dyn$c_walk_sd_prior,
                   paste(sprintf("%s=%s", strata, cs), collapse = " ")))
    } else if (osp_flag == 1L) {
      .say(paste0("    NOTE: no typed contacts in this window, so the combo-trip share c rests on its",
                  " walk prior and the OSP stream is a soft lower bound on f only.\n"))
    }
    if (dyn$CFI_n == 0 && osp_flag == 0L)
      .say("    NOTE: no classification rows in this window; f is its prior walk.\n")
  } else {
    .say(sprintf(paste0("  Crab fraction f: %d stratum/strata (%s); %d informed by egress classification",
                        " (>= %d boats), %d bounded below by OSP crab-only (>= %d boats), rest on set value %.2f.\n"),
                 K, mode, sum(use_eg), min_obs,
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
  }

  out <- c(base, list(
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
    osp_f_kappa_prior_mu = params$crab_fraction_osp_kappa_prior_mu %||% 20),
    dyn)
  attr(out, "f_strata") <- tibble(
    stratum = seq_len(K), label = strata,
    walk_prev = if (dyn_on) walk$prev else rep(0L, K), walk_gap = if (dyn_on) walk$gap else rep(1, K),
    contact_days = cfi_days, contacts = as.integer(eg$n_total), contacts_crabbing = as.integer(pmin(eg$n_crab, eg$n_total)),
    typed_days = cfc_days, typed_crabbing = cfc_crab_k, typed_combo = cfc_combo_k,
    osp_days = osp_days, osp_total = as.integer(osp_nt), osp_crab_only = as.integer(osp_nc))
  out
}


# Per-DAY point f for the design-based PE: the central value of the BSS f for each day's
# stratum, so the boat PE and boat BSS sit on the same crab-directed basis. Returns
# rep(1, nrow(days)) for shore / when f is off.
#
# LEGACY construction (crab_fraction_dynamic = FALSE). With the improvement-8 lower bound,
# E[f] is no longer a conjugate Beta mean, so it is computed by mirroring the Stan
# posterior on a 1-D grid over theta:
#     f_lower_hat = (1 + osp_crab) / (2 + osp_total)          Beta(1,1) posterior mean
#     p(theta) proportional to Beta(theta; a0, b0) x Binom(n_crab | n_total, f(theta))
#     E[f]        = f_lower_hat + (1 - f_lower_hat) * E[theta]
# Holding f_lower at its posterior mean is a deliberate approximation and a good one:
# the OSP denominators run to hundreds of boats per stratum, so f_lower's own sampling
# error is small next to theta's prior width. With no OSP rows this reduces EXACTLY to
# the previous Beta posterior mean.
#
# DYNAMIC construction (review item 1B). The design-based partner of the walk, deliberately
# without the walk's smoothing: a stratum with at least crab_fraction_min_obs classified
# boats takes the conjugate posterior mean (n_crab + a) / (n_total + a + b) under the weak
# Beta(set * kappa_pe, (1 - set) * kappa_pe) prior (kappa_pe = crab_fraction_pe_prior_kappa,
# default 1, one pseudo-contact; the legacy kappa = 20 would drag a 46-of-47 December to
# 0.78); a stratum below the minimum is interpolated on the LOGIT scale between the
# nearest informed strata of its walk chain (end values carried outward); a chain with no
# informed stratum takes the set value. The OSP stream does not enter the PE's f (it
# observes f(1 - c), and c is a model quantity).
crab_fraction_point_day <- function(is_boat, days, params) {
  D <- nrow(days)
  if (!isTRUE(is_boat) || !isTRUE(params$use_crab_fraction)) return(rep(1.0, D))
  cf <- crab_fraction_stan_data(is_shore = FALSE, days = days, params = params, quiet = TRUE)
  if (cf$crab_fraction_estimate == 0L) {
    fk <- as.numeric(cf$crab_fraction_value)                    # pinned per-stratum value
  } else if (identical(as.integer(cf$crab_fraction_dynamic), 1L)) {
    set_val <- pmin(pmax(as.numeric(params$crab_fraction_set %||% 0.3), 1e-4), 1 - 1e-4)
    kpe     <- as.numeric(params$crab_fraction_pe_prior_kappa %||% 1)
    min_obs <- as.numeric(params$crab_fraction_min_obs %||% 20L)
    nt <- as.numeric(cf$crab_fraction_n_total); nc <- as.numeric(cf$crab_fraction_n_crab)
    K  <- length(nt)
    a  <- set_val * kpe; b <- (1 - set_val) * kpe
    informed <- nt >= max(min_obs, 1)
    fk <- ifelse(informed, (nc + a) / (nt + a + b), NA_real_)
    walk <- crab_fraction_walk_structure(attr(cf, "f_strata")$label, params$crab_fraction_strata %||% "none")
    for (ch in unique(walk$chain)) {
      idx <- which(walk$chain == ch)
      inf <- idx[informed[idx]]; un <- idx[!informed[idx]]
      if (!length(un)) next
      if (!length(inf)) { fk[un] <- set_val; next }
      if (length(inf) == 1L) { fk[un] <- fk[inf]; next }
      lg <- stats::approx(x = walk$pos[inf], y = stats::qlogis(pmin(pmax(fk[inf], 1e-4), 1 - 1e-4)),
                          xout = walk$pos[un], rule = 2)$y
      fk[un] <- stats::plogis(lg)
    }
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

# ---------------------------------------------------------------------------
# crab_fraction_source_rows()  (review item 1, 2026-09-08)
#
# Assemble the crabbing-fraction classification rows the model reads
# (params$crab_fraction_rows: event_date, boats_crabbing, boats_total) from the sources
# params$crab_fraction_source names:
#   "interviews"  the sampler contacts (dwg$boat_contacts): every private boat approached
#                 at the launch, crabbing or not. A combo trip counts as CRABBING.
#   "ie"          the WPT/WBL egress classification columns of ingress_egress.xlsx
#                 (attr(ie_data, "crab_fraction_rows")); blank until the pilot delivers.
#   "both"        (default) the two bound together, one row per source per day.
# The two protocols see different boats (contacts happen during sampler shifts; an egress
# survey sees the whole day), so a day with both sources contributes two rows, not a merged
# one; the per-day beta-binomial in the dynamic model absorbs the between-row spread.
# ---------------------------------------------------------------------------
crab_fraction_source_rows <- function(dwg, ie_data, params, quiet = FALSE) {
  src <- tolower(params$crab_fraction_source %||% "both")
  if (!src %in% c("interviews", "ie", "both"))
    stop("params$crab_fraction_source must be interviews | ie | both (got '", src, "')", call. = FALSE)
  empty <- tibble(event_date = as.Date(character()), boats_crabbing = numeric(), boats_total = numeric(), source = character(),
                  boats_typed = numeric(), boats_crab_only = numeric(), boats_combo = numeric())
  from_int <- dwg$boat_contacts
  from_ie  <- attr(ie_data, "crab_fraction_rows")
  .col <- function(df, nm) if (nm %in% names(df)) as.numeric(df[[nm]]) else rep(NA_real_, nrow(df))
  parts <- list()
  if (src %in% c("interviews", "both") && !is.null(from_int) && nrow(from_int))
    parts$interviews <- from_int |>
      transmute(event_date = as.Date(event_date), boats_crabbing = as.numeric(boats_crabbing),
                boats_total = as.numeric(boats_total), source = "interviews",
                # 2026-09-09: the trip-type counts behind the combo-trip share c (typed rows only)
                boats_typed = .col(from_int, "boats_typed"), boats_crab_only = .col(from_int, "boats_crab_only"),
                boats_combo = .col(from_int, "boats_combo"))
  if (src %in% c("ie", "both") && !is.null(from_ie) && is.data.frame(from_ie) && nrow(from_ie))
    parts$ie <- from_ie |>
      transmute(event_date = as.Date(event_date), boats_crabbing = as.numeric(boats_crabbing),
                boats_total = as.numeric(boats_total), source = "ie",
                boats_typed = NA_real_, boats_crab_only = NA_real_, boats_combo = NA_real_) |>
      filter(is.finite(boats_total), boats_total > 0)
  rows <- if (length(parts)) bind_rows(parts) else empty
  if (!isTRUE(quiet)) {
    if (nrow(rows)) {
      by_src <- rows |> group_by(source) |>
        summarise(days = n(), total = sum(boats_total), crab = sum(boats_crabbing), .groups = "drop")
      cat(sprintf("  Crab-fraction classification rows (source = %s): %s\n", src,
                  paste(sprintf("%s %d days, %.0f boats, share %.3f", by_src$source, by_src$days,
                                by_src$total, by_src$crab / pmax(by_src$total, 1)), collapse = "; ")))
      if (any(is.finite(rows$boats_typed)) && sum(rows$boats_typed, na.rm = TRUE) > 0) {
        bc <- rows |> filter(is.finite(boats_typed), boats_typed > 0) |> mutate(m = format(event_date, "%Y-%m")) |>
          group_by(m) |> summarise(cr = sum(boats_crab_only + boats_combo), co = sum(boats_combo), .groups = "drop")
        cat(sprintf("    combo share of crabbing boats (typed contacts): %d combo of %d crabbing; by month: %s\n",
                    sum(bc$co), sum(bc$cr), paste(sprintf("%s %s", bc$m, ifelse(bc$cr > 0, sprintf("%.2f (n=%d)", bc$co / bc$cr, bc$cr), "-")), collapse = ", ")))
      }
      bm <- rows |> mutate(m = format(event_date, "%Y-%m")) |> group_by(m) |>
        summarise(t = sum(boats_total), c = sum(boats_crabbing), .groups = "drop")
      cat("    by month:", paste(sprintf("%s %.2f (n=%.0f)", bm$m, bm$c / pmax(bm$t, 1), bm$t), collapse = ", "), "\n")
    } else cat(sprintf("  Crab-fraction classification rows (source = %s): none; f rests on its prior.\n", src))
  }
  rows
}
