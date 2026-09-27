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
# bss_marine_hazard_covariates.R   (2026-09-25; shared by both drivers)
#
# Small Craft Advisories (and the gale / storm warnings that supersede them) for the Grays
# Harbor Bar and the coastal waters off Westport, and the USCG bar restrictions the samplers
# record on their survey form, as EFFORT covariates for the BSS. They ride on the existing
# K_open / X_open / B_open block (bss_opener_covariates.R): each selected covariate is one
# more column of X_open, so NEITHER STAN MODEL CHANGES and a run with
# marine_hazard_mode = "off" builds exactly the Stan data it built before this file existed.
#
# WHY THESE TWO, AND WHY THEY ARE NOT THE WEATHER MODULE THAT WAS REMOVED (A29)
#   The weather-tide module screened CONTINUOUS weather (wind, wave height, tide phase) and
#   excluded it on its own evidence. An advisory is a different object: it is the categorical
#   product the skipper actually reads at 5 a.m., it is issued for the bar the boats cross,
#   and it is archived to the minute for every day of the season, including the unsampled
#   days the effort process has to interpolate across. On the 2023-11-18 to 2026-09-08
#   launch trailer counts (2026-09-25, offline negative-binomial screen, 1,000 counts), an
#   SCA-or-higher in effect between 04:00 and 16:00 local carried a rate ratio of 0.25
#   [0.21, 0.30] after season, month and weekend/holiday; on the 2,117 dock gear counts 0.74
#   [0.67, 0.82]. Whether that association earns a term in the BSS is what
#   06_diagnostics/run_marine_hazard_batch_2026-09-25.R measures.
#
# THE TWO SOURCES, AND THE ASYMMETRY BETWEEN THEM
#   nws_sca_*          The NWS VTEC archive (04_input_files/nws_marine_hazards.xlsx, built by
#                      04_input_files/build_nws_marine_hazards.R from the Iowa Environmental
#                      Mesonet). Known for EVERY day in the window. Three definitions are
#                      offered: any zone (nws_sca_any), the bar zone only (nws_sca_bar), the
#                      coastal zone only (nws_sca_coastal); one of them per population at most.
#                      Since 2026-09-27 the any-zone flag is also offered split by season
#                      (nws_sca_any_winter + nws_sca_any_rest, THE SEASON SPLIT below), the pair
#                      counting as one definition.
#   bar_restriction    The samplers' "Bar Restrictions" tick in sampler_shifts.xlsx
#                      special_conditions. Observed on SAMPLED days only (and only from the
#                      day the option first appears on the form, detected from the first marine
#                      token in any survey's special_conditions: 2023-11-18 in the committed
#                      workbook; bar_restriction_field_start overrides). On an unsampled day it
#                      is IMPUTED: the expected value from a logistic regression of the observed
#                      tick on the archived advisory flags and a winter indicator, fitted on
#                      every observed day in the workbook (bar_restriction_impute = "nws"), or
#                      the observed rate (= "mean"). A PLUG-IN expected value p in (0, 1) is an
#                      approximation to the mixture the day really is: the effort model applies
#                      exp(B p) where the mixture mean is p exp(B) + (1 - p), and by Jensen the
#                      plug-in UNDERSTATES the mixture mean, by at most 3% for |B| = 0.5, 10%
#                      for |B| = 0.95 (the 2024-25 screen's bar effect) and 17% for |B| = 1.24,
#                      worst near p = 0.6 and vanishing at p = 0 or 1. Stated so a reader knows
#                      the unsampled-day contribution of this term rests on the NWS flags and
#                      on this approximation, not on an observation. Blank special_conditions
#                      on a sampled day is read as "no restriction".
#
# THE SCREEN (marine_hazard_mode = "auto") follows the opener screen in shape (the covariate's
# coefficient after day type and month, on the window's sampled days; Benjamini-Hochberg over
# the marine family, every candidate x population offered; threshold marine_hazard_auto_p on
# the ADJUSTED p) and departs from it in ONE respect: the model is a quasi-Poisson GLM with a
# log link, glm(daily count ~ day_type + month + covariate), not the opener screen's additive
# lm. The reason is scale. An advisory is hypothesized to REDUCE effort by a factor, the Stan
# effort process is log-linear (exp(... + X_open . B_open)), and the advisory days sit in the
# low-count winter months: an additive lm attributes most of a 22-to-1.4 trailer difference
# to the month term and reports p = 0.093 on 2024-25 (0.0006 log-link), while the log-link model reports the rate
# ratio the BSS would fit. The adj_estimate column is therefore a LOG rate ratio and
# rate_ratio = exp(adj_estimate). For bar_restriction the screen uses the OBSERVED tick only,
# never an imputed value. Read bss_opener_covariates.R on why a screen on the fitting data is
# a screen and not a verdict; the arbiter is the paired effort-stream elpd_loo against a
# covariate-free run, and the batch runner computes it.
#
# THE RESPONSE IS THE DAILY SUM OF THE COUNTS, as in diagnose_fishery_spillover(), and that
# choice has a known weakness the screen reports on itself: a day can carry one, two or three
# count sequences, and on 2024-25 the boat carried FEWER on advisory days (1.06 against 1.45
# per day), so part of a daily-sum difference is the number of counts, not the effort. The
# BSS observation model is per COUNT (each count ~ NB2 with that day's mean), so every row
# also carries the same GLM fitted per count with the day's covariates
# (rate_ratio_per_count, adj_p_per_count, n_counts). The auto decision uses the daily-sum
# adj_p, for continuity with the opener screen; the per-count columns are the sensitivity
# check, and on 2024-25 they agree on every decision (boat SCA 0.291 against 0.290, bar 0.325
# against 0.387, shore 0.972 against 0.905, the same three verdicts). A per-count GLM treats
# same-day counts as independent, which overstates its precision, so its p is reported and
# not selected on.
#
# THE SEASON SPLIT (2026-09-27; B41, for D32). The first block cross-validation (Section 1y)
# validated the boat SCA term out of sample on the spring-to-autumn weeks and could not see
# the winter, where the term moves the estimate most; a crude winter-only screen of the
# sampled trailer days put the winter rate ratio near 0.5 to 0.6 against 0.15 to 0.16 for the
# rest of the season and the fitted season-wide 0.31. Two further candidates are therefore
# offered, nws_sca_any_winter (the any-zone flag on days in marine_hazard_winter_months,
# default December to February) and nws_sca_any_rest (the same flag on every other day); they
# sum to nws_sca_any, so they are ONE NWS definition for the one-per-population rule (the pair
# may enter together; either may not enter beside nws_sca_any). They are offered, not in the
# shipped candidate lists, so `auto` and the 2026-09-25 rungs are unchanged; the ladder's M6
# rung enters them by name. December to February because that is the part of the all-gear
# window the OSP counts do not cover (they begin in March), so the "rest" coefficient is the
# one the block cross-validation can test, and because the split at March keeps both terms
# identifiable (52 and 50 flagged calendar days in the 2024-25 all-gear window, 36 and 22 of
# them sampled trailer days; a Dec-to-Mar / Apr-to-Sep split would leave the summer term eight
# sampled advisory days). TWO THINGS TO HOLD when the pair is on. (1) opener_design_matrix()
# drops a column with fewer than opener_min_days (10) flagged days in a fit's window, after
# which that window's winter advisory days are scored as ordinary days under the rest term;
# in the pot-closure window (September to November) that is intended and is the constant
# term for that sub-season, but a short or calm winter in an all-gear window would silently
# do the same, so read the "Effort day covariates" line the prep prints. (2) The winter's
# sampled counts are few and small (2024-25: 36 advisory days of five trailers or fewer), so
# a winter coefficient of the size the multi-season desk screen suggests (rate ratio about
# 0.43, 06_diagnostics/desk_sca_season_split_2026-09-27.R) is likely to sit astride zero in a
# single season's fit; rule 9 of the ladder runner says what that outcome means.
#
# THE CRABBING-FRACTION INTERACTION (boat). Boat effort is ALL private boats. If crabbing boats
# respond to an advisory differently from the finfish boats (a bay crabber does not cross the
# bar; a salmon boat must), f differs on advisory days, and the dynamic monthly f cannot see a
# within-month day effect. Not measured; recorded as an open item, and printed as a note when
# a boat covariate is active.
#
# CONTENTS
#   marine_hazard_events(params)                     read + parse the archive workbook
#   marine_hazard_flag_series(params, start, end)    per-date NWS flags for a window
#   bar_restriction_series(params, dates, nws_day)   the sampler tick, observed + imputed
#   marine_hazard_screen(dwg, flags, params)         the day-type + month adjusted screen
#   marine_hazard_select(screen, params)             off / auto / on / manual -> selection
#   marine_hazard_prepare(dwg, params, output_dir)   orchestration for both drivers
###############################################################################

.mh_labels <- c(
  nws_sca_any     = "NWS SCA or higher, any zone",
  nws_sca_bar     = "NWS SCA or higher, bar zone",
  nws_sca_coastal = "NWS SCA or higher, coastal zone",
  bar_restriction = "USCG bar restriction (sampler tick)",
  nws_sca_any_winter = "NWS SCA or higher, any zone, winter months",
  nws_sca_any_rest   = "NWS SCA or higher, any zone, the other months"
)
# the season split: one NWS definition in two columns (they sum to nws_sca_any)
.mh_split <- c("nws_sca_any_winter", "nws_sca_any_rest")
.mh_series <- c(shore = "Shore gear (effort)", private_boat = "Boat trailers (effort)")

# NULL-only (a zero-length candidate list is a deliberate "none", not a missing key). Used for
# the list-valued keys because the harness's `%||%` treats length-0 as missing and rlang's
# does not; the covariate lists must read the same under both.
.mh_or <- function(a, b) if (is.null(a)) b else a

.mh_defaults <- function(params) list(
  file   = params$marine_hazard_file   %||% "nws_marine_hazards.xlsx",
  sheet  = params$marine_hazard_sheet  %||% "data",
  zones  = params$marine_hazard_zones  %||% c(bar = "PZZ110", coastal = "PZZ156"),
  codes  = params$marine_hazard_codes  %||% c("SC.Y", "RB.Y", "SW.Y", "SI.Y", "GL.W", "SR.W", "SE.W", "HF.W"),
  window = params$marine_hazard_window %||% c(4, 16),
  tz     = params$marine_hazard_tz     %||% "America/Los_Angeles",
  winter = params$marine_hazard_winter_months %||% c(12L, 1L, 2L)   # B41: the months nws_sca_any_winter covers
)

# ---- time helpers -----------------------------------------------------------------------

# ISO-8601 UTC text ("2025-01-03T21:00:00Z"), "yyyy-mm-dd HH:MM:SS" text, or a datetime cell
# -> POSIXct UTC. Anything else NA.
.mh_parse_utc <- function(x) {
  if (inherits(x, "POSIXct")) return(as.POSIXct(format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = "UTC"))
  x <- trimws(as.character(x))
  out <- as.POSIXct(rep(NA_real_, length(x)), origin = "1970-01-01", tz = "UTC")
  iso <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?$", x)
  if (any(iso)) {
    y <- sub("Z$", "", x[iso]); y <- ifelse(nchar(y) == 16, paste0(y, ":00"), y)
    out[iso] <- as.POSIXct(y, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  }
  sp <- !iso & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}(:[0-9]{2})?$", x)
  if (any(sp)) {
    y <- x[sp]; y <- ifelse(nchar(y) == 16, paste0(y, ":00"), y)
    out[sp] <- as.POSIXct(y, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  out
}

# Local wall-clock `hour` (0..24, fractional allowed) on each of `dates`, DST-exact, as UTC.
# Hour 24 is midnight at the START of the next day. A wall-clock time that does not exist on
# a spring-forward day (02:00 to 02:59 in America/Los_Angeles) is moved FORWARD by the gap,
# so 02:30 becomes 03:30 PDT. Platforms disagree on what as.POSIXct() returns for such a
# time (NA on some, one hour EARLIER on glibc), so the result is checked against the wall
# clock that was asked for rather than against NA: any mismatch is rebuilt from local
# midnight (which always exists; DST changes at 02:00) plus the elapsed hours.
.mh_local_time <- function(dates, hour, tz) {
  dates <- as.Date(dates)
  h <- floor(hour + 1e-9); m <- as.integer(round((hour - h) * 60))
  if (m >= 60) { h <- h + 1; m <- 0L }
  d <- dates + (h %/% 24); h <- h %% 24
  want <- sprintf("%02d:%02d", as.integer(h), m)
  s <- sprintf("%s %s:00", format(d, "%Y-%m-%d"), want)
  out <- as.POSIXct(s, tz = tz, format = "%Y-%m-%d %H:%M:%S")
  bad <- is.na(out) | format(out, "%H:%M", tz = tz) != want
  if (any(bad)) {
    midnight <- as.POSIXct(sprintf("%s 00:00:00", format(d[bad], "%Y-%m-%d")), tz = tz, format = "%Y-%m-%d %H:%M:%S")
    out[bad] <- midnight + (h + m / 60) * 3600
  }
  out
}

# For each window [a_i, b_i) is any event [start_j, end_j) in effect? Half-open on both
# sides, so an event that STARTS exactly at b_i, or ENDS exactly at a_i, does not count.
# (NWS times fall on the hour and the default 04:00 / 16:00 PST edges are 12Z / 00Z.)
.mh_any_overlap <- function(a, b, ev_start, ev_end) {
  if (!length(ev_start) || !length(a)) return(rep(0L, length(a)))
  vapply(seq_along(a), function(i) as.integer(any(ev_start < b[i] & ev_end > a[i])), integer(1))
}

# ---- the archive -----------------------------------------------------------------------

# Reads 04_input_files/nws_marine_hazards.xlsx (sheet "data"). Returns a tibble with
# ugc, phenomena, significance, ps ("SC.Y"), start (POSIXct UTC), end, in_effect, and the
# provenance columns when present. attr(., "coverage") is one row per ugc with the pull window
# the builder recorded (pull_start / pull_end, ISO dates); when the workbook carries no pull
# window the coverage is INFERRED from the event dates and attr(., "coverage_inferred") is TRUE,
# which marine_hazard_flag_series() reports rather than trusts silently.
marine_hazard_events <- function(params, path = NULL) {
  mh <- .mh_defaults(params)
  # a name in 04_input_files/, or a path (absolute, or relative to the working directory),
  # the same rule as read_input_workbook()
  path <- path %||% (if (file.exists(mh$file)) mh$file else here::here("04_input_files", mh$file))
  if (!file.exists(path))
    stop("NWS marine hazard archive not found: ", path,
         "\n  Build it with `Rscript 04_input_files/build_nws_marine_hazards.R` (it pulls the",
         " Iowa Environmental Mesonet VTEC archive for the configured zones), point",
         " run_config$marine_hazard_file at it, or set marine_hazard_mode = \"off\".",
         call. = FALSE)
  d <- as.data.frame(read_input_workbook(path, sheet = mh$sheet))
  names(d) <- tolower(trimws(names(d)))
  need <- c("ugc", "phenomena", "significance", "start_utc", "end_utc")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("nws_marine_hazards workbook lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  ev <- tibble::tibble(
    ugc          = toupper(trimws(as.character(d$ugc))),
    phenomena    = toupper(trimws(as.character(d$phenomena))),
    significance = toupper(trimws(as.character(d$significance))),
    eventid      = if ("eventid" %in% names(d)) suppressWarnings(as.integer(d$eventid)) else NA_integer_,
    start        = .mh_parse_utc(d$start_utc),
    end          = .mh_parse_utc(d$end_utc),
    product_id   = if ("product_id" %in% names(d)) as.character(d$product_id) else NA_character_
  ) |>
    dplyr::mutate(ps = paste0(phenomena, ".", significance),
                  in_effect = !is.na(start) & !is.na(end) & end > start)
  n_bad <- sum(is.na(ev$start) | is.na(ev$end))
  if (n_bad) warning(sprintf("nws_marine_hazards: %d row(s) with an unparseable start/end time are ignored", n_bad), call. = FALSE)
  # coverage: the pull window per zone, or inferred from the events
  inferred <- !all(c("pull_start", "pull_end") %in% names(d))
  cov <- if (!inferred) {
    tibble::tibble(ugc = ev$ugc, pull_start = as.Date(substr(as.character(d$pull_start), 1, 10)),
                   pull_end = as.Date(substr(as.character(d$pull_end), 1, 10))) |>
      dplyr::group_by(ugc) |>
      dplyr::summarise(pull_start = min(pull_start, na.rm = TRUE), pull_end = max(pull_end, na.rm = TRUE), .groups = "drop")
  } else {
    ev |> dplyr::filter(!is.na(start), !is.na(end)) |> dplyr::group_by(ugc) |>
      dplyr::summarise(pull_start = as.Date(min(start)), pull_end = as.Date(max(end)), .groups = "drop")
  }
  attr(ev, "coverage") <- cov
  attr(ev, "coverage_inferred") <- inferred
  attr(ev, "path") <- path
  ev
}

# Per-date flags for [date_start, date_end] (default: the estimation window). Each flag is 1
# when an event with one of the configured codes was in effect, in the named zone(s), at any
# moment of the local window (params$marine_hazard_window, hours; c(0, 24) = any time that
# day). Stops when the archive does not cover the window for a configured zone, because a
# covariate that is silently 0 on the uncovered days is worse than no run.
marine_hazard_flag_series <- function(params, date_start = NULL, date_end = NULL, events = NULL,
                                      window = NULL, quiet = TRUE) {
  mh <- .mh_defaults(params)
  window <- window %||% mh$window
  if (length(window) != 2 || !is.numeric(window) || window[1] < 0 || window[2] > 24 || window[2] <= window[1])
    stop("marine_hazard_window must be c(start_hour, end_hour) with 0 <= start < end <= 24, e.g. c(4, 16)", call. = FALSE)
  date_start <- as.Date(date_start %||% params$est_date_start)
  date_end   <- as.Date(date_end   %||% params$est_date_end)
  if (is.na(date_start) || is.na(date_end) || date_end < date_start)
    stop("marine_hazard_flag_series(): a valid date_start / date_end (or est_date_start / est_date_end) is required", call. = FALSE)
  events <- events %||% marine_hazard_events(params)
  zones <- mh$zones
  if (is.null(names(zones)) || any(!nzchar(names(zones))))
    stop("marine_hazard_zones must be a NAMED character vector, e.g. c(bar = \"PZZ110\", coastal = \"PZZ156\")", call. = FALSE)
  zones <- stats::setNames(toupper(trimws(zones)), tolower(names(zones)))

  # coverage: every configured zone must cover the whole window
  cov <- attr(events, "coverage")
  for (z in zones) {
    r <- cov[cov$ugc == z, , drop = FALSE]
    if (!nrow(r) || is.na(r$pull_start) || is.na(r$pull_end) || r$pull_start > date_start || r$pull_end < date_end)
      stop(sprintf(paste0("NWS marine hazard archive does not cover zone %s over %s to %s (it covers %s to %s%s). ",
                          "Re-run `Rscript 04_input_files/build_nws_marine_hazards.R --end <date>` and commit the workbook, ",
                          "or set marine_hazard_mode = \"off\"."),
                   z, date_start, date_end,
                   if (nrow(r)) as.character(r$pull_start) else "?", if (nrow(r)) as.character(r$pull_end) else "?",
                   if (isTRUE(attr(events, "coverage_inferred"))) ", inferred from the events" else ""), call. = FALSE)
  }

  dates <- seq.Date(date_start, date_end, by = "day")
  a <- .mh_local_time(dates, window[1], mh$tz)
  b <- .mh_local_time(dates, window[2], mh$tz)
  ev <- events |> dplyr::filter(in_effect, ps %in% mh$codes)
  flag_for <- function(ugcs) { e <- ev[ev$ugc %in% ugcs, , drop = FALSE]; .mh_any_overlap(a, b, e$start, e$end) }
  wm <- suppressWarnings(as.numeric(mh$winter))
  if (!length(wm) || anyNA(wm) || any(wm != floor(wm)) || any(wm < 1 | wm > 12) || anyDuplicated(wm) || is.logical(mh$winter))
    stop("marine_hazard_winter_months must be distinct whole month numbers in 1..12, e.g. c(12, 1, 2)", call. = FALSE)
  wm <- as.integer(wm)
  is_winter <- as.integer(as.integer(format(dates, "%m")) %in% wm)
  out <- tibble::tibble(
    event_date      = dates,
    nws_sca_any     = flag_for(zones),
    nws_sca_bar     = if ("bar" %in% names(zones)) flag_for(zones[["bar"]]) else NA_integer_,
    nws_sca_coastal = if ("coastal" %in% names(zones)) flag_for(zones[["coastal"]]) else NA_integer_
  )
  # B41: the season split of the any-zone flag; the two columns sum to nws_sca_any
  out$nws_sca_any_winter <- out$nws_sca_any * is_winter
  out$nws_sca_any_rest   <- out$nws_sca_any * (1L - is_winter)
  attr(out, "definition") <- sprintf("%s in effect at any moment between %02d:%02d and %02d:%02d %s; zones %s",
                                     paste(mh$codes, collapse = "/"),
                                     floor(window[1]), round((window[1] %% 1) * 60), floor(window[2]), round((window[2] %% 1) * 60),
                                     mh$tz, paste(sprintf("%s = %s", names(zones), zones), collapse = ", "))
  attr(out, "n_events") <- nrow(ev)
  if (!isTRUE(quiet))
    cat(sprintf("  NWS marine hazards: %d event(s) with codes %s; flags on %d of %d days (any zone; %d in the winter months %s, %d in the rest)\n",
                nrow(ev), paste(mh$codes, collapse = "/"), sum(out$nws_sca_any), nrow(out),
                sum(out$nws_sca_any_winter), paste(wm, collapse = "/"), sum(out$nws_sca_any_rest)))
  out
}

# ---- the samplers' bar-restriction tick ------------------------------------------------

.mh_bar_token   <- "bar restriction"
.mh_form_tokens <- c("bar restriction", "small craft advisory", "gale warning")

# `dates`: the fit window's dates (Date). `nws_day`: per-date flags with window c(0, 24) for
# every date the imputation fit and the prediction need (built here over the archive's
# coverage when NULL). Returns one row per date in `dates`:
#   bar_restriction_obs     1 / 0 on a sampled day (any Grays Harbor survey that day ticked it
#                           / none did), NA on an unsampled day or before the option existed
#   bar_restriction         the covariate: the observation where there is one, else imputed
#   bar_restriction_source  "observed" | "imputed_nws" | "imputed_mean"
# attr(., "fit") is the imputation model's coefficient table (or NULL), attr(., "note") says
# what was done, attr(., "field_start") the first date the option appears on the form.
bar_restriction_series <- function(params, dates, nws_day = NULL, shifts = NULL, events = NULL, quiet = TRUE) {
  dates <- as.Date(dates)
  impute <- tolower(params$bar_restriction_impute %||% "nws")
  if (!impute %in% c("nws", "mean")) stop("bar_restriction_impute must be \"nws\" or \"mean\"", call. = FALSE)
  if (is.null(shifts)) {
    f <- here::here("04_input_files", params$sampler_shifts_file %||% "sampler_shifts.xlsx")
    if (!file.exists(f)) stop("sampler_shifts workbook not found: ", f, " (needed for bar_restriction)", call. = FALSE)
    shifts <- read_input_workbook(f, sheet = params$sampler_shifts_sheet %||% "data")
  }
  need <- c("date", "creel_location", "special_conditions")
  if (!all(need %in% names(shifts))) stop("sampler_shifts lacks column(s): ", paste(setdiff(need, names(shifts)), collapse = ", "), call. = FALSE)
  sh <- tibble::tibble(event_date = as.Date(as.character(shifts$date)),
                       loc = as.character(shifts$creel_location),
                       cond = tolower(trimws(as.character(shifts$special_conditions))))
  sh$cond[is.na(sh$cond)] <- ""
  # the option's first appearance ANYWHERE on the form (any site): before that date the tick
  # could not have been made, so those surveys are not observations of "no restriction"
  field_start <- params$bar_restriction_field_start %||% {
    hit <- sh$event_date[!is.na(sh$event_date) & grepl(paste(.mh_form_tokens, collapse = "|"), sh$cond)]
    if (length(hit)) min(hit) else as.Date(NA)
  }
  field_start <- as.Date(field_start)
  gh <- sh |> dplyr::filter(!is.na(event_date), loc == (params$gh_creel_location %||% "Grays Harbor"))
  if (!is.na(field_start)) gh <- gh |> dplyr::filter(event_date >= field_start)
  obs <- gh |> dplyr::group_by(event_date) |>
    dplyr::summarise(bar_restriction_obs = as.integer(any(grepl(.mh_bar_token, cond))), .groups = "drop")

  out <- tibble::tibble(event_date = dates) |> dplyr::left_join(obs, by = "event_date")
  out$bar_restriction <- as.numeric(out$bar_restriction_obs)
  out$bar_restriction_source <- ifelse(is.na(out$bar_restriction_obs), NA_character_, "observed")
  todo <- is.na(out$bar_restriction_obs)
  fit_tbl <- NULL; note <- character(0)
  rate_window <- if (any(!todo)) mean(out$bar_restriction_obs[!todo]) else NA_real_
  rate_all    <- if (nrow(obs)) mean(obs$bar_restriction_obs) else NA_real_
  fallback_rate <- if (is.finite(rate_window)) rate_window else rate_all

  if (any(todo)) {
    done <- FALSE
    if (impute == "nws") {
      res <- tryCatch({
        # day flags over every date the archive covers, so the imputation model is fitted on
        # EVERY observed day in the workbook (all seasons), not only the window's
        nd <- nws_day
        if (is.null(nd)) {
          events <- events %||% marine_hazard_events(params)
          cov <- attr(events, "coverage")
          mh <- .mh_defaults(params)
          cov <- cov[cov$ugc %in% toupper(mh$zones), , drop = FALSE]
          span <- c(max(c(cov$pull_start, min(obs$event_date, dates))), min(c(cov$pull_end, max(obs$event_date, dates))))
          nd <- marine_hazard_flag_series(params, span[1], span[2], events = events, window = c(0, 24))
        }
        .win <- function(d) as.integer(format(d, "%m") %in% c("10", "11", "12", "01", "02", "03"))
        tr <- obs |> dplyr::inner_join(nd, by = "event_date") |>
          dplyr::mutate(winter = .win(event_date))
        tr <- tr[stats::complete.cases(tr[, c("nws_sca_bar", "nws_sca_coastal")]), , drop = FALSE]
        if (nrow(tr) < 30 || length(unique(tr$bar_restriction_obs)) < 2)
          stop(sprintf("only %d observed day(s) with archive flags, or no variation", nrow(tr)))
        m <- suppressWarnings(stats::glm(bar_restriction_obs ~ nws_sca_bar + nws_sca_coastal + winter,
                                         family = stats::binomial(), data = tr))
        co <- stats::coef(summary(m))
        if (!isTRUE(m$converged) || any(!is.finite(co[, "Estimate"])) || any(abs(co[, "Estimate"]) > 10))
          stop("the logistic imputation did not converge or separated")
        pr <- tibble::tibble(event_date = out$event_date[todo]) |>
          dplyr::left_join(nd, by = "event_date") |> dplyr::mutate(winter = .win(event_date))
        p <- as.numeric(stats::predict(m, newdata = pr, type = "response"))
        list(p = p, coef = tibble::tibble(term = rownames(co), estimate = co[, "Estimate"], se = co[, "Std. Error"],
                                          p = co[, "Pr(>|z|)"], n_days = nrow(tr), n_restricted = sum(tr$bar_restriction_obs)))
      }, error = function(e) e)
      if (inherits(res, "error")) {
        note <- c(note, sprintf("bar_restriction: the NWS logistic imputation was not used (%s); unsampled days take the observed rate", conditionMessage(res)))
      } else {
        p <- res$p; ok <- is.finite(p)
        out$bar_restriction[todo][ok] <- p[ok]
        out$bar_restriction_source[todo][ok] <- "imputed_nws"
        if (any(!ok)) { out$bar_restriction[todo][!ok] <- fallback_rate; out$bar_restriction_source[todo][!ok] <- "imputed_mean" }
        fit_tbl <- res$coef; done <- TRUE
        .est <- function(term) { v <- fit_tbl$estimate[match(term, fit_tbl$term)]; if (length(v) && is.finite(v)) v else NA_real_ }
        note <- c(note, sprintf("bar_restriction: %d unsampled day(s) imputed from logit(P) = %.2f %+.2f*bar SCA %+.2f*coastal SCA %+.2f*winter, fitted on %d observed days (%d restricted)",
                                sum(todo), .est("(Intercept)"), .est("nws_sca_bar"), .est("nws_sca_coastal"), .est("winter"),
                                fit_tbl$n_days[1], fit_tbl$n_restricted[1]))
      }
    }
    if (!done) {
      if (!is.finite(fallback_rate)) stop("bar_restriction: no observed sampler tick exists to impute from", call. = FALSE)
      out$bar_restriction[todo] <- fallback_rate
      out$bar_restriction_source[todo] <- "imputed_mean"
      note <- c(note, sprintf("bar_restriction: %d unsampled day(s) set to the observed rate %.3f", sum(todo), fallback_rate))
    }
  }
  attr(out, "fit") <- fit_tbl
  attr(out, "note") <- note
  attr(out, "field_start") <- field_start
  attr(out, "n_observed") <- sum(!todo)
  if (!isTRUE(quiet)) {
    cat(sprintf("  Bar restrictions: option on the form from %s; %d observed day(s) in the window (%d restricted), %d imputed\n",
                as.character(field_start), sum(!todo), sum(out$bar_restriction_obs %in% 1), sum(todo)))
    for (n in note) cat("   ", n, "\n")
  }
  out
}

# ---- the screen -------------------------------------------------------------------------

# Candidate x population: the covariate's coefficient in a quasi-Poisson (log-link) GLM,
# glm(daily count ~ day_type + month + covariate), on the window's sampled days; the same
# adjustment set as diagnose_fishery_spillover()'s table, on the multiplicative scale the BSS
# uses (see the header). adj_estimate is the log rate ratio; rate_ratio = exp(adj_estimate);
# adj_p is the t-test on the quasi-likelihood dispersion. bar_restriction is screened on its
# OBSERVED days only. n_counts / rate_ratio_per_count / adj_p_per_count are the same GLM per
# count (the sensitivity check described in the header); nothing selects on them.
marine_hazard_screen <- function(dwg, flags, params) {
  cand <- list(shore        = .mh_or(params$marine_hazard_candidates_shore, c("nws_sca_any")),
               private_boat = .mh_or(params$marine_hazard_candidates_boat, c("nws_sca_any", "bar_restriction")))
  daily <- function(df) df |> dplyr::group_by(event_date) |>
    dplyr::summarise(value = sum(count_quantity), .groups = "drop")
  series <- list(shore = daily(dwg$shore_effort), private_boat = daily(dwg$boat_effort))
  raw    <- list(shore = dwg$shore_effort, private_boat = dwg$boat_effort)
  percount <- function(p, col) {
    out <- list(n = NA_integer_, rr = NA_real_, p = NA_real_)
    d <- raw[[p]] |> dplyr::left_join(flags, by = "event_date") |>
      dplyr::mutate(day_type = classify_day_type(event_date, params),
                    month = factor(format(event_date, "%Y-%m")))
    d$flag <- suppressWarnings(as.numeric(d[[col]]))
    d <- d[is.finite(d$count_quantity) & is.finite(d$flag) & d$flag %in% c(0, 1), , drop = FALSE]
    out$n <- nrow(d)
    if (nrow(d) < 10 || length(unique(d$flag)) < 2) return(out)
    fit <- tryCatch(suppressWarnings(stats::glm(count_quantity ~ day_type + month + flag, data = d,
                                                family = stats::quasipoisson(link = "log"))),
                    error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged)) return(out)
    co <- stats::coef(summary(fit))
    if (!"flag" %in% rownames(co) || !is.finite(co["flag", "Std. Error"])) return(out)
    out$rr <- exp(co["flag", "Estimate"]); out$p <- co["flag", "Pr(>|t|)"]
    out
  }
  rows <- list()
  for (p in names(cand)) {
    df <- series[[p]] |> dplyr::left_join(flags, by = "event_date") |>
      dplyr::mutate(day_type = classify_day_type(event_date, params),
                    month = factor(format(event_date, "%Y-%m")))
    for (cc in cand[[p]]) {
      base <- tibble::tibble(population = p, covariate = cc, label = unname(.mh_labels[cc]),
                             n_days = NA_integer_, n_flag = NA_integer_, mean_flag = NA_real_, mean_other = NA_real_,
                             adj_estimate = NA_real_, rate_ratio = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                             adj_p = NA_real_, n_counts = NA_integer_, rate_ratio_per_count = NA_real_,
                             adj_p_per_count = NA_real_, note = "")
      col <- if (identical(cc, "bar_restriction")) "bar_restriction_obs" else cc
      if (!col %in% names(df)) { base$note <- "flag absent"; rows[[length(rows) + 1]] <- base; next }
      d <- df; d$flag <- suppressWarnings(as.numeric(d[[col]]))
      d <- d[is.finite(d$value) & is.finite(d$flag) & d$flag %in% c(0, 1), , drop = FALSE]
      base$n_days <- nrow(d); base$n_flag <- sum(d$flag == 1)
      if (nrow(d) < 10 || length(unique(d$flag)) < 2) { base$note <- "covariate constant in the sampled days"; rows[[length(rows) + 1]] <- base; next }
      base$mean_flag <- mean(d$value[d$flag == 1]); base$mean_other <- mean(d$value[d$flag == 0])
      fit <- tryCatch(suppressWarnings(stats::glm(value ~ day_type + month + flag, data = d,
                                                  family = stats::quasipoisson(link = "log"))),
                      error = function(e) NULL)
      if (is.null(fit) || !isTRUE(fit$converged)) { base$note <- "model failed or did not converge"; rows[[length(rows) + 1]] <- base; next }
      co <- stats::coef(summary(fit))
      if (!"flag" %in% rownames(co) || is.na(co["flag", "Estimate"]) || !is.finite(co["flag", "Std. Error"])) {
        base$note <- "not identified (collinear with day-type/month)"; rows[[length(rows) + 1]] <- base; next }
      est <- co["flag", "Estimate"]; se <- co["flag", "Std. Error"]
      base$adj_estimate <- est; base$rate_ratio <- exp(est); base$se <- se
      base$ci_lo <- est - 1.96 * se; base$ci_hi <- est + 1.96 * se
      base$adj_p <- co["flag", "Pr(>|t|)"]
      # the sensitivity check: the same GLM per COUNT (the BSS observation unit), reported only
      pc <- percount(p, col)
      base$n_counts <- pc$n; base$rate_ratio_per_count <- pc$rr; base$adj_p_per_count <- pc$p
      rows[[length(rows) + 1]] <- base
    }
  }
  dplyr::bind_rows(rows)
}

# off / auto / on / manual -> which covariates enter each population's effort model, with a
# table saying why. `screen` is marine_hazard_screen()'s output (NULL tolerated: auto then
# selects nothing and says so).
marine_hazard_select <- function(screen, params) {
  mode <- tolower(params$marine_hazard_mode %||% "off")
  if (!mode %in% c("off", "auto", "on", "manual"))
    stop("marine_hazard_mode must be one of \"off\", \"auto\", \"on\", \"manual\" (got '", mode, "')", call. = FALSE)
  cand <- list(shore        = .mh_or(params$marine_hazard_candidates_shore, c("nws_sca_any")),
               private_boat = .mh_or(params$marine_hazard_candidates_boat, c("nws_sca_any", "bar_restriction")))
  unknown <- setdiff(unlist(cand), names(.mh_labels))
  if (length(unknown)) stop("unknown marine hazard candidate(s): ", paste(unknown, collapse = ", "),
                            "; the offered flags are ", paste(names(.mh_labels), collapse = ", "), call. = FALSE)
  thresh <- params$marine_hazard_auto_p %||% 0.05
  adj_method <- tolower(params$marine_hazard_auto_p_adjust %||% "BH")
  empty <- list(shore = character(0), private_boat = character(0), table = tibble::tibble(), mode = mode, note = character(0))
  if (identical(mode, "off")) { empty$note <- "marine_hazard_mode = 'off'; no marine hazard covariates."; return(empty) }

  # one NWS definition per population: the three nws_sca_* flags are near-collinear
  # B41: the season-split pair (nws_sca_any_winter + nws_sca_any_rest) is ONE definition and
  # may enter together; it may not enter beside nws_sca_any, which it sums to.
  .one_nws <- function(sel, p_by_name = NULL) {
    nws <- sel[startsWith(sel, "nws_")]
    fam <- ifelse(nws %in% .mh_split, "nws_sca_any_split", nws)
    if (length(unique(fam)) <= 1) return(list(keep = sel, dropped = character(0)))
    fam_p <- vapply(unique(fam), function(f) {
      m <- nws[fam == f]
      if (!is.null(p_by_name) && all(m %in% names(p_by_name)) && any(is.finite(p_by_name[m]))) min(p_by_name[m], na.rm = TRUE) else Inf
    }, numeric(1))
    best <- if (any(is.finite(fam_p))) names(fam_p)[which.min(fam_p)] else unique(fam)[1]
    keep1 <- nws[fam == best]
    list(keep = c(setdiff(sel, nws), keep1), dropped = setdiff(nws, keep1))
  }

  if (mode %in% c("on", "manual")) {
    want <- if (identical(mode, "on")) cand else
      list(shore        = intersect(.mh_or(params$marine_hazard_manual_shore, character(0)), names(.mh_labels)),
           private_boat = intersect(.mh_or(params$marine_hazard_manual_boat,  character(0)), names(.mh_labels)))
    sel <- list(); tbl <- list(); notes <- character(0)
    for (p in names(cand)) {
      w <- .mh_or(want[[p]], character(0))
      drop_nc <- setdiff(w, cand[[p]])
      if (length(drop_nc)) notes <- c(notes, sprintf("marine_hazard_manual_%s names non-candidate(s), ignored: %s",
                                                     if (p == "shore") "shore" else "boat", paste(drop_nc, collapse = ", ")))
      w <- intersect(w, cand[[p]])
      k <- .one_nws(w)
      sel[[p]] <- k$keep
      if (length(cand[[p]]))
        tbl[[p]] <- tibble::tibble(population = p, covariate = cand[[p]], label = unname(.mh_labels[cand[[p]]]),
                                   adj_estimate = NA_real_, p_raw = NA_real_, p_adj = NA_real_,
                                   selected = cand[[p]] %in% k$keep,
                                   reason = dplyr::case_when(cand[[p]] %in% k$keep ~ mode,
                                                             cand[[p]] %in% k$dropped ~ sprintf("one NWS definition per population (kept %s)", paste(k$keep[startsWith(k$keep, "nws_")], collapse = " + ")),
                                                             TRUE ~ "not named"))
      if (length(k$dropped)) notes <- c(notes, sprintf("%s: %s dropped, one NWS definition per population", p, paste(k$dropped, collapse = ", ")))
    }
    if (!is.null(screen) && nrow(screen)) {
      sc <- screen |> dplyr::select(population, covariate, adj_estimate, rate_ratio, p_raw = adj_p,
                                    dplyr::any_of("rate_ratio_per_count"))
      out_tbl <- dplyr::bind_rows(tbl) |> dplyr::select(-adj_estimate, -p_raw) |> dplyr::left_join(sc, by = c("population", "covariate"))
    } else out_tbl <- dplyr::bind_rows(tbl)
    # the same columns whatever the mode and whether the screen ran, so the report chunk and
    # the runner never have to ask which shape they were handed
    for (cc in c("adj_estimate", "rate_ratio", "p_raw", "p_adj", "rate_ratio_per_count"))
      if (nrow(out_tbl) && !cc %in% names(out_tbl)) out_tbl[[cc]] <- NA_real_
    return(list(shore = sel$shore, private_boat = sel$private_boat, table = out_tbl, mode = mode,
                note = c(sprintf("Marine hazard covariates forced (%s); no significance screen applied.", mode), notes)))
  }

  # --- auto -----------------------------------------------------------------
  if (is.null(screen) || !nrow(screen)) {
    empty$note <- "marine_hazard_mode = 'auto' but the screen produced no table; no marine hazard covariates selected."
    return(empty)
  }
  fam <- screen |> dplyr::filter(population %in% names(cand)) |>
    dplyr::rowwise() |> dplyr::filter(covariate %in% cand[[population]]) |> dplyr::ungroup()
  p_adj_method <- switch(adj_method, bh = "BH", bonferroni = "bonferroni", none = "none", adj_method)
  n_family <- nrow(fam)   # pinned to the family OFFERED, identified or not (see opener_select)
  fam$p_adj <- if (identical(p_adj_method, "none")) fam$adj_p else stats::p.adjust(fam$adj_p, method = p_adj_method, n = n_family)
  fam <- fam |> dplyr::mutate(selected = is.finite(p_adj) & p_adj < thresh,
                              reason = dplyr::case_when(!is.finite(p_adj) ~ paste0("no adjusted p (", note, ")"),
                                                        p_adj < thresh ~ sprintf("adjusted p %.4g < %.3g", p_adj, thresh),
                                                        TRUE ~ sprintf("adjusted p %.4g >= %.3g", p_adj, thresh)))
  sel <- list(); notes <- character(0)
  for (p in names(cand)) {
    rows <- fam[fam$population == p, , drop = FALSE]
    s <- rows$covariate[rows$selected]
    k <- .one_nws(s, stats::setNames(rows$p_adj, rows$covariate))
    sel[[p]] <- k$keep
    if (length(k$dropped)) {
      fam$selected[fam$population == p & fam$covariate %in% k$dropped] <- FALSE
      fam$reason[fam$population == p & fam$covariate %in% k$dropped] <-
        sprintf("cleared the screen but dropped: one NWS definition per population (kept %s)", paste(k$keep[startsWith(k$keep, "nws_")], collapse = " + "))
      notes <- c(notes, sprintf("%s: %s dropped, one NWS definition per population", p, paste(k$dropped, collapse = ", ")))
    }
  }
  if (!"rate_ratio_per_count" %in% names(fam)) fam$rate_ratio_per_count <- NA_real_   # a screen table from elsewhere
  tbl <- fam |> dplyr::transmute(population, covariate, label, adj_estimate, rate_ratio, p_raw = adj_p, p_adj, selected, reason,
                                 rate_ratio_per_count)
  note <- c(sprintf(paste("Auto marine hazard screen: %d effort test(s) in the multiplicity family (%d with an identifiable",
                          "estimate), adjustment '%s', threshold %.3g on the ADJUSTED p. A screen, not a verdict:",
                          "confirm a selected term against a covariate-free run's paired effort elpd_loo."),
                    n_family, sum(is.finite(fam$adj_p)), p_adj_method, thresh), notes)
  list(shore = sel$shore, private_boat = sel$private_boat, table = tbl, mode = mode, note = note)
}

# ---- orchestration ----------------------------------------------------------------------

# Builds the flags, runs the screen, applies the mode, and installs the result on `params`:
#   params$marine_hazard_selected   list(shore = chr, private_boat = chr), read by both preps
#   params$opener_flags             the per-date flag table the preps' opener_design_matrix()
#                                   reads, with the marine columns joined on (the opener columns
#                                   are untouched; a NULL opener table becomes the marine table)
# With marine_hazard_mode = "off" nothing is read and params changes only by the empty
# selection, so production is untouched. Any other mode STOPS on a missing or non-covering
# archive rather than running without the covariate the user asked for.
marine_hazard_prepare <- function(dwg, params, output_dir = NULL, quiet = FALSE) {
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  mode <- tolower(params$marine_hazard_mode %||% "off")
  empty_sel <- list(shore = character(0), private_boat = character(0))
  if (identical(mode, "off")) {
    params$marine_hazard_selected <- empty_sel
    .say("Marine hazard effort covariates: mode=off (nothing read; the K_open block is unchanged)\n")
    return(invisible(list(params = params, sel = c(empty_sel, list(table = tibble::tibble(), mode = mode,
                                                                    note = "marine_hazard_mode = 'off'")),
                          flags = NULL, screen = NULL, active = FALSE)))
  }
  events <- marine_hazard_events(params)
  cand_all <- unique(c(.mh_or(params$marine_hazard_candidates_shore, c("nws_sca_any")),
                       .mh_or(params$marine_hazard_candidates_boat,  c("nws_sca_any", "bar_restriction")),
                       .mh_or(params$marine_hazard_manual_shore, character(0)),
                       .mh_or(params$marine_hazard_manual_boat,  character(0))))
  nws <- marine_hazard_flag_series(params, events = events, quiet = quiet)
  flags <- nws
  bar <- NULL
  if ("bar_restriction" %in% cand_all) {
    bar <- bar_restriction_series(params, nws$event_date, events = events, quiet = quiet)
    flags <- flags |> dplyr::left_join(bar, by = "event_date")
  }
  screen <- tryCatch(marine_hazard_screen(dwg, flags, params),
                     error = function(e) { .say("  Marine hazard screen failed:", conditionMessage(e), "\n"); NULL })
  sel <- marine_hazard_select(screen, params)
  params$marine_hazard_selected <- list(shore = sel$shore, private_boat = sel$private_boat)
  # the preps read ONE per-date flag table; join the marine columns onto the opener columns
  mcols <- flags |> dplyr::select(event_date, dplyr::any_of(c("nws_sca_any", "nws_sca_bar", "nws_sca_coastal", "bar_restriction", .mh_split)))
  of <- params$opener_flags
  params$opener_flags <- if (is.null(of) || !nrow(of)) mcols else
    dplyr::full_join(of |> dplyr::select(-dplyr::any_of(names(mcols)[-1])), mcols, by = "event_date") |> dplyr::arrange(event_date)

  .say(sprintf("Marine hazard effort covariates: mode=%s\n  definition: %s\n  shore: %s\n  private boat: %s\n",
               sel$mode, attr(nws, "definition"),
               if (length(sel$shore)) paste(sel$shore, collapse = ", ") else "(none)",
               if (length(sel$private_boat)) paste(sel$private_boat, collapse = ", ") else "(none)"))
  for (n in sel$note) .say("  ", n, "\n", sep = "")
  if (length(sel$private_boat) && isTRUE(params$use_crab_fraction))
    .say(paste0("  NOTE: a BOAT marine hazard covariate is active. Boat effort is ALL private boats; if crabbing boats\n",
                "        respond to an advisory or a bar restriction differently from finfish boats, the crabbing\n",
                "        fraction f differs on those days and the monthly f cannot see it. Read the boat total with that\n",
                "        in mind (CHANGE_REGISTER D30).\n"))
  if (!is.null(output_dir) && dir.exists(output_dir)) {
    tryCatch({
      utils::write.csv(flags, file.path(output_dir, "marine_hazard_flags.csv"), row.names = FALSE)
      if (!is.null(screen) && nrow(screen)) utils::write.csv(screen, file.path(output_dir, "marine_hazard_screen.csv"), row.names = FALSE)
      if (!is.null(sel$table) && nrow(sel$table)) utils::write.csv(sel$table, file.path(output_dir, "marine_hazard_selection.csv"), row.names = FALSE)
      if (!is.null(bar) && !is.null(attr(bar, "fit"))) utils::write.csv(attr(bar, "fit"), file.path(output_dir, "marine_hazard_bar_imputation.csv"), row.names = FALSE)
    }, error = function(e) .say("  (marine hazard CSVs not written: ", conditionMessage(e), ")\n"))
  }
  invisible(list(params = params, sel = sel, flags = flags, screen = screen, bar = bar,
                 definition = attr(nws, "definition"), active = length(c(sel$shore, sel$private_boat)) > 0))
}
