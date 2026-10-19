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
# fetch_osp_boat_counts.R  (Phase 0: OSP boat-count ingest)
#
# Reads the WDFW Ocean Sampling Program (OSP) daily private-boat counts for the
# Westport boat launch (WBL) and returns a tidy per-day series shaped exactly
# like the trailer `boat_effort` object from fetch_crab_data(), so Phase 1 can
# add it as a SECOND boat effort-observation stream on the same latent process.
#
# WHAT THE FIELD IS
#   `WestportPrivateEffort` is a DAILY BOAT TOTAL: the count of ALL private boats
#   using the launch that day, across every fishery (confirmed with OSP). It is
#   NOT crab-specific (it needs the crabbing fraction f to become crab effort;
#   that is Phase 2) and it is NOT an instantaneous snapshot like the trailer
#   count (it is a full-day aggregate). OSP samplers are on the docks only during
#   the bottomfish season (~mid-March to mid-October), so the series is seasonal.
#
# NS vs ZERO (the rule that must not be gotten wrong)
#   - A date that is ABSENT from the file is NON-SAMPLED (NS): no observation, a
#     latent day the AR process imputes. We therefore never fabricate a 0 for a
#     missing date (in particular, OSP-dark winter days stay absent, not zero).
#   - A row with value 0 is an OBSERVED no-effort day: real data, kept as a 0
#     count (neg_binomial_2 handles it).
#   This function only ever returns OBSERVED rows (including observed zeros); the
#   model imputes everything absent.
#
# DE-DUPLICATION (defensive)
#   The 2024-25 delivery had every 2025 date duplicated with identical values (a
#   copy artifact). A date with >1 row is collapsed. If the duplicate values
#   AGREE it is lossless; if they DISAGREE they are treated as genuine sub-daily
#   counts and combined by params$osp_dupe_resolve (default "mean"), never
#   silently dropped.
#
# IMPROVEMENT 8 (2026-08-25): OPTIONAL CRAB-ONLY COLUMN
#   If the workbook carries params$osp_crab_only_col, this reader also emits a per-day
#   table of (osp_total, osp_crab_only) as attr(<result>, "osp_crab_rows"), which the
#   crabbing-fraction helper turns into a LOWER BOUND on f: a hard bound on the legacy f path,
#   and under the shipped DYNAMIC f a soft beta-binomial likelihood on f(1 - c), the crab-only
#   share (corrected 2026-09-28; see crab_fraction.R). The column is a count of
#   boats OSP labelled as crabbing ONLY: OSP does not record combo trips, so a boat that
#   crabs AND fishes something else is labelled by the OTHER fishery. That is why the
#   count bounds f from below rather than estimating it; see crab_fraction.R. The column
#   does not exist yet, so today the attribute is an empty tibble and every downstream f
#   feature is inert.
#
#   Rows where the crab-only count exceeds the daily total are a data error, not a combo
#   effect. They are clamped and COUNTED in the log rather than silently truncated.
#
# CONFIG (all optional; historical Westport defaults via %||%)
#   osp_boat_counts_file  = "WBL_boat_counts.xlsx"
#   osp_boat_counts_sheet = "Sheet1"
#   osp_effort_col        = "WestportPrivateEffort"
#   osp_crab_only_col     = c("WPTPrivateCrabOnly", "crabbing_only", "WestportCrabOnlyEffort")   # optional; first present wins
#   osp_crab_checked_col  = "WestportCrabClassified"      # optional; private boats SAMPLED that day
#   osp_sample_rate_col   = "WestportPrivateSampleRate"   # optional; the day's sampling rate (0.5 or 50)
#   osp_sampling_rate_source = "auto"    # auto | column | schedule | none (osp_sampling_rates.R)
#   osp_sampling_rates_file  = "osp_sampling_rates.xlsx"  # the manual's minimum-rate schedule
#   osp_crab_only_unit    = "count"      # count | fraction (of the boats sampled)
#   osp_crab_only_basis   = "sampled"    # sampled | expanded (already divided by the rate)
#   osp_absent_day        = "unclassified"  # unclassified | zero: a day with boats and NO crab-only value (D39)
#   osp_dupe_resolve      = "mean"    # mean | sum | max | first (only on DISAGREEING dupes)
#   est_date_start / est_date_end     # the estimation window (reused, already in run_config)
#
# RETURNS a tibble with the boat_effort schema:
#   event_date, count_sequence, count_quantity, section_num, count_type, population
# plus attr(., "osp_crab_rows") = tibble(event_date, osp_total (the boats SAMPLED, the binomial n),
#   osp_crab_only, osp_boat_total, osp_sample_rate, osp_rate_source) or empty; the drivers write it
#   as osp_crab_only_daily.csv; and attr(., "osp_absent_days") = tibble(event_date, osp_boat_total,
#   read_as): the days with boats and no crab-only value, how osp_absent_day read them (D39); the
#   drivers write it as osp_absent_days.csv.
#
# Requires: dplyr, tibble, readxl, here.
###############################################################################

fetch_osp_boat_counts <- function(params) {
  cat("\n  Reading OSP Westport boat-launch counts...\n")
  # D41 (2026-10-19): what a reader of the run must see goes through bss_warn() when the run-warnings
  # log is loaded (the drivers source the whole function library, so the note reaches the report's
  # last section and run_warnings.csv), and to the console in any case (a standalone call, the
  # harness). Until then the schedule-minimum note was a console line inside a hidden chunk.
  .osp_note <- function(msg) {
    if (exists("bss_warn", mode = "function")) bss_warn("OSP", msg, severity = "note")
    cat("  NOTE:", msg, "\n")
  }

  osp_file <- here::here("04_input_files", params$osp_boat_counts_file %||% "WBL_boat_counts.xlsx")
  empty <- tibble::tibble(
    event_date = as.Date(character()), count_sequence = integer(),
    count_quantity = numeric(), section_num = integer(),
    count_type = character(), population = character())

  if (!file.exists(osp_file)) {
    cat("  WARNING: OSP boat-count file not found at", osp_file, "\n")
    return(empty)
  }

  val_col  <- params$osp_effort_col    %||% "WestportPrivateEffort"
  raw <- read_input_workbook(osp_file, sheet = params$osp_boat_counts_sheet %||% "Sheet1")
  # B48 (2026-09-28): the crab-only column may be named by any of several candidates, the first
  # present wins. OSP will deliver it as "crabbing_only" (Matt, 2026-09-28); the older
  # placeholder name stays accepted.
  crab_cands <- params$osp_crab_only_col %||% c("WPTPrivateCrabOnly", "crabbing_only", "WestportCrabOnlyEffort")
  crab_col <- intersect(crab_cands, names(raw))[1]
  if (is.na(crab_col)) crab_col <- crab_cands[1]

  # Build event_date from Year/Month/Day, or from a pre-existing ISO `date` column.
  if (all(c("Year", "Month", "Day") %in% names(raw))) {
    raw <- raw |> dplyr::mutate(event_date = as.Date(sprintf(
      "%04d-%02d-%02d", as.integer(Year), as.integer(Month), as.integer(Day))))
  } else if ("date" %in% names(raw)) {
    raw <- raw |> dplyr::mutate(event_date = as.Date(date))
  } else {
    stop("fetch_osp_boat_counts: need Year/Month/Day or a `date` column in ", osp_file, call. = FALSE)
  }
  if (!val_col %in% names(raw))
    stop("fetch_osp_boat_counts: value column '", val_col, "' not found in ", osp_file, call. = FALSE)

  has_crab_col <- crab_col %in% names(raw)
  # 2026-09-28 (B46): OSP's crab-only data are awaited with their sampling frequency (whether
  # OSP classifies every boat or every Nth on busy days). Two settings are ready for it:
  #   osp_crab_checked_col  the number of boats OSP CLASSIFIED that day, when it is not every
  #                         boat; it becomes the binomial n of the crab-only share (the crab-only
  #                         count is out of the boats classified, not out of all boats). Absent:
  #                         every returning boat is taken as classified, as before.
  #   osp_crab_only_unit    "count" (boats; the default) or "fraction" (a 0-1 share of the boats
  #                         classified, converted to a count). A column of values all in [0, 1],
  #                         some fractional, against totals above 1 under "count" STOPS: rounding
  #                         a share to a count would read almost every day as 0 crab-only boats.
  chk_col  <- params$osp_crab_checked_col %||% "WestportCrabClassified"
  has_chk  <- chk_col %in% names(raw)
  # B48 (2026-09-28): OSP samples every k-th private boat at a rate fixed for the day (its
  # sampling manual), so the crab-only count is out of the boats SAMPLED. The day's rate, when
  # delivered, is read here; osp_resolve_sample_n() (osp_sampling_rates.R) turns it, the sampled
  # count above or the manual's schedule into the binomial n. osp_crab_only_basis says whether
  # the crab-only column counts SAMPLED boats ("sampled", the default) or has already been
  # expanded to all boats ("expanded", divided by the rate), in which case it is converted back.
  rate_col  <- params$osp_sample_rate_col %||% "WestportPrivateSampleRate"
  has_rate  <- rate_col %in% names(raw)
  crab_basis <- params$osp_crab_only_basis %||% "sampled"
  if (!crab_basis %in% c("sampled", "expanded"))
    stop("fetch_osp_boat_counts: osp_crab_only_basis must be \"sampled\" or \"expanded\" (got '", crab_basis, "').", call. = FALSE)
  crab_unit <- params$osp_crab_only_unit %||% "count"
  if (!crab_unit %in% c("count", "fraction"))
    stop("fetch_osp_boat_counts: osp_crab_only_unit must be \"count\" or \"fraction\" (got '", crab_unit, "').", call. = FALSE)
  osp <- raw |>
    dplyr::mutate(osp_boat_total = suppressWarnings(as.numeric(.data[[val_col]])),
                  osp_crab_only  = if (has_crab_col) suppressWarnings(as.numeric(.data[[crab_col]]))
                                   else NA_real_,
                  osp_checked    = if (has_chk) suppressWarnings(as.numeric(.data[[chk_col]])) else NA_real_,
                  osp_rate       = if (has_rate) osp_rate_as_fraction(.data[[rate_col]], rate_col) else NA_real_) |>
    dplyr::filter(!is.na(event_date), !is.na(osp_boat_total), osp_boat_total >= 0)  # keep observed zeros

  # --- Defensive de-duplication ---
  dup_dates <- osp |> dplyr::count(event_date) |> dplyr::filter(n > 1) |> dplyr::pull(event_date)
  if (length(dup_dates) > 0) {
    n_disagree <- osp |>
      dplyr::filter(event_date %in% dup_dates) |>
      dplyr::group_by(event_date) |>
      dplyr::summarise(k = dplyr::n_distinct(osp_boat_total), .groups = "drop") |>
      dplyr::filter(k > 1) |> nrow()
    resolve <- params$osp_dupe_resolve %||% "mean"
    agg <- switch(resolve,
                  mean = base::mean, sum = base::sum, max = base::max, first = dplyr::first,
                  stop("params$osp_dupe_resolve must be mean|sum|max|first (got '", resolve, "')", call. = FALSE))
    osp <- osp |>
      dplyr::group_by(event_date) |>
      dplyr::summarise(osp_boat_total = agg(osp_boat_total),
                       osp_crab_only  = if (all(is.na(osp_crab_only))) NA_real_
                                        else agg(osp_crab_only[!is.na(osp_crab_only)]),
                       osp_checked    = if (all(is.na(osp_checked))) NA_real_
                                        else agg(osp_checked[!is.na(osp_checked)]),
                       osp_rate       = if (all(is.na(osp_rate))) NA_real_
                                        else agg(osp_rate[!is.na(osp_rate)]),
                       .groups = "drop")
    cat(sprintf(paste0("  De-dup: %d date(s) had >1 row; %d had DISAGREEING values ",
                       "(combined by '%s'); the rest were identical copies collapsed losslessly.\n"),
                length(dup_dates), n_disagree, resolve))
  } else {
    osp <- osp |> dplyr::select(event_date, osp_boat_total, osp_crab_only, osp_checked, osp_rate)
  }

  # --- Restrict to the estimation window; OSP-dark days stay ABSENT (= NS/latent) ---
  d0 <- as.Date(params$est_date_start %||% "2024-09-16")
  d1 <- as.Date(params$est_date_end   %||% "2025-09-15")
  osp <- osp |> dplyr::filter(event_date >= d0, event_date <= d1) |> dplyr::arrange(event_date)

  out <- osp |>
    dplyr::transmute(
      event_date,
      count_sequence = 1L,                 # OSP is one daily total per day
      count_quantity = osp_boat_total,     # DAILY BOAT TOTAL (all private boats)
      section_num    = 1L,
      count_type     = "OSP Boat Count",
      population     = "private_boat")

  n_zero <- sum(out$count_quantity == 0)
  cat(sprintf("  OSP boat-count days in window: %d (%d observed zeros)%s\n",
              nrow(out), n_zero,
              if (nrow(out)) sprintf("; range %s to %s", min(out$event_date), max(out$event_date)) else ""))
  cat("  NS handling: absent dates are NON-SAMPLED (latent, imputed by the AR); observed zeros are kept as data.\n")

  # --- improvement 8: crab-only rows (the lower bound on f) -------------------
  crab_rows <- tibble::tibble(event_date = as.Date(character()),
                              osp_total = numeric(), osp_crab_only = numeric(), osp_boat_total = numeric(),
                              osp_sample_rate = numeric(), osp_rate_source = character())
  absent <- tibble::tibble(event_date = as.Date(character()), osp_boat_total = numeric(), read_as = character())
  if (has_crab_col) {
    v <- osp$osp_crab_only[is.finite(osp$osp_crab_only)]
    if (identical(crab_unit, "count") && length(v) && all(v <= 1) && any(v != round(v)) &&
        any(osp$osp_boat_total > 1, na.rm = TRUE))
      stop(sprintf(paste0("fetch_osp_boat_counts: '%s' holds values in [0, 1] with fractions, against daily totals above 1: ",
                          "it looks like a SHARE of the boats, not a count. Set osp_crab_only_unit = \"fraction\" in run_config.R ",
                          "(or deliver counts)."), crab_col), call. = FALSE)
    # D39 (2026-10-19): a day with boats and NO crab-only value. OSP's sheet appears to be built from
    # crab-trip records, which cannot list a zero day (17 of the 2024-25 window's 140 OSP days with
    # boats, 9 of them in September and October 2024; CHANGE_REGISTER D39). osp_absent_day says how
    # such a day is read:
    #   "unclassified" (the default, the reading since B64): the day is dropped from the crab-only
    #                  rows, so the share is fitted to zero-truncated data and reads HIGH if the
    #                  absent days are zero days;
    #   "zero"         the day is a classified day with 0 crab-only boats (n from the sampling rate).
    # The question is with OSP; until it answers, "zero" is the sensitivity, not the method.
    absent_mode <- params$osp_absent_day %||% "unclassified"
    if (!absent_mode %in% c("unclassified", "zero"))
      stop("fetch_osp_boat_counts: osp_absent_day must be \"unclassified\" or \"zero\" (got '", absent_mode, "').", call. = FALSE)
    absent <- osp |> dplyr::filter(is.na(osp_crab_only), osp_boat_total > 0) |>
      dplyr::transmute(event_date, osp_boat_total, read_as = absent_mode)
    if (nrow(absent) > 0) {
      .osp_note(sprintf(paste0("%d OSP day(s) with boats (%.0f returning in all) carry no crab-only value and are read as %s ",
                               "(osp_absent_day = \"%s\"; D39: if they are zero days the fitted crab-only share is biased up under ",
                               "\"unclassified\"; the question is with OSP; osp_absent_days.csv lists them)"),
                        nrow(absent), sum(absent$osp_boat_total),
                        if (identical(absent_mode, "zero")) "0 crab-only boats, n from the day's sampling rate" else "UNCLASSIFIED and dropped",
                        absent_mode))
      if (identical(absent_mode, "zero"))
        osp <- osp |> dplyr::mutate(osp_crab_only = dplyr::if_else(is.na(osp_crab_only) & osp_boat_total > 0, 0, osp_crab_only))
    }
    cr <- osp |> dplyr::filter(!is.na(osp_crab_only))
    sn <- osp_resolve_sample_n(cr$osp_boat_total, cr$osp_checked, cr$osp_rate, params)
    cr <- cr |>
      dplyr::mutate(osp_n = sn$n, osp_sample_rate = sn$rate, osp_rate_source = sn$source,
                    osp_crab_only = if (identical(crab_unit, "fraction")) round(osp_crab_only * osp_n)
                                    else if (identical(crab_basis, "expanded")) round(osp_crab_only * osp_sample_rate)
                                    else osp_crab_only) |>
      dplyr::filter(!is.na(osp_crab_only), osp_crab_only >= 0, osp_n > 0) |>
      dplyr::transmute(event_date,
                       osp_total      = osp_n,              # the binomial n: boats SAMPLED (classified)
                       osp_crab_only  = osp_crab_only,
                       osp_boat_total = osp_boat_total,
                       osp_sample_rate = osp_sample_rate,
                       osp_rate_source = osp_rate_source)
    if (nrow(cr)) {
      tb <- table(cr$osp_rate_source)
      cat(sprintf("  OSP sampling: the crab-only share is out of %.0f SAMPLED boats (of %.0f returning); n from %s.\n",
                  sum(cr$osp_total), sum(cr$osp_boat_total), paste(sprintf("%s on %d day(s)", names(tb), as.integer(tb)), collapse = ", ")))
      if (any(cr$osp_rate_source == "schedule (minimum rate)")) {
        cat(paste0("  NOTE: days on the schedule take the manual's MINIMUM rate; OSP may have sampled above it, so on\n",
                   "        those days n is a lower bound and the crab-only share an upper bound. Deliver the day's rate\n",
                   "        (osp_sample_rate_col) or the number sampled (osp_crab_checked_col) to remove this.\n"))
        # D41 (a), 2026-10-19: the same fact for the report and run_warnings.csv
        .osp_note(sprintf(paste0("OSP crab-only share: on %d of %d OSP day(s) n is the manual's MINIMUM sampling rate (no daily rate or ",
                                 "sampled count delivered), so there n is a lower bound and the share an UPPER bound on the truth (D2); ",
                                 "osp_sampling_rate_source = \"none\" (every boat sampled) is the bound in the other direction"),
                          sum(cr$osp_rate_source == "schedule (minimum rate)"), nrow(cr)))
      }
    }
    n_over <- sum(cr$osp_crab_only > cr$osp_total)
    if (n_over > 0) {
      cat(sprintf(paste0("  WARNING: %d OSP day(s) report MORE crab-only boats than boats sampled. ",
                         "That is a data error, not a combo-trip effect; clamped to the total. ",
                         "Check the delivery.\n"), n_over))
      cr$osp_crab_only <- pmin(cr$osp_crab_only, cr$osp_total)
    }
    crab_rows <- cr
    if (nrow(cr) > 0)
      cat(sprintf(paste0("  OSP crab-only classification: %d day(s), %.0f crab-only of %.0f sampled boats ",
                         "(raw share %.3f). This is a LOWER BOUND on the crabbing fraction f: OSP ",
                         "labels combo trips by the non-crab fishery.\n"),
                  nrow(cr), sum(cr$osp_crab_only), sum(cr$osp_total),
                  sum(cr$osp_crab_only) / max(sum(cr$osp_total), 1)))
  } else {
    cat(sprintf("  OSP crab-only column '%s' absent; the f lower bound is inert this run.\n", crab_col))
  }
  attr(out, "osp_crab_rows") <- crab_rows
  attr(out, "osp_absent_days") <- absent
  out
}
