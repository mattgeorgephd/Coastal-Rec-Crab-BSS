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
#   crabbing-fraction helper turns into a HARD LOWER BOUND on f. The column is a count of
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
#   osp_crab_only_col     = "WestportCrabOnlyEffort"   # optional; improvement 8
#   osp_dupe_resolve      = "mean"    # mean | sum | max | first (only on DISAGREEING dupes)
#   est_date_start / est_date_end     # the estimation window (reused, already in run_config)
#
# RETURNS a tibble with the boat_effort schema:
#   event_date, count_sequence, count_quantity, section_num, count_type, population
# plus attr(., "osp_crab_rows") = tibble(event_date, osp_total, osp_crab_only) or empty.
#
# Requires: dplyr, tibble, readxl, here.
###############################################################################

fetch_osp_boat_counts <- function(params) {
  cat("\n  Reading OSP Westport boat-launch counts...\n")

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
  crab_col <- params$osp_crab_only_col %||% "WestportCrabOnlyEffort"
  raw <- readxl::read_excel(osp_file, sheet = params$osp_boat_counts_sheet %||% "Sheet1")

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
  osp <- raw |>
    dplyr::mutate(osp_boat_total = suppressWarnings(as.numeric(.data[[val_col]])),
                  osp_crab_only  = if (has_crab_col) suppressWarnings(as.numeric(.data[[crab_col]]))
                                   else NA_real_) |>
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
                       .groups = "drop")
    cat(sprintf(paste0("  De-dup: %d date(s) had >1 row; %d had DISAGREEING values ",
                       "(combined by '%s'); the rest were identical copies collapsed losslessly.\n"),
                length(dup_dates), n_disagree, resolve))
  } else {
    osp <- osp |> dplyr::select(event_date, osp_boat_total, osp_crab_only)
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
                              osp_total = numeric(), osp_crab_only = numeric())
  if (has_crab_col) {
    cr <- osp |>
      dplyr::filter(!is.na(osp_crab_only), osp_crab_only >= 0, osp_boat_total > 0) |>
      dplyr::transmute(event_date,
                       osp_total     = osp_boat_total,
                       osp_crab_only = osp_crab_only)
    n_over <- sum(cr$osp_crab_only > cr$osp_total)
    if (n_over > 0) {
      cat(sprintf(paste0("  WARNING: %d OSP day(s) report MORE crab-only boats than total boats. ",
                         "That is a data error, not a combo-trip effect; clamped to the total. ",
                         "Check the delivery.\n"), n_over))
      cr$osp_crab_only <- pmin(cr$osp_crab_only, cr$osp_total)
    }
    crab_rows <- cr
    if (nrow(cr) > 0)
      cat(sprintf(paste0("  OSP crab-only classification: %d day(s), %.0f crab-only of %.0f total boats ",
                         "(raw share %.3f). This is a LOWER BOUND on the crabbing fraction f: OSP ",
                         "labels combo trips by the non-crab fishery.\n"),
                  nrow(cr), sum(cr$osp_crab_only), sum(cr$osp_total),
                  sum(cr$osp_crab_only) / max(sum(cr$osp_total), 1)))
  } else {
    cat(sprintf("  OSP crab-only column '%s' absent; the f lower bound is inert this run.\n", crab_col))
  }
  attr(out, "osp_crab_rows") <- crab_rows
  out
}
