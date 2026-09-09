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
# sampler_shifts.R  (2026-09-09)
#
# The sampler SHIFTS (04_input_files/sampler_shifts.xlsx, built from the survey export
# by 04_input_files/build_sampler_shifts.R) and what they say about the crabbing-fraction
# classification.
#
# THE QUESTION. The crabbing fraction f is read from the boats the samplers contacted at
# the launch, and a sampler is in port for a shift (about 9:45 to 15:50 at Westport in
# 2024-25). A boat returning outside the shift is never classified. If the trip-type mix
# of the boats returning after the shift differs from the mix during it (finfish boats
# returning later, say), the shift-time share misstates f. Nothing in the contacts
# themselves can settle that, because they only exist inside the shift; what CAN be
# measured now is (a) how much of the day's boat traffic the shift covers, from the boat
# I/E return profile, and (b) whether the trip-type mix drifts with the hour INSIDE the
# shift, from the contact clock times. OSP's all-day crabbing-only count, when it
# arrives, is the outside-the-shift check, and the return-time weighting hook below is
# where it would enter.
#
# fetch_sampler_shifts(params)          the shifts for this window and location
# diagnose_shift_coverage(...)           per-day coverage, contact hour by trip type,
#                                        written as shift_coverage_*.csv
#
# THE WEIGHTING HOOK (crab_fraction_shift_weighting, ships "none"). A return-time
# weighting re-weights each contact by the ratio (share of the day's returns in its hour,
# from an all-day count) / (share of the CONTACTS in that hour), so that the classified
# sample stands for the whole day's returns. It needs an all-day return profile that the
# I/E boat survey gives on 48 days and OSP will give on every day; and it cannot say
# anything about hours no sampler covered. It is documented here and not implemented,
# because with the 2024-25 within-shift mix flat across the hour the weighting would move
# nothing, and the hours outside the shift are exactly where it has no data.
###############################################################################

fetch_sampler_shifts <- function(params, quiet = FALSE) {
  f <- here("04_input_files", params$sampler_shifts_file %||% "sampler_shifts.xlsx")
  if (!file.exists(f)) {
    if (!isTRUE(quiet)) cat("  Sampler shifts: file not found (", basename(f), "); the shift diagnostic is skipped.\n")
    return(NULL)
  }
  sh <- readxl::read_excel(f, sheet = params$sampler_shifts_sheet %||% "data")
  need <- c("survey_id", "date", "creel_location", "check_in_hour", "check_out_hour", "shift_hours", "qc_flag")
  if (!all(need %in% names(sh))) stop("sampler_shifts.xlsx lacks column(s): ", paste(setdiff(need, names(sh)), collapse = ", "), call. = FALSE)
  sh <- sh |>
    mutate(event_date = as.Date(as.character(date)),
           qc_flag = ifelse(is.na(qc_flag), "", as.character(qc_flag))) |>
    filter(creel_location == (params$gh_creel_location %||% "Grays Harbor"))
  if ("season" %in% names(sh) && !is.null(params$season_filter)) sh <- sh |> filter(season %in% params$season_filter)
  if (!is.null(params$est_date_start) && !is.null(params$est_date_end))
    sh <- sh |> filter(event_date >= as.Date(params$est_date_start), event_date <= as.Date(params$est_date_end))
  n_flag <- sum(nzchar(sh$qc_flag))
  if (!isTRUE(quiet))
    cat(sprintf("  Sampler shifts: %d surveys on %d days in the window (%d flagged and held out of the hours); check-in median %s, check-out median %s\n",
                nrow(sh), n_distinct(sh$event_date), n_flag,
                .hhmm_of(median(sh$check_in_hour[!nzchar(sh$qc_flag)], na.rm = TRUE)),
                .hhmm_of(median(sh$check_out_hour[!nzchar(sh$qc_flag)], na.rm = TRUE))))
  attr(sh, "n_flagged") <- n_flag
  sh
}

.hhmm_of <- function(h) {
  h <- as.numeric(h)
  ifelse(is.finite(h), sprintf("%02d:%02d", floor(h), round((h - floor(h)) * 60)), "NA")
}

# The union of a day's shift windows as a set of [start, end] intervals (two surveys on one
# day are two sites, usually overlapping in time; the union is the hours SOMEONE was in port).
.shift_union <- function(t_in, t_out) {
  ok <- is.finite(t_in) & is.finite(t_out) & t_out >= t_in
  if (!any(ok)) return(matrix(numeric(0), 0, 2))
  iv <- cbind(t_in[ok], t_out[ok]); iv <- iv[order(iv[, 1]), , drop = FALSE]
  out <- iv[1, , drop = FALSE]
  if (nrow(iv) > 1) for (i in 2:nrow(iv)) {
    if (iv[i, 1] <= out[nrow(out), 2]) out[nrow(out), 2] <- max(out[nrow(out), 2], iv[i, 2])
    else out <- rbind(out, iv[i, ])
  }
  out
}
.hours_covered <- function(iv) if (!nrow(iv)) 0 else sum(iv[, 2] - iv[, 1])
# share of a return profile (hour, weight) falling inside the union of windows
.share_inside <- function(hour, w, iv) {
  if (!nrow(iv) || !length(hour) || sum(w) <= 0) return(NA_real_)
  inside <- vapply(hour, function(h) any(h >= iv[, 1] & h <= iv[, 2]), logical(1))
  sum(w[inside]) / sum(w)
}

# shifts   : fetch_sampler_shifts() output (this window)
# ie_data  : fetch_ie_data() output; attr(., "ie_boat_intervals") carries the boat I/E rows
# detail   : dwg$boat_contacts_detail (one row per contacted private boat)
diagnose_shift_coverage <- function(shifts, ie_data, detail, params, output_dir = NULL, quiet = FALSE) {
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  if (is.null(shifts) || !nrow(shifts)) { .say("  Shift coverage: no shifts in the window.\n"); return(invisible(NULL)) }
  ok_sh <- shifts |> filter(!nzchar(qc_flag), is.finite(check_in_hour), is.finite(check_out_hour))

  # --- the boat return profile: WBL when it has enough days, else every boat-I/E site ----
  bi <- attr(ie_data, "ie_boat_intervals")
  min_days <- as.integer(params$shift_coverage_ie_min_days %||% 10L)
  profile <- NULL; profile_src <- "none"
  if (!is.null(bi) && nrow(bi)) {
    wbl <- bi |> filter(location_name == (params$ie_boat_location %||% "WBL"))
    use <- if (n_distinct(wbl$event_date) >= min_days) { profile_src <- params$ie_boat_location %||% "WBL"; wbl }
           else { profile_src <- sprintf("all boat-I/E sites pooled (%s has %d days, floor %d)", params$ie_boat_location %||% "WBL", n_distinct(wbl$event_date), min_days); bi }
    profile <- use |> mutate(hr = floor(hour * 4) / 4) |> group_by(hr) |>
      summarise(returns = sum(boats_out), arrivals = sum(boats_in), .groups = "drop") |>
      mutate(return_share = returns / max(sum(returns), 1), arrival_share = arrivals / max(sum(arrivals), 1))
  }

  # --- per day: the union of the shift windows and the share of returns it covers -------
  by_day <- ok_sh |> group_by(event_date) |> group_modify(function(g, k) {
    iv <- .shift_union(g$check_in_hour, g$check_out_hour)
    tibble(n_surveys = nrow(g), first_check_in = min(g$check_in_hour), last_check_out = max(g$check_out_hour),
           hours_in_port = .hours_covered(iv),
           return_share_covered = if (is.null(profile)) NA_real_ else .share_inside(profile$hr, profile$returns, iv),
           arrival_share_covered = if (is.null(profile)) NA_real_ else .share_inside(profile$hr, profile$arrivals, iv))
  }) |> ungroup()

  # --- contacts: how many per day fell inside the union, and the hour by trip type -------
  cd <- NULL; ch <- NULL; crab_by_hour <- NULL
  if (!is.null(detail) && nrow(detail)) {
    areas <- params$crab_fraction_contact_areas %||% params$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch")
    d <- detail
    if ("creel_area" %in% names(d) && !identical(tolower(areas[1]), "all")) d <- d |> filter(is.na(creel_area) | creel_area %in% areas)
    d <- d |> filter(is.finite(contact_hour))
    if (nrow(d)) {
      d <- d |> left_join(by_day |> select(event_date, first_check_in, last_check_out), by = "event_date") |>
        mutate(inside = is.finite(first_check_in) & contact_hour >= first_check_in - 0.05 & contact_hour <= last_check_out + 0.05,
               cls = ifelse(is.na(trip_type_class), ifelse(is.finite(crabbers) & crabbers > 0, "crabbing (untyped)", "not crabbing (untyped)"), trip_type_class),
               hour_bin = cut(contact_hour, breaks = c(0, 8, 10, 12, 14, 16, 18, 24), right = FALSE,
                              labels = c("<08", "08-10", "10-12", "12-14", "14-16", "16-18", ">=18")))
      cd <- d |> group_by(cls) |>
        summarise(contacts = n(), inside_shift = sum(inside, na.rm = TRUE),
                  hour_q10 = quantile(contact_hour, 0.10), hour_q25 = quantile(contact_hour, 0.25), hour_median = median(contact_hour),
                  hour_q75 = quantile(contact_hour, 0.75), hour_q90 = quantile(contact_hour, 0.90), .groups = "drop")
      ch <- d |> count(hour_bin, cls) |> tidyr::pivot_wider(names_from = cls, values_from = n, values_fill = 0)
      # the within-shift drift of the crabbing share with the hour: a share by hour bin
      crab_by_hour <- d |> filter(!is.na(trip_type_class)) |>
        group_by(hour_bin) |> summarise(contacts = n(), crabbing = sum(trip_type_class %in% c("crab_only", "combo")),
                                        combo = sum(trip_type_class == "combo"), .groups = "drop") |>
        mutate(crabbing_share = crabbing / pmax(contacts, 1), combo_share_of_crabbing = ifelse(crabbing > 0, combo / crabbing, NA_real_))
    }
  }

  summ <- tibble(
    n_shift_days = nrow(by_day),
    median_first_check_in = median(by_day$first_check_in), median_last_check_out = median(by_day$last_check_out),
    median_hours_in_port = median(by_day$hours_in_port),
    return_profile_source = profile_src,
    median_return_share_covered = if (all(is.na(by_day$return_share_covered))) NA_real_ else median(by_day$return_share_covered, na.rm = TRUE),
    returns_before_first_check_in = if (is.null(profile)) NA_real_ else sum(profile$return_share[profile$hr < median(by_day$first_check_in)]),
    returns_after_last_check_out  = if (is.null(profile)) NA_real_ else sum(profile$return_share[profile$hr > median(by_day$last_check_out)]),
    contacts_inside_shift = if (is.null(cd)) NA_real_ else sum(cd$inside_shift) / max(sum(cd$contacts), 1),
    weighting = params$crab_fraction_shift_weighting %||% "none")

  .say(sprintf(paste0("  Shift coverage: %d shift days; in port %s to %s (medians), %.1f h; the boat return profile (%s)",
                      " puts %.0f%% of a day's returns inside the shift (%.0f%% before the first check-in, %.0f%% after the last check-out)\n"),
               summ$n_shift_days, .hhmm_of(summ$median_first_check_in), .hhmm_of(summ$median_last_check_out), summ$median_hours_in_port,
               profile_src, 100 * (summ$median_return_share_covered %||% NA), 100 * (summ$returns_before_first_check_in %||% NA),
               100 * (summ$returns_after_last_check_out %||% NA)))
  if (!is.null(cd)) {
    .say(sprintf("    contact hour by trip type (median [q25, q75]): %s\n",
                 paste(sprintf("%s %s [%s, %s] n=%d", cd$cls, .hhmm_of(cd$hour_median), .hhmm_of(cd$hour_q25), .hhmm_of(cd$hour_q75), cd$contacts), collapse = "; ")))
    if (!is.null(crab_by_hour) && nrow(crab_by_hour))
      .say(sprintf("    crabbing share of typed contacts by hour bin: %s\n",
                   paste(sprintf("%s %.2f (n=%d)", crab_by_hour$hour_bin, crab_by_hour$crabbing_share, crab_by_hour$contacts), collapse = ", ")))
  }
  if (!is.null(output_dir)) {
    utils::write.csv(by_day, file.path(output_dir, "shift_coverage_daily.csv"), row.names = FALSE)
    utils::write.csv(summ, file.path(output_dir, "shift_coverage_summary.csv"), row.names = FALSE)
    if (!is.null(profile)) utils::write.csv(profile, file.path(output_dir, "shift_coverage_return_profile.csv"), row.names = FALSE)
    if (!is.null(cd)) {
      utils::write.csv(cd, file.path(output_dir, "contact_hour_by_trip_type.csv"), row.names = FALSE)
      utils::write.csv(ch, file.path(output_dir, "contact_hour_bins_by_trip_type.csv"), row.names = FALSE)
      if (!is.null(crab_by_hour)) utils::write.csv(crab_by_hour, file.path(output_dir, "crabbing_share_by_contact_hour.csv"), row.names = FALSE)
    }
  }
  invisible(list(by_day = by_day, summary = summ, profile = profile, contact_hours = cd,
                 crabbing_by_hour = crab_by_hour))
}
