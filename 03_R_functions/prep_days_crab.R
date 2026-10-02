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
# prep_days_crab.R
#
# Build the per-day calendar (day index, week/month index, day type, day-type
# integer, and I/E-derived effective day length) for a date range. Extracted from
# the pooled and gear-resolved drivers so both share one implementation. The gear
# superset (which adds the Stan day_type_idx column) is used; the extra column is
# harmless for the pooled model, which does not read it. Auto-sourced by both
# drivers via the 03_R_functions walk.
#
# Signature takes `params` and derives weekends / holidays / period_pe / sections
# from it, so the function is pure at source time. Day-length assignment is
# delegated to bss_assign_day_length() in bss_day_length.R.
###############################################################################

prep_days_crab <- function(date_begin, date_end, params, L_eff_model = NULL) {
  # Derive day-typing inputs from the centralized config (single source of truth).
  weekends      <- params$days_wkend
  holiday_dates <- params$crabbing_holiday_dates
  period_pe     <- params$period_pe
  sections      <- params$sections
  date_begin <- as.Date(date_begin); date_end <- as.Date(date_end)
  days <- tibble(
    event_date = seq.Date(date_begin, date_end, by="day"),
    day = bss_weekday(event_date),
    day_type = case_when(
      event_date %in% holiday_dates ~ "holiday",
      day %in% weekends ~ "weekend", TRUE ~ "weekday"),
    # v5.1: Integer day type index for Stan (1=weekday, 2=weekend, 3=holiday)
    day_type_idx = case_when(
      day_type == "weekday" ~ 1L,
      day_type == "weekend" ~ 2L,
      day_type == "holiday" ~ 3L),
    day_type_num_weekend = as.integer(day_type %in% c("weekend","holiday")),
    day_type_num_holiday = as.integer(day_type == "holiday"),
    # 2026-09-28 (B46): weeks are ISO 8601 weeks keyed WITH their ISO year (%G-%V, Monday
    # start, as %W was). %W restarts at 0 every 1 January, so (a) the PE's (week x day_type)
    # strata pooled week N of one year with week N of the next whenever a window spans both
    # (a window opening before mid-September does, within a single season), and (b) the
    # New Year week was cut into two short periods (Mon 30 Dec to Sun 5 Jan 2025 became a
    # 2-day and a 5-day week, in the PE strata and the weekly AR alike). Month periods carry
    # their year for the same reason.
    week = as.numeric(format(event_date,"%V")),
    iso_year = as.numeric(format(event_date,"%G")),
    month = as.numeric(format(event_date,"%m")),
    year = as.numeric(format(event_date,"%Y")),
    # 2026-10-02: an if/else, not case_when(). period_pe is ONE value for the whole window, and
    # dplyr 1.2.0 (the renv.lock version) soft-deprecates a case_when() whose conditions are
    # all scalars while its outcomes are columns. The result is identical (harness section 95).
    period = if (isTRUE(period_pe == "month")) year * 100 + month else iso_year * 100 + week,
    day_index = as.integer(seq_along(event_date)),
    week_index = as.integer(factor(
      paste(iso_year, sprintf("%02d", week)),
      levels = unique(paste(iso_year, sprintf("%02d", week)))
    )),
    month_index = as.integer(factor(paste(year,sprintf("%02d",month)),
                  levels=unique(paste(year,sprintf("%02d",month))))),
    day_length = NA_real_
  )
  # Civil twilight + L_effective assignment. Shared with the pooled driver:
  # 03_R_functions/bss_day_length.R. Sets day_length, day_length_civil_twilight,
  # L_mu and L_prior_sigma. Falls back to civil twilight (clamped to
  # [day_length_min_hours, day_length_max_hours], default [9, 17]) only when
  # L_eff_model is NULL, i.e. no usable I/E data.
  days <- bss_assign_day_length(days, L_eff_model, params)

  for(s in sections) days[[paste0("open_section_",s)]] <- TRUE
  days
}
