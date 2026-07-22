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
    day = weekdays(event_date),
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
    week = as.numeric(format(event_date,"%W")),
    month = as.numeric(format(event_date,"%m")),
    year = as.numeric(format(event_date,"%Y")),
    period = case_when(
      period_pe == "month" ~ as.numeric(format(event_date, "%m")),
      period_pe == "week"  ~ as.numeric(format(event_date, "%W")),
      TRUE ~ as.numeric(format(event_date, "%W"))
    ),
    day_index = as.integer(seq_along(event_date)),
    week_index = as.integer(factor(
      paste(year, sprintf("%02d", week)),
      levels = unique(paste(year, sprintf("%02d", week)))
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
