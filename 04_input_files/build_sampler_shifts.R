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
# build_sampler_shifts.R  (2026-09-09; rebuilt on the season workbooks 2026-09-10)
#
# Builds 04_input_files/sampler_shifts.xlsx (sheet "data") from the "crab creel survey
# data" sheet of every per-season creel workbook in raw/: one row per sampler site-visit
# (a survey), with the check-in and check-out in 24-hour clock time.
#
# WHY. The crabbing-fraction classification comes from the boats the samplers contacted,
# and samplers are in port for a shift (about 9:45 to 15:50 at Westport in 2024-25), so
# a boat returning outside the shift is never classified. The shift windows are what a
# return-time weighting of the classification needs once OSP's all-day counts arrive,
# and until then they support the shift-coverage diagnostic
# (03_R_functions/sampler_shifts.R). The survey id is also the link between an interview
# or an effort count (survey_id = "S<n>") and the shift it was taken in.
#
# WHAT THE SEASON WORKBOOKS ADD (2026-09-10) over the 2026-09-09 survey export: the
# sampler's on-site conditions, recorded once per survey (tide stage in 2022-23 and
# 2023-24; rain, cloud, wind and wind direction every season; the weather station site
# from 2024-25), which are the covariate candidates for the effort process nearest to
# what a crabber saw when deciding to go (../README.md); and, in 2025-26, the Westport
# vessel tallies (commercial / private / charter), which build_comm_charter_tally.R reads.
#
# COLUMNS WRITTEN
#   survey_id        "S<n>", matching interview_combined.xlsx / effort_combined.xlsx
#   survey_num       the numeric ID
#   season           fishery season (the workbook's; checked against the date)
#   date             ISO text yyyy-mm-dd
#   day_of_week, holiday (0/1: the sampler's flag), creel_location, weather_location,
#   samplers, tide, rain, weather, wind, wind_direction, special_conditions, notes
#   check_in, check_out          "HH:MM" (24-hour)
#   check_in_hour, check_out_hour decimal hours (10.72 = 10:43)
#   shift_hours                  check_out_hour - check_in_hour; NA when flagged
#   qc_flag          "" | "missing_check_in" | "missing_check_out" |
#                    "check_out_before_check_in" | "duplicate_survey_id"; a flagged row
#                    keeps its raw values and gets NA hours, so nothing is silently
#                    corrected
#
# USAGE (from the repository root): Rscript 04_input_files/build_sampler_shifts.R
###############################################################################

source(file.path(if (dir.exists("04_input_files")) "04_input_files" else ".", "build_helpers.R"))
root <- build_root()
out_path <- file.path(root, "04_input_files", "sampler_shifts.xlsx")

wb <- season_workbooks(root)
say("Season workbooks: %s", paste(sprintf("%s (%s)", basename(wb$path), wb$season), collapse = ", "))

txt <- function(x) { x <- str_squish(as.character(x)); ifelse(is.na(x) | x == "", NA_character_, x) }

read_season <- function(path, season) {
  ss <- sheet_named(path, c("crab creel survey data"))
  if (is.null(ss)) stop("no 'crab creel survey data' sheet in ", basename(path))
  s <- read_sheet_text(path, ss)
  date_raw <- pick_col(s, "Date"); date_d <- parse_date_any(date_raw); report_unparsed(sprintf("%s date", season), date_raw, date_d)
  ci_raw <- pick_col(s, "Check In Time"); ci <- parse_clock_hours(ci_raw); report_unparsed(sprintf("%s check-in", season), ci_raw, ci)
  co_raw <- pick_col(s, "Check Out Time"); co <- parse_clock_hours(co_raw); report_unparsed(sprintf("%s check-out", season), co_raw, co)
  tibble(
    survey_num  = suppressWarnings(as.integer(round(num(pick_col(s, "ID"))))),
    season      = season,
    date_d      = date_d,
    date        = format(date_d, "%Y-%m-%d"),
    day_of_week = txt(pick_col(s, "Day of Week")),
    holiday     = { h <- num(pick_col(s, c("Holiday?", "Holiday"))); ifelse(is.na(h), NA_integer_, as.integer(h != 0)) },
    creel_location   = txt(pick_col(s, "Creel Location")),
    weather_location = txt(pick_col(s, "Weather Location")),
    samplers         = txt(pick_col(s, c("Sampler(s)", "Samplers"))),
    tide             = txt(pick_col(s, "Tide")),
    rain             = txt(pick_col(s, "Rain")),
    weather          = txt(pick_col(s, "Weather")),
    wind             = txt(pick_col(s, "Wind")),
    wind_direction   = txt(pick_col(s, "Wind Direction")),
    special_conditions = txt(pick_col(s, "Special Conditions")),
    check_in_hour  = ci,
    check_out_hour = co,
    notes = txt(pick_col(s, "Notes"))
  ) |>
    mutate(survey_id = ifelse(is.na(survey_num), NA_character_, paste0("S", survey_num)))
}

shifts <- pmap_dfr(list(wb$path, wb$season), read_season)

bad_season <- shifts |> filter(!is.na(date_d), season_of(date_d) != season)
if (nrow(bad_season)) stop(sprintf("%d survey row(s) dated outside their workbook's fishery season, e.g. %s in %s",
                                   nrow(bad_season), bad_season$date[1], bad_season$season[1]))

shifts <- shifts |>
  mutate(
    qc_flag = case_when(
      is.na(check_in_hour)                              ~ "missing_check_in",
      is.na(check_out_hour)                             ~ "missing_check_out",
      check_out_hour < check_in_hour                    ~ "check_out_before_check_in",
      TRUE ~ ""),
    .dup = !is.na(survey_id) & (duplicated(survey_id) | duplicated(survey_id, fromLast = TRUE)),
    qc_flag = ifelse(.dup, paste0(qc_flag, ifelse(nzchar(qc_flag), ";", ""), "duplicate_survey_id"), qc_flag),
    shift_hours = ifelse(qc_flag == "", check_out_hour - check_in_hour, NA_real_),
    check_in    = hhmm(check_in_hour),
    check_out   = hhmm(check_out_hour)
  ) |>
  select(survey_id, survey_num, season, date, day_of_week, holiday, creel_location, weather_location, samplers,
         tide, rain, weather, wind, wind_direction, special_conditions,
         check_in, check_out, check_in_hour, check_out_hour, shift_hours, qc_flag, notes) |>
  arrange(date, survey_num)

if (any(is.na(shifts$date))) warning(sprintf("%d row(s) with an unparseable date", sum(is.na(shifts$date))))
cat("\nSurveys by season x creel_location:\n"); print(table(shifts$season, shifts$creel_location))
cat("\nDate range by season:\n"); print(as.data.frame(shifts |> group_by(season) |> summarise(first = min(date), last = max(date), rows = n(), .groups = "drop")))
cat("\nQC flags:\n"); print(table(shifts$qc_flag, useNA = "ifany"))
fl <- shifts |> filter(qc_flag != "")
if (nrow(fl)) { cat("flagged rows:\n"); print(as.data.frame(fl |> select(survey_id, date, creel_location, check_in, check_out, qc_flag))) }
gh <- shifts |> filter(creel_location == "Grays Harbor", !is.na(shift_hours))
cat("\nGrays Harbor shifts by season (usable rows; medians):\n")
print(as.data.frame(gh |> group_by(season) |> summarise(n = n(), check_in = hhmm(median(check_in_hour)), check_out = hhmm(median(check_out_hour)), shift_h = round(median(shift_hours), 2), .groups = "drop")))
cat("\nSampler-flagged holidays (dates any survey flagged holiday = 1):\n")
print(as.data.frame(shifts |> filter(holiday %in% 1) |> distinct(season, date, day_of_week) |> arrange(date)))

if (file.exists(out_path)) {
  old <- suppressWarnings(readxl::read_excel(out_path, sheet = "data", guess_max = 1e5))
  cat("\nComparison with the previous workbook on survey_id:\n")
  compare_with_previous(old |> distinct(survey_id, .keep_all = TRUE), shifts |> distinct(survey_id, .keep_all = TRUE), "survey_id",
                        c("season", "date", "creel_location", "check_in_hour", "check_out_hour", "shift_hours", "qc_flag", "holiday"))
}

write_data_sheet(shifts, out_path)
