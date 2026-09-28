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
# build_effort_combined.R  (2026-09-10)
#
# Builds 04_input_files/effort_combined.xlsx (sheet "data") from the "Effort count" sheet
# of every per-season creel workbook in raw/. One row per instantaneous count at a site,
# every season. The 2026-07-16 workbook held 2024-25 only (3,256 rows); this one holds
# 2022-23 through 2025-26 (12,056 rows through 2026-09-08), which is what the staged
# two-season run (2023-24 + 2024-25) and a 2025-26 run were waiting for
# (CHANGE_REGISTER D8: "effort_combined.xlsx has NO 2023-24 rows").
#
# THE COUNT COLUMNS ACROSS SEASONS. The 2022-23 and 2023-24 sheets record
# "Truck/Boat Trailer Count" (the truck-and-trailer rigs in the launch lot) and a "Buoy
# Count"; from 2024-25 the protocol records "Boat Trailer Count" and, separately, a
# "Vehicle Count" (at Westport Boat Launch the vehicle count is zero on every 2024-25 and
# 2025-26 row, so the trailer count is the boat-effort count in every season). The
# builder maps both trailer headers to boat_trailer_count (the column the reader uses)
# and carries every other count under its own name; a count a season did not record is
# NA in that season.
#
# QC FLAGS (rows kept, held out of the model by the reader):
#   gear_count_from_interviews  During the L&I PFD stand-down of Feb 3 to 13, 2024 the
#       samplers could not walk the floats, and the Float 20 / Float 17-21 gear "counts"
#       on those days are the day's TOTAL gear from interviews (the notes say so: "Total
#       from Interviews", "gear count from interviews", "Based on interviews"), not an
#       instantaneous count. An instantaneous count times the turnover is the day's
#       gear-deployments; a day total already is, so feeding it to the count model would
#       inflate the shore effort by about the turnover (2.5x). The same note appears on
#       three later 2023-24 dock rows and on Tokeland rows. 20 rows in all, every one 2023-24.
#   missing_count_time          No clock time (the reader already drops these; flagged so
#       the loss is visible). 27 + 3 rows in 2022-23 / 2023-24, none later.
#   implausible_count_time      A clock time before 05:00 (three 2024-25 rows: "02:25",
#       "03:10" at Float 20, "04:52" at Mocrocks; 12-hour-clock slips, the earliest
#       check-in in four seasons is 04:50). INFORMATIONAL: the count value is real and
#       stays in the model; only its clock time is wrong, and the turnover derivation
#       already drops hours outside the I/E window. The reader drops only the flags in
#       run_config$effort_qc_drop (default: gear_count_from_interviews).
#
# COLUMNS: season, date (ISO text), survey_id ("S<n>"), creel_area (trimmed, the 2022-23
# names mapped as in build_helpers.R), count_time ("HH:MM:SS"), total_gear_count,
# boat_trailer_count, vehicle_count, boats_entering_marina, buoy_count, crabber_count,
# jetty_people_count, qc_flag. Notes stay in raw/.
#
# USAGE (from the repository root): Rscript 04_input_files/build_effort_combined.R
###############################################################################

source(file.path(if (dir.exists("04_input_files")) "04_input_files" else ".", "build_helpers.R"))
root <- build_root()
out_path <- file.path(root, "04_input_files", "effort_combined.xlsx")

wb <- season_workbooks(root)
say("Season workbooks: %s", paste(sprintf("%s (%s)", basename(wb$path), wb$season), collapse = ", "))

read_season <- function(path, season) {
  es <- sheet_named(path, c("Effort count"))
  if (is.null(es)) stop("no 'Effort count' sheet in ", basename(path))
  e <- read_sheet_text(path, es)
  date_raw <- pick_col(e, "Date"); date_d <- parse_date_any(date_raw); report_unparsed(sprintf("%s date", season), date_raw, date_d)
  ct_raw <- pick_col(e, "Count Time"); ct <- parse_clock_hours(ct_raw); report_unparsed(sprintf("%s count time", season), ct_raw, ct)
  notes <- str_squish(pick_col(e, "Notes"))
  tibble(
    season      = season,
    date_d      = date_d,
    date        = format(date_d, "%Y-%m-%d"),
    survey_id   = survey_id_of(pick_col(e, "ID")),
    creel_area  = harmonise_area(pick_col(e, c("Creel Area", "creel Area"))),
    count_time  = hhmmss(ct),
    total_gear_count      = num(pick_col(e, "Total Gear Count")),
    boat_trailer_count    = num(pick_col(e, c("Boat Trailer Count", "Truck/Boat Trailer Count"))),
    vehicle_count         = num(pick_col(e, "Vehicle Count")),
    boats_entering_marina = num(pick_col(e, "Boats Entering Marina Count")),
    buoy_count            = num(pick_col(e, "Buoy Count")),
    crabber_count         = num(pick_col(e, "Crabber Count")),
    jetty_people_count    = num(pick_col(e, "Jetty People Count")),
    .notes = notes
  ) |>
    mutate(
      qc_flag = case_when(
        !is.na(.notes) & str_detect(.notes, "(?i)from interviews|based on interviews|interviews/observ|total from interview") &
          !is.na(total_gear_count) & total_gear_count > 0 ~ "gear_count_from_interviews",
        is.na(count_time) ~ "missing_count_time",
        !is.na(ct) & ct < 5 ~ "implausible_count_time",
        TRUE ~ ""
      )
    )
}

combined <- pmap_dfr(list(wb$path, wb$season), read_season)

bad_season <- combined |> filter(!is.na(date_d), season_of(date_d) != season)
if (nrow(bad_season)) stop(sprintf("%d effort row(s) dated outside their workbook's fishery season, e.g. %s in %s",
                                   nrow(bad_season), bad_season$date[1], bad_season$season[1]))
if (any(is.na(combined$date))) say("  NOTE %d row(s) without a parseable date", sum(is.na(combined$date)))

flagged <- combined |> filter(qc_flag == "gear_count_from_interviews")
cat("\nRows flagged gear_count_from_interviews (held out of the shore counts by the reader):\n")
print(as.data.frame(flagged |> select(season, date, creel_area, count_time, total_gear_count, .notes)), right = FALSE)

combined <- combined |> arrange(season, date, survey_id, creel_area, count_time) |>
  select(season, date, survey_id, creel_area, count_time, total_gear_count, boat_trailer_count, vehicle_count,
         boats_entering_marina, buoy_count, crabber_count, jetty_people_count, qc_flag)

cat("\nRows by season:\n"); print(table(combined$season))
cat("\nDate range by season:\n"); print(as.data.frame(combined |> group_by(season) |> summarise(first = min(date, na.rm = TRUE), last = max(date, na.rm = TRUE), rows = n(), days = n_distinct(date), .groups = "drop")))
cat("\nqc_flag by season:\n"); print(table(combined$season, combined$qc_flag, useNA = "ifany"))
cat("\nWestport sites, counts per season (rows / days / mean trailer / mean gear):\n")
print(as.data.frame(combined |> filter(str_detect(creel_area, "^Westport|Ocean Shores")) |>
  group_by(season, creel_area) |> summarise(rows = n(), days = n_distinct(date), mean_trailers = round(mean(boat_trailer_count, na.rm = TRUE), 2),
                                            mean_gear = round(mean(total_gear_count, na.rm = TRUE), 2), .groups = "drop")))

if (file.exists(out_path)) {
  old <- suppressWarnings(readxl::read_excel(out_path, sheet = "data", guess_max = 1e5))
  cat("\nComparison with the previous workbook on (season, survey_id, creel_area, count_time):\n")
  keys <- c("season", "survey_id", "creel_area", "count_time")
  o <- old |> mutate(creel_area = harmonise_area(creel_area)) |> distinct(across(all_of(keys)), .keep_all = TRUE)
  n <- combined |> distinct(across(all_of(keys)), .keep_all = TRUE)
  compare_with_previous(o, n, keys, c("date", "total_gear_count", "boat_trailer_count"))
}

write_data_sheet(combined, out_path)
