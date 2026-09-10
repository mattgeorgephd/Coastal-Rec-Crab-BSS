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
# build_crabbing_holidays.R  (2026-09-10)
#
# Builds 04_input_files/crabbing_holidays.xlsx (sheet "data"): the crabbing-holiday
# calendar the day-typing reads (read_crabbing_holidays.R; a holiday is typed as a
# weekend day). The 2024-25 rows are the list that was hard-coded in the drivers and
# moved to the workbook on 2026-07-16; the other seasons apply THE SAME NAMED HOLIDAYS by
# rule, so a multi-season run types every season the same way:
#
#   Native American Heritage Day   the Friday after Thanksgiving
#   New Year's Eve / New Year's Day
#   Super Bowl Eve                 the Saturday before the Super Bowl (a Saturday, so a
#                                  no-op for day typing; kept because the 2024-25 list
#                                  and the I/E workbook's crabbing_holiday flag carry it)
#   Memorial Day weekend           Saturday, Sunday, Monday
#   Father's Day                   (a Sunday; no-op)
#   Independence Day               July 4, and the federal OBSERVED day when July 4 falls
#                                  on a weekend (Fri Jul 3, 2026: the samplers flagged it
#                                  and Float 20 carried 118 gear, 3.8x the July weekday
#                                  mean; Mon Jul 5, 2027)
#   Labor Day
#
# What the samplers' own Holiday? flag adds (sampler_shifts.xlsx, column holiday) and the
# effort counts say about it is in ../README.md: Veterans Day (observed), Thanksgiving
# Day, Juneteenth read 1.7 to 2.4x the month's weekday gear count at Float 20 on the days
# they were sampled; Christmas Eve, Presidents Day and MLK Day do not. They are NOT in
# this calendar; adding a name here is the way to change the rule.
#
# Super Bowl dates: LVII 2023-02-12, LVIII 2024-02-11, LIX 2025-02-09, LX 2026-02-08,
# LXI 2027-02-14 (scheduled).
#
# COLUMNS: season, date (ISO text), holiday_name, note.
# USAGE (from the repository root): Rscript 04_input_files/build_crabbing_holidays.R
###############################################################################

source(file.path(if (dir.exists("04_input_files")) "04_input_files" else ".", "build_helpers.R"))
root <- build_root()
out_path <- file.path(root, "04_input_files", "crabbing_holidays.xlsx")

nth_weekday <- function(year, month, wday, n) {        # wday: 1 = Monday ... 7 = Sunday
  d <- as.Date(sprintf("%d-%02d-01", year, month)); w <- as.integer(format(d, "%u"))
  d + ((wday - w) %% 7) + 7 * (n - 1)
}
last_weekday <- function(year, month, wday) {
  d <- as.Date(sprintf("%d-%02d-01", year, month + 1)) - 1; if (month == 12) d <- as.Date(sprintf("%d-12-31", year))
  w <- as.integer(format(d, "%u")); d - ((w - wday) %% 7)
}
observed <- function(d) { w <- format(d, "%u"); if (w == "6") d - 1 else if (w == "7") d + 1 else d }
super_bowl <- c("2023" = "2023-02-12", "2024" = "2024-02-11", "2025" = "2025-02-09", "2026" = "2026-02-08", "2027" = "2027-02-14")

season_rows <- function(y0) {                      # season y0-(y0+1): Sep 16 y0 .. Sep 15 y0+1
  y1 <- y0 + 1L; season <- sprintf("%d-%02d", y0, y1 %% 100)
  thanks <- nth_weekday(y0, 11, 4, 4); mem <- last_weekday(y1, 5, 1); jul4 <- as.Date(sprintf("%d-07-04", y1)); jul4o <- observed(jul4)
  rows <- tribble(~date, ~holiday_name, ~note,
    thanks + 1,                              "Native American Heritage Day",    "Friday after Thanksgiving",
    as.Date(sprintf("%d-12-31", y0)),        "New Year's Eve",                  NA_character_,
    as.Date(sprintf("%d-01-01", y1)),        "New Year's Day",                  NA_character_,
    as.Date(super_bowl[[as.character(y1)]]) - 1, "Super Bowl Eve",             "Saturday; no-op for day typing",
    mem - 2,                                 "Memorial Day weekend - Saturday", "no-op for day typing",
    mem - 1,                                 "Memorial Day weekend - Sunday",   "no-op for day typing",
    mem,                                     "Memorial Day",                    NA_character_,
    nth_weekday(y1, 6, 7, 3),                "Father's Day",                    "Sunday; no-op for day typing",
    jul4,                                    "Independence Day",                NA_character_,
    nth_weekday(y1, 9, 1, 1),                "Labor Day",                       NA_character_
  )
  if (jul4o != jul4) rows <- bind_rows(rows, tibble(date = jul4o, holiday_name = "Independence Day (observed)", note = "federal observed day"))
  rows |> mutate(season = season) |> arrange(date)
}

hol <- map_dfr(2022:2026, season_rows) |>
  mutate(date = format(date, "%Y-%m-%d")) |>
  select(season, date, holiday_name, note)

# the 2024-25 rows must reproduce the shipped list exactly
ref <- c("2024-11-29", "2024-12-31", "2025-01-01", "2025-02-08", "2025-05-24", "2025-05-25", "2025-05-26", "2025-06-15", "2025-07-04", "2025-09-01")
got <- hol$date[hol$season == "2024-25"]
if (!setequal(ref, got)) stop("2024-25 holidays differ from the shipped list: ", paste(setdiff(union(ref, got), intersect(ref, got)), collapse = ", "))

cat("Crabbing holidays by season:\n"); print(as.data.frame(hol), row.names = FALSE)
# cross-check against the samplers' flags where a shift workbook exists
sp <- file.path(root, "04_input_files", "sampler_shifts.xlsx")
if (file.exists(sp)) {
  sh <- suppressWarnings(readxl::read_excel(sp, sheet = "data", guess_max = 1e5))
  fl <- sh |> filter(holiday %in% 1) |> distinct(date) |> pull(date)
  cat("\nSampler-flagged dates NOT in the calendar (weekday flags are candidates for the rule):\n")
  miss <- setdiff(fl, hol$date); print(data.frame(date = miss, dow = weekdays(as.Date(miss))), row.names = FALSE)
  cat("\nCalendar dates the samplers did not flag on a sampled day (sampled but flag 0):\n")
  sampled0 <- sh |> filter(date %in% hol$date) |> group_by(date) |> summarise(flag = max(holiday, na.rm = TRUE), .groups = "drop") |> filter(flag == 0)
  print(as.data.frame(sampled0), row.names = FALSE)
}
write_data_sheet(hol, out_path)
