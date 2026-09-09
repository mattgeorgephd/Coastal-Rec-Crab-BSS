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
# build_sampler_shifts.R  (2026-09-09)
#
# Formats the creel SURVEY export (one row per sampler site-visit: check-in and check-out
# in 24-hour clock time) into 04_input_files/sampler_shifts.xlsx, sheet "data".
#
# WHY. The crabbing-fraction classification comes from the boats the samplers contacted,
# and samplers are in port for a shift (about 9:45 to 15:50 at Westport in 2024-25), so
# a boat returning outside the shift is never classified. The shift windows are what a
# return-time weighting of the classification needs once OSP's all-day counts arrive,
# and until then they support the shift-coverage diagnostic
# (03_R_functions/sampler_shifts.R). The survey id is also the link between an interview
# (interview_combined.xlsx: survey_id = "S<n>") and the shift it was taken in.
#
# COLUMNS WRITTEN
#   survey_id        "S<n>", matching interview_combined.xlsx / effort_combined.xlsx
#   survey_num       the export's numeric ID
#   season           fishery season by the Sep 16 boundary (a date on or after Sep 16 of
#                    year Y is "Y-(Y+1)"), the convention of every other workbook
#   date             ISO text yyyy-mm-dd
#   day_of_week, holiday (0/1), creel_location, samplers, special_conditions, notes
#   check_in, check_out          "HH:MM" (24-hour)
#   check_in_hour, check_out_hour decimal hours (10.72 = 10:43)
#   shift_hours                  check_out_hour - check_in_hour; NA when flagged
#   qc_flag          "" | "missing_check_out" | "check_out_before_check_in" |
#                    "duplicate_survey_id"; a flagged row keeps its raw strings and gets
#                    NA hours, so nothing is silently corrected
#
# USAGE (from the repository root):
#   Rscript 04_input_files/build_sampler_shifts.R [raw_survey.xlsx] [output.xlsx]
# Defaults: raw/surveydata20222026.xlsx -> sampler_shifts.xlsx.
###############################################################################

suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tibble); library(stringr); library(writexl) })

args <- commandArgs(trailingOnly = TRUE)
.root <- getwd()
if (!dir.exists(file.path(.root, "04_input_files")) && dir.exists(file.path(.root, "..", "04_input_files")))
  .root <- normalizePath(file.path(.root, ".."))
raw_path <- if (length(args) >= 1) args[1] else file.path(.root, "04_input_files", "raw", "surveydata20222026.xlsx")
out_path <- if (length(args) >= 2) args[2] else file.path(.root, "04_input_files", "sampler_shifts.xlsx")
stopifnot(file.exists(raw_path))

raw <- suppressWarnings(readxl::read_excel(raw_path, sheet = 1, .name_repair = "minimal"))
need <- c("ID", "Date", "Day of Week", "Holiday", "Samplers", "Creel Location", "Special Conditions",
          "Check In Time", "Check Out Time", "Notes")
miss <- setdiff(need, names(raw))
if (length(miss)) stop("raw survey export lacks column(s): ", paste(miss, collapse = ", "))

# "H:M:S" or "H:M" text in 24-hour clock time -> decimal hours; anything else NA.
.clock_hours <- function(x) {
  x <- str_squish(as.character(x))
  m <- str_match(x, "^(\\d{1,2}):(\\d{1,2})(?::(\\d{1,2}))?$")
  h <- suppressWarnings(as.numeric(m[, 2])); mi <- suppressWarnings(as.numeric(m[, 3])); s <- suppressWarnings(as.numeric(m[, 4]))
  s[is.na(s)] <- 0
  out <- h + mi / 60 + s / 3600
  out[is.na(h) | is.na(mi) | h > 24 | mi >= 60] <- NA_real_
  out
}
.hhmm <- function(hrs) ifelse(is.na(hrs), NA_character_, sprintf("%02d:%02d", floor(hrs), round((hrs - floor(hrs)) * 60)))
# The export's Date is M/D/YYYY text (an Excel date cell would come through as POSIXct).
.as_date <- function(x) {
  if (inherits(x, "POSIXct")) return(as.Date(x, tz = "UTC"))
  if (inherits(x, "Date")) return(x)
  d <- suppressWarnings(as.Date(as.character(x), format = "%m/%d/%Y"))
  d2 <- suppressWarnings(as.Date(as.character(x), format = "%Y-%m-%d"))
  d[is.na(d)] <- d2[is.na(d)]
  d
}
.season_of <- function(d) {
  y <- as.integer(format(d, "%Y")); after <- format(d, "%m%d") >= "0916"
  y0 <- ifelse(after, y, y - 1L)
  sprintf("%d-%02d", y0, (y0 + 1L) %% 100)
}

shifts <- raw |>
  transmute(
    survey_num  = suppressWarnings(as.integer(ID)),
    survey_id   = ifelse(is.na(survey_num), NA_character_, paste0("S", survey_num)),
    date_d      = .as_date(Date),
    season      = .season_of(date_d),
    date        = format(date_d, "%Y-%m-%d"),
    day_of_week = as.character(`Day of Week`),
    holiday     = suppressWarnings(as.integer(Holiday)),
    creel_location     = as.character(`Creel Location`),
    samplers           = as.character(Samplers),
    special_conditions = as.character(`Special Conditions`),
    check_in_raw  = as.character(`Check In Time`),
    check_out_raw = as.character(`Check Out Time`),
    check_in_hour  = .clock_hours(`Check In Time`),
    check_out_hour = .clock_hours(`Check Out Time`),
    notes = as.character(Notes)
  ) |>
  mutate(
    qc_flag = case_when(
      is.na(check_in_hour)                              ~ "missing_check_in",
      is.na(check_out_hour)                             ~ "missing_check_out",
      check_out_hour < check_in_hour                    ~ "check_out_before_check_in",
      TRUE ~ ""),
    qc_flag = ifelse(!is.na(survey_id) & duplicated(survey_id) | duplicated(survey_id, fromLast = TRUE) & !is.na(survey_id),
                     paste0(qc_flag, ifelse(nzchar(qc_flag), ";", ""), "duplicate_survey_id"), qc_flag),
    shift_hours    = ifelse(qc_flag == "", check_out_hour - check_in_hour, NA_real_),
    check_in       = .hhmm(check_in_hour),
    check_out      = .hhmm(check_out_hour)
  ) |>
  select(survey_id, survey_num, season, date, day_of_week, holiday, creel_location, samplers, special_conditions,
         check_in, check_out, check_in_hour, check_out_hour, shift_hours, qc_flag, notes) |>
  arrange(date, survey_num)

if (any(is.na(shifts$date))) warning(sprintf("%d row(s) with an unparseable date", sum(is.na(shifts$date))))
cat("Surveys by season x creel_location:\n"); print(table(shifts$season, shifts$creel_location))
cat("\nQC flags:\n"); print(table(shifts$qc_flag, useNA = "ifany"))
fl <- shifts |> filter(qc_flag != "")
if (nrow(fl)) { cat("flagged rows:\n"); print(as.data.frame(fl |> select(survey_id, date, creel_location, check_in, check_out, qc_flag))) }
gh <- shifts |> filter(creel_location == "Grays Harbor", !is.na(shift_hours))
cat(sprintf("\nGrays Harbor: %d surveys with a usable shift; check-in median %s, check-out median %s, shift hours median %.2f\n",
            nrow(gh), .hhmm(median(gh$check_in_hour)), .hhmm(median(gh$check_out_hour)), median(gh$shift_hours)))

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(list(data = shifts), out_path)
cat(sprintf("\nWrote %s: %d rows x %d columns.\n", out_path, nrow(shifts), ncol(shifts)))
