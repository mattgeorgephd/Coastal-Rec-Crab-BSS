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
# read_crabbing_holidays.R
#
# Single-source reader for the crabbing-holiday calendar. Holidays used to be a
# hardcoded as.Date(c(...)) vector in run_config.R (and a duplicated copy in the
# weather module); they now live in an editable workbook,
# 04_input_files/<params$crabbing_holidays_file> (default crabbing_holidays.xlsx,
# sheet params$crabbing_holidays_sheet, default "data"), so the season-to-season
# update is a spreadsheet edit and multiple seasons can coexist in one file.
#
# Expected columns (case-insensitive, trimmed): season, date, holiday_name.
#   season        matches params$season_filter (e.g. "2024-25")
#   date          ISO yyyy-mm-dd (or any format base::as.Date parses by default)
#   holiday_name  free text (documentation only; not read downstream)
#
# Returns a sorted, de-duplicated Date vector of the holidays for
# params$season_filter. Fails LOUDLY (stop) if the file, the required columns, or
# the requested season are missing, so a mistyped season can never silently blank
# out the holiday day-typing (that would misclassify high-effort holidays as
# ordinary weekdays and bias the effort expansion).
#
# Auto-sourced by both drivers via the 03_R_functions walk; requires readxl + here
# (both loaded in every driver's setup chunk before this is called).
###############################################################################

read_crabbing_holidays <- function(params) {
  file  <- params$crabbing_holidays_file  %||% "crabbing_holidays.xlsx"
  sheet <- params$crabbing_holidays_sheet %||% "data"
  path  <- here::here("04_input_files", file)

  if (!file.exists(path))
    stop("Crabbing-holiday workbook not found: ", path,
         "\n  Add the file, or point run_config$crabbing_holidays_file at it.",
         call. = FALSE)

  hol <- readxl::read_excel(path, sheet = sheet)
  names(hol) <- tolower(trimws(names(hol)))

  if (!all(c("season", "date") %in% names(hol)))
    stop("Crabbing-holiday workbook must have 'season' and 'date' columns; got: ",
         paste(names(hol), collapse = ", "), call. = FALSE)

  season <- as.character(params$season_filter)
  if (is.null(season) || !nzchar(season))
    stop("params$season_filter is unset; cannot select crabbing holidays.", call. = FALSE)

  hol_season <- hol[trimws(as.character(hol$season)) == season, , drop = FALSE]
  if (nrow(hol_season) == 0)
    stop("No crabbing holidays for season '", season, "' in ", path,
         ".\n  Add rows for this season (columns: season, date, holiday_name), ",
         "or fix season_filter.", call. = FALSE)

  dates <- as.Date(hol_season$date)
  if (any(is.na(dates)))
    stop("Unparseable date(s) in ", path, " for season '", season, "': ",
         paste(hol_season$date[is.na(dates)], collapse = ", "), call. = FALSE)

  sort(unique(dates))
}
