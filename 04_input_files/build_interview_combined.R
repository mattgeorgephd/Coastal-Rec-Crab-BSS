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
# build_interview_combined.R  (2026-09-09)
#
# Rebuilds 04_input_files/interview_combined.xlsx (sheet "data") from the raw creel
# database export 04_input_files/raw/interviewdata<YYYY><YYYY>.xlsx, so that the model
# workbook is DERIVED from the export by one documented rule rather than hand-edited.
#
# WHY. The 2026-09-09 export (2022-23 through 2025-26, 37,460 rows) carries three fields
# the pipeline now needs and the 2026-07-16 trimmed workbook did not:
#   Trip Type        crab-only / <fishery> & Crab (combo) / Finfish Only / Non Fishing Trip.
#                    The crabbing fraction f and the combo-trip share c are read from it
#                    (03_R_functions/fetch_crab_data.R, crab_fraction.R).
#   creel_area       which site the interview was taken at (the boat launch vs the marina
#                    or the docks), so the f classification can be restricted to the
#                    boats the trailer and OSP counts measure.
#   Interview Time   the clock time of the contact, which places every contact inside its
#                    sampler shift (sampler_shifts.xlsx) and gives the within-shift
#                    return-time distribution by trip type.
# It also fills gear_tampered ("Crabber states their pots were ran by someone else"),
# which the trimmed workbook carried as an all-blank column, so the reader's tampered
# filter now removes the flagged rows (133 across all locations and seasons).
#
# WHAT IS PRESERVED. The 17 legacy columns keep their names, order, types and cleaning
# (ISO-text dates, numeric counts, completed_trip as 0/1). Checked on 2026-09-09 against
# the previous workbook row by row (33,083 rows, same order): every difference is either
# floating-point noise on crabber_hours / gear_hours (identical to 6 decimals) or a
# correction that had been made in the source database since the July export (a few
# completed_trip and boat_type fills, seven crabbers / gear_type corrections outside
# Grays Harbor). The export is the newer truth.
#
# USAGE (from the repository root):
#   Rscript 04_input_files/build_interview_combined.R [raw_export.xlsx] [output.xlsx]
# Defaults: raw/interviewdata20222026.xlsx -> interview_combined.xlsx. Prints the row
# counts by season and location, the trip-type table for the Grays Harbor boat
# interviews, and a column-by-column comparison with the previous workbook when one
# exists. Re-run it whenever a new export arrives (NEW_SEASON_GUIDE.md, step 1).
###############################################################################

suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tibble); library(stringr); library(writexl) })

args <- commandArgs(trailingOnly = TRUE)
.root <- getwd()
if (!dir.exists(file.path(.root, "04_input_files")) && dir.exists(file.path(.root, "..", "04_input_files")))
  .root <- normalizePath(file.path(.root, ".."))
raw_path <- if (length(args) >= 1) args[1] else file.path(.root, "04_input_files", "raw", "interviewdata20222026.xlsx")
out_path <- if (length(args) >= 2) args[2] else file.path(.root, "04_input_files", "interview_combined.xlsx")
stopifnot(file.exists(raw_path))

num <- function(x) suppressWarnings(as.numeric(x))

# ---- the trip-type classes the model reads --------------------------------------
# crab_only     the boat crabbed and did nothing else
# combo         the boat crabbed AND fished another fishery ("<fishery> & Crab"); it is a
#               crabbing boat to f, and the boat OSP's crabbing-only column does not see
# other_fishery the boat fished but did not crab ("Finfish Only")
# non_fishing   a launch with no fishing at all (pleasure trip); counted by the trailer
#               and OSP totals, so it is a not-crabbing boat to f
trip_type_class <- function(x) {
  x <- str_squish(as.character(x))
  dplyr::case_when(
    is.na(x) | x == ""                       ~ NA_character_,
    str_detect(x, "(?i)^crab only$")         ~ "crab_only",
    str_detect(x, "(?i)&\\s*crab")           ~ "combo",
    str_detect(x, "(?i)finfish only|fishing only|other fisher") ~ "other_fishery",
    str_detect(x, "(?i)non.?fishing")        ~ "non_fishing",
    TRUE                                     ~ "unclassified"
  )
}

raw <- suppressWarnings(readxl::read_excel(raw_path, sheet = "data", .name_repair = "minimal"))
need <- c("season", "creel_location", "date", "ID", "Interview #", "Crabbing Mode", "Boat Type", "Crabbers",
          "Gear Type", "Number of Gear", "Number of Dungeness Crab Kept", "Number of Red Rock Crab Kept",
          "Number of Hours Fished", "Crabber Hours", "Gear Hours", "Completed Fishing Trip?",
          "Crabber states their pots were ran by someone else", "Trip Type", "creel_area", "Interview Time")
miss <- setdiff(need, names(raw))
if (length(miss)) stop("raw export lacks column(s): ", paste(miss, collapse = ", "))

# The export's date is an Excel datetime; readxl returns it as UTC POSIXct, so the day is
# taken in UTC (a local-time conversion would shift the day at the boundary).
.gt <- raw[["Crabber states their pots were ran by someone else"]]
combined <- raw |>
  transmute(
    season         = as.character(season),
    creel_location = as.character(creel_location),
    date           = format(as.Date(date, tz = "UTC"), "%Y-%m-%d"),
    survey_id      = as.character(ID),
    interview_num  = num(`Interview #`),
    crabbing_mode  = as.character(`Crabbing Mode`),
    boat_type      = as.character(`Boat Type`),
    crabbers       = num(Crabbers),
    gear_type      = as.character(`Gear Type`),
    number_of_gear = num(`Number of Gear`),        # the intended gear COUNT (see the README:
                                                   # "Number Of Gear" is a different, 0/1 field
                                                   # and is deliberately not used)
    dungeness_kept = num(`Number of Dungeness Crab Kept`),
    red_rock_kept  = num(`Number of Red Rock Crab Kept`),
    hours_fished   = num(`Number of Hours Fished`),
    crabber_hours  = num(`Crabber Hours`),
    gear_hours     = num(`Gear Hours`),
    completed_trip = num(`Completed Fishing Trip?`),
    gear_tampered  = ifelse(is.na(.gt), NA_real_, ifelse(as.logical(.gt), 1, 0)),
    # --- 2026-09-09 additions ---
    trip_type       = as.character(`Trip Type`),
    trip_type_class = trip_type_class(`Trip Type`),
    creel_area      = as.character(creel_area),
    interview_time  = ifelse(is.na(`Interview Time`), NA_character_, format(`Interview Time`, "%H:%M", tz = "UTC"))
  )
if (any(combined$trip_type_class %in% "unclassified"))
  warning(sprintf("%d trip-type label(s) not recognised: %s", sum(combined$trip_type_class %in% "unclassified"),
                  paste(unique(combined$trip_type[combined$trip_type_class %in% "unclassified"]), collapse = " | ")))

cat("Rows by season x creel_location:\n"); print(table(combined$season, combined$creel_location))
cat("\nGrays Harbor BOAT interviews: trip-type class by season:\n")
gb <- combined |> filter(creel_location == "Grays Harbor", tolower(crabbing_mode) == "boat")
print(table(gb$season, gb$trip_type_class, useNA = "ifany"))
cat("\ngear_tampered = 1 rows by season x location:\n"); print(table(combined$season[combined$gear_tampered %in% 1], combined$creel_location[combined$gear_tampered %in% 1]))

# ---- comparison with the previous workbook, when present ------------------------
if (file.exists(out_path)) {
  old <- suppressWarnings(readxl::read_excel(out_path, sheet = "data"))
  legacy <- c("season", "creel_location", "date", "survey_id", "interview_num", "crabbing_mode", "boat_type", "crabbers",
              "gear_type", "number_of_gear", "dungeness_kept", "red_rock_kept", "hours_fished", "crabber_hours",
              "gear_hours", "completed_trip", "gear_tampered")
  n <- min(nrow(old), nrow(combined))
  cat(sprintf("\nPrevious workbook: %d rows; rebuilt: %d rows. Column comparison on the first %d rows (tolerance 1e-6):\n",
              nrow(old), nrow(combined), n))
  for (cn in intersect(legacy, names(old))) {
    a <- old[[cn]][seq_len(n)]; b <- combined[[cn]][seq_len(n)]
    if (is.numeric(b)) { a <- num(a); diff <- !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-6)) }
    else { a <- as.character(a); diff <- !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)) }
    cat(sprintf("  %-16s rows differing: %d\n", cn, sum(diff)))
  }
}

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(list(data = combined), out_path)
cat(sprintf("\nWrote %s: %d rows x %d columns.\n", out_path, nrow(combined), ncol(combined)))
