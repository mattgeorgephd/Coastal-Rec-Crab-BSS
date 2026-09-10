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
# build_charter_trips.R  (2026-09-10)
#
# Builds 04_input_files/charter_trips.xlsx (sheet "data") from the "charter trips" sheet
# of every season workbook in raw/ that has one (2024-25 on): the CHARTER TRIP ROSTER,
# one row per charter crab trip the operators reported, with whether the samplers
# interviewed it ("interviewed"), did not ("missed"), or the trip did not sail
# ("canceled"). Ports: Westport, Ilwaco, Chinook, Cape D.
#
# WHY IT IS AN INPUT. The commercial/charter census (estimate_comm_charter.R) counts
# vessels from the daily tally, which exists only on the days a sampler was in port. The
# roster is a trip-level frame that does not depend on a sampler being present: on
# 2024-25 it lists 31 Westport charter trips in the census window that sailed, 8 of them
# (all "missed") on days WITHOUT a tally, which the exact-sum census (census_expansion =
# "none", "unsampled days = no operation") counts as zero. run_config$charter_frame =
# "roster" reads this file and takes, per day, the larger of the roster's trips and the
# tally's charter count (the union of the two frames); "tally" ignores it.
#
# The sheet also carries the operators' own summary (interviewed / missed / total trips,
# observed crab, an estimate of the unobserved crab as trips x mean per interviewed trip),
# which is the roster method by hand; those summary cells are not read.
#
# COLUMNS: season, port (the sheet's Location, trimmed), vessel (the sheet's Charter
# Boat, trimmed and title-cased so "Coho charters" and "Coho Charters" are one vessel),
# date (ISO text), status ("interviewed" | "missed" | "canceled"), contact ("phone" |
# "in person", 2025-26 notes), notes.
#
# USAGE (from the repository root): Rscript 04_input_files/build_charter_trips.R
###############################################################################

source(file.path(if (dir.exists("04_input_files")) "04_input_files" else ".", "build_helpers.R"))
root <- build_root()
out_path <- file.path(root, "04_input_files", "charter_trips.xlsx")

wb <- season_workbooks(root)
txt <- function(x) { x <- str_squish(as.character(x)); ifelse(is.na(x) | x == "", NA_character_, x) }

read_roster <- function(path, season) {
  cs <- sheet_named(path, c("charter trips"))
  if (is.null(cs)) { say("  %s: no charter trips sheet", season); return(NULL) }
  r <- read_sheet_text(path, cs)
  need <- c("Location", "Charter Boat", "Date", "Sampled")
  if (!all(vapply(need, function(a) has_col(r, a), logical(1)))) stop("charter trips sheet of ", basename(path), " lacks ", paste(need, collapse = "/"))
  d <- parse_date_any(pick_col(r, "Date"))
  status <- tolower(txt(pick_col(r, "Sampled")))
  notes <- txt(pick_col(r, c("Notes", "Unnamed: 4")))
  # the 2024-25 sheet has an unnamed fifth column with the occasional note
  if (all(is.na(notes)) && ncol(r) >= 5 && !has_col(r, "Notes")) notes <- txt(r[[5]])
  out <- tibble(
    season = season, date_d = d, date = format(d, "%Y-%m-%d"),
    port   = txt(pick_col(r, "Location")),
    vessel = str_to_title(txt(pick_col(r, "Charter Boat"))),
    status = case_when(is.na(status) ~ NA_character_,
                       str_detect(status, "^interview") ~ "interviewed",
                       str_detect(status, "^miss") ~ "missed",
                       str_detect(status, "^cancel") ~ "canceled",
                       TRUE ~ status),
    contact = case_when(is.na(notes) ~ NA_character_, str_detect(notes, "(?i)^phone") ~ "phone",
                        str_detect(notes, "(?i)^in person") ~ "in person", TRUE ~ NA_character_),
    notes = notes
  ) |> filter(!is.na(date_d), !is.na(port))
  bad <- out |> filter(!status %in% c("interviewed", "missed", "canceled"))
  if (nrow(bad)) say("  NOTE %s: %d row(s) with an unrecognised status kept verbatim: %s", season, nrow(bad), paste(unique(bad$status), collapse = " | "))
  if (any(season_of(out$date_d) != season)) stop("charter trip dated outside the ", season, " workbook's season")
  out
}

roster <- pmap(list(wb$path, wb$season), read_roster) |> compact() |> bind_rows() |>
  arrange(season, port, date, vessel) |> select(season, port, vessel, date, status, contact, notes)

cat("\nCharter trips by season x port x status:\n"); print(ftable(table(roster$season, roster$port, roster$status)))
cat("\nWestport trips by vessel and season:\n"); print(table(roster$vessel[roster$port == "Westport"], roster$season[roster$port == "Westport"]))
write_data_sheet(roster, out_path)
