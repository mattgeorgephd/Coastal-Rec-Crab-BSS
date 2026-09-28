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
# build_comm_charter_tally.R  (2026-09-10)
#
# Builds 04_input_files/wes_commercial_tally.xlsx (sheet "data"): the Westport daily
# vessel tally (private boats, commercial vessels fishing as recreational, charters) on
# the days a sampler was in port during the commercial-as-recreational window. Read by
# estimate_comm_charter() as the census frame of the commercial/charter component.
#
# SOURCES, per season workbook in raw/:
#   - the "wes commercial tally" sheet (2024-25 on): a two-row header ("tally counts" /
#     "actual interviews" over the column names), Date, then the three tallies and the
#     three interviewed counts. THE COLUMN ORDER DIFFERS BETWEEN SEASONS (2024-25:
#     Private, Commercial, Charter; 2025-26: Commercial, Private, Charter), so the
#     builder reads by header name, never by position.
#   - the survey sheet's "Commercial Boat Tally" / "Private Boat Tally" / "Charter Boat
#     Tally" columns (recorded per survey in the field from 2025-26): the same counts,
#     entered by the sampler. The two are reconciled per date; a disagreement is printed
#     and the tally sheet's value is kept (it carries the interviewed counts and the
#     notes); a date present only in the survey sheet is added from there.
# 2022-23 and 2023-24 have neither: no tally was kept before 2024-25, so those seasons
# have no census frame (the 33 Westport commercial-vessel interviews of 2023-24 have no
# vessel count to expand to; estimate_comm_charter() says so at run time).
#
# COLUMNS: season, date (ISO text), private_tally, commercial_tally, charter_tally,
# private_interviewed, commercial_interviewed, charter_interviewed, source ("tally sheet",
# "survey sheet", "both"), notes.
#
# The commercial-as-recreational window each tally covers (for run_config$census_windows):
#   2024-25  2024-12-03 to 2025-02-08 (47 days; the commercial fishery opened Feb 11, 2025)
#   2025-26  2025-12-01 to 2026-01-03 (23 days; Grays Harbor opened Jan 4, 2026, WDFW news
#            release 2025-12-29; the last south-area recreational trips ran Dec 27 and the
#            gear-set days are noted on the sheet)
#
# USAGE (from the repository root): Rscript 04_input_files/build_comm_charter_tally.R
###############################################################################

source(file.path(if (dir.exists("04_input_files")) "04_input_files" else ".", "build_helpers.R"))
root <- build_root()
out_path <- file.path(root, "04_input_files", "wes_commercial_tally.xlsx")

wb <- season_workbooks(root)
txt <- function(x) { x <- str_squish(as.character(x)); ifelse(is.na(x) | x == "", NA_character_, x) }

read_tally_sheet <- function(path, season) {
  ts <- sheet_named(path, c("wes commercial tally"))
  if (is.null(ts)) return(NULL)
  # find the header row: the first row whose first cell is "Date"
  raw <- suppressWarnings(readxl::read_excel(path, sheet = ts, col_types = "text", col_names = FALSE, .name_repair = "minimal"))
  hdr <- which(tolower(str_squish(as.character(raw[[1]]))) == "date")[1]
  if (is.na(hdr)) stop("no 'Date' header row in the tally sheet of ", basename(path))
  t <- read_sheet_text(path, ts, skip = hdr - 1)
  d_raw <- pick_col(t, "Date"); d <- parse_date_any(d_raw)
  t <- t[!is.na(d), , drop = FALSE]; d <- d[!is.na(d)]
  tibble(
    season = season, date_d = d, date = format(d, "%Y-%m-%d"),
    private_tally    = num(pick_col(t, c("Private"))),
    commercial_tally = num(pick_col(t, c("Commercial"))),
    charter_tally    = num(pick_col(t, c("Charter"))),
    private_interviewed    = num(pick_col(t, c("private int"))),
    commercial_interviewed = num(pick_col(t, c("commercial int"))),
    charter_interviewed    = num(pick_col(t, c("charter int"))),
    notes = txt(pick_col(t, "notes"))
  )
}

read_survey_tallies <- function(path, season) {
  ss <- sheet_named(path, c("crab creel survey data"))
  if (is.null(ss)) return(NULL)
  s <- read_sheet_text(path, ss)
  if (!has_col(s, "Commercial Boat Tally")) return(NULL)
  d <- parse_date_any(pick_col(s, "Date"))
  out <- tibble(
    season = season, date_d = d, date = format(d, "%Y-%m-%d"),
    creel_location = txt(pick_col(s, "Creel Location")),
    private_tally    = num(pick_col(s, "Private Boat Tally")),
    commercial_tally = num(pick_col(s, "Commercial Boat Tally")),
    charter_tally    = num(pick_col(s, "Charter Boat Tally"))
  ) |>
    filter(!is.na(private_tally) | !is.na(commercial_tally) | !is.na(charter_tally))
  if (nrow(out) && !all(out$creel_location %in% "Grays Harbor"))
    say("  NOTE %s: %d survey-sheet tally row(s) outside Grays Harbor are ignored", season, sum(!out$creel_location %in% "Grays Harbor"))
  out |> filter(creel_location %in% "Grays Harbor") |> select(-creel_location)
}

build_season <- function(path, season) {
  ts <- read_tally_sheet(path, season); sv <- read_survey_tallies(path, season)
  if (is.null(ts) && is.null(sv)) { say("  %s: no vessel tally in the workbook", season); return(NULL) }
  if (is.null(ts)) return(sv |> mutate(private_interviewed = NA_real_, commercial_interviewed = NA_real_, charter_interviewed = NA_real_, source = "survey sheet", notes = NA_character_))
  ts <- ts |> mutate(source = "tally sheet")
  if (!is.null(sv) && nrow(sv)) {
    j <- full_join(ts, sv, by = c("season", "date_d", "date"), suffix = c("", ".sv"))
    both <- j |> filter(!is.na(source), !is.na(private_tally.sv))
    dis <- both |> filter(private_tally != private_tally.sv | commercial_tally != commercial_tally.sv | charter_tally != charter_tally.sv)
    say("  %s: tally sheet %d days, survey-sheet tallies %d days, %d in both, %d disagreeing", season, nrow(ts), nrow(sv), nrow(both), nrow(dis))
    if (nrow(dis)) print(as.data.frame(dis |> select(date, private_tally, private_tally.sv, commercial_tally, commercial_tally.sv, charter_tally, charter_tally.sv)))
    only_sv <- j |> filter(is.na(source))
    if (nrow(only_sv)) {
      say("  %s: %d day(s) only in the survey sheet are added from it: %s", season, nrow(only_sv), paste(only_sv$date, collapse = ", "))
      add <- only_sv |> transmute(season, date_d, date, private_tally = private_tally.sv, commercial_tally = commercial_tally.sv, charter_tally = charter_tally.sv,
                                  private_interviewed = NA_real_, commercial_interviewed = NA_real_, charter_interviewed = NA_real_, notes = NA_character_, source = "survey sheet")
      ts <- bind_rows(ts, add)
    }
    ts$source[ts$date %in% both$date] <- "both"
  } else say("  %s: tally sheet %d days (no survey-sheet tallies)", season, nrow(ts))
  ts
}

tally <- pmap(list(wb$path, wb$season), build_season) |> compact() |> bind_rows() |>
  arrange(date) |>
  select(season, date, private_tally, commercial_tally, charter_tally, private_interviewed, commercial_interviewed, charter_interviewed, source, notes)

if (any(duplicated(tally$date))) stop("duplicated tally dates: ", paste(tally$date[duplicated(tally$date)], collapse = ", "))
cat("\nTally days by season (window, vessels):\n")
print(as.data.frame(tally |> group_by(season) |> summarise(days = n(), first = min(date), last = max(date), private = sum(private_tally, na.rm = TRUE),
                                                          commercial = sum(commercial_tally, na.rm = TRUE), charter = sum(charter_tally, na.rm = TRUE), .groups = "drop")))

if (file.exists(out_path)) {
  old <- suppressWarnings(readxl::read_excel(out_path, sheet = "data", guess_max = 1e5))
  cat("\nComparison with the previous workbook on date:\n")
  compare_with_previous(old, tally, "date", c("season", "private_tally", "commercial_tally", "charter_tally"))
}
write_data_sheet(tally, out_path)
