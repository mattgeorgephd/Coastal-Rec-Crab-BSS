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
# build_interview_combined.R  (2026-09-09; rebuilt on the season workbooks 2026-09-10)
#
# Builds 04_input_files/interview_combined.xlsx (sheet "data") from the "Harvest Creel"
# sheet of every per-season creel workbook in raw/ (see build_helpers.R for why the
# season workbooks replaced the pasted export). One row per interview, every season,
# through the last survey day in the newest workbook.
#
# WHAT IT FIXES AGAINST THE 2026-09-09 WORKBOOK
#   - 2023-24 gear count. The 2023-24 sheet heads the count "Number Of Gear" (capital O);
#     the pasted export had carried it as a second column under that header and the
#     2026-09-09 builder read only "Number of Gear", so number_of_gear was NA for the
#     whole 2023-24 season (13,629 rows) in both the July and the September workbooks.
#     The builder now takes the count from either header; 2023-24 has 13,629 counts,
#     0 to 23, matching the season workbook exactly (CHANGE_REGISTER D15, closed).
#   - Coverage. 2025-26 extends from 2026-08-01 (8,490 rows) to 2026-09-08 (9,992 rows).
#   - Gear vocabulary. The iForm's gear labels changed for 2024-25: "Collapsible trap or
#     ring" became "Ring Net", "Fishing rod with foldable trap" / "Star trap" became
#     "Trap (foldable, star)", "Fishing rod with snare" became "Snare", "Rake, net or
#     hands" became "Rake or Net" (the 2023-24 workbook's gear_key sheet). The gear-
#     resolved classifier (prep_bss_crab_gear.R) matches "ring net", "trap|star",
#     "snare", "pot" by regex, so on the old labels a ring net was classed as a TRAP.
#     gear_type now carries the 2024-25 vocabulary for every season (gear_type_raw keeps
#     the label as recorded); for 2024-25 and 2025-26 rows the two are identical.
#
# COLUMNS. The 21 columns of 2026-09-09 keep their names, order, types and cleaning
# (ISO-text dates, numeric counts, completed_trip 0/1, gear_tampered 0/1, interview_time
# "HH:MM"). Added on 2026-09-10, each with a modelling use recorded in ../README.md:
#   gear_type_raw, boat_name, bay_or_ocean, river_or_ocean (where the boat fished),
#   total_vehicles (2024-25 on), crab_released (0/1), dungeness_returned,
#   dungeness_returned_reason, red_rock_returned, red_rock_returned_reason.
# Notes (free text with vehicle and personal descriptions) stay in raw/.
#
# CLEANING RULES (applied identically to every season; every non-parsing value is NA
# and counted in the report):
#   - creel_area trimmed (the source carries "Cape Disappointment State Park Boat Launch "
#     with a trailing space) and the 2022-23 site names mapped to the 2023-24+ vocabulary
#     (area_map below; no Westport site is affected); creel_location from the survey
#     sheet by survey id.
#   - crabbing_mode / boat_type: case normalised ("boat" -> "Boat", "private" ->
#     "Private"); the source's "Commerical" spelling is kept (the reader matches it).
#   - completed_trip: 0/1 only; any other value (one "24" in 2023-24) is NA.
#   - a formula fragment left in a count cell ("3+N6847:N13622" in 2023-24 Dungeness
#     kept) is NA.
#   - season is the workbook's; a row whose date falls outside that fishery season
#     stops the build (a mis-filed row would otherwise change the season it is read in).
#
# USAGE (from the repository root): Rscript 04_input_files/build_interview_combined.R
# Prints the rows by season and location, the trip-type table for the Grays Harbor boat
# interviews, the gear-label harmonisation table, and a comparison with the previous
# workbook. Re-run when a season workbook is replaced (NEW_SEASON_GUIDE.md, step 1).
###############################################################################

source(file.path(if (dir.exists("04_input_files")) "04_input_files" else ".", "build_helpers.R"))
root <- build_root()
out_path <- file.path(root, "04_input_files", "interview_combined.xlsx")

# ---- trip-type classes (unchanged from 2026-09-09) --------------------------------------
# crab_only     the boat crabbed and did nothing else
# combo         the boat crabbed AND fished another fishery ("<fishery> & Crab")
# other_fishery the boat fished but did not crab ("Finfish Only")
# non_fishing   a launch with no fishing at all; counted by the trailer and OSP totals,
#               so a not-crabbing boat to f
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

# ---- gear vocabulary (2023-24 workbook, sheet gear_key; 2024-25 labels are the standard) --
gear_map <- c(
  "pot" = "Pot", "pots" = "Pot", "ring pot" = "Pot",
  "collapsible trap or ring" = "Ring Net", "ring net" = "Ring Net",
  "fishing rod with foldable trap" = "Trap (foldable, star)", "fishing rod with collapsible trap" = "Trap (foldable, star)",
  "star trap" = "Trap (foldable, star)", "trap (foldable, star)" = "Trap (foldable, star)",
  "fishing rod with snare" = "Snare", "snare" = "Snare",
  "slip ring pot" = "Slip Ring Pot",
  "rake, net or hands" = "Rake or Net", "rake, net or hand" = "Rake or Net", "rake or net" = "Rake or Net",
  "handline" = "Handline", "other" = "Other", "none" = "None", "unknown" = "Unknown", "n/a" = "Unknown", "multiple" = "Multiple"
)
# tokenise a comma-separated label list by matching the known labels (longest first, so
# "Trap (foldable, star)" and "Rake, net or hands" survive their own commas), map each,
# keep the source order; an unknown fragment is kept verbatim and reported.
gear_tokens <- function(s) {
  keys <- names(gear_map)[order(-nchar(names(gear_map)))]
  if (is.na(s) || !nzchar(str_squish(s))) return(character(0))
  rest <- str_squish(s); out <- character(0)
  while (nzchar(rest)) {
    rest <- sub("^[,\\s]+", "", rest, perl = TRUE); if (!nzchar(rest)) break
    low <- tolower(rest); hit <- NA_character_
    for (k in keys) if (startsWith(low, k) && (nchar(low) == nchar(k) || grepl("^[,\\s]", substr(low, nchar(k) + 1, nchar(k) + 1), perl = TRUE))) { hit <- k; break }
    if (!is.na(hit)) { out <- c(out, gear_map[[hit]]); rest <- substr(rest, nchar(hit) + 1, nchar(rest)) }
    else { frag <- sub(",.*$", "", rest); out <- c(out, str_squish(frag)); rest <- substr(rest, nchar(frag) + 1, nchar(rest)) }
  }
  out
}
harmonise_gear <- function(x) {
  vapply(as.character(x), function(s) { t <- gear_tokens(s); if (!length(t)) NA_character_ else paste(t, collapse = ", ") },
         character(1), USE.NAMES = FALSE)
}

# ---- read every season ----------------------------------------------------------------
wb <- season_workbooks(root)
say("Season workbooks: %s", paste(sprintf("%s (%s)", basename(wb$path), wb$season), collapse = ", "))

read_season <- function(path, season) {
  hs <- sheet_named(path, c("Harvest Creel")); ss <- sheet_named(path, c("crab creel survey data"))
  if (is.null(hs)) stop("no 'Harvest Creel' sheet in ", basename(path))
  h <- read_sheet_text(path, hs)
  loc <- if (!is.null(ss)) {
    s <- read_sheet_text(path, ss)
    tibble(survey_id = survey_id_of(pick_col(s, "ID")), creel_location = str_squish(pick_col(s, "Creel Location"))) |>
      filter(!is.na(survey_id)) |> distinct(survey_id, .keep_all = TRUE)
  } else tibble(survey_id = character(), creel_location = character())

  date_raw <- pick_col(h, "Date"); date_d <- parse_date_any(date_raw)
  report_unparsed(sprintf("%s date", season), date_raw, date_d)
  itime_raw <- pick_col(h, "Interview Time"); itime <- parse_clock_hours(itime_raw)
  report_unparsed(sprintf("%s interview time", season), itime_raw, itime)
  dk_raw <- pick_col(h, "Number of Dungeness Crab Kept"); dk <- num(dk_raw); report_unparsed(sprintf("%s dungeness kept", season), dk_raw, dk)
  ct_raw <- num(pick_col(h, "Completed Fishing Trip?")); ct <- ifelse(ct_raw %in% c(0, 1), ct_raw, NA_real_)
  if (any(!is.na(ct_raw) & !ct_raw %in% c(0, 1))) say("  NOTE %s completed trip: %d value(s) other than 0/1 set to NA (%s)", season,
                                                     sum(!is.na(ct_raw) & !ct_raw %in% c(0, 1)), paste(unique(ct_raw[!is.na(ct_raw) & !ct_raw %in% c(0, 1)]), collapse = ", "))
  gt_raw <- num(pick_col(h, "Crabber states their pots were ran by someone else"))
  cr_raw <- str_squish(pick_col(h, c("Crab Released?", "Crab released?")))
  rr_ret_raw <- pick_col(h, "Red Rock Returned"); rr_ret <- num(rr_ret_raw); report_unparsed(sprintf("%s red rock returned", season), rr_ret_raw, rr_ret)
  gear_raw <- str_squish(pick_col(h, "Gear Type")); gear_raw[gear_raw == ""] <- NA_character_

  tibble(
    season         = season,
    survey_id      = survey_id_of(pick_col(h, "ID")),
    date_d         = date_d,
    date           = format(date_d, "%Y-%m-%d"),
    interview_num  = num(pick_col(h, "Interview #")),
    crabbing_mode  = str_to_title(str_squish(pick_col(h, "Crabbing Mode"))),
    boat_type      = { b <- str_squish(pick_col(h, "Boat Type")); ifelse(tolower(b) == "private", "Private", b) },
    crabbers       = num(pick_col(h, "Crabbers")),
    gear_type_raw  = gear_raw,
    gear_type      = harmonise_gear(gear_raw),
    number_of_gear = num(pick_col(h, c("Number of Gear", "Number Of Gear"))),
    dungeness_kept = dk,
    red_rock_kept  = num(pick_col(h, "Number of Red Rock Crab Kept")),
    hours_fished   = num(pick_col(h, "Number of Hours Fished")),
    crabber_hours  = num(pick_col(h, "Crabber Hours")),
    gear_hours     = num(pick_col(h, "Gear Hours")),
    completed_trip = ct,
    gear_tampered  = ifelse(is.na(gt_raw), NA_real_, ifelse(gt_raw != 0, 1, 0)),
    trip_type      = { t <- str_squish(pick_col(h, "Trip Type")); ifelse(is.na(t) | t == "", NA_character_, t) },
    creel_area     = harmonise_area(pick_col(h, "Creel Area")),
    interview_time = hhmm(itime),
    boat_name      = { b <- str_squish(pick_col(h, "Boat Name")); ifelse(is.na(b) | b == "", NA_character_, b) },
    bay_or_ocean   = { b <- str_squish(pick_col(h, "Bay or Ocean")); ifelse(is.na(b) | b == "", NA_character_, b) },
    river_or_ocean = { b <- str_squish(pick_col(h, "River or Ocean")); ifelse(is.na(b) | b == "", NA_character_, b) },
    total_vehicles = num(pick_col(h, "Total Vehicles")),
    crab_released  = case_when(is.na(cr_raw) | cr_raw == "" ~ NA_real_, str_detect(cr_raw, "(?i)^y") ~ 1, str_detect(cr_raw, "(?i)^n") ~ 0, TRUE ~ NA_real_),
    dungeness_returned        = num(pick_col(h, "Dungeness Returned")),
    dungeness_returned_reason = { b <- str_squish(pick_col(h, "Dungeness Returned Reason")); ifelse(is.na(b) | b == "", NA_character_, b) },
    red_rock_returned         = rr_ret,
    red_rock_returned_reason  = { b <- str_squish(pick_col(h, "Red Rock Returned Reason")); ifelse(is.na(b) | b == "", NA_character_, b) }
  ) |>
    mutate(trip_type_class = trip_type_class(trip_type)) |>
    left_join(loc, by = "survey_id")
}

combined <- pmap_dfr(list(wb$path, wb$season), read_season)

# ---- integrity checks ----------------------------------------------------------------
bad_season <- combined |> filter(!is.na(date_d), season_of(date_d) != season)
if (nrow(bad_season)) stop(sprintf("%d row(s) dated outside their workbook's fishery season, e.g. %s in the %s workbook",
                                   nrow(bad_season), bad_season$date[1], bad_season$season[1]))
if (any(is.na(combined$date))) say("  NOTE %d row(s) without a parseable date", sum(is.na(combined$date)))
if (any(is.na(combined$creel_location))) say("  NOTE %d row(s) whose survey id has no survey-sheet row (creel_location NA)", sum(is.na(combined$creel_location)))
if (any(combined$trip_type_class %in% "unclassified"))
  warning(sprintf("%d trip-type label(s) not recognised: %s", sum(combined$trip_type_class %in% "unclassified"),
                  paste(unique(combined$trip_type[combined$trip_type_class %in% "unclassified"]), collapse = " | ")))
# gear labels the map did not recognise
tok_raw <- unlist(lapply(na.omit(unique(combined$gear_type_raw)), gear_tokens))
unknown_gear <- setdiff(unique(tok_raw), unique(unname(gear_map)))
if (length(unknown_gear)) say("  NOTE gear label fragment(s) not in the vocabulary, kept verbatim: %s", paste(unknown_gear, collapse = " | "))

# ---- the reader contract: the 21 columns of 2026-09-09 first, in their order ----------------
legacy <- c("season", "creel_location", "date", "survey_id", "interview_num", "crabbing_mode", "boat_type", "crabbers",
            "gear_type", "number_of_gear", "dungeness_kept", "red_rock_kept", "hours_fished", "crabber_hours",
            "gear_hours", "completed_trip", "gear_tampered", "trip_type", "trip_type_class", "creel_area", "interview_time")
added  <- c("gear_type_raw", "boat_name", "bay_or_ocean", "river_or_ocean", "total_vehicles", "crab_released",
            "dungeness_returned", "dungeness_returned_reason", "red_rock_returned", "red_rock_returned_reason")
combined <- combined |> arrange(season, date, survey_id, interview_num) |> select(all_of(c(legacy, added)))

# ---- report ----------------------------------------------------------------------------
cat("\nRows by season x creel_location:\n"); print(table(combined$season, combined$creel_location, useNA = "ifany"))
cat(sprintf("\nDate range by season:\n")); print(combined |> group_by(season) |> summarise(first = min(date), last = max(date), rows = n(), .groups = "drop") |> as.data.frame())
cat("\nnumber_of_gear present (non-NA) by season:\n"); print(tapply(!is.na(combined$number_of_gear), combined$season, sum))
cat("\nGrays Harbor BOAT interviews: trip-type class by season:\n")
gb <- combined |> filter(creel_location == "Grays Harbor", crabbing_mode == "Boat")
print(table(gb$season, gb$trip_type_class, useNA = "ifany"))
cat("\nGear labels harmonised (raw -> standard), single-label rows by season (multi-gear lists follow the same map):\n")
ch <- combined |> filter(!is.na(gear_type_raw), gear_type_raw != gear_type) |> count(season, gear_type_raw, gear_type) |> arrange(season, desc(n))
for (i in seq_len(nrow(ch))) if (!grepl(",", ch$gear_type[i])) say("  %s  %-36s -> %-24s %5d", ch$season[i], ch$gear_type_raw[i], ch$gear_type[i], ch$n[i])
say("  rows whose gear label changed: %s", paste(sprintf("%s=%d", names(table(ch$season)), tapply(ch$n, ch$season, sum)), collapse = ", "))
cat("\ngear_tampered = 1 rows by season x location:\n"); print(table(combined$season[combined$gear_tampered %in% 1], combined$creel_location[combined$gear_tampered %in% 1]))

# ---- comparison with the previous workbook, when present -----------------------------
if (file.exists(out_path)) {
  old <- suppressWarnings(readxl::read_excel(out_path, sheet = "data"))
  cat("\nComparison with the previous workbook on (season, survey_id, interview_num), rows with an interview number:\n")
  keys <- c("season", "survey_id", "interview_num")
  o <- old |> filter(!is.na(interview_num)) |> mutate(interview_num = num(interview_num), creel_area = str_squish(creel_area)) |> distinct(across(all_of(keys)), .keep_all = TRUE)
  n <- combined |> filter(!is.na(interview_num)) |> distinct(across(all_of(keys)), .keep_all = TRUE)
  compare_with_previous(o, n, keys, setdiff(legacy, keys))
}

write_data_sheet(combined, out_path)
