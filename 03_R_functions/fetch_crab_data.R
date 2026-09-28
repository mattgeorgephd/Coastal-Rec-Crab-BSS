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
# fetch_crab_data.R  (shared: pooled + gear-resolved drivers)
#
# Read and assemble a model's input data (effort counts, interviews, catch,
# commercial tally) from 04_input_files/ and classify interviews by population.
# Pure given params; auto-sourced by both drivers and called by BOTH. This single
# reader replaced the former fetch_crab_data / fetch_crab_data_v2 pair (2026-08-01):
# the two were identical except for one private-boat filter, now controlled by
# params$boat_require_gear_time. Gear-type CPUE classification is NOT done here; it
# is downstream in prep_bss_crab_gear.R, so this function is the same for both tracks.
#
# params$boat_require_gear_time (default TRUE): drop private_boat interviews with no
# positive gear_time_total. TRUE reproduces the historical pooled reader; the
# gear-resolved driver sets FALSE to reproduce the historical fetch_crab_data_v2.
#
# INPUTS (2026-07-16): all inputs are now .xlsx workbooks with a single "data"
# sheet (converted from CSV). Filenames, the sheet name, the interview creel
# location, and the effort site lists are all run_config parameters (with the
# historical values as %||% defaults), so a different port/site set or a renamed
# workbook can be run without editing this function. Dates are stored as ISO
# yyyy-mm-dd text and parsed with as.Date().
###############################################################################

fetch_crab_data <- function(params) {
  cat("Reading data...\n")
  in_sheet <- params$input_sheet %||% "data"
  # Drop private_boat interviews with no positive gear_time_total? TRUE = pooled
  # (historical fetch_crab_data); the gear-resolved driver sets FALSE (historical _v2).
  req_boat_gear_time <- isTRUE(params$boat_require_gear_time %||% TRUE)

  # 2026-09-09: season_filter accepts a character VECTOR so a multi-season span (e.g.
  # c("2024-25", "2025-26") with a matching est_date window) can run. A scalar behaves
  # exactly as before. The calendar indices are already span-safe (prep_days_crab builds
  # sequential year+week / year+month factors), and the one thing a multi-season span
  # CANNOT yet express is more than one pot-closure window; see build_subseasons.R and
  # NEW_SEASON_GUIDE.md.
  effort_raw <- read_input_workbook(params$effort_file %||% "effort_combined.xlsx", sheet = in_sheet) |>
    filter(season %in% params$season_filter) |> mutate(date = as.Date(date))
  # 2026-09-10: the effort workbook carries a qc_flag (04_input_files/build_effort_combined.R).
  # Rows whose flag is in params$effort_qc_drop are held out of the counts. The default drops
  # "gear_count_from_interviews": the Feb 2024 PFD stand-down rows, where the Float 20 / 17-21
  # "count" is the day's TOTAL gear from interviews (already a day's deployments, so feeding it
  # to the count model would inflate the shore effort by about the turnover). The other flags
  # (missing_count_time: dropped below anyway; implausible_count_time: a mistyped clock time
  # on a real count) are informational and stay.
  if ("qc_flag" %in% names(effort_raw)) {
    drop_flags <- params$effort_qc_drop %||% "gear_count_from_interviews"
    if (length(drop_flags) && identical(tolower(drop_flags[1]), "none")) drop_flags <- character(0)   # "none" keeps every row
    .fl <- ifelse(is.na(effort_raw$qc_flag), "", as.character(effort_raw$qc_flag))
    if (any(.fl %in% drop_flags))
      cat(sprintf("  Effort counts: %d row(s) held out by qc_flag (%s)\n", sum(.fl %in% drop_flags),
                  paste(unique(.fl[.fl %in% drop_flags]), collapse = ", ")))
    effort_raw <- effort_raw[!.fl %in% drop_flags, , drop = FALSE]
  }

  interview_raw <- read_input_workbook(params$interview_file %||% "interview_combined.xlsx", sheet = in_sheet) |>
    mutate(completed_trip = as.character(completed_trip)) |>
    filter(season %in% params$season_filter)

  # gear_tampered flags interviews where the crabber believes a third party
  # pulled their pots; it may be absent in older workbooks, so ensure it exists.
  if (!"gear_tampered" %in% names(interview_raw)) interview_raw$gear_tampered <- NA_real_

  gh_effort <- effort_raw |> filter(creel_area %in% (params$gh_effort_areas %||% c(
    "Westport Docks Float 20","Westport Docks Float 17-21",
    "Westport Boat Launch","Westport Marina","Westport Jetty",
    "Ocean Shores Boat Launch","Damon Point")))

  gh_interview <- interview_raw |> filter(creel_location == (params$gh_creel_location %||% "Grays Harbor"))

  gh_interview <- gh_interview |>
    mutate(
      event_date = as.Date(as.character(date), format = "%Y-%m-%d"),
      boat_type_clean = case_when(
        str_detect(boat_type, "(?i)commer") ~ "Commercial",
        str_detect(boat_type, "(?i)charter") ~ "Charter",
        str_detect(boat_type, "(?i)guide") ~ "Charter",
        str_detect(boat_type, "(?i)private") ~ "Private",
        TRUE ~ NA_character_
      ),
      population = case_when(
        boat_type_clean %in% c("Commercial","Charter") ~ "comm_charter",
        crabbing_mode == "Boat" & (is.na(boat_type_clean) | boat_type_clean == "Private") ~ "private_boat",
        crabbing_mode %in% c("Dock","Jetty","Beach") ~ "shore",
        TRUE ~ "shore"
      )
    )

  cat("  Grays Harbor interviews by population:\n")
  gh_interview |> count(population) |> mutate(l=sprintf("    %s: %d",population,n)) |>
    pull(l) |> walk(cat,"\n")

  # --- REVIEW ITEM 1 (2026-09-08): boat CONTACTS, before the crabbers > 0 filter -------
  # Samplers log every private boat they approach at the launch; a boat that was not
  # crabbing is recorded with crabbers = 0 (and no gear, no catch). Those rows have always
  # been dropped by the filter just below, but they are the crabbing-fraction classification
  # the boat effort model needs: per day, boats contacted and boats crabbing (a boat that
  # crabbed AND fished another fishery counts as crabbing here, because the sampler saw the
  # crab gear; OSP's crabbing-only column labels such a boat by the other fishery).
  # 2026-09-09: the workbook now carries the TRIP TYPE (crab_only / combo / other_fishery /
  # non_fishing) and the interview's creel_area, so the classification reads the trip type
  # where it exists and is restricted to the launch sites the trailer and OSP counts
  # measure (crab_fraction_contact_areas); the combo count per day is the observation of
  # the combo-trip share c. 2024-25 at the launches: 300 contacts, 143 crabbing (43 of them
  # combos), 50 non-fishing launches. Kept as dwg$boat_contacts; crab_fraction_source_rows()
  # feeds it to f and c.
  boat_contacts <- boat_contacts_from_interviews(gh_interview, params)
  # 2026-09-09: the per-contact detail (site, trip type, clock time, survey id) behind the
  # daily table, for the sampler-shift coverage diagnostic (sampler_shifts.R).
  boat_contacts_detail <- boat_contact_detail(gh_interview)
  cat(sprintf(paste0("  Private-boat contacts for the crabbing fraction (%s): %d boats on %d days, %d crabbing",
                     " (share %.3f); trip type on %d, of which %d crab-only, %d combo, %d other fishery, %d non-fishing\n"),
              attr(boat_contacts, "areas_label") %||% "all areas",
              sum(boat_contacts$boats_total), nrow(boat_contacts), sum(boat_contacts$boats_crabbing),
              if (sum(boat_contacts$boats_total) > 0) sum(boat_contacts$boats_crabbing) / sum(boat_contacts$boats_total) else NA_real_,
              sum(boat_contacts$boats_typed), sum(boat_contacts$boats_crab_only), sum(boat_contacts$boats_combo),
              sum(boat_contacts$boats_other), sum(boat_contacts$boats_nonfishing)))

  gh_interview <- gh_interview |>
    filter(!is.na(crabbers), as.numeric(crabbers) > 0) |>
    mutate(
      interview_id = paste0(survey_id,"_",interview_num),
      section_num = 1,
      angler_count = as.integer(crabbers),
      number_of_gear = as.numeric(number_of_gear),
      gear_tampered_num = suppressWarnings(as.numeric(gear_tampered)),
      hours_fished = as.numeric(hours_fished),
      crabber_hours_calc = as.numeric(crabber_hours),
      fishing_time_total = case_when(
        !is.na(crabber_hours_calc) & crabber_hours_calc > 0 ~ crabber_hours_calc,
        !is.na(hours_fished) & !is.na(angler_count) ~ hours_fished * angler_count,
        TRUE ~ NA_real_
      ),
      gear_time_total = case_when(
        !is.na(as.numeric(gear_hours)) & as.numeric(gear_hours) > 0 ~ as.numeric(gear_hours),
        !is.na(hours_fished) & !is.na(number_of_gear) ~ hours_fished * number_of_gear,
        TRUE ~ NA_real_
      ),
      dungeness_kept = replace_na(as.numeric(dungeness_kept), 0),
      red_rock_kept = replace_na(as.numeric(red_rock_kept), 0),
      trip_status = case_when(
        completed_trip=="1"~"Complete", completed_trip=="0"~"Incomplete", TRUE~NA_character_
      ),
      angler_final_int = case_when(population=="shore"~1L, population=="private_boat"~2L, TRUE~NA_integer_)
    ) |>
    # --- Exclude non-crabbing and gear-tampered interviews --------------------
    # number_of_gear == 0 means no crab gear was deployed, so the interview is
    # non-crabbing regardless of trip-completion status (complete, incomplete, or
    # blank completed_trip); drop all such rows. gear_tampered == 1 flags
    # interviews where the crabber believes a third party pulled their pots, so
    # catch and hours_fished are unreliable; drop those too. NA number_of_gear
    # (gear count not recorded) is left untouched.
    filter(is.na(number_of_gear) | number_of_gear != 0) |>
    filter(is.na(gear_tampered_num) | gear_tampered_num != 1)

  # --- Fishing-time filters, UNIT-AWARE (review item 8, 2026-09-08) ----------------
  # See apply_fishing_time_filters() below (a pure helper so the harness can test it).
  gh_interview <- apply_fishing_time_filters(gh_interview, params, req_boat_gear_time)

  # --- Duplicate interview_id repair (2026-09-28, B44) -----------------------------------
  # interview_id is survey_id + "_" + interview_num, and the catch table below is keyed on
  # it. Distinct interviews can share one (a blank interview_num: 2022-23 has survey days
  # where every interview reads "S1_NA"; or a number re-used within a survey), and then the
  # catch join gives EVERY row sharing the id the SUM of their catches, while the BSS preps'
  # distinct(interview_id) keeps only the first row. The repair makes the id unique per
  # row, before the catch table is built, and says which ids it touched.
  gh_interview <- repair_interview_ids(gh_interview)

  # --- SHORE EFFORT: Pair Float 20 + Float 17-21 gear counts ---
  dock_f20 <- params$shore_dock_float20 %||% "Westport Docks Float 20"
  dock_f17 <- params$shore_dock_float17 %||% "Westport Docks Float 17-21"
  dock_effort <- gh_effort |>
    filter(creel_area %in% c(dock_f20, dock_f17)) |>
    mutate(event_date = date,
           count_time_posix = as.POSIXct(paste(date,count_time), format="%Y-%m-%d %H:%M:%S", tz="America/Los_Angeles")) |>
    filter(!is.na(count_time_posix))

  f20 <- dock_effort |> filter(creel_area==dock_f20) |>
    arrange(event_date,count_time_posix) |> group_by(event_date) |>
    mutate(count_sequence=row_number()) |> ungroup()

  f17 <- dock_effort |> filter(creel_area==dock_f17)

  # 2026-09-28 (B45): a Float 20 count with no Float 17-21 count beside it is UNSAMPLED at
  # Float 17-21, not empty (see shore_dock_counts() below), and is filled by ratio.
  shore_effort <- shore_dock_counts(f20, f17, params) |>
    mutate(section_num=1, count_type="Gear Count", population="shore")

  # --- BOAT EFFORT ---
  boat_effort <- gh_effort |>
    filter(creel_area %in% (params$boat_launch_areas %||% c("Westport Boat Launch","Ocean Shores Boat Launch"))) |>
    mutate(event_date = date,
           count_time_posix = as.POSIXct(paste(date,count_time), format="%Y-%m-%d %H:%M:%S", tz="America/Los_Angeles"),
           count_quantity = as.numeric(boat_trailer_count)) |>
    filter(!is.na(count_time_posix), !is.na(count_quantity), count_quantity >= 0) |>
    arrange(event_date,count_time_posix) |> group_by(event_date) |>
    mutate(count_sequence=row_number(),
           count_hour = as.numeric(format(count_time_posix, "%H")) +
                        as.numeric(format(count_time_posix, "%M")) / 60) |> ungroup() |>
    select(event_date,count_sequence,count_quantity,count_hour) |>
    mutate(section_num=1, count_type="Trailer Count", population="private_boat")

  # --- COMMERCIAL TALLY ---
  comm_tally <- read_input_workbook(params$tally_file %||% "wes_commercial_tally.xlsx", sheet = in_sheet) |>
    mutate(date = as.Date(date))

  # 2026-09-10: the charter trip roster (charter_trips.xlsx), the trip-level census frame for
  # the charter part of the census; NULL when the file is absent (estimate_comm_charter falls
  # back to the tally frame and says so).
  charter_roster <- fetch_charter_roster(params)

  cat(sprintf("\n  Shore effort obs: %d (%d days)\n", nrow(shore_effort), n_distinct(shore_effort$event_date)))
  cat(sprintf("  Boat effort obs: %d (%d days)\n", nrow(boat_effort), n_distinct(boat_effort$event_date)))
  cat(sprintf("  Commercial tally days: %d\n", nrow(comm_tally)))

  # Dungeness only. The Red Rock catch rows (estimate_red_rock) were removed 2026-09-28
  # (B44). red_rock_kept stays on the interview frame: the unfished-zero-catch guard in
  # apply_fishing_time_filters() reads it as evidence that the gear was fished.
  catch <- gh_interview |> filter(dungeness_kept>0) |>
    transmute(interview_id,event_date,population,species="Dungeness",fate="Kept",
              fish_count=as.integer(dungeness_kept),catch_group="Dungeness_Kept")

  # 2026-09-09: say what the season selection captured, and stop loudly when it captured
  # nothing. See 03_R_functions/validate_season_window.R for why this exists.
  validate_season_window(effort_raw, gh_interview, params)

  return(list(
    shore_effort = shore_effort,
    boat_effort = boat_effort,
    interview = gh_interview,
    boat_contacts = boat_contacts,   # review item 1: per-day crabbing classification of contacted boats
    boat_contacts_detail = boat_contacts_detail,   # 2026-09-09: one row per contacted private boat
    catch = catch,
    comm_tally = comm_tally,
    charter_roster = charter_roster,   # 2026-09-10: the charter trip roster (NULL when absent)
    ll = tibble(centroid_lat=46.904, centroid_lon=-124.105)
  ))
}


###############################################################################
# fetch_charter_roster()  (2026-09-10)
#
# Reads 04_input_files/<params$charter_trips_file> (default charter_trips.xlsx, sheet
# "data"; built by 04_input_files/build_charter_trips.R from the season workbooks' "charter
# trips" sheets): one row per charter crab trip the operators reported, with status
# interviewed / missed / canceled. Filtered to params$charter_roster_port (default
# "Westport") and params$season_filter. Returns a tibble(season, port, vessel, date (Date),
# status, contact, notes) or NULL when the file is absent, so a workbook set without a
# roster runs exactly as before (the census uses the tally frame).
###############################################################################
fetch_charter_roster <- function(params, quiet = FALSE) {
  f <- here("04_input_files", params$charter_trips_file %||% "charter_trips.xlsx")
  if (!file.exists(f)) {
    if (!isTRUE(quiet)) cat("  Charter roster: no charter_trips.xlsx; the census uses the tally frame for charters.\n")
    return(NULL)
  }
  r <- read_input_workbook(f, sheet = params$charter_trips_sheet %||% "data")
  need <- c("season", "port", "vessel", "date", "status")
  if (!all(need %in% names(r))) stop("charter_trips.xlsx lacks column(s): ", paste(setdiff(need, names(r)), collapse = ", "), call. = FALSE)
  port <- params$charter_roster_port %||% "Westport"
  r <- r |>
    mutate(date = as.Date(as.character(date)), status = tolower(as.character(status))) |>
    filter(tolower(port) == tolower(!!port))
  if (!is.null(params$season_filter)) r <- r |> filter(season %in% params$season_filter)
  if (!isTRUE(quiet))
    cat(sprintf("  Charter roster (%s): %d trips on %d days; %d interviewed, %d missed, %d canceled\n", port, nrow(r),
                n_distinct(r$date), sum(r$status == "interviewed"), sum(r$status == "missed"), sum(r$status == "canceled")))
  r
}


###############################################################################
# apply_fishing_time_filters()  (review item 8, 2026-09-08)
#
# The historical reader dropped every interview with fishing_time_total below
# params$min_fishing_time (0.5 crabber-hours) and every private-boat interview with no
# positive gear_time_total. Both rules existed for the TIME-DENOMINATED effort units,
# where hours are the CPUE denominator and a near-zero denominator is a division hazard.
# Under the gear-DEPLOYMENT unit (production for shore since v7.7 and for the boat since
# v7.6) hours play no part in the likelihood, so the threshold only discarded data: on
# 2024-25 it removed 153 shore and 44 boat rows, most of them incomplete trips the
# incomplete-trip filter removes anyway, but also 18 shore and 14 boat complete or
# unlabelled trips of 0 to 0.4 h.
#
# The rule is now:
#   time-denominated unit (shore on crabber-hours / gear-hours) -> the historical
#       threshold, unchanged.
#   gear-deployment unit (boat always; shore in production; comm/charter) -> no time
#       threshold. ONE guard remains, params$drop_unfished_zero_catch (default TRUE): a
#       row with NO positive time in any of hours_fished / crabber_hours / gear_hours AND
#       zero catch is gear that was set and not yet fished (a boat "complete for today"
#       with its pots still soaking), not a fished deployment, whatever completed_trip
#       says. Rows with any positive time, however short, are real trips and stay.
# boat_require_gear_time (the old private-boat gear-time rule) is INERT under the
# deployment unit and is only reported; it fired on zero rows once the time threshold
# had run, and the boat has been on deployments since v7.6.
# The counts each rule removes are printed per population so a run log shows them.
#
# Pure: takes the classified interview frame (needs population, hours_fished,
# crabber_hours_calc, gear_hours, fishing_time_total, dungeness_kept, red_rock_kept)
# and returns it filtered. Effect on 2024-25 (measured 2026-09-08): shore 3,597 -> 3,739
# rows (complete-trip CPUE 0.979 -> 0.976), boat 162 -> 184 (3.26 -> 3.17),
# commercial/charter unchanged.
###############################################################################
###############################################################################
# shore_dock_counts()  (2026-09-28, B45)
#
# The shore gear count is Float 20 + Float 17-21, one row per Float 20 count. Float 17-21
# is counted only when a second sampler is on shift and has time (Matt, 2026-09-28: with
# one sampler, Float 20 is sampled and Float 17-21 is not; with two, one is posted at Float
# 20 and the other at the boat launch samples Float 17-21 when they can). So a Float 20
# count with no Float 17-21 count beside it is UNSAMPLED at Float 17-21, not empty.
#
# Until B45 every such count was given Float 17-21 = 0. On 2024-25 that is 70 of 197 shore
# count days where Float 17-21 was counted fewer times than Float 20, and 29 days where it
# was not counted at all (28 of them one-sampler days by sampler_shifts.xlsx); the mean
# daily shore gear count read about 9% low.
#
#   shore_f17_fill = "ratio" (shipped): each Float 20 count keeps the Float 17-21 count
#     paired to it (nearest in time, as before; two Float 17-21 counts paired to one Float
#     20 count are averaged, where they used to duplicate the row). An unpaired count gets
#     round(R_month x its own Float 20 count), R_month the ratio of sums of the time-paired
#     counts in its calendar month ("%Y-%m"), falling back to the ratio over every pair in
#     the data when the month has fewer than shore_f17_ratio_min_pairs pairs. Scaling the
#     concurrent Float 20 count carries the within-day pattern (Float 17-21's within-day CV
#     0.29 against Float 20's 0.34 on the 2024-25 days with repeat counts), which carrying
#     one Float 17-21 value across the day would not. Rounded, so the count stays an integer
#     for the NB2 likelihood.
#   shore_f17_fill = "zero": the pre-B45 behaviour, kept so the change can be priced.
#
# A Float 17-21 count on a day with no Float 20 count has no row to join to and is not
# used (Float 20 is the anchor, as it always was). Returns the shore effort rows with
# f20_gear, f17_gear, count_quantity, count_hour and f17_source ("observed", "ratio
# (month)", "ratio (all pairs)", "zero"); attr(, "f17_fill") carries the per-month ratios
# and counts, and a summary line is printed unless quiet = TRUE. Pure given its inputs.
###############################################################################
shore_dock_counts <- function(f20, f17, params = list(), quiet = FALSE) {
  fill <- params$shore_f17_fill %||% "ratio"
  if (!fill %in% c("ratio", "zero"))
    stop(sprintf("shore_dock_counts(): shore_f17_fill = '%s' is not one of ratio | zero.", fill), call. = FALSE)
  min_pairs <- as.numeric(params$shore_f17_ratio_min_pairs %||% 5)

  base <- f20 |>
    # 2026-09-08 (review item 2): keep the count's clock time (decimal hour) so the shore
    # turnover can be evaluated at the hours the counts were actually taken.
    mutate(count_hour = as.numeric(format(count_time_posix, "%H")) +
                        as.numeric(format(count_time_posix, "%M")) / 60) |>
    select(event_date, count_sequence, f20_gear = total_gear_count, count_hour)

  f17_paired <- if (nrow(f17) > 0 && nrow(f20) > 0) {
    f17 |>
      left_join(f20 |> select(event_date, count_sequence, f20_time = count_time_posix),
                by = "event_date", relationship = "many-to-many") |>
      filter(!is.na(f20_time)) |>
      mutate(time_diff = abs(as.numeric(difftime(count_time_posix, f20_time, units = "mins")))) |>
      group_by(event_date, survey_id, count_time) |> slice_min(time_diff, n = 1, with_ties = FALSE) |> ungroup() |>
      group_by(event_date, count_sequence) |>
      summarise(f17_gear = mean(as.numeric(total_gear_count)), n_f17 = n(), .groups = "drop")
  } else tibble(event_date = as.Date(character()), count_sequence = integer(), f17_gear = numeric(), n_f17 = integer())

  out <- base |>
    left_join(f17_paired, by = c("event_date", "count_sequence")) |>
    mutate(.observed = !is.na(f17_gear), .month = format(event_date, "%Y-%m"))

  pairs <- out |> filter(.observed)
  r_all <- if (nrow(pairs) && sum(pairs$f20_gear, na.rm = TRUE) > 0)
             sum(pairs$f17_gear, na.rm = TRUE) / sum(pairs$f20_gear, na.rm = TRUE) else NA_real_
  by_month <- pairs |> group_by(.month) |>
    summarise(n_pairs = n(), f17_sum = sum(f17_gear, na.rm = TRUE), f20_sum = sum(f20_gear, na.rm = TRUE), .groups = "drop") |>
    mutate(r_month = if_else(f20_sum > 0, f17_sum / f20_sum, NA_real_))

  if (identical(fill, "ratio") && is.finite(r_all)) {
    out <- out |>
      left_join(by_month |> select(.month, n_pairs, r_month), by = ".month") |>
      mutate(.use_month = is.finite(r_month) & !is.na(n_pairs) & n_pairs >= min_pairs,
             .r = if_else(.use_month, r_month, r_all),
             f17_source = case_when(.observed ~ "observed",
                                    .use_month ~ "ratio (month)",
                                    TRUE ~ "ratio (all pairs)"),
             f17_gear = if_else(.observed, round(f17_gear), round(.r * f20_gear)))
  } else {
    if (identical(fill, "ratio") && !quiet)
      cat("  *** Float 17-21: no Float 17-21 count is paired with a Float 20 count anywhere in the data, so there is no ratio to fill from; unpaired counts are 0 (shore_f17_fill = 'ratio' could not apply). ***\n")
    out <- out |> mutate(f17_source = if_else(.observed, "observed", "zero"),
                         f17_gear = if_else(.observed, round(f17_gear), 0))
  }
  out <- out |>
    mutate(count_quantity = f20_gear + f17_gear) |>
    select(event_date, count_sequence, f20_gear, f17_gear, count_quantity, count_hour, f17_source)

  attr(out, "f17_fill") <- list(fill = fill, r_all = r_all, by_month = by_month, min_pairs = min_pairs,
                                n_counts = nrow(out), n_observed = sum(out$f17_source == "observed"),
                                n_filled = sum(grepl("^ratio", out$f17_source)),
                                filled_gear = sum(out$f17_gear[grepl("^ratio", out$f17_source)]))
  if (!quiet) {
    k <- attr(out, "f17_fill")
    cat(sprintf(paste0("  Float 17-21 (shore_f17_fill = '%s'): %d of %d Float 20 counts carry an observed Float 17-21 count; ",
                       "%d unpaired filled by ratio (all-pairs ratio %.3f; month ratios %s), %s gear.\n"),
                fill, k$n_observed, k$n_counts, k$n_filled, if (is.finite(r_all)) r_all else NA_real_,
                if (nrow(by_month)) paste(sprintf("%.2f", range(by_month$r_month, na.rm = TRUE)), collapse = "-") else "none",
                format(round(k$filled_gear), big.mark = ",")))
  }
  out
}

###############################################################################
# repair_interview_ids()  (2026-09-28, B44)
#
# Makes df$interview_id unique per row. Rows that share an id keep it with a suffix
# "__r<k>" (k = 1, 2, ... in row order within the id); rows with a unique id are
# untouched, so a dataset without the defect passes through byte-identical. The repaired
# ids are recorded in attr(df, "interview_id_repairs") (one row per original id: the id,
# how many rows shared it, and their seasons) and printed unless quiet = TRUE.
#
# Why a suffix and not a de-duplication: every shared id found in 2022-23 to 2025-26 (14
# ids, 37 Grays Harbor rows) is two or more DIFFERENT interviews (different party size,
# gear, catch or time), never one interview entered twice, so dropping rows would lose
# real interviews. Pure; base R only.
###############################################################################
repair_interview_ids <- function(df, quiet = FALSE) {
  if (!is.data.frame(df) || !nrow(df) || !"interview_id" %in% names(df)) return(df)
  id  <- as.character(df$interview_id)
  dup <- id %in% id[duplicated(id)]
  if (!any(dup)) return(df)
  k <- stats::ave(seq_along(id), id, FUN = seq_along)
  seasons <- if ("season" %in% names(df)) as.character(df$season) else rep(NA_character_, length(id))
  shared  <- unique(id[dup])
  log <- data.frame(
    interview_id = shared,
    n_rows       = vapply(shared, function(x) sum(id == x), integer(1)),
    seasons      = vapply(shared, function(x) paste(unique(seasons[id == x]), collapse = ","), character(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  id[dup] <- paste0(id[dup], "__r", k[dup])
  df$interview_id <- id
  attr(df, "interview_id_repairs") <- log
  if (!isTRUE(quiet))
    cat(sprintf(paste0("  Interview ids: %d id(s) were shared by %d rows (distinct interviews under one ",
                       "survey_id + interview_num); each row now carries its own id, so its own catch.\n",
                       "    %s\n"),
                nrow(log), sum(log$n_rows),
                paste(sprintf("%s x%d (%s)", log$interview_id, log$n_rows, log$seasons), collapse = "; ")))
  df
}

apply_fishing_time_filters <- function(df, params, req_boat_gear_time = TRUE, quiet = FALSE) {
  .shore_unit <- params$shore_effort_unit %||% "crabber-hours"
  .time_units <- c("crabber-hours", "gear-hours")
  .min_time   <- as.numeric(params$min_fishing_time %||% 0.5)
  .zero_guard <- isTRUE(params$drop_unfished_zero_catch %||% TRUE)
  gh <- as.numeric(suppressWarnings(as.numeric(df$gear_hours %||% rep(NA_real_, nrow(df)))))
  df <- df |>
    mutate(
      .time_unit_pop = (population == "shore" & .shore_unit %in% .time_units),
      .any_pos_time  = (!is.na(hours_fished) & hours_fished > 0) |
                       (!is.na(crabber_hours_calc) & crabber_hours_calc > 0) |
                       (!is.na(gh) & gh > 0),
      .zero_catch    = (dungeness_kept <= 0) & (red_rock_kept <= 0),
      .drop_time     = .time_unit_pop & (is.na(fishing_time_total) | fishing_time_total < .min_time),
      .drop_unfished = !.time_unit_pop & .zero_guard & !.any_pos_time & .zero_catch
    )
  if (!isTRUE(quiet)) {
    .tab <- df |>
      group_by(population) |>
      summarise(n = n(), time_threshold = sum(.drop_time),
                unfished_zero_catch = sum(.drop_unfished), .groups = "drop")
    cat(sprintf("  Fishing-time filters (unit-aware; shore unit '%s'; boat_require_gear_time %s is inert under deployments):\n",
                .shore_unit, req_boat_gear_time))
    for (i in seq_len(nrow(.tab)))
      cat(sprintf("    %-13s %4d rows | dropped: time threshold %d, unfished zero-catch %d [%s]\n",
                  .tab$population[i], .tab$n[i], .tab$time_threshold[i], .tab$unfished_zero_catch[i],
                  if (.tab$population[i] == "shore" && .shore_unit %in% .time_units) "time unit" else "deployment unit"))
  }
  df |>
    filter(!.drop_time, !.drop_unfished) |>
    select(-.time_unit_pop, -.any_pos_time, -.zero_catch, -.drop_time, -.drop_unfished)
}


###############################################################################
# boat_contacts_from_interviews()  (review item 1, 2026-09-08; trip types 2026-09-09)
#
# Per-day classification of contacted private boats, from the CLASSIFIED but UNFILTERED
# interview frame (population assigned, crabbers still raw). Returns
# tibble(event_date, boats_total, boats_crabbing, boats_typed, boats_crab_only, boats_combo,
#        boats_other, boats_nonfishing), the shape crab_fraction.R consumes.
#
# CLASSIFICATION RULE, per row:
#   trip_type_class present (2026-09-09 workbook): crabbing = crab_only | combo;
#       not crabbing = other_fishery | non_fishing (a launch with no fishing is still a
#       boat in the trailer and OSP totals, so it is a not-crabbing boat to f).
#   trip_type_class absent (older rows): crabbing = crabbers > 0, combo unknown.
#   A row with neither a trip type nor a crabbers value is unclassified and counts nowhere.
# boats_typed / boats_crab_only / boats_combo count TYPED rows only; they are the data for
# the combo-trip share c (combo of crabbing, per day).
#
# SITE RESTRICTION. The trailer count and the OSP total measure boats launched at the
# ramp; a private boat interviewed at the marina or the docks is moored, not trailered,
# and belongs to neither count. params$crab_fraction_contact_areas (default: the
# boat_launch_areas, "all" for no restriction) selects the interviews by creel_area when
# the workbook carries that column. 2024-25: 300 contacts at the launches against 422
# everywhere (34 marina rows, all crab-only; 88 dock rows).
#
# Pure, so the harness tests it on a synthetic frame.
###############################################################################
boat_contacts_from_interviews <- function(gh_interview, params = list()) {
  empty <- tibble(event_date = as.Date(character()), boats_total = integer(), boats_crabbing = integer(),
                  boats_typed = integer(), boats_crab_only = integer(), boats_combo = integer(),
                  boats_other = integer(), boats_nonfishing = integer())
  if (is.null(gh_interview) || !nrow(gh_interview) ||
      !all(c("population", "event_date") %in% names(gh_interview)) ||
      !any(c("crabbers", "trip_type_class") %in% names(gh_interview))) return(empty)
  d <- gh_interview |> filter(population == "private_boat", !is.na(event_date))
  areas <- params$crab_fraction_contact_areas %||% params$boat_launch_areas %||%
    c("Westport Boat Launch", "Ocean Shores Boat Launch")
  areas_label <- "all areas"
  if ("creel_area" %in% names(d) && !identical(tolower(areas[1]), "all")) {
    d <- d |> filter(creel_area %in% areas)
    areas_label <- paste(areas, collapse = " + ")
  }
  cr <- if ("crabbers" %in% names(d)) suppressWarnings(as.numeric(d$crabbers)) else rep(NA_real_, nrow(d))
  tt <- if ("trip_type_class" %in% names(d)) as.character(d$trip_type_class) else rep(NA_character_, nrow(d))
  tt[!tt %in% c("crab_only", "combo", "other_fishery", "non_fishing")] <- NA_character_
  d$.typed <- !is.na(tt)
  d$.crab  <- ifelse(d$.typed, tt %in% c("crab_only", "combo"), !is.na(cr) & cr > 0)
  d$.class <- ifelse(d$.typed | !is.na(cr), TRUE, FALSE)      # classifiable at all
  d$.tt <- tt
  out <- d |>
    filter(.class) |>
    group_by(event_date) |>
    summarise(boats_total = n(), boats_crabbing = sum(.crab),
              boats_typed = sum(.typed),
              boats_crab_only = sum(.tt %in% "crab_only"), boats_combo = sum(.tt %in% "combo"),
              boats_other = sum(.tt %in% "other_fishery"), boats_nonfishing = sum(.tt %in% "non_fishing"),
              .groups = "drop") |>
    mutate(event_date = as.Date(event_date), across(-event_date, as.integer)) |>
    arrange(event_date)
  attr(out, "areas_label") <- areas_label
  out
}


###############################################################################
# boat_contact_detail()  (2026-09-09)
#
# One row per contacted PRIVATE boat (every area, before any filter): the date, the
# creel_area, the survey id (the link to sampler_shifts.xlsx), the trip-type class, the
# crabbers value, and the contact clock time as a decimal hour (from interview_time,
# "HH:MM"). Feeds the shift-coverage diagnostic; NULL columns degrade to NA.
###############################################################################
boat_contact_detail <- function(gh_interview) {
  empty <- tibble(event_date = as.Date(character()), creel_area = character(), survey_id = character(),
                  trip_type_class = character(), crabbers = numeric(), contact_hour = numeric())
  if (is.null(gh_interview) || !nrow(gh_interview) || !all(c("population", "event_date") %in% names(gh_interview))) return(empty)
  d <- gh_interview |> filter(population == "private_boat", !is.na(event_date))
  if (!nrow(d)) return(empty)
  .col <- function(nm, cast) if (nm %in% names(d)) cast(d[[nm]]) else rep(cast(NA), nrow(d))
  it <- .col("interview_time", as.character)
  hr <- suppressWarnings(as.numeric(sub(":.*$", "", it)) + as.numeric(sub("^[^:]*:", "", it)) / 60)
  tibble(event_date = as.Date(d$event_date),
         creel_area = .col("creel_area", as.character),
         survey_id  = .col("survey_id", as.character),
         trip_type_class = .col("trip_type_class", as.character),
         crabbers   = suppressWarnings(.col("crabbers", as.numeric)),
         contact_hour = ifelse(is.finite(hr) & hr >= 0 & hr < 24, hr, NA_real_))
}
