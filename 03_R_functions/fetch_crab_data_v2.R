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
# fetch_crab_data_v2.R  (gear-resolved driver)
#
# Read and assemble the gear-resolved model's input data (effort counts, interviews
# with weighted gear-type classification, catch, commercial tally). Extracted
# verbatim from the gear-resolved driver; pure given params. Auto-sourced by both
# drivers via the 03_R_functions walk (only the gear-resolved driver calls it).
#
# INPUTS (2026-07-16): all inputs are now .xlsx workbooks with a single "data"
# sheet. Filenames, the sheet name, the interview creel location, and the effort
# site lists are run_config parameters (historical values as %||% defaults). Dates
# are stored as ISO yyyy-mm-dd text and parsed with as.Date().
###############################################################################

fetch_crab_data_v2 <- function(params) {

  cat("Reading data...\n")
  in_sheet <- params$input_sheet %||% "data"

  # --- Read effort workbook (ISO yyyy-mm-dd dates) ---
  effort_raw <- readxl::read_excel(
      here("04_input_files", params$effort_file %||% "effort_combined.xlsx"), sheet = in_sheet) |>
    filter(season == params$season_filter) |> mutate(date = as.Date(date))

  # --- Read interview workbook ---
  interview_raw <- readxl::read_excel(
      here("04_input_files", params$interview_file %||% "interview_combined.xlsx"), sheet = in_sheet) |>
    mutate(completed_trip = as.character(completed_trip)) |>
    filter(season == params$season_filter)

  # --- Filter to Grays Harbor sites ---
  gh_effort <- effort_raw |> filter(creel_area %in% (params$gh_effort_areas %||% c(
    "Westport Docks Float 20","Westport Docks Float 17-21",
    "Westport Boat Launch","Westport Marina","Westport Jetty",
    "Ocean Shores Boat Launch","Damon Point")))

  gh_interview <- interview_raw |> filter(creel_location == (params$gh_creel_location %||% "Grays Harbor"))

  # --- Classify each interview into a population ---
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

  # --- Process interviews: compute fishing time, catch, trip status ---
  gh_interview <- gh_interview |>
    filter(!is.na(crabbers), as.numeric(crabbers) > 0) |>
    mutate(
      interview_id = paste0(survey_id,"_",interview_num),
      section_num = 1,
      angler_count = as.integer(crabbers),
      number_of_gear = as.numeric(number_of_gear),
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
    filter(!is.na(fishing_time_total), fishing_time_total >= params$min_fishing_time)

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

  if(nrow(f17) > 0) {
    f17_paired <- f17 |>
      left_join(f20 |> select(event_date,count_sequence,f20_time=count_time_posix),
                by="event_date", relationship="many-to-many") |>
      mutate(time_diff=abs(as.numeric(difftime(count_time_posix,f20_time,units="mins")))) |>
      group_by(event_date,survey_id,count_time) |> slice_min(time_diff,n=1,with_ties=FALSE) |> ungroup() |>
      select(event_date,count_sequence,f17_gear=total_gear_count)
  } else {
    f17_paired <- tibble(event_date=Date(),count_sequence=integer(),f17_gear=numeric())
  }

  shore_effort <- f20 |> select(event_date,count_sequence,f20_gear=total_gear_count) |>
    left_join(f17_paired, by=c("event_date","count_sequence")) |>
    mutate(f17_gear=replace_na(f17_gear,0), count_quantity=f20_gear+f17_gear,
           section_num=1, count_type="Gear Count", population="shore")

  # --- BOAT EFFORT: Trailer counts at boat launches ---
  boat_effort <- gh_effort |>
    filter(creel_area %in% (params$boat_launch_areas %||% c("Westport Boat Launch","Ocean Shores Boat Launch"))) |>
    mutate(event_date = date,
           count_time_posix = as.POSIXct(paste(date,count_time), format="%Y-%m-%d %H:%M:%S", tz="America/Los_Angeles"),
           count_quantity = as.numeric(boat_trailer_count)) |>
    filter(!is.na(count_time_posix), !is.na(count_quantity), count_quantity >= 0) |>
    arrange(event_date,count_time_posix) |> group_by(event_date) |>
    mutate(count_sequence=row_number()) |> ungroup() |>
    select(event_date,count_sequence,count_quantity) |>
    mutate(section_num=1, count_type="Trailer Count", population="private_boat")

  # --- COMMERCIAL TALLY ---
  comm_tally <- readxl::read_excel(
      here("04_input_files", params$tally_file %||% "wes_commercial_tally.xlsx"), sheet = in_sheet) |>
    mutate(date = as.Date(date))

  cat(sprintf("\n  Shore effort obs: %d (%d days)\n", nrow(shore_effort), n_distinct(shore_effort$event_date)))
  cat(sprintf("  Boat effort obs: %d (%d days)\n", nrow(boat_effort), n_distinct(boat_effort$event_date)))
  cat(sprintf("  Commercial tally days: %d\n", nrow(comm_tally)))

  # --- CATCH tables ---
  catch <- bind_rows(
    gh_interview |> filter(dungeness_kept>0) |>
      transmute(interview_id,event_date,population,species="Dungeness",fate="Kept",
                fish_count=as.integer(dungeness_kept),catch_group="Dungeness_Kept"),
    if(params$estimate_red_rock) {
      gh_interview |> filter(red_rock_kept>0) |>
        transmute(interview_id,event_date,population,species="Red_Rock",fate="Kept",
                  fish_count=as.integer(red_rock_kept),catch_group="Red_Rock_Kept")
    } else { tibble() }
  )

  return(list(
    shore_effort = shore_effort,
    boat_effort = boat_effort,
    interview = gh_interview,
    catch = catch,
    comm_tally = comm_tally,
    ll = tibble(centroid_lat=46.904, centroid_lon=-124.105)
  ))
}
