#!/usr/bin/env Rscript
###############################################################################
# gear_coverage_audit.R  --  GR-7 Phase 0: single-gear interview coverage audit.
#
# Answers the question that GATES the per-gear CPUE build (GR-7 / Option A1):
# for each population x sub-season, how many interviews report a SINGLE gear
# type (the only ones that, per the design decision, may feed a gear-specific
# CPUE), how many report MIXED (multiple) gear, and which gears therefore clear
# the estimability threshold (bss_min_gear_effective_n, default 15).
#
# This is a design-time audit, NOT part of a model run. It re-derives the gear
# classification and the standard interview filters directly from the input
# workbook, so it can be run on its own in seconds:
#
#     source("06_diagnostics/gear_coverage_audit.R")
#
# It prints a per-cell table and writes gear_coverage_audit.csv to the repo root.
# Output columns per (population, subseason, gear): single_gear_n (interviews
# reporting only that gear), frac_effective_n (the OLD fractional metric the
# collapse used), estimable (single_gear_n >= threshold). A "Mixed" row per cell
# reports the multi-gear count and whether Mixed itself clears the threshold.
#
# Governing rule (Matt, 2026-07-20): only single-gear interviews contribute to a
# gear-specific CPUE; multi-gear interviews form a "Mixed" gear. See
# 07_documentation/development_notes/GR-7-per-gear-CPUE-design.md.
###############################################################################

suppressPackageStartupMessages({ library(here); library(dplyr); library(stringr); library(readxl); library(tidyr) })

if (!exists("run_config")) source(here::here("run_config.R"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
p <- run_config

thr        <- p$bss_min_gear_effective_n %||% 15
season_f   <- p$season_filter            %||% "2024-25"
loc        <- p$gh_creel_location        %||% "Grays Harbor"
min_ft     <- p$min_fishing_time         %||% 0.5
d_start    <- as.Date(p$est_date_start   %||% "2024-09-16")
d_potclose <- as.Date(p$pot_closure_end  %||% "2024-11-30")
d_potopen  <- as.Date(p$pot_open_date    %||% "2024-12-01")
d_end      <- as.Date(p$est_date_end     %||% "2025-09-15")

int <- readxl::read_excel(here::here("04_input_files", p$interview_file %||% "interview_combined.xlsx"),
                          sheet = p$input_sheet %||% "data") |>
  mutate(completed_trip = as.character(completed_trip)) |>
  filter(season == season_f, creel_location == loc) |>
  mutate(gear_type = tidyr::replace_na(as.character(gear_type), ""))

int <- int |>
  mutate(
    event_date = as.Date(as.character(date)),
    boat_type_clean = case_when(
      str_detect(boat_type, "(?i)commer")   ~ "Commercial",
      str_detect(boat_type, "(?i)charter")  ~ "Charter",
      str_detect(boat_type, "(?i)guide")    ~ "Charter",
      str_detect(boat_type, "(?i)private")  ~ "Private",
      TRUE ~ NA_character_),
    population = case_when(
      boat_type_clean %in% c("Commercial", "Charter") ~ "comm_charter",
      crabbing_mode == "Boat" & (is.na(boat_type_clean) | boat_type_clean == "Private") ~ "private_boat",
      crabbing_mode %in% c("Dock", "Jetty", "Beach") ~ "shore",
      TRUE ~ "shore"),
    crabbers_n = suppressWarnings(as.numeric(crabbers)),
    ch = suppressWarnings(as.numeric(crabber_hours)),
    hf = suppressWarnings(as.numeric(hours_fished)),
    fishing_time_total = case_when(!is.na(ch) & ch > 0 ~ ch,
                                   !is.na(hf) ~ hf * crabbers_n, TRUE ~ NA_real_),
    trip_status = case_when(completed_trip == "1" ~ "Complete",
                            completed_trip == "0" ~ "Incomplete", TRUE ~ NA_character_)) |>
  filter(!is.na(crabbers_n), crabbers_n > 0,
         !is.na(fishing_time_total), fishing_time_total >= min_ft)

if (isTRUE(p$filter_incomplete_trips %||% TRUE))
  int <- int |> filter(is.na(trip_status) | trip_status == "Complete")   # keep Complete + NA

# Gear classification (identical regex to prep_bss_crab_gear.R:197-201)
int <- int |>
  mutate(
    has_pot      = as.integer(str_detect(gear_type, "(?i)\\bpot\\b") & !str_detect(gear_type, "(?i)\\bslip\\s*ring\\b")),
    has_ring_net = as.integer(str_detect(gear_type, "(?i)\\bring\\s*net\\b")),
    has_trap     = as.integer(str_detect(gear_type, "(?i)\\b(trap|star)\\b")),
    has_snare    = as.integer(str_detect(gear_type, "(?i)\\bsnare\\b")),
    subseason = case_when(
      event_date >= d_start   & event_date <= d_potclose ~ "pot_closure",
      event_date >= d_potopen & event_date <= d_end      ~ "all_gear",
      TRUE ~ NA_character_)) |>
  filter(!is.na(subseason))

# Regulatory exclusion: pots illegal in the pot-closure sub-season
int <- int |> mutate(has_pot = if_else(subseason == "pot_closure", 0L, has_pot))

gears  <- c(Pot = "has_pot", `Ring Net` = "has_ring_net", Trap = "has_trap", Snare = "has_snare")
int <- int |> mutate(n_types = has_pot + has_ring_net + has_trap + has_snare,
                     n_types_eff = pmax(n_types, 1L))

rows <- list()
cat(sprintf("\nGR-7 Phase 0 gear-coverage audit  |  season %s  |  single-gear threshold >= %d\n",
            season_f, thr))
cat(strrep("=", 78), "\n")
for (popn in c("shore", "private_boat")) {
  for (ss in c("pot_closure", "all_gear")) {
    sub <- int |> filter(population == popn, subseason == ss)
    if (nrow(sub) == 0) next
    n_single <- sum(sub$n_types == 1); n_mixed <- sum(sub$n_types > 1)
    cat(sprintf("\n%-13s / %-11s  (N=%d;  single-gear=%d, mixed=%d, mixed share=%.0f%%)\n",
                popn, ss, nrow(sub), n_single, n_mixed,
                100 * n_mixed / max(n_single + n_mixed, 1)))
    for (gt in names(gears)) {
      hc <- gears[[gt]]
      single <- sum(sub$n_types == 1 & sub[[hc]] == 1)
      frac   <- sum(sub[[hc]] / sub$n_types_eff)
      est    <- single >= thr
      cat(sprintf("    %-10s single-gear=%4d  frac_eff_n=%6.1f  %s\n",
                  gt, single, frac, if (est) "ESTIMABLE" else "-"))
      rows[[length(rows) + 1]] <- data.frame(population = popn, subseason = ss, gear = gt,
        single_gear_n = single, frac_effective_n = round(frac, 1), estimable = est)
    }
    rows[[length(rows) + 1]] <- data.frame(population = popn, subseason = ss, gear = "Mixed",
      single_gear_n = n_mixed, frac_effective_n = NA_real_, estimable = n_mixed >= thr)
    est_set <- names(gears)[vapply(names(gears), function(gt)
      sum(sub$n_types == 1 & sub[[gears[[gt]]]] == 1) >= thr, logical(1))]
    if (n_mixed >= thr) est_set <- c(est_set, "Mixed")
    cat(sprintf("    => estimable gear set: %s\n",
                if (length(est_set)) paste(est_set, collapse = ", ") else "NONE (stays G=1 'All')"))
  }
}
out <- dplyr::bind_rows(rows)
write.csv(out, here::here("gear_coverage_audit.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s\n", here::here("gear_coverage_audit.csv")))
