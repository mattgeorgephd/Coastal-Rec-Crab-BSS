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
# prep_population_summary.R
#
# Filter the data-warehouse bundle (dwg) to one population x sub-season and build
# its effort-index, interview, and catch-wide frames, plus the crabbers-per-gear
# and empirical R_G ratios. Extracted from the pooled and gear-resolved drivers so
# both share one implementation (the pooled superset, which also reports
# empirical_R_G, is used; it is inert for the gear model). Auto-sourced by both
# drivers via the 03_R_functions walk. catch_groups is derived from params.
###############################################################################

prep_population_summary <- function(dwg, population_name, date_start, date_end, params) {
  # Derive catch groups from the centralized config (single source of truth).
  catch_groups <- if (isTRUE(params$estimate_red_rock)) c("Dungeness_Kept", "Red_Rock_Kept") else "Dungeness_Kept"
  ds <- as.Date(date_start); de <- as.Date(date_end)
  summ <- list()

  if(population_name == "shore") {
    summ$effort_index <- dwg$shore_effort |> filter(between(event_date, ds, de))
  } else if(population_name == "private_boat") {
    summ$effort_index <- dwg$boat_effort |> filter(between(event_date, ds, de))
  }

  int_pop <- dwg$interview |>
    filter(population == population_name, between(event_date, ds, de))

  catch_wide <- dwg$catch |>
    filter(population == population_name) |>
    group_by(interview_id, catch_group) |>
    summarise(fish_count=sum(fish_count), .groups="drop") |>
    pivot_wider(names_from=catch_group, values_from=fish_count, values_fill=0)

  summ$interview <- int_pop |>
    left_join(catch_wide, by="interview_id") |>
    mutate(across(any_of(catch_groups), ~replace_na(.,0)))

  ratio_data <- summ$interview |> filter(!is.na(number_of_gear), number_of_gear>0, angler_count>0)
  summ$crabbers_per_gear <- if(nrow(ratio_data)>0) sum(ratio_data$angler_count)/sum(ratio_data$number_of_gear) else 1.0

  # Compute empirical gear-per-crabber ratio for R_G prior
  summ$empirical_R_G <- if(nrow(ratio_data)>0) sum(ratio_data$number_of_gear)/sum(ratio_data$angler_count) else 1.3

  cat(sprintf("\n  %s [%s to %s]: %d effort obs, %d interviews, crab/gear=%.2f, R_G_empirical=%.2f\n",
              population_name, ds, de, nrow(summ$effort_index), nrow(summ$interview),
              summ$crabbers_per_gear, summ$empirical_R_G))

  return(summ)
}
