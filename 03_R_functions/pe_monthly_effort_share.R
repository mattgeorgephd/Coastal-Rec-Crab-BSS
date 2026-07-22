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
# pe_monthly_effort_share.R  (shared by the pooled and gear-resolved drivers)
#
# Monthly effort SHARE for distributing a PE-fallback component's catch/effort across
# months (the piece the drivers' "7.8" / "7.8b" monthly-estimate blocks duplicated; T4.4).
# Each caller keeps its OWN draw accumulation and uncertainty handling, which legitimately
# differ (7.8 pools modes with point shares; 7.8b is per-mode with a lognormal effort-level
# scale; the gear block pools modes), so ONLY the share math is centralized here.
#
# The daily effort matches run_pe_*(): boat on the gear-DEPLOYMENT scale
# (mean_count * gear_per_group * tau_boat, day-length-free), shore on crabber-hours
# (mean_count * crabbers_per_gear * day_length). gear_per_group is recomputed from the
# same sub-season interview filter run_pe uses, so the split stays consistent with the PE.
# The share is normalized, so it changes only the across-month distribution, never the
# component total.
#
# na.rm defaults TRUE: a day in effort_index with no day_length/month_label match is
# dropped rather than nulling the whole month. This standardizes the one na.rm difference
# the two pooled blocks carried (behavior-neutral when every effort day is in-window), and
# for the gear driver it also gives the boat PE-fallback its correct deployment-scale
# effort instead of the shore formula that block used inline.
###############################################################################

pe_monthly_effort_share <- function(pop, summ_ss, days_ss, params, na.rm = TRUE) {
  is_boat_pe <- stringr::str_detect(pop, "private_boat")
  gpg_pe <- params$gear_per_group_default %||% 4.0
  tau_pe <- params$tau_boat_prior_mu %||% 1.2
  if (is_boat_pe) {
    rd_pe <- summ_ss$interview |>
      dplyr::filter(!is.na(number_of_gear), number_of_gear > 0, angler_count > 0)
    if (nrow(rd_pe) > 0) gpg_pe <- mean(rd_pe$number_of_gear)
  }
  pe_daily <- summ_ss$effort_index |>
    dplyr::filter(count_sequence <= params$bss_max_count_seq) |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(mean_count = mean(count_quantity), .groups = "drop") |>
    dplyr::left_join(dplyr::select(days_ss, event_date, day_length, month_label),
                     by = "event_date") |>
    dplyr::mutate(daily_effort = if (is_boat_pe) mean_count * gpg_pe * tau_pe
                                 else mean_count * summ_ss$crabbers_per_gear * day_length)
  pe_daily |>
    dplyr::group_by(month_label) |>
    dplyr::summarise(month_effort = sum(daily_effort, na.rm = na.rm), .groups = "drop") |>
    dplyr::mutate(share = month_effort / sum(month_effort, na.rm = na.rm))
}
