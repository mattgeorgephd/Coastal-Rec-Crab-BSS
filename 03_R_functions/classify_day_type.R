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
# classify_day_type.R
#
# Standalone day-type classifier: assign holiday / weekend / weekday to ANY date,
# independent of the estimation calendar (used for diagnostic plots that may show
# data outside the estimation window). Extracted from the pooled driver. weekends
# and holidays are read from params (params$days_wkend, params$crabbing_holiday_dates).
###############################################################################

classify_day_type <- function(dates, params) {
  weekends <- params$days_wkend
  holidays <- params$crabbing_holiday_dates
  case_when(
    dates %in% holidays ~ "holiday",
    weekdays(dates) %in% weekends ~ "weekend",
    TRUE ~ "weekday"
  )
}
