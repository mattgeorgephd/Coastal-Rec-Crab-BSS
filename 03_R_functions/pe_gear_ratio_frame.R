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
# ---------------------------------------------------------------------------
# pe_gear_ratio_frame.R
#
# ONE definition of which interviews the Point Estimator uses for the boat gear-per-group
# ratio, shared by run_pe_pooled() and run_pe_gear().
#
# THE INCONSISTENCY THIS CLOSES. The BSS learns R_G_boat from Gear_A_boat ~ Poisson(R_G_boat),
# whose frame (`intA`) descends from `int_d` AFTER filter_incomplete_trips is applied, so the
# BSS behaves like the four-arm diagnostic's "exclude" arm. Both Point Estimators built the
# same ratio from the UNFILTERED interview set, i.e. the "gear_only" arm. Two arms of one
# fused estimator disagreeing about which interviews count is a documented inconsistency
# (sensitivity_incomplete_trips.csv has reported it since 2026-08-25 as
# production_arm_bss = "exclude" against production_arm_pe = "gear_only") and it is the kind
# of thing a reviewer asks about first.
#
# THE SIZE, measured 2026-09-01 before deciding: the gear ratio moves 3.552 to 3.550 on the
# boat, 0.1%, and at most 0.7% on any component (shore all-gear, where n = 2,741 makes even
# that detectable at p = 0.0006 and where the shift is DOWNWARD, incomplete trips carrying
# FEWER gear, which is what a trip interrupted while still setting up looks like). So this is
# a CONSISTENCY fix, not an accuracy one. It is worth making precisely because it is small:
# there is no estimate to defend, only an inconsistency to remove.
#
# WHY A SHARED FUNCTION RATHER THAN THE SAME EDIT TWICE. The gear-resolved PE's comment
# already claimed it matched the BSS's R_G_boat while it did not. Two copies of a rule that
# is supposed to agree is how that happens; 07_documentation/CLAUDE.md makes shared helpers
# the convention for exactly this reason.
#
# ARMS
#   "match_bss"  (default) apply filter_incomplete_trips, matching the BSS frame.
#   "gear_only"  the pre-2026-09-02 behaviour: keep an interrupted trip's fully observed
#                gear count. Kept because it is a defensible position (the gear count IS
#                complete; only the catch is truncated) and because historical runs need to
#                stay reproducible.
#
# A caller that supplies its own frame is deliberately choosing one, which is what the
# four-arm diagnostic does, so nothing is applied on top of it.
# ---------------------------------------------------------------------------
pe_gear_ratio_frame <- function(interview, interview_gear = NULL, params = list(),
                                label = "", quiet = FALSE) {
  frame <- interview_gear %||% interview
  if (is.null(frame) || !nrow(frame)) return(frame)
  arm <- params$pe_gear_ratio_arm %||% "match_bss"
  if (!arm %in% c("match_bss", "gear_only"))
    stop("params$pe_gear_ratio_arm must be \"match_bss\" or \"gear_only\"; got \"", arm, "\".",
         call. = FALSE)
  if (!identical(arm, "match_bss") || !is.null(interview_gear) ||
      !isTRUE(params$filter_incomplete_trips) || !"trip_status" %in% names(frame))
    return(frame)
  n0 <- nrow(frame)
  # Missing status is KEPT, matching prep_bss_crab_pooled.R exactly: a blank completed_trip
  # may well be a complete trip, and dropping it would lose data the BSS keeps.
  out <- frame[frame$trip_status == "Complete" | is.na(frame$trip_status), , drop = FALSE]
  if (!isTRUE(quiet))
    cat(sprintf("  PE gear ratio arm 'match_bss'%s: %d -> %d interviews (matches the BSS R_G_boat frame)\n",
                if (nzchar(label)) paste0(" [", label, "]") else "", n0, nrow(out)))
  out
}
