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
# build_subseasons.R  (shared by the pooled and gear-resolved drivers)
#
# Build the list of within-season sub-seasons, split by the pot-closure window.
#
# The "pot closure" is the period when pots are not legal, so only non-pot gear
# (ring nets, snares, traps) may be used. It is NOT ring-net-only, so the internal
# label "ring_net_only" is kept only for output-filename continuity; the display_name
# is "Pot closure" and user-facing strings should use it. Outside the closure window,
# pots are allowed (all-gear).
#
# The closure window is given explicitly by params$pot_closure_start /
# params$pot_closure_end (added 2026-07-13). It is no longer assumed to start at the
# season start. This builder therefore handles the general case:
#   - an optional pre-closure all-gear period  [est_date_start, pot_closure_start - 1]
#   - the pot closure                          [pot_closure_start, pot_closure_end]
#   - an optional post-closure all-gear period [pot_closure_end + 1, est_date_end]
# A zero-length period is dropped.
#
# Backward compatibility: when pot_closure_start == est_date_start (the historical
# assumption, and the 2024-25 config), there is no pre-closure period and the result
# is the historical two-element list with the SAME internal names ("ring_net_only",
# "all_gear"), so keys, Stan fit labels, and output filenames are unchanged. Only
# when a mid-season closure is configured do the two all-gear periods take distinct
# names ("all_gear_pre" / "all_gear_post") to keep the population keys unique.
#
# Defaults: if params$pot_closure_start / _end are absent, they fall back to the old
# derivation (season start, and the day before params$pot_open_date), so an old config
# still runs. Each element also carries gear-model fields (period_bss, gear_exclude)
# that are inert for the pooled model.
#
# The model logic keys off gear_regime ("pot_closure" | "all_gear"), never off name.
###############################################################################

build_subseasons <- function(params) {
  est_start <- as.Date(params$est_date_start)
  est_end   <- as.Date(params$est_date_end)
  pc_start  <- as.Date(params$pot_closure_start %||% params$est_date_start)
  pc_end    <- as.Date(params$pot_closure_end   %||% (as.Date(params$pot_open_date) - 1))

  if (pc_start < est_start) pc_start <- est_start
  if (pc_end   > est_end)   pc_end   <- est_end
  if (pc_end < pc_start)
    stop("build_subseasons(): pot_closure_end is before pot_closure_start.", call. = FALSE)

  has_pre  <- pc_start > est_start
  has_post <- pc_end   < est_end

  allgear <- function(nm, s, e, disp) list(
    name = nm, display_name = disp, gear_regime = "all_gear",
    start = s, end = e, period_bss = "month", gear_exclude = character(0))
  closure <- list(
    name = "ring_net_only", display_name = "Pot closure", gear_regime = "pot_closure",
    start = pc_start, end = pc_end, period_bss = "biweekly", gear_exclude = c("Pot"))

  ss <- list()
  if (has_pre)
    ss <- c(ss, list(allgear(if (has_post) "all_gear_pre" else "all_gear",
                             est_start, pc_start - 1,
                             if (has_post) "All gear (pre-closure)" else "All gear")))
  ss <- c(ss, list(closure))
  if (has_post)
    ss <- c(ss, list(allgear(if (has_pre) "all_gear_post" else "all_gear",
                             pc_end + 1, est_end,
                             if (has_pre) "All gear (post-closure)" else "All gear")))
  ss
}
