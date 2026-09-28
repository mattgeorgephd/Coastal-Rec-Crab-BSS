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
# bss_rng.R  (2026-09-28, B46; shared by both drivers and the diagnostics)
#
# bss_with_seed(seed, expr): evaluate `expr` under set.seed(seed) and RESTORE the caller's
# random-number state afterwards (or leave none, if there was none). A NULL seed evaluates
# `expr` on the caller's stream, unseeded, and still restores nothing it did not touch.
#
# Why one helper. Several functions called set.seed() and left the global stream where their
# own draws had moved it (estimate_shore_turnover, bss_ppc_calibration, the run diagnostics'
# draw subsample, gear_share_bootstrap, the gear prep's interview subsample, the drivers'
# census draw), and two drew unseeded (save_bss_ppc_draws' thinning and the monthly PE-vs-BSS
# subsample). So any random draw made AFTER one of them depended on call order: turning a
# diagnostic on or off moved every later draw, the defect the 2026-09-07 ladder tripped over
# (bss_rung_adequacy.R) and B38 / B43 fixed one site at a time. Seeding here reproduces the
# values each site drew before (same seed, same sequence of draws inside the block); what
# changes is only that the stream AFTER the block is the caller's again.
###############################################################################
bss_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  if (!is.null(seed)) set.seed(as.integer(seed))
  expr   # forced here, after the seed
}
