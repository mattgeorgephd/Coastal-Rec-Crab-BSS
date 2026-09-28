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
# build_all_inputs.R  (2026-09-10)
#
# Rebuilds every model workbook in 04_input_files/ from the per-season creel workbooks
# in raw/ (<YYYY><YY>_rec_crab_harvest_data.xlsx), in dependency order. A new season is:
# drop its workbook into raw/, run this script from the repository root, check the
# builders' reports, commit the workbook and the rebuilt inputs together
# (07_documentation/NEW_SEASON_GUIDE.md, step 1).
#
#   Rscript 04_input_files/build_all_inputs.R
#
# Not rebuilt here (their sources are not the season workbooks): ingress_egress.xlsx
# (the I/E database export; hand-maintained, see ../README.md), WBL_boat_counts.xlsx
# (OSP), fishery_opener_dates.xlsx (the regulation calendar).
###############################################################################

builders <- c(
  "build_interview_combined.R",     # interviews, every season           -> interview_combined.xlsx
  "build_effort_combined.R",        # instantaneous effort counts        -> effort_combined.xlsx
  "build_sampler_shifts.R",         # sampler site-visits (shifts)       -> sampler_shifts.xlsx
  "build_comm_charter_tally.R",     # Westport daily vessel tally        -> wes_commercial_tally.xlsx
  "build_charter_trips.R",          # charter trip roster                -> charter_trips.xlsx
  "build_crabbing_holidays.R"       # holiday calendar (rule + check against the shifts) -> crabbing_holidays.xlsx
)
dir <- if (dir.exists("04_input_files")) "04_input_files" else "."
for (b in builders) {
  cat("\n", strrep("=", 100), "\n", b, "\n", strrep("=", 100), "\n", sep = "")
  st <- system2(file.path(R.home("bin"), "Rscript"), file.path(dir, b))
  if (!identical(st, 0L)) stop(b, " failed (exit status ", st, ")")
}
cat("\nAll inputs rebuilt.\n")
