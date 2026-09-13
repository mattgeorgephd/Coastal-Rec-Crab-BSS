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
# bss_superseded_runner.R  --  a batch runner whose question is CLOSED refuses to run
# -----------------------------------------------------------------------------
# THE PROBLEM THIS SOLVES, measured rather than assumed. 06_diagnostics/ holds eleven
# batch runners. A static audit on 2026-09-13 asked how many of the fifteen levers that
# define Method v2.0 each one PINS or SETS:
#
#   run_improvements_2026-09-08.R  12/15   (it has the WINDOW pin, added by B22)
#   run_improvement_plan_2026-08-27.R  5/15
#   run_stage5_2026-08-30.R            4/15
#   every other runner               0-2/15
#
# All eleven still parse. All eleven still run: sourced today in DRY_RUN they complete
# without error. That is what makes them dangerous rather than merely stale. A runner
# written before the dynamic f, the calibrated turnovers and the census split will fit its
# "pre-patch baseline" rung with TODAY's f, TODAY's tau and TODAY's census, because it
# never pins them, and it will then label the result as the baseline it is not. Nothing in
# its output says so. A broken runner announces itself; this class does not.
#
# WHAT THIS IS NOT. It is not a claim that the runner is worthless. Each one is the record
# of how a decision was reached, several are the ONLY provenance for their output folders
# (the 20260825 renders carry no run_parameters.txt at all, because the driver did not
# write one yet), each is referenced by 7 to 15 documents, and each is cited by the
# harness. Deleting them would dangle those references and destroy the provenance to save
# about 400 KB. So they stay, readable, and they refuse to FIT.
#
# THE OVERRIDE IS DELIBERATE. A person who wants to re-render a historical batch, having
# read why they should not trust it as a comparison, sets the override and proceeds. The
# point is to make that a decision instead of an accident.
#
#   I_KNOW_THIS_IS_SUPERSEDED <- TRUE
#   source("06_diagnostics/run_stage5_2026-08-30.R")
#
# Modelled on the hard stop run_estimation.R gained for the retired --weather flags: a
# retired path that silently does something reasonable is worse than one that stops.
###############################################################################

bss_superseded_runner <- function(runner, question, settled_by, what_would_happen,
                                  levers_pinned = NULL) {
  ovr <- tryCatch(get("I_KNOW_THIS_IS_SUPERSEDED", envir = globalenv()),
                  error = function(e) NULL)
  msg <- paste0(
    "\n", strrep("=", 78), "\n",
    " SUPERSEDED RUNNER: ", runner, "\n",
    strrep("=", 78), "\n",
    " THE QUESTION IT WAS BUILT TO ANSWER IS CLOSED.\n",
    "   question   : ", question, "\n",
    "   settled by : ", settled_by, "\n",
    if (!is.null(levers_pinned))
      paste0("   pins       : ", levers_pinned, " of the 15 Method v2.0 levers\n") else "",
    "\n WHAT WOULD HAPPEN IF YOU RAN IT ANYWAY:\n   ",
    what_would_happen, "\n",
    "\n THE FILE IS KEPT ON PURPOSE. It is the record of how that decision was reached, it\n",
    " is cited by the development notes and by the harness, and for the earliest batches it\n",
    " is the only surviving description of what its output folders were configured as.\n",
    " Reading it is the point; fitting it is not.\n",
    "\n TO RUN IT ANYWAY, having read the above:\n",
    "   I_KNOW_THIS_IS_SUPERSEDED <- TRUE\n",
    "   source(\"", runner, "\")\n",
    "\n THE CURRENT LADDER is 06_diagnostics/run_improvements_2026-09-08.R (pooled and the\n",
    " gear cross-check) and 06_diagnostics/run_gear_ar_zi_2026-09-13.R (D3 and D6). Both\n",
    " carry a WINDOW pin and a preflight that FAILS on any undeclared difference between\n",
    " rungs, which is the thing this runner lacks.\n",
    strrep("=", 78), "\n")
  if (isTRUE(ovr)) {
    cat(msg)
    cat(" I_KNOW_THIS_IS_SUPERSEDED is TRUE. Proceeding. Do not report a number from this\n",
        " run as a comparison against anything fitted after its own date.\n\n")
    return(invisible(FALSE))
  }
  stop(msg, call. = FALSE)
}
