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
# bss_gear_period.R  --  the gear track's AR period, PER POPULATION
# -----------------------------------------------------------------------------
# WHY THIS EXISTS, and it is a correction rather than a feature.
#
# The gear-resolved driver runs ar_adaptive = FALSE, so it passes
# fixed_resolution = <the sub-season's period> to bss_select_ar_resolution() and the
# per-population cap ar_max_resolution$gear_resolved is never consulted. period_bss is a
# property of the SUB-SEASON, and a sub-season contains BOTH populations. So on this track
# the shore all-gear fit and the boat all-gear fit were structurally locked to the SAME AR
# period, with no way to separate them.
#
# THAT IS WHAT BROKE THE FIRST D3 LADDER (run 2026-09-14). Moving
# gear_period_bss$all_gear from "month" to "weekly" was meant to move the shore fit, which
# is what CHANGE_REGISTER D3 asks about. It moved both:
#
#   shore all-gear   28,334 -> 29,266   (+3.3%)
#   boat  all-gear   45,374 -> 53,584   (+18.1%)   <- not what D3 asked about
#   port             93,278 -> 102,430  (+9.8%)
#
# so 90% of the port movement came from a component the question was not about, and the
# ladder's recommendation of "weekly" would, if adopted through that key, have shipped the
# WORST of the four configurations measured: +8.62% from the pooled track, against -0.35%
# for the shipped one. This is the same class of defect as the 2026-08-27 Stage C ar_force
# bug, which forced both boat sub-seasons to biweekly and "made its port total
# uninterpretable as the change it was supposed to isolate". One layer up, same shape.
#
# WHAT THE SEPARATION BUYS, measured from the committed draws of that run rather than
# predicted (the four fits are independent by design, so summing one rung's shore draws
# with another's boat draws is exactly the operation the driver performs within a run):
#
#   configuration                              port      vs pooled R4
#   gear shore monthly + boat monthly (shipped)  93,680      -0.35%
#   gear shore weekly  + boat weekly            102,120      +8.62%
#   gear shore weekly  + boat monthly            94,278      +0.28%   <- MATCHED
#
# The pooled track fits shore all-gear at WEEKLY and boat all-gear at MONTHLY, because it
# uses the per-population cap. +0.28% at the PORT is the closest the two tracks have ever
# come; the previous best was 0.08% on the shore component alone. The -1.39% gap that
# opened D3 is therefore not a shore-resolution artefact: it is the gear track being
# unable to express a per-population configuration at all.
#
# THE SHAPE. Either form is accepted, and the flat one is what ships, so nothing moves:
#
#   flat, per gear_regime (SHIPPED, and what build_subseasons() reads):
#     gear_period_bss = list(all_gear = "month", pot_closure = "biweekly")
#
#   per population, then per gear_regime (what D3 needs to adopt):
#     gear_period_bss = list(
#       shore        = list(all_gear = "weekly", pot_closure = "biweekly"),
#       private_boat = list(all_gear = "month",  pot_closure = "biweekly"))
#
# Resolution order: population -> gear_regime -> "default" -> the sub-season's own
# period_bss (which build_subseasons() already set from the flat form). Returning NULL
# means "no per-population opinion", and the caller keeps ss$period_bss, so a config that
# has never heard of this function behaves exactly as before.
#
# NOT ar_force. That key does the same thing and outranks everything, but it is documented
# as an experiment toggle whose production value is NULL, and it bypasses the cap silently.
# A ladder should use ar_force (06_diagnostics/run_gear_ar_zi_2026-09-13.R does); a shipped
# configuration should say what it means in a key a reader will find next to the other AR
# settings.
###############################################################################

bss_gear_period <- function(params, population_name, gear_regime = NULL) {
  g <- params$gear_period_bss
  if (is.null(g) || !is.list(g)) return(NULL)
  # per-population form: the top level names populations
  entry <- g[[population_name]]
  if (is.null(entry)) {
    # flat form. build_subseasons() has already applied it to ss$period_bss, so there is
    # no per-population opinion to add and the caller keeps what it has.
    return(NULL)
  }
  if (!is.list(entry)) return(as.character(entry))   # scalar: every sub-season of this population
  key <- if (!is.null(gear_regime)) gear_regime else "default"
  val <- entry[[key]] %||% entry[["default"]]
  if (is.null(val)) return(NULL)
  as.character(val)
}

# Is gear_period_bss in the per-population form? Used by the drivers to say so in the run
# log, because a silent per-population override is exactly the kind of thing that should
# not be silent.
bss_gear_period_is_per_pop <- function(params) {
  g <- params$gear_period_bss
  is.list(g) && any(c("shore", "private_boat") %in% names(g))
}
