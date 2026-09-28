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
# pe_monthly_split.R  (shared by the pooled and gear-resolved drivers and by
# save_run_diagnostics.R)
#
# The PE component total split into calendar months, for a component that reports its PE
# (the drivers' 7.8 / 7.8b monthly blocks) and for the PE column of monthly_pe_vs_bss.csv.
#
# HOW (2026-09-28, B44). run_pe_pooled() and run_pe_gear() build the estimate stratum by
# stratum: a (period x day_type) cell carries est_total = mean daily effort x its calendar
# days, and est_catch = est_total x that cell's OWN ratio-of-sums CPUE (or the empty-stratum
# fill). Both runners keep those frames (effort_strata, catch_strata[[cg]]). A cell's
# estimate is spread evenly over its own calendar days, and the days are summed by month,
# so the months add up to the component totals EXACTLY and each month's catch is carried at
# the CPUE of the cells that make it up.
#
# WHAT THIS REPLACED, AND WHY. pe_monthly_effort_share() (drivers) and
# .srd_monthly_share() (monthly_pe_vs_bss.csv) split the component's catch TOTAL by each
# month's share of raw effort on the SAMPLED days. Two defects followed:
#   1. one CPUE for every month. Catch = total x effort share assumes the catch rate is
#      constant across the sub-season; crab CPUE has a strong within-season gradient
#      (pe_effort_strata.R), so catch moved from high-CPUE months into low-CPUE ones;
#   2. sampled days only. A month's weight was its sampled effort, not its expanded
#      effort, so a thinly sampled month was under-weighted however many calendar days
#      (and imputed cells) it carried.
# Each was also a SECOND formula for the PE's daily effort, which had to be kept in step
# with run_pe_*() by hand: it drifted twice (day length, 2026-08-25; the crabbing fraction,
# 2026-09-12). Reading the PE's own strata removes that class of defect: there is no
# second formula left to drift.
###############################################################################

# pe_monthly_split()
#   pe_res        one run_pe_pooled() / run_pe_gear() result (needs effort_strata; uses
#                 catch_strata[[catch_group]] when present, else the catch is 0, which is
#                 what the component total is when the catch group has no column)
#   days_ss       the sub-season calendar the PE was run on (event_date, period, day_type,
#                 open_section_1), i.e. prep_days_crab() for that sub-season
#   catch_group   e.g. "Dungeness_Kept"
# Returns one row per month: month_label ("%Y-%m"), month_effort, month_catch. NULL when
# there are no strata.
pe_monthly_split <- function(pe_res, days_ss, catch_group = "Dungeness_Kept") {
  es <- pe_res$effort_strata
  if (is.null(es) || !nrow(es)) return(NULL)
  cal <- days_ss
  if ("open_section_1" %in% names(cal)) cal <- cal[cal$open_section_1 %in% TRUE, , drop = FALSE]
  cal <- data.frame(period = cal$period, day_type = cal$day_type,
                    month_label = format(as.Date(cal$event_date), "%Y-%m"), stringsAsFactors = FALSE)
  per_day <- function(frame, col) {
    x <- data.frame(section_num = frame$section_num, period = frame$period, day_type = frame$day_type,
                    v = ifelse(is.finite(frame[[col]]), frame[[col]], 0) / pmax(frame$n_total_days, 1),
                    stringsAsFactors = FALSE)
    z <- merge(x, cal, by = c("period", "day_type"))
    if (!nrow(z)) return(stats::setNames(numeric(0), character(0)))
    tapply(z$v, z$month_label, sum)
  }
  eff <- per_day(es, "est_total")
  cs  <- pe_res$catch_strata[[catch_group]]
  cat_ <- if (is.null(cs) || !nrow(cs)) stats::setNames(rep(0, length(eff)), names(eff)) else per_day(cs, "est_catch")
  months <- sort(unique(c(names(eff), names(cat_))))
  data.frame(month_label  = months,
             month_effort = as.numeric(ifelse(is.na(eff[months]), 0, eff[months])),
             month_catch  = as.numeric(ifelse(is.na(cat_[months]), 0, cat_[months])),
             stringsAsFactors = FALSE)
}
