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
# bss_ar_resolution.R
#
# Shared AR(1) temporal-resolution selector for the crab BSS models. Extracted
# from the pooled driver (v6.5 B1.2 cap + v7.4 ar_force toggle) and generalized
# so the gear-resolved track can call it too. Auto-sourced by both drivers via
# the 03_R_functions walk.
#
# WHY THIS EXISTS
#   The AR(1) latent process needs enough observations distributed across its
#   time steps to identify phi and sigma_eps. Too fine a resolution and the
#   process is unidentified (the pooled boat's daily AR diverged on ~100% of
#   iterations at a 289-state latent dimension); too coarse and real temporal
#   structure is smoothed away. The selector picks the finest resolution the
#   effort series supports, then applies a per-population cap.
#   References: Vehtari et al. (2021) for the diagnostics that exposed the
#   failure; Betancourt (2017) on divergences as a geometry symptom.
#
# TWO MODES
#   fixed_resolution = "<name>"  -> use that resolution verbatim (no data-driven
#       selection, no cap). This is the gear-resolved track's default: each
#       sub-season declares its own period_bss (biweekly for ring-net, monthly
#       for all-gear). Passing an explicit resolution is an explicit choice and
#       is NOT silently overridden by ar_max_resolution.
#   fixed_resolution = NULL      -> data-driven selection (daily / weekly /
#       monthly), then capped by params$ar_max_resolution[[population_name]].
#       This is the pooled track's behavior. The cap entry may be a scalar or a
#       per-gear_regime named list (see gear_regime below).
#
#   params$ar_force[[population_name]] overrides BOTH modes. It is an experiment
#   toggle; NULL (production) is a no-op.
#
#   gear_regime (optional) selects a per-sub-season cap when the population's
#   ar_max_resolution entry is a named list (e.g. list(all_gear = "daily",
#   pot_closure = "biweekly")). NULL, or a scalar cap entry, preserves the
#   original per-population behavior. Only consulted in the data-driven branch.
#
# RESOLUTIONS
#   "daily" > "weekly" > "biweekly" > "monthly" (finest to coarsest).
#   "biweekly" is reachable only via fixed_resolution or ar_force; the
#   data-driven branch never selects it, matching the pooled selector's
#   original three-way choice.
#
# RETURNS a list:
#   resolution     final resolution name (normalized)
#   P_n            number of AR periods
#   pvec           integer day -> period index, length nrow(days)
#   source         "fixed" | "adaptive" | "forced"
#   coverage       fraction of days with an effort observation
#   n_effort_days  distinct days with an effort observation
#   obs_per_week   mean effort observations per week
#
# REQUIRES days to carry day_index, week_index, month_index (prep_days_crab
# provides all three in both drivers), and eff_d to carry day_index.
###############################################################################

# Finest to coarsest. Used for the cap comparison: a resolution is reduced only
# when it is FINER than the population's cap.
.bss_res_rank <- c(monthly = 1L, biweekly = 2L, weekly = 3L, daily = 4L)

# Accept the drivers' historical spellings ("week"/"month") alongside the
# canonical adverbial names.
.bss_normalize_resolution <- function(x) {
  if (is.null(x) || is.na(x)) return(NA_character_)
  x <- tolower(trimws(as.character(x)))
  switch(x,
    "day"      = ,
    "daily"    = "daily",
    "week"     = ,
    "weekly"   = "weekly",
    "biweek"   = ,
    "biweekly" = "biweekly",
    "month"    = ,
    "monthly"  = "monthly",
    stop("Unrecognized AR resolution: '", x,
         "'. Expected one of daily, weekly, biweekly, monthly.", call. = FALSE)
  )
}

bss_select_ar_resolution <- function(days, eff_d, population_name, params,
                                     fixed_resolution = NULL,
                                     gear_regime = NULL,
                                     verbose = TRUE) {

  D <- nrow(days)

  # --- Effort-density diagnostics (reported in both modes) -------------------
  n_effort_days <- dplyr::n_distinct(eff_d$day_index)
  coverage      <- if (D > 0) n_effort_days / D else NA_real_
  n_weeks       <- max(days$week_index, na.rm = TRUE)
  obs_per_week  <- if (n_weeks > 0) n_effort_days / n_weeks else NA_real_

  # --- Choose the resolution -------------------------------------------------
  if (!is.null(fixed_resolution)) {
    # Explicit per-sub-season choice. Not subject to the cap: an explicit
    # instruction should not be silently overridden.
    ar_resolution <- .bss_normalize_resolution(fixed_resolution)
    sel_source    <- "fixed"

  } else {
    # Data-driven: finest resolution the effort series can support.
    min_cov  <- params$ar_daily_min_coverage  %||% 0.25
    min_obs  <- params$ar_daily_min_obs       %||% 20
    min_week <- params$ar_weekly_min_per_week %||% 1.5

    if (isTRUE(coverage >= min_cov) && n_effort_days >= min_obs) {
      ar_resolution <- "daily"
    } else if (isTRUE(obs_per_week >= min_week) && n_weeks >= 3) {
      ar_resolution <- "weekly"
    } else {
      ar_resolution <- "monthly"
    }
    sel_source <- "adaptive"

    # Per-population cap (v6.5 / B1.2). Reduces the latent AR dimension where the
    # effort series cannot identify a finer process. Only ever coarsens.
    # The cap entry may be a scalar resolution OR a named list keyed by gear_regime
    # ("all_gear" | "pot_closure") for a per-sub-season cap; a missing regime falls
    # back to a "default" key, then to "daily".
    cap_entry <- params$ar_max_resolution[[population_name]]
    raw_cap   <- if (is.list(cap_entry)) {
      regime_key <- if (!is.null(gear_regime)) gear_regime else "default"
      cap_entry[[regime_key]] %||% cap_entry[["default"]] %||% "daily"
    } else {
      cap_entry %||% "daily"
    }
    pop_cap <- .bss_normalize_resolution(raw_cap)
    if (.bss_res_rank[[ar_resolution]] > .bss_res_rank[[pop_cap]]) {
      ar_resolution <- pop_cap
      if (verbose) {
        cat(sprintf("  AR resolution capped to '%s' for %s/%s (ar_max_resolution)\n",
                    ar_resolution, population_name, gear_regime %||% "all"))
      }
    }
  }

  # --- Experiment override: bypasses both modes and the cap ------------------
  if (!is.null(params$ar_force) && !is.null(params$ar_force[[population_name]])) {
    ar_resolution <- .bss_normalize_resolution(params$ar_force[[population_name]])
    sel_source    <- "forced"
    if (verbose) {
      cat(sprintf("  AR resolution FORCED to '%s' for %s (ar_force experiment override)\n",
                  ar_resolution, population_name))
    }
  }

  # --- Map resolution -> period count and day -> period index ----------------
  if (ar_resolution == "daily") {
    P_n  <- D
    pvec <- as.integer(days$day_index)
  } else if (ar_resolution == "weekly") {
    pvec <- as.integer(days$week_index)
    P_n  <- max(pvec, na.rm = TRUE)
  } else if (ar_resolution == "biweekly") {
    pvec <- as.integer(ceiling(days$day_index / 14))
    P_n  <- max(pvec, na.rm = TRUE)
  } else {  # monthly
    pvec <- as.integer(days$month_index)
    P_n  <- max(pvec, na.rm = TRUE)
  }

  if (verbose) {
    cat(sprintf(
      "  AR resolution: %s [%s] (P_n=%d) | coverage=%.0f%% (%d/%d days), %.1f obs/week\n",
      ar_resolution, sel_source, P_n, coverage * 100, n_effort_days, D, obs_per_week))
  }

  list(
    resolution    = ar_resolution,
    P_n           = P_n,
    pvec          = pvec,
    source        = sel_source,
    coverage      = coverage,
    n_effort_days = n_effort_days,
    obs_per_week  = obs_per_week
  )
}


###############################################################################
# bss_ar_ladder()  --  improvement 7 (2026-08-25): the escalation ladder.
#
# WHAT PROBLEM THIS SOLVES
#   Until now, AR resolution was decided in two steps that both happen BEFORE any
#   sampling: a coverage rule (is the effort series dense enough for daily?) and a
#   per-population cap in run_config (ar_max_resolution). The cap is the part that
#   encodes sampler behaviour, and it was hand-tuned from earlier runs: shore
#   pot-closure was capped at biweekly because it funnelled at daily on Run 1, the boat
#   at monthly because it diverged on ~100% of iterations at daily. That works, but it
#   freezes a per-season empirical finding into config, and a fit that FAILS its gate is
#   demoted straight to the Point Estimator rather than being retried somewhere it could
#   succeed.
#
#   The ladder replaces the hand-tuning with a procedure: start at the finest rung, fit,
#   put the result through the SAME convergence gate that decides PE-vs-BSS, and if it
#   fails, coarsen one rung and refit. Stop at the first rung that passes. Every
#   component then reports at the finest resolution it can actually support, decided by
#   that run's own sampler behaviour.
#
# WHAT IT COSTS
#   Every rung is a real multi-hour MCMC fit. On the 2024-25 config the two capped
#   components will burn their known-bad rungs before settling where the caps already put
#   them. That is why params$ar_escalate ships FALSE and why the ladder is auditable
#   (ar_escalation_log.csv records every attempt, its gate verdict, and its runtime).
#
# WHAT IT IS NOT
#   It is not a substitute for the coverage rule. Coverage answers "can this series
#   identify a daily process at all"; the ladder answers "does the sampler survive it".
#   With ar_escalate_respect_cap = TRUE the ladder starts from the capped rung instead of
#   the top, which is the cheap variant: no extra fits when the caps are already right,
#   but it cannot discover that a cap is too coarse.
#
# DEGENERATE RUNGS ARE DROPPED
#   A rung that yields the same period count as a finer rung already on the ladder, or
#   fewer than min_periods periods, is removed: refitting an identical model wastes hours,
#   and a 1-2 period AR is not an AR. For a 76-day pot-closure window that removes
#   nothing (daily 76 / weekly ~11 / biweekly 6 / monthly 3); for a short window it can
#   collapse the ladder to a single rung, which the caller reports.
#
# RETURNS a character vector of resolutions, finest to coarsest, length >= 1.
###############################################################################
bss_ar_ladder <- function(days, eff_d, population_name, params,
                          fixed_resolution = NULL, gear_regime = NULL,
                          min_periods = 3L) {

  # ar_force is an experiment override and outranks everything, including the ladder.
  if (!is.null(params$ar_force) && !is.null(params$ar_force[[population_name]]))
    return(.bss_normalize_resolution(params$ar_force[[population_name]]))

  base_sel <- bss_select_ar_resolution(days, eff_d, population_name, params,
                                       fixed_resolution = fixed_resolution,
                                       gear_regime = gear_regime, verbose = FALSE)

  if (!isTRUE(params$ar_escalate)) return(base_sel$resolution)

  ladder <- vapply(params$ar_escalate_ladder %||% c("daily", "weekly", "biweekly", "monthly"),
                   .bss_normalize_resolution, character(1), USE.NAMES = FALSE)
  ladder <- ladder[order(-.bss_res_rank[ladder])]   # finest first

  # Start rung. Default: the top of the ladder (the point of the feature). With
  # ar_escalate_respect_cap the ladder starts no finer than what the cap + coverage rule
  # already chose.
  if (isTRUE(params$ar_escalate_respect_cap))
    ladder <- ladder[.bss_res_rank[ladder] <= .bss_res_rank[[base_sel$resolution]]]
  if (length(ladder) == 0) ladder <- base_sel$resolution

  # Drop rungs that are degenerate or duplicate an existing rung's period count.
  D <- nrow(days)
  p_of <- function(res) {
    if (res == "daily")    D
    else if (res == "weekly")   max(as.integer(days$week_index),  na.rm = TRUE)
    else if (res == "biweekly") max(as.integer(ceiling(days$day_index / 14)), na.rm = TRUE)
    else                        max(as.integer(days$month_index), na.rm = TRUE)
  }
  pn   <- vapply(ladder, p_of, numeric(1))
  keep <- !duplicated(pn) & (pn >= min_periods)
  # Never return an empty ladder: if every rung is degenerate, keep the coarsest.
  if (!any(keep)) keep[length(keep)] <- TRUE
  ladder <- ladder[keep]

  max_att <- params$ar_escalate_max_attempts %||% length(ladder)
  utils::head(ladder, max(1L, as.integer(max_att)))
}
