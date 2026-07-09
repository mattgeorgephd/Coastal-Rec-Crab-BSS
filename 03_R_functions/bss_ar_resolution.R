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
#       This is the pooled track's behavior.
#
#   params$ar_force[[population_name]] overrides BOTH modes. It is an experiment
#   toggle; NULL (production) is a no-op.
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
    pop_cap <- .bss_normalize_resolution(
      params$ar_max_resolution[[population_name]] %||% "daily"
    )
    if (.bss_res_rank[[ar_resolution]] > .bss_res_rank[[pop_cap]]) {
      ar_resolution <- pop_cap
      if (verbose) {
        cat(sprintf("  AR resolution capped to '%s' for %s (ar_max_resolution)\n",
                    ar_resolution, population_name))
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
