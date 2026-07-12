###############################################################################
# run_config.R  --  the single control surface for a production run.
#
# This is the ONE file you edit to run the estimation. Set the RUN SELECTION
# block (which model, and whether to run the weather module), then set the
# toggles below, then launch with either:
#
#     source("run_estimation.R")          # in RStudio (Source, not Knit)
#   or
#     Rscript run_estimation.R            # from a terminal, unattended
#     Rscript run_estimation.R --model gear_resolved --weather   # CLI override
#
# How the override works: run_estimation.R injects `run_config` (defined below)
# into the render environment of the chosen .Rmd. Each model does
# `params <- modifyList(params, run_config)`, so every key listed here OVERRIDES
# the model's in-file default. As of the 2026-07-11 consolidation this file is
# the SINGLE SOURCE OF TRUTH for every user-selectable toggle: the two model
# .Rmd files no longer carry their own copies of these keys, so there is nothing
# to keep in sync. Each .Rmd keeps only its own model-internal tuning (per-fit
# sampler settings, convergence-gate thresholds, AR-selector thresholds, and a
# few model constants), which you rarely touch and which legitimately differ
# between the two models.
#
# Standalone knit: each .Rmd sources THIS file automatically when `run_config`
# is not already present (the `if (!exists("run_config")) source(...)` guard in
# its setup chunk), so knitting a model .Rmd directly in RStudio uses exactly
# the same toggles as an orchestrated run. You never have to edit the .Rmd.
###############################################################################


# ============================ RUN SELECTION ================================ #
#            ^^^^ edit these two lines for a routine run ^^^^

model       <- "pooled"        # "pooled"  or  "gear_resolved"

run_weather <- FALSE           # TRUE also runs the weather-tide covariate
                               # module AFTER the model. Only valid with
                               # model = "pooled" (the weather module reuses the
                               # pooled run's in-memory objects). run_estimation.R
                               # will stop early if you set TRUE with
                               # "gear_resolved", before any multi-hour fit.

# =========================================================================== #


# ===================== USER TOGGLES (single source) ======================== #
# Everything a user changes season to season or to steer the model lives here.
# Keys are applied to whichever model runs; a key a model does not read is simply
# ignored by that model (harmless), so the model-specific toggles at the bottom
# can sit in this one shared list without affecting the other model.

run_config <- list(

  # --- Identifiers ---------------------------------------------------------
  # These unify the two models onto one set of strings. The committed gear-
  # resolved driver used "Rec Crab Grays Harbor Westport 2024-25" while pooled
  # used the string below; centralizing here makes both use one value. If either
  # string must be preserved verbatim as an output identifier, give that model
  # its own value in its .Rmd instead.
  project_name      = "Coastal Recreational Crab",
  fishery_name      = "Rec Crab Grays Harbor 2024-25",

  # --- Season window (the values you change most often) --------------------
  est_date_start    = "2024-09-16",   # first day of the estimation window
  est_date_end      = "2025-09-15",   # last day (commercial pot closure)
  season_filter     = "2024-25",

  # --- Regulatory / structural dates ---------------------------------------
  pot_open_date     = "2024-12-01",   # pots legal from this date (sub-season split)
  commercial_opener = "2025-01-01",   # (was malformed "2025-01-1" in gear-resolved)
  census_start_date = "2024-12-01",
  census_end_date   = "2025-02-08",

  # --- Catch groups --------------------------------------------------------
  estimate_red_rock = FALSE,          # TRUE adds Red_Rock_Kept alongside Dungeness

  # --- Day-type / stratification -------------------------------------------
  days_wkend        = c("Friday", "Saturday", "Sunday"),
  min_fishing_time  = 0.5,            # min crabber-hours to keep an interview
  period_pe         = "week",         # PE temporal stratum
  sections          = c(1),

  # --- Incomplete-trip filter (both models) --------------------------------
  # Incomplete trips (soak-time gear not yet retrieved) read systematically low
  # (about -20% CPUE for pots/traps), biasing CPUE and hence the harvest estimate
  # low. TRUE excludes them from CPUE estimation (PE and BSS), keeping Complete +
  # NA; FALSE keeps all trips (pre-filter behavior). Missing trip_status is kept.
  filter_incomplete_trips = TRUE,

  # --- Effort unit (both models) -------------------------------------------
  # Shore and boat CPUE denominators. As of pooled v7.7 / gear-resolved v5.5 both
  # components run on gear-DEPLOYMENTS: the pipeline's own linearity diagnostic
  # flags every time-denominated unit as invalid for pots (shore beta_h 0.57 for
  # crabber-hours, 0.73 for gear-hours, 1.05 for deployments; deployments is the
  # only harvest-unbiased unit). Routed through 03_R_functions/bss_effort_spec.R
  # so the BSS and PE always share a unit. Set shore_effort_unit = "crabber-hours"
  # to revert shore only.
  shore_effort_unit      = "gear-deployments",  # "crabber-hours" | "gear-hours" | "gear-deployments"
  tau_shore_prior_mu     = 1.7,       # shore deployment turnover (trips/gear-slot/day)
  tau_shore_prior_sigma  = 0.3,
  tau_boat_prior_mu      = 1.2,       # boat deployment turnover
  tau_boat_prior_sigma   = 0.3,
  gear_per_group_default = 4.0,       # PE fallback gear-per-boat-group when no interview records it

  # --- BSS run-level settings (NOT per-fit tuning) -------------------------
  bss_chains        = 4,
  bss_cores         = 4,
  bss_seed          = 20260619,       # fixed seed for reproducible fits
  bss_max_count_seq = 3,              # cap on count sequences per day

  # --- AR resolution experiment toggle -------------------------------------
  # PRODUCTION VALUE IS NULL. Forces a population's AR resolution, bypassing both
  # the data-driven selection and the per-population cap for the named population.
  # To run the boat daily-vs-weekly experiment, set list(private_boat = "daily").
  ar_force          = NULL,

  # --- Ingress/egress input + shore day length (both models) ---------------
  # Shore effort is expanded by the I/E-derived effective day length (~3.5-5 h),
  # not civil twilight (9-17 h). Fallback is automatic: regression -> grand mean
  # -> civil twilight (only when there is effectively no I/E data, with a
  # warning). Boats always use L = 24 h (gear soaks continuously).
  ie_data_file      = "ingress_egress.xlsx",
  ie_sheet          = "data",
  ie_shore_location = "WDF20",
  ie_boat_location  = "WBL",
  use_ie_day_length = TRUE,
  ie_min_obs_for_regression = 5,

  # Civil-twilight clamp. Binds only on the fallback rung and on the
  # day_length_civil_twilight diagnostic column.
  day_length_min_hours = 9.0,
  day_length_max_hours = 17.0,

  # --- Crabbing holidays (single source of truth for all three drivers) ----
  # High-effort non-weekend days treated as weekend for day-typing. Update this
  # ONE list each season; the models pull from here instead of their own copies.
  crabbing_holiday_dates = as.Date(c(
    "2024-11-29",  # Native American Heritage Day
    "2024-12-31",  # New Year's Eve
    "2025-01-01",  # New Year's Day
    "2025-02-08",  # Super Bowl Eve
    "2025-05-24",  # Memorial Day weekend - Saturday
    "2025-05-25",  # Memorial Day weekend - Sunday
    "2025-05-26",  # Memorial Day
    "2025-06-15",  # Father's Day
    "2025-07-04",  # Independence Day
    "2025-09-01"   # Labor Day
  )),

  # --- Model-specific toggles (centralized here; each is read only by its own
  #     model and ignored by the other, so they are safe to keep in one list) --
  collapse_mu_hier           = FALSE, # (pooled) collapse the single-cell mu-hierarchy
                                       #   (B1.7/POOL-4 experiment lever). FALSE = current
                                       #   hierarchy, posterior unchanged. Accepts a per-
                                       #   population named list, e.g. list(private_boat = TRUE).
  estimate_B1_C              = TRUE,   # (gear-resolved) weekend/holiday CPUE effect B1_C.
                                       #   TRUE matches the pooled model; FALSE drops B1_C
                                       #   from the likelihood (v5.4 behavior).
  ar_adaptive                = FALSE,  # (gear-resolved) FALSE preserves the fixed per-sub-
                                       #   season period_bss (biweekly ring-net, monthly all-
                                       #   gear) EXACTLY. TRUE hands AR choice to the data-driven
                                       #   selector; that is inference-changing, so validate first.
  loo_effort_unit_comparison = FALSE,  # (gear-resolved) TRUE restricts interviews to the common
                                       #   valid-denominator subset for a legitimate cross-unit
                                       #   elpd_loo comparison. FALSE for a production run of the
                                       #   chosen unit (the comparison is done).
  use_boat_ie                = TRUE    # (gear-resolved) use WBL boat I/E ingress counts to
                                       #   identify tau once enough days exist. IE_n = 0 is safe.
)

# =========================================================================== #
