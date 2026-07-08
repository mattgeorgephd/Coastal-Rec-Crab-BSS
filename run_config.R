###############################################################################
# run_config.R  --  the single control surface for a production run.
#
# This is the ONE file you edit to run the estimation. Set the RUN SELECTION
# block (which model, and whether to run the weather module), set the SEASON
# PARAMETERS, then launch with either:
#
#     source("run_estimation.R")          # in RStudio (Source, not Knit)
#   or
#     Rscript run_estimation.R            # from a terminal, unattended
#     Rscript run_estimation.R --model gear_resolved --weather   # CLI override
#
# How the override works: run_estimation.R injects `run_config` (defined below)
# into the render environment of the chosen .Rmd. Each model does
# `params <- modifyList(params, run_config)`, so the keys listed here OVERRIDE
# the model's in-file defaults, and every key you do NOT list is left at the
# model's own default (per-fit sampler tuning, convergence-gate thresholds,
# I/E flags, weather stations, and so on). That is deliberate: this file holds
# only what must be identical across models and what you change season to
# season. It is not a dump of every parameter.
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


# ======================= SEASON / RUN PARAMETERS =========================== #
# Change these when you move to a new season. Everything here is applied to the
# model (and, under Option A, to the weather module) as a run-level override.

run_config <- list(

  # --- Identifiers ---------------------------------------------------------
  # NOTE: these unify the two models. The committed gear-resolved driver used
  # "Rec Crab Grays Harbor Westport 2024-25" while pooled used the string below;
  # centralizing here makes both use one value. If either string feeds an output
  # identifier you need to preserve verbatim, drop these two keys and let each
  # model keep its own.
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

  # --- BSS run-level sampler settings (NOT per-fit tuning) ------------------
  bss_chains        = 4,
  bss_cores         = 4,
  bss_seed          = 20260619,       # fixed seed for reproducible fits
  bss_max_count_seq = 3,              # cap on count sequences per day

  # --- AR resolution experiment toggle -------------------------------------
  # PRODUCTION VALUE IS NULL. To run the experiment on purpose, 
  # set this to list(private_boat = "daily").
  ar_force          = NULL,

  # --- Ingress/egress input (pooled model only; harmless to the others) ----
  ie_data_file      = "ingress_egress.xlsx",
  ie_sheet          = "data",

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
  ))
)

# =========================================================================== #
