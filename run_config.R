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
# How it works: run_estimation.R injects `run_config` (defined below) into the
# render environment of the chosen .Rmd. As of the 2026-07-12 restructure (P5),
# run_config is the BASE parameter set and each model layers its own internal
# tuning on top: each .Rmd does `params <- modifyList(run_config, params_model)`,
# where `params_model` holds ONLY that model's specifics (Stan file, per-fit
# sampler settings, gate thresholds, AR-selector thresholds, model constants).
# The two key sets are disjoint, so the merge order carries intent only.
# This file is the SINGLE SOURCE OF TRUTH for every user-selectable toggle,
# including the AR resolution map below; the two .Rmd files no longer carry their
# own copies, so there is nothing to keep in sync. The AR map is per-model (the
# two models legitimately differ), so each driver selects its slice
# (run_config$ar_max_resolution$<model>) right after the merge.
#
# Standalone knit: each .Rmd sources THIS file automatically when `run_config`
# is not already present (the `if (!exists("run_config")) source(...)` guard in
# its setup chunk), so knitting a model .Rmd directly in RStudio uses exactly
# the same toggles as an orchestrated run. You never have to edit the .Rmd.
###############################################################################


# ============================ RUN SELECTION ================================ #
#            ^^^^ edit these two lines for a routine run ^^^^

model       <- "gear_resolved"        # "pooled"  or  "gear_resolved"

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

  # --- Run label -----------------------------------------------------------
  # Optional short label for THIS run's output subfolder, so same-day runs of a
  # model no longer overwrite each other. Blank = the driver auto-appends a HHMMSS
  # timestamp (folder like 05_output/<date>/pooled-CPUE-143022). Set a meaningful
  # string (e.g. "run5") for a named folder like pooled-CPUE-run5. run_rg_sweep.R
  # sets this per run automatically.
  run_tag           = "gear_resolved_G = TRUE",

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
  est_date_end      = "2025-09-15",   # last day
  season_filter     = "2024-25",

  # --- Regulatory / structural dates ---------------------------------------
  pot_open_date     = "2024-12-01",   # pots legal from this date (used for L_effective
                                      #   and as the default pot_closure_end + 1 day)
  # Pot-closure window: the period when pots are NOT legal (only non-pot gear, ring
  # nets/snares/traps). Given explicitly here rather than assumed to start at the
  # season start, so a future season whose start does not coincide with the closure
  # start is supported. Outside this window pots are allowed (all-gear). If a closure
  # starts after est_date_start or ends before est_date_end, the driver adds the
  # corresponding all-gear period(s) automatically (see 03_R_functions/build_subseasons.R).
  # Keep pot_open_date = pot_closure_end + 1.
  pot_closure_start = "2024-09-16",   # first day of the pot closure (here = season start)
  pot_closure_end   = "2024-11-30",   # last day pots are illegal (day before pots open)
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
  # PE empty-stratum CPUE fallback (item 2, 2026-07-13). A week x day-type stratum with
  # expanded effort but no surviving interviews: "pooled" (default) borrows the
  # population x sub-season ratio-of-sums CPUE; "zero" assigns it zero catch (old
  # behavior). "pooled" removes the sparse-stratum sign instability in the thin boat PE
  # (the incomplete-trip "anomaly", item 2); shore is dense so it barely moves.
  pe_empty_stratum  = "pooled",       # "pooled" | "zero"

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

  # --- R_G prior sensitivity (T1.3; OFF by default = data-driven) ----------
  # The pooled model's R_G prior (gear per crabber) is data-driven by default. To run
  # the R_G prior-sensitivity sweep (backlog T1.3 / critique 2), uncomment R_G_prior_mu
  # and set it to each of 1.0, ~1.28 (the empirical value), and 1.5 in turn, re-running
  # each; then compare the port totals. A tighter R_G_prior_sigma makes the prior bind
  # harder. Leave commented for production. (tau_shore/tau_boat above are the analogous
  # effort-expansion priors and are already live toggles.)
  # R_G_prior_mu    = 1.27,
  # R_G_prior_sigma = 0.3,

  # --- tau_boat prior sensitivity (GR-12; single-run projection, ON by default) ---
  # The boat catch is proportional to tau_boat (the boat deployment turnover) when tau
  # is prior-dominated, which it is whenever no in-window boat I/E days pin it (the
  # 2024-25 case: 0 boat I/E days, so the boat total rests on the tau_boat prior). With
  # diagnose_tau_sensitivity = TRUE the gear driver PROJECTS the boat and port totals
  # across tau_sensitivity_grid in a single run (no refit), via
  # 03_R_functions/diagnose_tau_boat_sensitivity.R; it writes tau_boat_sensitivity.csv
  # plus an on-page table and states whether tau was prior-dominated (projection exact)
  # or not (projection = upper bound). For the EXACT multi-refit check that stays valid
  # even when boat I/E informs tau, source 06_diagnostics/run_tau_sweep.R instead.
  # Default TRUE since the 2026-07-20 multi-refit sweep confirmed the projection
  # reproduces the exact result to ~0.2% (boat elasticity 1.00), so it is ~free.
  diagnose_tau_sensitivity = TRUE,
  tau_sensitivity_grid     = c(0.9, 1.0, 1.2, 1.5, 1.8),  # tau_boat_prior_mu values to project

  # --- BSS run-level settings (NOT per-fit tuning) -------------------------
  bss_chains        = 4,
  bss_cores         = 4,
  bss_seed          = 20260619,       # fixed seed for reproducible fits
  bss_max_count_seq = 3,              # cap on count sequences per day

  # --- AR resolution experiment toggle -------------------------------------
  # PRODUCTION VALUE IS NULL. Forces a population's AR resolution, bypassing both
  # the data-driven selection and the per-population cap for the named population.
  # To run the boat daily-vs-weekly experiment, set list(private_boat = "daily").
  # To test the boat monthly-AR reconciliation (open item in the model-state
  # review), set list(private_boat = "monthly").
  ar_force          = NULL,

  # --- AR resolution caps (per-model map; each driver selects its own slice) ----
  # Cap on the finest AR resolution the data-driven selector may choose, per
  # population. The two models legitimately differ, so both maps live here and each
  # driver reads run_config$ar_max_resolution$<model> just after it merges run_config:
  #   pooled:        runs adaptive AR, capped here (boat weakly informative -> weekly).
  #   gear_resolved: reads its map only when ar_adaptive = TRUE; production
  #                  gear-resolved is ar_adaptive = FALSE (fixed period_bss), so its
  #                  map is dormant. Coarser than pooled: the gear-resolved latent AR
  #                  is P_n x (G*S), ~4x the pooled dimension with 4 gear types.
  # A population absent from a map defaults to "daily" (no cap). A population's
  # entry may be a single resolution (applied to all its sub-seasons) OR a named
  # list keyed by gear_regime ("all_gear" | "pot_closure", from build_subseasons)
  # for a PER-SUB-SEASON cap; an unlisted regime falls back to a "default" key,
  # then to "daily". The pooled shore uses this to keep daily AR on the well-
  # sampled all-gear fit while capping the thin pot-closure (ring-net) fit at
  # biweekly: the pot-closure fit funnels at daily AR (~1,165 divergences on Run 1),
  # fails its gate, and falls back to PE; biweekly removes the funnel so it reports
  # BSS, matching the gear track's biweekly ring-net period_bss. all-gear is left
  # data-driven (daily) because it fits cleanly there. The gear_resolved map takes
  # the SAME per-sub-season structure; its values mirror the gear track's fixed
  # period_bss (monthly all-gear, biweekly ring-net). NOTE the gear map is dormant
  # in production: gear-resolved runs ar_adaptive = FALSE, so fixed_resolution =
  # period_bss bypasses the cap. It is consulted only in the ar_adaptive = TRUE
  # experiment, where it now agrees with the fixed periods instead of the old
  # blanket "weekly".
  ar_max_resolution = list(
    pooled        = list(shore = list(all_gear = "daily",   
                                      pot_closure = "biweekly"),
                         private_boat = "monthly"),
    gear_resolved = list(shore = list(all_gear = "monthly", 
                                      pot_closure = "biweekly"),
                         private_boat = "monthly")
  ),

  # --- Input workbooks (all .xlsx, single "data" sheet) --------------------
  # Every pipeline input is an .xlsx workbook (converted from CSV on 2026-07-16) with
  # one "data" sheet; dates are ISO yyyy-mm-dd text. Filenames are parameters so a
  # workbook can be swapped without touching the readers. The I/E, holiday, and
  # fishery-opener workbooks keep their own filename keys in the sections below.
  effort_file       = "effort_combined.xlsx",
  interview_file    = "interview_combined.xlsx",
  tally_file        = "wes_commercial_tally.xlsx",
  input_sheet       = "data",         # sheet name shared by the flat input workbooks

  # --- Input selection: which rows/sites each reader keeps -----------------
  # Surfaced here (mirroring the ingress/egress ie_shore_location / ie_boat_location
  # below) so a different port or site set can be run without editing fetch_crab_data*.
  # A key left unset falls back to the historical Grays Harbor / Westport value.
  gh_creel_location  = "Grays Harbor",                         # interview creel_location filter
  gh_effort_areas    = c("Westport Docks Float 20","Westport Docks Float 17-21",
                         "Westport Boat Launch","Westport Marina","Westport Jetty",
                         "Ocean Shores Boat Launch","Damon Point"),   # effort creel_area whitelist
  shore_dock_float20 = "Westport Docks Float 20",              # paired shore gear-count floats
  shore_dock_float17 = "Westport Docks Float 17-21",
  boat_launch_areas  = c("Westport Boat Launch","Ocean Shores Boat Launch"),  # boat-trailer count sites

  # --- Ingress/egress input + shore day length (both models) ---------------
  # Shore effort is expanded by the I/E-derived effective day length (~3.5-5 h),
  # not civil twilight (9-17 h). Fallback is automatic: regression -> grand mean
  # -> civil twilight (only when there is effectively no I/E data, with a
  # warning). Boats always use L = 24 h (gear soaks continuously).
  ie_data_file      = "ingress_egress.xlsx",
  ie_sheet          = "data",
  ie_shore_location = "WDF20",       # location_name kept as the SHORE I/E series
  ie_boat_location  = "WBL",         # location_name kept as the BOAT I/E series
  ie_filter_by_season = FALSE,       # FALSE pools all seasons of I/E for the L_effective
                                     #   regression (historical, current behavior). TRUE
                                     #   restricts to season_filter via the workbook's
                                     #   season column (now the fishery season label).
  use_ie_day_length = TRUE,
  ie_min_obs_for_regression = 5,
  # GR-8 (2026-07-13): minimum in-window I/E days required before the I/E likelihood is
  # allowed to bind, per component. Below this the stream is dropped and sigma_IE stays
  # prior-only (decoupled), removing the sparse-data sigma_IE funnel (only 2 in-window
  # I/E days in the shore ring-net / pot-closure fit). Set to 0 to always use whatever
  # I/E data exists. The sigma_IE prior itself is left as exponential(5) on purpose;
  # tightening it would push the shore all-gear sigma_IE (~1.07) down and force
  # possibly-unrepresentative I/E days to bind harder (see the shore-I/E diagnostic).
  ie_min_obs_shore = 3,   # shore components (both models)
  ie_min_obs_boat  = 2,   # boat components (gear-resolved; boat I/E identifies tau)

  # Civil-twilight clamp. Binds only on the fallback rung and on the
  # day_length_civil_twilight diagnostic column.
  day_length_min_hours = 9.0,
  day_length_max_hours = 17.0,

  # --- Crabbing holidays (now in an editable workbook, not hardcoded here) --
  # High-effort non-weekend days treated as weekend for day-typing. These used to be
  # a hardcoded as.Date(c(...)) vector here and were duplicated in the weather module;
  # they now live in ONE editable workbook, 04_input_files/crabbing_holidays.xlsx
  # (columns: season, date, holiday_name), so the season-to-season update is a
  # spreadsheet edit and multiple seasons coexist in one file. All three drivers read
  # it via 03_R_functions/read_crabbing_holidays.R, which filters to season_filter and
  # STOPS if the file, the required columns, or the requested season are missing (so a
  # mistyped season can never silently blank out holiday day-typing). Override the
  # name/sheet with these keys.
  crabbing_holidays_file  = "crabbing_holidays.xlsx",
  crabbing_holidays_sheet = "data",

  # --- Other-fishery opener dates (spillover DIAGNOSTIC; pooled report only) ----
  # One consolidated daily OPEN/CLOSED calendar for Marine Area 2 finfish and coastal
  # razor-clam digs, read by 03_R_functions/prep_fishery_events.R and used by
  # diagnose_fishery_spillover.R to test whether crab effort/CPUE differs on those dates
  # (candidate day categories, like the crabbing holidays above). DIAGNOSTIC ONLY: it
  # reports associations and changes no estimate. Set the toggle FALSE to skip it.
  # razor_nearby_beaches are the beaches closest to Grays Harbor; note that in the
  # 2024-25 data Twin Harbors is open on every listed dig day, so the nearby flag
  # coincides with "any dig" that season (the report surfaces this overlap).
  # As of 2026-07-16 the combined sheet below is the ONLY source: the former per-fishery
  # workbooks (MA2-fishing-dates*.xlsx, razor-clam-dig-dates*.xlsx) are retired and
  # prep_fishery_events STOPS if this file is absent (no silent fallback).
  run_fishery_spillover_diag = TRUE,
  razor_nearby_beaches = c("Twin Harbors", "Copalis", "Mocrocks"),
  fishery_opener_dates_file = "fishery_opener_dates.xlsx", # consolidated daily calendar; required.
  fishery_opener_sheet      = "data",

  # --- razor_dig SHORE-effort term (item 1, 2026-07-13) --------------------
  # Adds a razor-dig day-type effect to the SHORE effort model (a B3 * razor[d] term,
  # analogous to the holiday B2 effort term). "no" = off (production default); "yes" = on
  # for the shore fits; "auto" = on only if the spillover diagnostic's day-type/month-
  # adjusted shore-effort razor effect is significant at razor_dig_auto_p. Boat fits and
  # inactive runs pass razor = 0, so B3 stays decoupled (prior-only). Compare the shore-
  # effort elpd_loo against a "no" run to test the gain. RE-COMPILES the Stan model.
  razor_dig_mode   = "no",       # "no" | "yes" | "auto"  (Run 3 disqualified razor-dig B3: no elpd gain; keep off unless deliberately re-testing)
  razor_dig_auto_p = 0.05,     # auto-mode significance threshold (adjusted shore-effort p)

  # --- CPUE holiday + density terms (item 6, 2026-07-13) -------------------
  # B2_C (holiday CPUE effect, analogous to the effort B2) is ALWAYS on now (effort had
  # weekend + holiday terms, CPUE previously had only weekend). estimate_cpue_density adds
  # an OPTIONAL same-day-effort density-dependence term (gamma_C) to CPUE; it couples the
  # CPUE and effort processes, so it is OFF by default and should be validated on a test
  # fit first. RE-COMPILES the Stan model.
  estimate_cpue_density = FALSE,

  # --- Model-specific toggles (centralized here; each is read only by its own
  #     model and ignored by the other, so they are safe to keep in one list) --
  collapse_mu_hier           = FALSE, # (pooled) collapse the single-cell mu-hierarchy
                                       #   (B1.7/POOL-4 experiment lever). FALSE = current
                                       #   hierarchy, posterior unchanged. Accepts a per-
                                       #   population named list, e.g. list(private_boat = TRUE).
  estimate_B1_C              = TRUE,   # (gear-resolved) weekend/holiday CPUE effect B1_C.
                                       #   TRUE matches the pooled model; FALSE drops B1_C
                                       #   from the likelihood (v5.4 behavior).
  gear_resolved_G            = TRUE,  # (gear-resolved) GR-7 Phase 1. FALSE = production G = 1
                                       #   (gear split by PE apportionment). TRUE turns on genuine
                                       #   per-gear CPUE for SHORE fits (Option A1): only single-gear
                                       #   interviews feed a gear-specific CPUE, multi-gear trips form
                                       #   "Mixed", effort is split across gears by the pi_gear shares.
                                       #   Boat stays G = 1 (Pot-dominated; Phase 0). Changes shore
                                       #   inference, so validate by run. See
                                       #   07_documentation/development_notes/GR-7-per-gear-CPUE-design.md
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
