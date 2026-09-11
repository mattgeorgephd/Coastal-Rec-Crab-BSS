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
# run_config.R  --  the single control surface for a production run.
#
# ---------------------------------------------------------------------------
# RUNNING A NEW SEASON OR A NEW WINDOW? Start with 07_documentation/NEW_SEASON_GUIDE.md,
# which walks the whole workflow (naive run -> AR ladder -> read the diagnostics -> pin
# resolutions -> production run -> gear-track cross-check). The 2024-25 season was the
# DEVELOPMENT TEST SEASON; the architecture is built to run on any window the user
# selects (full season, part of a season, or a multi-season span via a vector
# season_filter). The short checklist of keys that MUST be revisited together:
#
#   the window       est_date_start, est_date_end
#   the data filter  season_filter (must match the season column of the workbooks;
#                    a stale value now stops the run with a plain message)
#   the calendar     pot_closure_start/end + pot_open_date; census_start_date/end;
#                    crabbing_holidays.xlsx rows; fishery_opener_dates.xlsx rows
#   the AR caps      ar_max_resolution (2024-25-derived; RE-DERIVE with the ladder,
#                    ar_escalate + ar_rung_adequacy, before trusting them on new data)
#   season-derived   every key tagged "SEASON-DERIVED" below (prior centers, floors,
#   priors/floors    the set crabbing fraction): revisit, do not assume
#
# Multi-season spans (2026-09-10, CHANGE_REGISTER A14): a span containing several pot
# closures is expressed with `pot_closures` (one entry per season) plus `census_windows`
# (a named per-season list); the report then adds season-level totals. See
# 07_documentation/NEW_SEASON_GUIDE.md section 7 for the four settings that must agree.
# ---------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------
# WHERE THIS WORK STANDS (context, 2026-09-04). Read this before treating any
# number in this repository as an estimate.
#
# NOTHING HAS BEEN PUBLISHED. WDFW has published no recreational Dungeness crab
# harvest estimate from this pipeline. The `main` branch holds the state of the
# model BEFORE two things happened: a meeting with the WDFW freshwater creel
# team, and confirmation from OSP that they can supply daily boat-count data.
# Everything on the `OSP-boat-count-incorporation` branch is work toward
# incorporating both. There is therefore NO published figure that a change here
# has to stay consistent with, and continuity with an earlier internal run is
# not by itself a reason to prefer one modelling choice over another.
#
# WHAT THE OSP DATA IS, and why the boat model is built the way it is. OSP will
# provide, per day: (a) the TOTAL number of vessels returning, and (b) the
# fraction of those that were CRABBING ONLY. (b) deliberately EXCLUDES combo
# trips that included crabbing alongside another fishery, so it is a LOWER BOUND
# on the vessels that did any crabbing, not the crabbing fraction itself. That
# is exactly what `osp_crab_lower` / `f_lower` exist for: the crab-only count
# bounds f from below while the combo-trip share stays on a prior. Do not wire
# the crab-only column in as if it were f.
#
# WHAT THE OSP DATA IS FOR. Two things, in order: improve the ACCURACY of the
# boat harvest estimate, and reduce its UNCERTAINTY. Before the OSP stream the
# boat rested on `tau_boat` ~ 1.2 from two I/E days, with catch proportional to
# it at elasticity 1.0; the dense OSP series identifies the turnover directly.
# The acknowledged limitation is that more BOAT INTERVIEWS are needed next
# season; no amount of boat counting fixes a thin CPUE sample.
#
# THE AR ESCALATION LADDER came out of the FW creel team meeting. The idea is
# that a component should report at the finest AR resolution its own sampler
# behaviour supports, decided per run rather than frozen into config. It is
# implemented as `ar_escalate` plus `ar_escalate_stop` / `ar_escalate_select`
# below, and it ships OFF so a production run is reproducible; turn it on
# deliberately.
# ---------------------------------------------------------------------------
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
  # 2026-09-08: renamed from "boat-count-validation-run", which had named every production
  # output folder since the OSP validation era and no longer described what a default run
  # is. Purely cosmetic: run_tag names the output folder, config_delta() ignores it, and no
  # fit or estimate depends on it. Historical folders keep their old names.
  run_tag           = "two-season-2023-25",

  # --- Identifiers ---------------------------------------------------------
  # These unify the two models onto one set of strings. The committed gear-
  # resolved driver used "Rec Crab Grays Harbor Westport 2024-25" while pooled
  # used the string below; centralizing here makes both use one value. If either
  # string must be preserved verbatim as an output identifier, give that model
  # its own value in its .Rmd instead.
  project_name      = "Coastal Recreational Crab",
  fishery_name      = "Rec Crab Grays Harbor 2024-25",

  # --- Season window (the values you change most often) --------------------
  # =========================================================================
  # ACTIVE RUN: TWO-SEASON SPAN, 2023-24 + 2024-25 (staged 2026-09-10).
  #
  # DATA PREREQUISITES (the 2026-09-10 rebuild of 04_input_files from the per-season
  # creel workbooks closed the first and third; the second cannot be closed):
  #   1. effort_combined.xlsx      NOW HOLDS 2022-23 through 2025-26 (2023-24: 4,116
  #      counts on 360 days; 20 Feb-2024 rows flagged and held out, see effort_qc_drop).
  #   2. wes_commercial_tally.xlsx has NO 2023-24 rows AND NONE CAN EXIST: no vessel tally
  #      was kept before 2024-25 and no charter roster either. The 2023-24 census component
  #      is therefore 0 and estimate_comm_charter() WARNS that the frame is missing (33
  #      Westport commercial-vessel interviews with no vessel count to expand to). A
  #      two-season port total is short by the 2023-24 commercial/charter catch.
  #   3. crabbing_holidays.xlsx    NOW HOLDS 2022-23 through 2026-27 (the 2024-25 named
  #      holidays applied by rule; build_crabbing_holidays.R).
  #   WBL_boat_counts.xlsx covers 2024-2025 only: no OSP for 2023-24 is SURVIVABLE
  #   (stream absent, boat trailer-only) but weakens the 2023-24 boat.
  #   interview_combined.xlsx: the 2023-24 gear count is present again (D15) and the
  #   2022-24 gear labels are harmonised to the 2024-25 vocabulary (a ring net was being
  #   classed as a trap by the gear-resolved regex); Grays Harbor 2023-24 shifts were
  #   shorter (median 4.1 h against 6.0 h in 2024-25), so the contact-based crabbing
  #   fraction covers less of the day in that season.
  #
  # SINGLE-SEASON ROLLBACK / 2025-26: swap the five values below for one of the
  # commented blocks that follow them.
  # =========================================================================
  est_date_start    = "2023-09-16",   # first day of the estimation window
  est_date_end      = "2025-09-15",   # last day
  # A vector runs a multi-season span (e.g. c("2024-25", "2025-26")) with a matching
  # est window; every listed season needs rows in the workbooks and a holiday calendar.
  # A stale value now stops the run loudly (validate_season_window.R).
  season_filter     = c("2023-24", "2024-25"),
  # -- rollback block (single-season 2024-25): --------------------------------
  #   est_date_start = "2024-09-16",  est_date_end = "2025-09-15",
  #   season_filter  = "2024-25",     pot_closures = NULL,
  #   census_windows = NULL,          run_tag = "production",
  # -- single-season 2025-26 (data through 2026-09-08; the season ends 2026-09-15): ----
  #   est_date_start = "2025-09-16",  est_date_end = "2026-09-15",
  #   season_filter  = "2025-26",     pot_closures = NULL,
  #   pot_closure_start = "2025-09-16", pot_closure_end = "2025-11-30", pot_open_date = "2025-12-01",
  #   census_windows = NULL,  census_start_date = "2025-12-01", census_end_date = "2026-01-03",
  #   commercial_opener = "2026-01-04", run_tag = "season-2025-26"
  #   (no OSP / WBL_boat_counts rows for 2025-26 yet: the boat stream runs on trailers + I/E;
  #   the tally has 23 days, the charter roster 11 Westport trips)
  # ---------------------------------------------------------------------------

  # --- Regulatory / structural dates ---------------------------------------
  # 2026-09-12 CORRECTION. The note that stood here said pot_open_date "is a SINGLE date
  # feeding the L_effective I/E regression split", a known multi-season approximation
  # (CHANGE_REGISTER A14). IT DID NOT. estimate_L_effective() took pot_open_date as an
  # argument and never referenced it; its only predictors are yday (quadratic) and day
  # type. The dead argument is gone and A14 is corrected. pot_open_date now has exactly
  # three consumers, all of them per-season and none of them a model term:
  #   - the "Pots open" vertical line on the season plots (pot_closure_markers);
  #   - the gear driver's monthly-panel vline;
  #   - build_subseasons()'s fallback for pot_closure_end when that is NULL.
  # SO A STALE VALUE IS COSMETIC, NOT SILENT-WRONG -- but it is stale here on any
  # single-season rollback: the value below belongs to 2023-24 because the shipped span
  # starts there. validate_season_window() now WARNS when it, or pot_closure_start/end,
  # falls outside the estimation window, which is what the 2026-09-11 ladder pin missed.
  # THE REAL multi-season approximation is that the L_effective regression pools every
  # season's I/E days into one yday -> L curve; it is documented at the function.
  pot_open_date     = "2023-12-01",   # pots legal from this date (plot markers, and the
                                      #   default pot_closure_end + 1 day)
  # Pot-closure window: the period when pots are NOT legal (only non-pot gear, ring
  # nets/snares/traps). Given explicitly here rather than assumed to start at the
  # season start, so a future season whose start does not coincide with the closure
  # start is supported. Outside this window pots are allowed (all-gear). If a closure
  # starts after est_date_start or ends before est_date_end, the driver adds the
  # corresponding all-gear period(s) automatically (see 03_R_functions/build_subseasons.R).
  # Keep pot_open_date = pot_closure_end + 1.
  # MULTI-SEASON CLOSURE CALENDAR (2026-09-10). Non-NULL outranks the scalar pair below:
  # one entry per pot closure in the span, each with its season label; the all-gear block
  # following a closure belongs to that season, and the report rolls fits up to season
  # totals (season_totals.csv + an HTML table on multi-season runs). Closure dates are
  # the same rule each season: season start through Nov 30, pots legal Dec 1.
  pot_closures = list(
    list(season = "2023-24", start = "2023-09-16", end = "2023-11-30"),
    list(season = "2024-25", start = "2024-09-16", end = "2024-11-30")
  ),
  pot_closure_start = "2024-09-16",   # scalar fallback, used only when pot_closures is NULL
  pot_closure_end   = "2024-11-30",   #   (here = season start .. day before pots open)
  commercial_opener = "2025-01-01",   # (was malformed "2025-01-1" in gear-resolved)
  # PER-SEASON CENSUS WINDOWS (2026-09-10). Non-NULL outranks the scalar pair below; a
  # named list season -> c(start, end), one commercial/charter window per season: pots
  # legal (Dec 1) through the last day before the coastal commercial fishery opened and
  # the vessels stopped fishing as recreational.
  #   2023-24  opened Feb 1, 2024 coastwide (WDFW bulletin, "Washington's coastal Dungeness
  #            crab commercial season opens Feb. 1"); no tally exists for this window.
  #   2024-25  the tally ran Dec 3 to Feb 8, 2025 (opener Feb 11).
  #   2025-26  Grays Harbor opened Jan 4, 2026 (WDFW news release 2025-12-29); the tally ran
  #            Dec 1 to Jan 3 with the gear-set days noted on the sheet.
  census_windows = list(
    "2023-24" = c("2023-12-01", "2024-01-31"),
    "2024-25" = c("2024-12-01", "2025-02-08")
  ),
  census_start_date = "2024-12-01",   # scalar fallback, used only when census_windows is NULL
  census_end_date   = "2025-02-08",
  # REVIEW ITEM 4 (2026-09-08, settled 2026-09-09): the commercial/charter tally is a census
  # on the days it was taken (47 of the 70 days of the 2024-25 window), and the samplers
  # are scheduled on the days the charter companies and the commercial vessels (fishing
  # as recreational) are confirmed to be operating, so the 23 days without a tally are
  # days with NO operation, not missed counts. census_expansion = "none" (shipped) makes
  # the component the exact sum over the tally days (7,884 crab on 2024-25) with the
  # unsampled days at zero; "day_type" is the 2026-09-08 expansion (unsampled days filled
  # with the sampled days' day-type mean, 11,753 on 2024-25, imputation SE 426), kept for
  # reproduction and for a season whose roster did not track operations.
  census_expansion = "none",
  # 2026-09-11 (Matt): the census-without-error statement applies to the COMMERCIAL boats.
  # The CHARTER vessels are not 100% sampled and need expansion, so the two populations now
  # have two estimators (03_R_functions/estimate_comm_charter.R):
  #   commercial  the tally is a census of the vessels; the component is the exact sum over
  #               the tally days. What is left is the per-vessel catch MEAN, a near-census
  #               sample mean (2024-25: 141 interviews on 164 vessel-trips, 86%, SE 123 =
  #               1.9%; 2025-26: 44 on 67, 66%, SE 146 = 5.5%), reported not carried.
  #   charter     an expansion over the charter TRIP frame (charter_frame below): trips x
  #               mean catch per interviewed trip, stratified by vessel, with the
  #               finite-population-corrected variance N^2 (1 - n/N) s^2 / n. 2024-25:
  #               34 trips, 20 interviewed (59%), 2,133 crab, SE 73 (3.4%).
  # census_uncertainty says which of those variances enters the port interval:
  #   "charter"  (shipped) the charter expansion only -- the commercial census is a constant
  #   "none"     neither (both reported in census_variance.csv)
  #   "sampling" all of it: the charter expansion, the commercial per-vessel mean, and the
  #              day-type imputation when census_expansion = "day_type" ("imputed_days" is
  #              still accepted as the 2026-09-08 name)
  census_uncertainty = "charter",
  # How the charter trips expand: "vessel" (shipped) strata by vessel, which matters because
  # the vessels differ (2024-25 Westport: Ultimate 70.2 crab/trip on 17 interviews, Outta
  # Line 30.7 on 3); "pooled" uses one charter mean for every trip (2,186, SE 118).
  charter_expansion = "vessel",

  # --- Catch groups --------------------------------------------------------
  estimate_red_rock = FALSE,          # TRUE adds Red_Rock_Kept alongside Dungeness

  # --- Day-type / stratification -------------------------------------------
  # WEEKEND DEFINITION (changed 2026-08-25, improvement 3). Friday was previously
  # typed as a weekend day; the 2024-25 data say it is a weekday. Within month, paired
  # across the 12 months carrying both day types:
  #   shore gear count  log(Fri/Sat) = -0.844 (Fri = 0.43x Sat), t = -8.88, p = 2.4e-6
  #                     log(Fri/Mon-Thu) = +0.190 (1.21x),        t =  1.81, p = 0.10
  #   boat trailers     log(Fri/Sat) = -0.507 (0.60x),            t = -3.81, p = 0.003
  #                     log(Fri/Mon-Thu) = -0.005 (0.99x),        t = -0.05, p = 0.96
  # So Friday is statistically indistinguishable from a weekday and 1.7-2.3x below
  # Saturday. Pooling it into "weekend" dragged the fitted weekend multiplier from
  # 2.34x to 1.74x (shore) and 1.81x to 1.42x (boat), and regressing log effort on
  # month + day-type shows moving Friday to the weekday stratum cuts the residual
  # variance the AR(1) must absorb by 15% (shore) / 5% (boat).
  # KNOWN RESIDUAL (documented, not fixed here): Sunday is sampled on only 6 of 50
  # shore days and 5 of 50 boat days, and the boat Sunday mean (3.4) sits BELOW the
  # weekday mean (6.9), so grouping Sunday with Saturday is not obviously right for
  # the boat stream. Splitting Sunday out needs a shift change in the field plan, not
  # a config edit; see PIPELINE_STATUS Section 6.
  days_wkend        = c("Saturday", "Sunday"),
  # --- Fishing-time filters (review item 8, 2026-09-08; unit-aware) ---------------
  # min_fishing_time applies ONLY when a population's effort unit is time-denominated
  # (shore on crabber-hours / gear-hours). Under the gear-deployment unit (production for
  # both populations) hours are not in the likelihood, so the threshold only discarded
  # data (153 shore and 44 boat rows on 2024-25, mostly incomplete trips but also 18
  # shore and 14 boat complete/unlabelled trips of 0 to 0.4 h). The one guard that
  # remains under deployments is drop_unfished_zero_catch: a row with NO positive time in
  # hours_fished, crabber_hours or gear_hours AND zero catch is gear set and not yet
  # fished (pots still soaking), not a fished deployment, whatever completed_trip says.
  # Measured effect (reader, 2024-25): shore 3,597 -> 3,739 rows, complete-trip CPUE
  # 0.979 -> 0.976; boat 162 -> 184 rows, 3.26 -> 3.17; commercial/charter unchanged.
  # See apply_fishing_time_filters() in 03_R_functions/fetch_crab_data.R.
  min_fishing_time  = 0.5,            # min crabber-hours to keep an interview (TIME units only)
  drop_unfished_zero_catch = TRUE,    # deployment units: drop no-time AND zero-catch rows
  period_pe         = "week",         # PE temporal stratum
  sections          = c(1),
  # =========================================================================
  # THE PE's THREE UNSAMPLED-CELL LEVERS (rebuilt 2026-09-12; one shared
  # implementation in 03_R_functions/pe_effort_strata.R, used by both PE runners)
  # -------------------------------------------------------------------------
  # The PE stratifies by (week x day_type) from the SAMPLED days and left-joins the
  # calendar. Roughly half of every component's calendar days therefore sit in a cell
  # that is either UNSAMPLED (no sampled day at all) or SINGLETON (exactly one), and the
  # three levers below say what happens to each. On 2024-25 that is not a corner case:
  #
  #   component            cells  unsampled        singleton        effort with no variance
  #   shore all-gear         93   21 (45 days)     37 (82 days)     54% (SE read 1.5%)
  #   boat all-gear          93   23 (48 days)     34 (77 days)     49% (SE read 5.5%)
  #   shore pot-closure      24    1 ( 1 day )     10 (19 days)     41%
  #   boat pot-closure       24    6 (11 days)      7 (14 days)     29%
  #
  # WHY THIS IS NOT PE-ONLY, contrary to the note that stood here until 2026-09-12. A
  # component whose convergence gate FAILS reports its PE point in the port total, as a
  # CONSTANT with no interval (see section 7.7 of either driver). The boat pot closure is
  # the component most likely to fall back. So these levers can reach the headline.
  # =========================================================================
  #
  # 1. THE CPUE OF A CELL WITH NO SURVIVING INTERVIEWS (item 2, 2026-07-13; extended
  #    2026-09-12). A cell can carry expanded effort and still have no interview left
  #    after the incomplete-trip filter.
  #      "local"  (SHIPPED 2026-09-12) that MONTH's ratio-of-sums, falling back to the
  #               sub-season's. Matches the scale of the effort fill below.
  #      "pooled" the sub-season-wide ratio-of-sums (the behaviour through 2026-09-11).
  #      "zero"   assign zero catch (the pre-2026-07-13 behaviour; reintroduces the
  #               sparse-stratum sign instability in the thin boat PE).
  #    MEASURED, and it is the LARGER lever on the boat, which the 2026-09-11 note missed:
  #    with the month-local effort fill below, boat all-gear PE catch is 42,841 under
  #    "pooled" and 37,018 under "local", -13.6%. The imputed boat cells are in the SUMMER
  #    months and the sub-season-pooled boat CPUE is dominated by the high-CPUE winter, so
  #    "pooled" was multiplying summer effort by a winter-inflated rate. The shore is
  #    almost unmoved (29,705 -> 29,737), so this is a boat correction. The shipped PE
  #    port total is 85,076, +17.8% on the historical 72,224, NOT the +26% reported on
  #    2026-09-11: that figure mixed a month-local effort fill with a season-pooled CPUE
  #    and so double-counted the seasonal gradient. The full decomposition, 2024-25:
  #      72,224  effort "zero",           CPUE "pooled"   (through 2026-09-11)
  #      90,861  effort "local_day_type", CPUE "pooled"   (+18,637 from the effort fill)
  #      85,076  effort "local_day_type", CPUE "local"    (-5,785 from the CPUE fill)
  pe_empty_stratum  = "local",        # "local" | "pooled" | "zero"
  #
  # 2. THE EFFORT OF A CELL WITH NO SAMPLED DAY (2026-08-25; default changed 2026-09-12).
  #      "local_day_type" (SHIPPED 2026-09-12, Matt) the same day type in the same MONTH,
  #               falling back to the same day type sub-season-wide and then to the
  #               sub-season mean. Local, so a February cell is filled at February rates
  #               rather than at a mean the summer dominates.
  #      "day_type" that sub-season's day-type mean. A sub-season spans a 20-fold seasonal
  #               swing, so this overfills a thin winter week and underfills a summer one.
  #      "zero"   the pre-2026-09-12 default: the cell contributes NOTHING.
  #    WHY "zero" IS NO LONGER THE DEFAULT. It is not a neutral choice; it is an
  #    assumption that 15.6% of the shore's days and 16.6% of the boat's had no fishing,
  #    on days that are disproportionately weekend and holiday days carrying 1.7-2.3x
  #    weekday effort. And its error is a BIAS, which no SE can carry: under "zero" the
  #    component's interval is silent about the omission, while under a fill the
  #    imputation is priced (lever 3). Every run now reports pe_zeroed_effort_bias, the
  #    effort the zeroed cells WOULD carry, so choosing "zero" states its own cost.
  #    2024-25 PE port total at the shipped CPUE fill: 72,224 ("zero"), 81,160
  #    ("day_type"), 85,076 ("local_day_type", shipped). The effort totals are 45,175 /
  #    53,846 / 56,793: "zero" omits 20% of the season's effort.
  #    Still worth reading against the post-ladder BSS rather than the stale
  #    20260904/pooled-CPUE-AD-A1-adopted reference (CHANGE_REGISTER D19); the ladder's R0
  #    desk rung prints all three so the comparison costs no MCMC.
  pe_empty_effort_stratum = "local_day_type",  # "local_day_type" | "day_type" | "zero"
  #
  # 3. WHAT AN UNSAMPLED OR SINGLETON CELL CONTRIBUTES TO THE SE (NEW 2026-09-12,
  #    CHANGE_REGISTER D21). The SE was sqrt(N^2 sd^2 / max(n,1)) with sd from the cell's
  #    own sampled days, and sd() of ONE observation is NA, replaced by 0. So a singleton
  #    cell contributed its full point estimate and no variance, and an IMPUTED cell did
  #    too: the shore all-gear effort SE was bit-identical (410) whether the fill added
  #    0, 6,240 or 7,947 effort units. Half of every component's reported effort had no
  #    uncertainty attached to it at all.
  #      "impute_aware" (SHIPPED) a singleton cell borrows its donor's spread with the
  #               divisor still 1 (the collapsed-stratum estimator standard for
  #               one-per-stratum designs: Hansen, Hurwitz & Madow 1953; Cochran 1977
  #               sec. 5A.12; Wolter 2007 ch. 2); an imputed cell carries the donor
  #               mean's own sampling variance PLUS a between-cell surrogate of one
  #               cell's worth of the donor's day-to-day spread, s^2 (1/n_donor + 1).
  #               Deliberately conservative: an imputed number should read wider than a
  #               measured one.
  #      "donor_mean_only" drops the between-cell term (s^2 / n_donor). The lower bound.
  #      "sampled_only"   the pre-2026-09-12 arithmetic, exactly. Kept so any earlier
  #               number stays reachable; every run reports effort_se_sampled_only
  #               alongside effort_se whatever this is set to.
  #    MEASURED on 2024-25 at the shipped fill, and the POINT ESTIMATE DOES NOT MOVE (all
  #    three settings give the same 85,076), so this lever is pure honesty about the
  #    interval: shore all-gear effort SE 410 -> 1,379 (1.5% -> 3.9% of the component),
  #    boat all-gear 573 -> 1,310 (5.5% -> 9.4%), shore pot-closure 238 -> 550, boat
  #    pot-closure 29 -> 79. "donor_mean_only" lands between (1,184 and 1,103 on the two
  #    all-gear components). Under the retired "zero" fill the same lever gives 1,099 and
  #    1,046, i.e. MOST of the understatement was never about imputation at all: it was
  #    the 37 and 34 singleton cells. The mean and the spread are
  #    borrowed from DIFFERENT donor levels on purpose (finest with any sampled day for
  #    the mean, finest with two or more for the spread), so this lever and lever 2 move
  #    independent things and the ladder can attribute them one at a time.
  #    NOT INCLUDED, and tracked separately as D22: the PE applies no finite-population
  #    correction anywhere, and the PE catch still has NO SE at all (only effort does).
  pe_variance       = "impute_aware", # "impute_aware" | "donor_mean_only" | "sampled_only"

  # --- Incomplete-trip filter (both models) --------------------------------
  # Incomplete trips (soak-time gear not yet retrieved) read systematically low
  # (about -20% CPUE for pots/traps), biasing CPUE and hence the harvest estimate
  # low. TRUE excludes them from CPUE estimation (PE and BSS), keeping Complete +
  # NA; FALSE keeps all trips (pre-filter behavior). Missing trip_status is kept.
  filter_incomplete_trips = TRUE,

  # --- Incomplete-trip TREATMENT diagnostic (improvement 5, 2026-08-25) -----
  # DIAGNOSTIC ONLY. The production estimator is unchanged (exclude, above); this
  # runs the fast design-based PE under four treatments of the incomplete trips and
  # reports what each would do, so the choice can be made on evidence rather than
  # argument. Output: sensitivity_incomplete_trips.csv + an on-page table.
  #
  #   exclude            drop them entirely from BOTH the CPUE estimate and the
  #                      gear-ratio estimates (R_G shore / gear-per-group boat).
  #                      This is what the pipeline does today.
  #   gear_only          drop their CATCH from the CPUE estimate but KEEP their gear
  #                      counts in the gear ratios. An interrupted trip's gear count
  #                      is fully observed; only its catch is truncated, so this is
  #                      the arm that answers "is option 2 already how the model
  #                      functions?" -- it is NOT: prep_bss_crab_*() derives intA
  #                      (the R_G / R_G_boat set) from the ALREADY-FILTERED interview
  #                      frame, so today an incomplete trip's gear count is discarded.
  #                      Note the PE and BSS currently DISAGREE on this: the boat PE's
  #                      gear_per_group is computed from the unfiltered interview set,
  #                      while the boat BSS learns R_G_boat from complete trips only.
  #   impute_mean_cpue   keep them, replacing their observed catch with (complete-trip
  #                      ratio-of-sums CPUE x their gear count). Expect ~0 movement:
  #                      at the pooled level a ratio-of-sums over the union returns the
  #                      complete-trip CPUE identically, so any movement is purely the
  #                      re-weighting of thin week x day-type strata toward the pooled
  #                      rate. Reported so that "no-op" is a measured result, not a claim.
  #   keep               no filter at all (the pre-v7.5 behavior).
  #
  # The diagnostic also reports the gear ratio (R_G shore / gear-per-group boat) under
  # each arm, since that -- not the CPUE -- is where the recovered interviews actually
  # bite, and a length-biased-sampling check (are intercepted trips gear-heavier?).
  diagnose_incomplete_trips = TRUE,
  incomplete_trip_arms      = c("exclude", "gear_only", "impute_mean_cpue", "keep"),

  # --- Effort unit (both models) -------------------------------------------
  # Shore and boat CPUE denominators. As of pooled v7.7 / gear-resolved v5.5 both
  # components run on gear-DEPLOYMENTS: the pipeline's own linearity diagnostic
  # flags every time-denominated unit as invalid for pots (shore beta_h 0.57 for
  # crabber-hours, 0.73 for gear-hours, 1.05 for deployments; deployments is the
  # only harvest-unbiased unit). Routed through 03_R_functions/bss_effort_spec.R
  # so the BSS and PE always share a unit. Set shore_effort_unit = "crabber-hours"
  # to revert shore only.
  shore_effort_unit      = "gear-deployments",  # "crabber-hours" | "gear-hours" | "gear-deployments"
  # --- SHORE TURNOVER PRIOR (review item 2, 2026-09-08) ------------------------
  # 1.7 is arrivals / PEAK presence over the WDF20 I/E days (SEASON-DERIVED, 2024-25).
  # The Stan effort likelihood calibrates lambda_E to the gear counts AT THE HOURS THEY
  # WERE TAKEN, and 81% of 2024-25 counts fall between 10:00 and 13:59, where presence
  # averages 0.69 of the daily peak (40 I/E days with the `time` column added
  # 2026-09-08). The multiplier those counts need is arrivals / presence-at-count-time,
  # about 2.5, and the shore component is proportional to it.
  # "derived" resolves the centre from THIS run's I/E time column and count hours
  # (03_R_functions/bss_day_length.R, estimate_shore_turnover; written to
  # shore_turnover_*.csv EVERY run whatever this key says), sets the log-SD to the
  # bootstrap log-SE (floored at tau_shore_prior_sigma_floor), and waives the shore
  # shared-turnover floor so ONE shared shore level carries the level uncertainty.
  # ADOPTED 2026-09-09 (Matt): the count hours were not in the I/E workbook when 1.7 was
  # derived, so the old turnover was the wrong quantity, not a competing estimate; both
  # keys ship "derived". Rung R4 of 06_diagnostics/run_improvements_2026-09-08.R still
  # isolates the change (R3 pins 1.7) so its size on the shore is on the record. To hold
  # the old value: set tau_shore_prior_mu = 1.7 and tau_shore_prior_sigma = 0.3.
  tau_shore_prior_mu     = "derived", # shore deployment turnover from the I/E time column (2.48 on 2024-25)
  tau_shore_prior_sigma  = "derived", # bootstrap log-SE of the level, floored below
  tau_shore_prior_sigma_floor = 0.10, # never claim better than 10% level precision from the prior
  tau_shore_prior_mu_fallback = 1.7,  # used by "derived" when the I/E time column is absent
  tau_shore_derive_min_days   = 10,   # minimum I/E days behind a derived centre
  # WHICH I/E DAYS THE DERIVED CENTRE USES (surfaced 2026-09-11 by the first full ladder
  # run; CHANGE_REGISTER D24). FALSE (shipped) pools EVERY I/E interval day in the
  # workbook. On the 2024-25 run that is 40 days spanning 2023-08 to 2026-08, of which 7
  # are inside the season and 4 inside the all-gear sub-season the resulting 2.477 scales
  # by 1.36 -- so the shipped value is a multi-season quantity, and it is the SAME number
  # whichever season you run. That was not documented and NEW_SEASON_GUIDE described it as
  # derived from the window; the boat side ("calibration") really is window-filtered, this
  # one never was.
  #
  # It stays pooled by default because 7 days is a thin basis for the second-largest mover
  # in the series, and because the derivation needs a diel profile (presence by hour) that
  # a handful of days estimates badly. TRUE restricts it to the estimation window so the
  # alternative can be priced. Every run now reports n_days_in_window, n_days_total and
  # tau_in_window in shore_turnover_summary.csv either way.
  #
  # THE IN-WINDOW CHECK THAT EXISTS TODAY is sigma_IE in the fitted output: the shore
  # likelihood is arrivals ~ lognormal(log(lambda_E * tau_shore), sigma_IE), so if the
  # pooled turnover over-expands this window's counts, sigma_IE grows. On the 2024-25 run
  # it went 0.373 -> 0.577 (+55%) when the derived centre was adopted, on FOUR in-window
  # days -- a direction, not a measurement. A paired gear count on every I/E day is the
  # field fix; see PIPELINE_STATUS section 1v.
  tau_shore_derive_window_only = FALSE,  # FALSE (shipped) = pool all I/E days | TRUE = window only
  # --- BOAT TURNOVER PRIOR (review item 3, 2026-09-08) --------------------------
  # "calibration" resolves the prior centre from THIS window's OSP/trailer overlap
  # (03_R_functions/bss_turnover_prior.R): the implied turnover of the
  # tau_boat_calibration_metric row of osp_trailer_overlap_calibration.csv, i.e. the
  # OSP daily boat total divided by the trailer snapshot on the paired days. The
  # mean-per-visit metric is the consistent choice because the Stan likelihood treats
  # every trailer count as an observation of lambda_E / R_G_boat (so the model's trailer
  # level is the per-visit mean) and the PE expands mean_count per day the same way.
  # The resolved number REPLACES this key at run time, so the PE (review item 5), the
  # shared-turnover prior and the sensitivity projection all read one value.
  # WHY. The previous centre, 1.2 from two WBL I/E days, was shrinking the adopted
  # tau_bar (2.60 [2.06, 3.25]) about 13% below its likelihood-only value (~2.97),
  # and the PE still expanded on 1.2, so the boat PE-vs-BSS gap was mostly two turnovers
  # disagreeing. A number here (e.g. 1.2) restores the historical behaviour exactly.
  # The prior is deliberately WIDE (log-SD 0.5): the overlap days are in the likelihood
  # too, and at 0.5 the prior carries ~6% of the posterior precision, so the double use
  # is harmless. Do not tighten below 0.3 without removing the overlap from one side.
  tau_boat_prior_mu      = "calibration",
  tau_boat_prior_sigma   = 0.5,
  tau_boat_calibration_metric = "trailer_mean_per_visit",  # | "trailer_max_per_day" | "trailer_sum_per_day"
  # SEASON-DERIVED fallback for a window with no OSP/trailer overlap (fewer than
  # tau_boat_calibration_min_pairs paired days): the 2024-25 max-per-day calibration.
  tau_boat_prior_mu_fallback  = 2.7,
  tau_boat_calibration_min_pairs = 3,

  # --- SHARED TURNOVER (improvement 2.1, 2026-08-27) -------------------------
  # FALSE (default) keeps the historical parameterization: L is D INDEPENDENT per-day draws,
  # each anchored on tau_*_prior_mu, with nothing pooling information across days. The
  # 2026-08-26 ladder showed what that costs. Because no shared parameter exists, an
  # observation stream covering a SUBSET of days cannot move the season-level turnover:
  # 4 shore I/E days out of 289 leave the median L at 1.6998 against a prior centre of
  # 1.7000, and 148 OSP days out of 289 leave the boat at 1.201 against 1.200 -- while the
  # OSP/trailer overlap calibration puts the real turnover at 2.0-3.0 and the free kappa_OSP
  # sat at 3.15. The ~2.5x conflict goes into the OSP overdispersion instead (r_OSP ~ 1.6),
  # and shows up independently as a boat trailer PIT mean of 0.42 against a nominal 0.50.
  #
  # TRUE replaces the D anchors with ONE estimated tau_bar and a fixed day-to-day spread, so
  # every observed day informs a common level. It is REFUSED (with a warning) whenever L is
  # not a constant turnover, e.g. under a time-denominated shore unit.
  #
  # THIS MOVES THE BOAT NUMBER. osp_scale_is_tau roughly doubles the private-boat harvest,
  # and this decides whether the size of that doubling comes from the data or from the prior.
  #
  # ADOPTED 2026-09-01 (was FALSE) on the evidence of the 2026-08-30 Stage 5 batch; the full
  # argument is in development_notes/stage5-batch-review-2026-08-31.md. Five independent
  # lines agree and none dissent:
  #   1. It improves the stream it should and leaves the other alone: trailer elpd -573.8 ->
  #      -559.1 (+14.7 nats) for 3.1 effective parameters, catch elpd -451.4 -> -451.5.
  #   2. It adds no unreliable LOO points (Pareto k > 0.7 stays at 0) and does not raise
  #      p_loo as a fraction of n_obs (8.0% either way).
  #   3. Every calibration statistic moves toward nominal: OSP PIT mean 0.579 -> 0.520
  #      (nominal 0.500), PIT sd 0.277 -> 0.285 (nominal 0.289), coverage_50 unchanged at
  #      0.508 (nominal 0.500); trailer PIT mean 0.424 -> 0.484.
  #   4. It replicates across the two independently parameterized CPUE tracks to 0.03% on
  #      tau_bar itself (pooled 2.5969, gear-resolved 2.5962) and 0.80% on the component.
  #   5. The posterior 2.597 [2.064, 3.249] straddles the INDEPENDENT OSP/trailer overlap
  #      calibration of 2.01-3.03, which is computed from different data (n = 61).
  # The contrary case that was tested and REJECTED is a daily boat AR, which reaches a
  # similar-looking gain by absorption: 27-48 effective parameters on 195 observations, the
  # catch stream 3.8 nats WORSE, 2 then 16 Pareto k above 0.7, sigma_r_OSP collapsing 0.806
  # -> 0.086, and OSP coverage_50 at 0.977 against a nominal 0.500. Do not adopt it.
  #
  # THE FEATURE IS BOAT-ONLY IN EFFECT, BY THE FLOOR BELOW, NOT BY A POPULATION SWITCH.
  # A global toggle moved the SHORE all-gear component +17.9% on 4 informed days out of 289,
  # with an interval containing its own prior centre and no replication in the gear track.
  shared_tau             = TRUE,
  # Fixed day-to-day log-scale spread of L around tau_bar. PINNED at the value every
  # adopted run used (0.15 = half of the old 0.3 prior SD) on 2026-09-08, because the boat
  # prior SD moved to 0.5 (above) and the NULL default (half the prior SD) would have
  # widened the day-to-day spread to 0.25 as a silent side effect of re-centring the prior.
  # May be a per-population list; the shore's between-day spread of the count-hour
  # turnover is ~0.37 (log-SD, 40 I/E days), the boat's 0.15 is the value every adopted
  # run used. Kept as a scalar (boat-only in effect) until the derived shore prior is
  # adopted, when list(shore = 0.35, private_boat = 0.15) is the recommended setting.
  shared_tau_sigma       = 0.15,
  # Minimum days that can INFORM L (I/E days, plus OSP days when osp_scale_is_tau = 1) before
  # a fit is allowed a shared level; below it the fit degrades to per-day draws with a printed
  # reason. Stated explicitly here rather than left to the helper default, because it is the
  # thing that makes the line above boat-only. Observed counts: shore all-gear 4, shore pot
  # closure 0, boat all-gear 130, boat pot closure 18. Any threshold in 5..18 separates the
  # shore from the boat; 15 sits just below 18 deliberately, to keep the boat pot closure,
  # whose tau_bar (1.873 pooled, 2.050 gear-resolved) is corroborated across tracks. Setting
  # 20 drops that component out of the feature and is worth running once as a sensitivity.
  # SEASON-DERIVED floor: 15 makes the shared turnover boat-only on 2024-25 (130
  # OSP-informed days vs 0 shore). A season with sparse OSP coverage can drop the boat
  # below this floor; the run prints the informed-day count next to this decision.
  # 2026-09-08: may also be a per-population list, list(shore = 0, private_boat = 15);
  # bss_resolve_tau_shore_prior() sets the shore entry to 0 when the shore prior is derived.
  shared_tau_min_obs     = 15,

  gear_per_group_default = 4.0,       # PE fallback gear-per-boat-group when no interview records it

  # --- Shore I/E OBSERVATION unit (improvement 1/2 fix, 2026-08-25) ---------
  # WHAT WAS WRONG. The shore I/E likelihood is  IE_obs ~ lognormal(log(lambda_E * L),
  # sigma_IE).  When shore ran on crabber-hours that was right: lambda_E is crabbers
  # present and L was L_effective in HOURS, so lambda_E * L is crabber-hours, which is
  # what ie_crabber_hours holds.  The v7.7 move to gear-deployments replaced L with
  # tau_shore (a dimensionless turnover, ~1.7) but left the OBSERVATION as crabber-hours,
  # so the model has since been comparing crabber-HOURS against a predicted
  # lambda_E * tau = crabber-TRIPS.  On a typical WDF20 survey day that is ~331 observed
  # against ~80 predicted, a ~4x scale mismatch that the gear-count stream then has to
  # fight; it is the most likely mechanism behind the unexplained shore all-gear
  # sigma_IE ~ 1.07 (backlog GR-9).  The BOAT stream was already fixed this way in F2
  # (observation = boat trips, predicted = groups x tau); this puts shore on the same
  # footing.
  #
  # "auto" (default, CORRECTED): the observation follows shore_effort_unit --
  #     gear-deployments -> IE_obs = crabber ARRIVALS (ie_trips), predicted lambda_E * tau_shore
  #     crabber-hours / gear-hours -> IE_obs = ie_crabber_hours,  predicted lambda_E * L_eff
  # "crabber_hours": force the historical (mismatched under deployments) behavior, for
  #     reproducing pre-2026-08-25 runs exactly.
  # THIS MOVES THE SHORE NUMBER. It is the headline inference change in this patch and
  # must be confirmed by a run against Run 6 / the OSP validation baseline.
  ie_shore_obs_unit      = "auto",    # "auto" | "crabber_hours"

  # --- OSP second effort stream (Phase 1; boat only) -----------------------
  # Adds the OSP daily boat-total series (fetch_osp_boat_counts) as a SECOND boat
  # effort observation on the same latent lambda_E, scaled by the OSP within-day
  # turnover kappa_OSP. Both tracks (pooled + gear-resolved) carry the stream. FALSE
  # reproduces the pre-Phase-1 boat EXACTLY. PRODUCTION DEFAULT = TRUE (adopted
  # 2026-07-31, validation batch): OSP contributes on its 148 operating days and the
  # model degrades to trailer-only outside them. Changes the boat effort posterior,
  # so validate by run.
  use_osp_boat_counts   = TRUE,
  # SEASON-DERIVED from the 2024-25 OSP/trailer overlap days:
  osp_scale_prior_mu    = 3.0,        # kappa_OSP prior center = OSP/trailer overlap ratio
                                      #   (1 / mean-per-visit origin slope ~0.33 -> ~3.0;
                                      #   validation posterior kappa_OSP = 3.15, 95% [2.50, 3.91])
  osp_scale_prior_sigma = 0.3,
  osp_effort_col        = "WestportPrivateEffort",  # OSP daily ALL-boat total column
  # OSP CRAB-ONLY column (improvement 8, 2026-08-25). Optional extra column in
  # WBL_boat_counts.xlsx holding the count of boats OSP labelled as crabbing ONLY.
  # Absent today (the field request is outstanding) -> every OSP feature below is inert
  # and f behaves exactly as it does now. See the use_osp_crab_lower block for the
  # combo-trip problem this column does and does not solve.
  osp_crab_only_col     = "WestportCrabOnlyEffort",

  # --- Crabbing fraction f (Phase 2; boat only) ----------------------------
  # f = share of private boats at the launch that are crabbing (vs other fisheries).
  # Trailer and OSP counts are ALL boats, so f converts all-boat effort to crab effort;
  # without f the model implicitly assumes f = 1 (biases the boat catch high). f is
  # applied in the BSS boat generated quantities ONLY (scales effort + catch), decoupled
  # from sampling. Hybrid: a Beta prior centered on the SET VALUE below, updated by the
  # crab-creel I/E crab-vs-total classification once it lands (columns ie_crab_col /
  # ie_total_col in ingress_egress.xlsx); until then f uses the set value.
  # NOTE: 0.3 is a STARTING set value chosen by Matt, NOT a measured estimate; it scales
  # the boat catch and effort to ~30% of the f = 1 value (validation Step 3/Step 4 confirm
  # the boat scales exactly linearly in f, CPUE invariant). Replace with the WBL egress
  # pilot's f_hat when it lands. The PE is f-adjusted too (Phase 2b), so the boat PE and
  # BSS both carry f. Set use_crab_fraction = FALSE to reproduce the pre-Phase-2 boat
  # (f = 1). PRODUCTION DEFAULT = TRUE (adopted 2026-07-31).
  # 2026-09-08 (review item 1): the classification data now EXIST (the sampler contacts,
  # crab_fraction_source) and f is a month-by-month random walk observed by them
  # (crab_fraction_dynamic, below); the set value is only the level prior's centre.
  use_crab_fraction         = TRUE,
  # SEASON-DERIVED placeholder (no supporting observations; see CHANGE_REGISTER D2):
  crab_fraction_set         = 0.3,    # set value = Beta prior mean and the thin-data fallback
  crab_fraction_prior_kappa = 20,     # Beta concentration (prior SD ~0.10 at mean 0.3)
  crab_fraction_fixed       = NA,     # a number PINS f exactly (no uncertainty; sensitivity)
  crab_fraction_min_obs     = 20,     # min classified boats before the I/E Binomial updates f
                                      # (legacy construction; under crab_fraction_dynamic it is the
                                      #  PE's interpolation trigger only, every day enters the BSS)
  ie_crab_col               = "boats_crabbing",  # ingress_egress column: crab-classified boats
  ie_total_col              = "boats_total",     # ingress_egress column: total classified boats

  # --- Crabbing-fraction stratification + OSP-tau (Phase 3) -----------------
  # crab_fraction_strata: "none" (one scalar f = Phase 2) | "month" | "day_type" |
  # "month_day_type". Per-stratum f (a Beta prior per stratum, updated by that stratum's
  # I/E classification counts). "none" reproduces Phase 2 exactly. Move to "month" once the
  # pilot classifies enough boats per month; a single annual f over-states the summer crab
  # share (summer is tuna/salmon-dominated), which is where the harvest is largest.
  # REVIEW ITEM 1 (2026-09-08): "month" (year-qualified, chronological) is now the shipped
  # stratification, because the classification data exist (crab_fraction_source below).
  crab_fraction_strata      = "month",
  # Where the classification rows come from: "interviews" = the sampler contacts (every
  # private boat approached at the launch, crabbing or not; 422 on 2024-25, share 0.49
  # pooled, 0.9+ in Dec-Feb, 0.26 in Jul-Sep; a combo trip counts as CRABBING), "ie" = the
  # WPT/WBL egress classification columns (blank until the pilot delivers), "both".
  # CAVEAT carried, not hidden: contacts happen during sampler shifts, so boats returning
  # outside the shift are not classified; if finfish boats return later than crab boats
  # the contact share overstates f, most in summer. OSP's all-day crabbing-only column,
  # when it arrives, is the check (crab_fraction_source_rows in crab_fraction.R).
  crab_fraction_source      = "both",
  # Improvement 4/8 interaction: "opener" and "month_opener" stratify f by whether a
  # competing fishery was open that day (params$opener_f_dates, built by the driver from
  # opener_f_flag below). This exists because an opener EFFORT covariate and a constant f
  # actively fight each other: an MA2-halibut term makes lambda_E track the ~+32-trailer
  # halibut surge more faithfully, and a constant f then converts that surge into crab
  # effort at the same 30% rate as a closed-halibut day, which is exactly the boats we
  # know are NOT crabbing. If you switch on a boat opener covariate (below), pair it with
  # an opener-aware f or you make the bias worse, not better.
  opener_f_flag             = "ma2_halibut_open",   # opener used by the "opener" f strata

  # --- OSP crab-only counts as a LOWER BOUND on f (improvement 8, 2026-08-25) ---
  # THE PROBLEM. OSP can report, per day, the number of boats that were ONLY crabbing.
  # It does not record combo trips: a boat that crabs AND fishes another fishery is
  # labelled by the non-crab fishery. So the OSP crab-only share is a LOWER BOUND on the
  # true crabbing fraction f, not an estimate of it, and OSP alone cannot separate
  # "few boats, all crabbing" from "many boats, some also crabbing".
  #
  # THE PARAMETERIZATION (hard lower bound). Per f stratum k:
  #     f_lower[k] ~ Beta(1,1)   updated by  osp_crab_only[k] ~ Binomial(osp_total[k], f_lower[k])
  #     f[k]       = f_lower[k] + (1 - f_lower[k]) * theta[k]
  # theta[k] in [0,1] is the share of the NOT-crab-labelled boats that were also crabbing
  # (the combo trips OSP cannot see). f can therefore never fall below what OSP directly
  # observed, which is the whole point of the bound.
  #
  # HOW IT DEGRADES (this is the "defaults back to interview data" behaviour):
  #   OSP + egress classification in a stratum -> f_lower from OSP, theta pinned by the
  #       egress Binomial n_crab ~ Binomial(n_total, f); both identified.
  #   OSP only                                 -> f_lower from OSP, theta on its combo
  #       prior; f is bounded below by data and the rest rests on the prior.
  #   Egress only / neither (winter, OSP dark) -> f_lower is PINNED TO 0 and theta's prior
  #       is set to the ordinary f prior Beta(set*kappa, (1-set)*kappa), so f reduces
  #       EXACTLY to today's Phase 2/3 behaviour. Backward compatibility is exact.
  # f still enters generated quantities only, so the boat stays exactly linear in f and
  # the model CPUE stays invariant (the validated Phase 2/3 property is preserved).
  # SHIPS OFF, deliberately. Every other behaviour-changing feature in the 2026-08-25
  # batch is opt-in, and this one has to be too: the moment WBL_boat_counts.xlsx gains a
  # crab-only column, flipping this to TRUE would change the boat harvest on the SAME run
  # that first reads the new data, with no baseline to compare against. Worse, at
  # crab_fraction_strata = "none" the whole window collapses to one stratum and f becomes
  # p_osp + (1 - p_osp) * crab_fraction_combo_share, so a large part of the move would come
  # from the combo-share PLACEHOLDER rather than from OSP's data. Turn it on deliberately,
  # in its own run, against a run with it off. Requires use_osp_boat_counts = TRUE.
  use_osp_crab_lower        = FALSE,
  crab_fraction_osp_min_obs = 20,     # min OSP-classified boats in a stratum before the bound binds
  # Beta-binomial concentration for the DAILY OSP crab-only shares. The bound is fitted
  # one observation PER DAY, not on the stratum sum: summing ~150 operating days at ~50
  # boats would give a Binomial n of ~7,500 and an f_lower SD near 0.005, which ~150
  # correlated daily observations do not support. kappa is estimated from the spread of
  # the daily shares, so the bound's precision is measured rather than assumed. 20 centres
  # it on a daily share SD of about 0.10 at a mean of 0.30; large kappa -> binomial.
  crab_fraction_osp_kappa_prior_mu = 20,
  # Aggregate BOTH f streams (egress classification and OSP crab-only) over the fit's own
  # day set, so the stratum theta is pinned on and the stratum f_lower is bounded by
  # describe the same dates. FALSE restores the older behaviour where the egress rows were
  # pooled across the whole workbook. Behaviour-neutral today (the egress columns are blank).
  crab_fraction_restrict_to_fit = TRUE,
  # theta's prior in an OSP-INFORMED stratum: the share of non-crab-labelled boats that
  # were also crabbing. 0.15 is a PLACEHOLDER with the same standing f = 0.3 had before
  # the pilot: it is a guess, not a measurement, and it is what turns the OSP lower bound
  # into a level. It is pinned by the WPT egress classification wherever both streams
  # cover the same stratum, which is the strongest argument for running the egress pilot
  # on days OSP is also in port.
  crab_fraction_combo_share = 0.15,
  crab_fraction_combo_kappa = 8,      # Beta concentration (prior SD ~0.12 at mean 0.15)

  # --- The DYNAMIC crabbing fraction (review item 1B, 2026-09-08) --------------------
  # TRUE replaces the per-stratum Beta/Binomial f above with a random walk on the log-odds
  # of f across the month strata in chronological order, observed PER DAY by the sampler
  # boat contacts (crab_fraction_source) and, when OSP's crabbing-only column arrives, by
  # the OSP counts through f x (1 - c), c the combo-trip share among crabbing boats:
  #     logit f[k] = logit f[k-1] + sigma_f * z[k],   sigma_f learned,  z ~ Student-t
  # A month with five contacts borrows from its neighbours instead of falling back to the
  # 0.30 placeholder; the per-day beta-binomial treats one day's contacted fleet as one
  # correlated observation, not n independent trials. f still enters generated quantities
  # only (effort and CPUE are untouched by construction). FALSE reproduces the legacy
  # construction exactly (the Phase A control rung, run_improvements_2026-09-08.R R3a).
  # The legacy path's crab_fraction_prior_kappa = 20 is why Phase A under-uses the winter
  # contacts: Beta(6, 14) is worth twenty boats at 0.30, so a 46-of-47 December comes out
  # near 0.78; the walk carries no such anchor beyond its first month.
  crab_fraction_dynamic     = TRUE,
  # Level prior of the FIRST stratum of each walk chain, on the logit scale, centred on
  # logit(crab_fraction_set). 1.5 is deliberately weak (95% of prior mass on f in roughly
  # 0.02-0.89); a tight prior at 0.30 would drag a near-1 winter month down by a third.
  crab_fraction_level_sd    = 1.5,
  # Half-normal scale of sigma_f, the SD of one month's step in logit f. The 2024-25
  # monthly contact shares move by RMS ~1.5 logits per month, part of it sampling noise;
  # with ~12 strata sigma_f is weakly identified and this prior matters. Read the
  # sigma_f_out posterior against it (bss_summary_*.csv / crab_fraction_strata_*.csv).
  crab_fraction_walk_sd_prior = 1.5,
  # Student-t degrees of freedom of the walk steps. The two regime jumps of a season (into
  # the winter, when the finfish fisheries close and nearly every boat out is crabbing,
  # and out of it) are real, not noise; heavy-tailed steps let them through without
  # loosening every quiet month. 0 (or negative) gives the Gaussian walk.
  crab_fraction_walk_df     = 4,
  # Beta-binomial concentration of the DAILY contact shares (lognormal centre, log-SD
  # 0.75). Most contact days hold 1-4 boats, where beta-binomial ~ binomial and this
  # rests on its prior; it bites on the busy days (20+ contacts).
  crab_fraction_contact_kappa_prior_mu = 20,
  # THE COMBO-TRIP SHARE c (2026-09-09). The interview workbook now carries every contacted
  # boat's TRIP TYPE (crab only / <fishery> & crab / finfish only / non fishing), so c, the
  # share of crabbing boats that also fished another fishery, is OBSERVED per day and gets
  # the same per-stratum logit walk as f (its own step SD). 2024-25 at the launches: 43
  # combos of 143 crabbing boats, none in Dec-Feb, half or more in the salmon and
  # bottomfish months. c does not move the boat total (f does); it is what OSP's
  # crabbing-only count reads low by, so f(1 - c) (f_lower_out) is the prediction of that
  # column, and when OSP delivers it the comparison is the check on the shift-time
  # contacts against an all-day count. The keys below are the walk's priors.
  crab_fraction_combo_level = 0.3,       # level prior centre (an anchored stratum's c)
  crab_fraction_combo_level_sd = 1.5,    # its logit SD (weak; the level is learned)
  crab_fraction_combo_walk_sd_prior = 1.5,  # half-normal scale of sigma_c
  crab_fraction_combo_kappa_prior_mu = 20,  # beta-binomial concentration of the daily combo shares
  # Which interviews are the f classification: the boats the trailer and OSP counts
  # measure are those launched at the ramp, and a private boat interviewed at the marina
  # or the docks is moored, not trailered (2024-25: 300 contacts at the launches, 34 at the
  # marina, all crab-only, 88 at the docks). NULL = boat_launch_areas; "all" = every
  # private-boat interview (the 2026-09-08 behaviour).
  crab_fraction_contact_areas = NULL,
  # SAMPLER SHIFTS (2026-09-09): 04_input_files/sampler_shifts.xlsx (built from the survey
  # export by build_sampler_shifts.R): check-in / check-out per survey. Read for the
  # shift-coverage diagnostic (shift_coverage_*.csv, contact_hour_by_trip_type.csv): what
  # share of a day's boat returns the shift covers and whether the trip-type mix drifts
  # with the hour inside it. crab_fraction_shift_weighting is the hook for a return-time
  # weighting of the classification once an all-day return profile exists (OSP); "none"
  # is the only implemented value (see sampler_shifts.R for why).
  sampler_shifts_file  = "sampler_shifts.xlsx",
  sampler_shifts_sheet = "data",
  shift_coverage_ie_min_days = 10,     # boat-I/E days the WBL series needs before it is used alone
  crab_fraction_shift_weighting = "none",
  # The PE's per-stratum shrinkage under the dynamic model: Beta(set * kappa, (1 - set) *
  # kappa) with ONE pseudo-contact, so an informed month is essentially its observed share;
  # months under crab_fraction_min_obs are interpolated on the logit scale from their
  # informed neighbours (the design-based partner of the walk, without its smoothing).
  crab_fraction_pe_prior_kappa = 1,
  # osp_scale_is_tau: makes the OSP mean use L (= tau_boat) as the within-day boat
  # turnover, so the dense OSP series (148 days) identifies tau_boat directly (closes
  # GR-12). FALSE keeps the free kappa_OSP OSP scale (Phase 1), which quarantines the
  # turnover question. PRODUCTION DEFAULT = TRUE (adopted 2026-07-31).
  # RESOLVED 2026-07-31: the crab-creel trailer count is an INSTANTANEOUS snapshot, so the
  # OSP/trailer ratio (~2.7-3.0; model kappa_OSP = 3.15) is REAL within-day boat turnover,
  # and the old tau_boat prior (~1.2) was a ~2x UNDER-count. TRUE therefore corrects that
  # bias and raises the boat effort/catch ~2x vs the free-kappa_OSP config (validation
  # Step 6 boat ~27.7k vs Step 3b ~13.8k). Caveat retained: the exact multiplier assumes
  # the trailer snapshot is timed representatively (not a fixed daily peak); see the
  # 2026-07-31 validation review. Requires use_osp_boat_counts = TRUE. Validate by run.
  osp_scale_is_tau          = TRUE,

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
  # 2026-09-08: grid re-centred on the overlap calibration (was 0.9 to 1.8 around the
  # retired 1.2 centre). With shared_tau the boat turnover is data-identified, so the
  # projection is an UPPER bound on prior sensitivity and the diagnostic says so.
  tau_sensitivity_grid     = c(2.0, 2.4, 2.7, 3.0, 3.4),  # tau_boat_prior_mu values to project

  # --- BSS run-level settings (NOT per-fit tuning) -------------------------
  bss_chains        = 4,
  bss_cores         = 4,
  bss_seed          = 20260619,       # fixed seed for reproducible fits
  bss_max_count_seq = 3,              # cap on count sequences per day

  # --- BSS data-sufficiency floor (improvement 6, 2026-08-25) --------------
  # LOWERED 20 -> 15. This was the only thing stopping the private-boat POT-CLOSURE
  # component from ever reaching the sampler: it carries 17 interviews, so at a floor of
  # 20 it was never fitted and entered the port total as a fixed PE point with no
  # interval. Together with the commercial/charter census (also a constant), roughly 18%
  # of the port total carried no uncertainty at all. At 15 that component now attempts a
  # BSS. Both thresholds previously lived in each driver's params_model, where they
  # silently overrode run_config; they now live HERE only.
  #
  # It is surgical: shore all-gear (2,741), shore pot-closure (856) and boat all-gear
  # (145) are all far above either value, so no other component's behaviour changes.
  #
  # HONEST EXPECTATION. 17 interviews over 76 days with ~44 effort counts is very thin.
  # The fit may well fail the convergence gate and fall back to PE anyway -- which is a
  # BETTER outcome than today, because the failure will then be a measured gate verdict
  # with per-fit diagnostics attached rather than a threshold decision made before any
  # data were looked at. Do not read a passing fit here as a precise estimate; read the
  # extrapolation-transparency and coverage rows alongside it.
  bss_min_interviews = 15,
  # SECOND floor, applied AFTER prep drops h <= 0 rows and (when filtering) incomplete
  # trips. The pre-existing guard counts UNFILTERED interviews, so it reads looser than
  # the count the likelihood actually sees; that was harmless when every component had
  # hundreds to spare and is not harmless at 15. NULL = reuse bss_min_interviews.
  bss_min_interviews_fitted = NULL,

  # --- AR resolution PIN (per fit) ------------------------------------------
  # PRODUCTION VALUE IS NULL. Sets a fit's AR resolution exactly, bypassing both the
  # data-driven selector and ar_max_resolution. Scope it per population or per
  # population x sub-season: list(private_boat = "monthly"),
  # list(shore = list(all_gear = "weekly")).
  #
  # 2026-09-09 REFRAME. This began life as an experiment lever, and for one-off
  # experiments it still is one. But in the NEW-SEASON workflow it is the sanctioned
  # second half of the ladder: the ladder (ar_escalate) finds the finest resolution each
  # fit supports, and ar_force is how the user PINS a fit to a chosen rung while deciding,
  # including a rung FINER than the data-driven selector would pick, which the cap alone
  # cannot express (caps only coarsen). Once a season's resolutions are settled, move them
  # into ar_max_resolution and return this to NULL, so production reads from the cap and
  # the config self-documents; the two routes were proven to build byte-identical Stan
  # data on 2026-09-07. Keep the 2026-08-28 lesson in mind: a per-population entry reaches
  # BOTH sub-seasons of that population; scope per sub-season when you mean one.
  ar_force          = NULL,

  # --- PE / BSS incomplete-trip arm alignment (2026-09-02) ---------------------
  # Which interviews the Point Estimator uses for the boat gear-per-group ratio.
  # "match_bss" (default) applies filter_incomplete_trips, matching the frame the BSS
  # learns R_G_boat from; "gear_only" is the pre-2026-09-02 behaviour, which kept an
  # interrupted trip's fully observed gear count. The two disagreed by 0.1% on the boat
  # (3.552 vs 3.550) and at most 0.7% on any component, so this is a consistency fix rather
  # than an accuracy one. See 03_R_functions/pe_gear_ratio_frame.R.
  pe_gear_ratio_arm = "match_bss",

  # --- Zero-inflated catch likelihood (2026-09-02; PROTOTYPE, ships OFF) -------
  # The 2026-09-01 zero-bin read found the NB2 places too little mass at zero on the SHORE
  # catch stream and only there: shore all-gear 676 observed zeros vs 605.4 expected
  # (z = +3.8, n = 1,649), shore pot closure 146 vs 118.5 (z = +2.9), every boat stream
  # inside |z| = 2.3. estimate_catch_zi = TRUE adds a structural-zero component theta_C to
  # the catch likelihood of the populations named below, so a single run carries the shore
  # fits as treatment and the boat fits as an untouched negative control.
  #
  # THE SEASON TOTAL IS SCALED BY (1 - theta_C) in generated quantities. lambda_C is fitted
  # to the non-inflated component and rises to absorb the zeros theta_C removes, so an
  # unscaled total would be inflated by 1/(1 - theta_C) purely by turning the feature on.
  # Judge it on elpd_loo for the catch stream and on the zero bin, not on the total.
  # Changing this forces a Stan recompile.
  # 2026-09-07: ADOPTED alongside the weekly shore AR (ar_max_resolution below). The two
  # changes are additive (+95 crab of interaction on a 21,500 component) and together give
  # the best-behaved fit this project has produced: ZERO Pareto k above 0.7, the only such
  # fit. Compared like for like at weekly, the ZINB halves both count bins on shore
  # all-gear (zero z +3.7 -> +2.0, one z -6.1 -> -3.3) and closes both on the pot-closure
  # replicate, for +11.6 nats at 2.30 paired SE.
  # POOLED ONLY. crab_bss_gear_resolved.stan has no theta_C and prep_bss_crab_gear.R never
  # emits zi_catch, so the gear track silently ignores this flag and fits plain NB2. That
  # makes the two-track cross-check compare UNLIKE shore catch likelihoods for the first
  # time. The effect is small (about -0.3% on the pooled shore component), so it explains
  # a sliver of the cross-track gap and none of a large one, but read the gap knowing it.
  estimate_catch_zi     = TRUE,
  catch_zi_populations  = c("shore"),
  # SEASON-DERIVED prior shape (chosen from the 2024-25 zero bin; weakly informative):
  zi_catch_prior_a      = 1,          # Beta(1, 9): mean 0.10, most mass below 0.25,
  zi_catch_prior_b      = 9,          # comfortably above the ~0.04 the zero bin implies
  #
  # 2026-09-04 RESULT OF THE PROTOTYPE (stage Z1, 2026-09-03 batch), after two corrections
  # to the way it was SCORED. Both corrections went the same way, toward the feature:
  #   elpd, paired: +14.8 nats at a paired SE of 5.5, i.e. 2.69 SE, clearing the stated
  #     2 SE bar. The batch had compared +14.8 against se_elpd_loo (about 46 nats), which
  #     is the SE of ONE model total, dominated by across-observation variation that is
  #     common to both models and cancels in the difference. See loo_elpd_paired.R.
  #   zero bin, under the MIXTURE: 676 observed against 637.9 expected, z = +2.0, from the
  #     NB2 baseline's 605.4 / z = +3.8. The batch computed 419.0 / z = +15.4 because the
  #     R-side PPC scored a ZINB fit as plain NB2. See zinb_ppc.R.
  #   replicate (shore pot closure, a separate fit on different data): theta_C 0.107
  #     against 0.179, elpd +9.1 at 1.96 SE, zero bin z +2.9 -> +0.8.
  # RENDERED RESULT (stage Z2, 2026-09-04 batch; supersedes the 2026-09-05 reading below).
  # Every objection is answered by rendered output:
  #   zero bin  676 observed vs 638.5 expected, z = +2.0, from the NB2's 605.5 / +3.8
  #   ONE bin   236 vs 285.0, z = -3.2, against the NB2's 236 vs 335.4, z = -6.1
  #   catch coverage_50 0.4645 (-2.9 SD) against the NB2's 0.4572 (-3.5 SD)
  #   elpd +14.8 at 2.69 PAIRED SE, replicated at 1.96 SE on the pot-closure fit
  # The 2026-09-05 note said the misfit had MIGRATED to y=1. That was wrong: it could not
  # be checked at the time because the NB2 baseline had no p_one column. The ZINB halves
  # BOTH ends. And the -42.0 nat elpd loss at y=1 is not a shape failure, it is the mass
  # tax a mixture levies on every non-zero observation: log(1 - 0.178) x 236 = -46.3. The
  # r_C tightening (0.947 -> 1.828) repays it at y >= 3, which carries 87% of the catch.
  #
  # SETTLED 2026-09-07 BY STAGE C1 (weekly + ZINB). Compared like for like at the SAME
  # resolution against L1, which is the comparison every earlier run lacked:
  #     shore all-gear  zero bin  z +3.7 -> +2.0     one bin  z -6.1 -> -3.3
  #     shore pot clos. zero bin  z +2.9 -> +0.7     one bin  z -3.5 -> -1.2
  #     elpd +11.6 nats at 2.30 PAIRED SE (it was +14.8 / 2.69 at daily)
  #     Pareto k > 0.7: 1 -> 0, the only fit in this project with none
  # The gain SHRINKS when the effort process is no longer overfitted, which is the right
  # direction and the expected size: it was not an artefact of daily absorbing structure
  # the catch likelihood should have explained. The two changes are ADDITIVE, +95 crab of
  # interaction on a 21,500 component.
  #
  # WHAT IT DOES NOT FIX. The all-gear one bin is still 3.3 SD out and the catch stream is
  # under-covered at every AR resolution: the data has more zeros AND fewer ones than the
  # mixture predicts, i.e. it is more bimodal than a zero-inflated NB can be. A hurdle
  # model, or a two-component NB with its own mean in the low regime, is the shape that
  # fits it; the remaining gain is bounded at roughly 6% of the catch, so that effort can
  # be judged before it is spent.
  #
  # RECOMMENDED: set this TRUE together with the weekly AR below, in ONE adoption render.

  # --- Persist the draws the PPC is computed from (2026-09-04) -----------------
  # Three separate defects in the DIAGNOSTIC arithmetic have each forced a full multi-hour
  # re-fit to correct a file the fit itself was never wrong about. Saving the draw objects
  # the PPC reads makes such a fix a recomputation instead of a re-run. Tens of MB per fit
  # at full draws; save_ppc_draws_max thins if that is too much for a given run.
  # See 03_R_functions/save_bss_ppc_draws.R for what is and is not covered.
  save_ppc_draws        = TRUE,
  save_ppc_draws_max    = NULL,       # NULL keeps every draw, so recomputation is EXACT

  # --- Per-rung adequacy inside the AR ladder (2026-09-06) --------------------
  # Only the fit the ladder KEEPS reaches write_bss_diagnostics(), so on the 2026-09-04
  # run the biweekly and monthly rungs cost 174 minutes and left no p_loo, no Pareto k
  # and no coverage. That is the evidence a ladder exists to produce: all three rungs
  # PASSED the convergence gate, so the gate separated nothing and the adequacy
  # statistics were the only thing that could. With this on, each rung's row in
  # ar_escalation_log.csv carries them. Costs seconds to a minute per rung against a
  # multi-hour fit; only ever computed when the ladder has more than one rung.
  ar_rung_adequacy      = TRUE,

  # --- SHORE ALL-GEAR AR: RESOLVED 2026-09-07, ADOPTION RENDER OUTSTANDING ---------
  # The pooled bracket is now measured across all four rungs (C2, 2026-09-04):
  #     rung   cov50 gear      cov50 catch    total p_loo  bad k   catch   PI rel
  #     daily  0.7010 (+7.1SD) 0.4572 (-3.5)     153.6      41    20,898  0.2974
  #     weekly 0.5659 (+2.3SD) 0.4433 (-4.6)      54.8       1    21,547  0.2530
  #   biweekly 0.5659 (+2.3SD) 0.4500 (-4.1)      34.2       0    21,383  0.2291
  #    monthly 0.5788 (+2.8SD) 0.4506 (-4.0)      22.2       0    20,771  0.2209
  # ONLY DAILY IS AN OUTLIER. The three coarse rungs sit within 0.5 sampling SD of each
  # other and span 3.6% in catch, so the evidence does not choose between them; the agreed
  # rule (finest rung that PASSES the gate) selects WEEKLY, and that is defensible here
  # because weekly is adequate on every statistic rather than merely passing a gate that
  # everything passes. NOTE the correction: the "monthly is also bad, cov50 0.035, -16.4 SD"
  # figure quoted from 2026-08-31 onward is from the GEAR-RESOLVED model, a different
  # likelihood with per-gear CPUE, and never applied to this one.
  # The catch stream is under-covered at EVERY resolution (-3.5 to -4.6 SD). That is a
  # likelihood problem, not a resolution problem; see estimate_catch_zi above.
  #
  # TO ADOPT: set ar_max_resolution$pooled$shore$all_gear to "weekly" (below) together with
  # estimate_catch_zi = TRUE, and render once. The gate on that render is that it must
  # reproduce 05_output/20260904/pooled-CPUE-LZ-C1-weekly-zi bit-identically, because C1
  # reached weekly through ar_force (an experiment lever that bypasses the cap) and
  # production must reach it through the cap: same resolution, different routing, and
  # nothing else may move.

  # --- SUPERSEDED by the block above, kept for the record ---------------------------
  # --- WHERE THE SHORE ALL-GEAR AR STANDS (measured 2026-09-04, NOT yet adopted) ---
  # ar_max_resolution$pooled leaves shore all-gear data-driven, which selects DAILY. The
  # ladder fitted it at daily / weekly / biweekly / monthly and ALL FOUR passed the
  # convergence gate, so the gate cannot choose between them. Adequacy can:
  #        rung   P_n   catch    p_loo/n_obs   bad Pareto k   gear coverage_50
  #        daily  289   20,898       0.352          41         0.701  (+7.1 SD)
  #        weekly  44   21,547       0.096           1         0.559  (+2.1 SD)
  #      biweekly  21   21,383        n/a          n/a           n/a
  #       monthly  10   20,771        n/a          n/a           n/a
  # The four rungs span only 3.7% in catch, so this is a smaller lever on the number than
  # feared; but daily is spending about one effective parameter per three observations and
  # is flagged miscalibrated, and weekly is not. elpd favours daily by 85.8 nats on the
  # effort stream and must NOT be believed: that is 1.08 nats per additional effective
  # parameter, and 41 of the daily fit's 311 gear observations have Pareto k above 0.7, so
  # PSIS-LOO has failed for 13% of the stream being cited. Moving to weekly takes the port
  # total from 71,450 to 72,122 (+0.94%). Biweekly and monthly have NO adequacy because
  # only the rung the ladder keeps reaches write_bss_diagnostics(); ar_rung_adequacy above
  # fixes that for the next run, which is what has to happen before this cap changes.

  # --- Sampler override escape hatch (EXPERIMENTS ONLY; production is NULL) -----
  # Each driver merges its own params_model ON TOP of run_config, so params_model
  # WINS every key it sets, including all the per-fit sampler settings. Setting
  # bss_iter_default here therefore does NOTHING on a normal merge: the driver
  # overwrites it and the run looks like it complied. This named list is applied
  # AFTER the merge (03_R_functions/bss_sampler_override.R), accepts sampler keys
  # only, prints every change, and errors on anything else rather than dropping it.
  #   e.g. list(bss_iter_default = 5000, bss_warmup_default = 2500)
  # to give the gear track 2,500 post-warmup draws instead of its default 1,000.
  bss_sampler_override = NULL,

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
  # SEASON-DERIVED, the WHOLE MAP: every cap below was settled on 2024-25 data (the boat
  # monthly cap by the 2026-07 reconciliation, the shore pot-closure biweekly cap by the
  # Run 6 funnel, the shore all-gear weekly cap by the 2026-09-04 ladder). On a NEW season
  # these are starting points, not answers: run the ladder (ar_escalate, with
  # ar_rung_adequacy = TRUE) and re-derive them. See NEW_SEASON_GUIDE.md.
  ar_max_resolution = list(
    # 2026-09-07 ADOPTION: pooled shore all_gear "daily" -> "weekly". See Section 1k of
    # PIPELINE_STATUS.md. The 2026-09-04 ladder showed daily is OVERFITTED (p_loo 35.2% of
    # n_obs, 41 Pareto k above 0.7, coverage_50 0.701 at +7.1 sampling SD, miscalibration
    # flag set); at weekly those become 9.6%, 1, 0.559 (+2.3 SD) and clear. ONLY daily is
    # an outlier: weekly, biweekly and monthly sit within 0.5 sampling SD of each other and
    # span 3.6% in catch, so the evidence does not choose among them and the agreed rule
    # does, reporting the finest rung that PASSES the gate.
    # NOTE the routing. Stage C1 reached weekly through ar_force, an experiment lever that
    # bypasses both the selector and the cap; this reaches it by letting the data-driven
    # selector pick daily and then coarsening. The two were verified to produce
    # byte-identical Stan data before this edit, and run_adoption_2026-09-07.R gates on the
    # fits being bit-identical to C1. If that gate FAILS, revert this line AND
    # estimate_catch_zi together: the routing would be at fault, not the resolution.
    pooled        = list(shore = list(all_gear = "weekly",
                                      pot_closure = "biweekly"),
                         private_boat = "monthly"),
    gear_resolved = list(shore = list(all_gear = "monthly",
                                      pot_closure = "biweekly"),
                         private_boat = "monthly")
  ),

  # --- AR RESOLUTION ESCALATION LADDER (improvement 7, 2026-08-25) ----------
  # FALSE (default) = the historical behaviour above: one fit per component at the
  # resolution the coverage rule picks, capped by ar_max_resolution, and a failed
  # convergence gate falls straight back to the Point Estimator.
  #
  # TRUE = try-and-escalate. Each component STARTS AT DAILY (the finest rung, ignoring
  # ar_max_resolution) and is refit one rung coarser each time it fails the gate, and
  # STOPS at the first rung that passes. Every component therefore reports at the finest
  # resolution it can actually support, decided by that run's own sampler behaviour
  # rather than by caps hand-tuned from earlier runs. If no rung passes, the component
  # falls back to PE exactly as before. The per-attempt gate verdicts are written to
  # ar_escalation_log.csv and rendered on-page, so the ladder is auditable.
  #
  # COST, STATED PLAINLY. Every rung is a real multi-hour MCMC fit. On the 2024-25
  # config the shore pot-closure fit funnels at daily (~1,165 divergences) and the boat
  # diverged on ~100% of iterations at daily, so those two will burn one and up to three
  # extra fits respectively before settling where ar_max_resolution already puts them.
  # Expect roughly 2-3x the wall clock of a capped run (a 3-6 h run becomes 8-16 h).
  # That is the price of the guarantee; it is why this ships OFF and is a toggle.
  #
  # ar_escalate_ladder is finest-to-coarsest and may be shortened (e.g. drop "weekly")
  # to trade guarantee for runtime. ar_escalate_max_attempts caps the number of fits per
  # component regardless of ladder length.
  # ar_escalate accepts THREE shapes, so a season escalates only the components that need
  # it. Every rung is a real multi-hour fit, and two of the four components already have a
  # known answer (the shore pot closure funnelled at daily on Run 1 with 1,165 divergences;
  # the boat diverged on ~100% of iterations at daily), so escalating those from the top
  # burns known-bad fits to rediscover the caps.
  #
  #   FALSE / TRUE                    off / on for every fit
  #   c("shore")                      on for the named POPULATIONS only
  #   list(shore = "all_gear")        on for named population x sub-season pairs
  #
  # For 2024-25 the component that actually needs it is shore all-gear, so
  #   ar_escalate = list(shore = "all_gear")
  # is the cheap diagnostic setting and TRUE is the exhaustive one.
  ar_escalate             = FALSE,
  ar_escalate_ladder      = c("daily", "weekly", "biweekly", "monthly"),
  ar_escalate_max_attempts = 4,

  # --- HOW THE LADDER STOPS, AND WHICH RUNG IS REPORTED (2026-09-04) -----------
  # ar_escalate_stop:
  #   "first_pass" (DEFAULT, and the behaviour that existed before this key) stop at the
  #      finest rung that passes the convergence gate and report it. This is the PRODUCTION
  #      rule agreed with the FW creel team: start at daily, and if a component cannot
  #      support daily, coarsen until it can. It costs one fit per failed rung and nothing
  #      extra when the first rung passes.
  #   "all_rungs" fit EVERY rung on the ladder, then report the one ar_escalate_select
  #      names. This is the DIAGNOSTIC mode, for honing in on the right resolution in a new
  #      season: it answers "how does the estimate depend on resolution" in ONE run and one
  #      output folder. It costs one fit per rung whether or not earlier rungs passed.
  #
  # ar_escalate_select, used only under "all_rungs":
  #   "first_pass"   report the FINEST rung that passed. Same answer as the production rule,
  #                  but with the whole ladder measured and logged alongside it.
  #   "narrowest_pi" among the rungs that passed, report the one with the narrowest RELATIVE
  #                  predictive interval (hi95 - lo95) / median.
  #
  # A WARNING ABOUT "narrowest_pi", because the instinct behind it is reasonable and the
  # data does not support it. On the 2026-08/09 private-boat 2x2, four cells on one track:
  #
  #     cell            catch    PI width   rel width   trailer coverage_50 (nominal 0.500)
  #     OFF x monthly  25,868      40,358      156.0%   0.538   calibrated
  #     ON  x monthly  31,008      47,671      153.7%   0.523   calibrated
  #     OFF x daily    37,359      54,769      146.6%   0.713   BROKEN
  #     ON  x daily    42,344      60,294      142.4%   0.744   BROKEN
  #
  # On RELATIVE width the two miscalibrated cells look the most precise, because a latent
  # process that absorbs observation noise reports a tighter interval than one that does
  # not. On ABSOLUTE width the calibrated cell wins, but only because catch and width are
  # nearly proportional here (the ratio is 0.64 in all four cells), so "narrowest absolute
  # interval" is close to "smallest estimate" and would systematically select the lowest
  # harvest number. Use precision to break a tie between rungs that are ALREADY adequate,
  # never to decide adequacy. That is why the default is "first_pass".
  #
  # WHAT THE GATE CANNOT SEE, stated because it bounds what this whole mechanism can do.
  # The ladder escalates on the CONVERGENCE gate, which asks whether a fit sampled. On the
  # 2026-08-31 production run every component passes it, including shore all-gear at daily,
  # so the production rule would change nothing this season. (That is exactly how it played
  # out: the 2026-09-04 ladder saw every rung pass the gate, adequacy alone separated them,
  # and weekly was adopted through ar_max_resolution rather than through this rule.) That same fit carries p_loo at
  # 35.2% of n_obs, 41 Pareto k above 0.7 and coverage_50 0.701 against a nominal 0.500,
  # all of which the gate is blind to and all of which model_adequacy.csv reports beside it.
  # If a future season wants the ladder to react to adequacy as well as to sampling, that is
  # a deliberate change to the gate, not a setting here.
  ar_escalate_stop        = "first_pass",
  ar_escalate_select      = "first_pass",
  # When escalating, should ar_max_resolution still apply as a floor on how FINE the
  # ladder may start? FALSE (default, and the point of the feature) starts at the top of
  # the ladder for every component. TRUE keeps the caps and escalates only from the
  # capped rung, which is the cheap variant: it costs no extra fits on a config whose
  # caps are already right, but it cannot discover that a cap is too coarse.
  ar_escalate_respect_cap = FALSE,

  # --- Input workbooks (all .xlsx, single "data" sheet) --------------------
  # Every pipeline input is an .xlsx workbook (converted from CSV on 2026-07-16) with
  # one "data" sheet; dates are ISO yyyy-mm-dd text. Filenames are parameters so a
  # workbook can be swapped without touching the readers. The I/E, holiday, and
  # fishery-opener workbooks keep their own filename keys in the sections below.
  effort_file       = "effort_combined.xlsx",
  interview_file    = "interview_combined.xlsx",
  tally_file        = "wes_commercial_tally.xlsx",
  input_sheet       = "data",         # sheet name shared by the flat input workbooks
  # 2026-09-10: every workbook is BUILT from the per-season creel workbooks in
  # 04_input_files/raw/ by the builders next to them (build_all_inputs.R runs the set), and
  # every one now spans 2022-23 through 2025-26 (through 2026-09-08). Readers go through
  # read_input_workbook() (guess_max over the whole column, since early seasons leave some
  # columns blank). The effort workbook carries a qc_flag; rows with a flag in this vector
  # are held out of the counts (build_effort_combined.R lists the flags and their counts):
  effort_qc_drop    = "gear_count_from_interviews",   # Feb 2024 PFD stand-down: day totals recorded as counts; "none" keeps every row
  # The charter trip roster (charter_trips.xlsx, 2024-25 on): the TRIP FRAME the charter
  # expansion runs on. charter_frame = "roster" (default) takes, per day, the larger of the
  # roster's sailed trips and the tally's charter column (the roster does not depend on a
  # sampler being in port: 8 of the 31 sailed Westport charter trips of 2024-25 fell on days
  # without a tally); "tally" uses the tally column alone, which reproduces the 2026-09-09
  # arithmetic (7,884 on 2024-25) and leaves the frame incomplete by construction. See
  # estimate_comm_charter.R.
  charter_frame       = "roster",
  charter_trips_file  = "charter_trips.xlsx",
  charter_trips_sheet = "data",
  charter_roster_port = "Westport",

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
  # EFFECTIVE DAY LENGTH, and what production actually uses (improvement 2).
  # The I/E surveys give two derived quantities per survey day, from the same 15-minute
  # presence series:
  #   L_effective = crabber-hours / peak crabbers present   (~5.3 h mean; "hours of
  #                 crabbing one snapshot crabber is worth")
  #   turnover    = crabber arrivals / peak crabbers present (~1.72 over 30 WDF20 days)
  # Their ratio is the implied trip length, 5.26 / 1.72 = 3.07 h, against an interview
  # mean trip length of 3.23 h -- two independent sources 5% apart, which is the free
  # consistency check on the method.
  # PRODUCTION USES THE TURNOVER, NOT THE DAY LENGTH. Since shore moved to
  # gear-deployments (v7.7) the expansion is E = lambda_E * R_G * tau_shore, so L in the
  # Stan model is tau_shore (~1.7, dimensionless), NOT L_effective in hours. L_effective
  # is still computed every run and drives the day_length column used by the diagnostics
  # and the civil-twilight comparison, but it is not on the estimation path unless
  # shore_effort_unit is set back to a time unit. The regression fallback ladder below
  # (regression -> grand mean -> civil twilight) governs L_effective only.
  # BOATS: L is tau_boat (deployment turnover), NOT 24 h. The old flat L = 24 gear-hours
  # was removed in POOL-3 / v7.6 when the boat moved to gear-deployments.
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

  # --- OTHER-FISHERY OPENER EFFORT COVARIATES (improvement 4, 2026-08-25) ---
  # Generalizes the old razor-dig-only B3 term. Any opener in the consolidated calendar
  # can now enter the EFFORT model of either population as a log-additive day covariate
  # (X_open / B_open in both Stan models), with per-population candidate lists. The old
  # razor_dig_mode below still works and is implemented on top of this machinery, so a
  # razor-only run is mathematically identical to the previous B3 term.
  #
  #   "off"    no opener covariates (the historical default; K_open = 0). Posterior-
  #            identical to a pre-2026-08-25 run, and bit-identical for the gear-resolved
  #            model; the POOLED model is not bit-identical, because the retired `real B3`
  #            was an always-declared sampled dimension whose removal shifts the RNG
  #            stream. Judge a baseline run on medians and intervals within Monte Carlo
  #            error, not by diffing CSVs.
  #   "auto"   include a candidate only if the spillover diagnostic's day-type + month
  #            ADJUSTED effect on THAT population's effort series clears opener_auto_p
  #            after the multiplicity adjustment below, AND the term is identifiable in
  #            the fit's own window (>= opener_min_days days on each side).
  #   "manual" include exactly opener_manual_shore / opener_manual_boat.
  #
  # WHY THE MULTIPLICITY ADJUSTMENT IS NOT OPTIONAL. The diagnostic runs 8 effort tests
  # (4 openers x 2 populations); at a raw alpha of 0.05 about one false positive is
  # expected by chance, and the project already has a worked example of exactly that
  # failure: razor-dig -> shore effort came in at p = 0.045, was fitted as B3 in Run 3,
  # and delivered NO predictive gain (elpd_loo within 1 SE on every stream). Under the
  # default Benjamini-Hochberg adjustment over the 8 effort tests that term correctly
  # does NOT survive, while MA2 halibut -> boat trailers (p = 3.6e-06) does. Selecting on
  # a p-value from the same data that fits the term also inflates the coefficient
  # (post-selection / winner's-curse bias), so treat "auto" as a screen and confirm any
  # selected term against a covariate-free run's elpd_loo before citing it.
  #
  # READ THIS BEFORE SWITCHING ON A BOAT OPENER. Boat effort is a count of ALL private
  # boats, so an MA2-halibut term is legitimate there (it explains real variance in the
  # all-boat series) -- but the halibut boats it captures are not crabbing, and a
  # CONSTANT crabbing fraction f then converts that surge to crab effort at the same rate
  # as any other day. Fitting the surge better while multiplying it by a flat f makes the
  # boat catch MORE biased on those days, not less. Pair any boat opener covariate with
  # crab_fraction_strata = "opener" or "month_opener" (see the crab-fraction block).
  opener_covariate_mode = "off",     # "off" | "auto" | "manual"
  opener_candidates_shore = c("razor_nearby_dig"),
  opener_candidates_boat  = c("ma2_halibut_open", "ma2_salmon_open", "ma2_bottomfish_open"),
  opener_manual_shore   = character(0),
  opener_manual_boat    = character(0),
  opener_auto_p         = 0.05,      # threshold on the ADJUSTED p (see opener_auto_p_adjust)
  opener_auto_p_adjust  = "BH",      # "BH" | "bonferroni" | "none"  (over the 8 effort tests)
  opener_min_days       = 10,        # min days OPEN and min days CLOSED inside a fit's window
  # CPUE is deliberately NOT offered an opener covariate: all eight opener-vs-catch-rate
  # tests in the 2024-25 diagnostic came back null (p = 0.29 to 0.90).

  # --- razor_dig SHORE-effort term (item 1, 2026-07-13; now an alias) -------
  # Retained compatibility switch, implemented as one column of the opener machinery
  # above. "no" = off (production default); "yes" = force razor_nearby_dig into the SHORE
  # effort model; "auto" = include it only if the adjusted shore-effort razor p clears
  # razor_dig_auto_p (RAW p, unadjusted -- the historical behaviour). Boat fits and
  # inactive runs get no razor column, so the term stays absent rather than decoupled.
  # Run 3 disqualified razor-dig B3 (no elpd gain); keep off unless deliberately
  # re-testing. Setting BOTH razor_dig_mode != "no" and a razor candidate under
  # opener_covariate_mode is de-duplicated (one column, no double counting).
  # RE-COMPILES the Stan model.
  razor_dig_mode   = "no",       # "no" | "yes" | "auto"
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
  gear_resolved_G            = FALSE, # (gear-resolved) GR-7 Phase 1. FALSE = production G = 1
                                       #   (gear split by PE apportionment). TRUE turns on genuine
                                       #   per-gear CPUE for SHORE fits (Option A1): only single-gear
                                       #   interviews feed a gear-specific CPUE, multi-gear trips form
                                       #   "Mixed", effort is split across gears by the pi_gear shares.
                                       #   Boat stays G = 1 (Pot-dominated; Phase 0). Changes shore
                                       #   inference, so validate by run. See
                                       #   07_documentation/development_notes/GR-7-per-gear-CPUE-design.md
  gear_share_dirichlet       = FALSE, # (gear-resolved) GR-7 Phase 2. FALSE = the gear split is the
                                       #   fixed pi_gear point estimate (Phase 1; output unchanged).
                                       #   TRUE replaces it with a Dirichlet posterior
                                       #   pi_gear ~ Dirichlet(alpha0_gear + weighted interview counts)
                                       #   per period x day_type, so share-sampling uncertainty widens
                                       #   the per-gear catch total C_sum_gear (C_sum is untouched).
                                       #   Requires gear_resolved_G = TRUE; inert (forced off) at G = 1.
  alpha0_gear                = 1.0,    # (gear-resolved) GR-7 Phase 2 Dirichlet concentration floor.
                                       #   1.0 = a flat (uniform) prior added to the counts; smaller
                                       #   (e.g. 0.5) lets sparsely-sampled period x day_type cells
                                       #   concentrate more on their observed gear, larger regularizes
                                       #   toward an even split. Only used when gear_share_dirichlet.
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
