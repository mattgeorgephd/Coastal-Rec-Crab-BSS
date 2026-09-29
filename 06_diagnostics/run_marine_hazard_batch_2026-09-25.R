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
# MARINE HAZARD EFFORT COVARIATES: SMALL CRAFT ADVISORIES AND BAR RESTRICTIONS (2026-09-25)
# -----------------------------------------------------------------------------
# WHY THIS RUN EXISTS
#   The samplers record "Small Craft Advisory" and "Bar Restrictions" on their survey form.
#   On 2026-09-25 the NWS VTEC archive was pulled for the Grays Harbor Bar (PZZ110) and the
#   coastal waters off Westport (PZZ156) and compared with those ticks: the sampler SCA tick
#   is a noisy copy of the archive (it misses about a third of the advisories in effect
#   during a shift), while the bar-restriction tick carries information no NWS product has
#   (a USCG closure of the bar to recreational vessels). Offline, on the 2023-24 to 2025-26
#   Westport Boat Launch trailer counts after season, month and day type, an SCA-or-higher
#   in the 04:00-16:00 window carried a rate ratio of about 0.25 and the bar tick about 0.67
#   beyond every archived NWS product. Those are associations on sampled days in a
#   quasi-Poisson regression. The question this run answers is different: does either EARN
#   a term in the BSS effort process, whose whole purpose is to carry effort across the
#   unsampled days?
#
#   The machinery landed with this file (03_R_functions/bss_marine_hazard_covariates.R,
#   run_config.R section 4.4b): the covariates ride on the existing K_open / X_open / B_open
#   block, so NEITHER STAN MODEL CHANGED and marine_hazard_mode = "off" (the shipped value
#   until 2026-09-27; since then the pre-adoption model) builds the Stan data it built before.
#   M1 below proves that empirically.
#
# THE RUNGS. One lever moves per rung; everything else is pinned by WINDOW below.
#   M0   desk, seconds. No fit. The archive covers the window; the flags and the screen
#        build; each rung differs from M1 in declared keys only; the Stan data of M4 differs
#        from M1's in K_open and X_open_flat ONLY (built for real, for the shore all-gear and
#        the boat all-gear components); the shipped config is reported (off until 2026-09-27,
#        the adopted method since).
#   M1   marine_hazard_mode = "off"      the baseline AND the inertness proof: must be
#        bit-identical to the committed R4 render (20260910/pooled-CPUE-IMP-R4-shore-tau-newf)
#   M2   SCA only    nws_sca_any forced on BOTH populations (manual)
#   M3   bar only    bar_restriction forced on the boat (manual)
#   M4   both        nws_sca_any on both populations + bar_restriction on the boat
#   M5   auto        marine_hazard_mode = "auto": what production would do with the screen
#   M6   split       (added 2026-09-27, B41, for D32; RENDERED 2026-09-27, Section 1z) the any-zone SCA flag on
#        the boat in TWO columns, nws_sca_any_winter (December to February) and
#        nws_sca_any_rest (every other day), no shore term. Judged AGAINST M2, the constant
#        term: rule 3 for each coefficient, rules 4 and 5 for the pair, and the block CV's
#        joint row (run_marine_block_cv_2026-09-26.R, M6 vs M2). What it asks: is the winter
#        effect the summer effect? Section 1y.3 found the term moves the boat estimate mostly
#        in the winter, where the block CV cannot see, and a crude winter-only screen put the
#        winter rate ratio near 0.5 to 0.6 against 0.15 to 0.16 for the rest of the season.
#
# THE DECISION RULE, STATED HERE BEFORE THE RUN so it cannot be fitted to the answer.
#   1. A rung is ELIGIBLE only if every fit in it passes the convergence gate, and every fit
#      used the SAME AR resolution as M1 (a covariate that only changes a resolution is not
#      a covariate result).
#   2. Adequacy must not degrade: shore all-gear and boat all-gear p_loo_frac <= 0.15 and
#      bad Pareto k <= 5% of the catch-stream n_obs, the same clauses as the D3 run.
#   3. A covariate is IDENTIFIED in a fit if the 95% interval of its B_open_out excludes 0.
#      A term that is not identified is not adopted, whatever the elpd says.
#   4. THE PAIRED EFFORT-STREAM elpd against M1 (gear stream for shore fits, trailer stream
#      for boat fits; loo_elpd_paired(), Vehtari et al. 2017 s3.3): a gain above +2 paired
#      SE on at least one fit where the term is active, and no fit worse than -2 SE, is
#      "earns its term ON THE SAMPLED DAYS". Between -2 and +2 SE is "no evidence either
#      way", and that is REVIEW, not FAIL.
#   5. The CATCH stream's paired elpd must sit within +-2 SE of M1: an effort covariate has
#      no business moving the catch fit. Outside that band is REVIEW.
#   6. bar_restriction is judged in M4 AGAINST M2, not against M1: on unsampled days it is an
#      expected value driven by the archived advisories, so against M1 it would be credited
#      with the SCA term's work. Its increment beyond the archive is what the offline screen
#      claimed (0.67) and what M4-vs-M2 measures.
#   7. THE PORT TOTAL IS NOT A CRITERION. It is reported at every rung because it is what a
#      reader will quote, and a covariate that moves it is the thing to investigate, not a
#      reason to prefer or reject the rung.
#   8. WHAT THIS RUN CANNOT MEASURE, stated up front. Observation-level LOO scores the
#      SAMPLED days. The covariate's value to the estimate is on the UNSAMPLED days, where
#      the AR process interpolates and an advisory flag is the only information that a day
#      was not an ordinary day; nothing in a PSIS-LOO on sampled observations tests that
#      ("LOO can't test what it never held out"). So this run can say a term is identified,
#      harmless to the catch fit and better on the sampled days; it cannot say the
#      unsampled-day interpolation improved. A leave-one-week-out block CV (CHANGE_REGISTER
#      D31) is the test that could, and it is not built. The recommendation block below
#      says which of the two it is resting on. [2026-09-26: built and run, B39, Section 1y.]
#   9. THE SEASON SPLIT (M6, stated 2026-09-27 before its render; the fit judged is the BOAT
#      ALL-GEAR fit: in the pot-closure window the winter column has no flagged day, is
#      dropped, and that fit is M2's). Read in this order, and stop at the first that applies:
#      (a) The winter coefficient is NOT identified (rule 3: its 95% interval covers 0).
#          Then this season's counts cannot place the winter effect, D32 stays open, and
#          NOTHING here is a verdict for or against the constant term: the constant term's
#          winter rests on an assumption the data cannot check, and the split's winter is a
#          free parameter (rule 3's own objection). The rest coefficient is still read: if it
#          is identified and the block CV's joint boat all-gear row for M6 against M2 is at
#          or above -2 paired SE, the split is not worse than the constant term where it can
#          be tested. The winter's treatment is then a modelling choice for the register
#          (D32), not a rule outcome. STATED BEFORE THE RENDER: this is the LIKELY outcome.
#          The module's own screen puts the 2024-25 winter term at 0.57 with a log-scale SE
#          of 0.8, and the multi-season desk screen (06_diagnostics/desk_sca_season_split_
#          2026-09-27.R; 794 sampled days over four seasons) puts the winter effect at 0.43
#          against 0.22 for the rest, so a coefficient near log(0.5) with the 2024-25 winter's
#          36 sampled advisory days of five trailers or fewer will very likely straddle zero.
#          M6's value in that case is the identified rest coefficient, the exact OSP check
#          on a rendered fit, and the recorded fact that the winter cannot be told apart.
#      (b) Both identified, and they DIFFER: the 95% interval of (winter - rest), taken as
#          independent normals from the two posterior means and SDs (bss_full_summary; the
#          two act on disjoint days, so the approximation errs a little conservative),
#          excludes 0. Then the split EARNS ITS PLACE over the constant term if the catch
#          stream is unmoved against M2 (rule 5) and the block CV's joint boat all-gear row
#          for M6 against M2 is at or above -2 paired SE (it can only see the rest
#          coefficient's weeks, so it is asked not to contradict, not to confirm; a gain
#          there is support). Rule 4 here, the paired trailer elpd against M2, is reported
#          and is not the arbiter (Section 1y.4).
#      (c) Both identified, and they do NOT differ (that interval covers 0). Then the
#          constant term stands and D32 is answered "no difference detectable in this
#          season's counts", which is not "the same size".
#      The port total is reported, never judged (rule 7).
#
# WHAT IT WRITES, merged by key so a partial re-run updates rather than truncates:
#   05_output/marine_hazard_2026-09-25_ladder.csv          per-rung totals, adequacy, B_open
#   05_output/marine_hazard_2026-09-25_verdicts.csv        every criterion, PASS/FAIL/REVIEW/INFO
#   05_output/marine_hazard_2026-09-25_recommendation.csv  SCA and bar, with the rule applied
#
# RUNTIME. The R4 render took 207 minutes (run_timings.csv: shore all-gear 150 min at
# weekly, boat all-gear 36 min at monthly). A covariate adds one parameter per active fit
# and does not change P_n, so budget the same per rung: about 3.5 h x 5 fitted rungs, 17 to
# 18 h (M1 to M5 rendered 2026-09-25/26; with RESUME they are read back in a minute and only
# M6, about 3.5 h, renders). STAGES below can drop rungs; M1 and M2 are the minimum that
# answers anything.
#
# WHERE THIS ENDED (2026-09-27, VALIDATION_CAMPAIGN Section 1z). Six rungs rendered. Matt's
# decision: the CONSTANT boat SCA term, the all-gear fit only (no shore term, not the bar
# tick, not the pot-closure fit), now the method of record in run_config.R (A30 ADOPTED).
# This runner keeps its rungs reproducible under the new run_config by setting the marine
# keys it needs in resolve_cfg(), and RESUME reads all six back in a minute. M0 also lays
# run_config.R's own marine keys over the window and builds all four production fits on the
# real inputs (row 4c): the term must sit in the boat all-gear fit and nowhere else.
#
# SHIPS DRY_RUN <- TRUE. Set it FALSE and source again to fit.
###############################################################################

DRY_RUN <- TRUE                    # TRUE prints the plan and runs M0; fits nothing
STAGES  <- c("M0", "M1", "M2", "M3", "M4", "M5", "M6")   # M6 added 2026-09-27 (B41); M1 to M5 RESUME from their folders
RESUME  <- TRUE                    # reuse a rung ONLY when its MH_STAGE.txt digest matches

# =========================================================================== #

.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) &&
    dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions")))
  stop("Run this from the repository root (or from 06_diagnostics/): 03_R_functions not found.")
.here <- function(...) file.path(.root, ...)
setwd(.root)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
banner <- function(msg) cat("\n", strrep("=", 78), "\n ", msg, "\n", strrep("=", 78), "\n", sep = "")
rule   <- function() cat(strrep("-", 78), "\n")
fmt <- function(x, d = 1) {
  if (length(x) == 0) return("NA")
  out <- formatC(suppressWarnings(as.numeric(x)), format = "f", digits = d, big.mark = ",")
  out[is.na(suppressWarnings(as.numeric(x)))] <- "NA"
  out
}
.num1 <- function(x) { v <- suppressWarnings(as.numeric(x)); if (length(v)) v[1] else NA_real_ }

if (!isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages({ library(here); library(rmarkdown) })
  # 2026-09-28 (B46): the shared loader (03_R_functions/bss_packages.R): renv.lock versions, a stop on a missing package.
  source(file.path(.root, "03_R_functions", "bss_packages.R")); bss_load_packages()
  rstan_options(auto_write = TRUE)
} else {
  suppressWarnings(suppressPackageStartupMessages(
    try({ library(dplyr); library(tidyr); library(readr); library(lubridate)
          library(readxl); library(here); library(purrr); library(stringr); library(tibble) },
        silent = TRUE)))
}
invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE),
                 function(f) source(f)))   # B46: a file that fails to source stops the runner (it was hidden by try())
source(.here("run_config.R"))
BASE <- run_config

POOLED_RMD  <- .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd")
POOLED_STAN <- .here("02_stan_models", "crab_bss_pooled.stan")
POOLED_PREP <- .here("03_R_functions", "prep_bss_crab_pooled.R")
MH_MODULE   <- .here("03_R_functions", "bss_marine_hazard_covariates.R")
PREFIX <- "pooled-CPUE-"
stopifnot(file.exists(POOLED_RMD), file.exists(POOLED_STAN), file.exists(POOLED_PREP), file.exists(MH_MODULE))

# THE REFERENCE the inertness control compares against: the authoritative R4 render, on the
# shipped configuration, made before this file and the covariate machinery existed.
REF_R4 <- "20260910/pooled-CPUE-IMP-R4-shore-tau-newf"

# the two fits the adequacy clauses read, and the fit each covariate is active in
FIT_SHORE <- "shore_all_gear_Dungeness_Kept"
FIT_BOAT  <- "private_boat_all_gear_Dungeness_Kept"
FITS_ALL  <- c("shore_ring_net_only_Dungeness_Kept", FIT_SHORE,
               "private_boat_ring_net_only_Dungeness_Kept", FIT_BOAT)
.stream_of <- function(fit) if (startsWith(fit, "shore")) "gear" else "trailer"

rd <- function(dir, f) {
  p <- file.path(dir %||% "", f); if (!file.exists(p)) return(NULL)
  tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}
.port_row <- function(x) {
  if (is.null(x) || !all(c("Estimate", "BSS_median") %in% names(x))) return(NULL)
  i <- which(grepl("^(Expected_)?Catch$", x$Estimate))
  if (!length(i)) return(NULL)
  as.list(x[i[1], , drop = FALSE])
}
.comp <- function(dir, key, col = "BSS_catch") {
  x <- rd(dir, "pe_vs_bss_comparison.csv"); if (is.null(x)) return(NA_real_)
  .num1(x[[col]][x$component == key])
}
.full <- function(dir, fit) {
  f <- file.path(dir %||% "", sprintf("bss_full_summary_%s.csv", fit))
  if (!file.exists(f)) return(NULL)
  tryCatch(utils::read.csv(f, row.names = 1, check.names = FALSE), error = function(e) NULL)
}
.adq <- function(dir, fit) {
  x <- rd(dir, "model_adequacy.csv"); if (is.null(x) || !"fit" %in% names(x)) return(NULL)
  i <- which(x$fit == fit); if (!length(i)) return(NULL)
  as.list(x[i[1], , drop = FALSE])
}
.arlog <- function(dir, fit) {
  x <- rd(dir, "ar_escalation_log.csv"); if (is.null(x) || !"fit" %in% names(x)) return(NULL)
  i <- which(x$fit == fit); if (!length(i)) return(NULL)
  as.list(x[i[length(i)], , drop = FALSE])
}
.gate <- function(dir) rd(dir, "convergence_report.csv")
.loo_summ <- function(dir, fit, stream) {
  x <- rd(dir, sprintf("loo_summary_%s.csv", fit)); if (is.null(x)) return(NULL)
  i <- which(x$stream == stream); if (!length(i)) return(NULL)
  as.list(x[i[1], , drop = FALSE])
}
# The covariate labels behind B_open_out[k] for one fit: opener_covariates_<fit>.csv, written
# by the driver whenever K_open > 0 (index, parameter, opener).
.cov_labels <- function(dir, fit) {
  x <- rd(dir, sprintf("opener_covariates_%s.csv", fit))
  if (is.null(x) || !all(c("parameter", "opener") %in% names(x))) return(NULL)
  stats::setNames(as.character(x$opener), as.character(x$parameter))
}
# B_open_out posterior for a named covariate in one fit: mean, 2.5%, 97.5%, and whether the
# interval excludes 0. NULL when the term is not in that fit.
.b_open <- function(dir, fit, covariate) {
  lab <- .cov_labels(dir, fit); fs <- .full(dir, fit)
  if (is.null(lab) || is.null(fs) || !covariate %in% lab) return(NULL)
  par <- names(lab)[match(covariate, lab)]
  if (!par %in% rownames(fs)) return(NULL)
  r <- fs[par, , drop = FALSE]
  list(parameter = par, mean = .num1(r[["mean"]]), lo = .num1(r[["2.5%"]]), hi = .num1(r[["97.5%"]]),
       rhat = .num1(r[["Rhat"]]), n_eff = .num1(r[["n_eff"]]),
       identified = isTRUE(is.finite(.num1(r[["2.5%"]])) && is.finite(.num1(r[["97.5%"]])) &&
                             (.num1(r[["2.5%"]]) > 0 || .num1(r[["97.5%"]]) < 0)))
}
.pct <- function(a, b) if (isTRUE(is.finite(a)) && isTRUE(is.finite(b)) && b != 0) 100 * (a - b) / b else NA_real_

V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = observed,
                                    threshold = threshold, verdict = verdict, why = why,
                                    stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# THE PIN. Every key that is NOT a declared delta is fixed here AT THE SHIPPED VALUE, so a
# rung differs from M1 in exactly the marine lever it names, and M1 differs from the
# committed R4 render in nothing but the keys that did not exist when R4 rendered.
# ---------------------------------------------------------------------------
WINDOW <- list(
  # the season, all nine per-season keys
  est_date_start = "2024-09-16", est_date_end = "2025-09-15", season_filter = "2024-25",
  pot_closures = NULL, census_windows = NULL,
  pot_closure_start = "2024-09-16", pot_closure_end = "2024-11-30", pot_open_date = "2024-12-01",
  census_start_date = "2024-12-01", census_end_date = "2025-02-08", commercial_opener = "2025-02-11",
  # the census and PE levers (they enter the port total)
  census_expansion = "none", census_uncertainty = "charter", charter_expansion = "vessel", charter_frame = "roster",
  pe_empty_effort_stratum = "local_day_type", pe_empty_stratum = "local", pe_variance = "impute_aware",
  # the turnover and f levers, at the shipped RESOLVERS (not resolved numbers)
  tau_shore_prior_mu = "derived", tau_shore_prior_sigma = "derived",
  tau_boat_prior_mu = "calibration", tau_boat_prior_sigma = 0.5,
  shared_tau = TRUE, shared_tau_sigma = 0.15, shared_tau_min_obs = 15, osp_scale_is_tau = TRUE,
  use_osp_boat_counts = TRUE, use_crab_fraction = TRUE,
  crab_fraction_strata = "month", crab_fraction_source = "both", crab_fraction_dynamic = TRUE,
  use_osp_crab_lower = FALSE,
  # the shore catch likelihood
  estimate_catch_zi = TRUE, catch_zi_populations = "shore", catch_zi_tracks = "pooled",
  zi_catch_prior_a = 1, zi_catch_prior_b = 9,
  # the OTHER day covariates stay off: the K_open block must carry the marine columns alone
  opener_covariate_mode = "off", razor_dig_mode = "no",
  # sampler and AR: one seed, the shipped caps, no ladder, no override
  bss_seed = 20260619, bss_chains = 4, bss_cores = 4, bss_sampler_override = NULL,
  # estimate_red_rock: retired 2026-09-28 (B44), kept in the pin so no stage digest moves
  ar_force = NULL, ar_escalate = FALSE, ar_rung_adequacy = TRUE, estimate_red_rock = FALSE,
  # the flag definition, pinned so a run_config edit cannot redefine the covariate mid-batch
  marine_hazard_candidates_shore = c("nws_sca_any"),
  marine_hazard_candidates_boat  = c("nws_sca_any", "bar_restriction"),
  marine_hazard_auto_p = 0.05, marine_hazard_auto_p_adjust = "BH",
  marine_hazard_file = "nws_marine_hazards.xlsx", marine_hazard_sheet = "data",
  marine_hazard_zones = c(bar = "PZZ110", coastal = "PZZ156"),
  marine_hazard_codes = c("SC.Y", "RB.Y", "SW.Y", "SI.Y", "GL.W", "SR.W", "SE.W", "HF.W"),
  marine_hazard_window = c(4, 16), marine_hazard_tz = "America/Los_Angeles",
  bar_restriction_impute = "nws", bar_restriction_field_start = NULL
)

# ---------------------------------------------------------------------------
# THE RUNGS. `delta` is the ONLY thing that may differ from the pinned configuration.
# ---------------------------------------------------------------------------
STAGE_DEFS <- list(
  M0 = list(tag = "MH-M0-desk", fit = FALSE, item = "prerequisites; no fit", delta = list()),
  M1 = list(tag = "MH-M1-off",  fit = TRUE,  item = "baseline: marine_hazard_mode = off (the pre-adoption model; the shipped value until 2026-09-27); must reproduce R4 bit for bit",
            delta = list(marine_hazard_mode = "off")),
  M2 = list(tag = "MH-M2-sca",  fit = TRUE,  item = "SCA only: nws_sca_any on shore AND boat",
            delta = list(marine_hazard_mode = "manual", marine_hazard_manual_shore = "nws_sca_any",
                         marine_hazard_manual_boat = "nws_sca_any")),
  M3 = list(tag = "MH-M3-bar",  fit = TRUE,  item = "bar only: bar_restriction on the boat",
            delta = list(marine_hazard_mode = "manual", marine_hazard_manual_shore = character(0),
                         marine_hazard_manual_boat = "bar_restriction")),
  M4 = list(tag = "MH-M4-both", fit = TRUE,  item = "both: nws_sca_any on shore and boat, bar_restriction on the boat",
            delta = list(marine_hazard_mode = "manual", marine_hazard_manual_shore = "nws_sca_any",
                         marine_hazard_manual_boat = c("nws_sca_any", "bar_restriction"))),
  M5 = list(tag = "MH-M5-auto", fit = TRUE,  item = "auto: the screen decides (what production would do with marine_hazard_mode = auto)",
            delta = list(marine_hazard_mode = "auto", marine_hazard_manual_shore = character(0),
                         marine_hazard_manual_boat = character(0))),
  # B41 (2026-09-27). The split pair must be OFFERED to be named under `manual`, so this rung
  # also widens the boat candidate list; the two extra keys are declared per stage (see
  # stage_digest() and M0's comparability check), so the M1 to M5 digests and folders stand.
  M6 = list(tag = "MH-M6-split", fit = TRUE,  item = "SCA split by season on the boat: nws_sca_any_winter (Dec-Feb) + nws_sca_any_rest; no shore term (judged against M2)",
            delta = list(marine_hazard_mode = "manual", marine_hazard_manual_shore = character(0),
                         marine_hazard_manual_boat = c("nws_sca_any_winter", "nws_sca_any_rest"),
                         marine_hazard_candidates_boat = c("nws_sca_any", "bar_restriction", "nws_sca_any_winter", "nws_sca_any_rest"),
                         marine_hazard_winter_months = c(12L, 1L, 2L)))
)
if (!all(STAGES %in% names(STAGE_DEFS)))
  stop("STAGES names a stage that does not exist: ", paste(setdiff(STAGES, names(STAGE_DEFS)), collapse = ", "))
DELTA_KEYS <- c("marine_hazard_mode", "marine_hazard_manual_shore", "marine_hazard_manual_boat")
# the keys a stage may differ from M1 in: the three levers above plus whatever its own delta
# names (M6 adds the widened candidate list and the winter months). Per stage, so that adding
# M6 changed no earlier stage's declared set or digest (the M1 to M5 folders RESUME as before).
declared_keys <- function(sid) unique(c(DELTA_KEYS, names(STAGE_DEFS[[sid]]$delta %||% list())))
# keys the driver ADDS to params at run time (data, not configuration) that config_delta()
# would otherwise report between two folders
RUNTIME_KEYS <- c("run_tag", "model", "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates",
                  "crab_fraction_rows", "osp_crab_rows", "tau_boat_prior_source",
                  "tau_boat_prior_calibration_table", "tau_shore_prior_source", "opener_flags",
                  "opener_selected", "razor_dig_active", "marine_hazard_selected")

resolve_cfg <- function(sid) {
  cfg <- BASE
  for (k in names(WINDOW)) cfg[[k]] <- WINDOW[[k]]
  cfg$marine_hazard_manual_shore <- character(0); cfg$marine_hazard_manual_boat <- character(0)
  # 2026-09-27: run_config.R now ships the ADOPTED method (mode "manual", boat nws_sca_any,
  # marine_hazard_gear_regimes = "all_gear"). Every rung here rendered with a selected term
  # applied to BOTH sub-season fits (the M2 pot-closure boat fit carries the SCA term: 1,233
  # against M1's 1,372), so that is set here, NOT in WINDOW, to keep the rungs what they were
  # and their digests, and RESUME, untouched. The mode is each rung's own delta.
  cfg$marine_hazard_gear_regimes <- c("pot_closure", "all_gear")
  d <- STAGE_DEFS[[sid]]$delta
  for (k in names(d)) cfg[[k]] <- d[[k]]
  cfg
}
# A digest over the delta keys and the pinned window, so RESUME cannot reuse a folder whose
# configuration is not this rung's (the 2026-09-10 lesson: reusing on a folder NAME mixed
# configurations into one ladder).
digest_or_hash <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(x))
  v <- utf8ToInt(x); h <- 2166136261
  for (b in v) { h <- bitwXor(h, b %% 256); h <- (h * 16777619) %% 2^32 }
  sprintf("%08x%08x", h %% 2^32, (sum(v) * 2654435761) %% 2^32)
}
stage_digest <- function(sid) {
  cfg <- resolve_cfg(sid)
  keys <- sort(unique(c(declared_keys(sid), names(WINDOW))))
  txt <- paste(vapply(keys, function(k)
    paste0(k, "=", paste(format(unlist(cfg[[k]] %||% "NULL")), collapse = "|")), character(1)),
    collapse = ";")
  substr(digest_or_hash(txt), 1, 8)
}
.code_fingerprint <- function() {
  .g <- function(paths) {
    fs <- sort(unlist(lapply(paths, function(p)
      if (dir.exists(p)) list.files(p, pattern = "[.](R|Rmd|stan)$", full.names = TRUE) else p)))
    txt <- unlist(lapply(fs, function(f) {
      l <- readLines(f, warn = FALSE)
      l <- sub("#.*$", "", l); l <- sub("//.*$", "", l); l[nzchar(trimws(l))]
    }))
    substr(digest_or_hash(paste(txt, collapse = "\n")), 1, 8)
  }
  sprintf("stan:%s drivers:%s fns:%s",
          .g(.here("02_stan_models")), .g(c(POOLED_RMD)), .g(.here("03_R_functions")))
}
# FINGERPRINT PAIRS DECLARED EQUIVALENT for RESUME, "<recorded> => <current>" with the reason,
# the same audit-trail idea as run_improvements_2026-09-08.R's CODE_EQUIVALENT: a pair is
# only defensible when the change cannot reach a fit. The key names BOTH ends, so the
# declaration lapses the moment the tree moves again (the harness asserts the current side).
CODE_EQUIVALENT_MH <- list(
  # The 2026-09-25/26 ladder (five rungs stamped stan:65b5adeb drivers:ae200663 fns:094c314f)
  # against the tree after the results patch: fns moved by bss_block_cv.R (new; called only
  # from the drivers' per-fit diagnostics loop AFTER every fit) and by bss_stan_fit()
  # rebuilding fit@sim$permutation from the Stan seed after each fit (B38: the permutation
  # governs only the ORDER of extract()'s draws, so no posterior can change); drivers moved by the block-CV
  # call in the diagnostics loop and by the seeded census draw in the port block (B38),
  # both downstream of every fit. The stan layer did not move. On 2026-09-27 (the block-CV
  # results patch, Section 1y) fns moved again, all of it downstream of every fit:
  # bss_block_cv.R (the joint table, the identical-fit floor, the per-week helper),
  # write_loo_diagnostics() writing the OSP stream's pointwise LOO beside the others, and
  # .bma_core() keeping the adequacy aggregate on gear / trailer / catch so
  # model_adequacy.csv stays comparable. No likelihood, prior or Stan datum is touched.
  "stan:65b5adeb drivers:ae200663 fns:094c314f => stan:65b5adeb drivers:a65be4bc fns:8022c96c" =
    "B38 (the draw permutation rebuilt from the Stan seed after each fit, and the seeded census draw) and B39 (block-CV diagnostics) landed after the run; neither reaches a likelihood, a prior or the Stan data, so every per-fit summary the verdicts read is the run's own. The port median is the one number B38 would change on a re-render (once, by the MC jitter it removes). SUPERSEDED 2026-09-27 by the entry below (fns moved again, in post-fit diagnostics only).",
  "stan:65b5adeb drivers:ae200663 fns:094c314f => stan:65b5adeb drivers:73ee167d fns:fd9f0d80" =
    paste("As above, plus the 2026-09-27 block-CV results patch (Section 1y): bss_block_cv.R gained the joint effort table, the identical-fit floor and the per-week helper; write_loo_diagnostics() now writes loo_pointwise_osp_*.csv; .bma_core() filters the adequacy aggregate to gear / trailer / catch (all post-fit diagnostics);",
          "and B41 in bss_marine_hazard_covariates.R, a PREP-layer change: two further candidate columns (nws_sca_any_winter, nws_sca_any_rest) and the one-definition rule treating that pair as one. Offered, not selected, by any of these five rungs' configurations:",
          "measured on the real 2024-25 inputs, M2 to M5 resolve to the same per-population selection, the same per-date values of every selected column, the same screen table and the same selection table before and after the change, and M1 (off) returns before the module reads anything, so every rung's Stan data is unchanged.",
          "And the 2026-09-27 ADOPTION (A30, B42): marine_hazard_terms_for() in the module and one call in each prep, confining a selected term to the sub-seasons in marine_hazard_gear_regimes; this runner sets that key to BOTH regimes in resolve_cfg(), which is what every rung here did,",
          "so the preps build the same X_open for M1 to M6 (checked on the real inputs: the boat pot-closure prep under both regimes carries the term as M2's did); the drivers moved by report prose and comments only (the marine_hazard_prepare() call is unchanged).",
          "And the 2026-09-28 render review (B43), after every fit: the fitted day-covariate table in both reports (bss_day_covariate_report.R), the census draws carried into the pooled season totals, and a seeded draw subsample in write_effort_overdispersion_diag().",
          "SUPERSEDED 2026-09-28 by B44, which is NOT inference-equivalent (see CODE_NOT_EQUIVALENT_MH below)."),
  # M6 rendered on 2026-09-27 from the tree at 4e23b15 (the block-CV results patch applied; its stamp
  # is that commit's fingerprint, recomputed and matched), so its recorded side is the second entry's
  # current side. The only fitting-layer change since is B42, covered by the last sentence above.
  "stan:65b5adeb drivers:a65be4bc fns:6f65d84a => stan:65b5adeb drivers:73ee167d fns:fd9f0d80" =
    paste("M6 rendered from the tree with B40 and B41 applied (4e23b15). Since then only B42 (A30, the 2026-09-27 adoption) touched the fitting layer: marine_hazard_terms_for() in the module and one call in each prep,",
          "confining a selected term to the sub-seasons in marine_hazard_gear_regimes; this runner pins that key to BOTH regimes in resolve_cfg(), under which the function returns the selection unchanged for every fit, so M6's Stan data is what it was",
          "(the M0 desk row 'M6 vs M2 differs in K_open and X_open_flat only, columns sum' is recomputed on the current tree at every run). The pooled driver moved by the section 3.7 prose only; the marine_hazard_prepare() call is unchanged.",
          "B43 (2026-09-28) is post-fit reporting and diagnostics only (the day-covariate table, the season totals' census draws, the overdispersion subsample's seed).",
          "SUPERSEDED 2026-09-28 by B44, which is NOT inference-equivalent (see CODE_NOT_EQUIVALENT_MH below).")
)
# EXAMINED AND DECLARED NOT EQUIVALENT (2026-09-28; B51 2026-09-29), as run_improvements_2026-09-08.R's
# CODE_NOT_EQUIVALENT: the current tree was checked against these rungs and B44 changes the
# shore all-gear Stan data, so no equivalence is written above. Record only; RESUME reports
# the code delta as REVIEW, as it did before this list existed.
CODE_NOT_EQUIVALENT_MH <- list(
  "stan:65b5adeb drivers:ae200663 fns:094c314f => stan:73a6cf27 drivers:be33b901 fns:a41ba28d" = "B44 (2026-09-28) is NOT inference-equivalent for these rungs, and is recorded here rather than declared equivalent. repair_interview_ids() in fetch_crab_data() gives each interview its own id: on 2024-25 one id (S6070_29, 2025-05-24, shore, Float 20) was shared by two different interviews, so the shore all-gear CPUE likelihood now carries both (catch 0 and 2) where distinct(interview_id) kept one carrying the pair's summed catch, and the shore all-gear PE catch falls 16 crab (29,737 to 29,721). Every rung's shore all-gear fit therefore differs, slightly, from a re-render at this tree: a RESUME still reuses the folders by digest and reports the code delta as REVIEW, which is the correct reading. The rest of B44 is post-fit or inert under these configurations: the PE monthly split (pe_monthly_split()), the PE effort SE's donor covariances, the gear driver's expected-catch totals (gear track only), and the removal of the estimate_red_rock switch (FALSE in every pin). And B45 (2026-09-28), which changes MORE: shore_dock_counts() fills a Float 20 count with no Float 17-21 count beside it by round(R_month x that Float 20 count) where it used 0 (shore_f17_fill = 'ratio'), which raises every shore fit's gear counts (2024-25: the mean daily shore gear count 40.9 to 44.6, +9.0%) and the shore PE with them; and the gear prep no longer drops interviews without a positive fishing time (none in 2024-25). And B46 (2026-09-28), which also reaches fits: prep_days_crab() keys the PE period and week_index by ISO week and ISO year (%V, %G) where it used %W, which on 2024-25 joins the week of 30 December 2024 to 5 January 2025 (a 2-day and a 5-day period before), so a fit at weekly AR resolution (the shore all-gear fit, where a rung fits it weekly) has one fewer AR period (43, not 44) and the PE port moves 88,819 to 88,758; the gear Stan's R_G prior centre is the season's interview ratio rather than the literal 1.3 (gear track only); both Stan files reject osp_crab_lower = 1 with a dynamic f and no dynamic c (unreachable under these pins). The rest of B46 is post-fit, reporting or infrastructure: Predictive_Catch rebuilt from the fitted ZINB (Expected_Catch unchanged), bss_with_seed() restoring the caller's RNG, the gate's NA verdict, the output folder, the package loader. And B48 (2026-09-28), inert under these pins: fetch_osp_boat_counts() resolves the OSP crab-only binomial n per day from the boats SAMPLED (osp_sampling_rates.R, new), and both drivers write osp_crab_only_daily.csv; with no crab-only column in WBL_boat_counts.xlsx the reader returns no crab rows and the same effort series, so no fit changes. And B51 (2026-09-29), NOT inference-equivalent for any rung either: every fit now starts its chains within init_r = 0.5 on the unconstrained scale (run_config$bss_init_r, passed by both drivers to bss_stan_fit()) where rstan's default radius of 2 applied, so a re-render at this tree changes every fit's draws. The posterior targeted is the same; the draws are not, and on a fit bound by a funnel (the shore all-gear fit, D33) the divergence count and the gate verdict can move with them. bss_sampler_override() now also accepts bss_init_r. And A31 (2026-09-29), NOT inference-equivalent: the pooled Stan collapses the single-section effort-level hierarchy (mu_hier_collapse_single = 'effort'; 'none' restores the old model exactly), which removes the unidentified sigma_mu_E * eps_mu_E pair from every S == 1 fit, i.e. every rung's four fits; the posterior of the reported totals moves by a small fraction of a posterior SD on the container refits, but no rung's draws are reproducible at this tree. B52 to B56 (the same day) are inert under these pins or post-fit: the census clipped to the estimation window (every pin's window contains its census), the turnover fallbacks (not reached: every pin resolves both centres), the day length computed only under a time unit (no pin uses one), run_pe_gear parity (inert on 2024-25), and the report.",
  "stan:65b5adeb drivers:a65be4bc fns:6f65d84a => stan:73a6cf27 drivers:be33b901 fns:a41ba28d" = "M6 as the five rungs: B42 was the only fitting-layer change between M6's tree and the adoption; B44 (2026-09-28) is NOT inference-equivalent for these rungs, and is recorded here rather than declared equivalent. repair_interview_ids() in fetch_crab_data() gives each interview its own id: on 2024-25 one id (S6070_29, 2025-05-24, shore, Float 20) was shared by two different interviews, so the shore all-gear CPUE likelihood now carries both (catch 0 and 2) where distinct(interview_id) kept one carrying the pair's summed catch, and the shore all-gear PE catch falls 16 crab (29,737 to 29,721). Every rung's shore all-gear fit therefore differs, slightly, from a re-render at this tree: a RESUME still reuses the folders by digest and reports the code delta as REVIEW, which is the correct reading. The rest of B44 is post-fit or inert under these configurations: the PE monthly split (pe_monthly_split()), the PE effort SE's donor covariances, the gear driver's expected-catch totals (gear track only), and the removal of the estimate_red_rock switch (FALSE in every pin). And B45 (2026-09-28), which changes MORE: shore_dock_counts() fills a Float 20 count with no Float 17-21 count beside it by round(R_month x that Float 20 count) where it used 0 (shore_f17_fill = 'ratio'), which raises every shore fit's gear counts (2024-25: the mean daily shore gear count 40.9 to 44.6, +9.0%) and the shore PE with them; and the gear prep no longer drops interviews without a positive fishing time (none in 2024-25). And B46 (2026-09-28), which also reaches fits: prep_days_crab() keys the PE period and week_index by ISO week and ISO year (%V, %G) where it used %W, which on 2024-25 joins the week of 30 December 2024 to 5 January 2025 (a 2-day and a 5-day period before), so a fit at weekly AR resolution (the shore all-gear fit, where a rung fits it weekly) has one fewer AR period (43, not 44) and the PE port moves 88,819 to 88,758; the gear Stan's R_G prior centre is the season's interview ratio rather than the literal 1.3 (gear track only); both Stan files reject osp_crab_lower = 1 with a dynamic f and no dynamic c (unreachable under these pins). The rest of B46 is post-fit, reporting or infrastructure: Predictive_Catch rebuilt from the fitted ZINB (Expected_Catch unchanged), bss_with_seed() restoring the caller's RNG, the gate's NA verdict, the output folder, the package loader. And B48 (2026-09-28), inert under these pins: fetch_osp_boat_counts() resolves the OSP crab-only binomial n per day from the boats SAMPLED (osp_sampling_rates.R, new), and both drivers write osp_crab_only_daily.csv; with no crab-only column in WBL_boat_counts.xlsx the reader returns no crab rows and the same effort series, so no fit changes. And B51 (2026-09-29), NOT inference-equivalent for any rung either: every fit now starts its chains within init_r = 0.5 on the unconstrained scale (run_config$bss_init_r, passed by both drivers to bss_stan_fit()) where rstan's default radius of 2 applied, so a re-render at this tree changes every fit's draws. The posterior targeted is the same; the draws are not, and on a fit bound by a funnel (the shore all-gear fit, D33) the divergence count and the gate verdict can move with them. bss_sampler_override() now also accepts bss_init_r. And A31 (2026-09-29), NOT inference-equivalent: the pooled Stan collapses the single-section effort-level hierarchy (mu_hier_collapse_single = 'effort'; 'none' restores the old model exactly), which removes the unidentified sigma_mu_E * eps_mu_E pair from every S == 1 fit, i.e. every rung's four fits; the posterior of the reported totals moves by a small fraction of a posterior SD on the container refits, but no rung's draws are reproducible at this tree. B52 to B56 (the same day) are inert under these pins or post-fit: the census clipped to the estimation window (every pin's window contains its census), the turnover fallbacks (not reached: every pin resolves both centres), the day length computed only under a time unit (no pin uses one), run_pe_gear parity (inert on 2024-25), and the report."
)
# 2026-09-28 (B46): the input workbooks and run_config.R's VALUES, which neither the stage digest
# (pinned keys only) nor the code fingerprint covers; stamped and checked on RESUME like the code
# line (reported, a REVIEW row, not enforced), as in run_improvements_2026-09-08.R.
.inputs_fingerprint <- function() {
  fs <- sort(list.files(.here("04_input_files"), pattern = "[.]xlsx$", full.names = TRUE, recursive = TRUE))
  rc <- tryCatch({ e <- new.env(); sys.source(.here("run_config.R"), envir = e)
                   paste(deparse(e$run_config[sort(names(e$run_config))]), collapse = "") },
                 error = function(e) "run_config.R unreadable")
  substr(digest_or_hash(paste(c(paste(basename(fs), unname(tools::md5sum(fs))), rc), collapse = "\n")), 1, 8)
}
.stage_stamp <- function(dir, sid) {
  writeLines(c(sprintf("stage: %s", sid),
               sprintf("digest: %s", stage_digest(sid)),
               sprintf("code: %s", .code_fingerprint()),
               sprintf("inputs: %s", .inputs_fingerprint()),
               sprintf("rstan: %s / StanHeaders %s",
                       as.character(utils::packageVersion("rstan")),
                       as.character(utils::packageVersion("StanHeaders"))),
               sprintf("written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               "", "# run_marine_hazard_batch_2026-09-25.R writes this after a rung renders. RESUME",
               "# reuses a folder ONLY when its digest matches the rung being requested."),
             file.path(dir, "MH_STAGE.txt"))
}
.stamp_field <- function(dir, key) {
  p <- file.path(dir %||% "", "MH_STAGE.txt")
  if (!file.exists(p)) return(NA_character_)
  l <- grep(paste0("^", key, ":"), readLines(p, warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_character_) else trimws(sub(paste0("^", key, ":"), "", l[1]))
}
find_outdir <- function(tag) {
  dirs <- list.dirs(.here("05_output"), recursive = FALSE)
  hits <- unlist(lapply(dirs, function(d)
    list.dirs(d, recursive = FALSE)[basename(list.dirs(d, recursive = FALSE)) == paste0(PREFIX, tag)]))
  if (!length(hits)) return(NA_character_)
  hits[order(basename(dirname(hits)), decreasing = TRUE)][1]
}

# ---------------------------------------------------------------------------
# M0: THE DESK STAGE. Everything that can be proven without MCMC is proven here.
# ---------------------------------------------------------------------------
stage_M0 <- function() {
  banner("M0  DESK: prerequisites, before any MCMC")
  q <- function(e) { s <- tempfile(); sink(s); on.exit(sink()); force(e) }

  # (1) the module and its entry points exist; the preps read the selection
  fns <- c("marine_hazard_events", "marine_hazard_flag_series", "bar_restriction_series",
           "marine_hazard_screen", "marine_hazard_select", "marine_hazard_prepare",
           "marine_hazard_terms_for")   # B42 (2026-09-27): the per-fit gate on the selection
  have <- vapply(fns, exists, logical(1))
  V1row("M0", "the marine hazard module is sourced and complete",
        sprintf("%d of %d functions present", sum(have), length(fns)), "all present",
        if (all(have)) "PASS" else "FAIL", "Nothing below can run without them.")
  for (f in c(POOLED_PREP, .here("03_R_functions", "prep_bss_crab_gear.R"))) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    # until 2026-09-27 the preps read params$marine_hazard_selected directly; since B42 they take the fit's
    # terms from marine_hazard_terms_for(params, population_name, gear_regime), which returns the population's
    # selection only for the sub-season gear regimes in marine_hazard_gear_regimes (this runner: both)
    wired <- grepl("marine_extra <- marine_hazard_terms_for(params, population_name, gear_regime)", src, fixed = TRUE) &&
             grepl("c(razor_extra, marine_extra)", src, fixed = TRUE)
    V1row("M0", sprintf("%s passes the marine selection into the K_open block", basename(f)),
          if (wired) "yes" else "NO",
          "takes marine_hazard_terms_for(params, population_name, gear_regime) and adds it to opener_design_matrix()'s extra",
          if (wired) "PASS" else "FAIL",
          "The driver installs the selection; the prep is where it becomes a column, for the fits whose gear regime is configured. Both tracks.")
  }
  st <- readLines(POOLED_STAN, warn = FALSE)
  V1row("M0", "the pooled Stan carries the K_open / X_open_flat / B_open block (unchanged)",
        paste(c("K_open", "X_open_flat", "B_open", "B_open_out")[vapply(c("K_open", "X_open_flat", "B_open", "B_open_out"),
              function(v) any(grepl(paste0("(^|[^A-Za-z0-9_])", v, "([^A-Za-z0-9_]|$)"), st)), logical(1))], collapse = ", "),
        "all four", if (all(vapply(c("K_open", "X_open_flat", "B_open", "B_open_out"),
                                  function(v) any(grepl(paste0("(^|[^A-Za-z0-9_])", v, "([^A-Za-z0-9_]|$)"), st)), logical(1)))) "PASS" else "FAIL",
        "The covariates are columns of this block; no Stan edit was made for them.")

  # (2) the shipped configuration. Until 2026-09-27 this row asserted "off" (the machinery
  #     BUILT, INERT until the run and its review said otherwise). The review said: adopt the
  #     constant boat term, all-gear fit only (Section 1z; A30 ADOPTED), so the row now
  #     REPORTS what run_config ships and PASSES on either the pre-adoption "off" or the
  #     adopted method; anything else is a configuration this ladder did not test.
  .ship_ok <- identical(BASE$marine_hazard_mode, "off") ||
    (identical(BASE$marine_hazard_mode, "manual") && identical(BASE$marine_hazard_manual_boat, "nws_sca_any") &&
     !length(BASE$marine_hazard_manual_shore %||% character(0)) && identical(BASE$marine_hazard_gear_regimes, "all_gear"))
  V1row("M0", "the SHIPPED config is the pre-adoption model (off) or the adopted method (manual: boat nws_sca_any, all-gear fit only, no shore term)",
        sprintf("marine_hazard_mode = %s; manual_boat = %s; manual_shore = %s; gear_regimes = %s", BASE$marine_hazard_mode %||% "NULL",
                paste(BASE$marine_hazard_manual_boat %||% character(0), collapse = "+"), paste(BASE$marine_hazard_manual_shore %||% character(0), collapse = "+"),
                paste(BASE$marine_hazard_gear_regimes %||% character(0), collapse = "+")),
        "off, or manual / nws_sca_any / none / all_gear", if (.ship_ok) "PASS" else "FAIL",
        "The rungs here set their own marine keys (resolve_cfg), so run_config's choice does not enter any rung; this row only says what production would do.")

  # (3) the archive covers the window; the flags and the screen build on the real data
  ok <- tryCatch({
    p <- resolve_cfg("M5"); p <- modifyList(p, list(bss_model_file = "crab_bss_pooled.stan", boat_require_gear_time = TRUE))
    p$ar_max_resolution <- p$ar_max_resolution$pooled
    p$crabbing_holiday_dates <- read_crabbing_holidays(p)
    ev <- marine_hazard_events(p); cov <- attr(ev, "coverage")
    V1row("M0", "the NWS archive covers the window for both zones",
          paste(sprintf("%s %s to %s", cov$ugc, cov$pull_start, cov$pull_end), collapse = "; "),
          sprintf("%s to %s", WINDOW$est_date_start, WINDOW$est_date_end),
          if (all(cov$pull_start <= as.Date(WINDOW$est_date_start)) && all(cov$pull_end >= as.Date(WINDOW$est_date_end)) &&
              all(c(WINDOW$marine_hazard_zones) %in% cov$ugc)) "PASS" else "FAIL",
          paste("marine_hazard_flag_series() stops on an uncovered window, so without this every",
                "covariate rung would error at the driver. The committed workbook is a live pull (2026-09-25)",
                "of the IEM archive (see build_nws_marine_hazards.R); re-pull before a window past its pull_end."))
    dwg <- q(fetch_crab_data(p))
    mh <- q(marine_hazard_prepare(dwg, p, output_dir = NULL, quiet = TRUE))
    fl <- mh$flags
    cat(sprintf("  definition: %s\n", mh$definition))
    cat(sprintf("  flag days in the window: nws_sca_any %d, nws_sca_bar %d, nws_sca_coastal %d of %d; bar restriction observed on %d days (%d restricted), imputed on %d\n",
                sum(fl$nws_sca_any), sum(fl$nws_sca_bar), sum(fl$nws_sca_coastal), nrow(fl),
                sum(!is.na(fl$bar_restriction_obs)), sum(fl$bar_restriction_obs %in% 1), sum(is.na(fl$bar_restriction_obs))))
    if (!is.null(mh$screen)) { cat("  the screen on this window:\n"); print(as.data.frame(mh$screen[, c("population", "covariate", "n_days", "n_flag", "rate_ratio", "adj_p")]), row.names = FALSE) }
    if (!is.null(mh$sel$table)) { cat("  auto selection:\n"); print(as.data.frame(mh$sel$table[, c("population", "covariate", "rate_ratio", "p_adj", "selected")]), row.names = FALSE) }
    V1row("M0", "the flags, the bar imputation and the screen build on the window's data",
          sprintf("SCA on %d of %d days; bar observed %d / imputed %d; screen rows %d; auto selects shore {%s} boat {%s}",
                  sum(fl$nws_sca_any), nrow(fl), sum(!is.na(fl$bar_restriction_obs)), sum(is.na(fl$bar_restriction_obs)),
                  if (is.null(mh$screen)) 0L else nrow(mh$screen),
                  paste(mh$sel$shore, collapse = ","), paste(mh$sel$private_boat, collapse = ",")),
          "no error; every flag in [0, 1]",
          if (all(fl$nws_sca_any %in% c(0, 1)) && all(is.finite(fl$bar_restriction)) &&
              all(fl$bar_restriction >= 0 & fl$bar_restriction <= 1)) "PASS" else "FAIL",
          "The M5 rung reads exactly this selection; the screen table is what its report will show.")
    V1row("M0", "REPORTED: what the auto screen selects on 2024-25",
          sprintf("shore: %s | boat: %s", if (length(mh$sel$shore)) paste(mh$sel$shore, collapse = ", ") else "(none)",
                  if (length(mh$sel$private_boat)) paste(mh$sel$private_boat, collapse = ", ") else "(none)"),
          "no threshold", "INFO",
          "Offline (2026-09-25) the boat SCA and bar terms cleared the BH screen and the shore SCA did not.")

    # (4) THE STAN DATA PROOF. M4 against M1, built for real for the two all-gear fits: the
    #     only entries that may differ are K_open and X_open_flat.
    ie <- q(fetch_ie_data(p))
    p$crab_fraction_rows <- q(crab_fraction_source_rows(dwg, ie, p))
    stt <- q(estimate_shore_turnover(attr(ie, "ie_intervals"), dwg$shore_effort, p))
    osp <- tryCatch(q(fetch_osp_boat_counts(p)), error = function(e) NULL)
    ov  <- tryCatch(q(diagnose_osp_trailer_overlap(osp, p, output_dir = NULL)), error = function(e) NULL)
    p$osp_crab_rows <- attr(osp, "osp_crab_rows")
    p <- q(bss_resolve_tau_boat_prior(p, ov)); p <- q(bss_resolve_tau_shore_prior(p, stt))
    Le <- if (isTRUE(p$use_ie_day_length) && nrow(ie) > 0) tryCatch(q(estimate_L_effective(ie, p)), error = function(e) NULL) else NULL
    sub <- build_subseasons(p)
    ss  <- sub[[which(vapply(sub, function(x) x$gear_regime == "all_gear", logical(1)))[1]]]
    days <- q(prep_days_crab(ss$start, ss$end, p, L_eff_model = Le))
    build <- function(sid, pop) {
      pp <- p; d <- STAGE_DEFS[[sid]]$delta; for (k in names(d)) pp[[k]] <- d[[k]]
      pp$marine_hazard_manual_shore <- d$marine_hazard_manual_shore %||% character(0)
      pp$marine_hazard_manual_boat  <- d$marine_hazard_manual_boat  %||% character(0)
      pp <- q(marine_hazard_prepare(dwg, pp, output_dir = NULL, quiet = TRUE))$params
      summ <- q(prep_population_summary(dwg, pop, ss$start, ss$end, pp))
      q(prep_bss_crab_pooled(days, summ, "Dungeness_Kept", pp, pop, gear_regime = ss$gear_regime, ie_data = ie))
    }
    expect <- list(shore = "nws_sca_any", private_boat = c("nws_sca_any", "bar_restriction"))
    for (pop in c("shore", "private_boat")) {
      a <- build("M1", pop); b <- build("M4", pop)
      keys <- union(names(a), names(b))
      diff <- keys[!vapply(keys, function(k) isTRUE(all.equal(a[[k]], b[[k]], tolerance = 0)), logical(1))]
      lab_b <- attr(b, "opener_labels")
      xb <- matrix(b$X_open_flat, nrow = b$D)
      ok_cols <- b$K_open == length(expect[[pop]]) && identical(lab_b, expect[[pop]])
      V1row("M0", sprintf("%s all-gear: M4's Stan data differs from M1's in K_open and X_open_flat ONLY", pop),
            sprintf("M1 K_open %d; M4 K_open %d [%s]; entries differing: %s", a$K_open, b$K_open,
                    paste(lab_b, collapse = ","), if (length(diff)) paste(diff, collapse = ", ") else "NONE"),
            "differing entries = {K_open, X_open_flat}; M1 K_open = 0; M4 labels as declared",
            if (setequal(diff, c("K_open", "X_open_flat")) && a$K_open == 0 && ok_cols) "PASS" else "FAIL",
            paste("This is the whole basis of the covariate being a covariate and nothing else: the",
                  "day set, the AR period, the counts, the interviews, the priors and the f data are",
                  "byte-identical between the two rungs, and the term enters through one block."))
      if (b$K_open > 0) {
        in01 <- all(xb >= 0 & xb <= 1)
        frac_days <- if ("bar_restriction" %in% lab_b) days$event_date[xb[, match("bar_restriction", lab_b)] %% 1 != 0] else as.Date(character(0))
        unsampled <- fl$event_date[is.na(fl$bar_restriction_obs)]
        frac_ok <- all(frac_days %in% unsampled)
        V1row("M0", sprintf("%s all-gear: the X_open columns are flags in [0, 1]; the bar column carries imputed days", pop),
              sprintf("D = %d; column sums %s; bar days with a fractional (imputed) value %d, all of them unsampled: %s", b$D,
                      paste(sprintf("%s = %.1f", lab_b, colSums(xb)), collapse = ", "), length(frac_days), frac_ok),
              "all in [0, 1]; every fractional bar value on a day with no observed tick",
              if (in01 && frac_ok) "PASS" else "FAIL",
              "An imputed expected value is the documented approximation (bss_marine_hazard_covariates.R header).")
      }
    }
    # (4b) B41: M6 against M2 for the boat all-gear fit. The split must enter as two flag
    #      columns that SUM to M2's one, and nothing else may differ.
    if ("M6" %in% names(STAGE_DEFS)) {
      a2 <- build("M2", "private_boat"); b6 <- build("M6", "private_boat")
      keys <- union(names(a2), names(b6))
      diff <- keys[!vapply(keys, function(k) isTRUE(all.equal(a2[[k]], b6[[k]], tolerance = 0)), logical(1))]
      lab6 <- attr(b6, "opener_labels"); x6 <- matrix(b6$X_open_flat, nrow = b6$D); x2 <- matrix(a2$X_open_flat, nrow = a2$D)
      sums_ok <- b6$K_open == 2L && identical(lab6, c("nws_sca_any_winter", "nws_sca_any_rest")) &&
                 a2$K_open == 1L && isTRUE(all.equal(rowSums(x6), x2[, 1], tolerance = 0))
      V1row("M0", "private_boat all-gear: M6's Stan data differs from M2's in K_open and X_open_flat ONLY, and its two columns sum to M2's one",
            sprintf("M2 K_open %d; M6 K_open %d [%s]; column sums %s; entries differing: %s", a2$K_open, b6$K_open, paste(lab6, collapse = ","),
                    paste(sprintf("%s = %.0f", lab6, colSums(x6)), collapse = ", "), if (length(diff)) paste(diff, collapse = ", ") else "NONE"),
            "differing entries = {K_open, X_open_flat}; K_open 1 -> 2; winter + rest = nws_sca_any on every day",
            if (setequal(diff, c("K_open", "X_open_flat")) && sums_ok) "PASS" else "FAIL",
            "The season split is the same flag in two columns; if it were anything else the rung would measure that too (rule 9).")
    }
    # (4c) THE SHIPPED CONFIGURATION, fit by fit, on the real inputs (2026-09-27, A30 ADOPTED / B42).
    #      This ladder pins marine_hazard_gear_regimes to BOTH regimes (its rungs applied a term to
    #      both fits); production restricts it. So here run_config.R's own marine keys are laid over
    #      the resolved window and all four pooled fits are built: the adopted method puts the term
    #      in the boat all-gear fit and nowhere else; the pre-adoption "off" puts it nowhere.
    ps <- p
    for (k in grep("^marine_hazard_|^bar_restriction_", names(BASE), value = TRUE)) ps[[k]] <- BASE[[k]]
    for (k in setdiff(grep("^marine_hazard_", names(ps), value = TRUE), names(BASE))) ps[[k]] <- NULL
    ps <- q(marine_hazard_prepare(dwg, ps, output_dir = NULL, quiet = TRUE))$params
    got <- list(); notes <- character(0)
    for (ss_ in sub) for (pop in c("shore", "private_boat")) {
      d_ <- q(prep_days_crab(ss_$start, ss_$end, ps, L_eff_model = Le))
      sm <- q(prep_population_summary(dwg, pop, ss_$start, ss_$end, ps))
      out <- utils::capture.output(bd <- prep_bss_crab_pooled(d_, sm, "Dungeness_Kept", ps, pop, gear_regime = ss_$gear_regime, ie_data = ie))
      notes <- c(notes, grep("NOT applied to this fit", out, value = TRUE))
      lab <- attr(bd, "opener_labels") %||% character(0)
      got[[sprintf("%s/%s", pop, ss_$gear_regime)]] <- sprintf("K_open %d%s", bd$K_open,
        if (bd$K_open > 0) sprintf(" [%s] flagged %d of %d days", paste(lab, collapse = ","), sum(matrix(bd$X_open_flat, nrow = bd$D)[, 1] > 0), bd$D) else "")
    }
    k1 <- vapply(got, function(x) as.integer(sub("^K_open (\\d+).*$", "\\1", x)), integer(1))
    adopted <- identical(BASE$marine_hazard_mode, "manual") && identical(BASE$marine_hazard_manual_boat, "nws_sca_any") &&
               !length(BASE$marine_hazard_manual_shore %||% character(0)) && identical(BASE$marine_hazard_gear_regimes, "all_gear")
    want <- if (adopted) c("shore/pot_closure" = 0L, "shore/all_gear" = 0L, "private_boat/pot_closure" = 0L, "private_boat/all_gear" = 1L)
            else if (identical(BASE$marine_hazard_mode, "off")) c("shore/pot_closure" = 0L, "shore/all_gear" = 0L, "private_boat/pot_closure" = 0L, "private_boat/all_gear" = 0L)
            else NULL
    ok4 <- !is.null(want) && all(names(want) %in% names(k1)) && identical(unname(k1[names(want)]), unname(want)) &&
           (!adopted || (grepl("[nws_sca_any]", got[["private_boat/all_gear"]], fixed = TRUE) && length(notes) == 1L && grepl("pot_closure", notes[1], fixed = TRUE)))
    V1row("M0", "the SHIPPED configuration on the real inputs, fit by fit: the term enters the boat all-gear fit and no other",
          paste(c(sprintf("%s: %s", names(got), unlist(got)), sprintf("withheld-notes printed: %d", length(notes))), collapse = "; "),
          if (adopted) "shore 0 / 0; boat pot-closure 0 (one note); boat all-gear 1 [nws_sca_any]" else if (is.null(want)) "no expectation for this configuration" else "K_open 0 in all four fits (off)",
          if (is.null(want)) "INFO" else if (ok4) "PASS" else "FAIL",
          paste("marine_hazard_terms_for() is the only thing between the per-population selection and the fit; this row is the",
                "production path (drivers pass ss$gear_regime) built on the season's data, not a unit test on a synthetic params."))
    TRUE
  }, error = function(e) {
    V1row("M0", "the flags, the screen and the Stan data proof were evaluated", conditionMessage(e), "no error", "ERROR",
          "Needs the data-reading packages and the input files; it is not a Stan check.")
    FALSE
  })

  # (5) COMPARABILITY. Each rung's resolved configuration must differ from M1's in declared keys only.
  base_cfg <- resolve_cfg("M1")
  for (sid in setdiff(names(STAGE_DEFS), c("M0", "M1"))) {
    cfg <- resolve_cfg(sid)
    ks  <- setdiff(union(names(base_cfg), names(cfg)), RUNTIME_KEYS)
    diff <- ks[!vapply(ks, function(k)
      identical(format(unlist(cfg[[k]] %||% "NULL")), format(unlist(base_cfg[[k]] %||% "NULL"))), logical(1))]
    undeclared <- setdiff(diff, declared_keys(sid))
    V1row(sid, "differs from M1 in DECLARED keys only",
          sprintf("differs in: %s", if (length(diff)) paste(diff, collapse = ", ") else "nothing"),
          sprintf("a subset of {%s}", paste(declared_keys(sid), collapse = ", ")),
          if (!length(undeclared)) "PASS" else "FAIL",
          paste("A rung that differs in an undeclared key measures that key as well as its own.",
                "Undeclared here:", if (length(undeclared)) paste(undeclared, collapse = ", ") else "none"))
  }
  # (6) the reference render for the inertness verdict
  V1row("M0", sprintf("the reference render %s is present", basename(REF_R4)),
        if (dir.exists(.here("05_output", REF_R4))) "present" else "MISSING", "present",
        if (dir.exists(.here("05_output", REF_R4))) "PASS" else "FAIL",
        "M1's bit-identity verdict reads it. Without it the inertness proof is structural only (M0's Stan-data row).")
  invisible(ok)
}

# ---------------------------------------------------------------------------
# RUNNING ONE RUNG
# ---------------------------------------------------------------------------
run_one <- function(sid) {
  st <- STAGE_DEFS[[sid]]
  rule()
  cat(sprintf("  %-3s %-12s  %s\n", sid, st$tag, st$item))
  cfg <- resolve_cfg(sid)
  cat(sprintf("       marine_hazard_mode = %-7s shore = {%s}  boat = {%s}  digest %s\n",
              cfg$marine_hazard_mode, paste(cfg$marine_hazard_manual_shore, collapse = ","),
              paste(cfg$marine_hazard_manual_boat, collapse = ","), stage_digest(sid)))
  existing <- find_outdir(st$tag)
  if (isTRUE(RESUME) && !is.na(existing) && file.exists(file.path(existing, "run_parameters.txt"))) {
    dg <- .stamp_field(existing, "digest")
    if (identical(dg, stage_digest(sid))) {
      cd <- .stamp_field(existing, "code")
      cat("  RESUME: output present at", basename(existing), "with a MATCHING digest - skipping the fit.\n")
      .inp <- .stamp_field(existing, "inputs")
      if (!is.na(.inp) && !identical(.inp, .inputs_fingerprint())) {
        cat("          *** THE INPUT WORKBOOKS OR run_config.R HAVE CHANGED SINCE THAT FIT; reused, cross-rung claims downgraded. ***\n")
        V1row(sid, "the reused fit was produced from DIFFERENT inputs than this run",
              sprintf("recorded %s; now %s", .inp, .inputs_fingerprint()), "the fingerprints match", "REVIEW",
              "RESUME matched the CONFIG digest (pinned keys); an input workbook or an unpinned run_config.R key changed since.")
      }
      if (!is.na(cd) && !identical(cd, .code_fingerprint())) {
        eq <- CODE_EQUIVALENT_MH[[paste(cd, .code_fingerprint(), sep = " => ")]]
        if (!is.null(eq)) {
          cat(sprintf("          code moved since that fit (%s -> %s), DECLARED EQUIVALENT: %s\n", cd, .code_fingerprint(), eq))
          V1row(sid, "the reused fit was produced by code declared EQUIVALENT to this run's",
                sprintf("recorded %s; now %s", cd, .code_fingerprint()), "a declared pair", "INFO", eq)
        } else {
          cat(sprintf("          *** THE CODE HAS CHANGED SINCE THAT FIT.\n              recorded %s\n              now      %s\n", cd, .code_fingerprint()))
          V1row(sid, "the reused fit was produced by DIFFERENT code than this run",
                sprintf("recorded %s; now %s", cd, .code_fingerprint()), "the fingerprints match", "REVIEW",
                "RESUME matched the CONFIG digest; any cross-rung claim involving this rung is indicative only.")
        }
      }
      return(existing)
    }
    cat(sprintf("  RESUME: %s exists but its digest %s does not match (%s). RE-RUNNING.\n",
                basename(existing), if (is.na(dg)) "is ABSENT" else dg, stage_digest(sid)))
  }
  if (isTRUE(DRY_RUN)) { cat("  DRY_RUN: not fitting.\n"); return(NA_character_) }
  cfg$model <- "pooled"; cfg$run_tag <- st$tag
  run_env <- new.env(parent = globalenv()); run_env$run_config <- cfg
  t0 <- Sys.time()
  html <- rmarkdown::render(POOLED_RMD, envir = run_env, quiet = FALSE)
  od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e) NA_character_)
  if (!is.na(od) && dir.exists(od) && file.exists(html) && normalizePath(dirname(html)) != normalizePath(od)) {
    if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE))) suppressWarnings(file.remove(html))
    else cat("  WARNING: could not move the rendered HTML; the NEXT rung will overwrite it.\n")
  }
  done <- find_outdir(st$tag)
  cat(sprintf("  %s finished in %.1f min -> %s\n", sid, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              if (is.na(done)) "OUTPUT FOLDER NOT FOUND" else basename(done)))
  if (is.na(done)) stop(sprintf("%s rendered but no folder named %s%s exists under 05_output", sid, PREFIX, st$tag))
  .stage_stamp(done, sid)
  done
}

# ---------------------------------------------------------------------------
# THE LADDER TABLE. One row per rung with every column the rule reads.
# ---------------------------------------------------------------------------
LAD <- list()
ladder_row <- function(sid, dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  pt <- .port_row(rd(dir, "port_total_Dungeness_Kept.csv")); gt <- .gate(dir)
  aqs <- .adq(dir, FIT_SHORE); aqb <- .adq(dir, FIT_BOAT)
  als <- .arlog(dir, FIT_SHORE); alb <- .arlog(dir, FIT_BOAT)
  bs <- .b_open(dir, FIT_SHORE, "nws_sca_any"); bb <- .b_open(dir, FIT_BOAT, "nws_sca_any"); bbar <- .b_open(dir, FIT_BOAT, "bar_restriction")
  bw <- .b_open(dir, FIT_BOAT, "nws_sca_any_winter"); br <- .b_open(dir, FIT_BOAT, "nws_sca_any_rest")   # B41 (M6)
  ls <- .loo_summ(dir, FIT_SHORE, "gear"); lb <- .loo_summ(dir, FIT_BOAT, "trailer")
  lsc <- .loo_summ(dir, FIT_SHORE, "catch"); lbc <- .loo_summ(dir, FIT_BOAT, "catch")
  nobs <- function(fit) { d <- rd(dir, sprintf("ppc_byobs_%s.csv", fit)); if (is.null(d)) NA_integer_ else sum(d$data_type == "catch", na.rm = TRUE) }
  LAD[[sid]] <<- data.frame(
    rung = sid, folder = basename(dir), mode = resolve_cfg(sid)$marine_hazard_mode,
    shore_cov = paste(.cov_labels(dir, FIT_SHORE) %||% character(0), collapse = "+"),
    boat_cov  = paste(.cov_labels(dir, FIT_BOAT) %||% character(0), collapse = "+"),
    shore_pc = .comp(dir, "shore (Pot closure)"), shore_ag = .comp(dir, "shore (All gear)"),
    boat_pc = .comp(dir, "private_boat (Pot closure)"), boat_ag = .comp(dir, "private_boat (All gear)"),
    shore_ag_effort = .comp(dir, "shore (All gear)", "BSS_effort"), boat_ag_effort = .comp(dir, "private_boat (All gear)", "BSS_effort"),
    census = .comp(dir, "comm_charter (census)", "PE_catch"),
    port = .num1(pt$BSS_median), port_lo95 = .num1(pt$BSS_lo95), port_hi95 = .num1(pt$BSS_hi95),
    n_fits_bss = if (is.null(gt)) NA_integer_ else sum(gt$method_selected == "BSS", na.rm = TRUE),
    n_fits = if (is.null(gt)) NA_integer_ else nrow(gt),
    gate_all_pass = if (is.null(gt)) NA else all(as.logical(gt$pass_convergence), na.rm = TRUE),
    shore_res = as.character(als$ar_resolution %||% NA), shore_P_n = .num1(als$P_n),
    boat_res = as.character(alb$ar_resolution %||% NA), boat_P_n = .num1(alb$P_n),
    shore_p_loo_frac = .num1(aqs$p_loo_frac), shore_pareto_bad = .num1(aqs$n_pareto_bad), shore_catch_n = nobs(FIT_SHORE),
    boat_p_loo_frac = .num1(aqb$p_loo_frac), boat_pareto_bad = .num1(aqb$n_pareto_bad), boat_catch_n = nobs(FIT_BOAT),
    shore_div = .num1(als$divergences), boat_div = .num1(alb$divergences),
    B_sca_shore = .num1(bs$mean), B_sca_shore_lo = .num1(bs$lo), B_sca_shore_hi = .num1(bs$hi),
    B_sca_boat = .num1(bb$mean), B_sca_boat_lo = .num1(bb$lo), B_sca_boat_hi = .num1(bb$hi),
    B_bar_boat = .num1(bbar$mean), B_bar_boat_lo = .num1(bbar$lo), B_bar_boat_hi = .num1(bbar$hi),
    B_sca_boat_winter = .num1(bw$mean), B_sca_boat_winter_lo = .num1(bw$lo), B_sca_boat_winter_hi = .num1(bw$hi),
    B_sca_boat_rest = .num1(br$mean), B_sca_boat_rest_lo = .num1(br$lo), B_sca_boat_rest_hi = .num1(br$hi),
    elpd_shore_gear = .num1(ls$elpd_loo), elpd_boat_trailer = .num1(lb$elpd_loo),
    elpd_shore_catch = .num1(lsc$elpd_loo), elpd_boat_catch = .num1(lbc$elpd_loo),
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# PER-RUNG VERDICTS
# ---------------------------------------------------------------------------
verdict_M1 <- function(dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  ref <- .here("05_output", REF_R4)
  if (!dir.exists(ref)) {
    V1row("M1", "bit-identity against the pre-covariate R4 render", "reference folder absent", REF_R4, "REVIEW",
          "Without it the inertness proof is M0's structural row only.")
    return(invisible(NULL))
  }
  # keys that legitimately differ between the R4 folder and this one: the block this patch
  # added (absent from R4's dump), the rung's own lever, and the NULL-vs-absent keys
  new_keys <- grep("^(marine_hazard_|bar_restriction_)", names(BASE), value = TRUE)
  # ar_force and bss_sampler_override are NULL in run_config and resolve_cfg()'s `cfg[[k]] <-`
  # DROPS a NULL-valued key, while R4's dump carries them as NULL: absent-vs-NULL, not a change
  ex <- tryCatch(fit_exactness(dir, ref, what = "the four fits",
                   expect_delta = unique(c(new_keys, DELTA_KEYS, RUNTIME_KEYS, "run_weather", "pot_closures",
                                           "census_windows", "ar_force", "bss_sampler_override",
                                           "tau_shore_derive_window_only", "gear_period_bss", "catch_zi_tracks"))),
                 error = function(e) NULL)
  ok <- !is.null(ex) && identical(ex$verdict, "PASS")
  V1row("M1", "THE MARINE CODE IS INERT WHEN OFF: M1 is bit-identical to the pre-covariate R4 render",
        if (is.null(ex)) "could not compare" else ex$observed,
        sprintf("every shared parameter row identical to %s", basename(REF_R4)),
        if (ok) "PASS" else "FAIL",
        paste("M1 ships R4's configuration and differs from it only in carrying the marine module",
              "(which reads nothing when off) and a prep that adds an empty vector to the K_open",
              "extras. If this FAILS, look at the unexpected-delta list and at rstan's version",
              "before reading any other rung: a baseline that is not the baseline measures the",
              "patch as well as the covariate."))
  cmp <- c("shore (Pot closure)", "shore (All gear)", "private_boat (Pot closure)", "private_boat (All gear)")
  dd <- vapply(cmp, function(k) .comp(dir, k), numeric(1)); rr <- vapply(cmp, function(k) .comp(ref, k), numeric(1))
  V1row("M1", "every BSS component reproduces the committed R4 figure exactly",
        paste(sprintf("%s %s vs %s", cmp, fmt(dd, 0), fmt(rr, 0)), collapse = "; "), "identical, to the crab",
        if (isTRUE(all(is.finite(dd)) && identical(dd, rr))) "PASS" else "FAIL",
        "The components come straight from the fits; the port adds the census as a random draw and is reported with a tolerance below.")
  p <- .num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_median)
  pr <- .num1(.port_row(rd(ref, "port_total_Dungeness_Kept.csv"))$BSS_median)
  V1row("M1", "the port total reproduces R4 within the census draw",
        sprintf("M1 %s vs R4 %s (%.4f%%)", fmt(p, 0), fmt(pr, 0), .pct(p, pr)), "within 0.3%",
        if (isTRUE(abs(.pct(p, pr)) < 0.3)) "PASS" else "FAIL",
        paste("The port median is assembled by pairing each component's posterior draws row by row, and",
              "rstan permutes draws with an UNSEEDED sample.int(), so bit-identical fits gave a port that",
              "moved 0.128% here (94,497 vs 94,376) and 'about 0.2%' in batch_verdict_helpers.R's record;",
              "the census draw adds a little more. The first version of this rule asked for 0.05% and",
              "FAILED a correct baseline (C row, 2026-09-26). 0.3% is the documented jitter with room;",
              "the components above are the exact test. B38 rebuilds the permutation from the Stan seed",
              "and seeds the census draw, so two post-B38 renders can be compared exactly; R4 predates it."))
}

# every fitted rung: the gate, the resolution, the adequacy (rules 1 and 2)
verdict_rung <- function(sid, dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  L <- LAD[[sid]]; B <- LAD[["M1"]]; if (is.null(L)) return(invisible(NULL))
  V1row(sid, "rule 1: every fit in the rung passes the convergence gate",
        if (is.na(L$gate_all_pass)) "NOT COMPUTABLE (convergence_report.csv missing)"
        else sprintf("%s of %s fits report BSS; gate all-pass = %s", fmt(L$n_fits_bss, 0), fmt(L$n_fits, 0), L$gate_all_pass),
        "all fits pass", if (is.na(L$gate_all_pass)) "REVIEW" else if (isTRUE(L$gate_all_pass)) "PASS" else "FAIL",
        "A rung with a PE fallback inside it reports a different estimator from the other rungs'. A missing statistic is REVIEW, never a verdict (B33).")
  if (!is.null(B) && !identical(sid, "M1")) {
    known <- !any(is.na(c(L$shore_res, B$shore_res, L$boat_res, B$boat_res, L$shore_P_n, B$shore_P_n, L$boat_P_n, B$boat_P_n)))
    same <- identical(L$shore_res, B$shore_res) && identical(L$boat_res, B$boat_res) &&
            isTRUE(L$shore_P_n == B$shore_P_n) && isTRUE(L$boat_P_n == B$boat_P_n)
    V1row(sid, "rule 1: the all-gear fits used the SAME AR resolution as M1",
          if (!known) "NOT COMPUTABLE (ar_escalation_log.csv missing in one of the two folders)" else
          sprintf("shore %s (P_n %s) vs M1 %s (%s); boat %s (P_n %s) vs M1 %s (%s)",
                  L$shore_res, fmt(L$shore_P_n, 0), B$shore_res, fmt(B$shore_P_n, 0),
                  L$boat_res, fmt(L$boat_P_n, 0), B$boat_res, fmt(B$boat_P_n, 0)),
          "identical resolutions", if (!known) "REVIEW" else if (same) "PASS" else "FAIL",
          "The caps are pinned, so a change here would mean the covariate moved the data-driven selector; the comparison would then be two things at once.")
  }
  for (side in c("shore", "boat")) {
    pf <- L[[paste0(side, "_p_loo_frac")]]; nb <- L[[paste0(side, "_pareto_bad")]]; nn <- L[[paste0(side, "_catch_n")]]
    share <- if (isTRUE(is.finite(nb)) && isTRUE(is.finite(nn)) && nn > 0) nb / nn else NA_real_
    known <- isTRUE(is.finite(pf)) && isTRUE(is.finite(share))
    V1row(sid, sprintf("rule 2: %s all-gear adequacy is not degraded", side),
          if (!known) sprintf("NOT COMPUTABLE (p_loo_frac %s from model_adequacy.csv; bad k %s of %s catch obs from ppc_byobs)", fmt(pf, 4), fmt(nb, 0), fmt(nn, 0)) else
          sprintf("p_loo_frac %s; %s bad k of %s catch obs (%s%%)%s", fmt(pf, 4), fmt(nb, 0), fmt(nn, 0), fmt(100 * share, 1),
                  if (!is.null(B) && !identical(sid, "M1")) sprintf("; M1 p_loo_frac %s", fmt(B[[paste0(side, "_p_loo_frac")]], 4)) else ""),
          "p_loo_frac <= 0.15 and bad k <= 5% of n_obs",
          if (!known) "REVIEW" else if (pf <= 0.15 && share <= 0.05) "PASS" else "FAIL",
          "The D3 run's calibration: the overfitted pooled daily fit sat at 0.352 with 41 bad k; the adopted weekly fit at 0.095 with 0.")
  }
}

# the covariate clauses (rules 3, 4, 5) for one covariate in one fit, `ctl` being the rung it
# is measured against (M1 for SCA; M2 for the bar's increment beyond the archive)
verdict_covariate <- function(sid, dir, ctl, dir_ctl, fit, covariate) {
  if (is.na(dir %||% NA) || is.na(dir_ctl %||% NA)) return(invisible(NULL))
  side <- if (startsWith(fit, "shore")) "shore" else "boat"
  b <- .b_open(dir, fit, covariate)
  V1row(sid, sprintf("rule 3: %s is IDENTIFIED in the %s all-gear fit (95%% interval excludes 0)", covariate, side),
        if (is.null(b)) "term not in this fit (dropped as unidentifiable in the window, or not selected)"
        else sprintf("%s = %s [%s, %s], rate ratio %s; Rhat %s, n_eff %s", b$parameter, fmt(b$mean, 3), fmt(b$lo, 3), fmt(b$hi, 3),
                     fmt(exp(b$mean), 3), fmt(b$rhat, 3), fmt(b$n_eff, 0)),
        "interval excludes 0", if (is.null(b)) "REVIEW" else if (isTRUE(b$identified)) "PASS" else "FAIL",
        paste("A term the data cannot place on one side of zero is a free parameter widening the level,",
              "which is what the razor-dig B3 was (elpd within 1 SE, port +0.6%)."))
  stream <- .stream_of(fit)
  pa <- file.path(dir_ctl, sprintf("loo_pointwise_%s_%s.csv", stream, fit))
  pb <- file.path(dir, sprintf("loo_pointwise_%s_%s.csv", stream, fit))
  el <- tryCatch(loo_elpd_paired(pa, pb, label = sprintf("%s -> %s (%s)", ctl, sid, covariate)), error = function(e) NULL)
  ratio <- if (is.null(el)) NA_real_ else el$ratio %||% NA_real_
  V1row(sid, sprintf("rule 4: %s, %s all-gear %s stream: paired elpd against %s", covariate, side, stream, ctl),
        if (is.null(el)) "NOT COMPUTABLE (loo_pointwise files missing or misaligned)" else loo_elpd_paired_str(el),
        "> +2 paired SE = better on the sampled days; within +-2 SE = no evidence; < -2 SE = worse",
        if (!isTRUE(is.finite(ratio))) "REVIEW" else if (ratio > 2) "PASS" else if (ratio < -2) "FAIL" else "REVIEW",
        paste("Paired, not naive (loo_elpd_paired.R). REVIEW between the bands is 'no evidence either way',",
              "not a failure. And rule 8: this scores the SAMPLED days only."))
  pa <- file.path(dir_ctl, sprintf("loo_pointwise_catch_%s.csv", fit)); pb <- file.path(dir, sprintf("loo_pointwise_catch_%s.csv", fit))
  ec <- tryCatch(loo_elpd_paired(pa, pb, label = sprintf("%s -> %s catch", ctl, sid)), error = function(e) NULL)
  rc <- if (is.null(ec)) NA_real_ else ec$ratio %||% NA_real_
  V1row(sid, sprintf("rule 5: %s, %s all-gear CATCH stream unmoved (paired elpd within +-2 SE of %s)", covariate, side, ctl),
        if (is.null(ec)) "NOT COMPUTABLE" else loo_elpd_paired_str(ec), "within +-2 paired SE",
        if (!isTRUE(is.finite(rc))) "REVIEW" else if (abs(rc) <= 2) "PASS" else "REVIEW",
        "An effort covariate has no business moving the catch fit; a move is something to understand, not a verdict.")
}

verdict_port <- function(sid, dir) {
  L <- LAD[[sid]]; B <- LAD[["M1"]]; if (is.null(L) || is.null(B) || identical(sid, "M1")) return(invisible(NULL))
  V1row(sid, "REPORTED, NOT A CRITERION: the estimate under the covariate(s)",
        sprintf("shore all-gear %s -> %s (%s%%); boat all-gear %s -> %s (%s%%); port %s -> %s (%s%%) [%s, %s]",
                fmt(B$shore_ag, 0), fmt(L$shore_ag, 0), fmt(.pct(L$shore_ag, B$shore_ag), 2),
                fmt(B$boat_ag, 0), fmt(L$boat_ag, 0), fmt(.pct(L$boat_ag, B$boat_ag), 2),
                fmt(B$port, 0), fmt(L$port, 0), fmt(.pct(L$port, B$port), 2), fmt(L$port_lo95, 0), fmt(L$port_hi95, 0)),
        "no threshold; rule 7", "INFO",
        "The direction to expect: a covariate that lowers effort on advisory days lowers the unsampled-day interpolation on those days, so the total moves DOWN if anything.")
}

verdict_M5 <- function(dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  sel <- rd(dir, "marine_hazard_selection.csv")
  V1row("M5", "REPORTED: what the auto screen selected in the rendered run",
        if (is.null(sel)) "marine_hazard_selection.csv missing" else
          paste(sprintf("%s/%s p_adj %s -> %s", sel$population, sel$covariate,
                        fmt(suppressWarnings(as.numeric(sel$p_adj)), 4), sel$selected), collapse = "; "),
        "no threshold", "INFO", "Compare with the forced rungs: auto should land on the terms the rule adopts, or say why not.")
}

# B41: the two season coefficients side by side, and whether they differ. The z below treats
# the two posteriors as independent (they act on disjoint days; the shared monthly state
# couples them a little), so it is a guide to read beside the intervals, not a criterion.
verdict_M6 <- function(dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  bw <- .b_open(dir, FIT_BOAT, "nws_sca_any_winter"); br <- .b_open(dir, FIT_BOAT, "nws_sca_any_rest")
  fs <- .full(dir, FIT_BOAT)
  sdw <- if (!is.null(bw) && !is.null(fs) && "sd" %in% names(fs)) .num1(fs[bw$parameter, "sd"]) else NA_real_
  sdr <- if (!is.null(br) && !is.null(fs) && "sd" %in% names(fs)) .num1(fs[br$parameter, "sd"]) else NA_real_
  z <- if (!is.null(bw) && !is.null(br) && isTRUE(is.finite(sdw)) && isTRUE(is.finite(sdr))) (bw$mean - br$mean) / sqrt(sdw^2 + sdr^2) else NA_real_
  V1row("M6", "REPORTED: the winter and rest SCA coefficients in the boat all-gear fit, and their difference",
        if (is.null(bw) || is.null(br)) sprintf("winter %s; rest %s", if (is.null(bw)) "absent (not in the fit)" else sprintf("%s [%s, %s]", fmt(bw$mean, 3), fmt(bw$lo, 3), fmt(bw$hi, 3)),
                                                if (is.null(br)) "absent (not in the fit)" else sprintf("%s [%s, %s]", fmt(br$mean, 3), fmt(br$lo, 3), fmt(br$hi, 3)))
        else sprintf("winter %s [%s, %s] (rate ratio %s); rest %s [%s, %s] (rate ratio %s); difference winter - rest %s, z %s (independence-approximate)",
                     fmt(bw$mean, 3), fmt(bw$lo, 3), fmt(bw$hi, 3), fmt(exp(bw$mean), 2), fmt(br$mean, 3), fmt(br$lo, 3), fmt(br$hi, 3), fmt(exp(br$mean), 2),
                     fmt(bw$mean - br$mean, 3), fmt(z, 2)),
        "no threshold; rule 9 reads the intervals and the block CV", "INFO",
        "Rule 9, in order: (a) winter not identified -> D32 open, no verdict on the constant term; (b) both identified and different (the difference's 95% interval excludes 0) -> the split earns its place if the catch is unmoved and the block CV's joint boat all-gear row is not below -2 SE; (c) both identified, no detectable difference -> the constant term stands.")
  branch <- if (is.null(bw) || !isTRUE(bw$identified)) "a" else if (is.null(br) || !isTRUE(br$identified)) "a-rest"
            else if (isTRUE(is.finite(z)) && abs(z) > qnorm(0.975)) "b" else "c"
  V1row("M6", "rule 9: which branch the season split lands in (a: winter not identified; b: both identified and different; c: both identified, no detectable difference)",
        switch(branch,
               a = "(a) the winter coefficient is not identified: this season's counts cannot place the winter effect; D32 stays open; no verdict on the constant term",
               `a-rest` = "the winter coefficient is identified and the REST is not: read the fit before anything else (this was not the expected shape)",
               b = sprintf("(b) both identified and different (z %s): the split earns its place if rule 5 holds and the block CV's joint boat all-gear row for M6 vs M2 is not below -2 SE", fmt(z, 2)),
               c = sprintf("(c) both identified, no detectable difference (z %s): the constant term stands", fmt(z, 2))),
        "stated in the header before the render", "INFO",
        "The block CV's row is in marine_hazard_2026-09-26_blockcv_pairs.csv after run_marine_block_cv_2026-09-26.R is re-sourced with M6 rendered.")
}

# ---------------------------------------------------------------------------
# THE RECOMMENDATION, from the ladder and the verdicts, by the rule in the header.
# ---------------------------------------------------------------------------
REC <- list()
recommend <- function() {
  banner("RECOMMENDATION: the SCA term and the bar-restriction term")
  if (!length(LAD)) { cat("  No fitted rungs to read. Nothing recommended.\n"); return(invisible(NULL)) }
  VV <- do.call(rbind, V)
  # the rows that bear on ONE covariate in ONE fit: rule 1 (the rung), rule 2 (that side's
  # adequacy), rules 3 to 5 (that covariate in that side's all-gear fit)
  # B33's lesson, applied here: a clause that is ABSENT (a rung or its control did not render)
  # or NOT COMPUTABLE (a file missing, the term not in the fit) is an open question, never a
  # verdict in either direction. Only when every one of rules 1 to 5 is present and evaluated
  # does the ladder get to say adopt, and rule 4's REVIEW inside +-2 SE (an evaluated "no
  # evidence either way") is named as such rather than lumped with "could not be evaluated".
  judge <- function(sid, covariate, side) {
    rows <- VV[VV$stage == sid & grepl("^rule ", VV$criterion), , drop = FALSE]
    keep <- grepl("^rule 1", rows$criterion) |
            (grepl("^rule 2", rows$criterion) & grepl(paste0(": ", side, " all-gear"), rows$criterion, fixed = TRUE)) |
            (grepl("^rule [345]", rows$criterion) & grepl(covariate, rows$criterion, fixed = TRUE) &
               grepl(paste0(side, " all-gear"), rows$criterion, fixed = TRUE))
    rows <- rows[keep, , drop = FALSE]
    if (!nrow(rows)) return(list(verdict = "not fitted", rows = rows))
    present <- vapply(1:5, function(n) any(grepl(paste0("^rule ", n), rows$criterion)), logical(1))
    if (!all(present))
      return(list(verdict = sprintf("open: rule(s) %s not evaluated (the rung or its control did not render)",
                                    paste(which(!present), collapse = ", ")), rows = rows))
    if (any(rows$verdict == "FAIL"))  return(list(verdict = "do not adopt", rows = rows))
    if (any(rows$verdict == "ERROR")) return(list(verdict = "open: a clause errored; read the verdicts table", rows = rows))
    not_computable <- grepl("^NOT COMPUTABLE|^term not in this fit", rows$observed)
    if (any(not_computable | (rows$verdict == "REVIEW" & grepl("^rule [123]", rows$criterion))))
      return(list(verdict = "open: a clause could not be evaluated (a file missing, or the term absent from the fit)", rows = rows))
    r4 <- rows[grepl("^rule 4", rows$criterion), , drop = FALSE]
    r5 <- rows[grepl("^rule 5", rows$criterion), , drop = FALSE]
    if (any(r5$verdict != "PASS"))
      return(list(verdict = "open: the catch stream moved (rule 5); understand why before adopting", rows = rows))
    list(verdict = if (all(r4$verdict == "PASS")) "adopt (sampled-day evidence)"
                   else "identified and harmless, no sampled-day gain (rule 4 within +-2 SE): adoption would rest on the mechanism, not on this run",
         rows = rows)
  }
  cat("\n  The ladder:\n")
  ladder <- do.call(rbind, LAD)
  print(ladder[, c("rung", "mode", "shore_cov", "boat_cov", "shore_ag", "boat_ag", "port", "gate_all_pass",
              "B_sca_shore", "B_sca_boat", "B_bar_boat", "B_sca_boat_winter", "B_sca_boat_rest", "elpd_shore_gear", "elpd_boat_trailer")], row.names = FALSE)
  items <- list(list(name = "SCA on the shore", sid = "M2", cov = "nws_sca_any", side = "shore"),
                list(name = "SCA on the boat",  sid = "M2", cov = "nws_sca_any", side = "boat"),
                list(name = "bar restriction on the boat (beyond the archive, M4 vs M2)", sid = "M4", cov = "bar_restriction", side = "boat"))
  if ("M6" %in% STAGES) items <- c(items, list(
    list(name = "SCA split by season on the boat, winter term (M6 vs M2; rule 9)", sid = "M6", cov = "nws_sca_any_winter", side = "boat"),
    list(name = "SCA split by season on the boat, rest-of-season term (M6 vs M2; rule 9)", sid = "M6", cov = "nws_sca_any_rest", side = "boat")))
  # B41: the M6 items are read by rule 9, not by judge()'s rules 1 to 5: an unidentified winter
  # coefficient is D32 restated, not "do not adopt", and a rule-4 PASS is not "adopt".
  judge_M6 <- function(covariate) {
    rows <- VV[VV$stage == "M6" & (grepl("^rule [1235]", VV$criterion) | grepl("^rule 9", VV$criterion)), , drop = FALSE]
    rows <- rows[!grepl("^rule 3", rows$criterion) | grepl(covariate, rows$criterion, fixed = TRUE), , drop = FALSE]
    rows <- rows[!grepl("^rule 5", rows$criterion) | grepl(covariate, rows$criterion, fixed = TRUE), , drop = FALSE]
    if (!nrow(rows)) return(list(verdict = "not fitted", rows = rows))
    br <- rows[grepl("^rule 9", rows$criterion), , drop = FALSE]
    if (!nrow(br)) return(list(verdict = "open: rule 9 not evaluated", rows = rows))
    o <- br$observed[1]
    v <- if (startsWith(o, "(a)")) "rule 9 (a): winter not identified; D32 open; no verdict on the constant term"
         else if (startsWith(o, "(b)")) paste("rule 9 (b): both identified and different; the split earns its place if rule 5 holds and the block CV's joint boat all-gear row",
                                              "(run_marine_block_cv_2026-09-26.R, M6 vs M2) is not below -2 SE: read that file")
         else if (startsWith(o, "(c)")) "rule 9 (c): both identified, no detectable difference; the constant term stands"
         else "open: an unexpected shape; read the fit"
    if (any(rows$verdict == "FAIL" & grepl("^rule [125]", rows$criterion))) v <- paste(v, "| BUT a rule 1, 2 or 5 row FAILED: read it first")
    list(verdict = v, rows = rows)
  }
  for (item in items) {
    j <- if (identical(item$sid, "M6")) judge_M6(item$cov) else judge(item$sid, item$cov, item$side)
    cat(sprintf("\n  %s: %s\n", item$name, toupper(j$verdict)))
    for (i in seq_len(nrow(j$rows))) cat(sprintf("     %-6s %s\n", j$rows$verdict[i], j$rows$criterion[i]))
    REC[[item$name]] <<- j$verdict
  }
  cat("\n  RULE 8, restated because it bounds everything above. Every elpd here is on the SAMPLED\n")
  cat("  days. The term's value to the ESTIMATE is on the unsampled days, and no clause above\n")
  cat("  measures that. 'adopt (sampled-day evidence)' means: identified, harmless to the catch\n")
  cat("  fit, better on the sampled days. Before the mode ships 'auto' or 'on' in run_config.R,\n")
  cat("  the leave-one-week-out block CV (CHANGE_REGISTER D31 -> B39) is the test that would show\n")
  cat("  the interpolation improved: 06_diagnostics/run_marine_block_cv_2026-09-26.R scores the\n")
  cat("  rendered rungs' saved draws (ppc_draws_<fit>.rds) without a refit. It ran 2026-09-26\n")
  cat("  (Section 1y): boat SCA +18.9 nats on the held-out OSP counts (3.9 SE), both boat streams\n")
  cat("  2.9 SE, the winter untested (D32). It ran again 2026-09-27 with M6: M6 vs M2 +1.2 at 0.6 SE.\n")
  adopted <- identical(BASE$marine_hazard_mode, "manual") && identical(BASE$marine_hazard_manual_boat, "nws_sca_any") &&
             !length(BASE$marine_hazard_manual_shore %||% character(0)) && identical(BASE$marine_hazard_gear_regimes, "all_gear")
  if (adopted) {
    cat("\n  ADOPTED 2026-09-27 (Matt's decision, Section 1z; A30): run_config.R section 2.10 ships manual,\n")
    cat("  boat nws_sca_any, no shore term, marine_hazard_gear_regimes = \"all_gear\". These rungs are the record\n")
    cat("  the decision was taken on. The production render of the shipped file is the authoritative run since 2026-09-28\n")
    cat("  (05_output/20260927/pooled-CPUE-canonical-2024-25, 96,118); the gear-resolved cross-check under it is the owed step.\n")
  } else {
    cat("\n  ADOPTION EDIT, if you take it: run_config.R section 2.10, marine_hazard_mode <- \"manual\" with\n")
    cat("  marine_hazard_manual_boat = \"nws_sca_any\", no shore term, marine_hazard_gear_regimes = \"all_gear\"\n")
    cat("  (what Section 1z adopted); \"auto\" (the screen) was M5 and adds the redundant bar column.\n")
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
banner(sprintf("MARINE HAZARD EFFORT COVARIATES  |  DRY_RUN = %s", DRY_RUN))
cat(sprintf("  stages: %s\n", paste(STAGES, collapse = ", ")))
cat(sprintf("  code   %s\n", .code_fingerprint()))
if (isTRUE(DRY_RUN))
  cat("\n  DRY RUN. M0 runs; no rung is fitted. Set DRY_RUN <- FALSE and source again.\n")

DIRS <- list()
for (sid in STAGES) {
  if (identical(sid, "M0")) { tryCatch(stage_M0(), error = function(e)
    V1row("M0", "the desk stage completed", conditionMessage(e), "no error", "ERROR",
          "M0 is pure desk work; an error here is a bug in this file, not a run failure.")); next }
  DIRS[[sid]] <- tryCatch(run_one(sid), error = function(e) {
    V1row(sid, "the rung completed", conditionMessage(e), "the render runs to completion", "ERROR",
          "The render errored. The other rungs are separate renders; RESUME = TRUE retries this one.")
    NA_character_ })
  if (!is.na(DIRS[[sid]] %||% NA)) tryCatch(ladder_row(sid, DIRS[[sid]]), error = function(e)
    V1row(sid, "the ladder row was extracted", conditionMessage(e), "no error", "ERROR", "The fit survived; only the extraction failed."))
}
fitted <- intersect(STAGES, names(STAGE_DEFS)[vapply(STAGE_DEFS, function(s) isTRUE(s$fit), logical(1))])
for (sid in fitted) tryCatch(verdict_rung(sid, DIRS[[sid]] %||% NA_character_), error = function(e)
  V1row(sid, "the rung verdicts were computed", conditionMessage(e), "no error", "ERROR", ""))
if ("M1" %in% STAGES) tryCatch(verdict_M1(DIRS$M1 %||% NA_character_), error = function(e)
  V1row("M1", "the bit-identity verdict was computed", conditionMessage(e), "no error", "ERROR", ""))
# the covariate clauses: SCA in M2 against M1 (both fits), bar in M3 against M1 and in M4 against M2
.cov_pairs <- list(
  list(sid = "M2", ctl = "M1", fit = FIT_SHORE, cov = "nws_sca_any"),
  list(sid = "M2", ctl = "M1", fit = FIT_BOAT,  cov = "nws_sca_any"),
  list(sid = "M3", ctl = "M1", fit = FIT_BOAT,  cov = "bar_restriction"),
  list(sid = "M4", ctl = "M2", fit = FIT_BOAT,  cov = "bar_restriction"),
  list(sid = "M4", ctl = "M1", fit = FIT_SHORE, cov = "nws_sca_any"),
  list(sid = "M4", ctl = "M1", fit = FIT_BOAT,  cov = "nws_sca_any"),
  # B41: the season split against the constant term
  list(sid = "M6", ctl = "M2", fit = FIT_BOAT,  cov = "nws_sca_any_winter"),
  list(sid = "M6", ctl = "M2", fit = FIT_BOAT,  cov = "nws_sca_any_rest"))
for (cp in .cov_pairs) if (cp$sid %in% STAGES && cp$ctl %in% STAGES)
  tryCatch(verdict_covariate(cp$sid, DIRS[[cp$sid]] %||% NA_character_, cp$ctl, DIRS[[cp$ctl]] %||% NA_character_, cp$fit, cp$cov),
           error = function(e) V1row(cp$sid, sprintf("the %s verdicts were computed", cp$cov), conditionMessage(e), "no error", "ERROR", ""))
for (sid in setdiff(fitted, "M1")) tryCatch(verdict_port(sid, DIRS[[sid]] %||% NA_character_), error = function(e) NULL)
if ("M5" %in% STAGES) tryCatch(verdict_M5(DIRS$M5 %||% NA_character_), error = function(e) NULL)
if ("M6" %in% STAGES) tryCatch(verdict_M6(DIRS$M6 %||% NA_character_), error = function(e) NULL)

if (length(V)) {
  banner("VERDICTS")
  VV <- do.call(rbind, V)
  for (i in seq_len(nrow(VV)))
    cat(sprintf("  %-6s %-4s %s\n         observed : %s\n         threshold: %s\n",
                VV$verdict[i], VV$stage[i], VV$criterion[i], VV$observed[i], VV$threshold[i]))
  cat(sprintf("\n  %d PASS  %d FAIL  %d REVIEW  %d INFO  %d ERROR\n",
              sum(VV$verdict == "PASS"), sum(VV$verdict == "FAIL"),
              sum(VV$verdict == "REVIEW"), sum(VV$verdict == "INFO"), sum(VV$verdict == "ERROR")))
}
tryCatch(recommend(), error = function(e) cat("  recommendation failed:", conditionMessage(e), "\n"))

if (!isTRUE(DRY_RUN)) {
  op <- .here("05_output")
  if (length(LAD)) merge_csv_by(do.call(rbind, LAD), file.path(op, "marine_hazard_2026-09-25_ladder.csv"), key = "rung")
  if (length(V))   merge_csv_by(do.call(rbind, V), file.path(op, "marine_hazard_2026-09-25_verdicts.csv"), key = c("stage", "criterion"))
  if (length(REC)) merge_csv_by(data.frame(item = names(REC), recommendation = unlist(REC),
                                           written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), code = .code_fingerprint(),
                                           stringsAsFactors = FALSE),
                                file.path(op, "marine_hazard_2026-09-25_recommendation.csv"), key = "item")
  cat("\n  wrote marine_hazard_2026-09-25_{ladder,verdicts,recommendation}.csv to 05_output/\n")
} else {
  cat("\n  DRY_RUN: nothing written. Set DRY_RUN <- FALSE and source again to start.\n")
}
