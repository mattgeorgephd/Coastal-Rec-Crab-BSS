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
# IMPROVEMENT LADDER 2026-09-08: the eight review items, one rung at a time
# -----------------------------------------------------------------------------
# WHY THIS RUN EXISTS
#   The 2026-09-08 patch series (items 1 to 8 of the 2026-09-06 branch review) changes
#   the boat turnover prior (item 3), the boat PE (5), the fishing-time filters (8), the
#   crabbing fraction f (1A, 1B), the gear split (6), the census reporting (4) and, OFF
#   by default, the shore turnover (2). Three of those push the port total the same way
#   (up), so they must be seen one at a time or their effects cannot be attributed. Each
#   rung below adds ONE change on top of the previous rung; every rung is a full pooled
#   render (four fits), so a rung is comparable to the last one and to the baseline.
#
# THE RUNGS (cumulative)
#   R0   desk, seconds: the PE-only and reporting changes (items 4, 5, 6) plus the data
#        the fitted rungs will read (contacts, calibration, shore turnover). No fit.
#   R1   pre-patch configuration + item 8 (the filters are code, always on)
#   R2   + item 3: tau_boat prior from the OSP/trailer calibration (shared_tau_sigma 0.15)
#   R3a  + item 1A: monthly f from the sampler contacts, LEGACY per-stratum construction
#   R3   + item 1B: the dynamic f (the shipped configuration)
#   R4   + item 2: tau_shore_prior_mu = "derived" (ADOPTED 2026-09-09: the shipped configuration)
#   R5   the gear-resolved cross-check on the configuration GEAR_FOLLOWS names
#
# 2026-09-09 (after the patches were applied): the census is exact over the tally days
# (census_expansion = "none", every rung), the contacts carry trip types and the combo
# share c has its own walk (R3a onward read the rebuilt interview workbook), and the
# derived shore turnover is adopted, so R4 is now the shipped configuration and R5
# follows it by default. R1 and R2 therefore compare to the baseline with the census
# 3,869 crab lower by design (11,821 -> 7,884 on 2024-25) and the tampered-gear filter
# live (13 Grays Harbor interviews of 2024-25); read the R1 verdict with that in mind.
#
# 2026-09-11: the census component is now a commercial CENSUS plus a charter EXPANSION
# (R1 and R2 keep the pre-patch tally frame, R4 the shipped roster frame: 7,884 -> 8,538
# on 2024-25), and the crabbing-holiday calendar gained Thanksgiving Day, Veterans Day
# (observed) and Juneteenth, which re-types three days of every season. Neither is a
# fitted-component change, but both move the PE, so R0 records them and the R1 verdict
# compares the FITTED components only.
#
# WHAT EACH RUNG MUST SHOW (verdict rows, written to 05_output/improvements_2026-09-08_verdicts.csv)
#   R1   shore fits within a few percent of the baseline (item 8 adds ~140 shore rows),
#        boat CPUE a few percent lower (22 zero-catch trips restored); nothing else moved.
#   R2   tau_bar near 2.97, boat all-gear +12 to 15%, shore fits BIT-IDENTICAL to R1.
#   R3a  f per month tracks the contact shares where n >= 20; thin months at 0.30.
#   R3   every non-f boat parameter within Monte Carlo error of R3a (the factorization
#        proof; bit-identity is impossible because the parameter vector changed), thin
#        months borrow from their neighbours, sigma_f away from zero, boat total moves by
#        roughly the catch-weighted change in f; shore BIT-IDENTICAL to R3a.
#   R4   shore effort x about 1.47 with the boat BIT-IDENTICAL to R3.
#   R5   port within 2% of the pooled rung it follows; tau_bar agrees.
#
# RUNTIME. About 4 h per pooled rung on 4 cores (the 2026-09-04 A1 render took 227 min),
# 0.5 h for the gear rung: 20 to 25 h for the full ladder. R3a can be dropped from STAGES
# if time is short; it exists so the dynamic model's smoothing is attributable. Every
# rung's Stan model is the current one (the dynamic f block recompiles once).
#
# WINDOW. Every rung fits the single 2024-25 season (the run_config.R rollback block) so
# the rungs compare to the 2026-09-04 A1 baseline; the live two-season configuration is
# not what this ladder measures (since 2026-09-10 it is no longer blocked on data: the
# 2023-24 effort counts, gear counts and holidays exist; only its 2023-24 census frame
# does not).
#
# 2026-09-10 (the inputs rebuilt from the season workbooks): every rung reads the rebuilt
# workbooks. For the 2024-25 season that changes nothing in the fitted components (the
# 2024-25 rows reproduce the earlier workbooks) and one thing in the census: the charter
# trip roster is now the charter frame (charter_frame = "roster"), which counts the 8
# sailed charter trips that fell on days without a tally, 7,884 -> 8,592 on 2024-25. The
# R0 row states both numbers; charter_frame = "tally" reproduces 7,884.
#
# =============================================================================
# 2026-09-12: SET UP FOR THE FULL RUN. Four changes, all of them about making the
# rungs comparable to each other rather than about adding a rung.
#
# 1. F_METHOD (new control). Matt 2026-09-12: "The full ladder rung run will only be for
#    the 2024-25 season and will use the new method for estimating f." The ladder as
#    designed does the opposite: R1 and R2 fit the RETIRED f (one scalar, Beta(6,14) at
#    0.30) precisely so that R3a and R3 can attribute the f change. Both readings are
#    defensible and they cost different amounts, so this is a switch and not a decision
#    taken here:
#      "new_throughout" (SHIPPED, per the instruction) every fitted rung carries the
#            dynamic monthly f. STAGES becomes R0, R1, R2, R4, R5: three pooled fits and
#            one gear fit, about 12 to 14 h. Every rung's port total is then a number that
#            could be cited, and the rung-to-rung deltas isolate the TURNOVER levers on a
#            common f. WHAT IT COSTS: the f change itself is never attributed by a fit,
#            and the factorization proof (that f touches generated quantities only, so
#            effort and CPUE posteriors are unchanged in distribution) is never run. Add
#            "R2f" to STAGES to recover the proof for one extra fit: R2f is R2 with the
#            f block rolled back, so R2-vs-R2f is the same test R3-vs-R3a was, at half
#            the cost and against the CURRENT turnover prior.
#      "ladder" the original design: R0, R1, R2, R3a, R3, R4, R5. Five pooled fits, about
#            20 to 25 h, and the f change is attributed twice over (legacy monthly, then
#            the walk). Two of the five rungs are then fitted on a retired f.
#    The mode is in every run_tag and every output filename, so the two cannot be mixed
#    in one verdicts file.
#
# 2. THE WINDOW PIN WAS INCOMPLETE. WINDOW pinned est_date_start, est_date_end,
#    season_filter, pot_closures and census_windows. It did NOT pin pot_open_date,
#    pot_closure_start/end, census_start_date/end or commercial_opener, all of which are
#    per-season and all of which the shipped two-season config sets for 2023-24. Every
#    rung would have run 2024-25 with pot_open_date = 2023-12-01. That one is cosmetic
#    (see run_config.R: plot markers and a NULL fallback, never a model term, and the
#    "L_effective regression split" it was documented as feeding never existed), but
#    pot_closure_start/end drive build_subseasons() and were right only by luck. All of
#    them are pinned below and validate_season_window() now warns if any is left behind.
#
# 3. THE PE's UNSAMPLED-CELL LEVERS ARE PINNED, and this is not cosmetic. A component
#    whose convergence gate FAILS reports its PE point in the port total as a constant,
#    so pe_empty_effort_stratum / pe_empty_stratum / pe_variance can reach the headline
#    of any rung where the thin boat pot closure falls back. run_config's defaults changed
#    on 2026-09-12 (zero -> local_day_type, pooled -> local, new impute_aware variance),
#    which moves the 2024-25 PE port total 72,224 -> 85,076. Pinning them in WINDOW means
#    a later edit to run_config cannot silently make one rung incomparable to the rest,
#    and the preflight now PROVES every rung resolved to the same values.
#
# 4. RESUME CANNOT REUSE A FOLDER BUILT UNDER A DIFFERENT CONFIG. It skipped any rung
#    whose output folder existed and contained run_parameters.txt, with no check that the
#    folder was built from the same delta. Since this file's deltas changed on 2026-09-11
#    and again today, that was a live hazard: a partially-run ladder would have been
#    silently completed with a mix of old and new rungs. Each finished rung now writes
#    IMP_STAGE.txt with a digest of its resolved configuration, and RESUME reuses a folder
#    only when the digest matches. A folder with no digest, or a stale one, is re-run.
# =============================================================================
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                    # TRUE: pre-flight + the R0 desk rung, nothing fitted. START HERE.
F_METHOD <- "new_throughout"       # "new_throughout" (shipped) | "ladder"  -- see note 1 above

# THE TWO-PASS PLAN (Matt 2026-09-13: "run the 4-rung version now and then follow up with
# the new R2f control rung"). Change ONE number between passes:
#
#   LADDER_PASS <- 1   R0, R1, R2, R4, R5.  Three pooled fits + one gear fit, ~12-14 h.
#                      These are the four rungs whose port totals are citable, so the
#                      answer you can act on arrives first.
#   LADDER_PASS <- 2   the same set PLUS R2f, the factorization control. Re-run with
#                      RESUME = TRUE and only R2f fits (~4 h): the pass-1 rungs are matched
#                      by their IMP_STAGE.txt config digest and skipped, and every verdict
#                      is recomputed from the folders already on disk.
#
# Verified before shipping: adding R2f does NOT change any other rung's digest, so pass 2
# reuses exactly the four pass-1 fits. WHAT TO WATCH: the digest covers the CONFIGURATION,
# not the code. If a Stan model, a driver or an R function changes between the two passes,
# R2f is fitted by different code than the R2 it is measured against; that is now detected
# (IMP_STAGE.txt records a code fingerprint), reported at the rung, and any cross-rung
# verdict resting on it is downgraded from PASS to REVIEW. So do not apply another patch
# between the two passes unless you mean to, and if you do, delete the affected rung
# folders so they re-fit.
# ALREADY RUN, 2026-09-11: all five fitted rungs completed in ONE pass with LADDER_PASS = 2
# and R4 fitted first. Results are committed under 05_output/20260910-11/*-newf/ and
# summarised in PIPELINE_STATUS section 1v (port total 94,376 [77,566, 118,602]; 14 PASS,
# 2 REVIEW, both diagnosed as D24 and D25). The two switches below are reset to their SAFE
# values so a clone does not start a 16-hour run by accident; RESUME will reuse every rung
# whose config digest still matches, so re-running costs minutes unless a delta changed.
LADDER_PASS <- 1                   # 1 = the four citable rungs; 2 = adds R2f. See above.

STAGES  <- if (identical(F_METHOD, "ladder")) {
             c("R0", "R1", "R2", "R3a", "R3", "R4", "R5")
           } else if (LADDER_PASS >= 2) {
             c("R0", "R4", "R1", "R2", "R2f", "R5")
           } else {
             c("R0", "R4", "R1", "R2", "R5")
           }
RESUME  <- TRUE                    # reuse a rung ONLY when its IMP_STAGE.txt digest matches
GEAR_FOLLOWS <- "R4"               # the pooled rung whose configuration R5 fits: "R4" (shipped) | "R3"

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
  load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
  install.lib <- load.lib[!load.lib %in% installed.packages()]
  for (lib in install.lib) install.packages(lib, dependencies = TRUE)
  invisible(sapply(load.lib, require, character.only = TRUE))
  rstan_options(auto_write = TRUE)
} else {
  suppressWarnings(suppressPackageStartupMessages(
    try({ library(dplyr); library(tidyr); library(readr); library(lubridate)
          library(readxl); library(here); library(purrr); library(stringr); library(tibble) },
        silent = TRUE)))
}
invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE),
                 function(f) try(source(f), silent = TRUE)))
source(.here("run_config.R"))
BASE <- run_config

model_rmd <- list(pooled        = .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
                  gear_resolved = .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd"))
prefix    <- list(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")
stopifnot(all(file.exists(unlist(model_rmd))))

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
.full_row <- function(dir, fit_pat, par) {
  f <- list.files(dir %||% "", pattern = paste0("^bss_full_summary_", fit_pat, ".*\\.csv$"), full.names = TRUE)
  if (!length(f)) return(NULL)
  x <- tryCatch(utils::read.csv(f[1], row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (is.null(x) || !par %in% rownames(x)) return(NULL)
  as.list(x[par, , drop = FALSE])
}
# THE f/c EXCLUSION LIST for fit_agreement(), in ONE place (2026-09-13).
#
# WHY IT IS A NAMED CONSTANT, AND WHAT IT ACTUALLY BUYS. The list was typed by hand in two
# places, and when the combo-share walk arrived on 2026-09-09 its parameters (z_c, sigma_c,
# cfc_kappa, and the generated sigma_c_out / cfc_kappa_out) were never added to either
# copy. crab_bss_pooled.stan declares sigma_f_out, cfi_kappa_out, sigma_c_out and
# cfc_kappa_out UNCONDITIONALLY and sets each to exactly 0.0 when its walk is off (lines
# 931-935, 966-970), so on an R2-vs-R2f or R3-vs-R3a comparison one side holds a real
# posterior and the other a hard zero.
#
# MEASURED, AND IT CORRECTS A STRONGER CLAIM MADE EARLIER THE SAME DAY. A 2-chain
# 300-iteration fit of the boat all-gear component under the R2f configuration reports
# those quantities as sd = 0 and therefore se_mean = NaN (n_eff and Rhat NaN too), NOT
# se_mean = 0. fit_agreement() skips any row whose combined se is not finite, so those
# rows are skipped whether or not this list names them. The rows that WOULD have been
# compared and would have produced an enormous z -- f_crab_out[*] (a real posterior on
# BOTH sides: 0.307 under R2f against the winter monthly values under R2), plus E_sum and
# C_expected_sum -- were already in the list before today.
#
# So closing the gap is DEFENSIVE HARDENING, not the bug first claimed: the masking rests
# on rstan::summary() reporting NaN rather than 0 for a zero-variance parameter, which is
# an implementation detail and not a property of the design. Worth more than the three
# added names is that harness section 60 now asserts the regex covers every f/c quantity
# BOTH Stan models report, does NOT cover the parameters the proof must actually compare,
# and pins the NaN-skip behaviour the analysis above rests on. The list maintains itself
# the next time an f output is added, and a future rstan that reports 0 instead of NaN
# cannot turn this into the spurious FAIL it currently is not.
#
# Covered: every f and c quantity (parameters, transformed parameters and generated
# quantities), plus the reported TOTALS the f multiplies (E, C, C_expected, lambda_Ctot),
# plus log_lik and lp__.
F_EXCLUDE <- paste0("^(f_crab|f_lower|f_theta|f_lower_param|z_f|sigma_f|cfi_kappa|",
                    "combo_c|z_c|sigma_c|cfc_kappa|eta_f|osp_f_kappa|",
                    "E\\[|E_sum|C\\[|C_sum|C_expected|lambda_Ctot|log_lik|lp__)")

V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = observed,
                                    threshold = threshold, verdict = verdict, why = why,
                                    stringsAsFactors = FALSE)
# 2026-09-13: a cross-rung claim (bit-identity, or agreement within Monte Carlo error) is
# only meaningful if the two folders were produced by the same code. In the two-pass plan
# they need not have been. This wraps such a row: it keeps the observation, appends what
# differs, and refuses to report PASS on a comparison whose premise is broken.
V1cross <- function(stage, criterion, observed, threshold, verdict, why, dir_a, dir_b) {
  note <- .comparability_note(dir_a, dir_b)
  if (nzchar(note) && identical(verdict, "PASS")) verdict <- "REVIEW"
  V1row(stage, criterion, paste0(observed, note), threshold, verdict,
        paste0(why, if (nzchar(note)) paste0("  COMPARABILITY: this claim assumes both rungs were fitted by the same",
                                             " code and the same rstan;", note, ". A PASS is downgraded to REVIEW",
                                             " because the premise, not the result, is what failed.") else ""))
}

# ---------------------------------------------------------------------------
# BASELINE, named beside each value: the 2026-09-04 authoritative render.
# ---------------------------------------------------------------------------
REF <- list(
  A1 = list(dir = "20260904/pooled-CPUE-AD-A1-adopted",
            port = 72027, lo95 = 53018, hi95 = 101364,
            shore_pc = 6320, shore_ag = 21489, boat_pc = 1018, boat_ag = 31008, census = 11821,
            tau_bar = 2.613, f_flat = 0.30),
  # the 2026-09-05 gear cross-check of that baseline
  G2 = list(dir = "20260905/gear-type-CPUE-model-AD-A2-gear-crosscheck"))

# ---------------------------------------------------------------------------
# CONFIGURATIONS. WINDOW pins the single season on every rung; each delta is the
# previous rung's plus one change. Keys set to NULL (shared_tau_sigma at R1) are kept
# as NULL by resolve_cfg(), which is what the pre-patch config had.
# ---------------------------------------------------------------------------
# WINDOW: the COMPLETE single-season 2024-25 pin plus the levers that must not differ
# between rungs. Everything here is applied to every rung, so a run_config edit between
# rungs cannot make one of them incomparable; the preflight proves it rung by rung.
WINDOW <- list(
  # -- the season, all nine per-season keys (2026-09-12: four of them were missing) ----
  est_date_start = "2024-09-16", est_date_end = "2025-09-15", season_filter = "2024-25",
  pot_closures = NULL, census_windows = NULL,
  pot_closure_start = "2024-09-16", pot_closure_end = "2024-11-30",
  pot_open_date = "2024-12-01",              # was inherited as 2023-12-01 from the 2-season config
  census_start_date = "2024-12-01", census_end_date = "2025-02-08",
  commercial_opener = "2025-02-11",          # the 2024-25 coastal commercial opener
  # -- the PE's unsampled-cell levers, PINNED (note 3) ------------------------------
  # These reach the port total of any rung whose gate fails, so they are held at the
  # 2026-09-12 shipped values rather than read from whatever run_config says at run time.
  pe_empty_effort_stratum = "local_day_type", pe_empty_stratum = "local",
  pe_variance = "impute_aware",
  # -- run-level, held fixed so no rung differs in a way nobody declared -------------
  run_weather = FALSE, bss_seed = 20260619, bss_chains = 4, bss_cores = 4,
  bss_sampler_override = NULL, ar_force = NULL, ar_escalate = FALSE,
  estimate_red_rock = FALSE)

# THE F BLOCK, named once so both modes read the same definition.
F_NEW <- list(crab_fraction_strata = "month", crab_fraction_source = "both", crab_fraction_dynamic = TRUE)
F_OLD <- list(crab_fraction_strata = "none",  crab_fraction_source = "ie",   crab_fraction_dynamic = FALSE)

D_R1  <- c(list(tau_boat_prior_mu = 1.2, tau_boat_prior_sigma = 0.3, shared_tau_sigma = NULL,
                tau_shore_prior_mu = 1.7, tau_shore_prior_sigma = 0.3, census_uncertainty = "none",
                # the pre-patch census: the whole component on the tally frame, nothing carried
                charter_frame = "tally", charter_expansion = "pooled"),
           if (identical(F_METHOD, "ladder")) F_OLD else F_NEW)
D_R2  <- modifyList(D_R1,  list(tau_boat_prior_mu = "calibration", tau_boat_prior_sigma = 0.5, shared_tau_sigma = 0.15))
# R2f (new_throughout only): R2 with the f block ROLLED BACK. R2-vs-R2f is the
# factorization proof -- the f terms enter generated quantities only, so every boat
# parameter that is not an f term must agree within Monte Carlo error and the shore must
# be bit-identical. It replaces the R3a/R3 pair at half the cost, and against the
# CURRENT turnover prior rather than the retired one.
D_R2f <- modifyList(D_R2,  F_OLD)
D_R3a <- modifyList(D_R2,  list(crab_fraction_strata = "month", crab_fraction_source = "both"))
D_R3  <- modifyList(D_R3a, list(crab_fraction_dynamic = TRUE))
# R4 is the shipped configuration: the derived shore turnover AND the 2026-09-11 census
# split (the commercial census exact, the charter expanded over the roster trip frame with
# its variance carried). The census is a PE-side quantity, so it does not disturb the
# fitted comparison between R3 and R4; it is here so R4 == run_config.
D_R4  <- modifyList(D_R3,  list(tau_shore_prior_mu = "derived", tau_shore_prior_sigma = "derived",
                                census_uncertainty = "charter", charter_frame = "roster", charter_expansion = "vessel"))

# The mode is in every tag and every output filename, so a "ladder" rung and a
# "new_throughout" rung can never land in the same folder or the same verdicts file.
.sfx <- if (identical(F_METHOD, "ladder")) "" else "-newf"
.tag <- function(x) paste0("IMP-", x, .sfx)
STAGE_DEFS <- list(
  # R0 reads the SHIPPED configuration (D_R4) so the desk rows describe what production
  # does, not what the R3 rung did: since 2026-09-11 that includes the census split (the
  # commercial census plus the charter expansion over the roster trip frame), and since
  # 2026-09-12 the PE's unsampled-cell levers.
  R0  = list(id = "R0",  model = "desk",   tag = .tag("R0-desk"),        delta = D_R4,
             headline = "desk: the PE fills and variance, the calibration turnover, the gear bootstrap, the census split, the f data, the holidays"),
  R1  = list(id = "R1",  model = "pooled", tag = .tag("R1-filters"),     delta = D_R1,
             headline = if (identical(F_METHOD, "ladder"))
               "pre-patch turnover priors + item 8 (unit-aware fishing-time filters)"
             else "pre-patch turnover priors + item 8 (filters) + the NEW f (dynamic monthly)"),
  R2  = list(id = "R2",  model = "pooled", tag = .tag("R2-tau-calib"),   delta = D_R2,
             headline = "+ item 3: boat turnover prior from the OSP/trailer calibration"),
  R2f = list(id = "R2f", model = "pooled", tag = .tag("R2f-legacy-f"),   delta = D_R2f,
             headline = "control: R2 with the f block rolled back (the factorization proof)"),
  R3a = list(id = "R3a", model = "pooled", tag = .tag("R3a-f-monthly"),  delta = D_R3a,
             headline = "+ item 1A: monthly f from the sampler contacts (legacy construction)"),
  R3  = list(id = "R3",  model = "pooled", tag = .tag("R3-f-dynamic"),   delta = D_R3,
             headline = "+ item 1B: the dynamic f"),
  R4  = list(id = "R4",  model = "pooled", tag = .tag("R4-shore-tau"),   delta = D_R4,
             headline = "+ item 2: shore turnover derived from the I/E time column (the shipped configuration)"),
  R5  = list(id = "R5",  model = "gear_resolved", tag = .tag("R5-gear-crosscheck"),
             delta = if (identical(GEAR_FOLLOWS, "R4")) D_R4 else D_R3,
             headline = sprintf("gear-resolved cross-check on the %s configuration", GEAR_FOLLOWS)))
if (!all(STAGES %in% names(STAGE_DEFS)))
  stop("STAGES names a stage that does not exist: ", paste(setdiff(STAGES, names(STAGE_DEFS)), collapse = ", "))
if (identical(F_METHOD, "ladder") && "R2f" %in% STAGES)
  stop("R2f is the new_throughout factorization control; under F_METHOD = 'ladder' that role is R3a vs R3.")
if (!identical(F_METHOD, "ladder") && any(c("R3a", "R3") %in% STAGES))
  stop("R3a / R3 are the 'ladder' f rungs; under F_METHOD = 'new_throughout' every rung already carries the new f (use R2f for the factorization proof).")

resolve_cfg <- function(sid) {
  cfg <- modifyList(BASE, WINDOW, keep.null = TRUE)
  modifyList(cfg, STAGE_DEFS[[sid]]$delta %||% list(), keep.null = TRUE)
}
# The KEYS any rung is allowed to differ in. Anything else differing between two rungs is
# a configuration leak, and the preflight fails on it rather than producing a comparison
# nobody can attribute.
DELTA_KEYS <- sort(unique(c(names(D_R1), names(D_R2), names(D_R2f), names(D_R3a), names(D_R3), names(D_R4),
                            "run_tag", "model", "tau_boat_prior_source", "tau_boat_prior_calibration_table",
                            "tau_shore_prior_source", "crab_fraction_rows", "osp_crab_rows",
                            "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates",
                            "ar_max_resolution", "shared_tau_min_obs", "tau_sensitivity_grid")))

# A digest of the resolved configuration, so RESUME can tell a folder built from THIS
# rung's config from one built before the deltas changed. Only the keys that can differ
# between rungs are hashed, plus the window, and the values are deparsed so a NULL and an
# absent key are distinguishable.
.cfg_fingerprint <- function(sid) {
  cfg <- resolve_cfg(sid)
  ks  <- sort(unique(c(DELTA_KEYS, names(WINDOW))))
  # DATA, not configuration: these are frames and date vectors read from the workbooks, so
  # they are identical across rungs by construction and would only make the digest depend
  # on row order.
  ks  <- setdiff(ks, c("crab_fraction_rows", "osp_crab_rows", "crabbing_holiday_dates",
                       "tau_boat_prior_calibration_table"))
  paste(c(paste0("model=", STAGE_DEFS[[sid]]$model %||% "?"),
          vapply(ks, function(k)
            paste0(k, "=", if (k %in% names(cfg)) paste(deparse(cfg[[k]]), collapse = " ") else "<absent>"),
            character(1))), collapse = "|")
}
# A plain position-weighted byte sum. No digest package is assumed (renv does not carry
# one) and cryptographic strength is not the point: the job is to notice that a folder was
# built from a DIFFERENT delta, and any collision-prone checksum does that.
stage_digest <- function(sid) {
  b <- as.integer(charToRaw(.cfg_fingerprint(sid)))
  sprintf("%s|%s|%08x", sid, F_METHOD,
          as.integer(sum(as.numeric(b) * seq_along(b)) %% 2147483647))
}
# THE CODE FINGERPRINT (2026-09-13). The config digest above says a folder was built
# from this rung's CONFIGURATION. It says nothing about the CODE. That gap matters for the
# two-pass plan Matt asked for -- fit R0/R1/R2/R4/R5 now, then come back and add R2f --
# because R2f's whole purpose is to be compared to R2, and if a Stan model, a driver or an
# R function changes between the two passes then R2f is fitted under different code than
# the R2 it is measured against, with nothing to flag it. Demonstrated before this was
# written: appending a line to crab_bss_pooled.stan left every digest unchanged and RESUME
# still reused R2.
#
# The fingerprint is three position-weighted byte sums, over the Stan models, the drivers
# and 03_R_functions, so a mismatch says WHICH layer moved. It is recorded, not enforced:
# re-fitting 12 h because a comment changed in an R function would be worse than the
# problem. What it does is (1) print the mismatch in the RESUME table and at the rung,
# (2) record a verdict row, and (3) make the cross-rung verdicts that depend on shared
# code (R2f's factorization proof, and every fit_exactness bit-identity claim) report
# REVIEW instead of PASS when the two folders disagree.
# COMMENTS AND BLANK LINES ARE STRIPPED BEFORE HASHING (2026-09-11). The first full run
# produced a record worth protecting, and the next documentation patch would otherwise have
# flagged every one of its folders as "code changed" and downgraded four PASS verdicts to
# REVIEW. An R comment cannot change a fit, so hashing it is a false positive generator --
# and a flag that cries wolf on doc edits is a flag nobody reads on the day it matters.
# This is a narrowing, not a weakening: any executable change still moves the hash.
.code_group <- function(dir_rel, pattern) {
  fs <- sort(list.files(.here(dir_rel), pattern = pattern, full.names = TRUE))
  if (!length(fs)) return("none")
  .strip <- function(f) {
    l <- readLines(f, warn = FALSE)
    l <- sub("(^|[^\\\\])#.*$", "\\1", l)     # trailing comments (Stan and R both use #; // handled below)
    l <- sub("//.*$", "", l)                  # Stan line comments
    l <- trimws(l)
    paste(l[nzchar(l)], collapse = "\n")
  }
  txt <- paste(vapply(fs, .strip, character(1)), collapse = "\n--\n")
  b <- as.integer(charToRaw(txt))
  sprintf("%08x", as.integer(sum(as.numeric(b) * seq_along(b)) %% 2147483647))
}
# PRE-B24 STAMPS, NORMALIZED FIRST. Until B24 (2026-09-11) the fingerprint hashed the RAW
# file text; it now strips comments and blank lines before hashing, so the current function
# can NEVER emit a pre-B24 string again and a layer-by-layer comparison against one is
# meaningless: every layer differs, including layers whose bytes never moved. The five
# 2026-09-11 rung folders carry pre-B24 stamps, so without this map they are permanently
# incomparable. Each entry maps a recorded pre-B24 stamp to code_fingerprint() COMPUTED ON
# THE SAME TREE. Verify a new entry, do not guess it:
#   git worktree add /tmp/wt <the commit the rung rendered from>   # then run code_fingerprint()
CODE_LEGACY <- list(
  # The five 2026-09-11 ladder rungs (R1, R2, R2f, R4, R5). Tree = 627a831, the last commit
  # before the run; 3609f1d touched only this runner, which is not a hashed layer. Checked
  # on that worktree: stan and drivers are bit-identical to HEAD and only 03_R_functions
  # has moved since (bss_day_length.R, D24).
  "stan:6e4aca2a drivers:3de636d0 fns:44079dee" = "stan:523f4e63 drivers:4c2ce454 fns:30ed14fb"
)
.code_norm <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(x)
  CODE_LEGACY[[x]] %||% x
}

# FINGERPRINT PAIRS DECLARED EQUIVALENT, with the reason, so a change that provably cannot
# alter a fit does not degrade an existing record. This is an AUDIT TRAIL, not an escape
# hatch: every entry names what changed and why the fits are unaffected, and an entry is only
# defensible when the default code path is identical. Anything you cannot write that sentence
# about belongs in a re-fit, not in this list.
#
# THE KEY NAMES BOTH ENDS, "<recorded> => <the fingerprint it is equivalent TO>", and that is
# load-bearing. Keyed on the recorded side alone -- as it was between 2026-09-11 and
# 2026-09-12 -- an entry NEVER EXPIRES: .code_delta short-circuits on the recorded string, so
# the declaration excuses that folder against whatever the tree later becomes. Measured, not
# argued: with the single-sided key, .code_delta("stan:6e4aca2a ...", a fingerprint with a
# changed Stan layer) returned "" -- a Stan edit would have gone unflagged on all five
# committed rungs, forever, and every cross-rung bit-identity claim would have kept its PASS.
# Naming both ends makes the declaration lapse the moment the current tree moves again, which
# is exactly when it should be re-examined.
CODE_EQUIVALENT <- list(
  # The 2026-09-11 full run (normalized by CODE_LEGACY above), against the tree as patched
  # the same day: the in-window I/E diagnostic in estimate_shore_turnover() (D24) and three
  # verdict-reporting fixes (D25, D26). The turnover derivation's DEFAULT path is
  # byte-identical: tau_shore_derive_window_only ships FALSE, under which `iv` is exactly the
  # frame the previous code built, and the new run_config key reaches no likelihood. The
  # verdict fixes are in this runner, which is not a hashed layer at all.
  "stan:523f4e63 drivers:4c2ce454 fns:30ed14fb => stan:523f4e63 drivers:4c2ce454 fns:0868556b" =
    "the 2026-09-11 ladder run vs the same tree plus D24's in-window turnover diagnostic; the default code path is byte-identical, so no fit is affected"
)
code_fingerprint <- function() {
  paste(sprintf("stan:%s", .code_group("02_stan_models", "\\.stan$")),
        sprintf("drivers:%s", .code_group("01_BSS_models", "\\.Rmd$")),
        sprintf("fns:%s", .code_group("03_R_functions", "\\.R$")), sep = " ")
}
.stage_stamp <- function(dir, sid) {
  writeLines(c(sprintf("stage: %s", sid), sprintf("F_METHOD: %s", F_METHOD),
               sprintf("digest: %s", stage_digest(sid)),
               sprintf("code: %s", code_fingerprint()),
               sprintf("rstan: %s / StanHeaders %s", utils::packageVersion("rstan"),
                       tryCatch(as.character(utils::packageVersion("StanHeaders")), error = function(e) "?")),
               sprintf("written: %s", format(Sys.time())),
               "", "# run_improvements_2026-09-08.R writes this after a rung renders. RESUME",
               "# reuses a folder ONLY when the `digest` line matches the stage it is resolving",
               "# now. The `code` line is NOT enforced: a mismatch is reported at the rung and in",
               "# the verdicts, and downgrades any cross-rung bit-identity or agreement claim that",
               "# rests on the two folders having been fitted by the same code."),
             file.path(dir, "IMP_STAGE.txt"))
}
.stamp_field <- function(dir, key) {
  p <- file.path(dir %||% "", "IMP_STAGE.txt")
  if (!file.exists(p)) return(NA_character_)
  l <- grep(paste0("^", key, ": "), readLines(p, warn = FALSE), value = TRUE)
  if (!length(l)) NA_character_ else sub(paste0("^", key, ": "), "", l[1])
}
.stage_digest_of <- function(dir) .stamp_field(dir, "digest")
.stage_code_of   <- function(dir) .stamp_field(dir, "code")
# Which layer(s) differ between a recorded fingerprint and another (or the current one).
.code_delta <- function(a, b = code_fingerprint()) {
  if (is.na(a %||% NA) || is.na(b %||% NA)) return(NA_character_)
  a <- .code_norm(a); b <- .code_norm(b)                       # pre-B24 stamps first
  if (identical(a, b)) return("")
  # declared equivalent; the key names BOTH ends, so the declaration lapses when either moves
  if (!is.null(CODE_EQUIVALENT[[paste(a, b, sep = " => ")]])) return("")
  pa <- strsplit(a, " ")[[1]]; pb <- strsplit(b, " ")[[1]]
  ka <- sub(":.*$", "", pa); kb <- sub(":.*$", "", pb)
  ks <- union(ka, kb)
  d <- ks[vapply(ks, function(k) !identical(pa[match(k, ka)], pb[match(k, kb)]), logical(1))]
  if (!length(d)) "" else paste(d, collapse = ", ")
}
# The one cross-pass risk the fingerprint cannot see: a package upgrade between passes.
# The driver already writes session_info.txt; this reads the one line that matters.
.stan_versions_of <- function(dir) {
  p <- file.path(dir %||% "", "session_info.txt")
  if (!file.exists(p)) return(NA_character_)
  l <- grep("^rstan .*StanHeaders", readLines(p, warn = FALSE), value = TRUE)
  if (!length(l)) NA_character_ else trimws(sub(";.*bss_seed.*$", "", l[1]))
}
# Used by every verdict that claims two folders are comparable. Returns "" when they are.
.comparability_note <- function(a, b) {
  cd <- .code_delta(.stage_code_of(a), .stage_code_of(b))
  vs <- if (!is.na(.stan_versions_of(a)) && !is.na(.stan_versions_of(b)) &&
            !identical(.stan_versions_of(a), .stan_versions_of(b)))
    sprintf("; rstan/StanHeaders differ (%s vs %s)", .stan_versions_of(a), .stan_versions_of(b)) else ""
  cc <- if (is.na(cd %||% NA)) "; one folder records no code fingerprint (pre-2026-09-13 run)"
        else if (nzchar(cd)) sprintf("; CODE DIFFERS between the two runs (%s)", cd) else ""
  paste0(cc, vs)
}

find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}
prev_pooled <- function(sid) {
  # the pooled rung this one is compared to; every optional rung is skipped when absent,
  # so R4 falls back to R2 under new_throughout and R3 falls back to R2 under "ladder".
  order <- c("R1", "R2", "R2f", "R3a", "R3", "R4")
  i <- match(sid, order); if (is.na(i) || i == 1) return(NA_character_)
  for (p in rev(order[seq_len(i - 1)])) if (p %in% STAGES && !identical(p, "R2f")) return(p)
  NA_character_
}

# ---------------------------------------------------------------------------
# PRE-FLIGHT
# ---------------------------------------------------------------------------
preflight <- function() {
  banner("PRE-FLIGHT")
  fails <- character(0)
  say <- function(ok, msg, detail = "") {
    cat(sprintf("  [%s] %s%s\n", if (ok) "ok" else "XX", msg,
                if (nzchar(detail)) paste0("  --  ", detail) else ""))
    if (!ok) fails <<- c(fails, msg)
  }
  say(identical(BASE$tau_boat_prior_mu, "calibration") && identical(BASE$crab_fraction_dynamic, TRUE) &&
        identical(BASE$crab_fraction_strata, "month") && identical(BASE$crab_fraction_source, "both") &&
        identical(BASE$tau_shore_prior_mu, "derived"),
      "run_config ships the R4 configuration (calibration turnover, dynamic monthly f from both sources, derived shore turnover)")
  say(identical(BASE$census_expansion, "none") && identical(BASE$census_uncertainty, "charter") &&
        identical(BASE$charter_frame, "roster") && identical(BASE$charter_expansion, "vessel"),
      paste("the commercial census is exact over the tally days and carries no error, and the charter part is an",
            "expansion whose SE enters the port interval (census_expansion none, census_uncertainty charter,",
            "charter_frame roster, charter_expansion vessel)"))
  say(is.null(BASE$ar_force) && !isTRUE(BASE$ar_escalate),
      "no experiment lever is active: ar_force NULL and ar_escalate off")
  for (m in c("crab_bss_pooled.stan", "crab_bss_gear_resolved.stan")) {
    src <- readLines(.here("02_stan_models", m), warn = FALSE)
    say(any(grepl("crab_fraction_dynamic", src, fixed = TRUE)) && any(grepl("cfi_kappa", src, fixed = TRUE)),
        sprintf("%s carries the dynamic-f block", m))
  }
  say(dir.exists(.here("05_output", REF$A1$dir)) &&
        file.exists(.here("05_output", REF$A1$dir, "pe_vs_bss_comparison.csv")),
      "the A1 baseline is on disk", REF$A1$dir)
  say(exists("fit_agreement") && exists("fit_exactness") && exists("merge_csv_by"),
      "the shared verdict helpers are sourced (batch_verdict_helpers.R)")
  say(GEAR_FOLLOWS %in% c("R3", "R4"), "GEAR_FOLLOWS names a pooled rung", GEAR_FOLLOWS)
  # The window is DELIBERATELY not run_config's. run_config ships the two-season span
  # (2023-24 + 2024-25); this ladder fits the single 2024-25 season so the rungs compare to
  # the 2026-09-04 A1 baseline and to each other, per Matt 2026-09-12. Said out loud here
  # so nobody reads "R4 == run_config" as including the window.
  say(!identical(BASE$season_filter, WINDOW$season_filter) || length(BASE$season_filter) == 1L,
      "the window is the ladder's pin, not run_config's",
      sprintf("ladder %s (%s to %s); run_config %s", paste(WINDOW$season_filter, collapse = "+"),
              WINDOW$est_date_start, WINDOW$est_date_end, paste(BASE$season_filter, collapse = "+")))
  say(identical(BASE$pe_empty_effort_stratum, "local_day_type") &&
        identical(BASE$pe_empty_stratum, "local") && identical(BASE$pe_variance, "impute_aware"),
      paste("run_config ships the 2026-09-12 PE unsampled-cell settings (local_day_type effort fill,",
            "local CPUE fill, impute_aware variance)"),
      sprintf("%s / %s / %s", BASE$pe_empty_effort_stratum %||% "?", BASE$pe_empty_stratum %||% "?",
              BASE$pe_variance %||% "?"))
  say(!identical(F_METHOD, "ladder") ||
        (identical(D_R1$crab_fraction_dynamic, FALSE) && identical(D_R1$crab_fraction_strata, "none")),
      "F_METHOD is consistent with D_R1's f block", sprintf("%s: R1 f = %s / %s / %s", F_METHOD,
        D_R1$crab_fraction_strata, D_R1$crab_fraction_source, D_R1$crab_fraction_dynamic))

  # ---- COMPARABILITY. Every rung must resolve to the SAME value for every key that is
  # not a declared delta, and each rung must differ from the previous one in the keys it
  # says it does and no others. This is the check the ladder's whole design rests on, and
  # until 2026-09-12 nothing verified it: WINDOW was applied to every rung and then each
  # rung's own values were compared only AFTER the fits, from run_parameters.txt, when the
  # MCMC had already been spent.
  .fit_stages <- setdiff(STAGES, "R0")
  if (length(.fit_stages) >= 2) {
    cfgs <- lapply(.fit_stages, resolve_cfg); names(cfgs) <- .fit_stages
    allk <- sort(unique(unlist(lapply(cfgs, names))))
    .same <- function(k) {
      v <- lapply(cfgs, function(c) c[[k]])
      all(vapply(v[-1], function(x) identical(x, v[[1]]), logical(1)))
    }
    leaked <- setdiff(allk[!vapply(allk, .same, logical(1))], DELTA_KEYS)
    say(!length(leaked),
        sprintf("no configuration leak: every non-delta key is identical across the %d fitted rungs", length(.fit_stages)),
        if (length(leaked)) paste("LEAKED:", paste(leaked, collapse = ", ")) else
          sprintf("%d keys compared", length(allk)))
    # the pins that reach the port total of a gate-failed component
    for (k in c("pe_empty_effort_stratum", "pe_empty_stratum", "pe_variance", "bss_seed",
                "est_date_start", "est_date_end", "season_filter", "pot_closure_start",
                "pot_closure_end", "pot_open_date", "census_start_date", "census_end_date"))
      say(.same(k), sprintf("pinned across every rung: %s", k),
          paste(deparse(cfgs[[1]][[k]]), collapse = " "))
    # Each rung differs from THE RUNG ITS VERDICT COMPARES IT TO, and only in declared
    # keys. The chain follows prev_pooled(), which skips the optional control rungs, so
    # R4-vs-R2 is checked under new_throughout and R2f-vs-R2 separately: comparing R4 to
    # R2f would report the f keys as a difference when nothing about f changed between the
    # rungs whose verdicts are actually read against each other.
    .pairs <- list()
    for (sid in intersect(c("R2", "R3a", "R3", "R4"), .fit_stages)) {
      pv <- prev_pooled(sid); if (!is.na(pv)) .pairs[[length(.pairs) + 1]] <- c(sid, pv)
    }
    if ("R2f" %in% .fit_stages && "R2" %in% .fit_stages) .pairs[[length(.pairs) + 1]] <- c("R2f", "R2")
    for (pr in .pairs) {
      a <- cfgs[[pr[2]]]; b <- cfgs[[pr[1]]]
      d <- union(names(a), names(b))
      d <- d[!vapply(d, function(k) identical(a[[k]], b[[k]]), logical(1))]
      say(all(d %in% DELTA_KEYS) && length(d) > 0,
          sprintf("%s differs from %s only in declared keys, and does differ", pr[1], pr[2]),
          sprintf("%d key(s): %s", length(d), paste(sort(d), collapse = ", ")))
    }
    # the manifest: what each rung WILL run, before any of it runs
    man <- do.call(rbind, lapply(.fit_stages, function(sid) {
      c1 <- cfgs[[sid]]
      data.frame(stage = sid, f_method = F_METHOD, model = STAGE_DEFS[[sid]]$model,
                 run_tag = STAGE_DEFS[[sid]]$tag, digest = stage_digest(sid),
                 tau_boat_prior_mu = paste(c1$tau_boat_prior_mu, collapse = ""),
                 tau_boat_prior_sigma = c1$tau_boat_prior_sigma %||% NA,
                 shared_tau_sigma = c1$shared_tau_sigma %||% NA,
                 tau_shore_prior_mu = paste(c1$tau_shore_prior_mu, collapse = ""),
                 f_strata = c1$crab_fraction_strata, f_source = c1$crab_fraction_source,
                 f_dynamic = c1$crab_fraction_dynamic,
                 census_uncertainty = c1$census_uncertainty, charter_frame = c1$charter_frame,
                 charter_expansion = c1$charter_expansion,
                 pe_effort_fill = c1$pe_empty_effort_stratum, pe_cpue_fill = c1$pe_empty_stratum,
                 pe_variance = c1$pe_variance,
                 headline = STAGE_DEFS[[sid]]$headline, stringsAsFactors = FALSE)
    }))
    mp <- .here("05_output", sprintf("improvements_2026-09-08_manifest%s.csv", .sfx))
    utils::write.csv(man, mp, row.names = FALSE)
    cat("\n  MANIFEST (written to ", basename(mp), "):\n", sep = "")
    print(man[, c("stage", "model", "tau_boat_prior_mu", "shared_tau_sigma", "tau_shore_prior_mu",
                  "f_strata", "f_dynamic", "census_uncertainty", "pe_effort_fill", "digest")], row.names = FALSE)
    # what an existing folder would do under RESUME
    if (isTRUE(RESUME)) {
      .reuse <- 0L
      for (sid in .fit_stages) {
        ex <- find_outdir(STAGE_DEFS[[sid]]$model, STAGE_DEFS[[sid]]$tag)
        if (!is.na(ex)) {
          dg <- .stage_digest_of(ex); keep <- identical(dg, stage_digest(sid))
          if (keep) .reuse <- .reuse + 1L
          cd <- if (keep) .code_delta(.stage_code_of(ex)) else NA_character_
          cat(sprintf("  RESUME: %s has an existing folder %s -- digest %s (%s)%s\n", sid, basename(ex),
                      if (is.na(dg)) "ABSENT" else dg,
                      if (keep) "matches: the fit will be SKIPPED"
                      else "does NOT match this stage's config: it will be RE-RUN",
                      if (!keep) "" else if (is.na(cd %||% NA)) "  [no code fingerprint recorded: pre-2026-09-13 run]"
                      else if (nzchar(cd)) sprintf("  *** CODE HAS CHANGED SINCE THAT FIT (%s): it is being REUSED anyway, and any cross-rung claim against it is downgraded ***", cd)
                      else if (!is.null(CODE_EQUIVALENT[[.stage_code_of(ex) %||% ""]]))
                        sprintf("  [code fingerprint declared EQUIVALENT: %s]", CODE_EQUIVALENT[[.stage_code_of(ex)]])
                      else "  [code fingerprint matches]"))
        }
      }
      if (.reuse > 0L)
        cat(sprintf("  RESUME will reuse %d of %d fitted rung(s); %s will be fitted.\n", .reuse, length(.fit_stages),
                    paste(setdiff(.fit_stages, .fit_stages[vapply(.fit_stages, function(sid) {
                      ex <- find_outdir(STAGE_DEFS[[sid]]$model, STAGE_DEFS[[sid]]$tag)
                      !is.na(ex) && identical(.stage_digest_of(ex), stage_digest(sid)) }, logical(1))]), collapse = ", ")))
    }
  }
  .self <- .here("06_diagnostics", "run_improvements_2026-09-08.R")
  .src <- if (file.exists(.self)) readLines(.self, warn = FALSE) else character(0)
  .src <- .src[!grepl("^\\s*#", .src)]
  .calls <- .src[grepl("rmarkdown::render(", .src, fixed = TRUE)]
  say(length(.calls) > 0 && !any(grepl("output_dir", .calls, fixed = TRUE)) &&
        any(grepl("cfg$run_tag <- st$tag", .src, fixed = TRUE)),
      "the tag goes INSIDE run_config and output_dir= is not passed to render()")
  say(any(grepl("file.copy(html", .src, fixed = TRUE)), "the rendered HTML is moved into the run folder")
  rule()
  if (length(fails)) {
    cat("  PRE-FLIGHT FAILED:\n"); for (f in fails) cat("   -", f, "\n")
    stop("Fix the above before running.")
  }
  cat("  pre-flight clean.\n")
}

# ---------------------------------------------------------------------------
# R0: the desk rung. Everything that changed without a fit, on the real inputs.
# ---------------------------------------------------------------------------
desk_R0 <- function() {
  banner("R0  desk: the PE fills and variance, the census split, and the data the fitted rungs read")
  out <- .here("05_output", format(Sys.Date(), "%Y%m%d"), "improvements-desk")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  q <- function(e) { s <- tempfile(); sink(s); on.exit(sink()); force(e) }
  ok <- tryCatch({
    p <- resolve_cfg("R0")
    p <- modifyList(p, list(bss_model_file = "crab_bss_pooled.stan", boat_require_gear_time = TRUE))
    p$ar_max_resolution <- p$ar_max_resolution$pooled
    p$crabbing_holiday_dates <- read_crabbing_holidays(p)
    dwg <- q(fetch_crab_data(p)); ie <- q(fetch_ie_data(p))
    p$crab_fraction_rows <- q(crab_fraction_source_rows(dwg, ie, p))
    st <- q(estimate_shore_turnover(attr(ie, "ie_intervals"), dwg$shore_effort, p))
    osp <- tryCatch(q(fetch_osp_boat_counts(p)), error = function(e) NULL)
    ov  <- tryCatch(q(diagnose_osp_trailer_overlap(osp, p, output_dir = out)), error = function(e) NULL)
    p$osp_crab_rows <- attr(osp, "osp_crab_rows")
    p <- q(bss_resolve_tau_boat_prior(p, ov))
    # item 5: the PE boat expands on the calibration turnover
    V1row("R0", "the boat PE expands the trailer count on the calibration turnover (item 5)",
          sprintf("tau_boat_prior_mu resolved to %.3f from '%s' (metric %s)", p$tau_boat_prior_mu,
                  p$tau_boat_prior_source %||% "?", p$tau_boat_calibration_metric %||% "?"),
          "a number near 3.0 from the overlap calibration, not 1.2",
          if (is.numeric(p$tau_boat_prior_mu) && p$tau_boat_prior_mu > 2 && grepl("calibration", p$tau_boat_prior_source %||% ""))
            "PASS" else "FAIL",
          paste("Until 2026-09-08 the PE expanded on the retired 1.2 while the BSS learned ~2.6 through",
                "tau_bar, so the boat PE-vs-BSS gap (3,709 vs 11,118 effort) was two turnovers",
                "disagreeing, not two estimators. The PE is the design-based cross-check and now is one."))
    # item 3: the prior centre the fitted rungs will use
    V1row("R0", "the boat turnover prior centre (item 3)",
          sprintf("lognormal(log %.3f, %.2f); shared_tau_sigma %s", p$tau_boat_prior_mu, p$tau_boat_prior_sigma %||% NA,
                  p$shared_tau_sigma %||% "NULL"),
          "centre 2.7-3.3 (the 61-day overlap gave 3.03 mean-per-visit), sigma 0.5", "INFO",
          paste("The previous prior lognormal(log 1.2, 0.3) sat 3+ prior SDs below the calibration and",
                "shrank tau_bar by 13% (posterior 2.60 vs likelihood-only 2.97). Re-centring removes",
                "the shrinkage; R2 measures what that does to the boat."))
    # item 1: the classification data
    bc <- dwg$boat_contacts
    bm <- p$crab_fraction_rows |> mutate(m = format(event_date, "%Y-%m")) |> group_by(m) |>
      summarise(t = sum(boats_total), c = sum(boats_crabbing), .groups = "drop")
    V1row("R0", "the crabbing-fraction classification stream exists (item 1)",
          sprintf("%d contacts on %d days, %d crabbing (share %.2f); by month: %s", sum(bc$boats_total), nrow(bc),
                  sum(bc$boats_crabbing), sum(bc$boats_crabbing) / max(sum(bc$boats_total), 1),
                  paste(sprintf("%s %.2f (n=%d)", bm$m, bm$c / pmax(bm$t, 1), bm$t), collapse = ", ")),
          "hundreds of contacts, strongly seasonal (0.9+ in winter, ~0.3 in summer)",
          if (sum(bc$boats_total) > 200) "PASS" else "REVIEW",
          paste("These rows were always in the interview workbook and were dropped by the reader's",
                "crabbers > 0 filter. They are the f data: R3a and R3 read them."))
    # item 2: the derived shore turnover
    V1row("R0", "the shore turnover at the count hours (item 2; R4 adopts it)",
          sprintf("tau_shore %.3f (peak-based %.3f) over %d I/E days, bootstrap log-SE %.3f, weekday %s / weekend %s",
                  st$tau, st$tau_peak %||% NA, st$n_days, st$log_se,
                  fmt(.num1(st$by_day_type$tau_ratio_of_sums[st$by_day_type$dt == "weekday"]), 2),
                  fmt(.num1(st$by_day_type$tau_ratio_of_sums[st$by_day_type$dt == "weekend"]), 2)),
          "about 2.5 against the 1.7 prior in use", "INFO",
          paste("The shore counts are taken at 9:00-14:00, when presence is ~0.7 of the daily peak, so",
                "the peak-based 1.7 under-expands the count. This is the largest single mover in the",
                "series (shore x ~1.47) and it ships OFF until R4 shows what it does."))
    # item 4: the census split
    cc <- q(estimate_comm_charter(dwg, p))
    utils::write.csv(cc$daily_full, file.path(out, "census_daily.csv"), row.names = FALSE)
    utils::write.csv(cc$variance_detail, file.path(out, "census_variance.csv"), row.names = FALSE)
    cc_tally <- q(estimate_comm_charter(dwg, modifyList(p, list(charter_frame = "tally", charter_expansion = "pooled"))))
    cc_pool  <- q(estimate_comm_charter(dwg, modifyList(p, list(charter_expansion = "pooled"))))
    V1row("R0", "the commercial CENSUS and the charter EXPANSION, separated (item 4, settled 2026-09-09; the charter corrected 2026-09-11)",
          sprintf(paste0("total %s = commercial %s (census, census_expansion = '%s': %.0f vessel-trips on %d tally days, %d unsampled calendar days %s; ",
                         "per-vessel mean over %d interviews, SE %s = %.1f%%, not carried) + charter %s (expansion over the %s trip frame, %s: %.0f trips, ",
                         "%d interviewed = %.0f%%, %s crab observed, SE %s = %.1f%%, CARRIED). Charter alternatives: %s pooled, %s on the tally frame. ",
                         "Baseline %s under the 2026-09-08 day-type expansion of the whole component."),
                  fmt(cc$Dungeness_Kept, 0), fmt(cc$commercial_dung, 0), cc$census_expansion %||% "none", cc$commercial_vessels,
                  sum(cc$daily_full$observed), sum(!cc$daily_full$observed),
                  if (identical(cc$census_expansion, "none")) "at zero" else sprintf("imputed at %s", fmt(cc$imputed_dung, 0)),
                  cc$commercial_interviews %||% NA_integer_,
                  fmt(cc$commercial_se, 0), 100 * cc$commercial_se / max(cc$commercial_dung, 1),
                  fmt(cc$charter_dung, 0), cc$charter_frame %||% "tally", cc$charter_expansion %||% "pooled", cc$charter_trips,
                  cc$charter_interviews, 100 * (cc$charter_sampled_frac %||% NA_real_), fmt(cc$charter_observed_dung, 0),
                  fmt(cc$charter_se, 0), 100 * cc$charter_se / max(cc$charter_dung, 1),
                  fmt(cc_pool$charter_dung, 0), fmt(cc_tally$charter_dung, 0), fmt(REF$A1$census, 0)),
          "the commercial part is the tally-day sum with no carried error; the charter part is an expansion whose SE enters the port interval",
          if (identical(cc$census_expansion, "none") && isTRUE(all.equal(cc$imputed_dung, 0)) &&
              identical(cc$census_uncertainty, "charter") && (cc$charter_se %||% 0) > 0) "PASS" else "REVIEW",
          paste("Matt 2026-09-09: samplers are scheduled on the days the commercial (recreational) vessels are",
                "confirmed to be operating, so a window day without a tally is a day with no commercial fishing.",
                "Matt 2026-09-11: that census-without-error statement applies to the COMMERCIAL boats only; the",
                "charter vessels are not 100% sampled and need expansion. The charter trip roster is the frame,",
                "and on 2024-25 a quarter of its sailed trips fell on days the tally never saw."))
    .cal <- seq(as.Date(p$est_date_start), as.Date(p$est_date_end), by = "day")
    .hol <- .cal %in% p$crabbing_holiday_dates
    .wke <- weekdays(.cal) %in% p$days_wkend
    V1row("R0", "the crabbing-holiday calendar gained Thanksgiving Day, Veterans Day (observed) and Juneteenth (2026-09-11)",
          sprintf("%d holiday(s) in the window, %d of them on a weekday (so they re-type a day); day types: %d weekday, %d weekend, %d holiday of %d",
                  sum(.hol), sum(.hol & !.wke), sum(!.hol & !.wke), sum(!.hol & .wke), sum(.hol), length(.cal)),
          "three more days per season carry the weekend + holiday effort effects (B1 + B2) in the BSS",
          "READ",
          paste("The three were added on the samplers' own Holiday? flag and the Float 20 gear count against the",
                "same month's weekday mean (Thanksgiving 2.4x, Veterans Day observed 2.3x, Juneteenth 1.7-2.3x).",
                "In the PE, an unsampled holiday in a week with no other sampled holiday lands in an EMPTY",
                "(week x day-type) stratum, which the shipped pe_empty_effort_stratum = 'zero' expands at ZERO",
                "effort: that is why the 2024-25 PE port total falls 0.5% when the three are added. The BSS",
                "imputes those days from B1 + B2 and does not have the problem. See CHANGE_REGISTER D19."))
    # 2026-09-09: the combo-trip share from the trip types, and the shift coverage
    fs0 <- attr(crab_fraction_stan_data(FALSE, tibble(event_date = seq(as.Date(p$est_date_start), as.Date(p$est_date_end), by = "day")), p, quiet = TRUE), "f_strata")
    V1row("R0", "the combo-trip share c is observed from the contacts' trip types (item 1, 2026-09-09)",
          sprintf("%d typed crabbing boats, %d combos (%.2f); by month: %s", sum(fs0$typed_crabbing), sum(fs0$typed_combo),
                  sum(fs0$typed_combo) / max(sum(fs0$typed_crabbing), 1),
                  paste(sprintf("%s %s", fs0$label, ifelse(fs0$typed_crabbing > 0, sprintf("%.2f (n=%d)", fs0$typed_combo / pmax(fs0$typed_crabbing, 1), fs0$typed_crabbing), "-")), collapse = ", ")),
          "none in Dec-Feb, half or more in the salmon and bottomfish months", if (sum(fs0$typed_crabbing) > 100) "PASS" else "REVIEW",
          paste("c does not move the boat total (f does); it is what OSP's crabbing-only count reads low by,",
                "so f(1 - c) is the prediction of that column and the comparison, when OSP delivers it, is the",
                "check of the shift-time contacts against an all-day count. In summer more than half the",
                "crabbing boats are combos, which is why OSP's column could never have identified f alone."))
    # which interviews are the classification: the launches (shipped), the launches plus the
    # docks, or every private-boat interview. The choice moves the summer f, so it is on
    # the record here until the protocol question (are dock-interviewed private boats
    # trailered or moored?) is answered.
    .area_sets <- list(launches = p$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch"),
                       launches_and_docks = c(p$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch"),
                                              p$shore_dock_float20 %||% "Westport Docks Float 20", p$shore_dock_float17 %||% "Westport Docks Float 17-21"),
                       all = "all")
    .gh_int <- tryCatch(q({ pp <- p; pp$crab_fraction_contact_areas <- "all"; fetch_crab_data(pp)$boat_contacts }), error = function(e) NULL)
    if (!is.null(.gh_int)) {
      .sens <- do.call(rbind, lapply(names(.area_sets), function(nm) {
        pp <- p; pp$crab_fraction_contact_areas <- .area_sets[[nm]]
        b <- q(fetch_crab_data(pp)$boat_contacts) |> mutate(m = format(event_date, "%Y-%m")) |> group_by(m) |>
          summarise(t = sum(boats_total), c = sum(boats_crabbing), .groups = "drop")
        data.frame(area_set = nm, month = b$m, contacts = b$t, crabbing_share = round(b$c / pmax(b$t, 1), 3))
      }))
      utils::write.csv(.sens, file.path(out, "contact_area_sensitivity.csv"), row.names = FALSE)
      .w <- tidyr::pivot_wider(.sens |> select(area_set, month, crabbing_share), names_from = area_set, values_from = crabbing_share)
      V1row("R0", "which interviews are the f classification: launch sites (shipped) vs launches + docks vs all",
            paste(sprintf("%s %s", .w$month, apply(.w[, -1], 1, function(r) paste(sprintf("%.2f", r), collapse = "/"))), collapse = "; "),
            "the summer months differ (docks are finfish-heavy, the marina is crab-only moorage); the protocol decides", "READ",
            paste("Order in each triple: launches / launches+docks / all. The trailer and OSP counts measure boats",
                  "launched at the ramp. A private boat interviewed at the docks is either a trailered boat that",
                  "stopped to unload (then it belongs in the classification and excluding it biases f up in",
                  "summer) or a moored boat (then it is outside the counted population). The marina rows are",
                  "crab-only moorage boats either way. crab_fraction_contact_areas selects the set."))
    }
    sh <- tryCatch(q(fetch_sampler_shifts(p)), error = function(e) NULL)
    sc <- tryCatch(q(diagnose_shift_coverage(sh, ie, dwg$boat_contacts_detail, p, output_dir = out)), error = function(e) NULL)
    if (!is.null(sc)) {
      cb <- sc$crabbing_by_hour
      V1row("R0", "sampler shift coverage of the day's boat returns (item 1, shift times)",
            sprintf("%d shift days; in port %s to %s (medians), %.1f h; %.0f%% of a day's boat returns inside the shift (%s), %.0f%% after the last check-out; crabbing share of typed contacts by hour: %s",
                    sc$summary$n_shift_days, .hhmm_of(sc$summary$median_first_check_in), .hhmm_of(sc$summary$median_last_check_out),
                    sc$summary$median_hours_in_port, 100 * sc$summary$median_return_share_covered, sc$summary$return_profile_source,
                    100 * sc$summary$returns_after_last_check_out,
                    if (is.null(cb)) "n/a" else paste(sprintf("%s %.2f (n=%d)", cb$hour_bin, cb$crabbing_share, cb$contacts), collapse = ", ")),
            "flat crabbing share across the hours inside the shift; the share of returns outside it is the exposure", "INFO",
            paste("The contacts classify the boats returning during the shift. Inside it the trip-type mix does",
                  "not drift with the hour (2024-25), which is the evidence available today; the returns outside",
                  "the shift are unclassified and only OSP's all-day count can say whether their mix differs.",
                  "shift_coverage_*.csv and contact_hour_by_trip_type.csv carry the detail."))
    }
    # ---- 2026-09-12: THE PE ITSELF, under every unsampled-cell setting. -------------
    # R0 is called "the PE and reporting changes" and until now never ran the PE. It costs
    # a minute per arm and it is the only way the D19 fill question can be settled without
    # spending 4 h of MCMC per arm; the shipped BSS reference for the comparison is R4,
    # which this same run produces.
    # The desk PE must be the DRIVER's PE, so it needs the same two resolutions the driver
    # does before section 4: the shore turnover prior (D_R4 sets tau_shore_prior_mu =
    # "derived", a STRING that run_pe_pooled cannot multiply) and the L_effective
    # regression. p already carries the resolved boat turnover from bss_resolve_tau_boat_prior
    # above.
    p <- q(bss_resolve_tau_shore_prior(p, st))
    .Lm <- if (isTRUE(p$use_ie_day_length) && nrow(ie) > 0)
             tryCatch(q(estimate_L_effective(ie, p)), error = function(e) NULL) else NULL
    .pe_arm <- function(pp) {
      subs <- build_subseasons(pp)
      out <- list()
      for (pop in c("shore", "private_boat")) for (ss in subs) {
        lab <- paste0(pop, "_", ss$name)
        sm <- q(prep_population_summary(dwg, pop, ss$start, ss$end, pp))
        dd <- q(prep_days_crab(ss$start, ss$end, pp, L_eff_model = .Lm))
        out[[lab]] <- tryCatch(q(run_pe_pooled(sm, dd, pp, lab)), error = function(e) NULL)
      }
      out
    }
    .arms <- list(
      "effort zero, CPUE pooled (through 2026-09-11)" = list(pe_empty_effort_stratum = "zero", pe_empty_stratum = "pooled"),
      "effort day_type, CPUE local"                    = list(pe_empty_effort_stratum = "day_type", pe_empty_stratum = "local"),
      "effort local_day_type, CPUE pooled"             = list(pe_empty_effort_stratum = "local_day_type", pe_empty_stratum = "pooled"),
      "effort local_day_type, CPUE local (SHIPPED)"    = list(pe_empty_effort_stratum = "local_day_type", pe_empty_stratum = "local"))
    .rows <- list(); .pe_ship <- NULL
    for (nm in names(.arms)) {
      pa <- modifyList(p, .arms[[nm]])
      res <- .pe_arm(pa)
      if (grepl("SHIPPED", nm)) .pe_ship <- res
      for (k in names(res)) {
        r <- res[[k]]; if (is.null(r)) next
        .rows[[length(.rows) + 1]] <- data.frame(
          arm = nm, component = k,
          effort_fill = r$pe_empty_effort_fill %||% NA, cpue_fill = pa$pe_empty_stratum %||% NA,
          variance = r$pe_variance %||% NA,
          effort = r$effort_total %||% NA, effort_se = r$effort_se %||% NA,
          effort_se_sampled_only = r$effort_se_sampled_only %||% NA,
          catch = r$Dungeness_Kept %||% NA, imputed_catch = r$imputed_Dungeness_Kept %||% NA,
          n_empty_strata = r$n_empty_effort_strata %||% NA, n_empty_days = r$n_empty_effort_days %||% NA,
          n_single_strata = r$n_single_effort_strata %||% NA, n_single_days = r$n_single_effort_days %||% NA,
          n_strata_total = r$n_effort_strata_total %||% NA,
          n_calendar_days = r$n_calendar_days %||% NA, zeroed_effort_bias = r$pe_zeroed_effort_bias %||% NA,
          stringsAsFactors = FALSE)
      }
    }
    .pe_tab <- do.call(rbind, .rows)
    utils::write.csv(.pe_tab, file.path(out, "pe_unsampled_cell_arms.csv"), row.names = FALSE)
    # The commercial/charter component is fill-invariant (it is a census plus an expansion,
    # not a stratified day expansion), so add it back as a constant: these totals are then
    # directly comparable to pe_port_summary.csv from any fitted rung.
    .cc_c <- cc$Dungeness_Kept %||% 0; .cc_e <- cc$effort_total %||% 0
    .tot <- .pe_tab |> group_by(arm) |> summarise(catch = sum(catch, na.rm = TRUE) + .cc_c,
                                                  effort = sum(effort, na.rm = TRUE) + .cc_e, .groups = "drop")
    .tot <- .tot[match(names(.arms), .tot$arm), ]          # the declared order, not alphabetical
    .base_c <- .num1(.tot$catch[grepl("^effort zero", .tot$arm)])
    V1row("R0", "the PE's unsampled-cell settings, all four arms on the real 2024-25 inputs (D19; the CPUE half found 2026-09-12)",
          paste(sprintf("%s: catch %s (%+.1f%% vs the retired arm), effort %s", .tot$arm, fmt(.tot$catch, 0),
                        100 * (.tot$catch - .base_c) / max(.base_c, 1), fmt(.tot$effort, 0)), collapse = "; "),
          "the shipped arm is the one the fitted rungs use; the spread between arms is the size of the D19 question",
          "READ",
          paste("The 2026-09-11 figure for this was +26% on the PE port total, and it was WRONG in a",
                "specific way: it paired a month-local EFFORT fill with a sub-season-pooled CPUE, so the",
                "seasonal gradient was counted once in the effort and once more in the rate. The imputed",
                "boat cells are in the summer and the pooled boat CPUE is dominated by the high-CPUE",
                "winter. Putting both halves on the month scale takes the boat all-gear PE from 42,841",
                "to 37,018 and the port total from 90,861 to 85,076, +17.8% on the retired 72,224.",
                "Read the arm spread against R4's BSS, not against 20260904 A1, which predates the",
                "derived shore turnover and the boat recentring."))
    if (!is.null(.pe_ship)) {
      .sh <- .pe_tab[grepl("SHIPPED", .pe_tab$arm), ]
      V1row("R0", "how much of the PE rests on 0 or 1 sampled days, and what that does to its SE (D21)",
            paste(sprintf("%s: %d of %d cells unsampled (%d days) and %d singleton (%d days) = %.0f%% of days; SE %s vs %s on the retired arithmetic",
                          .sh$component, .sh$n_empty_strata, .sh$n_strata_total,
                          .sh$n_empty_days, .sh$n_single_strata, .sh$n_single_days,
                          100 * (.sh$n_empty_days + .sh$n_single_days) / pmax(.sh$n_calendar_days, 1),
                          fmt(.sh$effort_se, 0), fmt(.sh$effort_se_sampled_only, 0)), collapse = "; "),
            "the SE must widen when half the days rest on one observation or none",
            if (all(.sh$effort_se >= .sh$effort_se_sampled_only, na.rm = TRUE) &&
                any(.sh$effort_se > .sh$effort_se_sampled_only * 1.5, na.rm = TRUE)) "PASS" else "REVIEW",
            paste("The retired SE was sqrt(N^2 sd^2 / max(n,1)) with sd from the cell's own sampled days,",
                  "and sd() of ONE observation is NA, replaced by 0. So 37 of 93 shore all-gear cells and",
                  "34 of 93 boat cells contributed their full point estimate and no variance, and so did",
                  "every imputed cell: the shore SE was bit-identical at 410 whether the fill added 0 or",
                  "7,947 effort units. This is a PE-side interval, but a gate-failed component reports its",
                  "PE point in the port total as a CONSTANT, which is the case where it reaches the",
                  "headline. pe_variance = 'sampled_only' restores the old arithmetic exactly."))
    }

    # item 6: the gear bootstrap on the boat all-gear interviews
    sub <- build_subseasons(p)
    ss  <- sub[[which(vapply(sub, function(x) x$gear_regime == "all_gear", logical(1)))]]
    iv  <- dwg$interview |> filter(population == "private_boat", event_date >= ss$start, event_date <= ss$end) |>
      transmute(gear_primary = gear_primary_class(gear_type), catch = as.numeric(dungeness_kept))
    shares <- gear_share_bootstrap(iv, n_boot = 2000, seed = 1L)
    pt <- gear_share_point(iv)
    qs <- apply(shares, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    V1row("R0", "gear shares by trip-level bootstrap (item 6)",
          sprintf("boat all-gear, %d trips: %s", nrow(iv),
                  paste(sprintf("%s %.3f [%.3f, %.3f]", names(pt), pt, qs[1, ], qs[2, ]), collapse = "; ")),
          "intervals several times wider than the retired Dirichlet on crab counts", "INFO",
          paste("Crab arrive in trips, not one at a time; the Dirichlet(catch + 0.5) treated every crab",
                "as an independent draw. The report now applies per-population x sub-season shares."))
    cat(sprintf("  R0 written to %s\n", out))
    TRUE
  }, error = function(e) {
    cat("  R0 desk rung unavailable:", conditionMessage(e), "\n")
    V1row("R0", "desk rung", paste("not evaluated:", conditionMessage(e)), "runs to completion", "REVIEW",
          "Needs the data-reading packages and the input files; it is not a Stan check.")
    NA
  })
  invisible(ok)
}

# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------
run_stage <- function(sid) {
  st  <- STAGE_DEFS[[sid]]; cfg <- resolve_cfg(sid)
  if (identical(st$model, "desk")) { desk_R0(); return(NA_character_) }
  banner(sprintf("%s  %s", sid, st$headline))
  cat("  run_tag:", st$tag, "\n")
  cat(sprintf("  config : tau_boat_prior_mu %s | shared_tau_sigma %s | f strata %s source %s dynamic %s | tau_shore_prior_mu %s | census %s\n",
              paste(cfg$tau_boat_prior_mu, collapse = ""), cfg$shared_tau_sigma %||% "NULL",
              cfg$crab_fraction_strata, cfg$crab_fraction_source, cfg$crab_fraction_dynamic,
              paste(cfg$tau_shore_prior_mu, collapse = ""), cfg$census_uncertainty))
  existing <- find_outdir(st$model, st$tag)
  if (isTRUE(RESUME) && !is.na(existing) && file.exists(file.path(existing, "run_parameters.txt"))) {
    # 2026-09-12: a matching DIGEST is now required, not just a finished-looking folder.
    # The deltas in this file changed on 2026-09-11 and again on 2026-09-12, so reusing a
    # folder on the strength of its name would have silently mixed configurations into one
    # ladder. run_parameters.txt is still checked first: it is written last by the driver,
    # so its presence is the completion marker.
    dg <- .stage_digest_of(existing)
    if (identical(dg, stage_digest(sid))) {
      .cd <- .code_delta(.stage_code_of(existing))
      cat("  RESUME: output present at", basename(existing), "with a MATCHING config digest - skipping the fit.\n")
      if (is.na(.cd %||% NA))
        cat("          (no code fingerprint recorded in that folder: it predates 2026-09-13.)\n")
      else if (nzchar(.cd)) {
        cat(sprintf(paste0("          *** THE CODE HAS CHANGED SINCE THAT FIT (%s). The folder is REUSED, because\n",
                           "          re-fitting hours of MCMC over a code edit that may be a comment is worse than\n",
                           "          the problem, but any verdict comparing this rung to another is downgraded. ***\n"), .cd))
        V1row(sid, "reused fit was produced by DIFFERENT code than this run",
              sprintf("folder %s; layer(s) changed: %s", basename(existing), .cd),
              "the reused folder's code fingerprint matches the current tree", "REVIEW",
              paste("RESUME matched this rung's CONFIG digest and reused the fit. The code fingerprint",
                    "(Stan models / drivers / 03_R_functions) does not match, so this rung was fitted by",
                    "a different version of the pipeline than the rungs fitted in this pass. Cross-rung",
                    "bit-identity and agreement claims involving it are not interpretable; either accept",
                    "them as indicative or delete the folder and re-fit."))
      }
      return(existing)
    }
    cat(sprintf(paste0("  RESUME: output present at %s but its config digest %s does not match this stage (%s).\n",
                       "          RE-RUNNING. The old folder is overwritten by the render.\n"),
                basename(existing), if (is.na(dg)) "is ABSENT (pre-2026-09-12 run)" else dg, stage_digest(sid)))
  }
  if (isTRUE(DRY_RUN)) { cat("  DRY_RUN: not fitting.\n"); return(NA_character_) }
  cfg$model <- st$model; cfg$run_tag <- st$tag; cfg$run_weather <- FALSE
  run_env <- new.env(parent = globalenv()); run_env$run_config <- cfg
  t0 <- Sys.time()
  html <- rmarkdown::render(model_rmd[[st$model]], envir = run_env, quiet = FALSE)
  od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e) NA_character_)
  if (!is.na(od) && dir.exists(od) && file.exists(html) &&
      normalizePath(dirname(html)) != normalizePath(od)) {
    if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE)))
      suppressWarnings(file.remove(html))
    else cat("  WARNING: could not move the rendered HTML; the NEXT stage will overwrite it.\n")
  }
  done <- find_outdir(st$model, st$tag)
  cat(sprintf("  %s finished in %.1f min -> %s\n", sid,
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              if (is.na(done)) "OUTPUT FOLDER NOT FOUND" else done))
  if (is.na(done)) stop(sprintf("%s rendered but no folder named %s%s exists under 05_output",
                                sid, prefix[[st$model]], st$tag))
  .stage_stamp(done, sid)   # the config digest RESUME will check on a later pass
  done
}

# ---------------------------------------------------------------------------
# VERDICTS. Every pooled rung records its components; the per-rung blocks read the
# change against the previous rung and the baseline.
# ---------------------------------------------------------------------------
LAD <- list()
ladder_row <- function(sid, dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  pt <- .port_row(rd(dir, "port_total_Dungeness_Kept.csv"))
  tb <- .full_row(dir, "private_boat_all_gear", "tau_bar_out")
  sf <- .full_row(dir, "private_boat_all_gear", "sigma_f_out")
  LAD[[sid]] <<- data.frame(
    rung = sid, folder = basename(dir),
    shore_pc = .comp(dir, "shore (Pot closure)"), shore_ag = .comp(dir, "shore (All gear)"),
    boat_pc = .comp(dir, "private_boat (Pot closure)"), boat_ag = .comp(dir, "private_boat (All gear)"),
    census = .comp(dir, "comm_charter (census)", "PE_catch"),
    port = .num1(pt$BSS_median), port_lo95 = .num1(pt$BSS_lo95), port_hi95 = .num1(pt$BSS_hi95),
    tau_bar = .num1(tb$mean), sigma_f = .num1(sf$mean),
    stringsAsFactors = FALSE)
}
.pct <- function(a, b) if (isTRUE(is.finite(a)) && isTRUE(is.finite(b)) && b != 0) 100 * (a - b) / b else NA_real_

verdict_R1 <- function(dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  base <- .here("05_output", REF$A1$dir)
  sa <- .comp(dir, "shore (All gear)"); sp <- .comp(dir, "shore (Pot closure)")
  ba <- .comp(dir, "private_boat (All gear)"); bp <- .comp(dir, "private_boat (Pot closure)")
  # 2026-09-12: under F_METHOD = "new_throughout" this rung also carries the NEW f, so the
  # BOAT is EXPECTED to move against the baseline and cannot be a pass/fail criterion. The
  # SHORE still can: f is boat-only (crab_fraction_stan_data() sets apply_cf = 0 for shore),
  # so the shore is the clean regression test against A1 in either mode.
  .newf <- !identical(F_METHOD, "ladder")
  V1row("R1", if (.newf)
          "item 8 + the new f: the SHORE against the baseline (the boat is expected to move with f)"
        else "item 8 alone (plus the 2026-09-09 workbook): the fitted components against the baseline",
        sprintf("shore pc %s (%+.1f%%), shore ag %s (%+.1f%%), boat pc %s (%+.1f%%), boat ag %s (%+.1f%%); census %s (baseline %s: exact over the tally days by design)",
                fmt(sp, 0), .pct(sp, REF$A1$shore_pc), fmt(sa, 0), .pct(sa, REF$A1$shore_ag),
                fmt(bp, 0), .pct(bp, REF$A1$boat_pc), fmt(ba, 0), .pct(ba, REF$A1$boat_ag),
                fmt(.comp(dir, "comm_charter (census)", "PE_catch"), 0), fmt(REF$A1$census, 0)),
        if (.newf) "both shore components within a few percent; the boat rises with the winter f (no threshold)"
        else "every fitted component within a few percent; the boat a few percent LOWER (zero-catch trips restored)",
        if (.newf) {
          if (all(abs(c(.pct(sp, REF$A1$shore_pc), .pct(sa, REF$A1$shore_ag))) < 5, na.rm = TRUE)) "PASS" else "REVIEW"
        } else if (all(abs(c(.pct(sp, REF$A1$shore_pc), .pct(sa, REF$A1$shore_ag))) < 5, na.rm = TRUE) &&
              isTRUE(.pct(ba, REF$A1$boat_ag) < 2)) "PASS" else "REVIEW",
        paste("The hours filter dropped 22 boat trips (14 with recorded catch) and ~140 shore rows whose",
              "time field was blank or under 0.5 h; under the deployment unit a set pot with zero time",
              "is real effort with real catch. Restoring them lowers the boat CPUE slightly (more zeros)",
              "and barely moves the shore. The rebuilt workbook (2026-09-09) also applies the tampered-gear",
              "filter (13 Grays Harbor interviews in 2024-25) and a few source corrections. The census is",
              "not compared: it is now the exact sum over the tally days (-3,869 by design).",
              if (.newf) paste("UNDER new_throughout the boat ALSO carries the dynamic monthly f, whose",
                               "winter months run 0.88 to 1.00 against the retired 0.30 anchor, so the boat",
                               "moving several tens of percent here is the design and not a regression. The",
                               "f change is not attributed by any fit in this mode; add R2f for the",
                               "factorization proof or run F_METHOD = 'ladder' for the full attribution.") else ""))
  b <- .full_row(dir, "private_boat_all_gear", "tau_bar_out")
  V1row("R1", "the boat turnover under the pre-patch prior (the R2 control)",
        sprintf("tau_bar %s (baseline %s)", fmt(.num1(b$mean), 4), fmt(REF$A1$tau_bar, 4)),
        "about 2.6: the calibration says ~3.0 and the 1.2 prior is pulling it down", "INFO",
        "This is the number R2 exists to move.")
}

verdict_R2 <- function(dir, prev) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  fe <- fit_exactness(dir, prev %||% "", pat = "shore", what = "shore fits vs R1",
                      expect_delta = c("tau_boat_prior_mu", "tau_boat_prior_sigma", "shared_tau_sigma", "run_tag", "model",
                                       "tau_boat_prior_source", "tau_boat_prior_calibration_table", "crabbing_holiday_dates",
                                       "opener_f_dates", "razor_dig_dates", "tau_sensitivity_grid"))
  V1cross("R2", "the shore did not move (the boat prior cannot reach it)", fe$observed,
        "shore fits bit-identical to R1", fe$verdict,
        "A FAIL means the turnover prior leaked into the shore fits, which it has no path to do.",
        dir, prev %||% "")
  tb <- .full_row(dir, "private_boat_all_gear", "tau_bar_out")
  ba <- .comp(dir, "private_boat (All gear)"); ba0 <- .comp(prev %||% "", "private_boat (All gear)")
  # 2026-09-11 CORRECTION, from the first full run. The threshold here was 8-20% on the
  # boat, reasoned from "effort = lambda_E x tau x f, so the boat moves in proportion to
  # tau_bar". THAT PROPORTIONALITY DOES NOT HOLD, and the run is what showed it: tau_bar
  # rose 14.0% (2.612 -> 2.977) and boat all-gear rose 4.4% (43,668 -> 45,604), because
  # the latent effort level absorbed most of it (mu_mu_E 2.7067 -> 2.6449 on the log
  # scale, a -6.0% multiplicative fall) while kappa_OSP held still (3.120 -> 3.136). The
  # boat effort scale is JOINTLY identified by the trailer counts, the OSP counts and the
  # catch likelihood, so re-centring the prior moves the reported component by roughly a
  # third of the prior move, not by all of it. That is a better property than the one the
  # old threshold assumed -- it means the OSP stream is doing work -- but it made a correct
  # run report REVIEW. The band is widened to 2-20% and the pass-through ratio is now
  # REPORTED, because the ratio is the quantity worth watching: a value near 1 would mean
  # the effort streams had stopped constraining the level.
  .tb0 <- .num1(.full_row(prev %||% "", "private_boat_all_gear", "tau_bar_out")$mean)
  .tpass <- if (isTRUE(is.finite(.tb0)) && .tb0 > 0 && isTRUE(is.finite(.pct(.num1(tb$mean), .tb0))) &&
                abs(.pct(.num1(tb$mean), .tb0)) > 1e-9) .pct(ba, ba0) / .pct(.num1(tb$mean), .tb0) else NA_real_
  V1row("R2", "tau_bar lands on the calibration, and how much of it reaches the boat",
        sprintf("tau_bar %s [%s, %s] (R1 %s, %+.1f%%); boat all-gear %s vs R1 %s (%+.1f%%); pass-through %s of the prior move",
                fmt(.num1(tb$mean), 3), fmt(.num1(tb[["2.5%"]]), 3), fmt(.num1(tb[["97.5%"]]), 3),
                fmt(.tb0, 3), .pct(.num1(tb$mean), .tb0), fmt(ba, 0), fmt(ba0, 0), .pct(ba, ba0),
                if (is.na(.tpass)) "NA" else sprintf("%.2f", .tpass)),
        "tau_bar 2.8-3.1; boat all-gear +2 to +20%; pass-through well below 1",
        if (isTRUE(.num1(tb$mean) > 2.75) && isTRUE(.num1(tb$mean) < 3.2) && isTRUE(.pct(ba, ba0) > 2) && isTRUE(.pct(ba, ba0) < 20))
          "PASS" else "REVIEW",
        paste("tau_bar is identified by the OSP/trailer overlap and the prior should no longer fight it:",
              "the likelihood-only value was 2.97 against the shrunk 2.60. What the boat COMPONENT does",
              "with that is a separate question, and the answer measured on 2026-09-11 is that it absorbs",
              "about two thirds of the move in the latent effort level rather than passing it through.",
              "A pass-through near 1 would mean the trailer, OSP and catch streams had stopped",
              "constraining the effort scale, which is the thing to watch here."))
  pt <- rd(dir, "ppc_calibration_private_boat_all_gear_Dungeness_Kept.csv")
  if (!is.null(pt))
    V1row("R2", "boat PPC still calibrated after the prior move",
          paste(utils::capture.output(print(utils::head(pt, 6))), collapse = " | "),
          "PIT means near 0.5, coverage near nominal", "READ",
          "The prior move should improve, not worsen, the OSP stream's fit: read the OSP/trailer rows.")
}

# ---------------------------------------------------------------------------
# R2f (new_throughout only): THE FACTORIZATION PROOF, at one fit instead of two.
# R2f is R2 with the f block rolled back, so the pair differs ONLY in f. The design
# claim is that f enters the boat GENERATED QUANTITIES and nothing else, so:
#   - the shore fits must be BIT-IDENTICAL (the shore has no f at all), and
#   - every boat parameter that is not an f term must agree within Monte Carlo error.
# Bit-identity is impossible on the boat: the parameter vector gained z_f, sigma_f and
# cfi_kappa, so every HMC trajectory differs. Agreement in DISTRIBUTION is what the design
# guarantees, and fit_agreement() is that test.
# ---------------------------------------------------------------------------
verdict_R2f <- function(dir, ref) {
  if (is.na(dir %||% NA) || is.na(ref %||% NA)) return(invisible(NULL))
  fe <- fit_exactness(dir, ref, pat = "shore", what = "shore fits, R2f vs R2",
                      expect_delta = c("crab_fraction_strata", "crab_fraction_source", "crab_fraction_dynamic",
                                       "crab_fraction_rows", "run_tag", "model",
                                       "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1cross("R2f", "the shore did not move when f was rolled back (f is boat-only)", fe$observed,
        "shore fits bit-identical to R2", fe$verdict,
        "crab_fraction_stan_data() sets apply_crab_fraction = 0 for shore; a FAIL means an f term reached the shore.",
        dir, ref)
  fa <- fit_agreement(dir, ref, pat = "private_boat",
                      exclude = F_EXCLUDE,
                      what = "boat non-f parameters, R2f vs R2")
  V1cross("R2f", "the f block leaves effort and CPUE untouched (factorization)", fa$observed,
        "max |z| under 5 and under 1% of rows above 3", fa$verdict,
        paste("This is the R3-vs-R3a test of the 2026-09-08 design, run as a single control rung",
              "against the CURRENT turnover prior instead of the retired one. A FAIL means an f term",
              "reached an effort or CPUE likelihood, which would make every boat number in the ladder",
              "attributable to two things at once."),
        dir, ref)
  ba <- .comp(dir, "private_boat (All gear)"); ba1 <- .comp(ref, "private_boat (All gear)")
  V1row("R2f", "what the new f is worth on the boat, at a fixed turnover",
        sprintf("boat all-gear with the new f %s vs the retired f %s (%+.1f%%); shore unchanged",
                fmt(ba1, 0), fmt(ba, 0), .pct(ba1, ba)),
        "informational: this is the f effect the 'ladder' mode spends two fits to attribute", "INFO",
        paste("The boat is linear in f day by day, so this ratio is the effort-weighted change in f.",
              "Read it against crab_fraction_strata_*.csv from the R2 folder, which carries the",
              "monthly f the walk produced and the contact counts behind each month."))
}

verdict_R3a <- function(dir, prev) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  fe <- fit_exactness(dir, prev %||% "", pat = "shore", what = "shore fits vs the previous rung",
                      expect_delta = c("crab_fraction_strata", "crab_fraction_source", "crab_fraction_rows", "run_tag", "model",
                                       "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1cross("R3a", "the shore did not move (f is boat-only)", fe$observed, "shore fits bit-identical", fe$verdict,
        "f enters the boat generated quantities only; the shore has no f.", dir, prev %||% "")
  fs <- rd(dir, "crab_fraction_strata_private_boat_all_gear_Dungeness_Kept.csv")
  if (!is.null(fs)) {
    sh <- fs$contacts_crabbing / pmax(fs$contacts, 1); inf <- fs$contacts >= 20
    V1row("R3a", "monthly f under the legacy construction: informed months track the contacts, thin months sit at 0.30",
          paste(sprintf("%s: share %.2f (n=%d) -> f %.2f [%.2f, %.2f]", fs$label, sh, fs$contacts, fs$f_median, fs$f_lo95, fs$f_hi95), collapse = "; "),
          "months with n >= 20 pulled toward their share (shrunk by the Beta(6,14) prior); others 0.30",
          if (all(abs(fs$f_median[!inf] - 0.30) < 0.05) && all(abs(fs$f_median[inf] - sh[inf]) < abs(0.30 - sh[inf]) + 1e-9)) "PASS" else "REVIEW",
          paste("This is Phase A: the existing machinery fed the contacts. Its kappa = 20 prior is worth",
                "twenty boats at 0.30, so a 46-of-47 December comes out near 0.78 and the four thin months",
                "(Nov, Feb, Mar, Jun) fall back to the placeholder. R3 exists to remove both defects."))
  }
  ba <- .comp(dir, "private_boat (All gear)"); ba0 <- .comp(prev %||% "", "private_boat (All gear)")
  V1row("R3a", "the boat under monthly f", sprintf("boat all-gear %s vs previous rung %s (%+.1f%%)", fmt(ba, 0), fmt(ba0, 0), .pct(ba, ba0)),
        "up: the winter months carry most of the boat catch and their f rises from 0.30", "INFO",
        "Read together with R3, which moves the same months further because it does not shrink them to 0.30.")
}

verdict_R3 <- function(dir, prev) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  fe <- fit_exactness(dir, prev %||% "", pat = "shore", what = "shore fits vs the previous rung",
                      expect_delta = c("crab_fraction_dynamic", "crab_fraction_strata", "crab_fraction_source", "crab_fraction_rows",
                                       "run_tag", "model", "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1cross("R3", "the shore did not move (f is boat-only)", fe$observed, "shore fits bit-identical", fe$verdict,
        "f enters the boat generated quantities only; the shore has no f.", dir, prev %||% "")
  # THE FACTORIZATION PROOF, numerically: every non-f boat parameter within MC error.
  fa <- fit_agreement(dir, prev %||% "", pat = "private_boat",
                      exclude = F_EXCLUDE,
                      what = "boat non-f parameters vs the previous rung")
  V1cross("R3", "the dynamic f leaves effort and CPUE untouched (factorization)", fa$observed,
        "max |z| under 5 and under 1% of rows above 3", fa$verdict,
        paste("Bit-identity is impossible here: the parameter vector changed (z_f, sigma_f, cfi_kappa),",
              "so every HMC trajectory differs. What the design guarantees is agreement in DISTRIBUTION,",
              "and this is that test. A FAIL means an f term reached an effort or CPUE likelihood."),
        dir, prev %||% "")
  fs <- rd(dir, "crab_fraction_strata_private_boat_all_gear_Dungeness_Kept.csv")
  fa_prev <- rd(prev %||% "", "crab_fraction_strata_private_boat_all_gear_Dungeness_Kept.csv")
  if (!is.null(fs)) {
    sh <- fs$contacts_crabbing / pmax(fs$contacts, 1); thin <- fs$contacts < 20
    V1row("R3", "monthly f under the walk: tracks the contacts, thin months borrow from their neighbours",
          paste(sprintf("%s: share %.2f (n=%d) -> f %.2f [%.2f, %.2f]", fs$label, sh, fs$contacts, fs$f_median, fs$f_lo95, fs$f_hi95), collapse = "; "),
          "informed months within ~0.05 of their share; thin months between their neighbours, wider intervals",
          if (all(abs(fs$f_median[!thin] - sh[!thin]) < 0.08) &&
                all((fs$f_hi95 - fs$f_lo95)[thin] >= (fs$f_hi95 - fs$f_lo95)[!thin][1] * 0.8)) "PASS" else "REVIEW",
          paste("Compare with R3a month by month: the same contacts, no 0.30 anchor. A thin month's f",
                "should sit between its neighbours with a visibly wider interval; an informed month's",
                "should sit on its share. The December shrinkage R3a showed should be gone."))
  }
  sf <- .full_row(dir, "private_boat_all_gear", "sigma_f_out")
  V1row("R3", "the walk's step SD is identified", sprintf("sigma_f %s [%s, %s]", fmt(.num1(sf$mean), 3), fmt(.num1(sf[["2.5%"]]), 3), fmt(.num1(sf[["97.5%"]]), 3)),
        "posterior away from zero and inside the half-normal(1.5) prior", if (isTRUE(.num1(sf[["2.5%"]]) > 0.1) && isTRUE(.num1(sf$mean) < 1.5)) "PASS" else "REVIEW",
        "The 2024-25 desk fit gave 0.66 [0.34, 1.38]; a posterior pinned near zero would mean the walk is flat.")
  ba <- .comp(dir, "private_boat (All gear)"); ba0 <- .comp(prev %||% "", "private_boat (All gear)")
  fw <- if (!is.null(fs) && !is.null(fa_prev) && nrow(fs) == nrow(fa_prev))
    sprintf("; catch-weighted f moved %.3f -> %.3f", stats::weighted.mean(fa_prev$f_median, fa_prev$contacts + 1), stats::weighted.mean(fs$f_median, fs$contacts + 1)) else ""
  V1row("R3", "the boat moves with the catch-weighted change in f", sprintf("boat all-gear %s vs previous rung %s (%+.1f%%)%s", fmt(ba, 0), fmt(ba0, 0), .pct(ba, ba0), fw),
        "roughly the change in the effort-weighted f (+15 to +25% against R3a, more against R2)", "READ",
        paste("The boat is linear in f day by day, so the component moves by the effort-weighted change",
              "in f, which the strata file lets you compute by hand. Anything not explained by f is a",
              "failure of the factorization row above."))
  V1row("R3", "the port total at the shipped configuration",
        sprintf("%s [%s, %s] (baseline %s [%s, %s])", fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_median), 0),
                fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_lo95), 0), fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_hi95), 0),
                fmt(REF$A1$port, 0), fmt(REF$A1$lo95, 0), fmt(REF$A1$hi95, 0)),
        "informational", "INFO",
        "The candidate for adoption if R4 is not adopted: items 8, 3, 1A, 1B in one number.")
}

verdict_R4 <- function(dir, prev, prev_id = "the previous rung") {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  # 2026-09-11, from the first full run: the three census keys were MISSING from
  # expect_delta, so this row reported "3 UNEXPECTED: census_uncertainty,
  # charter_expansion, charter_frame" on a PASS. They belong to D_R4 by design (the
  # 2026-09-11 census split ships with the derived shore turnover), so they are expected.
  # An UNEXPECTED flag that fires on a correct configuration is worse than no flag: it
  # trains the reader to ignore the one place the runner says "this comparison is not what
  # you think it is".
  fe <- fit_exactness(dir, prev %||% "", pat = "private_boat",
                      what = sprintf("boat fits vs %s", prev_id),
                      expect_delta = c("tau_shore_prior_mu", "tau_shore_prior_sigma", "tau_shore_prior_source", "shared_tau_min_obs",
                                       "census_uncertainty", "charter_frame", "charter_expansion",
                                       "run_tag", "model", "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1cross("R4", "the boat did not move (the shore prior cannot reach it)", fe$observed,
        sprintf("boat fits bit-identical to %s", prev_id), fe$verdict,
        "A FAIL means the shore turnover leaked into the boat fits, which it has no path to do.",
        dir, prev %||% "")
  sa <- .comp(dir, "shore (All gear)"); sa0 <- .comp(prev %||% "", "shore (All gear)")
  sp <- .comp(dir, "shore (Pot closure)"); sp0 <- .comp(prev %||% "", "shore (Pot closure)")
  ea <- .comp(dir, "shore (All gear)", "BSS_effort"); ea0 <- .comp(prev %||% "", "shore (All gear)", "BSS_effort")
  V1row("R4", "the shore scales by the ratio of the two turnovers",
        sprintf("shore all-gear effort %s vs %s (x%.3f); catch %s vs %s (x%.3f); pot closure catch %s vs %s (x%.3f)",
                fmt(ea, 0), fmt(ea0, 0), ea / ea0, fmt(sa, 0), fmt(sa0, 0), sa / sa0, fmt(sp, 0), fmt(sp0, 0), sp / sp0),
        "about x1.45 (2.48 / 1.7), effort and catch alike, both sub-seasons",
        if (isTRUE(abs(ea / ea0 - 2.48 / 1.7) < 0.12)) "PASS" else "REVIEW",
        paste("The shore is linear in tau_shore; the derived prior is 2.48 against 1.7 and the shared level",
              "carries it. A ratio well off 1.45 means the data (I/E trips) pulled the level, which is",
              "worth reading: the I/E stream is the only shore observation that can. Measured 2026-09-11:",
              "x1.356 on the all-gear effort, i.e. the data pulled the posterior turnover slightly BELOW",
              "the prior centre (2.394 against 2.477) and the tight derived prior (log SD 0.1) held the",
              "rest. Read this row together with the sigma_IE row below: they are two views of the same",
              "disagreement."))
  s1 <- .full_row(dir, "shore_all_gear", "sigma_IE_out"); s0 <- .full_row(prev %||% "", "shore_all_gear", "sigma_IE_out")
  # HOW MANY OBSERVATIONS IS THIS VERDICT RESTING ON? Added 2026-09-11, because the first
  # full run returned sigma_IE 0.577 against 0.373 (+55%, a REVIEW) and the number that
  # makes that interpretable was nowhere on the page: the shore ALL-GEAR fit has only FOUR
  # in-window I/E days (the shore turnover itself is derived from 40 days, 36 of them
  # outside this sub-season). A +55% move in a parameter estimated from four observations
  # is a direction, not a measurement, and the row has to say so.
  .ie_n <- tryCatch({
    x <- rd(dir, "L_effective_ie_detail.csv")
    pc <- rd(dir, "bss_period_coverage_shore_all_gear_Dungeness_Kept.csv")
    if (is.null(x) || is.null(pc)) NA_integer_ else {
      d1 <- min(as.Date(pc$date_start)); d2 <- max(as.Date(pc$date_end))
      sum(as.Date(x$event_date) >= d1 & as.Date(x$event_date) <= d2, na.rm = TRUE)
    }
  }, error = function(e) NA_integer_)
  .ie_tot <- tryCatch(nrow(rd(dir, "L_effective_ie_detail.csv")), error = function(e) NA_integer_)
  V1row("R4", "the I/E observation scale did not worsen",
        sprintf("sigma_IE %s (%s %s); the shore all-gear likelihood has %s in-window I/E day(s) of %s in the workbook",
                fmt(.num1(s1$mean), 3), prev_id, fmt(.num1(s0$mean), 3),
                if (is.na(.ie_n)) "?" else as.character(.ie_n),
                if (is.na(.ie_tot)) "?" else as.character(.ie_tot)),
        sprintf("not more than 5%% above %s's", prev_id),
        if (isTRUE(.num1(s1$mean) <= .num1(s0$mean) * 1.05)) "PASS" else "REVIEW",
        paste("The I/E stream is the shore's only direct observation of the turnover: the likelihood is",
              "arrivals ~ lognormal(log(lambda_E * tau_shore), sigma_IE), so if the derived tau over-expands",
              "the count, sigma_IE grows to absorb it. Measured 2026-09-11: 0.373 -> 0.577, and on the four",
              "in-window days the predicted arrivals went from 1.25x the observed to 1.64x. The direction is",
              "consistent and physically interpretable; the magnitude is not well identified (the two",
              "posteriors overlap heavily, [0.03, 0.84] against [0.29, 1.03]) because n = 4. The pot-closure",
              "component, which has its own I/E days, did NOT move (0.206 -> 0.202). See CHANGE_REGISTER D24."))
  V1row("R4", "the port total with the derived shore turnover",
        sprintf("%s [%s, %s] (%s %s; baseline %s)", fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_median), 0),
                fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_lo95), 0), fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_hi95), 0),
                prev_id, fmt(.num1(.port_row(rd(prev %||% "", "port_total_Dungeness_Kept.csv"))$BSS_median), 0), fmt(REF$A1$port, 0)),
        "informational; the decision (adopt, or hold for a field check) is Matt's", "INFO",
        "A paired gear count on every I/E day is the cheap field confirmation of the derived turnover.")
}

verdict_R5 <- function(dir, pooled_dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  gp <- .num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_median)
  pp <- .num1(.port_row(rd(pooled_dir %||% "", "port_total_Dungeness_Kept.csv"))$BSS_median)
  V1row("R5", sprintf("the two tracks agree on the %s configuration", GEAR_FOLLOWS),
        sprintf("gear port %s vs pooled %s (%+.2f%%); boat all-gear gear %s vs pooled %s; shore all-gear gear %s vs pooled %s",
                fmt(gp, 0), fmt(pp, 0), .pct(gp, pp), fmt(.comp(dir, "private_boat (All gear)"), 0), fmt(.comp(pooled_dir, "private_boat (All gear)"), 0),
                fmt(.comp(dir, "shore (All gear)"), 0), fmt(.comp(pooled_dir, "shore (All gear)"), 0)),
        "within about 2% at the port", if (isTRUE(abs(.pct(gp, pp)) <= 2)) "PASS" else "REVIEW",
        paste("Two independently parameterized CPUE structures fitting the same data. The gear track's shore",
              "all-gear fit sits at monthly and its catch likelihood is plain NB2, so a small standing gap",
              "is expected (about -0.8% on 2026-09-05); a large one after these patches points at the boat."))
  tg <- .full_row(dir, "private_boat_all_gear", "tau_bar_out"); tp <- .full_row(pooled_dir %||% "", "private_boat_all_gear", "tau_bar_out")
  V1row("R5", "the shared turnover agrees across tracks", sprintf("gear tau_bar %s vs pooled %s (%+.2f%%)", fmt(.num1(tg$mean), 4), fmt(.num1(tp$mean), 4), .pct(.num1(tg$mean), .num1(tp$mean))),
        "within a few percent", if (isTRUE(abs(.pct(.num1(tg$mean), .num1(tp$mean))) < 5)) "PASS" else "REVIEW",
        "tau_bar is a property of the OSP/trailer data; the two tracks estimated it to 0.03% on 2026-09-01.")
  fg <- rd(dir, "crab_fraction_strata_private_boat_all_gear_Dungeness_Kept.csv"); fp <- rd(pooled_dir %||% "", "crab_fraction_strata_private_boat_all_gear_Dungeness_Kept.csv")
  if (!is.null(fg) && !is.null(fp) && nrow(fg) == nrow(fp))
    V1row("R5", "the dynamic f agrees across tracks (it reads the same contacts)",
          sprintf("max |f_gear - f_pooled| over months %.3f", max(abs(fg$f_median - fp$f_median))), "under 0.03", if (max(abs(fg$f_median - fp$f_median)) < 0.03) "PASS" else "REVIEW",
          "The f block is identical in both Stan files and reads the same rows; a gap here is a data-plumbing bug, not a model difference.")
}

# ---------------------------------------------------------------------------
# DRIVE
# ---------------------------------------------------------------------------
banner(sprintf("IMPROVEMENT LADDER  2026-09-08   DRY_RUN=%s  F_METHOD=%s  stages: %s  (R5 follows %s)",
               DRY_RUN, F_METHOD, paste(STAGES, collapse = ", "), GEAR_FOLLOWS))
preflight()

# A FAILED RUNG MUST NOT DESTROY THE OTHERS (2026-09-13). This was a bare
#   for (sid in STAGES) dirs[[sid]] <- run_stage(sid)
# so an error inside any render -- a Stan failure, a full disk, one bad interview row --
# propagated out of the loop and killed the script before a single verdict was written.
# On a 12-to-14 h overnight run that means a failure in the last rung throws away the
# hours the earlier ones cost. The rungs are independent RENDERS (R4 does not read R2's
# output; only the VERDICTS compare folders), so the right behaviour is to record the
# failure, keep going, and let every verdict that can still be computed be computed.
dirs <- list(); .failed <- character(0); .t_all <- Sys.time()
for (sid in STAGES) {
  dirs[[sid]] <- tryCatch(run_stage(sid), error = function(e) {
    .failed <<- c(.failed, sid)
    cat(sprintf("\n  *** RUNG %s FAILED: %s\n", sid, conditionMessage(e)))
    cat("      Continuing with the remaining rungs. Every rung that finished keeps its\n")
    cat("      output and its verdicts; re-run with RESUME = TRUE to retry this one.\n")
    V1row(sid, "the rung did not complete", conditionMessage(e), "the render runs to completion", "ERROR",
          paste("The render errored. Nothing about the other rungs is affected: they are separate",
                "renders and only the verdict blocks read across folders. Re-running the script with",
                "RESUME = TRUE will reuse every completed rung by its config digest and retry only",
                "this one. Look in the rung's stan_console_*.log first if it is a sampler failure."))
    NA_character_
  })
  if (!isTRUE(DRY_RUN) && !identical(STAGE_DEFS[[sid]]$model, "desk"))
    cat(sprintf("  [elapsed %.1f h of the ladder so far]\n",
                as.numeric(difftime(Sys.time(), .t_all, units = "hours"))))
}
if (length(.failed))
  cat(sprintf("\n  %d rung(s) failed: %s. The verdicts below cover the rest.\n",
              length(.failed), paste(.failed, collapse = ", ")))

.safe <- function(sid, expr) tryCatch(force(expr), error = function(e) {
  cat(sprintf("\n  *** VERDICT BLOCK %s FAILED: %s\n", sid, conditionMessage(e)))
  cat("      The run itself is intact; re-run with RESUME = TRUE after fixing the block.\n")
  V1row(sid, "verdict block did not complete", conditionMessage(e), "the block runs to completion", "ERROR",
        "The fits are on disk and unaffected: this is a defect in the code that READS them; re-run with RESUME = TRUE.")
  invisible(NULL)
})
.dir_of <- function(sid) if (!is.na(sid %||% NA) && sid %in% names(dirs)) dirs[[sid]] else NA_character_
for (sid in intersect(c("R1", "R2", "R2f", "R3a", "R3", "R4"), STAGES)) .safe(sid, ladder_row(sid, dirs[[sid]]))
if ("R1"  %in% STAGES) .safe("R1",  verdict_R1(dirs$R1))
if ("R2"  %in% STAGES) .safe("R2",  verdict_R2(dirs$R2, .dir_of(prev_pooled("R2"))))
if ("R2f" %in% STAGES) .safe("R2f", verdict_R2f(dirs$R2f, .dir_of("R2")))
if ("R3a" %in% STAGES) .safe("R3a", verdict_R3a(dirs$R3a, .dir_of(prev_pooled("R3a"))))
if ("R3"  %in% STAGES) .safe("R3",  verdict_R3(dirs$R3, .dir_of(prev_pooled("R3"))))
if ("R4"  %in% STAGES) .safe("R4",  verdict_R4(dirs$R4, .dir_of(prev_pooled("R4")), prev_pooled("R4")))
if ("R5"  %in% STAGES) .safe("R5",  verdict_R5(dirs$R5, .dir_of(GEAR_FOLLOWS)))

if (length(LAD)) {
  lp <- .here("05_output", sprintf("improvements_2026-09-08_ladder%s.csv", .sfx))
  merge_csv_by(do.call(rbind, LAD), lp, "rung")
  banner("LADDER")
  print(do.call(rbind, LAD), row.names = FALSE)
  cat("\n  written to", lp, "\n")
}
# ---------------------------------------------------------------------------
# WHAT TO DO NEXT. Printed at the end of every non-dry run so the two-pass plan does not
# have to be remembered, and so a pass that ended early says so.
# ---------------------------------------------------------------------------
if (!isTRUE(DRY_RUN)) {
  .fitted <- setdiff(STAGES, "R0")
  .have <- vapply(.fitted, function(sid) !is.na(.dir_of(sid) %||% NA), logical(1))
  banner("NEXT")
  if (!all(.have)) {
    cat(sprintf("  %d of %d fitted rung(s) produced no output folder: %s%s\n", sum(!.have), length(.fitted),
                paste(.fitted[!.have], collapse = ", "),
                if (length(.failed)) sprintf("  (errored: %s)", paste(.failed, collapse = ", ")) else ""))
    cat("  Re-run this script as-is: RESUME = TRUE reuses everything that finished and re-fits the rest.\n")
  } else if (identical(F_METHOD, "ladder")) {
    cat("  The full 'ladder' set is on disk. Read the verdicts file, then the ladder table.\n")
  } else if (LADDER_PASS < 2) {
    cat(paste0("  PASS 1 COMPLETE: R1, R2, R4 and R5 are on disk and every rung's port total is citable.\n",
               "  To add the factorization control (R2f: R2 with the f block rolled back, ~4 h), set\n",
               "        LADDER_PASS <- 2\n",
               "  in the control block and re-run this script. RESUME will match the four fits above by\n",
               "  their config digest, skip them, fit R2f only, and recompute every verdict from disk.\n",
               "  Do not apply an unrelated patch in between: the digest covers the config, not the code,\n",
               "  and R2f is only interpretable against an R2 fitted by the same code (the run will say so\n",
               "  if they differ, and downgrade the affected verdicts).\n"))
  } else {
    cat(paste0("  PASS 2 COMPLETE: R2f is on disk. The rows to read are the two R2f verdicts -- the shore\n",
               "  must be bit-identical to R2 (f is boat-only) and every boat parameter that is not an f\n",
               "  term must agree within Monte Carlo error (the factorization proof). If either says\n",
               "  REVIEW for a COMPARABILITY reason, the fits came from different code or a different\n",
               "  rstan, not from a modelling failure.\n"))
  }
  cat(sprintf("  D19 is decidable now: compare the four PE arms in the R0 desk folder\n%s\n",
              "  (pe_unsampled_cell_arms.csv) against R4's BSS in pe_vs_bss_comparison.csv."))
}

if (length(V)) {
  vp <- .here("05_output", sprintf("improvements_2026-09-08_verdicts%s.csv", .sfx))
  merge_csv_by(do.call(rbind, V), vp, c("stage", "criterion"))
  banner("VERDICTS")
  for (r in seq_len(length(V))) with(V[[r]], cat(sprintf("  [%s] %-14s %s\n     obs: %s\n",
    stage, verdict, criterion, observed)))
  cat("\n  written to", vp, "\n")
}
