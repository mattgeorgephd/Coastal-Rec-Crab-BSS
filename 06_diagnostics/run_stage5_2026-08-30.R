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
###############################################################################
# run_stage5_2026-08-30.R
#
# Runs Stage 5 of development_notes/improvement-plan-2026-08-27.md end to end,
# unattended. Start it, walk away, come back to
#   05_output/stage5_2026-08-30_summary.csv
#   05_output/stage5_2026-08-30_verdicts.csv
#   05_output/stage5_2026-08-30_tau_by_ar.csv     <- the 2x2 the stage exists to fill in
#
# THE PROBLEM STAGE 5 EXISTS TO SOLVE, in one line: six configurations, all passing every
# convergence-gate criterion, span 25,883 to 37,392 on the private-boat all-gear component
# (44%) and 66,237 to 78,250 at the port (18%). The gate tests whether a fit SAMPLED. It was
# never meant to test whether the model is the right one, and the project has been reading
# it as though it were.
#
# THE ORDERING IS THE ARGUMENT. Read it once before starting.
#
#   S1  PE-fill comparison, minutes, no MCMC. Plan 5.1. This is a re-run: the 2026-08-28
#       stage A did not load the crabbing-holiday workbook, so it stratified 289 days into
#       85 week x day-type cells where the driver builds 92, and its zeroed-day counts were
#       41/44/8/0 against the driver's 44/47/9/0. The comparison was still fair (both arms
#       shared the stratification) but the absolute numbers were not production numbers.
#       The fix is in; this re-run must now reproduce the DRIVER's PE exactly, which is a
#       sharper criterion than the original stage A ever had.
#
#   S2  shared_tau ON, pooled, with the new informed-day floor. Plan 5.2. THE GATE.
#       In the 2026-08-29 batch the toggle was global and moved the SHORE all-gear component
#       +17.9% on the strength of 4 informed days out of 289, with an interval containing
#       the prior centre and no replication in the gear track. The floor
#       (shared_tau_min_obs, default 15) should now refuse the shore and allow the boat.
#       This run tests that as an EXACT statement, not a plausible one: see the gate below.
#
#   S3  shared_tau ON, gear track, boat all-gear fit given more draws AND tighter adaptation.
#       Plan 5.2 item 3. This DEPARTS from the plan's prescription; the departure is argued
#       at S3_SAMPLER in the control block and is the one place this script overrides the
#       written plan. Read that comment before running.
#
#   S4a shared_tau OFF, boat all-gear forced to a daily AR.   } Plan 5.3, the 2x2, and the
#   S4b shared_tau ON,  boat all-gear forced to a daily AR.   } run that settles Stage 5.
#       The plan calls this "one new run" because stage F already supplies an (off, daily)
#       cell. IT DOES NOT, and running only S4b would confound the answer: stage F reached
#       daily through ar_escalate, which moved ALL FOUR fits to daily (shore pot closure
#       biweekly->daily, shore all-gear already daily, boat pot closure monthly->daily, boat
#       all-gear monthly->daily) through a different code path with a different number of
#       fits per component. Two of its four components therefore differ from stage D for
#       reasons that have nothing to do with the boat's AR. S4a costs one run and buys a
#       cell that differs from S4b in exactly one toggle. That is the whole point of a 2x2.
#
#   S5  Boat pot-closure AR alignment, re-run. Plan 5.6. LAST because it is the least
#       decisive: the 2026-08-28 stage C answered the question (biweekly reconciles the two
#       tracks to 1.1%) but its ar_force was honoured per POPULATION, so it silently forced
#       the boat's ALL-GEAR sub-season to biweekly too and moved that component
#       25,883 -> 28,893 as a side effect. The nested form is now honoured per sub-season.
#       This re-run produces an uncontaminated port total for the change.
#
# WHAT IS NOT HERE, and why:
#   * Plan 5.4 (what ar_max_resolution is for) is a decision that reads S4a/S4b's output.
#     It cannot be run; it can only be made, and this batch is what makes it makeable.
#   * Plan 5.5 (adequacy beside the gate) is code, applied. Every stage here writes
#     model_adequacy.csv, and the archived 2026-08-29 cells were reconstructed from their
#     committed CSVs (annotate_model_adequacy_run) so the 2x2 quotes the same four
#     quantities for all four cells rather than for the new ones only.
#   * The remaining 5.6 items (tau_bar in prior_vs_posterior, expansion_ratios decoupled,
#     n_interviews_fitted on the gear track, the ie_obs_unit string) are code, applied, and
#     are checked by 06_diagnostics/test_improvements_2026-08-25.R rather than by a fit.
#
# HOW TO RUN
#     Rscript 06_diagnostics/run_stage5_2026-08-30.R      # from the repo root
#     source("06_diagnostics/run_stage5_2026-08-30.R")    # RStudio: Source
#
# START WITH DRY_RUN <- TRUE. It resolves and prints every stage's config, runs the
# pre-flight R-to-Stan data-contract check, runs the extractor self-test against the
# archived 2026-08-29 runs, and fits nothing. It needs no Stan toolchain. Then set
# DRY_RUN <- FALSE.
#
# RESUMABLE. A stage whose output folder already holds port_total_Dungeness_Kept.csv is
# skipped and re-extracted. Rows are appended to the summary CSV as each stage finishes, so
# a batch killed at hour 15 still leaves everything it learned.
#
# RUNTIME. Roughly 18-24 h on 4 cores. S1 is minutes; S4a and S4b are the long ones (a daily
# boat AR is 289 latent periods where monthly is 10).
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                  # TRUE: resolve, print, self-test, fit nothing. START HERE.
STAGES  <- c("S1", "S2", "S3", "S4a", "S4b", "S5")
RESUME  <- TRUE                  # skip a stage whose output folder already looks complete

# ---------------------------------------------------------------------------
# S3_SAMPLER -- READ THIS BEFORE RUNNING. It is a deliberate departure from the plan.
#
# Plan item 5.2 says: "Re-run the gear track at 2,500 iterations (it runs boat all-gear at
# 1,000, which is where its chain 3 stalled at the phi_E unit-root boundary)". Taken
# literally that prescription is very unlikely to work, and it is worth saying why rather
# than running it and reporting a null result.
#
# WHAT THE 2026-08-29 GEAR RUN ACTUALLY RECORDED for private_boat_all_gear:
#     tau_bar      n_eff   49   Rhat 1.081        554 divergences
#     r_OSP        n_eff   40   Rhat 1.077        gate: REJECTED, component fell back to PE
#     f_crab[1]    n_eff   25   Rhat 1.174
# n_eff 49 out of 4 x 1,000 = 4,000 post-warmup draws is an efficiency of 1.2%. If the only
# problem were draw count, reaching the gate's floor of 400 would need roughly 32,000 draws,
# EIGHT times the plan's 2,500-per-chain, not two and a half. And R-hat 1.08-1.17 does not
# fall with more draws when one chain is in a different place; more draws buys precision on
# a biased answer.
#
# WHAT THE COMPARISON WITH THE POOLED TRACK POINTS AT INSTEAD. The two tracks fit this
# component on the same data with the same latent structure, and the 2026-08-29 stage E runs
# differ like this:
#     pooled  boat all-gear:  adapt_delta 0.99, max_treedepth 13, 2,500 draws/chain ->  39 divergences, gate PASSED
#     gear    boat all-gear:  adapt_delta 0.90, max_treedepth 10, 1,000 draws/chain -> 554 divergences, gate REJECTED
# The pooled driver has an explicit per-fit override for this component
# (bss_iter_boat_allgear / bss_delta_boat_allgear = 0.99 / bss_treedepth_boat_allgear = 13).
# The gear driver has overrides for shore-all-gear and for both pot-closure fits, and NONE
# for boat all-gear, so that fit alone falls through to the track defaults. A 14x difference
# in divergences at 0.90 versus 0.99 adapt_delta is the ordinary signature of loose
# adaptation on a difficult posterior, not of too few draws.
#
# WHAT THIS SCRIPT DOES BY DEFAULT: raises draws AND adaptation together, to the pooled
# track's settings for the same component.
#
# THE COST, STATED: this changes three things at once, so if the gear fit converges you
# cannot say which of them did it. That is accepted deliberately, because the OBJECTIVE here
# is not attribution, it is obtaining a converged cross-track fit to compare against the
# pooled boat result. Attribution between draws and adaptation would need its own 2x2 and
# would answer a question about the sampler, not about the crab.
#
# IF YOU WANT THE PLAN AS WRITTEN, use the commented line instead and expect a rejected fit.
#
# WHY THIS NEEDS bss_sampler_override AT ALL: each driver merges its params_model ON TOP of
# run_config, so params_model wins every sampler key. Setting bss_iter_default in a stage
# delta does NOTHING; the driver overwrites it and the run looks like it complied. The
# override list is applied after the merge, accepts sampler keys only, prints every change,
# and errors on anything else (03_R_functions/bss_sampler_override.R). This is the same
# class of trap that left bss_min_interviews pinned at 20 until the 2026-08-25 audit.
#
# TARGETING, CHECKED: in the gear driver the per-fit branch is
#   shore & all_gear -> bss_*_shore_allgear ; pot_closure -> bss_*_ringnet ; else -> defaults
# and "else" is private_boat all-gear and nothing else. Raising the DEFAULTS therefore
# touches exactly the fit that failed, and leaves the other three fits bit-identical to the
# 2026-08-29 gear run. That is checked as a criterion, not assumed.
# ---------------------------------------------------------------------------
S3_SAMPLER <- list(bss_iter_default          = 5000,   # 2,500 post-warmup draws per chain
                   bss_warmup_default        = 2500,
                   bss_adapt_delta_default   = 0.99,   # pooled uses 0.99 for this component
                   bss_max_treedepth_default = 13)     # pooled uses 13
# S3_SAMPLER <- list(bss_iter_default = 5000, bss_warmup_default = 2500)   # the plan as written

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
fmt    <- function(x, d = 1) if (length(x) == 0 || is.na(x)) "NA" else
  formatC(as.numeric(x), format = "f", digits = d, big.mark = ",")

if (!isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages({ library(here); library(rmarkdown) })
  load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
  install.lib <- load.lib[!load.lib %in% installed.packages()]
  for (lib in install.lib) install.packages(lib, dependencies = TRUE)
  invisible(sapply(load.lib, require, character.only = TRUE))
  rstan_options(auto_write = TRUE)
  invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), source))
} else {
  # The retro annotator and the adequacy reader are pure R over CSVs; the dry run's
  # self-test uses them, so source just those two rather than the whole library (which
  # would pull rstan in).
  for (f in c("model_diagnostics.R", "bss_model_adequacy.R", "annotate_decoupled_run.R",
              "bss_sampler_override.R", "bss_ar_resolution.R"))
    try(sys.source(.here("03_R_functions", f), envir = globalenv()), silent = TRUE)
}

source(.here("run_config.R"))          # run_config = the SHIPPED config = ladder rung 4
BASE <- run_config

model_rmd <- list(
  pooled        = .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd"))
prefix <- list(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")
stopifnot(all(file.exists(unlist(model_rmd))))

# ---------------------------------------------------------------------------
# PRE-FLIGHT: the R-to-Stan data contract. A text parse, no packages, milliseconds. This is
# the check that would have saved the 2026-08-25 ladder's six lost runs, where five
# unforwarded Stan data fields made rstan return an EMPTY stanfit rather than raise, and the
# failure surfaced hours later as "dim(X) must have a positive length".
# ---------------------------------------------------------------------------
local({
  src_fn <- .here("03_R_functions", "bss_stan_fit.R")
  if (!file.exists(src_fn)) { warning("bss_stan_fit.R not found; pre-flight skipped."); return(invisible(NULL)) }
  e <- new.env(); sys.source(src_fn, envir = e)
  pairs <- list(pooled        = c("crab_bss_pooled.stan",        "prep_bss_crab_pooled.R"),
                gear_resolved = c("crab_bss_gear_resolved.stan", "prep_bss_crab_gear.R"))
  fail <- character(0)
  for (m in names(pairs)) {
    need <- e$bss_stan_data_names(.here("02_stan_models", pairs[[m]][1]))
    src  <- paste(readLines(.here("03_R_functions", pairs[[m]][2]), warn = FALSE), collapse = "\n")
    miss <- need[!vapply(need, function(v)
      grepl(paste0("(^|[^A-Za-z0-9_.])", v, "([^A-Za-z0-9_]|$)"), src, perl = TRUE), logical(1))]
    cat(sprintf("  pre-flight  %-14s %2d Stan data variables, %d not built in %s\n",
                m, length(need), length(miss), basename(pairs[[m]][2])))
    if (length(miss)) fail <- c(fail, sprintf("%s: %s", m, paste(miss, collapse = ", ")))
  }
  if (length(fail))
    stop("PRE-FLIGHT FAILED, the R-to-Stan data contract is broken. Every fit would return ",
         "an empty stanfit.\n  ", paste(fail, collapse = "\n  "), call. = FALSE)
})

# PRE-FLIGHT 2: the sampler-override hatch S3 depends on. If bss_apply_sampler_override is
# absent or not wired into both drivers, S3 would run at the driver's defaults and report
# that more draws did not help. Fail now, not at hour 6.
local({
  fn <- .here("03_R_functions", "bss_sampler_override.R")
  if (!file.exists(fn)) stop("03_R_functions/bss_sampler_override.R is missing; S3 cannot ",
                             "raise the gear track's sampler settings.", call. = FALSE)
  hooked <- vapply(unlist(model_rmd), function(p)
    any(grepl("bss_apply_sampler_override", readLines(p, warn = FALSE), fixed = TRUE)),
    logical(1))
  if (!all(hooked))
    stop("bss_apply_sampler_override is not called in: ",
         paste(basename(names(hooked)[!hooked]), collapse = ", "),
         ". A sampler override would be silently ignored.", call. = FALSE)
  cat(sprintf("  pre-flight  sampler-override hatch present and wired into %d/%d drivers\n",
              sum(hooked), length(hooked)))
})

# ---------------------------------------------------------------------------
# BASELINES. Hard-coded on purpose: a criterion you can edit after seeing the answer is not
# a criterion. Every value below is read from a committed run folder, named beside it.
#
# ON PORT TOTALS AND BIT-IDENTITY. Stage D reproduced rung 4's fits EXACTLY (every shared
# parameter row identical at full double precision) and still reported a port total of
# 66,094 against 66,237, a 0.2% gap. That is not a discrepancy: rstan::extract(permuted =
# TRUE) permutes draws, so a port total assembled by resampling component draws is
# RNG-sensitive even when the fits are bit-identical. Consequence for this script: EVERY
# exactness gate below compares per-fit summaries, NEVER port totals.
# ---------------------------------------------------------------------------
REF <- list(
  # Ladder rung 4, the shipped configuration.
  rung4 = list(dir = "20260826/pooled-CPUE-PV4-minint",
               port = 66237, shore_ag = 20898, shore_pc = 6331, boat_ag = 25868, boat_pc = 849,
               # PE points, as the driver computes them with the holiday workbook loaded.
               pe_shore_ag = 15755, pe_shore_pc = 5959, pe_boat_ag = 11176, pe_boat_pc = 400,
               # Zeroed-day counts from the driver's own pe_empty_effort_strata.csv.
               empty_days = c(shore_all_gear = 44, shore_ring_net_only = 0,
                              private_boat_all_gear = 47, private_boat_ring_net_only = 9)),
  gear5 = list(dir = "20260826/gear-type-CPUE-model-PV5-gear-ship",
               port = 65756, boat_ag = 25796, boat_pc = 743),
  # 2026-08-29 batch. D = (tau off, boat monthly); E = (tau on, boat monthly);
  # F = ar_escalate, which put ALL FOUR fits at daily.
  D = list(dir = "20260828/pooled-CPUE-IP-D-tau-off",
           port = 66094, shore_ag = 20898, shore_pc = 6331, boat_ag = 25868, boat_pc = 849,
           boat_ag_ar = "monthly", r_OSP = 1.617, r_OSP_hi95 = 2.208,
           pit_trailer = 0.4238, pit_osp = 0.5790, div_boat_ag = 97,
           p_loo_frac = 0.080, pareto_bad = 1, disp_neff = 9574),
  E = list(dir = "20260829/pooled-CPUE-IP-E-tau-on-pooled",
           port = 75619, shore_ag = 24629, shore_pc = 6122, boat_ag = 31008, boat_pc = 1018,
           boat_ag_ar = "monthly", tau_bar = 2.5969, tau_lo95 = 2.0644, tau_hi95 = 3.2488,
           r_OSP = 1.540, r_OSP_hi95 = 2.030,
           pit_trailer = 0.4858, pit_osp = 0.5211, div_boat_ag = 39,
           p_loo_frac = 0.080, pareto_bad = 1, disp_neff = 11125),
  Egear = list(dir = "20260829/gear-type-CPUE-model-IP-E-tau-on-gear_resolved",
               boat_ag_bss = 32689, boat_ag_method = "PE", tau_bar = 2.5515,
               tau_neff = 49, tau_Rhat = 1.0813, div_boat_ag = 554, disp_neff = 40),
  F = list(dir = "20260829/pooled-CPUE-IP-F-escalate",
           port = 78250, shore_ag = 20898, shore_pc = 6352, boat_ag = 37359, boat_pc = 1106,
           boat_ag_ar = "daily", r_OSP = 10.607, r_OSP_hi95 = 1405.5,
           pit_trailer = 0.4039, pit_osp = 0.5582, div_boat_ag = 154,
           p_loo_frac = 0.163, pareto_bad = 2, disp_neff = 307),
  # 2026-08-28 stage C, CONTAMINATED: ar_force was honoured per population, so both boat
  # sub-seasons went biweekly. Kept only so S5 can show what the contamination was worth.
  Cbad = list(dir = "20260828/pooled-CPUE-IP-C-boatpc-ar",
              port = 68853, boat_ag = 28893, boat_pc = 735))

# ---------------------------------------------------------------------------
# THE STAGES. `delta` is applied on top of the SHIPPED config, never cumulatively, so every
# result is attributable to one named change against rung 4.
# ---------------------------------------------------------------------------
stage <- function(id, tag, model, delta, headline, plan_item, kind = "run")
  list(id = id, tag = tag, model = model, delta = delta, headline = headline,
       plan_item = plan_item, kind = kind)

# The nested ar_force form. As of 2026-08-29 .bss_resolve_ar_force() honours a list keyed by
# gear_regime, and ERRORS on an unnamed list rather than silently forcing every sub-season,
# which is what contaminated stage C. It bypasses both the selector and ar_max_resolution.
AR_BOAT_AG_DAILY <- list(private_boat = list(all_gear    = "daily"))
AR_BOAT_PC_BIWK  <- list(private_boat = list(pot_closure = "biweekly"))

STAGE_DEFS <- list(
  S1  = stage("S1", "S5-1-pe-fill", "pooled", list(),
              "empty-stratum fill re-run on the fixed holiday path, PE only, no fitting",
              "5.1", kind = "pe_only"),
  S2  = stage("S2", "S5-2-tau-pooled", "pooled", list(shared_tau = TRUE),
              "GATE: the informed-day floor must refuse the shore and allow the boat",
              "5.2"),
  S3  = stage("S3", "S5-3-tau-gear", "gear_resolved",
              list(shared_tau = TRUE, bss_sampler_override = S3_SAMPLER),
              "cross-track check: the gear boat all-gear fit, given the pooled track's sampler",
              "5.2"),
  S4a = stage("S4a", "S5-4a-daily-tauoff", "pooled",
              list(shared_tau = FALSE, ar_force = AR_BOAT_AG_DAILY),
              "2x2 cell (tau OFF, boat daily): the matched control stage F could not provide",
              "5.3"),
  S4b = stage("S4b", "S5-4b-daily-tauon", "pooled",
              list(shared_tau = TRUE, ar_force = AR_BOAT_AG_DAILY),
              "2x2 cell (tau ON, boat daily): do the two mechanisms compound, or overlap?",
              "5.3"),
  S5  = stage("S5", "S5-5-boatpc-ar", "pooled", list(ar_force = AR_BOAT_PC_BIWK),
              "boat pot closure at biweekly, uncontaminated this time (per-sub-season ar_force)",
              "5.6"))

resolve_cfg <- function(sid) modifyList(BASE, STAGE_DEFS[[sid]]$delta, keep.null = TRUE)

# ---------------------------------------------------------------------------
# Extraction. Every reader is guarded: a missing or renamed file yields NA rather than
# aborting a batch that has already spent hours.
# ---------------------------------------------------------------------------
find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}
rd <- function(dir, f) {
  p <- file.path(dir, f)
  if (!file.exists(p)) return(NULL)
  tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL)
}
num1 <- function(x) { x <- suppressWarnings(as.numeric(x)); if (length(x)) x[1] else NA_real_ }
chr1 <- function(x) { if (length(x)) as.character(x)[1] else NA_character_ }

# structural_params_*.csv, any column. `estimate` (added 2026-08-30) is `median` with the
# decoupled rows blanked, so it is the column to read when the number is going into a table
# a reader will treat as an estimate; `median` is kept for the raw value.
sp_get <- function(dir, fitlab, par, col = "median") {
  d <- rd(dir, sprintf("structural_params_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"parameter" %in% names(d) || !col %in% names(d)) return(NA_real_)
  num1(d[[col]][d$parameter == par])
}
sp_flag <- function(dir, fitlab, par, col = "decoupled") {
  d <- rd(dir, sprintf("structural_params_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !col %in% names(d)) return(NA)
  v <- d[[col]][d$parameter == par]; if (length(v)) as.logical(v[1]) else NA
}
# tau_bar is reported as "tau_bar[1]" (it is a length-shared_tau vector); accept either.
sp_tau <- function(dir, fitlab, col) {
  v <- sp_get(dir, fitlab, "tau_bar[1]", col)
  if (is.na(v)) v <- sp_get(dir, fitlab, "tau_bar", col)
  v
}
pit_get <- function(dir, fitlab, stream) {
  d <- rd(dir, sprintf("ppc_calibration_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"data_type" %in% names(d)) return(NA_real_)
  num1(d$pit_mean[d$data_type == stream])
}
# model_adequacy.csv when the run wrote one; model_adequacy_reconstructed.csv when it is an
# archived run that has been retro-annotated. Same columns, same code path, different source.
adq_get <- function(dir, fitlab, col) {
  d <- rd(dir, "model_adequacy.csv") %||% rd(dir, "model_adequacy_reconstructed.csv")
  if (is.null(d) || !"fit" %in% names(d) || !col %in% names(d)) return(NA)
  v <- d[[col]][d$fit == paste0(fitlab, "_Dungeness_Kept")]
  if (!length(v)) NA else v[1]
}
fds_get <- function(dir, fitlab, col) {
  d <- rd(dir, "fit_data_summary.csv")
  if (is.null(d) || !"fit" %in% names(d) || !col %in% names(d)) return(NA)
  v <- d[[col]][d$fit == paste0(fitlab, "_Dungeness_Kept")]
  if (!length(v)) NA else v[1]
}

extract_run <- function(sid, model, run_tag, outdir, minutes, ok) {
  row <- data.frame(
    stage = sid, plan_item = STAGE_DEFS[[sid]]$plan_item %||% NA_character_,
    model = model, run_tag = run_tag, ok = ok, minutes = minutes,
    outdir = if (is.na(outdir)) "" else basename(outdir),
    port_BSS_catch = NA_real_, port_BSS_lo95 = NA_real_, port_BSS_hi95 = NA_real_,
    port_PE_catch = NA_real_,
    shore_ag_BSS = NA_real_, shore_pc_BSS = NA_real_,
    boat_ag_BSS = NA_real_, boat_pc_BSS = NA_real_,
    # bss_reported (added 2026-08-30) is TRUE only where method_selected == "BSS". Without
    # it, a BSS_catch column carries the value of a REJECTED fit with nothing marking it:
    # the 2026-08-29 gear run printed BSS_catch = 32,689 for a component that reported PE.
    shore_ag_reported = NA, boat_ag_reported = NA, boat_pc_reported = NA,
    shore_ag_method = NA_character_, boat_ag_method = NA_character_,
    boat_pc_method = NA_character_,
    boat_ag_ar = NA_character_, boat_pc_ar = NA_character_, shore_ag_ar = NA_character_,
    # The shared-turnover feature and the floor that now guards it.
    shore_ag_n_informed = NA_real_, boat_ag_n_informed = NA_real_, boat_pc_n_informed = NA_real_,
    shore_ag_shared_tau = NA_real_, boat_ag_shared_tau = NA_real_, boat_pc_shared_tau = NA_real_,
    boat_tau_bar = NA_real_, boat_tau_lo95 = NA_real_, boat_tau_hi95 = NA_real_,
    boat_tau_neff = NA_real_, boat_tau_Rhat = NA_real_,
    # The absorption signature: a latent process eating the observation error shows up as
    # r_OSP running away with a huge upper bound and its parent's n_eff collapsing.
    boat_r_OSP = NA_real_, boat_r_OSP_hi95 = NA_real_, boat_sigma_r_OSP_neff = NA_real_,
    boat_kappa_OSP = NA_real_, boat_kappa_decoupled = NA,
    pit_trailer = NA_real_, pit_osp = NA_real_, pit_catch = NA_real_,
    boat_ag_div = NA_real_, boat_ag_neff = NA_real_,
    # Model adequacy (plan 5.5): reported beside the gate, never gating.
    boat_p_loo_frac = NA_real_, boat_p_loo_stream = NA_character_,
    boat_pareto_bad = NA_real_, boat_disp_neff = NA_real_, boat_disp_par = NA_character_,
    fits_passed = NA_character_, ar_resolutions = NA_character_,
    empty_effort_note = NA_character_, stringsAsFactors = FALSE)
  if (is.na(outdir) || !dir.exists(outdir)) return(row)

  pt <- rd(outdir, "port_total_Dungeness_Kept.csv")
  if (!is.null(pt) && "Estimate" %in% names(pt)) {
    is_catch <- pt$Estimate %in% c("Expected_Catch", "Catch")
    row$port_BSS_catch <- num1(pt$BSS_median[is_catch])
    row$port_BSS_lo95  <- num1(pt$BSS_lo95[is_catch])
    row$port_BSS_hi95  <- num1(pt$BSS_hi95[is_catch])
    row$port_PE_catch  <- num1(pt$PE[is_catch])
  }
  pv <- rd(outdir, "pe_vs_bss_comparison.csv")
  if (!is.null(pv) && "component" %in% names(pv)) {
    g <- function(pat, col) if (col %in% names(pv)) {
      v <- pv[[col]][grepl(pat, pv$component)]; if (length(v)) v[1] else NA } else NA
    row$shore_ag_BSS <- num1(g("^shore \\(All gear\\)", "BSS_catch"))
    row$shore_pc_BSS <- num1(g("^shore \\(Pot closure\\)", "BSS_catch"))
    row$boat_ag_BSS  <- num1(g("^private_boat \\(All gear\\)", "BSS_catch"))
    row$boat_pc_BSS  <- num1(g("^private_boat \\(Pot closure\\)", "BSS_catch"))
    row$shore_ag_reported <- as.logical(g("^shore \\(All gear\\)", "bss_reported"))
    row$boat_ag_reported  <- as.logical(g("^private_boat \\(All gear\\)", "bss_reported"))
    row$boat_pc_reported  <- as.logical(g("^private_boat \\(Pot closure\\)", "bss_reported"))
  }
  cr <- rd(outdir, "convergence_report.csv")
  if (!is.null(cr) && "fit" %in% names(cr)) {
    get <- function(pat, col) { if (!col %in% names(cr)) return(NA)
      i <- grepl(pat, cr$fit); if (any(i)) cr[[col]][i][1] else NA }
    row$shore_ag_method <- chr1(get("^shore_all_gear", "method_selected"))
    row$boat_ag_method  <- chr1(get("^private_boat_all_gear", "method_selected"))
    row$boat_pc_method  <- chr1(get("^private_boat_ring_net_only", "method_selected"))
    row$shore_ag_ar     <- chr1(get("^shore_all_gear", "ar_resolution"))
    row$boat_ag_ar      <- chr1(get("^private_boat_all_gear", "ar_resolution"))
    row$boat_pc_ar      <- chr1(get("^private_boat_ring_net_only", "ar_resolution"))
    row$boat_ag_div     <- num1(get("^private_boat_all_gear", "divergences"))
    row$boat_ag_neff    <- num1(get("^private_boat_all_gear", "C_sum_neff"))
  }
  # The AR resolution is also written to fit_data_summary.csv; the gear track's convergence
  # report has no ar_resolution column, so fall back rather than reporting NA.
  if (is.na(row$boat_ag_ar)) row$boat_ag_ar <- chr1(fds_get(outdir, "private_boat_all_gear", "ar_resolution"))
  if (is.na(row$boat_pc_ar)) row$boat_pc_ar <- chr1(fds_get(outdir, "private_boat_ring_net_only", "ar_resolution"))

  for (nm in c("shore_all_gear", "private_boat_all_gear", "private_boat_ring_net_only")) {
    key <- c(shore_all_gear = "shore_ag", private_boat_all_gear = "boat_ag",
             private_boat_ring_net_only = "boat_pc")[[nm]]
    row[[paste0(key, "_n_informed")]] <- num1(fds_get(outdir, nm, "n_L_informed"))
    row[[paste0(key, "_shared_tau")]] <- num1(fds_get(outdir, nm, "shared_tau"))
  }

  bag <- "private_boat_all_gear"
  row$boat_tau_bar   <- sp_tau(outdir, bag, "median")
  row$boat_tau_lo95  <- sp_tau(outdir, bag, "lo95")
  row$boat_tau_hi95  <- sp_tau(outdir, bag, "hi95")
  row$boat_tau_neff  <- sp_tau(outdir, bag, "n_eff")
  row$boat_tau_Rhat  <- sp_tau(outdir, bag, "Rhat")
  row$boat_r_OSP     <- sp_get(outdir, bag, "r_OSP", "median")
  row$boat_r_OSP_hi95 <- sp_get(outdir, bag, "r_OSP", "hi95")
  row$boat_sigma_r_OSP_neff <- sp_get(outdir, bag, "sigma_r_OSP", "n_eff")
  row$boat_kappa_OSP <- sp_get(outdir, bag, "kappa_OSP", "median")
  row$boat_kappa_decoupled <- sp_flag(outdir, bag, "kappa_OSP")
  row$pit_trailer <- pit_get(outdir, bag, "trailer")
  row$pit_osp     <- pit_get(outdir, bag, "osp")
  row$pit_catch   <- pit_get(outdir, bag, "catch")
  row$boat_p_loo_frac   <- num1(adq_get(outdir, bag, "p_loo_frac"))
  row$boat_p_loo_stream <- chr1(adq_get(outdir, bag, "p_loo_worst_stream"))
  row$boat_pareto_bad   <- num1(adq_get(outdir, bag, "n_pareto_bad"))
  row$boat_disp_neff    <- num1(adq_get(outdir, bag, "disp_neff_min"))
  row$boat_disp_par     <- chr1(adq_get(outdir, bag, "disp_neff_min_par"))

  ss <- rd(outdir, "season_summary.csv")
  if (!is.null(ss) && "metric" %in% names(ss)) {
    row$fits_passed    <- chr1(ss$value[ss$metric == "BSS fits passed"])
    row$ar_resolutions <- chr1(ss$value[ss$metric == "AR temporal resolution"])
  }
  ee <- rd(outdir, "pe_empty_effort_strata.csv")
  if (!is.null(ee) && "component" %in% names(ee) && nrow(ee) > 0)
    row$empty_effort_note <- paste(sprintf("%s %s/%s", ee$component, ee$n_empty_days,
                                           ee$n_calendar_days), collapse = "; ")
  row
}

# ---------------------------------------------------------------------------
# EXACTNESS. Compare a run's FITS against a reference folder at full double precision.
#
# `fit_pattern` restricts the comparison to the fits whose label matches, which is what
# makes the S2 gate sharp. The informed-day floor is a claim about WHICH fits the feature
# reaches, so the right test is not "the shore total came back near 20,898" but "the shore
# fits are the same object they were with the feature globally off, and the boat fits are
# the same object they were with it globally on". Two exact statements instead of two
# tolerances, on the same run.
#
# Shared parameter ROWS only, because a feature legitimately ADDS reported rows (tau_bar,
# tau_bar_out; f_lower_out did the same in the 2026-08-25 batch). New rows are not a
# behaviour change; changed values are.
#
# VERIFIED PRECONDITION (2026-08-30, measured not assumed): extra elements in the list
# passed to rstan::stan(data = ) are inert. A two-chain fit run with, and without, three
# added elements including the new `.n_L_informed` returned bit-identical draw matrices
# (identical() TRUE on the full 600 x 3 matrix). So the dotted metadata that prep_bss_crab_*
# attaches to stan_data cannot be what moves a gate below.
# ---------------------------------------------------------------------------
fit_exactness <- function(new_dirname, ref_dir, fit_pattern = NULL, what = "fits") {
  if (is.na(new_dirname) || !nzchar(new_dirname))
    return(list(observed = "no output folder found", verdict = "REVIEW"))
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == new_dirname]
  new_dir <- if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NULL
  if (is.null(new_dir) || !dir.exists(new_dir) || !dir.exists(ref_dir))
    return(list(observed = "run or reference folder missing", verdict = "REVIEW"))

  keep <- function(fs) if (is.null(fit_pattern)) fs else fs[grepl(fit_pattern, fs)]
  fs <- keep(c(list.files(new_dir, pattern = "^bss_summary_.*\\.csv$"),
               list.files(new_dir, pattern = "^bss_full_summary_.*\\.csv$")))
  fs <- fs[file.exists(file.path(ref_dir, fs))]
  rows_ok <- TRUE; n_rows <- 0L; n_diff <- 0L; worst <- ""
  for (f in fs) {
    a <- tryCatch(utils::read.csv(file.path(ref_dir, f), row.names = 1, check.names = FALSE),
                  error = function(e) NULL)
    b <- tryCatch(utils::read.csv(file.path(new_dir, f), row.names = 1, check.names = FALSE),
                  error = function(e) NULL)
    if (is.null(a) || is.null(b)) { rows_ok <- FALSE; next }
    common <- intersect(rownames(a), rownames(b)); cols <- intersect(names(a), names(b))
    n_rows <- n_rows + length(common)
    cmp <- all.equal(a[common, cols], b[common, cols], tolerance = 0)
    if (!isTRUE(cmp)) {
      rows_ok <- FALSE; n_diff <- n_diff + 1L
      if (!nzchar(worst)) worst <- sprintf(" first difference in %s: %s", f, cmp[1])
    }
  }
  # Convergence-report rows for the same fits. Byte-comparing the whole file is only valid
  # when every fit is in scope; for a subset, compare the matching rows.
  cr_note <- ""
  a <- tryCatch(utils::read.csv(file.path(ref_dir, "convergence_report.csv"), stringsAsFactors = FALSE),
                error = function(e) NULL)
  b <- tryCatch(utils::read.csv(file.path(new_dir, "convergence_report.csv"), stringsAsFactors = FALSE),
                error = function(e) NULL)
  cr_same <- NA
  if (!is.null(a) && !is.null(b) && "fit" %in% names(a) && "fit" %in% names(b)) {
    ia <- if (is.null(fit_pattern)) rep(TRUE, nrow(a)) else grepl(fit_pattern, a$fit)
    ib <- if (is.null(fit_pattern)) rep(TRUE, nrow(b)) else grepl(fit_pattern, b$fit)
    cols <- intersect(names(a), names(b))
    cr_same <- isTRUE(all.equal(a[ia, cols], b[ib, cols], tolerance = 0, check.attributes = FALSE))
    cr_note <- sprintf("; convergence rows (%d) %s", sum(ia), if (isTRUE(cr_same)) "identical" else "DIFFER")
  }
  list(observed = sprintf("%s: %d shared parameter rows across %d summaries %s%s%s",
                          what, n_rows, length(fs),
                          if (rows_ok) "identical at full precision"
                          else sprintf("DIFFER in %d summaries", n_diff),
                          cr_note, worst),
       verdict = if (rows_ok && isTRUE(cr_same) && n_rows > 0) "PASS" else "FAIL")
}

# ---------------------------------------------------------------------------
# STAGE S1: the empty-stratum fill, PE only (plan 5.1).
#
# pe_empty_effort_stratum changes ONLY the Point Estimator, so a full render would spend
# ~4 h of MCMC to answer a question the PE settles in a minute. This builds the same inputs
# the driver builds, INCLUDING the crabbing-holiday workbook, and calls run_pe_pooled()
# twice per component.
#
# The holiday load is the 2026-08-27b fix and it is the reason this is a re-run: day-typing
# without holidays put 289 days into 85 week x day-type cells instead of 92 and shifted every
# zeroed-day count. The re-run's `zero` arm must therefore now reproduce the DRIVER's PE
# exactly, component by component, which is a criterion the original stage A could not state.
# ---------------------------------------------------------------------------
run_stage_S1 <- function() {
  banner("STAGE S1  |  empty-effort-stratum fill, PE only (plan item 5.1)")
  out_dir <- .here("05_output", format(Sys.Date(), "%Y%m%d"), "stage5-S1-pe-fill")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  params <- modifyList(BASE, list(bss_model_file = "crab_bss_pooled.stan"))
  params$crabbing_holiday_dates <- read_crabbing_holidays(params)   # <- the fix
  dwg     <- fetch_crab_data(params)
  ie_data <- fetch_ie_data(params)
  # Mirrors the driver exactly: the day-length model is built only when there is I/E data,
  # and NULL falls back to civil twilight the same way.
  L_eff   <- if (!is.null(ie_data) && nrow(ie_data) > 0)
               estimate_L_effective(ie_data, params$pot_open_date, params) else NULL
  subs    <- build_subseasons(params)

  rows <- list()
  for (fill in c("zero", "day_type")) {
    p <- modifyList(params, list(pe_empty_effort_stratum = fill))
    pe_all <- list()
    for (pop in c("shore", "private_boat")) for (ss in subs) {
      lab <- paste0(pop, "_", ss$name)
      days_ss <- prep_days_crab(ss$start, ss$end, p, L_eff_model = L_eff)
      summ_ss <- prep_population_summary(dwg, pop, ss$start, ss$end, p)
      r <- run_pe_pooled(summ_ss, days_ss, p, lab)
      pe_all[[lab]] <- r
      rows[[length(rows) + 1]] <- data.frame(
        fill = fill, component = lab,
        pe_effort = r$effort_total %||% NA_real_,
        pe_catch  = r[["Dungeness_Kept"]] %||% NA_real_,
        n_empty_strata  = r$n_empty_effort_strata %||% NA_integer_,
        n_strata_total  = r$n_effort_strata_total %||% NA_integer_,
        n_empty_days    = r$n_empty_effort_days   %||% NA_integer_,
        n_calendar_days = r$n_calendar_days       %||% NA_integer_,
        stringsAsFactors = FALSE)
    }
    write_pe_empty_stratum_report(pe_all, out_dir, p)
    if (identical(fill, "day_type"))
      file.rename(file.path(out_dir, "pe_empty_effort_strata.csv"),
                  file.path(out_dir, "pe_empty_effort_strata_day_type.csv"))
  }
  df <- do.call(rbind, rows)
  df$empty_day_fraction <- df$n_empty_days / pmax(df$n_calendar_days, 1)
  w <- merge(df[df$fill == "zero", ], df[df$fill == "day_type", ],
             by = "component", suffixes = c("_zero", "_daytype"))
  w$catch_change_pct  <- round(100 * (w$pe_catch_daytype  - w$pe_catch_zero)  / pmax(w$pe_catch_zero, 1), 1)
  w$effort_change_pct <- round(100 * (w$pe_effort_daytype - w$pe_effort_zero) / pmax(w$pe_effort_zero, 1), 1)
  utils::write.csv(df, file.path(out_dir, "pe_fill_comparison_long.csv"), row.names = FALSE)
  utils::write.csv(w,  file.path(out_dir, "pe_fill_comparison.csv"), row.names = FALSE)

  # Does the fixed stage now agree with the driver? Exact match expected on both.
  ref_pe   <- c(shore_all_gear = REF$rung4$pe_shore_ag, shore_ring_net_only = REF$rung4$pe_shore_pc,
                private_boat_all_gear = REF$rung4$pe_boat_ag,
                private_boat_ring_net_only = REF$rung4$pe_boat_pc)
  chk <- data.frame(component = w$component,
                    empty_days_here = w$n_empty_days_zero,
                    empty_days_driver = as.integer(REF$rung4$empty_days[w$component]),
                    pe_catch_here = round(w$pe_catch_zero),
                    pe_catch_driver = as.integer(ref_pe[w$component]), stringsAsFactors = FALSE)
  chk$days_match <- chk$empty_days_here == chk$empty_days_driver
  chk$pe_match   <- abs(chk$pe_catch_here - chk$pe_catch_driver) <= 1
  utils::write.csv(chk, file.path(out_dir, "pe_fill_vs_driver_check.csv"), row.names = FALSE)

  rule()
  cat("  component                    zeroed days   PE catch (zero -> day_type)   change\n")
  for (i in seq_len(nrow(w)))
    cat(sprintf("  %-28s %3d/%-3d %5.1f%%   %8s -> %-8s  %+6.1f%%\n",
                w$component[i], w$n_empty_days_zero[i], w$n_calendar_days_zero[i],
                100 * w$empty_day_fraction_zero[i],
                fmt(w$pe_catch_zero[i], 0), fmt(w$pe_catch_daytype[i], 0), w$catch_change_pct[i]))
  rule()
  cat("  agreement with the driver's own PE (the 2026-08-27b holiday fix):\n")
  for (i in seq_len(nrow(chk)))
    cat(sprintf("    %-28s days %2d vs %2d %-5s   PE %6d vs %6d %s\n",
                chk$component[i], chk$empty_days_here[i], chk$empty_days_driver[i],
                if (chk$days_match[i]) "ok" else "MISMATCH",
                chk$pe_catch_here[i], chk$pe_catch_driver[i],
                if (chk$pe_match[i]) "ok" else "MISMATCH"))
  cat("\n  Written to", basename(out_dir), "\n")
  invisible(list(w = w, chk = chk))
}

# ---------------------------------------------------------------------------
# The run loop.
# ---------------------------------------------------------------------------
sum_path <- .here("05_output", "stage5_2026-08-30_summary.csv")
ver_path <- .here("05_output", "stage5_2026-08-30_verdicts.csv")
tab_path <- .here("05_output", "stage5_2026-08-30_tau_by_ar.csv")
append_row <- function(r) utils::write.table(
  r, sum_path, sep = ",", row.names = FALSE, qmethod = "double",
  col.names = !file.exists(sum_path), append = file.exists(sum_path))

# ---- DRY RUN ---------------------------------------------------------------
if (isTRUE(DRY_RUN)) {
  banner(sprintf("DRY RUN  |  stages = %s  |  NOTHING WILL BE FITTED", paste(STAGES, collapse = ", ")))
  for (sid in STAGES) {
    S <- STAGE_DEFS[[sid]]
    rule()
    cat(sprintf("  STAGE %-3s (%s)  plan item %s  model: %s\n", S$id,
                if (identical(S$kind, "pe_only")) "PE only, no fitting" else "full run",
                S$plan_item, paste(S$model, collapse = " + ")))
    cat(sprintf("    %s\n", S$headline))
    cfg <- resolve_cfg(sid)
    keys <- names(S$delta)
    if (!length(keys)) cat("    delta: none (shipped config)\n")
    for (k in keys) {
      v <- cfg[[k]]
      cat(sprintf("    delta: %-24s %s\n", k,
                  if (is.null(v)) "NULL" else paste(utils::capture.output(str(v)), collapse = " ")))
    }
  }
  rule()
  banner("SELF-TEST (extractor + gates against the archived runs; no fitting)")
  st_ok <- 0; st_bad <- 0
  st <- function(nm, cond, extra = "") {
    if (isTRUE(cond)) { st_ok <<- st_ok + 1; cat(sprintf("  PASS  %s %s\n", nm, extra)) }
    else              { st_bad <<- st_bad + 1; cat(sprintf("  FAIL  %s %s\n", nm, extra)) } }

  dD <- .here("05_output", REF$D$dir); dE <- .here("05_output", REF$E$dir)
  dF <- .here("05_output", REF$F$dir); dG <- .here("05_output", REF$Egear$dir)

  if (dir.exists(dD)) {
    r <- extract_run("S2", "pooled", "ref", dD, NA_real_, TRUE)
    st("stage D (tau off, monthly) components parse",
       isTRUE(r$shore_ag_BSS == REF$D$shore_ag) && isTRUE(r$boat_ag_BSS == REF$D$boat_ag),
       sprintf("shore %s boat %s", fmt(r$shore_ag_BSS,0), fmt(r$boat_ag_BSS,0)))
    st("boat all-gear AR parses (the 2x2 axis)", identical(r$boat_ag_ar, REF$D$boat_ag_ar), r$boat_ag_ar)
    st("tau_bar absent when the feature is off", is.na(r$boat_tau_bar))
    st("kappa_OSP reads as decoupled under osp_scale_is_tau", isTRUE(r$boat_kappa_decoupled))
    st("adequacy reconstructed for the archived cell",
       isTRUE(abs(r$boat_p_loo_frac - REF$D$p_loo_frac) < 0.005),
       sprintf("p_loo_frac %s (%s)", r$boat_p_loo_frac, r$boat_p_loo_stream))
    st("n_L_informed absent in the archived runs (added 2026-08-30)", is.na(r$boat_ag_n_informed))
    ex <- fit_exactness(basename(dD), dD)
    st("exactness helper compares a folder to itself", identical(ex$verdict, "PASS"), ex$observed)
    exs <- fit_exactness(basename(dD), dD, "shore")
    st("exactness helper restricts to a fit subset", identical(exs$verdict, "PASS"), exs$observed)
    exd <- fit_exactness(basename(dE), dD, "private_boat_all_gear")
    st("exactness helper DETECTS a real difference (E vs D on the boat)",
       identical(exd$verdict, "FAIL"), substr(exd$observed, 1, 96))
    exsh <- fit_exactness(basename(dE), dD, "shore_all_gear")
    st("...and stage E's shore DID move, which is why the floor exists",
       identical(exsh$verdict, "FAIL"), substr(exsh$observed, 1, 72))
  } else st("stage D folder present", FALSE, dD)

  if (dir.exists(dE)) {
    r <- extract_run("S2", "pooled", "ref", dE, NA_real_, TRUE)
    st("stage E tau_bar parses with its interval",
       isTRUE(abs(r$boat_tau_bar - REF$E$tau_bar) < 0.005),
       sprintf("%s [%s, %s]", fmt(r$boat_tau_bar,3), fmt(r$boat_tau_lo95,3), fmt(r$boat_tau_hi95,3)))
    st("stage E PIT means parse", isTRUE(abs(r$pit_trailer - REF$E$pit_trailer) < 0.005) &&
       isTRUE(abs(r$pit_osp - REF$E$pit_osp) < 0.005),
       sprintf("trailer %s osp %s", fmt(r$pit_trailer,3), fmt(r$pit_osp,3)))
  } else st("stage E folder present", FALSE, dE)

  if (dir.exists(dF)) {
    r <- extract_run("S4a", "pooled", "ref", dF, NA_real_, TRUE)
    st("stage F r_OSP runaway parses (the absorption signature)",
       isTRUE(r$boat_r_OSP_hi95 > 1000) && isTRUE(r$boat_sigma_r_OSP_neff < 400),
       sprintf("r_OSP %s hi95 %s, sigma_r_OSP n_eff %s", fmt(r$boat_r_OSP,2),
               fmt(r$boat_r_OSP_hi95,0), fmt(r$boat_sigma_r_OSP_neff,0)))
    st("stage F adequacy separates from stage D",
       isTRUE(r$boat_p_loo_frac > 0.15) && isTRUE(r$boat_pareto_bad >= 2),
       sprintf("p_loo_frac %s (%s), Pareto k>0.7 %s", r$boat_p_loo_frac,
               r$boat_p_loo_stream, r$boat_pareto_bad))
    st("stage F ran the boat at daily", identical(r$boat_ag_ar, "daily"), r$boat_ag_ar)
    st("stage F moved the SHORE pot closure too, so it is not a matched control",
       !isTRUE(r$shore_pc_BSS == REF$D$shore_pc),
       sprintf("shore pot closure %s vs stage D %s", fmt(r$shore_pc_BSS,0), fmt(REF$D$shore_pc,0)))
  } else st("stage F folder present", FALSE, dF)

  if (dir.exists(dG)) {
    r <- extract_run("S3", "gear_resolved", "ref", dG, NA_real_, TRUE)
    st("gear boat all-gear reads as REJECTED",
       grepl("^PE", r$boat_ag_method %||% ""), r$boat_ag_method)
    st("gear tau_bar n_eff is the problem S3 addresses",
       isTRUE(r$boat_tau_neff < 400) && isTRUE(r$boat_tau_Rhat > 1.01),
       sprintf("tau_bar n_eff %s Rhat %s", fmt(r$boat_tau_neff,0), fmt(r$boat_tau_Rhat,4)))
    st("bss_reported is absent in the archived gear run, which is the defect it fixes",
       is.na(r$boat_ag_reported),
       sprintf("BSS_catch %s printed beside method %s", fmt(r$boat_ag_BSS,0), r$boat_ag_method))
  } else st("gear stage E folder present", FALSE, dG)

  st("S3 sampler override passes the whitelist",
     !inherits(try(bss_apply_sampler_override(list(), S3_SAMPLER, "dry", quiet = TRUE),
                   silent = TRUE), "try-error"))
  st("a structural key is REFUSED by the override, not silently dropped",
     inherits(try(bss_apply_sampler_override(list(), list(shared_tau = TRUE), "dry", quiet = TRUE),
                  silent = TRUE), "try-error"))
  st("the nested ar_force resolves per sub-season",
     identical(.bss_resolve_ar_force(list(ar_force = AR_BOAT_AG_DAILY), "private_boat", "all_gear"), "daily") &&
     is.null(.bss_resolve_ar_force(list(ar_force = AR_BOAT_AG_DAILY), "private_boat", "pot_closure")),
     "all_gear -> daily, pot_closure -> unforced")

  banner(sprintf("SELF-TEST: %d passed, %d failed", st_ok, st_bad))
  cat("\n  Set DRY_RUN <- FALSE and source again to start.\n\n")
} else {

# ---- REAL RUN --------------------------------------------------------------
banner(sprintf("STAGE 5 BATCH  |  stages = %s  |  start %s",
               paste(STAGES, collapse = ", "), format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("  Baselines: rung 4 (port 66,237) and the 2026-08-29 batch (D/E/F).\n")
cat("  Each stage is an ISOLATED change against rung 4, not cumulative.\n")

# Put the archived 2x2 cells on the same footing as the new ones before anything runs, so a
# batch that dies at hour 3 still leaves a comparable table behind.
banner("Reconstructing model adequacy for the archived cells (plan 5.5)")
for (d in c(REF$rung4$dir, REF$D$dir, REF$E$dir, REF$F$dir, REF$Egear$dir, REF$Cbad$dir))
  try(annotate_model_adequacy_run(.here("05_output", d)), silent = TRUE)
if (exists("annotate_decoupled_run"))
  for (d in c(REF$rung4$dir, REF$D$dir, REF$E$dir, REF$F$dir, REF$Egear$dir, REF$Cbad$dir))
    try(annotate_decoupled_run(.here("05_output", d)), silent = TRUE)

collected <- list()
skip_S4b <- FALSE
gate_S2  <- NULL

for (sid in STAGES) {
  S <- STAGE_DEFS[[sid]]

  if (identical(S$kind, "pe_only")) {
    t0 <- Sys.time()
    ok <- tryCatch({ assign("S1_RESULT", run_stage_S1(), envir = globalenv()); TRUE },
                   error = function(e) { message("*** STAGE S1 FAILED: ", conditionMessage(e)); FALSE })
    banner(sprintf("STAGE S1 -> %s (%.1f min)", if (ok) "OK" else "FAILED",
                   as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    next
  }

  if (identical(sid, "S4b") && isTRUE(skip_S4b)) {
    banner(paste("STAGE S4b SKIPPED: the S2 gate did not hold, so `shared_tau = TRUE` is not",
                 "yet a boat-only change and an (ON, daily) cell could not be attributed to",
                 "the AR resolution. S4a still ran and still gives you the matched (OFF,",
                 "daily) cell that stage F could not. Fix the floor, then run S4b alone."))
    next
  }

  m <- S$model
  run_tag <- S$tag
  folder  <- paste0(prefix[[m]], run_tag)

  done <- find_outdir(m, run_tag)
  if (isTRUE(RESUME) && !is.na(done) &&
      file.exists(file.path(done, "port_total_Dungeness_Kept.csv"))) {
    banner(sprintf("SKIP (already complete): stage %s / %s / %s", sid, m, folder))
    r <- extract_run(sid, m, run_tag, done, NA_real_, TRUE)
    collected[[length(collected)+1]] <- r; append_row(r)
    if (identical(sid, "S2")) gate_S2 <- NULL   # re-evaluated below from the folder
  } else {
    rc <- resolve_cfg(sid)
    rc$model <- m; rc$run_tag <- run_tag; rc$run_weather <- FALSE

    banner(sprintf("STAGE %s  %s  |  %s  |  start %s", sid, folder, S$headline,
                   format(Sys.time(), "%H:%M:%S")))
    cat(sprintf("  plan item %s | shared_tau=%s | shared_tau_min_obs=%s | ar_force=%s | sampler_override=%s\n",
                S$plan_item, rc$shared_tau %||% FALSE, rc$shared_tau_min_obs %||% 15,
                if (is.null(rc$ar_force)) "NULL" else paste(utils::capture.output(str(rc$ar_force)), collapse = " "),
                if (is.null(rc$bss_sampler_override)) "NULL"
                else paste(names(rc$bss_sampler_override), unlist(rc$bss_sampler_override),
                           sep = "=", collapse = " ")))

    run_env <- new.env(parent = globalenv())
    run_env$run_config <- rc
    t0 <- Sys.time()
    ok <- tryCatch({
      html <- rmarkdown::render(model_rmd[[m]], envir = run_env, quiet = FALSE)
      od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE),
                     error = function(e) NA_character_)
      if (!is.na(od) && dir.exists(od) &&
          normalizePath(dirname(html)) != normalizePath(od)) {
        if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE)))
          suppressWarnings(file.remove(html))
      }
      TRUE
    }, error = function(e) {
      message("*** STAGE ", sid, " / ", m, " FAILED: ", conditionMessage(e))
      # The render error is often NOT the cause: rstan returns an empty stanfit rather than
      # raising, so a chain failure surfaces downstream. Point at what the sampler said.
      od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE),
                     error = function(e2) NA_character_)
      if (!is.na(od) && dir.exists(od)) {
        logs <- list.files(od, pattern = "^stan_console_.*\\.log$", full.names = TRUE)
        if (length(logs)) {
          last <- logs[which.max(file.mtime(logs))]
          message("    Stan console for the last fit attempted (", basename(last), "):")
          message(paste0("      ", utils::tail(readLines(last, warn = FALSE), 8), collapse = "\n"))
        } else message("    No stan_console_*.log in ", od,
                       " -- failed before or during the first sampler call.")
      }
      FALSE })
    mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

    r <- extract_run(sid, m, run_tag, find_outdir(m, run_tag), mins, ok)
    collected[[length(collected)+1]] <- r
    tryCatch(append_row(r), error = function(e)
      message("  (summary row not written: ", conditionMessage(e), ")"))
    banner(sprintf("STAGE %s / %s -> %s (%s min)", sid, m, if (ok) "OK" else "FAILED", mins))
  }

  # ---- THE S2 GATE --------------------------------------------------------
  # Two exact statements about one run:
  #   the SHORE fits must be the object they were with the feature globally OFF (stage D),
  #   the BOAT  fits must be the object they were with the feature globally ON  (stage E).
  # Together those say the floor changed which fits the feature reaches and changed nothing
  # else. A tolerance-based "the shore came back near 20,898" would pass a run where the
  # shore moved for some other reason and happened to land nearby.
  if (identical(sid, "S2")) {
    nd <- basename(find_outdir(m, run_tag) %||% NA_character_)
    exsh <- fit_exactness(nd, .here("05_output", REF$D$dir), "shore", "shore fits vs stage D (tau OFF)")
    exbt <- fit_exactness(nd, .here("05_output", REF$E$dir), "private_boat", "boat fits vs stage E (tau ON)")
    gate_S2 <- list(shore = exsh, boat = exbt)
    cat("\n  GATE  informed-day floor\n")
    cat("        shore: ", exsh$verdict, " -- ", exsh$observed, "\n", sep = "")
    cat("        boat:  ", exbt$verdict, " -- ", exbt$observed, "\n", sep = "")
    if (!identical(exsh$verdict, "PASS") || !identical(exbt$verdict, "PASS")) {
      skip_S4b <- TRUE
      message("*** The floor did not do exactly what it claims. S4b will be skipped; S4a is ",
              "unaffected because it runs with shared_tau OFF.")
    }
  }
}

# ---- THE 2x2 -------------------------------------------------------------
# The table Stage 5 exists to produce. Archived cells are read from their committed folders;
# new cells from this batch. Every cell carries the four adequacy quantities plan 5.3 asks
# for, which is why the archived folders were retro-annotated above.
S <- if (length(collected)) do.call(rbind, collected) else NULL
g1 <- function(sid, col) { if (is.null(S)) return(NA)
  x <- S[[col]][S$stage == sid]; if (length(x)) x[1] else NA }

cell <- function(name, tau, ar, src, from_dir = NULL, sid = NULL) {
  r <- if (!is.null(from_dir)) {
    d <- .here("05_output", from_dir)
    if (dir.exists(d)) extract_run("S2", "pooled", "archived", d, NA_real_, TRUE) else NULL
  } else if (!is.null(sid) && !is.null(S) && any(S$stage == sid)) S[S$stage == sid, ][1, ] else NULL
  if (is.null(r)) return(NULL)
  data.frame(cell = name, shared_tau = tau, boat_ar = ar, source = src,
             boat_all_gear = r$boat_ag_BSS, reported = r$boat_ag_reported,
             port = r$port_BSS_catch,
             tau_bar = r$boat_tau_bar, r_OSP = r$boat_r_OSP, r_OSP_hi95 = r$boat_r_OSP_hi95,
             sigma_r_OSP_neff = r$boat_sigma_r_OSP_neff,
             p_loo_frac = r$boat_p_loo_frac, p_loo_stream = r$boat_p_loo_stream,
             pareto_bad = r$boat_pareto_bad, disp_neff = r$boat_disp_neff,
             pit_trailer = r$pit_trailer, pit_osp = r$pit_osp, pit_catch = r$pit_catch,
             divergences = r$boat_ag_div, stringsAsFactors = FALSE)
}
TAB <- do.call(rbind, Filter(Negate(is.null), list(
  cell("OFF x monthly", FALSE, "monthly", "archived: 2026-08-29 stage D", from_dir = REF$D$dir),
  cell("ON  x monthly", TRUE,  "monthly", "this batch: S2",               sid = "S2"),
  cell("OFF x daily",   FALSE, "daily",   "this batch: S4a",              sid = "S4a"),
  cell("ON  x daily",   TRUE,  "daily",   "this batch: S4b",              sid = "S4b"),
  cell("(unmatched) escalated, all four fits daily", FALSE, "daily",
       "archived: 2026-08-29 stage F, DIFFERENT code path", from_dir = REF$F$dir))))
if (!is.null(TAB)) utils::write.csv(TAB, tab_path, row.names = FALSE)

# ---- VERDICTS ------------------------------------------------------------
V <- list()
V1 <- function(stage, criterion, observed, threshold, verdict, why)
  data.frame(stage = stage, criterion = criterion, observed = observed,
             threshold = threshold, verdict = verdict, why = why, stringsAsFactors = FALSE)

# ---- S1
if (exists("S1_RESULT")) {
  chk <- S1_RESULT$chk; w <- S1_RESULT$w
  V[[length(V)+1]] <- V1("S1", "the fixed PE stage reproduces the driver's own Point Estimator",
    sprintf("zeroed days %s; PE %s",
            paste(sprintf("%s %d/%d", sub("_.*", "", chk$component), chk$empty_days_here,
                          chk$empty_days_driver), collapse = ", "),
            if (all(chk$pe_match)) "matches on all four components" else "MISMATCH"),
    "exact agreement on all four components, both quantities",
    if (all(chk$days_match) && all(chk$pe_match)) "PASS" else "FAIL",
    paste("The 2026-08-28 stage A did not load crabbing_holidays.xlsx, so it stratified 289",
          "days into 85 week x day-type cells where the driver builds 92, and reported",
          "41/44/8/0 zeroed days against the driver's 44/47/9/0 -- it even assigned the 44",
          "to the wrong component. The zero-vs-day_type COMPARISON was still fair, because",
          "both arms shared the stratification, but the absolute numbers were not",
          "production numbers and could not be quoted. A FAIL here means the holiday fix is",
          "incomplete and every fill number in the plan is still provisional."))
  V[[length(V)+1]] <- V1("S1", "what adopting day_type would cost",
    paste(sprintf("%s %+.1f%%", sub("_Dungeness.*", "", w$component), w$catch_change_pct),
          collapse = "; "),
    "decision, not a threshold", "DECIDE",
    paste("`zero` expands unsampled effort strata at zero effort, which asserts that",
          "unsampled days had no crabbers on them. The unsampled strata skew weekend and",
          "holiday and those days carry 1.7-2.3x weekday effort, so `zero` is biased DOWN",
          "by construction and only the size is in question. The counter-argument is",
          "continuity with every published number to date, and it is real. Note where this",
          "actually changes a REPORTED number: a component whose gate passed reports its",
          "BSS, so the fill moves only the cross-check; a component on PE fallback reports",
          "this figure directly."))
}

# ---- S2, the gate
if (!is.null(gate_S2)) {
  V[[length(V)+1]] <- V1("S2", "the informed-day floor REFUSES the shore",
    gate_S2$shore$observed, "shore fits bit-identical to stage D (shared_tau globally OFF)",
    gate_S2$shore$verdict,
    paste("The shore all-gear fit has 4 in-window I/E days out of 289 and its pot-closure",
          "fit has none. In the 2026-08-29 batch a global shared_tau moved that component",
          "+17.9% (20,898 -> 24,629) on those four observations, with a tau_bar interval",
          "still containing the 1.700 prior centre, no measurable improvement in fit, and",
          "no replication in the gear track, which put the same parameter at 1.681. The",
          "floor should now leave the shore untouched. Bit-identity, not a tolerance:",
          "if the shore moved at all, the floor is not the only thing that changed."))
  V[[length(V)+1]] <- V1("S2", "the informed-day floor ALLOWS the boat",
    gate_S2$boat$observed, "boat fits bit-identical to stage E (shared_tau globally ON)",
    gate_S2$boat$verdict,
    paste("The boat is why the feature exists: 130 OSP days inform its level once",
          "osp_scale_is_tau puts L into the OSP mean, tau_bar lands at 2.597 [2.064, 3.249]",
          "excluding the 1.200 prior centre, it agrees with the independent OSP/trailer",
          "overlap calibration (2.01-3.03), and it pulls the trailer and OSP PIT means",
          "toward nominal simultaneously and in opposite directions. The boat pot closure",
          "has 18 informed days and also clears the default floor of 15. A FAIL here means",
          "the floor is catching a fit it should not, or something else moved with it."))
  ni <- sprintf("shore all-gear %s informed / shared_tau=%s; boat all-gear %s / %s; boat pot closure %s / %s",
                fmt(g1("S2","shore_ag_n_informed"),0), g1("S2","shore_ag_shared_tau"),
                fmt(g1("S2","boat_ag_n_informed"),0),  g1("S2","boat_ag_shared_tau"),
                fmt(g1("S2","boat_pc_n_informed"),0),  g1("S2","boat_pc_shared_tau"))
  V[[length(V)+1]] <- V1("S2", "the floor's own bookkeeping agrees with its effect",
    ni, "shore 4 informed with shared_tau 0; boat 130 and 18 with shared_tau 1",
    if (isTRUE(g1("S2","shore_ag_shared_tau") == 0) &&
        isTRUE(g1("S2","boat_ag_shared_tau") == 1) &&
        isTRUE(g1("S2","boat_pc_shared_tau") == 1)) "PASS" else "REVIEW",
    paste("Direct evidence rather than inference from totals. n_L_informed and shared_tau",
          "are now written per fit to fit_data_summary.csv, so the decision the floor made",
          "is auditable from the output folder without re-reading the console. If the two",
          "exactness gates above pass but this disagrees, the totals coincided for another",
          "reason and both gates are lucky rather than right."))
}

# ---- S3, the cross-track check
if (!is.null(S) && any(S$stage == "S3")) {
  meth <- g1("S3","boat_ag_method"); tn <- g1("S3","boat_tau_neff"); tr <- g1("S3","boat_tau_Rhat")
  V[[length(V)+1]] <- V1("S3", "the gear boat all-gear fit converges with the pooled track's sampler",
    sprintf("method %s; tau_bar n_eff %s Rhat %s (was 49 / 1.0813); divergences %s (was 554)",
            meth, fmt(tn,0), fmt(tr,4), fmt(g1("S3","boat_ag_div"),0)),
    "method_selected == BSS",
    if (identical(meth, "BSS")) "PASS" else "FAIL",
    paste("The 2026-08-29 gear run rejected this fit at n_eff 25-49 and R-hat 1.08-1.17 with",
          "554 divergences, so the batch got no cross-track verdict on the shared turnover.",
          "This run gives that fit the pooled track's settings for the same component",
          "(2,500 draws/chain, adapt_delta 0.99, max_treedepth 13) where it previously ran",
          "at 1,000 / 0.90 / 10 -- the gear driver has per-fit overrides for shore-all-gear",
          "and both pot-closure fits and NONE for boat all-gear. THREE THINGS MOVED, so a",
          "PASS does not attribute the fix; it buys a converged fit to compare. A FAIL is",
          "the more interesting answer: it would mean the phi_E unit-root boundary is",
          "structural in the gear parameterization and no sampler setting reaches it."))
  V[[length(V)+1]] <- V1("S3", "cross-track agreement on the shared turnover",
    sprintf("gear tau_bar %s [%s, %s] vs pooled S2 %s (2026-08-29 pooled 2.597)",
            fmt(g1("S3","boat_tau_bar"),3), fmt(g1("S3","boat_tau_lo95"),3),
            fmt(g1("S3","boat_tau_hi95"),3), fmt(g1("S2","boat_tau_bar"),3)),
    "intervals overlap and both exclude the 1.200 prior centre",
    if (isTRUE(g1("S3","boat_tau_lo95") > 1.2) && isTRUE(g1("S2","boat_tau_lo95") > 1.2) &&
        isTRUE(g1("S3","boat_tau_lo95") <= g1("S2","boat_tau_hi95")) &&
        isTRUE(g1("S2","boat_tau_lo95") <= g1("S3","boat_tau_hi95"))) "PASS" else "REVIEW",
    paste("Two tracks with different CPUE structure fitting the same turnover on the same",
          "data. The rejected 2026-08-29 gear fit already put tau_bar at 2.552, close to",
          "the pooled 2.597, but a rejected fit is not evidence. If a converged gear fit",
          "lands in the same place, the shared turnover is a property of the data rather",
          "than of the pooled parameterization, which is what would justify adopting it."))
  V[[length(V)+1]] <- V1("S3", "the other three gear fits are untouched by the sampler change",
    (fit_exactness(basename(find_outdir("gear_resolved", STAGE_DEFS$S3$tag) %||% NA_character_),
                   .here("05_output", REF$Egear$dir), "shore",
                   "gear shore fits vs the 2026-08-29 gear run"))$observed,
    "shore fits bit-identical",
    (fit_exactness(basename(find_outdir("gear_resolved", STAGE_DEFS$S3$tag) %||% NA_character_),
                   .here("05_output", REF$Egear$dir), "shore"))$verdict,
    paste("The gear driver's per-fit branch is: shore & all_gear -> the *_shore_allgear",
          "settings; pot_closure -> the *_ringnet settings; everything else -> the track",
          "defaults. Only private_boat all-gear falls through to the defaults, so raising",
          "them should touch exactly one fit. This checks that claim instead of trusting it.",
          "A FAIL means the override reached further than intended and S3's comparison with",
          "the 2026-08-29 gear run is not clean."))
}

# ---- S4a / S4b, the 2x2
if (!is.null(S) && any(S$stage == "S4a")) {
  V[[length(V)+1]] <- V1("S4a", "a MATCHED (tau off, boat daily) cell, which stage F is not",
    sprintf("boat all-gear %s at %s AR (stage D monthly %s; stage F escalated %s); shore pot closure %s vs stage D %s",
            fmt(g1("S4a","boat_ag_BSS"),0), g1("S4a","boat_ag_ar"), fmt(REF$D$boat_ag,0),
            fmt(REF$F$boat_ag,0), fmt(g1("S4a","shore_pc_BSS"),0), fmt(REF$D$shore_pc,0)),
    "differs from stage D in the boat all-gear AR and nothing else",
    if (isTRUE(abs(g1("S4a","shore_pc_BSS") - REF$D$shore_pc) < 1) &&
        isTRUE(abs(g1("S4a","shore_ag_BSS") - REF$D$shore_ag) < 1)) "PASS" else "REVIEW",
    paste("The plan called the 2x2 'one new run' on the grounds that stage F supplies the",
          "(off, daily) cell. It does not. Stage F reached daily through ar_escalate, which",
          "moved ALL FOUR fits (shore pot closure biweekly->daily, boat pot closure",
          "monthly->daily, boat all-gear monthly->daily) through a different code path with",
          "a different number of fits per component; its shore pot closure came back 6,352",
          "against stage D's 6,331. Two of its four components therefore differ from the",
          "reference for reasons unrelated to the boat's AR. This cell forces ONE component",
          "and leaves the other three at their selected resolutions. If the two shore",
          "components do not come back at stage D's values, ar_force is leaking again and",
          "the 2x2 is contaminated in the same way stage C was."))
}
if (!is.null(S) && any(S$stage == "S4b")) {
  b_off_m <- REF$D$boat_ag; b_on_m <- g1("S2","boat_ag_BSS") %||% REF$E$boat_ag
  b_off_d <- g1("S4a","boat_ag_BSS"); b_on_d <- g1("S4b","boat_ag_BSS")
  add <- if (all(is.finite(c(b_off_m, b_on_m, b_off_d, b_on_d))))
    (b_on_d - b_off_d) - (b_on_m - b_off_m) else NA_real_
  V[[length(V)+1]] <- V1("S4b", "do the shared turnover and the daily AR compound, or overlap?",
    sprintf("2x2 boat all-gear: off/monthly %s, on/monthly %s, off/daily %s, on/daily %s; interaction %s",
            fmt(b_off_m,0), fmt(b_on_m,0), fmt(b_off_d,0), fmt(b_on_d,0), fmt(add,0)),
    "an interaction near zero means the two act on different things; a large positive one means both absorb",
    if (!is.finite(add)) "REVIEW" else if (abs(add) < 0.15 * b_off_m) "ADDITIVE" else "INTERACTING",
    paste("This is the run Stage 5 exists for. Stage E (shared tau at monthly) and stage F",
          "(per-day tau at daily) each raise the boat, by +20% and +44%, and had never been",
          "varied jointly; the gap between them is roughly 6,000 crab. THE WORKING",
          "HYPOTHESIS IS THAT E IS IDENTIFICATION AND F IS ABSORPTION, and the adequacy",
          "columns are how to tell: stage E improved the stream it should (trailer elpd",
          "+14.67, p_loo unchanged at 8.0% of n_obs, PIT bias 0.079 -> 0.021) while stage F",
          "bought a larger elpd gain by spending ~27 more effective parameters (p_loo",
          "8.0% -> 16.3%, Pareto k>0.7 1 -> 2), made the CATCH stream's elpd WORSE, and",
          "collapsed r_OSP to 10.61 with a 95% upper bound of 1,405 and sigma_r_OSP at",
          "n_eff 307. If on/daily lands near on/monthly, the shared level is doing the work",
          "and the daily AR is redundant. If they compound toward ~45,000, neither is",
          "identifying and both are absorbing observation error into the latent process."))
  V[[length(V)+1]] <- V1("S4b", "which cell is ADEQUATE, as distinct from which one SAMPLED",
    if (is.null(TAB)) "2x2 table not assembled" else
      paste(sprintf("%s: p_loo %s%% (%s), k>0.7 %s, disp n_eff %s, PIT trailer %s / osp %s",
                    TAB$cell, fmt(100*TAB$p_loo_frac,1), TAB$p_loo_stream, TAB$pareto_bad,
                    fmt(TAB$disp_neff,0), fmt(TAB$pit_trailer,3), fmt(TAB$pit_osp,3)),
            collapse = " | "),
    "reported, never gated (plan 5.5)",
    "INFO",
    paste("All four cells are expected to pass the convergence gate; six configurations",
          "already did while spanning 44% on this component. These columns are the ones",
          "that separate them, and they deliberately do NOT feed the gate: turning any of",
          "them into a hard gate is a modelling decision needing its own justification and",
          "its own run, not a default someone inherits. Note the dispersion n_eff column in",
          "particular -- the gate checks n_eff on C_sum and E_sum only, so sigma_r_OSP sat",
          "at 307 in stage F, below the gate's own floor of 400, without the gate noticing."))
}

# ---- S5
if (!is.null(S) && any(S$stage == "S5")) {
  pc <- g1("S5","boat_pc_BSS"); ag <- g1("S5","boat_ag_BSS")
  V[[length(V)+1]] <- V1("S5", "boat pot closure at biweekly reconciles the two tracks, uncontaminated",
    sprintf("pot closure %s at %s (pooled monthly %s, gear biweekly %s); all-gear %s at %s (must be rung 4's %s)",
            fmt(pc,0), g1("S5","boat_pc_ar"), fmt(REF$D$boat_pc,0), fmt(REF$gear5$boat_pc,0),
            fmt(ag,0), g1("S5","boat_ag_ar"), fmt(REF$D$boat_ag,0)),
    "pot closure within 5% of the gear track's 743 AND all-gear unmoved from 25,868",
    if (isTRUE(abs(pc - REF$gear5$boat_pc)/REF$gear5$boat_pc < 0.05) &&
        isTRUE(abs(ag - REF$D$boat_ag) < 1)) "PASS" else "REVIEW",
    paste("The 2026-08-28 stage C answered the reconciliation question (735 at biweekly",
          "against the gear track's 743, a 1.1% gap where the pooled monthly fit gave 850)",
          "but its ar_force was honoured per POPULATION, so as.character() on a one-element",
          "list silently forced the boat's ALL-GEAR sub-season to biweekly as well and moved",
          "that component 25,883 -> 28,893 -- a side effect worth 3,000 crab that was",
          "reported as a pot-closure result. The nested form is now resolved per sub-season",
          "and errors on an unnamed list. THE SECOND HALF OF THIS CRITERION IS THE POINT:",
          "the all-gear component must come back at rung 4's value. If it has moved, the",
          "leak is not fixed and stage C's port total stays unusable."))
}

VD <- if (length(V)) do.call(rbind, V) else NULL
if (!is.null(VD)) utils::write.csv(VD, ver_path, row.names = FALSE)

banner("BATCH SUMMARY")
if (!is.null(S)) {
  print(S[, c("stage","plan_item","model","ok","minutes","port_BSS_catch","boat_ag_BSS",
              "boat_ag_ar","boat_tau_bar","boat_ag_method","fits_passed")], row.names = FALSE)
  cat("\n  full summary -> ", sum_path, "\n", sep = "")
}
if (!is.null(TAB)) {
  rule(); cat("  TAU x AR, private-boat all-gear\n\n")
  cat(sprintf("  %-14s %-8s %10s %8s %9s %9s %7s\n",
              "cell", "reported", "catch", "tau_bar", "p_loo", "k>0.7", "r_OSP"))
  for (i in seq_len(nrow(TAB)))
    cat(sprintf("  %-14s %-8s %10s %8s %8s%% %9s %7s\n",
                substr(TAB$cell[i], 1, 14),
                if (is.na(TAB$reported[i])) "?" else if (isTRUE(TAB$reported[i])) "BSS" else "PE",
                fmt(TAB$boat_all_gear[i], 0), fmt(TAB$tau_bar[i], 3),
                fmt(100 * TAB$p_loo_frac[i], 1), fmt(TAB$pareto_bad[i], 0),
                fmt(TAB$r_OSP[i], 2)))
  cat("\n  2x2 -> ", tab_path, "\n", sep = "")
}
if (!is.null(VD)) {
  rule()
  for (i in seq_len(nrow(VD)))
    cat(sprintf("[%s] %-11s %s\n        %s\n", VD$stage[i], VD$verdict[i],
                VD$criterion[i], VD$observed[i]))
  cat("\n  verdicts -> ", ver_path, "\n", sep = "")
}
banner(sprintf("DONE  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("  Next: the S1 fill decision (5.1), then read the 2x2 and settle 5.4 -- whether\n")
cat("  ar_max_resolution encodes a live pathology or a preference.\n\n")
}
