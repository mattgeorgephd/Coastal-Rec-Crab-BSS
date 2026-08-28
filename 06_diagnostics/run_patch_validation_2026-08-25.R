#!/usr/bin/env Rscript
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
# run_patch_validation_2026-08-25.R
#
# The validation ladder for the 2026-08-25 improvement batch, run end to end with
# ONE change per rung and a pre-set pass criterion for each. Follows the house
# pattern of run_osp_validation.R (job matrix, resume, incremental summary CSV) and
# adds the part that batch lacked: the criteria are written down BEFORE the runs, in
# code, so each rung is judged rather than interpreted afterwards.
#
#   source("06_diagnostics/run_patch_validation_2026-08-25.R")        # RStudio: Source
#   Rscript 06_diagnostics/run_patch_validation_2026-08-25.R          # terminal
#
# START HERE:  set DRY_RUN <- TRUE and source it. That prints the fully resolved
# config for every rung and every criterion, and fits nothing. Read that before
# committing a day of compute.
#
# ============================================================================
# THE LADDER
# ============================================================================
#
#   1  baseline    every 2026-08-25 feature set to its PRE-PATCH value.
#                  pooled AND gear-resolved.
#   2  ie-unit     + ie_shore_obs_unit = "auto"      (the shore I/E unit fix)
#   3  weekend     + days_wkend = Sat/Sun
#   4  minint      + bss_min_interviews = 15         (rung 4 == the shipped config)
#   5  gear-ship   the shipped config on the gear-resolved track (cross-check)
#
# MODE = "cumulative" (default) means each rung adds one change to the previous one,
# so rung 4 IS the configuration you would publish and each rung's delta is read
# against the rung below it. MODE = "isolated" instead applies each change to the
# baseline alone. Cumulative is the default for one practical reason: it ends on a
# configuration you can actually adopt. Its one assumption is that the changes do not
# interact, and there IS a plausible interaction to watch here -- rungs 3 and 4 both
# act on the boat pot closure (the weekend change alters its PE strata, the interview
# floor decides whether it is fitted at all). If rung 4 surprises you, re-run rung 4
# in isolated mode before drawing a conclusion.
#
# ============================================================================
# WHAT THE BASELINE CAN AND CANNOT REVERT (read before judging rung 1)
# ============================================================================
#
# Revertible, and reverted by RUNG 1:
#   ie_shore_obs_unit, days_wkend, bss_min_interviews, bss_min_interviews_fitted,
#   ar_escalate, opener_covariate_mode, use_osp_crab_lower, pe_empty_effort_stratum,
#   crab_fraction_restrict_to_fit.
#
# NOT revertible, because they are corrections rather than toggles:
#   * the shore monthly PE effort share moved off the crabber-hours formula. Affects
#     only the monthly SPLIT of a component that reports PE (here: the boat pot
#     closure). Totals are untouched -- the share is normalised.
#   * the fishery-opener spillover diagnostic's CPUE denominator, the gear-resolved
#     daily-combined series, plot labels, the bss_L_effective column names, and the
#     extra rows in structural_params_*. All output-only.
#   * the retired `real B3` in the POOLED Stan model. Read the next paragraph.
#
# WHY RUNG 1 IS JUDGED DIFFERENTLY ON THE TWO TRACKS.
#   GEAR-RESOLVED is expected to be BIT-IDENTICAL. Every parameter this patch adds to
#   that model is zero-size at the default settings (K_open = 0, osp_crab_lower = 0),
#   and f_crab_param -> f_theta is a rename at the same position, so the unconstrained
#   parameter vector is unchanged and the fixed-seed RNG stream is preserved. Its port
#   total should match the reference to the last digit. That is the sharpest test in
#   this whole ladder and it is why the gear track is run first.
#   POOLED cannot be. The retired `real B3` was declared unconditionally with a proper
#   prior and, with razor off, entered no likelihood -- a genuine sampled dimension
#   that is now gone. Removing it shifts the HMC trajectories, so a fixed-seed rerun
#   will not diff to zero even though the posterior is unchanged. The pooled rung is
#   therefore judged in MONTE CARLO STANDARD ERRORS of the reference median, which the
#   script computes from that run's own interval width and effective sample size.
#
# ============================================================================
# COST
# ============================================================================
# Six model-runs. A pooled run is roughly 3-6 h on 4 cores and a gear-resolved run
# 2-4 h, so budget 16-32 h wall clock, plus a Stan recompile on the first run of each
# model (both .stan files changed). The script is RESUMABLE: a rung whose output
# folder already contains port_total_Dungeness_Kept.csv is skipped and its numbers are
# still extracted, so an interrupted batch can simply be re-sourced.
#
# ============================================================================
# OUTPUT
# ============================================================================
#   05_output/patch_validation_2026-08-25_summary.csv   one row per rung x model,
#       appended as each finishes, so a batch killed at hour 20 still leaves evidence.
#   05_output/patch_validation_2026-08-25_verdicts.csv  one row per criterion, with
#       the observed value, the threshold, and PASS / REVIEW / FAIL.
#   a console report at the end.
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- FALSE          # TRUE: resolve and print everything, fit nothing. START HERE.
MODE    <- "cumulative"   # "cumulative" (rung N = rung N-1 + one change) | "isolated"
RUNGS   <- 1:5            # subset to re-run part of the ladder, e.g. 2:4
RESUME  <- TRUE           # skip a rung whose output folder already looks complete

# =========================================================================== #

# Repo root, resolved from wherever this was launched. Deliberately NOT via here::here,
# so a DRY RUN needs no packages at all: you should be able to inspect the ladder on a
# machine with no Stan toolchain before committing a day of compute to it.
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) &&
    dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions")))
  stop("Run this from the repository root (or from 06_diagnostics/): 03_R_functions not found.")
.here <- function(...) file.path(.root, ...)

# The fitting dependencies are loaded ONLY for a real run.
if (!isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages({ library(here); library(rmarkdown) })
  load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
  install.lib <- load.lib[!load.lib %in% installed.packages()]
  for (lib in install.lib) install.packages(lib, dependencies = TRUE)
  invisible(sapply(load.lib, require, character.only = TRUE))
  rstan_options(auto_write = TRUE)
  invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), source))
}

source(.here("run_config.R"))          # defines run_config (the SHIPPED config)
BASE <- run_config

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
banner <- function(msg) cat("\n", strrep("=", 78), "\n ", msg, "\n", strrep("=", 78), "\n", sep = "")
rule   <- function() cat(strrep("-", 78), "\n")

model_rmd <- list(
  pooled        = .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd"))
prefix <- list(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")
stopifnot(all(file.exists(unlist(model_rmd))))

# ---------------------------------------------------------------------------
# PRE-FLIGHT: the R-to-Stan data contract.
#
# Added 2026-08-26 after the first attempt at this ladder lost all six model-runs to a
# defect this check finds in milliseconds. Five variables were declared in both .stan
# data blocks and never passed by the prep functions; Stan failed at data initialization
# on every fit, rstan returned an empty stanfit instead of raising, and the run died three
# layers downstream with "dim(X) must have a positive length". Nothing about that message
# named the cause, and the fits are hours apart, so the cost of NOT running this check is
# the whole batch.
#
# Deliberately runs in the DRY RUN too, and needs no packages: it is a text parse of the
# .stan data block against the prep source. bss_assert_stan_data() repeats the check
# exactly, per fit, immediately before each sampler call.
# ---------------------------------------------------------------------------
local({
  src_fn <- .here("03_R_functions", "bss_stan_fit.R")
  if (!file.exists(src_fn)) {
    warning("bss_stan_fit.R not found; skipping the Stan data-contract pre-flight.", call. = FALSE)
    return(invisible(NULL))
  }
  e <- new.env(); sys.source(src_fn, envir = e)
  pairs <- list(
    pooled        = c("crab_bss_pooled.stan",        "prep_bss_crab_pooled.R"),
    gear_resolved = c("crab_bss_gear_resolved.stan", "prep_bss_crab_gear.R"))
  fail <- character(0)
  for (m in names(pairs)) {
    stan_f <- .here("02_stan_models", pairs[[m]][1])
    prep_f <- .here("03_R_functions", pairs[[m]][2])
    if (!file.exists(stan_f) || !file.exists(prep_f)) next
    need <- e$bss_stan_data_names(stan_f)
    src  <- paste(readLines(prep_f, warn = FALSE), collapse = "\n")
    miss <- need[!vapply(need, function(v)
      grepl(paste0("(^|[^A-Za-z0-9_.])", v, "([^A-Za-z0-9_]|$)"), src, perl = TRUE),
      logical(1))]
    cat(sprintf("  pre-flight  %-14s %2d Stan data variables, %d not built in %s\n",
                m, length(need), length(miss), basename(prep_f)))
    if (length(miss))
      fail <- c(fail, sprintf("%s: %s (not built in %s)", m,
                              paste(miss, collapse = ", "), basename(prep_f)))
  }
  if (length(fail))
    stop("PRE-FLIGHT FAILED, the R-to-Stan data contract is broken. Every fit would ",
         "return an empty stanfit.\n  ", paste(fail, collapse = "\n  "),
         "\n  Fix the prep function(s) before starting the ladder; see the 2026-08-26 ",
         "entry in development_notes/PIPELINE_STATUS.md.", call. = FALSE)
})

# ---------------------------------------------------------------------------
# The pre-patch configuration. Every key here is a 2026-08-25 addition or change,
# set back to what the code did before the patch. NOTE bss_min_interviews_fitted = 0:
# the post-filter floor did not exist before, and leaving it NULL would make it inherit
# bss_min_interviews, which is NEW behaviour and would spoil the baseline.
# ---------------------------------------------------------------------------
PRE_PATCH <- list(
  ie_shore_obs_unit             = "crabber_hours",
  days_wkend                    = c("Friday", "Saturday", "Sunday"),
  bss_min_interviews            = 20,
  bss_min_interviews_fitted     = 0,
  ar_escalate                   = FALSE,
  opener_covariate_mode         = "off",
  razor_dig_mode                = "no",
  use_osp_crab_lower            = FALSE,
  pe_empty_effort_stratum       = "zero",
  crab_fraction_restrict_to_fit = FALSE
)

# ---------------------------------------------------------------------------
# Reference numbers, from the two runs this ladder is validated against. Hard-coded
# on purpose: a criterion you can edit after seeing the answer is not a criterion.
# ---------------------------------------------------------------------------
REF <- list(
  pooled = list(
    dir            = "20260804/pooled-CPUE-boat-count-validation-run",
    port_catch     = 67312, port_effort   = 39145,
    port_lo95      = 50601, port_hi95     = 93461,
    shore_ag_catch = 20665, shore_pc_catch = 6263,
    boat_ag_catch  = 27668, boat_pc_pe     = 351,
    shore_ag_sigma_IE = 1.0574,      # THE number the I/E unit fix should move
    shore_ag_B1       = 0.5087),     # weekend effort effect, log scale (1.66x)
  gear_resolved = list(
    dir        = "20260805/gear-type-CPUE-model-boat-count-validation-run",
    port_catch = 66461, port_effort = 38649)
)

# ---------------------------------------------------------------------------
# THE LADDER. `delta` is the one change this rung introduces.
# ---------------------------------------------------------------------------
rung <- function(id, tag, models, delta, headline) {
  list(id = id, tag = tag, models = models, delta = delta, headline = headline)
}
LADDER <- list(
  rung(1, "PV1-baseline", c("gear_resolved", "pooled"), list(),
       "pre-patch equivalent; gear-resolved must reproduce EXACTLY"),
  rung(2, "PV2-ie-unit",  "pooled", list(ie_shore_obs_unit = "auto"),
       "shore I/E observation follows the effort unit (crabber trips, not hours)"),
  rung(3, "PV3-weekend",  "pooled", list(days_wkend = c("Saturday", "Sunday")),
       "Friday moves to the weekday stratum"),
  rung(4, "PV4-minint",   "pooled", list(bss_min_interviews = 15,
                                          bss_min_interviews_fitted = NULL),
       "boat pot closure becomes eligible for a BSS fit; == the shipped config"),
  rung(5, "PV5-gear-ship", "gear_resolved",
       list(ie_shore_obs_unit = "auto", days_wkend = c("Saturday", "Sunday"),
            bss_min_interviews = 15, bss_min_interviews_fitted = NULL),
       "shipped config on the gear track; re-establishes the pooled/gear cross-check")
)

# Resolve a rung's full config under the selected MODE.
resolve_cfg <- function(i) {
  cfg <- modifyList(BASE, PRE_PATCH)
  if (identical(MODE, "cumulative")) {
    for (k in seq_len(i)) cfg <- modifyList(cfg, LADDER[[k]]$delta, keep.null = TRUE)
  } else {
    cfg <- modifyList(cfg, LADDER[[i]]$delta, keep.null = TRUE)
  }
  cfg
}

# ---------------------------------------------------------------------------
# Extraction. Every reader is guarded: a missing or renamed file yields NA rather
# than aborting a batch that has already spent hours.
# ---------------------------------------------------------------------------
find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}
rd <- function(outdir, f) {
  if (is.na(outdir)) return(NULL)
  tryCatch(read.csv(file.path(outdir, f), stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL, warning = function(w) NULL)
}
num1 <- function(x) { x <- suppressWarnings(as.numeric(x)); if (length(x) >= 1) x[1] else NA_real_ }
sp_get <- function(outdir, fit, par, col = "mean") {
  d <- rd(outdir, sprintf("structural_params_%s_Dungeness_Kept.csv", fit))
  if (is.null(d) || !"parameter" %in% names(d)) return(NA_real_)
  num1(d[[col]][d$parameter == par])
}

extract_run <- function(rung_id, model, run_tag, outdir, minutes, ok) {
  row <- data.frame(
    rung = rung_id, model = model, run_tag = run_tag, ok = ok, minutes = minutes,
    outdir = if (is.na(outdir)) "" else basename(outdir),
    port_BSS_catch = NA_real_, port_BSS_lo95 = NA_real_, port_BSS_hi95 = NA_real_,
    port_BSS_effort = NA_real_, port_PE_catch = NA_real_,
    shore_ag_BSS_catch = NA_real_, shore_pc_BSS_catch = NA_real_,
    boat_ag_BSS_catch = NA_real_, boat_ag_BSS_effort = NA_real_,
    boat_pc_method = NA_character_, boat_pc_PE_catch = NA_real_,
    shore_ag_sigma_IE = NA_real_, shore_ag_B1 = NA_real_, shore_ag_R_G = NA_real_,
    shore_ag_L_med = NA_real_, shore_ag_L_width = NA_real_,
    shore_ag_neff = NA_real_, fits_passed = NA_character_, ar_resolutions = NA_character_,
    # 2026-08-27: METHOD PER COMPONENT. Rung 2 of the first ladder reported a shore all-gear
    # BSS catch of 21,017 in this table from a fit the gate had REJECTED (2,957 divergences,
    # E_sum n_eff = 30, method "PE (convergence fail)"), and nothing in the row said so. Any
    # number here that came from a rejected fit has to arrive carrying that label.
    shore_ag_method = NA_character_, shore_pc_method = NA_character_,
    boat_ag_method = NA_character_,
    empty_effort_note = NA_character_, stringsAsFactors = FALSE)
  if (is.na(outdir) || !dir.exists(outdir)) return(row)

  pt <- rd(outdir, "port_total_Dungeness_Kept.csv")
  if (!is.null(pt) && "Estimate" %in% names(pt)) {
    is_catch  <- pt$Estimate %in% c("Expected_Catch", "Catch")
    is_effort <- grepl("^Effort", pt$Estimate)
    row$port_BSS_catch  <- num1(pt$BSS_median[is_catch])
    row$port_BSS_lo95   <- num1(pt$BSS_lo95[is_catch])
    row$port_BSS_hi95   <- num1(pt$BSS_hi95[is_catch])
    row$port_PE_catch   <- num1(pt$PE[is_catch])
    row$port_BSS_effort <- num1(pt$BSS_median[is_effort])
  }
  pv <- rd(outdir, "pe_vs_bss_comparison.csv")
  if (!is.null(pv) && "component" %in% names(pv)) {
    g <- function(pat, col) num1(pv[[col]][grepl(pat, pv$component)])
    row$shore_ag_BSS_catch <- g("^shore \\(All gear\\)",        "BSS_catch")
    row$shore_pc_BSS_catch <- g("^shore \\(Pot closure\\)",     "BSS_catch")
    row$boat_ag_BSS_catch  <- g("^private_boat \\(All gear\\)", "BSS_catch")
    row$boat_ag_BSS_effort <- g("^private_boat \\(All gear\\)", "BSS_effort")
    row$boat_pc_PE_catch   <- g("^private_boat \\(Pot closure\\)", "PE_catch")
  }
  cr <- rd(outdir, "convergence_report.csv")
  if (!is.null(cr) && "fit" %in% names(cr)) {
    meth <- function(pat) {
      i <- grepl(pat, cr$fit)
      if (any(i)) as.character(cr$method_selected[i][1]) else NA_character_
    }
    row$boat_pc_method  <- meth("^private_boat_ring_net_only")
    row$shore_ag_method <- meth("^shore_all_gear")
    row$shore_pc_method <- meth("^shore_ring_net_only")
    row$boat_ag_method  <- meth("^private_boat_all_gear")
    j <- grepl("^shore_all_gear", cr$fit)
    if (any(j)) row$shore_ag_neff <- num1(cr$C_sum_neff[j])
  }
  # 2026-08-27: the empty-effort-stratum audit, which run_pe_*() now writes to a CSV instead
  # of only cat()-ing it (the pooled PE chunk is results='hide', so on that track the report
  # reached nothing). The weekend redefinition changes which strata are empty, and an empty
  # stratum is expanded at ZERO effort under the shipped fill, so this belongs in the summary.
  ee <- rd(outdir, "pe_empty_effort_strata.csv")
  if (!is.null(ee) && "component" %in% names(ee) && nrow(ee) > 0) {
    row$empty_effort_note <- paste(sprintf("%s %s/%s d (%.0f%%)",
      ee$component, ee$n_empty_days, ee$n_calendar_days,
      100 * suppressWarnings(as.numeric(ee$empty_day_fraction))), collapse = "; ")
  }
  row$shore_ag_sigma_IE <- sp_get(outdir, "shore_all_gear", "sigma_IE")
  row$shore_ag_B1       <- sp_get(outdir, "shore_all_gear", "B1")
  row$shore_ag_R_G      <- sp_get(outdir, "shore_all_gear", "R_G")

  # L (the daily expansion factor = tau_shore in production). The column that carries
  # the prior centre was renamed on 2026-08-25 (L_prior_mu -> L_prior_center), so accept
  # either; the posterior columns are unchanged.
  Ld <- rd(outdir, "bss_L_effective_shore_all_gear_Dungeness_Kept.csv")
  if (!is.null(Ld) && "L_posterior_median" %in% names(Ld)) {
    row$shore_ag_L_med   <- num1(stats::median(suppressWarnings(as.numeric(Ld$L_posterior_median)), na.rm = TRUE))
    row$shore_ag_L_width <- num1(stats::median(
      suppressWarnings(as.numeric(Ld$L_posterior_hi95) - as.numeric(Ld$L_posterior_lo95)), na.rm = TRUE))
  }
  ss <- rd(outdir, "season_summary.csv")
  if (!is.null(ss) && "metric" %in% names(ss)) {
    row$fits_passed    <- as.character(ss$value[ss$metric == "BSS fits passed"])[1]
    row$ar_resolutions <- as.character(ss$value[ss$metric == "AR temporal resolution"])[1]
  }
  row
}

# ---------------------------------------------------------------------------
# CRITERIA. Written before the runs. Each returns rows of
# (rung, criterion, observed, threshold, verdict, why).
# PASS = the pre-set expectation held. REVIEW = look at it before going on.
# FAIL = the expectation did NOT hold, which is a finding, not necessarily a bug.
# ---------------------------------------------------------------------------
V <- function(rung, criterion, observed, threshold, verdict, why)
  data.frame(rung = rung, criterion = criterion, observed = observed,
             threshold = threshold, verdict = verdict, why = why, stringsAsFactors = FALSE)
fmt <- function(x, d = 1) if (is.na(x)) "NA" else formatC(x, format = "f", digits = d, big.mark = ",")

# Monte Carlo SE of a posterior median, from the run's own 95% interval and n_eff.
# sd ~ (hi - lo) / 3.92; MCSE(median) ~ 1.253 * sd / sqrt(n_eff).
mcse_median <- function(lo, hi, neff) {
  if (any(is.na(c(lo, hi, neff))) || neff <= 0) return(NA_real_)
  1.253 * ((hi - lo) / 3.92) / sqrt(neff)
}

# ---------------------------------------------------------------------------
# Bit-identity of the gear-resolved FITS against the reference folder (2026-08-27).
#
# What must match byte for byte: bss_summary_*.csv and convergence_report.csv. Both are pure
# functions of the fitted draws, so any difference is a behaviour change.
#
# What must NOT be required to match: bss_full_summary_*.csv gains rows for quantities the
# 2026-08-25 batch newly REPORTS (f_lower, f_lower_out, osp_f_kappa_out) and renames one
# (f_crab_param -> f_theta), so it is compared on shared rows only; and every file built from
# permuted or subsampled draws (bss_draws_summed_*, ppc_*, loo_*, port_total, monthly_*)
# shuffles even when the fits are identical, so it is excluded entirely.
#
# Returns list(observed, verdict) for V().
gear_fit_exactness <- function(new_dirname, ref_dir) {
  if (is.na(new_dirname) || !nzchar(new_dirname))
    return(list(observed = "no output folder found", verdict = "REVIEW"))
  # Resolve the basename back to a path the same way find_outdir() does: NEWEST first. A
  # failed earlier attempt at the same rung leaves a same-named folder under an earlier date
  # (the 2026-08-25 scratch folders are exactly that), and picking the first match instead of
  # the newest compares this run against an abandoned one.
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == new_dirname]
  new_dir <- if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NULL
  if (is.null(new_dir) || !dir.exists(new_dir) || !dir.exists(ref_dir))
    return(list(observed = "run or reference folder missing", verdict = "REVIEW"))

  # convergence_report.csv is the one file whose SCHEMA did not change, so it must match byte
  # for byte. bss_summary_*.csv must NOT: the 2026-08-25 batch adds reported quantities
  # (f_crab_out on the boat, f_lower_out / osp_f_kappa_out everywhere), so those files
  # legitimately gain rows and are compared row-wise below alongside the full summaries.
  exact_files <- "convergence_report.csv"
  exact_files <- exact_files[file.exists(file.path(ref_dir, exact_files)) &
                             file.exists(file.path(new_dir, exact_files))]
  if (!length(exact_files))
    return(list(observed = "no convergence report to compare", verdict = "REVIEW"))

  same <- vapply(exact_files, function(f)
    identical(readBin(file.path(ref_dir, f), "raw", file.size(file.path(ref_dir, f))),
              readBin(file.path(new_dir, f), "raw", file.size(file.path(new_dir, f)))),
    logical(1))

  # Shared-row comparison of every per-fit parameter summary, at full double precision.
  fs <- c(list.files(new_dir, pattern = "^bss_summary_.*\\.csv$"),
          list.files(new_dir, pattern = "^bss_full_summary_.*\\.csv$"))
  fs <- fs[file.exists(file.path(ref_dir, fs))]
  rows_ok <- TRUE; n_rows <- 0L
  for (f in fs) {
    a <- tryCatch(utils::read.csv(file.path(ref_dir, f), row.names = 1, check.names = FALSE),
                  error = function(e) NULL)
    b <- tryCatch(utils::read.csv(file.path(new_dir, f), row.names = 1, check.names = FALSE),
                  error = function(e) NULL)
    if (is.null(a) || is.null(b)) { rows_ok <- FALSE; next }
    common <- intersect(rownames(a), rownames(b))
    cols   <- intersect(names(a), names(b))
    n_rows <- n_rows + length(common)
    if (!isTRUE(all.equal(a[common, cols], b[common, cols], tolerance = 0))) rows_ok <- FALSE
  }

  obs <- sprintf("convergence_report byte-identical: %s; %d shared parameter rows across %d summaries %s",
                 if (all(same)) "yes" else "NO", n_rows, length(fs),
                 if (rows_ok) "identical at full precision" else "DIFFER")
  list(observed = obs, verdict = if (all(same) && rows_ok) "PASS" else "FAIL")
}

judge <- function(S) {
  get1 <- function(r, m, col) {
    x <- S[[col]][S$rung == r & S$model == m]
    if (length(x)) x[1] else NA
  }
  out <- list()

  # ---- RUNG 1: baseline reproduction --------------------------------------
  gp <- get1(1, "gear_resolved", "port_BSS_catch")
  if (!is.na(gp)) {
    d <- gp - REF$gear_resolved$port_catch
    # 2026-08-27, CRITERION REPAIRED. This used to demand delta == 0 on the port total and it
    # FAILED on the 2026-08-26 ladder at delta = -19 (0.03%). The model was bit-identical the
    # whole time; the criterion was testing the wrong object. The port total is not a function
    # of the fits alone: it is assembled from rstan::extract(permuted = TRUE) draws, and the
    # per-fit draw files come from an unseeded sample.int() in save_run_diagnostics.R, so only
    # ~490 of 2,000 draw indices were shared between two runs of the same fits. PIPELINE_STATUS
    # had already recorded that behaviour in the Run 6 note ("totals derived from them, and
    # RNG-sensitive downstream diagnostics shuffle") and this criterion contradicted it.
    #
    # Exactness is now tested where exactness exists: byte equality of the per-fit parameter
    # summaries and the convergence report against the reference folder. The port total is
    # judged on Monte Carlo error, exactly as the pooled track is.
    out[[length(out)+1]] <- V(1, "gear-resolved port total within Monte Carlo error of the reference",
      sprintf("%s (ref %s, delta %s)", fmt(gp,0), fmt(REF$gear_resolved$port_catch,0), fmt(d,0)),
      "< 0.5% of the reference",
      if (isTRUE(abs(d) < 0.005 * REF$gear_resolved$port_catch)) "PASS" else "FAIL",
      paste("The gear fits are expected to be BIT-identical (every parameter this patch adds",
            "there is zero-size at the default settings), but the port total is assembled from",
            "permuted and subsampled draws and moves by a fraction of one MCSE even when the",
            "fits do not. Exactness is tested by the next criterion; this one only catches a",
            "gross change. See the 2026-08-27 note in PIPELINE_STATUS."))

    # The real exactness test: the fitted parameters themselves.
    exact <- gear_fit_exactness(get1(1, "gear_resolved", "outdir"),
                                file.path(.here("05_output"), REF$gear_resolved$dir))
    out[[length(out)+1]] <- V(1, "gear-resolved FITS are bit-identical to the reference",
      exact$observed, "every shared file byte-identical",
      exact$verdict,
      paste("This is the sharpest test in the batch. The gear model never carried B3 and every",
            "parameter the 2026-08-25 batch adds to it is zero-size when the feature is off, so",
            "the unconstrained parameter vector and the fixed-seed RNG stream must be unchanged.",
            "bss_summary_*.csv and convergence_report.csv are pure functions of the fit, so they",
            "must match byte for byte; bss_full_summary_*.csv legitimately gains rows for the new",
            "reported quantities, so it is compared on SHARED rows only. Any difference here is a",
            "real behaviour change and must be found before any later rung is interpreted."))
  }
  pp <- get1(1, "pooled", "port_BSS_catch")
  if (!is.na(pp)) {
    se <- mcse_median(get1(1,"pooled","port_BSS_lo95"), get1(1,"pooled","port_BSS_hi95"),
                      get1(1,"pooled","shore_ag_neff"))
    d  <- pp - REF$pooled$port_catch
    n  <- if (is.na(se) || se <= 0) NA_real_ else abs(d) / se
    out[[length(out)+1]] <- V(1, "pooled port total within Monte Carlo error of the reference",
      sprintf("%s (ref %s, delta %s = %s MCSE)", fmt(pp,0), fmt(REF$pooled$port_catch,0),
              fmt(d,0), if (is.na(n)) "NA" else fmt(n,1)),
      "< 3 MCSE",
      if (is.na(n)) "REVIEW" else if (n < 3) "PASS" else if (n < 6) "REVIEW" else "FAIL",
      paste("The pooled model lost the retired `real B3` dimension, so the RNG stream shifted",
            "and a byte comparison is meaningless. MCSE is derived from this run's own interval",
            "width and n_eff. A delta beyond ~6 MCSE is a real change, not sampling noise."))
  }

  # ---- RUNG 2: the shore I/E unit fix -------------------------------------
  s2 <- get1(2, "pooled", "shore_ag_sigma_IE"); s1 <- get1(1, "pooled", "shore_ag_sigma_IE")
  if (!is.na(s2)) {
    out[[length(out)+1]] <- V(2, "shore all-gear sigma_IE falls sharply",
      sprintf("%s (rung 1: %s, reference: %s)", fmt(s2,3), fmt(s1,3), fmt(REF$pooled$shore_ag_sigma_IE,3)),
      "< 0.70",
      if (s2 < 0.70) "PASS" else if (s2 < 0.95) "REVIEW" else "FAIL",
      paste("sigma_IE ~ 1.07 was the only free parameter positioned to absorb a ~4x scale",
            "mismatch between crabber-HOURS observed and crabber-TRIPS predicted. If it does not",
            "fall, the unit mismatch was NOT the cause of GR-9 and that item stays open with a",
            "new hypothesis needed -- most likely the I/E days' representativeness, which the",
            "diagnostic already flags (mean percentile 0.27, p = 0.06)."))
  }
  b2 <- get1(2, "pooled", "boat_ag_BSS_catch"); b1 <- get1(1, "pooled", "boat_ag_BSS_catch")
  if (!is.na(b2) && !is.na(b1)) {
    out[[length(out)+1]] <- V(2, "boat is untouched by the shore I/E fix",
      sprintf("%s vs %s (delta %s)", fmt(b2,0), fmt(b1,0), fmt(b2-b1,0)), "delta == 0",
      if (abs(b2-b1) < 0.5) "PASS" else "FAIL",
      paste("The change is scoped to the SHORE observation column and the boat fit carries no",
            "I/E observations at all (IE_n = 0), so the boat must be identical. Any movement",
            "means the change leaked outside its scope."))
  }
  Lw2 <- get1(2, "pooled", "shore_ag_L_width"); Lw1 <- get1(1, "pooled", "shore_ag_L_width")
  if (!is.na(Lw2) && !is.na(Lw1) && Lw1 > 0) {
    r <- Lw2 / Lw1
    # 2026-08-27, CRITERION RETIRED TO INFO. This expected the fixed I/E stream to inform
    # tau_shore, and returned REVIEW at ratio 0.93 on the 2026-08-26 ladder. The prediction was
    # not merely unmet, it was structurally impossible: crab_bss_pooled.stan declares
    # `vector[D * estimate_L] L_raw` with L[d] = L_data[d] * exp(L_prior_sigma[d] * L_raw[d]),
    # so there is NO shared tau_shore parameter -- L is D independent per-day draws, each
    # anchored on the same prior. Four in-window I/E days can inform four days out of 289, and
    # the ratio reported here is a MEDIAN over the 285 days that carry no likelihood at all.
    # On the four days that do carry data the converged run gives 0.74, 0.69, 1.22, 1.15.
    #
    # It stays as INFO rather than being deleted: the per-day widths are still worth recording,
    # and this line is where the number goes when a shared-tau model makes the test meaningful.
    # See the "shared tau" item in development_notes/improvement-plan-2026-08-27.md.
    out[[length(out)+1]] <- V(2, "shore L posterior width (median over days; NOT a tau test)",
      sprintf("95%% width %s vs %s (ratio %s)", fmt(Lw2,3), fmt(Lw1,3), fmt(r,2)), "informational",
      "INFO",
      paste("L is a per-day vector of independent draws, not a shared tau_shore, so a handful of",
            "I/E days cannot narrow a season-level turnover and this median is dominated by the",
            "days with no observation. Recorded, not gated. Restore this as a PASS/FAIL criterion",
            "only once the model carries a shared tau with per-day deviations."))
  }
  # 2026-08-27, NEW. The first ladder's rung 2 lost its shore all-gear fit to the gate and no
  # criterion noticed: the summary reported the rejected fit's BSS catch, and the port total
  # silently carried the PE instead, which was the whole of the 67,323 -> 62,069 drop. Every
  # rung that is supposed to leave a component fitted must now say so explicitly.
  for (.r in 2:4) {
    m_ag <- get1(.r, "pooled", "shore_ag_method")
    if (!is.na(m_ag))
      out[[length(out)+1]] <- V(.r, "shore all-gear still reports a BSS fit",
        m_ag, "BSS",
        if (identical(m_ag, "BSS")) "PASS" else "FAIL",
        paste("A rung that changes an observation unit, a day-type definition or an interview",
              "floor must not silently cost a component its fit. When this FAILs, every other",
              "number for that rung is a PE substitution and must be read as such -- and the",
              "per-chain sampler_diagnostics_*.csv is the first place to look, because one",
              "stalled chain out of four looks exactly like a structural failure in the",
              "aggregate divergence count."))
  }

  s2c <- get1(2, "pooled", "shore_ag_BSS_catch"); s1c <- get1(1, "pooled", "shore_ag_BSS_catch")
  if (!is.na(s2c) && !is.na(s1c))
    out[[length(out)+1]] <- V(2, "shore all-gear catch moves (magnitude recorded, not gated)",
      sprintf("%s vs %s (%s%%)%s", fmt(s2c,0), fmt(s1c,0), fmt(100*(s2c-s1c)/s1c,1),
              if (identical(get1(2,"pooled","shore_ag_method"), "BSS")) ""
              else "  [FROM A REJECTED FIT -- not a result]"),
      "informational",
      "INFO", "This is the headline effect of the patch. The direction was not predictable in advance, because the old stream was pulling lambda_E UP against the gear counts.")

  # ---- RUNG 3: weekend = Sat/Sun ------------------------------------------
  B3 <- get1(3, "pooled", "shore_ag_B1"); B2r <- get1(2, "pooled", "shore_ag_B1")
  if (!is.na(B3) && !is.na(B2r))
    out[[length(out)+1]] <- V(3, "weekend effort effect B1 RISES once Friday leaves the stratum",
      sprintf("%s (exp %sx) vs %s (exp %sx)", fmt(B3,3), fmt(exp(B3),2), fmt(B2r,3), fmt(exp(B2r),2)),
      "B1 increases",
      if (B3 > B2r) "PASS" else "REVIEW",
      paste("Pooling Friday dragged the fitted multiplier from ~2.34x to ~1.74x on the season's",
            "own data, so removing it should push B1 back up. If it does not, the day-type",
            "signal is not what the marginal comparison suggested and the change needs a second look."))
  bp3 <- get1(3, "pooled", "boat_pc_PE_catch"); bp2 <- get1(2, "pooled", "boat_pc_PE_catch")
  if (!is.na(bp3) && !is.na(bp2))
    out[[length(out)+1]] <- V(3, "boat pot-closure PE (watch the empty-stratum side effect)",
      sprintf("%s vs %s (%s%%)", fmt(bp3,0), fmt(bp2,0),
              if (bp2 > 0) fmt(100*(bp3-bp2)/bp2,1) else "NA"), "informational",
      "INFO",
      paste("Expect a DOWNWARD move. The PE expands an unsampled week x day-type stratum at zero",
            "effort, and this change takes that component from 4 of 76 zeroed days to 9 of 76, all",
            "weekend or holiday days. The run log reports the count per component; if the drop is",
            "larger than you are willing to accept, re-run with pe_empty_effort_stratum = 'day_type'."))

  # ---- RUNG 4: interview floor 20 -> 15 -----------------------------------
  m4 <- get1(4, "pooled", "boat_pc_method")
  if (!is.na(m4))
    out[[length(out)+1]] <- V(4, "boat pot closure reaches the sampler",
      as.character(m4), "not 'PE (insufficient data)'",
      if (!grepl("insufficient", m4)) "PASS" else "FAIL",
      paste("At a floor of 20 this component (17 interviews) was never fitted and entered the",
            "port total as an interval-free point. Reaching the sampler is the pass condition;",
            "the gate verdict itself is informational. A 'PE (convergence fail)' here is a GOOD",
            "outcome relative to before -- the failure is now measured rather than assumed."))
  for (nm in list(c("shore_ag_BSS_catch","shore all-gear"), c("shore_pc_BSS_catch","shore pot closure"))) {
    a <- get1(4, "pooled", nm[1]); b <- get1(3, "pooled", nm[1])
    if (!is.na(a) && !is.na(b))
      out[[length(out)+1]] <- V(4, sprintf("%s untouched by the interview floor", nm[2]),
        sprintf("%s vs %s (delta %s)", fmt(a,0), fmt(b,0), fmt(a-b,0)), "delta == 0",
        if (abs(a-b) < 0.5) "PASS" else "FAIL",
        "Both shore components carry hundreds of interviews, far above either floor, so the change must not touch them. Movement means the floor is being applied somewhere it should not be.")
  }

  # ---- RUNG 5: the pooled / gear-resolved cross-check ---------------------
  g5 <- get1(5, "gear_resolved", "port_BSS_catch"); p4 <- get1(4, "pooled", "port_BSS_catch")
  if (!is.na(g5) && !is.na(p4)) {
    pct <- 100 * abs(g5 - p4) / p4
    out[[length(out)+1]] <- V(5, "pooled and gear-resolved port totals reconcile",
      sprintf("%s vs %s (%s%%)", fmt(g5,0), fmt(p4,0), fmt(pct,2)), "< 2%",
      if (pct < 2) "PASS" else if (pct < 5) "REVIEW" else "FAIL",
      paste("Two independent pipelines agreeing on the port total is the strongest internal",
            "cross-check the project has (they agreed to within 1-2% before this patch). A wider",
            "gap after the batch means one track picked up a change the other did not."))
    gb <- get1(5, "gear_resolved", "boat_ag_BSS_catch"); pb <- get1(4, "pooled", "boat_ag_BSS_catch")
    if (!is.na(gb) && !is.na(pb)) {
      pctb <- 100 * abs(gb - pb) / pb
      out[[length(out)+1]] <- V(5, "boat all-gear reconciles across tracks",
        sprintf("%s vs %s (%s%%)", fmt(gb,0), fmt(pb,0), fmt(pctb,2)), "< 3%",
        if (pctb < 3) "PASS" else if (pctb < 6) "REVIEW" else "FAIL",
        "The boat is the component most sensitive to AR resolution and effort-stream handling, so it is where a one-sided change would show first.")
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# DRY RUN: resolve and print the ladder, then stop.
# ---------------------------------------------------------------------------
show_cfg <- function(i) {
  cfg <- resolve_cfg(i)
  keys <- unique(c(names(PRE_PATCH), unlist(lapply(LADDER, function(z) names(z$delta)))))
  for (k in keys) {
    v <- cfg[[k]]
    changed <- !identical(v, modifyList(BASE, PRE_PATCH)[[k]])
    cat(sprintf("      %-30s %-28s %s\n", k,
                paste(as.character(v %||% "NULL"), collapse = ", "),
                if (changed) "  <- changed from the baseline" else ""))
  }
}
if (isTRUE(DRY_RUN)) {
  banner(sprintf("DRY RUN  |  mode = %s  |  rungs = %s  |  NOTHING WILL BE FITTED",
                 MODE, paste(RUNGS, collapse = ", ")))
  for (i in RUNGS) {
    L <- LADDER[[i]]
    rule(); cat(sprintf("  RUNG %d  %-14s  models: %s\n          %s\n",
                        L$id, L$tag, paste(L$models, collapse = " + "), L$headline))
    show_cfg(i)
  }
  rule()
  cat("\n  Reference runs this ladder is judged against:\n")
  for (m in names(REF)) {
    p <- .here("05_output", REF[[m]]$dir)
    cat(sprintf("    %-14s %s  [%s]\n", m, REF[[m]]$dir,
                if (dir.exists(p)) "found" else "MISSING -- rung 1 cannot be judged"))
  }
  cat(sprintf("\n  %d model-runs. Budget roughly 16-32 h wall clock plus a Stan recompile.\n",
              sum(vapply(LADDER[RUNGS], function(z) length(z$models), integer(1)))))

  # ---- SELF-TEST -----------------------------------------------------------
  # Prove the machinery works BEFORE it is trusted with a day of compute. Runs the
  # extractor against the two real reference folders (so a renamed column is caught
  # here rather than after rung 5), then pushes a synthetic clean scenario and a
  # synthetic broken scenario through the criteria to confirm they fire in BOTH
  # directions. A criterion that can only ever say PASS is not a criterion.
  banner("SELF-TEST (extractor + criteria; no fitting)")
  st_ok <- 0; st_bad <- 0
  st <- function(nm, cond, extra = "") {
    if (isTRUE(cond)) { st_ok <<- st_ok + 1; cat(sprintf("  PASS  %s %s\n", nm, extra)) }
    else              { st_bad <<- st_bad + 1; cat(sprintf("  FAIL  %s %s\n", nm, extra)) }
  }
  pd <- .here("05_output", REF$pooled$dir); gd <- .here("05_output", REF$gear_resolved$dir)
  if (dir.exists(pd) && dir.exists(gd)) {
    rp <- extract_run(1, "pooled", "ref", pd, NA_real_, TRUE)
    rg <- extract_run(1, "gear_resolved", "ref", gd, NA_real_, TRUE)
    st("pooled reference parses", isTRUE(rp$port_BSS_catch == REF$pooled$port_catch), rp$port_BSS_catch)
    st("gear reference parses (different header layout)",
       isTRUE(rg$port_BSS_catch == REF$gear_resolved$port_catch), rg$port_BSS_catch)
    st("component catches parse",
       isTRUE(rp$shore_ag_BSS_catch == REF$pooled$shore_ag_catch) &&
       isTRUE(rp$boat_ag_BSS_catch  == REF$pooled$boat_ag_catch))
    st("structural params parse",
       isTRUE(abs(rp$shore_ag_sigma_IE - REF$pooled$shore_ag_sigma_IE) < 1e-3) &&
       isTRUE(abs(rp$shore_ag_B1 - REF$pooled$shore_ag_B1) < 1e-3),
       sprintf("sigma_IE=%s B1=%s", fmt(rp$shore_ag_sigma_IE,3), fmt(rp$shore_ag_B1,3)))
    st("L posterior parses (column renamed 2026-08-25; both names accepted)",
       !is.na(rp$shore_ag_L_med) && !is.na(rp$shore_ag_L_width),
       sprintf("median=%s width=%s", fmt(rp$shore_ag_L_med,2), fmt(rp$shore_ag_L_width,2)))
    st("boat pot closure reads as unfitted in the reference",
       isTRUE(grepl("insufficient", rp$boat_pc_method)), rp$boat_pc_method)

    # 2026-08-27: the synthetic rows now carry shore_ag_method, because the criteria set gained
    # a per-rung "the component still reports a BSS fit" test. `outdir` is left as the REFERENCE
    # folder name so gear_fit_exactness() compares that folder against itself and returns PASS,
    # which is the correct clean-scenario answer and also exercises the comparison itself.
    mk <- function(rn, md, ...) { r <- rp; r$rung <- rn; r$model <- md
      r$shore_ag_method <- "BSS"; r$shore_pc_method <- "BSS"; r$boat_ag_method <- "BSS"
      v <- list(...); for (n in names(v)) r[[n]] <- v[[n]]; r }
    Sc <- rbind(
      mk(1,"gear_resolved", port_BSS_catch = REF$gear_resolved$port_catch,
                            outdir = basename(REF$gear_resolved$dir)),
      mk(1,"pooled"),
      mk(2,"pooled", shore_ag_sigma_IE = 0.32, shore_ag_L_width = rp$shore_ag_L_width * 0.55),
      mk(3,"pooled", shore_ag_sigma_IE = 0.32, shore_ag_B1 = 0.85),
      mk(4,"pooled", shore_ag_sigma_IE = 0.32, shore_ag_B1 = 0.85, boat_pc_method = "BSS"),
      mk(5,"gear_resolved", port_BSS_catch = REF$pooled$port_catch * 0.99))
    Vc <- judge(Sc)
    st("clean scenario produces no FAIL", !any(Vc$verdict == "FAIL"),
       paste(sum(Vc$verdict=="PASS"), "PASS /", sum(Vc$verdict=="REVIEW"), "REVIEW /",
             sum(Vc$verdict=="INFO"), "INFO"))

    Sb <- rbind(
      # A 5% port shift is now needed to trip the (deliberately loose) gear gross-change test;
      # bit-identity is tested separately, and by pointing outdir at a folder that does not
      # exist that criterion returns REVIEW rather than a spurious PASS.
      mk(1,"gear_resolved", port_BSS_catch = REF$gear_resolved$port_catch * 1.05,
                            outdir = "no-such-folder"),
      mk(1,"pooled",        port_BSS_catch = REF$pooled$port_catch + 9000),
      mk(2,"pooled", shore_ag_sigma_IE = 1.04, boat_ag_BSS_catch = rp$boat_ag_BSS_catch + 32,
                     shore_ag_method = "PE (convergence fail)"),
      mk(3,"pooled", shore_ag_B1 = 0.40),
      mk(4,"pooled", boat_pc_method = "PE (insufficient data)",
                     shore_ag_BSS_catch = rp$shore_ag_BSS_catch + 500),
      mk(5,"gear_resolved", port_BSS_catch = REF$pooled$port_catch * 0.6))
    Vb <- judge(Sb)
    want <- c("within Monte Carlo error of the reference","sigma_IE falls","boat is untouched",
              "reaches the sampler","shore all-gear untouched","port totals reconcile",
              "still reports a BSS fit")
    fired <- vapply(want, function(p) {
      v <- Vb$verdict[grepl(p, Vb$criterion, fixed = TRUE)]
      length(v) >= 1 && any(v == "FAIL") }, logical(1))
    st("broken scenario fires every hard criterion", all(fired),
       paste(names(fired)[!fired], collapse = "; "))
  } else {
    st("reference run folders present", FALSE,
       "-- rung 1 cannot be judged; check 05_output/20260804 and /20260805")
  }
  banner(sprintf("SELF-TEST: %d passed, %d failed", st_ok, st_bad))
  if (st_bad > 0)
    cat("  Fix the self-test before starting the batch: a broken extractor or a\n",
        "  criterion that cannot fail would let 20 hours of compute return nothing.\n", sep = "")

  cat("\n  Set DRY_RUN <- FALSE and source again to start.\n\n")
} else {

# ---------------------------------------------------------------------------
# THE BATCH
# ---------------------------------------------------------------------------
summary_path  <- .here("05_output", "patch_validation_2026-08-25_summary.csv")
verdict_path  <- .here("05_output", "patch_validation_2026-08-25_verdicts.csv")
append_row <- function(row) {
  hdr <- !file.exists(summary_path)
  write.table(row, summary_path, sep = ",", row.names = FALSE, col.names = hdr,
              append = !hdr, qmethod = "double")
}

banner(sprintf("PATCH VALIDATION LADDER  |  mode = %s  |  start %s",
               MODE, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("  Rungs: ", paste(RUNGS, collapse = ", "),
    "   Resume: ", RESUME, "\n", sep = "")

collected <- list()
for (i in RUNGS) {
  L   <- LADDER[[i]]
  cfg <- resolve_cfg(i)
  for (m in L$models) {
    run_tag <- L$tag
    folder  <- paste0(prefix[[m]], run_tag)

    done <- find_outdir(m, run_tag)
    if (isTRUE(RESUME) && !is.na(done) &&
        file.exists(file.path(done, "port_total_Dungeness_Kept.csv"))) {
      banner(sprintf("SKIP (already complete): rung %d / %s / %s", L$id, m, folder))
      r <- extract_run(L$id, m, run_tag, done, NA_real_, TRUE)
      collected[[length(collected)+1]] <- r; append_row(r); next
    }

    rc <- cfg
    rc$model <- m; rc$run_tag <- run_tag; rc$run_weather <- FALSE

    banner(sprintf("RUNG %d  %s  |  %s  |  start %s", L$id, folder, L$headline,
                   format(Sys.time(), "%H:%M:%S")))
    cat(sprintf("  ie_shore_obs_unit=%s  days_wkend=%s  bss_min_interviews=%s\n",
                rc$ie_shore_obs_unit, paste(rc$days_wkend, collapse = "/"),
                rc$bss_min_interviews))

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
      message("*** RUNG ", L$id, " / ", m, " FAILED: ", conditionMessage(e))
      # 2026-08-26: the render error is often NOT the cause. A Stan chain failure surfaces
      # downstream (rstan returns an empty stanfit rather than raising), so point the
      # operator at what the sampler itself said before they start guessing.
      od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE),
                     error = function(e2) NA_character_)
      if (!is.na(od) && dir.exists(od)) {
        logs <- list.files(od, pattern = "^stan_console_.*\\.log$", full.names = TRUE)
        if (length(logs)) {
          last <- logs[which.max(file.mtime(logs))]
          message("    Stan console for the last fit attempted (", basename(last), "):")
          message(paste0("      ", utils::tail(readLines(last, warn = FALSE), 8), collapse = "\n"))
        } else {
          message("    No stan_console_*.log in ", od,
                  " -- the run failed before or during the first sampler call.")
        }
      }
      FALSE })
    mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

    r <- extract_run(L$id, m, run_tag, find_outdir(m, run_tag), mins, ok)
    collected[[length(collected)+1]] <- r
    tryCatch(append_row(r), error = function(e)
      message("  (summary row not written: ", conditionMessage(e), ")"))
    banner(sprintf("RUNG %d / %s -> %s (%s min)", L$id, m, if (ok) "OK" else "FAILED", mins))
  }
}

S <- if (length(collected)) do.call(rbind, collected) else NULL

# ---------------------------------------------------------------------------
# THE REPORT
# ---------------------------------------------------------------------------
banner("LADDER SUMMARY")
if (!is.null(S)) {
  print(S[, c("rung","model","ok","minutes","port_BSS_catch","port_PE_catch",
              "shore_ag_BSS_catch","boat_ag_BSS_catch","shore_ag_sigma_IE",
              "shore_ag_B1","boat_pc_method","fits_passed")], row.names = FALSE)
}

banner("PRE-SET CRITERIA")
Vd <- tryCatch(judge(S), error = function(e) { message("  (criteria skipped: ",
                                                        conditionMessage(e), ")"); NULL })
if (!is.null(Vd) && nrow(Vd)) {
  for (k in seq_len(nrow(Vd))) {
    v <- Vd[k, ]
    cat(sprintf("\n  [%s]  rung %s: %s\n", v$verdict, v$rung, v$criterion))
    cat(sprintf("        observed : %s\n        threshold: %s\n", v$observed, v$threshold))
    cat(sprintf("        %s\n", paste(strwrap(v$why, width = 72, prefix = "        "),
                                       collapse = "\n")))
  }
  tryCatch(write.csv(Vd, verdict_path, row.names = FALSE), error = function(e) NULL)
  n_fail <- sum(Vd$verdict == "FAIL"); n_rev <- sum(Vd$verdict == "REVIEW")
  banner(sprintf("VERDICT: %d PASS, %d REVIEW, %d FAIL", sum(Vd$verdict == "PASS"), n_rev, n_fail))
  if (n_fail > 0)
    cat("\n  A FAIL is a finding, not automatically a bug. Read the 'why' text before\n",
        "  changing anything: two of these criteria are testing a HYPOTHESIS (that the\n",
        "  I/E unit mismatch caused GR-9, and that Friday is a weekday), and the\n",
        "  hypothesis being wrong is a legitimate result worth writing down.\n", sep = "")
} else {
  cat("  No criteria could be evaluated (no completed runs found).\n")
}

banner("NEXT")
cat("  summary : ", summary_path,  "\n", sep = "")
cat("  verdicts: ", verdict_path,  "\n", sep = "")
cat("\n  Then update PIPELINE_STATUS.md (Tier 1: the I/E fix validation item) with the\n",
    "  measured numbers, and commit the run folders:\n",
    "    git add 05_output 06_diagnostics/run_patch_validation_2026-08-25.R\n",
    "    git commit -m 'patch validation ladder 2026-08-25' && git push\n", sep = "")
}
