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
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program (see the LICENSE file);
# if not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
###############################################################################
# run_improvement_plan_2026-08-27.R
#
# Executes the runnable stages of development_notes/improvement-plan-2026-08-27.md
# end to end, unattended, cheapest and most decisive first. Start it, walk away, come
# back to 05_output/improvement_plan_2026-08-27_summary.csv and _verdicts.csv.
#
# THE ORDERING IS THE ARGUMENT, so it is worth reading once before starting.
#
#   A  PE-only comparison, minutes, no fitting. Answers plan item 1.4 (the empty-stratum
#      fill) without spending a single MCMC iteration on it, because that toggle changes
#      only the Point Estimator. A full render would have burned ~4 h to learn the same
#      thing.
#   B  Reseed rung 2, one pooled fit. Plan item 1.2. Cheapest decisive run in the plan: it
#      says whether the 2026-08-26 rung-2 failure was one stalled chain (3 of 4 chains
#      behaved) or the I/E unit fix genuinely destabilising the daily-AR shore fit. It runs
#      before stage F because its answer changes how urgent stage F is.
#   C  Boat pot-closure AR alignment. Plan item 3.1: the two tracks fit that component at
#      different resolutions and get 850 against 743.
#   D  shared_tau OFF, the GATE. Plan item 2.1's control. A fixed-seed run that must come
#      back bit-identical to rung 4 on every fitted parameter. If it does not, the model
#      edit is not behaviour-neutral and stage E is meaningless, so E is SKIPPED
#      automatically rather than run and misread.
#   E  shared_tau ON, pooled then gear. The substantive question: does the boat turnover
#      move toward the 2.0-3.0 the OSP/trailer overlap implies, or stay at its 1.2 prior?
#   F  ar_escalate on the shore all-gear component. Plan item 3.2, and deliberately LAST:
#      it is the longest (each failed rung is a whole extra fit) and the least decisive.
#
# WHAT IS NOT HERE, and why:
#   * Plan items 4.1, 4.2 and 4.4 are code and documentation changes, already applied.
#   * Plan item 4.3 (the empty-stratum sign puzzle) is answered by stage A's output rather
#     than by a run, so it is reported there.
#   * Plan item 2.2 (the OSP posterior-predictive check) is a code change that rides along;
#     every stage from B onward writes an `osp` row into ppc_calibration_*.csv, and the
#     summary carries its PIT mean next to the trailer's.
#
# HOW TO RUN
#     Rscript 06_diagnostics/run_improvement_plan_2026-08-27.R      # from the repo root
#     source("06_diagnostics/run_improvement_plan_2026-08-27.R")    # RStudio: Source
#
# START WITH DRY_RUN <- TRUE. It resolves every stage's config, prints it, runs the
# pre-flight data-contract check and the extractor self-test, and fits nothing. It needs no
# Stan toolchain. Then set DRY_RUN <- FALSE.
#
# RESUMABLE. A stage whose output folder already holds port_total_Dungeness_Kept.csv is
# skipped and re-extracted, so an interrupted batch picks up where it stopped. Results are
# appended to the summary CSV as each stage finishes, so a batch killed at hour 20 still
# leaves everything it learned.
#
# RUNTIME. Roughly 25-35 h on 4 cores for all stages, dominated by F. Stage A is minutes.
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                   # TRUE: resolve and print everything, fit nothing. START HERE.
#        ^^^^ reset to TRUE after the 2026-08-29 batch. RESUME skips completed stages, so
#        sourcing this with DRY_RUN <- FALSE again would re-extract and re-judge the existing
#        outputs and APPEND another six rows to the summary CSV rather than re-fitting.
STAGES  <- c("A", "B", "C", "D", "E", "F")
RESUME  <- TRUE                  # skip a stage whose output folder already looks complete
RESEED  <- 20260827L             # stage B's replacement bss_seed (anything != the shipped one)

# STAGES D and E carry a MODEL CHANGE (the shared-turnover parameterization added
# 2026-08-27). It is guarded, off by default, and stage D exists precisely to prove the OFF
# path is bit-identical before stage E is believed. If you would rather eyeball the Stan
# diff before committing compute to it, drop "D" and "E" from STAGES above and run them
# afterwards; nothing else in the batch depends on them.

# =========================================================================== #

# Repo root, resolved from wherever this was launched. Deliberately NOT via here::here, so a
# DRY RUN needs no packages at all.
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
  formatC(x, format = "f", digits = d, big.mark = ",")

if (!isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages({ library(here); library(rmarkdown) })
  load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
  install.lib <- load.lib[!load.lib %in% installed.packages()]
  for (lib in install.lib) install.packages(lib, dependencies = TRUE)
  invisible(sapply(load.lib, require, character.only = TRUE))
  rstan_options(auto_write = TRUE)
  invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), source))
}

source(.here("run_config.R"))          # defines run_config (the SHIPPED config = ladder rung 4)
BASE <- run_config

model_rmd <- list(
  pooled        = .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd"))
prefix <- list(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")
stopifnot(all(file.exists(unlist(model_rmd))))

# ---------------------------------------------------------------------------
# PRE-FLIGHT: the R-to-Stan data contract, the check that would have saved the
# 2026-08-25 ladder's six lost runs. A text parse, no packages, milliseconds. This batch
# adds four Stan data variables (the shared_tau block), so it earns its place again here.
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

# ---------------------------------------------------------------------------
# The 2026-08-26 rung-4 baseline. Hard-coded on purpose: a criterion you can edit after
# seeing the answer is not a criterion. These come from
# 05_output/20260826/pooled-CPUE-PV4-minint and .../gear-type-CPUE-model-PV5-gear-ship.
# ---------------------------------------------------------------------------
REF <- list(
  pooled = list(dir = "20260826/pooled-CPUE-PV4-minint",
                port_catch = 66237, port_lo95 = 50037, port_hi95 = 91210, port_effort = 40187,
                shore_ag = 20898, shore_pc = 6331, boat_ag = 25868, boat_pc = 849,
                # 849 is the EXPECTED-catch median as pe_vs_bss_comparison.csv reports it, which is
                # what extract_run() reads. bss_summary_*.csv shows 850 for C_sum, the PREDICTIVE
                # median; the one-crab gap is the Poisson observation noise, not a discrepancy.
                boat_pc_ar = "monthly", shore_ag_sigma_IE = 0.3004, shore_ag_B1 = 0.9068,
                shore_ag_neff = 6799, shore_ag_div = 333,
                trailer_pit = 0.42, boat_L_med = 1.201),
  gear   = list(dir = "20260826/gear-type-CPUE-model-PV5-gear-ship",
                port_catch = 65756, boat_ag = 25796, boat_pc = 743, boat_pc_ar = "biweekly"))

# ---------------------------------------------------------------------------
# THE STAGES. `delta` is applied on top of the SHIPPED config, never cumulatively: every
# stage is an isolated change against rung 4, so a result is attributable to one thing.
# ---------------------------------------------------------------------------
stage <- function(id, tag, model, delta, headline, plan_item, kind = "run")
  list(id = id, tag = tag, model = model, delta = delta, headline = headline,
       plan_item = plan_item, kind = kind)

STAGE_DEFS <- list(
  A = stage("A", "IP-A-pe-fill", "pooled", list(),
            "empty-stratum fill: 'zero' vs 'day_type', PE only, no fitting",
            "1.4", kind = "pe_only"),
  B = stage("B", "IP-B-reseed", "pooled",
            list(ie_shore_obs_unit = "auto", days_wkend = c("Friday","Saturday","Sunday"),
                 bss_min_interviews = 20, bss_min_interviews_fitted = 0,
                 bss_seed = RESEED),
            "ladder rung 2 on a different seed: stalled chain, or a real instability?",
            "1.2"),
  # 2026-08-27b: this delta was correct in intent and wrong in effect on the first run.
  # ar_force was per-POPULATION only, so the nested list silently forced BOTH boat
  # sub-seasons to biweekly and moved boat all-gear 25,883 -> 28,902 as a side effect. The
  # nested form is now honoured per sub-season (bss_ar_resolution.R), so re-running stage C
  # isolates the pot closure as intended.
  C = stage("C", "IP-C-boatpc-ar", "pooled",
            list(ar_force = list(private_boat = list(pot_closure = "biweekly"))),
            "boat pot closure at the GEAR track's biweekly AR (pooled ran monthly)",
            "3.1"),
  D = stage("D", "IP-D-tau-off", "pooled", list(shared_tau = FALSE),
            "GATE: the shared-turnover edit with the feature OFF must be bit-identical to rung 4",
            "2.1"),
  E = stage("E", "IP-E-tau-on", c("pooled", "gear_resolved"), list(shared_tau = TRUE),
            "shared turnover ON: does tau_boat move toward the 2.0-3.0 the OSP overlap implies?",
            "2.1"),
  F = stage("F", "IP-F-escalate", "pooled",
            list(ar_escalate = TRUE, ar_escalate_respect_cap = FALSE),
            "AR escalation ladder on the near-saturated shore all-gear daily fit",
            "3.2"))

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
sp_get <- function(dir, fitlab, par, col = "median") {
  d <- rd(dir, sprintf("structural_params_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"parameter" %in% names(d)) return(NA_real_)
  num1(d[[col]][d$parameter == par])
}
sp_decoupled <- function(dir, fitlab, par) {
  d <- rd(dir, sprintf("structural_params_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"decoupled" %in% names(d)) return(NA)
  v <- d$decoupled[d$parameter == par]
  if (length(v)) as.logical(v[1]) else NA
}
pit_get <- function(dir, fitlab, stream) {
  d <- rd(dir, sprintf("ppc_calibration_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"data_type" %in% names(d)) return(NA_real_)
  num1(d$pit_mean[d$data_type == stream])
}

extract_run <- function(sid, model, run_tag, outdir, minutes, ok) {
  row <- data.frame(
    stage = sid, plan_item = STAGE_DEFS[[sid]]$plan_item, model = model, run_tag = run_tag,
    ok = ok, minutes = minutes,
    outdir = if (is.na(outdir)) "" else basename(outdir),
    port_BSS_catch = NA_real_, port_BSS_lo95 = NA_real_, port_BSS_hi95 = NA_real_,
    port_BSS_effort = NA_real_, port_PE_catch = NA_real_,
    shore_ag_BSS = NA_real_, shore_pc_BSS = NA_real_,
    boat_ag_BSS = NA_real_, boat_pc_BSS = NA_real_,
    shore_ag_method = NA_character_, boat_ag_method = NA_character_,
    boat_pc_method = NA_character_, boat_pc_ar = NA_character_,
    shore_ag_sigma_IE = NA_real_, shore_ag_B1 = NA_real_, shore_ag_neff = NA_real_,
    shore_ag_div = NA_real_, shore_ag_min_chain_div = NA_real_, shore_ag_max_chain_div = NA_real_,
    boat_tau_bar = NA_real_, boat_L_med = NA_real_, boat_kappa_OSP = NA_real_,
    kappa_OSP_decoupled = NA, boat_r_OSP = NA_real_,
    pit_trailer = NA_real_, pit_osp = NA_real_,
    fits_passed = NA_character_, ar_resolutions = NA_character_,
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
    row$shore_ag_BSS <- g("^shore \\(All gear\\)", "BSS_catch")
    row$shore_pc_BSS <- g("^shore \\(Pot closure\\)", "BSS_catch")
    row$boat_ag_BSS  <- g("^private_boat \\(All gear\\)", "BSS_catch")
    row$boat_pc_BSS  <- g("^private_boat \\(Pot closure\\)", "BSS_catch")
  }
  cr <- rd(outdir, "convergence_report.csv")
  if (!is.null(cr) && "fit" %in% names(cr)) {
    get <- function(pat, col) { i <- grepl(pat, cr$fit); if (any(i)) cr[[col]][i][1] else NA }
    row$shore_ag_method <- as.character(get("^shore_all_gear", "method_selected"))
    row$boat_ag_method  <- as.character(get("^private_boat_all_gear", "method_selected"))
    row$boat_pc_method  <- as.character(get("^private_boat_ring_net_only", "method_selected"))
    row$boat_pc_ar      <- as.character(get("^private_boat_ring_net_only", "ar_resolution"))
    row$shore_ag_neff   <- num1(get("^shore_all_gear", "C_sum_neff"))
    row$shore_ag_div    <- num1(get("^shore_all_gear", "divergences"))
  }
  # Per-CHAIN divergences. The whole point of stage B: one stalled chain and a structural
  # failure are indistinguishable in the aggregate count and completely different findings.
  sdg <- rd(outdir, "sampler_diagnostics_shore_all_gear_Dungeness_Kept.csv")
  if (!is.null(sdg) && "divergent" %in% names(sdg)) {
    row$shore_ag_min_chain_div <- min(suppressWarnings(as.numeric(sdg$divergent)), na.rm = TRUE)
    row$shore_ag_max_chain_div <- max(suppressWarnings(as.numeric(sdg$divergent)), na.rm = TRUE)
  }
  row$shore_ag_sigma_IE <- sp_get(outdir, "shore_all_gear", "sigma_IE")
  row$shore_ag_B1       <- sp_get(outdir, "shore_all_gear", "B1")
  row$boat_tau_bar      <- sp_get(outdir, "private_boat_all_gear", "tau_bar[1]")
  if (is.na(row$boat_tau_bar)) row$boat_tau_bar <- sp_get(outdir, "private_boat_all_gear", "tau_bar")
  row$boat_kappa_OSP    <- sp_get(outdir, "private_boat_all_gear", "kappa_OSP")
  row$kappa_OSP_decoupled <- sp_decoupled(outdir, "private_boat_all_gear", "kappa_OSP")
  row$boat_r_OSP        <- sp_get(outdir, "private_boat_all_gear", "r_OSP")
  row$pit_trailer       <- pit_get(outdir, "private_boat_all_gear", "trailer")
  row$pit_osp           <- pit_get(outdir, "private_boat_all_gear", "osp")

  Ld <- rd(outdir, "bss_L_effective_private_boat_all_gear_Dungeness_Kept.csv")
  if (!is.null(Ld) && "L_posterior_median" %in% names(Ld))
    row$boat_L_med <- num1(stats::median(suppressWarnings(as.numeric(Ld$L_posterior_median)),
                                         na.rm = TRUE))
  ss <- rd(outdir, "season_summary.csv")
  if (!is.null(ss) && "metric" %in% names(ss)) {
    row$fits_passed    <- as.character(ss$value[ss$metric == "BSS fits passed"])[1]
    row$ar_resolutions <- as.character(ss$value[ss$metric == "AR temporal resolution"])[1]
  }
  ee <- rd(outdir, "pe_empty_effort_strata.csv")
  if (!is.null(ee) && "component" %in% names(ee) && nrow(ee) > 0)
    row$empty_effort_note <- paste(sprintf("%s %s/%s", ee$component, ee$n_empty_days,
                                           ee$n_calendar_days), collapse = "; ")
  row
}

# ---------------------------------------------------------------------------
# Bit-identity of a run's FITS against a reference folder (stage D's gate).
# Byte-compare convergence_report.csv, and compare every SHARED parameter row of the per-fit
# summaries at full double precision. Shared rows only, because the shared-turnover edit adds
# `tau_bar` / `tau_bar_out` to the reported set exactly as the 2026-08-25 batch added
# f_lower_out; new REPORTED rows are not a behaviour change, changed values are.
# ---------------------------------------------------------------------------
fit_exactness <- function(new_dirname, ref_dir) {
  if (is.na(new_dirname) || !nzchar(new_dirname))
    return(list(observed = "no output folder found", verdict = "REVIEW"))
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == new_dirname]
  new_dir <- if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NULL
  if (is.null(new_dir) || !dir.exists(new_dir) || !dir.exists(ref_dir))
    return(list(observed = "run or reference folder missing", verdict = "REVIEW"))

  cr_same <- if (file.exists(file.path(ref_dir, "convergence_report.csv")) &&
                 file.exists(file.path(new_dir, "convergence_report.csv"))) {
    a <- file.path(ref_dir, "convergence_report.csv"); b <- file.path(new_dir, "convergence_report.csv")
    identical(readBin(a, "raw", file.size(a)), readBin(b, "raw", file.size(b)))
  } else NA

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
    common <- intersect(rownames(a), rownames(b)); cols <- intersect(names(a), names(b))
    n_rows <- n_rows + length(common)
    if (!isTRUE(all.equal(a[common, cols], b[common, cols], tolerance = 0))) rows_ok <- FALSE
  }
  list(observed = sprintf("convergence_report byte-identical: %s; %d shared parameter rows across %d summaries %s",
                          if (isTRUE(cr_same)) "yes" else "NO", n_rows, length(fs),
                          if (rows_ok) "identical at full precision" else "DIFFER"),
       verdict = if (isTRUE(cr_same) && rows_ok && n_rows > 0) "PASS" else "FAIL")
}

# ---------------------------------------------------------------------------
# STAGE A: the empty-stratum fill, PE only.
#
# pe_empty_effort_stratum changes ONLY the Point Estimator, so a full render would spend
# ~4 h of MCMC to answer a question the PE settles in about a minute. This builds the same
# inputs the driver builds and calls run_pe_pooled() twice per component.
#
# It also answers plan item 4.3. The rung-3 rationale predicted the boat pot-closure PE
# would move DOWN when Friday left the weekend stratum; the count prediction was exactly
# right (4 of 76 zeroed days to 9 of 76) and the PE moved UP 14%. Running both fills side by
# side shows how much of each component's PE is currently being expanded at zero effort,
# which is the quantity that reasoning was really about.
# ---------------------------------------------------------------------------
run_stage_A <- function() {
  banner("STAGE A  |  empty-effort-stratum fill, PE only (plan item 1.4 / 4.3)")
  out_dir <- .here("05_output", format(Sys.Date(), "%Y%m%d"), "improvement-plan-A-pe-fill")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  params <- modifyList(BASE, list(bss_model_file = "crab_bss_pooled.stan"))
  # 2026-08-27b FIX. The driver loads the crabbing-holiday workbook into params at setup
  # (BSS-GH-pooled-CPUE-model.Rmd), and this stage did not, so its day-typing had no holidays
  # and it stratified 289 days into 85 week x day-type cells instead of the driver's 92. The
  # zeroed-day counts it reported were consequently 41/44/8/0 against the driver's 44/47/9/0.
  # The zero-vs-day_type COMPARISON was still fair (both arms used the same stratification),
  # but the absolute counts were not production numbers. Load them the way the driver does.
  params$crabbing_holiday_dates <- read_crabbing_holidays(params)
  dwg     <- fetch_crab_data(params)
  ie_data <- fetch_ie_data(params)
  # Mirrors the driver exactly (BSS-GH-pooled-CPUE-model.Rmd): the day-length model is only
  # built when there is I/E data, and NULL falls back to civil twilight the same way.
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
        n_empty_strata = r$n_empty_effort_strata %||% NA_integer_,
        n_strata_total = r$n_effort_strata_total %||% NA_integer_,
        n_empty_days   = r$n_empty_effort_days   %||% NA_integer_,
        n_calendar_days = r$n_calendar_days      %||% NA_integer_,
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

  rule()
  cat("  component                    zeroed days   PE catch (zero -> day_type)   change\n")
  for (i in seq_len(nrow(w)))
    cat(sprintf("  %-28s %3d/%-3d %5.1f%%   %8s -> %-8s  %+6.1f%%\n",
                w$component[i], w$n_empty_days_zero[i], w$n_calendar_days_zero[i],
                100 * w$empty_day_fraction_zero[i],
                fmt(w$pe_catch_zero[i], 0), fmt(w$pe_catch_daytype[i], 0), w$catch_change_pct[i]))
  rule()
  cat("  Written to", basename(out_dir), "\n")
  cat("  Read this next to pe_vs_bss_comparison.csv from the rung-4 run: where a component\n")
  cat("  reports its BSS, this changes only the cross-check; where it falls back to PE, it\n")
  cat("  changes the reported number.\n")
  invisible(w)
}

# ---------------------------------------------------------------------------
# The run loop.
# ---------------------------------------------------------------------------
sum_path <- .here("05_output", "improvement_plan_2026-08-27_summary.csv")
ver_path <- .here("05_output", "improvement_plan_2026-08-27_verdicts.csv")
append_row <- function(r) utils::write.table(
  r, sum_path, sep = ",", row.names = FALSE, qmethod = "double",
  col.names = !file.exists(sum_path), append = file.exists(sum_path))

# ---- DRY RUN ---------------------------------------------------------------
if (isTRUE(DRY_RUN)) {
  banner(sprintf("DRY RUN  |  stages = %s  |  NOTHING WILL BE FITTED", paste(STAGES, collapse = ", ")))
  for (sid in STAGES) {
    S <- STAGE_DEFS[[sid]]
    rule()
    cat(sprintf("  STAGE %s  (%s)  plan item %s  models: %s\n", S$id,
                if (identical(S$kind, "pe_only")) "PE only, no fitting" else "full run",
                S$plan_item, paste(S$model, collapse = " + ")))
    cat(sprintf("    %s\n", S$headline))
    cfg <- resolve_cfg(sid)
    keys <- names(S$delta)
    if (!length(keys)) cat("    delta: none (shipped config)\n")
    for (k in keys) {
      v <- cfg[[k]]
      cat(sprintf("    delta: %-28s %s\n", k,
                  if (is.null(v)) "NULL" else paste(utils::capture.output(str(v)), collapse = " ")))
    }
  }
  rule()
  banner("SELF-TEST (extractor against the rung-4 baseline; no fitting)")
  st_ok <- 0; st_bad <- 0
  st <- function(nm, cond, extra = "") {
    if (isTRUE(cond)) { st_ok <<- st_ok + 1; cat(sprintf("  PASS  %s %s\n", nm, extra)) }
    else              { st_bad <<- st_bad + 1; cat(sprintf("  FAIL  %s %s\n", nm, extra)) } }
  rdir <- .here("05_output", REF$pooled$dir)
  gdir <- .here("05_output", REF$gear$dir)
  if (dir.exists(rdir)) {
    r <- extract_run("D", "pooled", "ref", rdir, NA_real_, TRUE)
    st("rung-4 port total parses", isTRUE(r$port_BSS_catch == REF$pooled$port_catch), r$port_BSS_catch)
    st("rung-4 components parse",
       isTRUE(r$shore_ag_BSS == REF$pooled$shore_ag) && isTRUE(r$boat_pc_BSS == REF$pooled$boat_pc))
    st("rung-4 boat pot-closure AR parses",
       identical(r$boat_pc_ar, REF$pooled$boat_pc_ar), r$boat_pc_ar)
    st("per-chain divergences parse (stage B depends on these)",
       !is.na(r$shore_ag_min_chain_div) && !is.na(r$shore_ag_max_chain_div),
       sprintf("min=%s max=%s", r$shore_ag_min_chain_div, r$shore_ag_max_chain_div))
    st("boat L median parses", isTRUE(abs(r$boat_L_med - REF$pooled$boat_L_med) < 0.01), r$boat_L_med)
    st("trailer PIT parses",
       !is.na(r$pit_trailer) && abs(r$pit_trailer - REF$pooled$trailer_pit) < 0.05, fmt(r$pit_trailer, 3))
    st("tau_bar absent in the baseline (added by this batch)", is.na(r$boat_tau_bar))
    st("OSP PIT absent in the baseline (added by this batch)", is.na(r$pit_osp))
    ex <- fit_exactness(basename(rdir), rdir)
    st("exactness helper compares a folder to itself", identical(ex$verdict, "PASS"), ex$observed)
  } else st("rung-4 reference folder present", FALSE, rdir)
  if (dir.exists(gdir)) {
    g <- extract_run("E", "gear_resolved", "ref", gdir, NA_real_, TRUE)
    st("gear rung-5 parses", isTRUE(g$port_BSS_catch == REF$gear$port_catch), g$port_BSS_catch)
  } else st("gear rung-5 reference folder present", FALSE, gdir)
  banner(sprintf("SELF-TEST: %d passed, %d failed", st_ok, st_bad))
  cat("\n  Set DRY_RUN <- FALSE and source again to start.\n\n")
} else {

# ---- REAL RUN --------------------------------------------------------------
banner(sprintf("IMPROVEMENT PLAN BATCH  |  stages = %s  |  start %s",
               paste(STAGES, collapse = ", "), format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("  Baseline for every comparison: 2026-08-26 ladder rung 4 (port 66,237).\n")
cat("  Each stage is an ISOLATED change against that baseline, not cumulative.\n")

collected <- list()
skip_E <- FALSE

for (sid in STAGES) {
  S <- STAGE_DEFS[[sid]]

  if (identical(S$kind, "pe_only")) {
    t0 <- Sys.time()
    ok <- tryCatch({ run_stage_A(); TRUE },
                   error = function(e) { message("*** STAGE A FAILED: ", conditionMessage(e)); FALSE })
    banner(sprintf("STAGE A -> %s (%.1f min)", if (ok) "OK" else "FAILED",
                   as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    next
  }

  if (identical(sid, "E") && isTRUE(skip_E)) {
    banner("STAGE E SKIPPED: stage D did not reproduce the baseline, so the shared-turnover edit is NOT behaviour-neutral and an ON run could not be attributed to the feature. Fix D first.")
    next
  }

  for (m in S$model) {
    run_tag <- if (length(S$model) > 1) paste0(S$tag, "-", m) else S$tag
    folder  <- paste0(prefix[[m]], run_tag)

    done <- find_outdir(m, run_tag)
    if (isTRUE(RESUME) && !is.na(done) &&
        file.exists(file.path(done, "port_total_Dungeness_Kept.csv"))) {
      banner(sprintf("SKIP (already complete): stage %s / %s / %s", sid, m, folder))
      r <- extract_run(sid, m, run_tag, done, NA_real_, TRUE)
      collected[[length(collected)+1]] <- r; append_row(r); next
    }

    rc <- resolve_cfg(sid)
    rc$model <- m; rc$run_tag <- run_tag; rc$run_weather <- FALSE

    banner(sprintf("STAGE %s  %s  |  %s  |  start %s", sid, folder, S$headline,
                   format(Sys.time(), "%H:%M:%S")))
    cat(sprintf("  plan item %s | shared_tau=%s ar_escalate=%s bss_seed=%s\n",
                S$plan_item, rc$shared_tau %||% FALSE, rc$ar_escalate %||% FALSE,
                rc$bss_seed %||% NA))

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

    # Stage D is a GATE. If the OFF path did not reproduce the baseline, stage E cannot be
    # attributed to the feature, so do not spend hours producing a number nobody can read.
    if (identical(sid, "D")) {
      ex <- fit_exactness(basename(find_outdir(m, run_tag)),
                          .here("05_output", REF$pooled$dir))
      cat("\n  GATE  shared_tau OFF vs rung 4: ", ex$verdict, "\n        ", ex$observed, "\n", sep = "")
      if (!identical(ex$verdict, "PASS")) {
        skip_E <- TRUE
        message("*** The shared-turnover edit is NOT behaviour-neutral with the feature off. ",
                "Stage E will be skipped.")
      }
    }
  }
}

# ---- VERDICTS --------------------------------------------------------------
S <- if (length(collected)) do.call(rbind, collected) else NULL
V <- list()
V1 <- function(stage, criterion, observed, threshold, verdict, why)
  data.frame(stage = stage, criterion = criterion, observed = observed,
             threshold = threshold, verdict = verdict, why = why, stringsAsFactors = FALSE)
g1 <- function(sid, m, col) {
  if (is.null(S)) return(NA)
  x <- S[[col]][S$stage == sid & S$model == m]; if (length(x)) x[1] else NA
}

# STAGE B
bm <- g1("B", "pooled", "shore_ag_method")
if (!is.na(bm)) {
  mn <- g1("B","pooled","shore_ag_min_chain_div"); mx <- g1("B","pooled","shore_ag_max_chain_div")
  V[[length(V)+1]] <- V1("B", "reseeded rung 2: the shore all-gear fit survives",
    sprintf("%s (per-chain divergences %s to %s)", bm, mn, mx), "BSS",
    if (identical(bm, "BSS")) "PASS" else "FAIL",
    paste("On the 2026-08-26 ladder this fit was rejected at 2,957 divergences, but chains 1, 2",
          "and 4 recorded 163/140/154, in line with rung 1, while chain 3 recorded 2,500 of",
          "2,500 and never moved. PASS means it was a stalled chain and the I/E unit fix is",
          "exonerated; FAIL means the fix genuinely destabilises the daily-AR shore fit and",
          "plan item 3.2 becomes urgent rather than optional. Read the per-chain spread, not",
          "just the verdict: a wide min-to-max is the signature of the stall repeating."))
}

# STAGE C
cc <- g1("C","pooled","boat_pc_BSS"); car <- g1("C","pooled","boat_pc_ar")
if (!is.na(cc))
  V[[length(V)+1]] <- V1("C", "boat pot closure at biweekly AR reconciles the two tracks",
    sprintf("%s at %s (pooled monthly %s, gear biweekly %s)", fmt(cc,0), car,
            fmt(REF$pooled$boat_pc,0), fmt(REF$gear$boat_pc,0)),
    "within 5% of the gear track's 743",
    if (!is.na(cc) && abs(cc - REF$gear$boat_pc) / REF$gear$boat_pc < 0.05) "PASS" else "REVIEW",
    paste("The two tracks fit the same component at different AR resolutions and got 850",
          "against 743, a 13% gap with both passing the gate. If the pooled model at biweekly",
          "lands near 743, the gap was a configuration artefact and the maps should be",
          "aligned. If it does not, the difference is a genuine cross-track disagreement and",
          "belongs in the write-up rather than in the config."))

# STAGE D
dd <- g1("D","pooled","outdir")
if (!is.na(dd) && nzchar(dd)) {
  ex <- fit_exactness(dd, .here("05_output", REF$pooled$dir))
  V[[length(V)+1]] <- V1("D", "shared_tau OFF is bit-identical to rung 4",
    ex$observed, "every shared parameter row identical", ex$verdict,
    paste("tau_bar is declared vector<lower=0>[shared_tau], so it is zero-size when the",
          "feature is off and the unconstrained parameter vector is unchanged -- the same",
          "pattern that made the gear track bit-identical across the 2026-08-25 batch. This",
          "is the control for stage E. A FAIL here means the model edit is not",
          "behaviour-neutral and nothing in stage E can be attributed to the feature."))
}

# STAGE E
et <- g1("E","pooled","boat_tau_bar")
if (!is.na(et)) {
  eb <- g1("E","pooled","boat_ag_BSS"); ep <- g1("E","pooled","port_BSS_catch")
  V[[length(V)+1]] <- V1("E", "shared turnover: does tau_boat move off its prior?",
    sprintf("tau_bar = %s (prior centre 1.20; per-day L median was %s at rung 4)",
            fmt(et,3), fmt(REF$pooled$boat_L_med,3)),
    "informational; the OSP overlap implies 2.0-3.0",
    if (et > 1.5) "MOVED" else "HELD",
    paste("THE question this batch exists to answer. Under the per-day parameterization 148",
          "OSP days left the median at 1.201 against a prior centre of 1.200, while the",
          "OSP/trailer overlap calibration puts the real turnover at 2.0-3.0 and the Phase-1",
          "free kappa_OSP sat at 3.15. MOVED means the OSP data does identify the turnover",
          "once a shared level exists, and the boat number should be re-derived. HELD means",
          "the 1.2 prior survives contact with the data and the conflict lives somewhere",
          "else. Either answer settles it; the previous model could produce neither."))
  V[[length(V)+1]] <- V1("E", "what it does to the boat and the port total",
    sprintf("boat all-gear %s (rung 4 %s); port %s (rung 4 %s)",
            fmt(eb,0), fmt(REF$pooled$boat_ag,0), fmt(ep,0), fmt(REF$pooled$port_catch,0)),
    "informational", "INFO",
    paste("osp_scale_is_tau roughly doubles the private-boat harvest and the boat is about",
          "40% of the port total, so this is the size of the number that was resting on a",
          "prior. Do not publish either figure until the trailer and OSP PIT means below are",
          "read alongside it."))
  pt <- g1("E","pooled","pit_trailer"); po <- g1("E","pooled","pit_osp")
  V[[length(V)+1]] <- V1("E", "trailer PIT moves toward 0.50",
    sprintf("trailer %s (rung 4 %s), OSP %s", fmt(pt,3), fmt(REF$pooled$trailer_pit,2), fmt(po,3)),
    "trailer PIT closer to 0.50 than 0.42",
    if (!is.na(pt) && abs(pt - 0.5) < abs(REF$pooled$trailer_pit - 0.5)) "PASS" else "REVIEW",
    paste("The trailer PIT sat at 0.42 against a nominal 0.50 in all six ladder runs while",
          "every shore and catch stream sat at 0.498-0.516: the model predicted more trailers",
          "than were observed, which is what lambda_E being pulled up by a disagreeing OSP",
          "stream looks like. If a shared turnover lets the two streams reconcile, this is",
          "where it shows. The OSP column is new (plan item 2.2) and has no baseline."))
}

# STAGE F
fr <- g1("F","pooled","ar_resolutions")
if (!is.na(fr))
  V[[length(V)+1]] <- V1("F", "escalation ladder: where does shore all-gear settle?",
    sprintf("resolutions: %s | shore all-gear n_eff %s, divergences %s",
            fr, fmt(g1("F","pooled","shore_ag_neff"),0), fmt(g1("F","pooled","shore_ag_div"),0)),
    "informational", "INFO",
    paste("The pooled shore all-gear daily fit carries p_loo = 109 effective parameters on",
          "311 observations with 41 Pareto-k above 0.7, against the gear track's 11.3 and",
          "zero on the same data, and 143 of its 289 daily AR periods carry no observation.",
          "Read ar_escalation_log.csv for the rung-by-rung verdicts: if it settles coarser",
          "than daily, the coverage rule is choosing a resolution the sampler cannot support."))

VD <- if (length(V)) do.call(rbind, V) else NULL
if (!is.null(VD)) utils::write.csv(VD, ver_path, row.names = FALSE)

banner("BATCH SUMMARY")
if (!is.null(S)) {
  print(S[, c("stage","plan_item","model","ok","minutes","port_BSS_catch",
              "boat_ag_BSS","boat_pc_BSS","shore_ag_method","fits_passed")], row.names = FALSE)
  cat("\n  full summary -> ", sum_path, "\n", sep = "")
}
if (!is.null(VD)) {
  rule()
  for (i in seq_len(nrow(VD)))
    cat(sprintf("[%s] %-7s %s\n        %s\n", VD$stage[i], VD$verdict[i],
                VD$criterion[i], VD$observed[i]))
  cat("\n  verdicts -> ", ver_path, "\n", sep = "")
}
banner(sprintf("DONE  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("  Next: read the verdicts above, then update PIPELINE_STATUS.md and the plan.\n\n")
}
