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
# run_shore_ar_zi_2026-09-03.R
#
# The four next steps from the 2026-09-02 validation review. Start it, walk away, come back
# to
#   05_output/shore_ar_zi_2026-09-03_summary.csv
#   05_output/shore_ar_zi_2026-09-03_verdicts.csv
#   05_output/shore_ar_zi_2026-09-03_ar_ladder.csv    <- the shore AR ladder, the main result
#
# THE ORDERING IS THE ARGUMENT.
#
#   P1  PE incomplete-trip arm alignment, PE ONLY, minutes, no MCMC. A code change that
#       makes the Point Estimator use the same interview frame for the boat gear ratio that
#       the BSS uses for R_G_boat. It moves a number, so it gets measured rather than
#       asserted, and it costs a minute rather than a fit because it touches only the PE.
#
#   G1  Gear track re-run, about 30-40 minutes. The 2026-09-02 driver fix gave the gear
#       boat all-gear fit the per-fit sampler settings it had never had. Until this runs,
#       the cross-check that validates the pooled headline still reads 51,385 instead of
#       about 70,900. Cheapest decisive stage in the batch by a wide margin, so it goes
#       first among the fits.
#
#   Z0  ZINB code present, FEATURE OFF. Two jobs at once, which is why it is not merely a
#       control: (a) it must be bit-identical to the 2026-08-31 production run, proving the
#       Stan edit is behaviour-neutral when off, and (b) it IS the daily rung of the AR
#       ladder, fitted under the same Stan file as the other rungs. Without (b) the ladder
#       would compare three post-edit rungs against a pre-edit production run.
#
#   A1  shore all-gear at EVERY rung of the ladder, in ONE run, via the PRODUCTION toggle
#       (`ar_escalate` scoped to shore/all_gear, `ar_escalate_stop = "all_rungs"`). Each
#       rung's own estimate and 95% interval are logged and tabulated in the HTML report.
#       This replaces three separate forced runs: 7 fits instead of 12, one output folder
#       instead of three, and it exercises the mechanism a future season will actually use
#       rather than an experiment-only override.
#
#   Z1  ZINB ON for the shore. Last because it is the only stage whose value depends on
#       another stage passing (Z0's gate), and because its own negative control is free:
#       zi_catch is set per FIT, so the boat fits in this very run carry the feature off and
#       must come back bit-identical to Z0.
#
# WHAT THE LADDER CAN AND CANNOT DECIDE, stated up front because it bounds the whole stage.
# The production rule agreed with the FW creel team is: start at daily, and if a component
# fails the CONVERGENCE GATE, coarsen until one passes. On the 2026-08-31 production run
# EVERY component passes that gate, shore all-gear at daily included. So the production rule
# applied to this season's data changes nothing, and the ladder's value here is diagnostic:
# it measures how the estimate depends on resolution and whether the gate is the right thing
# to escalate on. The adequacy problems at daily (p_loo 35.2% of n_obs, 41 Pareto k above
# 0.7, coverage_50 0.701 against a nominal 0.500) are invisible to the gate by construction.
# If the ladder shows a coarser rung is materially better on those axes while daily keeps
# passing the gate, that is an argument about what the gate should test, and it needs a
# deliberate decision rather than a config edit.
#
# WHY THE SHORE AR LADDER IS THE MAIN EVENT. On 05_output/20260831/pooled-CPUE-VAL-1-adopted
# the shore all-gear fit carries p_loo = 35.2% of n_obs (109.4 effective parameters on 311
# gear observations), 41 Pareto k above 0.7, and coverage_50 = 0.701 against a nominal 0.500,
# which is +7.1 sampling SDs. The boat all-gear fit on the same run carries 8.0%, one bad k
# and 0.095. That is the signature the daily BOAT AR was rejected for in the Stage 5 batch,
# sitting unexamined on a component worth 29% of the port total. The gear-resolved track fits
# the same component at MONTHLY and is not better but differently wrong (coverage_50 0.035,
# -16.4 SDs), so both ends are bad and nobody has run what is between them. That cross-track
# comparison also differs in three things at once; this ladder differs in ONE.
#
# WHAT IS NOT HERE:
#   * Repository hygiene (untrack .Rproj.user, the 05_output policy) is a git operation and a
#     policy decision, not a run.
#   * The crabbing-fraction pilot, the OSP crab-only column, an external benchmark and the
#     Sunday split are data and field items.
#   * The weather module still forks a stale pre-deployment-scale engine and must be re-based
#     before any LOO from it means anything.
#
# HOW TO RUN
#     Rscript 06_diagnostics/run_shore_ar_zi_2026-09-03.R      # from the repo root
#     source("06_diagnostics/run_shore_ar_zi_2026-09-03.R")    # RStudio: Source
#
# START WITH DRY_RUN <- TRUE. It runs P1 for real (PE only, no Stan toolchain needed),
# resolves and prints every fitted stage, runs the R-to-Stan data-contract pre-flight and the
# extractor self-test, and fits nothing.
#
# RESUMABLE. A stage whose output folder already holds port_total_Dungeness_Kept.csv is
# skipped and re-extracted.
#
# RUNTIME. Roughly 17-20 h on 4 cores: G1 ~0.6 h, Z0 ~4.5 h, A1 ~7-9 h (four shore all-gear
# rungs plus three other components once each; coarser rungs are usually FASTER than daily,
# so this may come in well under), Z1 ~5 h. P1 is minutes. That is 5-7 h less than the
# three-forced-run design it replaces.
#
# ONE STAN RECOMPILE. The pooled model gains three data variables and one guarded parameter,
# so the first pooled stage pays a compile (a few minutes). rstan_options(auto_write = TRUE)
# caches it for the rest.
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- FALSE                  # TRUE: P1 runs, nothing is fitted. START HERE.
STAGES  <- c("P1", "G1", "Z0", "A1", "Z1")
RESUME  <- TRUE                  # skip a fitted stage whose output folder already looks complete

# Z0 is the gate for Z1: if the ZINB edit is not behaviour-neutral with the feature off,
# nothing in Z1 can be attributed to the feature. It does NOT gate the AR ladder, which
# shares the same Stan file across all four of its rungs and so stays internally valid
# either way; only the ladder's comparison against the pre-edit production run would be lost.
GATE_ON_Z0 <- TRUE

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
  .P1_ready <- TRUE
} else {
  # P1 is a Point-Estimator stage, so unlike the pure-CSV desk stages of the 2026-09-01 batch
  # it needs the whole helper library and the tidyverse, though NOT rstan and NOT a Stan
  # toolchain. Attach what is there and degrade gracefully: a dry run on a machine without
  # the data-reading packages should still print every stage's resolved config and run the
  # pre-flights, which is most of what a dry run is for.
  .P1_ready <- FALSE
  if ("P1" %in% STAGES) {
    .P1_ready <- tryCatch({
      suppressPackageStartupMessages({
        library(dplyr); library(tidyr); library(readr); library(lubridate)
        library(readxl); library(here); library(purrr); library(stringr); library(tibble)
      })
      invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), source))
      TRUE
    }, error = function(e) { message("  (P1 needs dplyr/lubridate/readxl and the helper ",
                                     "library: ", conditionMessage(e), ")"); FALSE })
  }
  if (!.P1_ready)
    for (f in c("model_diagnostics.R", "bss_model_adequacy.R", "annotate_decoupled_run.R",
                "bss_sampler_override.R", "bss_ar_resolution.R", "pe_gear_ratio_frame.R"))
      try(sys.source(.here("03_R_functions", f), envir = globalenv()), silent = TRUE)
}

source(.here("run_config.R"))
BASE <- run_config

model_rmd <- list(
  pooled        = .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd"))
prefix <- list(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")
stopifnot(all(file.exists(unlist(model_rmd))))
desk_dir <- .here("05_output", format(Sys.Date(), "%Y%m%d"), "shore-ar-zi-desk")

rd <- function(dir, f) {
  p <- file.path(dir, f); if (!file.exists(p)) return(NULL)
  tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}
num1 <- function(x) { x <- suppressWarnings(as.numeric(x)); if (length(x)) x[1] else NA_real_ }
chr1 <- function(x) if (length(x)) as.character(x)[1] else NA_character_

V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = observed,
                                    threshold = threshold, verdict = verdict, why = why,
                                    stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# BASELINES. Hard-coded from committed folders, named beside each value. Port totals are
# never used in an exactness gate: rstan::extract(permuted = TRUE) permutes draws, so a port
# total assembled by resampling component draws is RNG-sensitive even when the fits are
# bit-identical. Exactness is always tested on per-fit summaries.
# ---------------------------------------------------------------------------
REF <- list(
  # 2026-08-31 production, the adopted configuration. The daily-AR rung comes from here.
  PROD = list(dir = "20260831/pooled-CPUE-VAL-1-adopted",
              port = 71513, shore_ag = 20898, shore_pc = 6331, boat_ag = 31008, boat_pc = 1018,
              shore_ag_ar = "daily",
              # the adequacy that makes the shore ladder worth running
              shore_ag_p_loo_frac = 0.352, shore_ag_bad_k = 41, shore_ag_cov50 = 0.701,
              boat_ag_p_loo_frac = 0.080, boat_ag_bad_k = 1,
              # the zero bin the ZINB prototype targets
              shore_ag_obs_zeros = 676, shore_ag_exp_zeros = 605.4, shore_ag_zero_z = 3.8),
  # 2026-09-01 gear run at the OLD (broken) sampler settings, and the 2026-08-30 stage that
  # showed what the fix is worth. G1 must land on the S3 side of this.
  GEARBAD = list(dir = "20260901/gear-type-CPUE-model-VAL-2-gearG-phase1",
                 port = 51385, boat_ag_reported = 11176, divergences = 554,
                 tau_neff = 49, tau_Rhat = 1.0813),
  GEARFIX = list(dir = "20260829/gear-type-CPUE-model-S5-3-tau-gear",
                 port = 70886, boat_ag = 30760, boat_pc = 956, tau_bar = 2.5962,
                 divergences = 2, tau_neff = 19292),
  # The gear track fits shore all-gear at MONTHLY; the other end of the bracket.
  GEARMON = list(shore_ag_p_loo_frac = 0.036, shore_ag_bad_k = 0, shore_ag_cov50 = 0.035))

stage <- function(id, tag, model, delta, headline, item, kind = "run")
  list(id = id, tag = tag, model = model, delta = delta, headline = headline, item = item, kind = kind)

# The nested ar_force form, honoured per sub-season since 2026-08-29 and errors on an unnamed
# list. It bypasses both the data-driven selector and ar_max_resolution, which is what makes
# a forced ladder possible at all: shore all-gear is CAPPED at daily, so the cap alone cannot
# express weekly or coarser.
AR_SHORE_AG <- function(res) list(shore = list(all_gear = res))

STAGE_DEFS <- list(
  P1 = stage("P1", NA, "pooled", NULL,
             "PE / BSS incomplete-trip arm alignment, PE only, no fitting", "arm", "pe_only"),
  G1 = stage("G1", "SZ-G1-gearfix", "gear_resolved", list(),
             "gear track with the 2026-09-02 boat all-gear sampler fix: restore the cross-check",
             "gear fix"),
  Z0 = stage("Z0", "SZ-Z0-zi-off", "pooled", list(estimate_catch_zi = FALSE),
             "GATE and the ladder's DAILY rung: the ZINB edit must be bit-identical when off",
             "zinb"),
  # 2026-09-04: ONE stage replaces the three forced rungs A1/A2/A3.
  #
  # `ar_escalate_stop = "all_rungs"` scoped to shore all-gear fits that component at every
  # rung of the ladder INSIDE A SINGLE RENDER and logs each rung's own estimate and interval
  # to ar_escalation_log.csv, which the HTML report now tabulates. Three separate forced
  # runs fit 12 components to answer the same question; this fits 7 (4 shore all-gear rungs
  # plus the three untouched components once each), keeps every rung in one output folder,
  # and removes the need for a leak control, because the other three fits happen once and
  # cannot be reached by a per-sub-season escalation scope.
  #
  # It also EXERCISES THE PRODUCTION MECHANISM rather than a batch-only override. `ar_force`
  # is an experiment lever; `ar_escalate` is the toggle that will actually be used in future
  # seasons, so testing it here is worth more than testing three forced runs that no
  # production run will ever reproduce.
  A1 = stage("A1", "SZ-A1-shore-ladder", "pooled",
             list(ar_escalate = list(shore = "all_gear"),
                  ar_escalate_stop = "all_rungs",
                  ar_escalate_select = "first_pass",
                  ar_escalate_respect_cap = FALSE),
             "shore all-gear at EVERY rung in one run: the ladder, with each rung's estimate kept",
             "shore AR"),
  Z1 = stage("Z1", "SZ-Z1-zi-shore", "pooled",
             list(estimate_catch_zi = TRUE, catch_zi_populations = c("shore")),
             "ZINB on the SHORE catch likelihood; the boat fits are a free negative control",
             "zinb"))

resolve_cfg <- function(sid) modifyList(BASE, STAGE_DEFS[[sid]]$delta %||% list(), keep.null = TRUE)

find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}

# ---------------------------------------------------------------------------
# CONFIG DELTA, carried over from the 2026-09-01 batch. An exactness comparison is only
# meaningful when the two runs differ in ONE thing, and that was got wrong twice in three
# weeks by hand, both times producing a FAIL that looked like a code defect. `expect_delta`
# turns it into a check: any config key that differs but was not expected is reported inside
# the verdict.
# ---------------------------------------------------------------------------
.cfg_keys <- function(dir) {
  p <- file.path(dir, "run_parameters.txt")
  if (!file.exists(p)) return(NULL)
  kv <- list()
  for (l in readLines(p, warn = FALSE)) {
    m <- regmatches(l, regexec("^\\s*\\$ ([A-Za-z0-9_.]+)\\s*:(.*)$", l))[[1]]
    if (length(m) == 3) kv[[m[2]]] <- trimws(m[3])
  }
  kv
}
config_delta <- function(dir_a, dir_b) {
  a <- .cfg_keys(dir_a); b <- .cfg_keys(dir_b)
  if (is.null(a) || is.null(b)) return(NA_character_)
  ks <- setdiff(union(names(a), names(b)), c("run_tag", "model"))
  sort(ks[vapply(ks, function(k) !identical(a[[k]] %||% "<absent>", b[[k]] %||% "<absent>"), logical(1))])
}

fit_exactness <- function(new_dir, ref_dir, pat = NULL, what = "fits", expect_delta = NULL) {
  if (!dir.exists(new_dir %||% "") || !dir.exists(ref_dir %||% ""))
    return(list(observed = "run or reference folder missing", verdict = "REVIEW"))
  fs <- c(list.files(new_dir, pattern = "^bss_summary_.*\\.csv$"),
          list.files(new_dir, pattern = "^bss_full_summary_.*\\.csv$"))
  if (!is.null(pat)) fs <- fs[grepl(pat, fs)]
  fs <- fs[file.exists(file.path(ref_dir, fs))]
  n <- 0L; bad <- character(0)
  for (f in fs) {
    a <- tryCatch(utils::read.csv(file.path(ref_dir, f), row.names = 1, check.names = FALSE), error = function(e) NULL)
    b <- tryCatch(utils::read.csv(file.path(new_dir, f), row.names = 1, check.names = FALSE), error = function(e) NULL)
    if (is.null(a) || is.null(b)) { bad <- c(bad, f); next }
    cr <- intersect(rownames(a), rownames(b)); cc <- intersect(names(a), names(b))
    n <- n + length(cr)
    if (!isTRUE(all.equal(a[cr, cc], b[cr, cc], tolerance = 0))) bad <- c(bad, f)
  }
  cd <- config_delta(ref_dir, new_dir)
  extra <- if (length(cd) && !all(is.na(cd))) setdiff(cd, expect_delta %||% character(0)) else character(0)
  note <- if (!length(cd) || all(is.na(cd))) "" else
    sprintf("; configs differ in %d key(s): %s%s", length(cd), paste(cd, collapse = ", "),
            if (length(extra)) sprintf("  <-- %d UNEXPECTED: %s", length(extra), paste(extra, collapse = ", ")) else "")
  list(observed = sprintf("%s: %d shared parameter rows across %d summaries %s%s", what, n, length(fs),
                          if (!length(bad)) "identical at full precision"
                          else sprintf("DIFFER in %s", paste(basename(bad), collapse = ", ")), note),
       verdict = if (!length(bad) && n > 0) "PASS" else "FAIL", unexpected_delta = extra)
}

# ---------------------------------------------------------------------------
# Extraction.
# ---------------------------------------------------------------------------
sp_get <- function(dir, fitlab, par, col = "median") {
  d <- rd(dir, sprintf("structural_params_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"parameter" %in% names(d) || !col %in% names(d)) return(NA_real_)
  num1(d[[col]][d$parameter == par])
}
sp_idx <- function(dir, fitlab, par, col) {
  v <- sp_get(dir, fitlab, paste0(par, "[1]"), col); if (is.na(v)) v <- sp_get(dir, fitlab, par, col); v
}
adq_get <- function(dir, fitlab, col) {
  d <- rd(dir, "model_adequacy.csv") %||% rd(dir, "model_adequacy_reconstructed.csv")
  if (is.null(d) || !"fit" %in% names(d) || !col %in% names(d)) return(NA)
  v <- d[[col]][d$fit == paste0(fitlab, "_Dungeness_Kept")]; if (!length(v)) NA else v[1]
}
loo_get <- function(dir, fitlab, stream, col) {
  d <- rd(dir, sprintf("loo_summary_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"stream" %in% names(d) || !col %in% names(d)) return(NA_real_)
  num1(d[[col]][d$stream == stream])
}
# Coverage from the RANDOMIZED PIT. ppc_byobs has always used it; ppc_calibration only since
# 2026-09-01, and before that it over-covered small counts by construction.
cov50_byobs <- function(dir, fitlab, stream) {
  d <- rd(dir, sprintf("ppc_byobs_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !all(c("data_type", "in_50") %in% names(d))) return(NA_real_)
  y <- d$in_50[d$data_type == stream]; if (!length(y)) return(NA_real_)
  mean(as.logical(y), na.rm = TRUE)
}
zero_bin <- function(dir, fitlab, stream) {
  d <- rd(dir, sprintf("ppc_byobs_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !all(c("data_type", "p_zero", "observed") %in% names(d))) return(c(NA, NA, NA))
  y <- d[d$data_type == stream, , drop = FALSE]; if (!nrow(y)) return(c(NA, NA, NA))
  p <- suppressWarnings(as.numeric(y$p_zero)); p <- p[is.finite(p)]
  c(sum(y$observed == 0, na.rm = TRUE), sum(p), sqrt(sum(p * (1 - p))))
}

extract_run <- function(sid, model, run_tag, outdir, minutes, ok) {
  row <- data.frame(stage = sid, item = STAGE_DEFS[[sid]]$item %||% NA_character_,
    model = model, run_tag = run_tag, ok = ok, minutes = minutes,
    outdir = if (is.na(outdir)) "" else basename(outdir),
    port_BSS_catch = NA_real_, shore_ag_BSS = NA_real_, shore_pc_BSS = NA_real_,
    boat_ag_BSS = NA_real_, boat_pc_BSS = NA_real_,
    boat_ag_reported = NA, boat_ag_method = NA_character_,
    shore_ag_ar = NA_character_, boat_ag_ar = NA_character_,
    shore_ag_div = NA_real_,
    # the shore-ladder axes
    shore_ag_p_loo_frac = NA_real_, shore_ag_bad_k = NA_real_,
    shore_ag_cov50_gear = NA_real_, shore_ag_cov50_catch = NA_real_,
    shore_ag_reported = NA,
    shore_ag_elpd_gear = NA_real_, shore_ag_elpd_catch = NA_real_,
    shore_ag_se_elpd_gear = NA_real_, shore_ag_se_elpd_catch = NA_real_,
    shore_ag_ploo_gear = NA_real_, shore_ag_ploo_catch = NA_real_,
    shore_ag_phi_E = NA_real_, shore_ag_sigma_eps_E = NA_real_, shore_ag_disp_scale = NA_real_,
    # the ZINB axes
    shore_ag_theta_C = NA_real_, shore_ag_theta_lo95 = NA_real_, shore_ag_theta_hi95 = NA_real_,
    shore_ag_zero_obs = NA_real_, shore_ag_zero_exp = NA_real_, shore_ag_zero_z = NA_real_,
    shore_pc_theta_C = NA_real_,
    boat_tau_bar = NA_real_, boat_ag_div = NA_real_, boat_tau_neff = NA_real_,
    fits_passed = NA_character_, stringsAsFactors = FALSE)
  if (is.na(outdir) || !dir.exists(outdir)) return(row)

  pt <- rd(outdir, "port_total_Dungeness_Kept.csv")
  if (!is.null(pt) && "Estimate" %in% names(pt)) {
    i <- pt$Estimate %in% c("Expected_Catch", "Catch"); row$port_BSS_catch <- num1(pt$BSS_median[i]) }
  pv <- rd(outdir, "pe_vs_bss_comparison.csv")
  if (!is.null(pv) && "component" %in% names(pv)) {
    g <- function(pat, col) if (col %in% names(pv)) { v <- pv[[col]][grepl(pat, pv$component)]; if (length(v)) v[1] else NA } else NA
    row$shore_ag_BSS <- num1(g("^shore \\(All gear\\)", "BSS_catch"))
    row$shore_ag_reported <- as.logical(g("^shore \\(All gear\\)", "bss_reported"))
    row$shore_pc_BSS <- num1(g("^shore \\(Pot closure\\)", "BSS_catch"))
    row$boat_ag_BSS  <- num1(g("^private_boat \\(All gear\\)", "BSS_catch"))
    row$boat_pc_BSS  <- num1(g("^private_boat \\(Pot closure\\)", "BSS_catch"))
    row$boat_ag_reported <- as.logical(g("^private_boat \\(All gear\\)", "bss_reported"))
  }
  cr <- rd(outdir, "convergence_report.csv")
  if (!is.null(cr) && "fit" %in% names(cr)) {
    g <- function(pat, col) { if (!col %in% names(cr)) return(NA); i <- grepl(pat, cr$fit); if (any(i)) cr[[col]][i][1] else NA }
    row$shore_ag_ar     <- chr1(g("^shore_all_gear", "ar_resolution"))
    row$boat_ag_ar      <- chr1(g("^private_boat_all_gear", "ar_resolution"))
    row$shore_ag_div    <- num1(g("^shore_all_gear", "divergences"))
    row$boat_ag_div     <- num1(g("^private_boat_all_gear", "divergences"))
    row$boat_ag_method  <- chr1(g("^private_boat_all_gear", "method_selected"))
  }
  fds <- rd(outdir, "fit_data_summary.csv")
  if (!is.null(fds) && "fit" %in% names(fds) && "ar_resolution" %in% names(fds)) {
    if (is.na(row$shore_ag_ar)) row$shore_ag_ar <- chr1(fds$ar_resolution[fds$fit == "shore_all_gear_Dungeness_Kept"])
    if (is.na(row$boat_ag_ar))  row$boat_ag_ar  <- chr1(fds$ar_resolution[fds$fit == "private_boat_all_gear_Dungeness_Kept"])
  }

  sag <- "shore_all_gear"; bag <- "private_boat_all_gear"
  row$shore_ag_p_loo_frac  <- num1(adq_get(outdir, sag, "p_loo_frac"))
  row$shore_ag_bad_k       <- num1(adq_get(outdir, sag, "n_pareto_bad"))
  row$shore_ag_disp_scale  <- num1(adq_get(outdir, sag, "disp_scale_min"))
  row$shore_ag_cov50_gear  <- cov50_byobs(outdir, sag, "gear")
  row$shore_ag_cov50_catch <- cov50_byobs(outdir, sag, "catch")
  row$shore_ag_elpd_gear   <- loo_get(outdir, sag, "gear", "elpd_loo")
  row$shore_ag_elpd_catch  <- loo_get(outdir, sag, "catch", "elpd_loo")
  # loo_summary_*.csv carries se_elpd_loo. A bare "more than N nats" threshold is a rule of
  # thumb; an elpd difference is only evidence when it is large relative to its own SE.
  row$shore_ag_se_elpd_gear  <- loo_get(outdir, sag, "gear", "se_elpd_loo")
  row$shore_ag_se_elpd_catch <- loo_get(outdir, sag, "catch", "se_elpd_loo")
  row$shore_ag_ploo_gear   <- loo_get(outdir, sag, "gear", "p_loo")
  row$shore_ag_ploo_catch  <- loo_get(outdir, sag, "catch", "p_loo")
  row$shore_ag_phi_E       <- sp_get(outdir, sag, "phi_E")
  row$shore_ag_sigma_eps_E <- sp_get(outdir, sag, "sigma_eps_E")
  row$shore_ag_theta_C     <- sp_idx(outdir, sag, "theta_C", "median")
  row$shore_ag_theta_lo95  <- sp_idx(outdir, sag, "theta_C", "lo95")
  row$shore_ag_theta_hi95  <- sp_idx(outdir, sag, "theta_C", "hi95")
  row$shore_pc_theta_C     <- sp_idx(outdir, "shore_ring_net_only", "theta_C", "median")
  z <- zero_bin(outdir, sag, "catch")
  row$shore_ag_zero_obs <- z[1]; row$shore_ag_zero_exp <- z[2]
  row$shore_ag_zero_z   <- if (is.finite(z[3]) && z[3] > 0) (z[1] - z[2]) / z[3] else NA_real_
  row$boat_tau_bar   <- sp_idx(outdir, bag, "tau_bar", "median")
  row$boat_tau_neff  <- sp_idx(outdir, bag, "tau_bar", "n_eff")
  ss <- rd(outdir, "season_summary.csv")
  if (!is.null(ss) && "metric" %in% names(ss)) row$fits_passed <- chr1(ss$value[ss$metric == "BSS fits passed"])
  row
}

# ===========================================================================
# P1. The PE / BSS incomplete-trip arm alignment, PE ONLY.
#
# The BSS learns R_G_boat from a frame that HAS been incomplete-trip filtered
# (prep_bss_crab_pooled.R: intA descends from int_d, after filter_incomplete_trips). Both
# Point Estimators built the same ratio from the UNFILTERED set. sensitivity_incomplete_trips
# .csv has reported that as production_arm_bss = "exclude" against production_arm_pe =
# "gear_only" since 2026-08-25 without anyone acting on it.
#
# The fix is a shared helper (pe_gear_ratio_frame.R) rather than the same edit in two files,
# because the gear-resolved PE's comment ALREADY claimed it matched the BSS while it did not,
# which is what two copies of one rule produces.
#
# This stage runs both arms of the Point Estimator and reports the difference. It is a
# PE-only change, so it costs a minute rather than a fit, and it moves a reported number
# wherever a component falls back to PE, which makes measuring it non-optional.
# ===========================================================================
run_P1 <- function() {
  banner("STAGE P1  |  PE / BSS incomplete-trip arm alignment, PE only")
  dir.create(desk_dir, recursive = TRUE, showWarnings = FALSE)
  params <- modifyList(BASE, list(bss_model_file = "crab_bss_pooled.stan"))
  params$crabbing_holiday_dates <- read_crabbing_holidays(params)
  dwg     <- fetch_crab_data(params)
  ie_data <- fetch_ie_data(params)
  L_eff   <- if (!is.null(ie_data) && nrow(ie_data) > 0)
               estimate_L_effective(ie_data, params$pot_open_date, params) else NULL
  subs    <- build_subseasons(params)

  rows <- list()
  for (arm in c("gear_only", "match_bss")) {
    p <- modifyList(params, list(pe_gear_ratio_arm = arm))
    for (pop in c("shore", "private_boat")) for (ss in subs) {
      lab <- paste0(pop, "_", ss$name)
      days_ss <- prep_days_crab(ss$start, ss$end, p, L_eff_model = L_eff)
      summ_ss <- prep_population_summary(dwg, pop, ss$start, ss$end, p)
      r <- run_pe_pooled(summ_ss, days_ss, p, lab)
      rows[[length(rows) + 1]] <- data.frame(
        arm = arm, component = lab,
        pe_effort = r$effort_total %||% NA_real_,
        pe_catch  = r[["Dungeness_Kept"]] %||% NA_real_, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  w <- merge(df[df$arm == "gear_only", ], df[df$arm == "match_bss", ],
             by = "component", suffixes = c("_gear_only", "_match_bss"))
  w$catch_change_pct  <- round(100 * (w$pe_catch_match_bss  - w$pe_catch_gear_only)  / pmax(w$pe_catch_gear_only, 1), 3)
  w$effort_change_pct <- round(100 * (w$pe_effort_match_bss - w$pe_effort_gear_only) / pmax(w$pe_effort_gear_only, 1), 3)
  utils::write.csv(w, file.path(desk_dir, "P1_pe_arm_alignment.csv"), row.names = FALSE)

  rule()
  cat(sprintf("  %-28s %14s %14s %10s %10s\n", "component", "PE gear_only", "PE match_bss", "d catch", "d effort"))
  for (i in seq_len(nrow(w)))
    cat(sprintf("  %-28s %14s %14s %+9.3f%% %+9.3f%%\n", w$component[i],
                fmt(w$pe_catch_gear_only[i], 0), fmt(w$pe_catch_match_bss[i], 0),
                w$catch_change_pct[i], w$effort_change_pct[i]))
  rule()

  i_worst <- which.max(abs(w$catch_change_pct))
  worst <- w$catch_change_pct[i_worst]      # signed: the direction is part of the answer
  # The shore is expected to be untouched: the alignment changes only the BOAT gear-per-group
  # ratio, and the shore PE does not use one.
  sh <- w[grepl("^shore", w$component), , drop = FALSE]
  V1row("P1", "aligning the PE onto the BSS arm moves only the boat, and moves it little",
    sprintf("largest change %+.3f%% (%s); shore components %s",
            worst, w$component[i_worst],
            if (all(abs(sh$catch_change_pct) < 1e-6)) "unchanged, as they must be"
            else sprintf("MOVED by up to %+.3f%%", max(abs(sh$catch_change_pct)))),
    "boat only, and under 1%",
    if (all(abs(sh$catch_change_pct) < 1e-6) && abs(worst) < 1) "PASS" else "REVIEW",
    paste("Two arms of one fused estimator disagreeing about which interviews count is a",
          "documented inconsistency a reviewer asks about first, and the desk read on",
          "2026-09-01 sized it at 3.552 against 3.550 gear per boat group. The shore is the",
          "control here: the alignment touches the boat gear-per-group ratio only, and the",
          "shore PE has no such ratio, so any shore movement means the shared helper is",
          "reaching further than intended. A change under 1% on the boat is the expected",
          "result and is why this is a CONSISTENCY fix rather than an accuracy one; the",
          "reason to make it is that there is no estimate to defend, only an inconsistency",
          "to remove."))
  invisible(w)
}

# ---------------------------------------------------------------------------
# The run loop.
# ---------------------------------------------------------------------------
sum_path <- .here("05_output", "shore_ar_zi_2026-09-03_summary.csv")
ver_path <- .here("05_output", "shore_ar_zi_2026-09-03_verdicts.csv")
lad_path <- .here("05_output", "shore_ar_zi_2026-09-03_ar_ladder.csv")

# ---------------------------------------------------------------------------
# RESULT-FILE PERSISTENCE. Read this before changing it; two ways of getting it wrong have
# already destroyed real results in this project.
#
# 1. APPEND-ON-RESUME DUPLICATED ROWS. Every batch runner before this one appended to the
#    summary whenever the file existed, and the RESUME skip path appends too. RESUME is the
#    documented recovery for an interrupted multi-hour batch, so a resume wrote a SECOND row
#    for every completed stage. 05_output/osp_validation_summary.csv carries 28 rows for a
#    14-job matrix and patch_validation_2026-08-25_summary.csv carries 12 for 6 jobs, both
#    from exactly this.
#
# 2. UNGATED OVERWRITE LOST THE FITTED STAGES' VERDICTS. The verdicts file was rewritten
#    unconditionally from whatever the CURRENT invocation produced. A batch with desk stages
#    that run during a DRY RUN therefore truncates the file to desk rows the moment anyone
#    re-sources it, and the 2026-09-01 batch's own closing message told the reader to do
#    precisely that. 05_output/validation_2026-09-01_verdicts.csv now holds 11 desk rows and
#    none of the 7 fitted-stage criteria, while the summary still lists all five fits: the
#    runs happened, only their scored criteria were erased.
#
# The fix for both is the same: MERGE BY KEY, never append and never blind-overwrite. A
# re-run of one stage replaces that stage's rows and leaves every other row alone.
# ---------------------------------------------------------------------------
merge_csv_by <- function(new, path, key) {
  if (is.null(new) || !nrow(new)) return(invisible(NULL))
  old <- if (file.exists(path))
    tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e) NULL) else NULL
  if (!is.null(old) && nrow(old) && all(key %in% names(old)) && all(key %in% names(new))) {
    k_old <- do.call(paste, c(old[key], sep = "\r"))
    k_new <- do.call(paste, c(new[key], sep = "\r"))
    old <- old[!(k_old %in% k_new), , drop = FALSE]
    cols <- intersect(names(old), names(new))
    # A schema change between runs is possible (a new column added mid-project). Keep the
    # NEW schema and carry the old rows across on the shared columns rather than dropping
    # them, so an old row survives with blanks instead of vanishing.
    if (nrow(old)) {
      pad <- old[, cols, drop = FALSE]
      for (cn in setdiff(names(new), cols)) pad[[cn]] <- NA
      new <- rbind(pad[, names(new), drop = FALSE], new)
    }
  }
  utils::write.csv(new, path, row.names = FALSE)
  invisible(new)
}
append_row <- function(r) merge_csv_by(r, sum_path, "stage")

# PRE-FLIGHT 1: the R-to-Stan data contract. This batch adds THREE Stan data variables to the
# pooled model, so it earns its place: an unforwarded field makes rstan return an EMPTY
# stanfit rather than raise, and the 2026-08-25 ladder lost six runs to exactly that.
local({
  src_fn <- .here("03_R_functions", "bss_stan_fit.R")
  if (!file.exists(src_fn)) { warning("bss_stan_fit.R not found; pre-flight skipped."); return(invisible(NULL)) }
  e <- new.env(); sys.source(src_fn, envir = e)
  pairs <- list(pooled = c("crab_bss_pooled.stan", "prep_bss_crab_pooled.R"),
                gear_resolved = c("crab_bss_gear_resolved.stan", "prep_bss_crab_gear.R"))
  fail <- character(0)
  for (m in names(pairs)) {
    need <- e$bss_stan_data_names(.here("02_stan_models", pairs[[m]][1]))
    src  <- paste(readLines(.here("03_R_functions", pairs[[m]][2]), warn = FALSE), collapse = "\n")
    miss <- need[!vapply(need, function(v) grepl(paste0("(^|[^A-Za-z0-9_.])", v, "([^A-Za-z0-9_]|$)"), src, perl = TRUE), logical(1))]
    cat(sprintf("  pre-flight  %-14s %2d Stan data variables, %d not built in %s\n",
                m, length(need), length(miss), basename(pairs[[m]][2])))
    if (length(miss)) fail <- c(fail, sprintf("%s: %s", m, paste(miss, collapse = ", ")))
  }
  if (length(fail)) stop("PRE-FLIGHT FAILED, the R-to-Stan data contract is broken.\n  ",
                         paste(fail, collapse = "\n  "), call. = FALSE)
})

# PRE-FLIGHT 2: the two code changes this batch tests are actually present, and production
# still ships the ZINB off. A stage that silently tests the old code is worse than no stage.
local({
  ok_helper <- file.exists(.here("03_R_functions", "pe_gear_ratio_frame.R"))
  gear_rmd  <- paste(readLines(model_rmd$gear_resolved, warn = FALSE), collapse = "\n")
  ok_gear   <- grepl("bss_iter_boat_allgear", gear_rmd, fixed = TRUE)
  cat(sprintf("  pre-flight  PE arm helper present: %s | gear boat sampler fix present: %s | ships ZINB off: %s\n",
              ok_helper, ok_gear, isFALSE(BASE$estimate_catch_zi)))
  if (!ok_helper) stop("03_R_functions/pe_gear_ratio_frame.R is missing; P1 has nothing to test.", call. = FALSE)
  if (!ok_gear)   stop("The gear driver has no bss_iter_boat_allgear; G1 would reproduce the ",
                       "broken fit and look like a null result.", call. = FALSE)
  if (!isFALSE(BASE$estimate_catch_zi))
    stop("run_config.R ships estimate_catch_zi = TRUE, so Z0 would not be an OFF control.", call. = FALSE)
})

collected <- list(); gate_ok <- TRUE

for (sid in STAGES) {
  S <- STAGE_DEFS[[sid]]

  if (identical(S$kind, "pe_only")) {
    if (isTRUE(DRY_RUN) && !isTRUE(.P1_ready)) {
      banner("STAGE P1 DEFERRED: the packages it needs are not available in this session. It will run first when DRY_RUN <- FALSE, still without any fitting.")
      next
    }
    t0 <- Sys.time()
    tryCatch(run_P1(), error = function(e) message("*** STAGE P1 FAILED: ", conditionMessage(e)))
    cat(sprintf("\n  [P1 done in %.1f s]\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    next
  }
  if (isTRUE(DRY_RUN)) next

  if (!gate_ok && isTRUE(GATE_ON_Z0) && identical(sid, "Z1")) {
    banner(paste("STAGE Z1 SKIPPED: Z0 did not reproduce production with the feature off, so",
                 "the ZINB edit is not behaviour-neutral and nothing in Z1 could be attributed",
                 "to zero inflation. The AR ladder is unaffected: its four rungs share one Stan",
                 "file, so it stays internally valid; only its comparison against the pre-edit",
                 "production run would be lost."))
    next
  }

  m <- S$model; run_tag <- S$tag; folder <- paste0(prefix[[m]], run_tag)
  done <- find_outdir(m, run_tag)
  if (isTRUE(RESUME) && !is.na(done) && file.exists(file.path(done, "port_total_Dungeness_Kept.csv"))) {
    banner(sprintf("SKIP (already complete): %s / %s / %s", sid, m, folder))
    r <- extract_run(sid, m, run_tag, done, NA_real_, TRUE)
    collected[[length(collected) + 1]] <- r; append_row(r)
  } else {
    rc <- resolve_cfg(sid); rc$model <- m; rc$run_tag <- run_tag; rc$run_weather <- FALSE
    banner(sprintf("STAGE %s  %s  |  %s  |  start %s", sid, folder, S$headline, format(Sys.time(), "%H:%M:%S")))
    cat(sprintf("  item %s | ar_force=%s | estimate_catch_zi=%s (%s) | pe_gear_ratio_arm=%s\n",
                S$item,
                if (is.null(rc$ar_force)) "NULL" else paste(utils::capture.output(str(rc$ar_force)), collapse = " "),
                rc$estimate_catch_zi %||% FALSE, paste(rc$catch_zi_populations %||% "-", collapse = "+"),
                rc$pe_gear_ratio_arm %||% "match_bss"))
    run_env <- new.env(parent = globalenv()); run_env$run_config <- rc
    t0 <- Sys.time()
    ok <- tryCatch({
      html <- rmarkdown::render(model_rmd[[m]], envir = run_env, quiet = FALSE)
      od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e) NA_character_)
      if (!is.na(od) && dir.exists(od) && normalizePath(dirname(html)) != normalizePath(od)) {
        if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE))) suppressWarnings(file.remove(html))
      }
      TRUE
    }, error = function(e) {
      message("*** STAGE ", sid, " / ", m, " FAILED: ", conditionMessage(e))
      od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e2) NA_character_)
      if (!is.na(od) && dir.exists(od)) {
        logs <- list.files(od, pattern = "^stan_console_.*\\.log$", full.names = TRUE)
        if (length(logs)) { last <- logs[which.max(file.mtime(logs))]
          message("    Stan console (", basename(last), "):")
          message(paste0("      ", utils::tail(readLines(last, warn = FALSE), 8), collapse = "\n")) }
        else message("    No stan_console_*.log in ", od, " -- failed before the first sampler call.")
      }
      FALSE })
    mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    r <- extract_run(sid, m, run_tag, find_outdir(m, run_tag), mins, ok)
    collected[[length(collected) + 1]] <- r
    tryCatch(append_row(r), error = function(e) message("  (summary row not written: ", conditionMessage(e), ")"))
    banner(sprintf("STAGE %s / %s -> %s (%s min)", sid, m, if (ok) "OK" else "FAILED", mins))
  }

  if (identical(sid, "Z0")) {
    nd <- find_outdir(m, run_tag)
    ex <- fit_exactness(nd %||% "", .here("05_output", REF$PROD$dir), NULL,
                        "all four fits vs 2026-08-31 production",
                        expect_delta = c("estimate_catch_zi", "catch_zi_populations",
                                         "zi_catch_prior_a", "zi_catch_prior_b",
                                         "pe_gear_ratio_arm"))
    gate_ok <- identical(ex$verdict, "PASS")
    cat("\n  GATE  ZINB edit, feature OFF: ", ex$verdict, "\n        ", ex$observed, "\n", sep = "")
    if (!gate_ok) message("*** The ZINB edit is not behaviour-neutral with the feature off. Z1 will be skipped.")
  }
}

# ---- THE SHORE AR LADDER ---------------------------------------------------
# The table this batch exists to produce. Four rungs on ONE component, differing in one
# thing. Every rung carries the quantities that separated the boat cells in the Stage 5
# batch, so the two are read the same way.
S <- if (length(collected)) do.call(rbind, collected) else NULL
g1 <- function(sid, col) { if (is.null(S)) return(NA); x <- S[[col]][S$stage == sid]; if (length(x)) x[1] else NA }

lad_row <- function(name, res, sid = NULL, from_dir = NULL) {
  r <- if (!is.null(from_dir)) {
    d <- .here("05_output", from_dir)
    if (dir.exists(d)) extract_run("Z0", "pooled", "archived", d, NA_real_, TRUE) else NULL
  } else if (!is.null(sid) && !is.null(S) && any(S$stage == sid)) S[S$stage == sid, ][1, ] else NULL
  if (is.null(r)) return(NULL)
  # n for the coverage z. Read from the run rather than hard-coded: it is a property of the
  # DATA and so constant across rungs, but a hard-coded 311 would silently mis-scale every z
  # if the season window or the effort file ever changed.
  n_gear <- {
    d <- if (!is.null(from_dir)) .here("05_output", from_dir) else find_outdir("pooled", STAGE_DEFS[[sid]]$tag)
    x <- if (!is.na(d %||% NA) && dir.exists(d %||% ""))
      rd(d, "ppc_byobs_shore_all_gear_Dungeness_Kept.csv") else NULL
    if (!is.null(x) && "data_type" %in% names(x)) sum(x$data_type == "gear") else NA_integer_
  }
  data.frame(rung = name, ar = res,
             shore_all_gear = r$shore_ag_BSS, reported = r$shore_ag_reported,
             port = r$port_BSS_catch, n_gear = n_gear,
             p_loo_frac = r$shore_ag_p_loo_frac, pareto_bad = r$shore_ag_bad_k,
             p_loo_gear = r$shore_ag_ploo_gear, elpd_gear = r$shore_ag_elpd_gear,
             elpd_catch = r$shore_ag_elpd_catch,
             cov50_gear = r$shore_ag_cov50_gear, cov50_catch = r$shore_ag_cov50_catch,
             cov50_gear_sd = (r$shore_ag_cov50_gear - 0.5) / sqrt(0.25 / max(n_gear, 1)),
             se_elpd_gear = r$shore_ag_se_elpd_gear, se_elpd_catch = r$shore_ag_se_elpd_catch,
             phi_E = r$shore_ag_phi_E, sigma_eps_E = r$shore_ag_sigma_eps_E,
             disp_scale = r$shore_ag_disp_scale, divergences = r$shore_ag_div,
             stringsAsFactors = FALSE)
}
# The ladder now comes from ONE run's ar_escalation_log.csv, which carries a row per rung
# with that rung's own estimate and interval. Z0 supplies the production daily rung fitted
# under the same Stan file, so the two are directly comparable.
LAD <- local({
  d <- find_outdir("pooled", STAGE_DEFS$A1$tag)
  x <- if (!is.na(d %||% NA)) rd(d, "ar_escalation_log.csv") else NULL
  if (is.null(x) || !"fit" %in% names(x)) return(lad_row("daily (production)", "daily", sid = "Z0"))
  x <- x[grepl("^shore_all_gear", x$fit), , drop = FALSE]
  if (!nrow(x)) return(lad_row("daily (production)", "daily", sid = "Z0"))
  data.frame(rung = paste0(x$ar_resolution, ifelse(x$selected, " (REPORTED)", "")),
             ar = x$ar_resolution, shore_all_gear = x$catch_median,
             reported = x$selected, port = NA_real_, n_gear = NA_integer_,
             p_loo_frac = NA_real_, pareto_bad = NA_real_, p_loo_gear = NA_real_,
             elpd_gear = NA_real_, elpd_catch = NA_real_,
             cov50_gear = NA_real_, cov50_catch = NA_real_, cov50_gear_sd = NA_real_,
             se_elpd_gear = NA_real_, se_elpd_catch = NA_real_,
             phi_E = NA_real_, sigma_eps_E = NA_real_, disp_scale = NA_real_,
             divergences = x$divergences,
             lo95 = x$catch_lo95, hi95 = x$catch_hi95,
             pi_width = x$pi_width, pi_width_rel = x$pi_width_rel,
             passed_gate = x$pass_convergence, stringsAsFactors = FALSE)
})
# The pre-edit production run is a FALLBACK rung, used only when Z0 has not run yet. It is
# merged, never substituted: substituting it truncated this file to a single row on every
# dry run, and the ladder is the main result of the batch.
if (is.null(LAD) || !nrow(LAD))
  LAD <- lad_row("daily (pre-edit production)", "daily", from_dir = REF$PROD$dir)
if (!is.null(LAD)) LAD <- merge_csv_by(LAD, lad_path, "rung")

# ---- VERDICTS --------------------------------------------------------------
if (!is.null(S) && any(S$stage == "G1")) {
  meth <- g1("G1", "boat_ag_method"); div <- g1("G1", "boat_ag_div")
  V1row("G1", "the gear boat all-gear fit samples now that it has per-fit settings",
    sprintf("method %s; divergences %s (was %s); tau_bar n_eff %s (was %s); port %s (was %s, S3 with the experiment override %s)",
            meth, fmt(div, 0), fmt(REF$GEARBAD$divergences, 0), fmt(g1("G1","boat_tau_neff"), 0),
            fmt(REF$GEARBAD$tau_neff, 0), fmt(g1("G1","port_BSS_catch"), 0),
            fmt(REF$GEARBAD$port, 0), fmt(REF$GEARFIX$port, 0)),
    "method_selected == BSS", if (identical(meth, "BSS")) "PASS" else "FAIL",
    paste("Until this run the gear-resolved port total read 51,385 because the boat all-gear",
          "fit was rejected and the component fell back to PE, a 27% understatement of the",
          "cross-check that validates the pooled headline. The 2026-08-30 Stage 5 batch had",
          "already demonstrated the fix and left it in an EXPERIMENT-only override, so",
          "production never received it. A PASS here means the fix is in the driver where it",
          "belongs and the two tracks can be compared again. A FAIL would mean the driver",
          "edit does not reach the fit, which the pre-flight tries to catch first."))
  V1row("G1", "cross-track agreement is restored",
    sprintf("gear port %s vs pooled production %s (%.2f%%); gear boat all-gear %s vs pooled %s; tau_bar %s vs %s",
            fmt(g1("G1","port_BSS_catch"), 0), fmt(REF$PROD$port, 0),
            100 * (g1("G1","port_BSS_catch") - REF$PROD$port) / REF$PROD$port,
            fmt(g1("G1","boat_ag_BSS"), 0), fmt(REF$PROD$boat_ag, 0),
            fmt(g1("G1","boat_tau_bar"), 4), fmt(REF$GEARFIX$tau_bar, 4)),
    "within about 1% at the port, as the 2026-08-26 shipped pair was (0.73%)",
    if (!is.na(g1("G1","port_BSS_catch")) &&
        abs(g1("G1","port_BSS_catch") - REF$PROD$port) / REF$PROD$port < 0.02) "PASS" else "REVIEW",
    paste("Two independently parameterized CPUE structures fitting the same data. This is",
          "the strongest internal check the project has, and it has been unavailable since",
          "the shared turnover was adopted because one side could not sample. Read tau_bar",
          "alongside the totals: on 2026-08-30 the two tracks agreed on it to 0.03%, which",
          "is what makes the shared turnover a property of the data rather than of the",
          "pooled parameterization."))
}

if (!is.null(S) && any(S$stage == "Z0")) {
  nd <- find_outdir("pooled", STAGE_DEFS$Z0$tag)
  ex <- fit_exactness(nd %||% "", .here("05_output", REF$PROD$dir), NULL,
                      "all four fits vs 2026-08-31 production",
                      expect_delta = c("estimate_catch_zi", "catch_zi_populations",
                                       "zi_catch_prior_a", "zi_catch_prior_b", "pe_gear_ratio_arm"))
  V1row("Z0", "the ZINB edit is behaviour-neutral with the feature off",
    ex$observed, "every shared parameter row identical at full precision", ex$verdict,
    paste("theta_C is declared vector<lower=0,upper=1>[zi_catch], so it is zero-size when",
          "the feature is off and the unconstrained parameter vector is unchanged. That is",
          "the same guard pattern shared_tau and f_lower use and it has held both times, but",
          "the project's own rule is to validate by run rather than by reasoning: a change",
          "that looks inference-neutral on paper can still perturb the sampler geometry.",
          "This is also the ladder's daily rung, so a FAIL costs the ladder its comparison",
          "against the pre-edit production run as well as costing Z1 its interpretation."))
}

# ar_force LEAK CONTROL. On 2026-08-28 an ar_force keyed per POPULATION silently forced the
# boat's ALL-GEAR sub-season as well as its pot closure, moving that component 25,883 ->
# 28,893 and costing 3,025 crab in a result reported as a pot-closure finding. The nested
# form is fixed and V4 re-confirmed it, but this ladder forces a sub-season on the SHORE for
# the first time, so the same check is worth its zero cost here. Every fit is an independent
# Stan run with the same seed, so the three untouched fits must be BIT-IDENTICAL to Z0 in
# every rung, not merely close.
# SCOPE CONTROL. `ar_escalate = list(shore = "all_gear")` should reach exactly one fit. The
# other three are fitted once each and must be bit-identical to Z0. This is the same class
# of check that caught the 2026-08-28 ar_force leak, which cost 3,025 crab, and it is
# cheaper here because the scope is resolved per fit rather than per population.
if (!is.null(S) && any(S$stage == "A1")) {
  nd <- find_outdir("pooled", STAGE_DEFS$A1$tag); z0 <- find_outdir("pooled", STAGE_DEFS$Z0$tag)
  if (!is.na(nd %||% NA) && !is.na(z0 %||% NA)) {
    .esc_keys <- c("ar_escalate", "ar_escalate_stop", "ar_escalate_select", "ar_escalate_respect_cap")
    ex_pc <- fit_exactness(nd, z0, "shore_ring_net_only", "shore POT CLOSURE fits vs Z0",
                           expect_delta = .esc_keys)
    ex_bt <- fit_exactness(nd, z0, "private_boat", "BOAT fits vs Z0", expect_delta = .esc_keys)
    V1row("A1", "the escalation scope reaches shore all-gear and nothing else",
      sprintf("%s | %s", ex_pc$observed, ex_bt$observed),
      "the other three fits bit-identical to Z0",
      if (identical(ex_pc$verdict, "PASS") && identical(ex_bt$verdict, "PASS")) "PASS" else "FAIL",
      paste("`ar_escalate` accepts a per-population and per-sub-season scope as of",
            "2026-09-04, and this is its first use. A scope that reaches further than",
            "intended produces a number that looks like a finding: that is exactly what",
            "happened on 2026-08-28 when an ar_force keyed per population silently forced a",
            "second sub-season and moved a component 3,025 crab. The boat half also protects",
            "the headline, since the boat is about 43% of the port total and nothing in this",
            "ladder should move it."))
  }
}

if (!is.null(LAD) && nrow(LAD) > 1) {
  # A rung whose gate failed reports PE, and its BSS_catch column is then a number nobody
  # reported. Say so rather than ranking it alongside the others.
  .rej <- LAD$rung[!is.na(LAD$reported) & !LAD$reported]
  if (length(.rej))
    V1row("A1-A3", "every rung actually reported its BSS",
      sprintf("rung(s) whose shore all-gear fit was REJECTED and fell back to PE: %s",
              paste(.rej, collapse = ", ")),
      "all rungs report BSS", "REVIEW",
      paste("A rejected rung's BSS_catch is a value the run did not report, so it cannot be",
            "ranked against the others on catch. Its ADEQUACY columns are still readable and",
            "still informative: a rung that samples badly enough to fail the gate is telling",
            "you something about the resolution, which is the point of the ladder."))
  best_cov <- LAD$rung[which.min(abs(LAD$cov50_gear - 0.5))]
  best_elpd <- LAD$rung[which.max(LAD$elpd_gear)]
  V1row("A1-A3", "where does the shore all-gear effort process actually belong?",
    paste(sprintf("%s: p_loo %s%% of n_obs, k>0.7 %s, gear cov50 %s (%+.1f SD), gear elpd %s, catch elpd %s, shore catch %s",
                  LAD$rung, fmt(100 * LAD$p_loo_frac, 1), fmt(LAD$pareto_bad, 0),
                  fmt(LAD$cov50_gear, 3), LAD$cov50_gear_sd, fmt(LAD$elpd_gear, 1),
                  fmt(LAD$elpd_catch, 1), fmt(LAD$shore_all_gear, 0)), collapse = " | "),
    "best calibrated rung, read together with p_loo and the catch stream",
    sprintf("BEST COVERAGE: %s | BEST GEAR elpd: %s", best_cov, best_elpd),
    paste("THE QUESTION THIS BATCH EXISTS FOR. Production fits this component at daily and",
          "it is the worst-behaved fit in the production run: p_loo 35.2% of n_obs, 41",
          "Pareto k above 0.7, coverage_50 0.701 at +7.1 sampling SDs. The gear track fits",
          "it at monthly and gets coverage_50 0.035, -16.4 SDs, so both ends are bad. Read",
          "this ladder the way the Stage 5 boat 2x2 was read, and in this order: (1)",
          "COVERAGE first, because it is the statistic a reader can interpret without",
          "knowing what a PIT is and the one that separated the boat cells; (2) p_loo as a",
          "fraction of n_obs and the Pareto k count, because a rung that buys elpd by",
          "spending effective parameters at about one nat each is buying noise; (3) the",
          "CATCH stream's elpd, because the boat's daily AR improved the effort stream while",
          "making the catch stream worse, and catch is what the harvest estimate is made of.",
          "A rung that wins on all three is the answer. If the winner is not daily, the",
          "ar_max_resolution entry for shore all-gear should change and the shore component",
          "of the published estimate moves with it, so this is a number-moving decision and",
          "needs sign-off rather than a silent config edit."))
}

if (!is.null(S) && any(S$stage == "Z1")) {
  th <- g1("Z1","shore_ag_theta_C")
  e0 <- g1("Z0","shore_ag_elpd_catch"); e1 <- g1("Z1","shore_ag_elpd_catch")
  V1row("Z1", "does a structural-zero component earn its parameter on the shore catch stream?",
    sprintf("theta_C = %s [%s, %s]; shore pot closure theta_C = %s; catch elpd %s -> %s (%+.1f nats for 1 parameter)",
            fmt(th, 4), fmt(g1("Z1","shore_ag_theta_lo95"), 4), fmt(g1("Z1","shore_ag_theta_hi95"), 4),
            fmt(g1("Z1","shore_pc_theta_C"), 4), fmt(e0, 1), fmt(e1, 1), e1 - e0),
    "an elpd gain of at least 2 SE, and a theta_C interval excluding zero",
    { .se <- max(g1("Z1","shore_ag_se_elpd_catch") %||% NA, g1("Z0","shore_ag_se_elpd_catch") %||% NA, na.rm = TRUE)
      if (is.na(th) || is.na(e1) || is.na(e0)) "REVIEW"
      else if (!is.finite(.se) || .se <= 0) "REVIEW (no elpd SE)"
      else if ((e1 - e0) > 2 * .se && isTRUE(g1("Z1","shore_ag_theta_lo95") > 0.005)) "ADOPTABLE"
      else if ((e1 - e0) > .se) "MARGINAL" else "NOT WORTH IT" },
    paste("The 2026-09-01 zero bin put the shore all-gear catch stream at 676 observed zeros",
          "against 605.4 expected, z = +3.8 on n = 1,649, with the shore pot closure agreeing",
          "at z = +2.9 and every boat stream inside |z| = 2.3. So the target was chosen from",
          "data rather than from taste, and the available gain is bounded by roughly 70",
          "observations out of 1,649. THE THRESHOLD IS IN SE UNITS, not nats: loo_summary",
          "carries se_elpd_loo and an elpd difference is only evidence when it is large",
          "relative to its own standard error, which for a stream of this size is typically",
          "tens of nats. A gain above 2 SE is a clearly good trade for one parameter; between",
          "1 and 2 SE is a judgement call; below 1 SE is fitting noise. The shore POT CLOSURE",
          "fit is the internal replicate: an",
          "independent fit on different data that should land near the same theta_C if the",
          "effect is real."))
  V1row("Z1", "the zero bin closes on the treated stream",
    sprintf("shore all-gear zeros: %s observed vs %s expected (z = %s), was %s vs %s (z = %s)",
            fmt(g1("Z1","shore_ag_zero_obs"), 0), fmt(g1("Z1","shore_ag_zero_exp"), 1),
            fmt(g1("Z1","shore_ag_zero_z"), 1), fmt(REF$PROD$shore_ag_obs_zeros, 0),
            fmt(REF$PROD$shore_ag_exp_zeros, 1), fmt(REF$PROD$shore_ag_zero_z, 1)),
    "|z| under about 2", if (!is.na(g1("Z1","shore_ag_zero_z")) && abs(g1("Z1","shore_ag_zero_z")) < 2) "PASS" else "REVIEW",
    paste("The direct test that the feature does what it was added for. A theta_C with a",
          "good elpd gain but an unchanged zero bin would mean the parameter is absorbing",
          "something else. Note this is an IN-SAMPLE check on the quantity the parameter was",
          "fitted to, so it cannot by itself justify adoption; it is the necessary condition,",
          "and the elpd criterion above is the sufficient one."))
  ex <- fit_exactness(find_outdir("pooled", STAGE_DEFS$Z1$tag) %||% "",
                      find_outdir("pooled", STAGE_DEFS$Z0$tag) %||% "", "private_boat",
                      "boat fits vs Z0", expect_delta = c("estimate_catch_zi"))
  V1row("Z1", "the boat fits are an untouched negative control",
    ex$observed, "boat fits bit-identical to Z0", ex$verdict,
    paste("zi_catch is set per FIT from catch_zi_populations, so the boat fits in this very",
          "run carry the feature off. That makes the control free: no extra run, and the",
          "same RNG stream and the same everything else. A FAIL means the per-fit scoping",
          "leaks, which would make the shore result unattributable AND would have quietly",
          "changed the headline boat number."))
  V1row("Z1", "what it does to the reported shore total",
    sprintf("shore all-gear %s (Z0 %s, %+.1f%%); port %s (Z0 %s)",
            fmt(g1("Z1","shore_ag_BSS"), 0), fmt(g1("Z0","shore_ag_BSS"), 0),
            100 * (g1("Z1","shore_ag_BSS") - g1("Z0","shore_ag_BSS")) / g1("Z0","shore_ag_BSS"),
            fmt(g1("Z1","port_BSS_catch"), 0), fmt(g1("Z0","port_BSS_catch"), 0)),
    "informational, but READ THE WHY", "INFO",
    paste("The season total is scaled by (1 - theta_C) in generated quantities, deliberately.",
          "lambda_C is fitted to the NON-inflated component and rises to absorb the zeros",
          "theta_C removed, so reporting lambda_C * effort unscaled would inflate the total",
          "by 1/(1 - theta_C) purely as an artefact of turning the feature on. If this shows",
          "the shore total moving by roughly theta_C in either direction, check that scaling",
          "before believing it: a total that moves UP by about theta_C is the signature of",
          "the scaling having been dropped somewhere."))
}

VD <- if (length(V)) do.call(rbind, V) else NULL
# Merge by (stage, criterion): a dry run that scores only P1 must not erase the fitted
# stages' criteria. See the note at merge_csv_by().
if (!is.null(VD)) VD <- merge_csv_by(VD, ver_path, c("stage", "criterion"))

banner("BATCH SUMMARY")
if (!is.null(S)) {
  print(S[, c("stage","item","model","ok","minutes","port_BSS_catch","shore_ag_BSS",
              "shore_ag_ar","shore_ag_p_loo_frac","shore_ag_cov50_gear","fits_passed")], row.names = FALSE)
  cat("\n  full summary -> ", sum_path, "\n", sep = "")
}
if (!is.null(LAD) && nrow(LAD) > 1) {
  rule(); cat("  SHORE ALL-GEAR AR LADDER  (nominal coverage_50 = 0.500)\n\n")
  cat(sprintf("  %-20s %9s %7s %9s %10s %11s %11s\n", "rung", "catch", "p_loo", "k>0.7",
              "gear cov50", "gear elpd", "catch elpd"))
  for (i in seq_len(nrow(LAD)))
    cat(sprintf("  %-20s %9s %6.1f%% %9s %10s %11s %11s\n", LAD$rung[i],
                fmt(LAD$shore_all_gear[i], 0), 100 * LAD$p_loo_frac[i], fmt(LAD$pareto_bad[i], 0),
                sprintf("%.3f", LAD$cov50_gear[i]), fmt(LAD$elpd_gear[i], 1), fmt(LAD$elpd_catch[i], 1)))
  cat("\n  ladder -> ", lad_path, "\n", sep = "")
}
if (!is.null(VD)) {
  rule()
  for (i in seq_len(nrow(VD)))
    cat(sprintf("[%s] %-28s %s\n        %s\n", VD$stage[i], VD$verdict[i], VD$criterion[i], VD$observed[i]))
  cat("\n  verdicts -> ", ver_path, "\n", sep = "")
}
if (isTRUE(DRY_RUN)) {
  rule()
  banner("DRY RUN: P1 ran for real; the fitted stages were resolved and NOT run")
  for (sid in STAGES) {
    Sd <- STAGE_DEFS[[sid]]; if (identical(Sd$kind, "pe_only")) next
    cfg <- resolve_cfg(sid); ks <- names(Sd$delta %||% list())
    cat(sprintf("  %-3s %-24s %-14s %s\n", sid, Sd$tag, Sd$model, Sd$headline))
    if (!length(ks)) cat("      delta: none (production as shipped)\n")
    for (k in ks) cat(sprintf("      delta: %-24s %s\n", k,
      paste(utils::capture.output(str(cfg[[k]])), collapse = " ")))
  }
  cat("\n  Read the P1 result above, then set DRY_RUN <- FALSE.\n\n")
} else banner(sprintf("DONE  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
