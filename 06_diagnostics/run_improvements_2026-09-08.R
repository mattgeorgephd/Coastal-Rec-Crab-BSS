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
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                    # TRUE: pre-flight + the R0 desk rung, nothing fitted. START HERE.
STAGES  <- c("R0", "R1", "R2", "R3a", "R3", "R4", "R5")
RESUME  <- TRUE                    # skip a rung whose output folder already exists
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
V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = observed,
                                    threshold = threshold, verdict = verdict, why = why,
                                    stringsAsFactors = FALSE)

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
WINDOW <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15", season_filter = "2024-25",
               pot_closures = NULL, census_windows = NULL, run_weather = FALSE)
D_R1  <- list(tau_boat_prior_mu = 1.2, tau_boat_prior_sigma = 0.3, shared_tau_sigma = NULL,
              crab_fraction_strata = "none", crab_fraction_source = "ie", crab_fraction_dynamic = FALSE,
              tau_shore_prior_mu = 1.7, tau_shore_prior_sigma = 0.3, census_uncertainty = "none")
D_R2  <- modifyList(D_R1,  list(tau_boat_prior_mu = "calibration", tau_boat_prior_sigma = 0.5, shared_tau_sigma = 0.15))
D_R3a <- modifyList(D_R2,  list(crab_fraction_strata = "month", crab_fraction_source = "both"))
D_R3  <- modifyList(D_R3a, list(crab_fraction_dynamic = TRUE))
D_R4  <- modifyList(D_R3,  list(tau_shore_prior_mu = "derived", tau_shore_prior_sigma = "derived"))

STAGE_DEFS <- list(
  R0  = list(id = "R0",  model = "desk",   tag = "IMP-R0-desk",        delta = D_R3,
             headline = "desk: PE on the calibration turnover, gear bootstrap, census split, the f data"),
  R1  = list(id = "R1",  model = "pooled", tag = "IMP-R1-filters",     delta = D_R1,
             headline = "pre-patch configuration + item 8 (unit-aware fishing-time filters)"),
  R2  = list(id = "R2",  model = "pooled", tag = "IMP-R2-tau-calib",   delta = D_R2,
             headline = "+ item 3: boat turnover prior from the OSP/trailer calibration"),
  R3a = list(id = "R3a", model = "pooled", tag = "IMP-R3a-f-monthly",  delta = D_R3a,
             headline = "+ item 1A: monthly f from the sampler contacts (legacy construction)"),
  R3  = list(id = "R3",  model = "pooled", tag = "IMP-R3-f-dynamic",   delta = D_R3,
             headline = "+ item 1B: the dynamic f (the shipped configuration)"),
  R4  = list(id = "R4",  model = "pooled", tag = "IMP-R4-shore-tau",   delta = D_R4,
             headline = "+ item 2: shore turnover derived from the I/E time column (the shipped configuration)"),
  R5  = list(id = "R5",  model = "gear_resolved", tag = "IMP-R5-gear-crosscheck",
             delta = if (identical(GEAR_FOLLOWS, "R4")) D_R4 else D_R3,
             headline = sprintf("gear-resolved cross-check on the %s configuration", GEAR_FOLLOWS)))
resolve_cfg <- function(sid) {
  cfg <- modifyList(BASE, WINDOW, keep.null = TRUE)
  modifyList(cfg, STAGE_DEFS[[sid]]$delta %||% list(), keep.null = TRUE)
}
find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}
prev_pooled <- function(sid) {
  # the pooled rung this one is compared to (R3a is optional, so R3 falls back to R2)
  order <- c("R1", "R2", "R3a", "R3", "R4")
  i <- match(sid, order); if (is.na(i) || i == 1) return(NA_character_)
  for (p in rev(order[seq_len(i - 1)])) if (p %in% STAGES) return(p)
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
  say(identical(BASE$census_uncertainty, "none") && identical(BASE$census_expansion, "none"),
      "the census is exact over the tally days and stays a constant (census_expansion none, census_uncertainty none)")
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
  banner("R0  desk: the PE and reporting changes, and the data the fitted rungs read")
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
    cc_tally <- q(estimate_comm_charter(dwg, modifyList(p, list(charter_frame = "tally"))))
    V1row("R0", "the census is the exact sum over the tally days; unsampled days had no operation (item 4, settled 2026-09-09); the charter roster is the charter frame (2026-09-10)",
          sprintf(paste0("census_expansion = '%s': %s crab observed on %d tally days; %d unsampled calendar days %s; SE %s (%.1f%%, the per-vessel mean); ",
                         "total %s with charter_frame = '%s' (%s crab on %d charter-roster days without a tally; %s under the tally frame; baseline %s under the day-type expansion)"),
                  cc$census_expansion %||% "none", fmt(cc$observed_dung, 0), sum(cc$daily_full$observed),
                  sum(!cc$daily_full$observed), if (identical(cc$census_expansion, "none")) "at zero" else sprintf("imputed at %s", fmt(cc$imputed_dung, 0)),
                  fmt(cc$Dungeness_Kept_se, 0), 100 * cc$Dungeness_Kept_se / max(cc$Dungeness_Kept, 1), fmt(cc$Dungeness_Kept, 0),
                  cc$charter_frame %||% "tally", fmt(cc$charter_roster_dung %||% 0, 0), cc$n_roster_only_days %||% 0L, fmt(cc_tally$Dungeness_Kept, 0), fmt(REF$A1$census, 0)),
          "the total is the tally-day sum; the SE stays out of the port interval under census_uncertainty = none",
          if (identical(cc$census_expansion, "none") && isTRUE(all.equal(cc$imputed_dung, 0))) "PASS" else "REVIEW",
          paste("Samplers are scheduled on the days the charter and commercial (recreational) vessels are",
                "confirmed to be operating, so a window day without a tally is a day with no fishing, not a",
                "missed count. The 2026-09-08 day-type expansion (11,753 on 2024-25) filled those days with",
                "the sampled days' mean and is kept as census_expansion = 'day_type' for reproduction only."))
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
    cat("  RESUME: output already present at", basename(existing), "- skipping the fit.\n")
    return(existing)
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
  V1row("R1", "item 8 alone (plus the 2026-09-09 workbook): the fitted components against the baseline",
        sprintf("shore pc %s (%+.1f%%), shore ag %s (%+.1f%%), boat pc %s (%+.1f%%), boat ag %s (%+.1f%%); census %s (baseline %s: exact over the tally days by design)",
                fmt(sp, 0), .pct(sp, REF$A1$shore_pc), fmt(sa, 0), .pct(sa, REF$A1$shore_ag),
                fmt(bp, 0), .pct(bp, REF$A1$boat_pc), fmt(ba, 0), .pct(ba, REF$A1$boat_ag),
                fmt(.comp(dir, "comm_charter (census)", "PE_catch"), 0), fmt(REF$A1$census, 0)),
        "every fitted component within a few percent; the boat a few percent LOWER (zero-catch trips restored)",
        if (all(abs(c(.pct(sp, REF$A1$shore_pc), .pct(sa, REF$A1$shore_ag))) < 5, na.rm = TRUE) &&
              isTRUE(.pct(ba, REF$A1$boat_ag) < 2)) "PASS" else "REVIEW",
        paste("The hours filter dropped 22 boat trips (14 with recorded catch) and ~140 shore rows whose",
              "time field was blank or under 0.5 h; under the deployment unit a set pot with zero time",
              "is real effort with real catch. Restoring them lowers the boat CPUE slightly (more zeros)",
              "and barely moves the shore. The rebuilt workbook (2026-09-09) also applies the tampered-gear",
              "filter (13 Grays Harbor interviews in 2024-25) and a few source corrections. The census is",
              "not compared: it is now the exact sum over the tally days (-3,869 by design)."))
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
  V1row("R2", "the shore did not move (the boat prior cannot reach it)", fe$observed,
        "shore fits bit-identical to R1", fe$verdict,
        "A FAIL means the turnover prior leaked into the shore fits, which it has no path to do.")
  tb <- .full_row(dir, "private_boat_all_gear", "tau_bar_out")
  ba <- .comp(dir, "private_boat (All gear)"); ba0 <- .comp(prev %||% "", "private_boat (All gear)")
  V1row("R2", "tau_bar lands on the calibration and the boat rises with it",
        sprintf("tau_bar %s [%s, %s] (R1 %s); boat all-gear %s vs R1 %s (%+.1f%%)",
                fmt(.num1(tb$mean), 3), fmt(.num1(tb[["2.5%"]]), 3), fmt(.num1(tb[["97.5%"]]), 3),
                fmt(.num1(.full_row(prev %||% "", "private_boat_all_gear", "tau_bar_out")$mean), 3),
                fmt(ba, 0), fmt(ba0, 0), .pct(ba, ba0)),
        "tau_bar 2.8-3.1; boat all-gear +10 to +18%",
        if (isTRUE(.num1(tb$mean) > 2.75) && isTRUE(.num1(tb$mean) < 3.2) && isTRUE(.pct(ba, ba0) > 8) && isTRUE(.pct(ba, ba0) < 20))
          "PASS" else "REVIEW",
        paste("Effort = lambda_E x tau x f, so the boat moves in proportion to tau_bar. Likelihood-only tau",
              "was 2.97 against the shrunk posterior 2.60, i.e. +14%; the OSP series identifies tau_bar",
              "and the prior should no longer fight it. A move outside 8-20% means something else changed."))
  pt <- rd(dir, "ppc_calibration_private_boat_all_gear_Dungeness_Kept.csv")
  if (!is.null(pt))
    V1row("R2", "boat PPC still calibrated after the prior move",
          paste(utils::capture.output(print(utils::head(pt, 6))), collapse = " | "),
          "PIT means near 0.5, coverage near nominal", "READ",
          "The prior move should improve, not worsen, the OSP stream's fit: read the OSP/trailer rows.")
}

verdict_R3a <- function(dir, prev) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  fe <- fit_exactness(dir, prev %||% "", pat = "shore", what = "shore fits vs the previous rung",
                      expect_delta = c("crab_fraction_strata", "crab_fraction_source", "crab_fraction_rows", "run_tag", "model",
                                       "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1row("R3a", "the shore did not move (f is boat-only)", fe$observed, "shore fits bit-identical", fe$verdict,
        "f enters the boat generated quantities only; the shore has no f.")
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
  V1row("R3", "the shore did not move (f is boat-only)", fe$observed, "shore fits bit-identical", fe$verdict,
        "f enters the boat generated quantities only; the shore has no f.")
  # THE FACTORIZATION PROOF, numerically: every non-f boat parameter within MC error.
  fa <- fit_agreement(dir, prev %||% "", pat = "private_boat",
                      exclude = "^(f_crab|f_lower|f_theta|f_lower_param|z_f|sigma_f|cfi_kappa|combo_c|eta_f|osp_f_kappa|E\\[|E_sum|C\\[|C_sum|C_expected|lambda_Ctot|log_lik|lp__)",
                      what = "boat non-f parameters vs the previous rung")
  V1row("R3", "the dynamic f leaves effort and CPUE untouched (factorization)", fa$observed,
        "max |z| under 5 and under 1% of rows above 3", fa$verdict,
        paste("Bit-identity is impossible here: the parameter vector changed (z_f, sigma_f, cfi_kappa),",
              "so every HMC trajectory differs. What the design guarantees is agreement in DISTRIBUTION,",
              "and this is that test. A FAIL means an f term reached an effort or CPUE likelihood."))
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

verdict_R4 <- function(dir, prev) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  fe <- fit_exactness(dir, prev %||% "", pat = "private_boat", what = "boat fits vs R3",
                      expect_delta = c("tau_shore_prior_mu", "tau_shore_prior_sigma", "tau_shore_prior_source", "shared_tau_min_obs",
                                       "run_tag", "model", "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1row("R4", "the boat did not move (the shore prior cannot reach it)", fe$observed, "boat fits bit-identical to R3", fe$verdict,
        "A FAIL means the shore turnover leaked into the boat fits, which it has no path to do.")
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
              "worth reading: the I/E stream is the only shore observation that can."))
  s1 <- .full_row(dir, "shore_all_gear", "sigma_IE_out"); s0 <- .full_row(prev %||% "", "shore_all_gear", "sigma_IE_out")
  V1row("R4", "the I/E observation scale did not worsen", sprintf("sigma_IE %s (R3 %s)", fmt(.num1(s1$mean), 3), fmt(.num1(s0$mean), 3)),
        "not larger than R3's", if (isTRUE(.num1(s1$mean) <= .num1(s0$mean) * 1.05)) "PASS" else "REVIEW",
        "If the I/E trips disagree with the derived turnover, sigma_IE absorbs the disagreement; it should not grow.")
  V1row("R4", "the port total with the derived shore turnover",
        sprintf("%s [%s, %s] (R3 %s; baseline %s)", fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_median), 0),
                fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_lo95), 0), fmt(.num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_hi95), 0),
                fmt(.num1(.port_row(rd(prev %||% "", "port_total_Dungeness_Kept.csv"))$BSS_median), 0), fmt(REF$A1$port, 0)),
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
banner(sprintf("IMPROVEMENT LADDER  2026-09-08   DRY_RUN=%s  stages: %s  (R5 follows %s)", DRY_RUN,
               paste(STAGES, collapse = ", "), GEAR_FOLLOWS))
preflight()

dirs <- list()
for (sid in STAGES) dirs[[sid]] <- run_stage(sid)

.safe <- function(sid, expr) tryCatch(force(expr), error = function(e) {
  cat(sprintf("\n  *** VERDICT BLOCK %s FAILED: %s\n", sid, conditionMessage(e)))
  cat("      The run itself is intact; re-run with RESUME = TRUE after fixing the block.\n")
  V1row(sid, "verdict block did not complete", conditionMessage(e), "the block runs to completion", "ERROR",
        "The fits are on disk and unaffected: this is a defect in the code that READS them; re-run with RESUME = TRUE.")
  invisible(NULL)
})
.dir_of <- function(sid) if (!is.na(sid %||% NA) && sid %in% names(dirs)) dirs[[sid]] else NA_character_
for (sid in intersect(c("R1", "R2", "R3a", "R3", "R4"), STAGES)) .safe(sid, ladder_row(sid, dirs[[sid]]))
if ("R1"  %in% STAGES) .safe("R1",  verdict_R1(dirs$R1))
if ("R2"  %in% STAGES) .safe("R2",  verdict_R2(dirs$R2, .dir_of(prev_pooled("R2"))))
if ("R3a" %in% STAGES) .safe("R3a", verdict_R3a(dirs$R3a, .dir_of(prev_pooled("R3a"))))
if ("R3"  %in% STAGES) .safe("R3",  verdict_R3(dirs$R3, .dir_of(prev_pooled("R3"))))
if ("R4"  %in% STAGES) .safe("R4",  verdict_R4(dirs$R4, .dir_of(prev_pooled("R4"))))
if ("R5"  %in% STAGES) .safe("R5",  verdict_R5(dirs$R5, .dir_of(GEAR_FOLLOWS)))

if (length(LAD)) {
  lp <- .here("05_output", "improvements_2026-09-08_ladder.csv")
  merge_csv_by(do.call(rbind, LAD), lp, "rung")
  banner("LADDER")
  print(do.call(rbind, LAD), row.names = FALSE)
  cat("\n  written to", lp, "\n")
}
if (length(V)) {
  vp <- .here("05_output", "improvements_2026-09-08_verdicts.csv")
  merge_csv_by(do.call(rbind, V), vp, c("stage", "criterion"))
  banner("VERDICTS")
  for (r in seq_len(length(V))) with(V[[r]], cat(sprintf("  [%s] %-14s %s\n     obs: %s\n",
    stage, verdict, criterion, observed)))
  cat("\n  written to", vp, "\n")
}
