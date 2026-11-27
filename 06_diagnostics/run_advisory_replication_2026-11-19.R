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
# run_advisory_replication_2026-11-19.R
#
# THE BOAT ADVISORY TERM ON A SECOND SEASON (CHANGE_REGISTER B80, for D31 and D32;
# VALIDATION_CAMPAIGN Section 1z.21).
#
# WHY THIS RUN EXISTS
#   The NWS Small Craft Advisory flag on the private-boat all-gear effort process (A30) was
#   adopted on 2026-09-27 on the evidence of ONE season: identified at eight standard errors
#   on 2024-25, and better on that season's held-out weeks where the OSP counts are (the block
#   CV, D31). Two things were left open. The term has never been tested on held-out weeks of a
#   season it was not chosen on. And the winter, where the term moves the boat estimate most,
#   is where no held-out test on 2024-25 could see; the four-season desk screen says the
#   winter effect is about half the rest of the season's on the log scale, and rung M6 could
#   not place a winter coefficient on 2024-25 alone (D32).
#
#   2025-26 is a second season. Its boat all-gear fit has no OSP day yet (D18) and its
#   turnover is its prior (A41), but the advisory term sits on the PRESENCE process, which
#   the trailer counts observe, so the refused turnover does not touch it. Nothing rendered
#   here is an estimate of 2025-26 (D18); it is a covariate test.
#
# THE ARMS. Each is the pooled driver rendered on run_config.R AS SHIPPED, with the season's
# window and one delta on the marine hazard keys. Three per season:
#   const   the shipped term: nws_sca_any, one coefficient, the boat all-gear fit
#   off     marine_hazard_mode = "off": the model without the term
#   split   nws_sca_any_winter (December to February) + nws_sca_any_rest, the ladder's M6
#   S26-*   2025-26 (the window block of NEW_SEASON_GUIDE.md section 7.1)
#   S25-*   2024-25 (run_config.R as shipped; S25-const is the authoritative configuration
#           rendered on this container, reported against B70 and never an authoritative run)
# The shore fits do not carry the term, so with BSS_FIT_CACHE_DIR set they are sampled once
# per season and loaded from the cache in the other two arms.
#
# THE RULE, STATED HERE AND IN CHANGE_REGISTER B80 BEFORE ANY ARM WAS RENDERED OR READ.
#   Known before the rule and therefore reported, not judged: the committed coefficients
#   (B70 -1.146, SD 0.155; B76 -1.321, SD 0.186; 0.72 SE apart); the 2026-09-27 desk screen,
#   which had read the 2025-26 trailer days (winter 0.44 [0.24, 0.82], rest 0.20 [0.13, 0.28]
#   on that season); M6 on 2024-25 before A33 (winter -0.59 [-1.24, +0.06]).
#   R1  ELIGIBILITY. An arm is read only if every fit in it passes the convergence gate and
#       its boat all-gear fit ran at the AR resolution of the season's const arm. A pair with
#       an ineligible arm is REVIEW.
#   R2  D31 OUT OF SAMPLE, THE TRANSFER (desk_advisory_transfer_2026-11-19.R). Target 2025-26,
#       source 2024-25 alone: TC - M0 over the scored weeks above +2 paired SE is REPLICATED
#       OUT OF SAMPLE; below -2 is FAILED; between is no evidence. The pooled source and the
#       mirror (target 2024-25) are read the same way and reported beside it. F - TC above +2
#       SE says the target's own coefficient is detectably different from the carried one.
#   R3  D31 ON THE FITTED MODEL. Per season, const against off on the boat all-gear fit's
#       held-out weeks: the JOINT row where the fit has two effort streams, the trailer row
#       where it has one. Evaluable if at least 70% of the weeks are reliable (Pareto k <= 0.7
#       in both arms) and at least 8 remain. Above +2 paired SE is PASS, below -2 is FAIL,
#       between is REVIEW. STATED IN ADVANCE: on 2024-25 the trailer stream alone scored +7.9
#       at 1.2 SE and the OSP stream carried the evidence; the 2025-26 all-gear fit has only
#       the trailer stream, so REVIEW is the likely outcome there and would not be evidence
#       against the term. R2 is the instrument with power for that season.
#   R3a ADDENDUM, WRITTEN 2026-11-19 AFTER S26-const ALONE HAD BEEN READ AND BEFORE ANY OTHER ARM
#       FINISHED (so before any pair of arms existed). R3's instrument is not evaluable under the
#       current model, and the expectation stated for it below was wrong. Since the per-day
#       effort effect (A39) the PSIS leave-one-week-out is reliable on 12 of 42 joint weeks of
#       the authoritative run B70 (39 of 42 on B65, 38 of 42 on B67, the renders before it) and
#       on 12 of 41 trailer weeks of S26-const, whose table is B76's. A held-out week's trailer
#       counts are what identify those days' deviations uE, so conditional on uE the importance
#       ratios are heavy-tailed. The expectation of PASS on 2024-25 was written from the
#       2026-09-26 block CV without opening B70's committed table (a C row; D54).
#       R3 STAYS AS WRITTEN and will read REVIEW, not evaluable, on both seasons. R3a is the
#       SAME comparison (const minus off per season; the joint row where the fit has two effort
#       streams, the trailer row where it has one), the SAME evaluability (70% of weeks reliable
#       in both arms, at least 8) and the SAME thresholds (above +2 paired SE PASS, below -2
#       FAIL, between REVIEW), on the block CV with each held-out trailer day's deviation
#       INTEGRATED OUT of both the importance ratio and the score:
#         p(y_d | theta) = integral of prod_j NB2(y_dj | lambda_d / R_G_boat
#                          * exp(sigma_uE z - sigma_uE^2 / 2), r_E) N(z | 0, 1) dz
#       by quadrature on a fixed grid in z (the integrated form of Vehtari et al. 2016 and
#       Merkle, Furr and Rabe-Hesketh 2019); the OSP likelihood carries no day deviation and
#       enters unchanged. Preconditions, checked before a comparison is opened: (i) the
#       conditional log-likelihood rebuilt from the cached fit's draws reproduces the fit's own
#       log_lik_trailer to 1e-6; (ii) doubling the grid moves no week's leave-out elpd by more
#       than 0.01 nats. Wherever R4 (b) reads a block CV row it is read on both tables. R3a is
#       reported BESIDE R3, never in its place. EXPECTED: evaluable on both seasons; PASS on
#       2024-25, where the OSP stream is; REVIEW on 2025-26. Known when this was written:
#       S26-const's coefficient, totals and block table (equal to B76's committed ones) and
#       the transfer test; no off or split arm (desk_advisory_block_cv_2026-11-19.R).
#   R4  D32 PER SEASON (the ladder's rule 9, restated). On the split arm's boat all-gear fit:
#       (a) the winter coefficient's 95% interval covers 0: this season cannot place it;
#       (b) both identified and the 95% interval of winter - rest (independent normals)
#           excludes 0: the two differ on this season, and the split earns its place there
#           only if the catch stream is unmoved against const (paired elpd within 2 SE) and
#           the block CV row, split against const, is at or above -2 paired SE;
#       (c) both identified and that interval covers 0: no difference detectable.
#   R5  D32 OVER THE TWO SEASONS. d = winter - rest per season with SD sqrt(sd_w^2 + sd_r^2);
#       the precision-weighted mean of d over the eligible split arms and its z. |z| >= 2:
#       the two seasons together establish a difference (the sign says which effect is the
#       larger); 1 <= |z| < 2: suggestive; |z| < 1: none detectable. If the two seasons' d
#       differ from each other at |z| >= 2 the combination is reported with that warning and
#       is not read as established.
#   R6  D32 OUT OF SAMPLE (the transfer script). TS - TC on target 2025-26, source 2024-25:
#       above +2 SE the carried split predicts better than the carried constant; below -2
#       worse; between, the shape cannot be told apart out of sample at that level.
#   R7  THE ESTIMATE IS REPORTED, NEVER JUDGED: per season and arm the boat all-gear effort
#       and catch, December to February's boat catch, the port total.
#   R8  NOTHING HERE CHANGES THE SHIPPED TERM. The method changes on Matt's decision only.
#   EXPECTED, stated so a surprise is visible: R2 above +2 SE; R3 REVIEW on 2025-26 and PASS
#   on 2024-25; R4 (b) on 2025-26 and (a) on 2024-25; R5 at |z| >= 2 with the winter effect
#   the smaller reduction; R6 between -2 and +2.
#
# WHAT IT WRITES, merged by key: 05_output/advisory_replication_2026-11-19_{arms,pairs,
# monthly,verdicts}.csv, and a small evidence folder per arm under
# 05_output/advisory_replication_2026-11-19_arms/ (the run folders themselves are about
# 15 MB each and are not kept; the files the verdicts read are).
#
# HOW TO RUN. Ships DRY_RUN <- TRUE: sourcing it runs the desk stage and fits nothing.
#   Rscript 06_diagnostics/run_advisory_replication_2026-11-19.R --go --stage=S26-const
# renders one arm (a fresh process per arm keeps memory flat); --go with no --stage renders
# every arm in turn; --read reads the rendered arms and writes the tables without fitting.
# Set BSS_FIT_CACHE_DIR first. About 20 to 35 minutes per boat all-gear fit, plus the shore
# fits once per season.
###############################################################################
DRY_RUN <- TRUE                    # TRUE prints the plan and runs the desk stage; fits nothing
.args <- commandArgs(trailingOnly = TRUE)
if ("--go" %in% .args || identical(Sys.getenv("ADVISORY_REPLICATION_GO"), "1")) DRY_RUN <- FALSE
READ_ONLY <- "--read" %in% .args
RESUME <- TRUE                     # reuse an arm ONLY when its ADVISORY_STAGE.txt digest matches
ALL_STAGES <- c("S26-const", "S26-off", "S26-split", "S25-const", "S25-split", "S25-off")
.st <- sub("^--stage=", "", grep("^--stage=", .args, value = TRUE))
STAGES <- if (length(.st)) strsplit(.st[1], ",")[[1]] else ALL_STAGES
if (!all(STAGES %in% ALL_STAGES)) stop("unknown stage: ", paste(setdiff(STAGES, ALL_STAGES), collapse = ", "))

.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
.here <- function(...) file.path(.root, ...)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
banner <- function(msg) cat("\n", strrep("=", 78), "\n ", msg, "\n", strrep("=", 78), "\n", sep = "")
.num1 <- function(x) { v <- suppressWarnings(as.numeric(x)); if (length(v)) v[1] else NA_real_ }
FITTING <- !isTRUE(DRY_RUN) && !READ_ONLY
WRITING <- FITTING || READ_ONLY      # a plain dry run reads and prints, and writes nothing (D38)

if (FITTING) {
  suppressPackageStartupMessages({ library(here); library(rmarkdown) })
  source(.here("03_R_functions", "bss_packages.R")); bss_load_packages()      # B46: the shared loader
  rstan_options(auto_write = TRUE)
} else {
  suppressWarnings(suppressPackageStartupMessages(
    try({ library(dplyr); library(tidyr); library(lubridate); library(readxl); library(here); library(purrr); library(stringr); library(tibble) },
        silent = TRUE)))
}
invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), function(f) source(f)))
source(.here("run_config.R"))
BASE <- run_config
POOLED_RMD <- .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd")
PREFIX <- "pooled-CPUE-"
FIT_BOAT <- "private_boat_all_gear_Dungeness_Kept"
OUT <- function(kind) .here("05_output", sprintf("advisory_replication_2026-11-19_%s.csv", kind))
EVIDENCE <- .here("05_output", "advisory_replication_2026-11-19_arms")

# ---------------------------------------------------------------------------
# THE ARMS: the season's window, then one delta on the marine hazard keys.
# ---------------------------------------------------------------------------
SEASONS <- list(
  S26 = list(label = "2025-26", window = list(
    est_date_start = "2025-09-16", est_date_end = "2026-09-15", season_filter = "2025-26", pot_closures = NULL, census_windows = NULL,
    pot_closure_start = "2025-09-16", pot_closure_end = "2025-11-30", pot_open_date = "2025-12-01",
    census_start_date = "2025-12-01", census_end_date = "2026-01-03", commercial_opener = "2026-01-04")),
  S25 = list(label = "2024-25", window = list()))          # run_config.R as shipped
SPLIT_TERMS <- c("nws_sca_any_winter", "nws_sca_any_rest")
ARMS <- list(
  const = list(item = "the shipped term: nws_sca_any, one coefficient", delta = list()),
  off   = list(item = "no advisory term", delta = list(marine_hazard_mode = "off")),
  split = list(item = "the season split: nws_sca_any_winter + nws_sca_any_rest",
               delta = list(marine_hazard_manual_boat = SPLIT_TERMS,
                            marine_hazard_candidates_boat = unique(c(BASE$marine_hazard_candidates_boat, SPLIT_TERMS)))))
.season_of <- function(sid) sub("-.*$", "", sid)
.arm_of <- function(sid) sub("^[^-]*-", "", sid)
.tag <- function(sid) paste0("AR-", sid)
resolve_cfg <- function(sid) {
  cfg <- BASE
  # 2026-11-27 (B90): the six arms were rendered on 2026-11-19 with every boat day reporting the level;
  # run_config.R ships "carry_without_osp" since A42 (2026-11-26). Pinned here so a re-render of an arm
  # is the arm that was read (S25-const is B70, not B89).
  cfg$boat_day_effect_reporting <- "level"
  w <- SEASONS[[.season_of(sid)]]$window
  for (k in names(w)) cfg[k] <- list(w[[k]])
  d <- ARMS[[.arm_of(sid)]]$delta
  for (k in names(d)) cfg[k] <- list(d[[k]])
  cfg
}
stage_digest <- function(sid) {
  cfg <- resolve_cfg(sid); cfg$run_tag <- NULL; cfg$model <- NULL
  txt <- paste(vapply(sort(names(cfg)), function(k) paste0(k, "=", paste(format(unlist(cfg[[k]] %||% "NULL")), collapse = "|")), character(1)), collapse = ";")
  if (requireNamespace("digest", quietly = TRUE)) substr(digest::digest(txt), 1, 10) else sprintf("%010d", sum(utf8ToInt(txt)) %% 1e10)
}
.stamp_field <- function(dir, key) {
  f <- file.path(dir, "ADVISORY_STAGE.txt"); if (!file.exists(f)) return(NA_character_)
  l <- readLines(f, warn = FALSE); v <- sub(paste0("^", key, ": *"), "", grep(paste0("^", key, ":"), l, value = TRUE))
  if (length(v)) v[1] else NA_character_
}
find_outdir <- function(sid) {
  ds <- Sys.glob(.here("05_output", "*", paste0(PREFIX, .tag(sid), "*")))
  ds <- ds[dir.exists(ds) & file.exists(file.path(ds, "component_totals.csv"))]
  ok <- vapply(ds, function(d) identical(.stamp_field(d, "digest"), stage_digest(sid)), logical(1))
  if (any(ok)) return(ds[ok][length(ds[ok])])
  ev <- file.path(EVIDENCE, sid)                      # the kept evidence folder serves a later --read
  if (file.exists(file.path(ev, "component_totals.csv")) && identical(.stamp_field(ev, "digest"), stage_digest(sid))) return(ev)
  NA_character_
}

# ---------------------------------------------------------------------------
# THE DESK STAGE (always): what each arm is, and whether its data can carry it.
# ---------------------------------------------------------------------------
V <- list()
vrow <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = as.character(observed), threshold = threshold,
                                    verdict = verdict, why = why, stringsAsFactors = FALSE)
desk <- function() {
  banner("DESK STAGE: the arms, the flag and the sampled days behind each coefficient")
  ok_shipped <- identical(BASE$marine_hazard_mode, "manual") && identical(BASE$marine_hazard_manual_boat, "nws_sca_any") &&
    identical(as.character(BASE$marine_hazard_gear_regimes), "all_gear") && !length(BASE$marine_hazard_manual_shore)
  vrow("desk", "run_config.R ships the adopted term, so the const arm is the method of record", ok_shipped, "TRUE",
       if (ok_shipped) "PASS" else "FAIL", "A30: manual, boat nws_sca_any, all_gear only, no shore term")
  for (sn in names(SEASONS)) {
    cfg <- resolve_cfg(paste0(sn, "-const"))
    fl <- tryCatch(marine_hazard_flag_series(cfg, date_start = as.Date(cfg$est_date_start), date_end = as.Date(cfg$est_date_end)), error = function(e) e)
    if (inherits(fl, "error")) { vrow(sn, "the hazard workbook covers the window", conditionMessage(fl), "the flag series builds", "FAIL", "every production run stops on this"); next }
    ag <- fl[fl$event_date >= as.Date(cfg$pot_open_date), , drop = FALSE]
    e <- readxl::read_excel(.here("04_input_files", "effort_combined.xlsx"), sheet = "data")
    e$date <- as.Date(e$date)
    areas <- cfg$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch")
    sd_ <- unique(e$date[e$creel_area %in% areas & !is.na(e$boat_trailer_count) & e$date >= as.Date(cfg$pot_open_date) & e$date <= as.Date(cfg$est_date_end)])
    s <- ag[ag$event_date %in% sd_, , drop = FALSE]
    cat(sprintf("  %s all-gear window: %d days, %d advisory (%d winter, %d rest); sampled launch days %d, advisory %d (%d winter, %d rest)\n",
                SEASONS[[sn]]$label, nrow(ag), sum(ag$nws_sca_any), sum(ag$nws_sca_any_winter), sum(ag$nws_sca_any_rest),
                nrow(s), sum(s$nws_sca_any), sum(s$nws_sca_any_winter), sum(s$nws_sca_any_rest)))
    vrow(sn, "sampled advisory launch days behind each split coefficient (winter / rest)",
         sprintf("%d / %d of %d sampled days", sum(s$nws_sca_any_winter), sum(s$nws_sca_any_rest), nrow(s)), "at least 5 each",
         if (sum(s$nws_sca_any_winter) >= 5 && sum(s$nws_sca_any_rest) >= 5) "PASS" else "REVIEW", "a column with fewer is not identifiable in the window")
    vrow(sn, "advisory share of sampled against unsampled all-gear days",
         sprintf("%.0f%% of %d sampled, %.0f%% of %d unsampled", 100 * mean(s$nws_sca_any), nrow(s),
                 100 * mean(ag$nws_sca_any[!ag$event_date %in% sd_]), sum(!ag$event_date %in% sd_)), "reported", "INFO",
         "A30: this ratio sets the sign of the term's effect on the total")
  }
  for (sid in ALL_STAGES) {
    cfg <- resolve_cfg(sid)
    cat(sprintf("  %-9s %-58s mode %-6s boat {%s}  digest %s  %s\n", sid, ARMS[[.arm_of(sid)]]$item, cfg$marine_hazard_mode,
                paste(cfg$marine_hazard_manual_boat, collapse = ","), stage_digest(sid),
                { d <- find_outdir(sid); if (is.na(d)) "not rendered" else paste("rendered:", basename(dirname(d)), basename(d)) }))
  }
  cat(sprintf("  BSS_FIT_CACHE_DIR: %s\n", { v <- Sys.getenv("BSS_FIT_CACHE_DIR"); if (nzchar(v)) v else "NOT SET (every arm will sample the shore fits again)" }))
}

# ---------------------------------------------------------------------------
# ONE ARM
# ---------------------------------------------------------------------------
run_one <- function(sid) {
  banner(sprintf("%s  %s", sid, ARMS[[.arm_of(sid)]]$item))
  existing <- find_outdir(sid)
  if (isTRUE(RESUME) && !is.na(existing)) { cat("  RESUME: rendered with a matching digest at", existing, "\n"); return(existing) }
  if (!FITTING) { cat("  not fitting (DRY_RUN or --read).\n"); return(NA_character_) }
  cfg <- resolve_cfg(sid); cfg$model <- "pooled"; cfg$run_tag <- .tag(sid)
  run_env <- new.env(parent = globalenv()); run_env$run_config <- cfg
  t0 <- Sys.time()
  html <- rmarkdown::render(POOLED_RMD, envir = run_env, quiet = FALSE)
  od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e) NA_character_)
  if (is.na(od) || !dir.exists(od)) stop(sid, " rendered but its output folder was not found")
  if (file.exists(html) && normalizePath(dirname(html)) != normalizePath(od))
    if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE))) suppressWarnings(file.remove(html))
  writeLines(c(sprintf("stage: %s", sid), sprintf("tag: %s", .tag(sid)), sprintf("digest: %s", stage_digest(sid)),
               sprintf("git: %s", tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) NA)),
               sprintf("tree_clean_of_code: %s", tryCatch(!length(system("git status --porcelain -- 01_BSS_models 02_stan_models 03_R_functions run_config.R", intern = TRUE)), error = function(e) NA)),
               sprintf("platform: %s; rstan %s; StanHeaders %s", R.version.string, as.character(utils::packageVersion("rstan")), as.character(utils::packageVersion("StanHeaders"))),
               sprintf("minutes: %.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))),
               sprintf("rendered: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
             file.path(od, "ADVISORY_STAGE.txt"))
  cat(sprintf("  %s finished in %.1f min -> %s\n", sid, as.numeric(difftime(Sys.time(), t0, units = "mins")), od))
  od
}

# ---------------------------------------------------------------------------
# READING AN ARM
# ---------------------------------------------------------------------------
rd <- function(dir, f) { p <- file.path(dir %||% "", f); if (!file.exists(p)) return(NULL)
  tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL) }
KEEP <- c("ADVISORY_STAGE.txt", "run_parameters.txt", "run_warnings.csv", "run_timings.csv", "convergence_report.csv", "component_totals.csv",
          "pe_vs_bss_comparison.csv", "monthly_estimates.csv", "monthly_estimates_by_mode.csv", "fit_data_summary.csv", "model_adequacy.csv",
          "marine_hazard_selection.csv", "effort_day_covariates.csv",
          paste0(c("opener_covariates_", "structural_params_", "loo_block_trailer_", "loo_block_osp_", "loo_block_joint_", "loo_summary_",
                   "loo_pointwise_trailer_", "loo_pointwise_catch_", "loo_pointwise_osp_", "bss_daily_effort_", "bss_daily_catch_",
                   "sampler_diagnostics_", "advisory_coefficients_"), FIT_BOAT, ".csv"))
coef_rows <- function(dir) {
  kept <- rd(dir, paste0("advisory_coefficients_", FIT_BOAT, ".csv"))        # an evidence folder carries the extract, not the 630 KB summary
  if (!is.null(kept)) return(kept)
  lab <- rd(dir, sprintf("opener_covariates_%s.csv", FIT_BOAT))
  f <- file.path(dir, sprintf("bss_full_summary_%s.csv", FIT_BOAT))
  if (is.null(lab) || !file.exists(f)) return(NULL)
  fs <- utils::read.csv(f, row.names = 1, check.names = FALSE)
  out <- lapply(seq_len(nrow(lab)), function(i) { r <- fs[as.character(lab$parameter[i]), , drop = FALSE]
    data.frame(term = as.character(lab$opener[i]), mean = .num1(r[["mean"]]), sd = .num1(r[["sd"]]), q2.5 = .num1(r[["2.5%"]]), q50 = .num1(r[["50%"]]),
               q97.5 = .num1(r[["97.5%"]]), n_eff = .num1(r[["n_eff"]]), rhat = .num1(r[["Rhat"]]), stringsAsFactors = FALSE) })
  do.call(rbind, out)
}
keep_evidence <- function(sid, dir) {
  ev <- file.path(EVIDENCE, sid); if (normalizePath(dir) == normalizePath(ev, mustWork = FALSE)) return(invisible(ev))
  dir.create(ev, recursive = TRUE, showWarnings = FALSE)
  co <- coef_rows(dir)
  if (!is.null(co)) utils::write.csv(co, file.path(ev, paste0("advisory_coefficients_", FIT_BOAT, ".csv")), row.names = FALSE)
  for (f in KEEP) if (file.exists(file.path(dir, f))) file.copy(file.path(dir, f), file.path(ev, f), overwrite = TRUE)
  invisible(ev)
}
read_arm <- function(sid, dir) {
  cr <- rd(dir, "convergence_report.csv"); ct <- rd(dir, "component_totals.csv"); pv <- rd(dir, "pe_vs_bss_comparison.csv")
  mm <- rd(dir, "monthly_estimates_by_mode.csv"); co <- coef_rows(dir)
  # R1 is about FITS. A component the driver did not fit (too few interviews: its row carries no AR resolution and
  # "PE (insufficient data)") is not a failed fit; the first reading of S26-const counted the 2025-26 boat pot-closure
  # component, which has two interviews, as one (CHANGE_REGISTER C row, 2026-11-19).
  fitted <- if (is.null(cr)) logical(0) else !is.na(cr$ar_resolution)
  gate <- !is.null(cr) && any(fitted) && all(as.logical(cr$pass_convergence[fitted])) && isTRUE(fitted[cr$fit == FIT_BOAT])
  ar <- if (is.null(cr)) NA else cr$ar_resolution[cr$fit == FIT_BOAT][1]
  comp <- function(nm, col) if (is.null(ct)) NA_real_ else .num1(ct[[col]][ct$Component == nm])
  boat_m <- if (is.null(mm)) NULL else mm[mm$Mode == "Private Boat", , drop = FALSE]
  wm <- as.integer(BASE$marine_hazard_winter_months %||% c(12, 1, 2))
  winter <- if (is.null(boat_m)) NA_real_ else sum(boat_m$median[as.integer(substr(boat_m$Month, 6, 7)) %in% wm])
  term <- function(nm, col) if (is.null(co) || !nm %in% co$term) NA_real_ else co[[col]][co$term == nm][1]
  ck <- rd(file.path(EVIDENCE, sid), "blockint_checks.csv")          # desk_advisory_block_cv_2026-11-19.R (R3a); absent until it has run
  ckv <- function(col) if (is.null(ck) || !col %in% names(ck)) NA_real_ else .num1(ck[[col]])
  arm <- data.frame(
    stage = sid, season = SEASONS[[.season_of(sid)]]$label, arm = .arm_of(sid), folder = sub(paste0("^", .root, "/"), "", dir),
    platform = .stamp_field(dir, "platform"), git = .stamp_field(dir, "git"), minutes = .num1(.stamp_field(dir, "minutes")),
    gate_all_fits = gate, n_fits = if (is.null(cr)) NA_integer_ else sum(fitted), boat_all_gear_ar = ar,
    boat_all_gear_divergence_fraction = if (is.null(cr)) NA_real_ else .num1(cr$divergence_fraction[cr$fit == FIT_BOAT]),
    boat_all_gear_catch = comp("Private boat, All gear", "Median"), boat_all_gear_catch_lo = comp("Private boat, All gear", "lo"),
    boat_all_gear_catch_hi = comp("Private boat, All gear", "hi"),
    boat_all_gear_effort = if (is.null(pv)) NA_real_ else .num1(pv$BSS_effort[pv$component == "private_boat (All gear)"]),
    boat_winter_months_catch = winter, port_total = comp("PORT TOTAL", "Median"), port_lo = comp("PORT TOTAL", "lo"), port_hi = comp("PORT TOTAL", "hi"),
    b_const = term("nws_sca_any", "mean"), b_const_sd = term("nws_sca_any", "sd"), b_const_lo = term("nws_sca_any", "q2.5"), b_const_hi = term("nws_sca_any", "q97.5"),
    b_winter = term("nws_sca_any_winter", "mean"), b_winter_sd = term("nws_sca_any_winter", "sd"), b_winter_lo = term("nws_sca_any_winter", "q2.5"), b_winter_hi = term("nws_sca_any_winter", "q97.5"),
    b_rest = term("nws_sca_any_rest", "mean"), b_rest_sd = term("nws_sca_any_rest", "sd"), b_rest_lo = term("nws_sca_any_rest", "q2.5"), b_rest_hi = term("nws_sca_any_rest", "q97.5"),
    b_rhat_max = if (is.null(co)) NA_real_ else max(co$rhat), b_neff_min = if (is.null(co)) NA_real_ else min(co$n_eff),
    day_effect_sd = ckv("sigma_uE_mean"), day_effect_sd_sd = ckv("sigma_uE_sd"), day_effect_sd_lo = ckv("sigma_uE_q2.5"), day_effect_sd_hi = ckv("sigma_uE_q97.5"),
    block_weeks = ckv("weeks"), block_weeks_reliable_conditional = ckv("reliable_conditional"), block_weeks_reliable_integrated = ckv("reliable_integrated"),
    stringsAsFactors = FALSE)
  monthly <- if (is.null(boat_m)) NULL else data.frame(stage = sid, season = arm$season, arm = arm$arm, month = boat_m$Month, boat_catch_median = boat_m$median,
                                                       boat_catch_lo95 = boat_m$lo95, boat_catch_hi95 = boat_m$hi95, stringsAsFactors = FALSE)
  list(arm = arm, monthly = monthly)
}
block_tab <- function(dir, stream) rd(dir, sprintf("loo_block_%s_%s.csv", stream, FIT_BOAT))
pair_rows <- function(sn, a, b, dirs) {          # arm a minus arm b
  da <- dirs[[paste0(sn, "-", a)]]; db <- dirs[[paste0(sn, "-", b)]]
  if (is.null(da) || is.null(db) || is.na(da) || is.na(db)) return(NULL)
  out <- list()
  for (stream in c("joint", "trailer", "osp")) {
    ta <- block_tab(da, stream); tb <- block_tab(db, stream); if (is.null(ta) || is.null(tb)) next
    cmp <- bss_block_cv_compare(tb, ta)          # the helper returns its SECOND table minus its first; this row is a minus b
    out[[stream]] <- cbind(data.frame(season = SEASONS[[sn]]$label, pair = sprintf("%s - %s", a, b), statistic = "block CV, held-out weeks", stream = stream, stringsAsFactors = FALSE),
                           cmp[, c("n_blocks", "n_used", "reliable_share", "diff", "se", "ratio", "n_positive", "identical")])
  }
  for (stream in c("catch", "trailer")) {
    pa <- file.path(da, sprintf("loo_pointwise_%s_%s.csv", stream, FIT_BOAT)); pb <- file.path(db, sprintf("loo_pointwise_%s_%s.csv", stream, FIT_BOAT))
    if (!file.exists(pa) || !file.exists(pb)) next
    el <- tryCatch(loo_elpd_paired(pb, pa), error = function(e) NULL); if (is.null(el)) next     # likewise: second minus first
    out[[paste0("pt_", stream)]] <- data.frame(season = SEASONS[[sn]]$label, pair = sprintf("%s - %s", a, b), statistic = "pointwise LOO, paired", stream = stream,
      n_blocks = el$n %||% NA, n_used = el$n %||% NA, reliable_share = NA_real_, diff = el$elpd_diff, se = el$se_diff,
      ratio = if (isTRUE(el$se_diff > 0)) el$elpd_diff / el$se_diff else NA_real_, n_positive = NA_integer_,
      identical = isTRUE(abs(el$elpd_diff) < 1e-8), stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
band <- function(ratio, hi, mid, lo) if (!is.finite(ratio)) "REVIEW" else if (ratio > 2) hi else if (ratio < -2) lo else mid

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
banner(sprintf("ADVISORY TERM ON A SECOND SEASON (B80)  DRY_RUN = %s  READ_ONLY = %s  stages: %s", DRY_RUN, READ_ONLY, paste(STAGES, collapse = " ")))
desk()
DIRS <- list()
for (sid in STAGES) DIRS[[sid]] <- tryCatch(run_one(sid), error = function(e) { cat("  *** ", sid, " FAILED: ", conditionMessage(e), "\n", sep = ""); NA_character_ })
for (sid in setdiff(ALL_STAGES, STAGES)) DIRS[[sid]] <- find_outdir(sid)       # the arms rendered by earlier calls
have <- names(DIRS)[!vapply(DIRS, function(d) is.null(d) || is.na(d), logical(1))]
if (!length(have)) { cat("\nNo arm is rendered yet; the desk stage is all there is to report.\n")
} else {
  banner("READING THE RENDERED ARMS")
  R <- lapply(have, function(sid) read_arm(sid, DIRS[[sid]])); names(R) <- have
  arms <- do.call(rbind, lapply(R, `[[`, "arm")); monthly <- do.call(rbind, lapply(R, `[[`, "monthly"))
  if (WRITING) for (sid in have) keep_evidence(sid, DIRS[[sid]])
  # R1
  for (sn in names(SEASONS)) { ref <- arms$boat_all_gear_ar[arms$stage == paste0(sn, "-const")]
    for (sid in grep(paste0("^", sn, "-"), have, value = TRUE)) { a <- arms[arms$stage == sid, ]
      okk <- isTRUE(a$gate_all_fits) && (!length(ref) || identical(a$boat_all_gear_ar, ref))
      arms$eligible[arms$stage == sid] <- okk
      vrow(sid, "R1 eligibility: every fit passes the gate, the boat all-gear fit at the const arm's AR resolution",
           sprintf("gate %s; AR %s (const %s); boat all-gear divergent %.2f%%", a$gate_all_fits, a$boat_all_gear_ar, if (length(ref)) ref else "not rendered",
                   100 * a$boat_all_gear_divergence_fraction), "gate TRUE, same resolution", if (okk) "PASS" else "FAIL", "an ineligible arm is not read") } }
  elig <- function(sid) isTRUE(arms$eligible[arms$stage == sid])
  pairs <- do.call(rbind, lapply(names(SEASONS), function(sn) rbind(pair_rows(sn, "const", "off", DIRS), pair_rows(sn, "split", "const", DIRS), pair_rows(sn, "split", "off", DIRS))))
  # R3
  for (sn in names(SEASONS)) { if (is.null(pairs)) break
    p <- pairs[pairs$season == SEASONS[[sn]]$label & pairs$pair == "const - off" & pairs$statistic == "block CV, held-out weeks", , drop = FALSE]
    if (!nrow(p)) next
    pr <- if ("joint" %in% p$stream) p[p$stream == "joint", ] else p[p$stream == "trailer", ]
    ev <- isTRUE(pr$reliable_share >= 0.7) && isTRUE(pr$n_used >= 8) && elig(paste0(sn, "-const")) && elig(paste0(sn, "-off"))
    vrow(sn, sprintf("R3 const against off on held-out weeks (%s row)", pr$stream),
         sprintf("%+.1f nats, paired SE %.1f, ratio %+.2f; %d of %d weeks reliable, %d positive", pr$diff, pr$se, pr$ratio, pr$n_used, pr$n_blocks, pr$n_positive),
         "above +2 SE PASS; below -2 FAIL; evaluable at 70% and 8 weeks", if (!ev) "REVIEW" else band(pr$ratio, "PASS", "REVIEW", "FAIL"),
         if (!ev) "not evaluable (unreliable weeks or an ineligible arm)" else "the term on the season's held-out weeks") }
  # R3a (the addendum): the same comparison on the tables with the day deviation integrated out. The desk script
  # desk_advisory_block_cv_2026-11-19.R computes them and their preconditions; this runner only reads them.
  bi <- rd(dirname(OUT("blockint")), basename(OUT("blockint")))
  bi_row <- function(sn, pair) { if (is.null(bi)) return(NULL)
    x <- bi[bi$season == SEASONS[[sn]]$label & bi$pair == pair & grepl("integrated", bi$table) & bi$week_set == "all weeks", , drop = FALSE]
    if (!nrow(x)) NULL else if ("joint" %in% x$stream) x[x$stream == "joint", ] else x[x$stream == "trailer", ] }
  pre_ok <- function(sid) { k <- rd(file.path(EVIDENCE, sid), "blockint_checks.csv")
    !is.null(k) && isTRUE(as.logical(k$precondition_i)) && isTRUE(as.logical(k$precondition_ii)) }
  for (sn in names(SEASONS)) { pr <- bi_row(sn, "const - off"); if (is.null(pr)) next
    pre <- pre_ok(paste0(sn, "-const")) && pre_ok(paste0(sn, "-off"))
    ev <- pre && isTRUE(as.logical(pr$evaluable)) && elig(paste0(sn, "-const")) && elig(paste0(sn, "-off"))
    vrow(sn, sprintf("R3a const against off on held-out weeks, the day deviation integrated out (%s row)", pr$stream),
         sprintf("%+.1f nats, paired SE %.1f, ratio %+.2f; %d of %d weeks reliable, %d positive; preconditions (i) and (ii) %s", pr$diff, pr$se, pr$ratio,
                 pr$n_used, pr$n_weeks, pr$n_positive, if (pre) "met in both arms" else "NOT met"),
         "above +2 SE PASS; below -2 FAIL; evaluable at 70% and 8 weeks", if (!ev) "REVIEW" else band(pr$ratio, "PASS", "REVIEW", "FAIL"),
         if (!ev) "not evaluable (a precondition, unreliable weeks or an ineligible arm)" else "the addendum's instrument, reported beside R3") }
  # R4 and R5
  dd <- list()
  for (sn in names(SEASONS)) { sid <- paste0(sn, "-split"); if (!sid %in% have) next
    a <- arms[arms$stage == sid, ]; if (!is.finite(a$b_winter) || !is.finite(a$b_rest)) { vrow(sid, "R4 the split coefficients", "a split column is absent from the fit", "both present", "REVIEW", "the driver dropped an unidentifiable column"); next }
    idw <- a$b_winter_lo > 0 || a$b_winter_hi < 0; idr <- a$b_rest_lo > 0 || a$b_rest_hi < 0
    d <- a$b_winter - a$b_rest; s <- sqrt(a$b_winter_sd^2 + a$b_rest_sd^2); differ <- abs(d) > 1.96 * s
    branch <- if (!idw) "(a) winter not identified" else if (idr && differ) "(b) both identified, and they differ" else if (idr) "(c) both identified, no difference detectable" else "rest not identified"
    extra <- ""
    if (identical(branch, "(b) both identified, and they differ") && !is.null(pairs)) {
      pc <- pairs[pairs$season == a$season & pairs$pair == "split - const", , drop = FALSE]
      cat_ok <- { x <- pc[pc$statistic == "pointwise LOO, paired" & pc$stream == "catch", ]; nrow(x) == 1 && is.finite(x$ratio) && abs(x$ratio) <= 2 }
      blk <- { x <- pc[pc$statistic == "block CV, held-out weeks", ]; x <- if ("joint" %in% x$stream) x[x$stream == "joint", ] else x[x$stream == "trailer", ]; nrow(x) == 1 && is.finite(x$ratio) && x$ratio >= -2 }
      b4 <- bi_row(sn, "split - const"); blk_i <- !is.null(b4) && isTRUE(as.logical(b4$evaluable)) && is.finite(b4$ratio) && b4$ratio >= -2
      extra <- sprintf("; catch stream unmoved %s; block CV not below -2 SE: %s on the conditional table (R4 as written), %s on the integrated one (the addendum%s): the split %s its place on this season",
                       cat_ok, blk, blk_i, if (is.null(b4)) ", not computed yet" else sprintf(", %+.1f nats at %+.2f SE over %d of %d weeks", b4$diff, b4$ratio, b4$n_used, b4$n_weeks),
                       if (cat_ok && blk && blk_i) "EARNS" else if (cat_ok && blk_i) "EARNS on the integrated table only" else "does not earn")
    }
    vrow(sid, "R4 the season split on this season's boat all-gear fit",
         sprintf("winter %.2f [%.2f, %.2f], rest %.2f [%.2f, %.2f]; winter - rest %+.2f (SD %.2f, z %+.2f)%s", a$b_winter, a$b_winter_lo, a$b_winter_hi,
                 a$b_rest, a$b_rest_lo, a$b_rest_hi, d, s, d / s, extra), "branches (a), (b), (c)", if (elig(sid)) "INFO" else "REVIEW", branch)
    if (elig(sid)) dd[[sn]] <- c(d = d, s = s)
  }
  if (length(dd)) { d <- vapply(dd, `[[`, numeric(1), "d"); s <- vapply(dd, `[[`, numeric(1), "s"); w <- 1 / s^2
    m <- sum(w * d) / sum(w); se <- sqrt(1 / sum(w)); z <- m / se
    het <- if (length(d) == 2) (d[1] - d[2]) / sqrt(sum(s^2)) else NA_real_
    read <- if (length(d) < 2) "one season only: not the two-season statistic" else if (isTRUE(abs(het) >= 2)) "the two seasons disagree with each other: not read as established" else
      if (abs(z) >= 2) sprintf("ESTABLISHED over the two seasons: the winter effect is the %s reduction", if (m > 0) "smaller" else "larger") else if (abs(z) >= 1) "suggestive, not established" else "no difference detectable in two seasons"
    vrow("both seasons", "R5 winter - rest, precision-weighted over the eligible split arms",
         sprintf("%+.2f (SE %.2f, z %+.2f) from %s; the seasons differ from each other at z %+.2f", m, se, z, paste(sprintf("%s %+.2f (SD %.2f)", names(d), d, s), collapse = ", "), het),
         "|z| >= 2 established; 1 to 2 suggestive; below 1 none", "INFO", read) }
  # R7
  for (sn in names(SEASONS)) { g <- function(arm, col) { v <- arms[[col]][arms$stage == paste0(sn, "-", arm)]; if (length(v)) v else NA_real_ }
    if (is.finite(g("const", "boat_all_gear_catch")) && is.finite(g("off", "boat_all_gear_catch")))
      vrow(sn, "R7 the term's price: const against off", sprintf("boat all-gear catch %.0f against %.0f (%+.1f%%); port %.0f against %.0f (%+.1f%%)",
           g("const", "boat_all_gear_catch"), g("off", "boat_all_gear_catch"), 100 * (g("const", "boat_all_gear_catch") / g("off", "boat_all_gear_catch") - 1),
           g("const", "port_total"), g("off", "port_total"), 100 * (g("const", "port_total") / g("off", "port_total") - 1)), "reported, never judged", "INFO", "rule 7 of the ladder")
    if (is.finite(g("const", "boat_all_gear_catch")) && is.finite(g("split", "boat_all_gear_catch")))
      vrow(sn, "R7 the shape: split against const", sprintf("boat all-gear catch %.0f against %.0f (%+.1f%%); December to February boat catch %.0f against %.0f (%+.1f%%); port %.0f against %.0f (%+.1f%%)",
           g("split", "boat_all_gear_catch"), g("const", "boat_all_gear_catch"), 100 * (g("split", "boat_all_gear_catch") / g("const", "boat_all_gear_catch") - 1),
           g("split", "boat_winter_months_catch"), g("const", "boat_winter_months_catch"), 100 * (g("split", "boat_winter_months_catch") / g("const", "boat_winter_months_catch") - 1),
           g("split", "port_total"), g("const", "port_total"), 100 * (g("split", "port_total") / g("const", "port_total") - 1)), "reported, never judged", "INFO", "rule 7 of the ladder") }
  op <- options(width = 250)
  cat("\n== THE ARMS\n"); print(arms[, c("stage", "eligible", "boat_all_gear_ar", "boat_all_gear_catch", "boat_all_gear_effort", "boat_winter_months_catch", "port_total",
                                         "b_const", "b_const_sd", "b_winter", "b_winter_sd", "b_rest", "b_rest_sd", "minutes")], row.names = FALSE, digits = 4)
  if (!is.null(pairs)) { cat("\n== THE PAIRS\n"); print(pairs, row.names = FALSE, digits = 3) }
  options(op)
  if (WRITING) {
    merge_csv_by(arms, OUT("arms"), key = "stage")
    if (!is.null(monthly)) merge_csv_by(monthly, OUT("monthly"), key = c("stage", "month"))
    if (!is.null(pairs)) merge_csv_by(pairs, OUT("pairs"), key = c("season", "pair", "statistic", "stream"))
  }
}
if (length(V)) { vd <- do.call(rbind, V); cat("\n== VERDICTS\n"); op <- options(width = 250); print(vd[, c("stage", "criterion", "observed", "verdict", "why")], row.names = FALSE, right = FALSE); options(op)
  if (WRITING) merge_csv_by(vd, OUT("verdicts"), key = c("stage", "criterion")) else cat("  (a dry run writes nothing: D38)\n") }
cat("\nREAD BY THE RULE IN THIS FILE'S HEADER AND CHANGE_REGISTER B80. R2 and R6 are the transfer script's (desk_advisory_transfer_2026-11-19.R).\n")
