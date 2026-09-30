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
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program (see the LICENSE file);
# if not, see https://www.gnu.org/licenses/ .
# -----------------------------------------------------------------------------
###############################################################################
# THE AUTHORITATIVE RENDER AT THE A31 CODE, AND THE RUNS THAT CLOSE D29, D3 AND D6
# (2026-09-29, CHANGE_REGISTER B58)
# -----------------------------------------------------------------------------
# WHY THIS RUN EXISTS
#   main carries A31 (both single-section level hierarchies collapsed) and init_r = 0.5 (B51),
#   and neither has had a production render together: the only render at init_r = 0.5 (B57)
#   was of the pre-A31 code and failed its shore pot-closure gate. Three register items
#   also wait on renders, each named by the register as the way to close it:
#     D29  the boat all-gear AR period (monthly) has never been tested on either track, and
#          it is the largest sensitivity after moored boats (D14): +10% to +25% on a component
#          that is about half the port, measured on older code by the confounded D3 ladder.
#     D3   the gear track could not express a per-population AR period; B32 added the form.
#          One gear render at the pooled track's per-population periods closes it.
#     D6   the gear-track zero-inflated shore catch "earns its parameter"; adoption waits on
#          one render at the matched configuration.
#   This file renders all of them overnight, reads every folder back, applies the decision
#   rules written BELOW BEFORE THE RUN, and writes the verdicts.
#
# HOW TO RUN IT. In the RStudio CONSOLE, from the project (the route that works on Matt's
# machine, where Rscript is not on the Terminal's PATH):
#     source("06_diagnostics/run_authoritative_batch_2026-09-29.R")                              # dry run
#     Sys.setenv(BSS_BATCH_GO = "1"); source("06_diagnostics/run_authoritative_batch_2026-09-29.R")  # fits
#   Or, where Rscript and pandoc are on the PATH:
#     Rscript 06_diagnostics/run_authoritative_batch_2026-09-29.R [--go]
#   Every stage runs in an Rscript started by this session and loads packages from THIS
#   session's library (B59; see CHILD_ENV below), so the session you start it from decides the
#   package versions. Desk check S0 confirms a stage's process sees the same ones.
#   DO NOT edit DRY_RUN in this file to start it. The authoritative stage's manifest records
#   every TRACKED file that differs from the commit, so an edited runner would mark the
#   authoritative render as made from a modified tree. --go (or the environment variable)
#   starts it with the tree clean.
#   Re-running with --go after an interruption RESUMES: a stage whose folder carries a
#   matching AB_STAGE.txt digest is read back, not refitted. The digest covers the stage's
#   whole resolved configuration and the code fingerprint, so a changed run_config.R value or
#   a changed driver, Stan or function file refits the stage (conservative on purpose: an
#   authoritative number must come from the code and configuration it is filed under).
#
# THE STAGES (in run order; STAGES below can drop any; each is a separate R process, so a
# stage that errors is recorded and the batch goes on to the next)
#   S0   desk, a minute, no MCMC. The shipped configuration is the method of record; every
#        stage differs from it in its declared keys only; on the real 2024-25 inputs the
#        gear per-population periods resolve to the pooled track's caps fit by fit; the two
#        ladders resolve to weekly/biweekly/monthly for the boat all-gear fit and to one rung
#        for every other fit; the D6 stage builds zi_catch = 1 for the shore fits only; the
#        NWS workbook covers the window; the reference folders exist.
#   A    THE AUTHORITATIVE RENDER: run_estimation.R --model both, run_config.R exactly as
#        shipped, through the production orchestrator (manifest, cross-check). Both tracks.
#   D3   gear track with a PER-POPULATION AR period: the shore at the pooled track's caps
#        (all-gear weekly, pot closure biweekly), the boat at its shipped gear periods
#        (all-gear monthly, pot closure biweekly). Against the shipped gear track exactly ONE
#        fit's period moves, shore all-gear monthly -> weekly, which is the configuration the
#        register measured (+0.28% at the port against the pooled track, CHANGE_REGISTER D3)
#        and the adoption form documented in run_config.R section 5. (The gear boat
#        pot-closure period, biweekly against the pooled monthly, is left where it ships: no
#        register item asks about it, and moving it here would put a second, unjudged change
#        into D3's ADOPT.)
#   D6   D3 plus the zero-inflated shore catch on the gear track (catch_zi_tracks = both).
#        Judged against D3, which differs from it in that one key.
#   D29P POOLED: the escalation ladder on the boat all-gear fit ONLY (ar_escalate =
#        list(private_boat = "all_gear"), ar_escalate_stop = "all_rungs"), rungs weekly,
#        biweekly, monthly; the cap is ignored (ar_escalate_respect_cap = FALSE) because the
#        cap is the thing under test. Every rung's gate, estimate, adequacy (ar_rung_adequacy)
#        and, since this patch, its leave-one-week-out block CV (ladder_block_*.csv) are
#        written, not only the kept rung's.
#   D29G GEAR: the same ladder on the gear track, on D3's configuration (so the shore fits
#        are matched and only the boat all-gear period moves). Corroboration, never the
#        decider (rule D29-5).
#   R2   POOLED at init_r = 2, rstan's default: does the A31 model still need the small
#        radius? Judged against A, which differs from it in that one key.
#   D29D POOLED boat all-gear at DAILY (ar_force), the finest rung, last because it is the
#        slowest and the least likely to be adopted (the gear track failed adequacy there on
#        2026-09-14: p_loo 0.2751, 11 bad k). Its block CV is compared with the ladder's monthly rung.
#
# RUNTIME, from measured renders (your machine, 2026-09-28/29) and the A31 container refits:
#   A     pooled ~1.0-1.5 h (A31 shore fits were 5-7 min in the container, where the old
#         model took 30; B57 took 100 min for shore all-gear and 40 for boat all-gear on your
#         machine) + gear ~0.6 h (35.9 min on 2026-09-29)                      ~1.5-2 h
#   D3    gear, shore all-gear at weekly (P_n 43 against 10)                    ~0.7 h
#   D6    as D3                                                                 ~0.7 h
#   D29P  pooled, the boat all-gear fit three times (weekly and biweekly are
#         several times its monthly P_n of 10)                                  ~2.5-3.5 h
#   D29G  gear, the boat all-gear fit three times                               ~1.5 h
#   R2    as A's pooled half                                                    ~1-1.5 h
#   D29D  pooled with the boat all-gear fit at daily (P_n 289)                  ~2-4 h, unknown
#   TOTAL about 10-14 h. A to D29G (the answers the register asks for) are about 7-9 h.
#   Stop it any time; --go again resumes. Each stage also has a wall-clock limit (TIMEOUT_H
#   below; generous, several times the expected time) after which it is stopped and recorded
#   FAILED, so one hung fit cannot hold the rest of the night. If stage A is interrupted after
#   one model's folder was stamped, --go renders the other model only and the runner writes
#   the cross-check itself (a kill before the orchestrator's manifest was written leaves
#   nothing stamped, and A starts again).
#
# PREREQUISITES S0 CHECKS FIRST: rmarkdown must find pandoc (RStudio sets RSTUDIO_PANDOC in its
# Console, and the stages inherit it; a plain cmd / PowerShell window may not have pandoc on
# PATH), and a stage's R process must load this session's packages at the same versions (B59:
# the first overnight attempt, 2026-09-29, failed every stage because the stages activated the
# project's renv library, which is incomplete on Matt's machine, while the Console used the
# library the authoritative run was rendered with).
#
# ============================================================================
# THE DECISION RULES, STATED HERE BEFORE THE RUN so they cannot be fitted to the answer.
# ============================================================================
# A (can this render become the authoritative run?). Every clause is read from the files.
#   A1  every fit on both tracks passes the convergence gate (all 8 report BSS).
#   A2  the pooled shore all-gear fit is under 1% divergent (the A31 container refits were
#       0.05% to 0.19% at two seeds and two radii; the old model 1.6% to 23.3%).
#   A3  no pooled chain is stuck: no chain's mean step size is below one tenth of its fit's
#       median chain step size, and no fit has 5% or more of its iterations at maximum tree
#       depth. (B57's stuck chain: step 2.6e-5 against 0.006 to 0.014; 25% saturation.)
#   A4  each fitted component against the authoritative run of 2026-09-28, in units of this
#       render's posterior SD: within 0.2 PASS, beyond it REVIEW (never FAIL: A31 changes
#       the model as well as the draws, and a real move is recorded, not rejected).
#   A5  the gear cross-check is within run_config's 2% (the orchestrator's own verdict).
#   RESULT: A1, A2, A3 and A5 all PASS -> ELIGIBLE to become the authoritative run (the box in
#   PIPELINE_STATUS.md), with A4 read for what moved. Anything else -> NOT ELIGIBLE, and the
#   2026-09-28 run stays authoritative.
#
# R2 (does A31 still need init_r = 0.5?)
#   R2-1 every pooled fit at init_r = 2 passes the gate; R2-2 shore all-gear under 1%
#   divergent; R2-3 no stuck chain (A3's test); R2-4 each component within 0.2 posterior SD
#   of A. All four -> "the A31 model does not depend on the radius; 0.5 is a precaution, not
#   a crutch" (keep or revert it; either is defensible, and the method is the model).
#   Any failure -> "keep 0.5", with the failing clause.
#
# D3 (the gear track's AR period, per population). As the 2026-09-13 runner stated them:
#   D3-1 every fit in D3 passes the gate.
#   D3-2 the gear shore all-gear fit at weekly: worst-stream p_loo_frac <= 0.15, bad Pareto k
#        <= 5% of its catch-stream n_obs (the 2026-09-14 gear ladder: 0.0938, 0 bad k).
#   D3-3 every fit ran at the period asked for (shore all-gear weekly, shore pot closure
#        biweekly, boat all-gear monthly, boat pot closure biweekly). A fit that silently used
#        another period would make the comparison meaningless.
#   D3-4 REPORTED, NOT A CRITERION: the cross-track gap against A's pooled track, and the
#        shore all-gear component's move (monthly to weekly). Choosing the gear track's
#        period to shrink the gap would tune one estimate to another (the 2026-09-13 rule 5).
#   RESULT: D3-1 to D3-3 PASS -> ADOPT gear_period_bss in the per-population form above; D3
#   closes, and the shore half of the cross-check becomes like-for-like. Otherwise the
#   failing clause stands.
#
# D6 (the gear-track zero-inflated shore catch), as the 2026-09-13 runner stated them,
# D6 against D3:
#   D6-1 every fit in D6 passes the gate.
#   D6-2 the PAIRED elpd gain on the shore all-gear catch stream exceeds 2 paired SE
#        (loo_elpd_paired(); Vehtari et al. 2017, s3.3).
#   D6-3 both count bins improve: |z| of the zero bin and of the one bin both fall.
#   D6-4 the gain is not bought entirely at the zeros: the positive counts do not lose more
#        than the zeros gain.
#   D6-5 REPORTED: theta_C against its Beta(1, 9) prior mean; the pot-closure fit's pair;
#        the port total (never a criterion: the total is scaled by 1 - theta_C by design).
#   RESULT, read from the four clauses' verdicts: any FAIL -> DO NOT ADOPT; otherwise any
#   REVIEW (a clause that could not be computed, or only one count bin improving) ->
#   UNDECIDED, which blocks adoption without counting against it; otherwise (D6-1 to D6-4 all
#   PASS) -> ADOPT catch_zi_tracks = c("pooled", "gear_resolved"); D6 closes.
#
# D29 (the boat all-gear AR period). Rungs: weekly, biweekly, monthly on each track (D29P,
# D29G), and daily on the pooled track (D29D).
#   D29-0 COMPLETE: the pooled ladder log holds the boat all-gear fit at exactly weekly,
#         biweekly and monthly with the ladder engaged; otherwise the rule does not run (c).
#   D29-1 ELIGIBLE: the rung passes the convergence gate (an NA gate is not a pass).
#   D29-2 ADEQUATE: worst-stream p_loo_frac <= 0.15 and bad Pareto k <= 5% of the fit's
#         LOO observations over every stream (trailer, OSP and catch; n_pareto_bad counts
#         across all of them, so the denominator does too). An adequacy that could not be
#         computed is UNKNOWN, not inadequate, and sends the decision to (c).
#   D29-3 INTERPOLATION, the test an AR period is actually about: the leave-one-week-out
#         block CV of the boat's effort streams (trailer and OSP held out together, the
#         "joint" table; bss_block_cv.R), each finer rung paired week by week with the
#         LADDER'S MONTHLY rung (the daily stage too; A's monthly fit only if the ladder has
#         none) over the weeks reliable (Pareto k <= 0.7) in both. Gain > +2 paired SE:
#         BETTER. Below -2: WORSE. Between: NO EVIDENCE. Evaluable only if at least half of
#         the weeks both tables hold are reliable in both; otherwise NOT EVALUABLE.
#   D29-4 THE DECISION, pooled track only:
#         (a) a finer rung that is eligible, adequate and BETTER, and no rung finer than it
#             left undecided -> recommend the finest such rung for
#             ar_max_resolution$pooled$private_boat (all-gear), and the gear track to match;
#         (b) monthly eligible and adequate, and every finer eligible adequate rung WORSE or
#             NO EVIDENCE -> the data do not support moving the cap. NOT "monthly is right":
#             the span across the adequate rungs is the resolution uncertainty of this
#             component and is reported as such. The cap stays monthly unless Matt decides to
#             apply the shore precedent (finest adequate rung, 2026-09-07), which this rule
#             deliberately does not apply automatically, because the boat rungs differ by
#             +10% to +25% where the shore rungs differed by 3.6%;
#         (c) anything else (an incomplete ladder, an unknown adequacy or gate, a finer rung
#             NOT EVALUABLE on the block CV) -> REVIEW, with the reason. A rung that fails
#             the gate or is inadequate is DECIDED (out of the running), not undecided.
#   D29-5 NOT CRITERIA: agreement between the tracks (circular; the D3 rule's clause 5) and
#         the port total. The gear ladder corroborates or contradicts, and is reported.
#   D29-6 REPORTED: the ladder's own monthly rung against A's monthly boat fit, in posterior
#         SD. Near zero is the control that the ladder machinery changes nothing but the rung.
#
# WHAT IT WRITES (recomputed in full from the folders on every invocation):
#   05_output/authoritative_batch_2026-09-29_stages.csv          stage -> folder(s), status, digest
#   05_output/authoritative_batch_2026-09-29_components.csv      every stage's components
#   05_output/authoritative_batch_2026-09-29_d29_rungs.csv       the boat all-gear rungs, both tracks
#   05_output/authoritative_batch_2026-09-29_verdicts.csv        every clause, PASS/FAIL/REVIEW/INFO
#   05_output/authoritative_batch_2026-09-29_recommendation.csv  A, R2, D3, D6, D29
#   05_output/authoritative_batch_2026-09-29_logs/<stage>.log    each stage's R console
###############################################################################

DRY_RUN <- TRUE                    # ships TRUE; start with --go or BSS_BATCH_GO=1 (see above)
STAGES  <- c("S0", "A", "D3", "D6", "D29P", "D29G", "R2", "D29D")
RESUME  <- TRUE                    # reuse a stage ONLY when its AB_STAGE.txt digest matches
if ("--go" %in% commandArgs(trailingOnly = TRUE) || identical(Sys.getenv("BSS_BATCH_GO"), "1")) DRY_RUN <- FALSE

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

# The runner itself needs only base R, digest and the data-reading packages for S0; the
# drivers load everything else in their own R processes.
suppressWarnings(suppressPackageStartupMessages(
  try({ library(dplyr); library(tidyr); library(readr); library(lubridate)
        library(readxl); library(here); library(purrr); library(stringr); library(tibble) },
      silent = TRUE)))
invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), function(f) source(f)))
source(.here("run_config.R"))
BASE <- run_config
BASE_MODEL <- if (exists("model")) model else "both"

POOLED_RMD <- .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd")
GEAR_RMD   <- .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd")
ORCH       <- .here("run_estimation.R")
stopifnot(file.exists(POOLED_RMD), file.exists(GEAR_RMD), file.exists(ORCH))
OUT_PREFIX <- "authoritative_batch_2026-09-29"
LOG_DIR    <- .here("05_output", paste0(OUT_PREFIX, "_logs"))
RSCRIPT    <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

# EVERY STAGE'S R PROCESS USES THIS SESSION'S LIBRARY (2026-09-30, B59). A stage is a new
# Rscript started in the project folder, so the project's .Rprofile would activate renv in it
# and it would load packages from the renv library, whatever this session uses. On the first
# overnight attempt that library was incomplete and not the lockfile's (no tidyverse, rstan or
# loo; lubridate 1.9.5 against the lock's 1.9.1), while the Console that had rendered B50 and
# B57 used another: stage A spent 24 minutes compiling packages and died, every other stage
# died in seconds. So each stage is started with renv's autoloader off and with exactly this
# session's .libPaths() (the user and site libraries pointed at a path that does not exist, so
# R cannot add a library this session does not have). If this session runs under renv, the
# stages use the renv library; if not, they use the same library it does. Desk check S0 starts
# one such process and compares package versions before anything is fitted.
.nolib <- file.path(tempdir(), "no-library-here")
CHILD_ENV <- c(RENV_ACTIVATE_PROJECT = "FALSE",
               R_LIBS = paste(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), collapse = .Platform$path.sep),
               R_LIBS_USER = .nolib, R_LIBS_SITE = .nolib)
.with_child_env <- function(expr) {
  old <- Sys.getenv(names(CHILD_ENV), unset = NA)
  do.call(Sys.setenv, as.list(CHILD_ENV))
  on.exit(for (k in names(old)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k)),
          add = TRUE)
  force(expr)
}
# the packages whose versions must agree between this session and a stage's process: the
# required list, and the ones that decide how Stan compiles and what it computes
CHILD_PKGS <- unique(c(bss_required_packages, "StanHeaders", "Rcpp", "RcppEigen", "BH", "RcppParallel", "posterior"))

# The authoritative run this batch may replace, and whose components A4 reads against.
REF_AUTH_POOLED <- "20260928/pooled-CPUE-2024-25"
REF_AUTH_GEAR   <- "20260928/gear-type-CPUE-model-2024-25"

FITS <- c("shore_ring_net_only_Dungeness_Kept", "shore_all_gear_Dungeness_Kept",
          "private_boat_ring_net_only_Dungeness_Kept", "private_boat_all_gear_Dungeness_Kept")
FIT_SHORE_AG <- "shore_all_gear_Dungeness_Kept"
FIT_SHORE_PC <- "shore_ring_net_only_Dungeness_Kept"
FIT_BOAT_AG  <- "private_boat_all_gear_Dungeness_Kept"
COMP_KEY <- c(shore_ring_net_only_Dungeness_Kept = "shore (Pot closure)",
              shore_all_gear_Dungeness_Kept = "shore (All gear)",
              private_boat_ring_net_only_Dungeness_Kept = "private_boat (Pot closure)",
              private_boat_all_gear_Dungeness_Kept = "private_boat (All gear)")

# ---------------------------------------------------------------------------
# THE STAGES. Each is the shipped configuration plus a declared delta, nothing else.
# ---------------------------------------------------------------------------
# D3's per-population periods: the shore at the pooled caps, the boat at its shipped gear
# periods (see the header). Built from the configuration so it cannot drift from it.
MATCHED <- list(shore = BASE$ar_max_resolution$pooled$shore,
                private_boat = list(all_gear = BASE$gear_period_bss$all_gear %||% "month",
                                    pot_closure = BASE$gear_period_bss$pot_closure %||% "biweekly"))
LADDER <- list(ar_escalate = list(private_boat = "all_gear"), ar_escalate_stop = "all_rungs",
               ar_escalate_select = "first_pass", ar_escalate_ladder = c("weekly", "biweekly", "monthly"),
               ar_escalate_respect_cap = FALSE, ar_escalate_max_attempts = 3, ar_rung_adequacy = TRUE)
STAGE_DEFS <- list(
  S0   = list(model = NA, fit = FALSE, tag = NA, item = "desk: prerequisites, no MCMC", delta = list()),
  A    = list(model = "both", fit = TRUE, tag = BASE$run_tag, orchestrator = TRUE,
              item = "THE AUTHORITATIVE RENDER: run_estimation.R --model both, run_config.R as shipped", delta = list()),
  D3   = list(model = "gear_resolved", fit = TRUE, tag = "2024-25-AB-D3-gear-matched",
              item = "D3: gear track, per-population AR period (shore at the pooled caps; only shore all-gear moves)",
              delta = list(gear_period_bss = MATCHED)),
  D6   = list(model = "gear_resolved", fit = TRUE, tag = "2024-25-AB-D6-gear-matched-zi",
              item = "D6: D3 + the zero-inflated shore catch on the gear track",
              delta = list(gear_period_bss = MATCHED, catch_zi_tracks = c("pooled", "gear_resolved"))),
  D29P = list(model = "pooled", fit = TRUE, tag = "2024-25-AB-D29P-boat-ladder",
              item = "D29: pooled boat all-gear ladder weekly/biweekly/monthly (all rungs)",
              delta = LADDER),
  D29G = list(model = "gear_resolved", fit = TRUE, tag = "2024-25-AB-D29G-gear-boat-ladder",
              item = "D29: gear boat all-gear ladder weekly/biweekly/monthly, on D3's configuration",
              delta = c(list(gear_period_bss = MATCHED), LADDER)),
  R2   = list(model = "pooled", fit = TRUE, tag = "2024-25-AB-R2-init-r-2",
              item = "radius: the pooled A31 model at init_r = 2", delta = list(bss_init_r = 2)),
  D29D = list(model = "pooled", fit = TRUE, tag = "2024-25-AB-D29D-boat-daily",
              item = "D29: pooled boat all-gear at DAILY (ar_force)",
              delta = list(ar_force = list(private_boat = list(all_gear = "daily"))))
)
# Wall-clock limit per stage, in hours (system2(timeout = )); a stage that reaches it is
# stopped and recorded FAILED. Several times the expected time, so a slow fit is not killed.
TIMEOUT_H <- c(A = 10, D3 = 3, D6 = 3, D29P = 9, D29G = 5, R2 = 7, D29D = 10)
# The stages each desk check gates. A desk FAIL stops only these (and never the whole batch
# unless it is a prerequisite of every render: the method of record, the NWS archive, pandoc).
DEPENDS <- list(D3 = c("D3", "D6", "D29G"), D6 = "D6", D29P = "D29P", D29G = "D29G", R2 = "R2", D29D = "D29D")
if (!all(STAGES %in% names(STAGE_DEFS)))
  stop("STAGES names a stage that does not exist: ", paste(setdiff(STAGES, names(STAGE_DEFS)), collapse = ", "))
# keys the driver or orchestrator ADDS at run time (data, not configuration)
RUNTIME_KEYS <- c("run_tag", "model")

resolve_cfg <- function(sid) {
  cfg <- BASE
  d <- STAGE_DEFS[[sid]]$delta
  for (k in names(d)) cfg[k] <- list(d[[k]])     # list() so a NULL value is set, not dropped
  if (!is.na(STAGE_DEFS[[sid]]$tag %||% NA)) cfg$run_tag <- STAGE_DEFS[[sid]]$tag
  cfg
}

# ---------------------------------------------------------------------------
# DIGESTS AND STAMPS. RESUME may reuse a folder only when its stamp's digest matches.
# ---------------------------------------------------------------------------
digest_or_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The digest package is missing (it is in bss_required_packages): run renv::restore() first.", call. = FALSE)
  digest::digest(x)
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
  sprintf("stan:%s drivers:%s fns:%s orch:%s",
          .g(.here("02_stan_models")), .g(c(POOLED_RMD, GEAR_RMD)), .g(.here("03_R_functions")), .g(ORCH))
}
.inputs_fingerprint <- function() {
  fs <- sort(list.files(.here("04_input_files"), pattern = "[.]xlsx$", full.names = TRUE, recursive = TRUE))
  fs <- fs[!grepl("^~\\$", basename(fs))]     # Excel lock files of an open workbook are not inputs
  substr(digest_or_hash(paste(basename(fs), unname(tools::md5sum(fs)), collapse = "\n")), 1, 8)
}
.cfg_text <- function(cfg) {
  keys <- sort(setdiff(names(cfg), RUNTIME_KEYS))
  paste(vapply(keys, function(k) paste0(k, "=", paste(deparse(cfg[[k]]), collapse = "")), character(1)), collapse = ";")
}
CODE_FP <- .code_fingerprint(); INPUTS_FP <- .inputs_fingerprint()
stage_digest <- function(sid) {
  substr(digest_or_hash(paste(sid, STAGE_DEFS[[sid]]$model, .cfg_text(resolve_cfg(sid)), CODE_FP, INPUTS_FP,
                              sep = "\n")), 1, 12)
}
# NA when git did not answer (not installed, not a repository, a refused safe.directory):
# system2(stdout = TRUE) then returns what git printed, often nothing, with a status attribute,
# and "nothing printed" must never read as a clean tree.
.git <- function(args) {
  x <- tryCatch(suppressWarnings(system2("git", args, stdout = TRUE, stderr = FALSE)), error = function(e) NA_character_)
  if (!is.null(attr(x, "status")) || (length(x) == 1 && is.na(x))) NA_character_ else x
}
GIT_SHA <- { s <- .git(c("rev-parse", "--short", "HEAD")); if (length(s) == 1 && !is.na(s) && grepl("^[0-9a-f]{7,}$", s)) s else NA_character_ }
.stamp <- function(dir, sid, part, manifest = NA_character_) {
  writeLines(c(sprintf("stage: %s", sid), sprintf("part: %s", part), sprintf("digest: %s", stage_digest(sid)),
               sprintf("code: %s", CODE_FP), sprintf("inputs: %s", INPUTS_FP),
               sprintf("manifest: %s", if (is.na(manifest)) "none" else normalizePath(manifest, winslash = "/")),
               sprintf("git: %s", GIT_SHA %||% NA), sprintf("written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               "", "# run_authoritative_batch_2026-09-29.R writes this after a stage renders. RESUME reuses",
               "# a folder ONLY when its digest matches the stage being requested."),
             file.path(dir, "AB_STAGE.txt"))
}
.stamp_field <- function(dir, key) {
  p <- file.path(dir %||% "", "AB_STAGE.txt"); if (!file.exists(p)) return(NA_character_)
  l <- grep(paste0("^", key, ":"), readLines(p, warn = FALSE), value = TRUE)
  if (!length(l)) NA_character_ else trimws(sub(paste0("^", key, ":"), "", l[1]))
}
# Every run folder under 05_output/<date>/ that carries a stamp for this stage and part with
# the CURRENT digest; the newest wins.
find_stage_dir <- function(sid, part) {
  stamps <- Sys.glob(.here("05_output", "*", "*", "AB_STAGE.txt"))
  hits <- dirname(stamps)[vapply(dirname(stamps), function(d)
    identical(.stamp_field(d, "stage"), sid) && identical(.stamp_field(d, "part"), part) &&
      identical(.stamp_field(d, "digest"), stage_digest(sid)), logical(1))]
  if (!length(hits)) return(NA_character_)
  hits[order(file.mtime(file.path(hits, "AB_STAGE.txt")), decreasing = TRUE)][1]
}
PARTS <- function(sid) switch(STAGE_DEFS[[sid]]$model, both = c("pooled", "gear_resolved"), STAGE_DEFS[[sid]]$model)
PREFIX_OF <- c(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")

# ---------------------------------------------------------------------------
# VERDICT ROWS
# ---------------------------------------------------------------------------
V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = as.character(observed),
                                    threshold = threshold, verdict = verdict, why = why, stringsAsFactors = FALSE)
REC <- list()
rec_row <- function(item, recommendation, rule, evidence)
  REC[[length(REC) + 1]] <<- data.frame(item = item, recommendation = recommendation, rule = rule,
                                        evidence = evidence, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# S0: THE DESK STAGE
# ---------------------------------------------------------------------------
stage_S0 <- function() {
  banner("S0  DESK: prerequisites, before any MCMC")
  q <- function(e) { s <- tempfile(); sink(s); on.exit(sink()); force(e) }
  # the ladder's effort-day table, built exactly as both drivers build eff_d_ladder
  .eff_d <- function(dwg, pop, ss, d_, p) {
    sm <- q(prep_population_summary(dwg, pop, ss$start, ss$end, p))
    sm$effort_index |> dplyr::filter(count_sequence <= p$bss_max_count_seq) |>
      dplyr::left_join(d_ |> dplyr::select(event_date, day_index), by = "event_date") |> dplyr::filter(!is.na(day_index))
  }
  # (1) the shipped configuration is the method of record on main
  ship <- list(model = BASE_MODEL, bss_init_r = BASE$bss_init_r, mu_hier_collapse_single = BASE$mu_hier_collapse_single,
               marine_hazard_mode = BASE$marine_hazard_mode, marine_hazard_manual_boat = BASE$marine_hazard_manual_boat,
               marine_hazard_gear_regimes = BASE$marine_hazard_gear_regimes, catch_zi_tracks = BASE$catch_zi_tracks,
               ar_escalate = BASE$ar_escalate, ar_force = BASE$ar_force, bss_sampler_override = BASE$bss_sampler_override,
               season_filter = BASE$season_filter, est_date_start = BASE$est_date_start, est_date_end = BASE$est_date_end)
  want <- list(model = "both", bss_init_r = 0.5, mu_hier_collapse_single = "both", marine_hazard_mode = "manual",
               marine_hazard_manual_boat = "nws_sca_any", marine_hazard_gear_regimes = "all_gear", catch_zi_tracks = "pooled",
               ar_escalate = FALSE, ar_force = NULL, bss_sampler_override = NULL,
               season_filter = "2024-25", est_date_start = "2024-09-16", est_date_end = "2025-09-15")
  bad <- names(want)[!vapply(names(want), function(k) identical(ship[[k]], want[[k]]), logical(1))]
  V1row("S0", "run_config.R ships the method of record (the configuration stage A renders)",
        if (length(bad)) paste("differs:", paste(sprintf("%s = %s", bad, vapply(bad, function(k) paste(deparse(ship[[k]]), collapse = ""), character(1))), collapse = "; ")) else "as expected",
        paste(sprintf("%s = %s", names(want), vapply(want, function(v) paste(deparse(v), collapse = ""), character(1))), collapse = "; "),
        if (length(bad)) "FAIL" else "PASS",
        "Stage A renders run_config.R as shipped; if it is not the method of record, A is not an authoritative render.")
  # (2) each stage differs from A's configuration in its declared keys only
  base_cfg <- resolve_cfg("A")
  for (sid in setdiff(names(STAGE_DEFS), c("S0", "A"))) {
    cfg <- resolve_cfg(sid)
    ks <- setdiff(union(names(base_cfg), names(cfg)), RUNTIME_KEYS)
    diff <- ks[!vapply(ks, function(k) identical(cfg[[k]], base_cfg[[k]]), logical(1))]
    undeclared <- setdiff(diff, names(STAGE_DEFS[[sid]]$delta))
    V1row(sid, "differs from the shipped configuration in DECLARED keys only",
          sprintf("differs in: %s", if (length(diff)) paste(diff, collapse = ", ") else "nothing"),
          sprintf("a subset of {%s}", paste(names(STAGE_DEFS[[sid]]$delta), collapse = ", ")),
          if (!length(undeclared)) "PASS" else "FAIL",
          paste("A stage that differs in an undeclared key measures that key as well as its own. Undeclared:",
                if (length(undeclared)) paste(undeclared, collapse = ", ") else "none"))
  }
  V1row("D6", "differs from D3 in catch_zi_tracks ONLY",
        { a <- resolve_cfg("D3"); b <- resolve_cfg("D6"); ks <- setdiff(union(names(a), names(b)), RUNTIME_KEYS)
          paste(ks[!vapply(ks, function(k) identical(a[[k]], b[[k]]), logical(1))], collapse = ", ") },
        "catch_zi_tracks", { a <- resolve_cfg("D3"); b <- resolve_cfg("D6"); ks <- setdiff(union(names(a), names(b)), RUNTIME_KEYS)
          if (identical(ks[!vapply(ks, function(k) identical(a[[k]], b[[k]]), logical(1))], "catch_zi_tracks")) "PASS" else "FAIL" },
        "D6 is judged against D3; any second difference would be measured as ZI.")
  V1row("R2", "differs from A in bss_init_r ONLY", paste(names(STAGE_DEFS$R2$delta), collapse = ", "), "bss_init_r",
        if (identical(names(STAGE_DEFS$R2$delta), "bss_init_r")) "PASS" else "FAIL", "R2 is judged against A.")
  # (3) pandoc: every stage renders an .Rmd, and rmarkdown stops at once without it
  pd <- tryCatch(rmarkdown::pandoc_available("1.12.3"), error = function(e) FALSE)
  V1row("S0", "rmarkdown finds pandoc (every stage renders an .Rmd)",
        if (isTRUE(pd)) sprintf("pandoc %s at %s", as.character(rmarkdown::pandoc_version()), rmarkdown::pandoc_exec()) else "NOT FOUND",
        "pandoc >= 1.12.3 reachable", if (isTRUE(pd)) "PASS" else "FAIL",
        "Run from the RStudio Console (it sets RSTUDIO_PANDOC and the stages inherit it), set RSTUDIO_PANDOC to RStudio's pandoc folder, or install pandoc.")
  # (3b) a stage's R process sees THIS session's packages, loadable, at the same versions
  #      (B59: the first overnight attempt failed here, in every stage, after S0 had passed)
  tryCatch({
    cs <- file.path(tempdir(), "ab_child_packages.R"); co <- file.path(tempdir(), "ab_child_packages.txt")
    writeLines(c("a <- commandArgs(trailingOnly = TRUE); pk <- strsplit(a[1], ',', fixed = TRUE)[[1]]",
                 "v <- vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) 'MISSING'), '')",
                 "ld <- vapply(pk, function(p) isTRUE(suppressWarnings(requireNamespace(p, quietly = TRUE))), TRUE)",
                 "writeLines(c(paste0('renv\t', Sys.getenv('RENV_PROJECT')), paste0('libs\t', paste(normalizePath(.libPaths(), winslash = '/'), collapse = ';')),",
                 "             paste(pk, v, ld, sep = '\t')), a[2])"), cs)
    if (file.exists(co)) file.remove(co)
    stc <- .with_child_env(suppressWarnings(system2(RSCRIPT, c(shQuote(cs), shQuote(paste(CHILD_PKGS, collapse = ",")), shQuote(co)),
                                                    stdout = TRUE, stderr = TRUE, timeout = 300)))
    if (!file.exists(co)) stop("the test process wrote nothing: ", paste(utils::tail(stc, 5), collapse = " | "))
    cl <- strsplit(readLines(co, warn = FALSE), "\t", fixed = TRUE)
    kid <- do.call(rbind, lapply(cl[-(1:2)], function(x) data.frame(pkg = x[1], v = x[2], loads = identical(x[3], "TRUE"))))
    par_v <- vapply(CHILD_PKGS, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "MISSING"), "")
    need <- kid$pkg %in% bss_required_packages
    bad_load <- kid$pkg[need & !kid$loads]
    bad_ver  <- kid$pkg[kid$v != par_v[kid$pkg]]
    okk <- !length(bad_load) && !length(bad_ver)
    V1row("S0", "each stage's R process loads this session's packages at the same versions (renv autoloader off, this session's library)",
          sprintf("%s. Stage process: renv %s; libraries %s", if (okk) sprintf("all %d agree and load", nrow(kid)) else
                    paste(c(if (length(bad_load)) sprintf("do not load: %s", paste(bad_load, collapse = ", ")),
                            if (length(bad_ver)) sprintf("versions differ: %s", paste(sprintf("%s %s here, %s in the stage", bad_ver, par_v[bad_ver], kid$v[match(bad_ver, kid$pkg)]), collapse = "; "))),
                          collapse = "; "),
                  { rv <- cl[[1]][2]; if (!is.na(rv) && nzchar(rv)) sprintf("ACTIVE (%s)", rv) else "off" }, cl[[2]][2] %||% "?"),
          sprintf("%s load; %s same version", paste(bss_required_packages, collapse = ", "), paste(CHILD_PKGS, collapse = ", ")),
          if (okk) "PASS" else "FAIL",
          "Every stage renders in its own Rscript. If this fails, fix the library of THIS session (or run renv::restore() and restart R) before --go.")
  }, error = function(e) V1row("S0", "a stage's R process could be started to check its packages", conditionMessage(e), "no error", "FAIL",
                               "Every stage is such a process."))
  # (4) D3's per-population periods, fit by fit: the shore at the pooled caps, the boat at the
  #     shipped gear periods, so exactly one fit (shore all-gear) moves against the shipped gear track
  n <- function(x) .bss_normalize_resolution(x %||% NA)
  keys <- c("shore/pot_closure", "shore/all_gear", "private_boat/pot_closure", "private_boat/all_gear")
  got_per <- vapply(keys, function(k) { pr <- strsplit(k, "/")[[1]]; n(bss_gear_period(resolve_cfg("D3"), pr[1], pr[2])) }, character(1))
  pc <- BASE$ar_max_resolution$pooled$shore
  exp_per <- c("shore/pot_closure" = n(pc$pot_closure), "shore/all_gear" = n(pc$all_gear),
               "private_boat/pot_closure" = n(BASE$gear_period_bss$pot_closure), "private_boat/all_gear" = n(BASE$gear_period_bss$all_gear))
  ship_per <- c("shore/pot_closure" = n(BASE$gear_period_bss$pot_closure), "shore/all_gear" = n(BASE$gear_period_bss$all_gear),
                "private_boat/pot_closure" = n(BASE$gear_period_bss$pot_closure), "private_boat/all_gear" = n(BASE$gear_period_bss$all_gear))
  moved <- keys[got_per != ship_per]
  V1row("D3", "D3's gear periods: the shore at the pooled caps, the boat at its shipped periods; only shore all-gear moves",
        sprintf("%s; moved against the shipped gear track: %s", paste(sprintf("%s %s", keys, got_per), collapse = "; "),
                if (length(moved)) paste(moved, collapse = ", ") else "nothing"),
        sprintf("%s; moved: shore/all_gear", paste(sprintf("%s %s", keys, exp_per), collapse = "; ")),
        if (identical(unname(got_per), unname(exp_per)) && identical(moved, "shore/all_gear")) "PASS" else "FAIL",
        "bss_gear_period() is what the gear driver calls; D3 must move one fit, the one the register asks about.")
  # (5) the real inputs load (a prerequisite of every stage)
  inp <- tryCatch({
    p <- modifyList(resolve_cfg("D29P"), list(bss_model_file = "crab_bss_pooled.stan", boat_require_gear_time = TRUE))
    p$ar_max_resolution <- p$ar_max_resolution$pooled
    p$crabbing_holiday_dates <- read_crabbing_holidays(p)
    list(p = p, dwg = q(fetch_crab_data(p)), sub = build_subseasons(p))
  }, error = function(e) {
    V1row("S0", "the 2024-25 inputs load", conditionMessage(e), "no error", "FAIL", "Every stage reads them.")
    NULL
  })
  ok <- !is.null(inp)
  if (ok) {
    p <- inp$p; dwg <- inp$dwg; sub <- inp$sub
    # (5a) the pooled ladder
    tryCatch({
      lad <- list()
      for (ss in sub) for (pop in c("shore", "private_boat")) {
        d_ <- q(prep_days_crab(ss$start, ss$end, p, L_eff_model = NULL))
        lad[[sprintf("%s/%s", pop, ss$gear_regime)]] <- paste(bss_ar_ladder(d_, .eff_d(dwg, pop, ss, d_, p), pop, p, gear_regime = ss$gear_regime), collapse = ">")
      }
      V1row("D29P", "the pooled ladder: weekly>biweekly>monthly for the boat all-gear fit, ONE rung for every other fit",
            paste(sprintf("%s %s", names(lad), unlist(lad)), collapse = "; "), "boat all-gear weekly>biweekly>monthly; others a single rung",
            if (identical(lad[["private_boat/all_gear"]], "weekly>biweekly>monthly") &&
                all(!grepl(">", unlist(lad[setdiff(names(lad), "private_boat/all_gear")])))) "PASS" else "FAIL",
            "bss_ar_ladder() on the real 2024-25 days; a scoped ar_escalate must not reach the other fits.")
    }, error = function(e) V1row("D29P", "the pooled ladder resolves on the real inputs", conditionMessage(e), "no error", "FAIL", ""))
    # (5b) the gear ladder on D29G's configuration
    tryCatch({
      pg <- resolve_cfg("D29G"); pg$crabbing_holiday_dates <- p$crabbing_holiday_dates
      lg <- list()
      for (ss in sub) for (pop in c("shore", "private_boat")) {
        d_ <- q(prep_days_crab(ss$start, ss$end, pg, L_eff_model = NULL))
        per_ss <- bss_gear_period(pg, pop, ss$gear_regime) %||% ss$period_bss
        lg[[sprintf("%s/%s", pop, ss$gear_regime)]] <- paste(bss_ar_ladder(d_, .eff_d(dwg, pop, ss, d_, pg), pop, pg, fixed_resolution = per_ss,
                                                                          gear_regime = ss$gear_regime), collapse = ">")
      }
      lg_exp <- c(exp_per[c("shore/pot_closure", "shore/all_gear", "private_boat/pot_closure")], "private_boat/all_gear" = "weekly>biweekly>monthly")
      lg_got <- vapply(names(lg_exp), function(k) .bss_normalize_resolution_chain(lg[[k]]), character(1))
      V1row("D29G", "the gear ladder: weekly>biweekly>monthly for the boat all-gear fit, D3's period for every other fit",
            paste(sprintf("%s %s", names(lg_got), lg_got), collapse = "; "), paste(sprintf("%s %s", names(lg_exp), lg_exp), collapse = "; "),
            if (identical(unname(lg_got), unname(lg_exp))) "PASS" else "FAIL",
            "The gear driver passes fixed_resolution = the per-population period (ar_adaptive = FALSE); the ladder must override it for the boat all-gear fit only.")
    }, error = function(e) V1row("D29G", "the gear ladder resolves on the real inputs", conditionMessage(e), "no error", "FAIL", ""))
    # (5c) D6: the gear prep builds zi_catch = 1 for the shore all-gear fit and 0 for the boat all-gear fit
    tryCatch({
      ie <- q(fetch_ie_data(p))
      p6 <- modifyList(resolve_cfg("D6"), list(bss_model_file = "crab_bss_gear_resolved.stan", boat_require_gear_time = FALSE))
      p6$crabbing_holiday_dates <- p$crabbing_holiday_dates; p6$ar_max_resolution <- p6$ar_max_resolution$gear_resolved
      p6$crab_fraction_rows <- q(crab_fraction_source_rows(dwg, ie, p6))
      stt <- q(estimate_shore_turnover(attr(ie, "ie_intervals"), dwg$shore_effort, p6))
      osp <- tryCatch(q(fetch_osp_boat_counts(p6)), error = function(e) NULL)
      ov  <- tryCatch(q(diagnose_osp_trailer_overlap(osp, p6, output_dir = NULL)), error = function(e) NULL)
      p6$osp_crab_rows <- attr(osp, "osp_crab_rows")
      p6 <- q(bss_resolve_tau_boat_prior(p6, ov)); p6 <- q(bss_resolve_tau_shore_prior(p6, stt))
      zi <- list()
      for (pop in c("shore", "private_boat")) {
        ss <- sub[[which(vapply(sub, function(x) x$gear_regime == "all_gear", logical(1)))[1]]]
        d_ <- q(prep_days_crab(ss$start, ss$end, p6, L_eff_model = NULL))
        sm <- q(prep_population_summary(dwg, pop, ss$start, ss$end, p6))
        per_ss <- bss_gear_period(p6, pop, ss$gear_regime) %||% ss$period_bss
        bd <- q(prep_bss_crab_gear(d_, sm, "Dungeness_Kept", p6, pop, per_ss, gear_exclude = ss$gear_exclude %||% character(0),
                                   gear_regime = ss$gear_regime, ie_data = ie))
        zi[[pop]] <- bd$zi_catch %||% NA
      }
      V1row("D6", "the gear prep under D6 builds zi_catch = 1 for the shore all-gear fit and 0 for the boat all-gear fit",
            sprintf("shore %s; private_boat %s", zi$shore, zi$private_boat), "shore 1; private_boat 0",
            if (identical(as.integer(zi$shore), 1L) && identical(as.integer(zi$private_boat), 0L)) "PASS" else "FAIL",
            "catch_zi_populations = 'shore' scopes the likelihood; catch_zi_tracks is the only key D6 moves. (The pot-closure pair is scoped the same way.)")
    }, error = function(e) V1row("D6", "the D6 gear Stan data build on the real inputs", conditionMessage(e), "no error", "FAIL", ""))
    # (5d) the NWS archive covers the window for BOTH zones (every production run stops otherwise)
    tryCatch({
      ev <- marine_hazard_events(p); cov <- attr(ev, "coverage")
      zones <- unname(BASE$marine_hazard_zones)
      okc <- !is.null(cov) && nrow(cov) == length(zones) && all(zones %in% cov$ugc) &&
             all(cov$pull_start <= as.Date(BASE$est_date_start)) && all(cov$pull_end >= as.Date(BASE$est_date_end))
      V1row("S0", "the NWS archive covers the window for both zones",
            if (is.null(cov)) "no coverage table" else paste(sprintf("%s %s to %s", cov$ugc, cov$pull_start, cov$pull_end), collapse = "; "),
            sprintf("%s: %s to %s", paste(zones, collapse = " and "), BASE$est_date_start, BASE$est_date_end),
            if (okc) "PASS" else "FAIL", "The method of record reads nws_marine_hazards.xlsx and stops on an uncovered window.")
    }, error = function(e) V1row("S0", "the NWS archive was read", conditionMessage(e), "no error", "FAIL", ""))
  }
  # (6) the reference folders (A4 and D3-4 read them; nothing renders from them, so REVIEW, not FAIL)
  for (r in c(REF_AUTH_POOLED, REF_AUTH_GEAR))
    V1row("S0", sprintf("the reference render %s is present", r), if (dir.exists(.here("05_output", r))) "present" else "MISSING",
          "present", if (dir.exists(.here("05_output", r))) "PASS" else "REVIEW", "A4 and D3-4 read it; without it they report NA.")
  # (7) the per-rung block CV the D29 rule reads is wired into each driver (gates that track's ladder)
  for (f in c(POOLED_RMD, GEAR_RMD)) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    w <- grepl('prefix = "ladder_block"', src, fixed = TRUE) && grepl("cov50_trailer", src, fixed = TRUE)
    V1row(if (identical(f, POOLED_RMD)) "D29P" else "D29G", sprintf("%s writes every ladder rung's block CV (ladder_block_*) and trailer coverage", basename(f)),
          if (w) "yes" else "NO", "both present", if (w) "PASS" else "FAIL",
          "Without it the D29 rule's interpolation clause has nothing to read for the rungs the ladder did not keep.")
  }
  # (8) git: the commit, and whether the tracked tree is clean (the authoritative manifest records it)
  dirty <- .git(c("status", "--porcelain", "--untracked-files=no"))
  V1row("S0", "the tracked tree is clean (stage A's manifest records any tracked file that differs)",
        if (is.na(GIT_SHA) || (length(dirty) == 1 && is.na(dirty))) "git did not answer" else if (!length(dirty)) sprintf("clean at %s", GIT_SHA)
        else paste(trimws(substring(dirty, 4)), collapse = ", "),
        "clean", if (!is.na(GIT_SHA) && !(length(dirty) == 1 && is.na(dirty)) && !length(dirty)) "PASS" else "REVIEW",
        "Start with --go, not by editing DRY_RUN: an edited runner is a tracked file that differs. Not blocking; A6 reads the manifest.")
  invisible(ok)
}
.bss_normalize_resolution_chain <- function(x) {
  if (is.null(x) || is.na(x)) return(NA_character_)
  paste(vapply(strsplit(x, ">", fixed = TRUE)[[1]], .bss_normalize_resolution, character(1)), collapse = ">")
}

# ---------------------------------------------------------------------------
# RUNNING ONE STAGE, each in its own R process
# ---------------------------------------------------------------------------
.child_script <- function() {
  f <- file.path(tempdir(), "ab_render_child.R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "setwd(args[1]); cfg <- readRDS(args[2]); rmd <- args[3]; out_file <- args[4]",
    "run_env <- new.env(parent = globalenv()); run_env$run_config <- cfg",
    "html <- rmarkdown::render(rmd, envir = run_env, quiet = FALSE)",
    "od <- get('output_dir', envir = run_env, inherits = FALSE)",
    "if (dir.exists(od) && file.exists(html) && normalizePath(dirname(html)) != normalizePath(od)) {",
    "  if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE))) invisible(file.remove(html))",
    "}",
    "writeLines(normalizePath(od, winslash = '/'), out_file)"), f)
  f
}
# The run folder a manifest line names, mapped into THIS checkout's 05_output (the manifest
# holds the absolute path of the machine that rendered).
.map_outdir <- function(p) {
  p <- gsub("[\\\\/]+", "/", trimws(p)); parts <- strsplit(p, "/", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NA_character_)
  d <- .here("05_output", parts[length(parts) - 1], parts[length(parts)])
  if (dir.exists(d)) d else NA_character_
}
# The stage lines of run_estimation.R's manifest. Every line after "Stages:" that names a model
# is read (a FAILED line's error message can run over several lines, so the block is not
# assumed to end at the first line that does not look like a stage).
.parse_manifest <- function(path) {
  l <- readLines(path, warn = FALSE)
  i <- which(l == "Stages:"); if (!length(i) || i[1] >= length(l)) return(NULL)
  out <- list()
  for (x in l[(i[1] + 1):length(l)]) {
    if (grepl("^run_config \\(", x)) break
    if (!grepl("^  (pooled|gear_resolved)\\s", x)) next
    nm <- sub("^  (\\S+).*$", "\\1", x)
    failed <- grepl(" FAILED ", x, fixed = TRUE)
    path <- if (failed) sub("^.*\\(partial folder: (.*)\\)$", "\\1", x) else sub("^  \\S+\\s+\\S+ min\\s+(.*)$", "\\1", x)
    out[[nm]] <- list(failed = failed, outdir = if (!failed || grepl("partial folder", x)) .map_outdir(path) else NA_character_,
                      minutes = .num1(sub("^  \\S+\\s+([0-9.]+) min.*$", "\\1", x)))
  }
  out
}
STAGE_STATUS <- list()
# system2() with a wall-clock limit. Returns the exit status; 124 when the limit stopped it
# (R's own convention for system2(timeout =)). The limit stops the stage's R process; on
# Windows rstan's PSOCK workers exit when their master's connection closes.
.sys <- function(args, log, hours) {
  st <- tryCatch(.with_child_env(suppressWarnings(system2(RSCRIPT, args, stdout = log, stderr = log, timeout = round(3600 * hours)))),
                 error = function(e) { cat("  system2 error:", conditionMessage(e), "\n"); -1L })
  as.integer(st %||% -1L)
}
# The cross-check the orchestrator writes after a "both" render, rebuilt by the runner when
# stage A's two models were rendered by separate invocations (a resumed A). Same columns, same
# formula (run_estimation.R section 7), written beside the pooled folder.
.write_cross_check <- function(dp, dg) {
  rdp <- function(d) { pt <- utils::read.csv(file.path(d, "port_total_Dungeness_Kept.csv"), stringsAsFactors = FALSE)
    pt[pt$Estimate == "Expected_Catch", c("BSS_median", "BSS_lo95", "BSS_hi95")] }
  p <- rdp(dp); g <- rdp(dg); tol <- as.numeric(BASE$cross_check_tolerance %||% 0.02)
  rel <- (g$BSS_median - p$BSS_median) / p$BSS_median
  cc <- data.frame(pooled_median = p$BSS_median, pooled_lo95 = p$BSS_lo95, pooled_hi95 = p$BSS_hi95,
                   gear_median = g$BSS_median, gear_lo95 = g$BSS_lo95, gear_hi95 = g$BSS_hi95,
                   gear_minus_pooled_pct = round(100 * rel, 2), tolerance_pct = 100 * tol,
                   verdict = if (abs(rel) <= tol) "PASS" else "REVIEW",
                   pooled_folder = basename(dp), gear_folder = basename(dg))
  f <- file.path(dirname(dp), sprintf("cross_check_%s_batch-resumed.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))
  utils::write.csv(cc, f, row.names = FALSE)
  f
}
run_stage <- function(sid) {
  st <- STAGE_DEFS[[sid]]
  rule()
  cat(sprintf("  %-5s %s\n        digest %s\n", sid, st$item, stage_digest(sid)))
  have <- vapply(PARTS(sid), function(pt) find_stage_dir(sid, pt), character(1))
  if (isTRUE(RESUME) && all(!is.na(have))) {
    cat(sprintf("  RESUME: %s present with a MATCHING digest - skipping the fit.\n", paste(basename(have), collapse = ", ")))
    STAGE_STATUS[[sid]] <<- list(status = "reused", dirs = have, minutes = NA_real_)
    return(have)
  }
  # the parts still to render (a resumed "both" stage renders only the missing model)
  todo <- if (isTRUE(RESUME)) names(have)[is.na(have)] else names(have)
  if (isTRUE(DRY_RUN)) {
    cat(sprintf("  DRY_RUN: not fitting. Would render: %s%s\n", paste(todo, collapse = " + "),
                if (length(todo) < length(have)) sprintf(" (reusing %s)", paste(basename(have[!is.na(have)]), collapse = ", ")) else ""))
    STAGE_STATUS[[sid]] <<- list(status = "planned", dirs = have, minutes = NA_real_)
    return(have)
  }
  dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
  log <- file.path(LOG_DIR, sprintf("%s_%s.log", sid, format(Sys.time(), "%Y%m%d-%H%M%S")))
  hours <- unname(TIMEOUT_H[sid]) %||% 12
  t0 <- Sys.time()
  dirs <- if (isTRUE(RESUME)) have else stats::setNames(rep(NA_character_, length(have)), names(have))
  man <- NA_character_; note <- character(0)
  if (isTRUE(st$orchestrator)) {
    mdl <- if (length(todo) == length(have)) st$model else todo
    if (length(todo) < length(have))
      cat(sprintf("  PARTIAL RESUME: %s already rendered at this digest; rendering %s only.\n",
                  paste(basename(have[!is.na(have)]), collapse = ", "), mdl))
    before <- Sys.glob(.here("05_output", "*", "run_manifest_*.txt"))
    status <- .sys(c(shQuote(ORCH), "--model", mdl), log, hours)
    after <- setdiff(Sys.glob(.here("05_output", "*", "run_manifest_*.txt")), before)
    man <- if (length(after)) after[order(file.mtime(after), decreasing = TRUE)][1] else NA_character_
    mp <- if (!is.na(man)) .parse_manifest(man) else NULL
    for (pt in todo) {
      s <- mp[[pt]]
      d <- if (is.null(s) || isTRUE(s$failed)) NA_character_ else s$outdir %||% NA_character_
      dirs[[pt]] <- d
      if (!is.na(d)) .stamp(d, sid, pt, man)
    }
    # a resumed "both": the orchestrator rendered one model, so it wrote no cross-check
    if (identical(st$model, "both") && length(todo) < length(have) && all(!is.na(dirs)))
      tryCatch({ f <- .write_cross_check(dirs[["pooled"]], dirs[["gear_resolved"]]); note <- sprintf("cross-check written by the runner: %s", f) },
               error = function(e) note <<- sprintf("cross-check NOT written: %s", conditionMessage(e)))
  } else {
    cfg <- resolve_cfg(sid); cfg$model <- st$model
    cfg_file <- file.path(tempdir(), sprintf("ab_cfg_%s.rds", sid)); saveRDS(cfg, cfg_file)
    out_file <- file.path(tempdir(), sprintf("ab_out_%s.txt", sid)); if (file.exists(out_file)) file.remove(out_file)
    rmd <- if (identical(st$model, "pooled")) POOLED_RMD else GEAR_RMD
    status <- .sys(c(shQuote(.child_script()), shQuote(.root), shQuote(cfg_file), shQuote(rmd), shQuote(out_file)), log, hours)
    od <- if (file.exists(out_file)) .map_outdir(readLines(out_file, warn = FALSE)[1]) else NA_character_
    if (!identical(status, 0L)) od <- NA_character_
    dirs[[st$model]] <- od
    if (!is.na(od)) .stamp(od, sid, st$model)
  }
  ok <- all(!is.na(dirs)) && identical(status, 0L)
  STAGE_STATUS[[sid]] <<- list(status = if (ok) "rendered" else if (identical(status, 124L)) "TIMED OUT" else "FAILED",
                               dirs = dirs, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                               manifest = man, exit = status, log = log, note = note)
  s <- STAGE_STATUS[[sid]]
  cat(sprintf("  %s %s in %.1f min -> %s   (log: %s)\n", sid, s$status, s$minutes,
              paste(ifelse(is.na(s$dirs), "NO FOLDER", basename(s$dirs)), collapse = ", "), log))
  if (length(note)) cat("  ", note, "\n", sep = "")
  if (!ok)
    V1row(sid, "the stage rendered",
          sprintf("exit status %s%s; see %s", s$exit %||% "?", if (identical(status, 124L)) sprintf(" (stopped at the %s h limit)", hours) else "", log),
          "exit 0 and an output folder", "FAIL", "The batch continued with the next stage. --go again re-runs what is missing.")
  s$dirs
}

# ---------------------------------------------------------------------------
# READERS
# ---------------------------------------------------------------------------
rd <- function(dir, f) {
  p <- file.path(dir %||% "", f); if (is.na(dir %||% NA) || !file.exists(p)) return(NULL)
  tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}
.gate  <- function(dir) rd(dir, "convergence_report.csv")
.port  <- function(dir) { x <- rd(dir, "port_total_Dungeness_Kept.csv"); if (is.null(x)) return(NULL)
  i <- c(which(x$Estimate == "Expected_Catch"), which(x$Estimate == "Catch")); if (!length(i)) NULL else as.list(x[i[1], , drop = FALSE]) }
.draws <- function(dir, fit) { x <- rd(dir, sprintf("bss_draws_summed_%s.csv", fit)); if (is.null(x) || !"C_expected_sum" %in% names(x)) NULL else x$C_expected_sum }
.comp  <- function(dir, fit) {
  d <- .draws(dir, fit); g <- .gate(dir)
  m <- if (!is.null(g)) g$method_selected[g$fit == fit] else character(0)
  list(median = if (length(d)) stats::median(d) else NA_real_, sd = if (length(d) > 1) stats::sd(d) else NA_real_,
       lo = if (length(d)) unname(stats::quantile(d, 0.025)) else NA_real_, hi = if (length(d)) unname(stats::quantile(d, 0.975)) else NA_real_,
       method = if (length(m)) m[1] else NA_character_)
}
.adq <- function(dir, fit) { x <- rd(dir, "model_adequacy.csv"); if (is.null(x)) return(NULL)
  i <- which(x$fit == fit); if (!length(i)) NULL else as.list(x[i[1], , drop = FALSE]) }
.nobs <- function(dir, fit, streams) { x <- rd(dir, sprintf("loo_summary_%s.csv", fit)); if (is.null(x)) return(NA_real_)
  sum(x$n_obs[x$stream %in% streams], na.rm = TRUE) }
.samp <- function(dir, fit) rd(dir, sprintf("sampler_diagnostics_%s.csv", fit))
.stuck <- function(dir, fit) {
  s <- .samp(dir, fit); g <- .gate(dir); tp <- if (!is.null(g)) .num1(g$treedepth_pct[g$fit == fit]) else NA_real_
  if (is.null(s) || !"mean_stepsize" %in% names(s)) return(list(stuck = NA, why = "sampler_diagnostics missing", treedepth_pct = tp))
  med <- stats::median(s$mean_stepsize); r <- s$mean_stepsize / med
  list(stuck = any(r < 0.1, na.rm = TRUE) || isTRUE(tp >= 5),
       why = sprintf("step-size ratio to median %s; treedepth saturation %s%%", paste(fmt(r, 3), collapse = "/"), fmt(tp, 1)),
       treedepth_pct = tp)
}
.arlog <- function(dir) rd(dir, "ar_escalation_log.csv")
.bin <- function(dir, fit, k) {
  d <- rd(dir, sprintf("ppc_byobs_%s.csv", fit)); col <- if (k == 0L) "p_zero" else "p_one"
  if (is.null(d) || !all(c("data_type", col, "observed") %in% names(d))) return(NA_real_)
  y <- d[d$data_type == "catch", , drop = FALSE]; p <- suppressWarnings(as.numeric(y[[col]])); ok <- is.finite(p)
  obs <- sum(y$observed[ok] == k, na.rm = TRUE); ex <- sum(p[ok]); sdv <- sqrt(sum(p[ok] * (1 - p[ok])))
  if (is.finite(sdv) && sdv > 0) (obs - ex) / sdv else NA_real_
}
.theta <- function(dir, fit) { x <- rd(dir, sprintf("bss_full_summary_%s.csv", fit)); if (is.null(x)) return(NA_real_)
  nm <- names(x)[1]; i <- which(x[[nm]] %in% c("theta_C_out", "theta_C[1]")); if (!length(i)) NA_real_ else .num1(x$`50%`[i[1]]) }
.block <- function(dir, fname) { x <- rd(dir, fname); if (is.null(x) || !"elpd_block" %in% names(x)) NULL else x }

# ---------------------------------------------------------------------------
# THE VERDICTS
# ---------------------------------------------------------------------------
COMP <- list()
# What each stage's numbers ARE, printed beside them so a ladder's kept fit is never read as a
# reference estimate.
COMP_NOTE <- c(A = "the render of run_config.R as shipped (candidate authoritative run)",
               D3 = "gear track, per-population periods (the D3 candidate)",
               D6 = "D3 + zero-inflated shore catch on the gear track",
               D29P = "ladder stage: the boat all-gear fit is the rung the ladder KEPT (finest passing), not a reference; every rung is in the d29_rungs table",
               D29G = "ladder stage on D3's configuration: the boat all-gear fit is the KEPT rung, not a reference",
               R2 = "init_r = 2 (radius test), not a candidate estimate",
               D29D = "boat all-gear forced to daily (D29), not a candidate estimate")
record_components <- function(sid, dir, part) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  pt <- .port(dir); nt <- unname(COMP_NOTE[sid]) %||% ""
  for (f in FITS) { c1 <- .comp(dir, f)
    COMP[[length(COMP) + 1]] <<- data.frame(stage = sid, track = part, component = f, method = c1$method,
      median = round(c1$median), sd = round(c1$sd, 1), lo95 = round(c1$lo), hi95 = round(c1$hi), folder = basename(dir), note = nt,
      stringsAsFactors = FALSE) }
  COMP[[length(COMP) + 1]] <<- data.frame(stage = sid, track = part, component = "PORT (expected catch)", method = "BSS where the gate passed, else PE",
    median = round(.num1(pt$BSS_median)), sd = NA_real_, lo95 = round(.num1(pt$BSS_lo95)), hi95 = round(.num1(pt$BSS_hi95)), folder = basename(dir), note = nt,
    stringsAsFactors = FALSE)
}
.sd_units <- function(new, ref) if (isTRUE(is.finite(new$median)) && isTRUE(is.finite(ref$median)) && isTRUE(new$sd > 0)) (new$median - ref$median) / new$sd else NA_real_

# Three-valued clauses. TRUE -> PASS, FALSE -> FAIL, NA (a file missing, a number that could
# not be computed) -> REVIEW. A clause that could not be read never counts as a pass OR a fail.
.tri <- function(ok) if (isTRUE(ok)) "PASS" else if (identical(ok, FALSE)) "FAIL" else "REVIEW"
.worst <- function(vs) if (any(vs == "FAIL")) "FAIL" else if (any(vs == "REVIEW")) "REVIEW" else "PASS"
# every one of n fits reports BSS: NA when the report is missing, FALSE when a fit is missing
# from it (a component that was not fitted is not a BSS estimate) or any fit failed
.gate_all <- function(g, n = 4) {
  if (is.null(g)) return(NA)
  pc <- as.logical(g$pass_convergence)
  if (nrow(g) < n) return(FALSE)
  if (any(is.na(pc))) return(if (any(pc %in% FALSE)) FALSE else NA)
  all(pc)
}
.gate_str <- function(g) if (is.null(g)) "convergence_report.csv missing" else
  paste(sprintf("%s %s", sub("_Dungeness_Kept", "", g$fit), g$method_selected), collapse = ", ")
.stuck_ok <- function(st) { s <- vapply(st, function(x) if (is.na(x$stuck)) NA else !x$stuck, logical(1))
  if (any(s %in% FALSE)) FALSE else if (any(is.na(s))) NA else TRUE }
# The manifest that recorded this folder: the one its stamp names, else the newest manifest in
# its dated folder whose stage lines name it (a stamp written before the field existed).
.manifest_for <- function(dir) {
  m <- .stamp_field(dir, "manifest")
  if (!is.na(m) && !identical(m, "none") && file.exists(m)) return(m)
  ms <- Sys.glob(file.path(dirname(dir), "run_manifest_*.txt"))
  for (f in ms[order(file.mtime(ms), decreasing = TRUE)]) {
    mp <- tryCatch(.parse_manifest(f), error = function(e) NULL)
    if (any(vapply(mp, function(s) identical(basename(s$outdir %||% ""), basename(dir)) && !isTRUE(s$failed), logical(1)))) return(f)
  }
  NA_character_
}
.manifest_tree <- function(dir) {
  m <- .manifest_for(dir); if (is.na(m)) return(NA_character_)
  l <- grep("^git tree", readLines(m, warn = FALSE), value = TRUE)
  if (!length(l)) NA_character_ else trimws(sub("^git tree\\s*:", "", l[1]))
}

verdict_A <- function(dirs) {
  dp <- dirs[["pooled"]]; dg <- dirs[["gear_resolved"]]
  if (is.na(dp %||% NA) && is.na(dg %||% NA)) { rec_row("A", "NOT RENDERED", "A1-A5", "no folder"); return(invisible(NULL)) }
  gp <- .gate(dp); gg <- .gate(dg)
  a1 <- .tri({ x <- .gate_all(gp); y <- .gate_all(gg); if (identical(x, FALSE) || identical(y, FALSE)) FALSE else if (is.na(x) || is.na(y)) NA else TRUE })
  V1row("A", "A1 every fit on both tracks passes the convergence gate", sprintf("pooled %s; gear %s", .gate_str(gp), .gate_str(gg)),
        "8 of 8 BSS", a1, "A component on its PE is not an authoritative BSS estimate of that component.")
  dv <- if (!is.null(gp)) .num1(gp$divergence_fraction[gp$fit == FIT_SHORE_AG]) else NA_real_
  a2 <- .tri(if (is.na(dv)) NA else dv < 0.01)
  V1row("A", "A2 the pooled shore all-gear fit is under 1% divergent", fmt(100 * dv, 2), "< 1%", a2,
        "The A31 container refits: 0.05% to 0.19%. The old model: 1.6% to 23.3% (D33).")
  st <- lapply(FITS, function(f) .stuck(dp, f)); names(st) <- FITS
  a3 <- .tri(.stuck_ok(st))
  V1row("A", "A3 no pooled chain is stuck (step size and treedepth saturation)",
        paste(sprintf("%s: %s", sub("_Dungeness_Kept", "", names(st)), vapply(st, `[[`, "", "why")), collapse = " | "),
        "every chain's step >= 0.1 x its fit's median; saturation < 5%", a3,
        "B57's stuck chain: step 2.6e-5 against 0.006-0.014, 25% saturation, and a gate failure.")
  ref <- .here("05_output", REF_AUTH_POOLED)
  moves <- vapply(FITS, function(f) .sd_units(.comp(dp, f), .comp(ref, f)), numeric(1))
  a4 <- if (all(is.finite(moves)) && all(abs(moves) <= 0.2)) "PASS" else "REVIEW"
  V1row("A", "A4 each component against the 2026-09-28 authoritative run, in this render's posterior SD",
        paste(sprintf("%s %s -> %s (%s SD)", sub("_Dungeness_Kept", "", FITS), fmt(vapply(FITS, function(f) .comp(ref, f)$median, 0), 0),
                      fmt(vapply(FITS, function(f) .comp(dp, f)$median, 0), 0), fmt(moves, 2)), collapse = "; "),
        "|move| <= 0.2 SD", a4,
        "Never FAIL: A31 changes the model as well as the draws, so a move beyond 0.2 SD is recorded and explained, not rejected.")
  cc <- unique(unlist(lapply(unique(stats::na.omit(c(dp, dg))), function(d) Sys.glob(file.path(dirname(d), "cross_check_*.csv")))))
  ccx <- NULL
  for (f in cc[order(file.mtime(cc), decreasing = TRUE)]) { x <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(x) && identical(x$pooled_folder[1], basename(dp %||% "")) && identical(x$gear_folder[1], basename(dg %||% ""))) { ccx <- x; break } }
  a5 <- .tri(if (is.null(ccx)) NA else identical(ccx$verdict[1], "PASS"))
  V1row("A", "A5 the gear cross-check is within the configured tolerance",
        if (is.null(ccx)) "no cross_check csv names these two folders" else sprintf("pooled %s, gear %s: %s%% (%s)", fmt(ccx$pooled_median, 0), fmt(ccx$gear_median, 0), fmt(ccx$gear_minus_pooled_pct, 2), ccx$verdict),
        sprintf("|gap| <= %s%%", fmt(100 * (BASE$cross_check_tolerance %||% 0.02), 0)), a5,
        "The orchestrator's own verdict, read back (or the runner's, same formula, when A was resumed).")
  trees <- vapply(c(pooled = dp, gear_resolved = dg), function(d) if (is.na(d %||% NA)) NA_character_ else .manifest_tree(d), character(1))
  a6 <- if (all(!is.na(trees)) && all(grepl("^clean", trees))) "PASS" else "REVIEW"
  V1row("A", "A6 NOT BLOCKING: the manifest records a clean tree", paste(sprintf("%s: %s", names(trees), trees), collapse = " | "),
        "clean", a6, "A render from a modified tree can still be the method of record, but the difference must be named in the box.")
  vs <- c(A1 = a1, A2 = a2, A3 = a3, A5 = a5); res <- .worst(vs)
  pt <- .port(dp)
  rec_row("A", switch(res, PASS = "ELIGIBLE to become the authoritative run",
                      FAIL = "NOT ELIGIBLE: the 2026-09-28 run stays authoritative",
                      REVIEW = "NOT YET DECIDABLE: a clause could not be read; the 2026-09-28 run stays authoritative until it is"),
          "A1, A2, A3 and A5 PASS (A4 read for what moved; A6 noted)",
          sprintf("%s; port %s [%s, %s]; folders %s, %s%s", paste(sprintf("%s %s", names(vs), vs), collapse = ", "),
                  fmt(.num1(pt$BSS_median), 0), fmt(.num1(pt$BSS_lo95), 0), fmt(.num1(pt$BSS_hi95), 0),
                  basename(dp %||% "NA"), basename(dg %||% "NA"), if (a6 == "PASS") "" else "; A6: the tree was not clean or not recorded"))
}

verdict_R2 <- function(dir_r2, dir_a) {
  if (is.na(dir_r2 %||% NA)) { rec_row("R2", "NOT RENDERED", "R2-1..4", "no folder"); return(invisible(NULL)) }
  g <- .gate(dir_r2)
  r1 <- .tri(.gate_all(g))
  V1row("R2", "R2-1 every pooled fit at init_r = 2 passes the gate", .gate_str(g), "4 of 4 BSS", r1, "")
  dv <- if (!is.null(g)) .num1(g$divergence_fraction[g$fit == FIT_SHORE_AG]) else NA_real_
  r2 <- .tri(if (is.na(dv)) NA else dv < 0.01)
  V1row("R2", "R2-2 shore all-gear under 1% divergent at init_r = 2", fmt(100 * dv, 2), "< 1%", r2,
        "At init_r = 2 the old model left a chain stuck here (23.3%, D33).")
  st <- lapply(FITS, function(f) .stuck(dir_r2, f)); r3 <- .tri(.stuck_ok(st))
  V1row("R2", "R2-3 no stuck chain at init_r = 2", paste(vapply(st, `[[`, "", "why"), collapse = " | "), "as A3", r3, "")
  mv <- vapply(FITS, function(f) .sd_units(.comp(dir_r2, f), .comp(dir_a, f)), numeric(1))
  r4 <- .tri(if (!all(is.finite(mv))) NA else all(abs(mv) <= 0.2))
  V1row("R2", "R2-4 each component within 0.2 posterior SD of A", paste(sprintf("%s %s SD", sub("_Dungeness_Kept", "", FITS), fmt(mv, 2)), collapse = "; "),
        "|move| <= 0.2 SD", r4,
        "Same model and seed, another starting radius. With n_eff >= 400 the Monte Carlo SE of a difference is at most about 0.09 SD, so a spurious FAIL is possible; it errs toward keeping 0.5.")
  vs <- c(`R2-1` = r1, `R2-2` = r2, `R2-3` = r3, `R2-4` = r4); res <- .worst(vs)
  rec_row("R2", switch(res, PASS = "the A31 model does not depend on the radius: init_r = 0.5 is a precaution, not a crutch (keep or revert; the method is the model)",
                       FAIL = "keep init_r = 0.5: at the default radius the A31 model failed a clause",
                       REVIEW = "keep init_r = 0.5 (undecided: a clause could not be read)"),
          "R2-1..R2-4 all PASS",
          sprintf("%s; shore all-gear %s%% divergent; moves %s SD", paste(sprintf("%s %s", names(vs), vs), collapse = ", "), fmt(100 * dv, 2), paste(fmt(mv, 2), collapse = "/")))
}

verdict_D3 <- function(dir_d3, dir_a_pooled, dir_a_gear) {
  if (is.na(dir_d3 %||% NA)) { rec_row("D3", "NOT RENDERED", "D3-1..3", "no folder"); return(invisible(NULL)) }
  g <- .gate(dir_d3)
  c1 <- .tri(.gate_all(g))
  V1row("D3", "D3-1 every fit passes the gate", .gate_str(g), "4 of 4 BSS", c1, "")
  aq <- .adq(dir_d3, FIT_SHORE_AG); n_c <- .nobs(dir_d3, FIT_SHORE_AG, "catch")
  pl <- .num1(aq$p_loo_frac); nb <- .num1(aq$n_pareto_bad)
  c2 <- .tri(if (is.null(aq) || !is.finite(pl) || !is.finite(nb) || !is.finite(n_c) || n_c <= 0) NA else pl <= 0.15 && nb <= 0.05 * n_c)
  V1row("D3", "D3-2 the gear shore all-gear fit at weekly is adequate",
        if (is.null(aq)) "model_adequacy.csv missing" else sprintf("p_loo_frac %s (%s); bad k %s of %s catch obs", fmt(pl, 4), aq$p_loo_worst_stream, fmt(nb, 0), fmt(n_c, 0)),
        "p_loo_frac <= 0.15; bad k <= 5% of catch n_obs", c2,
        "The 2026-09-14 gear ladder: 0.0938 and 0 bad k at weekly; the pooled daily fit that was rejected: 0.352, 41.")
  n <- function(x) .bss_normalize_resolution(x %||% NA)
  pc <- BASE$ar_max_resolution$pooled$shore
  want <- c(shore_ring_net_only_Dungeness_Kept = n(pc$pot_closure), shore_all_gear_Dungeness_Kept = n(pc$all_gear),
            private_boat_ring_net_only_Dungeness_Kept = n(BASE$gear_period_bss$pot_closure),
            private_boat_all_gear_Dungeness_Kept = n(BASE$gear_period_bss$all_gear))
  got <- if (is.null(g)) rep(NA_character_, 4) else vapply(names(want), function(f) { r <- g$ar_resolution[g$fit == f]; if (length(r)) n(r[1]) else NA_character_ }, character(1))
  c3 <- .tri(if (any(is.na(got))) NA else identical(unname(got), unname(want)))
  V1row("D3", "D3-3 every fit ran at the period asked for", paste(sprintf("%s %s", sub("_Dungeness_Kept", "", names(want)), got), collapse = "; "),
        paste(sprintf("%s %s", sub("_Dungeness_Kept", "", names(want)), want), collapse = "; "), c3,
        "A fit that silently used another period would make the comparison meaningless (the 2026-09-14 lesson).")
  pd3 <- .num1(.port(dir_d3)$BSS_median); pa <- .num1(.port(dir_a_pooled)$BSS_median); pag <- .num1(.port(dir_a_gear)$BSS_median)
  V1row("D3", "D3-4 REPORTED, NOT A CRITERION: the cross-track gap at the matched periods",
        sprintf("gear matched %s vs pooled %s (%s%%); the shipped gear track %s (%s%%)", fmt(pd3, 0), fmt(pa, 0), fmt(100 * (pd3 - pa) / pa, 2), fmt(pag, 0), fmt(100 * (pag - pa) / pa, 2)),
        "no threshold", "INFO", "Choosing the gear period to shrink the gap would tune one estimate to another (rule 5, 2026-09-13).")
  for (f in FITS) {
    a <- .comp(dir_a_gear, f); b <- .comp(dir_d3, f)
    V1row("D3", sprintf("D3-4 REPORTED: %s, shipped gear track (A) -> D3%s", sub("_Dungeness_Kept", "", f),
                        if (f == FIT_SHORE_AG) " (the one fit whose period moves)" else " (period unchanged: a control)"),
          sprintf("%s -> %s (%s SD)", fmt(a$median, 0), fmt(b$median, 0), fmt(.sd_units(b, a), 2)), "no threshold", "INFO", "")
  }
  vs <- c(`D3-1` = c1, `D3-2` = c2, `D3-3` = c3); res <- .worst(vs)
  rec_row("D3", switch(res, PASS = sprintf("ADOPT: gear_period_bss per population (shore all-gear %s, shore pot closure %s; boat all-gear %s, boat pot closure %s, as shipped); D3 closes and the shore half of the cross-check becomes like-for-like",
                                           want[[2]], want[[1]], want[[4]], want[[3]]),
                       FAIL = "DO NOT ADOPT: a clause failed", REVIEW = "UNDECIDED: a clause could not be read"),
          "D3-1, D3-2, D3-3 PASS",
          sprintf("%s; port %s; shore all-gear p_loo_frac %s", paste(sprintf("%s %s", names(vs), vs), collapse = ", "), fmt(pd3, 0), fmt(pl, 4)))
}

verdict_D6 <- function(dir_d6, dir_d3) {
  if (is.na(dir_d6 %||% NA) || is.na(dir_d3 %||% NA)) { rec_row("D6", "NOT RENDERED (needs D3 and D6)", "D6-1..4", "no folder"); return(invisible(NULL)) }
  g <- .gate(dir_d6)
  c1 <- .tri(.gate_all(g))
  V1row("D6", "D6-1 every fit passes the gate", .gate_str(g), "4 of 4 BSS", c1, "")
  pr <- function(fit) loo_elpd_paired(file.path(dir_d3, sprintf("loo_pointwise_catch_%s.csv", fit)),
                                      file.path(dir_d6, sprintf("loo_pointwise_catch_%s.csv", fit)), label = sprintf("D3 -> D6 %s", fit))
  el <- tryCatch(pr(FIT_SHORE_AG), error = function(e) NULL)
  ratio <- if (is.null(el)) NA_real_ else .num1(el$ratio)
  c2 <- .tri(if (!is.finite(ratio)) NA else ratio > 2)
  V1row("D6", "D6-2 the paired catch elpd gain (shore all-gear, D6 minus D3) exceeds 2 paired SE",
        if (is.null(el)) "NOT COMPUTABLE" else loo_elpd_paired_str(el), "> 2 paired SE", c2,
        "The 2026-09-14 rung G5: +11.3 nats at 2.29 SE; the pooled adoption +11.6 at 2.30.")
  if (!is.null(el)) V1row("D6", "REPORTED: the elpd difference by observed count", loo_elpd_by_count_str(el), "no threshold", "INFO",
                          "Where the loss sits matters more than its sign.")
  z0a <- .bin(dir_d3, FIT_SHORE_AG, 0L); z0b <- .bin(dir_d6, FIT_SHORE_AG, 0L); z1a <- .bin(dir_d3, FIT_SHORE_AG, 1L); z1b <- .bin(dir_d6, FIT_SHORE_AG, 1L)
  imp0 <- abs(z0b) < abs(z0a); imp1 <- abs(z1b) < abs(z1a)
  c3 <- if (any(!is.finite(c(z0a, z0b, z1a, z1b)))) "REVIEW" else if (imp0 && imp1) "PASS" else if (imp0 || imp1) "REVIEW" else "FAIL"
  V1row("D6", "D6-3 both count bins improve", sprintf("zero z %s -> %s; one z %s -> %s", fmt(z0a, 2), fmt(z0b, 2), fmt(z1a, 2), fmt(z1b, 2)),
        "|z| falls in BOTH bins (one bin only: REVIEW; neither: FAIL)", c3, "The 2026-09-14 rung: zero 3.69 -> 2.03, one -6.13 -> -3.33.")
  zz <- el$zeros; pp <- el$positives
  c4 <- .tri(if (is.null(zz) || is.null(pp) || !is.finite(zz[["diff"]]) || !is.finite(pp[["diff"]])) NA
             else (pp[["diff"]] >= 0 || abs(pp[["diff"]]) < zz[["diff"]]))
  V1row("D6", "D6-4 the gain is not bought entirely at the zeros",
        if (c4 == "REVIEW") "NOT COMPUTABLE" else sprintf("zeros %s nats (n %s), positives %s (n %s)", fmt(zz[["diff"]], 1), fmt(zz[["n"]], 0), fmt(pp[["diff"]], 1), fmt(pp[["n"]], 0)),
        "positive-count loss < zero-count gain", c4, "")
  th <- .theta(dir_d6, FIT_SHORE_AG)
  V1row("D6", "D6-5 REPORTED: theta_C (shore all-gear) against its Beta(1, 9) prior mean 0.10", fmt(th, 4), "no threshold", "INFO",
        "2026-09-14: 0.170.")
  elp <- tryCatch(pr(FIT_SHORE_PC), error = function(e) NULL)
  V1row("D6", "D6-5 REPORTED: the shore pot-closure fit's pair", if (is.null(elp)) "NOT COMPUTABLE" else loo_elpd_paired_str(elp), "no threshold", "INFO", "")
  pa <- .num1(.port(dir_d3)$BSS_median); pb <- .num1(.port(dir_d6)$BSS_median)
  V1row("D6", "D6-5 REPORTED, NOT A CRITERION: the port total", sprintf("%s -> %s (%s%%)", fmt(pa, 0), fmt(pb, 0), fmt(100 * (pb - pa) / pa, 2)),
        "no threshold", "INFO", "The total is scaled by (1 - theta_C) by design.")
  vs <- c(`D6-1` = c1, `D6-2` = c2, `D6-3` = c3, `D6-4` = c4); res <- .worst(vs)
  rec_row("D6", switch(res, PASS = "ADOPT: catch_zi_tracks = c(\"pooled\", \"gear_resolved\"); D6 closes",
                       FAIL = "DO NOT ADOPT: a clause failed",
                       REVIEW = "UNDECIDED: a clause could not be computed or only one count bin improved (blocks adoption without counting against it)"),
          "D6-1 to D6-4 all PASS; any FAIL -> DO NOT ADOPT; otherwise any REVIEW -> UNDECIDED",
          sprintf("%s; %s", paste(sprintf("%s %s", names(vs), vs), collapse = ", "), if (is.null(el)) "no paired elpd" else loo_elpd_paired_str(el)))
}

D29 <- list()
d29_rows <- function(track, dir, stage) {
  lg <- .arlog(dir); if (is.null(lg)) return(invisible(NULL))
  lg <- lg[lg$fit == FIT_BOAT_AG, , drop = FALSE]; if (!nrow(lg)) return(invisible(NULL))
  n_all <- .nobs(dir, FIT_BOAT_AG, c("trailer", "osp", "catch", "gear"))
  for (i in seq_len(nrow(lg))) {
    res <- .bss_normalize_resolution(lg$ar_resolution[i])
    D29[[length(D29) + 1]] <<- data.frame(track = track, stage = stage, rung = res, P_n = lg$P_n[i],
      ladder = as.logical(lg$escalation_enabled[i]), selected = as.logical(lg$selected[i]),
      pass_gate = as.logical(lg$pass_convergence[i]), divergence_fraction = lg$divergence_fraction[i],
      catch_median = round(lg$catch_median[i]), catch_lo95 = round(lg$catch_lo95[i]), catch_hi95 = round(lg$catch_hi95[i]),
      p_loo_frac = lg$p_loo_frac[i], n_pareto_bad = lg$n_pareto_bad[i], n_obs_loo = n_all,
      cov50_trailer = if ("cov50_trailer" %in% names(lg)) lg$cov50_trailer[i] else NA_real_, cov50_catch = lg$cov50_catch[i],
      block_file = file.path(dir, sprintf("ladder_block_joint_%s_%s.csv", FIT_BOAT_AG, res)), stringsAsFactors = FALSE)
  }
}
d29_daily_row <- function(dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  g <- .gate(dir); if (is.null(g) || !any(g$fit == FIT_BOAT_AG)) return(invisible(NULL))
  aq <- .adq(dir, FIT_BOAT_AG); c1 <- .comp(dir, FIT_BOAT_AG)
  rung <- .bss_normalize_resolution(g$ar_resolution[g$fit == FIT_BOAT_AG][1])
  # ar_force must have reached the fit; a D29D fit at any other period is not the daily rung
  V1row("D29D", "the boat all-gear fit ran at daily (ar_force)", rung %||% "NA", "daily", .tri(if (is.na(rung)) NA else rung == "daily"),
        "A forced period that did not take would be read as the wrong rung.")
  if (!identical(rung, "daily")) return(invisible(NULL))
  D29[[length(D29) + 1]] <<- data.frame(track = "pooled", stage = "D29D", rung = rung,
    P_n = NA_real_, ladder = FALSE, selected = TRUE,
    pass_gate = as.logical(g$pass_convergence[g$fit == FIT_BOAT_AG][1]), divergence_fraction = g$divergence_fraction[g$fit == FIT_BOAT_AG][1],
    catch_median = round(c1$median), catch_lo95 = round(c1$lo), catch_hi95 = round(c1$hi),
    p_loo_frac = .num1(aq$p_loo_frac), n_pareto_bad = .num1(aq$n_pareto_bad), n_obs_loo = .nobs(dir, FIT_BOAT_AG, c("trailer", "osp", "catch", "gear")),
    cov50_trailer = NA_real_, cov50_catch = NA_real_,
    block_file = file.path(dir, sprintf("loo_block_joint_%s.csv", FIT_BOAT_AG)), stringsAsFactors = FALSE)
}
verdict_D29 <- function(dirs) {
  if (!length(D29)) { rec_row("D29", "NOT RENDERED", "D29-0..4", "no ladder folder"); return(invisible(NULL)) }
  R <- do.call(rbind, D29)
  # D29-2, three-valued: NA (not computed) is UNKNOWN, not inadequate
  R$adequate <- ifelse(!is.finite(R$p_loo_frac) | !is.finite(R$n_pareto_bad) | !is.finite(R$n_obs_loo) | R$n_obs_loo <= 0, NA,
                       R$p_loo_frac <= 0.15 & R$n_pareto_bad <= 0.05 * R$n_obs_loo)
  R$block_vs_monthly <- NA_character_; R$block_diff <- NA_real_; R$block_se <- NA_real_; R$block_weeks <- NA_integer_
  R$block_class <- ifelse(R$rung == "monthly", "reference", NA_character_)
  for (tr in unique(R$track)) {
    mref <- R[R$track == tr & R$rung == "monthly" & R$stage != "D29D", , drop = FALSE]
    # the partner of every finer rung, the daily stage's included, is the LADDER's monthly rung
    # (same fit, same seed); A's kept monthly fit stands in only if the ladder has none
    tab_m <- if (nrow(mref)) .block(dirname(mref$block_file[1]), basename(mref$block_file[1])) else NULL
    if (is.null(tab_m) && tr == "pooled" && !is.na(dirs$A[["pooled"]] %||% NA))
      tab_m <- .block(dirs$A[["pooled"]], sprintf("loo_block_joint_%s.csv", FIT_BOAT_AG))
    for (i in which(R$track == tr & R$rung != "monthly")) {
      tab_r <- .block(dirname(R$block_file[i]), basename(R$block_file[i]))
      if (is.null(tab_m) || is.null(tab_r)) { R$block_class[i] <- "NOT EVALUABLE (block file missing)"; next }
      cmp <- tryCatch(bss_block_cv_compare(tab_m, tab_r), error = function(e) NULL)
      if (is.null(cmp)) { R$block_class[i] <- "NOT EVALUABLE (compare failed)"; next }
      # bss_block_cv_compare(A, B) is B against A: diff = elpd(rung) - elpd(monthly), over the
      # weeks reliable in both; reliable_share = those weeks over the weeks both tables hold.
      R$block_diff[i] <- .num1(cmp$diff); R$block_se[i] <- .num1(cmp$se); R$block_weeks[i] <- as.integer(.num1(cmp$n_used))
      R$block_vs_monthly[i] <- tryCatch(bss_block_cv_str(cmp, "monthly", R$rung[i]), error = function(e) NA_character_)
      z <- .num1(cmp$ratio)
      R$block_class[i] <- if (isTRUE(cmp$identical)) "NO EVIDENCE (identical tables)"
                          else if (!isTRUE(.num1(cmp$reliable_share) >= 0.5)) "NOT EVALUABLE (under half the weeks reliable in both)"
                          else if (!is.finite(z)) "NOT EVALUABLE (no paired SE)"
                          else if (z > 2) "BETTER" else if (z < -2) "WORSE" else "NO EVIDENCE"
    }
  }
  # each rung's standing for D29-4: out (ineligible or inadequate), unknown, or its block class
  R$standing <- ifelse(R$pass_gate %in% FALSE, "out: fails the gate",
                ifelse(R$adequate %in% FALSE, "out: inadequate",
                ifelse(is.na(R$pass_gate) | is.na(R$adequate), "UNKNOWN (gate or adequacy not computed)",
                ifelse(R$rung == "monthly", "reference", R$block_class))))
  D29_TABLE <<- R
  for (i in seq_len(nrow(R)))
    V1row(R$stage[i], sprintf("D29 %s %s: gate / adequacy / interpolation vs monthly", R$track[i], R$rung[i]),
          sprintf("gate %s; div %s%%; catch %s [%s, %s]; p_loo_frac %s; bad k %s of %s; %s",
                  R$pass_gate[i], fmt(100 * R$divergence_fraction[i], 2), fmt(R$catch_median[i], 0), fmt(R$catch_lo95[i], 0), fmt(R$catch_hi95[i], 0),
                  fmt(R$p_loo_frac[i], 4), fmt(R$n_pareto_bad[i], 0), fmt(R$n_obs_loo[i], 0),
                  if (R$rung[i] == "monthly") "block: (reference)" else sprintf("block %s: %s", R$block_class[i], R$block_vs_monthly[i] %||% "")),
          "D29-1 pass; D29-2 p_loo_frac <= 0.15 and bad k <= 5% of the LOO observations; D29-3 > +2 paired SE",
          if (R$pass_gate[i] %in% FALSE || R$adequate[i] %in% FALSE) "FAIL" else if (is.na(R$pass_gate[i]) || is.na(R$adequate[i])) "REVIEW" else "INFO",
          "Per rung, per track. FAIL here means the rung is out of the running, not that the batch failed. The decision (D29-4) reads the pooled rows only.")
  # D29-6: the ladder's monthly rung against A's monthly boat all-gear fit, the control (pooled
  # only: A's gear boat fit is on the shipped gear configuration and D29G's base is D3)
  m <- R[R$track == "pooled" & R$rung == "monthly" & R$stage == "D29P", , drop = FALSE]; a <- dirs$A[["pooled"]]
  if (nrow(m) && !is.na(a %||% NA)) { ca <- .comp(a, FIT_BOAT_AG)
    V1row("D29P", "D29-6 REPORTED: the ladder's monthly rung against A's monthly boat all-gear fit",
          sprintf("%s vs %s (%s SD)", fmt(m$catch_median[1], 0), fmt(ca$median, 0), fmt((m$catch_median[1] - ca$median) / ca$sd, 2)),
          "near 0", "INFO", "The control that the ladder changes nothing but the rung.") }
  # D29-0: the pooled ladder is complete
  P <- R[R$track == "pooled", , drop = FALSE]; PL <- P[P$stage == "D29P", , drop = FALSE]
  complete <- nrow(PL) == 3 && setequal(PL$rung, c("weekly", "biweekly", "monthly")) && all(PL$ladder %in% TRUE)
  V1row("D29P", "D29-0 the pooled ladder is complete (weekly, biweekly, monthly, ladder engaged)",
        if (nrow(PL)) paste(sprintf("%s (ladder %s)", PL$rung, PL$ladder), collapse = "; ") else "no ladder rows",
        "exactly weekly, biweekly, monthly, each with escalation_enabled TRUE", .tri(complete), "An incomplete ladder cannot be read by D29-4.")
  rank <- c(daily = 4, weekly = 3, biweekly = 2, monthly = 1)
  mon <- P[P$rung == "monthly" & P$stage == "D29P", , drop = FALSE]
  finer <- P[P$rung != "monthly", , drop = FALSE]
  finer <- finer[order(-rank[finer$rung]), , drop = FALSE]
  decided <- function(s) grepl("^(out:|BETTER|WORSE|NO EVIDENCE)", s)
  ok_rows <- P$pass_gate %in% TRUE & P$adequate %in% TRUE
  span <- if (any(ok_rows)) range(P$catch_median[ok_rows]) else c(NA, NA)
  standing <- paste(sprintf("%s %s", P$rung, P$standing), collapse = "; ")
  if (!("daily" %in% P$rung)) standing <- paste0(standing, "; daily not rendered")
  better <- finer[finer$standing %in% "BETTER", , drop = FALSE]
  # a D29D that did not render (failed, timed out, or dropped from STAGES) is not a rung the
  # rule decided on; the recommendation says so rather than implying daily was weighed
  scope <- if ("daily" %in% P$rung) "" else " [scope: weekly to monthly; the daily rung was not rendered]"
  if (!complete) {
    rec_row("D29", "REVIEW: the pooled ladder is incomplete, so the rule does not run", "D29-0 / D29-4(c)", standing)
  } else if (nrow(better)) {
    b <- better[1, ]
    undecided_finer <- finer[rank[finer$rung] > rank[b$rung] & !decided(finer$standing), , drop = FALSE]
    if (nrow(undecided_finer))
      rec_row("D29", sprintf("REVIEW: %s interpolates BETTER than monthly, but the finer %s could not be decided", b$rung, paste(undecided_finer$rung, collapse = ", ")),
              "D29-4(c)", standing)
    else
      rec_row("D29", sprintf("MOVE the pooled boat all-gear cap to %s (ar_max_resolution$pooled$private_boat$all_gear), and the gear track to match (D3)%s", b$rung, scope),
              "D29-4(a): the finest eligible, adequate rung that interpolates BETTER than monthly",
              sprintf("%s: catch %s against monthly %s; block %s; adequate rungs span %s to %s. %s", b$rung, fmt(b$catch_median, 0), fmt(mon$catch_median[1], 0),
                      b$block_vs_monthly, fmt(span[1], 0), fmt(span[2], 0), standing))
  } else if (nrow(mon) && isTRUE(mon$pass_gate[1]) && isTRUE(mon$adequate[1]) && all(decided(finer$standing))) {
    rec_row("D29", paste0("KEEP monthly: the data do not support moving the cap. The span across the adequate rungs is this component's resolution uncertainty; applying the shore precedent (finest adequate rung) is Matt's decision, not the rule's", scope),
            "D29-4(b)", sprintf("adequate rungs span %s to %s. %s", fmt(span[1], 0), fmt(span[2], 0), standing))
  } else {
    rec_row("D29", "REVIEW: the rule could not decide", "D29-4(c)",
            sprintf("monthly gate %s, adequate %s. %s", mon$pass_gate[1] %||% NA, mon$adequate[1] %||% NA, standing))
  }
  G <- R[R$track == "gear_resolved", , drop = FALSE]
  if (nrow(G)) V1row("D29G", "D29-5 REPORTED, NOT A CRITERION: the gear ladder",
                     paste(sprintf("%s %s [%s, %s] (%s)", G$rung, fmt(G$catch_median, 0), fmt(G$catch_lo95, 0), fmt(G$catch_hi95, 0), G$standing), collapse = "; "),
                     "no threshold", "INFO", "Corroboration or contradiction of the pooled ladder; never the decider (circular).")
}
D29_TABLE <- NULL

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
banner(sprintf("AUTHORITATIVE BATCH (A31) + D29 / D3 / D6  |  DRY_RUN = %s  |  code %s  |  git %s", DRY_RUN, CODE_FP, GIT_SHA %||% "?"))
if (isTRUE(DRY_RUN))
  cat("\n  DRY RUN: S0 runs and every stage's plan is printed; nothing is fitted or written.\n",
      " Start with:  Rscript 06_diagnostics/run_authoritative_batch_2026-09-29.R --go\n", sep = "")
if ("S0" %in% STAGES) tryCatch(stage_S0(), error = function(e) V1row("S0", "the desk stage ran", conditionMessage(e), "no error", "FAIL", ""))
# A desk FAIL blocks only the stages it concerns (DEPENDS); a FAIL labelled S0 (the method of
# record, pandoc, the inputs, the NWS archive, the desk stage itself) blocks every stage.
.desk_fail <- unique(vapply(Filter(function(r) identical(r$verdict, "FAIL"), V), function(r) r$stage, character(1)))
BLOCKED <- if ("S0" %in% .desk_fail) setdiff(STAGES, "S0") else intersect(STAGES, unique(unlist(DEPENDS[intersect(.desk_fail, names(DEPENDS))])))
if (length(BLOCKED)) {
  cat(sprintf("\n  *** Desk check(s) FAILED for %s: %s %s NOT be fitted. Read the verdicts below. ***\n",
              paste(.desk_fail, collapse = ", "), paste(BLOCKED, collapse = ", "), if (isTRUE(DRY_RUN)) "would" else "will"))
  for (sid in BLOCKED) STAGE_STATUS[[sid]] <- list(status = "BLOCKED by a desk check", dirs = NULL, minutes = NA_real_)
}
DIRS <- list()
for (sid in setdiff(STAGES, c("S0", BLOCKED))) {
  DIRS[[sid]] <- tryCatch(run_stage(sid), error = function(e) {
    V1row(sid, "the stage ran", conditionMessage(e), "no error", "FAIL", "The batch continued with the next stage.")
    stats::setNames(rep(NA_character_, length(PARTS(sid))), PARTS(sid)) })
}
# Evaluate EVERY stage whose folder exists at the current digest, whether or not it ran now.
for (sid in setdiff(names(STAGE_DEFS), "S0"))
  if (is.null(DIRS[[sid]]) || all(is.na(DIRS[[sid]])))
    DIRS[[sid]] <- vapply(PARTS(sid), function(pt) find_stage_dir(sid, pt), character(1))
banner("EVALUATION")
for (sid in setdiff(names(STAGE_DEFS), "S0")) for (pt in PARTS(sid)) record_components(sid, DIRS[[sid]][[pt]], pt)
tryCatch(verdict_A(DIRS$A), error = function(e) V1row("A", "verdicts computed", conditionMessage(e), "no error", "REVIEW", ""))
tryCatch(verdict_R2(DIRS$R2[["pooled"]], DIRS$A[["pooled"]]), error = function(e) V1row("R2", "verdicts computed", conditionMessage(e), "no error", "REVIEW", ""))
tryCatch(verdict_D3(DIRS$D3[["gear_resolved"]], DIRS$A[["pooled"]], DIRS$A[["gear_resolved"]]), error = function(e) V1row("D3", "verdicts computed", conditionMessage(e), "no error", "REVIEW", ""))
tryCatch(verdict_D6(DIRS$D6[["gear_resolved"]], DIRS$D3[["gear_resolved"]]), error = function(e) V1row("D6", "verdicts computed", conditionMessage(e), "no error", "REVIEW", ""))
tryCatch({ d29_rows("pooled", DIRS$D29P[["pooled"]], "D29P"); d29_rows("gear_resolved", DIRS$D29G[["gear_resolved"]], "D29G")
           d29_daily_row(DIRS$D29D[["pooled"]]); verdict_D29(DIRS) },
         error = function(e) V1row("D29", "verdicts computed", conditionMessage(e), "no error", "REVIEW", ""))

VV <- if (length(V)) do.call(rbind, V) else NULL
RR <- if (length(REC)) do.call(rbind, REC) else NULL
CC <- if (length(COMP)) do.call(rbind, COMP) else NULL
SS <- do.call(rbind, lapply(setdiff(names(STAGE_DEFS), "S0"), function(sid) {
  s <- STAGE_STATUS[[sid]]; d <- DIRS[[sid]]
  data.frame(stage = sid, item = STAGE_DEFS[[sid]]$item, status = if (!is.null(s)) s$status else if (all(!is.na(d))) "present" else "absent",
             folders = paste(ifelse(is.na(d), "-", sub(paste0("^", .here("05_output"), "/?"), "", d)), collapse = "; "),
             minutes = round(s$minutes %||% NA_real_, 1), digest = stage_digest(sid), code = CODE_FP, git = GIT_SHA %||% NA,
             stringsAsFactors = FALSE) }))
if (!is.null(VV)) {
  cat("\n"); for (i in seq_len(nrow(VV)))
    cat(sprintf("  [%-6s] %-5s %s\n           observed: %s\n", VV$verdict[i], VV$stage[i], VV$criterion[i], VV$observed[i]))
}
banner("RECOMMENDATIONS")
if (!is.null(RR)) for (i in seq_len(nrow(RR))) cat(sprintf("  %-4s %s\n       rule: %s\n       evidence: %s\n\n", RR$item[i], RR$recommendation[i], RR$rule[i], RR$evidence[i]))
if (!isTRUE(DRY_RUN)) {
  w <- function(x, nm) if (!is.null(x)) utils::write.csv(x, .here("05_output", sprintf("%s_%s.csv", OUT_PREFIX, nm)), row.names = FALSE)
  w(SS, "stages"); w(CC, "components"); w(VV, "verdicts"); w(RR, "recommendation")
  if (!is.null(D29_TABLE)) w(D29_TABLE, "d29_rungs")
  cat(sprintf("\n  Written: 05_output/%s_{stages,components,verdicts,recommendation,d29_rungs}.csv\n", OUT_PREFIX))
} else {
  cat("\n  DRY_RUN: nothing written. Start with --go.\n")
}

