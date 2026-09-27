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
# run_marine_block_cv_2026-09-26.R
#
# THE TEST THE MARINE HAZARD LADDER COULD NOT RUN (CHANGE_REGISTER D31 / B39), computed
# AFTER THE FACT from the five committed rungs of run_marine_hazard_batch_2026-09-25.R,
# without a single refit.
#
# WHAT THE LADDER SAID (2026-09-26). The boat SCA term is identified at eight standard
# errors (B_open -1.16 [-1.45, -0.87], rate ratio 0.31), harmless to the catch fit, and
# gained +7.6 nats on the trailer stream at 1.46 PAIRED SE against the baseline: real,
# and short of the +2 SE the rule asked for. Pointwise LOO scores each sampled count with
# every other sampled count still in the fit. The covariate exists for the UNSAMPLED
# days, where the AR interpolates and the advisory flag is the only information that a
# day was not ordinary. This runner holds out one calendar WEEK of effort counts at a
# time, every effort stream of the fit together, so for those streams every day in it is
# unsampled, and scores the week's joint predictive density under the posterior without
# it (PSIS, bss_block_cv.R; what still anchors a held-out week is stated in R4).
#
# WHERE THE DRAWS COME FROM. Each rung saved ppc_draws_<fit>.rds (save_bss_ppc_draws();
# lambda_E_S, r_E, R_G, R_G_boat, L_out, r_OSP, kappa_OSP, plus the fit's stan_data).
# The effort log-likelihoods are rebuilt from those exactly as the Stan generated
# quantities form them, and the rebuild is CHECKED against the committed
# loo_pointwise_<stream>_<fit>.csv on every observation before anything is scored: a
# reconstruction that does not reproduce the run's own pointwise lpd to its rounding is
# refused, not averaged over. The .rds files are git-ignored, so this runs on the
# machine that rendered the rungs (or anywhere the folders were copied with them).
#
# THE DECISION RULE, STATED HERE BEFORE THE RUN.
#   R1. A fit is EVALUABLE only if its reconstructed pointwise lpd matches the committed
#       loo_pointwise file on every observation of every gear / trailer stream it has; a
#       gear or trailer stream WITHOUT a committed file is refused, not waved through. The
#       OSP stream had no committed pointwise file before B40 (2026-09-27): for those fits its
#       reconstructed posterior-mean count per observation is checked against the run's own
#       ppc_byobs_<fit>.csv fitted_mean (a 1,500-draw subsample, so to 10%), it is labelled
#       "approx-checked" in the output, and a fit whose OSP means miss that is refused too.
#       A fit rendered after B40 has loo_pointwise_osp_<fit>.csv and its OSP stream is
#       checked exactly, like the others.
#   R2. A comparison is EVALUABLE only if at least 70% of its weeks are RELIABLE (Pareto
#       k <= 0.7 in BOTH rungs) and at least 8 weeks remain. Otherwise the verdict is
#       REVIEW: "PSIS cannot carry this comparison", and the refit path is named. The
#       excluded weeks are part of the answer, not a footnote.
#   R3. Over the reliable weeks, the PAIRED held-out elpd difference (covariate rung minus
#       control): > +2 SE is PASS ("the covariate improves the interpolation on held-out
#       weeks"), < -2 SE is FAIL ("worse"), between is REVIEW ("no evidence either way").
#       Pairs: SCA in M2 against M1; bar alone in M3 against M1; bar beyond the archive in
#       M4 against M2; both in M4 against M1; auto (M5) against M1; and, once the ladder's M6
#       has rendered (B41, 2026-09-27), the season split against the constant term, M6
#       against M2, boat fits only (rule 9 of run_marine_hazard_batch_2026-09-25.R: the joint
#       boat row here is that rung's arbiter).
#   R4. THE LEAVE-OUT SET IS THE WHOLE WEEK'S EFFORT DATA. For the boat, a held-out week
#       removes its trailer counts AND its OSP counts together (35% of trailer days also
#       carry an OSP count), and each stream is then scored under those shared weights;
#       holding out one stream alone would leave the other anchoring lambda_E on the very
#       days being "predicted" (bss_block_cv.R header). Both boat streams are scored; the
#       OSP stream was never in the ladders' LOO (loo_summary_*.csv carries gear / trailer
#       / catch), so this is the first predictive read on the boat's second effort stream.
#       A term that helps one stream and hurts the other is a finding, not a pass. What
#       still anchors a held-out shore week: the four I/E observations (no log_lik_ie).
#   R5. The pot-closure fits are scored and reported; the adoption question is the
#       all-gear fits', where the harvest is.
#   R6. Nothing here is a total. The port is not read, not judged, not moved.
#
# WHAT THE RUN SAID (2026-09-26/27, c7e8cd5; read in VALIDATION_CAMPAIGN.md Section 1y).
#   R1 held on all 20 fits (gear / trailer to 5e-5; OSP means within 0.8 to 2.4% of the
#   1,500-draw ppc_byobs subsample). The boat all-gear fit was evaluable (34 of 37 trailer
#   weeks, 25 of 29 OSP weeks reliable); the shore fits were not (8 to 16 of 38 weeks at
#   weekly AR), as the header predicted. Boat SCA against the baseline: OSP +18.9 nats,
#   3.90 SE, PASS (21 of 25 weeks positive); trailer +7.9, 1.20 SE, REVIEW (22 of 34);
#   the two streams together +26.8 at 2.9 SE (the sum of the streams' weekly differences over
#   the 38 weeks reliable in every stream present, recomputed from the per-week tables; the
#   run did not carry the joint row). Bar beyond the archive: trailer -2.0 (-1.0 SE), OSP
#   -0.4: nothing. M5 == M4. Pot-closure boat trailer +4.1 at 2.66 SE on 8 of 11 weeks.
#   THE RULE AS PRE-COMMITTED WAS NOT MET (it asked each boat stream separately), AND IT WAS
#   UNDER-POWERED for the trailer stream: its paired SE over 34 weeks is 6.6 nats, so +2 SE
#   needed +13 nats, and the same effect scored 3.9 SE on the OSP counts and 1.6 SE on the
#   trailer counts over the same weeks (1.2 SE over all 34): the weaker instrument was asked to
#   clear the bar on its own. A C row.
#
# TWO AMENDMENTS, 2026-09-27, FOR EVERY RUN FROM HERE (the 2026-09-26 verdicts stand as
# written; this changes what the NEXT run is judged on, not what the last one said):
#   R3'. For a fit with more than one effort stream the JOINT row (the week's predictive
#        density of ALL its held-out effort observations together, bss_block_cv_fit()) is
#        the statistic R3 is applied to; the per-stream rows are reported beside it, and a
#        term that helps one stream and hurts the other is still a finding (R4).
#   R7.  Two rungs whose fits are the SAME model (a shore fit under a boat-only term) give
#        block tables equal to summation noise (1e-14 nats); that is reported INFO,
#        "identical fits", never PASS / FAIL / REVIEW. The 2026-09-26 run FAILED one such
#        row at "-2.31 SE" on a difference of 2e-14 nats (CHANGE_REGISTER C row).
#   Also written from here: marine_hazard_2026-09-26_blockcv_weeks.csv, the per-week
#   differences behind every pair, so the evidence can be sliced by season without
#   re-deriving it (the winter weeks are where the term moves the boat estimate and where
#   the trailer stream, the only winter stream, is uninformative: -0.6 +- 3.9 nats).
#
# WHAT IT WRITES
#   <rung folder>/loo_block_<stream>_<fit>.csv            per-week table, every rung
#                                                          (+ loo_block_joint_<fit>.csv for the boat)
#   05_output/marine_hazard_2026-09-26_blockcv.csv         per rung x fit x stream totals
#   05_output/marine_hazard_2026-09-26_blockcv_pairs.csv   the paired comparisons + verdicts
#   05_output/marine_hazard_2026-09-26_blockcv_weeks.csv   the per-week differences behind them
#   05_output/marine_hazard_2026-09-26_blockcv_verdicts.csv
#
# DRY_RUN <- TRUE ships: it inventories the rungs (folders, rds present, loo installed)
# and runs the reconstruction checks (R1) on the first fit that has its draws; it scores
# no block and writes nothing. Set FALSE and source again from the repository root to
# compute (minutes). Restore TRUE before committing: the harness asserts it.
###############################################################################

DRY_RUN <- TRUE
RUNGS   <- c(M1 = "MH-M1-off", M2 = "MH-M2-sca", M3 = "MH-M3-bar", M4 = "MH-M4-both", M5 = "MH-M5-auto",
             M6 = "MH-M6-split")   # M6 (B41, the season split) added 2026-09-27; absent until run_marine_hazard_batch renders it
K_MAX   <- 0.7          # PSIS reliability threshold
MIN_REL <- 0.70         # R2: reliable share needed for a comparison to be evaluable
MIN_WKS <- 8L           # R2: reliable weeks needed
OSP_TOL <- 0.10         # R1: relative tolerance of the OSP posterior-mean check against ppc_byobs (a 1,500-draw subsample)

# =========================================================================== #

.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions")))
  stop("Run this from the repository root (or from 06_diagnostics/): 03_R_functions not found.")
.here <- function(...) file.path(.root, ...)
setwd(.root)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
source(.here("03_R_functions", "bss_block_cv.R"))
source(.here("03_R_functions", "batch_verdict_helpers.R"))   # merge_csv_by
banner <- function(msg) cat("\n", strrep("=", 78), "\n ", msg, "\n", strrep("=", 78), "\n", sep = "")
PREFIX <- "pooled-CPUE-"
OUT_ROOT <- Sys.getenv("MH_BLOCKCV_ROOT", unset = .here("05_output"))   # tests point this at a fixture
FITS <- c(shore_all_gear = "shore_all_gear_Dungeness_Kept", boat_all_gear = "private_boat_all_gear_Dungeness_Kept",
          shore_pot_closure = "shore_ring_net_only_Dungeness_Kept", boat_pot_closure = "private_boat_ring_net_only_Dungeness_Kept")
STREAMS_OF <- list(shore_all_gear = c("gear"), shore_pot_closure = c("gear"),
                   boat_all_gear = c("trailer", "osp"), boat_pot_closure = c("trailer", "osp"))
# What gets a paired comparison: the reconstructed streams, plus the "joint" table
# bss_block_cv_fit() adds for a fit with more than one stream (R3').
SCORED_OF <- lapply(STREAMS_OF, function(s) if (length(s) > 1L) c(s, "joint") else s)
PAIRS <- list(list(b = "M2", a = "M1", what = "SCA (nws_sca_any) alone"),
              list(b = "M3", a = "M1", what = "bar restriction alone"),
              list(b = "M4", a = "M2", what = "bar restriction beyond the archive (M4 vs M2)"),
              list(b = "M4", a = "M1", what = "SCA + bar restriction"),
              list(b = "M5", a = "M1", what = "auto (what production would select)"),
              # B41: the season split against the constant term, boat fits only (M6 has no shore term)
              list(b = "M6", a = "M2", what = "SCA split by season beyond the constant term (M6 vs M2; rule 9 of the ladder)",
                   fits = c("boat_all_gear", "boat_pot_closure")))

find_outdir <- function(tag) {
  dirs <- list.dirs(OUT_ROOT, recursive = FALSE)
  hits <- unlist(lapply(dirs, function(d) { sub <- list.dirs(d, recursive = FALSE); sub[basename(sub) == paste0(PREFIX, tag)] }))
  if (!length(hits)) return(NA_character_)
  hits[order(basename(dirname(hits)), decreasing = TRUE)][1]
}
V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = observed, threshold = threshold,
                                    verdict = verdict, why = why, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Per rung x fit: the block tables for every effort stream, from the saved draws
# ---------------------------------------------------------------------------
# dates by day index: bss_daily_catch_<fit>.csv is one row per day of the sub-season, in order
.dates_by_day <- function(dir, fit) {
  f <- file.path(dir, sprintf("bss_daily_catch_%s.csv", fit))
  if (!file.exists(f)) return(NULL)
  as.Date(utils::read.csv(f, stringsAsFactors = FALSE)$event_date)
}
# The OSP stream's approximate check (R1): the reconstructed posterior-mean count per OSP
# observation against ppc_byobs_<fit>.csv's fitted_mean, which the run computed from a
# seeded 1,500-draw subsample of the same posterior (so agreement is to a few percent, and
# a formula error, forgetting L or R_G_boat, is a factor of two or more).
.osp_mean_check <- function(draws, sd, dir, fit, rec) {
  f <- file.path(dir, sprintf("ppc_byobs_%s.csv", fit))
  if (!file.exists(f)) return(list(ok = FALSE, skip = TRUE, why = "ppc_byobs_<fit>.csv absent; the OSP stream cannot be checked and is not scored"))
  pb <- utils::read.csv(f, stringsAsFactors = FALSE); pb <- pb[pb$data_type == "osp", , drop = FALSE]
  if (nrow(pb) != ncol(rec$ll)) return(list(ok = FALSE, why = sprintf("%d OSP rows in ppc_byobs, %d reconstructed", nrow(pb), ncol(rec$ll))))
  if (!identical(as.integer(pb$day_index), as.integer(rec$day))) return(list(ok = FALSE, why = "OSP day indices differ from ppc_byobs"))
  lam <- draws$lambda_E_S; idx <- seq_len(dim(lam)[1])
  mult <- if (!is.null(draws$R_G_boat)) 1 / as.numeric(draws$R_G_boat) else as.numeric(draws$R_T)
  mu <- vapply(seq_along(rec$day), function(i) {
    l <- rowSums(matrix(lam[idx, rec$section[i], rec$day[i], , drop = FALSE], nrow = length(idx))) * mult
    sc <- if (identical(as.integer(sd$osp_scale_is_tau %||% 0L), 1L)) as.matrix(draws$L_out)[idx, rec$day[i]] else as.numeric(draws$kappa_OSP)
    mean(l * sc) }, numeric(1))
  rel <- abs(mu - pb$fitted_mean) / pmax(pb$fitted_mean, 1e-6)
  list(ok = all(rel <= OSP_TOL), max_rel = max(rel), why = if (all(rel <= OSP_TOL)) "" else sprintf("max relative deviation %.1f%% > %.0f%%", 100 * max(rel), 100 * OSP_TOL))
}
# `score = FALSE` runs the checks only (the dry run). Returns list(ok, checks, tables, why).
block_tables <- function(dir, fit_key, write = TRUE, score = TRUE) {
  fit <- FITS[[fit_key]]
  rds <- file.path(dir, sprintf("ppc_draws_%s.rds", fit))
  if (!file.exists(rds)) return(list(ok = FALSE, why = sprintf("%s absent (git-ignored; present only where the rung rendered)", basename(rds)), checks = list()))
  obj <- readRDS(rds)
  draws <- obj$draws; sd <- obj$stan_data
  if (is.null(sd)) return(list(ok = FALSE, why = "the saved draws carry no stan_data", checks = list()))
  if (isTRUE(obj$thinned)) return(list(ok = FALSE, why = "the saved draws are THINNED; the pointwise check cannot pass and the block scores would not be the run's", checks = list()))
  dates <- .dates_by_day(dir, fit)
  if (is.null(dates)) return(list(ok = FALSE, why = "bss_daily_catch_<fit>.csv absent; no day-to-date map", checks = list()))
  if (length(dates) != (sd[["D"]] %||% length(dates))) return(list(ok = FALSE, why = sprintf("bss_daily_catch has %d rows; stan_data$D is %d", length(dates), sd[["D"]]), checks = list()))
  out <- list(ok = TRUE, tables = list(), checks = list(), why = "")
  recs <- list(); n_checked <- 0L
  for (sn in STREAMS_OF[[fit_key]]) {
    rec <- tryCatch(bss_effort_loglik(draws, sd, sn), error = function(e) e)
    if (inherits(rec, "error")) { out$ok <- FALSE; out$checks[[sn]] <- sprintf("%s: not reconstructed (%s)", sn, conditionMessage(rec)); break }
    if (!ncol(rec$ll)) { out$checks[[sn]] <- sprintf("%s: stream empty in this fit", sn); next }
    pw <- file.path(dir, sprintf("loo_pointwise_%s_%s.csv", sn, fit))
    # gear and trailer always have a committed pointwise file (refused without one); the OSP
    # stream has one only for fits rendered after B40 (2026-09-27), and is checked exactly
    # when it does and on its posterior mean against ppc_byobs when it does not (R1)
    if (sn %in% c("gear", "trailer") || file.exists(pw)) {
      ck <- bss_block_cv_check(rec$ll, pw)
      if (isTRUE(ck$ok)) {
        pwd <- utils::read.csv(pw, stringsAsFactors = FALSE)
        if (!identical(as.character(pwd$event_date), as.character(dates[rec$day]))) { ck$ok <- FALSE; ck$why <- "the pointwise file's dates do not match the day-to-date map" }
      }
      out$checks[[sn]] <- sprintf("%s: %d obs, %s", sn, ncol(rec$ll),
                                  if (isTRUE(ck$ok)) sprintf("max |lpd diff| %.2g -> MATCHES the committed pointwise file", ck$max_abs_diff) else paste("REFUSED:", ck$why))
      if (!isTRUE(ck$ok)) { out$ok <- FALSE; break }
      n_checked <- n_checked + 1L
    } else {
      ck <- .osp_mean_check(draws, sd, dir, fit, rec)
      if (isTRUE(ck$skip)) { out$checks[[sn]] <- sprintf("%s: %d obs, NOT SCORED (%s)", sn, ncol(rec$ll), ck$why); next }
      out$checks[[sn]] <- sprintf("%s: %d obs, no committed pointwise file; posterior-mean count against ppc_byobs %s", sn, ncol(rec$ll),
                                  if (isTRUE(ck$ok)) sprintf("within %.1f%% on every observation -> approx-checked", 100 * ck$max_rel) else paste("REFUSED:", ck$why))
      if (!isTRUE(ck$ok)) { out$ok <- FALSE; break }
    }
    recs[[sn]] <- list(ll = rec$ll, block = bss_block_weeks(dates[rec$day]))
  }
  if (!out$ok) { out$why <- "a stream failed its check (R1); nothing scored or written for this fit"; return(out) }
  if (n_checked == 0L) { out$ok <- FALSE; out$why <- "no gear or trailer stream could be checked against a committed file (R1); nothing scored"; return(out) }
  if (!isTRUE(score)) return(out)
  tabs <- tryCatch(bss_block_cv_fit(recs, k_max = K_MAX), error = function(e) e)
  if (inherits(tabs, "error")) { out$ok <- FALSE; out$why <- paste("block PSIS failed:", conditionMessage(tabs)); return(out) }
  for (sn in names(tabs)) {
    tab <- cbind(data_type = sn, tabs[[sn]], stringsAsFactors = FALSE)
    if (write) utils::write.csv(tab, file.path(dir, sprintf("loo_block_%s_%s.csv", sn, fit)), row.names = FALSE)
    out$tables[[sn]] <- tab
  }
  out
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
banner(sprintf("MARINE HAZARD COVARIATES: LEAVE-ONE-WEEK-OUT BLOCK CV  |  DRY_RUN = %s", DRY_RUN))
DIRS <- vapply(RUNGS, find_outdir, character(1))
for (r in names(RUNGS)) cat(sprintf("  %s  %-12s %s\n", r, RUNGS[[r]], if (is.na(DIRS[[r]])) "FOLDER NOT FOUND" else DIRS[[r]]))
cat(sprintf("  loo package: %s\n", if (requireNamespace("loo", quietly = TRUE)) "present" else "ABSENT (install.packages('loo'))"))
have_rds <- sapply(names(RUNGS), function(r) if (is.na(DIRS[[r]])) NA else
  sum(file.exists(file.path(DIRS[[r]], sprintf("ppc_draws_%s.rds", FITS)))))
cat(sprintf("  ppc_draws_<fit>.rds present per rung: %s (of 4 fits each)\n", paste(sprintf("%s=%s", names(have_rds), have_rds), collapse = ", ")))

if (isTRUE(DRY_RUN)) {
  cat("\n  DRY RUN. Inventory above; the reconstruction check runs on the first rung and fit that\n  has its draws, nothing is scored and nothing is written. Set DRY_RUN <- FALSE to compute.\n")
  done <- FALSE
  for (r in names(RUNGS)) for (fk in names(FITS)) {
    if (done || is.na(DIRS[[r]])) next
    if (!file.exists(file.path(DIRS[[r]], sprintf("ppc_draws_%s.rds", FITS[[fk]])))) next
    bt <- block_tables(DIRS[[r]], fk, write = FALSE, score = FALSE)
    cat(sprintf("\n  reconstruction check (R1), %s / %s: %s\n", r, fk, if (isTRUE(bt$ok)) "EVALUABLE" else paste("REFUSED:", bt$why)))
    for (ck in bt$checks) cat("   ", ck, "\n")
    done <- TRUE
  }
  if (!done) cat("\n  No rung has its draws here; run this where the rungs rendered.\n")
}

TAB <- list(); ROWS <- list(); PR <- list(); ROWSX <- list()
if (!isTRUE(DRY_RUN)) for (r in names(RUNGS)) {
  if (is.na(DIRS[[r]])) { V1row(r, "the rung folder is present", "absent", "present", "REVIEW", "Nothing to score."); next }
  for (fk in names(FITS)) {
    bt <- block_tables(DIRS[[r]], fk)
    if (!isTRUE(bt$ok)) {
      V1row(r, sprintf("R1: %s draws reconstruct the run's own pointwise lpd", fk),
            paste(c(bt$why, unlist(bt$checks)), collapse = "; "), "every gear / trailer stream matches its committed file; OSP means within 10% of ppc_byobs", "REVIEW",
            "Without the check a wrong reconstruction would be scored as a result. Not evaluated.")
      next
    }
    V1row(r, sprintf("R1: %s draws reconstruct the run's own pointwise lpd", fk), paste(unlist(bt$checks), collapse = "; "),
          "every gear / trailer stream matches its committed file; OSP means within 10% of ppc_byobs", "PASS",
          "The block scores below rest on the same likelihood the run reported.")
    check_of <- function(sn) {
      if (identical(sn, "joint")) {
        parts <- vapply(setdiff(names(bt$tables), "joint"), function(s) if (grepl("MATCHES", bt$checks[[s]] %||% "")) "exact" else "approx", character(1))
        return(paste0("joint (", paste(sprintf("%s %s", names(parts), parts), collapse = ", "), ")"))
      }
      if (grepl("MATCHES", bt$checks[[sn]] %||% "")) "exact" else "approx"
    }
    for (sn in names(bt$tables)) {
      t <- bt$tables[[sn]]; TAB[[paste(r, fk, sn)]] <- t
      ROWSX[[paste(r, fk, sn)]] <- if (identical(sn, "joint")) "joint" else check_of(sn)
      ROWS[[length(ROWS) + 1]] <- data.frame(rung = r, fit = fk, stream = sn, n_weeks = nrow(t), n_obs = sum(t$n_obs),
                                             n_reliable = sum(t$reliable), reliable_share = round(mean(t$reliable), 3),
                                             elpd_block_total = round(sum(t$elpd_block), 2),
                                             elpd_block_reliable = round(sum(t$elpd_block[t$reliable]), 2),
                                             lpd_insample_total = round(sum(t$lpd_block, na.rm = TRUE), 2), max_k = round(max(t$pareto_k), 2),
                                             leaveout_streams = t$leaveout_streams[1],
                                             check = check_of(sn), stringsAsFactors = FALSE)
    }
  }
}

if (!isTRUE(DRY_RUN)) banner("THE PAIRED COMPARISONS (R2, R3, R3', R7)")
WK <- list()
if (!isTRUE(DRY_RUN)) for (fk in names(FITS)) for (sn in SCORED_OF[[fk]]) for (p in PAIRS) {
  if (!is.null(p$fits) && !fk %in% p$fits) next
  ta <- TAB[[paste(p$a, fk, sn)]]; tb <- TAB[[paste(p$b, fk, sn)]]
  # a rung whose folder is not there yet (M6 before its render) gets its one "folder absent"
  # row above and no per-stream rows: a pending rung is not an open comparison
  if (is.na(DIRS[[p$a]] %||% NA) || is.na(DIRS[[p$b]] %||% NA)) next
  primary <- identical(sn, "joint") || length(SCORED_OF[[fk]]) == 1L      # R3': the row the rule is applied to
  approx_of <- function(rung, stream) { r <- ROWSX[[paste(rung, fk, stream)]]; !is.null(r) && identical(r, "approx") }
  suffix <- if (identical(sn, "osp") && (approx_of(p$a, "osp") || approx_of(p$b, "osp"))) " (approx-checked stream)"
            else if (identical(sn, "joint") && (approx_of(p$a, "osp") || approx_of(p$b, "osp"))) " (trailer exact, OSP approx-checked on at least one side)"
            else ""
  what_row <- sprintf("R3%s: %s, %s %s%s, held-out weeks against %s%s", if (identical(sn, "joint")) "'" else "", p$what, fk, sn,
                      if (identical(sn, "joint")) " (all effort streams together; the row R3 is applied to)" else " stream", p$a, suffix)
  if (is.null(ta) || is.null(tb)) {
    if (fk %in% c("shore_all_gear", "boat_all_gear") && (sn != "osp" || fk == "boat_all_gear")) {
      V1row(p$b, what_row,
            sprintf("not computable: %s", paste(c(if (is.null(ta)) p$a, if (is.null(tb)) p$b), "has no block table (draws absent or refused)", collapse = "; ")),
            "> +2 paired SE over the reliable weeks", "REVIEW", "A comparison with a missing side is open, not a verdict.")
      PR[[length(PR) + 1]] <- data.frame(pair = sprintf("%s vs %s", p$b, p$a), what = p$what, fit = fk, stream = sn,
                                         n_blocks = NA_integer_, n_used = NA_integer_, n_dropped = NA_integer_, reliable_share = NA_real_,
                                         diff = NA_real_, se = NA_real_, ratio = NA_real_, k_max = K_MAX, max_abs_diff = NA_real_,
                                         n_positive = NA_integer_, median_diff = NA_real_, identical = FALSE, primary = primary,
                                         evaluable = FALSE, verdict = "REVIEW",
                                         reading = "not computable: a side has no block table", stringsAsFactors = FALSE)
    }
    next
  }
  cmp <- bss_block_cv_compare(ta, tb, k_max = K_MAX)
  wk <- bss_block_cv_weeks(ta, tb, k_max = K_MAX)
  WK[[length(WK) + 1]] <- cbind(pair = sprintf("%s vs %s", p$b, p$a), fit = fk, stream = sn, wk, stringsAsFactors = FALSE)
  evaluable <- isTRUE(cmp$reliable_share >= MIN_REL) && isTRUE(cmp$n_used >= MIN_WKS)
  verdict <- if (isTRUE(cmp$identical)) "INFO"                                   # R7: the same fit twice is not a comparison
             else if (!evaluable) "REVIEW" else if (!is.finite(cmp$ratio)) "REVIEW" else if (cmp$ratio > 2) "PASS" else if (cmp$ratio < -2) "FAIL" else "REVIEW"
  reading <- if (isTRUE(cmp$identical)) sprintf("identical fits: %s and %s hold the same %s fit (max |diff| %.1e nats, summation order); nothing to compare (R7)", p$a, p$b, fk, cmp$max_abs_diff)
             else if (!evaluable) sprintf("PSIS cannot carry this comparison (%d of %d weeks reliable; need %.0f%% and %d): refit the excluded weeks or coarsen the block",
                                     cmp$n_used, cmp$n_blocks, 100 * MIN_REL, MIN_WKS)
             else if (verdict == "PASS") "the covariate improves the interpolation on held-out weeks"
             else if (verdict == "FAIL") "the covariate makes the held-out weeks WORSE"
             else "no evidence either way on held-out weeks"
  V1row(p$b, what_row, bss_block_cv_str(cmp, p$a, p$b), "> +2 paired SE over the reliable weeks", verdict, reading)
  PR[[length(PR) + 1]] <- data.frame(pair = sprintf("%s vs %s", p$b, p$a), what = p$what, fit = fk, stream = sn, cmp, primary = primary,
                                     evaluable = evaluable, verdict = verdict, reading = reading, stringsAsFactors = FALSE)
  cat(sprintf("  %-6s %s | %s | %s\n", verdict, fk, sn, bss_block_cv_str(cmp, p$a, p$b)))
}

if (!isTRUE(DRY_RUN) && length(V)) {
  banner("VERDICTS")
  VV <- do.call(rbind, V)
  for (i in seq_len(nrow(VV))) cat(sprintf("  %-6s %-3s %s\n         observed : %s\n", VV$verdict[i], VV$stage[i], VV$criterion[i], VV$observed[i]))
  cat(sprintf("\n  %d PASS  %d FAIL  %d REVIEW  %d INFO\n", sum(VV$verdict == "PASS"), sum(VV$verdict == "FAIL"), sum(VV$verdict == "REVIEW"), sum(VV$verdict == "INFO")))
}
if (!isTRUE(DRY_RUN)) {
  op <- OUT_ROOT
  if (length(ROWS)) merge_csv_by(do.call(rbind, ROWS), file.path(op, "marine_hazard_2026-09-26_blockcv.csv"), key = c("rung", "fit", "stream"))
  if (length(PR))   merge_csv_by(do.call(rbind, PR),   file.path(op, "marine_hazard_2026-09-26_blockcv_pairs.csv"), key = c("pair", "fit", "stream"))
  if (length(WK))   merge_csv_by(do.call(rbind, WK),   file.path(op, "marine_hazard_2026-09-26_blockcv_weeks.csv"), key = c("pair", "fit", "stream", "block"))
  if (length(V))    merge_csv_by(do.call(rbind, V),    file.path(op, "marine_hazard_2026-09-26_blockcv_verdicts.csv"), key = c("stage", "criterion"))
  cat(sprintf("\n  wrote marine_hazard_2026-09-26_blockcv{,_pairs,_weeks,_verdicts}.csv to %s\n", op))
  cat("\n  READ WITH R6 IN MIND: this scores the interpolation of held-out sampled weeks. It says\n")
  cat("  nothing about the port total, and a PASS here plus the ladder's identification is the case\n")
  cat("  for adopting a term; a REVIEW for want of reliable weeks is a case for refits, not a verdict.\n")
}
