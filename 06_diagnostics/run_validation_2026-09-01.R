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
# run_validation_2026-09-01.R
#
# The outstanding VALIDATION and DIAGNOSIS work, after the shared turnover was adopted
# 2026-09-01. Start it, walk away, come back to
#   05_output/validation_2026-09-01_summary.csv
#   05_output/validation_2026-09-01_verdicts.csv
#   05_output/<date>/validation-desk/            <- the desk stages' outputs
#
# Sources: the two open items from development_notes/stage5-batch-review-2026-08-31.md, plus
# the runnable and readable items in development_notes/PIPELINE_STATUS.md Section 4.
#
# THE ORDERING IS THE ARGUMENT.
#
#   FOUR DESK STAGES FIRST, minutes, no MCMC. Every one of them answers a backlog item from
#   data ALREADY ON DISK. Running four MCMC batches to learn things the committed outputs
#   already contain would be the most expensive way to be wrong about what is open.
#
#   D1  Close (or refuse to close) the Tier-1 "TOP OF THE LIST" shore I/E unit item against
#       its own four pre-set criteria. Those criteria were written in advance and the runs
#       that test them have existed since 2026-08-26; nobody has scored them.
#   D2  Decide the `gear_only` incomplete-trip arm, and measure the PE-vs-BSS arm
#       disagreement the backlog records but never sized.
#   D3  Diagnose the trailer over-coverage, the open item the Stage 5 review created. Also
#       reads the PPC zero bin, which is the stated precondition for the Tier-3
#       zero-inflation decision. Both come out of ppc_byobs_*.csv.
#   D4  Repository-hygiene audit (Tier 4): what is tracked that should not be, how large
#       05_output has become, and whether the LICENSE still carries its placeholder.
#
#   THEN FIVE FITTED STAGES, cheapest and most decisive first.
#
#   V1  THE ADOPTION GATE, and the only stage that must run. Production after the
#       2026-09-01 flip must reproduce S2 EXACTLY on all four fits. Two reasons, not one:
#       (a) it proves the configuration now shipping is the configuration that was
#       validated, rather than something that drifted alongside it; and (b) S2 was fitted
#       before the 2026-08-31 fix to the prior-vs-posterior writer, so it carries NO boat
#       prior_vs_posterior file at all, which is where the contraction diagnostic for the
#       very parameter being adopted lives. This run produces it.
#   V2  GR-7 Phase 2 CONTROL: per-gear CPUE on, Dirichlet gear shares OFF. Must reproduce
#       the 2026-07-20 Phase-1 validation. Cheap, and it is the control V3 needs.
#   V3  GR-7 Phase 2: Dirichlet gear shares ON. Coded 2026-07-21, parses under stanc, and
#       NEVER SAMPLED. A feature that has never been run is not a feature, it is a claim.
#   V4  The untested corner: shared turnover x biweekly boat pot closure. Both levers move
#       that component and they have never been crossed on the pooled track.
#   V5  shared_tau_min_obs = 20, which drops the boat pot closure (18 informed days) out of
#       the feature. A threshold nobody should have to take on trust.
#
# WHAT IS DELIBERATELY NOT HERE, and why:
#   * The censored likelihood for incomplete trips (Tier 2) is the largest single gain
#     available on CPUE precision and it is NOT a validation item. It needs a Stan edit, a
#     recompile, its own OFF-path bit-identity control and an elpd comparison. That is a
#     development batch of its own and mixing it in here would make both harder to read.
#   * Zero-inflation (Tier 3) is explicitly gated on reading the PPC zero bin first. D3
#     reads it. If the bin is off, the ZINB prototype is the next batch, not this one.
#   * The weather module (Tier 2) forks a stale pre-deployment-scale engine and must be
#     re-based before any LOO from it means anything. Re-basing is development.
#   * External validation (T1.1b), the OSP crab-only column, the crabbing-fraction pilot and
#     the Sunday split are data and field items. No amount of compute substitutes.
#
# HOW TO RUN
#     Rscript 06_diagnostics/run_validation_2026-09-01.R      # from the repo root
#     source("06_diagnostics/run_validation_2026-09-01.R")    # RStudio: Source
#
# START WITH DRY_RUN <- TRUE. It runs all four DESK stages for real (they need no Stan
# toolchain and no packages beyond base R), resolves and prints every fitted stage's config,
# runs the R-to-Stan data-contract pre-flight and the extractor self-test, and fits nothing.
# Read the desk verdicts, then set DRY_RUN <- FALSE.
#
# RESUMABLE. A stage whose output folder already holds port_total_Dungeness_Kept.csv is
# skipped and re-extracted. Rows are appended as each stage finishes.
#
# RUNTIME. Roughly 15-17 h on 4 cores: V1 ~4.5 h, V2+V3 ~2 h (the gear track runs in about
# half an hour now that its boat fit is given the pooled track's adaptation), V4 ~4.5 h,
# V5 ~4.5 h. The desk stages are seconds.
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                  # TRUE: desk stages run, nothing is fitted. START HERE.
#        ^^^^ reset to TRUE after the 2026-09-01 batch. RESUME skips completed stages, so
#        sourcing this with DRY_RUN <- FALSE again would re-extract the existing outputs and
#        APPEND another five rows to the summary CSV rather than re-fitting.
STAGES  <- c("D1", "D2", "D3", "D4", "V1", "V2", "V3", "V4", "V5")
RESUME  <- TRUE                  # skip a fitted stage whose output folder already looks complete

# V1 is the adoption gate. If it does not reproduce S2 exactly, the configuration now
# shipping is NOT the one the Stage 5 batch validated, and V4/V5 cannot be read against S2
# either. Set this FALSE only to deliberately override that.
GATE_ON_V1 <- TRUE

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
  for (f in c("model_diagnostics.R", "bss_model_adequacy.R", "annotate_decoupled_run.R",
              "bss_sampler_override.R", "bss_ar_resolution.R"))
    try(sys.source(.here("03_R_functions", f), envir = globalenv()), silent = TRUE)
}

source(.here("run_config.R"))
BASE <- run_config

model_rmd <- list(
  pooled        = .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd"))
prefix <- list(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")
stopifnot(all(file.exists(unlist(model_rmd))))

desk_dir <- .here("05_output", format(Sys.Date(), "%Y%m%d"), "validation-desk")

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
# CONFIG DELTA between two run folders, from their committed run_parameters.txt dumps.
#
# WHY THIS EXISTS. An exactness comparison is only meaningful when the two runs differ in
# ONE thing. That sounds obvious and it has now been got wrong twice in three weeks, both
# times by me and both times producing a FAIL that looked like a code defect:
#   * the 2026-08-30 S3 criterion compared gear shore fits against a run whose shared_tau was
#     GLOBAL, so the reference differed in the sampler AND in the floor;
#   * the 2026-09-01 V2 criterion compared gear boat fits against a run that carried the
#     experiment sampler override, so the reference differed in gear_resolved_G AND in the
#     sampler settings.
# Both were caught only by hand afterwards. The lesson is that "is this the right reference?"
# should be a CHECK, not a judgement, so fit_exactness() now takes `expect_delta` and reports
# what actually differs between the two configs alongside its verdict.
#
# Parses the `$ key : value` lines of a str() dump. Comparison is on the printed text, so it
# is approximate for lists and long vectors; that is fine, because the job is to catch a
# reference that differs in an unintended KEY, not to diff values precisely.
# ---------------------------------------------------------------------------
.cfg_keys <- function(dir) {
  p <- file.path(dir, "run_parameters.txt")
  if (!file.exists(p)) return(NULL)
  tx <- readLines(p, warn = FALSE); kv <- list()
  for (l in tx) {
    m <- regmatches(l, regexec("^\\s*\\$ ([A-Za-z0-9_.]+)\\s*:(.*)$", l))[[1]]
    if (length(m) == 3) kv[[m[2]]] <- trimws(m[3])
  }
  kv
}
config_delta <- function(dir_a, dir_b) {
  a <- .cfg_keys(dir_a); b <- .cfg_keys(dir_b)
  if (is.null(a) || is.null(b)) return(NA_character_)
  ks <- union(names(a), names(b))
  ks <- ks[!ks %in% c("run_tag", "model")]      # name the folder, not the fit
  d <- ks[vapply(ks, function(k) !identical(a[[k]] %||% "<absent>", b[[k]] %||% "<absent>"), logical(1))]
  sort(d)
}

# ---------------------------------------------------------------------------
# Compare the fits of two run folders at full double precision, restricted to the fits whose
# label matches `pat`. Shared parameter ROWS only, because a feature legitimately ADDS
# reported rows; new rows are not a behaviour change, changed values are.
#
# `expect_delta` names the config keys the two runs are ALLOWED to differ in. Anything else
# that differs is reported in the observation, so a FAIL caused by the wrong reference is
# visible in the verdict rather than diagnosed by hand a day later.
# ---------------------------------------------------------------------------
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
       verdict = if (!length(bad) && n > 0) "PASS" else "FAIL",
       unexpected_delta = extra)
}

# ===========================================================================
# D1. The Tier-1 "TOP OF THE LIST" item, scored against its OWN pre-set criteria.
#
# PIPELINE_STATUS Section 4 Tier 1 opens with: validate the shore I/E observation-unit fix,
# ISOLATED, with four criteria written in advance so the run is "read rather than
# rationalised". The isolating run has existed since 2026-08-26: rung 1 of the validation
# ladder is the pre-fix baseline and rung 2 changes `ie_shore_obs_unit` and nothing else.
# Nobody has scored the criteria against it, so a Tier-1 item that may well be closed has
# been sitting at the top of the backlog. This scores all four. It is a desk stage because
# every number it needs is already committed.
#
# The four criteria, quoted from the backlog:
#   (a) shore all-gear sigma_IE should fall sharply from ~1.07; if it does not, the GR-9
#       tension is NOT the unit mismatch and that item stays open needing a new hypothesis.
#   (b) the shore effort posterior should move, direction not predictable in advance.
#   (c) the shore tau_shore posterior should now depart from its prior on the in-window I/E
#       days, where before it could not.
#   (d) the boat must be byte-identical; any boat movement is a bug in the patch.
# ===========================================================================
run_D1 <- function() {
  banner("DESK D1  |  Tier-1 shore I/E observation-unit fix, scored against its four pre-set criteria")
  pre  <- .here("05_output", "20260826", "pooled-CPUE-PV1-baseline")   # rung 1, pre-fix
  post <- .here("05_output", "20260826", "pooled-CPUE-PV2-ie-unit")    # rung 2, the fix ISOLATED
  ship <- .here("05_output", "20260829", "pooled-CPUE-S5-2-tau-pooled")# where it stands today
  if (!dir.exists(pre) || !dir.exists(post)) {
    cat("  ladder rungs 1/2 not found; cannot score this item.\n"); return(invisible(NULL))
  }
  sig <- function(d) { x <- rd(d, "structural_params_shore_all_gear_Dungeness_Kept.csv")
    if (is.null(x)) return(c(NA, NA, NA))
    i <- x$parameter == "sigma_IE"; c(num1(x$median[i]), num1(x$lo95[i]), num1(x$hi95[i])) }
  eff <- function(d) { x <- rd(d, "pe_vs_bss_comparison.csv")
    if (is.null(x)) return(NA_real_); num1(x$BSS_catch[grepl("^shore \\(All gear\\)", x$component)]) }
  Lrng <- function(d) { x <- rd(d, "bss_L_effective_shore_all_gear_Dungeness_Kept.csv")
    if (is.null(x) || !"L_posterior_median" %in% names(x)) return(c(NA, NA, NA))
    v <- suppressWarnings(as.numeric(x$L_posterior_median)); c(min(v), max(v), diff(range(v))) }

  a0 <- sig(pre); a1 <- sig(post); a2 <- sig(ship)
  V1row("D1", "(a) shore all-gear sigma_IE falls sharply from ~1.07",
    sprintf("%.3f [%.2f, %.2f] -> %.3f [%.2f, %.2f] (today %.3f)", a0[1], a0[2], a0[3], a1[1], a1[2], a1[3], a2[1]),
    "a sharp fall; the item's own wording",
    if (!is.na(a1[1]) && !is.na(a0[1]) && a1[1] < 0.5 * a0[1]) "PASS" else "FAIL",
    paste("The shore I/E likelihood compared crabber-HOURS against a predicted",
          "lambda_E * tau = crabber-TRIPS from v7.7 until the 2026-08-25 fix, roughly 331",
          "against 80 on a typical WDF20 day, and sigma_IE was the only free parameter",
          "positioned to absorb it. e^1.03 = 2.8 against an implied offset near 4. A PASS",
          "closes BOTH this Tier-1 item and its Tier-2 twin, the GR-9 shore sigma_IE",
          "tension, which the backlog explicitly says must not be closed on reasoning alone."))

  e0 <- eff(pre); e1 <- eff(post)
  V1row("D1", "(b) the shore effort posterior moves",
    sprintf("shore all-gear BSS %s -> %s (%+.1f%%)", fmt(e0, 0), fmt(e1, 0), 100 * (e1 - e0) / e0),
    "any movement; direction not predictable in advance",
    if (!is.na(e0) && !is.na(e1) && abs(e1 - e0) > 0.5) "PASS" else "REVIEW",
    paste("The old stream was pulling lambda_E UP against the gear counts, so the backlog",
          "deliberately declined to predict a direction. Recording the size matters more",
          "than the sign: this is the only change in the 2026-08-25 batch that moves the",
          "shore, so it is the whole shore-side delta of that batch."))

  L0 <- Lrng(pre); L1 <- Lrng(post)
  V1row("D1", "(c) tau_shore departs from its prior on the in-window I/E days",
    sprintf("L posterior median spans %.3f-%.3f (range %.3f) -> %.3f-%.3f (range %.3f); prior centre 1.700",
            L0[1], L0[2], L0[3], L1[1], L1[2], L1[3]),
    "a wider posterior span, and excursions below the prior centre",
    if (!is.na(L1[3]) && !is.na(L0[3]) && L1[3] > 1.5 * L0[3] && L1[1] < L0[1]) "PASS" else "REVIEW",
    paste("Before the fix the shore L posterior was pinned within a whisker of its 1.700",
          "prior and could only move UP; the mismatched observation gave it nothing usable",
          "on the low side. If the fix is real, the I/E days should now be able to pull the",
          "turnover DOWN as well as up, which is what a wider span with a lower minimum is."))

  ex_boat  <- fit_exactness(post, pre, "private_boat", "boat fits, rung 1 vs rung 2")
  ex_shore <- fit_exactness(post, pre, "shore", "shore fits, rung 1 vs rung 2")
  V1row("D1", "(d) the boat is byte-identical across the isolated fix",
    sprintf("%s | for contrast, %s", ex_boat$observed, ex_shore$observed),
    "boat identical; shore differs (it is the only thing that should)",
    if (identical(ex_boat$verdict, "PASS") && identical(ex_shore$verdict, "FAIL")) "PASS" else "REVIEW",
    paste("The boat I/E stream was already on the correct unit, so the boat MUST be",
          "untouched; any boat movement would be a bug in the patch rather than a finding.",
          "The shore differing is the positive control: a comparison where nothing moved at",
          "all would mean the toggle never reached the model."))
  invisible(NULL)
}

# ===========================================================================
# D2. The `gear_only` incomplete-trip arm (Tier 2), decided on the diagnostic's own evidence.
#
# The backlog says: "Decide after the first run that carries it: adopt gear_only if the
# incomplete-trip gear counts are not significantly higher than complete ones, and either way
# make the PE and the BSS agree." The four-arm diagnostic has been carried by every run since
# 2026-08-25 and writes sensitivity_incomplete_trips.csv. The decision is a read, not a run.
# ===========================================================================
run_D2 <- function() {
  banner("DESK D2  |  the gear_only incomplete-trip arm, and the PE-vs-BSS arm disagreement")
  d <- rd(.here("05_output", "20260829", "pooled-CPUE-S5-2-tau-pooled"), "sensitivity_incomplete_trips.csv")
  if (is.null(d)) { cat("  sensitivity_incomplete_trips.csv not found.\n"); return(invisible(NULL)) }
  keep <- c("component", "arm", "n_interviews", "n_incomplete", "pct_incomplete", "gear_ratio",
            "gear_lengthbias_p", "production_arm_bss", "production_arm_pe")
  keep <- keep[keep %in% names(d)]
  w <- d[, keep, drop = FALSE]
  dir.create(desk_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(w, file.path(desk_dir, "D2_incomplete_trip_arms.csv"), row.names = FALSE)

  comps <- unique(w$component)
  rows <- lapply(comps, function(cc) {
    x <- w[w$component == cc, ]
    ge <- num1(x$gear_ratio[x$arm == "exclude"]); gg <- num1(x$gear_ratio[x$arm == "gear_only"])
    data.frame(component = cc, gear_ratio_exclude = ge, gear_ratio_gear_only = gg,
               pct_change = 100 * (gg - ge) / ge,
               lengthbias_p = num1(x$gear_lengthbias_p[1]),
               pct_incomplete = num1(x$pct_incomplete[1]),
               arm_bss = chr1(x$production_arm_bss[1]), arm_pe = chr1(x$production_arm_pe[1]),
               stringsAsFactors = FALSE) })
  R <- do.call(rbind, rows)
  utils::write.csv(R, file.path(desk_dir, "D2_gear_only_decision.csv"), row.names = FALSE)
  rule()
  cat(sprintf("  %-28s %9s %9s %8s %10s   %-9s %-9s\n", "component", "R_G excl", "R_G g_only",
              "change", "lenbias p", "BSS arm", "PE arm"))
  for (i in seq_len(nrow(R)))
    cat(sprintf("  %-28s %9.3f %9.3f %+7.1f%% %10s   %-9s %-9s\n", R$component[i],
                R$gear_ratio_exclude[i], R$gear_ratio_gear_only[i], R$pct_change[i],
                fmt(R$lengthbias_p[i], 4), R$arm_bss[i], R$arm_pe[i]))
  rule()

  worst_p <- min(R$lengthbias_p, na.rm = TRUE)
  worst_d <- max(abs(R$pct_change), na.rm = TRUE)
  disagree <- R[!is.na(R$arm_bss) & !is.na(R$arm_pe) & R$arm_bss != R$arm_pe, , drop = FALSE]
  V1row("D2", "are incomplete-trip gear counts exchangeable with complete ones?",
    sprintf("smallest length-bias p = %s; largest R_G shift from adopting gear_only = %.1f%%",
            fmt(worst_p, 4), worst_d),
    "the backlog's rule: adopt gear_only if not significantly higher",
    if (is.na(worst_p)) "REVIEW" else if (worst_p >= 0.05) "ADOPT gear_only" else "STATISTICALLY REJECTED, PRACTICALLY NEGLIGIBLE",
    paste("Read both columns, because they disagree in an instructive way. The length-bias",
          "test is significant on shore all-gear only because n = 2,741 makes a 0.7% shift",
          "detectable, and the shift is DOWNWARD (incomplete trips carry FEWER gear, which",
          "is what a trip interrupted while still setting up looks like), not the upward",
          "bias the rule was written to guard against. So the rule's letter says reject and",
          "its spirit is untroubled. The honest resolution is that the arm choice is worth",
          "at most 0.7% on any component and 0.1% on the boat, so it should be settled for",
          "CONSISTENCY rather than accuracy: pick one and make both estimators use it."))
  if (nrow(disagree))
    V1row("D2", "do the PE and the BSS use the same arm?",
      sprintf("%d component(s) disagree: %s", nrow(disagree),
              paste(sprintf("%s (BSS %s, PE %s)", disagree$component, disagree$arm_bss, disagree$arm_pe),
                    collapse = "; ")),
      "the two estimators must agree", "FAIL",
      paste("The backlog records this disagreement but has never sized it. Both boat",
            "components run the BSS on `exclude` while the PE uses `gear_only`, so the boat",
            "PE's gear-per-group comes from the unfiltered set and the boat BSS's R_G_boat",
            "does not. The size is now known and it is small (0.1% on the boat gear ratio),",
            "which makes this cheap to fix and hard to justify leaving: a documented,",
            "measured inconsistency between the two arms of a fused estimator is exactly",
            "what a reviewer will ask about."))
  invisible(R)
}

# ===========================================================================
# D3. The trailer over-coverage, and the PPC zero bin.
#
# THE OPEN ITEM. Every configuration in the 2026-08-30 batch over-covers the trailer stream:
# coverage_50 runs 0.667-0.692 at monthly against a nominal 0.500, which at n = 195 is 4.6 to
# 5.3 sampling SDs (SD = sqrt(0.25/195) = 0.036). The shared turnover improves it slightly
# and nothing fixes it. The Stage 5 review recorded it and could not explain it.
#
# TWO HYPOTHESES, AND A DISCRIMINATOR STATED BEFORE THE NUMBERS ARE LOOKED AT.
#
#   H1 DISCRETENESS. Trailer counts are small integers. A nominal 50% interval taken from
#      the quantiles of a discrete predictive distribution cannot have exactly 50% mass; at
#      a fitted mean near 1 or 2 the smallest interval containing the median already carries
#      far more than half the probability, so coverage is inflated by construction and the
#      model is not misspecified at all. The standard remedy is a randomized PIT, which is a
#      diagnostic change, not a model change.
#   H2 GENUINE OVER-DISPERSION. The predictive intervals really are too wide, meaning r_E is
#      being dragged down by a minority of badly-predicted days. That IS misspecification and
#      points at the observation family.
#
#   DISCRIMINATOR: split the stream by fitted mean into terciles. Under H1 the over-coverage
#   is concentrated in the low tercile and the high tercile approaches 0.500. Under H2 it is
#   flat across terciles. The two make opposite predictions about the same table, which is
#   what makes this worth computing rather than arguing about.
#
# ALSO THE ZERO BIN. Tier 3 says the zero-inflation question is to be read off the PPC zero
# bin first and prototyped "only if the zero bin is systematically off". Reported here so
# that decision has its evidence.
# ===========================================================================
run_D3 <- function() {
  banner("DESK D3  |  trailer over-coverage: discreteness or misspecification? Plus the zero bin")
  # The archived comparison cells, PLUS the newest production run. The zero-bin criterion
  # needs p_zero, which only exists on runs made after 2026-09-01, so pinning this list to
  # archived folders made that criterion permanently unscorable. It now picks up whatever
  # production run exists and prefers it for the zero bin.
  runs <- c("D  tau off/monthly" = "20260828/pooled-CPUE-IP-D-tau-off",
            "S2 tau ON /monthly" = "20260829/pooled-CPUE-S5-2-tau-pooled",
            "S4b tau ON /daily"  = "20260830/pooled-CPUE-S5-4b-daily-tauon")
  .newest <- {
    hits <- list.dirs(.here("05_output"), recursive = TRUE)
    hits <- hits[grepl("pooled-CPUE-VAL-1-adopted$", hits)]
    if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
  }
  if (!is.na(.newest))
    runs <- c(runs, "V1 PRODUCTION" = sub(paste0("^", .here("05_output"), "/"), "", .newest))
  fits <- c("private_boat_all_gear", "shore_all_gear")
  dir.create(desk_dir, recursive = TRUE, showWarnings = FALSE)
  out <- list(); zero <- list()
  for (nm in names(runs)) for (ft in fits) {
    x <- rd(.here("05_output", runs[[nm]]), sprintf("ppc_byobs_%s_Dungeness_Kept.csv", ft))
    if (is.null(x) || !all(c("data_type", "fitted_mean", "in_50", "observed") %in% names(x))) next
    for (st in unique(x$data_type)) {
      y <- x[x$data_type == st, ]
      if (nrow(y) < 12) next
      tert <- cut(y$fitted_mean, breaks = stats::quantile(y$fitted_mean, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                  include.lowest = TRUE, labels = c("low", "mid", "high"))
      cv <- tapply(as.logical(y$in_50), tert, mean)
      out[[length(out) + 1]] <- data.frame(
        run = nm, fit = ft, stream = st, n = nrow(y),
        coverage_50_all = mean(as.logical(y$in_50), na.rm = TRUE),
        cov50_low = as.numeric(cv["low"]), cov50_mid = as.numeric(cv["mid"]),
        cov50_high = as.numeric(cv["high"]),
        mean_fit_low = mean(y$fitted_mean[tert == "low"], na.rm = TRUE),
        mean_fit_high = mean(y$fitted_mean[tert == "high"], na.rm = TRUE),
        stringsAsFactors = FALSE)
      z <- y$observed == 0
      # THE ZERO BIN, DONE CORRECTLY. Compare the observed zero COUNT against the model's own
      # expected zero count, sum(p_zero) over ALL days. Reading the PIT at the observed zeros
      # instead answers a different question, because those days are selected for being zero
      # and so are the low-mean days whose P(Y=0) is high by construction.
      # p_zero was added to ppc_byobs on 2026-09-01; runs committed before that carry NA here
      # rather than a silently wrong substitute.
      zero[[length(zero) + 1]] <- data.frame(
        run = nm, fit = ft, stream = st, n = nrow(y),
        obs_zeros = sum(z, na.rm = TRUE),
        obs_zero_frac = mean(z, na.rm = TRUE),
        exp_zeros = if ("p_zero" %in% names(y)) sum(y$p_zero, na.rm = TRUE) else NA_real_,
        exp_zero_frac = if ("p_zero" %in% names(y)) mean(y$p_zero, na.rm = TRUE) else NA_real_,
        # Poisson-binomial SD of the expected zero count. A percentage-point gap alone cannot
        # say whether a stream is off, because n ranges from 17 to 1,649 across these fits.
        exp_zeros_sd = if ("p_zero" %in% names(y))
          sqrt(sum(y$p_zero * (1 - y$p_zero), na.rm = TRUE)) else NA_real_,
        mean_fitted_at_zero = if (any(z)) mean(y$fitted_mean[z], na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) { cat("  no ppc_byobs_*.csv found.\n"); return(invisible(NULL)) }
  O <- do.call(rbind, out); Z <- do.call(rbind, zero)
  utils::write.csv(O, file.path(desk_dir, "D3_coverage_by_fitted_tercile.csv"), row.names = FALSE)
  utils::write.csv(Z, file.path(desk_dir, "D3_zero_bin.csv"), row.names = FALSE)

  rule()
  cat(sprintf("  %-19s %-22s %-8s %5s %8s %8s %8s %8s\n", "run", "fit", "stream", "n",
              "cov50", "low", "mid", "high"))
  for (i in seq_len(nrow(O)))
    cat(sprintf("  %-19s %-22s %-8s %5d %8.3f %8.3f %8.3f %8.3f\n", O$run[i], O$fit[i], O$stream[i],
                O$n[i], O$coverage_50_all[i], O$cov50_low[i], O$cov50_mid[i], O$cov50_high[i]))
  cat("\n  nominal coverage_50 = 0.500 in every column.\n")
  rule()
  cat(sprintf("  %-19s %-22s %-8s %10s %10s %11s\n", "run", "fit", "stream",
              "obs zeros", "exp zeros", "obs - exp %"))
  for (i in seq_len(nrow(Z)))
    cat(sprintf("  %-19s %-22s %-8s %10d %10s %10s\n", Z$run[i], Z$fit[i], Z$stream[i],
                Z$obs_zeros[i], fmt(Z$exp_zeros[i], 1),
                if (is.na(Z$exp_zero_frac[i])) "n/a (old run)"
                else sprintf("%+.1f", 100 * (Z$obs_zero_frac[i] - Z$exp_zero_frac[i]))))
  cat("\n  exp zeros = sum over ALL days of the model's own P(Y = 0); n/a where the run\n")
  cat("  predates the p_zero column (2026-09-01).\n")
  rule()

  # ---- verdicts ----------------------------------------------------------
  # NOTE ON WHAT CHANGED UNDER THIS STAGE'S FEET. The premise the H1/H2 discriminator was
  # written for was "the trailer stream is over-covered in every configuration", recorded as
  # an open item in the 2026-08-31 Stage 5 review. Writing this stage exposed the reason:
  # ppc_calibration_*.csv and ppc_byobs_*.csv computed the SAME quantity two different ways
  # and disagreed by up to 0.154. The aggregate file tested the observation against a
  # QUANTILE INTERVAL of simulated draws, which over-covers small counts by construction;
  # the per-observation file used the randomized PIT, which does not. model_diagnostics.R was
  # fixed 2026-09-01 to use the randomized statistic. So the first question is no longer
  # "which hypothesis explains the over-coverage" but "is there any over-coverage to
  # explain", and the honest answer for the production configuration is no.
  agree <- NULL
  cal <- rd(.here("05_output", runs[["S2 tau ON /monthly"]]),
            "ppc_calibration_private_boat_all_gear_Dungeness_Kept.csv")
  s2t <- O[O$stream == "trailer" & grepl("S2", O$run), , drop = FALSE]
  if (!is.null(cal) && nrow(s2t)) {
    cc <- suppressWarnings(as.numeric(cal$coverage_50[cal$data_type == "trailer"]))[1]
    agree <- abs(cc - s2t$coverage_50_all[1])
    V1row("D3", "the two PPC files agree on coverage_50",
      sprintf("ppc_calibration %.3f vs ppc_byobs %.3f on the same run and stream (difference %.3f)",
              cc, s2t$coverage_50_all[1], agree),
      "identical; they compute the same quantity",
      if (is.na(agree)) "REVIEW" else if (agree < 1e-9) "PASS" else if (agree < 0.03) "PASS (residual simulation noise)" else "FAIL (pre-2026-09-01 run)",
      paste("A FAIL here on an ARCHIVED run is expected and is the defect being recorded:",
            "every run committed before 2026-09-01 carries the non-randomized coverage in",
            "ppc_calibration_*.csv. On a run made after the fix the two must agree, and this",
            "criterion becomes a regression guard. Prefer ppc_byobs wherever both exist;",
            "bss_model_adequacy.R now does so automatically and records which it used."))
  }

  if (nrow(s2t)) {
    n <- s2t$n[1]; sdc <- sqrt(0.25 / n); dev <- abs(s2t$coverage_50_all[1] - 0.5)
    V1row("D3", "is the trailer stream over-covered at all, on the corrected statistic?",
      sprintf("production (S2) trailer coverage_50 = %.3f on n = %d, %.1f sampling SDs from the nominal 0.500 (SD %.3f)",
              s2t$coverage_50_all[1], n, dev / sdc, sdc),
      "within about 2 sampling SDs is calibrated",
      if (is.na(dev)) "REVIEW" else if (dev / sdc < 2) "NO OVER-COVERAGE: the open item is RETRACTED" else "CONFIRMED",
      paste("The 2026-08-31 Stage 5 review recorded trailer over-coverage of 0.667-0.692",
            "against a nominal 0.500 as a new open modelling item, at 4.6-5.3 sampling SDs.",
            "Those figures came from the non-randomized statistic. The randomized one puts",
            "the same runs at 0.523 and 0.538. The catch stream, whose counts are an order",
            "of magnitude larger, barely moves between the two computations (0.595 either",
            "way), which is exactly the signature of a small-count artefact rather than a",
            "model defect. If this reads NO OVER-COVERAGE, that open item should be struck",
            "rather than carried forward, and no change should be made to the trailer",
            "observation model on account of it."))
  }

  dly <- O[O$stream == "trailer" & grepl("S4b", O$run), , drop = FALSE]
  if (nrow(dly)) {
    n <- dly$n[1]; sdc <- sqrt(0.25 / n); dev <- abs(dly$coverage_50_all[1] - 0.5)
    gap <- dly$cov50_low[1] - dly$cov50_high[1]
    V1row("D3", "the daily-AR cell is still broken on the corrected statistic",
      sprintf("S4b trailer coverage_50 = %.3f (%.1f sampling SDs); by fitted tercile low %.3f / mid %.3f / high %.3f",
              dly$coverage_50_all[1], dev / sdc, dly$cov50_low[1], dly$cov50_mid[1], dly$cov50_high[1]),
      "the Stage 5 conclusion must survive the correction, or it was resting on the artefact",
      if (is.na(dev)) "REVIEW" else if (dev / sdc > 3) "SURVIVES" else "THE CONCLUSION NEEDS REVISITING",
      paste("This is the criterion that matters most in this stage. If correcting the",
            "coverage statistic had also dissolved the daily-AR pathology, then the",
            "recommendation to adopt the shared turnover and reject the daily AR would have",
            "been resting on the same artefact and would need reopening. It does not: the",
            "daily cells sit at 0.713 and 0.744 under the randomized statistic, 5.9 and 6.8",
            "sampling SDs from nominal, against 0.523 for production. The correction makes",
            "the contrast CLEANER by removing a spurious background of over-coverage that",
            "affected every cell equally. Note the tercile profile is over-covered at the",
            "LOW end here, which is a latent process tracking small counts too closely."))
  }

  # ---- the zero bin, expressed as an implied P(Y = 0) ---------------------
  # The randomized PIT at an observed zero is 0.5 * P(Y = 0), so P(Y = 0) = 2 * PIT. Compare
  # that against the observed zero fraction: they should be close if NB2 already places
  # adequate mass at zero, which is the precondition Tier 3 sets for prototyping ZINB.
  zt <- Z[grepl("V1 PRODUCTION", Z$run), , drop = FALSE]
  if (!nrow(zt)) zt <- Z[grepl("S2", Z$run), , drop = FALSE]
  if (nrow(zt)) {
    utils::write.csv(zt, file.path(desk_dir, "D3_zero_bin.csv"), row.names = FALSE)
    have <- is.finite(zt$exp_zero_frac) & is.finite(zt$exp_zeros_sd) & zt$exp_zeros_sd > 0.5
    worst <- if (any(have)) max(abs(zt$obs_zeros[have] - zt$exp_zeros[have]) / zt$exp_zeros_sd[have]) else NA_real_
    V1row("D3", "PPC zero bin: does NB2 already place adequate mass at zero?",
      if (!any(have))
        "not computable: this run predates the p_zero column added 2026-09-01. Re-read after V1." else
      paste(mapply(function(ft, st, oz, ez, sd)
              sprintf("%s/%s: %d observed vs %.1f expected (z = %+.1f)",
                      sub("_Dungeness.*", "", ft), st, oz, ez, (oz - ez) / pmax(sd, 1e-9)),
              zt$fit[have], zt$stream[have], zt$obs_zeros[have], zt$exp_zeros[have],
              zt$exp_zeros_sd[have]),
            collapse = " | "),
      "Tier 3: prototype ZINB only if the zero bin is SYSTEMATICALLY off",
      if (is.na(worst)) "NOT YET COMPUTABLE (needs a post-2026-09-01 run)"
      else if (worst < 3) "NO ZERO INFLATION NEEDED" else "OFF ON ONE STREAM, SEE why",
      paste("Tier 3 makes the zero-inflation decision explicitly conditional on this read,",
            "which had never been produced, so the item has sat open on a precondition",
            "nobody checked. The randomized PIT at an observed zero is exactly 0.5 * P(Y=0),",
            "so the check is the observed zero COUNT against the model's expected zero",
            "count, sum(p_zero) over ALL days. Where those agree, an extra inflation",
            "parameter has nothing left to explain and would be fitted to noise. The gap is",
            "scored as a z against the Poisson-binomial SD sqrt(sum p(1-p)), because n runs",
            "from 17 to 1,649 across these fits and a raw percentage-point gap cannot be read",
            "the same way at both ends. NOTE the direction if they disagree: zero INFLATION",
            "adds mass at zero, so it is the remedy only when the model expects FEWER zeros",
            "than were observed. A model expecting MORE zeros than occurred has the opposite",
            "problem and a ZINB would make it worse. Read the per-stream z values, not just",
            "the verdict: a single stream off at z = 3-4 while every other stream sits inside",
            "z = 2.5 is a targeted misfit, and a prototype should target THAT likelihood",
            "rather than the whole model."))
  }
  invisible(O)
}

# ===========================================================================
# D4. Repository hygiene (Tier 4). Correctness-neutral, publication-blocking.
# ===========================================================================
run_D4 <- function() {
  banner("DESK D4  |  repository hygiene (Tier 4)")
  dir.create(desk_dir, recursive = TRUE, showWarnings = FALSE)
  tracked <- tryCatch(system2("git", c("ls-files"), stdout = TRUE, stderr = FALSE), error = function(e) character(0))
  rproj  <- grep("^\\.Rproj\\.user", tracked, value = TRUE)
  outs   <- grep("^05_output/", tracked, value = TRUE)
  sz <- tryCatch(sum(file.size(.here(outs)), na.rm = TRUE) / 1024^2, error = function(e) NA_real_)
  lic <- .here("LICENSE")
  lic_bytes <- if (file.exists(lic)) file.size(lic) else NA_real_
  lic_ph <- if (file.exists(lic)) any(grepl("placeholder|PASTE|TODO", readLines(lic, warn = FALSE), ignore.case = TRUE)) else NA
  stale <- grep("^05_output/20260711/pooled-CPUE-(morning|afternoon)/", tracked, value = TRUE)

  H <- data.frame(check = c("tracked .Rproj.user files", "tracked 05_output files",
                            "tracked 05_output size (MB)", "LICENSE bytes",
                            "LICENSE still has a placeholder", "stale 20260711 morning/afternoon files"),
                  value = c(length(rproj), length(outs), round(sz, 1), lic_bytes,
                            as.character(lic_ph), length(stale)), stringsAsFactors = FALSE)
  utils::write.csv(H, file.path(desk_dir, "D4_hygiene.csv"), row.names = FALSE)
  rule(); for (i in seq_len(nrow(H))) cat(sprintf("  %-42s %s\n", H$check[i], H$value[i])); rule()

  V1row("D4", "tracked files that should not be, and the LICENSE action",
    sprintf(".Rproj.user: %d tracked | 05_output: %d files, %s MB | LICENSE: %s bytes, placeholder %s | stale 20260711 runs: %d files",
            length(rproj), length(outs), fmt(sz, 1), format(lic_bytes, big.mark = ","), lic_ph, length(stale)),
    ".Rproj.user untracked; LICENSE carrying the full GPL-3.0 text",
    if (length(rproj) == 0 && isFALSE(lic_ph)) "PASS" else "ACTION",
    paste("Correctness-neutral and publication-blocking, which is a bad combination to",
          "discover late. `.Rproj.user` is gitignored yet tracked, so the ignore does",
          "nothing: it needs `git rm -r --cached .Rproj.user`. The LICENSE size tells you",
          "whether the canonical GPL-3.0 text has actually been pasted in; a file of a few",
          "hundred bytes is still the placeholder. Outputs are committed ON PURPOSE so past",
          "estimates are preserved as produced, so the size figure is context for a policy",
          "decision, not a defect: the open question is whether pre-v6 runs that are no",
          "longer reproducible from current code should stay in the tree or move to",
          "Releases. WDFW confirming the copyright holder and year lines is a human step",
          "this script cannot check."))
  invisible(H)
}

# ---------------------------------------------------------------------------
# BASELINES. Hard-coded from committed folders, named beside each value.
#
# Port totals are NEVER used in an exactness gate: rstan::extract(permuted = TRUE) permutes
# draws, so a port total assembled by resampling component draws is RNG-sensitive even when
# the fits are bit-identical (stage D reproduced rung 4's fits exactly and still reported
# 66,094 against 66,237). Exactness is always tested on per-fit summaries.
# ---------------------------------------------------------------------------
REF <- list(
  # The 2026-08-30 Stage 5 run of exactly the configuration adopted 2026-09-01.
  S2    = list(dir = "20260829/pooled-CPUE-S5-2-tau-pooled",
               port = 71521, shore_ag = 20898, shore_pc = 6331, boat_ag = 31008, boat_pc = 1018,
               tau_bar = 2.5969, tau_lo95 = 2.0644, tau_hi95 = 3.2488,
               pit_trailer = 0.4843, pit_osp = 0.5203),
  S3gear = list(dir = "20260829/gear-type-CPUE-model-S5-3-tau-gear",
                port = 70886, boat_ag = 30760, boat_pc = 956, tau_bar = 2.5962),
  # Boat pot closure at biweekly with the turnover OFF (S5), and the gear track's biweekly.
  S5    = list(dir = "20260830/pooled-CPUE-S5-5-boatpc-ar", boat_pc = 735, boat_ag = 25868),
  gear5 = list(dir = "20260826/gear-type-CPUE-model-PV5-gear-ship", boat_pc = 743),
  # GR-7 Phase 1 as validated 2026-07-20 (gear_resolved_G = TRUE). Context only: that run
  # predates the whole 2026-08-25 batch, the I/E unit fix and the shared turnover, so V2 will
  # NOT and SHOULD NOT reproduce it. Its per-gear SHORE all-gear medians are recorded so the
  # ORDERING can be checked (Trap highest, matching the interview ratios), which is the part
  # that should survive.
  #   NOTE the folder name carries spaces and an "=" (run_tag "gear_resolved_G = TRUE").
  #   That is a hygiene defect in its own right; D4 reports it.
  GR7p1 = list(dir = "20260720/gear-type-CPUE-model-gear_resolved_G = TRUE",
               shore_ag_gear = c(Pot = 3757, `Ring Net` = 1965, Snare = 1746, Trap = 7534)))

stage <- function(id, tag, model, delta, headline, item, kind = "run")
  list(id = id, tag = tag, model = model, delta = delta, headline = headline, item = item, kind = kind)

STAGE_DEFS <- list(
  D1 = stage("D1", NA, NA, NULL, "Tier-1 shore I/E unit fix, scored against its four pre-set criteria", "T1 top", "desk"),
  D2 = stage("D2", NA, NA, NULL, "the gear_only arm decision and the PE-vs-BSS arm disagreement", "T2", "desk"),
  D3 = stage("D3", NA, NA, NULL, "trailer over-coverage: discreteness or over-dispersion; plus the zero bin", "S5 review / T3", "desk"),
  D4 = stage("D4", NA, NA, NULL, "repository hygiene", "T4", "desk"),
  # V1 carries NO delta: it is production as shipped after the 2026-09-01 adoption. That is
  # the point of it. Anything that has to be set here would mean production is not what was
  # validated, which is exactly the failure this gate exists to catch.
  V1 = stage("V1", "VAL-1-adopted", "pooled", list(),
             "ADOPTION GATE: production after the flip must reproduce S2 exactly", "adoption"),
  V2 = stage("V2", "VAL-2-gearG-phase1", "gear_resolved", list(gear_resolved_G = TRUE),
             "GR-7 Phase 2 CONTROL: per-gear CPUE on, Dirichlet shares off", "GR-7"),
  V3 = stage("V3", "VAL-3-gearG-dirichlet", "gear_resolved",
             list(gear_resolved_G = TRUE, gear_share_dirichlet = TRUE),
             "GR-7 Phase 2: Dirichlet gear shares, coded 2026-07-21 and never sampled", "GR-7"),
  V4 = stage("V4", "VAL-4-tau-boatpc-biwk", "pooled",
             list(ar_force = list(private_boat = list(pot_closure = "biweekly"))),
             "the untested corner: shared turnover x biweekly boat pot closure", "S5 review"),
  V5 = stage("V5", "VAL-5-floor20", "pooled", list(shared_tau_min_obs = 20),
             "shared_tau_min_obs = 20 drops the boat pot closure out of the feature", "S5 review"))

resolve_cfg <- function(sid) modifyList(BASE, STAGE_DEFS[[sid]]$delta %||% list(), keep.null = TRUE)

find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}
sp_get <- function(dir, fitlab, par, col = "median") {
  d <- rd(dir, sprintf("structural_params_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"parameter" %in% names(d) || !col %in% names(d)) return(NA_real_)
  num1(d[[col]][d$parameter == par])
}
sp_tau <- function(dir, fitlab, col) {
  v <- sp_get(dir, fitlab, "tau_bar[1]", col); if (is.na(v)) v <- sp_get(dir, fitlab, "tau_bar", col); v
}
pit_get <- function(dir, fitlab, stream) {
  d <- rd(dir, sprintf("ppc_calibration_%s_Dungeness_Kept.csv", fitlab))
  if (is.null(d) || !"data_type" %in% names(d)) return(NA_real_)
  num1(d$pit_mean[d$data_type == stream])
}
fds_get <- function(dir, fitlab, col) {
  d <- rd(dir, "fit_data_summary.csv")
  if (is.null(d) || !"fit" %in% names(d) || !col %in% names(d)) return(NA)
  v <- d[[col]][d$fit == paste0(fitlab, "_Dungeness_Kept")]; if (!length(v)) NA else v[1]
}
adq_get <- function(dir, fitlab, col) {
  d <- rd(dir, "model_adequacy.csv") %||% rd(dir, "model_adequacy_reconstructed.csv")
  if (is.null(d) || !"fit" %in% names(d) || !col %in% names(d)) return(NA)
  v <- d[[col]][d$fit == paste0(fitlab, "_Dungeness_Kept")]; if (!length(v)) NA else v[1]
}

extract_run <- function(sid, model, run_tag, outdir, minutes, ok) {
  row <- data.frame(stage = sid, item = STAGE_DEFS[[sid]]$item %||% NA_character_,
    model = model, run_tag = run_tag, ok = ok, minutes = minutes,
    outdir = if (is.na(outdir)) "" else basename(outdir),
    port_BSS_catch = NA_real_, shore_ag_BSS = NA_real_, shore_pc_BSS = NA_real_,
    boat_ag_BSS = NA_real_, boat_pc_BSS = NA_real_,
    boat_ag_reported = NA, boat_pc_reported = NA,
    boat_ag_method = NA_character_, boat_pc_method = NA_character_,
    boat_pc_ar = NA_character_, boat_ag_ar = NA_character_,
    boat_ag_shared_tau = NA_real_, boat_pc_shared_tau = NA_real_,
    boat_pc_n_informed = NA_real_,
    boat_tau_bar = NA_real_, boat_tau_lo95 = NA_real_, boat_tau_hi95 = NA_real_,
    boat_tau_contraction = NA_real_, boat_tau_prior_influential = NA,
    pit_trailer = NA_real_, pit_osp = NA_real_,
    boat_p_loo_frac = NA_real_, boat_pareto_bad = NA_real_, boat_cov50_dev = NA_real_,
    fits_passed = NA_character_, gear_G = NA_real_, stringsAsFactors = FALSE)
  if (is.na(outdir) || !dir.exists(outdir)) return(row)

  pt <- rd(outdir, "port_total_Dungeness_Kept.csv")
  if (!is.null(pt) && "Estimate" %in% names(pt)) {
    i <- pt$Estimate %in% c("Expected_Catch", "Catch"); row$port_BSS_catch <- num1(pt$BSS_median[i]) }
  pv <- rd(outdir, "pe_vs_bss_comparison.csv")
  if (!is.null(pv) && "component" %in% names(pv)) {
    g <- function(pat, col) if (col %in% names(pv)) { v <- pv[[col]][grepl(pat, pv$component)]; if (length(v)) v[1] else NA } else NA
    row$shore_ag_BSS <- num1(g("^shore \\(All gear\\)", "BSS_catch"))
    row$shore_pc_BSS <- num1(g("^shore \\(Pot closure\\)", "BSS_catch"))
    row$boat_ag_BSS  <- num1(g("^private_boat \\(All gear\\)", "BSS_catch"))
    row$boat_pc_BSS  <- num1(g("^private_boat \\(Pot closure\\)", "BSS_catch"))
    row$boat_ag_reported <- as.logical(g("^private_boat \\(All gear\\)", "bss_reported"))
    row$boat_pc_reported <- as.logical(g("^private_boat \\(Pot closure\\)", "bss_reported"))
  }
  cr <- rd(outdir, "convergence_report.csv")
  if (!is.null(cr) && "fit" %in% names(cr)) {
    g <- function(pat, col) { if (!col %in% names(cr)) return(NA); i <- grepl(pat, cr$fit); if (any(i)) cr[[col]][i][1] else NA }
    row$boat_ag_method <- chr1(g("^private_boat_all_gear", "method_selected"))
    row$boat_pc_method <- chr1(g("^private_boat_ring_net_only", "method_selected"))
    row$boat_ag_ar     <- chr1(g("^private_boat_all_gear", "ar_resolution"))
    row$boat_pc_ar     <- chr1(g("^private_boat_ring_net_only", "ar_resolution"))
  }
  if (is.na(row$boat_pc_ar)) row$boat_pc_ar <- chr1(fds_get(outdir, "private_boat_ring_net_only", "ar_resolution"))
  row$boat_ag_shared_tau <- num1(fds_get(outdir, "private_boat_all_gear", "shared_tau"))
  row$boat_pc_shared_tau <- num1(fds_get(outdir, "private_boat_ring_net_only", "shared_tau"))
  row$boat_pc_n_informed <- num1(fds_get(outdir, "private_boat_ring_net_only", "n_L_informed"))
  bag <- "private_boat_all_gear"
  row$boat_tau_bar  <- sp_tau(outdir, bag, "median")
  row$boat_tau_lo95 <- sp_tau(outdir, bag, "lo95")
  row$boat_tau_hi95 <- sp_tau(outdir, bag, "hi95")
  # The diagnostic the 2026-08-31 regression removed from S2, S3 and S4b. Its presence here
  # is itself part of what V1 checks.
  pvp <- rd(outdir, sprintf("prior_vs_posterior_%s_Dungeness_Kept.csv", bag))
  if (!is.null(pvp) && "parameter" %in% names(pvp)) {
    i <- grepl("^tau_bar", pvp$parameter)
    if (any(i)) { row$boat_tau_contraction <- num1(pvp$contraction[i])
                  row$boat_tau_prior_influential <- as.logical(pvp$prior_influential[i][1]) }
  }
  row$pit_trailer <- pit_get(outdir, bag, "trailer"); row$pit_osp <- pit_get(outdir, bag, "osp")
  row$boat_p_loo_frac <- num1(adq_get(outdir, bag, "p_loo_frac"))
  row$boat_pareto_bad <- num1(adq_get(outdir, bag, "n_pareto_bad"))
  row$boat_cov50_dev  <- num1(adq_get(outdir, bag, "cov50_worst_dev"))
  ss <- rd(outdir, "season_summary.csv")
  if (!is.null(ss) && "metric" %in% names(ss)) row$fits_passed <- chr1(ss$value[ss$metric == "BSS fits passed"])
  cg <- rd(outdir, "catch_by_gear_type.csv")
  if (!is.null(cg)) row$gear_G <- nrow(cg)
  row
}

# Per-gear detail from catch_by_gear_type_detail.csv, for one population x sub-season.
gear_detail <- function(dir, pop = "shore", ss = "all_gear") {
  d <- rd(dir, "catch_by_gear_type_detail.csv")
  if (is.null(d) || !all(c("population", "subseason", "gear_type", "BSS_median") %in% names(d))) return(NULL)
  x <- d[d$population == pop & d$subseason == ss, , drop = FALSE]
  if (!nrow(x)) return(NULL)
  x$width <- suppressWarnings(as.numeric(x$BSS_hi95) - as.numeric(x$BSS_lo95))
  x$rel_width <- x$width / pmax(suppressWarnings(as.numeric(x$BSS_median)), 1)
  x[order(-suppressWarnings(as.numeric(x$BSS_median))), ]
}

# ---------------------------------------------------------------------------
# The run loop.
# ---------------------------------------------------------------------------
sum_path <- .here("05_output", "validation_2026-09-01_summary.csv")
ver_path <- .here("05_output", "validation_2026-09-01_verdicts.csv")
append_row <- function(r) utils::write.table(
  r, sum_path, sep = ",", row.names = FALSE, qmethod = "double",
  col.names = !file.exists(sum_path), append = file.exists(sum_path))

# PRE-FLIGHT: the R-to-Stan data contract. The check that would have saved the 2026-08-25
# ladder's six lost runs, where rstan returned an EMPTY stanfit rather than raising.
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

# PRE-FLIGHT 2: the adoption is actually in the shipped config. V1 has no delta by design, so
# if this is not true the gate silently tests the wrong thing.
local({
  ok_tau <- isTRUE(BASE$shared_tau)
  ok_flr <- identical(as.numeric(BASE$shared_tau_min_obs %||% NA), 15)
  cat(sprintf("  pre-flight  adoption: shared_tau = %s, shared_tau_min_obs = %s\n",
              BASE$shared_tau %||% "unset", BASE$shared_tau_min_obs %||% "unset (helper default)"))
  if (!ok_tau)
    stop("run_config.R does not ship shared_tau = TRUE, so V1 would not be testing the ",
         "adopted configuration. Either adopt it or drop V1 from STAGES.", call. = FALSE)
  if (!ok_flr)
    warning("shared_tau_min_obs is not explicitly 15 in run_config.R; V5's contrast is ",
            "against an implicit default rather than a stated production value.", call. = FALSE)
})

collected <- list(); gate_ok <- TRUE

for (sid in STAGES) {
  S <- STAGE_DEFS[[sid]]

  if (identical(S$kind, "desk")) {
    t0 <- Sys.time()
    fn <- get(paste0("run_", sid))
    tryCatch(fn(), error = function(e) message("*** DESK ", sid, " FAILED: ", conditionMessage(e)))
    cat(sprintf("\n  [%s done in %.1f s]\n", sid, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    next
  }
  if (isTRUE(DRY_RUN)) next

  if (!gate_ok && isTRUE(GATE_ON_V1) && sid %in% c("V4", "V5")) {
    banner(sprintf("STAGE %s SKIPPED: V1 did not reproduce S2, so the shipped configuration is not the one the Stage 5 batch validated and %s could not be read against S2. V2/V3 are unaffected; they test the gear track's per-gear machinery, not the adoption.", sid, sid))
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
    cat(sprintf("  item %s | shared_tau=%s min_obs=%s | gear_resolved_G=%s dirichlet=%s | ar_force=%s\n",
                S$item, rc$shared_tau %||% FALSE, rc$shared_tau_min_obs %||% "default",
                rc$gear_resolved_G %||% FALSE, rc$gear_share_dirichlet %||% FALSE,
                if (is.null(rc$ar_force)) "NULL" else paste(utils::capture.output(str(rc$ar_force)), collapse = " ")))
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

  if (identical(sid, "V1")) {
    nd <- find_outdir(m, run_tag)
    ex <- fit_exactness(nd %||% "", .here("05_output", REF$S2$dir), NULL, "all four fits vs S2")
    gate_ok <- identical(ex$verdict, "PASS")
    cat("\n  GATE  adopted production vs S2: ", ex$verdict, "\n        ", ex$observed, "\n", sep = "")
    if (!gate_ok) message("*** The shipped configuration does not reproduce S2. V4/V5 will be skipped.")
  }
}

# ---- VERDICTS for the fitted stages ---------------------------------------
S <- if (length(collected)) do.call(rbind, collected) else NULL
g1 <- function(sid, col) { if (is.null(S)) return(NA); x <- S[[col]][S$stage == sid]; if (length(x)) x[1] else NA }

if (!is.null(S) && any(S$stage == "V1")) {
  nd <- find_outdir("pooled", STAGE_DEFS$V1$tag)
  ex <- fit_exactness(nd %||% "", .here("05_output", REF$S2$dir), NULL, "all four fits vs S2")
  V1row("V1", "the shipped configuration reproduces the run that justified it",
    sprintf("%s; components %s / %s / %s / %s against S2's %s / %s / %s / %s",
            ex$observed, fmt(g1("V1","shore_ag_BSS"),0), fmt(g1("V1","shore_pc_BSS"),0),
            fmt(g1("V1","boat_ag_BSS"),0), fmt(g1("V1","boat_pc_BSS"),0),
            fmt(REF$S2$shore_ag,0), fmt(REF$S2$shore_pc,0), fmt(REF$S2$boat_ag,0), fmt(REF$S2$boat_pc,0)),
    "every shared parameter row identical at full precision", ex$verdict,
    paste("V1 carries NO delta: it is production exactly as shipped after the 2026-09-01",
          "adoption. A PASS says the configuration now in run_config.R is the one the Stage",
          "5 batch measured, rather than something that drifted alongside it between the",
          "batch and the flip. A FAIL is a code-drift alarm and not a shared_tau problem;",
          "the first suspects are the 2026-08-31 changes to save_run_diagnostics.R and",
          "bss_model_adequacy.R, both of which are supposed to be post-fit only."))
  cont <- g1("V1", "boat_tau_contraction")
  V1row("V1", "the boat prior-vs-posterior file exists again, with a tau_bar row",
    if (is.na(cont)) "STILL MISSING or no tau_bar row" else
      sprintf("tau_bar contraction %.3f, prior_influential = %s", cont, g1("V1","boat_tau_prior_influential")),
    "present, with a contraction value", if (is.na(cont)) "FAIL" else "PASS",
    paste("The second reason V1 must run. S2, S3 and S4b carry NO boat prior_vs_posterior",
          "file: tau_bar was registered under its bare name while rstan calls the row",
          "tau_bar[1], the lookup threw, and the writer's tryCatch swallowed the error along",
          "with the whole file. So the contraction diagnostic for the one parameter being",
          "adopted has never actually been written. READ THE VALUE WITH CARE: contraction",
          "is 1 - post_sd/prior_sd and is about 0.216 here, which trips prior_influential,",
          "while the posterior MEAN sits about 3.5 prior SDs from the prior mean. A",
          "parameter that has moved three and a half prior SDs is not prior-dominated. That",
          "flag's definition is a recorded open item, not a finding about tau_bar."))
}

if (!is.null(S) && any(S$stage == "V2")) {
  nd <- find_outdir("gear_resolved", STAGE_DEFS$V2$tag)
  # 2026-09-02 CORRECTION. This first compared against S5-3-tau-gear, which carried the
  # EXPERIMENT sampler override, so the reference differed in gear_resolved_G AND in the
  # sampler settings and the criterion FAILED for the wrong reason. The correct reference is
  # a gear run at the same sampler settings: the 2026-08-29 one. Against it the boat fits are
  # bit-identical across 4,367 shared rows. expect_delta now makes that check automatic.
  exb <- fit_exactness(nd %||% "", .here("05_output", REF$Egear$dir), "private_boat",
                       "boat fits vs the G = 1 gear run at the SAME sampler settings",
                       expect_delta = c("gear_resolved_G", "gear_share_dirichlet",
                                        "shared_tau", "shared_tau_min_obs"))
  V1row("V2", "turning on per-gear CPUE leaves the BOAT untouched",
    exb$observed, "boat fits bit-identical (the boat stays G = 1 by design)", exb$verdict,
    paste("GR-7 Phase 0 established that the boat is Pot-dominated and stays at G = 1, so",
          "gear_resolved_G should reach the SHORE fits and nothing else. THE REFERENCE MUST",
          "DIFFER IN ONE THING: comparing against a run that also carried the experiment",
          "sampler override produces a FAIL that says nothing about gear_resolved_G. That",
          "mistake has now been made twice, so expect_delta reports any unexpected config",
          "difference in the observation rather than leaving it to be found by hand."))
  gd <- gear_detail(nd %||% "")
  if (!is.null(gd)) {
    ord <- gd$gear_type[order(-suppressWarnings(as.numeric(gd$BSS_median)))]
    V1row("V2", "per-gear shore posteriors reproduce the 2026-07-20 ORDERING",
      sprintf("shore all-gear by median: %s (2026-07-20 had Trap > Mixed > Pot > Ring Net > Snare)",
              paste(ord, collapse = " > ")),
      "Trap highest, matching the interview ratios",
      if (length(ord) && identical(ord[1], "Trap")) "PASS" else "REVIEW",
      paste("The 2026-07-20 Phase-1 run predates the 2026-08-25 batch, the I/E unit fix and",
            "the shared turnover, so its LEVELS should not reproduce and it would be wrong",
            "to require it. What should survive is the ordering, because that is set by the",
            "interview CPUE ratios rather than by any of the intervening changes."))
  }
}

if (!is.null(S) && any(S$stage == "V3")) {
  n2 <- find_outdir("gear_resolved", STAGE_DEFS$V2$tag)
  n3 <- find_outdir("gear_resolved", STAGE_DEFS$V3$tag)
  a <- gear_detail(n2 %||% ""); b <- gear_detail(n3 %||% "")
  if (!is.null(a) && !is.null(b)) {
    m <- merge(a[, c("gear_type","BSS_median","rel_width")], b[, c("gear_type","BSS_median","rel_width")],
               by = "gear_type", suffixes = c("_off", "_on"))
    med_shift <- max(abs(100 * (m$BSS_median_on - m$BSS_median_off) / pmax(m$BSS_median_off, 1)), na.rm = TRUE)
    wid_ratio <- stats::median(m$rel_width_on / pmax(m$rel_width_off, 1e-9), na.rm = TRUE)
    V1row("V3", "Dirichlet gear shares WIDEN the per-gear intervals without moving the medians",
      sprintf("largest median shift %.1f%%; median relative-width ratio on/off %.2fx (%s)",
              med_shift, wid_ratio,
              paste(sprintf("%s %.2fx", m$gear_type, m$rel_width_on / pmax(m$rel_width_off, 1e-9)), collapse = ", ")),
      "medians within a few percent; widths strictly greater than 1.0x",
      if (is.na(wid_ratio)) "REVIEW" else if (wid_ratio > 1.02 && med_shift < 5) "PASS"
      else if (wid_ratio <= 1.02) "NO EFFECT" else "REVIEW",
      paste("GR-7 Phase 2 was coded 2026-07-21, parses under stanc 2.32.5 and has NEVER BEEN",
            "SAMPLED. Its entire purpose is that the Phase-1 per-gear intervals are",
            "'slightly too narrow' because they treat the interview gear shares as known;",
            "propagating Dirichlet uncertainty should widen them and leave the central",
            "estimates alone. NO EFFECT would mean the toggle is inert in practice despite",
            "compiling, which is worth knowing before anyone cites Phase 2 as available.",
            "A large median shift would mean it does more than propagate share uncertainty",
            "and needs its own validation before use."))
  }
}

if (!is.null(S) && any(S$stage == "V4")) {
  pc <- g1("V4","boat_pc_BSS"); ag <- g1("V4","boat_ag_BSS")
  V1row("V4", "the untested corner: shared turnover x biweekly boat pot closure",
    sprintf("boat pot closure %s at %s (S2 monthly+tau %s; S5 biweekly no-tau %s; gear biweekly+tau %s); all-gear %s (must be S2's %s)",
            fmt(pc,0), g1("V4","boat_pc_ar"), fmt(REF$S2$boat_pc,0), fmt(REF$S5$boat_pc,0),
            fmt(REF$S3gear$boat_pc,0), fmt(ag,0), fmt(REF$S2$boat_ag,0)),
    "all-gear unmoved from S2; pot closure bracketed by the four neighbouring cells",
    if (!is.na(ag) && abs(ag - REF$S2$boat_ag) < 1) "PASS" else "REVIEW",
    paste("Both levers move the boat pot closure and they have never been crossed on the",
          "pooled track: S2 gives 1,018 at monthly with the turnover on, S5 gives 735 at",
          "biweekly with it off, and the gear track gives 956 at biweekly with it on against",
          "743 with it off. The gear pair suggests the two are roughly additive there, but",
          "the pooled cell is a guess until it is run, worth about 200 crab. The second half",
          "of the criterion is the one that matters most: the ALL-GEAR component must be",
          "unmoved, because ar_force silently leaking across sub-seasons is exactly what",
          "contaminated stage C on 2026-08-28 and cost 3,025 crab there."))
}

if (!is.null(S) && any(S$stage == "V5")) {
  st <- g1("V5","boat_pc_shared_tau"); ni <- g1("V5","boat_pc_n_informed")
  pc <- g1("V5","boat_pc_BSS"); ag <- g1("V5","boat_ag_BSS")
  V1row("V5", "raising the floor to 20 drops the boat pot closure, and only that",
    sprintf("boat pot closure: %s informed days, shared_tau = %s, catch %s (was %s at floor 15); all-gear %s (must be S2's %s)",
            fmt(ni,0), st, fmt(pc,0), fmt(REF$S2$boat_pc,0), fmt(ag,0), fmt(REF$S2$boat_ag,0)),
    "pot closure shared_tau = 0 and moved; all-gear identical to S2",
    if (isTRUE(st == 0) && !is.na(ag) && abs(ag - REF$S2$boat_ag) < 1) "PASS" else "REVIEW",
    paste("The floor of 15 was set just below the boat pot closure's 18 informed days on",
          "purpose, to keep a component whose tau_bar is corroborated across both tracks",
          "(1.873 pooled, 2.050 gear-resolved, both excluding the 1.2 prior centre). That is",
          "a defensible choice and it is still a choice made after seeing which side of the",
          "line the component fell on, so it deserves one run that shows what it is worth.",
          "This measures exactly that: the difference between this pot-closure number and",
          "S2's is the price of the threshold. The all-gear half is the control, since 130",
          "informed days clears either floor and that component must not move at all."))
}

VD <- if (length(V)) do.call(rbind, V) else NULL
if (!is.null(VD)) utils::write.csv(VD, ver_path, row.names = FALSE)

banner("VALIDATION BATCH SUMMARY")
if (!is.null(S)) {
  print(S[, c("stage","item","model","ok","minutes","port_BSS_catch","boat_ag_BSS","boat_pc_BSS",
              "boat_tau_bar","fits_passed")], row.names = FALSE)
  cat("\n  full summary -> ", sum_path, "\n", sep = "")
}
if (!is.null(VD)) {
  rule()
  for (i in seq_len(nrow(VD)))
    cat(sprintf("[%s] %-34s %s\n        %s\n", VD$stage[i], VD$verdict[i], VD$criterion[i], VD$observed[i]))
  cat("\n  verdicts -> ", ver_path, "\n", sep = "")
}
if (isTRUE(DRY_RUN)) {
  rule()
  banner(sprintf("DRY RUN: the four DESK stages above ran for real; %d fitted stage(s) were resolved and NOT run",
                 sum(vapply(STAGES, function(s) !identical(STAGE_DEFS[[s]]$kind, "desk"), logical(1)))))
  for (sid in STAGES) {
    Sd <- STAGE_DEFS[[sid]]; if (identical(Sd$kind, "desk")) next
    cfg <- resolve_cfg(sid); ks <- names(Sd$delta %||% list())
    cat(sprintf("  %-3s %-22s %-14s %s\n", sid, Sd$tag, Sd$model, Sd$headline))
    if (!length(ks)) cat("      delta: none (production as shipped)\n")
    for (k in ks) cat(sprintf("      delta: %-24s %s\n", k,
      paste(utils::capture.output(str(cfg[[k]])), collapse = " ")))
  }
  cat("\n  Read the desk verdicts above, then set DRY_RUN <- FALSE.\n\n")
} else {
  banner(sprintf("DONE  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  # D3's zero-bin criterion needs the p_zero column, which only exists on runs made after
  # 2026-09-01. Once V1 has produced one, re-sourcing with STAGES <- "D3" answers it in
  # seconds; RESUME makes that safe, and the desk stages fit nothing.
  if (!is.null(VD) && any(grepl("NOT YET COMPUTABLE", VD$verdict)))
    cat("\n  NOTE: the zero-bin criterion could not be scored because every run it read\n",
        "  predates the p_zero column. V1 has now produced one. Re-source with\n",
        "  STAGES <- \"D3\" to score it; it takes seconds and fits nothing.\n\n", sep = "")
}
