###############################################################################
# BATCH: THE LADDER THAT DID NOT LADDER, AND THE ZINB DECISION RE-SCORED
# 2026-09-04. Follow-up to run_shore_ar_zi_2026-09-03.R.
# -----------------------------------------------------------------------------
# WHY THIS RUN EXISTS
#   The 2026-09-03 batch spent 23.6 h and returned five of eleven verdicts wrong. Three
#   one-line defects caused all five:
#
#   D1  THE LADDER NEVER LADDERED. The driver gated force_resolution on
#       isTRUE(params$ar_escalate). ar_escalate had become SCOPED that same day and the
#       stage passed list(shore = "all_gear"); isTRUE() on a list is FALSE, so
#       force_resolution stayed NULL while the ladder itself was correctly built with
#       four rungs and the loop correctly ran four attempts. Every attempt therefore
#       re-derived `daily` from the data-driven selector: four identical fits, P_n = 289,
#       catch_median 20,897.819 to three decimals on all four, 14.7 h of the 23.6 h batch.
#       Both A1-A3 verdicts are void; the question they were asked is still open.
#
#   D2  THE ZINB PPC WAS COMPUTED UNDER NB2. pit_block() and calib() scored the catch
#       stream with dnbinom/pnbinom while the fit was a mixture, so theta_C was excluded
#       from the expected zero count: 419.0 expected against 676 observed, z = +15.4, and
#       the batch recorded that the feature made the zero bin WORSE. Under the mixture it
#       is 637.9, z = +2.0, against the NB2 baseline's 605.4 / z = +3.8. The same defect
#       makes Z1's catch pit_mean, pit_sd, coverage_50 and flag_pit_bias unreadable.
#
#   D3  THE elpd CRITERION USED THE WRONG SE. +14.8 nats was compared against
#       se_elpd_loo (about 46 nats), which is the SE of ONE model's total and is dominated
#       by across-observation variation common to both models. The paired-difference SE is
#       5.5, so the gain is 2.69 SE and clears the stated 2 SE bar. NOT WORTH IT was wrong.
#
#   A fourth defect explains the Z0 caveat: run_parameters.txt was truncated by str()'s
#   list.len = 99 cap while run_config carries 120+ keys, so the "9 UNEXPECTED config
#   keys" were keys that crossed the cut, not a config difference.
#
#   All four are fixed. This batch re-runs only what the fixes actually change.
#
# WHAT THIS RUN DOES NOT DO
#   It does not re-litigate G1 or P1. Those two verdicts were computed correctly and are
#   unaffected by any of the four defects: the gear boat all-gear fit samples (2
#   divergences against 554, port 70,953 against 51,385), and the PE arm alignment moved
#   the boat by -0.059%. Nothing here touches them.
#
# RUNTIME, MEASURED. L1 took 5.4 h and Z2 3.7 h on the 2026-09-04 run, 9.1 h in total; the
# ladder rungs cost 95.6 / 91.1 / 83.3 min, so coarsening the AR barely reduces the cost and
# a ladder should be budgeted at one full fit per rung. With L1 and Z2 resumed from disk,
# the outstanding work (C1 + C2) is about 4 h + 5 h. The original estimate below was wrong.
#
# ORIGINAL ESTIMATE (kept because it was wrong by 2x and the reason is worth remembering):
# About 6-7 h with both switches at their defaults, against 23.6 h for the run it
# corrects: L1 about 2-3 h (the pooled shore all-gear fit took ~205 min at daily on this
# machine and the gear track fits it in 10 min at monthly, so weekly/biweekly/monthly
# should come in around 40-60, 25-40 and 15-25 min, plus ~60 min for the three untouched
# components), Z2 about 4 h (Z1 took 242 min). LADDER_INCLUDE_DAILY adds ~3.5 h;
# ZINB_RERENDER = FALSE removes ~4 h.
#
# NO STAN RECOMPILE unless ZINB_RERENDER is TRUE and the cached ZINB binary is gone; the
# Stan file is unchanged from 2026-09-03.
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                     # TRUE: desk stages run, nothing is fitted. START HERE.
# 2026-09-06: L1 and Z2 have RUN (results at 20260903/pooled-CPUE-LZ-*). Their stage
# definitions are kept so RESUME re-scores them without refitting. The work still to do is
# C1 (the candidate production configuration: weekly AR + ZINB together, never yet run) and
# C2 (the ladder again with per-rung adequacy, which the first ladder could not record).
# To re-score only, leave this as c("D0", "L1", "Z2") with RESUME = TRUE.
STAGES  <- c("D0", "L1", "Z2", "C1", "C2")
RESUME  <- TRUE

# ---- SWITCH 1: does the ladder refit its DAILY rung? -----------------------
# FALSE (default) runs weekly -> biweekly -> monthly and takes the daily rung from the
#   committed 2026-09-01 Z0 folder, whose config differs from L1's only in the ar_escalate
#   keys, none of which touch the model, the data or the seed. The daily rung is therefore
#   already known exactly: catch 20,897.819, 333 divergences, p_loo 35.2% of n_obs, 41
#   Pareto k above 0.7, coverage_50 0.701. Saves about 3.7 h.
# TRUE refits it, so all four rungs live in one folder with one config dump and one HTML
#   report, and the comparison needs no cross-folder claim. Costs about 3.7 h to reproduce
#   a number we hold. Choose TRUE if this ladder is going to be shown to anyone outside
#   the project, because "the daily rung is in a different folder" is a question a reviewer
#   will ask and a cross-folder answer is weaker than a same-folder one.
LADDER_INCLUDE_DAILY <- FALSE

# ---- SWITCH 2: re-render the ZINB shore fit with the corrected PPC? --------
# TRUE (default) refits Z1 under the corrected diagnostics, about 4 h. The elpd evidence
#   is already correct without it (Stan's log_lik always carried the mixture), but the
#   zero bin and the whole calibration block on disk are not, and the corrected numbers in
#   the 2026-09-03 verdicts file were recomputed OFFLINE from committed summaries using a
#   mean-product approximation, E[theta * p0] taken as E[theta] * E[p0]. That
#   approximation is very unlikely to move z by more than about 0.1, but "very unlikely"
#   is not a basis for adopting a change to the likelihood of a harvest estimate.
# FALSE decides ZINB on elpd alone and leaves ppc_* and model_adequacy wrong on disk for
#   that folder. Only choose FALSE if the ladder result is urgent and ZINB can wait.
ZINB_RERENDER <- TRUE

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
# 2026-09-06: VECTORISED. The scalar form used `if (is.na(x))`, which under R >= 4.2 is
# an ERROR on a length-3 input, not a warning; verdict_L1 passes whole ladder columns to
# it and the batch aborted after both stages had already been fitted. A formatter that
# only works on scalars has no business in a verdict block that summarises a table.
fmt <- function(x, d = 1) {
  if (length(x) == 0) return("NA")
  out <- formatC(suppressWarnings(as.numeric(x)), format = "f", digits = d, big.mark = ",")
  out[is.na(suppressWarnings(as.numeric(x)))] <- "NA"
  out
}

if (!isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages({ library(here); library(rmarkdown) })
  load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
  install.lib <- load.lib[!load.lib %in% installed.packages()]
  for (lib in install.lib) install.packages(lib, dependencies = TRUE)
  invisible(sapply(load.lib, require, character.only = TRUE))
  rstan_options(auto_write = TRUE)
}
invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE),
                 function(f) try(source(f), silent = TRUE)))
source(.here("run_config.R"))
BASE <- run_config

model_rmd <- list(pooled = .here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"))
prefix    <- list(pooled = "pooled-CPUE-")
stopifnot(all(file.exists(unlist(model_rmd))))
desk_dir  <- .here("05_output", format(Sys.Date(), "%Y%m%d"), "ladder-zinb-desk")

rd <- function(dir, f) {
  p <- file.path(dir, f); if (!file.exists(p)) return(NULL)
  tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}
V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = observed,
                                    threshold = threshold, verdict = verdict, why = why,
                                    stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# BASELINES, each named beside the folder it comes from.
# ---------------------------------------------------------------------------
REF <- list(
  # The daily rung. Stage Z0 of the 2026-09-03 batch: the ZINB Stan edit with the feature
  # OFF, proven bit-identical to 2026-08-31 production across 11,223 parameter rows.
  DAILY = list(dir = "20260901/pooled-CPUE-SZ-Z0-zi-off",
               shore_ag = 20897.819, lo95 = 18186.53, hi95 = 24402.44, pi_rel = 0.2974,
               divergences = 333, p_loo_frac = 0.3519, bad_k = 41, cov50_catch = 0.4572,
               port = 71450),
  # The other end of the bracket: the gear track fits this component at MONTHLY.
  GEARMON = list(p_loo_frac = 0.036, bad_k = 0, cov50 = 0.035),
  # The NB2 baseline for the ZINB comparison, same folder as DAILY.
  Z0 = list(dir = "20260901/pooled-CPUE-SZ-Z0-zi-off"),
  Z1 = list(dir = "20260902/pooled-CPUE-SZ-Z1-zi-shore",
            theta_ag = 0.1780, theta_pc = 0.1072,
            elpd_diff_ag = 14.8, se_paired_ag = 5.5, elpd_diff_pc = 9.1, se_paired_pc = 4.6,
            zero_mix_ag = 637.9, zero_z_ag = 2.0, zero_mix_pc = 138.1, zero_z_pc = 0.8))

LADDER <- if (isTRUE(LADDER_INCLUDE_DAILY)) c("daily", "weekly", "biweekly", "monthly") else
                                            c("weekly", "biweekly", "monthly")

STAGE_DEFS <- list(
  D0 = list(id = "D0", kind = "desk", model = NA, tag = NA, delta = NULL,
            headline = "desk: re-score the 2026-09-03 batch from committed files, no fitting"),
  # L1 is the A1 stage the last batch was supposed to run. Same production mechanism,
  # ar_escalate scoped to shore all-gear with all_rungs, now that force_resolution honours
  # the scoped form. ar_escalate_ladder drops the daily rung unless SWITCH 1 says otherwise.
  L1 = list(id = "L1", kind = "run", model = "pooled", tag = "LZ-L1-shore-ladder",
            delta = list(ar_escalate = list(shore = "all_gear"),
                         ar_escalate_stop = "all_rungs",
                         ar_escalate_select = "first_pass",
                         ar_escalate_respect_cap = FALSE,
                         ar_escalate_ladder = LADDER,
                         ar_escalate_max_attempts = length(LADDER),
                         estimate_catch_zi = FALSE),
            headline = "shore all-gear at every rung, for real this time"),
  Z2 = list(id = "Z2", kind = "run", model = "pooled", tag = "LZ-Z2-zi-shore",
            delta = list(estimate_catch_zi = TRUE, catch_zi_populations = c("shore")),
            headline = "ZINB shore catch, re-rendered under the CORRECTED PPC"),
  # ---- added 2026-09-06 after the L1/Z2 review -------------------------------
  # C1 is the CANDIDATE PRODUCTION CONFIGURATION. L1 moved the shore effort process to
  # weekly and Z2 changed the shore catch likelihood, in separate runs, and both look
  # right; neither has been run in the other's presence. Zero-inflation is a catch-stream
  # change and leaves the daily-AR overfitting untouched (Z2 still carries p_loo at 34.4%
  # of n_obs and 26 bad Pareto k), so the two are complementary rather than alternatives.
  # Read C1 against L1 to isolate the LIKELIHOOD at the new resolution, and against Z2 to
  # isolate the RESOLUTION under the new likelihood. The boat is untouched in both, which
  # makes it a free negative control for the third time.
  #
  # ar_force, not ar_escalate: this is not a ladder, it is one fit at a chosen resolution,
  # and ar_force is the lever that bypasses the data-driven selector and the cap.
  C1 = list(id = "C1", kind = "run", model = "pooled", tag = "LZ-C1-weekly-zi",
            delta = list(ar_force = list(shore = list(all_gear = "weekly")),
                         estimate_catch_zi = TRUE, catch_zi_populations = c("shore")),
            headline = "CANDIDATE: shore all-gear at WEEKLY with the ZINB catch likelihood"),
  # C2 is the missing half of the bracket. The 2026-09-04 ladder left biweekly and monthly
  # with no p_loo, no Pareto count and no coverage, so "weekly" currently means "the finest
  # rung that is clearly not overfitted" rather than "the best rung". With
  # ar_rung_adequacy = TRUE every rung reports them. Runs the SAME three coarse rungs; the
  # daily rung stays out for the same reason as before, and its adequacy is already known.
  C2 = list(id = "C2", kind = "run", model = "pooled", tag = "LZ-C2-ladder-adequacy",
            delta = list(ar_escalate = list(shore = "all_gear"),
                         ar_escalate_stop = "all_rungs",
                         ar_escalate_select = "first_pass",
                         ar_escalate_respect_cap = FALSE,
                         ar_escalate_ladder = c("weekly", "biweekly", "monthly"),
                         ar_escalate_max_attempts = 3L,
                         ar_rung_adequacy = TRUE,
                         estimate_catch_zi = FALSE),
            headline = "the ladder again, with adequacy recorded for EVERY rung"))
if (!isTRUE(ZINB_RERENDER)) STAGES <- setdiff(STAGES, "Z2")

resolve_cfg <- function(sid) modifyList(BASE, STAGE_DEFS[[sid]]$delta %||% list(), keep.null = TRUE)
find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}

# ---------------------------------------------------------------------------
# PRE-FLIGHT. Each check corresponds to one of the four defects, and each would have
# caught it before the run rather than after. A pre-flight that only restates the config
# is decoration; these read the SOURCE and assert the property that failed.
# ---------------------------------------------------------------------------
preflight <- function() {
  banner("PRE-FLIGHT")
  fails <- character(0)
  say <- function(ok, msg, detail = "") {
    cat(sprintf("  [%s] %s%s\n", if (ok) "ok" else "XX", msg,
                if (nzchar(detail)) paste0("  --  ", detail) else ""))
    if (!ok) fails <<- c(fails, msg)
  }
  # D1. No driver may test ar_escalate with isTRUE(); a scoped value reads as FALSE.
  for (drv in list.files(.here("01_BSS_models"), pattern = "\\.Rmd$", full.names = TRUE)) {
    src  <- readLines(drv, warn = FALSE)
    code <- src[!grepl("^\\s*#", src)]
    say(!any(grepl("isTRUE(params$ar_escalate)", code, fixed = TRUE)),
        sprintf("D1 %s: ar_escalate is resolved, never isTRUE()", basename(drv)))
    say(any(grepl("force_res <- if (.esc_on)", code, fixed = TRUE)),
        sprintf("D1 %s: force_resolution honours the resolved flag", basename(drv)))
    say(any(grepl(".kept_attempt[[label]] <- attempt", code, fixed = TRUE)),
        sprintf("D1 %s: the reported rung is tracked by attempt index", basename(drv)))
    say(any(grepl("str(params, list.len =", src, fixed = TRUE)),
        sprintf("D4 %s: config dump is untruncated", basename(drv)))
  }
  # D1 again, at the level that matters: this stage's OWN ladder must have distinct rungs.
  cfg <- resolve_cfg("L1")
  esc <- .bss_resolve_ar_escalate(cfg, "shore", "all_gear")
  say(isTRUE(esc), "D1 L1: ar_escalate resolves TRUE for shore/all_gear")
  say(!isTRUE(.bss_resolve_ar_escalate(cfg, "shore", "pot_closure")) &&
      !isTRUE(.bss_resolve_ar_escalate(cfg, "private_boat", "all_gear")),
      "D1 L1: the scope reaches shore all-gear and nothing else")
  say(!any(duplicated(LADDER)) && length(LADDER) >= 3,
      "D1 L1: the ladder carries distinct rungs", paste(LADDER, collapse = " -> "))
  # D2. The PPC must route through the shared mixture helper on the catch stream.
  for (f in c("save_run_diagnostics.R", "model_diagnostics.R")) {
    src <- readLines(.here("03_R_functions", f), warn = FALSE)
    say(any(grepl("bss_zi_pit(", src, fixed = TRUE)) &&
        any(grepl("bss_zi_theta_draws(fit, stan_data,", src, fixed = TRUE)),
        sprintf("D2 %s: catch PPC uses the mixture likelihood", f))
  }
  say(abs(bss_zi_p_zero(3.0, 0.8, 0.178) -
          (0.178 + 0.822 * dnbinom(0, size = 0.8, mu = 3.0))) < 1e-12,
      "D2: mixture p_zero arithmetic is correct")
  say(abs(bss_zi_pit(0, 3.0, 0.8, 0.178) -
          0.5 * (0.178 + 0.822 * dnbinom(0, size = 0.8, mu = 3.0))) < 1e-12,
      "D2: mixture PIT at y = 0 uses F(-1) = 0")
  # D3. The paired-SE helper is present and behaves.
  say(exists("loo_elpd_paired", mode = "function"), "D3: paired elpd helper is sourced")
  # Draw persistence, so the next diagnostic fix is not another re-fit.
  say(isTRUE(BASE$save_ppc_draws), "PPC draws will be persisted for every fit")
  # Folder naming. The driver reads run_config$run_tag and nothing else; simulate exactly
  # what run_stage() will hand it and assert the two fitted stages cannot collide.
  tags <- vapply(setdiff(STAGES, "D0"), function(sid) {
    cfg <- resolve_cfg(sid); cfg$run_tag <- STAGE_DEFS[[sid]]$tag
    if (is.null(cfg$run_tag) || !nzchar(cfg$run_tag)) "<blank>" else cfg$run_tag }, character(1))
  say(all(tags != "<blank>") && !any(duplicated(tags)) &&
      !any(tags == (BASE$run_tag %||% "")),
      "each fitted stage writes to its OWN folder, none to the run_config default",
      paste(names(tags), tags, sep = " -> ", collapse = "; "))
  .self <- .here("06_diagnostics", "run_ladder_zinb_2026-09-04.R")
  .src <- if (file.exists(.self)) readLines(.self, warn = FALSE) else character(0)
  .src <- .src[!grepl("^\\s*#", .src)]
  # Match the actual render CALL (rmarkdown::render(...)), which excludes this check's own
  # text; a literal-string self-match is how the first version of this line tripped itself.
  .calls <- .src[grepl("rmarkdown::render(", .src, fixed = TRUE)]
  say(length(.calls) > 0 && !any(grepl("output_dir", .calls, fixed = TRUE)) &&
      any(grepl("cfg$run_tag <- st$tag", .src, fixed = TRUE)),
      "run_stage() puts the tag INSIDE run_config and does not pass output_dir= to render()",
      paste(trimws(.calls), collapse = " ; "))
  say(any(grepl("file.copy(html, file.path(od, basename(html))", .src, fixed = TRUE)),
      "run_stage() MOVES the rendered HTML into the run folder (else the next stage overwrites it)")
  rule()
  if (length(fails)) {
    cat("  PRE-FLIGHT FAILED:\n"); for (f in fails) cat("   -", f, "\n")
    stop("Fix the above before running. Every one of these is a defect the 2026-09-03 batch shipped.")
  }
  cat("  pre-flight clean.\n")
}

# ---------------------------------------------------------------------------
# D0. Re-score the 2026-09-03 batch from committed files. Costs nothing and runs in a
# dry run, so the corrected numbers are reproducible by anyone with the repo.
# ---------------------------------------------------------------------------
stage_D0 <- function() {
  banner("D0  desk re-score of the 2026-09-03 batch")
  dir.create(desk_dir, recursive = TRUE, showWarnings = FALSE)
  z0 <- .here("05_output", REF$Z0$dir); z1 <- .here("05_output", REF$Z1$dir)

  # --- D3: the elpd decision, paired ---------------------------------------
  for (nm in c("shore_all_gear", "shore_ring_net_only")) {
    f <- sprintf("loo_pointwise_catch_%s_Dungeness_Kept.csv", nm)
    r <- loo_elpd_paired(file.path(z0, f), file.path(z1, f), nm)
    if (is.null(r)) next
    V1row("D0", sprintf("ZINB elpd on %s, split by COUNT SIZE", nm),
          loo_elpd_by_count_str(r), "read as a whole",
          if (isTRUE(r$by_count["1", "diff"] < -3 * r$by_count["1", "se"])) "MISFIT MIGRATED" else "READ",
          paste("2026-09-05. The two-way zero/positive split said 'positives -10.5' and that was read as",
                "the harvest-carrying counts getting worse. It is not: y=1 loses -42.0 at 16.7 SE, y=2 loses",
                "-6.6, and every bin from 3 up IMPROVES, with the 3+ bins carrying 87% of the catch. r_C",
                "doubles (0.95 -> 1.83): once theta absorbs structural zeros the NB2 no longer needs extreme",
                "overdispersion to reach zero, so it tightens, fitting the mid and upper counts better and",
                "the almost-zero count of 1 worse. The ZINB has MOVED the misfit from the 0 bin to the 1",
                "bin. Excess mass at 0 AND 1 relative to NB2 is the signature of a two-regime process",
                "(unsuccessful trips yielding 0-1 crab against successful ones), which a hurdle or a",
                "two-component mixture fits and a ZINB cannot."))
    V1row("D0", sprintf("ZINB elpd on %s, scored against the PAIRED difference SE", nm),
          loo_elpd_paired_str(r), "at least 2 SE",
          if (isTRUE(r$ratio >= 2)) "WORTH IT" else if (isTRUE(r$ratio >= 1)) "MARGINAL" else "NOT WORTH IT",
          paste("The 2026-09-03 batch compared the elpd difference against se_elpd_loo, the SE of ONE",
                "model's total. That quantity is dominated by variation ACROSS observations, which is",
                "common to both models and cancels in the difference; the model-comparison statistic is",
                "sd(d) * sqrt(n) on the paired per-observation differences (Vehtari, Gelman & Gabry 2017,",
                "Stat Comput 27:1413-1432, sec 3.3). READ THE ZERO/POSITIVE SPLIT, not just the total: a",
                "mixture that earns its elpd at the zeros while degrading the positive counts is a",
                "different object from one that fits better everywhere, and the harvest estimate is made",
                "of the positives."))
  }
  # --- D2: the zero bin under the mixture ----------------------------------
  # Recomputed from the committed per-observation NB2 p_zero and the posterior mean of
  # theta_C. This carries a mean-product approximation, E[theta * p0] = E[theta] * E[p0];
  # the exact value needs the joint draws, which is what stage Z2 produces.
  for (nm in c("shore_all_gear", "shore_ring_net_only")) {
    b <- rd(z1, sprintf("ppc_byobs_%s_Dungeness_Kept.csv", nm))
    a <- rd(z0, sprintf("ppc_byobs_%s_Dungeness_Kept.csv", nm))
    s <- rd(z1, sprintf("bss_full_summary_%s_Dungeness_Kept.csv", nm))
    if (is.null(a) || is.null(b) || is.null(s)) next
    bc <- b[b$data_type == "catch", ]; ac <- a[a$data_type == "catch", ]
    th <- suppressWarnings(as.numeric(s[s[[1]] == "theta_C_out", "mean"]))[1]
    if (!is.finite(th)) next
    obs0 <- sum(bc$observed == 0)
    pm   <- th + (1 - th) * bc$p_zero
    zmix <- (obs0 - sum(pm)) / sqrt(sum(pm * (1 - pm)))
    znb  <- (obs0 - sum(bc$p_zero)) / sqrt(sum(bc$p_zero * (1 - bc$p_zero)))
    z0z  <- (sum(ac$observed == 0) - sum(ac$p_zero)) / sqrt(sum(ac$p_zero * (1 - ac$p_zero)))
    V1row("D0", sprintf("ZINB zero bin on %s, scored under the MIXTURE", nm),
          sprintf(paste("theta_C %.4f; %d observed zeros vs %.1f expected, z = %+.1f",
                        "(as reported under NB2: %.1f, z = %+.1f; NB2 baseline: %.1f, z = %+.1f)"),
                  th, obs0, sum(pm), zmix, sum(bc$p_zero), znb, sum(ac$p_zero), z0z),
          "|z| under about 2", if (abs(zmix) < 2.5) "PASS" else "REVIEW",
          paste("pit_block() and calib() computed the catch stream's PIT, coverage and p_zero from a",
                "plain NB2 while the fit was a mixture, so theta_C was excluded from the model's own",
                "zero probability. Under the mixture P(Y=0) = theta + (1-theta)*NB2(0). elpd was never",
                "affected: Stan's log_lik carried the mixture from the start. NOTE the approximation:",
                "this uses E[theta]*E[p0] in place of E[theta*p0] because the committed files hold only",
                "the marginal means. Stage Z2 renders the exact value."))
  }
  # --- D1: what the ladder actually did ------------------------------------
  lg <- rd(.here("05_output", "20260902/pooled-CPUE-SZ-A1-shore-ladder"), "ar_escalation_log.csv")
  if (!is.null(lg)) {
    sa <- lg[grepl("^shore_all_gear", lg$fit), ]
    V1row("D0", "did the 2026-09-03 ladder produce distinct rungs?",
          sprintf("%d attempts, %d DISTINCT resolution(s): %s; catch_median values: %s; all pass_convergence = %s",
                  nrow(sa), length(unique(sa$ar_resolution)), paste(unique(sa$ar_resolution), collapse = ", "),
                  paste(unique(format(sa$catch_median, nsmall = 3)), collapse = ", "),
                  paste(unique(sa$pass_convergence), collapse = "/")),
          "one distinct resolution per attempt", if (length(unique(sa$ar_resolution)) == nrow(sa)) "PASS" else "FAIL",
          paste("The A1 stage's ar_escalate was list(shore = \"all_gear\") and the driver gated",
                "force_resolution on isTRUE(), which is FALSE for a list. The ladder was built correctly",
                "and the loop ran the right number of attempts; only the resolution failed to reach the",
                "data prep. Note pass_convergence = TRUE on every attempt, so the batch's separate claim",
                "that rungs were REJECTED and fell back to PE was also wrong: it inferred rejection from",
                "per-rung adequacy columns that only the REPORTED fit ever writes."))
  }
  # Same canonical path the final block uses, so a crash mid-batch leaves ONE verdicts
  # file rather than a dated duplicate that silently diverges from it.
  merge_csv_by(do.call(rbind, V), .here("05_output", "ladder_zinb_2026-09-04_verdicts.csv"),
               c("stage", "criterion"))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# RUN STAGES
# ---------------------------------------------------------------------------
run_stage <- function(sid) {
  st  <- STAGE_DEFS[[sid]]
  cfg <- resolve_cfg(sid)
  out <- .here("05_output", format(Sys.Date(), "%Y%m%d"), paste0(prefix[[st$model]], st$tag))
  banner(sprintf("%s  %s", sid, st$headline))
  cat("  run_tag:", st$tag, "\n  delta :\n")
  for (k in names(st$delta)) cat(sprintf("    %-28s %s\n", k,
    paste(utils::capture.output(str(st$delta[[k]])), collapse = " ")))
  if (sid == "L1")
    cat(sprintf("  ladder : %s   (daily rung %s)\n", paste(LADDER, collapse = " -> "),
                if (isTRUE(LADDER_INCLUDE_DAILY)) "REFITTED here"
                else sprintf("taken from %s: catch %.3f, PI rel %.4f", REF$DAILY$dir,
                             REF$DAILY$shore_ag, REF$DAILY$pi_rel)))
  existing <- find_outdir(st$model, st$tag)
  # Completion marker: run_parameters.txt is written by the driver's LAST chunk. Testing
  # the first shore summary instead would let a stage that crashed after the shore fit
  # (before the boat fits, the diagnostics and the draw persistence) resume as "complete".
  if (isTRUE(RESUME) && !is.na(existing) &&
      file.exists(file.path(existing, "run_parameters.txt")) &&
      file.exists(file.path(existing, "bss_summary_shore_all_gear_Dungeness_Kept.csv"))) {
    cat("  RESUME: output already present at", basename(existing), "- skipping the fit.\n")
    return(existing)
  }
  if (isTRUE(DRY_RUN)) { cat("  DRY_RUN: not fitting.\n"); return(NA_character_) }
  # 2026-09-05 DEFECT FIX, caught in the pre-run audit. The driver names its OWN output
  # folder from run_config$run_tag (BSS-GH-pooled-CPUE-model.Rmd line ~168) and ignores
  # any `output_dir` placed in the render environment. The first version of this runner
  # set run_tag as a bare variable and passed output_dir= to render(), neither of which
  # the driver reads, so BOTH stages would have written into
  # pooled-CPUE-boat-count-validation-run (the run_config.R default), Z2 overwriting the
  # L1 ladder it had just spent 2-3 h producing, and RESUME would never have found
  # either folder. The tag goes INSIDE the config, as run_shore_ar_zi_2026-09-03.R did.
  cfg$model <- st$model; cfg$run_tag <- st$tag; cfg$run_weather <- FALSE
  run_env <- new.env(parent = globalenv())
  run_env$run_config <- cfg
  t0 <- Sys.time()
  html <- rmarkdown::render(model_rmd[[st$model]], envir = run_env, quiet = FALSE)
  # 2026-09-06 DEFECT FIX. rmarkdown writes the rendered HTML NEXT TO THE .Rmd, not into
  # output_dir, and every driver names its CSV folder itself. The 2026-09-03 runner moved
  # the file afterwards; this one did not, so on the 2026-09-04 run stage L1 rendered its
  # report to 01_BSS_models/BSS-GH-pooled-CPUE-model.html and stage Z2 then OVERWROTE it.
  # The L1 HTML, which is where the AR ladder table is tabulated for a reader, is gone and
  # can only be regenerated by refitting; ar_escalation_log.csv preserved the numbers, so
  # the loss was cosmetic, but it also left a 6,000-line render artefact committed inside
  # 01_BSS_models/. Move it, and say so if the move fails.
  od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e) NA_character_)
  if (!is.na(od) && dir.exists(od) && file.exists(html) &&
      normalizePath(dirname(html)) != normalizePath(od)) {
    if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE)))
      suppressWarnings(file.remove(html))
    else cat(sprintf("  WARNING: could not move %s into %s; the NEXT stage will overwrite it.\n",
                     basename(html), basename(od)))
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
# VERDICTS FOR THE FITTED STAGES
# ---------------------------------------------------------------------------
verdict_L1 <- function(dir, stage_id = "L1") {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  lg <- rd(dir, "ar_escalation_log.csv"); ad <- rd(dir, "model_adequacy.csv")
  if (is.null(lg)) return(invisible(NULL))
  sa <- lg[grepl("^shore_all_gear", lg$fit), ]
  # THE control that the 2026-09-03 batch could not run: distinct rungs.
  V1row(stage_id, "the ladder produced one DISTINCT resolution per rung",
        sprintf("%d attempts: %s; P_n: %s", nrow(sa), paste(sa$ar_resolution, collapse = ", "),
                paste(sa$P_n, collapse = ", ")),
        "no two rungs share a resolution, and each carries its own adequacy",
        if (nrow(sa) != length(unique(sa$ar_resolution))) "FAIL"
        else if (!all(c("p_loo_frac", "n_pareto_bad") %in% names(sa)) ||
                 all(is.na(sa$p_loo_frac))) "PASS (adequacy missing)" else "PASS",
        paste("This is the assertion whose absence let the 2026-09-03 batch burn 14.7 h refitting one",
              "model four times. Identical P_n across rungs is the cheapest tell and is checked here",
              "as well as the resolution string."))
  # Rung comparison. The user's rule, from the FWC meeting: escalate on the BSS convergence
  # gate, report the finest rung that passes; if a margin is needed, prefer the narrowest
  # prediction interval. Read adequacy BESIDE that, never as a gate.
  # 2026-09-06: the adequacy columns are now per rung (bss_rung_adequacy.R). On the
  # 2026-09-04 run they existed for the REPORTED rung only, so the two coarser rungs cost
  # 174 minutes and left nothing to compare. Print them; they, not the gate, decide.
  has_adq <- all(c("p_loo_frac", "n_pareto_bad", "cov50_gear") %in% names(sa))
  rung <- paste(sprintf("%s: catch %s [%s, %s], PI rel %.4f, div %s, gate %s%s",
                        sa$ar_resolution, fmt(sa$catch_median), fmt(sa$catch_lo95), fmt(sa$catch_hi95),
                        sa$pi_width_rel, sa$divergences, sa$pass_convergence,
                        if (has_adq) sprintf(", p_loo %.1f%% of n_obs, bad k %s, cov50 gear %.3f / catch %.3f",
                                             100 * sa$p_loo_frac, sa$n_pareto_bad, sa$cov50_gear, sa$cov50_catch)
                        else ", adequacy NOT RECORDED for this rung"), collapse = " | ")
  if (!isTRUE(LADDER_INCLUDE_DAILY))
    rung <- paste0(sprintf("daily (from %s): catch %s [%s, %s], PI rel %.4f, div %d | ",
                           REF$DAILY$dir, fmt(REF$DAILY$shore_ag), fmt(REF$DAILY$lo95),
                           fmt(REF$DAILY$hi95), REF$DAILY$pi_rel, REF$DAILY$divergences), rung)
  V1row(stage_id, "where does the shore all-gear effort process actually belong?", rung,
        "finest rung that passes the gate; narrowest PI only as a tie-break",
        "READ",
        paste("THE QUESTION THIS BATCH EXISTS FOR, asked properly for the first time. Production fits",
              "this component at daily and it is the worst-behaved fit in the production run: p_loo",
              "35.2% of n_obs, 41 Pareto k above 0.7, coverage_50 0.701 at +7.1 sampling SDs. The gear",
              "track fits it at monthly and gets coverage_50 0.035 at -16.4 SDs. Both ends are bad, so",
              "the answer is probably between them and the ladder is how it is found. DECISION RULE, as",
              "agreed: escalate on the BSS convergence gate and report the finest rung that PASSES. Use",
              "the narrowest relative prediction interval only to break a tie, and note the caution from",
              "the Stage 5 boat 2x2, where the two MISCALIBRATED cells had the narrowest relative",
              "intervals: a narrow interval from a badly calibrated fit is false precision, so check",
              "cov50 and the Pareto k count on any rung a PI tie-break selects. If the answer is not",
              "daily, ar_max_resolution for shore all-gear should change; nothing is published, so this",
              "is a change to a working estimate, not to a released one."))
  if (!is.null(ad)) {
    r <- ad[grepl("shore_all_gear", ad$fit), ]
    if (nrow(r))
      V1row(stage_id, "adequacy of the REPORTED rung, read beside the gate and never as a gate",
            sprintf("p_loo %.1f%% of n_obs (daily was %.1f%%), Pareto k>0.7 %s (was %d), worst cov50 dev %.4f on %s, miscalibrated flag %s",
                    100 * r$p_loo_frac[1], 100 * REF$DAILY$p_loo_frac, r$n_pareto_bad[1], REF$DAILY$bad_k,
                    r$cov50_worst_dev[1], r$cov50_worst_stream[1], r$flag_miscalibrated[1]),
            "informational", "INFO",
            paste("The convergence gate answers 'did this fit sample'. It does not answer 'is this model",
                  "right'. p_loo at 35% of n_obs on the daily rung means the fit is spending about one",
                  "effective parameter per three observations, which is the signature of an AR that is",
                  "too fine for the data to identify, and it is invisible to the gate."))
  }
}

verdict_Z2 <- function(dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  z0 <- .here("05_output", REF$Z0$dir)
  # The 0 AND 1 bins, rendered under the corrected likelihood. Both are needed: the zero bin
  # alone passed on the 2026-09-03 prototype while the misfit migrated one bin over.
  zbin <- function(obs, p) (sum(obs) - sum(p)) / sqrt(sum(p * (1 - p)))
  for (nm in c("shore_all_gear", "shore_ring_net_only")) {
    b <- rd(dir, sprintf("ppc_byobs_%s_Dungeness_Kept.csv", nm))
    a <- rd(z0,  sprintf("ppc_byobs_%s_Dungeness_Kept.csv", nm))
    if (is.null(a) || is.null(b)) next
    bc <- b[b$data_type == "catch", ]; ac <- a[a$data_type == "catch", ]
    z0_ <- zbin(bc$observed == 0, bc$p_zero); zb0 <- zbin(ac$observed == 0, ac$p_zero)
    has1 <- "p_one" %in% names(bc)
    z1_ <- if (has1) zbin(bc$observed == 1, bc$p_one) else NA_real_
    # NB2 baseline for the one bin, recomputed from Z0 by the same rule if that file predates
    # p_one (it does): P(Y=1) needs the draws, so read it from the persisted ppc_draws if
    # present, otherwise report the baseline as unavailable rather than guessing.
    zb1 <- if ("p_one" %in% names(ac)) zbin(ac$observed == 1, ac$p_one) else NA_real_
    approx0 <- if (nm == "shore_all_gear") REF$Z1$zero_z_ag else REF$Z1$zero_z_pc
    V1row("Z2", sprintf("count bins 0 and 1 on %s, RENDERED under the mixture", nm),
          sprintf(paste("ZERO: %d obs vs %.1f exp, z = %+.1f (NB2 baseline z %+.1f; D0 approximation %+.1f).",
                        "ONE: %d obs vs %s exp, z = %s (NB2 baseline z %s)"),
                  sum(bc$observed == 0), sum(bc$p_zero), z0_, zb0, approx0,
                  sum(bc$observed == 1), if (has1) fmt(sum(bc$p_one)) else "NA",
                  if (has1) sprintf("%+.1f", z1_) else "NA",
                  if (is.finite(zb1)) sprintf("%+.1f", zb1) else "NA (Z0 predates p_one)"),
          "BOTH bins |z| under about 2.5; zero bin within 0.2 of the D0 approximation",
          if (is.finite(z1_)) { if (abs(z0_) < 2.5 && abs(z1_) < 2.5) "PASS" else "MISFIT MIGRATED" } else "REVIEW",
          paste("THE ADOPTION CRITERION, SET BEFORE THE RUN. The 2026-09-03 prototype halved the zero-bin",
                "z (+3.8 -> +2.0) and the elpd split by count size showed y=1 losing -42.0 nats at 16.7 SE:",
                "the ZINB moved the misfit rather than removing it. A count-bin table passes as a WHOLE or",
                "not at all. If the one bin fails here, the zero-inflation mechanism is the wrong shape for",
                "these data (excess at 0 AND 1 is a two-regime signature) and the next model to try is a",
                "hurdle or a two-component NB mixture, not a retuned theta prior. On the zero bin, a gap",
                "above 0.2 in z from D0's E[theta]*E[p0] approximation means the theta/p0 covariance is",
                "not negligible and that approximation must not be used again."))
  }
  # elpd by count size, from this run's own pointwise files against Z0.
  for (nm in c("shore_all_gear", "shore_ring_net_only")) {
    f <- sprintf("loo_pointwise_catch_%s_Dungeness_Kept.csv", nm)
    r <- loo_elpd_paired(file.path(z0, f), file.path(dir, f), nm)
    if (is.null(r)) next
    V1row("Z2", sprintf("ZINB elpd on %s by count size, and paired total", nm),
          paste(loo_elpd_by_count_str(r), " || TOTAL: ", loo_elpd_paired_str(r)),
          "3+ bins improve (they carry ~87% of catch); total >= 2 paired SE; y=1 loss is the price",
          if (isTRUE(r$ratio >= 2) && all(r$by_count[c("3-4", "5-8"), "diff"] > 0)) "WORTH IT (read the y=1 row)" else "READ",
          paste("Same fit as Z1 up to RNG, so this should reproduce D0's table; it is here so the",
                "adoption decision and the rendered PPC come from the same folder. WORTH IT here means",
                "the harvest-carrying counts fit better and the total clears 2 SE. It does not mean the",
                "model is right: read it with the count-bin row above."))
  }
  # Calibration, which was unreadable in Z1.
  cal <- rd(dir, "ppc_calibration_shore_all_gear_Dungeness_Kept.csv")
  cz  <- rd(z0,  "ppc_calibration_shore_all_gear_Dungeness_Kept.csv")
  if (!is.null(cal) && !is.null(cz)) {
    g <- function(x) x[x$data_type == "catch", ]
    V1row("Z2", "catch-stream calibration under the correct likelihood",
          sprintf("ZINB: cov50 %.3f, PIT mean %.4f, PIT sd %.4f | NB2: cov50 %.3f, PIT mean %.4f, PIT sd %.4f | Z1 as reported (WRONG): cov50 0.405, PIT mean 0.4325",
                  g(cal)$coverage_50, g(cal)$pit_mean, g(cal)$pit_sd,
                  g(cz)$coverage_50,  g(cz)$pit_mean,  g(cz)$pit_sd),
          "cov50 nearer 0.50; PIT mean not worse by more than 1 sampling SD",
          # 2026-09-06: the original form demanded BOTH statistics improve and fired REVIEW
          # on the 2026-09-04 run over a pit_mean difference of 0.0019, which is noise. The
          # sampling SD of a PIT mean on n observations is about 1/sqrt(12 n) = 0.0071 here,
          # so require only that pit_mean not DEGRADE beyond that.
          if (abs(g(cal)$coverage_50 - 0.5) <= abs(g(cz)$coverage_50 - 0.5) &&
              abs(g(cal)$pit_mean - 0.5) <= abs(g(cz)$pit_mean - 0.5) + 1 / sqrt(12 * g(cal)$n))
            "PASS" else "REVIEW",
          paste("Z1 reported cov50 0.405 and PIT mean 0.4325 for this stream and raised",
                "flag_pit_bias = TRUE on the strength of it. Both were computed under NB2 on a mixture",
                "fit: with theta absorbing the excess zeros lambda_C rises, and scoring the observed",
                "zeros against the risen NB2 mean alone drags the PIT down. This row is the first",
                "readable calibration for a ZINB fit in this project."))
  }
  # The negative control still has to hold.
  fe <- fit_exactness(dir, .here("05_output", REF$Z1$dir), pat = "private_boat",
                      what = "BOAT fits vs the 2026-09-03 Z1",
                      expect_delta = c("save_ppc_draws", "save_ppc_draws_max"))
  V1row("Z2", "the boat fits reproduce the 2026-09-03 Z1 run exactly", fe$observed,
        "bit-identical", fe$verdict,
        paste("zi_catch is set per FIT from catch_zi_populations, so the boat fits carry the feature",
              "off in both runs. The only code that changed between them is the R-side PPC arithmetic,",
              "which runs AFTER sampling and cannot touch a posterior. A FAIL here would mean the",
              "diagnostics patch reached the model, which it must not."))
  # Draw persistence, so this class of defect stops costing a re-fit.
  n_rds <- length(list.files(dir, pattern = "^ppc_draws_.*\\.rds$"))
  V1row("Z2", "PPC draws are persisted for every fit",
        sprintf("%d ppc_draws_*.rds present in %s", n_rds, basename(dir)),
        "one per fitted component, ON THE MACHINE THAT RAN THE FIT",
        if (n_rds >= 4) "PASS" else "NOT CHECKABLE HERE (gitignored)",
        paste("*.rds is gitignored, so these files stay on the machine that ran the fit and do",
              "NOT travel with a pushed result. A zero here on a cloned repo means the files were",
              "not committed, which is intended, not that they were not written; the run that",
              "produced them reported 4. Three diagnostic defects in five weeks have each forced a",
              "multi-hour re-fit to correct a file the fit itself was never wrong about, which is",
              "what these are for, but note that as of 2026-09-06 NO recompute path reads them",
              "back, so they are still write-only. See 03_R_functions/save_bss_ppc_draws.R."))
}

# C1: the candidate. Two isolations and one control, and the adequacy that decides.
verdict_C1 <- function(dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  L1d <- find_outdir("pooled", STAGE_DEFS$L1$tag); Z2d <- find_outdir("pooled", STAGE_DEFS$Z2$tag)
  ad <- rd(dir, "model_adequacy.csv")
  if (!is.null(ad)) {
    r <- ad[grepl("shore_all_gear", ad$fit), ]
    if (nrow(r))
      V1row("C1", "does the candidate keep the ladder's adequacy gain with the ZINB on?",
            sprintf(paste("p_loo %.1f%% of n_obs (daily+NB2 35.2%%, weekly+NB2 9.6%%, daily+ZINB 34.4%%);",
                          "Pareto k>0.7 %s (41 / 1 / 26); worst cov50 dev %.4f on %s; miscalibrated %s"),
                    100 * r$p_loo_frac[1], r$n_pareto_bad[1], r$cov50_worst_dev[1],
                    r$cov50_worst_stream[1], r$flag_miscalibrated[1]),
            "p_loo near the weekly figure and the miscalibration flag clear",
            if (isTRUE(r$p_loo_frac[1] < 0.15) && !isTRUE(as.logical(r$flag_miscalibrated[1]))) "PASS" else "REVIEW",
            paste("The whole point of running them together. Zero-inflation is a CATCH-stream change and",
                  "the AR resolution is an EFFORT-process change, so on the arithmetic they should not",
                  "interact and this should land on the weekly figures. If it does not, they DO interact",
                  "and neither result from the separate runs transfers to the combination."))
  }
  # isolate the likelihood: C1 vs L1 (both weekly, ZINB on/off)
  for (nm in c("shore_all_gear", "shore_ring_net_only")) {
    f <- sprintf("loo_pointwise_catch_%s_Dungeness_Kept.csv", nm)
    r <- loo_elpd_paired(file.path(L1d, f), file.path(dir, f), nm)
    if (is.null(r)) next
    V1row("C1", sprintf("the ZINB at WEEKLY: %s, elpd by count size", nm),
          paste(loo_elpd_by_count_str(r), " || TOTAL: ", loo_elpd_paired_str(r)),
          "at least 2 paired SE, and the 3+ bins improving",
          if (isTRUE(r$ratio >= 2)) "WORTH IT" else if (isTRUE(r$ratio >= 1)) "MARGINAL" else "NOT WORTH IT",
          paste("At daily the ZINB was worth +14.8 nats at 2.69 paired SE. The daily fit was overfitted,",
                "and an overfitted effort process can absorb structure the catch likelihood would",
                "otherwise have to explain, so the gain may be smaller at weekly. A gain that VANISHES",
                "at weekly would mean theta_C was compensating for the AR, not for structural zeros."))
  }
  # the zero and one bins at the new resolution
  zbin <- function(o, p) (sum(o) - sum(p)) / sqrt(sum(p * (1 - p)))
  for (nm in c("shore_all_gear", "shore_ring_net_only")) {
    b <- rd(dir, sprintf("ppc_byobs_%s_Dungeness_Kept.csv", nm))
    a <- rd(L1d, sprintf("ppc_byobs_%s_Dungeness_Kept.csv", nm))
    if (is.null(a) || is.null(b) || !("p_one" %in% names(b))) next
    bc <- b[b$data_type == "catch", ]; ac <- a[a$data_type == "catch", ]
    V1row("C1", sprintf("count bins 0 and 1 at WEEKLY: %s", nm),
          sprintf(paste("ZINB: zero %d vs %.1f (z %+.1f), one %d vs %.1f (z %+.1f) |",
                        "NB2 at the same resolution: zero %.1f (z %+.1f), one %.1f (z %+.1f)"),
                  sum(bc$observed == 0), sum(bc$p_zero), zbin(bc$observed == 0, bc$p_zero),
                  sum(bc$observed == 1), sum(bc$p_one),  zbin(bc$observed == 1, bc$p_one),
                  sum(ac$p_zero), zbin(ac$observed == 0, ac$p_zero),
                  sum(ac$p_one),  zbin(ac$observed == 1, ac$p_one)),
          "both bins improve on the NB2 at the SAME resolution, and both |z| under about 2.5",
          if (abs(zbin(bc$observed == 0, bc$p_zero)) < 2.5 &&
              abs(zbin(bc$observed == 1, bc$p_one)) < 2.5) "PASS" else "REVIEW",
          paste("THE CLEAN COMPARISON the 2026-09-04 batch could not make: same AR, same data, one",
                "likelihood difference. At daily the ZINB halved both bins against an NB2 baseline",
                "measured at a DIFFERENT resolution, because Z0 predates the p_one column."))
  }
  # the boat, untouched in every one of these runs
  fe <- fit_exactness(dir, L1d, pat = "private_boat", what = "BOAT fits vs L1",
                      # the ladder keys differ because L1 ran one and C1 does not, and the
                      # three *_dates entries are derived date lists rather than settings.
                      expect_delta = c("ar_force", "estimate_catch_zi", "catch_zi_populations",
                                       "ar_escalate", "ar_escalate_ladder", "ar_escalate_stop",
                                       "ar_escalate_select", "ar_escalate_max_attempts",
                                       "ar_escalate_respect_cap", "ar_rung_adequacy",
                                       "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1row("C1", "the boat fits are still an untouched negative control", fe$observed,
        "bit-identical to L1", fe$verdict,
        paste("Neither change is scoped to the boat: ar_force names shore/all_gear and",
              "catch_zi_populations names shore. The boat is 43% of the port total, so a FAIL here",
              "would mean the headline moved for a reason nobody chose. This control has held three",
              "times; it costs nothing and it is the reason a shore result stays attributable."))
  # what it does to the reported number
  pt <- rd(dir, "port_total_Dungeness_Kept.csv")
  if (!is.null(pt)) {
    row <- pt[pt$Estimate == "Expected_Catch", ]
    V1row("C1", "the candidate port total", 
          sprintf("%s [%s, %s]  (production 71,513; daily+NB2 71,450; weekly+NB2 72,122; daily+ZINB 71,287)",
                  fmt(row$BSS_median, 0), fmt(row$BSS_lo95, 0), fmt(row$BSS_hi95, 0)),
          "informational", "INFO",
          paste("The two changes push in opposite directions on the total: weekly AR added about +672",
                "and the ZINB removed about -163, so the candidate should land near 71,950. Nothing is",
                "published, so this is a change to a working estimate; it still needs to be a recorded",
                "decision with the adequacy table attached, not a config flip."))
  }
}

# ---------------------------------------------------------------------------
# DRIVE
# ---------------------------------------------------------------------------
banner(sprintf("LADDER + ZINB BATCH  2026-09-04   DRY_RUN=%s  stages: %s",
               DRY_RUN, paste(STAGES, collapse = ", ")))
cat(sprintf("  ladder rungs      : %s\n", paste(LADDER, collapse = " -> ")))
cat(sprintf("  daily rung        : %s\n", if (isTRUE(LADDER_INCLUDE_DAILY)) "refitted in L1" else
                                            paste("reused from", REF$DAILY$dir)))
cat(sprintf("  ZINB re-render    : %s\n", ZINB_RERENDER))
preflight()

dirs <- list()
for (sid in STAGES) {
  if (sid == "D0") { stage_D0(); next }
  dirs[[sid]] <- run_stage(sid)
}
if ("L1" %in% STAGES) verdict_L1(dirs$L1)
if ("Z2" %in% STAGES) verdict_Z2(dirs$Z2)
# 2026-09-07: pass the STAGE ID. Without it C2's ladder rows were written under
# stage "L1", giving the verdicts file two rows per criterion with the same key: the
# first ladder's reading and the second's, indistinguishable to a reader.
if ("C2" %in% STAGES) verdict_L1(dirs$C2, "C2")
if ("C1" %in% STAGES) verdict_C1(dirs$C1)

if (length(V)) {
  dir.create(desk_dir, recursive = TRUE, showWarnings = FALSE)
  vp <- .here("05_output", "ladder_zinb_2026-09-04_verdicts.csv")
  merge_csv_by(do.call(rbind, V), vp, c("stage", "criterion"))
  banner("VERDICTS")
  for (r in seq_len(length(V))) with(V[[r]], cat(sprintf("  [%s] %-14s %s\n     obs: %s\n",
    stage, verdict, criterion, observed)))
  cat("\n  written to", vp, "\n")
}
