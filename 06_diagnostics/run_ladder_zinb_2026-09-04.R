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
# RUNTIME. Roughly 3-6 h depending on the two switches below, against 23.6 h for the run
# it corrects. The saving comes from not refitting the daily rung, which stage Z0 already
# holds bit-identically.
#
# NO STAN RECOMPILE unless ZINB_RERENDER is TRUE and the cached ZINB binary is gone; the
# Stan file is unchanged from 2026-09-03.
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                   # TRUE: desk stages run, nothing is fitted. START HERE.
STAGES  <- c("D0", "L1", "Z2")
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
fmt    <- function(x, d = 1) if (length(x) == 0 || is.na(x)) "NA" else
  formatC(as.numeric(x), format = "f", digits = d, big.mark = ",")

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
            headline = "ZINB shore catch, re-rendered under the CORRECTED PPC"))
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
  if (isTRUE(RESUME) && !is.na(existing) &&
      file.exists(file.path(existing, "bss_summary_shore_all_gear_Dungeness_Kept.csv"))) {
    cat("  RESUME: output already present at", basename(existing), "- skipping the fit.\n")
    return(existing)
  }
  if (isTRUE(DRY_RUN)) { cat("  DRY_RUN: not fitting.\n"); return(NA_character_) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  run_env <- new.env(parent = globalenv())
  assign("run_config", cfg, envir = run_env)
  assign("output_dir", out, envir = run_env)
  assign("run_tag", st$tag, envir = run_env)
  t0 <- Sys.time()
  rmarkdown::render(model_rmd[[st$model]], output_dir = out, envir = run_env, quiet = FALSE)
  cat(sprintf("  %s finished in %.1f min\n", sid,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  out
}

# ---------------------------------------------------------------------------
# VERDICTS FOR THE FITTED STAGES
# ---------------------------------------------------------------------------
verdict_L1 <- function(dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  lg <- rd(dir, "ar_escalation_log.csv"); ad <- rd(dir, "model_adequacy.csv")
  if (is.null(lg)) return(invisible(NULL))
  sa <- lg[grepl("^shore_all_gear", lg$fit), ]
  # THE control that the 2026-09-03 batch could not run: distinct rungs.
  V1row("L1", "the ladder produced one DISTINCT resolution per rung",
        sprintf("%d attempts: %s; P_n: %s", nrow(sa), paste(sa$ar_resolution, collapse = ", "),
                paste(sa$P_n, collapse = ", ")),
        "no two rungs share a resolution",
        if (nrow(sa) == length(unique(sa$ar_resolution))) "PASS" else "FAIL",
        paste("This is the assertion whose absence let the 2026-09-03 batch burn 14.7 h refitting one",
              "model four times. Identical P_n across rungs is the cheapest tell and is checked here",
              "as well as the resolution string."))
  # Rung comparison. The user's rule, from the FWC meeting: escalate on the BSS convergence
  # gate, report the finest rung that passes; if a margin is needed, prefer the narrowest
  # prediction interval. Read adequacy BESIDE that, never as a gate.
  rung <- paste(sprintf("%s: catch %s [%s, %s], PI rel %.4f, div %s, gate %s",
                        sa$ar_resolution, fmt(sa$catch_median), fmt(sa$catch_lo95), fmt(sa$catch_hi95),
                        sa$pi_width_rel, sa$n_divergent %||% NA, sa$pass_convergence), collapse = " | ")
  if (!isTRUE(LADDER_INCLUDE_DAILY))
    rung <- paste0(sprintf("daily (from %s): catch %s [%s, %s], PI rel %.4f, div %d | ",
                           REF$DAILY$dir, fmt(REF$DAILY$shore_ag), fmt(REF$DAILY$lo95),
                           fmt(REF$DAILY$hi95), REF$DAILY$pi_rel, REF$DAILY$divergences), rung)
  V1row("L1", "where does the shore all-gear effort process actually belong?", rung,
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
      V1row("L1", "adequacy of the REPORTED rung, read beside the gate and never as a gate",
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
          "cov50 nearer 0.50 and PIT mean nearer 0.50 than the NB2 baseline",
          if (abs(g(cal)$coverage_50 - 0.5) <= abs(g(cz)$coverage_50 - 0.5) &&
              abs(g(cal)$pit_mean - 0.5) <= abs(g(cz)$pit_mean - 0.5)) "PASS" else "REVIEW",
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
        sprintf("%d ppc_draws_*.rds written", n_rds), "one per fitted component",
        if (n_rds >= 4) "PASS" else "REVIEW",
        paste("Three diagnostic defects in five weeks have each forced a multi-hour re-fit to correct a",
              "file the fit itself was never wrong about. With the draws on disk the next one is a",
              "recomputation. See 03_R_functions/save_bss_ppc_draws.R for what this does and does not cover."))
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

if (length(V)) {
  dir.create(desk_dir, recursive = TRUE, showWarnings = FALSE)
  vp <- .here("05_output", "ladder_zinb_2026-09-04_verdicts.csv")
  merge_csv_by(do.call(rbind, V), vp, c("stage", "criterion"))
  banner("VERDICTS")
  for (r in seq_len(length(V))) with(V[[r]], cat(sprintf("  [%s] %-14s %s\n     obs: %s\n",
    stage, verdict, criterion, observed)))
  cat("\n  written to", vp, "\n")
}
