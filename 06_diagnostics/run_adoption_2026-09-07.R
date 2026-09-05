###############################################################################
# ADOPTION RENDER: shore all-gear at WEEKLY with the zero-inflated shore catch
# 2026-09-07. Successor to run_ladder_zinb_2026-09-04.R, whose five stages have all run.
# -----------------------------------------------------------------------------
# WHY THIS RUN EXISTS
#   The 2x2 is complete (PIPELINE_STATUS.md Section 1k). Across daily/weekly x NB2/ZINB
#   the shore all-gear component reads 20,898 / 21,547 / 20,745 / 21,489, an interaction
#   of +95 crab, and the weekly+ZINB corner is the best-behaved fit this project has
#   produced: p_loo 9.5% of n_obs against 35.2%, ZERO Pareto k above 0.7 against 41,
#   coverage_50 deviation 0.069 against 0.204, miscalibration flag clear.
#
#   That corner was reached by stage C1 through `ar_force`, which is an EXPERIMENT lever:
#   it bypasses both the data-driven selector and `ar_max_resolution`. Production has to
#   reach the same place through the cap, by letting the selector pick daily and then
#   coarsening. Same resolution, different routing, and a routing change that silently
#   moved a number would be the worst possible way to adopt this.
#
#   So this run is not a new experiment. It is a ROUTING PROOF plus a re-baseline:
#     A1  render production as run_config now ships it, and require the four fits to be
#         BIT-IDENTICAL to C1. If they are, that folder becomes the authoritative run.
#     A2  the gear-resolved cross-check, which has not been re-measured since the shore
#         component moved. Two independently parameterized CPUE structures fitting the
#         same data is the strongest internal check the project has.
#
# WHAT WAS ALREADY CHANGED, AND HOW TO UNDO IT
#   run_config.R now ships ar_max_resolution$pooled$shore$all_gear = "weekly" and
#   estimate_catch_zi = TRUE. That is deliberate: the authoritative run should come out of
#   run_config as shipped, not out of a batch override, or its provenance is a footnote.
#   The cost of that choice is that a FAILED gate leaves run_config carrying a
#   configuration nobody validated, so if A1 fails, REVERT BOTH LINES TOGETHER and say so.
#   The desk pre-flight below refuses to start unless run_config resolves to exactly the
#   configuration C1 fitted, which is the cheap half of the check; the expensive half is
#   the bit-identity gate after the render.
#
# RUNTIME. A1 about 4 h (C1 took 227 min), A2 about 0.5 h (the gear track took 26 min).
# No Stan recompile: the model file is unchanged since 2026-09-02.
###############################################################################

# ============================ CONTROL BLOCK ================================ #
#            ^^^^ the only lines you normally edit ^^^^

DRY_RUN <- TRUE                    # TRUE: the desk pre-flight runs, nothing is fitted. START HERE.
STAGES  <- c("A1", "A2")
RESUME  <- TRUE

# A2 costs half an hour and answers a question a reviewer asks first ("do your two models
# still agree?"). Set FALSE only if the pooled result is needed urgently.
RUN_GEAR_CROSSCHECK <- TRUE

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

if (!isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages({ library(here); library(rmarkdown) })
  load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
  install.lib <- load.lib[!load.lib %in% installed.packages()]
  for (lib in install.lib) install.packages(lib, dependencies = TRUE)
  invisible(sapply(load.lib, require, character.only = TRUE))
  rstan_options(auto_write = TRUE)
} else {
  # The desk routing proof reads the input files and builds Stan data; it needs the
  # tidyverse verbs and readxl, but NOT rstan and not a Stan toolchain. Attach what is
  # there and degrade gracefully, so a dry run on a machine without them still prints
  # every stage's resolved config and runs the pre-flight, which is most of the point.
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

# 2026-09-08 DEFECT FIX. The two tracks label the port-total row DIFFERENTLY: the pooled
# driver writes "Expected_Catch" and the gear-resolved driver writes "Catch". verdict_A2
# looked for "Expected_Catch" in the gear file, got numeric(0), and `is.finite(numeric(0))
# && ...` is an ERROR in R (not FALSE), so the verdict block aborted AFTER 4.1 h of
# fitting and the run's verdicts were never written. Read the row by pattern, and return
# NA rather than a zero-length vector so a downstream `if` behaves.
.port_row <- function(x) {
  if (is.null(x) || !all(c("Estimate", "BSS_median") %in% names(x))) return(NULL)
  i <- which(grepl("^(Expected_)?Catch$", x$Estimate))
  if (!length(i)) return(NULL)
  as.list(x[i[1], , drop = FALSE])
}
.num1 <- function(x) { v <- suppressWarnings(as.numeric(x)); if (length(v)) v[1] else NA_real_ }

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
# BASELINES, named beside each value.
# ---------------------------------------------------------------------------
REF <- list(
  # The candidate, reached through ar_force. A1 must reproduce this exactly.
  C1 = list(dir = "20260904/pooled-CPUE-LZ-C1-weekly-zi",
            shore_ag = 21489, shore_pc = 6320, port = 72032, lo95 = 53044, hi95 = 101212,
            p_loo_frac = 0.0949, bad_k = 0, cov50_dev = 0.0691),
  # The run this supersedes.
  PROD = list(dir = "20260831/pooled-CPUE-VAL-1-adopted",
              port = 71513, lo95 = 52310, hi95 = 100686, shore_ag = 20898),
  # The gear track at its last good state: the boat all-gear fit sampling after the
  # 2026-09-02 driver fix. A2 is measured against this and against A1.
  GEAR = list(dir = "20260901/gear-type-CPUE-model-SZ-G1-gearfix",
              port = 70953, shore_ag = 20754, boat_ag = 30760, tau_bar = 2.5962,
              gap_vs_pooled_pct = -0.78))

STAGE_DEFS <- list(
  # NO DELTA. That is the point: the authoritative run must come out of run_config as
  # shipped. A batch override here would reintroduce exactly the provenance problem this
  # stage exists to remove.
  A1 = list(id = "A1", model = "pooled", tag = "AD-A1-adopted", delta = list(),
            headline = "PRODUCTION as run_config now ships it: weekly shore AR + ZINB shore catch"),
  A2 = list(id = "A2", model = "gear_resolved", tag = "AD-A2-gear-crosscheck", delta = list(),
            headline = "gear-resolved cross-check, first since the shore component moved"))
if (!isTRUE(RUN_GEAR_CROSSCHECK)) STAGES <- setdiff(STAGES, "A2")

resolve_cfg <- function(sid) modifyList(BASE, STAGE_DEFS[[sid]]$delta %||% list(), keep.null = TRUE)
find_outdir <- function(model, run_tag) {
  hits <- list.dirs(.here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}

# ---------------------------------------------------------------------------
# PRE-FLIGHT. The cheap half of the routing proof, plus the standing runner defects.
# ---------------------------------------------------------------------------
preflight <- function() {
  banner("PRE-FLIGHT")
  fails <- character(0)
  say <- function(ok, msg, detail = "") {
    cat(sprintf("  [%s] %s%s\n", if (ok) "ok" else "XX", msg,
                if (nzchar(detail)) paste0("  --  ", detail) else ""))
    if (!ok) fails <<- c(fails, msg)
  }
  # run_config must already BE the candidate configuration.
  say(identical(BASE$ar_max_resolution$pooled$shore$all_gear, "weekly"),
      "run_config caps pooled shore all-gear at weekly",
      BASE$ar_max_resolution$pooled$shore$all_gear %||% "<absent>")
  say(isTRUE(BASE$estimate_catch_zi) && identical(BASE$catch_zi_populations, "shore"),
      "run_config enables the ZINB catch likelihood, scoped to shore")
  say(is.null(BASE$ar_force) && !isTRUE(BASE$ar_escalate),
      "no experiment lever is active: ar_force NULL and ar_escalate off",
      sprintf("ar_force=%s ar_escalate=%s",
              if (is.null(BASE$ar_force)) "NULL" else "SET", BASE$ar_escalate %||% FALSE))
  # nothing else may have moved
  say(identical(BASE$ar_max_resolution$pooled$shore$pot_closure, "biweekly") &&
      identical(BASE$ar_max_resolution$pooled$private_boat, "monthly"),
      "the other pooled caps are untouched")
  say(identical(BASE$ar_max_resolution$gear_resolved$shore$all_gear, "monthly"),
      "the GEAR track's caps are untouched, so A2 measures the pooled change alone")
  # Stated, not asserted: it is a fact about the models, not something this run can fix.
  .gear_stan <- .here("02_stan_models", "crab_bss_gear_resolved.stan")
  .gear_has_zi <- file.exists(.gear_stan) &&
    any(grepl("zi_catch", readLines(.gear_stan, warn = FALSE), fixed = TRUE))
  cat(sprintf("  [--] NOTE: the gear-resolved Stan model %s a zero-inflation block, so A2's shore fits are\n       plain NB2 while A1's are ZINB. The cross-check now compares UNLIKE catch likelihoods;\n       the ZINB is worth about -0.3%% on the pooled shore component, so it explains a small part\n       of whatever gap A2 reports and none of a large one.\n",
              if (.gear_has_zi) "HAS" else "does NOT have"))
  say(isTRUE(BASE$shared_tau) && identical(BASE$pe_gear_ratio_arm, "match_bss"),
      "the 2026-09-01/02 adoptions are still in place (shared_tau, PE arm)")
  # the expensive-half prerequisite: C1 has to be on disk to gate against
  c1 <- .here("05_output", REF$C1$dir)
  say(dir.exists(c1) && file.exists(file.path(c1, "bss_summary_shore_all_gear_Dungeness_Kept.csv")),
      "the C1 reference run is on disk to gate against", REF$C1$dir)
  # standing runner defects, asserted here as well as in the harness
  .self <- .here("06_diagnostics", "run_adoption_2026-09-07.R")
  .src <- if (file.exists(.self)) readLines(.self, warn = FALSE) else character(0)
  .src <- .src[!grepl("^\\s*#", .src)]
  .calls <- .src[grepl("rmarkdown::render(", .src, fixed = TRUE)]
  say(length(.calls) > 0 && !any(grepl("output_dir", .calls, fixed = TRUE)) &&
      any(grepl("cfg$run_tag <- st$tag", .src, fixed = TRUE)),
      "the tag goes INSIDE run_config and output_dir= is not passed to render()")
  say(any(grepl("file.copy(html", .src, fixed = TRUE)),
      "the rendered HTML is moved into the run folder")
  rule()
  if (length(fails)) {
    cat("  PRE-FLIGHT FAILED:\n"); for (f in fails) cat("   -", f, "\n")
    stop("Fix the above before running.")
  }
  cat("  pre-flight clean.\n")
}

# ---------------------------------------------------------------------------
# DESK: the routing proof, run BEFORE spending four hours on it.
#
# ar_force and ar_max_resolution reach "weekly" by different code paths inside
# bss_select_ar_resolution(): the first takes the `fixed` branch and skips the cap, the
# second takes the `adaptive` branch, selects daily from the data, and is then coarsened.
# If those two produce different Stan data, the bit-identity gate cannot pass and there is
# no point starting the render. This costs seconds and needs no Stan toolchain.
# ---------------------------------------------------------------------------
desk_routing <- function() {
  banner("DESK  routing proof: ar_force vs ar_max_resolution must build the SAME Stan data")
  ok <- tryCatch({
    p <- modifyList(BASE, list(bss_model_file = "crab_bss_pooled.stan"))
    p$ar_max_resolution <- BASE$ar_max_resolution$pooled
    p$crabbing_holiday_dates <- read_crabbing_holidays(p)
    q <- function(e) { s <- tempfile(); sink(s); on.exit(sink()); force(e) }
    dwg <- q(fetch_crab_data(p)); ie <- q(fetch_ie_data(p))
    Le  <- q(if (!is.null(ie) && nrow(ie) > 0) estimate_L_effective(ie, p$pot_open_date, p) else NULL)
    sub <- build_subseasons(p)
    ss  <- sub[[which(vapply(sub, function(x) x$gear_regime == "all_gear", logical(1)))]]
    days <- q(prep_days_crab(ss$start, ss$end, p, L_eff_model = Le))
    summ <- q(prep_population_summary(dwg, "shore", ss$start, ss$end, p))
    # A: the C1 route (experiment lever).  B: the production route (the cap).
    pa <- modifyList(p, list(ar_force = list(shore = list(all_gear = "weekly"))))
    pb <- p
    pb$ar_max_resolution$shore$all_gear <- "daily"   # force the ADAPTIVE branch then cap
    pb$ar_max_resolution$shore$all_gear <- "weekly"
    mk <- function(x) q(prep_bss_crab_pooled(days, summ, "Dungeness_Kept", x, "shore",
                                             gear_regime = ss$gear_regime, ie_data = ie))
    a <- mk(pa); b <- mk(pb)
    keys <- union(names(a), names(b))
    diff <- keys[!vapply(keys, function(k)
      isTRUE(all.equal(a[[k]], b[[k]], tolerance = 0)), logical(1))]
    V1row("DESK", "the two routes to weekly build identical Stan data",
          sprintf("ar_force: %s, P_n %d | ar_max_resolution: %s, P_n %d | entries differing: %s",
                  attr(a, "ar_resolution"), a$P_n, attr(b, "ar_resolution"), b$P_n,
                  if (length(diff)) paste(diff, collapse = ", ") else "NONE"),
          "no stan_data entry differs, and both resolve to weekly",
          if (!length(diff) && identical(attr(a, "ar_resolution"), "weekly") &&
              identical(attr(b, "ar_resolution"), "weekly")) "PASS" else "FAIL",
          paste("The whole adoption rests on this. C1 reached weekly through ar_force, which takes the",
                "`fixed` branch of bss_select_ar_resolution() and is NOT subject to the cap; production",
                "reaches it through the adaptive branch, which selects daily from the effort density and",
                "is then coarsened by ar_max_resolution. Both end at weekly with P_n 44, but that has to",
                "be demonstrated rather than assumed, because a FAIL here means the four-hour render",
                "cannot possibly pass its bit-identity gate and should not be started."))
    cat(sprintf("  routes: ar_force -> %s (P_n %d) | cap -> %s (P_n %d) | differing entries: %s\n",
                attr(a, "ar_resolution"), a$P_n, attr(b, "ar_resolution"), b$P_n,
                if (length(diff)) paste(diff, collapse = ", ") else "NONE"))
    !length(diff)
  }, error = function(e) {
    cat("  desk routing proof unavailable:", conditionMessage(e), "\n")
    V1row("DESK", "the two routes to weekly build identical Stan data",
          paste("not evaluated:", conditionMessage(e)), "no stan_data entry differs", "REVIEW",
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
  banner(sprintf("%s  %s", sid, st$headline))
  cat("  run_tag:", st$tag, "\n")
  cat(sprintf("  config : shore all-gear cap %s | estimate_catch_zi %s | ar_force %s\n",
              cfg$ar_max_resolution[[st$model]]$shore$all_gear %||% "?",
              cfg$estimate_catch_zi %||% FALSE,
              if (is.null(cfg$ar_force)) "NULL" else "SET"))
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
# VERDICTS
# ---------------------------------------------------------------------------
verdict_A1 <- function(dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  c1 <- .here("05_output", REF$C1$dir)
  # THE GATE.
  fe <- fit_exactness(dir, c1, what = "all four fits vs C1",
                      expect_delta = c("ar_force", "ar_max_resolution", "run_tag", "model",
                                       "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1row("A1", "production routing reproduces the candidate EXACTLY", fe$observed,
        "every shared parameter row identical at full precision", fe$verdict,
        paste("C1 reached weekly through ar_force, which bypasses the selector and the cap; this run",
              "reaches it by letting the selector pick daily and then coarsening. Same resolution, same",
              "P_n, different code path. A PASS means the adoption is a routing change and nothing else,",
              "which is what makes it safe to call this the authoritative run. A FAIL means the cap and",
              "the override do not build the same model, and the right response is to REVERT both",
              "run_config lines rather than to prefer whichever number looks better."))
  # adequacy, stated once more so the authoritative run carries it
  ad <- rd(dir, "model_adequacy.csv")
  if (!is.null(ad)) {
    r <- ad[grepl("shore_all_gear", ad$fit), ]
    if (nrow(r))
      V1row("A1", "the adopted shore all-gear fit's adequacy",
            sprintf(paste("p_loo %.1f%% of n_obs (was 35.2%% at daily); Pareto k>0.7 %s (was 41);",
                          "worst cov50 dev %.4f on %s (was 0.2042); miscalibrated %s (was TRUE)"),
                    100 * r$p_loo_frac[1], r$n_pareto_bad[1], r$cov50_worst_dev[1],
                    r$cov50_worst_stream[1], r$flag_miscalibrated[1]),
            "matches C1: p_loo about 9.5%, zero bad Pareto k, flag clear",
            if (isTRUE(r$p_loo_frac[1] < 0.15) && !isTRUE(as.logical(r$flag_miscalibrated[1])))
              "PASS" else "REVIEW",
            paste("This is the number a reviewer should be shown next to the estimate. The convergence",
                  "gate says the fit sampled; these say the model is not spending an effective parameter",
                  "per three observations to do it. The catch stream remains under-covered at about -4",
                  "sampling SD at every AR resolution, which the ZINB improves and does not fix."))
  }
  pt <- rd(dir, "port_total_Dungeness_Kept.csv")
  if (!is.null(pt)) {
    row <- .port_row(pt)
    V1row("A1", "the new authoritative port total",
          sprintf("%s [%s, %s]  (superseded run %s [%s, %s]; C1 %s)",
                  fmt(row$BSS_median, 0), fmt(row$BSS_lo95, 0), fmt(row$BSS_hi95, 0),
                  fmt(REF$PROD$port, 0), fmt(REF$PROD$lo95, 0), fmt(REF$PROD$hi95, 0),
                  fmt(REF$C1$port, 0)),
          "within RNG noise of C1 (port totals move about 0.2% between bit-identical fits)",
          "INFO",
          paste("+0.73% on the superseded run, almost all of it the shore effort process no longer being",
                "over-imputed at a daily AR. NOTHING HAS BEEN PUBLISHED, so this supersedes an internal",
                "working number rather than a released one. Note that the port total is assembled by",
                "resampling component draws and rstan::extract(permuted = TRUE) permutes them, so it",
                "moves about 0.2% between bit-identical fits: judge the gate on the per-fit summaries",
                "above, never on this line."))
  }
}

verdict_A2 <- function(dir, pooled_dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  gp <- rd(dir, "port_total_Dungeness_Kept.csv"); pp <- rd(pooled_dir %||% "", "port_total_Dungeness_Kept.csv")
  gc <- rd(dir, "pe_vs_bss_comparison.csv");      pc <- rd(pooled_dir %||% "", "pe_vs_bss_comparison.csv")
  if (is.null(gp)) return(invisible(NULL))
  g_port <- .num1(.port_row(gp)$BSS_median)
  p_port <- .num1(.port_row(pp)$BSS_median)
  gap <- if (isTRUE(is.finite(p_port)) && isTRUE(is.finite(g_port))) 100 * (g_port - p_port) / p_port else NA_real_
  gs <- function(x, k) if (is.null(x)) NA_real_ else x$BSS_catch[x$component == k][1]
  V1row("A2", "the two tracks still agree after the pooled shore component moved",
        sprintf(paste("gear port %s vs pooled %s (%+.2f%%); gear shore all-gear %s vs pooled %s;",
                      "gear boat all-gear %s vs pooled %s. Last measured pair: %s vs 71,513 (-0.78%%)"),
                fmt(g_port, 0), fmt(p_port, 0), gap,
                fmt(gs(gc, "shore (All gear)"), 0), fmt(gs(pc, "shore (All gear)"), 0),
                fmt(gs(gc, "private_boat (All gear)"), 0), fmt(gs(pc, "private_boat (All gear)"), 0),
                fmt(REF$GEAR$port, 0)),
        "within about 2% at the port",
        if (isTRUE(is.finite(gap)) && abs(gap) <= 2) "PASS" else "REVIEW",
        paste("Two independently parameterized CPUE structures fitting the same data, and the strongest",
              "internal check the project has. EXPECT THE GAP TO CHANGE and read it carefully: the",
              "pooled shore all-gear component just moved from daily to weekly while the gear track",
              "still fits that component at MONTHLY (ar_max_resolution$gear_resolved is deliberately",
              "untouched), so the two tracks are now closer in resolution than they were and the gap",
              "should NARROW. If it widens instead, the pooled change did something the gear track",
              "disagrees with and that is worth more attention than the port total itself."))
  V1row("A2", "the two tracks now use DIFFERENT catch likelihoods on the shore",
        sprintf("gear-resolved Stan carries a zero-inflation block: %s. estimate_catch_zi is a POOLED-only flag.",
                if (any(grepl("zi_catch", readLines(.here("02_stan_models", "crab_bss_gear_resolved.stan"), warn = FALSE), fixed = TRUE))) "yes" else "NO"),
        "informational, but it changes how the gap is read", "INFO",
        paste("crab_bss_gear_resolved.stan has no theta_C and prep_bss_crab_gear.R never emits zi_catch, so",
              "the gear track ignores estimate_catch_zi entirely and fits plain NB2. As of this run the two",
              "tracks therefore differ in the shore CATCH LIKELIHOOD as well as in the AR resolution, which",
              "is the first time they have differed structurally on the shore. The ZINB is worth about",
              "-0.3% on the pooled shore component (21,547 to 21,489 at weekly), so it accounts for a small",
              "part of any gap and none of a large one. Porting the block to the gear model would restore",
              "symmetry and is a Stan edit plus a recompile; it is NOT needed to read this cross-check, but",
              "it should be recorded as the reason the two tracks are no longer like for like."))
  # THE comparison that says whether the tracks actually disagree. They fit the shore
  # all-gear component at DIFFERENT resolutions now (pooled weekly, gear monthly), so the
  # headline gap conflates a resolution difference with a structural one. The pooled
  # ladder fitted that component at monthly too, so the like-for-like number exists.
  c2 <- rd(.here("05_output", "20260904/pooled-CPUE-LZ-C2-ladder-adequacy"), "ar_escalation_log.csv")
  pm <- if (!is.null(c2)) .num1(c2$catch_median[c2$ar_resolution == "monthly" &
                                                grepl("^shore_all_gear", c2$fit)]) else NA_real_
  gsa <- .num1(gs(gc, "shore (All gear)"))
  V1row("A2", "do the tracks disagree, or are they just at different resolutions?",
        sprintf(paste("shore all-gear: pooled at the ADOPTED weekly %s vs gear at monthly %s (%+.2f%%);",
                      "pooled at MONTHLY (from the C2 ladder) %s vs gear at monthly %s (%+.2f%%)"),
                fmt(.num1(gs(pc, "shore (All gear)")), 0), fmt(gsa, 0),
                100 * (gsa - .num1(gs(pc, "shore (All gear)"))) / .num1(gs(pc, "shore (All gear)")),
                fmt(pm, 0), fmt(gsa, 0),
                if (isTRUE(is.finite(pm))) 100 * (gsa - pm) / pm else NA_real_),
        "the like-for-like gap is much smaller than the headline gap", "READ",
        paste("If the two tracks agree closely AT THE SAME RESOLUTION, the headline gap is a",
              "resolution difference and not a disagreement about the fishery, and the right response",
              "is to decide whether the GEAR track's cap should move too rather than to worry about",
              "the pooled estimate. Note the gear track's shore fits carry per-gear CPUE (G = 5), a",
              "thinner likelihood per gear type, so it may genuinely need a coarser AR than the pooled",
              "model does: that needs its own ladder, not an assumption copied across."))
  # The control that proves estimate_catch_zi really is inert on this track.
  fe <- fit_exactness(dir, .here("05_output", REF$GEAR$dir), what = "gear fits vs the G1 baseline",
                      expect_delta = c("ar_max_resolution", "estimate_catch_zi", "run_tag", "model",
                                       "save_ppc_draws", "save_ppc_draws_max", "ar_rung_adequacy",
                                       "crabbing_holiday_dates", "opener_f_dates", "razor_dig_dates"))
  V1row("A2", "the gear track did not move at all", fe$observed,
        "bit-identical to the 2026-09-01 G1 run", fe$verdict,
        paste("Nothing in this run's config reaches the gear track: its ar_max_resolution map is",
              "untouched and crab_bss_gear_resolved.stan has no zero-inflation block, so",
              "estimate_catch_zi is inert here. A FAIL would mean the pooled adoption leaked into the",
              "cross-check, which would make the cross-check worthless as an independent reading."))
  gt <- rd(dir, "bss_full_summary_private_boat_all_gear_Dungeness_Kept.csv")
  pt2 <- rd(pooled_dir %||% "", "bss_full_summary_private_boat_all_gear_Dungeness_Kept.csv")
  .tau <- function(x) if (is.null(x)) NA_real_ else .num1(x$mean[grepl("^tau_bar", x[[1]])])
  V1row("A2", "the shared turnover is still a property of the data, not of one parameterization",
        sprintf("gear tau_bar %s vs pooled %s (%+.2f%%); the 2026-09-01 pair agreed to 0.03%%",
                fmt(.tau(gt), 4), fmt(.tau(pt2), 4),
                100 * (.tau(gt) - .tau(pt2)) / .tau(pt2)),
        "close to the pooled value", "INFO",
        paste("tau_bar is estimated independently by the two tracks from different CPUE structures. Its",
              "agreement is what makes the shared-turnover adoption a finding about the fishery rather",
              "than an artefact of the pooled parameterization, and it is worth re-reading whenever",
              "either track's configuration changes."))
}

# ---------------------------------------------------------------------------
# DRIVE
# ---------------------------------------------------------------------------
banner(sprintf("ADOPTION RENDER  2026-09-07   DRY_RUN=%s  stages: %s", DRY_RUN,
               paste(STAGES, collapse = ", ")))
cat(sprintf("  run_config ships: shore all-gear cap %s | estimate_catch_zi %s | ar_force %s\n",
            BASE$ar_max_resolution$pooled$shore$all_gear %||% "?", BASE$estimate_catch_zi %||% FALSE,
            if (is.null(BASE$ar_force)) "NULL" else "SET"))
cat(sprintf("  gate: A1's four fits must be BIT-IDENTICAL to %s\n", REF$C1$dir))
preflight()
desk_routing()

dirs <- list()
for (sid in STAGES) dirs[[sid]] <- run_stage(sid)
# 2026-09-08: every verdict block is wrapped, and the verdicts are written whatever
# happens. This is the THIRD time a defect in a verdict block has aborted a batch AFTER
# all the fitting was done (fmt() on a vector, 2026-09-06; the port-row label, this run),
# and on both occasions the run's own verdicts were lost while the expensive part had
# succeeded. A verdict is a reading of a result; it must not be able to destroy one.
.safe <- function(sid, expr) tryCatch(force(expr), error = function(e) {
  cat(sprintf("\n  *** VERDICT BLOCK %s FAILED: %s\n", sid, conditionMessage(e)))
  cat("      The run itself is intact; re-run with RESUME = TRUE after fixing the block.\n")
  V1row(sid, "verdict block did not complete", conditionMessage(e),
        "the block runs to completion", "ERROR",
        paste("The fits are on disk and unaffected: this is a defect in the code that READS them.",
              "Fix the block and re-run with RESUME = TRUE, which re-scores from disk in seconds."))
  invisible(NULL)
})
if ("A1" %in% STAGES) .safe("A1", verdict_A1(dirs$A1))
if ("A2" %in% STAGES) .safe("A2", verdict_A2(dirs$A2, dirs$A1 %||% .here("05_output", REF$C1$dir)))

if (length(V)) {
  vp <- .here("05_output", "adoption_2026-09-07_verdicts.csv")
  merge_csv_by(do.call(rbind, V), vp, c("stage", "criterion"))
  banner("VERDICTS")
  for (r in seq_len(length(V))) with(V[[r]], cat(sprintf("  [%s] %-14s %s\n     obs: %s\n",
    stage, verdict, criterion, observed)))
  cat("\n  written to", vp, "\n")
}
