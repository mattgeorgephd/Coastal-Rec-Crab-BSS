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
# THE GEAR TRACK'S AR PERIOD (D3) AND ITS ZERO-INFLATED SHORE CATCH (D6)
# -----------------------------------------------------------------------------
# WHY THIS RUN EXISTS
#   Two register items are open, both about the same fit (gear-resolved, shore, all
#   gear), and neither can be settled by argument:
#
#   D3  The pooled track fits that component at WEEKLY; this track fits it at MONTHLY.
#       That difference is most of the -1.39% cross-track port gap: at a common monthly
#       resolution the two tracks agree to 0.08% (17 crab), the strongest agreement the
#       cross-check has produced. The open question is NOT "why do the tracks disagree",
#       which is answered; it is "what period should the GEAR track use on its own
#       merits". Do not assume the pooled answer transfers.
#
#       AND A CORRECTION TO THE REASON D3 HAS STAYED OPEN, found while building this file
#       and measured rather than argued. The register and adoption-review-2026-09-08.md
#       both justify not copying the pooled period across by saying the gear shore fits
#       carry per-gear CPUE at "G = 5", a thinner likelihood per gear type. THEY DO NOT.
#       gear_resolved_G ships FALSE and has since at least 2026-09-03; prep_bss_crab_gear
#       then collapses the Stan-facing gear dimension and prints, in the run log of the
#       committed R5 render itself, "NOTE: 4 gear types qualify; G = 1 this run (set
#       gear_resolved_G = TRUE for per-gear CPUE)". Built on the real 2024-25 shore
#       all-gear data on 2026-09-13: G = 1, D = 289, IntC = 1651. So the two tracks fit
#       the SAME latent dimension, P_n x S; the count in the claim is wrong as well as its
#       consequence (four gear types qualify, not five); and the caution becomes true only
#       if gear_resolved_G is switched on, which is GR-7 Phase 2 and a separate decision.
#       What is left between the two tracks is this period and the ZI block below.
#
#   D6  The pooled model has carried a zero-inflated shore catch likelihood since
#       2026-09-02; until 2026-09-13 this track's Stan file had no theta_C at all, so
#       estimate_catch_zi = TRUE was read and silently ignored here and the two tracks
#       compared UNLIKE shore catch likelihoods. The Stan port landed 2026-09-13 and
#       SHIPS OFF for this track (catch_zi_tracks = "pooled"), so nothing changed yet.
#       The open question is whether it earns its parameter on a PER-GEAR likelihood,
#       which is thinner than the pooled one it was adopted on.
#
# WHAT CHANGED IN THE CODE TO MAKE THIS RUNNABLE (2026-09-13, same patch as this file)
#   1. crab_bss_gear_resolved.stan gained the ZI block, line for line from the pooled
#      model: the zi_catch flag, theta_C declared vector<lower=0,upper=1>[zi_catch] so it
#      is ZERO-SIZE when off, the log_mix branch, theta_C_out / zi_scale, and zi_scale on
#      lambda_Ctot_S. prep_bss_crab_gear.R builds the three data variables. ALREADY MEASURED
#      INERT: the real 2024-25 shore all-gear stan_data (D = 289, IntC = 1651, G = 1,
#      zi_catch = 0) sampled at seed 20260619, 2 chains x 300 iterations, under the pre-edit
#      and post-edit models gives BIT-IDENTICAL draws on all 4,950 shared columns, the only
#      new columns being theta_C_out and zi_scale. G1 below re-proves it at production
#      iteration counts against the committed render, which is the version that matters.
#   2. run_config gained catch_zi_tracks (ships "pooled") and gear_period_bss (ships
#      list(all_gear = "month", pot_closure = "biweekly"), which are the exact literals
#      build_subseasons.R used before, so every committed run still reproduces).
#   3. D3 IS NOT THE KEY THE REGISTER NAMES. ar_max_resolution$gear_resolved is DORMANT:
#      production ships ar_adaptive = FALSE, so the gear driver passes
#      fixed_resolution = ss$period_bss and the cap is never consulted. Both run_config
#      and the gear .Rmd say so. Until this patch the real lever was a LITERAL inside
#      build_subseasons.R with no configuration surface at all, which is why "the gear
#      track's shore cap is still monthly" has sat open: changing the cap would have done
#      nothing. gear_period_bss is the lever, and this ladder moves it.
#
# THE RUNGS. One lever moves per rung; everything else is pinned by WINDOW below.
#   G0   desk, seconds. No fit. Proves the prerequisites before spending 3 h: the Stan
#        file parses and declares the ZI variables, the prep builds every declared
#        variable, each rung's resolved configuration differs from G1 ONLY in declared
#        keys, the reference folder for the bit-identity control exists, and the ZI
#        rung's configuration differs from its own control in exactly one key.
#   G1   gear_period_bss$all_gear = "month"     -- THE SHIPPED VALUE, and the control
#   G2   gear_period_bss$all_gear = "biweekly"
#   G3   gear_period_bss$all_gear = "weekly"    -- the pooled track's resolution
#   G4   gear_period_bss$all_gear = "daily"     -- completes the bracket
#   G5   ZI ON at ZI_AT (default "weekly"), catch_zi_tracks = c("pooled","gear_resolved")
#
#   G1 DOES DOUBLE DUTY. Its configuration is the shipped one, so it must come back
#   BIT-IDENTICAL to the committed 20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf
#   render. That is the proof that the D6 Stan edit is inert when off, and it costs
#   nothing because G1 is a rung the ladder needs anyway. If it is not bit-identical, the
#   Stan edit is NOT inert and every other rung in this ladder is uninterpretable, so the
#   verdict for G1 is reported first and the script says so.
#
# THE DECISION RULE FOR D3, STATED HERE BEFORE THE RUN so it cannot be fitted to the
# answer afterwards. This matters more than any single statistic below.
#   1. A rung is ELIGIBLE only if every fit in it passes the convergence gate.
#   2. Reject an eligible rung whose shore all-gear p_loo_frac exceeds 0.15, or whose
#      bad-Pareto-k count exceeds 5% of its catch-stream n_obs. Calibration: the POOLED
#      daily fit, which was rejected as overfitted, sat at p_loo 0.352 with 41 bad k; the
#      pooled weekly fit that replaced it sat at 0.095 with 0. The gear R5 monthly fit
#      sits at 0.0368 with 0.
#   3. DO NOT SELECT ON elpd. It is reported, and so is the bad-k share that makes it
#      unreliable, because on the pooled ladder elpd favoured the daily fit that every
#      other diagnostic said was overfitted, and 13% of that fit's gear-stream Pareto k
#      were above 0.7. An elpd computed from a fit whose importance weights are that
#      unstable is not evidence.
#   4. Among rungs surviving (1) and (2), prefer the FINEST resolution. The cap exists to
#      stop overfitting, not to prefer coarseness.
#   5. THE CROSS-TRACK GAP IS NOT A CRITERION. It is reported at every rung, and it is a
#      corroboration when the two tracks are at a common resolution. Choosing the gear
#      resolution because it makes the gap small would be circular: it would tune one
#      estimate to another estimate rather than to the data. The script refuses to let
#      the gap enter the recommendation, and says so in the output.
#
# THE DECISION RULE FOR D6, likewise stated first.
#   1. The ZI rung must pass the gate.
#   2. The PAIRED elpd gain on the catch stream must exceed 2 paired SE. Paired, not
#      naive: on the pooled Z1/Z0 pair the paired SE was 5.5 nats against a naive 46.7,
#      which is the difference between 2.69 SE and 0.32 SE (loo_elpd_paired.R; Vehtari
#      et al. 2017 s3.3).
#   3. Both count bins must improve. The pooled adoption halved both (zero z +3.7 -> +2.0,
#      one z -6.1 -> -3.3). A feature that fixes the zero bin by breaking the one bin is
#      not doing what it claims.
#   4. The elpd DECOMPOSITION must not show the whole gain bought at the zeros while the
#      positive counts degrade by more than the gain. loo_elpd_paired() splits it.
#   5. THE PORT TOTAL IS NOT A CRITERION. The season total is scaled by (1 - theta_C) in
#      generated quantities precisely so that turning the feature on does not inflate it,
#      so a total that barely moves is the design working, not weak evidence.
#
# RUNTIME, from measured evidence rather than a guess. run_timings.csv in the committed
# gear renders puts the shore all-gear fit at 10.0 min (20260911 R5) and 10.7 min
# (20260905 A2), both at monthly (P_n = 10), inside a whole gear render of about 32 min
# for four fits. On the POOLED track the same fit takes 96-165 min at weekly (P_n = 44)
# and 162 min at daily (P_n = 289), i.e. 6.6x the periods for about 1.5x the time, so
# time is strongly sub-linear in P_n. Applying that shape to a 10-minute monthly fit:
#   G1 month    ~32 min      G2 biweekly ~35 min      G3 weekly   ~40 min
#   G4 daily    ~55 min      G5 ZI       ~40 min                   TOTAL about 3.4 h
# The daily rung is the one most likely to overrun, and it is also the one the decision
# rule is least likely to select. STAGES below can drop it.
#
# WHAT IT WRITES, all merged by key so a partial re-run updates rather than truncates:
#   05_output/gear_ar_zi_2026-09-13_ladder.csv          per-rung totals and adequacy
#   05_output/gear_ar_zi_2026-09-13_verdicts.csv        every criterion, PASS/FAIL/REVIEW
#   05_output/gear_ar_zi_2026-09-13_recommendation.csv  D3 and D6, with the rule applied
# and a printed RECOMMENDATION block that names the rule it applied and the evidence.
#
# SHIPS DRY_RUN <- TRUE. Set it FALSE and source again to fit.
###############################################################################

DRY_RUN <- TRUE                    # TRUE prints the plan and the desk stage; fits nothing
STAGES  <- c("G0", "G1", "G3", "G2", "G4", "G5")   # G1 first after G0: it is the control
RESUME  <- TRUE                    # reuse a rung ONLY when its GEAR_STAGE.txt digest matches
ZI_AT   <- "weekly"                # the resolution G5 turns ZI on at; "d3" uses G3's value

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

GEAR_RMD <- .here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd")
GEAR_STAN <- .here("02_stan_models", "crab_bss_gear_resolved.stan")
GEAR_PREP <- .here("03_R_functions", "prep_bss_crab_gear.R")
PREFIX <- "gear-type-CPUE-model-"
stopifnot(file.exists(GEAR_RMD), file.exists(GEAR_STAN), file.exists(GEAR_PREP))

# THE REFERENCE the bit-identity control compares against: the committed R5 gear
# cross-check, rendered 2026-09-11 on the shipped configuration with the PRE-D6 Stan file.
REF_R5 <- "20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf"
# and the pooled rung whose shore all-gear component the cross-track gap is measured
# against. Reported only; see decision rule 5.
REF_POOLED <- "20260910/pooled-CPUE-IMP-R4-shore-tau-newf"

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
# The adequacy row for one fit. model_adequacy.csv is written per run by
# bss_model_adequacy(); every column the D3 rule reads lives there and nowhere else.
.adq <- function(dir, fit = "shore_all_gear_Dungeness_Kept") {
  x <- rd(dir, "model_adequacy.csv"); if (is.null(x) || !"fit" %in% names(x)) return(NULL)
  i <- which(x$fit == fit); if (!length(i)) return(NULL)
  as.list(x[i[1], , drop = FALSE])
}
# The AR resolution the fit ACTUALLY used, and its period count. A rung that silently
# fell back to a different resolution would corrupt the whole ladder without failing, so
# this is checked per rung rather than assumed from the configuration.
.arlog <- function(dir, fit = "shore_all_gear") {
  x <- rd(dir, "ar_escalation_log.csv"); if (is.null(x) || !"fit" %in% names(x)) return(NULL)
  i <- which(grepl(fit, x$fit)); if (!length(i)) return(NULL)
  as.list(x[i[1], , drop = FALSE])
}
.gate <- function(dir) {
  x <- rd(dir, "convergence_report.csv"); if (is.null(x)) return(NULL)
  x
}
# The count bins, from ppc_byobs. Returns observed zeros/ones, expected, and the z.
.bin <- function(dir, k = 0L, fitlab = "shore_all_gear_Dungeness_Kept", stream = "catch") {
  d <- rd(dir, sprintf("ppc_byobs_%s.csv", fitlab))
  col <- if (k == 0L) "p_zero" else "p_one"
  if (is.null(d) || !all(c("data_type", col, "observed") %in% names(d))) return(c(obs = NA, exp = NA, z = NA))
  y <- d[d$data_type == stream, , drop = FALSE]; if (!nrow(y)) return(c(obs = NA, exp = NA, z = NA))
  p <- suppressWarnings(as.numeric(y[[col]])); ok <- is.finite(p); p <- p[ok]
  obs <- sum(y$observed[ok] == k, na.rm = TRUE); ex <- sum(p); sd <- sqrt(sum(p * (1 - p)))
  c(obs = obs, exp = ex, z = if (is.finite(sd) && sd > 0) (obs - ex) / sd else NA_real_)
}

V <- list()
V1row <- function(stage, criterion, observed, threshold, verdict, why)
  V[[length(V) + 1]] <<- data.frame(stage = stage, criterion = criterion, observed = observed,
                                    threshold = threshold, verdict = verdict, why = why,
                                    stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# THE PIN. Every key that is NOT a declared delta is fixed here, so a rung differs from
# the control in exactly the lever it names. This is the B22 lesson from the pooled
# ladder: WINDOW there was missing four per-season keys, and nothing verified that the
# rungs were comparable until the preflight was written. A ladder whose rungs differ in
# an undeclared key is not a ladder, it is five runs.
# ---------------------------------------------------------------------------
WINDOW <- list(
  est_date_start = "2024-09-16", est_date_end = "2025-09-15", season_filter = "2024-25",
  pot_open_date = "2024-12-01", pot_closures = NULL, census_windows = NULL,
  pot_closure_start = "2024-09-16", pot_closure_end = "2024-11-30",
  commercial_opener = "2025-02-11",
  census_start_date = "2024-12-01", census_end_date = "2025-02-08",
  # the census and PE levers, pinned because they enter the PORT total the rungs report
  census_expansion = "none", census_uncertainty = "charter", charter_expansion = "vessel",
  charter_frame = "roster",
  pe_empty_effort_stratum = "local_day_type", pe_empty_stratum = "local",
  pe_variance = "impute_aware",
  # the turnover and f levers, pinned at the SHIPPED resolvers, not at resolved numbers:
  # resolving them here would freeze this window's values into the ladder and make the
  # rungs incomparable with any other run
  tau_shore_prior_mu = "derived", tau_boat_prior_mu = "calibration",
  shared_tau = TRUE, osp_scale_is_tau = TRUE,
  crab_fraction_strata = "month", crab_fraction_source = "both",
  # sampler: one seed, one iteration count, for every rung
  bss_seed = 20260619,
  # the AR machinery. ar_adaptive = FALSE is what makes gear_period_bss the lever at all;
  # ar_force NULL so nothing overrides it; ar_escalate FALSE so a rung cannot quietly
  # escalate to a resolution it was not asked for and report it as the one it was.
  ar_adaptive = FALSE, ar_force = NULL, ar_escalate = FALSE, ar_rung_adequacy = TRUE,
  # the pot-closure sub-season's period is NOT the lever and must not move
  gear_period_bss_pot_closure = "biweekly",
  # ZI: on for the pooled track only, which is the shipped state. G5 overrides.
  estimate_catch_zi = TRUE, catch_zi_populations = "shore", catch_zi_tracks = "pooled",
  zi_catch_prior_a = 1, zi_catch_prior_b = 9
)

# ---------------------------------------------------------------------------
# THE RUNGS. `delta` is the ONLY thing that may differ from the pinned configuration.
# ---------------------------------------------------------------------------
.pb <- function(ag) list(gear_period_bss = list(all_gear = ag, pot_closure = "biweekly"))

STAGE_DEFS <- list(
  G0 = list(tag = "GZ-G0-desk", fit = FALSE, res = NA_character_,
            item = "prerequisites; no fit",
            delta = list()),
  G1 = list(tag = "GZ-G1-month", fit = TRUE, res = "month",
            item = "D3 control: the SHIPPED period, and the bit-identity proof for the D6 Stan edit",
            delta = .pb("month")),
  G2 = list(tag = "GZ-G2-biweekly", fit = TRUE, res = "biweekly",
            item = "D3: biweekly",
            delta = .pb("biweekly")),
  G3 = list(tag = "GZ-G3-weekly", fit = TRUE, res = "weekly",
            item = "D3: weekly, the pooled track's resolution",
            delta = .pb("weekly")),
  G4 = list(tag = "GZ-G4-daily", fit = TRUE, res = "daily",
            item = "D3: daily, completing the bracket",
            delta = .pb("daily")),
  G5 = list(tag = "GZ-G5-zi", fit = TRUE, res = NA_character_,
            item = "D6: the zero-inflated shore catch likelihood ON for this track",
            delta = c(.pb("PLACEHOLDER"),
                      list(catch_zi_tracks = c("pooled", "gear_resolved"))))
)
# G5's resolution is resolved after ZI_AT so the pair is like for like.
.zi_res <- function() if (identical(ZI_AT, "d3")) STAGE_DEFS$G3$res else ZI_AT
STAGE_DEFS$G5$res <- .zi_res()
STAGE_DEFS$G5$delta$gear_period_bss <- list(all_gear = .zi_res(), pot_closure = "biweekly")
# the rung G5 is paired against: the AR rung at the same resolution, ZI off
.zi_control <- function() {
  r <- .zi_res()
  m <- c(month = "G1", monthly = "G1", biweekly = "G2", weekly = "G3", daily = "G4")
  m[[as.character(r)]] %||% "G3"
}

DELTA_KEYS <- unique(c("gear_period_bss", "catch_zi_tracks"))

resolve_cfg <- function(sid) {
  cfg <- BASE
  for (k in names(WINDOW)) if (!identical(k, "gear_period_bss_pot_closure")) cfg[[k]] <- WINDOW[[k]]
  cfg$gear_period_bss <- list(all_gear = "month",
                              pot_closure = WINDOW$gear_period_bss_pot_closure)
  d <- STAGE_DEFS[[sid]]$delta
  for (k in names(d)) cfg[[k]] <- d[[k]]
  cfg
}
# A digest over the delta keys and the pinned window, so RESUME cannot reuse a folder
# whose configuration is not this rung's. Same construction as the pooled ladder's
# stage_digest(), and for the same reason: reusing on the strength of a folder NAME
# silently mixed configurations into one ladder there on 2026-09-10.
stage_digest <- function(sid) {
  cfg <- resolve_cfg(sid)
  keys <- sort(unique(c(DELTA_KEYS, names(WINDOW)[names(WINDOW) != "gear_period_bss_pot_closure"])))
  txt <- paste(vapply(keys, function(k)
    paste0(k, "=", paste(format(unlist(cfg[[k]] %||% "NULL")), collapse = "|")), character(1)),
    collapse = ";")
  substr(digest_or_hash(txt), 1, 8)
}
# digest() is not a dependency of this repository, so hash with what is here. This is a
# CACHE KEY, not a security primitive: it only has to change when the configuration does.
digest_or_hash <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(x))
  v <- utf8ToInt(x); h <- 2166136261
  for (b in v) { h <- bitwXor(h, b %% 256); h <- (h * 16777619) %% 2^32 }
  sprintf("%08x%08x", h %% 2^32, (sum(v) * 2654435761) %% 2^32)
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
  sprintf("stan:%s drivers:%s fns:%s",
          .g(.here("02_stan_models")), .g(c(GEAR_RMD)), .g(.here("03_R_functions")))
}
.stage_stamp <- function(dir, sid) {
  writeLines(c(sprintf("stage: %s", sid),
               sprintf("resolution: %s", STAGE_DEFS[[sid]]$res %||% "NA"),
               sprintf("ZI_AT: %s", ZI_AT),
               sprintf("digest: %s", stage_digest(sid)),
               sprintf("code: %s", .code_fingerprint()),
               sprintf("rstan: %s / StanHeaders %s",
                       as.character(utils::packageVersion("rstan")),
                       as.character(utils::packageVersion("StanHeaders"))),
               sprintf("written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               "", "# run_gear_ar_zi_2026-09-13.R writes this after a rung renders. RESUME",
               "# reuses a folder ONLY when its digest matches the rung being requested."),
             file.path(dir, "GEAR_STAGE.txt"))
}
.stamp_field <- function(dir, key) {
  p <- file.path(dir %||% "", "GEAR_STAGE.txt")
  if (!file.exists(p)) return(NA_character_)
  l <- grep(paste0("^", key, ":"), readLines(p, warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_character_) else trimws(sub(paste0("^", key, ":"), "", l[1]))
}
find_outdir <- function(tag) {
  dirs <- list.dirs(.here("05_output"), recursive = FALSE)
  hits <- unlist(lapply(dirs, function(d)
    list.dirs(d, recursive = FALSE)[basename(list.dirs(d, recursive = FALSE)) == paste0(PREFIX, tag)]))
  if (!length(hits)) return(NA_character_)
  hits[order(basename(dirname(hits)), decreasing = TRUE)][1]
}

# ---------------------------------------------------------------------------
# G0: THE DESK STAGE. Everything that can be proven without MCMC is proven here, because
# the alternative is discovering it 3 h in. Nothing below fits anything.
# ---------------------------------------------------------------------------
stage_G0 <- function() {
  banner("G0  DESK: prerequisites, before any MCMC")

  # (1) the Stan file carries the ZI block, and the prep builds every variable it declares
  st <- readLines(GEAR_STAN, warn = FALSE)
  need <- c("zi_catch", "zi_catch_prior_a", "zi_catch_prior_b")
  decl <- all(vapply(need, function(v) any(grepl(paste0("(^|[^A-Za-z0-9_])", v, "\\s*;"), st)), logical(1)))
  V1row("G0", "the gear Stan declares the three ZI data variables",
        paste(need[vapply(need, function(v) any(grepl(paste0("(^|[^A-Za-z0-9_])", v, "\\s*;"), st)), logical(1))],
              collapse = ", "),
        "all three present", if (decl) "PASS" else "FAIL",
        "Without them nothing else in this file can settle D6.")
  V1row("G0", "theta_C is declared ZERO-SIZE when zi_catch = 0",
        paste(grep("vector<lower=0, upper=1>\\[zi_catch\\] theta_C", st, value = TRUE), collapse = ""),
        "vector<lower=0, upper=1>[zi_catch] theta_C",
        if (any(grepl("vector<lower=0, upper=1>\\[zi_catch\\] theta_C", st))) "PASS" else "FAIL",
        paste("This is the whole basis of the OFF path being bit-identical rather than merely",
              "similar: a zero-size parameter consumes no element of the unconstrained vector",
              "and no initialization draw. G1 proves it empirically against the committed R5 render."))
  V1row("G0", "the season total is scaled by zi_scale, so turning ZI on cannot inflate it",
        paste(length(grep("\\* zi_scale;", st)), "site(s)"), "1 site on lambda_Ctot_S",
        if (length(grep("f_crab\\[f_stratum\\[d\\]\\] \\* zi_scale;", st)) == 1) "PASS" else "FAIL",
        paste("lambda_C is fitted to the non-inflated component and rises to absorb the zeros",
              "theta_C removes, so an unscaled total would be inflated by 1/(1 - theta_C) purely",
              "by turning the feature on. The pooled model has done this since 2026-09-02."))
  if (exists("bss_stan_data_names")) {
    dn  <- bss_stan_data_names(GEAR_STAN)
    src <- paste(readLines(GEAR_PREP, warn = FALSE), collapse = "\n")
    miss <- dn[!vapply(dn, function(v)
      grepl(paste0("(^|[^A-Za-z0-9_.])", v, "([^A-Za-z0-9_]|$)"), src, perl = TRUE), logical(1))]
    V1row("G0", "every variable the gear Stan declares is built in prep_bss_crab_gear.R",
          sprintf("%d declared, %d missing", length(dn), length(miss)),
          "0 missing", if (!length(miss)) "PASS" else "FAIL",
          paste("The 2026-08-25 batch declared five variables that no prep built; Stan failed at",
                "data initialization on every fit, rstan returned an empty stanfit instead of",
                "raising, and every fit of every rung was lost to it. Missing:",
                if (length(miss)) paste(miss, collapse = ", ") else "none"))
  }

  # (2) COMPARABILITY. Each rung's resolved configuration must differ from G1's in
  #     declared keys only. This is what makes the ladder a ladder.
  base_cfg <- resolve_cfg("G1")
  for (sid in setdiff(names(STAGE_DEFS), c("G0", "G1"))) {
    cfg <- resolve_cfg(sid)
    ks  <- union(names(base_cfg), names(cfg))
    ks  <- setdiff(ks, c("run_tag", "model", "crabbing_holiday_dates", "opener_f_dates",
                         "razor_dig_dates", "crab_fraction_rows", "osp_crab_rows",
                         "tau_boat_prior_source", "tau_boat_prior_calibration_table"))
    diff <- ks[!vapply(ks, function(k)
      identical(format(unlist(base_cfg[[k]] %||% "NULL")), format(unlist(cfg[[k]] %||% "NULL"))),
      logical(1))]
    undeclared <- setdiff(diff, DELTA_KEYS)
    V1row(sid, "differs from G1 in DECLARED keys only",
          sprintf("differs in: %s", if (length(diff)) paste(diff, collapse = ", ") else "nothing"),
          sprintf("a subset of {%s}", paste(DELTA_KEYS, collapse = ", ")),
          if (!length(undeclared)) "PASS" else "FAIL",
          paste("A rung that differs in an undeclared key measures that key as well as its own,",
                "and nothing in the output would say so. Undeclared here:",
                if (length(undeclared)) paste(undeclared, collapse = ", ") else "none"))
  }

  # (3) each fitted rung's requested resolution must actually reach build_subseasons()
  for (sid in names(STAGE_DEFS)) {
    st_ <- STAGE_DEFS[[sid]]; if (!isTRUE(st_$fit)) next
    cfg <- resolve_cfg(sid)
    ss  <- tryCatch(build_subseasons(cfg), error = function(e) NULL)
    got <- if (is.null(ss)) NA_character_ else {
      ag <- Filter(function(x) identical(x$gear_regime, "all_gear"), ss)
      if (length(ag)) as.character(ag[[1]]$period_bss) else NA_character_
    }
    pc <- if (is.null(ss)) NA_character_ else {
      cl <- Filter(function(x) identical(x$gear_regime, "pot_closure"), ss)
      if (length(cl)) as.character(cl[[1]]$period_bss) else NA_character_
    }
    V1row(sid, "the requested all-gear period reaches build_subseasons(), and pot closure does not move",
          sprintf("all_gear = %s, pot_closure = %s", got, pc),
          sprintf("all_gear = %s, pot_closure = biweekly", st_$res),
          if (identical(got, st_$res) && identical(pc, "biweekly")) "PASS" else "FAIL",
          paste("gear_period_bss is read inside build_subseasons() and handed to the driver as",
                "fixed_resolution. If it did not arrive, the rung would fit at the default and",
                "report itself as the resolution it asked for, which is the one failure mode",
                "of this ladder that would not announce itself."))
  }

  # (4) the ZI pair is like for like in everything but ZI
  zc <- .zi_control()
  a <- resolve_cfg(zc); b <- resolve_cfg("G5")
  ks <- setdiff(union(names(a), names(b)), c("run_tag", "model"))
  diff <- ks[!vapply(ks, function(k)
    identical(format(unlist(a[[k]] %||% "NULL")), format(unlist(b[[k]] %||% "NULL"))), logical(1))]
  V1row("G5", sprintf("the D6 pair (%s vs G5) differs in catch_zi_tracks ONLY", zc),
        sprintf("differs in: %s", if (length(diff)) paste(diff, collapse = ", ") else "nothing"),
        "catch_zi_tracks", if (identical(diff, "catch_zi_tracks")) "PASS" else "FAIL",
        paste("A paired elpd is only a paired elpd if the two fits differ in one thing. The pooled",
              "adoption's Z1/Z0 pair was compared at a COMMON weekly resolution for exactly this",
              "reason, after an earlier comparison at different resolutions had to be retracted."))

  # (5) the reference folders the verdicts read must exist
  for (r in c(R5 = REF_R5, pooled = REF_POOLED)) {
    V1row("G0", sprintf("the reference render %s is present", basename(r)),
          if (dir.exists(.here("05_output", r))) "present" else "MISSING",
          "present", if (dir.exists(.here("05_output", r))) "PASS" else "FAIL",
          paste("G1's bit-identity verdict reads the R5 folder and the cross-track columns read",
                "the pooled R4 folder. Without them those verdicts are skipped, not wrong, but",
                "the D6 Stan edit then has no inertness proof."))
  }

  # (6) the shipped configuration must be UNCHANGED by this patch, which is the promise
  #     that landing the Stan port does not move any committed number
  V1row("G0", "the SHIPPED config still fits the gear track without ZI",
        sprintf("catch_zi_tracks = %s", paste(BASE$catch_zi_tracks %||% "NULL", collapse = ",")),
        "pooled", if (identical(as.character(BASE$catch_zi_tracks %||% ""), "pooled")) "PASS" else "FAIL",
        paste("The D6 Stan port is inert until this key names this track. If it shipped as both,",
              "the next production gear render would silently change and the committed R5",
              "cross-check figure quoted in PULL_REQUEST.md would go stale without anyone",
              "touching a number."))
  V1row("G0", "the SHIPPED gear all-gear period is still month",
        sprintf("gear_period_bss$all_gear = %s", BASE$gear_period_bss$all_gear %||% "NULL"),
        "month", if (identical(as.character(BASE$gear_period_bss$all_gear %||% ""), "month")) "PASS" else "FAIL",
        "Same promise for D3: the lever is exposed, not moved.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# RUNNING ONE RUNG
# ---------------------------------------------------------------------------
run_one <- function(sid) {
  st <- STAGE_DEFS[[sid]]
  rule()
  cat(sprintf("  %-3s %-16s  %s\n", sid, st$tag, st$item))
  cat(sprintf("       period_bss$all_gear = %-9s  catch_zi_tracks = %s  digest %s\n",
              st$res %||% "NA",
              paste(resolve_cfg(sid)$catch_zi_tracks, collapse = "+"), stage_digest(sid)))
  existing <- find_outdir(st$tag)
  if (isTRUE(RESUME) && !is.na(existing) && file.exists(file.path(existing, "run_parameters.txt"))) {
    dg <- .stamp_field(existing, "digest")
    if (identical(dg, stage_digest(sid))) {
      cd <- .stamp_field(existing, "code")
      cat("  RESUME: output present at", basename(existing), "with a MATCHING digest - skipping the fit.\n")
      if (!is.na(cd) && !identical(cd, .code_fingerprint())) {
        cat(sprintf("          *** THE CODE HAS CHANGED SINCE THAT FIT.\n              recorded %s\n              now      %s\n",
                    cd, .code_fingerprint()))
        V1row(sid, "the reused fit was produced by DIFFERENT code than this run",
              sprintf("recorded %s; now %s", cd, .code_fingerprint()),
              "the fingerprints match", "REVIEW",
              paste("RESUME matched the CONFIG digest and reused the fit, because re-fitting over",
                    "a code edit that may be a comment is worse than the problem. But any",
                    "cross-rung claim involving this rung is not interpretable; either accept it",
                    "as indicative or delete the folder and re-fit."))
      }
      return(existing)
    }
    cat(sprintf("  RESUME: %s exists but its digest %s does not match (%s). RE-RUNNING.\n",
                basename(existing), if (is.na(dg)) "is ABSENT" else dg, stage_digest(sid)))
  }
  if (isTRUE(DRY_RUN)) { cat("  DRY_RUN: not fitting.\n"); return(NA_character_) }
  cfg <- resolve_cfg(sid); cfg$model <- "gear_resolved"; cfg$run_tag <- st$tag
  run_env <- new.env(parent = globalenv()); run_env$run_config <- cfg
  t0 <- Sys.time()
  html <- rmarkdown::render(GEAR_RMD, envir = run_env, quiet = FALSE)
  od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e) NA_character_)
  if (!is.na(od) && dir.exists(od) && file.exists(html) &&
      normalizePath(dirname(html)) != normalizePath(od)) {
    if (isTRUE(file.copy(html, file.path(od, basename(html)), overwrite = TRUE)))
      suppressWarnings(file.remove(html))
    else cat("  WARNING: could not move the rendered HTML; the NEXT rung will overwrite it.\n")
  }
  done <- find_outdir(st$tag)
  cat(sprintf("  %s finished in %.1f min -> %s\n", sid,
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              if (is.na(done)) "OUTPUT FOLDER NOT FOUND" else basename(done)))
  if (is.na(done)) stop(sprintf("%s rendered but no folder named %s%s exists under 05_output",
                                sid, PREFIX, st$tag))
  .stage_stamp(done, sid)
  done
}

# ---------------------------------------------------------------------------
# THE LADDER TABLE. One row per rung, carrying every column the two decision rules read,
# so the recommendation is computed from a table a reader can check rather than from
# variables that only existed inside this script.
# ---------------------------------------------------------------------------
LAD <- list()
ladder_row <- function(sid, dir) {
  if (is.na(dir %||% NA) || !dir.exists(dir %||% "")) return(invisible(NULL))
  pt <- .port_row(rd(dir, "port_total_Dungeness_Kept.csv"))
  aq <- .adq(dir); al <- .arlog(dir); gt <- .gate(dir)
  b0 <- .bin(dir, 0L); b1 <- .bin(dir, 1L)
  th <- .full_row(dir, "shore_all_gear", "theta_C_out")
  nobs <- { d <- rd(dir, "ppc_byobs_shore_all_gear_Dungeness_Kept.csv")
            if (is.null(d)) NA_integer_ else sum(d$data_type == "catch", na.rm = TRUE) }
  LAD[[sid]] <<- data.frame(
    rung = sid, folder = basename(dir),
    requested_res = STAGE_DEFS[[sid]]$res %||% NA_character_,
    fitted_res = as.character(al$ar_resolution %||% NA), P_n = .num1(al$P_n),
    zi_tracks = paste(resolve_cfg(sid)$catch_zi_tracks, collapse = "+"),
    theta_C = .num1(th$mean),
    shore_pc = .comp(dir, "shore (Pot closure)"), shore_ag = .comp(dir, "shore (All gear)"),
    boat_pc = .comp(dir, "private_boat (Pot closure)"), boat_ag = .comp(dir, "private_boat (All gear)"),
    census = .comp(dir, "comm_charter (census)", "PE_catch"),
    port = .num1(pt$BSS_median), port_lo95 = .num1(pt$BSS_lo95), port_hi95 = .num1(pt$BSS_hi95),
    n_fits_bss = if (is.null(gt)) NA_integer_ else sum(gt$method_selected == "BSS", na.rm = TRUE),
    n_fits = if (is.null(gt)) NA_integer_ else nrow(gt),
    gate_all_pass = if (is.null(gt)) NA else all(as.logical(gt$pass_convergence), na.rm = TRUE),
    p_loo_frac = .num1(aq$p_loo_frac), p_loo_stream = as.character(aq$p_loo_worst_stream %||% NA),
    n_pareto_bad = .num1(aq$n_pareto_bad), catch_n_obs = nobs,
    pit_worst_bias = .num1(aq$pit_worst_bias), cov50_worst_dev = .num1(aq$cov50_worst_dev),
    disp_neff_min = .num1(aq$disp_neff_min),
    divergences = .num1(al$divergences), divergence_fraction = .num1(al$divergence_fraction),
    zero_obs = b0[["obs"]], zero_exp = b0[["exp"]], zero_z = b0[["z"]],
    one_obs = b1[["obs"]], one_exp = b1[["exp"]], one_z = b1[["z"]],
    stringsAsFactors = FALSE)
}
.pct <- function(a, b) if (isTRUE(is.finite(a)) && isTRUE(is.finite(b)) && b != 0) 100 * (a - b) / b else NA_real_

# ---------------------------------------------------------------------------
# PER-RUNG VERDICTS
# ---------------------------------------------------------------------------

# G1 is the control and its verdict is the one that gates every other reading in this
# file, so it is reported first and in full.
verdict_G1 <- function(dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  ref <- .here("05_output", REF_R5)
  if (!dir.exists(ref)) {
    V1row("G1", "bit-identity against the pre-D6 gear render", "reference folder absent",
          REF_R5, "REVIEW",
          paste("Without it the D6 Stan edit has no empirical inertness proof, only the",
                "structural argument that a zero-size parameter cannot change the",
                "unconstrained vector. Do not read the other rungs as if the edit were proven",
                "inert."))
    return(invisible(NULL))
  }
  ex <- tryCatch(fit_exactness(dir, ref, what = "the gear fits"), error = function(e) NULL)
  ok <- !is.null(ex) && isTRUE(ex$identical %||% FALSE)
  V1row("G1", "THE D6 STAN EDIT IS INERT WHEN OFF: G1 is bit-identical to the pre-edit render",
        if (is.null(ex)) "could not compare" else (ex$text %||% paste0("identical = ", ex$identical)),
        sprintf("every shared parameter row identical to %s", basename(REF_R5)),
        if (ok) "PASS" else "FAIL",
        paste("G1 ships the same configuration as the committed R5 cross-check and differs from",
              "it only in carrying a Stan file with the ZI block compiled in but switched off.",
              "theta_C is declared vector[zi_catch], so at zi_catch = 0 it is zero-size and",
              "consumes neither a parameter nor an initialization draw. If this FAILS, the edit",
              "is not inert, the 93,274 R5 figure in PULL_REQUEST.md is no longer reproducible,",
              "and EVERY rung in this ladder is measuring the edit as well as its own lever."))
  p  <- .num1(.port_row(rd(dir, "port_total_Dungeness_Kept.csv"))$BSS_median)
  pr <- .num1(.port_row(rd(ref, "port_total_Dungeness_Kept.csv"))$BSS_median)
  V1row("G1", "the port total reproduces the committed R5 figure",
        sprintf("G1 %s vs R5 %s (%.4f%%)", fmt(p, 0), fmt(pr, 0), .pct(p, pr)),
        "0.0000%", if (isTRUE(abs(.pct(p, pr)) < 1e-6)) "PASS" else "FAIL",
        "A weaker restatement of the row above, in the unit a reader will quote.")
}

# Every fitted rung: did it get the resolution it asked for, and did it pass the gate?
verdict_rung <- function(sid, dir) {
  if (is.na(dir %||% NA)) return(invisible(NULL))
  L <- LAD[[sid]]; if (is.null(L)) return(invisible(NULL))
  want <- STAGE_DEFS[[sid]]$res
  # the driver normalizes "month" to "monthly" on the way through bss_ar_resolution()
  norm <- function(x) sub("^month$", "monthly", as.character(x))
  V1row(sid, "the fit used the resolution the rung asked for",
        sprintf("requested %s, fitted %s (P_n = %s)", want, L$fitted_res, fmt(L$P_n, 0)),
        sprintf("fitted = %s", norm(want)),
        if (identical(norm(L$fitted_res), norm(want))) "PASS" else "FAIL",
        paste("Read from ar_escalation_log.csv, which records what the sampler was actually",
              "given. A rung that fell back would otherwise report itself as the resolution it",
              "requested and sit in the ladder as a duplicate of another rung."))
  V1row(sid, "every fit in the rung passes the convergence gate",
        sprintf("%s of %s fits report BSS; gate all-pass = %s",
                fmt(L$n_fits_bss, 0), fmt(L$n_fits, 0), L$gate_all_pass),
        "all fits pass", if (isTRUE(L$gate_all_pass)) "PASS" else "FAIL",
        paste("Decision rule 1. A rung whose fits do not sample is not a candidate resolution,",
              "whatever its adequacy statistics say: the gate is the single authority on method",
              "selection and a PE fallback inside a rung makes its port total a different",
              "estimator from the other rungs'."))
  bad_share <- if (isTRUE(is.finite(L$n_pareto_bad)) && isTRUE(is.finite(L$catch_n_obs)) && L$catch_n_obs > 0)
    L$n_pareto_bad / L$catch_n_obs else NA_real_
  V1row(sid, "adequacy: p_loo is a small fraction of n_obs and few Pareto k are bad",
        sprintf("p_loo_frac %s (worst stream %s); %s bad k of %s obs (%s%%); PIT bias %s; cov50 dev %s; min disp n_eff %s",
                fmt(L$p_loo_frac, 4), L$p_loo_stream, fmt(L$n_pareto_bad, 0), fmt(L$catch_n_obs, 0),
                fmt(100 * bad_share, 1), fmt(L$pit_worst_bias, 4), fmt(L$cov50_worst_dev, 4),
                fmt(L$disp_neff_min, 0)),
        "p_loo_frac <= 0.15 and bad k <= 5% of n_obs",
        if (isTRUE(L$p_loo_frac <= 0.15) && (is.na(bad_share) || isTRUE(bad_share <= 0.05))) "PASS" else "FAIL",
        paste("Decision rule 2, calibrated on the pooled ladder: the daily fit that was rejected",
              "as overfitted sat at p_loo 0.352 with 41 bad k; the weekly fit that replaced it at",
              "0.095 with 0; the gear R5 monthly fit at 0.0368 with 0. p_loo far above the number",
              "of parameters the model actually has is the signature of a latent process being",
              "fitted to individual observations."))
  # reported, never a criterion
  pl <- .num1(.comp(.here("05_output", REF_POOLED), "shore (All gear)"))
  V1row(sid, "REPORTED, NOT A CRITERION: the cross-track shore all-gear gap",
        sprintf("gear %s vs pooled R4 %s (%s%%); port %s vs %s (%s%%)",
                fmt(L$shore_ag, 0), fmt(pl, 0), fmt(.pct(L$shore_ag, pl), 2),
                fmt(L$port, 0),
                fmt(.num1(.port_row(rd(.here("05_output", REF_POOLED), "port_total_Dungeness_Kept.csv"))$BSS_median), 0),
                fmt(.pct(L$port, .num1(.port_row(rd(.here("05_output", REF_POOLED), "port_total_Dungeness_Kept.csv"))$BSS_median)), 2)),
        "no threshold; see decision rule 5", "INFO",
        paste("Decision rule 5. Agreement between the two tracks at a COMMON resolution is a",
              "corroboration of the fishery model; choosing this track's resolution so as to",
              "minimize the gap would tune one estimate to another estimate instead of to the",
              "data, and would make the corroboration circular. It is printed so the reader can",
              "see the gap close, and excluded from the recommendation by construction."))
}

# D6: the paired comparison, at a common resolution
verdict_G5 <- function(dir_zi, dir_ctl, ctl_name) {
  if (is.na(dir_zi %||% NA) || is.na(dir_ctl %||% NA)) return(invisible(NULL))
  lab <- "shore_all_gear_Dungeness_Kept"
  pa <- file.path(dir_ctl, sprintf("loo_pointwise_catch_%s.csv", lab))
  pb <- file.path(dir_zi,  sprintf("loo_pointwise_catch_%s.csv", lab))
  el <- tryCatch(loo_elpd_paired(pa, pb, label = sprintf("%s -> G5 (ZI on)", ctl_name)),
                 error = function(e) NULL)
  if (is.null(el)) {
    V1row("G5", "paired elpd on the shore all-gear catch stream", "could not compute",
          "gain > 2 paired SE", "REVIEW",
          paste("loo_pointwise_catch_*.csv is missing or the two runs' observation vectors do",
                "not align. Without it decision rule 2 cannot be applied and D6 stays open."))
  } else {
    gain <- el$diff %||% NA_real_; se <- el$se %||% NA_real_
    ratio <- if (isTRUE(is.finite(gain)) && isTRUE(is.finite(se)) && se > 0) gain / se else NA_real_
    V1row("G5", "D6 rule 2: the paired elpd gain exceeds 2 paired SE",
          sprintf("%s nats at %s paired SE (%s SE)", fmt(gain, 1), fmt(se, 2), fmt(ratio, 2)),
          "> 2 paired SE", if (isTRUE(ratio > 2)) "PASS" else if (isTRUE(ratio > 0)) "FAIL" else "FAIL",
          paste("Paired, not naive: on the pooled Z1/Z0 pair the paired SE was 5.5 nats against",
                "a naive 46.7, which is the difference between 2.69 SE and 0.32 SE. The pooled",
                "adoption earned +11.6 nats at 2.30 SE."))
    z <- el$zero %||% NULL; p <- el$positive %||% NULL
    V1row("G5", "D6 rule 4: the gain is not bought entirely at the zeros",
          if (is.null(z) || is.null(p)) "decomposition unavailable" else
            sprintf("zeros %s nats (n = %s), positives %s nats (n = %s)",
                    fmt(z[["diff"]], 1), fmt(z[["n"]], 0), fmt(p[["diff"]], 1), fmt(p[["n"]], 0)),
          "the positive-count loss is smaller than the zero-count gain",
          if (is.null(z) || is.null(p)) "REVIEW"
          else if (isTRUE(is.finite(z[["diff"]])) && isTRUE(is.finite(p[["diff"]])) &&
                   (p[["diff"]] > 0 || abs(p[["diff"]]) < z[["diff"]])) "PASS" else "FAIL",
          paste("A zero-inflation parameter that earns its elpd purely at the zeros while",
                "degrading the positive counts is doing something different from one that",
                "improves the fit everywhere, and the aggregate number cannot tell them apart.",
                "On the pooled pair the split was +25.2 at the zeros and -10.5 at the positives:",
                "net positive, but bought."))
  }
  La <- LAD[[ctl_name]]; Lb <- LAD[["G5"]]
  if (!is.null(La) && !is.null(Lb)) {
    V1row("G5", "D6 rule 3: both count bins improve",
          sprintf("zero z %s -> %s ; one z %s -> %s",
                  fmt(La$zero_z, 2), fmt(Lb$zero_z, 2), fmt(La$one_z, 2), fmt(Lb$one_z, 2)),
          "|z| falls in BOTH bins",
          if (isTRUE(abs(Lb$zero_z) < abs(La$zero_z)) && isTRUE(abs(Lb$one_z) < abs(La$one_z))) "PASS"
          else if (isTRUE(abs(Lb$zero_z) < abs(La$zero_z)) || isTRUE(abs(Lb$one_z) < abs(La$one_z))) "REVIEW"
          else "FAIL",
          paste("The pooled adoption halved both (zero +3.7 -> +2.0, one -6.1 -> -3.3). A feature",
                "that fixes the zero bin by breaking the one bin is moving the misfit, not",
                "removing it, which is what an earlier reading of the prototype got wrong and had",
                "to retract."))
    V1row("G5", "theta_C is identified and not pinned at its prior",
          sprintf("theta_C = %s (Beta(%s, %s) prior mean %s)", fmt(Lb$theta_C, 4),
                  fmt(BASE$zi_catch_prior_a, 0), fmt(BASE$zi_catch_prior_b, 0),
                  fmt(BASE$zi_catch_prior_a / (BASE$zi_catch_prior_a + BASE$zi_catch_prior_b), 3)),
          "distinguishable from the prior mean", if (isTRUE(is.finite(Lb$theta_C))) "INFO" else "REVIEW",
          paste("A theta_C sitting on its prior mean means the data did not inform it and the",
                "elpd gain, if any, came from somewhere else. kappa_OSP is the standing example",
                "of a parameter whose reported median is its prior."))
    V1row("G5", "REPORTED, NOT A CRITERION: the port total under ZI",
          sprintf("%s -> %s (%s%%); theta_C %s implies zi_scale %s",
                  fmt(La$port, 0), fmt(Lb$port, 0), fmt(.pct(Lb$port, La$port), 2),
                  fmt(Lb$theta_C, 4), fmt(1 - Lb$theta_C, 4)),
          "no threshold; see decision rule 5", "INFO",
          paste("The season total is multiplied by (1 - theta_C) in generated quantities precisely",
                "so that turning the feature on does not inflate it. A total that barely moves is",
                "the design working, not weak evidence, and a total that moves a lot would be the",
                "thing to investigate."))
  }
}

# ---------------------------------------------------------------------------
# THE RECOMMENDATION. Computed from the ladder table by applying the rule stated in the
# header, in that order, and printing which clause did the work. It is a RECOMMENDATION:
# it names the adoption edit and stops, because adopting is Matt's call.
# ---------------------------------------------------------------------------
REC <- list()
recommend <- function() {
  banner("RECOMMENDATION: D3 (the gear AR period) and D6 (the gear ZINB)")
  if (!length(LAD)) { cat("  No fitted rungs to read. Nothing recommended.\n"); return(invisible(NULL)) }
  T <- do.call(rbind, LAD)
  ar <- T[T$rung %in% c("G1", "G2", "G3", "G4"), , drop = FALSE]

  cat("\n  D3. THE RULE, APPLIED IN ORDER.\n")
  if (!nrow(ar)) { cat("     no AR rungs fitted.\n") } else {
    ar$bad_share <- ifelse(is.finite(ar$n_pareto_bad) & is.finite(ar$catch_n_obs) & ar$catch_n_obs > 0,
                           ar$n_pareto_bad / ar$catch_n_obs, NA_real_)
    ar$eligible  <- !is.na(ar$gate_all_pass) & ar$gate_all_pass
    ar$adequate  <- ar$eligible & is.finite(ar$p_loo_frac) & ar$p_loo_frac <= 0.15 &
                    (is.na(ar$bad_share) | ar$bad_share <= 0.05)
    rank <- c(monthly = 1L, month = 1L, biweekly = 2L, weekly = 3L, daily = 4L)
    ar$fineness <- rank[as.character(ar$requested_res)]
    cat(sprintf("     %-4s %-9s %5s  %-6s %-8s %9s %9s %8s %8s\n",
                "rung", "period", "P_n", "gate", "adequate", "p_loo", "bad k %", "shore_ag", "port"))
    for (i in order(ar$fineness)) cat(sprintf("     %-4s %-9s %5s  %-6s %-8s %9s %9s %8s %8s\n",
      ar$rung[i], ar$requested_res[i], fmt(ar$P_n[i], 0),
      if (isTRUE(ar$eligible[i])) "pass" else "FAIL",
      if (isTRUE(ar$adequate[i])) "yes" else "no",
      fmt(ar$p_loo_frac[i], 4), fmt(100 * ar$bad_share[i], 2),
      fmt(ar$shore_ag[i], 0), fmt(ar$port[i], 0)))
    surv <- ar[isTRUE(TRUE) & ar$adequate %in% TRUE, , drop = FALSE]
    if (!nrow(surv)) {
      cat("\n     NOTHING SURVIVES clauses 1 and 2. D3 STAYS OPEN.\n")
      cat("     That is a real answer, not a failed run: it would mean no resolution in the\n")
      cat("     bracket is both identified and not overfitted, and the next question is the\n")
      cat("     likelihood rather than the period.\n")
      REC$d3 <<- list(pick = NA_character_, why = "no rung passed the gate and the adequacy clauses")
    } else {
      pick <- surv[which.max(surv$fineness), , drop = FALSE]
      coarser <- ar[ar$fineness < pick$fineness & ar$adequate %in% TRUE, , drop = FALSE]
      cat(sprintf("\n     PICK: %s (%s, P_n = %s), by clause 4: the FINEST resolution surviving\n",
                  pick$rung, pick$requested_res, fmt(pick$P_n, 0)))
      cat("     clauses 1 and 2.\n")
      cat(sprintf("     Shore all-gear %s, port %s [%s, %s].\n",
                  fmt(pick$shore_ag, 0), fmt(pick$port, 0), fmt(pick$port_lo95, 0), fmt(pick$port_hi95, 0)))
      shipped <- ar[ar$rung == "G1", , drop = FALSE]
      if (nrow(shipped) && !identical(pick$rung, "G1"))
        cat(sprintf("     AGAINST THE SHIPPED PERIOD (%s): shore all-gear %s -> %s (%s%%), port %s -> %s (%s%%).\n",
                    shipped$requested_res, fmt(shipped$shore_ag, 0), fmt(pick$shore_ag, 0),
                    fmt(.pct(pick$shore_ag, shipped$shore_ag), 2),
                    fmt(shipped$port, 0), fmt(pick$port, 0), fmt(.pct(pick$port, shipped$port), 2)))
      if (identical(pick$rung, "G1"))
        cat("     This is the SHIPPED period. D3 closes with no change, and the -1.39% cross-track\n     gap is then a property of the two models, not a defect in this one.\n")
      else
        cat(sprintf("     ADOPTION EDIT, if you take it: run_config.R gear_period_bss$all_gear <- \"%s\".\n     One key. ar_max_resolution$gear_resolved stays dormant and is not the lever.\n",
                    pick$requested_res))
      if (nrow(coarser))
        cat(sprintf("     Also surviving, coarser: %s. Clause 4 prefers the finest; if you would rather\n     be conservative, %s is defensible and the reason is judgement, not evidence.\n",
                    paste(sprintf("%s (%s)", coarser$rung, coarser$requested_res), collapse = ", "),
                    coarser$rung[which.max(coarser$fineness)]))
      REC$d3 <<- list(pick = as.character(pick$requested_res),
                      why = sprintf("finest of {%s} surviving the gate and adequacy clauses",
                                    paste(surv$requested_res, collapse = ",")))
    }
    cat("\n     NOT USED IN THIS PICK, deliberately: elpd (clause 3) and the cross-track gap\n")
    cat("     (clause 5). Both are in the verdict table.\n")
  }

  cat("\n  D6. THE RULE, APPLIED IN ORDER.\n")
  zi <- T[T$rung == "G5", , drop = FALSE]
  ctl <- .zi_control(); ct <- T[T$rung == ctl, , drop = FALSE]
  if (!nrow(zi) || !nrow(ct)) {
    cat("     the ZI rung or its control was not fitted. D6 STAYS OPEN.\n")
    REC$d6 <<- list(pick = NA_character_, why = "the pair was not fitted")
  } else {
    vv <- do.call(rbind, V)
    g5 <- vv[vv$stage == "G5" & grepl("^D6 rule", vv$criterion), , drop = FALSE]
    for (i in seq_len(nrow(g5)))
      cat(sprintf("     %-6s %s\n              %s\n", g5$verdict[i], g5$criterion[i], g5$observed[i]))
    pass <- nrow(g5) > 0 && all(g5$verdict == "PASS")
    any_fail <- nrow(g5) > 0 && any(g5$verdict == "FAIL")
    cat(sprintf("\n     %s\n", if (pass)
      sprintf("ADOPT for this track: every D6 clause passes at %s resolution.", ct$requested_res)
      else if (any_fail)
        "DO NOT ADOPT for this track: at least one D6 clause fails."
      else "D6 STAYS OPEN: no clause fails outright but not all pass."))
    if (pass)
      cat("     ADOPTION EDIT: run_config.R catch_zi_tracks <- c(\"pooled\", \"gear_resolved\").\n     Note that this MOVES the committed R5 cross-check figure, so PULL_REQUEST.md and\n     the register's -1.17% gap have to be re-derived from a fresh gear render.\n")
    else
      cat("     The Stan port stays in the tree either way. It costs nothing when off (G1's\n     bit-identity verdict above is the proof), and it removes a real asymmetry: before\n     it, estimate_catch_zi = TRUE was read and silently ignored on this track.\n")
    REC$d6 <<- list(pick = if (pass) "adopt" else if (any_fail) "do not adopt" else "open",
                    why = sprintf("%d of %d D6 clauses pass at %s",
                                  sum(g5$verdict == "PASS"), nrow(g5), ct$requested_res))
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
banner(sprintf("GEAR AR PERIOD (D3) AND GEAR ZINB (D6)  |  DRY_RUN = %s  |  ZI_AT = %s",
               DRY_RUN, ZI_AT))
cat(sprintf("  stages: %s\n", paste(STAGES, collapse = ", ")))
cat(sprintf("  rstan %s / StanHeaders %s\n",
            as.character(utils::packageVersion("rstan")),
            as.character(utils::packageVersion("StanHeaders"))))
cat(sprintf("  code   %s\n", .code_fingerprint()))
cat(sprintf("  D6 pair: G5 (ZI on) against %s at %s resolution\n", .zi_control(), .zi_res()))
if (isTRUE(DRY_RUN))
  cat("\n  DRY RUN. G0 runs; no rung is fitted. Set DRY_RUN <- FALSE and source again.\n")

DIRS <- list()
# Every rung is wrapped, so an error inside one render (a Stan failure, a full disk, one
# bad interview row) costs that rung and not the ones already fitted. The pooled ladder
# lost an entire batch to an unwrapped verdict block on 2026-09-04.
for (sid in STAGES) {
  if (identical(sid, "G0")) { tryCatch(stage_G0(), error = function(e)
    V1row("G0", "the desk stage completed", conditionMessage(e), "no error", "ERROR",
          "G0 is pure desk work; an error here is a bug in this file, not a run failure.")); next }
  DIRS[[sid]] <- tryCatch(run_one(sid), error = function(e) {
    V1row(sid, "the rung completed", conditionMessage(e), "the render runs to completion", "ERROR",
          paste("The render errored. The other rungs are separate renders and are unaffected;",
                "re-running with RESUME = TRUE will skip whatever finished and retry this one."))
    NA_character_ })
  if (!is.na(DIRS[[sid]] %||% NA)) tryCatch(ladder_row(sid, DIRS[[sid]]), error = function(e)
    V1row(sid, "the ladder row was extracted", conditionMessage(e), "no error", "ERROR",
          "The fit survived; only the extraction failed. The folder is intact."))
}
for (sid in intersect(STAGES, c("G1", "G2", "G3", "G4", "G5")))
  tryCatch(verdict_rung(sid, DIRS[[sid]] %||% NA_character_), error = function(e)
    V1row(sid, "the rung verdicts were computed", conditionMessage(e), "no error", "ERROR", ""))
if ("G1" %in% STAGES) tryCatch(verdict_G1(DIRS$G1 %||% NA_character_), error = function(e)
  V1row("G1", "the bit-identity verdict was computed", conditionMessage(e), "no error", "ERROR", ""))
if ("G5" %in% STAGES) tryCatch(
  verdict_G5(DIRS$G5 %||% NA_character_, DIRS[[.zi_control()]] %||% NA_character_, .zi_control()),
  error = function(e) V1row("G5", "the D6 verdicts were computed", conditionMessage(e), "no error", "ERROR", ""))

if (length(LAD)) {
  banner("THE LADDER")
  T <- do.call(rbind, LAD)
  print(T[, c("rung", "requested_res", "fitted_res", "P_n", "zi_tracks", "shore_ag", "port",
              "p_loo_frac", "n_pareto_bad", "zero_z", "one_z", "theta_C")], row.names = FALSE)
}
if (length(V)) {
  banner("VERDICTS")
  VV <- do.call(rbind, V)
  for (i in seq_len(nrow(VV)))
    cat(sprintf("  %-6s %-4s %s\n         observed : %s\n         threshold: %s\n",
                VV$verdict[i], VV$stage[i], VV$criterion[i], VV$observed[i], VV$threshold[i]))
  cat(sprintf("\n  %d PASS  %d FAIL  %d REVIEW  %d INFO  %d ERROR\n",
              sum(VV$verdict == "PASS"), sum(VV$verdict == "FAIL"),
              sum(VV$verdict == "REVIEW"), sum(VV$verdict == "INFO"), sum(VV$verdict == "ERROR")))
}
tryCatch(recommend(), error = function(e) cat("  recommendation failed:", conditionMessage(e), "\n"))

if (!isTRUE(DRY_RUN)) {
  op <- .here("05_output")
  if (length(LAD)) merge_csv_by(do.call(rbind, LAD),
                                file.path(op, "gear_ar_zi_2026-09-13_ladder.csv"), key = "rung")
  if (length(V))   merge_csv_by(do.call(rbind, V),
                                file.path(op, "gear_ar_zi_2026-09-13_verdicts.csv"),
                                key = c("stage", "criterion"))
  rec <- data.frame(item = c("D3", "D6"),
                    recommendation = c(REC$d3$pick %||% NA_character_, REC$d6$pick %||% NA_character_),
                    why = c(REC$d3$why %||% NA_character_, REC$d6$why %||% NA_character_),
                    zi_at = ZI_AT, written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    code = .code_fingerprint(), stringsAsFactors = FALSE)
  merge_csv_by(rec, file.path(op, "gear_ar_zi_2026-09-13_recommendation.csv"), key = "item")
  cat("\n  wrote gear_ar_zi_2026-09-13_{ladder,verdicts,recommendation}.csv to 05_output/\n")
} else {
  cat("\n  DRY_RUN: nothing written. Set DRY_RUN <- FALSE and source again to start.\n")
}
