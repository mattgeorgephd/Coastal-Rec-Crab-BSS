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
# desk_advisory_block_cv_2026-11-19.R
#
# THE HELD-OUT-WEEK TEST WITH THE DAY DEVIATION INTEGRATED OUT (CHANGE_REGISTER B80's
# addendum R3a, and the prototype of the repair D54 asks for). Fits nothing: it reads
# posterior draws that B80's renders already sampled.
#
# WHY THIS EXISTS. The block CV (bss_block_cv.R, B39) approximates the posterior without a
# week by importance sampling on that week's joint likelihood. Since the per-day effort
# effect (A39) that likelihood is CONDITIONAL on one deviation uE per sampled day, and a
# held-out week's counts are exactly what identify its days' deviations: under the full
# posterior they are tight, under the leave-out posterior they are their prior, and the
# importance ratios between the two are heavy-tailed. On the authoritative run B70 12 of
# 42 joint weeks are reliable (39 of 42 on B65). The standard remedy is to integrate the
# latent out of both the ratio and the score, so the question asked of a held-out day is
# the one an unsampled day poses: what do the counts look like when the day's deviation is
# not known? This is integrated importance sampling: Li, Qiu, Zhang and Feng 2016, Statistics
# and Computing 26(4):881-897, who integrate out the latent variables tied to the held-out
# unit before the importance-sampling approximation and find it substantially better than the
# non-integrated form. See also Merkle, Furr and Rabe-Hesketh 2019, Psychometrika
# 84(3):802-829 (conditional against marginal likelihoods), and Vehtari, Mononen, Tolvanen,
# Sivula and Winther 2016, J. Mach. Learn. Res. 17(103):1-38 (leave-one-out approximations
# for latent variable models). The reliability bar is Vehtari, Simpson, Gelman, Yao and
# Gabry 2024, J. Mach. Learn. Res. 25(72):1-58.
#
# WHAT IT COMPUTES. For a sampled day d with trailer counts y_d1 .. y_dm,
#     p(y_d | theta) = integral of  prod_j NB2(y_dj | lambda_d / R_G_boat
#                                   * exp(sigma_uE z - sigma_uE^2 / 2), r_E)  N(z | 0, 1) dz,
# the day's JOINT density of its counts mixed over the lognormal day deviation, exactly the
# term crab_bss_pooled.stan carries (its lines for T_I), by the trapezoid rule on a fixed
# grid in z (step 0.025 over -8 to 8; the step-0.05 grid is every other node and is the
# doubling check). One column per sampled DAY replaces the conditional matrix's one column
# per count; the week's leave-out set, the PSIS smoothing, the reliability bar (Pareto k at
# or below 0.7) and the paired comparison are bss_block_cv.R's own, unchanged. The OSP
# likelihood carries no day deviation (crab_bss_pooled.stan, the OSP_I lines) and enters as
# the fit's own log_lik_osp. A fit whose term is off (A40's floor) has nothing to integrate
# and returns the conditional table.
#
# R3a's PRECONDITIONS, checked per arm and written to blockint_checks.csv before any pair is
# opened:
#   (i)  the conditional log-likelihood rebuilt here from the cached fit's draws (lambda_E_S,
#        R_G_boat, r_E, sigma_uE, uE_raw) reproduces the fit's own log_lik_trailer to 1e-6;
#   (ii) doubling the grid moves no week's leave-out elpd by more than 0.01 nats.
# Two further checks tie the arm to its run folder: the fit's pointwise lpd equals the
# folder's loo_pointwise_trailer file (so the cached fit IS that arm's fit), and the
# conditional block table rebuilt here equals the folder's loo_block table.
#
# WHAT IT NEEDS, AND WHY IT CANNOT BE RE-RUN FROM A CLONE. The arms' run folders
# (05_output/<date>/pooled-CPUE-AR-<stage>, about 15 MB each, not committed) for the Stan
# data, and the fit cache (BSS_FIT_CACHE_DIR, about 550 MB per boat fit, never committed)
# for sigma_uE and uE_raw, which ppc_draws_*.rds does not save (D54). What it writes IS
# committed, in the arm's evidence folder: loo_blockint_<stream>_<fit>.csv and
# blockint_checks.csv; `--read` then rebuilds the comparison tables from those alone.
#
# WHAT IT WRITES
#   05_output/advisory_replication_2026-11-19_arms/<stage>/loo_blockint_{trailer,osp,joint}_<fit>.csv
#   05_output/advisory_replication_2026-11-19_arms/<stage>/blockint_checks.csv
#   05_output/advisory_replication_2026-11-19_arms/<stage>/trailer_calibration.csv
#        the model's expected trailer count against the observed one, summed over the sampled
#        counts by day class (advisory or not, winter months or not). The expectation is the
#        MARGINAL mean lambda_d / R_G_boat, the quantity the reported boat effort is built
#        from (crab_bss_pooled.stan carries no day deviation into a boat day's effort). It is
#        reported, never judged (the rule's R7): it shows what an arm credits to a day class.
#   05_output/advisory_replication_2026-11-19_arms/<stage>/boat_reporting_whatif.csv
#        the arm's boat all-gear effort and expected catch as reported (a boat day keeps the
#        LEVEL: A39 treats the trailer counts' shared deviation as a snapshot nuisance, on
#        the premise that OSP's every-boat total pins the day), and what the SAME fit reports
#        if a sampled day's fitted deviation is carried into its effort and catch as it is
#        for the shore: on every sampled trailer day, and on the sampled trailer days with no
#        OSP count. Exact for the reporting rule, because generated quantities do not feed
#        back into the fit; it holds the day's crabbing share at its fitted value. Reported,
#        never judged.
#   05_output/advisory_replication_2026-11-19_blockint.csv        the pairs, by week set
#   05_output/advisory_replication_2026-11-19_blockint_weeks.csv  the per-week rows behind them
#   05_output/advisory_replication_2026-11-19_robustness.csv      for the comparisons the rule reads
#        (the transfer test's and the block CV's), what the paired standard error assumes and
#        what happens without it: the lag-1 autocorrelation of the weekly differences, a
#        Newey-West standard error (Bartlett weights, 4 lags), and a sign test. Reported, never
#        judged: the rule's thresholds are on the paired standard error.
#
# HOW TO RUN (from the repository root, after the arms are rendered)
#   BSS_FIT_CACHE_DIR=<dir> Rscript 06_diagnostics/desk_advisory_block_cv_2026-11-19.R
#   ... --stage=S26-const,S26-off     only those arms
#   ... --read                        no fit is opened: the pairs from the committed tables
#   ... --calibration-only            the calibration table alone (no quadrature; a minute per arm)
# About five minutes per arm on one core.
###############################################################################
.args <- commandArgs(trailingOnly = TRUE)
READ_ONLY <- "--read" %in% .args
CALIBRATION_ONLY <- "--calibration-only" %in% .args
ALL_STAGES <- c("S26-const", "S26-off", "S26-split", "S25-const", "S25-split", "S25-off")
.st <- sub("^--stage=", "", grep("^--stage=", .args, value = TRUE))
STAGES <- if (length(.st)) strsplit(.st[1], ",")[[1]] else ALL_STAGES
if (!all(STAGES %in% ALL_STAGES)) stop("unknown stage: ", paste(setdiff(STAGES, ALL_STAGES), collapse = ", "))

.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
.here <- function(...) file.path(.root, ...)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
setwd(.root)
suppressWarnings(suppressPackageStartupMessages(try({ library(dplyr); library(tidyr); library(lubridate); library(readxl); library(here); library(purrr); library(stringr); library(tibble) }, silent = TRUE)))
invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), function(f) source(f)))   # the block CV helpers and the advisory flag reader
source(.here("run_config.R"))

FIT_BOAT <- "private_boat_all_gear_Dungeness_Kept"
EVIDENCE <- .here("05_output", "advisory_replication_2026-11-19_arms")
OUT <- function(kind) .here("05_output", sprintf("advisory_replication_2026-11-19_%s.csv", kind))
SEASON <- c(S26 = "2025-26", S25 = "2024-25")
WINTER_MONTHS <- c(12L, 1L, 2L)          # run_config.R marine_hazard_winter_months, the split's own definition
GRID_H <- 0.025; GRID_ZMAX <- 8; K_MAX <- 0.7
TOL_LOGLIK <- 1e-6; TOL_GRID <- 0.01
rd <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL

# ---------------------------------------------------------------------------
# The integral. mu0 is [draws x n] (the count's expectation WITHOUT the day factor), y the
# counts, eday the day-effect index of each count, sig and size the draws of sigma_uE and
# r_E. Returns one column per sampled day, on the fine grid and on every other node of it.
# ---------------------------------------------------------------------------
.row_lse <- function(A) { m <- A[cbind(seq_len(nrow(A)), max.col(A, ties.method = "first"))]; m + log(rowSums(exp(A - m))) }
integrate_day_loglik <- function(mu0, y, eday, sig, size, h = GRID_H, zmax = GRID_ZMAX) {
  z <- seq(-zmax, zmax, by = h); K <- length(z); nd <- nrow(mu0)
  lw <- stats::dnorm(z, log = TRUE) + log(h)
  coarse <- seq(1L, K, by = 2L)
  shift <- outer(sig, z) - 0.5 * sig^2                      # [draws x K]: the log of the day factor at each node
  days <- sort(unique(eday))
  fine <- crs <- matrix(NA_real_, nd, length(days))
  for (j in seq_along(days)) {
    A <- matrix(0, nd, K)
    for (i in which(eday == days[j]))                       # size recycles down the rows: element [s, k] takes size[s]
      A <- A + stats::dnbinom(y[i], mu = exp(log(mu0[, i]) + shift), size = size, log = TRUE)
    A <- A + rep(lw, each = nd)
    fine[, j] <- .row_lse(A)
    crs[, j]  <- .row_lse(A[, coarse, drop = FALSE]) + log(2)
  }
  list(fine = fine, coarse = crs, days = days, n_nodes = K)
}

# ---------------------------------------------------------------------------
# One arm: find its fit in the cache by content, run the checks, write the tables.
# ---------------------------------------------------------------------------
run_dir <- function(sid) {
  ds <- Sys.glob(.here("05_output", "*", paste0("pooled-CPUE-AR-", sid, "*")))
  ds <- ds[file.exists(file.path(ds, "ADVISORY_STAGE.txt")) & file.exists(file.path(ds, sprintf("ppc_draws_%s.rds", FIT_BOAT)))]
  if (length(ds)) ds[length(ds)] else NA_character_
}
one_arm <- function(sid) {
  cat(sprintf("\n== %s ==\n", sid))
  dir <- run_dir(sid); ev <- file.path(EVIDENCE, sid)
  if (is.na(dir)) { cat("  no run folder with saved draws: skipped\n"); return(invisible(NULL)) }
  cache <- Sys.getenv("BSS_FIT_CACHE_DIR", unset = "")
  if (!nzchar(cache) || !dir.exists(cache)) stop("BSS_FIT_CACHE_DIR is not set or is not a directory: the fits' sigma_uE and uE_raw are only there")
  suppressPackageStartupMessages(library(rstan))
  sd_ <- readRDS(file.path(dir, sprintf("ppc_draws_%s.rds", FIT_BOAT)))$stan_data
  pw <- rd(file.path(dir, sprintf("loo_pointwise_trailer_%s.csv", FIT_BOAT)))
  po <- rd(file.path(dir, sprintf("loo_pointwise_osp_%s.csv", FIT_BOAT)))
  Tn <- as.integer(sd_[["T_n"]]); stopifnot(Tn > 0, nrow(pw) == Tn)
  # the arm's fit, by content: the candidate whose pointwise lpd IS the run folder's file
  cand <- Sys.glob(file.path(cache, sprintf("%s_att*_*.rds", FIT_BOAT)))
  cand <- cand[order(abs(as.numeric(file.mtime(cand)) - as.numeric(file.mtime(file.path(dir, sprintf("loo_pointwise_trailer_%s.csv", FIT_BOAT))))))]
  ex <- NULL; hit <- NA_character_; lpd_gap <- NA_real_
  for (f in cand) {
    fit <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(fit)) next
    ll <- tryCatch(as.matrix(rstan::extract(fit, pars = "log_lik_trailer")[[1]]), error = function(e) NULL)
    if (!is.null(ll) && ncol(ll) == Tn) {
      gap <- max(abs(apply(ll, 2, .bcv_lmeanexp) - pw$lpd))
      if (gap < 5e-4) {                                      # the file rounds lpd to four decimals
        pars <- c("lambda_E_S", "r_E", "R_G_boat", "log_lik_trailer", "E", "C_expected", if ((sd_[["OSP_n"]] %||% 0) > 0) "log_lik_osp",
                  if (identical(as.integer(sd_[["effort_day_effect"]] %||% 0L), 1L)) c("sigma_uE", "uE_raw"))
        ex <- rstan::extract(fit, pars = pars); hit <- f; lpd_gap <- gap
      }
    }
    rm(fit, ll); gc(verbose = FALSE)
    if (!is.null(ex)) break
  }
  if (is.null(ex)) { cat("  no cached fit reproduces this arm's pointwise file: skipped\n"); return(invisible(NULL)) }
  cat(sprintf("  fit: %s (pointwise lpd within %.1e of the run folder's file)\n", basename(hit), lpd_gap))
  nd <- nrow(ex$log_lik_trailer); G <- dim(ex$lambda_E_S)[4]
  y <- as.integer(sd_[["T_I"]])[seq_len(Tn)]; day <- as.integer(sd_[["day_T"]])[seq_len(Tn)]
  sec <- as.integer(sd_[["section_T"]] %||% rep(1L, Tn))[seq_len(Tn)]
  stopifnot(identical(y, as.integer(pw$observed)), identical(day, as.integer(pw$day_index)))
  mu0 <- vapply(seq_len(Tn), function(i) ex$lambda_E_S[, sec[i], day[i], G] / as.numeric(ex$R_G_boat), numeric(nd))
  size <- as.numeric(ex$r_E)
  has_ue <- !is.null(ex$sigma_uE)
  dir.create(ev, recursive = TRUE, showWarnings = FALSE)
  # the expected count against the observed one, by day class (reported, never judged)
  od <- as.Date(pw$event_date)
  fl <- marine_hazard_flag_series(run_config, date_start = min(od), date_end = max(od))
  adv <- fl$nws_sca_any[match(od, as.Date(fl$event_date))] == 1
  win <- as.integer(format(od, "%m")) %in% WINTER_MONTHS
  cls <- list(`all sampled counts` = rep(TRUE, Tn), `advisory days` = adv, `other days` = !adv, `advisory days, December to February` = adv & win,
              `advisory days, rest of the season` = adv & !win, `other days, December to February` = !adv & win, `other days, rest of the season` = !adv & !win)
  cal <- do.call(rbind, lapply(names(cls), function(nm) { k <- which(cls[[nm]]); if (!length(k)) return(NULL)
    e <- rowSums(mu0[, k, drop = FALSE]); o <- sum(y[k])
    data.frame(stage = sid, day_class = nm, counts = length(k), days = length(unique(od[k])), observed_trailers = o, expected_trailers_mean = mean(e),
               expected_trailers_median = stats::median(e), expected_lo95 = unname(stats::quantile(e, 0.025)), expected_hi95 = unname(stats::quantile(e, 0.975)),
               expected_over_observed = mean(e) / o, median_over_observed = stats::median(e) / o,
               ratio_lo95 = unname(stats::quantile(e, 0.025)) / o, ratio_hi95 = unname(stats::quantile(e, 0.975)) / o, stringsAsFactors = FALSE) }))
  utils::write.csv(cal, file.path(ev, "trailer_calibration.csv"), row.names = FALSE)
  cat(sprintf("  expected over observed trailers: all %.3f; advisory days %.3f; other days %.3f\n", cal$expected_over_observed[1], cal$expected_over_observed[2], cal$expected_over_observed[3]))
  # what the same fit reports if a sampled day's own deviation is carried (reported, never judged)
  if (has_ue && !is.null(ex$E) && !is.null(ex$C_expected)) {
    sg <- as.numeric(ex$sigma_uE); ed <- as.integer(sd_[["Eday_T"]])[seq_len(Tn)]; uu <- as.matrix(ex$uE_raw)
    Ed <- apply(ex$E, c(1, 3), sum); Cd <- apply(ex$C_expected, c(1, 3), sum)                 # [draws x days], summed over sections and gear groups
    sday <- sort(unique(day)); dev <- matrix(1, nd, ncol(Ed))
    for (dd in sday) dev[, dd] <- exp(sg * uu[, ed[match(dd, day)]] - 0.5 * sg^2)
    osp_days <- if ((sd_[["OSP_n"]] %||% 0) > 0) unique(as.integer(sd_[["day_OSP"]])[seq_len(sd_[["OSP_n"]])]) else integer(0)
    dev_no_osp <- dev; if (length(osp_days)) dev_no_osp[, osp_days] <- 1
    qs <- function(x) c(median = stats::median(x), lo95 = unname(stats::quantile(x, 0.025)), hi95 = unname(stats::quantile(x, 0.975)))
    wf <- function(label, dv, ndays) { e <- qs(rowSums(Ed * dv)); cc <- qs(rowSums(Cd * dv))
      data.frame(stage = sid, reporting = label, days_carrying_a_deviation = ndays, effort_median = e[["median"]], effort_lo95 = e[["lo95"]], effort_hi95 = e[["hi95"]],
                 expected_catch_median = cc[["median"]], expected_catch_lo95 = cc[["lo95"]], expected_catch_hi95 = cc[["hi95"]], stringsAsFactors = FALSE) }
    wi <- rbind(wf("as reported: every boat day keeps the level", matrix(1, nd, ncol(Ed)), 0L),
                wf("the deviation carried on every sampled trailer day", dev, length(sday)),
                wf("the deviation carried on sampled trailer days with no OSP count", dev_no_osp, length(setdiff(sday, osp_days))))
    wi$effort_over_reported <- wi$effort_median / wi$effort_median[1]; wi$catch_over_reported <- wi$expected_catch_median / wi$expected_catch_median[1]
    utils::write.csv(wi, file.path(ev, "boat_reporting_whatif.csv"), row.names = FALSE)
    cat(sprintf("  reporting what-if, expected catch: as reported %.0f; deviation carried on every sampled day %.0f (%.3f); on sampled days with no OSP count %.0f (%.3f)\n",
                wi$expected_catch_median[1], wi$expected_catch_median[2], wi$catch_over_reported[2], wi$expected_catch_median[3], wi$catch_over_reported[3]))
  }
  if (CALIBRATION_ONLY) return(invisible(cal))
  wk_obs <- bss_block_weeks(pw$event_date)
  streams_c <- list(trailer = list(ll = ex$log_lik_trailer, block = wk_obs))
  if (has_ue) {
    sig <- as.numeric(ex$sigma_uE); eday <- as.integer(sd_[["Eday_T"]])[seq_len(Tn)]
    ue <- as.matrix(ex$uE_raw)
    cond <- vapply(seq_len(Tn), function(i) stats::dnbinom(y[i], mu = mu0[, i] * exp(sig * ue[, eday[i]] - 0.5 * sig^2), size = size, log = TRUE), numeric(nd))
    gap_cond <- max(abs(cond - ex$log_lik_trailer))
    t0 <- Sys.time(); int <- integrate_day_loglik(mu0, y, eday, sig, size)
    cat(sprintf("  integrated %d counts on %d days over %d nodes in %.1f min\n", Tn, length(int$days), int$n_nodes, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    day_date <- as.Date(pw$event_date)[match(int$days, eday)]
    stopifnot(all(tapply(pw$event_date, eday, function(v) length(unique(v))) == 1))     # one date per day-effect index
    wk_day <- bss_block_weeks(day_date)
    streams_f <- list(trailer = list(ll = int$fine, block = wk_day)); streams_h <- list(trailer = list(ll = int$coarse, block = wk_day))
  } else {
    gap_cond <- max(abs(vapply(seq_len(Tn), function(i) stats::dnbinom(y[i], mu = mu0[, i], size = size, log = TRUE), numeric(nd)) - ex$log_lik_trailer))
    sig <- NA_real_; streams_f <- streams_c; streams_h <- streams_c
    cat("  the fit carries no day deviation: the integrated table is the conditional one\n")
  }
  if (!is.null(ex$log_lik_osp)) { stopifnot(!is.null(po), nrow(po) == ncol(ex$log_lik_osp))
    o <- list(ll = as.matrix(ex$log_lik_osp), block = bss_block_weeks(po$event_date))
    streams_c$osp <- o; streams_f$osp <- o; streams_h$osp <- o }
  tc <- bss_block_cv_fit(streams_c, k_max = K_MAX); tf <- bss_block_cv_fit(streams_f, k_max = K_MAX); th <- bss_block_cv_fit(streams_h, k_max = K_MAX)
  top <- if ("joint" %in% names(tf)) "joint" else "trailer"
  # the conditional table against the driver's own
  drv <- rd(file.path(dir, sprintf("loo_block_%s_%s.csv", top, FIT_BOAT)))
  gap_drv <- if (is.null(drv)) NA_real_ else { m <- merge(drv[, c("block", "elpd_block")], tc[[top]][, c("block", "elpd_block")], by = "block"); max(abs(m$elpd_block.x - m$elpd_block.y)) }
  both <- merge(tf[[top]][, c("block", "elpd_block", "reliable")], th[[top]][, c("block", "elpd_block", "reliable")], by = "block")
  gap_grid <- max(abs(both$elpd_block.x - both$elpd_block.y)[both$reliable.x & both$reliable.y])
  for (sn in names(tf)) utils::write.csv(cbind(data_type = sn, tf[[sn]], stringsAsFactors = FALSE),
                                         file.path(ev, sprintf("loo_blockint_%s_%s.csv", sn, FIT_BOAT)), row.names = FALSE)
  q <- function(x, p) if (all(is.na(x))) NA_real_ else unname(stats::quantile(x, p))
  chk <- data.frame(stage = sid, season = SEASON[[sub("-.*$", "", sid)]], arm = sub("^[^-]*-", "", sid), cached_fit = basename(hit), n_draws = nd,
    trailer_counts = Tn, trailer_days = length(unique(day)), osp_counts = as.integer(sd_[["OSP_n"]] %||% 0L), effort_day_effect = as.integer(has_ue),
    scored_row = top, pointwise_lpd_max_abs_gap = lpd_gap, conditional_loglik_max_abs_gap = gap_cond, precondition_i = isTRUE(gap_cond <= TOL_LOGLIK),
    conditional_block_table_max_abs_gap = gap_drv, grid_step = GRID_H, grid_doubling_max_abs_elpd_gap = gap_grid, precondition_ii = isTRUE(gap_grid <= TOL_GRID),
    weeks = nrow(tf[[top]]), reliable_conditional = sum(tc[[top]]$reliable), reliable_integrated = sum(tf[[top]]$reliable),
    median_k_conditional = stats::median(tc[[top]]$pareto_k), median_k_integrated = stats::median(tf[[top]]$pareto_k),
    elpd_reliable_integrated = sum(tf[[top]]$elpd_block[tf[[top]]$reliable]),
    sigma_uE_mean = mean(sig), sigma_uE_sd = stats::sd(sig), sigma_uE_q2.5 = q(sig, 0.025), sigma_uE_q50 = q(sig, 0.5), sigma_uE_q97.5 = q(sig, 0.975),
    r_E_mean = mean(size), stringsAsFactors = FALSE)
  utils::write.csv(chk, file.path(ev, "blockint_checks.csv"), row.names = FALSE)
  cat(sprintf("  (i) conditional log-likelihood gap %.2e (%s); (ii) grid doubling gap %.2e nats (%s); driver's table gap %.2e\n",
              gap_cond, if (chk$precondition_i) "met" else "NOT MET", gap_grid, if (chk$precondition_ii) "met" else "NOT MET", gap_drv))
  cat(sprintf("  %s row: %d weeks; reliable %d conditional, %d integrated (median k %.2f against %.2f); sigma_uE %.3f (SD %.3f)\n",
              top, chk$weeks, chk$reliable_conditional, chk$reliable_integrated, chk$median_k_conditional, chk$median_k_integrated, chk$sigma_uE_mean, chk$sigma_uE_sd))
  invisible(chk)
}

# ---------------------------------------------------------------------------
# The pairs, from the committed tables alone.
# ---------------------------------------------------------------------------
tab <- function(sid, stream, kind = "loo_blockint") rd(file.path(EVIDENCE, sid, sprintf("%s_%s_%s.csv", kind, stream, FIT_BOAT)))
# A WINTER WEEK, as the transfer script defines it: at least half of the week's sampled days fall in the
# split's winter months. Read from the dates of the season's const arm (every arm of a season holds the
# same observations).
winter_weeks <- function(sn) {
  d <- unique(as.Date(unlist(lapply(c("trailer", "osp"), function(st) {
    x <- rd(file.path(EVIDENCE, paste0(sn, "-const"), sprintf("loo_pointwise_%s_%s.csv", st, FIT_BOAT))); if (is.null(x)) NULL else x$event_date }))))
  if (!length(d)) return(NULL)
  sh <- tapply(as.integer(format(d, "%m")) %in% WINTER_MONTHS, bss_block_weeks(d), mean)
  names(sh)[sh >= 0.5]
}
pairs_all <- function() {
  P <- list(); W <- list()
  for (sn in names(SEASON)) for (pr in list(c("const", "off"), c("split", "const"), c("split", "off"))) for (kind in c("loo_blockint", "loo_block")) {
    ww <- winter_weeks(sn)
    for (stream in c("joint", "trailer", "osp")) {
      ta <- tab(paste0(sn, "-", pr[1]), stream, kind); tb <- tab(paste0(sn, "-", pr[2]), stream, kind)
      if (is.null(ta) || is.null(tb)) next
      w <- bss_block_cv_weeks(tb, ta, k_max = K_MAX)         # the helper returns its SECOND table minus its first: a minus b
      w$winter <- w$block %in% ww
      id <- data.frame(season = SEASON[[sn]], pair = sprintf("%s - %s", pr[1], pr[2]),
                       table = if (kind == "loo_blockint") "day deviation integrated out (R3a)" else "conditional on the day deviation (R3)",
                       stream = stream, stringsAsFactors = FALSE)
      W[[length(W) + 1]] <- cbind(id, w)
      for (set in c("all weeks", "winter weeks (December to February)", "rest of the season")) {
        keep <- if (set == "all weeks") rep(TRUE, nrow(w)) else if (grepl("^winter", set)) w$winter else !w$winter
        d <- w$diff[keep & w$used]; n <- length(d); se <- if (n >= 2) sqrt(n) * stats::sd(d) else NA_real_
        P[[length(P) + 1]] <- cbind(id, data.frame(week_set = set, n_weeks = sum(keep), n_used = n, reliable_share = if (sum(keep)) n / sum(keep) else NA_real_,
          diff = if (n) sum(d) else NA_real_, se = se, ratio = if (is.finite(se) && se > 0) sum(d) / se else NA_real_, n_positive = sum(d > 0),
          evaluable = set == "all weeks" && sum(keep) > 0 && n / sum(keep) >= 0.7 && n >= 8, stringsAsFactors = FALSE))
      }
    }
  }
  list(pairs = if (length(P)) do.call(rbind, P) else NULL, weeks = if (length(W)) do.call(rbind, W) else NULL)
}

# ---------------------------------------------------------------------------
# What the paired standard error assumes (weeks independent), and the readings without it.
# ---------------------------------------------------------------------------
robust_row <- function(d, label_cols) {
  n <- length(d); if (n < 6) return(NULL)
  m <- mean(d); e <- d - m; g <- function(l) sum(e[seq_len(n - l)] * e[seq_len(n - l) + l]) / n
  L <- 4L; lrv <- g(0) + 2 * sum(vapply(seq_len(L), function(l) (1 - l / (L + 1)) * g(l), numeric(1)))
  se_iid <- sqrt(n) * stats::sd(d); se_nw <- sqrt(max(lrv, 0) * n); pos <- sum(d > 0)
  cbind(label_cols, data.frame(n_weeks = n, sum = sum(d), se_paired = se_iid, ratio_paired = sum(d) / se_iid, lag1_autocorrelation = g(1) / g(0),
                               se_newey_west_4 = se_nw, ratio_newey_west_4 = if (se_nw > 0) sum(d) / se_nw else NA_real_, n_positive = pos,
                               sign_test_p_one_sided = stats::pbinom(pos - 1, n, 0.5, lower.tail = FALSE), stringsAsFactors = FALSE))
}
robustness <- function(weeks) {
  out <- list()
  tw <- rd(OUT("transfer_weeks"))
  if (!is.null(tw)) for (tg in c("2025-26", "2024-25")) { src <- setdiff(c("2025-26", "2024-25"), tg)
    x <- tw[tw$target == tg & tw$source == src & tw$scored %in% TRUE, ]; x <- x[order(x$week), ]
    for (cmp in list(c("TC", "M0", "R2: carried coefficient against no term"), c("TS", "TC", "R6: carried split against carried constant")))
      out[[length(out) + 1]] <- robust_row(x[[cmp[1]]] - x[[cmp[2]]], data.frame(instrument = "transfer test", comparison = cmp[3], target_season = tg, source = src, stringsAsFactors = FALSE)) }
  if (!is.null(weeks)) for (sn in unique(weeks$season)) for (pr in c("const - off", "split - const")) {
    x <- weeks[weeks$season == sn & weeks$pair == pr & grepl("integrated", weeks$table), ]; if (!nrow(x)) next
    x <- x[x$stream == (if ("joint" %in% x$stream) "joint" else "trailer") & x$used %in% TRUE, ]; x <- x[order(x$block), ]
    out[[length(out) + 1]] <- robust_row(x$diff, data.frame(instrument = "fitted model, held-out weeks (R3a table)", comparison = pr, target_season = sn, source = "its own fit", stringsAsFactors = FALSE)) }
  if (length(out)) do.call(rbind, out) else NULL
}

if (!READ_ONLY) for (sid in STAGES) tryCatch(one_arm(sid), error = function(e) cat(sprintf("  *** %s FAILED: %s\n", sid, conditionMessage(e))))
res <- pairs_all()
if (is.null(res$pairs)) { cat("\nNo pair of arms has its tables yet.\n")
} else {
  utils::write.csv(res$pairs, OUT("blockint"), row.names = FALSE); utils::write.csv(res$weeks, OUT("blockint_weeks"), row.names = FALSE)
  cat("\n== The pairs (a minus b, held-out weeks reliable in both arms) ==\n")
  s <- res$pairs[res$pairs$week_set == "all weeks", ]
  for (i in seq_len(nrow(s))) cat(sprintf("  %s  %-13s %-38s %-7s %+7.1f nats, SE %5.1f, ratio %+5.2f; %d of %d weeks, %d positive%s\n", s$season[i], s$pair[i], s$table[i], s$stream[i],
      s$diff[i], s$se[i], s$ratio[i], s$n_used[i], s$n_weeks[i], s$n_positive[i], if (isTRUE(s$evaluable[i])) "" else "  [not evaluable]"))
  cat(sprintf("\nWrote %s and %s\n", basename(OUT("blockint")), basename(OUT("blockint_weeks"))))
}
rb <- robustness(res$weeks)
if (!is.null(rb)) { utils::write.csv(rb, OUT("robustness"), row.names = FALSE)
  cat("\n== What the paired SE assumes, and the readings without it ==\n")
  for (i in seq_len(nrow(rb))) cat(sprintf("  %-40s %-44s %s: %+6.1f, paired %.2f SE, Newey-West %.2f SE, lag-1 r %+.2f, %d of %d positive (sign test p %.4f)\n", rb$instrument[i], rb$comparison[i],
      rb$target_season[i], rb$sum[i], rb$ratio_paired[i], rb$ratio_newey_west_4[i], rb$lag1_autocorrelation[i], rb$n_positive[i], rb$n_weeks[i], rb$sign_test_p_one_sided[i])) }
