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
# desk_boat_day_deviation_osp_2026-11-25.R
#
# DESK, MINUTES, NO MCMC (CHANGE_REGISTER B86, for D55).
#
# THE QUESTION. Since A39 a sampled trailer day's counts share one fitted deviation, and the
# reported boat effort ignores it: a boat day is credited with the level. A39's reason was
# that "OSP's every-boat total pins it; the trailer counts' shared deviation is a snapshot
# nuisance". In the model a day's boat effort is lambda_E[d] * L[d], and an OSP count observes
# exactly that product, so on an OSP day the premise holds: OSP's count sets L[d]. On a sampled
# trailer day WITHOUT an OSP count (every December to February day, 92 of the 143 sampled
# trailer days of 2024-25, every boat day of 2025-26 so far) L[d] is drawn from its prior around
# tau_bar, and the trailer deviation is dropped. Whether that is right depends on what the
# deviation is: more boats out that day (real effort, which the day's total should carry), or
# more trailers in the lot per boat (longer trips, the snapshot's timing: not more boats). The
# two cannot be told apart on a trailer-only day. They can on the days that have BOTH counts,
# because there OSP measures the boats.
#
# THE TEST. On the 2024-25 boat all-gear fit as shipped (the B80 arm S25-const, which reproduces
# the authoritative run B70 byte for byte), for every day with both a trailer count and an OSP
# count:
#   delta_d  the trailer day's fitted log deviation, the posterior mean of
#            sigma_uE * uE_raw[day] - sigma_uE^2 / 2. It is identified by the trailer counts alone:
#            the OSP likelihood carries no day deviation.
#   r_d      log(OSP count) - the posterior mean of log(lambda_E[d] * tau_bar / R_G_boat): the
#            day's actual boat total against what a day WITHOUT an OSP count is credited (the
#            level at the central turnover).
# Regress r_d on delta_d (ordinary least squares with an intercept; heteroscedasticity-robust
# HC3 standard error). Using posterior means of delta_d is regression calibration, so the
# slope is not attenuated by delta_d's own uncertainty.
# WHY 0.5 IS THE LINE. If r = a + b delta + e, then crediting a trailer-only day with the level
# misses by b delta + e, and carrying the deviation misses by (b - 1) delta + e. The second is
# the smaller in expected square exactly when b > 0.5. So the slope says which reporting rule
# predicts a day's actual boats better, and by how much.
#
# THE RULE, WRITTEN AND COMMITTED BEFORE r_d OR delta_d WAS COMPUTED ON ANY DAY.
#   Read on the slope b and its 95% interval (HC3):
#     lower bound above 0.5   the trailer deviation is mostly real boat effort; A39's premise
#                             fails on trailer-only days, and carrying the deviation there (as
#                             the shore does) predicts a day's boats better than the level.
#     upper bound below 0.5   the deviation is mostly snapshot; keeping the level is right.
#     covers 0.5              cannot tell: no change is supported by this test.
#   Fewer than 20 overlap days: too few to read.
#   Reported beside it, not read: the paired comparison of squared errors per day, level against
#   carried (mean difference and its SE); the slope on days with two or more trailer counts only;
#   a day bootstrap interval of b (2,000 resamples, seed 1); and a model-light version with
#   delta_d replaced by log(the day's mean trailer count) - log(the posterior mean of
#   lambda_E[d] / R_G_boat).
#   STATED IN ADVANCE: (1) the overlap days run March to September; the test says nothing direct
#   about December to February, where no OSP day exists. (2) OSP's count is every returning
#   private vessel, the trailer lot holds trailers: the same fleet, measured two ways. (3) Nothing
#   changes on this result; D55 is Matt's decision, and this is evidence for it.
#   EXPECTED: b positive; on which side of 0.5 I do not know.
#
# WRITES 05_output/boat_day_deviation_osp_2026-11-25.csv (the readings) and
#        05_output/boat_day_deviation_osp_2026-11-25_days.csv (the overlap days)
#
# HOW TO RUN (from the repository root; needs the S25-const run folder and the fit cache)
#   BSS_FIT_CACHE_DIR=<dir> RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_boat_day_deviation_osp_2026-11-25.R
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
suppressWarnings(suppressPackageStartupMessages(try({ library(dplyr); library(tidyr); library(lubridate); library(readxl); library(here); library(purrr); library(stringr); library(tibble) }, silent = TRUE)))
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
suppressPackageStartupMessages(library(rstan))

FIT <- "private_boat_all_gear_Dungeness_Kept"
DIR <- file.path("05_output", "20261119", "pooled-CPUE-AR-S25-const")
OUT <- "05_output/boat_day_deviation_osp_2026-11-25.csv"; OUT_DAYS <- "05_output/boat_day_deviation_osp_2026-11-25_days.csv"
MIN_DAYS <- 20L
CACHE <- Sys.getenv("BSS_FIT_CACHE_DIR", unset = "")
if (!nzchar(CACHE) || !dir.exists(CACHE)) stop("BSS_FIT_CACHE_DIR is not set or is not a directory")
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE)

sd_ <- readRDS(file.path(DIR, sprintf("ppc_draws_%s.rds", FIT)))$stan_data
pwT <- rd(file.path(DIR, sprintf("loo_pointwise_trailer_%s.csv", FIT))); pwO <- rd(file.path(DIR, sprintf("loo_pointwise_osp_%s.csv", FIT)))
ex <- NULL
for (f in Sys.glob(file.path(CACHE, sprintf("%s_att*_*.rds", FIT)))) {
  x <- readRDS(f); ll <- as.matrix(rstan::extract(x, pars = "log_lik_trailer")[[1]])
  if (ncol(ll) == nrow(pwT) && max(abs(apply(ll, 2, .bcv_lmeanexp) - pwT$lpd)) < 5e-4) {
    ex <- rstan::extract(x, pars = c("lambda_E_S", "R_G_boat", "tau_bar", "sigma_uE", "uE_raw")); cat("fit:", basename(f), "\n"); break }
  rm(x); gc(verbose = FALSE)
}
if (is.null(ex)) stop("no cached fit reproduces the S25-const boat all-gear fit")
stopifnot(identical(as.integer(sd_$effort_day_effect), 1L), identical(as.integer(sd_$shared_tau), 1L), identical(as.integer(sd_$osp_scale_is_tau), 1L))

nd <- length(ex$R_G_boat); G <- dim(ex$lambda_E_S)[4]
Tn <- as.integer(sd_$T_n); On <- as.integer(sd_$OSP_n)
tday <- as.integer(sd_$day_T)[seq_len(Tn)]; tsec <- as.integer(sd_$section_T %||% rep(1L, Tn))[seq_len(Tn)]; ted <- as.integer(sd_$Eday_T)[seq_len(Tn)]; ty <- as.integer(sd_$T_I)[seq_len(Tn)]
oday <- as.integer(sd_$day_OSP)[seq_len(On)]; osec <- as.integer(sd_$section_OSP %||% rep(1L, On))[seq_len(On)]; oy <- as.integer(sd_$OSP_I)[seq_len(On)]
both <- intersect(unique(tday), unique(oday))
cat(sprintf("trailer days %d, OSP days %d, days with both %d\n", length(unique(tday)), length(unique(oday)), length(both)))
sig <- as.numeric(ex$sigma_uE); uu <- as.matrix(ex$uE_raw); tb <- as.numeric(ex$tau_bar); rg <- as.numeric(ex$R_G_boat)

days <- do.call(rbind, lapply(both, function(d) {
  i <- which(tday == d); o <- which(oday == d)
  stopifnot(length(unique(ted[i])) == 1L, length(o) == 1L)
  delta <- mean(sig * uu[, ted[i[1]]] - 0.5 * sig^2)
  lam_t <- ex$lambda_E_S[, tsec[i[1]], d, G]; lam_o <- ex$lambda_E_S[, osec[o], d, G]
  data.frame(day_index = d, event_date = pwT$event_date[match(d, pwT$day_index)], trailer_counts = length(i), trailers_mean = mean(ty[i]),
             osp_count = oy[o], delta = delta,
             r = log(max(oy[o], 0.5)) - mean(log(lam_o * tb / rg)),
             r_modellight_x = log(max(mean(ty[i]), 0.5)) - mean(log(lam_t / rg)), stringsAsFactors = FALSE)
}))
days <- days[order(days$event_date), ]

hc3 <- function(y, x) {
  X <- cbind(1, x); XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y); e <- as.numeric(y - X %*% b)
  h <- rowSums((X %*% XtXi) * X); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  c(a = b[1], b = b[2], se = sqrt(V[2, 2]), n = length(y), r2 = 1 - sum(e^2) / sum((y - mean(y))^2))
}
reading_of <- function(lo, hi, n) if (n < MIN_DAYS) "too few overlap days to read" else
  if (lo > 0.5) "the trailer deviation is mostly real boat effort: carrying it predicts a trailer-only day's boats better than the level (A39's premise fails there)" else
  if (hi < 0.5) "the trailer deviation is mostly snapshot: keeping the level is right" else "cannot tell: the interval covers 0.5; no change is supported by this test"
res <- list()
row <- function(what, fit, read = FALSE, extra = "") {
  lo <- fit[["b"]] - 1.96 * fit[["se"]]; hi <- fit[["b"]] + 1.96 * fit[["se"]]
  res[[length(res) + 1L]] <<- data.frame(what = what, days = fit[["n"]], intercept = fit[["a"]], slope = fit[["b"]], se_hc3 = fit[["se"]], lo95 = lo, hi95 = hi,
                                         r_squared = fit[["r2"]], reading = if (read) reading_of(lo, hi, fit[["n"]]) else "", note = extra, stringsAsFactors = FALSE)
}
main <- hc3(days$r, days$delta); row("PRIMARY: r_d on delta_d, every overlap day", main, read = TRUE)
k2 <- days$trailer_counts >= 2; row("days with two or more trailer counts", hc3(days$r[k2], days$delta[k2]))
row("model-light: r_d on log(mean trailers) - log(lambda / R_G_boat)", hc3(days$r, days$r_modellight_x))
set.seed(1); bb <- replicate(2000, { j <- sample.int(nrow(days), replace = TRUE); hc3(days$r[j], days$delta[j])[["b"]] })
res[[length(res) + 1L]] <- data.frame(what = "day bootstrap of the primary slope (2,000, seed 1)", days = nrow(days), intercept = NA, slope = stats::median(bb), se_hc3 = stats::sd(bb),
                                      lo95 = unname(stats::quantile(bb, 0.025)), hi95 = unname(stats::quantile(bb, 0.975)), r_squared = NA, reading = "", note = "percentile interval", stringsAsFactors = FALSE)
# the two reporting rules as predictors of the day's boats, after the common offset (intercept) is removed
r0 <- days$r - mean(days$r); dl <- days$delta - mean(days$delta)
dsq <- (r0 - dl)^2 - r0^2                                   # carried minus level, per day
res[[length(res) + 1L]] <- data.frame(what = "squared error per day, carried minus level (negative favours carrying)", days = nrow(days), intercept = NA, slope = mean(dsq), se_hc3 = stats::sd(dsq) / sqrt(nrow(days)),
                                      lo95 = mean(dsq) - 1.96 * stats::sd(dsq) / sqrt(nrow(days)), hi95 = mean(dsq) + 1.96 * stats::sd(dsq) / sqrt(nrow(days)), r_squared = NA, reading = "",
                                      note = sprintf("mean squared error: level %.4f, carried %.4f", mean(r0^2), mean((r0 - dl)^2)), stringsAsFactors = FALSE)
out <- do.call(rbind, res)
op <- options(width = 250); print(out[, c("what", "days", "slope", "se_hc3", "lo95", "hi95", "r_squared")], row.names = FALSE, digits = 3); options(op)
cat("\nB86 (on the primary slope):", out$reading[1], "\n")
cat(sprintf("delta_d: mean %+.3f, SD %.3f across the overlap days; r_d: mean %+.3f, SD %.3f\n", mean(days$delta), stats::sd(days$delta), mean(days$r), stats::sd(days$r)))
utils::write.csv(out, OUT, row.names = FALSE); utils::write.csv(days, OUT_DAYS, row.names = FALSE)
cat("Wrote", OUT, "and", OUT_DAYS, "\n")
