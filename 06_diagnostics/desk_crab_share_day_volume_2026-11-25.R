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
# desk_crab_share_day_volume_2026-11-25.R
#
# DESK, MINUTES, NO MCMC (CHANGE_REGISTER B87; opens D56).
#
# THE QUESTION. A34 (for D40) makes the crabbing share fall with the day's boat volume,
# logit f_d = eta_f[month] + beta_fvol (log V_d - centre), so that "the share the totals apply is
# boat-weighted". A day's crab effort is V_d f(V_d), which is concave in log V_d (its slope in log V
# is 1 + beta_fvol (1 - f), between 0.24 and 1 on B70's beta_fvol of -0.76). The generated
# quantities evaluate it at the MODEL's day volume, lambda_E[d] L[d] / R_G_boat: the level, which
# varies across months, weekends, holidays and advisory days and, through L[d], by a fixed 0.15 on
# the log scale. Real boat days vary around that far more: the trailer-identified day deviation has
# SD sigma_uE (0.85 on B70), and B86 showed it is boat effort (slope 0.93 against OSP's counts). On
# a concave function the value at the mean volume exceeds the mean of the values, so the crab effort
# credited to a day at its level is above the day's expected crab effort: the boat-weighting A34
# brought covers the model's day-to-day variation and not the residual one, which is most of it.
# This prices that, exactly for the reporting (generated quantities do not feed back into the fit).
#
# THE COMPUTATION. On the authoritative configuration's boat all-gear fit (B80's S25-const, byte-
# identical to B70), for every draw and every day of the window, the factor applied to the day's
# effort and catch is changed from f(v_d) to the expectation over the day's volume,
#     F_d = integral of exp(delta) f(v_d exp(delta)) N(z | 0, 1) dz,   delta = s z - s^2 / 2,
# a mean-one lognormal day factor (so the expected BOATS are unchanged and only the share is
# weighted), on the trapezoid rule over z in [-8, 8] at step 0.05; v_d = lambda_E[d] L[d] / R_G_boat
# as the generated quantities form it, and f as they form it (1e-6 + (1 - 2e-6) inv_logit(...)).
# s is the residual day SD of log boat volume beyond what L[d] already carries: s^2 = s_day^2 -
# shared_tau_sigma^2, with s_day = sigma_uE draw by draw (PRIMARY), and 0.5 and 0.82 (the SD of the
# OSP days' log departure from the level in B86) as the sensitivity.
#
# THE RULE, WRITTEN AND COMMITTED BEFORE ANY F_d WAS COMPUTED. Reported, never judged (a pricing):
#   the boat all-gear expected catch (median of C_expected_sum) as reported, and with F_d on every
#   day (PRIMARY, s_day = sigma_uE), with F_d at s_day 0.5 and 0.82, and, beside them, with B86's
#   switch on in its volume-consistent form (a sampled trailer day without an OSP count carries its
#   own realised deviation, every other day takes F_d). The port effect is the boat change added to
#   B70's port median (a component change, as R7 of B80 reported).
#   Reading: if the PRIMARY change is larger than 1% of the boat all-gear catch, D56 is opened as an
#   open item that moves the number, a decision for Matt; below 1% it is opened as low.
#   STATED IN ADVANCE. (1) It assumes the day's residual log volume is normal with the fitted SD,
#   independent of the share's own residual; on OSP days the realised count is known and a model
#   could use it, which this does not. (2) It changes nothing; the fix is a model change (the share
#   integrated over the day volume in the generated quantities, or a day volume deviation shared by
#   the OSP and trailer streams), to be decided by Matt and rendered under a rule.
#   EXPECTED: a reduction of a few percent of the boat all-gear catch, larger in summer (low share,
#   steep in volume) than in winter (share near 0.9).
#
# WRITES 05_output/crab_share_day_volume_2026-11-25.csv, and _months.csv by month.
#
# HOW TO RUN (needs the S25-const run folder and the fit cache)
#   BSS_FIT_CACHE_DIR=<dir> RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_crab_share_day_volume_2026-11-25.R
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
setwd(.root)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
suppressWarnings(suppressPackageStartupMessages(try({ library(dplyr); library(tidyr); library(lubridate); library(readxl); library(here); library(purrr); library(stringr); library(tibble) }, silent = TRUE)))
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
suppressPackageStartupMessages(library(rstan))
FIT <- "private_boat_all_gear_Dungeness_Kept"
DIR <- file.path("05_output", "20261119", "pooled-CPUE-AR-S25-const")
PORT_B70 <- 86497
OUT <- "05_output/crab_share_day_volume_2026-11-25.csv"; OUT_M <- "05_output/crab_share_day_volume_2026-11-25_months.csv"
CACHE <- Sys.getenv("BSS_FIT_CACHE_DIR", unset = "")
if (!nzchar(CACHE) || !dir.exists(CACHE)) stop("BSS_FIT_CACHE_DIR is not set or is not a directory")

sd_ <- readRDS(file.path(DIR, sprintf("ppc_draws_%s.rds", FIT)))$stan_data
pw <- utils::read.csv(file.path(DIR, sprintf("loo_pointwise_trailer_%s.csv", FIT)))
fit <- NULL
for (f in Sys.glob(file.path(CACHE, sprintf("%s_att*_*.rds", FIT)))) {
  x <- readRDS(f); ll <- as.matrix(rstan::extract(x, pars = "log_lik_trailer")[[1]])
  if (ncol(ll) == nrow(pw) && max(abs(apply(ll, 2, .bcv_lmeanexp) - pw$lpd)) < 5e-4) { fit <- x; cat("fit:", basename(f), "\n"); break }
  rm(x); gc(verbose = FALSE)
}
if (is.null(fit)) stop("no cached fit reproduces S25-const")
mat <- function(p) as.matrix(fit, pars = p)
D <- as.integer(sd_$D); Gd <- as.integer(sd_$G %||% 1L)
lam <- mat("lambda_E_S"); li <- do.call(rbind, lapply(strsplit(sub("^lambda_E_S\\[(.*)\\]$", "\\1", colnames(lam)), ","), as.integer))
Lm <- mat("L_out"); rg <- mat("R_G_boat")[, 1]; etaf <- mat("eta_f"); bf <- mat("beta_fvol")[, 1]
sg <- mat("sigma_uE")[, 1]; uu <- mat("uE_raw")
nd <- length(rg)
byday <- function(p) { x <- mat(p); d <- as.integer(sub(sprintf("^%s\\[[0-9]+,([0-9]+),[0-9]+\\]$", p), "\\1", colnames(x))); sapply(seq_len(D), function(k) rowSums(x[, d == k, drop = FALSE])) }
Ed <- byday("E"); Cd <- byday("C_expected")
fstr <- as.integer(sd_$f_stratum); cen <- as.numeric(sd_$fvol_centre); tau_s <- as.numeric(sd_$shared_tau_sigma)
share <- function(v, k) 1e-6 + (1 - 2e-6) * stats::plogis(etaf[, k] + bf * (log(v) - cen[k]))
vday <- sapply(seq_len(D), function(d) pmax(rowSums(lam[, li[, 2] == d & li[, 3] == Gd, drop = FALSE]) / rg * Lm[, d], 1e-9))
f0 <- sapply(seq_len(D), function(d) share(vday[, d], fstr[d]))
z <- seq(-8, 8, by = 0.05); wz <- stats::dnorm(z) * 0.05
integ_ratio <- function(s) {          # [draws x days]: F_d / f(v_d) at residual SD s (a draws-length vector)
  s <- sqrt(pmax(s^2 - tau_s^2, 0))
  sapply(seq_len(D), function(d) {
    dl <- outer(s, z) - 0.5 * s^2                                 # [draws x K]
    fz <- 1e-6 + (1 - 2e-6) * stats::plogis(etaf[, fstr[d]] + bf * (log(vday[, d]) + dl - cen[fstr[d]]))
    as.numeric((exp(dl) * fz) %*% wz) / f0[, d] })
}
cat(sprintf("checks: the mean-one day factor integrates to %.8f at s = 0.85; the reported catch rebuilt from the days %.1f against C_expected_sum %.1f\n",
            sum(exp(0.85 * z - 0.85^2 / 2) * wz), stats::median(rowSums(Cd)), stats::median(mat("C_expected_sum")[, 1])))
base_C <- rowSums(Cd); base_E <- rowSums(Ed)
rows <- list(); months <- list()

dates <- seq(as.Date(min(pw$event_date)) - (pw$day_index[which.min(as.Date(pw$event_date))] - 1L), by = "day", length.out = D)
mon <- format(dates, "%Y-%m")
stopifnot(all(format(as.Date(pw$event_date), "%Y-%m") == mon[pw$day_index]))
one <- function(label, ratio) {
  Cn <- rowSums(Cd * ratio); En <- rowSums(Ed * ratio)
  rows[[length(rows) + 1L]] <<- data.frame(variant = label, boat_catch_median = stats::median(Cn), boat_catch_lo95 = unname(stats::quantile(Cn, 0.025)),
    boat_catch_hi95 = unname(stats::quantile(Cn, 0.975)), change_pct = 100 * (stats::median(Cn) / stats::median(base_C) - 1),
    effort_change_pct = 100 * (stats::median(En) / stats::median(base_E) - 1), port_change = stats::median(Cn) - stats::median(base_C),
    port_change_pct = 100 * (stats::median(Cn) - stats::median(base_C)) / PORT_B70, stringsAsFactors = FALSE)
  months[[length(months) + 1L]] <<- data.frame(variant = label, month = sort(unique(mon)),
    reported = sapply(sort(unique(mon)), function(m) stats::median(rowSums(Cd[, mon == m, drop = FALSE]))),
    priced = sapply(sort(unique(mon)), function(m) stats::median(rowSums((Cd * ratio)[, mon == m, drop = FALSE]))), stringsAsFactors = FALSE)
}
one("as reported (the share at the model's day volume)", matrix(1, nd, D))
rP <- integ_ratio(sg); one("PRIMARY: the share integrated over the day volume, s_day = sigma_uE", rP)
one("sensitivity: s_day = 0.5", integ_ratio(rep(0.5, nd)))
one("sensitivity: s_day = 0.82", integ_ratio(rep(0.82, nd)))
# beside it: B86's switch, volume-consistent, on the sampled trailer days without an OSP count; F_d elsewhere
Tn <- as.integer(sd_$T_n); tday <- as.integer(sd_$day_T)[seq_len(Tn)]; ted <- as.integer(sd_$Eday_T)[seq_len(Tn)]
osp_days <- unique(as.integer(sd_$day_OSP)[seq_len(sd_$OSP_n)]); carry <- setdiff(unique(tday), osp_days)
rS <- rP; rC <- matrix(1, nd, D)
for (d in carry) { df <- exp(sg * uu[, ted[match(d, tday)]] - 0.5 * sg^2)
  r <- df * share(vday[, d] * df, fstr[d]) / f0[, d]; rS[, d] <- r; rC[, d] <- r }
one("B86's switch alone (volume-consistent), every other day as reported", rC)
one("B86's switch on the carried days, the integrated share on every other day", rS)
res <- do.call(rbind, rows); mres <- do.call(rbind, months); mres$change_pct <- 100 * (mres$priced / mres$reported - 1)
op <- options(width = 220); print(res, row.names = FALSE, digits = 4)
print(mres[mres$variant == res$variant[2], c("month", "reported", "priced", "change_pct")], row.names = FALSE, digits = 4); options(op)
pc <- res$change_pct[2]
reading <- if (abs(pc) > 1) sprintf("the PRIMARY change is %.1f%% of the boat all-gear catch, above 1%%: D56 is opened as an open item that moves the number, a decision for Matt", pc) else
  sprintf("the PRIMARY change is %.1f%%, within 1%%: D56 is opened as low", pc)
cat("\nB87:", reading, "\n")
res$reading <- c("", reading, rep("", nrow(res) - 2L))
utils::write.csv(res, OUT, row.names = FALSE); utils::write.csv(mres, OUT_M, row.names = FALSE)
cat("Wrote", OUT, "and", OUT_M, "\n")
