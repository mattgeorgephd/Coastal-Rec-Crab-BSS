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
# desk_crab_share_boat_weighted_osp_2026-11-25.R
#
# DESK, MINUTES, NO MCMC (CHANGE_REGISTER B88, for D56; the model-light check of B87).
#
# THE QUESTION. B87 priced, through the fitted model, how far the crab effort the totals credit
# exceeds a day's expected crab effort once the crabbing share is weighted over the real
# day-to-day variation in boat volume. That rests on the model. The OSP days can answer the same
# question with almost none: on each OSP day, OSP counted the boats returning (V_d) and classified
# a sample of them (n_d sampled, y_d crab-only). D40's target, the share the totals need, is the
# BOAT-WEIGHTED crab-only share, sum_d V_d (y_d / n_d) / sum_d V_d. The Stan comment at f_volume
# says the share "as the totals apply it is BOAT-WEIGHTED". Here that is checked where it can be.
#
# THE COMPARISON, month by month over the 2024-25 OSP days in the all-gear window, on the
# authoritative configuration's boat all-gear fit (B80's S25-const, byte-identical to B70):
#   observed  sum_d V_d (y_d / n_d) / sum_d V_d                         (OSP's own data)
#   applied   sum_d v_d p_d / sum_d v_d, per draw, then its median     (what the totals apply)
#             v_d = lambda_E[d] L[d] / R_G_boat, the model's day volume as the generated
#             quantities form it; p_d = f(v_d) (1 - combo_c[month]), the crab-only share at it.
#   ratio     applied / observed, with a 95% interval from a day bootstrap of the observed side
#             (2,000 resamples within month, seed 1) combined with the posterior of the applied side.
#
# THE RULE, WRITTEN AND COMMITTED BEFORE THE COMPARISON WAS COMPUTED.
#   Over all the OSP days together (each month's ratio weighted by its boats):
#     the ratio's interval lies above 1   the totals apply a crab-only share above the boat-
#                                         weighted one on the days where it can be checked: B87's
#                                         direction is confirmed by the data, and the excess is its
#                                         model-light size there.
#     the interval covers 1               the data do not show the excess B87 priced; D56 is
#                                         reported as a model calculation the OSP days do not confirm.
#     the interval lies below 1           the opposite direction; reported as such.
#   Months with fewer than 5 OSP days are shown and not read on their own.
#   STATED IN ADVANCE. (1) Crab-only is a lower bound on crabbing boats (combo trips are classed by
#   their other fishery); the comparison is like for like because the applied side carries the same
#   (1 - combo_c). (2) n_d is the schedule's minimum sampled count (D2), as the model reads it.
#   (3) It covers March to September, where OSP samples; December to February is not checked.
#   EXPECTED: the ratio above 1 by several percent in the summer months.
#
# WRITES 05_output/crab_share_boat_weighted_osp_2026-11-25.csv
#
# HOW TO RUN (needs the S25-const run folder and the fit cache)
#   BSS_FIT_CACHE_DIR=<dir> RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_crab_share_boat_weighted_osp_2026-11-25.R
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
OSPD <- file.path("05_output", "20261019", "pooled-CPUE-2024-25", "osp_crab_only_daily.csv")
OUT <- "05_output/crab_share_boat_weighted_osp_2026-11-25.csv"
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
day1 <- as.Date(pw$event_date[1]) - (pw$day_index[1] - 1L)
dates <- day1 + seq_len(D) - 1L
stopifnot(all(dates[pw$day_index] == as.Date(pw$event_date)))
lam <- mat("lambda_E_S"); li <- do.call(rbind, lapply(strsplit(sub("^lambda_E_S\\[(.*)\\]$", "\\1", colnames(lam)), ","), as.integer))
Lm <- mat("L_out"); rg <- mat("R_G_boat")[, 1]; etaf <- mat("eta_f"); bf <- mat("beta_fvol")[, 1]; cc <- mat("combo_c")
fstr <- as.integer(sd_$f_stratum); cen <- as.numeric(sd_$fvol_centre)

o <- utils::read.csv(OSPD, stringsAsFactors = FALSE); o$event_date <- as.Date(o$event_date)
o <- o[o$event_date >= min(dates) & o$event_date <= max(dates) & o$osp_total > 0 & o$osp_boat_total > 0, ]
o$d <- match(o$event_date, dates); o$month <- format(o$event_date, "%Y-%m"); o$k <- fstr[o$d]
cat(sprintf("OSP days with a crab-only count in the window: %d (%s to %s)\n", nrow(o), min(o$event_date), max(o$event_date)))
v <- sapply(o$d, function(d) rowSums(lam[, li[, 2] == d & li[, 3] == Gd, drop = FALSE]) / rg * Lm[, d])     # [draws x days]
p <- sapply(seq_len(nrow(o)), function(j) { k <- o$k[j]; (1e-6 + (1 - 2e-6) * stats::plogis(etaf[, k] + bf * (log(pmax(v[, j], 1e-9)) - cen[k]))) * (1 - cc[, k]) })
y <- pmin(o$osp_crab_only, o$osp_total); n <- o$osp_total; V <- o$osp_boat_total
set.seed(1); B <- 2000L
months <- sort(unique(o$month))
rows <- lapply(c(months, "all"), function(m) {
  j <- if (m == "all") seq_len(nrow(o)) else which(o$month == m)
  obs <- sum(V[j] * y[j] / n[j]) / sum(V[j])
  day_avg <- mean(y[j] / n[j])
  app <- rowSums(v[, j, drop = FALSE] * p[, j, drop = FALSE]) / rowSums(v[, j, drop = FALSE])
  # the observed side's bootstrap, within month for "all"
  bs <- replicate(B, { jj <- if (m == "all") unlist(lapply(split(j, o$month[j]), function(q) q[sample.int(length(q), replace = TRUE)])) else j[sample.int(length(j), replace = TRUE)]
                       sum(V[jj] * y[jj] / n[jj]) / sum(V[jj]) })
  ratio <- app[sample.int(length(app), B, replace = TRUE)] / bs
  data.frame(month = m, osp_days = length(j), boats = sum(V[j]), observed_boat_weighted = obs, observed_day_average = day_avg,
             applied_median = stats::median(app), ratio_median = stats::median(ratio),
             ratio_lo95 = unname(stats::quantile(ratio, 0.025)), ratio_hi95 = unname(stats::quantile(ratio, 0.975)),
             readable = length(j) >= 5L || m == "all", stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)
a <- res[res$month == "all", ]
reading <- if (a$ratio_lo95 > 1) sprintf("the totals apply a crab-only share %.0f%% above the boat-weighted one on the OSP days [%.0f%%, %.0f%%]: B87's direction is confirmed by the data", 100 * (a$ratio_median - 1), 100 * (a$ratio_lo95 - 1), 100 * (a$ratio_hi95 - 1)) else
  if (a$ratio_hi95 < 1) "the totals apply a crab-only share BELOW the boat-weighted one on the OSP days: the opposite of B87's direction" else
  "the data do not show the excess B87 priced: D56 stays a model calculation the OSP days do not confirm"
res$reading <- ifelse(res$month == "all", reading, "")
op <- options(width = 220); print(res[, 1:10], row.names = FALSE, digits = 3); options(op)
cat("\nB88:", reading, "\n")
utils::write.csv(res, OUT, row.names = FALSE); cat("Wrote", OUT, "\n")
