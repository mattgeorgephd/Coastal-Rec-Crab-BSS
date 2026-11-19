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
# desk_advisory_crab_share_2026-11-19.R
#
# DESK, SECONDS, NO MCMC (CHANGE_REGISTER B81, for D30; a rider on B80).
#
# THE QUESTION. The boat advisory term is fitted to counts of ALL private boats, and the
# estimate needs CRABBING boats. D30 (2026-09-25) asked whether crabbing boats stay in port
# on an advisory day at a different rate from the rest of the fleet, and named what would
# settle it: OSP's classification on advisory and other days. That classification arrived on
# 2026-10-02 (B64), and since A34 the model's crabbing share varies with the day's boat
# volume, so the model already implies a differential: on the authoritative run's
# coefficients an advisory day has 0.32 times the boats and, through the volume term
# (beta -0.76), a higher crabbing share. Nobody has asked the OSP days whether the share on
# an advisory day is what that volume term predicts.
#
# THE DATA. The authoritative run's committed osp_crab_only_daily.csv
# (05_output/20261019/pooled-CPUE-2024-25): per OSP day, the boats returning (V), the boats
# OSP sampled (n, the schedule's minimum rate, D2) and the crab-only boats among them (y).
# Only the days on or after the pot opening are read, because that is where the term is
# applied; they run March to September (OSP does not sample the winter).
#
# THE RULE, WRITTEN BEFORE THE REGRESSION WAS RUN (only the table's header and three rows
# had been read, to learn its columns). R9, exploratory: it judges nothing in B80.
#   Two quasi-binomial logit regressions of y out of n, calendar month as a factor:
#     Ma   month + advisory                 the raw contrast
#     Mb   month + log(V) + advisory        beyond what the day's volume predicts, which is
#                                           the form of the model's own volume term (A34)
#   Mb's advisory coefficient, 95% interval from the quasi-binomial standard error:
#     excludes 0   the crab-only share on an advisory day differs from what the day's volume
#                  predicts. POSITIVE: crab-only boats are deterred less than the volume term
#                  allows, so the model still over-corrects crab effort on advisory days.
#                  NEGATIVE: the reverse.
#     covers 0     no differential beyond volume is detectable on these days; whatever
#                  differential there is, the volume term is carrying it.
#   Fewer than 10 advisory days among them: too few to read.
#   Ma is reported beside it. Both are repeated with n = V, the other bound of D2 (every boat
#   sampled), and the reading is Mb's on the schedule n; a sign change between the two bounds
#   is reported as such and the reading withdrawn.
#   STATED IN ADVANCE. Crab-only is a LOWER BOUND on crabbing boats (combo trips are
#   classed by their other fishery), so a shift could be the combo share's and not the
#   crabbing share's. There is no December to February OSP day, so this says nothing about
#   the winter, where D32 lives. The 2025-26 OSP days are all in the pot closure, where the
#   term is not applied: they are counted and not read. Nothing changes on this result.
#   EXPECTED: Ma positive (fewer boats, a higher crabbing share); Mb's interval covering 0.
#
# WRITES 05_output/advisory_replication_2026-11-19_crab_share.csv
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
suppressWarnings(suppressPackageStartupMessages(try({ library(dplyr); library(tidyr); library(lubridate); library(readxl); library(here); library(purrr); library(stringr); library(tibble) }, silent = TRUE)))
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
source("run_config.R")

SRC <- file.path("05_output", "20261019", "pooled-CPUE-2024-25", "osp_crab_only_daily.csv")
SRC26 <- file.path("05_output", "20261020", "pooled-CPUE-season-2025-26", "osp_crab_only_daily.csv")
OUT <- file.path("05_output", "advisory_replication_2026-11-19_crab_share.csv")
MIN_ADVISORY_DAYS <- 10L

d <- utils::read.csv(SRC, stringsAsFactors = FALSE)
d$event_date <- as.Date(d$event_date)
fl <- marine_hazard_flag_series(run_config, date_start = min(d$event_date), date_end = max(d$event_date))
d$advisory <- as.integer(fl$nws_sca_any[match(d$event_date, as.Date(fl$event_date))])
d <- d[d$event_date >= as.Date(run_config$pot_open_date) & d$osp_boat_total > 0 & d$osp_total > 0, , drop = FALSE]
d$month <- factor(format(d$event_date, "%Y-%m"))
cat(sprintf("2024-25 all-gear OSP days with a crab-only count: %d (%s to %s); advisory days %d; boats returning %d, sampled %d, crab-only %d\n",
            nrow(d), min(d$event_date), max(d$event_date), sum(d$advisory), sum(d$osp_boat_total), sum(d$osp_total), sum(d$osp_crab_only)))
by_month <- d |> group_by(month) |> summarise(days = n(), advisory_days = sum(advisory), .groups = "drop")
print(as.data.frame(by_month))

fit_one <- function(basis) {
  n <- if (basis == "schedule minimum") d$osp_total else d$osp_boat_total
  y <- pmin(d$osp_crab_only, n)
  dd <- data.frame(y = y, n = n, month = d$month, lv = log(d$osp_boat_total), advisory = d$advisory)
  row <- function(label, f) {
    m <- stats::glm(f, family = stats::quasibinomial(), data = dd); cf <- summary(m)$coefficients
    b <- cf["advisory", "Estimate"]; se <- cf["advisory", "Std. Error"]
    data.frame(season = "2024-25", model = label, n_basis = basis, days = nrow(dd), advisory_days = sum(dd$advisory),
               advisory_coef = b, se = se, lo95 = b - 1.96 * se, hi95 = b + 1.96 * se, odds_ratio = exp(b), or_lo95 = exp(b - 1.96 * se), or_hi95 = exp(b + 1.96 * se),
               log_volume_coef = if ("lv" %in% rownames(cf)) cf["lv", "Estimate"] else NA_real_, log_volume_se = if ("lv" %in% rownames(cf)) cf["lv", "Std. Error"] else NA_real_,
               dispersion = summary(m)$dispersion,
               share_advisory_days = sum(dd$y[dd$advisory == 1]) / sum(dd$n[dd$advisory == 1]), share_other_days = sum(dd$y[dd$advisory == 0]) / sum(dd$n[dd$advisory == 0]),
               stringsAsFactors = FALSE)
  }
  rbind(row("Ma: month + advisory", cbind(y, n - y) ~ month + advisory), row("Mb: month + log(V) + advisory", cbind(y, n - y) ~ month + lv + advisory))
}
res <- rbind(fit_one("schedule minimum"), fit_one("every boat sampled"))
mb <- res[res$model == "Mb: month + log(V) + advisory", ]
main <- mb[mb$n_basis == "schedule minimum", ]; other <- mb[mb$n_basis == "every boat sampled", ]
res$reading <- ""
reading <- if (main$advisory_days < MIN_ADVISORY_DAYS) "too few advisory days to read" else
  if (sign(main$advisory_coef) != sign(other$advisory_coef) && (other$lo95 > 0 || other$hi95 < 0 || main$lo95 > 0 || main$hi95 < 0)) "the two bounds of D2 disagree in sign: the reading is withdrawn" else
  if (main$lo95 > 0) "the crab-only share on an advisory day is HIGHER than the day's volume predicts: crab-only boats are deterred less than the volume term allows" else
  if (main$hi95 < 0) "the crab-only share on an advisory day is LOWER than the day's volume predicts" else
  "no differential beyond volume is detectable on these days"
res$reading[res$model == "Mb: month + log(V) + advisory" & res$n_basis == "schedule minimum"] <- reading
op <- options(width = 220); print(res[, c("model", "n_basis", "days", "advisory_days", "advisory_coef", "se", "lo95", "hi95", "odds_ratio", "log_volume_coef", "dispersion", "share_advisory_days", "share_other_days")], digits = 3, row.names = FALSE); options(op)
cat("\nR9 (Mb on the schedule n):", reading, "\n")
if (file.exists(SRC26)) { d26 <- utils::read.csv(SRC26, stringsAsFactors = FALSE)
  cat(sprintf("2025-26 OSP days with a crab-only count: %d, all in the pot closure (%s to %s): counted, not read.\n", nrow(d26), min(d26$event_date), max(d26$event_date))) }
utils::write.csv(res, OUT, row.names = FALSE)
cat("Wrote", OUT, "\n")
