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
# desk_ppc_day_effect_2026-11-27.R
#
# DESK, MINUTES, NO MCMC (CHANGE_REGISTER B90, for D58). Found in the housekeeping review of 2026-11-27.
#
# THE QUESTION. bss_ppc_calibration() (model_diagnostics.R), pit_block() (save_run_diagnostics.R) and
# .eod_decompose() (diagnose_effort_overdispersion.R) form an effort count's predictive mean as
# lambda_E[day] * m, m the stream's expansion multiplier. Since A39 and A40 (2026-10-18/19) the
# likelihood of a fit that carries the per-day effort effect has the mean
#   lambda_E[day] * m * exp(sigma_uE * uE_raw[eday] - sigma_uE^2 / 2),
# and none of the three carries that factor. On the authoritative run the boat all-gear fit's trailer
# stream reads coverage_50 0.359 and pit_mean 0.398 in ppc_calibration_*, and model_adequacy.csv raises
# flag_pit_bias on that fit; before A39 (B65, B67) the same stream read about 0.57 and no flag. Is the
# flag the model's or the diagnostic's?
#
# THE MEASUREMENT, on B89's two all-gear fits (the fits that carry the effect), from the saved
# ppc_draws_*.rds (every 5th draw), with the pipeline's own mid-PIT (bss_zi_pit) and statistics:
#   (a) as_computed   mu = lambda_E * m                          what the pipeline writes today
#   (b) conditional   mu = lambda_E * m * the day's fitted factor the likelihood's own mean
#   (c) integrated    mu = lambda_E * m * exp(sigma_uE z - sigma_uE^2 / 2), z ~ N(0, 1) fresh,
#                     20 replicates per draw (seed 1): the predictive of a count on a NEW day
#
# THE READING RULE, WRITTEN AND COMMITTED BEFORE THE SCRIPT WAS RUN.
#   CHECK   (a) must reproduce the committed ppc_calibration_* coverage_50 of the stream within 0.02
#           (the pipeline scores 400 seeded draws, this script 2,000), or nothing is read.
#   READ    with the pipeline's own thresholds (bss_model_adequacy(): |pit_mean - 0.5| > 0.05 is the
#           PIT-bias flag, |coverage_50 - 0.5| > 0.15 the miscalibration flag), per stream:
#     the diagnostic's   (a) raises a flag and (c) raises none: the flag comes from the omitted factor.
#     the model's        (c) raises the flag (a) raises: the per-day effect does not repair it.
#     mixed              anything else, reported as measured.
#   STATED IN ADVANCE. (1) (b) conditions on a deviation fitted to the same counts, so it is expected to
#   OVER-cover (coverage_50 above 0.5); it is shown and not read. (2) (c) treats the counts of one day as
#   separate new days, which leaves its expectation right and understates its sampling error; the
#   per-stream sampling SD printed is sqrt(0.25 / days), on days, not counts. (3) This judges the
#   diagnostic, not the model's adoption: A39 and A40 were adopted on their own rules. (4) Nothing in
#   the pipeline changes on this result; the repair is D58's and Matt's to choose.
#   EXPECTED: "the diagnostic's" on both streams.
#
# HOW TO RUN  RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_ppc_day_effect_2026-11-27.R [pooled run folder]
#             (needs the run folder's ppc_draws_*.rds, which are not committed; default B89's folder)
# WRITES      05_output/ppc_day_effect_2026-11-27.csv
###############################################################################
.args <- commandArgs(trailingOnly = TRUE)
DIR <- if (length(.args) >= 1L) .args[1] else "05_output/20261125/pooled-CPUE-2024-25-190951"
source("03_R_functions/zinb_ppc.R")
FITS <- c("shore_all_gear_Dungeness_Kept", "private_boat_all_gear_Dungeness_Kept")
stat <- function(mu, y, size) {
  pit <- vapply(seq_along(y), function(i) { m <- pmax(mu[, i], 1e-8); k <- is.finite(m) & is.finite(size); bss_zi_pit(y[i], m[k], size[k]) }, numeric(1))
  c(coverage_50 = mean(pit >= 0.25 & pit <= 0.75), coverage_95 = mean(pit >= 0.025 & pit <= 0.975), pit_mean = mean(pit), pit_sd = stats::sd(pit))
}
rows <- list()
for (fit in FITS) {
  f <- file.path(DIR, sprintf("ppc_draws_%s.rds", fit)); if (!file.exists(f)) stop("missing ", f)
  x <- readRDS(f); sd_ <- x$stan_data; dr <- x$draws
  if (!isTRUE(sd_[["effort_day_effect"]] == 1L)) stop(fit, ": the per-day effort effect is off in this fit")
  n <- nrow(dr$uE_raw); use <- seq(1L, n, by = 5L); nd <- length(use)
  lam <- dr$lambda_E_S; d <- dim(lam)
  lam <- matrix(if (length(d) == 4L) lam[use, 1, , 1] else if (length(d) == 3L) lam[use, , 1] else lam[use, ], nrow = nd)
  rE <- as.numeric(dr$r_E)[use]; s <- as.numeric(dr$sigma_uE)[use]; U <- dr$uE_raw[use, , drop = FALSE]
  committed <- utils::read.csv(file.path(DIR, sprintf("ppc_calibration_%s.csv", fit)), stringsAsFactors = FALSE)
  streams <- list()
  if (sd_$Gear_n > 0) streams$gear <- list(days = sd_$day_Gear, y = sd_$Gear_I, m = as.numeric(dr$R_G)[use], eday = sd_$Eday_Gear)
  if (sd_$T_n > 0 && !is.null(dr$R_G_boat)) streams$trailer <- list(days = sd_$day_T, y = sd_$T_I, m = 1 / as.numeric(dr$R_G_boat)[use], eday = sd_$Eday_T)
  for (st in names(streams)) {
    q <- streams[[st]]; mu0 <- lam[, q$days, drop = FALSE] * q$m
    fac <- vapply(q$eday, function(e) if (e > 0L) exp(s * U[, e] - 0.5 * s^2) else rep(1, nd), numeric(nd))
    set.seed(1); K <- 20L
    mu_int <- do.call(rbind, lapply(seq_len(K), function(k) mu0 * exp(s * matrix(stats::rnorm(nd * ncol(mu0)), nd) - 0.5 * s^2)))
    res <- rbind(as_computed = stat(mu0, q$y, rE), conditional = stat(mu0 * fac, q$y, rE), integrated = stat(mu_int, q$y, rep(rE, K)))
    days <- length(unique(q$days)); com <- committed$coverage_50[committed$data_type == st]
    flag <- function(r) paste(c(if (abs(r["pit_mean"] - 0.5) > 0.05) "pit_bias", if (abs(r["coverage_50"] - 0.5) > 0.15) "miscalibrated"), collapse = "+")
    fl <- apply(res, 1, flag); chk <- abs(res["as_computed", "coverage_50"] - com) <= 0.02
    reading <- if (!chk) "CHECK FAILED: (a) does not reproduce the committed coverage_50; nothing is read" else
      if (nzchar(fl["as_computed"]) && !nzchar(fl["integrated"])) "the diagnostic's: the flag comes from the omitted day factor" else
      if (nzchar(fl["as_computed"]) && nzchar(fl["integrated"]) && all(strsplit(fl["as_computed"], "+", fixed = TRUE)[[1]] %in% strsplit(fl["integrated"], "+", fixed = TRUE)[[1]])) "the model's: the per-day effect does not repair it" else "mixed: reported as measured"
    rows[[length(rows) + 1L]] <- data.frame(fit = fit, stream = st, n_counts = length(q$y), n_days = days, sampling_sd_days = sqrt(0.25 / days),
                                            committed_coverage_50 = com, mean_form = rownames(res), res, flags = ifelse(nzchar(fl), fl, "none"), reading = reading, row.names = NULL, stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
op <- options(width = 250); print(out[, setdiff(names(out), "reading")], row.names = FALSE, digits = 3); options(op)
for (k in unique(paste(out$fit, out$stream))) cat(k, ":", out$reading[paste(out$fit, out$stream) == k][1], "\n")
utils::write.csv(out, "05_output/ppc_day_effect_2026-11-27.csv", row.names = FALSE); cat("Wrote 05_output/ppc_day_effect_2026-11-27.csv\n")
