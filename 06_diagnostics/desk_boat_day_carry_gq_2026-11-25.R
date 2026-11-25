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
# desk_boat_day_carry_gq_2026-11-25.R
#
# DESK, MINUTES, NO MCMC (CHANGE_REGISTER B86, the D55 switch).
#
# WHAT IS CHECKED. run_config.R's boat_day_effect_reporting becomes the Stan data flag
# boat_day_carry (both programs; read only in transformed data, where it adds the sampled trailer
# days with no OSP count to eday_of_day, so the generated quantities carry their fitted day
# deviation). Generated quantities do not feed back into a fit, so the switch can be checked
# EXACTLY without sampling: the authoritative configuration's boat all-gear fit (the B80 arm
# S25-const, byte-identical to B70) is replayed through the new pooled program's generated
# quantities with rstan::gqs(), once with the flag 0 and once with 1.
#
# THE RULE, WRITTEN AND COMMITTED BEFORE THE REPLAY RAN.
#   (1) flag 0: C_expected_sum and E_sum equal the fit's own draws within 1e-9 relative, every draw
#       (the switch off is the shipped model).
#   (2) flag 1: the median of C_expected_sum equals B80's exact what-if for "the deviation carried on
#       sampled trailer days with no OSP count" (05_output/advisory_replication_2026-11-19_arms/
#       S25-const/boat_reporting_whatif.csv, 34,686) within 0.5 crab, and the median of E_sum its
#       effort column the same way. The what-if computed the same product in R; agreement says the
#       Stan code does what the register says it does.
#   (3) flag 1: no OSP day carries a deviation (the effort on every OSP day is unchanged, every draw).
#   ACCEPTED if (1) to (3) hold.
#
# READ 2026-11-25 (the switch as first built, 2442411): ACCEPTED. Flag 0 equals the fit to 1e-15 on
# every draw, flag 1 gives 34,686.44 against the what-if's 34,686.44 (-2.19%), and no OSP day moves.
# (The first read failed (1) and (3) on this script's draw order alone: gqs() follows as.matrix(),
# extract() permutes. Corrected and re-run; kept as 05_output/boat_day_carry_gq_2026-11-25_fixed_share.csv.)
#
# ADDENDUM (2b), WRITTEN AFTER THAT READ AND BEFORE THE STAN CHANGE IT CHECKS. The first form holds
# a carried day's crabbing share at the LEVEL's boat volume. A34's volume term was fitted on OSP
# days, where the volume is the day's actual boats, so on a carried day the consistent share is the
# one at the carried volume (level x the day factor); the share falls as the boats rise (beta_fvol
# -0.76 on B70), so holding it overstates how far crab effort moves with the deviation. The register
# had said so at D55 ("the day's volume in the crabbing share (A34) has to be decided with it"), and
# B80's what-if made the same simplification. The switch is changed so the generated quantities use
# vd * L[d] * df for the share (df is exactly 1.0 on every day that does not carry a deviation, so
# flag 0 stays bit-identical). Clause (2) is then replaced by:
#   (2b) flag 1: C_expected_sum and E_sum equal, within 1e-9 relative on every draw, an independent
#        recomputation in R from the fit's own draws (E, C_expected, L_out, lambda_E_S, R_G_boat,
#        eta_f, beta_fvol, sigma_uE, uE_raw and the Stan data): on each carried day the effort and
#        catch are multiplied by df * f(carried volume) / f(level volume).
#   Clauses (1) and (3) stand. The first form's what-if row is reported beside it, not judged.
#
# WRITES 05_output/boat_day_carry_gq_2026-11-25.csv
#
# HOW TO RUN (needs the S25-const run folder and the fit cache; compiles the pooled program once)
#   BSS_FIT_CACHE_DIR=<dir> RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_boat_day_carry_gq_2026-11-25.R
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
WHATIF <- "05_output/advisory_replication_2026-11-19_arms/S25-const/boat_reporting_whatif.csv"
OUT <- "05_output/boat_day_carry_gq_2026-11-25.csv"
CACHE <- Sys.getenv("BSS_FIT_CACHE_DIR", unset = "")
if (!nzchar(CACHE) || !dir.exists(CACHE)) stop("BSS_FIT_CACHE_DIR is not set or is not a directory")

sd_ <- readRDS(file.path(DIR, sprintf("ppc_draws_%s.rds", FIT)))$stan_data
sd_ <- sd_[!startsWith(names(sd_), ".")]
if (is.null(sd_[["sunday"]]) && !is.null(sd_[["sun"]])) { sd_[["sunday"]] <- sd_[["sun"]]; sd_[["sun"]] <- NULL }   # saved before B84
pw <- utils::read.csv(file.path(DIR, sprintf("loo_pointwise_trailer_%s.csv", FIT)))
fit <- NULL
for (f in Sys.glob(file.path(CACHE, sprintf("%s_att*_*.rds", FIT)))) {
  x <- readRDS(f); ll <- as.matrix(rstan::extract(x, pars = "log_lik_trailer")[[1]])
  if (ncol(ll) == nrow(pw) && max(abs(apply(ll, 2, .bcv_lmeanexp) - pw$lpd)) < 5e-4) { fit <- x; cat("fit:", basename(f), "\n"); break }
  rm(x); gc(verbose = FALSE)
}
if (is.null(fit)) stop("no cached fit reproduces S25-const")
m <- rstan::stan_model("02_stan_models/crab_bss_pooled.stan", auto_write = FALSE)
pars <- fit@model_pars[seq_len(which(fit@model_pars == "lp__") - 1L)]
# the parameters block only: gqs() takes draws of the parameters, by their constrained names
pblock <- names(rstan::get_inits(fit)[[1]])
dm <- as.matrix(fit, pars = pblock)
# THE DRAW ORDER (read 2026-11-25): gqs() returns its draws in the order of `dm`, which is as.matrix()'s
# (chain by chain, unpermuted), while rstan::extract() permutes them; the first read compared the two
# orders draw by draw and failed clauses (1) and (3) on that alone. Everything is now taken in
# as.matrix() order on both sides.
am <- function(f, p) as.matrix(f, pars = p)
day_sum <- function(f) { x <- am(f, "E"); d <- as.integer(sub("^E\\[[0-9]+,([0-9]+),[0-9]+\\]$", "\\1", colnames(x)))
  sapply(sort(unique(d)), function(k) rowSums(x[, d == k, drop = FALSE])) }
E0 <- list(C_expected_sum = am(fit, "C_expected_sum")[, 1], E_sum = am(fit, "E_sum")[, 1], Eday = day_sum(fit))
run <- function(flag) { d <- sd_; d$boat_day_carry <- as.integer(flag); rstan::gqs(m, data = d, draws = dm, seed = 1) }
g0 <- run(0); g1 <- run(1)
x0 <- list(C_expected_sum = am(g0, "C_expected_sum")[, 1], E_sum = am(g0, "E_sum")[, 1])
x1 <- list(C_expected_sum = am(g1, "C_expected_sum")[, 1], E_sum = am(g1, "E_sum")[, 1], Eday = day_sum(g1))
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
# (2b): the independent recomputation in R, every draw in as.matrix() order
mat <- function(p) as.matrix(fit, pars = p)
nd <- nrow(dm); D <- as.integer(sd_$D); S <- as.integer(sd_$S %||% 1L); Gd <- as.integer(sd_$G %||% 1L)
lam <- mat("lambda_E_S")                                   # columns lambda_E_S[s,d,g]
lam_idx <- do.call(rbind, lapply(strsplit(sub("^lambda_E_S\\[(.*)\\]$", "\\1", colnames(lam)), ","), as.integer))
Lm <- mat("L_out"); rg <- mat("R_G_boat")[, 1]; etaf <- mat("eta_f"); bf <- mat("beta_fvol")[, 1]
sg <- mat("sigma_uE")[, 1]; uu <- mat("uE_raw")
Tn <- as.integer(sd_$T_n); tday <- as.integer(sd_$day_T)[seq_len(Tn)]; ted <- as.integer(sd_$Eday_T)[seq_len(Tn)]
osp_days <- unique(as.integer(sd_$day_OSP)[seq_len(sd_$OSP_n)])
carry_days <- setdiff(unique(tday), osp_days)
fstr <- as.integer(sd_$f_stratum); cen <- as.numeric(sd_$fvol_centre)
share <- function(v, k, s) 1e-6 + (1 - 2e-6) * stats::plogis(etaf[s, k] + bf[s] * (log(v) - cen[k]))
Eday_fit <- day_sum(fit); Cday_fit <- { x <- mat("C_expected"); d <- as.integer(sub("^C_expected\\[[0-9]+,([0-9]+),[0-9]+\\]$", "\\1", colnames(x)))
  sapply(sort(unique(d)), function(k) rowSums(x[, d == k, drop = FALSE])) }
Er <- Eday_fit; Cr <- Cday_fit
for (d in carry_days) {
  ed <- ted[match(d, tday)]; df <- exp(sg * uu[, ed] - 0.5 * sg^2)
  vd <- rowSums(lam[, lam_idx[, 2] == d & lam_idx[, 3] == Gd, drop = FALSE]) / rg
  k <- fstr[d]; v0 <- pmax(vd * Lm[, d], 1e-9); v1 <- pmax(vd * Lm[, d] * df, 1e-9)
  ratio <- df * share(v1, k, seq_len(nd)) / share(v0, k, seq_len(nd))
  Er[, d] <- Er[, d] * ratio; Cr[, d] <- Cr[, d] * ratio
}
R_E <- rowSums(Er); R_C <- rowSums(Cr)
wi <- utils::read.csv(WHATIF); w <- wi[grepl("no OSP count", wi$reporting), ]
osp_days <- unique(as.integer(sd_$day_OSP)[seq_len(sd_$OSP_n)])
Eday0 <- E0$Eday; Eday1 <- x1$Eday
res <- data.frame(
  clause = c("1", "1", "2b", "2b", "3", "reported"),
  what = c("flag 0: C_expected_sum equals the fit's, every draw (max relative gap)", "flag 0: E_sum equals the fit's, every draw (max relative gap)",
           "flag 1: C_expected_sum equals the R recomputation, every draw (max relative gap)", "flag 1: E_sum equals the R recomputation, every draw (max relative gap)",
           "flag 1: effort on the OSP days unchanged, every draw (max relative gap)",
           "the first form's price (share held at the level's volume, B80's what-if), median C_expected_sum"),
  observed = c(rel(x0$C_expected_sum, E0$C_expected_sum), rel(x0$E_sum, E0$E_sum),
               rel(x1$C_expected_sum, R_C), rel(x1$E_sum, R_E), rel(Eday1[, osp_days], Eday0[, osp_days]), NA),
  value = c(NA, NA, stats::median(x1$C_expected_sum), stats::median(x1$E_sum), NA, w$expected_catch_median),
  stringsAsFactors = FALSE)
res$pass <- c(res$observed[1:2] <= 1e-9, res$observed[3:4] <= 1e-9, res$observed[5] <= 1e-9, NA)
res$verdict <- if (all(res$pass, na.rm = TRUE)) "ACCEPTED" else "NOT ACCEPTED"
cat(sprintf("flag 0 median C_expected_sum %.1f (fit %.1f); flag 1 %.1f (R recomputation %.1f), %.2f%%; the first form's what-if %.1f (%.2f%%)\n", stats::median(x0$C_expected_sum), stats::median(E0$C_expected_sum),
            stats::median(x1$C_expected_sum), stats::median(R_C), 100 * (stats::median(x1$C_expected_sum) / stats::median(x0$C_expected_sum) - 1),
            w$expected_catch_median, 100 * (w$expected_catch_median / stats::median(x0$C_expected_sum) - 1)))
cat(sprintf("flag 1 boat all-gear effort %.1f against %.1f (%.2f%%); carried days %d\n", stats::median(x1$E_sum), stats::median(x0$E_sum), 100 * (stats::median(x1$E_sum) / stats::median(x0$E_sum) - 1), length(carry_days)))
print(res, row.names = FALSE); cat("\nB86 switch:", res$verdict[1], "\n")
utils::write.csv(res, OUT, row.names = FALSE); cat("Wrote", OUT, "\n")
