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
E0 <- rstan::extract(fit, pars = c("C_expected_sum", "E_sum", "E"))
run <- function(flag) { d <- sd_; d$boat_day_carry <- as.integer(flag); rstan::gqs(m, data = d, draws = dm, seed = 1) }
g0 <- run(0); g1 <- run(1)
x0 <- rstan::extract(g0, pars = c("C_expected_sum", "E_sum")); x1 <- rstan::extract(g1, pars = c("C_expected_sum", "E_sum", "E"))
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
wi <- utils::read.csv(WHATIF); w <- wi[grepl("no OSP count", wi$reporting), ]
osp_days <- unique(as.integer(sd_$day_OSP)[seq_len(sd_$OSP_n)])
Eday0 <- apply(E0$E, c(1, 3), sum); Eday1 <- apply(x1$E, c(1, 3), sum)
res <- data.frame(
  clause = c("1", "1", "2", "2", "3"),
  what = c("flag 0: C_expected_sum equals the fit's, every draw (max relative gap)", "flag 0: E_sum equals the fit's, every draw (max relative gap)",
           "flag 1: median C_expected_sum against the what-if", "flag 1: median E_sum against the what-if", "flag 1: effort on the OSP days unchanged, every draw (max relative gap)"),
  observed = c(rel(x0$C_expected_sum, E0$C_expected_sum), rel(x0$E_sum, E0$E_sum),
               stats::median(x1$C_expected_sum) - w$expected_catch_median, stats::median(x1$E_sum) - w$effort_median,
               rel(Eday1[, osp_days], Eday0[, osp_days])),
  value = c(NA, NA, stats::median(x1$C_expected_sum), stats::median(x1$E_sum), NA),
  stringsAsFactors = FALSE)
res$pass <- c(res$observed[1:2] <= 1e-9, abs(res$observed[3:4]) <= 0.5, res$observed[5] <= 1e-9)
res$verdict <- if (all(res$pass)) "ACCEPTED" else "NOT ACCEPTED"
cat(sprintf("flag 0 median C_expected_sum %.1f (fit %.1f); flag 1 %.1f (what-if %.1f), %.2f%%\n", stats::median(x0$C_expected_sum), stats::median(E0$C_expected_sum),
            stats::median(x1$C_expected_sum), w$expected_catch_median, 100 * (stats::median(x1$C_expected_sum) / stats::median(x0$C_expected_sum) - 1)))
print(res, row.names = FALSE); cat("\nB86 switch:", res$verdict[1], "\n")
utils::write.csv(res, OUT, row.names = FALSE); cat("Wrote", OUT, "\n")
