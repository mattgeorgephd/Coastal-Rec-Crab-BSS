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
# run_block_cv_exact_refit_2026-11-19.R
#
# THE ACCEPTANCE TEST OF THE REPAIRED HELD-OUT-WEEK CHECK (CHANGE_REGISTER B83, for D54).
#
# WHY. B80's addendum R3a reads the fitted model on held-out weeks through an approximation:
# the posterior without a week is obtained by Pareto-smoothed importance sampling on the
# week's likelihood with each held-out day's deviation integrated out
# (desk_advisory_block_cv_2026-11-19.R). Its own diagnostics say it is reliable (Pareto k at
# or below 0.7 on 39 to 41 of 41 weeks). The direct test of an approximation to a refit is
# the refit.
#
# WHAT IT DOES. For one arm and a few weeks: remove the week's trailer counts from the arm's
# own Stan data (the list saved in ppc_draws_<fit>.rds, so the data are the driver's to the
# byte), re-index the day deviations, sample the boat all-gear fit again without them, and
# score the removed week under that posterior with each of its days' deviations integrated
# out by the same quadrature. That is the exact held-out log predictive density the
# approximation estimates. The catch and every other stream stay in, as they do in the
# approximation (it holds out the effort counts of the week, nothing else).
#
# THE RULE, WRITTEN BEFORE ANY REFIT WAS SAMPLED (R11).
#   Arm: S26-const, the shipped term on 2025-26 (one effort stream, so the test is clean).
#   Weeks, chosen by a rule and not by eye, from the arm's committed integrated table: the
#   week with the LARGEST Pareto k (the approximation's hardest week) and the week at the
#   MEDIAN k (a typical one).
#   Sampler: the pipeline's default budget (bss_iter_default / bss_warmup_default, 4 chains)
#   with the boat all-gear fit's adapt_delta and treedepth, bss_seed and bss_init_r. A refit
#   is read only if it has no divergent transition above 1% and R-hat below 1.01 on E_sum.
#   ACCEPTED if, on both weeks, the approximation is within 0.25 nats of the exact value.
#   Why 0.25: R3a sums weekly differences between two arms over about 40 weeks with a paired
#   SE of 6.7 nats; independent errors of 0.25 per week and arm would add about 2.2 nats of
#   noise to that sum, a third of its SE, and anything larger would begin to matter.
#   NOT ACCEPTED on either week: R3a's readings are reported with that warning, and D54's
#   repair goes back to an exact K-fold.
#   EXPECTED: accepted, the median-k week within 0.1 nats.
#
# HOW TO RUN. Ships DRY_RUN <- TRUE (prints the plan, samples nothing).
#   BSS_FIT_CACHE_DIR=<dir> Rscript 06_diagnostics/run_block_cv_exact_refit_2026-11-19.R --go
#   ... --stage=S26-off --weeks=2026-W19,2026-W05     another arm, named weeks
# Needs the arm's run folder (not committed). About ten minutes per week on four cores.
# WRITES 05_output/advisory_replication_2026-11-19_exact_refit.csv
###############################################################################
DRY_RUN <- TRUE                    # TRUE prints the plan and samples nothing
.args <- commandArgs(trailingOnly = TRUE)
if ("--go" %in% .args || identical(Sys.getenv("BLOCK_CV_EXACT_GO"), "1")) DRY_RUN <- FALSE
.arg <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), .args, value = TRUE)); if (length(v)) v[1] else default }
STAGE <- .arg("stage", "S26-const")
WEEKS <- { w <- .arg("weeks", ""); if (nzchar(w)) strsplit(w, ",")[[1]] else NULL }
TOL <- 0.25

.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
.here <- function(...) file.path(.root, ...)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
if (!DRY_RUN) { suppressPackageStartupMessages({ library(here) }); source(.here("03_R_functions", "bss_packages.R")); bss_load_packages(); rstan_options(auto_write = TRUE)
} else suppressWarnings(suppressPackageStartupMessages(try({ library(dplyr); library(tidyr); library(lubridate); library(readxl); library(here); library(purrr); library(stringr); library(tibble) }, silent = TRUE)))
invisible(lapply(list.files(.here("03_R_functions"), full.names = TRUE), function(f) source(f)))
source(.here("run_config.R"))

FIT_BOAT <- "private_boat_all_gear_Dungeness_Kept"
EVIDENCE <- .here("05_output", "advisory_replication_2026-11-19_arms")
OUT <- .here("05_output", "advisory_replication_2026-11-19_exact_refit.csv")
rd <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL
GRID_H <- 0.025; GRID_ZMAX <- 8
.row_lse <- function(A) { m <- A[cbind(seq_len(nrow(A)), max.col(A, ties.method = "first"))]; m + log(rowSums(exp(A - m))) }
.lmeanexp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
# the day's joint density of its counts with the day deviation integrated out, one value per draw (desk_advisory_block_cv's integral)
day_loglik <- function(mu0, y, sig, size) {
  z <- seq(-GRID_ZMAX, GRID_ZMAX, by = GRID_H); lw <- stats::dnorm(z, log = TRUE) + log(GRID_H); nd <- nrow(mu0)
  shift <- outer(sig, z) - 0.5 * sig^2; A <- matrix(0, nd, length(z))
  for (i in seq_along(y)) A <- A + stats::dnbinom(y[i], mu = exp(log(mu0[, i]) + shift), size = size, log = TRUE)
  .row_lse(A + rep(lw, each = nd))
}

int_tab <- rd(file.path(EVIDENCE, STAGE, sprintf("loo_blockint_trailer_%s.csv", FIT_BOAT)))
if (is.null(int_tab)) stop("the arm's integrated table is not there: run desk_advisory_block_cv_2026-11-19.R for ", STAGE, " first")
if (is.null(WEEKS)) { ok <- int_tab[int_tab$reliable %in% TRUE, ]; ok <- ok[order(ok$pareto_k), ]
  WEEKS <- c(ok$block[nrow(ok)], ok$block[ceiling(nrow(ok) / 2)]) }
cat(sprintf("EXACT REFIT OF HELD-OUT WEEKS (B83)  DRY_RUN = %s  arm %s  weeks %s\n", DRY_RUN, STAGE, paste(WEEKS, collapse = ", ")))
for (w in WEEKS) { r <- int_tab[int_tab$block == w, ]
  cat(sprintf("  %s: %d sampled days, Pareto k %.2f, approximated leave-out elpd %.3f\n", w, r$n_obs, r$pareto_k, r$elpd_block)) }
if (DRY_RUN) { cat("DRY RUN: nothing is sampled and nothing is written. --go samples.\n")
} else {
ds <- Sys.glob(.here("05_output", "*", paste0("pooled-CPUE-AR-", STAGE, "*")))
ds <- ds[file.exists(file.path(ds, sprintf("ppc_draws_%s.rds", FIT_BOAT)))]
if (!length(ds)) stop("the arm's run folder (with ppc_draws) is not on this machine")
dir <- ds[length(ds)]
sd0 <- readRDS(file.path(dir, sprintf("ppc_draws_%s.rds", FIT_BOAT)))$stan_data
pw <- rd(file.path(dir, sprintf("loo_pointwise_trailer_%s.csv", FIT_BOAT)))
Tn <- as.integer(sd0[["T_n"]]); stopifnot(nrow(pw) == Tn, identical(as.integer(sd0[["T_I"]])[seq_len(Tn)], as.integer(pw$observed)))
wk <- bss_block_weeks(pw$event_date)
P <- run_config
res <- list()
for (w in WEEKS) {
  out_i <- which(wk == w); keep <- setdiff(seq_len(Tn), out_i); stopifnot(length(out_i) > 0)
  sd1 <- sd0                                                      # the driver's data list, attributes and all, minus one week's trailer counts
  for (k in c("T_I", "day_T", "section_T", "Eday_T")) sd1[[k]] <- as.integer(sd0[[k]])[keep]
  old <- as.integer(sd0[["Eday_T"]])[keep]; map <- sort(unique(old)); sd1[["Eday_T"]] <- match(old, map)
  sd1[["T_n"]] <- length(keep); sd1[["n_Eday"]] <- length(map)
  for (k in c("T_I", "day_T", "section_T", "Eday_T")) dim(sd1[[k]]) <- length(keep)
  cat(sprintf("\n== %s: %d counts on %d days held out; %d counts on %d days kept ==\n", w, length(out_i), length(unique(sd0[["day_T"]][out_i])), length(keep), length(map)))
  t0 <- Sys.time()
  fit <- bss_stan_fit(file = .here("02_stan_models", "crab_bss_pooled.stan"), data = sd1, label = sprintf("blockcv_exact_%s_%s", STAGE, w),
                      chains = P$bss_chains %||% 4, iter = P$bss_iter_default %||% 3000, warmup = P$bss_warmup_default %||% 1500,
                      control = list(adapt_delta = P$bss_delta_boat_allgear %||% 0.99, max_treedepth = P$bss_treedepth_boat_allgear %||% 13),
                      seed = P$bss_seed, cores = P$bss_cores %||% 4, init_r = P$bss_init_r %||% 2,
                      pars = c("lambda_E_S", "R_G_boat", "r_E", "sigma_uE", "B_open", "E_sum", "C_expected_sum"), include = TRUE)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE); div <- mean(unlist(lapply(sp, function(x) x[, "divergent__"])))
  rhat <- tryCatch(rstan::summary(fit, pars = "E_sum")$summary[, "Rhat"], error = function(e) NA_real_)
  ex <- rstan::extract(fit, pars = c("lambda_E_S", "R_G_boat", "r_E", "sigma_uE")); nd <- length(ex$r_E); G <- dim(ex$lambda_E_S)[4]
  days <- unique(as.integer(sd0[["day_T"]])[out_i]); ll <- numeric(nd)
  for (dd in days) { ii <- out_i[as.integer(sd0[["day_T"]])[out_i] == dd]
    mu0 <- vapply(ii, function(i) ex$lambda_E_S[, as.integer(sd0[["section_T"]])[i], dd, G] / as.numeric(ex$R_G_boat), numeric(nd))
    ll <- ll + day_loglik(matrix(mu0, nrow = nd), as.integer(sd0[["T_I"]])[ii], as.numeric(ex$sigma_uE), as.numeric(ex$r_E)) }
  exact <- .lmeanexp(ll)
  pr <- exp(ll - max(ll)); mcse <- stats::sd(pr) / (mean(pr) * sqrt(nd))          # delta-method Monte Carlo SE of the exact value, draws taken as independent
  r <- int_tab[int_tab$block == w, ]; cnd <- rd(file.path(EVIDENCE, STAGE, sprintf("loo_block_trailer_%s.csv", FIT_BOAT))); rc <- cnd[cnd$block == w, ]
  ok_fit <- isTRUE(div <= 0.01) && isTRUE(rhat < 1.01)
  res[[w]] <- data.frame(stage = STAGE, week = w, days_held_out = length(days), counts_held_out = length(out_i), pareto_k_integrated = r$pareto_k,
    elpd_exact_refit = exact, exact_mc_se_naive = mcse, elpd_approximated = r$elpd_block, gap = r$elpd_block - exact, tolerance = TOL,
    refit_draws = nd, refit_divergent_fraction = div, refit_rhat_E_sum = unname(rhat), refit_minutes = mins, refit_read = ok_fit,
    accepted = ok_fit && abs(r$elpd_block - exact) <= TOL,
    pareto_k_conditional = if (nrow(rc)) rc$pareto_k else NA_real_, elpd_conditional_approximation = if (nrow(rc)) rc$elpd_block else NA_real_, stringsAsFactors = FALSE)
  cat(sprintf("  exact %.3f; approximated %.3f; gap %+.3f nats (tolerance %.2f); divergent %.2f%%, R-hat %.4f; %.1f min: %s\n", exact, r$elpd_block, r$elpd_block - exact, TOL,
              100 * div, rhat, mins, if (!ok_fit) "REFIT NOT READ" else if (abs(r$elpd_block - exact) <= TOL) "within tolerance" else "OUTSIDE tolerance"))
  rm(fit, ex); invisible(gc())
}
res <- do.call(rbind, res)
utils::write.csv(res, OUT, row.names = FALSE)
cat(sprintf("\nR11: %s. Wrote %s\n", if (all(res$accepted)) "ACCEPTED on every week refitted" else "NOT ACCEPTED on at least one week", basename(OUT)))
}
