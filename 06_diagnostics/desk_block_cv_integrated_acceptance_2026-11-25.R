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
# desk_block_cv_integrated_acceptance_2026-11-25.R
#
# DESK, NO MCMC: the acceptance of D54's repair as LIBRARY code (CHANGE_REGISTER B85).
#
# WHAT IS BEING ACCEPTED. write_block_cv_diagnostics() (03_R_functions/bss_block_cv.R), which
# every render calls once per fit, now scores the gear and trailer streams of a fit that carries
# the per-day effort effect with each sampled day's deviation integrated out, and writes the
# conditional tables beside them (loo_blockcond_*) with blockcv_checks_<fit>.csv. The method was
# prototyped as a desk script (desk_advisory_block_cv_2026-11-19.R, B80's R3a) and accepted
# against an exact refit of two held-out weeks (B83). This script checks that the library code
# IS that method, on the fits it was validated on, and runs it for the first time on the shore.
#
# THE INPUTS. The six B80 arms' run folders (05_output/20261119/pooled-CPUE-AR-<arm>, local, not
# committed) for each fit's Stan data and dates, and the fit cache (BSS_FIT_CACHE_DIR) for the
# stanfits; a fit is matched to its run folder by its pointwise log predictive density, as the
# prototype matched it.
#
# THE RULE, WRITTEN AND COMMITTED BEFORE THE LIBRARY CODE WAS RUN ON ANY FIT.
#   (a) On each of the six arms' boat all-gear fits, the library's loo_block_{trailer,osp,joint}
#       tables reproduce the prototype's committed loo_blockint_* tables
#       (05_output/advisory_replication_2026-11-19_arms/<arm>/): the same weeks, n_obs and
#       reliable flags, and elpd_block and pareto_k within 1e-9.
#   (b) On the same fits the library's loo_blockcond_* tables reproduce the drivers' own
#       loo_block_* tables in the run folders within 1e-9 (the conditional path is unchanged).
#   (c) On every fit read that carries the effect, both checks the library makes hold: (i) the
#       conditional log-likelihood rebuilt from mu0 and the fitted deviations equals the fit's
#       log_lik within 1e-6, and (ii) halving the grid moves no reliable week by 0.01 nats.
#   (d) The saved-draws path: save_bss_ppc_draws() now keeps sigma_uE and uE_raw, and from them
#       bss_effort_loglik() rebuilds the fit's own gear or trailer log_lik within 1e-9 (before
#       B85 it could not, D54's second consequence).
#   (e) On a fit WITHOUT the effect (the 2024-25 pot-closure fits, which A40's floor keeps
#       off), the library writes the conditional table only, equal to the driver's within 1e-9,
#       and no loo_blockcond_ or blockcv_checks_ file.
#   ACCEPTED if (a) to (e) hold on every fit they apply to. Not accepted: the library change
#   is not shipped, and the difference is found first.
#   THE SHORE, REPORTED AND NOT JUDGED, because no exact refit has been run on a shore fit: the
#   reliable weeks conditional and integrated on the shore all-gear and pot-closure fits of both
#   seasons. EXPECTED: the shore all-gear fit (weekly AR) recovers to about its pre-A39 state, 10
#   to 15 of 38 weeks on 2024-25, and no further, because under a weekly AR a held-out week still
#   removes its period's only anchor (campaign 1y.1); the 2025-26 shore pot-closure fit, 0 of 11
#   conditional, recovers to about the 2024-25 pot-closure fits' 8 or 9 of 11.
#
# WRITES 05_output/block_cv_integrated_acceptance_2026-11-25.csv
#
# HOW TO RUN (from the repository root)
#   BSS_FIT_CACHE_DIR=<dir> RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_block_cv_integrated_acceptance_2026-11-25.R
# About an hour on one core (the shore all-gear fits dominate).
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

ARMS <- c("S26-const", "S26-off", "S26-split", "S25-const", "S25-split", "S25-off")
SHORE_FROM <- c("S26-const", "S25-const")       # every arm of a season holds the same shore fits
FITS_SHORE <- c("shore_all_gear_Dungeness_Kept", "shore_ring_net_only_Dungeness_Kept")
FIT_BOAT <- "private_boat_all_gear_Dungeness_Kept"
EVIDENCE <- "05_output/advisory_replication_2026-11-19_arms"
OUT <- "05_output/block_cv_integrated_acceptance_2026-11-25.csv"
TOL <- 1e-9
CACHE <- Sys.getenv("BSS_FIT_CACHE_DIR", unset = "")
if (!nzchar(CACHE) || !dir.exists(CACHE)) stop("BSS_FIT_CACHE_DIR is not set or is not a directory")
rd <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL
run_dir <- function(sid) { d <- file.path("05_output", "20261119", paste0("pooled-CPUE-AR-", sid)); if (dir.exists(d)) d else NA_character_ }

rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# the dates of the sampled days, from the run folder's pointwise files (all the block CV reads)
days_from_pointwise <- function(dir, fit, D) {
  ev <- rep(as.Date(NA), D)
  for (st in c("gear", "trailer", "osp")) { p <- rd(file.path(dir, sprintf("loo_pointwise_%s_%s.csv", st, fit)))
    if (!is.null(p)) ev[p$day_index] <- as.Date(p$event_date) }
  data.frame(event_date = ev)
}
# the cached stanfit whose pointwise lpd is the run folder's
find_fit <- function(dir, fit, stream) {
  pw <- rd(file.path(dir, sprintf("loo_pointwise_%s_%s.csv", stream, fit))); if (is.null(pw)) return(NULL)
  for (f in Sys.glob(file.path(CACHE, sprintf("%s_att*_*.rds", fit)))) {
    x <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(x)) next
    ll <- tryCatch(as.matrix(rstan::extract(x, pars = paste0("log_lik_", stream))[[1]]), error = function(e) NULL)
    if (!is.null(ll) && ncol(ll) == nrow(pw) && max(abs(apply(ll, 2, .bcv_lmeanexp) - pw$lpd)) < 5e-4) return(list(fit = x, file = basename(f)))
    rm(x); gc(verbose = FALSE)
  }
  NULL
}
cmp_tab <- function(a, b) {     # same weeks, n_obs, reliable; max |diff| of elpd and k
  if (is.null(a) || is.null(b)) return(list(ok = FALSE, gap = NA_real_, why = "a table is missing"))
  m <- merge(a[, c("block", "n_obs", "elpd_block", "pareto_k", "reliable")], b[, c("block", "n_obs", "elpd_block", "pareto_k", "reliable")], by = "block")
  same <- nrow(m) == nrow(a) && nrow(m) == nrow(b) && all(m$n_obs.x == m$n_obs.y) && all(as.logical(m$reliable.x) == as.logical(m$reliable.y))
  gap <- max(abs(m$elpd_block.x - m$elpd_block.y), abs(m$pareto_k.x - m$pareto_k.y), na.rm = TRUE)
  list(ok = same && is.finite(gap) && gap <= TOL, gap = gap, why = if (same) "" else "weeks, n_obs or reliable flags differ")
}
# (d): the saved-draws path on this fit
saved_draws_check <- function(fit, sd_, label, stream) {
  td <- tempfile(); dir.create(td)
  f <- save_bss_ppc_draws(fit, sd_, label, td, quiet = TRUE)
  obj <- readRDS(f); unlink(td, recursive = TRUE)
  ll <- bss_effort_loglik(obj$draws, sd_, stream)$ll
  ref <- as.matrix(rstan::extract(fit, pars = paste0("log_lik_", stream))[[1]])
  c(has = all(c("sigma_uE", "uE_raw") %in% obj$pars), gap = max(abs(ll - ref)))
}

one_fit <- function(sid, fit_name, stream, boat) {
  dir <- run_dir(sid); if (is.na(dir)) { cat(sid, ": no run folder\n"); return(invisible()) }
  if (!file.exists(file.path(dir, sprintf("ppc_draws_%s.rds", fit_name)))) { cat(sprintf("%s %s: not fitted in this run (no saved draws)\n", sid, fit_name)); return(invisible()) }
  sd_ <- readRDS(file.path(dir, sprintf("ppc_draws_%s.rds", fit_name)))$stan_data
  hit <- find_fit(dir, fit_name, stream)
  if (is.null(hit)) { cat(sprintf("%s %s: no cached fit reproduces the run folder\n", sid, fit_name)); add(arm = sid, fit = fit_name, clause = "-", what = "fit found", observed = "no cached fit", pass = NA); return(invisible()) }
  has_ue <- identical(as.integer(sd_$effort_day_effect %||% 0L), 1L)
  td <- tempfile(); dir.create(td)
  t0 <- Sys.time()
  write_block_cv_diagnostics(hit$fit, sd_, days_from_pointwise(dir, fit_name, as.integer(sd_$D)), fit_name, td)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  streams <- sub(sprintf("^loo_block_(.*)_%s\\.csv$", fit_name), "\\1", list.files(td, pattern = sprintf("^loo_block_.*_%s\\.csv$", fit_name)))
  cat(sprintf("%s %s (%s): effect %s; streams %s; %.1f min\n", sid, fit_name, hit$file, if (has_ue) "on" else "off", paste(streams, collapse = "+"), mins))
  top <- if ("joint" %in% streams) "joint" else stream
  if (has_ue) {
    ck <- rd(file.path(td, sprintf("blockcv_checks_%s.csv", fit_name)))
    add(arm = sid, fit = fit_name, clause = "c", what = "checks (i) and (ii)",
        observed = if (is.null(ck)) "no checks file" else sprintf("(i) gap %.2e; (ii) grid gap %.2e; reliable %d conditional, %d integrated of %d", max(ck$conditional_loglik_max_abs_gap), ck$grid_halving_max_abs_elpd_gap[1], ck$reliable_conditional[1], ck$reliable_integrated[1], ck$weeks[1]),
        pass = !is.null(ck) && all(ck$check_i) && all(ck$check_ii))
    if (boat) for (s in c("trailer", "osp", "joint")) {
      r <- cmp_tab(rd(file.path(td, sprintf("loo_block_%s_%s.csv", s, fit_name))), rd(file.path(EVIDENCE, sid, sprintf("loo_blockint_%s_%s.csv", s, fit_name))))
      add(arm = sid, fit = fit_name, clause = "a", what = sprintf("%s table against the prototype's", s), observed = sprintf("max gap %.2e %s", r$gap, r$why), pass = r$ok)
    }
    for (s in streams) {
      r <- cmp_tab(rd(file.path(td, sprintf("loo_blockcond_%s_%s.csv", s, fit_name))), rd(file.path(dir, sprintf("loo_block_%s_%s.csv", s, fit_name))))
      add(arm = sid, fit = fit_name, clause = "b", what = sprintf("%s conditional table against the driver's", s), observed = sprintf("max gap %.2e %s", r$gap, r$why), pass = r$ok)
    }
    g <- saved_draws_check(hit$fit, sd_, fit_name, stream)
    add(arm = sid, fit = fit_name, clause = "d", what = "saved draws carry the effect and rebuild log_lik", observed = sprintf("sigma_uE and uE_raw saved: %s; max gap %.2e", as.logical(g[["has"]]), g[["gap"]]),
        pass = isTRUE(as.logical(g[["has"]])) && g[["gap"]] <= TOL)
    if (!boat) { ck <- rd(file.path(td, sprintf("blockcv_checks_%s.csv", fit_name)))
      add(arm = sid, fit = fit_name, clause = "shore", what = "reliable weeks, conditional and integrated (reported, not judged)",
          observed = if (is.null(ck)) "no checks" else sprintf("%d conditional, %d integrated of %d; median k %.2f and %.2f; sigma_uE %.2f", ck$reliable_conditional[1], ck$reliable_integrated[1], ck$weeks[1], ck$median_k_conditional[1], ck$median_k_integrated[1], ck$sigma_uE_mean[1]),
          pass = NA) }
  } else {
    extra <- list.files(td, pattern = "^(loo_blockcond_|blockcv_checks_)")
    for (s in streams) {
      r <- cmp_tab(rd(file.path(td, sprintf("loo_block_%s_%s.csv", s, fit_name))), rd(file.path(dir, sprintf("loo_block_%s_%s.csv", s, fit_name))))
      add(arm = sid, fit = fit_name, clause = "e", what = sprintf("%s table, fit without the effect, against the driver's", s), observed = sprintf("max gap %.2e %s; extra files %d", r$gap, r$why, length(extra)), pass = r$ok && !length(extra))
    }
    if (!boat) add(arm = sid, fit = fit_name, clause = "shore", what = "reliable weeks (the fit carries no effect)",
                   observed = { t <- rd(file.path(td, sprintf("loo_block_%s_%s.csv", top, fit_name))); sprintf("%d of %d", sum(t$reliable), nrow(t)) }, pass = NA)
  }
  unlink(td, recursive = TRUE); rm(hit); gc(verbose = FALSE)
}

for (sid in ARMS) one_fit(sid, FIT_BOAT, "trailer", boat = TRUE)
for (sid in SHORE_FROM) for (fn in FITS_SHORE) one_fit(sid, fn, "gear", boat = FALSE)
for (sid in SHORE_FROM) one_fit(sid, "private_boat_ring_net_only_Dungeness_Kept", "trailer", boat = TRUE)

res <- do.call(rbind, rows)
judged <- res[!is.na(res$pass), ]
verdict <- if (nrow(judged) && all(judged$pass) && all(c("a", "b", "c", "d") %in% judged$clause)) "ACCEPTED" else "NOT ACCEPTED"
res$verdict <- verdict
op <- options(width = 250); print(res, row.names = FALSE, right = FALSE); options(op)
cat("\nB85 rule:", verdict, "\n")
utils::write.csv(res, OUT, row.names = FALSE)
cat("Wrote", OUT, "\n")
