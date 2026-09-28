#!/usr/bin/env Rscript
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
# run_rg_sweep.R  --  T1.3 R_G prior-sensitivity sweep (pooled model).
#
# Renders three pooled runs back-to-back, one per R_G_prior_mu value, each into its OWN
# output folder so they do not overwrite each other:
#
#     05_output/<date>/pooled-CPUE-RG-sweep-1.00
#     05_output/<date>/pooled-CPUE-RG-sweep-1.28   (~ the 2024-25 empirical value)
#     05_output/<date>/pooled-CPUE-RG-sweep-1.50
#
# The R_G prior is data-driven by default (the season's interview ratio; both tracks since
# B46); this sweep OVERRIDES it via run_config$R_G_prior_mu, which both preps read first.
# R_G_prior_sigma stays at 0.3, the production value.
#
# SHIPS DRY_RUN <- TRUE (2026-09-28, B49). As committed until then it had no dry-run mode
# (sourcing it started three multi-hour renders) and its grid held two rungs while this
# header described three (the 1.00 rung was dropped mid-batch on 2026-07-14). The dry run
# prints each rung's resolved configuration, its output folder and the keys it changes
# against run_config.R, and fits nothing. Set DRY_RUN <- FALSE to render.
#
# Run it like run_estimation.R (Source, not Knit):
#     source("06_diagnostics/run_rg_sweep.R")
#
# COMPARE A NEW SWEEP AGAINST THE BOX AT THE TOP OF
# 07_documentation/development_notes/PIPELINE_STATUS.md, never against a number written into
# a runner comment (this header named a "current production total" three times and was wrong
# each time the box moved). The script ends by printing each rung's Expected_Catch median and
# interval from its own port_total_Dungeness_Kept.csv and writes them to
# 05_output/rg_sweep_<date>_summary.csv. This is a robustness study, not a correctness fix: if
# the port total is stable across the three priors, the estimate does not rest on the R_G
# prior (the 2026-07-14/15 sweep found every component within about 0.5%).
###############################################################################

DRY_RUN <- TRUE

suppressPackageStartupMessages({
  library(here)
  library(rmarkdown)
})

# 2026-09-28 (B46): the shared loader (03_R_functions/bss_packages.R): renv.lock versions, a stop on a missing package.
source(here::here("03_R_functions", "bss_packages.R")); bss_load_packages()
rstan_options(auto_write = TRUE)
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)

# Base configuration (defines run_config); the sweep overrides three fields per run.
source(here::here("run_config.R"))

model_rmd <- here::here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd")
stopifnot(file.exists(model_rmd))

# ---- Sweep grid ---------------------------------------------------------------
rg_grid  <- c(1.00, 1.28, 1.50)   # R_G_prior_mu values (1.28 ~ the 2024-25 empirical value)
rg_sigma <- 0.3                   # prior SD; tighter binds harder (0.3 = production)

banner <- function(msg) cat("\n", strrep("=", 74), "\n ", msg,
                            "\n", strrep("=", 74), "\n", sep = "")
rung_cfg <- function(rg) {
  cfg <- run_config
  cfg$R_G_prior_mu    <- rg
  cfg$R_G_prior_sigma <- rg_sigma
  cfg$run_tag         <- sprintf("RG-sweep-%s", formatC(rg, format = "f", digits = 2))
  cfg
}

if (isTRUE(DRY_RUN)) {
  banner(sprintf("R_G SWEEP, DRY RUN: %d rung(s), nothing is fitted", length(rg_grid)))
  for (rg in rg_grid) {
    cfg <- rung_cfg(rg)
    changed <- names(cfg)[!vapply(names(cfg), function(k) identical(cfg[[k]], run_config[[k]]), logical(1))]
    cat(sprintf("  R_G_prior_mu = %.2f  sigma = %.2f  ->  05_output/%s/pooled-CPUE-%s\n      keys changed against run_config.R: %s\n",
                rg, rg_sigma, format(Sys.Date(), "%Y%m%d"), cfg$run_tag, paste(changed, collapse = ", ")))
  }
  cat("\n  The pooled and gear preps read R_G_prior_mu before the season's empirical ratio:",
      any(grepl("params$R_G_prior_mu %||%", readLines(here::here("03_R_functions", "prep_bss_crab_pooled.R")), fixed = TRUE)) &&
      any(grepl("params$R_G_prior_mu %||%", readLines(here::here("03_R_functions", "prep_bss_crab_gear.R")), fixed = TRUE)), "\n")
  cat("  Set DRY_RUN <- FALSE to render the three rungs (about one production run each).\n")
} else {
  results <- list()
  for (rg in rg_grid) {
    cfg <- rung_cfg(rg)
    banner(sprintf("R_G SWEEP  |  R_G_prior_mu = %.2f  |  tag = %s  |  start %s",
                   rg, cfg$run_tag, format(Sys.time(), "%H:%M:%S")))
    t0 <- Sys.time()
    run_env <- new.env(parent = globalenv())
    run_env$run_config <- cfg
    # one rung failing must not lose the others: record it and move on
    html <- tryCatch(rmarkdown::render(model_rmd, envir = run_env, quiet = FALSE),
                     error = function(e) { message("  RUNG FAILED: ", conditionMessage(e)); NULL })
    outdir <- if (exists("output_dir", envir = run_env, inherits = FALSE))
      get("output_dir", envir = run_env, inherits = FALSE) else NA_character_
    # Relocate the rendered HTML into the driver's output folder, matching run_estimation.R.
    if (!is.null(html) && !is.na(outdir)) tryCatch({
      if (dir.exists(outdir) && normalizePath(dirname(html)) != normalizePath(outdir)) {
        dest <- file.path(outdir, basename(html))
        if (isTRUE(file.copy(html, dest, overwrite = TRUE))) suppressWarnings(file.remove(html))
      }
    }, error = function(e) message("  (note: could not relocate HTML: ", conditionMessage(e), ")"))
    mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    banner(sprintf("R_G_prior_mu = %.2f %s in %s min  ->  %s", rg,
                   if (is.null(html)) "FAILED" else "DONE", mins, outdir))
    results[[cfg$run_tag]] <- list(R_G_prior_mu = rg, outdir = outdir, minutes = mins, failed = is.null(html))
  }

  banner("R_G SWEEP COMPLETE")
  rows <- lapply(names(results), function(tag) {
    r <- results[[tag]]
    pt <- if (!is.na(r$outdir)) file.path(r$outdir, "port_total_Dungeness_Kept.csv") else ""
    d <- if (nzchar(pt) && file.exists(pt)) utils::read.csv(pt, stringsAsFactors = FALSE) else NULL
    e <- if (!is.null(d)) d[d[[2]] == "Expected_Catch", , drop = FALSE] else NULL
    data.frame(R_G_prior_mu = r$R_G_prior_mu, R_G_prior_sigma = rg_sigma, failed = r$failed, minutes = r$minutes,
               port_median = if (NROW(e)) e$BSS_median else NA, port_lo95 = if (NROW(e)) e$BSS_lo95 else NA,
               port_hi95 = if (NROW(e)) e$BSS_hi95 else NA, folder = r$outdir)
  })
  summ <- do.call(rbind, rows)
  print(summ, row.names = FALSE)
  out <- here::here("05_output", sprintf("rg_sweep_%s_summary.csv", format(Sys.Date(), "%Y%m%d")))
  utils::write.csv(summ, out, row.names = FALSE)
  cat("Written:", out, "\nCompare against the box at the top of PIPELINE_STATUS.md, on the same configuration.\n")
  if (any(summ$failed)) stop("R_G sweep: ", sum(summ$failed), " rung(s) failed; see above.", call. = FALSE)
}
