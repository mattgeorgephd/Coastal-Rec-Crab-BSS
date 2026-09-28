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
# run_tau_sweep.R  --  tau_boat prior-sensitivity sweep (gear-resolved model).
#
# The EXACT (multi-refit) companion to the single-run projection diagnostic
# (03_R_functions/diagnose_tau_boat_sensitivity.R). Where that projection scales
# the boat catch analytically under the assumption that tau is prior-dominated,
# this script actually re-fits the model at each tau_boat_prior_mu, so it is the
# gold-standard check that remains valid even when boat I/E days pin tau.
#
# It renders one gear-resolved run per tau_boat_prior_mu value, each into its own
# dated output folder so they do not overwrite each other:
#
#     05_output/<date>/gear-type-CPUE-model-tau-0.90
#     05_output/<date>/gear-type-CPUE-model-tau-1.20   (the RETIRED production prior)
#     05_output/<date>/gear-type-CPUE-model-tau-1.50
#
# tau_boat_prior_mu is the boat deployment turnover (trips per present group per
# day); it enters as L_data for the boat via 03_R_functions/bss_effort_spec.R and
# is already overridable from run_config. A tighter tau_boat_prior_sigma makes the
# prior bind harder. 0.3 was production until 2026-09-08; production now ships
# tau_boat_prior_mu = "calibration" (resolved per window, ~2.98 on 2024-25) with
# tau_boat_prior_sigma = 0.5, so NEITHER this grid nor this sigma is production any more.
#
# Run it like run_estimation.R (Source, not Knit):
#     source("06_diagnostics/run_tau_sweep.R")
#
# When it finishes, compare the boat all-gear catch and the port total across the
# STALE BASELINE WARNING (2026-09-03). The numbers on this line (43,314 boat all-gear, 82,957
# port) are from 2026-07 and are badly wrong under Method v2.0; the "current" figures this
# line quoted in their place went stale in turn (71,513 until 2026-09-12, 94,376 until
# 2026-09-28), so it quotes none. Read the box at the top of PIPELINE_STATUS.md instead of
# any number written into a runner comment. MORE IMPORTANT: this sweep is superseded. diagnose_tau_boat
# _sensitivity (run_config.R diagnose_tau_sensitivity = TRUE, a 5-point grid) reproduces its
# refits to ~0.2% on every run, and production adopted shared_tau = TRUE on 2026-09-01, which
# replaces the per-day tau_boat draws this sweep perturbs with a single estimated tau_bar of
# 2.597 [2.064, 3.249] -- an interval that excludes this script's entire 0.9/1.2/1.5 grid.
# folders against the production run (boat ~43,314; port ~82,957): read each
# gear-type-CPUE-model-tau-*/port_total_Dungeness_Kept.csv and
# pe_vs_bss_comparison.csv. If the boat scales ~linearly with tau_boat, the boat
# is prior-dominated and the field priority is boat I/E (egress-classification)
# coverage; if it barely moves, tau is data-identified and the prior is not the
# lever. This is a robustness study, not a correctness fix.
###############################################################################

# --- SUPERSEDED (2026-09-28, B46). This runner refuses to FIT; it is kept to be READ. ------
# Its own header already said "this sweep is superseded": the grid (0.9 / 1.2 / 1.5, SD 0.3)
# predates the OSP/trailer turnover calibration (~3.0, SD 0.5) and the shared turnover, and it
# printed a July 2026 "production" total as its reference. It was nonetheless LIVE and unguarded
# with no DRY_RUN, so sourcing it started three gear-resolved fits. The in-driver projection
# (diagnose_tau_sensitivity = TRUE, gear report section 6b) is the current tool.
if (!exists("bss_superseded_runner")) {
  .sr <- here::here("03_R_functions", "bss_superseded_runner.R")
  if (!file.exists(.sr)) stop("The superseded-runner guard was not found at ", .sr, "; refusing to run.", call. = FALSE)
  source(.sr)
}
bss_superseded_runner(
    runner   = "06_diagnostics/run_tau_sweep.R",
    question = "how does the boat total respond to the tau_boat prior centre",
    settled_by = paste("the turnover is calibrated from the OSP/trailer overlap (tau_boat_prior_mu = \"calibration\", ~3.0)",
                       "and shared across the boat fits (shared_tau, A2); the single-run projection in the gear report,",
                       "section 6b (diagnose_tau_sensitivity), answers the sensitivity question without refits"),
    what_would_happen = paste("Three gear-resolved fits at prior centres 0.9 / 1.2 / 1.5 with SD 0.3, none of them the",
                              "production prior, with no DRY_RUN; it would compare them against a July 2026 total."),
    levers_pinned = "0")

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

model_rmd <- here::here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd")
stopifnot(file.exists(model_rmd))

# ---- Sweep grid ---------------------------------------------------------------
tau_grid  <- c(0.9, 1.2, 1.5)   # tau_boat_prior_mu values. NONE of these is production:
                                #   1.2 was the centre until 2026-09-08, and production now
                                #   resolves "calibration" to ~2.98. See the warning above.
tau_sigma <- 0.3                # prior SD; tighter binds harder. Production ships 0.5.

banner <- function(msg) cat("\n", strrep("=", 74), "\n ", msg,
                            "\n", strrep("=", 74), "\n", sep = "")

results <- list()
for (tv in tau_grid) {
  tag <- sprintf("tau-%s", formatC(tv, format = "f", digits = 2))  # tau-0.90
  cfg <- run_config
  cfg$tau_boat_prior_mu    <- tv
  cfg$tau_boat_prior_sigma <- tau_sigma
  cfg$run_tag              <- tag
  cfg$diagnose_tau_sensitivity <- FALSE  # the sweep IS the exact check; no projection needed

  banner(sprintf("TAU_BOAT SWEEP  |  tau_boat_prior_mu = %.2f  |  tag = %s  |  start %s",
                 tv, tag, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()

  run_env <- new.env(parent = globalenv())
  run_env$run_config <- cfg

  html <- rmarkdown::render(model_rmd, envir = run_env, quiet = FALSE)

  # Relocate the rendered HTML into the driver's dated output folder, matching
  # run_estimation.R and run_rg_sweep.R.
  outdir <- if (exists("output_dir", envir = run_env, inherits = FALSE)) {
    get("output_dir", envir = run_env, inherits = FALSE)
  } else dirname(html)
  tryCatch({
    if (dir.exists(outdir) &&
        normalizePath(dirname(html)) != normalizePath(outdir)) {
      dest <- file.path(outdir, basename(html))
      if (isTRUE(file.copy(html, dest, overwrite = TRUE))) {
        suppressWarnings(file.remove(html))
      }
    }
  }, error = function(e) message("  (note: could not relocate HTML: ",
                                 conditionMessage(e), ")"))

  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  banner(sprintf("tau_boat_prior_mu = %.2f DONE in %s min  ->  %s", tv, mins, outdir))
  results[[tag]] <- list(tau_boat_prior_mu = tv, outdir = outdir, minutes = mins)
}

banner("TAU_BOAT SWEEP COMPLETE")
cat("Compare the boat all-gear catch and port total across the runs",
    "(the July 2026 reference: boat ~43,314, port ~82,957; NOT the current method of record):\n")
for (tag in names(results)) {
  r <- results[[tag]]
  cat(sprintf("  tau_boat_prior_mu = %.2f  ->  %s\n",
              r$tau_boat_prior_mu, file.path(r$outdir, "port_total_Dungeness_Kept.csv")))
}
