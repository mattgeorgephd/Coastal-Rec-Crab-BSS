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
# run_estimation.R  --  top-level orchestrator for the crab creel estimation.
#
# Renders the selected BSS model (pooled or gear-resolved) as a parameterized
# report, then, if requested, the weather-tide covariate module. All run-level
# settings come from run_config.R (edit that file, not this one).
#
# Run it either way:
#     source("run_estimation.R")                 # RStudio: Source (not Knit)
#     Rscript run_estimation.R                    # terminal / unattended
#     Rscript run_estimation.R --model gear_resolved
#     Rscript run_estimation.R --weather          # force weather on
#     Rscript run_estimation.R --no-weather       # force weather off
# CLI flags override the selection in run_config.R.
#
# Design notes (see the model documentation for the full rationale):
#   * The models stay as .Rmd reports; this script renders them via
#     rmarkdown::render(), so the full HTML diagnostic reports are preserved.
#   * Weather hand-off is Option A: the model and the weather module render into
#     the SAME environment (run_env), so the objects the model builds (dwg,
#     ie_data, L_eff_model, prep_days_crab, ...) satisfy the weather module's
#     `if(!exists("dwg"))` guard with no disk hand-off. Consequence: the pooled
#     model's memory is not freed before the weather module runs, and you cannot
#     re-run the weather module without re-running the model. If that becomes a
#     constraint, switch to the disk-bundle hand-off (documented as Phase 2).
###############################################################################

suppressPackageStartupMessages({
  library(here)
  library(rmarkdown)
})

load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)
rstan_options(auto_write = TRUE)
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)


# ---- 1. Load run configuration ------------------------------------------------
source(here::here("run_config.R"))     # defines: model, run_weather, run_config

# ---- 2. CLI overrides (optional) ----------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args)) {
  if ("--model" %in% .args) {
    .i <- which(.args == "--model")
    if (.i < length(.args)) model <- .args[.i + 1]
  }
  if ("--weather"    %in% .args) run_weather <- TRUE
  if ("--no-weather" %in% .args) run_weather <- FALSE
}

# ---- 3. Validate --------------------------------------------------------------
if (!model %in% c("pooled", "gear_resolved")) {
  stop("model must be 'pooled' or 'gear_resolved' (got '", model, "').",
       call. = FALSE)
}
if (isTRUE(run_weather) && model != "pooled") {
  stop("run_weather = TRUE is only valid with model = 'pooled'. The weather ",
       "module reuses the pooled run's objects. Set run_weather <- FALSE, or ",
       "model <- 'pooled'.", call. = FALSE)
}

model_rmd <- switch(model,
  pooled        = here::here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = here::here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd")
)
weather_rmd <- here::here("06_diagnostics",
                          "BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd")

stopifnot(file.exists(model_rmd))
if (isTRUE(run_weather)) stopifnot(file.exists(weather_rmd))

# ---- 4. Shared render environment (Option A hand-off) -------------------------
run_env <- new.env(parent = globalenv())
run_env$run_config <- run_config

# ---- 5. Helpers ---------------------------------------------------------------
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

banner <- function(msg) {
  cat("\n", strrep("=", 74), "\n ", msg, "\n", strrep("=", 74), "\n", sep = "")
}

# Render one report into run_env; co-locate its HTML with the CSVs the driver
# wrote (the driver sets its own `output_dir` inside run_env). Returns a small
# result list; the post-render file copy is protected so it can never turn a
# successful render into a reported failure.
render_stage <- function(rmd, label) {
  banner(sprintf("%s  |  %s  |  start %s",
                 label, basename(rmd), format(Sys.time(), "%H:%M:%S")))
  t0   <- Sys.time()
  html <- rmarkdown::render(rmd, envir = run_env, quiet = FALSE)
  outdir <- if (exists("output_dir", envir = run_env, inherits = FALSE)) {
    get("output_dir", envir = run_env, inherits = FALSE)
  } else {
    dirname(html)
  }
  # MOVE (not copy) the rendered HTML into the dated run folder, so no stale copy
  # is left in 01_BSS_models/. On any error the original render is kept as a
  # fallback, so relocation can never turn a successful render into a failure.
  final_html <- html
  tryCatch({
    if (dir.exists(outdir) &&
        normalizePath(dirname(html)) != normalizePath(outdir)) {
      dest <- file.path(outdir, basename(html))
      if (isTRUE(file.copy(html, dest, overwrite = TRUE))) {
        suppressWarnings(file.remove(html))
        final_html <- dest
      }
    }
  }, error = function(e) message("  (note: could not relocate HTML: ",
                                 conditionMessage(e), ")"))
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  banner(sprintf("%s DONE in %s min  ->  %s", label, mins, outdir))
  list(html = final_html, outdir = outdir, minutes = mins)
}

write_manifest <- function(stages, base_dir) {
  path <- file.path(base_dir, sprintf("run_manifest_%s.txt", run_stamp))
  con  <- file(path, "w")
  on.exit(close(con), add = TRUE)
  git_sha <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                      error = function(e) NA_character_)
  writeLines(c(
    "Run manifest",
    "============",
    paste("timestamp   :", run_stamp),
    paste("model       :", model),
    paste("run_weather :", run_weather),
    paste("git sha     :", if (length(git_sha)) git_sha else NA),
    "",
    "Stages:"), con)
  for (nm in names(stages)) {
    s <- stages[[nm]]
    if (is.null(s)) {
      writeLines(sprintf("  %-8s FAILED", nm), con)
    } else {
      writeLines(sprintf("  %-8s %6s min   %s", nm, s$minutes, s$outdir), con)
    }
  }
  writeLines(c("", "run_config (run-level overrides applied to the model):"), con)
  utils::capture.output(utils::str(run_config), file = con)
  writeLines(c("", "sessionInfo():"), con)
  utils::capture.output(print(utils::sessionInfo()), file = con)
  path
}

# ---- 6. Run -------------------------------------------------------------------
banner(sprintf("CRAB CREEL ESTIMATION  |  model = %s  |  weather = %s",
               model, run_weather))

stages <- list()

stages$model <- tryCatch(
  render_stage(model_rmd, sprintf("MODEL [%s]", model)),
  error = function(e) {
    message("\n*** MODEL render FAILED: ", conditionMessage(e), " ***")
    NULL
  })

if (is.null(stages$model)) {
  stop("Model stage failed; weather stage skipped. See console output above.",
       call. = FALSE)
}

if (isTRUE(run_weather)) {
  stages$weather <- tryCatch(
    render_stage(weather_rmd, "WEATHER covariates"),
    error = function(e) {
      message("\n*** WEATHER render FAILED (model outputs are intact): ",
              conditionMessage(e), " ***")
      NULL
    })
}

# ---- 7. Manifest --------------------------------------------------------------
# Base the manifest location on the last stage's output folder's PARENT, i.e.
# 05_output/<run_date>/, so it sits alongside the per-model subfolders.
last_outdir <- if (exists("output_dir", envir = run_env, inherits = FALSE)) {
  get("output_dir", envir = run_env, inherits = FALSE)
} else {
  here::here("05_output")
}
manifest_path <- tryCatch(write_manifest(stages, dirname(last_outdir)),
                          error = function(e) {
                            message("  (note: manifest not written: ",
                                    conditionMessage(e), ")"); NA_character_
                          })

banner(sprintf("ALL DONE  |  manifest: %s",
               if (is.na(manifest_path)) "(not written)" else manifest_path))
