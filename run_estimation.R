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
# report. All run-level settings come from run_config.R (edit that file, not
# this one).
#
# Run it either way:
#     source("run_estimation.R")                 # RStudio: Source (not Knit)
#     Rscript run_estimation.R                    # terminal / unattended
#     Rscript run_estimation.R --model gear_resolved
# The --model flag overrides the selection in run_config.R.
#
# Design notes (see the method documentation for the full rationale):
#   * The models stay as .Rmd reports; this script renders them via
#     rmarkdown::render(), so the full HTML diagnostic reports are preserved.
#   * Config reaches the driver through a shared environment, not rmarkdown
#     `params:`. This script builds run_env, sets run_env$run_config, and renders
#     into it; the driver then does params <- modifyList(run_config, params_model).
#
# REMOVED 2026-09-13: the weather-tide covariate module and its --weather /
#   --no-weather flags. The FWC creel team advised against using weather
#   covariates and weather on its own was not helpful; the module's own committed
#   conclusion was already EXCLUSION, and its Stan fork had drifted about 40 data
#   variables behind the production model. This orchestrator is now single-path.
#   The finding is kept at 07_documentation/WEATHER_COVARIATE_ANALYSIS.md and the
#   module's method document at 07_documentation/archive/. CHANGE_REGISTER A29.
###############################################################################

# 2026-09-28 (B46): one package list and loader (03_R_functions/bss_packages.R): renv.lock's
# versions when renv is active, a stop naming the package when one cannot be had.
# D37 (2026-09-29): the guard now runs FIRST. A library(here) above it used to fail on a
# missing package with R's generic error before this message could name the fix; here and
# rmarkdown are attached by bss_load_packages() with every other package.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is missing: run renv::restore() first.", call. = FALSE)
source(here::here("03_R_functions", "bss_packages.R")); bss_load_packages()
rstan_options(auto_write = TRUE)
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)


# ---- 1. Load run configuration ------------------------------------------------
source(here::here("run_config.R"))     # defines: model, run_config

# ---- 2. CLI overrides (optional) ----------------------------------------------
# 2026-09-28 (B46): parsed strictly. `--model X` and `--model=X` both work; X is pooled,
# gear_resolved or both. A flag this script does not know, or --model with no value, STOPS
# the run: they used to be ignored silently and the run went ahead as pooled.
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args)) {
  if (any(c("--weather", "--no-weather") %in% .args))
    stop("--weather / --no-weather were removed with the weather-tide module on 2026-09-13. ",
         "See 07_documentation/WEATHER_COVARIATE_ANALYSIS.md for why covariates are excluded.",
         call. = FALSE)
  .i <- 1L
  while (.i <= length(.args)) {
    a <- .args[.i]
    if (grepl("^--model=", a)) {
      model <- sub("^--model=", "", a)
    } else if (identical(a, "--model")) {
      if (.i == length(.args) || grepl("^--", .args[.i + 1]))
        stop("--model needs a value: pooled, gear_resolved or both.", call. = FALSE)
      model <- .args[.i + 1]; .i <- .i + 1L
    } else {
      stop(sprintf("Unknown argument '%s'. The only flag is --model (pooled | gear_resolved | both).", a), call. = FALSE)
    }
    .i <- .i + 1L
  }
}

# ---- 3. Validate --------------------------------------------------------------
# "both" renders the pooled model, then the gear-resolved cross-check, and compares their
# port totals (the gear track's purpose). Each renders in its OWN environment.
if (length(model) != 1 || !model %in% c("pooled", "gear_resolved", "both")) {
  stop("model must be 'pooled', 'gear_resolved' or 'both' (got '", paste(model, collapse = ","), "').",
       call. = FALSE)
}
models <- if (identical(model, "both")) c("pooled", "gear_resolved") else model
model_rmds <- c(
  pooled        = here::here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = here::here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd"))[models]
stopifnot(all(file.exists(model_rmds)))

# ---- 4. Render environments ----------------------------------------------------
# One environment per model, created in section 6; each gets run_config and writes its own
# output_dir there.

# ---- 5. Helpers ---------------------------------------------------------------
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

banner <- function(msg) {
  cat("\n", strrep("=", 74), "\n ", msg, "\n", strrep("=", 74), "\n", sep = "")
}

# Render one report into run_env; co-locate its HTML with the CSVs the driver
# wrote (the driver sets its own `output_dir` inside run_env). Returns a small
# result list; the post-render file copy is protected so it can never turn a
# successful render into a reported failure.
render_stage <- function(rmd, label, run_env) {
  banner(sprintf("%s  |  %s  |  start %s",
                 label, basename(rmd), format(Sys.time(), "%H:%M:%S")))
  t0   <- Sys.time()
  # a blank run_tag falls back to this timestamp (a stale one from an earlier standalone
  # knit in the same session must not name this run's folder)
  options(crab_run_tag = format(Sys.time(), "%H%M%S"))
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
  # system2(stderr = FALSE): outside a repository git's "fatal: not a git repository" would
  # otherwise reach the console at the end of a successful render and read as a failure.
  git_sha <- tryCatch(suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
                      error = function(e) NA_character_)
  # 2026-09-28 (B43): the sha names the COMMIT; whether the render ran on exactly that commit
  # is a separate fact, so the tracked files that differed from it are listed. (Untracked files
  # are not: run outputs and delivery folders would swamp it.) Read at manifest time, after the
  # render; an edit made during a 3-4 h render would show here too, which is the conservative side.
  git_dirty <- tryCatch(suppressWarnings(system2("git", c("status", "--porcelain", "--untracked-files=no"), stdout = TRUE, stderr = FALSE)),
                        error = function(e) structure(NA_character_, status = -1L))
  # system2(stdout = TRUE) does not error when git itself fails (not a repository, a refused
  # safe.directory): it returns what git printed, often nothing, with a non-zero "status"
  # attribute. Nothing printed must not read as "clean", so the status and the sha are both checked.
  git_ok <- is.null(attr(git_dirty, "status")) && !(length(git_dirty) == 1 && is.na(git_dirty)) &&
    length(git_sha) == 1 && !is.na(git_sha) && is.null(attr(git_sha, "status")) && grepl("^[0-9a-f]{7,}$", git_sha)
  tree <- if (!git_ok) "unknown (git did not answer for this folder; the sha above may not identify the code)" else
    if (!length(git_dirty)) "clean (the render ran on the committed tree)" else
      sprintf("%d tracked file(s) differ from %s: %s%s", length(git_dirty), if (length(git_sha)) git_sha else "HEAD",
              paste(utils::head(trimws(substring(git_dirty, 4)), 20), collapse = ", "),
              if (length(git_dirty) > 20) ", ..." else "")
  writeLines(c(
    "Run manifest",
    "============",
    paste("timestamp   :", run_stamp),
    paste("model       :", paste(get0("models", ifnotfound = model), collapse = " + ")),
    paste("git sha     :", if (length(git_sha)) git_sha else NA),
    paste("git tree    :", tree),
    "",
    "Stages:"), con)
  for (nm in names(stages)) {
    s <- stages[[nm]]
    if (isTRUE(s$failed)) {
      writeLines(sprintf("  %-14s FAILED after %s min: %s%s", nm, s$minutes %||% "?", s$error,
                         if (!is.na(s$outdir %||% NA)) paste0("   (partial folder: ", s$outdir, ")") else ""), con)
    } else {
      writeLines(sprintf("  %-14s %6s min   %s", nm, s$minutes, s$outdir), con)
    }
  }
  # 2026-09-28 (B43): EVERY key. str() stops at 99 list elements by default, and run_config has
  # carried more than 99 keys since 2026-08-25 (109 then, 177 now). No run_estimation.R render
  # happened between that date and 2026-09-27, so the first render of the method of record wrote
  # the first truncated manifest: 99 keys and "[list output truncated]", the 78 missing ones
  # including all of section 2.10, the term the run was made to confirm. The driver's own
  # run_parameters.txt had the same defect until 2026-09-04 and was fixed there, not here; it is
  # complete in every run folder since, including this one.
  writeLines(c("", sprintf("run_config (run-level overrides applied to the model; all %d keys):", length(run_config))), con)
  utils::capture.output(utils::str(run_config, list.len = length(run_config) + 1L, vec.len = 8), file = con)
  writeLines(c("", "sessionInfo():"), con)
  utils::capture.output(print(utils::sessionInfo()), file = con)
  path
}

# ---- 6. Run -------------------------------------------------------------------
banner(sprintf("CRAB CREEL ESTIMATION  |  model = %s", paste(models, collapse = " + ")))

# 2026-09-28 (B46): a failed render is RECORDED, not fatal here. The manifest (section 7) is
# written either way, so a multi-hour run that dies late still leaves its configuration, commit
# and sessionInfo beside its partial folder; the script then exits with an error (section 8).
# It used to stop() before the manifest, so the FAILED branch could never run.
stages <- list()
for (m in models) {
  run_env <- new.env(parent = globalenv())
  run_env$run_config <- run_config
  t0 <- Sys.time()
  stages[[m]] <- tryCatch(
    render_stage(model_rmds[[m]], sprintf("MODEL [%s]", m), run_env),
    error = function(e) {
      message("\n*** MODEL render FAILED [", m, "]: ", conditionMessage(e), " ***")
      list(failed = TRUE, error = conditionMessage(e),
           minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
           outdir = if (exists("output_dir", envir = run_env, inherits = FALSE))
                      get("output_dir", envir = run_env, inherits = FALSE) else NA_character_)
    })
}

# ---- 7. Manifest and the cross-check -----------------------------------------------
# The manifest sits in 05_output/<run_date>/, beside the per-model folders.
.outdirs <- vapply(stages, function(s) s$outdir %||% NA_character_, character(1))
.base <- if (any(!is.na(.outdirs))) dirname(.outdirs[!is.na(.outdirs)][1]) else
  here::here("05_output", format(Sys.Date(), "%Y%m%d"))
dir.create(.base, recursive = TRUE, showWarnings = FALSE)
manifest_path <- tryCatch(write_manifest(stages, .base),
                          error = function(e) {
                            message("  (note: manifest not written: ",
                                    conditionMessage(e), ")"); NA_character_
                          })

# Both tracks rendered: compare their port totals, like for like (both write Expected_Catch
# since B44). The pre-set criterion is agreement within run_config$cross_check_tolerance (2%).
cross_check_path <- NA_character_
if (length(models) == 2 && !any(vapply(stages, function(s) isTRUE(s$failed), logical(1)))) {
  cross_check_path <- tryCatch({
    rd <- function(m) {
      pt <- utils::read.csv(file.path(stages[[m]]$outdir, "port_total_Dungeness_Kept.csv"), stringsAsFactors = FALSE)
      pt[pt$Estimate == "Expected_Catch", c("BSS_median", "BSS_lo95", "BSS_hi95")]
    }
    p <- rd("pooled"); g <- rd("gear_resolved")
    tol <- as.numeric(run_config$cross_check_tolerance %||% 0.02)
    rel <- (g$BSS_median - p$BSS_median) / p$BSS_median
    cc <- data.frame(pooled_median = p$BSS_median, pooled_lo95 = p$BSS_lo95, pooled_hi95 = p$BSS_hi95,
                     gear_median = g$BSS_median, gear_lo95 = g$BSS_lo95, gear_hi95 = g$BSS_hi95,
                     gear_minus_pooled_pct = round(100 * rel, 2), tolerance_pct = 100 * tol,
                     verdict = if (abs(rel) <= tol) "PASS" else "REVIEW",
                     pooled_folder = basename(stages$pooled$outdir), gear_folder = basename(stages$gear_resolved$outdir))
    f <- file.path(.base, sprintf("cross_check_%s.csv", run_stamp))
    utils::write.csv(cc, f, row.names = FALSE)
    banner(sprintf("CROSS-CHECK  |  pooled %s  gear %s  (%+.2f%%, criterion %.0f%%): %s",
                   format(p$BSS_median, big.mark = ","), format(g$BSS_median, big.mark = ","),
                   100 * rel, 100 * tol, cc$verdict))
    f
  }, error = function(e) { message("  (note: cross-check not written: ", conditionMessage(e), ")"); NA_character_ })
}

banner(sprintf("ALL DONE  |  manifest: %s%s",
               if (is.na(manifest_path)) "(not written)" else manifest_path,
               if (is.na(cross_check_path)) "" else paste0("  |  cross-check: ", cross_check_path)))

# ---- 8. Exit status ---------------------------------------------------------------
.failed <- names(stages)[vapply(stages, function(s) isTRUE(s$failed), logical(1))]
if (length(.failed))
  stop(sprintf("Model stage(s) failed: %s. The manifest above records the error and the partial folder.",
               paste(.failed, collapse = ", ")), call. = FALSE)
