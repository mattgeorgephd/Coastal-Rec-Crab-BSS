#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# run_osp_validation.R  -- batch the OSP-boat-count validation ladder.
#
# Loops the validation steps across BOTH models, one output folder per
# (model x step), using the same render mechanism as run_estimation.R. Set it
# going (`Rscript run_osp_validation.R`), come back, and push all of 05_output
# together. Nothing here changes run_config.R on disk; each job builds its own
# run_config in memory and injects it into the render environment.
#
# OUTPUT FOLDERS: the driver names its own folder
#   05_output/<YYYYMMDD>/<model-prefix>-<run_tag>
# and this script sets run_tag = "OSPval-<step>", so every run lands in e.g.
#   05_output/20260730/pooled-CPUE-OSPval-step2-osp
#   05_output/20260730/gear-type-CPUE-model-OSPval-step2-osp
# (date + model prefix + OSP-val step id, all three, automatically).
#
# CROSS-STEP SUMMARY: after each run it appends the headline numbers to
#   05_output/osp_validation_summary.csv
# so you (and the reviewer) can compare all steps at a glance.
#
# RESUMABLE: a job whose output folder already has port_total_Dungeness_Kept.csv
# is skipped, so if the batch is interrupted you can just re-run this script.
#
# TRIMMING: comment out any line in the JOBS list below to skip that run. Each
# render is a full model fit (multi-hour); the full matrix is ~14 runs.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(here); library(rmarkdown) })

# ---- environment: mirror run_estimation.R so renders behave identically ------
load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for (lib in install.lib) install.packages(lib, dependencies = TRUE)
invisible(sapply(load.lib, require, character.only = TRUE))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
invisible(purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source))

# ---- base run_config (single source of truth), then per-job overrides --------
source(here::here("run_config.R"))          # defines run_config
BASE <- run_config

model_rmd <- c(
  pooled        = here::here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd"),
  gear_resolved = here::here("01_BSS_models", "BSS-GH-gear-type-CPUE-model.Rmd")
)
prefix <- c(pooled = "pooled-CPUE-", gear_resolved = "gear-type-CPUE-model-")

# Every job fully specifies the four toggles, so the committed run_config.R
# values do not matter. off = the behavior-neutral baseline.
off <- list(use_osp_boat_counts = FALSE, use_crab_fraction = FALSE,
            crab_fraction_strata = "none", crab_fraction_fixed = NA,
            osp_scale_is_tau = FALSE)
job <- function(step, model, ...) list(step = step, model = model,
                                        ov = modifyList(off, list(...)))

# ============================ JOB MATRIX ====================================
# Ordered by priority: if you stop early or a run fails, the most informative
# runs are already done. Comment out lines to trim.
JOBS <- list(
  # -- Step 1b: gear-resolved behavior-neutral baseline (pooled step 1 done) --
  job("step1-baseline", "gear_resolved"),

  # -- Step 2: OSP second effort stream on (f off) --
  job("step2-osp", "pooled",        use_osp_boat_counts = TRUE),
  job("step2-osp", "gear_resolved", use_osp_boat_counts = TRUE),

  # -- Step 3: crabbing fraction f = 0.3 scalar (OSP off, isolate f) --
  job("step3-f030", "pooled",        use_crab_fraction = TRUE),
  job("step3-f030", "gear_resolved", use_crab_fraction = TRUE),

  # -- Step 3b: OSP + f together (tau off) = the production candidate --
  job("step3b-ospf", "pooled",        use_osp_boat_counts = TRUE, use_crab_fraction = TRUE),
  job("step3b-ospf", "gear_resolved", use_osp_boat_counts = TRUE, use_crab_fraction = TRUE),

  # -- Step 4: f sensitivity sweep, pinned (pooled only; mechanism is model-agnostic) --
  job("step4-fix020", "pooled", use_crab_fraction = TRUE, crab_fraction_fixed = 0.2),
  job("step4-fix050", "pooled", use_crab_fraction = TRUE, crab_fraction_fixed = 0.5),
  job("step4-fix100", "pooled", use_crab_fraction = TRUE, crab_fraction_fixed = 1.0),

  # -- Step 5: time-varying f smoke, strata = month, no pilot data (both models) --
  job("step5-fmonth", "pooled",        use_crab_fraction = TRUE, crab_fraction_strata = "month"),
  job("step5-fmonth", "gear_resolved", use_crab_fraction = TRUE, crab_fraction_strata = "month"),

  # -- Step 6: OSP-informs-tau experiment (OSP + f + tau on) --
  job("step6-osptau", "pooled",        use_osp_boat_counts = TRUE, use_crab_fraction = TRUE, osp_scale_is_tau = TRUE),
  job("step6-osptau", "gear_resolved", use_osp_boat_counts = TRUE, use_crab_fraction = TRUE, osp_scale_is_tau = TRUE)
)
# ============================================================================

banner <- function(msg) cat("\n", strrep("=", 74), "\n ", msg, "\n", strrep("=", 74), "\n", sep = "")
`%||%` <- function(a, b) if (is.null(a)) b else a

# Locate the folder a job wrote (search all dates for <prefix><run_tag>).
find_outdir <- function(model, run_tag) {
  hits <- list.dirs(here::here("05_output"), recursive = TRUE)
  hits <- hits[basename(hits) == paste0(prefix[[model]], run_tag)]
  if (length(hits)) hits[order(file.mtime(hits), decreasing = TRUE)][1] else NA_character_
}

num_from <- function(df, pick) tryCatch(suppressWarnings(as.numeric(pick(df))), error = function(e) NA_real_)

# Pull the headline numbers from a completed run's output folder.
extract_summary <- function(step, model, run_tag, outdir, minutes, ok) {
  rd <- function(f) tryCatch(read.csv(file.path(outdir, f), stringsAsFactors = FALSE,
                                       check.names = FALSE), error = function(e) NULL)
  row <- data.frame(step = step, model = model, run_tag = run_tag, ok = ok,
                    minutes = minutes, outdir = if (is.na(outdir)) "" else basename(outdir),
                    port_BSS_catch = NA_real_, port_PE_catch = NA_real_,
                    boat_BSS_catch = NA_real_, boat_BSS_effort = NA_real_, boat_PE_catch = NA_real_,
                    boat_div_frac = NA_real_, boat_pass = NA, f_crab_mean = NA_real_,
                    kappa_OSP_mean = NA_real_, stringsAsFactors = FALSE)
  if (is.na(outdir) || !dir.exists(outdir)) return(row)

  pt <- rd("port_total_Dungeness_Kept.csv")
  if (!is.null(pt)) {
    row$port_BSS_catch <- num_from(pt, function(d) d$BSS_median[d$Estimate == "Expected_Catch"])
    row$port_PE_catch  <- num_from(pt, function(d) d$PE[d$Estimate == "Expected_Catch"])
  }
  pv <- rd("pe_vs_bss_comparison.csv")
  if (!is.null(pv)) {
    b <- pv[grepl("private_boat \\(All gear\\)", pv$component), , drop = FALSE]
    if (nrow(b)) { row$boat_BSS_catch <- suppressWarnings(as.numeric(b$BSS_catch[1]))
                   row$boat_BSS_effort <- suppressWarnings(as.numeric(b$BSS_effort[1]))
                   row$boat_PE_catch  <- suppressWarnings(as.numeric(b$PE_catch[1])) }
  }
  cr <- rd("convergence_report.csv")
  if (!is.null(cr)) {
    b <- cr[grepl("private_boat_all_gear", cr$fit), , drop = FALSE]
    if (nrow(b)) { row$boat_div_frac <- suppressWarnings(as.numeric(b$divergence_fraction[1]))
                   row$boat_pass <- b$pass_convergence[1] }
  }
  fs <- rd("bss_full_summary_private_boat_all_gear_Dungeness_Kept.csv")
  if (!is.null(fs) && ncol(fs) >= 2) {
    v <- fs[[1]]; mcol <- if ("mean" %in% names(fs)) "mean" else names(fs)[2]
    gm <- function(pat) { i <- grepl(pat, v); if (any(i)) mean(suppressWarnings(as.numeric(fs[[mcol]][i])), na.rm = TRUE) else NA_real_ }
    row$f_crab_mean    <- gm("^f_crab(_out)?\\[|^f_crab(_out)?$")
    row$kappa_OSP_mean <- gm("^kappa_OSP")
  }
  row
}

summary_path <- here::here("05_output", "osp_validation_summary.csv")
append_summary <- function(row) {
  hdr <- !file.exists(summary_path)
  write.table(row, summary_path, sep = ",", row.names = FALSE, col.names = hdr,
              append = !hdr, qmethod = "double")
}

# ---- run the matrix ----------------------------------------------------------
banner(sprintf("OSP VALIDATION BATCH  |  %d jobs  |  start %s",
               length(JOBS), format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

for (j in JOBS) {
  run_tag <- paste0("OSPval-", j$step)
  tag_folder <- paste0(prefix[[j$model]], run_tag)

  # resume: skip if already completed
  done_dir <- find_outdir(j$model, run_tag)
  if (!is.na(done_dir) && file.exists(file.path(done_dir, "port_total_Dungeness_Kept.csv"))) {
    banner(sprintf("SKIP (already done): %s / %s", j$model, run_tag))
    append_summary(extract_summary(j$step, j$model, run_tag, done_dir, NA_real_, TRUE))
    next
  }

  rc <- modifyList(BASE, j$ov)
  rc$model   <- j$model
  rc$run_tag <- run_tag
  rc$run_weather <- FALSE

  banner(sprintf("RUN %s / %s   [OSP=%s f=%s strata=%s fixed=%s tau=%s]",
                 j$model, tag_folder,
                 rc$use_osp_boat_counts, rc$use_crab_fraction, rc$crab_fraction_strata,
                 rc$crab_fraction_fixed, rc$osp_scale_is_tau))

  run_env <- new.env(parent = globalenv())
  run_env$run_config <- rc
  t0 <- Sys.time()
  ok <- tryCatch({
    html <- rmarkdown::render(model_rmd[[j$model]], envir = run_env, quiet = FALSE)
    # relocate the HTML next to the CSVs (mirrors run_estimation.R)
    od <- tryCatch(get("output_dir", envir = run_env, inherits = FALSE), error = function(e) NA_character_)
    if (!is.na(od) && dir.exists(od) && normalizePath(dirname(html)) != normalizePath(od)) {
      if (file.copy(html, file.path(od, basename(html)), overwrite = TRUE)) file.remove(html)
    }
    TRUE
  }, error = function(e) { message("*** FAILED: ", j$model, " / ", run_tag, " : ", conditionMessage(e)); FALSE })
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  outdir <- find_outdir(j$model, run_tag)
  append_summary(extract_summary(j$step, j$model, run_tag, outdir, mins, ok))
  banner(sprintf("%s / %s  ->  %s  (%s min)", j$model, run_tag,
                 if (ok) "OK" else "FAILED", mins))
}

banner(sprintf("BATCH DONE %s   summary: %s", format(Sys.time(), "%H:%M:%S"), summary_path))
cat("\nReview 05_output/osp_validation_summary.csv, then:\n",
    "  git add 05_output run_osp_validation.R && git commit -m 'OSP validation runs' && git push\n")
