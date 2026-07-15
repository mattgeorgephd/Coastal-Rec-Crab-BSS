#!/usr/bin/env Rscript
###############################################################################
# run_rg_sweep.R  --  T1.3 R_G prior-sensitivity sweep (pooled model).
#
# Renders three pooled runs back-to-back, one per R_G_prior_mu value, each into
# its OWN dated output folder so they do not overwrite each other:
#
#     05_output/<date>/pooled-CPUE-run5-RG-1.00
#     05_output/<date>/pooled-CPUE-run5-RG-1.28   (~ the empirical value)
#     05_output/<date>/pooled-CPUE-run5-RG-1.50
#
# The R_G prior is data-driven by default; this sweep OVERRIDES it via
# run_config$R_G_prior_mu (already wired in 03_R_functions/prep_bss_crab_pooled.R).
# A tighter R_G_prior_sigma makes the prior bind harder; 0.3 matches production.
#
# Run it like run_estimation.R (Source, not Knit):
#     source("run_rg_sweep.R")
#
# When it finishes, compare the port total across the three folders against the
# Run-1 baseline (83,035): read each pooled-CPUE-run5-RG-*/port_total_Dungeness_Kept.csv.
# This is a robustness study, not a correctness fix; if the port total is stable
# across the three priors, the boat does not rest heavily on the R_G prior.
###############################################################################

suppressPackageStartupMessages({
  library(here)
  library(rmarkdown)
})

load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here","readxl")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for (lib in install.lib) install.packages(lib, dependencies = TRUE)
invisible(sapply(load.lib, require, character.only = TRUE))
rstan_options(auto_write = TRUE)
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)

# Base configuration (defines run_config); the sweep overrides three fields per run.
source(here::here("run_config.R"))

model_rmd <- here::here("01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd")
stopifnot(file.exists(model_rmd))

# ---- Sweep grid ---------------------------------------------------------------
rg_grid  <- c(1.0, 1.28, 1.5)   # R_G_prior_mu values (1.28 ~ the empirical value)
rg_sigma <- 0.3                 # prior SD; tighter binds harder (0.3 = production)

banner <- function(msg) cat("\n", strrep("=", 74), "\n ", msg,
                            "\n", strrep("=", 74), "\n", sep = "")

results <- list()
for (rg in rg_grid) {
  tag <- sprintf("run5-RG-%s", formatC(rg, format = "f", digits = 2))  # run5-RG-1.00
  cfg <- run_config
  cfg$R_G_prior_mu    <- rg
  cfg$R_G_prior_sigma <- rg_sigma
  cfg$run_tag         <- tag

  banner(sprintf("R_G SWEEP  |  R_G_prior_mu = %.2f  |  tag = %s  |  start %s",
                 rg, tag, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()

  run_env <- new.env(parent = globalenv())
  run_env$run_config <- cfg

  html <- rmarkdown::render(model_rmd, envir = run_env, quiet = FALSE)

  # Relocate the rendered HTML into the driver's dated output folder (the driver
  # set output_dir inside run_env from cfg$run_tag), matching run_estimation.R.
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
  banner(sprintf("R_G_prior_mu = %.2f DONE in %s min  ->  %s", rg, mins, outdir))
  results[[tag]] <- list(R_G_prior_mu = rg, outdir = outdir, minutes = mins)
}

banner("R_G SWEEP COMPLETE")
cat("Compare the port total across the three runs against the Run-1 baseline (83,035):\n")
for (tag in names(results)) {
  r <- results[[tag]]
  cat(sprintf("  R_G_prior_mu = %.2f  ->  %s\n",
              r$R_G_prior_mu, file.path(r$outdir, "port_total_Dungeness_Kept.csv")))
}
