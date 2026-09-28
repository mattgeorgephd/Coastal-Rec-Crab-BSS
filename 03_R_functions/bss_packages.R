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
# bss_packages.R  (2026-09-28, B46; run_estimation.R and both drivers)
#
# ONE list of the packages a run needs and ONE way of getting them. Replaces three copies of
#   install.packages(lib, dependencies = TRUE); sapply(load.lib, require, character = TRUE)
# which (a) installed whatever CRAN held that day rather than the versions renv.lock pins,
# (b) had no repository set, so under Rscript (repos "@CRAN@") it failed outright, and
# (c) only WARNED when a package was missing, so the run went on and failed hours later with
# an unrelated-looking error. gt and patchwork were loaded and never used; they are gone.
#
# With the project's renv active (the .Rprofile does it; renv.lock pins every version), a
# missing package is restored from the lockfile. Without renv, a missing package is installed
# from CRAN (cloud.r-project.org unless a repository is already set). Either way a package that
# is still missing STOPS the run here, naming it.
###############################################################################
bss_required_packages <- c("tidyverse", "lubridate", "rstan", "here", "readxl", "rmarkdown",
                           "knitr", "loo", "suncalc", "digest")
bss_attached_packages <- c("tidyverse", "lubridate", "rstan", "here", "readxl")

bss_load_packages <- function(pkgs = bss_required_packages, attach = bss_attached_packages) {
  have <- function(p) vapply(p, requireNamespace, logical(1), quietly = TRUE)
  missing <- pkgs[!have(pkgs)]
  if (length(missing)) {
    renv_on <- nzchar(Sys.getenv("RENV_PROJECT")) && requireNamespace("renv", quietly = TRUE)
    if (renv_on) {
      message("Restoring from renv.lock: ", paste(missing, collapse = ", "))
      renv::restore(packages = missing, prompt = FALSE)
    } else {
      repos <- getOption("repos")
      if (is.null(repos) || !length(repos) || any(repos == "@CRAN@")) repos <- c(CRAN = "https://cloud.r-project.org")
      message("Installing from CRAN (renv is not active, so versions are NOT the pinned ones): ",
              paste(missing, collapse = ", "))
      utils::install.packages(missing, repos = repos)
    }
    still <- missing[!have(missing)]
    if (length(still))
      stop("Required package(s) could not be installed: ", paste(still, collapse = ", "),
           ". Run renv::restore() in the project (README, 'Setting up R'), then retry.", call. = FALSE)
  }
  for (p in attach) suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(TRUE)
}
