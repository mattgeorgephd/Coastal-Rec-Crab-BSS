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
# bss_output_dir.R  (2026-09-28, B46; both drivers)
#
# The run's output folder: 05_output/<YYYYMMDD>/<prefix><run_tag>. Three rules, each from a
# defect the fixed tag "canonical-2024-25" exposed:
#   1. A folder that already holds files is never written into again. Two same-day renders
#      under one tag used to share a folder, and the report's per-fit readers (every
#      <prefix>_*.csv, every plot_bss_*.png in it) then showed the earlier run's rows beside
#      the new one's whenever the sub-season set changed. The second run gets "-HHMMSS"
#      appended and says so. (A standalone RStudio knit is exempt: its knit hook creates the
#      empty folder first and passes the tag by option.)
#   2. A tag that names a season (yyyy-yy) the run does not model stops the run. A 2025-26
#      render under the shipped "canonical-2024-25" would have filed itself as 2024-25.
#   3. The tag is made path-safe: anything outside [A-Za-z0-9._-] becomes "-" (05_output
#      holds historical folders named with spaces and "=").
###############################################################################
bss_output_dir <- function(prefix, run_config = NULL, base = here::here("05_output"),
                           date = Sys.Date(), now = Sys.time()) {
  cfg_tag  <- if (is.list(run_config)) run_config$run_tag else NULL
  from_cfg <- !is.null(cfg_tag) && length(cfg_tag) == 1 && nzchar(cfg_tag)
  tag <- if (from_cfg) as.character(cfg_tag) else getOption("crab_run_tag", format(now, "%H%M%S"))
  clean <- gsub("[^A-Za-z0-9._-]+", "-", tag)
  if (!identical(clean, tag)) message(sprintf("run_tag '%s' is written as '%s' (path-safe characters only).", tag, clean))
  if (from_cfg) {
    named <- unique(regmatches(tag, gregexpr("[0-9]{4}-[0-9]{2}", tag))[[1]])
    sf <- sort(as.character(run_config$season_filter))
    # a span may be labelled by its first season's start year and last season's end year
    # ("two-season-2024-26" for c("2024-25", "2025-26"))
    span <- if (length(sf) > 1) paste0(substr(sf[1], 1, 4), "-", substr(sf[length(sf)], 6, 7)) else character(0)
    bad <- setdiff(named, c(sf, span))
    if (length(bad))
      stop(sprintf(paste0("run_tag '%s' names season %s, but season_filter is %s. Change run_tag in ",
                          "run_config.R with the season window (NEW_SEASON_GUIDE.md), or the run is filed ",
                          "under a season it does not estimate."),
                   tag, paste(bad, collapse = " and "), paste(run_config$season_filter, collapse = " + ")), call. = FALSE)
  }
  dir <- file.path(base, format(as.Date(date), "%Y%m%d"), paste0(prefix, clean))
  if (from_cfg && dir.exists(dir) && length(list.files(dir, all.files = FALSE, no.. = TRUE))) {
    dir2 <- paste0(dir, "-", format(now, "%H%M%S"))
    message(sprintf("Output folder %s already holds a run; this one is written to %s.", basename(dir), basename(dir2)))
    dir <- dir2
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}
