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
# bss_run_warnings.R  (2026-09-29; D36 and the report overhaul)
#
# ONE LOG OF THE CONDITIONS A READER OF THE REPORT MUST SEE. A run used to announce a
# fallback, a skipped diagnostic or a clipped frame with one cat() line in a chunk whose
# output is hidden (results = 'hide'), or with a warning() that knitr collects at the end of
# a chunk and an html_document may not show. The report then read as if nothing had
# happened. bss_warn() records the condition here AND raises it as an R warning; the report's
# closing section prints the whole log, and the drivers write it to run_warnings.csv.
#
#   bss_warn("turnover", "tau_shore prior: derivation unavailable ...")
#   bss_warn("OSP", "overlap diagnostic skipped: ...", severity = "note")
#   bss_warn_log()      # data.frame: time, source, severity, message
#
# State lives in a module-local environment (the bss_timers.R pattern), so no function
# depends on a driver global, and it resets each time 03_R_functions is sourced, which is
# once per render.
###############################################################################

.bss_run_warnings <- new.env(parent = emptyenv())
.bss_run_warnings$log <- data.frame(time = character(0), source = character(0),
                                    severity = character(0), message = character(0),
                                    stringsAsFactors = FALSE)

# severity: "warning" (the estimate may be affected; also raised as an R warning) or
# "note" (recorded for the reader, not raised).
bss_warn <- function(source, message, severity = c("warning", "note")) {
  severity <- match.arg(severity)
  msg <- paste(message, collapse = " ")
  .bss_run_warnings$log <- rbind(.bss_run_warnings$log,
    data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), source = as.character(source),
               severity = severity, message = msg, stringsAsFactors = FALSE))
  if (identical(severity, "warning")) warning(sprintf("[%s] %s", source, msg), call. = FALSE)
  invisible(msg)
}

bss_warn_log <- function() .bss_run_warnings$log

bss_warn_reset <- function() {
  .bss_run_warnings$log <- .bss_run_warnings$log[0, , drop = FALSE]
  invisible(NULL)
}
