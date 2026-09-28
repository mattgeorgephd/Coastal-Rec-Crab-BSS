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
# bss_day_covariate_report.R  (2026-09-28; CHANGE_REGISTER B43)
#
# THE EFFORT DAY COVARIATES AS FITTED, per fit, for the report.
#
# WHY THIS EXISTS. Since 2026-09-27 the method of record carries one day covariate: the
# NWS Small-Craft-Advisory flag on the private-boat ALL-GEAR effort process (A30), confined
# to that fit by marine_hazard_terms_for() (B42). The first render of the method of record
# (05_output/20260927/pooled-CPUE-canonical-2024-25) showed, in its HTML, the covariate's
# pre-fit SCREEN (a quasi-Poisson GLM on the sampled days, section 3.7 of the pooled report)
# and the configured scope ("applied to the sub-season gear regime(s): all_gear"), and never
# the FITTED coefficient or which fits actually carried the column: the preps print that per
# fit, but the fit loop's chunk hides its output, so it reached the console only. A reviewer
# reading the report could not see the one parameter the adoption added. This table is that
# record: every fit, the K_open columns it carries (or none), how many of its days each column
# flags, and the posterior of B_open_out with its rate ratio.
#
# CONTENTS
#   bss_day_covariate_labels(bss_data)            the K_open column labels of one fit's Stan
#                                                 data, under either prep's convention
#   bss_day_covariate_table(bss_all, params, ...) one row per fit x column (a "(none)" row for
#                                                 a fit that carries no column)
#
# Pure functions; nothing runs at source time. The draws are read through `extract_b_open`
# (default rstan::extract on B_open_out), so the harness can test the table without rstan.
###############################################################################

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# The pooled prep carries the labels as attr(bss_data, "opener_labels"); the gear prep as a
# comma-joined string in bss_data$.opener_labels. Either way the result is a character vector
# in K_open column order (character(0) when the fit carries no column).
bss_day_covariate_labels <- function(bss_data) {
  lab <- attr(bss_data, "opener_labels")
  if (is.null(lab) || !length(lab)) {
    s <- bss_data$.opener_labels %||% ""
    lab <- strsplit(paste(s, collapse = ","), ",", fixed = TRUE)[[1]]
  }
  lab <- as.character(lab)
  lab[!is.na(lab) & nzchar(lab)]
}

bss_day_covariate_table <- function(bss_all, params = list(), extract_b_open = NULL,
                                    digits = 3) {
  if (is.null(extract_b_open))
    extract_b_open <- function(fit) rstan::extract(fit, pars = "B_open_out")$B_open_out
  sel <- params$marine_hazard_selected %||% list()
  row <- function(fit, pop, ss, covariate, source, days_flagged, column_sum, D,
                  q = c(NA_real_, NA_real_, NA_real_)) {
    data.frame(fit = fit, population = pop %||% NA_character_, subseason = ss %||% NA_character_,
               covariate = covariate, source = source,
               days_flagged = days_flagged, column_sum = column_sum, fit_days = D,
               coef_median = round(q[1], digits), coef_lo95 = round(q[2], digits), coef_hi95 = round(q[3], digits),
               rate_ratio = if (all(is.finite(q))) sprintf("%.2f (%.2f-%.2f)", exp(q[1]), exp(q[2]), exp(q[3])) else NA_character_,
               stringsAsFactors = FALSE)
  }
  out <- list()
  for (label in names(bss_all)) {
    b <- bss_all[[label]]
    if (isTRUE(b$pe_fallback) || is.null(b$bss_data) || is.null(b$fit)) next
    bd <- b$bss_data
    lab <- bss_day_covariate_labels(bd)
    K <- as.integer(bd$K_open %||% 0L); D <- as.integer(bd$D %||% NA_integer_)
    if (K == 0L) {
      out[[label]] <- row(label, b$population, b$subseason, "(none)", NA_character_, NA_integer_, NA_real_, D)
      next
    }
    if (length(lab) != K)   # includes a fit that carries columns and has lost its labels: refused, not reported as "(none)"
      stop(sprintf("bss_day_covariate_table(): %s carries K_open = %d columns but %d labels (%s)",
                   label, K, length(lab), paste(lab, collapse = ", ")), call. = FALSE)
    X <- matrix(as.numeric(bd$X_open_flat), nrow = D)
    dr <- extract_b_open(b$fit)
    dr <- if (is.null(dim(dr))) matrix(dr, ncol = K) else matrix(dr, nrow = dim(dr)[1])
    if (ncol(dr) != K)
      stop(sprintf("bss_day_covariate_table(): %s: B_open_out has %d columns, K_open is %d", label, ncol(dr), K), call. = FALSE)
    marine <- sel[[b$population %||% ""]] %||% character(0)
    for (k in seq_len(K)) {
      q <- as.numeric(stats::quantile(dr[, k], c(0.5, 0.025, 0.975), names = FALSE))
      out[[paste(label, k)]] <- row(label, b$population, b$subseason, lab[k],
                                    if (lab[k] %in% marine) "marine hazard" else "other-fishery opener",
                                    sum(X[, k] == 1), round(sum(X[, k]), 1), D, q)
    }
  }
  if (!length(out))
    return(data.frame(fit = character(0), population = character(0), subseason = character(0),
                      covariate = character(0), source = character(0), days_flagged = integer(0),
                      column_sum = numeric(0), fit_days = integer(0), coef_median = numeric(0),
                      coef_lo95 = numeric(0), coef_hi95 = numeric(0), rate_ratio = character(0),
                      stringsAsFactors = FALSE))
  do.call(rbind, unname(out))
}
