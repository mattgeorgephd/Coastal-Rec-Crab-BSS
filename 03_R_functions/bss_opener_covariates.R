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
# bss_opener_covariates.R   (improvement 4, 2026-08-25; shared by both drivers)
#
# Turns the other-fishery opener calendar into EFFORT covariates for the BSS, with an
# automatic significance screen. Generalizes the razor-dig-only `B3` term into the
# K_open / X_open / B_open block now present in both Stan models.
#
# THE HONEST FRAMING, UP FRONT
#   "Include a covariate if it is significant" is model selection on the same data that
#   then fits the term. Two consequences the reader has to hold:
#
#   1. MULTIPLICITY. The spillover diagnostic runs eight EFFORT tests (four openers x
#      two populations). At a raw alpha of 0.05 about one false positive is expected by
#      chance alone. The default adjustment here is Benjamini-Hochberg over exactly those
#      eight, which is why opener_auto_p_adjust ships as "BH" and not "none".
#
#   2. POST-SELECTION BIAS. A coefficient selected because it was large is, conditional
#      on selection, biased away from zero, and its reported interval does not cover at
#      its nominal rate. So "auto" is a SCREEN, not a verdict. The arbiter is predictive:
#      refit with opener_covariate_mode = "off" and compare elpd_loo on the effort stream.
#
#   The project already has a worked example of why this matters. Razor dig -> shore
#   effort came in at an adjusted p of 0.045, was fitted as B3 in Run 3, and produced NO
#   predictive gain (elpd_loo within 1 SE of the B3-off run on every stream) while moving
#   the port total 0.6%. Under the default BH adjustment over the eight effort tests that
#   term does NOT survive the screen, while MA2 halibut -> boat trailers (p = 3.6e-06)
#   does. The default behaves correctly on the one case where the answer is known.
#
# THE BOAT / CRABBING-FRACTION INTERACTION (read before switching on a boat opener)
#   Boat effort is a count of ALL private boats, so an MA2-halibut term is a legitimate
#   description of that series: halibut days really do add about 32 trailers to a ramp
#   that averages 7. But those boats are not crabbing. A CONSTANT crabbing fraction f
#   converts the fitted surge into crab effort at the same rate as any other day, so
#   fitting the surge more faithfully while multiplying it by a flat f makes the boat
#   catch MORE biased on exactly those days, not less. Pair any boat opener covariate
#   with crab_fraction_strata = "opener" or "month_opener". This function emits a loud
#   warning when a boat opener is selected and f is not opener-aware.
#
# WHY EFFORT ONLY
#   All eight opener-vs-catch-rate tests in the 2024-25 diagnostic were null (p = 0.29 to
#   0.90). There is no evidence for a CPUE opener effect, so none is offered.
#
# CONTENTS
#   opener_flag_series(params)            per-date logical flags for every candidate
#   opener_select(spillover, params)      the screen: which openers, per population, why
#   opener_design_matrix(days, sel, pop)  the flattened D x K_open design + labels for one fit
###############################################################################

# Human-readable labels, shared with the spillover diagnostic's `Series` values.
.opener_labels <- c(
  ma2_salmon_open     = "MA2 salmon open",
  ma2_halibut_open    = "MA2 halibut open",
  ma2_bottomfish_open = "MA2 bottomfish open",
  razor_nearby_dig    = "Razor dig (nearby beaches)",
  razor_any_dig       = "Razor dig (any beach)"
)
.opener_series <- c(shore = "Shore gear (effort)", private_boat = "Boat trailers (effort)")


# Per-date logical flags for every opener the calendar carries. Returns a tibble keyed by
# event_date; a date absent from the sparse razor calendar means "no dig", not NA.
opener_flag_series <- function(params) {
  ev <- prep_fishery_events(params)
  dplyr::full_join(ev$ma2, ev$razor, by = "event_date") |>
    dplyr::mutate(dplyr::across(-event_date, ~ tidyr::replace_na(.x, FALSE))) |>
    dplyr::arrange(event_date)
}


# The screen. Returns a list with:
#   $shore / $private_boat  character vectors of selected opener flag names
#   $table                  one row per (population, candidate): raw p, adjusted p,
#                           estimate, whether it was selected, and why not if not
#   $mode                   the mode actually used
# `spillover` is the object returned by diagnose_fishery_spillover(); NULL is tolerated
# (auto mode then selects nothing and says so).
opener_select <- function(spillover, params) {
  mode <- tolower(params$opener_covariate_mode %||% "off")
  cand <- list(shore        = params$opener_candidates_shore %||% character(0),
               private_boat = params$opener_candidates_boat  %||% character(0))
  adj_method <- tolower(params$opener_auto_p_adjust %||% "BH")
  thresh     <- params$opener_auto_p %||% 0.05

  empty <- list(shore = character(0), private_boat = character(0),
                table = tibble::tibble(), mode = mode, note = character(0))

  if (identical(mode, "off")) {
    empty$note <- "opener_covariate_mode = 'off'; no opener covariates."
    return(empty)
  }

  if (identical(mode, "manual")) {
    # Intersect with the CANDIDATE list, not merely with the known flag names: the
    # report table below is built from the candidates, so a manually named non-candidate
    # would enter X_open and the fit with no row explaining it.
    .keep <- function(want, cand_p) {
      want <- intersect(want %||% character(0), names(.opener_labels))
      drop <- setdiff(want, cand_p)
      if (length(drop))
        warning("opener_manual_* names not in the candidate list, ignored: ",
                paste(drop, collapse = ", "), call. = FALSE)
      intersect(want, cand_p)
    }
    sel <- list(shore        = .keep(params$opener_manual_shore, cand$shore),
                private_boat = .keep(params$opener_manual_boat,  cand$private_boat))
    tbl <- dplyr::bind_rows(lapply(names(sel), function(p) {
      if (!length(cand[[p]])) return(NULL)
      tibble::tibble(population = p, opener = cand[[p]],
                     adj_estimate = NA_real_, p_raw = NA_real_, p_adj = NA_real_,
                     selected = cand[[p]] %in% sel[[p]], reason = "manual")
    }))
    return(list(shore = sel$shore, private_boat = sel$private_boat, table = tbl,
                mode = mode, note = "Manual opener selection; no significance screen applied."))
  }

  # --- auto -----------------------------------------------------------------
  if (is.null(spillover) || is.null(spillover$adjusted) || nrow(spillover$adjusted) == 0) {
    empty$note <- paste("opener_covariate_mode = 'auto' but the spillover diagnostic produced",
                        "no adjusted table (run_fishery_spillover_diag = FALSE?); no opener",
                        "covariates selected.")
    return(empty)
  }

  lab2flag <- stats::setNames(names(.opener_labels), unname(.opener_labels))
  eff <- spillover$adjusted |>
    dplyr::filter(.data$Series %in% unname(.opener_series)) |>
    dplyr::mutate(population = names(.opener_series)[match(.data$Series, .opener_series)],
                  opener     = unname(lab2flag[.data$Opener]))

  # Multiplicity adjustment over the EFFORT tests only. The CPUE tests are not part of
  # the selection family because CPUE is never offered a covariate, so including them
  # would make the adjustment conservative for no reason.
  p_adj_method <- switch(adj_method, bh = "BH", "bh" = "BH", bonferroni = "bonferroni", none = "none", adj_method)
  # n is pinned to the FULL family size. stats::p.adjust defaults n to length(p) AFTER
  # dropping NAs, so a term the spillover diagnostic could not identify (collinear with
  # day-type or month, which it warns is likely for salmon and bottomfish) would shrink
  # the family and make the screen quietly MORE permissive -- in exactly the case the
  # documentation uses as its worked example. The family is the effort tests offered,
  # identified or not.
  .n_family <- nrow(eff)
  eff$p_adj <- if (identical(tolower(p_adj_method), "none")) eff$adj_p
               else stats::p.adjust(eff$adj_p, method = p_adj_method, n = .n_family)

  tbl <- dplyr::bind_rows(lapply(names(cand), function(p) {
    if (!length(cand[[p]])) return(NULL)
    rows <- eff |> dplyr::filter(.data$population == p, .data$opener %in% cand[[p]])
    missing <- setdiff(cand[[p]], rows$opener)
    dplyr::bind_rows(
      rows |> dplyr::transmute(population = p, opener,
                               adj_estimate, p_raw = adj_p, p_adj,
                               selected = is.finite(p_adj) & p_adj < thresh,
                               reason = dplyr::case_when(
                                 !is.finite(p_adj)  ~ paste0("no adjusted p (", note, ")"),
                                 p_adj < thresh     ~ sprintf("adjusted p %.4g < %.3g", p_adj, thresh),
                                 TRUE               ~ sprintf("adjusted p %.4g >= %.3g", p_adj, thresh))),
      if (length(missing))
        tibble::tibble(population = p, opener = missing, adj_estimate = NA_real_,
                       p_raw = NA_real_, p_adj = NA_real_, selected = FALSE,
                       reason = "not tested by the spillover diagnostic")
    )
  }))

  sel <- lapply(names(cand), function(p) {
    if (is.null(tbl) || !nrow(tbl)) return(character(0))
    tbl$opener[tbl$population == p & tbl$selected]
  })
  names(sel) <- names(cand)

  note <- sprintf(paste("Auto opener screen: %d effort test(s) in the multiplicity family",
                        "(%d with an identifiable estimate), adjustment '%s', threshold %.3g on the",
                        "ADJUSTED p. Selection is a screen, not a verdict: confirm any selected term",
                        "against an opener-free run's effort elpd_loo."),
                  .n_family, sum(is.finite(eff$adj_p)), p_adj_method, thresh)

  list(shore = sel$shore, private_boat = sel$private_boat, table = tbl, mode = mode, note = note)
}


# Build the D x K_open design matrix for ONE fit. Drops any column that is not
# identifiable inside this fit's own window (fewer than params$opener_min_days days on
# either side of the indicator), because a near-constant column is a free parameter that
# does nothing but widen the level.
#
# `extra` lets the compatibility razor_dig_mode switch add its column without going
# through the screen; duplicates are removed, so razor can never be counted twice.
opener_design_matrix <- function(days, selected, flags, params, extra = character(0)) {
  cols <- unique(c(selected, extra))
  min_days <- params$opener_min_days %||% 10
  out <- list(); labels <- character(0); dropped <- character(0)

  if (length(cols) && !is.null(flags) && nrow(flags) > 0) {
    fl <- days |>
      dplyr::select(event_date) |>
      dplyr::left_join(flags, by = "event_date")
    for (cc in cols) {
      if (!cc %in% names(fl)) { dropped <- c(dropped, sprintf("%s (absent from calendar)", cc)); next }
      v <- as.numeric(tidyr::replace_na(fl[[cc]], FALSE))
      n_open <- sum(v == 1); n_shut <- sum(v == 0)
      if (n_open < min_days || n_shut < min_days) {
        dropped <- c(dropped, sprintf("%s (%d open / %d closed days < %d in this window)",
                                      cc, n_open, n_shut, min_days))
        next
      }
      out[[cc]] <- v; labels <- c(labels, cc)
    }
  }

  # Emitted FLAT (column-major, length D * K_open) to match the Stan data block: a
  # zero-length numeric is a marshalling path rstan already relies on throughout this
  # model, whereas a D x 0 matrix is a needless edge case. Stan rebuilds the matrix in
  # transformed data. No dimnames: the labels travel separately.
  X_flat <- if (length(out)) as.numeric(unlist(out, use.names = FALSE)) else numeric(0)

  list(K_open = length(labels), X_open = X_flat, labels = labels, dropped = dropped)
}
