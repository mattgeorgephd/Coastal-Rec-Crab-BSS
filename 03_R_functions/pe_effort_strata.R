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
# PE EFFORT STRATA: the fill and the variance, in ONE place  (2026-09-12)
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   run_pe_pooled() and run_pe_gear() each carried their own copy of the (week x
#   day-type) effort stratification, the empty-cell fill and the SE formula. Three
#   defects lived in both copies, and a fix applied to one would have been a fix
#   applied to half the pipeline:
#
#   1. THE SE IS COMPUTED OVER LESS THAN HALF THE EFFORT IT REPORTS. The SE was
#      sqrt(N^2 * sd^2 / max(n, 1)) with sd from the cell's own sampled days, and
#      sd() of ONE observation is NA, replaced by 0. On the rebuilt 2024-25 inputs
#      the shore all-gear component has 37 of 93 strata with exactly one sampled
#      day; they carry 82 of 289 calendar days and 14,816 of 27,345 effort units.
#      Together with the 21 empty cells, 54% of the reported effort contributed
#      ZERO variance, and the component's SE read 410 on 27,345 (1.5%) for a design
#      with ~50% day coverage and a within-cell CV of 0.2 to 0.7. The boat all-gear
#      figure was 49% of effort at zero variance, SE 5.5%.
#   2. AN IMPUTED CELL WAS FREE. Under pe_empty_effort_stratum = "day_type" or
#      "local_day_type" an empty cell gets a donor mean and sd_daily <- NA, so it
#      contributed its full point estimate and no variance at all. The measured
#      consequence: the shore all-gear effort SE was BIT-IDENTICAL (410) across
#      "zero" (27,345 effort), "day_type" (33,585) and "local_day_type" (35,292).
#      Imputing 29% more effort at zero cost is not a defensible artifact.
#   3. THE EFFORT FILL AND THE CPUE FILL WERE ON DIFFERENT SCALES. The effort fill
#      could be month-local while the empty-stratum CPUE fill (pe_empty_stratum)
#      was always the sub-season-wide ratio-of-sums. Crab CPUE has a strong
#      within-season gradient, so a February cell was getting February-rate effort
#      times a season-average CPUE.
#
# WHAT IT PROVIDES
#   pe_effort_donors()      the donor table: (month x day_type), then day_type, then
#                           the sub-season mean, each with n, mean and sd.
#   pe_stratum_months()     each (period x day_type) cell's month = the modal month
#                           of its CALENDAR days (not its sampled days, which an
#                           empty cell does not have).
#   pe_build_effort_strata() the stratification, the fill and the variance basis.
#   pe_empty_cpue_fill()    the matching CPUE fill for an empty stratum.
#
# THE VARIANCE, WRITTEN OUT
#   The stratum total is est_h = N_h * ybar_h and the estimator's variance is
#   Var(est_h) = N_h^2 * Var(ybar_h). What changes is Var(ybar_h):
#
#     n_h >= 2   s_h^2 / n_h from the cell's own sampled days. UNCHANGED, and this is
#                the only case the historical code got right. No finite-population
#                correction: the PE has never applied one, and adding it here would
#                mix an independent (smaller-SE) change into this one. See
#                CHANGE_REGISTER D22.
#     n_h == 1   s_h is undefined. The cell's spread is borrowed from its donor and
#                the divisor stays 1: Var = s_donor^2. This is the collapsed-stratum
#                estimator standard for one-per-stratum designs (Hansen, Hurwitz &
#                Madow 1953; Cochran 1977 sec. 5A.12; Wolter, Introduction to
#                Variance Estimation 2nd ed. 2007 ch. 2). It is conservative for a
#                cell whose true spread is below its donor's and anti-conservative
#                for one above it; both beat the zero it replaces.
#     n_h == 0   the estimate IS the donor mean, so the error has two parts: the
#      (filled)  donor mean's own sampling error, s_donor^2 / n_donor, and the
#                between-cell error, this cell's true mean not being the donor's.
#                The second cannot be estimated from a cell with no samples. The
#                surrogate used is one cell's worth of the donor's day-to-day
#                spread, giving Var = s_donor^2 * (1/n_donor + 1). Deliberately
#                conservative: an imputed number should be visibly wider than a
#                measured one. pe_variance = "donor_mean_only" drops the second
#                term if you want the lower bound instead.
#     n_h == 0   the point estimate is 0, and the error is a BIAS, not a variance.
#     (zeroed)   An SE around a zero would say "0 plus or minus something", implying
#                the truth could be negative and pretending the omission is random.
#                Zeroed cells therefore keep variance 0 and the omitted effort is
#                reported as a bias: pe_zeroed_effort_bias, the effort those cells
#                WOULD carry under "local_day_type". That is the honest reading of
#                pe_empty_effort_stratum = "zero", and the reason it is no longer the
#                shipped default (CHANGE_REGISTER D19).
#
#   Every run reports BOTH the new SE and effort_se_sampled_only, the historical
#   arithmetic, so no earlier number becomes unreachable.
###############################################################################

if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# pe_effort_donors(): the three donor levels, from the SAMPLED days only.
#   daily_effort  one row per sampled day, with est_daily_effort and day_type
#   day_month     event_date -> month, for the (month x day_type) level
# ---------------------------------------------------------------------------
pe_effort_donors <- function(daily_effort, day_month) {
  de <- daily_effort |> dplyr::left_join(day_month, by = "event_date")
  md <- de |> dplyr::group_by(.pe_month, day_type) |>
    dplyr::summarise(md_n = dplyr::n(),
                     md_mean = mean(est_daily_effort, na.rm = TRUE),
                     md_sd   = stats::sd(est_daily_effort, na.rm = TRUE), .groups = "drop")
  dt <- de |> dplyr::group_by(day_type) |>
    dplyr::summarise(dt_n = dplyr::n(),
                     dt_mean = mean(est_daily_effort, na.rm = TRUE),
                     dt_sd   = stats::sd(est_daily_effort, na.rm = TRUE), .groups = "drop")
  list(month_day_type = md, day_type = dt,
       all_n    = sum(is.finite(de$est_daily_effort)),
       all_mean = mean(de$est_daily_effort, na.rm = TRUE),
       all_sd   = stats::sd(de$est_daily_effort, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# pe_stratum_months(): the modal month of each cell's CALENDAR days.
# An empty cell has no sampled days, so its month must come from the calendar or the
# month-local fill has nothing to key on. `days` is the sub-season calendar with
# period, day_type and month.
# ---------------------------------------------------------------------------
pe_stratum_months <- function(days) {
  days |> dplyr::filter(.data$open_section_1) |>
    dplyr::group_by(period, day_type) |>
    dplyr::summarise(.pe_month = as.numeric(names(sort(table(month), decreasing = TRUE))[1]),
                     .groups = "drop")
}

# ---------------------------------------------------------------------------
# pe_build_effort_strata()
#   daily_effort  sampled days: event_date, section_num, period, day_type, est_daily_effort
#   days          the sub-season calendar: event_date, period, day_type, month, open_section_1
#   params        reads pe_empty_effort_stratum and pe_variance
# Returns the stratum frame with est_total, se_total, se_total_sampled_only and the
# per-cell provenance (fill_source, var_source), plus attributes carrying the run-level
# counts the diagnostics writer reports.
# ---------------------------------------------------------------------------
pe_build_effort_strata <- function(daily_effort, days, params) {
  fill <- params$pe_empty_effort_stratum %||% "local_day_type"
  vmode <- params$pe_variance %||% "impute_aware"
  if (!fill %in% c("zero", "day_type", "local_day_type"))
    stop(sprintf("pe_build_effort_strata(): pe_empty_effort_stratum = '%s' is not one of zero | day_type | local_day_type.", fill), call. = FALSE)
  if (!vmode %in% c("sampled_only", "impute_aware", "donor_mean_only"))
    stop(sprintf("pe_build_effort_strata(): pe_variance = '%s' is not one of sampled_only | impute_aware | donor_mean_only.", vmode), call. = FALSE)

  total_days_strat <- days |> dplyr::filter(.data$open_section_1) |>
    dplyr::group_by(period, day_type) |>
    dplyr::summarise(n_total_days = dplyr::n(), .groups = "drop")

  # Built from SAMPLED days, then full-joined to the calendar so a cell with calendar
  # days but no sampled day is VISIBLE rather than silently absent (2026-08-25).
  strat_sampled <- daily_effort |>
    dplyr::group_by(section_num, period, day_type) |>
    dplyr::summarise(mean_daily = mean(est_daily_effort, na.rm = TRUE),
                     sd_daily   = stats::sd(est_daily_effort, na.rm = TRUE),
                     n_sampled  = dplyr::n(), .groups = "drop")

  st <- strat_sampled |>
    dplyr::full_join(total_days_strat, by = c("period", "day_type")) |>
    dplyr::filter(!is.na(n_total_days)) |>
    dplyr::mutate(section_num = tidyr::replace_na(section_num, 1),
                  n_sampled   = tidyr::replace_na(as.integer(n_sampled), 0L),
                  empty_effort_stratum = n_sampled == 0L)

  day_month <- days |> dplyr::select(event_date, .pe_month = month)
  don  <- pe_effort_donors(daily_effort, day_month)
  cell <- pe_stratum_months(days)

  st <- st |>
    dplyr::left_join(cell, by = c("period", "day_type")) |>
    dplyr::left_join(don$month_day_type, by = c(".pe_month", "day_type")) |>
    dplyr::left_join(don$day_type,       by = "day_type") |>
    dplyr::mutate(md_n = tidyr::replace_na(as.numeric(md_n), 0),
                  dt_n = tidyr::replace_na(as.numeric(dt_n), 0))

  # THE MEAN AND THE SPREAD ARE BORROWED FROM DIFFERENT LEVELS, ON PURPOSE.
  #   the MEAN comes from the finest level with ANY sampled day, because that is the most
  #     local estimate of this cell's level and one observation is still an estimate;
  #   the SPREAD comes from the finest level with at least TWO, because sd() of a single
  #     observation does not exist.
  # Coupling them would force a choice between a local mean with no variance and a coarse
  # mean with one. Keeping them separate is also what lets pe_empty_effort_stratum and
  # pe_variance move independent things, so the ladder can attribute them one at a time,
  # and it makes the point estimate under "local_day_type" identical to the 2026-09-11
  # measurement (90,861 on the 2024-25 PE port total) rather than shifting with a
  # variance-driven fallback.
  .local_mean_ok <- if (identical(fill, "day_type")) rep(FALSE, nrow(st)) else (is.finite(st$md_mean) & st$md_n >= 1)
  st <- st |> dplyr::mutate(
    mean_level = dplyr::case_when(.local_mean_ok                  ~ "month x day_type",
                                  is.finite(dt_mean) & dt_n >= 1  ~ "day_type",
                                  TRUE                            ~ "sub-season mean"),
    donor_mean = dplyr::case_when(mean_level == "month x day_type" ~ md_mean,
                                  mean_level == "day_type"         ~ dt_mean,
                                  TRUE                             ~ don$all_mean),
    sd_level   = dplyr::case_when(is.finite(md_sd) & md_n >= 2 ~ "month x day_type",
                                  is.finite(dt_sd) & dt_n >= 2 ~ "day_type",
                                  TRUE                          ~ "sub-season"),
    donor_sd   = dplyr::case_when(sd_level == "month x day_type" ~ md_sd,
                                  sd_level == "day_type"         ~ dt_sd,
                                  TRUE                            ~ don$all_sd),
    donor_n    = dplyr::case_when(sd_level == "month x day_type" ~ md_n,
                                  sd_level == "day_type"         ~ dt_n,
                                  TRUE                            ~ as.numeric(don$all_n)),
    donor_level = mean_level)

  filled <- fill %in% c("day_type", "local_day_type")
  st <- st |> dplyr::mutate(
    imputed_effort_stratum = empty_effort_stratum & filled,
    mean_daily  = dplyr::case_when(!empty_effort_stratum ~ mean_daily,
                                   filled                ~ donor_mean,
                                   TRUE                  ~ 0),
    fill_source = dplyr::case_when(!empty_effort_stratum ~ "sampled",
                                   filled                ~ paste0("imputed from ", donor_level),
                                   TRUE                  ~ "zeroed (no sampled day)"))

  # ---- the variance basis -------------------------------------------------
  # var_sd / var_n are the (spread, divisor) pair the SE formula uses, so every case
  # above is one row of a table rather than a special case in a formula.
  st <- st |> dplyr::mutate(
    var_sd_hist = dplyr::if_else(empty_effort_stratum, NA_real_, sd_daily),
    var_n_hist  = pmax(n_sampled, 1),
    var_sd = var_sd_hist, var_n = as.numeric(var_n_hist),
    var_source = dplyr::if_else(n_sampled >= 2L, "own sampled days",
                                dplyr::if_else(empty_effort_stratum, "none (zero variance)",
                                               "none (single sampled day, sd undefined)")))
  if (vmode %in% c("impute_aware", "donor_mean_only")) {
    st <- st |> dplyr::mutate(
      # single sampled day: collapsed-stratum donor spread, divisor 1
      var_sd = dplyr::if_else(n_sampled == 1L, donor_sd, var_sd),
      var_n  = dplyr::if_else(n_sampled == 1L, 1, var_n),
      var_source = dplyr::if_else(n_sampled == 1L, paste0("collapsed stratum: ", sd_level), var_source),
      # imputed cell: donor-mean variance, plus (impute_aware only) a between-cell term
      var_sd = dplyr::if_else(imputed_effort_stratum, donor_sd, var_sd),
      var_n  = dplyr::if_else(imputed_effort_stratum,
                              if (identical(vmode, "impute_aware"))
                                1 / (1 / pmax(donor_n, 1) + 1) else pmax(donor_n, 1),
                              var_n),
      var_source = dplyr::if_else(imputed_effort_stratum,
                                  paste0(if (identical(vmode, "impute_aware")) "imputed (donor mean + between-cell): "
                                         else "imputed (donor mean only): ", sd_level),
                                  var_source))
  }

  st <- st |> dplyr::mutate(
    est_total = mean_daily * n_total_days,
    se_total  = sqrt((n_total_days^2) * tidyr::replace_na(var_sd^2, 0) / pmax(var_n, 1e-9)),
    se_total_sampled_only = sqrt((n_total_days^2) * tidyr::replace_na(var_sd_hist^2, 0) / var_n_hist),
    # what a zeroed cell omits, so "zero" reports a bias instead of a silent hole
    zeroed_effort_bias = dplyr::if_else(empty_effort_stratum & !filled, donor_mean * n_total_days, 0))

  attr(st, "pe_fill")  <- fill
  attr(st, "pe_variance") <- vmode
  attr(st, "counts") <- list(
    n_empty_strata      = sum(st$empty_effort_stratum),
    n_empty_days        = sum(st$n_total_days[st$empty_effort_stratum], na.rm = TRUE),
    n_single_strata     = sum(st$n_sampled == 1L),
    n_single_days       = sum(st$n_total_days[st$n_sampled == 1L], na.rm = TRUE),
    n_strata_total      = nrow(st),
    n_calendar_days     = sum(st$n_total_days, na.rm = TRUE),
    imputed_effort      = sum(st$est_total[st$imputed_effort_stratum], na.rm = TRUE),
    zeroed_effort_bias  = sum(st$zeroed_effort_bias, na.rm = TRUE))
  st
}

# ---------------------------------------------------------------------------
# pe_empty_cpue_fill(): the CPUE an empty stratum is expanded at.
#   daily_cpue  per sampled day: event_date, catch, hrs
#   strata      the frame from pe_build_effort_strata() (for each cell's .pe_month)
#   params      reads pe_empty_stratum: "local" (shipped) | "pooled" | "zero"
# Returns a numeric vector, one per row of `strata`, plus a `source` attribute.
# The point: under a MONTH-LOCAL effort fill, a sub-season-wide CPUE puts the two
# halves of the same imputed cell on different scales.
# ---------------------------------------------------------------------------
pe_empty_cpue_fill <- function(daily_cpue, strata, days, params) {
  mode <- params$pe_empty_stratum %||% "local"
  if (!mode %in% c("local", "pooled", "zero"))
    stop(sprintf("pe_empty_cpue_fill(): pe_empty_stratum = '%s' is not one of local | pooled | zero.", mode), call. = FALSE)
  pooled <- if (sum(daily_cpue$hrs, na.rm = TRUE) > 0)
    sum(daily_cpue$catch, na.rm = TRUE) / sum(daily_cpue$hrs, na.rm = TRUE) else 0
  if (identical(mode, "zero"))  return(structure(rep(0, nrow(strata)),      source = "zero"))
  if (identical(mode, "pooled")) return(structure(rep(pooled, nrow(strata)), source = "sub-season ratio-of-sums"))
  mc <- daily_cpue |>
    dplyr::left_join(days |> dplyr::select(event_date, .pe_month = month), by = "event_date") |>
    dplyr::group_by(.pe_month) |>
    dplyr::summarise(catch = sum(catch, na.rm = TRUE), hrs = sum(hrs, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(m_cpue = dplyr::if_else(hrs > 0, catch / hrs, NA_real_))
  v <- mc$m_cpue[match(strata$.pe_month, mc$.pe_month)]
  structure(dplyr::coalesce(v, pooled),
            source = sprintf("month ratio-of-sums (%d of %d cells), else sub-season",
                             sum(is.finite(v)), nrow(strata)))
}

# ---------------------------------------------------------------------------
# pe_effort_stratum_report(): the per-run console line and the numbers the
# diagnostics writer carries. Called by both PE runners so the wording cannot drift.
# ---------------------------------------------------------------------------
pe_effort_stratum_report <- function(st, population_name, params) {
  k <- attr(st, "counts"); fill <- attr(st, "pe_fill"); vmode <- attr(st, "pe_variance")
  tot <- sum(st$est_total, na.rm = TRUE)
  se  <- sqrt(sum(st$se_total^2, na.rm = TRUE))
  se0 <- sqrt(sum(st$se_total_sampled_only^2, na.rm = TRUE))
  cat(sprintf(paste0("  PE %s strata: %d of %d cells unsampled (%d of %d days, %.1f%%), %d cells with ONE ",
                     "sampled day (%d days); fill '%s', variance '%s'.\n"),
              population_name, k$n_empty_strata, k$n_strata_total, k$n_empty_days,
              k$n_calendar_days, 100 * k$n_empty_days / max(k$n_calendar_days, 1),
              k$n_single_strata, k$n_single_days, fill, vmode))
  if (k$imputed_effort > 0)
    cat(sprintf("    imputed effort %s of %s (%.1f%%); SE %s (%.1f%%) against %s (%.1f%%) on the sampled-only arithmetic.\n",
                format(round(k$imputed_effort), big.mark = ","), format(round(tot), big.mark = ","),
                100 * k$imputed_effort / max(tot, 1), format(round(se), big.mark = ","),
                100 * se / max(tot, 1), format(round(se0), big.mark = ","), 100 * se0 / max(tot, 1)))
  if (identical(fill, "zero") && k$n_empty_days > 0)
    cat(sprintf(paste0("    *** pe_empty_effort_stratum = 'zero': %d days are expanded at ZERO effort. That is a ",
                       "BIAS, not a variance, and no SE represents it. Those cells would carry about %s effort ",
                       "units under 'local_day_type' (%.1f%% of this component). ***\n"),
                k$n_empty_days, format(round(k$zeroed_effort_bias), big.mark = ","),
                100 * k$zeroed_effort_bias / max(tot, 1)))
  invisible(list(effort_total = tot, effort_se = se, effort_se_sampled_only = se0))
}
