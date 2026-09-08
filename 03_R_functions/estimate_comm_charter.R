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
# estimate_comm_charter.R
#
# Day-type-stratified census expansion of the commercial/charter vessel tally into
# a Dungeness (and optional red-rock) harvest total. Extracted from the pooled and
# gear-resolved drivers so both share one implementation; the gear superset (which
# adds red-rock support, all guarded by params$estimate_red_rock, and a stratified-
# expansion print) is used. With red-rock off the numeric result is identical to the
# pooled inline version. Auto-sourced by both drivers via the 03_R_functions walk.
#
# crabbing_holiday_dates is read from params (params$crabbing_holiday_dates).
#
# REVIEW ITEM 4 (2026-09-08): WHAT IS EXACT AND WHAT IS IMPUTED. The vessel tally on a
# sampled day is a census (charter coordination, vessel monitoring for the commercial
# fleet) and carries no error. But the tally covers only the days a sampler was in port:
# 47 of the 70 days of the 2024-25 window, the missing 23 being 9 Mondays, 8 Sundays,
# 3 Tuesdays, 2 Wednesdays and 1 Thursday (the roster, not closures). Those days are
# filled with the day-type MEAN of the sampled days, and they carry about 4,800 of the
# 11,821 crab. That part is an estimate with sampling error, whatever the tally days are.
# The result now separates the two (observed_dung / imputed_dung), writes a per-day
# table with an `observed` flag (census_daily.csv via the driver), and carries the
# imputation variance
#     Var(total) = sum_h (N_h - n_h)^2 * s_h^2 / n_h            (day-type strata h)
#                + sum_class V_class^2 * s_class^2 / n_class * (1 - n_class / N_class)
# where s_h^2 is the between-day variance of the sampled days' catch in stratum h, and
# the second term is the (near-census, hence tiny) sampling variance of the per-vessel
# catch means applied to the expanded vessel counts V_class. A stratum with one sampled
# day borrows the pooled between-day variance and is flagged.
# params$census_uncertainty: "none" (default, today's behaviour: Dungeness_Kept_se is
# reported but the port total treats the census as a constant) | "imputed_days" (the
# driver adds a normal draw with this SE to the port interval, clamped so the total never
# falls below the observed-exact part). The better fix is a data request: daily vessel
# counts for the unsampled days from the vessel-monitoring / charter records would make
# the tally a true 70-day census and bypass the expansion entirely.
###############################################################################

estimate_comm_charter <- function(dwg, params) {
  # -------------------------------------------------------------------------
  # 2026-09-10: PER-SEASON CENSUS WINDOWS (multi-season spans). The commercial season is
  # a per-season window, so a single census_start/end pair cannot describe a span, and
  # expanding one giant window across the gap between seasons would spread sampled-day
  # means over months with no fishery. params$census_windows, when non-NULL, is a NAMED
  # list season -> c(start, end); the census is estimated per window through this very
  # function (recursion with the scalar keys substituted, so the stratified expansion
  # logic exists once), the totals are summed, and result$by_season keeps each season's
  # own component for the season-summary table. NULL (the default) is the scalar path,
  # unchanged.
  # -------------------------------------------------------------------------
  cw <- params$census_windows
  if (!is.null(cw) && length(cw)) {
    if (is.null(names(cw)) || any(!nzchar(names(cw))))
      stop("estimate_comm_charter(): census_windows must be a NAMED list, season -> c(start, end).",
           call. = FALSE)
    per <- lapply(names(cw), function(sn) {
      ps <- params
      ps$census_windows   <- NULL
      ps$census_start_date <- as.character(cw[[sn]][1])
      ps$census_end_date   <- as.character(cw[[sn]][2])
      r <- estimate_comm_charter(dwg, ps)
      r$season <- sn
      r
    })
    names(per) <- names(cw)
    tot <- per[[1]]
    if (length(per) > 1) for (k in 2:length(per)) {
      tot$effort_total   <- (tot$effort_total   %||% 0) + (per[[k]]$effort_total   %||% 0)
      tot$Dungeness_Kept <- (tot$Dungeness_Kept %||% 0) + (per[[k]]$Dungeness_Kept %||% 0)
      if (!is.null(tot$Red_Rock_Kept) || !is.null(per[[k]]$Red_Rock_Kept))
        tot$Red_Rock_Kept <- (tot$Red_Rock_Kept %||% 0) + (per[[k]]$Red_Rock_Kept %||% 0)
      # review item 4: the observed/imputed split and the variance add across windows
      # (independent seasons), the per-day tables stack.
      tot$observed_dung  <- (tot$observed_dung  %||% 0) + (per[[k]]$observed_dung  %||% 0)
      tot$imputed_dung   <- (tot$imputed_dung   %||% 0) + (per[[k]]$imputed_dung   %||% 0)
      tot$Dungeness_Kept_var <- (tot$Dungeness_Kept_var %||% 0) + (per[[k]]$Dungeness_Kept_var %||% 0)
      if (!is.null(tot$daily_full) && !is.null(per[[k]]$daily_full))
        tot$daily_full <- dplyr::bind_rows(tot$daily_full, per[[k]]$daily_full)
      if (!is.null(tot$variance_detail) && !is.null(per[[k]]$variance_detail))
        tot$variance_detail <- dplyr::bind_rows(tot$variance_detail, per[[k]]$variance_detail)
    }
    tot$Dungeness_Kept_se <- sqrt(tot$Dungeness_Kept_var %||% 0)
    tot$season    <- NULL
    tot$by_season <- per
    return(tot)
  }

  # Holidays from the centralized config (single source of truth).
  crabbing_holiday_dates <- params$crabbing_holiday_dates
  cat("\n--- Commercial/Charter Census Estimation (Stratified) ---\n")

  census_start <- as.Date(params$census_start_date)
  census_end   <- as.Date(params$census_end_date)

  tally <- dwg$comm_tally |>
    filter(between(date, census_start, census_end))

  comm_int <- dwg$interview |>
    filter(population == "comm_charter",
           between(event_date, census_start, census_end))

  if(nrow(comm_int) == 0 || nrow(tally) == 0) {
    cat("  No commercial/charter data available.\n")
    result <- list(effort_total=0, Dungeness_Kept=0,
                   observed_dung = 0, imputed_dung = 0, Dungeness_Kept_var = 0, Dungeness_Kept_se = 0,
                   census_uncertainty = tolower(params$census_uncertainty %||% "none"))
    if(params$estimate_red_rock) result$Red_Rock_Kept <- 0
    return(result)
  }

  # T1.4 (2026-07-12): separate per-vessel catch means for commercial vs charter
  # vessels, each applied to its own tally column, instead of one mean across the
  # pooled tally. The two classes have materially different catch profiles (2024-25
  # census window: commercial ~35 vs charter ~51 Dungeness/vessel), so a single mean
  # biases the total whenever the commercial:charter mix in the tally differs from
  # the mix in the interviews. boat_type_clean is set in fetch_crab_data (Commercial,
  # Charter, with Guide folded into Charter). A class with no interviews in the
  # window falls back to the pooled mean. (Backlog T1.4 / critique 6.)
  pooled_mean_dung <- sum(comm_int$dungeness_kept) / nrow(comm_int)
  class_mean <- function(cls, col, fallback) {
    if (!"boat_type_clean" %in% names(comm_int)) return(fallback)  # degrade to pooled mean
    sub <- comm_int |> filter(boat_type_clean == cls)
    if (nrow(sub) > 0) sum(sub[[col]], na.rm = TRUE) / nrow(sub) else fallback
  }
  n_comm  <- sum(comm_int$boat_type_clean == "Commercial", na.rm = TRUE)
  n_char  <- sum(comm_int$boat_type_clean == "Charter",    na.rm = TRUE)
  md_comm <- class_mean("Commercial", "dungeness_kept", pooled_mean_dung)
  md_char <- class_mean("Charter",    "dungeness_kept", pooled_mean_dung)

  rr_str <- ""
  if(params$estimate_red_rock) {
    pooled_mean_rr <- sum(comm_int$red_rock_kept) / nrow(comm_int)
    mr_comm <- class_mean("Commercial", "red_rock_kept", pooled_mean_rr)
    mr_char <- class_mean("Charter",    "red_rock_kept", pooled_mean_rr)
    rr_str  <- sprintf(" | Red Rock/vessel: comm %.1f, charter %.1f", mr_comm, mr_char)
  }

  cat(sprintf("  Tally days: %d, Interviews: %d (commercial %d, charter %d)\n",
              nrow(tally), nrow(comm_int), n_comm, n_char))
  cat(sprintf("  Mean Dungeness/vessel: commercial %.1f, charter %.1f%s\n",
              md_comm, md_char, rr_str))

  daily_est <- tally |>
    mutate(
      total_comm_charter = commercial_tally + charter_tally,
      est_dung = commercial_tally * md_comm + charter_tally * md_char,
      day_of_week = weekdays(date),
      day_type = case_when(
        date %in% crabbing_holiday_dates ~ "weekend",
        day_of_week %in% params$days_wkend ~ "weekend",
        TRUE ~ "weekday"
      )
    )

  if(params$estimate_red_rock) {
    daily_est <- daily_est |> mutate(est_rr = commercial_tally * mr_comm + charter_tally * mr_char)
  }

  census_calendar <- tibble(
    date = seq.Date(census_start, census_end, by = "day"),
    day_of_week = weekdays(date),
    day_type = case_when(
      date %in% crabbing_holiday_dates ~ "weekend",
      day_of_week %in% params$days_wkend ~ "weekend",
      TRUE ~ "weekday"
    )
  )

  total_by_type <- census_calendar |> count(day_type, name = "n_total_days")
  sampled_by_type <- daily_est |> count(day_type, name = "n_sampled_days")

  strat <- total_by_type |>
    left_join(sampled_by_type, by = "day_type") |>
    mutate(n_sampled_days = replace_na(n_sampled_days, 0))

  strat_harvest <- daily_est |>
    group_by(day_type) |>
    summarise(
      mean_daily_dung = mean(est_dung),
      mean_daily_vessels = mean(total_comm_charter),
      .groups = "drop"
    )

  if(params$estimate_red_rock) {
    strat_harvest_rr <- daily_est |>
      group_by(day_type) |>
      summarise(mean_daily_rr = mean(est_rr), .groups = "drop")
    strat_harvest <- strat_harvest |> left_join(strat_harvest_rr, by = "day_type")
  }

  strat <- strat |>
    left_join(strat_harvest, by = "day_type")
  # Review item 4 (2026-09-08): a day type with NO sampled day used to expand to NA and
  # take the whole component with it (a short per-season window can hit this). It now
  # borrows the pooled sampled-day mean, and says so; its imputation variance below uses
  # the pooled between-day variance over all sampled days.
  .no_sample <- strat$n_sampled_days == 0
  if (any(.no_sample)) {
    strat$mean_daily_dung[.no_sample]    <- mean(daily_est$est_dung)
    strat$mean_daily_vessels[.no_sample] <- mean(daily_est$total_comm_charter)
    if (params$estimate_red_rock && "mean_daily_rr" %in% names(strat))
      strat$mean_daily_rr[.no_sample] <- mean(daily_est$est_rr)
    cat(sprintf("  NOTE: no sampled day in the %s stratum; its %d days take the pooled sampled-day mean.\n",
                paste(strat$day_type[.no_sample], collapse = "/"), sum(strat$n_total_days[.no_sample])))
  }
  strat <- strat |>
    mutate(
      est_total_dung = mean_daily_dung * n_total_days,
      est_total_vessels = mean_daily_vessels * n_total_days
    )

  cat("\n  Stratified expansion by day type:\n")
  cat(sprintf("    %-10s  Sampled  Total  Mean/day  Expanded\n", "Day Type"))
  for(i in 1:nrow(strat)) {
    cat(sprintf("    %-10s  %5d    %5d  %7.1f   %8.0f\n",
                strat$day_type[i], strat$n_sampled_days[i], strat$n_total_days[i],
                strat$mean_daily_dung[i], strat$est_total_dung[i]))
  }

  total_dung <- sum(strat$est_total_dung)
  total_vessels <- sum(strat$est_total_vessels)

  # T1.4: per-vessel-type split summary for the report (a sampled-day census, i.e.
  # before the day-type expansion that produces the headline total_dung).
  comm_tally_n <- sum(tally$commercial_tally, na.rm = TRUE)
  char_tally_n <- sum(tally$charter_tally,    na.rm = TRUE)
  vessel_type_detail <- tibble(
    vessel_type          = c("Commercial", "Charter"),
    n_interviews         = c(n_comm, n_char),
    mean_dung_per_vessel = c(md_comm, md_char),
    tally_vessels        = c(comm_tally_n, char_tally_n),
    sampled_catch        = c(comm_tally_n * md_comm, char_tally_n * md_char)
  )

  # --- REVIEW ITEM 4 (2026-09-08): observed vs imputed, and the imputation variance ------
  # Every calendar day gets a row: a sampled day carries its own (exact) tally estimate,
  # an unsampled day carries its stratum mean. The stratum variance is the between-day
  # variance of the SAMPLED days' catch, so Var(sum of imputed days in h) = (N_h - n_h)^2 *
  # s_h^2 / n_h treats the imputed days as N_h - n_h copies of one estimated mean, which is
  # what the expansion does. A stratum with a single sampled day cannot estimate s_h^2 and
  # borrows the pooled between-day variance (flagged in variance_detail).
  census_mode <- tolower(params$census_uncertainty %||% "none")
  if (!census_mode %in% c("none", "imputed_days"))
    stop("params$census_uncertainty must be 'none' or 'imputed_days' (got '", census_mode, "')", call. = FALSE)
  s2_pooled <- if (nrow(daily_est) > 1) stats::var(daily_est$est_dung) else 0
  var_by_type <- daily_est |>
    group_by(day_type) |>
    summarise(n_sampled_days = n(), s2_day = if (n() > 1) stats::var(est_dung) else NA_real_, .groups = "drop")
  variance_detail <- strat |>
    select(day_type, n_total_days, n_sampled_days, mean_daily_dung) |>
    left_join(var_by_type |> select(day_type, s2_day), by = "day_type") |>
    mutate(
      s2_source     = ifelse(is.na(s2_day) & n_sampled_days > 0, "pooled (1 sampled day)", ifelse(n_sampled_days == 0, "none (no sampled day)", "stratum")),
      s2_day        = ifelse(is.na(s2_day), s2_pooled, s2_day),
      n_imputed_days = n_total_days - n_sampled_days,
      imputed_dung  = mean_daily_dung * n_imputed_days,
      observed_dung = mean_daily_dung * n_sampled_days,
      var_imputed   = ifelse(n_sampled_days > 0, n_imputed_days^2 * s2_day / n_sampled_days,
                             n_imputed_days^2 * s2_pooled / max(nrow(daily_est), 1))
    )
  # per-vessel catch means: sample means of a near-census interview set (FPC applied)
  .mean_var <- function(cls, col) {
    if (!"boat_type_clean" %in% names(comm_int)) return(0)
    x <- comm_int[[col]][comm_int$boat_type_clean == cls]; x <- x[is.finite(x)]
    n <- length(x); if (n < 2) return(0)
    N <- if (cls == "Commercial") comm_tally_n else char_tally_n
    fpc <- if (is.finite(N) && N > n) (1 - n / N) else 0
    stats::var(x) / n * fpc
  }
  V_comm <- sum((strat$n_total_days / pmax(strat$n_sampled_days, 1)) *
                  vapply(strat$day_type, function(dt) sum(daily_est$commercial_tally[daily_est$day_type == dt], na.rm = TRUE), numeric(1)))
  V_char <- sum((strat$n_total_days / pmax(strat$n_sampled_days, 1)) *
                  vapply(strat$day_type, function(dt) sum(daily_est$charter_tally[daily_est$day_type == dt], na.rm = TRUE), numeric(1)))
  var_means <- V_comm^2 * .mean_var("Commercial", "dungeness_kept") + V_char^2 * .mean_var("Charter", "dungeness_kept")
  var_total <- sum(variance_detail$var_imputed) + var_means
  observed_dung <- sum(variance_detail$observed_dung)
  imputed_dung  <- sum(variance_detail$imputed_dung)
  daily_full <- census_calendar |>
    select(date, day_type) |>
    left_join(daily_est |> select(date, commercial_tally, charter_tally, total_comm_charter, est_dung), by = "date") |>
    left_join(strat |> select(day_type, mean_daily_dung, mean_daily_vessels), by = "day_type") |>
    mutate(observed = !is.na(est_dung),
           est_dung = ifelse(observed, est_dung, mean_daily_dung),
           est_vessels = ifelse(observed, total_comm_charter, mean_daily_vessels),
           source = ifelse(observed, "tally (census)", "imputed (day-type mean)")) |>
    select(date, day_type, observed, source, commercial_tally, charter_tally, est_vessels, est_dung)
  cat(sprintf(paste0("  Observed exactly on %d sampled days: %s crab; imputed on %d unsampled days: %s crab;",
                     " imputation SE %s (%.1f%% of the total); census_uncertainty = '%s'%s\n"),
              sum(strat$n_sampled_days), format(round(observed_dung), big.mark = ","),
              sum(variance_detail$n_imputed_days), format(round(imputed_dung), big.mark = ","),
              format(round(sqrt(var_total)), big.mark = ","), 100 * sqrt(var_total) / max(total_dung, 1), census_mode,
              if (census_mode == "none") " (reported, not carried into the port interval)" else " (carried into the port interval)"))
  if (any(variance_detail$s2_source != "stratum"))
    cat("    NOTE: ", paste(sprintf("%s stratum: %s", variance_detail$day_type[variance_detail$s2_source != "stratum"],
                                    variance_detail$s2_source[variance_detail$s2_source != "stratum"]), collapse = "; "), "\n")

  result <- list(
    Dungeness_Kept = total_dung,
    effort_total = total_vessels,
    daily_est = daily_est,
    strat_detail = strat,
    vessel_type_detail = vessel_type_detail,
    # review item 4
    observed_dung = observed_dung,
    imputed_dung = imputed_dung,
    Dungeness_Kept_var = var_total,
    Dungeness_Kept_se = sqrt(var_total),
    census_uncertainty = census_mode,
    daily_full = daily_full,
    variance_detail = variance_detail
  )

  if(params$estimate_red_rock) {
    strat <- strat |>
      mutate(est_total_rr = mean_daily_rr * n_total_days)   # mean_daily_rr joined above (pooled fallback applied)
    result$Red_Rock_Kept <- sum(strat$est_total_rr)
  }

  dung_str <- format(round(total_dung), big.mark=",")
  rr_out <- if(params$estimate_red_rock) sprintf(", Red Rock: %s", format(round(result$Red_Rock_Kept),big.mark=",")) else ""
  cat(sprintf("\n  Est Dungeness (stratified): %s%s\n", dung_str, rr_out))

  return(result)
}
