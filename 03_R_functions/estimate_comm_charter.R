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
###############################################################################

estimate_comm_charter <- function(dwg, params) {
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
    result <- list(effort_total=0, Dungeness_Kept=0)
    if(params$estimate_red_rock) result$Red_Rock_Kept <- 0
    return(result)
  }

  mean_dung_per_vessel <- sum(comm_int$dungeness_kept) / nrow(comm_int)
  rr_str <- ""
  if(params$estimate_red_rock) {
    mean_rr_per_vessel <- sum(comm_int$red_rock_kept) / nrow(comm_int)
    rr_str <- sprintf(", Mean Red Rock/vessel: %.1f", mean_rr_per_vessel)
  }

  cat(sprintf("  Tally days: %d, Interviews: %d\n", nrow(tally), nrow(comm_int)))
  cat(sprintf("  Mean Dungeness/vessel: %.1f%s\n", mean_dung_per_vessel, rr_str))

  daily_est <- tally |>
    mutate(
      total_comm_charter = commercial_tally + charter_tally,
      est_dung = total_comm_charter * mean_dung_per_vessel,
      day_of_week = weekdays(date),
      day_type = case_when(
        date %in% crabbing_holiday_dates ~ "weekend",
        day_of_week %in% params$days_wkend ~ "weekend",
        TRUE ~ "weekday"
      )
    )

  if(params$estimate_red_rock) {
    daily_est <- daily_est |> mutate(est_rr = total_comm_charter * mean_rr_per_vessel)
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
    left_join(strat_harvest, by = "day_type") |>
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

  result <- list(
    Dungeness_Kept = total_dung,
    effort_total = total_vessels,
    daily_est = daily_est,
    strat_detail = strat
  )

  if(params$estimate_red_rock) {
    strat <- strat |>
      left_join(daily_est |> group_by(day_type) |>
        summarise(mean_daily_rr = mean(est_rr), .groups = "drop"), by = "day_type") |>
      mutate(est_total_rr = mean_daily_rr * n_total_days)
    result$Red_Rock_Kept <- sum(strat$est_total_rr)
  }

  dung_str <- format(round(total_dung), big.mark=",")
  rr_out <- if(params$estimate_red_rock) sprintf(", Red Rock: %s", format(round(result$Red_Rock_Kept),big.mark=",")) else ""
  cat(sprintf("\n  Est Dungeness (stratified): %s%s\n", dung_str, rr_out))

  return(result)
}
