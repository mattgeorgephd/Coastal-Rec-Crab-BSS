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

  result <- list(
    Dungeness_Kept = total_dung,
    effort_total = total_vessels,
    daily_est = daily_est,
    strat_detail = strat,
    vessel_type_detail = vessel_type_detail
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
