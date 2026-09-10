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
# The commercial/charter component of the port total. TWO POPULATIONS, TWO ESTIMATORS
# since 2026-09-11 (Matt's correction: "a complete census without error applies to the
# commercial boats; the charter vessels are not 100% sampled and need expansion"):
#
#   COMMERCIAL vessels fishing as recreational -- A CENSUS. The daily tally counts the
#       vessels (vessel monitoring for the commercial fleet) and the samplers are
#       scheduled on the days those vessels are confirmed to operate, so a window day
#       without a tally is a day with NO OPERATION. The component is the exact sum over
#       the tally days, vessels x mean catch per vessel-trip. What is left is that MEAN,
#       a near-census sample mean (2024-25: 141 interviews on 164 vessel-trips, 86%),
#       reported with a finite-population-corrected SE and treated as a constant in the
#       port interval unless census_uncertainty = "sampling".
#   CHARTER vessels -- AN EXPANSION. They are not fully sampled, so the component is
#       N x (mean catch per interviewed trip) with a sampling variance. The frame N is
#       the charter TRIP roster the operators report (charter_trips.xlsx: one row per
#       trip, interviewed / missed / canceled), unioned per day with the tally's charter
#       column, because each frame is incomplete in its own way (the tally sees only the
#       days a sampler was in port; the roster only the operators who reported). The
#       expansion is stratified by VESSEL where the roster attributes the trip and the
#       vessel has interviews, since the vessels differ (2024-25 Westport: Ultimate 70.2
#       crab/trip on 17 interviews, Outta Line 30.7 on 3). 2024-25, Dec 1 to Feb 8:
#       34 trips in the frame, 20 interviewed (59%), 1,286 crab observed, 2,133 expanded,
#       SE 75. That SE IS carried into the port interval (census_uncertainty = "charter",
#       the shipped default).
#
# THE ESTIMATOR, per charter stratum (a vessel, or one pooled stratum):
#     est_v = N_v * m_v                      m_v = mean catch per interviewed trip
#     var_v = N_v^2 (1 - n_v / N_v) s_v^2 / n_v          (SRS with the FPC)
# A vessel with fewer than two interviews, or none, borrows the pooled charter mean and
# SD and the pooled n, and says so in charter_detail$mean_source. n_v >= N_v (every trip
# interviewed) gives zero: that stratum is itself a census.
#
# WHY NOT THE SPREADSHEET'S ARITHMETIC. The "charter trips" sheet estimates the unobserved
# trips per vessel by hand, sometimes at the observed mean per trip and sometimes at
# (crabbers x the daily limit); the two rules disagree by up to 45% per vessel and the
# result has no variance. On 2024-25 Westport the sheet's window total (~2,156 crab) and
# this estimator (2,133) agree to about 1%, so the choice does not move the number -- it
# gives the number a defensible variance. charter_expansion = "pooled" ignores the vessel
# strata (2,186 on 2024-25).
#
# HISTORY
#   2026-09-08 (review item 4): the tally covers only the days a sampler was in port (47
#       of 70 on 2024-25); the result separates observed from imputed, carries the
#       day-type imputation variance sum_h (N_h - n_h)^2 s_h^2 / n_h, and returns
#       daily_full (census_daily.csv) and variance_detail (census_variance.csv).
#   2026-09-09: the unsampled days had no operation, so census_expansion = "none"
#       (shipped) makes the commercial component exact; "day_type" keeps the 2026-09-08
#       expansion for reproduction and for a season whose roster did not track operations.
#   2026-09-10: the charter trip roster read as a second frame (charter_frame).
#   2026-09-11: the charter component became the expansion above; census_uncertainty
#       gained "charter" and ships as the default.
#
# CONFIGURATION
#   census_expansion   "none" (default) | "day_type"        the COMMERCIAL part
#   charter_frame      "roster" (default) | "tally"         the CHARTER frame
#   charter_expansion  "vessel" (default) | "pooled"        the CHARTER strata
#   census_uncertainty "charter" (default) | "none" | "sampling" ("imputed_days" is
#                      accepted as the 2026-09-08 name for "sampling")
#
# WITHOUT A ROSTER (charter_frame = "tally", or a season with no roster rows) the charter
# frame is the tally's charter column and follows census_expansion exactly as the
# commercial part does, which reproduces the pre-2026-09-10 arithmetic (7,884 on 2024-25);
# the interview-sampling variance is still carried. The frame is then incomplete by
# construction: days without a tally contribute no charter trips, and on 2024-25 the
# roster showed 8 of 31 sailed trips (26%) fell on exactly those days. The run says so.
#
# WHAT THE RESULT CARRIES. Dungeness_Kept_var / _se are the TOTAL variance of the
# component (reported, census_variance.csv). carried_var / carried_se are the part
# census_uncertainty puts into the port interval, and are what the drivers draw with.
# commercial_* and charter_* report the two components separately.
###############################################################################

estimate_comm_charter <- function(dwg, params) {
  # -------------------------------------------------------------------------
  # PER-SEASON CENSUS WINDOWS (2026-09-10, multi-season spans). The commercial season is
  # a per-season window, so one census_start/end pair cannot describe a span, and
  # expanding a single giant window across the gap between seasons would spread
  # sampled-day means over months with no fishery. params$census_windows, when non-NULL,
  # is a NAMED list season -> c(start, end): the component is estimated per window through
  # this very function (recursion with the scalar keys substituted, so the estimator
  # exists once), the totals and variances are summed, and by_season keeps each season's
  # own result for the season-summary table. NULL (the default) is the scalar path.
  # -------------------------------------------------------------------------
  cw <- params$census_windows
  if (!is.null(cw) && length(cw)) {
    if (is.null(names(cw)) || any(!nzchar(names(cw))))
      stop("estimate_comm_charter(): census_windows must be a NAMED list, season -> c(start, end).", call. = FALSE)
    per <- lapply(names(cw), function(sn) {
      ps <- params
      ps$census_windows    <- NULL
      ps$census_start_date <- as.character(cw[[sn]][1])
      ps$census_end_date   <- as.character(cw[[sn]][2])
      r <- estimate_comm_charter(dwg, ps); r$season <- sn; r
    })
    names(per) <- names(cw)
    tot <- per[[1]]
    sum_keys <- c("effort_total", "Dungeness_Kept", "Red_Rock_Kept", "observed_dung", "imputed_dung",
                  "Dungeness_Kept_var", "carried_var", "imputation_var",
                  "commercial_dung", "commercial_vessels", "commercial_var",
                  "charter_dung", "charter_observed_dung", "charter_var", "charter_trips", "charter_interviews", "commercial_interviews",
                  "charter_roster_dung", "n_unsampled_days", "n_roster_only_days", "n_roster_trips")
    bind_keys <- c("daily_full", "daily_est", "variance_detail", "roster_reconciliation", "charter_detail", "strat_detail")
    if (length(per) > 1) for (i in 2:length(per)) {
      for (k in sum_keys) if (!is.null(tot[[k]]) || !is.null(per[[i]][[k]]))
        tot[[k]] <- (tot[[k]] %||% 0) + (per[[i]][[k]] %||% 0)
      for (k in bind_keys) if (!is.null(per[[i]][[k]]))
        tot[[k]] <- if (is.null(tot[[k]])) per[[i]][[k]] else dplyr::bind_rows(tot[[k]], per[[i]][[k]])
    }
    tot$Dungeness_Kept_se <- sqrt(tot$Dungeness_Kept_var %||% 0)
    tot$carried_se        <- sqrt(tot$carried_var %||% 0)
    tot$commercial_se     <- sqrt(tot$commercial_var %||% 0)
    tot$charter_se        <- sqrt(tot$charter_var %||% 0)
    tot$charter_sampled_frac <- if ((tot$charter_trips %||% 0) > 0) (tot$charter_interviews %||% 0) / tot$charter_trips else NA_real_
    tot$charter_frame <- paste(unique(vapply(per, function(r) r$charter_frame %||% "tally", character(1))), collapse = "+")
    tot$season <- NULL
    tot$by_season <- per
    return(tot)
  }

  crabbing_holiday_dates <- params$crabbing_holiday_dates      # from the centralized config
  cat("\n--- Commercial/Charter Estimation (commercial census + charter expansion) ---\n")

  census_start <- as.Date(params$census_start_date)
  census_end   <- as.Date(params$census_end_date)

  census_expansion <- tolower(params$census_expansion %||% "none")
  if (!census_expansion %in% c("none", "day_type"))
    stop("params$census_expansion must be 'none' or 'day_type' (got '", census_expansion, "')", call. = FALSE)
  charter_frame_req <- tolower(params$charter_frame %||% "roster")
  if (!charter_frame_req %in% c("roster", "tally"))
    stop("params$charter_frame must be 'roster' or 'tally' (got '", charter_frame_req, "')", call. = FALSE)
  charter_expansion_req <- tolower(params$charter_expansion %||% "vessel")
  if (!charter_expansion_req %in% c("vessel", "pooled"))
    stop("params$charter_expansion must be 'vessel' or 'pooled' (got '", charter_expansion_req, "')", call. = FALSE)
  census_mode <- tolower(params$census_uncertainty %||% "charter")
  if (identical(census_mode, "imputed_days")) census_mode <- "sampling"     # the 2026-09-08 name
  if (!census_mode %in% c("none", "charter", "sampling"))
    stop("params$census_uncertainty must be 'none', 'charter' or 'sampling' (got '", census_mode, "')", call. = FALSE)

  tally <- dwg$comm_tally |>
    filter(between(date, census_start, census_end)) |>
    mutate(commercial_tally = replace_na(as.numeric(commercial_tally), 0),
           charter_tally    = replace_na(as.numeric(charter_tally), 0))
  comm_int <- dwg$interview |> filter(population == "comm_charter", between(event_date, census_start, census_end))

  roster <- dwg$charter_roster
  roster <- if (charter_frame_req == "roster" && !is.null(roster) && nrow(roster))
    roster |> filter(between(as.Date(date), census_start, census_end)) else NULL
  use_roster <- !is.null(roster) && nrow(roster) > 0
  charter_frame <- if (use_roster) "roster" else "tally"

  .day_type <- function(d) case_when(d %in% crabbing_holiday_dates ~ "weekend",
                                     weekdays(d) %in% params$days_wkend ~ "weekend", TRUE ~ "weekday")
  census_calendar <- tibble(date = seq.Date(census_start, census_end, by = "day")) |>
    mutate(day_of_week = weekdays(date), day_type = .day_type(date))

  if (nrow(comm_int) == 0 || (nrow(tally) == 0 && !use_roster)) {
    if (nrow(comm_int) > 0 && nrow(tally) == 0)
      warning(sprintf(paste0("estimate_comm_charter(): %d commercial/charter interview(s) in the window %s to %s but NO vessel tally ",
                             "rows (and no charter roster): the frame is missing for this window and the component is 0. A tally was ",
                             "kept from 2024-25 on; see 04_input_files/build_comm_charter_tally.R."),
                      nrow(comm_int), census_start, census_end), call. = FALSE)
    cat("  No commercial/charter data available.\n")
    out <- list(effort_total = 0, Dungeness_Kept = 0, observed_dung = 0, imputed_dung = 0,
                Dungeness_Kept_var = 0, Dungeness_Kept_se = 0, carried_var = 0, carried_se = 0, imputation_var = 0,
                census_uncertainty = census_mode, census_expansion = census_expansion, n_unsampled_days = 0L,
                charter_frame = charter_frame, charter_expansion = charter_expansion_req,
                commercial_dung = 0, commercial_vessels = 0, commercial_var = 0, commercial_se = 0,
                charter_dung = 0, charter_observed_dung = 0, charter_var = 0, charter_se = 0,
                charter_trips = 0, charter_interviews = 0L, charter_sampled_frac = NA_real_, commercial_interviews = 0L,
                charter_roster_dung = 0, n_roster_only_days = 0L, n_roster_trips = 0L)
    # the window's days still go into census_daily.csv, flagged, so the gap is documented
    # rather than absent (a 2023-24 window has interviews but no vessel tally at all)
    out$daily_full <- census_calendar |>
      transmute(date, day_type, observed = FALSE,
                source = if (nrow(comm_int) > 0) "no frame (no vessel tally or roster this season)" else "no commercial/charter data",
                commercial_tally = NA_real_, charter_tally = NA_real_, charter_trips = 0,
                est_vessels = 0, est_dung_char = 0, est_dung = 0)
    out$variance_detail <- tibble(day_type = character(), n_total_days = integer(), n_calendar_days = integer(),
                                  n_expand_days = integer(), n_sampled_days = integer(), n_imputed_days = integer(),
                                  mean_daily_dung = numeric(), s2_day = numeric(), s2_source = character(),
                                  n_unsampled_days = integer(), imputed_dung = numeric(), observed_dung = numeric(),
                                  var_imputed = numeric())
    if (isTRUE(params$estimate_red_rock)) out$Red_Rock_Kept <- 0
    return(out)
  }

  # ---- per-class catch means (T1.4, 2026-07-12) ------------------------------------
  # The classes have materially different catch profiles, so one pooled mean would bias
  # the total whenever the class mix in the frame differs from the mix in the interviews.
  # A class with no interviews falls back to the pooled mean.
  pooled_mean_dung <- sum(comm_int$dungeness_kept) / nrow(comm_int)
  .cls <- function(cls) if ("boat_type_clean" %in% names(comm_int)) comm_int |> filter(boat_type_clean == cls) else comm_int[0, , drop = FALSE]
  ci_comm <- .cls("Commercial"); ci_char <- .cls("Charter")
  n_comm <- nrow(ci_comm); n_char <- nrow(ci_char)
  .m <- function(d, col, fb) if (nrow(d)) sum(d[[col]], na.rm = TRUE) / nrow(d) else fb
  md_comm <- .m(ci_comm, "dungeness_kept", pooled_mean_dung)
  md_char <- .m(ci_char, "dungeness_kept", pooled_mean_dung)
  sd_char <- if (n_char > 1) stats::sd(ci_char$dungeness_kept) else NA_real_
  charter_obs <- if (n_char) sum(ci_char$dungeness_kept) else 0
  if (isTRUE(params$estimate_red_rock)) {
    pooled_mean_rr <- sum(comm_int$red_rock_kept) / nrow(comm_int)
    mr_comm <- .m(ci_comm, "red_rock_kept", pooled_mean_rr)
    mr_char <- .m(ci_char, "red_rock_kept", pooled_mean_rr)
  }
  cat(sprintf("  Tally days: %d, Interviews: %d (commercial %d, charter %d)\n", nrow(tally), nrow(comm_int), n_comm, n_char))
  cat(sprintf("  Mean Dungeness per vessel-trip: commercial %.1f, charter %.1f%s\n", md_comm, md_char,
              if (isTRUE(params$estimate_red_rock)) sprintf(" | Red Rock: comm %.1f, charter %.1f", mr_comm, mr_char) else ""))

  # =========================================================================
  # 1. THE CHARTER TRIP FRAME
  # roster frame: the operators' sailed trips per (date, vessel), unioned per day with the
  #     tally's charter column; the surplus on a day where the tally exceeds the roster is
  #     an unattributed trip (a vessel that did not report, or a roster gap).
  # tally frame:  the tally's charter column on the tally days, day-type expanded with the
  #     commercial part when census_expansion = "day_type" (the pre-2026-09-10 behaviour).
  # =========================================================================
  roster_recon <- NULL; n_roster_only_days <- 0L; n_roster_trips <- 0L; charter_roster_dung <- 0
  if (use_roster) {
    r_dv <- roster |>
      mutate(date = as.Date(date), status = tolower(as.character(status)),
             vessel = ifelse(is.na(vessel) | !nzchar(str_squish(as.character(vessel))), "(unnamed)", str_squish(as.character(vessel)))) |>
      filter(!status %in% "canceled") |>
      count(date, vessel, name = "trips")
    day_union <- full_join(tally |> select(date, charter_tally),
                           r_dv |> group_by(date) |> summarise(roster_trips = sum(trips), .groups = "drop"), by = "date") |>
      mutate(charter_tally = replace_na(charter_tally, 0), roster_trips = replace_na(roster_trips, 0),
             surplus = pmax(0, charter_tally - roster_trips), charter_n = roster_trips + surplus)
    frame_day <- bind_rows(r_dv,
                           day_union |> filter(surplus > 0) |> transmute(date, vessel = NA_character_, trips = surplus)) |>
      arrange(date, vessel)
    n_roster_trips <- as.integer(sum(r_dv$trips))
    only <- day_union |> filter(!date %in% tally$date)
    n_roster_only_days <- nrow(only)
    charter_roster_dung <- sum(only$charter_n) * md_char        # crab a tally-only frame would carry as zero
    roster_recon <- bind_rows(
      day_union |> filter(date %in% tally$date, roster_trips != charter_tally) |>
        transmute(date, tally_charter = charter_tally, roster_trips, charter_used = charter_n, note = "tally day; frames disagree"),
      only |> transmute(date, tally_charter = NA_real_, roster_trips, charter_used = charter_n, note = "no tally; roster trips")) |>
      arrange(date)
    cat(sprintf(paste0("  Charter frame: the trip roster, unioned with the tally. %d sailed roster trip(s); %d trip(s) on %d day(s) with ",
                       "no tally (a tally-only frame carries 0 there); %d tally day(s) where the frames disagree.\n"),
                n_roster_trips, as.integer(sum(only$charter_n)), n_roster_only_days,
                sum(grepl("frames disagree", roster_recon$note))))
    if (nrow(roster_recon)) {
      cat("    date         tally  roster  used  note\n")
      for (i in seq_len(nrow(roster_recon)))
        cat(sprintf("    %s   %5s  %6.0f  %4.0f  %s\n", roster_recon$date[i],
                    ifelse(is.na(roster_recon$tally_charter[i]), "-", sprintf("%.0f", roster_recon$tally_charter[i])),
                    roster_recon$roster_trips[i], roster_recon$charter_used[i], roster_recon$note[i]))
    }
  } else {
    frame_day <- tally |> filter(charter_tally > 0) |> transmute(date, vessel = NA_character_, trips = charter_tally)
    cat(sprintf(paste0("  Charter frame: the tally's charter column%s. The frame is incomplete by construction: a day without a tally ",
                       "carries no charter trip, and on 2024-25 the roster showed 8 of 31 sailed trips (26%%) fell on such days.\n"),
                if (charter_frame_req == "roster") " (no charter roster rows in this window)" else ""))
  }

  # =========================================================================
  # 2. THE COMMERCIAL CENSUS, and the day-type strata
  # est_dung_exp is what the day-type expansion acts on: the commercial estimate always,
  # plus the charter estimate when the charter rides on the tally frame. Splitting the two
  # and expanding each by day type gives the same total as expanding their sum (linearity),
  # so the tally-frame arithmetic is bit-for-bit the pre-2026-09-10 one.
  # =========================================================================
  ride_along <- !use_roster
  daily_est <- tally |>
    left_join(frame_day |> group_by(date) |> summarise(charter_n = sum(trips), .groups = "drop"), by = "date") |>
    mutate(charter_n = replace_na(charter_n, 0),
           total_comm_charter = commercial_tally + charter_n,
           est_dung_comm = commercial_tally * md_comm,
           est_dung_char = charter_n * md_char,
           est_dung_exp  = est_dung_comm + if (ride_along) est_dung_char else 0,
           est_dung      = est_dung_comm + est_dung_char,
           day_of_week   = weekdays(date),
           day_type      = .day_type(date))
  if (isTRUE(params$estimate_red_rock))
    daily_est <- daily_est |> mutate(est_rr_exp = commercial_tally * mr_comm + if (ride_along) charter_n * mr_char else 0)

  strat <- census_calendar |> count(day_type, name = "n_total_days") |>
    left_join(daily_est |> count(day_type, name = "n_sampled_days"), by = "day_type") |>
    mutate(n_sampled_days = replace_na(n_sampled_days, 0)) |>
    left_join(daily_est |> group_by(day_type) |>
                summarise(mean_daily_dung    = mean(est_dung_exp),
                          mean_daily_comm_v  = mean(commercial_tally),
                          mean_daily_char_v  = mean(charter_n),
                          mean_daily_vessels = mean(commercial_tally + if (ride_along) charter_n else 0),
                          .groups = "drop"), by = "day_type")
  if (isTRUE(params$estimate_red_rock))
    strat <- strat |> left_join(daily_est |> group_by(day_type) |> summarise(mean_daily_rr = mean(est_rr_exp), .groups = "drop"), by = "day_type")
  # Review item 4 (2026-09-08): a day type with NO sampled day used to expand to NA and
  # take the whole component with it (a short per-season window can hit this). It borrows
  # the pooled sampled-day mean and says so; its imputation variance uses the pooled one.
  .ns <- strat$n_sampled_days == 0
  if (any(.ns)) {
    strat$mean_daily_dung[.ns]    <- mean(daily_est$est_dung_exp)
    strat$mean_daily_comm_v[.ns]  <- mean(daily_est$commercial_tally)
    strat$mean_daily_char_v[.ns]  <- mean(daily_est$charter_n)
    strat$mean_daily_vessels[.ns] <- mean(daily_est$commercial_tally + if (ride_along) daily_est$charter_n else 0)
    if (isTRUE(params$estimate_red_rock)) strat$mean_daily_rr[.ns] <- mean(daily_est$est_rr_exp)
    cat(sprintf("  NOTE: no sampled day in the %s stratum; its %d day(s) take the pooled sampled-day mean.\n",
                paste(strat$day_type[.ns], collapse = "/"), sum(strat$n_total_days[.ns])))
  }
  strat <- strat |>
    mutate(n_calendar_days = n_total_days,
           n_expand_days   = if (census_expansion == "none") n_sampled_days else n_total_days,
           n_imputed_days  = n_expand_days - n_sampled_days,
           est_total_dung    = mean_daily_dung * n_expand_days,       # the expanded part (comm, + charter when riding along)
           est_comm_vessels  = mean_daily_comm_v * n_expand_days,
           est_char_vessels  = mean_daily_char_v * n_expand_days,
           est_total_vessels = mean_daily_vessels * n_expand_days)

  commercial_vessels <- sum(strat$est_comm_vessels)
  commercial_dung    <- commercial_vessels * md_comm
  commercial_observed <- sum(strat$mean_daily_comm_v * strat$n_sampled_days) * md_comm

  cat(sprintf("\n  Commercial vessels by day type (census_expansion = '%s'%s):\n", census_expansion,
              if (census_expansion == "none") "; unsampled days = no operation, so Total = Sampled" else ""))
  cat(sprintf("    %-10s  Sampled  Calendar  Expanded-to  Vessels/day  Vessels\n", "Day Type"))
  for (i in seq_len(nrow(strat)))
    cat(sprintf("    %-10s  %5d    %6d    %8d    %10.1f   %7.0f\n", strat$day_type[i], strat$n_sampled_days[i],
                strat$n_calendar_days[i], strat$n_expand_days[i], strat$mean_daily_comm_v[i], strat$est_comm_vessels[i]))

  # =========================================================================
  # 3. THE CHARTER EXPANSION, per stratum
  # =========================================================================
  iv <- if (n_char && "boat_name" %in% names(ci_char)) {
    ci_char |>
      mutate(vessel = ifelse(is.na(boat_name) | !nzchar(str_squish(as.character(boat_name))), "(unnamed)", str_squish(as.character(boat_name)))) |>
      group_by(vessel) |>
      summarise(n_v = n(), m_v = mean(dungeness_kept), s_v = if (n() > 1) stats::sd(dungeness_kept) else NA_real_,
                obs_v = sum(dungeness_kept), .groups = "drop")
  } else tibble(vessel = character(), n_v = integer(), m_v = numeric(), s_v = numeric(), obs_v = numeric())
  use_vessel <- charter_expansion_req == "vessel" && use_roster && nrow(iv) > 0

  # the strata: one row per rostered vessel (plus the unattributed surplus) when the
  # expansion is by vessel, otherwise ONE pooled stratum -- a pooled expansion has one
  # sample of size n_char, so splitting N across vessel rows and giving each the pooled n
  # would apply the finite-population correction stratum by stratum to a sample that was
  # never stratified, and understate the variance.
  fr_v <- if (!use_roster) tibble(vessel = NA_character_, N_v = sum(strat$est_char_vessels))
    else if (charter_expansion_req == "pooled") tibble(vessel = NA_character_, N_v = sum(frame_day$trips))
    else frame_day |> group_by(vessel) |> summarise(N_v = sum(trips), .groups = "drop")
  charter_detail <- fr_v |>
    left_join(iv, by = "vessel") |>
    mutate(n_v = replace_na(as.numeric(n_v), 0),
           mean_source = case_when(!use_vessel                ~ "pooled charter mean",
                                   is.na(vessel)              ~ "pooled charter mean (not on the roster)",
                                   n_v == 0                   ~ "pooled charter mean (vessel never interviewed)",
                                   n_v == 1                   ~ "vessel mean (1 interview; pooled SD)",
                                   TRUE                       ~ "vessel mean"),
           mean_per_trip = ifelse(grepl("^pooled", mean_source), md_char, m_v),
           s_use = ifelse(grepl("^pooled", mean_source) | is.na(s_v), sd_char, s_v),
           n_use = ifelse(grepl("^pooled", mean_source), n_char, n_v),
           obs_v = replace_na(obs_v, 0),
           est   = N_v * mean_per_trip,
           var   = ifelse(is.finite(s_use) & n_use > 0 & N_v > 0,
                          N_v^2 * pmax(0, 1 - n_use / pmax(N_v, n_use)) * s_use^2 / n_use, 0)) |>
    arrange(desc(N_v))
  charter_trips <- sum(charter_detail$N_v)
  charter_dung  <- sum(charter_detail$est)
  charter_var   <- sum(charter_detail$var)

  cat(sprintf("\n  Charter expansion (%s frame, %s; the vessels are NOT fully sampled): %.0f trip(s), %d interviewed (%.0f%%), %s crab observed -> %s crab; SE %s (%.1f%%)\n",
              charter_frame, if (use_vessel) "stratified by vessel" else "pooled mean", charter_trips, n_char,
              100 * n_char / max(charter_trips, 1), format(round(charter_obs), big.mark = ","),
              format(round(charter_dung), big.mark = ","), format(round(sqrt(charter_var)), big.mark = ","),
              100 * sqrt(charter_var) / max(charter_dung, 1)))
  for (i in seq_len(nrow(charter_detail)))
    cat(sprintf("    %-24s trips %6.1f  interviewed %3.0f  %6.1f crab/trip  est %7.0f  SE %5.0f  [%s]\n",
                ifelse(is.na(charter_detail$vessel[i]), "(not on the roster)", charter_detail$vessel[i]),
                charter_detail$N_v[i], charter_detail$n_v[i], charter_detail$mean_per_trip[i],
                charter_detail$est[i], sqrt(charter_detail$var[i]), charter_detail$mean_source[i]))

  total_dung    <- commercial_dung + charter_dung
  total_vessels <- commercial_vessels + charter_trips

  # =========================================================================
  # 4. VARIANCE: three sources, reported separately, carried per census_uncertainty
  # =========================================================================
  s2_pooled <- if (nrow(daily_est) > 1) stats::var(daily_est$est_dung_exp) else 0
  variance_detail <- strat |>
    select(day_type, n_total_days, n_calendar_days, n_expand_days, n_sampled_days, n_imputed_days, mean_daily_dung) |>
    left_join(daily_est |> group_by(day_type) |>
                summarise(s2_day = if (n() > 1) stats::var(est_dung_exp) else NA_real_, .groups = "drop"), by = "day_type") |>
    mutate(s2_source = ifelse(is.na(s2_day) & n_sampled_days > 0, "pooled (1 sampled day)",
                              ifelse(n_sampled_days == 0, "none (no sampled day)", "stratum")),
           s2_day    = ifelse(is.na(s2_day), s2_pooled, s2_day),
           n_unsampled_days = n_calendar_days - n_sampled_days,
           imputed_dung  = mean_daily_dung * n_imputed_days,
           observed_dung = mean_daily_dung * n_sampled_days,
           var_imputed   = ifelse(n_sampled_days > 0, n_imputed_days^2 * s2_day / n_sampled_days,
                                  n_imputed_days^2 * s2_pooled / max(nrow(daily_est), 1)))
  imputation_var <- sum(variance_detail$var_imputed)
  # the commercial per-vessel MEAN: a sample mean over a near-census of vessel-trips (FPC)
  commercial_var <- 0
  if (n_comm > 1) {
    fpc <- if (is.finite(commercial_vessels) && commercial_vessels > n_comm) (1 - n_comm / commercial_vessels) else 0
    commercial_var <- commercial_vessels^2 * stats::var(ci_comm$dungeness_kept) / n_comm * fpc
  }
  var_total   <- imputation_var + commercial_var + charter_var
  carried_var <- switch(census_mode, none = 0, charter = charter_var, sampling = var_total)
  imputed_dung <- sum(variance_detail$imputed_dung)
  # the floor a census draw can never fall below: the commercial census on the tally days
  # plus the charter crab actually observed on the interviewed trips
  observed_dung <- commercial_observed + charter_obs

  # ---- the per-day table over the whole calendar ----------------------------------
  char_day <- frame_day |>
    left_join(charter_detail |> filter(!is.na(vessel)) |> select(vessel, mean_per_trip), by = "vessel") |>
    mutate(mean_per_trip = replace_na(mean_per_trip, md_char)) |>
    group_by(date) |>
    summarise(charter_trips = sum(trips), est_dung_char = sum(trips * mean_per_trip), .groups = "drop")
  if (ride_along) char_day <- char_day |> mutate(est_dung_char = 0)     # it is inside the expanded commercial figure
  daily_full <- census_calendar |>
    select(date, day_type) |>
    left_join(daily_est |> select(date, commercial_tally, charter_tally, est_dung_exp), by = "date") |>
    left_join(char_day, by = "date") |>
    left_join(strat |> select(day_type, mean_daily_dung, mean_daily_vessels), by = "day_type") |>
    mutate(observed      = !is.na(est_dung_exp),
           charter_trips = replace_na(charter_trips, 0),
           est_dung_char = replace_na(est_dung_char, 0),
           est_dung_exp  = case_when(observed ~ est_dung_exp, census_expansion == "none" ~ 0, TRUE ~ mean_daily_dung),
           est_dung      = est_dung_exp + est_dung_char,
           est_vessels   = ifelse(observed, replace_na(commercial_tally, 0), if (census_expansion == "none") 0 else mean_daily_vessels) +
             if (ride_along) 0 else charter_trips,
           source = case_when(observed & charter_trips > 0 & use_roster ~ "tally (commercial census) + charter roster",
                              observed                                 ~ "tally (commercial census)",
                              charter_trips > 0 & census_expansion == "none" ~ "charter roster (no tally)",
                              charter_trips > 0                        ~ "charter roster (no tally) + imputed commercial",
                              census_expansion == "none"               ~ "no operation (unsampled day)",
                              TRUE                                     ~ "imputed commercial (day-type mean)")) |>
    select(date, day_type, observed, source, commercial_tally, charter_tally, charter_trips, est_vessels, est_dung_char, est_dung)

  cat(sprintf("\n  COMMERCIAL (census, exact on the tally days): %s crab, %.0f vessel-trips%s. Per-vessel mean over %d interviews: SE %s (%.1f%%), %s.\n",
              format(round(commercial_dung), big.mark = ","), commercial_vessels,
              if (census_expansion == "none") sprintf("; %d unsampled calendar day(s) = no operation", sum(variance_detail$n_unsampled_days))
              else sprintf("; %d unsampled day(s) imputed at the day-type mean (%s crab)", sum(variance_detail$n_unsampled_days), format(round(imputed_dung), big.mark = ",")),
              n_comm, format(round(sqrt(commercial_var)), big.mark = ","), 100 * sqrt(commercial_var) / max(commercial_dung, 1),
              if (census_mode == "sampling") "carried into the port interval" else "reported, not carried"))
  cat(sprintf("  CHARTER (expansion): %s crab. Expansion SE %s (%.1f%%), %s.\n",
              format(round(charter_dung), big.mark = ","), format(round(sqrt(charter_var)), big.mark = ","),
              100 * sqrt(charter_var) / max(charter_dung, 1),
              if (census_mode %in% c("charter", "sampling")) "carried into the port interval" else "reported, not carried"))
  if (any(variance_detail$s2_source != "stratum"))
    cat("    NOTE: ", paste(sprintf("%s stratum: %s", variance_detail$day_type[variance_detail$s2_source != "stratum"],
                                    variance_detail$s2_source[variance_detail$s2_source != "stratum"]), collapse = "; "), "\n")

  vessel_type_detail <- tibble(
    vessel_type = c("Commercial", "Charter"),
    estimator = c(sprintf("census (census_expansion = %s)", census_expansion),
                  sprintf("expansion (%s frame, %s)", charter_frame, if (use_vessel) "by vessel" else "pooled mean")),
    n_interviews = c(n_comm, n_char),
    mean_dung_per_vessel = c(md_comm, md_char),
    tally_vessels = c(commercial_vessels, charter_trips),
    sampled_frac = c(if (commercial_vessels > 0) n_comm / commercial_vessels else NA_real_,
                     if (charter_trips > 0) n_char / charter_trips else NA_real_),
    sampled_catch = c(commercial_dung, charter_dung),
    se = c(sqrt(commercial_var), sqrt(charter_var)))

  daily_est_out <- daily_full |>
    filter(est_dung > 0 | observed) |>
    transmute(date, day_of_week = weekdays(date), day_type, commercial_tally, charter_tally,
              charter_n = charter_trips, total_comm_charter = replace_na(commercial_tally, 0) + charter_trips, est_dung)

  result <- list(
    Dungeness_Kept = total_dung,
    effort_total = total_vessels,
    daily_est = daily_est_out,
    strat_detail = strat,
    vessel_type_detail = vessel_type_detail,
    observed_dung = observed_dung,
    imputed_dung = imputed_dung,
    Dungeness_Kept_var = var_total,
    Dungeness_Kept_se = sqrt(var_total),
    carried_var = carried_var,
    carried_se = sqrt(carried_var),
    imputation_var = imputation_var,
    census_uncertainty = census_mode,
    census_expansion = census_expansion,
    n_unsampled_days = sum(variance_detail$n_unsampled_days),
    daily_full = daily_full,
    variance_detail = variance_detail,
    # 2026-09-11: the two components, separately
    commercial_dung = commercial_dung, commercial_vessels = commercial_vessels,
    commercial_var = commercial_var, commercial_se = sqrt(commercial_var),
    commercial_interviews = n_comm,
    charter_dung = charter_dung, charter_observed_dung = charter_obs, charter_trips = charter_trips,
    charter_interviews = n_char, charter_sampled_frac = if (charter_trips > 0) n_char / charter_trips else NA_real_,
    charter_var = charter_var, charter_se = sqrt(charter_var),
    charter_frame = charter_frame, charter_expansion = if (use_vessel) "vessel" else "pooled",
    charter_detail = charter_detail, roster_reconciliation = roster_recon,
    charter_roster_dung = charter_roster_dung, n_roster_only_days = n_roster_only_days, n_roster_trips = n_roster_trips
  )
  if (isTRUE(params$estimate_red_rock))
    result$Red_Rock_Kept <- sum(strat$mean_daily_rr * strat$n_expand_days) + if (ride_along) 0 else charter_trips * mr_char

  cat(sprintf("\n  Est Dungeness (commercial census + charter expansion): %s%s | total SE %s; %s carried under census_uncertainty = '%s'\n",
              format(round(total_dung), big.mark = ","),
              if (isTRUE(params$estimate_red_rock)) sprintf(", Red Rock: %s", format(round(result$Red_Rock_Kept), big.mark = ",")) else "",
              format(round(sqrt(var_total)), big.mark = ","), format(round(sqrt(carried_var)), big.mark = ","), census_mode))

  return(result)
}
