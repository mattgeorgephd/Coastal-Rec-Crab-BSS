###############################################################################
# Crab Creel Estimation - Westport Docks Proof of Concept
# 2024-25 Season (Sep 16, 2024 - Sep 15, 2025)
#
# Adapted from WDFW CreelEstimates repository (FW creel)
#
# BUG FIXES from v1:
#   - Fixed interview CSV: number_of_gear column was empty (export bug)
#   - Fixed Stan O dimension: now correctly [D,S,G] (was [D,S])
#   - Added output folder structure: output/YYYYMMDD/
###############################################################################

# ===========================================================================
# 0. SETUP & PARAMETERS
# ===========================================================================

## Load Packages
load.lib<-c("tidyverse", "lubridate", "suncalc", "gt", "patchwork", "rstan", "here")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

rstan_options(auto_write = TRUE)

# Source the existing FW creel functions
purrr::walk(list.files(here("R_functions"), full.names = TRUE), source)

# --- Create output folder ---
run_date <- format(Sys.Date(), "%Y%m%d")
output_dir <- here("output", run_date)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat(sprintf("Output directory: %s\n", output_dir))

# --- Analysis parameters ---
params <- list(
  project_name         = "Coastal Recreational Crab",
  fishery_name         = "Rec Crab Grays Harbor Westport Docks 2024-25",
  est_date_start       = "2024-09-16",
  est_date_end         = "2025-09-15",
  season_filter        = "2024-25",
  est_catch_groups     = data.frame(
    species = c("Dungeness", "Red_Rock"),
    fate    = c("Kept", "Kept"),
    stringsAsFactors = FALSE
  ),
  study_design         = "Crab",
  days_wkend           = c("Friday", "Saturday", "Sunday"),
  min_fishing_time     = 0.5,
  period_pe            = "month",
  period_bss           = "week",
  sections             = c(1),
  bss_model_file_name  = "BSS_creel_model_02_2024-07-24.stan",
  model_used           = "Both models",
  data_grade           = "provisional",
  export               = "local"
)

# Crabbing holidays (HIGH-EFFORT holidays for 2024-25 season)
crabbing_holiday_dates <- as.Date(c(
  "2024-09-02",  # Labor Day
  "2024-11-29",  # Native American Heritage Day
  "2024-12-31",  # New Year's Eve
  "2025-01-01",  # New Year's Day
  "2025-02-08",  # Super Bowl Eve
  "2025-05-26",  # Memorial Day
  "2025-06-15"   # Father's Day
))

# ===========================================================================
# 1. DATA INGESTION
# ===========================================================================

fetch_crab_data <- function(params) {

  cat("Reading effort data...\n")
  effort_raw <- read_csv(
    here("input_files", "effort_combined.csv"),
    show_col_types = FALSE
  ) |>
    filter(season == params$season_filter) |>
    mutate(date = as.Date(date))

  cat("Reading interview data...\n")
  interview_raw <- read_csv(
    here("input_files", "interview_combined.csv"),
    show_col_types = FALSE,
    col_types = cols(
      date = col_date(), crabbers = col_double(),
      number_of_gear = col_double(), dungeness_kept = col_double(),
      red_rock_kept = col_double(), hours_fished = col_double(),
      crabber_hours = col_double(), gear_hours = col_double(),
      completed_trip = col_character(), total_vehicles = col_double(),
      crabbing_holiday = col_integer()
    )
  ) |>
    filter(season == params$season_filter)

  dock_areas <- c("Westport Docks Float 20", "Westport Docks Float 17-21")
  effort_docks <- effort_raw |> filter(creel_area %in% dock_areas)
  interview_docks <- interview_raw |> filter(creel_area %in% dock_areas)

  cat(sprintf("  Raw effort records at Westport Docks: %d\n", nrow(effort_docks)))
  cat(sprintf("  Raw interview records at Westport Docks: %d\n", nrow(interview_docks)))

  # --- EFFORT: assign count_sequence, pair Float 20 + 17-21 ---
  effort_parsed <- effort_docks |>
    mutate(
      event_date = date,
      count_time_posix = as.POSIXct(paste(date, count_time),
                                     format = "%Y-%m-%d %H:%M:%S",
                                     tz = "America/Los_Angeles")
    ) |>
    filter(!is.na(count_time_posix))

  f20 <- effort_parsed |>
    filter(creel_area == "Westport Docks Float 20") |>
    arrange(event_date, count_time_posix) |>
    group_by(event_date) |>
    mutate(count_sequence = row_number()) |>
    ungroup()

  f17 <- effort_parsed |>
    filter(creel_area == "Westport Docks Float 17-21")

  if (nrow(f17) > 0) {
    f17_paired <- f17 |>
      left_join(
        f20 |> select(event_date, count_sequence, f20_time = count_time_posix),
        by = "event_date", relationship = "many-to-many"
      ) |>
      mutate(time_diff = abs(as.numeric(difftime(count_time_posix, f20_time, units = "mins")))) |>
      group_by(event_date, survey_id, count_time) |>
      slice_min(time_diff, n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(event_date, count_sequence, f17_gear = total_gear_count)
  } else {
    f17_paired <- tibble(event_date = Date(), count_sequence = integer(), f17_gear = numeric())
  }

  effort_index <- f20 |>
    select(event_date, count_sequence, f20_gear = total_gear_count) |>
    left_join(f17_paired, by = c("event_date", "count_sequence")) |>
    mutate(
      f17_gear = replace_na(f17_gear, 0),
      count_quantity = f20_gear + f17_gear,
      section_num = 1, count_type = "Gear Count",
      water_body = "Grays Harbor",
      project_name = params$project_name,
      fishery_name = params$fishery_name,
      tie_in_indicator = 0L,
      no_count_reason = NA_character_,
      angler_final = "dock", angler_final_int = 1L
    )

  cat(sprintf("  Effort index (section totals): %d obs\n", nrow(effort_index)))
  cat(sprintf("  Days with effort: %d\n", n_distinct(effort_index$event_date)))
  cat(sprintf("  Max count_sequence (H): %d\n", max(effort_index$count_sequence)))
  effort_index |> count(event_date) |> count(n, name = "n_days") |>
    mutate(label = sprintf("    %d seq/day: %d days", n, n_days)) |>
    pull(label) |> walk(cat, "\n")

  # --- INTERVIEWS ---
  interview <- interview_docks |>
    filter(crabbing_mode == "Dock", !is.na(crabbers), as.numeric(crabbers) > 0) |>
    mutate(
      event_date = as.Date(date),
      water_body = "Grays Harbor",
      project_name = params$project_name,
      fishery_name = params$fishery_name,
      interview_id = paste0(survey_id, "_", interview_num),
      section_num = 1,
      angler_type = "dock", angler_final = "dock", angler_final_int = 1L,
      angler_count = as.integer(crabbers),
      person_count_final = as.integer(crabbers),
      number_of_gear = as.numeric(number_of_gear),
      hours_fished = as.numeric(hours_fished),
      crabber_hours_calc = as.numeric(crabber_hours),
      fishing_time_total = case_when(
        !is.na(crabber_hours_calc) & crabber_hours_calc > 0 ~ crabber_hours_calc,
        !is.na(hours_fished) & !is.na(angler_count) ~ hours_fished * angler_count,
        TRUE ~ NA_real_
      ),
      gear_hours_calc = as.numeric(gear_hours),
      gear_time_total = case_when(
        !is.na(gear_hours_calc) & gear_hours_calc > 0 ~ gear_hours_calc,
        !is.na(hours_fished) & !is.na(number_of_gear) ~ hours_fished * number_of_gear,
        TRUE ~ NA_real_
      ),
      dungeness_kept = replace_na(as.numeric(dungeness_kept), 0),
      red_rock_kept = replace_na(as.numeric(red_rock_kept), 0),
      trip_status = case_when(
        completed_trip == "1" ~ "Complete",
        completed_trip == "0" ~ "Incomplete",
        TRUE ~ NA_character_
      ),
      vehicle_count = as.integer(total_vehicles),
      trailer_count = 0L,
      crabbing_holiday = as.integer(crabbing_holiday)
    ) |>
    filter(!is.na(fishing_time_total), fishing_time_total >= params$min_fishing_time)

  cat(sprintf("  Interviews (filtered): %d\n", nrow(interview)))

  # --- CATCH ---
  catch <- bind_rows(
    interview |> filter(dungeness_kept > 0) |>
      transmute(interview_id, event_date,
                catch_id = paste0(interview_id, "_dung"),
                species = "Dungeness", fate = "Kept",
                fish_count = as.integer(dungeness_kept),
                catch_group = "Dungeness_Kept"),
    interview |> filter(red_rock_kept > 0) |>
      transmute(interview_id, event_date,
                catch_id = paste0(interview_id, "_rr"),
                species = "Red_Rock", fate = "Kept",
                fish_count = as.integer(red_rock_kept),
                catch_group = "Red_Rock_Kept")
  )

  return(list(
    effort = effort_index, interview = interview, catch = catch,
    ll = tibble(water_body_desc="Grays Harbor", centroid_lat=46.904, centroid_lon=-124.105),
    closures = tibble(fishery_name=character(), section_num=numeric(), event_date=Date()),
    fishery_manager = tibble(project_name=params$project_name, fishery_name=params$fishery_name,
                             section_num=1, p_census_bank=1.0, p_census_boat=1.0)
  ))
}

# ===========================================================================
# 2. PREP DAYS
# ===========================================================================

prep_days_crab <- function(date_begin, date_end, weekends, holiday_dates,
                           lat, long, period_pe, sections, ...) {
  date_begin <- as.Date(date_begin); date_end <- as.Date(date_end)
  days <- tibble(
    event_date = seq.Date(date_begin, date_end, by = "day"),
    day = weekdays(event_date),
    day_type = case_when(
      event_date %in% holiday_dates ~ "holiday",
      day %in% weekends             ~ "weekend",
      TRUE                          ~ "weekday"
    ),
    day_type_num_weekend = as.integer(day_type %in% c("weekend", "holiday")),
    day_type_num_holiday = as.integer(day_type == "holiday"),
    week = as.numeric(format(event_date, "%W")),
    month = as.numeric(format(event_date, "%m")),
    year = as.numeric(format(event_date, "%Y")),
    period = case_when(
      period_pe == "week" ~ as.numeric(format(event_date, "%W")),
      period_pe == "month" ~ as.numeric(format(event_date, "%m")),
      period_pe == "duration" ~ 0
    ),
    day_index = as.integer(seq_along(event_date)),
    week_index = as.integer(factor(
      paste(year, sprintf("%02d", week)),
      levels = unique(paste(year, sprintf("%02d", week)))
    )),
    month_index = as.integer(factor(
      paste(year, sprintf("%02d", month)),
      levels = unique(paste(year, sprintf("%02d", month)))
    )),
    day_length = if_else(month %in% 4:9, 10.0, 8.5)
  )
  for (s in sections) days[[paste0("open_section_", s)]] <- TRUE
  days |> mutate(fishery_name = params$fishery_name) |> relocate(fishery_name)
}

# ===========================================================================
# 3. DATA SUMMARIZATION
# ===========================================================================

prep_dwg_crab_summary <- function(dwg, params, days) {
  dwg_summ <- list()

  int_filtered <- dwg$interview |>
    filter(between(event_date, as.Date(params$est_date_start), as.Date(params$est_date_end)))

  catch_wide <- dwg$catch |>
    group_by(interview_id, catch_group) |>
    summarise(fish_count = sum(fish_count), .groups = "drop") |>
    pivot_wider(names_from = catch_group, values_from = fish_count, values_fill = 0)

  dwg_summ$interview <- int_filtered |>
    left_join(catch_wide, by = "interview_id") |>
    mutate(across(any_of(c("Dungeness_Kept", "Red_Rock_Kept")), ~replace_na(., 0)))

  dwg_summ$effort_index <- dwg$effort |>
    filter(between(event_date, as.Date(params$est_date_start), as.Date(params$est_date_end)))

  dwg_summ$effort_census <- dwg_summ$effort_index |>
    mutate(count_census = count_quantity, tie_in_indicator = 1L)

  dwg_summ$census_expan <- tibble(angler_final = "dock", section_num = 1, p_census = 1.0)

  # Crabbers-per-gear ratio
  ratio_data <- dwg_summ$interview |>
    filter(!is.na(number_of_gear), number_of_gear > 0, angler_count > 0)
  dwg_summ$crabbers_per_gear <- sum(ratio_data$angler_count) / sum(ratio_data$number_of_gear)

  int <- dwg_summ$interview; eff <- dwg_summ$effort_index
  cat("\n", strrep("=", 60), "\n DATA SUMMARY: Westport Docks 2024-25\n", strrep("=", 60), "\n")
  cat(sprintf("Days in season: %d\n", nrow(days)))
  cat(sprintf("Days sampled (effort): %d (%.0f%%)\n",
              n_distinct(eff$event_date), 100*n_distinct(eff$event_date)/nrow(days)))
  cat(sprintf("Days sampled (interviews): %d\n", n_distinct(int$event_date)))
  cat(sprintf("Total effort obs: %d | Total interviews: %d\n", nrow(eff), nrow(int)))
  cat(sprintf("Completed trips: %d (%.0f%%)\n",
              sum(int$trip_status == "Complete", na.rm=TRUE),
              100*mean(int$trip_status == "Complete", na.rm=TRUE)))
  cat(sprintf("Crabbers-per-gear ratio: %.2f\n", dwg_summ$crabbers_per_gear))

  if ("Dungeness_Kept" %in% names(int)) {
    cat(sprintf("Observed Dungeness kept: %d\n", sum(int$Dungeness_Kept)))
    comp <- int |> filter(trip_status == "Complete", fishing_time_total > 0)
    if (nrow(comp) > 0)
      cat(sprintf("Mean CPUE (completed): %.3f crab/crabber-hr\n",
                  sum(comp$Dungeness_Kept)/sum(comp$fishing_time_total)))
  }
  days |> count(day_type) |> print()

  return(dwg_summ)
}

# ===========================================================================
# 4. POINT ESTIMATE (PE) METHOD
# ===========================================================================

run_pe_crab <- function(dwg_summ, days, params) {
  cat("\n", strrep("=", 60), "\n PE METHOD\n", strrep("=", 60), "\n")
  results <- list()

  daily_effort <- dwg_summ$effort_index |>
    group_by(event_date, section_num) |>
    summarise(mean_gear_count = mean(count_quantity), n_counts = n(), .groups = "drop") |>
    mutate(est_crabbers = mean_gear_count * dwg_summ$crabbers_per_gear) |>
    left_join(days |> select(event_date, day_type, day_length, period, month_index), by = "event_date") |>
    mutate(
      est_daily_crabber_hours = est_crabbers * day_length,
      est_daily_gear_hours = mean_gear_count * day_length
    )
  results$daily_effort <- daily_effort

  total_days_strat <- days |> filter(open_section_1) |>
    group_by(period, day_type) |> summarise(n_total_days = n(), .groups = "drop")

  effort_strat <- daily_effort |>
    group_by(section_num, period, day_type) |>
    summarise(
      mean_daily_effort = mean(est_daily_crabber_hours, na.rm=TRUE),
      sd_daily_effort = sd(est_daily_crabber_hours, na.rm=TRUE),
      mean_daily_gear_effort = mean(est_daily_gear_hours, na.rm=TRUE),
      n_sampled = n(), .groups = "drop"
    ) |>
    left_join(total_days_strat, by = c("period", "day_type")) |>
    mutate(
      est_total_crabber_hours = mean_daily_effort * n_total_days,
      se_total_crabber_hours = sqrt((n_total_days^2) * replace_na(sd_daily_effort^2,0) / pmax(n_sampled,1)),
      est_total_gear_hours = mean_daily_gear_effort * n_total_days
    )
  results$effort_by_stratum <- effort_strat

  total_eff <- sum(effort_strat$est_total_crabber_hours, na.rm=TRUE)
  total_se <- sqrt(sum(effort_strat$se_total_crabber_hours^2, na.rm=TRUE))
  cat(sprintf("\nTotal crabber-hours: %s (SE: %s)\n",
              format(round(total_eff), big.mark=","), format(round(total_se), big.mark=",")))
  cat(sprintf("Total gear-hours: %s\n",
              format(round(sum(effort_strat$est_total_gear_hours, na.rm=TRUE)), big.mark=",")))

  for (cg in c("Dungeness_Kept", "Red_Rock_Kept")) {
    if (!cg %in% names(dwg_summ$interview)) next

    daily_cpue <- dwg_summ$interview |>
      group_by(event_date, section_num) |>
      summarise(
        total_catch = sum(.data[[cg]], na.rm=TRUE),
        total_crabber_hours = sum(fishing_time_total, na.rm=TRUE),
        total_gear_hours = sum(gear_time_total, na.rm=TRUE),
        n_int = n(), .groups = "drop"
      ) |>
      mutate(
        cpue_ch = if_else(total_crabber_hours > 0, total_catch/total_crabber_hours, 0),
        cpue_gh = if_else(total_gear_hours > 0, total_catch/total_gear_hours, 0)
      ) |>
      left_join(days |> select(event_date, day_type, period), by = "event_date")

    cpue_strat <- daily_cpue |>
      group_by(section_num, period, day_type) |>
      summarise(
        mean_cpue_ch = weighted.mean(cpue_ch, w=n_int, na.rm=TRUE),
        mean_cpue_gh = weighted.mean(cpue_gh, w=n_int, na.rm=TRUE),
        .groups = "drop"
      )

    catch_strat <- effort_strat |>
      left_join(cpue_strat, by = c("section_num", "period", "day_type")) |>
      mutate(
        est_catch_ch = est_total_crabber_hours * replace_na(mean_cpue_ch, 0),
        est_catch_gh = est_total_gear_hours * replace_na(mean_cpue_gh, 0)
      )

    cat(sprintf("\n%s:\n  Crabber-hr CPUE: %s crab\n  Gear-hr CPUE:    %s crab\n", cg,
                format(round(sum(catch_strat$est_catch_ch, na.rm=TRUE)), big.mark=","),
                format(round(sum(catch_strat$est_catch_gh, na.rm=TRUE)), big.mark=",")))
    results[[paste0("catch_", cg)]] <- catch_strat
  }
  return(results)
}

# ===========================================================================
# 5. BSS MODEL PREPARATION
# ===========================================================================
# FIX: O is now array(dim = c(D, S, G)) to match Stan declaration real O[D,S,G]

prep_inputs_bss_crab <- function(days, dwg_summ, est_catch_group, period, params) {
  cat(sprintf("\nPreparing BSS inputs: %s, period=%s\n", est_catch_group, period))

  eff <- dwg_summ$effort_index
  cen <- dwg_summ$effort_census
  int_cg <- dwg_summ$interview |>
    mutate(fish_count = .data[[est_catch_group]]) |>
    filter(!is.na(fishing_time_total), fishing_time_total > 0)

  D <- nrow(days); G <- 1L; S <- 1L; H <- max(eff$count_sequence)
  P_n <- case_when(
    tolower(period)=="day"~max(days$day_index), tolower(period)=="week"~max(days$week_index),
    tolower(period)=="month"~max(days$month_index), tolower(period)=="duration"~1L)
  pvec <- case_when(
    tolower(period)=="day"~days$day_index, tolower(period)=="week"~days$week_index,
    tolower(period)=="month"~days$month_index, tolower(period)=="duration"~rep(1L,D))

  eff_d <- eff |> left_join(days |> select(event_date,day_index), by="event_date") |> filter(!is.na(day_index))
  cen_d <- cen |> left_join(days |> select(event_date,day_index), by="event_date") |> filter(!is.na(day_index))
  int_d <- int_cg |> left_join(days |> select(event_date,day_index), by="event_date") |>
    filter(!is.na(day_index)) |> distinct(interview_id, .keep_all=TRUE)
  intA_d <- int_d |> filter(!is.na(number_of_gear), number_of_gear > 0, angler_count > 0)

  stan_list <- list(
    est_cg=est_catch_group, D=D, G=G, S=S, H=H, P_n=P_n, period=pvec,
    w = days$day_type_num_weekend, L = days$day_length,

    # FIX: O must be 3D array [D, S, G] to match Stan declaration
    O = array(1.0, dim = c(D, S, G)),

    V_n=0L, day_V=integer(0), section_V=integer(0), countnum_V=integer(0), V_I=integer(0),
    T_n=0L, day_T=integer(0), section_T=integer(0), countnum_T=integer(0), T_I=integer(0),
    A_n=nrow(eff_d), day_A=eff_d$day_index, gear_A=rep(1L,nrow(eff_d)),
    section_A=rep(1L,nrow(eff_d)), countnum_A=as.integer(eff_d$count_sequence),
    A_I=as.integer(eff_d$count_quantity),
    B_n=0L, day_B=integer(0), gear_B=integer(0), section_B=integer(0),
    countnum_B=integer(0), B_s=integer(0),
    E_n=nrow(cen_d), day_E=cen_d$day_index, gear_E=rep(1L,nrow(cen_d)),
    section_E=rep(1L,nrow(cen_d)), countnum_E=as.integer(cen_d$count_sequence),
    E_s=as.integer(cen_d$count_census),
    p_TI = matrix(1.0, nrow=G, ncol=S),
    IntC=nrow(int_d), day_IntC=int_d$day_index, gear_IntC=rep(1L,nrow(int_d)),
    section_IntC=rep(1L,nrow(int_d)), c=as.integer(int_d$fish_count), h=int_d$fishing_time_total,
    IntA=nrow(intA_d), day_IntA=intA_d$day_index, gear_IntA=rep(1L,nrow(intA_d)),
    section_IntA=rep(1L,nrow(intA_d)), A_A=as.integer(intA_d$angler_count),
    V_A=rep(0L,nrow(intA_d)), T_A=rep(0L,nrow(intA_d)), B_A=rep(0L,nrow(intA_d)),
    value_cauchyDF_sigma_eps_E=5, value_cauchyDF_sigma_r_E=5,
    value_cauchyDF_sigma_eps_C=5, value_cauchyDF_sigma_r_C=5,
    value_betashape_phi_E_scaled=2, value_betashape_phi_C_scaled=2,
    value_normal_sigma_B1=1,
    value_normal_sigma_omega_C_0=3, value_normal_sigma_omega_E_0=3,
    value_lognormal_sigma_b=0.5,
    value_normal_mu_mu_C=log(0.5), value_normal_sigma_mu_C=2,
    value_normal_mu_mu_E=log(25), value_normal_sigma_mu_E=2,
    value_cauchyDF_sigma_mu_C=5, value_cauchyDF_sigma_mu_E=5
  )

  cat(sprintf("  D=%d G=%d S=%d H=%d P_n=%d\n", D, G, S, H, P_n))
  cat(sprintf("  Gear index: %d | Census: %d | CPUE int: %d | Expan int: %d\n",
              stan_list$A_n, stan_list$E_n, stan_list$IntC, stan_list$IntA))
  return(stan_list)
}

# ===========================================================================
# 6. RUN THE ANALYSIS
# ===========================================================================

cat("\n", strrep("=", 60), "\n LOADING DATA\n", strrep("=", 60), "\n")
dwg <- fetch_crab_data(params)

days <- prep_days_crab(
  date_begin = params$est_date_start, date_end = params$est_date_end,
  weekends = params$days_wkend, holiday_dates = crabbing_holiday_dates,
  lat = dwg$ll$centroid_lat, long = dwg$ll$centroid_lon,
  period_pe = params$period_pe, sections = params$sections
)

dwg_summ <- prep_dwg_crab_summary(dwg, params, days)

# --- Diagnostic plots ---

p_effort <- dwg_summ$effort_index |>
  group_by(event_date) |>
  summarise(total_gear = sum(count_quantity), n_counts = n(), .groups = "drop") |>
  left_join(days |> select(event_date, day_type), by = "event_date") |>
  ggplot(aes(x = event_date, y = total_gear, color = day_type, size = n_counts)) +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(1,3), breaks = 1:3, name = "Counts/day") +
  labs(title = "Westport Docks: Daily Gear Count (section total)",
       subtitle = "Point size = count sequences per day (3-count protocol from ~Mar 2025)",
       x = "Date", y = "Total gear count", color = "Day type") +
  theme_bw()
ggsave(file.path(output_dir, "plot_effort_timeseries.png"), p_effort, width = 10, height = 6)
print(p_effort)

p_cpue <- dwg_summ$interview |>
  filter(fishing_time_total > 0, !is.na(Dungeness_Kept)) |>
  mutate(cpue = Dungeness_Kept / fishing_time_total) |>
  group_by(event_date) |>
  summarise(mean_cpue = mean(cpue), n_int = n(), .groups = "drop") |>
  left_join(days |> select(event_date, day_type), by = "event_date") |>
  ggplot(aes(x = event_date, y = mean_cpue, color = day_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, color = "grey30") +
  labs(title = "Westport Docks: Dungeness CPUE", x = "Date",
       y = "Mean CPUE (crab/crabber-hr)", color = "Day type") +
  theme_bw()
ggsave(file.path(output_dir, "plot_cpue_timeseries.png"), p_cpue, width = 10, height = 6)
print(p_cpue)

# Effort by month and day type
p_effort_month <- dwg_summ$effort_index |>
  group_by(event_date) |>
  summarise(total_gear = sum(count_quantity), .groups = "drop") |>
  left_join(days |> select(event_date, day_type, month, year), by = "event_date") |>
  mutate(month_label = factor(format(event_date, "%Y-%m"))) |>
  ggplot(aes(x = month_label, y = total_gear, fill = day_type)) +
  geom_boxplot() +
  labs(title = "Gear Counts by Month and Day Type", x = "Month", y = "Daily gear count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "plot_effort_by_month.png"), p_effort_month, width = 12, height = 6)
print(p_effort_month)

# --- Run PE ---
pe_results <- run_pe_crab(dwg_summ, days, params)

# Save PE effort plot
p_pe_effort <- pe_results$daily_effort |>
  ggplot(aes(x = event_date, y = est_daily_crabber_hours, color = day_type)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, color = "grey30") +
  labs(title = "PE: Estimated Daily Crabber-Hours at Westport Docks",
       x = "Date", y = "Estimated crabber-hours", color = "Day type") +
  theme_bw()
ggsave(file.path(output_dir, "plot_pe_daily_effort.png"), p_pe_effort, width = 10, height = 6)
print(p_pe_effort)

# --- Prepare BSS ---
bss_inputs_dung <- prep_inputs_bss_crab(days, dwg_summ, "Dungeness_Kept", params$period_bss, params)
bss_inputs_rr <- prep_inputs_bss_crab(days, dwg_summ, "Red_Rock_Kept", params$period_bss, params)

# --- Run BSS (Dungeness) ---
cat("\n", strrep("=", 60), "\n RUNNING BSS MODEL (Dungeness)\n", strrep("=", 60), "\n")
cat("This may take 10-30 minutes...\n")

fit_dung <- rstan::stan(
  file = here("stan_models", params$bss_model_file_name),
  data = bss_inputs_dung,
  chains = 4, iter = 2000, warmup = 1000, thin = 1,
  control = list(adapt_delta = 0.9, max_treedepth = 12),
  cores = 6
)

# --- Extract BSS results ---
bss_summary <- summary(fit_dung, pars = c("C_sum", "E_sum"))$summary
cat("\n--- BSS Summary (Dungeness) ---\n")
print(bss_summary)

# Save Stan fit summary
write.csv(as.data.frame(bss_summary), file.path(output_dir, "bss_summary_dungeness.csv"))

# Daily effort posteriors
E_draws <- rstan::extract(fit_dung, "E")$E  # [draws, S, D, G]
E_daily <- tibble(
  event_date = days$event_date,
  median = apply(E_draws[,1,,1], 2, median),
  lo95 = apply(E_draws[,1,,1], 2, quantile, 0.025),
  hi95 = apply(E_draws[,1,,1], 2, quantile, 0.975)
)

p_bss_effort <- E_daily |>
  left_join(days |> select(event_date, day_type), by = "event_date") |>
  ggplot(aes(x = event_date)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = median), color = "steelblue") +
  labs(title = "BSS: Estimated Daily Effort (Crabber-Hours)",
       subtitle = "Median with 95% credible interval",
       x = "Date", y = "Crabber-hours") +
  theme_bw()
ggsave(file.path(output_dir, "plot_bss_daily_effort.png"), p_bss_effort, width = 10, height = 6)
print(p_bss_effort)

# Daily catch posteriors
C_draws <- rstan::extract(fit_dung, "C")$C  # [draws, S, D, G]
C_daily <- tibble(
  event_date = days$event_date,
  median = apply(C_draws[,1,,1], 2, median),
  lo95 = apply(C_draws[,1,,1], 2, quantile, 0.025),
  hi95 = apply(C_draws[,1,,1], 2, quantile, 0.975)
)

p_bss_catch <- C_daily |>
  ggplot(aes(x = event_date)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, fill = "darkgreen") +
  geom_line(aes(y = median), color = "darkgreen") +
  labs(title = "BSS: Estimated Daily Dungeness Catch",
       subtitle = "Median with 95% credible interval",
       x = "Date", y = "Dungeness crab caught") +
  theme_bw()
ggsave(file.path(output_dir, "plot_bss_daily_catch_dungeness.png"), p_bss_catch, width = 10, height = 6)
print(p_bss_catch)

# Season totals posterior
C_sum_draws <- rstan::extract(fit_dung, "C_sum")$C_sum
E_sum_draws <- rstan::extract(fit_dung, "E_sum")$E_sum

p_bss_totals <- tibble(
  C_sum = C_sum_draws,
  E_sum = E_sum_draws
) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~param, scales = "free", labeller = labeller(
    param = c(C_sum = "Season Total Catch (Dungeness)", E_sum = "Season Total Effort (crabber-hrs)")
  )) +
  labs(title = "BSS: Posterior Distributions of Season Totals", x = "Value", y = "Density") +
  theme_bw()
ggsave(file.path(output_dir, "plot_bss_posteriors.png"), p_bss_totals, width = 12, height = 5)
print(p_bss_totals)

# ===========================================================================
# 7. COMBINED RESULTS SUMMARY
# ===========================================================================

cat("\n", strrep("=", 60), "\n")
cat(" FINAL RESULTS: Westport Docks 2024-25\n")
cat(strrep("=", 60), "\n\n")

# PE results
pe_effort <- sum(pe_results$effort_by_stratum$est_total_crabber_hours, na.rm=TRUE)
pe_dung_ch <- sum(pe_results$catch_Dungeness_Kept$est_catch_ch, na.rm=TRUE)
pe_dung_gh <- sum(pe_results$catch_Dungeness_Kept$est_catch_gh, na.rm=TRUE)
pe_rr_ch <- sum(pe_results$catch_Red_Rock_Kept$est_catch_ch, na.rm=TRUE)
pe_rr_gh <- sum(pe_results$catch_Red_Rock_Kept$est_catch_gh, na.rm=TRUE)

# BSS results
bss_E <- quantile(E_sum_draws, c(0.025, 0.5, 0.975))
bss_C <- quantile(C_sum_draws, c(0.025, 0.5, 0.975))

results_table <- tibble(
  Estimate = c("Total Effort (crabber-hrs)",
               "Dungeness Kept (crabber-hr CPUE)",
               "Dungeness Kept (gear-hr CPUE)",
               "Red Rock Kept (crabber-hr CPUE)",
               "Red Rock Kept (gear-hr CPUE)"),
  PE = c(round(pe_effort), round(pe_dung_ch), round(pe_dung_gh),
         round(pe_rr_ch), round(pe_rr_gh)),
  `BSS Median` = c(round(bss_E[2]), round(bss_C[2]), NA, NA, NA),
  `BSS 2.5%` = c(round(bss_E[1]), round(bss_C[1]), NA, NA, NA),
  `BSS 97.5%` = c(round(bss_E[3]), round(bss_C[3]), NA, NA, NA)
)

print(results_table)
write.csv(results_table, file.path(output_dir, "results_summary.csv"), row.names = FALSE)

# Save all PE stratum-level results
write.csv(pe_results$effort_by_stratum, file.path(output_dir, "pe_effort_by_stratum.csv"), row.names = FALSE)
write.csv(pe_results$catch_Dungeness_Kept, file.path(output_dir, "pe_catch_dungeness_by_stratum.csv"), row.names = FALSE)
write.csv(pe_results$catch_Red_Rock_Kept, file.path(output_dir, "pe_catch_redrock_by_stratum.csv"), row.names = FALSE)
write.csv(pe_results$daily_effort, file.path(output_dir, "pe_daily_effort.csv"), row.names = FALSE)

# Save BSS daily estimates
write.csv(E_daily, file.path(output_dir, "bss_daily_effort.csv"), row.names = FALSE)
write.csv(C_daily, file.path(output_dir, "bss_daily_catch_dungeness.csv"), row.names = FALSE)

# Save run parameters
writeLines(
  capture.output(str(params)),
  file.path(output_dir, "run_parameters.txt")
)

cat(sprintf("\nAll outputs saved to: %s\n", output_dir))
cat("Files:\n")
cat(paste(" ", list.files(output_dir), collapse = "\n"), "\n")
