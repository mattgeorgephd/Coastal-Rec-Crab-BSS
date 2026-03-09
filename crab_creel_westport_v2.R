###############################################################################
# Crab Creel Estimation - Westport / Grays Harbor
# 2024-25 Season | Framework v2
#
# Three population components estimated independently:
#   1. Shore crabbers (dock + jetty + beach) — BSS + PE
#   2. Private boat crabbers — BSS + PE
#   3. Commercial/charter boats — Census/tally
#
# Two gear-regime sub-seasons:
#   Sub-season 1: Ring-net only (Sep 16 – Nov 30, no pots)
#   Sub-season 2: All-gear (Dec 1 – Sep 15, pots allowed)
#
# Port total = sum of all components
###############################################################################

# ===========================================================================
# 0. SETUP
# ===========================================================================

load.lib <- c("tidyverse","lubridate","suncalc","gt","patchwork","rstan","here")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)
rstan_options(auto_write = TRUE)
purrr::walk(list.files(here("R_functions"), full.names = TRUE), source)

run_date <- format(Sys.Date(), "%Y%m%d")
output_dir <- here("output", run_date)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===========================================================================
# 0.1 PARAMETERS
# ===========================================================================

params <- list(
  project_name      = "Coastal Recreational Crab",
  fishery_name      = "Rec Crab Grays Harbor Westport 2024-25",
  est_date_start    = "2024-09-16",
  est_date_end      = "2025-09-15",
  season_filter     = "2024-25",

  # Gear restriction: pots prohibited Sep 16 – Nov 30
  pot_open_date     = "2024-12-01",

  # Commercial opener: commercial vessels stop rec fishing after this date
  commercial_opener = "2025-01-15",

  # Catch groups toggle: set to FALSE to skip Red Rock (saves ~50% BSS runtime)
  estimate_red_rock = FALSE,

  # Day type
  days_wkend        = c("Friday", "Saturday", "Sunday"),
  min_fishing_time  = 0.5,

  # Time strata
  period_pe         = "month",
  period_bss        = "month",
  sections          = c(1),

  # BSS settings
  bss_model_file    = "BSS_crab_model_01.stan",
  bss_chains        = 4,
  bss_iter          = 2000,
  bss_warmup        = 1000,
  bss_adapt_delta   = 0.9,
  bss_max_treedepth = 10,
  bss_cores         = 4,
  bss_max_interviews = 1000,  # set to NULL for full dataset
  bss_max_count_seq  = 3
)

# Crabbing holidays
crabbing_holiday_dates <- as.Date(c(
  "2024-09-02","2024-11-29","2024-12-31",
  "2025-01-01","2025-02-08","2025-05-26","2025-06-15"
))

# Derived: catch groups vector used throughout
catch_groups <- if(params$estimate_red_rock) {
  c("Dungeness_Kept", "Red_Rock_Kept")
} else {
  "Dungeness_Kept"
}

# Derived sub-season dates
subseason_1 <- list(
  name = "ring_net_only",
  start = as.Date(params$est_date_start),
  end   = as.Date(params$pot_open_date) - 1
)
subseason_2 <- list(
  name = "all_gear",
  start = as.Date(params$pot_open_date),
  end   = as.Date(params$est_date_end)
)

cat(sprintf("Sub-season 1 (ring-net): %s to %s (%d days)\n",
            subseason_1$start, subseason_1$end,
            as.integer(subseason_1$end - subseason_1$start + 1)))
cat(sprintf("Sub-season 2 (all-gear): %s to %s (%d days)\n",
            subseason_2$start, subseason_2$end,
            as.integer(subseason_2$end - subseason_2$start + 1)))
cat(sprintf("Catch groups: %s\n", paste(catch_groups, collapse = ", ")))

# ===========================================================================
# 1. DATA INGESTION — Classify by population
# ===========================================================================

fetch_crab_data_v2 <- function(params) {

  cat("Reading data...\n")
  effort_raw <- read_csv(here("input_files","effort_combined.csv"), show_col_types=FALSE) |>
    filter(season == params$season_filter) |> mutate(date = as.Date(date))

  interview_raw <- read_csv(here("input_files","interview_combined.csv"), show_col_types=FALSE,
    col_types = cols(date=col_date(), crabbers=col_double(), number_of_gear=col_double(),
      dungeness_kept=col_double(), red_rock_kept=col_double(), hours_fished=col_double(),
      crabber_hours=col_double(), gear_hours=col_double(), completed_trip=col_character(),
      total_vehicles=col_double(), crabbing_holiday=col_integer())) |>
    filter(season == params$season_filter)

  # Filter to Grays Harbor
  gh_effort <- effort_raw |> filter(creel_area %in% c(
    "Westport Docks Float 20","Westport Docks Float 17-21",
    "Westport Boat Launch","Westport Marina","Westport Jetty",
    "Ocean Shores Boat Launch","Damon Point"))

  gh_interview <- interview_raw |> filter(creel_location == "Grays Harbor")

  # --- CLASSIFY INTERVIEWS BY POPULATION ---
  gh_interview <- gh_interview |>
    mutate(
      event_date = as.Date(date),
      boat_type_clean = case_when(
        str_detect(boat_type, "(?i)commer") ~ "Commercial",
        str_detect(boat_type, "(?i)charter") ~ "Charter",
        str_detect(boat_type, "(?i)guide") ~ "Charter",
        str_detect(boat_type, "(?i)private") ~ "Private",
        TRUE ~ NA_character_
      ),
      population = case_when(
        boat_type_clean %in% c("Commercial","Charter") ~ "comm_charter",
        crabbing_mode == "Boat" & (is.na(boat_type_clean) | boat_type_clean == "Private") ~ "private_boat",
        crabbing_mode %in% c("Dock","Jetty","Beach") ~ "shore",
        TRUE ~ "shore"
      )
    )

  cat("  Grays Harbor interviews by population:\n")
  gh_interview |> count(population) |> mutate(l=sprintf("    %s: %d",population,n)) |>
    pull(l) |> walk(cat,"\n")

  # --- PROCESS INTERVIEWS (common transformations) ---
  gh_interview <- gh_interview |>
    filter(!is.na(crabbers), as.numeric(crabbers) > 0) |>
    mutate(
      interview_id = paste0(survey_id,"_",interview_num),
      section_num = 1,
      angler_count = as.integer(crabbers),
      number_of_gear = as.numeric(number_of_gear),
      hours_fished = as.numeric(hours_fished),
      crabber_hours_calc = as.numeric(crabber_hours),
      fishing_time_total = case_when(
        !is.na(crabber_hours_calc) & crabber_hours_calc > 0 ~ crabber_hours_calc,
        !is.na(hours_fished) & !is.na(angler_count) ~ hours_fished * angler_count,
        TRUE ~ NA_real_
      ),
      gear_time_total = case_when(
        !is.na(as.numeric(gear_hours)) & as.numeric(gear_hours) > 0 ~ as.numeric(gear_hours),
        !is.na(hours_fished) & !is.na(number_of_gear) ~ hours_fished * number_of_gear,
        TRUE ~ NA_real_
      ),
      dungeness_kept = replace_na(as.numeric(dungeness_kept), 0),
      red_rock_kept = replace_na(as.numeric(red_rock_kept), 0),
      trip_status = case_when(
        completed_trip=="1"~"Complete", completed_trip=="0"~"Incomplete", TRUE~NA_character_
      ),
      angler_final_int = case_when(population=="shore"~1L, population=="private_boat"~2L, TRUE~NA_integer_)
    ) |>
    filter(!is.na(fishing_time_total), fishing_time_total >= params$min_fishing_time)

  # --- SHORE EFFORT: Gear counts at docks (pair F20 + F17-21) ---
  dock_effort <- gh_effort |>
    filter(creel_area %in% c("Westport Docks Float 20","Westport Docks Float 17-21")) |>
    mutate(event_date = date,
           count_time_posix = as.POSIXct(paste(date,count_time), format="%Y-%m-%d %H:%M:%S", tz="America/Los_Angeles")) |>
    filter(!is.na(count_time_posix))

  f20 <- dock_effort |> filter(creel_area=="Westport Docks Float 20") |>
    arrange(event_date,count_time_posix) |> group_by(event_date) |>
    mutate(count_sequence=row_number()) |> ungroup()

  f17 <- dock_effort |> filter(creel_area=="Westport Docks Float 17-21")

  if(nrow(f17) > 0) {
    f17_paired <- f17 |>
      left_join(f20 |> select(event_date,count_sequence,f20_time=count_time_posix),
                by="event_date", relationship="many-to-many") |>
      mutate(time_diff=abs(as.numeric(difftime(count_time_posix,f20_time,units="mins")))) |>
      group_by(event_date,survey_id,count_time) |> slice_min(time_diff,n=1,with_ties=FALSE) |> ungroup() |>
      select(event_date,count_sequence,f17_gear=total_gear_count)
  } else {
    f17_paired <- tibble(event_date=Date(),count_sequence=integer(),f17_gear=numeric())
  }

  shore_effort <- f20 |> select(event_date,count_sequence,f20_gear=total_gear_count) |>
    left_join(f17_paired, by=c("event_date","count_sequence")) |>
    mutate(f17_gear=replace_na(f17_gear,0), count_quantity=f20_gear+f17_gear,
           section_num=1, count_type="Gear Count", population="shore")

  # --- BOAT EFFORT: Trailer counts at boat launches ---
  boat_effort <- gh_effort |>
    filter(creel_area %in% c("Westport Boat Launch","Ocean Shores Boat Launch")) |>
    mutate(event_date = date,
           count_time_posix = as.POSIXct(paste(date,count_time), format="%Y-%m-%d %H:%M:%S", tz="America/Los_Angeles"),
           count_quantity = as.numeric(boat_trailer_count)) |>
    filter(!is.na(count_time_posix), !is.na(count_quantity), count_quantity >= 0) |>
    arrange(event_date,count_time_posix) |> group_by(event_date) |>
    mutate(count_sequence=row_number()) |> ungroup() |>
    select(event_date,count_sequence,count_quantity) |>
    mutate(section_num=1, count_type="Trailer Count", population="private_boat")

  # --- COMMERCIAL TALLY ---
  comm_tally <- read_csv(here("input_files","wes_commercial_tally.csv"), show_col_types=FALSE) |>
    mutate(date = as.Date(date))

  cat(sprintf("\n  Shore effort obs: %d (%d days)\n", nrow(shore_effort), n_distinct(shore_effort$event_date)))
  cat(sprintf("  Boat effort obs: %d (%d days)\n", nrow(boat_effort), n_distinct(boat_effort$event_date)))
  cat(sprintf("  Commercial tally days: %d\n", nrow(comm_tally)))

  # --- CATCH tables ---
  catch <- bind_rows(
    gh_interview |> filter(dungeness_kept>0) |>
      transmute(interview_id,event_date,population,species="Dungeness",fate="Kept",
                fish_count=as.integer(dungeness_kept),catch_group="Dungeness_Kept"),
    if(params$estimate_red_rock) {
      gh_interview |> filter(red_rock_kept>0) |>
        transmute(interview_id,event_date,population,species="Red_Rock",fate="Kept",
                  fish_count=as.integer(red_rock_kept),catch_group="Red_Rock_Kept")
    } else { tibble() }
  )

  return(list(
    shore_effort = shore_effort,
    boat_effort = boat_effort,
    interview = gh_interview,
    catch = catch,
    comm_tally = comm_tally,
    ll = tibble(centroid_lat=46.904, centroid_lon=-124.105)
  ))
}

# ===========================================================================
# 2. PREP DAYS (works for any date range)
# ===========================================================================

prep_days_crab <- function(date_begin, date_end, weekends, holiday_dates, period_pe, sections) {
  date_begin <- as.Date(date_begin); date_end <- as.Date(date_end)
  days <- tibble(
    event_date = seq.Date(date_begin, date_end, by="day"),
    day = weekdays(event_date),
    day_type = case_when(
      event_date %in% holiday_dates ~ "holiday",
      day %in% weekends ~ "weekend", TRUE ~ "weekday"),
    day_type_num_weekend = as.integer(day_type %in% c("weekend","holiday")),
    week = as.numeric(format(event_date,"%W")),
    month = as.numeric(format(event_date,"%m")),
    year = as.numeric(format(event_date,"%Y")),
    period = if(period_pe=="month") as.numeric(format(event_date,"%m"))
             else as.numeric(format(event_date,"%W")),
    day_index = as.integer(seq_along(event_date)),
    month_index = as.integer(factor(paste(year,sprintf("%02d",month)),
                  levels=unique(paste(year,sprintf("%02d",month))))),
    day_length = if_else(month %in% 4:9, 10.0, 8.5)
  )
  for(s in sections) days[[paste0("open_section_",s)]] <- TRUE
  days
}

# ===========================================================================
# 3. SUMMARIZE DATA FOR A GIVEN POPULATION × SUB-SEASON
# ===========================================================================

prep_population_summary <- function(dwg, population_name, date_start, date_end, params) {
  ds <- as.Date(date_start); de <- as.Date(date_end)
  summ <- list()

  if(population_name == "shore") {
    summ$effort_index <- dwg$shore_effort |> filter(between(event_date, ds, de))
  } else if(population_name == "private_boat") {
    summ$effort_index <- dwg$boat_effort |> filter(between(event_date, ds, de))
  }

  int_pop <- dwg$interview |>
    filter(population == population_name, between(event_date, ds, de))

  catch_wide <- dwg$catch |>
    filter(population == population_name) |>
    group_by(interview_id, catch_group) |>
    summarise(fish_count=sum(fish_count), .groups="drop") |>
    pivot_wider(names_from=catch_group, values_from=fish_count, values_fill=0)

  summ$interview <- int_pop |>
    left_join(catch_wide, by="interview_id") |>
    mutate(across(any_of(catch_groups), ~replace_na(.,0)))

  ratio_data <- summ$interview |> filter(!is.na(number_of_gear), number_of_gear>0, angler_count>0)
  summ$crabbers_per_gear <- if(nrow(ratio_data)>0) sum(ratio_data$angler_count)/sum(ratio_data$number_of_gear) else 1.0

  cat(sprintf("\n  %s [%s to %s]: %d effort obs, %d interviews, crab/gear=%.2f\n",
              population_name, ds, de, nrow(summ$effort_index), nrow(summ$interview), summ$crabbers_per_gear))

  return(summ)
}

# ===========================================================================
# 4. PE METHOD (works for any population)
# ===========================================================================

run_pe <- function(summ, days, params, population_name) {
  results <- list()

  daily_effort <- summ$effort_index |>
    filter(count_sequence <= params$bss_max_count_seq) |>
    group_by(event_date, section_num) |>
    summarise(mean_count=mean(count_quantity), n_counts=n(), .groups="drop") |>
    mutate(est_crabbers = mean_count * summ$crabbers_per_gear) |>
    left_join(days |> select(event_date,day_type,day_length,period), by="event_date") |>
    mutate(est_daily_crabber_hours = est_crabbers * day_length,
           est_daily_gear_hours = mean_count * day_length)

  total_days_strat <- days |> filter(open_section_1) |>
    group_by(period, day_type) |> summarise(n_total_days=n(), .groups="drop")

  effort_strat <- daily_effort |>
    group_by(section_num, period, day_type) |>
    summarise(mean_daily=mean(est_daily_crabber_hours,na.rm=TRUE),
              sd_daily=sd(est_daily_crabber_hours,na.rm=TRUE),
              n_sampled=n(), .groups="drop") |>
    left_join(total_days_strat, by=c("period","day_type")) |>
    mutate(est_total = mean_daily * n_total_days,
           se_total = sqrt((n_total_days^2)*replace_na(sd_daily^2,0)/pmax(n_sampled,1)))

  results$effort_total <- sum(effort_strat$est_total, na.rm=TRUE)
  results$effort_se <- sqrt(sum(effort_strat$se_total^2, na.rm=TRUE))

  for(cg in catch_groups) {
    if(!cg %in% names(summ$interview)) { results[[cg]] <- 0; next }

    daily_cpue <- summ$interview |>
      group_by(event_date, section_num) |>
      summarise(catch=sum(.data[[cg]],na.rm=TRUE), hrs=sum(fishing_time_total,na.rm=TRUE),
                n_int=n(), .groups="drop") |>
      mutate(cpue=if_else(hrs>0, catch/hrs, 0)) |>
      left_join(days |> select(event_date,day_type,period), by="event_date")

    cpue_strat <- daily_cpue |>
      group_by(section_num, period, day_type) |>
      summarise(mean_cpue=weighted.mean(cpue,w=n_int,na.rm=TRUE), .groups="drop")

    catch_strat <- effort_strat |>
      left_join(cpue_strat, by=c("section_num","period","day_type")) |>
      mutate(est_catch = est_total * replace_na(mean_cpue, 0))

    results[[cg]] <- sum(catch_strat$est_catch, na.rm=TRUE)
  }

  rr_str <- if(params$estimate_red_rock) sprintf(", RR=%s", format(round(results$Red_Rock_Kept),big.mark=",")) else ""
  cat(sprintf("  PE %s: Effort=%s hrs, Dung=%s%s\n", population_name,
              format(round(results$effort_total),big.mark=","),
              format(round(results$Dungeness_Kept),big.mark=","), rr_str))

  return(results)
}

# ===========================================================================
# 5. BSS PREP (for crab-specific Stan model)
# ===========================================================================

prep_bss_crab <- function(days, summ, est_catch_group, params, population_name) {

  eff <- summ$effort_index |> filter(count_sequence <= params$bss_max_count_seq)
  D <- nrow(days); G <- 1L; S <- 1L; H <- max(c(eff$count_sequence, 1L))
  P_n <- max(days$month_index)
  pvec <- days$month_index

  eff_d <- eff |> left_join(days |> select(event_date,day_index), by="event_date") |> filter(!is.na(day_index))

  int_cg <- summ$interview |>
    mutate(fish_count = .data[[est_catch_group]]) |>
    filter(!is.na(fishing_time_total), fishing_time_total > 0)

  int_d <- int_cg |> left_join(days |> select(event_date,day_index), by="event_date") |>
    filter(!is.na(day_index)) |> distinct(interview_id, .keep_all=TRUE)

  intA <- int_d |> filter(!is.na(number_of_gear), number_of_gear>0, angler_count>0)

  # Subsample
  if(!is.null(params$bss_max_interviews) && nrow(int_d)>params$bss_max_interviews) {
    set.seed(42)
    int_d <- int_d |> slice_sample(n=params$bss_max_interviews)
    intA <- intA |> filter(interview_id %in% int_d$interview_id)
  }

  is_shore <- (population_name == "shore")

  list(
    D=D, G=G, S=S, H=H, P_n=P_n, period=pvec,
    w=days$day_type_num_weekend, L=days$day_length,
    O=array(1.0, dim=c(D,S,G)),

    Gear_n = if(is_shore) nrow(eff_d) else 0L,
    day_Gear = if(is_shore) eff_d$day_index else integer(0),
    section_Gear = if(is_shore) rep(1L,nrow(eff_d)) else integer(0),
    countnum_Gear = if(is_shore) as.integer(eff_d$count_sequence) else integer(0),
    Gear_I = if(is_shore) as.integer(eff_d$count_quantity) else integer(0),

    T_n = if(!is_shore) nrow(eff_d) else 0L,
    day_T = if(!is_shore) eff_d$day_index else integer(0),
    section_T = if(!is_shore) rep(1L,nrow(eff_d)) else integer(0),
    countnum_T = if(!is_shore) as.integer(eff_d$count_sequence) else integer(0),
    T_I = if(!is_shore) as.integer(eff_d$count_quantity) else integer(0),

    Crab_n=0L, day_Crab=integer(0), section_Crab=integer(0),
    countnum_Crab=integer(0), Crab_I=integer(0), p_I_crab=1.0,

    IntC=nrow(int_d), day_IntC=int_d$day_index, gear_IntC=rep(1L,nrow(int_d)),
    section_IntC=rep(1L,nrow(int_d)), c=as.integer(int_d$fish_count), h=int_d$fishing_time_total,

    IntA_gear = if(is_shore) nrow(intA) else 0L,
    Gear_A = if(is_shore) as.integer(intA$number_of_gear) else integer(0),
    A_A_gear = if(is_shore) as.integer(intA$angler_count) else integer(0),

    IntA_trailer = if(!is_shore) nrow(intA) else 0L,
    T_A_int = if(!is_shore) rep(1L, nrow(intA)) else integer(0),
    A_A_trailer = if(!is_shore) as.integer(intA$angler_count) else integer(0),

    value_cauchyDF_sigma_eps_E=5, value_cauchyDF_sigma_r_E=5,
    value_cauchyDF_sigma_eps_C=5, value_cauchyDF_sigma_r_C=5,
    value_betashape_phi_E_scaled=2, value_betashape_phi_C_scaled=2,
    value_normal_sigma_B1=1,
    value_normal_mu_mu_C=log(0.5), value_normal_sigma_mu_C=2,
    value_normal_mu_mu_E=if(is_shore) log(25) else log(10),
    value_normal_sigma_mu_E=2,
    value_cauchyDF_sigma_mu_C=5, value_cauchyDF_sigma_mu_E=5
  )
}

# ===========================================================================
# 6. COMMERCIAL/CHARTER CENSUS ESTIMATION
# ===========================================================================

estimate_comm_charter <- function(dwg, params) {
  cat("\n--- Commercial/Charter Census Estimation ---\n")

  tally <- dwg$comm_tally |>
    filter(date <= as.Date(params$commercial_opener))

  comm_int <- dwg$interview |>
    filter(population == "comm_charter",
           event_date <= as.Date(params$commercial_opener))

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
      est_dung = total_comm_charter * mean_dung_per_vessel
    )

  if(params$estimate_red_rock) {
    daily_est <- daily_est |> mutate(est_rr = total_comm_charter * mean_rr_per_vessel)
  }

  sampled_days <- nrow(daily_est)
  total_possible_days <- as.integer(
    min(as.Date(params$commercial_opener), as.Date(params$est_date_end)) -
    as.Date(params$est_date_start) + 1
  )
  expansion <- total_possible_days / sampled_days

  total_dung <- sum(daily_est$est_dung) * expansion
  result <- list(
    Dungeness_Kept = total_dung,
    effort_total = sum(daily_est$total_comm_charter) * expansion,
    daily_est = daily_est
  )

  if(params$estimate_red_rock) {
    result$Red_Rock_Kept <- sum(daily_est$est_rr) * expansion
  }

  dung_str <- format(round(total_dung), big.mark=",")
  rr_out <- if(params$estimate_red_rock) sprintf(", Red Rock: %s", format(round(result$Red_Rock_Kept),big.mark=",")) else ""
  cat(sprintf("  Sampled %d of %d possible days (expansion=%.1f)\n",
              sampled_days, total_possible_days, expansion))
  cat(sprintf("  Est Dungeness: %s%s\n", dung_str, rr_out))

  return(result)
}

# ===========================================================================
# 7. RUN EVERYTHING
# ===========================================================================

cat("\n", strrep("=",60), "\n LOADING DATA\n", strrep("=",60), "\n")
dwg <- fetch_crab_data_v2(params)

# --- PE for each population × sub-season ---
cat("\n", strrep("=",60), "\n PE ESTIMATES\n", strrep("=",60), "\n")

pe_all <- list()
for(pop in c("shore","private_boat")) {
  for(ss in list(subseason_1, subseason_2)) {
    label <- paste0(pop, "_", ss$name)
    days_ss <- prep_days_crab(ss$start, ss$end, params$days_wkend, crabbing_holiday_dates,
                              params$period_pe, params$sections)
    summ_ss <- prep_population_summary(dwg, pop, ss$start, ss$end, params)
    pe_all[[label]] <- run_pe(summ_ss, days_ss, params, label)
  }
}

# Commercial/charter
pe_all$comm_charter <- estimate_comm_charter(dwg, params)

# --- PE PORT TOTALS ---
cat("\n", strrep("=",60), "\n PE PORT TOTALS\n", strrep("=",60), "\n")

component_names <- c("Shore (ring-net)","Shore (all-gear)","Private Boat (ring-net)",
                     "Private Boat (all-gear)","Commercial/Charter","PORT TOTAL")
component_keys <- c("shore_ring_net_only","shore_all_gear",
                    "private_boat_ring_net_only","private_boat_all_gear","comm_charter")

pe_summary <- tibble(
  Component = component_names,
  Effort_hrs = c(sapply(component_keys, function(k) pe_all[[k]]$effort_total), NA),
  Dungeness = c(sapply(component_keys, function(k) pe_all[[k]]$Dungeness_Kept), NA)
)
pe_summary$Effort_hrs[6] <- sum(pe_summary$Effort_hrs[1:5], na.rm=TRUE)
pe_summary$Dungeness[6] <- sum(pe_summary$Dungeness[1:5], na.rm=TRUE)

if(params$estimate_red_rock) {
  pe_summary$Red_Rock <- c(
    sapply(component_keys, function(k) pe_all[[k]]$Red_Rock_Kept %||% 0), NA)
  pe_summary$Red_Rock[6] <- sum(pe_summary$Red_Rock[1:5], na.rm=TRUE)
}

print(pe_summary |> mutate(across(where(is.numeric), ~round(.))))
write.csv(pe_summary, file.path(output_dir, "pe_port_summary.csv"), row.names=FALSE)

# --- BSS RUNS ---
cat("\n", strrep("=",60), "\n BSS MODEL RUNS\n", strrep("=",60), "\n")

# Build list of BSS catch groups to run
bss_catch_groups <- catch_groups
cat(sprintf("Running BSS for: %s\n", paste(bss_catch_groups, collapse=", ")))

bss_all <- list()
for(cg in bss_catch_groups) {
  for(pop in c("shore","private_boat")) {
    for(ss in list(subseason_1, subseason_2)) {
      label <- paste0(pop, "_", ss$name, "_", cg)
      cat(sprintf("\n--- BSS: %s ---\n", label))

      days_ss <- prep_days_crab(ss$start, ss$end, params$days_wkend, crabbing_holiday_dates,
                                params$period_pe, params$sections)
      summ_ss <- prep_population_summary(dwg, pop, ss$start, ss$end, params)

      if(nrow(summ_ss$effort_index) < 5 || nrow(summ_ss$interview) < 10) {
        cat("  Insufficient data for BSS, using PE only.\n")
        bss_all[[label]] <- list(C_sum=NA, E_sum=NA, pe_fallback=TRUE,
                                 population=pop, subseason=ss$name, catch_group=cg)
        next
      }

      bss_data <- prep_bss_crab(days_ss, summ_ss, cg, params, pop)

      fit <- rstan::stan(
        file = here("stan_models", params$bss_model_file),
        data = bss_data,
        chains = params$bss_chains, iter = params$bss_iter, warmup = params$bss_warmup,
        control = list(adapt_delta=params$bss_adapt_delta, max_treedepth=params$bss_max_treedepth),
        cores = params$bss_cores
      )

      bss_summ <- summary(fit, pars=c("C_sum","E_sum"))$summary
      cat("  BSS Summary:\n"); print(bss_summ)
      write.csv(as.data.frame(bss_summ), file.path(output_dir, paste0("bss_summary_",label,".csv")))

      C_draws <- rstan::extract(fit, "C_sum")$C_sum
      E_draws <- rstan::extract(fit, "E_sum")$E_sum

      bss_all[[label]] <- list(
        C_sum = quantile(C_draws, c(0.025,0.5,0.975)),
        E_sum = quantile(E_draws, c(0.025,0.5,0.975)),
        C_draws = C_draws, E_draws = E_draws,
        fit = fit, pe_fallback = FALSE,
        population = pop, subseason = ss$name, catch_group = cg
      )
    }
  }
}

# --- COMBINED BSS PORT TOTAL (Dungeness) ---
cat("\n", strrep("=",60), "\n COMBINED RESULTS\n", strrep("=",60), "\n")

n_draws <- params$bss_chains * (params$bss_iter - params$bss_warmup)

for(cg in bss_catch_groups) {
  bss_C_total <- rep(0, n_draws)
  bss_E_total <- rep(0, n_draws)

  for(pop in c("shore","private_boat")) {
    for(ss in list(subseason_1, subseason_2)) {
      label <- paste0(pop, "_", ss$name, "_", cg)
      pe_label <- paste0(pop, "_", ss$name)
      b <- bss_all[[label]]

      if(!is.null(b$pe_fallback) && b$pe_fallback) {
        bss_C_total <- bss_C_total + (pe_all[[pe_label]][[cg]] %||% 0)
        bss_E_total <- bss_E_total + pe_all[[pe_label]]$effort_total
      } else {
        n <- min(n_draws, length(b$C_draws))
        bss_C_total[1:n] <- bss_C_total[1:n] + b$C_draws[1:n]
        bss_E_total[1:n] <- bss_E_total[1:n] + b$E_draws[1:n]
      }
    }
  }
  # Add commercial/charter as constant
  bss_C_total <- bss_C_total + (pe_all$comm_charter[[cg]] %||% 0)
  bss_E_total <- bss_E_total + pe_all$comm_charter$effort_total

  final_table <- tibble(
    Catch_Group = cg,
    Estimate = c("Effort (crabber-hrs)", "Catch"),
    PE = c(round(pe_summary$Effort_hrs[6]),
           round(if(cg=="Dungeness_Kept") pe_summary$Dungeness[6] else pe_summary$Red_Rock[6])),
    BSS_median = c(round(median(bss_E_total)), round(median(bss_C_total))),
    BSS_lo95 = c(round(quantile(bss_E_total,0.025)), round(quantile(bss_C_total,0.025))),
    BSS_hi95 = c(round(quantile(bss_E_total,0.975)), round(quantile(bss_C_total,0.975)))
  )

  cat(sprintf("\n--- Port Total: %s ---\n", cg))
  print(final_table)
  write.csv(final_table, file.path(output_dir, paste0("port_total_",cg,".csv")), row.names=FALSE)
}

# Save parameters
writeLines(capture.output(str(params)), file.path(output_dir, "run_parameters.txt"))

cat(sprintf("\nAll outputs saved to: %s\n", output_dir))
cat(paste(" ", list.files(output_dir), collapse="\n"), "\n")
