# ===========================================================================
# GEAR-HOURS REFORMULATION — Required R code changes
# 
# These diffs apply to BSS-GH-pooled-CPUE-model.Rmd
# Shore model is UNCHANGED. Only the boat (private_boat) code path changes.
# ===========================================================================


# ---------------------------------------------------------------------------
# CHANGE 1: prep_bss_crab() — Build Stan data list
#
# In the boat code path, three things change:
#   (a) h = gear_time_total instead of fishing_time_total
#   (b) L = 24 instead of day_length
#   (c) IntA data passes Gear_A_boat (gear count) instead of T_A_int/A_A_trailer
# ---------------------------------------------------------------------------

# CURRENT (line ~570 in Rmd):
#   h=int_d$fishing_time_total,

# NEW:
#   h = if(is_shore) int_d$fishing_time_total else int_d$gear_time_total,

# CURRENT (line ~534):
#   L=days$day_length,

# NEW:
#   L = if(is_shore) days$day_length else rep(24, D),

# CURRENT (lines ~575-578):
#   IntA_trailer = if(!is_shore) nrow(intA) else 0L,
#   T_A_int = if(!is_shore) rep(1L, nrow(intA)) else integer(0),
#   A_A_trailer = if(!is_shore) as.integer(intA$angler_count) else integer(0),

# NEW:
#   IntA_trailer = if(!is_shore) nrow(intA) else 0L,
#   Gear_A_boat = if(!is_shore) as.integer(intA$number_of_gear) else integer(0),

# CURRENT prior (line ~590):
#   value_normal_mu_mu_E=if(is_shore) log(25) else log(10),

# NEW (lambda_E is now gear count for boats; ~20 gear typical):
#   value_normal_mu_mu_E = if(is_shore) log(25) else log(20),

# CURRENT CPUE prior (line ~588):
#   value_normal_mu_mu_C=log(0.5), value_normal_sigma_mu_C=2,

# NEW (CPUE is now crab/gear-hr for boats, ~0.1; keep log(0.5) but could lower):
#   value_normal_mu_mu_C = if(is_shore) log(0.5) else log(0.1),


# ---------------------------------------------------------------------------
# CHANGE 2: prep_population_summary() — Compute gear_time_total for interviews
#
# The interview processing (lines ~228-250) already computes fishing_time_total
# and gear_time_total. Verify gear_time_total is retained on the interview tibble.
# ---------------------------------------------------------------------------

# EXISTING CODE (already present, just confirm it's retained):
#   gear_time_total = case_when(
#     !is.na(as.numeric(gear_hours)) & as.numeric(gear_hours) > 0 ~ as.numeric(gear_hours),
#     !is.na(hours_fished) & !is.na(number_of_gear) ~ hours_fished * number_of_gear,
#     TRUE ~ NA_real_
#   ),

# ADD a filter for boat interviews: gear_time_total must be valid
# After the existing filter line:
#   filter(!is.na(fishing_time_total), fishing_time_total >= params$min_fishing_time)
# Add for boats:
#   |> filter(if(population_name == "private_boat") !is.na(gear_time_total) & gear_time_total > 0 else TRUE)


# ---------------------------------------------------------------------------
# CHANGE 3: run_pe() — PE method for boats
#
# The PE for boats should also use gear-hours for consistency.
# ---------------------------------------------------------------------------

# CURRENT daily effort for boats (line ~438):
#   est_daily_crabber_hours = est_crabbers * day_length,

# For boats, this should be:
#   est_daily_gear_hours = mean_count * gear_per_trailer * 24

# Where gear_per_trailer = sum(number_of_gear) / n_groups from interviews.
# This replaces crabbers_per_gear with gear_per_trailer for boats.

# CURRENT CPUE denominator (line ~464):
#   hrs=sum(fishing_time_total, na.rm=TRUE)

# For boats:
#   hrs=sum(gear_time_total, na.rm=TRUE)


# Full replacement for run_pe() boat path:
run_pe_boat_gearhours <- function(summ, days, params, population_name = "private_boat") {
  results <- list()
  
  # Gear per trailer group (from interviews)
  ratio_data <- summ$interview |> 
    filter(!is.na(number_of_gear), number_of_gear > 0, angler_count > 0)
  gear_per_group <- if(nrow(ratio_data) > 0) {
    mean(ratio_data$number_of_gear)
  } else 4.0  # default
  
  # Daily effort: trailer count → gear in water → gear-hours (24 hrs/day)
  daily_effort <- summ$effort_index |>
    filter(count_sequence <= params$bss_max_count_seq) |>
    group_by(event_date, section_num) |>
    summarise(mean_count = mean(count_quantity), n_counts = n(), .groups = "drop") |>
    mutate(est_gear = mean_count * gear_per_group,
           est_daily_gear_hours = est_gear * 24) |>
    left_join(days |> select(event_date, day_type, day_length, period), by = "event_date")
  
  # Total days per stratum
  total_days_strat <- days |> filter(open_section_1) |>
    group_by(period, day_type) |> summarise(n_total_days = n(), .groups = "drop")
  
  # Stratum-level effort totals (in gear-hours)
  effort_strat <- daily_effort |>
    group_by(section_num, period, day_type) |>
    summarise(mean_daily = mean(est_daily_gear_hours, na.rm = TRUE),
              sd_daily = sd(est_daily_gear_hours, na.rm = TRUE),
              n_sampled = n(), .groups = "drop") |>
    left_join(total_days_strat, by = c("period", "day_type")) |>
    mutate(est_total = mean_daily * n_total_days,
           se_total = sqrt((n_total_days^2) * replace_na(sd_daily^2, 0) / pmax(n_sampled, 1)))
  
  results$effort_total <- sum(effort_strat$est_total, na.rm = TRUE)
  results$effort_se <- sqrt(sum(effort_strat$se_total^2, na.rm = TRUE))
  
  # CPUE using gear-hours denominator
  for (cg in catch_groups) {
    if (!cg %in% names(summ$interview)) { results[[cg]] <- 0; next }
    
    daily_cpue <- summ$interview |>
      filter(!is.na(gear_time_total), gear_time_total > 0) |>
      group_by(event_date, section_num) |>
      summarise(catch = sum(.data[[cg]], na.rm = TRUE),
                gear_hrs = sum(gear_time_total, na.rm = TRUE),
                n_int = n(), .groups = "drop") |>
      mutate(cpue = if_else(gear_hrs > 0, catch / gear_hrs, 0)) |>
      left_join(days |> select(event_date, day_type, period), by = "event_date")
    
    cpue_strat <- daily_cpue |>
      group_by(section_num, period, day_type) |>
      summarise(mean_cpue = weighted.mean(cpue, w = n_int, na.rm = TRUE), .groups = "drop")
    
    catch_strat <- effort_strat |>
      left_join(cpue_strat, by = c("section_num", "period", "day_type")) |>
      mutate(est_catch = est_total * replace_na(mean_cpue, 0))
    
    results[[cg]] <- sum(catch_strat$est_catch, na.rm = TRUE)
  }
  
  cat(sprintf("  PE %s (gear-hours): Effort=%s gear-hrs, Dung=%s\n", population_name,
              format(round(results$effort_total), big.mark = ","),
              format(round(results$Dungeness_Kept), big.mark = ",")))
  
  return(results)
}


# ---------------------------------------------------------------------------
# CHANGE 4: Output interpretation
#
# After this change:
#   - Shore E_sum = crabber-hours (unchanged)
#   - Boat E_sum = gear-hours (new unit)
#   - Shore CPUE = crab per crabber-hour
#   - Boat CPUE = crab per gear-hour
#
# For the port total, you have two options:
#   (a) Report shore and boat effort in their native units (recommended)
#   (b) Convert boat gear-hours to crabber-hours for a common unit:
#       boat_crabber_hours = boat_gear_hours / R_G_boat * crabbers_per_group
#       (but this introduces another ratio and defeats the purpose)
#
# Catch is always in crab, so C_sum is directly comparable across components.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# CHANGE 5: Stan model file reference
#
# Update the model file path:
#   CURRENT: bss_model_file = "crab_bss_pooled.stan"
#   NEW:     bss_model_file = "crab_bss_pooled_gearhours.stan"
# ---------------------------------------------------------------------------
