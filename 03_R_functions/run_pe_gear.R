###############################################################################
# run_pe_gear.R  (gear-resolved driver)
#
# Point Estimator (PE) for the gear-resolved model: stratified effort and catch
# expansion per population x sub-season, with the P0/P1/P2 fixes (explicit
# population argument, bss_effort_spec-driven effort unit, ratio-of-sums stratum
# CPUE, and a scale-consistency assertion). Extracted from the gear-resolved driver
# and renamed run_pe_gear. catch_groups is derived from params.
###############################################################################

run_pe_gear <- function(summ, days, params, population_name, population = NULL) {
  catch_groups <- if (isTRUE(params$estimate_red_rock)) c("Dungeness_Kept", "Red_Rock_Kept") else "Dungeness_Kept"
  results <- list()

  # P0 BUGFIX (2026-07): the caller passes `label` (e.g. "shore_all_gear") as
  # population_name, not the bare population. The old `is_shore_pe <-
  # (population_name == "shore")` was therefore ALWAYS FALSE, so every shore
  # component silently ran the boat branch: shore effort was expanded as
  # gear-deployments (gear x mean(number_of_gear) x tau) and its CPUE was
  # computed per gear rather than per crabber-hour. The PE was internally
  # consistent but on the wrong scale, and it was then compared against a BSS in
  # crabber-hours. Take the population explicitly, and fail loudly if absent.
  if (is.null(population)) {
    population <- if (grepl("^shore", population_name)) "shore"
                  else if (grepl("^private_boat", population_name)) "private_boat"
                  else stop("run_pe(): cannot infer population from '", population_name,
                            "'. Pass population = explicitly.", call. = FALSE)
  }
  stopifnot(population %in% c("shore", "private_boat"))

  # v5.2 Fix 1: Filter incomplete trips for CPUE estimation (PE method)
  interview_for_cpue <- summ$interview
  if(isTRUE(params$filter_incomplete_trips)) {
    interview_for_cpue <- interview_for_cpue |>
      filter(trip_status == "Complete" | is.na(trip_status))
  }

  # --- Effort expansion basis, by population (F2) ------------------------------
  # SHORE  : effort = crabbers x effective day length (hours). days$day_length is
  #          the I/E-derived L_effective set by bss_assign_day_length().
  # BOAT   : effort = gear DEPLOYMENTS = trailers x gear_per_group x tau.
  #          Boat catch is not linear in soak time (fitted exponent 0.133), so
  #          neither crabber-hours nor gear-hours is a stable effort unit; the
  #          deployment is. This puts the boat PE on the SAME scale as the boat
  #          BSS, so the BSS-PE gap is a model-disagreement diagnostic again
  #          rather than a unit artifact.
  is_shore_pe <- (population == "shore")

  if(is_shore_pe) {
    # P1: mirror the BSS's effort unit exactly, so the PE-BSS gap is a model
    # disagreement rather than a unit artifact, and so the assertion below holds.
    #   crabber-hours   : crabbers (= gear x crabbers_per_gear) x L_eff
    #   gear-hours      : gear                                   x L_eff
    #   gear-deployments: gear                                   x tau_shore
    eff_spec_pe <- bss_effort_spec(TRUE, days, params)
    effort_unit_pe <- eff_spec_pe$unit
    tau_shore_pe <- params$tau_shore_prior_mu %||% 1.7
    days <- days |> mutate(
      L_pe = if(effort_unit_pe == "gear-deployments") tau_shore_pe else day_length)
    gear_mult <- if(eff_spec_pe$effort_scale_gear == 1L) 1.0 else summ$crabbers_per_gear
    daily_effort <- summ$effort_index |>
      filter(count_sequence <= params$bss_max_count_seq) |>
      group_by(event_date, section_num) |>
      summarise(mean_count=mean(count_quantity), n_counts=n(), .groups="drop") |>
      mutate(est_units = mean_count * gear_mult) |>
      left_join(days |> select(event_date,day_type,L_pe,period), by="event_date") |>
      mutate(est_daily_effort = est_units * L_pe)
  } else {
    # Gear per boat group, from interviews where recorded (matches the BSS's
    # R_G_boat, which is learned from Gear_A_boat ~ Poisson(R_G_boat)).
    gpg_pe <- params$gear_per_group_default %||% 4.0
    ng_pe  <- suppressWarnings(as.numeric(summ$interview$number_of_gear))
    ng_pe  <- ng_pe[!is.na(ng_pe) & ng_pe > 0]
    if(length(ng_pe) > 0) gpg_pe <- mean(ng_pe)
    tau_pe <- params$tau_boat_prior_mu %||% 1.2
    effort_unit_pe <- "gear-deployments"
    cat(sprintf("  PE boat scale: gear_per_group=%.2f, tau=%.2f (deployments)\n", gpg_pe, tau_pe))
    days <- days |> mutate(L_pe = tau_pe)
    daily_effort <- summ$effort_index |>
      filter(count_sequence <= params$bss_max_count_seq) |>
      group_by(event_date, section_num) |>
      summarise(mean_count=mean(count_quantity), n_counts=n(), .groups="drop") |>
      left_join(days |> select(event_date,day_type,L_pe,period), by="event_date") |>
      mutate(est_daily_effort = mean_count * gpg_pe * L_pe)
  }
  results$effort_unit <- effort_unit_pe

  total_days_strat <- days |> filter(open_section_1) |>
    group_by(period, day_type) |> summarise(n_total_days=n(), .groups="drop")

  effort_strat <- daily_effort |>
    group_by(section_num, period, day_type) |>
    summarise(mean_daily=mean(est_daily_effort,na.rm=TRUE),
              sd_daily=sd(est_daily_effort,na.rm=TRUE),
              n_sampled=n(), .groups="drop") |>
    left_join(total_days_strat, by=c("period","day_type")) |>
    mutate(est_total = mean_daily * n_total_days,
           se_total = sqrt((n_total_days^2)*replace_na(sd_daily^2,0)/pmax(n_sampled,1)))

  results$effort_total <- sum(effort_strat$est_total, na.rm=TRUE)
  results$effort_se <- sqrt(sum(effort_strat$se_total^2, na.rm=TRUE))

  for(cg in catch_groups) {
    if(!cg %in% names(summ$interview)) { results[[cg]] <- 0; next }

    # F2: CPUE denominator matches the effort unit above. Ratio-of-sums within a
    # day (sum catch / sum denominator), which is the estimator consistent with a
    # Poisson/NB rate; do NOT switch this to a mean of per-interview ratios.
    # P1: CPUE denominator must be the same unit as the effort above.
    h_col_pe <- if(is_shore_pe) bss_effort_spec(TRUE, days, params)$h_col else "number_of_gear"
    daily_cpue <- interview_for_cpue |>
      mutate(.h_pe = suppressWarnings(as.numeric(.data[[h_col_pe]]))) |>
      # Numerator and denominator must come from the SAME interviews. Without this
      # filter, sum(.h_pe, na.rm = TRUE) silently drops rows from the denominator
      # while sum(catch) keeps them in the numerator. Currently inert (every boat
      # interview lacking number_of_gear has zero kept crab), but it is exactly the
      # failure mode that would bias CPUE upward the moment that stops being true.
      # This mirrors the CPUE denominator filter in prep_bss_crab().
      filter(is.finite(.h_pe), .h_pe > 0) |>
      group_by(event_date, section_num) |>
      summarise(catch = sum(.data[[cg]], na.rm=TRUE),
                hrs = sum(.h_pe, na.rm=TRUE),
                n_int=n(), .groups="drop") |>
      mutate(cpue=if_else(hrs>0, catch/hrs, 0)) |>
      left_join(days |> select(event_date,day_type,period), by="event_date")

    # P0: stratum CPUE is a RATIO-OF-SUMS within the stratum, sum(catch)/sum(h),
    # not a weighted mean of per-day ratios. The latter is a mean-of-ratios
    # estimator: with weekly strata and ~50% day coverage many strata rest on one
    # or two days, and a day with very little sampled effort produces an extreme
    # daily ratio that is then multiplied by the full stratum effort. The
    # ratio-of-sums is the estimator consistent with a Poisson/NB rate and is what
    # the BSS's lambda_C converges to as r_C -> Inf.
    cpue_strat <- daily_cpue |>
      group_by(section_num, period, day_type) |>
      summarise(catch_s = sum(catch, na.rm=TRUE),
                hrs_s   = sum(hrs,   na.rm=TRUE),
                n_int_s = sum(n_int, na.rm=TRUE), .groups="drop") |>
      mutate(mean_cpue = if_else(hrs_s > 0, catch_s / hrs_s, NA_real_))

    # Empty-stratum CPUE fallback (item 2, 2026-07-13; mirrors run_pe_pooled). An
    # effort-bearing stratum with no surviving interviews gets mean_cpue = NA; the old
    # replace_na(., 0) zeroed its catch (under-count, and the source of the boat PE's
    # sparse-stratum sign instability). params$pe_empty_stratum = "pooled" (default)
    # fills it with the population x sub-season ratio-of-sums CPUE; "zero" is the old behavior.
    pooled_cpue <- if (sum(daily_cpue$hrs, na.rm=TRUE) > 0)
                     sum(daily_cpue$catch, na.rm=TRUE) / sum(daily_cpue$hrs, na.rm=TRUE) else 0
    empty_fill  <- if (identical(params$pe_empty_stratum %||% "pooled", "zero")) 0 else pooled_cpue
    catch_strat <- effort_strat |>
      left_join(cpue_strat, by=c("section_num","period","day_type")) |>
      mutate(est_catch = est_total * replace_na(mean_cpue, empty_fill))

    results[[cg]] <- sum(catch_strat$est_catch, na.rm=TRUE)

    # P0: the PE's implied CPUE (catch / effort) must agree with the ratio-of-sums
    # over the interviews it was built from. If it does not, catch and effort are
    # on different scales, which is exactly the failure this fix removes.
    implied <- if (results$effort_total > 0) results[[cg]] / results$effort_total else NA_real_
    ros     <- if (sum(daily_cpue$hrs, na.rm=TRUE) > 0)
                 sum(daily_cpue$catch, na.rm=TRUE) / sum(daily_cpue$hrs, na.rm=TRUE) else NA_real_
    if (is.finite(implied) && is.finite(ros) && ros > 0) {
      rel <- implied / ros
      cat(sprintf("  PE check [%s / %s]: implied CPUE %.4f vs interview ratio-of-sums %.4f (%.2fx) [%s]\n",
                  population_name, cg, implied, ros, rel, effort_unit_pe))
      if (rel < 0.5 || rel > 2.0) {
        stop(sprintf(paste0("run_pe(): PE implied CPUE (%.4f) is %.2fx the interview ",
                            "ratio-of-sums (%.4f) for %s / %s. Catch and effort are not on ",
                            "the same scale."),
                     implied, rel, ros, population_name, cg), call. = FALSE)
      }
    }
  }

  rr_str <- if(params$estimate_red_rock) sprintf(", RR=%s", format(round(results$Red_Rock_Kept),big.mark=",")) else ""
  cat(sprintf("  PE %s: Effort=%s %s, Dung=%s%s\n", population_name,
              format(round(results$effort_total),big.mark=","), effort_unit_pe,
              format(round(results$Dungeness_Kept),big.mark=","), rr_str))

  return(results)
}
