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
# run_pe_pooled.R  (pooled-CPUE driver)
#
# Point Estimator (PE) for the pooled-CPUE model: stratified effort and catch
# expansion per population x sub-season. Extracted from the pooled driver and
# renamed run_pe_pooled (distinct from the gear-resolved run_pe_gear, since the
# whole 03_R_functions folder is sourced by both drivers). catch_groups is derived
# from params.
#
# 2026-07-11 shore-scale fix (finding F1): the shore effort expansion and CPUE
# denominator now flow through the shared 03_R_functions/bss_effort_spec.R, so the
# shore PE runs on the SAME effort unit as the shore BSS (gear-deployments in
# production, per params$shore_effort_unit). Before this fix the shore branch was
# still hard-coded to crabber-hours (est_crabbers * day_length; CPUE denominator
# fishing_time_total) even though pooled v7.7 moved the shore BSS onto
# gear-deployments, so the shore PE and shore BSS were on different units and the
# pe_vs_bss_comparison / monthly PE effort-share for shore were unit-inconsistent.
# NOTE: this MOVES the pooled shore PE number; confirm against a validation run.
# The stratum-CPUE estimator (weighted mean of daily ratios) is intentionally left
# unchanged; adopting the gear track's ratio-of-sums (P0) is a separate decision.
###############################################################################

run_pe_pooled <- function(summ, days, params, population_name) {
  catch_groups <- if (isTRUE(params$estimate_red_rock)) c("Dungeness_Kept", "Red_Rock_Kept") else "Dungeness_Kept"
  results <- list()
  is_boat <- str_detect(population_name, "private_boat")

  daily_effort <- summ$effort_index |>
    filter(count_sequence <= params$bss_max_count_seq) |>
    group_by(event_date, section_num) |>
    summarise(mean_count=mean(count_quantity), n_counts=n(), .groups="drop") |>
    mutate(est_crabbers = mean_count * summ$crabbers_per_gear) |>
    left_join(days |> select(event_date,day_type,day_length,period), by="event_date")

  if(is_boat) {
    # POOL-3: boat effort on the gear-DEPLOYMENT scale (matches the BSS via
    # bss_effort_spec()): gear_per_group * tau_boat deployments per day, replacing
    # the old flat gear-hours (gear_per_group * 24). Keeping PE and BSS on the same
    # unit makes the PE-vs-BSS gap a model disagreement, not a unit artifact.
    ratio_data <- summ$interview |>
      filter(!is.na(number_of_gear), number_of_gear > 0, angler_count > 0)
    gear_per_group <- if(nrow(ratio_data) > 0) mean(ratio_data$number_of_gear)
                      else (params$gear_per_group_default %||% 4.0)
    tau_boat_pe <- params$tau_boat_prior_mu %||% 1.2
    daily_effort <- daily_effort |>
      mutate(est_daily_effort = mean_count * gear_per_group * tau_boat_pe)
    effort_unit_pe <- "gear-deployments"
    cat(sprintf("  PE %s: gear_per_group=%.2f, tau=%.2f (gear-deployments)\n",
                population_name, gear_per_group, tau_boat_pe))
  } else {
    # POOL-7 (v7.7) shore-scale fix: shore effort on the unit set by
    # params$shore_effort_unit, via the shared bss_effort_spec(), so the PE matches
    # the shore BSS. gear-deployments (production): effort = gear_count * tau_shore.
    # crabber-hours (pre-v7.7): effort = crabbers * L_eff, which reproduces the old
    # est_crabbers * day_length exactly (gear_mult = crabbers_per_gear, L_pe = day_length).
    eff_spec_pe    <- bss_effort_spec(TRUE, days, params)
    effort_unit_pe <- eff_spec_pe$unit
    tau_shore_pe   <- params$tau_shore_prior_mu %||% 1.7
    gear_mult      <- if (eff_spec_pe$effort_scale_gear == 1L) 1.0 else summ$crabbers_per_gear
    daily_effort <- daily_effort |>
      mutate(L_pe = if (effort_unit_pe == "gear-deployments") tau_shore_pe else day_length,
             est_daily_effort = mean_count * gear_mult * L_pe)
  }

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

  # CPUE denominator matches the effort unit above. Boat: number_of_gear
  # (deployments). Shore: bss_effort_spec()$h_col (number_of_gear for deployments;
  # fishing_time_total for crabber-hours), so the shore PE CPUE denominator moves
  # in lockstep with the shore effort unit and the shore BSS.
  hrs_col <- if(is_boat) "number_of_gear" else bss_effort_spec(TRUE, days, params)$h_col

  # P0: accumulate the per-catch-group PE internal-consistency check (implied CPUE vs
  # interview ratio-of-sums) so the driver can render it in the report (Section 4).
  pe_check_rows <- list()
  for(cg in catch_groups) {
    if(!cg %in% names(summ$interview)) { results[[cg]] <- 0; next }
    # POOL-2: drop incomplete trips from the PE CPUE too, so the PE and BSS CPUE
    # definitions stay consistent. Toggle via params$filter_incomplete_trips;
    # missing trip_status is kept.
    int_cpue <- summ$interview |>
      filter(!is.na(.data[[hrs_col]]), .data[[hrs_col]] > 0)
    if(isTRUE(params$filter_incomplete_trips)) {
      int_cpue <- int_cpue |> filter(trip_status == "Complete" | is.na(trip_status))
    }
    daily_cpue <- int_cpue |>
      group_by(event_date, section_num) |>
      summarise(catch=sum(.data[[cg]],na.rm=TRUE), hrs=sum(.data[[hrs_col]],na.rm=TRUE),
                n_int=n(), .groups="drop") |>
      mutate(cpue=if_else(hrs>0, catch/hrs, 0)) |>
      left_join(days |> select(event_date,day_type,period), by="event_date")
    # P0 (2026-07-12): stratum CPUE is a RATIO-OF-SUMS within the stratum,
    # sum(catch)/sum(hrs), not a weighted mean of per-day ratios. The mean-of-ratios
    # (weighted.mean(cpue, w = n_int)) is unstable: with weekly strata and ~50% day
    # coverage many strata rest on one or two sampled days, and a day with very
    # little sampled effort produces an extreme daily ratio that is then multiplied
    # by the full stratum effort (pre-fix this made the shore PE implied CPUE 2.9x
    # its own ratio-of-sums). Ratio-of-sums is the estimator consistent with the
    # Poisson/NB rate and what the BSS lambda_C converges to as r_C -> Inf; it also
    # matches run_pe_gear, so the two tracks' PE cannot drift. (Backlog P0.)
    cpue_strat <- daily_cpue |>
      group_by(section_num, period, day_type) |>
      summarise(catch_s = sum(catch, na.rm=TRUE),
                hrs_s   = sum(hrs,   na.rm=TRUE), .groups="drop") |>
      mutate(mean_cpue = if_else(hrs_s > 0, catch_s / hrs_s, NA_real_))
    # Empty-stratum CPUE fallback (item 2, 2026-07-13). A stratum with expanded effort
    # but no surviving interviews (after the incomplete-trip filter) gets mean_cpue = NA.
    # The old behavior (replace_na(mean_cpue, 0)) assigned it ZERO catch, which
    # under-counts, and under thin BOAT sampling with weekly strata it made the boat PE
    # swing on a single empty cell (the incomplete-trip "anomaly": the boat filter effect
    # was a knife-edge because one week x day-type cell's only boat interviews were
    # incomplete and got zeroed). params$pe_empty_stratum: "pooled" (default) fills an
    # empty stratum with the population x sub-season ratio-of-sums CPUE (a sampled rate is
    # a better guess than zero and matches the P0 consistency target); "zero" restores the
    # old behavior. Shore is dense (few or no empty strata), so this mainly steadies the boat.
    pooled_cpue <- if (sum(daily_cpue$hrs, na.rm=TRUE) > 0)
                     sum(daily_cpue$catch, na.rm=TRUE) / sum(daily_cpue$hrs, na.rm=TRUE) else 0
    empty_fill  <- if (identical(params$pe_empty_stratum %||% "pooled", "zero")) 0 else pooled_cpue
    catch_strat <- effort_strat |>
      left_join(cpue_strat, by=c("section_num","period","day_type")) |>
      mutate(est_catch = est_total * replace_na(mean_cpue, empty_fill))
    results[[cg]] <- sum(catch_strat$est_catch, na.rm=TRUE)

    # P0: the PE's implied CPUE (catch / effort) must agree with the ratio-of-sums
    # over the interviews it was built from. A divergence beyond 2x means catch and
    # effort are on different scales, the failure this fix removes. This runs in the
    # PE section, before the multi-hour BSS fits, so it fails fast. Mirrors run_pe_gear.
    implied <- if (results$effort_total > 0) results[[cg]] / results$effort_total else NA_real_
    ros     <- if (sum(daily_cpue$hrs, na.rm=TRUE) > 0)
                 sum(daily_cpue$catch, na.rm=TRUE) / sum(daily_cpue$hrs, na.rm=TRUE) else NA_real_
    if (is.finite(implied) && is.finite(ros) && ros > 0) {
      rel <- implied / ros
      cat(sprintf("  PE check [%s / %s]: implied CPUE %.4f vs interview ratio-of-sums %.4f (%.2fx) [%s]\n",
                  population_name, cg, implied, ros, rel, effort_unit_pe))
      pe_check_rows[[cg]] <- tibble(catch_group = cg, effort_unit = effort_unit_pe,
                                    pe_implied_cpue = implied, interview_ros = ros, ratio = rel)
      if (rel < 0.5 || rel > 2.0) {
        stop(sprintf(paste0("run_pe_pooled(): PE implied CPUE (%.4f) is %.2fx the interview ",
                            "ratio-of-sums (%.4f) for %s / %s. Catch and effort are not on ",
                            "the same scale."),
                     implied, rel, ros, population_name, cg), call. = FALSE)
      }
    }
  }

  results$effort_unit <- effort_unit_pe
  results$pe_cpue_check <- if (length(pe_check_rows)) bind_rows(pe_check_rows) else NULL
  rr_str <- if(params$estimate_red_rock) sprintf(", RR=%s", format(round(results$Red_Rock_Kept),big.mark=",")) else ""
  cat(sprintf("  PE %s: Effort=%s %s, Dung=%s%s\n", population_name,
              format(round(results$effort_total),big.mark=","), effort_unit_pe,
              format(round(results$Dungeness_Kept),big.mark=","), rr_str))
  return(results)
}
