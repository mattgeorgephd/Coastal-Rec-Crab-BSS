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
# diagnose_incomplete_trips.R   (improvement 5, 2026-08-25; DIAGNOSTIC ONLY)
#
# An incomplete trip is one where the crabber had not finished when we spoke to them.
# Their gear is all in the CPUE denominator; only part of their catch is in the
# numerator, so the measured catch rate reads low: -21.2% for pots, -23.2% for traps,
# -19.9% for snares, +4.4% for ring nets (which are checked every few minutes and are
# therefore effectively exempt). Production drops them, which costs 36% of the sample
# (1,334 of 3,759 interviews on 2024-25) and raises the PE by +9.4% (shore pot closure),
# +11.7% (shore all gear) and +5.0% (boat all gear).
#
# THE STRUCTURAL POINT THAT DECIDES WHAT IS EVEN POSSIBLE HERE
#   Effort in this pipeline does NOT come from interviews. It comes from the gear counts
#   (Float 20 + 17-21) and the trailer/OSP boat counts. Interviews supply only two
#   things: the catch RATE, and the gear RATIOS (R_G, gear per crabber; R_G_boat / the PE
#   gear-per-group, gear per boat group). So no treatment of incomplete trips can add
#   effort to the expansion. Every option below moves the estimate through the catch rate
#   or the gear ratios, and nothing else.
#
# THE FOUR ARMS
#   exclude            Production today. Drop them from the CPUE estimate AND from the
#                      gear ratios.
#   gear_only          Drop their CATCH from the CPUE estimate, KEEP their gear counts in
#                      the gear ratios. An interrupted trip's gear count is fully
#                      observed -- the crabber has N pots out whether or not they are
#                      done -- so only the catch is truncated.
#   impute_mean_cpue   Keep them, replacing observed catch with (complete-trip
#                      ratio-of-sums CPUE x their gear count).
#   keep               No filter at all (pre-v7.5 behaviour).
#
# WHAT TO EXPECT, SO THE RESULT CAN BE READ AS CONFIRMATION RATHER THAN NEWS
#   * gear_only is the arm that can actually move something. The pipeline does NOT do
#     this today: prep_bss_crab_*() builds intA (the R_G / R_G_boat interview set) from
#     the ALREADY-FILTERED frame, so an incomplete trip's gear count is discarded. Note
#     the PE and the BSS currently DISAGREE on this point: run_pe_*() computes the boat
#     gear-per-group from the UNFILTERED interview set, so the boat PE already behaves
#     like gear_only while the boat BSS behaves like exclude. The `gear_ratio` columns
#     below make that visible.
#   * impute_mean_cpue should be close to a no-op ON THE CPUE, and that is arithmetic,
#     not luck: a ratio-of-sums over (complete + imputed) returns the complete-trip CPUE
#     identically, because the imputed numerator is by construction that CPUE times the
#     imputed denominator. Any movement is the re-weighting of thin week x day-type
#     strata toward the pooled rate, i.e. a shrinkage effect, not new information. It is
#     included so "no-op" is a measured result rather than a claim -- and as a warning
#     that in the BSS the same idea would be worse than a no-op, because it would add
#     fabricated observations that tighten the CPUE credible interval without adding
#     information.
#   * LENGTH-BIASED SAMPLING is the one thing that could undermine gear_only. Intercepted
#     trips are, by construction, over-representative of long trips; if long trips also
#     deploy more gear, the incomplete subset is gear-heavier and folding it into R_G
#     would bias the effort expansion UP. The `gear_ratio` column (per arm) and the
#     `gear_lengthbias_p` column (Welch t-test of incomplete vs complete gear counts,
#     one value per component) are that check. Read them before acting on the gear_only
#     arm: a small p with a HIGHER incomplete-trip gear count is the signature that would
#     disqualify it.
#
# WHAT IS NOT HERE
#   The statistically correct way to recover all 1,334 interviews is a CENSORED
#   likelihood: an incomplete trip's observed catch is a LOWER bound on its trip total,
#   so its contribution becomes neg_binomial_2_lccdf(c - 1 | lambda_C * h, r_C) instead of
#   the pmf. That is a Stan change with a real validation burden and it cannot be
#   evaluated by a design-based PE arm, so it is documented as a backlog item rather than
#   implemented here. See PIPELINE_STATUS Tier 2.
#
#   Also note what NONE of these arms fix: excluding incomplete trips does not remove
#   selection bias, it sharpens it. Field protocol already favours crabbers who have
#   finished, so early leavers are over-represented in what remains. There is no
#   diagnostic for that here and none in the pipeline.
#
# RETURNS a tibble, one row per (component x arm), written to
# sensitivity_incomplete_trips.csv by the caller.
###############################################################################

# Apply one arm to a population x sub-season interview frame. Returns the frame the CPUE
# estimate should see, plus the frame the gear ratios should see. `h_col` is the CPUE
# denominator for this population (from bss_effort_spec).
incomplete_trip_arm_frames <- function(interview, arm, h_col, catch_col = "Dungeness_Kept") {
  inc <- !is.na(interview$trip_status) & interview$trip_status == "Incomplete"
  valid_h <- is.finite(suppressWarnings(as.numeric(interview[[h_col]]))) &
             suppressWarnings(as.numeric(interview[[h_col]])) > 0

  complete_frame <- interview[!inc & valid_h, , drop = FALSE]

  cpue_frame <- switch(arm,
    "exclude"          = complete_frame,
    "gear_only"        = complete_frame,
    "keep"             = interview[valid_h, , drop = FALSE],
    "impute_mean_cpue" = {
      num <- sum(suppressWarnings(as.numeric(complete_frame[[catch_col]])), na.rm = TRUE)
      den <- sum(suppressWarnings(as.numeric(complete_frame[[h_col]])),     na.rm = TRUE)
      cpue_hat <- if (den > 0) num / den else 0
      imp <- interview[inc & valid_h, , drop = FALSE]
      if (nrow(imp) > 0)
        imp[[catch_col]] <- as.numeric(imp[[h_col]]) * cpue_hat
      dplyr::bind_rows(complete_frame, imp)
    },
    stop("Unknown incomplete-trip arm: ", arm, call. = FALSE))

  # The gear-ratio frame: which interviews are allowed to inform R_G / gear-per-group.
  # "exclude" is the only arm that throws away an incomplete trip's fully-observed gear
  # count.
  gear_frame <- if (identical(arm, "exclude")) complete_frame else interview[valid_h, , drop = FALSE]

  list(cpue = cpue_frame, gear = gear_frame, n_incomplete = sum(inc), n_total = nrow(interview))
}


# pe_fun: the track's own point estimator (run_pe_pooled or run_pe_gear), so each driver
# compares the arms against ITS production estimator rather than a re-implementation.
# Which arm each ESTIMATOR is actually running (2026-08-27).
#
# `filter_incomplete_trips = TRUE` drops incomplete trips from the CPUE likelihood and from
# the interview frame prep_bss_crab_*() derives R_G from, so the BSS runs `exclude`. The PE
# does not follow: run_pe_*() takes the boat gear-per-group from the UNFILTERED interview
# set, so the boat PE is already running `gear_only` while the boat BSS runs `exclude`. The
# shore PE counts gear directly and is unaffected either way. That asymmetry is improvement
# 5's central finding, and the 2026-08-26 ladder confirmed it numerically (the shipped boat
# PE equals the gear_only arm of this table), so it must be labelled rather than flattened.
incomplete_trip_production_arm <- function(is_shore, filter_incomplete_trips = TRUE) {
  if (!isTRUE(filter_incomplete_trips)) return(c(bss = "keep", pe = "keep"))
  c(bss = "exclude", pe = if (isTRUE(is_shore)) "exclude" else "gear_only")
}

diagnose_incomplete_trips <- function(dwg, subseasons, params, L_eff_model,
                                      catch_col = "Dungeness_Kept",
                                      pe_fun = run_pe_pooled) {
  arms <- params$incomplete_trip_arms %||% c("exclude", "gear_only", "impute_mean_cpue", "keep")
  rows <- list()

  for (pop in c("shore", "private_boat")) {
    for (ss in subseasons) {
      key <- paste0(pop, "_", ss$name)
      days_ss <- prep_days_crab(ss$start, ss$end, params, L_eff_model = L_eff_model)
      summ_ss <- prep_population_summary(dwg, pop, ss$start, ss$end, params)
      int     <- summ_ss$interview
      if (is.null(int) || nrow(int) == 0) next

      is_shore <- (pop == "shore")
      h_col <- if (is_shore) bss_effort_spec(TRUE, days_ss, params)$h_col else "number_of_gear"
      if (!h_col %in% names(int) || !catch_col %in% names(int)) next

      # Length-bias check on the gear counts themselves, once per component.
      ng   <- suppressWarnings(as.numeric(int$number_of_gear))
      incf <- !is.na(int$trip_status) & int$trip_status == "Incomplete"
      g_c  <- ng[!incf & is.finite(ng) & ng > 0]
      g_i  <- ng[ incf & is.finite(ng) & ng > 0]
      lb_p <- if (length(g_c) > 1 && length(g_i) > 1)
                tryCatch(stats::t.test(g_i, g_c)$p.value, error = function(e) NA_real_) else NA_real_

      for (arm in arms) {
        fr <- incomplete_trip_arm_frames(int, arm, h_col, catch_col)

        # Design-based PE under this arm. run_pe_*() is reused verbatim so the arm is
        # compared against the production estimator, not a re-implementation of it; the
        # arm is expressed by handing it a modified interview frame and switching its own
        # filter off (the arm has already applied whatever filtering it wants).
        # summ$interview drives the CPUE; summ$interview_gear (new, optional) drives the
        # gear ratio, so the gear_only arm can keep an incomplete trip's gear count while
        # still dropping its truncated catch. run_pe_*() falls back to summ$interview
        # when interview_gear is absent, so every other caller is unaffected.
        summ_arm <- summ_ss
        summ_arm$interview      <- fr$cpue
        summ_arm$interview_gear <- fr$gear
        p_arm <- modifyList(params, list(filter_incomplete_trips = FALSE))
        pe <- tryCatch(pe_fun(summ_arm, days_ss, p_arm, paste0(key, " [", arm, "]")),
                       error = function(e) NULL)

        # Numerator and denominator must come from the SAME rows, and from the same
        # subset prep_bss_crab_*() uses to build intA (number_of_gear > 0 AND
        # angler_count > 0), or this column is not comparable to the R_G the model
        # actually learns -- which is exactly the comparison this section asks for.
        .g <- suppressWarnings(as.numeric(fr$gear$number_of_gear))
        .a <- suppressWarnings(as.numeric(fr$gear$angler_count))
        keep_ga  <- is.finite(.g) & .g > 0 & is.finite(.a) & .a > 0
        gear_all <- .g[keep_ga]; ang <- .a[keep_ga]
        gear_ratio <- if (!length(gear_all)) NA_real_
                      else if (is_shore) sum(gear_all) / sum(ang)
                      else mean(gear_all)

        rows[[length(rows) + 1]] <- tibble::tibble(
          component        = key,
          component_disp   = paste0(gsub("_", " ", pop), " ", ss$display_name),
          arm              = arm,
          cpue_denominator = h_col,
          n_interviews     = fr$n_total,
          n_incomplete     = fr$n_incomplete,
          pct_incomplete   = round(100 * fr$n_incomplete / max(fr$n_total, 1), 1),
          n_cpue_used      = nrow(fr$cpue),
          n_gear_used      = length(gear_all),
          gear_ratio       = round(gear_ratio, 3),
          # NOTE on where the gear ratio bites. For the BOAT it multiplies the PE effort
          # directly (trailers x gear-per-group x tau), so the arm moves PE_effort and
          # PE_catch. For SHORE on gear-deployments it does NOT enter the PE at all (the
          # PE counts gear directly); it moves only the BSS, through R_G in
          # Gear_I ~ NB2(lambda_E * R_G, r_E) and E = lambda_E * R_G * tau. So a flat
          # shore PE across arms is expected and is not evidence that the arm is inert.
          gear_ratio_note  = if (is_shore) "gear per crabber (R_G; BSS only)"
                             else "gear per boat group (R_G_boat; PE + BSS)",
          is_shore_component = is_shore,
          PE_effort        = if (is.null(pe)) NA_real_ else round(pe$effort_total),
          PE_catch         = if (is.null(pe)) NA_real_ else round(pe[[catch_col]] %||% NA_real_),
          gear_lengthbias_p = round(lb_p, 4)
        )
      }
    }
  }

  if (!length(rows)) return(NULL)
  out <- dplyr::bind_rows(rows)

  # Express every arm against the production arm so the table reads as "what would change".
  base <- out |> dplyr::filter(.data$arm == "exclude") |>
    dplyr::select(component, base_catch = PE_catch, base_effort = PE_effort,
                  base_gear = gear_ratio)
  out |>
    dplyr::left_join(base, by = "component") |>
    dplyr::mutate(
      catch_vs_exclude_pct = ifelse(is.finite(base_catch) & base_catch > 0,
                                    round(100 * (PE_catch - base_catch) / base_catch, 1), NA_real_),
      gear_vs_exclude_pct  = ifelse(is.finite(base_gear) & base_gear > 0,
                                    round(100 * (gear_ratio - base_gear) / base_gear, 1), NA_real_),
      # 2026-08-27: production_arm is now PER ESTIMATOR, not one label for the whole table.
      # A single "exclude" was the honest label for the BSS and the wrong one for the boat PE,
      # and the 2026-08-26 ladder made that concrete: the shipped boat PE (3,565.75 effort /
      # 10,940.36 catch) equals this table's `gear_only` arm, not its `exclude` arm. That is not
      # a bug in the diagnostic, it is improvement 5's central finding -- run_pe_*() takes the
      # boat gear-per-group from the UNFILTERED interview set, so the boat PE already behaves
      # like gear_only while the boat BSS behaves like exclude. A table that labels both
      # "exclude" hides exactly the asymmetry it exists to surface.
      production_arm_bss = incomplete_trip_production_arm(
                             TRUE, params$filter_incomplete_trips)[["bss"]],
      production_arm_pe  = vapply(.data$is_shore_component, function(z)
                             incomplete_trip_production_arm(z, params$filter_incomplete_trips)[["pe"]],
                             character(1)),
      production_arm     = production_arm_bss) |>
    dplyr::select(-base_catch, -base_effort, -base_gear, -is_shore_component)
}
