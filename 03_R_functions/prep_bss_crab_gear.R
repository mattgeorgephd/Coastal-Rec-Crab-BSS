###############################################################################
# prep_bss_crab_gear.R  (gear-resolved driver)
#
# Build the Stan data list for crab_bss_gear_resolved.stan (per-gear CPUE machinery,
# weighted gear proportions, gear_exclude, effort-unit spec via bss_effort_spec).
# Extracted from the gear-resolved driver and renamed prep_bss_crab_gear. Pure given
# its arguments.
###############################################################################

prep_bss_crab_gear <- function(days, summ, est_catch_group, params, population_name, period_type,
                          gear_exclude = character(0), ie_data = NULL) {

  eff <- summ$effort_index |> filter(count_sequence <= params$bss_max_count_seq)
  D <- nrow(days); G <- 1L; S <- 1L

  # --- BSS period / AR temporal resolution ---------------------------------
  # Shared selector: 03_R_functions/bss_ar_resolution.R. eff_d is built first
  # because the selector reports effort coverage. Default (params$ar_adaptive =
  # FALSE) passes the sub-season's explicit period_bss, reproducing the v5.x
  # fixed-period behavior exactly. Setting ar_adaptive = TRUE lets the data
  # choose, subject to params$ar_max_resolution. params$ar_force overrides both.
  eff_d <- eff |> left_join(days |> select(event_date,day_index), by="event_date") |> filter(!is.na(day_index))

  ar_sel <- bss_select_ar_resolution(
    days             = days,
    eff_d            = eff_d,
    population_name  = population_name,
    params           = params,
    fixed_resolution = if (isTRUE(params$ar_adaptive)) NULL else period_type
  )
  ar_resolution <- ar_sel$resolution
  P_n  <- ar_sel$P_n
  pvec <- ar_sel$pvec

  # v5.3: is_shore is computed early so that the boat-specific
  # gear_time_total filter and the L/h branches below can use it.
  is_shore <- (population_name == "shore")

  int_cg <- summ$interview |>
    mutate(fish_count = .data[[est_catch_group]]) |>
    filter(!is.na(fishing_time_total), fishing_time_total > 0)

  int_d <- int_cg |>
    left_join(days |> select(event_date,day_index,day_type,day_type_idx), by="event_date") |>
    filter(!is.na(day_index)) |> distinct(interview_id, .keep_all=TRUE)

  # =========================================================================
  # v5.2 FIX 1: Incomplete trip filter
  # Crabbers still fishing at interview time have systematically lower CPUE
  # (-20% for pots/traps) because soak-time-dependent gear hasn't finished.
  # Ring nets are less affected (+4% bias) since they are checked frequently.
  # Including incomplete trips biases harvest estimates downward by ~7%.
  # =========================================================================
  n_before_trip_filter <- nrow(int_d)
  trip_status_counts <- int_d |> count(trip_status, name = "n")
  n_incomplete <- sum(int_d$trip_status == "Incomplete", na.rm = TRUE)
  n_complete <- sum(int_d$trip_status == "Complete", na.rm = TRUE)
  n_missing_status <- sum(is.na(int_d$trip_status))

  cat(sprintf("  Trip completion: %d complete, %d incomplete (%.0f%%), %d missing status\n",
              n_complete, n_incomplete,
              n_incomplete / max(n_complete + n_incomplete, 1) * 100,
              n_missing_status))

  # Report CPUE difference if both groups have data
  cpue_by_status <- int_d |>
    filter(fishing_time_total >= 0.5, fish_count >= 0) |>
    mutate(cpue = fish_count / fishing_time_total) |>
    group_by(trip_status) |>
    summarise(mean_cpue = mean(cpue), n = n(), .groups = "drop")
  if(nrow(cpue_by_status) >= 2) {
    cpue_comp <- cpue_by_status |> filter(trip_status == "Complete") |> pull(mean_cpue)
    cpue_inc <- cpue_by_status |> filter(trip_status == "Incomplete") |> pull(mean_cpue)
    if(length(cpue_comp) > 0 && length(cpue_inc) > 0 && cpue_comp > 0) {
      bias_pct <- (cpue_inc - cpue_comp) / cpue_comp * 100
      cat(sprintf("  CPUE: complete=%.3f, incomplete=%.3f (bias=%+.1f%%)\n",
                  cpue_comp, cpue_inc, bias_pct))
    }
  }

  if(isTRUE(params$filter_incomplete_trips)) {
    # Keep complete trips + missing status (which may be complete)
    int_d <- int_d |> filter(trip_status == "Complete" | is.na(trip_status))
    cat(sprintf("  Incomplete trip filter ON: %d → %d interviews (removed %d)\n",
                n_before_trip_filter, nrow(int_d), n_before_trip_filter - nrow(int_d)))
  } else {
    cat("  Incomplete trip filter OFF (v5.1 behavior)\n")
  }

  # --- 5b: I/E observations for this population and sub-season ---------------
  # Matched to day_index within the sub-season window. Shore only: the boat has no
  # shore I/E stream, so IE_n = 0 and sigma_IE is identified by its prior alone
  # (which is why that prior must be unconditional; see the Stan model's B1.6 note).
  # --- P1: effort-unit specification (shared with run_pe via bss_effort_spec) --
  eff_spec <- bss_effort_spec(is_shore, days, params)
  cat(sprintf("  Effort unit: %s (effort_scale_gear=%d)\n",
              eff_spec$unit, eff_spec$effort_scale_gear))

  # --- F2: I/E observation stream, population-aware ---------------------------
  #   shore: observation = crabber-hours; predicted = lambda_E * L        (hours)
  #   boat : observation = boat trips;    predicted = (lambda_E/R_G_boat) * tau
  # The boat stream is what identifies tau. It is inert this season (no WBL days
  # inside the window) but activates automatically as WBL I/E accumulates.
  ie_match <- tibble(day_index = integer(), ie_obs = numeric())
  ie_group_scale <- as.integer(!is_shore)

  if(!is.null(ie_data) && nrow(ie_data) > 0) {
    ie_pop  <- if(is_shore) "shore" else "private_boat"
    ie_want <- is_shore || isTRUE(params$use_boat_ie)
    if(ie_want) {
      ie_match <- ie_data |>
        filter(population == ie_pop) |>
        inner_join(days |> select(event_date, day_index), by = "event_date") |>
        filter(!is.na(day_index))
      ie_match <- if(is_shore) {
        ie_match |> mutate(ie_obs = ie_crabber_hours) |> filter(ie_obs > 0)
      } else {
        ie_match |> mutate(ie_obs = as.numeric(ie_trips)) |> filter(ie_obs > 0)
      }
      # Require a minimum number of boat I/E days before letting them inform tau.
      if(!is_shore && nrow(ie_match) < (params$ie_min_obs_boat %||% 2)) {
        if(nrow(ie_match) > 0)
          cat(sprintf("  Boat I/E: %d day(s) < ie_min_obs_boat; not used (tau stays prior-driven)\n",
                      nrow(ie_match)))
        ie_match <- ie_match[0, ]
      }
      # GR-8 (2026-07-13): symmetric shore guard. With only 1-2 in-window shore I/E days
      # the lognormal I/E likelihood barely informs L but lets sigma_IE (exponential(5),
      # mode at 0) shrink and stiffen against those few points, the funnel seen in the
      # sparse shore ring-net fit. Drop the stream so sigma_IE stays prior-only (decoupled).
      # Set ie_min_obs_shore = 0 to disable. Prior kept as exponential(5) on purpose
      # (tightening it would push the shore all-gear sigma_IE down; see item 4 / GR-9).
      if(is_shore && nrow(ie_match) > 0 && nrow(ie_match) < (params$ie_min_obs_shore %||% 3)) {
        cat(sprintf("  Shore I/E: %d day(s) < ie_min_obs_shore; not used (sigma_IE decoupled, GR-8)\n",
                    nrow(ie_match)))
        ie_match <- ie_match[0, ]
      }
    }
  }
  IE_n <- nrow(ie_match)
  cat(sprintf("  I/E observations: %d (%s)%s\n", IE_n,
              if(is_shore) "crabber-hours" else "boat trips -> tau",
              if(IE_n > 0) sprintf("  range %.0f-%.0f",
                                   min(ie_match$ie_obs), max(ie_match$ie_obs)) else ""))

  intA <- int_d |> filter(!is.na(number_of_gear), number_of_gear>0, angler_count>0)

  # v5.3: For boat fits, the CPUE denominator is gear-hours, not crabber-hours.
  # gear_time_total is computed in load_creel_data() from gear_hours (preferred)
  # or hours_fished * number_of_gear (fallback). Drop boat interviews lacking it
  # so the IntC likelihood (c ~ NegBin(lambda_C * h, r_C)) is well-defined.
  # P1: generic CPUE-denominator filter. The likelihood needs a finite, positive
  # h[a] for every interview, whatever effort unit is selected. Previously only
  # the boat had a filter, so a shore run under a non-default unit would have
  # pushed NA into stan_data.
  #
  # loo_effort_unit_comparison = TRUE additionally requires EVERY candidate
  # denominator to be valid, so the retained c[a] are identical across units and
  # elpd_loo is a valid model comparison. Set it TRUE for the comparison runs.
  # For 2024-25 shore this costs 2 of 2,869 interviews.
  h_cols_req <- if (isTRUE(params$loo_effort_unit_comparison))
                  bss_effort_h_candidates(is_shore) else eff_spec$h_col
  n_before_h <- nrow(int_d)
  for (hc in h_cols_req) {
    if (!hc %in% names(int_d)) stop("prep_bss_crab(): interview column '", hc,
                                    "' is missing; required for effort unit '",
                                    eff_spec$unit, "'.", call. = FALSE)
    v <- suppressWarnings(as.numeric(int_d[[hc]]))
    int_d <- int_d[is.finite(v) & v > 0, , drop = FALSE]
  }
  cat(sprintf("  CPUE denominator filter [%s%s]: %d -> %d interviews\n",
              eff_spec$unit,
              if (isTRUE(params$loo_effort_unit_comparison)) "; LOO-common subset" else "",
              n_before_h, nrow(int_d)))
  intA <- intA |> filter(interview_id %in% int_d$interview_id)

  # v5.1: bss_max_interviews = NULL means use full dataset (Issue 8)
  if(!is.null(params$bss_max_interviews) && nrow(int_d)>params$bss_max_interviews) {
    set.seed(42)
    int_d <- int_d |> slice_sample(n=params$bss_max_interviews)
    intA <- intA |> filter(interview_id %in% int_d$interview_id)
  }

  # v5.3: is_shore now defined earlier in the function (see comment above).

  # =========================================================================
  # v5.1 GEAR-TYPE CLASSIFICATION: Word-boundary regex + weighted assignment
  # =========================================================================
  # Detect each gear type with word-boundary patterns to avoid false matches
  # (e.g., "hotpot" would not match \bpot\b). The iForm uses controlled

  # vocabulary, but word boundaries make this robust to future data changes.

  int_d <- int_d |>
    mutate(
      has_pot      = as.integer(str_detect(gear_type, "(?i)\\bpot\\b") &
                                !str_detect(gear_type, "(?i)\\bslip\\s*ring\\b")),
      has_ring_net = as.integer(str_detect(gear_type, "(?i)\\bring\\s*net\\b")),
      has_trap     = as.integer(str_detect(gear_type, "(?i)\\b(trap|star)\\b")),
      has_snare    = as.integer(str_detect(gear_type, "(?i)\\bsnare\\b")),
      # Number of gear types reported per interview
      n_gear_types_reported = has_pot + has_ring_net + has_trap + has_snare,
      n_gear_types_reported = pmax(n_gear_types_reported, 1L)  # fallback for unrecognized
    )

  # Compute effective sample sizes (sum of fractional weights per gear type)
  gear_effective_n <- c(
    Pot       = sum(int_d$has_pot / int_d$n_gear_types_reported),
    `Ring Net` = sum(int_d$has_ring_net / int_d$n_gear_types_reported),
    Trap      = sum(int_d$has_trap / int_d$n_gear_types_reported),
    Snare     = sum(int_d$has_snare / int_d$n_gear_types_reported)
  )

  # =========================================================================
  # REGULATORY EXCLUSION: Remove gear types prohibited in this sub-season.
  # During the ring-net sub-season (Sep 16 – Nov 30), pots are illegal.
  # Any "Pot" detections are recording artifacts from multi-gear crabbers
  # naming gear they own rather than gear actually deployed. Attempting to
  # fit a Pot CPUE AR(1) process with <6 fractional interviews causes
  # catastrophic sampler failure: 800+ divergent transitions, overflow in
  # C_sum_gear, and contamination of other gear-type parameters.
  # =========================================================================
  if(length(gear_exclude) > 0) {
    excluded <- intersect(names(gear_effective_n), gear_exclude)
    if(length(excluded) > 0) {
      gear_effective_n[excluded] <- 0
      # Zero out the has_ flags for excluded types so their fractional
      # weight redistributes to remaining gear types in the weight matrix.
      # E.g., "Pot, Ring Net" becomes "Ring Net" (weight 1.0 instead of 0.5).
      has_columns_excl <- c(Pot="has_pot", `Ring Net`="has_ring_net",
                            Trap="has_trap", Snare="has_snare")
      for(gt in excluded) {
        if(gt %in% names(has_columns_excl)) {
          int_d[[has_columns_excl[gt]]] <- 0L
        }
      }
      # Recompute n_gear_types_reported after exclusion
      int_d <- int_d |>
        mutate(n_gear_types_reported = pmax(has_pot + has_ring_net + has_trap + has_snare, 1L))
      cat(sprintf("  Regulatory exclusion: %s removed (prohibited this sub-season)\n",
                  paste(excluded, collapse = ", ")))
    }
  }

  # Gear types with sufficient effective interviews are modeled independently.
  # Threshold is set conservatively because fractional weights from multi-gear
  # interviews can inflate effective N beyond what the data truly supports.
  min_gear_n <- params$bss_min_gear_effective_n %||% 15
  gear_types_present <- names(gear_effective_n[gear_effective_n >= min_gear_n])

  # If fewer than 2 gear types qualify, collapse to a single "All" category
  if(length(gear_types_present) < 2) {
    gear_type_labels <- "All"
    G_gear <- 1L
    gear_weight_matrix <- matrix(1.0, nrow = nrow(int_d), ncol = 1)
    colnames(gear_weight_matrix) <- "All"
    cat("  Insufficient gear type diversity — using single 'All' category\n")
  } else {
    gear_type_labels <- sort(gear_types_present)
    G_gear <- length(gear_type_labels)

    # Build weight matrix: each interview gets equal weight across its
    # reported gear types (e.g., "Pot, Ring Net" → 0.5 Pot + 0.5 Ring Net)
    has_columns <- c(Pot="has_pot", `Ring Net`="has_ring_net",
                     Trap="has_trap", Snare="has_snare")
    weight_cols <- list()
    for(gt in gear_type_labels) {
      weight_cols[[gt]] <- int_d[[has_columns[gt]]] / int_d$n_gear_types_reported
    }
    gear_weight_matrix <- do.call(cbind, weight_cols)

    # Interviews with no weight in any modeled type → assign to most common
    row_sums <- rowSums(gear_weight_matrix)
    zero_rows <- (row_sums == 0)
    if(any(zero_rows)) {
      most_common_idx <- which.max(colSums(gear_weight_matrix))
      gear_weight_matrix[zero_rows, most_common_idx] <- 1.0
      row_sums[zero_rows] <- 1.0
    }
    # Normalize rows to sum to 1
    gear_weight_matrix <- gear_weight_matrix / rowSums(gear_weight_matrix)
  }

  # --- Option B: the Stan model runs with G = 1 ------------------------------
  # crab_bss_gear_resolved.stan supports G > 1 structurally, but its only effort
  # observation touches g = 1 while E_sum / C_expected_sum sum over all g, so
  # raising G would add effort processes identified by their priors alone. Until
  # the gear-resolved reconstruction (Option A: effort shares in O[d,s,g] from
  # pi_gear_data, plus a multi-gear interview rule), collapse the Stan-facing gear
  # dimension to a single "All" category. The gear-composition diagnostics above
  # are retained for reporting.
  if(G_gear > 1) {
    cat(sprintf("  NOTE: %d gear types qualify (%s), but this Stan model runs with\n",
                G_gear, paste(gear_type_labels, collapse=", ")))
    cat("        G = 1. Collapsing to a single 'All' CPUE process. Gear-type catch\n")
    cat("        is reported by PE apportionment, not by the BSS. See Option A.\n")
    gear_type_labels <- "All"
    G_gear <- 1L
    gear_weight_matrix <- matrix(1.0, nrow = nrow(int_d), ncol = 1)
    colnames(gear_weight_matrix) <- "All"
  }

  cat(sprintf("  Gear types (%d): %s\n", G_gear, paste(gear_type_labels, collapse=", ")))
  cat(sprintf("  Effective N per gear type: %s\n",
    paste(sapply(1:G_gear, function(g) sprintf("%s=%.1f", gear_type_labels[g],
      colSums(gear_weight_matrix)[g])), collapse=", ")))

  # Report multi-gear interview fraction
  n_multi <- sum(int_d$n_gear_types_reported > 1)
  cat(sprintf("  Multi-gear interviews: %d of %d (%.0f%%)\n",
              n_multi, nrow(int_d), n_multi/max(nrow(int_d),1)*100))

  # =========================================================================
  # v5.1 EMPIRICAL pi_gear_data[P_n, 3, G_gear] with Laplace smoothing
  # =========================================================================
  # pi_gear is now DATA (not a Dirichlet parameter). Computed from interview
  # gear weights with Laplace smoothing (alpha=1) for each period × day_type.
  # Falls back to period-level, then sub-season-level if no interviews.

  n_day_types <- 3L  # 1=weekday, 2=weekend, 3=holiday
  alpha_smooth <- 1.0  # Laplace smoothing constant

  # Compute BSS period for each interview
  int_d$bss_period <- pvec[int_d$day_index]

  # Sub-season-level gear proportions (ultimate fallback)
  season_gear_sums <- colSums(gear_weight_matrix)
  season_pi <- (season_gear_sums + alpha_smooth) / (sum(season_gear_sums) + alpha_smooth * G_gear)

  pi_gear_data <- array(NA_real_, dim = c(P_n, n_day_types, G_gear))
  # v5.2 Fix B: Raw weighted counts for Dirichlet uncertainty propagation
  n_weighted_gear <- array(0, dim = c(P_n, n_day_types, G_gear))

  for(p in 1:P_n) {
    # Period-level proportions (intermediate fallback)
    mask_p <- (int_d$bss_period == p)
    if(any(mask_p)) {
      period_sums <- colSums(gear_weight_matrix[mask_p, , drop=FALSE])
      period_pi <- (period_sums + alpha_smooth) / (sum(period_sums) + alpha_smooth * G_gear)
    } else {
      period_pi <- season_pi
      period_sums <- season_gear_sums  # fallback raw counts for Dirichlet
    }

    for(dt in 1:n_day_types) {
      mask_pd <- mask_p & (int_d$day_type_idx == dt)
      if(any(mask_pd)) {
        dt_sums <- colSums(gear_weight_matrix[mask_pd, , drop=FALSE])
        pi_gear_data[p, dt, ] <- (dt_sums + alpha_smooth) / (sum(dt_sums) + alpha_smooth * G_gear)
        n_weighted_gear[p, dt, ] <- dt_sums  # raw counts for Stan Dirichlet
      } else {
        # Fallback: use period-level proportions and counts
        pi_gear_data[p, dt, ] <- period_pi
        n_weighted_gear[p, dt, ] <- if(any(mask_p)) period_sums else season_gear_sums
      }
    }
  }

  # --- Day type indices for Stan ---
  day_type_idx_vec <- days$day_type_idx

  # --- Sparse effort observation count ---
  n_gear <- if(is_shore) nrow(eff_d) else 0L
  n_trailer <- if(!is_shore) nrow(eff_d) else 0L
  n_effort_obs <- n_gear + n_trailer

  cat(sprintf("  Effort obs (NB-marginalized): %d (Gear_n + T_n)\n", n_effort_obs))

  # v5.3: Defensive check. Stan declares Gear_A_boat as int<lower=1>.
  # number_of_gear is stored as numeric (load_creel_data line ~305), so a
  # fractional value (e.g. 0.5) would pass the intA filter (> 0) but cast
  # to 0L and crash the Stan validator. Whole-number gear counts should
  # always hold in practice, but verify before building the data list.
  if(!is_shore && nrow(intA) > 0) {
    if(any(as.integer(intA$number_of_gear) < 1L)) {
      stop("prep_bss_crab: Gear_A_boat contains values < 1 after as.integer(). ",
           "Inspect intA$number_of_gear for non-integer or sub-1 entries.")
    }
  }

  stan_data <- list(
    D=D, G=G, S=S, P_n=P_n, period=pvec,
    w=days$day_type_num_weekend, holiday=days$day_type_num_holiday,
    # v5.3: Pots fish continuously while the trailer sits at the ramp,
    # so the daily fishing window for boats is 24 h. Shore crabbers
    # actively tend their gear, so the window remains day_length.
    # 5a/5b/F2/P1: per-day effort expansion factor L, a PARAMETER in both
    # populations (estimate_L = 1) with lognormal prior sd L_prior_sigma so its
    # uncertainty propagates. effort_scale_gear converts lambda_E (crabbers for
    # shore, gear for the boat) into the unit of h. See .effort_unit below.
    L_data        = eff_spec$L_data,
    estimate_L    = 1L,
    L_prior_sigma = eff_spec$L_prior_sigma,
    effort_scale_gear = eff_spec$effort_scale_gear,

    # F2: I/E observation stream. shore obs = crabber-hours, boat obs = trips.
    IE_n           = IE_n,
    day_IE         = if(IE_n > 0) ie_match$day_index else integer(0),
    section_IE     = if(IE_n > 0) rep(1L, IE_n) else integer(0),
    IE_obs         = if(IE_n > 0) ie_match$ie_obs else numeric(0),
    ie_group_scale = ie_group_scale,
    O=array(1.0, dim=c(D,S,G)),

    # Sparse effort observation count
    n_effort_obs = n_effort_obs,


    # Gear counts (shore only)
    Gear_n = n_gear,
    day_Gear = if(is_shore) eff_d$day_index else integer(0),
    section_Gear = if(is_shore) rep(1L,nrow(eff_d)) else integer(0),
    Gear_I = if(is_shore) as.integer(eff_d$count_quantity) else integer(0),

    # Trailer counts (boat only)
    T_n = n_trailer,
    day_T = if(!is_shore) eff_d$day_index else integer(0),
    section_T = if(!is_shore) rep(1L,nrow(eff_d)) else integer(0),
    T_I = if(!is_shore) as.integer(eff_d$count_quantity) else integer(0),

    # Direct crabber counts (reserved)
    Crab_n=0L, day_Crab=integer(0), section_Crab=integer(0),
    Crab_I=integer(0), p_I_crab=1.0,

    # Interview CPUE data (with gear weights — v5.1 Issue 4)
    # v5.3: For boats, h is gear-hours to match lambda_E (gear in water).
    # For shore, h remains crabber-hours.
    IntC=nrow(int_d), day_IntC=int_d$day_index, gear_IntC=rep(1L,nrow(int_d)),
    section_IntC=rep(1L,nrow(int_d)), c=as.integer(int_d$fish_count),
    # P1/F2: CPUE denominator, matched to the effort unit (see bss_effort_spec()).
    h = eff_spec$h_fun(int_d),

    # Gear-per-crabber expansion (shore only)
    IntA_gear = if(is_shore) nrow(intA) else 0L,
    Gear_A = if(is_shore) as.integer(intA$number_of_gear) else integer(0),
    A_A_gear = if(is_shore) as.integer(intA$angler_count) else integer(0),

    # Gear-per-boat-group expansion (boat only)
    # v5.3: Replaces T_A_int / A_A_trailer. The Stan model learns R_G_boat
    # (mean gear deployed per boat group) directly from observed
    # number_of_gear via Gear_A_boat[a] ~ Poisson(R_G_boat). The intA
    # filter already enforces number_of_gear > 0 and angler_count > 0,
    # which satisfies the Stan constraint int<lower=1> Gear_A_boat[].
    IntA_trailer = if(!is_shore) nrow(intA) else 0L,
    Gear_A_boat  = if(!is_shore) as.integer(intA$number_of_gear) else integer(0),



    # R_G estimation control (v5.1 Issue 5)
    # NOTE: estimate_R_G / R_G_fixed are not declared in the Stan model; R_G is
    # always a parameter with a proper lognormal prior. Retained for Option A.
    estimate_R_G = as.integer(is_shore),
    R_G_fixed = params$R_G_fixed_boat,

    # Hyperparameters (v5.1 Issue 7: Cauchy scales tightened to 2)
    value_cauchyDF_sigma_eps_E=2, value_cauchyDF_sigma_r_E=2,
    value_cauchyDF_sigma_eps_C=2, value_cauchyDF_sigma_r_C=2,
    value_betashape_phi_E_scaled=2, value_betashape_phi_C_scaled=2,
    value_normal_sigma_B1=1,
    value_normal_sigma_B2=1,
    value_normal_sigma_B1_C=1,
    use_B1_C = as.integer(isTRUE(params$estimate_B1_C)),
    value_normal_mu_mu_C=log(0.5), value_normal_sigma_mu_C=2,
    value_normal_mu_mu_E=if(is_shore) log(25) else log(10),
    value_normal_sigma_mu_E=2,
    value_cauchyDF_sigma_mu_C=2, value_cauchyDF_sigma_mu_E=2,

    # Metadata for output (not passed to Stan)
    .gear_type_labels = gear_type_labels,
    .ar_resolution    = ar_resolution,
    .effort_unit      = eff_spec$unit,
    .h_unit           = eff_spec$unit
  )

  # F4: minimal interview frame for the CPUE-estimator and saturation diagnostics.
  # Carried as an ATTRIBUTE, not a list element: rstan::stan(data=) tolerates the
  # dot-prefixed atomic vectors above, but a data frame inside the data list is a
  # needless risk. gear_time_total is retained for boats even though it no longer
  # enters the likelihood, because the saturation exponent needs it.
  # Safe accessor: a missing column yields NULL, and as.numeric(NULL) is length 0,
  # which would make tibble() error on recycling. Return all-NA of the right
  # length instead so a missing column degrades the diagnostic, not the run.
  .col_num <- function(df, nm) {
    if (!nm %in% names(df)) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(df[[nm]]))
  }
  attr(stan_data, "cpue_data") <- tibble(
    catch           = .col_num(int_d, "fish_count"),
    h               = as.numeric(eff_spec$h_fun(int_d)),
    number_of_gear  = .col_num(int_d, "number_of_gear"),
    gear_time_total = .col_num(int_d, "gear_time_total")
  )

  # F5: write_run_level_diagnostics() reads attr(stan_data, "ar_resolution")
  # (the pooled convention). Set it as well as the dot-prefixed element so both
  # readers work.
  attr(stan_data, "ar_resolution") <- ar_resolution
  attr(stan_data, "effort_unit")   <- eff_spec$unit

  # F4: fail in seconds, not after hours of sampling. C = E * lambda_C is catch
  # only if E and h carry the same unit.
  bss_assert_effort_units(stan_data$.effort_unit, stan_data$.h_unit, population_name)
  stan_data
}
