###############################################################################
# prep_bss_crab_pooled.R  (pooled-CPUE driver)
#
# Build the Stan data list for crab_bss_pooled.stan (adaptive AR resolution, data-
# driven R_G prior, L_effective, effort-unit spec via bss_effort_spec). Extracted
# from the pooled driver and renamed prep_bss_crab_pooled (distinct from the gear-
# resolved prep_bss_crab_gear). Pure given its arguments.
###############################################################################

prep_bss_crab_pooled <- function(days, summ, est_catch_group, params, population_name,
                          gear_regime = NULL, ie_data = NULL) {

  eff <- summ$effort_index |> filter(count_sequence <= params$bss_max_count_seq)
  D <- nrow(days); G <- 1L; S <- 1L

  eff_d <- eff |> left_join(days |> select(event_date,day_index), by="event_date") |> filter(!is.na(day_index))

  is_shore <- (population_name == "shore")

  # --- P1/POOL-3/POOL-7: effort-unit specification (shared with run_pe via the module
  # 03_R_functions/bss_effort_spec.R, so the BSS and PE can never land on different
  # scales). As of v7.7 BOTH components are on gear-deployments (h = number_of_gear):
  #   shore -> E_scale = R_G (lambda_E is crabbers -> gear), L = tau_shore   (POOL-7)
  #   boat  -> E_scale = 1   (lambda_E is already gear),     L = tau_boat    (POOL-1/3)
  # replacing the old crabber-hours (shore) and flat gear-hours L = 24 (boat) scales.
  # Set params$shore_effort_unit = "crabber-hours" to revert shore only.
  eff_spec <- bss_effort_spec(is_shore, days, params)
  cat(sprintf("  Effort unit: %s (effort_scale_gear=%d)\n",
              eff_spec$unit, eff_spec$effort_scale_gear))

  # --- Adaptive AR temporal resolution (POOL-6: shared selector) ------------
  # Resolution selection (the finest daily/weekly/monthly the effort series
  # supports), the per-population cap (params$ar_max_resolution), the
  # params$ar_force experiment override, and the day -> AR-period mapping now
  # live in the shared module 03_R_functions/bss_ar_resolution.R
  # (bss_select_ar_resolution), so the pooled and gear-resolved tracks use ONE
  # selector and cannot drift. Pooled is the DATA-DRIVEN caller:
  # fixed_resolution = NULL means 'pick from data density, then apply the cap'.
  # gear_regime lets the cap differ by sub-season (e.g. the pooled shore caps its
  # thin pot-closure fit at biweekly while the all-gear fit stays data-driven at
  # daily); NULL preserves the original per-population cap behavior.
  ar_sel        <- bss_select_ar_resolution(days, eff_d, population_name, params,
                                             fixed_resolution = NULL,
                                             gear_regime = gear_regime)
  ar_resolution <- ar_sel$resolution
  P_n           <- ar_sel$P_n
  pvec          <- ar_sel$pvec

  # --- Interview data (all interviews used, no subsampling) ---
  int_cg <- summ$interview |>
    mutate(fish_count = .data[[est_catch_group]])

  # P1/POOL-3/POOL-7: the CPUE denominator is the effort unit's h column (boat and,
  # as of v7.7, shore: number_of_gear = deployments; shore reverts to
  # fishing_time_total = crabber-hours if shore_effort_unit is changed back). Drop
  # interviews lacking a finite positive h so c[a] ~ NB2(lambda_C * h, r_C) is
  # well-defined, whatever effort unit is selected.
  h_vec_name <- eff_spec$h_col
  .hv <- suppressWarnings(as.numeric(int_cg[[h_vec_name]]))
  int_cg <- int_cg[is.finite(.hv) & .hv > 0, , drop = FALSE]
  cat(sprintf("  CPUE denominator: %s (%s), %d interviews with valid h\n",
              h_vec_name, eff_spec$unit, nrow(int_cg)))

  int_d <- int_cg |> left_join(days |> select(event_date,day_index), by="event_date") |>
    filter(!is.na(day_index)) |> distinct(interview_id, .keep_all=TRUE)

  # POOL-2: incomplete-trip filter for the BSS CPUE likelihood. Incomplete trips
  # (soak-time gear not yet retrieved) read ~-20% low for pots/traps, biasing the
  # shore CPUE, hence the publication shore estimate, low. Gear-resolved filters
  # these; pooled previously did not (trip_status was computed at L504 but never
  # used). Toggle with params$filter_incomplete_trips (default TRUE). Missing
  # status is KEPT (a blank completed_trip may be complete; dropping it loses data).
  # The filter also propagates to intA below (the gear-per-group / trailer set),
  # matching the gear-resolved structure.
  n_before_trip_filter <- nrow(int_d)
  if(isTRUE(params$filter_incomplete_trips)) {
    int_d <- int_d |> filter(trip_status == "Complete" | is.na(trip_status))
    cat(sprintf("  Incomplete-trip filter ON: %d -> %d interviews (removed %d)\n",
                n_before_trip_filter, nrow(int_d), n_before_trip_filter - nrow(int_d)))
  } else {
    cat("  Incomplete-trip filter OFF (retaining all trips)\n")
  }

  intA <- int_d |> filter(!is.na(number_of_gear), number_of_gear>0, angler_count>0)

  # POOL-1: Stan declares Gear_A_boat as int<lower=1>, fed from as.integer(number_of_gear).
  # number_of_gear is whole pots, but guard defensively: drop any boat interview whose
  # integer gear count is < 1 so the Gear_A_boat[a] ~ poisson(R_G_boat) term never
  # receives a 0 (which aborts the fit). Never triggers on clean data.
  if(!is_shore && nrow(intA) > 0) {
    .ng_int <- as.integer(suppressWarnings(as.numeric(intA$number_of_gear)))
    n_bad <- sum(!is.finite(.ng_int) | .ng_int < 1)
    if(n_bad > 0) {
      cat(sprintf("  POOL-1: dropping %d boat interview(s) with gear count < 1\n", n_bad))
      intA <- intA[is.finite(.ng_int) & .ng_int >= 1, , drop = FALSE]
    }
  }

  cat(sprintf("  Using ALL %d interviews (no subsampling)\n", nrow(int_d)))

  # --- Effort expansion factor L (from the effort spec) ---
  # P1/POOL-3: L is a PARAMETER in BOTH populations now (estimate_L = 1). Shore: the
  # I/E-derived effective day length in hours. Boat: tau_boat, the gear-deployment
  # turnover (~1.2), replacing the old flat L = 24 gear-hours (POOL-3). Its prior SD
  # comes from the effort spec so uncertainty propagates into the boat total.
  L_data_vec      <- eff_spec$L_data
  L_sigma_vec     <- eff_spec$L_prior_sigma
  estimate_L_flag <- 1L
  cat(sprintf("  L (%s): range [%.2f, %.2f], sigma [%.2f, %.2f], estimate_L=1\n",
              if(is_shore) "effective day length, hours" else "tau_boat, deployment turnover",
              min(L_data_vec), max(L_data_vec), min(L_sigma_vec), max(L_sigma_vec)))

  # --- Sparse effort observation counts ---
  n_gear <- if(is_shore) nrow(eff_d) else 0L
  n_trailer <- if(!is_shore) nrow(eff_d) else 0L
  n_effort_obs <- n_gear + n_trailer
  cat(sprintf("  Sparse eps_E_H: %d effort obs\n", n_effort_obs))

  # --- I/E observations (shore only) ---
  ie_match <- tibble(event_date = Date(), ie_crabber_hours = numeric())
  # The I/E stream is fed for SHORE only (predicted mean lambda_E * L = crabber-
  # hours). With POOL-1/POOL-3 the boat is now on the gear-deployment scale, so a
  # boat I/E observation would be boat TRIPS with predicted mean
  # (lambda_E / R_G_boat) * tau (see gear-resolved's ie_group_scale). Activating
  # that is a follow-up; the boat I/E is empty for the 2024-25 window anyway (no
  # WBL ingress days inside it), so leaving it off here is behavior-neutral.
  if(!is.null(ie_data) && nrow(ie_data) > 0 && is_shore) {
    ie_match <- ie_data |>
      filter(population == "shore") |>
      inner_join(days |> select(event_date, day_index), by = "event_date") |>
      filter(!is.na(day_index), ie_crabber_hours > 0)
  }

  IE_n <- nrow(ie_match)
  if(IE_n > 0) {
    cat(sprintf("  I/E observations: %d days (crabber-hrs range: %.0f-%.0f)\n",
                IE_n, min(ie_match$ie_crabber_hours), max(ie_match$ie_crabber_hours)))
  } else {
    cat("  I/E observations: 0\n")
  }

  # GR-8 (2026-07-13): drop the shore I/E stream when it has fewer than
  # params$ie_min_obs_shore in-window observations. With only 1-2 I/E days the lognormal
  # I/E likelihood barely informs the effective day length L, but it lets sigma_IE
  # (exponential(5), mode at 0) shrink and stiffen against those few points, a funnel
  # that drove divergences in the sparse shore pot-closure fit. Dropping the stream
  # leaves sigma_IE prior-only (decoupled, like the boat), removing the funnel at
  # negligible information loss. Set ie_min_obs_shore = 0 to disable. The prior is left
  # as exponential(5) ON PURPOSE: tightening it (e.g. lognormal(log(0.3), 0.5)) would
  # push the shore all-gear sigma_IE (~1.07) down and force possibly-unrepresentative
  # I/E days to bind harder; see the shore-I/E representativeness diagnostic (GR-9 / item 4).
  ie_min_obs_shore <- params$ie_min_obs_shore %||% 3L
  if (IE_n > 0 && IE_n < ie_min_obs_shore) {
    cat(sprintf("  I/E stream dropped: %d in-window obs < ie_min_obs_shore = %d; sigma_IE left decoupled (GR-8 guard).\n",
                IE_n, ie_min_obs_shore))
    IE_n <- 0L
  }

  # --- R_G prior (data-driven by default; overridable for prior sensitivity) ---
  # T1.3 (2026-07-12): params$R_G_prior_mu / R_G_prior_sigma override the data-driven
  # empirical R_G so the prior-sensitivity sweep (backlog T1.3 / critique 2) is a config
  # change, not a code edit. Leave both unset for production (the data-driven value). To
  # sweep, set params$R_G_prior_mu to 1.0, the empirical value (~1.28), and 1.5 in turn
  # and compare the port totals (a tighter R_G_prior_sigma makes the prior bind harder).
  R_G_empirical   <- params$R_G_prior_mu %||% summ$empirical_R_G
  R_G_prior_sigma <- params$R_G_prior_sigma %||% 0.3
  .rg_src <- if (is.null(params$R_G_prior_mu)) "data-driven" else "OVERRIDE (prior sensitivity)"
  cat(sprintf("  R_G prior = Lognormal(log(%.2f), %.2f) [%s]\n",
              R_G_empirical, R_G_prior_sigma, .rg_src))

  # Weakly-informative level priors (sigma 2, so the data dominate). Boat lambda_E
  # is gear in the water on the deployment scale; lambda_C is crab per deployment.
  # Matches crab_bss_gear_resolved.stan's boat priors.
  mu_E_prior <- if(is_shore) log(25) else log(10)
  mu_C_prior <- log(0.5)

  # POOL-4: resolve the collapse_mu_hier lever (global logical OR per-population
  # named list; default FALSE = current v6.8 hierarchy). Passed to Stan as an int.
  collapse_flag <- if(is.list(params$collapse_mu_hier)) {
    isTRUE(params$collapse_mu_hier[[population_name]])
  } else {
    isTRUE(params$collapse_mu_hier)
  }

  # item 1: razor-dig shore EFFORT covariate. Non-zero only for a shore fit when the razor
  # term is active (params$razor_dig_active, set by the driver from razor_dig_mode and, in
  # "auto" mode, the spillover diagnostic p-value); otherwise all zeros, so B3 stays
  # decoupled (prior-only, like sigma_IE at IE_n = 0). Dig dates come from the consolidated
  # opener calendar (params$razor_dig_dates, the nearby-beach dig days).
  razor_vec <- if (is_shore && isTRUE(params$razor_dig_active)) {
    as.numeric(days$event_date %in% (params$razor_dig_dates %||% as.Date(character(0))))
  } else rep(0, D)
  if (is_shore && isTRUE(params$razor_dig_active))
    cat(sprintf("  razor_dig effort term ACTIVE (%d razor days in window)\n", sum(razor_vec)))

  stan_data <- list(
    D=D, G=G, S=S,
    P_n=P_n, period=pvec,
    w=days$day_type_num_weekend, holiday=days$day_type_num_holiday,
    razor=razor_vec,
    O=array(1.0, dim=c(D,S,G)),
    collapse_mu_hier = as.integer(collapse_flag),   # POOL-4 lever (0 = v6.8 default)

    L_data = L_data_vec,
    estimate_L = estimate_L_flag,
    L_prior_sigma = L_sigma_vec,
    effort_scale_gear = as.integer(eff_spec$effort_scale_gear),  # P1/POOL-3

    n_effort_obs = n_effort_obs,

    Gear_n = n_gear,
    day_Gear = if(is_shore) eff_d$day_index else integer(0),
    section_Gear = if(is_shore) rep(1L,nrow(eff_d)) else integer(0),
    Gear_I = if(is_shore) as.integer(eff_d$count_quantity) else integer(0),

    T_n = n_trailer,
    day_T = if(!is_shore) eff_d$day_index else integer(0),
    section_T = if(!is_shore) rep(1L,nrow(eff_d)) else integer(0),
    T_I = if(!is_shore) as.integer(eff_d$count_quantity) else integer(0),

    Crab_n=0L, day_Crab=integer(0), section_Crab=integer(0),
    Crab_I=integer(0), p_I_crab=1.0,

    IntC=nrow(int_d), day_IntC=int_d$day_index, gear_IntC=rep(1L,nrow(int_d)),
    section_IntC=rep(1L,nrow(int_d)), c=as.integer(int_d$fish_count),
    # P1/POOL-3: CPUE denominator matched to the effort unit (see bss_effort_spec()).
    h = as.numeric(eff_spec$h_fun(int_d)),

    IntA_gear = if(is_shore) nrow(intA) else 0L,
    Gear_A = if(is_shore) as.integer(intA$number_of_gear) else integer(0),
    A_A_gear = if(is_shore) as.integer(intA$angler_count) else integer(0),

    # POOL-1: gear-per-boat-group interviews. R_G_boat is learned from observed
    # number_of_gear via Gear_A_boat[a] ~ poisson(R_G_boat). Replaces the degenerate
    # T_A_int (all ones, so bernoulli pinned R_T at 1) / A_A_trailer. intA is already
    # filtered to number_of_gear > 0, satisfying Stan's int<lower=1> Gear_A_boat[].
    IntA_trailer = if(!is_shore) nrow(intA) else 0L,
    Gear_A_boat  = if(!is_shore) as.integer(intA$number_of_gear) else integer(0),

    IE_n = IE_n,
    day_IE = if(IE_n > 0) ie_match$day_index else integer(0),
    section_IE = if(IE_n > 0) rep(1L, IE_n) else integer(0),
    IE_crabber_hours = if(IE_n > 0) ie_match$ie_crabber_hours else numeric(0),

    # Tightened Cauchy scales
    value_cauchyDF_sigma_eps_E=1, value_cauchyDF_sigma_r_E=1,
    value_cauchyDF_sigma_eps_C=1, value_cauchyDF_sigma_r_C=1,
    value_betashape_phi_E_scaled=2, value_betashape_phi_C_scaled=2,
    value_normal_sigma_B1=1, value_normal_sigma_B2=1,
    value_normal_sigma_B1_C=1,
    value_normal_sigma_B2_C=1,   # item 6a (holiday CPUE effect)
    value_normal_sigma_B3=1,     # item 1 (razor-dig effort effect)
    estimate_cpue_density=as.integer(isTRUE(params$estimate_cpue_density)),  # item 6b (off by default)
    log_E_ref=mu_E_prior,        # item 6b: center the density covariate at the effort level prior
    value_normal_mu_mu_C=mu_C_prior, value_normal_sigma_mu_C=2,
    value_normal_mu_mu_E=mu_E_prior, value_normal_sigma_mu_E=2,
    value_cauchyDF_sigma_mu_C=1, value_cauchyDF_sigma_mu_E=1,

    R_G_prior_mu = R_G_empirical,
    R_G_prior_sigma = R_G_prior_sigma
    # POOL-1: R_T_alpha / R_T_beta removed. R_G_boat carries a fixed lognormal prior
    # in the Stan model (log(4), 0.5), matching crab_bss_gear_resolved.stan.
  )

  # Store AR resolution for downstream reporting
  attr(stan_data, "ar_resolution") <- ar_resolution

  # POOL-5: attach the per-interview CPUE data + effort-unit tags so the shared CPUE
  # diagnostics (03_R_functions/bss_cpue_diagnostics.R) run per fit. The tags feed
  # bss_assert_effort_units() (effort E and the CPUE denominator h must share a unit).
  # With POOL-1/POOL-3 the boat is now gear-deployments/gear-deployments and shore is
  # crabber-hours/crabber-hours, both self-consistent, so the assertion passes; the
  # saturation/linearity CSVs then confirm the deployment scale is valid for pots
  # (boat catch is flat in soak HOURS but linear in deployments). cpue_data is an
  # ATTRIBUTE, not a stan_data list entry, so rstan never processes it (GR-13); the
  # two dot-prefixed unit tags are plain scalars, tolerated by rstan::stan().
  stan_data[[".effort_unit"]] <- eff_spec$unit   # E = lambda_E * E_scale * L, same unit as h
  stan_data[[".h_unit"]]      <- eff_spec$unit
  # gear_time_total is retained for the saturation diagnostic even though the boat
  # no longer uses it as the CPUE denominator (h is now number_of_gear).
  attr(stan_data, "cpue_data") <- tibble(
    catch           = as.numeric(int_d$fish_count),
    h               = as.numeric(eff_spec$h_fun(int_d)),
    number_of_gear  = suppressWarnings(as.numeric(int_d[["number_of_gear"]]  %||% rep(NA_real_, nrow(int_d)))),
    gear_time_total = suppressWarnings(as.numeric(int_d[["gear_time_total"]] %||% rep(NA_real_, nrow(int_d))))
  )

  stan_data
}
