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
# bss_day_length.R
#
# Shared shore day-length module for the crab BSS models. Both the pooled and
# gear-resolved drivers auto-source this folder, so both use ONE implementation
# and cannot drift.
#
# WHY SHORE DAY LENGTH IS NOT CIVIL TWILIGHT
#   Shore effort is counted at a peak instant and expanded by a day length L.
#   Civil twilight (dawn to dusk) at 46.9N runs 9 to 17 hours, but ingress/egress
#   surveys show the EFFECTIVE shore day, total crabber-hours divided by peak
#   crabbers present, averages only 3.5 to 5.0 hours. Using civil twilight
#   therefore overestimates shore effort by roughly 2x. L_effective corrects it.
#
#   Boats do NOT use a day length at all. STALE COMMENT CORRECTED 2026-08-25: the old
#   text here said "L = 24 by construction (the gear-hours formulation)", which stopped
#   being true at POOL-3 / v7.6 when the boat moved to gear-deployments. The boat's L is
#   tau_boat, the deployment turnover (~1.2), set by bss_effort_spec().
#
# WHAT PRODUCTION ACTUALLY EXPANDS ON (improvement 2, 2026-08-25)
#   Since the v7.7 shore unit move, SHORE also expands on a turnover, not on a day
#   length: E = lambda_E * R_G * tau_shore. Both quantities come from the same 15-minute
#   I/E presence series:
#       L_effective = crabber-hours / peak present   (~5.3 h mean)
#       turnover    = arrivals      / peak present   (~1.72 over 30 WDF20 days)
#   and their ratio is the implied trip length, 5.26 / 1.72 = 3.07 h against an interview
#   mean trip length of 3.23 h, which is the free consistency check on the method.
#   L_effective is still computed every run (it sets the day_length column the
#   diagnostics and the civil-twilight comparison use) but it is NOT on the estimation
#   path unless shore_effort_unit is set back to a time unit. The fallback ladder below
#   therefore governs a DIAGNOSTIC quantity in the production configuration.
#
#   Consequence for the I/E likelihood: because the predicted quantity is
#   lambda_E * tau = TRIPS, the observation must be the arrival count, not crabber-hours.
#   That pairing is now owned by bss_effort_spec()$ie_obs_col; see the note there.
#
# THE FALLBACK LADDER (automatic, no toggle)
#   estimate_L_effective() degrades WITHIN the same estimand rather than
#   switching estimands on a sample-size threshold:
#
#     n_ie >= params$ie_min_obs_for_regression   -> quadratic regression of
#           log(L_effective) on day-of-year + day type, with per-day prediction
#           uncertainty (L_sigma).
#     3 <= n_ie <  params$ie_min_obs_for_regression -> grand mean of
#           log(L_effective) with empirical SD.
#     n_ie <  3                                  -> returns NULL.
#
#   bss_assign_day_length() then falls back to civil twilight ONLY when the model
#   is NULL, i.e. when there is effectively no I/E data at all. It warns loudly
#   when it does, because that path reintroduces the ~2x shore bias.
#
#   Note the ladder's middle rung is what avoids a discontinuity: a season with
#   4 I/E surveys and one with 5 both estimate an effective day length (~3.5-5 h),
#   they differ only in how that estimate is smoothed. They do not jump between
#   3.5 h and 16 h.
#
# CAP
#   Civil twilight is clamped to [day_length_min_hours, day_length_max_hours],
#   default [9, 17]. 17 (not 16) because civil twilight at 46.9N reaches ~17 h at
#   the summer solstice; a 16 h cap slightly underestimated peak-season effort.
#   The cap only matters on the civil-twilight fallback rung, since whenever an
#   L_effective model exists day_length is overwritten by its prediction.
#
# CONTENTS
#   fetch_ie_data(params)                          lifted verbatim from the pooled
#   estimate_L_effective(ie_data, params)            driver (v7.4); the dead pot_open_date
#                                                    argument was removed 2026-09-12
#   bss_day_length_civil(dates, params)            civil-twilight helper
#   bss_assign_day_length(days, L_eff_model, params)  sets day_length, L_mu,
#                                                  L_prior_sigma on a days tibble
#
# Requires: dplyr, tibble, readxl, here, suncalc, lubridate (yday).
###############################################################################


# ===========================================================================
# 1.5 INGRESS/EGRESS DATA
# ===========================================================================

fetch_ie_data <- function(params) {
  cat("\n  Reading I/E data...\n")
  ie_file <- here("04_input_files", params$ie_data_file)
  if(!file.exists(ie_file)) {
    cat("  WARNING: I/E file not found at", ie_file, "\n")
    return(tibble(
      event_date = Date(), location_name = character(), population = character(),
      ie_crabber_hours = numeric(), n_intervals = integer(), max_present = numeric(),
      survey_hours = numeric(), total_arrivals = numeric(), total_departures = numeric(),
      mean_present = numeric(), L_effective = numeric()
    ))
  }

  ie_raw <- read_input_workbook(ie_file, sheet = params$ie_sheet) |>
    mutate(event_date = as.Date(date))

  # Optional: restrict I/E to the current fishery season. Default FALSE preserves the
  # historical pooling (the L_effective day-of-year regression intentionally uses all
  # seasons of I/E). The ingress_egress `season` column is now the fishery season label.
  if (isTRUE(params$ie_filter_by_season) && "season" %in% names(ie_raw))
    ie_raw <- ie_raw |> filter(season %in% params$season_filter)

  ie_shore <- ie_raw |>
    filter(location_name == params$ie_shore_location) |>
    mutate(
      crabbers_on = replace_na(as.numeric(crabbers_on), 0),
      crabbers_off = replace_na(as.numeric(crabbers_off), 0),
      crabber_flow = replace_na(as.numeric(crabber_flow), 0)
    ) |>
    group_by(event_date) |>
    summarise(
      location_name = first(location_name),
      ie_crabber_hours = sum(crabber_flow * 0.25),
      n_intervals = n(),
      max_present = max(crabber_flow),
      survey_hours = n() * 0.25,
      total_arrivals = sum(crabbers_on),
      total_departures = sum(crabbers_off),
      day_type = first(day_type),
      season = first(season),
      .groups = "drop"
    ) |>
    mutate(
      population = "shore",
      mean_present = ie_crabber_hours / survey_hours,
      L_effective = if_else(max_present > 0, ie_crabber_hours / max_present, 0),
      ie_trips    = total_arrivals,   # crabber arrivals; unused by the shore model
      ie_turnover = if_else(max_present > 0, total_arrivals / max_present, NA_real_)
    )

  ie_boat <- ie_raw |>
    filter(location_name == params$ie_boat_location) |>
    mutate(
      boats_in = replace_na(as.numeric(boats_in), 0),
      boats_out = replace_na(as.numeric(boats_out), 0),
      boat_flow = replace_na(as.numeric(boat_flow), 0)
    ) |>
    group_by(event_date) |>
    summarise(
      location_name = first(location_name),
      ie_crabber_hours = sum(boat_flow * 0.25),
      n_intervals = n(),
      max_present = max(boat_flow),
      survey_hours = n() * 0.25,
      total_arrivals = sum(boats_in),
      total_departures = sum(boats_out),
      day_type = first(day_type),
      season = first(season),
      .groups = "drop"
    ) |>
    mutate(
      population = "private_boat",
      mean_present = ie_crabber_hours / survey_hours,
      L_effective = if_else(max_present > 0, ie_crabber_hours / max_present, 0),
      # F2: the boat model expands on the DEPLOYMENT scale, so what it needs from
      # I/E is the number of boat trips (ingress count) and the turnover
      # tau = trips / peak boats present. ie_crabber_hours here is boat-hours and
      # is retained for diagnostics only; it is NOT the boat observation model's
      # predicted quantity (see the F2 note in crab_bss_gear_resolved.stan).
      ie_trips    = total_arrivals,
      ie_turnover = if_else(max_present > 0, total_arrivals / max_present, NA_real_)
    )

  ie_all <- bind_rows(ie_shore, ie_boat) |>
    filter(ie_crabber_hours > 0)

  # 2026-09-08 (review item 2): keep the SHORE interval rows, with their clock time, as an
  # attribute so estimate_shore_turnover() can evaluate presence at the hours the creel
  # counts were taken. The workbook gained a `time` column on 2026-09-08; without it the
  # attribute is NULL and the turnover derivation degrades to the peak-based value.
  ie_int <- NULL
  if ("time" %in% names(ie_raw)) {
    ie_int <- ie_raw |>
      filter(location_name == params$ie_shore_location) |>
      mutate(hour = .ie_hour_of(time),
             crabbers_on  = replace_na(as.numeric(crabbers_on), 0),
             crabber_flow = replace_na(as.numeric(crabber_flow), 0)) |>
      filter(is.finite(hour)) |>
      select(event_date, season, day_type, hour, crabbers_on, crabber_flow) |>
      arrange(event_date, hour)
  }
  attr(ie_all, "ie_intervals") <- ie_int
  # 2026-09-09: keep the BOAT interval rows of every boat-I/E site too (arrivals and
  # returns by clock time), for the sampler-shift coverage diagnostic
  # (03_R_functions/sampler_shifts.R): what share of a day's boat returns falls inside the
  # hours a sampler was in port to classify them.
  ie_boat_int <- NULL
  if ("time" %in% names(ie_raw) && all(c("boats_in", "boats_out") %in% names(ie_raw))) {
    ie_boat_int <- ie_raw |>
      filter(!is.na(boats_in) | !is.na(boats_out)) |>
      mutate(hour = .ie_hour_of(time),
             boats_in  = replace_na(as.numeric(boats_in), 0),
             boats_out = replace_na(as.numeric(boats_out), 0)) |>
      filter(is.finite(hour)) |>
      select(event_date, season, day_type, location_name, hour, boats_in, boats_out) |>
      arrange(location_name, event_date, hour)
  }
  attr(ie_all, "ie_boat_intervals") <- ie_boat_int

  cat(sprintf("  I/E survey days: %d shore (WDF20), %d boat (WBL)\n",
              sum(ie_all$population == "shore"),
              sum(ie_all$population == "private_boat")))
  if(nrow(ie_all) > 0)
    cat(sprintf("  Date range: %s to %s\n", min(ie_all$event_date), max(ie_all$event_date)))

  if(nrow(ie_all |> filter(population == "shore")) > 0) {
    shore_ie <- ie_all |> filter(population == "shore")
    cat(sprintf("  Shore L_effective: mean=%.1f hrs (range %.1f-%.1f)\n",
                mean(shore_ie$L_effective), min(shore_ie$L_effective), max(shore_ie$L_effective)))
  }

  # --- Phase 3: crabbing-fraction classification rows (crab-vs-total boats, per WBL row) ---
  # Optional columns on the boat I/E (WBL) rows: params$ie_crab_col / ie_total_col. Emitted
  # PER ROW so the crab-fraction helper (03_R_functions/crab_fraction.R) can aggregate by
  # stratum (month / day_type). Absent columns (the current state; the WBL classification
  # pilot is in progress) -> empty -> every stratum uses the set value. The preps/PE read
  # this via params$crab_fraction_rows (the driver lifts attr(ie_data, "crab_fraction_rows")).
  cf_crab_col  <- params$ie_crab_col  %||% "boats_crabbing"
  cf_total_col <- params$ie_total_col %||% "boats_total"
  cf_rows <- tibble(event_date = as.Date(character()),
                    boats_crabbing = numeric(), boats_total = numeric())
  wbl_raw <- ie_raw |> filter(location_name == params$ie_boat_location)
  if (nrow(wbl_raw) > 0 && all(c(cf_crab_col, cf_total_col) %in% names(wbl_raw))) {
    cf_rows <- wbl_raw |>
      transmute(event_date     = as.Date(date),
                boats_crabbing = suppressWarnings(as.numeric(.data[[cf_crab_col]])),
                boats_total    = suppressWarnings(as.numeric(.data[[cf_total_col]]))) |>
      filter(is.finite(boats_total))
    .nt <- sum(cf_rows$boats_total, na.rm = TRUE); .nc <- sum(cf_rows$boats_crabbing, na.rm = TRUE)
    cat(sprintf("  Crab-fraction I/E classification: %d WBL rows, %.0f total / %.0f crab boats (f_hat = %s)\n",
                nrow(cf_rows), .nt, .nc, if (.nt > 0) sprintf("%.2f", .nc / .nt) else "NA"))
  } else {
    cat("  Crab-fraction I/E classification: columns absent; f will use the set-value fallback.\n")
  }
  attr(ie_all, "crab_fraction_rows") <- cf_rows

  return(ie_all)
}


# ===========================================================================
# 1.6 L_EFFECTIVE MODEL - Regression with uncertainty propagation
#
# Fits a regression of log(L_effective) on day-of-year (quadratic) and day type.
# This captures the seasonal gradient WITHIN sub-seasons and provides per-day
# prediction uncertainty for propagation into the Stan model.
#
# 2026-09-12: the `pot_open_date` ARGUMENT WAS DEAD. It appeared in the signature and
# nowhere in the body, so there has never been a pots-open split in this regression; the
# only predictors are yday (quadratic) and day type. Both drivers passed
# params$pot_open_date into it, and run_config's comment on pot_open_date plus
# CHANGE_REGISTER A14 both described a multi-season approximation ("pot_open_date is a
# SINGLE date feeding the L_effective I/E regression split; with two seasons it is exact
# for one season only") that therefore never existed. The argument is removed rather than
# used, because a pots-open indicator is NOT the right fix either: the closure boundary
# is already inside the yday term, and adding a second, season-specific predictor to a
# model fitted on 40 I/E days would cost more than it buys. THE REAL multi-season
# approximation, now stated where it belongs: the regression POOLS every season's I/E
# days and assumes one yday -> L relationship across them. On the 2023-25 span that is
# 2023-24's shorter Grays Harbor shifts and 2024-25's being fitted as one curve. A
# per-season interaction is the fix if the two seasons' I/E day lengths diverge; check
# L_effective_ie_detail.csv (residuals by season) before assuming they do not.
#
# Returns: list with
#   $predict_fn: function(event_date, day_type) -> tibble(L_mu, L_sigma)
#   $model: the fitted lm object (NULL on the grand-mean rung)
#   $detail: per-I/E-day data
#   $n_obs, $method: "regression" or "grand_mean"
# ===========================================================================

estimate_L_effective <- function(ie_data, params) {
  cat("\n  Fitting L_effective regression from historical I/E data...\n")

  ie_shore <- ie_data |>
    filter(population == "shore", L_effective > 0)

  if(nrow(ie_shore) < 3) {
    cat("  WARNING: Fewer than 3 shore I/E days. Cannot fit L_effective model.\n")
    return(NULL)
  }

  # Prepare regression data
  ie_shore <- ie_shore |>
    mutate(
      yday = yday(event_date),
      day_type_group = if_else(day_type %in% c("Weekend", "weekend"), "weekend", "weekday"),
      log_L = log(L_effective)
    )

  n_ie <- nrow(ie_shore)
  cat(sprintf("  I/E observations for regression: %d\n", n_ie))

  if(n_ie >= params$ie_min_obs_for_regression) {
    # --- Fit quadratic regression on day-of-year + day type ---
    # log(L_effective) = b0 + b1*yday + b2*yday^2 + b3*weekend + error
    L_fit <- lm(log_L ~ poly(yday, 2) + day_type_group, data = ie_shore)

    cat("  L_effective regression summary:\n")
    cat(sprintf("    R² = %.3f, residual SE = %.3f (log scale)\n",
                summary(L_fit)$r.squared, sigma(L_fit)))
    cat(sprintf("    Coefficients:\n"))
    coefs <- coef(L_fit)
    for(nm in names(coefs)) {
      cat(sprintf("      %s: %.4f\n", nm, coefs[nm]))
    }

    sigma_resid <- sigma(L_fit)  # residual SD on log scale

    # Create prediction function
    predict_L <- function(event_dates, day_types) {
      newdata <- tibble(
        yday = yday(event_dates),
        day_type_group = if_else(day_types %in% c("weekend","holiday"), "weekend", "weekday")
      )
      pred <- predict(L_fit, newdata = newdata, se.fit = TRUE)

      tibble(
        L_mu = exp(pred$fit),  # predicted median on natural scale
        # Total prediction uncertainty = sqrt(regression SE² + residual variance)
        L_sigma = pmax(sqrt(pred$se.fit^2 + sigma_resid^2), 0.1),
        L_fit_only = exp(pred$fit),  # same as L_mu (for diagnostics)
        L_se_regression = pred$se.fit,
        L_se_total = pmax(sqrt(pred$se.fit^2 + sigma_resid^2), 0.1)
      )
    }

  } else {
    # --- Fallback: grand mean with empirical SD ---
    cat(sprintf("  Fewer than %d I/E obs - using grand mean fallback.\n",
                params$ie_min_obs_for_regression))
    grand_mean_log <- mean(ie_shore$log_L)
    grand_sd_log <- max(sd(ie_shore$log_L), 0.2)

    predict_L <- function(event_dates, day_types) {
      n <- length(event_dates)
      tibble(
        L_mu = rep(exp(grand_mean_log), n),
        L_sigma = rep(grand_sd_log, n),
        L_fit_only = rep(exp(grand_mean_log), n),
        L_se_regression = rep(grand_sd_log / sqrt(nrow(ie_shore)), n),
        L_se_total = rep(grand_sd_log, n)
      )
    }
    L_fit <- NULL
  }

  # Predict for I/E days (for diagnostics)
  ie_pred <- predict_L(ie_shore$event_date, ie_shore$day_type)
  ie_shore <- bind_cols(ie_shore, ie_pred |> rename(pred_L_mu = L_mu, pred_L_sigma = L_sigma))

  cat(sprintf("\n  Predicted L_effective range: %.1f-%.1f hrs\n",
              min(ie_shore$pred_L_mu), max(ie_shore$pred_L_mu)))
  cat(sprintf("  Prediction uncertainty (sigma on log scale): %.2f-%.2f\n",
              min(ie_shore$pred_L_sigma), max(ie_shore$pred_L_sigma)))

  result <- list(
    predict_fn = predict_L,
    model = L_fit,
    # 2026-08-27: carry the DEPLOYMENT-scale quantities alongside the hours-scale ones.
    # Production has expanded shore effort on a turnover since v7.7, and since 2026-08-25 the
    # shore I/E likelihood observes ie_trips, but this file reported only ie_crabber_hours and
    # L_effective (hours per crabber). That made the one file a reader would open to audit the
    # I/E unit change the one file that could not show it. turnover = arrivals / peak present
    # is the quantity `L` actually carries under gear-deployments; hours_per_trip is their
    # ratio, the method's free internal consistency check against the interview-reported
    # trip length.
    detail = ie_shore |>
      mutate(turnover      = if_else(max_present > 0, ie_trips / max_present, NA_real_),
             hours_per_trip = if_else(ie_trips > 0, ie_crabber_hours / ie_trips, NA_real_)) |>
      select(event_date, day_type, day_type_group, yday,
             ie_crabber_hours, ie_trips, max_present, L_effective, turnover, hours_per_trip,
             pred_L_mu, pred_L_sigma),
    n_obs = n_ie,
    method = if(!is.null(L_fit)) "regression" else "grand_mean"
  )

  return(result)
}


# ===========================================================================
# CIVIL TWILIGHT HELPER + DAY-LENGTH ASSIGNMENT
# ===========================================================================

# Dawn-to-dusk hours at the project centroid, clamped to the configured cap.
# Defaults match the historical Westport values so a caller that passes no
# location still reproduces prior behavior.
bss_day_length_civil <- function(dates, params = list()) {
  sun <- suncalc::getSunlightTimes(
    date = dates,
    lat  = params$centroid_lat %||% 46.904,
    lon  = params$centroid_lon %||% -124.105,
    tz   = params$local_tz     %||% "America/Los_Angeles"
  )
  dl <- as.numeric(difftime(sun$dusk, sun$dawn, units = "hours"))
  lo <- params$day_length_min_hours %||% 9.0
  hi <- params$day_length_max_hours %||% 17.0
  pmax(pmin(dl, hi), lo)
}

# Attach day_length, day_length_civil_twilight, L_mu and L_prior_sigma to a days
# tibble. `days` must already carry event_date and day_type.
#
#   day_length                shore day length used by the PE (L_effective when
#                             available, else civil twilight)
#   day_length_civil_twilight civil twilight, always retained for diagnostics
#   L_mu                      shore day-length point value handed to Stan
#   L_prior_sigma             log-scale uncertainty on L_mu (pooled passes this to
#                             Stan as a prior when estimate_L = 1; the
#                             gear-resolved model currently takes L_mu as data)
#
# L_eff_model = NULL means no usable I/E data: fall back to civil twilight and
# warn, because that path reintroduces the ~2x shore effort bias.
bss_assign_day_length <- function(days, L_eff_model, params = list()) {

  days$day_length <- bss_day_length_civil(days$event_date, params)
  days$day_length_civil_twilight <- days$day_length

  # Default rung: civil twilight with moderate (~30% log-scale) uncertainty.
  days$L_mu <- days$day_length
  days$L_prior_sigma <- rep(0.3, nrow(days))

  if (!is.null(L_eff_model)) {
    L_pred <- L_eff_model$predict_fn(days$event_date, days$day_type)
    days$L_mu <- L_pred$L_mu
    days$L_prior_sigma <- L_pred$L_sigma
    days$day_length <- L_pred$L_mu   # PE uses the predicted median
    attr(days, "l_source") <- L_eff_model$method %||% "L_effective"
  } else {
    attr(days, "l_source") <- "civil_twilight"
    cat("  WARNING: no L_effective model; shore day length falls back to civil",
        "twilight.\n           Shore effort is expected to be overestimated by",
        "roughly 2x on this path.\n")
  }

  days
}


# ===========================================================================
# SHORE TURNOVER FROM THE I/E TIME COLUMN  (review item 2, 2026-09-08)
#
# THE PROBLEM. tau_shore_prior_mu = 1.7 is arrivals / PEAK presence over the WDF20 I/E
# days. The Stan effort likelihood, Gear_I ~ NB2(lambda_E * R_G), calibrates lambda_E to
# the gear counts AT THE TIMES THEY WERE TAKEN, and the expansion E = lambda_E * R_G * tau
# is unbiased only if those counts sit at the daily peak. With the I/E `time` column the
# presence curves show that they do not: 81% of 2024-25 Float 20 counts fall between 10:00
# and 13:59, and presence in that window averages 0.69 of the daily peak (40 days). The
# multiplier the counts actually need is arrivals / presence-at-count-time, about 2.5.
#
# THE ESTIMATOR. For each I/E day d and each clock hour h, presence_d(h) / arrivals_d is the
# fraction of the day's crabber trips present at h (the diel profile). The season's count-
# time distribution w(h) comes from the shore effort counts themselves (count_hour, kept
# by fetch_crab_data since 2026-09-08; sequences <= bss_max_count_seq). Then
#     tau = sum_d arrivals_d / sum_d sum_h w(h) presence_d(h)         (ratio of sums)
# is the constant that makes count x tau unbiased for daily trips over the season, and it
# is the quantity the model's lambda_E * tau needs. Reported beside it: the geometric mean
# of the per-day ratios, the between-day log-SD (the natural day-to-day spread for the
# shared-turnover model), a day-resampling bootstrap log-SE for the level, the same by day
# type, and the profile itself. All go to shore_turnover_*.csv every run.
#
# CROSS-CHECK that justifies applying the profile to gear counts: on the six 2024-25 I/E
# days that also carry a Float 20 gear count, count / R_G matches I/E presence at the same
# instant to 3% (ratio 1.03), so the gear count is a presence snapshot.
#
# It is DERIVED every run and adopted only when params$tau_shore_prior_mu = "derived"
# (see bss_resolve_tau_shore_prior()), because it moves the shore component by the ratio of
# the two priors, about 1.47 on 2024-25, and must be validated by run first.
# ===========================================================================

.ie_hour_of <- function(x) {
  # readxl returns a time column as POSIXct (1899-12-31 base) or hms/difftime; a text
  # export gives "HH:MM:SS". Return the decimal hour, NA where unparseable.
  if (inherits(x, "POSIXct")) return(as.numeric(format(x, "%H")) + as.numeric(format(x, "%M")) / 60)
  if (inherits(x, "difftime") || inherits(x, "hms")) return(as.numeric(x, units = "hours") %% 24)
  xs <- as.character(x)
  hh <- suppressWarnings(as.numeric(sub("^\\s*(\\d{1,2}):(\\d{2}).*$", "\\1", xs)))
  mm <- suppressWarnings(as.numeric(sub("^\\s*(\\d{1,2}):(\\d{2}).*$", "\\2", xs)))
  ifelse(is.finite(hh) & is.finite(mm), hh + mm / 60, NA_real_)
}

estimate_shore_turnover <- function(ie_intervals, shore_effort, params = list(),
                                    n_boot = 2000, seed = 1L, quiet = FALSE) {
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  empty <- list(tau = NA_real_, tau_geomean = NA_real_, log_se = NA_real_, log_sd_days = NA_real_,
                tau_peak = NA_real_, n_days = 0L, method = "unavailable",
                profile = NULL, by_day = NULL, count_time_weights = NULL, by_day_type = NULL)
  if (is.null(ie_intervals) || !is.data.frame(ie_intervals) || !nrow(ie_intervals)) {
    .say("  Shore turnover: no I/E interval rows with a time column; derivation unavailable (peak-based prior stays).\n")
    return(empty)
  }
  # WHICH DAYS THE TURNOVER IS DERIVED FROM (2026-09-11). Until the first full ladder run
  # this function silently used EVERY I/E interval row in the workbook, and nothing said
  # so. On the 2024-25 run that is 40 days spanning 2023-08 to 2026-08, of which SEVEN are
  # inside the 2024-25 season and FOUR inside the all-gear sub-season the resulting 2.477
  # is applied to. It therefore returns the SAME number for every season, while
  # NEW_SEASON_GUIDE.md described it as derived "from the window's I/E time column and
  # count hours" -- a claim that was true of the boat side (`tau_boat_prior_mu =
  # "calibration"`, which reads the window-filtered overlap) and never true of this one.
  # The count-time WEIGHTS below are in-window (shore_effort is), so what shipped was a
  # hybrid: in-window count hours weighting an all-seasons diel profile.
  #
  # Pooling is defensible -- seven days is a thin basis for a quantity that scales the
  # whole shore component by 1.36 -- so it stays the DEFAULT. But it has to be a stated
  # choice with a number attached, not an accident, because this is the second-largest
  # mover in the whole improvement series. Both are now reported every run
  # (n_days_in_window, tau_in_window in shore_turnover_summary.csv) and
  # tau_shore_derive_window_only = TRUE restricts the derivation to the estimation window
  # so the alternative can be priced without editing code. See CHANGE_REGISTER D24.
  .ws <- suppressWarnings(as.Date(params$est_date_start %||% NA))
  .we <- suppressWarnings(as.Date(params$est_date_end   %||% NA))
  iv_all <- ie_intervals |>
    mutate(hbin = floor(hour)) |>
    group_by(event_date) |>
    mutate(arrivals = sum(crabbers_on), peak = max(crabber_flow)) |>
    ungroup() |>
    filter(arrivals > 0, peak > 0)
  .in_win <- if (is.na(.ws) || is.na(.we)) rep(TRUE, nrow(iv_all))
             else as.Date(iv_all$event_date) >= .ws & as.Date(iv_all$event_date) <= .we
  .n_win <- length(unique(iv_all$event_date[.in_win]))
  .n_all <- length(unique(iv_all$event_date))
  iv <- if (isTRUE(params$tau_shore_derive_window_only)) iv_all[.in_win, , drop = FALSE] else iv_all
  if (!nrow(iv)) { .say("  Shore turnover: no I/E day with arrivals; derivation unavailable.\n"); return(empty) }
  .say(sprintf(paste0("  Shore turnover derived from %d I/E day(s)%s; %d of the %d days in the workbook fall INSIDE the",
                      " estimation window (%s to %s).\n"),
               length(unique(iv$event_date)),
               if (isTRUE(params$tau_shore_derive_window_only)) " (window-only: tau_shore_derive_window_only = TRUE)" else " (ALL seasons pooled; the default)",
               .n_win, .n_all, if (is.na(.ws)) "?" else as.character(.ws), if (is.na(.we)) "?" else as.character(.we)))
  if (!isTRUE(params$tau_shore_derive_window_only) && .n_all > 0 && .n_win / .n_all < 0.5)
    .say(sprintf(paste0("  *** NOTE: only %.0f%% of the I/E days behind this turnover are inside the window it is",
                        " applied to. The shore component scales linearly in it; read sigma_IE in the fitted",
                        " output, which is the in-window check. ***\n"), 100 * .n_win / max(.n_all, 1)))

  # --- count-time weights: the hours the creel counts were actually taken ------------
  max_seq <- params$bss_max_count_seq %||% 3
  ch <- NULL
  if (!is.null(shore_effort) && "count_hour" %in% names(shore_effort))
    ch <- shore_effort |> filter(count_sequence <= max_seq, is.finite(count_hour)) |> pull(count_hour)
  if (is.null(ch) || !length(ch)) {
    .say("  Shore turnover: no count hours in the shore effort table; weighting the 10:00-14:00 window uniformly.\n")
    ch <- c(10.5, 11.5, 12.5, 13.5)
  }
  w <- tibble(hbin = floor(ch)) |> count(hbin, name = "n") |> mutate(weight = n / sum(n))

  # --- diel profile: presence / arrivals by hour bin, averaged over days ---------------
  # A day contributes a bin only if the survey covered it; bins outside every survey are
  # absent from the profile and get zero weight below (with a printed note).
  day_hour <- iv |>
    group_by(event_date, season, day_type, hbin, arrivals, peak) |>
    summarise(presence = mean(crabber_flow), .groups = "drop")
  profile <- day_hour |>
    mutate(pres_over_arr = presence / arrivals) |>
    group_by(hbin) |>
    summarise(n_days = n(), mean_presence_over_arrivals = mean(pres_over_arr),
              median_presence_over_arrivals = median(pres_over_arr), .groups = "drop") |>
    mutate(implied_turnover = 1 / mean_presence_over_arrivals)

  # --- per-day expected presence at the count times, then the ratio of sums ------------
  w_use <- w |> filter(hbin %in% day_hour$hbin)
  if (!nrow(w_use)) { .say("  Shore turnover: count hours fall outside every I/E survey; derivation unavailable.\n"); return(empty) }
  if (nrow(w_use) < nrow(w))
    .say(sprintf("  Shore turnover: %.0f%% of count hours fall outside the I/E survey window and are dropped from the weighting.\n",
                 100 * (1 - sum(w_use$n) / sum(w$n))))
  w_use <- w_use |> mutate(weight = n / sum(n))
  by_day <- day_hour |>
    inner_join(w_use |> select(hbin, weight), by = "hbin") |>
    group_by(event_date, season, day_type, arrivals, peak) |>
    # a day that lacks some weighted bins is renormalized over the bins it has
    summarise(presence_at_counts = sum(presence * weight) / sum(weight), .groups = "drop") |>
    mutate(tau_day = arrivals / presence_at_counts, tau_peak_day = arrivals / peak) |>
    filter(is.finite(tau_day), presence_at_counts > 0)

  tau_ros  <- sum(by_day$arrivals) / sum(by_day$presence_at_counts)
  tau_geo  <- exp(mean(log(by_day$tau_day)))
  lsd      <- stats::sd(log(by_day$tau_day))
  tau_peak <- exp(mean(log(by_day$tau_peak_day)))
  set.seed(seed)
  bs <- replicate(n_boot, { i <- sample.int(nrow(by_day), replace = TRUE)
                            sum(by_day$arrivals[i]) / sum(by_day$presence_at_counts[i]) })
  log_se <- stats::sd(log(bs))
  by_dt <- by_day |>
    mutate(dt = ifelse(tolower(day_type) %in% c("weekend", "holiday"), "weekend", "weekday")) |>
    group_by(dt) |>
    summarise(n_days = n(), tau_ratio_of_sums = sum(arrivals) / sum(presence_at_counts),
              tau_geomean = exp(mean(log(tau_day))), .groups = "drop")

  .say(sprintf(paste0("  Shore turnover from the I/E time column: %d days; arrivals / presence at the season's ",
                      "count hours = %.3f (ratio of sums; geometric mean of daily ratios %.3f; between-day ",
                      "log-SD %.2f; bootstrap log-SE %.3f). Peak-based value (the pre-2026-09-08 prior): %.3f.\n"),
               nrow(by_day), tau_ros, tau_geo, lsd, log_se, tau_peak))
  .say(sprintf("    presence at the count hours averages %.2f of the daily peak; weekday %.2f / weekend %.2f\n",
               mean(by_day$presence_at_counts / by_day$peak),
               by_dt$tau_ratio_of_sums[by_dt$dt == "weekday"] %||% NA_real_,
               by_dt$tau_ratio_of_sums[by_dt$dt == "weekend"] %||% NA_real_))

  # the IN-WINDOW subset, always computed and always reported, whatever the derivation used
  .bw <- if (is.na(.ws) || is.na(.we)) by_day
         else by_day[as.Date(by_day$event_date) >= .ws & as.Date(by_day$event_date) <= .we, , drop = FALSE]
  tau_win <- if (nrow(.bw) && sum(.bw$presence_at_counts) > 0)
    sum(.bw$arrivals) / sum(.bw$presence_at_counts) else NA_real_
  if (nrow(.bw) < nrow(by_day))
    .say(sprintf(paste0("    IN-WINDOW SUBSET: %d of %d day(s), ratio of sums %s against %.3f pooled.",
                        " The pooled value is what ships; tau_shore_derive_window_only = TRUE uses the subset.\n"),
                 nrow(.bw), nrow(by_day), if (is.na(tau_win)) "unavailable" else sprintf("%.3f", tau_win), tau_ros))

  list(tau = tau_ros, tau_geomean = tau_geo, log_se = log_se, log_sd_days = lsd,
       tau_peak = tau_peak, n_days = nrow(by_day), method = "count-time-weighted ratio of sums",
       n_days_in_window = nrow(.bw), n_days_total = .n_all, tau_in_window = tau_win,
       derived_window_only = isTRUE(params$tau_shore_derive_window_only),
       profile = profile, by_day = by_day, count_time_weights = w_use, by_day_type = by_dt)
}

# Write the derivation's tables into a run folder (guarded like every other writer).
write_shore_turnover <- function(st, output_dir) {
  if (is.null(st) || is.null(output_dir) || is.null(st$by_day)) return(invisible(NULL))
  tryCatch({
    utils::write.csv(st$profile, file.path(output_dir, "shore_turnover_profile.csv"), row.names = FALSE)
    utils::write.csv(dplyr::mutate(st$by_day, event_date = as.character(event_date)),
                     file.path(output_dir, "shore_turnover_by_day.csv"), row.names = FALSE)
    utils::write.csv(tibble::tibble(
      method = st$method, n_days = st$n_days, tau_ratio_of_sums = st$tau, tau_geomean = st$tau_geomean,
      log_se_bootstrap = st$log_se, log_sd_between_days = st$log_sd_days, tau_peak_based = st$tau_peak,
      weekday = st$by_day_type$tau_ratio_of_sums[st$by_day_type$dt == "weekday"] %||% NA_real_,
      weekend = st$by_day_type$tau_ratio_of_sums[st$by_day_type$dt == "weekend"] %||% NA_real_,
      # 2026-09-11 (D24): the shipped derivation POOLS every I/E day in the workbook, so
      # these three columns are how a reader sees what the number rests on. On the 2024-25
      # run n_days_in_window was 7 of 40 and the shore component scales linearly in tau.
      n_days_in_window = st$n_days_in_window %||% NA_integer_,
      n_days_total = st$n_days_total %||% NA_integer_,
      tau_in_window = st$tau_in_window %||% NA_real_,
      derived_window_only = isTRUE(st$derived_window_only)),
      file.path(output_dir, "shore_turnover_summary.csv"), row.names = FALSE)
  }, error = function(e) cat("  (shore turnover CSVs not written:", conditionMessage(e), ")\n"))
  invisible(NULL)
}

# Resolve tau_shore_prior_mu / _sigma = "derived" from the estimate above. The resolved
# numbers REPLACE the keys (as bss_resolve_tau_boat_prior does), so every consumer (the
# effort spec, run_pe_pooled / run_pe_gear, pe_monthly_effort_share, the gear driver's
# daily-combined series) reads one value. When the level is derived, the shared-turnover
# informed-day floor is waived for shore (params$shared_tau_min_obs$shore <- 0): the floor
# existed to stop four in-window I/E days dragging a level that rested on a 0.3 log-SD
# prior; a level anchored on 40 days with log-SE ~0.06 is not draggable, and the shared
# level is what carries the level uncertainty into the shore interval (289 independent
# per-day draws average it away).
bss_resolve_tau_shore_prior <- function(params, shore_turnover = NULL, quiet = FALSE) {
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  mu_raw <- params$tau_shore_prior_mu    %||% 1.7
  sd_raw <- params$tau_shore_prior_sigma %||% 0.3
  floor_sd <- as.numeric(params$tau_shore_prior_sigma_floor %||% 0.10)
  derived <- function(x) is.character(x) && identical(tolower(x), "derived")
  if (is.numeric(mu_raw) && is.numeric(sd_raw)) {
    params$tau_shore_prior_source <- "config (numeric)"
    return(params)
  }
  st <- shore_turnover
  ok <- !is.null(st) && is.finite(st$tau %||% NA_real_) && (st$n_days %||% 0) >= as.integer(params$tau_shore_derive_min_days %||% 10L)
  if (derived(mu_raw)) {
    if (!ok) {
      params$tau_shore_prior_mu <- as.numeric(params$tau_shore_prior_mu_fallback %||% 1.7)
      params$tau_shore_prior_source <- "fallback (turnover derivation unavailable)"
      .say(sprintf("  tau_shore prior: derivation unavailable; using the fallback centre %.2f.\n", params$tau_shore_prior_mu))
    } else {
      params$tau_shore_prior_mu <- st$tau
      params$tau_shore_prior_source <- sprintf("derived (I/E time column, %d days, %s)", st$n_days, st$method)
      if (is.list(params$shared_tau_min_obs) || is.numeric(params$shared_tau_min_obs)) {
        smo <- params$shared_tau_min_obs
        if (!is.list(smo)) smo <- list(shore = smo, private_boat = smo)
        smo$shore <- 0L
        params$shared_tau_min_obs <- smo
      } else params$shared_tau_min_obs <- list(shore = 0L, private_boat = 15L)
      .say(sprintf("  tau_shore prior centre RESOLVED from the I/E time column: %.3f (%d days); shore shared-turnover floor waived.\n",
                   st$tau, st$n_days))
    }
  } else if (!is.numeric(mu_raw)) {
    stop("params$tau_shore_prior_mu must be a positive number or \"derived\" (got ", deparse(mu_raw), ").", call. = FALSE)
  }
  if (derived(sd_raw)) {
    params$tau_shore_prior_sigma <- if (ok) max(st$log_se, floor_sd) else 0.3
    .say(sprintf("  tau_shore prior log-SD RESOLVED: %.3f (bootstrap log-SE %s, floor %.2f).\n",
                 params$tau_shore_prior_sigma, if (ok) sprintf("%.3f", st$log_se) else "n/a", floor_sd))
  } else if (!is.numeric(sd_raw)) {
    stop("params$tau_shore_prior_sigma must be a positive number or \"derived\".", call. = FALSE)
  }
  params
}
