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
#   Boats are unaffected: their gear soaks continuously, so L = 24 by
#   construction (the gear-hours formulation), set in each driver's prep_bss_crab.
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
#   estimate_L_effective(ie_data, pot_open, params)  driver (v7.4), unmodified
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

  ie_raw <- read_excel(ie_file, sheet = params$ie_sheet) |>
    mutate(event_date = as.Date(date))

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
      L_effective = if_else(max_present > 0, ie_crabber_hours / max_present, 0)
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
      L_effective = if_else(max_present > 0, ie_crabber_hours / max_present, 0)
    )

  ie_all <- bind_rows(ie_shore, ie_boat) |>
    filter(ie_crabber_hours > 0)

  cat(sprintf("  I/E survey days: %d shore (WDF20), %d boat (WBL)\n",
              sum(ie_all$population == "shore"),
              sum(ie_all$population == "private_boat")))
  if(nrow(ie_all) > 0)
    cat(sprintf("  Date range: %s to %s\n", min(ie_all$event_date), max(ie_all$event_date)))

  if(nrow(ie_all |> filter(population == "shore")) > 0) {
    shore_ie <- ie_all |> filter(population == "shore")
    cat(sprintf("  Shore L_effective: mean=%.1f hrs (range %.1f–%.1f)\n",
                mean(shore_ie$L_effective), min(shore_ie$L_effective), max(shore_ie$L_effective)))
  }

  return(ie_all)
}


# ===========================================================================
# 1.6 L_EFFECTIVE MODEL - Regression with uncertainty propagation
#
# Fits a regression of log(L_effective) on day-of-year (quadratic) and day type.
# This captures the seasonal gradient WITHIN sub-seasons and provides per-day
# prediction uncertainty for propagation into the Stan model.
#
# Returns: list with
#   $predict_fn: function(event_date, day_type) -> tibble(L_mu, L_sigma)
#   $model: the fitted lm object (NULL on the grand-mean rung)
#   $detail: per-I/E-day data
#   $n_obs, $method: "regression" or "grand_mean"
# ===========================================================================

estimate_L_effective <- function(ie_data, pot_open_date, params) {
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

  cat(sprintf("\n  Predicted L_effective range: %.1f–%.1f hrs\n",
              min(ie_shore$pred_L_mu), max(ie_shore$pred_L_mu)))
  cat(sprintf("  Prediction uncertainty (sigma on log scale): %.2f–%.2f\n",
              min(ie_shore$pred_L_sigma), max(ie_shore$pred_L_sigma)))

  result <- list(
    predict_fn = predict_L,
    model = L_fit,
    detail = ie_shore |> select(event_date, day_type, day_type_group, yday,
                                 ie_crabber_hours, max_present, L_effective,
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
