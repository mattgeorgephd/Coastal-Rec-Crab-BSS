###############################################################################
# bss_cpue_diagnostics.R
#
# F4: the diagnostics that would have caught the 2026-07-09 boat inflation.
# Shared by both drivers (auto-sourced from 03_R_functions/).
#
# THREE INDEPENDENT CHECKS, each of which flagged the same defect:
#
# 1. CPUE ESTIMATOR TRIAD  (bss_cpue_estimator_triad)
#    The likelihood c[a] ~ NB2(lambda_C * h[a], r_C) has a dispersion r_C that
#    does not scale with h. As r_C -> Inf the MLE for lambda_C approaches the
#    RATIO-OF-SUMS (sum catch / sum h); as r_C -> 0 the error becomes purely
#    multiplicative and it approaches the MEAN-OF-RATIOS (mean of catch_a / h_a).
#    With r_C ~ 0.8 and an h that spans three orders of magnitude, the fitted
#    lambda_C drifts far above ratio-of-sums. Reporting all three side by side
#    makes that drift visible:
#        boat 2024-25: ratio-of-sums 0.309, model 0.547, mean-of-ratios 1.704.
#    A model-implied value near mean-of-ratios is a warning, not a result.
#
# 2. SATURATION EXPONENT  (bss_saturation_exponent)
#    Fits log(catch per gear) ~ beta * log(hours per gear). The likelihood
#    assumes catch is proportional to the denominator, i.e. beta = 1. The boat
#    data give beta = 0.133: pots saturate and soak time barely matters. Any
#    beta well below 1 means the chosen effort unit is not a valid one.
#
# 3. UNIT ASSERTION  (bss_assert_effort_units)
#    C = E * lambda_C is only catch if E and h carry the same unit. This fails
#    loudly rather than silently producing a number. It is what would have caught
#    the pooled model, where E is group-hours and h is gear-hours.
#
# Reference: Betancourt (2017) on divergences as geometry; the estimator point is
# standard NB/quasi-likelihood weighting, see McCullagh & Nelder (1989) ch. 9.
###############################################################################


# Ratio-of-sums, mean-of-ratios, and the model-implied CPUE, side by side.
# `cpue_data` needs columns `catch` and `h`.
bss_cpue_estimator_triad <- function(cpue_data, model_lambda_C = NA_real_, label = "") {
  d <- cpue_data[is.finite(cpue_data$h) & cpue_data$h > 0 &
                 is.finite(cpue_data$catch), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)

  ratio_of_sums <- sum(d$catch) / sum(d$h)
  mean_of_ratios <- mean(d$catch / d$h)

  # Where does the model sit between the two? 0 = ratio-of-sums, 1 = mean-of-ratios.
  pos <- if (is.finite(model_lambda_C) && mean_of_ratios > ratio_of_sums) {
    (model_lambda_C - ratio_of_sums) / (mean_of_ratios - ratio_of_sums)
  } else NA_real_

  flag <- is.finite(pos) && pos > 0.25
  if (flag) {
    cat(sprintf(paste0("  WARNING: %s - model CPUE %.4f sits %.0f%% of the way from\n",
                       "           ratio-of-sums (%.4f) toward mean-of-ratios (%.4f).\n",
                       "           The NB dispersion is pulling lambda_C off the rate scale;\n",
                       "           check the saturation exponent before trusting the total.\n"),
                label, model_lambda_C, pos * 100, ratio_of_sums, mean_of_ratios))
  }

  tibble(
    fit = label,
    n_obs = nrow(d),
    h_min = min(d$h), h_median = median(d$h), h_max = max(d$h),
    cpue_ratio_of_sums = ratio_of_sums,
    cpue_model_implied = model_lambda_C,
    cpue_mean_of_ratios = mean_of_ratios,
    model_position_0rs_1mor = pos,
    estimator_drift_flag = flag
  )
}


# Saturation exponent: catch_per_gear ~ c * (hours_per_gear)^beta.
# beta = 1 is the linear-in-effort assumption the likelihood makes. Returns NULL
# when the fit lacks the columns (e.g. shore, where number_of_gear is absent).
bss_saturation_exponent <- function(cpue_data, label = "") {
  need <- c("number_of_gear", "gear_time_total", "catch")
  if (!all(need %in% names(cpue_data))) return(NULL)

  d <- cpue_data[is.finite(cpue_data$number_of_gear) & cpue_data$number_of_gear > 0 &
                 is.finite(cpue_data$gear_time_total) & cpue_data$gear_time_total > 0 &
                 is.finite(cpue_data$catch) & cpue_data$catch > 0, , drop = FALSE]
  if (nrow(d) < 30) return(NULL)

  hours_per_gear <- d$gear_time_total / d$number_of_gear
  catch_per_gear <- d$catch / d$number_of_gear
  ok <- is.finite(hours_per_gear) & hours_per_gear > 0 & is.finite(catch_per_gear) & catch_per_gear > 0
  if (sum(ok) < 30) return(NULL)

  fit <- stats::lm(log(catch_per_gear[ok]) ~ log(hours_per_gear[ok]))
  beta <- unname(stats::coef(fit)[2])
  se   <- unname(summary(fit)$coefficients[2, 2])

  # beta well below 1 means catch does not scale with the denominator.
  flag <- is.finite(beta) && (beta + 2 * se) < 0.75
  if (flag) {
    cat(sprintf(paste0("  WARNING: %s - saturation exponent beta = %.3f (SE %.3f).\n",
                       "           The likelihood assumes beta = 1 (catch proportional to h).\n",
                       "           Effort-hours is NOT a valid unit here; use deployments.\n"),
                label, beta, se))
  }

  tibble(
    fit = label, n_obs = sum(ok),
    beta = beta, beta_se = se,
    beta_lo95 = beta - 1.96 * se, beta_hi95 = beta + 1.96 * se,
    assumed_beta = 1.0,
    hours_per_gear_median = median(hours_per_gear[ok]),
    saturation_flag = flag
  )
}


# Generic test of the likelihood's core assumption, applicable to BOTH populations.
# c[a] ~ NB2(lambda_C * h[a], r_C) asserts log E[c] = log(lambda_C) + 1 * log(h),
# i.e. catch is PROPORTIONAL to the chosen denominator. Fit
#     glm(catch ~ log(h), family = quasipoisson)
# and read the coefficient on log(h): it should be 1. Anything well below 1 means
# catch saturates in h and h is not a valid effort unit. This is the check that
# generalizes the boat's soak-time finding to shore's crabber-hours, where the
# gear-specific test below does not apply.
bss_effort_linearity <- function(cpue_data, label = "") {
  d <- cpue_data[is.finite(cpue_data$h) & cpue_data$h > 0 &
                 is.finite(cpue_data$catch) & cpue_data$catch >= 0, , drop = FALSE]
  if (nrow(d) < 30) return(NULL)
  fit <- try(stats::glm(catch ~ log(h), family = stats::quasipoisson(link = "log"),
                        data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  co <- summary(fit)$coefficients
  if (!"log(h)" %in% rownames(co)) return(NULL)
  b  <- co["log(h)", "Estimate"]; se <- co["log(h)", "Std. Error"]
  flag <- is.finite(b) && is.finite(se) && (b + 1.96 * se) < 0.80

  if (flag) {
    cat(sprintf(paste0("  WARNING: %s - catch scales as h^%.3f (95%% CI %.3f-%.3f), not h^1.\n",
                       "           The likelihood assumes catch is proportional to the CPUE\n",
                       "           denominator. This effort unit is not valid; totals built on\n",
                       "           it will be biased. Check cpue_linearity_%s.csv.\n"),
                label, b, b - 1.96 * se, b + 1.96 * se, label))
  }
  tibble(fit = label, n_obs = nrow(d),
         beta_h = b, beta_h_se = se,
         beta_h_lo95 = b - 1.96 * se, beta_h_hi95 = b + 1.96 * se,
         assumed_beta_h = 1.0, linearity_flag = flag)
}


# C = E * lambda_C is catch only if E and h share a unit. Hard stop otherwise.
bss_assert_effort_units <- function(effort_unit, h_unit, label = "") {
  if (is.null(effort_unit) || is.null(h_unit)) return(invisible(TRUE))
  if (!identical(effort_unit, h_unit)) {
    stop(sprintf(paste0("Unit mismatch in %s: effort E is '%s' but the CPUE ",
                        "denominator h is '%s'. C = E * lambda_C is not catch. ",
                        "Fix the effort scale before reporting."),
                 label, effort_unit, h_unit), call. = FALSE)
  }
  invisible(TRUE)
}


# Writer. Pulls the model-implied CPUE as C_expected_sum / E_sum (both posterior
# medians), which is the quantity that actually multiplies effort into catch.
write_cpue_diagnostics <- function(b, label, output_dir) {
  sd_ <- b$bss_data
  if (is.null(sd_)) return(invisible(NULL))
  cpue_data <- attr(sd_, "cpue_data")
  if (is.null(cpue_data)) cpue_data <- sd_$.cpue_data
  if (is.null(cpue_data) || nrow(cpue_data) == 0) return(invisible(NULL))

  bss_assert_effort_units(sd_$.effort_unit, sd_$.h_unit, label)

  model_lambda_C <- NA_real_
  implied_per_gear_day <- NA_real_
  if (!is.null(b$fit)) {
    s <- tryCatch(rstan::summary(b$fit, pars = c("C_expected_sum", "E_sum"))$summary,
                  error = function(e) NULL)
    if (!is.null(s)) {
      Cx <- s["C_expected_sum", "50%"]; Ex <- s["E_sum", "50%"]
      if (is.finite(Ex) && Ex > 0) model_lambda_C <- Cx / Ex
      D <- sd_$D
      if (is.finite(Cx) && !is.null(D) && D > 0) implied_per_gear_day <- Cx / D
    }
  }

  triad <- bss_cpue_estimator_triad(cpue_data, model_lambda_C, label)
  if (!is.null(triad)) {
    triad$effort_unit <- sd_$.effort_unit %||% NA_character_
    triad$implied_catch_per_day <- implied_per_gear_day
    write.csv(triad, file.path(output_dir, sprintf("cpue_estimators_%s.csv", label)),
              row.names = FALSE)
  }

  sat <- bss_saturation_exponent(cpue_data, label)
  if (!is.null(sat)) {
    write.csv(sat, file.path(output_dir, sprintf("cpue_saturation_%s.csv", label)),
              row.names = FALSE)
  }

  lin <- bss_effort_linearity(cpue_data, label)
  if (!is.null(lin)) {
    lin$effort_unit <- sd_$.effort_unit %||% NA_character_
    write.csv(lin, file.path(output_dir, sprintf("cpue_linearity_%s.csv", label)),
              row.names = FALSE)
  }

  invisible(list(triad = triad, saturation = sat, linearity = lin))
}
