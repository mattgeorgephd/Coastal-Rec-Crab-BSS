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
# =============================================================================
# save_run_diagnostics.R
#
# Additive run outputs O1-O11 from ADDITIONAL_OUTPUTS_PROPOSAL.md. Persists
# quantities the run computes and holds in memory but did not previously save, so
# the workspace is rarely needed for later troubleshooting. Every writer is
# tryCatch-wrapped and writes only CSVs; nothing here changes an estimate.
#
# Two entry points, mirroring model_diagnostics.R:
#   write_fit_extended_diagnostics(fit, stan_data, days_ss, label, output_dir)
#       per-fit: O1 full parameter summary, O2 AR latent path, O3 AR-period data
#       coverage, O4 modeled daily CPUE, O5 per-observation PPC residuals,
#       O6 HMC sampler diagnostics (incl. E-BFMI), O8 summed-quantity draws,
#       O9 prior-vs-posterior.
#   write_run_level_diagnostics(bss_all, pe_all, gear_props, params, output_dir)
#       run-level: O7 monthly PE-vs-BSS by mode, O10 gear proportions,
#       O11 per-fit data summary.
#
# Sourced by the RMD's purrr::walk(list.files(here("03_R_functions"), ...), source).
# Base R + rstan:: / stats:: / utils:: only, so it does not depend on attached
# packages.
#
# NOTE on O7: the monthly PE share reproduces the v7.0 population-aware Fix-2 logic
# (boat: count-weighted, day-length-free; shore: count * day_length). That logic
# now lives in run_pe(), sections 7.8 and 7.8b, and here; it is a candidate for
# consolidation into one shared helper (review item T4.4).
# =============================================================================

# Null-coalescing fallback, only if not already provided by rlang/purrr.
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ---- shared extraction helpers ---------------------------------------------

# Reduce rstan's array[S] matrix[D,G] (returned as [iter, S, D, G]) to [iter, D]
# for the S=1, G=1 pooled model, guarding dropped size-1 dims. Same approach as
# the hardened get_lam in model_diagnostics.R.
.srd_get_DG <- function(arr, use) {
  d <- dim(arr)
  m <- if (length(d) == 4) arr[use, 1, , 1]
       else if (length(d) == 3) arr[use, , 1]
       else if (length(d) == 2) arr[use, ]
       else stop(sprintf("unexpected array dims: %s", paste(d, collapse = "x")))
  matrix(m, nrow = length(use))
}

# Reduce a [iter, P_n, G*S] matrix (omega) to [iter, P_n] for G*S = 1.
.srd_get_P <- function(arr, use) {
  d <- dim(arr)
  m <- if (length(d) == 3) arr[use, , 1]
       else if (length(d) == 2) arr[use, ]
       else stop(sprintf("unexpected omega dims: %s", paste(d, collapse = "x")))
  matrix(m, nrow = length(use))
}

.srd_q <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE, na.rm = TRUE)

# ---- per-fit writer ---------------------------------------------------------

write_fit_extended_diagnostics <- function(fit, stan_data, days_ss, label, output_dir,
                                           n_pit_draws = 1500, n_draw_save = 2000) {
  if (is.null(fit)) { cat(sprintf("  %s: no fit; skipped.\n", label)); return(invisible(NULL)) }
  ok <- function(tag, expr) tryCatch(expr, error = function(e) {
    cat(sprintf("    [save:%s] %s skipped: %s\n", tag, label, conditionMessage(e))); NULL })

  # day_index -> event_date (day_index = seq_along(event_date) in prep_days_crab)
  ev <- if (!is.null(days_ss) && "event_date" %in% names(days_ss)) as.Date(days_ss$event_date) else NULL
  daylen <- if (!is.null(days_ss) && "day_length" %in% names(days_ss)) as.numeric(days_ss$day_length) else NULL
  D <- stan_data$D

  # ---- O1. Full parameter posterior summary --------------------------------
  ok("O1", {
    s <- rstan::summary(fit)$summary
    df <- data.frame(parameter = rownames(s), s, row.names = NULL, check.names = FALSE)
    utils::write.csv(df, file.path(output_dir, sprintf("bss_full_summary_%s.csv", label)),
                     row.names = FALSE)
  })

  # ---- O6. HMC sampler diagnostics incl. E-BFMI ----------------------------
  ok("O6", {
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    ebfmi <- function(energy) {
      energy <- energy[is.finite(energy)]
      den <- sum((energy - mean(energy))^2)
      if (den <= 0) return(NA_real_)
      sum(diff(energy)^2) / den                       # Betancourt 2016 estimator
    }
    rows <- lapply(seq_along(sp), function(ci) {
      x <- sp[[ci]]
      md <- if ("treedepth__" %in% colnames(x)) max(x[, "treedepth__"]) else NA_real_
      data.frame(
        chain = ci,
        n_iter = nrow(x),
        divergent = if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"]) else NA_real_,
        mean_accept_stat = if ("accept_stat__" %in% colnames(x)) mean(x[, "accept_stat__"]) else NA_real_,
        mean_treedepth = if ("treedepth__" %in% colnames(x)) mean(x[, "treedepth__"]) else NA_real_,
        max_treedepth = md,
        mean_n_leapfrog = if ("n_leapfrog__" %in% colnames(x)) mean(x[, "n_leapfrog__"]) else NA_real_,
        mean_stepsize = if ("stepsize__" %in% colnames(x)) mean(x[, "stepsize__"]) else NA_real_,
        ebfmi = if ("energy__" %in% colnames(x)) ebfmi(x[, "energy__"]) else NA_real_
      )
    })
    df <- do.call(rbind, rows)
    df$ebfmi_low_flag <- is.finite(df$ebfmi) & df$ebfmi < 0.3   # < 0.3 flags energy problems
    utils::write.csv(df, file.path(output_dir, sprintf("sampler_diagnostics_%s.csv", label)),
                     row.names = FALSE)
    if (any(df$ebfmi_low_flag)) cat(sprintf("    [save:O6] %s: low E-BFMI on chain(s) %s\n",
                                             label, paste(df$chain[df$ebfmi_low_flag], collapse = ",")))
  })

  # Draws shared by O2/O3/O4/O5/O8 ------------------------------------------
  # Model-agnostic trailer expansion (pooled R_T vs gear-resolved R_G_boat).
  trailer_par <- bss_trailer_par(fit)
  ex <- ok("extract", rstan::extract(fit, pars = bss_extract_pars(fit,
                                     c("omega_E", "omega_C", "lambda_C_S",
                                       "lambda_E_S", "r_E", "r_C", "R_G",
                                       "C_sum", "C_expected_sum", "E_sum"))))
  if (is.null(ex)) return(invisible(NULL))
  ndraw <- length(ex$r_E)
  use_full <- seq_len(ndraw)
  use_pit  <- if (ndraw > n_pit_draws) sort(sample.int(ndraw, n_pit_draws)) else use_full

  # ---- O8. Summed-quantity posterior draws ---------------------------------
  ok("O8", {
    keep <- if (ndraw > n_draw_save) sort(sample.int(ndraw, n_draw_save)) else use_full
    df <- data.frame(draw = keep,
                     C_sum = as.numeric(ex$C_sum)[keep],
                     C_expected_sum = as.numeric(ex$C_expected_sum)[keep],
                     E_sum = as.numeric(ex$E_sum)[keep])
    utils::write.csv(df, file.path(output_dir, sprintf("bss_draws_summed_%s.csv", label)),
                     row.names = FALSE)
  })

  # ---- O2. AR latent path (per period) -------------------------------------
  period_of_day <- stan_data$period                     # length D, day_index -> period
  P_n <- stan_data$P_n
  ok("O2", {
    oE <- .srd_get_P(ex$omega_E, use_full)               # [iter, P_n]
    oC <- .srd_get_P(ex$omega_C, use_full)
    # date range each period spans
    per_start <- per_end <- rep(NA, P_n); per_ndays <- integer(P_n)
    if (!is.null(ev)) for (p in seq_len(P_n)) {
      dd <- which(period_of_day == p)
      per_ndays[p] <- length(dd)
      if (length(dd) > 0) { per_start[p] <- as.character(min(ev[dd])); per_end[p] <- as.character(max(ev[dd])) }
    }
    qE <- t(apply(oE, 2, .srd_q)); qC <- t(apply(oC, 2, .srd_q))
    df <- data.frame(period = seq_len(P_n), date_start = per_start, date_end = per_end,
                     n_days = per_ndays,
                     omega_E_median = qE[, 2], omega_E_lo95 = qE[, 1], omega_E_hi95 = qE[, 3],
                     omega_C_median = qC[, 2], omega_C_lo95 = qC[, 1], omega_C_hi95 = qC[, 3])
    utils::write.csv(df, file.path(output_dir, sprintf("bss_ar_path_%s.csv", label)),
                     row.names = FALSE)
  })

  # ---- O3. AR-period data coverage / extrapolation map ---------------------
  ok("O3", {
    eff_days <- if (!is.null(stan_data$Gear_n) && stan_data$Gear_n > 0) stan_data$day_Gear
                else if (!is.null(stan_data$T_n) && stan_data$T_n > 0) stan_data$day_T
                else integer(0)
    int_days <- if (!is.null(stan_data$IntC) && stan_data$IntC > 0) stan_data$day_IntC else integer(0)
    ie_days  <- if (!is.null(stan_data$IE_n) && stan_data$IE_n > 0) stan_data$day_IE else integer(0)
    per_eff <- tabulate(period_of_day[eff_days], nbins = P_n)
    per_int <- tabulate(period_of_day[int_days], nbins = P_n)
    per_ie  <- tabulate(period_of_day[ie_days],  nbins = P_n)
    # P4: day-level coverage, resolution-independent. n_effort_obs counts
    #     observations (can exceed n_days); these count distinct DAYS within each
    #     period carrying an observation, so a weekly period reads e.g. "2 of 7
    #     days sampled" and is comparable to a daily fit (n_days = 1). This is what
    #     makes the boat (weekly) and shore (daily) coverage directly comparable.
    ueff <- unique(eff_days); uint <- unique(int_days)
    per_days_eff <- integer(P_n); per_days_int <- integer(P_n)
    for (p in seq_len(P_n)) {
      dd <- which(period_of_day == p)
      per_days_eff[p] <- sum(dd %in% ueff)
      per_days_int[p] <- sum(dd %in% uint)
    }
    oE <- .srd_get_P(ex$omega_E, use_full); oC <- .srd_get_P(ex$omega_C, use_full)
    ciwE <- apply(oE, 2, function(z) diff(stats::quantile(z, c(0.025, 0.975), names = FALSE)))
    ciwC <- apply(oC, 2, function(z) diff(stats::quantile(z, c(0.025, 0.975), names = FALSE)))
    per_start <- per_end <- rep(NA, P_n); per_ndays <- integer(P_n)
    for (p in seq_len(P_n)) per_ndays[p] <- sum(period_of_day == p)
    if (!is.null(ev)) for (p in seq_len(P_n)) {
      dd <- which(period_of_day == p)
      if (length(dd) > 0) { per_start[p] <- as.character(min(ev[dd])); per_end[p] <- as.character(max(ev[dd])) }
    }
    df <- data.frame(period = seq_len(P_n), date_start = per_start, date_end = per_end,
                     n_days = per_ndays,
                     n_days_with_effort = per_days_eff, n_days_with_interview = per_days_int,
                     frac_days_effort = round(per_days_eff / pmax(per_ndays, 1), 3),
                     frac_days_interview = round(per_days_int / pmax(per_ndays, 1), 3),
                     n_effort_obs = per_eff, n_interviews = per_int, n_ie_obs = per_ie,
                     observed = (per_eff > 0 | per_int > 0),
                     omega_E_ci95_width = ciwE, omega_C_ci95_width = ciwC)
    utils::write.csv(df, file.path(output_dir, sprintf("bss_period_coverage_%s.csv", label)),
                     row.names = FALSE)
  })

  # ---- O4. Modeled daily CPUE (lambda_C) vs raw interview CPUE --------------
  ok("O4", {
    lamC <- .srd_get_DG(ex$lambda_C_S, use_full)         # [iter, D]
    qC <- t(apply(lamC, 2, .srd_q))                      # per day
    # raw interview CPUE per day from the fit's own interview data
    raw <- rep(NA_real_, D); nint <- integer(D)
    if (!is.null(stan_data$IntC) && stan_data$IntC > 0) {
      di <- stan_data$day_IntC; cpue_i <- stan_data$c / pmax(stan_data$h, 1e-8)
      ag <- tapply(cpue_i, di, mean); cnt <- tapply(cpue_i, di, length)
      idx <- as.integer(names(ag)); raw[idx] <- as.numeric(ag); nint[idx] <- as.integer(cnt)
    }
    df <- data.frame(day_index = seq_len(D),
                     event_date = if (!is.null(ev)) as.character(ev) else NA,
                     cpue_model_median = qC[, 2], cpue_model_lo95 = qC[, 1], cpue_model_hi95 = qC[, 3],
                     cpue_raw_interview = raw, n_interviews = nint)
    utils::write.csv(df, file.path(output_dir, sprintf("bss_daily_cpue_%s.csv", label)),
                     row.names = FALSE)
  })

  # ---- O5. Per-observation PPC residuals (exact randomized PIT) -------------
  # PIT_i = E_draws[ P(Y < y) + 0.5 P(Y = y) ] = mean over draws of
  #   pnbinom(y-1, size=r, mu) + 0.5 * dnbinom(y, size=r, mu). This is the exact
  #   expectation of the simulated PIT used by the aggregate PPC, with no
  #   simulation noise. in_50 / in_95 are the central-interval coverage flags
  #   (equivalent to PIT in [0.25,0.75] / [0.025,0.975] for the randomized PIT).
  ok("O5", {
    lamE <- .srd_get_DG(ex$lambda_E_S, use_pit)
    lamC <- .srd_get_DG(ex$lambda_C_S, use_pit)
    rE <- as.numeric(ex$r_E[use_pit]); rC <- as.numeric(ex$r_C[use_pit])
    RG <- as.numeric(ex$R_G[use_pit])
    RT <- bss_trailer_multiplier(ex, trailer_par, use_pit)   # R_T or 1/R_G_boat
    pit_block <- function(days, y, mu_mat, size) {
      no <- length(y); pit <- fit_mean <- rep(NA_real_, no)
      for (i in seq_len(no)) {
        mu <- pmax(mu_mat[, i], 1e-8); keep <- is.finite(mu) & is.finite(size)
        if (sum(keep) < 20) next
        mu <- mu[keep]; sz <- size[keep]
        pit[i] <- mean(stats::pnbinom(y[i] - 1, size = sz, mu = mu) +
                       0.5 * stats::dnbinom(y[i], size = sz, mu = mu))
        fit_mean[i] <- mean(mu)
      }
      data.frame(day_index = days,
                 event_date = if (!is.null(ev)) as.character(ev[days]) else NA,
                 observed = y, fitted_mean = round(fit_mean, 3), pit = round(pit, 4),
                 in_50 = pit >= 0.25 & pit <= 0.75,
                 in_95 = pit >= 0.025 & pit <= 0.975)
    }
    parts <- list()
    if (!is.null(stan_data$Gear_n) && stan_data$Gear_n > 0)
      parts$gear <- cbind(data_type = "gear",
                          pit_block(stan_data$day_Gear, stan_data$Gear_I,
                                    lamE[, stan_data$day_Gear, drop = FALSE] * RG, rE))
    if (!is.null(stan_data$T_n) && stan_data$T_n > 0 && !is.null(RT))
      parts$trailer <- cbind(data_type = "trailer",
                             pit_block(stan_data$day_T, stan_data$T_I,
                                       lamE[, stan_data$day_T, drop = FALSE] * RT, rE))
    if (!is.null(stan_data$IntC) && stan_data$IntC > 0) {
      muC <- lamC[, stan_data$day_IntC, drop = FALSE] * rep(stan_data$h, each = nrow(lamC))
      parts$catch <- cbind(data_type = "catch",
                           pit_block(stan_data$day_IntC, stan_data$c, muC, rC))
    }
    if (length(parts) > 0) {
      df <- do.call(rbind, parts)
      utils::write.csv(df, file.path(output_dir, sprintf("ppc_byobs_%s.csv", label)),
                       row.names = FALSE)
    }
  })

  # ---- O9. Prior vs posterior --------------------------------------------
  ok("O9", {
    sd_norm <- function(s) s
    beta_sd <- function(a, b) sqrt(a * b / ((a + b)^2 * (a + b + 1)))
    sd_p <- stan_data
    # Only summarize parameters this model declares. The pooled model carries R_T
    # and B1_C; the gear-resolved model carries R_G_boat and no weekend-CPUE
    # effect. Indexed names (mu_mu_E[1]) are checked on their base name.
    has_par <- function(p) sub("\\[.*$", "", p) %in% fit@model_pars
    b_phiE <- sd_p$value_betashape_phi_E_scaled; b_phiC <- sd_p$value_betashape_phi_C_scaled
    # Pooled passes the R_G prior through stan_data; the gear-resolved Stan model
    # hardcodes lognormal(log(1.3), 0.3), so fall back to those values.
    rg_mu <- sd_p$R_G_prior_mu %||% 1.3; rg_s <- sd_p$R_G_prior_sigma %||% 0.3
    lnorm_mean <- function(mu, s) mu * exp(s^2 / 2)
    lnorm_sd   <- function(mu, s) mu * exp(s^2 / 2) * sqrt(exp(s^2) - 1)
    prior_tbl <- list(
      R_G        = list(fam = sprintf("lognormal(log(%.3f), %.3f)", rg_mu, rg_s),
                        mean = lnorm_mean(rg_mu, rg_s),
                        sd = lnorm_sd(rg_mu, rg_s)),
      phi_E      = list(fam = sprintf("2*beta(%.1f,%.1f)-1", b_phiE, b_phiE), mean = 0,
                        sd = 2 / (2 * sqrt(2 * b_phiE + 1))),
      phi_C      = list(fam = sprintf("2*beta(%.1f,%.1f)-1", b_phiC, b_phiC), mean = 0,
                        sd = 2 / (2 * sqrt(2 * b_phiC + 1))),
      B1         = list(fam = sprintf("normal(0, %.1f)", sd_p$value_normal_sigma_B1), mean = 0, sd = sd_p$value_normal_sigma_B1),
      B2         = list(fam = sprintf("normal(0, %.1f)", sd_p$value_normal_sigma_B2), mean = 0, sd = sd_p$value_normal_sigma_B2),
      `mu_mu_E[1]` = list(fam = sprintf("normal(%.3f, %.1f)", sd_p$value_normal_mu_mu_E, sd_p$value_normal_sigma_mu_E), mean = sd_p$value_normal_mu_mu_E, sd = sd_p$value_normal_sigma_mu_E),
      `mu_mu_C[1]` = list(fam = sprintf("normal(%.3f, %.1f)", sd_p$value_normal_mu_mu_C, sd_p$value_normal_sigma_mu_C), mean = sd_p$value_normal_mu_mu_C, sd = sd_p$value_normal_sigma_mu_C),
      sigma_eps_E = list(fam = sprintf("half-cauchy(0, %.1f)", sd_p$value_cauchyDF_sigma_eps_E), mean = NA, sd = NA),
      sigma_eps_C = list(fam = sprintf("half-cauchy(0, %.1f)", sd_p$value_cauchyDF_sigma_eps_C), mean = NA, sd = NA),
      sigma_r_E  = list(fam = sprintf("half-cauchy(0, %.1f)", sd_p$value_cauchyDF_sigma_r_E), mean = NA, sd = NA),
      sigma_r_C  = list(fam = sprintf("half-cauchy(0, %.1f)", sd_p$value_cauchyDF_sigma_r_C), mean = NA, sd = NA),
      sigma_mu_E = list(fam = sprintf("half-cauchy(0, %.1f)", sd_p$value_cauchyDF_sigma_mu_E), mean = NA, sd = NA),
      sigma_mu_C = list(fam = sprintf("half-cauchy(0, %.1f)", sd_p$value_cauchyDF_sigma_mu_C), mean = NA, sd = NA)
    )

    # --- Model-specific priors, added only when the model declares them -------
    # Pooled: R_T ~ beta(alpha, beta) passed through stan_data; B1_C ~ normal.
    if (!is.null(sd_p$R_T_alpha) && !is.null(sd_p$R_T_beta)) {
      prior_tbl$R_T <- list(fam = sprintf("beta(%.1f, %.1f)", sd_p$R_T_alpha, sd_p$R_T_beta),
                            mean = sd_p$R_T_alpha / (sd_p$R_T_alpha + sd_p$R_T_beta),
                            sd = beta_sd(sd_p$R_T_alpha, sd_p$R_T_beta))
    }
    if (!is.null(sd_p$value_normal_sigma_B1_C)) {
      prior_tbl$B1_C <- list(fam = sprintf("normal(0, %.1f)", sd_p$value_normal_sigma_B1_C),
                             mean = 0, sd = sd_p$value_normal_sigma_B1_C)
    }
    # Gear-resolved: R_G_boat ~ lognormal(log(4), 0.5), hardcoded in the Stan model.
    if (has_par("R_G_boat")) {
      prior_tbl$R_G_boat <- list(fam = "lognormal(log(4.000), 0.500)",
                                 mean = lnorm_mean(4, 0.5), sd = lnorm_sd(4, 0.5))
    }

    pars <- names(prior_tbl)[vapply(names(prior_tbl), has_par, logical(1))]
    if (length(pars) == 0) return(NULL)
    post <- rstan::summary(fit, pars = pars)$summary

    rows <- lapply(pars, function(pn) {
      pr <- prior_tbl[[pn]]; po <- post[pn, ]
      contraction <- if (is.finite(pr$sd) && pr$sd > 0) 1 - po["sd"] / pr$sd else NA_real_
      data.frame(parameter = pn, prior = pr$fam,
                 prior_mean = pr$mean, prior_sd = pr$sd,
                 post_mean = po["mean"], post_sd = po["sd"], post_median = po["50%"],
                 post_lo95 = po["2.5%"], post_hi95 = po["97.5%"],
                 contraction = round(contraction, 3),
                 prior_influential = is.finite(contraction) & contraction < 0.5,
                 row.names = NULL)
    })
    df <- do.call(rbind, rows)
    utils::write.csv(df, file.path(output_dir, sprintf("prior_vs_posterior_%s.csv", label)),
                     row.names = FALSE)
  })

  invisible(TRUE)
}

# ---- O13: PSIS-LOO + Pareto-k influence (per fit) ---------------------------
# Consumes the Stan log_lik_gear / log_lik_trailer / log_lik_catch added to the
# pooled model. Enables the project's primary model-selection tool (PSIS-LOO),
# e.g. comparing AR resolutions for the boat across two runs by elpd_loo, and
# flags influential observations (high Pareto-k), e.g. whether a few sparse-month
# interviews drive a CPUE estimate. Caveat: this is conditional pointwise LOO
# given the latent AR path, not leave-future-out; treat elpd differences as a
# screen and use k-fold time-block CV for the definitive temporal comparison.
# Degrades gracefully: without the 'loo' package it still saves pointwise lpd and
# an across-draw influence proxy.
write_loo_diagnostics <- function(fit, stan_data, days_ss, label, output_dir) {
  if (is.null(fit)) { cat(sprintf("  %s: no fit; LOO skipped.\n", label)); return(invisible(NULL)) }
  ok <- function(tag, expr) tryCatch(expr, error = function(e) {
    cat(sprintf("    [save:%s] %s skipped: %s\n", tag, label, conditionMessage(e))); NULL })
  ev <- if (!is.null(days_ss) && "event_date" %in% names(days_ss)) as.Date(days_ss$event_date) else NULL
  has_loo <- requireNamespace("loo", quietly = TRUE)

  lpd_point <- function(ll) apply(ll, 2, function(col) {   # logsumexp-stable lpd
    m <- max(col); m + log(mean(exp(col - m)))
  })
  streams <- list(
    gear    = list(par = "log_lik_gear",    n = stan_data$Gear_n %||% 0, days = stan_data$day_Gear,  y = stan_data$Gear_I),
    trailer = list(par = "log_lik_trailer", n = stan_data$T_n %||% 0,    days = stan_data$day_T,     y = stan_data$T_I),
    catch   = list(par = "log_lik_catch",   n = stan_data$IntC %||% 0,   days = stan_data$day_IntC,  y = stan_data$c)
  )
  summ <- list()
  for (sn in names(streams)) {
    st <- streams[[sn]]
    if ((st$n %||% 0) == 0) next                            # stream absent in this fit
    res <- ok(paste0("O13:", sn), {
      ll <- as.matrix(rstan::extract(fit, pars = st$par)[[1]])   # [draws x n_obs]
      if (nrow(ll) == 0 || ncol(ll) == 0) return(NULL)
      n_obs <- ncol(ll)
      lpd <- lpd_point(ll)
      ll_sd <- apply(ll, 2, stats::sd)                     # across-draw sd: influence proxy
      dy <- st$days[seq_len(n_obs)]
      df <- data.frame(data_type = sn, obs_index = seq_len(n_obs), day_index = dy,
                       event_date = if (!is.null(ev)) as.character(ev[dy]) else NA,
                       observed = st$y[seq_len(n_obs)],
                       lpd = round(lpd, 4), loglik_sd = round(ll_sd, 4))
      elpd <- p_loo <- se_elpd <- NA_real_; n_khi <- NA_integer_
      if (has_loo) {
        lo <- loo::loo(ll)
        df$pareto_k <- round(lo$diagnostics$pareto_k, 3)
        df$elpd_loo <- round(lo$pointwise[, "elpd_loo"], 4)
        est <- lo$estimates
        elpd <- est["elpd_loo", "Estimate"]; se_elpd <- est["elpd_loo", "SE"]
        p_loo <- est["p_loo", "Estimate"]; n_khi <- sum(lo$diagnostics$pareto_k > 0.7)
      }
      utils::write.csv(df, file.path(output_dir, sprintf("loo_pointwise_%s_%s.csv", sn, label)),
                       row.names = FALSE)
      data.frame(stream = sn, n_obs = n_obs, elpd_loo = elpd, se_elpd_loo = se_elpd,
                 p_loo = p_loo, n_pareto_k_gt_0.7 = n_khi, sum_lpd = round(sum(lpd), 2))
    })
    if (!is.null(res)) summ[[sn]] <- res
  }
  if (length(summ) > 0) ok("O13:summary", {
    utils::write.csv(do.call(rbind, summ),
                     file.path(output_dir, sprintf("loo_summary_%s.csv", label)), row.names = FALSE)
  })
  if (!has_loo) cat(sprintf("    [save:O13] %s: 'loo' not installed; saved pointwise lpd + loglik_sd only (no Pareto-k/elpd). install.packages('loo') for PSIS-LOO.\n", label))
  invisible(TRUE)
}

# ---- run-level writers ------------------------------------------------------

# Population-aware monthly effort share (v7.0 Fix-2 logic; see header note).
# Boat: count-weighted (day-length-free). Shore: count * day_length.
.srd_monthly_share <- function(stan_data, days_ss, is_boat) {
  obs_days <- if (is_boat) stan_data$day_T else stan_data$day_Gear
  counts   <- if (is_boat) stan_data$T_I  else stan_data$Gear_I
  if (length(obs_days) == 0 || is.null(days_ss)) return(NULL)
  mc <- tapply(counts, obs_days, mean)
  di <- as.integer(names(mc))
  ev <- as.Date(days_ss$event_date)[di]
  dl <- if ("day_length" %in% names(days_ss)) as.numeric(days_ss$day_length)[di] else rep(1, length(di))
  w <- if (is_boat) as.numeric(mc) else as.numeric(mc) * dl
  mon <- format(ev, "%Y-%m")
  agg <- tapply(w, mon, sum)
  data.frame(month = names(agg), share = as.numeric(agg) / sum(agg, na.rm = TRUE),
             row.names = NULL)
}

# BSS monthly summed draws for one fit -> month, median, lo95, hi95 (catch & effort).
.srd_bss_monthly <- function(fit, days_ss, use_n = 2000) {
  ex <- rstan::extract(fit, pars = c("C_expected", "E"))
  nd <- dim(ex$C_expected)[1]
  use <- if (nd > use_n) sort(sample.int(nd, use_n)) else seq_len(nd)
  C <- .srd_get_DG(ex$C_expected, use); E <- .srd_get_DG(ex$E, use)   # [draws, D]
  ev <- as.Date(days_ss$event_date); mon <- format(ev, "%Y-%m"); um <- sort(unique(mon))
  out <- lapply(um, function(m) {
    cols <- which(mon == m)
    cs <- rowSums(C[, cols, drop = FALSE]); es <- rowSums(E[, cols, drop = FALSE])
    data.frame(month = m,
               BSS_catch_median = stats::median(cs), BSS_catch_lo95 = stats::quantile(cs, .025, names = FALSE),
               BSS_catch_hi95 = stats::quantile(cs, .975, names = FALSE),
               BSS_effort_median = stats::median(es), BSS_effort_lo95 = stats::quantile(es, .025, names = FALSE),
               BSS_effort_hi95 = stats::quantile(es, .975, names = FALSE))
  })
  do.call(rbind, out)
}

write_run_level_diagnostics <- function(bss_all, pe_all, gear_props, params, output_dir) {
  ok <- function(tag, expr) tryCatch(expr, error = function(e) {
    cat(sprintf("    [save:%s] skipped: %s\n", tag, conditionMessage(e))); NULL })

  # ---- O10. Gear proportions (with interview support) ----------------------
  ok("O10", {
    if (!is.null(gear_props)) {
      gp <- as.data.frame(gear_props)
      utils::write.csv(gp, file.path(output_dir, "gear_proportions.csv"), row.names = FALSE)
    }
  })

  # ---- O11. Per-fit data summary -------------------------------------------
  ok("O11", {
    rows <- lapply(names(bss_all), function(label) {
      b <- bss_all[[label]]
      sd_ <- b$bss_data; ds <- b$days_ss
      pe_label <- paste0(b$population, "_", b$subseason)
      n_eff_obs <- if (is.null(sd_)) NA else (if (!is.null(sd_$Gear_n) && sd_$Gear_n > 0) sd_$Gear_n else sd_$T_n)
      ev <- if (!is.null(ds) && "event_date" %in% names(ds)) as.Date(ds$event_date) else NULL
      pct_days_eff <- NA_real_; pct_days_int <- NA_real_; ints_per_mo <- NA_real_
      if (!is.null(sd_) && !is.null(ev)) {
        D <- sd_$D
        eff_days <- if (!is.null(sd_$Gear_n) && sd_$Gear_n > 0) sd_$day_Gear else sd_$day_T
        pct_days_eff <- round(100 * length(unique(eff_days)) / D, 1)
        pct_days_int <- round(100 * length(unique(sd_$day_IntC)) / D, 1)
        nmo <- length(unique(format(ev, "%Y-%m")))
        ints_per_mo <- if (nmo > 0) round((sd_$IntC %||% 0) / nmo, 1) else NA_real_
      }
      data.frame(
        fit = label, population = b$population %||% NA, subseason = b$subseason %||% NA,
        method = if (isTRUE(b$pe_fallback)) "PE (insufficient data)"
                 else if (isTRUE(b$use_bss)) "BSS" else if (!is.null(b$fit)) "PE (gate)" else "PE",
        ar_resolution = if (!is.null(sd_)) attr(sd_, "ar_resolution") %||% NA else NA,
        P_n = if (!is.null(sd_)) sd_$P_n else NA, D = if (!is.null(sd_)) sd_$D else NA,
        n_effort_obs = n_eff_obs,
        n_interviews = if (!is.null(sd_)) sd_$IntC else NA,
        n_ie_obs = if (!is.null(sd_)) sd_$IE_n else NA,
        date_start = if (!is.null(ev)) as.character(min(ev)) else NA,
        date_end = if (!is.null(ev)) as.character(max(ev)) else NA,
        pct_days_with_effort = pct_days_eff, pct_days_with_interview = pct_days_int,
        mean_interviews_per_month = ints_per_mo,
        pe_catch = pe_all[[pe_label]][[b$catch_group]] %||% NA,
        pe_effort = pe_all[[pe_label]]$effort_total %||% NA,
        row.names = NULL)
    })
    df <- do.call(rbind, rows)
    utils::write.csv(df, file.path(output_dir, "fit_data_summary.csv"), row.names = FALSE)
  })

  # ---- O7. Monthly PE vs BSS by mode ---------------------------------------
  ok("O7", {
    out <- list()
    for (label in names(bss_all)) {
      b <- bss_all[[label]]
      if (is.null(b$population)) next
      pop <- b$population; ds <- b$days_ss; sd_ <- b$bss_data
      pe_label <- paste0(pop, "_", b$subseason)
      pe_catch_tot <- pe_all[[pe_label]][[b$catch_group]] %||% NA
      pe_eff_tot   <- pe_all[[pe_label]]$effort_total %||% NA
      is_boat <- grepl("private_boat", pop)
      share <- if (!is.null(sd_)) .srd_monthly_share(sd_, ds, is_boat) else NULL
      pe_m <- if (!is.null(share)) data.frame(month = share$month,
                                              PE_catch = pe_catch_tot * share$share,
                                              PE_effort = pe_eff_tot * share$share) else NULL
      bss_m <- if (!is.null(b$fit) && isTRUE(b$use_bss)) .srd_bss_monthly(b$fit, ds) else NULL
      months <- sort(unique(c(if (!is.null(pe_m)) pe_m$month, if (!is.null(bss_m)) bss_m$month)))
      if (length(months) == 0) next
      df <- data.frame(mode = pop, subseason = b$subseason, month = months)
      df <- merge(df, pe_m, by = "month", all.x = TRUE)
      if (!is.null(bss_m)) df <- merge(df, bss_m, by = "month", all.x = TRUE)
      out[[label]] <- df
    }
    if (length(out) > 0) {
      allcols <- unique(unlist(lapply(out, names)))
      out <- lapply(out, function(d) { d[setdiff(allcols, names(d))] <- NA; d[allcols] })
      comb <- do.call(rbind, out)
      ordcols <- c("mode", "subseason", "month", "PE_catch", "BSS_catch_median",
                   "BSS_catch_lo95", "BSS_catch_hi95", "PE_effort", "BSS_effort_median",
                   "BSS_effort_lo95", "BSS_effort_hi95")
      ordcols <- intersect(ordcols, names(comb))
      comb <- comb[order(comb$mode, comb$subseason, comb$month), ordcols]
      # BSS columns stay NA (not 0) for sub-seasons with no fit (e.g. boat ring-net).
      utils::write.csv(comb, file.path(output_dir, "monthly_pe_vs_bss.csv"), row.names = FALSE)
    }
  })

  invisible(TRUE)
}
