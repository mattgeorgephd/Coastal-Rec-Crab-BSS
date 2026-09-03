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
# 03_R_functions/model_diagnostics.R
# Per-fit model-behavior diagnostics for the pooled BSS. Auto-sourced by the
# RMD (it sources everything in 03_R_functions/). All functions are read-only on
# the fit and safe to call on any stanfit from this model lineage.
#
# Produced per fit (see write_bss_diagnostics):
#   structural_params_<label>.csv   key scale/structural parameters with CI,
#                                   n_eff, R-hat (the tuning knobs)
#   divergence_localization_<label>.csv  where divergent draws sit vs the bulk,
#                                   plus the divergent-vs-bulk shift in the totals
#   ppc_calibration_<label>.csv     posterior predictive coverage + PIT for
#                                   effort counts and interview catches
#   ppc_pit_<label>.png             PIT histograms (uniform => well calibrated)
# =============================================================================

# --- 1. Structural / scale parameter summary --------------------------------
# Which parameters in a fit's curated table carry NO likelihood contribution (2026-08-27).
#
# Several parameters here are declared with a proper prior UNCONDITIONALLY, so that an unused
# one is decoupled-but-proper rather than an improper flat direction. That is the right
# modelling choice and a serious reporting hazard: a decoupled parameter reports its PRIOR, in
# the same columns and next to the same headings as an estimate, with nothing to separate the
# two. The 2026-08-26 ladder produced the worked example -- under the production
# `osp_scale_is_tau = TRUE` the OSP likelihood uses L instead of kappa_OSP, so kappa_OSP is
# inert and reports its lognormal(log 3, 0.3) prior exactly (median 3.008, 95% 1.63-5.40).
# Read as an estimate, that is a claim the model has MEASURED the boat turnover at 3.0. It has
# not; it has been told 3.0. The shore fits' r_OSP is worse: a 95% interval spanning roughly
# 0.0002 to 700.
#
# Deliberately conservative: a reason is given only where the Stan program's own guard
# demonstrably removes the term, so TRUE always means prior-only. FALSE does not by itself
# certify that a parameter is well identified.
#
# `parameters` may carry Stan indices (f_crab[1]); they are stripped before matching.
# Returns a character vector of reasons, NA where the parameter is coupled.
bss_decoupled_reasons <- function(parameters, stan_data = NULL) {
  reason <- rep(NA_character_, length(parameters))
  if (is.null(stan_data) || !length(parameters)) return(reason)
  g <- function(nm, default = NULL) { v <- stan_data[[nm]]; if (is.null(v)) default else v }
  base <- sub("\\[.*$", "", parameters)
  set  <- function(mask, why) reason[mask & is.na(reason)] <<- why

  ie_n    <- g("IE_n", 0L);   osp_n <- g("OSP_n", 0L);  t_n <- g("T_n", 0L)
  tau_sw  <- g("osp_scale_is_tau", 0L);                 k_open <- g("K_open", 0L)
  apply_f <- g("apply_crab_fraction", 0L); est_f <- g("crab_fraction_estimate", 0L)
  osp_lo  <- g("osp_crab_lower", 0L);      dens  <- g("estimate_cpue_density", 0L)
  w_sum   <- sum(as.numeric(g("w", 0)), na.rm = TRUE)
  h_sum   <- sum(as.numeric(g("holiday", 0)), na.rm = TRUE)

  if (ie_n == 0) set(base == "sigma_IE",
    "no I/E observations in this fit (IE_n = 0); sigma_IE is its prior")
  if (osp_n == 0) set(base %in% c("kappa_OSP", "sigma_r_OSP", "r_OSP"),
    "no OSP observations in this fit (OSP_n = 0)")
  else if (tau_sw == 1) set(base == "kappa_OSP",
    "osp_scale_is_tau = 1: the OSP mean uses L, not kappa_OSP; this is the prior")
  if (t_n == 0 && osp_n == 0) set(base %in% c("R_G_boat", "R_T"),
    "no boat count stream in this fit; the boat gear ratio is its prior")
  if (k_open == 0) set(base == "B_open", "K_open = 0: no opener covariate is active")
  if (dens == 0)   set(base == "gamma_C", "estimate_cpue_density = 0: the density term is inert")
  if (w_sum == 0)  set(base %in% c("B1", "B1_C"), "no weekend day in this fit's window")
  if (h_sum == 0)  set(base %in% c("B2", "B2_C"), "no holiday in this fit's window")
  if (apply_f == 0) set(base %in% c("f_crab", "f_lower"),
    "apply_crab_fraction = 0 (shore, or the feature is off): f is pinned at 1")
  else if (est_f == 0) set(base == "f_crab",
    "crab_fraction_estimate = 0: f is pinned at crab_fraction_value")
  if (osp_lo == 0) set(base == "f_lower",
    "osp_crab_lower = 0: the OSP lower bound is off and f_lower is pinned at 0")
  if (identical(as.integer(g("shared_tau", 0L)), 0L)) set(base %in% c("tau_bar", "tau_bar_out"),
    "shared_tau = 0: L is per-day independent draws and there is no shared turnover")
  # 2026-09-02: theta_C_out is written unconditionally so the reported parameter set does not
  # change shape between runs, which means it reports a hard 0.0 when the feature is off.
  # Without this rule a reader would see "theta_C_out = 0" and take it as an estimate that
  # zero inflation was tested and found absent, which is the opposite of what a zero means
  # there. This is the same reporting hazard kappa_OSP and f_lower carry.
  if (identical(as.integer(g("zi_catch", 0L)), 0L)) set(base %in% c("theta_C", "theta_C_out"),
    "zi_catch = 0: the catch likelihood is plain NB2 and theta_C is not in the model")
  reason
}

bss_structural_summary <- function(fit, stan_data = NULL, fit_method = NULL) {
  pars <- c("mu_mu_E", "mu_mu_C",
            "sigma_eps_E", "sigma_eps_C",
            "phi_E", "phi_C",
            "sigma_r_E", "sigma_r_C", "r_E", "r_C",
            "sigma_mu_E", "sigma_mu_C",
            "sigma_IE", "R_G", "R_T", "R_G_boat",
            "B1", "B2", "B1_C",
            # Added 2026-08-25. Closes the standing Tier-4 item "surface B2_C in the
            # curated report tables": the always-on holiday CPUE term, the opener effort
            # covariates (improvement 4), the density term when active, and the boat
            # scale/fraction parameters were all written only to bss_full_summary_*, so a
            # reader had to open the full summary to see terms that move the estimate.
            "B2_C", "gamma_C", "B_open",
            "kappa_OSP", "sigma_r_OSP", "r_OSP",
            "f_crab", "f_lower",
            # improvement 2.1 (2026-08-27): the shared turnover, when it exists.
            "tau_bar",
            # 2026-09-02: the zero-inflation probability, when the catch likelihood carries one.
            "theta_C")
  pars <- pars[pars %in% fit@model_pars]
  s <- summary(fit, pars = pars)$summary
  out <- data.frame(parameter = rownames(s),
                    mean   = s[, "mean"],
                    lo95   = s[, "2.5%"],
                    median = s[, "50%"],
                    hi95   = s[, "97.5%"],
                    n_eff  = round(s[, "n_eff"]),
                    Rhat   = round(s[, "Rhat"], 4),
                    row.names = NULL)

  # See bss_decoupled_reasons(): a parameter with no likelihood contribution in this fit
  # reports its PRIOR, and must not sit unlabelled beside parameters that were estimated.
  reason <- bss_decoupled_reasons(out$parameter, stan_data)
  out$decoupled        <- !is.na(reason)
  out$decoupled_reason <- reason

  # 2026-08-30: `estimate` is the column to quote, and it is NA for a decoupled parameter.
  #
  # The `decoupled` flag added on 2026-08-27 was necessary and not sufficient. A flag beside
  # a plausible-looking median still gets read as an estimate, and some of these medians are
  # not merely uninformative but actively misleading: the shore fits' `r_OSP` carried a
  # posterior MEAN of 2.5 million in the 2026-08-29 batch, and `kappa_OSP` reports a tidy
  # 3.0 in every fit of every run under the production `osp_scale_is_tau = TRUE`, which reads
  # exactly like a measured OSP scale and is its lognormal(log 3, 0.3) prior.
  #
  # The raw posterior summary stays: a decoupled parameter that does NOT match its prior is a
  # bug worth seeing. But `estimate` answers the question a reader is actually asking, and it
  # answers NA when there is no estimate to give.
  out$estimate <- ifelse(out$decoupled, NA_real_, out$median)

  # Provenance travels with the row. A parameter table lifted out of a folder whose fit the
  # gate REJECTED looks identical to one the gate accepted; the 2026-08-29 gear stage E run
  # is the worked example, where a rejected boat fit left tau_bar = 2.552 (n_eff 49) and
  # f_crab = 0.315 (n_eff 25, R-hat 1.17) sitting in a file with nothing to say so.
  out$fit_method <- fit_method %||% NA_character_

  out <- out[, c("parameter", "estimate", "median", "lo95", "hi95", "mean", "n_eff", "Rhat",
                 "decoupled", "decoupled_reason", "fit_method")]
  out
}

# --- 2. Divergence localization (where do divergences sit?) ------------------
bss_divergence_localization <- function(fit, candidate_pars = NULL) {
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- unlist(lapply(sp, function(x) x[, "divergent__"])) == 1
  n_div <- sum(divergent); n_tot <- length(divergent)
  res <- list(n_divergent = n_div, n_total = n_tot,
              frac_divergent = n_div / n_tot,
              distortion_C = NA_real_, distortion_E = NA_real_, table = NULL)
  if (n_div == 0 || n_div == n_tot) return(res)

  if (is.null(candidate_pars))
    candidate_pars <- c("sigma_eps_E", "sigma_eps_C", "sigma_r_E", "sigma_r_C",
                        "r_E", "r_C", "phi_E", "phi_C", "sigma_mu_E",
                        "sigma_mu_C", "sigma_IE", "R_T", "R_G_boat", "R_G")
  candidate_pars <- candidate_pars[candidate_pars %in% fit@model_pars]
  flat <- function(p) as.vector(as.array(fit, pars = p)[, , 1])  # chain-major

  rows <- lapply(candidate_pars, function(p) {
    x <- flat(p); xd <- x[divergent]; xn <- x[!divergent]; s <- stats::sd(x)
    smd <- if (is.finite(s) && s > 0) (mean(xd) - mean(xn)) / s else NA_real_
    data.frame(param = p, median_bulk = stats::median(xn),
               median_div = stats::median(xd), smd = smd, abs_smd = abs(smd))
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(-tab$abs_smd), c("param", "median_bulk", "median_div", "smd")]
  res$table <- tab

  for (q in c("C_expected_sum", "E_sum")) {
    if (q %in% fit@model_pars) {
      v <- flat(q); mb <- stats::median(v[!divergent]); md <- stats::median(v[divergent])
      sh <- abs(md - mb) / mb
      if (q == "C_expected_sum") res$distortion_C <- sh else res$distortion_E <- sh
    }
  }
  res
}

# --- 3. Posterior predictive calibration ------------------------------------
# For effort counts (NB(lambda_E * R, r_E)) and interview catches
# (NB(lambda_C * h, r_C)): draw from the posterior predictive for each
# observation, then report central-interval coverage and the PIT. Coverage near
# nominal and uniform PITs indicate the observation model fits. All in R from
# extracted quantities (no RNG added to Stan), capped at n_draws_use for speed.
bss_ppc_calibration <- function(fit, stan_data, n_draws_use = 400, seed = 1) {
  set.seed(seed)
  # Model-agnostic trailer expansion: current pooled (v7.6+) and gear-resolved both
  # carry R_G_boat; pre-v7.6 pooled fits carry R_T. bss_extract_pars() requests only
  # the one this fit declares.
  trailer_par <- bss_trailer_par(fit)
  ex <- rstan::extract(fit, pars = bss_extract_pars(fit, c("lambda_E_S", "lambda_C_S",
                                                           "r_E", "r_C", "R_G")))
  ndraw <- length(ex$r_E)
  use <- if (ndraw > n_draws_use) sort(sample.int(ndraw, n_draws_use)) else seq_len(ndraw)
  nd  <- length(use)
  # lambda_*_S is array[S] matrix[D,G]; rstan returns [iter, S, D, G]. get_lam
  # reduces it to [draws, D] for the S=1, G=1 pooled model and guards against
  # dropped size-1 dims so a shape change cannot abort the whole PPC.
  get_lam <- function(arr) {
    d <- dim(arr)
    m <- if (length(d) == 4) arr[use, 1, , 1]
         else if (length(d) == 3) arr[use, , 1]
         else if (length(d) == 2) arr[use, ]
         else stop(sprintf("unexpected lambda dims: %s", paste(d, collapse = "x")))
    matrix(m, nrow = nd)
  }
  lamE <- get_lam(ex$lambda_E_S)   # [draws, D]
  lamC <- get_lam(ex$lambda_C_S)

  # ROOT-CAUSE FIX (supersedes the earlier B1.7 attempt, which hardened lambda
  # extraction and the rnbinom NA path but did not touch the scalars and so left
  # the failure in place: every fit still aborted with "non-conformable arrays").
  # rstan::extract(permuted = TRUE) returns a SCALAR parameter (r_E, r_C, R_G,
  # R_T) as a 1-D ARRAY -- it carries a length-1 `dim` attribute, it is not a
  # bare vector. Single-bracket indexing of a 1-D array PRESERVES that dim
  # (dim(ex$R_G[use]) is length(use), not NULL), unlike a plain vector (no dim)
  # or an [iter,1] matrix (single-index drops the dim). The predictive means
  # below then evaluate `lamE[, day, drop=FALSE] * RG`, multiplying a 2-D array
  # by a 1-D array: both operands have a `dim`, the dims differ, and R throws
  # "non-conformable arrays". This fired on the gear branch (shore fits, R_G)
  # and the trailer branch (boat fit, R_T) -- exactly the fits that failed --
  # while the catch branch never did, because sweep() rebuilds STATS to match
  # dim(x) and cannot raise this error. as.numeric() strips the stray dim; the
  # multiply then recycles column-wise with the correct per-draw scaling
  # (verified: element [i, j] = lamE[i, day[j]] * scalar[i]). as.numeric() is a
  # no-op on an already-bare vector, so the fix is safe across rstan versions.
  # A two-line check on any fit: dim(rstan::extract(fit, "R_G")$R_G) returns the
  # iteration count (length-1 dim), not NULL.
  rE <- as.numeric(ex$r_E[use]); rC <- as.numeric(ex$r_C[use])
  RG <- as.numeric(ex$R_G[use])
  # RT is the per-draw MULTIPLIER m with mu_trailer = lambda_E * m:
  #   pooled        m = R_T
  #   gear-resolved m = 1 / R_G_boat
  # NULL when the fit declares no trailer expansion parameter.
  RT <- bss_trailer_multiplier(ex, trailer_par, use)

  # B1.7 fix: score each observation on its finite predictive draws only. An
  # extreme lambda draw (exp() overflow in a weakly-identified fit) yields a
  # non-finite mu, and rnbinom(mu = Inf) returns NA, which previously aborted
  # quantile() and the entire PPC. Non-finite mu and non-finite draws are now
  # dropped; an observation with < 20 usable draws is recorded NA and excluded.
  # 2026-09-04: `theta` is NULL for every stream except the catch stream of a zi_catch = 1
  # fit. NULL is the pre-2026-09-04 NB2 arithmetic unchanged, so ppc_calibration_*.csv is
  # bit-reproducible for every earlier run. The mixture arithmetic lives in
  # 03_R_functions/zinb_ppc.R, shared with pit_block() in save_run_diagnostics.R so the
  # aggregate and per-observation files cannot drift apart again.
  calib <- function(mu_mat, y, size_vec, theta = NULL) {
    nobs <- length(y); cov50 <- cov95 <- pit <- rep(NA_real_, nobs)
    for (i in seq_len(nobs)) {
      mu_i <- pmax(mu_mat[, i], 1e-8)
      keep <- is.finite(mu_i) & is.finite(size_vec)
      if (sum(keep) < 20) next
      mu_k <- mu_i[keep]; sz_k <- size_vec[keep]
      th_k <- if (is.null(theta)) NULL else theta[keep]
      # 2026-09-02: the EXACT expectation, matching ppc_byobs_*.csv, instead of simulating.
      # The simulated version is an unbiased estimate of the same quantity but carries Monte
      # Carlo noise, which left the two files disagreeing by 0.005-0.015 even after the
      # randomized-coverage fix. Same formula, no sampling, so the two now agree to floating
      # point and any future disagreement is a real defect rather than noise to be eyeballed.
      # The rnbinom draw is still taken, because the < 20 usable-draw guard below is what
      # protects against an exp() overflow in a weakly-identified fit producing a non-finite mu.
      yp <- stats::rnbinom(sum(keep), mu = mu_k, size = sz_k)
      yp <- yp[is.finite(yp)]
      if (length(yp) < 20) next
      pit[i]   <- bss_zi_pit(y[i], mu_k, sz_k, th_k)
      # 2026-09-01 FIX. Coverage is now read off the RANDOMIZED PIT, not off a quantile
      # interval of the simulated draws.
      #
      # WHAT WAS WRONG. `y[i] >= quantile(yp, .25) && y[i] <= quantile(yp, .75)` asks whether
      # the observation falls inside the central simulated interval. For CONTINUOUS data that
      # is a 50% interval. For small COUNTS it is not and cannot be: the quantiles snap to
      # integers, so the interval carries much more than half the probability mass and the
      # statistic over-covers by construction, with the inflation growing as the counts get
      # smaller. It is a property of the arithmetic, not of the model.
      #
      # WHAT IT COST. On the private-boat trailer stream (fitted means around 1-2 boats) the
      # two computations disagree by up to 0.154, and the 2026-08-31 Stage 5 review recorded
      # "the trailer stream is over-covered in EVERY configuration, 0.667-0.692 against a
      # nominal 0.500, 4.6-5.3 sampling SDs" as a new open modelling item. Under the
      # randomized statistic the same runs give 0.523 and 0.538, which is 0.6 and 1.1
      # sampling SDs: calibrated. The open item was an artefact of this line. The CATCH
      # stream, whose counts are larger, barely moves (0.595 either way), which is exactly
      # the signature of a discreteness effect rather than a model one.
      #
      # WHAT SURVIVES. The daily-AR cells still fail badly under the corrected statistic
      # (0.713 and 0.744 against 0.500, 5.9 and 6.8 sampling SDs), so the Stage 5 conclusion
      # is unchanged and now rests on a statistic that is not biased by count size.
      #
      # This matches ppc_byobs_*.csv, which has always used the randomized PIT, so the two
      # files now agree instead of disagreeing by up to 0.15 on the same quantity.
      cov50[i] <- pit[i] >= 0.25  && pit[i] <= 0.75
      cov95[i] <- pit[i] >= 0.025 && pit[i] <= 0.975
    }
    usable <- is.finite(pit)
    list(summary = data.frame(coverage_50 = mean(cov50, na.rm = TRUE),
                              coverage_95 = mean(cov95, na.rm = TRUE),
                              pit_mean = mean(pit, na.rm = TRUE),
                              pit_sd = stats::sd(pit, na.rm = TRUE),
                              n = sum(usable), n_obs = nobs),
         pit = pit[usable])
  }

  parts <- list()
  if (stan_data$Gear_n > 0)
    parts$gear    <- calib(lamE[, stan_data$day_Gear, drop = FALSE] * RG,
                           stan_data$Gear_I, rE)
  if (stan_data$T_n > 0 && !is.null(RT))
    parts$trailer <- calib(lamE[, stan_data$day_T, drop = FALSE] * RT,
                           stan_data$T_I, rE)
  # improvement 2.2 (2026-08-27): the OSP stream was the one observation stream with NO
  # posterior-predictive check, which is why the scale conflict it carries had to be
  # inferred from the TRAILER PIT rather than read off directly. The OSP mean mirrors the
  # Stan likelihood exactly: (lambda_E / R_G_boat) * (osp_scale_is_tau ? L : kappa_OSP),
  # with its own dispersion r_OSP. A PIT mean well below 0.5 here, alongside the trailer's,
  # is the direct read on whether lambda_E is being pulled between two streams that disagree.
  if ((stan_data$OSP_n %||% 0) > 0 && !is.null(RT)) {
    # NOTE the draw subsetting: lamE, RT and the dispersion vectors are all on the `use`
    # subset, so anything pulled fresh out of the fit has to be subset the same way or the
    # element-wise product silently recycles.
    osp_scale <- if (identical(as.integer(stan_data$osp_scale_is_tau %||% 0L), 1L)) {
      Lx <- try(rstan::extract(fit, pars = "L_out")$L_out, silent = TRUE)
      if (inherits(Lx, "try-error") || is.null(Lx)) NULL
      else Lx[use, stan_data$day_OSP, drop = FALSE]
    } else {
      kx <- try(rstan::extract(fit, pars = "kappa_OSP")$kappa_OSP, silent = TRUE)
      if (inherits(kx, "try-error") || is.null(kx)) NULL
      else matrix(as.numeric(kx)[use], nrow = nd, ncol = length(stan_data$day_OSP))
    }
    rO <- try(as.numeric(rstan::extract(fit, pars = "r_OSP")$r_OSP)[use], silent = TRUE)
    if (!is.null(osp_scale) && !inherits(rO, "try-error") && !is.null(rO) &&
        identical(dim(osp_scale), dim(lamE[, stan_data$day_OSP, drop = FALSE])))
      parts$osp <- calib(lamE[, stan_data$day_OSP, drop = FALSE] * RT * osp_scale,
                         stan_data$OSP_I, rO)
  }
  if (stan_data$IntC > 0)
    parts$catch   <- calib(sweep(lamC[, stan_data$day_IntC, drop = FALSE], 2,
                                 stan_data$h, "*"),
                           stan_data$c, rC,
                           bss_zi_theta_draws(fit, stan_data, use))
  if (length(parts) == 0) return(NULL)

  summ <- do.call(rbind, lapply(names(parts), function(nm)
    cbind(data_type = nm, parts[[nm]]$summary)))
  pit_long <- do.call(rbind, lapply(names(parts), function(nm)
    data.frame(data_type = nm, pit = parts[[nm]]$pit)))
  list(summary = summ, pit = pit_long)
}

# --- 4. Write all diagnostics for one fit -----------------------------------
write_bss_diagnostics <- function(fit, stan_data, label, output_dir, fit_method = NULL) {
  ok <- function(expr) tryCatch(expr, error = function(e) {
    cat(sprintf("    [diag] %s skipped: %s\n", label, conditionMessage(e))); NULL })

  ok({
    sp <- bss_structural_summary(fit, stan_data, fit_method)
    utils::write.csv(sp, file.path(output_dir, sprintf("structural_params_%s.csv", label)),
                     row.names = FALSE)
  })

  ok({
    dl <- bss_divergence_localization(fit)
    cat(sprintf("    [diag] %s: %d/%d divergent (%.1f%%); distortion C/E = %s/%s\n",
                label, dl$n_divergent, dl$n_total, 100 * dl$frac_divergent,
                ifelse(is.na(dl$distortion_C), "NA", sprintf("%.1f%%", 100 * dl$distortion_C)),
                ifelse(is.na(dl$distortion_E), "NA", sprintf("%.1f%%", 100 * dl$distortion_E))))
    if (!is.null(dl$table)) {
      dl$table$distortion_C <- dl$distortion_C
      dl$table$distortion_E <- dl$distortion_E
      utils::write.csv(dl$table,
                       file.path(output_dir, sprintf("divergence_localization_%s.csv", label)),
                       row.names = FALSE)
    }
  })

  ok({
    if (is.null(stan_data)) {
      cat("    [diag] PPC skipped: stan_data not stored with this fit\n")
    } else {
    ppc <- bss_ppc_calibration(fit, stan_data)
    if (!is.null(ppc)) {
      utils::write.csv(ppc$summary,
                       file.path(output_dir, sprintf("ppc_calibration_%s.csv", label)),
                       row.names = FALSE)
      cat("    [diag] PPC coverage (nominal 50/95):\n")
      for (k in seq_len(nrow(ppc$summary)))
        cat(sprintf("      %-8s 50%%=%.0f%%  95%%=%.0f%%  (n=%d)\n",
                    ppc$summary$data_type[k], 100 * ppc$summary$coverage_50[k],
                    100 * ppc$summary$coverage_95[k], ppc$summary$n[k]))
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        p <- ggplot2::ggplot(ppc$pit, ggplot2::aes(x = pit)) +
          ggplot2::geom_histogram(boundary = 0, bins = 10, fill = "steelblue",
                                  colour = "white") +
          ggplot2::geom_hline(yintercept = 0, colour = NA) +
          ggplot2::facet_wrap(~data_type, scales = "free_y") +
          ggplot2::labs(title = sprintf("PPC PIT: %s", label),
                        subtitle = "Uniform => calibrated; U-shape => underdispersed; hump => overdispersed",
                        x = "PIT", y = "count") +
          ggplot2::theme_bw()
        ggplot2::ggsave(file.path(output_dir, sprintf("ppc_pit_%s.png", label)),
                        p, width = 10, height = 4)
      }
    }
    }
  })
  invisible(NULL)
}
