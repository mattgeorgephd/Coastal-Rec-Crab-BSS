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
# bss_convergence_gate.R
#
# Shared scale-aware convergence gate (B1.8) for the crab BSS models. Extracted
# from the pooled driver so the pooled and gear-resolved tracks use one gate and
# do not drift. Auto-sourced by both drivers via the 03_R_functions walk.
#
# Two functions:
#   bss_compute_gate(b, label, ...) -> one-row tibble of gate diagnostics and the
#     pass/fail decision for a single fit.
#   bss_use_pe_for(b) -> TRUE if the fit should report PE (insufficient data OR
#     failed the gate), FALSE if it should report BSS.
#
# The gate decides PE-vs-BSS on THREE things beyond R-hat/n_eff:
#   1. divergence fraction below a hard backstop (max_div_fraction), and
#   2. the divergent draws do not move the reported totals by more than
#      max_impact_sd, measured in posterior-SD units (scale-invariant), on
#   3. the EXPECTED total (catch_par), not the Poisson-predictive realization.
#
# Design notes for dual-model use:
#   * catch_par / effort_par name the summed quantities to gate on. Both models
#     expose C_expected_sum and E_sum (gear-resolved got C_expected_sum in the
#     stage-1/2 Stan work); pass those.
#   * Thresholds are arguments (not globals) so each driver passes its own params.
#   * ar_resolution is an optional label for the report column (pooled passes its
#     adaptive AR resolution; gear-resolved uses fixed periods and passes NA).
#   * The function reads only b$fit, b$pe_fallback. It does not depend on any
#     model-specific object, so it is safe to call from either driver.
#
# Reference: Vehtari et al. (2021) for the R-hat/ESS thresholds; Betancourt
# (2017) for divergence-as-geometry-failure motivating the impact test.
###############################################################################

bss_compute_gate <- function(b, label,
                             catch_par        = "C_expected_sum",
                             effort_par       = "E_sum",
                             max_impact_sd    = 0.10,
                             max_div_fraction = 0.05,
                             rhat_threshold   = 1.01,
                             neff_threshold   = 400,
                             ar_resolution    = NA_character_) {

  # --- PE-only entry (insufficient data / no fit): nothing to gate. ----------
  if (is.null(b$fit) || (!is.null(b$pe_fallback) && b$pe_fallback)) {
    return(tibble(
      fit = label, ar_resolution = NA_character_,
      divergences = NA_integer_, divergence_fraction = NA_real_, treedepth_pct = NA_real_,
      C_sum_rhat = NA_real_, E_sum_rhat = NA_real_, C_sum_neff = NA_real_, E_sum_neff = NA_real_,
      impact_C_sd = NA_real_, impact_E_sd = NA_real_,
      distortion_C = NA_real_, distortion_E = NA_real_,
      pass_rhat = FALSE, pass_neff = FALSE, pass_div_fraction = FALSE, pass_impact = FALSE,
      pass_convergence = FALSE, method_selected = "PE (insufficient data)"
    ))
  }

  sp <- rstan::get_sampler_params(b$fit, inc_warmup = FALSE)
  n_divergent <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  n_treedepth <- sum(sapply(sp, function(x) sum(x[, "treedepth__"] >=
    b$fit@stan_args[[1]]$control$max_treedepth)))
  n_total  <- sum(sapply(sp, nrow))
  div_frac <- if (n_total > 0) n_divergent / n_total else NA_real_

  summ   <- summary(b$fit, pars = c(catch_par, effort_par))$summary
  rhat_C <- summ[catch_par,  "Rhat"];  rhat_E <- summ[effort_par, "Rhat"]
  neff_C <- summ[catch_par,  "n_eff"]; neff_E <- summ[effort_par, "n_eff"]

  pass_rhat <- max(rhat_C, rhat_E, na.rm = TRUE) < rhat_threshold
  pass_neff <- min(neff_C, neff_E, na.rm = TRUE) > neff_threshold

  # Divergent-draw flag aligned to the flattened draw vector. get_sampler_params
  # returns a per-chain list; unlist() concatenates chain-major, and
  # as.array(fit)[,,1] flattened with as.numeric() is also chain-major, so the
  # flag and the draws line up element-for-element.
  divergent_flag <- unlist(lapply(sp, function(x) x[, "divergent__"])) == 1
  flat_par <- function(p) as.numeric(as.array(b$fit, pars = p)[, , 1])

  # Scale-aware impact: |median(all) - median(non-divergent)| / sd(all). In units
  # of posterior SD, so it does not penalize a wide posterior; it asks whether
  # the divergences move the answer relative to how well it is pinned down.
  impact_sd <- function(p) {
    v  <- flat_par(p)
    vb <- v[!divergent_flag]
    s  <- sd(v)
    if (!is.finite(s) || s <= 0 || length(vb) == 0) return(NA_real_)
    abs(median(v) - median(vb)) / s
  }
  # Legacy level-distortion: reported for continuity, NOT gating. Scales with
  # posterior width, so it is a reported column only.
  level_distortion <- function(p) {
    v <- flat_par(p)
    if (sum(divergent_flag) == 0 || sum(!divergent_flag) == 0) return(NA_real_)
    mb <- median(v[!divergent_flag])
    if (!is.finite(mb) || mb == 0) return(NA_real_)
    abs(median(v[divergent_flag]) - mb) / mb
  }

  if (n_divergent > 0) {
    impact_C <- impact_sd(catch_par); impact_E <- impact_sd(effort_par)
    distortion_C <- level_distortion(catch_par); distortion_E <- level_distortion(effort_par)
  } else {
    impact_C <- 0; impact_E <- 0
    distortion_C <- 0; distortion_E <- 0
  }

  pass_impact <- is.finite(impact_C) && is.finite(impact_E) &&
                 max(impact_C, impact_E) < max_impact_sd
  pass_div_fraction <- is.finite(div_frac) && (div_frac < max_div_fraction)

  # Final gate: chains mixed (R-hat), enough effective draws (n_eff), divergence
  # fraction under the hard backstop, AND divergences do not move the totals.
  pass <- pass_rhat & pass_neff & pass_div_fraction & pass_impact

  # Soft warnings (do not affect pass/fail) so the run log still flags geometry.
  if (n_total > 0 && n_treedepth / n_total > 0.05) {
    cat(sprintf("  WARNING: %s - %.1f%% treedepth exceedances (consider increasing max_treedepth)\n",
                label, n_treedepth / n_total * 100))
  }
  if (n_divergent > 0) {
    cat(sprintf("  %s: %d/%d divergent (%.1f%%); impact C/E = %.4f/%.4f SD; level-distortion C/E = %.1f%%/%.1f%%\n",
                label, n_divergent, n_total, div_frac * 100,
                impact_C, impact_E, distortion_C * 100, distortion_E * 100))
  }
  if (!pass_div_fraction && is.finite(div_frac)) {
    cat(sprintf("  WARNING: %s - divergence fraction %.1f%% exceeds backstop %.0f%%; forcing PE.\n",
                label, div_frac * 100, max_div_fraction * 100))
  }
  if (pass_rhat && pass_neff && pass_div_fraction && !pass_impact) {
    cat(sprintf("  WARNING: %s - divergences move the totals %.3f SD (C) / %.3f SD (E), over %.2f SD; forcing PE.\n",
                label, impact_C, impact_E, max_impact_sd))
  }

  tibble(
    fit = label,
    ar_resolution = ar_resolution,
    divergences = n_divergent,
    divergence_fraction = round(div_frac, 4),
    treedepth_pct = round(n_treedepth / n_total * 100, 1),
    C_sum_rhat = round(rhat_C, 4),
    E_sum_rhat = round(rhat_E, 4),
    C_sum_neff = round(neff_C),
    E_sum_neff = round(neff_E),
    impact_C_sd = round(impact_C, 4),
    impact_E_sd = round(impact_E, 4),
    distortion_C = round(distortion_C, 4),
    distortion_E = round(distortion_E, 4),
    pass_rhat = pass_rhat,
    pass_neff = pass_neff,
    pass_div_fraction = pass_div_fraction,
    pass_impact = pass_impact,
    pass_convergence = pass,
    method_selected = if (pass) "BSS" else "PE (convergence fail)"
  )
}

# A fit reports PE either because it never fit (insufficient data) or because it
# fit but failed the scale-aware gate. Reported-estimate sections call this to
# choose PE vs BSS; diagnostic sections still key off b$pe_fallback directly so a
# fitted-but-gate-failed component keeps its per-fit diagnostics.
bss_use_pe_for <- function(b) isTRUE(b$pe_fallback) || !isTRUE(b$use_bss)
