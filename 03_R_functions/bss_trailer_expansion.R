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
# bss_trailer_expansion.R
#
# Model adapter for the shared diagnostics. Two trailer-expansion conventions exist
# across model versions:
#
#   crab_bss_pooled.stan (<= v7.5)   T_I[i] ~ NB2(lambda_E[d_i] * R_T,      r_E)
#   crab_bss_pooled.stan (v7.6+)     T_I[i] ~ NB2(lambda_E[d_i] / R_G_boat, r_E)
#   crab_bss_gear_resolved.stan      T_I[i] ~ NB2(lambda_E[d_i] / R_G_boat, r_E)
#
# R_T was trailers-per-crabber (in [0,1]); R_G_boat is gear-per-boat-group (~4).
# They are NOT a rename of one another: one multiplies, the other divides. Pooled
# v7.6 (POOL-1) adopted the gear-resolved / R_G_boat convention, so current fits of
# BOTH models carry R_G_boat; this adapter still handles R_T for pre-v7.6 fits.
#
# Every shared diagnostic that reconstructs the trailer predictive mean (PPC
# calibration, per-observation PIT, effort overdispersion decomposition) needs
# mu_trailer = lambda_E * m for a per-draw multiplier m. These helpers supply
# that m for whichever model produced the fit:
#
#     m = R_T            (pooled <= v7.5 fits)
#     m = 1 / R_G_boat   (pooled v7.6+ and gear-resolved fits)
#
# so every downstream `lamE * m` expression stays correct and unchanged.
###############################################################################

# Which trailer expansion parameter does this fit carry? NA if neither (e.g. a
# shore fit whose trailer stream is empty, or a future model).
bss_trailer_par <- function(fit) {
  mp <- fit@model_pars
  if ("R_T" %in% mp) "R_T"
  else if ("R_G_boat" %in% mp) "R_G_boat"
  else NA_character_
}

# Base parameter vector plus whichever trailer parameter exists. Use this to
# build the `pars` argument of rstan::extract() so extraction never requests a
# parameter the model does not declare (which errors).
bss_extract_pars <- function(fit, base) {
  tp <- bss_trailer_par(fit)
  if (is.na(tp)) base else c(base, tp)
}

# Per-draw MULTIPLIER m such that mu_trailer = lambda_E * m, given an extracted
# draws list `ex` and a draw index vector `idx`. Returns NULL when the fit has no
# trailer expansion parameter, so callers can skip the trailer stream cleanly.
#
# rstan::extract(permuted = TRUE) returns scalar parameters as 1-D arrays that
# retain a length-1 `dim`; as.numeric() strips it so column-wise recycling in
# `lamE * m` is correct (this is the same non-conformable-arrays trap documented
# at length in model_diagnostics.R).
bss_trailer_multiplier <- function(ex, par_name, idx) {
  if (is.null(par_name) || is.na(par_name) || is.null(ex[[par_name]])) return(NULL)
  v <- as.numeric(ex[[par_name]][idx])
  if (identical(par_name, "R_T")) v else 1 / v
}
