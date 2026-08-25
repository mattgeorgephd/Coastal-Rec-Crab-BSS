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
# bss_effort_spec.R
#
# P1: single source of truth for the effort unit of a population, used by BOTH
# prep_bss_crab() (which builds the Stan data) and run_pe() (which builds the
# point estimate). Keeping one function means the BSS and the PE can never end
# up on different scales, which is exactly the failure that made the 2026-07-09
# PE-versus-BSS comparison meaningless.
#
# WHY THIS EXISTS
#   c[a] ~ NB2(lambda_C * h[a], r_C) asserts catch is PROPORTIONAL to h. The
#   pipeline measures that assumption every run (cpue_linearity_*.csv). On the
#   2026-07-09 run:
#       boat, gear-deployments : beta_h = 0.754  (95% CI 0.468 to 1.039)  OK
#       shore, crabber-hours   : beta_h = 0.571  (95% CI 0.500 to 0.641)  FAIL
#       shore ring, crabber-hr : beta_h = 0.620  (95% CI 0.530 to 0.710)  FAIL
#   So crabber-hours is not a valid effort unit for shore. Rather than assume the
#   boat's answer transfers, offer the candidates and let PSIS-LOO decide.
#
# THE UNIT ALGEBRA
#   The gear-count likelihood is Gear_I ~ NB2(lambda_E * R_G, r_E), so for SHORE
#   lambda_E is CRABBERS and lambda_E * R_G is GEAR. For the BOAT the trailer
#   likelihood is T_I ~ NB2(lambda_E / R_G_boat, r_E), so lambda_E is GEAR.
#
#   The Stan model forms E = lambda_E * E_scale * L, with
#   E_scale = R_G when effort_scale_gear = 1 and 1 otherwise. Therefore:
#
#     population  unit               h                  E_scale  L            I/E observation
#     ----------  -----------------  -----------------  -------  -----------  ----------------
#     shore       crabber-hours      fishing_time_total    1     L_eff (hrs)  crabber-hours
#     shore       gear-hours         gear_time_total     R_G     L_eff (hrs)  crabber-hours
#     shore       gear-deployments   number_of_gear      R_G     tau_shore    crabber TRIPS
#     boat        gear-deployments   number_of_gear        1     tau_boat     boat trips
#
#   In every row E and h carry the same unit, which is the invariant that
#   bss_assert_effort_units() checks before sampling.
#
# THE I/E OBSERVATION COLUMN (improvement 1/2 fix, 2026-08-25)
#   The last column above is new, and it fixes a live unit error. The shore I/E likelihood
#   is  IE_obs ~ lognormal(log(lambda_E * L), sigma_IE).  Under crabber-hours that was
#   right: lambda_E is crabbers present and L is L_effective in HOURS, so lambda_E * L is
#   crabber-hours, matching ie_crabber_hours. The v7.7 move to gear-deployments replaced L
#   with tau_shore (a dimensionless turnover, ~1.7) but left the observation as
#   crabber-hours, so since v7.7 the model has been comparing crabber-HOURS against a
#   predicted lambda_E * tau = crabber-TRIPS. On a typical WDF20 day that is ~331 observed
#   against ~80 predicted, a ~4x mismatch the gear-count stream then has to absorb; it is
#   the most likely mechanism behind the unexplained shore all-gear sigma_IE ~ 1.07 (GR-9).
#   The BOAT stream was already correct (F2: observation = boat trips, predicted =
#   groups x tau). The spec now names the matching observation column per unit, and both
#   preps read it, so the two can no longer drift.
#
#   Both quantities come from the SAME I/E presence series, so nothing new is measured:
#       L_effective = crabber-hours / peak present   (~5.3 h)
#       turnover    = arrivals      / peak present   (~1.72)
#   and their ratio is the implied trip length, 3.07 h against an interview mean of
#   3.23 h. Production expands on the turnover, so the turnover is what the I/E stream
#   must observe.
#
#   params$ie_shore_obs_unit = "crabber_hours" forces the historical (mismatched under
#   deployments) behaviour, for reproducing pre-2026-08-25 runs.
#
# HOW TO CHOOSE
#   Run the model once per candidate value of params$shore_effort_unit and
#   compare elpd_loo on the CATCH stream in loo_summary_*.csv. The catch
#   observations c[a] are identical across runs; only h[a] and the model change,
#   so the comparison is a valid model comparison. Higher elpd_loo wins.
#
# RETURNS a list:
#   unit               character label, also used for the unit assertion
#   h_col              the interview column used as the CPUE denominator
#   h_fun(int_d)       numeric vector of h for the Stan data
#   effort_scale_gear  0L or 1L, passed to Stan
#   L_data             per-day expansion factor (hours, or a turnover)
#   L_prior_sigma      log-scale SD on L_data
###############################################################################

bss_effort_spec <- function(is_shore, days, params = list()) {

  D <- nrow(days)

  .num <- function(df, nm) {
    if (!nm %in% names(df)) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(df[[nm]]))
  }

  if (!is_shore) {
    # Boat: lambda_E is already gear, so no R_G conversion. L is the deployment
    # turnover tau_boat, identified by WBL I/E ingress counts when available.
    return(list(
      unit              = "gear-deployments",
      h_col             = "number_of_gear",
      h_fun             = function(int_d) .num(int_d, "number_of_gear"),
      effort_scale_gear = 0L,
      L_data            = rep(params$tau_boat_prior_mu    %||% 1.2, D),
      L_prior_sigma     = rep(params$tau_boat_prior_sigma %||% 0.3, D),
      ie_obs_col        = "ie_trips",           # boat ingress count (F2)
      ie_obs_unit       = "boat trips",
      L_unit            = "turnover (trips per present group per day)"
    ))
  }

  unit <- params$shore_effort_unit %||% "crabber-hours"
  allowed <- c("crabber-hours", "gear-hours", "gear-deployments")
  if (!unit %in% allowed) {
    stop("params$shore_effort_unit must be one of: ",
         paste(allowed, collapse = ", "), " (got '", unit, "')", call. = FALSE)
  }

  switch(unit,
    "crabber-hours" = list(
      unit              = unit,
      h_col             = "fishing_time_total",
      h_fun             = function(int_d) .num(int_d, "fishing_time_total"),
      effort_scale_gear = 0L,                 # lambda_E is already crabbers
      L_data            = days$L_mu,          # effective day length, hours
      L_prior_sigma     = days$L_prior_sigma,
      ie_obs_col        = "ie_crabber_hours",
      ie_obs_unit       = "crabber-hours",
      L_unit            = "effective day length (hours)"
    ),
    "gear-hours" = list(
      unit              = unit,
      h_col             = "gear_time_total",
      h_fun             = function(int_d) .num(int_d, "gear_time_total"),
      effort_scale_gear = 1L,                 # crabbers -> gear via R_G
      L_data            = days$L_mu,
      L_prior_sigma     = days$L_prior_sigma,
      ie_obs_col        = "ie_crabber_hours",
      ie_obs_unit       = "crabber-hours",
      L_unit            = "effective day length (hours)"
    ),
    "gear-deployments" = list(
      unit              = unit,
      h_col             = "number_of_gear",
      h_fun             = function(int_d) .num(int_d, "number_of_gear"),
      effort_scale_gear = 1L,                 # crabbers -> gear via R_G
      # tau_shore: trips per gear-slot per day. 30 WDF20 I/E days give
      # arrivals/peak = 1.72 (median 1.69, sd 0.45), and L_eff/tau = 3.06 h
      # against an interview mean trip length of 3.23 h.
      L_data            = rep(params$tau_shore_prior_mu    %||% 1.7, D),
      L_prior_sigma     = rep(params$tau_shore_prior_sigma %||% 0.3, D),
      # improvement 1/2 fix: under deployments the predicted I/E quantity is
      # lambda_E * tau_shore = crabber TRIPS, so the observation must be the crabber
      # ARRIVAL count, not crabber-hours. ie_shore_obs_unit = "crabber_hours" restores
      # the historical (mismatched) pairing.
      ie_obs_col        = if (identical(params$ie_shore_obs_unit %||% "auto", "crabber_hours"))
                            "ie_crabber_hours" else "ie_trips",
      ie_obs_unit       = if (identical(params$ie_shore_obs_unit %||% "auto", "crabber_hours"))
                            "crabber-hours (LEGACY, unit-mismatched under deployments)"
                          else "crabber trips",
      L_unit            = "turnover (trips per gear-slot per day)"
    )
  )
}


# Candidate CPUE denominators for a population. Used to build the common
# interview subset when comparing effort units by PSIS-LOO: elpd_loo is only
# comparable across models fitted to the SAME observations c[a], so every
# candidate denominator must be valid on every retained interview.
bss_effort_h_candidates <- function(is_shore) {
  if (is_shore) c("fishing_time_total", "gear_time_total", "number_of_gear")
  else          c("number_of_gear")
}
