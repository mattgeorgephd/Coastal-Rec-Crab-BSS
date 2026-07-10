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
#     population  unit               h                  E_scale  L
#     ----------  -----------------  -----------------  -------  ------------------
#     shore       crabber-hours      fishing_time_total    1     L_eff (hours)
#     shore       gear-hours         gear_time_total     R_G     L_eff (hours)
#     shore       gear-deployments   number_of_gear      R_G     tau_shore
#     boat        gear-deployments   number_of_gear        1     tau_boat
#
#   In every row E and h carry the same unit, which is the invariant that
#   bss_assert_effort_units() checks before sampling.
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
      L_prior_sigma     = rep(params$tau_boat_prior_sigma %||% 0.3, D)
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
      L_prior_sigma     = days$L_prior_sigma
    ),
    "gear-hours" = list(
      unit              = unit,
      h_col             = "gear_time_total",
      h_fun             = function(int_d) .num(int_d, "gear_time_total"),
      effort_scale_gear = 1L,                 # crabbers -> gear via R_G
      L_data            = days$L_mu,
      L_prior_sigma     = days$L_prior_sigma
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
      L_prior_sigma     = rep(params$tau_shore_prior_sigma %||% 0.3, D)
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
