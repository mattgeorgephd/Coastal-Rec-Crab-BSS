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
# bss_predictive_catch.R  (2026-09-28, B46; shared by both drivers)
#
# The PREDICTIVE season catch: what the realised harvest could have been given the fitted
# rates, drawn from the fitted OBSERVATION model rather than from a Poisson.
#
# WHY. Both Stan models draw the predictive daily catch as poisson_rng(lambda_Ctot), which
# carries neither the catch overdispersion r_C nor the zero-inflation theta_C that the
# likelihood fits, so Predictive_Catch was under-dispersed relative to the model's own
# observation layer (method document Section 14.8 and Section 20). Rebuilding it in R from
# the saved draws leaves the Stan code and every fit byte-identical.
#
# THE CONSTRUCTION. An interview is one party's trip: c_a ~ ZINB(theta, mean lambda_C * h_a,
# size r_C). A day's harvest is the sum over the n_d parties that fished it, n_d = E_d / hbar
# (E_d the day's effort in the fit's unit, hbar the mean h over the fitted interviews). Given
# the day's rate, the parties are independent, so the day total has
#   mean  mu_d = E_d * lambda_C * (1 - theta)            (= the expected catch C_expected)
#   var   mu_d * (1 + m_d / r_C + theta * m_d),  m_d = lambda_C * hbar (a party's NB mean)
# and it is drawn as NB2 with that mean and variance: size k_d = n_d (1 - theta) / (1/r_C + theta).
# With theta = 0 this is exactly the sum of n_d NB2(r_C) parties, NB2(n_d r_C). Days sum
# within a draw. The expected catch (the headline) is unchanged; only the predictive
# interval widens, by the party-level noise the Poisson left out.
#
# Returns one predictive season total per posterior draw, in rstan::extract()'s draw order
# (the order C_expected_sum / C_sum come in), or NULL when the fit lacks what it needs.
# Seeded, with the caller's RNG restored (bss_with_seed()).
###############################################################################
bss_predictive_catch <- function(fit, stan_data, seed = 1L) {
  if (is.null(fit) || is.null(stan_data)) return(NULL)
  pn <- fit@model_pars
  mu_par <- if ("C_expected" %in% pn) "C_expected" else if ("C_total" %in% pn) "C_total" else return(NULL)
  if (!all(c("E", "r_C") %in% pn)) return(NULL)
  ex <- rstan::extract(fit, pars = intersect(c(mu_par, "E", "r_C", "theta_C_out"), pn))
  # reduce [draws, S, D, G] (or [draws, D]) to [draws, D] by summing sections and gears
  to_dd <- function(x) {
    d <- dim(x)
    if (length(d) == 2) return(x)
    if (length(d) == 4) return(apply(x, c(1, 3), sum))
    NULL
  }
  mu <- to_dd(ex[[mu_par]]); E <- to_dd(ex$E)
  if (is.null(mu) || is.null(E) || !identical(dim(mu), dim(E))) return(NULL)
  hbar <- mean(as.numeric(stan_data$h), na.rm = TRUE)
  if (!is.finite(hbar) || hbar <= 0) return(NULL)
  r  <- as.numeric(ex$r_C)
  th <- if (!is.null(ex$theta_C_out)) as.numeric(ex$theta_C_out) else rep(0, length(r))
  bss_predictive_draws(mu, E, r, th, hbar, seed)
}

# The pure part: mu and E are [draws, days] matrices (expected catch and effort), r and th one
# value per draw (r_C, theta_C), hbar the mean gear per interviewed party. Returns one
# predictive season total per draw.
bss_predictive_draws <- function(mu, E, r, th, hbar, seed = 1L) {
  n_d  <- E / hbar
  size <- n_d * ((1 - th) / (1 / r + th))          # recycled down the draws (rows)
  ok   <- is.finite(mu) & mu > 0 & is.finite(size) & size > 0
  bss_with_seed(seed, {
    y <- matrix(0, nrow(mu), ncol(mu))
    y[ok] <- stats::rnbinom(sum(ok), size = size[ok], mu = mu[ok])
    rowSums(y)
  })
}

# One seed per fit, from the run seed and the fit's label, so two fits' predictive noise is
# never drawn from the same stream (the port total pairs the fits' draws row by row).
bss_predictive_seed <- function(base_seed, label) {
  as.integer((as.numeric(base_seed %||% 1L) + sum(utf8ToInt(as.character(label)))) %% .Machine$integer.max)
}
