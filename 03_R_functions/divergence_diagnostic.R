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
# divergence_diagnostic.R
# Localize the source of divergent transitions in a fitted pooled BSS, and test
# whether the divergences actually distort the reported totals.
#
# Two reparameterizations (B1.3 non-centered omega_0, B1.5 marginalized
# overdispersion) have not cleared the shore divergences. Rather than guess at a
# third, this finds which parameter the divergent draws cluster on (the funnel
# neck), per Betancourt (2017, arXiv:1701.02434).
#
# Usage (in the same R session right after the pooled run, fit objects still
# live in bss_all):
#
#   source(here::here("03_R_functions/divergence_diagnostic.R"))
#   diagnose_divergences(bss_all[["shore_ring_net_only_Dungeness_Kept"]]$fit)
#   diagnose_divergences(bss_all[["shore_all_gear_Dungeness_Kept"]]$fit)
#   diagnose_divergences(bss_all[["private_boat_all_gear_Dungeness_Kept"]]$fit)
#
# If the session was cleared, re-fit just shore_ring_net (the fastest, ~40 min)
# and pass that one fit. You do not need a full re-run.
# =============================================================================

diagnose_divergences <- function(fit, candidate_pars = NULL) {
  stopifnot(inherits(fit, "stanfit"))

  # --- Per-draw divergence flag, post-warmup, flattened chain-major -----------
  # get_sampler_params() returns one matrix per chain (chain-major). as.array()
  # returns [iter, chain, par]; as.vector() of the [,,1] slice is column-major =
  # chain-major, so the two align draw-for-draw.
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- unlist(lapply(sp, function(x) x[, "divergent__"])) == 1
  n_div <- sum(divergent); n_tot <- length(divergent)
  cat(sprintf("\nDivergences: %d / %d post-warmup draws (%.1f%%)\n",
              n_div, n_tot, 100 * n_div / n_tot))
  if (n_div == 0) { cat("No divergences to localize.\n"); return(invisible(NULL)) }
  if (n_div < 10) cat("  (fewer than 10; the ranking below is noisy)\n")

  # --- Candidate neck parameters (scalars that scale high-dim blocks) ---------
  if (is.null(candidate_pars))
    candidate_pars <- c("sigma_eps_E", "sigma_eps_C",   # AR innovation SDs
                        "sigma_r_E", "sigma_r_C",       # NB dispersion scale
                        "r_E", "r_C",                   # NB dispersion (= 1/sigma_r^2)
                        "phi_E", "phi_C",               # AR persistence (stationary var ~ 1/(1-phi^2))
                        "sigma_mu_E", "sigma_mu_C",     # hierarchical intercept SDs
                        "sigma_IE",                     # I/E lognormal SD (shore)
                        "R_T", "R_G")                   # expansion ratios (boundary check)
  candidate_pars <- candidate_pars[candidate_pars %in% fit@model_pars]

  flat <- function(p) as.vector(as.array(fit, pars = p)[, , 1])  # chain-major

  # --- Rank by how far divergent draws sit from the bulk ----------------------
  # Standardized mean difference on the natural scale, plus the raw medians so
  # the direction is obvious (e.g. divergent sigma_eps far BELOW bulk = the neck
  # is at small sigma_eps; divergent phi far ABOVE bulk = stationary-variance
  # blow-up near phi = 1).
  rows <- lapply(candidate_pars, function(p) {
    x  <- flat(p); xd <- x[divergent]; xn <- x[!divergent]
    s  <- stats::sd(x)
    smd <- if (is.finite(s) && s > 0) (mean(xd) - mean(xn)) / s else NA_real_
    data.frame(param = p,
               median_bulk = stats::median(xn),
               median_div  = stats::median(xd),
               smd         = smd, abs_smd = abs(smd))
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(-tab$abs_smd), ]
  cat("\nWhere do divergences sit, relative to the non-divergent bulk?\n")
  cat("Large |smd| => divergences concentrate at an extreme of this parameter,\n")
  cat("which is the funnel neck and the thing to reparameterize next.\n\n")
  print(tab[, c("param","median_bulk","median_div","smd")], row.names = FALSE, digits = 4)

  # --- Distortion check: do divergent draws move the reported totals? ---------
  # If divergent and non-divergent draws give the same C/E, the divergences are
  # not biasing the estimate, and using the (otherwise clean: Rhat 1.00,
  # n_eff > 4000, treedepth 0) BSS with a documented caveat is defensible
  # against falling back to the cruder PE.
  cat("\nDistortion check (does excluding divergent draws change the estimate?):\n")
  for (q in c("C_expected_sum", "E_sum")) {
    if (q %in% fit@model_pars) {
      v <- flat(q); mb <- stats::median(v[!divergent]); md <- stats::median(v[divergent])
      cat(sprintf("  %-15s bulk median = %10.0f   divergent median = %10.0f   (%+.1f%%)\n",
                  q, mb, md, 100 * (md - mb) / mb))
    }
  }
  cat("\nIf both shifts are small (say < a few %), the divergences are not\n")
  cat("distorting the totals; the gate, not the geometry, is what is forcing PE.\n")
  invisible(tab)
}
