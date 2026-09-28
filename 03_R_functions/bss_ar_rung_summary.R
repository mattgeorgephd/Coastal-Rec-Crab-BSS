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
# ---------------------------------------------------------------------------
# bss_ar_rung_summary()  --  2026-09-04.
#
# WHAT IT IS FOR. The escalation ladder (bss_ar_ladder / params$ar_escalate) refits a
# component at successively coarser AR resolutions until one passes the convergence gate.
# Until now the only thing kept from an attempt was its SAMPLER diagnostics: divergences,
# the divergence fraction, the impact test and the pass/fail. The ESTIMATE each rung
# produced was discarded the moment the loop moved on, because `fit` is overwritten.
#
# That is the wrong thing to throw away. The whole point of running a ladder is to see how
# the answer depends on the resolution, and a reader of the report cannot see that if only
# the surviving rung's number exists. It also makes the tie-break in
# params$ar_escalate_select unimplementable: you cannot prefer the rung with the narrowest
# predictive interval if you did not record any rung's interval.
#
# This extracts the season totals and their intervals from one attempt's fit, cheaply
# (three parameters, no draw extraction), so every rung's answer survives into
# ar_escalation_log.csv and into the rendered report.
#
# ON `pi_width_rel`, AND A WARNING ABOUT USING IT TO CHOOSE. It is (hi95 - lo95) / median,
# a coefficient-of-variation-like precision measure. It is reported because precision is a
# reasonable thing to want, and it is NOT the default selection rule because on this
# project's own data it points the wrong way. Measured on the 2026-08/09 private-boat 2x2,
# four cells on one track:
#
#     cell            catch    PI width   rel width   trailer coverage_50 (nominal 0.500)
#     OFF x monthly  25,868      40,358      156.0%   0.538   calibrated
#     ON  x monthly  31,008      47,671      153.7%   0.523   calibrated
#     OFF x daily    37,359      54,769      146.6%   0.713   BROKEN
#     ON  x daily    42,344      60,294      142.4%   0.744   BROKEN
#
# On RELATIVE width the two miscalibrated cells look the most precise. On ABSOLUTE width the
# calibrated one wins, but only because catch and width are nearly proportional here (the
# ratio is 0.64 in every cell), so "narrowest absolute interval" is very close to "smallest
# estimate" and would systematically select the lowest harvest number. Neither is a model
# selection criterion on its own. Use precision to break a tie between rungs that are
# already adequate, never to decide adequacy.
# ---------------------------------------------------------------------------
bss_ar_rung_summary <- function(fit, pars = c("C_expected_sum", "C_sum", "E_sum")) {
  out <- list(catch_median = NA_real_, catch_lo95 = NA_real_, catch_hi95 = NA_real_,
              pi_width = NA_real_, pi_width_rel = NA_real_, effort_median = NA_real_)
  if (is.null(fit)) return(out)
  have <- pars[pars %in% fit@model_pars]
  if (!length(have)) return(out)
  s <- tryCatch(rstan::summary(fit, pars = have, probs = c(0.025, 0.5, 0.975))$summary,
                error = function(e) NULL)
  if (is.null(s)) return(out)
  g <- function(p, col) if (p %in% rownames(s) && col %in% colnames(s)) s[p, col] else NA_real_
  # C_expected_sum is the expected-catch total the convergence gate itself is computed on
  # (catch_par = "C_expected_sum"), so the rung's reported number and its gate verdict are
  # about the same quantity. C_sum, the predictive total, differs by the Poisson observation
  # noise and is not what the gate scores.
  cm <- g("C_expected_sum", "50%"); if (is.na(cm)) cm <- g("C_sum", "50%")
  lo <- g("C_expected_sum", "2.5%"); if (is.na(lo)) lo <- g("C_sum", "2.5%")
  hi <- g("C_expected_sum", "97.5%"); if (is.na(hi)) hi <- g("C_sum", "97.5%")
  out$catch_median <- cm; out$catch_lo95 <- lo; out$catch_hi95 <- hi
  out$pi_width     <- if (is.finite(hi) && is.finite(lo)) hi - lo else NA_real_
  out$pi_width_rel <- if (is.finite(out$pi_width) && is.finite(cm) && cm > 0) out$pi_width / cm else NA_real_
  out$effort_median <- g("E_sum", "50%")
  out
}
