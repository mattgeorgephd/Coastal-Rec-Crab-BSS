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
# desk_a42_render_replay_2026-11-25.R
#
# DESK, MINUTES, NO MCMC (CHANGE_REGISTER B89). Written AFTER the A42 render was read, to find why its
# clauses (3) and (4) failed; it judges nothing. The A42 render's own boat all-gear fit (from its fit
# cache) is replayed through the generated quantities with the switch on and off: on must equal the
# render's own C_expected_sum draw for draw, and the two medians give the switch's effect on THOSE draws.
# Why the render's draws are not B70's although the switch acts only in the generated quantities: Stan
# draws the generated quantities' random numbers (poisson_rng) from the sampler's own stream, and a
# Poisson draw consumes a rate-dependent number of them, so a changed rate shifts every later transition.
# HOW TO RUN  RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_a42_render_replay_2026-11-25.R
#             (needs the render's run folder and its fit cache, /home/user/bss_fit_cache_a42 on the container)
# WRITES      05_output/a42_render_replay_2026-11-25.csv
###############################################################################
suppressPackageStartupMessages(library(rstan))
N <- "05_output/20261125/pooled-CPUE-2024-25-190951"
sd_ <- readRDS(file.path(N, "ppc_draws_private_boat_all_gear_Dungeness_Kept.rds"))$stan_data; sd_ <- sd_[!startsWith(names(sd_), ".")]
cat("boat_day_carry in the render's data:", sd_$boat_day_carry, "\n")
f <- Sys.glob("/home/user/bss_fit_cache_a42/private_boat_all_gear_Dungeness_Kept_att1_*.rds"); stopifnot(length(f) == 1)
fit <- readRDS(f)
own <- as.matrix(fit, pars = "C_expected_sum")[, 1]
m <- rstan::stan_model("02_stan_models/crab_bss_pooled.stan", auto_write = FALSE)
dm <- as.matrix(fit, pars = names(rstan::get_inits(fit)[[1]]))
g <- function(flag) { d <- sd_; d$boat_day_carry <- as.integer(flag); as.matrix(rstan::gqs(m, data = d, draws = dm, seed = 1), pars = "C_expected_sum")[, 1] }
x1 <- g(1); x0 <- g(0)
cat(sprintf("render's own C_expected_sum median %.2f; replay with the switch on %.2f (max relative gap per draw %.2e); with it off %.2f; the switch's effect on these draws %.2f%%\n",
            median(own), median(x1), max(abs(x1 - own) / own), median(x0), 100 * (median(x1) / median(x0) - 1)))
write.csv(data.frame(quantity = c("render_own_median", "replay_on_median", "replay_on_max_rel_gap", "replay_off_median", "switch_effect_pct"),
                     value = c(median(own), median(x1), max(abs(x1 - own) / own), median(x0), 100 * (median(x1) / median(x0) - 1))),
          "05_output/a42_render_replay_2026-11-25.csv", row.names = FALSE)
