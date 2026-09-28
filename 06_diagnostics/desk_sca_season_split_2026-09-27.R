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
# desk_sca_season_split_2026-09-27.R
#
# A DESK SCREEN, SECONDS, NO MCMC: is the boat SCA effect the same size in the winter as in
# the rest of the season? (CHANGE_REGISTER D32; VALIDATION_CAMPAIGN Section 1y.3; B41.)
#
# WHY. The 2024-25 block cross-validation validated the boat SCA term on held-out weeks where
# the OSP counts are (March to September) and could not see the winter, where the term moves
# the boat estimate most. The model fits ONE coefficient for the season. Within 2024-25 the
# sampled winter days are too few and too small to say whether that is right (the module's
# own screen: winter 0.57 with p 0.49, rest 0.24 with p 0.0008; an interaction test on the
# same days, p about 0.2). The one place with power to speak to it now is the whole trailer
# record: every sampled launch day in effort_combined.xlsx, four seasons, with the archived
# advisory flag on each day.
#
# WHAT. The daily sum of the launch trailer counts (Westport + Ocean Shores Boat Launch, the
# BSS's boat effort series), one row per sampled day, against season + calendar month +
# weekend, and the SCA-or-higher flag (marine_hazard_flag_series(), the production
# definition) either as ONE term or SPLIT into the winter months (marine_hazard_winter_months,
# default December to February) and the rest, plus the interaction form (sca + sca:winter)
# whose coefficient IS the test. Two families, because they disagree on precision when the
# counts are this overdispersed: negative binomial (MASS::glm.nb; the family of the 2026-09-25
# offline screen) and quasi-Poisson (the module's screen). Per-season splits are reported too.
# This is a screen on sampled days with a coarse adjustment set (calendar month, not the creel
# day type), not a fit; it says which shape the term should be tested in, not what to adopt.
#
# WHAT IT SAID (2026-09-27, 794 sampled days, 315 advisory days, 129 of them winter days):
#   NB  constant 0.26 [0.21, 0.31]; split winter 0.43 [0.29, 0.62], rest 0.22 [0.18, 0.27];
#       interaction 1.96 [1.27, 3.05], z 3.0; AIC 4228.9 constant vs 4221.9 split.
#   QP  constant 0.23 [0.17, 0.32]; split winter 0.43 [0.17, 1.06], rest 0.22 [0.15, 0.30];
#       interaction 1.97 [0.74, 5.22], z 1.4.
#   Per season (QP): 2023-24 winter 0.33 / rest 0.24; 2024-25 0.59 / 0.17; 2025-26 0.45 / 0.23.
#   Read: the winter effect is real (a reduction in every season and in both families) and
#   about half the rest-of-season effect on the log scale (log 0.43 = -0.84 against
#   log 0.22 = -1.51); the same direction in all three seasons that have a winter; the
#   difference is 3 SE under the negative binomial and 1.4 SE under the quasi-Poisson. A
#   season-constant coefficient of 0.26 to 0.31 therefore over-states the winter reduction
#   and under-states the summer one. Rule 9 of run_marine_hazard_batch_2026-09-25.R says
#   what the 2024-25 fit (rung M6) can and cannot add to this.
#
# WRITES 05_output/marine_hazard_2026-09-27_season_split_screen.csv (merged by key).
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tibble); library(tidyr); library(lubridate) })
if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS is needed for the negative-binomial fits (install.packages('MASS')).")
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))   # B46: a file that fails to source stops the runner (it was hidden by try())
source("run_config.R")
source("03_R_functions/batch_verdict_helpers.R")   # merge_csv_by

e <- read_excel("04_input_files/effort_combined.xlsx", sheet = "data") |> mutate(date = as.Date(date))
areas <- run_config$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch")
b <- e |> filter(creel_area %in% areas, !is.na(boat_trailer_count)) |>
  group_by(date, season) |> summarise(y = sum(boat_trailer_count), n_counts = n(), .groups = "drop")
rng <- range(b$date)
p <- modifyList(run_config, list(est_date_start = rng[1], est_date_end = rng[2]))
fl <- marine_hazard_flag_series(p, date_start = rng[1], date_end = rng[2])
wm <- run_config$marine_hazard_winter_months %||% c(12, 1, 2)
d <- b |> left_join(fl |> dplyr::select(event_date, nws_sca_any, nws_sca_any_winter, nws_sca_any_rest), by = c("date" = "event_date")) |>
  filter(!is.na(nws_sca_any)) |>
  mutate(month = factor(format(date, "%m")), wk = as.integer(wday(date) %in% c(1, 7)), season = factor(season),
         winter = as.integer(as.integer(format(date, "%m")) %in% wm), sca = nws_sca_any)
cat(sprintf("Boat launch trailer counts: %d sampled days over %s to %s (%s); %d advisory days, %d winter days, %d winter advisory days; winter months %s\n",
            nrow(d), rng[1], rng[2], paste(levels(d$season), collapse = ", "), sum(d$sca), sum(d$winter), sum(d$sca * d$winter), paste(wm, collapse = "/")))

row <- function(model, family, seasons, term, m, n_days, n_adv) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co) || !is.finite(co[term, 2])) return(NULL)
  est <- co[term, 1]; se <- co[term, 2]
  data.frame(model = model, family = family, seasons = seasons, term = term, estimate = est, se = se, z = est / se,
             rate_ratio = exp(est), lo95 = exp(est - 1.96 * se), hi95 = exp(est + 1.96 * se),
             n_days = n_days, n_advisory_days = n_adv, aic = if (family == "negbin") AIC(m) else NA_real_, stringsAsFactors = FALSE)
}
R <- list()
for (fam in c("negbin", "quasipoisson")) {
  fit <- function(f, dd) if (fam == "negbin") MASS::glm.nb(f, data = dd) else glm(f, data = dd, family = quasipoisson)
  m0 <- fit(y ~ season + month + wk + sca, d)
  m1 <- fit(y ~ season + month + wk + nws_sca_any_winter + nws_sca_any_rest, d)
  m2 <- fit(y ~ season + month + wk + sca + sca:winter, d)
  R <- c(R, list(row("constant", fam, "all", "sca", m0, nrow(d), sum(d$sca)),
                 row("split", fam, "all", "nws_sca_any_winter", m1, nrow(d), sum(d$nws_sca_any_winter)),
                 row("split", fam, "all", "nws_sca_any_rest", m1, nrow(d), sum(d$nws_sca_any_rest)),
                 row("interaction", fam, "all", "sca", m2, nrow(d), sum(d$sca)),
                 row("interaction", fam, "all", "sca:winter", m2, nrow(d), sum(d$sca * d$winter))))
  for (s in levels(d$season)) {
    ds <- d[d$season == s, ]
    if (sum(ds$nws_sca_any_winter) < 5 || sum(ds$nws_sca_any_rest) < 5) next
    m <- tryCatch(suppressWarnings(fit(y ~ month + wk + nws_sca_any_winter + nws_sca_any_rest, ds)), error = function(e) NULL)
    if (is.null(m)) next
    R <- c(R, list(row("split", fam, s, "nws_sca_any_winter", m, nrow(ds), sum(ds$nws_sca_any_winter)),
                   row("split", fam, s, "nws_sca_any_rest", m, nrow(ds), sum(ds$nws_sca_any_rest))))
  }
}
out <- do.call(rbind, Filter(Negate(is.null), R))
out$written <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
out$winter_months <- paste(wm, collapse = "/")
print(out[, c("model", "family", "seasons", "term", "rate_ratio", "lo95", "hi95", "z", "n_days", "n_advisory_days")], digits = 3, row.names = FALSE)
merge_csv_by(out, file.path("05_output", "marine_hazard_2026-09-27_season_split_screen.csv"), key = c("model", "family", "seasons", "term"))
cat("\nwrote 05_output/marine_hazard_2026-09-27_season_split_screen.csv\n")
cat("READ AS A SCREEN: sampled days, calendar-month adjustment, no creel day type, daily sums (not the per-count NB2 the BSS fits).\n")
cat("It says which SHAPE the term should be tested in (rule 9 of run_marine_hazard_batch_2026-09-25.R), not what to adopt.\n")
