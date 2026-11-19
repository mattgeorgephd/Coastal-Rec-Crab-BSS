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
# desk_advisory_fleet_composition_2026-11-19.R
#
# DESK, A MINUTE, NO MCMC (CHANGE_REGISTER B82, for D32 and D30; a rider on B80).
#
# THE QUESTION. D32 says the advisory effect on the launch trailer counts is weaker in
# December to February, and the split arm fits a separate coefficient for those months. A
# calendar boundary is a description, not a reason. A reason is on the table: in the winter
# nearly every boat at the launch is crabbing (the fitted crabbing share is 0.94 in December
# 2024 and 0.26 in August 2025), and B81 found the crab-only share of OSP's sampled boats at
# 29% on advisory days against 12% on other days. If crabbing boats are deterred less than
# the rest of the fleet, the effect on ALL boats must weaken as the crabbing share rises,
# month by month and not at a boundary, and two numbers (the crabbers' response and the
# other boats') would describe every month.
#
# THE DATA. The transfer script's series: the summed launch trailer count per sampled day
# with the archived advisory flag, here for the two in-scope seasons inside their all-gear
# windows (December to mid-September), 2024-25 and 2025-26. The crabbing share f of each
# season-month is the posterior median of the season's committed run (B70 and B76,
# crab_fraction_strata_private_boat_all_gear_*.csv); it is itself fitted from the samplers'
# boat contacts and, on 2024-25, OSP's crab-only counts, never from the trailer counts.
#
# THE MODELS, negative binomial with a level per season-month and a weekend term:
#   C  + advisory                              one coefficient (the shipped form)
#   W  + advisory + advisory x winter          the split (December to February against the rest)
#   F  + advisory + advisory x (f - mean f)    the effect as a log-linear function of the share
#   M  the advisory factor is  f exp(b_crab) + (1 - f) exp(b_other):  the crabbers' response
#      and the other boats', weighted by the month's share (maximum likelihood)
#
# THE RULE, WRITTEN BEFORE ANY OF THESE WAS FITTED. R10, exploratory: it judges nothing in
# B80. Known when it was written: the winter and rest coefficients of the desk regression on
# each season (so W improving on C is not news), the monthly shares, and B81's result.
#   (1) DOSE-RESPONSE. F's interaction coefficient with its 95% interval. Excludes 0 and
#       positive: the advisory effect weakens as the crabbing share rises. Covers 0: not
#       detectable on two seasons.
#   (2) BOUNDARY OR SHARE. W against F, which carry the same number of coefficients: the
#       difference in AIC, and the difference in the held-out log score with one ISO week
#       left out at a time, with its paired standard error over weeks. Beyond 2 SE either
#       way names the better description; inside, the two cannot be told apart.
#   (3) THE MIXTURE AGAINST AN INDEPENDENT STREAM. M's exp(b_crab) and exp(b_other) with
#       95% intervals, and its implied odds ratio of the crabbing share on an advisory day,
#       exp(b_crab - b_other), against the crab-only odds ratio OSP's boats show (B81's Ma:
#       2.82 [1.55, 5.15] on the schedule's sampled count, 3.58 [1.53, 8.37] if every boat
#       was sampled). AGREEMENT needs both: the mixture's interval excludes 1 (the trailer
#       counts alone separate the two responses) and it overlaps OSP's. An interval covering
#       1 says the trailer counts cannot separate them, whatever it overlaps: in the dry run
#       with the flag permuted (FLEET_COMPOSITION_DRY=1, the only run made before this rule
#       was committed) the interval was 0.43 to 7.4 and overlapped OSP's, which is why
#       overlap alone is not read. M's held-out score against C, W and F is reported beside it.
#   STATED IN ADVANCE. The share varies mostly BETWEEN winter and the rest, so (1) and (2)
#   lean on a few spring months (March to May) where the share is in between; two seasons may
#   not separate W from F. OSP's crab-only is a lower bound on crabbing boats and its days are
#   March to September. Everything here is a desk regression on sampled days with a plug-in
#   score, not the state-space model. Nothing changes on the result; a model form that used
#   the share would be new Stan code and a decision for Matt.
#   EXPECTED: (1) positive, its interval possibly covering 0; (2) inside 2 SE; (3) the
#   mixture's odds ratio at or above OSP's.
#
# WRITES 05_output/advisory_replication_2026-11-19_fleet_composition.csv (the readings) and
#        05_output/advisory_replication_2026-11-19_fleet_composition_months.csv (the advisory
#        rate ratio month by month, for the figure).
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tibble); library(tidyr); library(lubridate) })
if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS is needed for the negative-binomial fits (install.packages('MASS')).")
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
source("run_config.R")
OUT <- function(kind) file.path("05_output", sprintf("advisory_replication_2026-11-19_%s.csv", kind))

# ---- the series (the transfer script's), the all-gear windows, and the month's crabbing share ----
WINDOWS <- list(`2024-25` = as.Date(c("2024-12-01", "2025-09-15")), `2025-26` = as.Date(c("2025-12-01", "2026-09-15")))
FRUN <- c(`2024-25` = file.path("05_output", "20261019", "pooled-CPUE-2024-25"), `2025-26` = file.path("05_output", "20261020", "pooled-CPUE-season-2025-26"))
e <- read_excel("04_input_files/effort_combined.xlsx", sheet = "data") |> mutate(date = as.Date(date))
areas <- run_config$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch")
b <- e |> filter(creel_area %in% areas, !is.na(boat_trailer_count)) |> group_by(date) |> summarise(y = sum(boat_trailer_count), .groups = "drop")
rng <- range(b$date)
fl <- marine_hazard_flag_series(modifyList(run_config, list(est_date_start = rng[1], est_date_end = rng[2])), date_start = rng[1], date_end = rng[2])
wm <- as.integer(run_config$marine_hazard_winter_months %||% c(12, 1, 2))
fsh <- bind_rows(lapply(names(FRUN), function(s) utils::read.csv(file.path(FRUN[[s]], "crab_fraction_strata_private_boat_all_gear_Dungeness_Kept.csv"), stringsAsFactors = FALSE) |>
                          transmute(season = s, ym = label, f = f_median, f_lo = f_lo95, f_hi = f_hi95)))
d <- b |> left_join(fl |> dplyr::select(event_date, nws_sca_any), by = c("date" = "event_date")) |> filter(!is.na(nws_sca_any)) |>
  mutate(season = case_when(date >= WINDOWS[[1]][1] & date <= WINDOWS[[1]][2] ~ names(WINDOWS)[1], date >= WINDOWS[[2]][1] & date <= WINDOWS[[2]][2] ~ names(WINDOWS)[2], TRUE ~ NA_character_)) |>
  filter(!is.na(season)) |>
  mutate(ym = format(date, "%Y-%m"), sm = factor(paste(season, ym)), wk = as.integer(wday(date) %in% c(1, 7)), sca = as.numeric(nws_sca_any),
         winter = as.integer(as.integer(format(date, "%m")) %in% wm), week = paste(season, format(date, "%G-W%V"))) |>
  left_join(fsh, by = c("season", "ym"))
stopifnot(!anyNA(d$f))
# FLEET_COMPOSITION_DRY=1 checks the mechanics without reading the answer: the advisory flag is permuted
# within each season-month (a null data set) and nothing is written.
DRY <- identical(Sys.getenv("FLEET_COMPOSITION_DRY"), "1")
if (DRY) { set.seed(1); d <- d |> group_by(sm) |> mutate(sca = sample(sca)) |> ungroup(); cat("DRY RUN: the advisory flag is permuted within each season-month; nothing is written.\n") }
d$fc <- d$f - mean(d$f)
cat(sprintf("Sampled launch days in the two all-gear windows: %d (%d advisory); season-months %d; weeks %d; crabbing share %.2f to %.2f\n",
            nrow(d), sum(d$sca), nlevels(d$sm), length(unique(d$week)), min(d$f), max(d$f)))

# ---- the four models -----------------------------------------------------------------------------
nb <- function(f, dd) suppressWarnings(MASS::glm.nb(f, data = dd, control = stats::glm.control(maxit = 100)))
F_C <- y ~ sm + wk + sca; F_W <- y ~ sm + wk + sca + I(sca * winter); F_F <- y ~ sm + wk + sca + I(sca * fc)
fit_mix <- function(dd, start = NULL) {                 # M: log mu = X beta + sca * log(f e^bc + (1 - f) e^bo)
  X <- stats::model.matrix(~ sm + wk, dd)
  X <- X[, colSums(abs(X)) > 0, drop = FALSE]           # a held-out week can empty a season-month's column
  if (is.null(start)) { m0 <- nb(y ~ sm + wk + sca, dd); cf <- stats::coef(m0); start <- c(cf[colnames(X)], bc = unname(cf["sca"]), bo = unname(cf["sca"]) - 0.5, lth = log(m0$theta)) }
  start[is.na(start)] <- 0
  nll <- function(p) { k <- ncol(X); eta <- as.numeric(X %*% p[seq_len(k)]) + dd$sca * log(dd$f * exp(p[k + 1]) + (1 - dd$f) * exp(p[k + 2]))
    -sum(stats::dnbinom(dd$y, mu = exp(eta), size = exp(p[k + 3]), log = TRUE)) }
  o <- stats::optim(start, nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  list(par = o$par, nll = o$value, k = ncol(X), cols = colnames(X), conv = o$convergence, nllf = nll)
}
pred_mix <- function(m, new) { X <- stats::model.matrix(~ sm + wk, new)[, m$cols, drop = FALSE]; k <- m$k
  eta <- as.numeric(X %*% m$par[seq_len(k)]) + new$sca * log(new$f * exp(m$par[k + 1]) + (1 - new$f) * exp(m$par[k + 2]))
  sum(stats::dnbinom(new$y, mu = exp(eta), size = exp(m$par[k + 3]), log = TRUE)) }
mC <- nb(F_C, d); mW <- nb(F_W, d); mF <- nb(F_F, d); mM <- fit_mix(d)
H <- stats::optimHess(mM$par, mM$nllf); V <- tryCatch(solve(H), error = function(e) matrix(NA_real_, nrow(H), ncol(H)))
k <- mM$k; bc <- mM$par[k + 1]; bo <- mM$par[k + 2]; se_bc <- sqrt(V[k + 1, k + 1]); se_bo <- sqrt(V[k + 2, k + 2])
se_d <- sqrt(V[k + 1, k + 1] + V[k + 2, k + 2] - 2 * V[k + 1, k + 2])
aic <- c(C = stats::AIC(mC), W = stats::AIC(mW), F = stats::AIC(mF), M = 2 * mM$nll + 2 * (k + 3))

# ---- one ISO week left out at a time ----------------------------------------------------------------
lpd <- function(m, new) { mu <- exp(as.numeric(stats::predict(m, newdata = new, type = "link"))); sum(stats::dnbinom(new$y, mu = mu, size = m$theta, log = TRUE)) }
weeks <- unique(d$week); sc <- list()
for (w in weeks) { te <- d[d$week == w, ]; tr <- d[d$week != w, ]
  if (!all(te$sm %in% tr$sm)) { sc[[w]] <- tibble(week = w, scored = FALSE); next }      # a season-month with no other week cannot be predicted
  tr$sm <- droplevels(tr$sm); te$sm <- factor(te$sm, levels = levels(tr$sm))
  s <- tryCatch({ mm <- fit_mix(tr)
    c(C = lpd(nb(F_C, tr), te), W = lpd(nb(F_W, tr), te), F = lpd(nb(F_F, tr), te), M = pred_mix(mm, te)) }, error = function(e) rep(NA_real_, 4))
  sc[[w]] <- tibble(week = w, scored = all(is.finite(s)), C = s[1], W = s[2], F = s[3], M = s[4], n_days = nrow(te), n_advisory = sum(te$sca), winter_week = as.integer(mean(te$winter) >= 0.5)) }
sc <- bind_rows(sc); ok <- sc[sc$scored %in% TRUE, ]
pair <- function(a, bq) { x <- ok[[a]] - ok[[bq]]; n <- length(x); s <- sqrt(n) * stats::sd(x); c(diff = sum(x), se = s, ratio = sum(x) / s, n = n, pos = sum(x > 0)) }

# ---- the readings ----------------------------------------------------------------------------------
R <- list(); add <- function(reading, quantity, estimate, lo = NA_real_, hi = NA_real_, se = NA_real_, note = "")
  R[[length(R) + 1]] <<- data.frame(reading = reading, quantity = quantity, estimate = estimate, lo95 = lo, hi95 = hi, se = se, note = note, stringsAsFactors = FALSE)
cc <- summary(mC)$coefficients; cw <- summary(mW)$coefficients; cf <- summary(mF)$coefficients
ci <- function(tab, nm) c(tab[nm, 1], tab[nm, 1] - 1.96 * tab[nm, 2], tab[nm, 1] + 1.96 * tab[nm, 2], tab[nm, 2])
x <- ci(cc, "sca"); add("C", "advisory coefficient, one for the season", x[1], x[2], x[3], x[4], sprintf("%d days, %d advisory", nrow(d), sum(d$sca)))
x <- ci(cw, "sca"); add("W", "advisory coefficient, rest of the season", x[1], x[2], x[3], x[4])
x <- ci(cw, "I(sca * winter)"); add("W", "winter minus rest", x[1], x[2], x[3], x[4])
x <- ci(cf, "I(sca * fc)"); v1 <- if (x[2] > 0) "the advisory effect WEAKENS as the crabbing share rises" else if (x[3] < 0) "the advisory effect STRENGTHENS as the crabbing share rises" else "not detectable on two seasons"
add("(1) dose-response", "F: change in the advisory coefficient per unit of crabbing share", x[1], x[2], x[3], x[4], v1)
add("(1) dose-response", "F: implied rate ratio at a share of 0.25", exp(cf["sca", 1] + cf["I(sca * fc)", 1] * (0.25 - mean(d$f))))
add("(1) dose-response", "F: implied rate ratio at a share of 0.90", exp(cf["sca", 1] + cf["I(sca * fc)", 1] * (0.90 - mean(d$f))))
add("(2) boundary or share", "AIC, W minus F (negative favours the boundary)", aic[["W"]] - aic[["F"]], note = sprintf("AIC C %.1f, W %.1f, F %.1f, M %.1f", aic[["C"]], aic[["W"]], aic[["F"]], aic[["M"]]))
p <- pair("F", "W"); v2 <- if (p[["ratio"]] > 2) "the share describes it better" else if (p[["ratio"]] < -2) "the boundary describes it better" else "the two cannot be told apart"
add("(2) boundary or share", "held-out log score, F minus W", p[["diff"]], se = p[["se"]], note = sprintf("%+.2f SE over %d weeks, %d positive: %s", p[["ratio"]], p[["n"]], p[["pos"]], v2))
for (q in list(c("W", "C"), c("F", "C"), c("M", "C"), c("M", "W"), c("M", "F"))) { p <- pair(q[1], q[2])
  add("held-out log score", sprintf("%s minus %s", q[1], q[2]), p[["diff"]], se = p[["se"]], note = sprintf("%+.2f SE over %d weeks, %d positive", p[["ratio"]], p[["n"]], p[["pos"]])) }
add("(3) mixture", "exp(b_crab): crabbing boats on an advisory day, as a share of an ordinary day's", exp(bc), exp(bc - 1.96 * se_bc), exp(bc + 1.96 * se_bc), se_bc, sprintf("optim convergence %d", mM$conv))
add("(3) mixture", "exp(b_other): other boats on an advisory day, as a share of an ordinary day's", exp(bo), exp(bo - 1.96 * se_bo), exp(bo + 1.96 * se_bo), se_bo)
orr <- c(exp(bc - bo), exp(bc - bo - 1.96 * se_d), exp(bc - bo + 1.96 * se_d))
osp <- tryCatch(utils::read.csv(OUT("crab_share"), stringsAsFactors = FALSE), error = function(e) NULL)
oa <- if (is.null(osp)) NULL else osp[grepl("^Ma", osp$model), ]
ov <- if (is.null(oa)) NA else any(orr[2] <= oa$or_hi95 & orr[3] >= oa$or_lo95)
add("(3) mixture", "implied odds ratio of the crabbing share on an advisory day, exp(b_crab - b_other)", orr[1], orr[2], orr[3], se_d,
    if (is.null(oa)) "B81's table is not there" else sprintf("OSP crab-only odds ratio %.2f [%.2f, %.2f] (schedule n), %.2f [%.2f, %.2f] (every boat): %s",
       oa$odds_ratio[1], oa$or_lo95[1], oa$or_hi95[1], oa$odds_ratio[2], oa$or_lo95[2], oa$or_hi95[2],
       if (orr[2] <= 1 && orr[3] >= 1) "the mixture's interval covers 1: the trailer counts alone cannot separate the two responses" else
       if (isTRUE(ov)) "AGREEMENT: the interval excludes 1 and overlaps OSP's" else "the interval excludes 1 and does NOT overlap OSP's"))
res <- do.call(rbind, R)
op <- options(width = 230); print(res, digits = 3, row.names = FALSE, right = FALSE); options(op)
cat(sprintf("\nWeeks scored: %d of %d\n", nrow(ok), length(weeks)))
if (!DRY) utils::write.csv(res, OUT("fleet_composition"), row.names = FALSE)

# ---- the advisory rate ratio month by month, for the figure ---------------------------------------
mon <- d |> group_by(season, ym, f, f_lo, f_hi) |> summarise(days = n(), advisory_days = sum(sca), other_days = sum(1 - sca), mean_advisory = mean(y[sca == 1]), mean_other = mean(y[sca == 0]), .groups = "drop")
mon$log_rr <- NA_real_; mon$se_log_rr <- NA_real_
for (i in seq_len(nrow(mon))) { if (mon$advisory_days[i] < 3 || mon$other_days[i] < 3) next
  dd <- d[d$season == mon$season[i] & d$ym == mon$ym[i], ]
  m <- tryCatch(nb(y ~ wk + sca, dd), error = function(e) NULL); if (is.null(m)) next
  s <- summary(m)$coefficients; if ("sca" %in% rownames(s)) { mon$log_rr[i] <- s["sca", 1]; mon$se_log_rr[i] <- s["sca", 2] } }
mon$fit_F <- cf["sca", 1] + cf["I(sca * fc)", 1] * (mon$f - mean(d$f)); mon$fit_M <- log(mon$f * exp(bc) + (1 - mon$f) * exp(bo))
mon$fit_W <- cw["sca", 1] + cw["I(sca * winter)", 1] * as.integer(as.integer(substr(mon$ym, 6, 7)) %in% wm); mon$fit_C <- cc["sca", 1]
if (!DRY) { utils::write.csv(mon, OUT("fleet_composition_months"), row.names = FALSE)
  cat("Wrote", OUT("fleet_composition"), "and", OUT("fleet_composition_months"), "\n") }
