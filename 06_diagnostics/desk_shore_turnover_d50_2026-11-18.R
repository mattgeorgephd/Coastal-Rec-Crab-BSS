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
# desk_shore_turnover_d50_2026-11-18.R
#
# A DESK COMPARISON, SECONDS, NO MCMC: which of the two shore turnovers the I/E days give
# should be trusted? (CHANGE_REGISTER D50 and B78; VALIDATION_CAMPAIGN Section 1z.20.)
#
# WHY. On 2025-26 the same WDF20 ingress/egress days give two turnovers. The derivation
# (estimate_shore_turnover(): arrivals over presence at the count hours, both read off the
# I/E tally) says 2.595 on the 17 in-window days and 2.516 pooled over 40. The shore all-gear
# fit's tau_bar is 2.296 [1.98, 2.69] on 16 I/E observations, and the likelihood alone puts it
# near 2.1. The PE applies the first, the BSS the second, and the shore all-gear BSS sits
# 14.5% below the PE on effort (D50).
#
# WHAT THE CODE SAYS, before any number: the I/E observation is the WDF20 crabber ARRIVAL
# count (fetch_ie_data(), ie_trips, one location), and Stan compares it with lambda_E * L,
# where lambda_E is calibrated by Gear_I ~ NB2(lambda_E * R_G * day effect) to the SHORE gear
# count, Float 20 PLUS Float 17-21 (shore_dock_counts(); the unsampled Float 17-21 counts
# ratio-filled since B45). Nothing relates the one float to the two, and the I/E prediction
# carries no per-day effort effect although the count likelihood has since A39.
#
# THE IDENTITY. Per I/E day d the two readings differ by exactly lambda_d / Pbar_d (the
# fit's latent crabbers present over the tally's presence at the season's count hours), and
#     lambda_d / Pbar_d = F_smooth * F_cover * F_gear * F_time
#   F_smooth = lambda_d / ((G20 + G17) / R_G)   the fitted level over the day's own counts (H2)
#   F_cover  = (G20 + G17) / G20                two floats' gear over Float 20's (H1)
#   F_gear   = (G20 / R_G) / P_c                Float 20 gear as crabbers over the tally's
#                                               presence at the day's own count instants (H4)
#   F_time   = P_c / Pbar_d                     presence at the day's count instants over
#                                               presence at the season's count-hour weights (H3)
# G20, G17 are the day's mean counts (sequences <= bss_max_count_seq, as the prep takes them),
# R_G and lambda_d the fit's posterior medians. A day needs a same-day Float 20 count for the
# three count factors. The product F_gear * F_time does not depend on how P_c is measured;
# only the split between the two does.
#
# TWO MEASURES OF P_c, AND WHY. (A) AS FIRST COMPUTED, and the one the rule was first read
# on: the mean of the workbook's crabber_flow column over the count's hour bin, the binning
# the derivation uses. (B) THE REPAIR, made after that first table showed what (A) rests on:
# the workbook writes crabber_flow = 0 on each day's FIRST and LAST interval row whatever the
# running balance is (2026-04-03: 5 crabbers on at 07:45 and a flow of 0; 21 still present at
# 18:30 and a flow of 0), and on an I/E day the sampler takes the first gear count as the
# tally starts and the last as it ends, so those are exactly the rows the counts land on.
# (B) is the tally's own running balance, cumsum(on - off), interpolated to the count time.
# The derivation itself is barely touched by the zeros: the bins they fall in carry about 2%
# of the season's count-hour weight.
#
# THE READING RULES were written into CHANGE_REGISTER B78 (commit 8af6e2f) before this script
# produced a number: with g the mean log of lambda_d / Pbar_d over the fit's I/E days, a factor
# is MATERIAL at |f| >= 0.03 and ACCOUNTS FOR the gap at |f| >= |g| / 2 with g's sign. The
# script applies them mechanically, on (A) and on (B), and prints both.
#
# BESIDE THE RULE, labelled as such: the same factors as RATIOS OF SUMS (the scale the
# derivation is built on; a mean log weights a day with two crabbers present like a day with
# eighty), the gear per present crabber count by count for EVERY I/E day in the workbook, and
# a desk approximation of what the fit's tau_bar would be with the Float 20 share in the I/E
# prediction (prior and likelihood combined by precision on the log scale; it reproduces the
# fitted tau_bar, and it is an approximation, not a fit).
#
# WRITES (merged by key) to 05_output/, each named shore_turnover_d50_2026-11-18_<kind>.csv:
#   by_day     one row per in-window I/E day and fit
#   factors    the factor means per fit, the rule's verdicts, the desk approximation
#   by_count   every Float 20 gear count taken on a WDF20 I/E day, all seasons
#   by_season  gear per present crabber and the two-float ratio by season
# Fits nothing; reads the committed inputs through the pipeline's own readers and the
# committed run folders named in RUNS below.
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
suppressPackageStartupMessages({ library(here); library(readxl); library(dplyr); library(tibble); library(tidyr)
                                 library(purrr); library(stringr); library(lubridate) })
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
source("run_config.R")

# The two windows and the committed run each is read against. 2024-25 is run_config.R as
# shipped; 2025-26 is the block of NEW_SEASON_GUIDE.md section 7.1.
season_window <- function(y) list(
  est_date_start = sprintf("%d-09-16", y), est_date_end = sprintf("%d-09-15", y + 1), season_filter = sprintf("%d-%02d", y, (y + 1) %% 100),
  pot_closures = NULL, census_windows = NULL,
  pot_closure_start = sprintf("%d-09-16", y), pot_closure_end = sprintf("%d-11-30", y), pot_open_date = sprintf("%d-12-01", y),
  census_start_date = sprintf("%d-12-01", y), census_end_date = sprintf("%d-01-03", y + 1), commercial_opener = sprintf("%d-01-04", y + 1),
  run_tag = sprintf("season-%d-%02d", y, (y + 1) %% 100))
RUNS <- list(
  `2024-25` = list(window = list(), run = "05_output/20261019/pooled-CPUE-2024-25", label = "B70"),
  `2025-26` = list(window = season_window(2025), run = "05_output/20261020/pooled-CPUE-season-2025-26", label = "B76"))
RULE_MATERIAL <- 0.03                                    # CHANGE_REGISTER B78, written before the table
quietly <- function(expr) { r <- NULL; utils::capture.output(r <- expr); r }
rd  <- function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL
med <- function(sm, name) { r <- sm[sm$parameter == name, "50%"]; if (length(r)) r[1] else NA_real_ }
gm  <- function(x) { x <- x[is.finite(x) & x > 0]; if (length(x)) exp(mean(log(x))) else NA_real_ }
mlog <- function(x) { x <- x[is.finite(x) & x > 0]; if (length(x)) mean(log(x)) else NA_real_ }
selog <- function(x) { x <- x[is.finite(x) & x > 0]; if (length(x) > 1) stats::sd(log(x)) / sqrt(length(x)) else NA_real_ }
with_window <- function(w) { p <- run_config; for (k in names(w)) p[k] <- list(w[[k]])
  p$boat_require_gear_time <- TRUE                       # the pooled driver's params_model value
  p$crabbing_holiday_dates <- quietly(read_crabbing_holidays(p)); p }

# ---- 1. THE TALLY: every WDF20 interval row, with its own running balance -----------------
tally <- read_input_workbook(here("04_input_files", run_config$ie_data_file), sheet = run_config$ie_sheet) |>
  mutate(event_date = as.Date(date)) |> filter(location_name == run_config$ie_shore_location)
tally$hour <- .ie_hour_of(tally$time)
tally <- tally |> arrange(event_date, hour) |> group_by(event_date) |>
  mutate(on = replace_na(as.numeric(crabbers_on), 0), off = replace_na(as.numeric(crabbers_off), 0),
         flow_wb = replace_na(as.numeric(crabber_flow), 0), balance = cumsum(on - off)) |> ungroup()
tally_days <- tally |> group_by(event_date) |>
  summarise(season = first(season), survey_start = min(hour), survey_end = max(hour),
            arrivals = sum(on), departures = sum(off),
            on_first_interval = first(on),                # the headcount as the tally starts
            balance_at_end = last(balance),               # still present as it ends (< 0: more left than were seen arriving)
            first_row_zeroed = first(flow_wb) == 0 && first(balance) != 0,
            last_row_zeroed  = last(flow_wb) == 0 && last(balance) != 0,
            arrivals_after_18 = sum(on[hour >= 18]), .groups = "drop")
# presence at a count instant: (A) the workbook's flow, hour-bin mean; (B) the running balance, interpolated
presence_at <- function(d, h) {
  x <- tally[tally$event_date == d, , drop = FALSE]
  if (!nrow(x) || !is.finite(h)) return(c(A = NA_real_, B = NA_real_, pos = NA_real_))
  a <- x$flow_wb[floor(x$hour) == floor(h)]
  b <- if (h < min(x$hour) - 0.2501 || h > max(x$hour) + 0.2501) NA_real_
       else stats::approx(x$hour, x$balance, xout = min(max(h, min(x$hour)), max(x$hour)), ties = "ordered")$y
  pos <- if (!is.finite(b)) NA_real_ else if (h <= min(x$hour) + 0.2501) 1 else if (h >= max(x$hour) - 0.2501) 3 else 2
  c(A = if (length(a)) mean(a) else NA_real_, B = b, pos = pos)
}
add_presence <- function(cnt) {
  m <- t(mapply(presence_at, cnt$event_date, cnt$count_hour))
  cnt$P_A <- m[, "A"]; cnt$P_B <- m[, "B"]
  cnt$position <- c("tally start", "interior", "tally end")[m[, "pos"]]
  cnt
}

# ---- 2. THE TWO WINDOWS against their fits -----------------------------------------------
one_window <- function(season, spec) {
  p    <- with_window(spec$window)
  subs <- quietly(build_subseasons(p))
  dwg  <- quietly(fetch_crab_data(p))
  ie   <- quietly(fetch_ie_data(p))
  st   <- quietly(estimate_shore_turnover(attr(ie, "ie_intervals"), dwg$shore_effort, p))
  pr   <- quietly(bss_resolve_tau_shore_prior(p, st, quiet = TRUE))
  # the same derivation with the workbook's flow column replaced by the tally's running balance
  iv_b <- attr(ie, "ie_intervals") |> left_join(tally |> select(event_date, hour, balance), by = c("event_date", "hour")) |>
    mutate(crabber_flow = ifelse(is.finite(balance), balance, crabber_flow))
  st_b <- quietly(estimate_shore_turnover(iv_b, dwg$shore_effort, p))
  sched <- dwg$shore_effort |> filter(count_sequence <= p$bss_max_count_seq) |>
    mutate(ie_day = event_date %in% tally_days$event_date, off = count_hour < 10 | count_hour >= 15) |>
    group_by(ie_day) |> summarise(n = n(), off = sum(off), .groups = "drop")
  cnt  <- dwg$shore_effort |> filter(count_sequence <= p$bss_max_count_seq, event_date %in% st$by_day$event_date) |> add_presence()
  day_counts <- cnt |> group_by(event_date) |>
    summarise(n_counts = n(), n_counts_in_survey = sum(is.finite(P_B)),
              G20 = mean(f20_gear), G17 = mean(f17_gear), Gtot = mean(count_quantity),
              f17_observed = sum(f17_source == "observed"),
              P_c_A = if (any(is.finite(P_A))) mean(P_A, na.rm = TRUE) else NA_real_,
              P_c_B = if (any(is.finite(P_B))) mean(P_B, na.rm = TRUE) else NA_real_,
              counts_off_schedule = sum(count_hour < 10 | count_hour >= 15),
              count_hours = paste(sprintf("%.2f", count_hour), collapse = " "), .groups = "drop")
  out <- list()
  for (ss in subs) {
    key <- ss$name
    fs  <- rd(file.path(spec$run, sprintf("bss_full_summary_shore_%s_Dungeness_Kept.csv", key)))
    sp  <- rd(file.path(spec$run, sprintf("structural_params_shore_%s_Dungeness_Kept.csv", key)))
    fds <- rd(file.path(spec$run, "fit_data_summary.csv"))
    if (is.null(fs) || is.null(sp)) next
    days <- tibble(event_date = seq(as.Date(ss$start), as.Date(ss$end), by = "day")) |> mutate(day_index = row_number())
    spv <- function(nm) { v <- sp$median[sp$parameter == nm]; if (length(v)) v[1] else NA_real_ }
    sigma_uE <- med(fs, "sigma_uE[1]")                    # NA when the fit runs without the per-day effort effect (A40)
    # the per-day deviation the COUNT likelihood carries (A39): exp(sigma_uE u - sigma_uE^2 / 2), indexed by count day
    cdays <- sort(unique(dwg$shore_effort$event_date[dwg$shore_effort$count_sequence <= p$bss_max_count_seq &
                                                     dwg$shore_effort$event_date >= as.Date(ss$start) & dwg$shore_effort$event_date <= as.Date(ss$end)]))
    dev <- tibble(event_date = cdays,
                  day_effect = if (is.finite(sigma_uE)) exp(sigma_uE * vapply(seq_along(cdays), function(k) med(fs, sprintf("uE_raw[%d]", k)), numeric(1)) - 0.5 * sigma_uE^2)
                               else NA_real_)
    tab <- st$by_day |> mutate(event_date = as.Date(event_date)) |> inner_join(days, by = "event_date") |>
      mutate(lambda = vapply(day_index, function(d) med(fs, sprintf("lambda_E_S[1,%d,1]", d)), numeric(1)),
             L_day  = vapply(day_index, function(d) med(fs, sprintf("L[%d]", d)), numeric(1))) |>
      left_join(day_counts, by = "event_date") |> left_join(dev, by = "event_date") |>
      left_join(tally_days |> select(event_date, survey_start, survey_end, on_first_interval, balance_at_end), by = "event_date") |>
      mutate(season_window = season, run = spec$label, fit = sprintf("shore_%s", key),
             R_G = spv("R_G"), tau_bar = spv("tau_bar[1]"), sigma_IE = spv("sigma_IE"), sigma_uE = sigma_uE,
             n_ie_obs_fit = fds$n_ie_obs[fds$fit == sprintf("shore_%s_Dungeness_Kept", key)][1],
             n_count_days_fit = length(cdays), mean_log_day_effect_all_count_days = mlog(dev$day_effect),
             prior_mu = pr$tau_shore_prior_mu, prior_sd = pr$tau_shore_prior_sigma, shared_tau_sigma = p$shared_tau_sigma,
             log_sd_days_pooled = st$log_sd_days,   # between-day log-SD of the daily derivation ratios, all I/E days
             Pbar = presence_at_counts,
             tau_derivation = arrivals / Pbar,           # the derivation's reading on the day
             tau_likelihood = arrivals / lambda,         # what the I/E likelihood sees on the day
             F_total  = lambda / Pbar,
             F_smooth = lambda / (Gtot / R_G),
             F_cover  = Gtot / G20,
             F_gear_A = (G20 / R_G) / P_c_A, F_time_A = P_c_A / Pbar,
             F_gear_B = (G20 / R_G) / P_c_B, F_time_B = P_c_B / Pbar)
    out[[key]] <- tab
  }
  list(by_day = bind_rows(out), st = st, st_b = st_b, sched = sched, f17 = attr(dwg$shore_effort, "f17_fill"))
}
res <- lapply(names(RUNS), function(s) one_window(s, RUNS[[s]])); names(res) <- names(RUNS)
by_day <- bind_rows(lapply(res, `[[`, "by_day")) |>
  select(season_window, run, fit, event_date, day_type, arrivals, peak, Pbar, lambda, L_day, day_effect,
         n_counts, n_counts_in_survey, counts_off_schedule, f17_observed, count_hours, G20, G17, Gtot, P_c_A, P_c_B,
         R_G, tau_bar, sigma_IE, sigma_uE, n_ie_obs_fit, n_count_days_fit, mean_log_day_effect_all_count_days,
         prior_mu, prior_sd, shared_tau_sigma, log_sd_days_pooled, tau_derivation, tau_likelihood,
         F_total, F_smooth, F_cover, F_gear_A, F_time_A, F_gear_B, F_time_B,
         survey_start, survey_end, on_first_interval, balance_at_end)

# ---- 3. THE FACTOR MEANS, THE RULE, AND THE DESK APPROXIMATION ----------------------------
accounts <- function(f, g) is.finite(f) & is.finite(g) & sign(f) == sign(g) & abs(f) >= abs(g) / 2
factors <- by_day |> group_by(season_window, run, fit) |>
  summarise(
    n_ie_days = n(), n_ie_obs_fit = first(n_ie_obs_fit), n_with_counts = sum(is.finite(F_cover)),
    tau_bar_fit = first(tau_bar), R_G = first(R_G), sigma_IE = first(sigma_IE), sigma_uE = first(sigma_uE),
    prior_mu = first(prior_mu), prior_sd = first(prior_sd), shared_tau_sigma = first(shared_tau_sigma), log_sd_days_pooled = first(log_sd_days_pooled),
    tau_derivation_geomean = gm(tau_derivation), tau_derivation_ratio_of_sums = sum(arrivals) / sum(Pbar),
    tau_likelihood_geomean = gm(tau_likelihood), tau_likelihood_ratio_of_sums = sum(arrivals) / sum(lambda),
    # the rule's scale: mean logs over the I/E days
    g = mlog(F_total), g_se = selog(F_total),
    f_smooth = mlog(F_smooth), f_smooth_se = selog(F_smooth), f_cover = mlog(F_cover), f_cover_se = selog(F_cover),
    f_gear_A = mlog(F_gear_A), f_time_A = mlog(F_time_A),
    f_gear_B = mlog(F_gear_B), f_gear_B_se = selog(F_gear_B), f_time_B = mlog(F_time_B), f_time_B_se = selog(F_time_B),
    f_smooth_plus_time = mlog(F_smooth * F_time_B), f_smooth_plus_time_se = selog(F_smooth * F_time_B),
    f_all_but_cover = mlog(F_total / F_cover), f_all_but_cover_se = selog(F_total / F_cover),
    f_time_B_days_off_schedule = mlog(F_time_B[counts_off_schedule > 0]), f_time_B_days_on_schedule = mlog(F_time_B[counts_off_schedule == 0]),
    n_days_off_schedule = sum(counts_off_schedule > 0),
    # beside the rule: the same four as ratios of sums (their product is sum(lambda) / sum(Pbar) exactly)
    R_total = sum(lambda) / sum(Pbar), R_smooth = sum(lambda) / sum(Gtot / R_G), R_cover = sum(Gtot) / sum(G20),
    R_gear_B = sum(G20 / R_G) / sum(P_c_B), R_time_B = sum(P_c_B) / sum(Pbar),
    # the count schedule and the per-day effort effect on the I/E days
    counts_off_schedule = sum(counts_off_schedule), n_counts = sum(n_counts),
    mean_log_day_effect_ie_days = mlog(day_effect), mean_log_day_effect_all_count_days = first(mean_log_day_effect_all_count_days),
    on_first_interval_over_arrivals = sum(on_first_interval) / sum(arrivals),
    balance_at_end_over_arrivals = sum(pmax(balance_at_end, 0)) / sum(arrivals), .groups = "drop") |>
  mutate(
    material_smooth = abs(f_smooth) >= RULE_MATERIAL, material_cover = abs(f_cover) >= RULE_MATERIAL,
    material_gear_A = abs(f_gear_A) >= RULE_MATERIAL, material_time_A = abs(f_time_A) >= RULE_MATERIAL,
    material_gear_B = abs(f_gear_B) >= RULE_MATERIAL, material_time_B = abs(f_time_B) >= RULE_MATERIAL,
    accounts_smooth = accounts(f_smooth, g), accounts_cover = accounts(f_cover, g),
    accounts_gear_A = accounts(f_gear_A, g), accounts_time_A = accounts(f_time_A, g),
    accounts_gear_B = accounts(f_gear_B, g), accounts_time_B = accounts(f_time_B, g),
    exp_f_gear_A = exp(f_gear_A), exp_f_gear_B = exp(f_gear_B),
    # THE DESK APPROXIMATION of tau_bar: the prior and the I/E likelihood combined by precision on the log
    # scale; a day's log(A / lambda) has variance sigma_IE^2 + shared_tau_sigma^2 about log(tau_bar).
    w_prior = 1 / prior_sd^2, w_lik = ifelse(n_ie_obs_fit > 0, n_ie_obs_fit / (sigma_IE^2 + shared_tau_sigma^2), 0),
    tau_bar_desk = exp((w_prior * log(prior_mu) + w_lik * log(tau_likelihood_geomean)) / (w_prior + w_lik)),
    # ... with the Float 20 share of the day's gear in the I/E prediction (the likelihood's centre moves by f_cover)
    tau_bar_desk_cover = exp((w_prior * log(prior_mu) + w_lik * (log(tau_likelihood_geomean) + f_cover)) / (w_prior + w_lik)),
    tau_bar_desk_cover_over_fit = tau_bar_desk_cover / tau_bar_desk,   # the effort and the catch scale by this (desk over desk)
    tau_bar_prior_over_fit = prior_mu / tau_bar_fit,                   # ... and by this with the I/E stream dropped (the centre over the fitted value)
    # ... and with the day's OWN turnover (arrivals over the tally's presence at the season's count hours) observed on L[d]
    # instead (D51 option (d)): the likelihood's centre is the derivation's reading on the fit's own days, a day's log ratio
    # carrying the between-day log-SD of the pooled derivation
    w_own = ifelse(n_ie_obs_fit > 0, n_ie_obs_fit / log_sd_days_pooled^2, 0),
    tau_bar_desk_own = exp((w_prior * log(prior_mu) + w_own * log(tau_derivation_geomean)) / (w_prior + w_own)),
    tau_bar_desk_own_over_fit = tau_bar_desk_own / tau_bar_desk) |>
  select(-w_prior, -w_lik, -w_own)

# ---- 4. EVERY I/E DAY IN THE WORKBOOK, COUNT BY COUNT (the two factors that need no fit) --
all_counts <- bind_rows(lapply(2022:2025, function(y) {
  r <- tryCatch(quietly(fetch_crab_data(with_window(season_window(y)))), error = function(e) NULL)
  if (is.null(r)) NULL else r$shore_effort |> mutate(season = sprintf("%d-%02d", y, (y + 1) %% 100)) }))
hour_class <- function(h) ifelse(h < 10, "before 10:00", ifelse(h < 15, "10:00 to 15:00", "15:00 and after"))
by_count <- all_counts |> filter(event_date %in% tally_days$event_date) |> add_presence() |>
  transmute(season, event_date, count_sequence, count_hour = round(count_hour, 3), hour_class = hour_class(count_hour),
            f20_gear, f17_gear, f17_source, presence_workbook_bin = P_A, presence_balance = P_B, position,
            gear_per_present = ifelse(is.finite(P_B) & P_B > 0, f20_gear / P_B, NA_real_))
R_G_season <- factors |> filter(grepl("all_gear", fit)) |> transmute(season = season_window, R_G_all_gear_fit = R_G)
ros <- function(d) if (nrow(d) && sum(d$presence_balance) > 0) sum(d$f20_gear) / sum(d$presence_balance) else NA_real_
by_season <- by_count |> filter(is.finite(presence_balance)) |> group_by(season) |>
  group_modify(function(d, k) tibble(
    n_ie_days_with_count = n_distinct(d$event_date), n_counts = nrow(d),
    gear_per_present_all = ros(d),
    gear_per_present_10_to_15 = ros(d[d$hour_class == "10:00 to 15:00", ]), n_10_to_15 = sum(d$hour_class == "10:00 to 15:00"),
    gear_per_present_before_10 = ros(d[d$hour_class == "before 10:00", ]), n_before_10 = sum(d$hour_class == "before 10:00"),
    gear_per_present_15_on = ros(d[d$hour_class == "15:00 and after", ]), n_15_on = sum(d$hour_class == "15:00 and after"))) |>
  ungroup() |>
  left_join(by_count |> filter(f17_source == "observed") |> group_by(season) |>
              summarise(n_f17_observed_ie_days = n(), two_float_ratio_ie_days = sum(f20_gear + f17_gear) / sum(f20_gear), .groups = "drop"), by = "season") |>
  left_join(all_counts |> filter(f17_source == "observed") |> group_by(season) |>
              summarise(n_f17_observed_all_days = n(), two_float_ratio_all_days = sum(f20_gear + f17_gear) / sum(f20_gear), .groups = "drop"), by = "season") |>
  left_join(R_G_season, by = "season") |>
  mutate(F_gear_ratio_of_sums_10_to_15 = gear_per_present_10_to_15 / R_G_all_gear_fit,
         F_gear_ratio_of_sums_all = gear_per_present_all / R_G_all_gear_fit)
# IS THE TALLY ONE FLOAT'S OR TWO? If the WDF20 tally counted the crabbers of Floats 17-21 as well, its presence
# would be both floats' gear over R_G, and Float 20's gear per present crabber would read R_G over the two-float
# ratio, about 1.05, not R_G. A bootstrap over I/E days (the counts of a day stay together) says how far the
# measured value at the regular count hours is from that alternative.
boot_gpp <- function(s, n_boot = 4000, seed = 20261119L) {
  d <- by_count[by_count$season == s & is.finite(by_count$presence_balance) & by_count$hour_class == "10:00 to 15:00", , drop = FALSE]
  days <- unique(d$event_date)
  if (length(days) < 3) return(c(se = NA_real_, lo = NA_real_, hi = NA_real_))
  idx <- split(seq_len(nrow(d)), d$event_date)
  b <- bss_with_seed(seed, replicate(n_boot, { i <- unlist(idx[sample(length(idx), replace = TRUE)], use.names = FALSE)
                                                sum(d$f20_gear[i]) / sum(d$presence_balance[i]) }))
  c(se = stats::sd(b), lo = unname(stats::quantile(b, 0.025)), hi = unname(stats::quantile(b, 0.975)))
}
bt <- t(vapply(by_season$season, boot_gpp, numeric(3)))
by_season <- by_season |>
  mutate(gpp_10_to_15_boot_se = bt[, "se"], gpp_10_to_15_lo95 = bt[, "lo"], gpp_10_to_15_hi95 = bt[, "hi"],
         gpp_if_tally_covered_both_floats = R_G_all_gear_fit / two_float_ratio_all_days,
         se_above_both_floats_alternative = (gear_per_present_10_to_15 - gpp_if_tally_covered_both_floats) / gpp_10_to_15_boot_se)

# ---- 5. PRINT ----------------------------------------------------------------------------
op <- options(width = 240); on.exit(options(op), add = TRUE)
r2 <- function(x) round(x, 2); r3 <- function(x) round(x, 3)
cat(sprintf("\nTHE TALLY: %d WDF20 I/E days; the workbook's flow column is 0 on the first interval of %d days and on the last interval of %d days where the running balance is not.\n",
            nrow(tally_days), sum(tally_days$first_row_zeroed), sum(tally_days$last_row_zeroed)))
for (s in names(RUNS)) {
  k <- res[[s]]$f17; st <- res[[s]]$st
  cat(sprintf("\n== %s (%s): derivation %.3f pooled over %d I/E days, %.3f on the %d in-window days; Float 17-21 observed beside %d of %d Float 20 counts, all-pairs ratio %.3f\n",
              s, RUNS[[s]]$label, st$tau, st$n_days, st$tau_in_window, st$n_days_in_window, k$n_observed, k$n_counts, k$r_all))
  sb <- res[[s]]$st_b; sc <- res[[s]]$sched
  cat(sprintf("   on the tally's running balance in place of the workbook's flow column: %.3f pooled, %.3f in-window (%+.2f%% and %+.2f%%)\n",
              sb$tau, sb$tau_in_window, 100 * (sb$tau / st$tau - 1), 100 * (sb$tau_in_window / st$tau_in_window - 1)))
  cat(sprintf("   counts taken before 10:00 or at 15:00 and after: %d of %d on I/E days (%.0f%%), %d of %d on the other count days (%.0f%%)\n",
              sc$off[sc$ie_day], sc$n[sc$ie_day], 100 * sc$off[sc$ie_day] / sc$n[sc$ie_day], sc$off[!sc$ie_day], sc$n[!sc$ie_day], 100 * sc$off[!sc$ie_day] / sc$n[!sc$ie_day]))
  print(as.data.frame(by_day |> filter(season_window == s) |>
    transmute(fit = sub("shore_", "", fit), date = format(event_date), A = arrivals, Pbar = round(Pbar, 1), lambda = round(lambda, 1),
              G20 = round(G20, 1), G17 = round(G17, 1), f17obs = f17_observed, n = n_counts, off = counts_off_schedule,
              P_A = round(P_c_A, 1), P_B = round(P_c_B, 1), day_eff = r2(day_effect), tau_der = r2(tau_derivation), tau_lik = r2(tau_likelihood),
              F_tot = r2(F_total), F_smooth = r2(F_smooth), F_cover = r2(F_cover), F_gear_A = r2(F_gear_A), F_time_A = r2(F_time_A),
              F_gear_B = r2(F_gear_B), F_time_B = r2(F_time_B), hours = count_hours)), row.names = FALSE)
}
cat("\n== THE RULE'S SCALE: mean logs over the I/E days of each fit (se = between-day SD / sqrt(n)); A = as first computed, B = the repaired presence\n")
print(as.data.frame(factors |> transmute(season_window, fit = sub("shore_", "", fit), n = n_ie_days, n_fit = n_ie_obs_fit,
        g = r3(g), g_se = r3(g_se), f_smooth = r3(f_smooth), se_s = r3(f_smooth_se), f_cover = r3(f_cover), se_c = r3(f_cover_se),
        f_gear_A = r3(f_gear_A), f_time_A = r3(f_time_A), f_gear_B = r3(f_gear_B), se_gB = r3(f_gear_B_se), f_time_B = r3(f_time_B), se_tB = r3(f_time_B_se),
        `smooth+time` = r3(f_smooth_plus_time), se_st = r3(f_smooth_plus_time_se),
        all_but_cover = r3(f_all_but_cover), se_abc = r3(f_all_but_cover_se))), row.names = FALSE)
cat(sprintf("\n== THE RULE, applied mechanically (material at |f| >= %.2f; accounts for the gap at |f| >= |g| / 2 with g's sign); fits with an I/E stream only\n", RULE_MATERIAL))
flag <- function(m, a) ifelse(a, "ACCOUNTS", ifelse(m, "material", "-"))
print(as.data.frame(factors |> filter(n_ie_obs_fit > 0) |> transmute(season_window, fit = sub("shore_", "", fit), g = r3(g), half_g = r3(g / 2),
        F_smooth = flag(material_smooth, accounts_smooth), F_cover = flag(material_cover, accounts_cover),
        F_gear_A = flag(material_gear_A, accounts_gear_A), F_time_A = flag(material_time_A, accounts_time_A),
        F_gear_B = flag(material_gear_B, accounts_gear_B), F_time_B = flag(material_time_B, accounts_time_B),
        `exp(f_gear_A)` = r3(exp_f_gear_A), `exp(f_gear_B)` = r3(exp_f_gear_B))), row.names = FALSE)
cat("\n== BESIDE THE RULE: the same factors as ratios of sums (product = sum(lambda) / sum(Pbar)), and the two readings\n")
print(as.data.frame(factors |> transmute(season_window, fit = sub("shore_", "", fit), n = n_ie_days,
        tau_bar_fit = r3(tau_bar_fit), prior = r3(prior_mu), der_geomean = r3(tau_derivation_geomean), der_ros = r3(tau_derivation_ratio_of_sums),
        lik_geomean = r3(tau_likelihood_geomean), lik_ros = r3(tau_likelihood_ratio_of_sums),
        R_total = r3(R_total), R_smooth = r3(R_smooth), R_cover = r3(R_cover), R_gear_B = r3(R_gear_B), R_time_B = r3(R_time_B))), row.names = FALSE)
cat("\n== THE COUNT SCHEDULE AND THE PER-DAY EFFORT EFFECT on the I/E days; the survey's reach (H5, reported, not judged)\n")
print(as.data.frame(factors |> transmute(season_window, fit = sub("shore_", "", fit), counts = n_counts, off_schedule = counts_off_schedule,
        sigma_uE = r3(sigma_uE), mean_log_day_effect_ie = r3(mean_log_day_effect_ie_days), mean_log_day_effect_all = r3(mean_log_day_effect_all_count_days),
        days_off = n_days_off_schedule, f_time_B_off = r3(f_time_B_days_off_schedule), f_time_B_on = r3(f_time_B_days_on_schedule),
        on_first_over_A = r3(on_first_interval_over_arrivals), still_present_at_end_over_A = r3(balance_at_end_over_arrivals))), row.names = FALSE)
long <- tally_days |> filter(survey_end >= 20)
if (nrow(long)) cat(sprintf("   the %d day(s) the tally ran past 20:00 (%s): %s of %s arrivals came at 18:00 or later (%.0f%%)\n", nrow(long),
        paste(format(long$event_date), collapse = ", "), sum(long$arrivals_after_18), sum(long$arrivals), 100 * sum(long$arrivals_after_18) / sum(long$arrivals)))
cat("\n== DESK APPROXIMATION of tau_bar (precision-weighted prior and I/E likelihood, log scale; NOT a fit)\n")
print(as.data.frame(factors |> filter(n_ie_obs_fit > 0) |> transmute(season_window, fit = sub("shore_", "", fit), n_fit = n_ie_obs_fit, sigma_IE = r3(sigma_IE),
        prior = r3(prior_mu), prior_sd = r3(prior_sd), lik_centre = r3(tau_likelihood_geomean), tau_bar_fit = r3(tau_bar_fit), tau_bar_desk = r3(tau_bar_desk),
        with_float20_share = r3(tau_bar_desk_cover), effort_ratio_cover = r3(tau_bar_desk_cover_over_fit), effort_ratio_no_ie = r3(tau_bar_prior_over_fit),
        own_turnover = r3(tau_bar_desk_own), effort_ratio_own = r3(tau_bar_desk_own_over_fit), log_sd_days = r3(log_sd_days_pooled))), row.names = FALSE)
cat("\n== EVERY I/E DAY WITH A SAME-DAY FLOAT 20 COUNT, by season: gear per present crabber (ratio of sums, the running balance) and the two-float ratio (observed Float 17-21 counts)\n")
print(as.data.frame(by_season |> mutate(across(where(is.double), r3))), row.names = FALSE)
cat("\n== IS THE TALLY ONE FLOAT'S OR TWO? Float 20 gear per present crabber at 10:00 to 15:00 against R_G (one float) and R_G over the two-float ratio (both floats)\n")
print(as.data.frame(by_season |> filter(is.finite(R_G_all_gear_fit)) |> transmute(season, counts = n_10_to_15, measured = r3(gear_per_present_10_to_15),
        boot_se = r3(gpp_10_to_15_boot_se), lo95 = r3(gpp_10_to_15_lo95), hi95 = r3(gpp_10_to_15_hi95), if_one_float = r3(R_G_all_gear_fit),
        if_both_floats = r3(gpp_if_tally_covered_both_floats), se_above_both = r2(se_above_both_floats_alternative))), row.names = FALSE)

# ---- 6. WRITE ----------------------------------------------------------------------------
source("03_R_functions/batch_verdict_helpers.R")   # merge_csv_by
out_file <- function(kind) file.path("05_output", sprintf("shore_turnover_d50_2026-11-18_%s.csv", kind))
merge_csv_by(as.data.frame(by_day |> mutate(event_date = format(event_date))), out_file("by_day"), key = c("season_window", "fit", "event_date"))
merge_csv_by(as.data.frame(factors), out_file("factors"), key = c("season_window", "fit"))
merge_csv_by(as.data.frame(by_count |> mutate(event_date = format(event_date))), out_file("by_count"), key = c("event_date", "count_sequence"))
merge_csv_by(as.data.frame(by_season), out_file("by_season"), key = "season")
cat("\nwrote 05_output/shore_turnover_d50_2026-11-18_{by_day,factors,by_count,by_season}.csv\n")
cat("READ BY THE RULES OF CHANGE_REGISTER B78. The rule's verdicts are the table headed THE RULE; everything headed BESIDE THE RULE was added after the first table was seen.\n")
