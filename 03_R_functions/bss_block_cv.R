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
# bss_block_cv.R  (2026-09-26; CHANGE_REGISTER D31 -> B39)
#
# LEAVE-ONE-BLOCK-OUT CROSS-VALIDATION OF THE EFFORT STREAMS, from posterior draws.
#
# WHY THIS EXISTS. Observation-level PSIS-LOO (write_loo_diagnostics(), the ladders'
# rule 4) scores each SAMPLED count with the rest of the sampled counts still in the
# fit. A day covariate on the latent effort process earns its place on the UNSAMPLED
# days, where the AR process interpolates and the covariate is the only thing that says
# a day was not ordinary. Pointwise LOO never holds those days out, so it cannot measure
# that. The 2026-09-25 marine hazard ladder made this concrete: the boat SCA term was
# identified at eight standard errors and gained +7.6 nats on the trailer stream at 1.46
# paired SE, and the rule could say nothing about the interpolation.
#
# WHAT THIS DOES. Hold out a BLOCK, one calendar week of counts, and score the block's
# JOINT predictive density under the posterior without it. The leave-block-out posterior
# is approximated by Pareto-smoothed importance sampling on the block's joint
# log-likelihood (Vehtari, Gelman & Gabry 2017; the loo package), so no refit is needed,
# and the approximation's reliability is reported per block as its Pareto k. Blocks with
# k > k_max (0.7) are UNRELIABLE and are excluded from any comparison, and the share
# excluded is part of the answer: if it is large, the honest next step is an exact refit
# of those blocks, not a number.
#
# WHAT "HELD OUT" MEANS HERE, stated exactly because it is easy to overstate. The
# leave-out set for a week is EVERY effort observation of that week in the streams the
# caller passes: for the boat that is the trailer counts AND the OSP port counts (35% of
# the 2024-25 all-gear trailer days also carry an OSP count, and 24 of 37 trailer weeks
# do), scored stream by stream under the SAME importance weights. Holding out one stream
# alone would leave the other stream's counts anchoring lambda_E on the very days being
# "predicted", which is closer to pointwise LOO than to the interpolation test. What still
# anchors a held-out week: the shore fit's I/E observations (no log_lik_ie exists in the
# Stan models; four observations on 2024-25), and, under estimate_cpue_density = 1 only
# (off in production), the interviews through gamma_C. So a held-out week is unsampled
# for the effort streams that were passed, not for the whole model, and the header of a
# verdict should say which streams were.
#
# WHAT IT CANNOT DO. It tests interpolation between sampled weeks, the model's actual job
# on 2024-25 (79.8% of boat days unsampled, interspersed), not extrapolation past the last
# sampled day. And for a fit whose AR period equals the block (the shore's weekly AR with
# weekly blocks), a held-out week removes the period's only anchor and the importance
# ratios are expected to be heavy-tailed: expect many k > 0.7 there and read the reliable
# share before the difference.
#
# ENTRY POINTS
#   bss_effort_loglik(draws, stan_data, stream)   per-observation log-likelihood matrix
#       [draws x n] for "gear", "trailer" or "osp", rebuilt from lambda_E_S, r_E, R_G,
#       R_G_boat (or R_T), L_out, kappa_OSP, r_OSP EXACTLY as both Stan models form
#       log_lik_* (the gear-resolved model sums lambda_E_S over its gear groups; with
#       G = 1, the pooled model, the sum is the single column). Used when only the saved
#       ppc_draws_<fit>.rds is at hand (a committed run, no stanfit). With a live stanfit,
#       rstan::extract(fit, "log_lik_<stream>") is the same matrix and needs no
#       reconstruction. A reconstruction must be verified against
#       loo_pointwise_<stream>_<fit>.csv's lpd column before it is trusted
#       (bss_block_cv_check()).
#   bss_block_psis_loo(ll, block, k_max, weight_ll, weight_block)   one row per block of
#       the scored stream: n_obs, the in-sample joint lpd, the PSIS leave-block-out elpd,
#       the Pareto k, reliable. `weight_ll` / `weight_block` add the OTHER streams'
#       observations of each block to the leave-out set (the importance ratios) without
#       scoring them.
#   bss_block_cv_compare(tab_a, tab_b, k_max)      paired difference over the blocks both
#       tables have with k <= k_max, with its SE, the SE ratio, the excluded count, the
#       count of positive weeks and the median, and `identical` when the two tables differ
#       by summation noise only (2026-09-27; see the function).
#   bss_block_cv_weeks(tab_a, tab_b, k_max)        the per-week rows behind that comparison.
#   bss_block_cv_fit(streams, k_max)               every stream of one fit under the joint
#       leave-out set, plus, with more than one stream, the "joint" table: the week's joint
#       predictive density of ALL its held-out effort observations (2026-09-27).
#   write_block_cv_diagnostics(fit, ...)           the in-driver call: every effort stream
#       of one fit, each scored under the joint leave-out set of all of them, written as
#       loo_block_<stream>_<label>.csv (and loo_block_joint_<label>.csv for a fit with
#       more than one stream). On a fit that carries the per-day effort effect the tables
#       are the day-integrated ones (below, D54) and the conditional ones are written beside
#       them as loo_blockcond_<stream>_<label>.csv, with blockcv_checks_<label>.csv.
#   bss_block_integrate_days(mu0, y, eday, sig, size)  the day-integrated log-likelihood,
#       one column per sampled day (D54).
#   bss_block_day_streams(draws, stan_data, event_dates)  the gear and trailer streams of a
#       fit with the per-day effect, as day-integrated log-likelihoods with their weekly
#       labels and the checks; from rstan::extract() or from saved ppc draws (D54).
#
# THE PER-DAY EFFORT EFFECT (A39) AND WHY THE TABLE IS INTEGRATED ON SUCH A FIT (D54, B85;
# 2026-11-25). Since A39 a gear or trailer count's expectation carries its sampled day's own
# deviation, exp(sigma_uE * uE_raw[day] - sigma_uE^2 / 2), and log_lik_gear / log_lik_trailer
# are CONDITIONAL on it. A held-out week's counts are exactly what identify those days'
# deviations, so importance sampling on the conditional likelihood asks the posterior to stand
# in for one in which those deviations have reverted to their prior: the ratios are
# heavy-tailed and most weeks fail the Pareto-k bar (12 of 42 joint weeks on the authoritative
# boat all-gear fit, against 38 of 42 on the render before A39). It is not a sampler fault.
# The repair integrates each held-out day's deviation out of both the importance ratio and
# the score, which is the question an unsampled day poses (integrated importance sampling:
# Li, Qiu, Zhang and Feng 2016, Statistics and Computing 26(4), 881-897; the conditional
# against marginal distinction of Merkle, Furr and Rabe-Hesketh 2019, Psychometrika 84(3),
# 802-829). For a sampled day with counts y_1..y_m and count expectations mu0_j WITHOUT the
# day factor,
#     p(y | theta) = integral of prod_j NB2(y_j | mu0_j exp(sigma z - sigma^2 / 2), r_E) N(z | 0, 1) dz,
# computed by the trapezoid rule on z in [-8, 8] at step 0.025 (641 nodes) and checked against
# every other node. The OSP likelihood carries no day deviation and enters unchanged. The scored
# unit of the gear and trailer rows becomes the sampled DAY (n_obs counts days on those rows);
# the joint row adds the OSP counts. Prototyped in 06_diagnostics/desk_advisory_block_cv_2026-11-19.R
# (B80's R3a); accepted against an exact refit of two held-out weeks to 0.022 and 0.020 nats (B83);
# this library version reproduces that prototype's tables (B85). Two checks are made on every
# fit and written to blockcv_checks_<label>.csv: (i) the conditional log-likelihood rebuilt from
# mu0 and the fitted deviations equals the fit's own log_lik to 1e-6 (else the reconstruction is
# wrong, the integrated table is not trusted, and the conditional one is written as loo_block_*
# with a run warning); (ii) halving the grid moves no reliable week by 0.01 nats (else a run
# warning). WHAT IT DOES NOT FIX: a fit whose AR period equals the block (the shore all-gear
# fit's weekly AR) still loses its period's only anchor when a week is held out; there the
# integral returns the table to its state before A39 and no further.
# Blocks are labelled by bss_block_weeks(dates): the ISO week ("2025-W03"). A week is a
# choice: the resolution the shore AR runs at and a plausible gap between sampled boat
# days; pass another `block` vector to test another one.
#
# WHAT THE FIRST RUN SAID (2026-09-26/27, the five marine hazard rungs; Section 1y). The
# reconstruction matched every committed pointwise file to 5e-5. The boat all-gear fit at
# monthly AR was evaluable (34 of 37 trailer weeks, 25 of 29 OSP weeks reliable); the
# shore's weekly-AR fits were not (8 to 16 of 38), as the paragraph above predicted. The
# boat SCA term gained +18.9 nats on the held-out OSP counts (3.9 paired SE, 21 of 25
# weeks positive) and +7.9 on the trailer counts (1.2 SE; 1.6 SE over the weeks the OSP
# stream also covers), +26.8 on the two together (2.9 SE, the sum of the two streams' weekly
# differences over the 38 weeks reliable in every stream present, recomputed from the
# per-week tables; the exact joint row below is what the next run carries). The trailer's
# 13 winter weeks, the only held-out evidence for the months where the term moves the
# estimate most, summed to -0.6 +- 3.9: uninformative, not negative. The bar tick beyond
# the archive: -2.0 trailer, -0.4 OSP, nothing.
###############################################################################

# NULL-only, guarded, the same idiom as the sibling helpers: the drivers source this
# folder after library(tidyverse), so an unconditional definition would shadow rlang's.
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ISO week label per date, "YYYY-Www". A block is a calendar week, Monday to Sunday.
bss_block_weeks <- function(dates) format(as.Date(dates), "%G-W%V")

# `draws` is a list as save_bss_ppc_draws() stores it (lambda_E_S as [draws, S, D, G],
# scalars as vectors or 1-d arrays, L_out as [draws, D]) or as rstan::extract() returns
# it. `stan_data` is the fit's Stan data list. Returns a list:
#   ll     [draws x n] log-likelihood of every observation of the stream
#   y      the observations, day the day index, section the section index
# Mirrors the generated quantities of crab_bss_pooled.stan (G = 1) and
# crab_bss_gear_resolved.stan (sum over gear groups):
#   gear     Gear_I[i] ~ NB2(sum_g lambda_E_S[section][day, g] * R_G, r_E)
#   trailer  T_I[i]    ~ NB2(sum_g lambda_E_S[section][day, g] / R_G_boat, r_E)   (or * R_T)
#   osp      OSP_I[i]  ~ NB2((sum_g lambda_E_S[section][day, g] / R_G_boat)
#                            * (L[day] if osp_scale_is_tau else kappa_OSP), r_OSP)
bss_effort_loglik <- function(draws, stan_data, stream = c("gear", "trailer", "osp"),
                              draw_idx = NULL) {
  stream <- match.arg(stream)
  if (is.null(stan_data)) stop("stan_data is NULL", call. = FALSE)
  lam <- draws$lambda_E_S
  if (is.null(lam)) stop("draws carry no lambda_E_S", call. = FALSE)
  d <- dim(lam)
  if (is.null(d) || length(d) != 4L)
    stop(sprintf("lambda_E_S must be the [draws, S, D, G] array rstan::extract() returns (got dims %s)",
                 if (is.null(d)) "none" else paste(d, collapse = " x ")), call. = FALSE)
  # [[ ]] not $: a list's $ partial-matches, and stan_data$G would read Gear_n
  D <- stan_data[["D"]] %||% d[3]
  if (!identical(as.integer(d[3]), as.integer(D)))
    stop(sprintf("lambda_E_S has %d days; stan_data$D is %d", d[3], as.integer(D)), call. = FALSE)
  if (!is.null(stan_data[["S"]]) && !identical(as.integer(d[2]), as.integer(stan_data[["S"]])))
    stop(sprintf("lambda_E_S has %d sections; stan_data$S is %d", d[2], as.integer(stan_data[["S"]])), call. = FALSE)
  if (!is.null(stan_data[["G"]]) && !identical(as.integer(d[4]), as.integer(stan_data[["G"]])))
    stop(sprintf("lambda_E_S has %d gear groups; stan_data$G is %d", d[4], as.integer(stan_data[["G"]])), call. = FALSE)
  nd_all <- d[1]
  idx <- if (is.null(draw_idx)) seq_len(nd_all) else draw_idx
  num <- function(x) if (is.null(x)) NULL else as.numeric(x)[idx]
  spec <- switch(stream,
    gear    = list(n = stan_data$Gear_n %||% 0L, y = stan_data$Gear_I, day = stan_data$day_Gear, section = stan_data$section_Gear),
    trailer = list(n = stan_data$T_n %||% 0L, y = stan_data$T_I, day = stan_data$day_T, section = stan_data$section_T),
    osp     = list(n = stan_data$OSP_n %||% 0L, y = stan_data$OSP_I, day = stan_data$day_OSP, section = stan_data$section_OSP))
  n <- as.integer(spec$n)
  if (n == 0L) return(list(ll = matrix(numeric(0), nrow = length(idx), ncol = 0), y = integer(0),
                           day = integer(0), section = integer(0), stream = stream))
  sec <- if (is.null(spec$section)) rep(1L, n) else as.integer(spec$section)[seq_len(n)]
  day <- as.integer(spec$day)[seq_len(n)]; y <- as.integer(spec$y)[seq_len(n)]
  if (any(is.na(day)) || any(day < 1L | day > d[3])) stop(sprintf("%s day indices outside 1..%d", stream, d[3]), call. = FALSE)
  if (any(is.na(sec)) || any(sec < 1L | sec > d[2])) stop(sprintf("%s section indices outside 1..%d", stream, d[2]), call. = FALSE)
  if (stream == "gear") {
    mult <- num(draws$R_G); size <- num(draws$r_E)
    if (is.null(mult)) stop("gear stream needs R_G draws", call. = FALSE)
  } else {
    tp <- if (!is.null(draws$R_G_boat)) "R_G_boat" else if (!is.null(draws$R_T)) "R_T" else NA_character_
    if (is.na(tp)) stop("trailer / OSP stream needs R_G_boat (or R_T) draws", call. = FALSE)
    mult <- if (identical(tp, "R_T")) num(draws$R_T) else 1 / num(draws$R_G_boat)
    size <- if (stream == "trailer") num(draws$r_E) else num(draws$r_OSP)
    if (is.null(size)) stop(sprintf("%s stream needs %s draws", stream, if (stream == "trailer") "r_E" else "r_OSP"), call. = FALSE)
  }
  if (is.null(size)) stop(sprintf("%s stream needs its dispersion draws", stream), call. = FALSE)
  osp_L <- NULL; osp_k <- NULL
  if (stream == "osp") {
    if (identical(as.integer(stan_data$osp_scale_is_tau %||% 0L), 1L)) {
      osp_L <- draws$L_out
      if (is.null(osp_L)) stop("osp stream under osp_scale_is_tau needs L_out draws", call. = FALSE)
      osp_L <- as.matrix(osp_L)
      if (ncol(osp_L) != d[3]) stop("L_out does not have one column per day", call. = FALSE)
    } else {
      osp_k <- num(draws$kappa_OSP)
      if (is.null(osp_k)) stop("osp stream needs kappa_OSP draws", call. = FALSE)
    }
  }
  # D54 (2026-11-25): since A39 a gear or trailer count carries its sampled day's deviation. The
  # saved draws hold sigma_uE and uE_raw from B85 on; with them the rebuild is exact. Without
  # them (a run saved before B85) it is the count's likelihood at the LEVEL, which is not the
  # fit's log_lik, and bss_block_cv_check() will say so.
  dayf <- NULL
  if (stream != "osp" && identical(as.integer(stan_data$effort_day_effect %||% 0L), 1L) &&
      !is.null(draws$sigma_uE) && !is.null(draws$uE_raw) && length(draws$sigma_uE) > 0) {
    ed <- as.integer(if (stream == "gear") stan_data$Eday_Gear else stan_data$Eday_T)[seq_len(n)]
    sg <- as.numeric(draws$sigma_uE)[idx]; uu <- as.matrix(draws$uE_raw)[idx, , drop = FALSE]
    dayf <- function(i) exp(sg * uu[, ed[i]] - 0.5 * sg^2)
  }
  ll <- matrix(NA_real_, nrow = length(idx), ncol = n)
  for (i in seq_len(n)) {
    lam_i <- lam[idx, sec[i], day[i], , drop = FALSE]                 # [draws, 1, 1, G]
    lam_i <- rowSums(matrix(lam_i, nrow = length(idx)))              # sum over gear groups
    mu <- lam_i * mult
    if (!is.null(dayf)) mu <- mu * dayf(i)
    if (stream == "osp") mu <- mu * (if (!is.null(osp_L)) osp_L[idx, day[i]] else osp_k)
    ll[, i] <- stats::dnbinom(y[i], mu = pmax(mu, 1e-300), size = size, log = TRUE)
  }
  list(ll = ll, y = y, day = day, section = sec, stream = stream)
}

# log(mean(exp(x))) and log(sum(exp(x))), stable
.bcv_lmeanexp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
.bcv_lsumexp  <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# Self-check of a reconstruction against a run's committed loo_pointwise file: the same
# observations in the same order, and the pointwise lpd (log mean predictive density over
# draws) agreeing on EVERY observation. The file rounds lpd to 4 decimals (5e-5); the
# tolerance is an order of magnitude wider so floating-point summation order cannot fail
# it, and three orders of magnitude tighter than any formula error would produce. Returns
# list(ok, max_abs_diff, n, why).
bss_block_cv_check <- function(ll, pointwise_csv, tol = 5e-4) {
  bad <- function(why) list(ok = FALSE, max_abs_diff = NA_real_, n = ncol(ll), why = why)
  if (!file.exists(pointwise_csv)) return(bad("pointwise file absent"))
  pw <- utils::read.csv(pointwise_csv, stringsAsFactors = FALSE)
  if (!all(c("obs_index", "lpd") %in% names(pw))) return(bad("pointwise file lacks obs_index / lpd"))
  if (nrow(pw) != ncol(ll)) return(bad(sprintf("%d observations in the file, %d reconstructed", nrow(pw), ncol(ll))))
  if (!identical(as.integer(pw$obs_index), seq_len(ncol(ll)))) return(bad("obs_index is not 1..n in order"))
  lpd <- apply(ll, 2, .bcv_lmeanexp)
  dmax <- suppressWarnings(max(abs(lpd - as.numeric(pw$lpd))))
  if (!is.finite(dmax)) return(bad("non-finite lpd in the reconstruction or the file"))
  list(ok = isTRUE(dmax <= tol), max_abs_diff = dmax, n = ncol(ll), why = if (dmax <= tol) "" else sprintf("max |lpd diff| %.3g > %.0e", dmax, tol))
}

# Leave-one-block-out PSIS for ONE scored stream. `ll` [draws x n]; `block` length-n labels.
# `weight_ll` (a list of [draws x m_j] matrices) and `weight_block` (their labels) are the
# OTHER streams' observations: they join each block's leave-out set, so the importance
# ratios are those of the posterior without EVERY passed observation of the week, and the
# scored stream's block is then evaluated under those shared weights. Returns a data frame,
# one row per block of the scored stream (first-appearance order): block, n_obs (scored
# observations), n_leaveout (all observations removed for that block), lpd_block (the
# in-sample joint log predictive density of the scored observations), elpd_block (PSIS
# leave-block-out), p_eff = lpd - elpd, pareto_k, reliable (finite k <= k_max). A block
# with a non-finite log-likelihood anywhere in its leave-out set is scored NA with k = Inf
# rather than aborting the table. Draws are treated as exchangeable (r_eff = 1), the
# convention write_loo_diagnostics() uses for the pointwise LOO.
bss_block_psis_loo <- function(ll, block, k_max = 0.7, weight_ll = NULL, weight_block = NULL) {
  if (!requireNamespace("loo", quietly = TRUE)) stop("the 'loo' package is needed for block PSIS-LOO", call. = FALSE)
  block <- as.character(block)
  stopifnot(is.matrix(ll), length(block) == ncol(ll))
  if (!is.null(weight_ll)) {
    if (!is.list(weight_ll)) { weight_ll <- list(weight_ll); weight_block <- list(weight_block) }
    if (!is.list(weight_block) || length(weight_block) != length(weight_ll)) stop("weight_block must pair with weight_ll", call. = FALSE)
    for (j in seq_along(weight_ll)) {
      stopifnot(is.matrix(weight_ll[[j]]), nrow(weight_ll[[j]]) == nrow(ll), length(weight_block[[j]]) == ncol(weight_ll[[j]]))
      weight_block[[j]] <- as.character(weight_block[[j]])
    }
  }
  ub <- unique(block)
  rs <- function(m, j) if (length(j) == 1) m[, j] else rowSums(m[, j, drop = FALSE])
  scored <- matrix(vapply(ub, function(b) rs(ll, which(block == b)), numeric(nrow(ll))), nrow = nrow(ll))
  leave  <- scored
  n_leave <- as.integer(table(block)[ub])
  if (!is.null(weight_ll)) for (j in seq_along(weight_ll)) for (bi in seq_along(ub)) {
    jj <- which(weight_block[[j]] == ub[bi])
    if (length(jj)) { leave[, bi] <- leave[, bi] + rs(weight_ll[[j]], jj); n_leave[bi] <- n_leave[bi] + length(jj) }
  }
  finite <- apply(leave, 2, function(x) all(is.finite(x)))
  elpd <- rep(NA_real_, length(ub)); k <- rep(Inf, length(ub))
  if (any(finite)) {
    ps <- suppressWarnings(loo::psis(-leave[, finite, drop = FALSE], r_eff = rep(1, sum(finite))))
    lw <- stats::weights(ps, log = TRUE, normalize = TRUE)          # loo's method: log-normalised, [draws x blocks]
    k[finite] <- ps$diagnostics$pareto_k
    elpd[finite] <- vapply(seq_len(sum(finite)), function(b) .bcv_lsumexp(scored[, which(finite)[b]] + lw[, b]), numeric(1))
  }
  lpd <- apply(scored, 2, function(x) if (all(is.finite(x))) .bcv_lmeanexp(x) else NA_real_)
  data.frame(block = ub, n_obs = as.integer(table(block)[ub]), n_leaveout = n_leave, lpd_block = lpd, elpd_block = elpd,
             p_eff = lpd - elpd, pareto_k = k, reliable = is.finite(k) & k <= k_max, stringsAsFactors = FALSE)
}

# The per-week pairing behind a comparison: one row per block shared by the two tables,
# with both sides' elpd and Pareto k, the difference (B minus A) and whether the week counts
# (reliable in BOTH). This is what the runner writes out per pair so a reader can slice the
# held-out evidence by season, and what bss_block_cv_compare() sums.
bss_block_cv_weeks <- function(tab_a, tab_b, k_max = 0.7) {
  j <- merge(tab_a[, c("block", "n_obs", "elpd_block", "pareto_k")],
             tab_b[, c("block", "n_obs", "elpd_block", "pareto_k")], by = "block", suffixes = c("_a", "_b"))
  if (nrow(j) && any(j$n_obs_a != j$n_obs_b)) stop("the two block tables do not hold the same observations per block", call. = FALSE)
  j <- j[order(j$block), , drop = FALSE]
  j$diff <- j$elpd_block_b - j$elpd_block_a
  j$used <- is.finite(j$pareto_k_a) & is.finite(j$pareto_k_b) & j$pareto_k_a <= k_max & j$pareto_k_b <= k_max &
            is.finite(j$elpd_block_a) & is.finite(j$elpd_block_b)
  rownames(j) <- NULL
  j[, c("block", "n_obs_a", "elpd_block_a", "elpd_block_b", "diff", "pareto_k_a", "pareto_k_b", "used")]
}

# Paired comparison of two block tables of the SAME stream of the SAME fit under two
# configurations (same blocks, same observations): B against A. Only blocks reliable in
# BOTH count. One row: n_blocks, n_used, n_dropped, reliable_share, diff (sum over used
# blocks of elpd_B - elpd_A), se (sqrt(n_used) * sd of the per-block differences), ratio,
# and since 2026-09-27: max_abs_diff, n_positive, median_diff, identical.
#
# IDENTICAL TABLES (2026-09-27). Two rungs whose fits are the same model (a shore fit under a
# boat-only term, say) produce block tables equal to floating-point summation order, about
# 1e-14 nats per week. Summed and divided by a standard error of the same size, that noise
# read as "-2.31 SE" and the 2026-09-26 runner FAILED it (CHANGE_REGISTER C row). A
# difference smaller than `tol_identical` on EVERY used week is therefore reported as
# identical = TRUE with diff 0 and ratio NA: there is nothing to compare, and a verdict on
# it is a verdict on rounding. 1e-8 nats is 1e6 times the noise and 1e6 times below any
# difference a changed likelihood could produce over a week of counts.
bss_block_cv_compare <- function(tab_a, tab_b, k_max = 0.7, tol_identical = 1e-8) {
  j <- bss_block_cv_weeks(tab_a, tab_b, k_max = k_max)
  use <- j$used
  d <- j$diff[use]
  n_used <- sum(use)
  mad <- if (n_used) max(abs(d)) else NA_real_
  ident <- isTRUE(n_used > 0) && isTRUE(mad < tol_identical)
  se <- if (n_used >= 2 && !ident) sqrt(n_used) * stats::sd(d) else NA_real_
  data.frame(n_blocks = nrow(j), n_used = n_used, n_dropped = nrow(j) - n_used,
             reliable_share = if (nrow(j)) n_used / nrow(j) else NA_real_,
             diff = if (!n_used) NA_real_ else if (ident) 0 else sum(d), se = se,
             ratio = if (n_used >= 2 && !ident && is.finite(se) && se > 0) sum(d) / se else NA_real_,
             k_max = k_max, max_abs_diff = mad,
             n_positive = if (n_used) sum(d > tol_identical) else NA_integer_,
             median_diff = if (!n_used) NA_real_ else if (ident) 0 else stats::median(d),
             identical = ident, stringsAsFactors = FALSE)
}

# One line for a report or a verdict table.
bss_block_cv_str <- function(cmp, label_a = "A", label_b = "B") {
  if (is.null(cmp) || !nrow(cmp)) return("block CV not computable")
  if (isTRUE(cmp$identical))
    return(sprintf("%s -> %s: identical block tables on all %d used weeks (max |diff| %.1e nats): the two rungs hold the same fit",
                   label_a, label_b, cmp$n_used, cmp$max_abs_diff))
  sprintf("%s -> %s: %d of %d weeks reliable (k <= %.1f in both); held-out elpd diff %+.1f, SE %.1f (%s SE)%s",
          label_a, label_b, cmp$n_used, cmp$n_blocks, cmp$k_max %||% 0.7, cmp$diff, cmp$se,
          if (is.finite(cmp$ratio)) sprintf("%.2f", cmp$ratio) else "NA",
          if (!is.null(cmp$n_positive) && is.finite(cmp$n_positive)) sprintf("; %d of %d weeks positive, median %+.2f", cmp$n_positive, cmp$n_used, cmp$median_diff) else "")
}

# Every effort stream of one fit, each scored under the joint leave-out set of ALL of them,
# from a list of per-stream log-likelihood matrices and block labels:
#   streams = list(trailer = list(ll = <matrix>, block = <labels>), osp = list(...))
# Returns a named list of block tables. This is the one place the "other streams join the
# leave-out set" rule is applied, so the driver hook and the post-hoc runner cannot differ.
#
# THE JOINT TABLE (2026-09-27). When a fit has more than one effort stream the list also
# carries "joint": the week's log predictive density of EVERY held-out effort observation
# together, log E_w[ exp(sum of all streams' block log-likelihoods) ] under the same PSIS
# weights, one row per week that any stream observes. It is not the sum of the per-stream
# elpds: the two differ by log(1 + Cov(X, Y) / (E X . E Y)) with X and Y the streams' block
# likelihoods under the weights, a dependence term of either sign that runs to a few
# hundredths of a nat a week on the harness's conjugate check. It is the one number the
# question "did the term improve the interpolation of this week's boat effort?" has. The
# 2026-09-26 run scored the two boat streams separately, and the per-stream rule it was
# read under asked the trailer counts, the weaker instrument, to clear a bar on their own:
# on the same weeks the same effect scored 3.9 SE on the OSP counts and 1.6 SE on the
# trailer counts (Section 1y). The joint row is what that rule should have named. The
# per-stream tables are kept beside it: a term that helps one stream and hurts the other is
# still a finding, and the joint row is dominated by the stream with more information.
bss_block_cv_fit <- function(streams, k_max = 0.7) {
  streams <- streams[vapply(streams, function(s) !is.null(s$ll) && ncol(s$ll) > 0, logical(1))]
  out <- list()
  for (sn in names(streams)) {
    others <- setdiff(names(streams), sn)
    out[[sn]] <- bss_block_psis_loo(streams[[sn]]$ll, streams[[sn]]$block, k_max = k_max,
                                    weight_ll = if (length(others)) lapply(streams[others], `[[`, "ll") else NULL,
                                    weight_block = if (length(others)) lapply(streams[others], `[[`, "block") else NULL)
    out[[sn]]$leaveout_streams <- paste(names(streams), collapse = "+")
  }
  if (length(streams) > 1L) {
    nd <- vapply(streams, function(s) nrow(s$ll), numeric(1))
    if (length(unique(nd)) != 1L) stop("the streams' log-likelihood matrices do not share a draw count", call. = FALSE)
    ll_all <- do.call(cbind, lapply(streams, `[[`, "ll"))
    bl_all <- do.call(c, lapply(streams, function(s) as.character(s$block)))
    out[["joint"]] <- bss_block_psis_loo(ll_all, bl_all, k_max = k_max)
    out[["joint"]]$leaveout_streams <- paste(names(streams), collapse = "+")
    out[["joint"]] <- out[["joint"]][order(out[["joint"]]$block), , drop = FALSE]
    rownames(out[["joint"]]) <- NULL
  }
  out
}

# ---------------------------------------------------------------------------
# D54: the day-integrated likelihood. mu0 [draws x n] is each count's expectation WITHOUT the
# day factor, y the counts, eday each count's day-effect index, sig and size the draws of
# sigma_uE and r_E. Returns $fine and $coarse, [draws x days] (the second on every other node,
# the grid check), and $days, the day-effect indices in increasing order. The code is the
# prototype's (desk_advisory_block_cv_2026-11-19.R, integrate_day_loglik) line for line, so the
# two agree to the last bit (B85).
# ---------------------------------------------------------------------------
.bcv_row_lse <- function(A) { m <- A[cbind(seq_len(nrow(A)), max.col(A, ties.method = "first"))]; m + log(rowSums(exp(A - m))) }
bss_block_integrate_days <- function(mu0, y, eday, sig, size, h = 0.025, zmax = 8) {
  z <- seq(-zmax, zmax, by = h); K <- length(z); nd <- nrow(mu0)
  lw <- stats::dnorm(z, log = TRUE) + log(h)
  coarse <- seq(1L, K, by = 2L)
  shift <- outer(sig, z) - 0.5 * sig^2                      # [draws x K]: the log of the day factor at each node
  days <- sort(unique(eday))
  fine <- crs <- matrix(NA_real_, nd, length(days))
  for (j in seq_along(days)) {
    A <- matrix(0, nd, K)
    for (i in which(eday == days[j]))                       # size recycles down the rows: element [s, k] takes size[s]
      A <- A + stats::dnbinom(y[i], mu = exp(log(mu0[, i]) + shift), size = size, log = TRUE)
    A <- A + rep(lw, each = nd)
    fine[, j] <- .bcv_row_lse(A)
    crs[, j]  <- .bcv_row_lse(A[, coarse, drop = FALSE]) + log(2)
  }
  list(fine = fine, coarse = crs, days = days, n_nodes = K)
}

# The gear and trailer streams of a fit WITH the per-day effect, day-integrated. `draws` as
# rstan::extract() returns it (or a saved ppc draws list from B85 on): lambda_E_S, r_E, R_G and
# R_G_boat as the streams need, sigma_uE, uE_raw, and log_lik_<stream> when available (for check
# (i)). `event_dates` the date of each day index (length D; only the sampled days are read).
# Returns NULL when the fit carries no day effect, else a list with, per stream present,
#   $streams[[s]]  list(ll = [draws x days], block = weekly labels, n_counts = counts per day)
#   $coarse[[s]]   the same on the halved grid
#   $checks        one row per stream: counts, days, the check (i) gap, sigma_uE summary
bss_block_day_streams <- function(draws, stan_data, event_dates, h = 0.025, zmax = 8) {
  if (!identical(as.integer(stan_data$effort_day_effect %||% 0L), 1L) || as.integer(stan_data$n_Eday %||% 0L) == 0L) return(NULL)
  if (is.null(draws$sigma_uE) || is.null(draws$uE_raw) || !length(draws$sigma_uE))
    stop("the fit carries the per-day effect but the draws hold no sigma_uE / uE_raw", call. = FALSE)
  lam <- draws$lambda_E_S; nd <- dim(lam)[1]
  sig <- as.numeric(draws$sigma_uE); uu <- as.matrix(draws$uE_raw); size <- as.numeric(draws$r_E)
  ev <- as.Date(event_dates)
  out <- list(streams = list(), coarse = list(), checks = NULL)
  for (sn in c("gear", "trailer")) {
    n <- as.integer(if (sn == "gear") stan_data$Gear_n %||% 0L else stan_data$T_n %||% 0L)
    if (n == 0L) next
    y   <- as.integer(if (sn == "gear") stan_data$Gear_I else stan_data$T_I)[seq_len(n)]
    day <- as.integer(if (sn == "gear") stan_data$day_Gear else stan_data$day_T)[seq_len(n)]
    sec <- as.integer((if (sn == "gear") stan_data$section_Gear else stan_data$section_T) %||% rep(1L, n))[seq_len(n)]
    ed  <- as.integer(if (sn == "gear") stan_data$Eday_Gear else stan_data$Eday_T)[seq_len(n)]
    mult <- if (sn == "gear") as.numeric(draws$R_G) else 1 / as.numeric(draws$R_G_boat)
    if (!length(mult)) stop(sprintf("the %s stream needs %s draws", sn, if (sn == "gear") "R_G" else "R_G_boat"), call. = FALSE)
    mu0 <- vapply(seq_len(n), function(i) rowSums(matrix(lam[, sec[i], day[i], , drop = FALSE], nrow = nd)) * mult, numeric(nd))
    if (nd == 1L) mu0 <- matrix(mu0, nrow = 1L)
    # check (i): the conditional log-likelihood from mu0 and the fitted deviations against the fit's own
    llc <- draws[[paste0("log_lik_", sn)]]
    gap <- if (is.null(llc)) NA_real_ else {
      cond <- vapply(seq_len(n), function(i) stats::dnbinom(y[i], mu = mu0[, i] * exp(sig * uu[, ed[i]] - 0.5 * sig^2), size = size, log = TRUE), numeric(nd))
      max(abs(cond - as.matrix(llc)))
    }
    dd <- tapply(as.character(ev[day]), ed, function(v) length(unique(v)))
    if (any(dd != 1L)) stop(sprintf("a %s day-effect index spans more than one date", sn), call. = FALSE)
    int <- bss_block_integrate_days(mu0, y, ed, sig, size, h = h, zmax = zmax)
    day_date <- ev[day][match(int$days, ed)]
    wk <- bss_block_weeks(day_date)
    nc <- as.integer(table(factor(ed, levels = int$days)))
    out$streams[[sn]] <- list(ll = int$fine, block = wk, n_counts = nc)
    out$coarse[[sn]]  <- list(ll = int$coarse, block = wk, n_counts = nc)
    out$checks <- rbind(out$checks, data.frame(stream = sn, counts = n, days = length(int$days), grid_step = h, grid_nodes = int$n_nodes,
      conditional_loglik_max_abs_gap = gap, sigma_uE_mean = mean(sig), sigma_uE_q2.5 = unname(stats::quantile(sig, 0.025)),
      sigma_uE_q97.5 = unname(stats::quantile(sig, 0.975)), stringsAsFactors = FALSE))
  }
  out
}

# The in-driver convenience: block tables for every effort stream of one fit, from the
# stanfit's own log_lik_* generated quantities, written as loo_block_<stream>_<label>.csv.
# `days_ss` gives the dates (event_date by day index) the weekly blocks are cut from.
# Never stops a run: extraction is tryCatch-wrapped per stream and the scoring as a whole.
# `prefix` (2026-09-29, D29): the ladders write every RUNG's block CV as
# ladder_block_<stream>_<label>_<resolution>.csv, so the report's block-CV table, which reads
# the kept fit's loo_block_* files, is not mixed with rungs the run did not report.
write_block_cv_diagnostics <- function(fit, stan_data, days_ss, label, output_dir, k_max = 0.7,
                                       prefix = "loo_block", tol_loglik = 1e-6, tol_grid = 0.01) {
  if (is.null(fit) || !requireNamespace("loo", quietly = TRUE)) return(invisible(NULL))
  ev <- if (!is.null(days_ss) && "event_date" %in% names(days_ss)) as.Date(days_ss$event_date) else NULL
  if (is.null(ev)) return(invisible(NULL))
  spec <- list(gear    = list(par = "log_lik_gear",    n = stan_data$Gear_n %||% 0, days = stan_data$day_Gear),
               trailer = list(par = "log_lik_trailer", n = stan_data$T_n %||% 0,    days = stan_data$day_T),
               osp     = list(par = "log_lik_osp",     n = stan_data$OSP_n %||% 0,  days = stan_data$day_OSP))
  streams <- list()
  for (sn in names(spec)) {
    st <- spec[[sn]]
    if ((st$n %||% 0) == 0) next
    got <- tryCatch({
      ll <- as.matrix(rstan::extract(fit, pars = st$par)[[1]])
      if (nrow(ll) == 0 || ncol(ll) == 0) NULL else list(ll = ll, block = bss_block_weeks(ev[st$days[seq_len(ncol(ll))]]))
    }, error = function(e) { cat(sprintf("    [block CV] %s/%s not extracted: %s\n", label, sn, conditionMessage(e))); NULL })
    if (!is.null(got)) streams[[sn]] <- got
  }
  if (!length(streams)) return(invisible(NULL))
  # D54 (B85): on a fit with the per-day effect, the gear and trailer streams are scored with each
  # sampled day's deviation integrated out; the conditional tables are kept beside them.
  integ <- NULL
  if (identical(as.integer(stan_data$effort_day_effect %||% 0L), 1L) && as.integer(stan_data$n_Eday %||% 0L) > 0L &&
      any(c("gear", "trailer") %in% names(streams))) {
    integ <- tryCatch({
      pars <- c("lambda_E_S", "r_E", "sigma_uE", "uE_raw",
                if ("gear" %in% names(streams)) "R_G", if ("trailer" %in% names(streams)) "R_G_boat")
      dr <- rstan::extract(fit, pars = pars)
      for (sn in intersect(c("gear", "trailer"), names(streams))) dr[[paste0("log_lik_", sn)]] <- streams[[sn]]$ll
      bss_block_day_streams(dr, stan_data, ev)
    }, error = function(e) { cat(sprintf("    [block CV] %s: the day-integrated table was not built (%s); the conditional one is written\n", label, conditionMessage(e)))
      if (exists("bss_warn", mode = "function")) bss_warn("block CV", sprintf("%s: the day-integrated held-out-week table was not built (%s); loo_block_* is the conditional table, unusable on most weeks of a fit with the per-day effect (D54)", label, conditionMessage(e)))
      NULL })
  }
  tabs_cond <- tryCatch(bss_block_cv_fit(streams, k_max = k_max),
                        error = function(e) { cat(sprintf("    [block CV] %s skipped: %s\n", label, conditionMessage(e))); NULL })
  if (is.null(tabs_cond)) return(invisible(NULL))
  tabs <- tabs_cond; unit <- "count"
  if (!is.null(integ)) {
    ok <- tryCatch({
      s_int <- streams; s_crs <- streams
      for (sn in names(integ$streams)) { s_int[[sn]] <- integ$streams[[sn]]; s_crs[[sn]] <- integ$coarse[[sn]] }
      ti <- bss_block_cv_fit(s_int, k_max = k_max); th <- bss_block_cv_fit(s_crs, k_max = k_max)
      top <- if ("joint" %in% names(ti)) "joint" else names(integ$streams)[1]
      both <- merge(ti[[top]][, c("block", "elpd_block", "reliable")], th[[top]][, c("block", "elpd_block", "reliable")], by = "block")
      rel <- both$reliable.x & both$reliable.y
      gap_grid <- if (any(rel)) max(abs(both$elpd_block.x - both$elpd_block.y)[rel]) else NA_real_
      gap_ll <- suppressWarnings(max(integ$checks$conditional_loglik_max_abs_gap, na.rm = TRUE))
      chk <- cbind(fit = label, integ$checks, scored_row = top, weeks = nrow(ti[[top]]),
                   reliable_conditional = sum(tabs_cond[[top]]$reliable), reliable_integrated = sum(ti[[top]]$reliable),
                   median_k_conditional = stats::median(tabs_cond[[top]]$pareto_k), median_k_integrated = stats::median(ti[[top]]$pareto_k),
                   grid_halving_max_abs_elpd_gap = gap_grid, check_i = is.finite(gap_ll) && gap_ll <= tol_loglik,
                   check_ii = is.na(gap_grid) || gap_grid <= tol_grid, stringsAsFactors = FALSE)
      utils::write.csv(chk, file.path(output_dir, sprintf("blockcv_checks_%s.csv", label)), row.names = FALSE)
      if (!all(chk$check_i)) stop(sprintf("check (i) failed: the rebuilt conditional log-likelihood is %.2g from the fit's own", gap_ll))
      if (!all(chk$check_ii) && exists("bss_warn", mode = "function"))
        bss_warn("block CV", sprintf("%s: halving the quadrature grid moves a reliable week by %.3g nats (above %.2g); the day-integrated table is written, read it with that in mind", label, gap_grid, tol_grid))
      ti
    }, error = function(e) {
      cat(sprintf("    [block CV] %s: the day-integrated table is not used (%s)\n", label, conditionMessage(e)))
      if (exists("bss_warn", mode = "function")) bss_warn("block CV", sprintf("%s: the day-integrated held-out-week table is not used (%s); loo_block_* is the conditional table (D54)", label, conditionMessage(e)))
      NULL })
    if (!is.null(ok)) {
      tabs <- ok; unit <- "sampled day (gear and trailer rows: the day deviation integrated out, D54)"
      for (sn in names(tabs_cond)) {
        tab <- cbind(data_type = sn, tabs_cond[[sn]], scored_unit = "count (conditional on the day deviation)", stringsAsFactors = FALSE)
        utils::write.csv(tab, file.path(output_dir, sprintf("%scond_%s_%s.csv", prefix, sn, label)), row.names = FALSE)
      }
    }
  }
  for (sn in names(tabs)) {
    tab <- cbind(data_type = sn, tabs[[sn]], scored_unit = unit, stringsAsFactors = FALSE)
    utils::write.csv(tab, file.path(output_dir, sprintf("%s_%s_%s.csv", prefix, sn, label)), row.names = FALSE)
    cat(sprintf("    block CV %s/%s: %d weeks, %d reliable (k <= %.1f), held-out elpd %.1f; leave-out set %s%s\n",
                label, sn, nrow(tab), sum(tab$reliable), k_max, sum(tab$elpd_block[tab$reliable]), tab$leaveout_streams[1],
                if (identical(unit, "count")) "" else "; day-integrated"))
  }
  invisible(tabs)
}
