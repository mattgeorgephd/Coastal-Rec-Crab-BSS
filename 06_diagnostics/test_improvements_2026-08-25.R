# -----------------------------------------------------------------------------
# test_improvements_2026-08-25.R
#
# Standalone regression harness for the 2026-08-25 improvement batch. Deliberately
# dependency-light (dplyr / tibble / tidyr / stringr / purrr only, NO rstan, NO readxl),
# so it runs in seconds and can be used as a pre-flight check before committing to a
# multi-hour fit. It exercises the pure functions the batch added or changed, and asserts
# the SHIPPED defaults in run_config.R, which is the guard that stops a behaviour-changing
# toggle from drifting on unnoticed.
#
# Run from the repository root:   Rscript 06_diagnostics/test_improvements_2026-08-25.R
# Exits non-zero on any failure, so it can be wired into a pre-run check.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({library(dplyr); library(tibble); library(tidyr); library(stringr); library(purrr)})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Locate the repo root from wherever this was launched, so the harness tests the checkout
# it lives in rather than a hard-coded path.
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions")))
  stop("Run this from the repository root (or from 06_diagnostics/): 03_R_functions not found.")
here <- function(...) file.path(.root, ...)
setwd(.root)

for (f in c("bss_rng.R", "classify_day_type.R",   # 2026-09-28 (B46): bss_with_seed() and bss_weekday(), used by files below
            "bss_effort_spec.R","bss_ar_resolution.R","crab_fraction.R",
            "bss_opener_covariates.R","diagnose_incomplete_trips.R",
            "bss_stan_fit.R","save_run_diagnostics.R",
            "model_diagnostics.R","bss_model_adequacy.R",
            "annotate_decoupled_run.R",
            "bss_sampler_override.R",
            "pe_gear_ratio_frame.R",
            "bss_ar_rung_summary.R",
            "read_input_workbook.R")) source(file.path("03_R_functions", f))   # 2026-09-10: the workbook reader

ok <- 0; bad <- 0; skipped <- 0
# 2026-09-13: an ERROR while evaluating a condition is a FAIL, not the end of the run.
# `cond` is a promise, so an error raised while computing it surfaces here. Before this,
# any unexpected error in any one of a thousand conditions aborted the harness and threw
# away every result after it, which is the least informative possible outcome: you learn
# that something broke and nothing about the other 900 assertions. Found by
# negative-testing the census guard, where turning the missing-frame warning into a stop()
# killed the run at an exerciser 300 assertions before the assertion that was supposed to
# catch it.
chk <- function(nm, cond, extra="") {
  v <- tryCatch(isTRUE(cond), error = function(e) structure(FALSE, err = conditionMessage(e)))
  if (isTRUE(v)) { ok <<- ok+1; cat("PASS ", nm, extra, "\n") }
  else { bad <<- bad+1
         cat("FAIL ", nm, extra,
             if (!is.null(attr(v, "err"))) paste0("  [ERROR while evaluating: ",
                                                  substr(attr(v, "err"), 1, 220), "]") else "", "\n") } }

# D37 (2026-09-29): a check whose INPUT is absent is a SKIP, counted and printed apart. It
# used to be chk(name, TRUE, "skipped"), a PASS, or an if (dir.exists()) with no else, which
# ran nothing and said nothing; either way a checkout missing a committed run folder read as
# a clean harness. The summary line now reports the skips, and on a full checkout it is 0.
skp <- function(nm, why = "input absent") { skipped <<- skipped + 1; cat("SKIP ", nm, paste0("(", why, ")"), "\n") }

mkdays <- function(start, n) {
  d <- as.Date(start) + 0:(n-1)
  tibble(event_date = d, day_index = seq_len(n),
         week_index = as.integer(factor(format(d, "%Y-%W"), levels = unique(format(d, "%Y-%W")))),
         month_index = as.integer(factor(format(d, "%Y-%m"), levels = unique(format(d, "%Y-%m")))),
         day_type = ifelse(weekdays(d) %in% c("Saturday","Sunday"), "weekend", "weekday"))
}

P <- list(days_wkend=c("Saturday","Sunday"), shore_effort_unit="gear-deployments",
          tau_shore_prior_mu=1.7, tau_shore_prior_sigma=0.3,
          tau_boat_prior_mu=1.2, tau_boat_prior_sigma=0.3,
          use_crab_fraction=TRUE, crab_fraction_set=0.3, crab_fraction_prior_kappa=20,
          crab_fraction_min_obs=20, crab_fraction_strata="none",
          use_osp_boat_counts=TRUE, use_osp_crab_lower=TRUE,
          crab_fraction_osp_min_obs=20, crab_fraction_combo_share=0.15, crab_fraction_combo_kappa=8,
          opener_min_days=10)

days289 <- mkdays("2024-12-01", 289); days76 <- mkdays("2024-09-16", 76)

# ---------- 1. effort spec: I/E observation column follows the unit ----------
sp_s <- bss_effort_spec(TRUE, days289, P)
chk("shore deployments -> ie_trips", identical(sp_s$ie_obs_col, "ie_trips"), sp_s$ie_obs_col)
sp_s2 <- bss_effort_spec(TRUE, days289, modifyList(P, list(ie_shore_obs_unit="crabber_hours")))
chk("legacy override -> ie_crabber_hours", identical(sp_s2$ie_obs_col, "ie_crabber_hours"))
sp_h <- bss_effort_spec(TRUE, days289, modifyList(P, list(shore_effort_unit="crabber-hours")))
chk("shore crabber-hours -> ie_crabber_hours", identical(sp_h$ie_obs_col, "ie_crabber_hours"))
chk("boat -> ie_trips", identical(bss_effort_spec(FALSE, days289, P)$ie_obs_col, "ie_trips"))

# ---------- 2. AR ladder ----------
eff <- tibble(day_index = seq(1, 289, by = 2))
P0 <- modifyList(P, list(ar_escalate=FALSE, ar_max_resolution=list(shore=list(all_gear="daily", pot_closure="biweekly"), private_boat="monthly")))
chk("ladder off = single rung (shore all_gear)",
    identical(bss_ar_ladder(days289, eff, "shore", P0, gear_regime="all_gear"), "daily"))
chk("ladder off = single rung (boat capped monthly)",
    identical(bss_ar_ladder(days289, eff, "private_boat", P0, gear_regime="all_gear"), "monthly"))
P1 <- modifyList(P0, list(ar_escalate=TRUE, ar_escalate_ladder=c("daily","weekly","biweekly","monthly"),
                          ar_escalate_max_attempts=4, ar_escalate_respect_cap=FALSE))
L1 <- bss_ar_ladder(days289, eff, "private_boat", P1, gear_regime="all_gear")
chk("ladder on starts at daily, ignores cap", identical(L1, c("daily","weekly","biweekly","monthly")), paste(L1, collapse="->"))
L2 <- bss_ar_ladder(days289, eff, "private_boat", modifyList(P1, list(ar_escalate_respect_cap=TRUE)), gear_regime="all_gear")
chk("respect_cap starts at the cap", identical(L2, "monthly"), paste(L2, collapse="->"))
L3 <- bss_ar_ladder(days76, tibble(day_index=seq(1,76,by=2)), "shore", P1, gear_regime="pot_closure")
chk("76-day window keeps 4 distinct rungs", length(L3)==4 && L3[1]=="daily", paste(L3, collapse="->"))
L4 <- bss_ar_ladder(mkdays("2024-12-01", 20), tibble(day_index=1:10), "shore", P1)
chk("short window drops degenerate rungs", length(L4) >= 1 && !any(duplicated(L4)), paste(L4, collapse="->"))
L5 <- bss_ar_ladder(days289, eff, "private_boat", modifyList(P1, list(ar_force=list(private_boat="weekly"))))
chk("ar_force outranks the ladder", identical(L5, "weekly"))

# ---------- 3. crab fraction: exact back-compat, then the OSP lower bound ----------
Pn <- modifyList(P, list(use_osp_crab_lower=FALSE))
cf0 <- crab_fraction_stan_data(FALSE, days289, Pn, quiet=TRUE)
chk("no OSP: osp_crab_lower = 0", cf0$osp_crab_lower == 0L)
chk("no OSP: theta prior = the f prior", isTRUE(all.equal(as.numeric(cf0$crab_fraction_alpha0), 0.3*20)) &&
                                          isTRUE(all.equal(as.numeric(cf0$crab_fraction_beta0), 0.7*20)))
chk("shore is pinned to f = 1", crab_fraction_stan_data(TRUE, days289, P, quiet=TRUE)$apply_crab_fraction == 0L)
chk("point f with no data = set value",
    isTRUE(all.equal(unique(crab_fraction_point_day(TRUE, days289, Pn)), 0.3)))

osp_rows <- tibble(event_date = days289$event_date[1:200], osp_total = 50, osp_crab_only = 10)  # 20% crab-only
# NOTE: modifyList() recurses into list-like values, and a tibble IS a list, so
# modifyList(params, list(osp_crab_rows = <tibble>)) MERGES columns instead of replacing
# the frame. params carries two data-frame-valued keys (crab_fraction_rows,
# osp_crab_rows); set them by assignment, never through modifyList.
setp <- function(p, ...) { v <- list(...); for (nm in names(v)) p[[nm]] <- v[[nm]]; p }
Po <- setp(P, osp_crab_rows = osp_rows)
cf1 <- crab_fraction_stan_data(FALSE, days289, Po, quiet=TRUE)
chk("OSP present: bound switched on", cf1$osp_crab_lower == 1L)
chk("OSP present: counts aggregated", cf1$osp_f_n_total[1] == 10000 && cf1$osp_f_n_crab[1] == 2000,
    paste(cf1$osp_f_n_total[1], cf1$osp_f_n_crab[1]))
chk("OSP present: theta prior switches to combo",
    isTRUE(all.equal(as.numeric(cf1$crab_fraction_alpha0), 0.15*8)))
pf <- unique(crab_fraction_point_day(TRUE, days289, Po))
expected <- (1+2000)/(2+10000); expected <- expected + (1-expected)*0.15
chk("point f = lower + (1-lower)*theta", isTRUE(all.equal(pf, expected, tolerance=1e-6)),
    sprintf("got %.4f expected %.4f", pf, expected))
chk("point f is above the OSP lower bound", pf > 2000/10000)

# thin OSP stratum must not bind
osp_thin <- tibble(event_date = days289$event_date[1:2], osp_total = 5, osp_crab_only = 1)
cf2 <- crab_fraction_stan_data(FALSE, days289, setp(P, osp_crab_rows=osp_thin), quiet=TRUE)
chk("thin OSP stratum does not bind", cf2$osp_crab_lower == 0L && cf2$osp_f_n_total[1] == 0)

# egress + OSP together
eg_rows <- tibble(event_date = days289$event_date[1:50], boats_crabbing = 3, boats_total = 10)
cf3 <- crab_fraction_stan_data(FALSE, days289, setp(Po, crab_fraction_rows=eg_rows), quiet=TRUE)
chk("egress binomial still supplied", cf3$crab_fraction_n_total[1] == 500 && cf3$crab_fraction_n_crab[1] == 150)
pf3 <- unique(crab_fraction_point_day(TRUE, days289, setp(Po, crab_fraction_rows=eg_rows)))
chk("egress pulls point f toward 0.30", abs(pf3 - 0.30) < 0.03, sprintf("%.4f", pf3))

# per-month strata + opener strata
cfm <- crab_fraction_stan_data(FALSE, days289, setp(Po, crab_fraction_strata="month"), quiet=TRUE)
chk("month strata produce K = 10", cfm$n_f_strata == 10, cfm$n_f_strata)
cfo <- crab_fraction_stan_data(FALSE, days289, setp(Po,
        crab_fraction_strata="opener", opener_f_dates=days289$event_date[100:200]), quiet=TRUE)
chk("opener strata produce K = 2", cfo$n_f_strata == 2, cfo$n_f_strata)

# pinned f
cfp <- crab_fraction_stan_data(FALSE, days289, setp(Po, crab_fraction_fixed=0.4), quiet=TRUE)
chk("pinned f disables the bound", cfp$crab_fraction_estimate == 0L && cfp$osp_crab_lower == 0L)
chk("pinned point f = 0.4", isTRUE(all.equal(unique(crab_fraction_point_day(TRUE, days289, setp(Po, crab_fraction_fixed=0.4))), 0.4)))

# ---------- 4. opener design matrix ----------
flags <- tibble(event_date = days289$event_date,
                ma2_halibut_open = days289$event_date %in% days289$event_date[100:180],
                razor_nearby_dig = days289$event_date %in% days289$event_date[1:5])
o0 <- opener_design_matrix(days289, character(0), flags, P)
chk("K_open = 0 when nothing selected", o0$K_open == 0 && length(o0$X_open) == 0)
o1 <- opener_design_matrix(days289, "ma2_halibut_open", flags, P)
chk("selected opener -> 1 column", o1$K_open == 1 && length(o1$X_open) == 289)
chk("flat vector matches the flags", sum(o1$X_open) == 81, sum(o1$X_open))
o2 <- opener_design_matrix(days289, "razor_nearby_dig", flags, P)
chk("near-constant column dropped", o2$K_open == 0 && length(o2$dropped) == 1, o2$dropped)
o3 <- opener_design_matrix(days289, "ma2_halibut_open", flags, P, extra="ma2_halibut_open")
chk("duplicate column de-duplicated", o3$K_open == 1)
o4 <- opener_design_matrix(days289, character(0), flags, P, extra="nonexistent_flag")
chk("absent flag reported not fatal", o4$K_open == 0 && length(o4$dropped) == 1)

# ---------- 5. opener screen ----------
spill <- list(adjusted = tibble(
  Series = c(rep("Shore gear (effort)",4), rep("Boat trailers (effort)",4), rep("Shore CPUE",4)),
  Opener = rep(c("MA2 salmon open","MA2 halibut open","MA2 bottomfish open","Razor dig (nearby beaches)"), 3),
  adj_estimate = c(1,2,3,18, 5,31.7,2,1, 0,0,0,0),
  adj_p = c(0.77,0.40,0.55,0.045, 0.62,3.6e-6,0.30,0.51, .9,.9,.9,.9),
  note = ""))
Ps <- modifyList(P, list(opener_covariate_mode="auto", opener_auto_p=0.05, opener_auto_p_adjust="BH",
                         opener_candidates_shore="razor_nearby_dig",
                         opener_candidates_boat=c("ma2_halibut_open","ma2_salmon_open","ma2_bottomfish_open")))
sel <- opener_select(spill, Ps)
chk("BH screen keeps halibut->boat", identical(sel$private_boat, "ma2_halibut_open"), paste(sel$private_boat, collapse=","))
chk("BH screen drops razor->shore (Run 3 agreement)", length(sel$shore) == 0)
sel_raw <- opener_select(spill, modifyList(Ps, list(opener_auto_p_adjust="none")))
chk("unadjusted screen WOULD have kept razor", identical(sel_raw$shore, "razor_nearby_dig"))
chk("mode off selects nothing", length(opener_select(spill, modifyList(Ps, list(opener_covariate_mode="off")))$private_boat) == 0)
selm <- opener_select(spill, modifyList(Ps, list(opener_covariate_mode="manual", opener_manual_boat="ma2_salmon_open")))
chk("manual mode honours the list", identical(selm$private_boat, "ma2_salmon_open"))
chk("auto with no diagnostic degrades safely",
    length(opener_select(NULL, Ps)$private_boat) == 0 && length(opener_select(NULL, Ps)$note) > 0)

# ---------- 6. incomplete-trip arms ----------
set.seed(1)
iv <- tibble(trip_status = rep(c("Complete","Incomplete", NA), c(60,40,10)),
             number_of_gear = c(rep(4,60), rep(6,40), rep(4,10)),
             angler_count = 2,
             Dungeness_Kept = c(rep(8,60), rep(4,40), rep(8,10)))
fr_ex <- incomplete_trip_arm_frames(iv, "exclude", "number_of_gear")
fr_go <- incomplete_trip_arm_frames(iv, "gear_only", "number_of_gear")
fr_im <- incomplete_trip_arm_frames(iv, "impute_mean_cpue", "number_of_gear")
fr_kp <- incomplete_trip_arm_frames(iv, "keep", "number_of_gear")
chk("exclude drops incomplete from both frames", nrow(fr_ex$cpue)==70 && nrow(fr_ex$gear)==70)
chk("gear_only keeps gear, drops catch", nrow(fr_go$cpue)==70 && nrow(fr_go$gear)==110)
chk("keep retains everything", nrow(fr_kp$cpue)==110 && nrow(fr_kp$gear)==110)
ros <- function(d) sum(d$Dungeness_Kept)/sum(d$number_of_gear)
chk("impute_mean_cpue is a pooled no-op on CPUE",
    isTRUE(all.equal(ros(fr_im$cpue), ros(fr_ex$cpue))),
    sprintf("%.6f vs %.6f", ros(fr_im$cpue), ros(fr_ex$cpue)))
chk("gear_only raises the gear ratio here (length-bias signature)",
    mean(fr_go$gear$number_of_gear) > mean(fr_ex$gear$number_of_gear))

# ---------- 7. review fixes (2026-08-25, post-review) ----------
# 7a. the OSP bound now feeds DAILY beta-binomial rows, not one summed binomial
Po2 <- setp(P, osp_crab_rows = osp_rows, use_osp_crab_lower = TRUE)
cfd <- crab_fraction_stan_data(FALSE, days289, Po2, quiet=TRUE)
chk("daily OSP rows are emitted", cfd$OSPF_n == 200, cfd$OSPF_n)
chk("daily rows carry a valid stratum index",
    all(cfd$osp_f_stratum >= 1) && all(cfd$osp_f_stratum <= cfd$n_f_strata))
chk("daily totals/crab match the source rows",
    sum(cfd$osp_f_total) == 10000 && sum(cfd$osp_f_crab) == 2000)
chk("kappa prior passed through", cfd$osp_f_kappa_prior_mu == 20)
chk("bound off -> no daily rows",
    crab_fraction_stan_data(FALSE, days289, setp(Po2, use_osp_crab_lower=FALSE), quiet=TRUE)$OSPF_n == 0)
cfp2 <- crab_fraction_stan_data(FALSE, days289, setp(Po2, crab_fraction_fixed=0.4), quiet=TRUE)
chk("pinned f emits no daily rows", cfp2$OSPF_n == 0 && cfp2$osp_crab_lower == 0L)
# a stratum below the minimum must contribute no daily rows either
osp_mixed <- dplyr::bind_rows(osp_rows[1:100,], tibble(event_date = days289$event_date[250:251], osp_total = 3, osp_crab_only = 1))
cfm2 <- crab_fraction_stan_data(FALSE, days289, setp(Po2, osp_crab_rows = osp_mixed, crab_fraction_strata="month"), quiet=TRUE)
chk("thin strata contribute no daily rows",
    all(cfm2$osp_f_n_total[cfm2$osp_f_stratum] >= 20))
# 7b. min_obs = 0 must not hand a data-free stratum the combo prior (Stan pins f_lower=0)
cf0b <- crab_fraction_stan_data(FALSE, days289, setp(P,
          use_osp_crab_lower=TRUE, crab_fraction_osp_min_obs=0,
          osp_crab_rows = tibble(event_date=as.Date(character()), osp_total=numeric(), osp_crab_only=numeric())), quiet=TRUE)
chk("min_obs=0 with no OSP data keeps the ordinary f prior",
    isTRUE(all.equal(as.numeric(cf0b$crab_fraction_alpha0), 0.3*20)) && cf0b$osp_crab_lower == 0L)
# 7c. both f streams now share the fit's date window
eg_out <- tibble(event_date = days76$event_date, boats_crabbing = 5, boats_total = 20)
cfw <- crab_fraction_stan_data(FALSE, days289, setp(P, crab_fraction_rows = eg_out), quiet=TRUE)
chk("egress rows outside the fit window are excluded", cfw$crab_fraction_n_total[1] == 0)
cfw2 <- crab_fraction_stan_data(FALSE, days289, setp(P,
          crab_fraction_rows = eg_out, crab_fraction_restrict_to_fit = FALSE), quiet=TRUE)
chk("restrict_to_fit = FALSE restores season-wide pooling", cfw2$crab_fraction_n_total[1] > 0)
# 7d. opener screen: manual mode respects the candidate list; family size is pinned
selm2 <- opener_select(spill, modifyList(Ps, list(opener_covariate_mode="manual",
                                                  opener_manual_boat="razor_nearby_dig")))
chk("manual non-candidate is rejected", length(selm2$private_boat) == 0)
spill_na <- spill; spill_na$adjusted$adj_p[c(1,3)] <- NA_real_
sel_na <- opener_select(spill_na, Ps)
chk("NA tests do not shrink the multiplicity family",
    identical(sel_na$private_boat, "ma2_halibut_open") && length(sel_na$shore) == 0)
# 7e. shore gear ratio uses one row set for numerator and denominator
iv2 <- iv; iv2$angler_count[1:10] <- NA
fr2 <- incomplete_trip_arm_frames(iv2, "keep", "number_of_gear")
g2 <- suppressWarnings(as.numeric(fr2$gear$number_of_gear)); a2 <- suppressWarnings(as.numeric(fr2$gear$angler_count))
k2 <- is.finite(g2) & g2 > 0 & is.finite(a2) & a2 > 0
chk("gear ratio uses matched rows", isTRUE(all.equal(sum(g2[k2])/sum(a2[k2]), sum(g2[k2])/sum(a2[k2]))) && sum(k2) == 100, sum(k2))

# ---------- 8. the SHIPPED defaults in run_config.R ----------
# Regression guard: everything behaviour-changing in this batch must ship OFF except the
# three items that were explicitly requested as changes. If one of these flips, a routine
# run silently stops being comparable to the last one.
local({
  e <- new.env(); sys.source(file.path("run_config.R"), envir = e); rc <- e$run_config
  chk("shipped: weekend = Sat/Sun", identical(rc$days_wkend, c("Saturday","Sunday")))
  chk("shipped: bss_min_interviews = 15", identical(rc$bss_min_interviews, 15))
  chk("shipped: ie_shore_obs_unit = auto (the I/E unit FIX is on)", identical(rc$ie_shore_obs_unit, "auto"))
  chk("shipped: ar_escalate OFF", identical(rc$ar_escalate, FALSE))
  chk("shipped: opener_covariate_mode OFF", identical(rc$opener_covariate_mode, "off"))
  chk("shipped: razor_dig_mode no", identical(rc$razor_dig_mode, "no"))
  # 2026-10-04 (A33): ON, now that OSP delivered the crab-only column (B64); it needs the OSP
  # stream, the dynamic f and the combo walk (the Stan rejects osp_crab_lower without dynamic c)
  chk("shipped: use_osp_crab_lower ON (A33), with use_osp_boat_counts and the dynamic f it needs",
      identical(rc$use_osp_crab_lower, TRUE) && identical(rc$use_osp_boat_counts, TRUE) && identical(rc$crab_fraction_dynamic, TRUE))
  # 2026-09-10: the three unsampled-cell levers ship at their new values. "zero" is no
  # longer a neutral default; it is an assumption that 15.6% of the shore's days and 16.6%
  # of the boat's had no fishing, and its error is a bias no SE can carry.
  chk("shipped: pe_empty_effort_stratum = local_day_type (2026-09-10)", identical(rc$pe_empty_effort_stratum, "local_day_type"))
  chk("shipped: pe_empty_stratum = local, so the CPUE fill matches the effort fill's scale", identical(rc$pe_empty_stratum, "local"))
  chk("shipped: pe_variance = impute_aware, so an imputed or singleton cell is not free", identical(rc$pe_variance, "impute_aware"))
  chk("shipped: filter_incomplete_trips still TRUE (diagnostic only)", identical(rc$filter_incomplete_trips, TRUE))
  chk("shipped: crab_fraction_strata = month (review item 1, 2026-09-08)", identical(rc$crab_fraction_strata, "month"))
  chk("shipped: crab_fraction_dynamic = TRUE (review item 1B, 2026-09-08)", identical(rc$crab_fraction_dynamic, TRUE))
})


# ---------- 9. STAN DATA CONTRACT ----------
# THE TEST THAT WAS MISSING. The 2026-08-25 patch declared five new variables in both
# .stan data blocks (OSPF_n, osp_f_stratum, osp_f_total, osp_f_crab,
# osp_f_kappa_prior_mu) and forwarded three of the eight new crab_fraction_stan_data()
# fields out of the prep functions. Stan then failed at data initialization on every
# fit, rstan returned an EMPTY stanfit instead of raising, and the run died 300 lines
# downstream in apply() with "dim(X) must have a positive length". Every fit of every
# rung of the validation ladder was lost to it.
#
# Three layers of guard, cheapest first:
#   9a  the parser itself is sane (a silently-empty parse would make 9b/9c vacuous)
#   9b  every declared Stan data variable is at least MENTIONED in the prep that builds
#       that model's data list  --  a static check, no data files, runs in milliseconds
#   9c  every crab-fraction field the Stan models declare is returned by EVERY return
#       path of crab_fraction_stan_data()  --  exact, exercised on real return values
# The exact per-fit check lives in bss_assert_stan_data(), which bss_stan_fit() runs
# before every sampler call.
local({
  models <- c(pooled        = file.path("02_stan_models", "crab_bss_pooled.stan"),
              gear_resolved = file.path("02_stan_models", "crab_bss_gear_resolved.stan"))
  preps  <- c(pooled        = file.path("03_R_functions", "prep_bss_crab_pooled.R"),
              gear_resolved = file.path("03_R_functions", "prep_bss_crab_gear.R"))

  # 9a. parser sanity
  np <- bss_stan_data_names(models[["pooled"]])
  chk("stan data parser returns a plausible variable count", length(np) > 60, length(np))
  chk("stan data parser finds old-style array declarations", all(c("period","O") %in% np))
  chk("stan data parser finds the improvement-8 block",
      all(c("OSPF_n","osp_f_stratum","osp_f_total","osp_f_crab","osp_f_kappa_prior_mu") %in% np))

  # 9b. every declared variable is mentioned in its prep
  for (m in names(models)) {
    need <- bss_stan_data_names(models[[m]])
    src  <- paste(readLines(preps[[m]], warn = FALSE), collapse = "\n")
    miss <- need[!vapply(need, function(v)
      grepl(paste0("(^|[^A-Za-z0-9_.])", v, "([^A-Za-z0-9_]|$)"), src, perl = TRUE),
      logical(1))]
    chk(sprintf("%s: every Stan data variable is built in %s", m, basename(preps[[m]])),
        length(miss) == 0, if (length(miss)) paste("MISSING:", paste(miss, collapse=", ")) else "")
  }

  # 9c. every crab-fraction field, on every return path of crab_fraction_stan_data()
  cf_need <- intersect(bss_stan_data_names(models[["pooled"]]),
                       bss_stan_data_names(models[["gear_resolved"]]))
  # review item 1B (2026-09-08): the dynamic-f fields (f_walk_*, f_level_*, CFI_n, cfi_*,
  # combo_*) join the contract; they are declared unconditionally in both models.
  cf_need <- grep(paste0("^(apply_crab_fraction|crab_fraction_.+|n_f_strata|f_stratum|osp_crab_lower|osp_f_.+|OSPF_n|",
                         "f_walk_.+|f_level_.+|CFI_n|cfi_.+|combo_dynamic|c_level_.+|c_walk_.+|CFC_n|cfc_.+)$"),
                  cf_need, value = TRUE)
  chk("crab-fraction contract is non-trivial", length(cf_need) >= 38, length(cf_need))

  paths <- list(
    `shore / feature off` = crab_fraction_stan_data(TRUE,  days289, P,  quiet = TRUE),
    `boat, f pinned`      = crab_fraction_stan_data(FALSE, days289,
                              modifyList(P, list(crab_fraction_fixed = 0.35)), quiet = TRUE),
    `boat, f estimated`   = crab_fraction_stan_data(FALSE, days289, P,  quiet = TRUE))
  for (nm in names(paths)) {
    miss <- setdiff(cf_need, names(paths[[nm]]))
    chk(sprintf("crab_fraction_stan_data covers the contract (%s)", nm),
        length(miss) == 0, if (length(miss)) paste("MISSING:", paste(miss, collapse=", ")) else "")
  }

  # 9d. the guards themselves fire
  good <- as.list(setNames(rep(list(0L), length(np)), np))
  chk("bss_assert_stan_data passes a complete list",
      isTRUE(tryCatch(bss_assert_stan_data(good, models[["pooled"]], "unit"),
                      error = function(e) conditionMessage(e))))
  chk("bss_assert_stan_data rejects a missing variable",
      grepl("OSPF_n", tryCatch({bss_assert_stan_data(good[setdiff(np, "OSPF_n")],
                                                     models[["pooled"]], "unit"); ""},
                               error = function(e) conditionMessage(e)), fixed = TRUE))
  bad_na <- good; bad_na$L_data <- c(1, NA, 3)
  chk("bss_assert_stan_data rejects a non-finite value",
      grepl("L_data", tryCatch({bss_assert_stan_data(bad_na, models[["pooled"]], "unit"); ""},
                               error = function(e) conditionMessage(e)), fixed = TRUE))
  chk("bss_assert_fit_usable rejects a non-stanfit",
      grepl("EMPTY stanfit", tryCatch({bss_assert_fit_usable(NULL, "unit"); ""},
                                      error = function(e) conditionMessage(e)), fixed = TRUE))
})


# ---------- 10. THE 2026-08-27 POST-LADDER FIXES ----------
# Each block below is a defect the 2026-08-26 validation ladder exposed. The ladder itself
# is the regression test for the model; these are the regression tests for the things the
# ladder could not see because they live in the reporting layer.
local({

  # 10a. RETIRED 2026-09-28 (B44): .srd_monthly_share(), the second copy of the PE's
  # daily-effort formula whose shore day-length weighting these assertions pinned, was
  # removed with pe_monthly_effort_share(). monthly_pe_vs_bss.csv now takes its PE column
  # from pe_monthly_split(), which reads the PE's own strata (section 81), so no per-day
  # factor (day length, f) has a second formula to drift from. Assert the copy is gone.
  chk("monthly share: the retired second copy (.srd_monthly_share) is gone",
      !exists(".srd_monthly_share", mode = "function") &&
      !any(grepl("^\\.srd_monthly_share <- function", readLines("03_R_functions/save_run_diagnostics.R", warn = FALSE))))

  # 10b. The empty-effort-stratum report is a FILE, not a console line. The pooled driver's
  # PE chunk is results='hide', so the cat()-only version reached nothing on that track.
  td <- file.path(tempdir(), paste0("pe_empty_", as.integer(runif(1, 1e5, 1e6))))
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  pe_fake <- list(
    shore_all_gear = list(n_empty_effort_strata = 0L, n_empty_effort_days = 0L,
                          n_single_effort_strata = 30L, n_single_effort_days = 60L,
                          n_effort_strata_total = 84L, n_calendar_days = 289L,
                          pe_empty_effort_fill = "local_day_type", pe_variance = "impute_aware",
                          effort_total = 30000, effort_se = 1200, effort_se_sampled_only = 400,
                          pe_imputed_effort = 0, pe_zeroed_effort_bias = 0),
    private_boat_ring_net_only = list(n_empty_effort_strata = 9L, n_empty_effort_days = 9L,
                          n_single_effort_strata = 5L, n_single_effort_days = 8L,
                          n_effort_strata_total = 22L, n_calendar_days = 76L,
                          pe_empty_effort_fill = "zero", pe_variance = "impute_aware",
                          effort_total = 400, effort_se = 90, effort_se_sampled_only = 29,
                          pe_imputed_effort = 0, pe_zeroed_effort_bias = 63),
    comm_charter = list(effort_total = 283))
  rep_df <- write_pe_empty_stratum_report(pe_fake, td, list(pe_empty_stratum = "local", pe_variance = "impute_aware"))
  chk("empty-stratum report writes a file", file.exists(file.path(td, "pe_empty_effort_strata.csv")))
  chk("empty-stratum report excludes the census component",
      !"comm_charter" %in% rep_df$component && nrow(rep_df) == 2)
  chk("empty-stratum report computes the zeroed-day fraction",
      isTRUE(all.equal(rep_df$empty_day_fraction[rep_df$component == "private_boat_ring_net_only"],
                       9/76)))
  chk("empty-stratum report raises the >5%-at-zero flag",
      isTRUE(rep_df$exceeds_5pct_at_zero[rep_df$component == "private_boat_ring_net_only"]) &&
      isFALSE(rep_df$exceeds_5pct_at_zero[rep_df$component == "shore_all_gear"]))
  # 2026-09-10: the report now carries the SINGLETON cells, both SEs and the zeroing bias.
  # The singleton count was never reported anywhere before, and it is the LARGER half of
  # the variance understatement (SE 410 -> 1,099 on shore all-gear from singletons alone).
  chk("empty-stratum report carries the singleton cells and the thin-day share",
      identical(rep_df$n_single_strata[rep_df$component == "shore_all_gear"], 30L) &&
      isTRUE(all.equal(rep_df$thin_day_fraction[rep_df$component == "shore_all_gear"], 60/289)))
  chk("empty-stratum report carries BOTH SEs, so the historical number stays reachable",
      isTRUE(all.equal(rep_df$effort_se[rep_df$component == "shore_all_gear"], 1200)) &&
      isTRUE(all.equal(rep_df$effort_se_sampled_only[rep_df$component == "shore_all_gear"], 400)) &&
      isTRUE(all.equal(rep_df$effort_cv[rep_df$component == "shore_all_gear"], 1200/30000)))
  chk("empty-stratum report carries the zeroing BIAS for a component left at 'zero'",
      isTRUE(all.equal(rep_df$zeroed_effort_bias[rep_df$component == "private_boat_ring_net_only"], 63)) &&
      isTRUE(all.equal(rep_df$zeroed_effort_bias[rep_df$component == "shore_all_gear"], 0)))
  chk("empty-stratum report records all three levers per component",
      identical(rep_df$effort_fill[rep_df$component == "shore_all_gear"], "local_day_type") &&
      identical(unique(rep_df$cpue_fill), "local") && identical(unique(rep_df$variance), "impute_aware"))
  unlink(td, recursive = TRUE)

  # 10c. I/E observation PROVENANCE. bss_effort_spec() must name the column the shore
  # likelihood consumes, and the two settings must not resolve to the same column -- that is
  # the whole content of rung 2, and no run output recorded it before this batch.
  Pd <- list(shore_effort_unit = "gear-deployments", tau_shore_prior_mu = 1.7,
             tau_shore_prior_sigma = 0.3, tau_boat_prior_mu = 1.2, tau_boat_prior_sigma = 0.3)
  d6 <- data.frame(event_date = as.Date("2024-12-01") + 0:5, day_type = "weekday")
  a_col <- bss_effort_spec(TRUE, d6, modifyList(Pd, list(ie_shore_obs_unit = "auto")))$ie_obs_col
  l_col <- bss_effort_spec(TRUE, d6, modifyList(Pd, list(ie_shore_obs_unit = "crabber_hours")))$ie_obs_col
  chk("I/E provenance: auto and legacy name DIFFERENT columns",
      identical(a_col, "ie_trips") && identical(l_col, "ie_crabber_hours"))
  chk("I/E provenance: the boat spec names a column too",
      nzchar(bss_effort_spec(FALSE, d6, Pd)$ie_obs_col %||% ""))
})

# ---------- 11. per-estimator production arm ----------
# The single "exclude" label was honest for the BSS and wrong for the boat PE: run_pe_*()
# takes the boat gear-per-group from the UNFILTERED interview set, so the boat PE already
# behaves like gear_only. The 2026-08-26 ladder made that concrete -- the shipped boat PE
# (3,565.75 effort / 10,940.36 catch) equals this table's gear_only arm, not its exclude arm.
# A table that labels both "exclude" hides its own headline.
local({
  a_shore <- incomplete_trip_production_arm(TRUE,  TRUE)
  a_boat  <- incomplete_trip_production_arm(FALSE, TRUE)
  a_off   <- incomplete_trip_production_arm(FALSE, FALSE)
  chk("production arm: BSS is 'exclude' for both populations",
      identical(unname(a_shore[["bss"]]), "exclude") && identical(unname(a_boat[["bss"]]), "exclude"))
  chk("production arm: the SHORE PE matches the BSS",
      identical(unname(a_shore[["pe"]]), "exclude"))
  # 2026-09-28 (B46): under the shipped pe_gear_ratio_arm = "match_bss" the boat PE takes its
  # gear ratio from the filtered frame, so it runs "exclude"; "gear_only" is the old asymmetry.
  chk("production arm: the BOAT PE follows pe_gear_ratio_arm (exclude under match_bss, gear_only under gear_only)",
      identical(unname(a_boat[["pe"]]), "exclude") &&
      identical(unname(incomplete_trip_production_arm(FALSE, TRUE, "gear_only")[["pe"]]), "gear_only"))
  chk("production arm: with the filter off, every estimator is 'keep'",
      identical(unname(a_off[["bss"]]), "keep") && identical(unname(a_off[["pe"]]), "keep") &&
      identical(unname(incomplete_trip_production_arm(TRUE, FALSE)[["pe"]]), "keep"))
})

# ---------- 12. decoupled-parameter flag ----------
# structural_params_*.csv puts prior-only parameters in the same columns as estimated ones.
# The 2026-08-26 ladder's worked example: under the production osp_scale_is_tau = TRUE the OSP
# mean uses L, so kappa_OSP is inert and reports its lognormal(log 3, 0.3) prior EXACTLY
# (median 3.008, 95% 1.63-5.40) -- which reads as "the model measured the turnover at 3.0".
# It did not; it was told 3.0. Every such parameter must now arrive labelled.
local({
  shore <- list(IE_n = 0L, OSP_n = 0L, T_n = 0L, osp_scale_is_tau = 1L, K_open = 0L,
                apply_crab_fraction = 0L, crab_fraction_estimate = 0L, osp_crab_lower = 0L,
                estimate_cpue_density = 0L, w = c(1,0,1), holiday = c(0,0,0))
  boat  <- list(IE_n = 0L, OSP_n = 148L, T_n = 60L, osp_scale_is_tau = 1L, K_open = 0L,
                apply_crab_fraction = 1L, crab_fraction_estimate = 1L, osp_crab_lower = 0L,
                estimate_cpue_density = 0L, w = c(1,0,1), holiday = c(1,0,0))
  pars <- c("B1","B2","B2_C","gamma_C","sigma_IE","R_G","R_G_boat","kappa_OSP","r_OSP",
            "f_crab[1]","f_lower[1]","r_E")
  rs <- bss_decoupled_reasons(pars, shore); names(rs) <- pars
  rb <- bss_decoupled_reasons(pars, boat);  names(rb) <- pars

  chk("decoupled: shore sigma_IE flagged when IE_n = 0", !is.na(rs[["sigma_IE"]]))
  chk("decoupled: shore OSP parameters flagged when OSP_n = 0",
      !is.na(rs[["kappa_OSP"]]) && !is.na(rs[["r_OSP"]]) && !is.na(rs[["R_G_boat"]]))
  chk("decoupled: THE kappa_OSP CASE -- flagged on the boat under osp_scale_is_tau = 1",
      !is.na(rb[["kappa_OSP"]]) && grepl("osp_scale_is_tau", rb[["kappa_OSP"]], fixed = TRUE))
  chk("decoupled: r_OSP NOT flagged on a boat fit that has OSP data", is.na(rb[["r_OSP"]]))
  chk("decoupled: f_lower flagged while use_osp_crab_lower is off", !is.na(rb[["f_lower[1]"]]))
  chk("decoupled: f_crab flagged on shore, not on an estimating boat fit",
      !is.na(rs[["f_crab[1]"]]) && is.na(rb[["f_crab[1]"]]))
  chk("decoupled: holiday terms flagged only when the window has no holiday",
      !is.na(rs[["B2"]]) && !is.na(rs[["B2_C"]]) && is.na(rb[["B2"]]))
  chk("decoupled: genuinely estimated parameters are NOT flagged",
      is.na(rs[["B1"]]) && is.na(rs[["R_G"]]) && is.na(rs[["r_E"]]) && is.na(rb[["R_G_boat"]]))
  chk("decoupled: Stan indices are stripped before matching",
      identical(bss_decoupled_reasons("f_lower[3]", boat), bss_decoupled_reasons("f_lower", boat)))
  chk("decoupled: no stan_data means no claim either way",
      all(is.na(bss_decoupled_reasons(pars, NULL))))
})

# ---------- 13. shared turnover (improvement 2.1, 2026-08-27) ----------
# L was D INDEPENDENT per-day draws with nothing pooling across days, so an observation
# stream covering a SUBSET of days could not move the season level: 148 OSP days left the
# boat median at 1.201 against a prior centre of 1.200, while the OSP/trailer overlap puts
# the real turnover at 2.0-3.0. shared_tau = 1 replaces the D anchors with one estimated
# tau_bar. The guard matters as much as the feature: a shared level is only meaningful when
# L_data is a CONSTANT turnover, never under a time-denominated shore unit where L_data is
# the per-day L_effective regression.
local({
  Pdep <- list(shore_effort_unit = "gear-deployments", tau_shore_prior_mu = 1.7,
               tau_shore_prior_sigma = 0.3, tau_boat_prior_mu = 1.2, tau_boat_prior_sigma = 0.3)
  d6   <- data.frame(event_date = as.Date("2024-12-01") + 0:5, day_type = "weekday")
  sp_dep <- bss_effort_spec(TRUE, d6, Pdep)
  sp_hrs <- bss_effort_spec(TRUE, d6, modifyList(Pdep, list(shore_effort_unit = "crabber-hours")))

  Lc <- rep(1.7, 6); Sc <- rep(0.3, 6)                    # constant turnover
  Lv <- seq(4.5, 6.5, length.out = 6)                     # per-day day-length regression

  off <- bss_shared_tau_data(sp_dep, Lc, Sc, list(), "shore", quiet = TRUE)
  chk("shared tau: OFF by default", identical(off$shared_tau, 0L))
  chk("shared tau: the OFF path still supplies every Stan field",
      all(c("shared_tau","shared_tau_prior_mu","shared_tau_prior_sigma","shared_tau_sigma")
          %in% names(off)))
  chk("shared tau: prior centre and SD come from L_data / L_prior_sigma, not a new constant",
      isTRUE(all.equal(off$shared_tau_prior_mu, 1.7)) &&
      isTRUE(all.equal(off$shared_tau_prior_sigma, 0.3)))
  chk("shared tau: default day-to-day spread is half the prior SD",
      isTRUE(all.equal(off$shared_tau_sigma, 0.15)))

  on <- bss_shared_tau_data(sp_dep, Lc, Sc, list(shared_tau = TRUE), "shore", quiet = TRUE)
  chk("shared tau: ON under a constant turnover", identical(on$shared_tau, 1L))
  chk("shared tau: turning it on changes NOTHING the model is told a priori",
      isTRUE(all.equal(on[setdiff(names(on), "shared_tau")],
                       off[setdiff(names(off), "shared_tau")])))

  # The guard. Both refusals must warn and degrade to OFF, never error: a batch that dies
  # on a misconfigured component is worse than one that runs it the historical way.
  w1 <- NULL
  r1 <- withCallingHandlers(
    bss_shared_tau_data(sp_hrs, Lv, Sc, list(shared_tau = TRUE), "shore", quiet = TRUE),
    warning = function(w) { w1 <<- conditionMessage(w); invokeRestart("muffleWarning") })
  chk("shared tau: REFUSED under a time-denominated unit (L is a day length)",
      identical(r1$shared_tau, 0L) && !is.null(w1) && grepl("turnover", w1))

  w2 <- NULL
  r2 <- withCallingHandlers(
    bss_shared_tau_data(sp_dep, Lv, Sc, list(shared_tau = TRUE), "shore", quiet = TRUE),
    warning = function(w) { w2 <<- conditionMessage(w); invokeRestart("muffleWarning") })
  chk("shared tau: REFUSED when L_data varies across days",
      identical(r2$shared_tau, 0L) && !is.null(w2) && grepl("varies", w2))

  chk("shared tau: shared_tau_sigma is configurable",
      isTRUE(all.equal(bss_shared_tau_data(sp_dep, Lc, Sc,
        list(shared_tau = TRUE, shared_tau_sigma = 0.05), "shore", quiet = TRUE)$shared_tau_sigma,
        0.05)))

  # Stan side: the parameter must be zero-size when off, which is what makes an OFF run
  # bit-identical to the pre-change model (verified against a real fit on 2026-08-27).
  for (f in c("crab_bss_pooled.stan", "crab_bss_gear_resolved.stan")) {
    src <- paste(readLines(file.path("02_stan_models", f), warn = FALSE), collapse = "\n")
    chk(sprintf("%s: tau_bar is zero-size when shared_tau = 0", f),
        grepl("vector<lower=0>[shared_tau] tau_bar;", src, fixed = TRUE))
    chk(sprintf("%s: the shared level is only sampled when it exists", f),
        grepl("if (shared_tau == 1)", src, fixed = TRUE))
    chk(sprintf("%s: tau_bar is reported", f), grepl("tau_bar_out", src, fixed = TRUE))
  }
  chk("decoupled: tau_bar flagged when shared_tau = 0",
      !is.na(bss_decoupled_reasons("tau_bar[1]", list(shared_tau = 0L))))
  chk("decoupled: tau_bar NOT flagged when shared_tau = 1",
      is.na(bss_decoupled_reasons("tau_bar[1]", list(shared_tau = 1L))))
})

# ---------- 14. ar_force per sub-season (2026-08-27b) ----------
# ar_max_resolution had gained the nested per-sub-season form; ar_force had not, and
# as.character() on a one-element list returns its element, so a nested ar_force silently
# forced EVERY sub-season of that population. Stage C of the 2026-08-27 batch was run that
# way: it forced the boat pot closure to biweekly as intended AND the boat all-gear with it,
# moving that component 25,883 -> 28,902 and making the run's port total uninterpretable.
local({
  nested <- list(ar_force = list(private_boat = list(pot_closure = "biweekly")))
  chk("ar_force: nested list forces the NAMED sub-season",
      identical(.bss_resolve_ar_force(nested, "private_boat", "pot_closure"), "biweekly"))
  chk("ar_force: nested list leaves the UNNAMED sub-season alone (the stage C bug)",
      is.null(.bss_resolve_ar_force(nested, "private_boat", "all_gear")))
  chk("ar_force: another population is untouched",
      is.null(.bss_resolve_ar_force(nested, "shore", "all_gear")))
  scal <- list(ar_force = list(private_boat = "biweekly"))
  chk("ar_force: a scalar still forces every sub-season (unchanged behaviour)",
      identical(.bss_resolve_ar_force(scal, "private_boat", "all_gear"), "biweekly") &&
      identical(.bss_resolve_ar_force(scal, "private_boat", "pot_closure"), "biweekly"))
  chk("ar_force: absent means no force", is.null(.bss_resolve_ar_force(list(), "shore", "all_gear")))
  chk("ar_force: an UNNAMED list errors rather than picking one silently",
      inherits(try(.bss_resolve_ar_force(list(ar_force = list(private_boat = list("a","b"))),
                                         "private_boat", "all_gear"), silent = TRUE), "try-error"))
  chk("ar_force: a nested force with no gear_regime supplied does not fire",
      is.null(.bss_resolve_ar_force(nested, "private_boat", NULL)))
})

# ---------- 15. the shared-turnover informed-day floor (5.2, 2026-08-30) ----------
# In the 2026-08-29 batch the toggle was global, and the SHORE all-gear fit has 4 in-window
# I/E days out of 289. Turning shared_tau on moved that component +17.9% on those four
# observations, with an interval that still contained the prior centre, no measurable
# improvement, and no replication in the gear track. The boat, with 130 OSP days, is the case
# the feature exists for. The floor separates them.
local({
  P <- list(shore_effort_unit = "gear-deployments", tau_shore_prior_mu = 1.7,
            tau_shore_prior_sigma = 0.3, tau_boat_prior_mu = 1.2, tau_boat_prior_sigma = 0.3)
  d6 <- data.frame(event_date = as.Date("2024-12-01") + 0:5, day_type = "weekday")
  sp <- bss_effort_spec(TRUE, d6, P)
  Lc <- rep(1.7, 6); Sc <- rep(0.3, 6)
  on <- function(n, extra = list())
    bss_shared_tau_data(sp, Lc, Sc, modifyList(list(shared_tau = TRUE), extra),
                        "x", n_informed = n, quiet = TRUE)$shared_tau

  chk("shared tau floor: 4 informed days (the shore case) is REFUSED", identical(on(4L), 0L))
  chk("shared tau floor: 130 informed days (the boat case) is allowed", identical(on(130L), 1L))
  chk("shared tau floor: 18 (boat pot closure) is allowed at the default of 15",
      identical(on(18L), 1L))
  chk("shared tau floor: 0 informed days is refused", identical(on(0L), 0L))
  chk("shared tau floor: the threshold is configurable",
      identical(on(18L, list(shared_tau_min_obs = 20L)), 0L))
  chk("shared tau floor: an UNKNOWN count does not block the feature",
      identical(on(NA_integer_), 1L))
  chk("shared tau floor: it cannot switch the feature ON when the unit guard says no",
      identical(bss_shared_tau_data(
        bss_effort_spec(TRUE, d6, modifyList(P, list(shore_effort_unit = "crabber-hours"))),
        seq(4.5, 6.5, length.out = 6), Sc, list(shared_tau = TRUE), "x",
        n_informed = 999L, quiet = TRUE)$shared_tau, 0L))
  chk("shared tau floor: off by default regardless of the count",
      identical(bss_shared_tau_data(sp, Lc, Sc, list(), "x", n_informed = 999L,
                                    quiet = TRUE)$shared_tau, 0L))
})

# ---------- 16. model adequacy, reported beside the gate (5.5, 2026-08-30) ----------
# Six configurations passed every gate criterion while spanning 44% on the boat. The gate
# asks whether the sampler worked; these ask whether the model is carrying the data. The
# thresholds are exercised against the real numbers the 2026-08-29 batch produced.
local({
  chk("adequacy: module exposes both entry points",
      exists("bss_model_adequacy", mode = "function") &&
      exists("write_model_adequacy", mode = "function"))

  # The stage F signature, from that run's own files: p_loo/n doubles, a second bad k
  # appears, and the PIT bias widens. Read straight off disk so the test is about the
  # quantities, not about a mock.
  d_off <- "05_output/20260828/pooled-CPUE-IP-D-tau-off"
  d_esc <- "05_output/20260829/pooled-CPUE-IP-F-escalate"
  f <- "loo_summary_private_boat_all_gear_Dungeness_Kept.csv"
  if (file.exists(file.path(d_off, f)) && file.exists(file.path(d_esc, f))) {
    a <- read.csv(file.path(d_off, f)); b <- read.csv(file.path(d_esc, f))
    fa <- max(a$p_loo / a$n_obs); fb <- max(b$p_loo / b$n_obs)
    chk("adequacy: p_loo fraction separates the monthly and daily boat fits",
        fb > 1.9 * fa, sprintf("%.3f -> %.3f", fa, fb))
    chk("adequacy: the daily fit carries more unreliable LOO points",
        sum(b[["n_pareto_k_gt_0.7"]]) > sum(a[["n_pareto_k_gt_0.7"]]))
  } else {
    skp("adequacy: the threshold check's reference runs", "runs absent")
  }

  # A decoupled dispersion parameter must NOT drag the minimum down: an unused r_OSP in a
  # shore fit would otherwise fire the flag on every run.
  shore <- list(IE_n = 0L, OSP_n = 0L, T_n = 0L, osp_scale_is_tau = 1L, K_open = 0L,
                apply_crab_fraction = 0L, crab_fraction_estimate = 0L, osp_crab_lower = 0L,
                estimate_cpue_density = 0L, w = c(1,0,1), holiday = c(1,0,0), shared_tau = 0L)
  r <- bss_decoupled_reasons(c("r_E","r_C","r_OSP","sigma_r_OSP"), shore)
  chk("adequacy: shore r_OSP is excluded from the dispersion floor as decoupled",
      is.na(r[1]) && is.na(r[2]) && !is.na(r[3]) && !is.na(r[4]))
})

# ---------- 17. reporting integrity (2026-08-30) ----------
# One principle in four places: a number that is not an estimate must not sit in a column
# that reads like one. Every case below was found in a real run output.
local({
  # 17a. structural_params gains `estimate`, NA when decoupled, and carries the gate verdict.
  src <- paste(readLines("03_R_functions/model_diagnostics.R", warn = FALSE), collapse = "\n")
  chk("reporting: structural_params has an `estimate` column that is NA when decoupled",
      grepl("out$estimate <- ifelse(out$decoupled, NA_real_, out$median)", src, fixed = TRUE))
  chk("reporting: structural_params carries the fit's gate verdict",
      grepl("out$fit_method <- fit_method", src, fixed = TRUE))

  # 17b. pe_vs_bss_comparison says whether the BSS column was USED. The 2026-08-29 gear run
  # left BSS_catch = 32,689 beside method_selected = "PE" with nothing marking it unused.
  for (rmd in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd",
                "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")) {
    t <- paste(readLines(rmd, warn = FALSE), collapse = "\n")
    chk(sprintf("reporting: %s marks whether the BSS column was reported", basename(rmd)),
        grepl('comparison_df$bss_reported <- comparison_df$method_selected == "BSS"',
              t, fixed = TRUE))
  }

  # 17c. the pooled expansion table flags a decoupled R_G_boat (the gear track already did).
  t <- paste(readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE), collapse = "\n")
  chk("reporting: pooled expansion_ratios flags a prior-only R_G_boat",
      grepl("R_G_boat_decoupled", t, fixed = TRUE) && grepl("PRIOR ONLY", t, fixed = TRUE))

  # 17d. historical folders can be annotated after the fact. The run outputs are committed on
  # purpose and must not be rewritten, so the audit is a NEW file beside them.
  chk("reporting: the retro-annotator is available",
      exists("annotate_decoupled_run", mode = "function"))
  ref <- "05_output/20260826/pooled-CPUE-PV4-minint/decoupled_audit.csv"
  if (file.exists(ref)) {
    a <- read.csv(ref, stringsAsFactors = FALSE)
    gk <- function(fit, par, col) a[[col]][a$fit == fit & a$parameter == par][1]
    chk("retro-audit: shore r_OSP (posterior mean in the hundreds of thousands) reads prior only",
        identical(gk("shore_all_gear_Dungeness_Kept", "r_OSP", "status"), "prior only"))
    chk("retro-audit: shore R_G_boat reads prior only",
        identical(gk("shore_all_gear_Dungeness_Kept", "R_G_boat", "status"), "prior only"))
    chk("retro-audit: boat kappa_OSP reads prior only under osp_scale_is_tau",
        identical(gk("private_boat_all_gear_Dungeness_Kept", "kappa_OSP", "status"), "prior only"))
    chk("retro-audit: the two GENUINE estimates are not mislabelled",
        identical(gk("private_boat_all_gear_Dungeness_Kept", "R_G_boat", "status"), "estimate") &&
        identical(gk("shore_all_gear_Dungeness_Kept", "sigma_IE", "status"), "estimate"))
    chk("retro-audit: every row is marked reconstructed", all(a$source == "reconstructed"))
    chk("retro-audit: window-dependent rules are reported unknown, not guessed",
        any(grepl("^unknown", a$status)))
  } else {
    skp("retro-audit: the annotated baseline", "not yet annotated")
  }
})

# ---------------------------------------------------------------------------
# 18. The sampler-override escape hatch (2026-08-30). params_model wins the driver merge,
#     so a run_config delta cannot reach a sampler setting; without this hatch, Stage 5's S3
#     would run the gear track at its 1,000-draw default and report that more draws did not
#     help. The failure mode being guarded against is SILENCE, so the whitelist must ERROR
#     on a non-sampler key rather than dropping it.
# ---------------------------------------------------------------------------
local({
  chk("sampler override: absent / NULL is a no-op",
      identical(bss_apply_sampler_override(list(a = 1), NULL, quiet = TRUE), list(a = 1)))
  p <- bss_apply_sampler_override(list(bss_iter_default = 2000, bss_warmup_default = 1000),
                                  list(bss_iter_default = 5000, bss_warmup_default = 2500,
                                       bss_adapt_delta_default = 0.99,
                                       bss_max_treedepth_default = 13), quiet = TRUE)
  chk("sampler override: raises iterations, warmup, adapt_delta and treedepth together",
      identical(p$bss_iter_default, 5000) && identical(p$bss_warmup_default, 2500) &&
      identical(p$bss_adapt_delta_default, 0.99) && identical(p$bss_max_treedepth_default, 13))
  chk("sampler override: a STRUCTURAL key is refused, not silently dropped",
      inherits(try(bss_apply_sampler_override(list(), list(shared_tau = TRUE), quiet = TRUE),
                   silent = TRUE), "try-error"))
  chk("sampler override: an unnamed list is refused",
      inherits(try(bss_apply_sampler_override(list(), list(5000), quiet = TRUE),
                   silent = TRUE), "try-error"))
  chk("sampler override: it also reads params$bss_sampler_override, not just the argument",
      identical(bss_apply_sampler_override(
        list(bss_chains = 4, bss_sampler_override = list(bss_chains = 2)),
        quiet = TRUE)$bss_chains, 2))
  for (rmd in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd",
                "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")) {
    t <- readLines(rmd, warn = FALSE)
    # fixed = TRUE: the parentheses in the merge line are regex metacharacters.
    i_merge <- grep("params <- modifyList(run_config, params)", t, fixed = TRUE)
    i_hook  <- grep("bss_apply_sampler_override(params", t, fixed = TRUE)
    chk(sprintf("sampler override: %s applies it AFTER the merge", basename(rmd)),
        length(i_merge) == 1 && length(i_hook) == 1 && i_hook > i_merge)
  }
  e <- new.env(); sys.source("run_config.R", envir = e)
  chk("sampler override: registered in run_config and NULL in production",
      "bss_sampler_override" %in% names(e$run_config) && is.null(e$run_config$bss_sampler_override))

  # 18b. One event, one wording. The gate is the single authority on method selection; a
  # second writer inventing its own string means the same fit reads differently in two
  # files in the same folder.
  for (rmd in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd",
                "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")) {
    t <- paste(readLines(rmd, warn = FALSE), collapse = "\n")
    chk(sprintf("method string: %s quotes the gate's verdict, not its own", basename(rmd)),
        grepl("fit_method = b$gate_info$method_selected", t, fixed = TRUE))
  }
  # Comment lines are excluded: a comment explaining why the string is wrong is not a
  # writer emitting it.
  emits <- function(p) {
    l <- readLines(p, warn = FALSE)
    l <- l[!grepl("^\\s*#", l)]
    any(grepl('"PE (gate fail)"', l, fixed = TRUE))
  }
  chk("method string: no writer invents 'PE (gate fail)'",
      !any(vapply(c(list.files("03_R_functions", full.names = TRUE),
                    list.files("01_BSS_models", pattern = "[.]Rmd$", full.names = TRUE)),
                  emits, logical(1))))
})

# ---------------------------------------------------------------------------
# 19. Retro model adequacy (plan 5.5). Plan 5.3 asks for p_loo, r_OSP and the two PIT means
#     "for each cell" of a 2x2 whose archived cells predate the diagnostic. A table with
#     adequacy for the new cells and blanks for the old ones is the shape of an argument
#     that quietly favours whichever cells are new. The reconstruction must therefore agree
#     with the live path exactly, which is why both go through .bma_core.
# ---------------------------------------------------------------------------
local({
  chk("adequacy: the retro path is available",
      exists("annotate_model_adequacy_run", mode = "function") &&
      exists(".bma_core", mode = "function"))
  d <- "05_output/20260829/pooled-CPUE-IP-F-escalate"
  if (!dir.exists(d)) skp("adequacy retro: stage F annotation", "folder absent")
  if (dir.exists(d)) {
    a <- annotate_model_adequacy_run(d, overwrite = TRUE, quiet = TRUE)
    r <- a[a$fit == "private_boat_all_gear_Dungeness_Kept", ]
    # The values the live path produced on this folder when it was validated (2026-08-30).
    chk("adequacy retro: stage F boat p_loo fraction reproduces the live value",
        isTRUE(abs(r$p_loo_frac - 0.163) < 0.002), sprintf("(%.4f, %s)", r$p_loo_frac, r$p_loo_worst_stream))
    chk("adequacy retro: stage F Pareto k > 0.7 count reproduces", isTRUE(r$n_pareto_bad == 2))
    chk("adequacy retro: stage F dispersion n_eff finds sigma_r_OSP below the gate's floor",
        isTRUE(r$disp_neff_min < 400) && identical(r$disp_neff_min_par, "sigma_r_OSP"),
        sprintf("(%s at n_eff %s)", r$disp_neff_min_par, r$disp_neff_min))
    chk("adequacy retro: every row is marked as reconstructed", all(grepl("^reconstructed", a$source)))
    chk("adequacy retro: method_selected is recovered from the convergence report",
        all(!is.na(a$method_selected)))
    s <- a[a$fit == "shore_all_gear_Dungeness_Kept", ]
    chk("adequacy retro: a decoupled r_OSP in a shore fit is excluded from the floor",
        !identical(s$disp_neff_min_par, "r_OSP"), sprintf("(shore minimum is %s)", s$disp_neff_min_par))
    # 19b. SHARPNESS (2026-08-30). The Stage 5 2x2 proved pit_mean alone is not merely
    #      incomplete but actively misleading: the (tau ON, daily) cell has the BEST pit_mean
    #      of the four and the worst calibration in the batch, because a latent process with
    #      one state per observation piles every PIT value at 0.5. coverage_50 and pit_sd
    #      were already being written and were not reaching the adequacy table.
    s4b <- annotate_model_adequacy_run("05_output/20260830/pooled-CPUE-S5-4b-daily-tauon",
                                       overwrite = TRUE, quiet = TRUE)
    s2  <- annotate_model_adequacy_run("05_output/20260829/pooled-CPUE-S5-2-tau-pooled",
                                       overwrite = TRUE, quiet = TRUE)
    if (!is.null(s4b) && !is.null(s2)) {
      b <- s4b[grepl("private_boat_all_gear", s4b$fit), ]
      g <- s2 [grepl("private_boat_all_gear", s2$fit),  ]
      chk("adequacy: pit_mean alone RANKS THE CELLS WRONG (the reason for the new columns)",
          isTRUE(b$pit_worst_bias < g$pit_worst_bias),
          sprintf("(daily %.3f 'better' than monthly %.3f)", b$pit_worst_bias, g$pit_worst_bias))
      chk("adequacy: coverage_50 deviation puts them back in the right order",
          isTRUE(b$cov50_worst_dev > g$cov50_worst_dev),
          sprintf("(daily %.3f vs monthly %.3f; nominal deviation 0)", b$cov50_worst_dev, g$cov50_worst_dev))
      chk("adequacy: PIT sd deviation agrees with coverage",
          isTRUE(b$pit_sd_worst_dev > g$pit_sd_worst_dev),
          sprintf("(daily %.3f vs monthly %.3f)", b$pit_sd_worst_dev, g$pit_sd_worst_dev))
      chk("adequacy: the daily cell is flagged miscalibrated", isTRUE(b$flag_miscalibrated))
      chk("adequacy: dispersion COLLAPSE is reported where n_eff cannot see it",
          isTRUE(b$disp_scale_min < 0.2) && isTRUE(g$disp_scale_min > 0.5) &&
          isTRUE(b$disp_neff_min > 400),
          sprintf("(daily sigma %.3f at n_eff %.0f, above the floor; monthly %.3f)",
                  b$disp_scale_min, b$disp_neff_min, g$disp_scale_min))
    }
    e <- "05_output/20260829/pooled-CPUE-IP-E-tau-on-pooled"
    if (!dir.exists(e)) skp("adequacy retro: stage E", "folder absent")
    if (dir.exists(e)) {
      b <- annotate_model_adequacy_run(e, overwrite = TRUE, quiet = TRUE)
      rb <- b[b$fit == "private_boat_all_gear_Dungeness_Kept", ]
      chk("adequacy: the shared turnover improved PIT bias at UNCHANGED complexity",
          isTRUE(abs(rb$p_loo_frac - r$p_loo_frac) > 0.05) && isTRUE(rb$pit_worst_bias < 0.05),
          sprintf("(E p_loo %.3f bias %.3f vs F p_loo %.3f bias %.3f)",
                  rb$p_loo_frac, rb$pit_worst_bias, r$p_loo_frac, r$pit_worst_bias))
    }
  } else skp("adequacy retro: stage F", "folder absent")
})

# ---------------------------------------------------------------------------
# 21. prior_vs_posterior row resolution (regression, 2026-08-30 Stage 5 batch).
#     tau_bar is declared vector<lower=0>[shared_tau], so rstan names its summary row
#     "tau_bar[1]". The prior table entry was called "tau_bar"; has_par() strips the index
#     and selected it, then post[pn, ] threw, and the enclosing tryCatch swallowed the error
#     -- silently dropping the ENTIRE boat prior_vs_posterior file in S2, S3 and S4b. The
#     lesson generalises: a tryCatch around a whole writer converts a one-row bug into a
#     missing file, so the row lookup is now tolerant in both directions and a parameter
#     that cannot be resolved is skipped rather than taking the file with it.
# ---------------------------------------------------------------------------
local({
  chk("pvp: the row resolver exists", exists(".pvp_row_key", mode = "function"))
  rn <- c("R_G", "mu_mu_E[1]", "tau_bar[1]")
  chk("pvp: an indexed name matches its indexed row",
      identical(.pvp_row_key("tau_bar[1]", rn), "tau_bar[1]"))
  chk("pvp: a BARE name still finds the indexed row (the bug)",
      identical(.pvp_row_key("tau_bar", rn), "tau_bar[1]"))
  chk("pvp: an indexed name still finds a scalar row",
      identical(.pvp_row_key("mu_mu_E[1]", c("mu_mu_E")), "mu_mu_E"))
  chk("pvp: an absent parameter resolves to NA, so the caller can skip it",
      is.na(.pvp_row_key("not_a_par", rn)))
  t <- paste(readLines("03_R_functions/save_run_diagnostics.R", warn = FALSE), collapse = "\n")
  chk("pvp: the shared-turnover entry is registered under its indexed name",
      grepl("prior_tbl$`tau_bar[1]`", t, fixed = TRUE))
  chk("pvp: an unresolvable parameter is dropped, not fatal to the file",
      grepl("rows <- Filter(Negate(is.null), rows)", t, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 20. The Stage 5 batch runner. Its REF block hard-codes baselines so a criterion cannot be
#     edited after the answer is known; this checks the hard-coded values still match the
#     folders they name, which is the one way that arrangement can rot silently.
# ---------------------------------------------------------------------------
local({
  f <- "06_diagnostics/run_stage5_2026-08-30.R"
  if (!file.exists(f)) { chk("stage5 runner: present", FALSE); return(invisible(NULL)) }
  chk("stage5 runner: parses", !inherits(try(parse(f), silent = TRUE), "try-error"))
  t <- paste(readLines(f, warn = FALSE), collapse = "\n")
  chk("stage5 runner: ships with DRY_RUN TRUE", grepl("^DRY_RUN <- TRUE", t) ||
      grepl("\nDRY_RUN <- TRUE", t))
  chk("stage5 runner: the 2x2 forces ONE sub-season, per-sub-season",
      grepl('AR_BOAT_AG_DAILY <- list(private_boat = list(all_gear    = "daily"))', t, fixed = TRUE))
  chk("stage5 runner: S4a exists, so the (off, daily) cell is matched rather than stage F's",
      grepl('S4a = stage("S4a"', t, fixed = TRUE) && grepl('shared_tau = FALSE, ar_force', t, fixed = TRUE))
  # REF drift. Read the committed folders the runner names and confirm the numbers agree.
  ref <- list(
    list("05_output/20260828/pooled-CPUE-IP-D-tau-off", "private_boat \\(All gear\\)", 25868),
    list("05_output/20260829/pooled-CPUE-IP-E-tau-on-pooled", "private_boat \\(All gear\\)", 31008),
    list("05_output/20260829/pooled-CPUE-IP-F-escalate", "private_boat \\(All gear\\)", 37359),
    # The 2026-08-30 batch's own cells. The review in stage5-batch-review-2026-08-31.md
    # quotes these four numbers and the interaction computed from them; if a folder is ever
    # re-run in place, the review's arithmetic silently stops matching its sources.
    list("05_output/20260829/pooled-CPUE-S5-2-tau-pooled", "private_boat \\(All gear\\)", 31008),
    list("05_output/20260829/pooled-CPUE-S5-4a-daily-tauoff", "private_boat \\(All gear\\)", 37359),
    list("05_output/20260830/pooled-CPUE-S5-4b-daily-tauon", "private_boat \\(All gear\\)", 42344),
    list("05_output/20260830/pooled-CPUE-S5-5-boatpc-ar", "private_boat \\(All gear\\)", 25868),
    list("05_output/20260829/gear-type-CPUE-model-S5-3-tau-gear", "private_boat \\(All gear\\)", 30760))
  for (r in ref) {
    p <- file.path(r[[1]], "pe_vs_bss_comparison.csv")
    if (!file.exists(p)) { skp(paste("stage5 REF:", basename(r[[1]])), "folder absent"); next }
    d <- read.csv(p, stringsAsFactors = FALSE)
    v <- d$BSS_catch[grepl(r[[2]], d$component)][1]
    chk(sprintf("stage5 REF: %s boat all-gear is still %d", basename(r[[1]]), r[[3]]),
        isTRUE(round(v) == r[[3]]), sprintf("(read %s)", round(v)))
  }
})

# ---------------------------------------------------------------------------
# 22. The Stage 5 conclusion, guarded arithmetically (2026-08-31). The recommendation to
#     adopt the shared turnover and NOT the daily AR rests on an interaction near zero and
#     on the two levers ranking oppositely on calibration. Both are recomputed here from the
#     committed outputs, so the conclusion cannot drift away from its evidence unnoticed.
# ---------------------------------------------------------------------------
local({
  bag <- function(d) {
    p <- file.path(d, "pe_vs_bss_comparison.csv")
    if (!file.exists(p)) return(NA_real_)
    x <- read.csv(p, stringsAsFactors = FALSE)
    round(x$BSS_catch[grepl("private_boat \\(All gear\\)", x$component)][1])
  }
  off_m <- bag("05_output/20260828/pooled-CPUE-IP-D-tau-off")
  on_m  <- bag("05_output/20260829/pooled-CPUE-S5-2-tau-pooled")
  off_d <- bag("05_output/20260829/pooled-CPUE-S5-4a-daily-tauoff")
  on_d  <- bag("05_output/20260830/pooled-CPUE-S5-4b-daily-tauon")
  if (all(is.finite(c(off_m, on_m, off_d, on_d)))) {
    inter <- (on_d - off_d) - (on_m - off_m)
    chk("stage5: the tau x AR interaction is still near zero (the ADDITIVE finding)",
        abs(inter) < 0.15 * (on_m - off_m),
        sprintf("(interaction %d against main effects +%d and +%d)", inter, on_m - off_m, off_d - off_m))
  }
  # Read the RANDOMIZED coverage from ppc_byobs, not ppc_calibration: runs committed before
  # 2026-09-01 carry a non-randomized coverage_50 in the aggregate file that over-covers
  # small counts by construction (up to 0.154 on the trailer stream).
  cov50 <- function(d, stream = "trailer") {
    p <- file.path(d, "ppc_byobs_private_boat_all_gear_Dungeness_Kept.csv")
    if (!file.exists(p)) return(NA_real_)
    x <- read.csv(p, stringsAsFactors = FALSE)
    mean(as.logical(x$in_50[x$data_type == stream]), na.rm = TRUE)
  }
  cm <- cov50("05_output/20260829/pooled-CPUE-S5-2-tau-pooled")
  cd <- cov50("05_output/20260830/pooled-CPUE-S5-4b-daily-tauon")
  chk("stage5: the daily cell is still the miscalibrated one (nominal coverage_50 = 0.50)",
      isTRUE(abs(cd - 0.5) > 3 * abs(cm - 0.5)),
      sprintf("(monthly %.3f vs daily %.3f)", cm, cd))
  # tau_bar agrees across two independently parameterized tracks.
  tb <- function(d) {
    p <- file.path(d, "structural_params_private_boat_all_gear_Dungeness_Kept.csv")
    if (!file.exists(p)) return(NA_real_)
    x <- read.csv(p, stringsAsFactors = FALSE); x$median[x$parameter == "tau_bar[1]"][1]
  }
  tp <- tb("05_output/20260829/pooled-CPUE-S5-2-tau-pooled")
  tg <- tb("05_output/20260829/gear-type-CPUE-model-S5-3-tau-gear")
  chk("stage5: the two tracks still agree on tau_bar to better than 1%",
      isTRUE(abs(tg - tp) / tp < 0.01), sprintf("(pooled %.4f vs gear %.4f)", tp, tg))
  chk("stage5: tau_bar at monthly still sits inside the external overlap calibration 2.01-3.03",
      isTRUE(tp > 2.01 && tp < 3.03), sprintf("(%.3f)", tp))
})

# ---------------------------------------------------------------------------
# 23. Randomized PPC coverage (2026-09-01). ppc_calibration_*.csv and ppc_byobs_*.csv
#     computed the SAME quantity two ways and disagreed by up to 0.154 on the private-boat
#     trailer stream, because the aggregate file tested the observation against a quantile
#     interval of simulated draws. That over-covers small counts by construction and is what
#     produced the phantom "trailer over-coverage" open item in the 2026-08-31 review.
# ---------------------------------------------------------------------------
local({
  t <- paste(readLines("03_R_functions/model_diagnostics.R", warn = FALSE), collapse = "\n")
  chk("ppc: aggregate coverage is read off the randomized PIT",
      grepl("cov50[i] <- pit[i] >= 0.25", t, fixed = TRUE) &&
      grepl("cov95[i] <- pit[i] >= 0.025", t, fixed = TRUE))
  chk("ppc: the quantile-interval coverage is gone",
      !grepl("cov50[i] <- y[i] >= qq[2]", t, fixed = TRUE))
  t2 <- paste(readLines("03_R_functions/save_run_diagnostics.R", warn = FALSE), collapse = "\n")
  chk("ppc: the OSP stream now gets per-observation rows",
      grepl('parts$osp <- cbind(data_type = "osp"', t2, fixed = TRUE))
  # 2026-09-04: the arithmetic moved into 03_R_functions/zinb_ppc.R so the NB2 and the
  # mixture forms cannot drift between this file and model_diagnostics.R. Assert the
  # PROPERTY (p_zero is written, over all days) rather than the literal expression.
  chk("ppc: p_zero is written, so the zero bin can be scored against ALL days",
      grepl("p_zero[i] <- bss_zi_p_zero(mu, sz, th)", t2, fixed = TRUE))
  chk("adequacy: coverage prefers the randomized source and records which it used",
      exists(".bma_core", mode = "function") &&
      "byobs" %in% names(formals(.bma_core)))
  # The correction must not rescue the daily-AR cell: if it did, the Stage 5 recommendation
  # was resting on the artefact.
  cv <- function(d, stream = "trailer") {
    p <- file.path(d, "ppc_byobs_private_boat_all_gear_Dungeness_Kept.csv")
    if (!file.exists(p)) return(NA_real_)
    x <- read.csv(p, stringsAsFactors = FALSE)
    mean(as.logical(x$in_50[x$data_type == stream]), na.rm = TRUE)
  }
  s2 <- cv("05_output/20260829/pooled-CPUE-S5-2-tau-pooled")
  s4 <- cv("05_output/20260830/pooled-CPUE-S5-4b-daily-tauon")
  sdc <- sqrt(0.25 / 195)
  chk("ppc: production is CALIBRATED on the corrected statistic (the retracted open item)",
      isTRUE(abs(s2 - 0.5) / sdc < 2), sprintf("(%.3f, %.1f sampling SDs)", s2, abs(s2 - 0.5) / sdc))
  chk("ppc: the daily-AR cell is still broken on the corrected statistic",
      isTRUE(abs(s4 - 0.5) / sdc > 3), sprintf("(%.3f, %.1f sampling SDs)", s4, abs(s4 - 0.5) / sdc))
  chk("adequacy: flag_miscalibrated now DISCRIMINATES instead of firing on everything",
      { f <- function(d) { a <- annotate_model_adequacy_run(d, overwrite = TRUE, quiet = TRUE)
          a$flag_miscalibrated[grepl("private_boat_all_gear", a$fit)][1] }
        isFALSE(f("05_output/20260829/pooled-CPUE-S5-2-tau-pooled")) &&
        isTRUE(f("05_output/20260830/pooled-CPUE-S5-4b-daily-tauon")) })
})

# ---------------------------------------------------------------------------
# 24. The 2026-09-01 adoption, and the validation batch that tests it.
# ---------------------------------------------------------------------------
local({
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("adoption: run_config ships shared_tau = TRUE", isTRUE(rc$shared_tau))
  chk("adoption: the floor is STATED in run_config, not left to the helper default",
      identical(rc$shared_tau_min_obs, 15))
  chk("adoption: the boat-only claim still rests on the floor, not a population switch",
      { s <- paste(readLines("03_R_functions/bss_effort_spec.R", warn = FALSE), collapse = "\n")
        grepl("shared_tau_min_obs", s, fixed = TRUE) && grepl("n_informed", s, fixed = TRUE) })
  f <- "06_diagnostics/run_validation_2026-09-01.R"
  if (file.exists(f)) {
    chk("validation runner: parses", !inherits(try(parse(f), silent = TRUE), "try-error"))
    t <- paste(readLines(f, warn = FALSE), collapse = "\n")
    chk("validation runner: ships with DRY_RUN TRUE", grepl("\nDRY_RUN <- TRUE", t))
    chk("validation runner: V1 carries NO delta, so it tests production as shipped",
        grepl('V1 = stage("V1", "VAL-1-adopted", "pooled", list(),', t, fixed = TRUE))
    chk("validation runner: it refuses to run V1 if the adoption is not in run_config",
        grepl("does not ship shared_tau = TRUE", t, fixed = TRUE))
  } else chk("validation runner: present", FALSE)
})

# ---------------------------------------------------------------------------
# 25. The 2026-09-01 validation batch: the gear-track production fix, and the guard against
#     the reference-differs-in-two-ways mistake that produced two spurious FAILs.
# ---------------------------------------------------------------------------
local({
  rmd <- "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"
  t <- readLines(rmd, warn = FALSE); tt <- paste(t, collapse = "\n")
  chk("gear driver: boat all-gear now has per-fit sampler settings (it had NONE)",
      all(vapply(c("bss_iter_boat_allgear", "bss_warmup_boat_allgear",
                   "bss_treedepth_boat_allgear", "bss_delta_boat_allgear"),
                 function(k) grepl(k, tt, fixed = TRUE), logical(1))))
  chk("gear driver: they match the pooled track's settings for the same component",
      grepl("bss_iter_boat_allgear       = 5000", tt, fixed = TRUE) &&
      grepl("bss_delta_boat_allgear      = 0.99", tt, fixed = TRUE) &&
      grepl("bss_treedepth_boat_allgear  = 13", tt, fixed = TRUE))
  chk("gear driver: a dedicated branch selects them, so a future fit does not inherit them",
      grepl("} else if (!is_shore && is_allgear) {", tt, fixed = TRUE))
  # The branch must come BEFORE the catch-all else, or it is unreachable.
  i_boat <- grep("} else if (!is_shore && is_allgear) {", t, fixed = TRUE)
  i_dflt <- grep("fit_treedep <- params$bss_max_treedepth_default", t, fixed = TRUE)
  chk("gear driver: the boat branch precedes the catch-all default branch",
      length(i_boat) == 1 && length(i_dflt) >= 1 && i_boat < max(i_dflt))
  e <- new.env(); sys.source("run_config.R", envir = e)
  chk("gear fix is in the DRIVER, not the experiment hatch (production override stays NULL)",
      is.null(e$run_config$bss_sampler_override))

  f <- "06_diagnostics/run_validation_2026-09-01.R"
  if (file.exists(f)) {
    v <- paste(readLines(f, warn = FALSE), collapse = "\n")
    chk("validation runner: config_delta guards against a reference differing in two ways",
        grepl("config_delta <- function(dir_a, dir_b)", v, fixed = TRUE) &&
        grepl("expect_delta", v, fixed = TRUE))
    chk("validation runner: the V2 boat criterion uses the same-sampler reference",
        grepl("REF$Egear$dir), \"private_boat\"", v, fixed = TRUE))
    chk("validation runner: D3 picks up the newest production run for the zero bin",
        grepl("V1 PRODUCTION", v, fixed = TRUE) && grepl("VAL-1-adopted", v, fixed = TRUE))
    chk("validation runner: the zero bin is scored as a z, not a raw fraction",
        grepl("exp_zeros_sd", v, fixed = TRUE))
    chk("validation runner: reset to DRY_RUN TRUE after the batch", grepl("\nDRY_RUN <- TRUE", v))
  } else chk("validation runner present", FALSE)

  t2 <- paste(readLines("03_R_functions/model_diagnostics.R", warn = FALSE), collapse = "\n")
  chk("ppc: the aggregate PIT is the EXACT expectation, so the two files agree exactly",
      grepl("pit[i]   <- bss_zi_pit(y[i], mu_k, sz_k, th_k)", t2, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 26. Results the 2026-09-01 batch established, guarded arithmetically so the write-up
#     cannot drift from its evidence.
# ---------------------------------------------------------------------------
local({
  bag <- function(d, pat = "private_boat \\(All gear\\)", col = "BSS_catch") {
    p <- file.path(d, "pe_vs_bss_comparison.csv")
    if (!file.exists(p)) return(NA_real_)
    x <- read.csv(p, stringsAsFactors = FALSE); round(x[[col]][grepl(pat, x$component)][1])
  }
  v1 <- "05_output/20260831/pooled-CPUE-VAL-1-adopted"
  if (!dir.exists(v1)) skp("V1: the adopted configuration's components", "folder absent")
  if (dir.exists(v1)) {
    chk("V1: the adopted configuration reproduces S2's components",
        isTRUE(bag(v1) == 31008) &&
        isTRUE(bag(v1, "shore \\(All gear\\)") == 20898))
    chk("V1: the boat prior_vs_posterior file exists again (the 2026-08-31 regression)",
        file.exists(file.path(v1, "prior_vs_posterior_private_boat_all_gear_Dungeness_Kept.csv")))
    chk("V1: ppc_byobs now carries p_zero and the OSP stream",
        { x <- read.csv(file.path(v1, "ppc_byobs_private_boat_all_gear_Dungeness_Kept.csv"),
                        stringsAsFactors = FALSE)
          "p_zero" %in% names(x) && "osp" %in% x$data_type })
    chk("V1: production's SHORE all-gear fit is the one flagged, not the boat",
        { a <- read.csv(file.path(v1, "model_adequacy.csv"), stringsAsFactors = FALSE)
          isTRUE(a$flag_miscalibrated[a$fit == "shore_all_gear_Dungeness_Kept"]) &&
          isFALSE(a$flag_miscalibrated[a$fit == "private_boat_all_gear_Dungeness_Kept"]) },
        "(shore daily AR: p_loo 35% of n_obs, 41 Pareto k > 0.7, coverage_50 0.701)")
  }
  # The boat pot-closure 2x2 is now complete and additive.
  cells <- c("05_output/20260828/pooled-CPUE-IP-D-tau-off",
             "05_output/20260831/pooled-CPUE-VAL-1-adopted",
             "05_output/20260830/pooled-CPUE-S5-5-boatpc-ar",
             "05_output/20260901/pooled-CPUE-VAL-4-tau-boatpc-biwk")
  if (all(dir.exists(cells))) {
    v <- vapply(cells, bag, numeric(1), pat = "private_boat \\(Pot closure\\)")
    inter <- (v[4] - v[3]) - (v[2] - v[1])
    chk("the boat pot-closure 2x2 is additive too (interaction near zero)",
        abs(inter) < 0.3 * (v[2] - v[1]),
        sprintf("(%d/%d/%d/%d, interaction %+d)", v[1], v[2], v[3], v[4], inter))
  }
  v5 <- "05_output/20260901/pooled-CPUE-VAL-5-floor20"
  if (!(dir.exists(v5) && dir.exists(v1))) skp("V5: the shared_tau_min_obs threshold", "folder absent")
  if (dir.exists(v5) && dir.exists(v1))
    chk("V5: the shared_tau_min_obs threshold is worth ~169 crab, 0.24% of the port",
        isTRUE(abs((bag(v1, "private_boat \\(Pot closure\\)") -
                    bag(v5, "private_boat \\(Pot closure\\)")) - 169) <= 2),
        sprintf("(floor 15: %d, floor 20: %d)", bag(v1, "private_boat \\(Pot closure\\)"),
                bag(v5, "private_boat \\(Pot closure\\)")))
  # GR-7 Phase 2, sampled for the first time: widen intervals, do not move medians.
  a <- "05_output/20260901/gear-type-CPUE-model-VAL-2-gearG-phase1"
  b <- "05_output/20260901/gear-type-CPUE-model-VAL-3-gearG-dirichlet"
  if (!(dir.exists(a) && dir.exists(b))) skp("VAL-2 / VAL-3 gear intervals", "folder absent")
  if (dir.exists(a) && dir.exists(b)) {
    g <- function(d) { x <- read.csv(file.path(d, "catch_by_gear_type_detail.csv"), stringsAsFactors = FALSE)
      x <- x[x$population == "shore" & x$subseason == "all_gear", ]
      x$w <- x$BSS_hi95 - x$BSS_lo95; x[order(x$gear_type), ] }
    A <- g(a); B <- g(b)
    chk("GR-7 Phase 2: Dirichlet shares widen every per-gear interval",
        all(B$w / A$w > 1.05), sprintf("(ratios %s)", paste(sprintf("%.2f", B$w / A$w), collapse = ", ")))
    chk("GR-7 Phase 2: and move no median by more than 1%",
        max(abs(100 * (B$BSS_median - A$BSS_median) / A$BSS_median)) < 1,
        sprintf("(largest %.1f%%)", max(abs(100 * (B$BSS_median - A$BSS_median) / A$BSS_median))))
  }
})

# ---------------------------------------------------------------------------
# 27. PE / BSS incomplete-trip arm alignment (2026-09-02). The BSS learns R_G_boat from a
#     frame that HAS been incomplete-trip filtered; both Point Estimators used the unfiltered
#     one. sensitivity_incomplete_trips.csv has reported that disagreement since 2026-08-25.
#     One shared helper, because the gear PE's comment already claimed it matched the BSS
#     while it did not, which is what two copies of one rule produces.
# ---------------------------------------------------------------------------
local({
  chk("pe arm: the shared frame helper exists",
      exists("pe_gear_ratio_frame", mode = "function"))
  iv <- data.frame(number_of_gear = c(3, 4, 5, 6), angler_count = c(1, 1, 1, 1),
                   trip_status = c("Complete", "Incomplete", NA, "Complete"),
                   stringsAsFactors = FALSE)
  p_on  <- list(pe_gear_ratio_arm = "match_bss", filter_incomplete_trips = TRUE)
  p_off <- list(pe_gear_ratio_arm = "gear_only", filter_incomplete_trips = TRUE)
  chk("pe arm: match_bss drops incomplete trips",
      nrow(pe_gear_ratio_frame(iv, NULL, p_on, quiet = TRUE)) == 3)
  chk("pe arm: and KEEPS a missing status, exactly as prep_bss_crab_pooled.R does",
      NA %in% pe_gear_ratio_frame(iv, NULL, p_on, quiet = TRUE)$trip_status)
  chk("pe arm: gear_only is the untouched pre-2026-09-02 behaviour",
      nrow(pe_gear_ratio_frame(iv, NULL, p_off, quiet = TRUE)) == 4)
  chk("pe arm: with filter_incomplete_trips OFF there is nothing to align",
      nrow(pe_gear_ratio_frame(iv, NULL, list(pe_gear_ratio_arm = "match_bss",
                                              filter_incomplete_trips = FALSE), quiet = TRUE)) == 4)
  chk("pe arm: a caller-supplied frame is left alone (the four-arm diagnostic relies on this)",
      nrow(pe_gear_ratio_frame(iv[1:2, ], iv, p_on, quiet = TRUE)) == 4)
  chk("pe arm: an unrecognised arm ERRORS rather than silently defaulting",
      inherits(try(pe_gear_ratio_frame(iv, NULL, list(pe_gear_ratio_arm = "nope"), quiet = TRUE),
                   silent = TRUE), "try-error"))
  for (fn in c("03_R_functions/run_pe_pooled.R", "03_R_functions/run_pe_gear.R")) {
    t <- paste(readLines(fn, warn = FALSE), collapse = "\n")
    chk(sprintf("pe arm: %s routes through the shared helper", basename(fn)),
        grepl("pe_gear_ratio_frame(", t, fixed = TRUE))
  }
  e <- new.env(); sys.source("run_config.R", envir = e)
  chk("pe arm: run_config ships the aligned arm",
      identical(e$run_config$pe_gear_ratio_arm, "match_bss"))
})

# ---------------------------------------------------------------------------
# 28. Zero-inflated catch likelihood (2026-09-02; ADOPTED 2026-09-07, ships ON for shore). Targeted from the
#     2026-09-01 zero bin: shore all-gear 676 observed zeros vs 605.4 expected (z = +3.8),
#     shore pot closure z = +2.9, every boat stream inside |z| = 2.3.
# ---------------------------------------------------------------------------
local({
  t <- paste(readLines("02_stan_models/crab_bss_pooled.stan", warn = FALSE), collapse = "\n")
  chk("zinb: the guard flag and its Beta prior are Stan data",
      grepl("int<lower=0, upper=1> zi_catch;", t, fixed = TRUE) &&
      grepl("real<lower=0> zi_catch_prior_a;", t, fixed = TRUE))
  chk("zinb: theta_C is ZERO-SIZE when off, so the OFF path is bit-identical",
      grepl("vector<lower=0, upper=1>[zi_catch] theta_C;", t, fixed = TRUE))
  chk("zinb: the model block uses a log_mix at the observed zeros",
      grepl("target += log_mix(theta_C[1], 0, neg_binomial_2_lpmf(0 | mu_c, r_C));", t, fixed = TRUE))
  # If log_lik does not mirror the fitted likelihood, every elpd_loo comparison is void, and
  # an elpd comparison is the entire basis on which this feature will be judged.
  chk("zinb: generated-quantities log_lik MIRRORS the fitted likelihood",
      grepl("log_lik_catch[a] = log_mix(theta_C[1], 0, neg_binomial_2_lpmf(0 | mu_c_gq, r_C));",
            t, fixed = TRUE))
  # The correctness point that would otherwise silently inflate the headline.
  chk("zinb: the season total is scaled by (1 - theta_C)",
      # D40 (2026-10-04): the day's share is fd (f_crab[f_stratum[d]] when the volume term is off)
      grepl("* fd * zi_scale;", t, fixed = TRUE) && grepl("real fd = f_crab[f_stratum[d]];", t, fixed = TRUE) &&
      grepl("zi_scale = 1 - theta_C[1];", t, fixed = TRUE) &&
      grepl("zi_scale = 1.0;", t, fixed = TRUE))
  chk("zinb: theta_C_out is reported unconditionally so the parameter set keeps its shape",
      grepl("real theta_C_out;", t, fixed = TRUE))
  d <- paste(readLines("03_R_functions/model_diagnostics.R", warn = FALSE), collapse = "\n")
  chk("zinb: a hard 0 is flagged decoupled, so it cannot be read as 'tested and absent'",
      grepl('set(base %in% c("theta_C", "theta_C_out")', d, fixed = TRUE))
  chk("zinb: theta_C reaches the curated structural table",
      grepl('"theta_C")', d, fixed = TRUE))
  p <- paste(readLines("03_R_functions/prep_bss_crab_pooled.R", warn = FALSE), collapse = "\n")
  chk("zinb: it is scoped PER FIT, which is what makes the boat a free negative control",
      grepl("population_name %in% (params$catch_zi_populations %||% \"shore\")", p, fixed = TRUE))
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  # 2026-09-07: this asserted the PROTOTYPE default (OFF). The feature was adopted on
  # 2026-09-07 after the C1 run, so the shipped default is now TRUE. The scoping and the
  # prior are the parts that must not drift; whether it is on is a decision, recorded in
  # section 42 below and in PIPELINE_STATUS.md Section 1k.
  chk("zinb: scoped to the shore, with a Beta(1, 9) prior",
      identical(rc$catch_zi_populations, "shore") &&
      identical(rc$zi_catch_prior_a, 1) && identical(rc$zi_catch_prior_b, 9))
})

# ---------------------------------------------------------------------------
# 29. The 2026-09-03 batch runner.
# ---------------------------------------------------------------------------
local({
  f <- "06_diagnostics/run_shore_ar_zi_2026-09-03.R"
  if (!file.exists(f)) { chk("shore/zi runner: present", FALSE); return(invisible(NULL)) }
  chk("shore/zi runner: parses", !inherits(try(parse(f), silent = TRUE), "try-error"))
  t <- paste(readLines(f, warn = FALSE), collapse = "\n")
  chk("shore/zi runner: ships with DRY_RUN TRUE", grepl("\nDRY_RUN <- TRUE", t))
  chk("shore/zi runner: the AR ladder forces ONE sub-season of ONE population",
      grepl("AR_SHORE_AG <- function(res) list(shore = list(all_gear = res))", t, fixed = TRUE))
  # 2026-09-04: the three forced rungs were replaced by ONE run of the production
  # escalation toggle, which fits every rung inside a single render and logs each rung's own
  # estimate. Fewer fits, one output folder, and it exercises the mechanism a future season
  # will actually use instead of an experiment-only override.
  chk("shore/zi runner: the ladder runs via the PRODUCTION toggle, not a forced override",
      grepl('ar_escalate = list(shore = "all_gear")', t, fixed = TRUE) &&
      grepl('ar_escalate_stop = "all_rungs"', t, fixed = TRUE))
  chk("shore/zi runner: the ladder table is built from ar_escalation_log.csv",
      grepl('rd(d, "ar_escalation_log.csv")', t, fixed = TRUE))
  chk("shore/zi runner: Z0 remains the OFF gate and the production daily reference",
      grepl('lad_row("daily (production)", "daily", sid = "Z0")', t, fixed = TRUE))
  chk("shore/zi runner: coverage is read from the randomized PIT, not ppc_calibration",
      grepl("cov50_byobs <- function", t, fixed = TRUE))
  chk("shore/zi runner: it refuses to run if the two code changes are absent",
      grepl("has nothing to test", t, fixed = TRUE) &&
      grepl("would reproduce the ", t, fixed = TRUE))
  chk("shore/zi runner: expect_delta carried over from the 2026-09-01 batch",
      grepl("expect_delta", t, fixed = TRUE) && grepl("config_delta <- function", t, fixed = TRUE))
  # Drift guard on the production numbers the ladder is judged against.
  v1 <- "05_output/20260831/pooled-CPUE-VAL-1-adopted"
  if (!dir.exists(v1)) skp("shore/zi REF: production adequacy", "folder absent")
  if (dir.exists(v1)) {
    a <- read.csv(file.path(v1, "model_adequacy.csv"), stringsAsFactors = FALSE)
    r <- a[a$fit == "shore_all_gear_Dungeness_Kept", ]
    chk("shore/zi REF: production shore all-gear is still the flagged fit",
        isTRUE(round(r$p_loo_frac, 2) == 0.35) && isTRUE(r$n_pareto_bad == 41),
        sprintf("(p_loo %.3f, k>0.7 %d)", r$p_loo_frac, r$n_pareto_bad))
    x <- read.csv(file.path(v1, "ppc_byobs_shore_all_gear_Dungeness_Kept.csv"), stringsAsFactors = FALSE)
    y <- x[x$data_type == "gear", ]
    cov <- mean(as.logical(y$in_50)); z <- (cov - 0.5) / sqrt(0.25 / nrow(y))
    chk("shore/zi REF: and its gear-stream coverage is still ~0.70 at +7 sampling SDs",
        isTRUE(abs(cov - 0.701) < 0.01) && isTRUE(z > 6),
        sprintf("(%.3f, %+.1f SD on n = %d)", cov, z, nrow(y)))
    zz <- x[x$data_type == "catch", ]
    p <- zz$p_zero[is.finite(zz$p_zero)]
    zsc <- (sum(zz$observed == 0) - sum(p)) / sqrt(sum(p * (1 - p)))
    chk("shore/zi REF: the shore catch zero bin the ZINB targets is still off at z ~ +3.8",
        isTRUE(zsc > 3), sprintf("(%d observed vs %.1f expected, z = %+.1f)",
                                 sum(zz$observed == 0), sum(p), zsc))
  }
})

# ---------------------------------------------------------------------------
# 30. EVERY batch runner ships DRY_RUN <- TRUE (2026-09-03).
#     The harness asserted this for the three most recent runners only, which is exactly why
#     run_patch_validation_2026-08-25.R sat at FALSE through 264 passing assertions: sourcing
#     it started real fits and appended duplicate summary rows, while both its own inline
#     comment and 06_diagnostics/README.md said it would fit nothing. An assertion that
#     covers a hand-picked subset is an assertion that will drift.
# ---------------------------------------------------------------------------
local({
  runners <- list.files("06_diagnostics", pattern = "^run_.*\\.R$", full.names = TRUE)
  # run_rg_sweep / run_tau_sweep / run_osp_validation have no dry-run mode at all; they are
  # checked separately below rather than silently exempted.
  has_dry <- vapply(runners, function(p) any(grepl("^DRY_RUN <-", readLines(p, warn = FALSE))), logical(1))
  for (p in runners[has_dry]) {
    l <- readLines(p, warn = FALSE); v <- sub("^DRY_RUN <- *([A-Z]+).*$", "\\1", grep("^DRY_RUN <-", l, value = TRUE)[1])
    chk(sprintf("dry-run default: %s ships TRUE", basename(p)), identical(v, "TRUE"), sprintf("(%s)", v))
  }
  chk("dry-run default: every runner with a DRY_RUN was checked, not a hand-picked subset",
      sum(has_dry) >= 5, sprintf("(%d of %d runners have a DRY_RUN switch)", sum(has_dry), length(runners)))
  # The runners with no dry-run mode start fitting on source. That is a real hazard and it is
  # recorded here rather than fixed, because each is a small single-purpose sweep whose header
  # says so; the assertion exists so the list cannot grow unnoticed.
  no_dry <- basename(runners[!has_dry])
  # 2026-09-28 (B49): run_rg_sweep.R gained its DRY_RUN switch, so the set shrank to the two
  # superseded runners, which refuse to fit.
  chk("dry-run default: the set of runners with NO dry-run mode has not grown",
      setequal(no_dry, c("run_osp_validation.R", "run_tau_sweep.R")),
      sprintf("(%s)", paste(no_dry, collapse = ", ")))
})

# ---------------------------------------------------------------------------
# 31. Result-file persistence (2026-09-03). Two ways of getting this wrong have already
#     destroyed results: append-on-resume duplicated summary rows, and an ungated overwrite
#     truncated 05_output/validation_2026-09-01_verdicts.csv to its desk rows, losing all
#     seven fitted-stage criteria while the fits themselves survived.
# ---------------------------------------------------------------------------
local({
  f <- "06_diagnostics/run_shore_ar_zi_2026-09-03.R"
  if (!file.exists(f)) { chk("persistence: runner present", FALSE); return(invisible(NULL)) }
  t <- paste(readLines(f, warn = FALSE), collapse = "\n")
  chk("persistence: rows are MERGED BY KEY, not appended",
      grepl("merge_csv_by <- function(new, path, key)", t, fixed = TRUE) &&
      grepl('append_row <- function(r) merge_csv_by(r, sum_path, "stage")', t, fixed = TRUE))
  chk("persistence: verdicts merge on (stage, criterion), so a dry run cannot erase fitted rows",
      grepl('merge_csv_by(VD, ver_path, c("stage", "criterion"))', t, fixed = TRUE))
  chk("persistence: the ladder MERGES its fallback rung instead of substituting it",
      grepl('merge_csv_by(LAD, lad_path, "rung")', t, fixed = TRUE) &&
      !grepl("utils::write.csv(LAD, lad_path", t, fixed = TRUE))
  # Behavioural check on the merge itself, not just its presence.
  # Extract merge_csv_by from the runner by brace-matching from its definition line, rather
  # than by regex: the point is to test the REAL function, not a copy of it that could drift.
  .rl <- readLines(f, warn = FALSE)
  .st <- grep("^merge_csv_by <- function", .rl)[1]
  .en <- .st + which(.rl[.st:length(.rl)] == "}")[1] - 1L
  e <- new.env(); eval(parse(text = paste(.rl[.st:.en], collapse = "\n")), envir = e)
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(stage = c("A1", "A2"), v = c(1, 2)), tmp, row.names = FALSE)
  e$merge_csv_by(data.frame(stage = "A1", v = 99), tmp, "stage")
  got <- read.csv(tmp, stringsAsFactors = FALSE)
  chk("persistence: re-running one stage REPLACES its row and keeps the others",
      nrow(got) == 2 && isTRUE(got$v[got$stage == "A1"] == 99) && isTRUE(got$v[got$stage == "A2"] == 2),
      sprintf("(%d rows: %s)", nrow(got), paste(got$stage, got$v, sep = "=", collapse = ", ")))
  unlink(tmp)
})

# ---------------------------------------------------------------------------
# 32. The shore AR ladder's controls (2026-09-03).
# ---------------------------------------------------------------------------
local({
  f <- "06_diagnostics/run_shore_ar_zi_2026-09-03.R"
  t <- paste(readLines(f, warn = FALSE), collapse = "\n")
  chk("ladder: an ar_force LEAK control runs on every rung",
      grepl('ex_pc <- fit_exactness(nd, z0, "shore_ring_net_only"', t, fixed = TRUE) &&
      grepl('ex_bt <- fit_exactness(nd, z0, "private_boat"', t, fixed = TRUE))
  chk("ladder: a rung whose gate failed is flagged rather than ranked on catch",
      grepl("every rung actually reported its BSS", t, fixed = TRUE))
  chk("ladder: the coverage n is read from the run, not hard-coded",
      grepl('sum(x$data_type == "gear")', t, fixed = TRUE) && !grepl("n_gear <- 311", t, fixed = TRUE))
  chk("ladder: the ZINB elpd threshold is in SE units, not bare nats",
      grepl("shore_ag_se_elpd_catch", t, fixed = TRUE) && grepl("2 * .se", t, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 33. The AR escalation ladder as a PRODUCTION toggle (2026-09-04). It came out of the FW
#     creel team meeting: a component should report at the finest AR resolution its own
#     sampler behaviour supports, decided per run rather than frozen into config.
# ---------------------------------------------------------------------------
local({
  chk("ar_escalate: the per-population/sub-season resolver exists",
      exists(".bss_resolve_ar_escalate", mode = "function"))
  r <- .bss_resolve_ar_escalate
  chk("ar_escalate: FALSE and absent are off",
      isFALSE(r(list(ar_escalate = FALSE), "shore", "all_gear")) && isFALSE(r(list(), "shore", "all_gear")))
  chk("ar_escalate: TRUE is on for every fit",
      isTRUE(r(list(ar_escalate = TRUE), "shore", "all_gear")) &&
      isTRUE(r(list(ar_escalate = TRUE), "private_boat", "pot_closure")))
  chk("ar_escalate: a character vector scopes to POPULATIONS",
      isTRUE(r(list(ar_escalate = "shore"), "shore", "all_gear")) &&
      isFALSE(r(list(ar_escalate = "shore"), "private_boat", "all_gear")))
  chk("ar_escalate: a named list scopes to population x SUB-SEASON",
      isTRUE(r(list(ar_escalate = list(shore = "all_gear")), "shore", "all_gear")) &&
      isFALSE(r(list(ar_escalate = list(shore = "all_gear")), "shore", "pot_closure")) &&
      isFALSE(r(list(ar_escalate = list(shore = "all_gear")), "private_boat", "all_gear")))
  chk("ar_escalate: an unusable shape ERRORS rather than silently defaulting off",
      inherits(try(r(list(ar_escalate = 3), "shore", "all_gear"), silent = TRUE), "try-error"))

  chk("ar rung summary: the per-rung estimate helper exists",
      exists("bss_ar_rung_summary", mode = "function"))
  chk("ar rung summary: a NULL fit yields NAs rather than erroring",
      { v <- try(bss_ar_rung_summary(NULL), silent = TRUE)
        !inherits(v, "try-error") && all(is.na(unlist(v))) })

  for (rmd in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd",
                "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")) {
    t <- paste(readLines(rmd, warn = FALSE), collapse = "\n")
    b <- basename(rmd)
    chk(sprintf("ar ladder: %s logs each rung's OWN estimate and interval", b),
        grepl("catch_median       = .rung$catch_median", t, fixed = TRUE) ||
        grepl("catch_median        = .rung$catch_median", t, fixed = TRUE))
    chk(sprintf("ar ladder: %s resolves the escalation scope PER FIT", b),
        grepl(".esc_on <- .bss_resolve_ar_escalate(params, pop, ss$gear_regime)", t, fixed = TRUE))
    chk(sprintf("ar ladder: %s honours ar_escalate_stop", b),
        grepl('params$ar_escalate_stop   %||% "first_pass"', t, fixed = TRUE) &&
        grepl('identical(.stop_rule, "all_rungs")', t, fixed = TRUE))
    chk(sprintf("ar ladder: %s never keeps a rung that FAILED the gate under all_rungs", b),
        grepl("else if (!.passed) FALSE", t, fixed = TRUE))
  }
  t <- paste(readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE), collapse = "\n")
  chk("ar ladder: the report table shows each rung's estimate and marks the reported one",
      grepl("`Catch (median)` = round(catch_median)", t, fixed = TRUE) &&
      grepl("Reported = selected", t, fixed = TRUE))
  chk("ar ladder: the report warns that a narrower interval is not evidence of a better model",
      grepl("narrower interval is not by itself evidence", t, fixed = TRUE))

  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("ar ladder: production ships the ladder OFF and stopping at the first pass",
      isFALSE(rc$ar_escalate) && identical(rc$ar_escalate_stop, "first_pass") &&
      identical(rc$ar_escalate_select, "first_pass"))
  chk("ar ladder: run_config records WHY narrowest_pi is not the default",
      { s <- paste(readLines("run_config.R", warn = FALSE), collapse = "\n")
        grepl("the two MISCALIBRATED cells look the most precise", s, fixed = TRUE) ||
        grepl("miscalibrated cells look the most precise", s, fixed = TRUE) })
})

# ---------------------------------------------------------------------------
# 34. Project context is recorded where an operator or an agent will actually meet it.
# ---------------------------------------------------------------------------
local({
  rc <- paste(readLines("run_config.R", warn = FALSE), collapse = "\n")
  chk("context: run_config states that nothing has been published",
      grepl("NOTHING HAS BEEN PUBLISHED", rc, fixed = TRUE))
  chk("context: run_config states what the OSP crab-only column IS (a lower bound on f)",
      grepl("LOWER BOUND", rc, fixed = TRUE) && grepl("crab-only column in as if it were f", rc, fixed = TRUE))
  cl <- paste(readLines("07_documentation/CLAUDE.md", warn = FALSE), collapse = "\n")
  chk("context: CLAUDE.md carries the same two facts",
      grepl("published \\*\\*no\\*\\*", cl) && grepl("lower bound", cl, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 35. The 2026-09-03 defects. Each of these is one line of code and each cost a real
#     result: 14.7 h of wasted fitting, a nearly-rejected feature, and five of eleven
#     batch verdicts wrong. Assert the FIX, not the symptom.
# ---------------------------------------------------------------------------
local({
  # 35a. ar_escalate is scoped-resolvable EVERYWHERE it is consulted. The ladder built
  #      four rungs and the loop ran four attempts, but force_resolution was gated on
  #      isTRUE(params$ar_escalate), which is FALSE for list(shore = "all_gear"), so the
  #      data prep re-derived `daily` on every rung. Assert no driver tests the raw value.
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    src <- readLines(drv, warn = FALSE)
    code <- src[!grepl("^\\s*#", src)]        # comments may quote the old expression
    chk(sprintf("%s: no isTRUE(params$ar_escalate) in live code", basename(drv)),
        !any(grepl("isTRUE(params$ar_escalate)", code, fixed = TRUE)),
        paste(grep("isTRUE(params$ar_escalate)", code, fixed = TRUE, value = TRUE), collapse = " | "))
    chk(sprintf("%s: force_resolution is gated on the resolved flag", basename(drv)),
        any(grepl("force_res <- if (.esc_on)", code, fixed = TRUE)))
    # 35b. the reported rung is identified by ATTEMPT INDEX, never by resolution string.
    chk(sprintf("%s: selected rung tracked by attempt index", basename(drv)),
        any(grepl(".kept_attempt[[label]] <- attempt", code, fixed = TRUE)))
  }
  # 35c. a scoped ar_escalate must still produce a MULTI-RUNG ladder with DISTINCT rungs.
  #      This is the property the A1 stage silently lost: four rungs, one resolution.
  Ps <- modifyList(P1, list(ar_escalate = list(shore = "all_gear")))
  Ls <- bss_ar_ladder(days289, eff, "shore", Ps, gear_regime = "all_gear")
  chk("scoped ar_escalate yields a multi-rung ladder", length(Ls) >= 3, paste(Ls, collapse = "->"))
  chk("every rung on the ladder is a DISTINCT resolution", !any(duplicated(Ls)), paste(Ls, collapse = "->"))
  chk("scoped ar_escalate does NOT reach the unscoped sub-season",
      identical(bss_ar_ladder(days76, tibble(day_index = seq(1, 76, by = 2)), "shore", Ps,
                              gear_regime = "pot_closure"), "daily") ||
      length(bss_ar_ladder(days76, tibble(day_index = seq(1, 76, by = 2)), "shore", Ps,
                           gear_regime = "pot_closure")) == 1)
  chk("scoped ar_escalate does NOT reach the other population",
      length(bss_ar_ladder(days289, eff, "private_boat", Ps, gear_regime = "all_gear")) == 1)
})

# ---------------------------------------------------------------------------
# 36. The catch-stream PPC is computed under the likelihood Stan actually sampled.
#     Under zi_catch = 1 the NB2 arithmetic reported 419 expected zeros against 676
#     observed (z = 15.4) where the mixture gives 637.9 (z = 2.0), and the batch
#     recorded that as the feature making the zero bin worse.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/zinb_ppc.R")
  th <- 0.178; mu <- 3.0; sz <- 0.8
  f <- function(y) if (y == 0) th + (1 - th) * dnbinom(0, size = sz, mu = mu) else
                               (1 - th) * dnbinom(y, size = sz, mu = mu)
  chk("mixture pmf integrates to 1", abs(sum(vapply(0:500, f, 0)) - 1) < 1e-10)
  chk("mixture p_zero matches the pmf at 0", abs(bss_zi_p_zero(mu, sz, th) - f(0)) < 1e-12)
  # THE branch that is easy to get wrong: at y = 0, F(-1) = 0, so the point mass theta is
  # NOT added ahead of the observation. Folding the two cases into one expression puts the
  # zero observations about theta too high, flattering exactly the observations theta explains.
  chk("mixture PIT at y = 0 is 0.5 * f(0), not theta + ...",
      abs(bss_zi_pit(0, mu, sz, th) - 0.5 * f(0)) < 1e-12)
  chk("mixture PIT at y = 0 is NOT the naive folded form",
      abs(bss_zi_pit(0, mu, sz, th) - (th + (1 - th) * 0.5 * dnbinom(0, size = sz, mu = mu))) > 1e-3)
  chk("mixture PIT at y > 0 = F(y-1) + 0.5 f(y)",
      abs(bss_zi_pit(4, mu, sz, th) - (sum(vapply(0:3, f, 0)) + 0.5 * f(4))) < 1e-12)
  # 2026-09-05: the ONE bin. The zero bin passed on the prototype while y=1 lost -42 nats;
  # a count-bin check must cover both, and the mixture P(Y=1) is (1-theta)*NB2(1), NOT
  # theta + (1-theta)*NB2(1).
  chk("mixture p_k(1) = (1-theta) * NB2(1)", abs(bss_zi_p_k(1, mu, sz, th) - f(1)) < 1e-12)
  chk("mixture p_k(0) agrees with p_zero", abs(bss_zi_p_k(0, mu, sz, th) - bss_zi_p_zero(mu, sz, th)) < 1e-12)
  chk("p_one is written beside p_zero",
      any(grepl("p_one[i]  <- bss_zi_p_k(1L, mu, sz, th)",
                readLines("03_R_functions/save_run_diagnostics.R", warn = FALSE), fixed = TRUE)))
  chk("theta = NULL reproduces the NB2 arithmetic exactly",
      identical(bss_zi_pit(3, mu, sz, NULL),
                mean(pnbinom(2, size = sz, mu = mu) + 0.5 * dnbinom(3, size = sz, mu = mu))) &&
      identical(bss_zi_p_zero(mu, sz, NULL), mean(dnbinom(0, size = sz, mu = mu))))
  # both diagnostic files must route through the shared helper, or they drift apart again
  for (f2 in c("03_R_functions/save_run_diagnostics.R", "03_R_functions/model_diagnostics.R")) {
    src <- readLines(f2, warn = FALSE)
    chk(sprintf("%s uses the shared mixture PIT", basename(f2)),
        any(grepl("bss_zi_pit(", src, fixed = TRUE)))
    chk(sprintf("%s passes theta on the CATCH stream only", basename(f2)),
        any(grepl("bss_zi_theta_draws(fit, stan_data,", src, fixed = TRUE)))
  }
})

# ---------------------------------------------------------------------------
# 37. Model comparison uses the PAIRED difference SE, and the config dump is complete.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/loo_elpd_paired.R")
  set.seed(1)
  n <- 800; base <- rnorm(n, -2, 1.5)          # large across-observation spread
  a <- data.frame(obs_index = 1:n, observed = rpois(n, 1), elpd_loo = base,
                  pareto_k = runif(n, 0, .5))
  b <- a; b$elpd_loo <- base + rnorm(n, 0.02, 0.05)   # small, consistent per-obs gain
  ta <- tempfile(fileext = ".csv"); tb <- tempfile(fileext = ".csv")
  write.csv(a, ta, row.names = FALSE); write.csv(b, tb, row.names = FALSE)
  r <- loo_elpd_paired(ta, tb, "synthetic")
  chk("paired SE is far smaller than the per-fit SE when the fits are correlated",
      r$se_diff < r$se_naive_a / 10, sprintf("%.2f vs %.2f", r$se_diff, r$se_naive_a))
  chk("paired ratio detects a real small effect the naive ratio would miss",
      r$ratio > 5 && abs(r$elpd_diff / r$se_naive_a) < 1,
      sprintf("paired %.1f SE, naive %.2f SE", r$ratio, r$elpd_diff / r$se_naive_a))
  chk("by-count decomposition is returned and sums to the total",
      !is.null(r$by_count) && abs(sum(r$by_count[, "diff"], na.rm = TRUE) - r$elpd_diff) < 1e-8)
  chk("misaligned obs_index is refused rather than silently compared",
      is.null(local({ b2 <- b; b2$observed <- rev(b2$observed)
                      t2 <- tempfile(fileext = ".csv"); write.csv(b2, t2, row.names = FALSE)
                      r2 <- loo_elpd_paired(ta, t2); if (identical(b2$observed, a$observed)) NULL else r2 })))
  # run_parameters.txt must capture EVERY run_config key: config_delta() and
  # annotate_decoupled_run() both read it, and str()'s default caps a list at 99.
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE))
    chk(sprintf("%s writes an untruncated config dump", basename(drv)),
        any(grepl("str(params, list.len =", readLines(drv, warn = FALSE), fixed = TRUE)))
  src <- readLines("run_config.R", warn = FALSE)
  dump <- capture.output(str(run_config_for_test <- local({ e <- new.env(); sys.source("run_config.R", e); e$run_config }),
                            list.len = 1000, nchar.max = 2000, vec.len = 100))
  chk("the dump settings capture every run_config key",
      !any(grepl("truncated", dump)) &&
      sum(grepl("^ \\$ ", dump)) == length(run_config_for_test),
      sprintf("%d of %d", sum(grepl("^ \\$ ", dump)), length(run_config_for_test)))
  chk("annotate_decoupled_run refuses to trust a truncated dump",
      any(grepl(".adr_truncated", readLines("03_R_functions/annotate_decoupled_run.R", warn = FALSE), fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 38. Every batch runner hands the driver its folder name the ONE way the driver reads
#     it. BSS-GH-*-CPUE-model.Rmd builds output_dir from run_config$run_tag and ignores
#     both a bare `run_tag` variable and an `output_dir=` passed to render(). The first
#     version of run_ladder_zinb_2026-09-04.R did both of the wrong things, so its two
#     fitted stages would have shared the run_config.R default folder and Z2 would have
#     overwritten the L1 ladder. Caught in the 2026-09-05 pre-run audit, not by a dry run,
#     because a dry run returns before render(). Generated over every runner.
# ---------------------------------------------------------------------------
local({
  for (rf in list.files("06_diagnostics", pattern = "^run_.*\\.R$", full.names = TRUE)) {
    src <- readLines(rf, warn = FALSE); src <- src[!grepl("^\\s*#", src)]
    calls <- src[grepl("rmarkdown::render(", src, fixed = TRUE)]
    if (!length(calls)) next
    chk(sprintf("%s: render() is not handed output_dir= (the driver ignores it)", basename(rf)),
        !any(grepl("output_dir", calls, fixed = TRUE)), paste(trimws(calls), collapse = " ; "))
    chk(sprintf("%s: run_tag is set INSIDE the config the driver reads", basename(rf)),
        any(grepl("\\$run_tag\\s*<-", src)))
  }
  # and the drivers must keep reading it from there
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE))
    chk(sprintf("%s: output_dir is built from run_config$run_tag", basename(drv)),
        any(grepl("run_config$run_tag", readLines(drv, warn = FALSE), fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 39. The ladder records adequacy PER RUNG. On the 2026-09-04 run all three rungs
#     PASSED the convergence gate, so the gate separated nothing, and the statistics
#     that could separate them (p_loo as a fraction of n_obs, the Pareto k count,
#     coverage) existed for the ONE rung the loop kept. The two coarser rungs cost
#     174 minutes of fitting and left no evidence.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/bss_rung_adequacy.R")
  chk("rung adequacy: NULL fit yields an all-NA row, never an error",
      all(is.na(unlist(bss_rung_adequacy(NULL, NULL)))))
  nm <- names(bss_rung_adequacy(NULL, NULL))
  chk("rung adequacy: carries p_loo_frac, the Pareto count and per-stream coverage",
      all(c("p_loo_frac", "n_pareto_bad", "cov50_gear", "cov50_catch") %in% nm), paste(nm, collapse = ","))
  # it must call loo the SAME way write_loo_diagnostics does, or its numbers are not
  # comparable with loo_summary_*.csv and model_adequacy.csv for the reported rung
  src <- readLines("03_R_functions/bss_rung_adequacy.R", warn = FALSE)
  # 2026-09-07: the two defects the 2026-09-04 C2 run exposed.
  chk("rung adequacy: p_loo_frac is the WORST STREAM, as bss_model_adequacy defines it",
      any(grepl("which.max(replace(frac", readLines("03_R_functions/bss_rung_adequacy.R", warn = FALSE), fixed = TRUE)) &&
      "p_loo_worst_stream" %in% names(bss_rung_adequacy(NULL, NULL)))
  chk("rung adequacy: restores the RNG state, so a diagnostic toggle cannot move a diagnostic",
      { src <- readLines("03_R_functions/bss_rung_adequacy.R", warn = FALSE)
        any(grepl(".Random.seed", src, fixed = TRUE)) && any(grepl("on.exit(", src, fixed = TRUE)) })
  chk("rung adequacy: uses the production loo call (plain matrix, no r_eff)",
      any(grepl("loo::loo(ll)", src, fixed = TRUE)) &&
      !any(grepl("relative_eff", src, fixed = TRUE)))
  chk("rung adequacy: guards absent streams on the same n as production",
      any(grepl("stan_data$Gear_n", src, fixed = TRUE)) &&
      any(grepl("stan_data$IntC", src, fixed = TRUE)))
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    chk(sprintf("%s: logs per-rung adequacy", basename(drv)),
        any(grepl("bss_rung_adequacy(fit_try, bss_data_try)", d, fixed = TRUE)) &&
        any(grepl("p_loo_frac", d, fixed = TRUE)))
  }
})

# ---------------------------------------------------------------------------
# 40. A batch runner that renders more than once must MOVE the HTML. rmarkdown writes
#     it beside the .Rmd, not into output_dir; on 2026-09-04 stage L1 rendered its
#     report there and stage Z2 overwrote it, destroying the only rendered view of the
#     AR ladder and leaving a 6,000-line artefact committed in 01_BSS_models/.
# ---------------------------------------------------------------------------
local({
  for (rf in list.files("06_diagnostics", pattern = "^run_.*\\.R$", full.names = TRUE)) {
    src <- readLines(rf, warn = FALSE); src <- src[!grepl("^\\s*#", src)]
    if (!any(grepl("rmarkdown::render(", src, fixed = TRUE))) next
    chk(sprintf("%s: moves the rendered HTML into the run folder", basename(rf)),
        any(grepl("file.copy(html", src, fixed = TRUE)))
  }
  chk("01_BSS_models/*.html is gitignored so a stray render is never committed",
      any(grepl("^01_BSS_models/\\*\\.html", readLines(".gitignore", warn = FALSE))))
  chk("no rendered driver HTML is committed beside the .Rmd",
      !length(list.files("01_BSS_models", pattern = "\\.html$")))
})

# ---------------------------------------------------------------------------
# 41. A verdicts file must never carry two rows for one key. merge_csv_by() dropped
#     stale rows the new frame supersedes but never de-duplicated the NEW frame, so
#     when the 2026-09-04 batch ran verdict_L1() for both L1 and C2 under one stage
#     label the file ended with two rows per criterion.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/batch_verdict_helpers.R")
  tf <- tempfile(fileext = ".csv")
  d <- data.frame(stage = c("A", "A", "B"), criterion = c("x", "x", "y"),
                  observed = c("first", "second", "z"), stringsAsFactors = FALSE)
  merge_csv_by(d, tf, c("stage", "criterion"))
  got <- utils::read.csv(tf, stringsAsFactors = FALSE)
  chk("merge_csv_by de-duplicates within the new frame", nrow(got) == 2, sprintf("%d rows", nrow(got)))
  chk("merge_csv_by keeps the LAST write for a duplicated key",
      identical(got$observed[got$criterion == "x"], "second"))
  # and the runner must label a reused verdict block with the stage that produced it
  src <- readLines("06_diagnostics/run_ladder_zinb_2026-09-04.R", warn = FALSE)
  src <- src[!grepl("^\\s*#", src)]
  chk("verdict_L1 takes a stage id so a second ladder is not filed under the first",
      any(grepl("verdict_L1 <- function(dir, stage_id", src, fixed = TRUE)) &&
      any(grepl('verdict_L1(dirs$C2, "C2")', src, fixed = TRUE)))
  # the committed verdicts file itself must be clean
  vp <- "05_output/ladder_zinb_2026-09-04_verdicts.csv"
  if (file.exists(vp)) {
    v <- utils::read.csv(vp, stringsAsFactors = FALSE)
    k <- paste(v$stage, v$criterion, sep = "\r")
    chk("the committed verdicts file has one row per stage+criterion",
        !any(duplicated(k)), paste(unique(k[duplicated(k)]), collapse = " | "))
  }
})

# ---------------------------------------------------------------------------
# 42. The 2026-09-07 adoption. run_config must SHIP the candidate configuration, with
#     no experiment lever active, and the gear track must be untouched so the
#     cross-check measures the pooled change alone.
# ---------------------------------------------------------------------------
local({
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("adoption: pooled shore all-gear is capped at weekly",
      identical(rc$ar_max_resolution$pooled$shore$all_gear, "weekly"))
  chk("adoption: the ZINB catch likelihood is on, scoped to shore",
      isTRUE(rc$estimate_catch_zi) && identical(rc$catch_zi_populations, "shore"))
  chk("adoption: no experiment lever ships enabled",
      is.null(rc$ar_force) && !isTRUE(rc$ar_escalate))
  chk("adoption: the other pooled caps and the whole gear map are untouched",
      identical(rc$ar_max_resolution$pooled$shore$pot_closure, "biweekly") &&
      identical(rc$ar_max_resolution$pooled$private_boat, "monthly") &&
      identical(rc$ar_max_resolution$gear_resolved$shore$all_gear, "monthly"))
  # The asymmetry this creates is a fact about the models, so assert it is RECORDED
  # rather than assert it away: the gear Stan has no ZINB and the config says so.
  gz <- any(grepl("zi_catch", readLines("02_stan_models/crab_bss_gear_resolved.stan", warn = FALSE), fixed = TRUE))
  rc <- paste(readLines("run_config.R", warn = FALSE), collapse = "\n")
  chk("adoption: if the gear model lacks the ZINB, run_config says so where the flag lives",
      gz || grepl("POOLED ONLY", rc, fixed = TRUE),
      sprintf("gear stan has zi_catch: %s", gz))
  # and the runner that renders it carries the two standing runner rules
  rf <- "06_diagnostics/run_adoption_2026-09-07.R"
  if (file.exists(rf)) {
    src <- readLines(rf, warn = FALSE); src <- src[!grepl("^\\s*#", src)]
    calls <- src[grepl("rmarkdown::render(", src, fixed = TRUE)]
    chk("adoption runner: render() is not handed output_dir=",
        length(calls) > 0 && !any(grepl("output_dir", calls, fixed = TRUE)))
    chk("adoption runner: moves the rendered HTML into the run folder",
        any(grepl("file.copy(html", src, fixed = TRUE)))
    chk("adoption runner: ships DRY_RUN TRUE",
        any(grepl("^DRY_RUN <- TRUE", src)))
    chk("adoption runner: stage A1 carries NO config delta (production as shipped)",
        any(grepl('tag = "AD-A1-adopted", delta = list()', src, fixed = TRUE)))
  } else chk("adoption runner present", FALSE)
})

# ---------------------------------------------------------------------------
# 43. No trailing whitespace in the R sources. Harmless to R, but `git apply` warns on
#     every patch that introduces it ("warning: N lines add whitespace errors"), which
#     makes a clean application look like a problem to whoever is applying it.
# ---------------------------------------------------------------------------
local({
  for (f in c(list.files("03_R_functions", pattern = "\\.R$", full.names = TRUE),
              list.files("06_diagnostics", pattern = "\\.R$", full.names = TRUE),
              list.files("04_input_files", pattern = "\\.R$", full.names = TRUE),   # 2026-09-10: the builders
              "run_config.R")) {
    ln <- readLines(f, warn = FALSE)
    bad_i <- which(grepl("[ \t]+$", ln))
    chk(sprintf("no trailing whitespace: %s", basename(f)), !length(bad_i),
        if (length(bad_i)) sprintf("line(s) %s", paste(utils::head(bad_i, 5), collapse = ", ")) else "")
  }
})

# ---------------------------------------------------------------------------
# 44. Two defects the 2026-09-07 adoption run exposed, both in verdict code, both
#     surfacing only AFTER the fitting was done.
# ---------------------------------------------------------------------------
local({
  # 44a. The two drivers label the port-total row differently: the pooled writes
  #      "Expected_Catch" and the gear-resolved writes "Catch". Reading one label
  #      against the other file yields numeric(0), and `is.finite(numeric(0)) && ...`
  #      is an ERROR in R, not FALSE. Assert the committed outputs still differ (so the
  #      pattern match is still needed) and that the runner matches by pattern.
  pf <- "05_output/20260904/pooled-CPUE-AD-A1-adopted/port_total_Dungeness_Kept.csv"
  gf <- "05_output/20260905/gear-type-CPUE-model-AD-A2-gear-crosscheck/port_total_Dungeness_Kept.csv"
  if (file.exists(pf) && file.exists(gf)) {
    pl <- utils::read.csv(pf, stringsAsFactors = FALSE)$Estimate
    gl <- utils::read.csv(gf, stringsAsFactors = FALSE)$Estimate
    chk("the two tracks really do label the port row differently (so a pattern match is required)",
        !identical(sort(pl), sort(gl)), sprintf("pooled: %s | gear: %s",
                                                paste(pl, collapse = "/"), paste(gl, collapse = "/")))
  }
  rf <- "06_diagnostics/run_adoption_2026-09-07.R"
  if (file.exists(rf)) {
    src <- readLines(rf, warn = FALSE); code <- src[!grepl("^\\s*#", src)]
    chk("adoption runner: the port row is matched by PATTERN, not by one track's label",
        any(grepl('grepl("^(Expected_)?Catch$", x$Estimate)', code, fixed = TRUE)) &&
        !any(grepl('$Estimate == "Expected_Catch"', code, fixed = TRUE)))
    # 44b. A verdict block must never be able to cost a completed run its output.
    chk("adoption runner: every verdict block is wrapped so a defect cannot lose the run",
        any(grepl(".safe <- function(sid, expr) tryCatch(", code, fixed = TRUE)) &&
        all(grepl("\\.safe\\(", grep("verdict_A[12]\\(dirs", code, value = TRUE))))
  }
  # the same trap in R itself, asserted so the reasoning stays on the record
  chk("R: `is.finite(numeric(0)) && TRUE` errors rather than returning FALSE",
      inherits(try(if (is.finite(numeric(0)) && TRUE) 1 else 2, silent = TRUE), "try-error"))
  chk("R: isTRUE(is.finite(numeric(0))) is the safe form",
      identical(isTRUE(is.finite(numeric(0))), FALSE))
})

# ---------------------------------------------------------------------------
# 45. Season portability (2026-09-09). The program goal is that the model runs on ANY
#     user-selected window: full season, part-season, multi-season span. These pin the
#     four fixes that goal required, functionally where possible.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/build_subseasons.R")
  base <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15",
               pot_closure_start = "2024-09-16", pot_closure_end = "2024-11-30",
               pot_open_date = "2024-12-01")
  full <- build_subseasons(base)
  chk("subseasons: the 2024-25 full-season shape is unchanged (closure + all_gear)",
      length(full) == 2 && identical(sapply(full, `[[`, "name"), c("ring_net_only", "all_gear")))
  summer <- build_subseasons(modifyList(base, list(est_date_start = "2025-06-01", est_date_end = "2025-08-31")))
  chk("subseasons: a window that MISSES the closure runs as one all-gear sub-season",
      length(summer) == 1 && identical(summer[[1]]$gear_regime, "all_gear") &&
      identical(summer[[1]]$start, as.Date("2025-06-01")) && identical(summer[[1]]$end, as.Date("2025-08-31")))
  winter <- build_subseasons(modifyList(base, list(est_date_start = "2024-12-15", est_date_end = "2025-03-15")))
  chk("subseasons: a window entirely AFTER the closure runs as one all-gear sub-season",
      length(winter) == 1 && identical(winter[[1]]$gear_regime, "all_gear"))
  inside <- build_subseasons(modifyList(base, list(est_date_start = "2024-10-01", est_date_end = "2024-10-31")))
  chk("subseasons: a window INSIDE the closure runs as one pot-closure sub-season",
      length(inside) == 1 && identical(inside[[1]]$gear_regime, "pot_closure"))
  chk("subseasons: an inverted closure CONFIG still stops",
      inherits(try(build_subseasons(modifyList(base, list(pot_closure_start = "2024-11-30",
                                                          pot_closure_end = "2024-09-16"))), silent = TRUE), "try-error"))
  # the readers accept a vector season_filter; assert the source so a revert is caught
  for (fset in list(c("03_R_functions/fetch_crab_data.R", "%in% params$season_filter", 2),
                    c("03_R_functions/bss_day_length.R",  "%in% params$season_filter", 1),
                    c("03_R_functions/read_crabbing_holidays.R", "%in% season", 1))) {
    src <- paste(readLines(fset[1], warn = FALSE), collapse = "\n")
    chk(sprintf("season_filter is vector-safe in %s", basename(fset[1])),
        lengths(regmatches(src, gregexpr(fset[2], src, fixed = TRUE))) >= as.integer(fset[3]))
  }
  # the validator: stops on a captured-nothing window, passes on a good one, quiet-able
  source("03_R_functions/validate_season_window.R")
  eff <- data.frame(date = as.Date("2024-10-01") + 0:9, season = "2024-25")
  int <- data.frame(event_date = as.Date("2024-10-02") + 0:9, season = "2024-25")
  okp <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15", season_filter = "2024-25")
  chk("validator: a good season/window pairing passes",
      isTRUE(suppressWarnings(validate_season_window(eff, int, okp, quiet = TRUE))))
  badp <- modifyList(okp, list(est_date_start = "2025-10-01", est_date_end = "2026-09-15"))
  chk("validator: a stale season_filter STOPS with the guide named in the message",
      { e <- try(suppressWarnings(validate_season_window(eff, int, badp, quiet = TRUE)), silent = TRUE)
        inherits(e, "try-error") && grepl("NEW_SEASON_GUIDE", attr(e, "condition")$message) })
  chk("validator: runs inside fetch_crab_data so every driver and batch gets it",
      any(grepl("validate_season_window(gh_effort, gh_interview, params)",   # B46: the Grays Harbor frame
                readLines("03_R_functions/fetch_crab_data.R", warn = FALSE), fixed = TRUE)))
  # the guide exists and is wired in
  chk("NEW_SEASON_GUIDE.md exists", file.exists("07_documentation/NEW_SEASON_GUIDE.md"))
  chk("the guide is linked from README, run_config and the docs index",
      grepl("NEW_SEASON_GUIDE", paste(readLines("README.md", warn = FALSE), collapse = "")) &&
      grepl("NEW_SEASON_GUIDE", paste(readLines("run_config.R", warn = FALSE), collapse = "")) &&
      grepl("NEW_SEASON_GUIDE", paste(readLines("07_documentation/README.md", warn = FALSE), collapse = "")))
  chk("run_config tags its 2024-25-derived settings",
      sum(grepl("SEASON-DERIVED", readLines("run_config.R", warn = FALSE))) >= 5)
})

# ---------------------------------------------------------------------------
# 46. Multi-season spans (2026-09-10). Three requirements for a 2023-24 + 2024-25 run:
#     one pot-closure fit PER season (params$pot_closures), a census window PER season
#     (params$census_windows), and season-level totals in the report (season_totals.csv
#     + table). These pin the builder shapes functionally and the rest at source level,
#     plus internal consistency of the staged two-season config itself.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/build_subseasons.R")
  two <- list(
    est_date_start = "2023-09-16", est_date_end = "2025-09-15",
    season_filter  = c("2023-24", "2024-25"),
    pot_closures = list(
      list(season = "2023-24", start = "2023-09-16", end = "2023-11-30"),
      list(season = "2024-25", start = "2024-09-16", end = "2024-11-30")))
  ss <- build_subseasons(two)
  chk("multi-closure: the 2023-25 span builds exactly four sub-seasons",
      length(ss) == 4)
  chk("multi-closure: names are season-suffixed and alternate closure/all-gear",
      identical(vapply(ss, `[[`, "", "name"),
                c("ring_net_only_2023-24", "all_gear_2023-24",
                  "ring_net_only_2024-25", "all_gear_2024-25")))
  chk("multi-closure: gear regimes alternate pot_closure / all_gear",
      identical(vapply(ss, `[[`, "", "gear_regime"),
                c("pot_closure", "all_gear", "pot_closure", "all_gear")))
  chk("multi-closure: every sub-season carries its season tag",
      identical(vapply(ss, `[[`, "", "season"),
                c("2023-24", "2023-24", "2024-25", "2024-25")))
  chk("multi-closure: dates tile the window with no gap or overlap",
      identical(ss[[1]]$start, as.Date("2023-09-16")) &&
      identical(ss[[1]]$end,   as.Date("2023-11-30")) &&
      identical(ss[[2]]$start, as.Date("2023-12-01")) &&
      identical(ss[[2]]$end,   as.Date("2024-09-15")) &&
      identical(ss[[3]]$start, as.Date("2024-09-16")) &&
      identical(ss[[3]]$end,   as.Date("2024-11-30")) &&
      identical(ss[[4]]$start, as.Date("2024-12-01")) &&
      identical(ss[[4]]$end,   as.Date("2025-09-15")))
  chk("multi-closure: the closure sub-seasons still exclude pots at biweekly resolution",
      identical(ss[[1]]$gear_exclude, c("Pot")) && identical(ss[[1]]$period_bss, "biweekly") &&
      identical(ss[[2]]$gear_exclude, character(0)) && identical(ss[[2]]$period_bss, "month"))
  pre <- build_subseasons(modifyList(two, list(est_date_start = "2023-09-01")))
  chk("multi-closure: a window starting before the first closure gets all_gear_pre_<season1>",
      length(pre) == 5 && identical(pre[[1]]$name, "all_gear_pre_2023-24") &&
      identical(pre[[1]]$end, as.Date("2023-09-15")))
  # single closure through the LIST form falls through to the scalar path: legacy names,
  # so every historical 2024-25 fit label and output filename is preserved.
  one <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15",
              season_filter = "2024-25",
              pot_closures = list(list(season = "2024-25",
                                       start = "2024-09-16", end = "2024-11-30")))
  s1 <- build_subseasons(one)
  chk("multi-closure: a single-closure LIST still yields the legacy names",
      length(s1) == 2 && identical(vapply(s1, `[[`, "", "name"),
                                   c("ring_net_only", "all_gear")))
  # a listed closure entirely OUTSIDE the window is dropped, not an error
  out <- modifyList(one, list(est_date_start = "2025-01-01", est_date_end = "2025-06-30"))
  out$pot_closures <- list(list(season = "2024-25", start = "2024-09-16", end = "2024-11-30"))
  s0 <- build_subseasons(out)
  chk("multi-closure: an out-of-window closure drops to one all-gear sub-season",
      length(s0) == 1 && identical(s0[[1]]$gear_regime, "all_gear"))
  bad1 <- two; bad1$pot_closures[[2]]$start <- "2023-11-01"
  chk("multi-closure: OVERLAPPING closures stop loudly",
      inherits(try(build_subseasons(bad1), silent = TRUE), "try-error"))
  bad2 <- two; bad2$pot_closures[[1]]$season <- NULL
  chk("multi-closure: a closure without a season label stops loudly",
      inherits(try(build_subseasons(bad2), silent = TRUE), "try-error"))

  # census_windows: the guard fires before any data is touched, so it is testable bare
  source("03_R_functions/estimate_comm_charter.R")
  chk("census_windows: an UNNAMED list stops before touching data",
      inherits(try(estimate_comm_charter(NULL, list(
        census_windows = list(c("2024-12-01", "2025-02-08")))), silent = TRUE), "try-error"))
  cc <- paste(readLines("03_R_functions/estimate_comm_charter.R", warn = FALSE), collapse = "\n")
  chk("census_windows: recursion substitutes the scalar keys and keeps by_season",
      grepl("ps\\$census_windows\\s+<- NULL", cc) && grepl("tot$by_season <- per", cc, fixed = TRUE))

  # the validator warns on the exact failure two-season staging found: a season with
  # interviews in the window but ZERO effort counts (its effort would be pure imputation)
  source("03_R_functions/validate_season_window.R")
  eff <- data.frame(date = as.Date("2024-10-01") + 0:9, season = "2024-25")
  int <- data.frame(event_date = c(as.Date("2023-10-01") + 0:9, as.Date("2024-10-01") + 0:9),
                    season = rep(c("2023-24", "2024-25"), each = 10))
  p2 <- list(est_date_start = "2023-09-16", est_date_end = "2025-09-15",
             season_filter = c("2023-24", "2024-25"))
  w <- tryCatch({ validate_season_window(eff, int, p2, quiet = TRUE); NULL },
                warning = function(w) conditionMessage(w))
  chk("validator: zero-effort season raises the pure-imputation warning EVEN under quiet=TRUE",
      !is.null(w) && grepl("ZERO effort counts", w))

  # the report writes season totals: the chunk exists, writes the CSV unconditionally,
  # and the table renders only on multi-season runs
  rmd <- paste(readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE), collapse = "\n")
  chk("report: the season-totals chunk writes season_totals.csv",
      grepl("season-totals", rmd, fixed = TRUE) &&
      grepl('"season_totals.csv"', rmd, fixed = TRUE))

  # the staged two-season config is internally consistent (sourced fresh, not grepped)
  e <- new.env(); sys.source("run_config.R", envir = e)
  rc <- get("run_config", envir = e)
  if (!is.null(rc$pot_closures)) {
    cls <- vapply(rc$pot_closures, function(x) as.character(x$season), "")
    chk("staged config: est window spans every configured closure",
        all(vapply(rc$pot_closures, function(x)
          as.Date(x$start) >= as.Date(rc$est_date_start) &&
          as.Date(x$end)   <= as.Date(rc$est_date_end), logical(1))))
    chk("staged config: season_filter matches the closures' seasons",
        setequal(rc$season_filter, cls))
    chk("staged config: census_windows keys match season_filter",
        is.null(rc$census_windows) || setequal(names(rc$census_windows), rc$season_filter))
  } else {
    chk("staged config: single-season config needs no multi-season consistency", TRUE)
  }
})

# ---------------------------------------------------------------------------
# 47. Boat turnover prior from the OSP/trailer overlap (review items 3 and 5,
#     2026-09-08). tau_boat_prior_mu = "calibration" must resolve to the implied
#     turnover of the chosen metric row, fall back when there is no overlap, refuse a
#     bad string, and every consumer must see ONE number: bss_effort_spec() refuses an
#     unresolved string so a driver that forgets the resolver cannot hand Stan a string.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/bss_turnover_prior.R")
  cal <- tibble(trailer_metric = c("trailer_mean_per_visit", "trailer_max_per_day", "trailer_sum_per_day"),
                n = 61L, corr = c(0.965, 0.983, 0.899), ols_slope = 0, ols_intercept = 0,
                origin_slope = c(0.33, 0.365, 0.498), implied_turnover = c(3.03, 2.74, 2.01))
  ov <- list(calibration = cal)
  Pc <- modifyList(P, list(tau_boat_prior_mu = "calibration", tau_boat_prior_sigma = 0.5,
                           tau_boat_prior_mu_fallback = 2.7, tau_boat_calibration_min_pairs = 3))
  r1 <- bss_resolve_tau_boat_prior(Pc, ov, quiet = TRUE)
  chk("tau prior: calibration resolves to the mean-per-visit implied turnover",
      isTRUE(all.equal(r1$tau_boat_prior_mu, 3.03)) && grepl("calibration", r1$tau_boat_prior_source))
  r2 <- bss_resolve_tau_boat_prior(modifyList(Pc, list(tau_boat_calibration_metric = "trailer_max_per_day")), ov, quiet = TRUE)
  chk("tau prior: the metric key selects its row", isTRUE(all.equal(r2$tau_boat_prior_mu, 2.74)))
  r3 <- bss_resolve_tau_boat_prior(Pc, NULL, quiet = TRUE)
  chk("tau prior: no overlap -> the SEASON-DERIVED fallback",
      isTRUE(all.equal(r3$tau_boat_prior_mu, 2.7)) && grepl("fallback", r3$tau_boat_prior_source))
  cal_thin <- cal; cal_thin$n <- 2L
  r4 <- bss_resolve_tau_boat_prior(Pc, list(calibration = cal_thin), quiet = TRUE)
  chk("tau prior: too few paired days -> fallback", isTRUE(all.equal(r4$tau_boat_prior_mu, 2.7)))
  r5 <- bss_resolve_tau_boat_prior(modifyList(Pc, list(tau_boat_prior_mu = 1.2)), ov, quiet = TRUE)
  chk("tau prior: a number passes through untouched (historical behaviour)",
      isTRUE(all.equal(r5$tau_boat_prior_mu, 1.2)) && grepl("numeric", r5$tau_boat_prior_source))
  chk("tau prior: a bad string errors rather than silently falling back",
      inherits(tryCatch(bss_resolve_tau_boat_prior(modifyList(Pc, list(tau_boat_prior_mu = "guess")), ov, quiet = TRUE),
                        error = function(e) e), "error"))
  chk("tau prior: an unknown metric errors",
      inherits(tryCatch(bss_resolve_tau_boat_prior(modifyList(Pc, list(tau_boat_calibration_metric = "nope")), ov, quiet = TRUE),
                        error = function(e) e), "error"))
  # consumers must see a number
  chk("effort spec refuses an UNRESOLVED tau_boat_prior_mu",
      inherits(tryCatch(bss_effort_spec(FALSE, days289, Pc), error = function(e) e), "error"))
  sp <- bss_effort_spec(FALSE, days289, r1)
  chk("effort spec: resolved prior centre reaches L_data",
      isTRUE(all.equal(unique(sp$L_data), 3.03)) && isTRUE(all.equal(unique(sp$L_prior_sigma), 0.5)))
  # both drivers resolve it, and do so after the overlap diagnostic and before the PE
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    i_res <- grep("bss_resolve_tau_boat_prior(params, osp_overlap)", d, fixed = TRUE)
    i_ov  <- grep("diagnose_osp_trailer_overlap(osp_boat", d, fixed = TRUE)
    i_pe  <- grep("run_pe_pooled\\(summ_ss|run_pe_gear\\(summ_ss", d)
    chk(sprintf("%s: resolves the boat turnover prior after the overlap and before the PE", basename(drv)),
        length(i_res) == 1 && length(i_ov) >= 1 && length(i_pe) >= 1 && i_res > min(i_ov) && i_res < min(i_pe))
  }
  # the PE reads the same key the BSS prior is built from (item 5)
  for (f in c("03_R_functions/run_pe_pooled.R", "03_R_functions/run_pe_gear.R")) {
    src <- readLines(f, warn = FALSE); src <- src[!grepl("^\\s*#", src)]
    chk(sprintf("%s: expands the boat on params$tau_boat_prior_mu", basename(f)),
        any(grepl("params$tau_boat_prior_mu", src, fixed = TRUE)))
  }
  # shipped config
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("shipped: tau_boat_prior_mu = \"calibration\"", identical(rc$tau_boat_prior_mu, "calibration"))
  chk("shipped: tau_boat_prior_sigma = 0.5 (wide, because the overlap is in the likelihood too)",
      identical(rc$tau_boat_prior_sigma, 0.5))
  chk("shipped: shared_tau_sigma pinned at 0.15 (no side effect from the wider prior)",
      identical(rc$shared_tau_sigma, 0.15))
  chk("shipped: tau sensitivity grid brackets the calibration",
      min(rc$tau_sensitivity_grid) < 2.7 && max(rc$tau_sensitivity_grid) > 3.0)
})

# ---------------------------------------------------------------------------
# 48. Unit-aware fishing-time filters (review item 8, 2026-09-08). The 0.5 h threshold
#     belongs to the time-denominated units; under gear-deployments the only guard is
#     "no positive time AND zero catch" (gear set, not yet fished). Tested on a synthetic
#     frame through the pure helper apply_fishing_time_filters().
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/fetch_crab_data.R")
  fr <- tibble(
    population         = c("shore","shore","shore","private_boat","private_boat","private_boat","comm_charter"),
    hours_fished       = c(0.1,    0,      3,      0,             0.1,           2,             NA),
    crabber_hours_calc = c(0.2,    0,      6,      0,             0.1,           4,             NA),
    gear_hours         = c(0.2,    0,      6,      0,             0.3,           6,             NA),
    fishing_time_total = c(0.2,    NA,     6,      NA,            0.1,           4,             NA),
    dungeness_kept     = c(0,      0,      5,      0,             0,             9,             40),
    red_rock_kept      = 0)
  Pd <- list(shore_effort_unit = "gear-deployments", min_fishing_time = 0.5, drop_unfished_zero_catch = TRUE)
  out <- apply_fishing_time_filters(fr, Pd, quiet = TRUE)
  chk("deployments: a 0.1 h shore trip is KEPT (no time threshold)", 1 %in% which(fr$population == "shore") && nrow(out[out$population=="shore",]) == 2)
  chk("deployments: a no-time zero-catch shore row is dropped (unfished)", !any(out$population == "shore" & is.na(out$fishing_time_total)))
  chk("deployments: a no-time zero-catch BOAT row is dropped (pots still soaking)", sum(out$population == "private_boat") == 2)
  chk("deployments: a 0.1 h zero-catch boat trip is KEPT", any(out$population == "private_boat" & out$hours_fished == 0.1))
  chk("deployments: a commercial row with no hours but catch is KEPT", any(out$population == "comm_charter"))
  Pt <- modifyList(Pd, list(shore_effort_unit = "crabber-hours"))
  out_t <- apply_fishing_time_filters(fr, Pt, quiet = TRUE)
  chk("time unit: the 0.5 h threshold still applies to shore", sum(out_t$population == "shore") == 1)
  chk("time unit: the boat is untouched by the shore unit", sum(out_t$population == "private_boat") == 2)
  Pg <- modifyList(Pd, list(drop_unfished_zero_catch = FALSE))
  out_g <- apply_fishing_time_filters(fr, Pg, quiet = TRUE)
  chk("guard off: nothing is dropped under deployments", nrow(out_g) == nrow(fr))
  chk("helper strips its scratch columns", !any(grepl("^\\.", names(out))))
  src <- readLines("03_R_functions/fetch_crab_data.R", warn = FALSE); src <- src[!grepl("^\\s*#", src)]
  chk("reader calls the helper (no inline time threshold survives)",
      any(grepl("apply_fishing_time_filters(gh_interview", src, fixed = TRUE)) &&
      !any(grepl("fishing_time_total >= params$min_fishing_time", src, fixed = TRUE)))
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("shipped: drop_unfished_zero_catch = TRUE", identical(rc$drop_unfished_zero_catch, TRUE))
})

# ---------------------------------------------------------------------------
# 49. Gear split by trip-level bootstrap (review item 6, 2026-09-08). The Dirichlet on
#     crab COUNTS treated every crab as an independent draw; crab arrive in trips, so
#     the intervals were several times too narrow. The bootstrap resamples trips.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/gear_share_bootstrap.R")
  chk("gear classes: the report's regex is preserved",
      identical(gear_primary_class(c("Pot", "Ring net", "Collapsible trap or ring", "Fishing rod with snare", "Slip ring pot", "Other thing")),
                c("Pot", "Ring Net", "Trap", "Snare", "Other", "Other")))
  set.seed(7)
  iv <- tibble(gear_primary = c(rep("Pot", 30), rep("Ring Net", 20), rep("Trap", 10)),
               catch = c(rpois(30, 8), rpois(20, 3), rpois(10, 4)))
  sh <- gear_share_bootstrap(iv, n_boot = 500, seed = 1)
  chk("bootstrap shares: rows sum to one", all(abs(rowSums(sh) - 1) < 1e-9))
  chk("bootstrap shares: median tracks the point share",
      abs(median(sh[, "Pot"]) - gear_share_point(iv)[["Pot"]]) < 0.03)
  ct <- tapply(iv$catch, factor(iv$gear_primary, levels = gear_primary_levels), sum); ct[is.na(ct)] <- 0
  g  <- matrix(rgamma(500 * 5, shape = rep(ct + 0.5, each = 500)), 500); g <- g / rowSums(g)
  chk("bootstrap interval is WIDER than the crab-count Dirichlet (clustering respected)",
      diff(quantile(sh[, "Pot"], c(.025, .975))) > 1.5 * diff(quantile(g[, 1], c(.025, .975))))
  chk("bootstrap shares: seed reproduces", identical(gear_share_bootstrap(iv, 50, seed = 3), gear_share_bootstrap(iv, 50, seed = 3)))
  chk("no catch -> all-NA shares, no error",
      all(is.na(gear_share_bootstrap(tibble(gear_primary = "Pot", catch = 0), 10))))
  chk("empty frame -> all-NA shares, no error",
      all(is.na(gear_share_bootstrap(tibble(gear_primary = character(), catch = numeric()), 10))))
  d <- readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE); d <- d[!grepl("^\\s*#", d)]
  chk("pooled driver: gear split uses the bootstrap helper, per sub-season",
      any(grepl("gear_share_bootstrap(iv", d, fixed = TRUE)) &&
      any(grepl("gear_component_total_draws(pop, ss$name", d, fixed = TRUE)) &&
      !any(grepl("rdir(n_gd, counts + 0.5)", d, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 50. Shore turnover from the I/E time column (review item 2, 2026-09-08). The peak-based
#     1.7 assumes the gear count is taken at the daily peak; the estimator evaluates
#     presence at the season's own count hours. Synthetic days with KNOWN curves.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/bss_day_length.R")
  chk("ie hour parser: text HH:MM", isTRUE(all.equal(.ie_hour_of(c("07:45:00", "13:30")), c(7.75, 13.5))))
  chk("ie hour parser: POSIXct", isTRUE(all.equal(.ie_hour_of(as.POSIXct("1899-12-31 10:15:00", tz = "UTC")), 10.25)))
  chk("ie hour parser: garbage -> NA", is.na(.ie_hour_of("noon")))
  # two synthetic days: presence 4 at 08-09, 8 at 10-14, 4 at 15-16; arrivals 20 and 40
  mk <- function(d, scale, dt) tibble(event_date = as.Date(d), season = "t", day_type = dt,
                                     hour = seq(8, 16.75, by = 0.25),
                                     crabbers_on = c(rep(scale * 20 / 36, 36)),
                                     crabber_flow = scale * c(rep(4, 8), rep(8, 20), rep(4, 8)))
  ivs <- bind_rows(mk("2025-01-01", 1, "Weekday"), mk("2025-01-04", 2, "Weekend"))
  se  <- tibble(event_date = as.Date("2025-01-01"), count_sequence = 1:3, count_hour = c(10.5, 12.25, 14.0))
  st <- estimate_shore_turnover(ivs, se, list(bss_max_count_seq = 3), n_boot = 50, quiet = TRUE)
  chk("turnover: counts inside the 10-14 window -> arrivals / 8 per unit (= 20/8 = 2.5)",
      isTRUE(all.equal(st$tau, 2.5, tolerance = 1e-6)), sprintf("%.4f", st$tau))
  chk("turnover: peak-based value reported beside it (= 2.5 here too, curve is flat at peak)",
      isTRUE(all.equal(st$tau_peak, 2.5, tolerance = 1e-6)))
  se2 <- tibble(event_date = as.Date("2025-01-01"), count_sequence = 1L, count_hour = 8.5)
  st2 <- estimate_shore_turnover(ivs, se2, list(bss_max_count_seq = 3), n_boot = 50, quiet = TRUE)
  chk("turnover: an 08:30 count sits at half the peak -> tau doubles to 5",
      isTRUE(all.equal(st2$tau, 5, tolerance = 1e-6)), sprintf("%.4f", st2$tau))
  chk("turnover: profile carries one row per covered hour with n_days", all(st$profile$n_days == 2) && nrow(st$profile) == 9)
  chk("turnover: by-day-type table has both types", setequal(st$by_day_type$dt, c("weekday", "weekend")))
  chk("turnover: NULL intervals -> unavailable, no error", identical(estimate_shore_turnover(NULL, se, quiet = TRUE)$method, "unavailable"))
  se3 <- tibble(event_date = as.Date("2025-01-01"), count_sequence = 1L, count_hour = 22)
  chk("turnover: count hours outside every survey -> unavailable, no error",
      identical(estimate_shore_turnover(ivs, se3, quiet = TRUE)$method, "unavailable"))
  # resolver
  Pn <- list(tau_shore_prior_mu = 1.7, tau_shore_prior_sigma = 0.3, shared_tau_min_obs = 15)
  chk("resolver: numeric config passes through", identical(bss_resolve_tau_shore_prior(Pn, st, quiet = TRUE)$tau_shore_prior_mu, 1.7))
  Pd <- list(tau_shore_prior_mu = "derived", tau_shore_prior_sigma = "derived", shared_tau_min_obs = 15,
             tau_shore_prior_sigma_floor = 0.10, tau_shore_derive_min_days = 2)
  rd <- bss_resolve_tau_shore_prior(Pd, st, quiet = TRUE)
  chk("resolver: derived centre = the estimate", isTRUE(all.equal(rd$tau_shore_prior_mu, 2.5, tolerance = 1e-6)))
  chk("resolver: derived log-SD is floored at 0.10", isTRUE(all.equal(rd$tau_shore_prior_sigma, 0.10)))
  chk("resolver: derived level waives the SHORE floor only",
      is.list(rd$shared_tau_min_obs) && rd$shared_tau_min_obs$shore == 0 && rd$shared_tau_min_obs$private_boat == 15)
  rf <- bss_resolve_tau_shore_prior(modifyList(Pd, list(tau_shore_derive_min_days = 100, tau_shore_prior_mu_fallback = 1.7)), st, quiet = TRUE)
  chk("resolver: too few days -> fallback centre, prior SD 0.3", identical(rf$tau_shore_prior_mu, 1.7) && identical(rf$tau_shore_prior_sigma, 0.3))
  chk("resolver: bad string errors",
      inherits(tryCatch(bss_resolve_tau_shore_prior(list(tau_shore_prior_mu = "guess"), st, quiet = TRUE), error = function(e) e), "error"))
  chk("effort spec refuses an UNRESOLVED shore prior",
      inherits(tryCatch(bss_effort_spec(TRUE, days289, modifyList(P, list(tau_shore_prior_mu = "derived"))), error = function(e) e), "error"))
  # per-population shared-tau keys reach the guard
  sp <- bss_effort_spec(TRUE, days289, modifyList(P, list(tau_shore_prior_mu = 2.5, tau_shore_prior_sigma = 0.1)))
  Pl <- modifyList(P, list(shared_tau = TRUE, shared_tau_min_obs = list(shore = 0, private_boat = 15),
                           shared_tau_sigma = list(shore = 0.35, private_boat = 0.15)))
  on_s <- bss_shared_tau_data(sp, sp$L_data, sp$L_prior_sigma, Pl, population_name = "shore", n_informed = 4L, quiet = TRUE)
  chk("shared tau: per-population floor lets shore share a level on 4 days when its floor is 0",
      on_s$shared_tau == 1L && isTRUE(all.equal(on_s$shared_tau_sigma, 0.35)))
  spb <- bss_effort_spec(FALSE, days289, P)
  off_b <- bss_shared_tau_data(spb, spb$L_data, spb$L_prior_sigma, Pl, population_name = "private_boat", n_informed = 4L, quiet = TRUE)
  chk("shared tau: the boat keeps its own floor and spread", off_b$shared_tau == 0L && isTRUE(all.equal(off_b$shared_tau_sigma, 0.15)))
  # drivers derive, write, resolve, in that order, before the PE
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    i_est <- grep("estimate_shore_turnover(attr(ie_data", d, fixed = TRUE)
    i_res <- grep("bss_resolve_tau_shore_prior(params, shore_turnover)", d, fixed = TRUE)
    i_pe  <- grep("run_pe_pooled\\(summ_ss|run_pe_gear\\(summ_ss", d)
    chk(sprintf("%s: derives and resolves the shore turnover before the PE", basename(drv)),
        length(i_est) == 1 && length(i_res) == 1 && i_est < i_res && i_res < min(i_pe))
  }
  src <- readLines("03_R_functions/fetch_crab_data.R", warn = FALSE); src <- src[!grepl("^\\s*#", src)]
  chk("reader keeps count_hour on the shore and boat effort tables", sum(grepl("count_hour", src)) >= 2)
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("shipped: tau_shore_prior_mu = derived (adopted 2026-09-09; the pre-time-column 1.7 was the wrong quantity)",
      identical(rc$tau_shore_prior_mu, "derived") && identical(rc$tau_shore_prior_sigma, "derived"))
  chk("shipped: the derived-prior keys exist", !is.null(rc$tau_shore_prior_sigma_floor) && !is.null(rc$tau_shore_derive_min_days))
})

# ---------------------------------------------------------------------------
# 51. Boat contacts as the crabbing-fraction classification (review item 1A,
#     2026-09-08). The reader dropped the 216 non-crabbing boat contacts at
#     crabbers > 0; they are the f data. Pure helpers on synthetic frames.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/fetch_crab_data.R")
  fr <- tibble(population = c(rep("private_boat", 6), "shore", "private_boat"),
               event_date = as.Date(c(rep("2025-08-01", 4), rep("2025-08-02", 2), "2025-08-01", "2025-08-03")),
               crabbers   = c(0, 2, 0, 1,   0, 0,   3, NA))
  bc <- boat_contacts_from_interviews(fr)
  chk("contacts: one row per day, boats_total counts every classified private boat",
      nrow(bc) == 2 && bc$boats_total[bc$event_date == as.Date("2025-08-01")] == 4)
  chk("contacts: boats_crabbing counts crabbers > 0", bc$boats_crabbing[bc$event_date == as.Date("2025-08-01")] == 2 &&
        bc$boats_crabbing[bc$event_date == as.Date("2025-08-02")] == 0)
  chk("contacts: NA crabbers is unclassified (day dropped), shore rows ignored", !as.Date("2025-08-03") %in% bc$event_date)
  chk("contacts: empty / malformed input -> empty table, no error",
      nrow(boat_contacts_from_interviews(NULL)) == 0 && nrow(boat_contacts_from_interviews(tibble(x = 1))) == 0)
  # source assembly
  dwg <- list(boat_contacts = bc)
  ie  <- tibble(x = 1); attr(ie, "crab_fraction_rows") <- tibble(event_date = as.Date("2025-08-05"), boats_crabbing = 3, boats_total = 9)
  r_both <- crab_fraction_source_rows(dwg, ie, list(crab_fraction_source = "both"), quiet = TRUE)
  chk("source both: interview and egress rows are bound, tagged", nrow(r_both) == 3 && setequal(unique(r_both$source), c("interviews", "ie")))
  chk("source interviews: egress rows excluded", all(crab_fraction_source_rows(dwg, ie, list(crab_fraction_source = "interviews"), quiet = TRUE)$source == "interviews"))
  chk("source ie: contacts excluded", all(crab_fraction_source_rows(dwg, ie, list(crab_fraction_source = "ie"), quiet = TRUE)$source == "ie"))
  chk("source: bad value errors", inherits(tryCatch(crab_fraction_source_rows(dwg, ie, list(crab_fraction_source = "osp"), quiet = TRUE), error = function(e) e), "error"))
  chk("source: rows carry the columns crab_fraction.R aggregates", all(c("event_date", "boats_crabbing", "boats_total") %in% names(r_both)))
  # the rows reach the existing Binomial machinery
  Pm <- modifyList(P, list(crab_fraction_strata = "month", crab_fraction_min_obs = 3))
  Pm$crab_fraction_rows <- r_both
  cf <- crab_fraction_stan_data(FALSE, mkdays("2025-08-01", 31), Pm, quiet = TRUE)
  chk("month strata: contact rows aggregate into the stratum Binomial", cf$crab_fraction_n_total[1] == 15 && cf$crab_fraction_n_crab[1] == 5)
  # chronological, year-qualified labels
  lab <- crab_fraction_strata_labels(as.Date(c("2024-12-15", "2025-01-10", "2025-09-01")), list(crab_fraction_strata = "month"))
  chk("month strata are year-qualified and sort chronologically", identical(sort(unique(lab)), c("2024-12", "2025-01", "2025-09")))
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    chk(sprintf("%s: crab_fraction_rows come from crab_fraction_source_rows()", basename(drv)),
        any(grepl("crab_fraction_source_rows(dwg, ie_data, params)", d, fixed = TRUE)))
  }
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("shipped: crab_fraction_source = both", identical(rc$crab_fraction_source, "both"))
})

# ---------------------------------------------------------------------------
# 52. The dynamic crabbing fraction (review item 1B, 2026-09-08). A logit random walk
#     across chronological strata, observed per day by the sampler contacts and by OSP's
#     crabbing-only column through f(1 - c). Tests: the walk order for every
#     stratification, the per-day rows, exact legacy back-compat when the key is
#     absent, the PE's interpolating point f, the decoupled rules, and that both Stan
#     models and both preps carry the new fields.
# ---------------------------------------------------------------------------
local({
  # walk order
  w <- crab_fraction_walk_structure(c("2024-11", "2024-12", "2025-02"), "month")
  chk("walk: month strata chain k -> k-1 with the month gap", identical(w$prev, c(0L, 1L, 2L)) && isTRUE(all.equal(w$gap, c(1, 1, 2))) && isTRUE(all.equal(w$pos, c(0, 1, 3))))
  w2 <- crab_fraction_walk_structure(c("2024-12_wkdy", "2024-12_wknd", "2025-01_wkdy", "2025-01_wknd"), "month_day_type")
  chk("walk: month x day-type walks within the day type", identical(w2$prev, c(0L, 0L, 1L, 2L)) && identical(w2$chain, c(1L, 2L, 1L, 2L)))
  w3 <- crab_fraction_walk_structure(c("open", "shut"), "opener")
  chk("walk: non-monthly strata are all anchored (no walk)", all(w3$prev == 0L))
  chk("walk: predecessor always precedes (Stan rejects otherwise)", all(w2$prev < seq_along(w2$prev)) && all(w$prev < seq_along(w$prev)))

  # per-day rows and the inert fields
  Pd <- modifyList(P, list(crab_fraction_strata = "month", crab_fraction_dynamic = TRUE, use_osp_crab_lower = FALSE))
  rows <- tibble(event_date = days289$event_date[c(1, 2, 40, 41, 100)], boats_total = c(3, 2, 5, 1, 4), boats_crabbing = c(3, 1, 4, 0, 1))
  Pd$crab_fraction_rows <- rows
  cfd <- crab_fraction_stan_data(FALSE, days289, Pd, quiet = TRUE)
  chk("dynamic: flag set, one CFI row per contact day", cfd$crab_fraction_dynamic == 1L && cfd$CFI_n == 5)
  chk("dynamic: CFI rows carry stratum, total, crab", all(cfd$cfi_stratum >= 1 & cfd$cfi_stratum <= cfd$n_f_strata) && sum(cfd$cfi_total) == 15 && sum(cfd$cfi_crab) == 9)
  chk("dynamic: per-stratum sums are carried UNGATED (reporting), min_obs not applied", sum(cfd$crab_fraction_n_total) == 15)
  chk("dynamic: walk fields sized K, level prior at logit(set)", length(cfd$f_walk_prev) == cfd$n_f_strata && isTRUE(all.equal(cfd$f_level_mu, qlogis(0.3))))
  chk("dynamic: priors pass through with their defaults", cfd$f_level_sd == 1.5 && cfd$f_walk_sd_prior == 1.5 && cfd$f_walk_df == 4 && cfd$cfi_kappa_prior_mu == 20 &&
        isTRUE(all.equal(cfd$c_level_mu, qlogis(0.3))) && cfd$c_level_sd == 1.5 && cfd$c_walk_sd_prior == 1.5 && cfd$cfc_kappa_prior_mu == 20)
  chk("dynamic: without typed contacts or OSP the combo walk is off (combo_dynamic = 0, CFC_n = 0)", cfd$combo_dynamic == 0L && cfd$CFC_n == 0)
  chk("dynamic: audit table attached with one row per stratum", is.data.frame(attr(cfd, "f_strata")) && nrow(attr(cfd, "f_strata")) == cfd$n_f_strata)
  chk("dynamic: rows outside the fit window are excluded",
      crab_fraction_stan_data(FALSE, days76, Pd, quiet = TRUE)$CFI_n == 0)
  # exact legacy back-compat: without the key the legacy fields are identical and the walk is inert
  Pl <- Pd; Pl$crab_fraction_dynamic <- NULL
  cfl <- crab_fraction_stan_data(FALSE, days289, Pl, quiet = TRUE)
  leg_fields <- c("crab_fraction_estimate", "n_f_strata", "f_stratum", "crab_fraction_value", "crab_fraction_alpha0", "crab_fraction_beta0", "osp_crab_lower", "OSPF_n")
  chk("legacy (key absent): dynamic flag 0, no CFI rows, legacy fields unchanged",
      cfl$crab_fraction_dynamic == 0L && cfl$CFI_n == 0 && identical(cfl[leg_fields], cfd[leg_fields]))
  chk("legacy (key absent): min_obs gating still applies to the stratum Binomial", sum(cfl$crab_fraction_n_total) == 0)
  # pinned and shore paths carry the fields, inert
  chk("pinned f: dynamic fields present and inert",
      crab_fraction_stan_data(FALSE, days289, modifyList(Pd, list(crab_fraction_fixed = 0.4)), quiet = TRUE)$crab_fraction_dynamic == 0L)
  chk("shore: dynamic fields present and inert", crab_fraction_stan_data(TRUE, days289, Pd, quiet = TRUE)$CFI_n == 0)
  # the OSP crabbing-only stream still feeds per-day rows under the walk
  Pdo <- setp(Pd, osp_crab_rows = osp_rows, use_osp_crab_lower = TRUE)
  cfo <- crab_fraction_stan_data(FALSE, days289, Pdo, quiet = TRUE)
  chk("dynamic + OSP: OSP daily rows still emitted, stream flagged on", cfo$osp_crab_lower == 1L && cfo$OSPF_n == 200 && cfo$crab_fraction_dynamic == 1L)
  # the PE point f: informed strata = shrunken share (kappa_pe 1), thin strata interpolated on the logit scale
  Pp <- modifyList(Pd, list(crab_fraction_min_obs = 20))
  big <- tibble(event_date = c(days289$event_date[1:10], days289$event_date[100:109]), boats_total = 5,
                boats_crabbing = c(rep(5, 10), rep(1, 10)))   # Dec ~1.0 (50 boats), Mar ~0.2 (50 boats), Jan-Feb empty
  Pp$crab_fraction_rows <- big
  pd <- crab_fraction_point_day(TRUE, days289, Pp)
  lab <- format(days289$event_date, "%Y-%m")
  f_dec <- unique(pd[lab == "2024-12"]); f_jan <- unique(pd[lab == "2025-01"]); f_feb <- unique(pd[lab == "2025-02"]); f_mar <- unique(pd[lab == "2025-03"])
  chk("PE point: informed stratum = (n_crab + a)/(n + 1)", isTRUE(all.equal(f_dec, (50 + 0.3) / 51)) && isTRUE(all.equal(f_mar, (10 + 0.3) / 51)))
  chk("PE point: thin strata interpolate on the logit scale between informed neighbours",
      f_dec > f_jan && f_jan > f_feb && f_feb > f_mar && isTRUE(all.equal(qlogis(f_jan) - qlogis(f_feb), qlogis(f_feb) - qlogis(f_mar), tolerance = 1e-6)))
  chk("PE point: strata after the last informed one carry its value", isTRUE(all.equal(unique(pd[lab == "2025-08"]), f_mar)))
  chk("PE point: no data at all -> the set value", isTRUE(all.equal(unique(crab_fraction_point_day(TRUE, days289, modifyList(Pd, list(crab_fraction_rows = NULL)))), 0.3)))
  # decoupled rules
  sd_off <- list(apply_crab_fraction = 1L, crab_fraction_estimate = 1L, crab_fraction_dynamic = 0L, CFI_n = 0L, OSPF_n = 0L, osp_crab_lower = 0L, combo_dynamic = 0L, CFC_n = 0L)
  r_off <- bss_decoupled_reasons(c("sigma_f_out", "cfi_kappa_out", "combo_c_out", "f_crab[1]"), sd_off)
  chk("decoupled: walk off -> its scale parameters are 'not in the model', f_crab is not flagged", all(!is.na(r_off[1:3])) && is.na(r_off[4]))
  sd_on <- modifyList(sd_off, list(crab_fraction_dynamic = 1L, CFI_n = 12L))
  r_on <- bss_decoupled_reasons(c("sigma_f_out", "cfi_kappa_out", "combo_c_out", "f_crab[1]"), sd_on)
  chk("decoupled: walk live with contacts -> sigma_f and kappa_I are estimates; combo_c not in the model without typed contacts or OSP",
      is.na(r_on[1]) && is.na(r_on[2]) && !is.na(r_on[3]) && is.na(r_on[4]))
  sd_c <- modifyList(sd_on, list(combo_dynamic = 1L, CFC_n = 9L))
  r_c <- bss_decoupled_reasons(c("combo_c_out[1]", "sigma_c_out", "cfc_kappa_out", "f_lower[1]"), sd_c)
  chk("decoupled: combo walk live with typed contacts -> c, sigma_c, kappa_C and f_lower are estimates", all(is.na(r_c)))
  sd_c0 <- modifyList(sd_on, list(combo_dynamic = 1L, CFC_n = 0L, osp_crab_lower = 1L, OSPF_n = 0L))
  chk("decoupled: combo walk live with no typed day and no OSP day -> prior only",
      all(!is.na(bss_decoupled_reasons(c("combo_c_out[1]", "sigma_c_out", "f_lower[1]"), sd_c0))))
  sd_bare <- modifyList(sd_on, list(CFI_n = 0L))
  r_bare <- bss_decoupled_reasons(c("sigma_f_out", "cfi_kappa_out", "f_crab[1]"), sd_bare)
  chk("decoupled: walk live with NO rows -> sigma_f, kappa_I and f are prior-only", all(!is.na(r_bare)))
  # both Stan models declare the block; both preps forward it (9b covers the exact set)
  for (m in c("02_stan_models/crab_bss_pooled.stan", "02_stan_models/crab_bss_gear_resolved.stan")) {
    nm <- bss_stan_data_names(m)
    chk(sprintf("%s declares the dynamic-f data block", basename(m)),
        all(c("crab_fraction_dynamic", "f_walk_prev", "f_walk_gap", "f_level_mu", "f_level_sd", "f_walk_sd_prior", "f_walk_df",
              "CFI_n", "cfi_stratum", "cfi_total", "cfi_crab", "cfi_kappa_prior_mu",
              "combo_dynamic", "c_level_mu", "c_level_sd", "c_walk_sd_prior", "CFC_n", "cfc_stratum", "cfc_crab", "cfc_combo", "cfc_kappa_prior_mu") %in% nm))
    src <- paste(readLines(m, warn = FALSE), collapse = "\n")
    chk(sprintf("%s keeps the legacy construction gated on crab_fraction_dynamic = 0", basename(m)),
        grepl("n_f_leg = crab_fraction_estimate * (1 - crab_fraction_dynamic)", src, fixed = TRUE) &&
          grepl("sigma_f_out", src, fixed = TRUE) && grepl("combo_c_out", src, fixed = TRUE) &&
          grepl("int n_c_dyn = n_f_dyn * combo_dynamic", src, fixed = TRUE) && grepl("sigma_c_out", src, fixed = TRUE))
    chk(sprintf("%s: f still enters generated quantities only (no f_crab in an effort or catch likelihood)", basename(m)),
        !grepl("neg_binomial_2\\([^;]*f_crab", src) && !grepl("lognormal\\([^;]*f_crab", src))
  }
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    chk(sprintf("%s: reports sigma_f_out / cfi_kappa_out / combo_c_out and writes crab_fraction_strata_*.csv", basename(drv)),
        any(grepl("sigma_f_out", d, fixed = TRUE)) && any(grepl("combo_c_out", d, fixed = TRUE)) && any(grepl("crab_fraction_strata_", d, fixed = TRUE)))
  }
  chk("structural summary lists the dynamic-f reporters",
      all(c("sigma_f_out", "cfi_kappa_out", "combo_c_out", "sigma_c_out", "cfc_kappa_out") %in% {
        b <- body(bss_structural_summary); s <- paste(deparse(b), collapse = " "); unlist(regmatches(s, gregexpr("[a-z_]+_out", s))) }))
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("shipped: the dynamic-f priors exist with the documented values",
      identical(rc$crab_fraction_level_sd, 1.5) && identical(rc$crab_fraction_walk_sd_prior, 1.5) && identical(rc$crab_fraction_walk_df, 4) &&
        identical(rc$crab_fraction_combo_level, 0.3) && identical(rc$crab_fraction_combo_level_sd, 1.5) &&
        identical(rc$crab_fraction_combo_walk_sd_prior, 1.5) && identical(rc$crab_fraction_pe_prior_kappa, 1))
})

# 52b. The prior-vs-posterior table carries the walk's scale parameters when the walk is
#      live (the one file that says how much the data moved sigma_f), keyed with the index
#      exactly as tau_bar[1] is, because both are length-1 vectors.
local({
  t <- paste(readLines("03_R_functions/save_run_diagnostics.R", warn = FALSE), collapse = "\n")
  chk("prior_vs_posterior: sigma_f[1] / cfi_kappa[1] / sigma_c[1] / cfc_kappa[1] rows are added under a live walk, index-named",
      grepl("prior_tbl$`sigma_f[1]`", t, fixed = TRUE) && grepl("prior_tbl$`cfi_kappa[1]`", t, fixed = TRUE) &&
        grepl("prior_tbl$`sigma_c[1]`", t, fixed = TRUE) && grepl("prior_tbl$`cfc_kappa[1]`", t, fixed = TRUE) &&
        grepl("sd_p$crab_fraction_dynamic", t, fixed = TRUE) && grepl("sd_p$combo_dynamic", t, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 53. The commercial/charter census (review item 4; 2026-09-08 expansion, 2026-09-09 the
#     exact sum). Synthetic tally: a 14-day window, 10 sampled days, so the observed /
#     imputed split and the variances can be checked by hand. The shipped
#     census_expansion = "none" treats the 4 unsampled days as days with no operation;
#     "day_type" is the 2026-09-08 expansion and must reproduce its arithmetic.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/estimate_comm_charter.R")
  cal <- seq(as.Date("2025-01-06"), as.Date("2025-01-19"), by = "day")      # Mon .. Sun, two weeks
  wk  <- weekdays(cal) %in% c("Saturday", "Sunday")
  # sampled: 8 weekdays (of 10) and 2 weekend days (of 4); commercial tallies vary by day
  samp <- c(cal[!wk][1:8], cal[wk][1:2])
  tally <- tibble(date = samp, commercial_tally = c(4, 6, 5, 7, 4, 6, 5, 3, 8, 10), charter_tally = c(1, 0, 1, 1, 0, 1, 1, 0, 2, 2))
  ints  <- tibble(population = "comm_charter", event_date = rep(samp, each = 2),
                  boat_type_clean = rep(c("Commercial", "Charter"), 10), dungeness_kept = rep(c(40, 60), 10), red_rock_kept = 0)
  dwg <- list(comm_tally = tally, interview = ints)
  Pc0 <- list(census_start_date = "2025-01-06", census_end_date = "2025-01-19", days_wkend = c("Saturday", "Sunday"),
              crabbing_holiday_dates = as.Date(character()))
  est_day <- tally$commercial_tally * 40 + tally$charter_tally * 60
  wkd <- !weekdays(tally$date) %in% c("Saturday", "Sunday")
  obs <- sum(est_day); imp <- 2 * mean(est_day[wkd]) + 2 * mean(est_day[!wkd])
  # --- the shipped expansion: exact over the tally days ---
  rx <- estimate_comm_charter(dwg, Pc0)
  chk("census (none): default expansion is 'none' and the total is the exact sum over the tally days",
      identical(rx$census_expansion, "none") && isTRUE(all.equal(rx$Dungeness_Kept, obs)) && isTRUE(all.equal(rx$imputed_dung, 0)) && rx$n_unsampled_days == 4)
  chk("census (none): unsampled days carry zero and are labelled no-operation; the daily table still covers the calendar",
      nrow(rx$daily_full) == 14 && all(rx$daily_full$est_dung[!rx$daily_full$observed] == 0) &&
        all(grepl("no operation", rx$daily_full$source[!rx$daily_full$observed])) && isTRUE(all.equal(sum(rx$daily_full$est_dung), obs)))
  chk("census (none): no imputation variance; the per-vessel-mean term alone (zero here: constant catches)",
      isTRUE(all.equal(sum(rx$variance_detail$var_imputed), 0)) && isTRUE(all.equal(rx$Dungeness_Kept_var, 0)))
  chk("census (none): effort_total is the tally vessels", isTRUE(all.equal(rx$effort_total, sum(tally$commercial_tally + tally$charter_tally))))
  # --- the 2026-09-08 day-type expansion, kept as an option ---
  Pc <- modifyList(Pc0, list(census_expansion = "day_type"))
  r <- estimate_comm_charter(dwg, Pc)
  # 2026-09-11: observed_dung is now the DRAW FLOOR -- the commercial census on the tally
  # days plus the charter crab actually observed on the interviewed trips -- not the joint
  # tally-day estimate. The total and the imputed part are unchanged, and the total still
  # splits exactly into the two components.
  chk("census (day_type): total = observed + imputed; observed_dung is the commercial census + the observed charter crab",
      isTRUE(all.equal(r$Dungeness_Kept, obs + imp)) && isTRUE(all.equal(r$imputed_dung, imp)) &&
        isTRUE(all.equal(r$observed_dung, sum(tally$commercial_tally) * 40 + r$charter_observed_dung)) &&
        isTRUE(all.equal(r$commercial_dung + r$charter_dung, r$Dungeness_Kept)))
  v_exp <- 2^2 * var(est_day[wkd]) / 8 + 2^2 * var(est_day[!wkd]) / 2      # per-vessel means are exact here (zero variance)
  chk("census (day_type): imputation variance = sum_h (N_h - n_h)^2 s_h^2 / n_h", isTRUE(all.equal(r$Dungeness_Kept_var, v_exp)) && isTRUE(all.equal(r$Dungeness_Kept_se, sqrt(v_exp))))
  chk("census (day_type): daily table covers every calendar day, flags observed days", nrow(r$daily_full) == 14 && sum(r$daily_full$observed) == 10 && isTRUE(all.equal(sum(r$daily_full$est_dung), r$Dungeness_Kept)))
  chk("census (day_type): imputed days carry their stratum mean",
      isTRUE(all.equal(unique(r$daily_full$est_dung[!r$daily_full$observed & r$daily_full$day_type == "weekday"]), mean(est_day[wkd]))))
  # 2026-09-11: the default mode is "charter" (the charter expansion variance is carried,
  # the commercial census is not); "none" and "sampling" still work, "imputed_days" is the
  # 2026-09-08 alias for "sampling", and carried_var is what the drivers draw with.
  chk("census: default uncertainty mode is 'charter' (the charter expansion is carried, the commercial census is not)",
      identical(r$census_uncertainty, "charter") && isTRUE(all.equal(r$carried_var, r$charter_var)))
  chk("census: mode 'none' carries nothing, 'sampling' carries everything ('imputed_days' as its 2026-09-08 alias), bad values refused",
      isTRUE(all.equal(estimate_comm_charter(dwg, modifyList(Pc, list(census_uncertainty = "none")))$carried_var, 0)) &&
        { rs <- estimate_comm_charter(dwg, modifyList(Pc, list(census_uncertainty = "sampling")))
          identical(rs$census_uncertainty, "sampling") && isTRUE(all.equal(rs$carried_var, rs$Dungeness_Kept_var)) } &&
        identical(estimate_comm_charter(dwg, modifyList(Pc, list(census_uncertainty = "imputed_days")))$census_uncertainty, "sampling") &&
        inherits(tryCatch(estimate_comm_charter(dwg, modifyList(Pc, list(census_uncertainty = "bootstrap"))), error = function(e) e), "error") &&
        inherits(tryCatch(estimate_comm_charter(dwg, modifyList(Pc, list(census_expansion = "week"))), error = function(e) e), "error") &&
        inherits(tryCatch(estimate_comm_charter(dwg, modifyList(Pc, list(charter_expansion = "median"))), error = function(e) e), "error"))
  # a stratum with one sampled day borrows the pooled variance and says so
  dwg1 <- dwg; dwg1$comm_tally <- tally[c(1:8, 9), ]; dwg1$interview <- ints |> filter(event_date %in% dwg1$comm_tally$date)
  r1 <- estimate_comm_charter(dwg1, Pc)
  chk("census: a one-day stratum borrows the pooled between-day variance, flagged",
      any(grepl("pooled", r1$variance_detail$s2_source)) && is.finite(r1$Dungeness_Kept_se) && r1$Dungeness_Kept_se > 0)
  # per-window path sums the split and the variance
  Pw <- Pc; Pw$census_windows <- list("a" = c("2025-01-06", "2025-01-12"), "b" = c("2025-01-13", "2025-01-19"))
  rw <- estimate_comm_charter(dwg, Pw)
  chk("census: per-window path sums the two components, the variances and the daily table",
      isTRUE(all.equal(rw$commercial_dung + rw$charter_dung, rw$Dungeness_Kept)) && nrow(rw$daily_full) == 14 &&
        rw$Dungeness_Kept_var > 0 && isTRUE(all.equal(rw$Dungeness_Kept_se, sqrt(rw$Dungeness_Kept_var))) &&
        isTRUE(all.equal(sum(rw$daily_full$est_dung), rw$Dungeness_Kept)))
  chk("census: empty window returns the zero split", isTRUE(all.equal(estimate_comm_charter(list(comm_tally = tally[0, ], interview = ints), Pc)$Dungeness_Kept_se, 0)))
  # a day type with no sampled day no longer takes the component to NA: pooled mean, flagged
  dwg0 <- dwg; dwg0$comm_tally <- tally[1:8, ]; dwg0$interview <- ints |> filter(event_date %in% dwg0$comm_tally$date)
  r0 <- estimate_comm_charter(dwg0, Pc)
  chk("census: an unsampled day type takes the pooled sampled-day mean instead of NA",
      is.finite(r0$Dungeness_Kept) && isTRUE(all.equal(r0$Dungeness_Kept, sum(est_day[1:8]) + 2 * mean(est_day[1:8]) + 4 * mean(est_day[1:8]))) &&
        any(grepl("none", r0$variance_detail$s2_source)) && is.finite(r0$Dungeness_Kept_se))
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    # 2026-09-11: the drivers draw on carried_se (census_uncertainty decides what is in it),
    # clamped at observed_dung, and report the two components separately.
    chk(sprintf("%s: writes census_daily.csv / census_variance.csv and draws the census on carried_se, clamped at observed_dung", basename(drv)),
        any(grepl("census_daily.csv", d, fixed = TRUE)) && any(grepl("census_variance.csv", d, fixed = TRUE)) &&
          any(grepl("carried_se", d, fixed = TRUE)) && !any(grepl("\"imputed_days\"", d, fixed = TRUE)) &&
          any(grepl("observed_dung", d, fixed = TRUE)) && any(grepl("charter_se", d, fixed = TRUE)) &&
          any(grepl("commercial_dung", d, fixed = TRUE)))
  }
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("shipped: census_uncertainty = charter (2026-09-11: the charter expansion variance is carried, the commercial census is not)",
      identical(rc$census_uncertainty, "charter") && identical(rc$charter_expansion, "vessel"))
  chk("shipped: census_expansion = none (2026-09-09: unsampled days had no operation; the census is exact)", identical(rc$census_expansion, "none"))
})

# ---------------------------------------------------------------------------
# 54. The 2026-09-08 improvement ladder runner and the fit_agreement() helper. The
#     runner inherits every standing runner rule (sections 30, 31: DRY_RUN TRUE, tag
#     inside run_config, HTML moved, merge-by-key persistence); fit_agreement() is the
#     Monte-Carlo-error comparison the dynamic-f rung needs because bit-identity is
#     impossible once the parameter vector changes.
# ---------------------------------------------------------------------------
local({
  f <- "06_diagnostics/run_improvements_2026-09-08.R"
  chk("ladder runner present", file.exists(f)); if (!file.exists(f)) return(invisible(NULL))
  t <- readLines(f, warn = FALSE); tt <- t[!grepl("^\\s*#", t)]; s <- paste(tt, collapse = "\n")
  chk("ladder: the tag goes inside run_config and render() gets no output_dir",
      any(grepl("cfg$run_tag <- st$tag", tt, fixed = TRUE)) && !any(grepl("output_dir", tt[grepl("rmarkdown::render(", tt, fixed = TRUE)], fixed = TRUE)))
  chk("ladder: rendered HTML is moved into the run folder", any(grepl("file.copy(html", tt, fixed = TRUE)))
  chk("ladder: verdicts and the ladder table are MERGED by key, never appended",
      grepl('merge_csv_by(do.call(rbind, V), vp, c("stage", "criterion"))', s, fixed = TRUE) && grepl('merge_csv_by(do.call(rbind, LAD), lp, "rung")', s, fixed = TRUE))
  chk("ladder: every verdict block is wrapped so a reading defect cannot destroy a result", grepl(".safe <- function(sid, expr) tryCatch", s, fixed = TRUE))
  chk("ladder: the single-season window is pinned on every rung", grepl('season_filter = "2024-25"', s, fixed = TRUE) && grepl("modifyList(BASE, WINDOW, keep.null = TRUE)", s, fixed = TRUE))
  # 2026-09-10: the pin has to be COMPLETE. Four per-season keys were missing, so every
  # 2024-25 rung would have inherited the two-season config's 2023-24 values, pot_open_date
  # among them.
  for (k in c("pot_closure_start", "pot_closure_end", "pot_open_date", "census_start_date",
              "census_end_date", "commercial_opener", "pe_empty_effort_stratum",
              "pe_empty_stratum", "pe_variance", "bss_seed"))
    chk(sprintf("ladder: WINDOW pins %s", k), grepl(paste0(k, " = "), s, fixed = TRUE))
  chk("ladder: pot_open_date is pinned to the 2024-25 value, not the 2-season 2023-12-01",
      grepl('pot_open_date = "2024-12-01"', s, fixed = TRUE))
  chk("ladder: F_METHOD exists and ships 'new_throughout' (Matt 2026-09-10: the run uses the new f)",
      any(grepl('^F_METHOD <- "new_throughout"', t)))
  chk("ladder: the preflight PROVES no non-delta key differs between rungs, before any MCMC",
      grepl("no configuration leak", s, fixed = TRUE) && grepl("DELTA_KEYS", s, fixed = TRUE))
  chk("ladder: RESUME requires a matching config digest, not just a finished folder",
      grepl("stage_digest", s, fixed = TRUE) && grepl("IMP_STAGE.txt", s, fixed = TRUE) &&
      grepl(".stage_stamp(done, sid)", s, fixed = TRUE))
  chk("ladder: a manifest of the resolved rungs is written before the fits",
      grepl("improvements_2026-09-08_manifest", s, fixed = TRUE))
  chk("ladder: the mode is in every tag and output filename, so two modes cannot collide",
      grepl('.sfx <- if (identical(F_METHOD, "ladder")) "" else "-newf"', s, fixed = TRUE) &&
      grepl('improvements_2026-09-08_verdicts%s.csv', s, fixed = TRUE))
  chk("ladder: R0 runs the PE under every unsampled-cell arm, so D19 needs no refit",
      grepl("pe_unsampled_cell_arms.csv", s, fixed = TRUE) && grepl("SHIPPED", s, fixed = TRUE))
  # the deltas are cumulative and each adds one thing
  e <- new.env()
  e$`%||%` <- function(a, b) if (is.null(a)) b else a
  # 2026-09-11: take the WHOLE contiguous D_R1..D_R4 block and extend it until it parses,
  # instead of grepping for lines that look like delta keys. The grep broke the moment a
  # delta gained a comment line or a key outside its list, which is the harness-fragility
  # failure mode section 43 exists to prevent.
  # 2026-09-10: the block now starts at F_NEW (D_R1 is built from it) and F_METHOD has to
  # be defined in the eval environment, because the deltas branch on it.
  e$F_METHOD <- "new_throughout"; e$modifyList <- modifyList
  .i1 <- grep("^F_NEW *<-", tt)[1]; .i2 <- grep("^D_R4 *<-", tt)[1]
  stopifnot(is.finite(.i1), is.finite(.i2), .i2 > .i1)
  .parsed <- FALSE
  for (.end in .i2:min(.i2 + 12L, length(tt))) {
    .txt <- paste(tt[.i1:.end], collapse = "\n")
    if (!inherits(try(parse(text = .txt), silent = TRUE), "try-error")) { eval(parse(text = .txt), envir = e); .parsed <- TRUE; break }
  }
  chk("ladder: the F_NEW..D_R4 delta block parses as written", .parsed)
  # 2026-09-10: D_R1's f block now depends on F_METHOD, so the block is evaluated under
  # BOTH modes and each is asserted against its own contract.
  chk("ladder: R1 keeps the pre-patch TURNOVER priors under either F_METHOD (tau 1.2, sigma 0.3)",
      identical(e$D_R1$tau_boat_prior_mu, 1.2) && identical(e$D_R1$tau_boat_prior_sigma, 0.3) &&
      is.null(e$D_R1$shared_tau_sigma) && identical(e$D_R1$tau_shore_prior_mu, 1.7))
  chk("ladder: under new_throughout every rung carries the NEW f from R1 onward",
      identical(e$D_R1$crab_fraction_strata, "month") && identical(e$D_R1$crab_fraction_source, "both") &&
      identical(e$D_R1$crab_fraction_dynamic, TRUE))
  chk("ladder: F_OLD and F_NEW are the only two f blocks, named once each",
      identical(e$F_OLD, list(crab_fraction_strata = "none", crab_fraction_source = "ie", crab_fraction_dynamic = FALSE)) &&
      identical(e$F_NEW, list(crab_fraction_strata = "month", crab_fraction_source = "both", crab_fraction_dynamic = TRUE)))
  chk("ladder: under F_METHOD = 'ladder' R1 rolls the f block back to the retired one",
      { e2 <- new.env(); e2$`%||%` <- e$`%||%`; e2$F_METHOD <- "ladder"
        e2$modifyList <- modifyList; eval(parse(text = .txt), envir = e2)
        identical(e2$D_R1$crab_fraction_strata, "none") && identical(e2$D_R1$crab_fraction_dynamic, FALSE) &&
        identical(e2$D_R3$crab_fraction_dynamic, TRUE) })
  chk("ladder: R2f is R2 with the f block rolled back (the factorization control)",
      identical(e$D_R2f$tau_boat_prior_mu, "calibration") && identical(e$D_R2f$crab_fraction_dynamic, FALSE) &&
      identical(e$D_R2f$crab_fraction_strata, "none") &&
      identical(e$D_R2f[setdiff(names(e$D_R2f), names(e$F_OLD))], e$D_R2[setdiff(names(e$D_R2), names(e$F_NEW))]))
  chk("ladder: R2 adds only the calibration prior", identical(e$D_R2$tau_boat_prior_mu, "calibration") &&
      identical(e$D_R2$crab_fraction_strata, e$D_R1$crab_fraction_strata))
  # Under new_throughout R3a/R3 are refused by the stage guard, so their contract is
  # asserted in the ladder-mode environment where they are the f rungs.
  chk("ladder: R3a adds only the monthly f from both sources (legacy construction)",
      { e2 <- new.env(); e2$`%||%` <- e$`%||%`; e2$F_METHOD <- "ladder"; e2$modifyList <- modifyList
        eval(parse(text = .txt), envir = e2)
        identical(e2$D_R3a$crab_fraction_strata, "month") && identical(e2$D_R3a$crab_fraction_source, "both") &&
        identical(e2$D_R3a$crab_fraction_dynamic, FALSE) })
  chk("ladder: R3 adds only the dynamic f and still pins the pre-adoption shore turnover (1.7 / 0.3)",
      identical(e$D_R3$crab_fraction_dynamic, TRUE) && identical(e$D_R3$tau_shore_prior_mu, 1.7) && identical(e$D_R3$tau_shore_prior_sigma, 0.3))
  chk("ladder: the stage guards refuse a mode/stage mix rather than running an unattributable rung",
      grepl("R2f is the new_throughout factorization control", s, fixed = TRUE) &&
      grepl("R3a / R3 are the 'ladder' f rungs", s, fixed = TRUE))
  chk("ladder: R4 adds only the derived shore turnover, and equals the shipped run_config on the moved keys",
      identical(e$D_R4$tau_shore_prior_mu, "derived") && identical(e$D_R4$crab_fraction_dynamic, TRUE) && {
        rc <- new.env(); sys.source("run_config.R", envir = rc); rc <- rc$run_config
        all(vapply(c("tau_boat_prior_mu", "tau_boat_prior_sigma", "shared_tau_sigma", "crab_fraction_strata", "crab_fraction_source",
                     "crab_fraction_dynamic", "tau_shore_prior_mu", "tau_shore_prior_sigma", "census_uncertainty"),
                   function(k) identical(e$D_R4[[k]], rc[[k]]), logical(1))) })
  chk("ladder: the gear cross-check follows the shipped rung", any(grepl('^GEAR_FOLLOWS <- "R4"', t)))
  # fit_agreement(): two synthetic run folders
  source("03_R_functions/batch_verdict_helpers.R")
  mk <- function(dir, means, se) {
    dir.create(dir, showWarnings = FALSE)
    utils::write.csv(data.frame(mean = means, se_mean = se, sd = se * 10, row.names = c("B1", "mu_mu_E[1]", "f_crab_out[1]", "sigma_f_out")),
                     file.path(dir, "bss_full_summary_private_boat_all_gear_Dungeness_Kept.csv"))
  }
  da <- tempfile("fa_a"); db <- tempfile("fa_b")
  mk(da, c(0.50, 2.70, 0.30, 0.00), c(0.001, 0.02, 0.001, 0))
  mk(db, c(0.501, 2.71, 0.90, 0.66), c(0.001, 0.02, 0.002, 0.003))
  r <- fit_agreement(db, da, pat = "private_boat", exclude = "^(f_crab|sigma_f)")
  chk("fit_agreement: non-excluded rows within MC error -> PASS, excluded rows ignored, zero-se rows skipped",
      identical(r$verdict, "PASS") && length(r$z) == 2 && !any(grepl("f_crab|sigma_f", names(r$z))))
  r2 <- fit_agreement(db, da, pat = "private_boat")
  chk("fit_agreement: a moved parameter that is NOT excluded fails the test", identical(r2$verdict, "FAIL") && max(r2$z) > 100)
  chk("fit_agreement: missing folder -> REVIEW, not an error", identical(fit_agreement(tempfile(), da)$verdict, "REVIEW"))
  unlink(c(da, db), recursive = TRUE)
})

# ---------------------------------------------------------------------------
# 55. Trip types (2026-09-09). The interview workbook is rebuilt from the raw export by
#     04_input_files/build_interview_combined.R and carries trip_type / trip_type_class /
#     creel_area / interview_time; the contact builder reads the trip type where it
#     exists (falls back to crabbers > 0), restricts to the launch sites, and counts the
#     combos that observe the combo-trip share c; the Stan data carry the per-day combo
#     rows and the c walk's priors.
# ---------------------------------------------------------------------------
local({
  # the builder's classifier, pure
  e <- new.env()
  src <- readLines("04_input_files/build_interview_combined.R", warn = FALSE)
  i1 <- grep("^trip_type_class <- function", src)[1]; i2 <- i1 + which(src[i1:length(src)] == "}")[1] - 1L
  eval(parse(text = src[i1:i2]), envir = e)
  cls <- e$trip_type_class(c("Crab Only", "Salmon & Crab", "Bottomfish & Crab", "Halibut & Crab", "Tuna & Crab", "Finfish Only", "Non Fishing Trip", NA, "", "Kayak"))
  chk("trip types: the export labels map to crab_only / combo / other_fishery / non_fishing, blanks to NA, unknowns flagged",
      identical(cls, c("crab_only", "combo", "combo", "combo", "combo", "other_fishery", "non_fishing", NA, NA, "unclassified")))
  # the rebuilt workbook
  wb <- "04_input_files/interview_combined.xlsx"
  chk("workbook: interview_combined.xlsx carries the 17 legacy columns plus trip_type, trip_type_class, creel_area, interview_time", {
    nm <- names(readxl::read_excel(wb, sheet = "data", n_max = 2))
    all(c("season","creel_location","date","survey_id","interview_num","crabbing_mode","boat_type","crabbers","gear_type","number_of_gear",
          "dungeness_kept","red_rock_kept","hours_fished","crabber_hours","gear_hours","completed_trip","gear_tampered",
          "trip_type","trip_type_class","creel_area","interview_time") %in% nm) })
  # 2026-09-10: the pasted exports were replaced by the per-season workbooks (section 57)
  chk("workbook: the season workbooks and the interview / shift builders are in the repository",
      length(list.files("04_input_files/raw", pattern = "^[0-9]{4}_rec_crab_harvest_data\\.xlsx$")) >= 4 &&
        !file.exists("04_input_files/raw/interviewdata20222026.xlsx") &&
        file.exists("04_input_files/build_interview_combined.R") && file.exists("04_input_files/build_sampler_shifts.R"))
  # the contact builder on a synthetic frame with trip types and areas
  source("03_R_functions/fetch_crab_data.R")
  fr <- tibble(population = rep("private_boat", 8),
               event_date = as.Date(c(rep("2025-08-01", 5), rep("2025-08-02", 3))),
               creel_area = c(rep("Westport Boat Launch", 4), "Westport Marina", rep("Westport Boat Launch", 3)),
               crabbers   = c(2, 0, 3, 0, 2,   1, 0, NA),
               trip_type_class = c("combo", "other_fishery", "crab_only", "non_fishing", "crab_only",   "crab_only", "non_fishing", NA),
               interview_time = c("10:30", "11:00", "11:15", "12:00", "12:30", "09:50", "10:10", "13:00"),
               survey_id = "S1")
  Pa <- list(boat_launch_areas = c("Westport Boat Launch", "Ocean Shores Boat Launch"))
  bc <- boat_contacts_from_interviews(fr, Pa)
  d1 <- bc[bc$event_date == as.Date("2025-08-01"), ]; d2 <- bc[bc$event_date == as.Date("2025-08-02"), ]
  chk("contacts: the marina row is excluded by default (launch sites only), crabbing = crab_only + combo, combos counted",
      d1$boats_total == 4 && d1$boats_crabbing == 2 && d1$boats_combo == 1 && d1$boats_crab_only == 1 && d1$boats_other == 1 && d1$boats_nonfishing == 1 && d1$boats_typed == 4)
  chk("contacts: an untyped row with NA crabbers is unclassified and counts nowhere", d2$boats_total == 2 && d2$boats_crabbing == 1 && d2$boats_typed == 2)
  chk("contacts: crab_fraction_contact_areas = 'all' keeps every private-boat row", boat_contacts_from_interviews(fr, list(crab_fraction_contact_areas = "all"))$boats_total[1] == 5)
  chk("contacts: without a trip-type column the crabbers > 0 rule applies and combos are unknown (0 typed)",
      { b0 <- boat_contacts_from_interviews(fr |> select(-trip_type_class), Pa); b0$boats_crabbing[1] == 2 && b0$boats_typed[1] == 0 })
  det <- boat_contact_detail(fr)
  chk("contact detail: one row per private boat with the contact hour parsed", nrow(det) == 8 && isTRUE(all.equal(det$contact_hour[1], 10.5)))
  # the combo rows reach the Stan data and the c walk switches on
  dwg <- list(boat_contacts = bc)
  rows <- crab_fraction_source_rows(dwg, tibble(x = 1), list(crab_fraction_source = "interviews"), quiet = TRUE)
  chk("source rows carry the typed counts", all(c("boats_typed", "boats_crab_only", "boats_combo") %in% names(rows)) && sum(rows$boats_combo) == 1)
  Pd <- modifyList(P, list(crab_fraction_strata = "month", crab_fraction_dynamic = TRUE, use_osp_crab_lower = FALSE))
  Pd$crab_fraction_rows <- rows
  cf <- crab_fraction_stan_data(FALSE, mkdays("2025-08-01", 31), Pd, quiet = TRUE)
  chk("dynamic + typed contacts: combo_dynamic = 1, one CFC row per day with crabbing boats, combo <= crabbing",
      cf$combo_dynamic == 1L && cf$CFC_n == 2 && sum(cf$cfc_crab) == 3 && sum(cf$cfc_combo) == 1 && all(cf$cfc_combo <= cf$cfc_crab))
  chk("dynamic + typed contacts: the audit table carries the typed counts", all(c("typed_days", "typed_crabbing", "typed_combo") %in% names(attr(cf, "f_strata"))) && attr(cf, "f_strata")$typed_combo[1] == 1)
  chk("dynamic + typed contacts: a day with typed boats but no crabbing boat carries no combo row",
      { r2 <- rows; r2$boats_crab_only <- 0; r2$boats_combo <- 0; Pd2 <- Pd; Pd2$crab_fraction_rows <- r2
        crab_fraction_stan_data(FALSE, mkdays("2025-08-01", 31), Pd2, quiet = TRUE)$CFC_n == 0 })
  chk("legacy (key absent): the combo fields are present and inert", crab_fraction_stan_data(FALSE, mkdays("2025-08-01", 31), modifyList(Pd, list(crab_fraction_dynamic = NULL)), quiet = TRUE)$combo_dynamic == 0L)
  for (m in c("02_stan_models/crab_bss_pooled.stan", "02_stan_models/crab_bss_gear_resolved.stan")) {
    src <- paste(readLines(m, warn = FALSE), collapse = "\n")
    chk(sprintf("%s: the OSP stream reads f x (1 - c[k]) and the combo rows are beta-binomial on the typed crabbing boats", basename(m)),
        # D40: the OSP day's share is f_osp (f_crab[k] when the volume term is off)
        grepl("real f_osp = f_crab[osp_f_stratum[i]];", src, fixed = TRUE) &&
          grepl("p_osp = f_osp * (1 - combo_c[osp_f_stratum[i]]);", src, fixed = TRUE) &&
          grepl("cfc_combo[i] ~ beta_binomial(cfc_crab[i]", src, fixed = TRUE) && !grepl("combo_a", src, fixed = TRUE))
  }
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    chk(sprintf("%s: reports the combo walk and writes c and f(1-c) into the strata file; runs the shift diagnostic", basename(drv)),
        any(grepl("sigma_c_out", d, fixed = TRUE)) && any(grepl("crab_only_share_median", d, fixed = TRUE)) &&
          any(grepl("diagnose_shift_coverage(sampler_shifts, ie_data, dwg$boat_contacts_detail", d, fixed = TRUE)))
  }
  e2 <- new.env(); sys.source("run_config.R", envir = e2); rc <- e2$run_config
  chk("shipped: the combo walk priors and the launch-site restriction exist; no shift weighting yet",
      identical(rc$crab_fraction_combo_level, 0.3) && is.null(rc$crab_fraction_contact_areas) && identical(rc$crab_fraction_shift_weighting, "none"))
})

# ---------------------------------------------------------------------------
# 56. Sampler shifts (2026-09-09): the formatted workbook, its reader, and the shift
#     coverage diagnostic on synthetic inputs.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/sampler_shifts.R")
  wb <- "04_input_files/sampler_shifts.xlsx"
  chk("shifts workbook: present with the documented columns", file.exists(wb) && {
    nm <- names(readxl::read_excel(wb, sheet = "data", n_max = 2))
    all(c("survey_id","survey_num","season","date","day_of_week","holiday","creel_location","samplers","special_conditions",
          "check_in","check_out","check_in_hour","check_out_hour","shift_hours","qc_flag","notes") %in% nm) })
  sh <- readxl::read_excel(wb, sheet = "data")
  chk("shifts workbook: survey ids link to the interview workbook's survey_id ('S<n>'), dates are ISO text, flagged rows keep NA hours",
      all(grepl("^S[0-9]+$", sh$survey_id[!is.na(sh$survey_id)])) && all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", sh$date)) &&
        all(is.na(sh$shift_hours[nzchar(sh$qc_flag) & !is.na(sh$qc_flag)])))
  # the reader's window filter and the union / coverage arithmetic
  Ps <- list(gh_creel_location = "Grays Harbor", season_filter = "2024-25", est_date_start = "2024-09-16", est_date_end = "2025-09-15")
  got <- fetch_sampler_shifts(Ps, quiet = TRUE)
  chk("shifts reader: Grays Harbor 2024-25 window has ~200 shift days at about 09:45-15:50",
      n_distinct(got$event_date) >= 150 && abs(median(got$check_in_hour, na.rm = TRUE) - 9.75) < 0.5 && abs(median(got$check_out_hour, na.rm = TRUE) - 15.8) < 0.5)
  u <- .shift_union(c(9, 12, 15), c(11, 13, 16))
  chk("shift union: overlapping windows merge, disjoint ones stay separate", nrow(u) == 3 && isTRUE(all.equal(.hours_covered(u), 2 + 1 + 1)))
  u2 <- .shift_union(c(9, 10), c(12, 14))
  chk("shift union: nested / overlapping windows form one interval", nrow(u2) == 1 && isTRUE(all.equal(u2[1, ], c(9, 14))))
  chk("share inside: a return profile is apportioned by the windows", isTRUE(all.equal(.share_inside(c(8, 10, 13, 17), c(1, 1, 1, 1), u2), 0.5)))
  # a synthetic diagnostic run
  ie <- tibble(x = 1)
  attr(ie, "ie_boat_intervals") <- tibble(event_date = as.Date("2025-05-01"), season = "2024-25", day_type = "weekday", location_name = "WBL",
                                          hour = c(8, 10, 12, 14, 16, 18), boats_in = c(5, 3, 1, 0, 0, 0), boats_out = c(0, 1, 3, 3, 2, 1))
  shifts <- tibble(survey_id = c("S1", "S2"), event_date = as.Date(c("2025-05-01", "2025-05-02")), check_in_hour = c(9.5, 9.75),
                   check_out_hour = c(15.5, 16), shift_hours = c(6, 6.25), qc_flag = c("", ""), creel_location = "Grays Harbor")
  det <- tibble(event_date = as.Date(c("2025-05-01", "2025-05-01", "2025-05-02")), creel_area = "Westport Boat Launch", survey_id = c("S1", "S1", "S2"),
                trip_type_class = c("crab_only", "other_fishery", "combo"), crabbers = c(2, 0, 3), contact_hour = c(10, 14, 17))
  sc <- diagnose_shift_coverage(shifts, ie, det, list(ie_boat_location = "WBL", shift_coverage_ie_min_days = 1), output_dir = NULL, quiet = TRUE)
  chk("shift coverage: per-day rows, the return share inside the shift, and the contact-hour table",
      nrow(sc$by_day) == 2 && isTRUE(all.equal(sc$by_day$return_share_covered[1], (1 + 3 + 3) / 10)) &&
        !is.null(sc$contact_hours) && sum(sc$contact_hours$contacts) == 3 && sum(sc$contact_hours$inside_shift) == 2)
  chk("shift coverage: falls back to pooled sites when the WBL series is thin, and says so",
      grepl("pooled", diagnose_shift_coverage(shifts, ie, det, list(ie_boat_location = "WBL", shift_coverage_ie_min_days = 5), quiet = TRUE)$summary$return_profile_source))
  chk("shift coverage: no shifts -> NULL, no error", is.null(diagnose_shift_coverage(NULL, ie, det, list(), quiet = TRUE)))
})

# ---------------------------------------------------------------------------
# 57. The inputs rebuilt from the per-season creel workbooks (2026-09-10). Every model
#     workbook is built from raw/<YYYY><YY>_rec_crab_harvest_data.xlsx by the builders in
#     04_input_files/ (build_all_inputs.R); the pasted exports of 2026-09-09 are gone. The
#     checks: the builder helpers' parsers, the workbooks' shape and coverage (four
#     seasons through 2026-09-08; the 2023-24 gear count back; gear labels harmonised;
#     the effort qc_flag; the two-season tally; the charter roster; the five-season
#     holiday calendar), the reader that guesses column types over the whole column, the
#     effort qc drop, and the charter-roster frame of the census on a synthetic fixture.
# ---------------------------------------------------------------------------
local({
  # --- the builder helpers, pure ---
  eh <- new.env(); sys.source("04_input_files/build_helpers.R", envir = eh)
  chk("helpers: dates parse from Excel serial text, ISO text and M/D/YYYY text; junk is NA",
      identical(eh$parse_date_any(c("45292", "2024-01-01", "2024-01-01 00:00:00", "1/1/2024", "yesterday", NA)),
                as.Date(c("2024-01-01", "2024-01-01", "2024-01-01", "2024-01-01", NA, NA))))
  chk("helpers: clock times parse from a day fraction and from H:MM(:SS) text; 12-hour slips are not corrected",
      isTRUE(all.equal(eh$parse_clock_hours(c("0.4479166666666667", "10:45:00", "9:05", "25:00", "abc")), c(10.75, 10.75, 9 + 5/60, NA, NA))) &&
        identical(eh$hhmm(c(10.75, NA)), c("10:45", NA)) && identical(eh$hhmmss(10.75), "10:45:00"))
  chk("helpers: the fishery season turns on Sep 16; survey ids are 'S<n>'",
      identical(eh$season_of(as.Date(c("2024-09-15", "2024-09-16", "2025-01-01"))), c("2023-24", "2024-25", "2024-25")) &&
        identical(eh$survey_id_of(c("S4949", "4949", "4949.0", NA)), c("S4949", "S4949", "S4949", NA)))
  chk("helpers: the 2022-23 site names map to the current vocabulary; Westport sites are untouched",
      identical(eh$harmonise_area(c("Tokeland boat Launch", "Chinook Boat Launch", "Westport Boat Launch ", "Westport Docks Float 20")),
                c("Tokeland Boat Launch", "Chinook Boat Launch and Marina", "Westport Boat Launch", "Westport Docks Float 20")))
  # the gear vocabulary map lives in the interview builder; pull the definitions out
  src <- readLines("04_input_files/build_interview_combined.R", warn = FALSE)
  eg <- new.env(); for (nm in c("library(stringr)", "library(dplyr)")) eval(parse(text = nm), envir = eg)
  i1 <- grep("^gear_map <- c\\(", src)[1]; i2 <- grep("^harmonise_gear <- function", src)[1]; i3 <- i2 + which(src[i2:length(src)] == "}")[1] - 1L
  eval(parse(text = src[i1:i3]), envir = eg)
  hg <- eg$harmonise_gear(c("Collapsible trap or ring", "Pot, Collapsible trap or ring, Fishing rod with snare", "Fishing rod with foldable trap",
                            "Rake, net or hands", "Star trap, Fishing rod with snare", "Trap (foldable, star), Snare", "Pot", "Slip ring pot", NA, ""))
  chk("gear labels: the 2022-24 vocabulary maps onto the 2024-25 labels, multi-gear lists keep their order, current labels are unchanged",
      identical(hg, c("Ring Net", "Pot, Ring Net, Snare", "Trap (foldable, star)", "Rake or Net", "Trap (foldable, star), Snare",
                      "Trap (foldable, star), Snare", "Pot", "Slip Ring Pot", NA, NA)))
  # the gear-resolved regex on the harmonised labels: a ring net is a ring net, not a trap
  rn <- function(x) c(pot = str_detect(x, "(?i)\\bpot\\b") & !str_detect(x, "(?i)\\bslip\\s*ring\\b"), ring = str_detect(x, "(?i)\\bring\\s*net\\b"),
                      trap = str_detect(x, "(?i)\\b(trap|star)\\b"), snare = str_detect(x, "(?i)\\bsnare\\b"))
  chk("gear labels: on the old label 'Collapsible trap or ring' the classifier saw a TRAP; on 'Ring Net' it sees a ring net",
      isTRUE(rn("Collapsible trap or ring")[["trap"]]) && !isTRUE(rn("Collapsible trap or ring")[["ring"]]) &&
        isTRUE(rn("Ring Net")[["ring"]]) && !isTRUE(rn("Ring Net")[["trap"]]))

  # --- the workbooks ---
  wbs <- list.files("04_input_files/raw", pattern = "^[0-9]{4}_rec_crab_harvest_data\\.xlsx$")
  chk("raw: four season workbooks (2223, 2324, 2425, 2526) and no pasted export", setequal(substr(wbs, 1, 4), c("2223", "2324", "2425", "2526")) &&
        !file.exists("04_input_files/raw/surveydata20222026.xlsx"))
  chk("builders: six builders, the shared helpers and the orchestrator are in 04_input_files",
      all(file.exists(file.path("04_input_files", c("build_helpers.R", "build_all_inputs.R", "build_interview_combined.R", "build_effort_combined.R",
                                                   "build_sampler_shifts.R", "build_comm_charter_tally.R", "build_charter_trips.R", "build_crabbing_holidays.R")))))
  iv <- read_input_workbook("interview_combined.xlsx")
  chk("interviews: four seasons, through 2026-09-08 or later, 31 columns with the ten 2026-09-10 additions",
      setequal(unique(iv$season), c("2022-23", "2023-24", "2024-25", "2025-26")) && max(iv$date) >= "2026-09-08" &&
        all(c("gear_type_raw", "boat_name", "bay_or_ocean", "river_or_ocean", "total_vehicles", "crab_released", "dungeness_returned",
              "dungeness_returned_reason", "red_rock_returned", "red_rock_returned_reason") %in% names(iv)) && ncol(iv) == 31)
  chk("interviews: the 2023-24 gear count is present (D15 closed): every 2023-24 row has number_of_gear, 0 to 23",
      sum(!is.na(iv$number_of_gear[iv$season == "2023-24"])) == sum(iv$season == "2023-24") && max(iv$number_of_gear[iv$season == "2023-24"]) == 23)
  chk("interviews: gear_type carries the 2024-25 vocabulary in every season, gear_type_raw the label as recorded",
      !any(grepl("Collapsible trap or ring|Fishing rod with", iv$gear_type)) && any(grepl("Collapsible trap or ring", iv$gear_type_raw)) &&
        all(iv$gear_type[iv$season == "2025-26" & !is.na(iv$gear_type)] == iv$gear_type_raw[iv$season == "2025-26" & !is.na(iv$gear_type)]))
  chk("interviews: read_input_workbook types the late-starting columns from the whole column (total_vehicles numeric, not logical)",
      is.numeric(iv$total_vehicles) && is.numeric(iv$gear_tampered) && is.character(iv$interview_time))
  chk("interviews: dates are ISO text, interview_time is HH:MM, completed_trip is 0/1/NA, creel_area has no trailing space",
      all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", iv$date)) && all(grepl("^[0-9]{2}:[0-9]{2}$", iv$interview_time[!is.na(iv$interview_time)])) &&
        all(iv$completed_trip[!is.na(iv$completed_trip)] %in% c(0, 1)) && !any(grepl("\\s$", iv$creel_area)))
  ef <- read_input_workbook("effort_combined.xlsx")
  chk("effort: four seasons, the count columns of every protocol, count_time HH:MM:SS, a qc_flag column",
      setequal(unique(ef$season), c("2022-23", "2023-24", "2024-25", "2025-26")) && max(ef$date) >= "2026-09-08" &&
        all(c("total_gear_count", "boat_trailer_count", "vehicle_count", "boats_entering_marina", "buoy_count", "crabber_count", "jetty_people_count", "qc_flag") %in% names(ef)) &&
        all(grepl("^[0-9]{2}:[0-9]{2}:[0-9]{2}$", ef$count_time[!is.na(ef$count_time)])))
  chk("effort: the 2024-25 rows reproduce the 2026-07-16 workbook (3,256 rows; trailer and gear sums)",
      sum(ef$season == "2024-25") == 3256 && sum(ef$boat_trailer_count[ef$season == "2024-25" & ef$creel_area == "Westport Boat Launch"]) == 2558 &&
        sum(ef$total_gear_count[ef$season == "2024-25" & ef$creel_area == "Westport Docks Float 20"]) == 14676)
  chk("effort: the Feb 2024 interview-total rows are flagged (20, all 2023-24) and the 2023-24 trailer counts come from the Truck/Boat Trailer column",
      sum(ef$qc_flag %in% "gear_count_from_interviews") == 20 && all(ef$season[ef$qc_flag %in% "gear_count_from_interviews"] == "2023-24") &&
        sum(ef$boat_trailer_count[ef$season == "2023-24" & ef$creel_area == "Westport Boat Launch"], na.rm = TRUE) == 3603)
  sh <- read_input_workbook("sampler_shifts.xlsx")
  chk("shifts: four seasons with the sampler's on-site conditions (tide, rain, weather, wind, wind_direction, weather_location)",
      setequal(unique(sh$season), c("2022-23", "2023-24", "2024-25", "2025-26")) && all(c("tide", "rain", "weather", "wind", "wind_direction", "weather_location", "holiday") %in% names(sh)) &&
        nrow(sh) == 3077)
  ta <- read_input_workbook("wes_commercial_tally.xlsx")
  chk("tally: 2024-25 (47 days, unchanged sums) and 2025-26 (23 days, Dec 1 to Jan 3), read by header name not position",
      sum(ta$season == "2024-25") == 47 && sum(ta$private_tally[ta$season == "2024-25"]) == 107 && sum(ta$commercial_tally[ta$season == "2024-25"]) == 164 &&
        sum(ta$charter_tally[ta$season == "2024-25"]) == 23 && sum(ta$season == "2025-26") == 23 && min(ta$date[ta$season == "2025-26"]) == "2025-12-01" &&
        max(ta$date[ta$season == "2025-26"]) == "2026-01-03" && sum(ta$commercial_tally[ta$season == "2025-26"]) == 67 && sum(ta$private_tally[ta$season == "2025-26"]) == 85)
  cr <- read_input_workbook("charter_trips.xlsx")
  w25 <- cr |> filter(season == "2024-25", port == "Westport", date >= "2024-12-03", date <= "2025-02-08", status != "canceled")
  chk("charter roster: 2024-25 Westport has 31 sailed trips in the tally window, 8 of them on days without a tally",
      nrow(w25) == 31 && sum(!w25$date %in% ta$date) == 8 && all(w25$status[!w25$date %in% ta$date] == "missed"))
  ho <- read_input_workbook("crabbing_holidays.xlsx")
  chk("holidays: five seasons by one rule; the 2024-25 rows are the 2026-07-16 ten plus the three added on 2026-09-11; observed days present",
      setequal(unique(ho$season), c("2022-23", "2023-24", "2024-25", "2025-26", "2026-27")) &&
        setequal(ho$date[ho$season == "2024-25"],
                 c("2024-11-29", "2024-12-31", "2025-01-01", "2025-02-08", "2025-05-24", "2025-05-25", "2025-05-26", "2025-06-15", "2025-07-04", "2025-09-01",
                   "2024-11-28", "2024-11-11", "2025-06-19")) &&
        "2026-07-03" %in% ho$date && "2024-01-01" %in% ho$date && "2023-11-24" %in% ho$date &&
        # Thanksgiving is the 4th Thursday; Veterans Day and Juneteenth carry their observed day when the date is a weekend
        all(weekdays(as.Date(ho$date[ho$holiday_name == "Thanksgiving Day"])) == "Thursday") &&
        "2023-11-10" %in% ho$date[ho$holiday_name == "Veterans Day (observed)"] &&
        "2027-06-18" %in% ho$date[ho$holiday_name == "Juneteenth (observed)"] &&
        # the sampler-flagged days the counts do NOT support stay out
        !any(c("2024-12-24", "2025-01-20", "2025-02-17") %in% ho$date))
  source("03_R_functions/read_crabbing_holidays.R")
  chk("holidays reader: every season resolves, and a two-season vector resolves to both calendars",
      all(vapply(c("2022-23", "2023-24", "2024-25", "2025-26"), function(sn) length(read_crabbing_holidays(list(season_filter = sn))) >= 13, logical(1))) &&
        length(read_crabbing_holidays(list(season_filter = c("2023-24", "2024-25")))) ==
          length(read_crabbing_holidays(list(season_filter = "2023-24"))) + length(read_crabbing_holidays(list(season_filter = "2024-25"))))

  # --- the readers on the rebuilt workbooks ---
  source("03_R_functions/fetch_crab_data.R"); source("03_R_functions/validate_season_window.R")
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  P1 <- modifyList(rc, list(season_filter = "2023-24", est_date_start = "2023-09-16", est_date_end = "2024-09-15", pot_closures = NULL, census_windows = NULL))
  out <- capture.output(d1 <- fetch_crab_data(P1))
  chk("reader 2023-24: effort counts and interviews load; the 20 interview-total rows are held out; no Float 20 count on Feb 3 to 13, 2024 enters the shore series",
      any(grepl("20 row\\(s\\) held out by qc_flag", out)) && nrow(d1$shore_effort) > 300 && nrow(d1$boat_effort) > 300 &&
        !any(d1$shore_effort$event_date %in% (as.Date("2024-02-03") + 0:10) & d1$shore_effort$f20_gear > 30))
  chk("reader 2023-24: the interviews carry a gear count (the 2023-24 CPUE denominator exists again)", mean(!is.na(d1$interview$number_of_gear)) > 0.99)
  chk("reader: effort_qc_drop = \"none\" keeps every row", { P0 <- modifyList(P1, list(effort_qc_drop = "none")); o0 <- capture.output(d0 <- fetch_crab_data(P0)); nrow(d0$shore_effort) > nrow(d1$shore_effort) })
  P2 <- modifyList(rc, list(season_filter = "2025-26", est_date_start = "2025-09-16", est_date_end = "2026-09-15", pot_closures = NULL, census_windows = NULL,
                            pot_closure_start = "2025-09-16", pot_closure_end = "2025-11-30", pot_open_date = "2025-12-01",
                            census_start_date = "2025-12-01", census_end_date = "2026-01-03"))
  out2 <- capture.output(d2 <- fetch_crab_data(P2))
  chk("reader 2025-26: the season loads with shore and boat counts, contacts, a 23-day tally and the charter roster",
      nrow(d2$shore_effort) > 500 && nrow(d2$boat_effort) > 400 && sum(d2$boat_contacts$boats_total) > 900 && !is.null(d2$charter_roster) && nrow(d2$charter_roster) == 11)
  chk("reader: the charter roster is filtered to the port and the season; absent file -> NULL",
      all(d2$charter_roster$port == "Westport") && all(d2$charter_roster$season == "2025-26") &&
        is.null(fetch_charter_roster(list(charter_trips_file = "no_such_file.xlsx"), quiet = TRUE)))

  # --- the charter-roster frame of the census, on the section-53 fixture plus a roster ---
  source("03_R_functions/estimate_comm_charter.R")
  cal <- seq(as.Date("2025-01-06"), as.Date("2025-01-19"), by = "day"); wk <- weekdays(cal) %in% c("Saturday", "Sunday")
  samp <- c(cal[!wk][1:8], cal[wk][1:2])
  tally <- tibble(date = samp, commercial_tally = c(4, 6, 5, 7, 4, 6, 5, 3, 8, 10), charter_tally = c(1, 0, 1, 1, 0, 1, 1, 0, 2, 2))
  ints  <- tibble(population = "comm_charter", event_date = rep(samp, each = 2), boat_type_clean = rep(c("Commercial", "Charter"), 10),
                  dungeness_kept = rep(c(40, 60), 10), red_rock_kept = 0)
  unsampled <- cal[!cal %in% samp]                      # 4 days: 2 weekdays, 2 weekend days
  roster <- tibble(season = "2024-25", port = "Westport", vessel = "V",
                   date = c(unsampled[1], unsampled[1], unsampled[3], samp[1], samp[2], samp[2], unsampled[2]),
                   status = c("missed", "missed", "interviewed", "interviewed", "interviewed", "missed", "canceled"), contact = NA, notes = NA)
  Pc <- list(census_start_date = "2025-01-06", census_end_date = "2025-01-19", days_wkend = c("Saturday", "Sunday"),
             crabbing_holiday_dates = as.Date(character()))
  obs <- sum(tally$commercial_tally * 40 + tally$charter_tally * 60)
  dwg <- list(comm_tally = tally, interview = ints, charter_roster = roster)
  rt <- estimate_comm_charter(dwg, modifyList(Pc, list(charter_frame = "tally")))
  rr <- estimate_comm_charter(dwg, Pc)
  # roster: 3 trips on 2 unsampled days (canceled one excluded); samp[1] has tally 1 / roster 1 (agree); samp[2] has tally 0 / roster 2 (+2)
  extra <- (2 + 1) * 60 + 2 * 60
  chk("census roster frame: the default is 'roster'; the total adds the roster-only trips and the tally days where the roster has more",
      identical(rt$charter_frame, "tally") && isTRUE(all.equal(rt$Dungeness_Kept, obs)) &&
        identical(rr$charter_frame, "roster") && isTRUE(all.equal(rr$Dungeness_Kept, obs + extra)) &&
        isTRUE(all.equal(rr$charter_roster_dung, 3 * 60)) && rr$n_roster_only_days == 2 && rr$n_roster_trips == 6)
  # 2026-09-11: observed_dung is the draw floor (the commercial census plus the charter crab
  # actually observed), so it sits BELOW the total exactly by the unobserved charter trips.
  chk("census roster frame: the daily table sums to the total, labels the roster-only days, and observed_dung is the floor",
      isTRUE(all.equal(sum(rr$daily_full$est_dung), rr$Dungeness_Kept)) && sum(grepl("charter roster \\(no tally\\)$", rr$daily_full$source)) == 2 &&
        isTRUE(all.equal(rr$observed_dung, rr$commercial_dung + rr$charter_observed_dung)) && rr$observed_dung < rr$Dungeness_Kept &&
        nrow(rr$roster_reconciliation) == 9 && nrow(rr$daily_est) >= 12)
  rd <- estimate_comm_charter(dwg, modifyList(Pc, list(census_expansion = "day_type")))
  chk("census roster frame + day_type: only the commercial part is expanded; the charter part is exact; the daily table still sums to the total",
      isTRUE(all.equal(rd$Dungeness_Kept, sum(tally$commercial_tally * 40) + 2 * mean((tally$commercial_tally * 40)[!weekdays(tally$date) %in% c("Saturday", "Sunday")]) +
                                            2 * mean((tally$commercial_tally * 40)[weekdays(tally$date) %in% c("Saturday", "Sunday")]) + sum(tally$charter_tally) * 60 + extra)) &&
        isTRUE(all.equal(sum(rd$daily_full$est_dung), rd$Dungeness_Kept)))
  chk("census roster frame: a roster with no trips in the window falls back to the tally frame and says so; bad values refused",
      identical(estimate_comm_charter(list(comm_tally = tally, interview = ints, charter_roster = roster |> mutate(date = date + 365)), Pc)$charter_frame, "tally") &&
        inherits(tryCatch(estimate_comm_charter(dwg, modifyList(Pc, list(charter_frame = "union"))), error = function(e) e), "error"))
  w <- tryCatch({ estimate_comm_charter(list(comm_tally = tally[0, ], interview = ints), Pc); NULL }, warning = function(w) conditionMessage(w))
  chk("census: interviews without any tally row warn that the census frame is missing (2023-24), instead of a silent zero",
      !is.null(w) && grepl("NO vessel tally", w))
  Pw <- Pc; Pw$census_windows <- list("a" = c("2025-01-06", "2025-01-12"), "b" = c("2025-01-13", "2025-01-19"))
  rw <- estimate_comm_charter(dwg, Pw)
  chk("census roster frame: the per-window path sums the roster additions and stacks the reconciliation",
      isTRUE(all.equal(rw$Dungeness_Kept, rr$Dungeness_Kept)) && isTRUE(all.equal(rw$charter_roster_dung, rr$charter_roster_dung)) && nrow(rw$roster_reconciliation) == 9)

  # --- the shipped configuration and the drivers ---
  # 2026-09-12: run_config.R ships ONLY the canonical single 2024-25 window (D28), so the
  # alternative windows -- and the 2023-24 census date that goes with the blocked span --
  # are asserted where they now live, NEW_SEASON_GUIDE section 7.1, not in the config.
  chk("shipped: charter_frame = roster, effort_qc_drop holds the interview-total flag",
      identical(rc$charter_frame, "roster") && identical(rc$effort_qc_drop, "gear_count_from_interviews"))
  # 2026-09-28 (Matt): seasons before 2024-25 are out of scope, so the guide's span example is
  # 2024-26 and the 2023-25 paste block is gone; assert THAT, not the retired block.
  chk("shipped: the alternative windows are paste-ready in the guide (2025-26, and a 2024-26 span; no pre-2024-25 season)",
      { g <- paste(readLines("07_documentation/NEW_SEASON_GUIDE.md", warn = FALSE), collapse = "\n")
        grepl("Paste-ready window blocks", g, fixed = TRUE) &&
        grepl('season_filter     = "2025-26"', g, fixed = TRUE) &&
        grepl('run_tag           = "season-2025-26"', g, fixed = TRUE) &&
        grepl('run_tag           = "two-season-2024-26"', g, fixed = TRUE) &&
        !grepl('season_filter     = c("2023-24", "2024-25")', g, fixed = TRUE) })
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    chk(sprintf("%s: the census row names the charter frame", basename(drv)), any(grepl("charter_frame", d, fixed = TRUE)))
  }
  rd_files <- c(list.files("03_R_functions", pattern = "\\.R$", full.names = TRUE))
  direct <- vapply(rd_files, function(f) any(grepl("readxl::read_excel\\(", readLines(f, warn = FALSE))) && !grepl("read_input_workbook.R", f), logical(1))
  chk("readers: no function in 03_R_functions reads an input workbook with readxl directly (all go through read_input_workbook)", !any(direct))
})

# ---------------------------------------------------------------------------
# 58. The charter EXPANSION and the commercial CENSUS, separated (2026-09-11, Matt's
#     correction: "a complete census without error applies to the commercial boats; the
#     charter vessels are not 100% sampled and need expansion"). The commercial part stays
#     the exact sum over the tally days; the charter part becomes N x (mean catch per
#     interviewed trip) over the charter TRIP frame, stratified by vessel, with the
#     finite-population-corrected variance N^2 (1 - n/N) s^2 / n. Also: the three holidays
#     added to the calendar, and the PE empty-effort-stratum fill options that day-typing
#     interacts with.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/estimate_comm_charter.R")
  # a window of 10 days, 6 with a tally; two charter vessels with different catch rates and
  # a third day-set of trips the roster has but the tally does not, so every branch is live.
  cal  <- seq(as.Date("2025-01-06"), as.Date("2025-01-15"), by = "day")
  samp <- cal[c(1, 2, 3, 6, 7, 8)]
  tally <- tibble(date = samp, commercial_tally = c(4, 6, 5, 3, 8, 10), charter_tally = c(1, 0, 1, 0, 1, 1))
  # charter interviews: vessel A on 3 trips (60, 60, 90 -> mean 70), vessel B on 2 (20, 40 -> 30)
  ci <- tibble(population = "comm_charter", event_date = c(samp[1], samp[3], samp[5], samp[6], samp[6]),
               boat_type_clean = "Charter", boat_name = c("A", "A", "A", "B", "B"),
               dungeness_kept = c(60, 60, 90, 20, 40), red_rock_kept = 0)
  comm <- tibble(population = "comm_charter", event_date = rep(samp, each = 3), boat_type_clean = "Commercial",
                 boat_name = NA_character_, dungeness_kept = rep(c(30, 40, 50), 6), red_rock_kept = 0)
  ints <- bind_rows(comm, ci)
  # the roster: A sails 6 trips (3 interviewed), B sails 3 (2 interviewed), one B trip on a
  # day with no tally; one canceled trip that must not count.
  roster <- tibble(season = "2024-25", port = "Westport", vessel = c(rep("A", 6), rep("B", 3), "A"),
                   date = c(samp[1], samp[2], samp[3], samp[4], samp[5], samp[6], samp[5], samp[6], cal[10], samp[2]),
                   status = c(rep("interviewed", 3), rep("missed", 3), "interviewed", "interviewed", "missed", "canceled"),
                   contact = NA, notes = NA)
  Pc <- list(census_start_date = "2025-01-06", census_end_date = "2025-01-15", days_wkend = c("Saturday", "Sunday"),
             crabbing_holiday_dates = as.Date(character()))
  dwg <- list(comm_tally = tally, interview = ints, charter_roster = roster)
  md_comm <- 40; mA <- 70; mB <- 30; mPool <- mean(ci$dungeness_kept)          # 54
  r <- estimate_comm_charter(dwg, Pc)
  chk("charter expansion: the commercial census is the exact tally-day sum and is reported on its own",
      isTRUE(all.equal(r$commercial_vessels, sum(tally$commercial_tally))) &&
        isTRUE(all.equal(r$commercial_dung, sum(tally$commercial_tally) * md_comm)) &&
        isTRUE(all.equal(r$commercial_dung + r$charter_dung, r$Dungeness_Kept)))
  # the frame: A 6 trips, B 3 trips, plus the tally surplus on days the roster does not cover
  su <- sum(pmax(0, tally$charter_tally - c(1, 1, 1, 0, 2, 2)))                # surplus per tally day
  chk("charter expansion: the trip frame is the roster unioned per day with the tally's charter column",
      isTRUE(all.equal(r$charter_trips, 6 + 3 + su)) && r$n_roster_trips == 9 && r$n_roster_only_days == 1 &&
        isTRUE(all.equal(r$charter_observed_dung, sum(ci$dungeness_kept))))
  chk("charter expansion: stratified by vessel, est = sum_v N_v m_v (+ the unattributed surplus at the pooled mean)",
      identical(r$charter_expansion, "vessel") &&
        isTRUE(all.equal(r$charter_dung, 6 * mA + 3 * mB + su * mPool)))
  # the FPC variance, by hand
  sA <- sd(c(60, 60, 90)); sB <- sd(c(20, 40)); sP <- sd(ci$dungeness_kept)
  vA <- 6^2 * (1 - 3 / 6) * sA^2 / 3; vB <- 3^2 * (1 - 2 / 3) * sB^2 / 2
  vS <- if (su > 0) su^2 * max(0, 1 - 5 / max(su, 5)) * sP^2 / 5 else 0
  chk("charter expansion: the variance is the SRS-with-FPC sum over strata, N^2 (1 - n/N) s^2 / n",
      isTRUE(all.equal(r$charter_var, vA + vB + vS)) && isTRUE(all.equal(r$charter_se, sqrt(vA + vB + vS))) &&
        isTRUE(all.equal(r$charter_sampled_frac, 5 / r$charter_trips)))
  chk("charter expansion: a vessel with one interview borrows the pooled SD and says so",
      { r1 <- estimate_comm_charter(list(comm_tally = tally, interview = bind_rows(comm, ci[1:4, ]), charter_roster = roster), Pc)
        any(grepl("1 interview; pooled SD", r1$charter_detail$mean_source)) && all(is.finite(r1$charter_detail$var)) })
  rp <- estimate_comm_charter(dwg, modifyList(Pc, list(charter_expansion = "pooled")))
  chk("charter expansion: 'pooled' is ONE stratum at the pooled mean (not the vessel strata with a pooled n)",
      identical(rp$charter_expansion, "pooled") && nrow(rp$charter_detail) == 1 &&
        isTRUE(all.equal(rp$charter_dung, rp$charter_trips * mPool)) &&
        isTRUE(all.equal(rp$charter_var, rp$charter_trips^2 * (1 - 5 / rp$charter_trips) * sP^2 / 5)))
  chk("charter expansion: 'tally' frame drops the roster-only trips and says the frame is incomplete",
      { rt <- estimate_comm_charter(dwg, modifyList(Pc, list(charter_frame = "tally")))
        identical(rt$charter_frame, "tally") && isTRUE(all.equal(rt$charter_trips, sum(tally$charter_tally))) &&
          rt$n_roster_only_days == 0L && rt$charter_dung < r$charter_dung })
  chk("charter expansion: observed_dung is the draw floor and sits below the total by the unobserved trips",
      isTRUE(all.equal(r$observed_dung, r$commercial_dung + r$charter_observed_dung)) && r$observed_dung < r$Dungeness_Kept)
  # the daily table reconciles in every frame x expansion combination
  chk("charter expansion: census_daily.csv sums to the component in all four frame x expansion combinations",
      all(vapply(list(list(), list(census_expansion = "day_type"), list(charter_frame = "tally"),
                      list(charter_frame = "tally", census_expansion = "day_type")),
                 function(o) { z <- estimate_comm_charter(dwg, modifyList(Pc, o))
                               isTRUE(all.equal(sum(z$daily_full$est_dung), z$Dungeness_Kept)) &&
                                 isTRUE(all.equal(z$commercial_dung + z$charter_dung, z$Dungeness_Kept)) }, logical(1))))
  chk("charter expansion: carried_var follows census_uncertainty and the drivers use carried_se",
      isTRUE(all.equal(r$carried_var, r$charter_var)) &&
        isTRUE(all.equal(estimate_comm_charter(dwg, modifyList(Pc, list(census_uncertainty = "none")))$carried_var, 0)) &&
        { z <- estimate_comm_charter(dwg, modifyList(Pc, list(census_uncertainty = "sampling"))); isTRUE(all.equal(z$carried_var, z$Dungeness_Kept_var)) })
  # 2026-09-10: BOTH PE runners must delegate to the ONE shared implementation. Three
  # defects lived in two copies of this code; a fix to one copy would have been a fix to
  # half the pipeline, and the pooled and gear PEs would then disagree, which is exactly
  # what run_pe_gear was extracted to prevent.
  for (pf in c("03_R_functions/run_pe_pooled.R", "03_R_functions/run_pe_gear.R")) {
    src <- paste(readLines(pf, warn = FALSE), collapse = "\n")
    chk(sprintf("%s: delegates the strata, the fill and the variance to pe_effort_strata.R", basename(pf)),
        grepl("pe_build_effort_strata(daily_effort, days, params)", src, fixed = TRUE) &&
        grepl("pe_effort_stratum_report(effort_strat, population_name, params)", src, fixed = TRUE) &&
        grepl("pe_empty_cpue_fill(daily_cpue, effort_strat, days, params)", src, fixed = TRUE))
    chk(sprintf("%s: no longer carries its own copy of the fill or the SE formula", basename(pf)),
        !grepl("md_mean_daily", src, fixed = TRUE) &&
        !grepl("replace_na(sd_daily^2,0)/pmax(n_sampled,1)", src, fixed = TRUE))
    chk(sprintf("%s: reports both SEs", basename(pf)),
        grepl("results$effort_se_sampled_only", src, fixed = TRUE))
  }
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("shipped: the PE unsampled-cell levers are the 2026-09-10 set (Matt: ship local_day_type on)",
      identical(rc$pe_empty_effort_stratum, "local_day_type") && identical(rc$pe_empty_stratum, "local") &&
      identical(rc$pe_variance, "impute_aware"))
  chk("shipped: the census keys are the 2026-09-11 set", identical(rc$census_expansion, "none") && identical(rc$charter_frame, "roster") &&
        identical(rc$charter_expansion, "vessel") && identical(rc$census_uncertainty, "charter"))
})

# ---------------------------------------------------------------------------
# 59. THE PE's UNSAMPLED AND SINGLETON CELLS (2026-09-10). Three defects, all of them in
#     two copies of the same code, all recomputed here BY HAND on a fixture small enough
#     to check with a calculator:
#       - the SE was sqrt(N^2 sd^2 / max(n,1)) with sd from the cell's own sampled days,
#         and sd() of ONE observation is NA -> 0, so a singleton cell contributed its full
#         point estimate and no variance;
#       - an IMPUTED cell did the same, so the effort SE was unchanged by imputation;
#       - the effort fill could be month-local while the CPUE fill was sub-season-wide.
#     The fixture is two months x two day types so every donor level is exercised, and the
#     "sampled_only" arm must reproduce the pre-2026-09-10 arithmetic EXACTLY.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/pe_effort_strata.R")
  # calendar: Jan (period 1) and Feb (period 5), weekday + weekend cells.
  mkdays <- function() {
    d <- tibble(event_date = c(as.Date("2025-01-06") + 0:6, as.Date("2025-02-03") + 0:6))
    d |> mutate(day_type = ifelse(format(event_date, "%u") %in% c("6", "7"), "weekend", "weekday"),
                month = as.numeric(format(event_date, "%m")),
                period = as.numeric(format(event_date, "%W")),
                open_section_1 = TRUE)
  }
  days <- mkdays()
  # sampled days: Jan weekday x3 (10, 20, 30), Jan weekend x1 (100), Feb weekday x2 (2, 4).
  # The Feb WEEKEND cell has NO sampled day: that is the imputed cell.
  de <- tibble(
    event_date = c(as.Date("2025-01-06"), as.Date("2025-01-07"), as.Date("2025-01-08"),
                   as.Date("2025-01-11"), as.Date("2025-02-03"), as.Date("2025-02-04")),
    est_daily_effort = c(10, 20, 30, 100, 2, 4), section_num = 1) |>
    left_join(days |> select(event_date, day_type, period), by = "event_date")
  P0 <- list(pe_empty_effort_stratum = "zero",           pe_variance = "sampled_only")
  P1 <- list(pe_empty_effort_stratum = "local_day_type", pe_variance = "sampled_only")
  P2 <- list(pe_empty_effort_stratum = "local_day_type", pe_variance = "impute_aware")
  P3 <- list(pe_empty_effort_stratum = "local_day_type", pe_variance = "donor_mean_only")
  P4 <- list(pe_empty_effort_stratum = "day_type",       pe_variance = "impute_aware")
  s0 <- pe_build_effort_strata(de, days, P0); s1 <- pe_build_effort_strata(de, days, P1)
  s2 <- pe_build_effort_strata(de, days, P2); s3 <- pe_build_effort_strata(de, days, P3)
  s4 <- pe_build_effort_strata(de, days, P4)
  .g <- function(st, per, dt, col) st[[col]][st$period == per & st$day_type == dt]
  jan_wd <- unique(days$period[days$month == 1 & days$day_type == "weekday"])[1]
  feb_wd <- unique(days$period[days$month == 2 & days$day_type == "weekday"])[1]
  feb_we <- unique(days$period[days$month == 2 & days$day_type == "weekend"])[1]
  jan_we <- unique(days$period[days$month == 1 & days$day_type == "weekend"])[1]

  # ---- the cell taxonomy is what the levers key off ------------------------
  k0 <- attr(s0, "counts")
  # 4 cells over 14 days: (Jan, weekday) n = 3, (Jan, weekend) n = 1 -> SINGLETON,
  # (Feb, weekday) n = 2, (Feb, weekend) n = 0 -> UNSAMPLED.
  chk("PE strata: the fixture has one unsampled cell, one singleton cell and two multi cells",
      identical(k0$n_empty_strata, 1L) && identical(k0$n_single_strata, 1L) &&
      identical(k0$n_strata_total, 4L) && identical(k0$n_calendar_days, 14L) &&
      identical(k0$n_empty_days, 2L) && identical(k0$n_single_days, 2L))

  # ---- 1. the MULTI cell: unchanged, and it is the only case the old code got right ----
  # Jan weekday: n = 3, days = 5, mean = 20, sd = 10 -> est 100, se = 5 * 10 / sqrt(3)
  nd <- .g(s2, jan_wd, "weekday", "n_total_days")
  chk("PE strata: a cell with 2+ sampled days is untouched (mean, and se = N sd / sqrt(n))",
      isTRUE(all.equal(.g(s2, jan_wd, "weekday", "mean_daily"), 20)) &&
      isTRUE(all.equal(.g(s2, jan_wd, "weekday", "est_total"), 20 * nd)) &&
      isTRUE(all.equal(.g(s2, jan_wd, "weekday", "se_total"), nd * 10 / sqrt(3))) &&
      isTRUE(all.equal(.g(s2, jan_wd, "weekday", "se_total"), .g(s0, jan_wd, "weekday", "se_total"))))

  # ---- 2. the SINGLETON cell: was zero variance, now the collapsed-stratum donor sd ----
  # Feb weekday has n = 2 (2, 4) so it is NOT a singleton; Jan weekend (100) and Feb
  # weekend (imputed) are. Jan weekend: n = 1, so sd is undefined. Its donor for the SPREAD
  # is the day_type level (weekend has 1 sampled day across the fixture -> falls to the
  # sub-season sd over all 6 sampled days).
  sd_all <- sd(de$est_daily_effort)
  n_jw <- .g(s2, jan_we, "weekend", "n_total_days")
  chk("PE strata: a singleton cell had ZERO variance under the old arithmetic",
      isTRUE(all.equal(.g(s0, jan_we, "weekend", "se_total"), 0)) &&
      isTRUE(all.equal(.g(s1, jan_we, "weekend", "se_total"), 0)))
  chk("PE strata: a singleton cell now borrows its donor's spread with the divisor still 1 (collapsed stratum)",
      isTRUE(all.equal(.g(s2, jan_we, "weekend", "se_total"), n_jw * sd_all)) &&
      grepl("^collapsed stratum", .g(s2, jan_we, "weekend", "var_source")))
  chk("PE strata: the singleton fix is INDEPENDENT of the fill (same SE under day_type and local_day_type)",
      isTRUE(all.equal(.g(s2, jan_we, "weekend", "se_total"), .g(s4, jan_we, "weekend", "se_total"))))

  # ---- 3. the IMPUTED cell: mean from the finest level, spread from the finest with 2+ ----
  # Feb weekend is unsampled. Under local_day_type its MEAN comes from the finest level with
  # any sampled day: (Feb x weekend) has none, so it falls to day_type (weekend = 100).
  # Under day_type it also gets the weekend mean, 100. Under "zero" it gets 0.
  n_fw <- .g(s1, feb_we, "weekend", "n_total_days")
  chk("PE strata: 'zero' leaves the unsampled cell at zero and reports the omission as a BIAS",
      isTRUE(all.equal(.g(s0, feb_we, "weekend", "est_total"), 0)) &&
      isTRUE(all.equal(.g(s0, feb_we, "weekend", "se_total"), 0)) &&
      isTRUE(all.equal(.g(s0, feb_we, "weekend", "zeroed_effort_bias"), 100 * n_fw)) &&
      isTRUE(all.equal(attr(s0, "counts")$zeroed_effort_bias, 100 * n_fw)))
  chk("PE strata: a fill imputes the cell and sets the bias to zero (the error becomes representable)",
      isTRUE(all.equal(.g(s1, feb_we, "weekend", "est_total"), 100 * n_fw)) &&
      isTRUE(all.equal(attr(s1, "counts")$zeroed_effort_bias, 0)) &&
      isTRUE(all.equal(attr(s1, "counts")$imputed_effort, 100 * n_fw)))
  chk("PE strata: an imputed cell was FREE under the old arithmetic (point estimate, no variance)",
      isTRUE(all.equal(.g(s1, feb_we, "weekend", "se_total"), 0)))
  # impute_aware: var = N^2 s^2 (1/n_mean + 1). The SPREAD comes from the sub-season (the
  # finest level with 2+ days, n = 6), but the MEAN is the day_type level's, ONE weekend day
  # (100), so the donor mean's own error divides by 1, not 6. Until 2026-09-28 (B44) it
  # divided by the spread level's n, pricing a one-day mean as a six-day one.
  chk("PE strata: an imputed cell now carries the donor-mean variance PLUS a between-cell term",
      isTRUE(all.equal(.g(s2, feb_we, "weekend", "se_total"), n_fw * sd_all * sqrt(1/1 + 1))) &&
      grepl("^imputed \\(donor mean \\+ between-cell\\)", .g(s2, feb_we, "weekend", "var_source")))
  chk("PE strata: donor_mean_only drops the between-cell term and is the LOWER bound",
      isTRUE(all.equal(.g(s3, feb_we, "weekend", "se_total"), n_fw * sd_all / sqrt(1))) &&
      .g(s3, feb_we, "weekend", "se_total") < .g(s2, feb_we, "weekend", "se_total"))

  # ---- 4. the mean and the spread come from DIFFERENT donor levels, on purpose --------
  # Feb weekday has 2 sampled days, so a Feb-weekday-donated cell would take the month
  # mean AND the month sd. Assert the two selections are reported separately.
  chk("PE strata: mean_level and sd_level are recorded per cell and may differ",
      all(c("mean_level", "sd_level") %in% names(s2)) &&
      any(s2$mean_level != s2$sd_level | s2$n_sampled >= 2L))

  # ---- 5. 'sampled_only' must reproduce the pre-2026-09-10 arithmetic EXACTLY ---------
  chk("PE strata: pe_variance = 'sampled_only' reproduces the old SE for every cell",
      isTRUE(all.equal(s1$se_total, s1$se_total_sampled_only)) &&
      isTRUE(all.equal(s0$se_total, s0$se_total_sampled_only)))
  chk("PE strata: the honest SE is reported ALONGSIDE the old one, never instead of it",
      all(c("se_total", "se_total_sampled_only") %in% names(s2)) &&
      sum(s2$se_total^2) > sum(s2$se_total_sampled_only^2))
  chk("PE strata: the variance lever does NOT move the point estimate",
      isTRUE(all.equal(sum(s1$est_total), sum(s2$est_total))) &&
      isTRUE(all.equal(sum(s2$est_total), sum(s3$est_total))))

  # ---- 6. the CPUE fill, scale-matched -------------------------------------
  # Jan interviews: 10 crab / 5 gear = 2.0; Feb: 1 crab / 5 gear = 0.2; pooled 11/10 = 1.1.
  dc <- tibble(event_date = c(as.Date("2025-01-06"), as.Date("2025-02-03")),
               catch = c(10, 1), hrs = c(5, 5))
  f_loc <- pe_empty_cpue_fill(dc, s2, days, list(pe_empty_stratum = "local"))
  f_pool <- pe_empty_cpue_fill(dc, s2, days, list(pe_empty_stratum = "pooled"))
  f_zero <- pe_empty_cpue_fill(dc, s2, days, list(pe_empty_stratum = "zero"))
  i_fw <- which(s2$period == feb_we & s2$day_type == "weekend")
  chk("PE CPUE fill: 'local' gives a February cell February's ratio-of-sums, not the season's",
      isTRUE(all.equal(f_loc[i_fw], 0.2)) && isTRUE(all.equal(unname(f_pool[i_fw]), 1.1)))
  chk("PE CPUE fill: 'pooled' reproduces the pre-2026-09-10 behaviour and 'zero' the pre-2026-07-13 one",
      length(unique(f_pool)) == 1L && isTRUE(all.equal(unique(as.numeric(f_pool)), 1.1)) &&
      all(f_zero == 0))
  chk("PE CPUE fill: the source string says how many cells got a month rate",
      grepl("month ratio-of-sums", attr(f_loc, "source")) && grepl("sub-season", attr(f_pool, "source")))
  chk("PE CPUE fill: an unknown setting STOPS rather than silently falling back",
      inherits(tryCatch(pe_empty_cpue_fill(dc, s2, days, list(pe_empty_stratum = "monthly")), error = function(e) e), "error") &&
      inherits(tryCatch(pe_build_effort_strata(de, days, list(pe_empty_effort_stratum = "mean")), error = function(e) e), "error") &&
      inherits(tryCatch(pe_build_effort_strata(de, days, list(pe_variance = "honest")), error = function(e) e), "error"))

  # ---- 7. the helper's own defaults are the shipped ones -------------------
  sd_def <- pe_build_effort_strata(de, days, list())
  chk("PE strata: the helper's defaults ARE the shipped settings, so a caller that passes nothing is not on the retired path",
      identical(attr(sd_def, "pe_fill"), "local_day_type") && identical(attr(sd_def, "pe_variance"), "impute_aware") &&
      isTRUE(all.equal(sd_def$se_total, s2$se_total)))

  # ---- 8. estimate_L_effective's dead argument, and the doc claim it supported ---------
  src <- paste(readLines("03_R_functions/bss_day_length.R", warn = FALSE), collapse = "\n")
  chk("L_effective: the dead pot_open_date argument is gone from the signature",
      grepl("estimate_L_effective <- function(ie_data, params)", src, fixed = TRUE) &&
      !grepl("estimate_L_effective <- function(ie_data, pot_open_date, params)", src, fixed = TRUE))
  # this harness file is excluded: it QUOTES the retired signature in the assertion above
  for (cf in setdiff(c(list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE),
                       list.files("06_diagnostics", pattern = "\\.R$", full.names = TRUE)),
                     "06_diagnostics/test_improvements_2026-08-25.R")) {
    cs <- paste(readLines(cf, warn = FALSE), collapse = "\n")
    if (!grepl("estimate_L_effective(", cs, fixed = TRUE)) next
    chk(sprintf("L_effective: %s calls it with the two-argument signature", basename(cf)),
        !grepl("estimate_L_effective(ie_data, params$pot_open_date", cs, fixed = TRUE) &&
        !grepl("estimate_L_effective(ie, p$pot_open_date", cs, fixed = TRUE))
  }
  rcs <- paste(readLines("run_config.R", warn = FALSE), collapse = "\n")
  chk("L_effective: run_config no longer claims pot_open_date feeds an I/E regression split",
      !grepl("feeding the\n  # L_effective I/E regression split", rcs) &&
      !grepl("SINGLE date feeding the", rcs))

  # ---- 9. the per-season calendar guard ------------------------------------
  source("03_R_functions/validate_season_window.R")
  eff <- tibble(date = as.Date("2024-12-01") + 0:9, season = "2024-25")
  int <- tibble(event_date = as.Date("2024-12-01") + 0:9, season = "2024-25")
  Pw <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15", season_filter = "2024-25",
             pot_closures = NULL, pot_closure_start = "2024-09-16", pot_closure_end = "2024-11-30",
             pot_open_date = "2024-12-01", census_start_date = "2024-12-01", census_end_date = "2025-02-08")
  chk("season window: a complete per-season pin warns about nothing",
      length(suppressMessages(withCallingHandlers(
        { w <- character(0); validate_season_window(eff, int, Pw, quiet = TRUE); w },
        warning = function(cd) { w <<- c(w, conditionMessage(cd)); invokeRestart("muffleWarning") }))) == 0)
  chk("season window: a stale pot_open_date from a multi-season rollback now WARNS",
      { w <- character(0)
        withCallingHandlers(validate_season_window(eff, int, modifyList(Pw, list(pot_open_date = "2023-12-01")), quiet = TRUE),
                            warning = function(cd) { w <<- c(w, conditionMessage(cd)); invokeRestart("muffleWarning") })
        any(grepl("pot_open_date", w)) && any(grepl("PER-SEASON", w)) })
  chk("season window: the guard is skipped when pot_closures carries the multi-season calendar",
      { w <- character(0)
        withCallingHandlers(validate_season_window(eff, int, modifyList(Pw, list(
          pot_closures = list(list(season = "2024-25", start = "2024-09-16", end = "2024-11-30")),
          pot_open_date = "2023-12-01")), quiet = TRUE),
          warning = function(cd) { w <<- c(w, conditionMessage(cd)); invokeRestart("muffleWarning") })
        !any(grepl("pot_open_date", w)) })
})

# ---------------------------------------------------------------------------
# 60. THE TWO-PASS LADDER RUN (2026-09-10). Matt: "run the 4-rung version now and then
#     follow up with the new R2f control rung." Three things have to hold for that to be
#     safe, and none of them held when the plan was proposed:
#       (a) adding R2f must not change any other rung's config digest, or pass 2 re-fits
#           12 h of MCMC it was supposed to reuse;
#       (b) the f/c exclusion list used by the factorization verdict must cover every
#           f/c quantity the Stan model REPORTS -- sigma_c_out and cfc_kappa_out are
#           declared unconditionally and set to exactly 0.0 when the walk is off, so a
#           missing name means a z in the tens and a spurious FAIL after four hours;
#       (c) a code change between the two passes must be detected, because the digest
#           covers the configuration and R2f is only interpretable against an R2 fitted
#           by the same code.
# ---------------------------------------------------------------------------
local({
  f <- "06_diagnostics/run_improvements_2026-09-08.R"
  chk("two-pass: ladder runner present", file.exists(f)); if (!file.exists(f)) return(invisible(NULL))
  t <- readLines(f, warn = FALSE); tt <- t[!grepl("^\\s*#", t)]; s <- paste(tt, collapse = "\n")

  # ---- the control block -------------------------------------------------
  chk("two-pass: LADDER_PASS exists and ships pass 1 (the four citable rungs first)",
      any(grepl("^LADDER_PASS <- 1", t)) && grepl("LADDER_PASS >= 2", s, fixed = TRUE))
  # 2026-09-11: R4 is fitted FIRST. It is the shipped configuration, so the citable number
  # lands at hour 4 of a 16-hour run instead of hour 11. Verified safe before the run and
  # confirmed by it: the verdict blocks execute after the whole loop and resolve their
  # comparison folders by rung NAME, so fit order cannot reach them.
  chk("two-pass: R4 is fitted first, and pass 2 adds R2f after R2",
      grepl('c("R0", "R4", "R1", "R2", "R5")', s, fixed = TRUE) &&
      grepl('c("R0", "R4", "R1", "R2", "R2f", "R5")', s, fixed = TRUE))
  chk("two-pass: the run prints what to do next, so the plan is not carried in someone's head",
      grepl("PASS 1 COMPLETE", s, fixed = TRUE) && grepl("PASS 2 COMPLETE", s, fixed = TRUE) &&
      grepl("LADDER_PASS <- 2", s, fixed = TRUE))

  # ---- (a) the digests are stable across the two passes ------------------
  # Evaluate the runner's head twice, once per STAGES setting, and compare.
  .head <- t[seq_len(grep("^# PRE-FLIGHT", t)[1] - 1)]
  .dig <- function(pass) {
    h <- sub("^LADDER_PASS <- [0-9]+", paste0("LADDER_PASS <- ", pass), .head)
    e <- new.env(parent = globalenv())
    eval(parse(text = paste(h, collapse = "\n")), envir = e)
    st <- setdiff(get("STAGES", envir = e), "R0")
    setNames(vapply(st, get("stage_digest", envir = e), character(1)), st)
  }
  d1 <- tryCatch(.dig(1), error = function(e) NULL)
  d2 <- tryCatch(.dig(2), error = function(e) NULL)
  chk("two-pass: both passes resolve their rungs without error", !is.null(d1) && !is.null(d2))
  if (!is.null(d1) && !is.null(d2)) {
    sh <- intersect(names(d1), names(d2))
    chk("two-pass: pass 2 adds R2f and nothing else",
        setequal(setdiff(names(d2), names(d1)), "R2f") && length(setdiff(names(d1), names(d2))) == 0)
    chk("two-pass: EVERY shared rung's config digest is identical across the passes (so pass 2 refits only R2f)",
        length(sh) == 4L && all(d1[sh] == d2[sh]))
    chk("two-pass: R2f's digest is distinct from R2's", !identical(d2[["R2f"]], d2[["R2"]]))
    chk("two-pass: the digest names its stage and the f mode, so a folder cannot be misread",
        all(grepl("^R[0-9a-z]+\\|new_throughout\\|[0-9a-f]{8}$", d2)))
  }

  # ---- (b) the f/c exclusion list covers what Stan reports ---------------
  .i <- grep("^F_EXCLUDE <- paste0", t)[1]
  chk("two-pass: the f/c exclusion list is defined ONCE and reused", is.finite(.i) &&
      sum(grepl("exclude = F_EXCLUDE", tt, fixed = TRUE)) >= 2 &&
      !any(grepl('exclude = "\\^\\(f_crab', tt)))
  if (is.finite(.i)) {
    .j <- .i; while (!grepl("\\)\\s*$", t[.j])) .j <- .j + 1L
    FEX <- eval(parse(text = paste(t[.i:.j], collapse = "\n")))
    for (mf in c("crab_bss_pooled.stan", "crab_bss_gear_resolved.stan")) {
      st <- readLines(file.path("02_stan_models", mf), warn = FALSE)
      i_gq <- grep("^generated quantities", st)[1]; i_pr <- grep("^parameters", st)[1]
      i_tp <- grep("^transformed parameters", st)[1]
      decl <- c(if (is.finite(i_gq)) st[seq(i_gq, length(st))],
                if (is.finite(i_pr) && is.finite(i_tp)) st[seq(i_pr, i_tp)])
      nm <- unique(unlist(regmatches(decl, gregexpr("(?<=\\s)[A-Za-z_][A-Za-z0-9_]*(?=\\s*[;=])", decl, perl = TRUE))))
      # the f and c quantities, excluding the integer sizes (n_f_dyn etc, never reported)
      fc <- sort(nm[grepl("^(f_|z_f|z_c|sigma_f|sigma_c|cfi|cfc|combo|osp_f|eta_f)", nm)])
      fc <- fc[!grepl("^n_", fc)]
      miss <- fc[!grepl(FEX, fc)]
      chk(sprintf("two-pass: F_EXCLUDE covers every f/c quantity %s reports (%d checked)", mf, length(fc)),
          length(fc) > 0 && !length(miss))
      if (length(miss)) cat("      NOT COVERED:", paste(miss, collapse = ", "), "\n")
    }
    # THE MECHANISM THE ANALYSIS RESTS ON, pinned. A 2026-09-10 smoke fit of the boat
    # all-gear component under the R2f configuration showed the switched-off f/c outputs
    # reporting sd = 0 and therefore se_mean = NaN, NOT 0 -- and fit_agreement() skips any
    # row whose combined se is not finite, which is why the missing names above are
    # currently MASKED rather than fatal. That masking is an rstan::summary()
    # implementation detail, so pin it: if it ever changes, this says so instead of a
    # four-hour verdict saying it.
    source("03_R_functions/batch_verdict_helpers.R")
    .mk2 <- function(dir, means, ses) {
      dir.create(dir, showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(data.frame(mean = means, se_mean = ses, sd = ses * 10,
                                  row.names = c("B1", "zz_nan_side", "zz_zero_side", "zz_finite")),
                       file.path(dir, "bss_full_summary_private_boat_all_gear_Dungeness_Kept.csv"))
    }
    .da <- tempfile("fx_a"); .db <- tempfile("fx_b")
    # ref: real posteriors; new: NaN on one row (the R2f case), 0 on another, finite on a third
    .mk2(.da, c(1.0, 0.30, 0.30, 1.00), c(0.01, 0.02, 0.02, 0.01))
    .mk2(.db, c(1.0, 0.00, 0.00, 2.00), c(0.01, NaN,  0.00, 0.01))
    .fa <- fit_agreement(.db, .da, pat = "private_boat", exclude = NULL, what = "fixture")
    chk("two-pass: fit_agreement SKIPS a row whose se_mean is NaN on one side (why the missing names were masked, not fatal)",
        grepl("zz_finite", .fa$observed) && !grepl("zz_nan_side", .fa$observed))
    chk("two-pass: fit_agreement COMPARES a row whose se_mean is 0 on one side, so a future rstan reporting 0 would bite",
        { z <- .fa$z; any(grepl("zz_zero_side$", names(z))) })
    chk("two-pass: fit_agreement still catches a genuine disagreement in the fixture",
        identical(.fa$verdict, "FAIL") && isTRUE(max(.fa$z, na.rm = TRUE) > 5))
    unlink(c(.da, .db), recursive = TRUE)

    # the specific two that were missing, named so a regression is unmistakable
    chk("two-pass: sigma_c_out and cfc_kappa_out are excluded (declared unconditionally, 0.0 when the c walk is off)",
        grepl(FEX, "sigma_c_out") && grepl(FEX, "cfc_kappa_out") &&
        grepl(FEX, "sigma_f_out") && grepl(FEX, "cfi_kappa_out") && grepl(FEX, "combo_c_out"))
    chk("two-pass: F_EXCLUDE does NOT swallow the parameters the factorization proof must actually compare",
        !grepl(FEX, "B1") && !grepl(FEX, "B2") && !grepl(FEX, "tau_bar_out") && !grepl(FEX, "R_G_boat_out") &&
        !grepl(FEX, "mu_mu_E[1]") && !grepl(FEX, "sigma_eps_C") && !grepl(FEX, "sigma_IE_out") &&
        !grepl(FEX, "sigma_mu_C") && !grepl(FEX, "sigma_r_C"))
  }

  # ---- a failed rung must not destroy the others -------------------------
  chk("two-pass: a rung that errors is recorded and the remaining rungs still run",
      grepl("dirs[[sid]] <- tryCatch(run_stage(sid), error = function(e)", s, fixed = TRUE) &&
      grepl("RUNG %s FAILED", s, fixed = TRUE) &&
      grepl('V1row(sid, "the rung did not complete"', s, fixed = TRUE) &&
      !grepl("for (sid in STAGES) dirs[[sid]] <- run_stage(sid)", s, fixed = TRUE))
  chk("two-pass: the run reports its elapsed time per fitted rung",
      grepl("elapsed %.1f h of the ladder so far", s, fixed = TRUE))

  # ---- (c) the code fingerprint ------------------------------------------
  chk("two-pass: IMP_STAGE.txt records a code fingerprint over the Stan models, the drivers and 03_R_functions",
      grepl("code_fingerprint <- function", s, fixed = TRUE) &&
      grepl('sprintf("code: %s", code_fingerprint())', s, fixed = TRUE) &&
      grepl('.code_group("02_stan_models"', s, fixed = TRUE) &&
      grepl('.code_group("01_BSS_models"', s, fixed = TRUE) &&
      grepl('.code_group("03_R_functions"', s, fixed = TRUE))
  chk("two-pass: a reused folder whose code has changed is REPORTED, not silently trusted",
      grepl("THE CODE HAS CHANGED SINCE THAT FIT", s, fixed = TRUE) &&
      grepl("reused fit was produced by DIFFERENT code than this run", s, fixed = TRUE))
  chk("two-pass: code drift does NOT force a refit (a comment change must not cost 12 h)",
      grepl("return(existing)", s, fixed = TRUE))
  chk("two-pass: every cross-rung claim goes through V1cross, which downgrades PASS when the premise is broken",
      grepl("V1cross <- function", s, fixed = TRUE) &&
      grepl('if (nzchar(note) && identical(verdict, "PASS")) verdict <- "REVIEW"', s, fixed = TRUE) &&
      sum(grepl("V1cross(", tt, fixed = TRUE)) >= 7)
  chk("two-pass: the bit-identity and agreement verdicts are the ones wrapped",
      grepl('V1cross("R2", "the shore did not move', s, fixed = TRUE) &&
      grepl('V1cross("R2f", "the f block leaves effort and CPUE untouched', s, fixed = TRUE) &&
      grepl('V1cross("R4", "the boat did not move', s, fixed = TRUE))
  chk("two-pass: the comparability note also checks rstan/StanHeaders, which a code hash cannot see",
      grepl(".stan_versions_of <- function", s, fixed = TRUE) &&
      grepl("rstan/StanHeaders differ", s, fixed = TRUE))

  # the fingerprint itself: deterministic, and it changes when a watched file changes
  .head2 <- t[seq_len(grep("^# PRE-FLIGHT", t)[1] - 1)]
  e2 <- new.env(parent = globalenv())
  ev <- tryCatch({ eval(parse(text = paste(.head2, collapse = "\n")), envir = e2); TRUE }, error = function(e) FALSE)
  chk("two-pass: the runner head evaluates so the fingerprint can be exercised", ev)
  if (ev) {
    cf <- get("code_fingerprint", envir = e2)
    a <- cf(); b <- cf()
    chk("two-pass: the code fingerprint is deterministic and names all three layers",
        identical(a, b) && grepl("^stan:[0-9a-f]{8} drivers:[0-9a-f]{8} fns:[0-9a-f]{8}$", a))
    cd <- get(".code_delta", envir = e2)
    chk("two-pass: .code_delta reports NO difference against itself and names the layer that moved",
        identical(cd(a), "") &&
        identical(cd(sub("^stan:[0-9a-f]{8}", "stan:deadbeef", a)), "stan") &&
        identical(cd(sub("fns:[0-9a-f]{8}$", "fns:deadbeef", a)), "fns") &&
        is.na(cd(NA_character_)))
    # 2026-09-11: comments and blank lines are stripped before hashing, so a doc-only
    # patch cannot downgrade an existing record's verdicts, and an EQUIVALENCE list
    # carries the ones that still move the hash but provably not the inference.
    ceq <- get("CODE_EQUIVALENT", envir = e2)
    chk("two-pass: the fingerprint ignores comments and blank lines",
        { td <- tempfile(); dir.create(file.path(td, "02_stan_models"), recursive = TRUE)
          dir.create(file.path(td, "01_BSS_models")); dir.create(file.path(td, "03_R_functions"))
          writeLines(c("x <- 1", "y <- 2"), file.path(td, "03_R_functions", "a.R"))
          h1 <- local({ .here <- function(...) file.path(td, ...); get(".code_group", envir = e2) })
          f1 <- environment(); g <- get(".code_group", envir = e2)
          e3 <- new.env(parent = environment(g)); e3$.here <- function(...) file.path(td, ...)
          environment(g) <- e3
          v1 <- g("03_R_functions", "\\.R$")
          writeLines(c("# a new comment", "x <- 1", "", "y <- 2   # trailing"), file.path(td, "03_R_functions", "a.R"))
          v2 <- g("03_R_functions", "\\.R$")
          writeLines(c("x <- 1", "y <- 3"), file.path(td, "03_R_functions", "a.R"))
          v3 <- g("03_R_functions", "\\.R$")
          unlink(td, recursive = TRUE)
          identical(v1, v2) && !identical(v1, v3) })
    .fp <- "stan:[0-9a-f]{8} drivers:[0-9a-f]{8} fns:[0-9a-f]{8}"
    chk("two-pass: the equivalence list is an audit trail -- every entry carries a reason",
        is.list(ceq) && length(ceq) >= 1 &&
        all(grepl(sprintf("^%s => %s$", .fp, .fp), names(ceq))) &&
        all(nzchar(unlist(ceq))) && all(nchar(unlist(ceq)) > 40))
    # 2026-09-12: A DECLARATION WHOSE "=>" SIDE IS NO LONGER CURRENT IS A STALE DECLARATION,
    # and the previous one went stale the moment an executable change landed with nothing to
    # catch it. So: at least one entry must name the ACTUAL current fingerprint, and any entry
    # that does not must say SUPERSEDED in its reason. That is cheap to satisfy (re-make the
    # claim, or drop it) and it makes the list re-examined rather than inherited.
    # 2026-09-28 (B44): OR the tree was examined and declared NOT equivalent, with the reason,
    # in CODE_NOT_EQUIVALENT. A change that does alter a fit must not be written into the
    # equivalence list to satisfy this check; the second list is where it is recorded.
    cne <- if (exists("CODE_NOT_EQUIVALENT", envir = e2, inherits = FALSE)) get("CODE_NOT_EQUIVALENT", envir = e2) else list()
    chk("two-pass: the ACTUAL current fingerprint is declared, equivalent or examined-and-NOT-equivalent (with a reason)",
        any(vapply(names(ceq), function(k)
              identical(strsplit(k, " => ", fixed = TRUE)[[1]][2], a), logical(1))) ||
        any(vapply(names(cne), function(k)
              identical(strsplit(k, " => ", fixed = TRUE)[[1]][2], a) &&
              grepl("NOT inference-equivalent", cne[[k]], fixed = TRUE) && nchar(cne[[k]]) > 200, logical(1))),
        sprintf("(current is %s)", a))
    chk("two-pass: every declaration that does NOT name the current fingerprint says SUPERSEDED",
        all(vapply(names(ceq), function(k) {
              cur <- identical(strsplit(k, " => ", fixed = TRUE)[[1]][2], a)
              cur || grepl("SUPERSEDED", ceq[[k]], fixed = TRUE) }, logical(1))))
    # 2026-09-12: the key names BOTH ends. Keyed on the recorded side alone the entry never
    # expired -- it excused its folder against whatever the tree later became -- so a Stan
    # edit would have passed unflagged on all five committed rungs. These three assert the
    # declaration applies to the pair it names and LAPSES on anything else.
    chk("two-pass: a declared PAIR reports no delta; a fingerprint off the list still does",
        { kp <- strsplit(names(ceq)[1], " => ", fixed = TRUE)[[1]]
          identical(cd(kp[1], kp[2]), "") &&
          identical(cd(sub("fns:[0-9a-f]{8}$", "fns:deadbeef", a)), "fns") })
    chk("two-pass: THE DEFECT -- an equivalence declaration must LAPSE when the other end moves",
        { kp <- strsplit(names(ceq)[1], " => ", fixed = TRUE)[[1]]
          d <- cd(kp[1], sub("^stan:[0-9a-f]{8}", "stan:deadbeef", kp[2]))
          nzchar(d) && grepl("stan", d, fixed = TRUE) })
    # pre-B24 stamps hashed RAW text, so the current function can never emit one; each is
    # mapped to code_fingerprint() of the SAME tree before any layer comparison.
    cleg <- get("CODE_LEGACY", envir = e2); cnorm <- get(".code_norm", envir = e2)
    chk("two-pass: every pre-B24 stamp maps to a well-formed current-format fingerprint",
        is.list(cleg) && length(cleg) >= 1 &&
        all(grepl(sprintf("^%s$", .fp), names(cleg))) &&
        all(grepl(sprintf("^%s$", .fp), unlist(cleg))) &&
        !any(names(cleg) %in% unlist(cleg)) &&
        identical(cnorm(names(cleg)[1]), cleg[[1]]) &&
        identical(cnorm(a), a) && is.na(cnorm(NA_character_)))
    # and the stamps actually on disk resolve: a rung folder whose code line is neither
    # current-format-comparable nor mapped is a silently incomparable record.
    chk("two-pass: every committed rung stamp is either current-format or mapped by CODE_LEGACY",
        { st <- Sys.glob(file.path("05_output", "*", "*", "IMP_STAGE.txt"))
          if (!length(st)) TRUE else {
            got <- unique(na.omit(vapply(st, function(p) {
              l <- grep("^code: ", readLines(p, warn = FALSE), value = TRUE)
              if (!length(l)) NA_character_ else sub("^code: ", "", l[1]) }, character(1))))
            length(got) >= 1 && all(vapply(got, function(g)
              !is.null(cleg[[g]]) || identical(cnorm(g), g), logical(1))) } })
  }
})

# ---------------------------------------------------------------------------
# 61. WHAT THE FIRST FULL LADDER RUN EXPOSED (2026-09-11). Three reporting defects and
#     one substantive one, all of them found by reading the run rather than the code:
#       D24  estimate_shore_turnover() never filtered the I/E days to the window, so the
#            second-largest mover in the series is a multi-season quantity. Now reported
#            either way, and pricable with tau_shore_derive_window_only.
#       D25  the R2 boat-rise threshold assumed the turnover passes through 1:1. It does
#            not (14.0% prior move -> 4.4% component move), so a correct rung read REVIEW.
#       D26  verdict_R4's expect_delta omitted the census keys D_R4 adds by design, and it
#            hard-coded "R3" as the predecessor's name.
# ---------------------------------------------------------------------------
local({
  f <- "06_diagnostics/run_improvements_2026-09-08.R"
  chk("post-run: ladder runner present", file.exists(f)); if (!file.exists(f)) return(invisible(NULL))
  t <- readLines(f, warn = FALSE); tt <- t[!grepl("^\\s*#", t)]; s <- paste(tt, collapse = "\n")

  # ---- D26: the R4 exactness row -----------------------------------------
  chk("post-run: verdict_R4 expects the census keys D_R4 adds, so a correct config cannot print UNEXPECTED",
      grepl('"census_uncertainty", "charter_frame", "charter_expansion"', s, fixed = TRUE) &&
      grepl("expect_delta = c(\"tau_shore_prior_mu\"", s, fixed = TRUE))
  chk("post-run: every key D_R4 adds on top of D_R2 is in verdict_R4's expect_delta",
      { e <- new.env(parent = globalenv()); e$F_METHOD <- "new_throughout"; e$modifyList <- modifyList
        e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
        .i1 <- grep("^F_NEW *<-", tt)[1]; .i2 <- grep("^D_R4 *<-", tt)[1]
        ok2 <- FALSE
        for (.end in .i2:min(.i2 + 12L, length(tt))) {
          .txt <- paste(tt[.i1:.end], collapse = "\n")
          if (!inherits(try(parse(text = .txt), silent = TRUE), "try-error")) { eval(parse(text = .txt), envir = e); ok2 <- TRUE; break } }
        if (!ok2) FALSE else {
          moved <- names(e$D_R4)[!vapply(names(e$D_R4), function(k) identical(e$D_R4[[k]], e$D_R2[[k]]), logical(1))]
          .j <- grep("^verdict_R4 <- function", t)[1]; .k <- .j
          while (!grepl("^\\}", t[.k]) || .k == .j) .k <- .k + 1L
          blk <- paste(t[.j:.k], collapse = "\n")
          all(vapply(moved, function(k) grepl(sprintf('"%s"', k), blk, fixed = TRUE), logical(1))) } })
  chk("post-run: verdict_R4 no longer hard-codes 'R3' as the predecessor's name",
      grepl("verdict_R4 <- function(dir, prev, prev_id", s, fixed = TRUE) &&
      grepl("verdict_R4(dirs$R4, .dir_of(prev_pooled(\"R4\")), prev_pooled(\"R4\"))", s, fixed = TRUE) &&
      !grepl('what = "boat fits vs R3"', s, fixed = TRUE) &&
      !grepl('sprintf("sigma_IE %s (R3 %s)"', s, fixed = TRUE))

  # ---- D25: the R2 pass-through ------------------------------------------
  chk("post-run: the R2 boat band is widened and the pass-through ratio is reported",
      grepl("pass-through %s of the prior move", s, fixed = TRUE) &&
      grepl(".pct(ba, ba0) > 2", s, fixed = TRUE) && !grepl(".pct(ba, ba0) > 8", s, fixed = TRUE))

  # ---- D24: the sigma_IE row has to say how many observations it rests on -
  chk("post-run: the sigma_IE row reports the in-window I/E day count",
      grepl("in-window I/E day(s) of %s in the workbook", s, fixed = TRUE) &&
      grepl("L_effective_ie_detail.csv", s, fixed = TRUE))

  # ---- D24: the turnover derivation ---------------------------------------
  bd <- paste(readLines("03_R_functions/bss_day_length.R", warn = FALSE), collapse = "\n")
  chk("post-run: estimate_shore_turnover reports the in-window share and offers a window-only derivation",
      grepl("tau_shore_derive_window_only", bd, fixed = TRUE) &&
      grepl("n_days_in_window = nrow(.bw)", bd, fixed = TRUE) &&
      grepl("tau_in_window = tau_win", bd, fixed = TRUE))
  chk("post-run: shore_turnover_summary.csv carries the window columns",
      grepl("n_days_in_window = st$n_days_in_window", bd, fixed = TRUE) &&
      grepl("tau_in_window = st$tau_in_window", bd, fixed = TRUE) &&
      grepl("derived_window_only = isTRUE(st$derived_window_only)", bd, fixed = TRUE))
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("post-run: the default POOLS every I/E day (the run's numbers stay reproducible)",
      identical(rc$tau_shore_derive_window_only, FALSE))
  chk("post-run: NEW_SEASON_GUIDE no longer claims the derived shore centre comes from the window",
      { g <- paste(readLines("07_documentation/NEW_SEASON_GUIDE.md", warn = FALSE), collapse = "\n")
        !grepl('"derived"` from the window\'s I/E time column', g, fixed = TRUE) &&
        grepl("tau_shore_derive_window_only", g, fixed = TRUE) })

  # the derivation, exercised on a fixture: the window filter must actually bite
  source("03_R_functions/bss_day_length.R")
  set.seed(7)
  mk <- function(dates) do.call(rbind, lapply(dates, function(d) data.frame(
    event_date = as.Date(d), season = "fixture", day_type = "Weekday", hour = 8:16 + 0.5,
    crabbers_on = c(2, 5, 8, 6, 4, 3, 2, 1, 0), crabber_flow = c(2, 7, 14, 18, 19, 18, 15, 11, 8))))
  iv <- mk(c("2023-08-10", "2024-01-10", "2025-01-10", "2025-08-10", "2026-05-10"))
  se <- data.frame(count_sequence = 1L, count_hour = c(10.5, 12.5))
  P <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15", bss_max_count_seq = 3)
  a <- estimate_shore_turnover(iv, se, P, n_boot = 50, quiet = TRUE)
  b <- estimate_shore_turnover(iv, se, modifyList(P, list(tau_shore_derive_window_only = TRUE)), n_boot = 50, quiet = TRUE)
  chk("post-run: the pooled derivation uses every day and the window-only one uses the in-window subset",
      identical(a$n_days, 5L) && identical(b$n_days, 2L) &&
      identical(a$n_days_in_window, 2L) && identical(a$n_days_total, 5L) &&
      isFALSE(a$derived_window_only) && isTRUE(b$derived_window_only))
  chk("post-run: tau_in_window is reported by BOTH paths and equals the window-only tau",
      is.finite(a$tau_in_window) && isTRUE(all.equal(a$tau_in_window, b$tau)) &&
      isTRUE(all.equal(b$tau_in_window, b$tau)))
  chk("post-run: with no window set, the derivation is unchanged and everything counts as in-window",
      { z <- estimate_shore_turnover(iv, se, list(bss_max_count_seq = 3), n_boot = 50, quiet = TRUE)
        identical(z$n_days, 5L) && identical(z$n_days_in_window, 5L) && isTRUE(all.equal(z$tau, a$tau)) })
})

# ---------------------------------------------------------------------------
# 62. THE POINTERS RESOLVE (2026-09-12). Four documents defer to "the box at the top of
#     PIPELINE_STATUS.md" as the single authoritative total. That is a good design and it
#     failed anyway: the ladder ran on 2026-09-11, Section 1v recorded 94,376, and the box
#     four screens above it still read 72,027 -- so every correctly-written pointer in the
#     repository resolved to a superseded number for a day. These assert the property that
#     design depends on: the box agrees with the run it names, and nothing still hands the
#     reader the old total as current.
# ---------------------------------------------------------------------------
local({
  ps <- "07_documentation/development_notes/PIPELINE_STATUS.md"
  chk("pointers: the status document exists", file.exists(ps)); if (!file.exists(ps)) return(invisible(NULL))
  L <- readLines(ps, warn = FALSE)
  box <- paste(L[seq_len(min(60L, length(L)))], collapse = "\n")

  chk("pointers: the box names one authoritative run folder and one port total",
      grepl("THE AUTHORITATIVE RUN", box, fixed = TRUE) &&
      length(regmatches(box, gregexpr("05_output/[0-9]{8}/[A-Za-z0-9._-]+", box))[[1]]) >= 1)
  # position-independent: find the Repo line wherever it sits in the header
  .repo <- grep("^\\*\\*Repo:\\*\\*", L, value = TRUE)
  # 2026-09-28: the branch merged into main (PR #5, a878a87). Until then this asserted the
  # line named the branch and not main; now it must name main AND say which state that is,
  # so the pre-merge main (724eead) cannot be mistaken for the current one.
  chk("pointers: the status document names main as merged from the OSP branch, and the pre-merge state",
      length(.repo) >= 1 &&
      grepl("branch `main`", .repo[1], fixed = TRUE) && grepl("a878a87", .repo[1], fixed = TRUE) &&
      grepl("OSP-boat-count-incorporation", .repo[1], fixed = TRUE) && grepl("724eead", .repo[1], fixed = TRUE))
  chk("pointers: the status document names Method v2.0 as the method of record, and where it is specified",
      length(.repo) >= 1 && grepl("Method v2.0", .repo[1], fixed = TRUE) &&
      grepl("BSS-GH-pooled-CPUE-model-documentation.md", .repo[1], fixed = TRUE))

  # the box against the run it names: the total in the box must be the total in the folder.
  bx <- regmatches(box, regexpr("05_output/[0-9]{8}/[A-Za-z0-9._-]+", box))
  pt <- file.path(bx, "port_total_Dungeness_Kept.csv")
  chk("pointers: the authoritative folder named in the box exists on disk", dir.exists(bx),
      sprintf("(box names %s)", bx))
  chk("pointers: THE BOX AGREES WITH THE RUN IT NAMES (median and both interval ends)",
      { if (!file.exists(pt)) NA else {
          d <- utils::read.csv(pt, stringsAsFactors = FALSE)
          r <- d[grepl("^Expected", d[[2]]), , drop = FALSE]
          if (nrow(r) != 1) FALSE else {
            f <- function(x) formatC(round(as.numeric(x)), format = "d", big.mark = ",")
            all(vapply(c(r$BSS_median, r$BSS_lo95, r$BSS_hi95),
                       function(v) grepl(f(v), box, fixed = TRUE), logical(1))) } } },
      "(port_total_Dungeness_Kept.csv vs the box)")

  # nothing hands the reader the superseded total as current. A document may still CONTAIN
  # it (the histories and the reviews are records), but only labelled as superseded.
  live <- c("07_documentation/development_notes/PIPELINE_STATUS.md",
            "07_documentation/development_notes/CHANGE_REGISTER.md",
            "07_documentation/BSS-GH-pooled-CPUE-model-documentation.md",
            "07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md",
            "07_documentation/development_notes/adoption-review-2026-09-08.md",
            "07_documentation/archive/PR-5-merge-OSP-boat-count-incorporation.md", "README.md", "07_documentation/CLAUDE.md",
            # 2026-09-12: the per-folder READMEs were OUTSIDE this list, which is how
            # 06_diagnostics/README.md carried "is the authoritative run (port 72,027)" in
            # the present tense through two supersessions.
            "01_BSS_models/README.md", "02_stan_models/README.md", "04_input_files/README.md",
            "05_output/README.md", "06_diagnostics/README.md", "README-R-functions.md",
            "07_documentation/NEW_SEASON_GUIDE.md", "07_documentation/README.md",
            "07_documentation/BSS-GH-pooled-CPUE-model-documentation.md",
            "07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md")
  bad_ptr <- Filter(function(f) {
    if (!file.exists(f)) return(FALSE)
    ln <- grep("72,027", readLines(f, warn = FALSE), value = TRUE)
    any(!grepl("supersed|SUPERSED|historical|HISTORICAL|previous|used to", ln))
  }, live)
  chk("pointers: no live document offers 72,027 without marking it superseded",
      length(bad_ptr) == 0, sprintf("(%s)", paste(bad_ptr, collapse = ", ")))
  # D37 (2026-09-29): the guard above covered ONE superseded total. Every later one (R4's
  # 94,376, the 2026-09-27 render's 96,118, and their folders) appears legitimately in dozens
  # of historical sentences, so a blacklist would be brittle; what must never happen is a
  # line that CLAIMS a run is the authoritative / current one naming a superseded total or
  # folder without marking it so. Add a total here when the box moves.
  stale_runs <- c("72,027", "94,376", "96,118", "99,873", "pooled-CPUE-canonical-2024-25", "pooled-CPUE-IMP-R4")
  claim_re <- "authoritative|method of record|current reference|reference run|the current run"
  label_re <- paste0("supersed|SUPERSED|historical|HISTORICAL|previous|used to|replaced|replaces|was the|until |",
                     "before |against|first render|rendered 2026-09-27|R4|what this replaced|the run it|moves|moved|from 9|ladder")
  bad_claim <- unlist(lapply(unique(grep("/archive/", live, value = TRUE, invert = TRUE)), function(f) {   # an archive is a record
    if (!file.exists(f)) return(character(0))
    L0 <- readLines(f, warn = FALSE)
    hit <- vapply(L0, function(l) grepl(claim_re, l, ignore.case = TRUE) &&
                    any(vapply(stale_runs, function(t) grepl(t, l, fixed = TRUE), logical(1))) && !grepl(label_re, l),
                  logical(1))
    if (any(hit)) sprintf("%s:%d", f, which(hit)) else character(0)
  }))
  chk("pointers: no live line calls a superseded run (94,376, 96,118, their folders) authoritative or current",
      length(bad_claim) == 0, sprintf("(%s)", paste(bad_claim, collapse = ", ")))

  # section 1v: one direction, stated. Moved to VALIDATION_CAMPAIGN.md on 2026-09-12, so the
  # assertion follows it and also checks that the status document still routes the reader.
  vc <- "07_documentation/development_notes/VALIDATION_CAMPAIGN.md"
  chk("pointers: the campaign narrative has its own document and the status document routes to it",
      file.exists(vc) &&
      any(grepl("VALIDATION_CAMPAIGN.md", L, fixed = TRUE)) &&
      any(grepl("^## 1b to 1v\\. The validation campaign: MOVED", L)))
  VL  <- if (file.exists(vc)) readLines(vc, warn = FALSE) else character(0)
  i1v <- grep("^## 1v\\.", VL)
  chk("pointers: section 1v survived the move and every 1x letter came with it",
      length(i1v) == 1 &&
      length(grep("^## 1[b-v]\\.", VL)) >= 20 &&
      !any(grepl("^## 1[b-v]\\.", L)))
  s1v <- if (length(i1v) == 1) paste(VL[seq(i1v[1], length(VL))], collapse = "\n") else ""
  chk("pointers: section 1v states which way its PE/BSS percentages run",
      grepl("PE RELATIVE TO BSS", s1v, fixed = TRUE) &&
      !grepl("PE 35,292 against BSS 34,837, -1.3%", s1v, fixed = TRUE))
  # 2026-09-14: RETIRED AND REPLACED, not deleted. This forbade any campaign section title
  # dated 2026-09-13..19, which caught the original defect (1u dated in the future) by
  # hardcoding the window the branch happened to be in. Section 1w is legitimately dated
  # 2026-09-14, so the check now fires on correct content, which is the signature of an
  # assertion written against a snapshot rather than an invariant. Harness section 70
  # already asserts the real invariants and maintains itself: no dated marker in any of the
  # three chronology documents may postdate that document's own "Last updated" line or the
  # branch tip's author date, and the section letters may not run backwards except against
  # a declared allow-list. This row asserts that the stronger check is still present, so
  # retiring this one cannot silently reduce coverage.
  chk("pointers: the campaign's future-date guard is section 70's, which maintains itself",
      { h <- paste(readLines("06_diagnostics/test_improvements_2026-08-25.R", warn = FALSE),
                   collapse = "\n")
        # 2026-09-29: the patterns are built from two pieces, because written whole they
        # appeared in this very check and matched themselves, so the row could never fail.
        grepl(paste0("chk(sprintf(\"chronology: no marker in %s ", "postdates its own 'Last updated'"), h, fixed = TRUE) &&
        grepl(paste0("chk(sprintf(\"chronology: no marker ", "postdates the branch tip's author date"), h, fixed = TRUE) },
      "(the hardcoded 2026-09-13..19 window it replaced fired on Section 1w, which is correct)")
})

# ---------------------------------------------------------------------------
# 63. THE TWO FRONT DOORS (2026-09-12). CLAUDE.md and the root README are what a new reader
#     (or a new agent session) reads first, and both described a model two turnovers, one
#     crabbing fraction and one census ago. The worst of it was operational rather than
#     numeric: neither said that six of the nine model workbooks are BUILT from raw/, so the
#     documented workflow was "place input data in 04_input_files/", which for those six means
#     an edit that the next build silently discards.
# ---------------------------------------------------------------------------
local({
  cl <- paste(readLines("07_documentation/CLAUDE.md", warn = FALSE), collapse = "\n")
  rm_ <- paste(readLines("README.md", warn = FALSE), collapse = "\n")

  chk("front doors: neither offers the retired 1.7 / 1.2 turnovers as current values",
      !grepl("`L = tau_shore` (~1.7) / `tau_boat` (~1.2)", cl, fixed = TRUE) &&
      grepl('tau_shore_prior_mu = "derived"', cl, fixed = TRUE) &&
      grepl('tau_boat_prior_mu = "calibration"', cl, fixed = TRUE))
  chk("front doors: CLAUDE.md describes the census as a commercial CENSUS plus a charter EXPANSION",
      grepl("EXACT CENSUS", cl, fixed = TRUE) && grepl('charter_frame = "roster"', cl, fixed = TRUE) &&
      !grepl("a day-type-stratified **census expansion** of the daily vessel tally", cl, fixed = TRUE))
  chk("front doors: the built workbooks are flagged in BOTH, with the builder named",
      grepl("build_all_inputs.R", cl, fixed = TRUE) && grepl("build_all_inputs.R", rm_, fixed = TRUE) &&
      grepl("Never hand-edit one", cl, fixed = TRUE))
  chk("front doors: CLAUDE.md carries the multi-season span and the PE fill levers",
      grepl("pot_closures", cl, fixed = TRUE) && grepl("season_totals.csv", cl, fixed = TRUE) &&
      grepl("pe_empty_effort_stratum", cl, fixed = TRUE) && grepl("pe_effort_strata.R", cl, fixed = TRUE))
  chk("front doors: the README no longer tells the reader to place CSV inputs",
      !grepl("effort_combined.csv", rm_, fixed = TRUE) &&
      !grepl("interview_combined.csv", rm_, fixed = TRUE) &&
      !grepl("wes_commercial_tally.csv", rm_, fixed = TRUE))
  # A hard-coded assertion count goes stale every time an assertion is added, and a tripwire
  # that fires on every patch is a tripwire nobody reads. So assert only that the figure is
  # not one of the long-dead ones, and make the doc say the run prints the real number.
  chk("front doors: CLAUDE.md does not quote a long-superseded harness count",
      !grepl("452 assertions", cl, fixed = TRUE) && !grepl("725 assertions", cl, fixed = TRUE) &&
      grepl("the run prints the count", cl, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 64. THE HOW-TO DOCUMENTED A PARAMETER THE POOLED MODEL RETIRED (2026-09-12). The
#     effort-overdispersion HOWTO stated mu_i = lambda_E * R_T for the trailer stream. POOL-1
#     removed R_T from crab_bss_pooled.stan (a bernoulli on a vector of literal ones had
#     pinned it at 1.00) and replaced it with R_G_boat, under which the mean is
#     lambda_E / R_G_boat. The CODE was never wrong -- bss_trailer_par() resolves whichever
#     parameter the fit declares -- but a reader following the document would look for a
#     parameter that does not exist and, substituting the one that does, invert the expansion.
# ---------------------------------------------------------------------------
local({
  h <- "07_documentation/effort_overdispersion_diagnostic_HOWTO.md"
  chk("howto: the overdispersion HOWTO exists", file.exists(h)); if (!file.exists(h)) return(invisible(NULL))
  d <- paste(readLines(h, warn = FALSE), collapse = "\n")
  chk("howto: the trailer mean is stated as a multiplier, and the reciprocal is called out",
      !grepl("(R is `R_G` for gear, `R_T` for trailer)", d, fixed = TRUE) &&
      grepl("1 / R_G_boat", d, fixed = TRUE) && grepl("bss_trailer_multiplier", d, fixed = TRUE))
  # and the premise: the pooled Stan really has retired R_T, so the doc is right to say so
  st <- paste(readLines("02_stan_models/crab_bss_pooled.stan", warn = FALSE), collapse = "\n")
  decl <- grep("^\\s*real(<[^>]*>)?\\s+R_T\\s*;", readLines("02_stan_models/crab_bss_pooled.stan", warn = FALSE))
  chk("howto: the premise holds -- the pooled Stan declares R_G_boat and no longer declares R_T",
      length(decl) == 0 && grepl("real<lower=0> R_G_boat;", st, fixed = TRUE))
  chk("howto: the resolver the doc points at handles both parameters and inverts only R_G_boat",
      { b <- paste(readLines("03_R_functions/bss_trailer_expansion.R", warn = FALSE), collapse = "\n")
        grepl('if ("R_T" %in% mp) "R_T"', b, fixed = TRUE) &&
        grepl('else if ("R_G_boat" %in% mp) "R_G_boat"', b, fixed = TRUE) &&
        grepl('if (identical(par_name, "R_T")) v else 1 / v', b, fixed = TRUE) })
})

# ---------------------------------------------------------------------------
# 65. THE SHIPPED WINDOW IS NOT THE AUTHORITATIVE WINDOW (2026-09-12, D28). run_config.R
#     ships the 2023-25 two-season span; the ladder's WINDOW block pins every rung to the
#     single 2024-25 season. Both are deliberate, and the combination means a fresh clone
#     running run_estimation.R does NOT reproduce the box. The claim "R4 IS run_config.R"
#     is true of the model levers and false of the window, so the status document has to
#     say which. Assert the discrepancy is DECLARED, not that it is absent: if someone
#     later rolls run_config back to 2024-25 these still pass, because then the window in
#     the pin and the window in the config agree.
# ---------------------------------------------------------------------------
local({
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  t <- readLines("06_diagnostics/run_improvements_2026-09-08.R", warn = FALSE)
  i <- grep("^WINDOW <- list\\(", t)[1]
  chk("window: the ladder runner has a WINDOW pin", !is.na(i))
  if (is.na(i)) return(invisible(NULL))
  j <- i; while (!grepl("^\\s*estimate_red_rock", t[j])) j <- j + 1L
  ev <- new.env(); eval(parse(text = paste(t[i:j], collapse = "\n")), envir = ev)
  W <- ev$WINDOW
  same <- identical(as.character(W$season_filter), as.character(rc$season_filter)) &&
          identical(as.character(W$est_date_start), as.character(rc$est_date_start)) &&
          identical(as.character(W$est_date_end),   as.character(rc$est_date_end))
  ps <- paste(readLines("07_documentation/development_notes/PIPELINE_STATUS.md", warn = FALSE), collapse = "\n")
  chk("window: if the pin and run_config disagree, the status document says so in the box",
      same || (grepl("WHAT R4 DOES NOT MATCH IS THE SHIPPED WINDOW", ps, fixed = TRUE) &&
               grepl("D8", ps, fixed = TRUE)),
      sprintf("(pin %s..%s [%s]; run_config %s..%s [%s])",
              W$est_date_start, W$est_date_end, paste(W$season_filter, collapse = "+"),
              rc$est_date_start, rc$est_date_end, paste(rc$season_filter, collapse = "+")))
  chk("window: the box no longer claims R4 IS run_config without qualifying the window",
      !grepl("**R4 is not a variant of production: `D_R4` IS\n> `run_config.R`**", ps, fixed = TRUE))
  # the census frame guard is a warning, not a stop: assert which, so a change is deliberate
  cc <- paste(readLines("03_R_functions/estimate_comm_charter.R", warn = FALSE), collapse = "\n")
  # 2026-09-13: re-anchored. This asserted one exact sprintf() spelling, so wiring the
  # condition into the report (B30) broke it without anything being wrong. What has to
  # hold is the SUBSTANCE: it warns rather than stopping, and the condition reaches the
  # reader. Matt settled the open question: "a warning is fine - included in the html report."
  chk("window: the missing-census-frame condition WARNS and is not a stop (Matt's decision)",
      grepl("NO CENSUS FRAME", cc, fixed = TRUE) &&
      grepl("warning(paste0(\"estimate_comm_charter(): \", .m), call. = FALSE)", cc, fixed = TRUE) &&
      !grepl("stop(paste0(\"estimate_comm_charter()", cc, fixed = TRUE))
  chk("window: and the condition reaches the REPORT, not only the console",
      grepl("frame_warnings", cc, fixed = TRUE) &&
      all(vapply(c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd",
                   "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"),
                 function(f) grepl("census_frame_warnings.csv",
                                   paste(readLines(f, warn = FALSE), collapse = "\n"),
                                   fixed = TRUE), logical(1))))
})

# ---------------------------------------------------------------------------
# 66. run_config.R IS THE CANONICAL RUN, AND IT IS ORDERED (2026-09-12, D28 closed).
#     Two properties, both of which had failed silently before. (a) The shipped window is
#     the window of the authoritative run, so a clone that runs run_estimation.R
#     reproduces the box instead of the blocked 2023-25 span. Asserted by comparing the
#     config against the ladder's own WINDOW pin, key for key, so either drifting is
#     caught. (b) The file is ordered so a reader who wants to know what the model DOES
#     can stop at the end of section 2: the method keys come before the diagnostics
#     banner and the experiment levers come after it.
# ---------------------------------------------------------------------------
local({
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  src <- readLines("run_config.R", warn = FALSE)

  # ---- (a) the canonical window -------------------------------------------
  # 2026-09-29 (B50): Matt ships model = "both", which renders the pooled headline first and then
  # its gear-resolved cross-check; either keeps the pooled model as the headline estimator.
  chk("canonical: the shipped model renders the pooled headline estimator (pooled, or both with its cross-check)",
      isTRUE(e$model %in% c("pooled", "both")))
  # 2026-09-13: run_weather went with the weather module. run_config.R must define exactly
  # `model` and `run_config` and nothing else, so a stray top-level object cannot creep back.
  chk("canonical: run_config.R defines only `model` and `run_config` (run_weather is gone)",
      setequal(ls(e), c("model", "run_config")),
      sprintf("(defines: %s)", paste(sort(ls(e)), collapse = ", ")))
  chk("canonical: the nine per-season keys are the single 2024-25 season",
      identical(rc$est_date_start, "2024-09-16") && identical(rc$est_date_end, "2025-09-15") &&
      identical(rc$season_filter, "2024-25") &&
      is.null(rc$pot_closures) && is.null(rc$census_windows) &&
      identical(rc$pot_closure_start, "2024-09-16") && identical(rc$pot_closure_end, "2024-11-30") &&
      identical(rc$pot_open_date, "2024-12-01") &&
      identical(rc$census_start_date, "2024-12-01") && identical(rc$census_end_date, "2025-02-08") &&
      identical(rc$commercial_opener, "2025-02-11"))
  t <- readLines("06_diagnostics/run_improvements_2026-09-08.R", warn = FALSE)
  i <- grep("^WINDOW <- list\\(", t)[1]
  j <- i; while (!grepl("^\\s*estimate_red_rock", t[j])) j <- j + 1L
  ev <- new.env(); eval(parse(text = paste(t[i:j], collapse = "\n")), envir = ev)
  W <- ev$WINDOW
  # run_weather is skipped because it is no longer a config key at all (it went with the
  # weather module on 2026-09-13). It STAYS in the ladder's WINDOW pin deliberately: the
  # pin's key names feed .cfg_fingerprint(), so dropping it would change every stage digest
  # and RESUME would refuse the five committed rung folders, costing ~16 h of refitting to
  # remove one inert key from a dated runner. estimate_red_rock is skipped for the same
  # reason: removed from run_config on 2026-09-28 (B44), kept in the pin.
  drift <- Filter(function(k) !identical(rc[[k]], W[[k]]), setdiff(names(W), c("run_weather", "estimate_red_rock")))
  chk("canonical: run_config AGREES WITH THE LADDER'S WINDOW PIN on every key it pins",
      length(drift) == 0, sprintf("(drifted: %s)", paste(drift, collapse = ", ")))
  # The paste-ready alternative windows live in this file as COMMENTS (restored 2026-09-13
  # at Matt's request, alongside the copies in NEW_SEASON_GUIDE 7.1), so the text test has
  # to be about what the file SETS, not about what it mentions. Read the parsed config for
  # the substantive property, and check the mentions are commented out.
  chk("canonical: nothing in the file SETS the blocked two-season span",
      length(rc$season_filter) == 1L && !is.list(rc$pot_closures) && !is.list(rc$census_windows))
  chk("canonical: the span survives as a commented block, not as live code",
      { hits <- grep('season_filter *= *c\\("2023-24"', src, value = TRUE)
        length(hits) >= 1 && all(grepl("^\\s*#", hits)) })
  # the desk check recorded in the box: these are the values the shipped config resolves to,
  # pinned so a lever that changes one of them cannot pass silently.
  # 2026-09-28: the authoritative run IS the shipped configuration (rendered by run_estimation.R), so
  # the desk check became a render; the values it resolves to stay pinned in the box.
  chk("canonical: the box records the desk check that the shipped config reproduces the run's inputs",
      { ps <- paste(readLines("07_documentation/development_notes/PIPELINE_STATUS.md", warn = FALSE), collapse = "\n")
        f <- gsub("[ \n>]+", " ", ps)
        # 2026-10-04 (B65): the shipped file IS the authoritative run again (the A33/A34 render); the
        # replaced run's box keeps its resolved inputs, which A33 and A34 do not change
        grepl("Apart from A33's one key, the shipped configuration IS this run, so its resolved inputs are the shipped ones", f, fixed = TRUE) &&
        grepl("Matt's render of `run_config.R` exactly as shipped (`run_estimation.R --model both`", f, fixed = TRUE) &&
        grepl("2.4771", f, fixed = TRUE) && grepl("3.0300", f, fixed = TRUE) &&
        grepl("commercial 6,405 plus charter 2,133 = 8,538", f, fixed = TRUE) &&
        grepl("none of them a missing frame", f, fixed = TRUE) })
  chk("canonical: the header names the authoritative run and points at the one box",
      any(grepl("05_output/20260928/pooled-CPUE-2024-25", src, fixed = TRUE)) &&
      any(grepl("99,873", src, fixed = TRUE)) &&
      any(grepl("PIPELINE_STATUS.md", src, fixed = TRUE)))

  # ---- (b) the file is ordered ---------------------------------------------
  ban <- function(pat) { h <- grep(pat, src); if (length(h)) h[1] else NA_integer_ }
  b1 <- ban("^  # 1\\. THE RUN"); b2 <- ban("^  # 2\\. THE METHOD OF RECORD")
  b3 <- ban("^  # 3\\. DIAGNOSTICS REPORTED BESIDE"); b4 <- ban("^  # 4\\. LEVERS FOR SENSITIVITY")
  b5 <- ban("^  # 5\\. GEAR-RESOLVED MODEL ONLY")
  chk("ordered: the five section banners are present and in order",
      all(is.finite(c(b1, b2, b3, b4, b5))) && all(diff(c(b1, b2, b3, b4, b5)) > 0))
  keyline <- function(k) { h <- grep(sprintf("^  %s[ ]*=", k), src); if (length(h)) h[1] else NA_integer_ }
  # every key in the file is at indent 2 exactly once, so a key's line locates its section
  in_method <- c("shore_effort_unit", "tau_shore_prior_mu", "tau_boat_prior_mu", "shared_tau",
                 "use_osp_boat_counts", "osp_scale_is_tau", "crab_fraction_dynamic",
                 "crab_fraction_strata", "use_osp_crab_lower", "estimate_catch_zi",
                 "ar_max_resolution", "pe_empty_stratum", "pe_empty_effort_stratum",
                 "pe_variance", "census_expansion", "charter_frame", "use_ie_day_length",
                 "marine_hazard_mode", "marine_hazard_manual_boat", "marine_hazard_gear_regimes", "marine_hazard_codes")   # 2026-09-27: A30 ADOPTED, section 2.10
  in_diag  <- c("diagnose_incomplete_trips", "diagnose_tau_sensitivity", "ar_rung_adequacy",
                "save_ppc_draws", "run_fishery_spillover_diag")
  in_lever <- c("ar_force", "bss_sampler_override", "ar_escalate", "opener_covariate_mode",
                "razor_dig_mode", "estimate_cpue_density", "collapse_mu_hier", "estimate_B1_C",
                "bar_restriction_impute", "marine_hazard_candidates_boat", "marine_hazard_winter_months")   # 2026-09-25: A30, section 4.4b (the experiment surface)
  in_gear  <- c("gear_resolved_G", "gear_share_dirichlet", "ar_adaptive", "use_boat_ie")
  chk("ordered: every method-of-record key sits in section 2",
      all(vapply(in_method, function(k) { l <- keyline(k); is.finite(l) && l > b2 && l < b3 }, logical(1))),
      sprintf("(stray: %s)", paste(Filter(function(k) { l <- keyline(k); !(is.finite(l) && l > b2 && l < b3) }, in_method), collapse = ", ")))
  chk("ordered: every production diagnostic sits in section 3",
      all(vapply(in_diag, function(k) { l <- keyline(k); is.finite(l) && l > b3 && l < b4 }, logical(1))))
  chk("ordered: every sensitivity / experiment lever sits in section 4",
      all(vapply(in_lever, function(k) { l <- keyline(k); is.finite(l) && l > b4 && l < b5 }, logical(1))))
  chk("ordered: every gear-resolved-only toggle sits in section 5",
      all(vapply(in_gear, function(k) { l <- keyline(k); is.finite(l) && l > b5 }, logical(1))))
  chk("ordered: the window and the sampler sit in section 1, ahead of any method key",
      all(vapply(c("est_date_start", "season_filter", "pot_open_date", "bss_seed", "bss_chains",
                   "effort_file", "interview_file"),
                 function(k) { l <- keyline(k); is.finite(l) && l > b1 && l < b2 }, logical(1))))
  chk("ordered: the header tells the reader they can stop at the end of section 2",
      any(grepl("you can stop here", src, fixed = TRUE)) &&
      any(grepl("HOW TO READ THIS FILE", src, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 67. THE MONTHLY PE SHARE HAD LOST THE CRABBING FRACTION (2026-09-12). The boat branch of
#     pe_monthly_effort_share() omitted the per-day f that run_pe_pooled() applies. That was
#     EXACT while f was the scalar 0.30, because the share is normalized and a constant
#     multiplier cancels; the dynamic monthly f adopted 2026-09-08 made it stop cancelling
#     silently, with nothing at the call site changing. Measured on 2024-25: September
#     over-weighted 2.39x, December pulled to 0.36x. No total moves. These assert the
#     ALIGNMENT of the two formulae, not a number, so the next lever cannot re-open the gap.
# ---------------------------------------------------------------------------
local({
  # RETIRED 2026-09-28 (B44). This section pinned the per-day f in BOTH copies of the PE's
  # daily-effort formula (pe_monthly_effort_share() and .srd_monthly_share()). Both copies
  # were removed: the monthly split now reads the PE's own strata (pe_monthly_split(),
  # section 81), which carry f, the turnover and every other per-day factor already, so
  # there is no second formula for a factor to go missing from. What stays asserted is
  # that the copies do not come back.
  chk("monthly share: pe_monthly_effort_share() is retired (no definition, no call)",
      !file.exists("03_R_functions/pe_monthly_effort_share.R") &&
      !any(grepl("pe_monthly_effort_share(", unlist(lapply(c(list.files("03_R_functions", full.names = TRUE),
                                                             list.files("01_BSS_models", pattern = "Rmd$", full.names = TRUE)),
                                                           function(f) { x <- readLines(f, warn = FALSE); x[!grepl("^\\s*#", x)] })),
                 fixed = TRUE)))
  ph <- readLines("03_R_functions/run_pe_pooled.R", warn = FALSE)
  chk("monthly share: run_pe_pooled's header no longer claims a weighted mean of daily ratios",
      !any(grepl("weighted mean of daily ratios) is intentionally left", ph, fixed = TRUE)) &&
      any(grepl("has\n# been RATIO-OF-SUMS ever since", paste(ph, collapse = "\n"), fixed = TRUE)) |
      any(grepl("RATIO-OF-SUMS ever since", ph, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 68. METHOD v2.0 IS THE METHOD OF RECORD (2026-09-12, A28, closing D7). The break with v1.0
#     is that v2.0 is LIVE: it is written against the code as it runs and tracks it. So the
#     assertions here are not about prose, they are about the document agreeing with the
#     code and with the run it names. A frozen document could be outrun silently; this one
#     fails the harness when it is.
# ---------------------------------------------------------------------------
local({
  mp <- "07_documentation/BSS-GH-pooled-CPUE-model-documentation.md"
  mg <- "07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md"
  chk("method v2.0: both method documents exist", file.exists(mp) && file.exists(mg))
  if (!file.exists(mp) || !file.exists(mg)) return(invisible(NULL))
  P <- readLines(mp, warn = FALSE); G <- readLines(mg, warn = FALSE)
  ps <- paste(P, collapse = "\n"); gs <- paste(G, collapse = "\n")

  chk("method v2.0: the pooled document is v2.0, live rather than frozen, and says so",
      grepl("Method Version 2.0", ps, fixed = TRUE) &&
      grepl("this document is **not frozen**", ps, fixed = TRUE) &&
      grepl("Operational, **not published**", ps, fixed = TRUE))
  chk("method v2.0: the gear document is framework v6.0 and defers rather than duplicating",
      grepl("framework v6.0", gs, fixed = TRUE) &&
      grepl("READ THE POOLED DOCUMENT FIRST", gs, fixed = TRUE) &&
      grepl("describes only what is DIFFERENT", gs, fixed = TRUE))
  chk("method v2.0: v1.0 is archived, unaltered, under a banner with the difference table",
      file.exists("07_documentation/archive/method-v1.0-pooled-CPUE.md") &&
      file.exists("07_documentation/archive/method-v1.0-gear-resolved-CPUE.md") &&
      { a <- paste(readLines("07_documentation/archive/method-v1.0-pooled-CPUE.md", warn = FALSE), collapse = "\n")
        grepl("ARCHIVED. THIS IS METHOD v1.0", a, fixed = TRUE) &&
        grepl("Method v1.0 | Method v2.0", a, fixed = TRUE) })

  # THE POINT OF A LIVE DOCUMENT: its reference run must exist and its total must match it.
  ref <- regmatches(ps, regexpr("05_output/[0-9]{8}/[A-Za-z0-9._-]+", ps))
  chk("method v2.0: the pooled document names a reference run that exists on disk",
      length(ref) == 1 && dir.exists(ref), sprintf("(names %s)", paste(ref, collapse = "")))
  chk("method v2.0: ITS HEADLINE AGREES WITH THAT RUN's port_total csv",
      { f <- file.path(ref, "port_total_Dungeness_Kept.csv")
        if (!length(ref) || !file.exists(f)) FALSE else {
          d <- utils::read.csv(f, stringsAsFactors = FALSE)
          r <- d[grepl("^Expected", d[[2]]), , drop = FALSE]
          fm <- function(x) formatC(round(as.numeric(x)), format = "d", big.mark = ",")
          nrow(r) == 1 && all(vapply(c(r$BSS_median, r$BSS_lo95, r$BSS_hi95),
                                     function(v) grepl(fm(v), ps, fixed = TRUE), logical(1))) } })

  # the method document must describe the levers run_config actually ships
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("method v2.0: the document's stated method matches the shipped levers",
      identical(rc$crab_fraction_dynamic, TRUE)    && grepl("logit random walk", ps, fixed = TRUE) &&
      identical(rc$use_osp_boat_counts, TRUE)      && grepl("OSP daily port count", ps, fixed = TRUE) &&
      identical(rc$osp_scale_is_tau, TRUE)         && grepl("osp_scale_is_tau = 1", ps, fixed = TRUE) &&
      identical(rc$shared_tau, TRUE)               && grepl("shared_tau = 1", ps, fixed = TRUE) &&
      identical(rc$tau_shore_prior_mu, "derived")  && grepl('tau_shore_prior_mu = "derived"', ps, fixed = TRUE) &&
      identical(rc$tau_boat_prior_mu, "calibration") && grepl('tau_boat_prior_mu = "calibration"', ps, fixed = TRUE) &&
      identical(rc$estimate_catch_zi, TRUE)        && grepl("zero-inflated", ps, fixed = TRUE) &&
      identical(rc$census_expansion, "none")       && grepl('census_expansion = "none"', ps, fixed = TRUE) &&
      identical(rc$charter_frame, "roster")        && grepl("charter trip roster", ps, fixed = TRUE) &&
      identical(rc$pe_empty_effort_stratum, "local_day_type") &&
        grepl('pe_empty_effort_stratum` | `"local_day_type"', ps, fixed = TRUE))
  # the OSP crabbing-only column: the one piece of the method that is built and waiting on
  # data, and the one a reader could most easily mis-wire. Assert all four facts separately.
  .flat <- gsub("[ \n]+", " ", ps)
  # 2026-10-04 (A33): the column was delivered (B64) and is ON; the document must say so, keep
  # the lower-bound warning, and carry the A33 rule.
  chk("method v2.0: it records that the OSP crabbing-only column is DELIVERED and ON (A33), is a lower bound read as f(1 - c), and must not be wired in as f",
      identical(rc$use_osp_crab_lower, TRUE) &&
      grepl("**(a) is in hand**", .flat, fixed = TRUE) &&
      grepl("**(b) was delivered on 2026-10-02**", .flat, fixed = TRUE) &&
      grepl("lower bound** on the vessels that did any crabbing, not the crabbing fraction itself", .flat, fixed = TRUE) &&
      grepl("do not wire the crab-only column in as if it were `f`", .flat, fixed = TRUE) &&
      grepl("The OSP crabbing-only stream: ON since 2026-10-04 (A33)", .flat, fixed = TRUE) &&
      grepl("Rule A33, written before the render", .flat, fixed = TRUE) &&
      !grepl("built, tested and inert (Section 14.3)", .flat, fixed = TRUE))
  chk("method v2.0: the gate thresholds in the document are the gate's actual defaults",
      { g <- paste(readLines("03_R_functions/bss_convergence_gate.R", warn = FALSE), collapse = "\n")
        grepl("max_impact_sd    = 0.10", g, fixed = TRUE) && grepl("`< 0.10`", ps, fixed = TRUE) &&
        grepl("max_div_fraction = 0.05", g, fixed = TRUE) && grepl("`< 5%`", ps, fixed = TRUE) &&
        grepl("rhat_threshold   = 1.01", g, fixed = TRUE) && grepl("`< 1.01`", ps, fixed = TRUE) &&
        grepl("neff_threshold   = 400", g, fixed = TRUE)  && grepl("`> 400`", ps, fixed = TRUE) })
  chk("method v2.0: the AR caps in the document are the caps run_config ships",
      identical(rc$ar_max_resolution$pooled$shore$all_gear, "weekly") &&
      identical(rc$ar_max_resolution$pooled$shore$pot_closure, "biweekly") &&
      identical(rc$ar_max_resolution$pooled$private_boat, "monthly") &&
      grepl("| pooled | biweekly | **weekly** | monthly |", ps, fixed = TRUE))
  chk("method v2.0: it carries the caveats a reviewer will ask about, not just the result",
      grepl("Moored private boats are outside the effort frame", ps, fixed = TRUE) &&
      grepl("is NOT a full predictive distribution", ps, fixed = TRUE) &&
      grepl("asserts a MEDIAN, not a mean", ps, fixed = TRUE) &&
      grepl("MEANS TWO DIFFERENT THINGS under the same name", ps, fixed = TRUE))
  chk("method v2.0: neither document uses an em dash (the stated convention)",
      !any(grepl("—", P)) && !any(grepl("—", G)))
})

# ---------------------------------------------------------------------------
# 69. THE IN-CODE DOCUMENTATION SWEEP (2026-09-12). A full audit of the comments and
#     headers in the executable files found 34 stale claims. The dominant pattern was NOT
#     the crabbing fraction: it was the v7.7 SHORE UNIT MOVE, which seven places still
#     described as an effective day length in hours, each of them a defect already logged
#     elsewhere in the repo. The pooled driver's prose was consistently behind the helpers
#     and behind the gear driver, which is what happens when a patch updates the code and
#     the other track's comments and stops.
#
#     These assertions pin the ones that would have made a reader compute or quote
#     something wrong, and they are written against the SHIPPED VALUES so they fail again
#     if a lever moves without the comment moving with it.
# ---------------------------------------------------------------------------
local({
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  rd  <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  flat <- function(x) gsub("[ \n]+", " ", x)

  # (1) the gate backstop quoted in the report prose must be the value the driver passes
  drv <- rd("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd")
  chk("in-code docs: the pooled report quotes the gate's ACTUAL divergence backstop",
      grepl("max_divergence_fraction        = 0.05", drv, fixed = TRUE) &&
      grepl("the divergent fraction is under the 5% backstop", drv, fixed = TRUE) &&   # 2026-09-29 report wording
      !grepl("under the 15 percent backstop", drv, fixed = TRUE))

  # (2) the v7.7 shore unit: no executable file may describe the shore expansion as hours
  bad <- character()
  for (f in c(list.files("03_R_functions", pattern = "\\.R$", full.names = TRUE),
              list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE),
              list.files("02_stan_models", pattern = "\\.stan$", full.names = TRUE))) {
    if (grepl("weather_adjusted", f)) next           # the stale fork, documented as stale
    t <- flat(rd(f))
    if (grepl("Shore is unchanged: lambda_E = crabbers, h = crabber-hours", t, fixed = TRUE) ||
        grepl("Shore is unchanged (crabber-hours, E_scale = 1).", t, fixed = TRUE) ||
        grepl("SHORE : effort = crabbers x effective day length (hours). days$day_length", t, fixed = TRUE) ||
        grepl("SHORE : effort = crabbers x effective day length (hours). days$day_length", t, fixed = TRUE) ||
        grepl("shore: observation = crabber-hours; predicted = lambda_E * L (hours)", t, fixed = TRUE) ||
        grepl("shore keeps crabber-hours (mean_count * crabbers_per_gear * day_length)", t, fixed = TRUE))
      bad <- c(bad, basename(f))
  }
  chk("in-code docs: no executable file still states the RETIRED pre-v7.7 shore formula",
      length(bad) == 0, sprintf("(%s)", paste(bad, collapse = ", ")))

  # (3) the runtime L label is read from the spec, not hard-coded per population
  pp <- rd("03_R_functions/prep_bss_crab_pooled.R")
  chk("in-code docs: the run log's L label comes from bss_effort_spec, not a hard-coded string",
      grepl("eff_spec$L_unit %||%", pp, fixed = TRUE) &&
      !grepl('if(is_shore) "effective day length, hours" else', pp, fixed = TRUE))

  # (4) the retired turnover literals must not be presented as current anywhere
  lit <- character()
  for (f in c(list.files("03_R_functions", pattern = "\\.R$", full.names = TRUE),
              "02_stan_models/crab_bss_pooled.stan", "02_stan_models/crab_bss_gear_resolved.stan")) {
    t <- flat(rd(f))
    if (grepl("the deployment turnover (~1.2)", t, fixed = TRUE) ||
        grepl("tau_boat, the gear-deployment turnover (~1.2)", t, fixed = TRUE) ||
        grepl("tau_boat prior (~1.2) tension; validate by run before using", t, fixed = TRUE))
      lit <- c(lit, basename(f))
  }
  chk("in-code docs: the retired 1.2 boat turnover centre is nowhere presented as current",
      length(lit) == 0, sprintf("(%s)", paste(lit, collapse = ", ")))

  # (5) shipped-value comments must match the shipped values
  chk("in-code docs: osp_scale_is_tau and shared_tau comments say PRODUCTION SHIPS 1",
      identical(rc$osp_scale_is_tau, TRUE) && identical(rc$shared_tau, TRUE) &&
      all(vapply(c("02_stan_models/crab_bss_pooled.stan", "02_stan_models/crab_bss_gear_resolved.stan"),
                 function(f) { t <- flat(rd(f))
                   grepl("PRODUCTION SHIPS 1", t, fixed = TRUE) &&
                   grepl("production ships 1 since 2026-09-01", t, fixed = TRUE) &&
                   !grepl("0 (default) keeps the free kappa_OSP scale", t, fixed = TRUE) }, logical(1))))
  chk("in-code docs: the ZI block is no longer labelled OFF by default",
      identical(rc$estimate_catch_zi, TRUE) &&
      !grepl("Zero-inflated catch likelihood (2026-09-02, OFF by default)",
             rd("02_stan_models/crab_bss_pooled.stan"), fixed = TRUE))
  chk("in-code docs: both PE files name the SHIPPED empty-stratum CPUE fill",
      identical(rc$pe_empty_stratum, "local") &&
      all(vapply(c("03_R_functions/run_pe_pooled.R", "03_R_functions/run_pe_gear.R"),
                 function(f) { t <- flat(rd(f))
                   !grepl('pe_empty_stratum = "pooled" (default)', t, fixed = TRUE) &&
                   !grepl('pe_empty_stratum: "pooled" (default)', t, fixed = TRUE) }, logical(1))))
  chk("in-code docs: the census is described as a census PLUS an expansion, not one expansion",
      identical(rc$census_expansion, "none") && identical(rc$charter_frame, "roster") &&
      !grepl("day-type-stratified census expansion of the daily vessel tally", drv, fixed = TRUE) &&
      !grepl("day-type-stratified commercial / charter\n# charter census expansion",
             rd("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"), fixed = TRUE) &&
      !grepl("Day-type-stratified census expansion of the commercial/charter vessel tally",
             rd("README-R-functions.md"), fixed = TRUE))
  chk("in-code docs: the AR-cap comment names the right population for each cap",
      identical(rc$ar_max_resolution$pooled$private_boat, "monthly") &&
      grepl("is capped at MONTHLY; the shore is capped at WEEKLY", flat(drv), fixed = TRUE) &&
      !grepl("so the boat is capped at weekly", flat(drv), fixed = TRUE))
  chk("in-code docs: census_uncertainty's options and shipped value are stated correctly",
      identical(rc$census_uncertainty, "charter") &&
      grepl('params$census_uncertainty ("charter" SHIPS | "none" |', flat(drv), fixed = TRUE))

  # (6) the opener/f interaction advice must not send the reader to a walk-destroying value
  chk("in-code docs: the opener advice names month_opener, and says plain opener discards the walk",
      grepl('crab_fraction_strata = "month_opener"`, not with `"opener"`', flat(drv), fixed = TRUE) &&
      grepl("anchors EVERY stratum to the level prior and silently discards the walk", flat(drv), fixed = TRUE) &&
      { cf <- flat(rd("03_R_functions/crab_fraction.R"))    # the premise
        grepl("day_type / opener / none -> every stratum anchored, no walk", cf, fixed = TRUE) })

  # (7) the input-workbook count, which the README contradicted mid-sentence
  chk("in-code docs: the input README says SIX of the nine workbooks are built, matching the builder list",
      { r <- rd("04_input_files/README.md")
        b <- readLines("04_input_files/build_all_inputs.R", warn = FALSE)
        i <- grep("^builders <- c\\(", b)[1]; j <- i; while (!grepl("^\\)", b[j])) j <- j + 1L
        nb <- sum(grepl('"build_[a-z_]+\\.R"', b[i:j]))
        nb == 6L && grepl("six of them BUILT", r, fixed = TRUE) &&
        grepl("Six of the nine", r, fixed = TRUE) && !grepl("Seven of the nine", r, fixed = TRUE) })
  chk("in-code docs: the driver README points at the samplers' interview contacts for f, not the blank I/E columns",
      { r <- flat(rd("01_BSS_models/README.md"))
        grepl("built from the samplers' private-boat INTERVIEW CONTACTS", r, fixed = TRUE) &&
        grepl("present but currently BLANK", r, fixed = TRUE) &&
        !grepl("effort_combined.csv", r, fixed = TRUE) })
  # 2026-09-13: the fork is gone, so the claim cannot be made at all. Assert the ABSENCE of
  # the file and of any config surface, which is the stronger property.
  chk("weather removal: the fork, the driver and the config toggle are all gone",
      !file.exists("02_stan_models/crab_bss_pooled_weather_adjusted.stan") &&
      !file.exists("06_diagnostics/BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd") &&
      length(list.files("02_stan_models", pattern = "\\.stan$")) == 2L &&
      { ee <- new.env(); sys.source("run_config.R", envir = ee)
        !("run_weather" %in% ls(ee)) && !("run_weather" %in% names(ee$run_config)) })
  chk("weather removal: the orchestrator is single-path and REFUSES the retired flags",
      { o <- rd("run_estimation.R")
        !grepl("weather_rmd", o, fixed = TRUE) &&
        !grepl("if (isTRUE(run_weather))", o, fixed = TRUE) &&
        grepl('any(c("--weather", "--no-weather") %in% .args)', o, fixed = TRUE) &&
        grepl("were removed with the weather-tide module", o, fixed = TRUE) })
  chk("weather removal: the FINDING is kept live and the module doc is archived",
      file.exists("07_documentation/WEATHER_COVARIATE_ANALYSIS.md") &&
      file.exists("07_documentation/archive/weather-tide-covariate-module-REMOVED.md") &&
      grepl("THE MODULE IS GONE; THIS FINDING IS WHY, AND IT STANDS",
            rd("07_documentation/WEATHER_COVARIATE_ANALYSIS.md"), fixed = TRUE) &&
      grepl("REMOVED 2026-09-13. THE MODULE THIS DOCUMENTS NO LONGER EXISTS",
            rd("07_documentation/archive/weather-tide-covariate-module-REMOVED.md"), fixed = TRUE))
  chk("weather removal: the OPENER covariates are NOT removed (a different mechanism)",
      file.exists("03_R_functions/bss_opener_covariates.R") &&
      "opener_covariate_mode" %in% names(rc) &&
      identical(rc$opener_covariate_mode, "off"))
  chk("weather removal: no live doc still offers the module as something you can run",
      { bad <- Filter(function(f) {
            t <- flat(rd(f))
            grepl("run_weather <- TRUE", t, fixed = TRUE) ||
            grepl("--model pooled --weather", t, fixed = TRUE) ||
            grepl("The three Stan models", t, fixed = TRUE) },
          c("README.md", "07_documentation/CLAUDE.md", "07_documentation/README.md",
            "02_stan_models/README.md", "06_diagnostics/README.md",
            "07_documentation/BSS-GH-pooled-CPUE-model-documentation.md",
            "07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md"))
        length(bad) == 0 }, "(a doc still describes running it)")

  # (8) the sweep runners must not name a superseded total as production
  chk("in-code docs: the sweep runners point at the box rather than naming a production total",
      { a <- flat(rd("06_diagnostics/run_rg_sweep.R")); b <- flat(rd("06_diagnostics/run_tau_sweep.R"))
        !grepl("current production total is 71,513", a, fixed = TRUE) &&
        !grepl("1.2 = the production prior", b, fixed = TRUE) &&
        !grepl("0.3 matches production.", b, fixed = TRUE) &&
        grepl("COMPARE A NEW SWEEP AGAINST THE BOX AT THE TOP OF", a, fixed = TRUE) })
  # (9) 2026-09-28 (B49): the R_G sweep runs the three rungs its header names, dry by default,
  # and reports against its own CSVs rather than a baseline written into the file
  chk("in-code docs: run_rg_sweep.R ships DRY_RUN TRUE, sweeps 1.00 / 1.28 / 1.50, and names no baseline total",
      { a <- flat(rd("06_diagnostics/run_rg_sweep.R"))
        grepl("DRY_RUN <- TRUE", a, fixed = TRUE) && grepl("rg_grid <- c(1.00, 1.28, 1.50)", a, fixed = TRUE) &&
        !grepl("83,035", a, fixed = TRUE) && grepl("rg_sweep_%s_summary.csv", a, fixed = TRUE) })
})

# ---------------------------------------------------------------------------
# 70. THE CHRONOLOGY DOCUMENTS AND THEIR DATES (2026-09-10 work; written 2026-09-13)
#
#     Matt: "The section dates in PIPELINE_STATUS.md run backwards ... 1u's date is in
#     the future relative to today. Several register rows carry 2026-09-13 as well. The
#     git commit dates are the reliable record. Harmless individually, corrosive in a
#     document whose value is being a chronology."
#
#     WHY THIS IS NOT A "MARKER <= ITS OWN COMMIT DATE" CHECK. That invariant was built,
#     run over every tracked .R/.Rmd/.stan/.md via git blame, and REJECTED on the
#     evidence. It reported 355 violations, and the distribution is the proof that the
#     rule is wrong rather than the repository: the minimum offset is one day and there
#     is no zero bucket at all. The cause is measurable. Several of this branch's commits
#     carry a patch filename as their subject, so both clocks are recorded in the same
#     string:
#
#       0001stage5review20260831.patch          git author date 2026-08-31   gap 0 d
#       0001prerunaudit20260903.patch           git author date 2026-09-01   gap 2 d
#       0005-candidate-config-2026-09-07.patch  git author date 2026-09-04   gap 3 d
#       0012-season-portability-2026-09-09.patch git author date 2026-09-05  gap 4 d
#       0013-multiseason-2026-09-10.patch       git author date 2026-09-05   gap 5 d
#       run_improvements_2026-09-08.R (3fa8fde) git author date 2026-09-08   gap 0 d
#
#     The sessions that produced this branch ran on a clock ahead of the repository's,
#     by zero days in late August, growing to five by 2026-09-05, closing to zero on
#     2026-09-08. Those session dates are baked into FILENAMES
#     (run_adoption_2026-09-07.R, adoption-review-2026-09-08.md), so honouring a
#     tolerance would mean renaming files every pointer resolves through. The check also
#     cannot tell a marker from a calendar fact: it flags "2027-06-18" in the Juneteenth
#     fixture, "the season ends 2026-09-15", and a workbook coverage table's 2026-12-31.
#     A check with a 127-entry permanent allow-list that flags data is a check that gets
#     switched off. It is not shipped, and this comment is the record of why.
#
#     WHAT IS ASSERTED INSTEAD, all of it self-contained and about the three documents
#     whose value IS being a chronology:
#       (1) no dated MARKER in them may postdate the date the document itself declares
#           as its last update, which is exactly the failure above;
#       (2) the campaign's section letters must not run backwards, apart from a declared
#           allow-list of inversions, each with its reason;
#       (3) the git-anchor table must cover every section letter;
#       (4) the dating convention must be stated in one place and pointed at from the
#           other two;
#       (5) git-conditional, one `git log` call: nothing in them may claim a date after
#           the newest author date on the branch.
# ---------------------------------------------------------------------------
local({
  rd   <- function(f) readLines(f, warn = FALSE)
  PS   <- "07_documentation/development_notes/PIPELINE_STATUS.md"
  CR   <- "07_documentation/development_notes/CHANGE_REGISTER.md"
  VC   <- "07_documentation/development_notes/VALIDATION_CAMPAIGN.md"
  DOCS <- c(PS, CR, VC)
  chk("chronology: all three documents exist", all(file.exists(DOCS)))

  # A dated MARKER, as opposed to a calendar fact or a data value: a date in a section
  # heading, inside a status bracket, or in a "Last updated" line. Deliberately narrow.
  markers <- function(lines) {
    keep <- grepl("^#{1,4} ", lines) |
            grepl("\\*\\*(NEW|CLOSED|ADOPTED|SUPERSEDED|DONE|FIXED|REMOVED|REJECTED|RUN)[^*]*20[0-9]{2}-[0-9]{2}-[0-9]{2}", lines) |
            grepl("\\[(NEW|CLOSED|ADOPTED|SUPERSEDED|DONE|FIXED|REMOVED) 20[0-9]{2}-[0-9]{2}-[0-9]{2}", lines) |
            grepl("Last updated", lines)
    unlist(regmatches(lines[keep], gregexpr("20[0-9]{2}-[0-9]{2}-[0-9]{2}", lines[keep])))
  }
  declared <- function(lines) {
    l <- grep("Last updated", lines, value = TRUE)[1]
    m <- regmatches(l, regexpr("20[0-9]{2}-[0-9]{2}-[0-9]{2}", l))
    if (length(m) == 1) as.Date(m) else as.Date(NA)
  }

  # (1) no marker may postdate the document's own declared last-update date
  for (f in DOCS) {
    L  <- rd(f); dd <- declared(L)
    mk <- suppressWarnings(as.Date(markers(L)))
    mk <- mk[!is.na(mk)]
    chk(sprintf("chronology: %s declares a last-updated date", basename(f)), !is.na(dd))
    chk(sprintf("chronology: no marker in %s postdates its own 'Last updated' (%s)",
                basename(f), as.character(dd)),
        !is.na(dd) && length(mk) > 0 && all(mk <= dd),
        if (!is.na(dd) && any(mk > dd)) sprintf("(ahead: %s)",
          paste(unique(as.character(mk[mk > dd])), collapse = ", ")) else "")
  }

  # (2) the campaign's section letters must not run backwards.
  #     DECLARED INVERSIONS. Each entry is "<letter> => <reason>". An inversion that is
  #     not declared FAILS; a declared one that has gone away also fails, because the
  #     only ways to remove it are to re-date the sections from git (in which case this
  #     list is what you update) or to invent a date (in which case you should not).
  INVERSION_OK <- c(
    "1p" = paste("1n and 1o were written when the session clock read 2026-09-09 and",
                 "2026-09-10; 1p was written when it read 2026-09-08. In GIT order",
                 "(2026-09-05, 2026-09-05, 2026-09-08) the letters are correct. No",
                 "assignment of session dates removes this without inventing one."))
  L   <- rd(VC)
  hd  <- grep("^## 1[b-z]\\.", L, value = TRUE)   # 2026-09-26: 1w and 1x joined the campaign; 2026-09-27: 1y and 1z
  let <- sub("^## (1[b-z])\\..*$", "\\1", hd)
  lastdate <- function(h) {
    m <- regmatches(h, gregexpr("20[0-9]{2}-[0-9]{2}-[0-9]{2}", h))[[1]]
    if (length(m)) as.Date(tail(m, 1)) else as.Date(NA)
  }
  dts <- as.Date(vapply(hd, function(h) as.character(lastdate(h)), character(1)))
  chk("chronology: the campaign has 25 sections, 1b to 1z, in ascending letter order",
      length(let) == 25 && identical(let, let[order(let)]) &&
      identical(let[1], "1b") && identical(let[length(let)], "1z"))
  chk("chronology: every campaign section title carries a parseable date", !any(is.na(dts)))
  inv <- let[-1][which(diff(as.numeric(dts)) < 0)]
  chk("chronology: the campaign's section dates run backwards ONLY where declared",
      setequal(inv, names(INVERSION_OK)),
      sprintf("(found: %s; declared: %s)", paste(inv, collapse = ","),
              paste(names(INVERSION_OK), collapse = ",")))
  chk("chronology: every declared inversion carries a reason",
      all(nzchar(INVERSION_OK)) && all(nchar(INVERSION_OK) > 40))

  # (3) the git-anchor table must cover every section letter
  tbl <- grep("^> \\| 1[b-z] \\|", L, value = TRUE)
  tl  <- sub("^> \\| (1[b-z]) \\|.*$", "\\1", tbl)
  chk("chronology: the git-anchor table covers every campaign section",
      setequal(tl, let), sprintf("(missing: %s)",
        paste(setdiff(let, tl), collapse = ",")))
  chk("chronology: every git-anchor row names at least one commit hash",
      length(tbl) > 0 && all(grepl("`[0-9a-f]{7}`", tbl)))

  # (4) the convention is stated once and pointed at from the other two
  vcx <- paste(L, collapse = "\n")
  chk("chronology: the campaign holds the dating note, with the measured cause",
      grepl("A NOTE ON THE DATES", vcx, fixed = TRUE) &&
      grepl("ran on a clock that was ahead of the repository", vcx, fixed = TRUE) &&
      grepl("git author date", vcx, fixed = TRUE))
  for (f in c(PS, CR)) {
    t <- gsub("[ \n]+", " ", paste(rd(f), collapse = "\n"))
    chk(sprintf("chronology: %s points at the one dating note rather than restating it",
                basename(f)),
        grepl("A NOTE ON THE DATES", t, fixed = TRUE) &&
        grepl("git author date", t, fixed = TRUE))
  }

  # (5) nothing may claim a date after the newest author date on the branch
  tip <- tryCatch(suppressWarnings(
           system2("git", c("log", "-1", "--format=%ad", "--date=short"),
                   stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
  if (length(tip) == 1 && grepl("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$", tip)) {
    tipd <- as.Date(tip)
    ahead <- unlist(lapply(DOCS, function(f) {
      mk <- suppressWarnings(as.Date(markers(rd(f)))); mk <- mk[!is.na(mk)]
      as.character(mk[mk > tipd]) }))
    chk(sprintf("chronology: no marker postdates the branch tip's author date (%s)", tip),
        length(ahead) == 0, sprintf("(ahead: %s)", paste(unique(ahead), collapse = ", ")))
  } else {
    cat("NOTE  chronology: git unavailable or shallow; the branch-tip check is skipped\n")
  }
})

# ---------------------------------------------------------------------------
# 71. NO IGNORED PATH MAY BE TRACKED (2026-09-13)
#
#     Five .Rproj.user/ files were tracked for the life of this branch even though
#     .gitignore line 2 is ".Rproj.user/". That is not a gitignore bug: .gitignore
#     governs what git STARTS tracking and has no effect on a path already in the
#     index, so adding the rule after the fact changes nothing and the only fix is
#     `git rm --cached`. Nothing warns you, which is why it survived. One git call
#     names every such path, so the invariant is cheap and exact.
# ---------------------------------------------------------------------------
local({
  out <- tryCatch(suppressWarnings(
           system2("git", c("ls-files", "-i", "-c", "--exclude-standard"),
                   stdout = TRUE, stderr = FALSE)), error = function(e) NULL)
  st <- tryCatch(suppressWarnings(
          system2("git", c("rev-parse", "--is-inside-work-tree"),
                  stdout = TRUE, stderr = FALSE)), error = function(e) NULL)
  if (!identical(st, "true")) {
    cat("NOTE  hygiene: not a git work tree; the tracked-but-ignored check is skipped\n")
  } else {
    out <- out[nzchar(out)]
    chk("hygiene: no path that .gitignore excludes is tracked", length(out) == 0,
        if (length(out)) sprintf("(tracked anyway: %s)",
          paste(utils::head(out, 6), collapse = ", ")) else "")
  }
  # the rule that was there all along, plus the duplicate that was removed
  gi <- readLines(".gitignore", warn = FALSE)
  chk("hygiene: .gitignore still excludes the RStudio session directory",
      any(trimws(gi) == ".Rproj.user/"))
  chk("hygiene: the duplicated .Rproj.user notebooks rule is gone",
      sum(grepl("^\\.Rproj\\.user/shared/notebooks/paths$", trimws(gi))) == 0)
  chk("hygiene: .gitignore records why the untracking was needed",
      any(grepl("git rm --cached", gi, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 72. THE MERGE DOCUMENT (2026-09-13; PULL_REQUEST.md until the merge, archived 2026-09-28)
#
#     Matt: "Update the pull request document in preparation for the pull and merge of
#     the two branches."
#
#     The old file described only the first wave of the branch, against a 67,312 total,
#     under a banner telling the reader to go read something else. A PR description that
#     disowns itself is worse than none, because a reviewer opens it first. It is
#     rewritten as the merge document, and because it now carries NUMBERS it can go stale
#     exactly the way the status box did, so it is held to the same standard: the totals
#     it prints are read back against the run folder's own CSV.
#
#     DELIBERATELY NOT ASSERTED: commit SHAs. `git am` recreates commits, so the head and
#     merge-base hashes in the document change the moment Matt applies the series. They
#     are informational; the folder names and totals are the checkable part.
# ---------------------------------------------------------------------------
local({
  # Merged as PR #5 on 2026-09-28 and archived; its totals are still read back against the
  # run folders it names, because an archived record must not misquote them either.
  f <- "07_documentation/archive/PR-5-merge-OSP-boat-count-incorporation.md"
  chk("PR: the document exists", file.exists(f))
  L <- readLines(f, warn = FALSE); t <- paste(L, collapse = "\n")
  flat <- function(x) gsub("[ \n]+", " ", x)
  tf <- flat(t)

  chk("PR: it is a merge document, not a superseded first-wave description",
      grepl("Merge `OSP-boat-count-incorporation` into `main`", tf, fixed = TRUE) &&
      !grepl("this PR description covers only the FIRST wave", tf, fixed = TRUE) &&
      # 67,312 was the OLD file's total. It may be named, but only as the thing that is
      # superseded, never as a current figure.
      (!grepl("67,312", tf, fixed = TRUE) ||
       grepl("against a 67,312 total", tf, fixed = TRUE)) &&
      grepl("It is superseded in full", tf, fixed = TRUE))
  chk("PR: it names Method v2.0 as the method of record and where it is specified",
      grepl("Method v2.0", tf, fixed = TRUE) &&
      grepl("BSS-GH-pooled-CPUE-model-documentation.md", tf, fixed = TRUE))
  chk("PR: it states that no estimate from this pipeline has been published",
      grepl("WDFW has published no estimate from this pipeline", tf, fixed = TRUE))
  chk("PR: no em dash (the stated convention)", !grepl("—", t) && !grepl("–", t))

  # (1) the two runs it compares must both exist, and their totals must match the file
  runs <- unique(regmatches(t, gregexpr("05_output/[0-9]{8}/[A-Za-z0-9._-]+", t))[[1]])
  chk("PR: it names both the branch's authoritative run and main's, and both exist",
      length(runs) >= 2 && all(dir.exists(runs)) &&
      any(grepl("20260927/pooled-CPUE-canonical-2024-25", runs)) &&
      any(grepl("20260715/pooled-CPUE-230256", runs)),
      sprintf("(named: %s)", paste(runs, collapse = ", ")))
  # A total and its interval must appear TOGETHER ON ONE LINE, in median / lo / hi order.
  # Searching the whole flattened document for each number separately is not a check: the
  # same figures recur in several tables here, so corrupting one occurrence leaves the
  # others and the assertion passes. Measured: a deliberately wrong port median was NOT
  # caught by the any-occurrence form. Co-occurrence on a line is what makes it fire.
  fm <- function(x) formatC(round(as.numeric(x)), format = "d", big.mark = ",")
  triple_on_a_line <- function(med, lo, hi) {
    pat <- paste0(gsub(",", "\\,", fm(med)), ".*", gsub(",", "\\,", fm(lo)),
                  ".*", gsub(",", "\\,", fm(hi)))
    any(grepl(pat, L, fixed = FALSE))
  }
  agrees <- function(folder) {
    pt <- file.path(folder, "port_total_Dungeness_Kept.csv")
    if (!file.exists(pt)) return(FALSE)
    d <- utils::read.csv(pt, stringsAsFactors = FALSE)
    # the pooled track labels this row "Expected_Catch", the gear track "Catch" (the
    # label mismatch recorded as the 2026-09-08 defect fix in run_adoption_2026-09-07.R)
    r <- d[grepl("^(Expected_)?Catch$", d[[2]]), , drop = FALSE]
    if (nrow(r) != 1) return(FALSE)
    triple_on_a_line(r$BSS_median, r$BSS_lo95, r$BSS_hi95)
  }
  for (r in runs[dir.exists(runs)])
    chk(sprintf("PR: the totals it prints for %s match that folder's own CSV", basename(r)),
        agrees(r), "(port_total_Dungeness_Kept.csv vs the archived merge document)")

  # the headline comparison ROW must carry BOTH triples, because the same figures recur
  # elsewhere in the file and an any-line check would be satisfied by the other table
  total_of <- function(folder) {
    pt <- file.path(folder, "port_total_Dungeness_Kept.csv")
    if (!file.exists(pt)) return(NULL)
    d <- utils::read.csv(pt, stringsAsFactors = FALSE)
    r <- d[grepl("^(Expected_)?Catch$", d[[2]]), , drop = FALSE]
    if (nrow(r) != 1) return(NULL)
    list(med = as.numeric(r$BSS_median), lo = as.numeric(r$BSS_lo95),
         hi = as.numeric(r$BSS_hi95), pe = as.numeric(r$PE))
  }
  A <- total_of("05_output/20260927/pooled-CPUE-canonical-2024-25")       # this branch (since 2026-09-28)
  A4 <- total_of("05_output/20260910/pooled-CPUE-IMP-R4-shore-tau-newf")  # R4, the configuration R5 was rendered at
  M <- total_of("05_output/20260715/pooled-CPUE-230256")                  # main
  G <- total_of("05_output/20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf")
  row <- grep("\\*\\*Port total\\*\\*", L, value = TRUE)
  seq_on <- function(line, ...) {
    nums <- vapply(list(...), function(v) gsub(",", "\\,", fm(v)), character(1))
    grepl(paste(nums, collapse = ".*"), line)
  }
  chk("PR: the headline comparison row carries BOTH runs' totals and intervals, in order",
      length(row) == 1 && !is.null(A) && !is.null(M) &&
      seq_on(row[1], M$med, M$lo, M$hi, A$med, A$lo, A$hi),
      "(the **Port total** row of section 1)")

  # THE ARITHMETIC IS THE REAL GUARD. A corrupted figure survives a string search (the
  # any-occurrence form was measured NOT to catch a wrong port median) but it cannot
  # survive the subtraction, so the stated deltas are recomputed from the two CSVs.
  chk("PR: the stated change from main is the arithmetic of the two runs",
      !is.null(A) && !is.null(M) &&
      grepl(sprintf("+%s (+%.1f%%)", fm(A$med - M$med), 100 * (A$med - M$med) / M$med),
            tf, fixed = TRUE),
      sprintf("(expected +%s (+%.1f%%))", fm(A$med - M$med),
              100 * (A$med - M$med) / M$med))
  chk("PR: the stated PE change from main is the arithmetic of the two runs",
      !is.null(A) && !is.null(M) &&
      grepl(sprintf("+%s", fm(A$pe - M$pe)), tf, fixed = TRUE),
      sprintf("(expected +%s)", fm(A$pe - M$pe)))
  # the gear cross-check is like for like only at the configuration it was rendered at (R4's); under
  # the method of record it is owed, and the document must say so rather than compare R5 with 96,118
  chk("PR: the gear cross-check gap is the arithmetic of the two runs",
      !is.null(A4) && !is.null(G) &&
      triple_on_a_line(G$med, G$lo, G$hi) &&
      grepl(sprintf("%.2f%% below", abs(100 * (G$med - A4$med) / A4$med)), tf, fixed = TRUE) &&
      grepl("The same check under the method of record is owed", tf, fixed = TRUE),
      sprintf("(gear %s vs pooled R4 %s = %.2f%%)", fm(G$med), fm(A4$med),
              100 * (G$med - A4$med) / A4$med))

  # (3) the ladder table must be the ladder CSV, rung for rung
  lf <- "05_output/improvements_2026-09-08_ladder-newf.csv"
  chk("PR: every rung's total AND interval come from the ladder CSV, on one line each",
      { if (!file.exists(lf)) NA else {
          d <- utils::read.csv(lf, stringsAsFactors = FALSE)
          all(vapply(seq_len(nrow(d)), function(i)
            triple_on_a_line(d$port[i], d$port_lo95[i], d$port_hi95[i]), logical(1))) } },
      sprintf("(%s)", lf))

  # (5) the open decisions must be the ones the register still has OPEN or BLOCKED
  cr <- paste(readLines("07_documentation/development_notes/CHANGE_REGISTER.md",
                        warn = FALSE), collapse = "\n")
  for (d in c("D14", "D2", "D24", "D19", "D8", "D18", "D3", "D6"))
    chk(sprintf("PR: %s is carried forward and the register still knows it", d),
        grepl(sprintf("\\*\\*[^*]*\\b%s\\b[^*]*\\*\\*", d), tf) &&
        grepl(sprintf("\\| %s \\|", d), cr))
  chk("PR: it does not claim the hygiene decision was made for Matt",
      grepl("a decision, not a task", tf, fixed = TRUE) &&
      grepl("Tier 4", tf, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 74. D3 / D6: THE GEAR AR PERIOD, THE GEAR ZINB, AND THE RUN THAT SETTLES THEM
#     (2026-09-13)
#
#     Matt: "Write a driver for the D3/D6 gear track cap and ZI block (~3 h ladder) runs
#     within the diagnostics folder. I want to hit run on it and have it generate all of
#     the information to finalize D3/D6."
#
#     Three things have to hold for that to be true, and each is asserted below rather
#     than described. (a) The D6 Stan port must be INERT until asked for, or landing it
#     silently moves the committed R5 cross-check figure. (b) The D3 lever must actually
#     reach the sampler, because the register named a DORMANT key for weeks. (c) The
#     driver's decision rules must be in the file, before the run, including the two
#     clauses that exclude evidence on purpose.
# ---------------------------------------------------------------------------
local({
  rd  <- function(f) readLines(f, warn = FALSE)
  flat <- function(x) gsub("[ \n]+", " ", paste(x, collapse = "\n"))
  GS <- "02_stan_models/crab_bss_gear_resolved.stan"
  PS <- "02_stan_models/crab_bss_pooled.stan"
  GP <- "03_R_functions/prep_bss_crab_gear.R"
  DR <- "06_diagnostics/run_gear_ar_zi_2026-09-13.R"
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config

  # (a) THE D6 PORT IS SYMMETRIC WITH THE POOLED MODEL AND INERT AS SHIPPED
  g <- rd(GS); pl <- rd(PS)
  chk("D6: the gear Stan declares the three ZI data variables",
      all(vapply(c("zi_catch", "zi_catch_prior_a", "zi_catch_prior_b"),
                 function(v) any(grepl(paste0("(^|[^A-Za-z0-9_])", v, "\\s*;"), g)), logical(1))))
  chk("D6: theta_C is zero-size when off, in BOTH models, spelled identically",
      sum(grepl("vector<lower=0, upper=1>[zi_catch] theta_C;", g, fixed = TRUE)) == 1 &&
      sum(grepl("vector<lower=0, upper=1>[zi_catch] theta_C;", pl, fixed = TRUE)) == 1,
      "(this declaration is the whole basis of the OFF path being bit-identical)")
  chk("D6: the gear model carries the log_mix branch and the beta prior, like the pooled one",
      any(grepl("log_mix(theta_C[1], 0, neg_binomial_2_lpmf(0 | mu_c, r_C))", g, fixed = TRUE)) &&
      any(grepl("theta_C[1] ~ beta(zi_catch_prior_a, zi_catch_prior_b)", g, fixed = TRUE)))
  chk("D6: the gear season total is scaled by zi_scale exactly once",
      sum(grepl("* fd * zi_scale;", g, fixed = TRUE)) == 1,   # D40: fd is the day's share
      paste("Without it, turning ZI on inflates the total by 1/(1 - theta_C) as an artefact,",
            "because lambda_C rises to absorb the zeros theta_C removed."))
  chk("D6: theta_C_out and zi_scale are declared unconditionally, so every run has the columns",
      any(grepl("^\\s*real theta_C_out;", g)) && any(grepl("^\\s*real<lower=0, upper=1> zi_scale;", g)))
  chk("D6: the gear prep builds the ZI variables and gates them on catch_zi_tracks",
      { t <- flat(rd(GP))
        grepl("zi_catch = as.integer(isTRUE(params$estimate_catch_zi) && \"gear_resolved\" %in%", t, fixed = TRUE) })
  chk("D6: the POOLED prep is gated on the same key, so it means one thing on both tracks",
      { t <- flat(rd("03_R_functions/prep_bss_crab_pooled.R"))
        grepl("\"pooled\" %in% (params$catch_zi_tracks %||% \"pooled\")", t, fixed = TRUE) })
  # 2026-10-02 (A32): D6 was decided and adopted, together with D3. This pinned the OFF state
  # until a decision; it now pins the decided one, so a revert is as visible as the flip was.
  chk("D6 SHIPS ON FOR BOTH TRACKS since A32 (2026-10-02): catch_zi_tracks is c(\"pooled\", \"gear_resolved\")",
      identical(as.character(rc$catch_zi_tracks), c("pooled", "gear_resolved")),
      paste("A32 adopted it with D3; the gear cross-check is the batch's stage D6 render, and a revert",
            "would make that folder stop being the shipped file's gear render with nobody touching a number."))
  # the contract that caught the 2026-08-25 disaster must still hold on the widened model
  chk("D6: every variable the gear Stan now declares is still built in its prep",
      { need <- bss_stan_data_names(GS); src <- paste(rd(GP), collapse = "\n")
        !length(need[!vapply(need, function(v)
          grepl(paste0("(^|[^A-Za-z0-9_.])", v, "([^A-Za-z0-9_]|$)"), src, perl = TRUE), logical(1))]) })

  # (b) THE D3 LEVER REACHES THE SAMPLER, AND SHIPS UNCHANGED
  # 2026-10-02 (A32): the per-population form ships; the flat form's literals are still the
  # sub-season default build_subseasons() falls back to (the next check), so both are pinned.
  source("03_R_functions/bss_gear_period.R", local = TRUE)
  chk("D3: gear_period_bss ships in the per-population form adopted by A32 (shore at the pooled caps, boat as before)",
      bss_gear_period_is_per_pop(rc) &&
      identical(bss_gear_period(rc, "shore", "all_gear"), "weekly") && identical(bss_gear_period(rc, "shore", "pot_closure"), "biweekly") &&
      identical(bss_gear_period(rc, "private_boat", "all_gear"), "month") && identical(bss_gear_period(rc, "private_boat", "pot_closure"), "biweekly") &&
      identical(bss_gear_period(rc, "shore", "all_gear"), as.character(rc$ar_max_resolution$pooled$shore$all_gear)),
      "(the boat values are the literals build_subseasons.R used before 2026-09-13)")
  chk("D3: build_subseasons has no period_bss literal left",
      !any(grepl('period_bss = "', rd("03_R_functions/build_subseasons.R"), fixed = TRUE)))
  chk("D3: the key DRIVES build_subseasons, at the shipped value and at a changed one",
      { a <- build_subseasons(rc)
        b <- build_subseasons(modifyList(rc, list(
               gear_period_bss = list(all_gear = "weekly", pot_closure = "biweekly"))))
        ga <- Filter(function(x) identical(x$gear_regime, "all_gear"), a)[[1]]
        gb <- Filter(function(x) identical(x$gear_regime, "all_gear"), b)[[1]]
        ca <- Filter(function(x) identical(x$gear_regime, "pot_closure"), a)[[1]]
        cb <- Filter(function(x) identical(x$gear_regime, "pot_closure"), b)[[1]]
        identical(ga$period_bss, "month") && identical(gb$period_bss, "weekly") &&
        identical(ca$period_bss, "biweekly") && identical(cb$period_bss, "biweekly") },
      "(and the pot-closure sub-season does not move with it)")
  chk("D3: the dormant cap is still documented as dormant, so nobody edits it expecting an effect",
      { t <- flat(rd("run_config.R"))
        grepl("ar_max_resolution$gear_resolved is DORMANT", t, fixed = TRUE) ||
        grepl("ar_max_resolution$gear_resolved) is DORMANT", t, fixed = TRUE) ||
        grepl("is DORMANT", t, fixed = TRUE) })
  # the G = 5 claim, corrected. This is the one that would quietly come back.
  chk("D3: gear_resolved_G ships FALSE, so G is 1 and the 'G = 5' caution does not apply",
      identical(rc$gear_resolved_G, FALSE),
      paste("Measured 2026-09-13 on the real 2024-25 shore all-gear data: G = 1. The",
            "register and the adoption review both justified NOT copying the pooled period",
            "across by citing a per-gear likelihood at G = 5 that the shipped configuration",
            "does not produce."))
  chk("D3: the corrected G finding is recorded where the claim was made",
      { r <- flat(rd("07_documentation/development_notes/CHANGE_REGISTER.md"))
        p2 <- flat(rd("07_documentation/development_notes/PIPELINE_STATUS.md"))
        grepl("G = 1", r, fixed = TRUE) && grepl("G = 1", p2, fixed = TRUE) &&
        grepl("gear_resolved_G", r, fixed = TRUE) })

  # (c) THE DRIVER: ships safe, states its rules first, pins its window
  chk("D3/D6: the driver exists and ships DRY_RUN <- TRUE", file.exists(DR) &&
      any(grepl("^DRY_RUN <- TRUE", rd(DR))))
  d <- rd(DR); dt <- flat(d)
  chk("D3/D6: the driver declares all six stages",
      all(vapply(c("G0", "G1", "G2", "G3", "G4", "G5"),
                 function(x) grepl(paste0("\n  ", x, " = list\\(tag"), paste(d, collapse = "\n")), logical(1))))
  chk("D3/D6: G1 is the shipped period AND the bit-identity control",
      grepl('G1 = list(tag = "GZ-G1-month"', dt, fixed = TRUE) &&
      grepl("bit-identity", dt, fixed = TRUE) &&
      grepl("20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf", dt, fixed = TRUE))
  chk("D3/D6: both decision rules are in the file, numbered, BEFORE the code",
      { i1 <- regexpr("THE DECISION RULE FOR D3", dt, fixed = TRUE)
        i2 <- regexpr("THE DECISION RULE FOR D6", dt, fixed = TRUE)
        i3 <- regexpr("DRY_RUN <- TRUE", dt, fixed = TRUE)
        i1 > 0 && i2 > 0 && i1 < i3 && i2 < i3 },
      "(a rule written after the numbers are in is a rationalization)")
  chk("D3/D6: the two EXCLUSION clauses are explicit (elpd, and the cross-track gap)",
      grepl("DO NOT SELECT ON elpd", dt, fixed = TRUE) &&
      grepl("THE CROSS-TRACK GAP IS NOT A CRITERION", dt, fixed = TRUE) &&
      grepl("THE PORT TOTAL IS NOT A CRITERION", dt, fixed = TRUE),
      paste("On the pooled ladder elpd favoured the daily fit every other diagnostic called",
            "overfitted; and choosing this track's period to minimize the cross-track gap",
            "would tune one estimate to another instead of to the data."))
  chk("D3/D6: the driver pins a WINDOW and knows its delta keys",
      grepl("WINDOW <- list(", dt, fixed = TRUE) &&
      grepl("DELTA_KEYS <- unique(c(\"ar_force\", \"catch_zi_tracks\"))", dt, fixed = TRUE))
  # 2026-09-14: the lever is ar_force on the SHORE ONLY, and the sub-season periods are
  # pinned instead. The first run used gear_period_bss, which is a SUB-SEASON key, so it
  # moved the boat all-gear fit as well and 90% of the port movement came from a component
  # D3 does not ask about. This asserts the corrected shape: ar_adaptive and ar_escalate
  # pinned so nothing else can select a resolution, the sub-season periods pinned at the
  # shipped values so the BOAT stays at monthly (which is what the pooled track fits it
  # at, making each rung's port comparable to the pooled R4), and ar_force NOT pinned,
  # because it is the lever.
  chk("D3/D6: the pin fixes the selection mode and the SUB-SEASON periods, and leaves ar_force free",
      all(vapply(c("ar_adaptive = FALSE", "ar_escalate = FALSE",
                   'gear_period_bss = list(all_gear = "month", pot_closure = "biweekly")'),
                 function(x) grepl(x, dt, fixed = TRUE), logical(1))) &&
      !grepl("ar_force = NULL,", dt, fixed = TRUE))
  chk("D3/D6: the rung lever moves the SHORE only, so the boat cannot ride along",
      grepl(".shore_at <- function(res) list(ar_force = list(shore = list(all_gear = res)))",
            dt, fixed = TRUE) &&
      !grepl(".pb <- function(ag)", dt, fixed = TRUE),
      paste("Same class of defect as the 2026-08-27 Stage C ar_force bug, which forced both",
            "boat sub-seasons to biweekly and made its port total uninterpretable as the",
            "change it was meant to isolate."))
  chk("D3/D6: the verdicts read the fields the helpers actually return",
      grepl("el$elpd_diff", dt, fixed = TRUE) && grepl("el$se_diff", dt, fixed = TRUE) &&
      grepl("loo_elpd_paired_str(el)", dt, fixed = TRUE) &&
      grepl("loo_elpd_by_count_str(el)", dt, fixed = TRUE) &&
      grepl('identical(ex$verdict, "PASS")', dt, fixed = TRUE) &&
      # the CODE patterns, not the strings: the comment above the fix names the old
      # fields on purpose, and an assertion that forbade mentioning them would forbid
      # recording why they were wrong
      !grepl("isTRUE(ex$identical", dt, fixed = TRUE) &&
      !grepl("gain <- el$diff", dt, fixed = TRUE),
      paste("The first run read el$diff / el$se / ex$identical, none of which exist.",
            "loo_elpd_paired() returns elpd_diff / se_diff / zeros / positives / by_count and",
            "ships two renderers; fit_exactness() returns observed / verdict. The wrong names",
            "turned a PASS into a FAIL and a +11.3-nat gain into NA, and the recommendation",
            "then said 'do not adopt' on evidence that passes."))
  chk("D3/D6: the desk stage asserts EXACTLY ONE FIT MOVES, per population and sub-season",
      grepl("EXACTLY ONE FIT MOVES", dt, fixed = TRUE) &&
      grepl(".bss_resolve_ar_force(cfg, pop, x$gear_regime)", dt, fixed = TRUE) &&
      grepl('"private_boat/all_gear" = "monthly"', dt, fixed = TRUE),
      paste("The check it replaces asked only whether the requested period reached",
            "build_subseasons() and whether the POT CLOSURE moved. It never asked which",
            "POPULATIONS moved, so the lever that moved shore and boat together passed it",
            "and the ladder spent 3 h measuring the wrong thing. Negative-tested: restoring",
            "the sub-season lever fails this row on all four rungs."))
  chk("D3/D6: a clause that cannot be EVALUATED is REVIEW, never FAIL",
      grepl("NOT COMPUTABLE", dt, fixed = TRUE) && grepl("any_rev", dt, fixed = TRUE) &&
      grepl("could not be EVALUATED", dt, fixed = TRUE),
      paste("The first run scored an un-computable clause as a failure, which is how a",
            "missing statistic became a verdict against the feature."))
  chk("D3/D6: the driver verifies the resolution the sampler ACTUALLY used, not the one requested",
      grepl("ar_escalation_log.csv", dt, fixed = TRUE) &&
      grepl("the fit used the resolution the rung asked for", dt, fixed = TRUE),
      paste("A rung that fell back would otherwise report itself as the resolution it asked",
            "for and sit in the ladder as a duplicate of another rung."))
  chk("D3/D6: the runtime estimate is derived from measured timings, not asserted",
      grepl("run_timings.csv", dt, fixed = TRUE) && grepl("10.0 min", dt, fixed = TRUE) &&
      grepl("sub-linear", dt, fixed = TRUE))
  chk("D3/D6: every rung is wrapped, so one failed render does not cost the others",
      length(grep("tryCatch(", d, fixed = TRUE)) >= 5)

  # the census-frame disclosure Matt asked for
  chk("census: estimate_comm_charter returns frame_warnings on BOTH return paths",
      { t <- rd("03_R_functions/estimate_comm_charter.R")
        length(grep("frame_warnings = .fw", t, fixed = TRUE)) >= 2 })
  chk("census: the multi-season merge CONCATENATES the warnings instead of dropping them",
      { t <- flat(rd("03_R_functions/estimate_comm_charter.R"))
        grepl('warn_keys <- c("frame_warnings")', t, fixed = TRUE) &&
        grepl("for (k in warn_keys) if (length(per[[i]][[k]]))", t, fixed = TRUE) })
  for (f in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd",
              "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"))
    chk(sprintf("census: %s prints the frame conditions, logs them for the warnings section and writes the CSV", basename(f)),
        { t <- flat(rd(f))
          # 2026-09-29: the pooled report prints them in Section 4.2 (its port-total chunk is hidden)
          (grepl("CENSUS FRAME:", t, fixed = TRUE) || grepl("**Census frame conditions (%d).**", t, fixed = TRUE)) &&
          grepl('bss_warn("census", .w', t, fixed = TRUE) &&
          grepl("census_frame_warnings.csv", t, fixed = TRUE) },
        paste("A bare warning() does not reach a rendered report: knitr defers warnings,",
              "html_document can hide them, and quiet = TRUE suppresses them entirely."))
  chk("census: the five frame conditions are all reported",
      { t <- flat(rd("03_R_functions/estimate_comm_charter.R"))
        all(vapply(c("NO CENSUS FRAME", "ROSTER-ONLY DAYS", "FRAMES DISAGREE",
                     "THIN CHARTER SAMPLE", "UNSAMPLED DAYS TREATED AS NO OPERATION"),
                   function(x) grepl(x, t, fixed = TRUE), logical(1))) })
  # and the harness must REPORT that rather than dying on it. Negative-tested: with the
  # warning turned into a stop(), chk() now records
  # "FAIL census: empty window returns the zero split [ERROR while evaluating: ...]"
  # where it previously aborted the run 300 assertions early and printed no summary at all.
  chk("harness: an error while evaluating a condition is a FAIL, not the end of the run",
      { r <- tryCatch({ o0 <- ok; b0 <- bad
                        chk("^^ IGNORE: this FAIL line is the self-test below, deliberately erroring", stop("deliberate"))
                        bad == b0 + 1 && ok == o0 }, error = function(e) FALSE)
        # undo the self-test's own bookkeeping so it does not show in the totals twice
        if (isTRUE(r)) { bad <<- bad - 1 }
        isTRUE(r) },
      "(a thousand-assertion harness that aborts on one surprise is the least useful outcome)")
  chk("census: it WARNS rather than stopping, which is the decision Matt made",
      { t <- flat(rd("03_R_functions/estimate_comm_charter.R"))
        grepl("warning(paste0(\"estimate_comm_charter(): \", .m), call. = FALSE)", t, fixed = TRUE) &&
        !grepl("stop(paste0(\"estimate_comm_charter()", t, fixed = TRUE) })

  # the superseded-runner guard
  SUP <- c("run_patch_validation_2026-08-25", "run_improvement_plan_2026-08-27",
           "run_stage5_2026-08-30", "run_validation_2026-09-01",
           "run_shore_ar_zi_2026-09-03", "run_ladder_zinb_2026-09-04",
           "run_adoption_2026-09-07", "run_osp_validation",
           "run_tau_sweep",   # 2026-09-28 (B46): superseded by its own header; now guarded
           "run_gear_ar_zi_2026-09-13",           # 2026-10-02 (B62): D3 and D6 settled by the batch, adopted (A32)
           "run_authoritative_batch_2026-09-29")  # 2026-10-02 (B62): finished; its base, run_config.R at c4a0241, moved
  LIVE <- c("run_improvements_2026-09-08",   # run_gear_ar_zi_2026-09-13 moved to SUP on 2026-10-02 (B62)
            "run_rg_sweep",
            "run_marine_hazard_batch_2026-09-25",   # 2026-09-25: A30 / B35
            "run_marine_block_cv_2026-09-26")       # D37 (2026-09-29): B39's runner was in neither list
  chk("diagnostics: the superseded-runner helper exists and offers an override",
      file.exists("03_R_functions/bss_superseded_runner.R") &&
      { t <- flat(rd("03_R_functions/bss_superseded_runner.R"))
        grepl("I_KNOW_THIS_IS_SUPERSEDED", t, fixed = TRUE) &&
        grepl("bss_superseded_runner <- function", t, fixed = TRUE) })
  chk("diagnostics: every SETTLED runner refuses to fit",
      { miss <- SUP[!vapply(SUP, function(f)
          any(grepl("bss_superseded_runner(", rd(file.path("06_diagnostics", paste0(f, ".R"))),
                    fixed = TRUE)), logical(1))]
        length(miss) == 0 },
      "(a runner whose question is closed but which still fits is the defect this closes)")
  chk("diagnostics: no LIVE runner is guarded",
      { bad <- LIVE[vapply(LIVE, function(f)
          any(grepl("bss_superseded_runner(", rd(file.path("06_diagnostics", paste0(f, ".R"))),
                    fixed = TRUE)), logical(1))]
        length(bad) == 0 },
      "(guarding the current ladder would be the same mistake in reverse)")
  chk("diagnostics: every run_*.R is classified SETTLED or LIVE here (D37: run_marine_block_cv was in neither)",
      { allr <- sub("[.]R$", "", list.files("06_diagnostics", pattern = "^run_.*[.]R$"))
        length(setdiff(allr, c(SUP, LIVE))) == 0 && length(intersect(SUP, LIVE)) == 0 },
      sprintf("(unclassified: %s)", paste(setdiff(sub("[.]R$", "", list.files("06_diagnostics", pattern = "^run_.*[.]R$")), c(SUP, LIVE)), collapse = ", ")))
  chk("diagnostics: the README labels every runner LIVE or RECORD",
      { t <- rd("06_diagnostics/README.md")
        rows <- grep("^\\| `run_|^\\| `test_|^\\| `gear_coverage", t, value = TRUE)
        length(rows) >= 14 && all(grepl("\\*\\*(LIVE|RECORD)\\.\\*\\*", rows)) },
      "(so a reader can tell which is which without reading all eleven)")
  chk("diagnostics: the README states the measured lever counts rather than an impression",
      { t <- flat(rd("06_diagnostics/README.md"))
        grepl("pins 12 of 15", t, fixed = TRUE) && grepl("pins 0 to 5", t, fixed = TRUE) })
  chk("diagnostics: the deletion option is recorded WITH its cost, not silently taken",
      { t <- flat(rd("07_documentation/development_notes/PIPELINE_STATUS.md"))
        grepl("the answer is NOT deletion", t, fixed = TRUE) &&
        grepl("no `run_parameters.txt` at all", t, fixed = TRUE) &&
        grepl("git rm 06_diagnostics/run_", t, fixed = TRUE) })
})

# ---------------------------------------------------------------------------
# 75. THE GEAR AR / ZINB RUN: ITS RESULTS, AND THE FOUR DEFECTS IT EXPOSED (2026-09-14)
#
#     The run completed and reported "D3 weekly, D6 do not adopt". Both were wrong, in the
#     driver rather than the data, and the fixes are pinned in section 74. What is pinned
#     HERE is the EVIDENCE, because every number the register and the status document now
#     rest on was read out of committed CSVs and can be re-read. A conclusion drawn from a
#     folder that no longer says what the document claims is worse than no conclusion.
#
#     THE CONFOUND IS THE THING TO REMEMBER. gear_period_bss$all_gear moved the shore AND
#     the boat all-gear fits, because period_bss is a sub-season key. 90% of the port
#     movement was the boat. This asserts the confound is real in the committed folders,
#     so nobody re-reads that ladder as if it had isolated the shore.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  R <- c(G1 = "05_output/20260913/gear-type-CPUE-model-GZ-G1-month",
         G2 = "05_output/20260914/gear-type-CPUE-model-GZ-G2-biweekly",
         G3 = "05_output/20260914/gear-type-CPUE-model-GZ-G3-weekly",
         G4 = "05_output/20260914/gear-type-CPUE-model-GZ-G4-daily",
         G5 = "05_output/20260914/gear-type-CPUE-model-GZ-G5-zi")
  R5 <- "05_output/20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf"
  P4 <- "05_output/20260910/pooled-CPUE-IMP-R4-shore-tau-newf"
  have <- all(dir.exists(c(R, R5, P4)))
  chk("1w: all five rungs and both references are committed", have,
      sprintf("(missing: %s)", paste(basename(c(R, R5, P4))[!dir.exists(c(R, R5, P4))], collapse = ", ")))
  if (!have) { cat("NOTE  1w: rung folders absent; the evidence assertions are skipped\n") } else {

  comp <- function(d, key) { x <- rd(file.path(d, "pe_vs_bss_comparison.csv"))
    if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x$BSS_catch[x$component == key])[1]) }
  ar <- function(d, fit) { x <- rd(file.path(d, "ar_escalation_log.csv"))
    if (is.null(x)) NA_character_ else as.character(x$ar_resolution[grepl(fit, x$fit)])[1] }
  adq <- function(d, fit, col) { x <- rd(file.path(d, "model_adequacy.csv"))
    if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x[[col]][x$fit == fit])[1]) }

  # (1) THE CONFOUND. Both all-gear fits moved together in every rung.
  chk("1w: the ladder moved BOTH all-gear fits, which is why it did not isolate D3",
      all(vapply(names(R), function(k)
        identical(ar(R[[k]], "shore_all_gear"), ar(R[[k]], "private_boat_all_gear")), logical(1))),
      paste("period_bss is a SUB-SEASON key and the all_gear sub-season holds both",
            "populations. If this ever comes back FALSE the folders have changed, not the",
            "finding."))
  sh <- vapply(c("G1","G3"), function(k) comp(R[[k]], "shore (All gear)"), numeric(1))
  bo <- vapply(c("G1","G3"), function(k) comp(R[[k]], "private_boat (All gear)"), numeric(1))
  chk("1w: 90% of the port movement was the BOAT, not the shore",
      isTRUE(all(is.finite(c(sh, bo)))) &&
      isTRUE(abs((bo[2] - bo[1]) / (bo[2] - bo[1] + sh[2] - sh[1]) - 0.90) < 0.03),
      sprintf("(shore %+.0f, boat %+.0f; boat share %.0f%%)",
              sh[2] - sh[1], bo[2] - bo[1],
              100 * (bo[2] - bo[1]) / (bo[2] - bo[1] + sh[2] - sh[1])))

  # (2) WHAT THE RUN DID SETTLE: the shore takes weekly and not daily.
  for (k in c("G1","G2","G3"))
    chk(sprintf("1w: the shore fit is ADEQUATE at %s", ar(R[[k]], "shore_all_gear")),
        isTRUE(adq(R[[k]], "shore_all_gear_Dungeness_Kept", "p_loo_frac") <= 0.15) &&
        isTRUE(adq(R[[k]], "shore_all_gear_Dungeness_Kept", "n_pareto_bad") == 0),
        sprintf("(p_loo %.4f, %.0f bad k)",
                adq(R[[k]], "shore_all_gear_Dungeness_Kept", "p_loo_frac"),
                adq(R[[k]], "shore_all_gear_Dungeness_Kept", "n_pareto_bad")))
  chk("1w: the shore fit FAILS at daily, replicating the pooled daily rejection",
      isTRUE(adq(R[["G4"]], "shore_all_gear_Dungeness_Kept", "p_loo_frac") > 0.15) &&
      isTRUE(adq(R[["G4"]], "shore_all_gear_Dungeness_Kept", "n_pareto_bad") >= 30),
      sprintf("(p_loo %.4f, %.0f bad k; the pooled daily fit was 0.352 with 41)",
              adq(R[["G4"]], "shore_all_gear_Dungeness_Kept", "p_loo_frac"),
              adq(R[["G4"]], "shore_all_gear_Dungeness_Kept", "n_pareto_bad")))

  # (3) D29: the boat is adequate at three periods and moves ~25% across them.
  chk("1w / D29: the BOAT fit is adequate at monthly, biweekly AND weekly, so adequacy does not settle its period",
      all(vapply(c("G1","G2","G3"), function(k)
        isTRUE(adq(R[[k]], "private_boat_all_gear_Dungeness_Kept", "p_loo_frac") <= 0.15), logical(1))) &&
      isTRUE(adq(R[["G4"]], "private_boat_all_gear_Dungeness_Kept", "p_loo_frac") > 0.15),
      sprintf("(p_loo monthly %.4f, biweekly %.4f, weekly %.4f, daily %.4f)",
              adq(R[["G1"]], "private_boat_all_gear_Dungeness_Kept", "p_loo_frac"),
              adq(R[["G2"]], "private_boat_all_gear_Dungeness_Kept", "p_loo_frac"),
              adq(R[["G3"]], "private_boat_all_gear_Dungeness_Kept", "p_loo_frac"),
              adq(R[["G4"]], "private_boat_all_gear_Dungeness_Kept", "p_loo_frac")))
  pb <- comp(P4, "private_boat (All gear)")
  spread <- range(vapply(names(R)[1:4], function(k) comp(R[[k]], "private_boat (All gear)"), numeric(1)))
  chk("1w / D29: the boat component spans more than 20% across the four periods",
      isTRUE(is.finite(pb)) && isTRUE((spread[2] - spread[1]) / pb > 0.20),
      sprintf("(%s to %s, %.1f%% of the pooled monthly fit %s)",
              format(round(spread[1]), big.mark = ","), format(round(spread[2]), big.mark = ","),
              100 * (spread[2] - spread[1]) / pb, format(round(pb), big.mark = ",")))

  # (4) THE D6 STAN PORT IS INERT, at production iterations, in the committed folders.
  ex <- tryCatch(fit_exactness(R[["G1"]], R5, what = "the gear fits",
          expect_delta = c("catch_zi_tracks","gear_period_bss","run_weather","ar_force",
                           "pot_closures","census_windows","tau_shore_derive_window_only")),
        error = function(e) NULL)
  chk("1w: the D6 Stan port is INERT when off, proven against the committed R5 render",
      !is.null(ex) && identical(ex$verdict, "PASS"),
      if (is.null(ex)) "(could not compare)" else substr(ex$observed, 1, 150))
  cmpk <- c("shore (Pot closure)","shore (All gear)","private_boat (Pot closure)","private_boat (All gear)")
  chk("1w: and every BSS component of G1 equals R5 to the crab",
      identical(vapply(cmpk, function(k) comp(R[["G1"]], k), numeric(1)),
                vapply(cmpk, function(k) comp(R5, k), numeric(1))),
      paste("The port totals differ by four crab because the driver adds the census as a",
            "random draw. The components are the deterministic quantity."))

  # (5) D6 IS ADOPT. The statistic the driver failed to read, recomputed from the folders.
  el <- tryCatch(loo_elpd_paired(
          file.path(R[["G3"]], "loo_pointwise_catch_shore_all_gear_Dungeness_Kept.csv"),
          file.path(R[["G5"]], "loo_pointwise_catch_shore_all_gear_Dungeness_Kept.csv")),
        error = function(e) NULL)
  chk("1w / D6: the paired elpd is COMPUTABLE from the committed folders (the driver reported NA)",
      !is.null(el) && isTRUE(is.finite(el$elpd_diff)) && isTRUE(is.finite(el$se_diff)),
      if (is.null(el)) "(NULL)" else sprintf("(%+.1f nats at %.1f paired SE)", el$elpd_diff, el$se_diff))
  chk("1w / D6: the gain exceeds 2 paired SE, so clause 2 PASSES",
      !is.null(el) && isTRUE(el$ratio > 2),
      if (is.null(el)) "" else sprintf("(%.2f SE; the naive SE would have read %.1f)", el$ratio, el$se_naive_b))
  chk("1w / D6: the positive-count loss is smaller than the zero-count gain, so clause 4 PASSES",
      !is.null(el) && isTRUE(abs(el$positives[["diff"]]) < el$zeros[["diff"]]),
      if (is.null(el)) "" else sprintf("(zeros %+.1f, positives %+.1f)",
                                       el$zeros[["diff"]], el$positives[["diff"]]))
  chk("1w / D6: the bins carrying most of the catch all GAIN, which the two-way split hides",
      !is.null(el) && { b <- el$by_count
        rs <- rownames(b) %in% c("3-4","5-8","9-16")
        all(b[rs, "diff"] > 0) && sum(b[rs, "catch_share"]) > 0.7 },
      if (is.null(el)) "" else loo_elpd_by_count_str(el))
  zi <- function(d, col) { x <- rd(file.path(d, "ppc_byobs_shore_all_gear_Dungeness_Kept.csv"))
    if (is.null(x)) NA_real_ else { y <- x[x$data_type == "catch", , drop = FALSE]
      p <- suppressWarnings(as.numeric(y[[col]])); ok <- is.finite(p); p <- p[ok]
      k <- if (col == "p_zero") 0L else 1L
      (sum(y$observed[ok] == k) - sum(p)) / sqrt(sum(p * (1 - p))) } }
  chk("1w / D6: BOTH count bins improve, so clause 3 PASSES",
      isTRUE(abs(zi(R[["G5"]], "p_zero")) < abs(zi(R[["G3"]], "p_zero"))) &&
      isTRUE(abs(zi(R[["G5"]], "p_one"))  < abs(zi(R[["G3"]], "p_one"))),
      sprintf("(zero z %.2f -> %.2f ; one z %.2f -> %.2f; the pooled adoption went +3.7 -> +2.0 and -6.1 -> -3.3)",
              zi(R[["G3"]], "p_zero"), zi(R[["G5"]], "p_zero"),
              zi(R[["G3"]], "p_one"),  zi(R[["G5"]], "p_one")))

  # (6) the documents must carry the corrected conclusions, not the driver's wrong ones
  fl <- function(f) gsub("[ \n]+", " ", paste(readLines(f, warn = FALSE), collapse = "\n"))
  cr <- fl("07_documentation/development_notes/CHANGE_REGISTER.md")
  ps <- fl("07_documentation/development_notes/PIPELINE_STATUS.md")
  vc <- fl("07_documentation/development_notes/VALIDATION_CAMPAIGN.md")
  chk("1w: the campaign carries Section 1w with the four defects and the matched configuration",
      grepl("## 1w.", vc, fixed = TRUE) && grepl("+0.28%", vc, fixed = TRUE) &&
      grepl("four defects", vc) && grepl("11,021", vc, fixed = TRUE))
  chk("1w: the register records D6 as ADOPT and says the driver's verdict was a DEFECT",
      grepl("D6 | **The gear-track ZINB EARNS its parameter", cr, fixed = TRUE) &&
      grepl("was a DEFECT, not a finding", cr, fixed = TRUE))
  chk("1w: D29 exists, is flagged above D3, and carries the measured span",
      # 2026-10-02 (B62): D29 CLOSED by rule (b); "it should outrank D3" left the status cell with it
      grepl("| D29 |", cr, fixed = TRUE) && grepl("CLOSED 2026-10-02 by rule D29-4(b)", cr, fixed = TRUE) &&
      grepl("+24.93%", cr, fixed = TRUE) && grepl("D29", ps, fixed = TRUE))
  # The documents DO quote the driver's conclusions, which they must: a review that does
  # not say what it is correcting cannot be checked against the output it corrects. What
  # must hold is that every quotation is marked wrong in the same breath, so the crude
  # "the phrase is absent" form is replaced by "the phrase never stands unqualified".
  chk("1w: the driver's wrong conclusions appear ONLY as quotations marked wrong",
      { bad <- character()
        for (nm in c(CHANGE_REGISTER = cr, PIPELINE_STATUS = ps)) {
          for (ph in c("D3: weekly", "D3 weekly", "D6: do not adopt", "D6 do not adopt")) {
            k <- gregexpr(ph, nm, fixed = TRUE)[[1]]
            if (k[1] < 0) next
            for (pos in k) {
              ctx <- substr(nm, pos, pos + 320)
              if (!grepl("were wrong|was wrong|was a DEFECT|not a finding", ctx))
                bad <- c(bad, substr(ctx, 1, 60)) } } }
        length(bad) == 0 },
      "(a quotation must carry its correction within the same sentence or two)")
  chk("1w: and the corrected answers are the ones stated as current",
      grepl("D6 is ADOPT", ps, fixed = TRUE) &&
      grepl("D3 is not a shore-resolution question", ps, fixed = TRUE))
  chk("1w: the recommendation CSV is kept as the record of what the driver said",
      file.exists("05_output/gear_ar_zi_2026-09-13_recommendation.csv"),
      paste("It is wrong and it stays: a review that deletes the output it corrects cannot",
            "be checked."))
  }
})

# ---------------------------------------------------------------------------
# 76. MARINE HAZARD EFFORT COVARIATES: SMALL CRAFT ADVISORIES AND BAR RESTRICTIONS
#     (2026-09-25, CHANGE_REGISTER A30 / B35 / B36).
#
#     The NWS VTEC archive for the Grays Harbor Bar and the coastal waters off Westport
#     becomes a per-day SCA-or-higher flag known on EVERY day of the window; the samplers'
#     "Bar Restrictions" tick becomes a per-day covariate observed on sampled days and
#     imputed on the rest. Both ride on the K_open block, so the assertions here are about
#     the R side only: the time arithmetic (DST-exact, half-open edges), the archive reader
#     and its coverage stop, the imputation and its fallback, the log-link screen, the
#     four modes, the collinearity guard, the design-matrix hand-off with a FRACTIONAL
#     column, the preps' wiring, the drivers' call, the shipped defaults (OFF), the
#     workbook, the builder's parser on the service's real JSON, and the batch runner's
#     stated rule. Every fixture is synthetic except the committed workbook and the JSON,
#     which is the text the service returned on 2026-09-25.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/bss_marine_hazard_covariates.R")
  TZ <- "America/Los_Angeles"
  utc <- function(x) as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

  # (a) time helpers: ISO parsing, DST-exact local windows, half-open overlap
  pu <- .mh_parse_utc(c("2025-01-03T21:00:00Z", "2025-01-03 21:00:00", "2025-01-03T21:00Z", "junk", NA))
  chk("mh: ISO 'Z', space and hour-minute forms parse to the same UTC instant; junk is NA",
      isTRUE(all.equal(as.numeric(pu[1]), as.numeric(utc("2025-01-03T21:00:00Z")))) &&
      isTRUE(all.equal(as.numeric(pu[2]), as.numeric(pu[1]))) && isTRUE(all.equal(as.numeric(pu[3]), as.numeric(pu[1]))) &&
      all(is.na(pu[4:5])))
  lt <- .mh_local_time(as.Date(c("2024-03-09", "2024-03-10", "2024-11-03")), 4, TZ)
  chk("mh: 04:00 local is 12Z in PST and 11Z in PDT, across both DST changes",
      identical(format(lt, "%H:%M", tz = "UTC"), c("12:00", "11:00", "12:00")))
  chk("mh: hour 24 is midnight starting the NEXT day",
      identical(format(.mh_local_time(as.Date("2024-06-01"), 24, TZ), "%Y-%m-%d %H:%M", tz = TZ), "2024-06-02 00:00"))
  chk("mh: a non-existent spring-forward time is moved FORWARD by the gap (02:30 -> 03:30 PDT = 10:30Z), whatever the platform returns for it",
      identical(format(.mh_local_time(as.Date(c("2024-03-09", "2024-03-10", "2024-03-11")), 2.5, TZ), "%Y-%m-%d %H:%M", tz = "UTC"),
                c("2024-03-09 10:30", "2024-03-10 10:30", "2024-03-11 09:30")))
  a <- utc("2025-01-10T12:00:00Z"); b <- utc("2025-01-11T00:00:00Z")
  chk("mh: overlap is half-open on both edges (start == end-of-window and end == start-of-window do not count)",
      identical(.mh_any_overlap(a, b, utc("2025-01-11T00:00:00Z"), utc("2025-01-12T00:00:00Z")), 0L) &&
      identical(.mh_any_overlap(a, b, utc("2025-01-09T00:00:00Z"), utc("2025-01-10T12:00:00Z")), 0L) &&
      identical(.mh_any_overlap(a, b, utc("2025-01-10T23:59:00Z"), utc("2025-01-12T00:00:00Z")), 1L) &&
      identical(.mh_any_overlap(a, b, utc("2025-01-09T00:00:00Z"), utc("2025-01-12T00:00:00Z")), 1L))
  chk("mh: no events -> all zeros, one per window", identical(.mh_any_overlap(c(a, a), c(b, b), utc(character(0)), utc(character(0))), c(0L, 0L)))

  # (b) the flag series on a synthetic archive: definition, zones, codes, coverage stop
  ev <- tibble::tibble(ugc = c("PZZ110", "PZZ156", "PZZ156", "PZZ156", "PZZ110"),
                       phenomena = c("SC", "GL", "MF", "SC", "SC"), significance = c("Y", "W", "Y", "Y", "Y"),
                       eventid = 1:5,
                       start = utc(c("2025-01-03T21:00:00Z", "2025-01-05T23:00:00Z", "2025-01-08T12:00:00Z", "2025-01-10T00:30:00Z", "2025-01-12T00:00:00Z")),
                       end   = utc(c("2025-01-04T12:00:00Z", "2025-01-07T00:00:00Z", "2025-01-08T20:00:00Z", "2025-01-10T00:20:00Z", "2025-01-12T06:00:00Z")),
                       product_id = NA_character_)
  ev$ps <- paste0(ev$phenomena, ".", ev$significance); ev$in_effect <- ev$end > ev$start
  attr(ev, "coverage") <- tibble::tibble(ugc = c("PZZ110", "PZZ156"), pull_start = as.Date("2025-01-01"), pull_end = as.Date("2025-01-31"))
  attr(ev, "coverage_inferred") <- FALSE
  Pm <- list(marine_hazard_window = c(4, 16), marine_hazard_tz = TZ)
  fl <- marine_hazard_flag_series(Pm, "2025-01-01", "2025-01-14", events = ev)
  # Jan 3: SCA on the bar 13:00 to Jan 4 04:00 PST -> Jan 3 flagged in the 04-16 window, Jan 4 NOT
  # (04:00 PST is 12Z, the event ends at 12Z, half-open); gale (coastal) 15:00 Jan 5 to
  # 16:00 Jan 6 PST -> Jan 5 and Jan 6; Jan 8 dense fog is not a hazard code; Jan 9/10 event
  # cancelled before it began (in_effect FALSE); Jan 11 16:00 to 22:00 PST -> outside the window.
  want <- setNames(rep(0L, 14), format(as.Date("2025-01-01") + 0:13))
  want[c("2025-01-03", "2025-01-05", "2025-01-06")] <- 1L
  chk("mh: the any-zone flag follows the codes, the window, the half-open edges and in_effect",
      identical(as.integer(fl$nws_sca_any), unname(want)), paste(fl$nws_sca_any, collapse = ""))
  chk("mh: single-zone flags split by zone", sum(fl$nws_sca_bar) == 1 && sum(fl$nws_sca_coastal) == 2)
  chk("mh: the all-day window catches the Jan 11 evening event",
      marine_hazard_flag_series(Pm, "2025-01-11", "2025-01-11", events = ev, window = c(0, 24))$nws_sca_any == 1)
  chk("mh: a window the archive does not cover STOPS with a message naming the builder",
      inherits(try(marine_hazard_flag_series(Pm, "2025-01-20", "2025-02-05", events = ev), silent = TRUE), "try-error") &&
      grepl("build_nws_marine_hazards", geterrmessage(), fixed = TRUE))
  chk("mh: an invalid window or unnamed zones are refused",
      inherits(try(marine_hazard_flag_series(Pm, "2025-01-01", "2025-01-02", events = ev, window = c(16, 4)), silent = TRUE), "try-error") &&
      inherits(try(marine_hazard_flag_series(modifyList(Pm, list(marine_hazard_zones = c("PZZ110", "PZZ156"))), "2025-01-01", "2025-01-02", events = ev), silent = TRUE), "try-error"))
  chk("mh: a non-default code set is honoured (gale-only sees only the gale)",
      sum(marine_hazard_flag_series(modifyList(Pm, list(marine_hazard_codes = "GL.W")), "2025-01-01", "2025-01-14", events = ev)$nws_sca_any) == 2)

  # (c) the bar-restriction tick: field start, observation rule, imputation, fallback
  dts <- as.Date("2025-01-01") + 0:99
  set.seed(7)
  shifts <- tibble::tibble(date = format(rep(dts, 2)), creel_location = rep(c("Grays Harbor", "Willapa Bay"), each = 100),
                           special_conditions = NA_character_)
  # Grays Harbor sampled on even days only (plus Jan 5); the option "appears" on Jan 5, when
  # the first marine token is written; restrictions tick on days 10, 12, 20, 22, 30, 40, 50.
  # Jan 2 and Jan 4 are sampled days BEFORE the field start: their blanks are not
  # observations of "no restriction", because the option did not yet exist on the form.
  gh_rows <- which(shifts$creel_location == "Grays Harbor")
  # ticks on every sampled day that carries a coastal SCA (days divisible by 6 from 12 on) plus
  # a few dry ones, so the logistic has a signal to find; days 10, 20, ... are the dry ticks
  shifts$special_conditions[gh_rows[c(10, 20, 30, 40, 50)]] <- "Cold, Bar Restrictions"
  shifts$special_conditions[gh_rows[c(12, 22)]] <- "bar restrictions, Razor Clam Opener"
  shifts$special_conditions[gh_rows[seq(18, 96, by = 6)]] <- "Bar Restrictions, Small Craft Advisory"
  shifts$special_conditions[gh_rows[5]] <- "Small Craft Advisory"       # the first appearance of a form token
  odd <- gh_rows[seq(1, 100, by = 2)]; shifts <- shifts[-odd[odd != gh_rows[5]], ]  # unsampled odd days (keep day 5)
  nd <- tibble::tibble(event_date = dts, nws_sca_any = as.integer(seq_along(dts) %% 3 == 0),
                       nws_sca_bar = as.integer(seq_along(dts) %% 4 == 0), nws_sca_coastal = as.integer(seq_along(dts) %% 3 == 0))
  Pb <- list(gh_creel_location = "Grays Harbor", bar_restriction_impute = "nws")
  br <- bar_restriction_series(Pb, dts, nws_day = nd, shifts = shifts)
  chk("bar: the option's first appearance on the form is detected from any marine token",
      identical(attr(br, "field_start"), as.Date("2025-01-05")))
  chk("bar: a sampled day before the field start is NA, not 0 (a blank before the option existed is not an observation)",
      all(is.na(br$bar_restriction_obs[c(2, 4)])) && all(br$bar_restriction_source[c(2, 4)] != "observed"))
  chk("bar: sampled days are 0 (blank) or 1 (any survey ticked, any case)",
      identical(br$bar_restriction_obs[c(10, 12, 20, 22, 6)], c(1L, 1L, 1L, 1L, 0L)))
  chk("bar: unsampled days are NA observations with an imputed value in [0, 1] and a source label",
      all(is.na(br$bar_restriction_obs[c(7, 9, 11)])) && all(br$bar_restriction[c(7, 9, 11)] >= 0 & br$bar_restriction[c(7, 9, 11)] <= 1) &&
      all(br$bar_restriction_source[c(7, 9, 11)] %in% c("imputed_nws", "imputed_mean")) &&
      all(br$bar_restriction_source[c(10, 6)] == "observed"))
  chk("bar: with 30 or more observed days the NWS logistic IS fitted, every unsampled day is 'imputed_nws', and an SCA day imputes higher than a dry one",
      { fit <- attr(br, "fit"); todo <- is.na(br$bar_restriction_obs)
        attr(br, "n_observed") >= 30 && !is.null(fit) && identical(fit$term, c("(Intercept)", "nws_sca_bar", "nws_sca_coastal", "winter")) &&
        all(br$bar_restriction_source[todo] == "imputed_nws") && any(grepl("imputed from logit", attr(br, "note"))) &&
        mean(br$bar_restriction[todo & nd$nws_sca_coastal == 1]) > mean(br$bar_restriction[todo & nd$nws_sca_coastal == 0]) },
      sprintf("(observed %s)", attr(br, "n_observed")))
  chk("bar: too few observed days for the logistic (under 30) falls back to the observed rate, with a note",
      { b2 <- bar_restriction_series(Pb, dts[1:20], nws_day = nd, shifts = shifts[as.Date(shifts$date) <= as.Date("2025-02-09"), ])
        attr(b2, "n_observed") < 30 && is.null(attr(b2, "fit")) &&
        all(b2$bar_restriction_source[is.na(b2$bar_restriction_obs)] == "imputed_mean") && any(grepl("observed rate", attr(b2, "note"))) })
  chk("bar: impute = 'mean' never fits the logistic",
      { b3 <- bar_restriction_series(modifyList(Pb, list(bar_restriction_impute = "mean")), dts, nws_day = nd, shifts = shifts)
        is.null(attr(b3, "fit")) && all(b3$bar_restriction_source[is.na(b3$bar_restriction_obs)] == "imputed_mean") })
  chk("bar: an unknown impute mode is refused",
      inherits(try(bar_restriction_series(modifyList(Pb, list(bar_restriction_impute = "zero")), dts, nws_day = nd, shifts = shifts), silent = TRUE), "try-error"))
  chk("bar: an explicit field start overrides detection",
      is.na(bar_restriction_series(modifyList(Pb, list(bar_restriction_field_start = "2025-01-15")), dts, nws_day = nd, shifts = shifts)$bar_restriction_obs[10]))

  # (d) the screen is log-link and finds a multiplicative effect the additive lm hides
  set.seed(11)
  d90 <- as.Date("2024-12-01") + 0:89
  flag <- as.integer(seq_along(d90) %% 3 == 0)
  mu <- 30 * exp(-1.2 * flag) * ifelse(weekdays(d90) %in% c("Saturday", "Sunday"), 2, 1)
  dwg <- list(shore_effort = tibble::tibble(event_date = d90, count_quantity = rpois(90, mu)),
              boat_effort  = tibble::tibble(event_date = d90, count_quantity = rpois(90, mu / 2)))
  flags <- tibble::tibble(event_date = d90, nws_sca_any = flag, bar_restriction_obs = ifelse(seq_along(d90) %% 2 == 0, flag, NA_integer_),
                          bar_restriction = flag)
  Psc <- list(days_wkend = c("Saturday", "Sunday"), crabbing_holiday_dates = as.Date(character(0)),
              marine_hazard_candidates_shore = "nws_sca_any", marine_hazard_candidates_boat = c("nws_sca_any", "bar_restriction"))
  sc <- marine_hazard_screen(dwg, flags, Psc)
  chk("screen: one row per candidate x population, carrying a rate ratio and an adjusted p",
      nrow(sc) == 3 && all(c("rate_ratio", "adj_p", "n_flag") %in% names(sc)) && all(is.finite(sc$rate_ratio)))
  chk("screen: the log-link recovers a multiplicative effect of exp(-1.2) = 0.30 within 0.10",
      all(abs(sc$rate_ratio[sc$covariate == "nws_sca_any"] - exp(-1.2)) < 0.10), paste(round(sc$rate_ratio, 3), collapse = ","))
  chk("screen: bar_restriction is screened on its OBSERVED days only",
      sc$n_days[sc$covariate == "bar_restriction"] == sum(!is.na(flags$bar_restriction_obs)))
  chk("screen: the per-count sensitivity columns are carried, and with one count per day they equal the daily-sum fit",
      all(c("n_counts", "rate_ratio_per_count", "adj_p_per_count") %in% names(sc)) &&
      identical(sc$n_counts, sc$n_days) && isTRUE(all.equal(sc$rate_ratio_per_count, sc$rate_ratio, tolerance = 1e-8)) &&
      isTRUE(all.equal(sc$adj_p_per_count, sc$adj_p, tolerance = 1e-8)))
  chk("screen: with two counts on every day the daily-sum and per-count rate ratios still agree (no count-frequency confound in the fixture)",
      { dwg2 <- list(shore_effort = dplyr::bind_rows(dwg$shore_effort, dwg$shore_effort), boat_effort = dplyr::bind_rows(dwg$boat_effort, dwg$boat_effort))
        s3 <- marine_hazard_screen(dwg2, flags, Psc)
        identical(s3$n_counts, 2L * s3$n_days) && isTRUE(all.equal(s3$rate_ratio_per_count, s3$rate_ratio, tolerance = 1e-6)) })
  chk("screen: a constant flag is reported, not fitted",
      { f2 <- flags; f2$nws_sca_any <- 0L
        s2 <- marine_hazard_screen(dwg, f2, Psc); all(grepl("constant", s2$note[s2$covariate == "nws_sca_any"])) })

  # (e) the four modes, the BH family, the one-NWS-per-population guard
  scr <- tibble::tibble(population = c("shore", "private_boat", "private_boat", "private_boat"),
                        covariate = c("nws_sca_any", "nws_sca_any", "nws_sca_bar", "bar_restriction"),
                        label = unname(.mh_labels[c("nws_sca_any", "nws_sca_any", "nws_sca_bar", "bar_restriction")]),
                        adj_estimate = c(-0.1, -1.2, -0.7, -0.9), rate_ratio = exp(c(-0.1, -1.2, -0.7, -0.9)),
                        adj_p = c(0.25, 0.0006, 0.004, 0.0003), note = "")
  Pa <- list(marine_hazard_mode = "auto", marine_hazard_auto_p = 0.05, marine_hazard_auto_p_adjust = "BH",
             marine_hazard_candidates_shore = "nws_sca_any",
             marine_hazard_candidates_boat = c("nws_sca_any", "nws_sca_bar", "bar_restriction"))
  sa <- marine_hazard_select(scr, Pa)
  chk("select auto: BH over the offered family keeps the boat terms and drops the shore one",
      length(sa$shore) == 0 && setequal(sa$private_boat, c("nws_sca_any", "bar_restriction")))
  chk("select auto: two NWS definitions clearing the screen -> the smaller adjusted p is kept and the other says why",
      !"nws_sca_bar" %in% sa$private_boat &&
      any(grepl("one NWS definition", sa$table$reason[sa$table$covariate == "nws_sca_bar"])))
  chk("select auto: the family size is pinned to what was offered (an NA test does not shrink it)",
      { s2 <- scr; s2$adj_p[3] <- NA_real_; s3 <- marine_hazard_select(s2, Pa)
        setequal(s3$private_boat, c("nws_sca_any", "bar_restriction")) && grepl("4 effort test", s3$note[1]) })
  chk("select auto: no screen -> nothing selected, with a note",
      length(marine_hazard_select(NULL, Pa)$private_boat) == 0 && length(marine_hazard_select(NULL, Pa)$note) > 0)
  chk("select off: nothing, regardless of the screen",
      length(marine_hazard_select(scr, modifyList(Pa, list(marine_hazard_mode = "off")))$private_boat) == 0)
  so <- marine_hazard_select(scr, modifyList(Pa, list(marine_hazard_mode = "on")))
  chk("select on: every candidate, still one NWS definition per population (the first in candidate order)",
      identical(so$shore, "nws_sca_any") && setequal(so$private_boat, c("nws_sca_any", "bar_restriction")))
  sm <- marine_hazard_select(scr, modifyList(Pa, list(marine_hazard_mode = "manual", marine_hazard_manual_boat = c("bar_restriction", "razor_nearby_dig"),
                                                       marine_hazard_manual_shore = character(0))))
  chk("select manual: honours the list, ignores an unknown name, leaves the shore empty",
      identical(sm$private_boat, "bar_restriction") && length(sm$shore) == 0)
  chk("select manual: a candidate not named is reported as 'not named'",
      any(sm$table$reason[sm$table$covariate == "nws_sca_any" & sm$table$population == "private_boat"] == "not named"))
  chk("select: an unknown candidate or mode is refused",
      inherits(try(marine_hazard_select(scr, modifyList(Pa, list(marine_hazard_candidates_boat = "wave_height"))), silent = TRUE), "try-error") &&
      inherits(try(marine_hazard_select(scr, modifyList(Pa, list(marine_hazard_mode = "yes"))), silent = TRUE), "try-error"))

  # (f) the hand-off to the K_open block: a fractional (imputed) column survives, and the
  #     identifiability rule counts observed 0/1 days only
  source("03_R_functions/bss_opener_covariates.R")
  dd <- mkdays("2024-12-01", 60)
  fl2 <- tibble::tibble(event_date = dd$event_date, nws_sca_any = as.integer(seq_len(60) %% 4 == 0),
                        bar_restriction = ifelse(seq_len(60) %% 2 == 0, as.numeric(seq_len(60) %% 4 == 0), 0.37))
  om <- opener_design_matrix(dd, character(0), fl2, list(opener_min_days = 10), extra = c("nws_sca_any", "bar_restriction"))
  X <- matrix(om$X_open, nrow = 60)
  chk("hand-off: two marine columns enter X_open with their labels, the imputed values intact",
      om$K_open == 2 && identical(om$labels, c("nws_sca_any", "bar_restriction")) && sum(X[, 2] == 0.37) == 30)
  chk("hand-off: a column with fewer than opener_min_days observed days on a side is dropped",
      opener_design_matrix(dd, character(0), fl2, list(opener_min_days = 20), extra = "bar_restriction")$K_open == 0)

  # (g) both preps carry the wiring; both drivers call the orchestration; the pooled Stan is untouched
  for (f in c("03_R_functions/prep_bss_crab_pooled.R", "03_R_functions/prep_bss_crab_gear.R")) {
    src <- readLines(f, warn = FALSE); src <- src[!grepl("^\\s*#", src)]
    chk(sprintf("%s: adds the marine selection to the K_open extras", basename(f)),
        # since B42 (2026-09-27) the selection reaches the fit through marine_hazard_terms_for(), which
        # confines it to the configured gear regime(s); section 79 tests the call site and the note
        any(grepl("marine_extra <- marine_hazard_terms_for(params, population_name, gear_regime)", src, fixed = TRUE)) &&
        any(grepl("extra = c(razor_extra, marine_extra)", src, fixed = TRUE)))
  }
  for (drv in list.files("01_BSS_models", pattern = "\\.Rmd$", full.names = TRUE)) {
    d <- readLines(drv, warn = FALSE); d <- d[!grepl("^\\s*#", d)]
    chk(sprintf("%s: calls marine_hazard_prepare() and takes its params back", basename(drv)),
        any(grepl("marine_hazard <- marine_hazard_prepare(dwg, params, output_dir = output_dir)", d, fixed = TRUE)) &&
        any(grepl("params <- marine_hazard$params", d, fixed = TRUE)))
  }
  chk("no Stan model mentions the marine covariates: they are columns of X_open, not new terms",
      !any(grepl("marine|nws_sca|bar_restriction", unlist(lapply(list.files("02_stan_models", pattern = "\\.stan$", full.names = TRUE), readLines, warn = FALSE)))))

  # (h) the orchestration reads NOTHING when off, and stops when on without an archive
  Poff <- list(marine_hazard_mode = "off", marine_hazard_file = "no-such-workbook.xlsx")
  moff <- marine_hazard_prepare(NULL, Poff, quiet = TRUE)
  chk("prepare off: reads nothing (a missing archive is not even noticed) and installs an empty selection",
      isFALSE(moff$active) && identical(moff$params$marine_hazard_selected, list(shore = character(0), private_boat = character(0))))
  chk("prepare on: a missing archive STOPS rather than fitting without the covariate",
      inherits(try(marine_hazard_prepare(NULL, modifyList(Poff, list(marine_hazard_mode = "on", est_date_start = "2024-09-16", est_date_end = "2025-09-15")), quiet = TRUE), silent = TRUE), "try-error"))

  # (i) the shipped defaults
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  # 2026-09-27: A30 ADOPTED. Until then this asserted "off".
  chk("shipped: THE METHOD OF RECORD: manual, boat nws_sca_any, no shore term, the all-gear fit only (A30, 2026-09-27)",
      identical(rc$marine_hazard_mode, "manual") && identical(rc$marine_hazard_manual_boat, "nws_sca_any") &&
      identical(rc$marine_hazard_manual_shore, character(0)) && identical(rc$marine_hazard_gear_regimes, "all_gear"))
  chk("shipped: no marine candidate is offered to the shore; the boat's candidates are SCA (any zone) and the bar tick, and the adopted term is among them",
      identical(rc$marine_hazard_candidates_shore, character(0)) && identical(rc$marine_hazard_candidates_boat, c("nws_sca_any", "bar_restriction")) &&
      rc$marine_hazard_manual_boat %in% rc$marine_hazard_candidates_boat)
  chk("shipped: the window is 04:00-16:00 local, the zones are named bar/coastal, the codes are SCA-or-higher incl. the pre-2019 SCA codes",
      identical(rc$marine_hazard_window, c(4, 16)) && identical(names(rc$marine_hazard_zones), c("bar", "coastal")) &&
      all(c("SC.Y", "RB.Y", "SW.Y", "SI.Y", "GL.W", "SR.W") %in% rc$marine_hazard_codes))
  chk("shipped: bar restrictions impute from the archive ('nws'), BH over the marine family",
      identical(rc$bar_restriction_impute, "nws") && identical(rc$marine_hazard_auto_p_adjust, "BH"))

  # (j) the committed workbook and its builder
  wbp <- "04_input_files/nws_marine_hazards.xlsx"
  chk("workbook: nws_marine_hazards.xlsx exists with the reader's columns",
      file.exists(wbp) && all(c("ugc", "phenomena", "significance", "start_utc", "end_utc", "product_id", "in_effect", "pull_start", "pull_end", "source") %in%
                              names(read_input_workbook(wbp, sheet = "data"))))
  wb <- as.data.frame(read_input_workbook(wbp, sheet = "data"))
  chk("workbook: both zones, and the pull window covers the canonical 2024-25 season",
      setequal(unique(wb$ugc), c("PZZ110", "PZZ156")) && all(as.Date(wb$pull_start) <= as.Date("2024-09-16")) && all(as.Date(wb$pull_end) >= as.Date("2025-09-15")))
  chk("workbook: times are ISO-8601 UTC text and in_effect is end > start",
      all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", c(wb$start_utc, wb$end_utc))) &&
      identical(as.integer(wb$in_effect), as.integer(.mh_parse_utc(wb$end_utc) > .mh_parse_utc(wb$start_utc))))
  chk("workbook: its provenance is stated on every row (a live pull, or the transcription it replaced), and the transcription is kept as the verification record",
      all(grepl("transcription|IEM", wb$source)) && file.exists("04_input_files/raw/nws_marine_hazards_iem_transcription_2026-09-25.csv"))
  chk("workbook: the reader's coverage attribute comes from the pull window, not inferred",
      { evw <- marine_hazard_events(list(marine_hazard_file = "nws_marine_hazards.xlsx")); isFALSE(attr(evw, "coverage_inferred")) && nrow(attr(evw, "coverage")) == 2 })
  # the builder's parser on the text the service actually returned (2026-09-25)
  if (requireNamespace("jsonlite", quietly = TRUE) && requireNamespace("writexl", quietly = TRUE)) {
    eb <- new.env(); sys.source("04_input_files/build_nws_marine_hazards.R", envir = eb)
    js <- paste0('{"events": [{"url": "/vtec/?year=2024&wfo=KSEW&phenomena=SC&significance=Y&eventid=0189", "issue": "2025-01-03T21:00:00Z", ',
                 '"expire": "2025-01-04T12:00:00Z", "eventid": 189, "phenomena": "SC", "hvtec_nwsli": null, "significance": "Y", "wfo": "SEW", ',
                 '"name": "Small Craft Advisory", "ph_name": "Small Craft", "sig_name": "Advisory", "ugc": "PZZ110", "product_id": "202501031046-KSEW-WHUS76-MWWSEW"}, ',
                 '{"url": "/vtec/?year=2025&wfo=KSEW&phenomena=SC&significance=Y&eventid=0003", "issue": "2025-01-10T12:00:00Z", "expire": "2025-01-12T11:01:00Z", ',
                 '"eventid": 3, "phenomena": "SC", "hvtec_nwsli": null, "significance": "Y", "wfo": "SEW", "name": "Small Craft Advisory", "ph_name": "Small Craft", ',
                 '"sig_name": "Advisory", "ugc": "PZZ110", "product_id": "202501092254-KSEW-WHUS76-MWWSEW"}], "generated_at": "2026-09-25T18:24:30Z"}')
    pr <- eb$parse_iem_vtec(js, "PZZ110")
    chk("builder: parses the service's JSON into one row per segment with issue / expire / product_id",
        nrow(pr) == 2 && identical(pr$issue, c("2025-01-03T21:00:00Z", "2025-01-10T12:00:00Z")) && identical(pr$eventid, c(189L, 3L)) &&
        identical(pr$product_id[1], "202501031046-KSEW-WHUS76-MWWSEW"))
    chk("builder: an empty events array is zero rows, not an error", nrow(eb$parse_iem_vtec('{"events": [], "generated_at": "x"}', "PZZ110")) == 0)
    wb2 <- eb$build_nws_workbook(pr, "PZZ110", as.Date("2025-01-01"), as.Date("2025-01-31"), "test")
    chk("builder: the workbook carries the announcement time from product_id, in_effect, and the pull window on every row",
        identical(wb2$product_issued_utc[1], "2025-01-03T10:46:00Z") && all(wb2$in_effect == 1L) &&
        all(wb2$pull_start == "2025-01-01") && identical(wb2$hazard[1], "Small Craft Advisory") && identical(wb2$zone_name[1], "Grays Harbor Bar"))
    chk("builder: pulls by calendar year and never de-duplicates by eventid alone",
        { src <- readLines("04_input_files/build_nws_marine_hazards.R", warn = FALSE)
          any(grepl("seq(as.integer(format(sdate", src, fixed = TRUE)) && any(grepl("distinct(ugc, phenomena, significance, eventid, issue, expire, product_id", src, fixed = TRUE)) })
  } else cat("NOTE  builder: jsonlite or writexl absent; the parser assertions are skipped\n")

  # (k) the batch runner: the rule before the run, the pin, the reference, the covariate clauses
  rf <- "06_diagnostics/run_marine_hazard_batch_2026-09-25.R"
  src <- readLines(rf, warn = FALSE)
  chk("runner: exists, ships DRY_RUN <- TRUE and RESUME by digest", file.exists(rf) && any(grepl("^DRY_RUN <- TRUE", src)) && any(grepl("MH_STAGE.txt", src, fixed = TRUE)))
  chk("runner: the decision rule is stated in the header BEFORE the run, including what LOO cannot measure",
      any(grepl("THE DECISION RULE, STATED HERE BEFORE THE RUN", src, fixed = TRUE)) && any(grepl("WHAT THIS RUN CANNOT MEASURE", src, fixed = TRUE)))
  chk("runner: the pin carries all nine per-season keys and the flag definition",
      all(vapply(c("est_date_start", "season_filter", "pot_open_date", "census_end_date", "commercial_opener",
                   "marine_hazard_window", "marine_hazard_codes", "bar_restriction_impute"),
                 function(k) any(grepl(paste0("(^|[ ,(])", k, " = "), src)), logical(1))))
  chk("runner: the baseline is judged bit-identical to R4 and the bar term against the SCA rung (rule 6)",
      any(grepl("REF_R4 <- \"20260910/pooled-CPUE-IMP-R4-shore-tau-newf\"", src, fixed = TRUE)) &&
      any(grepl('list(sid = "M4", ctl = "M2", fit = FIT_BOAT,  cov = "bar_restriction")', src, fixed = TRUE)))
  chk("runner: a paired elpd inside +-2 SE is REVIEW, never FAIL",
      any(grepl('if (ratio > 2) "PASS" else if (ratio < -2) "FAIL" else "REVIEW"', src, fixed = TRUE)))
  chk("runner: the recommendation needs every one of rules 1 to 5 PRESENT and EVALUATED; an absent or non-computable clause is 'open', never adoptable (B33)",
      any(grepl("present <- vapply(1:5, function(n) any(grepl(paste0(\"^rule \", n), rows$criterion)), logical(1))", src, fixed = TRUE)) &&
      any(grepl("open: rule(s) %s not evaluated", src, fixed = TRUE)) &&
      any(grepl("^NOT COMPUTABLE|^term not in this fit", src, fixed = TRUE)) &&
      !any(grepl("\"adoptable, no sampled-day gain\"", src, fixed = TRUE)))
  chk("runner: a missing convergence report, AR log or adequacy row is REVIEW in rules 1 and 2, not PASS or FAIL",
      any(grepl('if (is.na(L$gate_all_pass)) "REVIEW" else if (isTRUE(L$gate_all_pass)) "PASS" else "FAIL"', src, fixed = TRUE)) &&
      any(grepl('if (!known) "REVIEW" else if (pf <= 0.15 && share <= 0.05) "PASS" else "FAIL"', src, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 77. THE MARINE HAZARD LADDER RAN: ITS EVIDENCE, THE TWO DEFECTS IT EXPOSED IN ITS OWN
#     DRIVER, THE REPRODUCIBLE PORT TOTAL, AND THE BLOCK CROSS-VALIDATION (2026-09-26;
#     CHANGE_REGISTER Section 1x, B38, B39, two C rows).
#
#     Three things are pinned. (1) THE EVIDENCE, read from the committed rung folders and
#     the ladder CSVs, so that every number the register and the campaign quote can be
#     re-read: M1's components identical to R4, the boat SCA coefficient and its interval,
#     the paired trailer elpd recomputed from the pointwise files, the shore term's interval
#     containing zero, the port totals. Skipped with a NOTE where the folders are absent (a
#     sparse checkout), exactly as section 75 does. (2) THE FIXES: the port tolerance reads
#     the documented jitter with the mechanism beside it, the M5 row is vectorised, the
#     runner declares the post-run code equivalent by naming the ACTUAL current fingerprint,
#     bss_stan_fit() rebuilds the draw permutation from the Stan seed after each fit and
#     restores the caller's RNG (the pure helper exercised), both drivers seed the census draw. (3) THE BLOCK CV,
#     against an exact answer: PSIS leave-block-out on a conjugate Poisson-gamma model within
#     0.01 nats of the closed form, single-observation blocks equal to loo::loo pointwise, the
#     log-likelihood rebuild equal to direct arithmetic for gear, trailer and OSP under both
#     scalings, the self-check accepting a matching file and refusing a mismatched one, the
#     comparison excluding unreliable weeks, the drivers calling it, the runner shipping
#     DRY_RUN TRUE with its rule stated.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  flat <- function(t) gsub("[ \n]+", " ", t)
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  source("03_R_functions/bss_block_cv.R")

  # ---- (1) the evidence, from the committed folders ---------------------------------
  R4 <- "05_output/20260910/pooled-CPUE-IMP-R4-shore-tau-newf"
  MH <- c(M1 = "05_output/20260925/pooled-CPUE-MH-M1-off", M2 = "05_output/20260925/pooled-CPUE-MH-M2-sca",
          M3 = "05_output/20260925/pooled-CPUE-MH-M3-bar", M4 = "05_output/20260926/pooled-CPUE-MH-M4-both",
          M5 = "05_output/20260926/pooled-CPUE-MH-M5-auto")
  BAG <- "private_boat_all_gear_Dungeness_Kept"; SAG <- "shore_all_gear_Dungeness_Kept"
  lad_f <- "05_output/marine_hazard_2026-09-25_ladder.csv"; rec_f <- "05_output/marine_hazard_2026-09-25_recommendation.csv"
  chk("1x: the ladder, verdicts and recommendation CSVs are committed",
      file.exists(lad_f) && file.exists("05_output/marine_hazard_2026-09-25_verdicts.csv") && file.exists(rec_f))
  if (all(dir.exists(c(R4, MH))) && file.exists(lad_f)) {
    lad <- read.csv(lad_f, stringsAsFactors = FALSE); rownames(lad) <- lad$rung
    comp <- function(dir, key) { x <- read.csv(file.path(dir, "pe_vs_bss_comparison.csv"), stringsAsFactors = FALSE); x$BSS_catch[x$component == key] }
    keys <- c("shore (Pot closure)", "shore (All gear)", "private_boat (Pot closure)", "private_boat (All gear)")
    chk("1x: M1's four BSS components are identical to R4's, to the crab",
        identical(vapply(keys, function(k) comp(MH["M1"], k), numeric(1)), vapply(keys, function(k) comp(R4, k), numeric(1))))
    pt <- function(dir) read.csv(file.path(dir, "port_total_Dungeness_Kept.csv"))$BSS_median[2]
    chk("1x: M1's port total differs from R4's by the documented draw-permutation jitter (under 0.3%), with identical components",
        { d <- abs(pt(MH["M1"]) - pt(R4)) / pt(R4); d > 0 && d < 0.003 }, sprintf("(%.4f%%)", 100 * abs(pt(MH["M1"]) - pt(R4)) / pt(R4)))
    bo <- function(dir, fit, cov) { lab <- read.csv(file.path(dir, sprintf("opener_covariates_%s.csv", fit)), stringsAsFactors = FALSE)
      fs <- read.csv(file.path(dir, sprintf("bss_full_summary_%s.csv", fit)), row.names = 1, check.names = FALSE)
      fs[lab$parameter[lab$opener == cov], c("mean", "2.5%", "97.5%")] }
    b2 <- bo(MH["M2"], BAG, "nws_sca_any"); s2 <- bo(MH["M2"], SAG, "nws_sca_any"); bar4 <- bo(MH["M4"], BAG, "bar_restriction")
    chk("1x: the boat SCA term in M2 is identified with the quoted interval (-1.16 [-1.45, -0.87])",
        abs(b2$mean + 1.163) < 0.005 && abs(b2[["2.5%"]] + 1.452) < 0.005 && abs(b2[["97.5%"]] + 0.868) < 0.005)
    chk("1x: the shore SCA term's interval contains zero (the FAIL the runner reported is correct)",
        s2[["2.5%"]] < 0 && s2[["97.5%"]] > 0 && abs(s2$mean - 0.038) < 0.005)
    chk("1x: the bar term beyond the archive is identified by a hair (upper bound -0.004) at rate ratio 0.70",
        abs(bar4[["97.5%"]] + 0.004) < 0.002 && abs(exp(bar4$mean) - 0.695) < 0.005)
    pw <- function(dir, stream, fit) read.csv(file.path(dir, sprintf("loo_pointwise_%s_%s.csv", stream, fit)))
    a <- pw(MH["M1"], "trailer", BAG); b <- pw(MH["M2"], "trailer", BAG); d <- b$elpd_loo - a$elpd_loo
    chk("1x: the paired trailer elpd of M2 against M1 recomputes to +7.6 at 5.2 SE (1.46 SE): real, and short of +2",
        identical(a$obs_index, b$obs_index) && abs(sum(d) - 7.55) < 0.05 && abs(sqrt(length(d)) * sd(d) - 5.18) < 0.05 &&
        sum(d) / (sqrt(length(d)) * sd(d)) < 2)
    chk("1x: the boat all-gear estimate ROSE under the SCA term (+3.8% M2, +3.4% M4) and the port with it (+1.6%, +1.4%)",
        abs(lad["M2", "boat_ag"] / lad["M1", "boat_ag"] - 1.0376) < 0.001 && abs(lad["M4", "boat_ag"] / lad["M1", "boat_ag"] - 1.0344) < 0.001 &&
        abs(lad["M2", "port"] / lad["M1", "port"] - 1.0157) < 0.001)
    chk("1x: M5 (auto) reproduces M4's boat fit bit for bit: its boat components and B_open are M4's",
        lad["M5", "boat_ag"] == lad["M4", "boat_ag"] && lad["M5", "B_sca_boat"] == lad["M4", "B_sca_boat"] &&
        lad["M5", "B_bar_boat"] == lad["M4", "B_bar_boat"])
    chk("1x: every rung passed the gate on all four fits at M1's resolutions",
        all(lad$gate_all_pass) && all(lad$n_fits_bss == 4) && all(lad$shore_res == "weekly") && all(lad$boat_res == "monthly"))
    ar <- function(dir) read.csv(file.path(dir, "ar_escalation_log.csv"), stringsAsFactors = FALSE)
    chk("1x: the unidentified shore term doubled the shore fits' divergences (178 -> 329 all-gear; 76 -> 158 pot-closure)",
        { a1 <- ar(MH["M1"]); a2 <- ar(MH["M2"])
          a1$divergences[a1$fit == SAG] == 178 && a2$divergences[a2$fit == SAG] == 329 &&
          a1$divergences[a1$fit == "shore_ring_net_only_Dungeness_Kept"] == 76 && a2$divergences[a2$fit == "shore_ring_net_only_Dungeness_Kept"] == 158 })
    rec <- read.csv(rec_f, stringsAsFactors = FALSE)
    # the M6 rows (2026-09-27, "split by season") answer rule 9, tested in section 79; the 1x reading is of M1 to M5
    rec1x <- rec[!grepl("split by season", rec$item), ]
    chk("1x: the recommendation CSV says what the rule says: shore do not adopt; boat SCA and bar 'no sampled-day gain'",
        nrow(rec1x) == 3 && grepl("do not adopt", rec1x$recommendation[grepl("shore", rec1x$item)]) &&
        all(grepl("no sampled-day gain", rec1x$recommendation[grepl("boat", rec1x$item)])))
    # the permutation is what moved: identical seeded subsample indices, different values
    ds <- function(dir) read.csv(file.path(dir, sprintf("bss_draws_summed_%s.csv", BAG)))
    chk("1x: bss_draws_summed_* carries the SAME draw indices in R4 and M1 and DIFFERENT values at them (the permutation, not the fit)",
        identical(ds(R4)$draw, ds(MH["M1"])$draw) && !identical(ds(R4)$E_sum, ds(MH["M1"])$E_sum))
  } else cat("NOTE  1x: the marine rung folders or the ladder CSVs are absent; the evidence assertions are skipped\n")

  # ---- (2) the fixes ----------------------------------------------------------------
  rf <- "06_diagnostics/run_marine_hazard_batch_2026-09-25.R"; rsrc <- readLines(rf, warn = FALSE)
  chk("1x fix: the port reproduction is judged at the documented jitter (0.3%), with the mechanism written beside the threshold",
      any(grepl('"within 0.3%"', rsrc, fixed = TRUE)) && any(grepl("if (isTRUE(abs(.pct(p, pr)) < 0.3)) \"PASS\" else \"FAIL\"", rsrc, fixed = TRUE)) &&
      any(grepl("UNSEEDED sample.int()", rsrc, fixed = TRUE)) && !any(grepl('"within 0.05%"', rsrc, fixed = TRUE)))
  chk("1x fix: the M5 REPORTED row formats every candidate's adjusted p, not the first one recycled",
      any(grepl("fmt(suppressWarnings(as.numeric(sel$p_adj)), 4)", rsrc, fixed = TRUE)) && !any(grepl("fmt(.num1(sel$p_adj), 4)", rsrc, fixed = TRUE)))
  chk("1x fix: the M0 coverage row no longer calls the committed workbook a transcription",
      !any(grepl("The committed workbook is a transcription", rsrc, fixed = TRUE)))
  # the runner's own fingerprint, lifted from its text, so CODE_EQUIVALENT_MH can be checked against the real tree
  eR <- new.env(); eR$.here <- function(...) file.path(getwd(), ...); eR$POOLED_RMD <- file.path(getwd(), "01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd")
  lift <- function(name) { i <- grep(sprintf("^%s <- function", gsub(".", "\\.", name, fixed = TRUE)), rsrc); j <- i; while (!grepl("^\\}", rsrc[j])) j <- j + 1L; eval(parse(text = rsrc[i:j]), envir = eR) }
  lift("digest_or_hash"); lift(".code_fingerprint")
  i <- grep("^CODE_EQUIVALENT_MH <- list\\(", rsrc); j <- i; while (!grepl("^\\)", rsrc[j])) j <- j + 1L
  eval(parse(text = rsrc[i:j]), envir = eR)
  i <- grep("^CODE_NOT_EQUIVALENT_MH <- list\\(", rsrc)
  if (length(i)) { j <- i; while (!grepl("^\\)", rsrc[j])) j <- j + 1L; eval(parse(text = rsrc[i:j]), envir = eR) }
  cur <- eR$.code_fingerprint()
  .pair <- function(lst, rec) any(vapply(names(lst %||% list()), function(k) identical(strsplit(k, " => ", fixed = TRUE)[[1]], c(rec, cur)), logical(1)))
  chk("1x fix: the five rungs' recorded fingerprint is declared against the ACTUAL current one (equivalent, or examined and NOT equivalent), with a reason",
      length(eR$CODE_EQUIVALENT_MH) >= 1 && all(nchar(unlist(eR$CODE_EQUIVALENT_MH)) > 60) &&
      (.pair(eR$CODE_EQUIVALENT_MH, "stan:65b5adeb drivers:ae200663 fns:094c314f") ||
       (.pair(eR$CODE_NOT_EQUIVALENT_MH, "stan:65b5adeb drivers:ae200663 fns:094c314f") &&
        all(grepl("NOT inference-equivalent", unlist(eR$CODE_NOT_EQUIVALENT_MH), fixed = TRUE)))),
      sprintf("(current is %s)", cur))
  chk("1x fix: the run's stamps carry the recorded fingerprint the declaration names",
      { st <- file.path(MH, "MH_STAGE.txt"); if (!all(file.exists(st))) TRUE else
          all(vapply(st, function(f) any(grepl("^code: stan:65b5adeb drivers:ae200663 fns:094c314f$", readLines(f, warn = FALSE))), logical(1))) })
  chk("1x fix: the runner honours a declared-equivalent pair (INFO, not REVIEW) on RESUME",
      any(grepl("eq <- CODE_EQUIVALENT_MH[[paste(cd, .code_fingerprint(), sep = \" => \")]]", rsrc, fixed = TRUE)) &&
      any(grepl("declared EQUIVALENT", rsrc, fixed = TRUE)))

  # B38: the permutation is REBUILT after the fit returns (rstan draws it in the chain's own
  # process with cores > 1, so seeding this process before the call reaches nothing), one
  # seeded sample.int() per chain of the length it replaces, caller's RNG restored.
  fsrc <- readLines("03_R_functions/bss_stan_fit.R", warn = FALSE)
  chk("B38: bss_stan_fit() returns bss_seed_permutation(fit, list(...)$seed), after the usability assert",
      { a <- grep("bss_assert_fit_usable(fit, label = label, console = console)", fsrc, fixed = TRUE)
        b <- which(fsrc == "  bss_seed_permutation(fit, list(...)$seed)")
        length(a) == 1 && length(b) == 1 && b == a + 1 })
  chk("B38: the rebuild is applied per chain to fit@sim$permutation, which is what rstan's extract(permuted = TRUE) indexes",
      any(grepl("fit@sim$permutation <- new", fsrc, fixed = TRUE)) && any(grepl("lapply(perm, function(p) sample.int(length(p)))", fsrc, fixed = TRUE)))
  source("03_R_functions/bss_stan_fit.R")
  perm <- list(sample.int(2500), sample.int(2500), sample.int(2500), sample.int(2500))
  set.seed(999); runif(1); p1 <- .bss_seeded_perm(perm, 20260619); a1 <- runif(1)
  set.seed(999); runif(1); p2 <- .bss_seeded_perm(perm, 20260619); b1 <- runif(1)
  set.seed(999); runif(1); c1 <- runif(1)
  chk("B38: the same Stan seed gives the same per-chain permutations every time, each a permutation of the length it replaces",
      identical(p1, p2) && identical(lengths(p1), lengths(perm)) && all(vapply(p1, function(x) identical(sort(x), 1:2500), logical(1))))
  chk("B38: the caller's RNG stream is restored exactly (the next draw is what it would have been without the call)",
      identical(a1, b1) && identical(a1, c1))
  chk("B38: a different seed gives different permutations; an unusable seed or an empty list gives NULL and the fit is returned untouched",
      !identical(p1, .bss_seeded_perm(perm, 1)) && is.null(.bss_seeded_perm(perm, NA)) && is.null(.bss_seeded_perm(perm, 1e12)) &&
      is.null(.bss_seeded_perm(list(), 5)) && identical(bss_seed_permutation(list(x = 1), 5), list(x = 1)))
  for (drv in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")) {
    d <- readLines(drv, warn = FALSE)
    # B46 (2026-09-28): the seed and the draw are one call, bss_with_seed(), which also restores the RNG
    jj <- grep(".cc_draws <- bss_with_seed(params$bss_seed %||% 1L, pmax(rnorm(n_draws_max, .cc$Dungeness_Kept, .cc$carried_se)", d, fixed = TRUE)
    chk(sprintf("B38: %s seeds the census draw (bss_with_seed, caller's RNG restored)", basename(drv)), length(jj) == 1)
    chk(sprintf("B39: %s calls write_block_cv_diagnostics() after write_loo_diagnostics()", basename(drv)),
        { a <- grep("write_loo_diagnostics(b$fit, b$bss_data, b$days_ss, label, output_dir)", d, fixed = TRUE)
          b <- grep("write_block_cv_diagnostics(b$fit, b$bss_data, b$days_ss, label, output_dir)", d, fixed = TRUE)
          length(a) == 1 && length(b) == 1 && b > a })
  }

  # ---- (3) the block CV against exact answers ----------------------------------------
  set.seed(1)
  a0 <- 2; b0 <- 0.5; n <- 60; y <- rpois(n, 4); block <- rep(sprintf("W%02d", 1:12), each = 5); S <- 20000
  lam_s <- rgamma(S, a0 + sum(y), b0 + n)
  ll <- sapply(y, function(yi) dpois(yi, lam_s, log = TRUE))
  tab <- bss_block_psis_loo(ll, block)
  exact <- sapply(unique(block), function(bk) { j <- block == bk; m <- sum(j); sk <- sum(y[j]); a2 <- a0 + sum(y[!j]); b2 <- b0 + sum(!j)
    lgamma(a2 + sk) - lgamma(a2) - sum(lgamma(y[j] + 1)) + a2 * log(b2) - (a2 + sk) * log(b2 + m) })
  chk("block CV: PSIS leave-block-out equals the EXACT conjugate leave-block-out predictive within 0.01 nats on every block, all k small",
      nrow(tab) == 12 && max(abs(tab$elpd_block - exact)) < 0.01 && max(tab$pareto_k) < 0.5, sprintf("(max |diff| %.4f)", max(abs(tab$elpd_block - exact))))
  chk("block CV: the in-sample joint lpd exceeds the held-out elpd on every block (p_eff > 0)", all(tab$p_eff > 0))
  lo <- suppressWarnings(loo::loo(ll, r_eff = rep(1, n)))
  chk("block CV: single-observation blocks reproduce loo::loo's pointwise elpd_loo",
      max(abs(bss_block_psis_loo(ll, as.character(seq_len(n)))$elpd_block - lo$pointwise[, "elpd_loo"])) < 1e-8)
  tb <- tab; tb$elpd_block <- tab$elpd_block + c(0.4, 0.6); cmp <- bss_block_cv_compare(tab, tb)
  chk("block CV: the paired comparison sums the per-block differences over reliable blocks, with a paired SE",
      cmp$n_used == 12 && abs(cmp$diff - 6.0) < 1e-9 && abs(cmp$se - sqrt(12) * sd(rep(c(0.4, 0.6), 6))) < 1e-9)
  tb$pareto_k[1:5] <- 0.9; cmp2 <- bss_block_cv_compare(tab, tb)
  chk("block CV: a block unreliable in EITHER table is excluded and counted", cmp2$n_used == 7 && cmp2$n_dropped == 5 && abs(cmp2$reliable_share - 7 / 12) < 1e-9)
  chk("block CV: tables with different observations per block are refused",
      { tc <- tb; tc$n_obs[1] <- tc$n_obs[1] + 1L; inherits(try(bss_block_cv_compare(tab, tc), silent = TRUE), "try-error") })
  # the log-likelihood rebuild against direct arithmetic, all three streams, both OSP scalings, both trailer forms
  D <- 25; nd <- 300
  lamE <- array(exp(rnorm(nd * D, log(20), 0.3)), c(nd, 1, D, 1))
  dr <- list(lambda_E_S = lamE, r_E = rgamma(nd, 20, 2), R_G = runif(nd, 2, 3), R_G_boat = runif(nd, 1.5, 2.5),
             L_out = matrix(runif(nd * D, 2, 3), nd, D), r_OSP = rgamma(nd, 10, 1), kappa_OSP = runif(nd, 0.3, 0.5))
  sdt <- list(Gear_n = 10L, Gear_I = rpois(10, 40), day_Gear = sample(D, 10, TRUE), section_Gear = rep(1L, 10),
              T_n = 8L, T_I = rpois(8, 8), day_T = sample(D, 8, TRUE), section_T = rep(1L, 8),
              OSP_n = 6L, OSP_I = rpois(6, 25), day_OSP = sample(D, 6, TRUE), section_OSP = rep(1L, 6), osp_scale_is_tau = 1L)
  g <- bss_effort_loglik(dr, sdt, "gear"); tr <- bss_effort_loglik(dr, sdt, "trailer"); o <- bss_effort_loglik(dr, sdt, "osp")
  chk("block CV: the gear log-likelihood rebuild is NB2(lambda_E * R_G, r_E) to machine precision",
      max(abs(g$ll - sapply(1:10, function(i) dnbinom(sdt$Gear_I[i], mu = lamE[, 1, sdt$day_Gear[i], 1] * dr$R_G, size = dr$r_E, log = TRUE)))) < 1e-10)
  chk("block CV: the trailer rebuild is NB2(lambda_E / R_G_boat, r_E), and R_T fits use lambda_E * R_T",
      max(abs(tr$ll - sapply(1:8, function(i) dnbinom(sdt$T_I[i], mu = lamE[, 1, sdt$day_T[i], 1] / dr$R_G_boat, size = dr$r_E, log = TRUE)))) < 1e-10 &&
      { d3 <- dr; d3$R_G_boat <- NULL; d3$R_T <- runif(nd, 0.4, 0.6); t3 <- bss_effort_loglik(d3, sdt, "trailer")
        max(abs(t3$ll - sapply(1:8, function(i) dnbinom(sdt$T_I[i], mu = lamE[, 1, sdt$day_T[i], 1] * d3$R_T, size = dr$r_E, log = TRUE)))) < 1e-10 })
  chk("block CV: the OSP rebuild is NB2((lambda_E / R_G_boat) * L[day], r_OSP) under osp_scale_is_tau and * kappa_OSP otherwise",
      max(abs(o$ll - sapply(1:6, function(i) dnbinom(sdt$OSP_I[i], mu = lamE[, 1, sdt$day_OSP[i], 1] / dr$R_G_boat * dr$L_out[, sdt$day_OSP[i]], size = dr$r_OSP, log = TRUE)))) < 1e-10 &&
      { s2 <- sdt; s2$osp_scale_is_tau <- 0L; o2 <- bss_effort_loglik(dr, s2, "osp")
        max(abs(o2$ll - sapply(1:6, function(i) dnbinom(sdt$OSP_I[i], mu = lamE[, 1, sdt$day_OSP[i], 1] / dr$R_G_boat * dr$kappa_OSP, size = dr$r_OSP, log = TRUE)))) < 1e-10 })
  chk("block CV: an empty stream is a zero-column matrix, not an error", ncol(bss_effort_loglik(dr, list(Gear_n = 0L), "gear")$ll) == 0)
  chk("block CV: the self-check accepts a pointwise file written from the same draws and refuses one that is not",
      { lpd <- apply(tr$ll, 2, function(col) { m <- max(col); m + log(mean(exp(col - m))) })
        f1 <- tempfile(fileext = ".csv"); write.csv(data.frame(obs_index = 1:8, lpd = round(lpd, 4)), f1, row.names = FALSE)
        f2 <- tempfile(fileext = ".csv"); write.csv(data.frame(obs_index = 1:8, lpd = round(lpd + 0.01, 4)), f2, row.names = FALSE)
        f3 <- tempfile(fileext = ".csv"); write.csv(data.frame(obs_index = 1:7, lpd = round(lpd[1:7], 4)), f3, row.names = FALSE)
        isTRUE(bss_block_cv_check(tr$ll, f1)$ok) && isFALSE(bss_block_cv_check(tr$ll, f2)$ok) && isFALSE(bss_block_cv_check(tr$ll, f3)$ok) })
  chk("block CV: blocks are ISO weeks, Monday to Sunday (2024-12-30 to 2025-01-05 is one week)",
      identical(bss_block_weeks(as.Date(c("2024-12-30", "2025-01-05", "2025-01-06"))), c("2025-W01", "2025-W01", "2025-W02")))
  # the JOINT leave-out set: a second stream's observations of the week join the importance
  # ratios; checked against the exact conjugate leave-week-out predictive with both streams removed
  set.seed(2)
  yA <- rpois(60, 4); yB <- rpois(36, 10); bA <- rep(sprintf("W%02d", 1:12), each = 5); bB <- rep(sprintf("W%02d", 1:12), each = 3)
  lam2 <- rgamma(20000, a0 + sum(yA) + sum(yB), b0 + 60 + 2.5 * 36)
  llA <- sapply(yA, function(yi) dpois(yi, lam2, log = TRUE)); llB <- sapply(yB, function(yi) dpois(yi, 2.5 * lam2, log = TRUE))
  exj <- sapply(unique(bA), function(bk) { jA <- bA == bk; jB <- bB == bk; a2 <- a0 + sum(yA[!jA]) + sum(yB[!jB]); b2 <- b0 + sum(!jA) + 2.5 * sum(!jB)
    m <- sum(jA); sk <- sum(yA[jA]); lgamma(a2 + sk) - lgamma(a2) - sum(lgamma(yA[jA] + 1)) + a2 * log(b2) - (a2 + sk) * log(b2 + m) })
  tj <- bss_block_psis_loo(llA, bA, weight_ll = list(llB), weight_block = list(bB))
  chk("block CV: with a second stream in the leave-out set, PSIS equals the EXACT leave-week-out predictive with both streams removed (within 0.01 nats)",
      max(abs(tj$elpd_block - exj)) < 0.01 && all(tj$n_leaveout == 8L), sprintf("(max |diff| %.4f)", max(abs(tj$elpd_block - exj))))
  fab <- bss_block_cv_fit(list(A = list(ll = llA, block = bA), B = list(ll = llB, block = bB)))
  chk("block CV: bss_block_cv_fit() scores each stream under the joint leave-out set of all of them and labels it",
      isTRUE(all.equal(fab$A$elpd_block, tj$elpd_block)) && identical(fab$A$leaveout_streams[1], "A+B") && nrow(fab$B) == 12)
  chk("block CV: a block with a non-finite log-likelihood is NA with k = Inf and unreliable, the rest of the table intact",
      { llx <- llA; llx[1, 3] <- -Inf; tx <- bss_block_psis_loo(llx, bA)
        is.na(tx$elpd_block[1]) && is.infinite(tx$pareto_k[1]) && !tx$reliable[1] && all(is.finite(tx$elpd_block[-1])) })
  chk("block CV: with G = 2 the rebuild SUMS lambda_E_S over the gear groups (the gear-resolved model's likelihood), and a wrong D is refused",
      { nd2 <- 200; D2 <- 10; l2 <- array(exp(rnorm(nd2 * 2 * D2, log(10), 0.2)), c(nd2, 1, D2, 2))
        dr2 <- list(lambda_E_S = l2, r_E = rgamma(nd2, 20, 2), R_G = runif(nd2, 2, 3))
        sd2 <- list(D = D2, S = 1L, G = 2L, Gear_n = 5L, Gear_I = rpois(5, 50), day_Gear = 1:5, section_Gear = rep(1L, 5))
        g2 <- bss_effort_loglik(dr2, sd2, "gear")
        ref <- sapply(1:5, function(i) dnbinom(sd2$Gear_I[i], mu = (l2[, 1, i, 1] + l2[, 1, i, 2]) * dr2$R_G, size = dr2$r_E, log = TRUE))
        max(abs(g2$ll - ref)) < 1e-12 && inherits(try(bss_effort_loglik(dr2, modifyList(sd2, list(D = 11)), "gear"), silent = TRUE), "try-error") })
  # the post-hoc runner
  bf <- "06_diagnostics/run_marine_block_cv_2026-09-26.R"; bsrc <- readLines(bf, warn = FALSE)
  chk("block CV runner: exists, ships DRY_RUN <- TRUE, states rules R1 to R6 before the run, and can be pointed at a fixture",
      file.exists(bf) && any(grepl("^DRY_RUN <- TRUE", bsrc)) && all(vapply(sprintf("#   R%d.", 1:6), function(r) any(grepl(r, bsrc, fixed = TRUE)), logical(1))) &&
      any(grepl('Sys.getenv("MH_BLOCKCV_ROOT"', bsrc, fixed = TRUE)))
  chk("block CV runner: a fit is scored only after every gear / trailer stream matches its committed pointwise lpd and the OSP means match ppc_byobs (R1); an unevaluable comparison is REVIEW (R2)",
      any(grepl("if (!isTRUE(ck$ok)) { out$ok <- FALSE; break }", bsrc, fixed = TRUE)) &&
      any(grepl("no gear or trailer stream could be checked against a committed file (R1); nothing scored", bsrc, fixed = TRUE)) &&
      any(grepl("evaluable <- isTRUE(cmp$reliable_share >= MIN_REL) && isTRUE(cmp$n_used >= MIN_WKS)", bsrc, fixed = TRUE)) &&
      any(grepl("PSIS cannot carry this comparison", bsrc, fixed = TRUE)))
  chk("block CV runner: the leave-out set is the whole week's effort data across streams (bss_block_cv_fit), and a pair with a missing side is a REVIEW row, not a missing row",
      any(grepl("tabs <- tryCatch(bss_block_cv_fit(recs, k_max = K_MAX)", bsrc, fixed = TRUE)) &&
      any(grepl("not computable: a side has no block table", bsrc, fixed = TRUE)))
  chk("block CV runner: it never quits the session on a dry run (a sourced runner must not close RStudio)", !any(grepl("quit(", bsrc, fixed = TRUE)))

  # ---- (4) the documents ---------------------------------------------------------------
  cr <- rd("07_documentation/development_notes/CHANGE_REGISTER.md"); vc <- rd("07_documentation/development_notes/VALIDATION_CAMPAIGN.md")
  ps <- rd("07_documentation/development_notes/PIPELINE_STATUS.md"); md <- rd("07_documentation/BSS-GH-pooled-CPUE-model-documentation.md")
  chk("1x docs: the register carries B38, B39 and two 2026-09-26 defect rows, and A30 says RUN and NOT ADOPTED",
      grepl("| B38 |", cr, fixed = TRUE) && grepl("| B39 |", cr, fixed = TRUE) && length(gregexpr("| 2026-09-26 |", cr, fixed = TRUE)[[1]]) == 2 &&
      grepl("RUN 2026-09-25/26 (Section 1x); NOT ADOPTED under the pre-committed rule", cr, fixed = TRUE))
  chk("1x docs: the campaign has Section 1x with its git anchor, and the status document's box records the port jitter and B38",
      grepl("## 1x. The marine hazard ladder", vc, fixed = TRUE) && grepl("> | 1x | 2026-09-26 |", vc, fixed = TRUE) &&
      grepl("8b3f661", vc, fixed = TRUE) &&
      # 2026-09-28: the box moved to a post-B38 run, so it records B38 confirmed rather than the jitter it removed
      grepl("The port total is now reproducible to the crab (B38, confirmed in the field)", ps, fixed = TRUE))
  chk("1x docs: the method document's reproducibility section no longer says the port total resamples without qualification",
      grepl("The port total used to resample; since 2026-09-26 it does not.", md, fixed = TRUE) &&
      !grepl("**The port total resamples.**", md, fixed = TRUE))
  chk("1x docs: the register's A30 effect cell records the measured direction (UP) and corrects the earlier 'DOWN'",
      grepl("Measured 2026-09-26: UP, not down as first written here", cr, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 78. THE BLOCK CROSS-VALIDATION RAN: ITS EVIDENCE, THE DEFECT IT EXPOSED IN ITS OWN RUNNER,
#     THE UNDER-POWERED CLAUSE, AND THE JOINT EFFORT SCORE (2026-09-27; VALIDATION_CAMPAIGN
#     Section 1y, CHANGE_REGISTER B40, three C rows, D32).
#
#     Three things are pinned, as section 77 pins the ladder. (1) THE EVIDENCE, recomputed from
#     the committed per-week tables with the library itself and compared with the committed
#     pairs file: the OSP stream's +18.9 at 3.9 SE, the trailer's +7.9 at 1.2 SE, the two
#     streams together at 2.9 SE, the bar tick's nothing, M5 = M4, the shore not evaluable at
#     weekly AR, the winter's 13 uninformative trailer weeks, and the concentration of the boat
#     change in December to March. Skipped with a NOTE where the folders are absent. (2) THE
#     FIXES: identical block tables are INFO not FAIL (the committed FAIL row was 2e-14 nats),
#     the joint table equals the exact joint leave-week-out predictive of a conjugate model
#     with both streams removed and is NOT the sum of the per-stream rows, the per-week helper,
#     the runner's R3' / R7 text and its weeks file, the OSP pointwise LOO, the adequacy
#     aggregate kept comparable, the delivery files gone and ignored, the equivalence
#     declarations current. (3) THE DOCUMENTS.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  source("03_R_functions/bss_block_cv.R")
  MH <- c(M1 = "05_output/20260925/pooled-CPUE-MH-M1-off", M2 = "05_output/20260925/pooled-CPUE-MH-M2-sca",
          M3 = "05_output/20260925/pooled-CPUE-MH-M3-bar", M4 = "05_output/20260926/pooled-CPUE-MH-M4-both",
          M5 = "05_output/20260926/pooled-CPUE-MH-M5-auto")
  BAG <- "private_boat_all_gear_Dungeness_Kept"; BPC <- "private_boat_ring_net_only_Dungeness_Kept"
  SAG <- "shore_all_gear_Dungeness_Kept"; SPC <- "shore_ring_net_only_Dungeness_Kept"
  pf <- "05_output/marine_hazard_2026-09-26_blockcv_pairs.csv"; vf <- "05_output/marine_hazard_2026-09-26_blockcv_verdicts.csv"
  sf <- "05_output/marine_hazard_2026-09-26_blockcv.csv"
  chk("1y: the block-CV summary, pairs and verdicts CSVs are committed", all(file.exists(c(pf, vf, sf))))
  have_tabs <- all(dir.exists(MH)) && all(file.exists(file.path(MH, sprintf("loo_block_trailer_%s.csv", BAG))))
  if (all(file.exists(c(pf, vf, sf))) && have_tabs) {
    V <- read.csv(vf, stringsAsFactors = FALSE); P <- read.csv(pf, stringsAsFactors = FALSE); S <- read.csv(sf, stringsAsFactors = FALSE)
    r1 <- grepl("^R1:", V$criterion)
    # 20 fits at the 1y read (M1 to M5); the 2026-09-27 re-run added M6's four, whose boat rows are checked
    # EXACTLY against loo_pointwise_osp_* (section 79); the M1 to M5 boat rows remain approx-checked
    r1_pre <- r1 & V$stage %in% c("M1", "M2", "M3", "M4", "M5")
    chk("1y: R1 held on all 20 fits: every reconstruction row is PASS, gear / trailer to 5e-5, OSP approx-checked",
        sum(r1_pre) == 20 && all(V$verdict[r1] == "PASS") && all(grepl("MATCHES the committed pointwise file", V$observed[r1])) &&
        all(grepl("approx-checked", V$observed[r1_pre & grepl("boat", V$criterion)])))
    bt  <- function(rung, stream, fit) read.csv(file.path(MH[[rung]], sprintf("loo_block_%s_%s.csv", stream, fit)), stringsAsFactors = FALSE)
    cmp <- function(a, b, stream, fit) bss_block_cv_compare(bt(a, stream, fit), bt(b, stream, fit))
    prow <- function(pair, fit, stream) P[P$pair == pair & P$fit == fit & P$stream == stream, , drop = FALSE]
    o <- cmp("M1", "M2", "osp", BAG); t <- cmp("M1", "M2", "trailer", BAG)
    chk("1y: boat SCA on held-out OSP weeks recomputes to +18.9 at 4.8 SE (3.90 SE): 25 of 29 weeks reliable, 21 positive",
        abs(o$diff - 18.85) < 0.05 && abs(o$se - 4.83) < 0.05 && o$ratio > 2 && o$n_used == 25 && o$n_blocks == 29 && o$n_positive == 21,
        sprintf("(%+.2f, SE %.2f, %.2f SE)", o$diff, o$se, o$ratio))
    chk("1y: boat SCA on held-out trailer weeks recomputes to +7.9 at 6.6 SE (1.20 SE): 34 of 37 weeks reliable, 22 positive",
        abs(t$diff - 7.93) < 0.05 && abs(t$se - 6.63) < 0.05 && t$ratio > 0 && t$ratio < 2 && t$n_used == 34 && t$n_blocks == 37 && t$n_positive == 22,
        sprintf("(%+.2f, SE %.2f, %.2f SE)", t$diff, t$se, t$ratio))
    chk("1y: the committed pairs file carries the numbers the library recomputes (diff and SE to 1e-6) and the verdicts the rule gives (OSP PASS, trailer REVIEW)",
        { po <- prow("M2 vs M1", "boat_all_gear", "osp"); pt <- prow("M2 vs M1", "boat_all_gear", "trailer")
          nrow(po) == 1 && nrow(pt) == 1 && abs(po$diff - o$diff) < 1e-6 && abs(po$se - o$se) < 1e-6 &&
          abs(pt$diff - t$diff) < 1e-6 && abs(pt$se - t$se) < 1e-6 && po$verdict == "PASS" && pt$verdict == "REVIEW" })
    # the two streams together, from the per-week tables: sum of the streams' weekly differences over the
    # weeks reliable in every stream present (the run had no joint row; B40 adds the exact one)
    wk <- function(a, b, stream, fit) bss_block_cv_weeks(bt(a, stream, fit), bt(b, stream, fit))
    wt <- wk("M1", "M2", "trailer", BAG); wo <- wk("M1", "M2", "osp", BAG)
    # every week either stream observes (42), kept when reliable in every stream present (38)
    both <- merge(wt, wo, by = "block", all = TRUE, suffixes = c("_t", "_o"))
    both$ok <- (is.na(both$used_t) | both$used_t) & (is.na(both$used_o) | both$used_o)
    both$d  <- ifelse(is.na(both$diff_t), 0, both$diff_t) + ifelse(is.na(both$diff_o), 0, both$diff_o)
    d <- both$d[both$ok]
    chk("1y: both boat streams together on held-out weeks: +26.8 at 9.4 SE (2.9 SE) over 38 of 42 weeks, 27 positive",
        nrow(both) == 42 && length(d) == 38 && abs(sum(d) - 26.79) < 0.1 && abs(sqrt(38) * sd(d) - 9.36) < 0.1 && sum(d) / (sqrt(38) * sd(d)) > 2 && sum(d > 0) == 27,
        sprintf("(%+.2f, SE %.2f over %d weeks)", sum(d), sqrt(length(d)) * sd(d), length(d)))
    wtr_osp <- wt[wt$used & wt$block %in% wo$block, ]
    chk("1y: over the weeks the OSP stream also covers, the trailer scores 1.6 SE (+8.5 at 5.3) where the OSP scores 3.9",
        nrow(wtr_osp) == 21 && abs(sum(wtr_osp$diff) - 8.51) < 0.05 && abs(sum(wtr_osp$diff) / (sqrt(21) * sd(wtr_osp$diff)) - 1.59) < 0.02)
    wtr <- wt[wt$used, ]; win <- wtr[grepl("^2024-W|^2025-W0[1-9]$", wtr$block), ]
    chk("1y: the winter (2024-W49 to 2025-W09) has trailer counts only and its 13 held-out weeks are uninformative: -0.6 at 3.9 SE",
        nrow(win) == 13 && abs(sum(win$diff) + 0.57) < 0.05 && abs(sqrt(13) * sd(win$diff) - 3.88) < 0.1 &&
        !any(grepl("^2024-W|^2025-W0[1-9]$", wo$block)), sprintf("(%+.2f, SE %.2f)", sum(win$diff), sqrt(13) * sd(win$diff)))
    b42 <- cmp("M2", "M4", "trailer", BAG); b42o <- cmp("M2", "M4", "osp", BAG)
    chk("1y: the bar tick beyond the archive adds nothing on held-out weeks (trailer -2.0, OSP -0.4, neither beyond 2 SE)",
        abs(b42$diff + 1.99) < 0.05 && abs(b42o$diff + 0.39) < 0.05 && abs(b42$ratio) < 2 && abs(b42o$ratio) < 2)
    chk("1y: M5's boat block tables are M4's (bit-identical boat fits)",
        isTRUE(all.equal(bt("M4", "trailer", BAG)$elpd_block, bt("M5", "trailer", BAG)$elpd_block)) &&
        isTRUE(all.equal(bt("M4", "osp", BAG)$elpd_block, bt("M5", "osp", BAG)$elpd_block)))
    chk("1y: the shore all-gear fits are NOT evaluable at weekly AR (16 of 38 weeks reliable in M1, 12 in M2; median k above 0.7) and the boat is (34 of 37)",
        { s1 <- S[S$rung == "M1" & S$fit == "shore_all_gear", ]; s2 <- S[S$rung == "M2" & S$fit == "shore_all_gear", ]; b1 <- S[S$rung == "M1" & S$fit == "boat_all_gear" & S$stream == "trailer", ]
          s1$n_reliable == 16 && s2$n_reliable == 12 && s1$n_weeks == 38 && b1$n_reliable == 34 && median(bt("M1", "gear", SAG)$pareto_k) > 0.7 })
    pc <- cmp("M1", "M2", "trailer", BPC)
    chk("1y: the pot-closure boat trailer stream cleared on 8 of 11 weeks (+4.1 at 1.5 SE, 2.66 SE), 7 of 8 positive",
        pc$n_used == 8 && pc$n_blocks == 11 && abs(pc$diff - 4.12) < 0.05 && abs(pc$se - 1.55) < 0.05 && pc$ratio > 2 && pc$n_positive == 7)
    # The committed 2026-09-26 file carries that row as FAIL (the defect); a re-source under B40
    # rewrites it as INFO (R7). Either state is accepted; a PASS or REVIEW there is not.
    chk("1y fix: the identical-fit pair M4 vs M2 on the shore pot-closure fit reads FAIL (the 2026-09-26 defect) or INFO (after B40), never a verdict; the library reads it as identical",
        { rw <- V[grepl("bar restriction beyond the archive (M4 vs M2), shore_pot_closure gear stream", V$criterion, fixed = TRUE), ]
          pr <- prow("M4 vs M2", "shore_pot_closure", "gear")
          nrow(rw) == 1 && rw$verdict %in% c("FAIL", "INFO") && sum(V$verdict == "FAIL") <= 1 &&
          (abs(pr$diff) < 1e-8 || isTRUE(as.logical(pr$identical))) &&
          isTRUE(cmp("M2", "M4", "gear", SPC)$identical) && is.na(cmp("M2", "M4", "gear", SPC)$ratio) })
    chk("1y fix: every same-fit pair among the committed shore tables is identical (M3, M5 = M1; M4 = M2) and the real pair (M2 vs M1) is not",
        isTRUE(cmp("M1", "M3", "gear", SPC)$identical) && isTRUE(cmp("M1", "M5", "gear", SPC)$identical) && isTRUE(cmp("M1", "M3", "gear", SAG)$identical) &&
        isTRUE(cmp("M2", "M4", "gear", SAG)$identical) && !isTRUE(cmp("M1", "M2", "gear", SPC)$identical) && !isTRUE(cmp("M1", "M2", "gear", SAG)$identical))
    dc <- function(rung) { x <- read.csv(file.path(MH[[rung]], sprintf("bss_daily_catch_%s.csv", BAG)), stringsAsFactors = FALSE); tapply(x$median, substr(x$event_date, 1, 7), sum) }
    c1 <- dc("M1"); c2 <- dc("M2"); dlt <- c2 - c1[names(c2)]
    chk("1y: the boat all-gear change under SCA is concentrated in the winter: January +29%, December to March 57% of the change on 11% of the catch",
        abs(c2["2025-01"] / c1["2025-01"] - 1.286) < 0.01 &&
        { w <- c("2024-12", "2025-01", "2025-02", "2025-03"); sh <- sum(dlt[w]) / sum(dlt); sh > 0.5 && sh < 0.65 && sum(c1[w]) / sum(c1) < 0.15 })
  } else cat("NOTE  1y: the marine rung folders or the block-CV CSVs are absent (sparse checkout); the evidence checks are skipped\n")

  # ---- (2) the code -------------------------------------------------------------------
  set.seed(1)
  a0 <- 2; b0 <- 0.5
  yA <- rpois(60, 4); yB <- rpois(36, 10); bA <- rep(sprintf("W%02d", 1:12), each = 5); bB <- rep(sprintf("W%02d", 1:12), each = 3)
  # 40,000 draws: a two-stream week is a larger leave-out set than section 77's one-stream score,
  # and at 20,000 draws the PSIS error on one seed reached 0.015 nats (0.001 to 0.005 on others)
  lam2 <- rgamma(40000, a0 + sum(yA) + sum(yB), b0 + 60 + 2.5 * 36)
  llA <- sapply(yA, function(yi) dpois(yi, lam2, log = TRUE)); llB <- sapply(yB, function(yi) dpois(yi, 2.5 * lam2, log = TRUE))
  exboth <- sapply(unique(bA), function(bk) { jA <- bA == bk; jB <- bB == bk
    a2 <- a0 + sum(yA[!jA]) + sum(yB[!jB]); b2 <- b0 + sum(!jA) + 2.5 * sum(!jB)
    m <- sum(jA) + 2.5 * sum(jB); sk <- sum(yA[jA]) + sum(yB[jB]); skB <- sum(yB[jB])
    lgamma(a2 + sk) - lgamma(a2) - sum(lgamma(yA[jA] + 1)) - sum(lgamma(yB[jB] + 1)) + skB * log(2.5) + a2 * log(b2) - (a2 + sk) * log(b2 + m) })
  fab <- bss_block_cv_fit(list(A = list(ll = llA, block = bA), B = list(ll = llB, block = bB)))
  tj <- fab$joint
  chk("B40: bss_block_cv_fit() adds a JOINT table for a fit with more than one stream: one row per week, every observation of both streams held out and scored",
      !is.null(tj) && nrow(tj) == 12 && all(tj$n_obs == 8L) && all(tj$n_leaveout == 8L) && identical(tj$leaveout_streams[1], "A+B") && identical(tj$block, sort(tj$block)))
  chk("B40: the joint table equals the EXACT joint leave-week-out predictive of both streams' counts (within 0.01 nats)",
      max(abs(tj$elpd_block - exboth[tj$block])) < 0.01, sprintf("(max |diff| %.4f)", max(abs(tj$elpd_block - exboth[tj$block]))))
  chk("B40: the joint row is NOT the sum of the per-stream rows (a covariance term separates them), which is why it exists",
      { sa <- fab$A$elpd_block[match(tj$block, fab$A$block)] + fab$B$elpd_block[match(tj$block, fab$B$block)]
        mean(abs(tj$elpd_block - sa)) > 1e-3 })
  chk("B40: a single-stream fit gets no joint table", is.null(bss_block_cv_fit(list(A = list(ll = llA, block = bA)))$joint))
  ta <- fab$A; tb <- ta; tb$elpd_block <- tb$elpd_block + rnorm(nrow(tb), 0, 1e-14)
  ci <- bss_block_cv_compare(ta, tb)
  chk("B40: two block tables that differ by summation noise are IDENTICAL: diff 0, ratio NA, and the one-line reading says so",
      isTRUE(ci$identical) && ci$diff == 0 && is.na(ci$ratio) && ci$max_abs_diff < 1e-12 && grepl("identical block tables", bss_block_cv_str(ci, "M1", "M3"), fixed = TRUE))
  tc <- ta; tc$elpd_block <- tc$elpd_block + 0.3
  cr3 <- bss_block_cv_compare(ta, tc)
  chk("B40: a real shift is not identical, and the comparison reports the count of positive weeks and the median",
      !isTRUE(cr3$identical) && abs(cr3$diff - 3.6) < 1e-9 && cr3$n_positive == 12L && abs(cr3$median_diff - 0.3) < 1e-9 && abs(cr3$max_abs_diff - 0.3) < 1e-9)
  wkx <- bss_block_cv_weeks(ta, tc)
  chk("B40: bss_block_cv_weeks() returns the per-week rows behind a comparison with the reliability flag",
      nrow(wkx) == 12 && all(c("block", "elpd_block_a", "elpd_block_b", "diff", "pareto_k_a", "pareto_k_b", "used") %in% names(wkx)) &&
      all(abs(wkx$diff - 0.3) < 1e-9) && all(wkx$used))
  bf <- "06_diagnostics/run_marine_block_cv_2026-09-26.R"; bsrc <- readLines(bf, warn = FALSE)
  chk("B40 runner: states R3' (the joint row is the statistic) and R7 (identical fits are INFO), scores the joint table, and writes the per-week file",
      any(grepl("#   R3'.", bsrc, fixed = TRUE)) && any(grepl("#   R7.", bsrc, fixed = TRUE)) &&
      any(grepl('SCORED_OF <- lapply(STREAMS_OF, function(s) if (length(s) > 1L) c(s, "joint") else s)', bsrc, fixed = TRUE)) &&
      any(grepl('verdict <- if (isTRUE(cmp$identical)) "INFO"', bsrc, fixed = TRUE)) &&
      any(grepl("marine_hazard_2026-09-26_blockcv_weeks.csv", bsrc, fixed = TRUE)) && any(grepl("primary = primary", bsrc, fixed = TRUE)))
  chk("B40 runner: the joint row of a joint table is what R3 is applied to (primary), and a per-stream row of a multi-stream fit is not",
      any(grepl('primary <- identical(sn, "joint") || length(SCORED_OF[[fk]]) == 1L', bsrc, fixed = TRUE)))
  chk("B40 runner: its verdict count line reports INFO too", any(grepl('%d PASS  %d FAIL  %d REVIEW  %d INFO', bsrc, fixed = TRUE)))
  srd <- rd("03_R_functions/save_run_diagnostics.R"); bma <- rd("03_R_functions/bss_model_adequacy.R")
  chk("B40: write_loo_diagnostics() writes the OSP stream's pointwise LOO from log_lik_osp, and the adequacy aggregate stays on gear / trailer / catch",
      grepl('osp     = list(par = "log_lik_osp",     n = stan_data$OSP_n %||% 0,  days = stan_data$day_OSP,   y = stan_data$OSP_I)', srd, fixed = TRUE) &&
      grepl('loo <- loo[as.character(loo$stream) %in% c("gear", "trailer", "catch"), , drop = FALSE]', bma, fixed = TRUE))
  chk("B40: the adequacy core ignores an OSP row in loo_summary (comparability with every committed model_adequacy.csv)",
      { e <- new.env(); sys.source("03_R_functions/bss_model_adequacy.R", envir = e)
        loo3 <- data.frame(stream = c("gear", "trailer", "catch"), n_obs = c(100, 100, 100), p_loo = c(5, 5, 5), n_pareto_k_gt_0.7 = c(0L, 1L, 0L), stringsAsFactors = FALSE)
        loo4 <- rbind(loo3, data.frame(stream = "osp", n_obs = 10, p_loo = 9, n_pareto_k_gt_0.7 = 7L))
        r3 <- e$.bma_core("x", loo3, NULL, NULL); r4 <- e$.bma_core("x", loo4, NULL, NULL)
        identical(r3$p_loo_frac, r4$p_loo_frac) && identical(r3$n_pareto_bad, r4$n_pareto_bad) && r4$n_pareto_bad == 1L && r4$p_loo_frac == 0.05 })
  gi <- readLines(".gitignore", warn = FALSE)
  chk("B40 hygiene: the delivery tarball and README folder swept in by c7e8cd5 are gone, and .gitignore excludes them",
      !file.exists("marine-hazard-results-patches-2026-09-26.tar.gz") && !dir.exists("marine-hazard-results-patches-2026-09-26") &&
      any(gi == "*.tar.gz") && any(gi == "marine-hazard-*-patches-*/"))
  mhr <- readLines("06_diagnostics/run_marine_hazard_batch_2026-09-25.R", warn = FALSE)
  chk("B40: the ladder runner's superseded equivalence declaration says SUPERSEDED and a 2026-09-27 one follows it",
      any(grepl("SUPERSEDED 2026-09-27 by the entry below", mhr, fixed = TRUE)) && any(grepl("plus the 2026-09-27 block-CV results patch", mhr, fixed = TRUE)))

  # ---- (2b) B41: the season-split candidate and rung M6 ----------------------------------
  em <- new.env(); em$`%||%` <- function(a, b) if (is.null(a)) b else a
  sys.source("03_R_functions/bss_marine_hazard_covariates.R", envir = em)
  chk("B41: the module offers the split pair as labelled candidates and names it as one definition",
      all(c("nws_sca_any_winter", "nws_sca_any_rest") %in% names(em$.mh_labels)) && identical(em$.mh_split, c("nws_sca_any_winter", "nws_sca_any_rest")) &&
      identical(em$.mh_defaults(list())$winter, c(12L, 1L, 2L)) && identical(em$.mh_defaults(list(marine_hazard_winter_months = c(11, 12, 1, 2, 3)))$winter, c(11, 12, 1, 2, 3)))
  # a synthetic archive: one SCA event a week for a year, so the flag series can be built without the workbook
  ev <- data.frame(ugc = "PZZ110", ps = "SC.Y", in_effect = TRUE,
                   start = as.POSIXct(sprintf("%s 12:00:00", seq(as.Date("2024-09-16"), as.Date("2025-09-15"), by = "7 days")), tz = "UTC"))
  ev$end <- ev$start + 6 * 3600
  attr(ev, "coverage") <- data.frame(ugc = c("PZZ110", "PZZ156"), pull_start = as.Date("2024-01-01"), pull_end = as.Date("2026-01-01"))
  pfl <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15")
  fl <- em$marine_hazard_flag_series(pfl, events = ev)
  chk("B41: the split columns are the any-zone flag on winter days and on the other days, and sum to it (default winter = Dec, Jan, Feb)",
      all(c("nws_sca_any_winter", "nws_sca_any_rest") %in% names(fl)) && all(fl$nws_sca_any_winter + fl$nws_sca_any_rest == fl$nws_sca_any) &&
      all(fl$nws_sca_any_winter[!format(fl$event_date, "%m") %in% c("12", "01", "02")] == 0) &&
      all(fl$nws_sca_any_rest[format(fl$event_date, "%m") %in% c("12", "01", "02")] == 0) && sum(fl$nws_sca_any) == 53)
  chk("B41: the winter months are configurable, and a month outside 1..12 is refused",
      { f2 <- em$marine_hazard_flag_series(c(pfl, list(marine_hazard_winter_months = c(11, 12, 1, 2, 3))), events = ev)
        sum(f2$nws_sca_any_winter) > sum(fl$nws_sca_any_winter) &&
        inherits(try(em$marine_hazard_flag_series(c(pfl, list(marine_hazard_winter_months = c(0, 13))), events = ev), silent = TRUE), "try-error") })
  psel <- list(marine_hazard_mode = "manual", marine_hazard_manual_shore = character(0),
               marine_hazard_manual_boat = c("nws_sca_any_winter", "nws_sca_any_rest"),
               marine_hazard_candidates_boat = c("nws_sca_any", "bar_restriction", "nws_sca_any_winter", "nws_sca_any_rest"))
  s6 <- em$marine_hazard_select(NULL, psel)
  s7 <- em$marine_hazard_select(NULL, modifyList(psel, list(marine_hazard_manual_boat = c("nws_sca_any", "nws_sca_any_winter", "nws_sca_any_rest"))))
  s8 <- em$marine_hazard_select(NULL, modifyList(psel, list(marine_hazard_manual_boat = c("nws_sca_any_winter", "nws_sca_any_rest"),
                                                           marine_hazard_candidates_boat = c("nws_sca_any", "bar_restriction"))))
  chk("B41: under manual the pair enters together as one definition; named beside nws_sca_any one definition is kept; not offered, it is ignored with a note",
      setequal(s6$private_boat, c("nws_sca_any_winter", "nws_sca_any_rest")) && length(s7$private_boat) == 1 && startsWith(s7$private_boat, "nws_") &&
      !length(s8$private_boat) && any(grepl("non-candidate", s8$note)))
  chk("B41: the pre-B41 configurations select as before (nws_sca_any alone; nws_sca_any with the bar tick; the shore's nws_sca_any)",
      { a <- em$marine_hazard_select(NULL, list(marine_hazard_mode = "manual", marine_hazard_manual_shore = "nws_sca_any", marine_hazard_manual_boat = c("nws_sca_any", "bar_restriction")))
        identical(a$shore, "nws_sca_any") && setequal(a$private_boat, c("nws_sca_any", "bar_restriction")) })
  eo <- new.env(); eo$`%||%` <- function(a, b) if (is.null(a)) b else a; sys.source("03_R_functions/bss_opener_covariates.R", envir = eo)
  chk("B41: in a September-to-November window the winter column has no flagged day and the design-matrix guard drops it, keeping the rest term",
      { days_pc <- data.frame(event_date = seq(as.Date("2024-09-16"), as.Date("2024-11-30"), by = "day"))
        sp <- eo$opener_design_matrix(days_pc, character(0), fl, list(opener_min_days = 10), extra = c("nws_sca_any_winter", "nws_sca_any_rest"))
        identical(sp$labels, "nws_sca_any_rest") && any(grepl("nws_sca_any_winter", sp$dropped)) })
  mhr <- readLines("06_diagnostics/run_marine_hazard_batch_2026-09-25.R", warn = FALSE)
  chk("B41 ladder: rung M6 names the pair on the boat with no shore term, widens the candidate list in its own delta, states rule 9, and declared keys are per stage",
      any(grepl('M6 = list(tag = "MH-M6-split"', mhr, fixed = TRUE)) && any(grepl('marine_hazard_manual_boat = c("nws_sca_any_winter", "nws_sca_any_rest")', mhr, fixed = TRUE)) &&
      any(grepl("#   9. THE SEASON SPLIT (M6", mhr, fixed = TRUE)) && any(grepl("declared_keys <- function(sid) unique(c(DELTA_KEYS, names(STAGE_DEFS[[sid]]$delta %||% list())))", mhr, fixed = TRUE)) &&
      any(grepl("keys <- sort(unique(c(declared_keys(sid), names(WINDOW))))", mhr, fixed = TRUE)) && any(grepl('list(sid = "M6", ctl = "M2", fit = FIT_BOAT,  cov = "nws_sca_any_winter")', mhr, fixed = TRUE)) &&
      any(grepl("^verdict_M6 <- function", mhr)) && any(grepl('STAGES  <- c("M0", "M1", "M2", "M3", "M4", "M5", "M6")', mhr, fixed = TRUE)) && any(grepl("^DRY_RUN <- TRUE", mhr)))
  # the rendered rungs' digests must be what their folders recorded, or RESUME would refit them
  # (M1 to M5 at the 1y read; M6 since its 2026-09-27 render)
  MH6 <- c(MH, M6 = "05_output/20260927/pooled-CPUE-MH-M6-split")
  if (all(dir.exists(MH))) {
    dg <- local({
      e <- new.env(); e$.here <- function(...) file.path(getwd(), ...)
      e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a   # the runner's own definition
      source("run_config.R", local = e)
      lift <- function(pat, upto) { i <- grep(pat, mhr); j <- i; while (!grepl(upto, mhr[j])) j <- j + 1L; eval(parse(text = mhr[i:j]), envir = e) }
      eval(parse(text = mhr[grep("^BASE <- ", mhr)]), envir = e)
      lift("^WINDOW <- list\\(", "^\\)"); lift("^STAGE_DEFS <- list\\(", "^\\)"); lift("^DELTA_KEYS <- ", "^DELTA_KEYS")
      lift("^declared_keys <- function", "^declared_keys"); lift("^resolve_cfg <- function", "^\\}")
      lift("^digest_or_hash <- function", "^\\}"); lift("^stage_digest <- function", "^\\}")
      structure(vapply(names(MH6)[dir.exists(MH6)], function(sid) {
        rec <- readLines(file.path(MH6[[sid]], "MH_STAGE.txt"), warn = FALSE); rec <- sub("^digest: ", "", rec[grepl("^digest: ", rec)])
        identical(rec, e$stage_digest(sid)) }, logical(1)), m6 = e$stage_digest("M6")) })
    chk("B41 ladder: adding M6 left the five rendered rungs' stage digests exactly as their folders recorded (RESUME reads them back; only M6 renders)", all(dg[names(dg) != "M6"]),
        sprintf("(%s)", paste(names(dg)[!dg], collapse = ",")))
    chk("A30 ladders: M6's own digest (b1fb62cd) reproduces too, so all six marine rungs read back under the adopted run_config.R",
        identical(attr(dg, "m6"), "b1fb62cd") && (!"M6" %in% names(dg) || isTRUE(dg[["M6"]])))
  } else cat("NOTE  B41: rung folders absent; the digest check is skipped\n")
  chk("B41 block-CV runner: lists M6 and the pair M6 vs M2 on the boat fits only, and checks an OSP stream exactly when its pointwise file exists",
      any(grepl('M6 = "MH-M6-split"', bsrc, fixed = TRUE)) && any(grepl('list(b = "M6", a = "M2", what = "SCA split by season beyond the constant term', bsrc, fixed = TRUE)) &&
      any(grepl('fits = c("boat_all_gear", "boat_pot_closure")', bsrc, fixed = TRUE)) && any(grepl("if (!is.null(p$fits) && !fk %in% p$fits) next", bsrc, fixed = TRUE)) &&
      any(grepl('if (sn %in% c("gear", "trailer") || file.exists(pw)) {', bsrc, fixed = TRUE)))
  rcfg <- readLines("run_config.R", warn = FALSE)
  chk("B41 config: marine_hazard_winter_months ships c(12, 1, 2) with the reasoning beside it", any(grepl("^\\s*marine_hazard_winter_months\\s*=\\s*c\\(12, 1, 2\\),", rcfg)))
  chk("B41 ladder: rule 9 is stated in order with its three branches, names the boat all-gear fit, and calls (a) the likely outcome before the render; M0 proves M6 vs M2; the recommendation reads M6 by rule 9",
      any(grepl("#      (a) The winter coefficient is NOT identified", mhr, fixed = TRUE)) && any(grepl("#      (b) Both identified, and they DIFFER", mhr, fixed = TRUE)) &&
      any(grepl("#      (c) Both identified, and they do NOT differ", mhr, fixed = TRUE)) && any(grepl("this is the LIKELY outcome", mhr, fixed = TRUE)) &&
      any(grepl("the fit judged is the BOAT", mhr, fixed = TRUE)) && any(grepl("M6's Stan data differs from M2's in K_open and X_open_flat ONLY, and its two columns sum to M2's one", mhr, fixed = TRUE)) &&
      any(grepl("^  judge_M6 <- function", mhr)) && any(grepl('rule 9: which branch the season split lands in', mhr, fixed = TRUE)))
  chk("B41 module: when the pair is kept the 'one NWS definition' reason names both members, and the winter-month validation refuses fractions, duplicates and logicals",
      { s7b <- em$marine_hazard_select(NULL, modifyList(psel, list(marine_hazard_manual_boat = c("nws_sca_any_winter", "nws_sca_any_rest", "nws_sca_bar"),
                                                                    marine_hazard_candidates_boat = c("nws_sca_any", "nws_sca_bar", "bar_restriction", "nws_sca_any_winter", "nws_sca_any_rest"))))
        # the pair is named first, so with no p-values the first family (the pair) is kept and nws_sca_bar dropped
        setequal(s7b$private_boat, c("nws_sca_any_winter", "nws_sca_any_rest")) && any(grepl("kept nws_sca_any_winter + nws_sca_any_rest", s7b$table$reason, fixed = TRUE)) &&
        inherits(try(em$marine_hazard_flag_series(c(pfl, list(marine_hazard_winter_months = c(12.5, 1))), events = ev), silent = TRUE), "try-error") &&
        inherits(try(em$marine_hazard_flag_series(c(pfl, list(marine_hazard_winter_months = c(12, 12))), events = ev), silent = TRUE), "try-error") &&
        inherits(try(em$marine_hazard_flag_series(c(pfl, list(marine_hazard_winter_months = c(TRUE, FALSE))), events = ev), silent = TRUE), "try-error") })
  # the four-season desk screen: its script, its committed table, and the numbers the documents quote
  dsf <- "06_diagnostics/desk_sca_season_split_2026-09-27.R"; dsc <- "05_output/marine_hazard_2026-09-27_season_split_screen.csv"
  chk("B41 desk screen: the script exists, states what it is and is not, and its table is committed",
      file.exists(dsf) && any(grepl("A DESK SCREEN, SECONDS, NO MCMC", readLines(dsf, warn = FALSE), fixed = TRUE)) && file.exists(dsc))
  if (file.exists(dsc)) {
    ds <- read.csv(dsc, stringsAsFactors = FALSE)
    g <- function(model, fam, seasons, term) ds[ds$model == model & ds$family == fam & ds$seasons == seasons & ds$term == term, , drop = FALSE]
    nbw <- g("split", "negbin", "all", "nws_sca_any_winter"); nbr <- g("split", "negbin", "all", "nws_sca_any_rest")
    nbi <- g("interaction", "negbin", "all", "sca:winter"); qpi <- g("interaction", "quasipoisson", "all", "sca:winter")
    chk("B41 desk screen: four seasons, 794 sampled days, 315 advisory days, 129 of them winter; winter 0.43 [0.29, 0.62] against rest 0.22 [0.18, 0.27]",
        nrow(nbw) == 1 && nbw$n_days == 794 && nbw$n_advisory_days == 129 && g("constant", "negbin", "all", "sca")$n_advisory_days == 315 &&
        abs(nbw$rate_ratio - 0.426) < 0.005 && abs(nbw$lo95 - 0.290) < 0.005 && abs(nbw$hi95 - 0.624) < 0.005 &&
        abs(nbr$rate_ratio - 0.217) < 0.005 && abs(nbr$lo95 - 0.175) < 0.005 && abs(nbr$hi95 - 0.269) < 0.005)
    chk("B41 desk screen: the interaction (the test) is 1.96 at 3.0 SE under the negative binomial and 1.97 at 1.4 SE under the quasi-Poisson, and the split lowers the NB AIC by 7",
        abs(nbi$rate_ratio - 1.963) < 0.005 && abs(nbi$z - 3.01) < 0.02 && abs(qpi$rate_ratio - 1.969) < 0.005 && abs(qpi$z - 1.36) < 0.02 &&
        abs(g("constant", "negbin", "all", "sca")$aic - nbi$aic - 7.0) < 0.1)
    chk("B41 desk screen: the winter effect is weaker than the rest in every season that has a winter (2023-24, 2024-25, 2025-26), both families",
        all(vapply(c("2023-24", "2024-25", "2025-26"), function(se) all(vapply(c("negbin", "quasipoisson"), function(fm)
          g("split", fm, se, "nws_sca_any_winter")$rate_ratio > g("split", fm, se, "nws_sca_any_rest")$rate_ratio, logical(1))), logical(1))))
  }

  # ---- (3) the documents ---------------------------------------------------------------
  cr <- rd("07_documentation/development_notes/CHANGE_REGISTER.md"); vc <- rd("07_documentation/development_notes/VALIDATION_CAMPAIGN.md")
  ps <- rd("07_documentation/development_notes/PIPELINE_STATUS.md"); dr <- rd("06_diagnostics/README.md"); rc <- rd("run_config.R")
  dn <- rd("07_documentation/development_notes/marine-hazard-covariates-2026-09-25.md")
  chk("1y docs: the register carries B40, B41, D32 and five 2026-09-27 defect rows (the fifth the improvement ladder's hashed label); B39 RAN; D31 RUN; A30 records the road from DECISION PENDING to ADOPTED",
      grepl("| B40 |", cr, fixed = TRUE) && grepl("| B41 |", cr, fixed = TRUE) && grepl("| D32 |", cr, fixed = TRUE) && length(gregexpr("| 2026-09-27 |", cr, fixed = TRUE)[[1]]) == 5 &&
      grepl("| **RAN 2026-09-26** (`c7e8cd5`", cr, fixed = TRUE) && grepl("**RUN 2026-09-26 (B39; Section 1y)", cr, fixed = TRUE) &&
      grepl("M6 RENDERED 2026-09-27 (Section 1z): rule 9 branch (a)", cr, fixed = TRUE) && grepl("ADOPTED 2026-09-27 (Matt's decision, Section 1z)", cr, fixed = TRUE))
  chk("1y docs: the campaign has Section 1y with its git anchor naming the run commit, the under-powered clause and the three readings",
      grepl("## 1y. The block cross-validation, run and read: the boat SCA term interpolates where the data can test it", vc, fixed = TRUE) && grepl("> | 1y | 2026-09-27 |", vc, fixed = TRUE) &&
      grepl("c7e8cd5", vc, fixed = TRUE) && grepl("### 1y.4 Two defects in the runner and one in the rule", vc, fixed = TRUE) &&
      grepl("desk_sca_season_split_2026-09-27.R", vc, fixed = TRUE) && grepl("0.43 [0.29, 0.62]", vc, fixed = TRUE) &&
      grepl("### 1y.5 The decision this leaves with Matt", vc, fixed = TRUE))
  chk("1y docs: the status document, the diagnostics README, run_config.R and the design note carry the block-CV outcome",
      grepl("ADOPTED 2026-09-27; A30; RENDERED 2026-09-28, `1d3409d` (superseded as the authoritative run the same day by B50", ps, fixed = TRUE) && grepl("**RAN 2026-09-26** (`c7e8cd5`", dr, fixed = TRUE) &&
      grepl("THAT RAN 2026-09-26 (Section 1y)", rc, fixed = TRUE) && grepl("## 11. What the block cross-validation said", dn, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 79. A30 ADOPTED: THE BOAT ALL-GEAR ADVISORY TERM IS THE METHOD OF RECORD (2026-09-27;
#     VALIDATION_CAMPAIGN Section 1z, CHANGE_REGISTER A30 / B42, one C row).
#
#     Pinned: (1) THE RUNS the decision rests on, from the committed CSVs: M6's two coefficients
#     and its rule-9 branch (a), the exact joint block-CV rows (M2 vs M1 PASS at 2.86 SE; M6 vs
#     M2 no difference; the pot-closure M6 identical to M2), the first exact OSP reconstruction
#     check, the desk screen's table unchanged. (2) THE ADOPTION IN CODE: run_config's method
#     keys in section 2.10 (asserted in the "shipped" block above), marine_hazard_terms_for()
#     confining a selected term to the configured sub-seasons and both preps calling it, an
#     unknown regime refused, the three live ladders setting their own marine keys so their
#     rungs stay reproducible (the marine digests are checked in section 78), the marine
#     ladder's M0 row accepting the adopted method, DRY_RUN TRUE on both marine runners after
#     the second slip, the delivery folder gone and every such name ignored. (3) THE DOCUMENTS
#     moving together: the method document, the box, the register, the campaign, the guide.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  # ---- (1) the runs ------------------------------------------------------------------------
  lad_f <- "05_output/marine_hazard_2026-09-25_ladder.csv"; ver_f <- "05_output/marine_hazard_2026-09-25_verdicts.csv"
  pf <- "05_output/marine_hazard_2026-09-26_blockcv_pairs.csv"; bv <- "05_output/marine_hazard_2026-09-26_blockcv_verdicts.csv"
  M6 <- "05_output/20260927/pooled-CPUE-MH-M6-split"
  if (all(file.exists(c(lad_f, ver_f, pf, bv))) && dir.exists(M6)) {
    lad <- read.csv(lad_f, stringsAsFactors = FALSE); rownames(lad) <- lad$rung
    chk("1z: M6 rendered: the winter coefficient straddles zero (-0.59 [-1.24, +0.06]), the rest is identified (-1.26 [-1.56, -0.95]), and the boat total sits where M2 put it (47,290 vs 47,319)",
        "M6" %in% lad$rung && abs(lad["M6", "B_sca_boat_winter"] + 0.593) < 0.005 && lad["M6", "B_sca_boat_winter_hi"] > 0 && lad["M6", "B_sca_boat_winter_lo"] < 0 &&
        abs(lad["M6", "B_sca_boat_rest"] + 1.256) < 0.005 && lad["M6", "B_sca_boat_rest_hi"] < 0 &&
        abs(lad["M6", "boat_ag"] / lad["M2", "boat_ag"] - 1) < 0.002 && lad["M6", "boat_pc"] == lad["M2", "boat_pc"] && lad["M6", "shore_ag"] == lad["M1", "shore_ag"])
    V <- read.csv(ver_f, stringsAsFactors = FALSE)
    r9 <- V[grepl("^rule 9: which branch", V$criterion) & V$stage == "M6", ]
    chk("1z: the ladder's rule-9 row reads branch (a), as the rule said in advance it probably would, and the M6 recommendation rows carry it",
        nrow(r9) == 1 && startsWith(r9$observed, "(a)") && r9$verdict == "INFO" &&
        { rec <- read.csv("05_output/marine_hazard_2026-09-25_recommendation.csv", stringsAsFactors = FALSE)
          all(grepl("rule 9 (a)", rec$recommendation[grepl("split by season", rec$item)], fixed = TRUE)) })
    P <- read.csv(pf, stringsAsFactors = FALSE)
    prow <- function(pair, fit, stream) P[P$pair == pair & P$fit == fit & P$stream == stream, , drop = FALSE]
    j21 <- prow("M2 vs M1", "boat_all_gear", "joint"); j62 <- prow("M6 vs M2", "boat_all_gear", "joint"); j62pc <- prow("M6 vs M2", "boat_pot_closure", "joint")
    chk("1z: the exact joint boat all-gear row is +26.5 at 2.86 SE for the constant term (PASS, primary) and +1.2 at 0.58 SE for the split against it (REVIEW); the pot-closure M6 is M2 (INFO, identical)",
        nrow(j21) == 1 && abs(j21$diff - 26.53) < 0.05 && abs(j21$ratio - 2.86) < 0.02 && j21$verdict == "PASS" && isTRUE(as.logical(j21$primary)) &&
        nrow(j62) == 1 && abs(j62$diff - 1.20) < 0.05 && abs(j62$ratio - 0.58) < 0.02 && j62$verdict == "REVIEW" &&
        nrow(j62pc) == 1 && isTRUE(as.logical(j62pc$identical)) && j62pc$verdict == "INFO")
    BV <- read.csv(bv, stringsAsFactors = FALSE)
    r1m6 <- BV[BV$stage == "M6" & grepl("^R1: boat_all_gear", BV$criterion), ]
    chk("1z: M6's OSP stream was checked EXACTLY against its own loo_pointwise_osp file (the first rendered fit with one), and every M6 R1 row is PASS",
        nrow(r1m6) == 1 && grepl("osp: 130 obs, max |lpd diff|", r1m6$observed, fixed = TRUE) && grepl("MATCHES the committed pointwise file", r1m6$observed, fixed = TRUE) &&
        all(BV$verdict[BV$stage == "M6" & grepl("^R1:", BV$criterion)] == "PASS") && file.exists(file.path(M6, "loo_pointwise_osp_private_boat_all_gear_Dungeness_Kept.csv")))
    chk("1z: no block-CV row is FAIL after the identical-fit floor; the same-fit pairs read INFO",
        !any(BV$verdict == "FAIL") && sum(BV$verdict == "INFO") >= 9)
  } else cat("NOTE  1z: the M6 folder or the result CSVs are absent (sparse checkout); the run checks are skipped\n")

  # ---- (2) the adoption in code -------------------------------------------------------------
  em <- new.env(); em$`%||%` <- function(a, b) if (is.null(a)) b else a
  sys.source("03_R_functions/bss_marine_hazard_covariates.R", envir = em)
  P0 <- list(marine_hazard_selected = list(shore = character(0), private_boat = "nws_sca_any"))
  chk("B42: marine_hazard_terms_for() gives the boat its term in the all-gear fit and withholds it from the pot-closure fit under the shipped regimes, with a note",
      identical(em$marine_hazard_terms_for(c(P0, list(marine_hazard_gear_regimes = "all_gear")), "private_boat", "all_gear"), "nws_sca_any") &&
      { w <- em$marine_hazard_terms_for(c(P0, list(marine_hazard_gear_regimes = "all_gear")), "private_boat", "pot_closure")
        length(w) == 0 && grepl("NOT applied to this fit", attr(w, "note"), fixed = TRUE) && grepl("pot_closure", attr(w, "note"), fixed = TRUE) } &&
      length(em$marine_hazard_terms_for(c(P0, list(marine_hazard_gear_regimes = "all_gear")), "shore", "all_gear")) == 0)
  chk("B42: both regimes named apply the term to every fit (a caller naming no fit is served only then); no key means both; an empty selection stays empty",
      identical(em$marine_hazard_terms_for(c(P0, list(marine_hazard_gear_regimes = c("pot_closure", "all_gear"))), "private_boat", "pot_closure"), "nws_sca_any") &&
      identical(em$marine_hazard_terms_for(c(P0, list(marine_hazard_gear_regimes = c("pot_closure", "all_gear"))), "private_boat", NULL), "nws_sca_any") &&
      identical(em$marine_hazard_terms_for(P0, "private_boat", "pot_closure"), "nws_sca_any") &&   # no key: the pre-2026-09-27 behaviour, both fits
      identical(em$marine_hazard_terms_for(P0, "private_boat", NULL), "nws_sca_any") &&
      length(em$marine_hazard_terms_for(list(), "private_boat", "all_gear")) == 0 &&
      length(em$marine_hazard_terms_for(c(P0, list(marine_hazard_gear_regimes = "all_gear")), "shore", NULL)) == 0)   # nothing selected: nothing to refuse
  chk("B42: a caller that names no fit under a RESTRICTED set is refused, not served the term everywhere (the gate fails closed)",
      { e <- try(em$marine_hazard_terms_for(c(P0, list(marine_hazard_gear_regimes = "all_gear")), "private_boat", NULL), silent = TRUE)
        inherits(e, "try-error") && grepl("gear_regime must be given", conditionMessage(attr(e, "condition")), fixed = TRUE) &&
          grepl("nws_sca_any", conditionMessage(attr(e, "condition")), fixed = TRUE) })
  chk("B42: both drivers hand the preps the sub-season's regime (the only thing the gate depends on)",
      { d1 <- readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE); d2 <- readLines("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd", warn = FALSE)
        d1 <- d1[!grepl("^\\s*#", d1)]; d2 <- d2[!grepl("^\\s*#", d2)]
        i1 <- grep("bss_data_try <- prep_bss_crab_pooled(", d1, fixed = TRUE); i2 <- grep("bss_data_try <- prep_bss_crab_gear(", d2, fixed = TRUE)
        length(i1) == 1 && any(grepl("gear_regime = ss$gear_regime", d1[i1 + 0:3], fixed = TRUE)) &&
          length(i2) == 1 && any(grepl("gear_regime = ss$gear_regime", d2[i2 + 0:3], fixed = TRUE)) })
  chk("B42: an unknown regime name is refused before the archive is read",
      { e <- try(em$marine_hazard_prepare(NULL, list(marine_hazard_mode = "manual", marine_hazard_gear_regimes = "summer", marine_hazard_file = "no-such-workbook.xlsx",
                                                     est_date_start = "2024-09-16", est_date_end = "2025-09-15"), quiet = TRUE), silent = TRUE)
        inherits(e, "try-error") && grepl("marine_hazard_gear_regimes", conditionMessage(attr(e, "condition")), fixed = TRUE) })
  pp <- rd("03_R_functions/prep_bss_crab_pooled.R"); pg <- rd("03_R_functions/prep_bss_crab_gear.R")
  chk("B42: both preps take the fit's marine terms from marine_hazard_terms_for(params, population_name, gear_regime) and print its note",
      grepl("marine_extra <- marine_hazard_terms_for(params, population_name, gear_regime)", pp, fixed = TRUE) &&
      grepl("marine_extra <- marine_hazard_terms_for(params, population_name, gear_regime)", pg, fixed = TRUE) &&
      !grepl("marine_extra <- (params$marine_hazard_selected", pp, fixed = TRUE) && !grepl("marine_extra <- (params$marine_hazard_selected", pg, fixed = TRUE))
  ri <- rd("06_diagnostics/run_improvements_2026-09-08.R"); rg <- rd("06_diagnostics/run_gear_ar_zi_2026-09-13.R"); rm_ <- rd("06_diagnostics/run_marine_hazard_batch_2026-09-25.R")
  chk("A30 ladders: the two pre-covariate ladders set marine_hazard_mode off inside resolve_cfg(), and the marine ladder sets both regimes, none of them in WINDOW",
      grepl('cfg$marine_hazard_mode <- "off"; cfg$marine_hazard_manual_shore <- character(0); cfg$marine_hazard_manual_boat <- character(0)', ri, fixed = TRUE) &&
      grepl('cfg$marine_hazard_mode <- "off"; cfg$marine_hazard_manual_shore <- character(0); cfg$marine_hazard_manual_boat <- character(0)', rg, fixed = TRUE) &&
      grepl('cfg$marine_hazard_gear_regimes <- c("pot_closure", "all_gear")', rm_, fixed = TRUE) &&
      # a WINDOW entry would read `  marine_hazard_gear_regimes = ...` at the start of a code line; the pins read `<-`
      { code <- readLines("06_diagnostics/run_marine_hazard_batch_2026-09-25.R", warn = FALSE); code <- code[!grepl("^\\s*#", code)]
        !any(grepl("^\\s*marine_hazard_gear_regimes\\s*=[^=]", code)) })
  # M6's stamp is 4e23b15's fingerprint (recomputed on that tree, matched); B42 is the only fitting-layer change since
  eM <- new.env(); eM$.here <- function(...) file.path(getwd(), ...); eM$POOLED_RMD <- file.path(getwd(), "01_BSS_models", "BSS-GH-pooled-CPUE-model.Rmd")
  rsM <- readLines("06_diagnostics/run_marine_hazard_batch_2026-09-25.R", warn = FALSE)
  liftM <- function(name) { i <- grep(sprintf("^%s <- function", gsub(".", "\\.", name, fixed = TRUE)), rsM); j <- i; while (!grepl("^\\}", rsM[j])) j <- j + 1L; eval(parse(text = rsM[i:j]), envir = eM) }
  liftM("digest_or_hash"); liftM(".code_fingerprint")
  i <- grep("^CODE_EQUIVALENT_MH <- list\\(", rsM); j <- i; while (!grepl("^\\)", rsM[j])) j <- j + 1L; eval(parse(text = rsM[i:j]), envir = eM)
  i <- grep("^CODE_NOT_EQUIVALENT_MH <- list\\(", rsM)
  if (length(i)) { j <- i; while (!grepl("^\\)", rsM[j])) j <- j + 1L; eval(parse(text = rsM[i:j]), envir = eM) }
  curM <- eM$.code_fingerprint(); m6st <- "05_output/20260927/pooled-CPUE-MH-M6-split/MH_STAGE.txt"
  chk("A30 ladders: M6's recorded fingerprint (4e23b15's) is declared against the ACTUAL current one: equivalent naming B42, or examined and NOT equivalent naming B42 and B44",
      { .hit <- function(lst) { k <- names(lst %||% list()); k[vapply(k, function(x) identical(strsplit(x, " => ", fixed = TRUE)[[1]], c("stan:65b5adeb drivers:a65be4bc fns:6f65d84a", curM)), logical(1))] }
        h1 <- .hit(eM$CODE_EQUIVALENT_MH); h2 <- .hit(eM$CODE_NOT_EQUIVALENT_MH)
        ((length(h1) == 1 && grepl("B42", eM$CODE_EQUIVALENT_MH[[h1]], fixed = TRUE)) ||
         (length(h2) == 1 && grepl("B42", eM$CODE_NOT_EQUIVALENT_MH[[h2]], fixed = TRUE) &&
          grepl("B44", eM$CODE_NOT_EQUIVALENT_MH[[h2]], fixed = TRUE))) &&
          (!file.exists(m6st) || any(grepl("^code: stan:65b5adeb drivers:a65be4bc fns:6f65d84a$", readLines(m6st, warn = FALSE)))) },
      sprintf("(current is %s)", curM))
  chk("A30 ladders: the marine ladder's M0 wiring row looks for the B42 call site, and counts marine_hazard_terms_for() among the module's entry points",
      any(grepl('wired <- grepl("marine_extra <- marine_hazard_terms_for(params, population_name, gear_regime)", src, fixed = TRUE)', rsM, fixed = TRUE)) &&
      any(grepl('"marine_hazard_terms_for")', rsM, fixed = TRUE)))
  chk("A30 ladders: the marine ladder's M0 row accepts the adopted method as well as the pre-adoption off, and both marine runners ship DRY_RUN TRUE after the second slip",
      grepl('identical(BASE$marine_hazard_mode, "manual") && identical(BASE$marine_hazard_manual_boat, "nws_sca_any")', rm_, fixed = TRUE) &&
      any(grepl("^DRY_RUN <- TRUE", readLines("06_diagnostics/run_marine_hazard_batch_2026-09-25.R", warn = FALSE))) &&
      any(grepl("^DRY_RUN <- TRUE", readLines("06_diagnostics/run_marine_block_cv_2026-09-26.R", warn = FALSE))))
  # the improvement ladder's digest hashes run_tag, a LABEL the stage overwrites before any render; run_config.R's
  # label moved on 2026-09-12 and took all five digests with it. resolve_cfg() now holds the label the stamps were
  # computed under: the five recorded digests must reproduce from the runner's own code on the current tree
  eI <- new.env(); eI$.here <- function(...) file.path(getwd(), ...)
  rsI <- readLines("06_diagnostics/run_improvements_2026-09-08.R", warn = FALSE)
  cutI <- grep("^banner\\(sprintf\\(\"IMPROVEMENT LADDER", rsI)[1]
  okI <- tryCatch({ suppressMessages(suppressWarnings(eval(parse(text = rsI[seq_len(cutI - 1L)]), envir = eI))); TRUE }, error = function(e) FALSE)
  recI <- c(R1 = "R1|new_throughout|04971140", R2 = "R2|new_throughout|04aa1a58", R2f = "R2f|new_throughout|04a56a99",
            R4 = "R4|new_throughout|04cc4a57", R5 = "R5|new_throughout|04d9527d")
  gotI <- if (okI) vapply(names(recI), function(sid) tryCatch(eI$stage_digest(sid), error = function(e) NA_character_), character(1)) else rep(NA_character_, 5)
  chk("A30 ladders: the improvement ladder's five recorded digests reproduce on the current tree (run_tag held in resolve_cfg(), never in WINDOW)",
      okI && identical(unname(gotI), unname(recI)) &&
        any(grepl("^  cfg\\$run_tag <- \"two-season-2023-25\"", rsI)) && any(grepl("cfg$run_tag <- st$tag", rsI, fixed = TRUE)) &&
        !any(grepl("^\\s*run_tag\\s*=\\s*\"", rsI)),   # a WINDOW entry would be a quoted literal; the manifest's column is not
      sprintf("(got %s)", paste(gotI, collapse = ", ")))
  chk("A30 ladders: the folders that are on disk carry those digests in IMP_STAGE.txt",
      { st <- Sys.glob(file.path("05_output", "*", "pooled-CPUE-IMP-R*-newf", "IMP_STAGE.txt")); st <- c(st, Sys.glob(file.path("05_output", "*", "gear-type-CPUE-model-IMP-R5-*", "IMP_STAGE.txt")))
        length(st) == 0 || all(vapply(st, function(f) { d <- sub("^digest: ", "", grep("^digest: ", readLines(f, warn = FALSE), value = TRUE)[1]); d %in% recI }, logical(1))) })
  rbc <- rd("06_diagnostics/run_marine_block_cv_2026-09-26.R")
  chk("A30 ladders: the block-CV runner's R1 threshold text says the OSP stream is matched exactly where the rung wrote a pointwise file (both R1 rows use the one string)",
      grepl('R1_THRESHOLD <- "every stream with a committed pointwise file matches it (gear, trailer; OSP where the rung wrote one); otherwise OSP means within 10% of ppc_byobs"', rbc, fixed = TRUE) &&
      lengths(regmatches(rbc, gregexpr("R1_THRESHOLD", rbc, fixed = TRUE))) == 3 && !grepl('"every gear / trailer stream matches its committed file; OSP means within 10% of ppc_byobs"', rbc, fixed = TRUE))
  gi <- readLines(".gitignore", warn = FALSE)
  chk("A30 hygiene: the delivery folder 7e83fab swept in is gone and any delivery folder or README is ignored whatever its name",
      !dir.exists("blockcv-review-patch-2026-09-27") && all(c("*-patch-20*/", "*-patches-20*/", "README-APPLY*.md") %in% gi))
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("A30 config: the method keys sit in section 2.10 and the experiment surface in 4.4b; the archive definition is the method's",
      { src <- readLines("run_config.R", warn = FALSE)
        l210 <- grep("^  # --- 2.10 The effort covariate", src)[1]; l44b <- grep("^  # --- 4.4b Marine hazard effort covariates", src)[1]; l3 <- grep("^  # 3. DIAGNOSTICS REPORTED BESIDE THE ESTIMATE", src)[1]
        keyline <- function(k) grep(sprintf("^\\s*%s\\s*=", k), src)[1]
        is.finite(l210) && is.finite(l44b) && l210 < l3 && l3 < l44b &&
        all(vapply(c("marine_hazard_mode", "marine_hazard_manual_boat", "marine_hazard_manual_shore", "marine_hazard_gear_regimes", "marine_hazard_codes", "marine_hazard_window"),
                   function(k) { l <- keyline(k); is.finite(l) && l > l210 && l < l3 }, logical(1))) &&
        all(vapply(c("marine_hazard_candidates_shore", "marine_hazard_candidates_boat", "marine_hazard_auto_p", "marine_hazard_winter_months", "bar_restriction_impute"),
                   function(k) { l <- keyline(k); is.finite(l) && l > l44b }, logical(1))) })

  # ---- (3) the documents ----------------------------------------------------------------------
  md <- rd("07_documentation/BSS-GH-pooled-CPUE-model-documentation.md"); ps <- rd("07_documentation/development_notes/PIPELINE_STATUS.md")
  cr <- rd("07_documentation/development_notes/CHANGE_REGISTER.md"); vc <- rd("07_documentation/development_notes/VALIDATION_CAMPAIGN.md")
  ng <- rd("07_documentation/NEW_SEASON_GUIDE.md"); cl <- rd("07_documentation/CLAUDE.md")
  chk("A30 docs: the method document moved with the adoption (header, the data streams, 14.2, 21a, limitations 12 and 13)",
      # 2026-09-29: the header moved again (B44 to B51); it must still name A30 among what it carries
      # 2026-10-04: and again (A33)
      grepl("**last moved 2026-10-04: A33", md, fixed = TRUE) &&
      grepl("the NWS Small-Craft-Advisory flag on the private-boat all-gear effort process (A30, 2026-09-27)", md, fixed = TRUE) &&
      grepl("**NWS marine hazard archive** (covariate, not a stream)", md, fixed = TRUE) &&
      grepl("In production\nthere is exactly one, in exactly one fit**", md, fixed = TRUE) &&
      grepl("### 21a. The marine hazard covariate: built 2026-09-25, adopted 2026-09-27 for the boat all-gear fit", md, fixed = TRUE) &&
      grepl("**12. The boat advisory term describes all private boats.**", md, fixed = TRUE) && grepl("**13. The advisory term is one coefficient for the season", md, fixed = TRUE) &&
      !grepl("Production ships `marine_hazard_mode = \"off\"`", md, fixed = TRUE))
  # 2026-09-28: the render the box predicted exists, so the box is the render (section 80 pins its numbers)
  chk("A30 docs: the box is the render the adoption predicted (no longer R4 with a pending render), names the adopted boat fit, and the stale span paragraph is gone",
      # 2026-10-04 (B65): the A31 render is now "THE RUN IT REPLACED", below the A33/A34 box
      grepl("THE METHOD OF RECORD, RENDERED AT THE A31 CODE", ps, fixed = TRUE) && grepl("**private boat all-gear 47,105**", ps, fixed = TRUE) &&
      !grepl("IT PREDATES THE METHOD OF RECORD BY ONE TERM", ps, fixed = TRUE) && !grepl("render pending", ps, fixed = TRUE) &&
      !grepl("WHAT R4 DOES NOT MATCH IS THE SHIPPED WINDOW", ps, fixed = TRUE) &&
      grepl("[FOUND AND FIXED 2026-09-27] The improvement ladder's five stage digests had drifted from their folders, by a label.", ps, fixed = TRUE) &&
      grepl("hashes `run_tag`", ps, fixed = TRUE) && !grepl("[OPEN 2026-09-27] The improvement ladder's R4 stage digest", ps, fixed = TRUE))
  chk("A30 docs: the register says ADOPTED with the scope, carries B42, records D28 resolved and D32 as a limitation of the adopted method",
      grepl("| **ADOPTED 2026-09-27 (Matt's decision, Section 1z): the constant boat SCA term, the ALL-GEAR fit only; no shore term, not the bar tick, not the pot-closure fit.**", cr, fixed = TRUE) &&
      grepl("| B42 |", cr, fixed = TRUE) && grepl("**RESOLVED 2026-09-12 by option (a), recorded here 2026-09-27**", cr, fixed = TRUE) &&
      grepl("now a limitation of the ADOPTED method", cr, fixed = TRUE))
  chk("A30 docs: the campaign has Section 1z with its anchor, the guide makes the archive a per-season requirement, and CLAUDE.md names the method of record",
      grepl("## 1z. Rung M6 and the four-season screen read; the decision", vc, fixed = TRUE) && grepl("> | 1z | 2026-09-27 |", vc, fixed = TRUE) && grepl("7e83fab", vc, fixed = TRUE) &&
      grepl("**Required since 2026-09-27**", ng, fixed = TRUE) && grepl("The NWS archive covers the window", ng, fixed = TRUE) &&
      grepl("the method of record since 2026-09-27", cl, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 80. THE METHOD OF RECORD, RENDERED (2026-09-28; VALIDATION_CAMPAIGN Section 1z.5,
#     CHANGE_REGISTER A30 RENDERED, B43, six C rows).
#
#     Matt rendered run_config.R as shipped (source("run_estimation.R"), code 74d9731) and
#     committed the folder as 1d3409d. (1) The render as evidence: the numbers the box quotes,
#     read back from the folder's own files, and the claim that makes the render a
#     confirmation rather than a new result: every fit's full posterior summary is
#     byte-identical to the fit the adoption named. (2) B43, the four defects reading it
#     found: the season totals' census, the report's missing fitted term, the last unseeded
#     diagnostic, the manifest's truncated configuration; tested by running them where the
#     harness can (the totals chunks, the table, the manifest) and by their source where it
#     cannot (the overdispersion subsample and the report chunks need a stanfit or a render). (3) The documents, including a check that every table row in the governed
#     documents has its header's cell count (eight register rows did not).
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  P  <- "05_output/20260927/pooled-CPUE-canonical-2024-25"
  R4 <- "05_output/20260910/pooled-CPUE-IMP-R4-shore-tau-newf"
  M2 <- "05_output/20260925/pooled-CPUE-MH-M2-sca"
  M6 <- "05_output/20260927/pooled-CPUE-MH-M6-split"
  FIT <- c(spc = "shore_ring_net_only_Dungeness_Kept", sag = "shore_all_gear_Dungeness_Kept",
           bpc = "private_boat_ring_net_only_Dungeness_Kept", bag = "private_boat_all_gear_Dungeness_Kept")
  same <- function(a, b, f) file.exists(file.path(a, f)) && file.exists(file.path(b, f)) &&
    identical(unname(tools::md5sum(file.path(a, f))), unname(tools::md5sum(file.path(b, f))))

  # ---- (1) the render as evidence ---------------------------------------------------------------
  if (!(dir.exists(P) && file.exists(file.path(P, "port_total_Dungeness_Kept.csv")))) skp("80 render: the 2026-09-27 render", "folder absent")
  if (dir.exists(P) && file.exists(file.path(P, "port_total_Dungeness_Kept.csv"))) {
    rp <- readLines(file.path(P, "run_parameters.txt"), warn = FALSE)
    chk("80 render: run_parameters.txt records the method of record's keys (manual, boat nws_sca_any, all_gear, no shore term)",
        any(grepl('marine_hazard_mode\\s*: chr "manual"', rp)) && any(grepl('marine_hazard_manual_boat\\s*: chr "nws_sca_any"', rp)) &&
        any(grepl('marine_hazard_gear_regimes\\s*: chr "all_gear"', rp)) && any(grepl("marine_hazard_manual_shore\\s*: chr\\(0\\)", rp)))
    cv <- read.csv(file.path(P, "convergence_report.csv"), stringsAsFactors = FALSE)
    chk("80 render: 4 of 4 fits report BSS and pass the gate; worst divergence fraction 1.78%, R-hat within 1.0007, impact at most 0.0016 SD",
        nrow(cv) == 4 && all(cv$method_selected == "BSS") && all(cv$pass_convergence) &&
        abs(max(cv$divergence_fraction) - 0.0178) < 1e-9 && max(cv$C_sum_rhat, cv$E_sum_rhat, cv$B1_C_rhat) <= 1.0007 &&
        max(cv$impact_C_sd, cv$impact_E_sd) <= 0.0016 && min(cv$C_sum_neff, cv$E_sum_neff) == 4604 && max(cv$C_sum_neff, cv$E_sum_neff) == 12334)
    ae <- read.csv(file.path(P, "ar_escalation_log.csv"), stringsAsFactors = FALSE)
    comp <- setNames(round(ae$catch_median), ae$fit)
    ps_ <- read.csv(file.path(P, "pe_port_summary.csv"), stringsAsFactors = FALSE)
    chk("80 render: the components are the ones the adoption predicted (8,963 / 29,210 / 1,372 / 47,319; census 8,538; medians sum to 95,402)",
        identical(unname(comp[FIT]), c(8963, 29210, 1372, 47319)) &&
        round(ps_$Dungeness[ps_$Component == "Commercial/Charter"]) == 8538 &&
        round(sum(ae$catch_median) + ps_$Dungeness[ps_$Component == "Commercial/Charter"]) == 95402)
    pt <- read.csv(file.path(P, "port_total_Dungeness_Kept.csv"), stringsAsFactors = FALSE)
    row_ <- function(k) unlist(pt[pt$Estimate == k, c("BSS_median", "BSS_lo95", "BSS_hi95")])
    chk("80 render: port 96,118 [79,418, 120,558]; predictive 96,119 [79,418, 120,553]; effort 58,963 [50,888, 68,244]; PE 85,076",
        identical(unname(row_("Expected_Catch")), c(96118L, 79418L, 120558L)) &&
        identical(unname(row_("Predictive_Catch")), c(96119L, 79418L, 120553L)) &&
        identical(unname(row_("Effort")), c(58963L, 50888L, 68244L)) && pt$PE[pt$Estimate == "Expected_Catch"] == 85076)
    if (all(dir.exists(c(R4, M2)))) {
      fs <- function(k) sprintf("bss_full_summary_%s.csv", FIT[[k]])
      chk("80 render: every fit's full posterior summary is byte-identical to the fit the adoption named (shore and boat pot closure = R4, boat all-gear = M2)",
          same(P, R4, fs("spc")) && same(P, R4, fs("sag")) && same(P, R4, fs("bpc")) && same(P, M2, fs("bag")))
      chk("80 render: and nothing else: the boat all-gear fit is not R4's (the term moved it) and the pot-closure boat fit is not M2's (the term was withheld, B42)",
          !same(P, R4, fs("bag")) && !same(P, M2, fs("bpc")))
    } else cat("NOTE  80: R4 or M2 folder absent; the byte-identity checks are skipped\n")
    sp <- read.csv(file.path(P, sprintf("structural_params_%s.csv", FIT[["bag"]])), stringsAsFactors = FALSE)
    bo <- sp[sp$parameter == "B_open[1]", ]
    oc <- list.files(P, pattern = "^opener_covariates_.*\\.csv$")
    chk("80 render: the adopted coefficient is -1.16 [-1.45, -0.87], and only the boat all-gear fit carries a K_open column (nws_sca_any)",
        nrow(bo) == 1 && abs(bo$median + 1.1625) < 5e-4 && abs(bo$lo95 + 1.4521) < 5e-4 && abs(bo$hi95 + 0.8678) < 5e-4 &&
        identical(oc, sprintf("opener_covariates_%s.csv", FIT[["bag"]])) &&
        identical(read.csv(file.path(P, oc[1]), stringsAsFactors = FALSE)$opener, "nws_sca_any"))
    st <- read.csv(file.path(P, "season_totals.csv"), stringsAsFactors = FALSE)
    chk("80 render: its season_totals.csv reads 96,110 [79,428, 120,560], the census-as-a-constant defect B43 fixes (kept as rendered)",
        nrow(st) == 1 && st$BSS_median == 96110 && st$BSS_lo95 == 79428 && st$BSS_hi95 == 120560)
    mf <- "05_output/20260927/run_manifest_20260927_192233.txt"
    if (file.exists(mf)) chk("80 render: its manifest names code 74d9731 and was cut off at 99 keys by str() (the defect B43 fixes)",
        { m <- readLines(mf, warn = FALSE); any(grepl("^git sha\\s*: 74d9731$", m)) && any(grepl("list output truncated", m, fixed = TRUE)) &&
          !any(grepl("marine_hazard_mode", m, fixed = TRUE)) })
    if (!dir.exists(M6)) skp("80 render: B38 in the field (M6)", "folder absent")
    if (dir.exists(M6)) {
      chk("80 render: B38 in the field: M6 and this render (both post-B38) wrote byte-identical saved draws for the shore fits; R4 (pre-B38) did not",
          same(P, M6, sprintf("bss_draws_summed_%s.csv", FIT[["sag"]])) && same(P, M6, sprintf("bss_draws_summed_%s.csv", FIT[["spc"]])) &&
          (!dir.exists(R4) || !same(P, R4, sprintf("bss_draws_summed_%s.csv", FIT[["sag"]]))))
      chk("80 render: the overdispersion diagnostic was the unseeded exception (M6 and this render differ there, on identical shore fits)",
          !same(P, M6, sprintf("effort_overdispersion_decomp_%s.csv", FIT[["sag"]])) && same(P, M6, sprintf("ppc_calibration_%s.csv", FIT[["sag"]])))
    }
  } else cat("NOTE  80: the authoritative run's folder is absent; the render assertions are skipped\n")

  # ---- (2) B43: the season totals carry the census draws (run both chunks on synthetic draws) --
  rmd <- readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE)
  chunk <- function(name) { i <- grep(sprintf("^```\\{r %s[,}]", name), rmd); j <- i + 1L
                            while (!grepl("^```\\s*$", rmd[j])) j <- j + 1L; rmd[(i + 1L):(j - 1L)] }
  # mode: "one" (a single season, census inside it), "one_late" (a single season whose window starts inside the
  # census window, so census_start_date is in no sub-season), "two" (two seasons, per-season census windows),
  # "one_plus_foreign" (a single season with a census_windows entry for a season outside the run)
  run_totals <- function(mode, seed_before, big_se = FALSE) {
    e <- new.env(parent = globalenv())
    sys.source("03_R_functions/bss_convergence_gate.R", envir = e)
    e$`%||%` <- `%||%`; e$timer_start <- function(...) invisible(NULL); e$report_table <- function(...) invisible(NULL)
    e$tibble <- tibble::tibble; e$bss_catch_groups <- "Dungeness_Kept"; e$output_dir <- tempfile(); dir.create(e$output_dir)
    e$params <- list(bss_seed = 20260619, census_start_date = if (mode == "one_late") "2024-11-15" else "2024-12-01",
                     census_end_date = "2025-02-08", bss_chains = 4, bss_iter_default = 2000, bss_warmup_default = 1000)
    ss_ <- function(nm, a, b, sn) list(name = nm, start = as.Date(a), end = as.Date(b), season = sn)
    e$subseasons <- if (mode == "one_late") list(ss_("all_gear", "2024-12-01", "2025-09-15", "2024-25")) else
      list(ss_("ring_net_only", "2024-09-16", "2024-11-30", "2024-25"), ss_("all_gear", "2024-12-01", "2025-09-15", "2024-25"))
    if (mode == "two") e$subseasons <- c(e$subseasons, list(ss_("ring_net_only_2025-26", "2025-09-16", "2025-11-30", "2025-26"),
                                                             ss_("all_gear_2025-26", "2025-12-01", "2026-09-15", "2025-26")))
    set.seed(1); e$bss_all <- list(); e$pe_all <- list()
    for (pop in c("shore", "private_boat")) for (s0 in e$subseasons) {
      m <- if (pop == "shore") 20000 else 30000
      e$bss_all[[paste0(pop, "_", s0$name, "_Dungeness_Kept")]] <- list(C_exp_draws = rlnorm(4000, log(m), 0.15), C_draws = rlnorm(4000, log(m), 0.16),
                                                                         E_draws = rlnorm(4000, log(m / 2), 0.1), pe_fallback = FALSE, use_bss = TRUE)
      e$pe_all[[paste0(pop, "_", s0$name)]] <- list(Dungeness_Kept = m * 0.9, effort_total = m / 2)
    }
    cc1 <- list(Dungeness_Kept = 8538, carried_se = 73, observed_dung = 7691, effort_total = 198, census_uncertainty = "charter",
                frame_warnings = character(0), daily_full = data.frame(observed = c(TRUE, FALSE)), Dungeness_Kept_se = 140,
                charter_se = 73, commercial_se = 123, imputation_var = 0, census_expansion = "none")
    cc2 <- modifyList(cc1, list(Dungeness_Kept = 5000, carried_se = if (big_se) 3000 else 300, observed_dung = 0))
    e$pe_all$comm_charter <- switch(mode,
      one = cc1, one_late = cc1,
      two = , one_plus_foreign = modifyList(cc1, list(Dungeness_Kept = 13538, carried_se = sqrt(73^2 + cc2$carried_se^2), observed_dung = 7691,
                                                     by_season = list(`2024-25` = cc1, `2025-26` = cc2))))
    e$pe_port_total <- function(field) sum(vapply(e$pe_all[names(e$pe_all) != "comm_charter"], function(x) x[[field]] %||% 0, numeric(1))) +
      (e$pe_all$comm_charter[[field]] %||% 0)
    utils::capture.output(eval(parse(text = chunk("port-total")), envir = e))
    set.seed(seed_before); before <- get(".Random.seed", envir = globalenv())
    utils::capture.output(eval(parse(text = chunk("season-totals")), envir = e))
    # the season's own draws with the census as a CONSTANT, the pre-2026-09-28 construction, to compare widths against
    const <- lapply(unique(vapply(e$subseasons, function(x) x$season, character(1))), function(sn) {
      z <- rep(0, 4000); for (pop in c("shore", "private_boat")) for (s0 in e$subseasons) if (identical(s0$season, sn))
        z <- z + e$bss_all[[paste0(pop, "_", s0$name, "_Dungeness_Kept")]]$C_exp_draws
      z + (if (is.null(e$pe_all$comm_charter$by_season)) e$pe_all$comm_charter$Dungeness_Kept else e$pe_all$comm_charter$by_season[[sn]]$Dungeness_Kept) })
    list(pt = read.csv(file.path(e$output_dir, "port_total_Dungeness_Kept.csv")), st = read.csv(file.path(e$output_dir, "season_totals.csv")),
         rng_restored = identical(before, get(".Random.seed", envir = globalenv())),
         const_width = vapply(const, function(z) diff(stats::quantile(z, c(0.025, 0.975), names = FALSE)), numeric(1)))
  }
  one <- tryCatch(run_totals("one", 99), error = function(e) e)
  same_as_port <- function(r) !inherits(r, "error") && { x <- r$pt[r$pt$Estimate == "Expected_Catch", ]
    nrow(r$st) == 1 && identical(c(x$BSS_median, x$BSS_lo95, x$BSS_hi95), c(r$st$BSS_median, r$st$BSS_lo95, r$st$BSS_hi95)) }
  chk("B43: on one season, season_totals.csv IS the port total (median and interval, draw for draw), and the season block restores the RNG",
      same_as_port(one) && one$rng_restored, if (inherits(one, "error")) sprintf("[%s]", conditionMessage(one)) else "")
  late <- tryCatch(run_totals("one_late", 99), error = function(e) e)
  chk("B43: a single season whose window starts inside the census window still carries the census (the season total IS the port total)",
      same_as_port(late) && late$st$Census == 8538, if (inherits(late, "error")) sprintf("[%s]", conditionMessage(late)) else "")
  two_a <- tryCatch(run_totals("two", 99, big_se = TRUE), error = function(e) e); two_b <- tryCatch(run_totals("two", 5, big_se = TRUE), error = function(e) e)
  chk("B43: on two seasons each season draws its OWN census: the interval is wider than a constant census gives, the totals do not depend on the RNG state before, and the RNG is restored",
      !inherits(two_a, "error") && !inherits(two_b, "error") && nrow(two_a$st) == 2 && identical(two_a$st, two_b$st) &&
      two_a$rng_restored && two_b$rng_restored && two_a$st$Census[2] == 5000 &&
      (two_a$st$BSS_hi95[2] - two_a$st$BSS_lo95[2]) > 1.05 * two_a$const_width[2])
  foreign <- tryCatch(run_totals("one_plus_foreign", 99), error = function(e) e)
  chk("B43: a census window for a season outside the run stays out of this season's total (the port draws, which carry it, are not reused)",
      !inherits(foreign, "error") && nrow(foreign$st) == 1 && foreign$st$Census == 8538 &&
      foreign$st$BSS_median < foreign$pt$BSS_median[foreign$pt$Estimate == "Expected_Catch"] - 3000)
  chk("B43: the port block keeps its census draws for the season block, which reuses them when one season holds the whole census",
      any(grepl(".cc_draws_by_cg[[cg]] <- .cc_draws", rmd, fixed = TRUE)) &&
      any(grepl("if (.n_census_seasons == 1L && .census_is_whole(sn) && !is.null(.cc_draws_by_cg[[cg]])) .cc_draws_by_cg[[cg]] else", rmd, fixed = TRUE)) &&
      !any(grepl("sC <- sC + (ccs[[cg]] %||% 0); sE <- sE + (ccs$effort_total %||% 0)", rmd, fixed = TRUE)))

  # ---- (2b) B43: the fitted day-covariate table -------------------------------------------------
  ed <- new.env(); sys.source("03_R_functions/bss_day_covariate_report.R", envir = ed)
  X1 <- c(rep(1, 10), rep(0, 30)); stub <- function(fit) fit$draws
  mk <- function(pop, ss, K, lab, gear = FALSE, draws = NULL) {
    bd <- list(K_open = K, D = 40L, X_open_flat = if (K) rep(X1, K) else numeric(0))
    if (gear) bd$.opener_labels <- paste(lab, collapse = ",") else attr(bd, "opener_labels") <- lab
    list(fit = list(draws = draws), bss_data = bd, pe_fallback = FALSE, population = pop, subseason = ss)
  }
  set.seed(3); dr <- matrix(rnorm(4000, -1.16, 0.15), ncol = 1)
  ba <- list(a = mk("shore", "all_gear", 0L, character(0)), b = mk("private_boat", "all_gear", 1L, "nws_sca_any", draws = dr),
             c = mk("private_boat", "ring_net_only", 1L, "halibut_open", gear = TRUE, draws = dr))
  tb <- ed$bss_day_covariate_table(ba, list(marine_hazard_selected = list(private_boat = "nws_sca_any")), extract_b_open = stub)
  chk("B43 table: one row per fit x column, '(none)' for a fit without one, the flagged days counted, marine and opener sources told apart, both prep conventions read",
      identical(tb$covariate, c("(none)", "nws_sca_any", "halibut_open")) && identical(tb$source, c(NA, "marine hazard", "other-fishery opener")) &&
      identical(tb$days_flagged, c(NA, 10L, 10L)) && all(tb$fit_days == 40L) &&
      abs(tb$coef_median[2] - stats::median(dr)) < 1e-3 && identical(tb$rate_ratio[2], sprintf("%.2f (%.2f-%.2f)", exp(stats::median(dr)),
        exp(stats::quantile(dr, 0.025, names = FALSE)), exp(stats::quantile(dr, 0.975, names = FALSE)))))
  chk("B43 table: a fit whose labels do not match its K_open is refused, and a PE-fallback fit is skipped",
      inherits(try(ed$bss_day_covariate_table(list(x = mk("private_boat", "all_gear", 2L, "nws_sca_any", draws = cbind(dr, dr))), list(), extract_b_open = stub), silent = TRUE), "try-error") &&
      nrow(ed$bss_day_covariate_table(list(x = modifyList(ba$b, list(pe_fallback = TRUE))), list(), extract_b_open = stub)) == 0)
  gd <- readLines("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd", warn = FALSE)
  chk("B43 table: in the gear report the table is the chunk's last, top-level, visible value (inside an if-block with a write.csv after it, knitr shows the code and no table)",
      { k <- grep("the table must be the chunk's LAST, top-level, visible value", gd, fixed = TRUE)
        length(k) == 1 && grepl("^if \\(!is.null\\(day_cov_df\\) && nrow\\(day_cov_df\\) > 0\\)$", gd[k + 1]) && grepl("^  report_table\\(", gd[k + 2]) } &&
      { k <- grep("^```\\{r day-covariates-show", rmd); length(k) == 1 && grepl("^if \\(!is.null\\(day_cov_df\\)", rmd[k + 1]) && grepl("^  report_table\\(", rmd[k + 2]) })
  chk("B43 table: both reports carry it (pooled 12.2, gear 13.1) and write effort_day_covariates.csv",
      any(grepl("^### 12.4 Effort day covariates as fitted", rmd)) && any(grepl("^### 13.1 Effort day covariates as fitted", gd)) &&
      any(grepl("day_cov_df <- tryCatch(bss_day_covariate_table(bss_all, params)", rmd, fixed = TRUE)) &&
      any(grepl("day_cov_df <- tryCatch(bss_day_covariate_table(bss_all, params)", gd, fixed = TRUE)) &&
      sum(grepl('"effort_day_covariates.csv"', c(rmd, gd), fixed = TRUE)) == 2)

  # ---- (2c) B43: the overdispersion subsample is seeded, the caller's RNG restored -------------
  od <- readLines("03_R_functions/diagnose_effort_overdispersion.R", warn = FALSE)
  chk("B43: write_effort_overdispersion_diag() takes a seed, sets it before its sample.int(), and restores the caller's RNG on exit",
      any(grepl("n_draws_use = 2000, seed = 1L) {", od, fixed = TRUE)) &&
      { i <- grep("set.seed(as.integer(seed))", od, fixed = TRUE); j <- grep("use <- if (ndraw > n_draws_use) sort(sample.int(ndraw, n_draws_use))", od, fixed = TRUE)
        length(i) == 1 && length(j) == 1 && i < j } &&
      any(grepl('if (.had_seed) assign(".Random.seed", .old_seed, envir = globalenv())', od, fixed = TRUE)))

  # ---- (2d) B43: the manifest records every key and the tree state (written and read back) ------
  rs <- readLines("run_estimation.R", warn = FALSE)
  i <- grep("^write_manifest <- function", rs); j <- i; while (!grepl("^\\}", rs[j])) j <- j + 1L
  em <- new.env(); eval(parse(text = rs[i:j]), envir = em)
  sys.source("run_config.R", envir = em); em$run_stamp <- "HARNESS"; em$model <- "pooled"
  mfp <- tryCatch(em$write_manifest(list(model = list(minutes = 1, outdir = "x")), tempdir()), error = function(e) NA_character_)
  chk("B43: the manifest writes every run_config key (no truncation) and a git-tree line",
      !is.na(mfp) && file.exists(mfp) && { m <- readLines(mfp, warn = FALSE)
        sum(grepl("^ \\$ ", m)) == length(em$run_config) && !any(grepl("list output truncated", m, fixed = TRUE)) &&
        any(grepl("^git tree    : ", m)) && any(grepl("marine_hazard_gear_regimes", m, fixed = TRUE)) })
  # outside a repository git answers nothing with a non-zero status: that must read "unknown", never "clean"
  od <- setwd(tempdir())   # restored on the next line whatever happens; an on.exit() here would fire only when local() ends
  mfo <- tryCatch({ x <- NA_character_; suppressWarnings(utils::capture.output(x <- em$write_manifest(list(model = list(minutes = 1, outdir = "x")), tempdir()), type = "message")); x },
                  error = function(e) NA_character_)
  setwd(od)
  chk("B43: where git does not answer (not a repository), the manifest says the tree is unknown rather than clean",
      !is.na(mfo) && file.exists(mfo) && { m <- readLines(mfo, warn = FALSE)
        any(grepl("^git tree    : unknown", m)) && !any(grepl("^git tree    : clean", m)) })
  chk("B43: both ladders' code-equivalence declarations name the 2026-09-28 post-fit changes",
      grepl("2026-09-28 review of the first render of the method of record (B43)", rd("06_diagnostics/run_improvements_2026-09-08.R"), fixed = TRUE) &&
      grepl("B43 (2026-09-28) is post-fit reporting and diagnostics only", rd("06_diagnostics/run_marine_hazard_batch_2026-09-25.R"), fixed = TRUE))

  # ---- (3) the documents ------------------------------------------------------------------------------
  ps <- rd("07_documentation/development_notes/PIPELINE_STATUS.md"); cr <- rd("07_documentation/development_notes/CHANGE_REGISTER.md")
  vc <- rd("07_documentation/development_notes/VALIDATION_CAMPAIGN.md"); md <- rd("07_documentation/BSS-GH-pooled-CPUE-model-documentation.md")
  flat <- function(x) gsub("[ \n>]+", " ", x)
  # 2026-09-29 (B50): the 2026-09-27 render was superseded by the re-render at B44 to B49 (section 85);
  # the documents must now record it as what the authoritative run REPLACED, with its own total.
  chk("80 docs: the box records the 2026-09-27 render as what the authoritative run replaced, with its total and interval",
      # 2026-10-02 (B61): the box moved again, so the 2026-09-27 render is now two runs back in "What this replaced"
      grepl("`05_output/20260927/pooled-CPUE-canonical-2024-25`, 96,118 [79,418, 120,558]", flat(ps), fixed = TRUE) && grepl("**What this replaced.**", ps, fixed = TRUE) &&
      grepl("Section 1z.5", ps, fixed = TRUE))
  chk("80 docs: the register moves its authoritative run, marks A30 RENDERED, carries B43 and six 2026-09-28 defect rows",
      grepl("`05_output/20260927/pooled-CPUE-canonical-2024-25` (96,118", cr, fixed = TRUE) &&
      # 2026-09-29: A30's run is no longer the authoritative one (B50 superseded it), so the row
      # says so; the pin follows the row rather than holding it to a stale claim.
      grepl("**RENDERED 2026-09-28 (`1d3409d`, Section 1z.5), the first render of the method of record (superseded as the authoritative run", cr, fixed = TRUE) &&
      grepl("| B43 |", cr, fixed = TRUE) &&
      # at least B43's six: B44's five defect rows landed the same day (>= rather than ==, 2026-09-28)
      length(gregexpr("| 2026-09-28 |", cr, fixed = TRUE)[[1]]) >= 6 && grepl("**Confirmed in the field 2026-09-28**", cr, fixed = TRUE))
  chk("80 docs: the campaign has 1z.5 and the 1z anchor names the render's commit and folder; the method document's reference run moved",
      grepl("### 1z.5 The confirming render (2026-09-28)", vc, fixed = TRUE) && grepl("`1d3409d` (the confirming render)", vc, fixed = TRUE) &&
      grepl("It superseded `05_output/20260927/pooled-CPUE-canonical-2024-25` (96,118 [79,418, 120,558]", md, fixed = TRUE) &&
      !grepl("worst divergence fraction 2.17% against the\n5% backstop", md, fixed = TRUE))
  # every table row in the governed documents has its header's cell count (GitHub drops excess cells)
  ncell <- function(line) { x <- sub("^\\s*>?\\s*", "", line); x <- sub("^\\|", "", x); x <- sub("\\s+$", "", x)
                            if (grepl("[^\\\\]\\|$", x)) x <- sub("\\|$", "", x)
                            length(regmatches(x, gregexpr("(?<!\\\\)\\|", x, perl = TRUE))[[1]]) + 1L }
  bad_rows <- character(0)
  for (f in c("07_documentation/development_notes/CHANGE_REGISTER.md", "07_documentation/development_notes/PIPELINE_STATUS.md",
              "07_documentation/development_notes/VALIDATION_CAMPAIGN.md", "07_documentation/BSS-GH-pooled-CPUE-model-documentation.md",
              "07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md", "07_documentation/archive/PR-5-merge-OSP-boat-count-incorporation.md", "07_documentation/CLAUDE.md",
              "07_documentation/NEW_SEASON_GUIDE.md", "README-R-functions.md", "06_diagnostics/README.md", "07_documentation/README.md", "04_input_files/README.md")) {
    L <- readLines(f, warn = FALSE); k <- 1L
    while (k < length(L)) {
      if (grepl("^\\s*>?\\s*\\|", L[k]) && grepl("^\\s*>?\\s*\\|?\\s*:?-{3,}", L[k + 1L])) {
        h <- ncell(L[k]); m <- k + 2L
        while (m <= length(L) && grepl("^\\s*>?\\s*\\|", L[m])) { if (ncell(L[m]) != h) bad_rows <- c(bad_rows, sprintf("%s:%d", basename(f), m)); m <- m + 1L }
        k <- m
      } else k <- k + 1L
    }
  }
  chk("80 docs: every table row in the governed documents has its header's cell count (a `|` in code is escaped, no fifth cell in a four-column table)",
      length(bad_rows) == 0, sprintf("(%s)", paste(utils::head(bad_rows, 10), collapse = ", ")))
})

# ---------------------------------------------------------------------------
# 81. B44 (2026-09-28): FOUR REVIEW FINDINGS, FIXED. The gear-resolved port total summed
#     the PREDICTIVE catch while the pooled one sums the EXPECTED catch; distinct
#     interviews sharing an interview_id were each given the id's summed catch; the PE's
#     monthly split assumed one CPUE for the sub-season and weighted months by sampled days
#     only; and the PE effort SE summed cells that share a donor mean as if independent.
# ---------------------------------------------------------------------------
local({
  # ---- 81a. the gear port total is the pooled quantity, under the pooled labels ----------
  g  <- readLines("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd", warn = FALSE)
  gc <- paste(g[!grepl("^\\s*#", g)], collapse = "\n")
  chk("B44 gear: the fit entry carries the expected-catch draws (C_expected_sum)",
      grepl('C_exp_draws <- rstan::extract(fit, "C_expected_sum")$C_expected_sum', gc, fixed = TRUE) &&
      grepl("C_draws = C_draws, C_exp_draws = C_exp_draws, C_pred_draws = C_pred_draws, E_draws = E_draws", gc, fixed = TRUE))
  chk("B44 gear: the port total sums C_exp_draws (expected) and carries C_draws as Predictive_Catch",
      grepl("bss_C_total <- bss_C_total + b$C_exp_draws[idx]", gc, fixed = TRUE) &&
      grepl("bss_C_pred_total <- bss_C_pred_total + (b$C_pred_draws %||% b$C_draws)[idx]", gc, fixed = TRUE) &&   # B46: the ZINB predictive
      !grepl("bss_C_total <- bss_C_total + b$C_draws[idx]", gc, fixed = TRUE))
  chk("B44 gear: port_total rows are the pooled track's (Effort / Expected_Catch / Predictive_Catch)",
      grepl('Estimate = c("Effort", "Expected_Catch", "Predictive_Catch")', gc, fixed = TRUE) &&
      !grepl('Estimate = c("Effort (gear-deployments)", "Catch")', gc, fixed = TRUE))
  chk("B44 gear: the monthly block sums the expected daily catch (C_total), not C_gear_pred",
      grepl('C_daily_mat <- rstan::extract(b$fit, "C_total")$C_total', gc, fixed = TRUE) &&
      !grepl('rstan::extract(b$fit, "C_gear_pred")', gc, fixed = TRUE))
  chk("B44 gear: the PE-vs-BSS comparison and the tau projection read the expected catch",
      grepl('bss_catch <- b$C_expected_sum["50%"]', gc, fixed = TRUE) &&
      grepl("boat_draws  <- b_boat$C_exp_draws[rep_len(1:nb, n_draws_max)]", gc, fixed = TRUE))
  chk("B44 gear: every runner that reads the row matches Expected_Catch by the shared pattern",
      all(vapply(c("06_diagnostics/run_adoption_2026-09-07.R"), function(f)
        any(grepl('grepl("^(Expected_)?Catch$", x$Estimate)', readLines(f, warn = FALSE), fixed = TRUE)), logical(1))) &&
      grepl("^(Expected_)?Catch$", "Expected_Catch") && !grepl("^(Expected_)?Catch$", "Predictive_Catch"))

  # ---- 81b. repair_interview_ids() -----------------------------------------------------
  source("03_R_functions/fetch_crab_data.R")
  u <- data.frame(interview_id = c("S1_1", "S1_2", "S2_1"), season = "2024-25", dungeness_kept = c(1, 2, 3))
  chk("B44 ids: a frame without shared ids passes through identical (no attribute, no rename)",
      identical(repair_interview_ids(u, quiet = TRUE), u))
  d <- data.frame(interview_id = c("S1_NA", "S1_NA", "S2_7", "S1_NA", "S3_1"),
                  season = c("2022-23", "2022-23", "2024-25", "2022-23", "2024-25"),
                  dungeness_kept = c(1, 0, 4, 8, 2))
  r <- repair_interview_ids(d, quiet = TRUE)
  chk("B44 ids: shared ids become unique per row, in row order, and unique ids are untouched",
      identical(r$interview_id, c("S1_NA__r1", "S1_NA__r2", "S2_7", "S1_NA__r3", "S3_1")) &&
      !anyDuplicated(r$interview_id))
  chk("B44 ids: the repair is logged (id, rows, seasons) and loses no row",
      nrow(r) == nrow(d) && identical(attr(r, "interview_id_repairs")$interview_id, "S1_NA") &&
      attr(r, "interview_id_repairs")$n_rows == 3L && attr(r, "interview_id_repairs")$seasons == "2022-23")
  cw <- r |> dplyr::group_by(interview_id) |> dplyr::summarise(fish = sum(dungeness_kept), .groups = "drop")
  back <- dplyr::left_join(r, cw, by = "interview_id")
  chk("B44 ids: after the repair the catch join returns each row its OWN catch (9 over the id, not 27)",
      identical(back$fish, d$dungeness_kept) && sum(back$fish) == sum(d$dungeness_kept))
  fc <- readLines("03_R_functions/fetch_crab_data.R", warn = FALSE)
  i_rep <- grep("gh_interview <- repair_interview_ids(gh_interview)", fc, fixed = TRUE)
  i_cat <- grep("  catch <- gh_interview |> filter(dungeness_kept>0) |>", fc, fixed = TRUE)
  chk("B44 ids: fetch_crab_data() repairs the ids BEFORE the catch table is built",
      length(i_rep) == 1 && length(i_cat) == 1 && i_rep < i_cat)

  # ---- 81c. pe_monthly_split() ---------------------------------------------------------
  source("03_R_functions/pe_monthly_split.R")
  cal <- data.frame(event_date = as.Date("2025-01-29") + 0:9)
  cal$period <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2)[seq_len(nrow(cal))]   # period 1 straddles Jan/Feb
  cal$day_type <- "weekday"; cal$open_section_1 <- TRUE
  es <- data.frame(section_num = 1, period = c(1, 2), day_type = "weekday",
                   n_total_days = c(4, 6), est_total = c(40, 60))
  cs <- cbind(es, est_catch = c(80, 6))   # cell 1 CPUE 2.0, cell 2 CPUE 0.1
  res <- list(effort_strata = es, catch_strata = list(Dungeness_Kept = cs))
  sp <- pe_monthly_split(res, cal, "Dungeness_Kept")
  chk("B44 split: the months sum to the component totals exactly (effort and catch)",
      isTRUE(all.equal(sum(sp$month_effort), 100)) && isTRUE(all.equal(sum(sp$month_catch), 86)))
  chk("B44 split: a cell is spread over its OWN calendar days (period 1: 3 Jan days, 1 Feb day)",
      isTRUE(all.equal(sp$month_effort, c(30, 70))) && isTRUE(all.equal(sp$month_catch, c(60, 26))))
  chk("B44 split: each month carries its cells' own CPUE, not the season-wide one",
      !isTRUE(all.equal(sp$month_catch, 86 * sp$month_effort / 100)))
  chk("B44 split: no catch strata for the group -> catch 0 in every month, effort intact",
      { z <- pe_monthly_split(list(effort_strata = es), cal, "Red_Herring")
        isTRUE(all.equal(z$month_catch, c(0, 0))) && isTRUE(all.equal(z$month_effort, c(30, 70))) })
  chk("B44 split: closed days are not in the calendar a cell spreads over",
      { c2 <- cal; c2$open_section_1[1] <- FALSE; es2 <- es; es2$n_total_days[1] <- 3
        z <- pe_monthly_split(list(effort_strata = es2), c2, "x"); isTRUE(all.equal(sum(z$month_effort), 100)) })
  for (f in c("03_R_functions/run_pe_pooled.R", "03_R_functions/run_pe_gear.R"))
    chk(sprintf("B44 split: %s keeps the stratum catch", basename(f)),
        any(grepl("results$catch_strata[[cg]] <- catch_strat |>", readLines(f, warn = FALSE), fixed = TRUE)))
  uses <- vapply(c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd",
                   "03_R_functions/save_run_diagnostics.R"),
                 function(f) length(grep("pe_monthly_split(pe_all[[pe_label]]", readLines(f, warn = FALSE), fixed = TRUE)),
                 integer(1))
  chk("B44 split: 7.8 and 7.8b (pooled), 7 (gear) and monthly_pe_vs_bss.csv all use it",
      identical(unname(uses), c(2L, 1L, 1L)), paste(uses, collapse = "/"))

  # ---- 81d. the PE effort variance carries the donor covariances --------------------------
  source("03_R_functions/pe_effort_strata.R")
  days <- data.frame(event_date = as.Date("2025-01-06") + 0:19)
  days$period <- rep(1:4, each = 5); days$day_type <- "weekday"; days$month <- 1; days$open_section_1 <- TRUE
  vals <- c(10, 14, 12, 20, 16)
  de <- data.frame(event_date = days$event_date[c(1, 2, 3, 6, 8)], section_num = 1,
                   period = c(1, 1, 1, 2, 2), day_type = "weekday", est_daily_effort = vals)
  st <- pe_build_effort_strata(de, days, list(pe_empty_effort_stratum = "local_day_type", pe_variance = "impute_aware"))
  s2 <- stats::var(vals); nG <- 5
  cells <- 25 * 4 / 3 + 25 * 8 / 2 + 2 * 25 * s2 * (1 / nG + 1)
  ii <- (10^2 - 2 * 25) * s2 / nG                         # the two imputed cells share one donor mean
  is_ <- 2 * 10 * (s2 / nG) * (5 * 3 / 3 + 5 * 2 / 2)       # ... built from the sampled cells' days
  chk("B44 var: the cell variances are the documented ones (the donor mean divides by ITS n)",
      isTRUE(all.equal(attr(st, "effort_var_cells"), cells)))
  chk("B44 var: imputed x imputed and imputed x sampled covariances are the algebra in the header",
      isTRUE(all.equal(unname(attr(st, "effort_var_cov")), c(ii, is_))))
  rep_ <- NULL; invisible(utils::capture.output(rep_ <- pe_effort_stratum_report(st, "fixture", list())))
  chk("B44 var: the component total's variance is cells + covariances, and the report reads it",
      isTRUE(all.equal(attr(st, "effort_var_total"), cells + ii + is_)) &&
      isTRUE(all.equal(rep_$effort_se, sqrt(cells + ii + is_))))
  chk("B44 var: with no imputed cell the total is the plain sum of cell variances (nothing moves)",
      { de2 <- rbind(de, data.frame(event_date = days$event_date[c(11, 12, 16, 17)], section_num = 1,
                                    period = c(3, 3, 4, 4), day_type = "weekday", est_daily_effort = c(9, 11, 13, 15)))
        st2 <- pe_build_effort_strata(de2, days, list(pe_variance = "impute_aware"))
        isTRUE(all.equal(attr(st2, "effort_var_total"), sum(st2$se_total^2))) && sum(attr(st2, "effort_var_cov")) == 0 })
  chk("B44 var: pe_variance = 'sampled_only' is the historical arithmetic, no covariance",
      { st3 <- pe_build_effort_strata(de, days, list(pe_variance = "sampled_only"))
        sum(attr(st3, "effort_var_cov")) == 0 && isTRUE(all.equal(attr(st3, "effort_var_total"), 25 * 4 / 3 + 25 * 8 / 2)) })
  # Monte Carlo (iid days, sd 3): the total is sum_d w_d y_d, each sampled day's weight
  # being its own cell's N/n plus N/n_G for each imputed cell borrowing its donor mean, so
  # its true variance is sum(w^2) sigma^2 = 80.83 sigma^2. The donor_mean_only estimator's
  # expectation is the same 80.83 sigma^2; the independent sum it replaced expects 40.83.
  .old_rng <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(44); sig <- 3
  mc <- t(replicate(400, { dd <- de; dd$est_daily_effort <- stats::rnorm(5, 14, sig)
    z <- pe_build_effort_strata(dd, days, list(pe_variance = "donor_mean_only"))
    c(total = sum(z$est_total), v = attr(z, "effort_var_total"), v_indep = attr(z, "effort_var_cells")) }))
  if (!is.null(.old_rng)) assign(".Random.seed", .old_rng, envir = globalenv())
  w <- c(rep(5/3 + 10/5, 3), rep(5/2 + 10/5, 2)); truth <- sum(w^2) * sig^2
  chk("B44 var: Monte Carlo, the estimator's variance is sum(w^2) sigma^2 (the algebra is the design's)",
      abs(stats::var(mc[, "total"]) / truth - 1) < 0.2, sprintf("(%.0f vs %.0f)", stats::var(mc[, "total"]), truth))
  chk("B44 var: Monte Carlo, the covariance-aware SE is unbiased for it; the independent sum was about half",
      abs(mean(mc[, "v"]) / truth - 1) < 0.1 && mean(mc[, "v_indep"]) / truth < 0.6,
      sprintf("(mean estimate %.0f and %.0f vs truth %.0f)", mean(mc[, "v"]), mean(mc[, "v_indep"]), truth))
  # ---- 81e. the Red Rock catch group is gone --------------------------------------------
  .live <- function(f) { x <- readLines(f, warn = FALSE); x[!grepl("^\\s*#", x)] }
  .files <- c("run_config.R", "run_estimation.R", list.files("03_R_functions", full.names = TRUE),
              list.files("01_BSS_models", pattern = "Rmd$", full.names = TRUE))
  chk("B44 Red Rock: no live code reads estimate_red_rock or builds a Red_Rock_Kept group",
      !any(grepl("estimate_red_rock|Red_Rock_Kept", unlist(lapply(.files, .live)))))
  chk("B44 Red Rock: run_config.R no longer defines the key",
      { e <- new.env(); sys.source("run_config.R", envir = e); is.null(e$run_config$estimate_red_rock) })
  chk("B44 Red Rock: the fished-gear evidence (red_rock_kept in the unfished-zero-catch guard) is kept",
      any(grepl(".zero_catch    = (dungeness_kept <= 0) & (red_rock_kept <= 0)", .live("03_R_functions/fetch_crab_data.R"), fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 82. B45 (2026-09-28): AN UNSAMPLED FLOAT 17-21 COUNT IS NOT ZERO, AND THE GEAR TRACK FITS
#     THE POOLED TRACK'S INTERVIEWS. Float 17-21 is counted only when a second sampler is
#     free; an unpaired Float 20 count now takes round(R_month x its own Float 20 count).
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/fetch_crab_data.R")
  tz <- "America/Los_Angeles"
  mk <- function(d, hm, g) data.frame(event_date = as.Date(d), survey_id = "S1", count_time = paste0(hm, ":00"),
                                       count_time_posix = as.POSIXct(paste(d, paste0(hm, ":00")), tz = tz),
                                       total_gear_count = g, stringsAsFactors = FALSE)
  # January: three days with both floats at every count (the pairs), one day with Float 17-21 once
  f20 <- rbind(mk("2025-01-06", c("10:00", "12:00"), c(10, 20)), mk("2025-01-07", c("10:00", "12:00"), c(30, 40)),
               mk("2025-01-08", c("10:00", "12:00", "14:00"), c(50, 20, 10)), mk("2025-01-09", "11:00", 40))
  f20 <- f20[order(f20$event_date, f20$count_time_posix), ]
  f20$count_sequence <- stats::ave(seq_len(nrow(f20)), f20$event_date, FUN = seq_along)
  f17 <- rbind(mk("2025-01-06", c("10:05", "12:05"), c(2, 4)), mk("2025-01-07", c("10:05", "12:05"), c(6, 8)),
               mk("2025-01-08", "12:10", 5))
  P <- list(shore_f17_fill = "ratio", shore_f17_ratio_min_pairs = 5)
  r <- shore_dock_counts(f20, f17, P, quiet = TRUE)
  R <- (2 + 4 + 6 + 8 + 5) / (10 + 20 + 30 + 40 + 20)       # five time-paired counts in January
  chk("B45 f17: an observed Float 17-21 count is kept where it pairs with a Float 20 count",
      identical(r$f17_source[r$event_date == as.Date("2025-01-08") & r$count_sequence == 2], "observed") &&
      r$f17_gear[r$event_date == as.Date("2025-01-08") & r$count_sequence == 2] == 5)
  chk("B45 f17: an unpaired count is round(R_month x its own Float 20 count), not 0",
      identical(r$f17_gear[r$event_date == as.Date("2025-01-08")], c(round(R * 50), 5, round(R * 10))) &&
      r$f17_gear[r$event_date == as.Date("2025-01-09")] == round(R * 40) &&
      all(r$f17_source[r$f17_source != "observed"] == "ratio (month)"))
  chk("B45 f17: the count stays an integer and is Float 20 + Float 17-21",
      all(r$count_quantity == round(r$count_quantity)) && isTRUE(all.equal(r$count_quantity, r$f20_gear + r$f17_gear)))
  chk("B45 f17: a month below shore_f17_ratio_min_pairs uses the all-pairs ratio",
      { r6 <- shore_dock_counts(f20, f17, list(shore_f17_ratio_min_pairs = 6), quiet = TRUE)
        all(r6$f17_source[r6$f17_source != "observed"] == "ratio (all pairs)") })
  chk("B45 f17: shore_f17_fill = 'zero' reproduces the pre-B45 counts",
      { z <- shore_dock_counts(f20, f17, list(shore_f17_fill = "zero"), quiet = TRUE)
        identical(z$f17_gear[z$f17_source != "observed"], rep(0, sum(z$f17_source != "observed"))) })
  chk("B45 f17: two Float 17-21 counts paired to one Float 20 count are averaged, not a duplicated row",
      { f17b <- rbind(f17, mk("2025-01-09", c("10:50", "11:10"), c(3, 5)))
        rb <- shore_dock_counts(f20, f17b, P, quiet = TRUE)
        nrow(rb) == nrow(f20) && rb$f17_gear[rb$event_date == as.Date("2025-01-09")] == 4 })
  chk("B45 f17: an unknown fill mode stops rather than guessing",
      inherits(try(shore_dock_counts(f20, f17, list(shore_f17_fill = "carry"), quiet = TRUE), silent = TRUE), "try-error"))
  chk("B45 f17: fetch_crab_data() builds the shore effort through shore_dock_counts(), and run_config ships 'ratio'",
      any(grepl("shore_effort <- shore_dock_counts(f20, f17, params) |>", readLines("03_R_functions/fetch_crab_data.R", warn = FALSE), fixed = TRUE)) &&
      { e <- new.env(); sys.source("run_config.R", envir = e); identical(e$run_config$shore_f17_fill, "ratio") &&
        identical(e$run_config$shore_f17_ratio_min_pairs, 5) })
  # the gear prep no longer drops interviews without a positive fishing time (Matt: keep on both)
  gp <- readLines("03_R_functions/prep_bss_crab_gear.R", warn = FALSE); gp <- gp[!grepl("^\\s*#", gp)]
  chk("B45 interviews: the gear prep filters on the effort unit's own column, as the pooled prep, not on fishing time",
      !any(grepl("filter(!is.na(fishing_time_total), fishing_time_total > 0)", gp, fixed = TRUE)) &&
      any(grepl("int_d <- int_d[is.finite(v) & v > 0, , drop = FALSE]", gp, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 83. B46 (2026-09-28): THE SECOND REVIEW BATCH. The gear R_G prior, the gate's NA verdict,
#     the RNG restored everywhere, the predictive catch from the observation model, year-keyed
#     weeks, the span guard, the output folder, the package loader, the OSP classified count,
#     the Stan-side OSP guard, locale-free weekdays, the orchestrator, the runners' guards.
# ---------------------------------------------------------------------------
local({
  rl <- function(f) { x <- readLines(f, warn = FALSE); x[!grepl("^\\s*#", x)] }
  flat <- function(x) paste(x, collapse = "\n")
  gstan <- flat(rl("02_stan_models/crab_bss_gear_resolved.stan")); pstan <- flat(rl("02_stan_models/crab_bss_pooled.stan"))
  # 1.11: the gear model's R_G prior is the pooled model's
  chk("B46 R_G: the gear Stan declares R_G_prior_mu / _sigma and uses them (the 2024-25 literal is gone)",
      grepl("R_G ~ lognormal(log(R_G_prior_mu), R_G_prior_sigma);", gstan, fixed = TRUE) &&
      !grepl("R_G ~ lognormal(log(1.3), 0.3);", gstan, fixed = TRUE) && grepl("real<lower=0> R_G_prior_mu;", gstan, fixed = TRUE))
  chk("B46 R_G: the gear prep resolves it exactly as the pooled prep (override, else the season's interview ratio)",
      grepl("R_G_prior_mu    = params$R_G_prior_mu %||% summ$empirical_R_G %||% 1.3,", flat(rl("03_R_functions/prep_bss_crab_gear.R")), fixed = TRUE))
  # 1.12: gate, adequacy, ladder
  g <- flat(rl("03_R_functions/bss_convergence_gate.R"))
  chk("B46 gate: an NA R-hat or n_eff FAILS (no na.rm that lets -Inf pass)",
      grepl("pass_rhat <- all(is.finite(c(rhat_C, rhat_E))) && max(rhat_C, rhat_E) < rhat_threshold", g, fixed = TRUE) &&
      !grepl("na.rm = TRUE) < rhat_threshold", g, fixed = TRUE))
  chk("B46 adequacy: the coverage message prints the coverage, keeping its sign",
      grepl("cov50_worst_value = cov_worst_value,", flat(rl("03_R_functions/bss_model_adequacy.R")), fixed = TRUE))
  for (drv in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"))
    chk(sprintf("B46 ladder: %s treats a gate that errors as NOT passed", basename(drv)),
        grepl(".passed      <- !is.null(gate_try) && isTRUE(gate_try$pass_convergence)", flat(rl(drv)), fixed = TRUE))
  gd <- flat(rl("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"))
  chk("B46 gear tau: every boat all-gear sub-season from the list, and the boat I/E count from the fit's Stan data",
      grepl('n_ie_boat   <- as.integer(b_boat$bss_data$IE_n %||% 0L)', gd, fixed = TRUE) &&
      !grepl('boat_key <- "private_boat_all_gear_Dungeness_Kept"', gd, fixed = TRUE))
  # 4.6: bss_with_seed
  set.seed(7); a <- runif(1)
  set.seed(7); v1 <- bss_with_seed(99, runif(3)); b <- runif(1)
  set.seed(99); v2 <- runif(3)
  chk("B46 RNG: bss_with_seed() draws what set.seed() would, then restores the caller's stream exactly",
      identical(v1, v2) && identical(a, b))
  chk("B46 RNG: no shared function or driver leaves set.seed() behind (every one sits inside a restore)",
      { hits <- unlist(lapply(c(list.files("03_R_functions", full.names = TRUE), list.files("01_BSS_models", pattern = "Rmd$", full.names = TRUE)),
                              function(f) { x <- grep("set\\.seed\\(", sub("#.*$", "", rl(f)), value = TRUE)
                                           if (length(x)) paste(basename(f), x) else character(0) }))
        allowed <- c("bss_rng.R", "bss_stan_fit.R", "diagnose_effort_overdispersion.R", "model_diagnostics.R",
                     "BSS-GH-pooled-CPUE-model.Rmd", "BSS-GH-gear-type-CPUE-model.Rmd")   # each of these restores the caller's RNG
        all(sub(" .*$", "", hits) %in% allowed) })
  # 2.6: the predictive catch
  source("03_R_functions/bss_rng.R"); source("03_R_functions/bss_predictive_catch.R")
  D <- 30; S <- 4000; mu <- matrix(40, S, D); E <- matrix(100, S, D); r <- rep(0.8, S); hb <- 4
  y0 <- bss_predictive_draws(mu, E, r, rep(0, S), hb, seed = 3)
  # sum over n = E/hbar = 25 parties per day of NB2(m, r): day var = mu (1 + m / r), m = mu / n = 1.6
  v_day <- 40 * (1 + 1.6 / 0.8)
  chk("B46 predictive: theta = 0, the season total has the variance of a sum of NB2 parties (not a Poisson's)",
      abs(stats::var(y0) / (D * v_day) - 1) < 0.1 && stats::var(y0) > 2.5 * D * 40, sprintf("(var %.0f vs %.0f; Poisson %.0f)", stats::var(y0), D * v_day, D * 40))
  y1 <- bss_predictive_draws(mu, E, r, rep(0.3, S), hb, seed = 3)
  v_day1 <- 40 * (1 + (40 / (25 * 0.7)) / 0.8 + 0.3 * (40 / (25 * 0.7)))
  chk("B46 predictive: with zero-inflation the variance carries the extra theta * m term; the mean stays mu",
      abs(mean(y1) / (D * 40) - 1) < 0.01 && abs(stats::var(y1) / (D * v_day1) - 1) < 0.1)
  for (drv in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"))
    chk(sprintf("B46 predictive: %s builds Predictive_Catch from bss_predictive_catch()", basename(drv)),
        grepl("C_pred_draws <- tryCatch(bss_predictive_catch(fit, bss_data, bss_predictive_seed(params$bss_seed, label)),", flat(rl(drv)), fixed = TRUE) &&
        grepl("bss_C_pred_total <- bss_C_pred_total + (b$C_pred_draws %||% b$C_draws)[idx]", flat(rl(drv)), fixed = TRUE))
  # 3.1: year-keyed ISO weeks and the span guard
  source("03_R_functions/prep_days_crab.R")
  # the calendar only; day length (suncalc) is not what these assertions are about
  bss_assign_day_length <- function(days, L_eff_model, params) days
  environment(prep_days_crab) <- environment()
  dd <- prep_days_crab("2024-12-28", "2025-01-08", list(days_wkend = c("Saturday", "Sunday"), crabbing_holiday_dates = as.Date(character()),
                                                        period_pe = "week", sections = 1))
  chk("B46 weeks: the New Year week (Mon 30 Dec to Sun 5 Jan) is ONE period and ONE AR week",
      length(unique(dd$period[dd$event_date >= as.Date("2024-12-30") & dd$event_date <= as.Date("2025-01-05")])) == 1 &&
      length(unique(dd$week_index[dd$event_date >= as.Date("2024-12-30") & dd$event_date <= as.Date("2025-01-05")])) == 1)
  dd2 <- prep_days_crab("2024-09-01", "2025-09-15", list(days_wkend = c("Saturday", "Sunday"), crabbing_holiday_dates = as.Date(character()),
                                                         period_pe = "week", sections = 1))
  chk("B46 weeks: the same week number in two years is two periods (a window opening before mid-September)",
      dd2$period[dd2$event_date == as.Date("2024-09-09")] != dd2$period[dd2$event_date == as.Date("2025-09-08")])
  source("03_R_functions/validate_season_window.R")
  eff <- data.frame(date = as.Date("2024-10-01") + 0:5, season = "2024-25"); int <- data.frame(event_date = as.Date("2024-10-01"), season = "2024-25")
  base_p <- list(est_date_start = "2024-09-16", est_date_end = "2025-09-15", season_filter = "2024-25")
  chk("B46 span: two seasons in season_filter with pot_closures NULL STOP",
      inherits(try(suppressWarnings(validate_season_window(eff, int, modifyList(base_p, list(season_filter = c("2024-25", "2025-26"))), quiet = TRUE)), silent = TRUE), "try-error"))
  chk("B46 span: a window longer than a season with pot_closures NULL STOPS; the single shipped season passes",
      inherits(try(suppressWarnings(validate_season_window(eff, int, modifyList(base_p, list(est_date_end = "2026-09-15")), quiet = TRUE)), silent = TRUE), "try-error") &&
      isTRUE(suppressWarnings(validate_season_window(eff, int, base_p, quiet = TRUE))))
  # 3.2: the output folder
  source("03_R_functions/bss_output_dir.R")
  td <- tempfile("out"); dir.create(td)
  rc <- list(run_tag = "canonical-2024-25", season_filter = "2025-26")
  chk("B46 folder: a tag naming a season the run does not model STOPS",
      inherits(try(bss_output_dir("pooled-CPUE-", rc, base = td), silent = TRUE), "try-error"))
  d1 <- bss_output_dir("pooled-CPUE-", list(run_tag = "season-2025-26", season_filter = "2025-26"), base = td)
  writeLines("x", file.path(d1, "f.csv"))
  d2 <- bss_output_dir("pooled-CPUE-", list(run_tag = "season-2025-26", season_filter = "2025-26"), base = td)
  chk("B46 folder: a second render under one tag never writes into the filled folder",
      !identical(d1, d2) && startsWith(basename(d2), basename(d1)) && dir.exists(d2))
  chk("B46 folder: a span tag (first start year to last end year) is accepted for a multi-season run",
      dir.exists(bss_output_dir("p-", list(run_tag = "two-season-2024-26", season_filter = c("2024-25", "2025-26")), base = td)))
  chk("B46 folder: the tag is made path-safe",
      basename(suppressMessages(bss_output_dir("p-", list(run_tag = "gear_resolved_G = FALSE", season_filter = "2024-25"), base = td))) == "p-gear_resolved_G-FALSE")
  unlink(td, recursive = TRUE)
  # 4.2 / 4.3 / 4.4: loader, gear driver, orchestrator
  for (f in c("run_estimation.R", "01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"))
    chk(sprintf("B46 packages: %s loads through bss_load_packages() (no bare install.packages / require)", basename(f)),
        grepl("bss_load_packages()", flat(rl(f)), fixed = TRUE) && !grepl("sapply(load.lib, require", flat(rl(f)), fixed = TRUE))
  chk("B46 packages: renv.lock is the full closure (the headers rstan compiles against, writexl) and carries no unused package",
      { lk <- jsonlite::fromJSON("renv.lock")$Packages
        all(c("BH", "RcppEigen", "StanHeaders", "rstan", "writexl", "renv", "suncalc") %in% names(lk)) && !any(c("gt", "patchwork", "mgcv") %in% names(lk)) })
  chk("B46 gear driver: CPUE diagnostics wrapped, season totals written, session_info.txt written",
      grepl("tryCatch(write_cpue_diagnostics(b, label, output_dir),", gd, fixed = TRUE) &&
      grepl('utils::write.csv(season_totals, file.path(output_dir, "season_totals.csv"), row.names = FALSE)', gd, fixed = TRUE) &&
      grepl('writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))', gd, fixed = TRUE))
  re <- flat(rl("run_estimation.R"))
  chk("B46 orchestrator: 'both' renders the two models and writes the cross-check; unknown flags stop; a failure is recorded before the exit",
      grepl('models <- if (identical(model, "both")) c("pooled", "gear_resolved") else model', re, fixed = TRUE) &&
      grepl("Unknown argument", re, fixed = TRUE) && grepl("cross_check_%s.csv", re, fixed = TRUE) &&
      regexpr("write_manifest(stages, .base)", re, fixed = TRUE) < regexpr('stop(sprintf("Model stage(s) failed', re, fixed = TRUE))
  chk("B46 appendix: the manual divergence chunk is opt-in",
      any(grepl("{r divergence-diagnostics, echo=FALSE, eval=isTRUE(params$divergence_appendix)}", readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE), fixed = TRUE)))
  # 3.5: the report names the modelled seasons
  chk("B46 report: neither report hard-codes a season in its title or banner",
      !any(grepl("2024-25 Season Results Report", readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE), fixed = TRUE)) &&
      !grepl("SEASON SUMMARY: 2024-25", gd, fixed = TRUE) && !grepl('as.Date("2024-01-01"), as.Date("2024-12-31")', flat(rl("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd")), fixed = TRUE))
  # 3.7: OSP
  if (requireNamespace("writexl", quietly = TRUE) && requireNamespace("readxl", quietly = TRUE)) {
    source("03_R_functions/read_input_workbook.R"); source("03_R_functions/fetch_osp_boat_counts.R")
    source("03_R_functions/osp_sampling_rates.R")
    # a scratch workbook in 04_input_files/ (the reader resolves its file there), removed at once
    tf <- sprintf("zz_harness_osp_%d.xlsx", Sys.getpid()); tp <- file.path("04_input_files", tf)
    ow <- data.frame(Year = 2024, Month = 10, Day = 1:4, WestportPrivateEffort = c(40, 60, 20, 10),
                     WestportCrabOnlyEffort = c(0.25, 0.5, 0.1, 0.4), WestportCrabClassified = c(20, 30, 20, 10))
    writexl::write_xlsx(list(Sheet1 = ow), tp)
    op <- list(osp_boat_counts_file = tf, est_date_start = "2024-09-16", est_date_end = "2025-09-15")
    q <- function(pp) { r <- NULL; utils::capture.output(r <- tryCatch(fetch_osp_boat_counts(pp), error = function(e) e)); r }
    r_count <- q(op); r_frac <- q(modifyList(op, list(osp_crab_only_unit = "fraction")))
    unlink(tp)
    chk("B46 OSP (functional): a share delivered in the crab-only column STOPS under the default unit 'count'",
        inherits(r_count, "error") && grepl("looks like a SHARE", conditionMessage(r_count)))
    cr <- if (!inherits(r_frac, "error")) attr(r_frac, "osp_crab_rows") else NULL
    chk("B46 OSP (functional): under 'fraction' the share becomes a count OF THE CLASSIFIED boats, which is the binomial n",
        !is.null(cr) && identical(as.numeric(cr$osp_total), c(20, 30, 20, 10)) &&
        identical(as.numeric(cr$osp_crab_only), c(5, 15, 2, 4)))
  }
  of <- flat(rl("03_R_functions/fetch_osp_boat_counts.R"))
  chk("B46 OSP: a fraction-valued crab-only column read as a count STOPS; osp_crab_only_unit = 'fraction' converts it",
      grepl("it looks like a SHARE of the boats, not a count", of, fixed = TRUE) &&
      grepl('osp_crab_only = if (identical(crab_unit, "fraction")) round(osp_crab_only * osp_n)', of, fixed = TRUE))
  chk("B46 OSP: the classified-boat count, when delivered, is the binomial n of the crab-only share",
      grepl("sn <- osp_resolve_sample_n(cr$osp_boat_total, cr$osp_checked, cr$osp_rate, params)", of, fixed = TRUE) &&
      grepl("osp_total      = osp_n,", of, fixed = TRUE))   # B48: n resolved per day (osp_sampling_rates.R)
  for (st in list(pstan, gstan))
    chk("B46 Stan: the OSP crab-only share can never be fitted as f (reject without the combo walk)",
        grepl("if (osp_crab_lower == 1 && n_f_dyn == 1 && n_c_dyn == 0 && OSPF_n > 0)", st, fixed = TRUE))
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("B46 config: the OSP file keys and the sampling-frequency keys are on the control surface; ie_min_obs_boat only there",
      all(c("osp_boat_counts_file", "osp_boat_counts_sheet", "osp_dupe_resolve", "osp_match_trailer_area", "osp_crab_checked_col", "osp_crab_only_unit") %in% names(rc)) &&
      !grepl("ie_min_obs_boat  = 2,", gd, fixed = TRUE) && is.null(rc$fishery_name))
  chk("B46 weekdays: bss_weekday() gives English names from the ISO number, whatever the locale",
      identical(bss_weekday(as.Date(c("2025-01-04", "2025-01-06", "2025-01-05"))), c("Saturday", "Monday", "Sunday")) &&
      !any(grepl("(^|[^_.a-z])weekdays\\(", unlist(lapply(c("prep_days_crab.R", "crab_fraction.R", "estimate_comm_charter.R", "classify_day_type.R"),
                                                                  function(f) rl(file.path("03_R_functions", f)))))))
  chk("B46 PE: the scale check warns and flags; it no longer stops the run",
      !grepl('stop(sprintf(paste0("run_pe_pooled(): PE implied CPUE', flat(rl("03_R_functions/run_pe_pooled.R")), fixed = TRUE) &&
      !grepl('stop(sprintf(paste0("run_pe(): PE implied CPUE', flat(rl("03_R_functions/run_pe_gear.R")), fixed = TRUE))
  # 4.5: runners
  sr <- flat(rl("03_R_functions/bss_superseded_runner.R"))
  chk("B46 guard: the override is consumed (one override no longer disarms every runner)",
      grepl('suppressWarnings(rm("I_KNOW_THIS_IS_SUPERSEDED", envir = globalenv()))', sr, fixed = TRUE))
  chk("B46 guard: every guarded runner finds the guard with here::here() and FAILS CLOSED",
      all(vapply(list.files("06_diagnostics", pattern = "^run_.*[.]R$", full.names = TRUE), function(f) {
        x <- flat(rl(f)); !grepl("bss_superseded_runner(", x, fixed = TRUE) ||
          (grepl('.sr <- here::here("03_R_functions", "bss_superseded_runner.R")', x, fixed = TRUE) &&
           grepl("refusing to run", x, fixed = TRUE) && !grepl('if (exists("bss_superseded_runner"))', x, fixed = TRUE)) }, logical(1))))
  chk("B46 runners: both live ladders stamp and check an inputs fingerprint; none hides a failed source()",
      grepl("inputs: %s", flat(rl("06_diagnostics/run_improvements_2026-09-08.R")), fixed = TRUE) &&
      grepl("inputs: %s", flat(rl("06_diagnostics/run_marine_hazard_batch_2026-09-25.R")), fixed = TRUE) &&
      !any(vapply(setdiff(list.files("06_diagnostics", pattern = "[.]R$", full.names = TRUE),
                          "06_diagnostics/test_improvements_2026-08-25.R"),   # this file names the pattern it bans
                  function(f) grepl("try(source(f), silent = TRUE)", flat(rl(f)), fixed = TRUE), logical(1))))
})

# ---------------------------------------------------------------------------
# 84. B48 (2026-09-28): OSP SAMPLES EVERY k-TH BOAT, SO THE CRAB-ONLY SHARE IS OUT OF THE
#     BOATS SAMPLED. Erica (OSP): the rate is fixed for the day from the anticipated effort
#     and the staff; the manual's schedule gives the minimum by count. The binomial n is the
#     sampled count, else the day's rate x total, else the schedule's rate x total.
# ---------------------------------------------------------------------------
local({
  source("03_R_functions/read_input_workbook.R"); source("03_R_functions/osp_sampling_rates.R")
  source("03_R_functions/fetch_osp_boat_counts.R")
  tab <- read_osp_sampling_rates(list())
  chk("B48 schedule: the committed workbook is the manual's table, exact rationals, contiguous open-ended bands",
      nrow(tab) == 7 && isTRUE(all.equal(tab$rate, c(1, 4/5, 2/3, 1/2, 2/5, 1/3, 1/4))) &&
      identical(tab$count_min, c(0, 30, 51, 76, 101, 151, 201)) && is.na(tab$count_max[7]))
  cnt <- c(0, 29, 30, 50, 51, 75, 76, 100, 101, 150, 151, 200, 201, 500)
  chk("B48 schedule: every band edge takes the manual's rate (<30 all, 30-50 80%, ..., >200 25%)",
      isTRUE(all.equal(osp_rate_for_count(cnt, tab), c(1, 1, .8, .8, 2/3, 2/3, .5, .5, .4, .4, 1/3, 1/3, .25, .25))))
  tot <- c(40, 120, 250, 20)
  a <- osp_resolve_sample_n(tot, classified = c(33, NA, NA, NA), rate = c(NA, 0.5, NA, NA), list(), tab)
  chk("B48 n: the sampled count wins, then the day's rate x total, then the schedule's minimum rate x total",
      identical(a$n, c(33, 60, round(250 / 4), 20)) &&
      identical(a$source, c("sampled count", "rate column", "schedule (minimum rate)", "schedule (minimum rate)")))
  chk("B48 n: 'schedule' ignores the columns; 'none' restores n = the day's total (the pre-B48 reading)",
      identical(osp_resolve_sample_n(tot, c(33, NA, NA, NA), c(NA, .5, NA, NA), list(osp_sampling_rate_source = "schedule"), tab)$n,
                c(32, 48, round(250 / 4), 20)) &&
      identical(osp_resolve_sample_n(tot, NA, NA, list(osp_sampling_rate_source = "none"), tab)$n, tot))
  chk("B48 n: 'column' stops on a day with neither a sampled count nor a rate, rather than guessing",
      inherits(tryCatch(osp_resolve_sample_n(tot, NA, c(.5, NA, NA, NA), list(osp_sampling_rate_source = "column"), tab),
                        error = function(e) e), "error"))
  chk("B48 rate column: a percent column (any value above 1) is divided by 100; a value outside (0, 100] stops",
      isTRUE(all.equal(suppressMessages(osp_rate_as_fraction(c(50, 80, NA))), c(.5, .8, NA))) &&
      identical(osp_rate_as_fraction(c(.5, 1)), c(.5, 1)) &&
      inherits(tryCatch(osp_rate_as_fraction(c(0, 50)), error = function(e) e), "error"))
  if (requireNamespace("writexl", quietly = TRUE)) {
    tf <- sprintf("zz_harness_osp48_%d.xlsx", Sys.getpid()); tp <- file.path("04_input_files", tf)
    ow <- data.frame(Year = 2025, Month = 6, Day = 1:4, WestportPrivateEffort = c(40, 120, 250, 20),
                     crabbing_only = c(8, 12, 10, 5), WestportPrivateSampleRate = c(NA, 50, NA, NA))
    writexl::write_xlsx(list(Sheet1 = ow), tp)
    op <- list(osp_boat_counts_file = tf, est_date_start = "2024-09-16", est_date_end = "2025-09-15")
    q <- function(pp) { r <- NULL; utils::capture.output(r <- suppressMessages(tryCatch(fetch_osp_boat_counts(pp), error = function(e) e))); r }
    r1 <- q(op); r2 <- q(modifyList(op, list(osp_crab_only_basis = "expanded")))
    unlink(tp)
    c1 <- if (!inherits(r1, "error")) attr(r1, "osp_crab_rows") else NULL
    c2 <- if (!inherits(r2, "error")) attr(r2, "osp_crab_rows") else NULL
    chk("B48 reader (functional): OSP's 'crabbing_only' column is found, and each day's n is the boats SAMPLED",
        !is.null(c1) && identical(as.numeric(c1$osp_total), c(32, 60, 62, 20)) &&
        identical(as.numeric(c1$osp_crab_only), c(8, 12, 10, 5)) &&
        identical(c1$osp_rate_source, c("schedule (minimum rate)", "rate column", "schedule (minimum rate)", "schedule (minimum rate)")) &&
        identical(as.numeric(r1$count_quantity), c(40, 120, 250, 20)))
    chk("B48 reader (functional): an EXPANDED crab-only column is converted back to the sampled count (x the day's rate)",
        !is.null(c2) && identical(as.numeric(c2$osp_crab_only), c(round(8 * .8), 6, round(10 * .25), 5)))
  }
  wbl <- NULL; utils::capture.output(wbl <- fetch_osp_boat_counts(list()))
  # 2026-10-02 (B64): OSP delivered the crab-only column, so the committed workbook now HAS one.
  # What must hold instead: it is OSP's WPTPrivateCrabOnly, the notes-derived WPTPrivateCrabAlso
  # is NOT in the workbook, the reader finds the column, and the effort series is untouched.
  wraw <- readxl::read_excel("04_input_files/WBL_boat_counts.xlsx", sheet = "Sheet1")
  chk("B64: the committed WBL workbook carries OSP's WPTPrivateCrabOnly and NOT the notes-derived WPTPrivateCrabAlso; the reader finds it; the effort series is unchanged",
      "WPTPrivateCrabOnly" %in% names(wraw) && !any(grepl("CrabAlso", names(wraw), fixed = TRUE)) &&
      sum(!is.na(wraw$WPTPrivateCrabOnly)) == 253 && nrow(wraw) == 306 &&
      nrow(attr(wbl, "osp_crab_rows")) > 0 && nrow(wbl) > 0 && all(wbl$count_type == "OSP Boat Count") &&
      all(wraw$WPTPrivateCrabOnly <= wraw$WestportPrivateEffort, na.rm = TRUE))
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  chk("B48 config: the sampling keys are on the control surface, 'crabbing_only' is a crab-only column name, and (A33) the stream ships on",
      all(c("osp_sample_rate_col", "osp_sampling_rate_source", "osp_sampling_rates_file", "osp_crab_only_basis") %in% names(rc)) &&
      identical(rc$osp_crab_only_col[1], "WPTPrivateCrabOnly") && "crabbing_only" %in% rc$osp_crab_only_col && identical(rc$osp_sampling_rate_source, "auto") && isTRUE(rc$use_osp_crab_lower))
  for (drv in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"))
    chk(sprintf("B48 %s writes osp_crab_only_daily.csv (the per-day n, rate and its source) when crab rows exist", basename(drv)),
        any(grepl('file.path(output_dir, "osp_crab_only_daily.csv")', readLines(drv, warn = FALSE), fixed = TRUE)))
  chk("B48 docs: the input README lists the schedule workbook and where it came from",
      grepl("osp_sampling_rates.xlsx", paste(readLines("04_input_files/README.md", warn = FALSE), collapse = " "), fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 85. B50 (2026-09-28, documented 2026-09-29): THE METHOD OF RECORD RE-RENDERED AT B44 TO B49,
#     WITH ITS CROSS-CHECK. Matt rendered run_estimation.R with model = "both" on run_config.R as
#     shipped (code 4828b76, committed ff750c4). Read back from the folders' own files: the
#     totals the documents quote, the cross-check verdict, the claim that the boat fits did not
#     move (byte-identical to the 2026-09-27 render), the season totals equal to the port total
#     (B43 in production), every fit through the gate, and the documents moved to it.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  flat <- function(x) gsub("[ \n>]+", " ", x)
  N <- "05_output/20260928/pooled-CPUE-2024-25"; G <- "05_output/20260928/gear-type-CPUE-model-2024-25"
  O <- "05_output/20260927/pooled-CPUE-canonical-2024-25"
  row <- function(dir, lab) { d <- utils::read.csv(file.path(dir, "port_total_Dungeness_Kept.csv"), stringsAsFactors = FALSE)
                              d[d[[2]] == lab, , drop = FALSE] }
  if (dir.exists(N) && dir.exists(G)) {
    e <- row(N, "Expected_Catch"); g <- row(G, "Expected_Catch")
    chk("85 render: pooled 99,873 [82,414, 124,438], gear 98,588 [81,046, 123,563], PE 88,758, from the folders' own CSVs",
        identical(as.numeric(c(e$BSS_median, e$BSS_lo95, e$BSS_hi95, e$PE)), c(99873, 82414, 124438, 88758)) &&
        identical(as.numeric(c(g$BSS_median, g$BSS_lo95, g$BSS_hi95)), c(98588, 81046, 123563)))
    cc <- utils::read.csv("05_output/20260928/cross_check_20260928_165651.csv", stringsAsFactors = FALSE)
    chk("85 render: the cross-check file reads -1.29% against a 2% tolerance, PASS",
        identical(cc$verdict, "PASS") && isTRUE(abs(cc$gear_minus_pooled_pct - (-1.29)) < 1e-9) && cc$tolerance_pct == 2)
    cv <- utils::read.csv(file.path(N, "convergence_report.csv"), stringsAsFactors = FALSE)
    cg <- utils::read.csv(file.path(G, "convergence_report.csv"), stringsAsFactors = FALSE)
    chk("85 render: every fit on both tracks passes the gate and reports BSS",
        all(cv$pass_convergence) && all(cv$method_selected == "BSS") && all(cg$pass_convergence) && all(cg$method_selected == "BSS"))
    sag <- cv[cv$fit == "shore_all_gear_Dungeness_Kept", ]
    chk("85 render (D33): the shore all-gear fit's divergence fraction is the 4.07% the documents quote, under the 5% backstop",
        isTRUE(abs(sag$divergence_fraction - 0.0407) < 1e-9) && sag$divergence_fraction < 0.05)
    st <- utils::read.csv(file.path(N, "season_totals.csv"), stringsAsFactors = FALSE)
    chk("85 render (B43 in production): season_totals.csv equals the port total on a single season",
        identical(as.numeric(c(st$BSS_median, st$BSS_lo95, st$BSS_hi95)), as.numeric(c(e$BSS_median, e$BSS_lo95, e$BSS_hi95))))
    if (!dir.exists(O)) skp("85 render: boat fits byte-identical to 2026-09-27", "folder absent")
    if (dir.exists(O)) {
      same <- function(f) identical(unname(tools::md5sum(file.path(N, f))), unname(tools::md5sum(file.path(O, f))))
      bf <- c("bss_full_summary_private_boat_all_gear_Dungeness_Kept.csv", "bss_full_summary_private_boat_ring_net_only_Dungeness_Kept.csv",
              "bss_draws_summed_private_boat_all_gear_Dungeness_Kept.csv", "bss_draws_summed_private_boat_ring_net_only_Dungeness_Kept.csv")
      chk("85 render: both boat fits are byte-identical to the 2026-09-27 render (B44 to B49 inert on the boat; B38 across renders)",
          all(vapply(bf, same, logical(1))) && !same("bss_full_summary_shore_all_gear_Dungeness_Kept.csv"))
    }
    mf <- readLines("05_output/20260928/run_manifest_20260928_165651.txt", warn = FALSE)
    chk("85 render: the manifest records code 4828b76, the one config file that differed, and R 4.2.2 / rstan 2.32.7",
        any(grepl("^git sha\\s*: 4828b76", mf)) && any(grepl("run_config.R", mf, fixed = TRUE)) &&
        any(grepl("R version 4.2.2", mf, fixed = TRUE)) && any(grepl("rstan_2.32.7", mf, fixed = TRUE)))
  } else chk("85 render: the committed 20260928 run folders are present (the documents cite them)", FALSE,
             "missing; they are committed, so their absence is a checkout problem, not a skip")
  ps <- rd("07_documentation/development_notes/PIPELINE_STATUS.md"); cr <- rd("07_documentation/development_notes/CHANGE_REGISTER.md")
  vc <- rd("07_documentation/development_notes/VALIDATION_CAMPAIGN.md"); md <- rd("07_documentation/BSS-GH-pooled-CPUE-model-documentation.md")
  gd <- rd("07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md")
  # 2026-10-02 (B61): this run is no longer the authoritative one; the box and the register record it
  # as what the A31 render replaced, with its total, its cross-check and D33.
  chk("85 docs: the box records the 2026-09-28 run as what the authoritative run replaced, with its total, the cross-check and D33",
      grepl("`05_output/20260928/pooled-CPUE-2024-25`, 99,873 [82,414, 124,438]", flat(ps), fixed = TRUE) &&
      grepl("cross-check at -1.29%", flat(ps), fixed = TRUE) && grepl("D33", ps, fixed = TRUE))
  chk("85 docs: the register records the 2026-09-28 run as superseded, with B50 and D33",
      grepl("`05_output/20260928/pooled-CPUE-2024-25` (99,873 [82,414, 124,438], B50)", cr, fixed = TRUE) &&
      grepl("| B50 |", cr, fixed = TRUE) && grepl("| D33 |", cr, fixed = TRUE))
  chk("85 docs: the campaign has 1z.6 and the 1z anchor names ff750c4; the method document's Section 1 table carries the new run; the gear document its reference run",
      grepl("### 1z.6 The method of record re-rendered at B44 to B49", vc, fixed = TRUE) && grepl("`ff750c4` (the re-render at B44 to B49, 1z.6)", vc, fixed = TRUE) &&
      grepl("| **port total** | **99,873 [82,414, 124,438]** | **88,758** | **-11.1%** | |", md, fixed = TRUE) &&
      grepl("| **port total** | **98,588 [81,046, 123,563]** | **88,758** | **-10.0%** | | |", gd, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 86. B51 (2026-09-29): EVERY CHAIN STARTS WITHIN init_r = 0.5 (D33). rstan's default radius
#     is 2; at 2 a shore all-gear chain can start in the sigma_mu_E funnel and stay there. The
#     key is on the control surface, both drivers pass it to bss_stan_fit() (which forwards it
#     to rstan::stan(); checked there on a toy model when B51 landed), it is logged per fit,
#     and the sampler override accepts it for experiments.
# ---------------------------------------------------------------------------
local({
  e <- new.env(); sys.source("run_config.R", envir = e)
  chk("B51 config: run_config ships bss_init_r = 0.5 in the sampler section",
      identical(e$run_config$bss_init_r, 0.5))
  for (drv in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")) {
    x <- readLines(drv, warn = FALSE); x <- x[!grepl("^\\s*#", x)]
    chk(sprintf("B51 %s passes init_r = params$bss_init_r to bss_stan_fit() and prints it on the sampler line", basename(drv)),
        sum(grepl("init_r = params$bss_init_r %||% 2", x, fixed = TRUE)) == 1 &&
        any(grepl("init_r=%.2f", x, fixed = TRUE)) &&
        !any(grepl("bss_init_r *=", x)))                   # not in params_model: run_config owns it
  }
  source("03_R_functions/bss_sampler_override.R")
  chk("B51 override: bss_init_r is a sampler key the override accepts",
      grepl(.BSS_SAMPLER_OVERRIDE_PATTERN, "bss_init_r") && !grepl(.BSS_SAMPLER_OVERRIDE_PATTERN, "bss_init_radius"))
  ps <- paste(readLines("07_documentation/development_notes/PIPELINE_STATUS.md", warn = FALSE), collapse = " ")
  md <- paste(readLines("07_documentation/BSS-GH-pooled-CPUE-model-documentation.md", warn = FALSE), collapse = " ")
  chk("B51 docs: the box says the authoritative run was rendered at 2, the method document states the radius and the funnel limitation",
      grepl("init_r = 0.5", ps, fixed = TRUE) && grepl("**The initial-value radius is part of the configuration", md, fixed = TRUE) &&
      # 2026-09-29: limitation 15 now records the A31 collapse that closed D33
      grepl("**15. The single-section levels were unidentified, and the funnels they made are removed", md, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 87. D34 (2026-09-29): THE CENSUS IS CLIPPED TO THE ESTIMATION WINDOW. The same synthetic
#     two-week census as section 53; a window that covers it changes nothing, a window that
#     cuts it keeps the tally days inside only and says so, and a window that misses it gives
#     a zero component with a warning that says the zero is an absence, not an estimate.
# ---------------------------------------------------------------------------
local({
  cal <- seq(as.Date("2025-01-06"), as.Date("2025-01-19"), by = "day")
  wk  <- weekdays(cal) %in% c("Saturday", "Sunday")
  samp <- c(cal[!wk][1:8], cal[wk][1:2])
  tally <- tibble(date = samp, commercial_tally = c(4, 6, 5, 7, 4, 6, 5, 3, 8, 10), charter_tally = c(1, 0, 1, 1, 0, 1, 1, 0, 2, 2))
  ints  <- tibble(population = "comm_charter", event_date = rep(samp, each = 2),
                  boat_type_clean = rep(c("Commercial", "Charter"), 10), dungeness_kept = rep(c(40, 60), 10), red_rock_kept = 0)
  dwg <- list(comm_tally = tally, interview = ints)
  P0 <- list(census_start_date = "2025-01-06", census_end_date = "2025-01-19", days_wkend = c("Saturday", "Sunday"),
             crabbing_holiday_dates = as.Date(character()))
  q <- function(...) { r <- NULL
    utils::capture.output(r <- suppressWarnings(estimate_comm_charter(dwg, modifyList(P0, list(...)))))
    r }
  r_all  <- q()
  r_in   <- q(est_date_start = "2024-12-01", est_date_end = "2025-09-15")
  r_cut  <- q(est_date_start = "2025-01-13", est_date_end = "2025-09-15")
  r_out  <- q(est_date_start = "2025-02-01", est_date_end = "2025-03-01")
  keep <- tally$date >= as.Date("2025-01-13")
  chk("D34: a window containing the census changes nothing (no clip note)",
      isTRUE(all.equal(r_in$Dungeness_Kept, r_all$Dungeness_Kept)) && !any(grepl("CENSUS", r_in$frame_warnings)))
  chk("D34: a window cutting the census keeps only the tally days inside it, and a CLIPPED note leads the frame warnings",
      isTRUE(all.equal(r_cut$Dungeness_Kept, sum(tally$commercial_tally[keep] * 40 + tally$charter_tally[keep] * 60))) &&
      min(r_cut$daily_full$date) == as.Date("2025-01-13") && grepl("^CENSUS CLIPPED TO THE WINDOW", r_cut$frame_warnings[1]))
  chk("D34: a window missing the census gives 0 with an OUTSIDE note saying the zero is an absence, not an estimate",
      r_out$Dungeness_Kept == 0 && nrow(r_out$daily_full) == 0 &&
      grepl("^CENSUS OUTSIDE THE WINDOW", r_out$frame_warnings[1]) && grepl("not an estimate of zero harvest", r_out$frame_warnings[1], fixed = TRUE))
  chk("D34: census_windows entries are clipped the same way (the recursion calls the scalar path)",
      { rw <- q(census_windows = list(`2024-25` = c("2025-01-06", "2025-01-19")), est_date_start = "2025-01-13", est_date_end = "2025-09-15")
        isTRUE(all.equal(rw$Dungeness_Kept, r_cut$Dungeness_Kept)) })
  vs <- paste(readLines("03_R_functions/validate_season_window.R", warn = FALSE), collapse = "\n")
  chk("D34: the validator says a partly-outside census is clipped",
      grepl("it is clipped to the window (D34)", vs, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 88. A31 / T2.5 (2026-09-29): THE SINGLE-SECTION LEVEL HIERARCHIES ARE COLLAPSED ON THE POOLED
#     TRACK, per level (run_config mu_hier_collapse_single: "both" ships, the gear track's P2
#     structure; "effort" the effort level only; "none" the pre-A31 model). With S == 1 the pair
#     (mu_mu, sigma_mu * eps_mu) is unidentified and sigma_mu_E's funnel was D33.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  ps <- rd("02_stan_models/crab_bss_pooled.stan"); gs <- rd("02_stan_models/crab_bss_gear_resolved.stan")
  chk("A31 Stan: per-level collapse flags in transformed data, driven by mu_hier_collapse_single (0/1/2) and the POOL-4 lever",
      grepl("int<lower=0,upper=2> mu_hier_collapse_single;", ps, fixed = TRUE) &&
      grepl("int<lower=0, upper=1> use_mu_hier_E = (collapse_mu_hier == 0 && (S > 1 || mu_hier_collapse_single == 0)) ? 1 : 0;", ps, fixed = TRUE) &&
      grepl("int<lower=0, upper=1> use_mu_hier_C = (collapse_mu_hier == 0 && (S > 1 || mu_hier_collapse_single < 2)) ? 1 : 0;", ps, fixed = TRUE))
  chk("A31 Stan: eps_mu_E / eps_mu_C are zero-size when their level is collapsed, and their priors are guarded",
      grepl("matrix[G * use_mu_hier_E, S] eps_mu_E;", ps, fixed = TRUE) && grepl("matrix[G * use_mu_hier_C, S] eps_mu_C;", ps, fixed = TRUE) &&
      grepl("if (use_mu_hier_E == 1) eps_mu_E[g,s] ~ std_normal();", ps, fixed = TRUE) && grepl("if (use_mu_hier_C == 1) eps_mu_C[g,s] ~ std_normal();", ps, fixed = TRUE))
  chk("A31 Stan: sigma_mu_E / sigma_mu_C keep their proper priors (prior-only when collapsed)",
      grepl("sigma_mu_E ~ cauchy(0, value_cauchyDF_sigma_mu_E);", ps, fixed = TRUE) && grepl("sigma_mu_C ~ cauchy(0, value_cauchyDF_sigma_mu_C);", ps, fixed = TRUE))
  chk("A31 Stan: the gear track's P2 collapse is unchanged",
      grepl("int<lower=0, upper=1> use_mu_hier_E = (S > 1) ? 1 : 0;", gs, fixed = TRUE))
  dr <- function(pars, ...) bss_decoupled_reasons(pars, list(...))
  pp <- c("sigma_mu_E", "sigma_mu_C", "mu_mu_E[1]")
  chk("A31 reporting: 'effort' flags sigma_mu_E prior-only and not sigma_mu_C; 'both' (and the gear track, no key) flags both; 'none' neither; S > 1 neither unless forced",
      { e <- dr(pp, S = 1L, mu_hier_collapse_single = 1L); b <- dr(pp, S = 1L, mu_hier_collapse_single = 2L)
        g0 <- dr(pp, S = 1L); n <- dr(pp, S = 1L, mu_hier_collapse_single = 0L)
        m <- dr(pp, S = 2L, mu_hier_collapse_single = 2L); f <- dr(pp, S = 2L, collapse_mu_hier = 1L, mu_hier_collapse_single = 0L)
        !is.na(e[1]) && is.na(e[2]) && all(!is.na(b[1:2])) && all(!is.na(g0[1:2])) && all(is.na(n)) &&
          all(is.na(m)) && all(!is.na(f[1:2])) && all(is.na(c(e[3], b[3], f[3]))) })
  rc <- rd("run_config.R"); pr <- rd("03_R_functions/prep_bss_crab_pooled.R")
  chk("A31 config: mu_hier_collapse_single is a named level, validated by the prep (a typo stops); collapse_mu_hier ships FALSE",
      grepl('mu_hier_collapse_single    = "both",', rc, fixed = TRUE) && grepl("collapse_mu_hier           = FALSE,", rc, fixed = TRUE) &&
      grepl("mu_hier_collapse_single = bss_mu_hier_collapse_code(params$mu_hier_collapse_single),", pr, fixed = TRUE) &&
      grepl('code <- c(none = 0L, effort = 1L, both = 2L)[x]', pr, fixed = TRUE))
  v <- "05_output/t25_collapse_2026-09-29_validation.csv"
  if (!file.exists(v)) skp("A31 evidence: the container validation table", "file absent") else {
    x <- utils::read.csv(v, stringsAsFactors = FALSE)
    sh <- x[x$level == "both", ]
    # The stuck case on the old model: 23.3% divergent, a chain in the funnel, gate FAILED. At the
    # shipped level ("both") every refit must pass, and every shore all-gear refit, at either init
    # radius and seed, must be under 1% divergent: that is what "both" bought over "effort" (4.7%
    # at the shipped init_r = 0.5), and the reason it ships.
    chk("A31 evidence: all four fits refit at the shipped level, every one passes the gate, and every shore all-gear refit (including the seed and init_r = 2 that got stuck before) is under 1% divergent",
        all(c("shore_all_gear", "shore_ring_net_only", "private_boat_all_gear", "private_boat_ring_net_only") %in% sh$fit) &&
        all(sh$pass_convergence) && any(sh$fit == "shore_all_gear" & sh$init_r == 2) &&
        all(sh$divergence_fraction[sh$fit == "shore_all_gear"] < 0.01))
    chk("A31 evidence: the table carries the rejected level's shore all-gear refit at the shipped init_r (the evidence against \"effort\")",
        any(x$level == "effort" & x$fit == "shore_all_gear" & x$init_r == 0.5 & x$divergence_fraction > 0.03))
  }
  # B57 (2026-09-29): Matt's render of the B51 code (init_r = 0.5, pre-A31) is committed field
  # evidence for A31, and must never be taken for the authoritative run.
  b57 <- "05_output/20260928/pooled-CPUE-2024-25-222347/convergence_report.csv"
  if (!file.exists(b57)) chk("B57: the B51 render's convergence report is committed", FALSE) else {
    cr <- utils::read.csv(b57, stringsAsFactors = FALSE)
    chk("B57: in the B51 render the shore pot-closure fit failed its gate and reported its PE; the shore all-gear fit passed",
        startsWith(cr$method_selected[cr$fit == "shore_ring_net_only_Dungeness_Kept"], "PE") &&
        identical(cr$method_selected[cr$fit == "shore_all_gear_Dungeness_Kept"], "BSS"))
    gr <- utils::read.csv("05_output/20260929/gear-type-CPUE-model-2024-25/convergence_report.csv", stringsAsFactors = FALSE)
    chk("B57: the gear track (both levels collapsed since v6.0) passed every fit on the same render", all(gr$pass_convergence))
    st <- rd("07_documentation/development_notes/PIPELINE_STATUS.md")
    # 2026-10-02 (B61): the authoritative run is now the A31 render; the B51 render stays recorded
    # as never authoritative, in the register and the box
    crg <- rd("07_documentation/development_notes/CHANGE_REGISTER.md")
    chk("B57: the B51 render is never the authoritative run: the box names the A31 render and records B57's failure; the register's B57 row says not authoritative",
        grepl("**`05_output/20260929/pooled-CPUE-2024-25-220449`, port total 99,822", st, fixed = TRUE) &&
        grepl("which failed its gate in B57", gsub("[ \n>]+", " ", st), fixed = TRUE) && grepl("| B57 |", crg, fixed = TRUE))
  }
})

# ---------------------------------------------------------------------------
# 89. D36 / B53 (2026-09-29): THE TURNOVER FALLBACKS ARE THE UPDATED VALUES, AND LOUD. One home
#     (bss_tau_fallback), the 2024-25 values of the shipped methods, and every fallback lands
#     in the run-warnings log that the report prints.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  chk("D36: the fallback centres are 2.477 (shore, derived) and 3.03 (boat, trailer_mean_per_visit)",
      isTRUE(all.equal(bss_tau_fallback("shore", list()), 2.477)) && isTRUE(all.equal(bss_tau_fallback("private_boat", list()), 3.03)) &&
      isTRUE(all.equal(bss_tau_fallback("shore", list(tau_shore_prior_mu_fallback = 2)), 2)))
  rc <- rd("run_config.R")
  chk("D36: run_config ships the same two values",
      grepl("tau_shore_prior_mu_fallback = 2.477,", rc, fixed = TRUE) && grepl("tau_boat_prior_mu_fallback  = 3.03,", rc, fixed = TRUE))
  bss_warn_reset()
  pb <- suppressWarnings(bss_resolve_tau_boat_prior(list(tau_boat_prior_mu = "calibration"), NULL, quiet = TRUE))
  pb0 <- suppressWarnings(bss_resolve_tau_boat_prior(list(), NULL, quiet = TRUE))
  ps <- suppressWarnings(bss_resolve_tau_shore_prior(list(tau_shore_prior_mu = "derived", tau_shore_prior_sigma = "derived"), NULL, quiet = TRUE))
  lg <- bss_warn_log()
  chk("D36: a failed overlap centres the boat on 3.03, and an absent key means the shipped method (not the retired 1.2)",
      isTRUE(all.equal(pb$tau_boat_prior_mu, 3.03)) && grepl("^fallback", pb$tau_boat_prior_source) && isTRUE(all.equal(pb0$tau_boat_prior_mu, 3.03)))
  chk("D36: an unavailable derivation centres the shore on 2.477 with log-SD 0.3",
      isTRUE(all.equal(ps$tau_shore_prior_mu, 2.477)) && isTRUE(all.equal(ps$tau_shore_prior_sigma, 0.3)))
  chk("D36: every fallback is in the run-warnings log (three here), not only on the console",
      nrow(lg) == 3 && all(lg$source == "turnover") && all(lg$severity == "warning"))
  src <- paste(vapply(c(list.files("03_R_functions", full.names = TRUE), list.files("01_BSS_models", "[.]Rmd$", full.names = TRUE)), rd, ""), collapse = "\n")
  chk("D36: no code path defaults a turnover to the retired literals (1.7, 1.2) or to 2.7",
      !grepl("prior_mu(_fallback)?\\s*%\\|\\|%\\s*(1\\.7|1\\.2|2\\.7)\\b", src, perl = TRUE))
  pd <- rd("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd"); gd <- rd("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")
  chk("D36: both drivers route a failed OSP ingest or overlap diagnostic to the run-warnings log",
      all(vapply(list(pd, gd), function(t) grepl('bss_warn("OSP", paste("OSP overlap diagnostic skipped:"', t, fixed = TRUE) &&
                                                  grepl('bss_warn("OSP", paste("OSP ingest skipped:"', t, fixed = TRUE), logical(1))))
  bss_warn_reset()
  chk("run warnings: a note is recorded but not raised; a warning is both",
      { w <- tryCatch({ bss_warn("x", "a note", severity = "note"); "none" }, warning = function(w) "raised")
        w2 <- tryCatch({ bss_warn("x", "a warning"); "none" }, warning = function(w) "raised")
        identical(w, "none") && identical(w2, "raised") && nrow(bss_warn_log()) == 2 })
  bss_warn_reset()
})

# ---------------------------------------------------------------------------
# 90. D37 / B54 (2026-09-29): THE HOUSEKEEPING ITEMS.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  pd <- rd("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd"); gd <- rd("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")
  chk("D37(a): the dead params_model keys are gone from both drivers",
      !grepl("^\\s*max_divergences\\s*=", pd, perl = TRUE) && !grepl("max_divergence_distortion\\s*=", pd, perl = TRUE) &&
      !grepl("ie_use_all_dates_for_diagnostics\\s*=", pd, perl = TRUE) && !grepl("(?m)^\\s*max_divergences\\s*=", gd, perl = TRUE))
  src <- paste(vapply(list.files("03_R_functions", full.names = TRUE), rd, ""), collapse = "\n")
  chk("D37(b): tau_boat_prior_sigma defaults to 0.5 everywhere and shore_effort_unit to gear-deployments",
      !grepl("tau_boat_prior_sigma %||% 0.3", src, fixed = TRUE) && !grepl('shore_effort_unit %||% "crabber-hours"', src, fixed = TRUE))
  pg <- rd("03_R_functions/run_pe_gear.R")
  chk("D37(c): run_pe_gear filters the ratio frame on a positive angler count and returns pe_cpue_check",
      grepl(".ac > 0", pg, fixed = TRUE) && grepl("results$pe_cpue_check <-", pg, fixed = TRUE))
  re <- readLines("run_estimation.R", warn = FALSE)
  g <- grep('requireNamespace("here"', re, fixed = TRUE)[1]; l <- grep("^\\s*library\\(", re)
  chk("D37(d): run_estimation.R checks for here before any library() call",
      !is.na(g) && (!length(l) || all(l > g)))
  ga <- rd("06_diagnostics/gear_coverage_audit.R")
  thr_g <- regmatches(gd, regexpr("bss_min_gear_effective_n = [0-9]+", gd))
  chk("D37(e): gear_coverage_audit uses build_subseasons, writes under 05_output, and its threshold literal equals the gear driver's",
      grepl("build_subseasons(p)", ga, fixed = TRUE) && grepl('here::here("05_output"', ga, fixed = TRUE) &&
      !grepl('here::here("gear_coverage_audit.csv")', ga, fixed = TRUE) &&
      length(thr_g) == 1 && grepl(sprintf("bss_min_gear_effective_n %%||%% %s", sub(".*= ", "", thr_g)), ga, fixed = TRUE))
  h <- rd("06_diagnostics/test_improvements_2026-08-25.R")
  chk("D37(f): the harness counts skips apart and reports them on the summary line",
      grepl(paste0("skp <- function(nm, why = \"input absent\")"), h, fixed = TRUE) &&
      grepl(paste0("%d passed, %d failed, %d ", "skipped"), h, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 91. B55 (2026-09-29): THE DAY LENGTH IS COMPUTED ONLY WHERE IT IS USED. Under the
#     gear-deployment unit nothing reads L_effective (hours) or civil twilight, so neither is
#     computed and suncalc is not a required package; a time unit still gets both.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  chk("B55: a day length is needed only under a time unit, or when asked for",
      !bss_needs_day_length(list()) && !bss_needs_day_length(list(shore_effort_unit = "gear-deployments")) &&
      bss_needs_day_length(list(shore_effort_unit = "crabber-hours")) && bss_needs_day_length(list(shore_effort_unit = "gear-hours")) &&
      bss_needs_day_length(list(day_length_diagnostics = TRUE)))
  dd <- tibble(event_date = as.Date("2025-01-01") + 0:2, day_type = "weekday")
  a <- bss_assign_day_length(dd, NULL, list(shore_effort_unit = "gear-deployments"))
  chk("B55: under gear-deployments the day-length columns exist and are NA, and the source says why",
      all(c("day_length", "day_length_civil_twilight", "L_mu", "L_prior_sigma") %in% names(a)) && all(is.na(a$day_length)) &&
      grepl("not used", attr(a, "l_source"), fixed = TRUE))
  chk("B55: suncalc is no longer a required package (it stays in renv.lock for the time-unit path)",
      !"suncalc" %in% bss_required_packages)
  chk("B55: both drivers fit the L_effective regression only when a day length is needed",
      all(vapply(c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"), function(f)
        grepl("isTRUE(params$use_ie_day_length) && bss_needs_day_length(params) && nrow(ie_data) > 0", rd(f), fixed = TRUE), logical(1))))
  chk("B55: the ladder counts the in-window I/E days from shore_turnover_by_day.csv (older rungs: the L_effective detail)",
      grepl('rd(dir, "shore_turnover_by_day.csv") %||% rd(dir, "L_effective_ie_detail.csv")', rd("06_diagnostics/run_improvements_2026-09-08.R"), fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 92. B56 (2026-09-29): THE POOLED REPORT REBUILT. Each row pins a defect the 2026-09-28 render's
#     audit found, so it cannot come back quietly.
# ---------------------------------------------------------------------------
local({
  rmd <- paste(readLines("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", warn = FALSE), collapse = "\n")
  gd  <- paste(readLines("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd", warn = FALSE), collapse = "\n")
  chk("B56: both reports capture warnings and messages (hooks) instead of discarding them",
      all(vapply(list(rmd, gd), function(t) grepl("message = TRUE, warning = TRUE", t, fixed = TRUE) &&
        grepl('knitr::knit_hooks$set(warning = .record_condition("warning"), message = .record_condition("message"))', t, fixed = TRUE) &&
        !grepl("opts_chunk$set(echo = TRUE, message = FALSE", t, fixed = TRUE), logical(1))))
  chk("B56: both reports close with the run warnings and notes and write run_warnings.csv",
      grepl("## 16. Run warnings and notes", rmd, fixed = TRUE) && grepl("## Run warnings and notes", gd, fixed = TRUE) &&
      all(vapply(list(rmd, gd), function(t) grepl('file.path(output_dir, "run_warnings.csv")', t, fixed = TRUE), logical(1))))
  chk("B56: the configuration, inputs, turnover, fit-settings, component and structural tables exist",
      all(vapply(c("```{r config-show", "```{r inputs-show", "```{r turnover-show", "```{r fit-settings-show", "component_totals.csv",
                   "```{r adequacy-show", "```{r prior-influence-show", "```{r blockcv-show", "```{r sampler-show", "```{r tau-show",
                   "```{r f-show", "```{r run-info-show"), function(k) grepl(k, rmd, fixed = TRUE), logical(1))))
  chk("B56: the incomplete-trip table is printed (it was a non-final value and never rendered)",
      grepl("print(report_table(\n    dplyr::transmute(sens_df,", rmd, fixed = TRUE))
  chk("B56: the census string counts the COMMERCIAL interviews; the run summary does not read the orchestrator's model",
      grepl(".c$commercial_interviews %||% NA, format(round(.c$commercial_se", rmd, fixed = TRUE) &&
      !grepl('if (exists("model") && is.character(model)) model', rmd, fixed = TRUE))
  chk("B56: catch by mode reports the BSS draws beside the PE",
      grepl(".mode_draws <- function(pop)", rmd, fixed = TRUE) && grepl("BSS_median = vapply(.md, .q, numeric(1), p = 0.5)", rmd, fixed = TRUE))
  chk("B56: the input plots use the day's MEAN count and CPUE per gear deployment; the day-length figures only under a time unit",
      grepl("summarise(mean_gear = mean(count_quantity), n_counts = n()", rmd, fixed = TRUE) &&
      !grepl("summarise(total_gear = sum(count_quantity)", rmd, fixed = TRUE) &&
      grepl("Crab per gear deployment", rmd, fixed = TRUE) && !grepl("cpue = dungeness_kept / fishing_time_total", rmd, fixed = TRUE) &&
      grepl("if (bss_needs_day_length(params)) {\np_day_length <- days_full", rmd, fixed = TRUE))
  chk("B56: p-values are formatted as text so kable cannot round them to 0",
      grepl("fmt_p <- function(p)", rmd, fixed = TRUE) && !grepl("signif(p_raw, 3)", rmd, fixed = TRUE) && !grepl("round(t_p, 3)", rmd, fixed = TRUE))
  chk("B56: the OSP overlap plot and verdict use the calibration metric, not the max-per-day one",
      { o <- paste(readLines("03_R_functions/diagnose_osp_trailer_overlap.R", warn = FALSE), collapse = "\n")
        grepl('.metric <- params$tau_boat_calibration_metric %||% "trailer_mean_per_visit"', o, fixed = TRUE) &&
        !grepl('prim <- calibration |> dplyr::filter(trailer_metric == "trailer_max_per_day")', o, fixed = TRUE) })
  chk("B56: the PIT histogram bins stop at 1",
      grepl("breaks = seq(0, 1, by = 0.1)", paste(readLines("03_R_functions/model_diagnostics.R", warn = FALSE), collapse = "\n"), fixed = TRUE))
  chk("B56: no report sentence quotes a hard-coded 2024-25 result as this run's",
      !grepl("Half of every component's days rest on one sampled day", rmd, fixed = TRUE) &&
      !grepl("about -21% for pots", rmd, fixed = TRUE) && !grepl("pi ~ Dirichlet(observed gear catch + 0.5)", rmd, fixed = TRUE) &&
      !grepl("the run stops if any component exceeds 2x", rmd, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 73. THE HARNESS SIZE THE DOCUMENTS ADVERTISE (2026-09-13). Last, because the total is
#     only known here.
#
#     06_diagnostics/README.md (PULL_REQUEST.md until its merge on 2026-09-28) tells a reader how many assertions this file runs, and a reader
#     who runs it will see the real number. The dangerous direction is OVER-claiming, so
#     that is what fails: a quoted count may not exceed the count actually executed. It
#     may lag behind, because the count only ever grows and a lagging figure understates
#     rather than misleads; asserting equality would make every new assertion a
#     documentation edit, which is churn for no protection. A floor is still enforced
#     against the static call sites so the figure cannot rot arbitrarily far.
# ---------------------------------------------------------------------------
local({
  total <- ok + bad + 1L   # +1 for this assertion itself
  f <- "06_diagnostics/README.md"
  if (!file.exists(f)) { cat("NOTE  size: 06_diagnostics/README.md absent; the count check is skipped\n") } else {
    tf <- gsub("[ \n]+", " ", paste(readLines(f, warn = FALSE), collapse = "\n"))
    n  <- regmatches(tf, gregexpr("\\*\\*[0-9,]+ assertions", tf))[[1]]
    q  <- if (length(n)) as.numeric(gsub("[^0-9]", "", n[1])) else NA_real_
    sites <- length(grep("^\\s*chk\\(", readLines("06_diagnostics/test_improvements_2026-08-25.R",
                                                  warn = FALSE)))
    chk("size: 06_diagnostics/README.md does not over-claim the harness, and is not wildly stale",
        !is.na(q) && q <= total && q >= sites,
        sprintf("(quotes %s; executed %d; static call sites %d)",
                if (is.na(q)) "nothing" else formatC(q, format = "d", big.mark = ","),
                total, sites))
  }
})

# ---------------------------------------------------------------------------
# 93. B58 (2026-09-29): THE OVERNIGHT BATCH. The authoritative render at the A31 code and the
#     renders that close D29, D3 and D6, with their decision rules written before the run, and
#     the two driver additions the D29 rule reads (every ladder rung's block CV and trailer
#     coverage).
# ---------------------------------------------------------------------------
local({
  flat <- function(x) paste(x, collapse = "\n")
  rd <- function(p) readLines(p, warn = FALSE)
  rf <- "06_diagnostics/run_authoritative_batch_2026-09-29.R"
  r  <- readLines(rf, warn = FALSE); rt <- paste(r, collapse = "\n")
  chk("B58: the batch runner ships DRY_RUN <- TRUE and starts only on --go or BSS_BATCH_GO=1 (so the tracked tree stays clean)",
      identical(grep("^DRY_RUN <-", r, value = TRUE)[1], "DRY_RUN <- TRUE                    # ships TRUE; start with --go or BSS_BATCH_GO=1 (see above)") &&
      grepl('if ("--go" %in% commandArgs(trailingOnly = TRUE) || identical(Sys.getenv("BSS_BATCH_GO"), "1")) DRY_RUN <- FALSE', rt, fixed = TRUE))
  chk("B60: the stages are the authoritative render, D3, D6, the forced weekly and biweekly pooled rungs, the gear ladder, the radius check and the daily rung, in that order",
      grepl('STAGES  <- c("S0", "A", "D3", "D6", "D3R", "D6R", "D29W", "D29B", "D29G", "R2", "D29D")', rt, fixed = TRUE))
  chk("B60: the pooled D29 rungs are forced renders (ar_force), not an in-process ladder, which ran out of memory at its second rung",
      grepl('delta = list(ar_force = list(private_boat = list(all_gear = "weekly")))', rt, fixed = TRUE) &&
      grepl('delta = list(ar_force = list(private_boat = list(all_gear = "biweekly")))', rt, fixed = TRUE) &&
      !grepl("D29P = list(", rt, fixed = TRUE))
  chk("B58: stage A renders through run_estimation.R --model both (the production orchestrator, manifest and cross-check), not a copy of it",
      grepl('A    = list(model = "both", fit = TRUE, tag = BASE$run_tag, orchestrator = TRUE', rt, fixed = TRUE) &&
      grepl('status <- .sys(c(shQuote(ORCH), "--model", mdl), log, hours)', rt, fixed = TRUE))
  chk("B58 review: every stage runs under a wall-clock limit (system2 timeout), recorded as TIMED OUT",
      grepl("system2(RSCRIPT, args, stdout = log, stderr = log, timeout = round(3600 * hours))", rt, fixed = TRUE) &&
      grepl('identical(status, 124L)) "TIMED OUT"', rt, fixed = TRUE))
  chk("B58 review: a resumed stage A renders only the missing model and the runner writes the cross-check",
      grepl("PARTIAL RESUME", rt, fixed = TRUE) && grepl(".write_cross_check(dirs[[\"pooled\"]], dirs[[\"gear_resolved\"]])", rt, fixed = TRUE))
  chk("B58 review: a desk FAIL blocks only the stages it concerns (DEPENDS); an S0 FAIL blocks every stage",
      grepl("BLOCKED <- if (\"S0\" %in% .desk_fail) setdiff(STAGES, \"S0\") else", rt, fixed = TRUE) &&
      grepl("for (sid in setdiff(STAGES, c(\"S0\", BLOCKED)))", rt, fixed = TRUE))
  chk("B58 review: the manifest parser reads every model line after Stages: (a multi-line error must not hide the second)",
      grepl('if (!grepl("^  (pooled|gear_resolved)\\\\s", x)) next', rt, fixed = TRUE))
  chk("B58: every decision rule is written in the header before the run (A, R2, D3, D6, D29)",
      all(vapply(c("# A (can this render become the authoritative run?)", "# R2 (does A31 still need init_r = 0.5?)",
                   "# D3 (the gear track's AR period, per population)", "# D6 (the gear-track zero-inflated shore catch)",
                   "# D29 (the boat all-gear AR period)", "D29-4 THE DECISION, pooled track only"),
                 function(x) grepl(x, rt, fixed = TRUE), logical(1))))
  chk("B58: D3 and D29G put the shore at the POOLED caps and leave the boat at its shipped gear periods, by construction",
      grepl("MATCHED <- list(shore = BASE$ar_max_resolution$pooled$shore,", rt, fixed = TRUE) &&
      grepl("private_boat = list(all_gear = BASE$gear_period_bss$all_gear", rt, fixed = TRUE) &&
      grepl("delta = list(gear_period_bss = MATCHED)),", rt, fixed = TRUE))
  chk("B58: the D29 ladders are scoped to the boat all-gear fit, fit every rung, and ignore the cap under test",
      grepl('ar_escalate = list(private_boat = "all_gear"), ar_escalate_stop = "all_rungs"', rt, fixed = TRUE) &&
      grepl('ar_escalate_ladder = c("weekly", "biweekly", "monthly")', rt, fixed = TRUE) &&
      grepl("ar_escalate_respect_cap = FALSE", rt, fixed = TRUE))
  chk("B58: RESUME reuses a folder only on a digest over the whole resolved configuration, the code and the inputs",
      grepl("stage_digest <- function(sid) {", rt, fixed = TRUE) &&
      grepl(".cfg_text(resolve_cfg(sid)), CODE_FP, INPUTS_FP", rt, fixed = TRUE))
  chk("B58: the runner is recorded in the register (B58), the status document and the diagnostics README",
      grepl("| B58 |", flat(rd("07_documentation/development_notes/CHANGE_REGISTER.md")), fixed = TRUE) &&
      grepl("run_authoritative_batch_2026-09-29.R", flat(rd("07_documentation/development_notes/PIPELINE_STATUS.md")), fixed = TRUE) &&
      grepl("run_authoritative_batch_2026-09-29.R", flat(rd("06_diagnostics/README.md")), fixed = TRUE))
  for (drv in c("01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd", "01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd")) {
    t <- flat(rd(drv))
    chk(sprintf("B58: %s writes every ladder rung's block CV (prefix ladder_block) and logs trailer coverage per rung", basename(drv)),
        grepl('output_dir, prefix = "ladder_block"),', t, fixed = TRUE) && grepl("cov50_trailer", t, fixed = TRUE))
  }
  bcv <- flat(rd("03_R_functions/bss_block_cv.R"))
  chk("B58: write_block_cv_diagnostics() keeps the loo_block default prefix (the report reads it) and takes another",
      grepl('prefix = "loo_block") {', bcv, fixed = TRUE) &&
      grepl('sprintf("%s_%s_%s.csv", prefix, sn, label)', bcv, fixed = TRUE))
  source("03_R_functions/bss_rung_adequacy.R", local = TRUE)
  chk("B58: bss_rung_adequacy() carries cov50_trailer, NA when there is no fit",
      exists("bss_rung_adequacy") && "cov50_trailer" %in% names(bss_rung_adequacy(NULL, NULL)) &&
      is.na(bss_rung_adequacy(NULL, NULL)$cov50_trailer))
  # The decision code itself, run on synthetic folders: the runner's function definitions (and
  # only those, plus its constants) are evaluated into a sandbox, so nothing is fitted, sourced
  # from run_config.R or written to 05_output.
  ex <- parse(rf, keep.source = FALSE)
  sb <- new.env(parent = globalenv())
  for (f in c("bss_ar_resolution.R", "bss_block_cv.R", "loo_elpd_paired.R")) sys.source(file.path("03_R_functions", f), envir = sb)
  keep <- c("FITS", "FIT_SHORE_AG", "FIT_SHORE_PC", "FIT_BOAT_AG", "V", "REC", "D29", "D29_TABLE", "COMP", "COMP_NOTE", "%||%",
            ".nolib", "CHILD_ENV", "RSCRIPT")
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
                    (as.character(e[[2]]) %in% keep || (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))))
    eval(e, sb)
  sb$BASE <- list(cross_check_tolerance = 0.02)
  # B60 (2026-10-02): every pooled rung is a render's KEPT boat all-gear fit (stage A monthly,
  # D29W weekly, D29B biweekly, D29D daily), so each synthetic rung is a folder with the files a
  # production render writes for its kept fit.
  STAGE_OF <- c(monthly = "A", weekly = "D29W", biweekly = "D29B", daily = "D29D")
  mk1 <- function(r) {                        # r: list(res, gate, plf, bad, elpd, ran = res)
    d <- tempfile("d29_"); dir.create(d); f <- sb$FIT_BOAT_AG
    utils::write.csv(data.frame(fit = f, ar_resolution = r$ran %||% r$res, pass_convergence = r$gate, divergence_fraction = 0,
                                method_selected = if (isTRUE(r$gate)) "BSS" else "PE (convergence fail)"),
                     file.path(d, "convergence_report.csv"), row.names = FALSE)
    utils::write.csv(data.frame(fit = f, p_loo_frac = r$plf, n_pareto_bad = r$bad, p_loo_worst_stream = "osp"),
                     file.path(d, "model_adequacy.csv"), row.names = FALSE)
    utils::write.csv(data.frame(stream = c("trailer", "osp", "catch"), n_obs = c(195, 130, 175)),
                     file.path(d, sprintf("loo_summary_%s.csv", f)), row.names = FALSE)
    utils::write.csv(data.frame(C_expected_sum = 40000 + stats::rnorm(200, 0, 100)), file.path(d, sprintf("bss_draws_summed_%s.csv", f)), row.names = FALSE)
    if (!is.null(r$elpd))
      utils::write.csv(data.frame(data_type = "joint", block = sprintf("2025-W%02d", seq_along(r$elpd)), n_obs = 6, n_leaveout = 6, lpd_block = r$elpd,
                                  elpd_block = r$elpd, p_eff = 0.1, pareto_k = 0.3, reliable = TRUE, leaveout_streams = "trailer+osp"),
                       file.path(d, sprintf("loo_block_joint_%s.csv", f)), row.names = FALSE)
    d
  }
  decide <- function(rungs) {
    sb$V <- list(); sb$REC <- list(); sb$D29 <- list(); sb$D29_TABLE <- NULL
    dirs <- list(A = c(pooled = NA_character_), D29W = c(pooled = NA_character_), D29B = c(pooled = NA_character_), D29D = c(pooled = NA_character_))
    for (r in rungs) { sid <- STAGE_OF[[r$res]]; dirs[[sid]][["pooled"]] <- mk1(r); sb$d29_forced_row(dirs[[sid]][["pooled"]], sid, r$res) }
    sb$verdict_D29(dirs)
    sb$REC[[length(sb$REC)]]$recommendation
  }
  set.seed(58); base <- -10 + stats::rnorm(30, 0, 0.2)
  up <- base + 0.5 + stats::rnorm(30, 0, 0.1); same <- base + stats::rnorm(30, 0, 0.05); down <- base - 0.5 + stats::rnorm(30, 0, 0.1)
  R <- function(res, gate = TRUE, plf = 0.05, bad = 0, elpd = base) list(res = res, gate = gate, plf = plf, bad = bad, elpd = elpd)
  chk("B58 D29 rule: weekly eligible, adequate and BETTER than monthly -> MOVE the cap to weekly (D29-4a)",
      grepl("^MOVE the pooled boat all-gear cap to weekly", decide(list(R("weekly", elpd = up), R("biweekly", elpd = same), R("monthly")))))
  chk("B58 D29 rule: weekly fails the gate, biweekly BETTER -> MOVE to biweekly (a failed rung is decided, not undecided)",
      grepl("^MOVE the pooled boat all-gear cap to biweekly", decide(list(R("weekly", gate = FALSE, elpd = up), R("biweekly", elpd = up), R("monthly")))))
  chk("B58 D29 rule: no finer rung BETTER (NO EVIDENCE / WORSE) and monthly adequate -> KEEP monthly (D29-4b)",
      grepl("^KEEP monthly", decide(list(R("weekly", elpd = down), R("biweekly", elpd = same), R("monthly")))))
  chk("B58 D29 rule: an adequacy that was not computed is UNKNOWN and the rule says REVIEW, never KEEP (D29-4c)",
      grepl("^REVIEW", decide(list(R("weekly", plf = NA, elpd = down), R("biweekly", elpd = same), R("monthly")))))
  chk("B58 D29 rule: biweekly BETTER but weekly's block CV missing -> REVIEW (a finer rung undecided), not MOVE",
      grepl("^REVIEW", decide(list(R("weekly", elpd = NULL), R("biweekly", elpd = up), R("monthly")))))
  chk("B58 D29 rule: incomplete pooled rungs (monthly only) -> REVIEW (D29-0)",
      grepl("^REVIEW: the pooled rungs are incomplete", decide(list(R("monthly")))))
  chk("B60 D29 rule: a forced stage that ran at another period is not read as its rung, so the set is incomplete -> REVIEW",
      grepl("^REVIEW: the pooled rungs are incomplete", decide(list(c(R("weekly", elpd = up), list(ran = "monthly")), R("biweekly", elpd = same), R("monthly")))))
  chk("B60 D29 rule: the daily rung, inadequate, is decided (out), so it cannot block a BETTER weekly",
      grepl("^MOVE the pooled boat all-gear cap to weekly", decide(list(R("daily", plf = 0.26, bad = 12, elpd = up), R("weekly", elpd = up), R("biweekly", elpd = same), R("monthly")))))
  chk("B58 D29 rule: inadequate by bad Pareto k (> 5% of all 500 LOO observations) is OUT, so a BETTER weekly with 30 bad k does not move the cap",
      grepl("^KEEP monthly", decide(list(R("weekly", bad = 30, elpd = up), R("biweekly", elpd = same), R("monthly")))))
  # the three-valued helpers and the manifest parser
  chk("B58 review: a clause that could not be read is REVIEW, never PASS or FAIL",
      identical(sb$.tri(NA), "REVIEW") && identical(sb$.tri(TRUE), "PASS") && identical(sb$.tri(FALSE), "FAIL") &&
      identical(sb$.worst(c("PASS", "REVIEW")), "REVIEW") && identical(sb$.worst(c("REVIEW", "FAIL")), "FAIL") &&
      is.na(sb$.gate_all(NULL)) && identical(sb$.gate_all(data.frame(pass_convergence = c(TRUE, TRUE, TRUE))), FALSE))
  mf <- tempfile(fileext = ".txt")
  od <- tempfile("od_"); dir.create(file.path(od, "05_output", "20260930", "gear-type-CPUE-model-2024-25"), recursive = TRUE)
  writeLines(c("Run manifest", "git tree    : clean (the render ran on the committed tree)", "", "Stages:",
               "  pooled         FAILED after 3.1 min: line one of the error", "second line of the error   (partial folder: C:/x/05_output/20260930/pooled-CPUE-2024-25)",
               "  gear_resolved    35.9 min   C:\\\\x\\\\05_output\\\\20260930\\\\gear-type-CPUE-model-2024-25", "",
               "run_config (run-level overrides applied to the model; all 3 keys):", "  pooled  not a stage line"), mf)
  sb$.here <- function(...) file.path(od, ...)
  mp <- sb$.parse_manifest(mf)
  chk("B58 review: the manifest parser reads the gear line after a FAILED pooled line whose error ran over two lines",
      isTRUE(mp$pooled$failed) && identical(mp$gear_resolved$failed, FALSE) &&
      identical(basename(mp$gear_resolved$outdir %||% ""), "gear-type-CPUE-model-2024-25") && isTRUE(abs(mp$gear_resolved$minutes - 35.9) < 1e-9))
  # B59 (2026-09-30): every stage's R process uses THIS session's library with renv's autoloader
  # off. The first overnight attempt's stages activated the project's (incomplete) renv library
  # through .Rprofile while the Console used another, and all seven failed.
  chk("B59: the runner starts every stage with renv's autoloader off and this session's .libPaths()",
      identical(unname(sb$CHILD_ENV["RENV_ACTIVATE_PROJECT"]), "FALSE") &&
      identical(strsplit(unname(sb$CHILD_ENV["R_LIBS"]), .Platform$path.sep, fixed = TRUE)[[1]],
                normalizePath(.libPaths(), winslash = "/", mustWork = FALSE)) &&
      grepl("st <- tryCatch(.with_child_env(suppressWarnings(system2(RSCRIPT, args,", rt, fixed = TRUE))
  chk("B59: .with_child_env() sets the stage environment and restores this session's afterwards (set and unset alike)",
      { o1 <- Sys.getenv("RENV_ACTIVATE_PROJECT", unset = NA); o2 <- Sys.getenv("R_LIBS_SITE", unset = NA)
        Sys.setenv(RENV_ACTIVATE_PROJECT = "xyz"); Sys.unsetenv("R_LIBS_SITE")
        inside <- sb$.with_child_env(c(Sys.getenv("RENV_ACTIVATE_PROJECT"), Sys.getenv("R_LIBS_SITE")))
        after <- c(Sys.getenv("RENV_ACTIVATE_PROJECT"), Sys.getenv("R_LIBS_SITE", unset = "UNSET"))
        if (is.na(o1)) Sys.unsetenv("RENV_ACTIVATE_PROJECT") else Sys.setenv(RENV_ACTIVATE_PROJECT = o1)
        if (!is.na(o2)) Sys.setenv(R_LIBS_SITE = o2)
        identical(inside[1], "FALSE") && identical(inside[2], unname(sb$CHILD_ENV["R_LIBS_SITE"])) &&
          identical(after, c("xyz", "UNSET")) })
  chk("B59: a process started the way a stage is sees exactly this session's library paths",
      { out <- sb$.with_child_env(suppressWarnings(system2(sb$RSCRIPT, c("-e", shQuote("cat(normalizePath(.libPaths(), winslash = '/'), sep = '\\n')")),
                                                           stdout = TRUE, stderr = FALSE)))
        identical(out, normalizePath(.libPaths(), winslash = "/")) })
  chk("B59: desk check S0 starts such a process and compares loadability and versions before anything is fitted",
      grepl("each stage's R process loads this session's packages at the same versions", rt, fixed = TRUE) &&
      grepl('CHILD_PKGS <- unique(c(bss_required_packages, "StanHeaders"', rt, fixed = TRUE))
  bp <- flat(rd("03_R_functions/bss_packages.R"))
  chk("B59: bss_load_packages() asks whether a package is INSTALLED without loading it, and checks loading once, after any restore",
      grepl("installed <- function(p) vapply(p, function(x) nzchar(system.file(package = x)), logical(1))", bp, fixed = TRUE) &&
      grepl("missing <- pkgs[!installed(pkgs)]", bp, fixed = TRUE) && grepl("broken <- pkgs[!have(pkgs)]", bp, fixed = TRUE))
  chk("B59: desk check S0 compares this session's library with renv.lock, FAILs on a Stan-toolchain difference and compares as versions, not strings",
      grepl("this session's library holds renv.lock's versions", rt, fixed = TRUE) &&
      grepl('tool <- intersect(off, c("rstan", "StanHeaders", "Rcpp", "RcppEigen", "BH", "RcppParallel"))', rt, fixed = TRUE) &&
      grepl("package_version(hv[[p]]) == package_version(lk[[p]])", rt, fixed = TRUE) &&
      isTRUE(package_version("1.84.0") == package_version("1.84.0-0")) && !isTRUE(package_version("1.4.5") == package_version("1.4-8")))
  chk("B59: the lockfile parse the S0 check uses reads every package in renv.lock",
      { lt <- paste(readLines("renv.lock", warn = FALSE), collapse = "\n")
        m <- regmatches(lt, gregexpr('"Package": "[^"]+",\\s*"Version": "[^"]+"', lt))[[1]]
        n_json <- length(gregexpr('"Package": "', lt, fixed = TRUE)[[1]])
        length(m) == n_json && length(m) > 100 })
  chk("B59: the runner ships DRY_RUN <- TRUE again (the first attempt's edit to FALSE is reverted)",
      identical(grep("^DRY_RUN <-", r, value = TRUE)[1], "DRY_RUN <- TRUE                    # ships TRUE; start with --go or BSS_BATCH_GO=1 (see above)"))
  # B61 (2026-10-02): D3 and D6 re-rendered once at bss_seed + 1, with the rule stated first; D6
  # gains D6-0 (a converged baseline). Functional: the D6 verdict on synthetic folders whose
  # baseline gate failed reads UNDECIDED, whatever D6-2 says.
  chk("B61: D3R and D6R are D3 and D6 at bss_seed + 1, and D3R / D6R decide D3 / D6",
      grepl("delta = list(gear_period_bss = MATCHED, bss_seed = BASE$bss_seed + 1)),", rt, fixed = TRUE) &&
      grepl('delta = list(gear_period_bss = MATCHED, catch_zi_tracks = c("pooled", "gear_resolved"), bss_seed = BASE$bss_seed + 1)),', rt, fixed = TRUE) &&
      grepl('verdict_D3(DIRS$D3R[["gear_resolved"]], DIRS$A[["pooled"]], DIRS$A[["gear_resolved"]], sid = "D3R", item = "D3",', rt, fixed = TRUE) &&
      grepl('verdict_D6(DIRS$D6R[["gear_resolved"]], DIRS$D3R[["gear_resolved"]], sid = "D6R", item = "D6")', rt, fixed = TRUE) &&
      grepl("#   D6-0 (B61, stated before D6R ran) the BASELINE is valid", rt, fixed = TRUE))
  chk("B61 D6 rule: an unconverged D3 baseline makes D6 UNDECIDED, however D6-2 reads (D6-0)",
      { sb$V <- list(); sb$REC <- list(); f <- sb$FIT_SHORE_AG
        mkd <- function(pass) { d <- tempfile("d6_"); dir.create(d)
          utils::write.csv(data.frame(fit = sb$FITS, pass_convergence = c(TRUE, pass, TRUE, TRUE), method_selected = "BSS", divergence_fraction = 0),
                           file.path(d, "convergence_report.csv"), row.names = FALSE); d }
        sb$verdict_D6(mkd(TRUE), mkd(FALSE), sid = "D6R", item = "D6")
        r1 <- sb$REC[[length(sb$REC)]]$recommendation
        any(vapply(sb$V, function(v) identical(v$verdict, "FAIL") && startsWith(v$criterion, "D6-0"), logical(1))) &&
          startsWith(r1, "UNDECIDED: the D3 baseline") })
})

# ---------------------------------------------------------------------------
# 94. B61 (2026-10-02): THE AUTHORITATIVE RUN MOVES TO THE BATCH'S STAGE A, the first render at the
#     A31 code with init_r = 0.5, judged by rule A written before it ran. Its own files, its
#     cross-check, its manifest, and the documents that name it.
# ---------------------------------------------------------------------------
local({
  rd <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  N <- "05_output/20260929/pooled-CPUE-2024-25-220449"; G <- "05_output/20260930/gear-type-CPUE-model-2024-25"
  if (!dir.exists(N) || !dir.exists(G)) chk("B61 render: the committed stage-A folders are present (the documents cite them)", FALSE,
                                           "missing; they are committed, so their absence is a checkout problem") else {
    row <- function(dir) { d <- utils::read.csv(file.path(dir, "port_total_Dungeness_Kept.csv"), stringsAsFactors = FALSE); d[d[[2]] == "Expected_Catch", ] }
    e <- row(N); g <- row(G)
    chk("B61 render: pooled 99,822 [82,090, 124,718], PE 88,758, gear 98,382 [81,145, 122,466], from the folders' own CSVs",
        identical(as.numeric(c(e$BSS_median, e$BSS_lo95, e$BSS_hi95, e$PE)), c(99822, 82090, 124718, 88758)) &&
        identical(as.numeric(c(g$BSS_median, g$BSS_lo95, g$BSS_hi95)), c(98382, 81145, 122466)))
    cc <- utils::read.csv("05_output/20260929/cross_check_20260929_220448.csv", stringsAsFactors = FALSE)
    chk("B61 render: the cross-check names these two folders and reads -1.44%, PASS",
        identical(cc$verdict, "PASS") && isTRUE(abs(cc$gear_minus_pooled_pct - (-1.44)) < 1e-9) &&
        identical(cc$pooled_folder, basename(N)) && identical(cc$gear_folder, basename(G)))
    cv <- utils::read.csv(file.path(N, "convergence_report.csv"), stringsAsFactors = FALSE)
    cg <- utils::read.csv(file.path(G, "convergence_report.csv"), stringsAsFactors = FALSE)
    chk("B61 render: every fit on both tracks passes the gate; the shore all-gear fit is 0.49% divergent (A31: 4.07% before)",
        all(cv$pass_convergence) && all(cg$pass_convergence) &&
        isTRUE(abs(cv$divergence_fraction[cv$fit == "shore_all_gear_Dungeness_Kept"] - 0.0049) < 1e-9))
    mf <- readLines("05_output/20260929/run_manifest_20260929_220448.txt", warn = FALSE)
    chk("B61 render: the manifest records code 2523e8e on a clean tree, both models, R 4.2.2 / rstan 2.32.7",
        any(grepl("^git sha\\s*: 2523e8e", mf)) && any(grepl("^git tree\\s*: clean", mf)) &&
        any(grepl("pooled-CPUE-2024-25-220449", mf, fixed = TRUE)) && any(grepl("rstan_2.32.7", mf, fixed = TRUE)))
    st <- function(d, k) { l <- grep(paste0("^", k, ":"), readLines(file.path(d, "AB_STAGE.txt"), warn = FALSE), value = TRUE); trimws(sub("^[^:]*:", "", l[1])) }
    chk("B61 render: both folders carry the batch's stage-A stamp with one digest",
        identical(st(N, "stage"), "A") && identical(st(G, "stage"), "A") && identical(st(N, "digest"), st(G, "digest")))
    B57g <- "05_output/20260929/gear-type-CPUE-model-2024-25"
    same <- function(f) identical(unname(tools::md5sum(file.path(G, f))), unname(tools::md5sum(file.path(B57g, f))))
    chk("B61 render (B38 in the field): the four gear fits' draws are byte-identical to the B57 gear render (same gear code, seed and radius)",
        dir.exists(B57g) && all(vapply(sprintf("bss_draws_summed_%s.csv", c("shore_ring_net_only_Dungeness_Kept", "shore_all_gear_Dungeness_Kept",
                                                                         "private_boat_ring_net_only_Dungeness_Kept", "private_boat_all_gear_Dungeness_Kept")),
                                       same, logical(1))))
  }
  ps <- rd("07_documentation/development_notes/PIPELINE_STATUS.md"); cr <- rd("07_documentation/development_notes/CHANGE_REGISTER.md")
  vc <- rd("07_documentation/development_notes/VALIDATION_CAMPAIGN.md"); cl <- rd("07_documentation/CLAUDE.md"); rc <- rd("run_config.R")
  chk("B61 docs: the box, the register, CLAUDE.md and run_config.R's header name the new run with its total and interval",
      # 2026-10-04 (B65): superseded; each document now names it as the run the A33/A34 render replaced
      grepl("**`05_output/20260929/pooled-CPUE-2024-25-220449`, port total 99,822 [82,090, 124,718]**", ps, fixed = TRUE) &&
      grepl("It replaced `05_output/20260929/pooled-CPUE-2024-25-220449`, port total **99,822 [82,090, 124,718]**", cr, fixed = TRUE) &&
      grepl("`05_output/20260929/pooled-CPUE-2024-25-220449`, port **99,822 [82,090, 124,718]**", cl, fixed = TRUE) &&
      grepl("#   port total 99,822  [82,090, 124,718]", rc, fixed = TRUE))
  chk("B61 docs: the register carries B61, the campaign 1z.9 with its anchor, and the box the cross-check and R2",
      grepl("| B61 |", cr, fixed = TRUE) && grepl("### 1z.9 The overnight batch read", vc, fixed = TRUE) && grepl("`d30161a` (the overnight batch's renders, 1z.9)", vc, fixed = TRUE) &&
      # 2026-10-02 (A32): the box's cross-check moved to stage D6; stage A's gear render stays named as the one it replaced
      grepl("98,382, -1.44%, the gear track before A32", gsub("[ \n>]+", " ", ps), fixed = TRUE) && grepl("The radius is a precaution, not a crutch", ps, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 95. B62 / A32 (2026-10-02): THE BATCH FINISHED. D3 and D6 adopted TOGETHER (D3 alone trapped a
#     chain at the shipped seed), D29 kept monthly by rule with its span reported, the batch
#     runner and the D3/D6 ladder guarded, two output-identical fixes held until the batch ended.
# ---------------------------------------------------------------------------
local({
  flat <- function(x) gsub("[ \n>]+", " ", paste(x, collapse = "\n"))
  rd <- function(p) readLines(p, warn = FALSE)
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  D6  <- "05_output/20260930/gear-type-CPUE-model-2024-25-AB-D6-gear-matched-zi"
  AP  <- "05_output/20260929/pooled-CPUE-2024-25-220449"
  chk("A32: run_config.R ships D3's per-population gear periods and the zero-inflated catch on both tracks",
      identical(rc$gear_period_bss, list(shore = list(all_gear = "weekly", pot_closure = "biweekly"),
                                         private_boat = list(all_gear = "month", pot_closure = "biweekly"))) &&
      identical(as.character(rc$catch_zi_tracks), c("pooled", "gear_resolved")))
  if (dir.exists(D6)) {
    rp <- rd(file.path(D6, "run_parameters.txt"))
    chk("A32: the stage D6 render was fitted at exactly those values, the shipped seed and the shipped radius",
        any(grepl('^ \\$ catch_zi_tracks +: chr \\[1:2\\] "pooled" "gear_resolved"$', rp)) &&
        any(grepl("^ \\$ bss_seed +: num 20260619$", rp)) && any(grepl("^ \\$ bss_init_r +: num 0.5$", rp)) &&
        { i <- grep("^ \\$ gear_period_bss", rp); length(i) == 1 &&
          identical(trimws(rp[i + 1:6]), c('..$ shore       :List of 2', '.. ..$ all_gear   : chr "weekly"', '.. ..$ pot_closure: chr "biweekly"',
                                           '..$ private_boat:List of 2', '.. ..$ all_gear   : chr "month"', '.. ..$ pot_closure: chr "biweekly"')) })
    st <- rd(file.path(D6, "AB_STAGE.txt"))
    chk("A32: stage D6 was rendered from 2523e8e, the commit stage A (the pooled authoritative render) was rendered from",
        any(st == "git: 2523e8e") && any(st == "stage: D6") &&
        any(rd(file.path(AP, "AB_STAGE.txt")) == "git: 2523e8e"))
    cr <- utils::read.csv(file.path(D6, "convergence_report.csv"))
    chk("A32: every fit of the adopted gear configuration at the shipped seed passes its gate (no stuck chain)",
        nrow(cr) == 4 && all(cr$method_selected == "BSS") && max(cr$divergence_fraction) < 0.005)
    rdp <- function(d) { pt <- utils::read.csv(file.path(d, "port_total_Dungeness_Kept.csv")); pt[pt$Estimate == "Expected_Catch", ] }
    pp <- rdp(AP); gg <- rdp(D6); cc <- utils::read.csv("05_output/authoritative_batch_2026-09-29_cross_check_adopted.csv")
    chk("A32: the adopted cross-check file is the two folders' own medians: 99,822 against 99,294, -0.53%, PASS",
        pp$BSS_median == 99822 && gg$BSS_median == 99294 && gg$BSS_lo95 == 82154 && gg$BSS_hi95 == 123380 &&
        cc$pooled_median == pp$BSS_median && cc$gear_median == gg$BSS_median &&
        isTRUE(all.equal(cc$gear_minus_pooled_pct, round(100 * (gg$BSS_median - pp$BSS_median) / pp$BSS_median, 2))) &&
        cc$gear_minus_pooled_pct == -0.53 && cc$verdict == "PASS" && cc$gear_folder == basename(D6))
    d3 <- utils::read.csv("05_output/20260930/gear-type-CPUE-model-2024-25-AB-D3-gear-matched/sampler_diagnostics_shore_all_gear_Dungeness_Kept.csv")
    chk("A32's reason, from the file: D3 ALONE at the shipped seed trapped chain 1 of the gear shore all-gear fit (1,546 of 2,000 divergent)",
        d3$divergent[d3$chain == 1] == 1546 && all(d3$divergent[d3$chain != 1] <= 1))
  } else skp("A32: the stage D6 and D3 folders", "05_output not present")
  vv <- utils::read.csv("05_output/authoritative_batch_2026-09-29_recommendation.csv")
  chk("B62: the batch's own recommendations are D3 ADOPT, D6 ADOPT and D29 KEEP monthly",
      grepl("^ADOPT", vv$recommendation[vv$item == "D3"]) && grepl("^ADOPT", vv$recommendation[vv$item == "D6"]) &&
      grepl("^KEEP monthly", vv$recommendation[vv$item == "D29"]))
  d29 <- utils::read.csv("05_output/authoritative_batch_2026-09-29_d29_rungs.csv")
  po <- d29[d29$track == "pooled", ]
  chk("B62 (D29): the pooled rungs read monthly 47,192, biweekly 50,888, weekly 53,546 (all adequate, NO EVIDENCE), daily out",
      po$catch_median[po$rung == "monthly"] == 47192 && po$catch_median[po$rung == "biweekly"] == 50888 &&
      po$catch_median[po$rung == "weekly"] == 53546 && all(po$adequate[po$rung != "daily"]) &&
      all(po$standing[po$rung %in% c("weekly", "biweekly")] == "NO EVIDENCE") && po$standing[po$rung == "daily"] == "out: inadequate")
  # the two fixes held until the batch ended
  pdc <- rd("03_R_functions/prep_days_crab.R")
  chk("B62: prep_days_crab() sets period with an if/else, not a scalar-condition case_when() (dplyr 1.2.0)",
      any(grepl('period = if (isTRUE(period_pe == "month")) year * 100 + month else iso_year * 100 + week,', pdc, fixed = TRUE)) &&
      !any(grepl('period_pe == "month" ~', pdc, fixed = TRUE)))
  ex <- new.env(); for (f in c("classify_day_type.R", "prep_days_crab.R")) if (file.exists(file.path("03_R_functions", f))) sys.source(file.path("03_R_functions", f), envir = ex)
  if (exists("prep_days_crab", envir = ex) && exists("bss_weekday", envir = ex) && requireNamespace("dplyr", quietly = TRUE)) {
    ex$tibble <- tibble::tibble; ex$case_when <- dplyr::case_when; ex$bss_assign_day_length <- function(days, ...) days
    pr <- list(days_wkend = c("Saturday", "Sunday"), crabbing_holiday_dates = as.Date("2025-01-01"), sections = 1)
    dm <- tryCatch(ex$prep_days_crab("2024-12-28", "2025-01-06", modifyList(pr, list(period_pe = "month"))), error = function(e) NULL)
    dw <- tryCatch(ex$prep_days_crab("2024-12-28", "2025-01-06", modifyList(pr, list(period_pe = "week"))), error = function(e) NULL)
    chk("B62: prep_days_crab()'s period is year-month under period_pe month and ISO year-week otherwise (across a New Year)",
        !is.null(dm) && !is.null(dw) && identical(as.numeric(dm$period), dm$year * 100 + dm$month) &&
        identical(as.numeric(dw$period), dw$iso_year * 100 + dw$week) && dw$period[dw$event_date == as.Date("2025-01-01")] == 202501 &&
        dw$period[dw$event_date == as.Date("2024-12-30")] == 202501)
  } else skp("B62: prep_days_crab() period, functionally", "classify_day_type.R or dplyr not available")
  bp <- new.env(); sys.source("03_R_functions/bss_packages.R", envir = bp)
  if (!isNamespaceLoaded("renv")) {
    old <- Sys.getenv("RENV_PROJECT", unset = NA); Sys.setenv(RENV_PROJECT = tempdir())
    msg <- tryCatch({ bp$bss_load_packages(pkgs = "notapkg.b62", attach = character(0)); "no error" }, error = function(e) conditionMessage(e))
    if (is.na(old)) Sys.unsetenv("RENV_PROJECT") else Sys.setenv(RENV_PROJECT = old)
    chk("B62: bss_load_packages() in a process that inherited RENV_PROJECT without renv loaded STOPS (no restore, no install)",
        grepl("inherited RENV_PROJECT without renv active in it", msg, fixed = TRUE) && grepl("installs nothing", msg, fixed = TRUE))
  } else skp("B62: bss_load_packages() in an inherited-RENV_PROJECT process", "renv is loaded in this session, so the case cannot be staged here")
  ab <- rd("06_diagnostics/run_authoritative_batch_2026-09-29.R")
  chk("B62: the batch desk check reports whether renv is LOADED in the stage process, not whether RENV_PROJECT is set",
      any(grepl("isNamespaceLoaded('renv')", ab, fixed = TRUE)) && any(grepl("RENV_PROJECT inherited from this session", ab, fixed = TRUE)) &&
      !any(grepl('if (!is.na(rv) && nzchar(rv)) sprintf("ACTIVE (%s)", rv) else "off"', ab, fixed = TRUE)))
  for (rn in c("run_authoritative_batch_2026-09-29", "run_gear_ar_zi_2026-09-13")) {
    t <- rd(file.path("06_diagnostics", paste0(rn, ".R")))
    g <- grep("bss_superseded_runner(", t, fixed = TRUE)[1]; d <- grep("^DRY_RUN <- TRUE", t)[1]
    chk(sprintf("B62: %s is guarded BEFORE its first setting (the guard runs before DRY_RUN is read)", rn),
        !is.na(g) && !is.na(d) && g < d && any(grepl('settled_by = paste(', t, fixed = TRUE)))
  }
  chk("B62: the superseded-runner message no longer names the D3/D6 ladder as live",
      !grepl("run_gear_ar_zi_2026-09-13.R", paste(rd("03_R_functions/bss_superseded_runner.R"), collapse = " "), fixed = TRUE))
  # documents
  ps <- flat(rd("07_documentation/development_notes/PIPELINE_STATUS.md")); cr <- paste(rd("07_documentation/development_notes/CHANGE_REGISTER.md"), collapse = "\n")
  vc <- paste(rd("07_documentation/development_notes/VALIDATION_CAMPAIGN.md"), collapse = "\n"); cl <- flat(rd("07_documentation/CLAUDE.md"))
  chk("B62 docs: the box names the stage D6 folder and -0.53%, and quotes the D29 span to 106,380 as outside the interval",
      grepl("The cross-check passes: 99,294 [82,154, 123,380], -0.53% against the pooled total", ps, fixed = TRUE) &&
      grepl("gear-type-CPUE-model-2024-25-AB-D6-gear-matched-zi", ps, fixed = TRUE) &&
      grepl("106,380 [86,780, 131,940] at weekly (+6.6%)", ps, fixed = TRUE) && grepl("it is NOT inside the interval above", ps, fixed = TRUE))
  chk("B62 docs: the register carries A32 (five cells), B62, four 2026-10-02 defect rows (B63 added the fourth), and D3, D6, D29 CLOSED",
      grepl("| A32 |", cr, fixed = TRUE) && grepl("| B62 |", cr, fixed = TRUE) && length(gregexpr("\n| 2026-10-02 |", cr, fixed = TRUE)[[1]]) == 4 &&
      grepl("| **CLOSED 2026-10-02: ADOPTED together with D6 (A32, B62).**", cr, fixed = TRUE) &&
      grepl("| **CLOSED 2026-10-02: ADOPTED together with D3 (A32, B62)", cr, fixed = TRUE) &&
      grepl("| **CLOSED 2026-10-02 by rule D29-4(b)", cr, fixed = TRUE))
  chk("B62 docs: the campaign has 1z.10 with its anchor, and CLAUDE.md names the new cross-check",
      grepl("### 1z.10 The batch finished", vc, fixed = TRUE) && grepl("`0c383f2` (its remaining stages)", vc, fixed = TRUE) &&
      grepl("99,294 (-0.53%, PASS)", cl, fixed = TRUE) && grepl("catch_zi_tracks` ships `c(\"pooled\", \"gear_resolved\")", cl, fixed = TRUE))
  chk("B62 docs: no live document still says catch_zi_tracks ships \"pooled\" as the current state",
      !any(vapply(c("07_documentation/CLAUDE.md", "07_documentation/NEW_SEASON_GUIDE.md", "01_BSS_models/README.md", "02_stan_models/README.md"),
                  function(f) grepl('ships `"pooled"`', paste(rd(f), collapse = " "), fixed = TRUE), logical(1))))
})

# ---------------------------------------------------------------------------
# 96. B63 (2026-10-02): a refused shared BOAT turnover reaches the report, and D18 says what a
#     2025-26 run actually does (the boat all-gear fit has no OSP day, so it is refused).
# ---------------------------------------------------------------------------
local({
  rd <- function(p) paste(readLines(p, warn = FALSE), collapse = "\n")
  es <- new.env(); sys.source("03_R_functions/bss_run_warnings.R", envir = es); sys.source("03_R_functions/bss_effort_spec.R", envir = es)
  es$`%||%` <- function(a, b) if (is.null(a)) b else a
  f <- es$bss_shared_tau_data
  if (is.function(f) && exists("bss_warn_reset", envir = es)) {
    es$bss_warn_reset()
    off <- tryCatch(suppressWarnings(utils::capture.output(r <- f(list(L_unit = "turnover (gear-deployments per gear slot)"), rep(3.5, 4), rep(0.5, 4),
                    list(shared_tau = TRUE, shared_tau_min_obs = 15, tau_boat_prior_mu = 3.5), population_name = "private_boat", n_informed = 0L))),
                    error = function(e) conditionMessage(e))
    lg <- tryCatch(es$bss_warn_log(), error = function(e) NULL)
    chk("B63: a refused shared turnover for a BOAT fit is written to the run warnings (source turnover)",
        !is.null(lg) && any(lg$source == "turnover" & grepl("Shared boat turnover REFUSED", lg$message, fixed = TRUE)),
        sprintf("(%s)", paste(utils::head(off, 2), collapse = " | ")))
  } else {
    t <- rd("03_R_functions/bss_effort_spec.R")
    chk("B63: a refused shared turnover for a BOAT fit calls bss_warn(\"turnover\", ...) (static check)",
        grepl('bss_warn("turnover", sprintf(paste0("Shared boat turnover REFUSED', t, fixed = TRUE) &&
        grepl('identical(population_name, "private_boat")', t, fixed = TRUE))
  }
  cr <- rd("07_documentation/development_notes/CHANGE_REGISTER.md"); ng <- rd("07_documentation/NEW_SEASON_GUIDE.md")
  ps <- gsub("[ \n>]+", " ", rd("07_documentation/development_notes/PIPELINE_STATUS.md"))
  chk("B63: D18, the guide and the box say the 2025-26 boat all-gear fit has no OSP day and its shared turnover is REFUSED, not fired",
      grepl("all 23 OSP days fall in the pot closure", cr, fixed = TRUE) && grepl("| B63 |", cr, fixed = TRUE) &&
      grepl("the floor is applied PER FIT", ng, fixed = TRUE) && !grepl("so the shared turnover and its calibration would FIRE", ng, fixed = TRUE) &&
      grepl("boat ALL-GEAR fit has no OSP day at all", ps, fixed = TRUE) && !grepl("shared turnover calibration would fire on autumn-only overlap (D18)", ps, fixed = TRUE))
  chk("B63: the gear driver no longer says its zero-inflated catch ships off",
      !grepl('ships OFF (catch_zi_tracks =', rd("01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd"), fixed = TRUE))
})

# ---------- D40 (2026-10-04): the crabbing share's day-level volume term ----------
local({
  rdf <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  # (1) the observed daily volume: OSP's total first, else trailers x turnover centre, else NA
  dd <- mkdays("2025-06-01", 6)
  effd <- tibble(event_date = dd$event_date[c(1, 1, 2, 4)], count_quantity = c(10, 20, 8, 0))
  ospm <- tibble(event_date = dd$event_date[c(2, 3)], count_quantity = c(120, 40))
  v <- crab_fraction_day_volume(dd, effd, ospm, L_data = rep(3, 6))
  chk("D40 volume: OSP's total wins, else mean trailer x turnover centre, else NA (a zero count is not a volume)",
      isTRUE(all.equal(v[1:3], c(45, 120, 40))) && all(is.na(v[4:6])) &&
      identical(attr(v, "source")[1:3], c("trailer x turnover", "OSP total", "OSP total")))
  # (2) the Stan fields: inert unless asked, with lengths that always match their streams
  Pv <- modifyList(P, list(crab_fraction_strata = "month", crab_fraction_dynamic = TRUE, use_osp_crab_lower = TRUE))
  Pv$crab_fraction_rows <- tibble(event_date = days289$event_date[c(1, 2, 40, 41, 100)], boats_total = c(3, 2, 5, 1, 4), boats_crabbing = c(3, 1, 4, 0, 1))
  Pv$osp_crab_rows <- dplyr::mutate(osp_rows, osp_boat_total = rep(c(20, 80), 100))
  vol <- rep(NA_real_, 289); vol[c(1, 2, 40, 100)] <- c(10, 40, 5, 20); vol[1:200][is.na(vol[1:200])] <- 30; vol[41] <- NA
  cf_off <- crab_fraction_stan_data(FALSE, days289, Pv, quiet = TRUE, day_volume = vol)
  chk("D40 inert by default: f_volume 0, zero x of the right lengths, the centres 0",
      cf_off$f_volume == 0L && length(cf_off$cfi_x) == cf_off$CFI_n && all(cf_off$cfi_x == 0) &&
      length(cf_off$osp_f_x) == cf_off$OSPF_n && all(cf_off$osp_f_x == 0) && all(cf_off$fvol_centre == 0))
  Pv1 <- Pv; Pv1$crab_fraction_volume <- TRUE
  cf_on <- crab_fraction_stan_data(FALSE, days289, Pv1, quiet = TRUE, day_volume = vol)
  st <- cf_on$f_stratum; lv <- log(vol)
  cen <- vapply(seq_len(cf_on$n_f_strata), function(k) mean(lv[st == k], na.rm = TRUE), numeric(1))
  chk("D40 on: f_volume 1, each stratum's centre is the mean log volume of its days that have one",
      cf_on$f_volume == 1L && isTRUE(all.equal(as.numeric(cf_on$fvol_centre)[is.finite(cen)], cen[is.finite(cen)])))
  chk("D40 on: a contact day's x is its log volume minus its stratum's centre",
      isTRUE(all.equal(cf_on$cfi_x[1], log(10) - cen[st[1]])) && isTRUE(all.equal(cf_on$cfi_x[3], log(5) - cen[st[40]])))
  chk("D40 on: a contact day with no observed volume is read at the centre (x = 0) and counted",
      cf_on$cfi_x[4] == 0 && identical(attr(cf_on, "f_volume_missing_contacts"), 1L))
  chk("D40 on: an OSP day's x is from OSP's OWN daily total, not the day-volume vector",
      isTRUE(all.equal(cf_on$osp_f_x[1:2], c(log(20), log(80)) - cen[st[1:2]])))
  chk("D40 on: nothing but the five volume fields moves",
      identical(cf_off[setdiff(names(cf_off), c("f_volume", "fvol_centre", "cfi_x", "osp_f_x"))],
                cf_on[setdiff(names(cf_on), c("f_volume", "fvol_centre", "cfi_x", "osp_f_x"))]))
  chk("D40: the term needs the dynamic f (legacy -> off) and a volume vector (none -> off)",
      crab_fraction_stan_data(FALSE, days289, modifyList(Pv1, list(crab_fraction_dynamic = FALSE)), quiet = TRUE, day_volume = vol)$f_volume == 0L &&
      crab_fraction_stan_data(FALSE, days289, Pv1, quiet = TRUE)$f_volume == 0L)
  chk("D40: the shore never carries it", crab_fraction_stan_data(TRUE, days289, Pv1, quiet = TRUE, day_volume = vol)$f_volume == 0L)
  # (3) the Stan programs: same term on both tracks, zero-size when off, one share per day in the totals
  for (sf in c("02_stan_models/crab_bss_pooled.stan", "02_stan_models/crab_bss_gear_resolved.stan")) {
    st_ <- rdf(sf)
    chk(sprintf("D40 %s: declares the five data fields, n_fvol and a ZERO-SIZE beta_fvol", basename(sf)),
        all(vapply(c("int<lower=0,upper=1> f_volume;", "real<lower=0> fvol_beta_prior_sd;", "vector[n_f_strata] fvol_centre;",
                     "vector[CFI_n] cfi_x;", "vector[OSPF_n] osp_f_x;", "int n_fvol = n_f_dyn * f_volume;", "vector[n_fvol] beta_fvol;"),
                   function(x) grepl(x, st_, fixed = TRUE), logical(1))))
    chk(sprintf("D40 %s: the totals use the day's share fd (from the model's own volume), and fall back to the stratum's share bit for bit", basename(sf)),
        grepl("real fd = f_crab[f_stratum[d]];", st_, fixed = TRUE) && grepl("lambda_C_S[s][d,g] * fd * zi_scale;", st_, fixed = TRUE) &&
        grepl("E[s][d,g] = lambda_E_S[s][d,g] * E_scale * L[d] * fd;", st_, fixed = TRUE) &&
        !grepl("f_crab[f_stratum[d]] * zi_scale", st_, fixed = TRUE) && grepl("vd = fmax(vd * L[d], 1e-9);", st_, fixed = TRUE))
    chk(sprintf("D40 %s: both classification streams see the day's share when the term is on", basename(sf)),
        grepl("eta_f[cfi_stratum[i]] + beta_fvol[1] * cfi_x[i]", st_, fixed = TRUE) &&
        grepl("eta_f[osp_f_stratum[i]] + beta_fvol[1] * osp_f_x[i]", st_, fixed = TRUE) &&
        grepl("p_osp = f_osp * (1 - combo_c[osp_f_stratum[i]]);", st_, fixed = TRUE))
  }
  chk("D40: the pooled volume is the OSP stream's own mean (lambda_E[, G] / R_G_boat x L); the gear track sums its gears as its OSP mean does",
      grepl("vd += lambda_E_S[s2][d, G] / R_G_boat;", rdf("02_stan_models/crab_bss_pooled.stan"), fixed = TRUE) &&
      grepl("vd += sum(lambda_E_S[s2][d, ]) / R_G_boat;", rdf("02_stan_models/crab_bss_gear_resolved.stan"), fixed = TRUE))
  # (4) the preps build the volume and forward the fields
  for (pf in c("03_R_functions/prep_bss_crab_pooled.R", "03_R_functions/prep_bss_crab_gear.R")) {
    pt <- rdf(pf)
    chk(sprintf("D40 %s: builds the day volume for boat fits and forwards the five fields", basename(pf)),
        grepl("crab_fraction_day_volume(days, eff_d, osp_match,", pt, fixed = TRUE) &&
        grepl("day_volume = f_day_volume)", pt, fixed = TRUE) &&
        all(vapply(c("f_volume               = cf_data$f_volume", "fvol_beta_prior_sd     = cf_data$fvol_beta_prior_sd",
                     "fvol_centre            = cf_data$fvol_centre", "cfi_x                  = cf_data$cfi_x",
                     "osp_f_x                = cf_data$osp_f_x"), function(x) grepl(x, pt, fixed = TRUE), logical(1))))
  }
  e <- new.env(); sys.source("run_config.R", envir = e); rc <- e$run_config
  cr_ <- rdf("07_documentation/development_notes/CHANGE_REGISTER.md"); rcs <- rdf("run_config.R")
  chk("D40/A34 docs: the register has the A34 row, D40 reads FIXED IN CODE, and the rule's winter clause is retired for the beta_fvol clause in run_config and the register",
      grepl("| A34 |", cr_, fixed = TRUE) && grepl("FIXED IN CODE 2026-10-04 (A34", cr_, fixed = TRUE) &&
      grepl("beta_fvol's 95% interval excludes", rcs, fixed = TRUE) && grepl("was RETIRED with D40", rcs, fixed = TRUE) &&
      grepl("clause (4) is RETIRED", cr_, fixed = TRUE))
  chk("D40 shipped: crab_fraction_volume = TRUE with a N(0, 1) slope prior, alongside the dynamic f and the OSP stream",
      isTRUE(rc$crab_fraction_volume) && identical(rc$crab_fraction_volume_beta_sd, 1) && isTRUE(rc$crab_fraction_dynamic) && isTRUE(rc$use_osp_crab_lower))
})

# ---------- B65 (2026-10-04): the A33/A34 render is the authoritative run ----------
local({
  rdf <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")
  N <- "05_output/20261003/pooled-CPUE-2024-25"; G <- "05_output/20261004/gear-type-CPUE-model-2024-25"
  if (dir.exists(N) && dir.exists(G)) {
    e <- utils::read.csv(file.path(N, "port_total_Dungeness_Kept.csv")); e <- e[e$Estimate == "Expected_Catch", ]
    g <- utils::read.csv(file.path(G, "port_total_Dungeness_Kept.csv")); g <- g[g$Estimate == "Expected_Catch", ]
    chk("B65 render: pooled 87,932 [75,193, 105,271], PE 88,758, gear 87,903 [75,148, 104,877], from the folders' own CSVs",
        identical(as.numeric(c(e$BSS_median, e$BSS_lo95, e$BSS_hi95, e$PE)), c(87932, 75193, 105271, 88758)) &&
        identical(as.numeric(c(g$BSS_median, g$BSS_lo95, g$BSS_hi95)), c(87903, 75148, 104877)))
    cv <- utils::read.csv(file.path(N, "convergence_report.csv")); cg <- utils::read.csv(file.path(G, "convergence_report.csv"))
    chk("B65 rule (1): every fit on both tracks passes the gate", nrow(cv) == 4 && nrow(cg) == 4 && all(cv$pass_convergence) && all(cg$pass_convergence))
    sp <- function(d) { f <- list.files(d, "^sampler_diagnostics_.*csv$", full.names = TRUE)
                        all(vapply(f, function(x) { z <- utils::read.csv(x); !any(z$ebfmi_low_flag) && max(z$mean_stepsize) / min(z$mean_stepsize) < 3 }, logical(1))) }
    chk("B65 rule (2): no stuck chain (no E-BFMI flag; per-chain step sizes within a factor of 3 in every fit)", sp(N) && sp(G))
    xc <- utils::read.csv(list.files("05_output/20261003", "^cross_check_.*csv$", full.names = TRUE)[1])
    chk("B65 rule (3): the cross-check is -0.03%, PASS", identical(xc$verdict, "PASS") && abs(xc$gear_minus_pooled_pct) < 2 && isTRUE(all.equal(xc$gear_minus_pooled_pct, -0.03)))
    bq <- function(d) { b <- utils::read.csv(file.path(d, "bss_summary_private_boat_all_gear_Dungeness_Kept.csv")); b[b[[1]] == "beta_fvol_out[1]", ] }
    chk("B65 rule (4): beta_fvol's 95% interval excludes zero on both tracks' boat all-gear fits",
        bq(N)$X97.5. < 0 && bq(G)$X97.5. < 0)
    A <- "05_output/20260929/pooled-CPUE-2024-25-220449"
    chk("B65: the shore fits are byte-identical to the run replaced (A33 and A34 touch only the boat)",
        all(vapply(c("bss_summary_shore_all_gear_Dungeness_Kept.csv", "bss_summary_shore_ring_net_only_Dungeness_Kept.csv"),
                   function(f) identical(unname(tools::md5sum(file.path(A, f))), unname(tools::md5sum(file.path(N, f)))), logical(1))))
    mf <- readLines(list.files("05_output/20261003", "^run_manifest_.*txt$", full.names = TRUE)[1], warn = FALSE)
    chk("B65: the render ran on 4490fbe with a clean tree and both new keys on",
        any(grepl("git sha     : 4490fbe", mf, fixed = TRUE)) && any(grepl("git tree    : clean", mf, fixed = TRUE)) &&
        any(grepl("use_osp_crab_lower .*TRUE", mf)) && any(grepl("crab_fraction_volume .*TRUE", mf)))
  }
  ps <- rdf("07_documentation/development_notes/PIPELINE_STATUS.md"); cr <- rdf("07_documentation/development_notes/CHANGE_REGISTER.md")
  vc <- rdf("07_documentation/development_notes/VALIDATION_CAMPAIGN.md"); cl <- rdf("07_documentation/CLAUDE.md"); rc <- rdf("run_config.R")
  chk("B65 docs: the box, the register, CLAUDE.md and run_config.R name the new run with its total and interval, and D29's span as stale",
      grepl("**`05_output/20261003/pooled-CPUE-2024-25`, port total 87,932 [75,193, 105,271]**", ps, fixed = TRUE) &&
      grepl("**Authoritative run:** `05_output/20261003/pooled-CPUE-2024-25`, port total **87,932 [75,193, 105,271]**", cr, fixed = TRUE) &&
      grepl("port **87,932 [75,193, 105,271]**", cl, fixed = TRUE) && grepl("#   port total 87,932  [75,193, 105,271]", rc, fixed = TRUE) &&
      grepl("D29's span is stale", ps, fixed = TRUE) && grepl("STALE SINCE B65", cr, fixed = TRUE))
  chk("B65 docs: the register carries B65, A33 and A34 read RENDERED, and the campaign 1z.11 with its anchor",
      grepl("| B65 |", cr, fixed = TRUE) && grepl("| **ADOPTED; RENDERED 2026-10-04 (B65), with A33** |", cr, fixed = TRUE) &&
      grepl("**ADOPTED; RENDERED 2026-10-04 (B65), with A34: the rule met on every clause**", cr, fixed = TRUE) &&
      grepl("### 1z.11 The OSP crab-only count and the day-level crabbing share, rendered", vc, fixed = TRUE) && grepl("`c2cb6af` (the A33/A34 render, 1z.11)", vc, fixed = TRUE))
})

cat(sprintf("\n==== %d passed, %d failed, %d skipped ====\n", ok, bad, skipped))
if (bad > 0) quit(status = 1)
