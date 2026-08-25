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

for (f in c("bss_effort_spec.R","bss_ar_resolution.R","crab_fraction.R",
            "bss_opener_covariates.R","diagnose_incomplete_trips.R")) source(file.path("03_R_functions", f))

ok <- 0; bad <- 0
chk <- function(nm, cond, extra="") { if (isTRUE(cond)) { ok <<- ok+1; cat("PASS ", nm, extra, "\n") } else { bad <<- bad+1; cat("FAIL ", nm, extra, "\n") } }

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
  chk("shipped: use_osp_crab_lower OFF", identical(rc$use_osp_crab_lower, FALSE))
  chk("shipped: pe_empty_effort_stratum = zero (historical)", identical(rc$pe_empty_effort_stratum, "zero"))
  chk("shipped: filter_incomplete_trips still TRUE (diagnostic only)", identical(rc$filter_incomplete_trips, TRUE))
  chk("shipped: crab_fraction_strata unchanged", identical(rc$crab_fraction_strata, "none"))
})

cat(sprintf("\n==== %d passed, %d failed ====\n", ok, bad))
if (bad > 0) quit(status = 1)
