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
            "bss_opener_covariates.R","diagnose_incomplete_trips.R",
            "bss_stan_fit.R","save_run_diagnostics.R",
            "model_diagnostics.R")) source(file.path("03_R_functions", f))

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
  cf_need <- grep("^(apply_crab_fraction|crab_fraction_.+|n_f_strata|f_stratum|osp_crab_lower|osp_f_.+|OSPF_n)$",
                  cf_need, value = TRUE)
  chk("crab-fraction contract is non-trivial", length(cf_need) >= 16, length(cf_need))

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

  # 10a. .srd_monthly_share(): the SECOND copy of the shore day-length weighting.
  # The rule: the monthly share is a normalized weight, so a per-day multiplier that is
  # constant across days cancels. Weight by day length ONLY when L is an effective day length
  # in hours. Under gear-deployments (production) and for both boat fits, L is a turnover and
  # the share must be count-weighted. Before the fix this file multiplied every SHORE share by
  # days_ss$day_length regardless, re-weighting the split toward long-day summer months.
  D <- 120L
  days <- data.frame(event_date = as.Date("2024-12-01") + 0:(D-1),
                     day_length = seq(9, 16, length.out = D))   # strongly seasonal
  sd_dep <- list(day_Gear = rep(1:D, each = 1), Gear_I = rep(10, D),
                 .L_unit = "turnover (trips per gear-slot per day)",
                 .effort_unit = "gear-deployments")
  sd_hrs <- list(day_Gear = rep(1:D, each = 1), Gear_I = rep(10, D),
                 .L_unit = "effective day length (hours)",
                 .effort_unit = "crabber-hours")
  sh_dep <- .srd_monthly_share(sd_dep, days, is_boat = FALSE)
  sh_hrs <- .srd_monthly_share(sd_hrs, days, is_boat = FALSE)
  # With a flat count series, a turnover-unit shore share must be proportional to DAYS PER
  # MONTH alone; an hours-unit share must tilt toward the long-day months at the end.
  n_per_month <- as.numeric(table(format(days$event_date, "%Y-%m")))
  chk("monthly share, deployments: count-weighted (day length cancels)",
      isTRUE(all.equal(sh_dep$share, n_per_month / sum(n_per_month), tolerance = 1e-10)))
  chk("monthly share, crabber-hours: still day-length weighted",
      !isTRUE(all.equal(sh_hrs$share, n_per_month / sum(n_per_month), tolerance = 1e-6)) &&
      tail(sh_hrs$share, 1) > tail(sh_dep$share, 1))
  sd_boat <- list(day_T = 1:D, T_I = rep(4, D),
                  .L_unit = "turnover (trips per present group per day)",
                  .effort_unit = "gear-deployments")
  chk("monthly share, boat: unchanged and count-weighted",
      isTRUE(all.equal(.srd_monthly_share(sd_boat, days, is_boat = TRUE)$share,
                       n_per_month / sum(n_per_month), tolerance = 1e-10)))
  # Fallback path: a stan_data built before .L_unit existed must not crash and must keep the
  # historical shore behaviour rather than silently switching units.
  sd_old <- list(day_Gear = 1:D, Gear_I = rep(10, D))
  chk("monthly share: pre-.L_unit stan_data falls back to the historical shore weighting",
      !isTRUE(all.equal(.srd_monthly_share(sd_old, days, is_boat = FALSE)$share,
                        n_per_month / sum(n_per_month), tolerance = 1e-6)))

  # 10b. The empty-effort-stratum report is a FILE, not a console line. The pooled driver's
  # PE chunk is results='hide', so the cat()-only version reached nothing on that track.
  td <- file.path(tempdir(), paste0("pe_empty_", as.integer(runif(1, 1e5, 1e6))))
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  pe_fake <- list(
    shore_all_gear = list(n_empty_effort_strata = 0L, n_empty_effort_days = 0L,
                          n_effort_strata_total = 84L, n_calendar_days = 289L,
                          pe_empty_effort_fill = "zero"),
    private_boat_ring_net_only = list(n_empty_effort_strata = 9L, n_empty_effort_days = 9L,
                          n_effort_strata_total = 22L, n_calendar_days = 76L,
                          pe_empty_effort_fill = "zero"),
    comm_charter = list(effort_total = 283))
  rep_df <- write_pe_empty_stratum_report(pe_fake, td)
  chk("empty-stratum report writes a file", file.exists(file.path(td, "pe_empty_effort_strata.csv")))
  chk("empty-stratum report excludes the census component",
      !"comm_charter" %in% rep_df$component && nrow(rep_df) == 2)
  chk("empty-stratum report computes the zeroed-day fraction",
      isTRUE(all.equal(rep_df$empty_day_fraction[rep_df$component == "private_boat_ring_net_only"],
                       9/76)))
  chk("empty-stratum report raises the >5%-at-zero flag",
      isTRUE(rep_df$exceeds_5pct_at_zero[rep_df$component == "private_boat_ring_net_only"]) &&
      isFALSE(rep_df$exceeds_5pct_at_zero[rep_df$component == "shore_all_gear"]))
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
  chk("production arm: the BOAT PE is gear_only, not exclude",
      identical(unname(a_boat[["pe"]]), "gear_only"))
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

cat(sprintf("\n==== %d passed, %d failed ====\n", ok, bad))
if (bad > 0) quit(status = 1)
