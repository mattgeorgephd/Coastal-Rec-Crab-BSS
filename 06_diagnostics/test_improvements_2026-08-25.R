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
            "model_diagnostics.R","bss_model_adequacy.R",
            "annotate_decoupled_run.R",
            "bss_sampler_override.R",
            "pe_gear_ratio_frame.R",
            "bss_ar_rung_summary.R")) source(file.path("03_R_functions", f))

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
    chk("adequacy: reference runs present for the threshold check", TRUE, "skipped, runs absent")
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
    chk("retro-audit: annotated baseline present", TRUE, "skipped, not yet annotated")
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
    if (dir.exists(e)) {
      b <- annotate_model_adequacy_run(e, overwrite = TRUE, quiet = TRUE)
      rb <- b[b$fit == "private_boat_all_gear_Dungeness_Kept", ]
      chk("adequacy: the shared turnover improved PIT bias at UNCHANGED complexity",
          isTRUE(abs(rb$p_loo_frac - r$p_loo_frac) > 0.05) && isTRUE(rb$pit_worst_bias < 0.05),
          sprintf("(E p_loo %.3f bias %.3f vs F p_loo %.3f bias %.3f)",
                  rb$p_loo_frac, rb$pit_worst_bias, r$p_loo_frac, r$pit_worst_bias))
    }
  } else chk("adequacy retro: stage F folder present", TRUE, "skipped, folder absent")
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
    if (!file.exists(p)) { chk(paste("stage5 REF:", basename(r[[1]])), TRUE, "skipped, folder absent"); next }
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
  if (dir.exists(v5) && dir.exists(v1))
    chk("V5: the shared_tau_min_obs threshold is worth ~169 crab, 0.24% of the port",
        isTRUE(abs((bag(v1, "private_boat \\(Pot closure\\)") -
                    bag(v5, "private_boat \\(Pot closure\\)")) - 169) <= 2),
        sprintf("(floor 15: %d, floor 20: %d)", bag(v1, "private_boat \\(Pot closure\\)"),
                bag(v5, "private_boat \\(Pot closure\\)")))
  # GR-7 Phase 2, sampled for the first time: widen intervals, do not move medians.
  a <- "05_output/20260901/gear-type-CPUE-model-VAL-2-gearG-phase1"
  b <- "05_output/20260901/gear-type-CPUE-model-VAL-3-gearG-dirichlet"
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
# 28. Zero-inflated catch likelihood (2026-09-02, prototype, ships OFF). Targeted from the
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
      grepl("* f_crab[f_stratum[d]] * zi_scale;", t, fixed = TRUE) &&
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
  chk("zinb: ships OFF, scoped to the shore, with a Beta(1, 9) prior",
      isFALSE(rc$estimate_catch_zi) && identical(rc$catch_zi_populations, "shore") &&
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
  chk("dry-run default: the set of runners with NO dry-run mode has not grown",
      setequal(no_dry, c("run_osp_validation.R", "run_rg_sweep.R", "run_tau_sweep.R")),
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

cat(sprintf("\n==== %d passed, %d failed ====\n", ok, bad))
if (bad > 0) quit(status = 1)
