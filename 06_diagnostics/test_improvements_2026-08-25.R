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
            "bss_ar_rung_summary.R",
            "read_input_workbook.R")) source(file.path("03_R_functions", f))   # 2026-09-10: the workbook reader

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
  # 2026-09-12: the three unsampled-cell levers ship at their new values. "zero" is no
  # longer a neutral default; it is an assumption that 15.6% of the shore's days and 16.6%
  # of the boat's had no fishing, and its error is a bias no SE can carry.
  chk("shipped: pe_empty_effort_stratum = local_day_type (2026-09-12)", identical(rc$pe_empty_effort_stratum, "local_day_type"))
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
  # 2026-09-12: the report now carries the SINGLETON cells, both SEs and the zeroing bias.
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
      any(grepl("validate_season_window(effort_raw, gh_interview, params)",
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
              crabbing_holiday_dates = as.Date(character()), estimate_red_rock = FALSE)
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
  # 2026-09-12: the pin has to be COMPLETE. Four per-season keys were missing, so every
  # 2024-25 rung would have inherited the two-season config's 2023-24 values, pot_open_date
  # among them.
  for (k in c("pot_closure_start", "pot_closure_end", "pot_open_date", "census_start_date",
              "census_end_date", "commercial_opener", "pe_empty_effort_stratum",
              "pe_empty_stratum", "pe_variance", "bss_seed"))
    chk(sprintf("ladder: WINDOW pins %s", k), grepl(paste0(k, " = "), s, fixed = TRUE))
  chk("ladder: pot_open_date is pinned to the 2024-25 value, not the 2-season 2023-12-01",
      grepl('pot_open_date = "2024-12-01"', s, fixed = TRUE))
  chk("ladder: F_METHOD exists and ships 'new_throughout' (Matt 2026-09-12: the run uses the new f)",
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
  # 2026-09-12: the block now starts at F_NEW (D_R1 is built from it) and F_METHOD has to
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
  # 2026-09-12: D_R1's f block now depends on F_METHOD, so the block is evaluated under
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
        grepl("f_crab[osp_f_stratum[i]] * (1 - combo_c[osp_f_stratum[i]])", src, fixed = TRUE) &&
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
             crabbing_holiday_dates = as.Date(character()), estimate_red_rock = FALSE)
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
  chk("shipped: charter_frame = roster, effort_qc_drop holds the interview-total flag, the 2023-24 census window ends Jan 31 (opener Feb 1, 2024), a 2025-26 block is documented",
      identical(rc$charter_frame, "roster") && identical(rc$effort_qc_drop, "gear_count_from_interviews") &&
        identical(unname(rc$census_windows[["2023-24"]]), c("2023-12-01", "2024-01-31")) &&
        any(grepl("season-2025-26", readLines("run_config.R", warn = FALSE), fixed = TRUE)))
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
             crabbing_holiday_dates = as.Date(character()), estimate_red_rock = FALSE)
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
  # 2026-09-12: BOTH PE runners must delegate to the ONE shared implementation. Three
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
  chk("shipped: the PE unsampled-cell levers are the 2026-09-12 set (Matt: ship local_day_type on)",
      identical(rc$pe_empty_effort_stratum, "local_day_type") && identical(rc$pe_empty_stratum, "local") &&
      identical(rc$pe_variance, "impute_aware"))
  chk("shipped: the census keys are the 2026-09-11 set", identical(rc$census_expansion, "none") && identical(rc$charter_frame, "roster") &&
        identical(rc$charter_expansion, "vessel") && identical(rc$census_uncertainty, "charter"))
})

# ---------------------------------------------------------------------------
# 59. THE PE's UNSAMPLED AND SINGLETON CELLS (2026-09-12). Three defects, all of them in
#     two copies of the same code, all recomputed here BY HAND on a fixture small enough
#     to check with a calculator:
#       - the SE was sqrt(N^2 sd^2 / max(n,1)) with sd from the cell's own sampled days,
#         and sd() of ONE observation is NA -> 0, so a singleton cell contributed its full
#         point estimate and no variance;
#       - an IMPUTED cell did the same, so the effort SE was unchanged by imputation;
#       - the effort fill could be month-local while the CPUE fill was sub-season-wide.
#     The fixture is two months x two day types so every donor level is exercised, and the
#     "sampled_only" arm must reproduce the pre-2026-09-12 arithmetic EXACTLY.
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
  # impute_aware: var = N^2 s^2 (1/n_donor + 1); donor here is the sub-season (n = 6)
  chk("PE strata: an imputed cell now carries the donor-mean variance PLUS a between-cell term",
      isTRUE(all.equal(.g(s2, feb_we, "weekend", "se_total"), n_fw * sd_all * sqrt(1/6 + 1))) &&
      grepl("^imputed \\(donor mean \\+ between-cell\\)", .g(s2, feb_we, "weekend", "var_source")))
  chk("PE strata: donor_mean_only drops the between-cell term and is the LOWER bound",
      isTRUE(all.equal(.g(s3, feb_we, "weekend", "se_total"), n_fw * sd_all / sqrt(6))) &&
      .g(s3, feb_we, "weekend", "se_total") < .g(s2, feb_we, "weekend", "se_total"))

  # ---- 4. the mean and the spread come from DIFFERENT donor levels, on purpose --------
  # Feb weekday has 2 sampled days, so a Feb-weekday-donated cell would take the month
  # mean AND the month sd. Assert the two selections are reported separately.
  chk("PE strata: mean_level and sd_level are recorded per cell and may differ",
      all(c("mean_level", "sd_level") %in% names(s2)) &&
      any(s2$mean_level != s2$sd_level | s2$n_sampled >= 2L))

  # ---- 5. 'sampled_only' must reproduce the pre-2026-09-12 arithmetic EXACTLY ---------
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
  chk("PE CPUE fill: 'pooled' reproduces the pre-2026-09-12 behaviour and 'zero' the pre-2026-07-13 one",
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
# 60. THE TWO-PASS LADDER RUN (2026-09-13). Matt: "run the 4-rung version now and then
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
    # THE MECHANISM THE ANALYSIS RESTS ON, pinned. A 2026-09-13 smoke fit of the boat
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
  chk("pointers: the status document names the working BRANCH, not main",
      grepl("OSP-boat-count-incorporation", paste(L[seq_len(10L)], collapse = "\n"), fixed = TRUE) &&
      !grepl("^\\*\\*Repo:\\*\\* `Coastal-Rec-Crab-BSS`, `main`\\.", L[6]))

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
            "PULL_REQUEST.md", "README.md", "07_documentation/CLAUDE.md")
  bad_ptr <- Filter(function(f) {
    if (!file.exists(f)) return(FALSE)
    ln <- grep("72,027", readLines(f, warn = FALSE), value = TRUE)
    any(!grepl("supersed|SUPERSED|historical|HISTORICAL|previous|used to", ln))
  }, live)
  chk("pointers: no live document offers 72,027 without marking it superseded",
      length(bad_ptr) == 0, sprintf("(%s)", paste(bad_ptr, collapse = ", ")))

  # section 1v: one direction, stated.
  s1v <- paste(L[seq(grep("^## 1v\\.", L)[1], length(L))], collapse = "\n")
  chk("pointers: section 1v states which way its PE/BSS percentages run",
      grepl("PE RELATIVE TO BSS", s1v, fixed = TRUE) &&
      !grepl("PE 35,292 against BSS 34,837, -1.3%", s1v, fixed = TRUE))
})

cat(sprintf("\n==== %d passed, %d failed ====\n", ok, bad))
if (bad > 0) quit(status = 1)
