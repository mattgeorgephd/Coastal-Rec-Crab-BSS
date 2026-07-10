# Outstanding Issues: Coastal Rec Crab BSS Pipeline

**Last updated:** 2026-07-10
**Repo state at audit:** commit `564071b` ("7/9 gear resolved run"), plus the uncommitted P0/P1/P2 changes delivered 2026-07-10 (a run using those was in flight at audit time).
**Scope:** three pipelines (pooled CPUE, gear-resolved CPUE, weather/tide covariate module), plus the orchestrator and repo hygiene.

**Update 2026-07-10 (v7.5, branch `pooled-CPUE-fixes`):** POOL-6, POOL-5, POOL-4, and POOL-2 are implemented (see the per-item **Resolution status** block below the pooled table).

**Update 2026-07-10 (after the gear-resolved 7/10 run, commit `d1968f3`):** the gear-resolved section below now carries a per-item **Status after the 2026-07-10 run** block (GR-7..GR-15 verdicts plus new items GR-16/GR-17). That run also re-opens the POOL-1/POOL-3 decision: it validates the boat deployment scale and shows the pooled boat is ~30% too high, so POOL-1 + POOL-3 are **reassessed from held to RECOMMENDED** (see their reassessment note below the pooled table).

This file is a living backlog. Each issue has an ID, a severity, the concrete evidence (file and line where possible), and enough mechanism to act on it without re-deriving it. Severity is about consequence for a publishable harvest estimate, not effort.

---

## Operational warnings (read before editing during a run)

- **Never edit `02_stan_models/*.stan` while a run is in progress.** `rstan::stan(file = ...)` is called once per fit inside the loop and re-reads the file each time, so an edit mid-run makes later fits compile a different model than earlier ones, and `stan_data` may no longer match. Editing `.Rmd` (knitr-parsed at start) and `03_R_functions/*.R` (sourced at start) is safe.

- **P2 (`use_mu_hier`) is on watch for the gear-resolved run.** Pooled attempted the identical change (v6.9 / B1.7) and reverted it because the shore all-gear fit hung; see GR-10 and POOL-4 below. Gear-resolved is expected to be fine because its shore AR is monthly/biweekly (P_n = 10 / 6), not daily (P_n = 289), but if `shore_all_gear` saturates `max_treedepth` or does not finish, P2 is the first suspect.

---

## Pooled model (publication priority)

**File:** `01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd`, `02_stan_models/crab_bss_pooled.stan`

Note two things pooled does NOT have, to avoid re-investigating them:
- Pooled's `run_pe` uses `str_detect(population_name, "private_boat")`, not `==`, so it is immune to the argument bug that hit gear-resolved (P0). This is why pooled's PE never showed that symptom.
- Pooled's only guarded prior is `R_T`, declared `real<lower=0,upper=1>`. A flat prior on bounded support is proper, so pooled does NOT have the `R_G_boat` improper-prior bug that F1 fixed in gear-resolved.

| ID | Severity | Issue |
|----|----------|-------|
| POOL-1 | High | **Dimensional inconsistency in the trailer expansion.** `T_A_int = rep(1L, nrow(intA))` is a vector of literal ones, so `T_A_int[a] ~ bernoulli(R_T)` contributes `R_T^n` to the likelihood; combined with `beta(5,1)` this pins `R_T` at 1.00 (posterior 0.98-1.00). `R_T` is asked to be both "probability a group has a trailer" (approximately 1) and the trailer-to-`lambda_E` expansion factor. Once pinned, `lambda_E` is groups while `h = gear_time_total` is gear-hours, so `E_sum` is group-hours mislabelled as gear-hours (off by gear-per-group, approximately 4). Pooled's boat number is right only through this compensating error. Fix: pass a genuine trailer indicator, or drop the Bernoulli term. Needs its own before/after validation because it moves the publication boat number. Evidence: Stan L112, L167, L334; `expansion_ratios.csv` shows `R_T` at the bound. |
| POOL-2 | High | **No incomplete-trip filter anywhere.** `trip_status` appears exactly once in the driver (its definition); `filter_incomplete_trips` appears zero times. 39.8% of shore all-gear interviews are incomplete, with a measured CPUE bias of about -20% (gear-dependent: Pot -21%, Trap -23.2%, Snare -19.9%, Ring Net +4.4%). Pooled's shore CPUE, hence the publication-priority estimates, is biased low. Gear-resolved filters these in four places; pooled in none. |
| POOL-3 | High | **Boat `L = 24` (flat gear-hours).** `L_data_vec <- rep(24.0, D)` at driver L868. The saturation analysis (`crab_per_gear ~ h^0.133`, from 1,532 boat interviews) shows time-denominated effort is invalid for pots. Same defect F2 fixed in gear-resolved. Should move to the gear-deployment scale with `L = tau`. |
| POOL-4 | Medium | **`sigma_mu_E` / `sigma_mu_C` funnel at G*S = 1**, known and currently unfixed. B1.7 (v6.9, 2026-06-21) collapsed the single-cell hierarchy and it cleared the boat funnel offline, but in production the shore all-gear fit (daily AR, 289 days, ~50% unobserved) hung beyond 24 h: removing the decoupled level term forced the level to reconcile directly against the high-dimensional AR ridge, saturating `max_treedepth`. Reverted to v6.8 structure. The durable fix is a more informative boat effort series, not parameter surgery. Evidence: Stan L48-60, L234. |
| POOL-5 | Medium | **No F4 diagnostics wired.** `write_cpue_diagnostics`, `bss_effort_spec`, and the linearity/saturation/estimator-triad checks are absent (0 call sites). Pooled cannot currently detect the effort-unit defect that these found in gear-resolved. This is the highest information-per-hour item: diagnostic-only, changes no estimates, and would immediately show whether pooled's boat and shore units are as broken as gear-resolved's were. |
| POOL-6 | Medium | **R-layer duplication.** Pooled defines its own `compute_gate` (L1487), `use_pe_for` (L1616), and inline AR-resolution logic (L811-829), while gear-resolved uses the shared `03_R_functions/` versions. Two implementations of the same logic, guaranteed to drift. Port pooled onto the shared modules and delete the inline copies. |

### Resolution status — v7.5 (2026-07-10, branch `pooled-CPUE-fixes`)

Four of the six pooled items are implemented; two are held by decision. Every change carries inline `POOL-<n>` comments. The pipeline could not be run offline (a full pooled run is a multi-hour rstan job), so the changes were validated by R static parse checks and an adversarial code review, **not** by a before/after harvest comparison. Re-run the pooled model to refresh the numbers and emit the new diagnostics.

- **POOL-6 — DONE (behavior-preserving).** Deleted the inline `compute_gate()`, `use_pe_for()`, and adaptive-AR block in `prep_bss_crab`; the driver now calls the shared `bss_compute_gate()`, `bss_use_pe_for()`, and `bss_select_ar_resolution(..., fixed_resolution = NULL)`. Identical gating quantities (`C_expected_sum`, `E_sum`) and thresholds (`params$max_divergence_total_impact_sd`, `params$max_divergence_fraction`); the report-only `B1_C_rhat` column is re-attached at the call site so `convergence_report.csv` keeps its schema. No estimate change.
- **POOL-5 — DONE (diagnostic-only).** `prep_bss_crab` attaches `attr(stan_data, "cpue_data")` and the `.effort_unit` / `.h_unit` tags; the extended-diagnostics loop calls `write_cpue_diagnostics()` (tryCatch-wrapped). Writes `cpue_estimators_`, `cpue_saturation_`, `cpue_linearity_<label>.csv` per fit. Unit labels are set equal per population (shore crabber-hours; boat gear-hours) so `bss_assert_effort_units()` passes and the run is not halted; the saturation/linearity CSVs carry the real signal (catch is roughly flat in soak time for pots). No estimate change.
- **POOL-4 — DONE (lever, default off).** Added `int<lower=0,upper=1> collapse_mu_hier;` to `crab_bss_pooled.stan` and guarded the `mu_E`/`mu_C` computation; the driver passes it from `params$collapse_mu_hier` (default `FALSE`; accepts a per-population named list). `collapse_mu_hier = 0` reproduces the v6.8 hierarchy exactly (default posterior unchanged); `= 1` collapses the single-cell level (the B1.7 experiment) and can now be tested per population from config without editing Stan mid-run. The durable funnel fix remains a more informative boat effort series; this is only a safe lever. Note: the `= 1` path is syntax/review-validated, not run-validated.
- **POOL-2 — DONE (default ON, toggleable). CHANGES THE SHORE NUMBER.** Added `params$filter_incomplete_trips` (default `TRUE`) and applied `trip_status == "Complete" | is.na(trip_status)` to the PE CPUE (`run_pe`) and the BSS CPUE (`prep_bss_crab`, which propagates to the gear-per-group `intA` set), matching gear-resolved. New section 7.2b writes `sensitivity_incomplete_trips.csv` (PE catch filter-off vs on) to quantify the move. Removing incomplete trips (~40% of shore all-gear interviews, ~-20% CPUE bias) raises the shore CPUE, so the shore harvest estimate rises. Set to `FALSE` for pre-POOL-2 behavior.

**POOL-1 + POOL-3 — held 2026-07-10, now REASSESSED as RECOMMENDED after the 2026-07-10 gear-resolved run.** One entangled fix that moves the publication boat number (pooled boat BSS all-gear ~56,266 in the 20260708 run). It was held pending a validation run; the gear-resolved 7/10 run supplies that validation and shows the pooled boat is ~30% too high (reassessment below). Confirmed evidence and the implementation spec:

- *Evidence.* In `05_output/20260708/pooled-CPUE`, boat `R_T` posterior is mean 0.993 / median 0.995 / 95% CI [0.975, 1.000] (pinned at the bound), exactly as POOL-1 predicts. With `R_T ~ 1`, `lambda_E` counts boat GROUPS while `L = 24` and `h = gear_time_total` are gear-hours, so `E_sum` is group-hours mislabelled as gear-hours; the boat total is right only through the compensating cancellation of "groups (~4x too few)" against "L = 24 h vs a true ~6 h soak (~4x too many)".
- *The fix (port the boat onto the gear-resolved structure).* In `crab_bss_pooled.stan`: replace parameter `R_T` with `R_G_boat ~ lognormal(log(4), 0.5)`; change the trailer likelihood to `T_I[i] ~ neg_binomial_2(lambda_E / R_G_boat, r_E)` (so `lambda_E` becomes GEAR and `lambda_E / R_G_boat` is boat groups); replace `T_A_int[a] ~ bernoulli(R_T)` with `Gear_A_boat[a] ~ poisson(R_G_boat)`; rename the data `T_A_int` / `A_A_trailer` to `Gear_A_boat` (gear per interviewed group); update the trailer `log_lik` and the `R_T_out` generated quantity; and optionally activate the boat I/E stream on the group scale (`ie_pred / R_G_boat`). In the driver `prep_bss_crab` (boat branch): `h = number_of_gear` (deployments), `L_data = rep(tau_boat_prior_mu ~ 1.2, D)`, `estimate_L = 1`, `L_prior_sigma = rep(tau_boat_prior_sigma ~ 0.3, D)`, `Gear_A_boat = as.integer(intA$number_of_gear)`, and `.effort_unit = .h_unit = "gear-deployments"`. This mirrors F1 + F2 on `crab_bss_gear_resolved.stan`; `03_R_functions/bss_trailer_expansion.R` already handles both conventions. Land it behind a boat-mode flag defaulting to current behavior, run pooled with the flag off then on, and adjudicate with the POOL-5 `cpue_saturation_*` / `cpue_linearity_*` diagnostics plus the boat catch total before adopting it as the default.
- *Reassessment (2026-07-10, after the gear-resolved 7/10 run) — now RECOMMENDED and de-risked.* The gear-resolved run supplies the validation that was missing when this was held. (1) The boat CPUE-saturation diagnostic is model-independent proof the pooled unit is wrong: `cpue_saturation_private_boat_all_gear` gives `beta = 0.249` (95% CI 0.13-0.36), i.e. boat catch barely rises with soak time, so gear-hours is not a valid denominator for pots. (2) The correctly-scaled gear-resolved boat (gear-deployments, `beta_h = 0.75`, model CPUE at/below ratio-of-sums, no drift) reports BSS all-gear catch **43,314**, against the pooled boat's **56,266** on the broken gear-hours scale: the pooled boat is **~13,000 / ~30% too high**, so the compensating error does NOT cancel. (3) The port is now low-effort and low-risk: `03_R_functions/bss_effort_spec.R` plus the `effort_scale_gear` / `R_G_boat` / `tau` machinery are shared modules, and gear-resolved runs the exact target structure cleanly (boat all-gear 3/3 gates pass, monthly AR, R-hat ~1.0, 21 divergences). **Recommendation:** implement POOL-1 + POOL-3 as one paired fix reusing `bss_effort_spec.R`; expect the pooled boat to fall ~25% toward ~43k; move the pooled boat to monthly AR (as gear-resolved does) and keep the POOL-4 `collapse_mu_hier` lever ready for the single-cell funnel; confirm with one pooled run (now a confirmation, not an exploration). Caveat: pooled and gear-resolved also differ in AR resolution and mu-structure, so the exact pooled magnitude is fixed only by that run, but a ~30% overestimate on a publication boat number is too large to leave, and the saturation diagnostic proves the unit itself is invalid. The related shore-unit question (GR-16) affects the pooled shore the same way and should be settled by the P1 LOO comparison in the same session.

---

## Gear-resolved model

**File:** `01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd`, `02_stan_models/crab_bss_gear_resolved.stan`

The P0/P1/P2 fixes (population argument, ratio-of-sums strata, PE assertion, configurable effort unit with LOO comparison, mu-hierarchy collapse) were delivered 2026-07-10 and are not listed here as open. The items below remain.

| ID | Severity | Issue |
|----|----------|-------|
| GR-7 | Medium | **Option A (genuine per-gear CPUE) unbuilt.** `G = 1`, so the per-gear machinery in the Stan model is inert and gear-type catch comes from PE apportionment. Requires effort shares in the `O[d,s,g]` offset fed by `pi_gear_data`, plus a rule for multi-gear interviews (the gear-weight matrix gives fractional weights, but `c[a] ~ NB2(lambda_C[gear_IntC[a]] * h, r_C)` needs one gear per interview). Cannot simply raise `G`: the only effort observation touches `g = 1` while `E_sum`/`C_expected_sum` sum over all `g`, so unobserved effort processes would enter the totals identified by priors alone. |
| GR-8 | Medium | **`sigma_IE` funnel in `shore_ring_net`.** Largest divergence-localization SMD, -0.420, with only 2 in-window I/E days. `exponential(5)` has its mode at zero, so `sigma_IE` can shrink and stiffen the lognormal likelihood. Add an `ie_min_obs_shore` guard (>= 3), or tighten the prior to something like `lognormal(log(0.3), 0.5)` that keeps mass off zero without a hard bound. |
| GR-9 | Medium | **`shore_all_gear` `sigma_IE` = 0.603 (0.334-1.083), unexplained.** The 4 I/E observations and the effort counts disagree by roughly 60% on the log scale. Investigate whether the I/E days are peak-count days and whether a monthly-AR `lambda_E` can represent a day-specific total. |
| GR-10 | Medium | **No `collapse_mu_hier` rollback lever.** P2's `use_mu_hier` is computed in `transformed data` as `G * S > 1`, so it cannot be toggled from config. Given the pooled B1.7 hang (POOL-4), make it a data flag defaulting to current behavior so it can be turned off without editing Stan mid-investigation. |
| GR-11 | Medium | **Gate policy: hard 5% divergence backstop overrides the scale-aware impact test.** On the 2026-07-09 run, `shore_ring_net` passed R-hat, `n_eff`, and the impact test (divergences moved totals by 0.004 SD) but failed only `pass_div_fraction` (7.5% vs 5%), so a good BSS estimate (6,930, matching pooled's 6,733) was replaced by a bad PE estimate (12,174). Keep the backstop on principle (divergences indicate biased exploration), but only after the PE is validated. Revisit the 5% threshold. |
| GR-12 | Medium | **`tau_boat = 1.2` rests on 2 WBL I/E days** (turnover 1.00 and 1.29). Boat catch is directly proportional to it. The boat I/E stream is built and activates automatically once `ie_min_obs_boat` (2) days fall inside a window; the egress-classification pilot is the critical path. State this dependence in any co-management number. |
| GR-13 | Low | **Dead `stan_data` entries and fragile metadata.** `estimate_R_G` and `R_G_fixed` are passed but not declared (retained deliberately for Option A). Dot-prefixed metadata (`.effort_unit`, `.h_unit`, etc.) is passed into `rstan::stan(data=)` and tolerated; `.cpue_data` was moved to an attribute for this reason, but the remaining dot-vectors are fragile. |
| GR-14 | Low (by design) | **Latent AR carries only 4-13% of effort variance;** NB observation overdispersion carries 81-91%. Consequence of the coarse fixed AR (`ar_adaptive = FALSE`, biweekly/monthly). Day-level variance lands in `r_E`. Not a bug, but the reason the AR contributes little. |
| GR-15 | Low | **`private_boat_ring_net` never fits** (insufficient data), always falls back to PE. |

### Status after the 2026-07-10 gear-resolved run (commit `d1968f3`)

The 7/10 run's code changes were P0/P1/P2 (population-argument fix + ratio-of-sums strata + a PE assertion; the configurable effort-unit module `bss_effort_spec.R`; and the `mu`-hierarchy collapse at `G*S == 1`). None of GR-7..GR-15's own proposed fixes were implemented, so the movement in the outputs is a side-effect of P0/P2, not a targeted GR fix. Per-item status against the run (comparing `05_output/20260710` vs `05_output/20260709`):

- **GR-7 — STILL OPEN.** `G = 1`; `catch_by_gear_type.csv` has `has_BSS = TRUE` only for `"All"` (0.844); Pot/Ring Net/Snare/Trap remain PE-apportioned. No `pi_gear_data` gear-share offset in the Stan diff.
- **GR-8 — PARTIALLY ADDRESSED (symptom only).** `shore_ring_net` now converges and reports BSS: divergences 598 -> 216, fraction 7.47% -> 2.7%, `pass_div_fraction` FALSE -> TRUE, catch PE 12,174 -> BSS 6,936. But the fix did not land: no `ie_min_obs_shore` guard, `sigma_IE ~ exponential(5)` unchanged, and `sigma_IE` is still the #1 divergence driver (SMD -0.444, marginally worse). The gain came from P2 collapsing the separate `sigma_mu` funnel, which pulled total divergences under the 5% backstop.
- **GR-9 — STILL OPEN.** `shore_all_gear` `sigma_IE` median 0.600 (0.331-1.100) vs 0.603 before: unchanged, still unexplained.
- **GR-10 — STILL OPEN.** P2's collapse is still hardcoded as `use_mu_hier = (G*S > 1)` in `transformed data`; not a config flag. (The pooled side got exactly this lever as POOL-4's `collapse_mu_hier`; gear-resolved still needs the same treatment.)
- **GR-11 — STILL OPEN as policy.** `max_divergence_fraction` is still 0.05 and `bss_convergence_gate.R` was not in the commit. The 7/9 symptom disappeared incidentally: P2 dropped ring-net divergences below 5% and P0 corrected the PE (ring-net PE 12,174 -> 6,152, now agreeing with BSS 6,936). A fit at 5.1% with negligible impact would still be force-flipped to PE.
- **GR-12 — STILL OPEN.** `tau_boat` still 1.2; boat I/E did not activate (`n_ie_obs = 0`, boat `sigma_IE` at its prior). Boat catch still rests entirely on the `tau_boat = 1.2` prior.
- **GR-13 — STILL OPEN.** `estimate_R_G` / `R_G_fixed` still passed-but-undeclared; the dot-metadata vectors still passed into `rstan::stan(data=)`. No cleanup in the commit.
- **GR-14 — STILL TRUE (by design).** `effort_overdispersion_decomp`: latent share 0.043 / 0.065 / 0.133; NB-overdispersion share 0.91 / 0.90 / 0.80. Unchanged.
- **GR-15 — STILL OPEN.** `private_boat_ring_net_only` is `PE (insufficient data)` (~17 interviews Sep-Nov). Unchanged.

New items the run surfaces:

| ID | Severity | Issue |
|----|----------|-------|
| GR-16 | High | **Shore runs on the effort unit the pipeline itself flags as invalid, and the P1 LOO decision was never run.** `shore_effort_unit = "crabber-hours"` and `loo_effort_unit_comparison = FALSE`, yet `cpue_linearity_shore_all_gear` `beta_h = 0.571` (CI 0.500-0.641, flag TRUE) and `cpue_saturation` `beta = 0.217` (flag TRUE); the estimator triad sits ~3/4 toward mean-of-ratios (`estimator_drift_flag = TRUE`). The P1 machinery (three candidate units, `bss_effort_h_candidates`, elpd_loo) is built but unexercised, so the production shore totals (19,684 all_gear + 6,936 ring_net) rest on a unit flagged as wrong. Next run: set `loo_effort_unit_comparison = TRUE` and compare crabber-hours vs gear-hours vs gear-deployments by `loo_summary_*`. This is the shore analogue of the boat's POOL-1/POOL-3 and affects the pooled shore too. |
| GR-17 | Low | **Cosmetic: collapsed `sigma_mu_E/C` posteriors look non-converged but are inert.** Post-P2 they sample their half-Cauchy prior (mean up to ~21-30, hi95 ~57) and enter no likelihood. Harmless, but drop them from `structural_params` when `use_mu_hier == 0`, or flag them, so a reviewer does not misread them as divergent. |

**Run headline (task 1 comparison).** Convergence 2/3 -> 3/3 fits pass. The shore ring-net estimate moved from an inflated PE 12,174 to BSS 6,936 (matching the pooled-scale expectation ~6,733). Shore all-gear BSS stable at ~19,684; boat all-gear BSS stable at ~43,314 on validated gear-deployments. PE totals dropped sharply from the P0 fix (shore PE 45,112 -> 19,767; PE total 80,160 -> 54,598; the old PE ran on the wrong scale), but the gated BSS numbers, not the PE, are what the port total reports.

---

## Weather / tide covariate module

**File:** `06_diagnostics/BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd` (2,225 lines, 27 inline function defs), `02_stan_models/crab_bss_pooled_weather_adjusted.stan`

| ID | Severity | Issue |
|----|----------|-------|
| WX-16 | High | **It forks the engine.** Private copies of `prep_bss_crab_augmented` (L1291) and `bss_convergence_check` (L1668), among others. The header claims "pooled-model v6.9 parity"; pooled is at v7.4 plus the shared-module refactor. Every fix to pooled has to be manually re-applied here or the two diverge silently. |
| WX-17 | High | **Its Stan model is stale.** `crab_bss_pooled_weather_adjusted.stan` was last touched Jul 8 and has received none of B1.8, F1, F4, P1, or P2. It needs the same guarded-prior and G*S = 1 funnel audit as the other two models. |
| WX-18 | Medium | **Not run since the day-length module was extracted.** It consumes `dwg`, `ie_data`, and `L_eff_model` from the pooled render environment (the Option A shared-environment hand-off). Those objects still exist, but the coupling is untested since the extraction. |
| WX-19 | Medium | **Leave-one-week-out block CV outstanding.** Standard PSIS-LOO cannot test genuine week-long sampling gaps. `kfold_time_block_cv` exists (L1506) but the backlog notes it is not yet the right block structure for the covariate question. |
| WX-20 | Medium | **Multi-year covariate pooling needs a protocol fixed effect.** The pre-2024-25 protocol (single peak count) and the 2024-25 protocol (three randomized counts) define different estimands; naive pooling is invalid. |

---

## Orchestrator, reporting, repo hygiene

| ID | Severity | Issue |
|----|----------|-------|
| ORCH-21 | Medium | **Weather re-run requires a full pooled refit.** The Option A shared-environment hand-off means pooled's memory is not freed before the weather module runs, and the weather module cannot be re-run without re-running the approximately 3 h pooled fit. The Phase-2 disk-bundle design (persist `dwg`/`ie_data`/`L_eff_model`) is specified but not built. |
| ORCH-22 | Low | **`season_summary.csv` mixes presentations.** A PE-based mode breakdown alongside a gate-combined total. |
| ORCH-23 | Low | **`port_total` Effort row sums incompatible units** (shore crabber-hours + boat gear-deployments). Flagged in the file via `effort_units_note`; arguably should not be summed at all. |
| ORCH-24 | Low | **`catch_by_gear_type.csv` mixes sources.** A boat "All" BSS column alongside shore gear types from PE apportionment. Given the PE and apportionment caveats, doubly suspect. |
| ORCH-25 | High (for publication) | **`renv.lock` absent.** No dependency pinning for a method intended for publication. Roughly a five-minute job with a large reproducibility payoff. |
| ORCH-26 | Medium | **Output management.** 53 MB and 714 git-tracked files under `05_output/`, including pre-v6 runs that cannot be reproduced from current code. Policy undecided (prune vs archive vs stop tracking). |
| ORCH-27 | Low | **Gear-resolved does not log rstan/StanHeaders versions** the way pooled does (pooled L2317). The orchestrator manifest covers this via `sessionInfo()`, so it is belt-and-suspenders. |

---

## Cross-cutting principle worth keeping visible

For trap and pot gear, catch per unit soak time is not a stable parameter, so any effort unit denominated in time will be unstable. This drove F2 (boat) and is the open question for shore (P1, decided by LOO). The pipeline now measures it automatically every run via `cpue_linearity_*.csv` and `cpue_saturation_*.csv`, so it can be checked each season rather than rediscovered. Any new gear type or population should be run through those diagnostics before its totals are trusted.

---

## Suggested next actions

POOL-5, POOL-6, POOL-4, and POOL-2 are DONE (v7.5, branch `pooled-CPUE-fixes`). Remaining, in priority order:

1. **POOL-1 + POOL-3 (now recommended)** — port the pooled boat onto the `R_G_boat` / `tau` deployment scale, reusing the shared `bss_effort_spec.R`. The 7/10 gear-resolved run validates the target structure and shows the pooled boat is ~30% (~13k) too high; fully specified in the reassessment note above. Do this together with the v7.5 re-run.
2. **Re-run pooled under v7.5** to refresh the published totals (POOL-2 raises the shore number) and to emit the new `cpue_saturation_*` / `cpue_linearity_*` diagnostics; read `sensitivity_incomplete_trips.csv` to see the shore move. Combine with (1) so the boat and shore corrections land in one validated run.
3. **GR-16 / P1 shore effort unit** — set `loo_effort_unit_comparison = TRUE` for gear-resolved and pick the shore unit by `elpd_loo`; the same question affects the pooled shore. Highest information-per-hour of the open items.
4. **ORCH-25** — create `renv.lock` (about a five-minute job with a large reproducibility payoff for a publication method).
5. **Gear-resolved GR items** — GR-8 (`ie_min_obs_shore` guard / `sigma_IE` prior), GR-10 (make the `mu`-collapse a config flag, as POOL-4 did for pooled), GR-11 (gate-policy revisit) remain open; the 7/10 run cleared GR-8's symptom only as a P2 side-effect.
