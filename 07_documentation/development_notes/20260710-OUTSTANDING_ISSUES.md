# Outstanding Issues: Coastal Rec Crab BSS Pipeline

**Last updated:** 2026-07-10
**Repo state at audit:** commit `564071b` ("7/9 gear resolved run"), plus the uncommitted P0/P1/P2 changes delivered 2026-07-10 (a run using those was in flight at audit time).
**Scope:** three pipelines (pooled CPUE, gear-resolved CPUE, weather/tide covariate module), plus the orchestrator and repo hygiene.

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

## Suggested next actions (all safe during a live run: `.Rmd` and `03_R_functions/` only)

1. **POOL-5** — wire `write_cpue_diagnostics` into pooled. Diagnostic-only, changes no estimate, and will likely reorder the rest of this list by revealing whether pooled's units are broken.
2. **POOL-6** — port pooled onto the shared gate and AR selector; delete inline copies.
3. **ORCH-25** — create `renv.lock`.
4. **POOL-1** — the `T_A_int` fix, specified but held until sign-off because it moves the publication boat number.
