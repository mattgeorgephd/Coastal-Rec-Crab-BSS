# Grays Harbor Crab Harvest Estimation: Gear-Resolved CPUE Model

## Development History

**Companion to:** `BSS-GH-gear-type-CPUE-model-documentation.md` (the published gear-resolved method reference).
**Scope:** the full version-by-version change log of the gear-resolved-CPUE pipeline and its Stan model `crab_bss_gear_resolved.stan`, plus the detailed working notes from the two largest episodes: the boat effort-scale correction (fix-markers F1 and F2) and the point-estimator population-and-estimator fix (P0). The current production state is framework **v5.5** (2026-07-11), in which both the shore and boat components run on the gear-deployment effort unit and `loo_effort_unit_comparison` is turned off for production.
**Convention:** no em dashes.

This file is the provenance record for the gear-resolved model. The published method document summarizes this history in one screen and refers here for the detail. Entries are newest-first. The version log traces the framework from the initial gear-resolved release (v5.0) through the empirical-proportion and data-alignment work (v5.1 to v5.4) to the shore effort-unit resolution (v5.5). The run-driven Stan fixes that do not carry a framework tag are recorded in their own section below the version log, and the outstanding backlog closes the log.

The gear-resolved track branched from the shared pooled/gear-resolved sequence at v5: versions v1 through v4 are one shared milestone line (documented in `README.md`), and since v5 the two tracks carry independent change logs. The gear-resolved track imports the pooled convergence work as parity ports rather than re-deriving it (the B1.x markers in the Stan file below), and adds its own hard-won equivalents, F1 and F2, which the pooled track then ported back as POOL-1 and POOL-3. The two tracks share the effort-specification module `03_R_functions/bss_effort_spec.R` and the convergence-gate and AR-selection helpers, so they cannot drift onto different scales or gate policies.

**A note on the two numbering systems.** The gear-resolved track carries a v5.x "framework" series in the RMD header (`BSS-GH-gear-type-CPUE-model.Rmd`) and in the documentation change log. That series is distinct from the Stan model's own history: `crab_bss_gear_resolved.stan` carries no vX.Y tag. The documentation's "Stan v3.2" label (attached to framework v5.2) is stale; the Stan file is materially past v3.2. It now carries the pooled-track parity ports B1.3, B1.5, B1.6, B1.8, and B1.9, plus the run-driven fixes F1, F2, P1, and P2. Read the framework tag for what the R pipeline does and the fix-marker section below for what the Stan model does; the "Stan vX.Y" labels in the version log are historical and should not be trusted as the current model state.

**Header status.** As of v5.6 (2026-07-12) the RMD banner, the method-document header, and this history are all in sync at v5.6. Earlier releases carried a header lag (the RMD banner and method-doc header read v5.4 while the body was v5.5); that is resolved.

------------------------------------------------------------------------

## Version log

### v5.6 (2026-07-12), run_config base-parameter restructure (P5), parity with pooled v7.9

Config-architecture parity with the pooled track's v7.9. No Stan change, no estimate change, and no gear-resolved run was required.

-   **`run_config.R` is the base parameter set.** The driver merge was inverted to `params <- modifyList(run_config, params_model)`, so `run_config` is the base and this model layers only its internal tuning (Stan file, per-fit sampler settings, gate thresholds) on top. The two key sets are disjoint, so behavior is identical; the inversion encodes intent (run_config is authoritative). The gear-resolved AR resolution cap map now lives in `run_config.R` as `ar_max_resolution$gear_resolved = list(shore = "weekly", private_boat = "monthly")`, selected after the merge. That map is dormant in production (`ar_adaptive = FALSE`, fixed per-sub-season `period_bss`); it binds only if `ar_adaptive` is turned on.
-   The divergence-fraction backstop was already 0.05 here, so the pooled P7 tightening (0.15 to 0.05) brought the pooled model into line with this one rather than changing anything in the gear-resolved track.
-   Files changed: `01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd` (merge inversion, AR map read from `run_config`, header to v5.6); this history.

### v5.5 (2026-07-11), Shore moved to the gear-deployment effort unit

The shore component is moved onto the gear-deployment effort unit (`shore_effort_unit = "gear-deployments"`), matching the boat, so both fitted populations now run on gear-deployments. The unit was chosen from the 2026-07-10 shore LOO comparison (shore all-gear, n=1649):

-   **Linearity, not fit, decided it.** Gear-deployments is the only shore unit whose CPUE-linearity slope `beta_h` covers 1 (1.05, 95% CI 0.94-1.15), against crabber-hours 0.57 (0.50-0.64) and gear-hours 0.73 (0.67-0.80). A `beta_h` that covers 1 means catch scales one-for-one with the effort unit, so the constant-CPUE season expansion is unbiased; a `beta_h` below 1 means the unit over-counts effort at the high end and the expansion is biased. Gear-deployments is therefore the only harvest-unbiased shore unit.
-   **No estimator-triad drift.** Gear-deployments is also the only shore unit with no drift across the three CPUE estimators: ratio-of-sums 0.87, mean-of-ratios 0.85, and model 0.85 crab/deployment all agree. When the three disagree, the unit is mis-specified.
-   **The LOO edge was a red herring.** Gear-hours had a marginally higher catch-stream `elpd_loo` (-3131, versus -3190 for deployments and -3228 for crabber-hours), but that edge is `lambda_C` absorbing the sub-linearity into the fitted rate, which is exactly what biases the season expansion. Predictive fit and harvest-unbiasedness diverge here, and the unit was chosen on harvest-unbiasedness.
-   Resolves backlog item GR-16. Sets `loo_effort_unit_comparison = FALSE` for production (the comparison is done). Harvest-unbiasedness is the operative criterion because the reported number is a constant-CPUE season expansion: if `beta_h` is below 1 the fitted rate is pulled toward the high-effort interviews and the expansion to unsampled days over-counts, so a unit with `beta_h` covering 1 is the one whose season total is not systematically biased.
-   **Retroactive record.** This entry also, belatedly, records that the boat had already moved to gear-deployments (`h = number_of_gear`, `L = tau_boat`) via the shared `03_R_functions/bss_effort_spec.R` plus the Stan `effort_scale_gear` / `E_scale` machinery. That is the F2 and P1 work; it predates v5.5 but was never carried in this framework log. Documentation Section 6.2 is corrected accordingly.
-   Files changed: gear RMD (`shore_effort_unit`, `loo_effort_unit_comparison`); gear Stan (via the shared effort module); documentation (Section 6.2, Section 13).

### v5.4, R-hat convergence gate tightened to 1.01

-   Tightened the convergence gate from R-hat < 1.05 to R-hat < 1.01 (Vehtari et al. 2021), in step with the pooled track's v6.4 so both gates use one threshold. Outcome-neutral for the 2024-25 season: every passing fit has R-hat near 1.00, and every failing fit fails on divergences or n_eff, not on an R-hat between 1.01 and 1.05.
-   Files changed: gear RMD (convergence gate); documentation.

### v5.3, R-side data prep aligned with the (then) gear-hours Stan formulation

-   Aligned the R-side data preparation with the gear-hours boat formulation the Stan model then expected. For boat fits the model expects `Gear_A_boat[IntA_trailer]` (observed gear per boat group, replacing `T_A_int` / `A_A_trailer`), `h = gear-hours`, and `L[d] = 24`. The v5.2 data list still passed the obsolete `T_A_int` / `A_A_trailer` variables and used crabber-hours / day_length for boats, which produced the runtime error `Gear_A_boat[IntA_trailer] is missing` and, if patched alone, an approximately 2x downward bias in boat catch. `prep_bss_crab()` now branches by `is_shore` for `L`, `h`, and the gear-expansion block; boat interviews lacking a valid `gear_time_total` are dropped before the Stan data list is built.
-   NOTE: this gear-hours (`L = 24`) boat formulation was itself later superseded by the gear-deployment scale (see fix-marker F2 and framework v5.5). Documentation Section 6.2 now calls the v5.3 gear-hours choice an over-correction: it fixed the missing-variable crash but replaced a wrong unit (crabber-hours) with another wrong unit (gear-hours), and the right unit is gear-deployments.
-   Files changed: gear RMD (`prep_bss_crab`, boat data list); documentation.

### v5.2 (Stan v3.2), Expected-value catch, incomplete-trip filter, divergence gate

-   **Expected-value daily catch `C_gear`** (no Poisson noise) is emitted for the CPUE trajectories, alongside a separate `C_gear_pred` (Dirichlet plus Poisson) for the season totals and prediction intervals. Trajectories read the smooth expected value; totals carry full predictive noise.
-   **Dirichlet-sampled `pi_gear`** in generated quantities propagates gear-proportion uncertainty into the per-gear catch. Raw `n_weighted_gear` is passed as data.
-   **Incomplete-trip filter** (`filter_incomplete_trips`, default TRUE) removes the roughly -20% CPUE bias for pots and traps, whose soak-time gear is often not yet retrieved at interview.
-   **Divergent transitions added to the convergence pass/fail** (`max_divergences`, default 5), matching the pooled v6.1 gate. (Reconciled 2026-07-13 to the code and the method documentation, both of which use 5; this entry previously said 0. A3.)
-   **Regulatory gear exclusion per sub-season** (`gear_exclude`) prevents fitting a CPUE for a gear type prohibited in that sub-season.
-   **Minimum effective-N per gear raised 3 -> 15** (`bss_min_gear_effective_n`) as a guard against phantom multi-gear mentions.
-   Trip-completion sensitivity output added.
-   Files changed: gear Stan (v3.2); gear RMD.

### v5.1 (Stan v3.1), Empirical gear proportions and per-gear process error

-   **Weighted fractional gear assignment** for multi-gear interviews (equal split across mentioned gear), replacing the Pot > Ring Net > Trap > Snare priority hierarchy that had assigned each interview to a single gear.
-   **`pi_gear` becomes DATA** (empirical proportions with Laplace smoothing) rather than a sampled parameter; the categorical likelihood is removed. `pi_gear` varies by period x day_type.
-   **Per-gear-type process-error SD** `sigma_eps_C_gear[G_gear]`, so each gear's CPUE process has its own volatility.
-   **`R_G` conditionally estimated** (shore) or fixed at 1.3 (boat).
-   **Cauchy prior scales tightened 5 -> 2.**
-   `bss_max_interviews` removed (the model now fits the full dataset); biweekly BSS periods added for ring-net; extended convergence checks; word-boundary regex for gear classification; day-length cap raised 16h -> 17h (46.9N summer).
-   Files changed: gear Stan (v3.1); gear RMD.

### v5.0 (Stan v3.0), Initial gear-resolved model

-   Initial gear-resolved release. Per-gear-type AR(1) CPUE processes with a shared `phi` and a shared `sigma`; Dirichlet-distributed `pi_gear` per period; the priority hierarchy for multi-gear interviews (later replaced in v5.1); the holiday effect `B2` in the effort process; sparse per-observation overdispersion.
-   Files changed: gear Stan (v3.0); gear RMD (branched from the shared pooled/gear sequence).

### Shared-sequence milestones (v1 to v4)

Versions v1 through v4 predate the gear-resolved branch and are a single shared milestone sequence, documented in `README.md`: v1 the single-population dock-only prototype, v2 the shared bug fixes (CSV columns, Stan dimensions, output folders), v3 the three-population two-sub-season pooled model with convergence tuning, and v4 the dawn/dusk day length, stat-week PE, census dates, and team review. The gear-resolved model forked at v5.0; everything above this line is gear-resolved-specific.

### Run-driven fix-markers (F and P series)

The Stan file `crab_bss_gear_resolved.stan` carries no vX.Y tag; it records its changes as fix-markers in the header comment block. The 20260709, 20260710, and 20260711 runs drove the markers below. These are the changes that the framework v5.x series does not name, and they include the single largest quantitative move in the model's history (F2).

-   **F1 (= B1.10), proper unconditional prior on `R_G_boat`.** Before: `R_G_boat` was declared `real<lower=0>` with no prior and no likelihood in the shore fits (it is a boat-only parameter). The `<lower=0>` constraint adds a Jacobian to the target, so with nothing else touching the parameter the log-posterior increased without bound, `exp(z)` plus the Jacobian drove it to +inf, `R_G_boat` ran to about 1e307, and 97.6% of shore transitions diverged. After: a proper unconditional prior `R_G_boat ~ lognormal(log(4), 0.5)`, applied in every fit. Shore_ring_net divergences 7818/8000 (97.7%) -> 598/8000 (7.5%); shore_all_gear 7807/8000 (97.6%) -> 213/8000 (2.7%); R-hat 1.0002-1.0026, n_eff 2937-9313. See Appendix A.
-   **F2, boat effort measured in gear DEPLOYMENTS, not gear-hours.** Boat catch is not linear in soak time. Binned by hours-per-gear, crab per gear-HOUR falls 43x (2.88 -> 0.068) across a 64x soak range, while crab per gear per trip rises only 1.8x; a log-log fit over 1,532 interviews gives `crab_per_gear ~ h^0.133`. Gear-hours forces that 43x per-hour gradient onto one constant CPUE, biasing the season expansion. The fix: `h = number_of_gear`, `L = tau` (turnover, a parameter with a lognormal prior), `E = lambda_E * tau`, and `lambda_C` becomes crab per deployment (stable, 4.0-7.6). Boat catch 155,038 -> 43,268; `R_G_boat = 3.554` (3.244-3.874); boat CPUE 2.930 crab/deployment; `beta_h` 0.754 (now covers 1). See Appendix A. This scale is what the pooled track later ported as POOL-1 and POOL-3.
-   **F4, CPUE diagnostics wired.** `write_cpue_diagnostics`, `bss_effort_spec`, and the linearity / saturation / estimator-triad checks plus the effort-unit assertion. The linearity check regresses log-catch on log-effort to estimate `beta_h`; the saturation check bins by effort and reports the per-unit rate gradient; the estimator-triad check compares ratio-of-sums, mean-of-ratios, and the model CPUE. These are the diagnostics that surfaced the shore effort-unit question resolved in v5.5 and are the model-independent evidence that gear-hours is invalid for pots (the F2 finding).
-   **F5, reporting.** `ar_resolution` populated; `divergence_fraction` / `impact` / `pass` columns added to the convergence report; `expansion_ratios` reports `R_G_boat` (replacing the removed `R_T`); `port_total` composition reported so the assembled total is auditable against its components. These are report-only and change no estimate.
-   **P0 (2026-07-10, commit d1968f3), PE population-argument bugfix and ratio-of-sums strata.** `run_pe` received the component label (for example `"shore_all_gear"`) in the slot where it expected the population, so every shore component silently ran the boat branch. Combined with a switch to ratio-of-sums PE strata (stratum CPUE = `sum(catch)/sum(h)`, replacing the unstable weighted mean-of-ratios) and a permanent PE implied-CPUE-vs-ratio-of-sums assertion. The PE had been internally inconsistent: implied CPUE was 2.88x the ratio-of-sums for shore_all_gear, 2.61x for shore_ring_net, and 0.55x for private_boat. Effect: shore PE 45,112 -> 19,767; PE total 80,160 -> 54,598. See Appendix B.
-   **P1 (2026-07-10), configurable effort unit with LOO comparison.** `bss_effort_spec.R`, `shore_effort_unit`, and `loo_effort_unit_comparison`; the Stan `effort_scale_gear` / `E_scale` machinery so `E` always carries `h`'s unit, which makes the cross-unit `elpd_loo` comparison valid. This is the lever the v5.5 shore decision was read off.
-   **P2 (2026-07-10), mu-hierarchy collapse at `G*S == 1`.** `use_mu_hier = (G*S > 1) ? 1 : 0` drops the redundant single-cell hierarchical level that funnels. The pooled track attempted the identical change (B1.7 in v6.9) and reverted it, because its shore all-gear DAILY AR (289 states) hung for more than 24 hours once the decoupled level was removed. The gear-resolved track is safe because its shore AR is monthly or biweekly, not daily, so removing the level does not force the effort level to reconcile against a long daily AR. P2 produced two new backlog items: GR-17 (the collapsed `sigma_mu` posteriors look non-converged but are inert) and GR-10 (the collapse is hardcoded and needs a config lever). The pooled track carries the equivalent lever as POOL-4's `collapse_mu_hier` (default off, so its production posterior is unchanged); the gear-resolved track applies the collapse unconditionally at `G*S == 1` and still needs the same config treatment.
-   **Parity ports already in the Stan file.** B1.3 non-centered AR(1) initial state; B1.5 effort overdispersion marginalized to NB2; B1.6 `sigma_IE` proper unconditional prior; B1.8 `C_expected_sum` deterministic gate total; B1.9 weekend CPUE effect `B1_C` (`estimate_B1_C`, default TRUE). These carry the pooled convergence-work conclusions into the gear-resolved model without re-deriving them.

### Outstanding backlog (GR-7 to GR-17)

The GR-series is maintained in `07_documentation/development_notes/20260710-OUTSTANDING_ISSUES.md`, which carries a per-item "Status after the 2026-07-10 run" block. An important caveat when reading the status column: the 2026-07-10 run's only code changes were P0, P1, and P2, so where a GR item's symptom cleared, the movement is usually a side-effect of P0 or P2 rather than a targeted GR fix. GR-8 is the clearest case: its funnel symptom cleared because P2 collapsed a separate `sigma_mu` level and pulled total divergences under the gate backstop, while the item's own proposed fix (the `ie_min_obs` guard and prior retune) never landed.

| ID | Item | Status |
| --- | --- | --- |
| GR-7 | Option A / genuine per-gear CPUE unbuilt: `G = 1`, the per-gear machinery is inert, gear catch is PE-apportioned rather than modeled | OPEN |
| GR-8 | `sigma_IE` funnel in shore_ring_net | PARTIAL. Symptom cleared by the P2 side-effect (divergences 598 -> 216, fraction 7.47% -> 2.7%), but the actual `ie_min_obs` guard / prior retune never landed; `sigma_IE` is still the #1 divergence driver (SMD -0.444) |
| GR-9 | shore_all_gear `sigma_IE = 0.603` unexplained | OPEN |
| GR-10 | No `collapse_mu_hier` config lever (`use_mu_hier` is hardcoded) | OPEN |
| GR-11 | Gate policy: the hard 5% divergence backstop overrides the scale-aware impact test | OPEN |
| GR-12 | `tau_boat = 1.2` rests on only 2 WBL I/E days | OPEN |
| GR-13 | Dead `stan_data` entries plus dot-prefixed metadata passed into `rstan::stan(data=)` | OPEN (low) |
| GR-14 | Latent AR carries only 4-13% of effort variance; NB overdispersion carries 81-91% | TRUE BY DESIGN (coarse fixed AR) |
| GR-15 | `private_boat_ring_net` never fits (about 17 interviews) | OPEN (low) |
| GR-16 | Shore effort unit | RESOLVED 2026-07-11 (see v5.5) |
| GR-17 | Collapsed `sigma_mu` posteriors look non-converged but are inert | NEW / OPEN (low) |

------------------------------------------------------------------------

## Appendix A: F1 and F2 working notes, the boat improper prior and the deployment-scale correction

These are the consolidated working notes behind fix-markers F1 and F2, the two changes that made the shore fits converge and corrected the boat effort scale. Together they are the single largest quantitative change in the gear-resolved model's history, and F2 is the reference implementation the pooled track later ported as POOL-1 and POOL-3.

**F1, the improper direction. What went wrong.** `R_G_boat` (gear per boat group) is a boat-only parameter, used only in the boat trailer-expansion likelihood. It was declared `real<lower=0>` with no prior. In a shore fit it appears in no likelihood at all. A `<lower=0>` declaration means Stan samples an unconstrained `z` and maps `R_G_boat = exp(z)`, adding the Jacobian `log|dR/dz| = z` to the target. With no prior and no likelihood, the target along that direction is just the Jacobian `+z`, which is unbounded above, so the sampler pushed `z` and hence `R_G_boat` toward positive infinity. `R_G_boat` ran to about 1e307 (the floating-point ceiling), the runaway direction wrecked mass-matrix adaptation, and 97.6% of shore transitions diverged. This is structurally the same bug as the pooled B1.6 `sigma_IE` improper direction (a constrained parameter with a Jacobian that no prior or likelihood ever touches), which is why it is filed as B1.10 in the shared marker sequence.

**F1, the fix and outcome.** Give `R_G_boat` a proper unconditional prior `R_G_boat ~ lognormal(log(4), 0.5)`, applied in every fit rather than only when the population is boat. In a shore fit it now enters no likelihood but is proper, so it sits near 4 and cannot run away; in a boat fit the likelihood updates it as before. The effect on the shore fits was decisive: shore_ring_net divergences 7818/8000 (97.7%) -> 598/8000 (7.5%), shore_all_gear 7807/8000 (97.6%) -> 213/8000 (2.7%), with R-hat 1.0002-1.0026 and n_eff 2937-9313. The shore components became reportable BSS fits rather than near-total divergence failures forced to PE.

**F2, the wrong unit. Why gear-hours is invalid for pots.** With the shore fits converging, the remaining error was the boat effort unit. The boat CPUE denominator was gear-hours (`h = gear_time_total`, `L = 24`), which assumes boat catch is linear in soak time. It is not. Binned by hours-per-gear, crab per gear-HOUR falls 43x (from 2.88 to 0.068) across a 64x soak range, while crab per gear per trip rises only 1.8x. A log-log fit over 1,532 interviews gives `crab_per_gear ~ h^0.133`, an exponent near zero. The mechanism is saturation: a pot or trap catches most of what it will catch soon after it is set, and the per-hour rate collapses as it soaks. Gear-hours forces a 43x per-hour gradient onto a single constant CPUE, so a fit that matches short soaks over-predicts long ones and the constant-CPUE season expansion inherits the bias.

**F2, the fix and outcome.** Measure boat effort in gear DEPLOYMENTS, not gear-hours: `h = number_of_gear` (the count of gear set), `L = tau` (turnover, a parameter with a lognormal prior, roughly how many times gear is pulled and reset per day), `E = lambda_E * tau`, and `lambda_C` becomes crab per deployment, which is stable across the soak range (4.0-7.6). The reported quantities moved substantially: boat catch 155,038 -> 43,268 (about a 3.6x reduction), `R_G_boat = 3.554` (3.244-3.874), boat CPUE 2.930 crab/deployment, and the linearity slope `beta_h = 0.754` now covers 1 where gear-hours had held it far from 1. This is the number that supplied the pooled track its missing validation and moved the pooled boat from gear-hours (about 56,266) onto deployments (about 43,314). The two fixes reinforce each other: F1 made the shore fits reportable and F2 put the boat on the harvest-unbiased scale, so both fitted populations now sit on BSS posteriors that the gate can accept.

**The standing lesson.** F1 is the same class of bug as the pooled B1.6 `sigma_IE` improper direction, and it was caught the same way, by a divergence-localization diagnostic that pointed at a single parameter pinned at the floating-point ceiling. That both tracks independently grew an improper boat-only direction in a shore fit is the argument for the shared parity ports and the shared effort module: a fix derived once should be carried, not re-discovered. F2 is the mirror-image lesson on the R side. A denominator that looks like a reasonable effort measure (gear-hours) can be quantitatively wrong by a factor of three on the harvest total, and the only reliable test is an empirical one (the saturation and linearity diagnostics of F4), not a plausibility argument. Neither fix was adopted on reasoning alone; both were confirmed by a run and by the before/after numbers recorded above.

------------------------------------------------------------------------

## Appendix B: P0 working notes, the PE population bug and the ratio-of-sums estimator

These are the consolidated working notes behind fix-marker P0 (2026-07-10, commit d1968f3). P0 is a point-estimator fix, not a Stan change, but it matters because the PE is both the fallback estimate when a BSS fit is gated out and the independent cross-check on every BSS fit. A badly biased PE corrupts the gate comparison and, when a component falls back, ships a wrong number.

**The population-argument bug.** `run_pe()` is called once per component. It received the component label (for example `"shore_all_gear"`) in the argument slot where it expected the population (`"shore"` or `"private_boat"`). Downstream, the branch that selects the shore-versus-boat effort and CPUE construction tested that argument against the boat population, and since `"shore_all_gear"` is not `"private_boat"` the test never matched for shore, so every shore component silently fell through to the boat branch and was computed with the boat's effort construction and denominator. The symptom was not a crash; it was a plausible-looking but wrong shore PE. The fix passes the population explicitly, separately from the label.

**The estimator inconsistency.** Independently, the PE stratum CPUE was a weighted mean-of-ratios, the average of per-interview `catch/effort` ratios. That estimator is unstable when effort varies across interviews within a stratum, because a low-effort interview with a few crab produces a large ratio that the mean over-weights relative to its share of the catch. The result was an internally inconsistent PE: the implied CPUE (PE catch divided by PE effort) was 2.88x the ratio-of-sums for shore_all_gear, 2.61x for shore_ring_net, and 0.55x for private_boat. The fix replaces mean-of-ratios with ratio-of-sums (stratum CPUE = `sum(catch)/sum(h)`), the harvest-consistent estimator whose implied CPUE equals itself by construction, and adds a permanent assertion that the PE implied CPUE equals the ratio-of-sums within tolerance, so this class of inconsistency cannot silently return.

**Effect and significance.** Together the two changes dropped the shore PE from 45,112 to 19,767 and the PE total from 80,160 to 54,598. Because the corrected PE is the cross-check the BSS is read against, P0 is also what made the shore effort-unit question legible: once the PE and the ratio-of-sums agreed, the linearity and estimator-triad diagnostics (F4) could show cleanly that no shore unit except gear-deployments produced a consistent, unbiased rate, which is the finding that GR-16 and framework v5.5 later resolved. P0 restored the PE as a trustworthy fallback and a trustworthy cross-check in one commit.

**Why the assertion is permanent.** The population-argument bug was silent: it produced a plausible number, not an error, so it survived until the ratio-of-sums cross-check exposed a 2.88x discrepancy that no correct estimator could produce. The permanent implied-CPUE-vs-ratio-of-sums assertion is the guard against that class of failure returning under a later refactor. It halts the run rather than shipping a PE whose implied CPUE and stratum CPUE disagree, which is the signature of either a mis-routed population or a re-introduced mean-of-ratios. This is the point-estimator analogue of the effort-unit assertion added in F4: a cheap invariant that converts a silent, plausible-looking bias into a loud, immediate failure.

------------------------------------------------------------------------

## Appendix C: Numbers of record

The two runs that bracket the F and P series. These are the component-level numbers the totals were assembled from; they are preserved here because the fit objects are not saved and the run outputs are the only durable record.

**20260709 run (commit 564071b).** Before the P0 PE fix and before the shore effort-unit decision.

| Component | PE | BSS | Reported basis |
| --- | --- | --- | --- |
| shore_ring_net | 12,174 | 6,930 (pooled BSS 6,733) | gated PE |
| shore_all_gear | 32,937 | 19,680 | gated BSS |
| private_boat_all | 22,070 | 43,268 | gated BSS |
| comm / charter | - | - | 12,007 census |

Reported port total 88,238 (95% CI 75,656-106,266), using the PE for shore_ring_net. The BSS-reconstituted interim total (substituting the shore ring-net BSS for its inflated PE) was about 82,856.

**20260710 run (commit d1968f3).** After P0, F1, and F2.

-   Convergence improved from 2/3 to 3/3 fits passing.
-   Shore ring-net: the inflated PE 12,174 was replaced by the BSS 6,936 (the fit now passes the gate).
-   Shore all-gear: BSS about 19,684.
-   Boat all-gear: BSS about 43,314, on validated deployments (the F2 scale).
-   PE totals dropped from the P0 fix: shore PE 45,112 -> 19,767; PE total 80,160 -> 54,598.

The 20260710 run is the first in which all three fitted components report on their BSS posteriors on the corrected effort scales, with the PE serving as a consistent cross-check rather than an inflated fallback. Framework v5.5 (2026-07-11) then moved the shore component onto gear-deployments to match the boat and closed GR-16.

------------------------------------------------------------------------

## Current production configuration

For a reader running the model today, the production state as of framework v5.5 (2026-07-11) is:

-   **Effort units.** Shore and boat both on gear-deployments (`h = number_of_gear`, `L = tau`); `loo_effort_unit_comparison = FALSE`. Both units have a linearity `beta_h` covering 1, so the season expansion is harvest-unbiased.
-   **Fitted components.** Shore ring-net, shore all-gear, and private boat all-gear are fit as BSS and, in the latest run, all pass the convergence gate; commercial and charter is a census tally. `private_boat_ring_net` never fits and falls back to PE (GR-15).
-   **Gates and filters.** R-hat < 1.01, n_eff > 400, a divergence-fraction backstop at 5%, and the scale-aware impact test (GR-11 notes the backstop can still override the impact test). `filter_incomplete_trips = TRUE`. `filter` and gate helpers are the shared modules, so pooled and gear-resolved cannot drift.
-   **Known open dependencies on the reported number.** The boat catch rests on the `tau_boat = 1.2` prior until the boat I/E stream activates (GR-12), and per-gear catch is PE-apportioned rather than modeled because `G = 1` (GR-7). State both in any co-management number.
-   **The two header lags** (method-doc header `v5.4`, RMD banner `Framework v5.4`) are cosmetic and do not reflect the running code, which is v5.5.
