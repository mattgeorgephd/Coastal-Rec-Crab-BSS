# Coastal Rec Crab BSS: Pipeline Status and Backlog

- **Last updated:** 2026-07-14
- **Maintainer note:** this is the single living status document for the pipeline. It replaces the six scattered development notes listed in Section 8, reconciling their issue IDs so nothing is lost. Update this file as work lands; do not re-fork it into per-session notes.

**Repo:** `Coastal-Rec-Crab-BSS`, `main`. **Method of record:** Method v1.0 (frozen against pooled code v7.4); the code has advanced past v7.9 with the Tier-2 batch (nine items, 2026-07-13), and the authoritative run is now **Run 1** (`20260713/pooled-CPUE-run1`). The 2024-25 reference numbers in the method documents remain pre-refresh; they should be regenerated to the Run 1 totals (port total 83,035; Section 1).

---

## 1. Current state at a glance

| Component | Version | State |
|---|---|---|
| Pooled CPUE model (`crab_bss_pooled.stan`, `BSS-GH-pooled-CPUE-model.Rmd`) | pipeline code **v7.9 + Tier-2 batch** (2026-07-13) | Production, validated by Run 1 (`20260713`). Both fitted components on gear-deployments; boat on monthly AR by the resolution map (`ar_force = NULL`). |
| Gear-resolved model (`crab_bss_gear_resolved.stan`, `BSS-GH-gear-type-CPUE-model.Rmd`) | framework **v5.6** (2026-07-12) | Production, validated by Run 2 (`20260713/gear-type-CPUE-model-run2`, 3/3 fits pass). Shore + boat on gear-deployments; `G = 1` (per-gear CPUE inert, gear split from PE apportionment, see GR-7). |
| Weather/tide covariate module (`crab_bss_pooled_weather_adjusted.stan`, `06_diagnostics/...Rmd`) | ~v6.9 parity, **stale** | Not production. Forked engine, pre-deployment-scale Stan, not re-run since the day-length extraction. Do not cite its boat number. |
| Config surface (`run_config.R`) | base-params architecture (P5) | `run_config` is the base parameter set; each driver layers model-specific tuning via `modifyList(run_config, params_model)`. |

**Authoritative run:** `05_output/20260713/pooled-CPUE-run1` (**Run 1**, pooled default config, 2026-07-13). Boat on monthly AR selected by the per-population resolution map (`ar_force = NULL`); gate at fraction 0.05 / impact 0.10 SD. Run 1 folds the nine-item Tier-2 + parity batch (below) on top of the 7/12 Tier-1 baseline and **validates** it: the two passing BSS fits are essentially identical to 7/12 (shore all-gear 20,655 -> 20,607; boat all-gear 43,480 -> 43,180), so the nine items are behavior-neutral on the reported total (port total 83,825 -> 83,035, -0.9%). The one behavior change is the boat point-estimator (empty-stratum fallback; see the headline), a cross-check gain rather than a change to the reported number. The prior authoritative run was `20260712/pooled-CPUE-morning`.

**Update 2026-07-13 (Tier-2 + parity batch; nine items landed after the 7/12 run, across two applied patches):**

- (2) PE empty-stratum CPUE fallback (`pe_empty_stratum = "pooled"`, both PE tracks).
- (3) Dirichlet `pi_gear` x BSS-total-draws intervals for the pooled gear breakdown.
- (4) shore-I/E representativeness diagnostic (I/E days are representative, so the shore all-gear `sigma_IE` ~1.07 is sparse-data scatter, not peak-day bias).
- (5) `ie_min_obs_shore >= 3` guard that drops the I/E stream and decouples `sigma_IE` for the thin shore pot-closure / ring-net fits (GR-8; prior deliberately left `exponential(5)` per item 4).
- (7) `max_divergences` doc reconciled to 5.
- (8) stale gear-resolved crabber-hours labels fixed to gear-deployments.
- (9) PE monthly-share consolidated into `pe_monthly_effort_share.R` (also fixes the gear boat PE-fallback to use the deployment formula), residual-globals verified absent.
- (6) **always-on** holiday CPUE term `B2_C` plus an optional same-day-effort density term `gamma_C` (`estimate_cpue_density`, off by default).
- (1) razor-dig shore-effort term `B3` (`razor_dig_mode = "no"|"yes"|"auto"`, off by default) with the two opener workbooks consolidated into `04_input_files/fishery_opener_dates.csv`.

The Stan model gained `B2_C` / `B3` / `gamma_C`, all confirmed present and sampling in Run 1. **Run 1 has landed and validates the batch** (2026-07-13; the authoritative run above). What actually moved: only the boat empty-stratum PE (boat PE +65%). `B2_C` came out indistinguishable from zero for both reported fits, so it did not move the total; `B3` and `gamma_C` correctly sampled their priors (inert when off). One regression surfaced: the thin shore pot-closure fit lost convergence and now falls back to PE (see the convergence note and Tier 2). **Validation runs:** Run 1 (pooled) and **Run 2 (gear-resolved) both landed 2026-07-13 and cross-validate** (Section 1, fact 1); the gear-resolved track reports 3/3 fits and its ring-net fit passes at biweekly AR where the pooled one failed at daily. **Run 3 (razor-dig `"yes"`, pooled) landed 2026-07-14:** activating the shore-effort razor term `B3` buys no predictive gain (elpd_loo within 1 SE of B3-off on every stream) and moves the port total only +0.6%, so it stays off by default; its one surprise is that it stabilizes the thin ring-net fit (Tier 2). **Remaining experiments:** density (`gamma_C`) and the R_G prior sweep.

**Headline numbers, Run 1 (2026-07-13 pooled default, Dungeness kept 2024-25):**

Reported port total (gated): BSS median **83,035** (95% CI 70,527 - 102,302; predictive 83,004, CI 70,461 - 102,438). Effort BSS median 44,156 (CI 40,194 - 49,336). Essentially unchanged from the 7/12 run (83,825), which confirms the estimate.

| Fitted component | BSS catch (median) | Gate | 7/12 BSS |
|---|---:|---|---:|
| Shore all-gear (Dec 1 - Sep 15) | 20,607 | BSS | 20,655 |
| Private boat all-gear | 43,180 | BSS | 43,480 |
| Shore pot-closure (Sep 16 - Nov 30) | 6,215 | **PE (convergence fail)** | 6,587 (was BSS) |
| Private boat pot-closure | n/a | PE (no data) | n/a |
| Commercial/charter | 11,986 | Census | 11,986 |

**Boat PE reconciliation (the one number that moved).** The empty-stratum pooled fallback (item 2) lifted the boat mode PE from 22,823 to **37,638** (+65%), raising implied boat CPUE 1.83 -> 3.02 (interview ratio-of-sums 3.36) and narrowing the boat PE-vs-BSS gap from 1.9x to 1.15x. The reported boat is still the BSS (43,180) in both runs, so this is a cross-check confidence gain, not a change to the published number. PE by mode: shore 21,533 (unchanged), boat 37,638, commercial/charter 11,986; PE port total 71,157 (was 56,342).

**Convergence (Run 1, 2026-07-13):** the two reported fits pass cleanly (shore all-gear daily 3.0% divergences, impact 0.001 SD; boat all-gear monthly 3.1%, impact 0.020 SD; R-hat 1.00, n_eff > 5,000). **Regression:** shore pot-closure (ring-net) went from 150 divergences / PASS on 7/12 to **1,165 divergences (19.4% fraction) / FAIL** on Run 1 (impact 0.11 SD, over the 0.10 gate), so it now reports PE. The most likely cause is the always-on `B2_C` destabilizing an already-marginal thin fit running at daily AR (the pot-closure sub-season carries holidays); the gate caught it and the fallback PE is defensible (its `beta_h` = 0.92 covers 1). Boat pot-closure remains PE (no data). See Tier 2 for the fix.

**Gear-resolved Run 2 (2026-07-13, independent cross-check).** The gear-resolved default run reports a port-total BSS of **82,995** (95% CI 70,238 - 101,111), within 0.05% of pooled Run 1 (83,035), and all 3 fits pass. Its ring-net (shore pot-closure) fit passes cleanly at **biweekly** AR (22 divergences, 0.27%) where the pooled ring-net failed at daily AR, which pins the pooled regression on AR resolution, not `B2_C` (Tier 2). Diagnostics match the pooled run (all `beta_h` cover 1; PIT means ~0.50; boat 79.9% extrapolated). Caveat: the gear-resolved gear breakdown is degenerate under `G = 1` (`catch_by_gear_type.csv` lumps 84% of catch into "All"); the pooled Dirichlet split is the better gear breakdown (GR-7).

**The three facts that most shape how to read the numbers:**
1. **The boat is cross-model consistent (RE-CONFIRMED on the Run 1 / Run 2 pair).** The two independent default-config pipelines now agree on the port total to 0.05% (pooled 83,035 vs gear-resolved 82,995) and on the boat all-gear BSS to 0.3% (pooled 43,180 vs gear 43,314); boat effort matches to ~1% (14,707 vs 14,805) and shore all-gear to ~2% (20,607 vs 20,134). This is the strongest internal cross-check available. The earlier ~26% pooled-vs-gear boat gap was entirely the weekly-vs-monthly AR resolution, not the effort unit; both tracks now run boat at monthly.
2. **The BSS sits above the PE, but the boat PE gap is now mostly closed (Run 1).** On 7/12, ratio-of-sums (P0) left the boat PE at 21,910 (implied CPUE 1.84), far below the interview ratio-of-sums (3.36) and the BSS latent CPUE (~2.95): P0 fixes the WITHIN-stratum estimator (it made the SHORE PE consistent, implied 0.865 vs RoS 0.873) but not the ACROSS-stratum effort weighting. Run 1's empty-stratum pooled fallback (item 2) targets exactly that, stopping the zeroing of effort-bearing but interview-empty boat strata: the boat PE rose to 37,638 (implied CPUE 3.02), closing ~75% of the gap to the interview ratio-of-sums, and the boat PE-vs-BSS ratio fell from 1.9x to 1.15x. The reported boat remains the BSS (43,180), cross-checked against the gear-resolved boat (Run 2 pending); the residual ~15% is genuine effort/CPUE imputation over the 80%-extrapolated boat series (fact 3), not a PE estimator artifact.
3. **The harvest is summer-dominated and lean on extrapolation, now quantified (T1.1a).** The boat all-gear BSS catch (the single largest component, ~43k) is **79.8% extrapolated** on Run 1: only 66 of 289 all-gear days (23%) carry a boat interview, and the AR process supplies the other 80% of the catch. Shore all-gear is 48.5% extrapolated (145/289 sampled), shore pot-closure 30.9% (51/76). This is the largest standing threat to the headline number's credibility, and it is the residual behind the boat's 15% PE-vs-BSS gap in fact 2. Report the boat number with this caveat attached. Run 2 reproduces the profile on the gear track (boat 79.9%, shore all-gear 52.5%, shore pot-closure 31.8%).

---

## 2. Repository map

- **`01_BSS_models/`**, the two production driver `.Rmd` (pooled, gear-resolved) and their rendered `.html`. Pooled is v7.9, gear is v5.6. The `-old.Rmd` snapshots were removed 2026-07-12.
- **`02_stan_models/`**, `crab_bss_pooled.stan` (v7.6+; `R_G_boat`, `effort_scale_gear`/`E_scale`, `collapse_mu_hier` lever), `crab_bss_gear_resolved.stan` (has the non-centered `omega_0`, see D4), `crab_bss_pooled_weather_adjusted.stan` (stale).
- **`03_R_functions/`**, 25 shared modules. Post-refactor: extracted drivers (`fetch_crab_data*`, `run_pe_pooled`/`run_pe_gear`, `prep_bss_crab_pooled`/`_gear`, `prep_days_crab`, `prep_population_summary`, `estimate_comm_charter`, `classify_day_type`, `bss_timers`) plus the diagnostic/engine modules (`bss_ar_resolution`, `bss_convergence_gate`, `bss_cpue_diagnostics`, `bss_effort_spec`, `bss_trailer_expansion`, `bss_day_length`, `diagnose_effort_overdispersion`, `divergence_diagnostic`, `model_diagnostics`, `save_run_diagnostics`).
- **`04_input_files/`**, `effort_combined.csv`, `interview_combined.csv`, `wes_commercial_tally.csv`, `ingress_egress.xlsx`, plus (2026-07-13) `MA2-fishing-dates-2023-2026.xlsx` and `razor-clam-dig-dates-2021-2025.xlsx` for the fishery-opener spillover diagnostic.
- **`05_output/`**, dated run folders. Committing full run outputs is a known hygiene problem (T4.2/ORCH-26) and keeps growing with each run. The current authoritative run is `20260713/pooled-CPUE-run1` (Run 1).
- **`06_diagnostics/`**, the weather/tide covariate module (stale).
- **`07_documentation/`**, method docs, development histories, READMEs, and `development_notes/` (this file plus the historical notes).
- **Root**, `run_config.R` (control surface), `run_estimation.R` (orchestrator), `README.md`, `README-R-functions.md`.

---

## 3. Completed work (reconciled across all prior notes)

Historical IDs are preserved in parentheses so the older notes remain traceable. Legend of ID schemes: **B1.x/A/B/C/D** from `PLANNED_IMPROVEMENTS.md`; **T1.x-T4.x** from `CODE_IMPROVEMENTS_REVIEW_v7.0.md`; **O1-O12** from `ADDITIONAL_OUTPUTS_PROPOSAL.md`; **F1-F5/P0-P3** from `pipeline_state_review_20260709.md`; **POOL-*/GR-*/WX-*/ORCH-*** from `20260710-OUTSTANDING_ISSUES.md`; **critique 1-11** from the 2026-03-31 model critique; **P4/P5/P7** from the 2026-07-12 review.

**Convergence and sampler geometry.**

- BSS convergence achieved (B1, B1.1-B1.6). Boat all-gear diverged on ~100% of iterations under daily AR (a funnel from a 289-state latent on a thin trailer series); fixed by the per-population AR cap to weekly (B1.2/v6.5), non-centering the AR initial state `omega_0` (B1.3/v6.6), and boat sampler tuning (v6.2). Effort overdispersion marginalized to negative binomial (B1.5/v6.7). Unconditional `sigma_IE` prior (B1.6/v6.8).
- Non-centered `omega_0` ported to the gear-resolved Stan (**D4 DONE**; `omega_E_0_raw`/`omega_C_0_raw` present in `crab_bss_gear_resolved.stan`).
- `R_G_boat` improper-prior blocker fixed in gear-resolved (F1): an unbounded `real<lower=0>` with no prior/likelihood drove the log posterior to +inf; now `lognormal(log(4), 0.5)`. Ported to pooled as POOL-1.
- Scale-aware convergence gate (B1.8/v7.0): pass/fail on the divergent draws' impact on the reported totals in posterior-SD units, not a raw count. The 15% fraction backstop was tightened to **0.05** (P7, 2026-07-12) to match the gear-resolved gate; the impact test (threshold 0.10 SD) remains primary.

**Effort unit and scale (the dominant correctness theme).**

- Boat moved off invalid time-denominated effort onto gear-deployments (F2 gear-resolved; POOL-1 + POOL-3 pooled/v7.6). The saturation diagnostic proved catch is sub-linear in soak time for pots (`beta` ~0.13-0.27), so any time unit is invalid; `R_T` (pinned at ~1) was replaced by `R_G_boat`, and `L = 24` gear-hours by `L = tau_boat`.
- Shore moved onto gear-deployments (GR-16 / POOL-7 / v7.7), chosen by a three-unit LOO comparison: gear-deployments is the only shore unit whose catch-vs-effort elasticity covers 1 (`beta_h = 1.05` vs 0.57 crabber-hours, 0.73 gear-hours). Applied to both pipelines. This is the resolution of critique-1's day-length/snapshot concern for the CPUE denominator (the level still uses a gear-count snapshot x turnover; see T1.1).
- Shore PE aligned to the same unit (v7.8): `run_pe_pooled` now takes the shore denominator from `bss_effort_spec.R`; the shore PE effort dropped from 42,541 (crabber-hours) to 18,104 (deployments), removing a PE-vs-BSS unit inconsistency.
- **Boat AR resolution resolved (BOAT_RESOLUTION_EXPERIMENT closed).** The experiment asked daily-vs-weekly; the answer is **monthly**. On monthly AR the pooled boat reconciles to the gear-resolved boat (Section 1, fact 1), and the sparse-month effort over-imputation collapses (June boat BSS/PE fell 8.7x to 4.1x; July 3.7x to 2.1x).

**CPUE structure.**

- Weekend/holiday CPUE effect `B1_C` added and confirmed to keep on (critique 4 RESOLVED, T3.3 base). Shore CPUE ~46% lower on weekend/holiday (`B1_C ~ -0.6` shore all-gear), indistinguishable from zero for boat (soaking gear does not care about day type). Physically sensible crowding signal.

**Filters.**

- Incomplete-trip filter added to the pooled track (POOL-2 / B7 / T2.3), matching gear-resolved; ~40% of shore all-gear interviews are incomplete with a measured ~-20% CPUE bias, so filtering raises the shore estimate. Toggle `filter_incomplete_trips` (default TRUE), with a `sensitivity_incomplete_trips.csv` output (now rendered in the report).

**Diagnostics and persisted outputs.**

- The output catalog O1-O11 is implemented (`save_run_diagnostics.R`): full parameter summaries, AR latent path, period coverage map, modeled daily CPUE, per-observation PPC residuals, sampler diagnostics (E-BFMI), summed-quantity draws, prior-vs-posterior, gear proportions, monthly PE-vs-BSS by mode, per-fit data summary. O12 (`.rds`) deliberately not committed.
- CPUE effort-unit diagnostics wired (POOL-5): `cpue_estimators_/cpue_saturation_/cpue_linearity_*.csv`, the estimator triad + saturation + linearity that surfaces an invalid effort unit automatically each run.
- Posterior predictive checks run and calibrated (B2 UNBLOCKED and reported): 95% coverage ~0.95-0.98, PIT means ~0.5; the 50% coverage is wide for shore catch / boat trailer (consistent with `r_C ~ 0.75`), i.e. the effort/catch predictive is mildly over-dispersed.
- Effort overdispersion decomposition (T1.5): NB observation overdispersion dominates (81-91%); the latent AR carries 4-13%, a consequence of the coarse AR, not a bug.
- **Report display pass (2026-07-12):** the pooled `.Rmd` now renders the incomplete-trip sensitivity, PPC calibration, effort-overdispersion decomposition, the CPUE validity triad, per-fit data coverage, PSIS-LOO, and top divergence drivers as on-page tables; the wide convergence table is curated. Sections that previously only wrote CSVs now show results.

**Refactor and config.**

- Behavior-preserving function extraction into `03_R_functions/` and centralization of all user toggles in `run_config.R` (v7.8). Config restructured so `run_config` is the base parameter set and each model layers its tuning (P5/v7.9); the per-model AR resolution map lives in `run_config.R`.
- **Pot-closure sub-season made explicit and renamed (2026-07-13).** The pot-closure window (when pots are illegal and only non-pot gear is legal) is now set explicitly in `run_config.R` via `pot_closure_start` / `pot_closure_end`, instead of being assumed to run from the season start to `pot_open_date - 1`, so a future season whose start does not coincide with the closure is supported. A shared builder `03_R_functions/build_subseasons.R` derives the sub-seasons for both drivers and handles the general case (optional pre/post all-gear periods for a mid-season closure, keyed off `gear_regime` rather than the sub-season name). The sub-season was renamed "ring-net only" to "pot closure": a **display-only** rename, since non-pot gear other than ring nets is legal, so the old name was a misnomer; the internal key stays `ring_net_only` for output-filename continuity, so keys and filenames are unchanged and the 2024-25 config produces bit-identical sub-seasons. Season plots in both drivers now draw vertical lines at the closure start and the pots-open date. (The report's downstream aggregation was initially left hardcoded to the two historical sub-seasons behind a fail-fast guard; that has since been generalized, see the next bullet.)
- **Aggregation generalized over the sub-season set (2026-07-13).** The PE port summary, catch-by-mode, final harvest table, and gear-allocation totals in both drivers now sum over whatever sub-seasons `build_subseasons` produced, via semantic helpers (`pe_pop_sum` / `pe_port_total` / `pe_field`) instead of hardcoded `pe_all$..._ring_net_only + pe_all$..._all_gear` references and fixed `pe_summary` row indices (`[6]`, `[1:5]`, `[3]+[4]`). Verified numerically identical to the prior hardcoded code for the 2024-25 two-sub-season config (independently reviewed; `pe_port_total` uses `na.rm=TRUE` to match the baseline exactly), and the fail-fast guard is downgraded to a soft note, so a mid-season pot closure now aggregates end-to-end. Confirmed behavior-neutral on the 7/12 run.
- **Fishery-opener spillover diagnostic added (2026-07-13; pooled report Section 3.5).** New `prep_fishery_events.R` + `diagnose_fishery_spillover.R` read the MA2 finfish and razor-clam-dig calendars and test whether crab effort/CPUE differs on opener days, raw and adjusted for day-type + month, across shore/boat effort and shore/boat CPUE. DIAGNOSTIC ONLY (no model change). See Tier 2 for the 7/12 result and the pending decision.

**PE (partial).**

- Shore PE unit alignment done (v7.8) and ratio-of-sums (P0) done + confirmed 2026-07-12: the shore PE is now internally consistent (implied CPUE 0.865 vs RoS 0.873). The boat PE remains structurally low (its effort expansion weights low-CPUE strata; implied CPUE 1.84 vs interview RoS 3.36), which P0 does not fix, so the boat is reported from the BSS. See Section 4 / fact 2.

**Repo cleanup (partial).**

- `-old` driver snapshots and the refactor tarball removed (2026-07-12).

---

## 4. Outstanding work (prioritized by publication impact)

### Tier 1: can move the headline number, or a reviewer will demand it

Status note: **P0, T1.1a, and T1.4 are CONFIRMED (7/12) and re-confirmed on Run 1 (7/13).** T1.3's prior-sensitivity harness is in place but the sweep runs are pending; external validation (T1.1b) still needs benchmark data. Run 1 also closed two of the 7/12 findings: the boat PE-vs-BSS gap (empty-stratum pooled fallback, item 2) and the boat incomplete-trip-filter anomaly (now resolved; Tier 2).

- **P0 [DONE + CONFIRMED 2026-07-12]: pooled PE on ratio-of-sums.** `run_pe_pooled` computes stratum CPUE as within-stratum `sum(catch)/sum(hrs)`, replacing the unstable `weighted.mean(daily_ratios, w = n_int)`, with an implied-CPUE-vs-ratio-of-sums guard that fails fast before the multi-hour fits. Mirrors `run_pe_gear`. **Run result:** the SHORE PE is now internally consistent (implied CPUE 0.865 vs interview ratio-of-sums 0.873; `cpue_estimators_*` drift flag FALSE). The BOAT PE, however, stayed at 21,910 (implied CPUE 1.84), NOT the ~40k an earlier version of this note predicted: ratio-of-sums fixes the within-stratum estimator but not the across-stratum effort weighting, and the boat's expanded effort concentrates on low-CPUE strata, so its effort-weighted CPUE remains far below the interview ratio-of-sums (3.36). The guard passes only marginally for the boat (0.55x vs the 0.5 floor). Net: P0 closed the shore inconsistency and is correct to keep, but the boat PE-vs-BSS gap is a real imputation difference, not an estimator artifact P0 could remove. Renders in Section 4.4. **Update (Run 1):** the across-stratum weighting P0 could not fix is what the empty-stratum pooled fallback (item 2) addresses, and it mostly closes the boat gap: boat PE 21,910 -> 37,638 (implied CPUE 3.02). See Section 1, fact 2.
- **T1.1(a) [DONE + CONFIRMED 2026-07-12]: summer-extrapolation transparency.** `extrapolation_transparency.csv` (RMD Section 11.4) decomposes each fitted component's BSS catch into interview-days vs extrapolated (no-interview) days. **Run result (Run 1):** boat all-gear **79.8% extrapolated** (66/289 days sampled), shore all-gear 48.5% (145/289), shore pot-closure 30.9% (51/76). The summer-dominance caveat is now a number, and it is large for the boat, the biggest component.
- **T1.1(b) [OPEN, needs external data]: external validation (B4).** Still requires a benchmark the repo does not contain: dockside census totals, a known-total year, or the legacy estimator. The strongest internal cross-checks now in place are the two independent BSS pipelines agreeing on the boat (43.5k) and the per-component PE-vs-BSS comparison (fair once P0's run lands); neither is a true external check. Pursue a benchmark before publication.
- **T1.3 [HARNESS DONE 2026-07-12, sweep runs pending]: R_G prior sensitivity.** The R_G prior is now overridable via `params$R_G_prior_mu` / `R_G_prior_sigma` (data-driven when unset), with a commented sweep block in `run_config.R`. Remaining work is the runs: set `R_G_prior_mu` to 1.0, ~1.28 (empirical), and 1.5 in turn and tabulate the port totals, then extend to the `L_effective` prior (`tau_shore`/`tau_boat` are already live toggles). `prior_vs_posterior_*.csv` (O9) shows the prior pull. **Effort: the runs (multi-hour each).**
- **T1.4 [DONE + CONFIRMED 2026-07-12]: commercial/charter vessel split.** `estimate_comm_charter` applies separate per-vessel Dungeness (and red-rock) means for commercial vs charter vessels (from `boat_type_clean`, Guide folded into Charter) to the matching tally columns, with a fallback to the pooled mean for a class with no interviews. **Run result:** the split moved the census from 12,007 (single pooled mean, 7/11) to 11,986 (7/12), a small net change this season because the tally mix happens to sit near the pooled mean, but the mechanism is now correct for seasons where the commercial:charter mix diverges. Applies to both tracks (shared estimator). Renders in Section 4.2. Still pair with a census uncertainty interval (T3.5).

### Tier 2: accuracy, cross-track parity, or credibility

- **[DONE, Run 1] Make monthly the boat AR production default.** Run 1 ran with `ar_force = NULL` and boat monthly selected by the per-population resolution map, so monthly is now the production default by policy, not an override. Confirmed in the Run 1 convergence report (boat `ar_resolution = monthly`).
- **Fishery-opener spillover: decide whether to add an opener day-category (reconfirmed on Run 1).** Adjusted for day-type + month (the valid test), the only crab-specific signal is **nearby razor digs -> shore effort: +18 gear/day, 95% CI [0.6, 35.4], p = 0.045** (the raw comparison was -25%, a Simpson reversal driven by digs falling in low-effort months; adjusting flips the sign). Halibut-open associates with +32 boat trailers (p < 0.001), but the boat response is launch trailer counts, which include halibut boats that are not crabbing, so it is not clean evidence of crab spillover. Salmon's large raw effect (+1,486% boat trailers) vanishes after adjustment (salmon = summer). No opener moved CPUE (catch rate) for either mode. Caveats before acting: single season, 16 tests (~1 false positive expected at alpha 0.05), and this season nearby-dig is collinear with any-dig (Twin Harbors open on all 109 dig days). Recommended: treat the razor-dig-> shore-effort signal as promising-but-not-decisive; if pursued, add a `razor_dig` day-type indicator to the SHORE effort model (analogous to the holiday term) and test its `elpd_loo` gain, rather than adding it blindly. **[TESTED, Run 3, closed] Verdict:** Run 3 set `razor_dig_mode = "yes"`, activating `B3`. The term is weak (shore all-gear 0.14 [-0.14, 0.43], crosses zero; shore pot-closure 0.35 [-0.04, 0.72], positive but the 95% CI just includes zero) and gives **no predictive improvement** (elpd_loo within 1 SE of the B3-off Run 1 on every stream: shore all-gear gear -1258.9 vs -1257.6; catch identical). It moves the port total only +0.6%, mostly because the ring-net fit flips to BSS (next item), not a direct effect. So the razor-dig term does not earn a place in the production estimate; **default stays `razor_dig_mode = "no"`.** The diagnostic signal was real but small, as the p = 0.045 already implied. **Effort: done.**
- **[RESOLVED, Run 1] Boat incomplete-trip filter direction.** On 7/12, dropping incomplete trips LOWERED the boat all-gear PE by 9.1% (24,093 -> 21,910), the wrong direction versus the soak-time rationale. Run 1 shows that was an artifact of the sparse-stratum zeroing the empty-stratum pooled fallback (item 2) fixed: the filter now RAISES the boat all-gear PE by +5.1% (34,699 -> 36,468), the same direction as shore (+9 to +12%). Resolved; no separate data check needed.
- **[DONE, item 3] T2.2 / B5: propagate `pi_gear` uncertainty in the pooled gear breakdown.** Run 1's `catch_by_gear_type.csv` now carries Dirichlet(`pi_gear`) x BSS-total-draws intervals: Pot 59,556 [48,521-76,906], Trap 11,441 [10,065-13,202], Ring Net 10,298 [8,949-12,229], Snare 1,987 [1,717-2,344], Other 410, replacing the too-narrow point-proportion intervals.
- **T2.4 / C1-C2 + WX module: weather covariate evaluation, then integration; rehabilitate the module.** The GAM screen ran for all four cells but the committed PSIS-LOO is only shore ring-net (marginal, ~2.4 SE, below the 4-SE bar). The module (WX-16..20) forks a stale v6.9-parity engine on a pre-deployment-scale Stan and has not been re-run; it must be re-based on the shared modules and current Stan before its LOO or any covariate promotion (C2) is trustworthy. **Effort: medium.**
- **T2.5: durable single-cell mu-hierarchy funnel fix.** With `G*S = 1`, `mu = mu_mu + sigma_mu * eps_mu` is over-parameterized and is the residual funnel the divergences sit on (`sigma_mu_*` is the top divergence driver each run). The `collapse_mu_hier` lever (POOL-4) exists as a safe off-by-default toggle; the naive full collapse (B1.7) was reverted because it forced the level to reconcile against the daily AR ridge. The durable fix is a careful partial reparameterization (single decoupled level absorber), a design task to prototype on shore all-gear, not a mechanical edit. Also make the gear-resolved collapse a config flag (GR-10). **Effort: medium.**
- **Shore `sigma_IE` effort-model tension (GR-9 and pooled analogue).** In the pooled monthly run shore all-gear `sigma_IE` = 1.07 (worse than the gear-resolved 0.60): the I/E day-length/turnover stream and the effort counts disagree by ~e^1 at the tails. Investigate whether the I/E observation days are unrepresentative (peak-count days) before trusting the shore effort posterior to three figures. **Effort: medium.**
- **[DONE, item 5] GR-8: `sigma_IE` guard for thin shore fits.** The `ie_min_obs_shore >= 3` guard is live. In Run 1 it dropped the I/E stream and decoupled `sigma_IE` to its prior for the thin shore pot-closure fit (`sigma_IE` median 0.10, exponential(5)-consistent), while shore all-gear keeps its stream (`sigma_IE` 1.03). The prior was deliberately left `exponential(5)` per item 4. The guard did its job; the pot-closure fit still fails on the AR funnel, a separate issue (new item below).
- **[Run 1 regression, diagnosis settled by Run 2] Shore pot-closure (ring-net) daily-AR funnel.** The thin shore pot-closure fit went from 150 divergences / PASS on 7/12 to 1,165 divergences (19.4% fraction) / FAIL on pooled Run 1, so it reports PE (impact 0.11 SD, over the 0.10 gate). **Run 2 settles the cause:** the gear-resolved track runs the same ring-net fit at biweekly AR and it passes cleanly (22 divergences, 0.27%) with `B2_C` still always-on, so the destabilizer is the daily-AR funnel on a thin series, not `B2_C`. Fix for the pooled track: adopt the gear track's coarser shore AR caps in the pooled `ar_max_resolution` map (biweekly for shore pot-closure; the gear track also caps shore all-gear at monthly, which incidentally calibrates the shore effort predictive better, coverage_50 0.57 vs the pooled daily 0.76). The gate keeps the pooled headline safe in the meantime (PE fallback; `beta_h` = 0.92 covers 1). **Run 3 reveals a second, unexpected fix:** with `B3` (razor-dig) active the same ring-net fit passes at daily AR (171 divergences, down from 1,165) and reports BSS, so the missing razor-dig covariate was itself part of what strained the pot-closure funnel. Two independent fixes now exist, coarser AR (clean, no added parameter, does not perturb the estimate) or `B3` on (physically motivated but adds a weak covariate with no elpd gain); **recommend the coarser-AR fix for production.** The regression is now fully understood, not a mystery. **Effort: low (config change + a confirming re-run).**

### Tier 3: refinements

- **T3.1 / B6 / critique 5: sub-season boundary edge effect.** The hard Dec 1 split fits two independent BSS with no cross-boundary anchor, ballooning the interval at the seam. Either document as a known artifact or move to a single regime-switching model with a `pot_available` indicator. **Effort: medium (doc) to high (model).**
- **[DONE, item 6] T3.3 / critique 4 extension: holiday CPUE term and optional density-dependence.** The holiday CPUE term `B2_C` is now always-on and a same-day-effort density term `gamma_C` is available (`estimate_cpue_density`, off by default). Run 1 result: `B2_C` is indistinguishable from zero for both reported fits (shore -0.18 [-0.58, 0.21], boat -0.19 [-0.85, 0.62]), so there is no holiday-specific CPUE effect beyond the weekend term `B1_C`. `gamma_C` remains available to test density-dependence when wanted.
- **GR-7: genuine per-gear CPUE (Option A).** `G = 1`, so per-gear catch comes from PE apportionment and `catch_by_gear_type.csv` mixes a boat "All" BSS column with PE-apportioned shore gear types (ORCH-24). A real build (effort shares in the offset, a multi-gear rule), not a fix; not to be started until the boat total is stable. Related: critique 7 (priority hierarchy) is resolved in the gear track via weighted pseudo-likelihood; critique 8 (shared `phi_C_gear`) remains a low-priority gear-track item. **Run 2 confirms the current state is not useful:** `catch_by_gear_type.csv` puts 84% of catch in "All" and the Pot/Ring Net/Snare/Trap split is degenerate, so the pooled Dirichlet breakdown (item 3) is materially better and should be the cited gear split until GR-7 is built. **Effort: high.**
- **GR-12: `tau_boat` identification.** Boat catch is proportional to `tau_boat = 1.2`, whose prior rests on 2 WBL I/E days; the boat I/E stream is built and activates once >=2 days fall in a window. The egress-classification pilot is the critical path. State this dependence in any co-management number. **Effort: medium (pilot).**
- **Zero-inflation / ZINB (critique 9 / B2 follow-on): decide on evidence.** The catch likelihood is NB2 with no zero inflation; a high zero fraction (snare, ring-net) may stretch `r_C`. Read the PPC zero-bin and per-observation residuals (O5); prototype ZINB (a variant exists in the freshwater lineage) on the worst-calibrated component only if the zero bin is systematically off, compare by PSIS-LOO. **Effort: low to read; medium if warranted.**
- **T3.2 / critique 11: shore soak-time basis.** Largely mooted now that shore is on gear-deployments rather than crabber-hours, but a one-line confirmation in the docs is worth it. **Effort: trivial.**

### Tier 4: reproducibility and repository hygiene (correctness-neutral, publication-blocking)

- **renv.lock (A4 / T4.1 / ORCH-25).** Still absent. A publication method with a multi-hour Stan pipeline needs a pinned dependency set; ~5 minutes with a large reproducibility payoff. Add a top-level runner too (`run_estimation.R` already serves this).
- **Output and repo hygiene (A1/A2 / T4.2 / ORCH-26).** `05_output/` is 1,017 tracked files / ~82 MB and grows each run, including pre-v6 runs not reproducible from current code; `.Rproj.user/` (5 files) is tracked despite being gitignored. Untrack `.Rproj.user`; decide an `output/` policy (gitignore + one canonical reference run, or Releases). Delete the stale `20260711/pooled-CPUE-morning` and `-afternoon` runs.
- **[DONE, item 7] A3: gear-resolved `max_divergences` doc/code reconciled.** Code and docs now agree on 5.
- **[DONE, item 8] Stale gear-resolved output labels.** The crabber-hours labels (PE plot title, `effort_units_note`) are fixed to gear-deployments.
- **[DONE, item 9] Global-variable coupling (T4.3), duplicated PE monthly-share blocks (T4.4).** The refactor moved functions onto `params`; item 9 consolidated the two PE monthly-share blocks into `pe_monthly_effort_share.R` and verified no residual globals.
- **Surface `B2_C` in the curated report tables.** The always-on holiday CPUE term is written only to `bss_full_summary_*`; it is absent from the curated `structural_params_*` and the on-page convergence/structural tables. Add `B2_C` (and, when active, `B3` / `gamma_C`) to the curated view so the term is visible without opening the full summary. **Effort: trivial.**

---

## 5. Cross-cutting principles and lessons

- **Time-denominated effort is invalid for pot and trap gear.** Catch is sub-linear in soak time (saturation `beta` ~0.13-0.27), so crabber-hours and gear-hours both fail the linearity test; gear-deployments is the harvest-unbiased unit. The pipeline now measures this every run (`cpue_linearity_*.csv`, `cpue_saturation_*.csv`), so any new gear type or population must be run through those before its totals are trusted.
- **AR resolution is a real inference lever for sparse series, not a tuning knob.** Finer AR over-imputes effort in thinly sampled periods; the boat's weekly-vs-monthly choice moved its catch ~20% and was the entire pooled-vs-gear-resolved discrepancy. Cap AR at what the effort series can identify.
- **The harvest is summer-extrapolation-dominated.** The largest uncertainty and the largest BSS-vs-PE divergence both live in the thinly sampled months; no code refinement manufactures information that was not sampled. The single highest-value field change for 2025-26 is more summer interview/effort coverage.
- **BSS-vs-PE gaps are a diagnostic, not a defect, but demand a mechanism.** The two estimate the same harvest by different routes (design-based stratum expansion vs a latent process that imputes unobserved days and propagates uncertainty). A gap is expected where sampling under-covers effort. Report the BSS (now cross-model stable for the boat) with the PE as a caveated cross-check; the boat gap is mostly the PE undercounting (P0), the shore gap is an effort-imputation tied to `sigma_IE`.
- **Validate by run, never by reasoning alone (the "pin" lesson).** A change that looks inference-neutral on paper (the v7.3 pin; the B1.7 collapse) can perturb a delicate mass-matrix geometry. Isolate each change, compare against a confirmed baseline, and decide on pre-set criteria. Every "done pending run" item is not done until the run confirms it.

---

## 6. Field-protocol items for 2025-26 (not code)

These require sampling-design changes, not code, and should feed the field plan: **jetty and beach direct effort counts** (critique 3 / D3; shore effort currently proxied by dock gear counts); **minimum summer coverage targets by mode** (motivated by T1.1, the estimate's weakest point); **interview-time vs effort-count-time representativeness checks** (critique 10); and confirming effort-count and incomplete-trip protocols. Spatial heterogeneity within the port (Float 20 vs 17-21, critique 11) is a design choice, not currently modeled.

---

## 7. Historical facts worth keeping (from the 2026-03-31 critique)

The critique framed the framework as a sound adaptation of the WDFW freshwater creel BSS to a different fishery, with the biggest risk upstream of the statistics in the effort-to-crabber conversion and day-length expansion. Its 11 issues map to current items as: 1 -> effort unit (largely resolved via deployments) + T1.1; 2 -> T1.3; 3 -> D3 (field); 4 -> RESOLVED (B1_C); 5 -> T3.1/B6; 6 -> T1.4; 7 -> resolved in gear track; 8 -> gear-track low priority; 9 -> ZINB decision; 10 -> field + representativeness; 11 -> weather (T2.4) + shore soak-time (mostly mooted by deployments).

---

## 8. Disposition of the source notes

This file consolidates and supersedes the following. Recommend deleting them from `development_notes/` once this file is confirmed, keeping only the two that retain unique reference value.

| File | Disposition |
|---|---|
| `PLANNED_IMPROVEMENTS.md` | Superseded (its A/B/C/D backlog is reconciled in Sections 3-4). Delete. |
| `CODE_IMPROVEMENTS_REVIEW_v7.0.md` | Superseded (its T1-T4 tiers are carried into Section 4). Delete, or keep for its detailed evidence appendix if wanted. |
| `ADDITIONAL_OUTPUTS_PROPOSAL.md` | Superseded (O1-O11 implemented; see Section 3). Delete. |
| `pipeline_state_review_20260709.md` | Superseded (F1-F5 done, P0-P3 carried forward). Delete. |
| `BOAT_RESOLUTION_EXPERIMENT.md` | Concluded (monthly AR chosen; Section 3). Delete, or keep as the experiment's provenance. |
| `CHANGES-2026-07-11.md` | Point-in-time refactor changelog; retained in the pooled development-history doc. Delete from here. |
| `20260710-OUTSTANDING_ISSUES.md` | Superseded by this file (POOL-*/GR-*/WX-*/ORCH-* reconciled in Sections 3-4). Delete. |
| `20260331-model-critique.docx` | **Keep.** The original external critique is primary-source provenance; Section 7 summarizes it but the docx should remain. |

The full version-by-version change log lives in the two development-history documents (`BSS-GH-pooled-CPUE-model-development-history.md`, `BSS-GH-gear-type-CPUE-model-development-history.md`); this file is the forward-looking status and backlog, not a changelog.
