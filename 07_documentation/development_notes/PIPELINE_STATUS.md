# Coastal Rec Crab BSS: Pipeline Status and Backlog

**Last updated:** 2026-07-13
**Maintainer note:** this is the single living status document for the pipeline. It replaces the six scattered development notes listed in Section 8, reconciling their issue IDs so nothing is lost. Update this file as work lands; do not re-fork it into per-session notes.

**Repo:** `Coastal-Rec-Crab-BSS`, `main`. **Method of record:** Method v1.0 (frozen against pooled code v7.4); the code has advanced to v7.9 with effort-unit, filtering, and config corrections that move published totals, so the 2024-25 reference numbers in the method documents are pre-refresh until regenerated.

---

## 1. Current state at a glance

| Component | Version | State |
|---|---|---|
| Pooled CPUE model (`crab_bss_pooled.stan`, `BSS-GH-pooled-CPUE-model.Rmd`) | pipeline code **v7.9** (2026-07-12) | Production. Both fitted components on gear-deployments; boat reconciled on monthly AR. |
| Gear-resolved model (`crab_bss_gear_resolved.stan`, `BSS-GH-gear-type-CPUE-model.Rmd`) | framework **v5.6** (2026-07-12) | Production. Shore + boat on gear-deployments; `G = 1` (per-gear CPUE inert, gear split from PE apportionment). |
| Weather/tide covariate module (`crab_bss_pooled_weather_adjusted.stan`, `06_diagnostics/...Rmd`) | ~v6.9 parity, **stale** | Not production. Forked engine, pre-deployment-scale Stan, not re-run since the day-length extraction. Do not cite its boat number. |
| Config surface (`run_config.R`) | base-params architecture (P5) | `run_config` is the base parameter set; each driver layers model-specific tuning via `modifyList(run_config, params_model)`. |

**Authoritative run:** `05_output/20260711/pooled-CPUE` (pooled, boat forced to monthly AR via `ar_force`, gate at 0.05). The sibling dirs `-morning` (v7.7) and `-afternoon` (weekly boat, `ar_force = NULL`, gate 0.15) are stale and should be deleted to avoid confusion.

**Headline numbers, that run (Dungeness kept, 2024-25):**

| Component | PE | BSS (median) | BSS vs PE | Method used |
|---|---:|---:|---:|---|
| Shore pot-closure (Sep 16 - Nov 30) | 5,884 | 6,586 | +11.9% | BSS |
| Shore all-gear (Dec 1 - Sep 15) | 15,732 | 20,654 | +31.3% | BSS |
| Private boat all-gear | 22,070 | 43,475 | +97% | BSS |
| Private boat pot-closure | 971 | n/a | n/a | PE (insufficient data) |
| Commercial/charter | 12,007 | n/a | n/a | Census |
| **Port total** | **56,665** | **83,914** | **+48%** | gated combination |

**Convergence:** all three fitted components pass (shore pot-closure daily 2.5% divergences; shore all-gear daily 2.23%; boat all-gear monthly 0.63%); R-hat near 1.00, n_eff well above 400. Boat pot-closure is the only PE fallback and it is data-limited (~17 interviews), not a geometry failure.

**The three facts that most shape how to read the numbers:**
1. **The boat is now cross-model consistent.** On monthly AR the pooled boat (catch 43,475, effort 14,716, CPUE 2.95) matches the independent gear-resolved boat (43,314 / 14,805 / 2.93). The weekly-vs-monthly AR resolution, not the effort unit, was the entire source of the earlier ~26% pooled-vs-gear-resolved boat gap.
2. **The BSS sits systematically above the PE**, most sharply for the boat (+97%). Most of the boat gap is a PE artifact (its stratified mean-of-daily-ratios undercounts in thin strata; see P0), not BSS over-reach; the shore gap is an effort-imputation difference (CPUE agrees). See Section 5.
3. **The harvest is summer-dominated and lean on extrapolation.** July-September and the sparse shoulder months carry wide intervals; even on monthly AR the June boat BSS is ~4x the design-based PE. This is the largest standing threat to the headline number's credibility (T1.1).

---

## 2. Repository map

- **`01_BSS_models/`**, the two production driver `.Rmd` (pooled, gear-resolved) and their rendered `.html`. Pooled is v7.9, gear is v5.6. The `-old.Rmd` snapshots were removed 2026-07-12.
- **`02_stan_models/`**, `crab_bss_pooled.stan` (v7.6+; `R_G_boat`, `effort_scale_gear`/`E_scale`, `collapse_mu_hier` lever), `crab_bss_gear_resolved.stan` (has the non-centered `omega_0`, see D4), `crab_bss_pooled_weather_adjusted.stan` (stale).
- **`03_R_functions/`**, 25 shared modules. Post-refactor: extracted drivers (`fetch_crab_data*`, `run_pe_pooled`/`run_pe_gear`, `prep_bss_crab_pooled`/`_gear`, `prep_days_crab`, `prep_population_summary`, `estimate_comm_charter`, `classify_day_type`, `bss_timers`) plus the diagnostic/engine modules (`bss_ar_resolution`, `bss_convergence_gate`, `bss_cpue_diagnostics`, `bss_effort_spec`, `bss_trailer_expansion`, `bss_day_length`, `diagnose_effort_overdispersion`, `divergence_diagnostic`, `model_diagnostics`, `save_run_diagnostics`).
- **`04_input_files/`**, `effort_combined.csv`, `interview_combined.csv`, `wes_commercial_tally.csv`, `ingress_egress.xlsx`.
- **`05_output/`**, dated run folders. **1,017 git-tracked files (~82 MB)**; a known hygiene problem (T4.2/ORCH-26). The current authoritative run is `20260711/pooled-CPUE`.
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
- **Pot-closure sub-season made explicit and renamed (2026-07-13).** The pot-closure window (when pots are illegal and only non-pot gear is legal) is now set explicitly in `run_config.R` via `pot_closure_start` / `pot_closure_end`, instead of being assumed to run from the season start to `pot_open_date - 1`, so a future season whose start does not coincide with the closure is supported. A shared builder `03_R_functions/build_subseasons.R` derives the sub-seasons for both drivers and handles the general case (optional pre/post all-gear periods for a mid-season closure, keyed off `gear_regime` rather than the sub-season name). The sub-season was renamed "ring-net only" to "pot closure": a **display-only** rename, since non-pot gear other than ring nets is legal, so the old name was a misnomer; the internal key stays `ring_net_only` for output-filename continuity, so keys and filenames are unchanged and the 2024-25 config produces bit-identical sub-seasons. Season plots in both drivers now draw vertical lines at the closure start and the pots-open date. **Known limitation:** the builder is general, but the report's downstream aggregation (PE port summary `component_keys` / `component_names`, catch-by-mode, and the final-table PE-fallback rows) is still hardcoded to the two historical sub-seasons per population; a fail-fast guard in each driver stops the run if a mid-season closure is configured (sub-season name set != {ring_net_only, all_gear}), rather than silently under-reporting the port total. Generalizing those aggregators to sum over all sub-seasons per population is the follow-up (relates to T3.1).

**PE (partial).**
- Shore PE unit alignment done (v7.8, above). Ratio-of-sums adoption (P0) is NOT done and is the top open PE item (Section 4).

**Repo cleanup (partial).**
- `-old` driver snapshots and the refactor tarball removed (2026-07-12).

---

## 4. Outstanding work (prioritized by publication impact)

### Tier 1: can move the headline number, or a reviewer will demand it

Status note: the four Tier 1 items were worked on 2026-07-12. Three are implemented in code (P0, T1.1a, T1.4) and need a confirming run; the prior-sensitivity harness (T1.3) is in place but the sweep runs are pending; external validation (T1.1b) still needs benchmark data.

- **P0 [DONE 2026-07-12, confirm on next run]: pooled PE on ratio-of-sums.** `run_pe_pooled` now computes stratum CPUE as within-stratum `sum(catch)/sum(hrs)`, replacing the unstable `weighted.mean(daily_ratios, w = n_int)`, and carries an implied-CPUE-vs-ratio-of-sums guard that fails fast (in the PE section, before the multi-hour BSS fits) if catch and effort drift onto different scales. This mirrors `run_pe_gear`, so the two tracks' PE cannot diverge. Expected effect on the next run: the boat PE moves from 22,070 toward ~40k (near the BSS and the gear-resolved boat), and the shore PE implied CPUE drops from ~2.9x its ratio-of-sums toward 1x, making the PE a fair cross-check. The estimator change and its per-component internal-consistency check now render in the knitted report (Section 4.4).
- **T1.1(a) [DONE 2026-07-12]: summer-extrapolation transparency.** A new report section (`extrapolation_transparency.csv` plus an on-page table, RMD Section 11.4) decomposes each fitted component's BSS catch into interview-days vs extrapolated (no-interview) days, and reports the share of the summed BSS catch that rests on extrapolated time. Turns the summer-dominance caveat into a number.
- **T1.1(b) [OPEN, needs external data]: external validation (B4).** Still requires a benchmark the repo does not contain: dockside census totals, a known-total year, or the legacy estimator. The strongest internal cross-checks now in place are the two independent BSS pipelines agreeing on the boat (43.5k) and the per-component PE-vs-BSS comparison (fair once P0's run lands); neither is a true external check. Pursue a benchmark before publication.
- **T1.3 [HARNESS DONE 2026-07-12, sweep runs pending]: R_G prior sensitivity.** The R_G prior is now overridable via `params$R_G_prior_mu` / `R_G_prior_sigma` (data-driven when unset), with a commented sweep block in `run_config.R`. Remaining work is the runs: set `R_G_prior_mu` to 1.0, ~1.28 (empirical), and 1.5 in turn and tabulate the port totals, then extend to the `L_effective` prior (`tau_shore`/`tau_boat` are already live toggles). `prior_vs_posterior_*.csv` (O9) shows the prior pull. **Effort: the runs (multi-hour each).**
- **T1.4 [DONE 2026-07-12, confirm on next run]: commercial/charter vessel split.** `estimate_comm_charter` now applies separate per-vessel Dungeness (and red-rock) means for commercial vs charter vessels (from `boat_type_clean`, Guide folded into Charter) to the matching tally columns, with a fallback to the pooled mean for a class with no interviews. In the 2024-25 census window commercial is ~35 and charter ~51 crab/vessel, so the split corrects a real bias. Applies to both tracks (shared estimator). The per-vessel-type split renders in the knitted report (Section 4.2). Pair with a census uncertainty interval (T3.5).

### Tier 2: accuracy, cross-track parity, or credibility

- **Make monthly the boat AR production default.** Monthly is currently reached via `ar_force = list(private_boat = "monthly")`, an experiment override. Set `ar_max_resolution$pooled$private_boat = "monthly"` in `run_config.R` so the production run selects it by policy (it is better-identified for the thin trailer series and now agrees with the independent pipeline), and set `ar_force = NULL`. **Effort: trivial.**
- **T2.2 / B5: propagate `pi_gear` uncertainty in the pooled gear breakdown, or scope it as approximate.** `catch_by_gear_type.csv` applies point interview proportions post-hoc, so its gear-type intervals are too narrow. Either mirror the gear-resolved Dirichlet propagation or add a one-paragraph scoping note. **Effort: low (note) to medium (propagation).**
- **T2.4 / C1-C2 + WX module: weather covariate evaluation, then integration; rehabilitate the module.** The GAM screen ran for all four cells but the committed PSIS-LOO is only shore ring-net (marginal, ~2.4 SE, below the 4-SE bar). The module (WX-16..20) forks a stale v6.9-parity engine on a pre-deployment-scale Stan and has not been re-run; it must be re-based on the shared modules and current Stan before its LOO or any covariate promotion (C2) is trustworthy. **Effort: medium.**
- **T2.5: durable single-cell mu-hierarchy funnel fix.** With `G*S = 1`, `mu = mu_mu + sigma_mu * eps_mu` is over-parameterized and is the residual funnel the divergences sit on (`sigma_mu_*` is the top divergence driver each run). The `collapse_mu_hier` lever (POOL-4) exists as a safe off-by-default toggle; the naive full collapse (B1.7) was reverted because it forced the level to reconcile against the daily AR ridge. The durable fix is a careful partial reparameterization (single decoupled level absorber), a design task to prototype on shore all-gear, not a mechanical edit. Also make the gear-resolved collapse a config flag (GR-10). **Effort: medium.**
- **Shore `sigma_IE` effort-model tension (GR-9 and pooled analogue).** In the pooled monthly run shore all-gear `sigma_IE` = 1.07 (worse than the gear-resolved 0.60): the I/E day-length/turnover stream and the effort counts disagree by ~e^1 at the tails. Investigate whether the I/E observation days are unrepresentative (peak-count days) before trusting the shore effort posterior to three figures. **Effort: medium.**
- **GR-8: `sigma_IE` funnel in `shore_ring_net`.** `exponential(5)` has its mode at zero so `sigma_IE` can shrink and stiffen the likelihood on the 2 in-window I/E days; add an `ie_min_obs_shore >= 3` guard or tighten the prior to `lognormal(log(0.3), 0.5)`. **Effort: low.**

### Tier 3: refinements

- **T3.1 / B6 / critique 5: sub-season boundary edge effect.** The hard Dec 1 split fits two independent BSS with no cross-boundary anchor, ballooning the interval at the seam. Either document as a known artifact or move to a single regime-switching model with a `pot_available` indicator. **Effort: medium (doc) to high (model).**
- **T3.3 / critique 4 extension: continuous density-dependence and a holiday CPUE term.** CPUE as a function of same-day effort rather than only a weekend dummy; effort has a holiday term (`B2`), CPUE does not. Lower priority since the dominant day-type signal is already modeled. **Effort: medium.**
- **GR-7: genuine per-gear CPUE (Option A).** `G = 1`, so per-gear catch comes from PE apportionment and `catch_by_gear_type.csv` mixes a boat "All" BSS column with PE-apportioned shore gear types (ORCH-24). A real build (effort shares in the offset, a multi-gear rule), not a fix; not to be started until the boat total is stable. Related: critique 7 (priority hierarchy) is resolved in the gear track via weighted pseudo-likelihood; critique 8 (shared `phi_C_gear`) remains a low-priority gear-track item. **Effort: high.**
- **GR-12: `tau_boat` identification.** Boat catch is proportional to `tau_boat = 1.2`, whose prior rests on 2 WBL I/E days; the boat I/E stream is built and activates once >=2 days fall in a window. The egress-classification pilot is the critical path. State this dependence in any co-management number. **Effort: medium (pilot).**
- **Zero-inflation / ZINB (critique 9 / B2 follow-on): decide on evidence.** The catch likelihood is NB2 with no zero inflation; a high zero fraction (snare, ring-net) may stretch `r_C`. Read the PPC zero-bin and per-observation residuals (O5); prototype ZINB (a variant exists in the freshwater lineage) on the worst-calibrated component only if the zero bin is systematically off, compare by PSIS-LOO. **Effort: low to read; medium if warranted.**
- **T3.2 / critique 11: shore soak-time basis.** Largely mooted now that shore is on gear-deployments rather than crabber-hours, but a one-line confirmation in the docs is worth it. **Effort: trivial.**

### Tier 4: reproducibility and repository hygiene (correctness-neutral, publication-blocking)

- **renv.lock (A4 / T4.1 / ORCH-25).** Still absent. A publication method with a multi-hour Stan pipeline needs a pinned dependency set; ~5 minutes with a large reproducibility payoff. Add a top-level runner too (`run_estimation.R` already serves this).
- **Output and repo hygiene (A1/A2 / T4.2 / ORCH-26).** `05_output/` is 1,017 tracked files / ~82 MB and grows each run, including pre-v6 runs not reproducible from current code; `.Rproj.user/` (5 files) is tracked despite being gitignored. Untrack `.Rproj.user`; decide an `output/` policy (gitignore + one canonical reference run, or Releases). Delete the stale `20260711/pooled-CPUE-morning` and `-afternoon` runs.
- **A3: gear-resolved `max_divergences` doc/code mismatch.** Code uses 5; older docs said default 0. Reconcile the doc. **Effort: trivial.**
- **Stale gear-resolved output labels.** After the v5.5 shore->deployments move, the PE plot title "Estimated Daily Crabber-Hours" and the port-total `effort_units_note` still say crabber-hours; one-line fixes.
- **Global-variable coupling (T4.3), duplicated PE monthly-share blocks (T4.4).** Largely addressed by the refactor (functions now take `params`); verify no residual globals and that the two PE monthly-share blocks are one helper.

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
