# Grays Harbor Recreational Dungeness Crab Harvest Estimation

## Method Version 1.0: Pooled CPUE Model

**Author:** Matthew George, Ph.D.
**Contact:** matthew.george@dfw.wa.gov
**Agency:** Washington Department of Fish and Wildlife (WDFW)
**Status:** Published, operational. The method of record for estimating recreational Dungeness crab harvest at Westport / Grays Harbor.
**Method version:** 1.0. "Method v1.0" is the frozen method label: the model structure, the estimators, and the design decisions in this document. It is distinct from the pipeline "code" revision (v7.x), which tracks implementation changes. Method v1.0 was first frozen against code v7.4.
**Current pipeline code:** v7.9 + Tier-2 batch (2026-07-13). Since v7.4 the code has advanced in ways that move published totals: v7.5 added the incomplete-trip filter (raises the shore estimate); v7.6 moved the private boat onto gear-deployments; v7.7 moved the shore BSS onto gear-deployments (moves the shore total); v7.8 completed the shore PE onto the same unit and refactored the code; v7.9 restructured `run_config.R` and tightened the pooled divergence backstop (no totals move). The 2026-07-11 v7.8 run refreshed the totals. One correction to an earlier expectation: the v7.6 boat move proved catch-neutral (on weekly AR the boat holds near 54,481, NOT the ~-25% originally predicted); it corrects the boat's effort unit and interpretation but does not move the boat harvest. What does move the boat is the AR resolution: on monthly boat AR (`05_output/20260713/pooled-CPUE-run1`) the boat reconciles to ~43,180, matching the gear-resolved pipeline, so monthly is the recommended production boat AR cap (see Sections 16 and 19). The method is unchanged by these; they are effort-unit and filtering corrections, not a re-derivation. Numbers quoted throughout that predate the 2026-07-11 run should be read as pre-refresh.
**Reference calibration season:** 2024-25.

------------------------------------------------------------------------

## How to read this document

This is the single authoritative reference for the pooled-CPUE harvest estimation pipeline (`01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd`), the Stan model it fits (`02_stan_models/crab_bss_pooled.stan`), the diagnostics it runs, and the inputs and outputs it uses. It is written for two audiences and is split into three parts:

- **Part I (Sections 1-6): For everyone.** Plain-language description of what the method does, the fishery, the data, how the estimate is built, and where it is valid. No statistics background required.
- **Part II (Sections 7-13): Running it next season.** The operational guide: prerequisites, step-by-step run, how to judge whether a season's estimate is trustworthy, the output catalog, the diagnostics, reproducibility, and the conditions under which the method stops applying.
- **Part III (Sections 14-20): Technical reference.** The full model specification, design rationale, limitations, glossary, and references.

The development history (how the model reached v1.0, the full change log from v3 through the current code version, and the convergence-debugging narrative) has been moved out of this document to keep it focused on the published method. It lives in `BSS-GH-pooled-CPUE-model-development-history.md`. Section 19 gives a one-screen summary and points there.

A note on naming: this is the "pooled" model because it uses a single catch-rate (CPUE) process shared across gear types. A separate model, documented in `BSS-GH-gear-type-CPUE-model-documentation.md`, instead estimates a CPUE process per gear type. The pooled model is the published v1 because it is the more robust of the two and answers the primary management question (total harvest with defensible uncertainty); the gear-resolved model is the alternative when modeled gear-type catch with full uncertainty is required.

---

## PART I: FOR EVERYONE

------------------------------------------------------------------------

### 1. What this method produces

This framework estimates the total recreational Dungeness crab (*Metacarcinus magister*) harvest at Westport and the greater Grays Harbor area for a season. It combines four kinds of field observations (gear counts at the docks, trailer counts at the boat launch, dockside crabber interviews, and ingress/egress surveys) with a statistical model that fills in the days when no sampling occurred.

Each run produces:

- A total Dungeness crab harvest estimate for the port, with a 95% credible interval (a range that has a 95% probability of containing the true harvest, given the data and model).
- Monthly harvest trends showing when crabbing pressure peaks and how it changes through the season.
- Breakdowns by crabbing mode (shore, private boat, commercial/charter) and an approximate breakdown by gear type (pot, ring net, trap, snare).
- A weekend catch-rate effect (whether weekend catch rates differ from weekday rates).
- Daily estimates of "effective day length" at the docks when ingress/egress data are available.

**How confident are we?** The framework runs two independent estimation methods, a simple average-based approach (the Point Estimator, PE) and a Bayesian time-series model (the Bayesian State-Space model, BSS), then compares them. When the two agree and the BSS passes its convergence checks, confidence is high. The output includes formal diagnostics and a side-by-side comparison so a reviewer can judge reliability. Section 9 explains how to read those checks.

------------------------------------------------------------------------

### 2. The fishery and study area

The recreational Dungeness crab fishery at Westport is one of the highest-volume recreational crabbing operations on the Washington coast. Crabbers use four main gear types: crab pots (highest catch rate), ring nets, foldable/star traps, and snares. WDFW rules prohibit pots from late September through November, which creates a structural break in both effort and catch rates and is the reason the season is split into two sub-seasons (Section 5).

Commercial Dungeness crab vessels also crab recreationally before the commercial season opens, under the same daily limits as private boats. Their harvest is tracked separately through a vessel tally at the marina.

Westport sits on the south side of the Grays Harbor estuary. Recreational crabbing occurs from public docks (Floats 17-21), a jetty, beaches, a public boat launch, and the commercial marina. The "shore" component pools dock, jetty, and beach crabbing.

------------------------------------------------------------------------

### 3. The four data streams

| Stream | What is collected | What it tells the model |
|---|---|---|
| **Effort counts** | Instantaneous point-in-time counts of crab gear at the docks and boat trailers at the launch, by field surveyors | The primary indicator of how much crabbing activity is happening |
| **Crabber interviews** | Dockside trip-level records: group size, gear deployed and type, hours fished, crab kept, trip status | Catch rate (CPUE) and the mix of gear in use |
| **Commercial/charter tally** | Daily count of commercial and charter vessels at the marina during the recreational pre-season | The commercial/charter component of harvest, via expansion |
| **Ingress/egress (I/E) surveys** | All-day surveys recording crabber arrivals and departures every 15 minutes | A direct measurement of crabber-hours that calibrates the gear-count pathway |

The four input files that carry these streams are listed in Section 7; their exact schema and known quirks are documented in `04_input_files/README.md`.

------------------------------------------------------------------------

### 4. How the estimate is built

The core problem is that field crews cannot sample every day. Both methods solve the same problem (estimate harvest on unsampled days), but differently.

**The Point Estimator (PE): a simple average.** For each stat-week by day-type group, the PE averages the daily harvest on sampled days and multiplies by the number of days in that group (Pollock et al. 1994; Hahn et al. 2000). It is transparent and assumption-light, but it cannot fill a group that had zero samples and it produces no uncertainty bounds.

**The Bayesian State-Space model (BSS): a time-series curve.** The BSS fits a smooth curve through the daily effort and catch-rate data using a statistical time-series process, then uses that curve to estimate every day in the season, including unsampled days (Conn 2002; Staton et al. 2017). It accounts for the fact that adjacent days are correlated, fills gaps with honest uncertainty that grows the further a day is from the nearest observation, and produces credible intervals. It is more complex, takes roughly half an hour to a few hours per fit, and must be checked for convergence.

**Combining them.** For each population component, the framework checks the BSS fit against formal convergence criteria (Section 9). If the fit passes, its estimate is used; if not, the PE estimate is used as a fallback. The two are reported side by side so a reviewer can see where they agree and where they differ.

The headline harvest number uses the BSS posterior expected catch (the model's best estimate of the average catch) rather than a single simulated draw, following the standard distinction between estimation and prediction in hierarchical models (Gelman et al. 2013, Ch. 7).

------------------------------------------------------------------------

### 5. The three population components and the two sub-seasons

The harvest is built from three components, estimated separately and summed:

1. **Shore crabbers** (dock + jetty + beach). Effort is indicated by gear counts at the docks.
2. **Private boat crabbers.** Effort is indicated by trailer counts at the boat launch.
3. **Commercial/charter vessels** crabbing recreationally pre-season. Estimated by expanding the marina vessel tally.

The season is split into two **sub-seasons**, defined by the pot closure, and each is estimated independently:

- **Pot closure** (Sep 16 to Nov 30): pots prohibited; non-pot gear (ring nets, snares, traps) allowed. Formerly labeled "ring-net only", a misnomer since gear other than ring nets is also legal during the closure. The internal key stays `ring_net_only` for output-filename continuity; the reports display "Pot closure".
- **All-gear** (Dec 1 to Sep 15): pots allowed.

The closure window is set explicitly in `run_config.R` via `pot_closure_start` and `pot_closure_end` (added 2026-07-13), so a future season whose start does not coincide with the closure start is supported. The shared builder `03_R_functions/build_subseasons.R` derives the sub-seasons from that window and adds pre/post all-gear periods automatically if the closure falls mid-season. Season plots mark the closure start and the pots-open date with vertical lines.

------------------------------------------------------------------------

### 6. Where this method is valid

Method v1.0 is calibrated to the Westport / Grays Harbor fishery as sampled in the 2024-25 season. It is designed to be re-run in future seasons **provided the fishing location, the input data streams, and the sampling design remain the same.** Section 13 sets out, in detail, which assumptions are baked in and the specific conditions under which the method must be re-derived rather than re-run. In short: a different port, a change in how effort counts are taken (for example, reverting from randomized counts to a single peak-time count), or a structural change in who participates would each require revisiting the method, not just feeding it new data.

---

## PART II: RUNNING IT NEXT SEASON

------------------------------------------------------------------------

### 7. Prerequisites and repository layout

**Software.** R 4.2 or later, with rstan 2.32 or later and a working C++ toolchain (rstan compiles the model), plus the packages tidyverse, lubridate, suncalc, gt, patchwork, here, and readxl.

**Repository layout.** The pipeline relies on the numbered stage folders (`01_BSS_models` through `05_output`) and on `here::here()`, which anchors all file paths to the repository root. You do not edit paths to run it; you edit one config file and place the input files in the right folder:

| Location | What it holds | Your job |
|---|---|---|
| `run_config.R` (repo root) | The single control surface: as of the 2026-07-11 consolidation this is the single source of truth for every user-selectable toggle (season window, structural dates, catch groups, day-typing, effort unit, filters, I/E settings, holidays, and model-behavior toggles) | This is the one file you edit season to season |
| `run_estimation.R` (repo root) | The orchestrator that injects `run_config` and renders the chosen driver | Launch it (`source()` in RStudio, or `Rscript`); nothing to edit |
| `01_BSS_models/` | The driver `BSS-GH-pooled-CPUE-model.Rmd` | Run via `run_estimation.R`, or knit standalone (its setup chunk sources `run_config.R` automatically) |
| `02_stan_models/` | `crab_bss_pooled.stan` | Leave in place; the driver compiles it |
| `03_R_functions/` | The helper functions, sourced wholesale by the driver's setup chunk. The 2026-07-11 refactor extracted the driver's inline helpers into standalone files here, both shared (`bss_effort_spec.R`, `bss_timers.R`, `classify_day_type.R`, `prep_days_crab.R`, `prep_population_summary.R`, `estimate_comm_charter.R`) and pooled-specific (`fetch_crab_data.R`, `run_pe_pooled.R`, `prep_bss_crab_pooled.R`) | Leave in place; auto-sourced |
| `04_input_files/` | The four input files | Replace these with the new season's data, same names and schema |
| `05_output/` | Dated run folders | The run writes here; nothing to place |

**The four input files** (see `04_input_files/README.md` for schema and quirks):

- `effort_combined.csv` (effort counts)
- `interview_combined.csv` (crabber interviews)
- `wes_commercial_tally.csv` (commercial/charter tally)
- `ingress_egress.xlsx` (I/E surveys; named through the `ie_data_file` / `ie_sheet` parameters)

------------------------------------------------------------------------

### 8. Step-by-step: running a season

1. **Place the new season's data** in `04_input_files/`, keeping the four filenames and their column schemas unchanged. Honor the schema quirks in the input-folder README (the interview gear column maps from column N; re-export the effort CSV with full quoting; dates are M/D/YYYY; the "Commerical" boat-type spelling is matched by regex).
2. **Edit `run_config.R` (repo root). This is the ONE file you edit.** As of the 2026-07-11 consolidation it is the single source of truth for every user-selectable toggle; the driver's `params` chunk holds only this model's internal tuning (the Stan file, per-fit sampler settings, gate thresholds, AR-selector thresholds), which rarely changes and legitimately differs from the gear-resolved model. The keys you normally set:
   - `est_date_start`, `est_date_end`: the season window. The driver fits each sub-season inside this window.
   - `pot_open_date`, `pot_closure_start`, `pot_closure_end`, `commercial_opener`, `census_start_date`, `census_end_date`: the regulatory / structural dates. `pot_closure_start` / `pot_closure_end` bound the pot-closure sub-season explicitly (keep `pot_open_date` = `pot_closure_end` + 1).
   - `shore_effort_unit` (default `gear-deployments`), `filter_incomplete_trips` (default on), and the `tau_shore` / `tau_boat` turnover priors: the effort-unit and CPUE-denominator controls (Sections 11, 14, 15).
   - `ie_data_file`, `ie_sheet`: the I/E workbook and sheet (defaults `ingress_egress.xlsx`, `data`).
   - `crabbing_holiday_dates`: update this one list each season.
   - `bss_seed`: the RNG seed (default 20260619). Leave fixed for reproducibility; change only if a pathological seed is ever suspected.
3. **Run the pipeline.** Use `source("run_estimation.R")` in RStudio (Source, not Knit) or `Rscript run_estimation.R` from a terminal; the orchestrator injects `run_config` into the driver and renders it. You can also knit `01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd` directly, because its setup chunk sources `run_config.R` automatically when `run_config` is not already present, so a standalone knit uses exactly the same toggles. Each population by sub-season is fit independently. Expect a total runtime of roughly 3 to 6 hours on a 4-core machine, depending on AR resolution and sub-season length.
4. **Check convergence** for each fit using `convergence_report.csv` and the rules in Section 9. A fit that fails falls back to PE automatically.
5. **Read the outputs** from `05_output/YYYYMMDD/pooled-CPUE/` (Section 10), starting with `port_total_Dungeness_Kept.csv` and `pe_vs_bss_comparison.csv`.

------------------------------------------------------------------------

### 9. Judging whether a season's estimate is trustworthy

This is the most important section for an operator. For each BSS fit, the framework monitors four diagnostics, reported per fit in `convergence_report.csv`: rank-normalized split-R-hat and bulk effective sample size (n_eff) for the seasonal totals `C_expected_sum` and `E_sum` (Vehtari et al. 2021), the number of divergent transitions, and the percentage of iterations that saturate the sampler's tree depth.

A fit **passes**, and its BSS estimate is used for that component, when all of the following hold; otherwise the PE estimate is used:

- **R-hat < 1.01** for both totals. R-hat near 1.00 means the independent sampler chains agree.
- **n_eff > 400** for both totals. This is the effective number of independent posterior samples.
- **Divergent fraction below 0.05** (the hard backstop). Above this rate, the sampler's geometry is untrustworthy and the fit is rejected regardless of anything else.
- **Divergences do not move the answer.** The shift the divergent draws induce in each total, measured in units of that total's posterior standard deviation (`|median(all) - median(bulk)| / sd(all)`), is below 0.10 SD.

Why divergences are in the gate: a sampler that cannot accurately integrate its trajectory is not faithfully exploring the target distribution, and can bias the posterior even when R-hat and n_eff look fine (Betancourt 2017).

Why the impact criterion is measured in standard deviations and not as a percentage of the estimate: a percentage-of-level threshold penalizes a component for having a wide posterior (being weakly identified) rather than for being biased. The SD-normalized criterion asks the question the gate exists to answer (do the divergences move the answer relative to how well the answer is pinned down) and is invariant to how wide the posterior is. This matters directly for the private boat, whose posterior is genuinely wide; see Section 16.

**Reading the comparison.** `pe_vs_bss_comparison.csv` shows PE and BSS effort and catch by component with the selected method. Large PE-vs-BSS gaps are not automatically errors; they can reflect a real disagreement between the design-based expansion and the model's reconciliation against interview data (the private boat is the standing example, Section 16). The convergence report and the comparison now always agree on which method was used, because the gate decision is computed once per fit and consumed by every downstream summary (this was a v7.0 fix; see Section 19).

One unit-consistency correction, from the 2026-07-11 refactor, bears on this comparison for shore. Before the refactor the pooled shore PE branch still computed effort as crabber-hours (`est_crabbers * day_length`) and CPUE per crabber-hour, even though v7.7 had moved the shore BSS to gear-deployments, so the shore PE and shore BSS were on different units and the shore row of `pe_vs_bss_comparison.csv` and the monthly PE effort-share for shore were unit-inconsistent. The shore PE branch now flows through the shared effort-unit spec (`03_R_functions/bss_effort_spec.R`), the same one the shore BSS uses, so shore PE matches the shore BSS unit. This moves the pooled shore PE number and, like the v7.5 to v7.7 changes, is pre-refresh: it requires a validation run before the shore PE figures are cited (see the header caveat).

------------------------------------------------------------------------

### 10. Output catalog

Each run writes to `05_output/YYYYMMDD/pooled-CPUE/`. Files tagged with a population follow the pattern `<metric>_<population>_Dungeness_Kept.{csv,png}`, where population is one of `shore_ring_net_only`, `shore_all_gear`, or `private_boat_all_gear`. The commercial/charter component has no separate BSS file; it enters the port total by census expansion.

**Headline estimates**

| File | Contents |
|---|---|
| `port_total_Dungeness_Kept.csv` | Combined PE + BSS port total (expected and predictive) |
| `pe_port_summary.csv` | PE estimates by component and port total |
| `pe_vs_bss_comparison.csv` | PE vs BSS effort and catch by component, with the selected method |
| `monthly_estimates.csv` | Monthly catch and effort with credible intervals |
| `monthly_estimates_by_mode.csv` | Monthly catch by crabbing mode with 95% intervals |
| `catch_by_mode.csv` | Catch by crabbing mode (shore, boat, commercial) |
| `catch_by_gear_type.csv` | Approximate catch by gear type (proportional allocation) |
| `season_summary.csv` | Season totals roll-up |
| `sensitivity_incomplete_trips.csv` | PE catch, effort and gear ratio under four incomplete-trip treatments (exclude / gear_only / impute_mean_cpue / keep), each expressed against the production `exclude` arm, plus a length-bias test on the gear counts (2026-08-25; was filter-off vs filter-on) |
| `ar_escalation_log.csv` | One row per BSS fit attempt: AR resolution tried, period count, divergences, divergence fraction, impact in posterior SD, and the gate verdict that decided whether to escalate (2026-08-25). With `ar_escalate = FALSE` this is one row per fit |

**Convergence and structure (per fit)**

| File | Contents |
|---|---|
| `convergence_report.csv` | Per-fit R-hat, n_eff, divergent count and fraction, tree-depth, AR resolution, the SD-normalized divergence impact (`impact_C_sd`, `impact_E_sd`; the gating criterion), and the retained level-distortion (`distortion_C`, `distortion_E`; reported only, no longer gating) |
| `structural_params_<label>.csv` | Posterior summary of scale/structural parameters (sigma_eps, phi, r, sigma_mu, sigma_IE, R_G, R_G_boat) with CI, n_eff, R-hat |
| `divergence_localization_<label>.csv` | Where divergent draws sit relative to the bulk, per parameter |
| `sampler_diagnostics_<label>.csv` | HMC sampler diagnostics including E-BFMI |
| `prior_vs_posterior_<label>.csv` | Prior vs posterior comparison per fit |

**Posterior predictive and cross-validation (per fit)**

| File | Contents |
|---|---|
| `ppc_calibration_<label>.csv`, `ppc_pit_<label>.png` | Posterior predictive coverage and PIT for effort counts and interview catches |
| `ppc_byobs_<label>.csv` | Per-observation PPC residuals (exact randomized PIT) |
| `effort_overdispersion_decomp_<label>.csv`, `effort_overdispersion_byobs_<label>.csv` | Effort-variance decomposition (Section 11) |
| `loo_summary_<label>.csv`, `loo_pointwise_*_<label>.csv` | PSIS-LOO summaries and pointwise contributions by likelihood component |
| `cpue_estimators_<label>.csv`, `cpue_saturation_<label>.csv`, `cpue_linearity_<label>.csv` | CPUE effort-unit checks (v7.5): estimator triad (ratio-of-sums vs model-implied vs mean-of-ratios), saturation exponent, and effort linearity; flag when catch does not scale with the chosen effort denominator (Section 11) |

**Effort, day length, and parameters**

| File | Contents |
|---|---|
| `effort_cpue_multipliers.csv` | B1, B2, B1_C posteriors |
| `expansion_ratios.csv` | R_G, R_G_boat posteriors |
| `bss_L_effective_<label>.csv` | Daily effective-day-length posteriors (prior, median, 95% CI) |
| `L_effective_ie_detail.csv`, `ie_analysis.csv` | I/E regression predictions vs observed, and I/E validation |
| `bss_daily_effort_<label>.csv`, `bss_daily_cpue_<label>.csv`, `bss_daily_catch_<label>.csv` | Posterior daily series |
| `bss_summary_<label>.csv`, `bss_full_summary_<label>.csv`, `bss_ar_path_<label>.csv`, `bss_period_coverage_<label>.csv`, `bss_draws_summed_<label>.csv` | Per-fit summaries and the AR path/coverage |

**Plots and metadata**

Plots (`plot_*.png`) cover the daily series, posteriors, monthly catch (total and by mode), the L_effective regression, and the day-length comparison. `run_parameters.txt` and `session_info.txt` record the exact parameters and the R/package/Stan session and seed for the run.

A complete, categorized listing (including how older runs differ and how the weather-tide module's outputs look) is in `05_output/README.md`.

------------------------------------------------------------------------

### 11. Diagnostics: what each one answers

The diagnostics are additive (each is wrapped so a failure cannot break a run) and are written every run. The three that an operator should be able to read:

**Posterior predictive checks (PPC).** `ppc_calibration_<label>.csv` and `ppc_pit_<label>.png` ask whether the model's predictions are calibrated against the actual effort counts and interview catches. A well-calibrated model has PIT values spread uniformly; a central hump means the predictive is too wide (over-dispersed), and 50% coverage above the nominal 0.50 says the same. In the reference run the effort predictive is somewhat over-dispersed (gear/trailer 50% coverage around 0.63 to 0.75), which is what the next diagnostic dissects.

**Effort over-dispersion decomposition.** `effort_overdispersion_decomp_<label>.csv` splits each effort observation's predictive variance into three additive parts via the law of total variance, so the lever behind any over-dispersion is identified before any prior or model change:

```text
Var(Y) = E[mu]            (Poisson floor: irreducible, not a lever)
       + E[mu^2 / r_E]    (NB observation over-dispersion: controlled by the r_E / sigma_r_E prior)
       + Var(mu)          (latent process + parameter uncertainty: controlled by sigma_eps_E)
```

The `lever` column reports the verdict. The decision rule: if the NB-overdispersion share dominates, the lever is the `r_E` / `sigma_r_E` prior (the cheaper, exact change); if the latent share dominates, the lever is the AR innovation scale, which is a more delicate change (the boat tends to show a larger latent share). The analytic decomposition was checked against a brute-force Monte Carlo predictive variance and matches within Monte Carlo noise. Two standing cautions apply: any such correction is a prior/inference change that needs a guarded test run, and tightening the effort dispersion narrows the reported intervals (including the headline summer intervals), which is a change to reported uncertainty and needs explicit sign-off. The target is calibration (50% coverage near 0.50), not zero over-dispersion; some over-dispersion is real.

**PSIS-LOO.** `loo_summary_<label>.csv` reports out-of-sample predictive performance (expected log predictive density, `elpd_loo`) and the Pareto-k influence diagnostic per likelihood component (gear/trailer/catch). This is the basis for principled model comparison; it is what was used to evaluate, and reject, weather covariates (Section 17).

**CPUE effort-unit checks (v7.5).** `cpue_estimators_<label>.csv`, `cpue_saturation_<label>.csv`, and `cpue_linearity_<label>.csv` test the likelihood's core assumption that catch is proportional to the chosen effort denominator `h`. The estimator triad reports the model-implied CPUE (`C_expected_sum / E_sum`) against the ratio-of-sums and the mean-of-ratios; a model sitting near the mean-of-ratios is a warning that the negative-binomial dispersion is pulling `lambda_C` off the rate scale. The saturation exponent fits `catch_per_gear ~ (hours_per_gear)^beta` (boat only) and the linearity check fits `glm(catch ~ log(h))`; the likelihood assumes `beta = 1`, so a value well below 1 means the effort unit is not valid (for pots, catch is nearly flat in soak time). The run also asserts that effort `E` and the CPUE denominator `h` carry the same unit. These are diagnostic only; they are what would surface a boat or shore effort-unit defect before its total is trusted.

**Why the deployment is the effort unit (saturation).** Binned by soak time, crab per gear-HOUR falls about 43-fold across the range of soak durations, while crab per gear per trip rises only about 1.8-fold; a log-log fit gives catch per gear scaling as soak-hours to the power ~0.13. In plain terms, soak time barely matters, so a set pot is a set pot whether it soaked two hours or eight. That makes the deployment the unit on which catch-per-unit-effort is a stable rate: roughly 4 to 7 crab per gear-deployment, steady across soak times. (Terminology, corrected 2026-08-25: a deployment is **a piece of gear a crabber had in the water on that trip**, which is what `number_of_gear` records. It is not a "pot lift". A crabber who sets four pots and checks each three times contributes four deployments, not twelve; the repeat use of a gear slot across the day is carried separately by the turnover `tau`, at the slot level. The earlier "pot lift" wording implied effort scales with checks, which this model does not assume and the data do not support.) A stable rate is exactly what the harvest method needs, because harvest is effort multiplied by that rate, and the multiplication is only unbiased if the rate does not drift with the effort denominator. On the 2026-07-10 shore comparison this is exactly what the linearity diagnostic shows: gear-deployments is the only shore unit whose `beta_h` covers 1 (1.05, 95% CI 0.94 to 1.15), while crabber-hours (0.57) and gear-hours (0.73) both fall well short, so v7.7 moves shore, like the boat before it (v7.6), onto the deployment scale.

**Incomplete-trip filter and its treatment diagnostic (v7.5; extended 2026-08-25).** CPUE is computed from completed trips only (`filter_incomplete_trips`, default on): incomplete trips have soak-time gear that has not finished and read systematically low (about -21% pots, -23% traps, -20% snares; ring nets are effectively exempt at +4%, being checked every few minutes). Missing trip status is kept (a blank `completed_trip` may still be a complete trip). Dropping them costs about 36% of the interview sample.

`sensitivity_incomplete_trips.csv` now reports **four** treatments rather than filter-on vs filter-off, and the production estimator is unchanged by it. The framing that decides what any treatment can do: **effort in this pipeline does not come from interviews.** It comes from the dock gear counts and the trailer/OSP boat counts. Interviews supply only the catch rate and the gear ratios (`R_G`; `R_G_boat` and the PE's gear-per-group), so no treatment can add effort to the expansion, and each arm moves the estimate through one of those two channels and nothing else.

| Arm | What it does | What to expect |
|---|---|---|
| `exclude` | production: drop them from the catch rate **and** the gear ratios | the reference row |
| `gear_only` | drop the truncated catch, keep the fully observed gear count in the gear ratios | the arm with real content; see below |
| `impute_mean_cpue` | keep them with catch replaced by (complete-trip rate x their gear count) | close to a no-op on the catch rate, by arithmetic: a ratio-of-sums over complete-plus-imputed returns the complete-trip rate identically. Any movement is thin-stratum shrinkage |
| `keep` | no filter (pre-v7.5) | the low-biased comparison |

The finding worth acting on is `gear_only`. An interrupted trip's gear count is fully observed — the crabber has N pots out whether or not they are done — and only the catch is truncated, yet `prep_bss_crab_*()` builds the `R_G` / `R_G_boat` interview set from the **already-filtered** frame, so today that gear count is discarded along with the catch. The PE and the BSS also disagree about this: `run_pe_*()` computes the boat gear-per-group from the unfiltered set, so the boat PE already behaves like `gear_only` while the boat BSS behaves like `exclude`. The diagnostic reports the gear ratio under each arm so that inconsistency is visible, and reports a Welch test of incomplete-vs-complete gear counts, because the one thing that would disqualify `gear_only` is length-biased sampling: intercepted trips over-represent long trips by construction, and if long trips also deploy more gear then folding them into `R_G` would bias the effort expansion up.

Two things none of these arms fix, stated so they are not mistaken for solved. First, the statistically correct way to recover all 1,334 dropped interviews is a **censored likelihood**, treating an incomplete trip's observed catch as a lower bound (`neg_binomial_2_lccdf`) rather than an observation; that is a Stan change with a validation burden, it cannot be evaluated by a design-based arm, and it sits in the backlog. Second, excluding incomplete trips does not remove selection bias, it sharpens it: field protocol already favours crabbers who have finished, so early leavers are over-represented in what remains. There is no diagnostic for that and none on any tier.

**Config levers.** Two experiment toggles default to production behavior and are documented in the driver's `params`: `collapse_mu_hier` (default off) collapses the single-cell mu-hierarchy per population for the funnel investigation, and `ar_force` (default null) forces a population's AR resolution. Both leave the default posterior unchanged.

------------------------------------------------------------------------

### 12. Reproducibility

The Stan fits take a fixed RNG seed (`bss_seed`, default 20260619, set in `run_config.R`), passed to `rstan::stan()`. rstan seeds each chain from `bss_seed + chain_id`, so the chains still differ (R-hat remains meaningful) while run-to-run variation is removed. Package and Stan versions and the seed are written to `session_info.txt` with each output set, and the resolved run configuration is written to `run_parameters.txt`. Because every user-selectable toggle now lives in the single `run_config.R` file (the single source of truth applied to the driver as an override), a run is fully specified by that one file plus the input data and the recorded session, whether it was launched through `run_estimation.R` or knit standalone. Change `bss_seed` only if a pathological seed is ever suspected. Expected runtime is 3 to 6 hours on a 4-core machine.

------------------------------------------------------------------------

### 13. Scope: when this method applies, and when it must be re-derived

Method v1.0 is built for one fishery under one sampling design. It can be re-run season after season as long as the following hold. Where one breaks, the method must be revisited, not merely re-fed.

**Assumptions that allow a straight re-run:**

- **Same location.** Westport / Grays Harbor access points (docks Floats 17-21, the jetty, beaches, the boat launch, the marina). The gear-per-crabber prior `R_G`, the gear-per-boat-group prior `R_G_boat`, and the effective-day-length regression are all calibrated to this site.
- **Same input streams, same schema.** The four input files in the same form (Section 7).
- **Same sampling design.** Instantaneous effort counts, dockside interviews, the commercial tally, and I/E surveys, collected as in 2024-25. The 2024-25 protocol of three randomized effort counts per day is the design the effort expansion assumes (it measures mean daily effort).
- **Same sub-season structure.** Pot closure Sep 16 to Nov 30 (non-pot gear only); all-gear Dec 1 to Sep 15, tied to the pot closure. The closure window is set explicitly via `pot_closure_start` / `pot_closure_end`.

**Conditions that require re-derivation, not just new data:**

- **A different port.** `R_G`, `R_G_boat`, and the L_effective regression would have to be re-estimated from that port's I/E and interview data; the access-point structure differs.
- **A change in the effort-count protocol.** Reverting to a single peak-time count per day measures a different quantity (peak, not mean daily effort) and would bias the effort level high. Mixing protocols across years is a genuine confound, addressable only with a protocol fixed effect and a peak-to-mean calibration (the multi-year question in Section 17).
- **A structural change in participation.** For example, a change in how commercial/charter vessels participate pre-season, or the opening of a new major access point (a jetty effort count, currently absent), would change what the components represent.
- **A season with large sampling gaps coinciding with anomalous weather.** The routine model interpolates gaps with its time-series process and deliberately excludes weather (Section 17). A season with extended unsampled stretches under unusual conditions is the one case where the shelved parsimonious weather-effort contingency (Section 17) should be considered, evaluated by leave-one-week-out block cross-validation.

---

## PART III: TECHNICAL REFERENCE

------------------------------------------------------------------------

### 14. Model specification (`crab_bss_pooled.stan`)

#### 14.1 Effort process

```text
log(lambda_E[d]) = mu_E + omega_E[period(d)] + B1 * w[d] + B2 * holiday[d] + X_open[d] . B_open
```

`X_open` is a `D x K_open` matrix of other-fishery opener indicators and `B_open` their log-additive effects (2026-08-25). `K_open = 0` is the production default and removes the term entirely; it generalizes the retired razor-dig-only `B3`. Section 15 covers how columns are selected and the one interaction that must be respected before switching a boat opener on.

**Weekend definition (changed 2026-08-25).** `w[d]` marks Saturday and Sunday. Friday was previously included and is not a weekend day in this fishery: within month, paired across the twelve months carrying both day types, Friday runs at 0.43x Saturday for the shore gear count (t = -8.88, p = 2.4e-06) and 0.60x for boat trailers (t = -3.81, p = 0.003), while sitting statistically on top of the Monday-Thursday mean (1.21x, p = 0.10 shore; 0.99x, p = 0.96 boat). Pooling it dragged the fitted weekend multiplier from 2.34x to 1.74x (shore) and 1.81x to 1.42x (boat) and left 15% more residual variance (shore) for the AR to absorb. A holiday sets both `w` and `holiday`, so a holiday multiplier is the combined effect.

The temporal deviation `omega_E` evolves as an AR(1) process:

```text
omega_E[p] = phi_E * omega_E[p-1] + sigma_eps_E * epsilon[p-1]
```

where `period(d)` maps day `d` to its AR period index. At daily resolution `period(d) = d` and the number of periods `P_n = D`; at weekly or monthly resolution `period(d)` maps to the week or month index and `P_n` is the number of weeks or months. Innovations `epsilon` are standard normal (non-centered parameterization for efficient HMC; Papaspiliopoulos et al. 2007). The AR(1) initial state is non-centered: `omega_E_0` is a raw standard normal scaled in the transformed-parameters block by the stationary standard deviation `sigma_eps_E / sqrt(1 - phi_E^2)`, so the process starts from its stationary distribution without a centered funnel (Harvey 1989; Betancourt and Girolami 2015).

**Adaptive temporal resolution.** The AR resolution is selected automatically per fit from effort-data density:

| Resolution | Condition | Rationale |
|---|---|---|
| Daily | >= 25% of days sampled AND >= 20 effort days | Dense data supports day-level smoothing with proper uncertainty scaling by distance from the nearest observation (Staton et al. 2017) |
| Weekly | >= 1.5 effort obs per week AND >= 3 weeks | Moderate data; weekly states smooth 3-5 day gaps without under-identifying the AR |
| Monthly | Fallback for sparse data | Few AR parameters; robust with limited observations |

This applies the finest resolution the data can identify and falls back gracefully when it cannot (Conn 2002; Sullivan 2003).

**Escalation ladder (`ar_escalate`, added 2026-08-25, ships OFF).** The coverage rule above answers "can this series identify a daily process at all". It does not answer "does the sampler survive it", and that second question was previously answered by hand: the per-population caps in `ar_max_resolution` are an empirical finding from earlier runs frozen into config, and a fit that failed its convergence gate was demoted straight to the Point Estimator rather than retried anywhere it might succeed. With `ar_escalate = TRUE` each component instead starts at the finest rung (daily, ignoring the caps), is put through the **same** convergence gate that decides PE-vs-BSS, and on failure is refit one rung coarser, stopping at the first rung that passes. Every component then reports at the finest resolution its own sampler behaviour supports, and `ar_escalation_log.csv` records each attempt with the verdict that triggered the next one. If no rung passes, the component falls back to PE exactly as before.

The cost is real and is why it ships off: every rung is a multi-hour MCMC fit, and on the 2024-25 configuration the two capped components burn their known-bad rungs (shore pot-closure funnels at daily with about 1,165 divergences; the boat diverged on roughly 100% of iterations at daily) before settling where the caps already put them. Expect roughly two to three times the wall clock. `ar_escalate_respect_cap = TRUE` is the cheap variant: it escalates from the capped rung instead of the top, costing no extra fits when the caps are right but unable to discover that a cap is too coarse. Rungs that would duplicate a finer rung's period count, or leave fewer than three AR periods, are dropped from the ladder.

The data-driven choice is additionally capped per population via `ar_max_resolution`: the boat fit is capped at monthly regardless of coverage, because the trailer-count series cannot identify a daily latent process even when its coverage exceeds the daily threshold (coverage counts how many days carry an observation, not how strongly each observation constrains the latent process). Shore is uncapped at daily, where it converges with n_eff above 2000; the thin shore pot-closure fit, however, fails convergence at daily AR (Run 1), and the top open fix is to adopt a coarser shore AR, as the gear-resolved track does. An `ar_force` parameter can override both the data-driven rule and the cap for a single population, used for the boat daily-vs-weekly resolution experiment, which concluded monthly (Sections 16 and 19); it defaults to `NULL` (production behavior).

#### 14.2 CPUE process

```text
log(lambda_C[d]) = mu_C + omega_C[period(d)] + B1_C * w[d] + B2_C * holiday[d]
```

`B1_C` allows weekend CPUE to differ from weekday CPUE, motivated by evidence that weekend/holiday crabber populations at tourist-accessible ports include more novice participants (Thomson 1991; Pollock et al. 1997). In the reference data `B1_C` is about -0.25 to -0.30 for shore crabbers (weekend crabbers catch roughly 21-26% fewer crab per unit effort than weekday regulars), consistent with the novice-dilution hypothesis. `B1_C` is a multiplier on the catch rate, so this ratio is effort-unit-independent; the point value is nonetheless pre-refresh (Section 19). This is a single pooled CPUE process; gear-type catch is apportioned afterward from interview proportions (the gear-resolved model is the alternative that models per-gear CPUE).

#### 14.3 Observation models

- Gear counts (shore): `Gear_I ~ NegBinomial2(lambda_E[d] * R_G, r_E)`
- Trailer counts (boats): `T_I ~ NegBinomial2(lambda_E[d] / R_G_boat, r_E)` (POOL-1; lambda_E is gear, lambda_E / R_G_boat is boat groups)
- Interview catch: `c ~ NegBinomial2(lambda_C[d] * h, r_C)`, where the CPUE denominator `h` is set per component by `03_R_functions/bss_effort_spec.R`. As of v7.6 (boat) and v7.7 (shore) both components run on gear-deployments, so `h` = `number_of_gear` (the older forms, crabber-hours for shore and gear-hours for the boat, are the revert options; see Sections 14.5 and 15). The model carries `E_scale` so that the seasonal effort `E = lambda_E * E_scale * L` always shares `h`'s unit (shore: `E_scale = R_G`, converting the crabber-scale `lambda_E` to gear; boat: `lambda_E` is already gear via `R_G_boat`, so `E_scale = 1`).
- I/E crabber-hours: `IE_crabber_hours ~ Lognormal(log(lambda_E[d] * L[d]), sigma_IE)` (the I/E survey measures crabber-hours as ground truth regardless of the chosen CPUE unit)

The negative binomial accommodates the overdispersion typical of recreational trip-level catch (Maunder and Punt 2004).

#### 14.4 Effort overdispersion (marginalized)

Each effort count is negative binomial with shape `r_E`. This was originally written as a Gamma-Poisson mixture with an explicit per-observation latent multiplier `eps_E_H_obs ~ Gamma(r_E, r_E)`. Because the Gamma-Poisson mixture integrates exactly to the negative binomial (Hilbe 2011), the latent multipliers are marginalized analytically and the negative binomial is written directly. The change is inference-preserving (the marginal likelihood is identical), removes a high-dimensional centered latent block from the sampler, and makes the model block consistent with the `log_lik` block. The data field `n_effort_obs` is retained as an unused field to keep the R prep interface stable.

#### 14.5 I/E integration and effective day length

On I/E survey days, the observed I/E quantity enters as a direct lognormal observation of `lambda_E * L`, a second independent constraint on the latent effort state that bypasses `R_G` and day-length assumptions and calibrates the gear-count pathway against the I/E ground truth (Robson 1991; Pollock et al. 1994).

**Which I/E quantity is observed depends on the effort unit (corrected 2026-08-25).** The predicted mean is `lambda_E * L`, so the observation has to carry the same unit as that product:

| Shore effort unit | `L` is | predicted `lambda_E * L` | observation |
|---|---|---|---|
| gear-deployments (production) | `tau_shore`, a turnover (~1.7) | crabber **trips** | crabber arrivals (`ie_trips`) |
| crabber-hours / gear-hours | `L_effective`, hours (~5.3) | crabber-**hours** | `ie_crabber_hours` |

This was a live defect between v7.7 and 2026-08-25. The shore move to gear-deployments replaced `L` with the turnover but left the observation as crabber-hours, so the model was comparing crabber-hours against a predicted crabber-trip count. On a typical Float 20 survey day that is roughly 331 observed against roughly 80 predicted, a four-fold scale mismatch the gear-count stream then had to absorb, and it is the most plausible mechanism behind the otherwise unexplained shore all-gear `sigma_IE` of about 1.07 (backlog GR-9). The boat stream had already been put on the matching pair (boat trips against groups x turnover) by F2. `bss_effort_spec()` now owns the pairing for both, so the two cannot drift again; `ie_shore_obs_unit = "crabber_hours"` reproduces the old behaviour for comparison runs. **This changes the shore posterior and must be confirmed by a run.** When no I/E data are available (`IE_n = 0`), the I/E likelihood contributes nothing and the effort and catch posteriors are unchanged. The prior on the I/E scale `sigma_IE ~ exponential(5)` is applied unconditionally (not only inside the `IE_n > 0` branch), so that with no I/E data `sigma_IE` is still proper rather than an improper flat direction; because `sigma_IE` is decoupled from effort and catch, this leaves those posteriors unchanged.

When `estimate_L = 1`, the daily effort-expansion factor `L[d]` is a parameter with a non-centered lognormal prior:

```text
L[d] = L_data[d] * exp(L_sigma[d] * L_raw[d]),    L_raw ~ Normal(0, 1)
```

What `L_data` represents depends on the effort unit, which `03_R_functions/bss_effort_spec.R` sets per component:

- **Production (gear-deployments): `L_data` = the deployment turnover `tau`** (trips per gear-slot per day), prior-centered on the `tau_*` values in `run_config.R`, with `tau_shore` about 1.7 (shore, since v7.7) and `tau_boat` about 1.2 (boat, since v7.6). This replaced the boat's old flat 24-hour soak (`L = 24`) and shore's effective-day-length-in-hours expansion, because catch does not scale with soak time for pots (Sections 11, 15).
- **Crabber-hours revert (shore only): `L_data` = the effective day length in hours, `L_mu`,** from a regression of log effective day length on day-of-year (quadratic) and day type, fit from the I/E data:

```text
log(L_effective) = b0 + b1 * yday + b2 * yday^2 + b3 * weekend + e
```

The quadratic captures the seasonal arc. Effective day length at the docks averages about 3.5 to 5.5 hours, substantially shorter than civil twilight (9 to 16 hours), because crabbers rotate through the dock rather than occupying it all day.

**How the two I/E-derived quantities relate, and which one production uses.** Both come from the same 15-minute presence series on a survey day:

```text
L_effective = crabber-hours / peak crabbers present     (mean ~5.26 h over the WDF20 days)
turnover    = crabber arrivals / peak crabbers present  (mean 1.72 over 30 WDF20 days)
```

Their ratio is the implied mean trip length, `5.26 / 1.72 = 3.07` hours, against an interview-reported mean trip length of 3.23 hours: two independent measurements 5% apart, which is the method's free internal consistency check. Since v7.7 the production expansion is `E = lambda_E * R_G * tau_shore`, so **the turnover is what the estimate rests on and the day length in hours is a diagnostic**. `L_effective` is still computed every run and sets the `day_length` column used by the diagnostics and the civil-twilight comparison, but it leaves the estimation path unless `shore_effort_unit` is set back to a time unit. The `bss_L_effective_*.csv` columns were renamed accordingly on 2026-08-25 (`L_prior_center` plus an explicit `L_unit`, with the hours regression kept in its own column), because the old `L_prior_mu` column reported the hours regression even when the fit's `L` was a turnover.

The known residual on this measurement is coverage, not method: only six shore I/E surveys fall inside the season and they skew toward low-effort days (mean percentile 0.27, rank-test p = 0.06), against the roughly 40 a season Pollock et al. (1997) would imply. That is the cheapest precision left on the table and it is a field-plan item, not a code item.

On I/E days `L[d]` is further constrained by the I/E likelihood; on other days it is informed by its prior and its uncertainty propagates into effort and catch (Pollock et al. 1994; Hartill et al. 2012). The I/E crabber-hour observations continue to calibrate the shore effort level under either unit. See Sections 11, 15, 16, and 19.

#### 14.6 Key parameters and priors

| Parameter | Description | Prior |
|---|---|---|
| B1 | Weekend effort multiplier (log) | Normal(0, 1) |
| B2 | Holiday effort multiplier (log) | Normal(0, 1) |
| B1_C | Weekend CPUE effect (log) | Normal(0, 1) |
| B2_C | Holiday CPUE effect (log) | Normal(0, sigma) |
| R_G | Gear per crabber | Lognormal(log(R_G_empirical), 0.3), data-driven |
| R_G_boat | Gear per boat group | Lognormal(log 4, 0.5) |
| phi_E, phi_C | AR(1) autocorrelation | Beta(2, 2) rescaled to [-1, 1] |
| r_E, r_C | Overdispersion | Half-Cauchy(0, 1) |
| sigma_IE | I/E measurement error (log) | Exponential(5) |
| L[d] | Deployment turnover `tau` (production); effective day length under a time unit | Lognormal centred on `tau_*_prior_mu`, log-SD `tau_*_prior_sigma` |
| B_open | Other-fishery opener effort covariates (log) | Normal(0, 1); absent when `K_open = 0` |
| kappa_OSP | OSP-to-trailer scale (within-day boat turnover) | Lognormal(log 3.0, 0.3); inert when `osp_scale_is_tau = 1` |
| theta (`f_theta`) | Combo-trip share among non-crab-labelled boats | Beta(`combo_share * kappa`, ...) in an OSP-informed stratum; otherwise the ordinary f prior |
| f_lower | OSP-observed crab-only share (lower bound on f) | Beta(1, 1), updated by the OSP Binomial |
| f | Crabbing fraction = `f_lower + (1 - f_lower) * theta` | derived; see Section 16 |

Prior rationale: `R_G` is centered on the empirical gear-per-crabber ratio in the relevant population by sub-season, eliminating prior-posterior conflict; `R_G_boat` Lognormal(log 4, 0.5) is centered on ~4 gear per boat group (POOL-1; replaces the old R_T Beta(5, 1), which was pinned at 1 by a degenerate bernoulli term); the Half-Cauchy(0, 1) variance priors are weakly informative (Gelman 2006).

#### 14.7 Generated quantities

- `C_expected[d] = lambda_E[d] * L[d] * lambda_C[d]`: the posterior expected daily catch, E[C | data], the quantity used for harvest estimation.
- `C[d] = Poisson_rng(C_expected[d])`: a predictive draw including Poisson sampling noise, reported separately for prediction intervals.

For seasonal totals the Poisson noise largely averages out, so the two are similar; for daily or monthly breakdowns the difference can be material. Pointwise `log_lik` for the gear, trailer, and catch streams is also produced, enabling PSIS-LOO.

------------------------------------------------------------------------

### 15. Design decisions and their rationale

- **Two sub-seasons** are estimated independently because the pot closure creates a structural break in both effort and catch rates; pooling across it would blur two different regimes.
- **Gear-deployments for both components (v7.6 boat, v7.7 shore).** Effort is denominated in gear-deployments, not soak-hours or crabber-hours. Earlier versions used gear-hours for the boat on the reasoning that soak time is the fishing-time measure for crab gear, but the saturation diagnostic overturned that: catch per gear-hour falls roughly 43-fold across soak durations while catch per deployment is nearly flat (catch per gear scales as about soak-hours^0.13 for the boat and ^0.22 for shore), so any time-denominated unit violates the likelihood's proportionality assumption (`beta_h = 1`). The deployment is the unit on which CPUE is a stable rate (about 4 to 7 crab per deployment), which is what keeps harvest = effort x CPUE unbiased. Crabber-hours for shore failed the same test (`beta_h = 0.57`), so v7.7 moves shore onto deployments too. See Sections 11, 16, and 19.
  - **Where the unit reaches, audited 2026-08-25.** The unit is not only the CPUE denominator; it has to hold in every place effort is formed, or the diagnostics and the monthly tables quietly contradict the totals. Four residues from before the v7.6/v7.7 move were found and fixed in the same pass: the shore I/E observation (Section 14.5, the substantive one), the shore branch of the monthly PE effort share (`pe_monthly_effort_share.R`, which still multiplied by the seasonal day-length curve and therefore re-weighted the monthly split toward long-day summer months; totals were unaffected because the share is normalised), the fishery-opener spillover diagnostic's CPUE (measured per crabber-hour, the denominator this pipeline's own linearity test rejects), and the gear-resolved daily-combined series (which applied the shore crabber-hours formula to the boat as well). The BSS effort plot axis and the `bss_L_effective_*.csv` column names were corrected at the same time.
  - **Deployments are the right unit; they are not a perfect one.** Every component covers `beta_h = 1` on deployments and none does on either time unit, so the unit choice is settled. But catch is not literally flat in soak time: the saturation exponent is 0.22 for shore all-gear, 0.27 for the pot closure and 0.25 for the boat, and all three intervals exclude zero as well as one. All three still raise the saturation flag on every production run. A soak-time term would close the remaining gap; the deployment unit gets most of the way.
- **Other-fishery opener effort covariates (`opener_covariate_mode`, added 2026-08-25, ships `"off"`).** Any opener in the consolidated calendar can enter either population's effort model as a log-additive day covariate, replacing the razor-dig-only `B3`. In `"auto"` a candidate is included only if the spillover diagnostic's day-type + month adjusted effect on that population's effort series clears `opener_auto_p` **after a multiplicity adjustment across the eight effort tests** (Benjamini-Hochberg by default), and only if the indicator is identifiable inside the fit's own window. Two cautions are structural, not stylistic. First, selecting on a p-value from the data that then fits the term inflates the coefficient and understates its interval, so `"auto"` is a screen and the arbiter is the effort-stream `elpd_loo` against an opener-free run. Second, the project already has the worked example: razor-dig on shore effort came in at p = 0.045, was fitted as `B3` in Run 3, and gave no predictive gain on any stream; under the default adjustment it does not survive the screen, while MA2 halibut on boat trailers (p = 3.6e-06) does. **Before switching on a boat opener, read the interaction with `f` in Section 16.**
- **L_effective from I/E, not civil twilight,** on the crabber-hours revert unit: the dock activity curve is peaked (crabbers rotate through), so civil twilight overstates the time gear is actively fished by roughly a factor of two. On the production gear-deployment unit the shore daily expansion is the turnover `tau_shore` instead (Section 14.5).
- **A weekend CPUE effect (B1_C)** because weekend crabber composition differs (more novices), which the catch-rate process should be allowed to reflect.
- **Expected catch, not a predictive draw, as the headline** because harvest estimation wants E[C | data], the estimation quantity, not a single noisy prediction.

------------------------------------------------------------------------

### 16. Limitations and the private-boat caveat

**General limitations.**

- I/E coverage is limited; expanding toward roughly 40 days per season would improve the L_effective regression (Pollock et al. 1997 recommend at least 3 I/E days per month by day-type stratum).
- The L_effective regression is a quadratic in day-of-year; with more data a GAM could capture non-monotonic patterns.
- Gear-type breakdowns from this pooled model are approximate (proportional allocation, not modeled). The gear-resolved model is the alternative when modeled gear catch with uncertainty is needed.
- There are no jetty effort counts; beach crabbing is unmeasured within the pooled shore count.
- The weekend CPUE effect is constant across the season; a time-varying effect may be warranted if tourist composition shifts seasonally.
- The adaptive AR selection is rule-based; a formal LOO/WAIC model comparison could provide principled resolution selection (Vehtari et al. 2017).

**The private-boat pot closure, and the interview floor (changed 2026-08-25).** The boat pot-closure component carries 17 interviews. Under the previous floor of 20 it was never fitted: it entered the port total as a fixed Point Estimate with no interval, which together with the commercial/charter census meant roughly 18% of the port total carried no uncertainty at all. `bss_min_interviews` is now 15, so the component attempts a BSS. Two things to hold about that. First, it is surgical: shore all-gear (2,741), shore pot-closure (856) and boat all-gear (145) are far above either value, so no other component's behaviour changes. Second, 17 interviews over 76 days is genuinely thin and the fit may well fail the convergence gate and fall back to PE anyway — which is a better outcome than before, because the failure is then a measured gate verdict with per-fit diagnostics attached rather than a threshold decision taken before any data were looked at. A second floor, `bss_min_interviews_fitted`, was added at the same time and is applied **after** the pipeline drops zero-denominator rows and incomplete trips; the pre-existing guard counts unfiltered interviews and so reads looser than the count the likelihood actually sees, which was harmless with hundreds to spare and is not harmless at 15. Both floors now live only in `run_config.R`; they previously sat in each driver's `params_model`, which is merged **on top of** `run_config` and therefore silently overrode it.

**The private boat.** The private-boat all-gear fit rests on a thin, weakly informative trailer-count series, and is the component most prone to wide posteriors and to PE-vs-BSS disagreement. A sequence of fixes was applied to make the fit converge: dedicated sampler tuning, a per-population AR cap to weekly, non-centering of the AR initial state, an unconditional `sigma_IE` prior (the boat has no I/E data, so this removed an improper flat direction that had inflated divergences), and finally the scale-aware convergence gate. Under that gate the boat all-gear component is reported on its BSS posterior, with a wide 95% interval and a catch CV around 27%, rather than the narrower PE point. The wide interval is not a defect to be hidden by substituting the PE point; it is the effort and CPUE data honestly reporting their own uncertainty for a component identified by only a handful of interviews per month.

The boat effort unit was corrected in v7.6 (POOL-1 + POOL-3): the boat now runs on gear-deployments, matching the gear-resolved model. The degenerate trailer expansion (`R_T` pinned at about 1) is replaced by `T_I ~ NB2(lambda_E / R_G_boat, r_E)`, with gear-per-group learned as `R_G_boat` from interviews, and the CPUE denominator becomes `number_of_gear` (deployments) with `L = tau_boat`, replacing the old expansion that assigned every trailer-day a full gear-per-group times a flat 24-hour soak. This is justified because catch is not linear in soak time for pots (the saturation diagnostic gives catch per gear scaling as about soak-hours^0.25), so any time-denominated effort unit is invalid. On this scale the gear-resolved model reports a boat total about 25% below the old pooled gear-hours figure (43,314 vs 56,266). The 2026-07-11 run is now complete, and it shows the pooled boat did NOT follow: on gear-deployments the pooled boat holds at 54,481 (median), essentially unchanged (about -3.2%), because within one pipeline the unit change is catch-neutral (rescaling `lambda_C` and the effort process inversely conserves `E x CPUE`). The correction is therefore an interpretive fix (it removes the `R_T`-pinned-at-1 degeneracy and makes the saturation/linearity diagnostics valid), not a boat-harvest change. Because the boat PE and BSS now share the `bss_effort_spec.R` unit, the remaining PE-vs-BSS gap is a genuine model-vs-design disagreement rather than a unit artifact. The reconciling run has now been done. On the weekly boat AR the pooled boat (54,481) sat above both the gear-resolved boat (43,314) and the design anchor (PE effort times interview ratio-of-sums, about 40,000); on **monthly** boat AR (`05_output/20260713/pooled-CPUE-run1`) the pooled boat falls to **43,180** (median) and matches the independent gear-resolved boat on all three quantities: catch 43,180 vs 43,314, effort 14,707 vs 14,805 deployments, model CPUE ~2.95 vs 2.93 crab/deployment. The weekly-vs-monthly AR resolution, not the effort unit, was the entire source of the pooled-vs-gear-resolved boat gap: monthly AR cannot chase sparse-week effort spikes, so the sparse-month over-imputation collapses (June boat BSS/PE fell from 8.7x to 4.1x). Monthly should be adopted as the production boat AR cap. Run 1's empty-stratum pooled PE fallback raised the boat PE to 37,638 (implied CPUE 3.02), closing the boat PE-vs-BSS gap to about 15% (1.15x, from 1.9x); the reported boat is the BSS at 43,180, cross-checked against the independent gear-resolved boat (43,314).

The data-thinness half of the caveat is now being addressed by the OSP boat-count work (branch `OSP-boat-count-incorporation`, 2026-07-31), which supplies exactly the more informative effort series the caveat called for, rather than further parameter surgery on the shared model.

**OSP boat counts (the second effort stream).** The Ocean Shores Patrol Westport Boat Launch (WBL) daily private-boat total (`04_input_files/WBL_boat_counts.xlsx`, read by `fetch_osp_boat_counts`, toggle `use_osp_boat_counts`, production ON) is added as a SECOND boat effort observation on the same latent `lambda_E`, scaled by an OSP within-day turnover `kappa_OSP`: `OSP_I[i] ~ neg_binomial_2((lambda_E / R_G_boat) * kappa_OSP, r_OSP)`. It contributes on its 148 operating days (mid-March to mid-October); on non-operating days the model runs on trailer counts alone, so the pre-OSP behavior is the emergent fallback with no special-casing. Validation (Step 2) shows OSP raises the boat effort about 10% and tightens it modestly (boat-effort CV 12.2% to 11.2%), with the larger value coming from 87 "OSP only" days that carry a direct effort observation the trailer series lacked.

**Crabbing fraction f (all-boat to crab-effort conversion).** Trailer counts AND OSP counts are ALL private boats at the launch, not just crabbers, so at f = 1 the model counts every launched boat as a crabber and biases the boat catch high in a mixed salmon/tuna/crab fishery. `use_crab_fraction` (production ON) multiplies the boat effort and catch by f in the BSS boat generated quantities ONLY, decoupled from sampling, so it cannot affect convergence and leaves CPUE invariant (validation Step 3/Step 4 confirm the boat scales exactly linearly in f: 0.2/0.5/1.0 give 8.7k/21.6k/43.2k, with f = 1 reproducing the pre-f baseline). f is a Beta prior centered on `crab_fraction_set` (production 0.3), to be updated by the crab-creel ingress/egress crab-vs-total classification (`ingress_egress.xlsx` columns `boats_crabbing` / `boats_total`) once the WBL egress pilot lands. **0.3 is a placeholder chosen by Matt, not a measured estimate**, and it is the single largest lever on the boat harvest; replace it with the pilot's f_hat.

**OSP-informs-tau (`osp_scale_is_tau`, production ON) and the turnover resolution.** GR-12 (boat catch rests on a `tau_boat` anchored on two winter I/E days) is closed by letting the dense OSP series identify the boat turnover directly: with `osp_scale_is_tau = TRUE` the OSP mean uses `L` (= `tau_boat`) in place of the free `kappa_OSP`. The Phase 0 overlap and the model agree that the OSP/trailer within-day turnover is about 2.7 to 3.0 (validation posterior `kappa_OSP` = 3.15, 95% [2.50, 3.91]), while the old `tau_boat` prior was about 1.2. **Resolved 2026-07-31: the crab-creel trailer count is an instantaneous snapshot**, so that 2.7 to 3.0 is real within-day turnover and `tau_boat` about 1.2 was a roughly 2x under-count of boat effort; turning the toggle on corrects it and roughly doubles the boat (validation Step 6 boat about 27.7k vs Step 3b about 13.8k). One caveat is retained honestly: the exact multiplier assumes the trailer snapshot is timed representatively (not fixed at a daily peak); a peak-timed snapshot would understate turnover, so 2.7 is a floor rather than a ceiling.

**OSP crab-only counts as a hard lower bound on f (`use_osp_crab_lower`, added 2026-08-25).** OSP can report, per day, how many boats were labelled as crabbing **only**. It cannot see combo trips: a boat that crabs and also fishes another fishery is labelled by the other fishery. That makes the OSP crab-only share a **lower bound** on f, never an estimate of it, and treating it as an estimate would bias the boat harvest down by exactly the combo-trip rate, on precisely the halibut and tuna days where that rate is highest. The bound is therefore imposed by construction rather than by convention. Per f stratum k:

```text
f_lower[k] ~ Beta(1,1)              osp_crab_only[k] ~ Binomial(osp_total[k], f_lower[k])
f[k]       = f_lower[k] + (1 - f_lower[k]) * theta[k]
```

`theta[k]` is the share of the not-crab-labelled boats that were also crabbing. State the identifiability plainly: **OSP alone identifies `f_lower`, never `theta`** — the counts cannot distinguish "few crab boats" from "many combo trips". `theta` rests on its prior (`crab_fraction_combo_share`, a placeholder with exactly the standing that f = 0.3 had before the pilot) until the WPT/WBL egress classification covers the same stratum, at which point the egress Binomial `n_crab ~ Binomial(n_total, f)` pins f and therefore theta. That is the operational argument for scheduling egress classification days on days OSP is also in port, rather than treating the two as independent efforts.

The degradation is exact, not approximate. In a stratum with no OSP classification, `f_lower` is pinned to 0, `f = theta`, and the R side hands Stan the ordinary f prior for theta instead of the combo prior, so the posterior for f is identical to the pre-2026-08-25 model — "defaults back to interview data when OSP is not in port" is a property of the parameterization, not a special case in the code. Both new likelihood terms sit on **observed boat counts**, never on the latent effort, so f still enters generated quantities only: the boat remains exactly linear in f and the model CPUE remains invariant, preserving the validated Phase 2/3 behaviour.

**`use_osp_crab_lower` ships FALSE.** Every other behaviour-changing feature in the 2026-08-25 batch is opt-in and this one has to be as well: with a TRUE default the boat harvest would move on the same run that first read the new column, with no baseline to compare against, and at the production `crab_fraction_strata = "none"` the whole window collapses to one stratum, so `f` becomes `p_osp + (1 - p_osp) * crab_fraction_combo_share` and a large part of the move would come from the combo-share placeholder rather than from OSP's data. Turn it on deliberately, in its own run. Beyond the toggle it is also inert until the crab-only column exists in `WBL_boat_counts.xlsx` (`osp_crab_only_col`), and a stratum with fewer than `crab_fraction_osp_min_obs` classified boats does not bind.

**The opener-covariate interaction (read before switching on a boat opener).** Boat effort counts all private boats, so an MA2-halibut effort covariate legitimately describes that series: halibut days add about 32 trailers to a ramp averaging 7. But those boats are not crabbing. A **constant** f converts the fitted surge into crab effort at the same rate as a closed-fishery day, so fitting the surge more faithfully while multiplying it by a flat f makes the boat catch **more** biased on those days, not less. The two features must move together: pair any boat opener covariate with `crab_fraction_strata = "opener"` or `"month_opener"` (which key off `opener_f_flag`, default `ma2_halibut_open`). Both drivers warn when a boat opener is active and f is not opener-aware.

**Production configuration (adopted 2026-07-31; `use_osp_crab_lower` added 2026-08-25 and inert without the crab-only column):** `use_osp_boat_counts = TRUE`, `use_crab_fraction = TRUE` with `crab_fraction_set = 0.3`, `osp_scale_is_tau = TRUE`. Confirmed by the production run `05_output/20260804/pooled-CPUE-boat-count-validation-run` (all fits converged): the private-boat harvest is **27,684** Dungeness kept (BSS median; 95% CI 11,432 to 53,617) and the **port total is 67,312** (95% CI 50,601 to 93,461; PE port 44,810). That is down from about 43k at f = 1 via the 0.3 crab share, then up about 2x via the OSP turnover correction, and the run reproduces the pre-merge validation batch (Step 6) fit-for-fit. The full evidence, the run matrix, and the two open decisions (the pilot f_hat, and confirming the snapshot timing) are in `development_notes/osp-validation-review-2026-07-31.md`, with the change log in the development-history document.

------------------------------------------------------------------------

### 17. Weather and tide covariates: evaluated and excluded

A weather-and-tide covariate module (`06_diagnostics/`) was built and run on the 2024-25 season to test whether weather improves the estimate. The conclusion, documented in full in `WEATHER_COVARIATE_ANALYSIS.md`, is that **covariates are excluded for all three components** under a pre-committed 4.0-SE PSIS-LOO improvement margin. No component cleared the margin; shore all-gear was a tie, and shore ring-net and the boat were meaningfully worse out-of-sample with covariates.

The instructive findings:

- **Weather drives effort, not CPUE.** Boat effort is suppressed by wave height and rain and raised by temperature and tide range; shore effort is suppressed by rain. CPUE is essentially weather-flat. These are descriptive effort-dynamics results worth reporting (and they bear on sampling design, since effort is predictably low on rough or rainy days), but they do not improve harvest prediction.
- **False precision.** The covariate models produced narrower credible intervals while predicting held-out data worse. Narrower-but-worse is the signature of overfitting, which is exactly why selection here uses LOO and not interval width.
- **Why significant effects do not help.** On a sampled day, the effort count already encodes the weather effect (a rough day shows a low count whether or not weather is in the model). Weather could only add value by improving interpolation on unsampled days, and the AR(1) process already does that interpolation; on routine data, weather does not beat it.
- **The one exception, kept on the shelf.** If a season has a stretch with no effort counts (a true sampling gap), the AR's redundancy with the counts disappears, and a parsimonious weather-effort model (a few strong drivers, not the full screened set) could correct the AR's naive interpolation, most valuably when the unsampled period's weather is anomalous. The observation-level LOO used here cannot test that scenario (it never holds out a full week); the correct test is leave-one-week-out block cross-validation. Weather is therefore not part of the routine production model, but is worth keeping as a contingency for grounding effort across sampling gaps.

A known reconciliation item: the covariate module's absolute boat effort is well below the main pipeline at identical CPUE, because the module predates the current gear-deployment expansion (v7.6 boat, v7.7 shore). This does not affect the covariate decision (the LOO is computed on the latent fit and is invariant to the post-hoc expansion) or the production estimate (the main pipeline is authoritative), but the module's boat magnitudes should not be read as harvest until reconciled.

------------------------------------------------------------------------

### 18. Glossary

| Term | Meaning |
|---|---|
| BSS | Bayesian State-Space model |
| PE | Point Estimator |
| CPUE | Catch Per Unit Effort |
| AR(1) | First-order autoregressive process |
| P_n | Number of AR periods (= D for daily, fewer for weekly/monthly) |
| period(d) | Mapping from day d to its AR period index |
| R_G | Gear-per-crabber ratio |
| R_G_boat | Gear per boat group |
| B1 / B2 | Weekend / holiday effort multipliers (log) |
| B1_C | Weekend CPUE multiplier; exp(B1_C) = weekend/weekday CPUE ratio |
| C_expected | Expected daily catch (no Poisson noise); E[C \| data] |
| C | Predictive daily catch draw (includes Poisson noise) |
| L_effective | Effective day length (hours): I/E crabber-hours divided by peak crabbers present. The daily effort expansion `L` when a component runs on the crabber-hours unit; on the production gear-deployment unit it is superseded by `tau` (Sections 14.5, 15) |
| L_mu, L_sigma | Regression-predicted median and uncertainty for L_effective |
| tau (tau_shore, tau_boat) | Gear-deployment turnover: trips per gear-slot per day; the daily effort expansion `L` on the production gear-deployment scale (tau_shore about 1.7, tau_boat about 1.2) |
| gear-deployment | The production effort unit for both components: one piece of gear a crabber had in the water on that trip, `h` = `number_of_gear`, replacing crabber-hours (shore) and gear-hours (boat) as of v7.7 / v7.6. **Not** a "pot lift": repeat checks of the same gear slot are not extra deployments; slot re-use across the day is carried by the turnover `tau` |
| turnover (`tau`) | Trips per gear-slot per day (arrivals / peak present from the I/E surveys; `tau_shore` about 1.7, `tau_boat` about 1.2). This is the daily expansion factor `L` in production, not a day length |
| L_effective | Crabber-hours divided by peak crabbers present on an I/E day (about 5.3 h). A **diagnostic** since v7.7: production expands on the turnover instead. `L_effective / tau` = 3.07 h implied trip length, against 3.23 h reported in interviews |
| K_open / B_open | The other-fishery opener effort covariates and their log effects; `K_open = 0` (production) removes the term |
| f | Crabbing fraction: the share of private boats that are crabbing. `f = f_lower + (1 - f_lower) * theta` |
| f_lower | The OSP-observed crab-**only** share, a hard lower bound on f (OSP labels combo trips by the non-crab fishery) |
| theta | The share of the not-crab-labelled boats that were also crabbing (the combo trips OSP cannot see) |
| I/E | Ingress/Egress survey |
| sigma_IE | I/E measurement error on the log scale |
| n_eff | Bulk effective sample size |
| R-hat | Rank-normalized split potential scale reduction factor (convergence) |
| PSIS-LOO | Pareto-smoothed importance-sampling leave-one-out cross-validation |
| PIT | Probability integral transform (posterior predictive calibration) |

------------------------------------------------------------------------

### 19. Development history (summary)

**2026-08-25 improvement batch (branch `OSP-boat-count-incorporation`; not yet run).** Eight items plus two defects found during the review. In brief, and in the order they matter to the number:

1. *(defect, inference-changing)* The shore I/E likelihood had been comparing crabber-hours against a predicted crabber-**trip** count ever since the v7.7 shore unit move, a roughly four-fold scale mismatch and the most plausible cause of the unexplained shore `sigma_IE` of about 1.07. The observation now follows the effort unit (Section 14.5). **This moves the shore number.**
2. *(defect)* The shore branch of the monthly PE effort share still used the crabber-hours formula with the seasonal day-length curve, re-weighting the monthly split toward long-day months. Component and port totals were unaffected (the share is normalised); the monthly split was not.
3. Weekend redefined to Saturday and Sunday; Friday is a weekday in this fishery on the season's own data (Section 14.1).
4. `bss_min_interviews` 20 to 15, so the boat pot closure can attempt a BSS instead of entering as an interval-free point, plus a post-filter floor and the removal of the driver-level override that silently beat `run_config` (Section 16).
5. Other-fishery openers generalised into selectable effort covariates with a multiplicity-controlled automatic screen (Section 15), replacing the razor-dig-only term.
6. AR resolution escalation ladder: start at the finest rung and coarsen only on a gate failure, so each component reports at the finest resolution its own sampler behaviour supports (Section 14.1). Ships off; it costs roughly two to three times the wall clock.
7. OSP crab-only counts wired in as a **hard lower bound** on the crabbing fraction f, with the combo-trip share as the separate, egress-identified quantity (Section 16).
8. A four-arm incomplete-trip treatment diagnostic (exclude / gear-only / mean-CPUE imputation / keep). The production estimator is unchanged; the finding worth acting on is that an interrupted trip's fully observed **gear count** is currently discarded along with its truncated catch, and that the PE and the BSS already disagree about this.

Terminology, unit-residue and curated-table fixes rode along; the full list is in `development_notes/PIPELINE_STATUS.md`. **None of this has been confirmed by a run.** Item 1 is the headline change and should be isolated in its own run against the 2026-08-04 production baseline before anything else is switched on.

Method v1.0 is the frozen method label; it was first frozen against pipeline code **v7.4**. The current pipeline code is **v7.9 + Tier-2 batch** (2026-07-13). The method itself is unchanged since v1.0, but several code revisions each move published totals: v7.5 added the incomplete-trip filter, v7.6 moved the boat onto gear-deployments (catch-neutral, as the 2026-07-11 run showed), v7.7 moved the shore BSS onto gear-deployments, and v7.8 moved the shore PE onto the same unit (v7.9 is a config/gate refinement that moves nothing). The 2026-07-11 v7.8 run refreshed the totals; reference numbers below that predate it should be read as pre-refresh. The model began as an adaptation of the WDFW freshwater-creel state-space framework and was hardened over a sequence of versions in response to a 2026-03-31 model critique and an extended convergence-debugging effort focused on the private boat. The arc in one screen:

- **v3-v5:** shared state-space milestones; the gear-resolved track branched at v5.
- **v6.0:** post-critique modeling upgrades (adaptive AR resolution; L_effective as an estimated parameter; the B1_C weekend CPUE effect; direct I/E integration; data-driven R_G; sparse effort overdispersion; expected and predictive catch both reported).
- **v6.1-v6.6:** the convergence gate gained divergence awareness; boat sampler tuning, a per-population AR cap to weekly, and non-centering of the AR initial state addressed boat non-convergence.
- **v6.7-v6.8:** effort overdispersion marginalized to negative binomial (inference-preserving); an unconditional `sigma_IE` prior fixed an improper direction specific to the boat; per-fit model diagnostics, a fixed seed, and session capture were added.
- **v6.9-v6.9.1:** a single-cell scale collapse (B1.7) was attempted and reverted after it hung the shore all-gear fit (the standing lesson: the durable boat fix is a better effort series, not parameter surgery); PPC calibration was hardened; monthly catch by mode was added.
- **v7.0:** the scale-aware convergence gate (impact measured in posterior standard deviations, not as a percentage of level), which moved the boat onto its BSS posterior and made the gate control the selection rather than merely label it; a PE monthly effort-share fix; a PPC extraction fix.
- **v7.1-v7.4:** the effort over-dispersion decomposition diagnostic; an extended set of persisted per-fit outputs (the O-series); pointwise `log_lik` enabling PSIS-LOO on the pooled model; and the `ar_force` experiment toggle (a tight-pin attempt in v7.3 was reverted in v7.4 after it tipped the shore funnel into failure).
- **v7.5 (2026-07-10):** the pooled backlog fixes POOL-2/4/5/6. The R layer was de-duplicated onto the shared gate and AR selector (POOL-6, behavior-preserving); the CPUE effort-unit diagnostics were wired in (POOL-5); a `collapse_mu_hier` lever was added for the funnel investigation (POOL-4, default off); and the incomplete-trip filter was added (POOL-2, default on), which raises the shore estimate, so a re-run is needed to refresh these numbers. The boat-structure items POOL-1 and POOL-3 were held for a validated session because they move the publication boat number. See `development_notes/PIPELINE_STATUS.md`.
- **v7.6 (2026-07-10):** POOL-1 + POOL-3. The private boat is moved onto the gear-deployment scale, matching the gear-resolved model: `R_T` (pinned at ~1) is replaced by `R_G_boat` with `T_I ~ NB2(lambda_E / R_G_boat)` and `Gear_A_boat ~ poisson(R_G_boat)`, and the boat CPUE denominator becomes `number_of_gear` with `L = tau_boat` instead of gear-hours with `L = 24`. This resolves the private-boat effort-unit caveat in Section 16. The 2026-07-11 run showed this is catch-neutral (the boat held at 54,481, about -3.2%, NOT the ~-25% originally predicted toward the gear-resolved boat); it corrects the boat's unit and interpretation but does not move the boat harvest. Shore is unchanged.
- **v7.7 (2026-07-11):** Shore moved onto the gear-deployment scale, so both components now share one effort unit (gear-deployments) with the gear-resolved model. The shore CPUE denominator becomes `number_of_gear` with `E_scale = R_G` and `L = tau_shore` (~1.7 turnover), replacing crabber-hours. This settles the shore half of the effort-unit question (backlog GR-16) using the 2026-07-10 shore LOO comparison (shore all_gear, n = 1649): gear-deployments is the only shore unit whose linearity coefficient covers 1 (`beta_h = 1.05`, 95% CI 0.94 to 1.15, flag off), against crabber-hours (0.57) and gear-hours (0.73), and the only one with no estimator-triad drift (ratio-of-sums 0.87 ~= mean-of-ratios 0.85 ~= model-implied 0.85 crab per deployment). gear-hours had a marginally better catch-stream `elpd_loo` (-3131 vs -3190 for deployments), but that predictive edge comes from the CPUE process absorbing the sub-linearity, which is what biases the season expansion; the choice therefore prioritizes harvest-unbiasedness over marginal predictive fit. There is no Stan change (the v7.6 `effort_scale_gear` / `E_scale` machinery already supports shore), but the shore BSS publication number moves.
- **v7.8 (2026-07-11):** Behavior-preserving code refactor (helper functions extracted to `03_R_functions/`, all user toggles centralized in `run_config.R`) plus the shore-PE completion fix: the pooled `run_pe` shore branch now reads its unit from `bss_effort_spec.R`, so the shore PE runs on gear-deployments to match the shore BSS (previously it was left on crabber-hours). The 2026-07-11 v7.8 run confirmed the fix (shore all-gear PE effort moved from 42,541 to 18,104; BSS unchanged) and refreshed the totals.
- **OSP boat-count incorporation (branch `OSP-boat-count-incorporation`, 2026-07-31; not yet merged to Method v1.0).** Three boat-effort additions, validated by a 14-run batch and detailed in `development_notes/osp-validation-review-2026-07-31.md` and the development-history document: the OSP Westport Boat Launch daily boat-total as a second effort stream on `lambda_E` (scaled by `kappa_OSP`); a crabbing fraction f that converts all-boat effort to crab effort in the boat generated quantities (production 0.3, a placeholder pending the WBL egress pilot; the boat scales exactly linearly in f, CPUE invariant); and `osp_scale_is_tau`, which lets the dense OSP series identify the boat turnover and closes GR-12. Because the crab-creel trailer count was confirmed to be an instantaneous snapshot, the OSP-implied turnover (~2.7 to 3.0; posterior `kappa_OSP` = 3.15) is real and the old `tau_boat` ~1.2 was a ~2x under-count, so production `osp_scale_is_tau = TRUE`. Production private-boat harvest 27,684 Dungeness kept, port total 67,312, confirmed by the run `05_output/20260804/pooled-CPUE-boat-count-validation-run`. Two interview data fixes ride along: non-crabbing interviews (any `number_of_gear = 0` row, regardless of trip status) and gear-tampered interviews (`gear_tampered = 1`) are dropped.
- **v7.9 (2026-07-12):** Config/gate refinement, no estimate change. `run_config.R` is now the base parameter set with each model layering its own tuning on top (`params <- modifyList(run_config, params_model)`), and the per-model AR resolution map moved into `run_config.R`. The pooled divergence-fraction backstop was tightened from 0.15 to 0.05 to match the gear-resolved model (the scale-aware impact test remains the primary gate; the 2026-07-11 fits all sit under 5%, so no gating decision changes).
- **2026-07-12 (later):** Two follow-ons. (a) The boat monthly-AR reconciliation run: on monthly boat AR (`05_output/20260713/pooled-CPUE-run1`) the pooled boat falls from the weekly 54,481 to 43,180, matching the independent gear-resolved boat (43,314) on catch, effort, and CPUE. The AR resolution, not the effort unit, was the source of the pooled-vs-gear-resolved boat gap; monthly should be the production boat AR cap. (b) A report display pass: the driver now renders the incomplete-trip sensitivity, PPC calibration, effort over-dispersion, the CPUE validity triad, per-fit coverage, PSIS-LOO, and top divergence drivers on the page, and curates the wide convergence table. The forward-looking backlog is consolidated in `development_notes/PIPELINE_STATUS.md`.

The full change log, with the per-version rationale, the divergence-diagnostic narrative, and the detailed B1.5 / B1.6 working notes, is in **`BSS-GH-pooled-CPUE-model-development-history.md`**.

------------------------------------------------------------------------

### 20. References

Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo. *arXiv preprint* arXiv:1701.02434.

Betancourt, M. & Girolami, M. (2015). Hamiltonian Monte Carlo for hierarchical models. *In:* Current Trends in Bayesian Methodology with Applications. CRC Press.

Conn, P.B. (2002). Bayesian methods for estimating recreational angler effort, catch rates, and total catch using creel survey data. Ph.D. Dissertation, University of Wisconsin-Madison.

Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models. *Bayesian Analysis*, 1(3), 515-534.

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., & Rubin, D.B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

Hahn, P.K.J., Brooks, L., & Hartill, B.W. (2000). Strategies and procedures for estimating catch and effort in freshwater fisheries. *In:* Inland Fisheries Management in North America (2nd ed.), American Fisheries Society.

Hartill, B.W., Cryer, M., Lyle, J.M., Rees, E.B., Ryan, K.L., Steffe, A.S., Taylor, S.M., West, L., & Wise, B.S. (2012). Scale- and context-dependent selection of recreational harvest estimation methods. *North American Journal of Fisheries Management*, 32(1), 109-123.

Harvey, A.C. (1989). *Forecasting, Structural Time Series Models and the Kalman Filter*. Cambridge University Press.

Hilbe, J.M. (2011). *Negative Binomial Regression* (2nd ed.). Cambridge University Press.

Maunder, M.N. & Punt, A.E. (2004). Standardizing catch and effort data: a review of recent approaches. *Fisheries Research*, 70(2-3), 141-159.

Papaspiliopoulos, O., Roberts, G.O., & Skold, M. (2007). A general framework for the parametrization of hierarchical models. *Statistical Science*, 22(1), 59-73.

Pollock, K.H., Jones, C.M., & Brown, T.L. (1994). *Angler Survey Methods and Their Applications in Fisheries Management*. American Fisheries Society Special Publication 25.

Pollock, K.H., Hoenig, J.M., Jones, C.M., Robson, D.S., & Greene, C.J. (1997). Catch rate estimation for roving and access point surveys. *North American Journal of Fisheries Management*, 17(1), 11-19.

Robson, D.S. (1991). The roving creel survey. *American Fisheries Society Symposium*, 12, 137-148.

Staton, B.A., Catalano, M.J., Connors, B.M., Coggins, L.G., Jones, M.L., Walters, C.J., Fleischman, S.J., & Beardsall, J.W. (2017). Evaluation of methods for spawner-recruit analysis in mixed-stock Pacific salmon fisheries. *Canadian Journal of Fisheries and Aquatic Sciences*, 74(7), 1108-1122.

Sullivan, M.G. (2003). Active management of walleye fisheries in Alberta. *North American Journal of Fisheries Management*, 23(4), 1343-1358.

Thomson, C.J. (1991). Effects of the avidity bias on survey estimates of fishing effort and economic value. *American Fisheries Society Symposium*, 12, 356-366.

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413-1432.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.C. (2021). Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.
