# Grays Harbor Recreational Dungeness Crab Harvest Estimation

## Method Version 1.0: Pooled CPUE Model

**Author:** Matthew George, Ph.D.
**Contact:** matthew.george@dfw.wa.gov
**Agency:** Washington Department of Fish and Wildlife (WDFW)
**Status:** Published, operational. The method of record for estimating recreational Dungeness crab harvest at Westport / Grays Harbor.
**Method version:** 1.0 (corresponds to pipeline code version v7.4; reference calibration season 2024-25).

------------------------------------------------------------------------

### How to read this document

This is the single authoritative reference for the pooled-CPUE harvest estimation pipeline (`01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd`), the Stan model it fits (`02_stan_models/crab_bss_pooled.stan`), the diagnostics it runs, and the inputs and outputs it uses. It is written for two audiences and is split into three parts:

- **Part I (Sections 1-6): For everyone.** Plain-language description of what the method does, the fishery, the data, how the estimate is built, and where it is valid. No statistics background required.
- **Part II (Sections 7-13): Running it next season.** The operational guide: prerequisites, step-by-step run, how to judge whether a season's estimate is trustworthy, the output catalog, the diagnostics, reproducibility, and the conditions under which the method stops applying.
- **Part III (Sections 14-20): Technical reference.** The full model specification, design rationale, limitations, glossary, and references.

The development history (how the model reached v1.0, the full change log from v3 through v7.4, and the convergence-debugging narrative) has been moved out of this document to keep it focused on the published method. It lives in `BSS-GH-pooled-CPUE-model-development-history.md`. Section 19 gives a one-screen summary and points there.

A note on naming: this is the "pooled" model because it uses a single catch-rate (CPUE) process shared across gear types. A separate model, documented in `BSS-GH-gear-type-CPUE-model-documentation.md`, instead estimates a CPUE process per gear type. The pooled model is the published v1 because it is the more robust of the two and answers the primary management question (total harvest with defensible uncertainty); the gear-resolved model is the alternative when modeled gear-type catch with full uncertainty is required.

================================================================

# PART I: FOR EVERYONE

------------------------------------------------------------------------

## 1. What this method produces

This framework estimates the total recreational Dungeness crab (*Metacarcinus magister*) harvest at Westport and the greater Grays Harbor area for a season. It combines four kinds of field observations (gear counts at the docks, trailer counts at the boat launch, dockside crabber interviews, and ingress/egress surveys) with a statistical model that fills in the days when no sampling occurred.

Each run produces:

- A total Dungeness crab harvest estimate for the port, with a 95% credible interval (a range that has a 95% probability of containing the true harvest, given the data and model).
- Monthly harvest trends showing when crabbing pressure peaks and how it changes through the season.
- Breakdowns by crabbing mode (shore, private boat, commercial/charter) and an approximate breakdown by gear type (pot, ring net, trap, snare).
- A weekend catch-rate effect (whether weekend catch rates differ from weekday rates).
- Daily estimates of "effective day length" at the docks when ingress/egress data are available.

**How confident are we?** The framework runs two independent estimation methods, a simple average-based approach (the Point Estimator, PE) and a Bayesian time-series model (the Bayesian State-Space model, BSS), then compares them. When the two agree and the BSS passes its convergence checks, confidence is high. The output includes formal diagnostics and a side-by-side comparison so a reviewer can judge reliability. Section 9 explains how to read those checks.

------------------------------------------------------------------------

## 2. The fishery and study area

The recreational Dungeness crab fishery at Westport is one of the highest-volume recreational crabbing operations on the Washington coast. Crabbers use four main gear types: crab pots (highest catch rate), ring nets, foldable/star traps, and snares. WDFW rules prohibit pots from late September through November, which creates a structural break in both effort and catch rates and is the reason the season is split into two sub-seasons (Section 5).

Commercial Dungeness crab vessels also crab recreationally before the commercial season opens, under the same daily limits as private boats. Their harvest is tracked separately through a vessel tally at the marina.

Westport sits on the south side of the Grays Harbor estuary. Recreational crabbing occurs from public docks (Floats 17-21), a jetty, beaches, a public boat launch, and the commercial marina. The "shore" component pools dock, jetty, and beach crabbing.

------------------------------------------------------------------------

## 3. The four data streams

| Stream | What is collected | What it tells the model |
|---|---|---|
| **Effort counts** | Instantaneous point-in-time counts of crab gear at the docks and boat trailers at the launch, by field surveyors | The primary indicator of how much crabbing activity is happening |
| **Crabber interviews** | Dockside trip-level records: group size, gear deployed and type, hours fished, crab kept, trip status | Catch rate (CPUE) and the mix of gear in use |
| **Commercial/charter tally** | Daily count of commercial and charter vessels at the marina during the recreational pre-season | The commercial/charter component of harvest, via expansion |
| **Ingress/egress (I/E) surveys** | All-day surveys recording crabber arrivals and departures every 15 minutes | A direct measurement of crabber-hours that calibrates the gear-count pathway |

The four input files that carry these streams are listed in Section 7; their exact schema and known quirks are documented in `04_input_files/README.md`.

------------------------------------------------------------------------

## 4. How the estimate is built

The core problem is that field crews cannot sample every day. Both methods solve the same problem (estimate harvest on unsampled days), but differently.

**The Point Estimator (PE): a simple average.** For each stat-week by day-type group, the PE averages the daily harvest on sampled days and multiplies by the number of days in that group (Pollock et al. 1994; Hahn et al. 2000). It is transparent and assumption-light, but it cannot fill a group that had zero samples and it produces no uncertainty bounds.

**The Bayesian State-Space model (BSS): a time-series curve.** The BSS fits a smooth curve through the daily effort and catch-rate data using a statistical time-series process, then uses that curve to estimate every day in the season, including unsampled days (Conn 2002; Staton et al. 2017). It accounts for the fact that adjacent days are correlated, fills gaps with honest uncertainty that grows the further a day is from the nearest observation, and produces credible intervals. It is more complex, takes roughly half an hour to a few hours per fit, and must be checked for convergence.

**Combining them.** For each population component, the framework checks the BSS fit against formal convergence criteria (Section 9). If the fit passes, its estimate is used; if not, the PE estimate is used as a fallback. The two are reported side by side so a reviewer can see where they agree and where they differ.

The headline harvest number uses the BSS posterior expected catch (the model's best estimate of the average catch) rather than a single simulated draw, following the standard distinction between estimation and prediction in hierarchical models (Gelman et al. 2013, Ch. 7).

------------------------------------------------------------------------

## 5. The three population components and the two sub-seasons

The harvest is built from three components, estimated separately and summed:

1. **Shore crabbers** (dock + jetty + beach). Effort is indicated by gear counts at the docks.
2. **Private boat crabbers.** Effort is indicated by trailer counts at the boat launch.
3. **Commercial/charter vessels** crabbing recreationally pre-season. Estimated by expanding the marina vessel tally.

The season is split into two **sub-seasons**, defined by the pot closure, and each is estimated independently:

- **Ring-net only** (Sep 16 to Nov 30): pots prohibited.
- **All-gear** (Dec 1 to Sep 15): pots allowed.

------------------------------------------------------------------------

## 6. Where this method is valid

Method v1.0 is calibrated to the Westport / Grays Harbor fishery as sampled in the 2024-25 season. It is designed to be re-run in future seasons **provided the fishing location, the input data streams, and the sampling design remain the same.** Section 13 sets out, in detail, which assumptions are baked in and the specific conditions under which the method must be re-derived rather than re-run. In short: a different port, a change in how effort counts are taken (for example, reverting from randomized counts to a single peak-time count), or a structural change in who participates would each require revisiting the method, not just feeding it new data.

================================================================

# PART II: RUNNING IT NEXT SEASON

------------------------------------------------------------------------

## 7. Prerequisites and repository layout

**Software.** R 4.2 or later, with rstan 2.32 or later and a working C++ toolchain (rstan compiles the model), plus the packages tidyverse, lubridate, suncalc, gt, patchwork, here, and readxl.

**Repository layout.** The pipeline relies on the numbered stage folders and on `here::here()`, which anchors all file paths to the repository root. You do not edit paths to run it; you place files in the right folders:

| Folder | What it holds | Your job |
|---|---|---|
| `01_BSS_models/` | The driver `BSS-GH-pooled-CPUE-model.Rmd` | This is the file you run |
| `02_stan_models/` | `crab_bss_pooled.stan` | Leave in place; the driver compiles it |
| `03_R_functions/` | Helper functions (auto-sourced) | Leave in place |
| `04_input_files/` | The four input files | Replace these with the new season's data, same names and schema |
| `05_output/` | Dated run folders | The run writes here; nothing to place |

**The four input files** (see `04_input_files/README.md` for schema and quirks):

- `effort_combined.csv` (effort counts)
- `interview_combined.csv` (crabber interviews)
- `wes_commercial_tally.csv` (commercial/charter tally)
- `ingress_egress.xlsx` (I/E surveys; named through the `ie_data_file` / `ie_sheet` parameters)

------------------------------------------------------------------------

## 8. Step-by-step: running a season

1. **Place the new season's data** in `04_input_files/`, keeping the four filenames and their column schemas unchanged. Honor the schema quirks in the input-folder README (the interview gear column maps from column N; re-export the effort CSV with full quoting; dates are M/D/YYYY; the "Commerical" boat-type spelling is matched by regex).
2. **Open `01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd`** and set the run parameters in the `params` chunk. The ones you normally touch:
   - `est_date_start`, `est_date_end`: the season window. The driver fits each sub-season inside this window.
   - `ie_data_file`, `ie_sheet`: the I/E workbook and sheet (defaults `ingress_egress.xlsx`, `data`).
   - `bss_seed`: the RNG seed (default 20260619). Leave fixed for reproducibility; change only if a pathological seed is ever suspected.
   - Sampler and gate controls (Sections 9 and 14) have sensible defaults and rarely need changing for the standard fishery.
3. **Run / knit the driver.** Each population by sub-season is fit independently. Expect a total runtime of roughly 3 to 6 hours on a 4-core machine, depending on AR resolution and sub-season length.
4. **Check convergence** for each fit using `convergence_report.csv` and the rules in Section 9. A fit that fails falls back to PE automatically.
5. **Read the outputs** from `05_output/YYYYMMDD/pooled-CPUE/` (Section 10), starting with `port_total_Dungeness_Kept.csv` and `pe_vs_bss_comparison.csv`.

------------------------------------------------------------------------

## 9. Judging whether a season's estimate is trustworthy

This is the most important section for an operator. For each BSS fit, the framework monitors four diagnostics, reported per fit in `convergence_report.csv`: rank-normalized split-R-hat and bulk effective sample size (n_eff) for the seasonal totals `C_expected_sum` and `E_sum` (Vehtari et al. 2021), the number of divergent transitions, and the percentage of iterations that saturate the sampler's tree depth.

A fit **passes**, and its BSS estimate is used for that component, when all of the following hold; otherwise the PE estimate is used:

- **R-hat < 1.01** for both totals. R-hat near 1.00 means the independent sampler chains agree.
- **n_eff > 400** for both totals. This is the effective number of independent posterior samples.
- **Divergent fraction below 0.15** (the hard backstop). Above this rate, the sampler's geometry is untrustworthy and the fit is rejected regardless of anything else.
- **Divergences do not move the answer.** The shift the divergent draws induce in each total, measured in units of that total's posterior standard deviation (`|median(all) - median(bulk)| / sd(all)`), is below 0.10 SD.

Why divergences are in the gate: a sampler that cannot accurately integrate its trajectory is not faithfully exploring the target distribution, and can bias the posterior even when R-hat and n_eff look fine (Betancourt 2017).

Why the impact criterion is measured in standard deviations and not as a percentage of the estimate: a percentage-of-level threshold penalizes a component for having a wide posterior (being weakly identified) rather than for being biased. The SD-normalized criterion asks the question the gate exists to answer (do the divergences move the answer relative to how well the answer is pinned down) and is invariant to how wide the posterior is. This matters directly for the private boat, whose posterior is genuinely wide; see Section 16.

**Reading the comparison.** `pe_vs_bss_comparison.csv` shows PE and BSS effort and catch by component with the selected method. Large PE-vs-BSS gaps are not automatically errors; they can reflect a real disagreement between the design-based expansion and the model's reconciliation against interview data (the private boat is the standing example, Section 16). The convergence report and the comparison now always agree on which method was used, because the gate decision is computed once per fit and consumed by every downstream summary (this was a v7.0 fix; see Section 19).

------------------------------------------------------------------------

## 10. Output catalog

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
| `sensitivity_incomplete_trips.csv` | PE catch with the incomplete-trip filter off vs on, and the % change (v7.5; the harvest impact of `filter_incomplete_trips`) |

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

## 11. Diagnostics: what each one answers

The diagnostics are additive (each is wrapped so a failure cannot break a run) and are written every run. The three that an operator should be able to read:

**Posterior predictive checks (PPC).** `ppc_calibration_<label>.csv` and `ppc_pit_<label>.png` ask whether the model's predictions are calibrated against the actual effort counts and interview catches. A well-calibrated model has PIT values spread uniformly; a central hump means the predictive is too wide (over-dispersed), and 50% coverage above the nominal 0.50 says the same. In the reference run the effort predictive is somewhat over-dispersed (gear/trailer 50% coverage around 0.63 to 0.75), which is what the next diagnostic dissects.

**Effort over-dispersion decomposition.** `effort_overdispersion_decomp_<label>.csv` splits each effort observation's predictive variance into three additive parts via the law of total variance, so the lever behind any over-dispersion is identified before any prior or model change:

```
Var(Y) = E[mu]            (Poisson floor: irreducible, not a lever)
       + E[mu^2 / r_E]    (NB observation over-dispersion: controlled by the r_E / sigma_r_E prior)
       + Var(mu)          (latent process + parameter uncertainty: controlled by sigma_eps_E)
```

The `lever` column reports the verdict. The decision rule: if the NB-overdispersion share dominates, the lever is the `r_E` / `sigma_r_E` prior (the cheaper, exact change); if the latent share dominates, the lever is the AR innovation scale, which is a more delicate change (the boat tends to show a larger latent share). The analytic decomposition was checked against a brute-force Monte Carlo predictive variance and matches within Monte Carlo noise. Two standing cautions apply: any such correction is a prior/inference change that needs a guarded test run, and tightening the effort dispersion narrows the reported intervals (including the headline summer intervals), which is a change to reported uncertainty and needs explicit sign-off. The target is calibration (50% coverage near 0.50), not zero over-dispersion; some over-dispersion is real.

**PSIS-LOO.** `loo_summary_<label>.csv` reports out-of-sample predictive performance (expected log predictive density, `elpd_loo`) and the Pareto-k influence diagnostic per likelihood component (gear/trailer/catch). This is the basis for principled model comparison; it is what was used to evaluate, and reject, weather covariates (Section 17).

**CPUE effort-unit checks (v7.5).** `cpue_estimators_<label>.csv`, `cpue_saturation_<label>.csv`, and `cpue_linearity_<label>.csv` test the likelihood's core assumption that catch is proportional to the chosen effort denominator `h`. The estimator triad reports the model-implied CPUE (`C_expected_sum / E_sum`) against the ratio-of-sums and the mean-of-ratios; a model sitting near the mean-of-ratios is a warning that the negative-binomial dispersion is pulling `lambda_C` off the rate scale. The saturation exponent fits `catch_per_gear ~ (hours_per_gear)^beta` (boat only) and the linearity check fits `glm(catch ~ log(h))`; the likelihood assumes `beta = 1`, so a value well below 1 means the effort unit is not valid (for pots, catch is nearly flat in soak time). The run also asserts that effort `E` and the CPUE denominator `h` carry the same unit. These are diagnostic only; they are what would surface a boat or shore effort-unit defect before its total is trusted.

**Why the deployment is the effort unit (saturation).** Binned by soak time, crab per gear-HOUR falls about 43-fold across the range of soak durations, while crab per gear per trip rises only about 1.8-fold; a log-log fit gives catch per gear scaling as soak-hours to the power ~0.13. In plain terms, soak time barely matters, so a pot lift is a pot lift whether it soaked two hours or eight. That makes the deployment the unit on which catch-per-unit-effort is a stable rate: roughly 4 to 7 crab per pot lift, steady across soak times. A stable rate is exactly what the harvest method needs, because harvest is effort multiplied by that rate, and the multiplication is only unbiased if the rate does not drift with the effort denominator. On the 2026-07-10 shore comparison this is exactly what the linearity diagnostic shows: gear-deployments is the only shore unit whose `beta_h` covers 1 (1.05, 95% CI 0.94 to 1.15), while crabber-hours (0.57) and gear-hours (0.73) both fall well short, so v7.7 moves shore, like the boat before it (v7.6), onto the deployment scale.

**Incomplete-trip filter and its sensitivity (v7.5).** CPUE is computed from completed trips only (`filter_incomplete_trips`, default on): incomplete trips have soak-time gear that has not finished and read systematically low (about -20% for pots and traps). `sensitivity_incomplete_trips.csv` reports the PE catch with the filter off vs on so the harvest impact is explicit each run. Missing trip status is kept (a blank `completed_trip` may still be a complete trip).

**Config levers.** Two experiment toggles default to production behavior and are documented in the driver's `params`: `collapse_mu_hier` (default off) collapses the single-cell mu-hierarchy per population for the funnel investigation, and `ar_force` (default null) forces a population's AR resolution. Both leave the default posterior unchanged.

------------------------------------------------------------------------

## 12. Reproducibility

The Stan fits take a fixed RNG seed (`bss_seed`, default 20260619), passed to `rstan::stan()`. rstan seeds each chain from `bss_seed + chain_id`, so the chains still differ (R-hat remains meaningful) while run-to-run variation is removed. Package and Stan versions and the seed are written to `session_info.txt` with each output set. Change `bss_seed` only if a pathological seed is ever suspected. Expected runtime is 3 to 6 hours on a 4-core machine.

------------------------------------------------------------------------

## 13. Scope: when this method applies, and when it must be re-derived

Method v1.0 is built for one fishery under one sampling design. It can be re-run season after season as long as the following hold. Where one breaks, the method must be revisited, not merely re-fed.

**Assumptions that allow a straight re-run:**

- **Same location.** Westport / Grays Harbor access points (docks Floats 17-21, the jetty, beaches, the boat launch, the marina). The gear-per-crabber prior `R_G`, the gear-per-boat-group prior `R_G_boat`, and the effective-day-length regression are all calibrated to this site.
- **Same input streams, same schema.** The four input files in the same form (Section 7).
- **Same sampling design.** Instantaneous effort counts, dockside interviews, the commercial tally, and I/E surveys, collected as in 2024-25. The 2024-25 protocol of three randomized effort counts per day is the design the gear-hours expansion assumes (it measures mean daily effort).
- **Same sub-season structure.** Ring-net only Sep 16 to Nov 30; all-gear Dec 1 to Sep 15, tied to the pot closure.

**Conditions that require re-derivation, not just new data:**

- **A different port.** `R_G`, `R_G_boat`, and the L_effective regression would have to be re-estimated from that port's I/E and interview data; the access-point structure differs.
- **A change in the effort-count protocol.** Reverting to a single peak-time count per day measures a different quantity (peak, not mean daily effort) and would bias the effort level high. Mixing protocols across years is a genuine confound, addressable only with a protocol fixed effect and a peak-to-mean calibration (the multi-year question in Section 17).
- **A structural change in participation.** For example, a change in how commercial/charter vessels participate pre-season, or the opening of a new major access point (a jetty effort count, currently absent), would change what the components represent.
- **A season with large sampling gaps coinciding with anomalous weather.** The routine model interpolates gaps with its time-series process and deliberately excludes weather (Section 17). A season with extended unsampled stretches under unusual conditions is the one case where the shelved parsimonious weather-effort contingency (Section 17) should be considered, evaluated by leave-one-week-out block cross-validation.

================================================================

# PART III: TECHNICAL REFERENCE

------------------------------------------------------------------------

## 14. Model specification (`crab_bss_pooled.stan`)

### 14.1 Effort process

```
log(lambda_E[d]) = mu_E + omega_E[period(d)] + B1 * w[d] + B2 * holiday[d]
```

The temporal deviation `omega_E` evolves as an AR(1) process:

```
omega_E[p] = phi_E * omega_E[p-1] + sigma_eps_E * epsilon[p-1]
```

where `period(d)` maps day `d` to its AR period index. At daily resolution `period(d) = d` and the number of periods `P_n = D`; at weekly or monthly resolution `period(d)` maps to the week or month index and `P_n` is the number of weeks or months. Innovations `epsilon` are standard normal (non-centered parameterization for efficient HMC; Papaspiliopoulos et al. 2007). The AR(1) initial state is non-centered: `omega_E_0` is a raw standard normal scaled in the transformed-parameters block by the stationary standard deviation `sigma_eps_E / sqrt(1 - phi_E^2)`, so the process starts from its stationary distribution without a centered funnel (Harvey 1989; Betancourt and Girolami 2015).

**Adaptive temporal resolution.** The AR resolution is selected automatically per fit from effort-data density:

| Resolution | Condition | Rationale |
|---|---|---|
| Daily | >= 25% of days sampled AND >= 20 effort days | Dense data supports day-level smoothing with proper uncertainty scaling by distance from the nearest observation (Staton et al. 2017) |
| Weekly | >= 1.5 effort obs per week AND >= 3 weeks | Moderate data; weekly states smooth 3-5 day gaps without under-identifying the AR |
| Monthly | Fallback for sparse data | Few AR parameters; robust with limited observations |

This applies the finest resolution the data can identify and falls back gracefully when it cannot (Conn 2002; Sullivan 2003). The data-driven choice is additionally capped per population via `ar_max_resolution`: the boat fit is capped at weekly regardless of coverage, because the trailer-count series cannot identify a daily latent process even when its coverage exceeds the daily threshold (coverage counts how many days carry an observation, not how strongly each observation constrains the latent process). Shore is uncapped at daily, where it converges with n_eff above 2000. An `ar_force` parameter can override both the data-driven rule and the cap for a single population, used for the boat daily-vs-weekly resolution experiment; it defaults to `NULL` (production behavior).

### 14.2 CPUE process

```
log(lambda_C[d]) = mu_C + omega_C[period(d)] + B1_C * w[d]
```

`B1_C` allows weekend CPUE to differ from weekday CPUE, motivated by evidence that weekend/holiday crabber populations at tourist-accessible ports include more novice participants (Thomson 1991; Pollock et al. 1997). In the reference data `B1_C` is about -0.25 to -0.30 for shore crabbers (weekend crabbers catch roughly 21-26% fewer crab per crabber-hour than weekday regulars), consistent with the novice-dilution hypothesis. This is a single pooled CPUE process; gear-type catch is apportioned afterward from interview proportions (the gear-resolved model is the alternative that models per-gear CPUE).

### 14.3 Observation models

- Gear counts (shore): `Gear_I ~ NegBinomial2(lambda_E[d] * R_G, r_E)`
- Trailer counts (boats): `T_I ~ NegBinomial2(lambda_E[d] / R_G_boat, r_E)` (POOL-1; lambda_E is gear, lambda_E / R_G_boat is boat groups)
- Interview catch: `c ~ NegBinomial2(lambda_C[d] * h, r_C)`, where `h` = crabber-hours (shore) or gear-hours (boats)
- I/E crabber-hours: `IE_crabber_hours ~ Lognormal(log(lambda_E[d] * L[d]), sigma_IE)`

The negative binomial accommodates the overdispersion typical of recreational trip-level catch (Maunder and Punt 2004).

### 14.4 Effort overdispersion (marginalized)

Each effort count is negative binomial with shape `r_E`. This was originally written as a Gamma-Poisson mixture with an explicit per-observation latent multiplier `eps_E_H_obs ~ Gamma(r_E, r_E)`. Because the Gamma-Poisson mixture integrates exactly to the negative binomial (Hilbe 2011), the latent multipliers are marginalized analytically and the negative binomial is written directly. The change is inference-preserving (the marginal likelihood is identical), removes a high-dimensional centered latent block from the sampler, and makes the model block consistent with the `log_lik` block. The data field `n_effort_obs` is retained as an unused field to keep the R prep interface stable.

### 14.5 I/E integration and effective day length

On I/E survey days, observed crabber-hours enter as a direct lognormal observation of `lambda_E * L`, a second independent constraint on the latent effort state that bypasses `R_G` and day-length assumptions and calibrates the gear-count pathway against the I/E ground truth (Robson 1991; Pollock et al. 1994). When no I/E data are available (`IE_n = 0`), the I/E likelihood contributes nothing and the effort and catch posteriors are unchanged. The prior on the I/E scale `sigma_IE ~ exponential(5)` is applied unconditionally (not only inside the `IE_n > 0` branch), so that with no I/E data `sigma_IE` is still proper rather than an improper flat direction; because `sigma_IE` is decoupled from effort and catch, this leaves those posteriors unchanged.

When `estimate_L = 1` (shore fits), effective day length `L[d]` is a parameter with a non-centered lognormal prior:

```
L[d] = L_mu[d] * exp(L_sigma[d] * L_raw[d]),    L_raw ~ Normal(0, 1)
```

`L_mu` and `L_sigma` come from a regression of log effective day length on day-of-year (quadratic) and day type, fit from the I/E data:

```
log(L_effective) = b0 + b1 * yday + b2 * yday^2 + b3 * weekend + e
```

The quadratic captures the seasonal arc. Effective day length at the docks averages about 3.5 to 5.5 hours, substantially shorter than civil twilight (9 to 16 hours), because crabbers rotate through the dock rather than occupying it all day. On I/E days `L[d]` is further constrained by the I/E likelihood; on other days it is informed by the regression prior and its uncertainty propagates into effort and catch (Pollock et al. 1994; Hartill et al. 2012). For boats, the effort unit is gear-deployments as of v7.6 (the daily expansion is `L = tau_boat`, the deployment turnover, rather than a fixed 24-hour soak), because catch does not scale with soak time for pots; as of v7.7 shore uses the same gear-deployment unit. See Sections 11, 16, and 19.

### 14.6 Key parameters and priors

| Parameter | Description | Prior |
|---|---|---|
| B1 | Weekend effort multiplier (log) | Normal(0, 1) |
| B2 | Holiday effort multiplier (log) | Normal(0, 1) |
| B1_C | Weekend CPUE effect (log) | Normal(0, 1) |
| R_G | Gear per crabber | Lognormal(log(R_G_empirical), 0.3), data-driven |
| R_G_boat | Gear per boat group | Lognormal(log 4, 0.5) |
| phi_E, phi_C | AR(1) autocorrelation | Beta(2, 2) rescaled to [-1, 1] |
| r_E, r_C | Overdispersion | Half-Cauchy(0, 1) |
| sigma_IE | I/E measurement error (log) | Exponential(5) |
| L[d] | Effective day length (shore) | Lognormal from the I/E regression |

Prior rationale: `R_G` is centered on the empirical gear-per-crabber ratio in the relevant population by sub-season, eliminating prior-posterior conflict; `R_G_boat` Lognormal(log 4, 0.5) is centered on ~4 gear per boat group (POOL-1; replaces the old R_T Beta(5, 1), which was pinned at 1 by a degenerate bernoulli term); the Half-Cauchy(0, 1) variance priors are weakly informative (Gelman 2006).

### 14.7 Generated quantities

- `C_expected[d] = lambda_E[d] * L[d] * lambda_C[d]`: the posterior expected daily catch, E[C | data], the quantity used for harvest estimation.
- `C[d] = Poisson_rng(C_expected[d])`: a predictive draw including Poisson sampling noise, reported separately for prediction intervals.

For seasonal totals the Poisson noise largely averages out, so the two are similar; for daily or monthly breakdowns the difference can be material. Pointwise `log_lik` for the gear, trailer, and catch streams is also produced, enabling PSIS-LOO.

------------------------------------------------------------------------

## 15. Design decisions and their rationale

- **Two sub-seasons** are estimated independently because the pot closure creates a structural break in both effort and catch rates; pooling across it would blur two different regimes.
- **Gear-deployments for both components (v7.6 boat, v7.7 shore).** Effort is denominated in gear-deployments (pot lifts), not soak-hours or crabber-hours. Earlier versions used gear-hours for the boat on the reasoning that soak time is the fishing-time measure for crab gear, but the saturation diagnostic overturned that: catch per gear-hour falls roughly 43-fold across soak durations while catch per pot lift is nearly flat (catch per gear scales as about soak-hours^0.13 for the boat and ^0.22 for shore), so any time-denominated unit violates the likelihood's proportionality assumption (`beta_h = 1`). The deployment is the unit on which CPUE is a stable rate (about 4 to 7 crab per pot lift), which is what keeps harvest = effort x CPUE unbiased. Crabber-hours for shore failed the same test (`beta_h = 0.57`), so v7.7 moves shore onto deployments too. See Sections 11, 16, and 19.
- **L_effective from I/E** rather than civil twilight, because the dock activity curve is peaked (crabbers rotate through), so civil twilight overstates the time gear is actively fished by roughly a factor of two.
- **A weekend CPUE effect (B1_C)** because weekend crabber composition differs (more novices), which the catch-rate process should be allowed to reflect.
- **Expected catch, not a predictive draw, as the headline** because harvest estimation wants E[C | data], the estimation quantity, not a single noisy prediction.

------------------------------------------------------------------------

## 16. Limitations and the private-boat caveat

**General limitations.**

- I/E coverage is limited; expanding toward roughly 40 days per season would improve the L_effective regression (Pollock et al. 1997 recommend at least 3 I/E days per month by day-type stratum).
- The L_effective regression is a quadratic in day-of-year; with more data a GAM could capture non-monotonic patterns.
- Gear-type breakdowns from this pooled model are approximate (proportional allocation, not modeled). The gear-resolved model is the alternative when modeled gear catch with uncertainty is needed.
- There are no jetty effort counts; beach crabbing is unmeasured within the pooled shore count.
- The weekend CPUE effect is constant across the season; a time-varying effect may be warranted if tourist composition shifts seasonally.
- The adaptive AR selection is rule-based; a formal LOO/WAIC model comparison could provide principled resolution selection (Vehtari et al. 2017).

**The private boat.** The private-boat all-gear fit rests on a thin, weakly informative trailer-count series, and is the component most prone to wide posteriors and to PE-vs-BSS disagreement. A sequence of fixes was applied to make the fit converge: dedicated sampler tuning, a per-population AR cap to weekly, non-centering of the AR initial state, an unconditional `sigma_IE` prior (the boat has no I/E data, so this removed an improper flat direction that had inflated divergences), and finally the scale-aware convergence gate. Under that gate the boat all-gear component is reported on its BSS posterior, with a wide 95% interval and a catch CV around 27%, rather than the narrower PE point. The wide interval is not a defect to be hidden by substituting the PE point; it is the effort and CPUE data honestly reporting their own uncertainty for a component identified by only a handful of interviews per month. The PE-vs-BSS gap for the boat is a real disagreement between the design-based trailer expansion (which assigns every trailer-day a full gear-per-group times 24 hours) and the model's reconciliation against interview-reported soak times (which include zero-hour deployment-day interviews and real soak durations), not an arithmetic error; the BSS is the better-reconciled of the two. The durable fix is a more informative effort series (for example, access-point or camera exit counts), not further parameter surgery on the shared model.

**v7.6 update (POOL-1 + POOL-3).** The effort-UNIT half of this caveat has now been fixed. The boat is moved onto the gear-deployment scale (Section 19): the trailer expansion becomes `T_I ~ NB2(lambda_E / R_G_boat)` with gear-per-group learned as `R_G_boat` from interviews, and the CPUE denominator becomes `number_of_gear` (deployments) with `L = tau_boat`, replacing the gear-per-group-times-24-hours expansion described above. This is justified because catch is not linear in soak time for pots (the saturation diagnostic gives catch per gear ~ soak-hours^0.25), so any time-denominated effort unit is invalid; the gear-resolved model on this scale reports a boat total ~25% below the old gear-hours figure. The DATA-THINNESS half of the caveat remains: the boat still rests on a sparse trailer-count series, and `tau_boat` leans on its prior until boat ingress/egress counts accumulate, so the wide posterior and the "more informative effort series" recommendation still stand. The v7.6 boat total must be confirmed by a run.

------------------------------------------------------------------------

## 17. Weather and tide covariates: evaluated and excluded

A weather-and-tide covariate module (`06_diagnostics/`) was built and run on the 2024-25 season to test whether weather improves the estimate. The conclusion, documented in full in `WEATHER_COVARIATE_ANALYSIS.md`, is that **covariates are excluded for all three components** under a pre-committed 4.0-SE PSIS-LOO improvement margin. No component cleared the margin; shore all-gear was a tie, and shore ring-net and the boat were meaningfully worse out-of-sample with covariates.

The instructive findings:

- **Weather drives effort, not CPUE.** Boat effort is suppressed by wave height and rain and raised by temperature and tide range; shore effort is suppressed by rain. CPUE is essentially weather-flat. These are descriptive effort-dynamics results worth reporting (and they bear on sampling design, since effort is predictably low on rough or rainy days), but they do not improve harvest prediction.
- **False precision.** The covariate models produced narrower credible intervals while predicting held-out data worse. Narrower-but-worse is the signature of overfitting, which is exactly why selection here uses LOO and not interval width.
- **Why significant effects do not help.** On a sampled day, the effort count already encodes the weather effect (a rough day shows a low count whether or not weather is in the model). Weather could only add value by improving interpolation on unsampled days, and the AR(1) process already does that interpolation; on routine data, weather does not beat it.
- **The one exception, kept on the shelf.** If a season has a stretch with no effort counts (a true sampling gap), the AR's redundancy with the counts disappears, and a parsimonious weather-effort model (a few strong drivers, not the full screened set) could correct the AR's naive interpolation, most valuably when the unsampled period's weather is anomalous. The observation-level LOO used here cannot test that scenario (it never holds out a full week); the correct test is leave-one-week-out block cross-validation. Weather is therefore not part of the routine production model, but is worth keeping as a contingency for grounding effort across sampling gaps.

A known reconciliation item: the covariate module's absolute boat effort is well below the main pipeline at identical CPUE, because the module predates the current gear-hours expansion. This does not affect the covariate decision (the LOO is computed on the latent fit and is invariant to the post-hoc expansion) or the production estimate (the main pipeline is authoritative), but the module's boat magnitudes should not be read as harvest until reconciled.

------------------------------------------------------------------------

## 18. Glossary

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
| C_expected | Expected daily catch (no Poisson noise); E[C | data] |
| C | Predictive daily catch draw (includes Poisson noise) |
| L_effective | Effective day length: I/E crabber-hours divided by peak crabbers present |
| L_mu, L_sigma | Regression-predicted median and uncertainty for L_effective |
| I/E | Ingress/Egress survey |
| sigma_IE | I/E measurement error on the log scale |
| n_eff | Bulk effective sample size |
| R-hat | Rank-normalized split potential scale reduction factor (convergence) |
| PSIS-LOO | Pareto-smoothed importance-sampling leave-one-out cross-validation |
| PIT | Probability integral transform (posterior predictive calibration) |

------------------------------------------------------------------------

## 19. Development history (summary)

Method v1.0 corresponds to pipeline code **v7.4**. The model began as an adaptation of the WDFW freshwater-creel state-space framework and was hardened over a sequence of versions in response to a 2026-03-31 model critique and an extended convergence-debugging effort focused on the private boat. The arc in one screen:

- **v3-v5:** shared state-space milestones; the gear-resolved track branched at v5.
- **v6.0:** post-critique modeling upgrades (adaptive AR resolution; L_effective as an estimated parameter; the B1_C weekend CPUE effect; direct I/E integration; data-driven R_G; sparse effort overdispersion; expected and predictive catch both reported).
- **v6.1-v6.6:** the convergence gate gained divergence awareness; boat sampler tuning, a per-population AR cap to weekly, and non-centering of the AR initial state addressed boat non-convergence.
- **v6.7-v6.8:** effort overdispersion marginalized to negative binomial (inference-preserving); an unconditional `sigma_IE` prior fixed an improper direction specific to the boat; per-fit model diagnostics, a fixed seed, and session capture were added.
- **v6.9-v6.9.1:** a single-cell scale collapse (B1.7) was attempted and reverted after it hung the shore all-gear fit (the standing lesson: the durable boat fix is a better effort series, not parameter surgery); PPC calibration was hardened; monthly catch by mode was added.
- **v7.0:** the scale-aware convergence gate (impact measured in posterior standard deviations, not as a percentage of level), which moved the boat onto its BSS posterior and made the gate control the selection rather than merely label it; a PE monthly effort-share fix; a PPC extraction fix.
- **v7.1-v7.4:** the effort over-dispersion decomposition diagnostic; an extended set of persisted per-fit outputs (the O-series); pointwise `log_lik` enabling PSIS-LOO on the pooled model; and the `ar_force` experiment toggle (a tight-pin attempt in v7.3 was reverted in v7.4 after it tipped the shore funnel into failure).
- **v7.5 (2026-07-10):** the pooled backlog fixes POOL-2/4/5/6. The R layer was de-duplicated onto the shared gate and AR selector (POOL-6, behavior-preserving); the CPUE effort-unit diagnostics were wired in (POOL-5); a `collapse_mu_hier` lever was added for the funnel investigation (POOL-4, default off); and the incomplete-trip filter was added (POOL-2, default on), which raises the shore estimate, so a re-run is needed to refresh these numbers. The boat-structure items POOL-1 and POOL-3 were held for a validated session because they move the publication boat number. See `development_notes/20260710-OUTSTANDING_ISSUES.md`.
- **v7.6 (2026-07-10):** POOL-1 + POOL-3. The private boat is moved onto the gear-deployment scale, matching the gear-resolved model: `R_T` (pinned at ~1) is replaced by `R_G_boat` with `T_I ~ NB2(lambda_E / R_G_boat)` and `Gear_A_boat ~ poisson(R_G_boat)`, and the boat CPUE denominator becomes `number_of_gear` with `L = tau_boat` instead of gear-hours with `L = 24`. This resolves the private-boat effort-unit caveat in Section 16 and moves the publication boat total (expected ~-25%, toward the gear-resolved boat), so v7.6 must be re-run before the totals are trusted. Shore is unchanged.
- **v7.7 (2026-07-11):** Shore moved onto the gear-deployment scale, so both components now share one effort unit (gear-deployments) with the gear-resolved model. The shore CPUE denominator becomes `number_of_gear` with `E_scale = R_G` and `L = tau_shore` (~1.7 turnover), replacing crabber-hours. This settles the shore half of the effort-unit question (backlog GR-16) using the 2026-07-10 shore LOO comparison (shore all_gear, n = 1649): gear-deployments is the only shore unit whose linearity coefficient covers 1 (`beta_h = 1.05`, 95% CI 0.94 to 1.15, flag off), against crabber-hours (0.57) and gear-hours (0.73), and the only one with no estimator-triad drift (ratio-of-sums 0.87 ~= mean-of-ratios 0.85 ~= model-implied 0.85 crab per deployment). gear-hours had a marginally better catch-stream `elpd_loo` (-3131 vs -3190 for deployments), but that predictive edge comes from the CPUE process absorbing the sub-linearity, which is what biases the season expansion; the choice therefore prioritizes harvest-unbiasedness over marginal predictive fit. There is no Stan change (the v7.6 `effort_scale_gear` / `E_scale` machinery already supports shore), but the shore publication number moves, so v7.7 must be re-run before the totals are trusted.

The full change log, with the per-version rationale, the divergence-diagnostic narrative, and the detailed B1.5 / B1.6 working notes, is in **`BSS-GH-pooled-CPUE-model-development-history.md`**.

------------------------------------------------------------------------

## 20. References

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
