# BSS-GH Weather & Tide Covariate Module — Technical Documentation

**Author:** Matt George (WDFW), with implementation support
**File:** `BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd` + `crab_bss_pooled_weather_adjusted.stan`
**Version:** 0.2.0
**Date:** 2026-06-17
**Companion to:** `BSS-GH-pooled-CPUE-model.Rmd` and `BSS-GH-pooled-CPUE-model-documentation.md`

---

## 1. Overview

This module extends the Bayesian State-Space (BSS) pooled CPUE model for the Grays Harbor recreational Dungeness crab fishery by testing and optionally incorporating environmental covariates, specifically **tide features** and **weather/sea-state features**, into the daily log-linear predictors for effort (`lambda_E`) and CPUE (`lambda_C`).

### 1.1 Hypotheses tested

| ID | Hypothesis | Test vehicle |
|----|-----------|--------------|
| H1 | Boat departure times cluster on rising/flood tide, after accounting for daylight availability. | Binomial test of interview-level departure times against tide-record-derived expected flood fraction during daylight. |
| H2 | Daily boat effort (trailer counts) varies with tide features after accounting for weekend, holiday, and seasonality. | Daily GAM screen with quasipoisson link; LOO-PSIS confirmation in BSS. |
| H3 | Daily shore CPUE varies with tide features (tide range, daytime high-tide timing). | Daily GAM screen on catch ~ tide + offset(log_hours); LOO-PSIS confirmation in BSS. |
| H4 | Daily boat and shore effort depend on sea state (wind speed, wave height, precipitation). | Daily GAM screen; LOO-PSIS confirmation. |

### 1.2 Design philosophy

Three principles drove the design:

1. **Reusability.** The module runs for any date range each season without code changes. It inherits `est_date_start` and `est_date_end` from the main model parameters.
2. **Resilience.** Every data source has ranked backup stations. If the primary station returns no data for part or all of the date range, the module fails over to the next closest station automatically. A station-availability log is written each run.
3. **Conservatism about inclusion.** Covariates are only added to the BSS if they survive three filters: (a) pass a loose GAM screen at p < 0.10, (b) produce a PSIS-LOO ELPD improvement greater than 4·SE over the baseline, and (c) posterior sign of `gamma` is consistent with the biological hypothesis (`>80%` of posterior mass on the expected side of zero).

### 1.3 Integration flow

```
  BSS-GH-pooled-CPUE-model.Rmd          (existing; baseline)
          |
          |  produces dwg, ie_data, days, summaries
          v
  BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd   (this module)
          |
          |--> fetch_tide_with_fallback()      [Westport -> Toke Point -> La Push]
          |--> fetch_buoy_with_fallback()      [46211 -> 46029 -> 46041 -> 46087]
          |--> fetch_asos_with_fallback()      [KHQM -> KAST -> KOLM -> GSOD]
          |
          |--> compute daily covariates
          |--> boat departure clustering test  [hypothesis H1]
          |--> Layer B GAM screen               [H2, H3, H4]
          |--> fit baseline BSS + augmented BSS (crab_bss_pooled_weather_adjusted.stan)
          |--> PSIS-LOO comparison (w/ k-fold fallback)
          |--> decision rule: include covariates per (population × sub-season)
          |
          v
  final_bss_estimates_with_covariate_selection.csv
```

---

## 2. Data sources

### 2.1 Tide: NOAA CO-OPS (Center for Operational Oceanographic Products and Services)

| Station ID | Name | Role | Coordinates |
|-----------|------|------|-------------|
| 9441102   | Westport, WA | Primary | 46.9043, -124.1051 |
| 9440910   | Toke Point, WA | Backup (outer coast, same tidal regime, ~25 mi SSW) | 46.7076, -123.9672 |
| 9442396   | La Push, WA | Deep backup (north, phase-shifted) | 47.9133, -124.6367 |

- **API endpoint (data):** `https://api.tidesandcurrents.noaa.gov/api/prod/datagetter`
- **API endpoint (metadata):** `https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{id}.json`
- **Product preference:** `water_level` (verified observations) first, `predictions` as gap-fill and final fallback.
- **Datum:** MLLW (Mean Lower Low Water).
- **Units:** metric (meters).
- **Interval:** hourly.
- **Timezone:** LST/LDT (local standard + daylight saving). This avoids the mismatch between UTC tide data and local-time interview timestamps.
- **Request chunking:** CO-OPS verified data is limited to 31 days per request. The module paginates transparently.

### 2.2 Offshore wind/wave: NOAA NDBC (National Data Buoy Center)

| Buoy ID | Name | Role | Coordinates |
|---------|------|------|-------------|
| 46211   | Grays Harbor | Primary | 46.857, -124.244 |
| 46029   | Columbia River Bar | Backup #1 | 46.144, -124.510 |
| 46041   | Cape Elizabeth, WA | Backup #2 | 47.353, -124.731 |
| 46087   | Neah Bay | Deep backup | 48.494, -124.728 |

- **Historical endpoint:** `https://www.ndbc.noaa.gov/view_text_file.php?filename={id}h{YYYY}.txt.gz&dir=data/historical/stdmet/`
- **Realtime endpoint (current year):** `https://www.ndbc.noaa.gov/data/realtime2/{ID}.txt`
- **Fields used:** `WDIR` (wind direction), `WSPD` (wind speed m/s), `GST` (gust m/s), `WVHT` (significant wave height m), `DPD`, `APD`, `MWD`, `ATMP`, `WTMP`, `PRES`.
- **Known data gaps:** Buoy 46211 has experienced multi-month outages in prior years. The fallback chain is the primary reliability mechanism.

### 2.3 Land weather: ASOS via Iowa State IEM, GSOD fallback

| Station | Role | Source |
|---------|------|--------|
| KHQM (Hoquiam / Bowerman Field) | Primary ASOS | Iowa State IEM archive |
| KAST (Astoria, OR) | Backup ASOS #1 | Iowa State IEM |
| KOLM (Olympia, WA) | Backup ASOS #2 | Iowa State IEM |
| 727923-24227 (KHQM), 727910-24227 (KAST) | GSOD daily fallback | NCEI v1 API |

- **ASOS endpoint:** `https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py`
  - Variables requested: `tmpc, dwpc, sknt, gust, drct, p01m, vsby, skyc1`
  - Timezone: `Etc/UTC`, converted to `America/Los_Angeles` post-fetch.
  - Output format: `onlycomma`, with 'M' as missing sentinel.
- **GSOD endpoint:** `https://www.ncei.noaa.gov/access/services/data/v1`
  - Dataset: `global-summary-of-the-day`
  - Variables: `TEMP`, `WDSP` (wind speed knots), `MXSPD`, `GUST`, `PRCP` (precip inches).
  - Daily only; used only if all ASOS stations fail.

### 2.4 Caching

All fetched data is cached as RDS files under `cache/weather_tide/`, keyed by date range:
- `tide_{start}_{end}.rds`
- `ndbc_{start}_{end}.rds`
- `asos_{start}_{end}.rds`

Re-running the module uses cached files. Set `params$force_refetch <- TRUE` to override and pull fresh data, which you should do whenever extending the date range or if a primary station's verified data is updated post-hoc.

---

## 3. Key definitions

### 3.1 Tide-level features (daily summaries from hourly observations)

| Feature | Formula | Rationale |
|---------|---------|-----------|
| `tide_max_m` | `max(v)` over 24h | Maximum tide height; indicates diurnal spring-tide magnitude. |
| `tide_min_m` | `min(v)` over 24h | |
| `tide_range_m` | `tide_max_m - tide_min_m` | Spring/neap proxy. Large range = spring tide = strong currents. |
| `daytime_tide_max_m` | `max(v)` during daylight only | Whether a spring high falls within fishable hours matters more than nighttime highs. |
| `daytime_tide_min_m` | `min(v)` during daylight | |
| `daytime_tide_range_m` | `daytime_max - daytime_min` | Daytime-specific spring/neap proxy. |

### 3.2 Tide-phase features (daily summaries)

| Feature | Formula | Rationale |
|---------|---------|-----------|
| `daytime_high_hour` | Hour of max tide during daylight | Morning vs evening high affects when crabbers go out. |
| `n_high_tides_daytime` | Count of local maxima during daylight | 1 vs 2 daytime highs per day (reflects lunar cycle). |
| `flood_hours_daytime` | Hours during daylight with rising tide | Proxy for "good crabbing window" under flood-tide hypothesis. |
| `tidal_energy` | `sqrt(sum(dv^2))` | Proxy for current energy; large energy = strong spring currents = possibly worse crabbing. |

### 3.3 Interview-derived daily features (boat population only)

| Feature | Formula |
|---------|---------|
| `prop_depart_on_flood` | Fraction of interviewed boat trips with rising tide at `fishing_start_time`. |
| `prop_depart_near_high` | Fraction with departure within ±2h of nearest daytime high tide. |
| `mean_depart_tide_height_m` | Mean tide height at interview-level departure times for that day. |

Note: these use interview-level timing information (sub-daily) aggregated to daily covariates for BSS input. The sub-daily information enters only via aggregation.

### 3.4 Weather features (daily summaries)

| Feature | Source | Formula |
|---------|--------|---------|
| `wind_ms_buoy_mean` | NDBC | Daily mean wind speed from offshore buoy. |
| `gust_ms_buoy_max` | NDBC | Daily max gust. |
| `wave_h_m_mean`, `wave_h_m_max` | NDBC | Significant wave height daily stats. |
| `sst_c_mean` | NDBC | Sea surface temperature daily mean. |
| `pressure_hpa_mean` | NDBC | Barometric pressure daily mean. |
| `temp_c_min`, `temp_c_max` | ASOS | Land-station daily temperature range. |
| `precip_mm_total` | ASOS | Daily precipitation sum. |
| `wind_ms_asos_mean` | ASOS | Land-station daily mean wind (redundant with buoy but serves as fallback). |
| `visibility_mi_min` | ASOS | Daily minimum visibility; proxy for fog. |

### 3.5 BSS structural covariates

- `K_E`, `K_C`: integer, number of covariates on effort and CPUE linear predictors respectively.
- `X_E[D × K_E]`, `X_C[D × K_C]`: standardized (mean 0, SD 1) daily covariate matrices.
- `gamma_E[K_E]`, `gamma_C[K_C]`: posterior coefficients on the log scale.
- `prior_sd_gamma`: prior SD on each gamma; default 0.35.
- `log_lik[N_obs]`: pointwise log-likelihood in `generated quantities`, for PSIS-LOO.

### 3.6 Model comparison metrics

- `elpd_diff`: Difference in Expected Log Pointwise Predictive Density between baseline and covariate model.
- `se_diff`: Standard error of `elpd_diff` (computed by `loo::loo_compare`).
- `pareto_k`: Diagnostic for PSIS-LOO reliability per observation. `k > 0.7` is unreliable.
- **Inclusion rule:** `elpd_diff > 4 × se_diff` AND sign consistency AND `< 10%` of observations have `k > 0.7`. If Pareto-k threshold exceeded, fall back to 5-fold time-block CV.

### 3.7 Sub-seasons

- **Sub-season 1 (ring_net_only):** `est_date_start` to `pot_open_date - 1` (2024-09-16 to 2024-11-30 in the 2024-25 data).
- **Sub-season 2 (all_gear):** `pot_open_date` to `est_date_end` (2024-12-01 to 2025-09-15).

Covariate selection is run **per sub-season × population** because the biology differs between ring-net and pot crabbing.

---

## 4. Code structure

### 4.1 Section map of `BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd`

| Section | Purpose |
|---------|---------|
| 0       | Setup, parameters, station priority lists, cache setup. |
| 1       | NOAA CO-OPS tide fetch with metadata check, chunked pagination, verified + prediction gap fill, station fallback. |
| 2       | NDBC buoy historical + realtime fetch with station fallback. |
| 3       | ASOS (Iowa State IEM) fetch + GSOD fallback. |
| 4       | Daily covariate construction: tide features, weather features, interview-derived features. |
| 4.5     | Boat departure clustering test (binomial test of H1). |
| 5       | Layer B GAM screen for effort (per population) and CPUE (per population). |
| 6       | Covariate matrix assembly with VIF pruning. |
| 7       | BSS baseline + augmented fit wrapper (`prep_bss_crab_augmented`). |
| 8       | PSIS-LOO with k-fold time-block CV fallback. |
| 9       | Main execution loop over (population × sub-season × catch group). |
|         | Step 7: Inclusion decisions per (population × sub-season). |
|         | Step 8: Uncertainty comparison (CI widths with/without covariates). |
|         | Step 9: Final estimates using selected model per fit. |
|         | Step 10: Diagnostic plots. |

### 4.2 Key functions

| Function | Purpose |
|----------|---------|
| `coops_station_metadata(id)` | Fetch station metadata for availability check. |
| `coops_fetch_chunk(...)` | 31-day-max water-level or predictions fetch. |
| `coops_fetch_range(...)` | Paginated fetch for arbitrary date range. |
| `fetch_tide_with_fallback(stations, start, end, cache_dir, force_refetch)` | Full pipeline: try primary, check coverage, fall back, fill gaps with predictions. |
| `ndbc_fetch_historical(buoy, year)` | Parse NDBC fixed-width text files. |
| `ndbc_fetch_realtime(buoy)` | Pull current-year realtime feed. |
| `fetch_buoy_with_fallback(...)` | NDBC fallback pipeline. |
| `iem_asos_fetch(station, start, end)` | ASOS via Iowa State IEM. |
| `gsod_fetch(station, start, end)` | GSOD daily fallback. |
| `fetch_asos_with_fallback(...)` | ASOS + GSOD fallback pipeline. |
| `build_sun_calendar(dates, lat, lon, tz)` | Generate dawn/dusk/daylight-hours per date via `suncalc`. |
| `compute_tide_daily_features(tide_obs, sun_cal)` | Produce all daily tide features. |
| `compute_weather_daily_features(buoy_obs, asos, gsod, sun_cal)` | Produce all daily weather features. |
| `compute_interview_tide_features(interview_df, tide_obs)` | Attach tide height, rising/falling, near-high to each interview; aggregate to daily. |
| `departure_clustering_test(per_interview, tide_obs, sun_cal)` | Binomial test of H1. |
| `gam_screen_effort(effort_daily, covariates, population)` | mgcv::gam screen for effort covariates. |
| `gam_screen_cpue(interview_df, covariates, population)` | mgcv::gam screen for CPUE covariates. |
| `build_covariate_matrix(features, covariate_df, date_vec)` | Standardize + VIF-prune into a Stan-ready matrix. |
| `prep_bss_crab_augmented(...)` | Build `stan_data` list with `K_E`, `K_C`, `X_E`, `X_C`, `prior_sd_gamma`. |
| `bss_loo_compare(fit_base, fit_cov, ...)` | PSIS-LOO with k-fold fallback. |
| `kfold_time_block_cv(...)` | Refit both models over contiguous time blocks, compute held-out predictive log-density. |
| `mask_stan_data`, `extract_stan_data`, `compute_lpd_held_out` | Helpers for k-fold CV. |

### 4.3 Augmented Stan model

The augmented Stan model is `crab_bss_pooled_weather_adjusted.stan` (in `stan_models/`). It extends the pooled model with covariate design matrices `X_E`, `X_C` and their dimensions `K_E`, `K_C`, adding linear predictors on `mu_E` and `mu_C`. It collapses to the baseline pooled model when `K_E = K_C = 0`, so this one file serves both the baseline fit (built by zeroing `K_E`, `K_C`, `X_E`, `X_C`) and the covariate-augmented fit. There is no separate patch file: the augmented model is maintained as a complete standalone Stan file.

---

## 5. Libraries used

| Package | Purpose | Version tested |
|---------|---------|----------------|
| `tidyverse` | Data manipulation, plotting | ≥ 2.0 |
| `lubridate` | Date/time handling, timezone conversion | ≥ 1.9 |
| `suncalc` | Sunrise/sunset calculation per date and location | ≥ 0.5 |
| `mgcv` | GAM fits for Layer B screen | ≥ 1.9 (base R) |
| `loo` | PSIS-LOO and k-fold CV | ≥ 2.6 |
| `rstan` | Stan model fitting | ≥ 2.26 |
| `here` | Project-relative paths | ≥ 1.0 |
| `readxl` | Legacy Excel reads if interview data is sourced from xlsx | ≥ 1.4 |
| `httr` | HTTP calls to CO-OPS, NDBC, IEM, NCEI | ≥ 1.4 |
| `jsonlite` | JSON parsing for CO-OPS and NCEI responses | ≥ 1.8 |
| `geosphere` | Distance calculations for proximity-ranked station fallback | ≥ 1.5 |
| `patchwork` | Composing diagnostic plots | ≥ 1.1 |

All dependencies are CRAN-available and auto-installed by the setup block if missing.

---

## 6. Data import

### 6.1 Internal objects expected at runtime

The module expects these objects to exist in the R session (from a prior run of `BSS-GH-pooled-CPUE-model.Rmd`):

| Object | Type | Source |
|--------|------|--------|
| `dwg` | list | Main RMD `fetch_crab_data()`; must contain `$interview`, `$shore_effort`, `$boat_effort`. |
| `ie_data` | tibble or `NULL` | Main RMD `fetch_ie_data()`; I/E anchor observations. |
| `prep_days_crab()` | function | Main RMD; builds per-day calendar with day/week/month indices. |
| `prep_population_summary()` | function | Main RMD; population-specific summary for BSS prep. |
| `L_eff_model` | optional | Main RMD; day-length regression, used for shore sunrise-to-sunset calc. |

If any of these are missing, the module halts with an explicit instruction to run the main RMD first.

### 6.2 External data fetched at runtime

| Source | Data | Typical volume for a full year |
|--------|------|-------------------------------|
| CO-OPS | Hourly water level (verified + predictions) | ~8,760 rows |
| NDBC   | Hourly offshore met (wind, wave, SST, pressure) | ~8,000 rows (some gaps common) |
| ASOS   | Sub-hourly land met (observations typically every 20 minutes) | ~25,000-30,000 rows |
| GSOD   | Daily summaries | ~365 rows |

Total external data volume per season: roughly 50 MB uncompressed, 5-10 MB cached (RDS).

### 6.3 Failure modes and behavior

| Failure | Behavior |
|---------|----------|
| CO-OPS primary station unavailable | Automatic fallback to Toke Point, then La Push. Station log records which was used. |
| CO-OPS all verified unavailable | Fall through to predictions from primary station. |
| CO-OPS all sources fail entirely | Returns empty tide observations; covariate construction returns NA features; BSS runs baseline only. |
| NDBC primary buoy unavailable | Fall through chain 46211 → 46029 → 46041 → 46087. |
| NDBC all unavailable | Weather covariates use ASOS-derived features only (no wave data). |
| ASOS all unavailable | Fall through to GSOD daily summaries (no sub-daily resolution). |
| ASOS + GSOD both fail | Weather covariates return NA; BSS runs baseline only. |

All failure modes produce a written station_log CSV and continue execution. No partial failure halts the pipeline.

---

## 7. Outputs

All outputs are written to `output/{run_date}/pooled-CPUE-covariates/`.

### 7.1 Station logs and metadata

| File | Contents |
|------|----------|
| `tide_station_log.csv` | Which CO-OPS station was used, status, n_obs. |
| `ndbc_station_log.csv` | Which NDBC buoy was used, status, n_obs. |
| `asos_station_log.csv` | Which ASOS/GSOD station was used, status, n_obs. |

### 7.2 Covariates

| File | Contents |
|------|----------|
| `daily_covariates.csv` | One row per date in the estimation window, all tide + weather + interview-derived features. |

### 7.3 Hypothesis tests

| File | Contents |
|------|----------|
| `boat_departure_clustering_test.csv` | Binomial test result for H1 (observed vs expected flood-departure fraction). |
| `boat_departure_clustering_test_hourly_robustness.csv` | H1 test rerun with departure times rounded to nearest hour. |

### 7.4 GAM screen results

| File | Contents |
|------|----------|
| `gam_shore_eff_smooths.csv` | Layer B GAM smooth terms, edf, p-values for shore effort. |
| `gam_boat_eff_smooths.csv` | Same for boat effort. |
| `gam_shore_cpue_smooths.csv` | Same for shore CPUE. |
| `gam_boat_cpue_smooths.csv` | Same for boat CPUE. |
| `layer_b_selected_features.csv` | Features passing the Layer B p-threshold for each (population × component). |

### 7.5 BSS comparison results

| File | Contents |
|------|----------|
| `loo_comparison_{label}.csv` | ELPD-LOO difference with SE per (population × sub-season × catch group). |
| `kfold_cv_{label}.csv` | Per-fold held-out log predictive density (only when k-fold fallback triggered). |
| `pareto_k_{label}.csv` | Pareto-k values per observation for both baseline and covariate fits. |

### 7.6 Decision and final estimates

| File | Contents |
|------|----------|
| `covariate_inclusion_decisions.csv` | Per-fit decision: include covariates or not, with elpd_diff, SE, reason, sign-check notes. |
| `ci_width_comparison.csv` | 95% CI widths for `C_expected_sum` and `E_sum`, baseline vs covariates, percentage change. |
| `final_bss_estimates_with_covariate_selection.csv` | Final season estimates per (population × sub-season × catch group), using selected model. |
| `covariate_effects_posterior.csv` | Posterior median, 95% CI for each `gamma_E` and `gamma_C`. |
| `covariate_effects_forest.png` | Forest plot of `gamma` posteriors across fits and features. |

### 7.7 File naming convention

`{label}` in output filenames follows: `{population}_{sub_season}_{catch_group}`
- Examples: `shore_ring_net_only_Dungeness_Kept`, `private_boat_all_gear_Dungeness_Kept`.

---

## 8. Future improvements

Items identified during the initial build that warrant future consideration:

### 8.1 Data quality concerns to revisit

1. **Interview `fishing_start_time` quality.** The hypothesis test for H1 depends on self-reported departure times, which are known to be rounded. The module includes a robustness check (coarsening to nearest hour). If the hourly-rounded and exact-time results diverge, it signals the test is being driven by apparent timing precision rather than true tide-phase effects.
2. **Sub-season heterogeneity in covariate effects.** The current design fits covariates independently per sub-season, but does not share information across seasons. A multilevel extension with a shared prior on `gamma` across sub-seasons could stabilize estimates in data-poor sub-seasons.
3. **Per-gear CPUE covariates.** Currently `gamma_C` is one coefficient per covariate across all gear types. For the all_gear sub-season, pot CPUE and ring-net CPUE may respond differently to tide. Extending `gamma_C` to be gear-specific would test this, at the cost of doubling parameters.

### 8.2 Statistical refinements

4. **Non-linear effects via splines in Stan.** The current implementation uses linear effects of standardized covariates on log(lambda). Non-linear effects (e.g., effort maximized at moderate wave heights but dropping at extremes) could be modeled with a small basis spline. This would match the flexibility of the mgcv::gam screen in the final BSS.
5. **Interaction with day-type.** The effect of wind on boat effort may differ between weekends and weekdays (weekend crabbers are less likely to cancel). Adding `gamma_E * w[d]` interactions would test this.
6. **Lagged weather effects.** Yesterday's wind may suppress today's effort if conditions are forecast to persist. A 1-day lag on key weather covariates would be a simple extension.
7. **Formal temporal block LFO-CV.** The current fallback CV uses 5 random-ordered time blocks. A more statistically appropriate approach for time series is Bürkner/Gabry/Vehtari (2020) leave-future-out (LFO) CV, which predicts forward in time. The compute cost is high but justified for a production-grade comparison.

### 8.3 Operational improvements

8. **Annual station availability pre-check.** Run a once-a-year job to verify all ranked stations are still operational before the season starts. Currently this is done at runtime, which means a failed primary station triggers the fallback every run.
9. **Tide station distance precomputation.** The `geosphere::distHaversine` ranking of backup stations by proximity is computed each run; could be cached per project.
10. **Integration with main RMD output format.** The `final_bss_estimates_with_covariate_selection.csv` currently produces a distinct format. A future iteration should produce a drop-in replacement for the main model's output files so the covariate-selected results flow downstream without format conversion.
11. **Diagnostic report knitting.** Auto-knit an HTML or PDF diagnostic report summarizing station availability, covariate selection, LOO results, and effect-size plots. Currently results are in CSVs, which are machine-readable but not human-glanceable.

### 8.4 Known limitations flagged in the code

- **Station backup ordering for tide.** Toke Point is in Willapa Bay, not directly on the Pacific; its tide phase is similar to Westport but amplitude differs. Using Toke Point as a substitute height may introduce bias; phase-based features (rising/falling, hour of high) are more robust across stations than absolute height.
- **Mixed-mode shore aggregation.** Per project decision, dock/jetty/beach are all coded as "shore". If shore effort or CPUE response to tide differs meaningfully by sub-mode, the aggregate signal may be muted. Revisit if jetty/beach interview counts grow.
- **GSOD temporal resolution.** When all ASOS stations fail and the module falls back to GSOD, only daily summary statistics are available. Covariates derived from sub-daily aggregation (e.g., `visibility_mi_min` during daylight) are NA.
- **Production pooled Stan model reconciliation.** The augmented model `crab_bss_pooled_weather_adjusted.stan` preserves the pooled model's likelihood and adds covariate blocks on `mu_E` and `mu_C`, collapsing to the baseline when `K_E = K_C = 0`. It must be kept in step with `crab_bss_pooled.stan`: the two share the same effort and CPUE process code, so any change to the pooled model's likelihood or parameter naming has to be mirrored in the augmented file. A `log_lik`-level cross-check of the two files (and confirmation of `p_I_shore` vs `p_TI` naming) is still outstanding (see Planned patches).

---

## 9. Version history

| Version | Date       | Notes |
|---------|------------|-------|
| 0.2.0   | 2026-06-17 | Made the module runnable against the actual crab data and brought its BSS engine to pooled-model v6.5 parity. **Schema fix:** `compute_interview_tide_features()` filtered on `fishing_start_time`, which the crab interview export does not have (it records `interview_time`, an Excel serial datetime, and `hours_fished`). It now derives a departure timestamp as `interview_time - hours_fished` for completed trips, controlled by `use_interview_departure_proxy` (default TRUE), and degrades gracefully to empty departure features when no start time and no usable interview time are available. The CPUE GAM screen does not use departure features and the effort screen drops them when absent, so the screens run either way. The proxy assumes the interview is taken at trip end; the rigorous alternative (interval-level I/E departure data) is noted as future work. **BSS parity:** the module keeps its own copy of the prep and fitting code, so the v6.1 divergence-aware gate, the v6.2 boat sampler tuning (`adapt_delta` 0.99 / `max_treedepth` 13 for boat all-gear, 14 / 0.95 for shore all-gear, 12 / 0.95 for ring-net), the v6.4 R-hat threshold of 1.01, and the v6.5 per-population AR resolution cap (boat capped at weekly) were all ported. Both the baseline and covariate fits now use the same per-fit tuning, and a new `bss_convergence_check()` records divergences, R-hat, n_eff, and treedepth per fit to `convergence_report.csv`; a component is flagged when either fit fails the gate, since its LOO comparison is then unreliable. The module does not fall back to PE (it is a comparison tool, not an estimator). |
| 0.1.2   | 2026-06-17 | Removed the dead `bss_model_file_baseline` parameter (defined but never used). Both the baseline and covariate fits load the augmented model `crab_bss_pooled_weather_adjusted.stan`; the baseline is built by zeroing the covariate blocks (`K_E = K_C = 0`) so that baseline and covariate `log_lik` are produced by the same model and stay comparable for PSIS-LOO. Added an inline comment to prevent re-introduction. (The redundant baseline refit itself, recomputing a baseline the main pooled run already produced, is tracked as a separate item.) |
| 0.1.1   | 2026-06-17 | Reference and file reconciliation. Corrected documentation and code-comment references from the never-committed patch file (`crab_bss_pooled_covariates_patch.md`) and the non-existent `crab_bss_pooled_covariates.stan` to the standalone augmented model `crab_bss_pooled_weather_adjusted.stan`, which the module already loads for both the baseline (K=0) and covariate fits. Corrected the documented module filename to `BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd` and removed the stale, misfiled duplicate `stan_models/BSS-GH-weather-tide-covariates.Rmd`. The module now runs from a clean clone using the root Rmd. The `log_lik`-level reconciliation and `p_I_shore` vs `p_TI` naming check originally planned for 0.1.1 are not part of this change and remain open (now tracked as 0.1.2). |
| 0.1.0   | 2026-04-22 | Initial build. Tide/weather fetch with fallback, daily covariate construction, Layer B GAM screen, BSS augmented fits, PSIS-LOO comparison with k-fold time-block CV fallback, decision rule, and final-estimate selection. Covers 2024-25 season data out of the box; parameterized for annual reuse. Companion Stan patch (`crab_bss_pooled_covariates_patch.md`) supplies minimum-surgical additions to `crab_bss_pooled.stan`. |

### 9.1 Planned patches

- **0.1.3 (planned):** Cross-check the augmented `crab_bss_pooled_weather_adjusted.stan` against the production `crab_bss_pooled.stan` at the `log_lik` level; reconcile observation likelihoods and confirm `p_I_shore` vs `p_TI` naming.
- **0.2.0 (planned):** Add sub-season-shared priors on `gamma` (multilevel extension per item 2 in future improvements).
- **0.3.0 (planned):** Non-linear covariate effects via small basis splines per item 4.

---

## 10. Acknowledgments and references

### 10.1 Methodology references

- Vehtari, A., Gelman, A., & Gabry, J. (2017). *Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC*. Statistics and Computing, 27, 1413-1432.
- Bürkner, P.-C., Gabry, J., & Vehtari, A. (2020). *Approximate leave-future-out cross-validation for Bayesian time series models*. Journal of Statistical Computation and Simulation, 90(14), 2499-2523.
- Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R (2nd ed.)*. CRC Press.

### 10.2 Data source attributions

- NOAA Center for Operational Oceanographic Products and Services (CO-OPS). Water level observations at Westport, WA (station 9441102). https://tidesandcurrents.noaa.gov/
- NOAA National Data Buoy Center (NDBC). Standard meteorological observations, station 46211 (Grays Harbor). https://www.ndbc.noaa.gov/
- Iowa Environmental Mesonet (IEM), Iowa State University. ASOS archive for KHQM. https://mesonet.agron.iastate.edu/
- NOAA National Centers for Environmental Information (NCEI). Global Summary of the Day. https://www.ncei.noaa.gov/

### 10.3 R package attributions

See Section 5 for full package list. Key packages: `rstan` (Stan Development Team), `loo` (Vehtari, Gabry, et al.), `mgcv` (Wood).
