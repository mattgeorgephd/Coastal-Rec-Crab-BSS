# Recreational Crab Creel Estimation, Grays Harbor / Westport

**Agency:** Washington Department of Fish and Wildlife (WDFW)
**Lead:** Matt George
**Status:** 2024-25 season, first full implementation. Pooled and gear-resolved are the production models; the weather-tide covariate work is an experimental module.

---

## What This Does

Estimates total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area by combining two statistical approaches:

- **Point Estimator (PE):** Stratified expansion. Averages sampled days within a stat-week by day-type stratum and expands to unsampled days. Fast, transparent, and used as the fallback when a BSS fit does not converge.
- **Bayesian State-Space (BSS):** A time-series model that treats daily effort and catch rate as latent quantities evolving smoothly over time. Fills temporal gaps, propagates uncertainty, and produces credible intervals.

Three crabbing populations are estimated independently and summed for the port total:

1. **Shore crabbers** (dock + jetty + beach), effort from gear counts, BSS + PE.
2. **Private boat crabbers**, effort from trailer counts, BSS + PE.
3. **Commercial/charter vessels**, effort from a daily vessel tally, census expansion.

Both BSS models share the same effort model, the PE estimator, the I/E (ingress/egress) handling, and the modular R pipeline in `R_functions/`. They differ in how catch-per-unit-effort (CPUE) is modeled.

---

## Three Approaches

### 1. Pooled CPUE Model (production)

All gear types (pots, ring nets, traps, snares) share a single CPUE process; gear-type catch breakdowns are derived after estimation by applying interview-based proportions to the total. This is the simpler of the two production models and the one to use for a single headline harvest number.

Despite the name, the pooled model is not minimal: it includes adaptive AR(1) temporal resolution (daily/weekly/monthly, selected per fit from effort-data density), a weekend CPUE effect (`B1_C`), effective day length (`L_effective`) estimated as a parameter from the I/E regression, direct I/E crabber-hour integration, a data-driven `R_G` prior, and a divergence-aware convergence gate.

| File | Description |
|---|---|
| `BSS-GH-pooled-CPUE-model.Rmd` | R analysis script |
| `BSS-GH-pooled-CPUE-model-documentation.md` | Technical documentation |
| `stan_models/crab_bss_pooled.stan` | Stan model (single CPUE process) |

### 2. Gear-Resolved CPUE Model (production)

Each gear type gets its own CPUE process with shared AR(1) dynamics, so gear-type catch estimates carry posterior uncertainty directly from the model. Also includes a separate holiday effort effect (`B2`), day-type stratified commercial/charter census expansion, an incomplete-trip filter, and explicit regulatory gear exclusions per sub-season. Use this when you need gear-type catch estimates with uncertainty.

| File | Description |
|---|---|
| `BSS-GH-gear-type-CPUE-model.Rmd` | R analysis script |
| `BSS-GH-gear-type-CPUE-model-documentation.md` | Technical documentation |
| `stan_models/crab_bss_gear_resolved.stan` | Stan model (per-gear CPUE processes) |

### 3. Weather & Tide Covariate Module (experimental)

Tests whether tide and weather covariates (tide phase/range, daytime high-tide timing, wind, wave height) improve the prediction of effort and CPUE, after accounting for weekend/holiday effects. It screens candidate covariates with daily GAMs, fits a covariate-augmented BSS alongside the baseline, and compares them with PSIS-LOO (with a k-fold time-block CV fallback).

This module is **experimental and not a production estimator on its own.** It is currently layered on the pooled model only and shares the pooled pipeline. The augmented Stan model `crab_bss_pooled_weather_adjusted.stan` adds covariate blocks on `mu_E` and `mu_C` and collapses to the baseline pooled model when `K_E = K_C = 0`, so one file serves both the baseline and augmented fits. **The intent is that any covariate shown to help is folded directly into the pooled and gear-resolved models once verified, rather than maintained as a separate track.**

| File | Description |
|---|---|
| `BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd` | R analysis and integration script |
| `BSS-GH-weather-tide-covariates-documentation.md` | Technical documentation |
| `stan_models/crab_bss_pooled_weather_adjusted.stan` | Augmented Stan model (baseline at K=0) |

### Which Should I Use?

Use the **pooled model** for simplicity and a single headline harvest number. Use the **gear-resolved model** when you need gear-type catch estimates with uncertainty, the holiday effect, or the stratified census expansion. The **weather-tide module** is a diagnostic / research tool for deciding whether environmental covariates are worth adding to the production models; it is not used to produce the official harvest estimate.

---

## Stan Models

| File | Used by | Description |
|---|---|---|
| `crab_bss_pooled.stan` | Pooled model | Single pooled CPUE process |
| `crab_bss_gear_resolved.stan` | Gear-resolved model | Per-gear-type CPUE processes, `B2` holiday effect |
| `crab_bss_pooled_weather_adjusted.stan` | Weather-tide module | Pooled model plus covariate blocks; collapses to baseline at `K_E = K_C = 0` |

The `.Rmd` files select their Stan model via the `bss_model_file` (or `bss_model_file_covariates`) parameter. Earlier prototype names (`BSS_crab_model_01/02/03.stan`) are retired.

---

## Quick Start

1. Clone this repository.
2. Place input data in `input_files/`:
   - `effort_combined.csv` (effort counts; re-exported with `QUOTE_ALL`)
   - `interview_combined.csv` (interviews; dates in M/D/YYYY format)
   - `wes_commercial_tally.csv` (daily vessel tally)
   - `ingress_egress.xlsx` (I/E surveys; used for `L_effective` and the temporal correction)
3. Open the desired `.Rmd` file.
4. Set `est_date_start` and `est_date_end` in the parameters block.
5. Run all chunks (or Knit). Output is written to `output/YYYYMMDD/<model>/`.

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here, readxl. The weather-tide module additionally requires mgcv, loo, httr, jsonlite, and geosphere, and reaches NOAA CO-OPS, NDBC, and Iowa State IEM/GSOD endpoints at runtime (results are cached locally under `cache/`).

---

## Output Files

Each run writes to `output/YYYYMMDD/<model>/`. Both models produce PE and combined PE+BSS port totals, monthly estimates, catch by mode and gear type, a per-fit `convergence_report.csv`, a `pe_vs_bss_comparison.csv`, daily BSS effort/catch series, and `run_parameters.txt`. The pooled model additionally writes the I/E and `L_effective` diagnostics (`ie_analysis.csv`, `bss_L_effective_*.csv`, `L_effective_ie_detail.csv`). The gear-resolved model additionally writes gear-type catch with posterior uncertainty (`catch_by_gear_type_detail.csv`) and the monthly/area/mode breakdowns. See each model's documentation for the exact file list.

---

## Season Structure

The 2024-25 season (Sep 16, 2024 to Sep 15, 2025) is split into two independent sub-seasons at the pot-open date (Dec 1):

- **Ring-net only** (Sep 16 to Nov 30, 76 days): ring nets, snares, foldable traps only.
- **All-gear** (Dec 1 to Sep 15, 289 days): all gear including pots.

Each sub-season gets its own BSS fit per population. The split prevents the model from bridging the structural break in effort and CPUE when pots become legal.

---

## Known Issues and Data Notes

- **Interview CSV column mapping:** the `number_of_gear` column maps from column N (not column W) in the raw iForm export, due to a duplicate field name.
- **Effort CSV quoting:** re-export with `QUOTE_ALL` to handle commas in the notes field.
- **Interview dates:** M/D/YYYY (`col_date(format="%m/%d/%Y")`).
- **Boat type typo:** iForm exports "Commerical" (one 'm'), handled by regex.
- **Windows MAX_PATH:** with OneDrive and long paths the output directory may exceed 260 characters; the code detects this and falls back to a short path.
- **Boat all-gear convergence:** the private boat all-gear BSS fit is prone to non-convergence (sparse trailer-count effort series); v6.2 adds dedicated sampler tuning, and the boat component falls back to PE when the fit fails the gate. See the pooled model documentation.

---

## Development History

Versions through v5 are a single shared milestone sequence. Since v5 the pooled and gear-resolved tracks have been versioned independently in their own documentation change logs, and the weather-tide module has its own version line.

| Version | Track | Key Changes |
|---|---|---|
| v1 | shared | Single-population dock-only prototype |
| v2 | shared | Bug fixes (CSV columns, Stan dimensions, output folder) |
| v3 | pooled | Three populations, two sub-seasons, convergence tuning |
| v4 | pooled | Dawn/dusk day length, stat-week PE, census dates, team review |
| v5.0 to v5.4 | gear-resolved | Per-gear CPUE processes, `B2` holiday effect, stratified census, incomplete-trip filter, regulatory gear exclusions, divergence-aware gate, gear-hours boat formulation, R-hat gate tightened to 1.01 |
| v6.0 | pooled | Post-critique upgrades: adaptive AR(1), `L_effective` from I/E, `B1_C` day-type CPUE effect, data-driven `R_G` prior, sparse overdispersion |
| v6.1 | pooled | Divergence-aware convergence gate (`max_divergences`, treedepth warnings, per-fit tuning) |
| v6.2 | pooled | Boat all-gear sampler tuning |
| v6.3 | pooled | Documentation corrections |
| v6.4 | pooled | R-hat convergence threshold tightened to 1.01 (gear-resolved track tightened in step, v5.4) |
| 0.1.0 to 0.1.1 | weather-tide module | Initial build (tide/weather fetch, GAM screen, augmented BSS, PSIS-LOO comparison); reference and file reconciliation |

See each model's documentation change log for details.
