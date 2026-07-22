# Recreational Crab Creel Estimation, Grays Harbor / Westport

- **Agency:** Washington Department of Fish and Wildlife (WDFW)
- **Lead:** Matt George
- **Status:** 2024-25 season, first full implementation. Pooled and gear-resolved are the production models; the weather-tide covariate work is an experimental module.

---

## What This Does

Estimates total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area by combining two statistical approaches:

- **Point Estimator (PE):** Stratified expansion. Averages sampled days within a stat-week by day-type stratum and expands to unsampled days. Fast, transparent, and used as the fallback when a BSS fit does not converge.
- **Bayesian State-Space (BSS):** A time-series model that treats daily effort and catch rate as latent quantities evolving smoothly over time. Fills temporal gaps, propagates uncertainty, and produces credible intervals.

Three crabbing populations are estimated independently and summed for the port total:

1. **Shore crabbers** (dock + jetty + beach), effort from gear counts, BSS + PE.
2. **Private boat crabbers**, effort from trailer counts, BSS + PE.
3. **Commercial/charter vessels**, effort from a daily vessel tally, census expansion.

Both BSS models share the same effort model, the PE estimator, the I/E (ingress/egress) handling, the modular R pipeline in `03_R_functions/`, and a single run configuration in `run_config.R`. They differ in how catch-per-unit-effort (CPUE) is modeled.

---

## Repository Layout

The repository is organized into numbered stage folders that follow the order of the pipeline, from analysis driver through model code, helper functions, inputs, outputs, diagnostics, and documentation. Most folders have their own `README.md` with a file inventory.

```text
Coastal-Rec-Crab-BSS/
├── 01_BSS_models/      Production .Rmd analysis drivers (pooled, gear-resolved)
├── 02_stan_models/     Stan model code (.stan) called by the drivers
├── 03_R_functions/     Modular R helpers, auto-sourced by every driver
├── 04_input_files/     Raw season inputs (effort, interviews, tally, I/E)
├── 05_output/          Per-run outputs, one dated folder per run (YYYYMMDD)
├── 06_diagnostics/     Experimental / research .Rmd (weather-tide covariates)
├── 07_documentation/   Technical docs, change logs, equations, instructions
├── run_config.R        Single control surface: user toggles and per-model settings
├── run_estimation.R    Run orchestrator (sources run_config.R)
├── README-R-functions.md   Inventory of the 03_R_functions/ helper library
├── README.md           This file
├── .gitignore
└── Coastal-Rec-Crab-BSS.Rproj
```

| Folder | Contents | README |
|---|---|---|
| `01_BSS_models/` | The two production analysis drivers (`*-pooled-CPUE-model.Rmd`, `*-gear-type-CPUE-model.Rmd`) | [01_BSS_models/README.md](01_BSS_models/README.md) |
| `02_stan_models/` | The three Stan models (pooled, gear-resolved, weather-adjusted) | [02_stan_models/README.md](02_stan_models/README.md) |
| `03_R_functions/` | All R helper functions; the drivers source the whole folder via `purrr::walk` | [README-R-functions.md](README-R-functions.md) |
| `04_input_files/` | `effort_combined.csv`, `interview_combined.csv`, `wes_commercial_tally.csv`, `ingress_egress.xlsx` | [04_input_files/README.md](04_input_files/README.md) |
| `05_output/` | Dated run folders, each with a per-model subfolder of CSVs and plots | [05_output/README.md](05_output/README.md) |
| `06_diagnostics/` | The experimental weather-tide covariate driver | [06_diagnostics/README.md](06_diagnostics/README.md) |
| `07_documentation/` | Per-model documentation, change logs, the rendered equations/landing pages, and the WDFW instruction docs | [07_documentation/README.md](07_documentation/README.md) |

### How paths work (important when moving files)

Every file read or written by the drivers is resolved with `here::here()`, which anchors paths to the repository root (located via the `.Rproj` / `.git` sentinels), **not** to the location of the `.Rmd`. Consequences:

- An `.Rmd` can sit in any subfolder (the drivers live in `01_BSS_models/` and `06_diagnostics/`) and still resolve `here("04_input_files", ...)` correctly, because `here()` walks up to the repo root regardless of the knit working directory.
- The directory **names** inside `here(...)` must match the folder names on disk. The numbered reorganization therefore required updating every `here("R_functions"/"stan_models"/"input_files"/"output", ...)` call to its numbered equivalent (`03_R_functions`, `02_stan_models`, `04_input_files`, `05_output`). If a stage folder is ever renamed again, update the corresponding string in the drivers.
- The weather-tide module also writes a runtime cache to `here("cache", "weather_tide")` at the repo root. This is regenerable and git-ignored, so it is intentionally **not** part of the numbered stage scheme.

---

## Three Approaches

### 1. Pooled CPUE Model (production)

All gear types (pots, ring nets, traps, snares) share a single CPUE process; gear-type catch breakdowns are derived after estimation by applying interview-based proportions to the total. This is the simpler of the two production models and the one to use for a single headline harvest number.

Despite the name, the pooled model is not minimal: it includes adaptive AR(1) temporal resolution (daily/weekly/monthly, selected per fit from effort-data density), a weekend CPUE effect (`B1_C`), effective day length (`L_effective`) estimated as a parameter from the I/E regression, direct I/E crabber-hour integration, a data-driven `R_G` prior, and a divergence-aware convergence gate.

| File | Description |
|---|---|
| `01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd` | R analysis script |
| `07_documentation/BSS-GH-pooled-CPUE-model-documentation.md` | Technical documentation |
| `02_stan_models/crab_bss_pooled.stan` | Stan model (single CPUE process) |

### 2. Gear-Resolved CPUE Model (production)

Each gear type gets its own CPUE process with shared AR(1) dynamics, so gear-type catch estimates carry posterior uncertainty directly from the model. Also includes a separate holiday effort effect (`B2`), day-type stratified commercial/charter census expansion, an incomplete-trip filter, and explicit regulatory gear exclusions per sub-season. Use this when you need gear-type catch estimates with uncertainty.

| File | Description |
|---|---|
| `01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd` | R analysis script |
| `07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md` | Technical documentation |
| `02_stan_models/crab_bss_gear_resolved.stan` | Stan model (per-gear CPUE processes) |

### 3. Weather & Tide Covariate Module (experimental)

Tests whether tide and weather covariates (tide phase/range, daytime high-tide timing, wind, wave height) improve the prediction of effort and CPUE, after accounting for weekend/holiday effects. It screens candidate covariates with daily GAMs, fits a covariate-augmented BSS alongside the baseline, and compares them with PSIS-LOO (with a k-fold time-block CV fallback).

This module is **experimental and not a production estimator on its own.** It is currently layered on the pooled model only and shares the pooled pipeline. The augmented Stan model `crab_bss_pooled_weather_adjusted.stan` adds covariate blocks on `mu_E` and `mu_C` and collapses to the baseline pooled model when `K_E = K_C = 0`, so one file serves both the baseline and augmented fits. **The intent is that any covariate shown to help is folded directly into the pooled and gear-resolved models once verified, rather than maintained as a separate track.**

| File | Description |
|---|---|
| `06_diagnostics/BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd` | R analysis and integration script |
| `07_documentation/BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md` | Technical documentation |
| `02_stan_models/crab_bss_pooled_weather_adjusted.stan` | Augmented Stan model (baseline at K=0) |

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
2. Place input data in `04_input_files/`:
   - `effort_combined.csv` (effort counts; re-exported with `QUOTE_ALL`)
   - `interview_combined.csv` (interviews; dates in M/D/YYYY format)
   - `wes_commercial_tally.csv` (daily vessel tally)
   - `ingress_egress.xlsx` (I/E surveys; used for `L_effective` and the temporal correction)
3. Edit `run_config.R`: choose the `model` ("pooled" or "gear_resolved"), set the season window (`est_date_start`, `est_date_end`), and set any other toggles. As of the 2026-07-11 consolidation, `run_config.R` is the single control surface for a run; you do not edit the `.Rmd` files for a routine run.
4. Launch the run with `source("run_estimation.R")` in RStudio (Source, not Knit) or `Rscript run_estimation.R` from a terminal. You can still knit a model `.Rmd` directly; it sources `run_config.R` automatically when `run_config` is not already present.
5. Output is written to `05_output/YYYYMMDD/<model>-<run_tag>/`.

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here, readxl. The weather-tide module additionally requires mgcv, loo, httr, jsonlite, and geosphere, and reaches NOAA CO-OPS, NDBC, and Iowa State IEM/GSOD endpoints at runtime (results are cached locally under `cache/`).

---

## Output Files

Each run writes to `05_output/YYYYMMDD/<model>-<run_tag>/`. Both models produce PE and combined PE+BSS port totals, monthly estimates, catch by mode and gear type, a per-fit `convergence_report.csv`, a `pe_vs_bss_comparison.csv`, daily BSS effort/catch series, and `run_parameters.txt`. The pooled model additionally writes the I/E and `L_effective` diagnostics (`ie_analysis.csv`, `bss_L_effective_*.csv`, `L_effective_ie_detail.csv`). The gear-resolved model additionally writes gear-type catch with posterior uncertainty (`catch_by_gear_type_detail.csv`) and the monthly/area/mode breakdowns. See each model's documentation for the exact file list.

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
- **Boat all-gear convergence:** the private boat all-gear BSS fit was historically prone to non-convergence (sparse trailer-count effort series). Dedicated sampler tuning (v6.2), the scale-aware convergence gate (v7.0), and moving the boat onto the gear-deployment effort scale (v7.6) have largely resolved this; the boat now typically reports its BSS posterior and falls back to PE only if a fit fails the gate. See the pooled model documentation.

---

## Development History

Versions through v5 are a single shared milestone sequence. Since v5 the pooled and gear-resolved tracks have been versioned independently, each with its own detailed development-history document; the weather-tide module has its own version line. The table below is a one-line-per-milestone summary. See the two development-history documents for the full change log with working notes:

- `07_documentation/BSS-GH-pooled-CPUE-model-development-history.md` (pooled, through v7.9)
- `07_documentation/BSS-GH-gear-type-CPUE-model-development-history.md` (gear-resolved, through v5.6)

For the current state and the prioritized backlog (what is done and what remains), see the single living status document: `07_documentation/development_notes/PIPELINE_STATUS.md`.

| Version | Track | Key Changes |
|---|---|---|
| v1 to v2 | shared | Single-population dock-only prototype; bug fixes (CSV columns, Stan dimensions, output folders) |
| v3 to v4 | pooled | Three populations, two sub-seasons, convergence tuning; dawn/dusk day length, stat-week PE, census dates, team review |
| v5.0 to v5.5 | gear-resolved | Per-gear CPUE processes, `B2` holiday effect, stratified census, incomplete-trip filter, regulatory gear exclusions; empirical `pi_gear`; divergence-aware then R-hat < 1.01 gate; boat and shore moved onto the gear-deployment effort scale (v5.5) |
| v6.0 to v7.4 | pooled | Post-critique upgrades (adaptive AR(1), `L_effective` from I/E, `B1_C`, data-driven `R_G`); the convergence-debugging arc (divergence gate, boat tuning, non-centered AR, marginalized NB, scale-aware gate); extended diagnostics and PSIS-LOO. Method v1.0 = code v7.4 |
| v7.5 to v7.8 | pooled | Backlog fixes (incomplete-trip filter, CPUE diagnostics, `collapse_mu_hier` lever); boat (v7.6) then shore (v7.7) moved onto the gear-deployment effort scale; behavior-preserving repository refactor and the shore-PE completion fix (v7.8) |
| 0.1.0 to 0.1.1 | weather-tide module | Initial build (tide/weather fetch, GAM screen, augmented BSS, PSIS-LOO comparison); reference and file reconciliation |

See each model's development-history document for details.

---

## License and attribution

Coastal-Rec-Crab-BSS is free software, licensed under the **GNU General Public License, version 3 (GPL-3.0)**; see [`LICENSE`](LICENSE). Copyright (C) 2024-2026 Washington Department of Fish and Wildlife.

This pipeline is a **derivative work of [CreelEstimates](https://github.com/dfw-wa/CreelEstimates)**, the WDFW freshwater creel estimation framework (also GPL-3.0). The BSS and PE methodology, the project structure, and substantial portions of the R and Stan code originate in CreelEstimates and remain copyright their authors under GPL-3.0; the adaptation to recreational Dungeness crab is by WDFW. See [`NOTICE`](NOTICE) for the full attribution. Every source file carries a GPL-3.0 header with the copyright and this attribution. If you use or redistribute this software, retain the license, the per-file headers, and this attribution, and cite CreelEstimates as the upstream source.
