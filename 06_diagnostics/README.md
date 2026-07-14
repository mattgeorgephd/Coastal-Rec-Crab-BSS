# 06_diagnostics

The **experimental weather-tide covariate module**. This folder holds the one research driver that tests whether environmental covariates improve the effort and CPUE predictions of the production pooled model. It is **not a production estimator** and does not produce the official harvest number.

For the production models, see [`01_BSS_models/`](../01_BSS_models/README.md); for the project overview, the [root README](../README.md).

## Files

| File | Role |
|---|---|
| `BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd` | The covariate module driver (module v0.2.x). Layered on the **pooled** model only. Screens candidate tide/weather covariates with daily GAMs, fits a covariate-augmented BSS alongside the baseline, and compares them with PSIS-LOO (with a leave-one-week-out block-CV fallback for true sampling gaps). |
| `README.md` | This file. |

The augmented Stan model it fits, `crab_bss_pooled_weather_adjusted.stan`, lives in [`02_stan_models/`](../02_stan_models/README.md); it adds covariate blocks (`gamma_E` on `mu_E`, `gamma_C` on `mu_C`) and **collapses exactly to `crab_bss_pooled.stan` when `K_E = K_C = 0`**, so one file serves both the baseline and augmented fits.

## How it runs

- As part of a full run: set `run_weather <- TRUE` in `run_config.R` (or pass `--weather`) with `model = "pooled"`, and `run_estimation.R` renders this module **after** the pooled model. It is only valid with the pooled model — the module reuses the pooled run's in-memory objects (`dwg`, `ie_data`, `L_eff_model`, …) via the shared render environment ("Option A" hand-off), so it cannot run without a preceding pooled run in the same session.
- Standalone: knit `BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd` directly; its setup chunk auto-sources `run_config.R` when `run_config` is not already present.

Outputs land in `05_output/<YYYYMMDD>/pooled-CPUE-covariates/` (paired `*_baseline` / `*_covariates` files, GAM smooths, and the covariate-vs-baseline LOO comparison). At runtime the module reaches NOAA CO-OPS, NDBC, and Iowa State IEM/GSOD endpoints for tide and weather series; results are cached under `cache/weather_tide/` at the repo root (regenerable and git-ignored).

## Status

**Experimental and currently stale.** Per `07_documentation/development_notes/PIPELINE_STATUS.md`, the module tracks an older pooled engine (~v6.9 parity, pre-deployment-scale) and has not been re-run against the current Stan models, so its boat-magnitude outputs must not be cited. Its committed conclusion is that weather/tide covariates are **excluded** — see the decision record in `07_documentation/WEATHER_COVARIATE_ANALYSIS.md` and Section 17 of the pooled-model documentation. The intent is that any covariate later shown to help is folded directly into the production pooled and gear-resolved models, rather than maintained as a separate track.

## Documentation

Technical write-up: `07_documentation/BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md`. Decision record (exclusion finding): `07_documentation/WEATHER_COVARIATE_ANALYSIS.md`.
