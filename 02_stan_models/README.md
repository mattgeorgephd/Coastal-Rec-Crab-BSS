# 02_stan_models

Stan model code for the Bayesian State-Space (BSS) estimator. These are called by the drivers in `01_BSS_models/` and `06_diagnostics/` via `rstan::stan(file = here("02_stan_models", <model_file>), ...)`. The driver passes only the filename; the folder is supplied by the `here()` call.

All three models share the same core architecture: an adaptive-resolution AR(1) process for effort and CPUE over `P_n` periods (daily / weekly / monthly, selected in R from effort-data density and mapped to days by `period[d]`), a sparse per-observation effort overdispersion term (`eps_E_H_obs`, one per actual count), I/E-anchored effort integration, and dual reporting of expected catch plus posterior predictive draws. They differ in how CPUE is modeled and in a few effort-side effects.

## Files

| File | Used by | CPUE structure | Distinguishing features |
|---|---|---|---|
| `crab_bss_pooled.stan` | pooled driver | **Single** pooled CPUE process across all gear | `B1_C` weekend CPUE effect; `L_effective` estimated as a parameter (lognormal prior, shore only) from the I/E regression; data-driven `R_G` prior; informative `R_T` Beta prior; non-centered AR(1) initial states |
| `crab_bss_gear_resolved.stan` | gear-resolved driver | **Per-gear** CPUE process (`mu_C` is `[G,S]`, `omega_C` runs over `G*S` with a Cholesky correlation), so gear-type catch carries posterior uncertainty | `B2` holiday effort effect separate from the `B1` weekend effect; gear-deployment effort formulation (`R_G_boat`; boat day length `L = tau` deployment turnover) |
| `crab_bss_pooled_weather_adjusted.stan` | weather-tide module (`06_diagnostics/`) | Pooled, plus covariate blocks `gamma_E`/`gamma_C` on `mu_E`/`mu_C` | Adds `log_lik` for PSIS-LOO; **collapses to `crab_bss_pooled.stan` when `K_E = K_C = 0`** (zero-column covariate matrices), so one file serves both the baseline and augmented fits |

## Effort unit: gear-deployments

Both models measure effort in **gear-deployments** (pot lifts), not a time-denominated unit. Interviews show crab catch is sub-linear in soak time (`crab_per_gear ~ h^0.13`), so crabber-hours and gear-hours are invalid CPUE denominators for pot/trap gear; an earlier gear-hours formulation with `L = 24` inflated boat catch by roughly 2x. For boats, `lambda_E` is gear in the water, `h = number_of_gear` (deployments), day length `L = tau` (a deployment-turnover *parameter*, identified by boat I/E ingress counts when available), and trailer counts map through `R_G_boat`; `E = lambda_E * tau` is gear-deployments per day. Shore likewise runs on gear-deployments (deployment turnover `tau_shore`). The single effort-unit contract lives in `03_R_functions/bss_effort_spec.R`, read by both the BSS prep and the Point Estimator so they always share a unit. See the `F2` header block in `crab_bss_gear_resolved.stan` for the full rationale.

## Selecting a model

The driver's `bss_model_file` (or `bss_model_file_covariates`) parameter names the file. Earlier prototype models (`BSS_creel_model_02_*.stan`, `BSS_crab_model_01/02/03.stan`) from the freshwater-creel lineage are retired and are not in this folder.

## Compiled artifacts

On first build, rstan writes a compiled `*.rds` next to each `*.stan`. These are machine-local and git-ignored (see `.gitignore`); they regenerate automatically when the `.stan` source changes.

## Documentation

Per-model technical documentation is in `07_documentation/` (`BSS-GH-pooled-CPUE-model-documentation.md`, `BSS-GH-gear-type-CPUE-model-documentation.md`, `BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md`).
