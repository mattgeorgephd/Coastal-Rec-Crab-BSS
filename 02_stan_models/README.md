# 02_stan_models

Stan model code for the Bayesian State-Space (BSS) estimator. These are called by the drivers in `01_BSS_models/` and `06_diagnostics/` via `rstan::stan(file = here("02_stan_models", <model_file>), ...)`. The driver passes only the filename; the folder is supplied by the `here()` call.

All three models share the same core architecture: an adaptive-resolution AR(1) process for effort and CPUE over `P_n` periods (daily / weekly / monthly, selected in R from effort-data density and mapped to days by `period[d]`), a sparse per-observation effort overdispersion term (`eps_E_H_obs`, one per actual count), I/E-anchored effort integration, and dual reporting of expected catch plus posterior predictive draws. They differ in how CPUE is modeled and in a few effort-side effects.

## Files

| File | Used by | CPUE structure | Distinguishing features |
|---|---|---|---|
| `crab_bss_pooled.stan` | pooled driver | **Single** pooled CPUE process across all gear | `B1_C` weekend CPUE effect; `L_effective` estimated as a parameter (lognormal prior, shore only) from the I/E regression; data-driven `R_G` prior; informative `R_T` Beta prior; non-centered AR(1) initial states |
| `crab_bss_gear_resolved.stan` | gear-resolved driver | **Per-gear** CPUE process (`mu_C` is `[G,S]`, `omega_C` runs over `G*S` with a Cholesky correlation), so gear-type catch carries posterior uncertainty | `B2` holiday effort effect separate from the `B1` weekend effect; gear-hours boat formulation (`R_G_boat`, `L = 24` for boats) |
| `crab_bss_pooled_weather_adjusted.stan` | weather-tide module (`06_diagnostics/`) | Pooled, plus covariate blocks `gamma_E`/`gamma_C` on `mu_E`/`mu_C` | Adds `log_lik` for PSIS-LOO; **collapses to `crab_bss_pooled.stan` when `K_E = K_C = 0`** (zero-column covariate matrices), so one file serves both the baseline and augmented fits |

## Gear-hours formulation (boats)

The gear-resolved model uses gear-hours for boat populations rather than crabber-hours: `lambda_E` is gear in the water, the CPUE denominator `h` is gear-hours, day length `L = 24` (pots fish continuously while the trailer is at the ramp), and trailer counts map through `R_G_boat` (gear per boat group). This corrects a unit mismatch that previously underestimated boat catch by roughly 2x. The shore side is unchanged (crabbers, crabber-hours, civil-twilight day length).

## Selecting a model

The driver's `bss_model_file` (or `bss_model_file_covariates`) parameter names the file. Earlier prototype models (`BSS_creel_model_02_*.stan`, `BSS_crab_model_01/02/03.stan`) from the freshwater-creel lineage are retired and are not in this folder.

## Known issue: stale header in `crab_bss_gear_resolved.stan`

> **Flag (documentation only, not a path issue).** The header comment block inside `crab_bss_gear_resolved.stan` still reads "Pooled CPUE Crab Creel Model ... Single CPUE process shared across all gear types" and carries an old alternate filename `(crab_bss_pooled_gearhours.stan)`. This is a leftover from the pooled-gearhours prototype the file was derived from and now **contradicts the model's own code**, which resolves CPUE by gear (`matrix[G,S] mu_C`, `mu_mu_C[G]`, per-gear `omega_C` over `G*S`, and a gear-indexed interview likelihood via `gear_IntC`). The code, the root README, and the model documentation all agree the model is gear-resolved; only this header comment is wrong. Recommend updating the header comment so the file's self-description matches its behavior. No functional impact.

## Compiled artifacts

On first build, rstan writes a compiled `*.rds` next to each `*.stan`. These are machine-local and git-ignored (see `.gitignore`); they regenerate automatically when the `.stan` source changes.

## Documentation

Per-model technical documentation is in `07_documentation/` (`BSS-GH-pooled-CPUE-model-documentation.md`, `BSS-GH-gear-type-CPUE-model-documentation.md`, `BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md`).
