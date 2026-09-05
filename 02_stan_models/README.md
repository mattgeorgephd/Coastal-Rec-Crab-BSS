# 02_stan_models

Stan model code for the Bayesian State-Space (BSS) estimator. These are called by the drivers in `01_BSS_models/` and `06_diagnostics/` via `rstan::stan(file = here("02_stan_models", <model_file>), ...)`. The driver passes only the filename; the folder is supplied by the `here()` call.

All three models share the same core architecture: an adaptive-resolution AR(1) process for effort and CPUE over `P_n` periods (daily / weekly / biweekly / monthly; selected in R from effort-data density, then coarsened by the per-population `ar_max_resolution` cap, or forced by the ladder/`ar_force`), a sparse per-observation effort overdispersion term (`eps_E_H_obs`, one per actual count), I/E-anchored effort integration, and dual reporting of expected catch plus posterior predictive draws. They differ in how CPUE is modeled and in a few effort-side effects.

## Files

| File | Used by | CPUE structure | Distinguishing features |
|---|---|---|---|
| `crab_bss_pooled.stan` | pooled driver | **Single** pooled CPUE process across all gear | `B1_C` weekend and `B2_C` holiday CPUE effects; gear-deployment effort (`R_G_boat` trailer mapping; `R_T` is the retired legacy expansion); **shared turnover** `tau_bar` (`shared_tau`, one estimated turnover with a fixed day spread, replacing per-day draws); **zero-inflated catch likelihood** (`zi_catch` / `theta_C`, per-fit via `catch_zi_populations`, season total scaled by `1 - theta_C` in generated quantities); other-fishery opener effort covariates (`K_open`/`X_open`/`B_open`, off in production); OSP second boat-effort likelihood (scaled by `L[day]` when `osp_scale_is_tau`, else `kappa_OSP`); crabbing fraction `f_crab` with the optional OSP crab-only lower bound (`f_lower`), scaling boat generated quantities only; non-centered AR(1) initial states; per-stream `log_lik` for PSIS-LOO |
| `crab_bss_gear_resolved.stan` | gear-resolved driver | **Per-gear** CPUE process (`mu_C` is `[G,S]`, `omega_C` runs over `G*S` with a Cholesky correlation), so gear-type catch carries posterior uncertainty | `B2` holiday effort effect separate from the `B1` weekend effect; gear-deployment effort formulation (`R_G_boat`; boat day length `L = tau` deployment turnover); OSP second boat-effort likelihood (scaled by `kappa_OSP`, or `tau_boat` when `osp_scale_is_tau`); `f_crab` scaling in the boat generated quantities only (decoupled from sampling, CPUE invariant) |
| `crab_bss_pooled_weather_adjusted.stan` | weather-tide module (`06_diagnostics/`) | Pooled, plus covariate blocks `gamma_E`/`gamma_C` on `mu_E`/`mu_C` | Adds `log_lik` for PSIS-LOO; **collapses to `crab_bss_pooled.stan` when `K_E = K_C = 0`** (zero-column covariate matrices), so one file serves both the baseline and augmented fits |

## Effort unit: gear-deployments

Both models measure effort in **gear-deployments** (pot lifts), not a time-denominated unit. Interviews show crab catch is sub-linear in soak time (`crab_per_gear ~ h^0.13`), so crabber-hours and gear-hours are invalid CPUE denominators for pot/trap gear; an earlier gear-hours formulation with `L = 24` inflated boat catch by roughly 2x. For boats, `lambda_E` is gear in the water, `h = number_of_gear` (deployments), day length `L = tau` (a deployment-turnover *parameter*, identified by boat I/E ingress counts when available), and trailer counts map through `R_G_boat`; `E = lambda_E * tau` is gear-deployments per day. Shore likewise runs on gear-deployments (deployment turnover `tau_shore`). The single effort-unit contract lives in `03_R_functions/bss_effort_spec.R`, read by both the BSS prep and the Point Estimator so they always share a unit. See the `F2` header block in `crab_bss_gear_resolved.stan` for the full rationale.

## OSP boat counts and the crabbing fraction f

Two boat-side features are opt-in via `run_config` and behavior-neutral when off. When `use_osp_boat_counts = TRUE`, the OSP daily private-boat total enters as a second effort observation on the boat latent `lambda_E`, scaled by an OSP within-day turnover `kappa_OSP` (or by `tau_boat` when `osp_scale_is_tau = TRUE`, so the dense OSP series identifies boat turnover). When `use_crab_fraction = TRUE`, a crabbing fraction `f_crab` (the trailer and OSP counts are all boats, including salmon and tuna trips) scales boat effort and catch in the boat generated quantities only, leaving the sampling model and the CPUE process unchanged. With both toggles off, each model reproduces its prior trailer-only behavior.

## The zero-inflated catch likelihood is POOLED-ONLY

As of 2026-09-07 the pooled model's catch likelihood is a zero-inflated negative binomial on the populations named by `catch_zi_populations` (shore in production): `P(0) = theta_C + (1 - theta_C) * NB2(0)`, with Stan's `log_lik` carrying the mixture and the season total scaled by `(1 - theta_C)` in generated quantities so enabling the feature is not itself an inflation. **`crab_bss_gear_resolved.stan` has no zero-inflation block** and its data prep never emits `zi_catch`, so the gear track silently ignores `estimate_catch_zi` and fits plain NB2. The two tracks therefore differ in the shore catch likelihood (worth about -0.3% on the pooled shore component); porting the block is a Stan edit plus a recompile, tracked as D6 in `07_documentation/development_notes/CHANGE_REGISTER.md`. Changing `estimate_catch_zi` (like `razor_dig_mode` and `estimate_cpue_density`) forces a Stan recompile.

## Selecting a model

The driver's `bss_model_file` (or `bss_model_file_covariates`) parameter names the file. Earlier prototype models (`BSS_creel_model_02_*.stan`, `BSS_crab_model_01/02/03.stan`) from the freshwater-creel lineage are retired and are not in this folder.

## Compiled artifacts

On first build, rstan writes a compiled `*.rds` next to each `*.stan`. These are machine-local and git-ignored (see `.gitignore`); they regenerate automatically when the `.stan` source changes.

## Documentation

Per-model technical documentation is in `07_documentation/` (`BSS-GH-pooled-CPUE-model-documentation.md`, `BSS-GH-gear-type-CPUE-model-documentation.md`, `BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md`).
