# 02_stan_models

Stan model code for the Bayesian State-Space (BSS) estimator. These are called by the drivers in `01_BSS_models/` (and by the batch runners in `06_diagnostics/`, which render those drivers) via `rstan::stan(file = here("02_stan_models", <model_file>), ...)`. The driver passes only the filename; the folder is supplied by the `here()` call.

Both models share the same core architecture: an adaptive-resolution AR(1) process for effort and CPUE over `P_n` periods (daily / weekly / biweekly / monthly; selected in R from effort-data density, then coarsened by the per-population `ar_max_resolution` cap, or forced by the ladder/`ar_force`), effort overdispersion marginalized as NB2 (`r_E`; the per-observation `eps_E_H_obs` latents were removed in B1.5), I/E-anchored effort integration, and dual reporting of expected catch plus posterior predictive draws. They differ in how CPUE is modeled and in a few effort-side effects.

## Files

| File | Used by | CPUE structure | Distinguishing features |
|---|---|---|---|
| `crab_bss_pooled.stan` | pooled driver | **Single** pooled CPUE process across all gear | `B1_C` weekend and `B2_C` holiday CPUE effects; gear-deployment effort (`R_G_boat` trailer mapping; `R_T` is the retired legacy expansion); **shared turnover** `tau_bar` (`shared_tau`, one estimated turnover with a fixed day spread, replacing per-day draws); **zero-inflated catch likelihood** (`zi_catch` / `theta_C`, per-fit via `catch_zi_populations`, season total scaled by `1 - theta_C` in generated quantities); other-fishery opener effort covariates (`K_open`/`X_open`/`B_open`, off in production); OSP second boat-effort likelihood (scaled by `L[day]` when `osp_scale_is_tau`, else `kappa_OSP`); crabbing fraction `f_crab` as a **dynamic logit random walk** across month strata observed per day by the sampler contacts and OSP's crabbing-only counts (`crab_fraction_dynamic`; the legacy per-stratum Beta/Binomial with the OSP lower bound `f_lower` stays as the gated control), scaling boat generated quantities only; non-centered AR(1) initial states; per-stream `log_lik` for PSIS-LOO |
| `crab_bss_gear_resolved.stan` | gear-resolved driver | **Per-gear** CPUE process (`mu_C` is `[G,S]`, `omega_C` runs over `G*S` with a Cholesky correlation), so gear-type catch carries posterior uncertainty | `B2` holiday effort effect separate from the `B1` weekend effect; gear-deployment effort formulation (`R_G_boat`; boat day length `L = tau` deployment turnover); OSP second boat-effort likelihood (scaled by `kappa_OSP`, or `tau_boat` when `osp_scale_is_tau`); `f_crab` scaling in the boat generated quantities only (decoupled from sampling, CPUE invariant), with the same dynamic-f block as the pooled model; since 2026-09-13 the same zero-inflated catch block (`zi_catch` / `theta_C`, the D6 port), off as shipped because `catch_zi_tracks` ships `"pooled"`; since B46 (2026-09-28) the shore `R_G` prior is data (`R_G_prior_mu`, `R_G_prior_sigma`, resolved as in the pooled prep) rather than the literal `lognormal(log(1.3), 0.3)` |

## Effort unit: gear-deployments

**Both models' `C` / `C_sum` generated quantities are a Poisson draw** on the expected rate and carry neither `r_C` nor `theta_C`. Since B46 (2026-09-28) the drivers do not report them: `Predictive_Catch` is rebuilt after the fit from the fitted zero-inflated NB (`03_R_functions/bss_predictive_catch.R`), which leaves the Stan code, and so every fit, unchanged. Both models also `reject()` `osp_crab_lower = 1` with a dynamic `f` and no dynamic `c`, which the R preps already refuse.

Both models measure effort in **gear-deployments** (a piece of gear a crabber had in the water on the trip, `number_of_gear`; NOT a pot lift: repeat checks of one gear slot are not extra deployments), not a time-denominated unit. Interviews show crab catch is sub-linear in soak time (`crab_per_gear ~ h^0.13`), so crabber-hours and gear-hours are invalid CPUE denominators for pot/trap gear; an earlier gear-hours formulation with `L = 24` inflated boat catch by roughly 2x. For boats, `lambda_E` is gear in the water, `h = number_of_gear` (deployments), day length `L = tau` (a deployment-turnover *parameter*, identified by boat I/E ingress counts when available), and trailer counts map through `R_G_boat`; `E = lambda_E * tau` is gear-deployments per day. Shore likewise runs on gear-deployments (deployment turnover `tau_shore`). The single effort-unit contract lives in `03_R_functions/bss_effort_spec.R`, read by both the BSS prep and the Point Estimator so they always share a unit. See the `F2` header block in `crab_bss_gear_resolved.stan` for the full rationale.

## OSP boat counts and the crabbing fraction f

Two boat-side features are opt-in via `run_config` and behavior-neutral when off. When `use_osp_boat_counts = TRUE`, the OSP daily private-boat total enters as a second effort observation on the boat latent `lambda_E`, scaled by an OSP within-day turnover `kappa_OSP` (or by `tau_boat` when `osp_scale_is_tau = TRUE`, so the dense OSP series identifies boat turnover). When `use_crab_fraction = TRUE`, a crabbing fraction `f_crab` (the trailer and OSP counts are all boats, including salmon and tuna trips) scales boat effort and catch in the boat generated quantities only, leaving the sampling model and the CPUE process unchanged. With both toggles off, each model reproduces its prior trailer-only behavior.

**The dynamic f (2026-09-08, review item 1B).** Under `crab_fraction_dynamic = TRUE` (the shipped value) both models replace the per-stratum Beta/Binomial `f` with a random walk on the log-odds of `f` across the month strata in chronological order (`f_walk_prev` / `f_walk_gap` encode the order; the first stratum of a chain sits on a weak `N(f_level_mu, f_level_sd)` level prior; the step SD `sigma_f` is half-normal and learned; `f_walk_df > 0` gives Student-t steps so the two seasonal regime jumps pass without loosening every quiet month). It is observed per DAY by the sampler boat contacts (`CFI_n` rows: every private boat approached at the launch, crabbing or not, a combo trip counted as crabbing; beta-binomial with concentration `cfi_kappa`) and, when `WBL_boat_counts.xlsx` carries OSP's crabbing-only column, by the OSP counts through `f x (1 - combo_c)`, `combo_c` the combo-trip share among crabbing boats (identified only where both streams cover a stratum). Both denominators are observed boat counts, so `f` stays out of the effort and CPUE likelihoods and the joint posterior factorizes exactly as before. The Phase 3 / improvement 8 construction stays in both files, gated on `crab_fraction_dynamic = 0`, as the Phase A control. Reported: `f_crab_out[k]` (labelled per stratum in `crab_fraction_strata_<fit>.csv`), `sigma_f_out`, `cfi_kappa_out`, `combo_c_out[k]`, `sigma_c_out`, `cfc_kappa_out` (0.0 when not in the model; see the decoupled flags).

**The combo-trip share is observed (2026-09-09).** The interview workbook now carries every contacted boat's trip type, so `c[k]`, the share of crabbing boats that also fished another fishery, is a per-day beta-binomial observation (`CFC_n` rows: combos of the typed crabbing boats) and carries its own logit random walk over the same strata and chains as `f` (`combo_dynamic`, `z_c`, `sigma_c`, `cfc_kappa`), replacing the 2026-09-08 scalar `combo_c ~ Beta`. The OSP stream reads `f x (1 - c[k])`, and `f_lower_out[k] = f(1 - c)` is the implied crab-only share, the prediction of OSP's column. On the 2024-25 boat all-gear fit `c` runs from about 0.02 in Dec-Feb to 0.4-0.8 in the salmon and bottomfish months (`sigma_c` about 1.0), which is why OSP's crab-only count alone could never have identified `f` in summer.

## The zero-inflated catch likelihood: in both models, on the pooled track only as shipped

As of 2026-09-07 the pooled model's catch likelihood is a zero-inflated negative binomial on the populations named by `catch_zi_populations` (shore in production): `P(0) = theta_C + (1 - theta_C) * NB2(0)`, with Stan's `log_lik` carrying the mixture and the season total scaled by `(1 - theta_C)` in generated quantities so enabling the feature is not itself an inflation. **`crab_bss_gear_resolved.stan` has carried the same block since 2026-09-13** (the D6 port; until then it had none and the gear track silently ignored `estimate_catch_zi`). Which TRACKS fit it is `catch_zi_tracks` in `run_config.R`, which ships `"pooled"`, so `prep_bss_crab_gear()` emits `zi_catch = 0` and the gear track still fits plain NB2 as shipped; the OFF path is bit-identical to the pre-port model (11,021 parameter rows). The two tracks therefore still differ in the shore catch likelihood as shipped (worth about -0.3% on the pooled shore component). D6 in `07_documentation/development_notes/CHANGE_REGISTER.md` says adopt, pending one render at the matched configuration. Changing `estimate_catch_zi` (like `razor_dig_mode` and `estimate_cpue_density`) forces a Stan recompile.

## The level hierarchy with one section (A31, 2026-09-29)

Each model writes the effort and CPUE levels as a hierarchy over sections,
`mu[g,s] = mu_mu[g] + eps_mu[g,s] * sigma_mu`. With ONE section, which is every production fit,
only the sum enters the likelihood, so `sigma_mu` is unidentified and the product is a funnel;
the pooled shore all-gear fit's divergences sat on the effort one (D33). The gear-resolved
model has collapsed both levels at `S == 1` since v6.0 (its P2 block). **The pooled model
collapses them per level since 2026-09-29**, from the data int `mu_hier_collapse_single`
(0 none, the pre-A31 model exactly; 1 the effort level, shipped; 2 both), set by
`run_config$mu_hier_collapse_single`, so no recompile. A collapsed level's `eps_mu` is zero-size
and its `sigma_mu` keeps a proper prior that enters nothing; it reports its prior and
`bss_decoupled_reasons()` flags it. `collapse_mu_hier = 1` still forces both collapses on a
multi-section fit. The container validation, and why `"effort"` rather than `"both"` ships, are
in CHANGE_REGISTER A31.

## Selecting a model

The driver's `bss_model_file` parameter names the file (the `bss_model_file_covariates` key belonged to the weather-tide driver and went with it on 2026-09-13). Earlier prototype models (`BSS_creel_model_02_*.stan`, `BSS_crab_model_01/02/03.stan`) from the freshwater-creel lineage are retired and are not in this folder.

## Compiled artifacts

On first build, rstan writes a compiled `*.rds` next to each `*.stan`. These are machine-local and git-ignored (see `.gitignore`); they regenerate automatically when the `.stan` source changes.

## Documentation

Per-model technical documentation is in `07_documentation/` (`BSS-GH-pooled-CPUE-model-documentation.md` (Method v2.0) and `BSS-GH-gear-type-CPUE-model-documentation.md`). The removed weather-tide module's method document is archived as `07_documentation/archive/weather-tide-covariate-module-REMOVED.md`.

## Removed 2026-09-13

`crab_bss_pooled_weather_adjusted.stan`, the weather-tide module's fork of the pooled model,
was deleted with the module (CHANGE_REGISTER A29). It had drifted about **40 data variables**
behind `crab_bss_pooled.stan` (missing the whole `shared_tau` block, the OSP effort stream, the
`crab_fraction_*` block and the zero-inflation flags), so the claim that it collapses to the
pooled model at `K_E = K_C = 0` had stopped being true on 2026-09-02. **There are now two Stan
models and both are production.** The covariate finding is kept at
`07_documentation/WEATHER_COVARIATE_ANALYSIS.md`.
