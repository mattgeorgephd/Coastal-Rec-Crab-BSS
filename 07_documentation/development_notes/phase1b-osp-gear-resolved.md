# Phase 1b: OSP boat counts as a second boat effort stream (gear-resolved)

**Applies after:** the Phase 1 patch (`phase1-osp-second-stream.patch`), which itself applies after Phase 0.
**Scope:** the gear-resolved track (`crab_bss_gear_resolved.stan`, `prep_bss_crab_gear.R`) plus the gear driver's diagnostic wiring. Mirrors Phase 1 (pooled) so both production models fit the boat with the OSP stream.
**Status:** behind the same `use_osp_boat_counts` toggle (default FALSE) added in Phase 1. FALSE reproduces the pre-Phase-1b gear-resolved boat EXACTLY. R/Stan not executed in the authoring environment; validate by run.

## What it does

Adds the OSP daily boat-total series as a second observation of the same latent boat effort, adapted to the gear-resolved effort structure (effort is one level split across gears by the `O` share, and the boat observation is against the total `sum_g lambda_E`):

```
trailer:  T_I[i]   ~ NB2( sum_g lambda_E[day]/R_G_boat,               r_E   )
OSP:      OSP_I[i] ~ NB2( (sum_g lambda_E[day]/R_G_boat) * kappa_OSP, r_OSP )
```

At `G = 1` (the production path) `sum_g lambda_E = lambda_E`, so this is identical in form to the pooled Phase 1 stream. `kappa_OSP` is the OSP-to-trailer scale (within-day boat turnover), lognormal prior centered on `osp_scale_prior_mu` (~3, the Phase 0 overlap ratio); own overdispersion `r_OSP`; proper priors unconditionally so a shore fit or an OSP-off boat fit (`OSP_n = 0`) samples them (decoupled), exactly like `R_G_boat` / `sigma_IE`.

## Files changed

- `02_stan_models/crab_bss_gear_resolved.stan`: OSP data block, `osp_scale_prior_mu/sigma`, params `kappa_OSP` + `sigma_r_OSP`, transformed `r_OSP`, model priors + likelihood loop (uses `sum(lambda_E_S[...])`), and `kappa_OSP_out` + `log_lik_osp` in generated quantities.
- `03_R_functions/prep_bss_crab_gear.R`: fetches OSP internally for boat fits (guarded by the toggle), joins to the fit day set, passes `OSP_n/day_OSP/section_OSP/OSP_I` + the kappa prior into `stan_data`. Mirrors `prep_bss_crab_pooled`.
- `01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd`: adds the guarded OSP ingest + overlap-diagnostic call after data load (parity with the pooled driver, since a gear-resolved run does not execute the pooled driver).

No `run_config.R` change: `use_osp_boat_counts`, `osp_scale_prior_mu`, `osp_scale_prior_sigma` were added by the Phase 1 patch and are model-agnostic (each driver reads them).

## Validation gate

Run gear-resolved with `use_osp_boat_counts = TRUE` vs. the FALSE baseline: boat effort intervals should tighten in summer, `kappa_OSP` ~2.5 to 3.5 and tighter than its prior, convergence still passing, boat `C_sum` moving only modestly (densification, not the Phase 2 `f` correction). Cross-check that the pooled and gear-resolved boat still reconcile with OSP on (they agreed at ~43.5k without it; PIPELINE_STATUS fact 1).

## Still open (backlog)

- Phase 2: crabbing fraction `f` (hybrid Beta prior from the I/E classification), the standing Gap A correction (implicit `f = 1` today), applied to both tracks.
- Phase 3: time-varying `f` + let OSP inform `tau_boat`.
