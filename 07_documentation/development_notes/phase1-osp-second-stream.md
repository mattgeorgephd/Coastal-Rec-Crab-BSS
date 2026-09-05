# Phase 1: OSP boat counts as a second boat effort stream (pooled)

**Status 2026-09-08: EXECUTED and in production.** The OSP stream ships on, and since 2026-09-01 it identifies the shared boat turnover `tau_bar` (`osp_scale_is_tau = TRUE`). See `CHANGE_REGISTER.md` A2/A7.

**Applies after:** the Phase 0 patch (`phase0-osp-ingest.patch`).
**Scope:** pooled model only (`crab_bss_pooled.stan`, `prep_bss_crab_pooled.R`). The gear-resolved mirror (`crab_bss_gear_resolved.stan`, `prep_bss_crab_gear.R`) is the identical pattern and is a follow-up (Phase 1b).
**Status:** behind `use_osp_boat_counts` (default FALSE). FALSE reproduces the pre-Phase-1 boat EXACTLY. R/Stan were not executed in the authoring environment; validate by run.

## What it does

Adds the OSP daily boat-total series as a SECOND observation of the same latent boat effort `lambda_E`, alongside the trailer counts:

```
trailer:  T_I[i]   ~ NB2( lambda_E[day]/R_G_boat,               r_E   )   # snapshot of crab boats present
OSP:      OSP_I[i] ~ NB2( (lambda_E[day]/R_G_boat) * kappa_OSP, r_OSP )   # daily total of all private boats
```

`kappa_OSP` is the OSP-to-trailer scale, i.e. the within-day boat turnover (a full-day boat total vs. an instantaneous snapshot). Its lognormal prior is centered on `osp_scale_prior_mu` (~3.0, the OSP/trailer overlap ratio measured by the Phase 0 diagnostic, `1 / mean-per-visit origin slope ~ 0.33`). It is identified by the ~61 overlap days where both streams see the same `lambda_E`; the dense summer OSP series then tightens boat effort on the exact months that dominate the harvest (T1.1), while the trailer stream carries the OSP-dark winter with no swap logic.

`kappa_OSP` and `sigma_r_OSP` carry proper priors unconditionally, so a shore fit or an OSP-off boat fit (`OSP_n = 0`) simply samples them (decoupled, proper), exactly like `R_G_boat` / `sigma_IE`. `OSP_I` enters its own overdispersion `r_OSP` (OSP is a denser, cleaner measurement than the thin trailer series, so it should not be forced to share `r_E`).

## Files changed

- `02_stan_models/crab_bss_pooled.stan`: OSP data block (`OSP_n/day_OSP/section_OSP/OSP_I`), `osp_scale_prior_mu/sigma`, params `kappa_OSP` + `sigma_r_OSP`, transformed `r_OSP`, model priors + likelihood loop, and `kappa_OSP_out` + `log_lik_osp` in generated quantities (PSIS-LOO on the OSP stream).
- `03_R_functions/prep_bss_crab_pooled.R`: fetches OSP internally (boat fits only, guarded by the toggle), joins to the fit's day set, and passes `OSP_n/day_OSP/section_OSP/OSP_I` + the kappa prior into `stan_data`. `OSP_I` is rounded to integer.
- `run_config.R`: `use_osp_boat_counts` (FALSE), `osp_scale_prior_mu` (3.0), `osp_scale_prior_sigma` (0.3).

Note: because the prep fetches OSP itself, no `.Rmd` change is needed for Phase 1 (the Phase 0 patch already added the diagnostic call to the driver).

## What it does NOT do

- No crabbing fraction `f` yet (Phase 2). In Phase 1, `lambda_E` remains the all-boat-scale quantity the model already uses; OSP and trailer are both all-boat, so `kappa_OSP` is pure turnover with no `f` in it. Phase 2 introduces `f` to convert `lambda_E` to crab effort and correct Gap A.
- Not wired into the gear-resolved track (Phase 1b).

## Validation gate (per the plan, §7)

Run pooled with `use_osp_boat_counts = TRUE` vs. the FALSE baseline and confirm: boat effort credible intervals SHRINK in the summer months; `kappa_OSP` posterior is sensible (~2.5 to 3.5) and tighter than its prior; convergence still passes (divergences, R-hat, E-BFMI); the boat `C_sum` moves only modestly (effort densification, not an `f` correction); and `log_lik_osp` gives a clean PSIS-LOO with no extreme Pareto-k. Compare against the confirmed 43,475 boat baseline.
