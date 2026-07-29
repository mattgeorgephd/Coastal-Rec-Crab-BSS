# Phase 3: time-varying f (by stratum) + optional OSP-informs-tau

**Applies after:** the Phase 2b patch. **Scope:** both tracks. **Status:** behavior-neutral by default: `crab_fraction_strata = "none"` reproduces the Phase 2 scalar f exactly, and `osp_scale_is_tau = FALSE` reproduces the Phase 1 OSP scale exactly. Opt into each feature per run. R/Stan not executed here; validate by run.

## Part 1: time-varying f by stratum (the main deliverable)

Phase 2 applied a single scalar `f`. But the crab share of boats very plausibly varies within the season (lower in the summer tuna/salmon peak, higher in the fall/late shoulder), and the harvest is summer-dominated, so a single annual `f` over-states summer crab effort, re-introducing a bias in the place the model is weakest.

Phase 3 makes `f` PER STRATUM. `crab_fraction_strata` selects the stratification: `"none"` (one scalar = Phase 2), `"month"`, `"day_type"`, or `"month_day_type"`. Each stratum gets its own `Beta(set*kappa, (1-set)*kappa)` prior, updated by that stratum's I/E crab-vs-total classification counts (Binomial) once it clears `crab_fraction_min_obs`. `n_f_strata = 1` reproduces Phase 2 exactly, so the default path is unchanged.

The Stan `f_crab` scalar became `vector[n_f_strata] f_crab`, and each day scales by `f_crab[f_stratum[d]]` in the generated quantities (still decoupled from sampling). `f_stratum[d]` and the per-stratum priors/counts come from `crab_fraction_stan_data()`, which labels every day and every WBL classification row by one date-based rule (`crab_fraction_strata_labels()`), so the day set and the classification rows are stratified consistently. The PE takes a per-day point `f` (`crab_fraction_point_day()`, the per-stratum Beta posterior mean).

`fetch_ie_data()` now emits the classification PER ROW (`attr(ie_data, "crab_fraction_rows")`, a tibble of `event_date, boats_crabbing, boats_total`) so the helper can aggregate by any stratification; the driver lifts it into `params$crab_fraction_rows`. Empty (the pilot is in progress) -> every stratum uses the set value, so `"month"` with no data equals the scalar 0.3, just reported per month. **Field note:** for time-varying `f` to actually bind, the WBL egress pilot must classify enough boats PER STRATUM (per month, or per month x day-type), not just in aggregate.

## Part 2: OSP-informs-tau (opt-in; carries a real tension)

`osp_scale_is_tau` (default FALSE) makes the OSP observation mean use `L` (= `tau_boat`) as the within-day turnover instead of the free `kappa_OSP`, so the dense OSP series (148 days) would identify the boat turnover directly and close GR-12 (today `tau_boat` rests on 2 winter WBL I/E days + its prior).

**The tension, stated plainly:** the Phase 0 overlap put the OSP/trailer ratio at ~2.7 (both series are all-boat, so the ratio is the within-day boat turnover), while `tau_boat`'s prior is ~1.2. If those are the same quantity, turning this on pulls `tau_boat` from 1.2 toward ~2.7 and roughly DOUBLES the boat effort/catch, partly offsetting the Phase 2 `f` reduction. Alternatively the ratio may be a trailer-snapshot-timing artifact (if the crab-creel trailer count is taken near the daily peak rather than at a random instant, OSP/trailer overstates turnover). These are opposite conclusions with a ~2x swing on the boat, so this MUST be resolved by a run, not shipped on by default.

Safe first step (no toggle): run Phase 1 with OSP on and read `kappa_OSP_out`. If its posterior sits near ~2.7 and far from `tau_boat`, that is the flag to investigate the snapshot-timing question before enabling `osp_scale_is_tau`.

## Files changed

- `02_stan_models/crab_bss_pooled.stan`, `crab_bss_gear_resolved.stan`: `f_crab` scalar -> `vector[n_f_strata]` with `f_stratum[D]` and per-stratum priors/counts; the `osp_scale_is_tau` conditional in the OSP likelihood + log_lik.
- `03_R_functions/crab_fraction.R`: rewritten for per-stratum (`crab_fraction_strata_labels`, per-stratum `crab_fraction_stan_data(is_shore, days, params)`, per-day `crab_fraction_point_day`).
- `03_R_functions/bss_day_length.R`: `fetch_ie_data()` emits per-row classification.
- `prep_bss_crab_pooled.R` / `_gear.R`: pass `n_f_strata` / `f_stratum` / `osp_scale_is_tau`.
- `run_pe_pooled.R` / `run_pe_gear.R`: per-day point `f`.
- both drivers: lift `crab_fraction_rows` into `params`.
- `run_config.R`: `crab_fraction_strata` ("none"), `osp_scale_is_tau` (FALSE).

## Validation gate

- `crab_fraction_strata = "none"` + `osp_scale_is_tau = FALSE` must reproduce the Phase 2b run bit-for-bit (default path).
- `crab_fraction_strata = "month"` with no pilot data: identical totals to the scalar 0.3 run, but `f_crab_out` is a length-(#months) vector all ~0.3.
- With pilot classification loaded: per-month `f_crab_out` tracks each month's `n_crab/n_total` shrunk toward 0.3; boat summer months move relative to shoulders.
- `osp_scale_is_tau = TRUE` (OSP on): expect `tau_boat` / boat effort to move materially (see the tension above); compare against the FALSE run and against `kappa_OSP_out`.
