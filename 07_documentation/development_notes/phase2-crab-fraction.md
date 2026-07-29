# Phase 2: directed-crabbing fraction f (both tracks)

**Applies after:** the Phase 1b patch. **Scope:** both production models. **Status:** ships ON (`use_crab_fraction = TRUE`, `crab_fraction_set = 0.3`) per Matt's choice. `use_crab_fraction = FALSE` reproduces the pre-Phase-2 boat exactly. R/Stan not executed in the authoring environment; validate by run.

## What f is and why it matters

`f` = the share of private boats at the Westport launch that are crabbing (vs. targeting tuna/salmon/halibut/bottomfish). The effort series the model uses (trailer counts, and the Phase 1 OSP boat totals) count ALL private boats, so crab effort = `f` x all-boat effort. The pre-Phase-2 model has no `f`, i.e. it implicitly sets `f = 1` (every boat crabbing), which biases the boat catch high. This is the standing "Gap A" from the plan. Phase 2 introduces `f`; it is the phase that actually moves the boat number, and it moves it DOWN.

**With the shipped set value f = 0.3, the boat BSS catch and effort scale to ~30% of their f = 1 values** (CPUE is unchanged, since f cancels in catch/effort). 0.3 is a starting value chosen pending the pilot, not a measured estimate.

## Design (hybrid; matches R_G / tau_boat / L_effective)

`f` is a parameter with an informative Beta prior centered on the set value, applied ONLY in the Stan generated quantities (scaling boat `E` and `lambda_Ctot` -> `C`). It is decoupled from the effort/CPUE sampling, so it adds no identifiability or convergence risk (f does not appear in any effort or catch likelihood; the catch interviews remain crab-boat CPUE, and trailer/OSP remain all-boat effort).

Set-value fallback is automatic:

- **I/E classification thin/absent** (`n_total < crab_fraction_min_obs`, the current state): `f ~ Beta(set*kappa, (1-set)*kappa)`, prior only, centered on the set value. Uncertainty is carried.
- **Enough classified boats:** a `Binomial(n_total, f)` likelihood on the crab-vs-total counts updates `f`, so data take over as the pilot accumulates.
- **`crab_fraction_fixed` set:** `f` is pinned to that value (hard set value, no uncertainty), for a sensitivity sweep.

`apply_crab_fraction = 0` (shore fits, or `use_crab_fraction = FALSE`) pins `f = 1` in Stan (behavior-neutral).

## Data hook (I/E classification -> f)

`fetch_ie_data()` now emits the classification counts as `attr(ie_data, "crab_fraction_counts")`, read from two OPTIONAL columns on the boat-I/E (WBL) rows of `ingress_egress.xlsx`: `boats_crabbing` and `boats_total` (names configurable via `ie_crab_col` / `ie_total_col`). They do not exist yet (the WBL egress-classification pilot is in progress), so the reader reports them absent and f uses the set-value fallback. When the pilot lands, add those two columns to the WBL rows and f updates automatically; no code change.

## Files changed

- `02_stan_models/crab_bss_pooled.stan`, `crab_bss_gear_resolved.stan`: `f_crab` (parameter or pinned) with a Beta prior + optional Binomial classification likelihood; applied to boat `E` and `C` in generated quantities; `f_crab_out` reported.
- `03_R_functions/crab_fraction.R` (new): `crab_fraction_stan_data()`, the shared helper (set value, prior, fallback, fixed-pin), used by both preps.
- `03_R_functions/bss_day_length.R`: `fetch_ie_data()` emits the classification counts attribute.
- `03_R_functions/prep_bss_crab_pooled.R`, `prep_bss_crab_gear.R`: call the helper (boat fits) and pass the 7 `crab_fraction_*` keys to `stan_data`.
- `run_config.R`: `use_crab_fraction` (TRUE), `crab_fraction_set` (0.3), `crab_fraction_prior_kappa` (20), `crab_fraction_fixed` (NA), `crab_fraction_min_obs` (20), `ie_crab_col` / `ie_total_col`.

## Known limitations / backlog

- **PE not yet f-adjusted (Phase 2b).** `run_pe_pooled` / `run_pe_gear` still produce all-boat catch, so with f < 1 the boat PE will read HIGH relative to the BSS (the reverse of today). Apply the same f to the boat PE for a fair cross-check. Small follow-up.
- **Single scalar f (Phase 3).** Crab share very plausibly varies by month/day-type (lower in the summer tuna peak). A single annual f over-states summer crab effort. Phase 3 makes f time-varying once the pilot supports per-stratum counts; the pilot should classify by month and day-type.
- **The set value is load-bearing.** Until the pilot lands, the boat number scales linearly with `crab_fraction_set`. Treat the f = 0.3 run as illustrative, and run a sensitivity sweep over `crab_fraction_fixed` in {0.2, 0.3, 0.5, 1.0}.

## Validation gate

Run each model with `use_crab_fraction = TRUE`: boat `C_sum` and `E_sum` drop to ~0.3x their FALSE-baseline values, model CPUE (crab/deployment) is UNCHANGED (f cancels), convergence is unchanged (f is decoupled; identical effort/CPUE posterior), and `f_crab_out` ~ Beta(6, 14) (median ~0.29) in the prior-only state. Then a `crab_fraction_fixed` sweep to show the linear f-sensitivity of the boat total.
