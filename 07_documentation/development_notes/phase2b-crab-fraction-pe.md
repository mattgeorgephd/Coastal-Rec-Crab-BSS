# Phase 2b: crabbing fraction f applied to the Point Estimate

**Applies after:** the Phase 2 patch. **Scope:** the design-based PE (`run_pe_pooled.R`, `run_pe_gear.R`) plus a one-line driver hand-off. Small follow-up to Phase 2; no Stan change. **Status:** inherits the Phase 2 toggles (`use_crab_fraction`, `crab_fraction_set = 0.3`); FALSE reproduces the pre-Phase-2 boat PE exactly. R not executed here; validate by run.

## Why

Phase 2 applied `f` to the BSS boat only. With `f < 1` the BSS boat drops to ~f of its old value while the PE boat stayed all-boat, so the boat PE would read HIGH relative to the BSS (the reverse of today) and the `pe_vs_bss_comparison` cross-check would be an f artifact, not a model disagreement. Phase 2b puts the boat PE on the same crab-directed basis. It also keeps the PE FALLBACK consistent: when a boat BSS fit fails its convergence gate and the pipeline reports the PE instead, that PE is now f-adjusted too.

## What it does

Both `run_pe_pooled()` and `run_pe_gear()` multiply the boat `est_daily_effort` by a POINT fraction `f_pe = crab_fraction_point(is_boat, params)`. Because the PE catch is `effort x ratio-of-sums CPUE`, scaling the boat effort flows into the boat PE effort, its SE, and the boat PE catch, all by `f_pe`. Shore and commercial/charter are untouched (`f_pe = 1`).

`crab_fraction_point()` (new, in `03_R_functions/crab_fraction.R`) returns the CENTRAL value of the BSS `f`, i.e. the Beta posterior mean `(alpha0 + n_crab) / (alpha0 + beta0 + n_total)`. So:

- I/E classification thin/absent (now): `f_pe = crab_fraction_set` (0.3), matching the BSS prior-only median.
- classified boats accumulate: `f_pe` is the same shrunken posterior mean the BSS centers on.
- `crab_fraction_fixed` set: `f_pe` = that pinned value.

The I/E counts reach the PE via `params$crab_fraction_counts`, which each driver now lifts from `attr(ie_data, "crab_fraction_counts")` right after `fetch_ie_data()`, so the PE and the BSS read the SAME counts.

## The consistency check still holds (f cancels)

`run_pe_*` asserts the PE implied CPUE (catch / effort) matches the interview ratio-of-sums within 2x. Scaling both boat catch and boat effort by the same `f_pe` leaves implied CPUE unchanged, so the assertion is unaffected and cannot false-trip on f.

## Files changed

- `03_R_functions/crab_fraction.R`: add `crab_fraction_point()`.
- `03_R_functions/run_pe_pooled.R`, `run_pe_gear.R`: compute `f_pe` and apply it to the boat `est_daily_effort` (+ f in the log line).
- `01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd`, `BSS-GH-gear-type-CPUE-model.Rmd`: `params$crab_fraction_counts <- attr(ie_data, "crab_fraction_counts")` after `fetch_ie_data()`.

No `run_config` change (reuses the Phase 2 keys).

## Validation gate

Run each model with `use_crab_fraction = TRUE`: the boat PE effort and catch drop to ~f_pe (0.3) of their FALSE-baseline values, the PE implied-CPUE check still passes (f cancels), and `pe_vs_bss_comparison` for the boat is back to a like-for-like (both crab-directed) comparison. Shore and comm/charter PE unchanged.

## Still open

- **Phase 3:** time-varying f (month / day-type) once the pilot classifies by stratum; the PE would then take a per-stratum `f_pe` rather than one scalar.
- Any OTHER PE-derived boat breakdown computed outside `run_pe_*` (e.g. a separate monthly PE table, if present in a driver) should carry the same `f_pe`; the headline PE totals and the PE fallback are handled here.
