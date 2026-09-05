# Incorporate OSP boat counts: second effort stream, crabbing fraction f, OSP-informs-tau

**Branch:** `OSP-boat-count-incorporation` -> `main`
**Convention:** no em dashes.

> ## STATUS 2026-09-08: this PR description covers only the FIRST wave of the branch
>
> It was written 2026-08-05, when the branch held the OSP integration alone, and its
> "Confirmed harvest numbers" (67,312 pooled) are three adoptions out of date. The branch
> has since also adopted, each validated by its own gated run: the **shared boat turnover**
> `tau_bar` (2026-09-01; the OSP series identifies one estimated turnover, roughly doubling
> the boat component relative to the pre-OSP baseline), the **weekly shore all-gear AR**
> (2026-09-07; settled by an escalation ladder that showed daily was overfitted, p_loo
> 35.2% -> 9.5% of n_obs, 41 -> 0 bad Pareto k), and the **zero-inflated shore catch
> likelihood** (2026-09-07), plus the diagnostics and validation machinery in
> `06_diagnostics/` (452-assertion harness, dated batch runners, per-rung ladder adequacy,
> posterior draw persistence).
>
> **The branch as it stands is summarized in two places, which a reviewer should read
> instead of this file:** `07_documentation/development_notes/CHANGE_REGISTER.md` (every
> change, its status, evidence, and effect on the number) and
> `07_documentation/development_notes/PIPELINE_STATUS.md` (the narrative, Sections 1b-1m,
> with the authoritative run in the box at the top: `20260904/pooled-CPUE-AD-A1-adopted`,
> port total 72,027 [53,018, 101,364]). The text below is kept as the record of the
> original OSP scope and its validation.

## Summary

Adds the OSP Westport Boat Launch (WBL) daily boat-count series to the private-boat effort model, converts all-boat effort to crab effort with a crabbing fraction `f`, and lets the dense OSP series identify the within-day boat turnover (`tau_boat`), resolving the long-standing GR-12 dependence on two winter I/E days. Two interview data fixes and a reader-deduplication refactor ride along. The production config is now the committed default, and both models have been re-run to confirm the numbers.

This is a boat-effort change; shore and commercial/charter are untouched.

## What changed

- **OSP second boat-effort stream (Phase 1, both models).** `04_input_files/WBL_boat_counts.xlsx`, read by `03_R_functions/fetch_osp_boat_counts.R`, enters both Stan models as a second observation on the latent boat effort `lambda_E`, scaled by the OSP within-day turnover `kappa_OSP`. Calibration/coverage diagnostics in `03_R_functions/diagnose_osp_trailer_overlap.R`. Toggle `use_osp_boat_counts`. OSP covers 148 days (mid-March to mid-October); outside that window the model runs trailer-only, so the pre-OSP behavior is the emergent fallback.
- **Crabbing fraction f (Phase 2/2b/3, boat only).** `03_R_functions/crab_fraction.R`. Trailer and OSP counts are all boats, so `f` converts all-boat effort to crab effort in the boat generated quantities only, decoupled from sampling, so it never affects convergence and leaves CPUE invariant. Beta prior centered on `crab_fraction_set` (production 0.3), updatable per-stratum by the WBL egress classification once it lands. The Point Estimator is f-adjusted too.
- **OSP-informs-tau (`osp_scale_is_tau`, Phase 3), production ON.** The dense OSP series identifies `tau_boat` directly. The crab-creel trailer count was confirmed to be an instantaneous snapshot, so the OSP/trailer within-day turnover (posterior `kappa_OSP` = 3.15, 95% [2.50, 3.91]) is real and the old `tau_boat` ~1.2 was a roughly 2x under-count; turning this on corrects that.
- **Two interview filters (both readers).** Non-crabbing interviews (`number_of_gear == 0`, any trip-completion status) and gear-tampered interviews (a new `gear_tampered == 1` column) are dropped. `gear_tampered` is all blank today, so it removes nothing yet.
- **Reader de-duplication.** `fetch_crab_data.R` and `fetch_crab_data_v2.R` were near-identical (one differing private-boat filter). They are now one shared `fetch_crab_data(params)` called by both drivers, with the differing filter behind `params$boat_require_gear_time` (TRUE pooled, FALSE gear). `fetch_crab_data_v2.R` is deleted.
- **`renv.lock`.** A top-level lockfile pins R 4.2.2 and the ~99 CRAN packages used to produce the estimates (see Reproducibility below).

## Production configuration (committed default)

`use_osp_boat_counts = TRUE`, `use_crab_fraction = TRUE` (`crab_fraction_set = 0.3`), `osp_scale_is_tau = TRUE`.

Note this default diverges from the frozen Method v1.0 (f = 1, no OSP); it is the intended production change, and the method-doc reference numbers have been refreshed to it.

## Confirmed harvest numbers (2024-25, Dungeness kept)

| Model | Private boat (BSS median) | Port total (BSS median) | PE port | Confirmation run |
|---|---|---|---|---|
| Pooled | 27,684 | 67,312 (95% CI 50,601 to 93,461) | 44,810 | `05_output/20260804/pooled-CPUE-boat-count-validation-run` |
| Gear-resolved | 27,524 | 66,461 (95% CI 50,259 to 91,309) | 44,810 | `05_output/20260805/gear-type-CPUE-model-boat-count-validation-run` |

The two independent pipelines reconcile to 0.6% on the boat and 1.3% on the port, the strongest internal cross-check available.

## Validation evidence

- **14-run validation batch** (`05_output/osp_validation_summary.csv`; review in `07_documentation/development_notes/osp-validation-review-2026-07-31.md`). All runs converged. Key checks: the boat catch is exactly linear in `f` (0.2/0.5/1.0 give 8.7k/21.6k/43.2k, CPUE invariant), `f = 1` reproduces the pre-f baseline to 0.1%, shore is byte-identical across every toggle (no leakage), and both models reconcile at every step.
- **Two production confirmation runs** (the table above). Both converge, and each reproduces its pre-merge validation Step 6 fit-for-fit (interview counts and every component `C_sum` byte-identical), which certifies that the reader merge and the interview filters are behavior-neutral on the current data. The two `boat_require_gear_time` branches (TRUE pooled, FALSE gear) are each exercised by one of the runs.

## Known limitations and next steps (not blockers)

- **The crabbing fraction `f = 0.3` is a placeholder,** the largest single lever on the boat harvest, pending the WBL egress crab-vs-total classification pilot. The boat scales linearly in `f`, so the Step 4 sweep in the review is the current sensitivity band. Time-varying `f` by month is coded and validated as plumbing; it binds once the pilot classifies enough boats per stratum.
- **The `osp_scale_is_tau` correction assumes a representative trailer snapshot** (not a fixed daily peak). A peak-timed snapshot would make the ~2x turnover correction a floor rather than a point estimate. Confirming the snapshot-timing protocol tightens the multiplier; it does not change its direction. See PIPELINE_STATUS Section 6.

## Notes for reviewers

- **Diff size.** The code and documentation change is about 50 files; the branch-vs-main diff is much larger because committed `05_output/` run artifacts (the validation batch and the two confirmation runs) are included per the repository's standing commit-outputs policy. Review the non-`05_output` files for the substance.
- **`renv.lock` provenance.** It was reconstructed from the confirmation run's `session_info.txt` (the exact environment used), not from `renv::snapshot()`. It therefore carries no package hashes and covers the two production models' observed dependencies; the experimental weather module additionally needs `geosphere`. Run `renv::restore()` to validate, or `renv::init()` + `renv::snapshot()` locally to regenerate the fully canonical lock.

## Reproducibility

R 4.2.2, rstan 2.32.7 (fixed `bss_seed = 20260619`). Restore the environment with `renv::restore()` against the committed `renv.lock`, place inputs in `04_input_files/`, set the model in `run_config.R`, and run `Rscript run_estimation.R`.
