# 03_R_functions

Helper functions for the estimation pipeline. There is no entry point here; the drivers in `01_BSS_models/` and `06_diagnostics/` load this folder wholesale in their setup chunk:

```r
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)
```

Because the whole folder is sourced by BOTH drivers, every `.R` file here must define functions only, with no side effects at source time (no reads, writes, or plotting at the top level), and no assumption about sourcing order. Two consequences follow, and both are load-bearing:

1. **Pure functions.** A helper takes everything it needs as arguments (or derives it from a passed `params` list). It does not read driver globals such as `catch_groups` or `crabbing_holiday_dates`. Config is passed, not captured.
2. **No name collisions.** Since both drivers source every file, two files cannot define the same function name with different bodies. Functions that are genuinely shared keep one name and one definition; functions that differ between the pooled and gear-resolved tracks carry distinct names (the `_pooled` / `_gear` suffix), so both can coexist in the sourced environment and each driver calls its own.

For the project overview and the PE-vs-BSS split these functions implement, see the [root README](README.md).

## Function groups

The 2026-07-11 refactor pulled the per-driver data-prep, PE, BSS-prep, and utility functions out of the two `.Rmd` drivers and into this folder (previously they were inline in each `.Rmd`). The files now fall into four groups.

### Shared driver modules (pooled + gear-resolved)

One implementation each, called by both production drivers, so the two tracks cannot drift.

| File | Public function(s) | Role |
|---|---|---|
| `bss_convergence_gate.R` | `bss_compute_gate`, `bss_use_pe_for` | Scale-aware convergence gate (B1.8): the per-fit PE-vs-BSS pass/fail decision and its report row. Thresholds are arguments, so each driver passes its own. |
| `bss_ar_resolution.R` | `bss_select_ar_resolution` | Adaptive AR(1) temporal-resolution selector (daily / weekly / biweekly / monthly), with the per-population cap and the `ar_force` override. |
| `bss_cpue_diagnostics.R` | `write_cpue_diagnostics` (+ `bss_cpue_estimator_triad`, `bss_saturation_exponent`, `bss_effort_linearity`, `bss_assert_effort_units`) | Per-fit CPUE effort-unit diagnostics (estimator triad, saturation exponent, linearity slope) and the effort-unit assertion. |
| `bss_effort_spec.R` | `bss_effort_spec`, `bss_effort_h_candidates` | The single effort-unit specification (crabber-hours / gear-hours / gear-deployments) that both the BSS prep and the PE read, so effort and CPUE always share a unit. |
| `bss_trailer_expansion.R` | `bss_trailer_par`, `bss_extract_pars`, `bss_trailer_multiplier` | Boat trailer-expansion adapter that abstracts the `R_T` (legacy) vs `R_G_boat` (current) split so downstream expansion code is unchanged. |
| `bss_day_length.R` | `fetch_ie_data`, `estimate_L_effective`, `bss_day_length_civil`, `bss_assign_day_length` | I/E ingest and the effective-day-length (`L_effective`) model with the civil-twilight fallback ladder. |
| `bss_timers.R` | `timer_start`, `timer_stop`, `bss_timer_log` | Section timers for the end-of-run timing summary. State lives in a module-local environment (reset on source), so no driver global is needed. |
| `prep_days_crab.R` | `prep_days_crab` | Builds the per-day calendar (indices, day type, day-type integer, effective day length). Takes `params` and derives day-typing inputs from it. |
| `prep_population_summary.R` | `prep_population_summary` | Filters the data bundle to one population x sub-season and builds its effort, interview, and catch frames plus the gear/crabber ratios. |
| `estimate_comm_charter.R` | `estimate_comm_charter` | Day-type-stratified census expansion of the commercial/charter vessel tally (with optional red-rock, guarded by `params$estimate_red_rock`). |
| `classify_day_type.R` | `classify_day_type` | Standalone day-type classifier for any date, used by diagnostic plots outside the estimation window. |

### Pooled-CPUE driver functions

Called only by `BSS-GH-pooled-CPUE-model.Rmd`; named to avoid a collision with the gear-resolved equivalents.

| File | Function | Role |
|---|---|---|
| `fetch_crab_data.R` | `fetch_crab_data` | Read and assemble the pooled model's inputs and classify interviews by population. |
| `run_pe_pooled.R` | `run_pe_pooled` | Pooled Point Estimator (stratified effort and catch). The shore branch reads its effort unit and CPUE denominator from `bss_effort_spec`, so the shore PE matches the shore BSS (2026-07-11 fix). |
| `prep_bss_crab_pooled.R` | `prep_bss_crab_pooled` | Build the Stan data list for `crab_bss_pooled.stan`. |

### Gear-resolved driver functions

Called only by `BSS-GH-gear-type-CPUE-model.Rmd`.

| File | Function | Role |
|---|---|---|
| `fetch_crab_data_v2.R` | `fetch_crab_data_v2` | Read and assemble the gear-resolved model's inputs, with weighted gear-type classification of interviews. |
| `run_pe_gear.R` | `run_pe_gear` | Gear-resolved Point Estimator, with the P0/P1/P2 fixes (explicit population argument, `bss_effort_spec` effort unit, ratio-of-sums stratum CPUE, and a scale-consistency assertion). |
| `prep_bss_crab_gear.R` | `prep_bss_crab_gear` | Build the Stan data list for `crab_bss_gear_resolved.stan`. |

### Crab-specific BSS diagnostics

Per-fit and per-run diagnostic writers, all `tryCatch`-wrapped so one fit cannot abort a run. They write the convergence, divergence, over-dispersion, PPC, and extended-output files catalogued in [05_output/README.md](05_output/README.md).

| File | Public function(s) | Role |
|---|---|---|
| `model_diagnostics.R` | `bss_structural_summary`, `bss_divergence_localization`, `bss_ppc_calibration`, `write_bss_diagnostics` | Structural-parameter summary, divergence localization, and posterior-predictive calibration per fit. |
| `diagnose_effort_overdispersion.R` | `write_effort_overdispersion_diag` | Law-of-total-variance decomposition of each effort-count predictive variance (Poisson floor / NB over-dispersion / latent process). |
| `divergence_diagnostic.R` | `diagnose_divergences` | Interactive, console divergence funnel-neck ranking (run by hand post-fit). |
| `save_run_diagnostics.R` | `write_fit_extended_diagnostics`, `write_loo_diagnostics`, `write_run_level_diagnostics` | The extended per-fit output series (O1-O13) and the run-level PE-vs-BSS and gear-proportion summaries. |

## Configuration and paths

User-selectable toggles live in [`run_config.R`](run_config.R), the single source of truth; these functions read them through the `params` list the driver passes. Path-writing helpers (the diagnostic writers above) resolve outputs through `here("05_output", ...)`. If a stage folder is renamed, check those writers plus the drivers (see "How paths work" in the root README).
