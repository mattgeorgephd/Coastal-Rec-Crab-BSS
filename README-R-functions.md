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

The 2026-07-11 refactor pulled the per-driver data-prep, PE, BSS-prep, and utility functions out of the two `.Rmd` drivers and into this folder (previously they were inline in each `.Rmd`). The files now fall into five groups.

### Shared driver modules (pooled + gear-resolved)

One implementation each, called by both production drivers, so the two tracks cannot drift.

| File | Public function(s) | Role |
|---|---|---|
| `bss_convergence_gate.R` | `bss_compute_gate`, `bss_use_pe_for` | Scale-aware convergence gate (B1.8): the per-fit PE-vs-BSS pass/fail decision and its report row. Thresholds are arguments, so each driver passes its own. |
| `bss_ar_resolution.R` | `bss_select_ar_resolution`, `bss_ar_ladder` | Adaptive AR(1) temporal-resolution selector (daily / weekly / biweekly / monthly), with the per-population cap and the `ar_force` override. `bss_ar_ladder` (2026-08-25) returns the opt-in escalation ladder: with `ar_escalate = TRUE` a fit starts at the finest rung and coarsens only when it fails the convergence gate, so each component reports at the finest resolution its own sampler behaviour supports. Ships off; every rung is a real multi-hour fit. |
| `bss_cpue_diagnostics.R` | `write_cpue_diagnostics` (+ `bss_cpue_estimator_triad`, `bss_saturation_exponent`, `bss_effort_linearity`, `bss_assert_effort_units`) | Per-fit CPUE effort-unit diagnostics (estimator triad, saturation exponent, linearity slope) and the effort-unit assertion. |
| `bss_stan_fit.R` | `bss_stan_fit`, `bss_assert_stan_data`, `bss_assert_fit_usable`, `bss_stan_data_names` | Guarded wrapper around `rstan::stan()`, called by both drivers (2026-08-26). On the happy path it is `rstan::stan()` with the same arguments, seed and RNG stream, so fits are unchanged. Before sampling it parses the `data` block of the `.stan` file in use and refuses to run when a declared variable is missing from the data list, non-numeric, or non-finite; after sampling it rejects an empty `stanfit` and quotes what Stan said. It exists because rstan reports a data-initialization failure on the MESSAGE stream and then returns an empty `stanfit` rather than raising, which on 2026-08-25 cost six model-runs and surfaced only as `dim(X) must have a positive length` three layers downstream. `log_file` keeps the per-fit Stan console in `output_dir` so a `results='hide'` chunk stays diagnosable. |
| `bss_effort_spec.R` | `bss_effort_spec`, `bss_effort_h_candidates` | The single effort-unit specification (crabber-hours / gear-hours / gear-deployments) that both the BSS prep and the PE read, so effort and CPUE always share a unit. As of 2026-08-25 it also names the matching I/E OBSERVATION column (`ie_obs_col`): under gear-deployments the predicted `lambda_E * L` is crabber TRIPS, so the observation is the arrival count, not crabber-hours. Route any new consumer of effort through this function rather than re-deriving the formula. |
| `bss_trailer_expansion.R` | `bss_trailer_par`, `bss_extract_pars`, `bss_trailer_multiplier` | Boat trailer-expansion adapter that abstracts the `R_T` (legacy) vs `R_G_boat` (current) split so downstream expansion code is unchanged. |
| `bss_day_length.R` | `fetch_ie_data`, `estimate_L_effective`, `bss_day_length_civil`, `bss_assign_day_length` | I/E ingest and the effective-day-length (`L_effective`) model with the civil-twilight fallback ladder. |
| `bss_timers.R` | `timer_start`, `timer_stop`, `bss_timer_log` | Section timers for the end-of-run timing summary. State lives in a module-local environment (reset on source), so no driver global is needed. |
| `build_subseasons.R` | `build_subseasons` | Derives the pot-closure sub-seasons (the non-pot vs all-gear split) that both drivers fit separately, keeping the internal `ring_net_only` key for output-filename continuity. |
| `prep_days_crab.R` | `prep_days_crab` | Builds the per-day calendar (indices, day type, day-type integer, effective day length). Takes `params` and derives day-typing inputs from it. |
| `prep_population_summary.R` | `prep_population_summary` | Filters the data bundle to one population x sub-season and builds its effort, interview, and catch frames plus the gear/crabber ratios. |
| `estimate_comm_charter.R` | `estimate_comm_charter` | Day-type-stratified census expansion of the commercial/charter vessel tally (with optional red-rock, guarded by `params$estimate_red_rock`). |
| `pe_monthly_effort_share.R` | `pe_monthly_effort_share` | Consolidated PE monthly effort-share helper: the shared share math for spreading a PE-fallback component's catch and effort across months, with each driver keeping its own draw accumulation and uncertainty handling. Its shore branch was still on the pre-v7.7 crabber-hours formula (and therefore on the seasonal day-length curve) until 2026-08-25; it now reads `bss_effort_spec`. |
| `bss_opener_covariates.R` | `opener_flag_series`, `opener_select`, `opener_design_matrix` | Other-fishery opener EFFORT covariates (2026-08-25), generalizing the retired razor-dig-only `B3` into the `K_open` / `X_open` / `B_open` block in both Stan models. `opener_select` is the automatic screen: adjusted p-values from the spillover diagnostic, with a multiplicity correction across the eight effort tests (BH by default). Treat it as a screen, not a verdict; the arbiter is the effort-stream `elpd_loo`. |
| `diagnose_incomplete_trips.R` | `diagnose_incomplete_trips`, `incomplete_trip_arm_frames` | The four-arm incomplete-trip treatment diagnostic (exclude / gear_only / impute_mean_cpue / keep) run through each track's own point estimator. DIAGNOSTIC ONLY; production remains `exclude`. |
| `classify_day_type.R` | `classify_day_type` | Standalone day-type classifier for any date, used by diagnostic plots outside the estimation window. |
| `read_crabbing_holidays.R` | `read_crabbing_holidays` | Single-source reader for the crabbing-holiday calendar (`04_input_files/crabbing_holidays.xlsx`), filtered to `season_filter`. Replaces the hardcoded `crabbing_holiday_dates` vector formerly in `run_config.R`; stops loudly if the file, columns, or season are missing. Each driver calls it once and sets `params$crabbing_holiday_dates`, so downstream consumers are unchanged. |

### Pooled-CPUE driver functions

Called only by `BSS-GH-pooled-CPUE-model.Rmd`; named to avoid a collision with the gear-resolved equivalents.

| File | Function | Role |
|---|---|---|
| `fetch_crab_data.R` | `fetch_crab_data` | Shared reader for BOTH tracks (merged 2026-08-01 from the former fetch_crab_data / _v2 pair; the one differing private-boat filter is now the `boat_require_gear_time` param, TRUE for pooled, FALSE for gear-resolved). Reads and assembles inputs and classifies interviews by population; drops non-crabbing (`number_of_gear == 0`, any trip status) and gear-tampered (`gear_tampered == 1`) interviews. Per-gear CPUE classification is downstream in `prep_bss_crab_gear.R`, not here. |
| `run_pe_pooled.R` | `run_pe_pooled` | Pooled Point Estimator (stratified effort and catch). The shore branch reads its effort unit and CPUE denominator from `bss_effort_spec`, so the shore PE matches the shore BSS (2026-07-11 fix). |
| `prep_bss_crab_pooled.R` | `prep_bss_crab_pooled` | Build the Stan data list for `crab_bss_pooled.stan`. |

### Gear-resolved driver functions

Called only by `BSS-GH-gear-type-CPUE-model.Rmd`.

| File | Function | Role |
|---|---|---|
| `run_pe_gear.R` | `run_pe_gear` | Gear-resolved Point Estimator, with the P0/P1/P2 fixes (explicit population argument, `bss_effort_spec` effort unit, ratio-of-sums stratum CPUE, and a scale-consistency assertion). |
| `prep_bss_crab_gear.R` | `prep_bss_crab_gear` | Build the Stan data list for `crab_bss_gear_resolved.stan`. |

### OSP boat counts and crabbing fraction

Feature functions from the OSP boat-count incorporation branch, used by both production drivers. The first reads the OSP daily boat-total as a second boat-effort stream; the second calibrates it against the trailer counts; the third builds the crabbing fraction `f`.

| File | Public function(s) | Role |
|---|---|---|
| `fetch_osp_boat_counts.R` | `fetch_osp_boat_counts` | Reads `04_input_files/WBL_boat_counts.xlsx` as the OSP daily boat-total stream (the second boat-effort observation on `lambda_E`), gated by `run_config$use_osp_boat_counts`. Also emits the optional per-day crab-ONLY counts as `attr(., "osp_crab_rows")` (2026-08-25), which feed the hard lower bound on `f`. |
| `diagnose_osp_trailer_overlap.R` | `diagnose_osp_trailer_overlap` | Calibrates the OSP boat-total against the trailer counts on overlapping days and writes the OSP-vs-trailer calibration, pairs, plot, and coverage audit. |
| `crab_fraction.R` | `crab_fraction_strata_labels`, `crab_fraction_stan_data`, `crab_fraction_point_day` | Builds the crabbing fraction `f` (per-stratum / per-day) that scales boat effort and catch in the BSS boat generated quantities and the PE. As of 2026-08-25 `f = f_lower + (1 - f_lower) * theta`, where `f_lower` is the OSP-observed crab-ONLY share (a hard lower bound: OSP labels combo trips by the non-crab fishery) and `theta` is the combo share among the not-crab-labelled boats, identified only by the egress classification. With no OSP counts `f_lower = 0` and the behaviour is bit-identical to the previous Beta-prior `f`. |

### Crab-specific BSS diagnostics

Per-fit and per-run diagnostic writers, mostly `tryCatch`-wrapped so one fit cannot abort a run, plus two report-embedded pooled diagnostics (fishery-opener spillover and shore I/E representativeness) and the opener-calendar prep they use. The per-fit writers produce the convergence, divergence, over-dispersion, PPC, and extended-output files catalogued in [05_output/README.md](05_output/README.md).

| File | Public function(s) | Role |
|---|---|---|
| `model_diagnostics.R` | `bss_structural_summary`, `bss_decoupled_reasons`, `bss_divergence_localization`, `bss_ppc_calibration` | Structural-parameter summary, divergence localization, and posterior-predictive calibration per fit. **As of 2026-09-01 the PPC coverage is read off the RANDOMIZED PIT** rather than a quantile interval of simulated draws, which over-covers small counts by construction and disagreed with `ppc_byobs_*.csv` by up to 0.154 on the private-boat trailer stream. As of 2026-08-27 `structural_params_*.csv` carries a `decoupled` flag and its reason: several parameters here are declared with a proper prior unconditionally, so an unused one reports its PRIOR in the same columns as an estimate. The 2026-08-26 ladder's worked example is `kappa_OSP`, inert under the production `osp_scale_is_tau = TRUE` and reporting `lognormal(log 3, 0.3)` exactly. As of 2026-08-30 the same file also carries `fit_method` and an `estimate` column, which is `median` with the decoupled rows blanked, so a prior cannot be read out of an estimate column by a reader who did not check the flag. |
| `bss_model_adequacy.R` | `bss_model_adequacy`, `write_model_adequacy`, `annotate_model_adequacy_run` | Model-adequacy diagnostics reported BESIDE the convergence gate, never gating it (improvement-plan 5.5): `p_loo` as a fraction of `n_obs`, the Pareto k > 0.7 count, the worst PIT bias, and the smallest `n_eff` among the observation-model dispersion parameters the gate never inspects (decoupled ones excluded). Writes `model_adequacy.csv`, quoting the convergence gate's own `method_selected` rather than a second wording for it. **As of 2026-08-31 it also carries `cov50_worst_dev`, `pit_sd_worst_dev` and `flag_miscalibrated`,** because the 2026-08-30 Stage 5 batch showed `pit_mean` alone ranks configurations BACKWARDS: a latent process with one state per observation interpolates the data and piles every PIT value at 0.5, so the worst-calibrated cell in that batch (98% of observations inside a nominal 50% interval) had the best PIT mean. `disp_scale_min` reports the smallest dispersion SCALE, unflagged, because `disp_neff_min` sees sampling efficiency and not collapse. `annotate_model_adequacy_run()` rebuilds the same table for an ARCHIVED run folder from its committed CSVs into `model_adequacy_reconstructed.csv`; both paths share one body so the two are directly comparable. Exists because six configurations passed every gate criterion while spanning 44% on the boat component. |
| `annotate_decoupled_run.R` | `annotate_decoupled_run` | Retro-fits the `decoupled` flag onto a run folder that predates it, reconstructing each rule from the fit label, `fit_data_summary.csv` and `run_parameters.txt` into `decoupled_audit.csv`. Committed outputs are never rewritten; every row is marked `source = "reconstructed"`, and a rule that cannot be reconstructed reports `unknown` rather than a guess. |
| `bss_sampler_override.R` | `bss_apply_sampler_override` | The one sanctioned way for `run_config` to reach a per-fit sampler setting. Both drivers merge `params_model` ON TOP of `run_config`, so `params_model` wins every sampler key and a `run_config` delta on `bss_iter_default` silently does nothing. This named list is applied AFTER the merge, is restricted to sampler keys by pattern, prints every change, and ERRORS on a non-sampler key rather than dropping it. Experiments only; production value is `NULL`. |
| `diagnose_effort_overdispersion.R` | `write_effort_overdispersion_diag` | Law-of-total-variance decomposition of each effort-count predictive variance (Poisson floor / NB over-dispersion / latent process). |
| `divergence_diagnostic.R` | `diagnose_divergences` | Interactive, console divergence funnel-neck ranking (run by hand post-fit). |
| `save_run_diagnostics.R` | `write_fit_extended_diagnostics`, `write_loo_diagnostics`, `write_run_level_diagnostics`, `write_pe_empty_stratum_report` | The extended per-fit output series (O1-O13) and the run-level PE-vs-BSS and gear-proportion summaries. 2026-08-27: `write_pe_empty_stratum_report()` writes `pe_empty_effort_strata.csv` from both drivers, because the counts were previously only `cat()`-ed and the pooled PE chunk is `results='hide'`. `.srd_monthly_share()` no longer applies the shore day-length weighting when `L` is a turnover; that was the second, unfixed copy of improvement 1's defect (B). `fit_data_summary.csv` gained I/E and effort-unit provenance. |
| `prep_fishery_events.R` | `prep_fishery_events` | Reads the consolidated MA2 finfish and razor-clam-dig opener calendar (`fishery_opener_dates.xlsx`, now the only source) into per-date OPEN/CLOSED flags for the spillover diagnostic and the opener effort covariates. Stops if the workbook is absent; the former per-fishery fallback workbooks are retired. |
| `diagnose_fishery_spillover.R` | `diagnose_fishery_spillover` | Tests whether crab effort and CPUE differ on other-fishery opener days (pooled report Section 3.5), raw and adjusted for day type and month. Its adjusted table is also the input to the opener-covariate screen. CPUE is measured per gear-deployment as of 2026-08-25 (it was per crabber-hour, the denominator the pipeline's own linearity test rejects). |
| `diagnose_ie_representativeness.R` | `diagnose_ie_representativeness` | Shore I/E representativeness diagnostic: whether the I/E observation days are unrepresentative (measured on peak-effort days) or just sparse (pooled report). |

## Configuration and paths

User-selectable toggles live in [`run_config.R`](run_config.R), the single source of truth; these functions read them through the `params` list the driver passes. Path-writing helpers (the diagnostic writers above) resolve outputs through `here("05_output", ...)`. If a stage folder is renamed, check those writers plus the drivers (see "How paths work" in the root README).
