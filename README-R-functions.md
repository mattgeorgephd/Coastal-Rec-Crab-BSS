# 03_R_functions

Helper functions for the estimation pipeline. There is no entry point here; the drivers in `01_BSS_models/` and `06_diagnostics/` load this folder wholesale in their setup chunk:

```r
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)
```

Because the whole folder is sourced, every `.R` file here is expected to define functions only and to have no side effects at source time (no reads, writes, or plotting at the top level). Sourcing order is not guaranteed, so a function in one file may not assume another file has already run anything beyond defining its functions.

For the project overview and the PE-vs-BSS split these functions implement, see the [root README](../README.md).

## Function groups

The files fall into eight groups (the shared `bss_*` driver modules were centralized in the v6.5 to v7.5 refactors so the pooled and gear-resolved tracks share one implementation). Several in the data-fetch and PE groups carry over from the WDFW freshwater-creel codebase this project was forked from; they are retained because the Point Estimator (PE) path and the optional database export reuse that machinery.

| Group | Files | Role |
|---|---|---|
| Database / ETL (creel lineage) | `establish_db_con`, `fetch_db_table`, `fetch_dwg`, `write_db_tables`, `confirm_db_upload`, `export_estimates`, `prep_export`, `JSON_conversion`, `generate_analysis_lut`, `map_data_grade`, `transform_estimates` | Connect to the warehouse, pull raw tables, and push or serialize finished estimates. Used by the PE/export path; not required for a local BSS-only run. |
| PE data prep | `prep_days`, `prep_dwg_census_expan`, `prep_dwg_effort_census`, `prep_dwg_effort_index`, `prep_dwg_interview_angler_types`, `prep_dwg_interview_catch`, `prep_dwg_interview_fishing_time` | Reshape raw effort-count and interview pulls into the per-day, per-stratum frames the PE consumes. |
| PE input assembly | `prep_inputs_pe_df`, `prep_inputs_pe_ang_hrs`, `prep_inputs_pe_days_total`, `prep_inputs_pe_int_ang_per_object`, `prep_inputs_pe_paired_census_index_counts`, `prep_inputs_pe_daily_cpue_catch_est` | Build the specific input objects the PE estimators need (angler-hours, day totals, paired census/index counts, daily CPUE/catch). |
| PE estimation | `est_pe_effort`, `est_pe_catch`, `process_estimates_pe` | Compute the Point Estimator effort and catch, then collate into the PE result tables used for the convergence-gate fallback and the PE-vs-BSS comparison. |
| BSS | `prep_inputs_bss`, `fit_bss`, `get_bss_overview`, `get_bss_effort_daily`, `get_bss_cpue_daily`, `get_bss_catch_daily`, `process_estimates_bss` | Assemble the Stan data list, extract the fitted posterior into daily effort/CPUE/catch series, and collate BSS results. (See the `fit_bss` flag below.) |
| Plotting | `plot_census_index_counts`, `plot_est_pe_effort`, `plot_est_pe_catch`, `plot_inputs_pe_census_vs_index`, `plot_inputs_pe_cpue_period`, `plot_inputs_pe_index_effort_counts` | Input-diagnostic and PE-result plots written to the run's output folder. |
| Shared driver modules (pooled + gear-resolved) | `bss_convergence_gate`, `bss_ar_resolution`, `bss_cpue_diagnostics`, `bss_trailer_expansion`, `bss_day_length` | One implementation each of the scale-aware convergence gate and PE-vs-BSS selector (`bss_compute_gate` / `bss_use_pe_for`), the adaptive AR-resolution selector (`bss_select_ar_resolution`), the CPUE effort-unit diagnostics (`write_cpue_diagnostics`), the trailer-expansion adapter, and the I/E effective-day-length model. Both production drivers call these so the two tracks cannot drift; as of pooled v7.5 the pooled driver uses the shared gate, AR selector, and CPUE diagnostics rather than inline copies (POOL-5 / POOL-6). |
| Crab-specific BSS diagnostics | `model_diagnostics`, `diagnose_effort_overdispersion`, `divergence_diagnostic`, `save_run_diagnostics` | Convergence reporting, effort-overdispersion decomposition, divergence localization, and the per-run diagnostic bundle. These write the convergence/divergence/overdispersion files catalogued in [05_output/README.md](../05_output/README.md). |

## Path handling

The four diagnostic files in the last group write their outputs through `here("05_output", ...)` and were updated when the stage folders were renumbered. The rest of the functions either take paths as arguments from the driver or operate purely in memory, so they did not need path edits. If a stage folder is renamed again, check these four files plus the drivers (see "How paths work" in the root README).

## Known issue: `fit_bss.R` is dead code

> **Flag (not a path issue).** `fit_bss.R` has **zero call sites** anywhere in the repo; the production drivers fit the model with a direct `rstan::stan(file = here("02_stan_models", params$bss_model_file), ...)` call instead. On top of being unused, its default argument still points at the old, un-prefixed folder **and** a retired model that no longer exists:
> ```r
> model_file_name = here::here("stan_models/BSS_creel_model_02_2024-04-03.stan")
> ```
> Neither `stan_models/` nor that `.stan` file exists in the reorganized tree. This path was deliberately left unedited: "fixing" it would make a non-functional helper look usable. Recommend either deleting `fit_bss.R` or repairing it to call the current `02_stan_models/` models if you intend to revive a reusable fit wrapper. No impact on current runs.
