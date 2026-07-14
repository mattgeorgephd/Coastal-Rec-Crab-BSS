# 05_output

Run results. Each driver writes a complete set of estimates, diagnostics, and plots here. Nothing in this folder is hand-authored; it is all regenerated from `04_input_files/` plus the chosen model and parameters, so any run can be reproduced by re-knitting the matching driver.

For the meaning of the populations, sub-seasons, and the PE-vs-BSS gate referenced throughout, see the [root README](../README.md).

## Folder structure

```text
05_output/
  <run_date>/                 # YYYYMMDD, set by run_date <- format(Sys.Date(), "%Y%m%d")
    <model>/                  # one subfolder per driver run that day
      <output files>
```

`output_dir <- here("05_output", run_date, "<model>")` is set once in each driver's setup chunk; every write lands under it. A folder is created fresh on the date you knit, so older dated folders are immutable records of past runs. Because the pipeline has grown over time, **older runs legitimately contain fewer or differently named files than recent ones**; that is history, not breakage.

The `<model>` subfolder name is whatever the driver set, and it has changed across versions. You will see, for example, `pooled-CPUE-model`, `pooled-CPUE`, `gear-type-CPUE-model`, and `pooled-CPUE-covariates`. The `*-covariates` subfolders are produced by the weather-tide module in `06_diagnostics/`, not by a production driver.

## File naming convention

Most per-population files follow:

```text
<metric>_<population>_<species>_<fate>.{csv,png}
```

- `population` is one of `shore_ring_net_only`, `shore_all_gear`, `private_boat_all_gear`.
- `species` is `Dungeness`; `fate` is `Kept`.

So `bss_daily_catch_private_boat_all_gear_Dungeness_Kept.csv` is the daily catch series for the private-boat all-gear component. Port-level and cross-population files (totals, monthly, comparisons, metadata) drop the population tag.

The commercial/charter component has **no** per-population BSS file; it enters the port total through census expansion of `wes_commercial_tally.csv`.

## File catalog (representative production run)

The groups below reflect a recent pooled production run. A gear-resolved run produces the same families plus the per-gear extras noted at the end.

| Group | Files (per population unless noted) | What it is |
|---|---|---|
| BSS daily series | `bss_daily_effort_*`, `bss_daily_cpue_*`, `bss_daily_catch_*` | Posterior daily effort, CPUE, and catch with intervals. |
| BSS summaries | `bss_summary_*`, `bss_full_summary_*`, `bss_draws_summed_*`, `bss_ar_path_*`, `bss_period_coverage_*`, `bss_L_effective_*` (shore only) | Season totals, full parameter summary, summed posterior draws, the AR(1) path, period coverage, and the estimated effective day length. |
| Port / monthly / mode | `port_total_Dungeness_Kept.csv`, `monthly_estimates.csv`, `monthly_estimates_by_mode.csv`, `catch_by_mode.csv`, `season_summary.csv`, `expansion_ratios.csv`, `effort_cpue_multipliers.csv`, `fit_data_summary.csv` | The headline harvest number and its monthly and by-mode breakdowns. |
| Gear | `catch_by_gear_type.csv`, `gear_proportions.csv` | Gear-type catch. In a pooled run these are derived after estimation from interview proportions; in a gear-resolved run they carry posterior uncertainty (see extras below). |
| PE vs BSS | `pe_port_summary.csv`, `pe_vs_bss_comparison.csv`, `monthly_pe_vs_bss.csv` | Point Estimator results and the side-by-side reconciliation used by the convergence gate. |
| I/E | `ie_analysis.csv`, `L_effective_ie_detail.csv` | Ingress/egress summary and the per-observation detail behind the `L_effective` regression. |
| Convergence and structure | `convergence_report.csv` (one file, all populations), `divergence_localization_*`, `sampler_diagnostics_*`, `structural_params_*`, `prior_vs_posterior_*` | R-hat / divergence / treedepth reporting, where divergences land, sampler behavior, structural parameters, and prior-vs-posterior overlap. |
| Effort overdispersion | `effort_overdispersion_byobs_*`, `effort_overdispersion_decomp_*` | The sparse per-observation effort overdispersion, by observation and decomposed. |
| Cross-validation (LOO) | `loo_summary_*`, `loo_pointwise_catch_*`, `loo_pointwise_gear_*` (shore), `loo_pointwise_trailer_*` (boat) | PSIS-LOO summaries and pointwise contributions by likelihood component. |
| Posterior predictive checks | `ppc_byobs_*`, `ppc_calibration_*`, `ppc_pit_*.png` | Observation-level PPC, calibration tables, and PIT histograms. |
| Plots | `plot_bss_effort_*`, `plot_bss_catch_*`, `plot_bss_posteriors_*`, `plot_cpue_timeseries.png`, `plot_shore_effort_timeseries.png`, `plot_boat_effort_timeseries.png`, `plot_effort_by_month.png`, `plot_monthly_catch.png`, `plot_monthly_catch_by_mode.png`, `plot_monthly_catch_by_mode_facet.png`, `plot_catch_by_gear_type.png`, `plot_L_effective_regression.png`, `plot_day_length_comparison.png` | Figure versions of the series and summaries. |
| Run metadata | `run_parameters.txt`, `session_info.txt` | The exact parameters and the R/package session for that run. |

**Gear-resolved extras.** A `gear-type-CPUE-model` run adds per-gear detail such as `catch_by_gear_type_detail.csv` and `sensitivity_incomplete_by_gear.csv`, and its gear-type catch comes from the per-gear CPUE process rather than post-hoc proportions.

**Weather-tide (covariates) runs.** A `*-covariates` subfolder looks different by design: paired `*_baseline` / `*_covariates` files, `loo_comparison_*` and `pareto_k_*`, GAM smooths (`gam_*_smooths.csv`), covariate effect and inclusion tables, `ci_width_comparison.csv`, `daily_covariates.csv`, and data-source logs (`asos_station_log.csv`, `ndbc_station_log.csv`, `tide_station_log.csv`). These come from `06_diagnostics/`; see [06_diagnostics/README.md](../06_diagnostics/README.md) and the conclusion in `07_documentation/WEATHER_COVARIATE_ANALYSIS.md`.

## Why outputs are committed

These results are checked into the repository so a given seasonal estimate and its diagnostics are preserved exactly as produced, without needing to re-run. The trade-off is repository size; this folder is the bulk of it. Per-run R workspaces (`*.RData`) and compiled model caches (`*.rds`) are **not** committed; they are git-ignored because they are machine-local and regenerate on demand (see `.gitignore`). If repository size becomes a problem, the alternative is to stop committing dated runs and keep only the inputs plus a tagged driver, accepting that past estimates would then have to be reproduced rather than read.
