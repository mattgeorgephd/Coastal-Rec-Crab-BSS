# FWC-estimation-method

# Crab Creel Estimation - Westport Docks 2024-25 (v2 - Bug Fixes)

## Bugs Fixed from v1

| Bug | Symptom | Fix |
|-----|---------|-----|
| Interview CSV export: duplicate column name overwrite | `Crabbers-per-gear ratio: NaN`, PE crabber-hours = 0, IntA = 0 | Re-exported CSV mapping columns by letter position instead of name |
| Stan O dimension | `dims declared=(365,1,1); dims found=(365,1)` | Changed `O = matrix(...)` to `O = array(1.0, dim = c(D, S, G))` |
| No output folder | Results not saved | Added `output/YYYYMMDD/` with all plots, CSVs, and parameters |

## Files to Replace

1. **`interview_combined.csv`** → Replace in `input_files/` (number_of_gear column now populated)
2. **`crab_creel_westport_docks.R`** → Replace in project root (fixed O dimension + output folder)

The `effort_combined.csv` is unchanged.

## What the Script Now Produces (in output/YYYYMMDD/)

**Plots:**
- `plot_effort_timeseries.png` — Daily gear counts with day type colors
- `plot_cpue_timeseries.png` — Daily Dungeness CPUE trend
- `plot_effort_by_month.png` — Gear counts by month and day type (boxplots)
- `plot_pe_daily_effort.png` — PE estimated daily crabber-hours
- `plot_bss_daily_effort.png` — BSS estimated daily effort with 95% CI
- `plot_bss_daily_catch_dungeness.png` — BSS daily catch with 95% CI
- `plot_bss_posteriors.png` — Posterior distributions of season totals

**Data:**
- `results_summary.csv` — PE vs BSS comparison table
- `pe_effort_by_stratum.csv` — PE effort by month × day_type
- `pe_catch_dungeness_by_stratum.csv` — PE catch by stratum
- `pe_catch_redrock_by_stratum.csv` — PE catch by stratum
- `pe_daily_effort.csv` — PE daily effort estimates
- `bss_summary_dungeness.csv` — Stan summary statistics
- `bss_daily_effort.csv` — BSS daily effort (median + 95% CI)
- `bss_daily_catch_dungeness.csv` — BSS daily catch (median + 95% CI)
- `run_parameters.txt` — Analysis parameters used

## Expected Output After Fix

With the `number_of_gear` column now populated:
- `Crabbers-per-gear ratio` should be ~0.7–1.0 (a real number, not NaN)
- PE crabber-hours will be non-zero
- PE Dungeness catch via crabber-hr CPUE will be non-zero
- `IntA` (expansion interviews) will be >0, informing the gear-to-crabber relationship in BSS
- Stan model should compile and sample (O dimension now correct)

## Stan Model Notes

The BSS model will run but with some caveats:
- V_n=0, T_n=0, B_n=0: vehicle/trailer/boat loops don't execute
- R_V, R_T, R_B, b parameters still get sampled from priors (uninformed) — harmless but wasteful
- If convergence issues: increase `adapt_delta` to 0.95, `max_treedepth` to 14
