# Recreational Crab Creel Estimation — Grays Harbor / Westport

**Agency:** Washington Department of Fish and Wildlife (WDFW)  
**Lead:** Matt George  
**Status:** Proof of concept — 2024-25 season  

---

## What This Does

Estimates total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area by combining two statistical approaches:

- **Point Estimator (PE):** Simple stratified expansion — averages sampled days to estimate unsampled days. Fast, transparent, no modeling assumptions beyond "sampled days represent unsampled days."
- **Bayesian State-Space (BSS):** Fits a time-series model that treats daily effort and catch rate as hidden quantities evolving smoothly over time. Fills temporal gaps, propagates uncertainty, and produces credible intervals.

Three crabbing populations are estimated independently and summed for the port total:

1. **Shore crabbers** (dock + jetty + beach) — effort from gear counts, BSS + PE
2. **Private boat crabbers** — effort from trailer counts, BSS + PE
3. **Commercial/charter vessels** — effort from daily vessel tally, census expansion

---

## Two Model Variants

The framework includes two BSS model specifications. Both share the same effort model, data pipeline, and PE estimator. They differ in how catch-per-unit-effort (CPUE) is modeled:

### Pooled CPUE Model

All gear types (pots, ring nets, traps, snares) share a single CPUE process. Gear-type catch breakdowns are derived after estimation by applying interview-based proportions to the total. Simpler, faster, fewer parameters.

| File | Description |
|---|---|
| `BSS-GH-pooled-CPUE-model.Rmd` | R analysis script |
| `BSS-GH-pooled-CPUE-model-documentation.md` | Technical documentation |
| `stan_models/crab_bss_pooled.stan` | Stan model (single CPUE process) |

### Gear-Resolved CPUE Model

Each gear type gets its own CPUE process with independent AR(1) dynamics. Gear-type catch estimates carry posterior uncertainty directly from the model. Also includes a separate holiday effort effect (B2) and day-type stratified commercial/charter census expansion.

| File | Description |
|---|---|
| `BSS-GH-gear-type-CPUE-model.Rmd` | R analysis script |
| `BSS-GH-gear-type-CPUE-model-documentation.md` | Technical documentation |
| `stan_models/crab_bss_gear_resolved.stan` | Stan model (per-gear CPUE processes) |

### Which Should I Use?

Use the **pooled model** if you want simplicity, faster runtime (~3 hrs), and a single headline harvest number. Use the **gear-resolved model** if you need gear-type catch estimates with uncertainty, want to track whether pot CPUE differs from ring net CPUE over the season, or need the holiday effect and stratified census improvements (~4-5 hrs).

---

## Stan Model Naming

| Old Name | New Name | Description |
|---|---|---|
| `BSS_crab_model_01.stan` | (retired) | Original crab adaptation, single CPUE, no R_T guard |
| `BSS_crab_model_02.stan` | `crab_bss_pooled.stan` | Pooled CPUE, R_T guard, jetty reserved |
| `BSS_crab_model_03.stan` | `crab_bss_gear_resolved.stan` | Gear-type CPUE (G_gear processes), B2 holiday effect, Dirichlet gear proportions |

Both new names are accepted by the Rmd files via the `bss_model_file` parameter. Update the parameter to match whichever filename you use.

---

## Quick Start

1. Clone this repository
2. Place input data in `input_files/`:
   - `effort_combined.csv` (effort counts, re-exported with QUOTE_ALL)
   - `interview_combined.csv` (interviews, dates in M/D/YYYY format)
   - `wes_commercial_tally.csv` (daily vessel tally)
3. Place the Stan model file in `stan_models/`
4. Open the desired `.Rmd` file
5. Set `est_date_start` and `est_date_end` in the parameters section
6. Run all chunks (or Knit)

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here

---

## Output Files

Both models produce output in `output/YYYYMMDD/`. The gear-resolved model produces additional gear-type breakdowns:

### Both Models
| File | Contents |
|---|---|
| `pe_port_summary.csv` | PE estimates by component and port total |
| `port_total_Dungeness_Kept.csv` | Combined PE + BSS port total |
| `bss_summary_{label}.csv` | Stan convergence diagnostics per fit |
| `bss_daily_effort_{label}.csv` | Daily BSS effort (median + 95% CI) |
| `bss_daily_catch_{label}.csv` | Daily BSS catch (median + 95% CI) |
| `run_parameters.txt` | All parameters for reproducibility |

### Gear-Resolved Model Only
| File | Contents |
|---|---|
| `monthly_by_population.csv` | Month × Population with PE, BSS, and Combined estimates |
| `monthly_port_totals.csv` | Port-level monthly totals |
| `monthly_by_mode.csv` | Monthly catch by crabbing mode |
| `catch_by_gear_type_detail.csv` | Gear-type catch with BSS posterior uncertainty |
| `catch_by_gear_type.csv` | Port-level catch by gear type |
| `monthly_by_area.csv` | Monthly catch by creel area |
| `wide_summary.csv` | All dimensions in one table |

---

## Season Structure

The 2024-25 season (Sep 16, 2024 – Sep 15, 2025) is split into two independent sub-seasons at the pot-open date (Dec 1):

- **Ring-net only** (Sep 16 – Nov 30, 76 days): Ring nets, snares, foldable traps only
- **All-gear** (Dec 1 – Sep 15, 289 days): All gear including pots

Each sub-season gets its own BSS model fit. The split prevents the model from trying to bridge the structural break in effort and CPUE when pots become legal.

---

## Known Issues and Data Notes

- **Interview CSV column bug:** The `number_of_gear` column must be mapped from column N (not column W) in the raw iForm export due to a duplicate field name
- **Effort CSV quoting:** Re-exported with `QUOTE_ALL` to handle commas in the notes field
- **Interview dates:** Format is M/D/YYYY (read with `col_date(format="%m/%d/%Y")`)
- **Boat type typo:** iForm exports "Commerical" (one 'm') — handled by regex
- **Windows MAX_PATH:** If using OneDrive with long paths, the output directory may exceed 260 characters. The code detects this and falls back to a short path

---

## Development History

| Version | Model | Key Changes |
|---|---|---|
| v1 | — | Initial single-population dock-only prototype |
| v2 | — | Bug fixes (CSV columns, Stan O dimension, output folder) |
| v3 | Pooled | Three populations, two sub-seasons, convergence tuning |
| v4 | Pooled | Dawn/dusk day length, stat-week PE, census dates, team review |
| v5 | Gear-Resolved | Per-gear CPUE processes, B2 holiday effect, stratified census, comprehensive outputs |
