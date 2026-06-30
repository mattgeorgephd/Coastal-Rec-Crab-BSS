# 01_BSS_models

Production analysis drivers. These are the R Markdown files you actually run to produce a seasonal harvest estimate. Each driver wires together the inputs (`04_input_files/`), the helper functions (`03_R_functions/`), and a Stan model (`02_stan_models/`), then writes a full set of estimates and diagnostics to a dated folder under `05_output/`.

For the project-level overview (what the estimate is, the PE vs. BSS split, the three crabbing populations, the two sub-seasons), see the [root README](../README.md).

## Files

| File | Stan model used | Purpose |
|---|---|---|
| `BSS-GH-pooled-CPUE-model.Rmd` | `02_stan_models/crab_bss_pooled.stan` | Production driver with a single pooled CPUE process. Use for a headline harvest number. Gear-type catch is derived after estimation from interview proportions. |
| `BSS-GH-gear-type-CPUE-model.Rmd` | `02_stan_models/crab_bss_gear_resolved.stan` | Production driver with a per-gear CPUE process (shared AR(1) dynamics). Use when you need gear-type catch with posterior uncertainty, the `B2` holiday effort effect, or the stratified commercial/charter census expansion. |

The experimental weather-tide covariate driver lives separately in `06_diagnostics/`, not here, because it is a research/diagnostic tool rather than a production estimator.

## How a driver runs

1. **Setup chunk.** Loads packages, sets `run_date <- format(Sys.Date(), "%Y%m%d")`, sources every file in `03_R_functions/` with `purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)`, and sets `output_dir <- here("05_output", run_date, "<model>")`.
2. **Parameters.** A `params <- list(...)` chunk holds the run configuration. The ones you normally touch:
   - `est_date_start`, `est_date_end`: the season (or sub-season) window.
   - `bss_model_file`: the Stan filename (just the name; the folder `02_stan_models/` is supplied by the `here()` call at the fit step).
   - `ie_data_file`, `ie_sheet`: the I/E workbook and sheet (pooled driver only).
   - Sampler controls (`adapt_delta`, `max_treedepth`, iterations) and the convergence-gate thresholds.
3. **Data prep.** Reads `effort_combined.csv`, `interview_combined.csv`, `wes_commercial_tally.csv` from `04_input_files/`, plus `ingress_egress.xlsx` (pooled).
4. **Fit.** Each population × sub-season is fit independently by calling `rstan::stan(file = here("02_stan_models", params$bss_model_file), ...)`.
5. **Convergence gate.** Each fit is checked (R-hat, divergences, treedepth). A fit that fails the gate falls back to its PE estimate for that population; the boat all-gear fit is the usual fallback case.
6. **Outputs.** Per-population daily series, port and monthly totals, PE-vs-BSS comparison, and a large set of diagnostics are written to `output_dir`. See [05_output/README.md](../05_output/README.md) for the file catalog.

## Path handling

All reads and writes use `here::here()`, which resolves to the repository root regardless of where the `.Rmd` sits or what the knit working directory is. That is why these drivers can live in `01_BSS_models/` while still reading from `04_input_files/` and writing to `05_output/`. If a stage folder is renamed, update the matching directory string inside the `here(...)` calls in these drivers (see "How paths work" in the root README).

## Companion documentation

Technical write-ups and change logs for each driver are in `07_documentation/`:

- `BSS-GH-pooled-CPUE-model-documentation.md`
- `BSS-GH-gear-type-CPUE-model-documentation.md`
