# 01_BSS_models

Production analysis drivers. These are the R Markdown files you actually run to produce a seasonal harvest estimate. Each driver wires together the inputs (`04_input_files/`), the helper functions (`03_R_functions/`), and a Stan model (`02_stan_models/`), then writes a full set of estimates and diagnostics to a dated folder under `05_output/`.

For the project-level overview (what the estimate is, the PE vs. BSS split, the three crabbing populations, the two sub-seasons), see the [root README](../README.md).

## Files

| File | Stan model used | Purpose |
|---|---|---|
| `BSS-GH-pooled-CPUE-model.Rmd` | `02_stan_models/crab_bss_pooled.stan` | Production driver with a single pooled CPUE process. Use for a headline harvest number. Gear-type catch is derived after estimation from interview proportions. |
| `BSS-GH-gear-type-CPUE-model.Rmd` | `02_stan_models/crab_bss_gear_resolved.stan` | Production driver with a per-gear CPUE process (shared AR(1) dynamics). Use when you need gear-type catch with posterior uncertainty, the `B2` holiday effort effect, or the stratified commercial/charter census expansion. |

`06_diagnostics/` holds the regression harness, the dated batch runners that validate changes to these drivers, and the experimental weather-tide covariate driver; none of them are production estimators. **For where the model currently stands** (the authoritative run, adopted configuration, open items), read the box at the top of `07_documentation/development_notes/PIPELINE_STATUS.md` and the change register `07_documentation/development_notes/CHANGE_REGISTER.md`.

## How a driver runs

1. **Setup chunk.** Loads packages, sets `run_date <- format(Sys.Date(), "%Y%m%d")`, sources every file in `03_R_functions/` with `purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)`, and sets `output_dir <- here("05_output", run_date, "<model>")`.
2. **Configuration.** User-selectable toggles (season window, structural dates, catch groups, effort unit, filters, I/E settings, holidays, and the model-behavior levers) live in `run_config.R` at the repository root, the single control surface for a run. The setup chunk sources `run_config.R` automatically when the orchestrator has not already defined `run_config`, then merges as `params <- modifyList(run_config, params_model)`: `run_config` is the BASE and the driver's own `params_model` list is layered ON TOP, so **a key present in both is silently decided by the driver**. That direction is load-bearing: `bss_min_interviews` once sat in `params_model` and silently overrode `run_config` until 2026-08-25, and `bss_sampler_override` exists as the one sanctioned way for `run_config` to reach a sampler key (it is applied AFTER the merge). `params_model` holds only this model's internal tuning: `bss_model_file` (the Stan filename), the per-fit sampler controls (`adapt_delta`, `max_treedepth`, iterations), the convergence-gate thresholds, and the AR-selector thresholds. For a routine run you edit `run_config.R`, not the `.Rmd`.
3. **Data prep.** Reads `effort_combined.csv`, `interview_combined.csv`, `wes_commercial_tally.csv` from `04_input_files/`, plus `ingress_egress.xlsx` (pooled). When `use_osp_boat_counts = TRUE`, both drivers also read `WBL_boat_counts.xlsx` (OSP daily boat totals) as a second boat-effort stream, and build the crabbing fraction f from the `boats_crabbing` / `boats_total` columns in `ingress_egress.xlsx`.
4. **Fit.** Each population x sub-season is fit independently by calling `rstan::stan(file = here("02_stan_models", params$bss_model_file), ...)`.
5. **Convergence gate, adequacy, and the optional AR ladder.** Each fit is checked by the scale-aware convergence gate (R-hat, n_eff, divergence fraction, and the SD-normalized divergence impact on the reported totals); a fit that fails falls back to its PE estimate. As of the 2026-09-07 configuration all four fitted components pass. The gate answers "did this fit sample", never "is this model right", so **model adequacy** (`p_loo` fraction, Pareto k count, randomized-PIT coverage, `flag_miscalibrated`) is reported BESIDE it in `model_adequacy.csv` and never gates. With `ar_escalate` on (off by default; every rung is a multi-hour fit) a component is fitted at a ladder of AR resolutions with per-rung adequacy logged to `ar_escalation_log.csv`; that machinery is how the shore all-gear cap was settled at weekly.
6. **Outputs.** Per-population daily series, port and monthly totals, PE-vs-BSS comparison, and a large set of diagnostics are written to `output_dir`. See [05_output/README.md](../05_output/README.md) for the file catalog.

## Path handling

All reads and writes use `here::here()`, which resolves to the repository root regardless of where the `.Rmd` sits or what the knit working directory is. That is why these drivers can live in `01_BSS_models/` while still reading from `04_input_files/` and writing to `05_output/`. If a stage folder is renamed, update the matching directory string inside the `here(...)` calls in these drivers (see "How paths work" in the root README).

## Companion documentation

Technical write-ups and change logs for each driver are in `07_documentation/`:

- `BSS-GH-pooled-CPUE-model-documentation.md`
- `BSS-GH-gear-type-CPUE-model-documentation.md`
