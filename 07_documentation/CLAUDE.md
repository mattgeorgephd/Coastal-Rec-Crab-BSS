# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

> ### PROJECT CONTEXT: nothing here has been published
>
> WDFW has published **no** recreational Dungeness crab harvest estimate from this pipeline.
> The `main` branch is the state of the model *before* a meeting with the WDFW freshwater
> creel team and *before* OSP confirmed they can supply daily boat counts; the
> `OSP-boat-count-incorporation` branch is the work incorporating both. **There is no
> published figure that a change has to stay consistent with**, so continuity with an
> earlier internal run is not on its own a reason to prefer one modelling choice over
> another. Judge changes on the evidence, and record what moved.
>
> **The OSP data** is, per day, the total number of returning vessels and the fraction that
> were *crabbing only*. The second deliberately excludes combo trips that also crabbed, so
> it is a **lower bound** on crabbing vessels, not the crabbing fraction `f`. That is what
> `osp_crab_lower` / `f_lower` are for; do not wire it in as if it were `f`. Its purpose is
> to improve the accuracy of the boat estimate and reduce its uncertainty. The acknowledged
> limitation is that more **boat interviews** are needed next season, which no amount of
> boat counting fixes.
>
> **The authoritative run** and its totals live in one place:
> `07_documentation/development_notes/PIPELINE_STATUS.md`, in the box at the top.


A WDFW recreational Dungeness crab creel-estimation pipeline for Grays Harbor / Westport (R + Stan). It estimates total seasonal harvest by fusing a design-based **Point Estimator (PE)** with a **Bayesian State-Space (BSS)** time-series model, across three crabbing populations, over two gear-regime sub-seasons. There are two production models (pooled and gear-resolved CPUE) and one experimental module (weather-tide covariates).

The repository is organized as a numbered stage pipeline: `01_BSS_models/` (drivers) → `02_stan_models/` (Stan code) → `03_R_functions/` (shared helpers) → `04_input_files/` (raw data) → `05_output/` (dated runs) → `06_diagnostics/` (experimental) → `07_documentation/` (reference layer). Most folders have their own `README.md` with a file inventory; this file covers what those don't, the cross-cutting architecture and the conventions that will bite you.

## Running the estimation

**One command runs everything.** Edit `run_config.R`, then launch:

```r
source("run_estimation.R")        # RStudio: use Source, NOT Knit
```
```sh
Rscript run_estimation.R                          # terminal / unattended
Rscript run_estimation.R --model gear_resolved    # override the model
Rscript run_estimation.R --model pooled --weather # also run the weather module
Rscript run_estimation.R --no-weather             # force weather off
```

CLI flags (`--model`, `--weather`, `--no-weather`) override `run_config.R` for that run.

- **In RStudio you must _Source_ `run_estimation.R`, not _Knit_ it**, knitting would try to render the script itself.
- A single model `.Rmd` can also be knit standalone; its setup chunk auto-sources `run_config.R` when `run_config` isn't already defined, so it uses identical toggles. You never edit the `.Rmd` for a routine run.
- **Runtime is long** (~3-6 h on 4 cores for a full pooled run, real MCMC over many fits).
- **Requirements:** R 4.2+, rstan 2.32+ (this uses **rstan, not cmdstanr**), plus tidyverse, lubridate, suncalc, gt, patchwork, here, readxl. `run_estimation.R` auto-`install.packages()` anything missing, so a fresh machine's first run may trigger a long compile. The weather module additionally needs mgcv, loo, httr, jsonlite, geosphere and reaches NOAA CO-OPS / NDBC / Iowa IEM endpoints at runtime (cached under `cache/`, git-ignored).

## `run_config.R` is the single control surface

This is the **one file you edit** for a routine run (do not edit `run_estimation.R`, the `.Rmd` drivers, or the `.stan` files). For a season re-run you typically change only:

- `model`, `"pooled"` or `"gear_resolved"`
- `run_weather`, `TRUE`/`FALSE` (**only valid with `model = "pooled"`**; the orchestrator hard-stops early otherwise, before any multi-hour fit, because the weather module reuses the pooled run's in-memory objects)
- the season window: `est_date_start`, `est_date_end`, and the structural dates (`pot_closure_start/end`, `pot_open_date`, `census_*`)
- the holiday list, now the `04_input_files/crabbing_holidays.xlsx` workbook (was the `crabbing_holiday_dates` config vector; read by `read_crabbing_holidays.R`), updated once per season

**Config keys added since this guide was first written**, all in `run_config.R` with comments at the definition: `shared_tau` (**TRUE in production since 2026-09-01**) and `shared_tau_sigma` / `shared_tau_min_obs` (the informed-day floor, 15, which is what makes the shared turnover boat-only); `estimate_catch_zi` / `catch_zi_populations` / `zi_catch_prior_a` / `zi_catch_prior_b` (the zero-inflated catch prototype, ships off); `pe_gear_ratio_arm` (PE/BSS incomplete-trip arm alignment, ships `"match_bss"`); `bss_sampler_override` (below); `bss_min_interviews_fitted`; `ar_escalate_ladder` / `ar_escalate_max_attempts` / `ar_escalate_respect_cap`; `crabbing_holidays_file` / `_sheet`.

`run_config` is a flat list; a key a given model doesn't read is silently ignored, which is why model-specific toggles (`collapse_mu_hier`, `estimate_B1_C`, `ar_adaptive`, `use_boat_ie`, …) all live in the one shared list. **Watch the merge direction:** each driver does `params <- modifyList(run_config, params_model)`, so `params_model` WINS. A key present in both is silently decided by the driver, not by `run_config`; that is exactly how `bss_min_interviews` sat at 20 while `run_config` appeared to own it (fixed 2026-08-25). Before adding a key to `run_config`, grep both `params_model` blocks for it. **Three toggles force a Stan recompile** when changed: `razor_dig_mode`, `estimate_cpue_density` and `estimate_catch_zi`. Several sensitivity levers ship commented-out or `NULL` for production (`R_G_prior_mu/sigma`, `ar_force`, `bss_sampler_override`). **`bss_sampler_override` is the one sanctioned way past the merge direction described above**: it is applied AFTER the merge, accepts sampler keys only, prints every change and errors on anything else (`03_R_functions/bss_sampler_override.R`). Use it rather than trying to set a sampler key in `run_config` and wondering why nothing changed. `bss_seed` is fixed for reproducibility, leave it fixed.

## How a run is wired (orchestration)

`run_estimation.R` sources the entire `03_R_functions/` library, sources `run_config.R`, applies CLI overrides, then renders the chosen driver `.Rmd`. Two non-obvious mechanics:

- **Config is injected via a shared environment, not rmarkdown `params:`.** The orchestrator builds `run_env <- new.env(parent = globalenv())`, sets `run_env$run_config`, and calls `rmarkdown::render(rmd, envir = run_env)`. Each driver then does `params <- modifyList(run_config, params_model)`, where `params_model` holds **only** that model's internal tuning (Stan filename, per-fit sampler settings, gate/AR thresholds). `run_config` and `params_model` are *intended* to be disjoint; when a key lands in both, the driver's value wins silently. That happened with `bss_min_interviews` (fixed 2026-08-25), so treat "intended" as a rule to check, not a guarantee.
- **The weather module shares that same `run_env`** ("Option A" hand-off): the pooled model's objects (`dwg`, `ie_data`, `L_eff_model`, …) satisfy the weather module's `if(!exists(...))` guards with no disk hand-off. Consequence: you cannot re-run the weather module without re-running the pooled model.

The **driver** (not the orchestrator) creates the output folder: `run_date <- format(Sys.Date(), "%Y%m%d")`, `output_dir <- here("05_output", run_date, "<model>-<run_tag>")`. The `run_tag` is set in `run_config.R` (or by `run_rg_sweep.R`) and names the folder (e.g. `pooled-CPUE-run1`); when it is blank the driver appends an `HHMMSS` timestamp instead, so same-day re-runs of the same model land in **distinct** folders and no longer overwrite each other. The orchestrator reads `output_dir` back out of `run_env`, moves the rendered HTML into it, and writes `run_manifest_<timestamp>.txt` one level up in `05_output/<run_date>/` (recording model, `run_weather`, git SHA, per-stage timing, a `str()` dump of `run_config`, and full `sessionInfo()`). The dated folder is **today's system date at render time**.

## The estimation architecture (big picture)

**Three populations, estimated independently and summed into a port total:**

1. **Shore** (dock + jetty + beach), effort from gear counts; BSS + PE.
2. **Private boat**, effort from trailer counts; BSS + PE.
3. **Commercial/charter**, **not modeled**; a day-type-stratified **census expansion** of the daily vessel tally (`estimate_comm_charter.R`). No BSS fit, no per-population output file. This asymmetry is intentional.

**PE and BSS are fused per fit by the convergence gate** (`bss_convergence_gate.R`), which is the single authority on method selection. A fit reports its BSS posterior only if it passes **all** of: R-hat < 1.01, n_eff > 400, divergent fraction < 0.05, **and** an SD-normalized divergence-impact test (< 0.10 posterior SD, does the divergence *move* the answer, per Betancourt 2017). Otherwise that population × sub-season contributes its **PE point** instead. PE and BSS are always reported side-by-side (`pe_vs_bss_comparison.csv`, `convergence_report.csv`).

- **`bss_use_pe_for(b)` is the only correct way to ask "should this component report PE?"** in a totals section, it combines pre-fit data-sufficiency (`pe_fallback`) with the gate result (`use_bss`). Diagnostic sections deliberately check `b$pe_fallback` *alone* so that a fitted-but-gate-failed component keeps its per-fit diagnostics, do not "fix" those to use `bss_use_pe_for`.

**Two sub-seasons, fit separately, split at the pot-open date** (`build_subseasons.R`): a pot-closure sub-season (non-pot gear only) and an all-gear sub-season. The split stops the model from bridging the structural break when pots become legal; totals sum over sub-seasons. Note the deliberate key/display mismatch: the pot-closure sub-season's internal key is `ring_net_only` (kept for output-filename continuity) even though it's displayed as "Pot closure".

**Pooled vs gear-resolved differ only in how CPUE is modeled.** Pooled uses one CPUE process and allocates gear-type catch after estimation via Dirichlet-propagated interview shares. Gear-resolved is *written* for a per-gear CPUE process, but **as currently driven it runs with `G = 1`, so that machinery is inert** and gear-type catch is PE-apportioned. Do not raise `G > 1` without adding per-gear effort shares; only gear 1 is observed in the effort stream.

## `03_R_functions/` contract (load-bearing)

Every driver sources the **entire folder** wholesale in its setup chunk:
```r
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)
```
Because both production drivers source every file, two rules are structural, not stylistic:

1. **Pure functions, zero source-time side effects.** No top-level reads/writes/plots/`library()` calls, no assumption about sourcing order. Config is *passed* through `params`, never captured from a driver global. (`bss_timers.R` keeps state in a module-local environment precisely to avoid needing a global.)
2. **No name collisions.** Shared functions keep one name and one body. Functions that differ between tracks carry a `_pooled` / `_gear` suffix (`run_pe_pooled` vs `run_pe_gear`; `prep_bss_crab_pooled` vs `prep_bss_crab_gear`). **Adding a track-specific helper without a suffix silently overwrites the other track's function** depending on alphabetical source order. (`fetch_crab_data` was merged 2026-08-01 into one shared reader parameterized by `boat_require_gear_time`, retiring the former `_v2` variant; prefer this pattern, a shared function with a per-model param, over a name split when the two versions are nearly identical.)

Key shared helpers to know: `bss_model_adequacy.R` (model adequacy reported BESIDE the gate, never gating it), `bss_sampler_override.R`, `pe_gear_ratio_frame.R`, `annotate_decoupled_run.R`, `fetch_osp_boat_counts.R`, `diagnose_osp_trailer_overlap.R`, `bss_effort_spec.R` (single source of the effort unit AND the matching I/E observation column, read by both PE and BSS prep so they can't drift), `bss_convergence_gate.R` (the gate), `bss_ar_resolution.R` (adaptive AR selector + `bss_ar_ladder()`, the opt-in escalation ladder), `bss_day_length.R` (I/E ingest + `L_effective`), `crab_fraction.R` (the crabbing fraction `f`, its per-stratum priors, and the OSP crab-only lower bound), `bss_opener_covariates.R` (opener effort covariates and their multiplicity-controlled screen), `diagnose_incomplete_trips.R` (the four-arm treatment diagnostic), `build_subseasons.R`, `prep_days_crab.R`, `prep_population_summary.R`.

## Cross-cutting invariants

- **Paths: everything resolves through `here::here()`**, anchored to the repo root via the `.Rproj` / `.git` sentinels, *not* the `.Rmd` location, so a driver knits correctly from any working directory. But the directory **names** inside `here(...)` must match the on-disk folder names exactly; renaming a numbered stage folder requires updating every `here()` string in the drivers and diagnostic writers.
- **Effort unit is gear-deployments** for both shore and boat (as of v7.7 / v7.6). Time-denominated units (crabber-hours, gear-hours) are invalid for pot/trap gear because catch is sub-linear in soak time; the pipeline re-measures this every run (`cpue_linearity_*` / `cpue_saturation_*` CSVs). Any new gear type must pass those before its totals are trusted. A "gear-deployment" is a piece of gear a crabber had in the water on that trip (`number_of_gear`); **not** a pot lift; repeat checks of one gear slot are not extra deployments, and slot re-use across the day is carried separately by `tau`.
- **The daily expansion factor `L` is a TURNOVER, not a day length.** `L = tau_shore` (~1.7) / `tau_boat` (~1.2). `L_effective` in hours (~5.3) is computed every run but is a diagnostic; it only becomes `L` if `shore_effort_unit` is set back to a time unit. Anything that consumes `L`; including the I/E observation stream; must carry the matching unit. `bss_effort_spec()` is the single source for that pairing (`h_col`, `L_data`, and as of 2026-08-25 `ie_obs_col`); if you add a consumer of effort, route it through that function rather than re-deriving the formula. Four separate places had drifted onto the pre-v7.7 crabber-hours formula before the 2026-08-25 audit, and one of them (the shore I/E observation) was materially wrong.
- **AR temporal resolution is an inference lever, not a tuning knob.** It's selected per fit from data density (daily/weekly/monthly) then capped per population via `ar_max_resolution`; too-fine AR is unidentified (the pooled boat's daily AR once diverged ~100%). As of 2026-08-25 an opt-in escalation ladder (`ar_escalate`, default FALSE) can instead start at the finest rung and coarsen only on a convergence-gate failure; it costs one extra multi-hour fit per failed rung, which is why it ships off.

## Input data quirks (real, not bugs, do not "fix")

From `04_input_files/` (see its README). As of the 2026-07-16 input migration every model input is an **`.xlsx` workbook with a single `data` sheet** (previously a mix of CSVs), read via `readxl::read_excel`, with dates stored as ISO `yyyy-mm-dd` text parsed by `as.Date()`. Quirks that survive re-export and are matched in code:
- Interview `number_of_gear` maps from **column N, not W** (duplicate iForm field name).
- The commercial `boat_type` is the typo **"Commerical"** (one m), matched by regex, don't correct the spelling without updating the matcher.
- Windows/OneDrive long paths can exceed MAX_PATH; the code detects this and falls back to a short path.
- Historical note (pre-migration, no longer applies): the CSVs needed `QUOTE_ALL` for the notes field, and carried **M/D/YYYY** dates that silently parsed to `NA` under the default reader; the xlsx conversion fixed both, and the notes / unused columns were dropped in the same modeling-only strip.

`ingress_egress.xlsx` is read by **both** drivers (via `fetch_ie_data` in `bss_day_length.R`, called from the pooled and gear-resolved drivers alike), feeding the shore `L_effective` day-length model; the earlier "pooled driver only" claim was wrong (input-audit F1). The one genuinely pooled-only input is the opener calendar `fishery_opener_dates.xlsx`, read solely by the pooled report's spillover diagnostic, so it changes no estimate.

## Outputs

Each run writes to `05_output/<YYYYMMDD>/<model>-<run_tag>/`, where the `<model>` prefix is `pooled-CPUE`, `gear-type-CPUE-model`, or `pooled-CPUE-covariates` and `run_tag` (from `run_config.R`, or `run_rg_sweep.R`) names the folder, e.g. `pooled-CPUE-run1`; a blank `run_tag` falls back to an `HHMMSS` timestamp, so same-day runs of a model no longer overwrite. Per-population files follow `<metric>_<population>_<species>_<fate>.{csv,png}`; port/monthly/comparison files drop the population tag. **Outputs are committed to git on purpose** so past estimates are preserved as-produced; only per-run `*.RData` workspaces and compiled `*.rds` Stan caches are git-ignored. A given dated folder may be a *partial* run; use a recent complete run as the reference catalog, not any single folder.

## Working conventions

- **Validate by run, never by reasoning alone.** A change that looks inference-neutral on paper can perturb the sampler geometry. Isolate each change, compare against a confirmed baseline run against pre-set criteria; an item is not "done" until a run confirms it. Changes that *narrow* reported uncertainty (e.g. tightening dispersion priors) need explicit sign-off because they move the headline intervals.
- **Diagnostics are additive and `tryCatch`-wrapped** so one failing fit can't abort a run.
- **`07_documentation/development_notes/PIPELINE_STATUS.md` is the single living status + backlog document** (forward-looking; last updated 2026-09-03). **The authoritative run and its port total are in a box at the top of that file**; do not take an estimate from anywhere else in the repo without checking it there first, because several documents still quote superseded totals. Read it first for current state and update it in place; don't fork per-session notes. The two `*-development-history.md` docs are the backward-looking, newest-first change logs.
- **Method version ≠ code version.** "Method v1.0" is the frozen method-of-record; the pipeline code (v7.x) moves faster and has changed published totals since the freeze (v7.5 incomplete-trip filter, v7.6/v7.7 gear-deployment switch, v7.8 refactor). Numbers predating the 2026-07-11 refresh are labeled "pre-refresh", don't cite them as final.
- `06_diagnostics/` holds **the validation harness and the batch runners**, not only the weather module: `test_improvements_2026-08-25.R` is a dependency-light regression harness (279 assertions, seconds, no rstan) that should be run before committing to any long fit, and `run_*.R` are the dated experiment batches. Every batch runner ships `DRY_RUN <- TRUE`; that is asserted by the harness. The weather-tide module in the same folder is **experimental and currently stale** (its Stan fork is missing ~40 data variables the production model now declares); its committed conclusion is that covariates are excluded, and it is not part of the official harvest estimate.
