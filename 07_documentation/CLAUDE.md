# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository. It lives in `07_documentation/`; the root `CLAUDE.md` imports it with `@07_documentation/CLAUDE.md`, so Claude Code loads it from the repository root too. Edit this file, not the root one.

## What this is

> ### PROJECT CONTEXT: nothing here has been published
>
> WDFW has published **no** recreational Dungeness crab harvest estimate from this pipeline.
> Until 2026-09-28 the `main` branch was the state of the model *before* a meeting with the
> WDFW freshwater creel team and *before* OSP confirmed they can supply daily boat counts
> (`724eead`), and the `OSP-boat-count-incorporation` branch was the work incorporating both.
> That branch merged into `main` on 2026-09-28 (pull request #5, merge commit `a878a87`), and
> pull requests #6 (2026-09-28) and #7 (2026-09-29, merge `3a55a2e`) followed, so
> `main` now carries Method v2.0; the merge description is archived at
> `07_documentation/archive/PR-5-merge-OSP-boat-count-incorporation.md`. **There is no
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
> **THE PROGRAM GOAL (stated by Matt, 2026-09-09).** The OSP branch was to supersede `main`,
> and did (2026-09-28).
> The 2024-25 season was the DEVELOPMENT TEST SEASON, not the target: the model must run
> on any season the user selects, a part of a season, or a multi-season span, and on a
> naive new window it must hand the user the information to tune the fit (the AR ladder
> finds the minimum resolution at which each fit passes the gate; per-rung adequacy says
> what a finer or coarser choice gains and loses). Do not optimize for 2024-25; build and
> document for the novel-season workflow. That workflow is
> `07_documentation/NEW_SEASON_GUIDE.md`, and the 2024-25-derived settings are tagged
> SEASON-DERIVED in `run_config.R`.
>
> **The authoritative run** and its totals live in one place:
> `07_documentation/development_notes/PIPELINE_STATUS.md`, in the box at the top.
> **Every change since the OSP work began, and its status** (ADOPTED / BUILT, INERT / OPEN /
> REJECTED / BLOCKED) is tabulated in
> `07_documentation/development_notes/CHANGE_REGISTER.md`. The adopted configuration includes
> the shared boat turnover, the WEEKLY shore all-gear AR cap and the zero-inflated shore
> catch likelihood from 2026-09-07 (on BOTH tracks since 2026-10-02, A32: `catch_zi_tracks` ships `c("pooled", "gear_resolved")`, adopted together with the gear track's per-population AR periods, D3), and then
> four changes that between them moved the port total **+31%** in the 2026-09-11 ladder: a
> **dynamic monthly crabbing fraction `f`** (a per-stratum logit random walk observed by the
> sampler boat contacts, replacing a flat 0.30 that was wrong by roughly 3x in the winter
> months that carry most of the boat catch), the **boat turnover recentred on the OSP/trailer
> calibration**, the **shore turnover derived from the I/E `time` column**, and the census
> split into a **commercial CENSUS plus a charter EXPANSION**; and, since 2026-09-27, the **NWS
> Small-Craft-Advisory day flag on the private-boat all-gear effort process** (A30, one
> season-constant term, `marine_hazard_mode = "manual"`; boat all-gear 47,319 against 45,604; first
> rendered 2026-09-27 at port 96,118). **The authoritative run** is the method of record rendered at
> the A31 code with `init_r = 0.5`, stage A of `06_diagnostics/run_authoritative_batch_2026-09-29.R`
> (B58 to B62): `05_output/20260929/pooled-CPUE-2024-25-220449`, port **99,822 [82,090, 124,718]**,
> rendered by Matt from `2523e8e` on a clean tree with every package at its `renv.lock` version,
> and judged by rule A, written in the runner's header before it ran. Every component is within 0.02
> posterior SD of the 2026-09-28 run it replaced (99,873, B50); what A31 changed is the sampler (the
> shore all-gear fit 0.49% divergent, against 4.07%; D33 closed). Its gear-resolved cross-check is
> `05_output/20260930/gear-type-CPUE-model-2024-25-AB-D6-gear-matched-zi` at **99,294 (-0.53%, PASS)**
> since A32 (2026-10-02): the same batch's stage D6, which is the shipped file's gear render once D3
> (the gear shore fits at the pooled periods) and D6 (the gear-track zero-inflated catch) were
> adopted TOGETHER. D3 alone trapped a chain at the shipped `bss_seed`; with D6 it passed at both
> seeds tried. The same batch closed D29 by rule: the boat all-gear fit stays monthly, and its
> adequate rungs span 47,105 to 53,651 (full-posterior medians; 47,192 to 53,546 on the batch's
> 2,000 saved draws), so the port reads 99,822 at monthly and 106,380 at weekly
> (+6.6%), reported as resolution uncertainty, not inside the interval. R2 showed the A31 model needs
> no small radius. The batch is finished and its runner now refuses to fit (B62).
> If a number in this file disagrees with the box, the box wins.


A WDFW recreational Dungeness crab creel-estimation pipeline for Grays Harbor / Westport (R + Stan). It estimates total seasonal harvest by fusing a design-based **Point Estimator (PE)** with a **Bayesian State-Space (BSS)** time-series model, across three crabbing populations, over two gear-regime sub-seasons. There are two models: pooled CPUE (the headline estimator) and gear-resolved CPUE (its cross-check). The weather-tide covariate module was removed 2026-09-13 (A29); its finding, exclusion, stands.

The repository is organized as a numbered stage pipeline: `01_BSS_models/` (drivers) → `02_stan_models/` (Stan code) → `03_R_functions/` (shared helpers) → `04_input_files/` (raw data) → `05_output/` (dated runs) → `06_diagnostics/` (harness and batch runners) → `07_documentation/` (reference layer). Most folders have their own `README.md` with a file inventory; this file covers what those don't, the cross-cutting architecture and the conventions that will bite you.

## Running the estimation

**One command runs everything.** Edit `run_config.R`, then launch:

```r
source("run_estimation.R")        # RStudio: use Source, NOT Knit
```
```sh
Rscript run_estimation.R                          # terminal / unattended
Rscript run_estimation.R --model gear_resolved    # override the model
Rscript run_estimation.R --model both             # pooled headline, then the gear cross-check
```

`--model` (also `--model=X`) is the only CLI flag and it overrides `run_config.R` for that run; it takes `pooled`, `gear_resolved` or `both`, and since B46 (2026-09-28) ANY other argument stops the run. `both` renders each model in a fresh environment, writes `cross_check_<timestamp>.csv` (the two `Expected_Catch` medians and PASS or REVIEW against `cross_check_tolerance`, 0.02), and a failed stage is recorded in the manifest before the orchestrator exits non-zero. `--weather` and
`--no-weather` were removed with the weather-tide module on 2026-09-13 (A29); passing either
now stops the run with a message rather than being silently ignored.

- **In RStudio you must _Source_ `run_estimation.R`, not _Knit_ it**, knitting would try to render the script itself.
- A single model `.Rmd` can also be knit standalone; its setup chunk auto-sources `run_config.R` when `run_config` isn't already defined, so it uses identical toggles. You never edit the `.Rmd` for a routine run.
- **Runtime is long** (~3-6 h on 4 cores for a full pooled run, real MCMC over many fits).
- **Requirements:** R 4.2+, rstan 2.32+ (this uses **rstan, not cmdstanr**) and a C++ toolchain. Packages load through ONE function, `bss_load_packages()` (`03_R_functions/bss_packages.R`; the list is `bss_required_packages`: tidyverse, lubridate, rstan, here, readxl, rmarkdown, knitr, loo, digest; suncalc left it on 2026-09-29, B55, because only a time-denominated shore effort unit reads the civil-twilight day length it computes, and it stays in `renv.lock` for that path). `renv.lock` pins the full closure (120 packages, B46); the committed `.Rprofile` activates renv, so `renv::restore()` once from a fresh clone (README, "Setting up R"). With renv active a missing package is restored from the lock; without it, installed from CRAN with a message that the versions are not pinned. Set `RENV_ACTIVATE_PROJECT=FALSE` to run against a site library (the harness in a container). Do not add a package with `library()` / `install.packages()` in a driver or runner: add it to `bss_required_packages` and the lock. Since the weather-tide module was removed (2026-09-13) a run needs NO network access: mgcv, httr, jsonlite and geosphere, and the NOAA CO-OPS / NDBC / Iowa IEM calls, went with it. That still holds after the marine hazard covariates (2026-09-25): a run reads the committed `04_input_files/nws_marine_hazards.xlsx`; only its builder, `04_input_files/build_nws_marine_hazards.R`, calls the IEM service (and needs `jsonlite`), and it is run by hand, not by a run and not by `build_all_inputs.R`. Since 2026-09-27 EVERY production run reads that workbook (the covariate is the method of record) and STOPS if it does not cover the window, so rebuilding it is a per-season step (`NEW_SEASON_GUIDE.md`, section 0).

## `run_config.R` is the single control surface

This is the **one file you edit** for a routine run (do not edit `run_estimation.R`, the `.Rmd` drivers, or the `.stan` files). For a season re-run you typically change only:

- `model`, `"pooled"`, `"gear_resolved"` or `"both"`
- `run_tag`: letters, digits, `-` and `_` (anything else is replaced); a season token in it (`2025-26`) must be in `season_filter` or be the span label, or the run stops; a tag that already names a non-empty folder gets `-HHMMSS` appended (`bss_output_dir()`, B46)
- the season window: `est_date_start`, `est_date_end`, and the structural dates (`pot_closure_start/end`, `pot_open_date`, `census_*`)
- the holiday list, now the `04_input_files/crabbing_holidays.xlsx` workbook (was the `crabbing_holiday_dates` config vector; read by `read_crabbing_holidays.R`). It is BUILT by `04_input_files/build_crabbing_holidays.R`: for a new season extend the builder's season/year range and rebuild, never hand-edit the workbook (Matt updates the holiday input for new seasons as needed)

**Config keys added since this guide was first written**, all in `run_config.R` with comments at the definition: `shared_tau` (**TRUE in production since 2026-09-01**) and `shared_tau_sigma` / `shared_tau_min_obs` (the informed-day floor, 15, which is what makes the shared turnover boat-only); `estimate_catch_zi` / `catch_zi_populations` / `zi_catch_prior_a` / `zi_catch_prior_b` (the zero-inflated shore catch likelihood, ADOPTED 2026-09-07: `estimate_catch_zi = TRUE`, `catch_zi_populations = "shore"`; `catch_zi_tracks` decides which tracks fit it and ships `c("pooled", "gear_resolved")` since 2026-10-02, A32); `gear_period_bss` (the gear track's AR period per population; since A32 the shore at the pooled track's periods, the boat as before); `pe_gear_ratio_arm` (PE/BSS incomplete-trip arm alignment, ships `"match_bss"`); `bss_sampler_override` (below); `bss_min_interviews_fitted`; `ar_escalate_ladder` / `ar_escalate_max_attempts` / `ar_escalate_respect_cap`; `crabbing_holidays_file` / `_sheet`; `mu_hier_collapse_single` (A31, 2026-09-29: `"both"` ships; which single-section level hierarchies the pooled Stan collapses; Stan data, no recompile); `day_length_diagnostics` (B55: computes `L_effective` and civil twilight under the deployment unit anyway; ships FALSE); and the marine hazard effort covariate (CHANGE_REGISTER A30): **the method of record since 2026-09-27** is `marine_hazard_mode = "manual"`, `marine_hazard_manual_boat = "nws_sca_any"`, `marine_hazard_manual_shore = character(0)` and `marine_hazard_gear_regimes = "all_gear"` (section 2.10, with the flag definition `marine_hazard_file` / `_sheet` / `_zones` / `_codes` / `_window` / `_tz`: what an SCA-or-higher day IS, the codes, the two NWS zones, the local 04:00 to 16:00 window). One season-constant term on the private-boat ALL-GEAR effort process; no shore term, not the samplers' bar-restriction tick, not the pot-closure fit. The experiment surface stays in section 4.4b: `marine_hazard_candidates_shore` / `_boat` (what `"auto"` / `"on"` consider; a manual term must be offered there), `marine_hazard_auto_p` / `_p_adjust`, `marine_hazard_winter_months` (the season-split candidates, B41), `bar_restriction_impute` / `_field_start`. Under `"off"` nothing is read and both preps build the pre-2026-09-27 Stan data, measured. The evidence is campaign Sections 1x to 1z; the limitations are D30 and D32. Do not change the term's scope in passing: `marine_hazard_terms_for()` is where a selected term is confined to its sub-seasons.

`run_config` is a flat list; a key a given model doesn't read is silently ignored, which is why model-specific toggles (`collapse_mu_hier`, `estimate_B1_C`, `ar_adaptive`, `use_boat_ie`, …) all live in the one shared list. **Watch the merge direction:** each driver does `params <- modifyList(run_config, params_model)`, so `params_model` WINS. A key present in both is silently decided by the driver, not by `run_config`; that is exactly how `bss_min_interviews` sat at 20 while `run_config` appeared to own it (fixed 2026-08-25). Before adding a key to `run_config`, grep both `params_model` blocks for it. **Three toggles force a Stan recompile** when changed: `razor_dig_mode`, `estimate_cpue_density` and `estimate_catch_zi`. Several sensitivity levers ship commented-out or `NULL` for production (`R_G_prior_mu/sigma`, `ar_force`, `bss_sampler_override`). **`bss_sampler_override` is the one sanctioned way past the merge direction described above**: it is applied AFTER the merge, accepts sampler keys only, prints every change and errors on anything else (`03_R_functions/bss_sampler_override.R`). Use it rather than trying to set a sampler key in `run_config` and wondering why nothing changed. `bss_seed` is fixed for reproducibility, leave it fixed. Since 2026-09-26 (B38) `bss_stan_fit()` also rebuilds from it the draw permutation rstan had drawn from R's RNG in each chain's process, so the port total, which pairs component draws row by row, is reproducible between identical runs; before B38 it drifted about 0.1 to 0.2% between bit-identical fits, which is why R4's 94,376 will not be hit exactly by a post-B38 render of the same configuration. Two post-B38 renders on one machine have written byte-identical saved draws (the boat fits of 2026-09-27 and 2026-09-28); across platforms or compilers Stan's floating point, and so the last digits, may differ. **`bss_init_r` (0.5, B51, 2026-09-29) is also fixed and should stay fixed**: it is the radius rstan draws each chain's starting point from (default 2), set small because a shore all-gear chain started in the `sigma_mu_E` funnel can stay stuck there (D33). The same day A31 removed that funnel (`mu_hier_collapse_single = "both"` collapses the unidentified single-section effort and CPUE level hierarchies; `"effort"` alone left 4.7% divergences at the shipped radius; `"none"` restores the old model exactly), so the radius is now a precaution; keep it. The authoritative run (99,822, B61) is the first render at the A31 code with `init_r = 0.5`, and it is the reproducible reference: a render of `run_config.R` as shipped on the same machine reproduces it. The same batch rendered the model at 2 (R2): every gate passed and every component sat within 0.08 posterior SD, so the radius is a precaution, not something the result depends on. The first render at 0.5, Matt's of the pre-A31 code (`05_output/20260928/pooled-CPUE-2024-25-222347`, B57), failed the shore pot-closure fit's gate with a chain stuck in the level funnels and was never the authoritative run.

## How a run is wired (orchestration)

`run_estimation.R` sources the entire `03_R_functions/` library, sources `run_config.R`, applies CLI overrides, then renders the chosen driver `.Rmd` (or both, in turn, each in its own fresh `run_env`). Two non-obvious mechanics:

- **Config is injected via a shared environment, not rmarkdown `params:`.** The orchestrator builds `run_env <- new.env(parent = globalenv())`, sets `run_env$run_config`, and calls `rmarkdown::render(rmd, envir = run_env)`. Each driver then does `params <- modifyList(run_config, params_model)`, where `params_model` holds **only** that model's internal tuning (Stan filename, per-fit sampler settings, gate/AR thresholds). `run_config` and `params_model` are *intended* to be disjoint; when a key lands in both, the driver's value wins silently. That happened with `bss_min_interviews` (fixed 2026-08-25), so treat "intended" as a rule to check, not a guarantee.
- **`run_env` is one per model stage.** It used to be shared with the weather-tide module ("Option A" hand-off), which read the pooled driver's in-memory objects (`dwg`, `ie_data`, `L_eff_model`, …) rather than anything on disk. That module was removed on 2026-09-13 (A29), so the coupling is gone; since B46 `model = "both"` renders the two production drivers in turn, each in its own fresh `run_env`, so neither can read the other's objects. Note that the coupling was one-directional: no production driver ever referenced the weather module, which is why the removal touched no fit.

The **driver** (not the orchestrator) creates the output folder through `bss_output_dir("<model>-", run_config)` (B46), which resolves to `here("05_output", format(Sys.Date(), "%Y%m%d"), "<model>-<run_tag>")` after the `run_tag` checks above. The `run_tag` is set in `run_config.R` (or by `run_rg_sweep.R`) and names the folder (e.g. `pooled-CPUE-run1`); when it is blank the driver appends an `HHMMSS` timestamp instead, so same-day re-runs of the same model land in **distinct** folders and no longer overwrite each other. The orchestrator reads `output_dir` back out of `run_env`, moves the rendered HTML into it, and writes `run_manifest_<timestamp>.txt` one level up in `05_output/<run_date>/` (recording model, git SHA and whether the working tree matched it, per-stage timing, a `str()` dump of EVERY `run_config` key, and full `sessionInfo()`; until 2026-09-28 the dump stopped at 99 keys and there was no tree line, B43). The run folder's own `run_parameters.txt` has recorded the full resolved parameters since 2026-09-04 (it had the same 99-key truncation until then); it is the record to read. The dated folder is **today's system date at render time**.

## The estimation architecture (big picture)

**Three populations, estimated independently and summed into a port total:**

1. **Shore** (dock + jetty + beach), effort from gear counts; BSS + PE.
2. **Private boat**, effort from trailer counts; BSS + PE.
3. **Commercial/charter**, **not modeled**; **two different estimators in one component** (`estimate_comm_charter.R`), split 2026-09-11 on Matt's answer that the census premise holds for the commercial vessels only. The **commercial** part is an EXACT CENSUS summed over the tally days (`census_expansion = "none"`: the unsampled days had no commercial operation by the sampling design, so there is nothing to expand to). The **charter** part is an EXPANSION over the charter trip roster stratified by vessel (`charter_frame = "roster"`, `charter_expansion = "vessel"`), because charter trips are not fully sampled (34 trips, 20 interviewed in 2024-25), with SRS-with-FPC variance `N^2 (1 - n/N) s^2 / n`; `census_uncertainty = "charter"` carries that SE into the port interval. The residual it does NOT carry: the commercial per-vessel catch MEAN is itself a sample (141 of 164 vessel-trips in 2024-25, 44 of 67 in 2025-26), so "census without error" is a statement about the COUNT, not the catch (D16). No BSS fit, no per-population output file. This asymmetry is intentional. **This FPC is the only one in the pipeline**; the PE applies none anywhere (D22), so the two estimators are knowingly inconsistent on that point.

**PE and BSS are fused per fit by the convergence gate** (`bss_convergence_gate.R`), which is the single authority on method selection. A fit reports its BSS posterior only if it passes **all** of: R-hat < 1.01, n_eff > 400, divergent fraction < 0.05, **and** an SD-normalized divergence-impact test (< 0.10 posterior SD, does the divergence *move* the answer, per Betancourt 2017). Otherwise that population × sub-season contributes its **PE point** instead. PE and BSS are always reported side-by-side (`pe_vs_bss_comparison.csv`, `convergence_report.csv`).

- **`bss_use_pe_for(b)` is the only correct way to ask "should this component report PE?"** in a totals section, it combines pre-fit data-sufficiency (`pe_fallback`) with the gate result (`use_bss`). Diagnostic sections deliberately check `b$pe_fallback` *alone* so that a fitted-but-gate-failed component keeps its per-fit diagnostics, do not "fix" those to use `bss_use_pe_for`.

**Sub-seasons, fit separately, split at the pot-open date** (`build_subseasons.R`): a pot-closure sub-season (non-pot gear only) and an all-gear sub-season. The split stops the model from bridging the structural break when pots become legal; totals sum over sub-seasons. Note the deliberate key/display mismatch: the pot-closure sub-season's internal key is `ring_net_only` (kept for output-filename continuity) even though it's displayed as "Pot closure".

**It is no longer two.** Since 2026-09-10 (A14) a run may span MULTIPLE SEASONS: `pot_closures` takes one closure window per season and `build_subseasons()` emits `ring_net_only_<season>` / `all_gear_<season>` per closure (plus `all_gear_pre_<season1>` for any window before the first), with per-season census windows in `census_windows` and `season_totals.csv` written every run. A config with 0 or 1 closure falls through to the byte-identical scalar path, so single-season runs are unchanged. Two things to know before using it: nothing is shared across seasons except the pooled shore turnover derivation and the priors (each sub-season is an independent fit), and that derivation pools every I/E day in the workbook into one diel profile (D24), so check `shore_turnover_by_day.csv` by season before trusting a span whose seasons behave differently. (The `yday -> L` I/E regression pooled the same way; since B55 it runs only under a time-denominated shore unit.) **Pre-2024-25 seasons are out of scope (Matt, 2026-09-28).** 2022-23 and 2023-24 were exploratory: insufficient sampling coverage and not the current protocol (2023-24 has one peak count per day), and they are not to be run with this modelling approach. So the 2023-25 span still commented in `run_config.R` is not a run to make (D8, closed as out of scope; it was also blocked, with no 2023-24 vessel tally or charter roster). Multi-season spans remain supported mechanically from 2024-25 onward (e.g. 2024-26 once the 2025-26 data, including its OSP record, are complete).

**Pooled vs gear-resolved differ only in how CPUE is modeled.** Pooled uses one CPUE process and allocates gear-type catch after estimation by trip-level bootstrap interview shares per population and sub-season (review item 6, 2026-09-08; the Dirichlet it replaced treated each crab as independent). Gear-resolved is *written* for a per-gear CPUE process, but **as currently driven it runs with `G = 1`, so that machinery is inert** and gear-type catch is PE-apportioned. Do not raise `G > 1` without adding per-gear effort shares; only gear 1 is observed in the effort stream.

## `03_R_functions/` contract (load-bearing)

Every driver sources the **entire folder** wholesale in its setup chunk:
```r
purrr::walk(list.files(here("03_R_functions"), full.names = TRUE), source)
```
Because both production drivers source every file, two rules are structural, not stylistic:

1. **Pure functions, zero source-time side effects.** No top-level reads/writes/plots/`library()` calls, no assumption about sourcing order. Config is *passed* through `params`, never captured from a driver global. (`bss_timers.R` keeps state in a module-local environment precisely to avoid needing a global.)
2. **No name collisions.** Shared functions keep one name and one body. Functions that differ between tracks carry a `_pooled` / `_gear` suffix (`run_pe_pooled` vs `run_pe_gear`; `prep_bss_crab_pooled` vs `prep_bss_crab_gear`). **Adding a track-specific helper without a suffix silently overwrites the other track's function** depending on alphabetical source order. (`fetch_crab_data` was merged 2026-08-01 into one shared reader parameterized by `boat_require_gear_time`, retiring the former `_v2` variant; prefer this pattern, a shared function with a per-model param, over a name split when the two versions are nearly identical.)

Key shared helpers to know: `bss_model_adequacy.R` (model adequacy reported BESIDE the gate, never gating it), `bss_run_warnings.R` (`bss_warn()`: the run-level log of fallbacks, skips and frame conditions that both reports print in their closing "Run warnings and notes" section and write to `run_warnings.csv`; use it, not a bare `cat()`, for anything a reader must see, since most compute chunks are `results='hide'`), `bss_sampler_override.R`, `pe_gear_ratio_frame.R`, `annotate_decoupled_run.R`, `fetch_osp_boat_counts.R`, `diagnose_osp_trailer_overlap.R`, `bss_effort_spec.R` (single source of the effort unit AND the matching I/E observation column, read by both PE and BSS prep so they can't drift), `bss_convergence_gate.R` (the gate), `bss_ar_resolution.R` (adaptive AR selector + `bss_ar_ladder()`, the opt-in escalation ladder), `bss_day_length.R` (I/E ingest + `L_effective`), `crab_fraction.R` (the crabbing fraction `f`, its per-stratum priors, and the OSP crab-only lower bound), `bss_opener_covariates.R` (opener effort covariates and their multiplicity-controlled screen), `bss_marine_hazard_covariates.R` (the NWS Small Craft Advisory flag as an effort covariate on the same block: the boat all-gear term is the method of record since 2026-09-27, A30; the bar-restriction tick and the season split are offered for experiments; `marine_hazard_terms_for()` confines a selected term to the sub-seasons in `marine_hazard_gear_regimes`), `diagnose_incomplete_trips.R` (the four-arm treatment diagnostic), `build_subseasons.R`, `prep_days_crab.R`, `prep_population_summary.R`.

## Cross-cutting invariants

- **Paths: everything resolves through `here::here()`**, anchored to the repo root via the `.Rproj` / `.git` sentinels, *not* the `.Rmd` location, so a driver knits correctly from any working directory. But the directory **names** inside `here(...)` must match the on-disk folder names exactly; renaming a numbered stage folder requires updating every `here()` string in the drivers and diagnostic writers.
- **Effort unit is gear-deployments** for both shore and boat (as of v7.7 / v7.6). Time-denominated units (crabber-hours, gear-hours) are invalid for pot/trap gear because catch is sub-linear in soak time; the pipeline re-measures this every run (`cpue_linearity_*` / `cpue_saturation_*` CSVs). Any new gear type must pass those before its totals are trusted. A "gear-deployment" is a piece of gear a crabber had in the water on that trip (`number_of_gear`); **not** a pot lift; repeat checks of one gear slot are not extra deployments, and slot re-use across the day is carried separately by `tau`.
- **The daily expansion factor `L` is a TURNOVER, not a day length.** `L = tau_shore` / `tau_boat`. **Both centres are DERIVED FROM DATA, not set, and both moved in 2026-09:** `tau_shore_prior_mu = "derived"` reads the I/E `time` column and gives **2.477** on 2024-25 (the retired literal 1.7 was arrivals over PEAK presence, which is the wrong quantity); `tau_boat_prior_mu = "calibration"` reads the OSP/trailer overlap, and with `shared_tau = TRUE` the fitted `tau_bar` is **3.120** in the posterior mean (median 3.102 [2.49, 3.86]) on the authoritative run (3.118 on the 2026-09-28 run; 2.977 on R4, before the boat advisory-day term). The literals `1.7` and `1.2` survive in the code only as `%||%` fallbacks and in the ladder's R1 rung, which exists to price exactly this change; do not read either as a current value. Open item D24: the derived shore centre pools ALL 40 I/E days, and this window's own 6 give 2.225, 10.2% lower, worth about 3,500 crab. `L_effective` in hours (~5.3), the civil-twilight day length and the I/E day-length regression are computed ONLY when `shore_effort_unit` is a time unit (or `day_length_diagnostics = TRUE`), since 2026-09-29 (B55): under gear-deployments no estimate reads them, so they are not computed, the day-length columns are NA and `L_effective_ie_detail.csv` is not written. Anything that consumes `L`; including the I/E observation stream; must carry the matching unit. `bss_effort_spec()` is the single source for that pairing (`h_col`, `L_data`, and as of 2026-08-25 `ie_obs_col`); if you add a consumer of effort, route it through that function rather than re-deriving the formula. Four separate places had drifted onto the pre-v7.7 crabber-hours formula before the 2026-08-25 audit, and one of them (the shore I/E observation) was materially wrong.
- **The PE prices its unsampled cells, and the choice is a real lever.** The PE builds (week x day-type) strata from sampled days and left-joins the calendar, so a cell with calendar days but no sampled day used to contribute ZERO effort: 45 of 289 shore all-gear days and 48 of 289 boat days, concentrated in the high months. `pe_effort_strata.R` is the single shared home of the stratification, the fill and the variance, called by both `run_pe_pooled()` and `run_pe_gear()`, so the two tracks cannot drift. Shipped: `pe_empty_effort_stratum = "local_day_type"` (borrow the month-local day-type mean), `pe_empty_stratum = "local"` (the CPUE fill), `pe_variance = "impute_aware"`. **Mean and spread are borrowed from DIFFERENT donor levels on purpose** (the finest level with any sampled day for the mean, the finest with at least two for the spread); coupling them shifts the point estimate, which is how the fill was first got wrong. The effort fill is vindicated (shore PE within 1.3% of shore BSS, against -21.5% at `zero`); the CPUE fill is NOT settled and the evidence runs against the shipped `"local"` (D19), so treat it as open.
- **AR temporal resolution is an inference lever, not a tuning knob.** It's selected per fit from data density (daily/weekly/monthly) then capped per population via `ar_max_resolution`; too-fine AR is unidentified (the pooled boat's daily AR once diverged ~100%). As of 2026-08-25 an opt-in escalation ladder (`ar_escalate`, default FALSE) can instead start at the finest rung and coarsen only on a convergence-gate failure; it costs one extra multi-hour fit per failed rung, which is why it ships off.

## Input data quirks (real, not bugs, do not "fix")

From `04_input_files/` (see its README).

**SIX OF THE NINE MODEL AND DIAGNOSTIC WORKBOOKS ARE BUILT, NOT AUTHORED. Never hand-edit one.** Since the 2026-09-10 rebuild, `interview_combined.xlsx`, `effort_combined.xlsx`, `sampler_shifts.xlsx`, `wes_commercial_tally.xlsx`, `charter_trips.xlsx` and `crabbing_holidays.xlsx` are generated from the per-season creel workbooks in `04_input_files/raw/<YYYY><YY>_rec_crab_harvest_data.xlsx` by `04_input_files/build_all_inputs.R` (which runs the six `build_*.R` scripts in dependency order); a hand edit is silently discarded the next time anyone adds a season. Fix the raw workbook or the builder, then rebuild and commit the raw file and the rebuilt inputs together. The three inputs that are NOT built, because their sources are not the season workbooks, are `ingress_egress.xlsx` (the I/E database export, hand-maintained), `WBL_boat_counts.xlsx` (OSP) and `fishery_opener_dates.xlsx` (the regulation calendar). A tenth workbook, `nws_marine_hazards.xlsx` (2026-09-25), is built by its OWN builder from the NWS VTEC archive (`build_nws_marine_hazards.R`, network), is not in `build_all_inputs.R`, and since 2026-09-27 is read by EVERY production run (the boat all-gear advisory term is the method of record); the committed copy is a live pull (2026-09-25; the transcription it was verified against is in `04_input_files/raw/`). Rebuild it before a new season's run or the run stops. An eleventh, `osp_sampling_rates.xlsx` (2026-09-28, B48), is OSP's minimum sampling-rate schedule from its sampling manual, hand-maintained; it is read only when `WBL_boat_counts.xlsx` carries a crab-only column, to turn the day's count into the number of boats OSP SAMPLED (the binomial n of the crab-only share; `osp_sampling_rates.R`). OSP samples every k-th private boat at a rate set for the day, so the crab-only count must never be divided by the day's total.

As of the 2026-07-16 input migration every model input is an **`.xlsx` workbook with a single `data` sheet** (previously a mix of CSVs), read via `readxl::read_excel`, with dates stored as ISO `yyyy-mm-dd` text parsed by `as.Date()`. Quirks that survive re-export and are matched in code:
- **Float 17-21 is counted only when a second sampler is free** (Matt, 2026-09-28): one sampler covers Float 20 only; with two, the second works the boat launch and counts Floats 17-21 when they can. A Float 20 count with no Float 17-21 count beside it is therefore UNSAMPLED at 17-21, not empty. `shore_dock_counts()` fills it by the month's time-paired ratio (`shore_f17_fill = "ratio"`, B45); until 2026-09-28 it was 0, which read the 2024-25 shore gear count about 9% low. Do not "fix" it back to a zero.
- Interview `number_of_gear` maps from **column N, not W** (duplicate iForm field name).
- The commercial `boat_type` is the typo **"Commerical"** (one m), matched by regex, don't correct the spelling without updating the matcher.
- Windows/OneDrive long paths can exceed MAX_PATH; the code detects this and falls back to a short path.
- Historical note (pre-migration, no longer applies): the CSVs needed `QUOTE_ALL` for the notes field, and carried **M/D/YYYY** dates that silently parsed to `NA` under the default reader; the xlsx conversion fixed both, and the notes / unused columns were dropped in the same modeling-only strip.

`ingress_egress.xlsx` is read by **both** drivers (via `fetch_ie_data` in `bss_day_length.R`, called from the pooled and gear-resolved drivers alike), feeding the shore turnover derivation (`shore_turnover_*.csv`) and the shore I/E observation stream (and the `L_effective` day-length model only under a time unit, B55); the earlier "pooled driver only" claim was wrong (input-audit F1). The one genuinely pooled-only input is the opener calendar `fishery_opener_dates.xlsx`, read solely by the pooled report's spillover diagnostic, so it changes no estimate.

## Outputs

Each run writes to `05_output/<YYYYMMDD>/<model>-<run_tag>/`, where the `<model>` prefix is `pooled-CPUE` or `gear-type-CPUE-model` (the `pooled-CPUE-covariates` prefix belonged to the weather-tide module removed 2026-09-13, A29; no code writes it now) and `run_tag` (from `run_config.R`, or `run_rg_sweep.R`) names the folder, e.g. `pooled-CPUE-run1`; a blank `run_tag` falls back to an `HHMMSS` timestamp, so same-day runs of a model no longer overwrite. Per-population files follow `<metric>_<population>_<species>_<fate>.{csv,png}`; port/monthly/comparison files drop the population tag. **Outputs are committed to git on purpose** so past estimates are preserved as-produced; only per-run `*.RData` workspaces and compiled `*.rds` Stan caches are git-ignored. A given dated folder may be a *partial* run; use a recent complete run as the reference catalog, not any single folder.

## Working conventions

- **Validate by run, never by reasoning alone.** A change that looks inference-neutral on paper can perturb the sampler geometry. Isolate each change, compare against a confirmed baseline run against pre-set criteria; an item is not "done" until a run confirms it. Changes that *narrow* reported uncertainty (e.g. tightening dispersion priors) need explicit sign-off because they move the headline intervals.
- **Diagnostics are additive and `tryCatch`-wrapped** so one failing fit can't abort a run.
- **DATING A MARKER YOU WRITE.** A dated label in this repository (a section title, a `[NEW ...]`
  bracket, a register row, a `# 2026-mm-dd:` comment) names the **working session**, not the
  commit. The two are not the same: the sessions that produced this branch ran on a clock that
  was ahead of the repository's by up to five days before 2026-09-08, and the drift is baked
  into filenames (`run_adoption_2026-09-07.R` was committed 2026-09-04), so it cannot be
  normalized away. Two rules follow. **(a) When you write a marker into
  `PIPELINE_STATUS.md`, `CHANGE_REGISTER.md` or `VALIDATION_CAMPAIGN.md`, use the git author
  date of the commit that lands the work, and give the hash when it differs from the session
  date.** Those three documents exist to be a chronology and harness section 70 enforces it:
  no marker in them may postdate the document's own `Last updated` line or the branch tip, the
  campaign's section letters may not run backwards except where declared, and every section
  must appear in the git-anchor table. **(b) Everywhere else, the session date is fine and is
  what the convention means**; do not start a repository-wide re-dating. The measurement, the
  evidence and the reconciliation live in one place: "A NOTE ON THE DATES" at the top of
  `VALIDATION_CAMPAIGN.md`.
- **FIVE DOCUMENTS, FIVE JOBS. Know which one you are in before you write to it** (the layout settled 2026-09-12):

| document | job | when to write to it |
|---|---|---|
| `07_documentation/BSS-GH-pooled-CPUE-model-documentation.md` | **Method v2.0**: what the model IS. Specification, and the limitations | an adopted change moves the method |
| `07_documentation/development_notes/PIPELINE_STATUS.md` | the CURRENT state and the backlog. **The authoritative run and its total are in the box at the top** | state changes, or an item opens or closes |
| `07_documentation/development_notes/CHANGE_REGISTER.md` | the tabular register: every change, its status, its evidence, its effect on the number; every defect and what it cost | every change, and every defect, without exception |
| `07_documentation/development_notes/VALIDATION_CAMPAIGN.md` | HOW it got there: the dated run-by-run narrative, Sections 1b to 1z (1b to 1v were the status document's until 2026-09-12). Every "Section 1x" reference resolves here | a batch runs and is reviewed |
| the two `*-development-history.md` files | the backward-looking, newest-first version logs | a version lands |

  **Do not take an estimate from anywhere in this repository without checking the box at the top of `PIPELINE_STATUS.md` first**, including from the method document, whose numbers come from a named reference run that may have been superseded. Read the status document first for current state and update it in place; do not fork per-session notes.
- **Method version ≠ code version, and since 2026-09-12 the method version is LIVE.** **Method v2.0** is the method of record and it tracks the working model: when an adopted change moves the method, `BSS-GH-pooled-CPUE-model-documentation.md` moves with it. That is a deliberate break with **Method v1.0**, which was frozen against pooled code v7.4 and which the code then outran in nine separate ways; v1.0 is archived at `07_documentation/archive/method-v1.0-pooled-CPUE.md` with the difference table, because several current decisions are only legible as departures from it. The code version (v7.x, and framework v6.0 on the gear track) still moves independently and faster. Numbers predating the 2026-07-11 refresh are labeled "pre-refresh"; do not cite them as final.
- `06_diagnostics/` holds **the validation harness and the batch runners**, and since 2026-09-13 nothing else (the weather-tide module was removed, A29): `test_improvements_2026-08-25.R` is a dependency-light regression harness (about **1,488** assertions as of 2026-10-02 and growing, a couple of minutes, no rstan; the run prints the count, so take that over this figure) that should be run before committing to any long fit, and `run_*.R` are the dated experiment batches. Every batch runner that HAS a `DRY_RUN` switch ships it `TRUE`; that is asserted by the harness. Four scripts have none: `run_osp_validation.R` and `run_tau_sweep.R` refuse to fit (superseded; B31, and B46 for the tau sweep; every guarded runner fails closed if the guard file is missing), and `desk_sca_season_split_2026-09-27.R` and `gear_coverage_audit.R` fit nothing. (`run_rg_sweep.R` gained a `DRY_RUN` switch, shipped `TRUE`, and its three-rung grid on 2026-09-28, B49.) `run_authoritative_batch_2026-09-29.R` (B58) ships `DRY_RUN <- TRUE` and was STARTED with `--go` (or `BSS_BATCH_GO=1`), never by editing the file, because its stage A was the authoritative render; it finished on 2026-10-02 and, like `run_gear_ar_zi_2026-09-13.R` (whose question, D3 and D6, it settled), now refuses to fit (B62): its stages are deltas from `run_config.R` as shipped at `c4a0241`, which A32 moved. (An earlier version of this bullet still described the weather-tide module as living in this folder; it was removed on 2026-09-13, and the sentence went stale with it. Corrected 2026-09-25.)
