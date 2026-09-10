# Running a new season (or any new window)

**Audience:** the analyst pointing this pipeline at data it has never seen.
**Framing:** the 2024-25 season was the DEVELOPMENT TEST SEASON. Every cap, floor, prior center and sampler setting in the shipped config was derived on it. The architecture is built to run on any window you select: a full season, part of one, or a multi-season span. This guide is the workflow that takes you from a naive first run to a defensible estimate, using the diagnostics the pipeline writes to make each tuning decision a measured one instead of a guess.

The one-sentence version: **configure the season, run naively with the AR ladder on, read what the ladder and the adequacy files tell you, pin each fit's resolution deliberately, then produce.**

---

## 0. What a season needs before anything runs

All inputs live in `04_input_files/` and carry a `season` column, so multiple seasons coexist in one workbook and a new season is rows added, not files replaced. Since 2026-09-10 seven of the nine workbooks are BUILT from the per-season creel workbooks: **drop the season's creel workbook into `04_input_files/raw/` as `<YYYY><YY>_rec_crab_harvest_data.xlsx` (e.g. `2627_...`) and run `Rscript 04_input_files/build_all_inputs.R` from the repository root**; read each builder's report (rows by season, the flagged rows, the comparison with the previous workbook) before committing the workbook and the rebuilt inputs together.

| Input | New-season action |
|---|---|
| `interview_combined.xlsx` | **Built** by `build_interview_combined.R` from every season workbook's `Harvest Creel` sheet. The `Trip Type` labels must classify (crab only / `<fishery> & Crab` / finfish only / non fishing) and the gear labels must be in the vocabulary map; an unrecognised value is printed. A row dated outside its workbook's fishery season stops the build. |
| `effort_combined.xlsx` | **Built** by `build_effort_combined.R` from the `Effort count` sheets. Check the `qc_flag` table it prints: a count recorded from interviews rather than observed must carry the note the builder recognises ("from interviews", "based on interviews", "total from interviews"), or add the pattern. |
| `sampler_shifts.xlsx` | **Built** by `build_sampler_shifts.R` from the `crab creel survey data` sheets; flagged rows (a check-out typed on a 12-hour clock, a missing check-out, a duplicated survey id) are listed and held out of the hours. |
| `wes_commercial_tally.xlsx` | **Built** by `build_comm_charter_tally.R` from the `wes commercial tally` sheet and the survey sheet's vessel-tally columns (reconciled per day; read by header name because the column order has differed between seasons). Then set the season's `census_windows` entry: pots legal (Dec 1) to the day before the commercial opener. |
| `charter_trips.xlsx` | **Built** by `build_charter_trips.R` from the `charter trips` sheet (interviewed / missed / canceled). Under `charter_frame = "roster"` it is the frame the charter component is EXPANDED over (the charter vessels are not fully sampled), so the roster's completeness sets the charter estimate: check that every operator's trips are on the sheet. A season without the sheet falls back to the tally frame, which counts no charter trip on a day without a tally, and the run says so. |
| `crabbing_holidays.xlsx` | **Built** by `build_crabbing_holidays.R` by rule (the named holidays, including Thanksgiving Day, Veterans Day and Juneteenth since 2026-09-11, each with its federal observed day; Super Bowl dates are a table in the builder, extend it). **Every season in `season_filter` needs rows**: a missing season STOPS the run by design. The builder prints the sampler-flagged dates the rule does not cover, which is where the next candidate would come from. |
| `ingress_egress.xlsx` | Hand-maintained from the I/E database export (put the export in `raw/`; keep the `time` column). **Zero I/E rows is survivable**: `L_effective` falls back down the civil-twilight ladder, and the run says so. |
| `WBL_boat_counts.xlsx` | Add OSP daily boat totals if available. **Absent dates are non-sampled, not zero**; with no OSP rows in the window the stream is simply absent and the boat runs trailer-only, and the shared boat turnover is refused for want of overlap days (2025-26 as of 2026-09-10). |
| `fishery_opener_dates.xlsx` | Add the season's opener rows if you use the opener diagnostics (production covariate mode is off). |

## 1. The per-season config checklist

Everything is in `run_config.R`; the keys below must move TOGETHER. A checklist version lives at the top of that file.

1. **The window:** `est_date_start`, `est_date_end`. Any span; it does not need to be a whole season.
2. **The data filter:** `season_filter`. A stale value used to produce a wall of empty fits with no explanation; it now stops the run with a plain message (`validate_season_window.R`), and every run prints what the season/window selection actually captured. A character vector selects a multi-season span.
3. **The closure calendar:** `pot_closure_start`, `pot_closure_end`, `pot_open_date`. A window that does not intersect the closure yields a single all-gear sub-season (fixed 2026-09-09); a window inside the closure yields a single pot-closure sub-season; a mid-window closure yields three. A span containing several closures uses `pot_closures` instead, one entry per season (added 2026-09-10; see section 7). **Since 2026-09-12 `validate_season_window()` warns when any of these three, or the two census dates, falls outside the estimation window** -- the failure mode is a single-season ROLLBACK from a multi-season span, where these are exactly the keys left behind. `pot_open_date` is the least dangerous of them (its only consumers are the "Pots open" plot marker and a fallback for a NULL `pot_closure_end`; the "L_effective I/E regression split" it was once documented as feeding never existed) and `pot_closure_start`/`end` are the most, because `build_subseasons()` splits the season on them.
4. **The census window:** `census_start_date`, `census_end_date` (the commercial/charter tally span; independent of the estimation window and easy to forget).
5. **The AR caps:** `ar_max_resolution`. These are 2024-25 answers and the subject of section 3. Treat them as starting points.
6. **Every key tagged `SEASON-DERIVED`** in `run_config.R`: the `kappa_OSP` prior center (3.0, from the 2024-25 overlap days), the ZI prior shape (Beta(1,9), from the 2024-25 zero bin), `shared_tau_min_obs` (15; check the printed OSP-informed-day count against it), `tau_boat_prior_mu_fallback` (2.7) and `tau_shore_prior_mu_fallback` (1.7). Since 2026-09-08/09 the two turnover prior centres resolve per run (`tau_boat_prior_mu = "calibration"` from the window's OSP/trailer overlap, `tau_shore_prior_mu = "derived"` from the window's I/E time column and count hours), the crabbing fraction f and the combo share c are read from the season's contacts and trip types, and the census is the exact sum over the tally days; the fallbacks apply only when a season lacks the data behind them, and the run says so.
7. **The PE's unsampled-cell levers:** `pe_empty_effort_stratum`, `pe_empty_stratum`, `pe_variance`. Not season-derived, but read the per-component numbers in `pe_empty_effort_strata.csv` for the new season before citing anything: with weekly strata and ~50% day coverage, roughly half of a component's calendar days rest on one sampled day or none, and a thinner season makes it worse. These levers change the PE point and its SE, and a component whose gate FAILS reports its PE point in the port total as a constant, so on a thin season they can reach the headline.
8. **`run_tag`:** name the run so you can find the folder.

Worked values for **2025-26** (data through 2026-09-08; the season ends 2026-09-15), in the commented block under `season_filter` in `run_config.R`: window 2025-09-16 to 2026-09-15, `season_filter = "2025-26"`, closure 2025-09-16 to 2025-11-30 with pots legal 2025-12-01, census window 2025-12-01 to 2026-01-03 (Grays Harbor's commercial fishery opened Jan 4, 2026), `commercial_opener = "2026-01-04"`. As of 2026-09-10 the season has no OSP rows and no WBL boat I/E days, so the boat runs on the trailer counts with per-day turnover draws (a short smoke fit of the boat all-gear component sampled in 4 minutes with 37 of 2,000 divergences, the known funnel of the per-day-tau geometry; expect to pin the boat on the OSP file when it arrives).

Also know that the drivers' own `params_model` blocks carry per-fit sampler settings (iterations, treedepth, adapt_delta keyed by fit name) tuned on 2024-25 geometry. They are conservative, so they rarely need touching, but a new season's pathological fit is tuned there, not in `run_config` (`bss_sampler_override` is the sanctioned route from the config side).

## 2. The naive first run

Run the harness first; it is seconds and pins the shipped invariants:

```r
Rscript 06_diagnostics/test_improvements_2026-08-25.R
```

Then configure the ladder for discovery and run the pooled driver once:

```r
ar_escalate             <- TRUE        # every fit climbs the ladder
ar_escalate_stop        <- "first_pass"
ar_rung_adequacy        <- TRUE        # per-rung p_loo / Pareto k / coverage (the decisive columns)
ar_escalate_respect_cap <- FALSE       # ignore the 2024-25 caps; that is the point
```

**Cost, stated plainly:** every rung is a full MCMC fit, and coarsening barely reduces per-fit cost (the 2026-09-04 ladder's rungs cost 96/91/83 minutes each). Budget one full fit per rung per component. Scope the ladder (`ar_escalate = list(shore = "all_gear")`, or a character vector of populations) when you already trust some components.

**What "first_pass" answers:** the ladder starts each fit at daily and coarsens until the convergence gate passes, so the reported rung is the FINEST resolution at which that fit samples, which is exactly the "minimum AR necessary to pass" question a naive season asks. `ar_escalate_stop = "all_rungs"` instead fits every rung and is the full-diagnostic mode, at full cost.

## 3. Reading the run: the decision files, in order

1. **The console/season check** (printed at read time): did the season/window selection capture the data you think it did?
2. **`season_summary.csv` and `convergence_report.csv`:** which fits attempted, which passed, which fell back to PE and why. A PE fallback on a thin component is a measured outcome, not a failure.
3. **`ar_escalation_log.csv`:** one row per rung per fit: resolution, `P_n`, divergences, the gate verdict, each rung's own catch estimate and interval, and (with `ar_rung_adequacy`) per-rung `p_loo`, Pareto k count and coverage. The `selected` flag marks the reported rung.
4. **`model_adequacy.csv`:** the reported fit's adequacy beside the gate. **The gate answers "did this fit sample"; it never answers "is this model right"**, and on a dense season every rung can pass the gate, in which case adequacy is the only thing that separates them.
5. **`pe_vs_bss_comparison.csv`:** the design-based cross-check per component.

**How to choose a rung when several pass** (the 2024-25 ladder, labeled as the worked example, not the answer): all four rungs passed the gate and spanned only 3.7% in catch, so the gate and the estimate decided nothing. What decided it: at daily, `p_loo` was 35% of `n_obs` (one effective parameter per three observations), 41 Pareto k above 0.7 (so its own elpd was unreliable), and effort coverage was +7.1 sampling SD; at weekly all three were clean. Two cautions that generalize: **a narrower prediction interval from a coarser latent process is not precision you earned** (the miscalibrated cells in the 2026-08 boat 2x2 had the narrowest relative intervals), and **an elpd advantage means nothing when the fit's Pareto k have failed**. Prefer the finest rung whose adequacy is clean; use interval width only to break ties between rungs that are ALREADY adequate.

## 4. Pinning the resolutions

Two levers, used in sequence:

- **While deciding:** `ar_force` pins a fit to an exact rung, finer or coarser than the selector would pick (caps can only coarsen, so a pin FINER than the data-driven choice is only expressible here). Scope it per population x sub-season: `ar_force = list(shore = list(all_gear = "weekly"))`. A per-population entry reaches BOTH of that population's sub-seasons; that scoping mistake once moved a component 3,025 crab.
- **Once settled:** move the answer into `ar_max_resolution` and return `ar_force` to NULL, so production reads from the cap and the config self-documents. The two routes were proven to build byte-identical Stan data (2026-09-07), and the adoption pattern to copy is `run_adoption_2026-09-07.R`: render once with the new caps and require bit-identity with the pinned run.

## 5. The production run and the cross-check

With the ladder off and the caps set, `source("run_estimation.R")` is the production run. Then run the gear-resolved track once (about half an hour) and read the two tracks together, with the lesson the 2026-09-07 cross-check taught: **compare the tracks at the same resolution before calling a gap a disagreement.** The two tracks agreed on the shore component to 0.08% at a common resolution while sitting 3.4% apart as configured, because the resolution difference dominates. Also note the gear track's shore catch is plain NB2 (no zero-inflation block), a small (~0.3%) structural asymmetry.

## 6. Part-season windows

Supported. What changes:

- The sub-season list adapts (section 1.3): one sub-season when the window misses or sits inside the closure; `all_gear_pre`/`all_gear_post` names when a closure falls mid-window.
- Components with no data in the window fall to their floors (`bss_min_interviews`, post-filter floor) and report PE or empty, with the reason in `convergence_report.csv`.
- The census component is empty unless the census window overlaps yours (the validator warns).
- Expect coarser feasible AR: fewer days means fewer periods, and the ladder's degenerate-rung dropping will shorten itself automatically.

## 7. Multi-season spans

Supported as of 2026-09-10 (CHANGE_REGISTER A14). A span takes four settings that must agree:

1. `est_date_start` / `est_date_end` spanning the whole range, and `season_filter` as a vector, e.g. `c("2023-24", "2024-25")`.
2. `pot_closures`: a list with one `list(season =, start =, end =)` entry per closure in the span. This outranks the scalar `pot_closure_start/end` pair and yields one pot-closure sub-season per season (`ring_net_only_<season>`, biweekly, pots excluded) with the following all-gear block as `all_gear_<season>`; a window starting before the first closure gets `all_gear_pre_<season1>`. Overlapping closures, an unlabeled entry, or end-before-start stop loudly. With zero or one in-window closure the legacy scalar path runs byte-for-byte, so single-season names, fit labels, and filenames never change.
3. `census_windows`: a NAMED list, season to `c(start, end)`, for the commercial/charter census; the census is expanded per window and summed, and each season's own component is kept for the season table. An unnamed list stops.
4. The data. Every season in the span needs rows in EVERY workbook of section 0: effort counts, interviews, ingress/egress, holidays, opener dates, and the tally. The validator prints per-season capture and warns hard on the asymmetric case (interviews present, zero effort counts: that season's effort would be pure imputation, not estimation); the holiday reader stops on a missing season.

The report then writes `season_totals.csv` (always) and renders a season-summary table (when the span has more than one season). Monthly figures and fit diagnostics already cover the full span: month labels are year-qualified (`%Y-%m`) and the calendar indices are span-safe (sequential year-week and year-month factors; no aliasing).

Two things to know about a span. `pot_open_date` is still a single scalar feeding the ingress/egress `L_effective` regression split, so the split is exact only for the first season (tracked in A14). And a span is NOT one statistical process: `build_subseasons()` emits one sub-season per closure and the driver fits each one independently (its own effort, CPUE and turnover parameters), so nothing is shared across seasons except the I/E day-length regression (pooled across seasons unless `ie_filter_by_season = TRUE`) and the config priors. A season's total from a span therefore equals a standalone single-season run of the same data up to those two inputs and Monte Carlo error; the value of a span is one run, one report and one set of figures, not borrowed strength between seasons. (Corrected 2026-09-08; the first version of this section claimed a shared process.)

## 8. Failure modes, and what each one means

| Symptom | Meaning | Fix |
|---|---|---|
| Run stops at "estimation window contains NO effort counts and NO interviews" | `season_filter` does not match the season rows in the window | Update `season_filter` with the window |
| Run stops at "No crabbing holidays for season(s)" | The holiday calendar has no rows for a requested season | Add rows to `crabbing_holidays.xlsx` |
| "OSP boat-count days in window: 0" | No OSP coverage; stream absent, boat runs trailer-only, `f` lower bound inert | Expected when OSP did not operate; nothing to fix |
| "I/E observations: 0" | `L_effective` regression has no data | Civil-twilight fallback engages; check the day-length plot |
| A fit reports "PE fallback" with post-filter interviews below the floor | Too thin to fit | A measured outcome; the PE carries the component |
| Every ladder rung fails the gate for one fit | The component cannot support ANY AR on this data | PE fallback is the answer this season; collect more data |
| Gate passes but `flag_miscalibrated` is TRUE / `p_loo` is a large fraction of `n_obs` | Sampled fine, model too flexible at this resolution | Coarsen via the ladder evidence (section 3) |

## 9. What this guide does not cover

Site generality (the readers filter to Westport / Grays Harbor locations via `gh_effort_areas` / `gh_creel_location`; a new PORT is a data-mapping exercise, not a season), the frozen Method v1.0 documents (historical; see their banners), and the open modelling items, which live with their evidence in `CHANGE_REGISTER.md`.
