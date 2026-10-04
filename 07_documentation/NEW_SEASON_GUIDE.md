# Running a new season (or any new window)

**Audience:** the analyst pointing this pipeline at data it has never seen.
**Scope of seasons (Matt, 2026-09-28):** this workflow is for 2024-25 and later seasons. The seasons before 2024-25 (2022-23, 2023-24) were exploratory, with insufficient sampling coverage and not the current protocol (2023-24 has one peak count per day), and they are **not to be run with this modelling approach**, alone or in a span (CHANGE_REGISTER D8). Their rows stay in the workbooks as a record.
**Framing:** the 2024-25 season was the DEVELOPMENT TEST SEASON. Every cap, floor, prior center and sampler setting in the shipped config was derived on it. The architecture is built to run on any window you select: a full season, part of one, or a multi-season span. This guide is the workflow that takes you from a naive first run to a defensible estimate, using the diagnostics the pipeline writes to make each tuning decision a measured one instead of a guess.

The one-sentence version: **configure the season, run naively with the AR ladder on, read what the ladder and the adequacy files tell you, pin each fit's resolution deliberately, then produce.**

---

## 0. What a season needs before anything runs

All inputs live in `04_input_files/`, and all but three carry a `season` column, so multiple seasons coexist in one workbook and a new season is rows added, not files replaced. The three exceptions are `WBL_boat_counts.xlsx` (OSP's own layout, below), `nws_marine_hazards.xlsx` (a calendar of hazard events) and `osp_sampling_rates.xlsx` (OSP's sampling schedule, not per season). Since 2026-09-10 six of the nine model and diagnostic workbooks are BUILT from the per-season creel workbooks (the tenth, the NWS marine hazard archive, has its own builder, needs the network, and since 2026-09-27 is read by every production run; last row of the table): **drop the season's creel workbook into `04_input_files/raw/` as `<YYYY><YY>_rec_crab_harvest_data.xlsx` (e.g. `2627_...`) and run `Rscript 04_input_files/build_all_inputs.R` from the repository root**; read each builder's report (rows by season, the flagged rows, the comparison with the previous workbook) before committing the workbook and the rebuilt inputs together.

| Input | New-season action |
|---|---|
| `interview_combined.xlsx` | **Built** by `build_interview_combined.R` from every season workbook's `Harvest Creel` sheet. The `Trip Type` labels must classify (crab only / `<fishery> & Crab` / finfish only / non fishing) and the gear labels must be in the vocabulary map; an unrecognised value is printed. A row dated outside its workbook's fishery season stops the build. |
| `effort_combined.xlsx` | **Built** by `build_effort_combined.R` from the `Effort count` sheets. Check the `qc_flag` table it prints: a count recorded from interviews rather than observed must carry the note the builder recognises ("from interviews", "based on interviews", "total from interviews"), or add the pattern. |
| `sampler_shifts.xlsx` | **Built** by `build_sampler_shifts.R` from the `crab creel survey data` sheets; flagged rows (a check-out typed on a 12-hour clock, a missing check-out, a duplicated survey id) are listed and held out of the hours. |
| `wes_commercial_tally.xlsx` | **Built** by `build_comm_charter_tally.R` from the `wes commercial tally` sheet and the survey sheet's vessel-tally columns (reconciled per day; read by header name because the column order has differed between seasons). Then set `census_start_date` / `census_end_date` (or, in a span, that season's `census_windows` entry): pots legal (Dec 1) to the day before the commercial opener. |
| `charter_trips.xlsx` | **Built** by `build_charter_trips.R` from the `charter trips` sheet (interviewed / missed / canceled). Under `charter_frame = "roster"` it is the frame the charter component is EXPANDED over (the charter vessels are not fully sampled), so the roster's completeness sets the charter estimate: check that every operator's trips are on the sheet. A season without the sheet falls back to the tally frame, which counts no charter trip on a day without a tally, and the run says so. |
| `crabbing_holidays.xlsx` | **Built** by `build_crabbing_holidays.R` by rule (the named holidays, including Thanksgiving Day, Veterans Day and Juneteenth since 2026-09-11, each with its federal observed day; Super Bowl dates are a table in the builder, extend it). The rule covers the seasons 2022-23 through 2026-27 (`map_dfr(2022:2026, season_rows)` in the builder): for a later season extend that year range and rebuild; **never hand-edit the workbook** (Matt updates the holiday input for new seasons as needed). **Every season in `season_filter` needs rows**: a missing season STOPS the run by design. The builder prints the sampler-flagged dates the rule does not cover, which is where the next candidate would come from. |
| `ingress_egress.xlsx` | Hand-maintained from the I/E database export (put the export in `raw/`; keep the `time` column). **Zero I/E rows is survivable**: the shore turnover derivation falls back to `tau_shore_prior_mu_fallback` (2.477, the 2024-25 derivation) with log-SD 0.3, and the run records it in the report's warnings section and `run_warnings.csv` (D36). (`L_effective` in hours and the civil-twilight day length are computed only under a time-denominated `shore_effort_unit` since 2026-09-29, B55.) |
| `WBL_boat_counts.xlsx` | Hand-maintained from OSP's deliveries, in OSP's own layout, NOT the `data`/`season` layout of the other workbooks: one sheet, `Sheet1` (`osp_boat_counts_sheet`), with integer columns `Year`, `Month`, `Day` and `WestportPrivateEffort` (the day's private-boat total), plus the crab-only column (`WPTPrivateCrabOnly`, as OSP delivered it on 2026-10-02 for 2024-03 to 2025-10, from its `crabbing trips` sheet; `crabbing_only` and `WestportCrabOnlyEffort` are also read, first present wins; never carry OSP's `WPTPrivateCrabAlso`, which comes from the form's free-text notes and was not consistently sampled) and the day's sampling rate or boats sampled (`WestportPrivateSampleRate` / `WestportCrabClassified`, B48). To add a season, append one row per OSP day in those columns; no season column is needed, since the reader builds the date and filters to the window. **Absent dates are non-sampled, not zero**; an observed 0 is data. With no OSP rows in the window the stream is simply absent and the boat runs trailer-only, and the shared boat turnover is refused for want of overlap days. **2025-26 as of 2026-09-28: PARTIAL**, 23 OSP days, 17 Sep to 18 Oct 2025 (14 paired with Westport trailer counts); the rest is awaited from OSP. Corrected 2026-10-02: **the floor is applied PER FIT** (`prep_bss_crab_pooled()` counts the days that can inform L in that fit), and all 23 OSP days fall in the POT-CLOSURE sub-season (to 30 Nov). So the boat ALL-GEAR fit, about 47% of the port on 2024-25, has NO OSP day: no OSP effort stream, trailer counts only, and the shared turnover REFUSED (per-day draws around the calibration centre, the structure `shared_tau` was adopted to remove; since 2026-10-02 the refusal is written to `run_warnings.csv`). The calibration centre itself comes from 14 autumn pairs (3.50, against 2024-25's 3.03 from 61). The boat pot-closure fit has 2 interviews, below the floor of 15, and reports its PE. A citable 2025-26 boat figure needs OSP's March to September 2026 days (CHANGE_REGISTER D18). |
| `osp_sampling_rates.xlsx` | Nothing per season (2026-09-28, B48): OSP's minimum sampling-rate schedule by daily count, from its sampling manual, read only when `WBL_boat_counts.xlsx` carries a crab-only column. Update it only if OSP changes its manual. The day's actual rate (or the number of boats sampled), when OSP delivers it as a column, takes precedence over the schedule. |
| `fishery_opener_dates.xlsx` | Add the season's opener rows. The production covariate mode is off, but the spillover diagnostic runs by default (`run_fishery_spillover_diag = TRUE`) and stops if the workbook is absent; as of 2026-10-02 the workbook (last edited 2026-07-17) has no summer-2026 salmon or halibut openers, so complete it for 2025-26. |
| `nws_marine_hazards.xlsx` | **Required since 2026-09-27**: the method of record carries the NWS Small-Craft-Advisory flag on the boat all-gear effort process (`marine_hazard_mode = "manual"`, CHANGE_REGISTER A30 ADOPTED), so every production run READS this archive. It must COVER the window, or the run stops naming this step: `Rscript 04_input_files/build_nws_marine_hazards.R --end <the day after the window>` re-pulls both NWS zones from the IEM VTEC service (network; not part of `build_all_inputs.R`). Leave `--start` at its default, 2008-01-01: the committed archive starts there, and a later start (the builder's help text shows `--start 2023-01-01` as an example) would drop the earlier years that the four-season desk screen (`desk_sca_season_split_2026-09-27.R`, from 2022-23) reads. The bar-restriction tick needs nothing extra: it is read from `sampler_shifts.xlsx`, which the standard rebuild produces. |

## 1. The per-season config checklist

Everything is in `run_config.R`; the keys below must move TOGETHER. A checklist version lives at the top of that file.

1. **The window:** `est_date_start`, `est_date_end`. Any span; it does not need to be a whole season. Change `run_tag` with it: a season token in the tag that is not in `season_filter` stops the run (B46).
2. **The data filter:** `season_filter`. A stale value used to produce a wall of empty fits with no explanation; it now stops the run with a plain message (`validate_season_window.R`), and every run prints what the season/window selection actually captured. A character vector selects a multi-season span.
3. **The closure calendar:** `pot_closure_start`, `pot_closure_end`, `pot_open_date`. A window that does not intersect the closure yields a single all-gear sub-season (fixed 2026-09-09); a window inside the closure yields a single pot-closure sub-season; a mid-window closure yields three. A span containing several closures uses `pot_closures` instead, one entry per season (added 2026-09-10; see section 7). **Since 2026-09-12 `validate_season_window()` warns when any of these three, or the two census dates, falls outside the estimation window** -- the failure mode is a single-season ROLLBACK from a multi-season span, where these are exactly the keys left behind. `pot_open_date` is the least dangerous of them (its only consumers are the "Pots open" plot marker and a fallback for a NULL `pot_closure_end`; the "L_effective I/E regression split" it was once documented as feeding never existed) and `pot_closure_start`/`end` are the most, because `build_subseasons()` splits the season on them.
4. **The census window:** `census_start_date`, `census_end_date` (the commercial/charter tally span; independent of the estimation window and easy to forget), and `commercial_opener`, which `validate_season_window()` checks the census end against.
5. **The AR caps:** `ar_max_resolution`. These are 2024-25 answers and the subject of section 3. Treat them as starting points. **5b.** `gear_period_bss` (section 5): the gear track's periods; its shore entries mirror the pooled shore caps (A32), so change them together.
6. **Every key tagged `SEASON-DERIVED`** in `run_config.R`: the `kappa_OSP` prior center (3.0, from the 2024-25 overlap days), the ZI prior shape (Beta(1,9), from the 2024-25 zero bin), `shared_tau_min_obs` (15; check the printed OSP-informed-day count against it), `tau_boat_prior_mu_fallback` (3.03, the 2024-25 calibration on the shipped metric; 2.7 until 2026-09-29) and `tau_shore_prior_mu_fallback` (2.477, the 2024-25 derivation; 1.7 until 2026-09-29, D36). A fallback is used only when the per-run resolution below cannot run, and it is recorded in the report's warnings section. Since 2026-09-08/09 the two turnover prior centres resolve per run (`tau_boat_prior_mu = "calibration"` from the window's OSP/trailer overlap, `tau_shore_prior_mu = "derived"` from the I/E time column and the window's count hours -- **note, corrected 2026-09-11: the count HOURS are the window's but the diel presence PROFILE pools every I/E day in the workbook, so the derived centre is a multi-season quantity and is the same number whichever season you run. On 2024-25 only 6 of its 40 days are in the season, and the window's own 6 give 2.225 against the pooled 2.477, a gap of 0.83 SE that the existing I/E data cannot resolve either way. `tau_shore_derive_window_only = TRUE` restricts it; `shore_turnover_summary.csv` reports both either way. See CHANGE_REGISTER D24**), the crabbing fraction f and the combo share c are read from the season's contacts and trip types, and the census is the exact sum over the tally days; the fallbacks apply only when a season lacks the data behind them, and the run says so.
7. **The PE's unsampled-cell levers:** `pe_empty_effort_stratum`, `pe_empty_stratum`, `pe_variance`. Not season-derived, but read the per-component numbers in `pe_empty_effort_strata.csv` for the new season before citing anything: with weekly strata and ~50% day coverage, roughly half of a component's calendar days rest on one sampled day or none, and a thinner season makes it worse. These levers change the PE point and its SE, and a component whose gate FAILS reports its PE point in the port total as a constant, so on a thin season they can reach the headline.
8. **The NWS archive covers the window** (section 0, last row; the committed archive, pulled 2026-09-25, already covers 2025-26 to 2026-09-15, so a 2025-26 run needs no re-pull and a 2026-27 run does): `Rscript 04_input_files/build_nws_marine_hazards.R --end <the day after est_date_end>` (the default `--start`, 2008-01-01, keeps the whole archive) before the first run of a season, and commit the workbook. The flag definition itself (`marine_hazard_codes`, `_zones`, `_window`, `_tz`, section 2.10 of `run_config.R`) is part of the method and does not change per season; `marine_hazard_winter_months` (section 4.4b) matters only to the season-split experiment.
9. **`run_tag`:** name the run so you can find the folder. Letters, digits, `-` and `_` only. A season token in the tag must match `season_filter` (a `2025-26` tag on a `2024-25` run stops; a span label such as `2024-26` is accepted for a two-season run), and re-using a tag on the same day appends `-HHMMSS` rather than writing into the earlier run's folder (B46).

Worked values for **2025-26** (data through 2026-09-08; the season ends 2026-09-15), also written out as a paste-ready block in section 7.1 below (they used to sit commented out in `run_config.R`; since 2026-09-12 that file carries only the canonical 2024-25 window, CHANGE_REGISTER D28): window 2025-09-16 to 2026-09-15, `season_filter = "2025-26"`, closure 2025-09-16 to 2025-11-30 with pots legal 2025-12-01, census window 2025-12-01 to 2026-01-03 (Grays Harbor's commercial fishery opened Jan 4, 2026), `commercial_opener = "2026-01-04"`. As of 2026-09-28 the season's OSP record is PARTIAL: 23 OSP days, 17 Sep to 18 Oct 2025 (14 paired with trailer counts), and one WBL boat I/E day (2025-12-19); the rest is awaited from OSP. (Until 2026-09-28 this paragraph said the season had no OSP rows and no WBL boat I/E days, which was wrong.) Corrected 2026-10-02 (this said the shared turnover would fire): **the floor is applied PER FIT** (`prep_bss_crab_pooled()` counts the days that can inform L in that fit), and all 23 OSP days fall in the POT-CLOSURE sub-season (to 30 Nov). So the boat ALL-GEAR fit, about 47% of the port on 2024-25, has NO OSP day: no OSP effort stream, trailer counts only, and the shared turnover REFUSED (per-day draws around the calibration centre, the structure `shared_tau` was adopted to remove; since 2026-10-02 the refusal is written to `run_warnings.csv`). The calibration centre itself comes from 14 autumn pairs (3.50, against 2024-25's 3.03 from 61). The boat pot-closure fit has 2 interviews, below the floor of 15, and reports its PE. A citable 2025-26 boat figure needs OSP's March to September 2026 days (CHANGE_REGISTER D18). The raw workbook and every built input also stop at 2026-09-08 and `ingress_egress.xlsx` at 2026-08-28, so drop the final season export into `raw/`, rerun `build_all_inputs.R` and update the I/E workbook before a production run (section 1).

Also know that the drivers' own `params_model` blocks carry per-fit sampler settings (iterations, treedepth, adapt_delta keyed by fit name) tuned on 2024-25 geometry. They are conservative, so they rarely need touching, but a new season's pathological fit is tuned there, not in `run_config` (`bss_sampler_override` is the sanctioned route from the config side).

## 2. The naive first run

Run the harness first, on the configuration AS SHIPPED, before editing `run_config.R`, from a shell at the repository root; it takes about two minutes, pins the shipped invariants and prints its assertion count and every FAIL. After you change the window, about four checks that pin the canonical 2024-25 window FAIL by design (`canonical: the nine per-season keys are the single 2024-25 season`, the two WINDOW-pin agreement checks, and for a span the blocked-span check), and turning the ladder on fails `shipped: ar_escalate OFF`; any other FAIL is real:

```sh
Rscript 06_diagnostics/test_improvements_2026-08-25.R
```

Then configure the ladder for discovery and run the pooled driver once (`Rscript run_estimation.R --model pooled`: `run_config.R` ships `model = "both"`, and the gear driver honours `ar_escalate` too, so a plain `source("run_estimation.R")` would ladder both tracks and roughly double the cost). These are ELEMENTS of the `run_config <- list(...)` in `run_config.R`, not top-level assignments: edit the existing entries in place (a top-level `ar_escalate <- TRUE` would never reach the run). Only the first differs from the shipped value:

```r
# inside run_config <- list( ... ) in run_config.R
  ar_escalate             = TRUE,          # every fit climbs the ladder (ships FALSE)
  ar_escalate_stop        = "first_pass",  # (the shipped value)
  ar_rung_adequacy        = TRUE,          # per-rung p_loo / Pareto k / coverage, the decisive columns (shipped)
  ar_escalate_respect_cap = FALSE,         # ignore the 2024-25 caps; that is the point (shipped)
```

**Cost, stated plainly:** every rung is a full MCMC fit, and coarsening barely reduces per-fit cost (the 2026-09-04 ladder's rungs cost 96/91/83 minutes each). Budget one full fit per rung per component. Scope the ladder (`ar_escalate = list(shore = "all_gear")`, or a character vector of populations) when you already trust some components.

**What "first_pass" answers:** the ladder starts each fit at daily and coarsens until the convergence gate passes, so the reported rung is the FINEST resolution at which that fit samples, which is exactly the "minimum AR necessary to pass" question a naive season asks. `ar_escalate_stop = "all_rungs"` instead fits every rung and is the full-diagnostic mode, at full cost.

## 3. Reading the run: the decision files, in order

0. **The report's Section 16, "Run warnings and notes" (`run_warnings.csv`), FIRST** (2026-09-29): every fallback (a turnover prior on its 2024-25 fallback, a component on its PE), every skipped diagnostic, every census frame condition (a census clipped to, or outside, the window) and every R warning the run raised, in one table. A naive season is where these fire; nothing below is safe to read until you know which did.
1. **The console/season check** (printed at read time; the report's Section 2 inputs table since 2026-09-29): did the season/window selection capture the data you think it did? Section 3.1 then shows the turnover centres the run used and whether either came from a fallback.
2. **`season_summary.csv` and `convergence_report.csv`:** which fits attempted, which passed, which fell back to PE and why. A PE fallback on a thin component is a measured outcome, not a failure.
3. **`ar_escalation_log.csv`:** one row per rung per fit: resolution, `P_n`, divergences, the gate verdict, each rung's own catch estimate and interval, and (with `ar_rung_adequacy`) per-rung `p_loo`, Pareto k count and coverage. The `selected` flag marks the reported rung.
4. **`model_adequacy.csv`:** the reported fit's adequacy beside the gate. **The gate answers "did this fit sample"; it never answers "is this model right"**, and on a dense season every rung can pass the gate, in which case adequacy is the only thing that separates them.
5. **`pe_vs_bss_comparison.csv`:** the design-based cross-check per component.

**How to choose a rung when several pass** (the 2024-25 ladder, labeled as the worked example, not the answer): all four rungs passed the gate and spanned only 3.7% in catch, so the gate and the estimate decided nothing. What decided it: at daily, `p_loo` was 35% of `n_obs` (one effective parameter per three observations), 41 Pareto k above 0.7 (so its own elpd was unreliable), and effort coverage was +7.1 sampling SD; at weekly all three were clean. Two cautions that generalize: **a narrower prediction interval from a coarser latent process is not precision you earned** (the miscalibrated cells in the 2026-08 boat 2x2 had the narrowest relative intervals), and **an elpd advantage means nothing when the fit's Pareto k have failed**. Prefer the finest rung whose adequacy is clean; use interval width only to break ties between rungs that are ALREADY adequate.

## 4. Pinning the resolutions

Two levers, used in sequence:

- **While deciding:** `ar_force` pins a fit to an exact rung, finer or coarser than the selector would pick (caps can only coarsen, so a pin FINER than the data-driven choice is only expressible here). Scope it per population x sub-season: `ar_force = list(shore = list(all_gear = "weekly"))`. A per-population entry reaches BOTH of that population's sub-seasons; that scoping mistake once moved a component 3,025 crab.
- **Once settled:** move the answer into `ar_max_resolution` and set the matching shore entries of `gear_period_bss` (section 5 of `run_config.R`; the gear track's operative periods, which mirror the pooled shore caps since A32, or the cross-check stops being like-for-like), and return `ar_force` to NULL, so production reads from the cap and the config self-documents. The two routes were proven to build byte-identical Stan data (2026-09-07), and the adoption pattern to copy is `run_adoption_2026-09-07.R`: render once with the new caps and require bit-identity with the pinned run.

## 5. The production run and the cross-check

With the ladder off and the caps set, `source("run_estimation.R")` is the production run: with `model = "both"` (shipped) it renders the pooled headline, then the gear cross-check (about 34 minutes of sections at its shipped periods, the batch's stage D6), and writes `cross_check_<timestamp>.csv` against the 2% criterion (`--model gear_resolved` re-runs the cross-check alone), and read the two tracks together, with the lesson the 2026-09-07 cross-check taught: **compare the tracks at the same resolution before calling a gap a disagreement.** The two tracks agreed on the shore component to 0.08% at a common resolution while sitting 3.4% apart as configured, because the resolution difference dominates. Since 2026-10-02 (CHANGE_REGISTER A32) the gear track ships at the pooled track's shore periods (`gear_period_bss` per population) and with the zero-inflated shore catch (`catch_zi_tracks = c("pooled", "gear_resolved")`), so the two tracks differ only in how CPUE is modelled, plus the boat pot-closure period (biweekly on the gear track, monthly on the pooled; about 1,300 crab in 2024-25). **On a new window, check the gear shore all-gear fit's chains** (`sampler_diagnostics_shore_all_gear_*.csv`): on 2024-25 the weekly fit without the zero inflation trapped a chain at the shipped seed, and with it passed at both seeds tried. A trapped chain there is D3 reopening, not bad luck; report the cross-check as not like-for-like and say so, rather than changing `bss_seed` until it passes.

**Report the boat all-gear resolution span beside the total (D29).** On 2024-25 the boat all-gear fit was adequate and gate-passing at monthly, biweekly and weekly, the held-out weeks could not choose between them, and the component rose from 47,192 to 53,546 (+6.6% at the port) as the period got finer. Monthly ships; the span is resolution uncertainty outside the posterior interval. On a new season the span needs the ladder in its full-diagnostic mode for the boat all-gear fit (section 2: `ar_escalate = list(private_boat = "all_gear")`, `ar_escalate_stop = "all_rungs"`; one full fit per rung). `"first_pass"` stops at the finest rung that passes the gate and gives no span. On 2024-25 the pooled ladder holding every rung in one process ran out of memory (CHANGE_REGISTER B60); forced renders, one rung per run through `ar_force`, are the reliable route.

## 6. Part-season windows

Supported. What changes:

- The sub-season list adapts (section 1.3): one sub-season when the window misses or sits inside the closure; `all_gear_pre`/`all_gear_post` names when a closure falls mid-window.
- Components with no data in the window fall to their floors (`bss_min_interviews`, post-filter floor) and report PE or empty, with the reason in `convergence_report.csv`.
- **The census is clipped to your window** (D34, fixed 2026-09-29, B52). `estimate_comm_charter()` runs the tally, the commercial/charter interviews and the charter roster over the intersection of the census window (or each `census_windows` entry) with `est_date_start` to `est_date_end`. A window that cuts the census covers only its part (on 2024-25 a window from 2025-01-01 takes 3,756 of the 8,538 crab), and one that misses it carries a zero component that the report labels as an absence, not an estimate; either case leads `census_frame_warnings.csv` and the report's warnings section. Before this fix the whole census was added whatever the window. You still set the census dates to the season's census; you no longer have to trim them to the window.
- Expect coarser feasible AR: fewer days means fewer periods, and the ladder's degenerate-rung dropping will shorten itself automatically.

## 7. Multi-season spans

Supported as of 2026-09-10 (CHANGE_REGISTER A14), **for seasons from 2024-25 onward**: the pre-2024-25 seasons are out of scope (top of this guide; D8), so a span must not include 2022-23 or 2023-24. A span takes four settings that must agree:

1. `est_date_start` / `est_date_end` spanning the whole range, and `season_filter` as a vector, e.g. `c("2024-25", "2025-26")`.
2. `pot_closures`: a list with one `list(season =, start =, end =)` entry per closure in the span. This outranks the scalar `pot_closure_start/end` pair and yields one pot-closure sub-season per season (`ring_net_only_<season>`, biweekly, pots excluded) with the following all-gear block as `all_gear_<season>`; a window starting before the first closure gets `all_gear_pre_<season1>`. Overlapping closures, an unlabeled entry, or end-before-start stop loudly. With zero or one in-window closure the legacy scalar path runs byte-for-byte, so single-season names, fit labels, and filenames never change.
3. `census_windows`: a NAMED list, season to `c(start, end)`, for the commercial/charter census; the census is expanded per window and summed, and each season's own component is kept for the season table. An unnamed list stops.
4. The data. Every season in the span needs rows in EVERY workbook of section 0: effort counts, interviews, ingress/egress, holidays, opener dates, and the tally. The validator prints per-season capture and warns hard on the asymmetric case (interviews present, zero effort counts: that season's effort would be pure imputation, not estimation); the holiday reader stops on a missing season.

The report then writes `season_totals.csv` (always) and renders a season-summary table (when the span has more than one season). Monthly figures and fit diagnostics already cover the full span: month labels are year-qualified (`%Y-%m`) and the calendar indices are span-safe (sequential year-week and year-month factors; no aliasing).

Two things to know about a span. The first is NOT what this section used to say: `pot_open_date` was described as a single scalar feeding the ingress/egress `L_effective` regression split, exact only for the first season. It never did that; `estimate_L_effective()` took it as an argument and never referenced it, and the dead argument is gone (corrected 2026-09-12, CHANGE_REGISTER A14). **The real approximation is that the derived shore turnover POOLS every season's I/E days into one diel profile (D24)**, so check `shore_turnover_by_day.csv` by season before trusting a span whose seasons behave differently. (Under a time-denominated shore unit the I/E day-length regression pools the same way, into one `yday -> L` curve; its residuals by season are in `L_effective_ie_detail.csv`, and `ie_filter_by_season = TRUE` is the per-season alternative. Since B55 neither is computed under the shipped gear-deployment unit.) And a span is NOT one statistical process: `build_subseasons()` emits one sub-season per closure and the driver fits each one independently (its own effort, CPUE and turnover parameters), so nothing is shared across seasons except the pooled shore turnover derivation (and, under a time unit only, the I/E day-length regression) and the config priors. A season's total from a span therefore equals a standalone single-season run of the same data up to those two inputs and Monte Carlo error; the value of a span is one run, one report and one set of figures, not borrowed strength between seasons. (Corrected 2026-09-08; the first version of this section claimed a shared process.)

### 7.1 Paste-ready window blocks

`run_config.R` ships the canonical single 2024-25 season and nothing else, so the alternative
windows live here. Each block replaces the nine per-season keys in section 1.2 of that file as a
set; change one, change them all.

**The canonical run (shipped; no edit needed).** The window of the authoritative run,
`05_output/20261003/pooled-CPUE-2024-25`, port total 87,932 [75,193, 105,271], with its gear
cross-check `05_output/20261004/gear-type-CPUE-model-2024-25` (87,903, -0.03%), which are the renders of
`run_config.R` exactly as shipped, `model <- "both"` (B65; and of the runs before it, 99,822 on 2026-09-29
before A33 and A34, 99,873 on 2026-09-28, 96,118 on 2026-09-27 and R4, 94,376). The window keys below are unchanged:

```r
  est_date_start    = "2024-09-16",   est_date_end      = "2025-09-15",
  season_filter     = "2024-25",
  pot_closures      = NULL,
  pot_closure_start = "2024-09-16",   pot_closure_end   = "2024-11-30",
  pot_open_date     = "2024-12-01",
  census_windows    = NULL,
  census_start_date = "2024-12-01",   census_end_date   = "2025-02-08",
  commercial_opener = "2025-02-11",
  run_tag           = "2024-25",
```

**Single season 2025-26** (data through 2026-09-08; the season ends 2026-09-15). Grays Harbor's
commercial fishery opened Jan 4, 2026 (WDFW news release 2025-12-29), so the census window ends
Jan 3:

```r
  est_date_start    = "2025-09-16",   est_date_end      = "2026-09-15",
  season_filter     = "2025-26",
  pot_closures      = NULL,
  pot_closure_start = "2025-09-16",   pot_closure_end   = "2025-11-30",
  pot_open_date     = "2025-12-01",
  census_windows    = NULL,
  census_start_date = "2025-12-01",   census_end_date   = "2026-01-03",
  commercial_opener = "2026-01-04",
  run_tag           = "season-2025-26",
```

Known before you run it (corrected 2026-09-28; this paragraph said no OSP rows existed):
the 2025-26 OSP record is PARTIAL. `WBL_boat_counts.xlsx` holds 23 OSP days, 17 Sep to 18 Oct
2025, 14 of them paired with Westport trailer counts, and one WBL boat I/E day (2025-12-19); the
rest is awaited from OSP (CHANGE_REGISTER D18). Corrected 2026-10-02 (this said the 23 days
would make the shared turnover fire): the floor is applied per fit and all 23 OSP days are in the
pot closure, so the boat ALL-GEAR fit has no OSP stream and its shared turnover is REFUSED
(per-day draws around a centre of 3.50 from 14 autumn pairs; the refusal is in
`run_warnings.csv`), and the boat pot-closure fit (2 interviews) reports its PE. Do not cite a
2025-26 boat or port figure before the full OSP record is in. A smoke fit of the boat all-gear
component on 2026-09-10 (before A30, A31 and B44 to B62) sampled in 4 minutes with 37 of 2,000
divergences. The tally has 23 days and the charter roster 11 Westport trips.

**A two-season span, 2024-26: the mechanical example** (not yet a run to make: wait until the
2025-26 data, its OSP record included, are complete). The multi-season keys, per season in
`pot_closures` and `census_windows`, with the scalar keys on the last season:

```r
  est_date_start    = "2024-09-16",   est_date_end      = "2026-09-15",
  season_filter     = c("2024-25", "2025-26"),
  pot_closures      = list(
    list(season = "2024-25", start = "2024-09-16", end = "2024-11-30"),
    list(season = "2025-26", start = "2025-09-16", end = "2025-11-30")
  ),
  pot_closure_start = "2025-09-16",   pot_closure_end   = "2025-11-30",
  pot_open_date     = "2025-12-01",
  census_windows    = list(
    "2024-25" = c("2024-12-01", "2025-02-08"),
    "2025-26" = c("2025-12-01", "2026-01-03")
  ),
  census_start_date = "2025-12-01",   census_end_date   = "2026-01-03",
  commercial_opener = "2026-01-04",
  run_tag           = "two-season-2024-26",
```

**The 2023-25 span is OUT OF SCOPE (Matt, 2026-09-28; CHANGE_REGISTER D8, closed).** This guide
staged it from 2026-09-10 (and `run_config.R` shipped it until 2026-09-12, D28; a commented copy
is still in that file). It is not to be run: 2023-24 is a pre-protocol exploratory season (one
peak count per day, insufficient coverage), and it was also blocked on data, since no vessel tally
or charter roster was kept for 2023-24, so its census component would have been 0. The 2022-23
and 2023-24 rows remain in the workbooks as a record, and the four-season desk screen for D32
reads their trailer counts descriptively; neither is a model fit.

## 8. Failure modes, and what each one means

| Symptom | Meaning | Fix |
|---|---|---|
| Run stops at "estimation window contains NO effort counts and NO interviews" | `season_filter` does not match the season rows in the window | Update `season_filter` with the window |
| Run stops at "No crabbing holidays for season(s)" | The holiday calendar has no rows for a requested season | Extend the season/year range in `04_input_files/build_crabbing_holidays.R` and rebuild (`Rscript 04_input_files/build_all_inputs.R`). Never hand-edit `crabbing_holidays.xlsx`: it is built, and a hand edit is discarded by the next rebuild |
| Run stops because the window spans more than one season but `pot_closures` is NULL | `season_filter` names several seasons, or the window is longer than 366 days, with the single-season closure keys | Set `pot_closures` and `census_windows` (section 7), or narrow the window |
| Run stops in `bss_output_dir()` on the run tag | `run_tag` carries a season token (`2024-25`) that is not in `season_filter` | Change `run_tag` with the window |
| `run_warnings.csv`: "Shared boat turnover REFUSED for this fit" | That boat fit has fewer than `shared_tau_min_obs` days that can inform the turnover (OSP + boat I/E); it falls back to per-day draws | Expected on a fit with no OSP days (2025-26's all-gear fit until OSP's spring/summer days arrive, D18); get the OSP record before citing the boat |
| "OSP boat-count days in window: 0" | No OSP coverage; stream absent, boat runs trailer-only, `f` lower bound inert | Expected when OSP did not operate; nothing to fix |
| "tau_shore prior: the turnover derivation is unavailable" | Fewer than `tau_shore_derive_min_days` I/E days with the `time` column | The shore turnover is centred on the 2.477 fallback (log-SD 0.3); the report's warnings section says so. Add I/E days, or set `tau_shore_prior_mu` to a number you can defend |
| "tau_boat prior: no OSP/trailer overlap calibration" | No OSP file, a failed overlap diagnostic, or fewer than 3 paired days | The boat turnover is centred on the 3.03 fallback (the 2024-25 calibration); check `osp_trailer_overlap_calibration.csv` and the warnings section |
| A fit reports "PE fallback" with post-filter interviews below the floor | Too thin to fit | A measured outcome; the PE carries the component |
| Every ladder rung fails the gate for one fit | The component cannot support ANY AR on this data | PE fallback is the answer this season; collect more data |
| Gate passes but `flag_miscalibrated` is TRUE / `p_loo` is a large fraction of `n_obs` | Sampled fine, model too flexible at this resolution | Coarsen via the ladder evidence (section 3) |
| Run stops at "NWS marine hazard archive does not cover zone ... over ..." | the committed archive ends before the window does (every production run reads it since 2026-09-27) | Re-pull with `build_nws_marine_hazards.R` (section 0) and commit the workbook. Setting `marine_hazard_mode = "off"` runs the pre-2026-09-27 model, which is not the method of record |
| "bar_restriction: the NWS logistic imputation was not used (...)" in the run log | Fewer than 30 observed sampler-tick days with archive flags, no variation, or a separated fit | Not a stop: unsampled days take the observed rate (`bar_restriction_source = imputed_mean` in `marine_hazard_flags.csv`, and no `marine_hazard_bar_imputation.csv` is written); read the note before citing the term |

## 9. What this guide does not cover

Site generality (the readers filter to Westport / Grays Harbor locations via `gh_effort_areas` / `gh_creel_location`; a new PORT is a data-mapping exercise, not a season), the frozen Method v1.0 documents (historical; see their banners), and the open modelling items, which live with their evidence in `CHANGE_REGISTER.md`.
