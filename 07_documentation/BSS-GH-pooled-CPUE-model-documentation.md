# Grays Harbor Recreational Dungeness Crab Harvest Estimation

## Method Version 2.0: Pooled CPUE Model

**Author:** Matthew George, Ph.D.
**Contact:** matthew.george@dfw.wa.gov
**Agency:** Washington Department of Fish and Wildlife (WDFW)
**Status:** Operational, **not published**. This is the internal method of record for estimating recreational Dungeness crab harvest at Westport / Grays Harbor. WDFW has released no estimate from this pipeline; "method of record" means the method the working model implements, and there is no external figure that a change here has to stay consistent with.
**Method version:** 2.0, adopted 2026-09-12. Method v1.0 (frozen against pooled code v7.4) is archived at `archive/method-v1.0-pooled-CPUE.md`, with a table of the nine places the two methods differ.
**Reference season:** 2024-25, the development test season. The pipeline runs on any window: a full season, part of one, or a multi-season span.
**Reference run:** `05_output/20260910/pooled-CPUE-IMP-R4-shore-tau-newf`.
**Convention:** no em dashes.

> ### THE NUMBERS IN THIS DOCUMENT, AND WHERE THE CURRENT ONES LIVE
>
> Every figure quoted here comes from the reference run named above, and it was the
> authoritative run when this document was written. **It may not be the authoritative run
> when you read it.** The authoritative run and its total live in exactly one place, the box
> at the top of `development_notes/PIPELINE_STATUS.md`. Check there before quoting a number
> from anywhere in this repository, including this file.
>
> Unlike Method v1.0, this document is **not frozen**. It describes the model that
> `run_config.R` ships. When an adopted change moves the method, this file moves with it and
> the change is logged in `development_notes/CHANGE_REGISTER.md`. The two
> `*-development-history.md` files are the backward-looking version logs.

---

## How to read this document

Three parts, and you probably want one of them.

- **Part I is for everyone**: what the estimate is, what data it rests on, how the pieces fit
  together, and what changed from Method v1.0. No equations.
- **Part II is for whoever runs it**: prerequisites, the step-by-step, how to tell whether a
  season's estimate is trustworthy, and what every output file is.
- **Part III is the technical reference**: the full generative specification, the
  design-based estimators, the gate, the priors, and the limitations. Equations, parameter
  names and file line references.

Terms are defined in the glossary (Section 22). Two that trip people up immediately:
**effort** here is measured in **gear-deployments**, not hours; and the **turnover** `tau` is
a dimensionless count of trips per gear slot per day, not a day length.

---

# PART I: FOR EVERYONE

## 1. What this method produces

One number with an interval, for a chosen window, broken into its parts.

For the 2024-25 season the reference run gives a **port total of 94,376 crab, 95% credible
interval [77,566, 118,602]**, assembled from five components:

| component | Bayesian estimate | design-based (PE) | PE relative to BSS | share of port |
|---|---:|---:|---:|---:|
| shore, pot closure | 8,963 | 8,591 | -4.2% | 9.5% |
| shore, all gear | 29,210 | 29,737 | +1.8% | 30.9% |
| private boat, pot closure | 1,372 | 1,192 | -13.1% | 1.5% |
| private boat, all gear | 45,604 | 37,018 | -18.8% | 48.3% |
| commercial + charter vessels | 8,538 | (the same; not modelled) | n/a | 9.0% |
| **port total** | **94,376 [77,566, 118,602]** | **85,076** | **-9.9%** | |

Two things about that table to get right straight away.

**The component medians sum to 93,687, not to 94,376.** That is correct and not a rounding
error: the port total is the median of the summed posterior draws, which is not the sum of
the component medians. Every interval in this document is a posterior quantile, so intervals
do not add either.

**Every percentage in this document is stated PE relative to BSS**, so a negative number
means the design-based estimate sits below the Bayesian one. The output file
`pe_vs_bss_comparison.csv` reports the opposite direction in its `effort_diff_pct` and
`catch_diff_pct` columns, and the column names do not say so. Check the direction before
quoting either.

Alongside the total, a run produces: monthly estimates, catch split by gear type, a
side-by-side design-based comparison for every component, per-fit convergence and adequacy
diagnostics, and the daily latent effort and catch-rate series with intervals.

## 2. The fishery and study area

Recreational Dungeness crab fishing at Westport and in the greater Grays Harbor area, taken
by three distinct groups of people who have to be estimated separately because they are
counted differently:

- **Shore crabbers** fish from the docks (Float 20 and Floats 17-21), the jetty and the
  beach, mostly with ring nets, snares and folding traps, and with pots once pots are legal.
- **Private boat crabbers** launch trailered boats at the Westport Boat Launch and Ocean
  Shores Boat Launch. Most private boat effort is pots.
- **Commercial and charter vessels** fish recreationally out of Westport before the coastal
  commercial season opens, and are counted vessel by vessel rather than modelled.

The season runs mid-September to mid-September. Pots are illegal from the season start until
Nov 30 and legal from Dec 1, which is a structural break in both effort and catch rate, and
the estimate is built on either side of it separately (Section 5).

## 3. The data streams

Seven streams reach the estimate. The first four are counts and interviews the creel program
collects; the fifth is new in Method v2.0; the last two are frames for the vessel component.

| stream | what it is | what it identifies |
|---|---|---|
| **Shore gear counts** | instantaneous counts of gear in the water at the dock, jetty and beach, up to three per day | shore effort |
| **Boat trailer counts** | instantaneous counts of trailers in the censused launch lot | boat effort |
| **OSP daily port count** | the total number of private vessels returning to the Westport Boat Launch each day, supplied by Oregon State Police | boat effort, and the within-day boat turnover |
| **Interviews** | creel samplers' completed-trip interviews: gear count, hours, catch by species and fate, trip type | catch per unit effort, the gear-per-group ratios, and the crabbing fraction |
| **Ingress / egress (I/E) surveys** | timed arrival and departure counts at the dock, with the hour recorded | the shore turnover, and the effective day length |
| **Commercial / charter vessel tally** | the daily count of commercial and charter vessels landing, on the days a sampler was present | the commercial census |
| **Charter trip roster** | every charter trip, marked interviewed / missed / cancelled | the frame the charter component is expanded over |

**The OSP stream is the reason Method v2.0 exists, and it is half-delivered.** OSP will
provide, per day, (a) the total number of vessels returning and (b) the fraction of those
that were **crabbing only**. **(a) is in hand** and is what the boat effort model now runs
on. **(b) is still outstanding.** Note carefully what (b) is when it arrives: it
deliberately excludes combo trips that crabbed alongside another fishery, so it is a **lower
bound** on the vessels that did any crabbing, not the crabbing fraction itself. The machinery
that would use it as a bound is built, tested and inert (Section 14.3); do not wire the
crab-only column in as if it were `f`.

**The interview stream is what makes the crabbing fraction possible.** Every private boat a
sampler contacted at a launch site is recorded with whether it was crabbing and what kind of
trip it was. On 2024-25 that is **300 boats on 107 days, of which 143 were crabbing (a raw
share of 0.477)**, and of the 300, **100 crab-only, 43 combo, 107 another fishery, 50 not
fishing**. Those contacts are what the monthly crabbing fraction is fitted to. Method v1.0
had no such data and used a flat 0.30.

## 4. How the estimate is built

Two estimators run on the same data, and the estimate is the fusion of them, decided per fit.

**The Point Estimator (PE)** is a classical stratified expansion. Sampled days are averaged
within a (stat-week x day-type) stratum and expanded to that stratum's calendar days; catch
comes from the stratum's ratio-of-sums catch rate. It is fast, transparent, makes no
assumption about how effort evolves in time, and it is what a component reports when its
Bayesian fit is not trustworthy.

**The Bayesian state-space model (BSS)** treats daily effort and daily catch rate as latent
quantities evolving smoothly through an AR(1) process, observed through the count and
interview streams. It fills the unsampled days from the temporal structure rather than from a
stratum mean, propagates every uncertainty into the interval, and can use streams the PE
cannot (the OSP counts, the I/E surveys, the sampler contacts).

**The convergence gate decides between them, per fit, on pre-set criteria** (Section 17). A
fit reports its Bayesian posterior only if the sampler demonstrably worked: R-hat below 1.01,
effective sample size above 400, divergent transitions below 5% of draws, and a
divergence-impact test showing that the divergences do not move the answer. Otherwise that
component reports its PE point. In the reference run **all four fits passed and no component
fell back**.

The two are always reported side by side, which is the single most useful validation in the
pipeline: they share almost no machinery, so where they agree, both are more credible. On the
reference run the shore all-gear components agree to **1.3% on effort and 1.8% on catch**,
which is the strongest cross-estimator agreement the project has produced. The boat does not
agree as well (-18.8%), and Section 20 says what is known about why.

## 5. The three population components and the sub-seasons

**Three populations, estimated independently and summed.** Shore and private boat each get a
PE and a BSS fit. The commercial/charter component is **not modelled**: it is a vessel-by-
vessel census plus an expansion (Section 16), with no fit and no per-population output file.
That asymmetry is intentional; it reflects that those vessels are enumerated, not sampled.

**Two sub-seasons, fit separately, split at the pot-open date.** A pot-closure sub-season
(non-pot gear only) and an all-gear sub-season. The split exists because pots becoming legal
is a structural break that one latent process should not bridge, and totals sum across it.
For 2024-25: pot closure Sep 16 to Nov 30 (76 days) and all gear Dec 1 to Sep 15 (289 days).

One naming quirk to know: the pot-closure sub-season's internal key is `ring_net_only`, kept
for output-filename continuity, and it is displayed as "Pot closure".

**A window may span several seasons.** `pot_closures` takes one closure per season and the
report adds season-level totals. Nothing is shared across seasons except the pooled I/E
day-length regression and the config priors; each sub-season is an independent fit.

## 6. What changed from Method v1.0, and what each change was worth

Method v1.0's total for 2024-25, now superseded, was 72,027 [53,018, 101,364]. Method v2.0's is 94,376
[77,566, 118,602]: **+31.0% on the median, with the interval tightening from 67% of the
median to 43%**, and the PE-vs-BSS port gap closing from 37% to 10%.

That is a large move and it was not accepted on faith. It was produced by an **improvement
ladder**: a sequence of runs in which each change is switched on alone, against the same
data, so its effect is measured rather than argued. The four movers sum to the whole within
**24 crab**:

| change | effect on the port | what it replaced |
|---|---:|---|
| the dynamic monthly crabbing fraction `f` | **+11,963** | a flat `f = 0.30` |
| the shore turnover, derived from the I/E `time` column | **+10,327** | `tau_shore` fixed at 1.7 |
| the boat turnover, recentred on the OSP/trailer calibration | **+2,621** | `tau_boat` fixed at 1.2 |
| the census split: a commercial census plus a charter expansion | **-3,283** | one day-type census expansion |

**The crabbing fraction is the largest mover and the best supported.** The retired flat 0.30
was wrong by roughly 3x in the winter months that carry most of the boat catch. Fitted
monthly, with the sampler contacts observing it:

| month | contacts | raw share | fitted `f` | | month | contacts | raw share | fitted `f` |
|---|---:|---:|---:|---|---|---:|---:|---:|
| 2024-12 | 31 | 0.97 | 0.92 | | 2025-05 | 17 | 0.47 | 0.51 |
| 2025-01 | 26 | 0.89 | 0.89 | | 2025-06 | 13 | 0.46 | 0.47 |
| 2025-02 | **5** | 1.00 | **0.87** | | 2025-07 | 17 | 0.47 | 0.42 |
| 2025-03 | 13 | 0.69 | 0.75 | | 2025-08 | 44 | 0.27 | 0.30 |
| 2025-04 | 26 | 0.62 | 0.62 | | 2025-09 | 29 | 0.14 | 0.20 |

Every informed month (20 or more contacts) tracks its own raw share within 0.06. The five
thin months carry visibly wider intervals (mean width 0.360 against 0.249). February is the
clearest demonstration that the walk is doing the right thing: 5 contacts, all 5 crabbing,
and it reports 0.87 rather than 1.00, shrinking toward January (0.89) and March (0.75),
where the retired Beta(6,14) prior would have dragged it toward 0.30.

**And `f` was proved to enter in exactly one place.** A control rung re-ran the model with
the `f` block rolled back and nothing else changed. Over 1,957 shared boat parameters the
largest standardised difference was 3.20, with 2 rows above 3 (0.10%): `f` enters the boat
generated quantities and touches neither the effort nor the CPUE posterior. The boat total is
therefore exactly linear in `f`, which is what makes the +11,963 attributable.

**Two design guarantees were measured in the same ladder rather than assumed.** The shore
components are bit-identical across the rungs that change only boat quantities (6,270 shared
parameter rows at full precision), and the boat components are bit-identical across the rung
that changes only the shore turnover (4,170 rows). A change that should not reach a component
demonstrably does not.

**The cross-check held.** The gear-resolved model, an independent implementation with a
per-gear CPUE structure, read 93,274 on the same configuration: **-1.17%**, inside the
pre-set 2% criterion. The shared turnover agreed to 0.02% and the monthly `f` to 0.002 across
the two parameterizations.

## 7. Where this method is valid

**It is valid for the fishery it was built on**: recreational Dungeness crab at Westport and
Grays Harbor, with the current sampling design (instantaneous gear and trailer counts,
completed-trip interviews at the launch and dock sites, an I/E survey with the hour recorded,
a vessel tally, a charter roster, and the OSP port count).

**It is not a general creel estimator you can point at another fishery.** Four things are
specific:

1. **The effort unit is gear-deployments.** Catch is sub-linear in soak time for pot and trap
   gear, so time-denominated units fail the pipeline's own linearity test. A fishery whose
   catch is linear in time needs a different unit and a different `L`.
2. **The turnover structure assumes gear slots turn over within a day.** That is what `tau`
   is, and it is estimated from the I/E and OSP series. A fishery without that structure has
   no `tau`.
3. **The crabbing fraction exists because boat counts are all-boat counts.** A fishery where
   the effort count is already directed at the target species does not need `f`, and applying
   one would be wrong.
4. **Every cap, floor and prior centre in `run_config.R` was derived on 2024-25** and is
   tagged `SEASON-DERIVED`. On a new season they are starting points, not answers; the
   workflow for re-deriving them is `NEW_SEASON_GUIDE.md`.

**It must be re-derived, not merely re-run, if:** the sampling design changes (different
count protocol, different sites, a different shift structure); the gear vocabulary changes
(see the 2022-24 ring-net label problem, CHANGE_REGISTER D20); or the effort frame changes.
On that last point, note now that the current frame is known to be incomplete: private boats
moored at the marina floats are interviewed and their catch rates enter the CPUE, but their
effort is in neither the trailer count nor the OSP ramp total (Section 20).

---

# PART II: RUNNING IT

## 8. Prerequisites and repository layout

**Software.** R 4.2 or later, and **rstan** 2.32 or later (this pipeline uses rstan, not
cmdstanr). Plus tidyverse, lubridate, suncalc, gt, patchwork, here, readxl, loo.
`run_estimation.R` installs anything missing, so a fresh machine's first run may trigger a
long Stan compile. The reference run was made with rstan 2.32.7 / StanHeaders 2.32.10.

**Runtime.** About 3 to 6 hours on 4 cores for a full pooled run: four real MCMC fits plus
the diagnostics. The gear-resolved cross-check is another 3 hours or so.

**Repository layout**, as a numbered pipeline:

| folder | contents |
|---|---|
| `01_BSS_models/` | the two production driver `.Rmd` reports |
| `02_stan_models/` | the Stan models; this method uses `crab_bss_pooled.stan` |
| `03_R_functions/` | the shared helper library, sourced whole by every driver |
| `04_input_files/` | the nine input workbooks, and the builders that generate six of them |
| `05_output/` | one dated folder per run |
| `06_diagnostics/` | the regression harness, the dated batch runners, the weather module |
| `07_documentation/` | this file and the rest of the reference layer |
| `run_config.R` | **the one file you edit** |
| `run_estimation.R` | the orchestrator |

**The inputs, and which are built.** Every input is an `.xlsx` workbook with a single `data`
sheet and ISO `yyyy-mm-dd` dates. **Six of the nine are BUILT** from the per-season creel
workbooks in `04_input_files/raw/` by `04_input_files/build_all_inputs.R`, and must not be
hand-edited: `interview_combined.xlsx`, `effort_combined.xlsx`, `sampler_shifts.xlsx`,
`wes_commercial_tally.xlsx`, `charter_trips.xlsx`, `crabbing_holidays.xlsx`. An edit to one of
those is silently discarded the next time anyone adds a season; fix the raw workbook or the
builder. The three that are maintained by hand, because their sources are not the season
workbooks, are `ingress_egress.xlsx` (the I/E database export), `WBL_boat_counts.xlsx` (the
OSP port counts) and `fishery_opener_dates.xlsx` (the regulation calendar, used only by a
diagnostic).

**Data quirks that are real and are matched in code.** Interview `number_of_gear` maps from
column N, not W, because of a duplicate iForm field name. The commercial `boat_type` is the
typo "Commerical", one m, matched by regex; do not correct the spelling without updating the
matcher. Windows and OneDrive long paths can exceed MAX_PATH, which the code detects and
works around.

## 9. Step-by-step: running a season

**One command runs everything**, after you edit one file.

```r
source("run_estimation.R")        # RStudio: use Source, NOT Knit
```
```sh
Rscript run_estimation.R                          # terminal / unattended
Rscript run_estimation.R --model gear_resolved    # the cross-check
```

1. **Add the season's data.** Drop the creel workbook into `04_input_files/raw/` as
   `<YYYY><YY>_rec_crab_harvest_data.xlsx` and run
   `Rscript 04_input_files/build_all_inputs.R`. Read each builder's report before committing.
2. **Edit `run_config.R`, section 1.2**: the nine per-season keys, as a set. Paste-ready
   blocks for the canonical 2024-25 window, for 2025-26, and for the 2023-25 span are in
   `NEW_SEASON_GUIDE.md` section 7.1.
3. **Run the harness first.** `Rscript 06_diagnostics/test_improvements_2026-08-25.R` takes
   seconds, needs no rstan, and pins every shipped invariant. It has caught a committed
   `DRY_RUN <- FALSE` five separate times.
4. **Run naively, with the AR ladder on**, if the season is new. Then read
   `ar_escalation_log.csv` and `model_adequacy.csv` and pin each fit's resolution
   deliberately. The full workflow, including what each diagnostic decides, is
   `NEW_SEASON_GUIDE.md`.
5. **Produce**, then run the gear-resolved cross-check and compare the port totals. The
   pre-set criterion is agreement within 2%.

**As shipped, `run_config.R` is the canonical run**: the single 2024-25 season, pooled model,
Method v2.0 throughout. It reproduces the reference run named at the top of this document.

**Do not edit the `.Rmd` drivers or the `.stan` files for a routine run.** `run_config.R` is
the single control surface, and it is ordered so that section 1 is the run, section 2 is the
method, and sections 3 to 5 are diagnostics and levers you will rarely touch. Three toggles
force a Stan recompile when changed: `razor_dig_mode`, `estimate_cpue_density` and
`estimate_catch_zi`.

## 10. Judging whether a season's estimate is trustworthy

Read these five things, in this order. The first two are pass/fail; the last three are
judgement.

**1. Did every component report BSS?** `convergence_report.csv`, column `method_selected`. A
component that reads `PE (convergence fail)` or `PE (insufficient data)` contributed a point
estimate with no interval to the port total, which widens nothing and hides everything. On
the reference run all four fits reported BSS.

**2. Is anything flagged in `model_adequacy.csv`?** Five flags, and each means something
different. `flag_overparameterised` (`p_loo` above 25% of observations) says the fit is
spending effective parameters at a rate the data cannot support. `flag_miscalibrated`
(50% coverage off by more than 0.15) says the intervals are the wrong width.
`flag_loo_unreliable`, `flag_pit_bias` and `flag_dispersion_neff` are the others. **On the
reference run nothing is flagged**, `p_loo` runs 7 to 18% of observations, and the only
`flag_loo_unreliable` cases are one bad Pareto k each on shore all-gear and boat pot closure.

**3. How much of the estimate is extrapolated, and from what?** `fit_data_summary.csv` gives
the sampled-day fraction per component. On 2024-25 the shore effort series covers 197 of 365
days and the boat 186 of 365, but the **boat CPUE** rests on far less: roughly 200 usable
private-boat interviews across the season, on about a quarter of the days. The boat is half
the port total. No modelling choice manufactures uncollected information.

**4. Do the PE and the BSS agree?** `pe_vs_bss_comparison.csv`. They share almost nothing, so
agreement is real evidence. Two cautions. The annual agreement is not monthly agreement:
on the reference run the shore annual totals agree to 1.8% while the PE over-allocates
January to March by about 2.2x and under-allocates June and July, and the errors cancel. And
the comparison depends on a lever: the boat PE sits 18.8% below the boat BSS under the
shipped `pe_empty_stratum = "local"` and 6.1% below under `"pooled"` (Section 20).

**5. Is the PE resting on thin cells?** `pe_empty_effort_strata.csv`. With weekly strata and
about 50% day coverage, **roughly 44% of the shore all-gear component's calendar days and 43%
of the boat's sit in a stratum with one sampled day or none**. That is priced (Section 15),
but on a thinner season it gets worse, and if a component's gate then fails, those cells
reach the headline.

## 11. Output catalog

Each run writes to `05_output/<YYYYMMDD>/<model>-<run_tag>/`. Outputs are committed to git on
purpose, so past estimates are preserved as produced; only the per-run `*.RData` workspaces
and compiled Stan caches are ignored. A given dated folder may be a partial run; use a recent
complete run as the reference catalog. `05_output/README.md` is the full per-file inventory;
what follows is what to read first.

**The estimate**

| file | contents |
|---|---|
| `port_total_<species>_<fate>.csv` | the port total: PE, BSS median, and the 95% interval, for expected catch and for predictive catch |
| `season_totals.csv` | per-season totals (always written; the season table renders on a multi-season span) |
| `pe_port_summary.csv` | the design-based port total, component by component |
| `pe_vs_bss_comparison.csv` | the side-by-side, per component, with the census row's construction spelled out |
| `monthly_pe_vs_bss.csv` | the same comparison split by month |
| `catch_by_gear_type*.csv` | catch apportioned to gear types from interview shares |

**The fits**

| file | contents |
|---|---|
| `convergence_report.csv` | the gate, per fit: every criterion, its value, its verdict, and `method_selected` |
| `model_adequacy.csv` | the adequacy statistics and the five flags, reported BESIDE the gate |
| `bss_summary_*.csv`, `bss_full_summary_*.csv` | posterior summaries per fit |
| `bss_daily_effort_*.csv`, `bss_daily_catch_*.csv`, `bss_daily_cpue_*.csv` | the latent daily series with intervals |
| `bss_draws_summed_*.csv` | the summed draws each component contributes to the port |
| `ar_escalation_log.csv` | one row per fit attempt: resolution, gate verdict, that rung's own estimate and interval, and per-rung adequacy |

**The method-of-record quantities**

| file | contents |
|---|---|
| `crab_fraction_strata_<fit>.csv` | every `f` stratum: its label, the contacts that informed it, the combo share, and the posterior median and interval |
| `shore_turnover_summary.csv` | the derived shore turnover, its bootstrap SE, the in-window day count and the window-only alternative |
| `osp_trailer_overlap_calibration.csv`, `osp_trailer_overlap_pairs.csv` | the boat turnover calibration and the paired overlap days |
| `census_daily.csv`, `census_variance.csv` | every census day flagged observed or no-operation, the two components separately, and the per-stratum variance |
| `pe_empty_effort_strata.csv` | how much of each PE component rests on thin cells, and what the three levers did |
| `L_effective_ie_detail.csv`, `bss_L_effective_*.csv` | the I/E day-length regression and the fitted turnover |

**Provenance**

`run_parameters.txt` (the exact config, untruncated), `session_info.txt` (R and package
versions), and `run_manifest_<timestamp>.txt` one level up (model, git SHA, per-stage timing,
a `str()` dump of the config, and full `sessionInfo()`).

## 12. Diagnostics: what each one answers

Every diagnostic in this pipeline is `tryCatch`-wrapped, so a failing diagnostic cannot abort
a multi-hour fit, and none of them changes an estimate.

| diagnostic | the question it answers |
|---|---|
| `cpue_linearity_*`, `cpue_saturation_*` | is catch linear in the effort unit? This is what justifies gear-deployments and rejects crabber-hours; it is re-measured every run, and any new gear type must pass it before its totals are trusted |
| `model_adequacy.csv` | is the model carrying the data, as distinct from did the sampler work |
| `ar_escalation_log.csv` | what would a finer or coarser AR resolution have given, and would it have passed |
| `sensitivity_incomplete_trips.csv` | what would the four treatments of interrupted trips do to the estimate and to the gear ratios |
| `effort_overdispersion_decomp_*.csv` | is the effort variance Poisson floor, observation overdispersion, or latent process? This says which lever would help |
| `osp_coverage_audit.csv` | which days the OSP stream covers, and how its counts compare with the trailer counts on the overlap days |
| `shift_coverage_*.csv`, `contact_hour_*.csv` | what fraction of a day's boat returns the sampler shifts cover, and whether the trip-type mix drifts with the hour inside them |
| `cpue_saturation_*`, `ppc_calibration_*`, `ppc_byobs_*` | posterior predictive checks: does the fitted observation model reproduce the observed count distribution |
| `pe_empty_effort_strata.csv` | how much of the design-based estimate is imputed, and what it would be under the other levers |
| `tau_sensitivity_*` | how the boat component moves across a grid of turnover prior centres |
| `fishery_opener_spillover_*` | do other fisheries' openers coincide with effort surges that a constant crabbing fraction would mis-convert |

## 13. Reproducibility

`bss_seed` is fixed and should stay fixed. A fixed-seed re-run of the same configuration on
the same data reproduces the same fits, which is what makes the ladder's bit-identity claims
possible.

Two caveats that matter when comparing runs.

**A configuration change that resizes the parameter vector breaks bit-identity legitimately.**
Most optional features in this model use a zero-size-when-off declaration precisely so that
switching them off leaves the unconstrained parameter vector unchanged and a fixed-seed rerun
reproduces byte for byte (Section 14.9 lists which). Some features do not, and for those a
baseline reproduction has to be judged on medians and intervals within Monte Carlo error.

**The port total resamples.** Component draws are permuted when the port is assembled, so the
port total moves by about 0.2% between bit-identical fits. The gate is judged on per-fit
posterior summaries, never on the port line.

Every run records its git SHA, its full config and its session info (Section 11). The
improvement ladder additionally records a three-layer **code fingerprint** (Stan models,
drivers, R functions), so a cross-run comparison that rests on two folders having been fitted
by the same code can be checked rather than assumed.

---

# PART III: TECHNICAL REFERENCE

Everything below is written against `02_stan_models/crab_bss_pooled.stan` and the helper
library as they stand. Parameter names are the Stan file's own. Where a statement holds only
under the shipped configuration, the lever is named.

## 14. Model specification (`crab_bss_pooled.stan`)

### 14.1 Notation and dimensions

`D` days in the sub-season, `G = 1` gear group, `S = 1` section, `P_n` AR periods. The Stan
code carries `[g]` and `[s]` loops throughout and the R prep sets both to 1, so every such
loop is a single cell; the indices are kept below because the code keeps them.

One fit is one **population x sub-season**. Four fits per pooled run: shore pot closure,
shore all gear, private boat pot closure, private boat all gear.

**`lambda_E_S` does not have a fixed unit, and the Stan file does not declare one.** Its unit
is set entirely by which observation stream R feeds:

| population | `lambda_E_S` is | because the effort likelihood is |
|---|---|---|
| shore | crabbers | `Gear_I ~ NB2(lambda_E * R_G, r_E)` |
| private boat | gear units in the water | `T_I ~ NB2(lambda_E / R_G_boat, r_E)` |

The authority for that pairing is `03_R_functions/bss_effort_spec.R`, which is the single
source of the effort unit AND the matching I/E observation column, read by both the PE and
the BSS prep so the two cannot drift onto different scales.

### 14.2 The effort process

**The latent intensity.**

```
lambda_E_S[s][d,g] = exp( mu_E[g,s] + omega_E[period[d], gs]
                          + B1 * w[d] + B2 * holiday[d]
                          + X_open[d] . B_open )         (the last term only if K_open > 0)
```

**The level** is two-tier: `mu_E[g,s] = mu_mu_E[g] + eps_mu_E[g,s] * sigma_mu_E`, collapsing
to `mu_mu_E[g]` under `collapse_mu_hier = 1` (off in production).

**The AR(1)** is on `omega_E`, indexed by PERIOD, not by day, so every day in a period shares
one deviation exactly:

```
omega_E[1]   = sigma_eps_E / sqrt(1 - phi_E^2) * omega_E_0_raw        (non-centred, stationary)
omega_E[p]   = phi_E * omega_E[p-1] + sigma_eps_E * eps_E[p-1]        p = 2..P_n
phi_E        = 2 * phi_E_scaled - 1,   phi_E_scaled ~ Beta(2,2)
eps_E        ~ std_normal()
```

`period[d]` and `P_n` come from the resolution (Section 18): daily gives `P_n = D` and
`period[d] = d`; weekly, biweekly and monthly give the corresponding calendar index. **The
Stan code is unchanged across resolutions**; only the index and `P_n` change.

**The day-type effects are NESTED, not mutually exclusive, and this is easy to get wrong.**
The R prep builds `w[d] = 1` on weekends **and** holidays, and `holiday[d] = 1` only on
holidays. So the multiplicative day-type effect is `1` on a weekday, `exp(B1)` on a weekend,
and `exp(B1 + B2)` on a holiday: **`B2` is an increment on top of the weekend effect, not a
separate level.** The same nesting applies to `B1_C` and `B2_C` on CPUE. This follows from
`prep_days_crab.R`, not from anything in the Stan file.

**Expected effort**, the quantity that becomes the estimate:

```
E[s][d,g] = lambda_E_S[s][d,g] * E_scale * L[d] * f_crab[f_stratum[d]]
E_scale   = R_G   if effort_scale_gear == 1   (shore, gear-deployments)
          = 1     otherwise                   (boat)
```

so **E = latent intensity x expansion x turnover x crabbing fraction**, and `E_sum` is its
sum over days. `E` is **crab-directed** effort: `f` is inside it.

**The four streams that inform effort.**

| stream | likelihood | note |
|---|---|---|
| shore gear count | `Gear_I[i] ~ NB2(lambda_E[day,1] * R_G, r_E)` | `R_G` is a MULTIPLIER |
| boat trailer count | `T_I[i] ~ NB2(lambda_E[day,G] / R_G_boat, r_E)` | `R_G_boat` is a **DIVISOR**: `lambda_E` is gear in the water, `R_G_boat` is gear per boat group, so the quotient is boat groups, i.e. trailers. `L` does NOT appear here |
| OSP port count | `OSP_I[i] ~ NB2((lambda_E[day,G] / R_G_boat) * L[day], r_OSP)` | under `osp_scale_is_tau = 1`, production. Its own dispersion `r_OSP`, not `r_E`. `OSP_I` is the daily total of ALL private boats, crabbing or not |
| ingress / egress | `IE_obs[i] ~ lognormal(log(lambda_E[day,1] * L[day]), sigma_IE)` | shore only; see the caution below |

`r_E` is shared between the gear and trailer streams. Two interview streams inform the
expansion ratios indirectly: `Gear_A[a] ~ Poisson(A_A_gear[a] * R_G)` and
`Gear_A_boat[a] ~ Poisson(R_G_boat)`. Because `R_G` is also `E_scale` for shore, the
gear-per-crabber interviews enter shore effort twice, once through the count mean and once
through the expansion.

> **The I/E likelihood asserts a MEDIAN, not a mean.** The first argument of a lognormal is
> its location, so `lambda_E * L` is the median of observed arrivals and the mean is
> `lambda_E * L * exp(sigma_IE^2 / 2)`. That is a genuine asymmetry with the three count
> streams, whose stated means are means. On the reference run's shore all-gear fit
> `sigma_IE = 0.577`, so the implied mean exceeds the median by 1.18x. The Stan file does not
> comment on this, and it is worth keeping in view whenever the I/E stream is used to argue
> about the level of shore effort.

**Effort overdispersion is marginalized.** The earlier form was
`Poisson(lambda * eps * R)` with `eps ~ Gamma(r_E, r_E)`; the gamma-Poisson marginal is
exactly `NB2(lambda * R, r_E)`, so the per-observation latents no longer exist. The
dispersion parameters are reparameterized: `r_E = 1/sigma_r_E^2`, and likewise `r_OSP`,
`r_C`, with the `sigma_r_*` being what is sampled. `NB2(mu, r)` has variance `mu + mu^2/r`.

### 14.3 The crabbing fraction `f` and the combo share `c`

This is the block that defines Method v2.0, so it is specified in full.

**The walk.** `f` is a logit random walk over strata, with year-month strata in production:

```
eta_f[k] = f_level_mu + f_level_sd * z_f[k]                      if k is anchored
         = eta_f[prev[k]] + sigma_f * sqrt(gap[k]) * z_f[k]      otherwise
f_crab[k] = inv_logit(eta_f[k])                                  (squeezed off 0 and 1 by 1e-6)
```

- **Level prior:** `f_level_mu = logit(crab_fraction_set)` with `crab_fraction_set = 0.3`,
  `f_level_sd = 1.5`. The 0.3 is now only the anchored stratum's prior CENTRE, not a value
  the model uses; on 2024-25 no month's posterior sits near it.
- **Innovation SD:** `sigma_f ~ half-normal(0, 1.5)`. On the reference run
  `sigma_f = 0.701`, away from zero, so the walk is identified rather than collapsing to a
  constant.
- **Step scaling:** `sqrt(gap[k])`, Brownian, so a two-month gap takes a step of SD
  `sigma_f * sqrt(2)`.
- **Innovation family:** `z_f ~ student_t(f_walk_df, 0, 1)` when `f_walk_df > 0`, else
  standard normal. **Shipped `f_walk_df = 4`.** Note precisely: `z_f` is the same vector for
  the anchored level deviation and for the innovations, so under the shipped setting the
  LEVEL prior is `f_level_mu + f_level_sd * t_4`, not a normal. The file's own prose writes
  the two lines separately and then gives one distribution for `z_f`, which does not
  distinguish the two roles. This is the code's behaviour, stated without any claim about
  what was intended.
- **Walk order:** built by `crab_fraction_walk_structure()`. Monthly strata chain each month
  to the latest earlier month in the same class. **Under `crab_fraction_strata = "day_type"`,
  `"opener"` or `"none"` there is no walk at all**: every stratum is anchored, an independent
  draw from the level prior, and `sigma_f` enters no likelihood.
- `transformed data` **rejects** a walk whose predecessor does not precede it, rather than
  producing a silently wrong chain.

**The observation model for the sampler contacts**, which is what makes the walk informative:

```
cfi_crab[i] ~ beta_binomial(cfi_total[i],
                            f_crab[stratum[i]] * cfi_kappa,
                            (1 - f_crab[stratum[i]]) * cfi_kappa)
```

**One observation per CONTACT DAY**, not per boat and not per stratum sum. That choice is
load-bearing: a binomial on the stratum sum would treat roughly 7,500 boat-days as
independent trials and return an `f` posterior SD near 0.005, which would be a false
precision, not a result. `cfi_kappa ~ lognormal(log(20), 0.75)`, with the log-SD hard-coded
in the Stan file. The stream observes `f` itself, because a combo trip counts as crabbing
when the sampler saw the crab gear. **The denominators are observed boat counts, never the
latent effort**, which is what keeps `f` out of the effort and CPUE likelihoods.

**The combo-trip share `c`** is a second walk with the same order, the same gaps and the same
innovation family, but its own step SD `sigma_c ~ half-normal(0, 1.5)`, its own level prior
centred at 0.3, and its own concentration `cfc_kappa ~ lognormal(log(20), 0.75)`:

```
cfc_combo[i] ~ beta_binomial(cfc_crab[i],
                             combo_c[stratum[i]] * cfc_kappa,
                             (1 - combo_c[stratum[i]]) * cfc_kappa)
```

where `cfc_crab[i]` is the day's typed crabbing boats and `cfc_combo[i]` how many of those
were combo trips. The `c` walk is live whenever typed contacts exist or the OSP crab-only
stream is on; typed counts come from the interview workbook's `trip_type_class` column.

**`f` enters generated quantities and nothing else.** It appears in the model block only on
observed boat counts, and in the generated quantities at `E`, `lambda_Ctot_S` and `f_crab_out`.
It appears nowhere in `lambda_E_S` or `lambda_C_S`. Therefore **the boat total is exactly
linear in `f`, and the model's CPUE is invariant to it** - which is the property the R2-vs-R2f
control rung measured rather than assumed (Section 6).

**The OSP crabbing-only lower bound: built, tested, and inert.** Its parameters are
zero-size unless `use_osp_crab_lower = TRUE` **and** the legacy (non-dynamic) construction is
in use, and it additionally requires the workbook to carry the crab-only column and at least
one stratum to clear `crab_fraction_osp_min_obs`. Today none of those holds. When it is
active, the legacy form is a genuine hard bound: `f = f_lower + (1 - f_lower) * theta`, so
`f` can never fall below the crab-only share OSP observed directly.

> **`f_lower_out` MEANS TWO DIFFERENT THINGS under the same name, and this is a cross-run
> hazard.** Under the legacy construction it is the hard lower bound just described. Under
> the dynamic construction, which is production, it is **not a bound at all**: it is the
> derived quantity `f * (1 - c)`, the model's prediction of what OSP's crabbing-only column
> should read, and it is pinned to 0 when the `c` walk is off. Any comparison of
> `f_lower_out` across runs must first establish which construction was live. This is also
> the column that will be the direct check on the shift-time contacts when OSP delivers the
> crab-only data.

### 14.4 The turnovers

`L[d]` is the daily expansion factor and it is a **TURNOVER**, not a day length: trips per
gear slot per day, dimensionless. `L_effective` in hours is computed every run but is a
diagnostic; it becomes `L` only if the shore effort unit is set back to a time unit.

```
L[d] = tau_bar * exp(shared_tau_sigma * L_raw[d])     if shared_tau == 1
     = L_data[d] * exp(L_prior_sigma[d] * L_raw[d])   otherwise
L_raw ~ std_normal()
```

**What `shared_tau` fixes.** With `shared_tau = 0` there are `D` independent per-day draws,
each anchored on its own prior centre, with nothing pooling information across days. That
sounds harmless and is not: measured, shore with 4 in-window I/E days out of 289 produced a
season median `L` of 1.6998 against a prior centre of 1.7000, and the boat with 148 OSP days
gave 1.201 against a centre of 1.200. **The prior, not the data, was setting the turnover.**
With `shared_tau = 1` there is one estimated level `tau_bar` with per-day lognormal
deviations of FIXED spread. The spread is fixed deliberately: `shared_tau_sigma` is data, not
a parameter, so the day-to-day scatter cannot trade off against the level on a series where
most days are unobserved.

R refuses `shared_tau = 1` when `L_data` varies across days (a time-denominated unit), or
when fewer than `shared_tau_min_obs = 15` days can inform `L`. **That floor is what makes the
feature boat-only in effect**: on 2024-25 the shore fits have 4 informed days and the boat
has 148.

**Where the two centres come from, and this is the second-largest change in Method v2.0.**

| | shipped setting | 2024-25 prior centre | fitted | what it replaced |
|---|---|---|---|---|
| shore | `tau_shore_prior_mu = "derived"` | **2.477**, log-SD 0.10 (floored) | posterior 2.394 | a literal 1.7 |
| boat | `tau_boat_prior_mu = "calibration"` | **3.030**, SD 0.50, from 61 paired overlap days | `tau_bar` **2.977** | a literal 1.2 |

Both resolve per run and both are printed with their source, so a run says where its turnover
came from rather than leaving it to be looked up: `derived (I/E time column, 40 days,
count-time-weighted ratio of sums)` and `calibration (trailer_mean_per_visit, n = 61)`.

The shore centre is derived from the I/E `time` column: a diel presence profile, evaluated at
the hours the creel counts were actually taken, giving a count-time-weighted ratio of
arrivals to presence. The retired 1.7 was arrivals over PEAK presence, which is a different
quantity, not a competing estimate of the same one. The boat centre is the implied turnover
of the OSP daily total divided by the trailer snapshot on the paired overlap days.

**`L` enters the likelihood in exactly two places**: the I/E lognormal location, and the OSP
mean under `osp_scale_is_tau = 1`. It does not enter the gear or trailer means. That is
precisely the definition of an "informed day" behind the `shared_tau_min_obs` floor.

**One tension the code carries and does not settle.** `kappa_OSP`, the free OSP scale under
`osp_scale_is_tau = 0`, sat near 2.7 while the `tau_boat` prior centred near 1.2, and the
Stan file's own comments record that the roughly 2.5x conflict was being absorbed by the OSP
overdispersion (`r_OSP` near 1.6) rather than moving `L`, showing up independently as a boat
trailer PIT mean of 0.42 against a nominal 0.50. Production runs with
`osp_scale_is_tau = 1`, under which `L` IS the OSP turnover and `kappa_OSP` drops out of the
likelihood entirely while still being sampled from its prior and reported. The file does not
claim the tension is resolved, and neither does this document.

### 14.5 The CPUE process

```
lambda_C_S[s][d,g] = exp( mu_C[g,s] + omega_C[period[d], gs]
                          + B1_C * w[d] + B2_C * holiday[d]
                          + gamma_C * (log lambda_E - log_E_ref) )   (last term off in production)
```

Identical in form to the effort process: a two-tier level, a non-centred stationary AR(1)
initial state, `phi_C = 2 * phi_C_scaled - 1` with `phi_C_scaled ~ Beta(2,2)`, standard-normal
innovations, and the same `period[d]` index and `P_n`. The day-type effects nest the same way
(Section 14.2). The density term `gamma_C` is the rejected same-day-effort interaction; it
ships off.

**Expected catch is formed twice, consistently.** Per interview, in the likelihood:
`mu_c = lambda_C * h[a]`, with `h[a]` that interview's own effort denominator. Per day, for
reporting:

```
lambda_Ctot_S[s][d,g] = E[s][d,g] * lambda_C_S[s][d,g] * zi_scale
```

**The `zi_scale` factor is essential.** Under the zero-inflated mixture `lambda_C` is fitted
to the non-inflated component and rises to absorb the structural zeros, so reporting
`lambda_C * E` unscaled would inflate the season total by `1/(1 - theta_C)`. `zi_scale` is
`1 - theta_C` when the mixture is on and 1 when it is off.

### 14.6 Observation models

| stream | distribution | mean (or location) | dispersion | zero-inflation |
|---|---|---|---|---|
| shore gear count | NB2 | `lambda_E[day,1] * R_G` | `r_E` | none |
| boat trailer count | NB2 | `lambda_E[day,G] / R_G_boat` | `r_E` | none |
| OSP port count | NB2 | `(lambda_E[day,G] / R_G_boat) * L[day]` | `r_OSP` | none |
| I/E effort | lognormal | location `log(lambda_E[day,1] * L[day])` | `sigma_IE` | n/a (continuous) |
| interview catch | NB2, or a two-component ZINB mixture | `lambda_C[day,gear] * h[a]` | `r_C` | `theta_C`, shore only |
| gear per crabber | Poisson | `A_A_gear[a] * R_G` | - | none |
| gear per boat group | Poisson | `R_G_boat` | - | none |
| sampler contacts | beta-binomial | `f_crab[k]`, concentration `cfi_kappa` | - | none |
| typed combo contacts | beta-binomial | `combo_c[k]`, concentration `cfc_kappa` | - | none |

**The zero-inflated catch block, exactly.**

```
if (c[a] == 0)  target += log_mix(theta_C, 0, NB2_lpmf(0 | mu_c, r_C));
else            target += log1m(theta_C) + NB2_lpmf(c[a] | mu_c, r_C);
```

It applies **only to the interview catch stream**. No effort stream is zero-inflated. The
switch is **per fit, not per run**: with `estimate_catch_zi = TRUE` and
`catch_zi_populations = c("shore")`, the mixture is on for the shore fits and off for the
boat fits **within the same run**, which makes the boat fits a deliberate untouched negative
control. `theta_C` is declared zero-size when off, so an off run is bit-identical rather than
merely similar. Prior `theta_C ~ Beta(1, 9)`, which is SEASON-DERIVED from the 2024-25 zero
bin.

The evidence that put the mixture in the method, compared like for like at the same
resolution: shore all-gear zero bin z from +3.7 to +2.0, one bin z from -6.1 to -3.3; shore
pot closure zero bin +2.9 to +0.7, one bin -3.5 to -1.2; elpd +11.6 nats at 2.30 paired SE;
Pareto k above 0.7 from 1 to 0, the only fit in this project with none. What it does **not**
fix is in Section 20.

The pointwise log-likelihood generated for PSIS-LOO mirrors the mixture exactly, because a
mismatched `log_lik` would invalidate every `elpd_loo` comparison. **PSIS-LOO on this model
covers the gear, trailer, OSP and catch streams only**; there is no `log_lik` for the I/E
stream, the two gear-ratio interview streams, or any crabbing-fraction stream.

### 14.7 Priors

| parameter | prior | hyperparameter source |
|---|---|---|
| `mu_mu_E[G]` | `normal(mu, 2)` | `mu = log(25)` shore, `log(10)` boat |
| `sigma_mu_E`, `sigma_mu_C` | half-Cauchy(0, 1) | hard-coded scale 1 |
| `eps_mu_E`, `eps_mu_C` | `std_normal()` | - |
| `mu_mu_C[G]` | `normal(log(0.5), 2)` | - |
| `phi_E_scaled`, `phi_C_scaled` | `Beta(2, 2)` on (0,1), mapped to `phi = 2p - 1` | symmetric about 0 with a mild interior mode |
| `sigma_eps_E`, `sigma_eps_C` | half-Cauchy(0, 1) | - |
| `sigma_r_E`, `sigma_r_C`, `sigma_r_OSP` | half-Cauchy(0, 1) | `sigma_r_OSP` reuses the effort hyperparameter |
| `eps_E`, `eps_C`, `omega_*_0_raw`, `L_raw` | `std_normal()` | - |
| `B1`, `B2`, `B1_C`, `B2_C`, `B_open` | `normal(0, 1)` | - |
| `gamma_C` | `normal(0, 1)` | hard-coded; off in production |
| `R_G` | `lognormal(log(mu), 0.3)` | `mu` = the run's empirical gear-per-crabber ratio |
| `R_G_boat` | `lognormal(log(4), 0.5)` | hard-coded |
| `kappa_OSP` | `lognormal(log(3.0), 0.3)` | SEASON-DERIVED from the 2024-25 overlap days |
| `sigma_IE` | `exponential(5)`, **unconditional** | hard-coded |
| `tau_bar` | `lognormal(log(median L_data), median L_prior_sigma)` | the resolved turnover centre |
| `theta_C` | `Beta(1, 9)` | SEASON-DERIVED |
| `z_f`, `z_c` | `student_t(4, 0, 1)` | `f_walk_df = 4` |
| `sigma_f`, `sigma_c` | half-normal(0, 1.5) | - |
| `cfi_kappa`, `cfc_kappa` | `lognormal(log(20), 0.75)` | log-SD hard-coded |

**Two naming traps in the data block.** The five `value_cauchyDF_*` names read as degrees of
freedom and are used as the **scale** of a Cauchy centred at zero; combined with the
`<lower=0>` declarations each is a **half-Cauchy(0, 1)**, not a Cauchy with one degree of
freedom. And `value_betashape_phi_*_scaled` is one value used for **both** Beta shapes, giving
the symmetric `Beta(2,2)`.

`sigma_IE` is placed unconditionally on purpose. It previously had no prior when `IE_n = 0`,
an improper flat direction that drifted to about 1e307 and was the boat's dominant divergence
source.

### 14.8 Generated quantities

| quantity | what it is |
|---|---|
| `E`, `E_sum` | crab-directed effort per day and its total: `lambda_E * E_scale * L * f` |
| `lambda_Ctot_S`, `C_expected`, `C_expected_sum` | expected catch per day and its total: `E * lambda_C * zi_scale` |
| `C`, `C_sum` | posterior predictive catch draws and their total |
| `L_out`, `tau_bar_out` | the fitted daily turnover and its shared level |
| `f_crab_out`, `f_lower_out`, `combo_c_out` | the per-stratum crabbing fraction, the implied crab-only share, the combo share |
| `sigma_f_out`, `cfi_kappa_out`, `sigma_c_out`, `cfc_kappa_out` | the walk and contact-model hyperparameters |
| `theta_C_out`, `zi_scale` | the zero-inflation probability and the reporting multiplier |
| `R_G_out`, `R_G_boat_out`, `kappa_OSP_out`, `sigma_IE_out` | the expansion ratios and scales |
| `B1_C_out`, `B2_C_out`, `B_open_out`, `gamma_C_out` | the CPUE effects |
| `log_lik_gear`, `log_lik_trailer`, `log_lik_osp`, `log_lik_catch` | pointwise log-likelihood for PSIS-LOO |

> **`C_sum` is NOT a full predictive distribution for the season total, and the reported
> `Predictive_Catch` interval should be read with that in mind.** The predictive draw is
> `C = poisson_rng(lambda_Ctot_S)`: Poisson on the expected rate. It therefore propagates all
> posterior parameter uncertainty but carries **neither** the catch likelihood's NB2
> overdispersion `r_C` **nor** the zero-inflation mixture, so it is under-dispersed relative
> to the fitted observation model. On the reference run this is why `Predictive_Catch`
> (94,386 [77,519, 118,776]) is barely wider than `Expected_Catch` (94,376 [77,566, 118,602]),
> when a true predictive interval would be materially wider. The headline interval this
> document quotes is the **expected-catch** interval, which is the right quantity for a
> seasonal harvest total; the predictive line should not be cited as a prediction interval
> for an observed catch. The Stan file does not comment on this.

### 14.9 Declared-but-inert parameters, and dead code

This section exists because it matters for cross-run comparison: several parameters are
sampled, reported, and mean nothing in a given run.

**Sampled from their prior, entering no likelihood.** Each of these is a genuine sampled
dimension whose posterior IS its prior, reported in the output as if it were an estimate.

| parameter | inert when | reported as |
|---|---|---|
| `kappa_OSP` | `osp_scale_is_tau = 1` (**production**) or no OSP days | `kappa_OSP_out` |
| `sigma_r_OSP` / `r_OSP` | no OSP days | `sigma_r_OSP` |
| `R_G_boat` | every shore fit | `R_G_boat_out` |
| `R_G` | every boat fit | `R_G_out` |
| `sigma_IE` | every boat fit, and any shore fit below the I/E day floor | `sigma_IE_out` |
| `gamma_C` | `estimate_cpue_density = 0` (**production**) | `gamma_C_out` |
| `B2`, `B2_C` | a window containing no holiday | `B2_C_out` |
| `Lcorr_E`, `Lcorr_C` | always, at `G*S = 1` | `Omega_E`, `Omega_C` (the constant 1) |

The unconditional priors on `R_G`, `R_G_boat` and `sigma_IE` are deliberate and must not be
moved inside a guard: a `real<lower=0>` with no prior is improper, which is what the
`sigma_IE` history above records. `model_diagnostics.R` flags a decoupled parameter so a
prior-only posterior is never read as an estimate.

**Zero-size when off**, so switching the feature off leaves the unconstrained parameter
vector unchanged and a fixed-seed rerun reproduces byte for byte: `theta_C`, `f_theta`,
`f_lower_param`, `z_f` / `sigma_f` / `cfi_kappa`, `z_c` / `sigma_c` / `cfc_kappa`,
`osp_f_kappa`, `tau_bar`, `B_open`, `L_raw`.

**Unreachable under the shipped configuration.** The entire legacy scalar-`f` construction
(`crab_fraction_dynamic = TRUE` makes it zero-size); the whole OSP crab-only block
(`use_osp_crab_lower = FALSE`); the `kappa_OSP` arm of the OSP mean
(`osp_scale_is_tau = TRUE`); and the `estimate_L = 0` branch, which the prep never selects.

**Dead code, kept and named so nobody re-derives it by accident.** The `Crab_*` census stream
(`Crab_n`, `day_Crab`, `section_Crab`, `Crab_I`, `p_I_crab`) is declared in the data block and
referenced nowhere else in the program; it contributes no likelihood term and cannot be
activated without editing the Stan file. `n_effort_obs` is retained for R-interface
compatibility after the overdispersion marginalization. `osp_f_n_crab` is carried for
reporting only.

**The boat I/E path needs a code change, not a flag.** The prep gates the I/E stream on
`is_shore`, so `IE_n = 0` for every boat fit. The Stan I/E likelihood as written uses
`lambda_E * L` with no `/ R_G_boat`, which would be the wrong mean for the boat, whose
`lambda_E` is gear rather than groups. Activating the boat I/E stream therefore requires
editing the likelihood, and `use_boat_ie` alone will not do it correctly.

## 15. The Point Estimator

One call per population x sub-season. The stratum is **(stat-week x day-type)** in production
(`period_pe = "week"`, three day types: weekday, weekend, holiday; `sections = c(1)` collapses
the section dimension).

**Effort.** Sampled days are averaged within a cell and expanded to that cell's calendar days:

```
ybar_h = mean( daily effort )  over the cell's sampled days
est_h  = ybar_h * N_h                                 N_h = the cell's calendar days
effort_total = SUM_h est_h
```

Daily effort matches the BSS's unit exactly, read from the same `bss_effort_spec()`:

```
boat:   mean_count(d) * gear_per_group * tau_boat * f_crab(d)
shore:  mean_count(d) * tau_shore                             (gear-deployments; no f)
```

**Catch.** The within-stratum catch rate is a **RATIO OF SUMS**, `sum(catch)/sum(gear)` over
the cell's sampled days, not a weighted mean of per-day ratios. That is not a stylistic
choice: with weekly strata and about 50% day coverage many cells rest on one or two sampled
days, and a day with very little sampled effort produces an extreme daily ratio that then
multiplies the full cell effort. Before the fix that made the shore PE's implied CPUE 2.9x
its own ratio-of-sums. Ratio-of-sums is also the estimator the BSS `lambda_C` converges to as
`r_C` grows, and it is what the gear-resolved PE uses, so the two tracks cannot drift.

```
est_catch_h = est_h * cpue_h        cpue_h = sum(catch) / sum(gear) over the cell's sampled days
catch_total = SUM_h est_catch_h
```

A hard internal-consistency stop guards the whole thing: the implied CPUE
(`catch_total / effort_total`) must sit within a factor of 2 of the interview ratio-of-sums,
or the run stops rather than reporting.

**Thin cells, which is where the interesting behaviour is.** A cell with calendar days but no
sampled day is retained with `n = 0`, and a cell resting on exactly one sampled day has no
sample SD. On 2024-25 that is not rare: **21 of 93 shore all-gear cells are unsampled and 37
rest on one day, covering 44% of the component's calendar days**; for the boat, 23 unsampled
and 34 singletons, 43% of days. Three levers price them.

| lever | shipped | what it decides |
|---|---|---|
| `pe_empty_effort_stratum` | `"local_day_type"` | the MEAN an unsampled cell is filled with |
| `pe_empty_stratum` | `"local"` | the CPUE an unsampled cell is expanded at |
| `pe_variance` | `"impute_aware"` | what an unsampled or singleton cell contributes to the SE |

**The mean and the spread are borrowed from DIFFERENT donor levels, on purpose.** Donors are
computed over sampled days at three nested levels: month x day-type, day-type, and
sub-season. The **mean** comes from the finest level with at least ONE sampled day; the
**spread** from the finest level with at least TWO. Coupling them (requiring two for both)
moves the point estimate, which is how the fill was first got wrong.

**The variance, case by case.** `Var(est_h) = N_h^2 * Var(ybar_h)`:

| case | `Var(ybar_h)` | rationale |
|---|---|---|
| `n_h >= 2` | `s_h^2 / n_h` | the cell's own sampled days |
| `n_h == 1` | `s_donor^2` | a collapsed stratum: the donor spread with divisor 1. `sd()` of one observation is NA, and replacing it with 0 is what made a singleton cell contribute its full point estimate and no variance at all |
| `n_h == 0`, filled | `s_donor^2 * (1/n_donor + 1)` | the donor mean's own error plus one between-cell deviation, which is what imputing a cell actually costs |
| `n_h == 0`, zeroed | `0`, and a **bias** is reported instead | the point estimate is 0 and the error is a bias, not a variance. An SE around a zero would imply the truth could be negative |

That last row is why the retired `"zero"` fill is retired: on 2024-25 it priced 45 of 289
shore all-gear days and 48 of 289 boat days at zero effort, with the empty cells sitting in
the high months, and reported no uncertainty for any of it. Under `"zero"` the shore PE sat
21.5% below the shore BSS; under the shipped `"local_day_type"` it sits **1.3% above**, which
is what the fill was adopted to achieve.

**Two things the PE does NOT do, both deliberate and both open.** It applies **no
finite-population correction anywhere**, although a cell's `n_h` sampled days are drawn
without replacement from its `N_h` calendar days, so `Var(ybar_h)` should carry `(1 - n/N)`.
Adding one would REDUCE the SE and would have mixed an opposite-signed change into the
variance work, making neither attributable. And **the PE produces no catch variance at all**,
only an effort SE. In the headline port total a PE component therefore enters as a constant;
the monthly-by-mode table borrows the effort relative SE as a lognormal multiplier, which
means the PE catch inherits a spread rather than having one. Both are tracked as
CHANGE_REGISTER D22. Note the inconsistency this creates with Section 16, where the charter
expansion DOES apply an FPC.

## 16. The commercial census and the charter expansion

This component is two different estimators under one heading, split on the design fact that
the census premise holds for the commercial vessels and not for the charter vessels.

**The commercial part is an exact census.** The frame is the daily vessel tally, on the days
it was taken. Under the shipped `census_expansion = "none"` the estimate is exactly

```
commercial_dung = SUM over tally days of ( commercial_tally(d) * mean crab per commercial vessel )
```

with unsampled calendar days contributing zero. That is not an omission: samplers are
scheduled on the days those vessels are confirmed to operate, so a window day with no tally
is a day with no operation. On 2024-25 that is **6,405 crab from 164 vessel-trips on 47 tally
days, with 23 unsampled days that had no operation**.

**The charter part is an expansion**, because charter trips are not fully sampled. The frame
is the charter trip roster, taking per day the larger of the roster's trips and the tally's
charter count, with any excess carried as an unattributed row. The expansion is stratified by
**vessel**, and the variance is simple-random-sampling with a finite-population correction:

```
est_v = N_v * mean_per_trip_v
var_v = N_v^2 * (1 - n_v / N_v) * s_v^2 / n_v
charter_dung = SUM_v est_v          charter_var = SUM_v var_v
```

A vessel with one interview takes the vessel mean with the pooled SD and divisor 1, which is
the same collapsed-stratum logic as the PE's singleton cell. A vessel never interviewed, and
the unattributed surplus row, take the pooled charter mean. Pooling to a single stratum is
what `charter_expansion = "pooled"` does, and it is the honest fallback: splitting `N` across
vessel rows while giving each the pooled `n` would apply the FPC stratum by stratum to a
sample that was never stratified, and understate the variance.

On 2024-25: **34 trips, 20 interviewed (59%), 1,286 crab observed, expanding to 2,133 crab
with an SE of 73**, including 514 crab on 8 days that had roster trips but no tally.

**What reaches the port interval.** `census_uncertainty` decides. Under the shipped
`"charter"` setting, the charter expansion variance is carried and the commercial census
enters as a constant: the port assembly draws the component from a normal centred on the
total with that SE, truncated below at the catch actually observed (the commercial census
plus the charter crab seen). The total variance is reported in `census_variance.csv` either
way.

**The residual this component reports but does not carry.** "A census without error" is a
statement about the COUNT, not about the catch. The commercial per-vessel catch MEAN is
itself a sample: 141 of 164 vessel-trips in 2024-25 (86%, SE 123, which is 1.9%) and 44 of 67
in 2025-26 (66%, SE 146, 5.5%). `census_uncertainty = "sampling"` carries it.

**When the frame is missing** the component is zero, not NA, and it says so: a `warning()`
naming the interview count and the window, a per-day table flagging every day as
`no frame (no vessel tally or roster this season)`, and a zero-row variance table with the
full schema so downstream binds do not break. This is the 2023-24 situation, where 33
Westport commercial-vessel interviews exist with no vessel count to expand to. A `warning()`
four hours into a knit is easy to miss; read the census section of the report before quoting
a total from a span that includes such a season.

## 17. The convergence gate and model adequacy

**The gate is the single authority on method selection**, and it asks one question: did the
sampler work? Four criteria, all of which must pass:

| criterion | threshold |
|---|---|
| `pass_rhat` | `max(R-hat)` over the catch and effort totals `< 1.01` |
| `pass_neff` | `min(n_eff)` over the same `> 400` |
| `pass_div_fraction` | divergent transitions `< 5%` of post-warmup draws |
| `pass_impact` | the divergence-impact test `< 0.10` posterior SD |

**The divergence-impact test is the one worth understanding.** Divergences are not
automatically fatal; what matters is whether they move the answer. The test is

```
impact = | median(all draws) - median(non-divergent draws) | / sd(all draws)
```

computed on the catch total and the effort total, in posterior-SD units, after Betancourt
(2017). A fit with a few hundred divergences that do not shift the posterior location passes;
a fit with fewer that do, fails. When there are no divergences the test returns zero rather
than NA, so a clean fit passes rather than being unjudgeable.

Reported beside the gate but **not gating**: the legacy level-distortion statistic (kept for
continuity; it scales with posterior width) and treedepth saturation, which raises a soft
warning above 5% and explicitly does not affect the verdict.

On the reference run every fit passed with room: worst divergence fraction 2.17% against the
5% backstop, treedepth saturation 0%, every R-hat within 1.0007, `n_eff` from 4,600 to
21,200 against the 400 floor, and divergence impact at most **0.003** posterior SD against
the 0.10 threshold.

**`bss_use_pe_for(b)` is the only correct way to ask whether a component should report PE**
in a totals section. It combines the pre-fit data-sufficiency flag with the post-fit gate
verdict, and its shape is fail-safe: a missing verdict reads as PE. Diagnostic sections
deliberately check the data-sufficiency flag alone, so a fitted-but-gate-failed component
keeps its per-fit diagnostics; do not "fix" those to use the combined helper.

**Model adequacy is a different question and is never allowed to gate.** The gate asks
whether the sampler worked; adequacy asks whether the model is carrying the data. A fit can
pass the gate and still be flagged. The statistics, from PSIS-LOO and randomized-PIT
posterior predictive checks:

| statistic | flag threshold | what it means |
|---|---|---|
| `p_loo_frac` | `> 0.25` | effective parameters as a fraction of observations: the fit is spending structure the data cannot support |
| `n_pareto_bad` | `> 0` | PSIS-LOO is unreliable for some observations |
| `pit_worst_bias` | `> 0.05` | the predictive distribution is systematically off-centre |
| `cov50_worst_dev` | `> 0.15` | the 50% intervals are the wrong width |
| `disp_neff_min` | `< 400` | a dispersion parameter is poorly sampled |

Turning any of these into a hard gate is a modelling decision that would need its own
justification and its own run; it is not a default anyone should inherit. One statistic,
`disp_scale_min`, is reported and deliberately NOT flagged, because flagging it would need a
defensible reference for "too small" and the half-Cauchy priors on those scales have no
finite SD, so the contraction diagnostic is undefined.

## 18. AR temporal resolution

**Resolution is an inference lever, not a tuning knob.** It is selected per fit from effort
data density, then capped per population.

```
daily    if coverage >= 0.25 and effort days >= 20
weekly   else if effort days per week >= 1.5 and weeks >= 3
monthly  otherwise
```

then coarsened by `ar_max_resolution` if the cap is finer. The cap **only ever coarsens**; a
cap finer than the data-driven pick is ignored. Shipped caps:

| | shore pot closure | shore all gear | private boat |
|---|---|---|---|
| pooled | biweekly | **weekly** | monthly |
| gear-resolved | biweekly | monthly | monthly |

`"biweekly"` is unreachable from the data-driven branch; it arrives only from a cap or from
`ar_force`.

**Why the shore all-gear cap is weekly, and why that decision needed a ladder.** The
data-driven selector picks daily for that component, and **all four rungs passed the
convergence gate**, so the gate could not choose between them. Adequacy could:

| rung | periods | catch | `p_loo`/`n_obs` | bad Pareto k | gear coverage_50 |
|---|---:|---:|---:|---:|---|
| daily | 289 | 20,898 | 0.352 | 41 | 0.701 (+7.1 SD) |
| weekly | 44 | 21,547 | 0.096 | 1 | 0.559 (+2.1 SD) |
| biweekly | 21 | 21,383 | - | - | - |
| monthly | 10 | 20,771 | - | - | - |

The four rungs span only 3.7% in catch, so this is a smaller lever on the number than it
looks. But daily was spending about one effective parameter per three observations and was
flagged miscalibrated, and weekly was not. **`elpd` favoured daily by 85.8 nats on the effort
stream and must not be believed**: that is 1.08 nats per additional effective parameter, and
41 of the daily fit's 311 gear observations had Pareto k above 0.7, so PSIS-LOO had failed
for 13% of the stream being cited. Only daily is an outlier; weekly, biweekly and monthly sit
within 0.5 sampling SD of each other, so the evidence does not choose among them and the
agreed rule does: **report the finest rung that PASSES the gate.**

**The escalation ladder** implements that rule per run rather than freezing it in config.
With `ar_escalate` on, a component starts at the finest rung and coarsens only on a gate
failure, stopping at the first pass. It costs one extra multi-hour fit per failed rung, which
is why it ships off; it is the tool for a new season, and every attempt is logged with its
gate verdict, its own estimate and interval, and its per-rung adequacy. Degenerate rungs are
pruned (a resolution that yields the same period count as a coarser one, or fewer than three
periods: a one- or two-period AR is not an AR).

---

## 19. Design decisions and their rationale

The decisions that would look arbitrary without the reason, each with the evidence that
settled it.

**Gear-deployments, not crabber-hours.** Catch is sub-linear in soak time for pot and trap
gear (measured saturation 0.13 to 0.27), so a time-denominated effort unit fails the
pipeline's own linearity test while gear-deployments passes it (`beta_h` covers 1). This is
re-measured every run. A "gear-deployment" is a piece of gear a crabber had in the water on
that trip; it is **not** a pot lift, repeat checks of one gear slot are not extra
deployments, and slot re-use across the day is carried separately by `tau`.

**One CPUE process, with gear-type catch apportioned afterwards.** The pooled model fits a
single latent catch rate and splits the total to gear types using Dirichlet-propagated
interview shares. The alternative, a genuine per-gear CPUE process, is what the gear-resolved
model is for, and it is the cross-check rather than the headline because the per-gear
likelihood is thinner and the two agree to 1.17% at the port.

**The AR(1) is on a period index, not on days.** That is what makes resolution a lever at
all, and it means the same Stan program serves every resolution.

**A beta-binomial per contact day, not a binomial on the stratum sum.** Stated in Section
14.3 and repeated here because it is the single most consequential choice in the `f` block: a
binomial on the sum would have returned a posterior SD near 0.005 and called it precision.

**`f` in generated quantities only.** It would have been easier to put the crabbing fraction
inside `lambda_E`. Keeping it out means the boat total is exactly linear in `f`, the CPUE
posterior is invariant to it, and the +11,963 it moved is attributable to it alone, which a
control rung then measured.

**The turnover level is shared and its day-to-day spread is fixed.** Two hundred
eighty-nine independent per-day draws reproduced their own prior centre to four decimal
places. Fixing the spread as data is what stops it trading off against the level on a series
where most days are unobserved.

**Zero-inflation on shore only, in the same run as an unmodified boat.** The boat fits are a
deliberate negative control, and they are bit-identical to a run with the feature off, which
is what makes them one.

**Adequacy beside the gate, never gating.** Two different questions deserve two different
answers, visibly. Collapsing them would let a well-sampled but overparameterised fit be
reported as clean, or a slightly ragged but adequate fit be discarded.

**Outputs are committed to git.** Past estimates are preserved as produced, so a claim about
what a run said can be checked rather than recalled.

**Validate by run, never by reasoning alone.** A change that looks inference-neutral on paper
can perturb the sampler geometry. Every adopted change in Section 6 was isolated in a ladder
rung against pre-set criteria, and an item is not "done" until a run confirms it. Changes that
NARROW reported uncertainty need explicit sign-off, because they move the headline interval.

## 20. Limitations

Ordered by how much they could move or break the number, which is not the order in which they
are easiest to measure.

**1. Moored private boats are outside the effort frame entirely, and this is unbounded.**
Private boats interviewed at the marina (34 in 2024-25, all crab-only) and at the docks (88
in 2024-25; 221 dock boat interviews in 2025-26) are moored at the floats. **Their catch
rates enter the CPUE; their effort is in neither the trailer count nor the OSP ramp total.**
The private boat component is half the port total, and this omission is one-directional:
correcting it can only add. It is not sized, and it cannot be sized from the desk. What would
size it is a Westport moorage count, slips occupied or a marina boat count, on enough days to
build a ratio to the ramp total. The 2022-24 effort protocol had a "Boats Entering Marina
Count" column, recorded at Ilwaco (2,209 in 2023-24) and never at Westport, so the field form
already supports it. **If a reviewer finds one thing in this method, it will be this.**

**2. The boat CPUE denominator is thin, and the boat is half the port.** Roughly 200 usable
private-boat interviews across the 2024-25 season, on about a quarter of the days. The OSP
stream fixed the boat's EFFORT problem and does nothing for its CPUE problem; no amount of
boat counting fixes a thin catch-rate sample, and no modelling choice manufactures
uncollected information. More boat interviews is the acknowledged limitation of this method.

**3. The derived shore turnover is a multi-season quantity presented as a per-window one.**
The shipped 2.477 pools all 40 I/E days in the workbook, spanning 2023-08 to 2026-08, of
which **6** fall inside the 2024-25 season. The window's own 6 days give **2.225**, 10.2%
lower, which is worth about **3,500 crab (3.7%)** on the port total. Two things should be said
about it together. The DIRECTION of the change from 1.7 is well supported: 1.7 was the wrong
quantity, and both candidates sit far above it. But the MAGNITUDE is not resolvable from the
I/E data that exists: with a between-day log-SD of 0.344, the gap between a 6-day subset and
the pooled 40 is about **0.83 standard errors**, and the window value sits **1.07 prior SDs**
from the shipped centre, so the shipped prior already covers it, which is why the fitted
posterior landed at 2.394 between the two. `tau_shore_derive_window_only = TRUE` prices the
alternative, and the field fix that would actually settle it is a paired gear count on every
I/E day. Tracked as CHANGE_REGISTER D24.

**4. The PE's unsampled-cell CPUE fill is not settled, and it moves the validation rather
than the estimate.** The boat PE sits 18.8% below the boat BSS under the shipped
`pe_empty_stratum = "local"` and 6.1% below under `"pooled"`, and the PE's own internal target
points the same way: the boat interview ratio-of-sums is 3.276 crab per deployment while the
PE's implied CPUE is 2.650 under `local` (0.81x) against 3.067 under `pooled` (0.94x), with
the BSS at 2.928 (0.89x). Against that, the theoretical case for `local` is sound: a
month-local effort fill multiplied by a season-pooled rate counts the seasonal gradient
twice, and a correct month-local fill SHOULD pull the effort-weighted implied CPUE below an
interview ratio-of-sums that is not weighted by the calendar. So the open question is whether
0.81x is too far. Because no component fell back to PE in the reference run, this affects the
cross-check that the estimate is judged by, not the estimate. Tracked as D19.

**5. The catch stream is under-covered at every AR resolution.** Randomized-PIT 50% coverage
runs 3.5 to 4.6 sampling SD low on the catch stream regardless of resolution, so it is a
likelihood problem, not a resolution problem. The zero-inflated mixture halved both count
bins and did not close the shore all-gear one bin, which is still 3.3 SD out: the data has
more zeros AND fewer ones than a zero-inflated NB can produce, i.e. it is more bimodal than
that mixture can be. A hurdle model, or a two-component NB with its own mean in the low
regime, is the shape that would fit it. **The remaining gain is bounded at roughly 6% of the
catch**, which is stated so the effort can be judged before it is spent. Tracked as D4.

**6. `Predictive_Catch` is not a prediction interval.** See the box in Section 14.8: the
predictive draw is Poisson on the expected rate and carries neither `r_C` nor the
zero-inflation, so it is under-dispersed relative to the fitted observation model.

**7. The annual PE-BSS agreement is not monthly agreement.** On the reference run the shore
all-gear annual totals agree to 1.8% while the PE over-allocates January to March by about
2.2x and under-allocates June and July, and the errors cancel. Do not cite the annual
agreement as month-by-month agreement. (The boat half of this comparison was additionally
distorted until 2026-09-12 by a reporting defect: the monthly share omitted the per-day
crabbing fraction, which was exact while `f` was a constant and stopped being exact when `f`
became a monthly walk. Fixed; the boat monthly comparison needs re-reading on the next
render.)

**8. The design-based layer is inconsistent on finite-population corrections, and the PE has
no catch variance.** Section 15 states both. Neither affects the headline, because the
headline is the BSS posterior; both affect what the PE column can be asked to support.

**9. The contact stream sees boats during sampler shifts only.** The shifts cover about 75%
of a day's boat returns, and the trip-type mix does not drift with the hour inside them on
2024-25, which is measured rather than assumed. The roughly 15% of returns after the last
check-out are unclassified. **Only OSP's all-day crabbing-only count can say whether their
mix differs**, and when it arrives the comparison is direct: `f_lower_out` is the model's
prediction of that column. Tracked as D11.

**10. The gear-resolved cross-check still differs in two ways.** Its shore all-gear cap is
monthly rather than weekly, and its Stan has no zero-inflation block, so the two tracks
differ in a resolution and in a likelihood. At a common resolution they agree on shore
all-gear to 0.08%. Closing the gap needs a gear-track ladder, about 3 hours. Tracked as D3
and D6.

**11. Pre-2024-25 gear labels are one option short.** There was no ring-net option on the
2022-23 or 2023-24 forms, which offered "Collapsible trap or ring" and, separately, "FISHING
ROD WITH foldable trap" and "Star trap", so a hand-held foldable trap had no home but the "or
ring" option. The builder implements the 2023-24 workbook's own `gear_key` mapping to "Ring
net", which is one-to-many in truth. Affects the gear-resolved 2022-24 fits only. Tracked as
D20.

## 21. Weather and tide covariates: evaluated and excluded

Tide phase and range, daytime high-tide timing, wind and wave height were screened as
candidate covariates on both effort and catch rate, after accounting for weekend and holiday
effects, and compared against the baseline with PSIS-LOO under a **pre-committed** margin.

**They are excluded.** No candidate cleared the margin, and the result is worth stating in
its own terms: the apparent covariate signal was **false precision**, an artefact of
comparing a covariate model's `elpd_loo` against the SE of one model total rather than
against the paired SE of the difference. The module that produced this conclusion is in
`06_diagnostics/`, its analysis is `WEATHER_COVARIATE_ANALYSIS.md`, and it is now **stale**:
its Stan fork is missing roughly 40 data variables the production model declares, so do not
cite any harvest number from it. Its conclusion, exclusion, stands.

The related **other-fishery opener** covariates are a separate mechanism and are also off:
they are built, tested and inert, and the razor-dig predecessor bought no predictive gain
(elpd within 1 SE) while moving the port +0.6%. One thing to know before ever switching a
BOAT opener covariate on: boat effort is a count of ALL private vessels, so an opener term
makes the latent effort track a surge of boats that are going fishing for something else, and
a crabbing fraction that is not opener-aware then converts that surge to crab effort at the
same rate as a closed-opener day. Pair the two or the bias gets worse, not better. CPUE is
deliberately offered no opener covariate at all: all eight opener-versus-catch-rate tests on
2024-25 came back null (p 0.29 to 0.90).

## 22. Glossary

| term | meaning |
|---|---|
| **BSS** | Bayesian state-space model: the time-series estimator, `crab_bss_pooled.stan` |
| **PE** | Point Estimator: the design-based stratified expansion, `run_pe_pooled()` |
| **the gate** | the convergence gate, which decides PE vs BSS per fit (Section 17) |
| **adequacy** | the model-adequacy statistics, reported beside the gate and never gating |
| **gear-deployment** | a piece of gear a crabber had in the water on that trip; the effort unit. Not a pot lift |
| **turnover, `tau`** | trips per gear slot per day, dimensionless. The daily expansion factor `L` |
| **`L_effective`** | the effective day length in hours, from the I/E regression. A DIAGNOSTIC under the production effort unit, not the expansion factor |
| **crabbing fraction, `f`** | the share of counted boats that were crabbing. A per-stratum logit random walk in Method v2.0 |
| **combo share, `c`** | among crabbing boats, the share whose trip also targeted another fishery |
| **`f_lower`** | under the dynamic `f`, the derived quantity `f(1-c)`: the model's prediction of OSP's crabbing-only column. Under the legacy construction, a hard lower bound on `f`. Same name, two meanings |
| **sub-season** | pot closure or all gear, split at the pot-open date. Each is an independent fit |
| **`ring_net_only`** | the internal key for the pot-closure sub-season, kept for filename continuity. Displayed as "Pot closure" |
| **stratum (PE)** | a (stat-week x day-type) cell |
| **stratum (`f`)** | a year-month, in production |
| **thin cell** | a PE stratum with one sampled day or none |
| **donor level** | the nested level (month x day-type, day-type, sub-season) an unsampled cell borrows from |
| **AR resolution** | the period index the AR(1) runs on: daily, weekly, biweekly or monthly |
| **the ladder** | the improvement ladder: a sequence of runs isolating one change per rung. Also the AR escalation ladder, which is a different thing |
| **rung** | one run of the improvement ladder, or one resolution attempt of the AR ladder |
| **`p_loo`** | the effective number of parameters from PSIS-LOO; read as a fraction of observations |
| **SEASON-DERIVED** | a `run_config.R` tag: this value was derived on 2024-25 and is a starting point on new data |

## 23. Development history (summary)

Newest first. The full version-by-version log with working notes is
`BSS-GH-pooled-CPUE-model-development-history.md`; every change on the current branch with
its status and evidence is `development_notes/CHANGE_REGISTER.md`.

- **Method v2.0 (2026-09-12).** The method of record becomes the model that runs: the dynamic
  monthly crabbing fraction from the coastal crab samplers' interview contacts, the OSP daily
  port count as a second boat effort stream, both turnovers data-derived, the zero-inflated
  shore catch likelihood, the weekly shore all-gear AR, the commercial census plus charter
  expansion, and the PE's month-local unsampled-cell fill. Method v1.0 is archived.
- **2026-09-11.** The first full improvement ladder ran. Port total 72,027 (superseded) to **94,376
  [77,566, 118,602]**, with the four movers summing to the whole within 24 crab, the `f`
  factorization proved, shore and boat bit-identity measured, and the gear-resolved
  cross-check at -1.17%.
- **2026-09-08 to 2026-09-10.** The dynamic `f` and the combo share `c`; both turnover priors
  resolved per window; the inputs rebuilt from the per-season creel workbooks; multi-season
  spans; the PE's unsampled-cell levers; the commercial/charter split.
- **2026-09-07 / 2026-09-08.** The weekly shore all-gear AR and the zero-inflated shore catch
  likelihood adopted together and gate-confirmed, the adoption render reproducing its
  candidate bit-identically across 10,253 shared parameter rows.
- **2026-09-01 / 2026-09-02.** The shared boat turnover `tau_bar` adopted. The gear-track boat
  sampler fix, which restored a cross-check that had never worked.
- **2026-08-25 to 2026-08-31.** The improvement batch: the shore I/E observation-unit fix,
  model adequacy reported beside the gate, the AR escalation ladder, opener covariates, the
  regression harness.
- **2026-07-31.** The OSP boat-count branch opens: the OSP second effort stream, the crabbing
  fraction as a scalar, and OSP-identifies-tau.
- **v7.5 to v7.9 (2026-07).** The incomplete-trip filter; the boat and then the shore moved
  onto the gear-deployment effort scale; the repository refactor and the shore PE completion
  fix; `run_config.R` restructured.
- **v6.0 to v7.4.** The post-critique arc: adaptive AR(1), `L_effective` from the I/E
  surveys, the weekend CPUE effect, a data-driven `R_G` prior, the convergence-debugging
  sequence (divergence gate, boat tuning, non-centred AR, marginalized NB, scale-aware
  gate), extended diagnostics and PSIS-LOO. **Method v1.0 = code v7.4.**
- **v1 to v5.** Single-population dock prototype through three populations and two
  sub-seasons.

## 24. References

Betancourt, M. (2017). *A Conceptual Introduction to Hamiltonian Monte Carlo.*
arXiv:1701.02434. The divergence-impact reasoning behind the gate.

Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley. Stratified expansion, the
finite-population correction, and collapsed strata (section 5A.12).

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., and Rubin, D. B. (2013).
*Bayesian Data Analysis*, 3rd ed. CRC Press.

Hansen, M. H., Hurwitz, W. N., and Madow, W. G. (1953). *Sample Survey Methods and Theory.*
Wiley. The donor-variance treatment of an imputed cell.

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model evaluation using
leave-one-out cross-validation and WAIC. *Statistics and Computing* 27, 1413-1432. PSIS-LOO
and the Pareto k diagnostic.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Bürkner, P.-C. (2021). Rank-
normalization, folding, and localization: an improved R-hat for assessing convergence of
MCMC. *Bayesian Analysis* 16(2), 667-718.

Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed. Springer, ch. 2.
Imputation variance.

**Upstream.** This pipeline is a derivative work of
[CreelEstimates](https://github.com/dfw-wa/CreelEstimates), the WDFW freshwater creel
estimation framework (GPL-3.0). The BSS and PE methodology, the project structure, and
substantial portions of the R and Stan code originate there and remain copyright their
authors; the adaptation to recreational Dungeness crab is by WDFW. See `NOTICE`.
