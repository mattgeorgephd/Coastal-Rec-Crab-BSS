# Marine hazard effort covariates: NWS Small Craft Advisories and the USCG bar-restriction tick

**Date:** 2026-09-25. **Status:** BUILT, INERT (`marine_hazard_mode = "off"` ships). **Register:** A30 (the method), B35 (the runner), B36 (the archive), D30 and D31 (the two open items). **Target branch:** `OSP-boat-count-incorporation` (delivered as a patch series on top of 3954a29).

**What was asked.** Incorporate Small Craft Advisories and bar restrictions into the estimation as a co-factor that can help with estimation, the way the razor-clam digs and the other fishery openers were incorporated: an auto option in the config that tests the correlation and applies it, with a toggle; and a batch runner that measures the impact of SCA and bar restrictions separately and in combination against the baseline.

**What was built, in one paragraph.** Two per-day covariates for the EFFORT process, entering both Stan models through the `K_open / X_open / B_open` block the opener covariates already use, so no Stan file changed and the shipped configuration builds Stan data identical to the pre-patch preps (measured, Section 6). `nws_sca_any` is 1 on a day when an NWS Small Craft Advisory or a higher warning was in effect for the Grays Harbor Bar or the coastal waters off Westport at any moment between 04:00 and 16:00 local time, read from an archive of the NWS's own VTEC events, so it is known on every day of the window. `bar_restriction` is the samplers' "Bar Restrictions" tick, observed on sampled days and imputed on the rest from the archived advisories. `marine_hazard_mode` is `off` / `auto` / `on` / `manual`; `auto` screens each candidate per population with a log-link GLM on the sampled days and BH over the marine family. `06_diagnostics/run_marine_hazard_batch_2026-09-25.R` is the ladder (baseline, SCA, bar, both, auto) with its decision rule stated before the run.

**What was not done, stated first.** The ladder has NOT been run (about 17 h of MCMC). Nothing here is evidence that the covariates improve the estimate; the sampled-day association is strong for the boat and absent for the shore, and Section 7 says why that is not the question that decides adoption.

---

## 1. Why the K_open block, and not a new term or the old weather fork

The pooled effort process is log-linear in its day covariates:

```
lambda_E[d] = exp( mu_E + omega_E[period[d]] + B1 w[d] + B2 holiday[d] + X_open[d] . B_open )
```

with `X_open` a `D x K_open` matrix the R prep hands over and `B_open ~ normal(0, 1)` per column (`value_normal_sigma_B_open = 1` in `prep_bss_crab_pooled.R`). That block was built for the other-fishery openers (A9): `opener_design_matrix()` takes a `selected` list plus an `extra` vector of column names, looks each up in `params$opener_flags`, drops any column with fewer than `opener_min_days` days on either side, and emits the matrix flat with its labels. The marine covariates are two more names in `extra` and two more columns in `opener_flags`. Every property the openers have, these inherit: identifiability guard, labels in `opener_covariates_<fit>.csv`, `B_open_out[k]` in the summaries, the `structural_params` table.

The alternative, a dedicated `B_sca` term, would have meant a Stan edit in two models, a new data variable, a new prior key, a new column in every summary writer, and a bit-identity proof for the edit when off. It buys nothing the block does not already give. The weather module that was removed (A29) was a separate Stan FORK that drifted 40 data variables behind production; that is the failure mode the block avoids.

One property of the block matters for these covariates in particular and is stated in `opener_design_matrix()`'s own terms: the columns are REAL, not integer. `X_open` is `vector[D*K_open]`, so an imputed bar-restriction day carries a probability in (0, 1), not a rounded 0/1, and `opener_min_days` counts only the observed 0 and 1 days when deciding identifiability. The harness pins both (section 76(f)).

## 2. The NWS flag

**Source.** `04_input_files/nws_marine_hazards.xlsx` (a live pull, 2026-09-25, 2008-01-01 to the pull day; see Section 8, item 5, for the transcription it replaced), one row per VTEC hazard SEGMENT for zones PZZ110 (Grays Harbor Bar) and PZZ156 (coastal waters from Point Grenville to Cape Shoalwater out 10 nm), with `start_utc` (when the hazard takes effect), `end_utc`, the product id (which carries the announcement time), `in_effect` (a segment cancelled before it began is kept with FALSE) and the requested pull window on every row. Built by `04_input_files/build_nws_marine_hazards.R` from the Iowa Environmental Mesonet's archive of NWS products (`vtec_events_byugc.py`, one call per zone per calendar year). It is NOT in `build_all_inputs.R`, which stays offline; a model run reads only the committed workbook.

**Definition.** `nws_sca_any[d] = 1` if any row with a code in `marine_hazard_codes` and `in_effect` overlaps the half-open interval `[04:00, 16:00)` of local day `d` in `America/Los_Angeles`. The codes are SC.Y (Small Craft Advisory), the three pre-2019-12-03 advisories NWS merged into it (RB.Y rough bar, SW.Y hazardous seas, SI.Y winds; Service Change Notice 19-83), and the warnings that supersede an SCA (GL.W gale, SR.W storm, SE.W hazardous seas, HF.W hurricane force). Watches, dense fog and special marine warnings are carried in the workbook and not counted. `nws_sca_bar` and `nws_sca_coastal` are the single-zone variants, offered as alternative candidates; at most one NWS definition enters a population's model, because the three are near-collinear.

**Time arithmetic, because it is where these things go wrong.** The window edges are computed as local wall-clock instants and converted to UTC per day, so 04:00 is 12Z in PST and 11Z in PDT; the harness checks both DST transitions of 2024. Overlap is half-open at both ends: an event ending exactly at 04:00 does not flag the day, and one starting exactly at 16:00 does not either. NWS winter products often start and end exactly on 04:00, 10:00 and 16:00 PST, so the edge rule matters and is written down. The R series was checked against an independent Python computation on every day of both windows (0 mismatches).

**Coverage is enforced.** The absence of a row is not evidence that nothing was in effect; it is evidence only within the window that was pulled. `marine_hazard_flag_series()` STOPS if the estimation window falls outside the pull window of either zone, naming the builder. A future season needs a re-pull before the mode is switched on.

**Why 04:00 to 16:00.** The launch decision is made before dawn and the counts are taken through early afternoon. On 2023-11-18 to 2026-09-08 (1,000 launch trailer counts; negative-binomial GLM after season, month and weekend/holiday) the 04:00-16:00, all-day and 06:00-14:00 definitions fit about equally (AIC within 4 of each other, all about 200 below the no-flag model; rate ratios 0.25, 0.28, 0.23) and the 04:00-10:00 and 10:00-16:00 halves fit worse (by 8 and 11). The choice is not sensitive within that family; it would be a mistake to take a half-day window.

## 3. The bar-restriction tick

**What it is.** A USCG closure of the bar to recreational vessels (33 CFR 165.1325), which no NWS product records. The samplers have a "Bar Restrictions" option in the survey's special-conditions field. In the earlier analysis of the samplers' records against the archive (project note of 2026-09-25) the SCA tick added nothing once the archive's timing was in the model (rate ratio 0.98), while the bar-restriction tick still did (0.67 [0.49, 0.90]); that is why the SCA comes from the archive and the restriction from the samplers.

**Observation rule.** On a sampled Grays Harbor day, `bar_restriction_obs = 1` if any survey that day ticked it, 0 if none did. Days before the option existed on the form are NOT observations of "no restriction": the field start is detected as the first day any survey carries a marine token (bar restriction, small craft advisory or gale warning), 2023-11-18 in the committed workbook, and `bar_restriction_field_start` overrides it. Blank special-conditions on a sampled day after that is read as no restriction, which is the same reading the SCA agreement analysis had to make and carries the same weakness (samplers under-record; the earlier note measured SCA sensitivity at 0.73 to 0.78 after 2023-24).

**Imputation on unsampled days.** `bar_restriction_impute = "nws"` fits, on EVERY observed day in the shifts workbook that the archive covers (716 days, 235 restricted), a logistic regression of the tick on the bar-zone SCA, the coastal SCA and a winter (October to March) indicator, all-day definitions:

```
logit P(restriction) = -2.25 + 1.29 * bar SCA + 1.38 * coastal SCA + 0.50 * winter
                       (SE 0.19, 0.24, 0.25, 0.21)
```

and assigns each unsampled day its fitted probability. On the 165 unsampled 2024-25 days that is a mean of 0.54 on advisory days and 0.15 on the rest. Under 30 usable days, no variation, a separated fit or a non-converged one falls back to the observed rate with a note (`"mean"` forces that). The fit is written to `marine_hazard_bar_imputation.csv`.

**The approximation this introduces, quantified.** The effort model applies `exp(B p)` on an imputed day, where the day is really a mixture with mean `p exp(B) + (1 - p)`. By Jensen the plug-in UNDERSTATES the mixture mean: by at most 3% at |B| = 0.5, 10% at |B| = 0.95 (the 2024-25 bar effect) and 17% at |B| = 1.24, worst near p = 0.6 and vanishing at p = 0 or 1. So on imputed days a bar term pulls interpolated effort slightly lower than a fully observed one would. The alternative, treating the tick as missing data inside Stan, would need a Bernoulli latent per unsampled day (marginalised, since Stan has no discrete parameters) and is more than the term is worth before the ladder has run. Recorded here rather than hidden.

## 4. The screen and the modes

**Modes.** `off` reads nothing, not even the archive, and installs an empty selection. `auto` builds the flags, runs the screen and selects. `on` enters every candidate (one NWS definition per population, the first in candidate order). `manual` enters exactly `marine_hazard_manual_shore` / `_boat`, ignoring unknown names with a note. Any mode but `off` stops on a missing or non-covering archive rather than fitting without the covariate that was asked for.

**The screen.** Per candidate and population, `glm(daily count ~ day_type + month + covariate, family = quasipoisson(link = "log"))` on the window's sampled days, the same adjustment set as `diagnose_fishery_spillover()` (the opener screen) and the same response, the daily sum of the counts. It departs from the opener screen in the LINK: that screen is an additive `lm`. The effort process is multiplicative and the advisory days sit in the low-count winter months, so an additive model attributes most of the difference to the month term. Measured on 2024-25 with the same days and adjustment set: boat SCA p 0.093 additive against 0.0006 log-link; bar 0.099 against 0.0003; shore 0.27 against 0.25. The reported `rate_ratio` is `exp(coefficient)`, the number the BSS would fit. The bar tick is screened on its OBSERVED days only.

**A weakness of the daily sum, reported rather than fixed.** A day carries one to three count sequences, and on 2024-25 the boat carried fewer on advisory days (1.06 against 1.45 per day), so part of a daily-sum difference is the number of counts. The BSS observation model is per count, so every screen row also carries the same GLM fitted per count (`rate_ratio_per_count`, `adj_p_per_count`, `n_counts`) and the report shows both. On 2024-25 they agree on every decision (boat SCA 0.291 against 0.290; bar 0.325 against 0.387; shore 0.972 against 0.905). Nothing selects on the per-count columns: a per-count GLM treats same-day counts as independent and overstates its precision.

**Selection under `auto`.** The multiplicity family is every candidate x population offered (3 as shipped), `n` pinned to the family so an unidentifiable test does not shrink it, Benjamini-Hochberg by default, threshold 0.05 on the adjusted p. If two NWS definitions clear the screen in one population, the smaller adjusted p is kept and the table says why. The selection is then subject to `opener_design_matrix()`'s identifiability rule in each fit's window.

## 5. What the 2024-25 data say (desk, no MCMC)

The window is 2024-09-16 to 2025-09-15, the R4 configuration. SCA-or-higher was in effect in the 04:00-16:00 window on 140 of 365 days (80 on the bar zone; every bar day is also a coastal day), concentrated in October to March (11 to 24 days a month) and rare in July to September (3, 2, 0). The bar tick was observed on 200 sampled days (86 restricted) and imputed on 165.

| Population | Covariate | Sampled days (flagged) | Raw mean, flagged / other | Rate ratio (daily sum) | Rate ratio (per count) | p raw | p BH | Selected |
|---|---|---|---|---|---|---|---|---|
| Boat trailers | SCA or higher, any zone | 186 (77) | 1.4 / 22.5 | **0.290** [0.15, 0.58] | 0.291 | 0.00056 | 0.00083 | yes |
| Boat trailers | Bar restriction (observed) | 186 (81) | 2.3 / 22.6 | **0.387** [0.23, 0.64] | 0.325 | 0.00032 | 0.00083 | yes |
| Shore gear | SCA or higher, any zone | 197 (82) | 53.0 / 101.0 | 0.905 [0.77, 1.07] | 0.972 | 0.25 | 0.25 | no |

The intervals are `exp(estimate +/- 1.96 SE)`. Three things to read out of the table rather than past it.

- **The raw boat difference (22.5 against 1.4 trailers) is mostly season.** Advisories are a winter phenomenon and so is the absence of boats; the month term takes most of it, and what remains is a factor of about 3. That factor is the covariate's claim.
- **The two boat flags overlap and are not the same thing.** Of the 86 observed restriction days, 63 had an SCA in the window; of the 114 observed non-restriction days, 21 did. Twenty-three restriction days had no advisory in the window and 21 advisory days had no restriction. That overlap is why the runner judges the bar term against the SCA rung (M4 against M2) rather than against the baseline, and why the two terms in one model will be posterior-correlated.
- **The shore does not respond within the season, though it does across seasons.** On 2023-11-18 to 2026-09-08 the dock gear counts give an SCA rate ratio of 0.74 [0.67, 0.82] after season, month and weekend (negative-binomial, 2,117 counts). Within 2024-25 alone, after day type and month, it is 0.905 with p 0.25. A season fit sees the season's data, so the screen does not select it, and `on` would enter it for the runner's M2 and M4 to measure.

## 6. What is proven about the patch as shipped

- `marine_hazard_mode = "off"` reads nothing: `marine_hazard_prepare()` returns before touching the archive and installs `marine_hazard_selected = list(shore = character(0), private_boat = character(0))`. The harness runs it with a non-existent archive path and it does not notice.
- Both preps then call `opener_design_matrix(..., extra = c(razor_extra, character(0)))`, which is `razor_extra`; the function takes `unique(c(selected, extra))`, so nothing can move.
- MEASURED on the real 2024-25 all-gear sub-season (D = 289), for the shore and the boat, for `prep_bss_crab_pooled()` and `prep_bss_crab_gear()`: the `stan_data` list built by the patched prep under `off` is `identical()` entry by entry, attributes included, to the list built by the 3954a29 prep from the same inputs, and to the 3954a29 prep fed a params list with every `marine_hazard_*` / `bar_restriction_*` key removed (the pre-patch `run_config` shape). K_open = 0 in every case. The pair is declared in `CODE_EQUIVALENT` with this measurement; M1 of the runner is the live bit-identity check against the R4 render.
- The runner's desk stage builds M4's Stan data against M1's for both all-gear fits and finds them different in `K_open` and `X_open_flat` ONLY (dry run, 16 PASS / 0 FAIL).
- No Stan file was touched; the harness asserts no Stan model mentions the covariates.
- Harness section 76: about 60 assertions over the time helpers, the flag series on a synthetic archive (definition, zones, codes, window, half-open edges, in_effect, coverage stop), the bar series (field start, observation rule, imputation, fallback, override), the screen (log-link recovers exp(-1.2) within 0.10; per-count columns), the four modes and the BH family, the hand-off with a fractional column, the preps' and drivers' wiring, the shipped defaults, the workbook, the builder's parser on the service's real JSON, and the runner's stated rule.

## 7. The runner and its rule

`06_diagnostics/run_marine_hazard_batch_2026-09-25.R`, pooled track, R4 configuration pinned key for key (`WINDOW`), one lever per rung:

| Stage | Mode | Shore | Boat | Purpose |
|---|---|---|---|---|
| M0 | desk | | | prerequisites; the Stan-data proof; no fit |
| M1 | off | | | the baseline; bit-identity against `20260910/pooled-CPUE-IMP-R4-shore-tau-newf` |
| M2 | manual | SCA | SCA | the advisory alone, both populations |
| M3 | manual | | bar | the restriction alone, boat |
| M4 | manual | SCA | SCA + bar | both |
| M5 | auto | (screen) | (screen) | what production would do with `marine_hazard_mode = "auto"` |

The rule, in the header before the run: (1) a rung is eligible only if every fit passes the convergence gate at the SAME AR resolution as M1; (2) adequacy must not degrade (shore and boat all-gear p_loo fraction <= 0.15, bad Pareto k <= 5% of the catch-stream n_obs); (3) a term is identified in a fit if the 95% interval of its `B_open_out` excludes 0, and an unidentified term is not adopted whatever the elpd says; (4) the PAIRED effort-stream `elpd_loo` against the control (gear stream for shore fits, trailer stream for boat fits, `loo_elpd_paired()`) must gain more than +2 paired SE on at least one fit where the term is active with no fit worse than -2 SE; between the two it is REVIEW, not FAIL; (5) the catch stream's paired elpd must sit within +/-2 SE, since an effort covariate has no business moving the catch fit; (6) the bar term is judged M4 against M2, not against M1, because on unsampled days it is driven by the archived advisories and would otherwise be credited with the SCA term's work; (7) the port total is reported at every rung and is not a criterion; (8) what the run cannot measure is stated up front (next section). Outputs: `05_output/marine_hazard_2026-09-25_{ladder,verdicts,recommendation}.csv`. Budget about 3.5 h per fitted rung from the R4 timings, 17 to 18 h for five; `STAGES` can drop rungs and M1 plus M2 is the minimum that answers anything. `DRY_RUN <- TRUE` ships.

**Adoption edit, if the run earns it:** `run_config.R` `marine_hazard_mode <- "auto"` (the screen decides each season) or `"manual"` with the adopted terms named. `auto` is the option that was asked for; `manual` is the option that a reviewer can read.

## 8. What this cannot show, and what would

1. **LOO scores the sampled days; the covariate exists for the unsampled ones (D31).** The boat is 79.8% extrapolated. A day covariate on the latent process earns its place by improving the interpolation on days with no count, and `elpd_loo` on the trailer stream can only score days that have one. Rule (4) therefore certifies "better on the sampled days", which is necessary and not sufficient. The diagnostic that would close the gap is a leave-one-WEEK-out block cross-validation on the effort stream: hold out each sampled week, refit (or approximate by PSIS with moment matching), and score the held-out counts under the AR alone against the AR plus the term. It does not exist. Until it does, adoption rests on the sampled-day evidence plus the mechanism, and the register says so.
2. **A boat term describes all private boats (D30).** Boat effort is trailers; `f` is monthly. If crabbing boats stay in port on an advisory day at a different rate from finfish boats, the crabbing share on those days differs from the month's and the monthly `f` cannot see it. The direction is not obvious (crabbers work pots near the entrance; finfish boats go further out). The egress classification (`boats_crabbing` / `boats_total`) on advisory and non-advisory days would settle it. The driver prints this note whenever a boat term is active.
3. **The bar tick on an unsampled day is an imputed probability**, so that column's unsampled-day contribution rests on the NWS flags and on the Jensen approximation of Section 3. If the ladder adopts the bar term, the honest reading is "the bar restriction as predicted by the advisories", and the archive alone might do the same work: M4 against M2 is where that shows.
4. **The OSP stream has no LOO stream.** `loo_summary_<fit>.csv` carries `catch`, `gear` and `trailer`; the boat's second effort stream (OSP port counts) is not scored, so a covariate's effect on the OSP likelihood goes unmeasured by rule (4). Adding an OSP pointwise log-likelihood to the generated quantities is a Stan edit and was kept out of this patch on purpose.
5. **The archive was a transcription for half a day (B36).** The environment that built the machinery could not reach the IEM host; the 617 rows were transcribed twice from the service's pages and agree row for row, and the first workbook was built from that CSV with `--from-csv`. The live pull that replaced it (590b5c3, 2008-01-01 to 2026-09-25, 4,387 rows) agrees with the transcription on 615 of 617 rows, the two differences being a product id and an event still in progress on the pull day; the 2024-25 flags, imputation fit and screen are byte-identical under the two. What remains true: `pull_end` is the pull day, so a window past it needs a re-pull before the reader will run.
6. **The screen is on one season.** `auto` re-screens every run on that run's window, which is what was asked for and also means the selected set can change from season to season. A `manual` selection is the reproducible one; the report table records which was used.
7. **Not measured: interactions with the openers.** `opener_covariate_mode` and `marine_hazard_mode` can both be on; the columns simply co-exist in `X_open`. No run has had both.

## 9. Files

- `03_R_functions/bss_marine_hazard_covariates.R` (new): the archive reader, the flag series, the bar series, the screen, the selection and the orchestration; header documents each.
- `03_R_functions/prep_bss_crab_pooled.R`, `prep_bss_crab_gear.R`: one line reading the selection into `extra`; two console labels.
- `01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd` (section 3.7, the decision chunk and the screen table), `BSS-GH-gear-type-CPUE-model.Rmd` (the decision chunk).
- `run_config.R` section 4.4b: the keys, shipped `off`.
- `04_input_files/build_nws_marine_hazards.R` (new), `04_input_files/nws_marine_hazards.xlsx` (new), `04_input_files/raw/nws_marine_hazards_iem_transcription_2026-09-25.csv` (new).
- `06_diagnostics/run_marine_hazard_batch_2026-09-25.R` (new); harness section 76; `run_improvements_2026-09-08.R` `CODE_EQUIVALENT` entry.
- Documentation: this note; CHANGE_REGISTER A30 / B35 / B36 / D30 / D31; PIPELINE_STATUS; method document 14.2 and 21a; `04_input_files/README.md`, `raw/README.md`, `README-R-functions.md`, `06_diagnostics/README.md`, `NEW_SEASON_GUIDE.md`, `CLAUDE.md`, `WEATHER_COVARIATE_ANALYSIS.md`, root README.
