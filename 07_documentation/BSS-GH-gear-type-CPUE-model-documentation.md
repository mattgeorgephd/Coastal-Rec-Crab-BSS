# Grays Harbor Recreational Dungeness Crab Harvest Estimation

## Gear-Resolved CPUE Model, framework v6.0

**Author:** Matthew George, Ph.D.
**Contact:** matthew.george@dfw.wa.gov
**Agency:** Washington Department of Fish and Wildlife (WDFW)
**Status:** Operational, **not published**. This is the **CROSS-CHECK** to the pooled model, not the headline estimator, and that is a design decision rather than a ranking of quality (Section 2).
**Framework version:** 6.0, adopted 2026-09-12, when the method of record moved to Method v2.0 and this track was brought onto the same configuration. Framework v5.6, frozen against the same era as Method v1.0, is archived at `archive/method-v1.0-gear-resolved-CPUE.md`.
**Reference run:** `05_output/20260928/gear-type-CPUE-model-2024-25`, rendered in the same `run_estimation.R` call as the pooled reference run (`model = "both"`, 2026-09-28, committed `ff750c4`): 98,588 [81,046, 123,563], **-1.29%** against the pooled 99,873, inside the 2% criterion. It superseded ladder rung R5 (`05_output/20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf`, 93,274, at R4's configuration and on the predictive catch).
**Convention:** no em dashes.

> ### READ THE POOLED DOCUMENT FIRST
>
> **This document describes only what is DIFFERENT.** The two models share the effort
> process, the crabbing fraction, the turnovers, the Point Estimator, the commercial census
> and charter expansion, the convergence gate, the adequacy reporting, the AR-resolution
> machinery, the input workbooks and the whole `03_R_functions/` library. All of that is
> specified in **`BSS-GH-pooled-CPUE-model-documentation.md`** (Method v2.0) and is not
> repeated here. What differs is how CPUE is modelled, and a short list of consequences.
>
> The authoritative run and its total live in one place, the box at the top of
> `development_notes/PIPELINE_STATUS.md`. The numbers below are from the reference run named
> above and may be superseded when you read this.

---

## 1. What this model is for

**Two questions the pooled model cannot answer as well.**

**Gear-type catch with its own uncertainty.** The pooled model fits one catch rate and splits
the total to gear types afterwards, using trip-level bootstrap interview shares (A17). That gives a
gear breakdown whose uncertainty is the uncertainty of the shares, not of the gear-specific
catch rates. With `gear_resolved_G = TRUE` the shore fits carry a genuine per-gear CPUE
process, so a per-gear catch estimate carries posterior uncertainty from the model.

**An independent check on the port total.** The two models share the effort side and differ
in the catch side, which makes their agreement informative. **Under the method of record** (the
re-render at B44 to B49, 2026-09-28, `run_estimation.R` with `model = "both"`) the gear track
reads **98,588 [81,046, 123,563]** against the pooled **99,873 [82,414, 124,438]**: **-1.29%**,
inside the pre-set 2% criterion (`05_output/20260928/cross_check_20260928_165651.csv`, PASS).
It is the first like-for-like check under the method of record: both tracks now sum the
expected catch (B44), and both carry the boat advisory term and the Float 17-21 fill. The
pooled boat all-gear turnover `tau_bar` is 3.118 in the mean. On R4's configuration (the
2026-09-12 method, before the boat term) the gear track had read 93,274 against the pooled
94,376, -1.17%, with the shared turnover agreeing to 0.02% and the monthly crabbing fraction to
0.002; that figure summed the predictive catch (see Section 3).

That agreement is the single most useful external validation the project has, because the two
implementations were written separately and reconcile through no shared catch code.

## 2. Why this is the cross-check and not the headline

Three reasons, none of which is "the pooled model is better".

**The per-gear likelihood is thinner.** Splitting the catch stream by gear splits the
interviews with it, so each gear's catch rate rests on fewer observations. That is the right
trade when you want per-gear catch and the wrong one when you want one seasonal total.

**It runs at `G = 1` by default, where the machinery is inert.** With `gear_resolved_G = FALSE`
(the shipped default) the gear split is PE-apportioned exactly as in the pooled model, so the
default gear-resolved run is a structurally simpler cross-check rather than a per-gear model.
**Do not raise `G > 1` without adding per-gear effort shares**: only gear 1 is observed in the
effort stream, so a naive `G > 1` would apportion the whole effort series to every gear.

**Its shore all-gear AR period is still monthly (`gear_period_bss`), and its zero-inflation
block ships off.** The block was ported from the pooled model on 2026-09-13 (D6) and is
switched by `catch_zi_tracks`, which ships `"pooled"`, so as shipped this track fits plain NB2.
So the two tracks currently differ in a resolution AND in a likelihood, which is why the port
gap is -1.29% under the method of record (-1.17% at R4's configuration) rather than smaller. At a COMMON resolution the two agree on shore all-gear to
**0.08%** (20,771 against 20,754), which is the number that shows the gap is a resolution
difference and not a disagreement. The gear-track ladder that was to close it ran on
2026-09-13/14 (`run_gear_ar_zi_2026-09-13.R`, campaign Section 1w): D6 says ADOPT the ZINB
(+11.3 nats at 2.29 paired SE, as on the pooled track), and D3 was re-scoped, because the
flat `gear_period_bss` moves the shore and boat all-gear fits together; at shore weekly, boat
monthly (the per-population form, `bss_gear_period()`) the tracks agree to +0.28% at the port.
Each is one ~35 minute render from closed. Tracked as CHANGE_REGISTER D3 and D6 (and D29, the
boat all-gear period that ladder exposed).

## 3. The reference run

`05_output/20260928/gear-type-CPUE-model-2024-25`, rendered in the same call as the pooled
reference run and on the same configuration (Matt, 2026-09-28, `ff750c4`; R 4.2.2, rstan 2.32.7).

| component | BSS median | PE | PE relative to BSS | AR resolution | divergences |
|---|---:|---:|---:|---|---:|
| shore, pot closure | 9,605 | 9,126 | -5.0% | biweekly | 15 |
| shore, all gear | 31,425 | 32,960 | +4.9% | monthly | 14 |
| private boat, pot closure | 1,287 | 1,192 | -7.4% | biweekly | 10 |
| private boat, all gear | 46,832 | 36,941 | -21.1% | monthly | **0** |
| commercial + charter | 8,538 | (the same) | n/a | not modelled | - |
| **port total** | **98,588 [81,046, 123,563]** | **88,758** | **-10.0%** | | |

Every fit passes the gate (at most 0.19% divergences; R-hat within 1.0000 on the summed catch
and effort). Against the pooled track, component by component: shore pot closure +1.5%, shore
all-gear -2.8% (monthly AR here against the pooled weekly, the D3 difference), boat pot closure
-6.2% (biweekly against monthly), boat all-gear -1.0%. The previous reference run,
`05_output/20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf` (93,274 [76,537, 117,227],
R4's configuration), predates B44, B45 and the boat term.

**Read the port total with B44 (2026-09-28).** Until B44 this track's port total, and every
component median above, summed the PREDICTIVE catch (`C_sum`: a Poisson draw per day on top
of the expected catch) under the row label `Catch`; the pooled headline is the EXPECTED catch
(`C_expected_sum`, `Expected_Catch`). The medians differ by little, the intervals by more. From
B44 this driver writes the pooled track's three rows (`Effort`, `Expected_Catch`,
`Predictive_Catch`) from the same quantities, so a cross-check rendered at or after B44 compares
like with like; the run above is rendered after it.

All four fits passed the gate. Sampler health is clean: every R-hat on the summed catch and
effort within 1.0000, `n_eff` from 7,981 to 10,618 against the 400 floor, treedepth saturation
0%, and the worst divergence fraction 0.19%, well inside the 5% backstop. Percentages are PE relative to BSS, as in the pooled
document.

**The boat all-gear fit sampling with zero divergences is worth noting**, because for most of
this project's history the gear track could not fit its boat component at all: it produced
554 divergences and a `tau_bar` effective sample size of 49. The 2026-09-02 fix moved the
per-fit sampler settings out of an experiment override and into the driver, and restored a
cross-check that had never worked.

## 4. How the CPUE process differs

Everything in the pooled document's Section 14.2 (the effort process), 14.3 (the crabbing
fraction), 14.4 (the turnovers) and 14.6 (the observation models, except catch) applies
unchanged. What changes is the catch side.

**Per-gear latent catch rate with shared AR(1) dynamics.** Instead of one `lambda_C_S`
process, there are `G` of them. They share the AR(1) temporal structure, which is what makes
the model estimable at all on a split interview stream: the gears are allowed to differ in
level and to move together in time, rather than each carrying an independent time series.

**An interview contributes to its own gear's likelihood.** The catch likelihood is indexed by
the interview's gear, so a ring-net interview informs the ring-net catch rate and not the pot
rate. At `G = 1` every interview lands in the single gear column and the likelihood is
identical in form to the pooled one.

**Per-sub-season regulatory gear exclusions are explicit.** The pot-closure sub-season
excludes pots by rule, so its `G` is one smaller than the all-gear sub-season's: with
`gear_resolved_G = TRUE`, shore all-gear runs `G = 5` and shore pot closure `G = 4`. The
exclusion is a config-declared regulatory fact, not something inferred from the interviews.

**Gear-share uncertainty (`gear_share_dirichlet`) is built and off.** It propagates
interview-share uncertainty into the per-gear intervals. It parses, it is byte-identical to
the previous behaviour when off, and it is forced off at `G = 1`. It has been sampled in
validation stages and never adopted, because at `G = 1` there is nothing for it to do.

**The holiday effort effect `B2`, and no separate holiday CPUE term.** The effort process
carries the same holiday term `B2` as the pooled model; the difference is on the CPUE side,
where this track folds holidays into `B1_C` and the pooled model has `B2_C`. As in the pooled model, note that the day-type indicators NEST: the weekend
indicator is 1 on weekends and on holidays, so `B2` is an increment on top of the weekend
effect rather than a separate level.

**The zero-inflation block is present and ships off.** Since 2026-09-13 (the D6 port)
`crab_bss_gear_resolved.stan` carries the pooled model's `zi_catch` / `theta_C` block line for
line. `catch_zi_tracks` decides which tracks fit it and ships `"pooled"`, so the gear prep
emits `zi_catch = 0` and the shore catch likelihood here is plain NB2 as shipped; the OFF path
is bit-identical to the pre-port model (11,021 parameter rows). This is the likelihood half of
the cross-track gap (-1.29% under the method of record), worth about -0.3%. Adding `"gear_resolved"` to `catch_zi_tracks` is the D6
adoption, which the evidence supports pending one render at the matched configuration.

**The shore gear-count expansion prior `R_G` is the pooled one (B46, 2026-09-28).** Until
B46 the gear Stan hard-coded `R_G ~ lognormal(log(1.3), 0.3)`, the 2024-25 interview ratio,
while the pooled Stan takes the centre from the season being fitted. The gear prep now passes
`R_G_prior_mu` and `R_G_prior_sigma` resolved exactly as the pooled prep does (the
`run_config.R` override, else the season's empirical interview ratio, about 1.28 on 2024-25,
else 1.3), so a season other than 2024-25 no longer inherits 2024-25's prior. The T1.3 sweep
(pooled document, PIPELINE_STATUS) found every component invariant to this prior.

## 5. Toggles this model reads that the pooled model ignores

All in `run_config.R` section 5.

| key | shipped | what it does |
|---|---|---|
| `gear_resolved_G` | `FALSE` | `TRUE` gives the SHORE fits a genuine per-gear CPUE process (all-gear `G = 5`, pot closure `G = 4`); the boat stays `G = 1`. `FALSE` apportions the gear split from interview shares, as the pooled model does |
| `gear_share_dirichlet` | `FALSE` | propagates gear-share uncertainty into the per-gear intervals. Forced off at `G = 1` |
| `alpha0_gear` | see config | the Dirichlet concentration for the gear shares |
| `gear_period_bss` | `list(all_gear = "month", pot_closure = "biweekly")` | the gear track's AR period: per sub-season (the flat form, shipped), or per population then per sub-season (`list(shore = list(...), private_boat = list(...))`, resolved by `bss_gear_period()`, B32). This, not `ar_max_resolution$gear_resolved`, sets the gear fits' periods |
| `catch_zi_tracks` | `"pooled"` | (section 2, read by both tracks) which tracks fit the zero-inflated catch likelihood; add `"gear_resolved"` to fit it here (D6) |
| `ar_adaptive` | `FALSE` | `FALSE` preserves the fixed per-sub-season period from `gear_period_bss` exactly (biweekly pot closure, monthly all gear; `period_bss` is the sub-season field it fills). `TRUE` hands the AR choice to the data-driven selector, which is inference-changing: validate first |
| `loo_effort_unit_comparison` | `FALSE` | `TRUE` restricts interviews to the common valid-denominator subset so a cross-unit `elpd_loo` comparison is legitimate. The comparison is done; `FALSE` for production |
| `use_boat_ie` | `TRUE` | use the WBL boat I/E ingress counts to identify the turnover once enough days exist. `IE_n = 0` is safe, and today there are 2 WBL days, so the stream is effectively absent |

**One caution on `ar_adaptive`.** The gear track's shipped resolutions come from the fixed
`gear_period_bss` path, not from the pooled track's `ar_max_resolution` ladder. Do not copy the
pooled caps across on the strength of symmetry alone. (The reason once given, that the gear
shore fits carry a thinner per-gear likelihood at `G = 5`, holds only with
`gear_resolved_G = TRUE`; as shipped the shore fits are `G = 1`. The question got its own
ladder on 2026-09-13/14, Section 1w, and D3 carries what it found.)

## 6. Running it

Identical to the pooled model, with one flag:

```sh
Rscript run_estimation.R --model gear_resolved   # this model alone
Rscript run_estimation.R --model both            # the pooled headline, then this cross-check
```

or set `model <- "gear_resolved"` (or `"both"`) in the RUN SELECTION block at the top of
`run_config.R`. An unknown flag stops the orchestrator rather than being ignored (B46). The
weather-tide module, which was pooled-only anyway, was removed on 2026-09-13 (A29).

**Run it AFTER a pooled run, on the same configuration, and compare the port totals** (the
`Expected_Catch` rows of the two `port_total_*.csv`; both tracks write the same rows since B44).
That is what it is for, and `--model both` does it in one call: each model renders in a fresh
environment and the orchestrator writes `cross_check_<timestamp>.csv` at the run-date level
(the two `Expected_Catch` medians, their relative difference, and PASS or REVIEW against
`cross_check_tolerance`, shipped 0.02). The pre-set criterion is agreement within 2%. If the gap exceeds it, the
first thing to check is whether the two tracks are at the same AR resolution, because that
has explained every gap so far.

Outputs land in `05_output/<YYYYMMDD>/gear-type-CPUE-model-<run_tag>/` and follow the same
catalog as the pooled model (pooled document Section 11), plus:

| file | contents |
|---|---|
| `catch_by_gear_type_detail.csv` | per-gear catch with posterior uncertainty when `G > 1`, PE-apportioned otherwise |
| the monthly / area / mode breakdowns | the gear track's own splits |

## 7. Version history, and a note on two numbering systems

Newest first. The full log is `BSS-GH-gear-type-CPUE-model-development-history.md`.

- **v6.0 (2026-09-12).** Brought onto Method v2.0: the dynamic crabbing fraction, both
  data-derived turnovers, the commercial census plus charter expansion, and the PE's
  month-local unsampled-cell fill, all of which are shared code and arrived with it. Ran as
  ladder rung R5 and reconciled to the pooled track at -1.17%. Still differs in the shore AR
  period and in the zero-inflation block, which was absent then and was ported 2026-09-13 but
  ships off (`catch_zi_tracks = "pooled"`).
- **2026-09-02.** The boat all-gear sampler settings moved from an experiment override into
  the driver. Divergences 554 to 2, `tau_bar` `n_eff` 49 to 19,292, gear port 51,385 to
  70,953. This restored a cross-check that had never worked.
- **v5.6 (2026-07-12).** `run_config.R` becomes the base parameter set, in parity with the
  pooled track's v7.9.
- **v5.5.** Both shore and boat on the gear-deployment effort unit;
  `loo_effort_unit_comparison` turned off for production.
- **v5.1 to v5.4.** Empirical gear proportions and the data-alignment work.
- **v5.0.** The initial gear-resolved release: per-gear CPUE processes, the `B2` holiday
  effect, the stratified census, the incomplete-trip filter, and the regulatory gear
  exclusions.
- **v1 to v4.** Shared with the pooled track; see the root `README.md`.

**The two numbering systems.** This track carries a v5.x / v6.x "framework" series in the
`.Rmd` header and in this document. That series is distinct from the Stan model's own
history: `crab_bss_gear_resolved.stan` carries no `vX.Y` tag, and the "Stan v3.2" label that
appears in older documentation is **stale**. The Stan file is materially past it: it carries
the pooled-track parity ports B1.3, B1.5, B1.6, B1.8 and B1.9, plus the run-driven fixes F1,
F2, P1 and P2. Read the framework tag for what the R pipeline does and the fix-marker section
of the development history for what the Stan model does; do not trust a "Stan vX.Y" label as
a statement of the current model state.

**On the shared library.** The two tracks share `bss_effort_spec.R` (the single source of the
effort unit and the matching I/E observation column), `bss_convergence_gate.R`,
`bss_ar_resolution.R`, `pe_effort_strata.R`, `estimate_comm_charter.R`, `crab_fraction.R` and
the rest of `03_R_functions/`, so they cannot drift onto different effort scales, gate
policies, PE stratifications or crabbing fractions. Functions that legitimately differ carry
a `_pooled` / `_gear` suffix. **Adding a track-specific helper without a suffix silently
overwrites the other track's function**, because every driver sources the whole folder.
