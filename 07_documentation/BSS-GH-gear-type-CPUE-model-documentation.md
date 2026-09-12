# Grays Harbor Recreational Dungeness Crab Harvest Estimation

## Gear-Resolved CPUE Model, framework v6.0

**Author:** Matthew George, Ph.D.
**Contact:** matthew.george@dfw.wa.gov
**Agency:** Washington Department of Fish and Wildlife (WDFW)
**Status:** Operational, **not published**. This is the **CROSS-CHECK** to the pooled model, not the headline estimator, and that is a design decision rather than a ranking of quality (Section 2).
**Framework version:** 6.0, adopted 2026-09-12, when the method of record moved to Method v2.0 and this track was brought onto the same configuration. Framework v5.6, frozen against the same era as Method v1.0, is archived at `archive/method-v1.0-gear-resolved-CPUE.md`.
**Reference run:** `05_output/20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf`, ladder rung R5 on the pooled reference run's configuration.
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
the total to gear types afterwards, using Dirichlet-propagated interview shares. That gives a
gear breakdown whose uncertainty is the uncertainty of the shares, not of the gear-specific
catch rates. With `gear_resolved_G = TRUE` the shore fits carry a genuine per-gear CPUE
process, so a per-gear catch estimate carries posterior uncertainty from the model.

**An independent check on the port total.** The two models share the effort side and differ
in the catch side, which makes their agreement informative. On the reference configuration
the gear track reads **93,274 [76,537, 117,227]** against the pooled **94,376 [77,566,
118,602]**: **-1.17%**, inside the pre-set 2% criterion. The shared turnover agrees to 0.02%
across the two parameterizations and the monthly crabbing fraction to 0.002.

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

**Its shore all-gear AR cap is still monthly, and its Stan has no zero-inflation block.** So
the two tracks currently differ in a resolution AND in a likelihood, which is why the port
gap is -1.17% rather than smaller. At a COMMON resolution the two agree on shore all-gear to
**0.08%** (20,771 against 20,754), which is the number that shows the gap is a resolution
difference and not a disagreement. Closing it needs a gear-track AR ladder, about 3 hours of
fitting. Tracked as CHANGE_REGISTER D3 and D6.

## 3. The reference run

`05_output/20260911/gear-type-CPUE-model-IMP-R5-gear-crosscheck-newf`, on the same
configuration as the pooled reference run.

| component | BSS median | PE | PE relative to BSS | AR resolution | divergences |
|---|---:|---:|---:|---|---:|
| shore, pot closure | 9,036 | 8,591 | -4.9% | biweekly | 19 |
| shore, all gear | 28,334 | 29,737 | +4.9% | monthly | 8 |
| private boat, pot closure | 1,295 | 1,192 | -8.0% | biweekly | 24 |
| private boat, all gear | 45,374 | 37,018 | -18.4% | monthly | **0** |
| commercial + charter | 8,538 | (the same) | n/a | not modelled | - |
| **port total** | **93,274 [76,537, 117,227]** | **85,076** | **-8.8%** | | |

All four fits passed the gate. Sampler health is clean: every R-hat within 1.0011, `n_eff`
from 8,008 to 11,715 against the 400 floor, treedepth saturation 0%, and the worst divergence
fraction well inside the 5% backstop. Percentages are PE relative to BSS, as in the pooled
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

**A separate holiday effort effect `B2`.** The gear-resolved effort process carries its own
holiday term. As in the pooled model, note that the day-type indicators NEST: the weekend
indicator is 1 on weekends and on holidays, so `B2` is an increment on top of the weekend
effect rather than a separate level.

**No zero-inflation block.** `crab_bss_gear_resolved.stan` has no `theta_C`, so
`estimate_catch_zi` is inert on this track and the shore catch likelihood is plain NB2. This
is the likelihood half of the -1.17% gap, worth about -0.3%. Restoring symmetry is a Stan
edit plus a recompile (D6).

## 5. Toggles this model reads that the pooled model ignores

All in `run_config.R` section 5.

| key | shipped | what it does |
|---|---|---|
| `gear_resolved_G` | `FALSE` | `TRUE` gives the SHORE fits a genuine per-gear CPUE process (all-gear `G = 5`, pot closure `G = 4`); the boat stays `G = 1`. `FALSE` apportions the gear split from interview shares, as the pooled model does |
| `gear_share_dirichlet` | `FALSE` | propagates gear-share uncertainty into the per-gear intervals. Forced off at `G = 1` |
| `alpha0_gear` | see config | the Dirichlet concentration for the gear shares |
| `ar_adaptive` | `FALSE` | `FALSE` preserves the fixed per-sub-season `period_bss` exactly (biweekly pot closure, monthly all gear). `TRUE` hands the AR choice to the data-driven selector, which is inference-changing: validate first |
| `loo_effort_unit_comparison` | `FALSE` | `TRUE` restricts interviews to the common valid-denominator subset so a cross-unit `elpd_loo` comparison is legitimate. The comparison is done; `FALSE` for production |
| `use_boat_ie` | `TRUE` | use the WBL boat I/E ingress counts to identify the turnover once enough days exist. `IE_n = 0` is safe, and today there are 2 WBL days, so the stream is effectively absent |

**One caution on `ar_adaptive`.** The gear track's shipped resolutions come from the fixed
`period_bss` path, not from the pooled track's `ar_max_resolution` ladder. Do not copy the
pooled caps across: the gear shore fits carry a thinner per-gear likelihood and may genuinely
need a coarser AR. That question needs its own ladder, which the machinery now supports with
per-rung adequacy.

## 6. Running it

Identical to the pooled model, with one flag:

```sh
Rscript run_estimation.R --model gear_resolved
```

or set `model <- "gear_resolved"` in the RUN SELECTION block at the top of `run_config.R`.
The weather-tide module is **not** valid with this model; the orchestrator hard-stops early
if you ask for both, before any multi-hour fit.

**Run it AFTER a pooled run, on the same configuration, and compare the port totals.** That
is what it is for. The pre-set criterion is agreement within 2%. If the gap exceeds it, the
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
  cap and the absent zero-inflation block.
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
