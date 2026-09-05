# 07_documentation

Written documentation for the project: how the models work, what decisions were made and why, change history, and the forward-looking status and backlog. None of this is executed by a run; it is the reference layer.

For the one-paragraph project overview and quick start, see the [root README](../README.md).

## Start here

| File | Role |
|---|---|
| `NEW_SEASON_GUIDE.md` | **How to run the model on a new season, a part-season window, or a multi-season span**: the naive-run -> ladder -> pin-resolutions -> production workflow, the per-season config checklist, and the failure-mode table. The 2024-25 season was the development test season; this is the document for every season after it. |
| `development_notes/PIPELINE_STATUS.md` | **The single living status document**: current state, what is done, and the prioritized backlog. Read this first to see where the pipeline is. |

## Start here, part 2

| File | Role |
|---|---|
| `development_notes/CHANGE_REGISTER.md` | **Every change on the branch and its status** (ADOPTED / BUILT, INERT / OPEN / REJECTED / BLOCKED), the evidence, the effect on the number, the defects found and what each cost. The tabular companion to `PIPELINE_STATUS.md`. |

## Method documentation (Method v1.0, FROZEN)

These record **Method v1.0**, frozen against pooled code v7.4. The method description (estimators, likelihoods, gate criteria, expansion structure) is what they are for; **every number in them is historical**, and since 2026-09-07 the working model also differs from them in a likelihood (the zero-inflated shore catch) and a resolution (the weekly shore AR). Each carries a SUPERSEDED-NUMBERS banner at the top pointing at the authoritative run. Regenerating them against the current model is open documentation debt (CHANGE_REGISTER item D7).

| File | Describes |
|---|---|
| `BSS-GH-pooled-CPUE-model-documentation.md` | The pooled-CPUE production model, as frozen (Method v1.0 / code v7.4; code is now well past v7.9). |
| `BSS-GH-gear-type-CPUE-model-documentation.md` | The gear-resolved production model, as frozen (framework v5.6). |
| `BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md` | The weather-tide covariate module (stale; not production; its conclusion, exclusion, stands). |

## Development histories (the version-by-version change log)

| File | Describes |
|---|---|
| `BSS-GH-pooled-CPUE-model-development-history.md` | Full change log for the pooled model and its Stan file, newest first. |
| `BSS-GH-gear-type-CPUE-model-development-history.md` | Full change log for the gear-resolved model. |

The method documents summarize the history in one screen and point here for detail. `PIPELINE_STATUS.md` is the forward-looking backlog, not a changelog; these histories are the backward-looking record.

## Decision records and how-tos

| File | Content |
|---|---|
| `WEATHER_COVARIATE_ANALYSIS.md` | The finding that weather/tide covariates are excluded under the pre-committed PSIS-LOO margin (the false-precision result). Pairs with `06_diagnostics/`. Note: superseded in part by the deployment-scale move; the module itself is stale (see `PIPELINE_STATUS.md`, T2.4). |
| `effort_overdispersion_diagnostic_HOWTO.md` | How to read and run the effort-overdispersion diagnostic. |

## development_notes/

`PIPELINE_STATUS.md` (above) consolidated and superseded the generic historical working notes that used to live here. The folder now holds the living status document, the original critique, and the design and validation notes for the larger feature branches:

The 2026-08/09 validation campaign is a CHAIN of review documents, one per batch, newest first; each supersedes parts of the one before and `PIPELINE_STATUS.md` Sections 1b-1m are the running summary:

| File | Status |
|---|---|
| `adoption-review-2026-09-08.md` | **The current endpoint**: the adoption gate passed, `AD-A1-adopted` became the authoritative run, and the cross-track gap was shown to be a resolution difference (0.08% agreement at a common resolution). |
| `candidate-config-review-2026-09-07.md` | The completed 2x2 (AR resolution x catch likelihood); C1 identified as the candidate; the gear-track "monthly is also bad" figure retracted as a cross-model misattribution. |
| `ladder-zinb-review-2026-09-06.md` | The ladder's first real run (daily is overfitted; the estimate barely moves) and the ZINB re-scored under the corrected PPC. |
| `shore-ar-zi-review-2026-09-04.md` | Review of the 2026-09-03 batch: three one-line defects, five of eleven verdicts wrong, all corrected with an audit trail. |
| `validation-batch-review-2026-09-02.md` | Review of the 2026-09-01 validation batch: the shared-turnover adoption confirmed, the gear-track boat fit that had never worked, and the shore all-gear AR promoted to Tier 1. |
| `stage5-batch-review-2026-08-31.md` | Review of the 2026-08-30 Stage 5 batch: the tau x AR 2x2, and the retracted trailer over-coverage item. |
| `improvement-batch-review-2026-08-29.md` | Review of the 2026-08-27 improvement batch. |
| `improvement-plan-2026-08-27.md` | The sequenced follow-up plan (Stages 0-5), with status marks carried forward. |
| `ladder-validation-review-2026-08-27.md` | Review of the 2026-08-26 five-rung validation ladder. |
| `20260331-model-critique.docx` | Keep. The original external critique (primary source). |
| `PIPELINE_STATUS.md` | The single living status document (also linked under "Start here" above). |
| `GR-7-per-gear-CPUE-design.md` | Design note for the gear-resolved per-gear CPUE work. |

**OSP boat-count and crabbing-fraction design + validation notes.** `WBL-boat-count-plan.md` (plus its rendered `.html`), `phase1-osp-second-stream.md`, `phase1b-osp-gear-resolved.md`, `phase2-crab-fraction.md`, `phase2b-crab-fraction-pe.md`, `phase3-time-varying-f-and-osp-tau.md`, `osp-validation-review-2026-07-31.md`, and `osp_trailer_overlap.png`.

## Change register

`development_notes/CHANGE_REGISTER.md` lists every change made on the `OSP-boat-count-incorporation` branch with its status (ADOPTED / BUILT, INERT / OPEN / REJECTED / BLOCKED), the evidence for it, and its effect on the reported number, plus the defects found and what each cost. Start there for a view of where the work stands; `development_notes/PIPELINE_STATUS.md` is the narrative version with the run-by-run detail.
