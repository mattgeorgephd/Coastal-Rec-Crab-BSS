# 07_documentation

Written documentation for the project: how the models work, what decisions were made and why, change history, and the forward-looking status and backlog. None of this is executed by a run; it is the reference layer.

For the one-paragraph project overview and quick start, see the [root README](../README.md).

## Start here

| File | Role |
|---|---|
| `development_notes/PIPELINE_STATUS.md` | **The single living status document**: current state, what is done, and the prioritized backlog. Read this first to see where the pipeline is. |

## Current model documentation (the method of record)

These track the live models in `02_stan_models/` and the drivers in `01_BSS_models/` and `06_diagnostics/`:

| File | Describes |
|---|---|
| `BSS-GH-pooled-CPUE-model-documentation.md` | The pooled-CPUE production model (pipeline code v7.9). |
| `BSS-GH-gear-type-CPUE-model-documentation.md` | The gear-resolved production model (framework v5.6). |
| `BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md` | The weather-tide covariate module (currently stale; not production). |

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
| `documentation_tables.xlsx` | Source tables that feed the write-ups. |

## development_notes/

`PIPELINE_STATUS.md` (above) consolidates and supersedes the historical notes in this folder. The historical notes are retained for provenance until confirmed for deletion (see the disposition table at the end of `PIPELINE_STATUS.md`):

| File | Status |
|---|---|
| `20260331-model-critique.docx` | Keep. The original external critique (primary source). |
| `20260710-OUTSTANDING_ISSUES.md` | Superseded by `PIPELINE_STATUS.md`. |
| `pipeline_state_review_20260709.md` | Superseded (F1-F5 done, P0-P3 carried forward). |
| `CODE_IMPROVEMENTS_REVIEW_v7.0.md` | Superseded (T1-T4 carried forward); detailed evidence appendix retained if wanted. |
| `PLANNED_IMPROVEMENTS.md` | Superseded (A/B/C/D backlog reconciled). |
| `ADDITIONAL_OUTPUTS_PROPOSAL.md` | Superseded (O1-O11 implemented). |
| `BOAT_RESOLUTION_EXPERIMENT.md` | Concluded (monthly AR chosen). |
| `CHANGES-2026-07-11.md` | Point-in-time refactor changelog (also recorded in the pooled development history). |
