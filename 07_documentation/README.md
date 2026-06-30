# 07_documentation

Written documentation for the project: how the models work, what decisions were made and why, change history, proposals, and the rendered documentation site. None of this is executed by a run; it is the reference layer.

For the one-paragraph project overview, see the [root README](../README.md).

## Current model documentation

These track the three live models in `02_stan_models/` and the drivers in `01_BSS_models/` and `06_diagnostics/`:

| File | Describes |
|---|---|
| `BSS-GH-pooled-CPUE-model-documentation.md` | The pooled-CPUE production model. |
| `BSS-GH-gear-type-CPUE-model-documentation.md` | The gear-resolved production model. |
| `BSS-GH-pooled-CPUE-weather-tide-covariates-documentation.md` | The weather-tide covariate module. |

## Rendered documentation site

| File | Role |
|---|---|
| `docs_index.Rmd` / `docs_index.html` | Landing page of the built doc site. |
| `docs_equations.Rmd` / `docs_equations.html` | Model equations. |
| `documentation_tables.xlsx` | Source tables that feed the write-ups. |

The `.html` files are the built output of the `.Rmd` sources; re-knit the `.Rmd` to regenerate them. These `.Rmd` files do not read the renumbered stage folders, so the folder reorganization did not require any edits to them.

## Decision records and how-tos

| File | Content |
|---|---|
| `WEATHER_COVARIATE_ANALYSIS.md` | The conclusion that weather/tide covariates are excluded for all three components under the pre-committed 4.0-SE PSIS-LOO margin (the false-precision finding). Pairs with `06_diagnostics/`. |
| `effort_overdispersion_diagnostic_HOWTO.md` | How to read and run the effort-overdispersion diagnostic. |
| `diagnostics_and_reproducibility_notes.md` | Notes on diagnostics and reproducing runs. |

## Change logs and experiments

| File | Content |
|---|---|
| `B1.5_change_notes.md`, `B1.6_change_notes.md` | Versioned change notes for those model revisions. |
| `BOAT_RESOLUTION_EXPERIMENT.md` | The daily-vs-weekly AR temporal-resolution experiment for boat effort. |

## Proposals and reviews

| File | Content |
|---|---|
| `PLANNED_IMPROVEMENTS.md` | Roadmap of intended changes. |
| `ADDITIONAL_OUTPUTS_PROPOSAL.md` | Proposed additions to the output set. |
| `CODE_IMPROVEMENTS_REVIEW_v7.0.md` | Internal code-quality review. |
| `20260331-model-critique.docx` | External model critique (Word). |

## Creel-lineage carryover

This project was forked from the WDFW freshwater-creel framework, and some of that documentation rode along. It is kept for reference but does not describe the current crab models:

| Item | Note |
|---|---|
| `Instructions for using Creel Estimates.docx` | Freshwater-creel user guide. |
| `Instructions for using the Creel Schedule Generator.docx` | Freshwater-creel scheduler guide. |
| `DRAFT_FreshwaterCreel rep_file sub structure mockup.docx` | Draft report-file structure mockup from the creel framework. |
| `FWC_bss_docs/` | Legacy reference subfolder (see below). |

### `FWC_bss_docs/`

Older reference material predating this project, retained for lineage:

- `02_Creel_Models_2022-01-20.csv` (2022 creel model list)
- `Skagit creel model - definitions_updated 2021-01-12.xlsx` (Skagit model definitions)
- `overview of time-series model_2019-02-28.pptx` (2019 time-series model overview)

These describe the antecedent state-space creel models the current BSS approach grew out of; treat them as historical context, not as documentation of the code in this repo.
