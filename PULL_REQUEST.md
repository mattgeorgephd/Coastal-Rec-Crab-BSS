# Merge `OSP-boat-count-incorporation` into `main`: Method v2.0 becomes the method of record

**Branch:** `OSP-boat-count-incorporation` → `main`
**Head:** `7e59ca4` · **Merge base:** `724eead` (2026-07-28, "WBL count data")
**Scope:** 85 commits · 135 files outside `05_output/` · +44,164 / −5,001 lines · 75 files added, 8 removed
**`main` has not moved since the merge base** (0 commits), so this merges without conflict.

> **What this changes for anyone reading `main`.** `main` today is the pre-FW-creel-meeting,
> pre-OSP pipeline: **Method v1.0**, one effort stream, literal turnover constants
> (`tau_boat_prior_mu = 1.2`, `tau_shore_prior_mu = 1.7`), a constant crabbing fraction, and an
> imputed commercial/charter census. It has no `osp_scale_is_tau`, no `shared_tau`, no
> `crab_fraction_strata` and no `census_expansion` key at all. After this merge the method of
> record is **Method v2.0**, specified in
> `07_documentation/BSS-GH-pooled-CPUE-model-documentation.md`, and it describes what the code
> actually runs rather than a frozen earlier version.

---

## 1. The number, before and after

`main`'s authoritative run is `05_output/20260715/pooled-CPUE-230256` (Run 6, 2026-07-15).
This branch's is `05_output/20260910/pooled-CPUE-IMP-R4-shore-tau-newf`, rung **R4** of the
improvement ladder, rendered overnight 2026-09-10 / 2026-09-11.

| 2024-25 Dungeness kept | `main` (Run 6) | this branch (R4) | change |
|---|---|---|---|
| Shore, pot closure | 6,275 BSS | **8,963** BSS | +2,688 |
| Shore, all gear | 20,608 BSS | **29,210** BSS | +8,602 |
| Private boat, pot closure | 1,170 **PE** (gate failed) | **1,372** BSS | +202, and it now fits |
| Private boat, all gear | 43,221 BSS | **45,604** BSS | +2,383 |
| Commercial / charter | 11,986 imputed census | **8,538** census + expansion | **−3,448** |
| **Port total** | **83,488** [70,866, 102,897] | **94,376** [77,566, 118,602] | **+10,888 (+13.0%)** |
| Port total, PE | 71,157 | 85,076 | +13,919 |
| Effort (gear-deployments) | 44,484 | 58,557 | +14,073 |

Two features of that table matter more than the headline.

**Every BSS fit now passes the gate.** On `main`, four of five components reported BSS and the
boat pot-closure component fell back to its PE point. Here **4 of 4 fits report BSS**, with the
worst R-hat 1.0007, the smallest relevant n_eff 4,604, the largest divergent fraction 1.78%
against a 5% backstop, and 0% treedepth saturation on every fit.

**The census went DOWN, and that is a correction, not a loss.** `main`'s 11,986 imputed
commercial/charter catch over unsampled days. The branch splits the component: commercial is an
**exact census** over the 47 tally days (6,405 crab from 164 vessel-trips; the 23 unsampled days
are recorded as *no operation*, which is what the tally says happened, rather than imputed), and
charter is an **expansion over the vessel roster** (34 trips, 20 interviewed, 1,286 crab
observed → 2,133 crab, SRS-with-FPC SE 73). Only the charter SE is carried into the port
interval (`census_uncertainty = "charter"`), because the commercial side is a census and has no
sampling variance.

## 2. What moves the number, attributed by its own ladder

The branch does not ask anyone to take +10,888 on trust. `06_diagnostics/run_improvements_2026-09-08.R`
fits the same season on successive rungs, each differing from the previous one only in declared
keys, with comparability proven before the MCMC rather than argued after it:

| rung | shore pc | shore ag | boat pc | boat ag | census | port | 95% CI |
|---|---|---|---|---|---|---|---|
| R1 pre-patch turnover priors + unit-aware filters + dynamic `f` | 6,343 | 21,502 | 1,070 | 43,668 | 7,884 | 80,715 | 66,042 to 103,228 |
| R2 + the calibration `tau_bar` prior | 6,343 | 21,502 | 1,372 | 45,604 | 7,884 | 83,336 | 67,624 to 107,227 |
| R2f control: R2 with the `f` block rolled back | 6,343 | 21,502 | 1,290 | 33,641 | 7,884 | 71,179 | 50,627 to 104,722 |
| **R4 + derived shore turnover + the census split (shipped)** | **8,963** | **29,210** | **1,372** | **45,604** | **8,538** | **94,376** | **77,566 to 118,602** |

R2 → R2f is the control that isolates the dynamic crabbing fraction: rolling `f` back moves the
boat all-gear component from 45,604 to 33,641 and nothing else, which is the design guarantee
holding exactly, because `f` enters only the generated quantities and the boat total is linear in
it. R2 → R4 is the derived shore turnover plus the census split, and it is the largest single
step (+11,040).

**The gear-resolved track is an independent cross-check, not a second opinion to average.** Rung
R5 re-estimates the same season with a per-gear CPUE structure and returns
**93,274** [76,537, 117,227], which is **1.17% below** the pooled R4. Two model structures on
the same data landing within 1.2% is the strongest corroboration the project has produced.

## 3. What changed, by mechanism

**A second effort stream, and what it identifies.** OSP total port boat counts
(`WBL_boat_counts.xlsx`) enter as a second observation of boat effort alongside the trailer
counts. With `osp_scale_is_tau = 1` the OSP series identifies the boat **deployment turnover**
`tau_bar` rather than a free scale, so the two streams constrain one quantity instead of
competing. `tau_boat_prior_mu` is no longer the literal 1.2: it is `"calibration"`, resolving to
a prior centre of **3.0300** derived from 61 paired OSP/trailer overlap days, and the fitted
`tau_bar` is **2.977**. The retired 1.2 and 1.7 survive only as `%||%` fallbacks and as the
ladder's R1 rung, which is what makes R1 a measurement rather than a memory.

**The crabbing fraction `f` became dynamic.** It was a constant 0.30. It is now a per-stratum
(year-month) logit random walk with Brownian `sqrt(gap)` steps, `z_f ~ student_t(4, 0, 1)`,
`sigma_f ~ half-normal(0, 1.5)`, observed by `cfi_crab[i] ~ beta_binomial(cfi_total[i], f *
cfi_kappa, (1-f) * cfi_kappa)` at one observation per contact day, from the **coastal crab
samplers' interview contacts** (not the blank I/E columns). The combined-trip share `c` walks on
its own. This is the R2 → R2f delta above.

**The shore turnover stopped being a guess.** `tau_shore_prior_mu = "derived"` resolves to
**2.4771** from the I/E time column at the count hours, with the prior log-SD floored at 0.10.

**The commercial/charter component became a census plus an expansion**, as described in §1.

**Effort is gear-deployments throughout**, both tracks, so turnovers are dimensionless and the
PE expands on the same turnover the BSS fits. No shore expansion anywhere is in hours.

**Multi-season spans are supported** (`pot_closures` takes one closure window per season). The
shipped configuration is the single 2024-25 season; the 2023-25 span is staged as a commented
block and **blocked on data**, because there is no 2023-24 vessel tally or charter roster (D8).

**The weather-tide covariate module is removed** (A29). The FWC advised against its use and
weather alone was not informative. Its Stan fork, its driver and the `run_weather` toggle are
gone; `WEATHER_COVARIATE_ANALYSIS.md` keeps the exclusion finding live, because the finding is
the reason, and the module document is archived under a banner.

## 4. Evidence a reviewer can run

| check | how | current result |
|---|---|---|
| Regression harness | `Rscript 06_diagnostics/test_improvements_2026-08-25.R` | **964 assertions, 0 failing**, no rstan needed, seconds |
| Shipped config reproduces the authoritative run's inputs | desk-checked 2026-09-12 | sub-seasons `ring_net_only [2024-09-16..2024-11-30]` / `all_gear [2024-12-01..2025-09-15]`, tau_shore 2.4771, tau_boat 3.0300, census 6,405 + 2,133 = 8,538 SE 73, **zero warnings** |
| Convergence gate | `convergence_report.csv` in the run folder | 4/4 BSS; worst R-hat 1.0007, min n_eff 4,604, max divergent fraction 1.78% vs the 5% backstop |
| PE / BSS agreement | `pe_vs_bss_comparison.csv` | shore within 2 to 4%; boat all-gear PE runs 19% below BSS, which is the turnover and `f` treatment and is expected |
| Cross-model | ladder rung R5 | gear-resolved 93,274, **−1.17%** from pooled R4 |
| Every rung comparable | `improvements_2026-09-08_manifest-newf.csv` | preflight FAILS on any non-delta key differing between rungs; `WINDOW` pins all nine per-season keys, the three PE unsampled-cell levers and the sampler seed; a three-layer code fingerprint is recorded per rung and `V1cross()` downgrades any verdict rendered under mismatched code |

The harness is the fastest way in. It asserts the **shipped defaults** in `run_config.R`, not
just that the functions work, which is the guard that stops a behaviour-changing toggle from
drifting on unnoticed.

## 5. Open decisions this merge carries into `main`

None of these blocks the merge; all of them are recorded with what they are worth, so `main`
inherits the questions rather than a silence.

| id | question | what it is worth |
|---|---|---|
| **D14** | Moored boats are outside the effort frame entirely. The trailer count counts trailers; a boat kept in a slip leaves no trailer. | **The largest unbounded risk in the estimate.** Needs a Westport moorage count. Cannot be bounded from existing data. |
| **D2** | The fraction of boats OSP records as *crabbing-only*. The column is built and inert, awaiting the OSP field. | Would replace part of what `f` currently infers with a direct observation. |
| **D24** | The derived shore turnover 2.4771 is a multi-season quantity used for one window; this window's own I/E data put it ~10% lower. | ~3,500 crab. **The existing data cannot settle it**: the gap is 0.1073 in log space against an approximate SE of 0.1295, **z = 0.83**, and the window value sits 1.07 prior SDs from the shipped centre. A window-only refit would move the number without evidence that it should. |
| **D19** | The PE's unsampled-cell fill. Shipped as `local_day_type` / `local` / `impute_aware`. | Moves the PE cross-check by up to 18% (72,224 / 81,160 / 90,861 / 85,076 across the four arms); **does not move the BSS estimate**. |
| **D3 / D6** | The gear track's monthly shore cap and its absent ZI block. | ~3 h of ladder to settle. Affects only the cross-check. |
| **D8 / D18** | The 2023-25 span (no 2023-24 tally or roster) and OSP counts stopping 2025-10-18. | Blocks the multi-season run and a 2025-26 boat fit on OSP. Data requests, not code. |
| open | Should `estimate_comm_charter()` **stop** rather than warn when the census frame is missing? | Currently warns and returns 0, which is how a silently-zero census component could ship. |

## 6. Not in scope, and stated plainly

- **WDFW has published no estimate from this pipeline.** Every "adopted" on this branch means
  adopted into the working model, not released. There is no published figure any change here has
  to stay consistent with.
- **The gear-resolved track runs with `G = 1`**, so its per-gear CPUE machinery is inert and
  gear-type catch is PE-apportioned. Do not raise `G > 1` without adding per-gear effort shares;
  only gear 1 is observed in the effort stream.
- **Repository hygiene is a decision, not a task, and it is deliberately left open.**
  `.Rproj.user/` is untracked in this branch. `05_output/` is 9,942 tracked files and 835 MB, of
  which the PNG plots and rendered HTML are **218 MB, 64% of the whole history, with no reader in
  the codebase**. Three options with their costs are in Tier 4 of
  `07_documentation/development_notes/PIPELINE_STATUS.md`. `origin/main` already carries 82 MB of
  the same blobs, so the aggressive option would have to be applied to `main` as well. Nothing is
  near a GitHub limit, so this is a clone-size question, not a blocked-push one.

## 7. Reproducing the authoritative run

```
# one command, from the repository root
source("run_estimation.R")
```

`run_config.R` ships the canonical configuration: `model <- "pooled"`, the single 2024-25 season,
and every default the authoritative run used. It is organized in five sections, method-affecting
levers first, with rarely-used and diagnostic-only toggles moved to the end. Paste-ready
commented blocks for the 2025-26 season and the blocked 2023-25 span sit in section 1.2.

Environment: R 4.2.2, rstan 2.32.7 / StanHeaders 2.32.10, pinned in the top-level `renv.lock`
(~99 CRAN packages from the confirmation run's `session_info.txt`). Expect roughly 12 to 14 h for
the five-rung ladder; a single production run is a fraction of that.

## 8. Where to read what, after the merge

| document | job |
|---|---|
| `07_documentation/BSS-GH-pooled-CPUE-model-documentation.md` | **Method v2.0**: what the model is. Specification and the eleven limitations, ordered by risk |
| `07_documentation/BSS-GH-gear-type-CPUE-model-documentation.md` | framework v6.0, the cross-check track; defers everything shared |
| `07_documentation/development_notes/PIPELINE_STATUS.md` | current state and backlog. **The authoritative run and its total are in the box at the top, and that box is the only place to take a number from** |
| `07_documentation/development_notes/CHANGE_REGISTER.md` | every change, its status, its evidence, its effect on the number; every defect and what it cost |
| `07_documentation/development_notes/VALIDATION_CAMPAIGN.md` | the dated run-by-run narrative, Sections 1b to 1v: why each decision was taken and what the run said. Also holds the git-anchor table that reconciles every section date against the commit that landed it |
| `07_documentation/archive/` | Method v1.0 for both tracks, frozen as a historical artefact, and the removed weather module |
| `07_documentation/CLAUDE.md` | the working conventions, including how to date a marker and which of the five documents to write to |

The earlier version of this file described only the first wave of the branch (the OSP integration
alone, against a 67,312 total). It is superseded in full; its record survives in git history and,
in narrative form, in `VALIDATION_CAMPAIGN.md`.

---

**Recommended merge:** no squash. The 85 commits are the audit trail; each one passes the harness
on its own, which is what makes the branch bisectable, and several commit messages are the only
record of a defect and why the fix is shaped the way it is.
