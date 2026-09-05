# Change register: the OSP boat-count branch

**Last updated:** 2026-09-09
**Branch:** `OSP-boat-count-incorporation`. **`main` is the pre-FW-creel-meeting, pre-OSP state.**
**Authoritative run:** `05_output/20260904/pooled-CPUE-AD-A1-adopted`, port total **72,027 [53,018, 101,364]**, 4 of 4 components fitted.

**WDFW has published no estimate from this pipeline.** Every "adopted" below means adopted into the working model on this branch, not released. There is no published figure a change has to stay consistent with.

Status vocabulary: **ADOPTED** (shipping in `run_config.R`, validated by a run) / **BUILT, INERT** (coded and tested, off by default, waiting on data or a decision) / **OPEN** (not resolved) / **REJECTED** (tried, not kept) / **BLOCKED** (waiting on something outside the code).

---

## A. Model and estimator changes

| # | Change | Status | Evidence | Effect on the number |
|---|---|---|---|---|
| A1 | **Gear-deployments as the effort unit** (`h = number_of_gear`, `L = tau`), shore and boat, both tracks | **ADOPTED** (v7.7) | `cpue_linearity_*.csv`: `beta_h` covers 1 on both units; catch is sub-linear in soak time (saturation 0.13-0.27), so crabber-hours fails the linearity test | Structural; the basis of every number since |
| A2 | **Shared turnover** `shared_tau`, boat-only via `shared_tau_min_obs = 15` | **ADOPTED** 2026-09-01 | Stage S2 bit-identity gate; `tau_bar` agrees across the two independent tracks to **0.05%** (2.6128 vs 2.6114) as of 2026-09-08 | Replaced 289 per-day draws with one estimated turnover |
| A3 | **PE / BSS incomplete-trip arm alignment** (`pe_gear_ratio_arm = "match_bss"`) | **ADOPTED** 2026-09-02 | Stage P1: largest change **-0.059%** on the boat PE, shore unchanged (the control) | A consistency fix, not an accuracy one |
| A4 | **Shore all-gear AR: daily -> weekly** (`ar_max_resolution$pooled$shore$all_gear`) | **ADOPTED** 2026-09-07 | Ladder (C2) fitted all four rungs; adequacy went `p_loo` 35.2% -> 9.5% of `n_obs`, Pareto k>0.7 **41 -> 0**, `cov50` dev 0.204 -> 0.066, miscalibration flag cleared. Routing proven identical to the candidate across **10,253 parameter rows** | Shore all-gear 20,898 -> 21,489; port **+0.72%** |
| A5 | **Zero-inflated shore catch likelihood** (`estimate_catch_zi = TRUE`, scoped to shore) | **ADOPTED** 2026-09-07 | At the same resolution: zero bin z +3.7 -> +2.0, one bin z -6.1 -> -3.3, pot-closure replicate closes both; elpd **+11.6 nats at 2.30 paired SE** | Shore all-gear -58 crab at weekly (-0.27%) |
| A6 | **Gear-track boat all-gear sampler settings** moved from an experiment override into the driver | **ADOPTED** 2026-09-02 | Divergences **554 -> 2**, `tau_bar` n_eff 49 -> 19,292, gear port 51,385 -> 70,953 | Restored a cross-check that had been unavailable |
| A7 | **OSP boat-count stream** (`osp_scale_is_tau = TRUE`, `kappa_OSP` decoupled) | **ADOPTED** | 130 in-window days ingested; the OSP mean uses `L[day]` rather than a free scale | The reason this branch exists |
| A8 | **OSP crabbing-only lower bound on `f`** (`use_osp_crab_lower`, `f_lower`) | **BUILT, INERT** | The column `WestportCrabOnlyEffort` is absent from the workbook, so the bound never binds and `f` sits on the set value 0.30 | None until the data arrives |
| A9 | **Other-fishery opener effort covariates** (`K_open` / `X_open` / `B_open`) | **BUILT, INERT** | `opener_covariate_mode = "off"`; the razor-dig predecessor `B3` bought no predictive gain (elpd within 1 SE) and moved the port +0.6% | None |
| A10 | **Same-day-effort density term** `gamma_C` (`estimate_cpue_density`) | **REJECTED** | Run 4 was pathological (one chain ~17x slower, treedepth saturation); `gamma_C` indistinguishable from zero where it did fit | None |
| A11 | **GR-7 Phase 2 Dirichlet gear shares** (`gear_share_dirichlet`) | **BUILT, INERT** | Coded 2026-07-21, parses under stanc 2.32.5, byte-identical to Phase 1 when off; sampled in stages V2/V3 | None at `G = 1` |
| A12 | **Weather / tide covariate module** | **REJECTED and STALE** | Excluded on its own evidence; the fork is now missing ~40 data variables the production model declares | None; do not cite its boat number |

| A13 | **Season portability** (2026-09-09): non-intersecting closure windows yield a single all-gear sub-season instead of an error; `season_filter` accepts a vector for multi-season spans; `validate_season_window()` stops loudly on a stale season/window pairing; SEASON-DERIVED tags + NEW SEASON CHECKLIST in `run_config.R`; `ar_force` reframed as the per-fit resolution pin in the new-season workflow | **ADOPTED** 2026-09-09 | Functional tests in the harness (window cases, vector filter, validator stop); the program goal recorded in `CLAUDE.md` and Section 1n | None on 2024-25 outputs (full-season configs unaffected); part-season windows go from ERROR to running |

## B. Diagnostics and infrastructure

| # | Change | Status | Why it exists |
|---|---|---|---|
| B1 | **Scale-aware convergence gate** (`bss_convergence_gate.R`) | **ADOPTED** | Decides PE vs BSS per fit on R-hat, n_eff, divergence fraction and SD-normalized divergence impact. It answers "did this fit sample", never "is this model right" |
| B2 | **Model adequacy reported beside the gate** (`bss_model_adequacy.R`) | **ADOPTED** | Six configurations once passed every gate criterion while spanning 44% on the boat. `p_loo_frac`, Pareto k count, `cov50_worst_dev`, `pit_sd_worst_dev`, `flag_miscalibrated` |
| B3 | **Randomized-PIT coverage** in both the aggregate and per-observation PPC | **ADOPTED** 2026-09-01/02 | The quantile-interval form over-covers small counts by construction and produced a phantom "trailer over-coverage" finding |
| B4 | **Decoupled-parameter flag** on `structural_params_*.csv` | **ADOPTED** | Several parameters carry a proper prior unconditionally, so an unused one reports its PRIOR in the same columns as an estimate |
| B5 | **AR escalation ladder** (`ar_escalate`, scoped per population and sub-season) | **ADOPTED as a toggle**, off by default | The FW creel team's idea. Every rung is a multi-hour fit, so it ships off; it settled A4 |
| B6 | **Per-rung adequacy** (`bss_rung_adequacy.R`, `ar_rung_adequacy = TRUE`) | **ADOPTED** 2026-09-06 | The first ladder fitted two rungs for 174 minutes and recorded nothing about either, because only the kept fit reaches `write_bss_diagnostics()` |
| B7 | **Paired-difference elpd** (`loo_elpd_paired.R`) | **ADOPTED** 2026-09-05 | A model comparison needs the SE of the paired difference, not of either total. On the ZINB that was 5.5 nats against a naive 46, turning an apparent non-effect into a 2.69 SE effect |
| B8 | **Mixture-aware PPC** (`zinb_ppc.R`) | **ADOPTED** 2026-09-05 | The R-side PPC scored a ZINB fit as plain NB2 and reported the opposite of the truth on the zero bin |
| B9 | **Posterior draw persistence** (`save_bss_ppc_draws.R`, `save_ppc_draws = TRUE`) | **ADOPTED, HALF-BUILT** | Three diagnostic defects in five weeks each forced a multi-hour re-fit. The draws are written; **nothing reads them back yet**, so the feature is still write-only |
| B10 | **Untruncated config dump** (`str(params, list.len = 1000, ...)`) | **ADOPTED** 2026-09-04 | `str()` capped at 99 keys while `run_config` carried 120+, so every `run_parameters.txt` before that date is missing ~21 keys and every config comparison built on one was unreliable |
| B11 | **Shared batch helpers** (`batch_verdict_helpers.R`) | **ADOPTED** 2026-09-04 | Four runners had four copies. `merge_csv_by` exists because an ungated overwrite once destroyed a batch's verdicts; `config_delta` because an exactness comparison was got wrong by hand twice in three weeks |
| B12 | **Regression harness** (`test_improvements_2026-08-25.R`) | **ADOPTED**, **452 assertions** | Runs in seconds without rstan. Sections 35-44 read driver SOURCE and assert the property that failed, because most defects here have been in code that validates, not code that models |

## C. Defects found and fixed, with what each cost

| Date | Defect | Cost | Fixed by |
|---|---|---|---|
| 2026-08-25 | rstan reports a data-init failure on the message stream and returns an empty `stanfit` rather than raising | six model-runs | `bss_stan_fit.R` data-contract guard |
| 2026-08-28 | `ar_force` keyed per population silently forced a second sub-season | 3,025 crab on a component | per-sub-season resolution and a scope control in every batch |
| 2026-08-31 | PPC coverage from a quantile interval over-covers small counts | a phantom open item, retracted | B3 |
| 2026-09-01 | `tau_bar` registered under its bare name while rstan writes `tau_bar[1]`; the lookup threw inside a `tryCatch` | whole `prior_vs_posterior` files silently dropped | row-key helper |
| 2026-09-02 | The gear boat all-gear fit had never sampled in production | a 27% understatement of the cross-check | A6 |
| 2026-09-03 | `force_res` gated on `isTRUE(params$ar_escalate)`, which is FALSE for the scoped form | **14.7 h**: four ladder rungs all refit `daily` | resolve through the one helper; harness asserts no driver contains the expression |
| 2026-09-03 | The ZINB PPC computed under NB2 | the feature was nearly rejected on a diagnostic that reported the opposite of the truth | B8 |
| 2026-09-03 | The elpd criterion used `se_elpd_loo` | "NOT WORTH IT" on a 2.69 SE effect | B7 |
| 2026-09-04 | `run_parameters.txt` truncated at 99 keys | one false caveat, one latent missing flag | B10 |
| 2026-09-06 | The runner never moved the rendered HTML | **stage L1's report destroyed** by the next stage; numbers survived in CSV | move it, gitignore `01_BSS_models/*.html`, assert over every runner |
| 2026-09-06 | `bss_rung_adequacy` reseeded inside the fitting loop | a diagnostic-only toggle moved a committed diagnostic on a **bit-identical** fit | save and restore `.Random.seed` |
| 2026-09-06 | `p_loo_frac` summed streams where `model_adequacy` takes the worst | two statistics under one name in one run | match the existing definition |
| 2026-09-08 | The two drivers label the port-total row differently; `is.finite(numeric(0)) && ...` errors | the adoption run's verdicts lost after 4.1 h of fitting | pattern match; every verdict block wrapped |
| recurring | `DRY_RUN <- FALSE` committed | a fresh clone starts real fits | asserted since 2026-09-03; caught every time, **five occurrences** |

## D. Open items, ordered by how much they could move the number

| # | Item | Status | What would settle it |
|---|---|---|---|
| D1 | **The boat is 43% of the port and 79.8% extrapolated** (66 sampled days of 289) | **OPEN, not a code problem** | More boat interviews in 2025-26. No modelling change manufactures uncollected information |
| D2 | **OSP crabbing-only daily counts** | **BLOCKED on a data request** | One column in `WBL_boat_counts.xlsx`. Machinery built and inert (A8). Note it is a LOWER BOUND on crabbing vessels, not `f` |
| D3 | **The gear track's shore cap is still monthly** | **OPEN** | Now the largest single contributor to the -1.39% cross-track gap. A gear-track ladder, about 3 h. At the same resolution the tracks agree to **0.08%** |
| D4 | **The catch stream is under-covered at every AR resolution** (-3.5 to -4.6 SD) | **OPEN, bounded** | A hurdle or two-component NB mixture. The ZINB halved both count bins and did not close the shore all-gear one bin (still 3.3 SD). Gain bounded at ~6% of the catch |
| D5 | **`ppc_draws_*.rds` are write-only** | **HALF-BUILT** | A recompute script. Code, no run |
| D6 | **The ZINB is pooled-only** | **OPEN, low priority** | The two tracks now differ in the shore catch likelihood (~-0.3%). A Stan edit plus a recompile restores symmetry |
| D8 | **One closure window per run** | **OPEN, by design for now** | A multi-season span containing two pot closures cannot be expressed; run per season and combine. Generalization = a closures table feeding `build_subseasons()` N ordered windows |
| D9 | **New-season workflow documented** | **DONE 2026-09-09** | `07_documentation/NEW_SEASON_GUIDE.md`: naive run -> ladder -> per-rung adequacy -> pin (`ar_force`) -> cap -> production -> cross-check |
| D7 | **Method of record is frozen at v1.0** (pooled code v7.4); the code is past v7.9 | **OPEN, documentation debt** | The two method documents carry pre-refresh 2024-25 reference numbers that have never been regenerated |
