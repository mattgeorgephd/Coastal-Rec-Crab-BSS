# Improvement plan, 2026-08-27

- **Written after:** the 2026-08-26 validation ladder (six model-runs, ~15 h) and the review of its results.
- **Branch:** `OSP-boat-count-incorporation`, ladder results at `ec9a52d`.
- **Companion:** `PIPELINE_STATUS.md` (state and backlog), `ladder-validation-review-2026-08-27.md` (the review this plan comes from).
- **Convention:** no em dashes.
- **STATUS 2026-08-29:** stages A-F have RUN (`06_diagnostics/run_improvement_plan_2026-08-27.R`, ~22 h). Results and what they changed are marked per item below; the new work they generated is in the new Stage 5. Read `improvement-batch-review-2026-08-29.md` alongside this.

This is a sequenced plan, not a wish list. Every item states what it changes, what "done" looks like as a check somebody else can run, and roughly what it costs. Items in Stage 0 are already applied on this branch; the rest are ordered so that each one's result is interpretable before the next begins.

The organising judgement: **the 2026-08-25 batch is validated and should be adopted, and the largest open question in the project is not in that batch at all.** It is `osp_scale_is_tau`, which predates it, roughly doubles the private-boat harvest, and does not do what the documentation says it does.

---

## Stage 0. Applied in this patch (2026-08-27)

Six defects the ladder exposed. None of them changes a fitted number; all of them change what a future run can be audited on, which is why they go first.

| # | Change | Files | Why now |
|---|---|---|---|
| 0.1 | The **second copy** of the shore day-length weighting is fixed. `.srd_monthly_share()` now weights by day length only when the fit's `L` is a day length in hours, and count-weights when `L` is a turnover. | `save_run_diagnostics.R` | Improvement 1's defect (B) was fixed in `pe_monthly_effort_share.R` and missed here, so `monthly_pe_vs_bss.csv` kept the pre-v7.7 split. The ladder proved it: that file's shore PE effort was numerically unchanged from the pre-fix reference. |
| 0.2 | The **empty-effort-stratum audit is a file**, `pe_empty_effort_strata.csv`, written by both drivers. | `run_pe_pooled.R`, `run_pe_gear.R`, `save_run_diagnostics.R`, both drivers | The counts were only `cat()`-ed and the pooled PE chunk is `results='hide'`, so on that track the report reached nobody. See the finding in Stage 1.4, which this fix produced within minutes of existing. |
| 0.3 | **Unit provenance** in `fit_data_summary.csv`: `ie_obs_col`, `ie_obs_unit`, `effort_unit`, `L_unit`, `n_interviews_fitted`. `ie_analysis.csv` and `L_effective_ie_detail.csv` gain `ie_trips`, `turnover` and `hours_per_trip`. | `bss_day_length.R`, both preps, `save_run_diagnostics.R`, pooled driver | Rung 2 changed which column the shore I/E likelihood consumes and no run output recorded it. Auditing that change meant recomputing `ie_trips` from the raw workbook. |
| 0.4 | **Decoupled parameters are labelled.** `structural_params_*.csv` gains `decoupled` and `decoupled_reason`, from a new `bss_decoupled_reasons()`. | `model_diagnostics.R` | Under the production `osp_scale_is_tau = TRUE`, `kappa_OSP` reports its prior exactly (median 3.008, 95% 1.63-5.40) in the same columns as estimated parameters. Read as an estimate, that is a claim the model measured the boat turnover at 3.0. It did not; it was told 3.0. |
| 0.5 | **`production_arm` is per estimator**, via `incomplete_trip_production_arm()`. | `diagnose_incomplete_trips.R` | The single "exclude" label was right for the BSS and wrong for the boat PE. The shipped boat PE equals this table's `gear_only` arm, which is improvement 5's central finding, and one flat label hid it. |
| 0.6 | **Skipped CPUE diagnostics say so.** A saturation or linearity regression below its 30-interview minimum now writes a `status = "skipped"` row instead of no file. | `bss_cpue_diagnostics.R` | With the floor at 15, the boat pot closure fits on 17 interviews and two of its three CPUE diagnostics are legitimately unavailable. Two missing files among a hundred read as an error. |

Also in Stage 0, in the validation runner:

- **The gear exactness criterion is repaired.** It demanded `delta == 0` on the port total and FAILED at delta = -19 while the model was bit-identical the whole time. The port total is assembled from `rstan::extract(permuted = TRUE)` draws and an unseeded `sample.int()`, so it moves by a fraction of one MCSE even when the fits do not; `PIPELINE_STATUS` had already recorded that. Exactness is now tested where it exists: byte equality of `convergence_report.csv` plus full-precision equality of every shared parameter row across `bss_summary_*` and `bss_full_summary_*`. Re-judging the 2026-08-26 outputs, it now reports **9,895 shared parameter rows identical at full precision** and passes.
- **The tau-narrowing criterion is retired to INFO.** See Stage 2.1: the prediction was structurally impossible, not merely unmet.
- **A per-rung convergence criterion is added.** Rung 2 lost its shore all-gear fit to the gate and nothing noticed; the summary reported the rejected fit's numbers and the port total silently carried the PE instead. `method_selected` per component is now carried into the summary, and a rung that costs a component its fit FAILs.
- **`empty_effort_note` is populated** from 0.2 instead of being a dead column.

Re-judging the existing ladder outputs with the repaired criteria gives **13 PASS, 1 FAIL, 3 INFO**, where the single FAIL is the rung-2 convergence failure that the original set missed entirely.

**Acceptance:** `Rscript 06_diagnostics/test_improvements_2026-08-25.R` reports 104 assertions passing; sourcing the runner with `DRY_RUN <- TRUE` reports pre-flight clean and self-test 8/8.

---

## Stage 1. Close out the 2026-08-25 batch

### 1.1 Adopt the rung-4 configuration
**Do:** nothing to the code. `run_config.R` already ships it.
**Numbers:** port total 66,237 (95% CI 50,037 - 91,210), effort 40,187, 4 of 4 components fitted.
**Basis:** the baseline reproduces on both tracks, each change moved for its predicted reason, and the two independent pipelines reconcile to 0.73%.
**Done when:** `PIPELINE_STATUS.md` names rung 4 as the authoritative run (done in this patch).
**Cost:** none.

### 1.2 Re-run rung 2 with a different `bss_seed`  <!-- DONE 2026-08-29 -->
**RESULT (stage B): answered, with a sharper conclusion than the question expected.** The fit survived on the new seed (488 divergences, 4.88%, worst chain 350) so the 2026-08-26 failure was seed-dependent. But 4.88% against a 5% backstop and `E_sum` n_eff of 924 against the shipped config's 5,914 make this the worst fit that passed anywhere in the project. Two seeds, two struggling chains. The I/E unit fix alone DOES strain the shore all-gear daily fit; the weekend fix is what makes the shipped configuration comfortable. The clean isolation is now available and confirms the I/E result three ways: `sigma_IE` 0.238 (stage B) / 0.300 (rung 3) / 0.272 (rung 4) against 1.068 before.

**Why:** rung 2 is the only rung without a clean isolation. Its shore all-gear fit was rejected by the gate at 2,957 divergences, but the per-chain record says chains 1, 2 and 4 recorded 163 / 140 / 154 divergences, in line with rung 1's 385 across four, while **chain 3 recorded 2,500 of 2,500** with mean accept-stat 0.053 and mean treedepth 1.21. It never moved. One stalled chain and a structural pathology look identical in an aggregate divergence count and are completely different findings.

**Do:**
1. In `run_patch_validation_2026-08-25.R` set `RUNGS <- 2` and `MODE <- "isolated"`.
2. Add `bss_seed = <a different integer>` to the rung-2 delta, or set it in `run_config.R` for the one run.
3. Run. It is one pooled model-run.

**Done when:** either (a) 4 of 4 chains behave and the fit passes the gate, in which case the I/E unit fix is exonerated and rung 2 becomes the clean isolation the ladder lacked; or (b) it stalls again, in which case the fix genuinely destabilises the daily-AR shore fit and Stage 3.2 becomes urgent rather than optional.
**Cost:** ~3 h.
**Risk if skipped:** the headline `sigma_IE` result rests on rung 3, which also carries the weekend change. The two are almost certainly separable (rung 2's 0.345 and rung 3's 0.300 sit close together, far from the reference's 1.057), but "almost certainly" is not the standard this ladder was built to meet.

### 1.3 Correct the documentation
Applied in this patch: the tau-informativeness claim retracted, the extrapolation figures regenerated, the hours-to-trips mismatch recorded as 3.2x rather than 4x, and the ladder result written into both development histories.

### 1.4 Decide what to do about the empty-effort strata  <!-- RUN 2026-08-29, DECISION OPEN -->
**RESULT (stage A): `day_type` raises every affected component's PE by 20-27%** (boat all-gear +21.2%, boat pot closure +26.8%, shore all-gear +20.7%, shore pot closure unchanged), roughly +13% at the port level. The `zero` fill has been biasing the PE down by about an eighth. Caveat: stage A did not load the holiday workbook, so its absolute stratum counts (41/44/8/0) are not the driver's (44/47/9/0); the comparison is sound because both arms shared the stratification. Fixed in the 2026-08-29 patch. **The decision is still open, and it is now a real one** - see Stage 5.1.

**New finding, from fix 0.2 within minutes of it existing.** The zeroed-day counts on the shipped pooled config are much larger than the one number that had been quoted:

| Component | Empty strata | Zeroed days | Share |
|---|---:|---:|---:|
| shore pot closure | 0 / 23 | 0 / 76 | 0.0% |
| shore all-gear | 20 / 92 | 44 / 289 | **15.2%** |
| private boat pot closure | 5 / 23 | 9 / 76 | **11.8%** |
| private boat all-gear | 22 / 92 | 47 / 289 | **16.3%** |

Three of four components sit above the 5% warning threshold at the default `zero` fill. The figure previously on the record, 4 of 76 rising to 9 of 76, was the boat pot closure alone; the two all-gear components were never reported on the pooled track at all, because of the defect fixed in 0.2.

**What it does and does not affect.** These are PE numbers. Where the BSS passes the gate the PE is a cross-check, not the reported number, so the headline port total is not biased by this. It matters in three places: the boat pot closure entered the port total as a PE point in rungs 1-3; the PE-vs-BSS comparison is the project's main sanity check and is being run against a PE that is biased down; and any component that falls back to PE in a future run inherits the bias directly.

**Do:** run the shipped config once with `pe_empty_effort_stratum = "day_type"` and compare `pe_port_summary.csv` and `pe_vs_bss_comparison.csv` against rung 4. This is a PE-only change, so it does not need a full refit to interpret, but it does need a run.
**Done when:** the PE-vs-BSS ratios are recorded under both fills and a default is chosen deliberately rather than inherited.
**Cost:** ~4 h for a pooled run.
**Note:** the unsampled strata skew weekend and holiday, and those days carry 1.7 to 2.3x weekday effort, so the bias is worse than the day count suggests. The `day_type` fill exists precisely for this and has never been exercised in production.

---

## Stage 2. Model corrections

### 2.1 A shared `tau` with per-day deviations  <!-- RUN 2026-08-29: BOAT CONFIRMED, SHORE MUST BE SUPPRESSED -->
**RESULT (stages D + E). The gate PASSED** (11,217 shared parameter rows identical at full double precision), so the edit is behaviour-neutral when off. **With it on, the boat turnover moves off its prior and four independent lines agree**: `tau_bar` = 2.597 [2.064, 3.249] against the overlap calibration's 2.01-3.03; the trailer and OSP PIT means converge on 0.50 simultaneously and in opposite directions (bias down 81% and 73%); boat trailer `elpd_loo` +14.67 at paired z = +2.66 with `p_loo` staying low; and convergence improved (divergences 97 -> 39). Boat all-gear catch 25,883 -> 31,012, port 66,237 -> 75,619. Two qualifications: the posterior sits above the prior's 95% envelope, so 2.597 is shrunk rather than prior-driven; and `n_ie_obs = 0` for this fit, so the turnover is identified indirectly through the count streams and must never be called measured. **The SHORE result must not be carried forward** - it moved +17.9% on 4 of 289 observations, its interval still contains the 1.7 prior centre, it improved nothing measurable, and the gear track put the same parameter at 1.681. See Stage 5.2. The gear track's boat fit failed on a single stuck chain (`phi_E` at the unit-root boundary, 544 of 554 divergences in chain 3, 1,000 iterations against pooled's 2,500) and gives no verdict; its boat POT-CLOSURE fit passed at `tau_bar` 2.050 [1.415, 2.930], weaker but real support.

**Status: the model change is IMPLEMENTED and verified; the runs are stages D and E of `06_diagnostics/run_improvement_plan_2026-08-27.R`.** `shared_tau` ships FALSE. `tau_bar` is declared `vector<lower=0>[shared_tau]`, so it is zero-size when the feature is off and the unconstrained parameter vector is unchanged; a fixed-seed fit on real data reproduces the pre-change model at full double precision (verified 2026-08-27 on the shore pot-closure data: `C_sum`, `E_sum`, `B1`, `B2`, `B1_C`, `R_G`, `sigma_IE` and `L_out` all identical to the last digit, zero `tau_bar` draws, `tau_bar_out` exactly 0). With it on the model samples cleanly and the per-day `L` collapses onto one level as designed. `bss_shared_tau_data()` refuses the toggle, with a warning rather than an error, whenever `L` is not a constant turnover.

**The defect.** `crab_bss_pooled.stan` declares

```stan
vector[D * estimate_L] L_raw;
...
L[d] = L_data[d] * exp(L_prior_sigma[d] * L_raw[d]);   L_raw ~ std_normal();
```

There is no shared `tau_shore` or `tau_boat` parameter. `L` is `D` independent per-day draws, each anchored on the same prior centre with log-SD 0.3. Nothing pools information across days, so no amount of observation on a subset of days can move the season-level turnover.

**What that costs, concretely.**

- Shore: four in-window I/E days can inform four days out of 289. The season median `L` posterior is 1.6998 against a prior centre of 1.7000. The claim that the fixed I/E stream would "become genuinely informative about `tau_shore`" cannot be true as the model is written.
- Boat: **this is the expensive one.** With `osp_scale_is_tau = TRUE` the OSP likelihood uses `L[day]` as the turnover, and 148 days of OSP data move it almost not at all: posterior median 1.201 against a prior centre of 1.200, with only 26 of 289 days narrowing below 90% of the prior width. Meanwhile `osp_trailer_overlap_calibration.csv` puts the implied turnover at 2.01 to 3.03 and the Phase-1 free `kappa_OSP` posterior sat at 3.15. That roughly 2.5x conflict is absorbed by the OSP overdispersion instead: `r_OSP` lands at 1.58 to 1.62, a negative binomial loose enough to treat a systematic scale error as noise.

**Independent corroboration.** The boat trailer PIT mean is 0.420 to 0.425 in all six ladder runs, against a nominal 0.50, while every shore and catch stream sits at 0.498 to 0.516. The model predicts more trailers than were observed. That is the signature of `lambda_E` being pulled up by OSP against a trailer stream that disagrees with it, which is exactly what a scale conflict with nowhere to go looks like.

**Why it matters.** `osp_scale_is_tau` roughly doubles the private-boat harvest, and the boat is about 40% of the port total. As implemented, the size of that doubling is set by the `tau_boat` prior rather than by the OSP data. That is not the same claim as "the dense OSP series identifies the boat turnover", which is what the notes say.

**Do:**
1. Add `real<lower=0> tau_bar;` with the existing lognormal prior, and re-express `L[d] = tau_bar * exp(sigma_tau * L_raw[d])` with `sigma_tau` a small estimated or fixed day-to-day scale. Both models.
2. Keep `estimate_L = 0` reproducing the current fixed-`L` behaviour exactly, and keep the per-day deviation so no flexibility is lost.
3. Guard it behind a config toggle, off by default, exactly as every other behaviour-changing item in the 2026-08-25 batch was.

**Done when:**
- With the toggle off, a fixed-seed run is bit-identical to rung 4 on `bss_summary_*` and `convergence_report.csv` (the exactness test from Stage 0 now does this for you).
- With it on, the boat `tau_bar` posterior either moves toward the 2.0 to 3.0 the overlap calibration implies, or stays at 1.2 with the interval to justify it. Either answer settles the question; the current model can produce neither.
- The trailer PIT mean moves toward 0.50, or does not, and we learn which.

**Cost:** a day of Stan work, then two pooled runs plus one gear run, so roughly 12 to 16 h of compute.
**Risk:** this changes the boat number. It should be run and reviewed before any external presentation of the boat figure, not after.

### 2.2 A posterior-predictive check on the OSP stream  <!-- IMPLEMENTED 2026-08-27 -->
**Status: IMPLEMENTED.** `ppc_calibration_*.csv` now carries an `osp` row for every boat fit, using the same mean the Stan likelihood uses, and the batch reports its PIT mean beside the trailer's.

`model_diagnostics.R` builds PPC calibration for the gear, trailer and catch streams and **not** for OSP, which is why the scale conflict in 2.1 had to be inferred from the trailer PIT rather than read directly.
**Do:** add an `osp` part to the `calib()` block, using `(lambda_E / R_G_boat) * (osp_scale_is_tau ? L : kappa_OSP)` as the mean and `r_OSP` as the dispersion.
**Done when:** `ppc_calibration_*.csv` carries an `osp` row for every boat fit, and its PIT mean is reported alongside the trailer's.
**Cost:** an hour, plus whatever run it rides along with.

---

## Stage 3. Sampler and resolution

### 3.1 Give the boat pot closure one AR resolution  <!-- ANSWERED 2026-08-29 -->
**RESULT (stage C): the 13% cross-track gap was a configuration artefact.** Forced to biweekly, the pooled boat pot closure gives 735 against the gear track's 743, 1.1% apart. Align the two `ar_max_resolution` maps. The run was contaminated by an `ar_force` defect this review found (it ignored the sub-season key and forced both boat sub-seasons, moving boat all-gear 25,883 -> 28,902), so stage C's port total is not interpretable; the pot-closure conclusion stands. Defect fixed in the 2026-08-29 patch.

The two tracks fit the same component at different resolutions and get different answers: pooled at monthly (`P_n = 3`) gives 850, gear at biweekly (`P_n = 6`) gives 743, a 13% gap, and their over-dispersion decompositions differ substantially (`r_E` 1.13 against 2.39). Both pass the gate.
**Do:** reconcile the two `ar_max_resolution` maps for `private_boat` / `pot_closure`, or let Stage 3.2 choose.
**Done when:** the two tracks agree on the resolution, and the residual gap is reported as a genuine cross-track difference rather than a configuration artefact.
**Cost:** two runs.

### 3.2 Exercise `ar_escalate` on the shore all-gear component  <!-- RUN 2026-08-29: ANSWER IS UNCOMFORTABLE -->
**RESULT (stage F): every component passes at DAILY on the first attempt, and the daily boat fit is absorbing rather than identifying.** The ladder never escalated. Shore pot closure at daily gives 111 divergences where Run 6 recorded 1,165 - the funnel that justified the biweekly cap is gone, so `ar_max_resolution` is now the only thing holding production at monthly/biweekly and it is doing real work on the number. But the boat at daily raises catch 44% (25,883 -> 37,392) while `p_loo` on its trailer stream goes 4.9 -> 31.8 with 2 bad Pareto k, the CATCH stream's elpd gets worse (-451.4 -> -455.3), `r_OSP` collapses 1.62 -> 10.61 with a 95% upper bound of 1,405, and its parent `sigma_r_OSP` drops to n_eff 307 (below the gate's own floor, which the gate does not check on that parameter). The entire increase lands in the extrapolated portion. See Stage 5.3.

Improvement 7 shipped `FALSE` and has never been run. The fit it was designed for is the pooled shore all-gear at daily AR, which is the most strained fit in the project:

- `p_loo` = 109 effective parameters on 311 gear observations, 41 Pareto-k above 0.7 (max 1.16). The gear track's monthly-AR fit of the same data: `p_loo` 11.3, zero high-k.
- 143 of the 289 daily AR periods carry no observation at all.
- `elpd_loo` is better at daily (-1,251 against -1,353), but at k > 1 that margin is not a trustworthy model comparison.
- It is the fit that stalled a chain in rung 2, and historically the one that wobbles across the gate.

**Do:** run the shipped config with `ar_escalate = TRUE` and `ar_escalate_respect_cap = FALSE`, and read `ar_escalation_log.csv`.
**Done when:** we know whether the shore all-gear component settles at daily, weekly or biweekly on its own sampler behaviour, and what that costs in `elpd_loo` and in the reported catch.
**Cost:** the escalation ladder refits each rung, so expect 2 to 3x a normal run, roughly 8 to 16 h. Order it after Stage 1.2, since a reseeded rung 2 may change what "strained" means here.

---

## Stage 4. Reporting and hygiene

### 4.1 Seed the diagnostic subsampling  <!-- IMPLEMENTED 2026-08-27 -->
**Status: IMPLEMENTED.** `write_fit_extended_diagnostics()` takes a `seed` argument, wired to `params$bss_seed` from both drivers. The port total remains RNG-sensitive through `rstan::extract(permuted = TRUE)` regardless, which is now stated in the header.

`save_run_diagnostics.R` draws `sample.int()` twice without a seed, for the saved draw subset and the PIT subset, so `bss_draws_summed_*`, `ppc_*` and everything built on them shuffle between two runs of identical fits. That is what turned the Stage 0 exactness criterion into a false alarm.
**Do:** seed from `params$bss_seed` at the top of `write_fit_extended_diagnostics()`, and note in the header that the port total remains RNG-sensitive through `rstan::extract(permuted = TRUE)` regardless.
**Cost:** minutes. Do it with the next code change, not on its own.

### 4.2 Make the commercial/charter census's day-typing explicit  <!-- IMPLEMENTED 2026-08-27 -->
**Status: IMPLEMENTED.** The method label is now "Census (day-type stratified)" in both drivers.

It moves 11,986 to 11,821 under the weekend redefinition, because `estimate_comm_charter.R` stratifies on `params$days_wkend`. That is defensible, but a column labelled "Census" in `pe_vs_bss_comparison.csv` that responds to a day-typing toggle will mislead somebody eventually.
**Do:** relabel it "Census (day-type stratified)" and note the dependency in the method document.

### 4.3 Reconcile the empty-stratum prediction with its sign
The recorded rationale for rung 3 said "expect a downward move" in the boat pot-closure PE. The count prediction was exactly right (4 of 76 to 9 of 76) and the PE moved **up** 14%, from 351 to 400. Something in the reasoning is wrong even though the arithmetic was right. Worth ten minutes with `pe_effort_strata.csv` before the same rationale is reused.

### 4.4 Refresh the documented extrapolation figures
The notes carry 79.8 / 48.5 / 30.9%. Observed across the six runs: boat all-gear 81.8 to 82.0%, shore pot closure 31.8 to 33.3%, shore all-gear 48.5% at rung 1 rising to 51.0% shipped and 54.2% on the gear track. The boat pot closure, at **76.2%**, is undocumented and is the most extrapolated stratum in the set. Applied in this patch.

---

---

## Stage 5. What the 2026-08-29 batch generated

The batch answered its questions. In answering them it turned a suspicion into a measurement, and that measurement is now the organising problem.

**The private-boat all-gear estimate ranges over 25,883 to 37,392 - a 44% spread, 18% at the port total - across six configurations that ALL pass the convergence gate.** The 2026-08-25 batch, by comparison, moved the port total 1.6%. The gate tests whether a fit sampled, not whether the model is the right one, and it was never meant to do the latter. Everything below exists to narrow that spread on evidence rather than on default settings.

### 5.1 Adopt a fill for the empty effort strata  <!-- decision, then one run -->
**Where it stands.** `day_type` raises the PE 20-27% per affected component, about 13% at the port. Three of four components sit above the 5% warning threshold at the shipped `zero` fill. The unsampled strata skew weekend and holiday, and those days carry 1.7-2.3x weekday effort, so `zero` is biased down by construction and the only question is by how much.

**Do:** re-run stage A on the fixed script (it now loads the holiday workbook, so the counts will be the driver's). Then decide. This is a PE-only change, so it costs minutes, not a fit.
**Recommendation, stated so it can be argued with:** adopt `day_type`. `zero` is not a neutral default; it is an assumption that unsampled weekend days had no crabbers on them, which is false. The counter-argument is continuity with every published number to date, and it is a real cost.
**Done when:** the fill is chosen deliberately and `PIPELINE_STATUS.md` records which, with the PE-vs-BSS ratios under both.
**Cost:** minutes to re-run, plus the decision.

### 5.2 Make `shared_tau` boat-only, then re-run it  <!-- code, then 2 runs -->
**Why.** The boat result is well supported by four independent lines. The shore result is not supported by anything: 4 of 289 informed days, an interval containing the prior centre, no measurable improvement, and no replication in the gear track. One global toggle produced both.

**Do:**
1. Gate `bss_shared_tau_data()` on a minimum informed-day count, `shared_tau_min_obs`, defaulting to something like 20. A fit with fewer degrades to the per-day parameterization with a warning, exactly as the existing turnover guard does. The boat qualifies through its 148 OSP days; the shore's 4 I/E days do not.
2. Re-run pooled with `shared_tau = TRUE` and confirm the shore returns to 20,898 while the boat holds near 31,012.
3. Re-run the gear track at 2,500 iterations (it runs boat all-gear at 1,000, which is where its chain 3 stalled at the `phi_E` unit-root boundary) and get the cross-track verdict the batch could not produce.

**Done when:** the shore is unmoved, the boat is unchanged from stage E, and the gear track either corroborates the boat or contradicts it on a fit that passed.
**Cost:** an hour of code, then ~4 h pooled + ~0.5 h gear.
**Do not adopt `shared_tau` in production before this.** It currently moves the shore on four observations.

### 5.3 Separate AR resolution from turnover  <!-- 2 runs, the important one -->
**Why.** Stage E (shared tau at monthly) and stage F (per-day tau at daily) both raise the boat, by +20% and +44%, and they have never been varied jointly. They may be measuring the same thing two ways, or they may compound. Right now nobody can say which, and the difference between them is roughly 6,000 crab.

The evidence already separates them in quality. Stage E improves the stream it should (trailer elpd +14.67, `p_loo` 4.9 -> 8.0, both PIT means converging on nominal). Stage F buys a larger elpd gain by spending ~27 effective parameters, makes the CATCH stream worse, and collapses `r_OSP` to 10.61 with a 95% upper bound of 1,405 - a latent process eating the observation error. **The working hypothesis is that stage E is identification and stage F is absorption.** It is a hypothesis, and it is testable.

**Do:** run the 2x2. `shared_tau` in {off, on} x boat AR in {monthly, daily}. Two of the four cells already exist (stage D and stage F are the off row; stage E is on/monthly), so this is **one new run**: `shared_tau = TRUE` with the boat forced to daily.
**Read:** if shared tau at daily lands near shared tau at monthly, the shared level is doing the work and the daily AR is redundant. If they compound to something near 45,000, neither is identifying and both are absorbing.
**Done when:** the boat estimate is attributed to one mechanism, with `p_loo`, `r_OSP` and the two PIT means quoted for each cell.
**Cost:** ~4 h for the one missing cell.

### 5.4 Decide what `ar_max_resolution` is for  <!-- decision -->
**Why.** The caps were introduced because specific fits funneled at finer resolutions. Stage F shows those funnels are gone: shore pot closure at daily gives 111 divergences where Run 6 recorded 1,165. The caps are now the only thing holding production at monthly and biweekly, and they are worth 44% on the boat.

A cap justified by a pathology that no longer exists is not a cap, it is a preference. Either re-derive it from current diagnostics, or state plainly that it encodes a judgement about how much daily structure the boat series can support and defend that judgement on its own terms. **The `p_loo` evidence in 5.3 is the right basis** - a resolution that triples effective parameters and makes the catch stream worse is not the one to pick, whatever the gate says.
**Done when:** each cap entry cites either a current diagnostic or an explicit modelling judgement.
**Cost:** none beyond 5.3's output.

### 5.5 Give the gate something to say about model choice  <!-- design, then code -->
**Why.** Six configurations, all passing, spanning 44% on the boat. The gate answers "did this fit sample?" and the project has been reading it as "is this fit good?". Nothing currently stops a future run from adopting the daily boat AR on the strength of four green ticks.

**Do:** add the model-adequacy diagnostics the project already computes to the per-fit gate report, without necessarily gating on them at first: `p_loo` as a fraction of `n_obs`, the count of Pareto k > 0.7, the PIT mean per stream, and a flag when an observation-model dispersion parameter has n_eff below the gate's own floor (`sigma_r_OSP` at 307 in stage F would have tripped this). Report them beside `pass_convergence` so the two questions are visibly different.
**Done when:** `convergence_report.csv` carries both, and the distinction between "sampled" and "adequate" is documented.
**Cost:** half a day, no new runs; the quantities are already in `loo_summary_*` and `ppc_calibration_*`.

### 5.6 Smaller items the review turned up
- **Align the boat pot-closure AR maps** across the two tracks (plan 3.1 is answered: biweekly reconciles them to 1.1%).
- **`tau_bar` is missing from `prior_vs_posterior_*.csv`**, so the one parameter this batch existed to evaluate has no contraction or `prior_influential` diagnostic. Add it.
- **`expansion_ratios.csv` prints decoupled shore `R_G_boat` values without the flag** that `structural_params_*.csv` now carries. Propagate `decoupled` to that file.
- **`n_interviews_fitted` is NA for every gear-track fit** while the pooled track populates it.
- **The gear track labels boat `ie_obs_unit` as "boat trips -> tau"** where pooled says "boat trips". Cosmetic, but they should agree.
- **Re-run stage C** on the fixed `ar_force` to get an uncontaminated port total for the boat pot-closure change.

## Sequencing

Stages A-F have run. What remains, cheapest and most decisive first:

```
5.1  empty-stratum fill ....................  minutes (re-run) + a decision
   |
5.2  shared_tau boat-only, then re-run ......  1 h code + ~4.5 h
   |     <- do this before anything else touches the boat number
5.3  the 2x2: tau x AR resolution ..........  ~4 h (one new cell)
   |     <- the run that settles what the boat estimate actually rests on
5.4  what ar_max_resolution is for .........  decision, no new compute
5.5  gate reports model adequacy ...........  half a day, no runs
5.6  smaller items ..........................  alongside whatever runs next
```

Roughly 9 h of compute and a day of code. Note that 5.2 gates 5.3: running the 2x2 with a
`shared_tau` that also moves the shore on four observations would confound it again.

**One thing to hold on to while working through this.** The 2026-08-26 ladder established that the 2026-08-25 batch is behaviour-neutral where it claims to be and moves the port total 1.6%. The 2026-08-29 batch established that ordinary modelling choices move it 18%, and the boat component 44%, without tripping a single gate. The shipped configuration is not wrong; it is one point in a space the project had not measured. Until 5.3 is done, the honest way to report the boat is with the configuration named beside it. And `crab_fraction_set = 0.3` is still a placeholder and still the single largest lever on that component - a field question, not a code one.
