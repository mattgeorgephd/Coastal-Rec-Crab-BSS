# Stage 5 batch review: what the 2x2 settled, and what it exposed

> **SUPERSEDED 2026-09-01.** The authoritative run this review crowns has been superseded; the current one is in the box at the top of `PIPELINE_STATUS.md` (`20260904/pooled-CPUE-AD-A1-adopted`, 72,027, as of 2026-09-08). This document is the record of what was true when written.

- **Batch:** `06_diagnostics/run_stage5_2026-08-30.R`, six stages, branch `OSP-boat-count-incorporation`.
- **Landed:** `bd8e298` (the patch, byte-identical to what was delivered) and `d9ed8ff` (the results).
- **Runtime:** 24.0 h of fitting (S2 262 min, S3 28.5 min, S4a 327 min, S4b 535 min, S5 249 min) plus minutes for S1.
- **Verdicts as the batch scored itself:** 11 of 12 criteria PASS / ADDITIVE / INFO / DECIDE, 1 FAIL.
- **Harness:** 199 assertions, 0 failing, including guards that recompute the interaction, the calibration ordering and the cross-track `tau_bar` agreement from the committed outputs.
- **Convention:** no em dashes.

---

## 1. The headline

**The 2x2 is decisive, and it does not say what the point estimates alone would suggest.**

| private-boat all-gear | boat AR monthly | boat AR daily |
|---|---|---|
| **shared turnover OFF** | 25,868 | 37,359 |
| **shared turnover ON** | 31,008 | 42,344 |

The two levers are almost exactly **additive**: the interaction is **-155 crab** against main effects of +5,140 and +11,491, that is 3% of the smaller main effect. They are not two measurements of one thing. They act on different parts of the model, and either could be adopted without the other.

Additivity settles the arithmetic and adjudicates nothing. What adjudicates is calibration and complexity, and on those the four cells are not close:

| cell | trailer elpd | catch elpd | p_loo (trailer) | Pareto k>0.7 | OSP coverage_50 | OSP PIT sd | smallest sigma_r_* |
|---|---:|---:|---:|---:|---:|---:|---:|
| OFF x monthly | -573.8 | -451.4 | 4.9 | 0 | 0.508 | 0.277 | 0.768 |
| **ON x monthly** | **-559.1** | **-451.5** | **8.0** | **0** | **0.508** | **0.285** | **0.767** |
| OFF x daily | -540.4 | -455.3 | 31.8 | 2 | 0.931 | 0.160 | 0.307 |
| ON x daily | -473.4 | -455.1 | 52.5 | 16 | 0.977 | 0.138 | 0.086 |

*Nominal for a calibrated model: coverage_50 = 0.500, PIT sd = 1/sqrt(12) = 0.289. Trailer n = 195, OSP n = 130.*

**The shared turnover buys 14.7 nats of trailer elpd for 3.1 effective parameters (4.7 nats each), costs the catch stream 0.1 nats, adds no unreliable LOO points, and moves every calibration statistic toward nominal.**

**The daily boat AR buys elpd by spending 27 to 48 effective parameters on 195 observations, makes the catch stream 3.8 nats WORSE in both tau conditions, generates 2 then 16 Pareto k above 0.7, and destroys calibration.**

That is identification versus absorption, and it is no longer a hypothesis.

---

## 2. The trap the batch walked into, and why the new diagnostics exist

Read `pit_mean` alone and the ranking inverts:

```
                pit_mean   pit_sd   coverage_50      (OSP stream, private-boat all-gear)
ON x monthly       0.520    0.285        0.508
ON x DAILY         0.495    0.138        0.977
nominal            0.500    0.289        0.500
```

On the statistic the adequacy table shipped with, **the worst-calibrated fit in the batch is the best-calibrated cell.** Its PIT mean is near-perfect precisely because a latent process carrying one state per observation interpolates the data and piles every PIT value at 0.5. Ninety-eight percent of OSP observations fall inside a nominal fifty-percent predictive interval.

Both correcting statistics were already being computed into `ppc_calibration_*.csv` and neither reached the adequacy table. They do now (`cov50_worst_dev`, `pit_sd_worst_dev`, `flag_miscalibrated`), and the ordering they produce is the ordering every other line of evidence produces.

A second gap in the same table: **`disp_neff_min` asks whether a dispersion parameter SAMPLED well, which turns out to be a different question from whether it COLLAPSED.**

```
                sigma_r_OSP posterior median      its n_eff      gate floor
OFF x monthly                       0.786            12,463             400
ON  x monthly                       0.806            16,658             400
OFF x daily                         0.307               307             400   <- flagged
ON  x daily                         0.086               657             400   <- NOT flagged
```

The cell that squeezed the OSP observation error to a ninth of its monthly value passes the n_eff flag comfortably. `disp_scale_min` now reports the value. It is deliberately **not** flagged: the half-Cauchy priors on these scales have no finite SD, so the contraction diagnostic is undefined for them, and any fixed cut-off would be a number chosen after seeing these five runs. Compare the column across configurations, which is what it is for.

---

## 3. Stage by stage

### S1, the empty-stratum fill (plan 5.1). PASS, then a decision.

The holiday-loading fix is confirmed complete. The re-run reproduces the driver **exactly** on both quantities and all four components: zeroed days 44/0/47/9, PE catch 15,755 / 5,959 / 11,176 / 400. The 2026-08-28 stage A reported 41/44/8/0 and had even assigned the 44 to the wrong component.

| component | zeroed days | PE at `zero` | PE at `day_type` | change |
|---|---|---:|---:|---:|
| shore all-gear | 44 / 289 | 15,755 | 19,365 | **+22.9%** |
| shore pot closure | 0 / 76 | 5,959 | 5,959 | 0.0% |
| private boat all-gear | 47 / 289 | 11,176 | 13,750 | **+23.0%** |
| private boat pot closure | 9 / 76 | 400 | 514 | **+28.5%** |

At the port that is a PE of 45,112 going to **51,409, +14.0%**. Where a component reports its BSS the fill moves only the cross-check; where it falls back to PE it moves the reported number.

### S2, the informed-day floor (plan 5.2). PASS on all three criteria, and they are exact.

- **Shore fits bit-identical to stage D** (shared turnover globally OFF): 7,225 shared parameter rows at full double precision, convergence rows identical.
- **Boat fits bit-identical to stage E** (globally ON): 3,998 shared rows, identical.
- The floor's own bookkeeping agrees: shore all-gear 4 informed days with `shared_tau = 0`; boat all-gear 130 and boat pot closure 18, both with `shared_tau = 1`.

The feature is now boat-only by construction rather than by hope. Port total **71,521** [52,350, 100,759].

### S3, the gear-track cross-check (plan 5.2). PASS, and the departure from the plan is vindicated.

The plan prescribed more draws. The batch raised draws **and** adaptation to the pooled track's settings for the same component. The result is not marginal:

| | 2026-08-29 gear | S3 |
|---|---:|---:|
| `tau_bar` n_eff | 49 | **19,292** |
| `tau_bar` R-hat | 1.0813 | **0.9997** |
| divergences | 554 | **2** |
| gate | REJECTED | **BSS** |

Draws rose 2.5x; n_eff rose 394x. Draw count was not the binding constraint, adaptation was, and the plan's prescription alone would have produced another rejected fit.

**Cross-track agreement is now the strongest evidence in the project for the shared turnover:**

| | pooled S2 | gear S3 | gap |
|---|---:|---:|---:|
| `tau_bar` | 2.5969 | 2.5962 | **0.03%** |
| boat all-gear | 31,008 | 30,760 | -0.80% |
| shore all-gear | 20,898 | 20,754 | -0.69% |
| shore pot closure | 6,331 | 6,332 | +0.02% |
| boat pot closure | 1,018 | 956 | -6.09% |
| port total | 71,521 | 70,886 | **-0.89%** |

Two independently parameterized CPUE structures put the same turnover in the same place to three decimal places, and all four gear fits now report BSS.

**The one FAIL in the batch is the criterion's, not the code's.** It compared S3's shore fits against the 2026-08-29 gear run, which had `shared_tau` GLOBAL, so its shore fits carried a `tau_bar` the floor now correctly refuses; the reference differed in two ways at once. That is exactly the mistake the script calls out in stage F, committed in the script. Against the correct control, the 2026-08-26 shipped gear run with the turnover off for the shore in both, **S3's shore fits are bit-identical across 6,463 shared parameter rows.** The sampler override touched exactly the one fit it was aimed at. The criterion has been corrected in the runner.

### S4a, the matched (off, daily) cell (plan 5.3). PASS, and it returns a clean null.

Boat all-gear **37,359, identical to stage F's 37,359**, with shore all-gear 20,898 and shore pot closure 6,331 exactly matching stage D. So stage F's boat number was NOT contaminated by its other three components moving: fits are independent Stan runs, and `ar_escalate` settled the boat at daily on its first attempt, which is the same single fit a forced daily produces.

The insurance was worth buying and it paid out a null. Worth stating both halves: the objection to stage F as a control was methodologically right, and empirically it did not bite. What S4a adds beyond the number is the confirmation, on 20,898 and 6,331, that `ar_force` no longer leaks across sub-seasons.

### S4b, the (on, daily) corner (plan 5.3). ADDITIVE, and the least trustworthy fit in the batch.

Boat all-gear **42,344**, port **82,738**, `tau_bar` **3.232 [2.87, 3.63]**.

That `tau_bar` is the sharpest single argument against this cell. The independent OSP/trailer overlap calibration (`osp_trailer_overlap_calibration.csv`, n = 61, correlations 0.90-0.98) implies a turnover of **2.01 / 2.74 / 3.03** by three trailer metrics. At monthly the posterior sits at 2.597 [2.06, 3.25], straddling the middle of that range. At daily it moves to 3.232 with a **lower** bound of 2.87, so the whole 95% interval sits at or above the top of the external range, and it got *narrower* while doing so. A finer AR on the effort process has no business making the turnover more precisely identified; that precision is the latent process squeezing out the residual that turnover was being asked to explain.

The rest of the cell agrees. Both observation-error scales collapse (`sigma_r_OSP` 0.806 to 0.086, `sigma_r_E` 0.934 to 0.306), both observation dispersions run away (`r_OSP` 1.54 to 133.7 with a 95% upper bound of 43,875; `r_E` 1.15 to 10.6), and `phi_E` sits at 0.886 [0.82, 0.93], a near-random-walk latent process at daily resolution. A near-random-walk with the observation noise removed is a curve through the data.

**A mechanism worth naming, and its limits.** The reported effort is a sum of exponentiated latent states. The fitted stationary log-variance of that process, `sigma_eps_E^2 / (1 - phi_E^2)`, rises from roughly 1.9-2.3 in the monthly cells to 3.0-3.7 in the daily ones. Under a lognormal mean, that direction of change raises the expected sum for an unchanged level, which is a route by which the boat total can rise without any more crabbers existing. This is a mechanism **consistent with** the numbers, not a decomposition of them: the level `mu_mu_E` is refit simultaneously and these outputs do not separate the two contributions. Settling it would need a posterior decomposition of `E_sum` into level and variance contributions, which is one small addition to the diagnostics rather than a run.

### S5, the boat pot-closure AR (plan 5.6). PASS on both halves.

Pot closure **735 at biweekly** against the gear track's 743, a 1.1% gap where the pooled monthly fit gives 849. All-gear **25,868 at monthly**, exactly rung 4's value.

The second half is the one that mattered. The 2026-08-28 stage C reported the same 735, but its `ar_force` was honoured per population, so it silently forced the boat's all-gear sub-season to biweekly as well and moved that component to 28,893. **The contamination was worth 3,025 crab and 2,865 at the port** (68,853 against this run's 65,988). The pot-closure finding survives it unchanged; the port total it was reported alongside did not.

---

## 4. Where the numbers now stand

| run | configuration | port total | boat all-gear |
|---|---|---:|---:|
| S5 | shipped + boat pot closure biweekly | 65,988 | 25,868 |
| rung 4 | **shipped (authoritative)** | **66,237** | **25,868** |
| D | shipped, turnover machinery off | 66,094 | 25,868 |
| **S2** | **+ shared turnover, boat-only** | **71,521** | **31,008** |
| S3 | as S2, gear-resolved track | 70,886 | 30,760 |
| S4a | + boat AR daily | 77,486 | 37,359 |
| F | + escalate all fits to daily | 78,250 | 37,359 |
| S4b | + both | 82,738 | 42,344 |

Port spread **65,988 to 82,738, +25.4%**; boat component **25,868 to 42,344, +63.7%**. Both wider than the 2026-08-29 batch's 18% and 44%, and **every one of these fits passes every convergence-gate criterion with room to spare** (divergence impacts 0.0001 to 0.0056 posterior SD against a 0.10 threshold; `C_sum` n_eff 4,839 to 19,239 against a floor of 400). The gate is not being fooled by near-misses. It is answering a different question, confidently.

---

## 5. What this justifies changing, and what it does not

**Adopt the shared turnover, boat-only.** Five independent lines now agree and none dissent: it improves the stream it should and leaves the catch stream alone; it costs 3.1 effective parameters and no unreliable LOO points; it moves PIT mean, PIT sd and coverage all toward nominal; it replicates across two tracks to 0.03% on the parameter and 0.8% on the component; and it agrees with an external overlap calibration computed from different data. The pre-run objection (it moved the shore on four observations) is answered exactly by S2's bit-identity gate. **Recommended authoritative run: S2, port total 71,521 [52,350, 100,759].**

**Do not adopt the daily boat AR.** It fails on every axis that distinguishes a better model from a more flexible one.

**Two things this does not settle, stated so they are not quietly assumed.**

1. ~~*The trailer stream is over-covered in every configuration.*~~ **RETRACTED 2026-09-01. This was an artefact of the diagnostic, not a property of the model.** `ppc_calibration_*.csv` computed `coverage_50` by testing the observation against a QUANTILE INTERVAL of the simulated draws. For small counts that interval cannot carry 50% of the mass, so it over-covers by construction, and `ppc_byobs_*.csv` had always computed the same quantity correctly from the randomized PIT. The two disagreed by up to 0.154 on the private-boat trailer stream and were never compared. On the correct statistic the same runs give **0.523 (production) and 0.538**, which is 0.6 and 1.1 sampling SDs from nominal: calibrated. The catch stream, whose counts are an order of magnitude larger, is 0.595 either way, exactly the signature of a small-count artefact. `model_diagnostics.R` was fixed the same day so the two files agree.

   **The Stage 5 conclusion survives the correction and is sharper for it.** The daily-AR cells sit at 0.713 and 0.744 under the randomized statistic, 5.9 and 6.8 sampling SDs from nominal, against 0.523 for production. Removing the spurious background of over-coverage that affected every cell equally makes the contrast cleaner, and `flag_miscalibrated` now fires on the two daily cells and stays quiet on the two monthly ones, where before it fired on all five. A flag that fires on everything carries no information; at the time that was rationalised rather than treated as the tell it was.
2. *`crab_fraction_set = 0.3` is still a placeholder and still the largest single lever on this component.* A field question, not a code one.

---

## 6. Defects this review found

- **REGRESSION, fixed.** `tau_bar` was registered in the prior-vs-posterior table under its bare name while rstan summarises it as `tau_bar[1]` (it is declared `vector<lower=0>[shared_tau]`). `has_par()` strips the index and selected it, the `post[pn, ]` lookup then threw, and the writer's enclosing `tryCatch` swallowed the error, **silently dropping the entire boat `prior_vs_posterior_*.csv` in S2, S3 and S4b.** The entry is now named `tau_bar[1]` like `mu_mu_E[1]`, the row lookup tolerates a scalar/vector mismatch in either direction, and an unresolvable parameter is dropped rather than taking the file with it. The generalisable lesson: a `tryCatch` around a whole writer converts a one-row bug into a missing file.

  The contraction it should have reported, computed by hand: posterior sd 0.302 against a prior sd of 0.385, contraction **0.216**, so `prior_influential` would read TRUE. **That flag is misleading here and the definition is worth revisiting.** It measures variance reduction only. The posterior mean sits at 2.613 against a prior mean of 1.255, which is **3.5 prior SDs away**. A parameter that has moved three and a half prior SDs is not prior-dominated whatever happened to its variance.

- **`pit_worst_bias` alone ranks the cells backwards** (section 2). Fixed by adding `cov50_worst_dev`, `pit_sd_worst_dev` and `flag_miscalibrated`.
- **`disp_neff_min` cannot see a dispersion collapse** (section 2). `disp_scale_min` added, reported unflagged.
- **The S3 criterion compared against a reference differing in two ways** (section 3). Corrected in the runner.
- **`DRY_RUN` was left FALSE** after the batch. Reset to TRUE with the standing comment; `RESUME` would otherwise append five duplicate rows to the summary CSV on a re-source. The summary was checked and carries no duplicates.

---

## 7. Untested combinations, named so nobody assumes them

- **Shared turnover x biweekly boat pot closure.** S2 puts the pot closure at 1,018 (monthly AR, turnover on); S5 puts it at 735 (biweekly, turnover off). Both levers move that component and they have never been crossed on the pooled track. The gear track's S3 gives 956 at biweekly with the turnover on, against 743 at biweekly with it off, so the two look roughly additive there too, but the pooled cell is a guess until it is run. It is a small component; the honest cost of that guess is about 200 crab.
- **`day_type` fill x any of the above.** The fill changes only the PE, so it cannot move a component that reports its BSS. It moves the reported number only where a component falls back to PE, and in every Stage 5 run all four components reported BSS. It is therefore a decision about robustness and about the PE cross-check, not about the current headline.
- **`shared_tau_min_obs = 20`**, which would drop the boat pot closure (18 informed days) out of the feature. One run, and worth doing once as a sensitivity check on a threshold nobody should have to take on trust.
