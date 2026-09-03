# Review: the 2026-09-03 shore-AR / ZINB batch

**Date:** 2026-09-04
**Batch:** `06_diagnostics/run_shore_ar_zi_2026-09-03.R`, commit `daf299a` on `OSP-boat-count-incorporation`
**Runtime:** 23.6 h of fitting across five stages
**Verdict on the batch:** two stages are sound, three defects invalidated five of eleven verdicts, and the stage the batch existed for did not run at all.

---

## 1. Headline

| stage | model | min | port | shore all-gear | AR | boat all-gear | method |
|---|---|---:|---:|---:|---|---:|---|
| G1 | gear | 26.5 | 70,953 | 20,754 | monthly | 30,760 | BSS |
| Z0 | pooled | 263.9 | 71,450 | 20,898 | daily | 31,008 | BSS |
| A1 | pooled | 883.0 | 71,535 | 20,898 | daily | 31,008 | BSS |
| Z1 | pooled | 242.4 | 71,326 | 20,745 | daily | 31,008 | BSS |

**Nothing here changes the authoritative run.** `20260831/pooled-CPUE-VAL-1-adopted`, port 71,513, still stands. G1 restores a cross-check that had been unavailable since the shared turnover was adopted; P1 removes a documented inconsistency worth 0.059% on the boat; A1 produced nothing; Z1's evidence is real but was scored wrongly in both directions.

## 2. What actually worked

**P1, the PE/BSS incomplete-trip arm alignment. PASS.** Largest change -0.059% on `private_boat_all_gear`; the shore components did not move, which is the control (the shore PE has no gear-per-group ratio, so any shore movement would mean the shared helper reached further than intended). This was always a consistency fix, not an accuracy one: two arms of one fused estimator no longer disagree about which interviews count.

**G1, the gear-track boat all-gear sampler fix. PASS twice.** The fit that had never worked in production now samples: 2 divergences against 554, `tau_bar` n_eff 19,292 against 49, and the gear port total 70,953 against 51,385. That 51,385 was a 27% understatement of the cross-check that validates the pooled headline, caused by the boat component being rejected by the gate and falling back to PE. Cross-track agreement is restored at -0.78% on the port, against 0.73% for the 2026-08-26 shipped pair, and the two independently parameterized tracks agree on `tau_bar` to four decimal places (2.5962 against 2.5962). This is the strongest internal check the project has and it is working again.

**Z0, the ZINB behaviour-neutrality gate. PASS.** All four fits reproduce 2026-08-31 production across 11,223 shared parameter rows at full precision, confirming that `vector<lower=0,upper=1>[zi_catch] theta_C` is genuinely zero-size when the feature is off. The batch attached a caveat naming 9 "UNEXPECTED" config keys; that caveat was wrong and is resolved in section 4.

## 3. The stage that did not run

`ar_escalation_log.csv` for stage A1:

| fit | attempt | planned | AR resolution | P_n | pass gate | catch median | PI rel |
|---|---:|---:|---|---:|---|---:|---:|
| shore_all_gear | 1 | 4 | daily | 289 | TRUE | 20,897.819 | 0.2974 |
| shore_all_gear | 2 | 4 | daily | 289 | TRUE | 20,897.819 | 0.2974 |
| shore_all_gear | 3 | 4 | daily | 289 | TRUE | 20,897.819 | 0.2974 |
| shore_all_gear | 4 | 4 | daily | 289 | TRUE | 20,897.819 | 0.2974 |

One model, refitted four times, 14.7 h of the 23.6 h batch.

**Root cause (D1).** `01_BSS_models/BSS-GH-pooled-CPUE-model.Rmd` line 1202:

```r
force_res <- if (isTRUE(params$ar_escalate)) ar_ladder[attempt] else NULL
```

`ar_escalate` became SCOPED on 2026-09-04 (`.bss_resolve_ar_escalate()` accepts `TRUE`, a character vector of populations, or `list(shore = "all_gear")`), and this stage passed the list form. `isTRUE()` on a list is `FALSE`, so `force_res` stayed `NULL` on every attempt while `.esc_on` twelve lines above resolved correctly to `TRUE`. The consequence is specific and nasty: the ladder was BUILT correctly with four rungs, the loop ran the correct number of attempts, the log recorded `n_attempts_planned = 4`, and every observable except the resolution string said the ladder was working. Only the resolution failed to reach `prep_bss_crab_pooled()`, which then re-derived `daily` from the data-driven selector each time.

I introduced the scoped form and updated `.esc_on`, `bss_ar_ladder()` and the continuation logic in the same patch, and missed this one line. The identical defect was latent in the gear driver at line 819.

**A second, smaller defect in the same log (D1b).** `selected = TRUE` was written on attempt 4, but `.keep` had retained attempt 1. The post-loop marker matched on the resolution string and took the last hit, which mislabels the reported rung whenever two rungs share a resolution, i.e. exactly the state a broken ladder leaves the log in. The gear driver never set `selected` at all; it wrote a column of `FALSE`.

**A third error, in the batch's reading of its own log (D1c).** The A1-A3 verdict read "rung(s) whose shore all-gear fit was REJECTED and fell back to PE: daily, daily, daily". No rung was rejected; `pass_convergence = TRUE` on all four rows. The batch inferred rejection from per-rung adequacy columns that only the REPORTED fit ever writes, so absent diagnostics were read as a failed gate.

Both A1-A3 verdicts are **VOID**. The question they were asked is still open, and it is the most consequential open question in the batch.

## 4. Three diagnostics computed under the wrong model

### D2. The ZINB posterior-predictive check was computed as NB2

`pit_block()` in `save_run_diagnostics.R` and `calib()` in `model_diagnostics.R` computed the catch stream's PIT, coverage flags and `p_zero` from `dnbinom`/`pnbinom` alone. The fit was a mixture. `theta_C` was excluded entirely.

| shore all-gear catch stream | observed zeros | expected | z |
|---|---:|---:|---:|
| NB2 baseline (Z0) | 676 | 605.5 | +3.8 |
| **as the batch reported it (NB2 arithmetic on a ZINB fit)** | 676 | 419.0 | **+15.4** |
| **under the mixture** | 676 | 637.9 | **+2.0** |

The batch recorded that the feature made the zero bin dramatically worse. It halves it. The shore pot-closure replicate closes outright: 146 against 138.1, z = +0.8, from 118.8 / z = +2.9.

The same defect makes Z1's catch-stream `pit_mean` (0.4325 against the baseline's 0.4984), `pit_sd`, `coverage_50` (0.405 against 0.457) and the `flag_pit_bias = TRUE` it raised on that fit unreadable. The mechanism is easy to state: with `theta` absorbing the excess zeros, `lambda_C` rises to fit the non-inflated component, and scoring the observed zeros against the risen NB2 mean alone drags the PIT down. `elpd_loo` and every Pareto k are unaffected, because Stan's `log_lik` carried the mixture from the start.

The mixture PIT is **not** the obvious expression. The randomized PIT is `F(y-1) + 0.5 f(y)`, and under the mixture `F(-1) = 0` while `F(y)` for `y >= 0` carries the whole point mass:

```
y = 0 :  0.5 * (theta + (1 - theta) * NB2_pmf(0))
y > 0 :  theta + (1 - theta) * (NB2_cdf(y-1) + 0.5 * NB2_pmf(y))
```

Folding the two into one expression puts the zero observations about `theta` too high, flattering exactly the observations the parameter was introduced to explain. The harness now asserts against that specific mistake.

### D3. The elpd criterion used the wrong standard error

The batch compared the elpd difference (+14.8 nats) against `se_elpd_loo` from each fit's own `loo_summary` (about 46 nats) and recorded NOT WORTH IT.

`se_elpd_loo` is the standard error of ONE model's total elpd. It is dominated by variation ACROSS observations, i.e. by the fact that some interviews are easier to predict than others. That variation is common to both models and cancels in the difference. The model-comparison statistic is the standard error of the PAIRED difference, `sd(d) * sqrt(n)` (Vehtari, Gelman & Gabry 2017, *Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC*, Statistics and Computing 27:1413-1432, section 3.3; this is what the `loo` package's `loo_compare` reports as `se_diff`).

| stream | elpd A -> B | diff | paired SE | ratio | naive SE |
|---|---|---:|---:|---:|---:|
| shore all-gear | -3,162.5 -> -3,147.7 | +14.8 | 5.5 | **2.69 SE** | 46.3 |
| shore pot closure | -1,450.0 -> -1,440.9 | +9.1 | 4.6 | 1.96 SE | 21.8 |

The gain clears the batch's own stated 2 SE bar, and the independent replicate on a different data subset lands at 1.96 SE with `theta_C` 0.107 against 0.179. The verdict should have been WORTH IT.

**Read the decomposition before adopting anything.** The gain is bought entirely at the zeros:

| stream | zeros | positives |
|---|---|---|
| shore all-gear | n = 676, **+25.2** (11.6 SE) | n = 973, **-10.5** (-2.1 SE) |
| shore pot closure | n = 146, +19.0 (8.9 SE) | n = 481, -10.0 (-2.6 SE) |

A mixture that improves the zeros by degrading the positive counts is a different object from one that fits better everywhere, and the harvest estimate is made of the positives. Pareto k above 0.7 on the shore all-gear fit did fall 41 -> 26, which is a real improvement computed from the correct likelihood.

### D4. `run_parameters.txt` has been silently truncated in every run

`str()` caps a list at `list.len = 99`; `run_config` carries 120 keys. Every `run_parameters.txt` in this repo written before 2026-09-04 is 105 lines ending in `[list output truncated]` and is missing roughly 21 keys.

That file is not decoration. `config_delta()` in the batch runners compares two runs key by key from it, and `annotate_decoupled_run.R` reads run flags out of it. A truncated dump produces false positives (a key that crossed the cut between two runs reads as differing) and false negatives (a key past entry 99 in both dumps is never compared at all).

The Z0 caveat naming 9 "UNEXPECTED" keys was exactly this: `ar_escalate_stop` is present in the new dump and absent from the old; `opener_covariate_mode` is present in the old and absent from the new. Nothing about the config differed. Both are 105 lines.

There is one live consequence beyond the false caveat. `annotate_decoupled_run.R` reads `opener_covariate_mode` to decide whether `B_open` is decoupled, and that key fell past the cut in every folder from 2026-09-01 onward. Post-hoc annotation of those folders would have left `B_open` unflagged, the same class of silent missing-flag as the unflagged `R_G_boat` raised on 2026-08-25. The in-driver annotation reads `params` directly and is not affected, so no committed `structural_params_*.csv` is wrong.

## 5. What has been fixed

| id | fix | file |
|---|---|---|
| D1 | `force_res` resolves the scoped `ar_escalate` via `.esc_on` | both drivers |
| D1b | the reported rung is tracked by ATTEMPT INDEX, not by resolution string; the gear driver now marks `selected` at all | both drivers |
| D1c | the report-table gate reads the per-fit `escalation_enabled` column instead of `isTRUE()` | pooled driver |
| D2 | mixture PIT / `p_zero` in one shared helper both diagnostic files call | new `03_R_functions/zinb_ppc.R` |
| D3 | paired-difference SE with a zero/positive decomposition | new `03_R_functions/loo_elpd_paired.R` |
| D4 | `str(params, list.len = 1000, ...)`; 120 of 120 keys captured. `annotate_decoupled_run()` warns on a truncated dump; `config_delta()` refuses to raise UNEXPECTED against one | both drivers, `annotate_decoupled_run.R`, new `03_R_functions/batch_verdict_helpers.R` |
| - | the injected `[A1] PASS fake fitted criterion` test row is removed from the committed verdicts | `05_output/shore_ar_zi_2026-09-03_verdicts.csv` |
| - | `run_shore_ar_zi_2026-09-03.R` shipped `DRY_RUN <- FALSE`; reset to TRUE | `06_diagnostics/` |
| new | posterior draws the PPC reads are persisted per fit, so the next diagnostic fix is a recomputation and not a re-fit | new `03_R_functions/save_bss_ppc_draws.R`, `save_ppc_draws = TRUE` |

The committed verdicts file now carries `verdict_2026_09_04` and `correction_2026_09_04` columns beside the originals, so the supersession is auditable rather than rewritten.

The harness (`06_diagnostics/test_improvements_2026-08-25.R`) is at **330 assertions, 0 failures**, including source-level assertions that no driver tests `ar_escalate` with `isTRUE()`, that a scoped ladder yields distinct rungs, that the mixture PIT uses `F(-1) = 0` at `y = 0` and is NOT the naive folded form, that `theta = NULL` reproduces the NB2 arithmetic exactly, that the paired SE detects an effect the naive SE misses, and that the config dump captures every key.

## 6. Why this happened, and the one structural change

Three of the four defects share a shape: **a change was validated by the mechanism it was about to break.** The ladder was tested by a harness that checked `bss_ar_ladder()`'s output and never that the resolution reached the data prep. The ZINB feature was validated by a zero-bin diagnostic written before the mixture existed. The config comparison was validated by reading the file it was silently truncating.

The structural change is `save_ppc_draws`. Three defects in five weeks (the quantile-interval coverage defect on 2026-08-31, the `tau_bar` row-key defect on 2026-09-01, this one) were all in code that runs AFTER sampling, on a fit that was never wrong, and each cost a multi-hour re-fit to correct. Nothing persisted a stanfit or any subset of its draws. With the draw objects on disk, tens of MB per fit, the next such fix is a recomputation.

## 7. Where the two open questions stand

**The shore all-gear AR resolution: OPEN, unchanged since 2026-08-31.** Production fits this component, 29% of the port total, at a daily AR with p_loo at 35.2% of `n_obs`, 41 Pareto k above 0.7, and `coverage_50` 0.701 (+7.1 sampling SDs). The gear track fits it at monthly and gets `coverage_50` 0.035 (-16.4 SDs). Both ends are bad, so the answer is probably between them, and 14.7 h of fitting did not move this one inch.

**ZINB on the shore catch stream: the evidence is positive but the decision needs rendered diagnostics.** +14.8 nats at 2.69 paired SE, replicated at 1.96 SE on an independent fit, the zero bin halving from z = +3.8 to +2.0, Pareto bad k 41 -> 26. Against it: the gain is entirely at the zeros while the positives get worse by 2.1 SE, and the reported shore total moves -0.7%. The corrected numbers above were recomputed OFFLINE from committed marginal means and carry a mean-product approximation in the zero bin (`E[theta * p0]` taken as `E[theta] * E[p0]`). That approximation is very unlikely to move z by more than about 0.1, but "very unlikely" is not a basis for changing the likelihood of a harvest estimate. `estimate_catch_zi` stays FALSE until a run with the corrected PPC renders the exact values.

## 8. Next steps

`06_diagnostics/run_ladder_zinb_2026-09-04.R`, roughly 3-6 h against the 23.6 h it corrects.

- **D0 (desk, free, runs in a dry run).** Re-derives every corrected number above from committed files, so the corrections are reproducible by anyone with the repo rather than asserted in a document.
- **L1 (about 2-3 h).** The shore all-gear ladder, for real. `ar_escalate = list(shore = "all_gear")` with `all_rungs`, exercising the production toggle rather than a batch-only override. Two switches are exposed at the top of the file with their tradeoffs written out:
  - `LADDER_INCLUDE_DAILY` defaults FALSE, taking the daily rung from Z0 (whose config differs only in `ar_escalate` keys, none of which touch the model, data or seed) and saving 3.7 h. Set it TRUE if the ladder will be shown outside the project, because a cross-folder claim is weaker than a same-folder one.
  - Decision rule, as agreed from the FWC discussion: escalate on the BSS convergence gate, report the finest rung that PASSES, use the narrowest relative prediction interval only as a tie-break. Carry the Stage 5 caution: on the boat 2x2 the two MISCALIBRATED cells had the narrowest relative intervals, so any rung a PI tie-break selects needs `cov50` and the Pareto k count checked before it is believed.
- **Z2 (about 4 h, `ZINB_RERENDER`).** Re-renders the ZINB shore fit under the corrected PPC so the adoption decision rests on rendered output, checks the rendered zero bin against D0's approximation (a gap above 0.2 in z would mean the theta/`p0` covariance is not negligible), produces the first readable calibration for a ZINB fit in this project, and confirms the boat negative control still reproduces bit-for-bit.

Both fitted stages will write `ppc_draws_*.rds`, so the next diagnostic defect costs a recomputation rather than a batch.

---

## 9. Correction, 2026-09-05: where the ZINB actually loses

Section 4's "the positives get worse" reading was misleading and is corrected here rather than rewritten above.

| shore all-gear, elpd change by count | y=0 | y=1 | y=2 | y=3-4 | y=5-8 | y=9-16 | y=17+ |
|---|---:|---:|---:|---:|---:|---:|---:|
| n | 676 | 236 | 161 | 239 | 230 | 91 | 16 |
| diff (nats) | **+25.2** | **-42.0** | -6.6 | **+14.3** | **+18.1** | **+5.6** | +0.1 |
| SE ratio | 11.6 | **-16.7** | -4.5 | 10.9 | 11.8 | 4.6 | 0.2 |
| share of catch | 0% | 6% | 8% | 20% | 34% | 25% | 8% |

Every bin from 3 upward improves, and those bins carry 87% of the shore catch. The loss is at y=1, and it is enormous. `r_C` doubles (0.95 -> 1.83): once `theta` absorbs structural zeros the NB2 no longer needs extreme overdispersion to reach zero, so it tightens, fitting the harvest-carrying counts better and the almost-zero count of 1 worse. The pot-closure replicate has the same shape (y=1 -20.4; `r_C` 1.68 -> 2.86).

So the ZINB has **moved** the misfit from the 0 bin to the 1 bin, not removed it. Excess mass at both 0 and 1 relative to NB2 is the signature of a two-regime process, unsuccessful trips yielding 0-1 crab against successful ones, which a hurdle or a two-component NB mixture fits and a ZINB by construction cannot. The zero-bin check alone, which is all stage Z2 originally carried, would have passed while recording none of this.

Changes: `ppc_byobs_*.csv` carries `p_one` beside `p_zero` (mixture-aware, `bss_zi_p_k()`); `loo_elpd_paired()` returns the by-count table; stage Z2's criterion is pre-set as a count-bin table that passes as a whole or not at all, with the named alternative if the one bin fails. Harness 335.
