# Review: the 2026-09-04 ladder / ZINB batch

**Date:** 2026-09-06
**Batch:** `06_diagnostics/run_ladder_zinb_2026-09-04.R`, results at `dabe4a5` on `OSP-boat-count-incorporation`
**Runtime:** 9.1 h of fitting (L1 323.6 min, Z2 223.4 min)
**Verdict:** both stages did what they were built to do. The ladder ran for real and the shore all-gear AR question is now answerable; the ZINB re-render corrected two of my own readings, one of which was wrong in the ZINB's favour.

---

## 1. The ladder, at last

| rung | P_n | divergences | frac | catch (median) | 95% CI | PI rel | effort | gate |
|---|---:|---:|---:|---:|---|---:|---:|---|
| **daily** (from Z0) | 289 | 333 | 0.033 | 20,898 | [18,187, 24,402] | 0.2974 | 24,997 | PASS |
| **weekly** (reported) | 44 | 313 | 0.0313 | 21,547 | [19,076, 24,528] | 0.2530 | 25,567 | PASS |
| **biweekly** | 21 | 212 | 0.0212 | 21,383 | [19,143, 24,041] | 0.2291 | 25,258 | PASS |
| **monthly** | 10 | 152 | 0.0152 | 20,771 | [18,628, 23,217] | 0.2209 | 25,092 | PASS |

Four distinct resolutions, four distinct `P_n`, four distinct estimates. The mechanism works.

**Finding 1: the AR resolution is not a large lever on the shore point estimate.** The four rungs span 20,771 to 21,547, a range of 777 crab, **3.7%** of the component and about 1% of the port total. Compare the Stage 5 boat 2x2, which spanned 44%. The worry that motivated this whole line of work, that the daily AR might be materially distorting a component worth 29% of the port, is not borne out at the level of the point estimate.

**Finding 2: the convergence gate separated nothing.** All four rungs pass. Under the agreed production rule, "report the finest rung that PASSES", the answer is simply the finest rung run, whatever its adequacy. That is the limit recorded in Section 1h of `PIPELINE_STATUS.md`, now demonstrated rather than predicted.

**Finding 3: the adequacy statistics separate the rungs decisively, and they say daily is overfitted.**

| statistic | daily | weekly |
|---|---:|---:|
| `p_loo` as a fraction of `n_obs` | **0.352** | **0.096** |
| `p_loo`, gear stream | 109.4 | 29.7 |
| Pareto k > 0.7 | **41** of 311 | **1** |
| gear `coverage_50` (nominal 0.500) | 0.701 (**+7.1 SD**) | 0.559 (**+2.1 SD**) |
| catch `coverage_50` | 0.457 (-3.5 SD) | 0.441 (-4.8 SD) |
| `flag_miscalibrated` | **TRUE** | **FALSE** |

The daily fit was spending about one effective parameter per three observations on the effort stream. At weekly it spends 80 fewer, its coverage moves from +7.1 sampling SDs to +2.1, and the miscalibration flag clears.

**Finding 4: `elpd` favours daily, and should not be believed.**

| stream | daily -> weekly | paired SE | ratio |
|---|---:|---:|---:|
| gear | -1251.3 -> -1337.1 (**-85.8**) | 11.7 | -7.34 SE |
| catch | -3162.5 -> -3170.0 (-7.4) | 3.6 | -2.08 SE |

Daily predicts better on both streams. Two things disqualify that as an argument for daily:

1. The gear-stream gain is 85.8 nats for 79.7 additional effective parameters, i.e. **1.08 nats per parameter**. The batch's own written criterion is that "a rung that buys elpd by spending effective parameters at about one nat each is buying noise". This is that, almost exactly.
2. **41 of the daily fit's 311 gear observations have Pareto k above 0.7.** PSIS-LOO has failed for 13% of the stream, so the daily fit's own `elpd_loo` is not a reliable estimate of anything. The weekly fit has 1. Comparing the two on elpd favours the fit whose elpd cannot be trusted.

**What this does NOT establish.** Only two of the four rungs have adequacy at all. Biweekly and monthly were fitted, cost 174 minutes, and left no `p_loo`, no Pareto count and no coverage, because only the fit the ladder KEEPS reaches `write_bss_diagnostics()`. So the ladder shows that daily is too fine and weekly is much better behaved; it does not show where coarsening stops helping. The gear track fits this component at monthly and gets `coverage_50` 0.035 (-16.4 SD), badly under-covering, so the turnover is somewhere between weekly and monthly and remains unlocated. **Fixed for future runs**: `bss_rung_adequacy()` now writes those statistics into every rung's row of `ar_escalation_log.csv`.

**Effect on the reported number.** Port total 72,122 [53,228, 101,553] at weekly against 71,450 at daily in the same configuration, **+0.94%**, with the lower bound rising from 52,398 to 53,228. Against the current authoritative run (71,513) it is +0.85%.

**Runtime.** My estimate of 2-3 h for L1 was wrong; it took 5.4 h. The rungs cost 95.6, 91.1 and 83.3 minutes: **coarsening the AR barely reduces the cost**, so the latent AR dimension is not what makes these fits slow. Any future ladder should be budgeted at roughly one full fit per rung.

## 2. The ZINB, re-rendered, and two corrections to my own readings

Z2 reproduces Z1 exactly (`theta_C` 0.1780, `C_expected_sum` 20,872.6 in both), so the re-render control passes and everything below is attributable to the corrected diagnostics rather than to a different fit.

**Correction A: the offline approximation was fine.** The zero bin renders at 676 observed against **638.5** expected, z = **+2.0**. My offline `E[theta]*E[p0]` approximation gave 637.9 and z = +2.0: accurate to 0.6 crab and to the first decimal in z. The pot-closure replicate likewise (138.6 / +0.7 against 138.1 / +0.8). That approximation can be trusted in future, which matters because it costs nothing.

**Correction B: "the misfit migrated to y = 1" was wrong.** That reading came from the elpd split and could not be checked, because `Z0` predates the `p_one` column. This run gives the comparison:

| fit | likelihood | AR | y=1 observed | expected | z |
|---|---|---|---:|---:|---:|
| L1 shore all-gear | NB2 | weekly | 236 | **335.4** | **-6.1** |
| Z2 shore all-gear | ZINB | daily | 236 | **285.0** | **-3.2** |

The NB2 over-predicts ones by 6.1 SD; the ZINB by 3.2. **The ZINB roughly halves the one-bin misfit as well as the zero-bin misfit** (+3.8 to +2.0). It is not moving the problem, it is reducing both ends of it. Caveat: the two fits differ in AR resolution as well as likelihood, so this is not a clean like-for-like; a daily NB2 fit with `p_one` does not exist because Z0 predates the column and `LADDER_INCLUDE_DAILY` was FALSE.

**Why the y = 1 elpd loss is not evidence against the model.** The loss is -42.0 nats over 236 observations. A zero-inflated mixture multiplies every non-zero probability by (1 - theta): log(1 - 0.178) x 236 = **-46.3 nats**. The observed loss is that mass tax, not a failure of shape. The tax falls on all 973 non-zero observations (-190.7 nats if nothing else changed), and the tightening of `r_C` from 0.947 to 1.828 repays it with interest at y >= 3, which is where 87% of the catch is. The aggregate is +14.8 nats at 2.69 paired SE.

**Calibration, readable for the first time.**

| run | catch `coverage_50` | PIT mean |
|---|---:|---:|
| Z0, NB2 | 0.4572 (-3.5 SD) | 0.4984 |
| Z1, **computed under the wrong likelihood** | 0.4051 (-7.7 SD) | 0.4325 |
| Z2, ZINB correctly scored | **0.4645 (-2.9 SD)** | 0.4965 |

Under the correct arithmetic the ZINB's catch calibration is slightly **better** than the NB2's, where the 2026-09-03 batch reported it as much worse and raised `flag_pit_bias` on the strength of it.

**What the ZINB does not do.** Z2's shore all-gear adequacy is still `p_loo_frac` 0.3437, 26 bad Pareto k, `flag_miscalibrated` TRUE, because Z2 runs at daily AR. Zero-inflation is a catch-stream change and does not touch the effort-process overfitting; that is what the ladder addresses. **The two have never been run together.**

## 3. Defects in the run

| what | consequence | status |
|---|---|---|
| The runner never moved the rendered HTML. `rmarkdown::render()` writes it beside the `.Rmd`; the 2026-09-03 runner moved it afterwards and this one did not. L1 rendered its report to `01_BSS_models/BSS-GH-pooled-CPUE-model.html` and **Z2 overwrote it**. | **L1's HTML report is destroyed** and can only be regenerated by refitting (5.4 h). The AR ladder table it tabulates is the artefact the FW creel discussion specifically asked to be preserved. The numbers survive in `ar_escalation_log.csv`, so the loss is cosmetic. A 6,062-line render artefact was also committed inside `01_BSS_models/`. | Fixed: the runner moves it and warns if the move fails; `01_BSS_models/*.html` is gitignored; the stray file is removed; the harness asserts both over every runner. |
| Only the kept rung got adequacy. | Two rungs, 174 minutes, no evidence. | Fixed: `03_R_functions/bss_rung_adequacy.R`, `ar_rung_adequacy = TRUE`. |
| `fmt()` was scalar-only (`if (length(x) == 0 \|\| is.na(x))`) and `verdict_L1` passes it whole ladder columns. | On **R >= 4.3 this is an error, not a warning**, and the batch aborts in the verdict block AFTER 9.1 h of fitting. It survived here only because the run used an older R. | Fixed: vectorised. |
| `DRY_RUN <- FALSE` committed. | Same class as 2026-08-25 and 2026-09-03; a fresh clone starts real fits. | Reset to TRUE. |
| `ppc_draws_*.rds` are gitignored. | The 4 files written by each stage stay on the machine that ran the fit and did not travel with the push. Not a defect, but as of today **no recompute path reads them back**, so they remain write-only. | Verdict now says so instead of reporting a misleading zero. |

Harness: **372 assertions, 0 failing.**

## 4. Where the two decisions stand

**Shore all-gear AR resolution.** The evidence supports moving off daily. Every adequacy statistic improves sharply at weekly, the point estimate moves +3.1% on the component and +0.94% on the port, and the only statistic favouring daily is an elpd whose own reliability diagnostic has failed on 13% of the stream. What is missing is the other half of the bracket: biweekly and monthly have no adequacy, so "weekly" is currently "the finest rung that is clearly not overfitted", not "the best rung". A short re-run with `ar_rung_adequacy` on would settle it, and would cost about 4.5 h for the three coarser rungs.

**ZINB on the shore catch stream.** Every objection I raised has now been answered by rendered output: the zero bin halves, the one bin halves, the catch calibration improves, the elpd gain is 2.69 paired SE with an independent replicate at 1.96, and the y=1 elpd loss is the arithmetic of the mixing weight rather than a shape failure. `estimate_catch_zi` still ships FALSE for one reason only: it has never been run at the AR resolution the ladder now points to, and adopting a likelihood change and a resolution change in one step would leave neither attributable.

## 5. Next steps

**One run settles both.** `L1` at weekly and `Z2` at daily each moved the shore component in the same direction for different reasons, and neither has been tested in the other's presence.

1. **C1, the combined candidate** (about 4.5 h): shore all-gear at weekly with `estimate_catch_zi = TRUE` on shore, boat untouched as the negative control. This is the candidate production configuration. Read it against L1 (weekly, NB2) to isolate the ZINB at the new resolution, and against Z2 (daily, ZINB) to isolate the resolution under the new likelihood.
2. **C2, the missing half of the bracket** (about 4.5 h, can share the run): the ladder again with `ar_rung_adequacy = TRUE` so biweekly and monthly finally report `p_loo`, Pareto k and coverage. Without this the choice of weekly rests on it being the finest non-overfitted rung rather than on it being the best one.
3. **Then re-baseline.** If C1 holds, one clean production render supersedes `20260831/pooled-CPUE-VAL-1-adopted` and carries both changes, with `ar_max_resolution$pooled$shore$all_gear` set to whatever C2 supports.
4. **Write the recompute path** for `ppc_draws_*.rds`, so the next diagnostic defect is a recomputation rather than a re-fit. It is the only part of that feature still missing, and the feature has now paid for itself twice on paper and zero times in practice.

Unchanged and still open: the OSP crabbing-only column (a data request to OSP), more boat interviews next season (field protocol), and the stale weather module.
