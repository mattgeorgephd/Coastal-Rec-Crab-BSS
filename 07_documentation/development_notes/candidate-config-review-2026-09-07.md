# Review: the C1 / C2 run, and the candidate production configuration

**Date:** 2026-09-07
**Batch:** stages C1 and C2 of `06_diagnostics/run_ladder_zinb_2026-09-04.R`, results at `542bd06`
**Runtime:** 9.3 h (C1 3.8 h, C2 5.5 h)
**Verdict:** the 2x2 is complete and it is clean. C1 (shore all-gear at weekly, ZINB on the shore catch) is the best-behaved fit this project has produced, and the two changes do not interact.

---

## 1. The 2x2

| configuration | shore all-gear | port total | `p_loo_frac` | Pareto k>0.7 | `cov50` dev | miscalibrated |
|---|---:|---:|---:|---:|---:|---|
| daily + NB2 (production) | 20,898 | 71,450 | 0.3519 | **41** | 0.2042 | **TRUE** |
| weekly + NB2 (L1) | 21,547 | 72,122 | 0.0956 | 1 | 0.0691 | FALSE |
| daily + ZINB (Z2) | 20,745 | 71,287 | 0.3437 | 26 | 0.2010 | **TRUE** |
| **weekly + ZINB (C1)** | **21,489** | **72,032** | **0.0949** | **0** | **0.0691** | FALSE |

**The effects are additive.** Interaction on the shore all-gear component: (21,489 - 21,547) - (20,745 - 20,898) = **+95 crab**, 0.4% of the component. The AR resolution moves the effort process and the ZINB moves the catch likelihood, and on this data they do not talk to each other, which is what the arithmetic predicted and what makes the separate results from L1 and Z2 transfer to the combination.

**The division of labour is equally clean.** The AR resolution does essentially all the work on `p_loo` and coverage (0.352 to 0.095; 0.204 to 0.069), and the ZINB does the remaining work on the Pareto diagnostics (41 to 26 at daily, 1 to **0** at weekly). C1 is the only fit in this project with **zero** Pareto k above 0.7.

**Port total 72,032 [53,044, 101,212]**, +0.73% on the current authoritative 71,513, with the lower bound rising from 52,310 to 53,044.

## 2. A correction I have been carrying since 2026-08-31

I have repeatedly written that the shore all-gear AR is bracketed by two bad ends: daily at `coverage_50` 0.701 (+7.1 SD) and monthly at 0.035 (-16.4 SD). **The second figure is from the gear-resolved model, not the pooled one.** It is a different model with a different likelihood (per-gear CPUE, `G = 5` on the shore), and it has no bearing on where the pooled model's AR belongs.

The pooled model's own bracket, measured by C2:

| rung | `cov50` gear | `cov50` catch | total `p_loo` | Pareto k>0.7 | catch | PI rel |
|---|---:|---:|---:|---:|---:|---:|
| daily | 0.7010 (**+7.1 SD**) | 0.4572 (-3.5 SD) | 153.6 | 41 | 20,898 | 0.2974 |
| weekly | 0.5659 (+2.3 SD) | 0.4433 (-4.6 SD) | 54.8 | 1 | 21,547 | 0.2530 |
| biweekly | 0.5659 (+2.3 SD) | 0.4500 (-4.1 SD) | 34.2 | 0 | 21,383 | 0.2291 |
| monthly | 0.5788 (+2.8 SD) | 0.4506 (-4.0 SD) | 22.2 | 0 | 20,771 | 0.2209 |

**Only daily is an outlier.** Weekly, biweekly and monthly sit within 0.5 sampling SD of each other on effort coverage and span 776 crab, 3.6% of the component. There is no bad far end in the pooled model, and the ladder's message is narrower and cleaner than I have been stating it: *daily is overfitted; everything coarser is adequate.*

**Which means the choice among the three coarse rungs is not made by the evidence.** `p_loo` keeps falling with coarseness (parsimony favours monthly), the catch stream's coverage improves very slightly with coarseness, the effort stream's degrades very slightly, and the prediction interval narrows monotonically, which is the tie-break the Stage 5 boat 2x2 showed to be a trap. Under the rule agreed at the FW creel meeting, "escalate on the convergence gate and report the finest rung that PASSES", the answer is **weekly**, and that rule is defensible here precisely because weekly is adequate on every statistic rather than merely passing a gate that everything passes.

**One statistic is not fixed by any rung.** The catch stream is under-covered at every resolution, -3.5 to -4.6 SD, and coarsening the effort process makes it marginally worse. That is a catch-likelihood problem, and it is what the ZINB addresses.

## 3. The ZINB, finally compared like for like

C1 against L1 is the comparison every previous run was missing: same resolution, same data, one likelihood difference.

| stream | bin | NB2 (L1) | ZINB (C1) |
|---|---|---|---|
| shore all-gear | zero | 607.2 expected, **z = +3.7** | 638.4, **z = +2.0** |
| shore all-gear | one | 335.4 expected, **z = -6.1** | 285.7, **z = -3.3** |
| shore pot closure | zero | 118.4, z = +2.9 | 138.4, **z = +0.7** |
| shore pot closure | one | 114.0, z = -3.5 | 91.1, **z = -1.2** |

The ZINB halves both bins on all-gear and closes both on the pot closure. The 2026-09-05 reading that "the misfit migrated to y = 1" is definitively refuted; the 2026-09-06 correction is confirmed on a clean comparison.

`elpd` at weekly: **+11.6 nats at 2.30 paired SE** (it was +14.8 at 2.69 SE at daily). The gain shrinks slightly when the effort process is no longer overfitted, which is the right direction and the expected size, and it still clears the 2 SE bar. So the gain was not an artefact of the daily overfitting absorbing structure the catch likelihood should have explained.

**What still does not close.** The shore all-gear one bin remains at z = -3.3, so the batch's own criterion (both bins under 2.5) records REVIEW. The data has more zeros AND fewer ones than even the mixture predicts: it is more bimodal than a zero-inflated NB can be, which is the two-regime signature. A hurdle model, or a two-component NB mixture with its own mean in the low regime, is the shape that fits it. The available further gain is bounded: the one bin is 6% of the catch and the zero bin is 0%.

## 4. Four defects, all in code I wrote in the previous patch

| what | consequence | status |
|---|---|---|
| `bss_rung_adequacy()` reported `p_loo_frac` as the SUM over streams divided by summed `n_obs`. `bss_model_adequacy.R` reports the **worst stream's** fraction. | The C2 ladder printed "2.8% of n_obs" for a rung whose `model_adequacy.csv` figure is 9.6%: two different statistics under one name, in one run, being compared with each other. My README claim that they were "directly comparable" was false. | Fixed; reports the worst stream and names it. Verified against a real fit. |
| `bss_ppc_calibration()` calls `set.seed()` unconditionally, and both it and `save_run_diagnostics()` take a random subset of draws. Calling the rung helper inside the fitting loop shifted the global RNG for every diagnostic afterwards. | C2's kept fit is **bit-identical to L1's across 10,251 parameter rows** and still reported gear `coverage_50` 0.5659 against L1's 0.5595. A diagnostic-only toggle moved a committed diagnostic. Bit-comparability is the tool that has caught five defects in this project. | Fixed: the RNG state is saved and restored. Verified. |
| `verdict_L1(dirs$C2)` wrote C2's ladder rows under the stage label `"L1"`, and `merge_csv_by()` de-duplicated stale rows against the new frame but never within it. | The verdicts file carried **two rows per criterion**, the first ladder's reading and the second's, indistinguishable to a reader. | Both fixed: the verdict block takes a stage id, and `merge_csv_by()` de-duplicates within the new frame (last write wins). The committed file is relabelled and annotated. |
| `DRY_RUN <- FALSE` committed. | Fourth occurrence (2026-08-25, 09-03, 09-06, now). | Reset. The harness has asserted this since 2026-09-03 and catches it every time; what it cannot do is stop the commit. |

The C2 per-rung numbers in the committed run carry the superseded `p_loo` definition and the RNG shift. They are annotated in the verdicts file rather than silently corrected, because they cannot be recomputed from what was written. **Read those three rungs against each other, not against `model_adequacy.csv`.**

Harness: **378 assertions, 0 failing.**

## 5. Recommendation

**Adopt C1.** Shore all-gear at weekly AR, `estimate_catch_zi = TRUE` scoped to shore. It is the best-behaved fit in the project on every adequacy statistic, the two changes are independent, and the evidence for each was collected separately and then confirmed in combination.

The adoption is a config change, not a code change:

```r
ar_max_resolution$pooled$shore$all_gear  <- "weekly"   # currently data-driven -> daily
estimate_catch_zi                        <- TRUE
catch_zi_populations                     <- c("shore")
```

C1 itself used `ar_force`, which is an experiment lever that bypasses the cap and the selector; production must express the same thing through `ar_max_resolution`. That is a one-line difference and it needs one confirming render, which doubles as the new authoritative run.

**What a reviewer will ask, and the honest answers.**

- *Why not monthly, which is more parsimonious and has a narrower interval?* Because narrower intervals from coarser latent processes are what the Stage 5 boat 2x2 showed to be false precision, and because the agreed rule is the finest adequate rung. The three coarse rungs differ by 3.6% and nothing in the data distinguishes them.
- *Why is the port total higher than the last one?* +0.73%, and almost all of it is the shore effort process no longer being over-imputed at daily. Nothing has been published, so this supersedes an internal working number, not a released one.
- *Is the catch likelihood right now?* Better, not right. The shore all-gear one bin is still 3.3 SD out and the catch stream is under-covered at every resolution. A hurdle model is the named next step and its available gain is bounded by about 6% of the catch.

## 6. Next steps

1. **P1, the adoption render** (about 4 h). `ar_max_resolution$pooled$shore$all_gear = "weekly"` plus `estimate_catch_zi = TRUE`, no `ar_force`. Gate: it must reproduce C1's four fits bit-identically, because expressing the same resolution through the cap rather than the override is a routing change and nothing else. If it does, it becomes the authoritative run and supersedes `20260831/pooled-CPUE-VAL-1-adopted`.
2. **P2, the gear-track cross-check** (about 0.5 h, can share the run). The two-track agreement is the strongest internal check the project has and it has not been re-measured since the shore component moved. Expect the gap to change, since the gear track fits this component at monthly.
3. **Optional, low priority: re-run C2** with the fixed helper so the three coarse rungs carry comparable adequacy. It would not change the decision, since the rule selects the finest adequate rung and weekly is adequate; it would only close the record.
4. **The hurdle likelihood** as the next modelling step, with the gain bounded in advance at roughly 6% of the catch so the effort can be judged before it is spent.

Standing and unchanged: the OSP crabbing-only column (a data request), more boat interviews next season (field protocol, and the acknowledged limitation of the whole OSP approach), and the stale weather module.
