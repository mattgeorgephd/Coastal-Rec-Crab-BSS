# Review: the adoption render, and a new authoritative run

**Date:** 2026-09-08
**Batch:** `06_diagnostics/run_adoption_2026-09-07.R`, results at `4ca59e6`
**Runtime:** 4.1 h (A1 220.4 min, A2 26.4 min)
**Verdict:** the gate passes. `05_output/20260904/pooled-CPUE-AD-A1-adopted` is the new authoritative run.

---

## 1. The gate

**A1's four fits are identical to stage C1's across 10,253 shared parameter rows, at full precision.**

That is the whole point of the run. C1 reached the weekly shore AR through `ar_force`, an experiment lever that takes the `fixed` branch of `bss_select_ar_resolution()` and skips the cap. A1 reaches it through `ar_max_resolution`, which lets the data-driven selector pick daily and then coarsens. Different code paths, same model, and now demonstrated rather than assumed: the desk stage proved the two build byte-identical Stan data before the render started, and the render proved the posteriors match.

So the adoption is a routing change and nothing else, which is what makes the folder callable authoritative.

## 2. The new authoritative run

**`05_output/20260904/pooled-CPUE-AD-A1-adopted`, port total 72,027 [53,018, 101,364]**, 4 of 4 components fitted, superseding `20260831/pooled-CPUE-VAL-1-adopted` (71,513). That is **+0.72%**, essentially all of it the shore effort process no longer being over-imputed at a daily AR.

C1 read 72,032 on the same fits. The 5-crab difference is the port assembly resampling component draws under `rstan::extract(permuted = TRUE)`, which permutes them: port totals move about 0.2% between bit-identical fits, which is why the gate is judged on per-fit summaries and never on this line.

| component | catch | method |
|---|---:|---|
| shore all-gear | 21,489 | BSS |
| shore pot closure | 6,320 | BSS |
| private boat all-gear | 31,008 | BSS |
| private boat pot closure | 1,018 | BSS |
| commercial / charter | 11,821 | census |

**Adequacy, and this is the first run where it is clean across the board:**

| fit | `p_loo_frac` | Pareto k>0.7 | `cov50` dev | miscalibrated |
|---|---:|---:|---:|---|
| shore all-gear | 0.0949 | **0** | 0.0659 | FALSE |
| shore pot closure | 0.1129 | 1 | 0.0769 | FALSE |
| private boat all-gear | 0.0801 | 1 | 0.0954 | FALSE |
| private boat pot closure | 0.2120 | 1 | 0.1136 | FALSE |

**No fit carries the miscalibration flag.** The shore all-gear component, which was the standing caveat on every headline since 2026-08-31 (p_loo 35.2% of `n_obs`, 41 Pareto k above 0.7, `coverage_50` 0.701 at +7.1 sampling SD, flag set), is now the best-behaved fit in the run.

## 3. The cross-check, and what it actually says

The gear-resolved track reports **71,026 [52,498, 100,038]**, so the two tracks are **-1.39%** apart at the port, against **-0.78%** for the last measured pair. **The gap widened**, which is the opposite of what I predicted.

The prediction was wrong for an interesting reason, and the number that explains it is this:

| comparison | pooled | gear | gap |
|---|---:|---:|---:|
| shore all-gear, as configured (pooled weekly vs gear monthly) | 21,489 | 20,754 | **-3.42%** |
| shore all-gear, **both at monthly** | 20,771 | 20,754 | **-0.08%** |

**At the same resolution the two tracks agree to 0.08%.** Two independently parameterized CPUE structures, one pooled and one carrying per-gear CPUE with `G = 5` on the shore, fitting the same data and landing within 17 crab of each other. That is the strongest agreement this cross-check has ever produced.

So the -1.39% port gap is a **resolution difference, not a disagreement about the fishery**. The pooled shore component moved to weekly; the gear track still fits it at monthly because `ar_max_resolution$gear_resolved` was deliberately left untouched so the cross-check would measure the pooled change alone.

`tau_bar` tells the same story: pooled 2.6128, gear 2.6114, **0.05% apart**, each estimated independently from a different CPUE structure. The shared-turnover adoption remains a property of the data rather than of the pooled parameterization.

**Two controls held.** The gear track is bit-identical to its own 2026-09-01 baseline on all four components, confirming that `estimate_catch_zi` is genuinely inert there. And that inertness is itself worth stating: `crab_bss_gear_resolved.stan` has no `theta_C` and `prep_bss_crab_gear.R` never emits `zi_catch`, so **the two tracks now differ in the shore catch likelihood as well as the AR resolution**, for the first time. That is worth about -0.3% on the pooled shore component, so it explains a sliver of the gap and none of the rest.

**The open question this raises.** Should the gear track's shore cap move to weekly too, so the cross-check is like for like? Probably, but not on the strength of this. Its shore fits carry per-gear CPUE, a thinner likelihood per gear type, so it may genuinely need a coarser AR than the pooled model does. That needs its own ladder on the gear track, not a cap copied across, and the ladder machinery now exists and records per-rung adequacy.

## 4. One defect, and the third instance of its shape

**`verdict_A2` crashed after 4.1 h of fitting and the run's verdicts were never written.** The pooled driver labels the port-total row `"Expected_Catch"`; the gear-resolved driver labels it `"Catch"`. Looking for the pooled label in the gear file returns `numeric(0)`, and `is.finite(numeric(0)) && ...` is an **error** in R rather than FALSE, so the block aborted before the write. The committed verdicts file held only the DESK row from my dry run.

Nothing was lost from the run itself: both folders are complete, and re-scoring from disk takes seconds. But this is the **third time** a defect in verdict code has aborted a batch after all the expensive work succeeded (`fmt()` on a vector, 2026-09-06; the port-row label, this run), and on two of those occasions the run's own verdicts were lost.

Fixed three ways: the port row is matched by pattern rather than by one track's label; every verdict block is now wrapped so a defect records itself as an `ERROR` row and the verdicts are written regardless; and the harness asserts both, plus asserts that the two tracks really do still use different labels, so the pattern match cannot be "simplified" away.

Harness: **452 assertions, 0 failing.**

## 5. What is now settled, and what is not

**Settled.** The shore all-gear AR resolution (weekly, adopted, routed through the cap and proven identical to the candidate). The zero-inflated shore catch likelihood (adopted). The shared turnover (adopted 2026-09-01, re-confirmed here at 0.05% cross-track). The PE/BSS incomplete-trip arm alignment (adopted 2026-09-02). The gear-track boat all-gear sampler settings (fixed 2026-09-02, and that fit still samples).

**Not settled, in order of how much it could move the number.**

1. **The boat is 43% of the port total and rests on a thin series.** 79.8% of the boat all-gear catch is extrapolated from 66 sampled days out of 289. No modelling change manufactures information that was not collected; this is the standing case for more boat interviews.
2. **The OSP crabbing-only column has still not arrived.** The machinery to use it as a hard lower bound on `f` is built, tested and inert. Without it, `f` sits on a set value of 0.30 for every stratum.
3. **The gear track's shore cap.** Now the largest single contributor to the cross-track gap, and answerable with a gear-track ladder.
4. **The catch stream is under-covered at every AR resolution** (-3.5 to -4.6 sampling SD). The ZINB improves it and does not fix it: the shore all-gear one bin is still 3.3 SD out because the data is more bimodal than a zero-inflated NB can be. A hurdle or two-component mixture is the shape that fits, with the gain bounded in advance at roughly 6% of the catch.
5. **`ppc_draws_*.rds` are written but nothing reads them back.** The recompute path is the last piece of a feature that has now paid for itself twice on paper and zero times in practice.

## 6. Next steps

There is no compelling reason to run anything immediately. The estimate is in a defensible state and the remaining items are either data requests or bounded improvements. In rough order of value per hour:

1. **Nothing, for a moment.** The authoritative run is current, its adequacy is clean, and the cross-check is understood. If a number is needed for the FW creel team or for co-management discussion, it is ready now, with the caveats in section 5 attached.
2. **The gear-track shore ladder** (about 3 h). `ar_escalate = list(shore = "all_gear")` on the gear track with `ar_rung_adequacy = TRUE`. It would settle the largest remaining cross-track difference and it exercises the ladder on the track it has never run on.
3. **The `ppc_draws` recompute path** (code, no run). Turns the next diagnostic defect from a multi-hour re-fit into a recomputation.
4. **The hurdle likelihood** (a Stan edit, a recompile, and one run). Bounded gain, and the bound is known before the effort is spent.
5. **Port the ZINB to the gear model** (a Stan edit and a recompile) if and when symmetry of the cross-check matters more than it does today.

Standing and unchanged: the OSP data request, boat interviews for the 2025-26 season, and the stale weather module.
