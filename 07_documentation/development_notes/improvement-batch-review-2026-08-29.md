# Improvement-batch review, 2026-08-29

- **Reviewed:** `06_diagnostics/run_improvement_plan_2026-08-27.R`, six stages, ~22 h, results at `b6ec3b6`.
- **Baseline for every comparison:** 2026-08-26 ladder rung 4 (`05_output/20260826/pooled-CPUE-PV4-minint`), port total 66,237.
- **Companions:** `PIPELINE_STATUS.md` (state), `improvement-plan-2026-08-27.md` Stage 5 (what to do next).
- **Convention:** no em dashes.

## Bottom line

Every stage answered its question. Together they exposed a bigger one.

**The private-boat all-gear estimate ranges over 25,883 to 37,392 across six configurations that all pass
the convergence gate.** That is a 44% spread on the component and 18% at the port total. Every one of those
fits has R-hat below 1.01, n_eff above 400, a divergence fraction below 5% and divergence impact below
0.10 SD. For scale, the entire 2026-08-25 improvement batch moved the port total 1.6%.

The gate is not malfunctioning. It answers "did this fit sample?" and the project has been reading it as
"is this model right?". Those are different questions and only one of them is currently asked.

| Configuration | AR | boat all-gear | 95% CI | width | `r_OSP` | port |
|---|---|---:|---|---:|---:|---:|
| **shipped (rung 4)** | monthly | **25,883** | 10,400 - 50,675 | 40,275 | 1.62 | **66,237** |
| D, `shared_tau` off (gate) | monthly | 25,883 | 10,400 - 50,675 | 40,275 | 1.62 | 66,094 |
| B, rung-2 config reseeded | monthly | 27,560 | | 42,921 | 1.58 | 67,106 |
| C, boat AR forced biweekly | biweekly | 28,902 | | 45,569 | 1.68 | 68,853 |
| E, `shared_tau` on | monthly | 31,012 | | 47,604 | 1.54 | 75,619 |
| F, `ar_escalate` (settles daily) | daily | **37,392** | 15,592 - 70,368 | 54,777 | **10.61** | **78,250** |

## Stage by stage

### A. Empty-stratum fill (plan 1.4 / 4.3). PE only, minutes.

`day_type` against `zero`:

| component | zeroed days | PE catch `zero` | PE catch `day_type` | change |
|---|---:|---:|---:|---:|
| boat all-gear | 44/289 | 11,059 | 13,404 | **+21.2%** |
| boat pot closure | 8/76 | 403 | 511 | **+26.8%** |
| shore all-gear | 41/289 | 16,027 | 19,339 | **+20.7%** |
| shore pot closure | 0/76 | 5,891 | 5,891 | 0.0% |

About +13% at the port level. The `zero` fill has been biasing the PE down by roughly an eighth. It matters
in three places: the PE is the project's main sanity check on the BSS; the boat pot closure entered the
port total as a PE point before the interview floor changed; and any future PE fallback inherits it
directly. The unsampled strata skew weekend and holiday, and those days carry 1.7-2.3x weekday effort, so
`zero` is biased down by construction.

**Caveat, a defect in the batch and not in the result.** Stage A did not load the crabbing-holiday workbook
the way the driver does, so it stratified 289 days into 85 week x day-type cells against the driver's 92,
and reported 41/44/8/0 zeroed days against the driver's 44/47/9/0. Both arms shared the stratification, so
the 20-27% comparison is sound; the absolute counts are not production numbers. Fixed.

### B. Reseeded rung 2 (plan 1.2). The verdict says PASS and undersells it.

| run | config | divergences | fraction | worst chain | `E_sum` n_eff | verdict |
|---|---|---:|---:|---:|---:|---|
| 2026-08-26 rung 2 | I/E fix, Fri/Sat/Sun | 2,957 | 29.6% | **2,500 of 2,500** | 30 | REJECTED |
| **B (new seed)** | I/E fix, Fri/Sat/Sun | **488** | **4.88%** | **350** | **924** | passed, barely |
| shipped (rung 4) | I/E fix + Sat/Sun | 333 | 3.33% | 152 | 5,914 | comfortable |

Two seeds, two struggling chains: one fatal, one survivable at 4.88% against a 5% backstop. So the 2026-08-26
failure was seed-dependent, but "bad luck, move on" is the wrong conclusion. **The I/E unit fix on its own
strains the shore all-gear daily fit, and the weekend fix is what makes the shipped configuration
comfortable** - exactly what its rationale predicted, since it cuts the residual variance the AR must
absorb by about 15%.

The isolation the ladder lacked now exists, and it confirms the I/E result three ways. Shore all-gear
`sigma_IE`: **0.238** (B, the clean isolation), 0.300 (rung 3), 0.272 (rung 4), against 1.068 before the fix.

### C. Boat pot-closure AR (plan 3.1). Answered, on a contaminated run.

Forced to the gear track's biweekly, the pooled boat pot closure gives **735** against the gear track's
**743**, 1.1% apart, where the two tracks had disagreed 13% at monthly (849) versus biweekly (743).
**The gap was a configuration artefact.** Align the two `ar_max_resolution` maps.

The run was contaminated by a defect this review found. `ar_force` was per-POPULATION only while
`ar_max_resolution` had already gained the nested per-sub-season form, and `as.character()` on a
one-element list returns its element, so `ar_force = list(private_boat = list(pot_closure = "biweekly"))`
silently forced BOTH boat sub-seasons. Boat all-gear moved 25,883 -> 28,902 as a side effect, so stage C's
port total is not interpretable as the change it was meant to isolate. Fixed: `ar_force` now takes the same
shapes as `ar_max_resolution` and errors on an unnamed list.

### D. The gate (plan 2.1's control). PASSED.

`shared_tau = FALSE` reproduces rung 4 with **11,217 shared parameter rows identical at full double
precision across 8 summaries** and a byte-identical `convergence_report.csv`. The shared-turnover edit is
behaviour-neutral when off, so stage E is attributable to the feature and nothing else.

### E. Shared turnover (plan 2.1). The boat result is good. The shore result is not.

**Boat: `tau_bar` = 2.597 [2.064, 3.249]** against a prior centred on 1.20, excluding the prior centre at
about 2.6 prior SD. Four independent lines agree:

1. **Independent calibration.** `osp_trailer_overlap_calibration.csv`, unchanged across every run, implies
   a turnover of 2.01 / 2.74 / 3.03 by three trailer metrics (n = 61, correlation 0.90-0.98). The posterior
   sits inside that range.
2. **Both boat streams reconcile at once, in opposite directions.** Trailer PIT mean 0.4238 -> 0.4858 (bias
   down 81%), OSP PIT mean 0.5790 -> 0.5211 (bias down 73%). The trailer stream was under-predicted and OSP
   over-predicted by nearly equal amounts; one shared level fixes both. That is hard to do by accident.
3. **LOO improves on the stream it should.** Boat trailer `elpd_loo` -573.8 -> -559.1, a paired gain of
   +14.67 at z = +2.66, `p_loo` 4.9 -> 8.0, zero bad Pareto k. Every other stream is a wash.
4. **Convergence improved:** boat all-gear divergences 97 -> 39, `tau_bar` n_eff 16,750, R-hat 0.9997.

Two qualifications. The posterior lies almost entirely **above** the prior's 95% envelope (prior hi95 2.160
< posterior lo95 2.064), so the prior is pulling `tau_bar` down and 2.597 is shrunk, not prior-driven. And
this fit has `n_ie_obs = 0`: there are no direct turnover observations, so `tau_bar` is identified
indirectly through the trailer and OSP count streams. Legitimate inference; never call it measured.

**Shore: do not carry it forward.** Shore all-gear moved 20,898 -> 24,629 (+17.9%), not boat spillover
(separate fits, no shared parameters) but its own `tau_bar` moving 1.700 -> 2.080. Backed by:
only **4 of 289 shore days** carry any I/E observation; the interval **[1.358, 2.979] still contains the
1.7 prior centre**; it bought no measurable improvement (gear PIT 0.5005 -> 0.5012, catch elpd -0.18); and
**the gear track put the same parameter at 1.681**, no move at all. A 3,731-crab increase from a parameter
indistinguishable from its prior, improving nothing, not replicating. The toggle is global; it needs an
informed-day floor.

**The gear track's boat fit failed** (554 divergences, R-hat 1.32, n_eff 15) and gives no verdict. The
failure is diagnosable and is not about `tau_bar`: 544 of 554 divergences are in chain 3 alone (chains
1/2/4: 0/2/8), the divergence localization implicates `phi_E` at the AR unit-root boundary (SMD +1.26),
`tau_bar` does not appear in that table, and the gear track runs this fit at 1,000 iterations against the
pooled track's 2,500. Its boat POT-CLOSURE fit did pass, at `tau_bar` 2.050 [1.415, 2.930], also excluding
1.2 - weaker but real cross-track support.

### F. The escalation ladder (plan 3.2). The answer is uncomfortable.

**Every component passes at DAILY on the first attempt.** `ar_escalation_log.csv` records `attempt = 1`
for all four fits. Shore pot closure at daily gives 111 divergences (1.85%) where the Run 6 note recorded
1,165 (19.4%): **the funnel that justified the biweekly cap is gone**, presumably retired by the
gear-deployment unit, the weekend fix and the I/E unit fix between them. So `ar_max_resolution` is now the
only thing holding production at monthly and biweekly, and it is worth 44% on the boat.

But the daily boat fit is absorbing, not identifying:

- Boat trailer **`p_loo` 4.9 -> 31.8** with 2 bad Pareto k, for a +33 elpd gain. About 27 extra effective
  parameters on 195 observations.
- **The catch stream gets worse**: elpd -451.4 -> -455.3. The stream the harvest estimate is about does not
  improve.
- **`r_OSP` collapses 1.62 -> 10.61**, 95% CI [3.57, **1,405**], and its parent `sigma_r_OSP` drops to
  n_eff 307 with R-hat 1.014, below the gate's own n_eff floor - on a parameter the gate does not check.
  A 289-state latent process has swallowed the OSP observation error.
- The whole increase lands in the extrapolated portion: boat extrapolated catch 19,257 -> 26,895 on the
  same 66 of 289 sampled days.

Treedepth was checked and is not the issue (`treedepth_pct = 0` against `bss_treedepth_boat_allgear = 13`).

## Numbers that must not be quoted

- **`pe_vs_bss_comparison.csv` in the gear stage E run shows `BSS_catch = 32,689` for private_boat (All
  gear) while `method_selected = "PE"`.** That is the rejected fit. The gate did work downstream, and
  `port_total_Dungeness_Kept.csv` substitutes the PE, but the rejected value sits unmarked in the
  comparison CSV. Same for that fit's `tau_bar` 2.552 (n_eff 49), `f_crab` 0.315 (n_eff 25, R-hat 1.17).
- **`kappa_OSP` is `decoupled = TRUE` in all four fits of every run.** Its ~3.0 median is its prior. Under
  the production `osp_scale_is_tau = TRUE` the OSP mean uses `L`, not `kappa_OSP`.
- **Shore `r_OSP` posterior mean is 2.5 million** in one fit. Flagged decoupled. The median is meaningless too.
- **`R_G_boat` ~ 4.0 in the shore fits is decoupled.** `expansion_ratios.csv` still prints it without the flag.
- **The 2026-08-26 baseline `structural_params_*` files have no `decoupled` column at all**, so every trap
  above is unmarked there.

## What this changes

Nothing in the shipped configuration, deliberately. Rung 4 stands as the authoritative run. Stage E is the
strongest candidate to change it and should not be adopted on one pooled run with a failed gear
cross-check and a shore side effect that has to be suppressed first.

Next steps are Stage 5 of `improvement-plan-2026-08-27.md`. In order: adopt a fill for the empty strata
(minutes); make `shared_tau` boat-only and re-run both tracks (~4.5 h); run the 2x2 of shared tau against
AR resolution, which is one new cell and settles what the boat estimate actually rests on (~4 h); then
decide what `ar_max_resolution` is for, and give the gate something to say about model adequacy rather
than only about sampling.
