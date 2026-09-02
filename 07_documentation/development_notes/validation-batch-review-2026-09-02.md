# Validation batch review: an adoption confirmed, a production fit that never worked, and where the daily-AR argument points next

- **Batch:** `06_diagnostics/run_validation_2026-09-01.R`, four desk stages and five fitted, branch `OSP-boat-count-incorporation`.
- **Landed:** `ccf2c0e`. Code and config byte-identical to `28ef9c9`; the only difference is the `DRY_RUN <- FALSE` needed to run it.
- **Runtime:** 14.1 h of fitting (V1 257 min, V2 26 min, V3 31 min, V4 259 min, V5 271 min) plus seconds for the desk stages.
- **Verdicts:** 19 criteria, 17 resolving as intended, 1 FAIL that is a real production defect, 1 FAIL that is my criterion.
- **Harness:** 233 assertions, 0 failing.
- **Convention:** no em dashes.

---

## 1. The adoption is confirmed

**V1 reproduces S2 exactly: 11,223 shared parameter rows across 8 summaries, identical at full double precision, with all four components matching (20,898 / 6,331 / 31,008 / 1,018).** The configuration now shipping in `run_config.R` is the configuration the Stage 5 batch measured, not something that drifted alongside it. The port total reads 71,513 against S2's 71,521, an 8-crab gap that is the documented `rstan::extract(permuted = TRUE)` draw permutation and not a discrepancy.

The second reason V1 had to run also paid out: **the boat `prior_vs_posterior_*.csv` files exist again**, and the file-set difference between V1 and S2 is now exactly those two files. `tau_bar` contraction is 0.216, and `prior_influential` reads TRUE, which is the open item about that flag's definition rather than a finding about `tau_bar`: contraction measures variance reduction only, and this posterior sits about 3.5 prior SDs from the prior mean.

**`05_output/20260831/pooled-CPUE-VAL-1-adopted` is the authoritative run.** Port total **71,513** [52,346, 100,742].

---

## 2. The production defect: the gear track has never been able to fit its boat component

This is the most important thing in the batch and it was not what the batch was looking for.

| run | sampler for boat all-gear | divergences | `tau_bar` n_eff | R-hat | gate | boat all-gear reported | port |
|---|---|---:|---:|---:|---|---:|---:|
| 2026-08-29 gear | production defaults | 554 | 49 | 1.0813 | REJECTED | 11,176 (PE) | - |
| **S3** (2026-08-30) | **experiment override** | **2** | **19,292** | **0.9997** | **BSS** | **30,760** | **70,886** |
| V2 (2026-09-01) | production defaults | 554 | 49 | 1.0813 | REJECTED | 11,176 (PE) | 51,385 |
| V3 (2026-09-01) | production defaults | 554 | 49 | 1.0813 | REJECTED | 11,176 (PE) | 51,350 |

The Stage 5 batch found the fix and put it in an **experiment-only** escape hatch (`bss_sampler_override`), so production never received it. The gear driver had explicit per-fit sampler settings for shore all-gear and for both pot-closure fits, and **none at all** for boat all-gear, which therefore fell through to the track defaults (2,000/1,000 draws, `adapt_delta` 0.90, `max_treedepth` 10) while the pooled track fits the same component on the same data at 5,000/2,500, 0.99 and 13.

**The cost: the gear-resolved port total reads 51,385 instead of about 70,900, a 27% understatement of the cross-check that validates the pooled headline.** Anyone comparing the two tracks in the last month was comparing a working pooled fit against a rejected gear fit.

**Fixed permanently 2026-09-02**: the four settings now live in the gear driver's `params_model` with their own branch (`else if (!is_shore && is_allgear)`), so they select for that component and nothing inherits them by accident. `bss_sampler_override` stays NULL in production, where it belongs.

**What I got wrong, stated plainly.** I demonstrated this fix on 2026-08-30, wrote it up as vindicating a departure from the plan, and left production without it. A fix that lives only in an experiment flag is not a fix. The batch caught it only because V2 and V3 happened to exercise the gear track at production settings.

---

## 3. The other FAIL is mine, and it is the same mistake twice

V2's criterion "turning on per-gear CPUE leaves the boat untouched" FAILED. It is not a code defect:

- against **S3**: boat fits DIFFER (4,367 shared rows)
- against the **2026-08-29 gear run at the same sampler settings**: boat fits are **bit-identical** across 4,367 shared rows

S3 carried the experiment sampler override, so the reference differed in `gear_resolved_G` *and* in the sampler. `gear_resolved_G` reaches the shore fits and nothing else, exactly as GR-7 Phase 0 says.

**This is the second time in three weeks I have compared against a reference differing in two ways**, the first being the S3 shore criterion on 2026-08-30, and both times it produced a FAIL that looked like a code defect and took manual diagnosis to clear. Judgement is evidently not reliable here, so it is now a check: `fit_exactness()` takes `expect_delta`, parses both runs' `run_parameters.txt`, and reports any config key that differs but was not expected, inside the verdict itself.

---

## 4. GR-7 Phase 2 works, and has now been sampled for the first time

Coded 2026-07-21, parsing under stanc since, never run. V3 against V2, shore all-gear:

| gear | median off | median on | shift | interval width off | width on | ratio |
|---|---:|---:|---:|---:|---:|---:|
| Trap | 7,602 | 7,584 | -0.2% | 2,442 | 2,879 | 1.18x |
| Mixed | 5,503 | 5,501 | -0.0% | 1,656 | 1,868 | 1.13x |
| Pot | 3,823 | 3,815 | -0.2% | 1,324 | 1,512 | 1.14x |
| Ring Net | 2,120 | 2,106 | -0.7% | 947 | 1,057 | 1.12x |
| Snare | 1,824 | 1,821 | -0.2% | 927 | 1,139 | 1.23x |

**Exactly what it was designed to do:** every interval widens (1.12x to 1.23x), no median moves as much as 1%, and the per-gear medians still sum to the component total (20,872 against 20,827). The Phase-1 intervals were described as "slightly too narrow" because they treated the interview gear shares as known; propagating Dirichlet uncertainty widens them by 12-23% and leaves the central estimates alone.

V2 also reproduces the 2026-07-20 per-gear ORDERING (Trap > Mixed > Pot > Ring Net > Snare) despite everything that has changed since, which is the part that should survive because it is set by the interview CPUE ratios.

**Phase 2 is validated.** Whether `gear_resolved_G` becomes the gear-model default is now a decision rather than an unknown, and it should be taken together with the sampler fix above, since neither track's gear breakdown is trustworthy while the boat fit is rejected.

---

## 5. The boat pot-closure 2x2 completes, and is additive too

| boat pot closure | AR monthly | AR biweekly |
|---|---:|---:|
| turnover OFF | 849 | 735 |
| turnover ON | 1,018 | **939** |

Interaction **+35 crab** against main effects of +169 and -114. The same additivity the all-gear 2x2 showed. V4 also confirms the `ar_force` leak stays closed: the all-gear component came back at 31,008, S2's value exactly, where the 2026-08-28 leak had moved it 3,025 crab.

**V5 prices the floor.** Raising `shared_tau_min_obs` from 15 to 20 drops the boat pot closure (18 informed days) out of the feature, and it reverts to 849, the per-day value. The threshold is therefore worth **169 crab, 0.24% of the port**. That is the honest answer to a threshold set just below the component it keeps: the choice is defensible and it is also nearly inconsequential, which is the best possible outcome for a number chosen after seeing which side of the line a component fell on. The all-gear component (130 informed days) is unmoved at 31,008 either way, as it must be.

---

## 6. What the desk stages retired, with no MCMC

- **Tier 1, "TOP OF THE LIST", the shore I/E observation-unit fix: CLOSED on all four of its own pre-set criteria.** (a) `sigma_IE` 1.034 [0.64, 1.68] to 0.341 [0.02, 0.71], 0.272 today; (b) shore all-gear BSS 20,708 to 21,017, +1.5%; (c) the shore `L` posterior span widens 0.354 to 0.873 and reaches BELOW the 1.700 prior centre, where before it could only go up; (d) the boat bit-identical across 3,158 shared rows while the shore differs. The isolating run has existed since 2026-08-26 and had never been scored. **This closes the Tier-2 GR-9 `sigma_IE` tension with it**, which the backlog said must not be closed on reasoning alone.
- **Tier 2, the `gear_only` arm.** Adopting it shifts `R_G` by at most 0.7% (shore all-gear) and 0.1% (boat). The length-bias test is significant on shore all-gear (p = 0.0006) only because n = 2,741 makes a 0.7% shift detectable, and the shift is DOWNWARD, not the upward bias the rule was written against. The live defect is that the PE and the BSS use **different arms on both boat components** (BSS `exclude`, PE `gear_only`), which is cheap to fix and hard to defend leaving.
- **Tier 3, zero inflation: the precondition is now met, and it points at one stream.** Scored as observed zero count against the model's own expected count over ALL days, with a Poisson-binomial z:

  | fit / stream | n | observed zeros | expected | z |
  |---|---:|---:|---:|---:|
  | boat all-gear / trailer | 195 | 39 | 43.7 | -0.9 |
  | boat all-gear / osp | 130 | 7 | 4.1 | +1.5 |
  | boat all-gear / catch | 131 | 12 | 6.4 | +2.3 |
  | **shore all-gear / catch** | **1,649** | **676** | **605.4** | **+3.8** |
  | shore pot closure / catch | 627 | 146 | 118.5 | +2.9 |

  Every boat stream is fine. The two **shore catch** streams under-predict zeros in the same direction at z = 3.8 and 2.9, about 71 and 27 extra zero-catch interviews. That is a real but targeted misfit: a ZINB prototype, if built, should target the shore catch likelihood, not the whole model, and its available gain is bounded by roughly 70 observations out of 1,649.
- **Tier 4, hygiene.** `.Rproj.user` still tracked (5 files); `05_output` at 6,837 files and 565 MB, including 222 files from the stale 20260711 morning/afternoon runs; LICENSE 35,823 bytes with no placeholder, so the GPL-3.0 paste is done and only WDFW's confirmation of the copyright lines remains.

---

## 7. The finding the batch was not looking for: production's SHORE all-gear fit is the badly-behaved one

The adequacy table on the production run (V1) reads:

| fit | p_loo / n_obs | Pareto k > 0.7 | coverage_50 dev | smallest sigma_r | miscalibrated |
|---|---:|---:|---:|---:|---|
| shore pot closure | 11.4% | 1 | 0.077 | 0.506 | no |
| **shore all-gear** | **35.2%** | **41** | **0.201** | **0.234** | **YES** |
| boat pot closure | 21.2% | 1 | 0.136 | 0.667 | no |
| boat all-gear | 8.0% | 1 | 0.095 | 0.767 | no |

**The component the Stage 5 batch cleared is fine. The one nobody tested is not.** Shore all-gear runs at a DAILY AR in production, by `ar_max_resolution$pooled$shore$all_gear = "daily"`, and it shows the same signature the daily boat AR was rejected for: 109.4 effective parameters on 311 gear observations, 41 Pareto k above 0.7, and `coverage_50` at 0.701 against a nominal 0.500, which is **+7.1 sampling SDs**. That last number needs no cross-track comparison to interpret.

The gear track fits the same shore component at MONTHLY, and it is not better, it is differently wrong:

| shore all-gear | gear stream p_loo | Pareto k > 0.7 | gear stream coverage_50 |
|---|---:|---:|---:|
| pooled, DAILY (production) | 109.4 / 311 | 41 | 0.701 (+7.1 SD) |
| gear-resolved, MONTHLY | 11.3 / 311 | 0 | **0.035 (-16.4 SD)** |

Daily over-fits: it buys 101 nats of gear elpd for 98 effective parameters, about 1.0 nats each, which is the rate a pure-noise parameter buys, and it costs 5 nats on the catch stream. Monthly cannot track the series at all: only 3.5% of observations fall inside a nominal 50% interval.

**Caveat, stated because it limits the claim.** That table is a CROSS-TRACK comparison, so three things differ at once: the AR resolution, the CPUE structure, and `gear_resolved_G`. It is suggestive of where the answer lies, not a measurement of it. The pooled daily number alone is sound, and it is enough to say production's shore all-gear component is miscalibrated.

**Neither end is right and the answer is between them.** Shore all-gear is 20,898 of a 71,513 port total, 29%. Nobody has run it at weekly or biweekly. Plan 5.4 explicitly set the shore cap aside as untested; it now has evidence and deserves the same treatment the boat got.

---

## 8. Two diagnostic fixes that rode along

- **The two PPC files now agree exactly, not approximately.** After the 2026-09-01 randomized-coverage fix they still differed by 0.005 to 0.015, because the aggregate file estimated the PIT by simulation while the per-observation file computed the exact expectation. `calib()` now uses the same closed form. Any future disagreement between the two is a defect rather than noise to be eyeballed.
- **D3 reads the newest production run.** It was pinned to archived folders, none of which carry `p_zero`, so its zero-bin criterion was permanently unscorable. It now picks up `VAL-1-adopted` and scores it.

---

## 9. Where the numbers stand

| run | configuration | port | boat all-gear | boat pot closure |
|---|---|---:|---:|---:|
| rung 4 | pre-adoption shipped | 66,237 | 25,868 | 849 |
| **V1** | **production, adopted** | **71,513** | **31,008** | **1,018** |
| V4 | + boat pot closure biweekly | 71,366 | 31,008 | 939 |
| V5 | + floor 20 | 71,205 | 31,008 | 849 |
| S3 | gear-resolved cross-check | 70,886 | 30,760 | 956 |
| V2/V3 | gear-resolved, boat fit REJECTED | 51,385 | 11,176 (PE) | 956 |

The pooled and gear-resolved tracks agree to 0.89% once the gear boat fit can sample. With the 2026-09-02 driver fix, that is now the production path rather than an experiment.
