# Gear-Resolved BSS Pipeline: State Review and Improvement Plan

**Run reviewed:** `05_output/20260709/gear-type-CPUE-model/` (commit `564071b`, "7/9 gear resolved run")
**Comparison run:** commit `8101d65` (pre-F1/F2), and pooled `20260627`
**Scope of this review:** all 107 output files, plus code and raw-input verification

---

## 0. Executive summary

The three fixes shipped last session (F1 `R_G_boat`, F2 boat deployment scale, F5 reporting) all worked, and the F4 diagnostics paid for themselves immediately by exposing a defect I had not predicted. Two of three fitted components now report BSS instead of falling back to the point estimator.

Two new findings dominate everything else:

1. **The shore effort unit is invalid.** Catch scales as `h^0.57`, not `h^1`, against crabber-hours. This is the same failure mode as the boat's gear-hours, one step less extreme, and shore is the majority of the harvest.
2. **The Point Estimator is internally inconsistent, and the convergence gate falls back to it.** Its implied CPUE is 2.9x the ratio-of-sums computed on its own interviews. In this run the gate discarded a good shore BSS estimate in favour of a bad PE estimate.

**The port total of 88,238 should not be reported.** A defensible interim figure, using BSS wherever a fit exists, is approximately **82,900**.

---

## 1. What is now demonstrably fixed

### 1.1 F1: the `R_G_boat` improper prior (blocker)

| fit | before | after | gate |
|---|---|---|---|
| shore_ring_net | 7,818 / 8,000 (97.7%) | **598 / 8,000 (7.5%)** | FAIL (divergence fraction only) |
| shore_all_gear | 7,807 / 8,000 (97.6%) | **213 / 8,000 (2.7%)** | **PASS → BSS** |
| private_boat | 109 / 4,000 (2.7%) | **83 / 4,000 (2.1%)** | **PASS → BSS** |

R-hat is now 1.0002 to 1.0026 and `n_eff` 2,937 to 9,313 across all reported quantities. The diagnosis was exactly right: an unbounded `real<lower=0>` parameter with no prior and no likelihood, whose Jacobian `+z` drove the log posterior to `+inf`.

### 1.2 F2: boat effort on the gear-deployment scale

- Boat catch: **155,038 → 43,268**
- `R_G_boat` = 3.554 (3.244 to 3.874), well identified
- Boat CPUE = 2.930 crab per deployment
- `beta_h` (linearity) = **0.754, 95% CI 0.468 to 1.039**, which now covers 1.0

Cross-check against Grays Harbor interviews only: ratio-of-sums 2.555 crab per deployment; the BSS's own filtered subset gives 3.357. The model's 2.930 sits between them. A data-only expansion (14.4 trips/day x 3.51 gear/trip x 2.555 crab/deployment x 289 days) gives ~37,300, so the model is 1.16x a purely empirical expansion. That is a defensible position.

### 1.3 F5 and earlier work, all verified in outputs

- `ar_resolution` is populated (`biweekly`, `monthly`) in `fit_data_summary.csv` and `convergence_report.csv`
- `convergence_report.csv` now carries `divergence_fraction`, `impact_C_sd`, `impact_E_sd`, and the four `pass_*` flags
- `expansion_ratios.csv` correctly reports boat `R_G` as "prior-only (no gear counts in boat fits)", with CI 0.722 to 2.299 confirming it
- `port_total` carries `composition` and `effort_units_note`
- **B1.6 confirmed empirically:** boat `sigma_IE` median = 0.1403 against `exponential(5)` median = `ln2/5` = 0.1386. Prior-only sampling to three digits, proving `IE_n = 0` and that the shore-gate and the unconditional prior both work.

### 1.4 The `B1_C` question is answered

| fit | `B1_C` median | 95% CI |
|---|---|---|
| shore_all_gear | **-0.612** | -0.751 to -0.473 |
| shore_ring_net | **-0.335** | -0.512 to -0.154 |
| private_boat | -0.028 | -0.323 to 0.270 |

Weekend and holiday CPUE is about 46% lower for shore (`exp(-0.612) = 0.54`) and indistinguishable from zero for the boat. That is physically sensible: shore crowding depresses catch rates, soaking gear does not care what day it is. `B1_C` should stay on.

---

## 2. Correction to my earlier analysis

**My "boat catch should be ~90,700" claim was wrong by a factor of about 2.4, and I want that on the record.**

I computed that target from all boat interviews in the 2024-25 season. Those span three water bodies: Columbia River (1,679), Willapa Bay (965), Grays Harbor (607). The pipeline correctly restricts to Grays Harbor. Restricted to Grays Harbor private boats:

| quantity | all water bodies (my error) | Grays Harbor (correct) |
|---|---|---|
| crab per trip | 21.75 | **8.97** |
| gear per trip | 4.46 | **3.51** |
| crab per deployment | 4.881 | **2.555** |
| implied season catch | 90,735 | **~37,300** |

The model's 43,268 is therefore reasonable, not "2x low" as my earlier framing would have implied. The lesson is that any hand-built cross-check must apply the same population filters as the pipeline.

---

## 3. New finding A: the shore effort unit is invalid

`cpue_linearity_*.csv` tests the likelihood's own assumption. `c[a] ~ NB2(lambda_C * h[a], r_C)` asserts `log E[c] = log(lambda_C) + 1 * log(h)`.

| component | effort unit | `beta_h` | 95% CI | flagged |
|---|---|---|---|---|
| private_boat | gear-deployments | 0.754 | 0.468 to 1.039 | no |
| **shore_all_gear** | crabber-hours | **0.571** | 0.500 to 0.641 | **YES** |
| **shore_ring_net** | crabber-hours | **0.620** | 0.530 to 0.710 | **YES** |

`cpue_saturation_*.csv` (catch per gear against hours per gear) gives `beta` = 0.22 to 0.27 for all three components, including shore.

Shore's `h` spans 0.5 to 102 crabber-hours, a 200-fold range, with a median of 6. Forcing a 200-fold denominator range through a single proportional constant is exactly what inflated the boat.

The estimator triad confirms the consequence:

| component | ratio-of-sums | model implied | mean-of-ratios | model position |
|---|---|---|---|---|
| shore_all_gear | 0.3025 | 0.3657 | 0.3866 | **75% toward MoR** |
| shore_ring_net | 0.5018 | 0.6231 | 0.6560 | **79% toward MoR** |
| private_boat | 3.357 | 2.930 | 3.566 | below RoS (healthy) |

Both shore fits triggered `estimator_drift_flag`. With `r_C = 0.748` the NB2 behaves close to a multiplicative-error model and drags `lambda_C` toward mean-of-ratios.

**Consequence:** shore BSS catch is plausibly on the order of 20% high, and `lambda_C` for shore is not a stable parameter across seasons whose trip-length mix differs.

**Nuance worth preserving:** ring nets are actively pulled, so their catch genuinely should scale with time more than a soaking pot does. `beta_h` = 0.620 for ring net against 0.571 for all-gear is consistent with that. Neither is close to 1.

---

## 4. New finding B: the Point Estimator is internally inconsistent (highest priority)

The PE's implied CPUE (`PE_catch / PE_effort`) against the ratio-of-sums on its own interviews:

| component | PE catch | PE effort | implied CPUE | interview RoS | ratio |
|---|---|---|---|---|---|
| shore_all_gear | 32,937 | 37,825 | 0.871 | 0.3025 | **2.88x** |
| shore_ring_net | 12,174 | 9,281 | 1.312 | 0.5018 | **2.61x** |
| private_boat | 22,070 | 11,886 | 1.857 | 3.357 | **0.55x** |

The direction differs by component, so this is not a unit error.

**Cross-run evidence.** Between the two runs the shore day-length basis changed. PE effort fell 3.07x (116,220 → 37,825) while PE catch fell only 1.11x (36,578 → 32,937). Since `est_catch = est_total * mean_cpue`, catch must scale with effort at fixed CPUE. It did not. Catch and effort are not being computed on the same basis.

**Independent verification.** Computing directly from raw Grays Harbor shore interviews (Dec 1 to Sep 15, complete trips): ratio-of-sums 0.2747, mean of daily ratios 0.332, `n_int`-weighted 0.304, mean across 70 weekly x day-type strata 0.320. Nothing approaches 0.871. `fishing_time_total` is confirmed to be genuine crabber-hours (`crabber_hours = crabbers * hours_fished`, correlation 1.0000), so the denominator itself is correct.

**Most likely mechanism.** `mean_cpue` is `weighted.mean(daily_cpue, w = n_int)` per `(section, period, day_type)` stratum, where `daily_cpue` is a per-day ratio. With `period_pe = "week"` and roughly 50% day coverage, many strata rest on one or two sampled days, and days with very small sampled effort (minimum observed: 1.0 crabber-hour) produce extreme daily ratios that are then multiplied by the full stratum effort. This is a mean-of-ratios-of-sums estimator, unstable in both directions: heavy-tailed daily CPUE inflates shore, zero-catch days deflate the boat.

**Why this is the top priority.** The convergence gate routes failing components to the PE. In this run `shore_ring_net` failed the gate **only** on the divergence-fraction backstop (7.47% against a 5% threshold), while the scale-aware impact test said the divergences move the totals by 0.0044 SD (catch) and 0.0077 SD (effort), which is nothing. The gate therefore replaced a BSS estimate of **6,930** (pooled independently gives 6,733, agreement within 2.9%) with a PE estimate of **12,174**. The safety mechanism made the answer worse.

**Diagnostic that settles it in seconds:**

```r
# inside run_pe, after catch_strat is built
implied <- sum(catch_strat$est_catch) / sum(effort_strat$est_total)
ros     <- sum(daily_cpue$catch)     / sum(daily_cpue$hrs)
stopifnot(abs(implied / ros - 1) < 0.15)   # these must agree
```

**Fixes:**
1. Compute stratum CPUE as a ratio-of-sums within stratum (`sum(catch) / sum(hrs)`), not a weighted mean of daily ratios.
2. Coarsen `period_pe` for CPUE (month rather than week), or pool strata with fewer than N sampled days.
3. Add the assertion above as a permanent diagnostic.

---

## 5. Residual convergence issues

### 5.1 `sigma_IE` funnel (shore_ring_net, 7.5% divergences)

`divergence_localization_shore_ring_net_only.csv` puts `sigma_IE` at the largest absolute SMD, **-0.420** (bulk median 0.1167, divergent median 0.0495). With only 2 in-window I/E days, and `exponential(5)` placing its mode at zero, `sigma_IE` can shrink toward 0; as it does, the lognormal likelihood becomes razor-sharp on `lambda_E * L` and the geometry stiffens.

Options, in order of preference:
1. Tighten the prior to something like `lognormal(log(0.3), 0.5)`, which keeps mass away from zero without a hard bound.
2. Add `ie_min_obs_shore` (>= 3) before enabling the stream, mirroring `ie_min_obs_boat`.
3. A hard lower bound (`real<lower=0.05> sigma_IE`) works but is a hack.

Separately: `shore_all_gear` has `sigma_IE` = 0.603 (0.334 to 1.083). The 4 I/E observations and the effort counts disagree by roughly 60% on the log scale. That deserves its own investigation. Are the I/E days peak-count days, and is `lambda_E` on a monthly AR able to represent a day-specific total?

### 5.2 Non-identified hierarchical mean layer at `G = 1`

With `G = 1` and `S = 1`, `mu_E = mu_mu_E + sigma_mu_E * eps_mu_E` has a single element, so `mu_mu_E` and `sigma_mu_E * eps_mu_E` are jointly unidentified: only their sum enters the likelihood.

Evidence: `sigma_mu_E` median 1.13 to 1.20 with `hi95` around 7.5 (pure prior), `mu_mu_E` 95% CI spanning about 5 log units, and `sigma_mu_E` appearing in the divergence localization of every fit (SMD 0.163 shore ring, -0.147 boat).

**Fix:** when `G * S == 1`, drop the layer (`mu_E = mu_mu_E`, `sigma_mu_* = 0`). This removes a redundant funnel at zero inferential cost. It applies to `sigma_mu_C` identically.

### 5.3 Gate policy: the divergence backstop versus the impact test

`shore_ring_net` passes R-hat, passes `n_eff`, passes the impact test, and fails only `pass_div_fraction`. The scale-aware gate (B1.8) was built precisely to distinguish "divergences that matter" from "divergences that do not", and then a hard 5% count-based backstop overrides it.

I still favour keeping the backstop on principle: divergences indicate biased exploration whether or not the summary statistic moved (Betancourt 2017). **But that argument only holds if the fallback estimator is trustworthy, and right now it is not.** Fix the PE first, then revisit whether 5% is the right threshold.

---

## 6. Model fit quality (the good news)

- **PPC:** 95% coverage 0.951 to 0.985 across every stream; PIT means 0.495 to 0.523. Well calibrated in the tails.
- **PPC caveat:** 50% coverage is too wide for shore catch (0.737) and boat trailer (0.692), consistent with `r_C = 0.748` producing NB intervals wider than the data warrant.
- **LOO:** exactly one Pareto k above 0.7 (boat catch, of 131 observations). `p_loo` 5.8 to 12.3 against 52 to 1,651 observations. No pathologies.
- **Effort overdispersion:** NB observation overdispersion dominates (81% to 91%); the latent AR carries only 4% to 13%. Variance-to-mean ratios of 13 to 16. This is a structural consequence of a coarse (biweekly/monthly) AR that cannot track daily variation, so day-level variance lands in `r_E`. It is a design consequence, not a bug, but it is the reason the AR contributes so little.

---

## 7. Where the numbers actually stand

| component | PE | BSS | pooled BSS | gate selected | status |
|---|---|---|---|---|---|
| shore_ring_net | 12,174 | **6,930** | 6,733 | PE | PE is wrong; BSS matches pooled |
| shore_all_gear | 32,937 | **19,680** | 19,116 | BSS | BSS credible, possibly ~20% high |
| private_boat_ring | 971 | n/a | n/a | PE | insufficient data |
| private_boat_all | 22,070 | **43,268** | 74,025 | BSS | plausible; 1.16x a data-only expansion |
| comm/charter | 12,007 | census | census | census | fine |

Reported port total: **88,238** (95% CI 75,656 to 106,266). This uses the PE for `shore_ring_net`.

Reconstituted using BSS wherever a fit exists: **6,930 + 19,680 + 971 + 43,268 + 12,007 = 82,856**.

Neither should be published until Section 4 is resolved.

Two dependencies to state explicitly in any co-management conversation:
- Boat catch is proportional to `tau`, whose prior mean (1.2) rests on **two** WBL I/E days (turnover 1.00 and 1.29). Its uncertainty now propagates, but the evidence base is thin. `IE_n = 0` this season.
- Shore catch rests on an effort unit that the pipeline's own diagnostic flags as invalid.

---

## 8. Prioritized improvements

### P0. Blocking any reported number

1. **Fix `run_pe` stratum CPUE.** Ratio-of-sums within stratum; add the implied-CPUE assertion. (Section 4)
2. **Revisit the gate's PE fallback.** Do not fall back to an estimator that has not been validated against the interview ratio-of-sums.

### P1. Validity of the shore estimate (shore is the majority of catch)

3. **Resolve the shore effort unit.** Test a deployment-style denominator (`number_of_gear`, or gear-lifts where recorded) against crabber-hours, and adjudicate with PSIS-LOO on the catch stream. `log_lik_catch` already exists, so this is a two-run comparison.
4. **`sigma_IE` prior and minimum-observation guard**, and investigate the 60% I/E-versus-effort disagreement in `shore_all_gear`.

### P2. Model correctness, cheap to implement

5. **Drop the hierarchical mean layer when `G * S == 1`.** Removes a redundant funnel.
6. **`tau` identification.** The boat I/E path is built and inert. It activates automatically once two or more WBL days fall inside an estimation window. The egress-classification pilot is the critical path for the boat estimate.
7. **Pooled model.** `T_A_int = rep(1L, ...)` makes `T_A_int ~ bernoulli(R_T)` a likelihood of `R_T^n`, which with `beta(5,1)` pins `R_T` at 1.00 (0.98 to 1.00) and forces `lambda_E` to be groups while `h` remains gear-hours. Pooled's boat `E_sum` is group-hours mislabelled as gear-hours; its boat number is right only through compensating errors. Pooled also filters incomplete trips nowhere, while 39.8% of shore all-gear interviews are incomplete with a measured -20% CPUE bias.

### P3. Housekeeping and known gaps

8. `season_summary.csv` reports a PE-based mode breakdown alongside a gate-combined total. Inconsistent presentation.
9. `port_total` Effort row sums crabber-hours and gear-deployments. Flagged in the file, but it should probably not be summed at all.
10. **Option A** (genuine gear resolution) remains unbuilt. `G = 1`; the per-gear machinery in the Stan model is inert. Gear-type catch comes from PE apportionment, which given Section 4 is now doubly suspect.
11. `renv.lock` is absent. 37 MB of git-tracked outputs, some pre-v6 and unreproducible.

---

## 9. Recommended sequence

**Session 1 (P0).** Fix `run_pe`, add the assertion, decide gate policy. Re-run. Success criterion: PE implied CPUE agrees with interview ratio-of-sums to within about 15% for every component, and the PE-versus-BSS gap becomes a model-disagreement diagnostic rather than an artifact.

**Session 2 (P1).** Shore effort unit. Run the model twice, once with `h = crabber-hours` and once with a deployment denominator, and compare `elpd_loo` on the catch stream. Decide on evidence, not on symmetry with the boat.

**Session 3 (P2).** The `G*S == 1` mu layer, the `sigma_IE` prior, and then the pooled `T_A_int` fix with its own before-and-after comparison.

The single most important structural insight from this run generalizes beyond the boat: **for trap and pot gear, catch per unit soak time is not a stable parameter, and any effort unit denominated in time will be unstable.** The pipeline now measures this automatically (`cpue_linearity_*.csv`, `cpue_saturation_*.csv`), which means it can be checked every season rather than rediscovered.
