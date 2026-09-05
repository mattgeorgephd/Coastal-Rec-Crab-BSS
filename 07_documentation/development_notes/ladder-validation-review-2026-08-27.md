# Crab BSS ladder review, 2026-08-27

> **SUPERSEDED 2026-09-01.** The authoritative run this review crowns has been superseded; the current one is in the box at the top of `PIPELINE_STATUS.md` (`20260904/pooled-CPUE-AD-A1-adopted`, 72,027, as of 2026-09-08). This document is the record of what was true when written.

Branch `OSP-boat-count-incorporation` @ `ec9a52d`. Six model-runs, ~15 h wall clock.
References: `20260804/pooled-CPUE-boat-count-validation-run`, `20260805/gear-type-CPUE-model-boat-count-validation-run`.

## Bottom line

The ladder ran clean and the batch is sound. Ten of twelve criteria pass; the one FAIL and the one REVIEW
are both mis-specified criteria, not code defects. The shipped config moves the port total DOWN 1.1-1.6%
against the pre-patch production runs, and essentially all of it is the boat all-gear component responding
to the weekend redefinition, partly offset by the newly fitted boat pot closure. Three real defects
surfaced, all mine, all small. One pre-existing issue is larger than anything in this batch:
`osp_scale_is_tau` does not do what the documentation says it does.

## The ladder

| Rung | Track | Change | Port BSS | 95% CI | Effort | Fits | Status |
|---|---|---|---:|---:|---:|---:|---|
| ref | pooled | pre-patch production | 67,312 | 50,601-93,461 | 39,145 | 3/3 | - |
| 1 | pooled | pre-patch equivalent | 67,323 | 50,432-93,793 | 39,156 | 3/3 | reproduces |
| 2 | pooled | + shore I/E unit fix | 62,069 | 45,483-88,517 | 33,079 | 2/3 | one chain stalled |
| 3 | pooled | + weekend = Sat/Sun | 65,493 | 49,529-90,501 | 40,016 | 3/3 | as predicted |
| **4** | **pooled** | **+ bss_min_interviews 15 (SHIPPED)** | **66,237** | **50,037-91,210** | **40,187** | **4/4** | **clean win** |
| ref | gear | pre-patch production | 66,461 | 50,259-91,309 | 38,649 | 3/3 | - |
| 1 | gear | pre-patch equivalent | 66,442 | 50,337-91,740 | 38,676 | 3/3 | bit-identical fits |
| **5** | **gear** | **shipped config** | **65,756** | **49,660-91,476** | **40,178** | **4/4** | **cross-check holds** |

Components (Dungeness kept):

| Component | ref pooled | R1 | R2 | R3 | R4 shipped | R5 gear |
|---|---:|---:|---:|---:|---:|---:|
| Shore all-gear | 20,665 | 20,708 | 21,017 (rejected) | 20,898 | 20,898 | 20,754 |
| Shore pot closure | 6,263 | 6,255 | 6,255 | 6,331 | 6,331 | 6,332 |
| Boat all-gear | 27,668 | 27,750 | 27,750 | 25,868 | 25,868 | 25,796 |
| Boat pot closure | 351 PE | 351 PE | 351 PE | 400 PE | 850 BSS | 743 BSS |
| Commercial/charter | 11,986 | 11,986 | 11,986 | 11,821 | 11,821 | 11,821 |

Decomposition of the -1,075 pooled move: boat all-gear -1,882 (weekend), boat pot closure +744
(interview floor), shore roughly +190 net, commercial census -165.

## The two non-PASS verdicts

### FAIL, rung 1 gear bit-identity: false alarm

The gear model IS bit-identical. All 1,508 / 4,935 / 3,431 parameter rows shared with the reference,
across the three fits, match at full double precision. `convergence_report.csv` is byte-identical (22/19/6
divergences, n_eff 8,478/8,316/5,328). The only differences are three new inert reporters (`f_lower[1]`,
`f_lower_out[1]`, `osp_f_kappa_out`, all exactly 0 with 0 SD) and a pure rename `f_crab_param[1]` ->
`f_theta[1]`, identical to every digit.

The port total moves by 19 (~0.1 MCSE) because it is assembled from `rstan::extract(permuted = TRUE)`
draws and an unseeded `sample.int()` at `save_run_diagnostics.R:145`; only 490 of 2,000 draw indices are
common. PIPELINE_STATUS already recorded that behaviour in the Run 6 note. My criterion tested the wrong
object. Fix: gate exactness on `bss_summary_*.csv` + `convergence_report.csv` byte equality.

### REVIEW, rung 2 tau narrowing: retract the prediction

`crab_bss_pooled.stan:370,418` declares `vector[D * estimate_L] L_raw` with
`L[d] = L_data[d] * exp(L_prior_sigma[d] * L_raw[d])`, `L_raw ~ std_normal()`. There is NO shared
`tau_shore` parameter; L is 289 independent per-day draws each anchored on 1.7 with log-SD 0.3. Four I/E
days can inform four days out of 289. The reported 0.93 ratio is a median over 285 prior-only days; on the
four informed days the converged run gives 0.74, 0.69, 1.22, 1.15. Season median L = 1.6998 vs prior 1.7000.

The 2026-08-25 claim that the shore I/E stream "becomes genuinely informative about tau_shore" is wrong as
the model is written and should be retracted. Making it true is a model change, not a config change.

## Rung by rung

**Rung 1.** Pooled 67,323 vs 67,312 (delta 11, ~0.1 MCSE), components within 0.2%. PE side byte-identical
to the reference. Gear bit-identical. The baseline reproduces.

**Rung 2, shore I/E unit fix.** Effect confirmed and large: shore all-gear sigma_IE 1.068 -> 0.300 pooled
(-72%), 0.871 -> 0.598 gear (-31%). Arithmetic corroboration: geometric-mean hours:trips ratio on the four
in-window I/E days is 3.44; log(3.44) = 1.24 against a pre-fix sigma_IE of 1.03-1.07. The mismatch was
~3.2x, not the ~4x the notes estimated.

But the isolating fit failed: 2,957 divergences (29.6%), E_sum n_eff 30, R-hat 1.049, rejected to
"PE (convergence fail)". That substitution is the entire rung-2 port drop. Per-chain: chains 1/2/4 at
163/140/154 divergences (in line with rung 1's 385 across four), chain 3 at **2,500 of 2,500**, mean accept
0.053, mean treedepth 1.21 - it never moved. One stalled chain, not a global pathology.

RECOMMENDED: re-run rung 2 with a different `bss_seed` (~3 h). Also: the summary CSV reports rung 2's
shore all-gear BSS numbers with no indication the fit was rejected; the extractor should carry
`method_selected` per component.

**Rung 3, weekend = Sat/Sun.** B1 rises as predicted, slightly past the forecast: shore all-gear
1.67x -> 2.48x (B1 0.511 -> 0.907), shore pot closure 1.61x -> 2.02x, boat all-gear 1.56x -> 1.90x. Shore
pot-closure r_E rises 3.22 -> 3.90 and predictive VMR falls 14.9 -> 12.5, i.e. less residual for the AR.

This is the rung that moves the number: boat all-gear -6.8% (27,750 -> 25,868), unpredicted and the largest
single move in the ladder. Monthly pattern is summer down / shoulder up (Jul -370, Aug -496, Jun -180 vs
Sep +19, Oct +40, Nov +55, Apr +128).

Side effects: (a) the empty-stratum prediction was exactly right in count (boat pot closure 4/76 -> 9/76)
and backwards in sign - its PE moved UP 14% (351 -> 400) while the recorded rationale said "expect a
downward move"; the other two components moved down (48->47, 47->44 days). (b) The commercial/charter
"census" is not invariant - it moves 11,986 -> 11,821 because `estimate_comm_charter.R:93,107` stratifies
on `days_wkend`. (c) The BH opener screen is stable: one covariate survives in every run (boat trailers x
MA2 halibut, q = 1.1e-04); razor dig on shore effort does not (q = 0.12-0.18), reproducing the documented
worked example. K_open = 0 everywhere.

**Rung 4, interview floor 20 -> 15.** The boat pot closure reaches the sampler and PASSES the gate cleanly:
35 divergences (0.58%), R-hat 1.0003, n_eff 5,207, max treedepth 8. BSS median 850 (95% CI 274-2,693)
against a PE of 400. Isolation is perfect: rungs 3 and 4 are numerically identical on the three shared fits
(same per-chain divergence counts 110/333/97, same elpd_loo to five decimals).

Two cautions: it is the most extrapolated stratum in the set (13 of 76 days sampled, 76.2% of catch
extrapolated; B2, B1_C, B2_C all at their priors), and the two tracks disagree by 13% because they fit it
at different AR resolutions (pooled monthly P=3 -> 850; gear biweekly P=6 -> 743).

**Rung 5, cross-track.** Pooled 66,237 vs gear 65,756 (0.73%, threshold 2%); boat all-gear 25,868 vs
25,796 (0.28%, threshold 3%). Both tracks 4/4. Slightly tighter than the 1-2% before the batch.

## Defects found (all from the 2026-08-25 batch)

- **D1** `save_run_diagnostics.R:437` - a second, unfixed copy of the shore day-length weighting
  (`w <- if (is_boat) mc else mc * dl`). Defect (B) was fixed in `pe_monthly_effort_share.R` only. This copy
  feeds `monthly_pe_vs_bss.csv`, whose shore PE effort is numerically unchanged from the pre-fix reference.
  The fixed copy also never executed in this ladder (no shore PE fallback in any run).
- **D2** `run_pe_pooled.R:144` + `BSS-GH-pooled-CPUE-model.Rmd:903` - the empty-stratum count is computed
  and only `cat()`ed, and the pooled PE chunk is `results='hide'`. The claim that "every run now reports the
  zeroed-stratum count" is true only of the gear track; the 4/76 -> 9/76 numbers above came from the gear HTML.
- **D3** `run_patch_validation_2026-08-25.R:306` - `empty_effort_note` declared, never assigned; NA on all rows.
- **D4** `bss_day_length.R` - `ie_analysis.csv` and `L_effective_ie_detail.csv` are byte-identical across all
  eight runs and still report only `ie_crabber_hours`. You cannot audit from the outputs which observation
  column the model consumed, i.e. the exact thing rung 2 changed.
- **D5** `diagnose_incomplete_trips.R` - `production_arm` is a static "exclude" label, but the shipped boat PE
  (3,565.75 / 10,940.36) equals the `gear_only` arm. The diagnostic is right; the label should be per component.

Not broken after all: the `opener_covariates_*.csv` writer exists in both drivers and correctly wrote
nothing at K_open = 0. The `stan_console_*.log` files are absent by design (at cores > 1 a successful fit
reports nothing in the parent process).

## Pre-existing issues surfaced

**HIGH - `osp_scale_is_tau` pins the OSP scale at a prior the data contradicts.** With the toggle on, the
OSP mean uses `L[day]`, which is the same per-day independent construction as above: 289 draws anchored on
1.2 with log-SD 0.3, no shared parameter. Posterior median of L across 289 days is 1.201 against a prior
centre of 1.200; only 26 days narrow below 90% of prior width. 148 days of OSP data barely move it.

Meanwhile the data say the turnover is 2.0-3.0 (`osp_trailer_overlap_calibration.csv`: implied turnovers
3.03 / 2.74 / 2.01; Phase-1 free kappa_OSP posterior 3.15). That ~2.5x conflict is absorbed by the OSP
overdispersion instead: r_OSP lands at 1.58-1.62, loose enough to treat a systematic scale error as noise.
Corroborating: the boat trailer PIT mean is 0.420-0.425 in all six runs against a nominal 0.50, while every
shore and catch stream sits at 0.498-0.516 - the model predicts more trailers than observed, the signature
of lambda_E pulled up by OSP against a disagreeing trailer stream.

Reporting hazard: in this configuration kappa_OSP is decoupled and reports its prior exactly (median 3.008,
CI 1.63-5.40), sitting in `structural_params_*.csv` with nothing marking it inert. r_OSP in shore fits is
worse - 95% CI spanning roughly 0.0002 to 700.

This toggle roughly doubles the private-boat harvest, a quarter of the port total, and as implemented the
size of that doubling is set by the tau_boat prior rather than by OSP. Suggested: refit with a shared
`tau_boat` parameter carrying per-day deviations, plus a posterior-predictive check on the OSP stream.

**MEDIUM - the pooled shore all-gear daily-AR fit is near-saturated.** p_loo = 109 effective parameters on
311 gear observations, 41 Pareto-k > 0.7 (max 1.16). The gear track's monthly-AR fit of the same data:
p_loo 11.3, zero high-k. 143 of 289 daily AR periods carry no observation. Its elpd_loo is better
(-1,251 vs -1,353) but at k > 1 that margin is not a trustworthy model comparison. This is also the fit
that stalled a chain in rung 2. `ar_escalate` (improvement 7) was built for exactly this and shipped FALSE.

**LOW - documented extrapolation figures do not match any run.** Documented 79.8 / 48.5 / 30.9%. Observed:
boat all-gear 81.8-82.0%, shore pot closure 31.8-33.3%, shore all-gear 48.5% only at rung 1, rising to
51.0% shipped pooled and 54.2% gear. Boat pot closure at 76.2% is undocumented and the most extrapolated
stratum in the set.

Smaller: ring-net incomplete-trip bias flips sign with the unit change (-21.3% on deployments vs +4.4% on
legacy crabber-hours); nominal-50% interval coverage runs 0.55-0.84 for every fit (conservative
throughout); the estimator triad, linearity and saturation diagnostics are bit-identical across all six
runs and discriminate nothing between rungs; no cpue_linearity/cpue_saturation file is written for the boat
pot closure even though the fit now exists.

## Status of these findings

All six defects above are FIXED as of 2026-08-27, in the same patch that added this file. The
empty-stratum fix (D2) immediately produced a finding that had been invisible on the pooled
track: on the shipped config the zeroed-day shares are shore all-gear 44/289 (15.2%), boat pot
closure 9/76 (11.8%), boat all-gear 47/289 (16.3%), shore pot closure 0/76. Three of four
components sit above the 5% warning threshold at the default `zero` fill, and the
previously-quoted "4 of 76 rising to 9 of 76" was the boat pot closure alone.

The two mis-specified criteria are repaired in the validation runner, and re-judging the same
outputs now gives 13 PASS, 1 FAIL, 3 INFO - the single FAIL being the rung-2 convergence
failure the original criteria set missed.

The sequenced follow-up work is in `improvement-plan-2026-08-27.md`.

## Recommended next steps

1. Adopt the rung-4 config. Port total 66,237 (95% CI 50,037-91,210).
2. Re-run rung 2 with a different `bss_seed` (~3 h) - the one rung without a clean isolation.
3. Fix D1 and D2 (few lines each); D3-D5 can ride along.
4. Repair the two criteria in the validation runner; add `method_selected` per component to the extractor.
5. Open the `osp_scale_is_tau` question properly - largest open item in the project, predates the batch.
6. Try `ar_escalate = TRUE` on shore all-gear; align the boat pot-closure AR map across tracks.
7. Correct the documentation: retract the tau-informativeness claim, regenerate the extrapolation figures,
   record the hours-to-trips mismatch as 3.2x rather than 4x.

## What this ladder does and does not establish

It establishes that the batch is behaviour-neutral where it claims to be, that each change moves the number
in the direction and roughly the magnitude predicted, and that the two independent tracks still agree. It
does not establish that the shipped port total is closer to the truth than the previous one - only that it
is built on fewer known errors. The boat remains ~82% extrapolated and the f = 0.3 placeholder is still the
single largest lever on it.
