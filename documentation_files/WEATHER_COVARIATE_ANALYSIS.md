# Weather and Tide Covariate Analysis: Results and Interpretation

**Date:** 2026-06-25
**Run:** `output/20260625/pooled-CPUE-covariates`
**Model:** `stan_models/crab_bss_pooled_weather_adjusted.stan`, driver `BSS-GH-pooled-CPUE-weather-tide-covariates.Rmd`
**Season:** 2024-25, Grays Harbor / Westport
**No em dashes per project convention.**

---

## 1. Purpose

Test whether daily weather and tide covariates improve the Bayesian harvest model, by augmenting the effort and CPUE process equations with covariate terms and comparing the augmented model to the baseline by PSIS-LOO. The question is predictive: do covariates improve out-of-sample prediction enough to justify inclusion in the production estimate.

## 2. Data sources

- **Tides:** NOAA CO-OPS station 9441102 (Westport), verified, 90,470 observations.
- **Waves, wind, sea-surface temperature:** NDBC buoy 46211, 17,173 observations.
- **Air temperature, precipitation, wind, visibility, pressure:** KHQM ASOS (Hoquiam), 8,726 observations.
- Fallback hierarchy NOAA to NDBC to ASOS to GSOD, with per-source status logged.

Candidate features (~25): tide level and range (full-day and daytime-only), daytime high-tide hour, flood hours, number of daytime high tides, a tidal-energy index; wind and gust; wave height mean and max; sea-surface temperature; barometric pressure; air temperature min and max; total precipitation; visibility; and three departure-timing behavioral features (proportion departing on flood, near high, mean departure tide height).

## 3. Method

Two layers:

1. **GAM screening (Layer A/B).** Generalized additive models with smooths on the candidate features identify which have a non-flat association with effort and CPUE, per population. The screened feature sets (`layer_b_selected_features.csv`) were, for shore effort, daytime tide range, tidal energy, and precipitation; for shore CPUE, daytime tide range, daytime high hour, and wave height; for boat effort, tide range, daytime tide range, number of daytime high tides, tidal energy, wave height, precipitation, and max temperature (7 features); for boat CPUE, eight tide, wave, SST, and precipitation features.

2. **Bayesian model plus PSIS-LOO (selection).** The screened features enter the effort and CPUE log-rate equations as linear covariate terms. The augmented model is compared to baseline by `loo::loo_compare`, and a feature set is included only if it improves elpd by more than 4.0 times the standard error of the difference. Pareto-k is checked for LOO reliability.

## 4. Results

### 4.1 The decision: covariates excluded for all three components

| Component | elpd_diff (covariates vs baseline) | SE | margin | Decision |
|---|---|---|---|---|
| shore ring_net | -11.4 | 4.7 | 2.4 SE worse | exclude |
| shore all_gear | -0.03 | 3.3 | tied | exclude |
| private boat all_gear | -26.6 | 9.5 | 2.8 SE worse | exclude |

No component clears the 4.0-SE improvement margin. shore all_gear is a dead tie; shore ring and the boat are meaningfully worse out-of-sample with covariates. The LOO is reliable for the comparison: the covariate and baseline models have nearly identical Pareto-k profiles (shore all_gear 57 vs 52 observations with k > 0.7, boat 4 vs 5), so the differences are trustworthy, including the negative boat result.

### 4.2 The coefficients: weather drives effort, not CPUE

Several individual coefficients are physically sensible and have 95% credible intervals excluding zero, and they are concentrated in the **effort** equations:

- **Boat effort:** wave height -0.65 (-0.89, -0.41), max temperature +0.68 (+0.41, +0.95), tide range +0.37, number of daytime high tides +0.20, precipitation -0.27. Rough water and rain suppress boat effort; warm weather and bigger tides raise it.
- **Shore effort:** precipitation -0.25 (ring) and -0.22 (all_gear); daytime tide range +0.37 (ring). Rain suppresses shore effort; bigger daytime tides raise it.
- **CPUE:** nearly flat. Only shore all_gear daytime high hour (+0.16) and a marginal boat SST (-0.29, CI just excluding zero) are distinguishable from zero. Catch-per-unit-effort is essentially weather-insensitive.

### 4.3 Convergence held

Adding covariates did not destabilize the fits. shore all_gear, the funnel-prone component, had fewer divergences with covariates (444 to 394), the boat was stable (57 to 66), and all Rhat are ~1.00. The covariate parameters did not aggravate the `sigma_mu_E` geometry. (The module's `pass_convergence` column reads FALSE for all fits including baselines that pass the production gate, so it uses a stricter flag than the main pipeline and does not indicate a regression.)

### 4.4 A false-precision warning

The covariate models produce **narrower** credible intervals (shore ring catch -12%, effort -16%; boat catch -16%) while predicting held-out data **worse** by LOO. Narrower-but-worse is the signature of overfitting and overconfidence, not a real precision gain. This is why model selection here uses LOO and not interval width: judged on interval width alone, the tighter boat catch CI would wrongly argue for inclusion. The final estimates correctly use the baseline (`model_used = "baseline"`), so the production numbers are not contaminated.

## 5. Interpretation: why significant effects do not improve prediction

The result looks paradoxical, significant weather-effort coefficients yet no predictive value, but it is explained by what the model already observes. Effort is estimated from the effort counts (gear and trailer counts on sampled days). On any sampled day, the count is a direct observation of effort, so it already encodes the effect of weather: the count is low on a rainy or rough day whether or not weather is in the model. Weather can therefore only add value by improving the **interpolation of effort on days with no count**, and the AR(1) latent process already performs that interpolation from neighboring days. Weather competes with the AR for the interpolation and, on routine data, does not beat it. CPUE is the one place covariates would have to earn their keep on their own, because catchability is not directly observed the way effort counts observe effort, but CPUE turns out to be weather-flat. So covariates have no routine predictive role, and on the sparse boat (15 features on 145 interviews and 195 effort observations) they overfit, which is why the boat is the worst at -26.6.

## 6. Consistency note

The covariate module's boat estimate (catch 12,400, effort 19,358) is ~4.5x below the v7.4 production boat (catch 56,266, effort 87,359), with the CPUE identical (0.64), so it is purely an effort-scale difference. This points to the boat gear-hours expansion (trailer count times gear-per-group times 24 hours) not being applied identically in the two pipelines; the covariate Rmd predates the current model. This does not affect the covariate decision (the LOO is computed on the latent fit and is invariant to the post-hoc expansion, and the baseline-versus-covariate comparison is internally consistent) or the production estimate (the main pipeline is authoritative and covariates were excluded). It should be reconciled so the module's estimates match and the 12,400 is not mistaken for a boat harvest.

## 7. Outstanding questions

### 7.1 Would covariates help ground effort on an unsampled period?

This is the one scenario where the section-5 redundancy argument breaks, and the answer is a qualified yes.

The reason covariates fail on routine data is that the effort counts already encode the weather effect. But if a week has **no counts at all** (a sampling gap), that redundancy disappears. The only information about that week's effort is then the AR interpolation from neighboring weeks plus, potentially, weather. The AR carries the **level** (persistence from neighbors) but cannot know that the unsampled week had, say, a major storm; a wave-height covariate (-0.65 for the boat) would pull the estimate down for exactly that week. So covariates provide information that is independent of the absent counts, and they are most valuable when the unsampled week's weather is **anomalous**, because that is when the AR's naive interpolation is most wrong. For a typical unsampled week the AR alone is adequate; for an unusual one, weather is a correction.

Two important qualifications:

- **The LOO in this run does not test this scenario and so does not rule it out.** PSIS-LOO leaves out one observation at a time, and a single held-out effort count is easily predicted by the AR plus the other counts in the same week. The aggregate LOO therefore cannot see the unsampled-week benefit, because it never holds out a full week. The correct test for the sampling-gap use case is **leave-one-week-out block cross-validation** (hold out every observation in a contiguous week, predict the week from the rest), which is the same time-block CV designed for the boat resolution test. A covariate set that loses on leave-one-observation-out LOO can still win on leave-one-week-out CV, and that is the experiment that would answer this question.

- **The vehicle should be parsimonious, not the 15-feature model.** For grounding gaps you want the few strong, well-identified effort drivers (boat wave height and temperature; shore precipitation), not the full screened set that overfits. A model with two or three robust covariates per population is the right tool, both because it will not overfit and because the coefficients you apply to an unsampled week must themselves be trustworthy.

**Bottom line:** weather is not worth including in the routine production model, but it is worth keeping a parsimonious weather-effort model on the shelf as a contingency for grounding effort when a period goes unsampled, evaluated by leave-one-week-out block CV rather than by the observation-level LOO used here. It is insurance against sampling gaps, especially gaps that coincide with anomalous weather, not a routine improvement.

### 7.2 Would expanding to the full multi-year dataset help resolve the story?

It would help the covariate question specifically, but only if the 2024-25 protocol change is handled, and that change is a serious confound, not a detail.

**The benefit.** More seasons mean more weather variation and more data, which would tighten the weather-effort coefficients and could either confirm or refute the effects that the single-season analysis can only weakly resolve. If the goal is robust weather-effort slopes (which is exactly what 7.1 needs to ground unsampled periods), multi-year estimation is the way to get them.

**The confound.** Before 2024-25, the protocol was a single effort count per day at the estimated peak effort time; in 2024-25 it became three randomized counts per day. These measure different quantities. A single peak count is an upward-biased estimate of a different estimand (peak effort, not mean daily effort), while the randomized 2024-25 counts give unbiased mean daily effort, which is what the gear-hours expansion assumes. Pooling the years naively would feed peak counts into a model that treats counts as observations of mean daily effort, biasing the effort level high in the early years. So the early years cannot simply be added.

**Why the slopes can still be pooled, with one condition.** The effort model is log-linear: `log(effort) = level + B_weather * weather + ...`. If the peak-to-mean ratio is roughly constant (peak is about a fixed multiple of the mean), then `log(peak) = log(constant) + log(mean)`, and the weather **slopes** `B_weather` are identical for peak counts and mean counts; the protocol difference is a pure level shift. So the early years can contribute to estimating the weather slopes if the model includes a **protocol fixed effect** (a separate intercept or level for the peak-count years). The condition is that the peak-to-mean ratio not itself depend on weather. If effort is more concentrated in time on calm, pleasant days (so peak/mean is higher then), there is a protocol-by-weather interaction and the simple fixed effect is not enough; that interaction would have to be modeled or the slopes would be biased.

**The 2024-25 data is the calibration key.** The randomized three-count protocol is what makes it possible to characterize the within-day effort distribution and therefore to estimate the peak-to-mean relationship needed to bring the earlier years onto a common footing. The new protocol does not just give better current data; it provides the information to correct the old data. This is the natural bridge for any multi-year extension, and it is also independently valuable for the main harvest model, since a peak-to-mean correction is what the earlier years need to be usable for production estimates at all.

**Bottom line:** expanding to all years is worthwhile for robustly estimating the weather-effort slopes, which is the quantity that matters for the sampling-gap use case in 7.1, but it requires a protocol fixed effect at minimum, a check on whether the peak-to-mean ratio is weather-dependent, and use of the 2024-25 within-day data to calibrate the correction. For the production harvest estimate, pooling the early years is a larger undertaking because the peak counts need that full peak-to-mean correction before they can be trusted as observations of mean daily effort.

## 8. Recommendations

1. **Accept the exclusion for routine estimation and report it as a clean negative result.** Covariates were screened with a GAM and tested with PSIS-LOO, and weather does not improve harvest prediction because effort is already measured by counts and CPUE is weather-insensitive.
2. **Report the effort-driver coefficients descriptively.** Boat effort is strongly suppressed by waves and raised by temperature, shore effort by rain. These are effort-dynamics findings for the discussion and bear on sampling design, since effort is predictably low on rough or rainy days.
3. **Reconcile the boat effort expansion** between the covariate module and the v7.4 main pipeline.
4. **If the sampling-gap use case (7.1) matters, build a parsimonious weather-effort model and evaluate it by leave-one-week-out block CV**, not observation-level LOO, using only the strong effort drivers. This is the targeted experiment that the current run cannot substitute for.
5. **A multi-year weather analysis (7.2) is the way to get robust slopes for (4), but design it with a protocol fixed effect, test the peak-to-mean ratio for weather dependence, and use the 2024-25 within-day data to calibrate the early-year correction.** Treat this as a distinct study, not a re-run of the single-season analysis.

## 9. Version history

- 2026-06-25: initial documentation of the 2024-25 single-season weather/tide covariate run and the two outstanding questions (unsampled-period grounding; multi-year expansion under the protocol change).
