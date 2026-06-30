# Boat AR Resolution Experiment: Daily vs Weekly

**Date:** 2026-06-24
**Companion to:** the pooled model (v7.4), `ADDITIONAL_OUTPUTS_PROPOSAL.md`, `CODE_IMPROVEMENTS_REVIEW_v7.0.md`.
**No em dashes per project convention.**

---

## The question

The private boat is 56% of the port harvest and rests on a thin data series (145 interviews, 24.6% of days; 195 effort observations, 49.5% of days). It runs on a weekly AR by a deliberate cap (v6.5/B1). In the production run, the boat August catch (the peak month, 26,005 of the boat's 63,864 PE catch) is where BSS and PE diverge most: BSS is 0.58 of PE, and the monthly effort columns (O7) show PE puts a large August effort peak (96,795, ~40% of the boat's annual effort) that the weekly latent does not track (25,568, essentially equal to July). The hypothesis is that the weekly AR cannot represent the August effort spike, so that variance is absorbed as negative-binomial observation noise (consistent with the T1.5 finding, boat NB share 0.83) and the peak is smoothed out of the estimate.

The experiment tests whether a finer (daily) AR recovers the August peak, and whether the data supports it.

## The constraint that shapes the experiment

The model shares one AR period grid between effort (`omega_E`) and CPUE (`omega_C`); they are not independently resolved (RMD lines ~975-981). So "daily" means daily for both. The boat's effort can plausibly support daily (49.5% coverage, ~143 effort-days, both above the daily thresholds of 0.25 coverage and 20 observations), which is why the data-driven logic already selects daily and only the cap holds it at weekly. But the boat's CPUE cannot: 145 interviews over ~71 distinct days means a 289-period daily `omega_C` is identified on only a quarter of its periods, with the rest floating on the AR prior.

This is the central tension. A clean "daily effort, weekly CPUE" is the parameterization the data actually argues for, but it requires separate grids, which is a model change. The daily experiment below is therefore diagnostic of two things at once: whether daily effort recovers the August peak, and how badly daily CPUE floats. The outcome tells us whether to pursue the split-grid model change, not just whether to flip the cap.

## How to run it

The toggle is in place (v7.4). To run the boat at daily resolution:

1. Set `ar_force = list(private_boat = "daily")` in the params list (it defaults to `NULL`, which is production behavior and leaves shore and the baseline untouched). This bypasses both the data-driven selection and the cap, and forces only the boat. The run logs `AR resolution FORCED to 'daily' for private_boat`.
2. Re-run the pooled pipeline exactly as for the v7.4 baseline. Shore ring and shore all_gear are unaffected (their resolution is unchanged), so this isolates the boat.
3. Compare the daily-boat outputs against the current weekly baseline (`output/20260624/`) using the metrics below. All of them come from diagnostics already produced (O1-O13); no new code is needed for Stage 1.

Keep this as an isolated experiment: change only `ar_force`, compare against the confirmed baseline, decide on the criteria below before looking. This is the discipline the pin episode taught.

## Stage 1 metrics, all from existing outputs

1. **Does daily converge?** `convergence_report.csv` (divergences, pass flags) and O1 `bss_full_summary_private_boat_all_gear` for `sigma_mu_E`, `sigma_eps_E`, `phi_E` Rhat and n_eff; O6 `sampler_diagnostics` for E-BFMI. A daily boat has 289 AR periods versus 44; if it cannot identify `phi_E` / `sigma_eps_E`, expect Rhat above 1.01 and divergences climbing. Convergence failure here is itself an answer: the data does not support daily.

2. **Does the latent identify, and where does it float?** O2 `bss_ar_path_private_boat_all_gear` and O3 `bss_period_coverage`. Read the `omega_E` versus `omega_C` 95% CI widths by period. The prediction is that `omega_E` (effort) identifies reasonably on sampled days while `omega_C` (CPUE) CI widths balloon on the ~75% of days with no interview. If `omega_C` floats as expected, that is the evidence that daily CPUE is unsupported and motivates the split grid. Use the new day-level coverage columns to separate "period observed" from "most days in period unsampled."

3. **Does the August peak recover?** O7 `monthly_pe_vs_bss.csv`, boat rows. Compare daily versus weekly for August (2025-08): does `BSS_effort_median` rise from ~25,568 toward the PE 96,795, and does `BSS_catch_median` rise from ~15,067 toward the PE 26,005? Also check the off-peak months do not inflate further (June was already 3.1x PE on the weekly fit). The signature of success is the seasonal shape sharpening (peak up, troughs not worse), not a uniform level shift.

4. **Predictive comparison on the reliable stream.** O13 `loo_summary_private_boat_all_gear`, the `catch` row. Compare `elpd_loo` daily versus weekly. The catch LOO is reliable for the boat (zero Pareto-k above 0.7 in the baseline), so this is a valid leave-one-interview-out comparison. Higher `elpd_loo` predicts held-out interviews better. Read the `p_loo` column too: if daily's `p_loo` jumps far above weekly's, the finer model is buying its fit with effective parameters, which is the overfitting signal.

5. **Estimate movement.** `port_total` and `catch_by_mode`. Does the boat catch and the port total move, and in which direction. The boat is 56% of the total, so a real August recovery would raise the port total above 100,420.

## Pass / fail and the interpretation tree

Decide before looking:

- **Daily wins** if it converges cleanly (boat `sigma_eps_E` / `phi_E` Rhat below 1.01, divergences not worse than weekly), the August effort and catch rise toward PE, and the catch `elpd_loo` improves by more than its SE without a large `p_loo` jump. Action: lift the cap for the boat in production, or, if CPUE `omega_C` floated badly even though effort helped, pursue the split-grid model change (daily effort, weekly CPUE) as the principled version.

- **Daily loses** if it fails to converge, or `omega_C` floats so badly that CPUE becomes noise, or August does not move, or `elpd_loo` is worse or within its SE while `p_loo` balloons. Action: the weekly cap is confirmed, and the August discrepancy is not a resolution artifact the boat data can fix. The durable levers are then the split grid only if effort identification clearly improved, or more boat summer interviews.

- **Ambiguous** if August moves the right way but `elpd_loo` is within its SE (the conditional leave-one-out cannot cleanly separate the two resolutions). This is the case where the rigorous Stage 2 below is warranted.

## Stage 2, time-block cross-validation (build only if Stage 1 is ambiguous)

Leave-one-interview-out LOO conditional on the latent path slightly favors the finer resolution, because a more flexible latent fits each held-out interview better without fully paying for temporal overfitting. The clean tie-breaker is leave-future / held-out-block prediction. Design, to be implemented and unit-tested only if Stage 1 lands ambiguous:

1. Partition the boat season into K contiguous time blocks (K = 6, roughly monthly, to put a full block in the summer peak).
2. For each block k and each resolution (daily, weekly): refit the boat with the interviews in block k removed from the catch likelihood (`IntC`, `day_IntC`, `c`, `h`, `gear_IntC`), keeping the effort data and the latent so the held-out days are informed by effort plus neighboring interviews.
3. Score the held-out block by the posterior predictive log-density of its removed interviews, computed from the refitted `lambda_C_S` and `r_C` exactly as O13 computes pointwise lpd (logsumexp over draws of `neg_binomial_2_lpmf(c_i | lambda_C[day_i] * h_i, r_C)`).
4. Sum the held-out block scores per resolution. The higher total is the better-predicting resolution, and because the held-out interviews were never in the fit, this penalizes a too-flexible latent that overfits within block.

This is 2K boat refits (the boat is the sparse, fast fit, so this is tractable). I will write and unit-test the masking and the held-out scoring before running it, rather than ship it speculatively.

## Why this staging

Stage 1 reuses metrics that already exist and the catch LOO that O13 showed is reliable for the boat, so it is a config flip and a comparison, with no untested code. It will very likely resolve the question on its own, because the predicted outcome (effort helps August, CPUE floats) is exactly the pattern that motivates the split-grid change regardless of the LOO margin. Stage 2 exists for the narrow case where the predictive comparison is too close to call, and it is built only then, tested first. That sequencing keeps every change isolated and verifiable, which is the lesson from the pin.
