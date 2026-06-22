# Additional Run Outputs to Persist: a Prioritized Catalog

**Companion to:** `CODE_IMPROVEMENTS_REVIEW_v7.0.md` and the pooled model (v7.1).
**Date:** 2026-06-22
**Purpose:** identify quantities the run computes and holds in memory but does not save, so they can be persisted (CSV or plot) and you rarely need to keep the R workspace to diagnose. Each item is tied to a planned-improvement ID (T1.x and the critique issues) where relevant, and ranked by value.
**No em dashes per project convention.**

---

## How this was scoped

I compared what the run persists (every `write.csv`, `ggsave`, `writeLines`, plus the `model_diagnostics.R` and v7.1 effort-overdispersion writes) against what it computes and keeps in memory (the per-fit `stanfit` objects and their full posterior, the AR latent path, the per-observation fitted values, the summed-quantity draws, and the PE and data objects). The gaps below are quantities already in the workspace after a run; persisting them is additive and changes no estimate. Every candidate is a small CSV unless noted, so storage is not the binding constraint; the only heavy option is O12 (saving the fit objects), called out separately.

A note on the trade-off you are managing: more saved artifacts improve diagnosability but also grow `output/`, which is already a flagged repo-hygiene problem (review item T4.2, `output/` is 410 tracked files). The Tier A set below is small CSVs; O12 is the one that would bloat the repo, so it is recommended off-git.

---

## Tier A: build before the next run (high leverage for the planned improvements and for general diagnosis)

### O1. Full parameter posterior summary per fit `bss_full_summary_<label>.csv`
**What.** The complete `rstan::summary(fit)$summary`: every parameter (the AR path `omega_E`/`omega_C`, the innovations `eps_E`/`eps_C`, `mu_E`/`mu_C`, all scales, `R_G`, `R_T`, the derived quantities) with mean, sd, the 2.5/25/50/75/97.5 quantiles, n_eff, and R-hat. The current `structural_params_<label>.csv` is a curated ~19-parameter subset; this is everything.
**Why.** The single highest-value "never need the workspace" artifact. It captures every parameter's posterior and convergence in one file, so any future question about a parameter you did not pre-select is answerable from disk. Directly supports T2.5 (the `sigma_mu_E` identification), the effort-overdispersion work (the `sigma_r_E`/`r_E` posterior in context), and any convergence post-mortem.
**Format.** CSV. **Effort.** Trivial. **Size.** Small (a few hundred rows for shore all-gear).

### O2. AR latent path per fit `bss_ar_path_<label>.csv`
**What.** The posterior of the AR(1) deviations `omega_E` and `omega_C` per AR period (period index, the date range the period spans, median and 95 percent CI for each), plus the period-level `lambda_E` and `lambda_C`.
**Why.** The latent AR process is the heart of the state-space model and is currently invisible in the outputs. It is needed to see where the AR is anchored by data versus extrapolating (T1.1, the summer-extrapolation transparency), to characterize the sub-season boundary edge effect (T3.1 / B6), and to diagnose AR dynamics generally. The daily catch and effort CSVs show `C` and `E` but not the underlying process that generated them.
**Format.** CSV, optionally a companion plot of the path with its CI. **Effort.** Low. **Size.** Small (P_n rows; 289 for shore all-gear daily).

### O3. AR-period data coverage and extrapolation map per fit `bss_period_coverage_<label>.csv`
**What.** Per AR period: the number of effort observations, the number of interviews, an "observed vs projected" flag, and the `omega` posterior CI width as a projection-uncertainty proxy. Built from `bss_data$period`, `day_Gear`/`day_T`, and `day_IntC`.
**Why.** The most direct possible support for T1.1. It quantifies, per period, how much of the estimate rests on sampled versus unsampled time, which is the headline transparency the paper needs (the harvest is dominated by sparsely sampled summer months). It turns the summer-extrapolation argument from a narrative into a table.
**Format.** CSV. **Effort.** Low. **Size.** Small.

### O4. Modeled daily CPUE trajectory per fit `bss_daily_cpue_<label>.csv`
**What.** The modeled CPUE rate `lambda_C` per day with 95 percent CI, next to the raw interview CPUE for that day and the interview count. The interview likelihood is `c ~ NB(lambda_C * h, r_C)`, so `lambda_C` is the model's CPUE.
**Why.** The current `plot_cpue_timeseries.png` plots raw interview CPUE (`dungeness_kept / fishing_time_total`), not the modeled `lambda_C`. The model's AR-smoothed, day-type-adjusted CPUE is invisible in the outputs. It is needed for the weather-on-CPUE work (C1/C2, which acts on `lambda_C`) and to see the `B1_C` weekend CPUE effect in the daily series.
**Format.** CSV, optionally a modeled-vs-raw overlay plot. **Effort.** Low. **Size.** Small.

### O5. Per-observation PPC residuals per fit `ppc_byobs_<label>.csv`
**What.** For each effort count and each interview catch: the observed value, the fitted mean (and CI), the PIT value, and the day and gear. The PPC already computes these internally in `bss_ppc_calibration`; this exposes them instead of discarding them after the aggregate.
**Why.** The PPC currently saves only aggregate coverage and PIT moments. Per-observation residuals enable the residual analyses the review calls for: which gears and days fit poorly, the bag-limit upper-tail misfit (critique issue 9, seen in the shore catch PIT), the incomplete-trip CPUE bias (T2.3), and interview representativeness (critique issue 10). It is the difference between knowing the model is mildly miscalibrated and knowing where.
**Format.** CSV. **Effort.** Low to medium (recompute or expose from the PPC path). **Size.** Small (observation counts in the low thousands).

### O6. HMC sampler diagnostics per fit `sampler_diagnostics_<label>.csv`
**What.** Per chain: mean `accept_stat`, mean and max treedepth, mean `n_leapfrog`, divergent count, and E-BFMI (the energy-based fraction of missing information). From `get_sampler_params`, which the run already calls for divergences.
**Why.** `convergence_report.csv` records divergences and treedepth but not E-BFMI, which is the key energy diagnostic for the AR-funnel pathologies this model has fought (Betancourt 2017); a low E-BFMI flags exactly the slow mixing the boat and B1.7 episodes involved. Cheap insurance for any future convergence post-mortem.
**Format.** CSV. **Effort.** Low. **Size.** Tiny (one row per chain).

---

## Tier B: strong secondary value

### O7. Monthly PE-vs-BSS by component `monthly_pe_vs_bss.csv`
**What.** Per month and component: PE catch and effort next to BSS catch and effort. The current `pe_vs_bss_comparison.csv` is by sub-season only.
**Why.** T1.1 needs the by-month comparison to locate and explain the shore (+50 percent) versus boat (-12 percent) BSS-vs-PE divergence, which the run output and the review both flag. A by-sub-season table cannot show that the divergence concentrates in summer; a by-month table can.
**Format.** CSV, optionally an overlay plot. **Effort.** Low to medium. **Size.** Small.

### O8. Summed-quantity posterior draws per fit `bss_draws_summed_<label>.csv`
**What.** The draws of `C_sum`, `C_expected_sum`, and `E_sum` (optionally thinned to ~2,000). These live in `bss_all` in memory but are never written.
**Why.** Lets you recompute any quantile, HDI, or the gate impact metric offline without re-fitting. It is the cheapest way to make the gate analysis and any interval re-derivation reproducible from disk, which is squarely the workspace-independence you asked for.
**Format.** CSV. **Effort.** Trivial. **Size.** Moderate but small in absolute terms (a few thousand rows by three columns).

### O9. Prior-vs-posterior comparison per fit `prior_vs_posterior_<label>.csv`
**What.** For the key parameters (`R_G`, `R_T`, the `sigma_eps`, `sigma_r`, `phi`, `sigma_mu`, `B1`, `B2`, `B1_C`), the prior family with its mean and sd next to the posterior mean, sd, and CI, plus a "prior pull" flag when the posterior sits near the prior.
**Why.** Supports T1.3 (prior sensitivity) and surfaces weakly-identified parameters pulled by their priors (for example `R_T` floats at its prior for shore fits, where there are no trailer counts). It makes the prior sensitivity a reading rather than a re-derivation.
**Format.** CSV. **Effort.** Medium (the priors must be encoded alongside the extraction). **Size.** Small.

### O10. Gear-proportion table `gear_proportions.csv`
**What.** The `pi_gear` point proportions per period and day-type and gear, with the interview counts behind them.
**Why.** Supports B5 (the pooled gear breakdown carries no proportion uncertainty). `catch_by_gear_type.csv` saves the gear catch but not the proportions or their support, so the approximation cannot currently be audited from disk.
**Format.** CSV. **Effort.** Low. **Size.** Small.

### O11. Per-fit input and data summary `fit_data_summary.csv`
**What.** One row per fit: population, sub-season, AR resolution, P_n, D, number of gear observations, number of trailer observations, number of interviews, number of I/E observations, the observed date range, the percent of days with an effort observation, the percent of days with an interview, and mean interviews per month.
**Why.** Reproducibility plus a compact, per-fit view of the data density that T1.1 is about. `season_summary.csv` has aggregate counts but not the per-fit breakdown that explains why the boat ring-net falls to PE or why summer is uncertain.
**Format.** CSV. **Effort.** Low. **Size.** Tiny.

---

## Tier C: optional, and the one heavy item

### O12. Saved fit objects per fit `fit_<label>.rds`
**What.** The full `stanfit` object per fit.
**Why and the trade-off.** This is the ultimate "diagnose anything later," because you can re-extract any quantity without re-running. But it is the one expensive option: a shore all-gear `stanfit` with 10,000 draws over several hundred parameters can run to 100-plus MB, and committing it would worsen the `output/` bloat the review already flags (T4.2). With O1, O2, O5, and O8 saved as CSV, almost everything you would reach into the `.rds` for is already on disk, so the `.rds` is largely redundant for the planned improvements. Recommendation: if you want it, save it locally or to a GitHub Release rather than committing it, or save a thinned draws object (a few hundred draws) instead of the full fit.
**Format.** RDS. **Effort.** Trivial to write, but see the size trade-off. **Size.** Large.

---

## Summary and recommendation

| ID | Output | Primary improvement served | Effort | Size |
|---|---|---|---|---|
| O1 | Full parameter summary per fit | General; T2.5; convergence | Trivial | Small |
| O2 | AR latent path per fit | T1.1; T3.1/B6 | Low | Small |
| O3 | AR-period data coverage map | T1.1 | Low | Small |
| O4 | Modeled daily CPUE per fit | C1/C2; B1_C | Low | Small |
| O5 | Per-observation PPC residuals | T2.3; issues 9, 10 | Low-med | Small |
| O6 | HMC sampler diagnostics (E-BFMI) | General; T2.5 | Low | Tiny |
| O7 | Monthly PE-vs-BSS by component | T1.1 | Low-med | Small |
| O8 | Summed-quantity draws per fit | Gate re-derivation; general | Trivial | Moderate |
| O9 | Prior-vs-posterior per fit | T1.3 | Medium | Small |
| O10 | Gear-proportion table | B5 | Low | Small |
| O11 | Per-fit data summary | Reproducibility; T1.1 | Low | Tiny |
| O12 | Saved fit objects (.rds) | Everything, but heavy | Trivial | Large |

**Recommended build before the next run: O1 through O6, plus O8 and O11.** That set is all small CSVs and one tiny table, it covers the workspace-independence you asked for (full parameter summary, AR path, modeled CPUE, per-observation residuals, sampler diagnostics, and the summed-quantity draws), and it directly equips the three highest Tier-1 improvements (T1.1 summer extrapolation via O2/O3/O7/O11, the effort and residual work via O5/O6, and general diagnosis via O1/O8). I would hold O7, O9, and O10 only if you want them now (O7 is the strongest of the three for T1.1), and I would not commit O12.

**Implementation approach.** One additive file, `R_functions/save_run_diagnostics.R`, holding a per-fit writer (O1 to O6, O8) and a run-level writer (O7, O11), each `tryCatch`-wrapped, called from a single chunk after section 7.13, exactly mirroring the `model_diagnostics.R` and effort-overdispersion patterns already in the repo. This keeps the RMD change to a few lines and the logic in a sourced, testable file. It writes only CSVs and changes nothing in the model, so it is safe to add before the next run.

Tell me the scope (the recommended set, or adjust it) and I will build the module and the chunk, with the same parse-check and math-verification discipline as the effort-overdispersion diagnostic.
