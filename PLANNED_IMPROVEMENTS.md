# Planned Improvements

A living backlog for the crab creel estimation framework. We check items off as we complete them. Each item records **what** it is and **why** it is needed so the rationale is not lost between sessions.

**Primary near-term goal:** clean the repository and bring the **pooled** method to publication quality, ideally with weather impacts on CPUE incorporated.

**Current focus order:** B1 (converge the BSS) -> C1 (identify useful weather covariates) -> C2 (fold them into the production pooled CPUE).

**Status legend:** `[ ]` not started, `[~]` in progress, `[x]` done. "Done" means the code/doc change is made; a change that needs a model run to validate is marked done with a note that validation is pending.

---

## A. Repo cleanliness (quick wins)

- [ ] **A1. Untrack `.Rproj.user/`.** *What:* it is gitignored but still tracked from before the ignore existed. *Why:* IDE state should not be version-controlled. *Action:* `git rm -r --cached .Rproj.user`. Effort: trivial.
- [ ] **A2. Decide an `output/` policy.** *What:* 14 MB and ~160 files of dated run artifacts (PNGs, CSVs) are committed and growing. *Why:* a publication repo should be reproducible from code, not bloated with every run. *Options:* gitignore `output/` and keep one canonical reference run, or move runs to GitHub Releases. Effort: low (one decision).
- [ ] **A3. Reconcile a gear-resolved doc/code mismatch.** *What:* the gear-resolved v5.2 changelog and Section 8.3 say `max_divergences` default is 0, but the code uses 5. *Why:* doc must match code for a reviewer. Effort: low.
- [ ] **A4. Reproducibility scaffolding.** *What:* add an `renv` lockfile (pinned package versions) and a top-level `run_all.R` or Makefile. *Why:* a methods paper must reproduce from a clean clone. Effort: medium.
- [ ] **A5. Data-shipping decision.** *What:* decide whether `input_files/` should carry real survey data publicly. *Why:* size and any confidentiality; the usual pattern is a small synthetic example plus a data-access note. Effort: low-medium.

---

## B. Pooled method, publication quality

- [~] **B1. Make the BSS demonstrably converge, or rigorously characterize the PE fallback.**
  *What:* the boat all-gear pooled fit fails (near-total divergence) and falls back to PE; for a methods paper the BSS must converge on the components that matter, or the fallback must be rigorously characterized.
  *Why:* the paper's contribution is that the BSS fills gaps and propagates uncertainty; silent degradation to PE on a key component undermines that claim.
  *Status / sub-steps:*
  - [x] **B1.1 Diagnose.** Boat all-gear used daily AR and diverged on ~100% of iterations (treedepth 0, n_eff 76, R-hat 1.07): a funnel from a 289-state daily latent process that the weakly-informative trailer-count series cannot identify. Shore uses the same daily AR but converges (n_eff > 2000) because gear counts are far more informative.
  - [x] **B1.2 Coarsen boat AR (v6.5).** Added `ar_max_resolution` per-population cap; boat capped at weekly, cutting the latent dimension ~7x. *Validation pending a model run.*
  - [ ] **B1.3 Non-center the AR initial state `omega_0` (Stan).** The `omega_0 ~ normal(0, sqrt(sigma^2/(1 - phi^2)))` prior is a centered funnel; non-centering it improves geometry for all fits (the shore fits also carry 400-800 divergences and shore all-gear saturates treedepth at 91.6%). Do this if B1.2 plus the v6.2 tuning does not converge the boat fit, and for cleaner geometry generally. Effort: medium (Stan reparameterization, recompile).
  - [ ] **B1.4 Characterize the fallback.** If a component still cannot converge, document the PE-vs-BSS agreement and the conditions under which PE is used, so the fallback is a defensible part of the method rather than a silent failure.

- [ ] **B2. Posterior predictive checks (PPCs).** *What:* check that the effort and CPUE observation models reproduce key features of the data (zero fraction, overdispersion) by gear and day-type. *Why:* standard expectation for a Bayesian estimation paper; demonstrates the likelihood fits. Effort: medium.

- [ ] **B3. Prior sensitivity table.** *What:* vary the `R_G` prior, the `L_effective` regression prior, `R_T`, and the half-Cauchy scales; report how the harvest estimate moves. *Why:* shows robustness and answers the R_G-sensitivity recommendation from the 2026-03-31 critique that the data-driven prior only partly addressed. Effort: medium.

- [ ] **B4. External validation.** *What:* compare BSS/PE to an independent benchmark (dockside census totals, a known-total year, or the legacy estimator). *Why:* even one validation point materially strengthens the paper. Effort: medium, depends on data availability.

- [ ] **B5. Propagate pi_gear uncertainty in the pooled gear breakdown.** *What:* the pooled model applies point interview proportions post-hoc, so its gear-type catch carries no proportion uncertainty (the gear-resolved model does this properly). *Why:* either propagate it or scope the pooled gear breakdown explicitly as approximate, so the reported gear-type CIs are honest. Effort: low (doc) to medium (code).

- [ ] **B6. Sub-season boundary edge effect.** *What:* the hard Dec 1 split inflates the credible interval at the boundary (critique issue 5, still open). *Why:* either document it as a known artifact or move to a single regime-switching model. Effort: medium-high.

- [ ] **B7. Incomplete-trip handling in the pooled track.** *What:* the pooled Rmd derives `trip_status` but the explicit filter and sensitivity output exist only in the gear-resolved track. *Why:* confirm the pooled CPUE handles incomplete trips (a documented ~20% downward bias) and report a sensitivity, matching the gear-resolved treatment. Effort: low-medium.

---

## C. Weather impacts on CPUE (priority alongside B)

- [~] **C1. Run the weather-tide module and read the PSIS-LOO results.** *What:* determine which tide/weather covariates actually improve CPUE (and effort) for 2024-25. *Why:* prerequisite for integration; we should only add covariates that demonstrably help. *Status:* the module (v0.2.0) is now runnable against the crab data (the `fishing_start_time` schema mismatch is fixed via a departure proxy) and its BSS engine is at pooled v6.5 parity, so its fits and LOO comparison are trustworthy. Next action: run the GAM screen for the early read, then the BSS comparison once pooled v6.5 confirms the boat fit converges at weekly AR. Effort: low (analysis, needs a run).

- [ ] **C2. Fold supported covariates into the production pooled CPUE.** *What:* move the verified covariates from the separate `crab_bss_pooled_weather_adjusted.stan` into `crab_bss_pooled.stan` and the main Rmd. *Why:* makes weather-on-CPUE a permanent part of the pooled method rather than an experimental side file; this is the stated end goal. Effort: medium. Depends on C1.

- [ ] **C3. Reproducible environmental inputs.** *What:* the fetch caches are gitignored, so a fresh clone re-fetches from NOAA/NDBC/IEM, which can drift. *Why:* for publication, archive the exact processed daily covariate table used (commit it or snapshot the cache) so results reproduce even if the APIs change. Effort: low.

- [ ] **C4. Non-linear covariate effects (splines).** *What:* allow non-linear tide/weather effects on CPUE (the module's planned 0.3.0). *Why:* tide effects on CPUE are plausibly non-monotonic. Effort: medium. Do after C1 confirms linear effects are real.

- [ ] **C5. Weather in the gear-resolved CPUE (optional).** *What:* a gear-resolved augmented model so weather effects can differ by gear (e.g., ring-net CPUE more tide-sensitive than pot). *Why:* parity and possibly better fit, but lower priority given the pooled focus. Effort: medium-high.

- [ ] **C6. Rigorous boat-departure timing via interval-level I/E data.** *What:* the boat departure-on-flood test currently uses an interview-time proxy (departure = `interview_time - hours_fished` for completed trips, which assumes the interview is at trip end). Replace it with actual departure events from the I/E surveys, weighting each survey interval by its observed departure count and the tide phase during that interval. *Why:* removes the interview-timing assumption and uses the correct data source for departure timing. Requires exposing interval-level I/E data (`fetch_ie_data` currently aggregates to daily totals and discards per-interval timestamps) and confirming the raw `ingress_egress.xlsx` carries interval times. Effort: medium. Only needed if the departure-on-flood hypothesis is pursued seriously.

---

## D. Carried over / separate items

- [ ] **D1. Resolve the redundant baseline refit in the covariate module.** *What:* the module recomputes a baseline fit the main pooled run already produced (it needs a baseline for the PSIS-LOO comparison, but recomputes rather than reusing). *Why:* a workflow inefficiency; resolving it needs a save/load handoff between the two Rmds or merging them. Effort: medium.
- [ ] **D2. Module log_lik reconciliation (module 0.1.3, planned).** *What:* cross-check `crab_bss_pooled_weather_adjusted.stan` against `crab_bss_pooled.stan` at the `log_lik` level; confirm `p_I_shore` vs `p_TI` naming. *Why:* ensures the baseline-vs-covariate LOO comparison is valid. Effort: low-medium.
- [ ] **D3. Jetty/beach effort counts (critique issue 3).** *What:* shore effort lumps dock, jetty, and beach. *Why:* improves shore effort accuracy; requires a field-protocol change for 2025-26, not a code change. Effort: out of code scope, sampling-dependent.

---

## Completed

Recent changes already merged (most recent first):

- **Module v0.2.0** Schema fix (interview departure proxy replaces the missing `fishing_start_time`) so the module runs on the crab data, plus the BSS engine brought to pooled v6.5 parity: divergence gate, boat sampler tuning, per-population AR cap, and a per-fit convergence report. Unblocks C1.
- **Pooled v6.5 / B1.2** Per-population AR resolution cap; boat capped at weekly.
- **Pooled v6.4 + gear-resolved v5.4** R-hat convergence threshold tightened from 1.05 to 1.01 (Vehtari et al. 2021), consistently across both tracks.
- **Pooled v6.3** Documentation corrections (output-file listing, Vehtari citation).
- **Module v0.1.2** Removed the dead `bss_model_file_baseline` parameter.
- **Module v0.1.1** Weather-tide reference and file reconciliation; removed misfiled duplicate; README brought current to three approaches; `.gitignore` extended (rds, cache).
- **Pooled v6.2** Boat all-gear sampler tuning (adapt_delta 0.99, max_treedepth 13).
- **Pooled v6.1 + gear-resolved (existing v5.2 gate)** Divergence-aware convergence gate ported into the pooled track so both tracks use one standard.
