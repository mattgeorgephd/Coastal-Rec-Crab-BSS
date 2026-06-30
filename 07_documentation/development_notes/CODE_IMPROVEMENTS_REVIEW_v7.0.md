# Outstanding Code Improvements: Comprehensive Impact-Ranked Review

**Repository:** `FWC-estimation-method` (pooled CPUE track, primary)
**Review date:** 2026-06-22
**Baseline reviewed:** code at v6.9.1 / Stan v6.9, plus the v7.0 changes delivered this session (scale-aware gate, PE monthly day-length fix, PPC extraction fix) which are pending their first validating run.
**Basis of review:** the full repo (Stan models, R functions, the three Rmd tracks), the documentation folder (`PLANNED_IMPROVEMENTS.md`, the B1.5/B1.6 change notes, the diagnostics-and-reproducibility notes, the 2026-03-31 model critique, and all three documentation markdown files), and the committed run outputs from `output/20260331` through `output/20260621` including convergence reports, PE-vs-BSS comparisons, monthly estimates, structural-parameter tables, divergence-localization tables, and the covariate-module LOO and GAM outputs.

This document supersedes `PLANNED_IMPROVEMENTS.md` as the current ranked backlog. It cross-references the existing backlog IDs (A1 to D4) and the 2026-03-31 critique issue numbers (1 to 11) so nothing is lost, marks what v7.0 changed, and adds items found in this review. Appendix A is a one-row-per-item status reconciliation of every prior ID and critique issue.

---

## How to read this

**Impact tiers** rank by effect on the harvest estimate's accuracy and on the credibility of the methods paper, weighted by the probability the issue materially matters for the 2024-25 result. Tier 1 items can change the headline number or are blocking requirements a reviewer will demand. Tier 4 items are correctness-neutral cleanup.

**Status legend:** `OPEN` not started; `PARTIAL` started, not complete; `UNBLOCKED` a dependency just cleared so it is now actionable; `NEW` first raised here; `RESOLVED (pending run)` fixed in code but awaiting a validating run.

**Effort:** trivial (minutes), low (under a day), medium (a few days), high (a week or more or needs new data).

**A note on what v7.0 already cleared.** The previous Tier 1 blocker in `PLANNED_IMPROVEMENTS.md` was B1: "make the BSS demonstrably converge or rigorously characterize the PE fallback," because the boat all-gear fit was reported as failing. That is now substantially resolved. By the v6.9 run the boat all-gear fit was converging cleanly (59 of 10,000 divergent, 0.6 percent; R-hat 1.001; n_eff 5,266 on catch, 9,520 on effort), and v7.0's scale-aware gate moves it onto the BSS posterior. The only remaining PE fallback is the boat ring-net sub-season, which is data-limited (17 interviews, below the 20-interview threshold), not a geometry failure. B1.1 through B1.6 are done; B1.4 is satisfied for the one component still on PE. v7.0 also unblocked B2 (posterior predictive checks), which could not run because of the PPC extraction bug. The tiers below reflect this post-v7.0 state.

---

## Tier 1: Critical for a defensible published estimate

These are the items most likely to change the reported number or to be non-negotiable requirements for publication. Ordered by priority.

### T1.1 Validate and make transparent the summer-extrapolation dominance of the harvest estimate
**(NEW; sharpens and partly subsumes B4 external validation and the 2025-26 coverage items)**

**What.** The reported harvest is dominated by the three most thinly sampled months, and rests heavily on the AR(1) process extrapolating effort and CPUE into them. The paper must quantify how much of the total is extrapolation rather than expansion of observed data, characterize the resulting uncertainty honestly, and validate the extrapolated portion against any available benchmark.

**Why it matters most.** This is the single largest threat to the credibility of the headline number. From `output/20260621/pooled-CPUE/monthly_estimates.csv`, July (19,795) and August (17,806) alone are 37.4 percent of the 100,448-crab port total, and July through September are 45.7 percent. These months carry 3 to 4-fold credible intervals (July 10,139 to 42,262; August 9,610 to 32,734), far wider than the winter months (December 9,326 to 11,586). Per the project record, summer boat interview density is roughly 6 to 16 per month. The AR(1) process, by design, resists a late-season CPUE decline when interview density is low, so it projects a high summer effort and CPUE that the data only weakly constrain. The contribution this makes is exactly the BSS-vs-PE divergence visible in `pe_vs_bss_comparison.csv`: shore all-gear BSS catch is +50.1 percent over PE (19,116 vs 12,733) and effort +30.8 percent, while the boat all-gear BSS is -11.9 percent under PE. These two large gaps point in opposite directions, which is itself a result a reviewer will demand an explanation for. A paper whose central estimate is 40-plus percent driven by months with a handful of interviews, without an explicit accounting of that dependence, will not survive review.

**Concrete actions.**
1. Add an output that decomposes the port total into the share coming from months (or sub-season strata) above versus below an interview-density threshold, so the fraction of harvest resting on sparse data is reported, not implicit.
2. Add a sampled-vs-extrapolated diagnostic per fit: the share of effort-days and interview-days that fall in each AR period, so the reader can see which periods are anchored by data and which are projected.
3. Produce the mechanistic explanation of the shore (+50 percent) versus boat (-12 percent) BSS-vs-PE divergence. The likely story is that for shore the BSS fills unsampled summer days with effort the PE never expands into, while for the boat the trailer-count series and the gear-hours correction pull the BSS below the PE's trailer expansion. Confirm this against the daily effort and CPUE series (`bss_daily_effort_*`, `bss_daily_catch_*`) and write it up.
4. Pursue at least one external validation point (B4): dockside census totals, a known-total year, or the legacy estimator. Even one anchor materially strengthens the paper, and it is the only thing that can adjudicate whether the summer extrapolation is right.

**Effort.** Medium for the transparency diagnostics and the divergence write-up; medium-to-high for external validation depending on data availability.
**Depends on.** Nothing for the diagnostics; benchmark data for validation.

### T1.2 Run and report posterior predictive checks; assess zero-inflation
**(backlog B2, UNBLOCKED by v7.0; folds in critique issue 9)**

**What.** Now that the PPC extraction bug is fixed, run the posterior predictive checks for effort counts and interview catches by gear and day-type, report central-interval coverage and PIT calibration, and use them to decide whether the negative-binomial catch likelihood needs a zero-inflated form.

**Why it matters.** A Bayesian estimation paper is expected to demonstrate that its observation models reproduce the data's key features (zero fraction, overdispersion). Nothing in the framework currently tests this. The run outputs confirm the gap: `output/20260621/pooled-CPUE/` contains `structural_params_*` and `divergence_localization_*` for every fit but no `ppc_calibration_*` or `ppc_pit_*` files, because `bss_ppc_calibration` aborted on every fit (the non-conformable-arrays bug v7.0 fixes). The PPC is also the correct vehicle for critique issue 9 (zero-inflation): the catch likelihood is `neg_binomial_2(lambda_C * h, r_C)` (no zero inflation), and the 2026-03-31 critique flagged that a high zero fraction, especially for snare and ring-net gear, may stretch `r_C`. A ZINB variant already exists in the freshwater lineage. The PIT histograms will show directly whether the NB is mis-calibrated at zero; if so, ZINB is the documented next step. Decide this on evidence rather than assume.

**Concrete actions.** Confirm `ppc_calibration_*` and `ppc_pit_*` now write on the v7.0 run; read the PIT shapes and coverage; report them as a standard output; if the zero bin is systematically off, prototype ZINB on the worst-calibrated component and compare via PSIS-LOO.

**Effort.** Low to read and report (the code now runs); medium if ZINB is warranted.
**Depends on.** The v7.0 PPC fix (delivered, pending run).

### T1.3 Prior sensitivity analysis, with R_G as the priority
**(backlog B3; critique issues 1 and 2)**

**What.** Vary the consequential priors and report how the harvest estimate moves: the `R_G` prior first, then the `L_effective` regression prior, `R_T`, and the half-Cauchy scales.

**Why it matters.** `R_G` (gear per crabber, modeled as `lognormal(log(R_G_prior_mu), R_G_prior_sigma)`) and the day-length expansion are the two multiplicative factors the 2026-03-31 critique identified as the most consequential assumptions in the entire framework, sitting upstream of both PE and BSS. The data-driven prior introduced for `R_G` only partly addresses the critique's recommendation; the structural-parameter table shows `R_G` is tightly identified at 1.271 (95 percent CI 1.243 to 1.299), but that is the posterior under one prior, and the critique explicitly asked for the estimate's movement under `R_G` fixed at 1.0, 1.27, and 1.5. A prior sensitivity table is a standard robustness expectation and the direct answer to a recommendation that is still formally open. The `L_effective` prior matters for the same reason: it sets the day-length expansion that multiplies through everything.

**Effort.** Medium (re-runs across a small grid of priors; can reuse the existing fit machinery).
**Depends on.** Nothing.

### T1.4 Commercial/charter: separate per-vessel-type catch means
**(NEW, confirming critique issue 6)**

**What.** Replace the single mean-catch-per-vessel applied to all vessel-days with separate commercial and charter per-vessel means, applied to the respective tally columns.

**Why it matters.** The commercial/charter component is 12,007 crab, about 12 percent of the port total, and is currently computed with a known structural simplification. In `estimate_comm_charter` (Rmd around line 1013) the estimator is `mean_dung_per_vessel <- sum(comm_int$dungeness_kept) / nrow(comm_int)` applied to `total_comm_charter = commercial_tally + charter_tally`. This lumps two vessel classes with very different catch profiles (a commercial crab vessel running many pots recreationally versus a six-pack charter). If the commercial-to-charter mix shifts across the season, the constant mean introduces bias. The fix is fully supported by the data: the interview file carries a `boat_type` field with 478 "Commercial Crab Vessel," 119 "Charter," and 10 "Guide" records, and the tally file (`wes_commercial_tally.csv`) already splits `commercial_tally` and `charter_tally`. Computing two means and applying them to the matching tally columns is a small, well-scoped change that could move the 12,007 figure meaningfully. Fold the "Guide" records into charter and normalize the data typo ("Commerical").

**Effort.** Low to medium (confirm the boat_type-to-tally-column mapping and sample sizes per stratum, then split the mean).
**Depends on.** Nothing.

---

## Tier 2: High value for accuracy or for closing a parity gap

### T2.1 Resolve the `catch_by_mode.csv` shore PE/BSS inconsistency
**(carried from v7.0 delivery; pre-existing)**

**What.** `catch_by_mode.csv` reports the shore catch on a PE basis (about 18,803) while the port total and monthly outputs report shore on BSS (about 25,849: ring-net 6,733 plus all-gear 19,116). Reconcile so all reported shore figures use the selected method.

**Why it matters.** This is an internal inconsistency in the published outputs: two committed files disagree on the shore number by about 7,000 crab. It shares a root with the by-mode method selection and should be fixed alongside it. v7.0 makes it more visible, because the boat all-gear row in the comparison and convergence outputs now correctly reads BSS, which draws attention to the remaining mode-table mismatch. A reviewer or a downstream user reading `catch_by_mode.csv` against `port_total` will see the discrepancy immediately.

**Effort.** Low.
**Depends on.** Best done after the v7.0 run confirms the new selection logic.

### T2.2 Propagate gear-proportion (pi_gear) uncertainty in the pooled gear breakdown
**(backlog B5; gear-resolved parity)**

**What.** The pooled model produces `catch_by_gear_type.csv` by applying point interview proportions post-hoc, so its gear-type catch carries no proportion uncertainty. Either propagate it (as the gear-resolved track does via Dirichlet sampling in generated quantities) or scope the pooled gear breakdown explicitly as approximate.

**Why it matters.** The gear-resolved track already solved this (its documentation records "pi_gear does not contribute posterior uncertainty" as addressed in v5.2 via Dirichlet sampling). The pooled track's gear-type credible intervals are therefore too narrow, which is a correctness issue if those intervals are reported. The cheap honest fix is a one-paragraph scoping note that the pooled gear breakdown is a point decomposition; the full fix mirrors the gear-resolved Dirichlet propagation.

**Effort.** Low for the scoping note, medium for full propagation.
**Depends on.** Nothing.

### T2.3 Explicit incomplete-trip filter and sensitivity in the pooled track
**(backlog B7; gear-resolved parity)**

**What.** The pooled Rmd derives `trip_status` but the explicit incomplete-trip filter and its sensitivity diagnostic exist only in the gear-resolved track. Add the filter and a sensitivity output to the pooled track.

**Why it matters.** Incomplete trips carry a documented downward CPUE bias of roughly 20 percent (gear-dependent), and the gear-resolved record notes the incomplete fraction is about 35 percent for shore crabbers and about 9 percent for the boat. If the pooled CPUE is not filtering incomplete trips the same way, its shore CPUE (and therefore the shore harvest that dominates the +50 percent BSS-vs-PE gap in T1.1) could be biased low or be inconsistent with the gear-resolved track. At minimum the two tracks should handle this identically and the pooled track should report the same sensitivity.

**Effort.** Low to medium.
**Depends on.** Nothing.

### T2.4 Complete the weather-covariate evaluation across all components
**(backlog C1, PARTIAL)**

**What.** Finish the PSIS-LOO baseline-vs-covariate comparison for shore all-gear and the boat, not only shore ring-net, then apply the decision rule per population and sub-season.

**Why it matters.** The covariate module ran the GAM screen and produced selected feature lists for all four population-by-component cells (`layer_b_selected_features.csv`), but the committed LOO output is only `loo_comparison_shore_ring_net_only_Dungeness_Kept.csv`. That single result is marginal: covariates improve ELPD by 11.4 with a standard error of 4.74, about 2.4 SE, which does not clear the module's own 4-SE inclusion bar. The component that matters most for harvest, shore all-gear, was not evaluated, and it is precisely the component with the summer-extrapolation problem in T1.1, where a real CPUE covariate (the GAM screen finds daytime tide range, tidal energy, and wave height significant for shore CPUE) could reduce the extrapolation uncertainty. C1 is the gate for C2 (promotion into production), so it blocks the stated end goal of weather-on-CPUE.

**Effort.** Low to medium (analysis, needs the LOO runs on the remaining components).
**Depends on.** The fits converging, which they now do.

### T2.5 Residual single-cell hierarchy funnel: a careful partial reparameterization
**(NEW, building on B1.6's optional note and the reverted B1.7)**

**What.** The pooled model carries `mu_E[g,s] = mu_mu_E[g] + eps_mu_E[g,s] * sigma_mu_E` with G = 1 and S = 1, which is three parameters for one identified intercept. This redundancy is the residual funnel the boat's divergences sit on. A careful partial reparameterization could remove it without the failure mode that reverted B1.7.

**Why it matters.** The divergence-localization table for the boat (`divergence_localization_private_boat_all_gear_Dungeness_Kept.csv`) shows the top divergent-vs-bulk correlate is `sigma_mu_E` at standardized separation 0.36, the largest of any parameter; the structural table shows `sigma_mu_E` only weakly identified (median 0.686, 95 percent CI 0.035 to 3.80), the signature of the redundancy. The naive full collapse (set `mu_E = mu_mu_E`, drop `eps_mu` and `sigma_mu`) was tried as B1.7 and reverted because it forced the overall level to reconcile directly against the high-dimensional daily AR for shore all-gear, producing a long thin ridge that saturated max_treedepth and would not finish in 24 hours. The lesson is that the level term must remain decoupled from the AR. A reparameterization that keeps a single decoupled level absorber (for example, drop `mu_mu` and carry the level entirely in `eps_mu * sigma_mu`, or impose a sum-to-zero constraint) could remove the funnel without re-creating the ridge.

This is now lower urgency than before, because v7.0's gate accepts the boat with its mild funnel on the evidence that the divergences do not move the totals. But removing the funnel honestly is preferable to accepting it via the gate, and it would benefit the gear-resolved track too. Treat it as a design task, not a mechanical edit, given that this exact area has already broken one run. Prototype and test on shore all-gear before any production run.

**Effort.** Medium (design plus a guarded test run).
**Depends on.** v7.0 validated first (so the baseline is stable).

---

## Tier 3: Medium value, accuracy or robustness refinements

### T3.1 Sub-season boundary edge effect
**(backlog B6; critique issue 5)**

**What.** The hard December 1 split fits two fully independent BSS models, so the AR(1) has no anchor across the boundary and the credible interval balloons at the seam (visible as the late-November spike in the daily catch plots). Either document it as a known artifact or move to a single regime-switching model with a `pot_available` indicator that shifts the CPUE baseline while keeping the effort process continuous.

**Why it matters.** It inflates uncertainty near the boundary and prevents information sharing across a date where effort does not actually jump discontinuously. The documentation fix is cheap; the model fix is a meaningful structural change.

**Effort.** Medium for documentation, medium-to-high for the regime-switching model.

### T3.2 Shore CPUE: soak-time versus crabber-hours sensitivity
**(NEW, from critique issue 11 applied to shore)**

**What.** The boat now uses gear-soak-hours (L = 24, `h = gear_time_total`), but shore still uses crabber-hours and day-length (confirmed in the v7.0 Fix 2: shore PE effort is `crabbers_per_gear * day_length`, and the Stan CPUE denominator for shore is `fishing_time_total`). For passive gear the catch is generated by soak time, not presence time. Run a sensitivity that re-expresses shore CPUE on a gear-soak-hours basis and report the harvest movement.

**Why it matters.** The 2026-03-31 critique flagged this structural mismatch for passive gear. It is plausibly small for shore because trap gear dominates shore harvest and traps are checked frequently (short soaks, so crabber-hours and gear-hours track each other), but that is a hypothesis, not a measurement. A sensitivity either retires the concern or reveals a bias in the dominant shore component.

**Effort.** Medium.

### T3.3 Continuous density-dependence in CPUE (beyond the day-type effect)
**(refines critique issue 4, which is otherwise RESOLVED)**

**What.** The model already captures a weekend CPUE effect (`B1_C * w[d]`), and it is real and large: the shore all-gear structural table gives `B1_C = -0.43` (95 percent CI -0.59 to -0.27), an exp(-0.43) of about 0.65, so weekend CPUE is roughly 35 percent lower, while weekend effort is higher (`B1 = 0.51`, about 1.67-fold). That is the crowding-plus-composition signature the critique hypothesized, captured at day-type resolution. The remaining refinement is a continuous density-dependence term, CPUE as a function of same-day effort rather than only a weekend dummy, and a holiday CPUE effect (effort has `B2` for holidays; CPUE has no holiday term).

**Why it matters.** Lower priority because the dominant day-type signal is already modeled. A continuous term would capture within-day-type crowding gradients and could slightly sharpen the summer estimates where effort varies most. Mark critique issue 4 itself as resolved; this is the incremental extension.

**Effort.** Medium.

### T3.4 Rigorous boat-departure timing via interval-level I/E data
**(backlog C6)**

**What.** The boat departure-on-flood-tide test (hypothesis H1 in the covariate module) uses an interview-time proxy (departure equals interview time minus hours fished for completed trips, which assumes the interview is at trip end). Replace it with actual departure events from the I/E surveys, weighting each survey interval by its observed departure count and the tide phase during that interval. This requires exposing interval-level I/E data, since `fetch_ie_data` currently aggregates to daily totals and discards per-interval timestamps, and confirming the raw `ingress_egress.xlsx` carries interval times.

**Why it matters.** Removes a timing assumption and uses the correct data source. Only needed if the departure-on-flood hypothesis is pursued seriously.

**Effort.** Medium.

### T3.5 Commercial/charter census uncertainty quantification
**(gear-resolved planned item 6; affects both tracks)**

**What.** The commercial/charter component is reported as a point with no interval in the port total. The by-mode plot applies a catch-per-vessel CV, but the headline census number carries no formal uncertainty. Add a bootstrap or parametric interval (per-vessel catch sampling variance plus tally coverage).

**Why it matters.** It is the one component with no propagated uncertainty in the port total, so the port-total interval is slightly too narrow. Smaller effect than the BSS components because census coverage is high, but it is a gap a careful reviewer will note. Pairs naturally with T1.4 (the vessel-type split), since both touch the same estimator.

**Effort.** Medium.

---

## Tier 4: Repository, reproducibility, and code-quality cleanup

These do not change the estimate but are required for a publication-quality, reproducible repository.

### T4.1 Reproducibility scaffolding
**(backlog A4)**

**What.** There is no `renv.lock`, no `run_all.R`, and no Makefile (confirmed absent). Add a pinned-package lockfile and a single top-level entry point so the analysis reproduces from a clean clone.

**Why it matters.** A methods paper must reproduce from a clean clone with pinned dependencies. The diagnostics notes already added a fixed seed and session-info capture, which is the right direction; an `renv` lockfile and a one-command runner complete it.

**Effort.** Medium.

### T4.2 Output and repo hygiene
**(backlog A1, A2)**

**What.** `output/` is committed and large: 410 tracked files across eight dated run folders, totaling roughly 30 MB and growing each run. `.Rproj.user/` is still tracked (5 files) despite being gitignored. Untrack `.Rproj.user` (`git rm -r --cached .Rproj.user`) and decide an `output/` policy: gitignore it and keep one canonical reference run, or move runs to GitHub Releases.

**Why it matters.** A publication repo should be reproducible from code, not bloated with every run's artifacts, and IDE state should not be version-controlled. This is the most visible repo-cleanliness issue.

**Effort.** Trivial for `.Rproj.user`; low (one decision) for `output/`.

### T4.3 Global-variable coupling in functions
**(NEW, code quality)**

**What.** Several functions read top-level globals rather than receiving them as arguments. `estimate_comm_charter(dwg, params)` uses `crabbing_holiday_dates` (defined at top level, around Rmd line 226) inside its body; `run_pe` uses the global `catch_groups` (around line 786). `prep_days_crab` does it correctly (it takes `holidays` as a parameter), which highlights the inconsistency.

**Why it matters.** Implicit global dependencies make functions harder to test, reuse, and reason about, and they are a latent source of silent errors if a global is renamed or scoped differently. For publication-quality code, pass these explicitly. Low risk to change, and it improves the auditability the methods paper will be judged on.

**Effort.** Low.

### T4.4 Consolidate the duplicated PE monthly-share blocks
**(carried from v7.0 delivery)**

**What.** Sections 7.8 and 7.8b contain duplicated PE monthly-effort-share code (with a minor `na.rm` difference). v7.0 fixed the day-length bug identically in both rather than refactoring, to keep the change small. Consolidate into one helper.

**Why it matters.** Duplicated logic is a drift risk; the two copies must be kept in sync by hand, as the day-length bug itself demonstrated (it existed in both). One helper removes that risk.

**Effort.** Low.
**Depends on.** Best after the v7.0 run validates the current behavior.

### T4.5 Gear-resolved track maintenance items
**(backlog A3, D1, D2, D4)**

**What.** Four gear-resolved or module items: (A3) the gear-resolved documentation and Section 8.3 say `max_divergences` defaults to 0, but the code uses 5 (confirmed: `max_divergences = 5` in the gear-resolved Rmd, while the doc comment says default 0), a doc-code mismatch to reconcile. (D4) the gear-resolved Stan model does not carry the non-centered `omega_0` reparameterization (confirmed: no `omega_E_0_raw` in `crab_bss_gear_resolved.stan`), so it almost certainly retains the centered funnel that B1.3 fixed in the pooled track, and its last run had every component fall back to PE; port B1.3 once the pooled version is validated. (D1) the covariate module recomputes a baseline fit the main pooled run already produced; resolve via a save/load handoff or by merging the two Rmds. (D2) cross-check `crab_bss_pooled_weather_adjusted.stan` against `crab_bss_pooled.stan` at the `log_lik` level and confirm the `p_I_shore` versus `p_TI` naming so the baseline-vs-covariate LOO is valid.

**Why it matters.** These keep the secondary tracks honest and the covariate comparison valid. D4 in particular is the same convergence fix the pooled track benefited from; the gear-resolved track is currently all-PE and would likely benefit the same way. Lower priority than the pooled-track items above because the pooled track is the publication focus, but D4 and D2 gate the gear-resolved and covariate results respectively.

**Effort.** Low for A3 and D4 (mechanical); medium for D1; low-medium for D2.

### T4.6 Field-protocol items for 2025-26 (not code)
**(backlog D3 and the critique's sampling recommendations; critique issues 3, 10, 11)**

**What.** Several improvements require a field-protocol change for 2025-26, not a code change, and should feed the sampling-design plan rather than this backlog: jetty and beach direct effort counts (critique issue 3, backlog D3); minimum coverage targets by mode with a specific emphasis on the summer months that dominate the estimate (directly motivated by T1.1); interview-time versus effort-count-time representativeness checks (critique issue 10); and effort-count and incomplete-trip protocols. The single highest-value sampling change implied by this review is increasing summer interview and effort coverage, because the harvest estimate's largest uncertainty and its BSS-vs-PE divergence both live in July through September.

**Why it matters.** Code refinements cannot manufacture information that was never sampled. The durable fix for the summer-extrapolation problem in T1.1, and for the boat component generally, is more informative summer effort and interview data. This belongs in the 2025-26 sampling design, flagged here so the link to the estimate's weakest point is explicit.

**Effort.** Out of code scope; design-dependent.

---

## Appendix A: Status reconciliation of prior IDs and critique issues

| Prior ID / Issue | Description | Status after v7.0 and this review | Mapped to |
|---|---|---|---|
| A1 | Untrack `.Rproj.user/` | OPEN (5 files still tracked) | T4.2 |
| A2 | `output/` policy | OPEN (410 files tracked) | T4.2 |
| A3 | Gear-resolved `max_divergences` doc/code mismatch | OPEN (doc says 0, code uses 5) | T4.5 |
| A4 | Reproducibility scaffolding (renv, runner) | OPEN (none present) | T4.1 |
| A5 | Data-shipping decision | OPEN | T4.1 area |
| B1 | Make BSS converge or characterize PE fallback | DONE pending run (boat converges at 0.6 percent; v7.0 gate moves it to BSS; only boat ring-net remains on PE, data-limited) | T1.1 (validation) |
| B1.1 to B1.3 | Diagnose, coarsen boat AR, non-center omega_0 | DONE | |
| B1.4 | Characterize the fallback | SATISFIED for the one remaining PE component | |
| B1.5 | Marginalize effort overdispersion to NB | DONE (v6.7) | |
| B1.6 | sigma_IE proper prior; distortion-aware gate | DONE (v6.8); gate now superseded by v7.0 scale-aware gate | |
| B2 | Posterior predictive checks | UNBLOCKED by v7.0 PPC fix; now Tier 1 | T1.2 |
| B3 | Prior sensitivity table | OPEN | T1.3 |
| B4 | External validation | OPEN | T1.1 |
| B5 | Propagate pi_gear uncertainty (pooled) | OPEN | T2.2 |
| B6 | Sub-season boundary edge effect | OPEN | T3.1 |
| B7 | Incomplete-trip handling (pooled) | OPEN | T2.3 |
| C1 | Run weather module, read LOO | PARTIAL (GAM screen done; LOO only for shore ring-net, marginal) | T2.4 |
| C2 | Fold supported covariates into production | OPEN (blocked by C1) | T2.4 |
| C3 | Reproducible environmental inputs | OPEN | T4.1 area |
| C4 | Non-linear covariate effects (splines) | OPEN | T2.4 follow-on |
| C5 | Weather in gear-resolved CPUE | OPEN (low priority) | |
| C6 | Interval-level I/E for departure timing | OPEN | T3.4 |
| D1 | Redundant baseline refit in covariate module | OPEN | T4.5 |
| D2 | Module log_lik reconciliation | OPEN | T4.5 |
| D3 | Jetty/beach effort counts | OPEN (sampling, not code) | T4.6 |
| D4 | Port non-centered omega_0 to gear-resolved | OPEN (confirmed absent in gear-resolved Stan) | T4.5 |
| Critique 1 | Effort-to-crabber and day-length expansion | PARTIAL (boats use gear-hours; shore still snapshot times day-length, partly mitigated by I/E-informed L_effective) | T1.1, T1.3, T3.2 |
| Critique 2 | R_G ratio sensitivity | OPEN (data-driven prior only partial) | T1.3 |
| Critique 3 | Jetty/beach lumping | OPEN (sampling) | T4.6 |
| Critique 4 | No day-type CPUE effect | RESOLVED (B1_C implemented; weekend CPUE effect is real at -0.43); continuous density-dependence remains | T3.3 |
| Critique 5 | Sub-season hard boundary | OPEN | T3.1 |
| Critique 6 | Commercial/charter single per-vessel mean | OPEN (confirmed; boat_type and tally split support the fix) | T1.4 |
| Critique 7 | Gear-resolved priority hierarchy | RESOLVED in gear-resolved track (weighted pseudo-likelihood, pi_gear empirical) | |
| Critique 8 | Shared phi_C_gear | OPEN (gear-resolved, low priority) | |
| Critique 9 | Zero-inflation / ZINB | OPEN but now testable via PPC | T1.2 |
| Critique 10 | Interview representativeness | OPEN (incomplete-trip filter partly addresses; representativeness check is sampling) | T2.3, T4.6 |
| Critique 11 | Spatial heterogeneity, weather, soak time | PARTIAL (weather module exists; boat soak-time done; shore soak-time and spatial open) | T2.4, T3.2 |

---

## Appendix B: Evidence cited

All figures above are read from the committed repository, not from memory:

- Monthly harvest concentration and credible intervals: `output/20260621/pooled-CPUE/monthly_estimates.csv`.
- PE-vs-BSS divergences (shore +50.1 percent, boat -11.9 percent) and method labels: `output/20260621/pooled-CPUE/pe_vs_bss_comparison.csv`.
- Port total (100,448 BSS expected catch) and effort: `output/20260621/pooled-CPUE/port_total_Dungeness_Kept.csv`.
- Boat convergence (59 divergent, R-hat, n_eff) and the gate verdict: `output/20260621/pooled-CPUE/convergence_report.csv`.
- Boat funnel on `sigma_mu_E` (smd 0.36): `output/20260621/pooled-CPUE/divergence_localization_private_boat_all_gear_Dungeness_Kept.csv`.
- Weekend effort and CPUE effects (B1 = 0.51, B1_C = -0.43) and `sigma_mu_E` identification: `output/20260621/pooled-CPUE/structural_params_shore_all_gear_Dungeness_Kept.csv`.
- Absence of PPC output files pre-v7.0: file listing of `output/20260621/pooled-CPUE/`.
- Covariate LOO marginal result and selected features: `output/20260621/pooled-CPUE-covariates/loo_comparison_shore_ring_net_only_Dungeness_Kept.csv` and `layer_b_selected_features.csv`.
- Commercial/charter estimator and the single per-vessel mean: `BSS-GH-pooled-CPUE-model.Rmd`, `estimate_comm_charter`.
- Vessel-type support for the split: `input_files/interview_combined.csv` (`boat_type`) and `input_files/wes_commercial_tally.csv` (split tally columns).
- Catch likelihood (no zero inflation) and the single-cell hierarchy: `stan_models/crab_bss_pooled.stan`.
- Repository hygiene (tracked `.Rproj.user`, tracked `output/`, no renv/runner) and gear-resolved funnel absence: `git ls-files` and `stan_models/crab_bss_gear_resolved.stan`.
