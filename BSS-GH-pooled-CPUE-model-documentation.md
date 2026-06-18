# Recreational Dungeness Crab Harvest Estimation --- Grays Harbor

## Pooled CPUE Model: Technical Documentation

**Author:** Matthew George, Ph.D.\
**Contact:** [matthew.george\@dfw.wa.gov](mailto:matthew.george@dfw.wa.gov){.email}\
**Agency:** Washington Department of Fish and Wildlife (WDFW)\
**Status:** Operational --- Annual estimation framework

------------------------------------------------------------------------

## 1. Summary for Decision-Makers

This framework estimates the total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area. It combines four types of field observations --- gear counts at the docks, trailer counts at the boat launch, dockside crabber interviews, and ingress/egress surveys --- with a statistical model that fills in the days when no sampling occurred.

**What this model produces:**

-   A total Dungeness crab harvest estimate for the port with a 95% credible interval.
-   Monthly harvest trends showing when crabbing pressure peaks and how it changes through the season.
-   Breakdowns by crabbing mode (shore, private boat, commercial/charter) and gear type (pot, ring net, trap, snare).
-   A CPUE day-type effect (B1_C) quantifying whether weekend catch rates differ from weekday catch rates.
-   Daily posterior estimates of effective day length (L_effective) when I/E data is available.

**How confident are we?** The framework runs two independent estimation methods --- a simple average-based approach (PE) and a Bayesian time-series model (BSS) --- then compares them. When the two methods agree closely and the BSS model converges properly, confidence is high. The output includes formal convergence diagnostics and a comparison table so reviewers can assess reliability.

------------------------------------------------------------------------

## 2. Understanding the Two Estimation Methods

### 2.1 Point Estimator (PE) --- The Simple Average

The PE method computes the average daily harvest for each stat-week × day-type stratum, then multiplies by the total number of days in that stratum to estimate harvest on unsampled days (Pollock et al. 1994; Hahn et al. 2000).

**Strengths:** Easy to understand. Transparent. No modeling assumptions beyond representativeness within strata.

**Weaknesses:** Cannot fill temporal gaps with zero samples. No uncertainty bounds. Treats each stratum independently.

### 2.2 Bayesian State-Space Model (BSS)

The BSS method fits a smooth curve through the daily effort and catch rate data using an autoregressive (AR(1)) state-space model, then uses that curve to estimate every day in the season --- including days with no field sampling. The approach follows the Bayesian creel survey framework developed by Conn (2002) and extended by Staton et al. (2017), where latent daily effort and CPUE processes evolve as first-order autoregressive time series with observation error.

**Adaptive temporal resolution:** The AR(1) process resolution is selected automatically based on effort data density for each population × sub-season fit:

| Resolution  | Condition                                | Rationale                                                                                                                              |
|-------------|------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| **Daily**   | ≥25% of days sampled AND ≥20 effort days | Dense data supports day-level smoothing; proper uncertainty scaling by temporal distance from nearest observation (Staton et al. 2017) |
| **Weekly**  | ≥1.5 effort obs per week AND ≥3 weeks    | Moderate data; weekly states smooth over 3--5 day gaps without under-identifying the AR dynamics                                       |
| **Monthly** | Fallback for sparse data                 | Conservative; few AR parameters to estimate, robust with limited observations                                                          |

This adaptive approach resolves a tension identified in the creel survey literature: daily AR processes provide the smoothest interpolation and most honest uncertainty quantification (Conn 2002; Sullivan 2003), but require sufficient observation density to identify the autocorrelation parameters. When data is too sparse --- as with boat trailer counts during winter months --- a daily AR with D ≈ 60 latent states but only 10 observations creates a poorly-identified posterior with difficult HMC geometry. The adaptive rule applies the finest resolution the data can support while falling back gracefully when it cannot.

On days with ingress/egress surveys, the BSS receives a direct crabber-hours observation in addition to the gear count. The effective day length on those days is informed by both the I/E observation and a regression prior, with the posterior narrowing around the I/E-anchored value.

**Strengths:** Fills temporal gaps with proper uncertainty scaling. Accounts for temporal autocorrelation. Produces rigorous uncertainty bounds. I/E anchor points constrain the effort trajectory.

**Weaknesses:** More complex. Requires 30--120 minutes per model fit. Must be checked for convergence.

### 2.3 Combined Best Estimate

The framework checks convergence using rank-normalized split-R-hat and bulk effective sample size (Vehtari et al. 2021) together with divergent transitions (Betancourt 2017), evaluated on C_expected_sum and E_sum. When all criteria pass, the BSS expected catch estimate is preferred; otherwise the PE is used. Section 8 gives the exact thresholds.

The primary harvest estimate uses C_expected (the posterior expected catch, E[C\|data]) rather than the Poisson predictive draw, following the distinction between estimation and prediction in hierarchical models (Gelman et al. 2013, Ch. 7). The predictive distribution (C_sum) is reported separately for computing prediction intervals.

------------------------------------------------------------------------

## 3. The Recreational Crab Fishery at Grays Harbor

### 3.1 Fishery Overview

The recreational Dungeness crab (*Metacarcinus magister*) fishery at Westport and the greater Grays Harbor area is one of the highest-volume recreational crabbing operations on the Washington coast. Recreational crabbers use four primary gear types: crab pots (highest CPUE), ring nets, foldable/star traps, and snares. WDFW regulations restrict pot use to December through September, creating a structural break in both effort and catch rates.

Commercial Dungeness crab vessels also participate in the recreational fishery before the commercial season opens, crabbing recreationally under the same daily limits as private boats. Their harvest is tracked separately through a vessel tally system at the marina.

### 3.2 Study Area

Westport is a small coastal town on the south side of the Grays Harbor estuary. Recreational crabbing occurs from multiple access points: public docks (Floats 17--21), a jetty, beaches, a public boat launch, and a commercial marina.

------------------------------------------------------------------------

## 4. Data Sources and Field Collection

### 4.1 Effort Counts

Instantaneous point-in-time counts of gear or trailers conducted by field surveyors. The primary input for estimating total crabbing activity.

| Site                                  | What Is Counted       | Role in Model                 |
|---------------------------------------|-----------------------|-------------------------------|
| Westport Docks Float 20 + Float 17-21 | Crab gear in water    | Shore effort indicator        |
| Westport Boat Launch                  | Boat trailers at ramp | Private boat effort indicator |

### 4.2 Crabber Interviews

Dockside interviews recording trip-level information: group size, gear deployed, gear type, hours fished, crab kept, trip status. The CPUE denominator is crabber-hours for shore interviews and gear-hours for boat interviews (Pollock et al. 1994).

### 4.3 Commercial/Charter Vessel Tally

Daily count of commercial and charter boats at Westport Marina during the recreational pre-season period. Combined with mean catch per vessel from interviews for stratified expansion.

### 4.4 Ingress/Egress (I/E) Surveys

All-day surveys recording crabber arrivals and departures every 15 minutes. Daily crabber-hours equals the area under the crabber-presence curve (sum of present × 0.25 hr). This provides a direct measurement of daily effort that bypasses the gear-count-to-effort conversion chain. The approach follows the bus-route and access-point survey methods described by Robson (1991) and Pollock et al. (1997).

------------------------------------------------------------------------

## 5. Key Design Decisions

### 5.1 Sub-Seasons

**Ring-net only** (Sep 16 -- Nov 30): Pots prohibited. **All-gear** (Dec 1 -- Sep 15): Pots allowed. Each sub-season is estimated independently.

### 5.2 Day Length and Effort Units

**Shore crabbers --- L_effective regression:** Rather than using civil twilight (9--16 hours) as a day-length proxy, the framework estimates an empirical "effective day length" from I/E survey data. L_effective captures the peaked activity curve at the docks --- crabbers rotate through rather than occupying the dock all day. Analysis of I/E data shows L_effective averaging 3.5--5.5 hours, substantially shorter than civil twilight.

The L_effective model fits a regression of log(L_effective) on day-of-year (quadratic) and day type:

```         
log(L_effective) = β₀ + β₁ × yday + β₂ × yday² + β₃ × weekend + ε
```

The quadratic captures the seasonal arc in effective day length. For each day in the estimation period, the regression produces a predicted median (L_mu) and total prediction uncertainty (L_sigma) on the log scale. These enter the Stan model as a lognormal prior on the day-length parameter, propagating L_effective uncertainty into the catch estimate (Pollock et al. 1994 §4.3; Hartill et al. 2012).

**Private boats:** Day length is fixed at 24 hours because boat gear (primarily pots) soaks continuously. The effort unit is gear-hours.

### 5.3 Day Type

Weekday (Mon--Thu), weekend (Fri--Sun), or holiday. Separate B1 (weekend) and B2 (holiday) effort effects in the BSS. The CPUE process includes B1_C, a weekend CPUE effect that allows catch rate to differ by day type --- an extension motivated by evidence that weekend crabber populations include a higher proportion of less-experienced participants (Thomson 1991; Pollock et al. 1997).

------------------------------------------------------------------------

## 6. Population Components

### 6.1 Shore Crabbers (Dock + Jetty + Beach)

Effort indicator: gear counts at the docks. Conversion: (Gear counted ÷ R_G) × L_effective = crabber-hours. On I/E days, crabber-hours are observed directly.

### 6.2 Private Boat Crabbers

Effort indicator: trailer counts at boat launches. Conversion: trailer_count × gear_per_group × 24 = gear-hours. CPUE uses gear-hours as the denominator.

### 6.3 Commercial/Charter Vessels

Estimated via day-type-stratified census expansion from the vessel tally.

------------------------------------------------------------------------

## 7. The Pooled CPUE Model (Technical)

### 7.1 Effort Process

```         
log(lambda_E[d]) = mu_E + omega_E[period(d)] + B1 × w[d] + B2 × holiday[d]
```

The temporal deviation omega_E evolves as an AR(1) process:

```         
omega_E[p] = phi_E × omega_E[p-1] + sigma_eps_E × epsilon[p-1]
```

where `period(d)` maps day d to its AR period index. When AR resolution is daily, period(d) = d and P_n = D. When weekly or monthly, period(d) maps to the corresponding week or month index, and P_n equals the number of weeks or months. The innovations epsilon are standard normal (non-centered parameterization for efficient HMC sampling; Papaspiliopoulos et al. 2007).

The stationary initial state prior is omega_E_0 \~ Normal(0, sigma_eps_E / sqrt(1 - phi_E²)), ensuring the AR process starts from its stationary distribution rather than requiring a burn-in period (Harvey 1989, Ch. 3).

### 7.2 CPUE Process

```         
log(lambda_C[d]) = mu_C + omega_C[period(d)] + B1_C × w[d]
```

B1_C allows weekend CPUE to differ from weekday CPUE. This is motivated by the observation that weekend/holiday crabber populations at tourist-accessible ports include more novice participants with potentially different catch rates (Thomson 1991; Pollock et al. 1997). Empirical estimates from Grays Harbor show B1_C ≈ -0.25 to -0.30 for shore crabbers (weekenders catch 21--26% fewer crab per crabber-hour than weekday regulars), consistent with the novice-dilution hypothesis.

### 7.3 Observation Models

-   **Gear counts (shore):** `Gear_I ~ Poisson(lambda_E[d] × eps_E_H × R_G)`
-   **Trailer counts (boats):** `T_I ~ Poisson(lambda_E[d] × eps_E_H × R_T)`
-   **Interview catch:** `c ~ NegBin(lambda_C[d] × h, r_C)` where h = crabber-hours (shore) or gear-hours (boats)
-   **I/E crabber-hours:** `IE_crabber_hours ~ Lognormal(log(lambda_E[d] × L[d]), sigma_IE)`

The negative binomial catch likelihood accommodates overdispersion in individual catch data, which is common in recreational fisheries where trip-level catch rates are highly variable (Maunder & Punt 2004).

### 7.4 Sparse Overdispersion

Within-day overdispersion parameters `eps_E_H_obs` are allocated only for actual effort observations (not all D × count_sequence slots), eliminating 64--77% of effort parameters. Each eps_E_H_obs \~ Gamma(r_E, r_E) with mean 1, following the Gamma-Poisson mixture parameterization of the negative binomial (Hilbe 2011).

### 7.5 I/E Integration

On I/E survey days, observed crabber-hours enter as a direct lognormal observation of lambda_E × L. This provides a second, independent constraint on the latent effort state that bypasses R_G and day-length assumptions. The dual-observation design (gear count + I/E on paired days) calibrates the gear-count pathway against the I/E ground truth, analogous to paired census-index designs in roving creel surveys (Robson 1991; Pollock et al. 1994).

When no I/E data is available (IE_n = 0), the I/E likelihood contributes nothing and the model reverts to gear-count-only estimation with no change in behavior.

### 7.6 L_effective as a Parameter

When `estimate_L = 1` (shore fits), effective day length L[d] is a parameter with non-centered lognormal prior:

```         
L[d] = L_mu[d] × exp(L_sigma[d] × L_raw[d]),    L_raw ~ Normal(0, 1)
```

where L_mu and L_sigma come from the I/E regression. On I/E days, L[d] is additionally constrained by the I/E likelihood, producing a tighter posterior. On non-I/E days, L[d] is informed only by the regression prior, and its uncertainty propagates into the effort and catch estimates.

### 7.7 Key Parameters

| Parameter    | Description                     | Prior                                              |
|--------------|---------------------------------|----------------------------------------------------|
| B1           | Weekend effort multiplier (log) | Normal(0, 1)                                       |
| B2           | Holiday effort multiplier (log) | Normal(0, 1)                                       |
| B1_C         | Weekend CPUE effect (log)       | Normal(0, 1)                                       |
| R_G          | Gear per crabber                | Lognormal(log(R_G_empirical), 0.3) --- data-driven |
| R_T          | Trailers per boat group         | Beta(5, 1)                                         |
| phi_E, phi_C | AR(1) autocorrelation           | Beta(2,2) rescaled [-1,1]                          |
| r_E, r_C     | Overdispersion                  | Half-Cauchy(0, 1)                                  |
| sigma_IE     | I/E measurement error (log)     | Exponential(5)                                     |
| L[d]         | Effective day length (shore)    | Lognormal from regression                          |

**Prior rationale notes:**

-   **R_G**: Prior center computed from the empirical gear-per-crabber ratio in interview data for the relevant population × sub-season, eliminating prior-posterior conflict.
-   **R_T**: Beta(5, 1) concentrates mass near 1 (most boat groups bring one trailer), replacing the uninformative Beta(0.5, 0.5).
-   **Half-Cauchy(0, 1)**: Weakly informative variance priors following Gelman (2006) and the Stan development team recommendations. Scale of 1 is appropriate for variance components in a well-characterized recreational fishery.

### 7.8 Generated Quantities

The model reports two catch quantities:

-   **C_expected[d]** = lambda_E[d] × L[d] × lambda_C[d] --- the posterior expected daily catch with no Poisson sampling noise. This is E[C\|data], the natural quantity for harvest estimation.
-   **C[d]** = Poisson_rng(C_expected[d]) --- a predictive draw including Poisson sampling variability.

For seasonal totals, Poisson noise largely averages out (CLT), so the two distributions are similar. For daily or monthly breakdowns, the difference can be substantial.

**Stan model file:** `crab_bss_pooled.stan`

------------------------------------------------------------------------

## 8. Convergence and Model Selection

For each BSS fit the framework monitors four diagnostics, reported per fit in `convergence_report.csv`: rank-normalized split-R-hat and bulk effective sample size for `C_expected_sum` and `E_sum` (Vehtari et al. 2021), the number of divergent transitions, and the percentage of post-warmup iterations that saturate `max_treedepth`.

A fit **passes**, and its BSS estimate is preferred, when all of the following hold; otherwise the PE estimate is used for that component:

-   R-hat \< 1.01 for `C_expected_sum` and `E_sum`.
-   n_eff \> 400 for `C_expected_sum` and `E_sum`.
-   Divergent transitions ≤ `max_divergences` (default 5).

Divergent transitions are part of the gate because they can bias the posterior even when R-hat and n_eff look satisfactory; a sampler that cannot integrate the Hamiltonian trajectory accurately is not exploring the target distribution, regardless of how well the chains agree (Betancourt 2017). This is the same standard applied in the gear-resolved track, so both models now use one convergence gate. Treedepth saturation above 5% raises a warning rather than a hard failure: it signals truncated trajectories that reduce effective sample size, and is addressed by raising `max_treedepth` for the affected fit. Per-fit `max_treedepth` and `adapt_delta` overrides are set for the shore all-gear fit (14, 0.95), the boat all-gear fit (13, 0.99), and the ring-net fits (12, 0.95), with all other fits using the defaults (10, 0.9). The report also records the AR resolution used for each fit, so reviewers can judge whether the selected resolution was appropriate.

> Note on the R-hat threshold: the gate uses R-hat \< 1.01, following Vehtari et al. (2021), who developed the rank-normalized R-hat and recommend this threshold (the ESS \> 400 criterion is from the same source). This was tightened from a historical 1.05 in v6.4; the change does not affect the 2024-25 results, because every passing fit has R-hat near 1.00 and every failing fit fails on divergent transitions or n_eff. The gear-resolved track uses the same 1.01 threshold.

------------------------------------------------------------------------

## 9. Output Files

Each run produces output in `output/YYYYMMDD/`:

| File                            | Contents                                               |
|---------------------------------|--------------------------------------------------------|
| `pe_port_summary.csv`           | PE estimates by component and port total               |
| `port_total_Dungeness_Kept.csv` | Combined PE + BSS port total (expected and predictive) |
| `monthly_estimates.csv`         | Monthly catch and effort with credible intervals       |
| `catch_by_mode.csv`             | Catch by crabbing mode (shore, boat, commercial)       |
| `catch_by_gear_type.csv`        | Catch by gear type (proportional allocation)           |
| `convergence_report.csv`        | Diagnostics per BSS fit including AR resolution        |
| `effort_cpue_multipliers.csv`   | B1, B2, B1_C posteriors                                |
| `expansion_ratios.csv`          | R_G, R_T posteriors                                    |
| `ie_analysis.csv`               | I/E validation with f_temporal                         |
| `bss_L_effective_{label}.csv`   | Daily L posteriors (prior, median, 95% CI)             |
| `L_effective_ie_detail.csv`     | Per-I/E-day regression predictions vs observed         |
| `pe_vs_bss_comparison.csv`     | PE vs BSS effort and catch by component, with the selected method |

------------------------------------------------------------------------

## 10. Limitations and Future Directions

-   Limited I/E coverage; expand to \~40 days per season for better L_effective regression (Pollock et al. 1997 recommend ≥3 I/E days per month × day-type stratum).
-   L_effective regression uses a quadratic in day-of-year; with more data, a GAM could capture non-monotonic patterns.
-   No weather covariates in effort or CPUE processes; NOAA buoy data could improve prediction on unsampled days (Conn 2002 included temperature in the Kenai River creel model).
-   Gear-type breakdowns are approximate (proportional allocation, not modeled separately).
-   No jetty effort counts. Beach crabbing unmeasured.
-   The B1_C effect is constant across the season; a time-varying weekend CPUE effect may be warranted if tourist composition shifts seasonally.
-   The adaptive AR selection is rule-based; a formal model comparison (LOO-CV or WAIC; Vehtari et al. 2017) could provide principled resolution selection.
-   The private boat all-gear BSS fit is prone to non-convergence. The trailer-count effort series is sparse and the latent effort process is weakly identified, producing many divergent transitions (5998 in the 2026-04-08 run, with no treedepth saturation, which points to a step-size rather than a trajectory-length problem). v6.2 adds dedicated sampler tuning (adapt_delta 0.99, max_treedepth 13, more iterations) as an attempt to achieve convergence; when the fit still fails the gate, the boat estimate uses PE. If tuning does not succeed, the likely structural levers are an informative prior on the boat effort process or forcing a coarser AR resolution for the boat fit (the adaptive rule selected daily AR despite the sparse effort series).

------------------------------------------------------------------------

## 11. Glossary

| Term              | Meaning                                                             |
|-------------------|---------------------------------------------------------------------|
| **BSS**           | Bayesian State-Space model                                          |
| **PE**            | Point Estimator                                                     |
| **CPUE**          | Catch Per Unit Effort                                               |
| **AR(1)**         | First-order autoregressive process                                  |
| **P_n**           | Number of AR periods (= D for daily, fewer for weekly/monthly)      |
| **period(d)**     | Mapping from day d to its AR period index                           |
| **R_G**           | Gear-per-crabber ratio                                              |
| **B1_C**          | Weekend CPUE multiplier; exp(B1_C) = weekend/weekday CPUE ratio     |
| **C_expected**    | Expected daily catch (no Poisson noise)                             |
| **L_effective**   | Effective day length: I/E crabber-hours ÷ peak crabbers present     |
| **L_mu, L_sigma** | Regression-predicted median and uncertainty for L_effective         |
| **I/E**           | Ingress/Egress survey                                               |
| **f_temporal**    | Temporal correction factor: I/E crabber-hours ÷ gear-count estimate |
| **sigma_IE**      | I/E measurement error on log scale                                  |

------------------------------------------------------------------------

## 12. Reproducibility

1.  Clone the repository. Place input files in `input_files/`.
2.  Place `crab_bss_pooled.stan` in `stan_models/`.
3.  Open `BSS-GH-pooled-CPUE-model.Rmd`, update parameters, run.

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here, readxl.\
**Expected runtime:** 3--6 hours on a 4-core machine (varies with AR resolution and sub-season length).

------------------------------------------------------------------------

## 13. References

Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo. *arXiv preprint* arXiv:1701.02434.

Conn, P.B. (2002). Bayesian methods for estimating recreational angler effort, catch rates, and total catch using creel survey data. Ph.D. Dissertation, University of Wisconsin-Madison.

Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models. *Bayesian Analysis*, 1(3), 515--534.

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., & Rubin, D.B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

Hahn, P.K.J., Brooks, L., & Hartill, B.W. (2000). Strategies and procedures for estimating catch and effort in freshwater fisheries. *In:* Inland Fisheries Management in North America (2nd ed.), American Fisheries Society.

Hartill, B.W., Cryer, M., Lyle, J.M., Rees, E.B., Ryan, K.L., Steffe, A.S., Taylor, S.M., West, L., & Wise, B.S. (2012). Scale- and context-dependent selection of recreational harvest estimation methods: the Australasian experience. *North American Journal of Fisheries Management*, 32(1), 109--123.

Harvey, A.C. (1989). *Forecasting, Structural Time Series Models and the Kalman Filter*. Cambridge University Press.

Hilbe, J.M. (2011). *Negative Binomial Regression* (2nd ed.). Cambridge University Press.

Maunder, M.N. & Punt, A.E. (2004). Standardizing catch and effort data: a review of recent approaches. *Fisheries Research*, 70(2--3), 141--159.

Papaspiliopoulos, O., Roberts, G.O., & Sköld, M. (2007). A general framework for the parametrization of hierarchical models. *Statistical Science*, 22(1), 59--73.

Pollock, K.H., Jones, C.M., & Brown, T.L. (1994). *Angler Survey Methods and Their Applications in Fisheries Management*. American Fisheries Society Special Publication 25.

Pollock, K.H., Hoenig, J.M., Jones, C.M., Robson, D.S., & Greene, C.J. (1997). Catch rate estimation for roving and access point surveys. *North American Journal of Fisheries Management*, 17(1), 11--19.

Robson, D.S. (1991). The roving creel survey. *American Fisheries Society Symposium*, 12, 137--148.

Staton, B.A., Catalano, M.J., Connors, B.M., Coggins, L.G., Jones, M.L., Walters, C.J., Fleischman, S.J., & Beardsall, J.W. (2017). Evaluation of methods for spawner-recruit analysis in mixed-stock Pacific salmon fisheries. *Canadian Journal of Fisheries and Aquatic Sciences*, 74(7), 1108--1122.

Sullivan, M.G. (2003). Active management of walleye fisheries in Alberta: dilemmas of managing recovering fisheries. *North American Journal of Fisheries Management*, 23(4), 1343--1358.

Thomson, C.J. (1991). Effects of the avidity bias on survey estimates of fishing effort and economic value. *American Fisheries Society Symposium*, 12, 356--366.

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413--1432.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.C. (2021). Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667--718.

------------------------------------------------------------------------

## 14. Version History

Versions continue the shared milestone sequence used in `README.md` (which documents v1--v5). The pooled and gear-resolved tracks have interleaved since v5; the gear-resolved documentation maintains its own v5.x change log. If a different numbering scheme is preferred, these entries can be renumbered.

### v6.4 (2026-06-17), R-hat threshold tightened to 1.01

-   Tightened the convergence gate from R-hat < 1.05 to R-hat < 1.01 for `C_expected_sum` and `E_sum`, matching Vehtari et al. (2021), the source of the rank-normalized R-hat and the ESS > 400 criterion. The same change was applied to the gear-resolved track (v5.4) so both gates use one threshold.
-   Outcome-neutral for the 2024-25 season: every passing fit has R-hat near 1.00, and every failing fit fails on divergent transitions or n_eff, not on an R-hat between 1.01 and 1.05.
-   Files changed: `BSS-GH-pooled-CPUE-model.Rmd` (convergence gate); this documentation (Sections 8, 14).

### v6.3 (2026-06-17), Documentation corrections (no code change)

-   Corrected the Section 9 output-file listing: replaced `daily_combined_estimate.csv` (produced by the gear-resolved model, not the pooled pipeline) with `pe_vs_bss_comparison.csv`, which the pooled run actually writes.
-   Corrected the Vehtari citation in Sections 2.2 and 8. The rank-normalized R-hat and the ESS > 400 criterion are attributed to Vehtari et al. (2021); the note now states accurately that Vehtari recommend R-hat < 1.01, while the gate retains its operational R-hat < 1.05 pending a cross-track threshold decision. Section 2.2 also now reflects the divergence criterion added in v6.1.
-   No change to the Stan model, the R pipeline, or the convergence gate logic.
-   Files changed: this documentation (Sections 2.2, 8, 9, 14).

### v6.2 (2026-06-17), Boat all-gear sampler tuning

-   **Dedicated sampler settings for the private boat all-gear fit.** Added `bss_iter_boat_allgear` (5000), `bss_warmup_boat_allgear` (2500), `bss_treedepth_boat_allgear` (13), and `bss_delta_boat_allgear` (0.99), applied through a new boat all-gear branch in the per-fit Stan control. The boat all-gear fit failed convergence in the 2026-04-08 run (5998 divergent transitions, R-hat 1.074, n_eff 76) under the defaults (adapt_delta 0.9, max_treedepth 10) and fell back to PE. Because `treedepth_pct` was 0, the divergences came from step size rather than truncated trajectories, so `adapt_delta` (raised to 0.99) is the primary lever and `max_treedepth` (raised to 13) is a buffer against the longer trajectories a higher `adapt_delta` produces.
-   **Per-fit tuning consolidated.** `iter`, `warmup`, `max_treedepth`, and `adapt_delta` are now set together in one `if/else` block (matching the gear-resolved pattern), rather than `iter`/`warmup` in separate ternaries. Behavior is unchanged for the shore all-gear and ring-net fits.
-   **This is an attempt, not a guarantee.** If the boat all-gear fit still exceeds `max_divergences` or fails R-hat / n_eff under the new settings, it falls back to PE exactly as before. n_eff = 76 was far below the 400 threshold, and additional iterations help only if the geometry fix removes the autocorrelation; if n_eff (not divergences) remains the binding failure, the structural levers in Section 10 apply. Expect a longer runtime for this fit (higher adapt_delta and more iterations mean more leapfrog steps per iteration).
-   Files changed: `BSS-GH-pooled-CPUE-model.Rmd` (parameters, per-fit Stan control); this documentation (Sections 8, 10, 14).

### v6.1 (2026-06-17), Divergence-aware convergence gate

-   **Convergence gate now includes divergent transitions.** A BSS fit passes only when R-hat, n_eff, **and** divergent transitions (≤ `max_divergences`, default 5) all pass; otherwise it falls back to PE. Previously the gate checked only R-hat and n_eff, so the shore all-gear fit in the 2026-04-08 run was reported as a clean BSS pass despite 842 divergent transitions and 91.6% treedepth saturation. Divergences can bias the posterior even when R-hat and n_eff are satisfactory (Betancourt 2017). This matches the gear-resolved v5.2 gate, so both tracks now apply one standard.
-   **Treedepth and divergence warnings.** The run log now warns when treedepth saturation exceeds 5% or when any divergent transitions are detected, even if the fit still passes.
-   **Per-fit `max_treedepth` / `adapt_delta` overrides.** Added so that fits needing deeper trajectories actually use them: shore all-gear (`max_treedepth` 14, `adapt_delta` 0.95) and ring-net (`max_treedepth` 12, `adapt_delta` 0.95). Any fit without an override uses the defaults (`max_treedepth` 10, `adapt_delta` 0.9). This reduces the divergences and treedepth truncation that previously affected the dense shore all-gear geometry.
-   **Expected effect on the 2024-25 estimate.** On re-run, the shore all-gear fit will either converge cleanly under the deeper-tree / higher-delta settings (preferred outcome) or fail the divergence gate and fall back to PE. Either way the reported BSS estimate will no longer be a divergence-contaminated pass. The private boat all-gear fit already failed on R-hat and continues to use PE; its dedicated tuning is handled as a separate change (boat refit).
-   Files changed: `BSS-GH-pooled-CPUE-model.Rmd` (parameters, per-fit Stan control, convergence report); this documentation (Sections 8, 13, 14).

### v6.0 (pre-change baseline, current GitHub state), Post-critique modeling upgrades

Documentation catch-up entry: records pooled-model features already present in the repository but not previously captured in a change log. These were developed in response to the 2026-03-31 model critique:

-   Adaptive AR(1) temporal resolution (daily / weekly / monthly) selected per fit from effort-data density.
-   `L_effective` estimated as a parameter with a lognormal prior from the I/E regression, propagating effective-day-length uncertainty into catch (addresses the critique's primary concern: the snapshot-times-day-length effort expansion).
-   `B1_C` day-type CPUE effect (weekend vs weekday catch rate).
-   Direct I/E crabber-hour integration as lognormal anchor points.
-   Data-driven `R_G` prior centered on the empirical gear-per-crabber ratio.
-   Sparse per-observation overdispersion (`eps_E_H` allocated per observation).
-   Expected (`C_expected`) and predictive (`C`) catch both reported in generated quantities.

### v3-v4 (pooled), v5 (gear-resolved branched)

Earlier shared-sequence milestones are documented in the `README.md` development-history table.
