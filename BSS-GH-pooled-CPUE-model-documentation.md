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

The framework checks convergence using R-hat (\< 1.05) and effective sample size (\> 400) for both C_expected_sum and E_sum (Vehtari et al. 2021). If both criteria are met, the BSS expected catch estimate is preferred. Otherwise, the PE is used.

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

R-hat \< 1.05 AND n_eff \> 400 for C_expected_sum and E_sum (Vehtari et al. 2021). The convergence report records the AR resolution used for each fit, enabling reviewers to assess whether the selected resolution was appropriate.

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
| `daily_combined_estimate.csv`   | Daily PE + BSS estimates with method flag              |

------------------------------------------------------------------------

## 10. Limitations and Future Directions

-   Limited I/E coverage; expand to \~40 days per season for better L_effective regression (Pollock et al. 1997 recommend ≥3 I/E days per month × day-type stratum).
-   L_effective regression uses a quadratic in day-of-year; with more data, a GAM could capture non-monotonic patterns.
-   No weather covariates in effort or CPUE processes; NOAA buoy data could improve prediction on unsampled days (Conn 2002 included temperature in the Kenai River creel model).
-   Gear-type breakdowns are approximate (proportional allocation, not modeled separately).
-   No jetty effort counts. Beach crabbing unmeasured.
-   The B1_C effect is constant across the season; a time-varying weekend CPUE effect may be warranted if tourist composition shifts seasonally.
-   The adaptive AR selection is rule-based; a formal model comparison (LOO-CV or WAIC; Vehtari et al. 2017) could provide principled resolution selection.

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
