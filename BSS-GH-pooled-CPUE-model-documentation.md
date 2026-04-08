# Recreational Dungeness Crab Harvest Estimation — Grays Harbor
## Pooled CPUE Model: Technical Documentation

**Author:** Matthew George, Ph.D.  
**Contact:** matthew.george@dfw.wa.gov  
**Agency:** Washington Department of Fish and Wildlife (WDFW)  
**Status:** Operational — Annual estimation framework  

---

## 1. Summary for Decision-Makers

This framework estimates the total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area. It combines four types of field observations — gear counts at the docks, trailer counts at the boat launch, dockside crabber interviews, and ingress/egress surveys — with a statistical model that fills in the days when no sampling occurred. The result is a season-long harvest estimate for the port, broken down by month, crabbing mode, and gear type, with a measure of how uncertain the estimate is.

**What this model produces:**

- A total Dungeness crab harvest estimate for the port with a 95% credible interval (e.g., "an estimated 70,000 crab were harvested, with 95% probability that the true number falls between 62,000 and 79,000").
- Monthly harvest trends showing when crabbing pressure peaks and how it changes through the season.
- Breakdowns by crabbing mode (shore, private boat, commercial/charter) showing which user group contributes what share of the total.
- Gear-type harvest proportions (pots, ring nets, traps, snares) based on what crabbers report in interviews.
- A temporal correction factor (f_temporal) quantifying how well gear-count-based effort estimates match directly-observed crabber-hours from I/E surveys.

**What this model does not do:** It does not model each gear type's catch rate independently. Gear-type breakdowns are approximate, derived from proportions observed in interviews and applied to the total catch estimate after estimation. For gear-type-specific estimates with their own uncertainty bounds, see the companion Gear-Resolved CPUE Model.

**How confident are we?** The framework runs two independent estimation methods — a simple average-based approach (PE) and a Bayesian time-series model (BSS) — then compares them. When the two methods agree closely and the BSS model converges properly, confidence is high. On days with ingress/egress surveys, the effort estimate is anchored by direct observation rather than relying on the gear-count conversion, providing an independent calibration check. The output includes formal convergence diagnostics and a comparison table so reviewers can assess reliability.

---

## 2. Understanding the Two Estimation Methods

The framework produces two sets of numbers because no single method is ideal in all situations.

### 2.1 Point Estimator (PE) — The Simple Average

The PE method works like this: look at the days when surveyors were in the field, compute the average daily harvest for each week × day-type (weekday vs weekend) combination, then multiply by the total number of days in that combination to estimate harvest on unsampled days.

**Strengths:** Easy to understand. Transparent calculations. No modeling assumptions beyond "sampled days represent unsampled days within a stratum."

**Weaknesses:** If the sampled days happen to be unusually busy or quiet, the estimate is biased. Cannot produce estimates for time periods with zero samples. Does not produce uncertainty bounds. Treats each stratum independently — a busy Saturday provides no information about the following Tuesday.

### 2.2 Bayesian State-Space Model (BSS) — The Time-Series Model

The BSS method fits a smooth curve through the daily effort and catch rate data, then uses that curve to estimate every day in the season — including days with no field sampling. It works by assuming that daily effort and catch rates evolve smoothly over time: yesterday's effort level is informative about today's. The model produces a full probability distribution of plausible values (a "posterior distribution"), from which we extract a median estimate and a 95% credible interval.

On days with ingress/egress surveys, the BSS receives a direct crabber-hours observation in addition to the gear count. This dual-observation design lets the model calibrate the gear-count-to-effort mapping within the estimation, rather than relying on external assumptions.

**Strengths:** Fills temporal gaps by borrowing information from neighboring days. Accounts for temporal autocorrelation. Produces rigorous uncertainty bounds. I/E anchor points constrain the effort trajectory.

**Weaknesses:** More complex. Requires 30–60 minutes per model fit. Must be checked for convergence — if the model doesn't converge, the estimates are unreliable.

### 2.3 Combined Best Estimate

The framework automatically checks whether each BSS fit converged properly using two standard diagnostics: R-hat (should be below 1.05) and effective sample size (should exceed 400). If both criteria are met, the BSS estimate is preferred. If convergence fails, the PE is used as a robust fallback. The output files include a convergence report showing which method was selected for each component and why.

---

## 3. The Recreational Crab Fishery at Grays Harbor

### 3.1 Fishery Overview

The recreational Dungeness crab (*Metacarcinus magister*) fishery at Westport and the greater Grays Harbor area is one of the highest-volume recreational crabbing operations on the Washington coast. Thousands of crabbers participate annually, using a variety of gear types from shore-based locations and private boats. The fishery operates year-round, though effort peaks during summer months and around major holidays.

Recreational crabbers at Grays Harbor use four primary gear types: crab pots (the highest-catch gear), ring nets, foldable/star traps, and snares. WDFW regulations restrict pot use to a specific portion of the season (typically December through September), creating a structural break in both effort and catch rates when pots become available.

Commercial Dungeness crab vessels also participate in the recreational fishery before the commercial season opens — these large vessels crab recreationally under the same daily limits as private boats but tend to have much higher catch rates per vessel. Their harvest is tracked separately through a vessel tally system at the marina.

### 3.2 Why Estimate Harvest?

WDFW needs accurate recreational harvest estimates to monitor total removals against sustainable yield targets, evaluate the effectiveness of gear restrictions and season structures, inform allocation decisions between recreational and commercial sectors, and track long-term trends in recreational participation and catch rates.

### 3.3 Study Area

Westport is a small coastal town on the south side of the Grays Harbor estuary. Recreational crabbing occurs from multiple access points within a few square miles: a system of public docks (Floats 17–21), a jetty, beaches, a public boat launch, and a commercial marina. Each access point serves a different crabbing mode and requires a different type of effort measurement.

---

## 4. Data Sources and Field Collection

The estimation framework uses four types of field data.

### 4.1 Effort Counts

**What they are:** Instantaneous point-in-time counts of gear or trailers, conducted by a field surveyor who visits a site and records the number of crabbing indicators visible at that moment.

**Why they matter:** Effort counts are the primary input for estimating total crabbing activity on any given day. Since not every crabber is interviewed, effort counts provide an independent measure of how many people are actively crabbing.

**How they are collected:** Field staff conduct effort counts during scheduled survey events. Protocol calls for multiple counts per day at standardized times. The number of within-day counts matters because the BSS model uses replication to estimate within-day variability.

**Sites and count types:**

| Site | What Is Counted | Why This Indicator | Role in Model |
|---|---|---|---|
| Westport Docks Float 20 | Crab gear in the water | Each piece of gear = an active crabber | Primary shore effort indicator |
| Westport Docks Float 17-21 | Crab gear in the water | Paired with Float 20 for full dock coverage | Summed with Float 20 for section total |
| Westport Boat Launch | Boat trailers at ramp | Each trailer = a boat group out crabbing | Private boat effort indicator |
| Ocean Shores Boat Launch | Boat trailers at ramp | Secondary launch on north shore | Supplementary boat effort |
| Westport Jetty | Crabbers (future) | Direct person count | Reserved for future use |

**Assumptions:**
- A point-in-time count reflects relative effort on that day.
- The time of day introduces noise but not systematic bias.
- The surveyor's count is accurate — all deployed gear is visible and countable.

### 4.2 Crabber Interviews

**What they are:** Dockside interviews conducted by WDFW field staff who approach crabbers and record trip-level information.

**Why they matter:** Interviews provide the catch rate (crab per unit effort) and the gear-per-crabber ratio (R_G) — quantities that effort counts alone cannot measure.

**What each interview records:** Number of crabbers in the group, number of gear units deployed, gear type(s) used, hours fished, crabber-hours, gear-hours, Dungeness crab kept, Red Rock crab kept, trip completion status, crabbing mode, and boat type.

**CPUE denominator by population:** Shore interviews use crabber-hours (crabbers × hours fished). Boat interviews use gear-hours (gear units × hours fished), because boat gear soaks continuously and the relevant effort unit is gear deployment time, not crabber presence time.

**Population classification:**

| Crabbing Mode | Boat Type | Population |
|---|---|---|
| Dock, Jetty, or Beach | Any | **Shore** |
| Boat | Private or blank | **Private Boat** |
| Boat | Commercial, Charter, or Guide | **Commercial/Charter** |

**Assumptions:**
- Interviewed crabbers are representative of all crabbers in their population.
- Crabbers accurately report catch, hours fished, and group size.

### 4.3 Commercial/Charter Vessel Tally

A daily count of commercial crab vessels and charter boats at Westport Marina during the period when these vessels participate in the recreational fishery. Combined with mean catch per vessel from interviews for stratified expansion.

### 4.4 Ingress/Egress (I/E) Surveys

**What they are:** All-day surveys where a field staff member records the number of crabbers arriving at and departing from a site every 15 minutes. At the boat launch, boat arrivals and departures are tracked instead.

**Why they matter:** I/E surveys provide a direct measurement of daily crabber-hours that bypasses the conversion chain used by gear counts (gear count → R_G division → day-length multiplication). A gear count is a snapshot at one moment; an I/E survey tracks the full activity curve across the day.

**How crabber-hours are computed:** The running total (`crabber_flow`) tracks crabbers currently present at each 15-minute interval. Daily crabber-hours equals the sum of `crabber_flow × 0.25` across all intervals — the area under the crabber-presence curve.

**How I/E data enters the model:** On I/E survey days, the BSS receives the observed crabber-hours as a direct lognormal observation of lambda_E × L, with measurement error sigma_IE. This dual-observation design (gear count + I/E on the same day) lets the model calibrate the gear-count pathway against the I/E ground truth.

**What the I/E data reveals:** The effective day length (crabber-hours ÷ peak crabbers present) averages approximately 3.5–5.0 hours at Float 20, substantially shorter than the civil-twilight-based day length (9–16 hours). Crabbers rotate through the dock in a peaked activity curve, not a flat all-day presence. The f_temporal correction factor varies with the time of day the gear count is taken.

**Sites:** WDF20 (Westport Docks Float 20) and WBL (Westport Boat Launch).

---

## 5. Key Design Decisions

### 5.1 Sub-Seasons

**Ring-net only** (Sep 16 – Nov 30): Crab pots prohibited. **All-gear** (Dec 1 – Sep 15): Pots allowed. Each sub-season is estimated independently.

### 5.2 Day Length and Effort Units

**Shore crabbers — I/E-derived L_effective:** Rather than using civil twilight (9–16 hours), the framework estimates an empirical "effective day length" from all available historical I/E data. L_effective is defined as total crabber-hours ÷ peak crabbers present on each I/E survey day. It captures the fact that crabbers rotate through the dock in a peaked activity curve — the peak count at any moment substantially exceeds the average count across the day.

The L_effective model groups I/E days into three sub-season periods aligned with the pot closure:

| Sub-season period | Months | Rationale |
|---|---|---|
| **ring_net** | Sep 16 – Nov 30 | Pots prohibited. Shorter, more active ring-net trips with frequent gear checking. |
| **allgear_winter** | Dec 1 – Feb 28 | Pots open. Higher CPUE. Pot crabbers may "set and leave" for hours, creating the largest divergence between gear counts and crabber presence. |
| **allgear_spring_summer** | Mar 1 – Sep 15 | Longer days, more tourist crabbers, potentially different turnover patterns. |

Within each sub-season period, L_effective is estimated separately for weekdays and weekends (with holidays pooled as weekends). For cells with fewer than 2 I/E observations, the grand mean across all I/E days is used as a conservative fallback.

The original civil twilight day length is retained as `day_length_civil_twilight` for reference and for any downstream analysis that needs the astronomical value. The outputs include `L_effective_lookup.csv` (the 6-cell lookup table) and `L_effective_ie_detail.csv` (per-day I/E measurements) for transparency.

**Why this matters:** Analysis of I/E data at Float 20 shows L_effective averaging 3.5–5.0 hours — roughly half the civil twilight value. Using civil twilight systematically overestimates daily shore effort (and therefore harvest). The L_effective model corrects this using direct field observations. As the I/E dataset grows over future seasons, the cell estimates will become more precise and the model can be refined (e.g., monthly rather than seasonal grouping).

**When no I/E data is available:** If the I/E file is missing or `use_ie_day_length = FALSE`, the framework reverts to civil twilight day length with no change in behavior.

**Private boats:** Day length is fixed at **24 hours** because boat crab gear (primarily pots) soaks continuously. The boat effort unit is **gear-hours** (gear units × hours deployed; CPUE is crab per gear-hour).

### 5.3 Day Type

Weekday (Mon–Thu), weekend (Fri–Sun), or holiday. Separate B1 (weekend) and B2 (holiday) effort effects in the BSS.

---

## 6. Population Components

### 6.1 Shore Crabbers (Dock + Jetty + Beach)

Effort indicator: gear counts at the docks. Conversion: `(Gear counted ÷ R_G) × day_length = crabber-hours`. On I/E days, crabber-hours are observed directly.

### 6.2 Private Boat Crabbers

Effort indicator: trailer counts at boat launches. Conversion: `trailer_count × gear_per_group × 24 = gear-hours`. CPUE uses gear-hours as the denominator.

### 6.3 Commercial/Charter Vessels

Estimated via day-type stratified census expansion from the vessel tally.

---

## 7. The Pooled CPUE Model (Technical)

### 7.1 Effort Process

```
log(lambda_E[d]) = mu_E + omega_E[period[d]] + B1 × w[d] + B2 × holiday[d]
```

AR(1) structure: `omega_E[p] = phi_E × omega_E[p-1] + innovation`.

### 7.2 CPUE Process

```
log(lambda_C[d]) = mu_C + omega_C[period[d]]
```

### 7.3 Observation Models

- **Gear counts (shore):** `Gear_I ~ Poisson(lambda_E × eps_E_H_obs × R_G)`
- **Trailer counts (boats):** `T_I ~ Poisson(lambda_E × eps_E_H_obs × R_T)`
- **Interview catch:** `c ~ NegBin(lambda_C × h, r_C)` where h = crabber-hours (shore) or gear-hours (boats)
- **I/E crabber-hours:** `IE_crabber_hours ~ Lognormal(log(lambda_E × L), sigma_IE)`

### 7.4 Sparse Overdispersion

`eps_E_H_obs` allocated only for actual effort observations. Eliminates 64–77% of effort parameters.

### 7.5 I/E Integration (Option 2) — How Ingress/Egress Data Enters the Model

#### The problem the I/E data solves

On a typical survey day, the model sees a gear count — say 50 pieces of gear at the dock at 10:30am. To turn this into daily crabber-hours (which is what gets multiplied by CPUE to produce catch), the model has to do three things: divide by R_G (~1.27 gear per crabber) to get ~39 crabbers present at that moment, then multiply by day length (~10 hours) to get ~390 crabber-hours. That day-length multiplication assumes those 39 crabbers were there all day. They weren't — crabbers rotate through. I/E data from February 10, 2024 showed 117 crabbers at peak but only 581 true crabber-hours, meaning an effective day length of about 5 hours, not 10.

#### What the model "sees" without I/E

The BSS has a latent effort variable `lambda_E[d]` for each day — the model's belief about how much crabbing happened. On days with gear counts, the model observes:

```
Gear_count ~ Poisson(lambda_E[d] × eps_overdispersion × R_G)
```

This tells the model about lambda_E, but only through the lens of R_G and overdispersion. The model then computes daily effort as `E[d] = lambda_E[d] × L[d]` (day length), and daily catch as `C[d] = lambda_E[d] × L[d] × lambda_C[d]`. If L[d] is wrong, both the effort and catch estimates inherit that bias.

#### What the model "sees" with I/E

On I/E days, the model gets a second observation of the same latent effort state:

```
IE_crabber_hours ~ Lognormal(log(lambda_E[d] × L[d]), sigma_IE)
```

This is a direct observation of `E[d]` — the quantity the model is trying to estimate — with some measurement noise (sigma_IE). It doesn't go through R_G, doesn't go through overdispersion, and doesn't depend on when during the day the gear count was taken. The I/E surveyor tracked every arrival and departure across the full day, so the resulting crabber-hours *is* the daily effort.

#### How the two observations interact on paired days

On a day with both a gear count and an I/E survey, the model has two independent signals constraining lambda_E[d]:

1. The gear count says: "lambda_E is approximately Gear_count / (R_G × eps)" — but this is noisy because R_G is uncertain and eps adds overdispersion.

2. The I/E says: "lambda_E × L is approximately IE_crabber_hours" — this is much more precise (sigma_IE is small, ~0.2 on log scale, meaning ~±20% measurement error).

The model resolves these two signals through the posterior. If the gear-count pathway consistently implies higher effort than the I/E observations, the model can adjust by: (a) shifting R_G upward (more gear per crabber means the same gear count implies fewer crabbers), (b) adjusting the effort process level downward, or (c) absorbing it through the overdispersion terms. The key point is that the model figures out the calibration internally rather than requiring an externally imposed correction factor.

#### What happens on non-I/E days

On the ~95% of days without I/E surveys, the model only has gear counts. But the I/E anchor points have already constrained the effort trajectory — the AR(1) process smoothly interpolates between them. And the posterior for R_G has been informed by the dual-observation days, so the gear-count conversion is now calibrated. The I/E data effectively "teaches" the model how to interpret gear counts more accurately, even on days when no I/E survey was conducted.

#### Why lognormal for the I/E likelihood

The I/E measurement has multiplicative error — if the surveyor misses 10% of arrivals, that's a proportional undercount regardless of whether 20 or 200 crabbers came through. Lognormal naturally handles this: `sigma_IE = 0.2` means the true value is within about ±20% of the observation, with the uncertainty being proportional to the magnitude.

#### What sigma_IE tells you

The posterior for sigma_IE is itself informative. If it comes back very small (~0.05–0.10), the I/E measurements are highly precise and the model trusts them almost exactly. If it's larger (~0.3+), there's substantial discrepancy between what the I/E predicts and what the rest of the model expects — which could indicate that the I/E survey window didn't cover the full crabbing day, or that the gear-count pathway has a systematic bias the model is struggling to reconcile.

#### Analogy

Think of it like having two thermometers. The gear count is a cheap thermometer that reads through a window — indirect, noisy, and potentially biased by when you look. The I/E survey is a precision thermometer placed directly in the room. On days when you have both, the model can figure out how much the cheap thermometer is off. Then on days when you only have the cheap thermometer, it applies that learned correction.

#### When no I/E data is available

When `IE_n = 0`, the I/E likelihood contributes nothing and the model reverts to the gear-count-only pathway with no change in behavior. This makes the I/E integration fully backward-compatible — it improves estimates when I/E data is available but doesn't alter the model structure when it isn't.

### 7.6 Key Parameters

| Parameter | Description | Prior |
|---|---|---|
| B1 | Weekend effort multiplier (log scale) | Normal(0, 1) |
| B2 | Holiday effort multiplier | Normal(0, 1) |
| R_G | Gear per crabber (shore) | Lognormal(log(1.3), 0.3) |
| R_T | Trailers per boat group | Beta(0.5, 0.5) |
| phi_E, phi_C | AR(1) autocorrelation | Beta(2,2) rescaled [-1,1] |
| r_E, r_C | Overdispersion | Half-Cauchy(5) |
| sigma_IE | I/E measurement error (log scale) | Exponential(5) |

**Stan model file:** `crab_bss_pooled.stan`

---

## 8. Convergence and Model Selection

R-hat < 1.05 AND n_eff > 400 for C_sum and E_sum. The `convergence_report.csv` output records diagnostics for every fit.

---

## 9. Output Files

16+ CSV files, 10+ plots, and run parameters. Key additions: `ie_analysis.csv` (I/E validation with f_temporal and L_effective), `L_effective_lookup.csv` (sub-season × day-type predicted day lengths), `L_effective_ie_detail.csv` (per-day I/E measurements), `plot_L_effective_ie.png` (historical I/E L_effective with model predictions), and `plot_day_length_comparison.png` (civil twilight vs I/E day length across the season).

---

## 10. Limitations and Future Directions

- Limited I/E coverage; expand to ~40 days per season with month × day-type stratification.
- L_effective model has thin cells (some sub-season × day-type combinations have < 2 I/E days); as the I/E dataset grows, finer temporal resolution (monthly grouping) will become feasible.
- No jetty effort counts. Beach crabbing unmeasured.
- Gear-type breakdowns are approximate (proportional allocation, not modeled).
- The L_effective model assumes the structural relationship between I/E crabber-hours and peak count is stable across years within a sub-season × day-type cell.
- Add per-gear catch recording in interviews for future gear-resolved analysis.

---

## 11. Glossary

| Term | Meaning |
|---|---|
| **BSS** | Bayesian State-Space model |
| **PE** | Point Estimator |
| **CPUE** | Catch Per Unit Effort |
| **Crabber-hour** | One person crabbing for one hour |
| **Gear-hour** | One piece of crab gear deployed for one hour |
| **R_G** | Gear-per-crabber ratio |
| **I/E** | Ingress/Egress survey |
| **f_temporal** | Temporal correction factor: I/E crabber-hours ÷ gear-count-derived estimate |
| **L_effective** | Effective day length: I/E crabber-hours ÷ peak crabbers present |
| **L_effective model** | Sub-season × day-type regression predicting L_effective from historical I/E data |
| **sigma_IE** | I/E measurement error on log scale |
| **AR(1)** | First-order autoregressive process |
| **Credible interval** | Bayesian range containing the true value with stated probability |

---

## 12. Reproducibility

1. Clone the repository. Place input files in `input_files/`: `effort_combined.csv`, `interview_combined.csv`, `wes_commercial_tally.csv`, `Matt_Ingress-Egress_Compilation.xlsx`.
2. Place `crab_bss_pooled.stan` in `stan_models/`.
3. Open `BSS-GH-pooled-CPUE-model.Rmd`, update parameters, run.

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here, readxl.  
**Expected runtime:** 3–4 hours on a 4-core machine.
