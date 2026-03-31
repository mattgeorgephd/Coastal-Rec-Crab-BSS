# Recreational Crab Harvest Estimation — Grays Harbor
## Pooled CPUE Model Documentation

**Authors:** Matt George, with analytical development support  
**Agency:** Washington Department of Fish and Wildlife (WDFW)  
**Date:** March 2026  
**Status:** Proof of concept — 2024-25 season  

---

## 1. Summary for Decision-Makers

This framework estimates the total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area. It combines field observations (gear counts, trailer counts, and crabber interviews) with a statistical model that fills in the days when no sampling occurred.

**Key results (2024-25, through Jan 30):**
- Estimated total harvest: ~32,000–34,000 Dungeness crab
- Shore crabbers (docks, jetty, beach) contribute ~51% of the catch
- Commercial/charter vessels contribute ~30% despite low vessel-hours
- Private boats contribute ~18%
- Pots account for ~57% of catch, traps ~24%, ring nets ~15%, snares ~4%

**What this model does well:** Produces defensible total harvest estimates with uncertainty bounds. When the model says "32,000 ± 3,000 crab," that means we're 95% confident the true harvest falls in that range.

**What this model does not do:** It does not estimate catch separately for each gear type within the statistical model. Gear-type breakdowns are approximate, based on interview proportions applied after estimation. For gear-type-specific estimates with their own uncertainty bounds, see the Gear-Resolved CPUE Model.

---

## 2. Understanding the Two Estimation Methods

Before diving into technical details, it helps to understand *why* this framework produces two sets of numbers (PE and BSS) and what each one means.

### Point Estimator (PE) — The Simple Average

The PE method works like this: look at sampled days, compute the average daily harvest for each week × day-type combination, then multiply by the total number of days in that combination to estimate the unsampled days.

**Strength:** Easy to understand, no modeling assumptions.  
**Weakness:** If the sampled days happen to be unusually busy or quiet, the estimate is biased. Can't estimate harvest for weeks with zero samples.

### Bayesian State-Space Model (BSS) — The Time-Series Model

The BSS method fits a smooth curve through the daily effort and catch rate data, then uses that curve to estimate every day — including days with no sampling. It produces not just a single number but a range of plausible values (a "credible interval").

**Strength:** Fills temporal gaps, accounts for temporal autocorrelation (yesterday's effort predicts today's), produces uncertainty bounds.  
**Weakness:** More complex, slower to run, requires convergence checking.

### Combined Best Estimate

The framework checks whether the BSS model converged properly (R-hat < 1.05 and effective sample size > 400). If it did, the BSS estimate is preferred. If not, the PE is used as a robust fallback.

---

## 3. Glossary

| Term | Meaning |
|---|---|
| **BSS** | Bayesian State-Space model — a time-series approach that estimates hidden daily quantities |
| **PE** | Point Estimator — a simpler stratified-expansion approach |
| **CPUE** | Catch Per Unit Effort — the average number of crab caught per crabber-hour |
| **Crabber-hour** | One person crabbing for one hour (the basic unit of effort) |
| **Effort** | The total amount of crabbing activity, measured in crabber-hours |
| **Gear count** | A snapshot count of crab gear deployed in the water at a point in time |
| **R_G** | Gear-per-crabber ratio — how many pieces of gear each crabber deploys (~1.27) |
| **R_T** | Trailer-per-boat-group ratio — fraction of boat groups that have a trailer at the ramp |
| **AR(1)** | First-order autoregressive process — a model where today's value depends on yesterday's |
| **Credible interval** | The Bayesian equivalent of a confidence interval — a range containing the true value with 95% probability |
| **R-hat** | A convergence diagnostic — should be below 1.05 for reliable results |
| **n_eff** | Effective sample size — how many independent samples the MCMC chain produced |
| **Posterior** | The probability distribution of a parameter after seeing the data |
| **Latent** | A hidden quantity that is not directly observed but inferred from data |
| **Sub-season** | One of two periods: ring-net only (Sep–Nov) or all-gear (Dec–Sep) |
| **Day type** | Weekday, weekend, or holiday — effort differs systematically across these |

---

## 4. Study Area and Data Sources

### 4.1 Westport / Grays Harbor

Westport is a small coastal town on the south side of the Grays Harbor estuary in Washington State. It is one of the highest-volume recreational crabbing ports on the Pacific coast. Crabbing occurs year-round from multiple access points within a few square miles.

### 4.2 Effort Data

Effort is measured through **instantaneous counts** — a field surveyor visits a site and counts the gear or trailers visible at that moment.

| Site | What Is Counted | 2024-25 Records |
|---|---|---|
| Westport Docks Float 20 | Crab gear in the water | 348 |
| Westport Docks Float 17-21 | Crab gear in the water | 169 |
| Westport Boat Launch | Boat trailers at ramp | 226 |
| Ocean Shores Boat Launch | Boat trailers at ramp | 4 |
| Westport Marina | Vessels (limited) | 5 |

Float 20 and Float 17-21 counts are paired by time and summed to produce a section-level gear total.

### 4.3 Interview Data

Creel interviews capture trip-level information from crabbers: group size, gear deployed, hours fished, crab kept by species, and crabbing mode. In 2024-25, 4,359 interviews were conducted. Each is classified into a population:

| crabbing_mode | boat_type | Population |
|---|---|---|
| Dock, Jetty, Beach | any | Shore |
| Boat | Commercial/Charter/Guide | Commercial/Charter |
| Boat | Private or blank | Private Boat |

### 4.4 Commercial/Charter Vessel Tally

A daily tally of commercial and charter vessels was maintained at Westport Marina from December 1, 2024 through February 8, 2025 (47 tally days within a 70-day period).

---

## 5. Season Structure

The 2024-25 season (Sep 16, 2024 – Sep 15, 2025) is split into two sub-seasons at the pot-open date (Dec 1):

| Sub-season | Dates | Duration | Gear Allowed |
|---|---|---|---|
| Ring-net only | Sep 16 – Nov 30 | 76 days | Ring nets, snares, foldable traps |
| All-gear | Dec 1 – Sep 15 | 289 days | All gear including pots |

Each sub-season is estimated independently. This prevents the model from trying to bridge the large discontinuity in effort and CPUE that occurs when pots become legal.

### Day Length

Day length is computed daily from civil twilight (dawn to dusk) using the `suncalc` R package, capped between 9 and 16 hours. This is more accurate than fixed seasonal windows and accounts for crabbers who set gear before dawn and retrieve after dusk.

### Day Type

Each day is classified as weekday (Mon–Thu), weekend (Fri–Sun), or holiday. In the BSS model, effort on these day types is modeled with two separate effects: B1 (weekend boost) and B2 (additional holiday boost beyond B1). This captures the pattern that holidays generate substantially higher effort than regular weekends.

---

## 6. Population Components

The framework estimates three populations independently:

### 6.1 Shore Crabbers

People crabbing from the Westport docks, jetty, or beaches. Effort is measured by gear counts at the docks.

**How gear counts become crabber-hours:** The number of gear in the water is divided by the gear-per-crabber ratio (R_G ≈ 1.27, meaning each crabber deploys about 1.27 pieces of gear), then multiplied by day length.

In equation form: `Gear count ÷ R_G × day length = crabber-hours`

### 6.2 Private Boat Crabbers

People who trailer their boats to the boat launch. Effort is measured by trailer counts.

**How trailer counts become crabber-hours:** Each trailer represents one boat group out crabbing. The trailer-per-group ratio R_T converts trailers to boat groups.

Boat crabbers have substantially higher CPUE (~4.8 crab/trip vs ~1.1 for shore), justifying separate estimation.

### 6.3 Commercial/Charter Vessels

Large vessels moored at Westport Marina. Effort is measured by the daily vessel tally. Estimation uses a day-type stratified census expansion: weekday and weekend/holiday harvest rates are computed separately and expanded to the full census period (Dec 1 – Feb 8), preventing bias from uneven sampling across day types.

---

## 7. The Pooled CPUE Model

### 7.1 How It Works (Non-Technical)

The model assumes that each day has a "true" effort level and a "true" catch rate that we can't directly observe. On days when we have data, the observations pin down these hidden quantities. On days without data, the model fills in the gaps by assuming effort and catch rates change smoothly from day to day. The result is a complete daily time series with uncertainty bounds.

All gear types share a single catch rate process — the model doesn't distinguish whether crab were caught by pots, ring nets, or traps.

### 7.2 How It Works (Technical)

**Effort process:**
```
log(lambda_E[d]) = mu_E + omega_E[period[d]] + B1 × w[d] + B2 × holiday[d]
```
In words: daily effort on the log scale equals a baseline level (`mu_E`), plus a smooth deviation that changes each period (`omega_E`, following an AR(1) process), plus a weekend boost (`B1`), plus an additional holiday boost (`B2`).

**CPUE process:**
```
log(lambda_C[d]) = mu_C + omega_C[period[d]]
```
In words: daily catch rate follows a similar smooth process but without a weekend effect (CPUE is assumed constant across day types within a period).

**Observation models:**
- Gear counts: `count ~ Poisson(crabbers × gear_per_crabber × time_of_day_noise)`
- Trailer counts: `count ~ Poisson(boat_groups × trailers_per_group × time_of_day_noise)`
- Interview catch: `crab_caught ~ NegBin(catch_rate × hours_fished, overdispersion)`

**Daily catch:** `C[d] = Poisson_draw(crabbers[d] × day_length[d] × catch_rate[d])`

### 7.3 Stan Model File

The Stan model is `crab_bss_pooled.stan`. Key parameters:

| Parameter | Description | Prior |
|---|---|---|
| B1 | Weekend effort multiplier (log scale) | Normal(0, 1) |
| B2 | Additional holiday effort multiplier | Normal(0, 1) |
| R_G | Gear per crabber | Lognormal(log(1.3), 0.3) |
| R_T | Trailers per boat group | Beta(0.5, 0.5), guarded |
| phi_E, phi_C | AR(1) autocorrelation | Beta(2, 2), rescaled to [-1, 1] |
| r_E, r_C | Overdispersion (effort, CPUE) | Half-Cauchy |

---

## 8. Convergence and Model Selection

For each BSS fit, convergence is assessed using R-hat (< 1.05) and effective sample size (> 400). Fits passing both thresholds use the BSS estimate; otherwise PE is used as a fallback.

---

## 9. Output Files

Each run produces output in `output/YYYYMMDD/`:

### Tables
| File | Contents |
|---|---|
| `pe_port_summary.csv` | PE estimates by component and port total |
| `port_total_Dungeness_Kept.csv` | Combined PE + BSS port total |
| `monthly_by_population.csv` | Month × Population with PE, BSS, and Combined estimates |
| `monthly_port_totals.csv` | Port-level monthly totals |
| `monthly_by_mode.csv` | Monthly catch by crabbing mode |
| `catch_by_gear_type.csv` | Catch by gear type (proportional allocation from interviews) |
| `monthly_by_area.csv` | Monthly catch by creel area |
| `wide_summary.csv` | All dimensions in one table |
| `bss_summary_{label}.csv` | Stan convergence diagnostics per fit |
| `run_parameters.txt` | Full parameter list for reproducibility |

---

## 10. Key Assumptions

1. **Gear counts represent relative effort** across days.
2. **Gear-per-crabber ratio is constant** within a sub-season.
3. **Trailer counts represent boat crabber effort.**
4. **Day length from civil twilight is a reasonable proxy** for available fishing hours.
5. **Interviews are representative** of all crabbers in the population.
6. **Released crab mortality is not estimated.**
7. **The AR(1) process captures temporal dynamics** — effort and CPUE change smoothly.
8. **Within-day overdispersion is gamma-distributed.**
9. **Holidays and weekends are pooled as "high-effort days"** in the effort process, with B2 capturing the additional holiday boost.
10. **Mean catch per commercial vessel is constant** within each day-type stratum across sampled and unsampled days.

---

## 11. Limitations and Planned Improvements

### Current Limitations
- No jetty effort counts (28 interviews contribute to CPUE only)
- Single-count days in ring-net sub-season limit within-day overdispersion estimation
- Beach crabbing unmeasured (1 interview in 2024-25)
- **Gear-type breakdowns are approximate** — derived from interview proportions, not modeled

### Planned Improvements
1. Add Westport Jetty effort counts
2. Restructure `eps_E_H` allocation for data-informed days only
3. Expand to other ports (Willapa Bay, Columbia River)

---

## 12. Reproducibility

1. Clone the `FWC-estimation-method` repository
2. Place input files in `input_files/`
3. Place `crab_bss_pooled.stan` in `stan_models/`
4. Open `BSS-GH-pooled-CPUE-model.Rmd`
5. Set parameters and run

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here  
**Expected runtime:** ~3-4 hours on a 4-core machine
