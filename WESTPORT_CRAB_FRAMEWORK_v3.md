# Westport Crab Creel Estimation Framework
## Technical Documentation

**Authors:** Matt George, with analytical development support  
**Agency:** Washington Department of Fish and Wildlife (WDFW)  
**Date:** March 2026  
**Status:** Proof of concept — 2024-25 season  

---

## 1. Purpose

This document describes the statistical framework used to estimate total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area. The framework uses a Bayesian state-space (BSS) model alongside a traditional point estimator (PE), adapted from WDFW's freshwater creel estimation system for the distinct characteristics of the ocean recreational crab fishery.

### Why a purpose-built framework?

Recreational crabbing at Westport differs from river-based fisheries in several important ways:

- Effort is measured by counting **gear in the water** (crab pots, ring nets, snares), not people or vehicles.
- Multiple **crabbing modes** coexist at the same port — dock crabbers, jetty crabbers, private boat crabbers, and commercial/charter vessels all contribute to total harvest.
- Each mode has a **different effort indicator** (gear counts for docks, trailer counts for boat launches, vessel tallies for the marina).
- **Gear restrictions** change mid-season: crab pots are prohibited from September 16 through November 30, creating a structural break in both effort and catch rates.
- **Commercial crab vessels** participate in the recreational fishery before the commercial season opens, requiring a separate estimation approach.

The framework addresses each of these differences while preserving the core statistical machinery (AR(1) latent effort and CPUE processes, negative-binomial catch likelihood, gamma within-day overdispersion) that makes the BSS model effective for filling temporal gaps in creel data.

---

## 2. Study Area and Data Sources

### 2.1 Westport / Grays Harbor

Westport is a small coastal town on the south side of the Grays Harbor estuary in Washington State. It is one of the highest-volume recreational crabbing ports on the Pacific coast. Crabbing occurs year-round from multiple access points within a few square miles.

### 2.2 Effort Data

Effort is measured through **instantaneous counts** conducted by field staff during creel surveys. The type of count varies by location:

| Site | What Is Counted | Count Type | 2024-25 Records | Days Sampled |
|---|---|---|---|---|
| Westport Docks Float 20 | Crab gear in the water | Gear count | 348 | 186 |
| Westport Docks Float 17-21 | Crab gear in the water | Gear count | 169 | 162 |
| Westport Boat Launch | Boat trailers at ramp | Trailer count | 226 | 174 |
| Ocean Shores Boat Launch | Boat trailers at ramp | Trailer count | 4 | 4 |
| Westport Marina | Vessels (limited) | Gear count | 5 | 3 |
| Westport Jetty | None in 2024-25 | — | 0 | 0 |

**Float 20 and Float 17-21** are paired: Float 20 is the primary count location; Float 17-21 counts are matched to the nearest-in-time Float 20 count and summed to produce a section-level gear total. This pairing accounts for the spatial distribution of dock crabbers across two adjacent float areas.

**Effort count protocol:** Prior to approximately March 2025, surveyors conducted a single effort count per day. Starting in March 2025, the protocol shifted to three counts per day at standardized times. This transition matters because the BSS model uses within-day replication (multiple counts on the same day) to estimate effort overdispersion — a key component of the state-space structure. Sub-seasons falling entirely before this transition (the ring-net sub-season, Sep–Nov) have limited within-day replication.

### 2.3 Interview Data

Creel interviews are conducted by field staff who approach crabbers and record trip-level information. Each interview captures the number of crabbers in the group, gear deployed, hours fished, crab kept (by species), trip completion status, and crabbing mode. In 2024-25, 4,359 interviews were conducted across all Grays Harbor sites. After filtering to valid trips (crabbers > 0, fishing time >= 0.5 hours), approximately 3,925 interviews are available for analysis.

Each interview is classified into a **population** based on two fields:

- `crabbing_mode`: Dock, Jetty, Beach, or Boat
- `boat_type`: Private, Commercial Crab Vessel, Charter, or Guide

The classification rules are:

| crabbing_mode | boat_type | Assigned Population | Rationale |
|---|---|---|---|
| Dock | any | Shore | Crabbing from the docks |
| Jetty | any | Shore | Crabbing from the jetty |
| Beach | any | Shore | Crabbing from beaches |
| Boat | Commercial or Charter or Guide | Commercial/Charter | Vessels moored at marina |
| Boat | Private or blank | Private Boat | Trailered boats from boat launch |

**"Boat" interviews at the docks:** Approximately 104 interviews at Float 20 and Float 17-21 have `crabbing_mode = "Boat"` with `boat_type = "Private"`. These are private boat crabbers who launched from the boat launch but pulled up to the docks for offloading or interviewing. Their **effort** is captured by the trailer counts at the boat launch (not the dock gear counts), so they are correctly classified as private boat population for both effort and CPUE purposes.

**Marina interviews:** The 312 interviews at Westport Marina are split by `boat_type`:
- 141 Commercial Crab Vessel -> Commercial/Charter population
- 26 Charter -> Commercial/Charter population
- 34 Private -> Private Boat population
- 111 Dock mode -> Shore population

The `boat_type` field contains a known typo in the iForm export: "Commerical" (one 'm') rather than "Commercial". The classification code handles this with a case-insensitive regex match.

### 2.4 Commercial/Charter Vessel Tally

A separate daily tally of vessels was maintained at Westport Marina from December 3, 2024 through February 8, 2025 (47 days). This tally recorded the number of private, commercial, and charter vessels observed each day, along with the number actually interviewed. This data source is used for the commercial/charter census estimation (Section 5.3) and is not fed into the BSS model.

### 2.5 Data Files

| File | Location | Contents |
|---|---|---|
| `effort_combined.csv` | `input_files/` | All effort count records, all seasons. Re-exported with `QUOTE_ALL` to handle commas in the notes field that otherwise corrupt CSV parsing. |
| `interview_combined.csv` | `input_files/` | All interview records, all seasons. Dates in M/D/YYYY format (read with `col_date(format="%m/%d/%Y")`). The `number_of_gear` column must be mapped from column N (not column W) in the raw iForm export due to a duplicate field name bug. |
| `wes_commercial_tally.csv` | `input_files/` | Daily vessel tally from the "wes commercial tally" sheet of the 2024-25 seasonal workbook. |

---

## 3. Season Structure

### 3.1 Full Season

The 2024-25 crab season runs from **September 16, 2024 through September 15, 2025** (365 days). This annual cycle corresponds to the state recreational crabbing regulation year.

### 3.2 Gear-Regime Sub-Seasons

Washington State prohibits the use of crab pots in the recreational fishery from September 16 through November 30 each year. During this period, only ring nets, snares, foldable traps, and handlines are permitted. Pots become legal on December 1.

This gear restriction creates a **structural break** in both effort intensity and catch rates. When pots become available, gear counts at the docks increase dramatically (from ~20-50 in Oct-Nov to 100+ in the following months), and CPUE shifts because pots are more effective than ring nets.

The framework splits the season into **two independent sub-seasons**:

| Sub-season | Dates | Duration | Gear Allowed |
|---|---|---|---|
| Ring-net only | Sep 16 – Nov 30 | 76 days | Ring nets, snares, foldable traps, handlines |
| All-gear (pot open) | Dec 1 – Sep 15 | 289 days | All gear including pots |

Each sub-season is estimated independently with its own BSS model fit (or PE estimate). This approach:

- Avoids forcing the AR(1) process to bridge a discontinuity in effort and CPUE
- Allows each gear regime to have its own baseline levels
- Keeps each model fit smaller and faster
- Makes the gear-regime effect explicit in the results

The `pot_open_date` parameter controls the split and can be adjusted if regulations change.

### 3.3 Commercial/Charter Active Period

Commercial crab vessels and charter boats participate in the recreational crab fishery primarily before the **commercial crab season opener**. Once the commercial season opens, these vessels shift to commercial crabbing and their recreational harvest becomes negligible.

The commercial opener date varies by year and is set as a parameter:

```r
commercial_opener = "2025-01-15"   # 2024-25 season
```

The commercial/charter estimation is truncated at this date.

### 3.4 Day Type Classification

Each day in the season is classified into one of three types:

- **Weekday:** Monday through Thursday
- **Weekend:** Friday, Saturday, Sunday
- **Holiday:** Specific high-effort dates designated as crabbing holidays

The 2024-25 holidays are:

| Date | Holiday |
|---|---|
| Sep 2, 2024 | Labor Day |
| Nov 29, 2024 | Native American Heritage Day |
| Dec 31, 2024 | New Year's Eve |
| Jan 1, 2025 | New Year's Day |
| Feb 8, 2025 | Super Bowl Eve |
| May 26, 2025 | Memorial Day |
| Jun 15, 2025 | Father's Day |

In the BSS model, weekends and holidays are combined into a single binary indicator (`w[d] = 1` for weekend/holiday, `0` for weekday) that enters the effort process as a fixed effect (`B1 * w[d]`). In the PE method, holidays are treated as their own stratum, separate from both weekdays and weekends.

### 3.5 Day Length

Day length determines how many hours of crabbing effort are possible each day. The framework uses **fixed seasonal windows**:

| Months | Day Length | Rationale |
|---|---|---|
| April – September | 10.0 hours | Longer summer days |
| October – March | 8.5 hours | Shorter winter days |

This simplification is appropriate because dock crabbing does not strictly follow sunrise/sunset — crabbers often deploy gear before dawn and retrieve after dusk, and the docks are lit. Day length primarily serves as a scaling factor to convert instantaneous effort counts (crabbers present at a point in time) to daily crabber-hours.

---

## 4. Population Components

The framework estimates harvest independently for three crabbing populations and sums the results for the port total:

```
Total Harvest = Shore Crabbers + Private Boat Crabbers + Commercial/Charter Vessels
```

### 4.1 Shore Crabbers

**In plain language:** People who crab from the Westport docks, the jetty, or nearby beaches. This is the largest group by interview count and the primary focus of the creel survey.

**Population includes:** All interviews with `crabbing_mode` of Dock, Jetty, or Beach, plus Dock-mode interviews at the Marina.

**2024-25 data:**
- Ring-net sub-season: 52 effort observations, ~856 interviews
- All-gear sub-season: 313 effort observations, ~2,743 interviews

**Effort observation model:** Gear counts at the docks measure crabbing equipment deployed in the water. The latent quantity of interest is the number of *crabbers* (people), not gear units. The relationship is mediated by the parameter **R_G** (gear per crabber):

```
Observed gear count ~ Poisson(lambda_E x eps_H x R_G)
```

Where:
- `lambda_E` is the latent true crabber count at that point in time
- `eps_H` is within-day overdispersion (gamma-distributed; accounts for effort fluctuating throughout the day)
- `R_G` is the gear-per-crabber ratio (~1.27 from interview data)

R_G is estimated jointly from two sources: the gear count observations themselves and interview expansion data (`Gear_A ~ Poisson(A_A x R_G)`), where both gear units and crabber count are observed per group. With 3,000+ expansion interviews informing R_G, its posterior is tightly constrained. The prior is weakly informative: `R_G ~ Lognormal(log(1.3), 0.3)`.

**CPUE model:**
```
crab_caught ~ NegBin(lambda_C x crabber_hours, r_C)
```

The negative binomial distribution handles the substantial overdispersion in crab catch (many groups catch zero while some catch the daily limit).

### 4.2 Private Boat Crabbers

**In plain language:** People who trailer their own boats to the boat launch, motor into Grays Harbor, and crab from their boats. They are interviewed either at the boat launch upon return or at the docks if they pull up there.

**Population includes:** All interviews with `crabbing_mode = "Boat"` and `boat_type = "Private"` (or blank), regardless of interview location.

**2024-25 data:**
- Ring-net sub-season: 44 effort observations, ~86 interviews
- All-gear sub-season: 196 effort observations, ~336 interviews

**Effort observation model:** Trailer counts at the boat launch measure the number of boat trailers parked at the ramp, each corresponding to one boat group currently out crabbing:

```
Observed trailer count ~ Poisson(lambda_E x eps_H x R_T)
```

R_T (trailers per boat group) is bounded between 0 and 1 and receives a `Beta(0.5, 0.5)` prior, evaluated only when trailer data is present (`T_n > 0` guard in the Stan model).

**Key difference from shore:** Boat crabbers have substantially higher CPUE than shore crabbers (~4.8 crab/trip vs ~1.1 crab/trip), reflecting the ability of boats to access deeper water and deploy gear more efficiently. This justifies running separate CPUE processes for each population.

**Data limitation:** The ring-net sub-season has only ~17-86 interviews depending on processing filters. When interview count falls below `bss_min_interviews` (default: 20), the framework falls back to PE for that stratum.

### 4.3 Commercial/Charter Vessels

**In plain language:** Large fishing vessels — commercial crab boats and charter boats — moored at Westport Marina that participate in recreational crabbing before the commercial season opens. These boats are much more efficient than private recreational boats and catch a disproportionate amount of crab relative to their numbers.

**Population includes:** All interviews with `boat_type` matching "Commercial Crab Vessel", "Charter", or "Guide".

**2024-25 data:**
- 47 tally days (Dec 3, 2024 – Feb 8, 2025)
- ~164 interviews (after filtering)
- Active period: Sep 16, 2024 through Jan 15, 2025 (commercial opener)

**Why not BSS?** The Marina has only 5 gear count records across 3 days — far too sparse for a state-space model. Vessel tallies are also inherently different from gear or trailer counts (counting discrete large vessels rather than individual pieces of equipment).

**Estimation method:** Direct census expansion:
```
Daily harvest = (commercial vessels + charter vessels) x mean crab per vessel
Season harvest = sum(daily harvest on sampled days) x (total possible days / sampled days)
```

The total possible days runs from the season start through the commercial opener.

---

## 5. Statistical Models

### 5.1 Point Estimator (PE)

**In plain language:** The PE method calculates a simple average of daily effort for each month x day-type combination, multiplies by the total number of days in that stratum to get total effort, and estimates catch by multiplying effort by average CPUE.

**Technical description:**

For each stratum defined by period (month) and day type (weekday/weekend/holiday):

1. **Daily effort:** On each sampled day, the mean gear count across count sequences is converted to estimated crabbers using the gear-per-crabber ratio, then multiplied by day length to get crabber-hours.

2. **Stratum total effort:**
```
E_stratum = mean(daily_crabber_hours) x N_total_days_in_stratum
SE_stratum = sqrt(N^2 x var(daily_crabber_hours) / n_sampled)
```

3. **CPUE:** Daily CPUE is calculated as total crab caught / total crabber-hours across all interviews on that day. Stratum-level CPUE is the interview-weighted mean of daily CPUEs.

4. **Catch:**
```
C_stratum = E_stratum x mean_CPUE_stratum
```

5. **Season total:** Sum across all strata.

The PE method is simple and transparent. Its main limitation is that it cannot estimate effort or catch on unsampled days except through the stratum-mean assumption — if a particular sampled weekend was unusually busy, the PE treats all weekends in that month as equally busy.

### 5.2 Bayesian State-Space Model (BSS)

**In plain language:** The BSS model treats daily effort and CPUE as hidden quantities that evolve smoothly over time, following an autoregressive process. On days with data, the model uses observations to pin down the hidden quantities. On days without data, it fills in the gaps using the temporal pattern learned from observed days. The result is a complete time series of daily effort and catch with uncertainty bounds.

**Technical description:**

The BSS model estimates two latent processes:

**Effort process:** Daily effort follows a log-normal AR(1) process with a weekend/holiday effect:

```
log(lambda_E[d]) = mu_E + omega_E[period[d]] + B1 x w[d]
```

Where:
- `mu_E` is the season-long effort intercept (log scale)
- `omega_E[p]` is the period-level residual: `omega_E[p] = phi_E x omega_E[p-1] + epsilon_E[p]`
- `B1` is the weekend/holiday effect on effort
- `phi_E` is the AR(1) coefficient (-1 to 1; positive values mean effort is autocorrelated across periods)

**CPUE process:** Daily CPUE follows a similar AR(1) structure without the weekend effect:

```
log(lambda_C[d]) = mu_C + omega_C[period[d]]
```

**Within-day overdispersion:** The instantaneous effort at count time `i` on day `d` is:

```
lambda_E_instantaneous[d,i] = lambda_E[d] x eps_E_H[d,i]
```

Where `eps_E_H ~ Gamma(r_E, r_E)` with mean 1. This accounts for the fact that a single snapshot of effort at one time of day may not represent the daily average. Multiple counts per day help inform `r_E`; with only one count per day, `r_E` is essentially unidentified and samples from its prior.

**Observation models** (how latent quantities connect to observed data):

| Data Type | Likelihood | Link to Latent Process |
|---|---|---|
| Gear counts (shore) | `Gear_I ~ Poisson(lambda_E x eps_H x R_G)` | Gear = crabbers x gear-per-crabber |
| Trailer counts (boat) | `T_I ~ Poisson(lambda_E x eps_H x R_T)` | Trailers = boat groups x trailers-per-group |
| Direct crabber counts (future) | `Crab_I ~ Poisson(lambda_E x eps_H x p_I)` | Direct observation with visibility fraction |
| Interview CPUE | `c ~ NegBin(lambda_C x h, r_C)` | Crab caught = CPUE x hours fished |
| Interview gear expansion | `Gear_A ~ Poisson(A_A x R_G)` | Group gear = group size x gear-per-crabber |
| Interview trailer expansion | `T_A ~ Bernoulli(R_T)` | Group has trailer with probability R_T |

**Generated quantities:** After fitting, the model generates realized daily catch and effort:

```
E[d] = lambda_E[d] x L[d]                              # daily crabber-hours
C[d] = Poisson_rng(lambda_E[d] x L[d] x lambda_C[d])   # daily catch (stochastic draw)
```

Season totals (`E_sum`, `C_sum`) are the sums across all days.

### 5.3 Stan Model Specification

The Stan model file is `BSS_crab_model_01.stan`. Its design is tailored for the crab fishery's observation structure:

**Effort observation types supported:**
- Gear counts (docks) with `R_G` expansion parameter — lognormal prior, unbounded >0 because crabbers deploy >1 gear unit each
- Trailer counts (boat launch) with `R_T` expansion parameter — beta prior, bounded 0-1, guarded by `if (T_n > 0)` so it is only evaluated when trailer data is present
- Direct crabber counts (jetty, future) with `p_I_crab` visibility parameter — passed as data, not estimated
- Overflow protection in generated quantities — Poisson_rng rates capped at 1e9 to prevent crashes during warmup

### 5.4 Period Stratification

The AR(1) process operates at the **period** level, not the daily level. The default period is **month**, giving P_n = 3 for the ring-net sub-season (Sep, Oct, Nov) and P_n = 10 for the all-gear sub-season (Dec through Sep).

Monthly periods balance temporal resolution (capturing seasonal trends) with computational tractability (fewer AR(1) states = faster fitting).

---

## 6. Convergence Tuning

The BSS model's convergence behavior varies across the four population x sub-season combinations, reflecting differences in data density and parameter identifiability. The framework implements fit-specific tuning to address these differences.

### 6.1 Convergence Diagnostics

| Fit | D (days) | Effort Obs | Interviews | Counts/Day | Primary Challenge |
|---|---|---|---|---|---|
| Shore ring-net | 76 | 52 | ~856 | Mostly 1 | Few multi-count days; `r_E` poorly identified |
| Shore all-gear | 289 | 313 | ~2,743 | Mix of 1 and 3 | Large parameter space (867 eps_E_H) |
| Boat ring-net | 76 | 44 | ~17-86 | 1 | Too few interviews for reliable CPUE |
| Boat all-gear | 289 | 196 | ~336 | 1 | Moderate data; trailer counts noisier than gear |

### 6.2 Higher Max Treedepth for Shore Fits

Shore fits have D x H `eps_E_H` gamma overdispersion parameters (228 for ring-net, 867 for all-gear). Most have no data informing them — only days with effort counts contribute likelihood. This creates a high-dimensional space with many flat directions that the NUTS sampler struggles to navigate within the default 2^10 = 1,024 leapfrog step budget.

Shore fits use `max_treedepth = 14`, allowing up to 2^14 = 16,384 leapfrog steps per transition. Boat fits retain `max_treedepth = 10` because their parameter space is smaller.

### 6.3 Guarded R_T Prior

The `R_T` parameter (trailers per boat group) receives a `Beta(0.5, 0.5)` prior, which is U-shaped with infinite density at 0 and 1. When no trailer data is present (shore-only fits where `T_n = 0`), the sampler encounters sharp curvature at these boundaries, causing divergent transitions.

The Stan model wraps the R_T prior in a conditional:
```stan
if (T_n > 0) {
  R_T ~ beta(0.5, 0.5);
}
```

When `T_n = 0`, R_T contributes nothing to the log-posterior and is effectively ignored.

### 6.4 PE Fallback for Sparse Strata

A minimum interview threshold is enforced:
```r
bss_min_interviews = 20
```

When a population x sub-season combination has fewer than 20 interviews, the BSS model is skipped and the PE estimate is used instead. Twenty interviews is roughly the minimum needed for a stratum-level CPUE estimate with reasonable precision; below this the BSS's latent CPUE process is dominated by its prior rather than data.

This primarily affects the boat ring-net sub-season (~17 interviews after processing).

### 6.5 Increased Iterations for Shore All-Gear

The shore all-gear fit (289 days) has the largest parameter space and requires more samples for adequate effective sample size (ESS) on the season-total summaries:

```r
bss_iter_shore_allgear = 4000    # (default for other fits: 2000)
bss_warmup_shore_allgear = 2000  # (default: 1000)
```

This produces 8,000 post-warmup draws. When combining posterior draws across fits for the port total, shorter draw vectors are recycled using `rep_len()`. This is statistically valid because draws from independent fits are independent samples from their respective posteriors.

---

## 7. Estimation Workflow

```
                    SET PARAMETERS
                  Season dates, pot open date,
                  commercial opener, holidays,
                  BSS tuning, catch group toggle
                         |
                         v
                  LOAD & CLASSIFY DATA
                  fetch_crab_data_v2()
                  -> shore effort (gear counts)
                  -> boat effort (trailer counts)
                  -> interviews (classified by pop)
                  -> catch, commercial tally
                         |
                         v
                   DIAGNOSTIC PLOTS
                  Gear & trailer timeseries,
                  CPUE by population,
                  monthly boxplot, vessel tally
                         |
            +------------+-------------+
            |            |             |
            v            v             v
         SHORE        BOAT      COMM/CHARTER
         PE+BSS       PE+BSS    Census/tally
         For each     For each
         sub-season:  sub-season:   Truncated at
          - PE          - PE        commercial
          - BSS*        - BSS*      opener date
            |            |             |
            v            v             v
                 COMBINE ESTIMATES
               Sum posterior draws
               PE fallback for sparse strata
               Comm/charter as constant
               -> Port total with uncertainty
                         |
                         v
                       OUTPUT
                  CSVs, plots, parameters

         * BSS skipped when < 20 interviews
```

### 7.1 BSS Fits

With `estimate_red_rock = FALSE` (Dungeness only):

| Fit | Population | Sub-season | Method | Iterations | Treedepth | Est. Runtime |
|---|---|---|---|---|---|---|
| 1 | Shore | Ring-net | BSS | 2000/1000 | 14 | ~30 min |
| 2 | Shore | All-gear | BSS | 4000/2000 | 14 | ~2 hr |
| 3 | Boat | Ring-net | **PE only** | — | — | instant |
| 4 | Boat | All-gear | BSS | 2000/1000 | 10 | ~30 min |
| 5 | Comm/Charter | Full season | Census | — | — | instant |

**Total estimated runtime: ~3 hours** on a 4-core machine (dominated by the Shore all-gear fit).

With `estimate_red_rock = TRUE`, fits 1-4 are duplicated for Red Rock Kept, approximately doubling total runtime.

---

## 8. Output Files

Each run produces the following files in `output/YYYYMMDD/`:

### Tables
| File | Contents |
|---|---|
| `pe_port_summary.csv` | PE estimates by component and port total |
| `port_total_Dungeness_Kept.csv` | Combined PE + BSS port total |
| `bss_summary_{label}.csv` | Stan summary (mean, SE, SD, quantiles, n_eff, R-hat) per BSS fit |
| `bss_daily_effort_{label}.csv` | Daily effort posteriors (median, 95% CI) per BSS fit |
| `bss_daily_catch_{label}.csv` | Daily catch posteriors per BSS fit |
| `run_parameters.txt` | Full parameter list for reproducibility |

### Plots
| File | Description |
|---|---|
| `plot_shore_effort_timeseries.png` | Gear counts over time, colored by day type, with pot-open reference line |
| `plot_boat_effort_timeseries.png` | Trailer counts over time |
| `plot_cpue_timeseries.png` | CPUE by population (faceted), with LOESS smoother |
| `plot_effort_by_month.png` | Monthly gear count boxplots by day type |
| `plot_commercial_tally.png` | Stacked bar chart of daily vessel tallies |
| `plot_pe_daily_effort.png` | PE daily crabber-hours by population (faceted) |
| `plot_bss_effort_{label}.png` | BSS daily effort ribbon (median + 95% CI) per fit |
| `plot_bss_catch_{label}.png` | BSS daily catch ribbon per fit |
| `plot_bss_posteriors_{label}.png` | Posterior density plots for season totals per fit |

---

## 9. Key Assumptions

### 9.1 Assumptions Common to PE and BSS

1. **Gear counts represent relative effort.** A gear count of 50 on day A and 25 on day B implies approximately twice as much crabbing effort on day A. This assumes gear deployment duration is roughly constant across days.

2. **The gear-per-crabber ratio (R_G) is constant across the sub-season.** In reality, this ratio may vary (e.g., experienced crabbers deploy more gear). The model estimates a single R_G from all interviews in the sub-season.

3. **Trailer counts represent boat crabber effort.** Not all trailers at the boat launch correspond to crabbers — some may be for finfish anglers. The classification relies on interview data to establish CPUE for boat-mode crabbers, which implicitly accounts for the mix.

4. **Day length is a reasonable proxy for available fishing hours.** The fixed day length serves as a scaling factor; any systematic bias is absorbed into the effort intercept.

5. **Interviews are representative of all crabbers in the population.** If creel clerks preferentially interview completed-trip crabbers (who tend to have higher catch), CPUE estimates may be biased upward. Trip status is recorded and can be used for sensitivity analysis.

6. **Crab mortality from released catch is not estimated.** Only kept crab are counted. Handling mortality of released crab (sub-legal, female) is not included.

### 9.2 Additional BSS Assumptions

7. **The AR(1) process captures temporal dynamics.** Effort and CPUE change smoothly between periods, with mean reversion. Sudden disruptions (storms, regulatory changes mid-season) appear as large process errors rather than being modeled explicitly.

8. **Within-day overdispersion is gamma-distributed.** The `eps_E_H` parameters assume effort within a day follows a gamma distribution around the daily mean.

9. **Poisson observation error.** Effort counts and catch counts are Poisson conditional on the latent rates, implying variance equals the mean.

### 9.3 Commercial/Charter Assumptions

10. **Mean catch per vessel is constant across sampled and unsampled days.** Given that the tally covers 47 of approximately 122 possible days (39%) spanning the full active period, this is reasonable for a first estimate.

11. **Commercial/charter recreational harvest is negligible after the commercial opener.**

---

## 10. Interpreting PE vs BSS Differences

It is normal and expected for PE and BSS estimates to differ. The 2024-25 proof-of-concept run produced:

| | PE | BSS Median | BSS 95% CI |
|---|---|---|---|
| Effort (crabber-hrs) | 150,232 | 151,586 | 141,175 – 163,566 |
| Dungeness Caught | 74,237 | 90,061 | 81,045 – 104,387 |

**Effort agreement (~1% difference)** validates that the observation models (gear counts -> crabbers, trailers -> boat groups) are functioning correctly.

**Catch disagreement (~21% higher for BSS)** is expected because:

1. **BSS interpolates unsampled days.** PE assigns zero catch to strata with no interview data. BSS uses the latent CPUE process to estimate catch on unsampled days.

2. **BSS accounts for weekend peaks.** High-effort weekend days may be undersampled. The BSS's weekend effect and temporal smoothing capture this.

3. **BSS propagates uncertainty.** The BSS 95% CI contains the PE estimate at its lower tail, indicating the methods are statistically compatible.

The BSS estimate should be preferred when convergence is adequate (R-hat < 1.05, ESS > 400, few divergences). When convergence is poor, the PE provides a robust fallback.

---

## 11. Known Limitations and Future Work

### 11.1 Current Limitations

- **No jetty effort counts.** The 28 jetty interviews contribute to shore CPUE but there is no independent effort data for the jetty.

- **Single-count days in ring-net sub-season.** The within-day overdispersion parameter `r_E` is poorly identified during this period.

- **Marina effort data is sparse.** The commercial/charter census relies entirely on the tally sheet.

- **No ingress/egress data.** Digital camera data from dock access points has been collected but not yet digitized.

- **Beach crabbing is unmeasured.** Only 1 beach interview in 2024-25.

### 11.2 Planned Improvements

1. **Digitize ingress/egress camera data** to create continuous daily crabber counts.

2. **Add Westport Jetty effort counts** using the `Crab_n`/`Crab_I` inputs already in the Stan model.

3. **Improve commercial/charter estimation** by stratifying the census expansion by day type.

4. **Restructure `eps_E_H` allocation** so only days with observations receive overdispersion parameters.

5. **Expand to other ports:** Willapa Bay (Tokeland, Chinook) and Columbia River.

6. **Add a holiday effect** as a separate covariate (`B2 x holiday[d]`) in the BSS effort process.

---

## 12. Reproducibility

To reproduce the 2024-25 estimates:

1. Clone the `FWC-estimation-method` repository.
2. Place input files in `input_files/`:
   - `effort_combined.csv` (re-exported with proper quoting)
   - `interview_combined.csv` (dates in M/D/YYYY; number_of_gear from column N)
   - `wes_commercial_tally.csv`
3. Place `BSS_crab_model_01.stan` in `stan_models/`.
4. Open `FWC-estimation-method-wesport-v3.Rmd`.
5. Verify parameters in Section 0.1 match the desired settings.
6. Knit or run all chunks.

**Software requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, here, patchwork, gt, suncalc.

**Expected runtime:** ~3 hours on a 4-core machine (dominated by Shore all-gear BSS fit with 4,000 iterations).

**Random seed:** Interview subsampling uses `set.seed(42)` for reproducibility. Stan's internal RNG seed is not fixed; results will vary slightly between runs due to MCMC stochasticity. Port-total estimates should be stable to within ~2-3% across runs when convergence is adequate.
