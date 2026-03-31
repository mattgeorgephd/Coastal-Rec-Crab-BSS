# Recreational Crab Harvest Estimation — Grays Harbor
## Gear-Resolved CPUE Model Documentation

**Authors:** Matt George, with analytical development support  
**Agency:** Washington Department of Fish and Wildlife (WDFW)  
**Date:** March 2026  
**Status:** Proof of concept — 2024-25 season  

---

## 1. Summary for Decision-Makers

This framework estimates the total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area, broken down by gear type (pots, ring nets, traps, snares). Unlike the Pooled CPUE Model which treats all gear types as a single group, this model estimates a separate catch rate for each gear type, producing gear-specific harvest estimates with their own uncertainty bounds.

**Why this matters for management:** Gear regulations are a primary management tool for the crab fishery. If managers are considering restricting or expanding pot use, they need to know how much harvest each gear type contributes — and how confident we are in that number. This model provides that directly, rather than relying on approximate proportional allocations.

**Key results (2024-25, through Jan 30):**
- Estimated total harvest: ~32,000–34,000 Dungeness crab (similar to pooled model)
- Each gear type now has its own posterior distribution, so we can say things like "Pots contributed 18,000 crab (95% CI: 15,000–22,000)" rather than just "Pots were ~57% of the total"
- Pot CPUE is estimated at ~X crab/crabber-hr; Ring Net at ~Y; these rates can be tracked over the season to detect depletion or seasonal patterns

**Trade-off:** This model takes ~50% longer to run and has more parameters to monitor for convergence. The total harvest estimate is similar to the pooled model; the advantage is in the gear-type decomposition.

---

## 2. Understanding the Two Estimation Methods

### Point Estimator (PE) — The Simple Average

The PE method averages sampled days within each stat-week × day-type combination and expands to unsampled days. Fast and transparent, but can't fill temporal gaps or provide uncertainty bounds.

### Bayesian State-Space Model (BSS) — The Time-Series Model

The BSS fits a smooth curve through daily effort and catch rate data, estimating every day including unsampled ones. In this model, each gear type gets its own catch rate curve, so pot CPUE and ring net CPUE can follow different seasonal trajectories.

### Combined Best Estimate

Uses BSS when convergence diagnostics pass (R-hat < 1.05, n_eff > 400); otherwise falls back to PE.

---

## 3. Glossary

| Term | Meaning |
|---|---|
| **BSS** | Bayesian State-Space model — a time-series approach that estimates hidden daily quantities |
| **PE** | Point Estimator — a simpler stratified-expansion approach |
| **CPUE** | Catch Per Unit Effort — the average number of crab caught per crabber-hour |
| **Crabber-hour** | One person crabbing for one hour (the basic unit of effort) |
| **Gear count** | A snapshot count of crab gear deployed in the water at a point in time |
| **R_G** | Gear-per-crabber ratio — how many pieces of gear each crabber deploys (~1.27) |
| **R_T** | Trailer-per-boat-group ratio — fraction of boat groups that have a trailer at the ramp |
| **AR(1)** | First-order autoregressive process — a model where today's value depends on yesterday's |
| **Credible interval** | The Bayesian equivalent of a confidence interval — a range containing the true value with 95% probability |
| **R-hat** | A convergence diagnostic — should be below 1.05 for reliable results |
| **n_eff** | Effective sample size — how many independent samples the MCMC chain produced |
| **G_gear** | Number of gear types modeled in a sub-season (3 for ring-net, 4 for all-gear) |
| **pi_gear** | The proportion of crabbers using each gear type in a given period |
| **lambda_C_gear** | Gear-type-specific catch rate (crab per crabber-hour for a given gear type) |
| **r_C_gear** | Per-gear-type overdispersion — how variable catch is around the mean for each gear type |
| **Dirichlet** | A probability distribution over proportions that must sum to 1 |
| **Sub-season** | One of two periods: ring-net only (Sep–Nov) or all-gear (Dec–Sep) |

---

## 4. Study Area and Data Sources

### 4.1 Westport / Grays Harbor

Westport is a small coastal town on the south side of the Grays Harbor estuary in Washington State, and one of the highest-volume recreational crabbing ports on the Pacific coast.

### 4.2 Effort Data

| Site | What Is Counted | 2024-25 Records |
|---|---|---|
| Westport Docks Float 20 | Crab gear in the water | 348 |
| Westport Docks Float 17-21 | Crab gear in the water | 169 |
| Westport Boat Launch | Boat trailers at ramp | 226 |
| Ocean Shores Boat Launch | Boat trailers at ramp | 4 |
| Westport Marina | Vessels (limited) | 5 |

### 4.3 Interview Data

4,359 interviews conducted across all Grays Harbor sites. Each classified into Shore, Private Boat, or Commercial/Charter based on crabbing mode and boat type.

**Gear type classification:** Each interview is assigned to a primary gear type using a priority hierarchy: **Pot > Ring Net > Trap > Snare**. About 30% of interviews list multiple gear types (e.g., "Pot, Ring Net"); these are assigned to the highest-priority type. This is a simplification — future work could split catch proportionally across reported gear types.

**Interview counts by gear type (2024-25):**

| Gear Type | Interviews | Notes |
|---|---|---|
| Pot | ~683 | Dominant in all-gear sub-season; illegal in ring-net |
| Ring Net | ~676 | Used year-round |
| Trap (foldable/star) | ~909 | Most common overall |
| Snare | ~302 | Lower CPUE |

### 4.4 Commercial/Charter Vessel Tally

Daily vessel tally maintained Dec 1, 2024 – Feb 8, 2025 (47 tally days within a 70-day census period).

---

## 5. Season Structure

Same as the pooled model: two sub-seasons split at the pot-open date (Dec 1).

| Sub-season | Dates | Duration | Gear Types Modeled |
|---|---|---|---|
| Ring-net only | Sep 16 – Nov 30 | 76 days | 3 (Ring Net, Trap, Snare) |
| All-gear | Dec 1 – Sep 15 | 289 days | 4 (Pot, Ring Net, Trap, Snare) |

The number of gear types modeled (G_gear) adapts to what's available. Since pots are illegal during the ring-net sub-season, the model fits 3 gear types for that period and 4 for the all-gear period.

### Day Length

Computed daily from civil twilight (dawn to dusk) via `suncalc`, capped 9–16 hours.

### Day Type and Holiday Effect

Each day is weekday, weekend, or holiday. The effort equation includes two separate effects: B1 (weekend boost) and B2 (additional holiday boost), so holidays are not just treated as "another weekend" — they can have their own, typically larger, effort multiplier.

---

## 6. Population Components

Three populations estimated independently and summed:

### 6.1 Shore Crabbers
Dock + jetty + beach crabbers. Effort from gear counts. Largest group (~3,600 interviews).

### 6.2 Private Boat Crabbers
Trailered boats from boat launch. Effort from trailer counts. Higher CPUE than shore (~4.8 vs ~1.1 crab/trip).

### 6.3 Commercial/Charter Vessels
Moored at Westport Marina. Day-type stratified census expansion from the vessel tally. Weekday and weekend/holiday harvest rates are computed separately and expanded to the full census period, preventing bias from uneven sampling.

---

## 7. The Gear-Resolved CPUE Model

### 7.1 How It Works (Non-Technical)

This model is identical to the Pooled CPUE Model on the effort side — it estimates the same total number of crabbers present each day. The difference is on the catch side: instead of one catch rate for all gear types combined, it estimates a separate catch rate for each gear type.

Think of it this way: the model knows that on a given day, there are 100 crabbers present. Of those 100, about 40 are using pots, 25 ring nets, 25 traps, and 10 snares (these proportions come from interview data and are allowed to change over the season). Each gear type has its own catch rate — pot users catch more per hour than ring net users. The daily catch for each gear type is: (number of crabbers using that gear) × (hours available) × (catch rate for that gear type).

Because the model estimates each gear type's catch rate as a time series with uncertainty, the resulting gear-type harvest estimates have proper credible intervals — not just proportional allocations from interviews.

### 7.2 How It Works (Technical)

**Effort process (shared, same as pooled):**
```
log(lambda_E[d]) = mu_E + omega_E[period[d]] + B1 × w[d] + B2 × holiday[d]
```

**Gear-type CPUE processes (one per gear type):**
```
log(lambda_C_gear[d, g]) = mu_C_gear[g] + omega_C_gear[period[d], g]
```
In words: each gear type `g` has its own baseline catch rate (`mu_C_gear[g]`) and its own smooth temporal deviation (`omega_C_gear`). All gear types share the same AR(1) coefficient (`phi_C_gear`) and process error SD (`sigma_eps_C_gear`), but have independent trajectories.

**Gear-type proportions:**
```
pi_gear[period] ~ Dirichlet(alpha)
```
In words: the fraction of crabbers using each gear type varies across periods and is estimated from a Dirichlet distribution. The prior concentration `alpha` is weakly informative, centered on observed interview proportions.

**Observation models:**
- Effort observations: same as pooled (gear counts → Poisson with R_G expansion)
- Interview gear assignment: `gear_type ~ Categorical(pi_gear[period])` — each interview's gear type informs the proportions
- Interview catch: `crab_caught ~ NegBin(lambda_C_gear[gear_type] × hours, r_C_gear[gear_type])` — catch depends on the gear-specific CPUE and the gear-specific overdispersion

**Per-gear-type overdispersion:** Each gear type has its own negative binomial overdispersion parameter `r_C_gear[g]`. This allows pot catches to have a different variance structure than ring net catches, which is biologically plausible — pots tend to produce more consistent catches while ring net catches are more variable.

**Daily catch by gear type:**
```
C_gear[d, g] = Poisson_draw(crabbers[d] × day_length[d] × pi_gear[g] × lambda_C_gear[d, g])
```
In words: daily catch for gear type g equals crabbers present × hours × fraction using that gear × catch rate for that gear.

### 7.3 Stan Model File

The Stan model is `crab_bss_gear_resolved.stan`. Key parameters beyond the pooled model:

| Parameter | Description | Prior |
|---|---|---|
| mu_C_gear[G_gear] | CPUE intercept per gear type (log scale) | Normal(log(0.5), 2) |
| omega_C_gear[P_n, G_gear] | CPUE AR(1) residuals per gear type | AR(1) with shared phi_C_gear |
| phi_C_gear | Shared CPUE AR(1) coefficient | Beta(2,2) rescaled |
| sigma_eps_C_gear | Shared CPUE process error SD | Half-Cauchy(5) |
| r_C_gear[G_gear] | NegBin overdispersion per gear type | Half-Cauchy(5) |
| pi_gear[P_n] | Gear-type proportions per period | Dirichlet(alpha) |
| B1 | Weekend effort multiplier | Normal(0, 1) |
| B2 | Additional holiday effort multiplier | Normal(0, 1) |

---

## 8. Convergence and Model Selection

Same as pooled: R-hat < 1.05 AND n_eff > 400 for C_sum and E_sum. The gear-resolved model has more parameters and may require closer monitoring — watch `C_sum_gear` convergence per gear type.

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
| `catch_by_gear_type_detail.csv` | **Gear-type catch with BSS posterior uncertainty (median + 95% CI)** |
| `catch_by_gear_type.csv` | Port-level catch by gear type |
| `monthly_by_area.csv` | Monthly catch by creel area |
| `wide_summary.csv` | All dimensions in one table |
| `bss_summary_{label}.csv` | Stan convergence diagnostics including `C_sum_gear` and `r_C_gear_out` |
| `run_parameters.txt` | Full parameter list for reproducibility |

---

## 10. Key Assumptions

### Shared with Pooled Model (1–8)
1. Gear counts represent relative effort across days
2. Gear-per-crabber ratio is constant within a sub-season
3. Trailer counts represent boat crabber effort
4. Day length from civil twilight is a reasonable proxy
5. Interviews are representative of the population
6. Released crab mortality is not estimated
7. AR(1) captures temporal dynamics
8. Within-day overdispersion is gamma-distributed

### Gear-Type Specific (9–13)
9. **Multi-gear interviews are assigned to a single primary type** using a priority hierarchy (Pot > Ring Net > Trap > Snare). This may undercount secondary gear types.
10. **All gear types share the same AR(1) autocorrelation coefficient.** The temporal persistence of CPUE is assumed similar across gear types.
11. **Gear-type proportions vary by period but not by day type.** If the gear mix differs on weekends vs weekdays, this is not captured.
12. **Each gear type has its own overdispersion parameter.** This allows different variance structures across gear types — pots may produce more consistent catches than ring nets.
13. **The Dirichlet prior is weakly informative**, centered on observed interview proportions × 10. The data dominates the prior when interview counts are large.

### Commercial/Charter (14–15)
14. Mean catch per vessel is constant within each day-type stratum
15. Commercial/charter recreational harvest is negligible after the census period

---

## 11. Limitations and Planned Improvements

### Current Limitations
- No jetty effort counts
- Single-count days in ring-net sub-season
- Beach crabbing unmeasured
- Multi-gear interview assignment is deterministic (Pot > Ring Net > Trap > Snare)
- Shared `phi_C_gear` assumes similar autocorrelation across gear types

### Planned Improvements
1. Add Westport Jetty effort counts
2. Restructure `eps_E_H` for data-informed days only
3. Expand to other ports
4. Independent `phi_C_gear` per gear type if data supports it
5. Full effort decomposition by gear type (separate effort processes per gear)
6. Proportional multi-gear interview assignment (split catch across reported gear types)

---

## 12. Reproducibility

1. Clone the `FWC-estimation-method` repository
2. Place input files in `input_files/`
3. Place `crab_bss_gear_resolved.stan` in `stan_models/`
4. Open `BSS-GH-gear-type-CPUE-model.Rmd`
5. Set parameters and run

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here  
**Expected runtime:** ~4-5 hours on a 4-core machine
