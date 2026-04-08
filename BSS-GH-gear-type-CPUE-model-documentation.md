# Recreational Dungeness Crab Harvest Estimation — Grays Harbor
## Gear-Resolved CPUE Model: Technical Documentation

**Author:** Matthew George, Ph.D.  
**Contact:** matthew.george@dfw.wa.gov  
**Agency:** Washington Department of Fish and Wildlife (WDFW)  
**Status:** Operational — Annual estimation framework  
**Model version:** v5.2 (Stan model v3.2)

---

## 1. Summary for Decision-Makers

This framework estimates the total recreational Dungeness crab harvest at Westport and the greater Grays Harbor area, broken down by gear type (pots, ring nets, traps, snares). Unlike the companion Pooled CPUE Model which treats all gear types as a single group, this model estimates a separate catch rate for each gear type, producing gear-specific harvest numbers with their own uncertainty bounds.

**Why this matters for management:** Gear regulations are a primary management tool for the crab fishery. If managers are considering restricting pot use, extending the ring-net-only period, or evaluating snare effectiveness, they need to know how much harvest each gear type contributes — and how confident we are in that number. This model provides that directly.

**What this model produces:**

- A total Dungeness crab harvest estimate with a 95% credible interval, identical in structure to the Pooled model.
- Per-gear-type harvest estimates with their own posterior uncertainty (e.g., "Pots contributed 18,000 crab, 95% CI: 15,000–22,000"). These credible intervals account for uncertainty in effort, catch rate, and gear-type proportions.
- Gear-type CPUE trajectories over the season — enabling questions like "is pot CPUE declining through summer?" or "how does ring net efficiency compare to traps?"
- Time-varying gear-type proportions showing how the mix of gear in use shifts across the season and between day types.
- Monthly harvest trends, mode breakdowns, and all standard outputs.
- A sensitivity analysis quantifying how incomplete trip filtering affects harvest estimates by population, sub-season, and gear type — providing transparency about data quality decisions.

**Trade-off:** This model takes ~50% longer to run and has more parameters to monitor for convergence. The total harvest estimate is similar to the Pooled model; the advantage is in the gear-type decomposition with formal uncertainty.

---

## 2. Understanding the Two Estimation Methods

### 2.1 Point Estimator (PE) — The Simple Average

The PE method computes the average daily harvest for each week × day-type (weekday vs weekend) combination from sampled days, then expands to unsampled days by multiplying by the total number of days in that combination.

**Strengths:** Easy to understand. Transparent calculations. No modeling assumptions beyond "sampled days represent unsampled days within a stratum."

**Weaknesses:** Biased if sampled days are unrepresentative. Cannot estimate for unsampled time periods. No uncertainty bounds. Treats each stratum independently.

### 2.2 Bayesian State-Space Model (BSS) — The Time-Series Model

The BSS fits a smooth curve through daily effort and catch rate data, estimating every day including unsampled ones. In this model, each gear type gets its own catch rate curve with its own temporal variability, so pot CPUE and ring net CPUE can follow different seasonal trajectories with different degrees of smoothness. The model produces full posterior distributions for all quantities, enabling rigorous uncertainty quantification.

**Strengths:** Fills temporal gaps. Produces uncertainty bounds. Gear-type decomposition with posterior uncertainty. Accounts for temporal autocorrelation. Per-gear-type process variability.

**Weaknesses:** More complex. Longer runtime. Must be checked for convergence.

### 2.3 Combined Best Estimate

BSS is preferred when all convergence diagnostics pass: R-hat < 1.05 AND n_eff > 400 for all monitored parameters (including per-gear-type quantities), AND no divergent transitions during sampling. Divergent transitions indicate that the sampler has encountered regions of the posterior with difficult geometry, which can bias parameter estimates even when R-hat and n_eff appear satisfactory. When any diagnostic fails, PE is used as the fallback. The `convergence_report.csv` output documents the decision and all diagnostic values for every fit.

---

## 3. The Recreational Crab Fishery at Grays Harbor

### 3.1 Fishery Overview

The recreational Dungeness crab (*Metacarcinus magister*) fishery at Westport and the greater Grays Harbor area is one of the highest-volume recreational crabbing operations on the Washington coast. Thousands of crabbers participate annually using four primary gear types: crab pots, ring nets, foldable/star traps, and snares. The fishery operates year-round, though effort peaks during summer months and around major holidays.

WDFW regulations restrict pot use to a specific portion of the season (typically December through September), creating a structural break in effort and catch rates. This gear restriction is the primary reason the season is split into two independent sub-seasons for estimation.

Commercial Dungeness crab vessels also participate in the recreational fishery before the commercial season opens. These large vessels crab recreationally under the same daily limits as private boats but tend to have much higher catch rates per vessel. Their harvest is tracked separately.

### 3.2 Why Estimate Harvest by Gear Type?

Beyond the total harvest estimate, gear-specific information supports several management needs: evaluating whether proposed gear restrictions would meaningfully reduce harvest, monitoring whether specific gear types are becoming more or less efficient over time (potentially indicating changes in crab abundance or behavior), and understanding how the recreational fleet adapts when gear regulations change.

### 3.3 Study Area

Westport is a small coastal town on the south side of the Grays Harbor estuary. Recreational crabbing occurs from public docks (Floats 17–21), a jetty, beaches, a public boat launch, and a commercial marina. Each access point serves a different crabbing mode and requires a different type of effort measurement.

---

## 4. Data Sources and Field Collection

The framework requires three types of field data collected by WDFW field staff using Apple iPads running the WDFW iForm "Crab Creel Survey" application.

### 4.1 Effort Counts

**What they are:** Instantaneous point-in-time counts of gear or trailers at specific sites.

**Why they matter:** Effort counts measure total crabbing activity independently of interviews. The model combines effort counts with interview catch rates to estimate total harvest.

**How they are collected:** Field staff visit sites and record the number of crabbing indicators visible. Protocol calls for multiple counts per day at standardized times. The number of within-day counts affects the model's ability to estimate within-day variability.

**Sites:**

| Site | What Is Counted | Role in Model |
|---|---|---|
| Westport Docks Float 20 | Crab gear in the water | Primary shore effort indicator |
| Westport Docks Float 17-21 | Crab gear in the water | Summed with Float 20 for section total |
| Westport Boat Launch | Boat trailers at ramp | Private boat effort indicator |
| Ocean Shores Boat Launch | Boat trailers at ramp | Supplementary boat effort |
| Westport Jetty | Crabbers (future) | Reserved for future use |

Float 20 and Float 17-21 counts are paired by time and summed because dock crabbers distribute across both float areas.

**Assumptions:**
- A point-in-time count reflects relative effort on that day.
- Time-of-day introduces noise but not systematic bias.
- The surveyor's count is accurate.

### 4.2 Crabber Interviews

**What they are:** Dockside interviews recording trip-level information from crabbers.

**Why they matter:** Interviews provide the catch rate (crab per crabber-hour), the gear-per-crabber ratio (R_G), and — critically for this model — the gear type(s) each crabber used. The gear-type field is the foundation of the gear-resolved CPUE decomposition.

**What each interview records:** Number of crabbers, number of gear units, gear type(s) (select all that apply from: pot, ring net, trap, snare), hours fished, crabber-hours, Dungeness crab kept, Red Rock crab kept, trip completion status, crabbing mode (Dock, Boat, Jetty, Beach), and boat type (Private, Commercial, Charter, Guide).

**Gear-type recording:** The iForm allows multiple gear type selections per interview. About 28% of valid Grays Harbor interviews report multiple gear types (e.g., "Pot, Ring Net"). The model handles this through weighted fractional assignment across reported gear types (see Section 7).

**Population classification:**

| Crabbing Mode | Boat Type | Population |
|---|---|---|
| Dock, Jetty, or Beach | Any | **Shore** |
| Boat | Private or blank | **Private Boat** |
| Boat | Commercial, Charter, or Guide | **Commercial/Charter** |

#### 4.2.1 Trip Completion Filtering

Field staff record whether the crabber has completed their trip at the time of interview. Only completed-trip interviews are used for CPUE estimation. This filtering is controlled by the `filter_incomplete_trips` parameter (default: `TRUE`).

**Why this matters:** Crabbers interviewed mid-trip have systematically lower catch than they would at the end of the trip, because their gear has not finished soaking. In the 2024-25 Grays Harbor data, approximately 35% of shore interviews are incomplete trips, and these show a mean CPUE that is 20% lower than completed trips. This downward bias is not uniform across gear types — it is largest for soak-time-dependent gear (pots: -21%, traps: -23%) and negligible for ring nets (+4%), which are checked and pulled frequently throughout a trip. Including incomplete trips without adjustment would underestimate total harvest by approximately 7% for shore crabbers and, more consequentially for this model, would differentially suppress pot and trap catch rates relative to ring net, distorting the gear-type harvest decomposition.

**How it works:** Interviews with `trip_status == "Incomplete"` are excluded from both BSS and PE CPUE estimation. Interviews with missing trip status are retained, as they may represent completed trips with incomplete metadata. The framework logs the number of interviews filtered, the CPUE difference between complete and incomplete groups, and produces a sensitivity analysis (`sensitivity_incomplete_trips.csv`) comparing harvest estimates with and without the filter for each population × sub-season combination.

**Effort counts are not affected.** The incomplete trip filter applies only to interviews used for catch-rate estimation, not to effort counts (gear counts, trailer counts), which measure activity at a point in time regardless of trip status.

**Assumptions:**
- Interviewed crabbers are representative of all crabbers in their population, conditional on trip completion.
- Crabbers accurately report catch, hours, group size, and gear types used.
- The gear type(s) reported reflect what was actually deployed.
- Completed trips provide an unbiased estimate of the average crab-per-crabber-hour achieved by crabbers on that day. Incomplete trips represent a left-censored observation of catch and are excluded to avoid systematic downward bias.

### 4.3 Commercial/Charter Vessel Tally

**What it is:** Daily counts of commercial crab vessels and charter boats at Westport Marina during the period when these vessels participate in the recreational fishery.

**Why it matters:** Commercial vessels have much higher per-vessel catch and are difficult to interview comprehensively. The tally provides a census-like measure of vessel-days for stratified expansion.

**Assumptions:**
- The tally captures all participating vessels on each sampled day.
- Mean harvest per vessel from interviews applies to uninterviewed vessels.
- The tally is stratified by day type to prevent sampling bias.

### 4.4 Input Files

| File | Contents | Notes |
|---|---|---|
| `effort_combined.csv` | Shore gear counts and boat trailer counts | Re-exported with QUOTE_ALL |
| `interview_combined.csv` | Crabber interviews with gear type field | Dates in M/D/YYYY |
| `wes_commercial_tally.csv` | Daily vessel counts | One row per tally day |

---

## 5. Season Structure

### 5.1 Sub-Season Split

| Sub-season | Typical Dates | Gear Allowed | Gear Types Modeled | BSS Period Type |
|---|---|---|---|---|
| Ring-net only | Sep 16 – Nov 30 | Ring nets, snares, foldable traps | 3 (Ring Net, Trap, Snare) | Biweekly |
| All-gear | Dec 1 – Sep 15 | All gear including crab pots | 4 (Pot, Ring Net, Trap, Snare) | Monthly |

Each sub-season is estimated independently. The number of gear types modeled (`G_gear`) is determined by two filters applied in sequence: regulatory exclusions and a minimum data threshold.

**Regulatory exclusion.** Each sub-season definition carries a `gear_exclude` list specifying gear types that are prohibited by regulation during that period. For the ring-net sub-season, pots are excluded. This is enforced structurally — regardless of what interview data contains, the model will not attempt to fit a CPUE process for a prohibited gear type. Any interviews that mention an excluded gear type (typically multi-gear crabbers naming gear they own rather than what they deployed) have their fractional weight redistributed to the remaining gear types they reported. For example, an interview recorded as "Pot, Ring Net" during the ring-net sub-season contributes 100% of its weight to Ring Net rather than 50% to each.

This exclusion is necessary because even a small number of phantom gear mentions (5–6 fractional interviews in a 76-day sub-season) can cause catastrophic model failure. With near-zero data to constrain the AR(1) CPUE process, the sampler produces hundreds of divergent transitions, extreme parameter values, and overflow in the catch generated quantities. The regulatory exclusion prevents this by encoding known fishery structure rather than relying solely on a sample-size heuristic.

**Minimum effective-N threshold.** After regulatory exclusions, each remaining gear type must have at least `bss_min_gear_effective_n` effective interviews (default: 15) — computed as the sum of fractional gear weights — to be modeled independently. Gear types below this threshold are collapsed. The threshold is set conservatively because fractional weights from multi-gear interviews can inflate the effective sample count beyond what the data truly supports for independent CPUE estimation. If fewer than two gear types qualify, the model falls back to a single "All" category.

The ring-net sub-season uses biweekly BSS periods instead of monthly. With only ~76 days and 2.5 calendar months, monthly periods provide too few AR(1) transitions for effective temporal smoothing. Biweekly periods yield ~5–6 periods, giving the AR(1) process more room to capture temporal variation.

**Assumption:** The pot-open date creates a clean structural break.

### 5.2 Day Length

Computed daily from civil twilight (dawn to dusk) via `suncalc` at Westport coordinates (46.904°N, 124.105°W), capped 9–17 hours.

At 46.9°N latitude, civil twilight reaches approximately 17 hours at the summer solstice. The 17-hour cap ensures the available crabbing window is not underestimated during peak summer months.

**Assumption:** Civil twilight is a reasonable proxy for the crabbing window.

### 5.3 Day Type and Holiday Effect

Each day is classified as weekday (Mon–Thu), weekend (Fri–Sun), or holiday. The BSS model includes two separate effort effects: B1 (weekend boost) and B2 (additional holiday boost beyond B1). This captures the pattern that holidays generate substantially higher effort than regular weekends.

Gear-type proportions vary by day type (weekday/weekend/holiday) in addition to period, allowing the model to capture differences in gear mix between weekdays and weekends (see Section 7.5).

**Assumption:** Day-type classification is the same across all populations. Effort effects are estimated from data, not fixed.

---

## 6. Population Components

Three populations estimated independently and summed for the port total.

### 6.1 Shore Crabbers (Dock + Jetty + Beach)

Effort from gear counts. Conversion: `(Gear counted ÷ R_G) × day_length = crabber-hours`. Largest population by interview count. R_G (gear per crabber) is estimated from interview data.

**Assumptions:** Dock gear count is proportional to total shore effort. R_G is constant within a sub-season. Jetty/beach effort captured via dock count as spatial proxy.

### 6.2 Private Boat Crabbers

Effort from trailer counts. Boat crabbers have 2–5× higher CPUE than shore crabbers.

R_G is fixed at 1.3 for boat fits (not estimated) because no gear-count data informs this parameter for boats. Allowing R_G to float at its prior in boat fits would waste parameter space and potentially interact with other parameters.

**R_T assumption (trailers per group):** All boat-launched crabbing groups are assumed to have exactly one trailer. This drives R_T toward 1.0, which is the expected value for a boat launch. The `T_A_int` variable is set to 1 for all boat interviews, reflecting this structural assumption rather than missing data.

**Assumptions:** Every crabbing boat has a trailer. Trailer counts capture all private boat crabbing. Boat interviews are representative.

### 6.3 Commercial/Charter Vessels

Day-type stratified census expansion from the vessel tally. Weekday and weekend harvest rates expanded independently.

**Assumptions:** Mean catch per vessel is constant within each day-type stratum. All harvest occurs during the census period.

---

## 7. The Gear-Resolved CPUE Model

### 7.1 How It Works (Non-Technical)

This model estimates the same total number of crabbers present each day as the Pooled model — the effort side is identical. The difference is on the catch side: instead of one catch rate for all gear types combined, it estimates a separate catch rate for each gear type, each with its own degree of temporal variability.

Think of it this way: the model knows that on a given day there are 100 crabbers present. Of those 100, about 40 are using pots, 25 ring nets, 25 traps, and 10 snares (these proportions come from interview data and vary by period and day type). Each gear type has its own catch rate — pot users catch more per hour than snare users — and the catch rate for each type evolves at its own pace. The daily catch for each gear type is: (crabbers using that gear) × (hours available) × (catch rate for that gear type).

Because the model estimates each gear type's catch rate as a time series with uncertainty, the gear-type harvest estimates have proper credible intervals — not just proportional allocations.

### 7.2 Weighted Gear-Type Classification of Interviews

Before gear-type classification, interviews are filtered for quality: only completed trips with at least 0.5 crabber-hours of fishing time are included (see Section 4.2.1). Gear types prohibited by regulation during the sub-season are then excluded: any mention of a prohibited type (e.g., "Pot" during the ring-net sub-season) is removed from the interview's gear list, and the fractional weight is redistributed to the remaining types (see Section 5.1). The gear-type classification and weighting described below operate on this filtered and regulation-adjusted set.

Each interview contributes to gear-type CPUE estimation in proportion to how many *eligible* gear types the crabber reported. Gear types are detected using word-boundary regular expressions (`\bpot\b`, `\bring\s*net\b`, `\b(trap|star)\b`, `\bsnare\b`) for robust matching.

For single-gear interviews (~70%), the interview contributes 100% to that gear type's CPUE. For multi-gear interviews (~30%), the interview's catch is split equally across the reported gear types.

**Example:** A crabber reports "Pot, Ring Net" with 6 crab kept during the all-gear sub-season. This interview contributes to the Pot CPUE likelihood with weight 0.5 and to the Ring Net CPUE likelihood with weight 0.5. If the same interview occurred during the ring-net sub-season (when pots are prohibited), the Pot mention would be excluded and the interview would contribute 100% to Ring Net.

**Why weighted assignment?** A naive approach would assign each multi-gear interview entirely to a single gear type (e.g., the first listed or the one assumed to have the highest catch rate). This would systematically inflate CPUE for the favored type and deflate it for others. With approximately 28% of valid interviews reporting multiple gear types, such a bias would meaningfully distort the gear-type harvest decomposition. Equal weighting across reported gear types is a reasonable approximation when per-gear catch totals are not recorded.

**Assumptions:**
- Equal weighting is a reasonable proxy when per-gear-type catch is not recorded. (If a crabber used 2 pots and 4 ring nets, the true split would be ~33%/67%, but we don't have per-type catch data.)
- At least 15 effective interviews (sum of fractional weights ≥ 15) per gear type are needed for independent modeling. Below this threshold, the gear type is collapsed. This conservative threshold accounts for the inflation of effective sample size by fractional weights from multi-gear interviews.
- Interviews with unrecognized gear types are assigned to the most common gear type.

### 7.3 Effort Process (Shared)

```
log(lambda_E[d]) = mu_E + omega_E[period[d]] + B1 × w[d] + B2 × holiday[d]
```

Identical to the Pooled model. AR(1) temporal evolution with weekend and holiday effects.

**Assumptions:** Effort evolves log-normally. The AR(1) coefficient captures persistence. Weekend and holiday effects are multiplicative.

### 7.4 Gear-Type CPUE Processes

```
log(lambda_C_gear[d, g]) = mu_C_gear[g] + omega_C_gear[period[d], g]
```

Each gear type `g` has its own baseline catch rate (`mu_C_gear[g]`) and its own AR(1) temporal deviation (`omega_C_gear`).

Each gear type has its own innovation standard deviation (`sigma_eps_C_gear[g]`), allowing pot CPUE to evolve smoothly while snare CPUE can be more volatile — consistent with the biological expectation that pots (which soak and integrate over time) produce more stable catch rates than snares (which are highly dependent on tide, weather, and skill). All gear types share the same AR(1) coefficient (`phi_C_gear`).

**Assumptions:**
- Each gear type's CPUE evolves smoothly and independently.
- All gear types share the same degree of temporal autocorrelation (phi_C_gear) but have independent process variability (sigma_eps_C_gear[g]).
- CPUE does not vary by day type within a period for any gear type.

### 7.5 Gear-Type Proportions and Their Uncertainty

The fraction of crabbers using each gear type (`pi_gear`) is computed from interview data as empirical proportions, stratified by period and day type (weekday/weekend/holiday). Laplace smoothing (alpha = 1) is applied to prevent zero proportions and provide gentle regularization.

```
pi_gear[period, day_type, g] = (n_weighted[period, day_type, g] + 1) / 
                                (N_weighted[period, day_type] + G_gear)
```

When a period × day_type combination has no interviews, the model falls back to the period-level proportions, then to the sub-season-level proportions.

These point-estimate proportions are used in the CPUE likelihood (model block) to allocate crabbers across gear types. However, `pi_gear` carries genuine sampling uncertainty — especially in sparse strata where a period × day_type cell may contain only a handful of interviews, or where the fallback hierarchy has been triggered. Ignoring this uncertainty would produce gear-type harvest credible intervals that are too narrow, because they would reflect only effort and CPUE uncertainty while treating the gear-type allocation as perfectly known.

**How gear-proportion uncertainty is propagated:** In the generated quantities block, where daily catches are computed for posterior summaries, the model draws gear-type proportions from a Dirichlet posterior on each MCMC iteration. The concentration parameters are the raw weighted interview counts (`n_weighted_gear`) plus the Laplace smoothing constant (alpha = 1), following the conjugate Dirichlet-multinomial relationship:

```
pi_draw[g] ~ Dirichlet(n_weighted[period, day_type, g] + 1,  for g = 1..G_gear)
```

In strata with many interviews, the Dirichlet draws are tightly concentrated around the point estimate, and the gear-proportion uncertainty contributes negligibly to the total. In sparse strata — particularly holiday cells and early/late-season periods — the Dirichlet draws vary substantially, widening the credible intervals on gear-type harvest estimates to honestly reflect the limited data.

The framework produces two sets of daily catch quantities to serve different purposes. The expected-value version (`C_gear`) uses the point-estimate `pi_gear` for smooth daily trajectories in plots. The prediction version (`C_gear_pred`) uses the Dirichlet-sampled proportions for season totals and credible intervals. Season-level summaries (`C_sum`, `C_sum_gear`) are always derived from the prediction version so that reported uncertainty bounds account for all quantified sources of variability.

**Assumptions:**
- Gear proportions vary by period and day type. In sparse strata, the model borrows strength from the parent period or sub-season via the fallback hierarchy, and the Dirichlet sampling reflects this reduced information through wider credible intervals.
- Laplace smoothing is adequate for regularization. With sufficient interview data, both the smoothing and the Dirichlet sampling uncertainty have negligible effect.
- The weighted gear assignment from Section 7.2 provides the basis for these proportions.

### 7.6 Observation Models

- Gear counts: `Gear_I ~ Poisson(lambda_E × eps_E_H_obs × R_G)` (sparse overdispersion)
- Trailer counts: `T_I ~ Poisson(lambda_E × eps_E_H_obs × R_T)`
- Interview catch (weighted): `target += w_g × NegBin_lpmf(c | lambda_C_gear[g] × hours, r_C_gear[g])` for each gear type `g` with weight `w_g > 0`

The interview catch likelihood is a weighted pseudo-likelihood. Each interview contributes to multiple gear types' CPUE in proportion to `gear_weights[interview, gear_type]`. For single-gear interviews, this is equivalent to a standard hard-assignment approach. For multi-gear interviews, the full observed catch is evaluated against each reported gear type's CPUE, scaled by the fractional weight. This allows multi-gear interviews to inform all of their reported gear types rather than being arbitrarily assigned to a single type.

### 7.7 Per-Gear-Type Overdispersion

Each gear type has its own negative binomial overdispersion parameter `r_C_gear[g]`. This allows pot catches to have a different variance structure than ring net catches, which is biologically plausible — pots tend to produce more consistent catches while ring net and snare catches are more variable.

**Assumption:** Each gear type has its own catch variance structure.

### 7.8 Daily Catch by Gear Type

The model produces two complementary representations of daily gear-type catch, each serving a different analytical purpose.

**Expected daily catch** (`C_gear`) is the deterministic product of the latent effort, day length, gear-type proportion, and gear-type CPUE for each MCMC iteration:

```
C_gear[d, g] = lambda_E[d] × L[d] × pi_gear[period[d], day_type[d], g] × lambda_C_gear[d, g]
```

This quantity uses the point-estimate `pi_gear` and involves no additional stochastic sampling. Because the posterior already captures uncertainty in effort (`lambda_E`) and catch rate (`lambda_C_gear`), the distribution of `C_gear` across MCMC iterations produces credible intervals that reflect process-level uncertainty. These smooth trajectories are used for daily time-series plots and for comparing BSS to PE on a day-by-day basis.

**Predicted daily catch** (`C_gear_pred`) adds two additional sources of variability on top of the expected value:

```
C_gear_pred[d, g] = Poisson_draw(lambda_E[d] × L[d] × pi_draw[g] × lambda_C_gear[d, g])
```

where `pi_draw` is a Dirichlet sample from the gear-proportion posterior for that day's period × day_type stratum (see Section 7.5). This quantity represents a plausible realized daily catch — the kind of count that would actually be observed if we could enumerate every crab — and reflects both the uncertainty in how many crabbers chose each gear type and the inherent count variability (Poisson noise) around the expected rate. Season totals (`C_sum`, `C_sum_gear`) and their credible intervals are computed from `C_gear_pred`, ensuring that reported uncertainty bounds account for all quantified sources of variability: effort, CPUE, gear-type allocation, and count noise.

**Why both are needed:** Using only the expected value for season totals would understate uncertainty by omitting gear-proportion sampling variability, which can be substantial in sparse strata. Using only the predicted value for daily plots would produce noisy trajectories that obscure the underlying trend. The dual representation provides honest uncertainty for management summaries while preserving clear visual diagnostics.

### 7.9 Sparse Overdispersion

Within-day effort overdispersion parameters are allocated only for actual observations (not every possible day × count-sequence slot), eliminating 64–77% of effort parameters and dramatically improving convergence.

### 7.10 Key Parameters

| Parameter | Description | Prior |
|---|---|---|
| B1 | Weekend effort multiplier (log scale) | Normal(0, 1) |
| B2 | Additional holiday effort multiplier | Normal(0, 1) |
| R_G | Gear per crabber (shore: estimated, boat: fixed at 1.3) | Lognormal(log(1.3), 0.3) or fixed |
| R_T | Trailers per boat group | Beta(0.5, 0.5), guarded |
| phi_E | Effort AR(1) autocorrelation | Beta(2,2) rescaled |
| phi_C_gear | CPUE AR(1) autocorrelation (shared across gear types) | Beta(2,2) rescaled |
| mu_C_gear[G_gear] | CPUE intercept per gear type (log scale) | Normal(log(0.5), 2) |
| sigma_eps_C_gear[G_gear] | Per-gear-type CPUE process error SD | Half-Cauchy(2) |
| r_C_gear[G_gear] | NegBin overdispersion per gear type | Half-Cauchy(2) |
| pi_gear[P_n, 3, G_gear] | Gear-type proportions per period × day_type | **DATA** (empirical, Laplace smoothed) |
| n_weighted_gear[P_n, 3, G_gear] | Raw fractional interview counts per stratum | **DATA** (Dirichlet concentration parameters) |

**Prior scale rationale:** All half-Cauchy priors use scale=2. For log-scale processes, sigma=2 implies the quantity could change by a factor of `exp(2) ≈ 7.4` between periods, which accommodates genuine seasonal variation while preventing implausible magnitudes in data-sparse strata.

**Stan model file:** `crab_bss_gear_resolved.stan` (v3.2)

---

## 8. Convergence and Model Selection

Convergence is assessed on all key parameters, not just the aggregates. The following criteria must all be satisfied for BSS to be selected over PE for a given fit:

**Parameter diagnostics:**

- **Aggregates:** R-hat < 1.05 AND n_eff > 400 for `C_sum` and `E_sum`
- **Per-gear-type catch:** R-hat < 1.05 AND n_eff > 400 for all `C_sum_gear[g]`
- **Per-gear-type overdispersion:** R-hat < 1.05 for all `r_C_gear[g]`
- **Per-gear-type process error:** R-hat < 1.05 for all `sigma_eps_C_gear[g]`

**Sampler diagnostics:**

- **Divergent transitions:** The number of divergent transitions must not exceed the `max_divergences` threshold (default: 0, meaning any divergence causes the fit to fail). Divergent transitions indicate regions of the posterior where the Hamiltonian Monte Carlo sampler cannot maintain accurate trajectories, potentially biasing parameter estimates. Even when R-hat and n_eff appear satisfactory, divergences signal that the reported posteriors may not fully explore the target distribution.
- **Maximum treedepth:** While not a hard pass/fail criterion, the framework issues a warning when more than 5% of iterations hit the maximum treedepth limit. Frequent treedepth saturation suggests that the sampler is working harder than expected to explore the posterior, which can reduce effective sample size and slow runtime. Increasing `max_treedepth` (default: 10) may resolve this at the cost of longer runtime.

All criteria must pass for BSS to be selected. If the aggregates pass but gear-type parameters do not, a warning is issued — this indicates that the total harvest estimate may be reliable but the gear-type decomposition should be interpreted with caution.

The `convergence_report.csv` output includes `divergences`, `treedepth_pct`, `max_gear_rhat`, `min_gear_neff`, `max_r_C_rhat`, and `max_sigma_C_rhat` columns for transparent monitoring.

---

## 9. Output Files

All outputs from the Pooled model, plus gear-type-specific tables:

| File | Contents |
|---|---|
| `catch_by_gear_type_detail.csv` | Gear-type catch with BSS posterior uncertainty (median + 95% CI) |
| `catch_by_gear_type.csv` | Port-level catch by gear type |
| `convergence_report.csv` | Per-fit diagnostics including gear-type parameters and divergence counts |
| `daily_combined_estimate.csv` | Daily PE, BSS, and combined estimates |
| `data_coverage.csv` | Month × population sampling coverage |
| `effort_multipliers.csv` | B1/B2 as human-readable multipliers |
| `expansion_ratios.csv` | R_G and R_T posteriors (with estimation status) |
| `interview_cpue_summary.csv` | Monthly CPUE from raw interviews |
| `pe_vs_bss_comparison.csv` | Component-level PE vs BSS comparison |
| `season_summary.csv` | One-table season summary |
| `sensitivity_incomplete_trips.csv` | PE harvest estimates with and without incomplete trip filter, by component |
| `sensitivity_incomplete_by_gear.csv` | Per-gear-type CPUE comparison between complete and incomplete trips |

Total: 19 CSV files, 10+ plots, and run parameters.

---

## 10. Limitations and Future Directions

### Current Limitations

- No jetty effort counts. Beach crabbing unmeasured. Single-count days in ring-net sub-season.
- Multi-gear interview assignment uses equal weighting across reported gear types; true catch-per-gear-type data would improve precision. The weighted pseudo-likelihood used for multi-gear interviews evaluates the full observed catch against each gear type's CPUE independently (scaled by fractional weight). For gear types where a large fraction of contributing interviews are multi-gear (~50% for Ring Net and Snare), this may narrow credible intervals slightly beyond what the independent information warrants. A latent catch allocation model would address this but adds complexity.
- Shared `phi_C_gear` assumes similar temporal autocorrelation across gear types (per-gear-type process error SD partially addresses this).
- Commercial/charter harvest estimation lacks formal uncertainty quantification.
- The incomplete trip filter reduces CPUE sample size by approximately 35% for shore crabbers. The sensitivity analysis output quantifies the resulting harvest difference. In populations where the incomplete fraction is low (private boat: ~9%), the filter has minimal impact.

### Addressed in v5.2

1. ~~Incomplete trips included without adjustment~~ → Completed-trip filter with sensitivity diagnostics.
2. ~~pi_gear does not contribute posterior uncertainty~~ → Dirichlet sampling in generated quantities propagates gear-proportion uncertainty into season totals.
3. ~~Poisson noise in generated quantities inflates season-total variance~~ → Expected-value daily catch for trajectory plots; Dirichlet + Poisson prediction for season totals.
4. ~~Convergence checks ignore divergent transitions~~ → Divergence count included in pass/fail criterion.
5. ~~Prohibited gear types can be modeled if mentioned in interviews~~ → Regulatory exclusion per sub-season prevents fitting CPUE processes for illegal gear types. Minimum effective-N threshold raised from 3 to 15 as a safety net against phantom gear types from multi-gear interview noise.

### Addressed in v5.1

1. ~~Multi-gear interview assignment is deterministic~~ → Weighted fractional assignment.
2. ~~Single shared process error SD for all gear types~~ → Per-gear-type `sigma_eps_C_gear`.
3. ~~Gear proportions do not vary by day type~~ → pi_gear varies by period × day_type.
4. ~~R_G floats at prior for boat fits~~ → Fixed at 1.3 for boat fits.
5. ~~Convergence checks only on aggregates~~ → Extended to all gear-type parameters.
6. ~~Cauchy priors excessively diffuse~~ → Tightened from scale=5 to scale=2.
7. ~~Day length cap underestimates summer~~ → Raised from 16h to 17h.
8. ~~Gear regex fragile~~ → Word-boundary patterns.

### Planned Improvements

1. Add Westport Jetty direct effort counts.
2. Expand to other Washington coastal ports.
3. Independent `phi_C_gear` per gear type if data supports it.
4. Collect per-gear-type catch data in interviews to replace equal weighting.
5. Full effort decomposition by gear type (separate effort processes per gear).
6. Bootstrap or parametric uncertainty for commercial/charter census component.
7. Latent catch allocation model for multi-gear interviews to replace the weighted pseudo-likelihood.

---

## 11. Glossary

| Term | Meaning |
|---|---|
| **BSS** | Bayesian State-Space model |
| **PE** | Point Estimator |
| **CPUE** | Catch Per Unit Effort (crab per crabber-hour) |
| **Crabber-hour** | One person crabbing for one hour |
| **R_G** | Gear-per-crabber ratio (estimated for shore, fixed for boat) |
| **R_T** | Trailer-per-boat-group ratio |
| **AR(1)** | First-order autoregressive process |
| **Credible interval** | Bayesian range containing the true value with stated probability |
| **R-hat** | Convergence diagnostic comparing MCMC chains |
| **n_eff** | Effective sample size from MCMC |
| **G_gear** | Number of gear types modeled in a sub-season |
| **pi_gear** | Proportion of crabbers using each gear type per period × day_type (data-derived point estimate) |
| **pi_draw** | Dirichlet-sampled gear-type proportions used in prediction quantities to propagate allocation uncertainty |
| **lambda_C_gear** | Gear-type-specific catch rate |
| **sigma_eps_C_gear** | Per-gear-type CPUE process error standard deviation |
| **r_C_gear** | Per-gear-type overdispersion |
| **gear_weights** | Fractional gear-type assignment per interview (rows sum to 1) |
| **n_weighted_gear** | Raw fractional interview counts per period × day_type × gear type, used as Dirichlet concentration parameters |
| **C_gear** | Expected daily catch by gear type (no Poisson noise); used for trajectory plots |
| **C_gear_pred** | Predicted daily catch by gear type (Dirichlet-sampled pi + Poisson draw); used for season totals |
| **Laplace smoothing** | Adding a constant (alpha=1) to counts before normalizing to proportions |
| **Overdispersion** | Extra variability beyond the base distribution |
| **Divergent transition** | MCMC diagnostic indicating difficult posterior geometry; can bias parameter estimates even when R-hat appears satisfactory |
| **Incomplete trip** | Interview conducted before the crabber has finished fishing; excluded from CPUE estimation due to systematic downward bias in catch rate |
| **gear_exclude** | Per-sub-season list of gear types prohibited by regulation; prevents fitting CPUE processes for illegal gear and redistributes fractional interview weights to remaining types |
| **bss_min_gear_effective_n** | Minimum sum of fractional gear weights required for a gear type to be modeled independently (default: 15); prevents near-empty CPUE processes from destabilizing the sampler |

---

## 12. Reproducibility

1. Clone the `FWC-estimation-method` repository.
2. Place input files in `input_files/`.
3. Place `crab_bss_gear_resolved.stan` (v3.2) in `stan_models/`.
4. Open `BSS-GH-gear-type-CPUE-model.Rmd`, update parameters, run.

**Requirements:** R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here.  
**Expected runtime:** 4–5 hours on a 4-core machine.

---

## 13. Change Log

### v5.2 (current)
- Incomplete trip filter: excludes interviews with `trip_status == "Incomplete"` from CPUE estimation to remove systematic downward bias (~20% for pots/traps)
- Dirichlet-sampled gear-type proportions in generated quantities propagate pi_gear uncertainty into season totals and credible intervals
- Expected-value daily catch (C_gear) for smooth trajectory plots; Dirichlet + Poisson predicted catch (C_gear_pred) for season totals
- Divergent transitions included in convergence pass/fail criterion (max_divergences parameter, default 0)
- Treedepth exceedance warnings at >5% of iterations
- Regulatory gear exclusion: per-sub-season `gear_exclude` prevents fitting CPUE processes for prohibited gear types (e.g., Pot during ring-net sub-season); interview weights redistributed to remaining types
- Minimum effective-N threshold raised from 3 to 15 (`bss_min_gear_effective_n`) to prevent phantom gear types from multi-gear interview noise
- Sensitivity analysis output comparing harvest estimates with and without incomplete trip filter
- Per-gear-type CPUE comparison between complete and incomplete trips in diagnostics
- Raw weighted interview counts (n_weighted_gear) passed to Stan as Dirichlet concentration parameters

### v5.1
- Weighted gear-type assignment for multi-gear interviews (equal split across reported types)
- Removed categorical likelihood; pi_gear now data-derived with Laplace smoothing
- pi_gear varies by period × day_type (weekday/weekend/holiday)
- Per-gear-type process error SD (sigma_eps_C_gear[G_gear])
- R_G fixed at 1.3 for boat fits (not estimated)
- Half-Cauchy prior scales tightened from 5 to 2
- Biweekly BSS periods for ring-net sub-season
- Full interview dataset used (no subsampling cap)
- Extended convergence checks (C_sum_gear, r_C_gear, sigma_eps_C_gear)
- Word-boundary regex for gear-type classification
- Day length cap raised from 16h to 17h
- Backward-compatible C[S][D,G] replaced with C_total[D]

### v5.0
- Initial gear-resolved CPUE model
- Per-gear-type AR(1) CPUE processes with shared phi and sigma
- Dirichlet-distributed pi_gear per period
- Priority hierarchy for multi-gear interviews
- Holiday effect B2 in effort process
- Sparse overdispersion parameterization
