# Grays Harbor Recreational Dungeness Crab Harvest Estimation

## Method Version 1.0: Gear-Resolved CPUE Model

**Author:** Matthew George, Ph.D.
**Contact:** matthew.george@dfw.wa.gov
**Agency:** Washington Department of Fish and Wildlife (WDFW)
**Status:** Operational. The method of record for gear-type harvest decomposition at Westport / Grays Harbor, and the companion to the published pooled-CPUE model (`BSS-GH-pooled-CPUE-model-documentation.md`). Use the pooled model for the headline total; use this model when gear-type catch structure is the management question.
**Method version:** 1.0 (corresponds to the gear-resolved pipeline framework code v5.5; reference calibration season 2024-25).
**Code version:** framework v5.5, 2026-07-11.

------------------------------------------------------------------------

### How to read this document

This is the single authoritative reference for the gear-resolved harvest estimation pipeline (`01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd`), the Stan model it fits (`02_stan_models/crab_bss_gear_resolved.stan`), the diagnostics it runs, and the inputs and outputs it uses. It is written for two audiences and is split into three parts:

- **Part I (Sections 1-6): For everyone.** Plain-language description of what the method does, the fishery, the data, how the estimate is built, and where it is valid. No statistics background required.
- **Part II (Sections 7-13): Running it next season.** The operational guide: prerequisites, the single-file run, how to judge whether a season's estimate is trustworthy, the output catalog, the diagnostics, reproducibility, and the conditions under which the method stops applying.
- **Part III (Sections 14-20): Technical reference.** The full model specification, design rationale, limitations, glossary, and references.

The development history (how the model reached v5.5, the full version-by-version change log, the run-driven Stan fix markers, and the boat effort-scale and point-estimator working notes) has been moved out of this document to keep it focused on the published method. It lives in `BSS-GH-gear-type-CPUE-model-development-history.md`. Section 19 gives a one-screen summary and points there.

A note on naming: this is the "gear-resolved" model because it is written to estimate a catch-rate (CPUE) process per gear type. The companion "pooled" model, documented in `BSS-GH-pooled-CPUE-model-documentation.md`, uses a single CPUE process shared across gear types and is the published v1, because it is the more robust of the two and answers the primary management question (total harvest with defensible uncertainty). The gear-resolved model is the alternative when the gear-type decomposition is what is needed. An important current caveat, stated up front and in full in Section 16: the model is written for a per-gear CPUE process, but the production configuration runs with a single gear dimension (G = 1), so the per-gear machinery is inert and gear-type catch is apportioned from interview proportions rather than independently modeled. The fully modeled per-gear path is specified but not yet built.

================================================================

# PART I: FOR EVERYONE

------------------------------------------------------------------------

## 1. What this method produces

This framework estimates the total recreational Dungeness crab (*Metacarcinus magister*) harvest at Westport and the greater Grays Harbor area for a season, and organizes that harvest by gear type (pots, ring nets, foldable/star traps, snares). It combines four kinds of field observations (gear counts at the docks, trailer counts at the boat launch, dockside crabber interviews, and ingress/egress surveys) with a statistical model that fills in the days when no sampling occurred.

Each run produces:

- A total Dungeness crab harvest estimate for the port, with a 95% credible interval (a range that has a 95% probability of containing the true harvest, given the data and model). This total is structurally identical to the pooled model's total.
- A gear-type breakdown of harvest with its own uncertainty, a season-long CPUE trajectory per gear type, and time-varying gear-type proportions showing how the mix of gear in use shifts across the season and between day types.
- Breakdowns by crabbing mode (shore, private boat, commercial/charter) and monthly harvest trends.
- A weekend/holiday catch-rate effect (whether weekend and holiday catch rates differ from weekday rates).
- A sensitivity analysis quantifying how the incomplete-trip filter affects harvest by population, sub-season, and gear type.

**Why this matters for management.** Gear regulations are a primary management tool for the crab fishery. If managers are weighing a pot restriction, an extension of the ring-net-only period, or an evaluation of snare effectiveness, they need to know how much harvest each gear type contributes and how confident the estimate is.

**A caveat on the current per-gear numbers.** In the production configuration the gear dimension is set to one, so the headline gear-type catch is apportioned from interview proportions rather than carried through an independent per-gear CPUE posterior. The credible intervals on the total are full; the per-gear split is an apportionment, not a modeled quantity with its own posterior. The architecture for the modeled per-gear path exists in the Stan model (Part III) and is tracked as the next build (Section 16, GR-7). Read the per-gear split as a proportional allocation until that path is activated.

**How confident are we?** The framework runs two independent estimation methods, a simple average-based approach (the Point Estimator, PE) and a Bayesian time-series model (the Bayesian State-Space model, BSS), then compares them. When the two agree and the BSS passes its convergence checks, confidence is high. The output includes formal diagnostics and a side-by-side comparison so a reviewer can judge reliability. Section 9 explains how to read those checks.

------------------------------------------------------------------------

## 2. The fishery and study area

The recreational Dungeness crab fishery at Westport is one of the highest-volume recreational crabbing operations on the Washington coast. Thousands of crabbers participate annually using four primary gear types: crab pots (highest catch rate), ring nets, foldable/star traps, and snares. The fishery operates year-round, though effort peaks during summer months and around major holidays. WDFW rules prohibit pots from late September through November, which creates a structural break in both effort and catch rates and is the reason the season is split into two sub-seasons (Section 5).

Commercial Dungeness crab vessels also crab recreationally before the commercial season opens, under the same daily limits as private boats, but tend to have much higher catch rates per vessel. Their harvest is tracked separately through a vessel tally at the marina.

Westport sits on the south side of the Grays Harbor estuary. Recreational crabbing occurs from public docks (Floats 17-21), a jetty, beaches, a public boat launch, and the commercial marina. The "shore" component pools dock, jetty, and beach crabbing. Each access point serves a different crabbing mode and requires a different type of effort measurement.

**Why estimate harvest by gear type.** Beyond the total, gear-specific information supports several management needs: evaluating whether a proposed gear restriction would meaningfully reduce harvest, monitoring whether a gear type is becoming more or less efficient over time (which can indicate changes in crab abundance or behavior), and understanding how the recreational fleet adapts when gear regulations change.

------------------------------------------------------------------------

## 3. The data streams and field collection

The framework requires four data streams. Three are collected by WDFW field staff using Apple iPads running the WDFW iForm "Crab Creel Survey" application; the fourth is the all-day ingress/egress survey.

| Stream | What is collected | What it tells the model |
|---|---|---|
| **Effort counts** | Instantaneous point-in-time counts of crab gear at the docks and boat trailers at the launch | The primary indicator of how much crabbing activity is happening |
| **Crabber interviews** | Dockside trip-level records: group size, gear deployed and type, hours fished, crab kept, trip status | Catch rate (CPUE), the gear mix in use, and the gear-per-group ratios |
| **Commercial/charter tally** | Daily count of commercial and charter vessels at the marina during the recreational pre-season | The commercial/charter component of harvest, via expansion |
| **Ingress/egress (I/E) surveys** | All-day surveys recording crabber arrivals and departures | A direct measurement that identifies the deployment turnover (tau) used to expand effort |

### 3.1 Effort counts

Field staff visit sites and record the number of crabbing indicators visible. Protocol calls for multiple counts per day at standardized times; the number of within-day counts affects the model's ability to estimate within-day variability.

| Site | What is counted | Role in model |
|---|---|---|
| Westport Docks Float 20 | Crab gear in the water | Primary shore effort indicator |
| Westport Docks Float 17-21 | Crab gear in the water | Summed with Float 20 for the section total |
| Westport Boat Launch | Boat trailers at ramp | Private boat effort indicator |
| Ocean Shores Boat Launch | Boat trailers at ramp | Supplementary boat effort |
| Westport Jetty | Crabbers (future) | Reserved for future use |

Float 20 and Float 17-21 counts are paired by time and summed, because dock crabbers distribute across both float areas.

**Assumptions.** A point-in-time count reflects relative effort on that day; time-of-day introduces noise but not systematic bias; the surveyor's count is accurate.

### 3.2 Crabber interviews

Interviews record trip-level information: number of crabbers, number of gear units, gear type(s) (select all that apply from pot, ring net, trap, snare), hours fished, crabber-hours, Dungeness crab kept, Red Rock crab kept, trip completion status, crabbing mode (Dock, Boat, Jetty, Beach), and boat type (Private, Commercial, Charter, Guide). They provide the catch rate, the gear-per-group ratios (R_G for shore, R_G_boat for boats), and, critically for this model, the gear type(s) each crabber used.

**Gear-type recording.** The iForm allows multiple gear-type selections per interview. About 28% of valid Grays Harbor interviews report multiple gear types (for example "Pot, Ring Net"). The model handles this through weighted fractional assignment across the reported gear types (Section 14.4).

**Population classification.**

| Crabbing mode | Boat type | Population |
|---|---|---|
| Dock, Jetty, or Beach | Any | **Shore** |
| Boat | Private or blank | **Private Boat** |
| Boat | Commercial, Charter, or Guide | **Commercial/Charter** |

**Trip-completion filtering.** Field staff record whether the crabber has finished the trip at interview. Only completed-trip interviews are used for CPUE estimation, controlled by `filter_incomplete_trips` (default TRUE). Crabbers interviewed mid-trip have systematically lower catch than they would at trip end, because their gear has not finished soaking. In the 2024-25 Grays Harbor data, roughly 35% of shore interviews are incomplete trips, and these show a mean CPUE about 20% lower than completed trips. The bias is not uniform across gear: it is largest for soak-time-dependent gear (pots about -21%, traps about -23%) and negligible for ring nets (about +4%), which are checked and pulled frequently. Including incomplete trips unadjusted would underestimate shore harvest by roughly 7% and, more consequentially for this model, would differentially suppress pot and trap catch rates relative to ring net, distorting the gear-type decomposition. Interviews with missing trip status are kept (a blank status may be a completed trip with incomplete metadata). The effort counts are not affected by this filter; it applies only to the interviews used for catch-rate estimation. The run writes `sensitivity_incomplete_trips.csv` (harvest with the filter off vs on) each season.

**Assumptions.** Interviewed crabbers are representative of their population, conditional on trip completion; crabbers accurately report catch, hours, group size, and gear types; the gear type(s) reported reflect what was actually deployed; completed trips give an unbiased estimate of the day's average catch per unit effort, and incomplete trips are a left-censored observation excluded to avoid a downward bias.

### 3.3 Commercial/charter vessel tally

Daily counts of commercial crab vessels and charter boats at Westport Marina during the pre-season window when these vessels crab recreationally. Commercial vessels have much higher per-vessel catch and are hard to interview comprehensively, so the tally provides a census-like measure of vessel-days for a day-type-stratified expansion.

**Assumptions.** The tally captures all participating vessels on each sampled day; mean harvest per vessel from interviews applies to uninterviewed vessels; the tally is stratified by day type to prevent sampling bias.

### 3.4 Ingress/egress (I/E) surveys

All-day surveys that record crabber arrivals and departures. As of the shared day-length module the gear-resolved pipeline uses the same I/E workbook the pooled model uses, to identify the deployment turnover `tau` that expands effort (Section 14.2). When enough I/E days fall inside a fit's window, they constrain `tau`; when they do not, `tau` rests on its prior (a standing dependence for the boat; Section 16, GR-12).

The four input files that carry these streams are listed in Section 7; their exact schema and known quirks are documented in `04_input_files/README.md`.

------------------------------------------------------------------------

## 4. How the estimate is built

Field crews cannot sample every day. Both methods solve the same problem (estimate harvest on unsampled days), but differently.

**The Point Estimator (PE): a simple average.** For each stat-week by day-type group, the PE averages the daily harvest on sampled days and multiplies by the number of days in that group (Pollock et al. 1994; Hahn et al. 2000). Stratum catch rates are formed as ratio-of-sums (total catch over total effort), the harvest-consistent estimator. The PE is transparent and assumption-light, but it cannot fill a group that had zero samples and it produces no uncertainty bounds.

**The Bayesian State-Space model (BSS): a time-series curve.** The BSS fits a smooth curve through the daily effort and catch-rate data using a statistical time-series process, then uses that curve to estimate every day in the season, including unsampled days (Conn 2002; Staton et al. 2017). In this model, the catch side is written so that each gear type can carry its own catch-rate curve; the effort side is identical to the pooled model. The BSS accounts for correlation between adjacent days, fills gaps with honest uncertainty that grows the further a day is from the nearest observation, and produces credible intervals. It is more complex, takes roughly four to five hours per full run, and must be checked for convergence.

**Combining them.** For each population component, the framework checks the BSS fit against formal convergence criteria (Section 9). If the fit passes, its estimate is used; if not, the PE estimate is used as a fallback. The two are reported side by side so a reviewer can see where they agree and where they differ. The headline number uses the BSS posterior expected catch (the model's best estimate of the average catch) rather than a single simulated draw, following the standard distinction between estimation and prediction in hierarchical models (Gelman et al. 2013, Ch. 7).

------------------------------------------------------------------------

## 5. The three population components and the two sub-seasons

The harvest is built from three components, estimated separately and summed:

1. **Shore crabbers** (dock + jetty + beach). Effort is indicated by gear counts at the docks.
2. **Private boat crabbers.** Effort is indicated by trailer counts at the boat launch.
3. **Commercial/charter vessels** crabbing recreationally pre-season. Estimated by a day-type-stratified expansion of the marina vessel tally (about 12,007 crab in the 2024-25 reference runs; this component is a census and does not carry a BSS fit).

The season is split into two **sub-seasons**, defined by the pot closure, and each is estimated independently:

| Sub-season | Typical dates | Gear allowed | Gear types modeled | BSS period type |
|---|---|---|---|---|
| Ring-net only | Sep 16 to Nov 30 | Ring nets, snares, foldable traps | 3 (Ring Net, Trap, Snare) | Biweekly |
| All-gear | Dec 1 to Sep 15 | All gear including crab pots | 4 (Pot, Ring Net, Trap, Snare) | Monthly |

The number of gear types set up for modeling in a sub-season is determined by two filters applied in sequence.

**Regulatory exclusion.** Each sub-season carries a `gear_exclude` list of gear prohibited by regulation during that period. For the ring-net sub-season, pots are excluded. This is enforced structurally: regardless of interview contents, the model will not set up a CPUE process for a prohibited gear type. Any interview that mentions an excluded gear type (typically a multi-gear crabber naming gear they own rather than what they deployed) has its fractional weight redistributed to the remaining reported types. An interview recorded as "Pot, Ring Net" during the ring-net sub-season contributes 100% of its weight to Ring Net rather than 50% to each. This matters because even a handful of phantom gear mentions (5 to 6 fractional interviews in a 76-day sub-season) can, with near-zero data to constrain a CPUE process, drive hundreds of divergent transitions and overflow in the generated catch. Encoding the known fishery structure prevents that.

**Minimum effective-N threshold.** After regulatory exclusions, each remaining gear type must have at least `bss_min_gear_effective_n` effective interviews (default 15), computed as the sum of fractional gear weights, to be set up independently. The threshold is conservative because fractional weights from multi-gear interviews can inflate the effective count. If fewer than two gear types qualify, the model falls back to a single "All" category.

The ring-net sub-season uses biweekly BSS periods rather than monthly: with only about 76 days and 2.5 calendar months, monthly periods give too few AR(1) transitions for effective temporal smoothing, while biweekly periods yield about 5 to 6 periods.

**Day type and the holiday effect.** Each day is classified as weekday (Mon-Thu), weekend (Fri-Sun), or holiday. The effort process includes two separate effects: a weekend boost (B1) and an additional holiday boost beyond the weekend (B2), because holidays generate substantially higher effort than regular weekends. Gear-type proportions vary by day type as well as by period. The single season list of crabbing holidays lives in `run_config.R`.

**Assumption.** The pot-open date creates a clean structural break, and the day-type classification is the same across all populations.

------------------------------------------------------------------------

## 6. Where this method is valid

Method v1.0 is calibrated to the Westport / Grays Harbor fishery as sampled in the 2024-25 season. It is designed to be re-run in future seasons **provided the fishing location, the input data streams, and the sampling design remain the same.** Section 13 sets out which assumptions are baked in and the specific conditions under which the method must be re-derived rather than re-run. In short: a different port, a change in how effort counts are taken, or a structural change in who participates would each require revisiting the method, not just feeding it new data.

================================================================

# PART II: RUNNING IT NEXT SEASON

------------------------------------------------------------------------

## 7. Prerequisites and repository layout

**Software.** R 4.2 or later, with rstan 2.32 or later and a working C++ toolchain (rstan compiles the model), plus the packages tidyverse, lubridate, suncalc, gt, patchwork, here, and readxl. The gear-resolved driver now reads the I/E workbook through the shared day-length module, so `readxl` is required even though earlier gear-resolved runs did not need it.

**Repository.** The code lives in the `Coastal-Rec-Crab-BSS` repository, organized into numbered stage folders. All reads and writes go through `here::here()`, which anchors paths to the repository root, so you do not edit paths to run it; you place files in the right folders.

| Folder | What it holds | Your job |
|---|---|---|
| `01_BSS_models/` | The driver `BSS-GH-gear-type-CPUE-model.Rmd` | This is the file you knit |
| `02_stan_models/` | `crab_bss_gear_resolved.stan` | Leave in place; the driver compiles it |
| `03_R_functions/` | Helper functions (auto-sourced wholesale) | Leave in place |
| `04_input_files/` | The four input files | Replace with the new season's data, same names and schema |
| `05_output/` | Dated run folders | The run writes here; nothing to place |

**The four input files** (see `04_input_files/README.md` for schema and quirks):

- `effort_combined.csv` (effort counts; re-exported with full quoting)
- `interview_combined.csv` (crabber interviews with the gear-type field; dates M/D/YYYY)
- `wes_commercial_tally.csv` (commercial/charter tally; one row per tally day)
- `ingress_egress.xlsx` (I/E surveys; named through the `ie_data_file` / `ie_sheet` settings)

**One file you edit, and where tuning lives.** As of the 2026-07-11 consolidation, every user-selectable toggle (season window, structural dates, catch groups, day-typing, effort unit, incomplete-trip filter, I/E settings, holidays, and the model-behavior levers) lives in `run_config.R`, the single source of truth. The two model drivers no longer carry their own copies of these keys, so there is nothing to keep in sync. Each driver keeps only its own model-internal tuning in its `params` block (the Stan filename, per-fit sampler settings, the convergence-gate thresholds, the AR-selector thresholds, and a few gear-model constants), which you rarely touch and which legitimately differs from the pooled model. When the driver is knit standalone, its setup chunk sources `run_config.R` automatically (an `if (!exists("run_config")) source(...)` guard), so a standalone knit uses exactly the same toggles as an orchestrated run.

**Code organization (2026-07-11 refactor).** The helper functions that used to be defined inline in the driver were extracted into `03_R_functions/` so both tracks share one implementation and cannot drift. Shared helpers include `prep_days_crab.R` (calendar and day-typing), `prep_population_summary.R` (population by sub-season filtering), `estimate_comm_charter.R` (the census expansion), and `bss_timers.R` (run timing). Gear-specific helpers include `fetch_crab_data_v2.R` (input assembly with weighted gear classification), `run_pe_gear.R` (the gear Point Estimator), and `prep_bss_crab_gear.R` (the Stan data list). The effort-unit specification (`bss_effort_spec.R`), the scale-aware convergence gate, the AR-resolution selector, the CPUE effort-unit diagnostics, and the I/E day-length model are all shared modules used by both drivers.

------------------------------------------------------------------------

## 8. Step-by-step: running a season

1. **Place the new season's data** in `04_input_files/`, keeping the four filenames and column schemas unchanged. Honor the schema quirks in the input-folder README (re-export the effort CSV with full quoting; interview dates are M/D/YYYY; the commercial boat-type spelling is matched by a case-insensitive regex).
2. **Edit `run_config.R`**, the one file you touch. Set `model <- "gear_resolved"`, then the values you change most often: `est_date_start` and `est_date_end` (the season window), `season_filter`, the structural dates (`pot_open_date`, `census_start_date`, `census_end_date`), and the season's `crabbing_holiday_dates`. Leave `bss_seed` fixed (default 20260619) for reproducibility. The effort unit is set here too (`shore_effort_unit`, defaulting to `"gear-deployments"`; Section 14.2). Do not run the weather module with this model: `run_weather` is only valid with the pooled model, and the orchestrator stops early if you set it TRUE here.
3. **Knit the driver** `01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd`, or launch through the orchestrator. Each population by sub-season is fit independently. Expect a total runtime of roughly four to five hours on a 4-core machine (about 50% longer than the pooled model, because there are more parameters to sample and monitor).
4. **Check convergence** for each fit using `convergence_report.csv` and the rules in Section 9. A fit that fails falls back to PE automatically.
5. **Read the outputs** from `05_output/YYYYMMDD/gear-type-CPUE-model/` (Section 10), starting with the port total, `pe_vs_bss_comparison.csv`, and `catch_by_gear_type.csv`.

------------------------------------------------------------------------

## 9. Judging whether a season's estimate is trustworthy

This is the most important section for an operator. For each BSS fit, the framework applies the shared scale-aware convergence gate, reported per fit in `convergence_report.csv`. A fit **passes**, and its BSS estimate is used for that component, when all of the following hold; otherwise the PE estimate is used:

- **R-hat < 1.01** for the seasonal totals `C_expected_sum` (equivalently `C_sum`) and `E_sum`, following Vehtari et al. (2021). R-hat near 1.00 means the independent sampler chains agree.
- **n_eff > 400** for both totals. This is the effective number of independent posterior samples.
- **Divergent fraction below 0.05** (the hard backstop, `max_divergence_fraction`). Above this rate the sampler's geometry is untrustworthy and the fit is rejected regardless of anything else.
- **Divergences do not move the answer.** The shift the divergent draws induce in each total, measured in units of that total's posterior standard deviation, is below 0.10 SD (`max_divergence_total_impact_sd`).

A count-based divergence criterion (`max_divergences`, default 5) is retained for reporting but is superseded for the PE/BSS decision by the scale-aware gate above.

**Per-gear monitoring.** Because this model is written to resolve catch by gear type, convergence is also reported on the per-gear quantities: R-hat and n_eff for each `C_sum_gear[g]`, and R-hat for the per-gear overdispersion and process-error scales. If the aggregates pass but a gear-type quantity does not, the run issues a warning: the total may be reliable while the gear-type decomposition should be read with caution. The report carries `divergences`, `treedepth_pct`, `max_gear_rhat`, `min_gear_neff`, `max_r_C_rhat`, and `max_sigma_C_rhat` for transparent monitoring. In the current G = 1 configuration these per-gear columns describe a single-gear posterior (Section 16, GR-7).

Why divergences are in the gate: a sampler that cannot accurately integrate its trajectory is not faithfully exploring the target distribution, and can bias the posterior even when R-hat and n_eff look fine (Betancourt 2017). Why the impact criterion is measured in standard deviations and not as a percentage of the estimate: a percentage-of-level threshold penalizes a component for having a wide posterior rather than for being biased; the SD-normalized criterion asks whether the divergences move the answer relative to how well the answer is pinned down, and is invariant to posterior width.

**A standing gate caveat (GR-11).** The hard 5% divergence backstop can override the scale-aware impact test. On the 2026-07-09 reference run, `shore_ring_net` passed R-hat, n_eff, and the impact test (divergences moved the totals by only 0.004 SD) but failed the backstop alone (7.5% vs 5%), so a good BSS estimate (about 6,930, matching the pooled-scale expectation of about 6,733) was replaced by an inflated PE estimate. The backstop is kept on principle, but the 5% threshold and its precedence over the impact test are an open policy question. A later run brought that same fit below the backstop, so it now reports on its BSS posterior.

**Reading the comparison.** `pe_vs_bss_comparison.csv` shows PE and BSS effort and catch by component with the selected method. Large PE-vs-BSS gaps are not automatically errors; they can reflect a real disagreement between the design-based expansion and the model's reconciliation against interview data. The gate decision is computed once per fit by the shared gate module and consumed by every downstream summary, so the convergence report and the comparison always agree on which method was used.

------------------------------------------------------------------------

## 10. Output catalog

Each run writes to `05_output/YYYYMMDD/gear-type-CPUE-model/`. Files tagged with a population follow the pattern `<metric>_<population>_Dungeness_Kept.{csv,png}`, where population is one of `shore_ring_net_only`, `shore_all_gear`, or `private_boat_all_gear`. The commercial/charter component has no separate BSS file; it enters the port total by census expansion.

**Headline estimates and gear-type detail**

| File | Contents |
|---|---|
| `port_total_Dungeness_Kept.csv` | Combined PE + BSS port total (expected and predictive), with its component composition |
| `pe_port_summary.csv`, `pe_vs_bss_comparison.csv`, `monthly_pe_vs_bss.csv` | PE estimates and the side-by-side PE vs BSS reconciliation used by the gate |
| `catch_by_gear_type.csv` | Port-level catch by gear type |
| `catch_by_gear_type_detail.csv` | Gear-type catch with BSS posterior uncertainty (median + 95% CI) |
| `monthly_estimates.csv`, `monthly_estimates_by_mode.csv`, `catch_by_mode.csv` | Monthly catch and effort, and catch by crabbing mode, with credible intervals |
| `season_summary.csv` | One-table season roll-up |
| `sensitivity_incomplete_trips.csv` | PE catch with the incomplete-trip filter off vs on, by component |
| `sensitivity_incomplete_by_gear.csv` | Per-gear-type CPUE, complete vs incomplete trips |

**Convergence, structure, and expansion (per fit unless noted)**

| File | Contents |
|---|---|
| `convergence_report.csv` | One file, all populations: R-hat, n_eff, divergent count and fraction, tree-depth, AR resolution, the per-gear R-hat/n_eff columns, and the SD-normalized divergence impact (the gating criterion) |
| `structural_params_<label>.csv` | Posterior summary of scale/structural parameters (sigma_eps, phi, r, sigma_mu, sigma_IE, R_G, R_G_boat, tau) with CI, n_eff, R-hat |
| `expansion_ratios.csv` | R_G and R_G_boat posteriors, with the R_G_boat expansion (replacing the removed R_T) |
| `divergence_localization_<label>.csv`, `sampler_diagnostics_<label>.csv`, `prior_vs_posterior_<label>.csv` | Where divergent draws sit per parameter, HMC sampler diagnostics including E-BFMI, and prior-vs-posterior overlap |
| `effort_cpue_multipliers.csv` | B1, B2, and the B1_C weekend/holiday CPUE effect as human-readable multipliers |

**Posterior predictive, cross-validation, and the CPUE effort-unit checks (per fit)**

| File | Contents |
|---|---|
| `ppc_calibration_<label>.csv`, `ppc_pit_<label>.png`, `ppc_byobs_<label>.csv` | Posterior predictive coverage, PIT histograms, and per-observation residuals for effort counts and interview catches |
| `effort_overdispersion_decomp_<label>.csv`, `effort_overdispersion_byobs_<label>.csv` | Effort-variance decomposition (Section 11) |
| `loo_summary_<label>.csv`, `loo_pointwise_*_<label>.csv` | PSIS-LOO summaries and pointwise contributions by likelihood component (gear/trailer/catch) |
| `cpue_estimators_<label>.csv`, `cpue_saturation_<label>.csv`, `cpue_linearity_<label>.csv` | The CPUE effort-unit checks central to this model: estimator triad, saturation exponent, and effort linearity (Section 11) |

**Daily series, plots, and metadata**

Per-fit daily series (`bss_daily_effort_*`, `bss_daily_cpue_*`, `bss_daily_catch_*`), per-fit summaries and the AR path (`bss_summary_*`, `bss_full_summary_*`, `bss_ar_path_*`, `bss_period_coverage_*`, `bss_draws_summed_*`, `bss_L_effective_*` for shore), the I/E detail (`ie_analysis.csv`, `L_effective_ie_detail.csv`), a family of plots (daily series, posteriors, monthly catch total and by mode, catch by gear type, the L_effective regression and day-length comparison), and `run_parameters.txt` / `session_info.txt` recording the exact parameters and the R/package/Stan session and seed. A complete, categorized listing is in `05_output/README.md`.

------------------------------------------------------------------------

## 11. Diagnostics: what each one answers

The diagnostics are additive (each is wrapped so a failure cannot break a run) and are written every run. The CPUE effort-unit checks are the ones most central to this model, because the harvest total is a constant-CPUE season expansion and is only unbiased if catch scales with the chosen effort unit.

**Posterior predictive checks (PPC).** `ppc_calibration_<label>.csv` and `ppc_pit_<label>.png` ask whether the model's predictions are calibrated against the actual effort counts and interview catches. A well-calibrated model has PIT values spread uniformly; a central hump means the predictive is too wide (over-dispersed), and 50% coverage above the nominal 0.50 says the same.

**Effort over-dispersion decomposition.** `effort_overdispersion_decomp_<label>.csv` splits each effort observation's predictive variance into three additive parts via the law of total variance:

```
Var(Y) = E[mu]            (Poisson floor: irreducible, not a lever)
       + E[mu^2 / r_E]    (NB observation over-dispersion: controlled by the r_E / sigma_r_E prior)
       + Var(mu)          (latent process + parameter uncertainty: controlled by sigma_eps_E)
```

The `lever` column reports the verdict: if the NB-overdispersion share dominates, the lever is the `r_E` / `sigma_r_E` prior; if the latent share dominates, the lever is the AR innovation scale. In this model the coarse fixed AR is by design (`ar_adaptive = FALSE`, biweekly/monthly periods), so the latent AR carries only a small share of the effort variance and the NB observation term carries most of it. In the reference run the latent share is about 4 to 13% and the NB-overdispersion share about 81 to 91% across the three fits (GR-14). This is expected under a coarse AR, not a defect; day-level variance simply lands in `r_E` rather than in the latent path. It is the reason the AR contributes little to the effort fit here.

**PSIS-LOO.** `loo_summary_<label>.csv` reports out-of-sample predictive performance (expected log predictive density, `elpd_loo`) and the Pareto-k influence diagnostic per likelihood component (gear/trailer/catch). This is the basis for principled model comparison, and it is the machinery the shore effort-unit decision was read off (Section 15).

**CPUE effort-unit checks.** `cpue_linearity_<label>.csv`, `cpue_saturation_<label>.csv`, and `cpue_estimators_<label>.csv` test the likelihood's core assumption that catch is proportional to the chosen effort denominator `h`. The linearity check fits `glm(catch ~ log(h))` to estimate a scaling exponent `beta_h`; the likelihood assumes `beta_h = 1`, so a value well below 1 means the effort unit over-counts effort at the high end and the season expansion is biased. The saturation check bins catch per gear by soak time and reports the per-unit rate gradient. The estimator triad compares the model-implied CPUE against the ratio-of-sums and the mean-of-ratios; drift toward the mean-of-ratios warns that the negative-binomial dispersion is pulling the fitted rate off the rate scale. The run also asserts that effort `E` and the CPUE denominator `h` carry the same unit before sampling. These are the checks that surfaced the effort-unit corrections for both components.

**Why the deployment is the effort unit (saturation).** Binned by soak time, crab per gear-HOUR falls about 43-fold across the range of soak durations, while crab per gear per trip rises only about 1.8-fold; a log-log fit over the boat interviews gives catch per gear scaling as soak-hours to the power about 0.13. In plain terms, soak time barely matters, so a pot lift is a pot lift whether it soaked two hours or eight. That makes the deployment the unit on which catch-per-unit-effort is a stable rate: roughly 4 to 7 crab per pot lift, steady across soak times. A stable rate is exactly what the harvest method needs, because harvest is effort multiplied by that rate, and the multiplication is only unbiased if the rate does not drift with the effort denominator. On the 2026-07-10 shore comparison, gear-deployments is the only shore unit whose `beta_h` covers 1 (1.05, 95% CI 0.94 to 1.15), while crabber-hours (0.57) and gear-hours (0.73) both fall well short.

**Incomplete-trip sensitivity.** `sensitivity_incomplete_trips.csv` and `sensitivity_incomplete_by_gear.csv` report the harvest and per-gear CPUE impact of `filter_incomplete_trips` (default on) so the data-quality decision is explicit each run.

------------------------------------------------------------------------

## 12. Reproducibility

The Stan fits take a fixed RNG seed (`bss_seed`, default 20260619), passed to `rstan::stan()`. rstan seeds each chain from `bss_seed + chain_id`, so the chains still differ (R-hat remains meaningful) while run-to-run variation is removed. Package and Stan versions and the seed are written to `session_info.txt` with each output set. Change `bss_seed` only if a pathological seed is ever suspected.

To reproduce a run: clone `Coastal-Rec-Crab-BSS`; place the four input files in `04_input_files/` (their names and schema unchanged); leave `crab_bss_gear_resolved.stan` in `02_stan_models/` and the helpers in `03_R_functions/`; set the season toggles in `run_config.R`; and knit `01_BSS_models/BSS-GH-gear-type-CPUE-model.Rmd` (or launch through the orchestrator). You do not edit paths (`here::here()` resolves them) and you do not edit the driver (it sources `run_config.R` and holds only model-internal tuning). Requirements: R 4.2+, rstan 2.32+, tidyverse, lubridate, suncalc, gt, patchwork, here, readxl. Expected runtime is about four to five hours on a 4-core machine.

------------------------------------------------------------------------

## 13. Scope: when this method applies, and when it must be re-derived

Method v1.0 is built for one fishery under one sampling design. It can be re-run season after season as long as the following hold. Where one breaks, the method must be revisited, not merely re-fed.

**Assumptions that allow a straight re-run:**

- **Same location.** Westport / Grays Harbor access points (docks Floats 17-21, the jetty, beaches, the boat launch, the marina). The gear-per-crabber prior `R_G`, the gear-per-boat-group prior `R_G_boat`, the deployment turnovers `tau_shore` and `tau_boat`, and the effective-day-length regression are all calibrated to this site.
- **Same input streams, same schema.** The four input files in the same form (Section 7).
- **Same sampling design.** Instantaneous effort counts, dockside interviews, the commercial tally, and I/E surveys, collected as in 2024-25. The 2024-25 protocol of multiple randomized effort counts per day is the design the effort expansion assumes (it measures mean daily effort).
- **Same sub-season structure.** Ring-net only Sep 16 to Nov 30, all-gear Dec 1 to Sep 15, tied to the pot closure, with the same `gear_exclude` regulation.
- **Same gear regulations.** The set of legal gear types, and which are prohibited when, is what `gear_exclude` and the sub-season split encode.

**Conditions that require re-derivation, not just new data:**

- **A different port.** `R_G`, `R_G_boat`, the `tau` turnovers, and the L_effective regression would have to be re-estimated from that port's I/E and interview data; the access-point structure differs.
- **A change in the effort-count protocol.** Reverting to a single peak-time count per day measures a different quantity (peak, not mean daily effort) and would bias the effort level. Mixing protocols across years is a genuine confound, addressable only with a protocol fixed effect and a peak-to-mean calibration.
- **A change in gear regulations.** A new legal gear type, or a change in which gear is prohibited when, changes the sub-season structure and the `gear_exclude` logic, and would require re-checking the minimum effective-N and the gear classification.
- **A structural change in participation.** For example, a change in how commercial/charter vessels participate pre-season, or the opening of a new major access point (a jetty effort count, currently absent), would change what the components represent.
- **A move to genuinely modeled per-gear catch.** Activating the per-gear CPUE path (raising the gear dimension above one; Section 16, GR-7) is a re-derivation, not a re-run: it requires an effort-share offset and a rule for multi-gear interviews, and it changes what the gear-type numbers mean.

================================================================

# PART III: TECHNICAL REFERENCE

------------------------------------------------------------------------

## 14. Model specification (`crab_bss_gear_resolved.stan`)

### 14.1 Effort process (shared with the pooled model)

```
log(lambda_E[d]) = mu_E + omega_E[period(d)] + B1 * w[d] + B2 * holiday[d]
```

The temporal deviation `omega_E` evolves as an AR(1) process:

```
omega_E[p] = phi_E * omega_E[p-1] + sigma_eps_E * epsilon[p-1]
```

where `period(d)` maps day `d` to its AR period index (biweekly for the ring-net sub-season, monthly for all-gear). Innovations are standard normal under a non-centered parameterization for efficient HMC (Papaspiliopoulos et al. 2007). The AR(1) initial state is non-centered: the initial deviation is a raw standard normal scaled by the stationary standard deviation `sigma_eps_E / sqrt(1 - phi_E^2)`, so the process starts from its stationary distribution without a centered funnel (Harvey 1989; Betancourt and Girolami 2015; parity port B1.3). `B1` is the weekend effort boost and `B2` is the additional holiday boost beyond a weekend.

**Adaptive resolution is off by default.** The AR resolution is fixed per sub-season (biweekly ring-net, monthly all-gear) with `ar_adaptive = FALSE`. The data-driven selector and per-population caps (shared `bss_ar_resolution.R`) are available by setting `ar_adaptive = TRUE`, but that is inference-changing and should be validated first. The fixed coarse AR is why the latent process carries only a small share of effort variance (Section 11, GR-14).

### 14.2 Population effort construction and the gear-deployment effort unit

As of framework v5.5 **both fitted components run on gear-deployments** (pot lifts), routed through the shared module `03_R_functions/bss_effort_spec.R` so the BSS and PE always share a unit. The CPUE denominator is the number of gear deployments (`h = number_of_gear`), and the daily effort expansion is `L = tau`, the deployment turnover (how many times gear is pulled and reset per day). Effort in deployment units is formed as `E = lambda_E * E_scale * L`, and the run asserts that `E` and `h` carry the same unit before sampling.

**Shore.** `lambda_E` is the number of crabbers, related to the dock gear count through the gear-per-crabber ratio `R_G` (learned from interviews). The effort scale is `E_scale = R_G` and the expansion is `L = tau_shore` (prior mean about 1.7), so `E = lambda_E * R_G * tau_shore` is in gear-deployments, matching `h`. This replaces the earlier shore construction, in which effort was crabber-hours formed as `(gear counted / R_G) * day_length`. Shore no longer expands by a civil-twilight day length under the deployment unit.

**Private boat.** `lambda_E` is the number of gear units in the water, related to trailer counts through `lambda_E / R_G_boat` (boat groups). The effort scale is `E_scale = 1` (`lambda_E` is already gear) and the expansion is `L = tau_boat` (prior mean about 1.2), so `E = lambda_E * tau_boat` is in gear-deployments. `R_G_boat`, the gear per boat group, is learned from observed `number_of_gear` in boat interviews.

**Commercial/charter.** A day-type-stratified census expansion of the vessel tally (weekday and weekend harvest rates expanded independently), not a BSS fit.

**Why deployments, and what it supersedes.** Pots and traps do not fish harder the longer they soak; a pot lift returns roughly the same catch whether it soaked two hours or eight (Section 11). Two earlier boat formulations were superseded. The crabber-hours-with-day-length version underestimated boat catch by roughly a factor of two, because pots keep fishing while the party is away but a civil-twilight day length is only 9 to 17 hours. The gear-hours version (`L = 24`) fixed that underestimate but over-corrected: with catch nearly flat in soak time, a fixed 24-hour denominator inflates effort for gear that is checked and re-set several times a day. The deployment scale removes both biases by denominating effort in pot lifts, on which CPUE is a stable rate. The same reasoning moved shore off crabber-hours (Section 15). The effort unit is set by `shore_effort_unit` in `run_config.R` (default `"gear-deployments"`; set to `"crabber-hours"` to revert shore only), with the turnover priors `tau_shore_prior_*` and `tau_boat_prior_*`. `tau` is identified by I/E ingress counts when enough days fall in the window, and otherwise rests on its prior (Section 16, GR-12).

**Assumptions.** Dock gear count is proportional to total shore effort; `R_G` and `R_G_boat` are constant within a sub-season; every crabbing boat trailers a vehicle and trailer counts capture all private boat crabbing; jetty and beach effort is captured via the dock count as a spatial proxy; mean catch per commercial vessel is constant within each day-type stratum.

### 14.3 Gear-type CPUE processes

```
log(lambda_C_gear[d, g]) = mu_C_gear[g] + omega_C_gear[period(d), g] + B1_C * w[d]
```

Each gear type `g` has its own baseline catch rate `mu_C_gear[g]` and its own AR(1) temporal deviation `omega_C_gear`. The design gives each gear its own process-error scale `sigma_eps_C_gear[g]` and its own negative-binomial overdispersion `r_C_gear[g]`, so pot CPUE can evolve smoothly (pots integrate over soak time) while snare CPUE can be more volatile (snares depend on tide, weather, and skill), and pot catches can carry a different variance structure than ring-net catches. All gear types share the AR(1) autocorrelation `phi_C_gear`. The initial AR state is non-centered (B1.3), as for effort.

`B1_C` is the weekend/holiday CPUE effect, controlled by `estimate_B1_C` (default TRUE, parity port B1.9; matching the pooled model). Setting it FALSE drops `B1_C` from the CPUE likelihood, reproducing the earlier v5.4 behavior. The effect is motivated by evidence that weekend and holiday crabber populations at tourist-accessible ports include more novice participants (Thomson 1991; Pollock et al. 1997), so weekend catch rates can differ from weekday rates.

**A caveat that governs this whole subsection (GR-7).** The Stan model is written for a per-gear CPUE process (per-gear `mu_C_gear`, per-gear scales, an AR that runs over the gear-by-section blocks with a correlation across them, and an interview catch likelihood indexed by gear type). As currently driven the R pipeline passes a single gear dimension (G = 1), so all of this per-gear machinery is inert: the model behaves as a pooled-CPUE model with a gear dimension of length one, and gear-type catch is apportioned after the fact. The per-gear structure below is therefore the model's design, not a currently active decomposition. See Section 16.

**Assumptions.** Each gear type's CPUE evolves smoothly; gear types share the degree of temporal autocorrelation but have independent process variability; CPUE does not vary by day type within a period beyond the shared `B1_C` weekend/holiday effect.

### 14.4 Weighted gear-type classification of interviews

Before classification, interviews are filtered for quality: only completed trips with at least `min_fishing_time` (default 0.5) crabber-hours are included. Gear types prohibited by regulation during the sub-season are then removed from each interview's gear list and the fractional weight is redistributed to the remaining reported types (Section 5). Gear types are detected with word-boundary regular expressions (`\bpot\b`, `\bring\s*net\b`, `\b(trap|star)\b`, `\bsnare\b`).

Each interview contributes in proportion to how many eligible gear types the crabber reported. For single-gear interviews (about 70%), the interview contributes fully to that gear type. For multi-gear interviews (about 30%), the interview's catch is split equally across the reported types. A crabber reporting "Pot, Ring Net" with 6 crab kept in the all-gear sub-season contributes weight 0.5 to Pot and 0.5 to Ring Net; the same interview in the ring-net sub-season contributes fully to Ring Net (Pot excluded). Equal weighting is a reasonable approximation when per-gear catch totals are not recorded; a naive first-listed or highest-rate assignment would systematically inflate one type and deflate the others, and with about 28% of interviews multi-gear that bias would meaningfully distort the decomposition.

These weights build the `gear_weights` matrix and, aggregated, the `pi_gear` proportions (Section 14.5). In the per-gear CPUE design the catch likelihood is a weighted pseudo-likelihood in which each interview informs all of its reported gear types; as driven with G = 1 the catch likelihood collapses to a single pooled term (every interview maps to the one gear dimension).

**Assumptions.** Equal weighting is a reasonable proxy when per-gear catch is not recorded; at least 15 effective interviews per gear type are needed for independent modeling, below which the type is collapsed; interviews with unrecognized gear types are assigned to the most common type.

### 14.5 Gear-type proportions and their uncertainty

The fraction of crabbers using each gear type, `pi_gear`, is computed from interview data as empirical proportions stratified by period and day type, with Laplace smoothing (alpha = 1) to prevent zero proportions:

```
pi_gear[period, day_type, g] = (n_weighted[period, day_type, g] + 1) / (N_weighted[period, day_type] + G_gear)
```

When a period by day-type cell has no interviews, the model falls back to period-level, then sub-season-level, proportions. These point-estimate proportions are used to allocate crabbers across gear types. But `pi_gear` carries genuine sampling uncertainty, especially in sparse strata (holiday cells, early or late season), and ignoring it would make the gear-type intervals too narrow.

**Propagation.** In the generated quantities block, gear-type proportions are drawn from a Dirichlet posterior on each MCMC iteration, using the conjugate Dirichlet-multinomial relationship (Gelman et al. 2013, Ch. 3) with the raw weighted counts plus the Laplace constant as concentration parameters:

```
pi_draw[g] ~ Dirichlet(n_weighted[period, day_type, g] + 1,  for g = 1..G_gear)
```

In data-rich strata the draws are tightly concentrated around the point estimate; in sparse strata they vary substantially, widening the gear-type intervals to reflect the limited data.

**Assumptions.** Gear proportions vary by period and day type; in sparse strata the model borrows strength through the fallback hierarchy and the Dirichlet sampling reflects the reduced information through wider intervals; Laplace smoothing is adequate regularization.

### 14.6 Observation models

- Gear counts (shore): `Gear_I ~ NegBinomial2(lambda_E * R_G, r_E)`, where `lambda_E` is crabbers and `R_G` is gear per crabber.
- Trailer counts (boat): `T_I ~ NegBinomial2(lambda_E / R_G_boat, r_E)`, where `lambda_E` is gear in the water and `lambda_E / R_G_boat` is boat groups.
- Gear per crabber (shore): `Gear_A ~ Poisson(A_A_gear * R_G)`, learning `R_G` from interview gear counts.
- Gear per boat group (boat): `Gear_A_boat ~ Poisson(R_G_boat)`, learning `R_G_boat` from interview gear counts.
- Interview catch: `c ~ NegBinomial2(lambda_C_gear[g] * h, r_C)`, where `h = number_of_gear` (deployments). In the per-gear design each interview contributes to its reported gear types with fractional weight; as driven with G = 1 this is a single pooled term.
- I/E ingress: an ingress count constrains the turnover `tau` through a lognormal observation of the predicted trips; for boats the prediction is on the group scale (`lambda_E / R_G_boat`).

The negative binomial accommodates the overdispersion typical of recreational trip-level catch and count data (Maunder and Punt 2004; Hilbe 2011).

### 14.7 Effort overdispersion (marginalized)

Each effort count is negative binomial with shape `r_E`. This was originally written as a Poisson-Gamma mixture with an explicit per-observation latent multiplier (a sparse parameterization that allocated a multiplier only for actual observations). Because the Gamma-Poisson mixture integrates exactly to the negative binomial (Hilbe 2011), the latent multipliers are marginalized analytically and the negative binomial is written directly (parity port B1.5). The change is inference-preserving (identical mean and variance), removes a high-dimensional latent block from the sampler, and makes the model block consistent with the `log_lik` block. The former per-observation count field is retained only for R-interface compatibility.

### 14.8 I/E integration and effective day length

On I/E survey days, an observed ingress count enters as a lognormal constraint on the predicted trips (`lambda_E * L` on the shore/deployment scale, `lambda_E / R_G_boat * tau` on the boat group scale), a second independent constraint on the latent effort state that identifies the turnover `tau`. When no I/E data fall in a fit's window, the I/E likelihood contributes nothing and `tau` rests on its prior. The prior on the I/E scale `sigma_IE` is applied unconditionally (parity port B1.6), so that with no I/E data `sigma_IE` is proper rather than an improper flat direction; because it is decoupled from effort and catch in that case, this leaves those posteriors unchanged.

For shore, an effective-day-length regression (log effective day length on a day-of-year quadratic and day type, fit from the I/E data) supplies the historical `L_effective` and remains available as a diagnostic; under the gear-deployment unit the shore expansion is the turnover `tau_shore` rather than an hours-based day length. Civil twilight (computed daily from `suncalc` at the Westport centroid, 46.904 N, 124.105 W, clamped to 9 to 17 hours) is retained only as a diagnostic column and as the last-resort fallback rung when there is effectively no I/E data.

### 14.9 Daily catch by gear type

The model produces two complementary representations of daily gear-type catch.

**Expected daily catch** (`C_gear`), the deterministic product of the latent effort, the expansion, the gear-type proportion, and the gear-type CPUE for each MCMC iteration:

```
C_gear[d, g] = lambda_E[d] * E_scale * L[d] * pi_gear[period(d), day_type(d), g] * lambda_C_gear[d, g]
```

This uses the point-estimate `pi_gear` and adds no stochastic sampling; its spread across iterations reflects effort and CPUE process uncertainty. It is used for the smooth daily trajectories in plots and for day-by-day BSS-vs-PE comparison.

**Predicted daily catch** (`C_gear_pred`), which adds gear-proportion sampling and count noise on top of the expected value:

```
C_gear_pred[d, g] = Poisson_draw(lambda_E[d] * E_scale * L[d] * pi_draw[g] * lambda_C_gear[d, g])
```

where `pi_draw` is a Dirichlet sample for that day's period by day-type stratum (Section 14.5). Season totals (`C_sum`, `C_sum_gear`) and their credible intervals are computed from `C_gear_pred`, so the reported bounds account for effort, CPUE, gear-type allocation, and count noise. Using only the expected value for totals would understate uncertainty by omitting allocation variance; using only the predicted value for daily plots would produce noisy trajectories that obscure the trend.

### 14.10 Key parameters and priors

| Parameter | Description | Prior |
|---|---|---|
| B1 | Weekend effort multiplier (log) | Normal(0, 1) |
| B2 | Additional holiday effort multiplier (log) | Normal(0, 1) |
| B1_C | Weekend/holiday CPUE effect (log), when `estimate_B1_C` is TRUE | Normal(0, 1) |
| R_G | Gear per crabber (shore: learned; boat: not used) | Lognormal(log(R_G_empirical), 0.3) |
| R_G_boat | Gear per boat group | Lognormal(log 4, 0.5), unconditional |
| tau_shore, tau_boat | Deployment turnover (shore about 1.7, boat about 1.2) | Lognormal from `tau_*_prior_*` |
| phi_E, phi_C_gear | AR(1) autocorrelation (phi_C_gear shared across gears) | Beta(2, 2) rescaled to [-1, 1] |
| r_E, r_C_gear | Overdispersion (effort; per-gear catch by design) | Half-Cauchy(0, 2) via sigma parameterization |
| mu_C_gear[g] | CPUE intercept per gear type (log) | Normal(log 0.5, 2) |
| sigma_eps_C_gear[g] | Per-gear CPUE process-error SD | Half-Cauchy(0, 2) |
| sigma_IE | I/E measurement error (log) | Exponential(5), unconditional |
| pi_gear, n_weighted_gear | Gear proportions and their Dirichlet concentrations | DATA (empirical, Laplace smoothed) |

Prior rationale. `R_G` is centered on the empirical gear-per-crabber ratio in the relevant population by sub-season, eliminating prior-posterior conflict. `R_G_boat` is given a proper unconditional Lognormal(log 4, 0.5) prior, centered on about 4 gear per group; this prior is applied in every fit, including shore fits where `R_G_boat` enters no likelihood, because a bounded-below parameter with no prior and no likelihood carries an unbounded Jacobian direction that otherwise runs away (fix marker F1). The half-Cauchy(0, 2) variance priors are weakly informative; on the log scale, scale 2 implies a plausible factor of about `exp(2) = 7.4` change between periods, which accommodates genuine seasonal variation while preventing implausible magnitudes in data-sparse strata (Gelman 2006).

### 14.11 Generated quantities and Stan file status

- `C_expected_sum` is the deterministic seasonal total used by the convergence gate (parity port B1.8), summed from `lambda_E * E_scale * L * lambda_C_gear` over days and gear.
- `C[d, g] = Poisson_rng(...)` gives predictive daily draws including Poisson noise, reported separately for prediction intervals; `C_gear` and `C_gear_pred` are the trajectory and total quantities of Section 14.9.
- Pointwise `log_lik` for the gear, trailer, and catch streams is produced, enabling PSIS-LOO.
- A single-cell hierarchical mean layer is collapsed automatically when the gear-by-section dimension is one (`use_mu_hier`, fix marker P2), removing an unidentified funnel; its now-inert scale parameters are a cosmetic reporting item (Section 16, GR-17).

**A note on the Stan file version.** The Stan model file `crab_bss_gear_resolved.stan` carries no vX.Y tag, and the older documentation label "Stan v3.2" is stale: the file is materially past that point. It now carries the pooled-track parity ports B1.3 (non-centered AR initial state), B1.5 (effort overdispersion marginalized to NB2), B1.6 (unconditional `sigma_IE` prior), B1.8 (the deterministic gate total), and B1.9 (the `B1_C` weekend/holiday CPUE effect), together with the run-driven fixes F1 (the unconditional `R_G_boat` prior), F2 (the gear-deployment effort scale), P1 (the configurable effort unit with `effort_scale_gear` / `E_scale`), and P2 (the single-cell mu-hierarchy collapse). Read the framework version tag for what the R pipeline does and the fix-marker record in the development history for what the Stan model does; do not attach a version number to the Stan file. The full marker record is in `BSS-GH-gear-type-CPUE-model-development-history.md`.

------------------------------------------------------------------------

## 15. Design decisions and their rationale

- **Two sub-seasons** are estimated independently because the pot closure creates a structural break in both effort and catch rates; pooling across it would blur two different regimes.
- **A per-gear CPUE architecture** is the point of this model: gear regulations are a primary management lever, and a modeled per-gear catch rate with its own uncertainty is what supports regulatory questions that a proportional split cannot. The architecture is built even though the production configuration runs it inert at G = 1 (Section 16).
- **Regulatory gear exclusion and a minimum effective-N** encode known fishery structure so the sampler is never asked to fit a CPUE process for a prohibited or near-empty gear type, which is a catastrophic-failure mode rather than a mild one.
- **Weighted fractional gear assignment** for multi-gear interviews, rather than a first-listed or highest-rate assignment, because about 28% of interviews are multi-gear and a hard assignment would systematically bias the decomposition.
- **Dirichlet propagation of gear proportions** so the gear-type intervals reflect allocation uncertainty in sparse strata, not just effort and CPUE uncertainty.
- **Gear-deployments for both components.** Effort is denominated in gear-deployments (pot lifts), not soak-hours or crabber-hours. The saturation diagnostic shows catch per gear-hour falls about 43-fold across soak durations while catch per pot lift is nearly flat, so any time-denominated unit violates the likelihood's proportionality assumption (`beta_h = 1`). The deployment is the unit on which CPUE is a stable rate, which keeps harvest = effort x CPUE unbiased. Crabber-hours failed the same test for shore (`beta_h = 0.57`), so shore moved onto deployments too (framework v5.5). The choice prioritizes harvest-unbiasedness over marginal predictive fit: gear-hours had a slightly better catch-stream `elpd_loo`, but that edge comes from the CPUE process absorbing the sub-linearity, which is exactly what biases the season expansion.
- **A weekend/holiday CPUE effect (B1_C)** because weekend and holiday crabber composition differs (more novices), which the catch-rate process should be allowed to reflect; default on, matching the pooled model.
- **Expected catch, not a predictive draw, as the headline** because harvest estimation wants E[C | data], the estimation quantity, not a single noisy prediction.
- **Shared modules and parity ports** so the gear-resolved and pooled tracks cannot drift onto different effort scales, gate policies, or convergence fixes: the convergence gate, the AR selector, the CPUE diagnostics, the effort specification, and the I/E day-length model are one implementation each, and the pooled convergence work is imported as the B1.x markers rather than re-derived.

------------------------------------------------------------------------

## 16. Limitations and future directions

**The standing caveats, stated honestly.**

- **G = 1, so per-gear catch is apportioned, not modeled (GR-7).** The Stan model supports a genuine per-gear CPUE process, but the R driver passes a single gear dimension, so the per-gear machinery is inert and gear-type catch comes from PE apportionment of interview proportions. The gear dimension cannot simply be raised: the only effort observation touches one gear dimension while the seasonal totals sum over all gear, so unmodeled effort processes would enter the totals identified by priors alone. The fully modeled path ("Option A") needs an effort-share offset fed by `pi_gear_data` and a rule for multi-gear interviews (the gear-weight matrix gives fractional weights, but a per-interview catch likelihood needs one gear per interview). Until it is built, read the per-gear split as an apportionment.
- **Boat catch rests on a thin `tau_boat` (GR-12).** `tau_boat = 1.2` is anchored on only two within-window boat I/E days (turnover 1.00 and 1.29). Boat catch is directly proportional to it. The boat I/E stream activates automatically once enough days fall inside a window (`ie_min_obs_boat`, default 2), but until it does, the boat total leans on this prior. State this dependence in any co-management number.
- **The gate backstop can override the impact test (GR-11).** The hard 5% divergence-fraction backstop can force a good BSS fit to fall back to PE even when the divergences demonstrably do not move the totals (Section 9). The backstop is kept on principle, but its 5% threshold and its precedence are an open policy question, to be revisited only after the PE fallback is validated.
- **`private_boat_ring_net` never fits (GR-15).** With about 17 interviews across September to November, this component has insufficient data and always falls back to PE.
- **Effort variance sits mostly in the observation term (GR-14).** Under the coarse fixed AR the latent process carries only about 4 to 13% of effort variance and the NB observation term about 81 to 91%. This is by design, not a defect, but it is why the AR contributes little to the effort fit.
- **A collapsed-layer cosmetic (GR-17).** When the single-cell mu-hierarchy is collapsed (P2), its scale parameters sample their prior and enter no likelihood; they can look non-converged in `structural_params` but are inert. They should be dropped or flagged when the collapse is active so a reviewer does not misread them.

**General limitations.**

- Multi-gear interview assignment uses equal weighting; true per-gear catch data would improve precision. A latent catch-allocation model would replace the weighted pseudo-likelihood but adds complexity.
- The shared `phi_C_gear` assumes similar temporal autocorrelation across gear types (the per-gear process-error scale partially addresses this).
- The commercial/charter census lacks formal uncertainty quantification.
- There are no jetty effort counts; beach crabbing is unmeasured within the pooled shore count, and the ring-net sub-season has single-count days.
- I/E coverage is limited; expanding it would sharpen both the `tau` turnovers and the shore day-length regression (Pollock et al. 1997 recommend at least 3 I/E days per month by day-type stratum).

**Future directions.**

- Build the modeled per-gear CPUE path (Option A): the effort-share offset and the multi-gear rule that GR-7 requires.
- Add Westport Jetty direct effort counts and expand the framework to other Washington coastal ports (re-deriving the site-specific priors).
- Allow an independent `phi_C_gear` per gear type if the data support it.
- Collect per-gear catch in interviews to replace equal weighting.
- Add bootstrap or parametric uncertainty for the commercial/charter census.
- Make the mu-hierarchy collapse a configuration lever (GR-10), as the pooled track did with `collapse_mu_hier`.

------------------------------------------------------------------------

## 17. Weather and tide covariates: not used here

A weather-and-tide covariate module exists in the pooled track (`06_diagnostics/`) as a research diagnostic. It was built and run on the 2024-25 season to test whether weather improves the estimate, and the conclusion, documented in `07_documentation/WEATHER_COVARIATE_ANALYSIS.md`, was to exclude covariates for all components under a pre-committed PSIS-LOO improvement margin: weather drives effort (not CPUE), and on routine data the AR(1) process already does the gap interpolation that weather could otherwise help with, so the covariate models were narrower but predicted held-out data worse. The gear-resolved model does not use the weather-tide module, and `run_weather` is only valid with the pooled model. The module is retained as a pooled-track contingency for seasons with extended sampling gaps under anomalous weather, evaluated by leave-one-week-out block cross-validation.

------------------------------------------------------------------------

## 18. Glossary

| Term | Meaning |
|---|---|
| BSS | Bayesian State-Space model |
| PE | Point Estimator |
| CPUE | Catch Per Unit Effort; as of v5.5 the denominator is gear-deployments (crab per gear deployment) for both shore and boat |
| gear-deployments | Count of gear set (pot lifts); the current CPUE denominator `h = number_of_gear` and the effort unit for both components |
| tau (tau_shore, tau_boat) | Deployment turnover: trips per present group per day; the daily effort expansion `L` (shore about 1.7, boat about 1.2), a parameter |
| E_scale | Converts `lambda_E` to the unit of `h`: R_G for shore (crabbers to gear), 1 for boat (already gear) |
| crabber-hours | One person crabbing for one hour; the pre-v5.5 shore CPUE denominator, superseded by gear-deployments |
| gear-hours | Total time gear spent in the water; a superseded boat denominator, invalid for pots (catch is flat in soak time) |
| AR(1) | First-order autoregressive process |
| period(d) | Mapping from day d to its AR period index (biweekly ring-net, monthly all-gear) |
| R_G | Gear-per-crabber ratio (learned for shore) |
| R_G_boat | Gear-per-boat-group ratio (learned for boat); replaces the retired R_T |
| Gear_A_boat | Observed gear per boat-interview group; the data that informs R_G_boat |
| B1 / B2 | Weekend / additional-holiday effort multipliers (log) |
| B1_C | Weekend/holiday CPUE effect (log); controlled by `estimate_B1_C` (default TRUE) |
| G_gear | Number of gear types set up in a sub-season (as driven, the fitted gear dimension is 1) |
| pi_gear | Gear-type proportions per period by day_type (data-derived point estimate, Laplace smoothed) |
| pi_draw | Dirichlet-sampled gear proportions used in prediction quantities to propagate allocation uncertainty |
| n_weighted_gear | Raw fractional interview counts per stratum; the Dirichlet concentration parameters |
| gear_weights | Fractional gear-type assignment per interview (rows sum to 1) |
| lambda_C_gear | Gear-type-specific catch rate (crab per deployment) |
| sigma_eps_C_gear | Per-gear CPUE process-error standard deviation |
| r_C_gear | Per-gear catch overdispersion |
| C_gear | Expected daily catch by gear type (no Poisson noise); used for trajectory plots |
| C_gear_pred | Predicted daily catch by gear type (Dirichlet-sampled pi + Poisson draw); used for season totals |
| C_expected_sum | Deterministic seasonal total used by the convergence gate; E[C \| data] |
| gear_exclude | Per-sub-season list of gear prohibited by regulation; prevents fitting a CPUE for illegal gear and redistributes fractional weights |
| bss_min_gear_effective_n | Minimum sum of fractional gear weights for a gear type to be modeled independently (default 15) |
| Incomplete trip | Interview taken before the crabber finished; excluded from CPUE (systematic downward bias) |
| Divergent transition | HMC diagnostic indicating difficult posterior geometry; can bias estimates even when R-hat looks fine |
| n_eff | Bulk effective sample size |
| R-hat | Rank-normalized split potential scale reduction factor (convergence) |
| PSIS-LOO | Pareto-smoothed importance-sampling leave-one-out cross-validation |
| PIT | Probability integral transform (posterior predictive calibration) |
| I/E | Ingress/Egress survey |
| sigma_IE | I/E measurement error on the log scale |

------------------------------------------------------------------------

## 19. Development history (summary)

Method v1.0 corresponds to gear-resolved framework code **v5.5** (2026-07-11). The model began as an adaptation of the WDFW freshwater-creel state-space framework, shared a single milestone line with the pooled model through v4, and branched at v5.0 into an independent gear-resolved track. The arc in one screen:

- **v1-v4 (shared line):** the single-population prototype, shared bug fixes, the three-population two-sub-season pooled model with convergence tuning, and the day-length, stat-week PE, and census-date work.
- **v5.0:** the initial gear-resolved release: per-gear AR(1) CPUE processes with shared phi and sigma, Dirichlet `pi_gear` per period, the multi-gear priority hierarchy (later replaced), the B2 holiday effort effect, and sparse overdispersion.
- **v5.1-v5.2:** empirical proportions and per-gear structure. Weighted fractional gear assignment replaced the priority hierarchy; `pi_gear` became data with Laplace smoothing; per-gear process-error was added; then the incomplete-trip filter, Dirichlet-sampled `pi_gear` in generated quantities, `C_gear` vs `C_gear_pred`, the divergence gate, the regulatory `gear_exclude`, and the minimum effective-N raise from 3 to 15.
- **v5.3-v5.4:** the R-side data prep was aligned with the boat formulation the Stan model then expected (later superseded by the deployment scale), and the R-hat gate was tightened from 1.05 to 1.01 (Vehtari et al. 2021), in step with the pooled track.
- **v5.5 (2026-07-11):** shore was moved onto the gear-deployment effort unit to match the boat, so both fitted components now share one effort unit with each other and with the pooled track. Chosen from the 2026-07-10 shore LOO comparison as the only harvest-unbiased shore unit (`beta_h = 1.05`, covering 1, against 0.57 for crabber-hours and 0.73 for gear-hours), on harvest-unbiasedness rather than marginal predictive fit. Resolves backlog GR-16 and sets `loo_effort_unit_comparison = FALSE` for production. This entry also retroactively records that the boat had already moved to gear-deployments via the shared effort module (fix markers F1, F2, P1).

The run-driven Stan fix markers (F1 the unconditional `R_G_boat` prior; F2 the gear-deployment scale, the single largest quantitative move in the model's history; F4 the CPUE diagnostics; F5 reporting; P0 the point-estimator population and ratio-of-sums fix; P1 the configurable effort unit; P2 the mu-hierarchy collapse) and the pooled parity ports (B1.3, B1.5, B1.6, B1.8, B1.9) are recorded, with their before/after numbers and the F1/F2 and P0 working notes, in **`BSS-GH-gear-type-CPUE-model-development-history.md`**. The outstanding backlog (GR-7 to GR-17) lives in `07_documentation/development_notes/20260710-OUTSTANDING_ISSUES.md`.

**Reference-run figures, and a note on staleness.** The last full assembled totals on record are from the 2026-07-09 and 2026-07-10 reference runs, before the v5.5 shore deployment move. Because v5.5 changes the shore effort unit, the shore components must be refreshed by a v5.5 run before they are cited as current.

| Component | Reference figure | Basis and status |
|---|---|---|
| shore_ring_net | about 6,936 (BSS) | crabber-hours scale; PRE-REFRESH under v5.5 |
| shore_all_gear | about 19,684 (BSS) | crabber-hours scale; PRE-REFRESH under v5.5 |
| private_boat_all_gear | about 43,314 (BSS) | gear-deployment scale; current, but rests on the `tau_boat = 1.2` prior (GR-12) |
| private_boat_ring_net | PE fallback | insufficient data (GR-15) |
| commercial/charter | about 12,007 | day-type census; scale-independent |

The 2026-07-10 run was the first in which all three fitted components reported on their BSS posteriors on corrected effort scales with the PE as a consistent cross-check. The v5.5 shore move then requires a confirming run before a refreshed port total is published.

------------------------------------------------------------------------

## 20. References

Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo. *arXiv preprint* arXiv:1701.02434.

Betancourt, M. & Girolami, M. (2015). Hamiltonian Monte Carlo for hierarchical models. *In:* Current Trends in Bayesian Methodology with Applications. CRC Press.

Conn, P.B. (2002). Bayesian methods for estimating recreational angler effort, catch rates, and total catch using creel survey data. Ph.D. Dissertation, University of Wisconsin-Madison.

Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models. *Bayesian Analysis*, 1(3), 515-534.

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., & Rubin, D.B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

Hahn, P.K.J., Brooks, L., & Hartill, B.W. (2000). Strategies and procedures for estimating catch and effort in freshwater fisheries. *In:* Inland Fisheries Management in North America (2nd ed.), American Fisheries Society.

Hartill, B.W., Cryer, M., Lyle, J.M., Rees, E.B., Ryan, K.L., Steffe, A.S., Taylor, S.M., West, L., & Wise, B.S. (2012). Scale- and context-dependent selection of recreational harvest estimation methods. *North American Journal of Fisheries Management*, 32(1), 109-123.

Harvey, A.C. (1989). *Forecasting, Structural Time Series Models and the Kalman Filter*. Cambridge University Press.

Hilbe, J.M. (2011). *Negative Binomial Regression* (2nd ed.). Cambridge University Press.

Maunder, M.N. & Punt, A.E. (2004). Standardizing catch and effort data: a review of recent approaches. *Fisheries Research*, 70(2-3), 141-159.

Papaspiliopoulos, O., Roberts, G.O., & Skold, M. (2007). A general framework for the parametrization of hierarchical models. *Statistical Science*, 22(1), 59-73.

Pollock, K.H., Jones, C.M., & Brown, T.L. (1994). *Angler Survey Methods and Their Applications in Fisheries Management*. American Fisheries Society Special Publication 25.

Pollock, K.H., Hoenig, J.M., Jones, C.M., Robson, D.S., & Greene, C.J. (1997). Catch rate estimation for roving and access point surveys. *North American Journal of Fisheries Management*, 17(1), 11-19.

Robson, D.S. (1991). The roving creel survey. *American Fisheries Society Symposium*, 12, 137-148.

Staton, B.A., Catalano, M.J., Connors, B.M., Coggins, L.G., Jones, M.L., Walters, C.J., Fleischman, S.J., & Beardsall, J.W. (2017). Evaluation of methods for spawner-recruit analysis in mixed-stock Pacific salmon fisheries. *Canadian Journal of Fisheries and Aquatic Sciences*, 74(7), 1108-1122.

Sullivan, M.G. (2003). Active management of walleye fisheries in Alberta. *North American Journal of Fisheries Management*, 23(4), 1343-1358.

Thomson, C.J. (1991). Effects of the avidity bias on survey estimates of fishing effort and economic value. *American Fisheries Society Symposium*, 12, 356-366.

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413-1432.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.C. (2021). Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.
