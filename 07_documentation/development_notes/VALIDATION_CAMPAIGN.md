# The validation campaign, 2026-08 to 2026-09: the run-by-run record

- **What this is.** The dated, run-by-run narrative of the validation and improvement
  campaign that produced Method v2.0. It was Sections 1b through 1v of
  `PIPELINE_STATUS.md` until 2026-09-12, and it is moved here unaltered so that the status
  document can be what its name says: the current state and the backlog.
- **Cross-references to "Section 1x" anywhere in this repository resolve HERE.** The section
  letters and titles are unchanged.
- **Companions.** `PIPELINE_STATUS.md` is the living status and backlog.
  `CHANGE_REGISTER.md` is the tabular register of every change, its status and its evidence.
  The two `../BSS-GH-*-development-history.md` files are the version-by-version change logs.
  This file is the middle layer: why each decision was taken, on what run, and what the run
  actually said.

> ### A NOTE ON THE DATES, because they do not all agree
>
> **The section LETTERS are the order. The parenthetical dates are working-session dates and
> they drifted from the git commit dates in the 1r-to-1v range**, where several sessions'
> work landed in the same day's commits. Where a date matters, the authorities are the git
> log and the output folders, not a section title: by git, the work in Sections 1r, 1s, 1t
> and 1u all landed on 2026-09-10, and the full ladder run of Section 1v rendered overnight
> on 2026-09-10 / 2026-09-11 (`IMP_STAGE.txt` in each rung folder carries the exact
> timestamp). Two titles that read 2026-09-12 and 2026-09-13 have been corrected to
> 2026-09-10 for that reason; the rest are left as written. The same drift appears in dated
> comments inside `06_diagnostics/run_improvements_2026-09-08.R` and the harness, where
> "2026-09-13" marks work the git log dates 2026-09-10. Those are left alone: they are
> internally consistent with each other and several are load-bearing strings in assertions.
>
> **Every total below the first section is historical.** The authoritative run lives in one
> place, the box at the top of `PIPELINE_STATUS.md`.

---

## 1b. Stage 5 prerequisites (2026-08-30): what landed before the next batch

Everything below is code and reporting, applied and covered by the harness (**264 assertions** as of 2026-09-03, 0 failing). No estimate has changed; the batch that tests these is `06_diagnostics/run_stage5_2026-08-30.R`, which ships with `DRY_RUN <- TRUE`.

**1. The shared-turnover informed-day floor (plan 5.2).** `bss_shared_tau_data()` now refuses `shared_tau` for a fit with fewer than `shared_tau_min_obs` (default **15**) days that can inform `L`, and falls back to the per-day parameterization with a printed reason. Informed days = I/E days, plus OSP days when `osp_scale_is_tau = 1` puts `L` into the OSP mean. The observed counts are 4 (shore all-gear), 0 (shore pot closure), 130 (boat all-gear), 18 (boat pot closure), so any threshold in 5..18 separates the shore from the boat; 15 sits just below 18 deliberately, to keep the cross-track-corroborated boat pot closure. `shared_tau` remains **off in production** and stays off until the Stage 5 gate confirms the floor does exactly what it claims.

**2. Model adequacy beside the gate (plan 5.5).** `03_R_functions/bss_model_adequacy.R` writes `model_adequacy.csv` per run: `p_loo` as a fraction of `n_obs` (worst stream), the count of Pareto k > 0.7, the worst PIT bias, and the smallest `n_eff` among the observation-model dispersion parameters the gate never looks at. **It does not gate.** `pass_convergence` keeps its exact meaning; these columns sit beside it so "did this fit sample?" and "is this model adequate?" are visibly different questions. On the 2026-08-29 cells the separation is stark: stage E carries `p_loo` at 8.0% of `n_obs` with PIT bias 0.021, stage F carries 16.3% with bias 0.096, two bad Pareto k, and `sigma_r_OSP` at n_eff 307, below the gate's own floor of 400, with the gate passing.

**3. Retro annotation, so archived runs are comparable rather than blank.** `annotate_model_adequacy_run()` and `annotate_decoupled_run()` rebuild both tables for a committed run folder from its own CSVs, writing `model_adequacy_reconstructed.csv` and `decoupled_audit.csv` **beside** the originals; committed outputs are never rewritten. Both mark every row `source = "reconstructed"`. This exists because plan 5.3 asks for adequacy "for each cell" of a 2x2 whose archived cells predate the diagnostic, and a table with numbers for the new cells and blanks for the old ones is the shape of an argument that quietly favours whichever cells are new.

**4. Four reporting defects, all fixed.**

| Defect | Fix |
|---|---|
| `pe_vs_bss_comparison.csv` printed `BSS_catch = 32,689` for a component whose `method_selected` was `PE`; the rejected fit's `tau_bar` (n_eff 49) and `f_crab` (n_eff 25, R-hat 1.17) sat unmarked in `structural_params_*` | both drivers now write `bss_reported`; `structural_params_*` carries `fit_method` and an `estimate` column that is `median` with the decoupled rows blanked |
| `kappa_OSP` is decoupled in every fit under `osp_scale_is_tau = 1` (its ~3.0 median is its prior) and shore `r_OSP` has a posterior mean of 2.5 million | `bss_decoupled_reasons()` flags both with a stated reason; the new `estimate` column is `NA` for them so a prior cannot be read as a result |
| `expansion_ratios.csv` printed a decoupled shore `R_G_boat` (~4.0) with no flag | pooled table gained `R_G_boat_decoupled` and a `[PRIOR ONLY: no boat count stream in this fit]` suffix; the gear track already said `not applicable (shore fit)` |
| `structural_params_*` and `model_adequacy.csv` invented `"PE (gate fail)"` where `convergence_report.csv` says `"PE (convergence fail)"`, so one fit read two ways in one folder | all three writers now quote `gate_info$method_selected`, the gate being the single authority on method selection |
| the 2026-08-26 baselines have no `decoupled` column at all | `annotate_decoupled_run()` reconstructs it from the fit label, `fit_data_summary.csv` and `run_parameters.txt`, marking window-dependent rules `unknown` rather than guessing |

**5. `bss_sampler_override`, an explicit escape hatch (new).** Each driver merges `params_model` **on top of** `run_config`, so `params_model` wins every sampler key. Setting `bss_iter_default` in a `run_config` delta therefore does nothing and the run looks like it complied; the same class of trap that left `bss_min_interviews` pinned at 20 until the 2026-08-25 audit, and one that would have silently voided plan 5.2's request to re-run the gear track at more draws. `bss_sampler_override` is a named list applied **after** the merge, restricted to sampler keys, printing every change and **erroring** on anything else rather than dropping it. Production value is `NULL`.

**6. `fit_data_summary.csv` now records the turnover decision** (`n_L_informed`, `shared_tau`) per fit, so the floor's decision is auditable from the output folder rather than from the console.

**One thing measured rather than assumed (2026-08-30).** Extra elements in the list passed to `rstan::stan(data = )` are inert: a two-chain fit run with, and without, three added elements including the new `.n_L_informed` returned bit-identical draw matrices. The dotted metadata `prep_bss_crab_*` attaches to `stan_data` therefore cannot perturb any bit-identity gate.

---

## 1c. The Stage 5 batch (2026-08-31): the 2x2 is settled

Full review: `development_notes/stage5-batch-review-2026-08-31.md`. Six stages, 24 h of fitting, 11 of 12 criteria passing; the one FAIL was a criterion comparing against a reference that differed in two ways, not a code defect. Harness now **199 assertions**, including guards that recompute the interaction, the calibration ordering and the cross-track `tau_bar` agreement from the committed outputs, so the review's conclusion cannot drift from its evidence unnoticed.

### The finding

| private-boat all-gear | boat AR monthly | boat AR daily |
|---|---|---|
| **shared turnover OFF** | 25,868 | 37,359 |
| **shared turnover ON** | **31,008** | 42,344 |

The two levers are **additive** (interaction -155 against main effects of +5,140 and +11,491). They are not two measurements of one thing, and additivity adjudicates nothing. Calibration and complexity do, and they are unanimous:

| cell | trailer elpd | catch elpd | p_loo | k>0.7 | OSP coverage_50 | smallest sigma_r_* |
|---|---:|---:|---:|---:|---:|---:|
| OFF x monthly | -573.8 | -451.4 | 4.9 | 0 | 0.508 | 0.768 |
| **ON x monthly** | **-559.1** | **-451.5** | **8.0** | **0** | **0.508** | **0.767** |
| OFF x daily | -540.4 | -455.3 | 31.8 | 2 | 0.931 | 0.307 |
| ON x daily | -473.4 | -455.1 | 52.5 | 16 | 0.977 | 0.086 |

*Nominal coverage_50 = 0.500.* The shared turnover buys 14.7 nats of trailer elpd for 3.1 effective parameters, costs the catch stream 0.1, and moves every calibration statistic toward nominal. The daily AR buys elpd by spending 27-48 effective parameters on 195 observations, makes the catch stream 3.8 nats worse in both tau conditions, and destroys calibration. **Identification versus absorption, no longer a hypothesis.**

### Cross-track corroboration, the strongest the project has produced

With the gear track's boat all-gear fit finally converged (S3: n_eff 49 to 19,292, R-hat 1.081 to 0.9997, 554 divergences to 2, after raising adaptation as well as draws), the two tracks agree on `tau_bar` to **0.03%** (2.5969 vs 2.5962), on the boat component to 0.80%, and on the port total to 0.89%. The independent OSP/trailer overlap calibration puts the turnover at 2.01-3.03; the monthly posterior is 2.597 [2.06, 3.25], straddling it. The daily cell's 3.232 [2.87, 3.63] sits at or above its top.

### Recommended change to the authoritative run

**Adopt the shared turnover, boat-only** (`shared_tau = TRUE` with `shared_tau_min_obs = 15`). New authoritative run **S2, `05_output/20260829/pooled-CPUE-S5-2-tau-pooled`, port total 71,521 [52,350, 100,759]**, cross-checked by S3 at 70,886. **Do not adopt the daily boat AR.** *(Adopted 2026-09-01: `run_config.R` ships `shared_tau = TRUE` with `shared_tau_min_obs = 15`. Confirmed 2026-09-02 by V1. See Sections 1d and 1e.)*

### Two open items this batch created

1. ~~**The trailer stream is over-covered in every configuration.**~~ **RETRACTED 2026-09-01, see Section 1d.** It was a diagnostic artefact: `ppc_calibration_*.csv` computed `coverage_50` from a quantile interval of simulated draws, which over-covers small counts by construction. On the randomized PIT the production run is 0.523 against a nominal 0.500, 0.6 sampling SDs: calibrated.
2. **`prior_influential` measures variance reduction only.** `tau_bar` contracts 0.216 (flagged prior-influential) while its posterior mean sits **3.5 prior SDs** from the prior mean. The definition needs revisiting before that flag is quoted.

### Diagnostic defects the batch exposed, all fixed

- **`pit_worst_bias` alone ranked the cells backwards.** The (ON, daily) cell has the best PIT mean in the batch (0.495) and 98% of OSP observations inside a nominal 50% interval: a latent process with one state per observation interpolates the data and piles every PIT value at 0.5. `cov50_worst_dev`, `pit_sd_worst_dev` and `flag_miscalibrated` added; both statistics were already in `ppc_calibration_*.csv` and neither reached the table.
- **`disp_neff_min` cannot see a dispersion collapse.** `sigma_r_OSP` fell 0.806 to 0.086 while its n_eff stayed at 657, above the gate's floor. `disp_scale_min` added, reported **unflagged**: the half-Cauchy priors have no finite SD so contraction is undefined, and any cut-off would be chosen after seeing these runs.
- **REGRESSION, fixed.** `tau_bar` was registered in the prior-vs-posterior table under its bare name while rstan names the row `tau_bar[1]`; the lookup threw, the writer's `tryCatch` swallowed it, and the entire boat `prior_vs_posterior_*.csv` vanished from S2, S3 and S4b. Entry renamed, lookup made tolerant in both directions, unresolvable parameters now dropped rather than taking the file with them.

### Also settled

- **Plan 5.1 (the fill).** The re-run reproduces the driver exactly on zeroed days (44/0/47/9) and PE catch. `day_type` raises the PE 22.9-28.5% per component, **+14.0% at the port**. Decision still open; it moves a reported number only where a component falls back to PE, and all four reported BSS in every Stage 5 run.
- **Plan 5.6 (stage C re-run).** Boat pot closure 735 at biweekly (gear track 743, 1.1% apart) with the all-gear component back at rung 4's 25,868. The 2026-08-28 `ar_force` leak was worth **3,025 crab** on that component and 2,865 at the port.

---

## 1d. Adoption and correction (2026-09-01)

### The shared turnover is ADOPTED, boat-only

`run_config.R` now ships `shared_tau = TRUE` with `shared_tau_min_obs = 15` stated explicitly rather than left to the helper default. The feature is boat-only **by the informed-day floor, not by a population switch**: shore all-gear has 4 informed days out of 289 and shore pot closure has 0, so both degrade to per-day draws with a printed reason; boat all-gear (130) and boat pot closure (18) take the shared level. The five lines of evidence are recorded at the toggle itself and argued in `stage5-batch-review-2026-08-31.md`.

**Recommended authoritative run: S2, `05_output/20260829/pooled-CPUE-S5-2-tau-pooled`, port total 71,521 [52,350, 100,759]** *(adopted 2026-09-01 and confirmed 2026-09-02 by V1, which reproduced it across 11,223 shared parameter rows; the production run is now `20260831/pooled-CPUE-VAL-1-adopted` at 71,513, the 8-crab gap being the documented draw-permutation effect)*, cross-checked by the gear track at 70,886 (0.89% apart). It supersedes rung 4 (66,237). **Confirmation still outstanding:** stage V1 of `run_validation_2026-09-01.R` re-runs production as shipped and must reproduce S2 bit-for-bit. A config-key comparison shows no substantive difference, but a text diff of a `str()` dump is not a bit-identity proof, and S2 additionally predates the prior-vs-posterior fix so it carries no boat `prior_vs_posterior` file at all.

### A diagnostic defect that corrected the 2026-08-31 review

`ppc_calibration_*.csv` and `ppc_byobs_*.csv` computed `coverage_50` two different ways and disagreed by **up to 0.154** on the private-boat trailer stream. The aggregate file tested the observation against a quantile interval of simulated draws; for small counts that interval cannot carry 50% of the mass, so it over-covers by construction. The per-observation file had always used the randomized PIT. `model_diagnostics.R` now uses the randomized statistic too, so the two agree; `bss_model_adequacy.R` prefers `ppc_byobs` wherever it exists and records which source it used, which keeps archived runs comparable.

| stream (production run) | non-randomized (shipped until 2026-09-01) | randomized (correct) |
|---|---:|---:|
| private-boat trailer, n = 195 | 0.667 | **0.523** (0.6 sampling SDs) |
| private-boat catch, n = 131 | 0.595 | 0.595 |
| trailer, daily-AR cell | 0.836 | **0.744** (6.8 sampling SDs) |

Two consequences: the trailer over-coverage item is retracted, and the Stage 5 recommendation survives on a cleaner contrast. Two smaller fixes rode along: the OSP stream now gets per-observation PPC rows (it was the one stream whose calibration could not be cross-checked, and it is where the daily-AR pathology is most extreme), and `p_zero` is now written per observation so the zero bin can be scored against all days rather than against the selected zero days.

### What the four desk stages of the validation batch retired, with no MCMC

- **Tier 1, "TOP OF THE LIST", the shore I/E observation-unit fix: CLOSED, all four pre-set criteria met.** The isolating run has existed since 2026-08-26 (ladder rung 1 to rung 2) and nobody had scored it. (a) shore all-gear `sigma_IE` 1.034 [0.64, 1.68] to 0.341 [0.02, 0.71], and 0.272 today; (b) shore all-gear BSS 20,708 to 21,017, +1.5%; (c) the shore `L` posterior span widens 0.354 to 0.873 and can now go BELOW the 1.700 prior centre, where before it could only go up; (d) the boat is bit-identical across 3,158 shared parameter rows while the shore differs, which is exactly the pattern the criterion demanded. **This also closes the Tier-2 GR-9 `sigma_IE` tension**, which the backlog said must not be closed on reasoning alone.
- **Tier 2, the `gear_only` incomplete-trip arm: decidable, and the decision is about consistency rather than accuracy.** Adopting it shifts `R_G` by at most 0.7% (shore all-gear) and 0.1% (boat). The length-bias test is significant on shore all-gear (p = 0.0006) only because n = 2,741 makes a 0.7% shift detectable, and the shift is DOWNWARD, not the upward bias the rule was written against. **The live defect is that the PE and the BSS disagree**: both boat components run the BSS on `exclude` and the PE on `gear_only`.
- **Tier 3, zero inflation: the precondition is now computable** (it was not, before `p_zero`). D3 scores it against all days once a post-2026-09-01 run exists.
- **Tier 4, hygiene: `.Rproj.user` is still tracked (5 files); `05_output` is 6,837 tracked files and 563 MB; the stale 20260711 morning/afternoon runs are 222 files. The LICENSE is 35,149 bytes with no placeholder, so the GPL-3.0 paste is done** and only WDFW's confirmation of the copyright lines remains.

---

## 1e. The validation batch (2026-09-02): adoption confirmed, and a production fit that never worked

Full review: `development_notes/validation-batch-review-2026-09-02.md`. Nineteen criteria, 14.1 h of fitting. Harness **233 assertions**.

### Adoption CONFIRMED

V1 reproduces S2 across **11,223 shared parameter rows, identical at full double precision**, all four components matching. The shipped configuration is the one the Stage 5 batch measured. **Authoritative run: `05_output/20260831/pooled-CPUE-VAL-1-adopted`, port total 71,513 [52,346, 100,742].** The boat `prior_vs_posterior` files exist again, closing the 2026-08-31 regression.

### THE PRODUCTION DEFECT: the gear track has never been able to fit its boat component

The gear driver had per-fit sampler settings for shore all-gear and both pot-closure fits and **none at all** for boat all-gear, which fell through to the track defaults (2,000/1,000 draws, `adapt_delta` 0.90, `max_treedepth` 10) while the pooled track fits the same component on the same data at 5,000/2,500, 0.99, 13.

| | divergences | `tau_bar` n_eff | R-hat | gate | port |
|---|---:|---:|---:|---|---:|
| gear, production defaults (V2, V3, and 2026-08-29) | 554 | 49 | 1.0813 | REJECTED, falls back to PE | **51,385** |
| gear, with the settings (S3) | 2 | 19,292 | 0.9997 | BSS | **70,886** |

**A 27% understatement of the cross-check that validates the pooled headline.** The Stage 5 batch found the fix and left it in an experiment-only escape hatch, so production never received it; the validation batch reproduced the failure bit-identically. **Fixed 2026-09-02**: the four settings now live in the gear driver's `params_model` with their own branch. `bss_sampler_override` stays NULL in production.

### The finding the batch was not looking for: production's SHORE all-gear fit is the flagged one

| V1 fit | p_loo / n_obs | Pareto k > 0.7 | coverage_50 dev | miscalibrated |
|---|---:|---:|---:|---|
| **shore all-gear (DAILY AR)** | **35.2%** | **41** | **0.201** | **YES** |
| boat all-gear | 8.0% | 1 | 0.095 | no |

Shore all-gear carries 109.4 effective parameters on 311 gear observations and `coverage_50` 0.701 against a nominal 0.500, **+7.1 sampling SDs**: the same signature the daily boat AR was rejected for, on a component that is 29% of the port total and has been at daily in production all along. The gear track's monthly fit of the same component is not better but differently wrong (coverage 0.035, -16.4 SD), so neither end is right and nobody has run weekly or biweekly. That cross-track table differs in three things at once and is suggestive only; the pooled daily number alone is sound.

### Other results

- **GR-7 Phase 2 SAMPLED for the first time and works**: every per-gear interval widens 1.12x-1.23x, no median moves as much as 1%, per-gear medians still sum to the component total. Validated; whether `gear_resolved_G` becomes the gear default is now a decision.
- **The boat pot-closure 2x2 completes and is additive**: 849 / 1,018 / 735 / **939**, interaction +35 crab. `ar_force` leak stays closed (all-gear back at 31,008 exactly).
- **The `shared_tau_min_obs` threshold is worth 169 crab, 0.24% of the port** (V5). Defensible and nearly inconsequential.
- **Tier 1 "TOP OF THE LIST" shore I/E fix CLOSED** on all four pre-set criteria, and the Tier-2 GR-9 tension with it. Details in Section 1d.
- **Tier 3 zero inflation: precondition met, and it targets ONE stream.** Every boat stream is fine (|z| <= 2.3); the two shore CATCH streams under-predict zeros at z = +3.8 and +2.9, about 71 and 27 extra zero-catch interviews. A prototype should target the shore catch likelihood, not the whole model.
- **`gear_only` arm**: worth at most 0.7% anywhere. The live defect is that the PE and BSS use different arms on both boat components.

### A process failure worth recording

Two of the last three batches contained a criterion that compared against a reference differing in **two** things, both times producing a FAIL that looked like a code defect and needed manual diagnosis (the S3 shore criterion on 2026-08-30, the V2 boat criterion on 2026-09-01). `fit_exactness()` now takes `expect_delta`, parses both runs' `run_parameters.txt`, and reports any unexpected config difference inside the verdict. Judgement was not reliable here; it is now a check.

---

## 1f. The next batch (2026-09-03): shore AR ladder, gear cross-check, arm alignment, ZINB prototype

`06_diagnostics/run_shore_ar_zi_2026-09-03.R`, one PE-only stage and six fitted, roughly 23-26 h. Ships `DRY_RUN <- TRUE`. Harness **279 assertions**.

### Two code changes, both applied and both testable

**1. PE / BSS incomplete-trip arm alignment (Tier 2, closed).** The BSS learns `R_G_boat` from a frame that HAS been incomplete-trip filtered; both Point Estimators used the unfiltered one. New shared helper `03_R_functions/pe_gear_ratio_frame.R`, called by both PEs, with `pe_gear_ratio_arm = "match_bss"` shipping as the default and `"gear_only"` retained for the four-arm diagnostic and historical reproduction. **Measured (stage P1, PE only, seconds): the boat PE moves 11,176 to 11,170, -0.059%, and the shore is unchanged as it must be.** A consistency fix, not an accuracy one; the reason to make it is that there is no estimate to defend, only an inconsistency a reviewer asks about first. One helper rather than two edits because the gear PE's comment already claimed it matched the BSS while it did not.

**2. Zero-inflated catch likelihood (Tier 3, prototype, ships OFF).** `crab_bss_pooled.stan` gains `zi_catch`, a Beta prior and `vector<lower=0,upper=1>[zi_catch] theta_C`, so the OFF path is zero-size and bit-identical. Targeted from the 2026-09-01 zero bin rather than from taste: shore all-gear 676 observed zeros against 605.4 expected (z = +3.8, n = 1,649), shore pot closure z = +2.9, every boat stream inside |z| = 2.3.

> **The season total is scaled by (1 - theta_C) in generated quantities, deliberately.** `lambda_C` is fitted to the non-inflated component and rises to absorb the zeros `theta_C` removes, so reporting `lambda_C * effort` unscaled would inflate the total by 1/(1 - theta_C) purely by turning the feature on. Judge the prototype on the catch stream's `elpd_loo` and on the zero bin, never on the total.

It is scoped **per fit** (`catch_zi_populations`, default `"shore"`), so one run carries the shore fits as treatment and the boat fits as an untouched negative control at no extra cost. The model compiles to C++ in 3.1 min.

### The stages

```
P1  PE arm alignment ..............  seconds, PE only        DONE in the dry run: -0.059% on the boat
 |
G1  gear track with the driver fix   ~35 min                 restores the cross-check (51,385 -> ~70,900)
 |
Z0  ZINB code present, OFF ........  ~4.5 h    GATE, and the ladder's DAILY rung
 |
A1  shore all-gear AR weekly ......  ~4-5 h    }  THE MAIN QUESTION
A2  shore all-gear AR biweekly ....  ~4-5 h    }  four rungs, one component, one thing changing
A3  shore all-gear AR monthly .....  ~4-5 h    }
 |
Z1  ZINB ON, shore only ...........  ~5 h      boat fits are a free negative control
```

**Read the ladder in this order**, the same order that settled the boat 2x2: coverage_50 from the randomized PIT first, then `p_loo` as a fraction of `n_obs` with the Pareto k count, then the CATCH stream's elpd. Production sits at daily with p_loo 35.2% of `n_obs`, 41 Pareto k above 0.7 and coverage 0.701 (+7.1 sampling SDs); the gear track's monthly fit of the same component gives 0.035 (-16.4 SDs). Both ends are bad and the middle has never been run. If the winner is not daily, the `ar_max_resolution` entry changes and the shore component of the published estimate moves with it, so it is a number-moving decision needing sign-off rather than a silent config edit.

**Z0 does two jobs**, which is why it is not merely a control: it proves the Stan edit is behaviour-neutral when off, and it IS the ladder's daily rung, fitted under the same Stan file as the other three. Without the second job the ladder would compare three post-edit rungs against a pre-edit production run.

---

## 1g. Pre-run audit (2026-09-03): three defects that would have corrupted the next batch

Found by a full sweep of `06_diagnostics/` and `07_documentation/` before starting the 2026-09-03 batch. All three are fixed; the first two had already caused real damage.

**1. An ungated overwrite erased the fitted stages' verdicts.** Every batch runner rewrote its verdicts file from whatever the CURRENT invocation produced. A batch with desk stages that run during a dry run therefore truncates the file to desk rows the moment anyone re-sources it, and the 2026-09-01 batch's own closing message told the reader to do exactly that. `validation_2026-09-01_verdicts.csv` was left holding 11 desk rows and none of the 7 fitted-stage criteria, while the summary still listed all five fits: the runs happened, only their scored criteria were erased. **Restored 2026-09-03** from `validation-batch-review-2026-09-02.md`, which had quoted them in full.

**2. Append-on-resume duplicated summary rows.** `append_row()` appended whenever the file existed, and the `RESUME` skip path appends too. Since RESUME is the documented recovery for an interrupted multi-hour batch, a resume wrote a second row per completed stage: `osp_validation_summary.csv` carries 28 rows for a 14-job matrix and `patch_validation_2026-08-25_summary.csv` 12 for 6 jobs.

Both are fixed by one change in `run_shore_ar_zi_2026-09-03.R`: **`merge_csv_by()` merges by key** (summary by `stage`, verdicts by `stage` + `criterion`, ladder by `rung`), so re-running a stage replaces its rows and leaves the rest alone. The AR ladder's pre-edit fallback rung is now merged rather than substituted, which had been truncating that file to a single row on every dry run.

**3. `run_patch_validation_2026-08-25.R` shipped `DRY_RUN <- FALSE`.** Sourcing it started real multi-hour fits, while its own inline comment said "fit nothing" and `06_diagnostics/README.md` told the reader to set it TRUE first. It sat that way from 2026-08-25 through 279 passing assertions because the harness checked the three most recent runners by name. **The assertion is now generated over every `run_*.R`**, and a second assertion pins the set of runners that have no dry-run mode at all so it cannot grow unnoticed.

### Documentation corrections

`PIPELINE_STATUS.md` gave **six different answers** for the authoritative run across its own length (Run 6 at `:6`, rung 4 at `:20`, S2 at `:367`, VAL-1-adopted at `:422`, Run 6 again in the repo map), and stated both that `run_config.R` ships `shared_tau = FALSE` and that it ships TRUE. There is now **one box at the top of the file** carrying the authoritative run, and every superseded reference is marked as historical. Also corrected: the LICENSE placeholder claim (done), the `05_output` size (6,837 files / 565 MB, not 1,017 / 82 MB), the helper-module count (37, not 33), the `04_input_files` list (all `.xlsx` since the 2026-07-16 migration), and the `06_diagnostics` description (the validation harness, not just the weather module).

`CLAUDE.md` never mentioned `shared_tau` at all, the single production-changing decision of the last month. It now lists every config key added since it was written, corrects "two toggles force a recompile" to three, points at `bss_sampler_override` as the sanctioned way past the merge-direction trap it already warned about, and describes `06_diagnostics/` as the harness it has become.

Both **method documents** now carry a superseded-numbers banner: they present 67,312 (pooled) and 66,461 (gear) as production, the gear one with no staleness hedge of any kind, and both still instruct an operator to edit `crabbing_holiday_dates`, a config key that no longer exists.

`effort_overdispersion_diagnostic_HOWTO.md` instructed the reader to act on `coverage_50` from `ppc_calibration_*.csv`, which is the non-randomized statistic corrected on 2026-09-01. It now carries a warning and the `ppc_byobs` replacement recipe.

**Not fixed, recorded instead:** three rendered `.html` files are content-stale by up to seven weeks; `PIPELINE_STATUS.md` Section 4 still carries ~14 items marked DONE inline rather than moving them to Section 3; five superseded source notes referenced by the ID legend no longer exist, so the `T1.x` / `GR-9` / `POOL-1` identifiers are untraceable; and `07_documentation/README.md` indexes only three development notes.

---

## 1h. The AR escalation ladder as a production toggle (2026-09-04)

Came out of the FW creel team meeting. Harness **303 assertions**.

### Project context, recorded because it changes how everything above should be read

**WDFW has published no estimate from this pipeline.** `main` is the state of the model before the FW creel meeting and before OSP confirmed they can supply boat counts; this branch is the work incorporating both. **There is no published figure a change has to stay consistent with**, so continuity with an earlier internal run is not on its own a reason to prefer one modelling choice over another. This is now stated at the top of `run_config.R` and `CLAUDE.md`.

**The OSP data** is, per day, the total returning vessels and the fraction that were *crabbing only*, which deliberately EXCLUDES combo trips that also crabbed. It is therefore a **lower bound on crabbing vessels, not the crabbing fraction `f`**, which is exactly what `osp_crab_lower` / `f_lower` were built for. Its purpose is to improve the accuracy of the boat estimate and reduce its uncertainty; the acknowledged limitation is that more **boat interviews** are needed next season, which boat counting cannot fix.

### What the toggle does, and its one important limit

`ar_escalate` already existed. What it lacked, and now has:

- **A per-population and per-sub-season scope.** `TRUE`, `c("shore")`, or `list(shore = "all_gear")`. Every rung is a multi-hour fit, and two components already have known answers (shore pot closure funnelled at daily with 1,165 divergences; the boat diverged on ~100% of iterations), so escalating those from the top burns known-bad fits.
- **Each rung's own estimate.** The log kept sampler diagnostics only; the estimate was discarded when `fit` was overwritten. `ar_escalation_log.csv` now carries `catch_median`, `catch_lo95`, `catch_hi95`, `pi_width`, `pi_width_rel`, `effort_median` and a `selected` flag per rung, and the HTML report tabulates them.
- **`ar_escalate_stop`**: `"first_pass"` (default, and byte-identical to the previous behaviour) stops at the finest rung that passes the gate, which is the production rule. `"all_rungs"` fits every rung and reports the one `ar_escalate_select` names, which is the diagnostic mode for honing in on a resolution in a new season.

> **THE LIMIT, stated because it bounds what the mechanism can do.** The ladder escalates on the CONVERGENCE gate. On the 2026-08-31 production run **every component passes that gate, shore all-gear at daily included**, so the production rule applied to this season changes nothing. That same fit carries p_loo at 35.2% of `n_obs`, 41 Pareto k above 0.7 and `coverage_50` 0.701 against a nominal 0.500. The gate is blind to all of it by construction. Making the ladder react to adequacy is a deliberate change to what the gate tests, not a setting.

### On selecting by the narrowest prediction interval

Tested against the boat 2x2, four cells on one track:

| cell | catch | PI width | rel width | trailer coverage_50 |
|---|---:|---:|---:|---:|
| OFF x monthly | 25,868 | 40,358 | **156.0%** | 0.538 calibrated |
| ON x monthly | 31,008 | 47,671 | 153.7% | 0.523 calibrated |
| OFF x daily | 37,359 | 54,769 | 146.6% | 0.713 broken |
| ON x daily | 42,344 | 60,294 | **142.4%** | 0.744 broken |

On **relative** width the two miscalibrated cells look the most precise, because a latent process that absorbs observation noise reports a tighter interval. On **absolute** width the calibrated cell wins, but catch and width are nearly proportional here (ratio 0.64 in all four cells), so "narrowest absolute interval" is close to "smallest estimate" and would systematically select the lowest harvest number. `"narrowest_pi"` is available and is **not** the default; use precision to break a tie between rungs that are already adequate, never to decide adequacy.

### Batch simplification

The 2026-09-03 batch replaced its three forced rungs (A1/A2/A3) with **one** run of the production toggle: 7 fits instead of 12, one output folder instead of three, and it exercises the mechanism a future season will use rather than an experiment-only `ar_force`. Runtime drops from 23-26 h to 17-20 h.

---

## 1i. The 2026-09-03 batch: two good stages, three defects, and a ladder that did not ladder (reviewed 2026-09-04)

Full review: `07_documentation/development_notes/shore-ar-zi-review-2026-09-04.md`. Harness **355 assertions, 0 failures**.

**The authoritative run does not move.** `20260831/pooled-CPUE-VAL-1-adopted`, port 71,513, still stands.

| stage | model | min | port | shore all-gear | AR | boat all-gear | outcome |
|---|---|---:|---:|---:|---|---:|---|
| P1 | PE only | ~0 | - | - | - | - | PASS, -0.059% on the boat |
| G1 | gear | 26.5 | 70,953 | 20,754 | monthly | 30,760 | PASS x2 |
| Z0 | pooled | 263.9 | 71,450 | 20,898 | daily | 31,008 | PASS |
| A1 | pooled | 883.0 | 71,535 | 20,898 | daily | 31,008 | **VOID** |
| Z1 | pooled | 242.4 | 71,326 | 20,745 | daily | 31,008 | re-scored |

### What worked

**G1 restores the cross-track check.** The gear-track boat all-gear fit had never sampled in production: 554 divergences, gate rejection, PE fallback, and a gear port total of 51,385, a 27% understatement of the number that cross-validates the pooled headline. With the per-fit sampler settings moved from an experiment override into the driver, it now runs at **2 divergences**, `tau_bar` n_eff **19,292** against 49, and a port total of **70,953**, which is **-0.78%** from the pooled production run (the 2026-08-26 shipped pair was 0.73% apart). The two independently parameterized tracks agree on `tau_bar` to four decimals, 2.5962 against 2.5962, which is what makes the shared turnover a property of the data rather than of the pooled parameterization.

**P1 closes the incomplete-trip inconsistency.** The PE and BSS arms of one fused estimator no longer disagree about which interviews count. Largest movement -0.059% on the boat, shore unchanged, which is the control.

**Z0 confirms the ZINB Stan edit is inference-neutral when off.** 11,223 shared parameter rows identical at full precision.

### The four defects

| id | what | cost |
|---|---|---|
| **D1** | `force_res <- if (isTRUE(params$ar_escalate))` in both drivers. `ar_escalate` became SCOPED the same day; `isTRUE()` on `list(shore = "all_gear")` is FALSE, so the resolution never reached the data prep while `.esc_on` twelve lines above resolved correctly. All four rungs ran at `daily`, `P_n = 289`, catch 20,897.819 to three decimals. | **14.7 h of 23.6 h; both A1-A3 verdicts VOID** |
| **D1b** | `selected` was marked by matching the resolution string and taking the last hit, so the log named attempt 4 while attempt 1 was kept. The gear driver never set it at all. | log unreadable exactly when it matters |
| **D2** | `pit_block()` and `calib()` scored the catch stream as plain NB2 on a ZINB fit, excluding `theta_C`. Reported **419.0** expected zeros against 676 observed, z = **+15.4**, i.e. "the feature made the zero bin worse". Under the mixture it is **637.9, z = +2.0**, against the NB2 baseline's 605.5 / +3.8: the bin halves. Also makes Z1's catch `pit_mean`, `pit_sd`, `coverage_50` and `flag_pit_bias` unreadable. `elpd` unaffected; Stan's `log_lik` always carried the mixture. | the feature was nearly rejected on it |
| **D3** | the elpd criterion compared +14.8 nats against `se_elpd_loo` (~46), the SE of ONE model's total, which is dominated by across-observation variation common to both models and cancels in the difference. The paired-difference SE is **5.5**, so the gain is **2.69 SE** and clears the stated 2 SE bar. | "NOT WORTH IT" was wrong |
| **D4** | `run_parameters.txt` is written by `str(params)`, which caps a list at `list.len = 99`; `run_config` carries 120 keys. **Every dump in the repo before 2026-09-04 is truncated.** `config_delta()` and `annotate_decoupled_run()` both read it, so it yields false positives (the Z0 "9 UNEXPECTED keys" caveat, on a correct config) and false negatives (a key past entry 99 in both dumps is never compared). `opener_covariate_mode` fell past the cut from 2026-09-01 on, so post-hoc annotation of those folders would leave `B_open` unflagged. | one false caveat; one latent missing-flag |

**Five of eleven verdicts were wrong.** The committed `05_output/shore_ar_zi_2026-09-03_verdicts.csv` now carries `verdict_2026_09_04` and `correction_2026_09_04` beside the originals so the supersession is auditable rather than rewritten. An injected `[A1] PASS fake fitted criterion` test row that was never removed is deleted, and `run_shore_ar_zi_2026-09-03.R` shipped `DRY_RUN <- FALSE` and is reset.

### The ZINB decision, re-scored

| stream | elpd diff | paired SE | ratio | zeros | positives | zero bin z |
|---|---:|---:|---:|---|---|---|
| shore all-gear | +14.8 | 5.5 | **2.69** | n=676, +25.2 (11.6 SE) | n=973, **-10.5 (-2.1 SE)** | +3.8 -> **+2.0** |
| shore pot closure | +9.1 | 4.6 | 1.96 | n=146, +19.0 (8.9 SE) | n=481, **-10.0 (-2.6 SE)** | +2.9 -> **+0.8** |

Pareto k>0.7 on the shore all-gear fit fell **41 -> 26**. `theta_C` 0.179 on all-gear, 0.107 on the independent pot-closure replicate.

**CORRECTED 2026-09-05: the "positives get worse" reading above was misleading.** Split by count size on shore all-gear, the elpd change is y=0 **+25.2**, y=1 **-42.0 (16.7 SE)**, y=2 -6.6, y=3-4 **+14.3**, y=5-8 **+18.1**, y=9-16 **+5.6**, y=17+ +0.1. The 3+ bins carry **87% of the shore catch** and every one of them improves; the loss is concentrated at y=1. Mechanism: `r_C` doubles (0.95 -> 1.83) because once `theta` absorbs structural zeros the NB2 no longer needs extreme overdispersion to reach zero, so it tightens, fitting the harvest-carrying counts better and the almost-zero count of 1 worse. The pot-closure replicate shows the same shape (y=1 -20.4, `r_C` 1.68 -> 2.86). **The ZINB has moved the misfit from the 0 bin to the 1 bin, not removed it.** Excess mass at both 0 and 1 relative to NB2 is the signature of a two-regime process (unsuccessful trips yielding 0-1 crab against successful ones), which a hurdle or two-component NB mixture fits and a ZINB cannot.

**`estimate_catch_zi` stays FALSE**, now for one reason rather than two: the decision needs a rendered count-bin table (0 AND 1, not the zero bin alone, which passed on the prototype while the misfit migrated) from a run under the corrected PPC. Stage Z2's criterion is pre-set in the runner: both bins |z| under 2.5 or it fails as a whole, and if the one bin fails the next model is a hurdle or two-component mixture, not a retuned `theta` prior. The offline `E[theta]*E[p0]` approximation is checked against the rendered value at the same time. `ppc_byobs_*.csv` now carries `p_one` beside `p_zero`, and `loo_elpd_paired()` returns the by-count table, so this read is reproducible from any pair of runs.

### The structural change: posterior draws are now persisted

Three defects in five weeks (the quantile-interval coverage defect 2026-08-31, the `tau_bar` row-key defect 2026-09-01, D2 above) were all in code that runs AFTER sampling, on a fit that was never wrong, and each cost a multi-hour re-fit to correct a file. Nothing persisted a stanfit or any subset of its draws. `save_ppc_draws = TRUE` now writes `ppc_draws_<fit>.rds` per fit, tens of MB, carrying the draw objects every R-side predictive statistic reads plus `stan_data`. The next such fix is a recomputation. See `03_R_functions/save_bss_ppc_draws.R` for what it does and does not cover.

### Why this keeps happening

Three of the four defects share a shape: **a change was validated by the mechanism it was about to break.** The ladder was tested by a harness that checked `bss_ar_ladder()`'s output and never that the resolution reached the data prep. The ZINB feature was validated by a zero-bin diagnostic written before the mixture existed. The config comparison was validated by reading the file it was silently truncating. The new assertions target that shape: they read the driver SOURCE and assert the property that failed, not the helper's return value.

### Next: `06_diagnostics/run_ladder_zinb_2026-09-04.R`, about 6-7 h against the 23.6 h it corrects

**Pre-run audit 2026-09-05** (review doc section 10): traced by execution rather than by reading. `prep_bss_crab_pooled()` returns `P_n` 289/44/21/10 for NULL/weekly/biweekly/monthly under the L1 config; the draw-persistence and theta-extraction code reproduce the live PIT exactly from a saved `.rds` on a real ZINB stanfit; nothing the batch injects is overridden by the driver merge. One defect found and fixed that a dry run cannot reach: the runner set `run_tag` where the driver does not read it, so L1 and Z2 would have shared the default output folder and **Z2 would have overwritten the L1 ladder**. Now asserted over every runner (harness 355).

- **D0**, desk, free, runs in a dry run: re-derives every corrected number above from committed files.
- **L1**, ~2-3 h: the shore all-gear ladder for real, via the production `ar_escalate` toggle. `LADDER_INCLUDE_DAILY` defaults FALSE and reuses Z0's daily rung (config differs only in `ar_escalate` keys, none of which touch model, data or seed), saving 3.7 h; set TRUE if the ladder will be shown outside the project. Decision rule as agreed: escalate on the convergence gate, report the finest rung that PASSES, narrowest relative PI only as a tie-break, and check `cov50` and Pareto k on any rung the tie-break selects.
- **Z2**, ~4 h, switchable: re-render the ZINB shore fit under the corrected PPC, so adoption rests on rendered output; also checks the rendered zero bin against D0's approximation and confirms the boat negative control reproduces bit-for-bit.

---

## 1j. The 2026-09-04 ladder / ZINB batch: the ladder ran, and it says daily is overfitted (reviewed 2026-09-06)

Full review: `07_documentation/development_notes/ladder-zinb-review-2026-09-06.md`. Harness **372 assertions, 0 failures**. 9.1 h of fitting (L1 5.4 h, Z2 3.7 h).

### The ladder, for real this time

| rung | P_n | divergences | catch | 95% CI | PI rel | gate |
|---|---:|---:|---:|---|---:|---|
| daily (from Z0) | 289 | 333 | 20,898 | [18,187, 24,402] | 0.2974 | PASS |
| **weekly (reported)** | 44 | 313 | **21,547** | [19,076, 24,528] | 0.2530 | PASS |
| biweekly | 21 | 212 | 21,383 | [19,143, 24,041] | 0.2291 | PASS |
| monthly | 10 | 152 | 20,771 | [18,628, 23,217] | 0.2209 | PASS |

Four distinct resolutions, four distinct `P_n`, four distinct estimates: the mechanism works.

**The AR resolution is NOT a large lever on the shore point estimate.** The four rungs span 777 crab, **3.7%** of the component and about 1% of the port total, against 44% for the Stage 5 boat 2x2. The worry that motivated this line of work is not borne out at the level of the estimate.

**The gate separated nothing.** All four rungs pass it, so the production rule ("report the finest rung that PASSES") returns the finest rung RUN whatever its adequacy. That is the limit stated in Section 1h, now demonstrated.

**The adequacy statistics separate them decisively, and they say daily is overfitted.**

| statistic | daily | weekly |
|---|---:|---:|
| `p_loo` / `n_obs` | **0.352** | **0.096** |
| `p_loo`, gear stream | 109.4 | 29.7 |
| Pareto k > 0.7 | **41** of 311 | **1** |
| gear `coverage_50` | 0.701 (**+7.1 SD**) | 0.559 (**+2.1 SD**) |
| catch `coverage_50` | 0.457 (-3.5 SD) | 0.441 (-4.8 SD) |
| `flag_miscalibrated` | **TRUE** | **FALSE** |

`elpd` favours daily (gear -85.8 nats at a paired SE of 11.7; catch -7.4 at 3.6) and **should not be believed**: the gear gain is 85.8 nats for 79.7 extra effective parameters, i.e. 1.08 nats each, which is this project's own written definition of buying noise; and **41 of the daily fit's 311 gear observations have Pareto k above 0.7**, so PSIS-LOO has failed for 13% of the stream and the daily fit's own elpd is not a reliable estimate of anything. Weekly has 1.

**Port total at weekly: 72,122 [53,228, 101,553]**, +0.94% on Z0's 71,450 in the same configuration and +0.85% on the current authoritative 71,513, with the lower bound rising.

**What it does NOT establish.** Only the REPORTED rung got adequacy, because only the kept fit reaches `write_bss_diagnostics()`. Biweekly and monthly cost 174 minutes and left no `p_loo`, no Pareto count, no coverage. The gear track fits this component at monthly and gets `coverage_50` 0.035 (-16.4 SD), so the turnover is between weekly and monthly and is still unlocated. Fixed for future runs: `03_R_functions/bss_rung_adequacy.R`, `ar_rung_adequacy = TRUE`, writes them into every rung's row.

**Runtime.** The rungs cost 95.6, 91.1 and 83.3 minutes. **Coarsening the AR barely reduces the cost**, so budget a ladder at one full fit per rung; my 2-3 h estimate for L1 was wrong and it took 5.4 h.

### The ZINB, re-rendered: two of my own readings were wrong, one in the ZINB's favour

Z2 reproduces Z1 exactly (`theta_C` 0.1780, `C_expected_sum` 20,872.6), so the re-render control passes.

- **The offline approximation was accurate.** Zero bin renders at 676 vs **638.5**, z = **+2.0**; the `E[theta]*E[p0]` approximation gave 637.9 / +2.0. Trustworthy in future.
- **"The misfit migrated to y = 1" was WRONG.** With `p_one` now written: the NB2 (L1, weekly) over-predicts ones at 236 vs **335.4, z = -6.1**; the ZINB (Z2, daily) at 236 vs **285.0, z = -3.2**. The ZINB roughly **halves both** the zero-bin and the one-bin misfit. Caveat: the two differ in AR as well as likelihood, because Z0 predates the `p_one` column.
- **The y=1 elpd loss is arithmetic, not shape.** -42.0 nats over 236 observations against log(1 - 0.178) x 236 = **-46.3**: it is the mass tax a mixture levies on every non-zero observation. The `r_C` tightening (0.947 -> 1.828) repays it at y >= 3, where 87% of the catch is.
- **Catch calibration, readable for the first time:** ZINB `coverage_50` **0.4645 (-2.9 SD)** against NB2's 0.4572 (-3.5 SD). The 2026-09-03 batch reported 0.4051 (-7.7 SD) and raised `flag_pit_bias` on it; that was the wrong-likelihood artefact.
- **It does not fix the AR.** Z2 still shows `p_loo_frac` 0.3437, 26 bad k, `flag_miscalibrated` TRUE, because Z2 runs at daily. **The two changes have never been run together.**

### Defects in the run

- **L1's HTML report was destroyed.** `rmarkdown::render()` writes beside the `.Rmd`; the 2026-09-03 runner moved it afterwards and this one did not, so L1 rendered there and **Z2 overwrote it**. The rendered ladder table, the artefact the FW creel discussion asked to preserve, is gone and needs a refit to regenerate; `ar_escalation_log.csv` kept the numbers, so the loss is cosmetic. A 6,062-line render artefact was also committed inside `01_BSS_models/`. Fixed: the runner moves it and warns on failure, `01_BSS_models/*.html` is gitignored, the stray file is deleted, and the harness asserts both over every runner.
- **`fmt()` was scalar-only** and `verdict_L1` hands it whole ladder columns. On **R >= 4.3 that is an error**, so the batch would have aborted in the verdict block after 9.1 h of fitting; it survived only because the run used an older R. Vectorised.
- **`DRY_RUN <- FALSE` committed again** (third time: 2026-08-25, 2026-09-03, now). Reset.
- **`ppc_draws_*.rds` are gitignored**, so the 4 files per stage stayed on the machine that ran the fit. Intended, but as of today **no recompute path reads them back**, so the feature is still write-only. The verdict now says so instead of reporting a misleading zero.

### Where the two decisions stand, and what settles them

**Shore all-gear AR**: the evidence supports moving off daily. What is missing is the other half of the bracket, so "weekly" today means "the finest rung that is clearly not overfitted", not "the best rung".

**ZINB**: every objection I raised has been answered by rendered output. It ships FALSE for one reason only, that it has never been run at the resolution the ladder points to, and adopting a likelihood change and a resolution change together would leave neither attributable.

**Next: one run settles both.** `C1`, shore all-gear at weekly WITH `estimate_catch_zi = TRUE`, boat untouched as the control, about 4.5 h; read against L1 to isolate the likelihood and against Z2 to isolate the resolution. Share the run with `C2`, the ladder again with `ar_rung_adequacy = TRUE` so biweekly and monthly finally report adequacy. If C1 holds, one clean production render supersedes `20260831/pooled-CPUE-VAL-1-adopted` carrying both changes.

---

## 1k. The C1 / C2 run (2026-09-04): the 2x2 is complete and there is a candidate configuration (reviewed 2026-09-07)

Full review: `07_documentation/development_notes/candidate-config-review-2026-09-07.md`. Harness **378 assertions, 0 failures**. 9.3 h of fitting (C1 3.8 h, C2 5.5 h).

### The 2x2

| configuration | shore all-gear | port total | `p_loo_frac` | Pareto k>0.7 | `cov50` dev | miscalibrated |
|---|---:|---:|---:|---:|---:|---|
| daily + NB2 (production) | 20,898 | 71,450 | 0.3519 | **41** | 0.2042 | **TRUE** |
| weekly + NB2 (L1) | 21,547 | 72,122 | 0.0956 | 1 | 0.0691 | FALSE |
| daily + ZINB (Z2) | 20,745 | 71,287 | 0.3437 | 26 | 0.2010 | **TRUE** |
| **weekly + ZINB (C1)** | **21,489** | **72,032** | **0.0949** | **0** | **0.0691** | FALSE |

**The effects are additive**: the interaction on the shore component is **+95 crab**, 0.4%. The AR resolution does all the work on `p_loo` and coverage; the ZINB does the rest on the Pareto diagnostics. **C1 is the only fit in this project with zero Pareto k above 0.7.** Port total **72,032 [53,044, 101,212]**, +0.73% on the current authoritative 71,513, with the lower bound rising.

### A correction carried since 2026-08-31

I have repeatedly framed the shore all-gear AR as bracketed by two bad ends, daily at `coverage_50` 0.701 (+7.1 SD) and monthly at 0.035 (-16.4 SD). **The second figure is from the GEAR-RESOLVED model**, a different likelihood with per-gear CPUE, and says nothing about where the pooled model's AR belongs. The pooled model's own bracket, measured by C2:

| rung | `cov50` gear | `cov50` catch | total `p_loo` | bad k | catch | PI rel |
|---|---:|---:|---:|---:|---:|---:|
| daily | 0.7010 (**+7.1 SD**) | 0.4572 (-3.5 SD) | 153.6 | 41 | 20,898 | 0.2974 |
| weekly | 0.5659 (+2.3 SD) | 0.4433 (-4.6 SD) | 54.8 | 1 | 21,547 | 0.2530 |
| biweekly | 0.5659 (+2.3 SD) | 0.4500 (-4.1 SD) | 34.2 | 0 | 21,383 | 0.2291 |
| monthly | 0.5788 (+2.8 SD) | 0.4506 (-4.0 SD) | 22.2 | 0 | 20,771 | 0.2209 |

**Only daily is an outlier.** The three coarse rungs sit within 0.5 sampling SD of each other and span 3.6% in catch, so **the choice among them is not made by the evidence**; it is made by the agreed rule (the finest rung that passes), which selects weekly and is defensible here because weekly is adequate on every statistic rather than merely passing a gate everything passes. Separately, the catch stream is under-covered at EVERY resolution (-3.5 to -4.6 SD): that is a likelihood problem, not a resolution problem, and it is what the ZINB addresses.

### The ZINB, compared like for like at last

C1 against L1 is the comparison every earlier run lacked: same resolution, same data, one likelihood difference.

| stream | bin | NB2 | ZINB |
|---|---|---|---|
| shore all-gear | zero | z = **+3.7** | z = **+2.0** |
| shore all-gear | one | z = **-6.1** | z = **-3.3** |
| shore pot closure | zero | z = +2.9 | **+0.7** |
| shore pot closure | one | z = -3.5 | **-1.2** |

Both bins halve on all-gear and close on the pot closure. `elpd` +11.6 nats at **2.30 paired SE** (it was +14.8 / 2.69 at daily): the gain shrinks when the effort process is no longer overfitted, which is the right direction, and still clears 2 SE. **The all-gear one bin does not close** (z = -3.3): the data is more bimodal than a zero-inflated NB can be, and a hurdle or two-component mixture is the shape that fits it, with the remaining gain bounded at about 6% of the catch.

### Four defects, all in code from the previous patch

- **`bss_rung_adequacy()` reported the wrong `p_loo_frac`**: the sum over streams / summed `n_obs`, where `bss_model_adequacy.R` reports the WORST stream's fraction. The C2 ladder printed 2.8% for a rung whose `model_adequacy.csv` figure is 9.6%. Fixed and verified against a real fit.
- **It shifted the global RNG.** `bss_ppc_calibration()` reseeds, and both it and `save_run_diagnostics()` subsample draws, so calling the helper inside the fitting loop moved every diagnostic computed afterwards: C2's kept fit is **bit-identical to L1's across 10,251 parameter rows** and still reported gear `coverage_50` 0.5659 against 0.5595. A diagnostic-only toggle moved a committed diagnostic. Fixed with a save/restore.
- **C2's ladder rows were filed under the stage label `L1`** and `merge_csv_by()` de-duplicated only against the OLD file, so the verdicts carried two rows per criterion. Both fixed; the committed file is relabelled and annotated.
- **`DRY_RUN <- FALSE` committed** for the fourth time.

The C2 per-rung numbers carry the superseded `p_loo` definition and the RNG shift and cannot be recomputed from what was written, so they are annotated rather than silently corrected. Read those rungs against each other, not against `model_adequacy.csv`.

### Recommendation: adopt C1

`ar_max_resolution$pooled$shore$all_gear = "weekly"` and `estimate_catch_zi = TRUE` scoped to shore. C1 used `ar_force`, an experiment lever; production expresses the same thing through the cap, which is a one-line difference needing one confirming render that doubles as the new authoritative run.

**Next: P1, the adoption render** (about 4 h), gated on reproducing C1's four fits bit-identically, since routing the same resolution through the cap rather than the override should change nothing. **P2**, the gear-track cross-check, can share the run: the two-track agreement is the strongest internal check the project has and has not been re-measured since the shore component moved.

---

## 1l. Adoption (2026-09-07): run_config now ships the candidate, pending one confirming render

`run_config.R` has been changed:

```r
ar_max_resolution$pooled$shore$all_gear  : "daily"  ->  "weekly"
estimate_catch_zi                        :  FALSE   ->  TRUE     (scoped to shore)
```

Nothing else moved: the shore pot-closure cap stays biweekly, the boat stays monthly, and **the entire gear-resolved map is untouched**, so the cross-check measures the pooled change alone.

**Why the change ships before the render.** The authoritative run should come out of `run_config` as shipped, not out of a batch override, or its provenance is a footnote. The cost is that a failed gate leaves `run_config` carrying an unvalidated configuration, so **if A1 fails, revert both lines together**; the routing would be at fault, not the resolution.

**The routing proof, run before spending four hours on it.** Stage C1 reached weekly through `ar_force`, which takes the `fixed` branch of `bss_select_ar_resolution()` and skips the cap; production reaches it through the adaptive branch, which selects daily from the effort density and is then coarsened. **Verified: both build byte-identical Stan data, weekly with `P_n` 44, no entry differing.** The desk stage of the runner repeats this check every time and refuses to start the render if it fails.

**One asymmetry this creates, stated because a reviewer will find it.** `crab_bss_gear_resolved.stan` has no `theta_C` and `prep_bss_crab_gear.R` never emits `zi_catch`, so the gear track silently ignores `estimate_catch_zi` and fits plain NB2. **The two tracks now differ in the shore CATCH LIKELIHOOD as well as the AR resolution**, for the first time. The ZINB is worth about -0.3% on the pooled shore component, so it explains a sliver of any cross-track gap and none of a large one. Porting the block to the gear model would restore symmetry; it is a Stan edit plus a recompile and is not needed to read the cross-check.

### The run: `06_diagnostics/run_adoption_2026-09-07.R`, about 4.5 h

- **A1** (about 4 h), production as `run_config` ships it, **no config delta at all**. Gate: the four fits must be **bit-identical to `20260904/pooled-CPUE-LZ-C1-weekly-zi`**. A PASS means the adoption is a routing change and nothing else, which is what makes the folder callable authoritative. Judge the gate on the per-fit summaries, never on the port total: `rstan::extract(permuted = TRUE)` permutes draws, so a port total moves about 0.2% between bit-identical fits.
- **A2** (about 0.5 h), the gear-resolved cross-check, not re-measured since the shore component moved. **Expect the gap to NARROW**: the gear track fits shore all-gear at monthly, so the two tracks are now closer in resolution than they were (last measured pair -0.78%). If it widens, the pooled change did something the gear track disagrees with, and that is worth more attention than the port total.

Harness **391 assertions, 0 failures**, including a section that asserts `run_config` ships the adopted configuration with no experiment lever active and the gear map untouched.

---

## 1m. The adoption render (2026-09-07): the gate passed and the authoritative run moves (reviewed 2026-09-08)

Full review: `07_documentation/development_notes/adoption-review-2026-09-08.md`. Harness **452 assertions, 0 failures**. 4.1 h of fitting (A1 220.4 min, A2 26.4 min).

### The gate

**A1's four fits are identical to stage C1's across 10,253 shared parameter rows, at full precision.** C1 reached the weekly shore AR through `ar_force`, which takes the `fixed` branch of `bss_select_ar_resolution()` and skips the cap; A1 reaches it through `ar_max_resolution`, which selects daily from the effort density and then coarsens. Different code paths, same posterior. The desk stage proved the two build byte-identical Stan data before the render started; the render proved the fits match. **The adoption is a routing change and nothing else**, which is what makes the folder callable authoritative.

### The new authoritative run

**Port total 72,027 [53,018, 101,364]** (itself superseded 2026-09-12 by the ladder's R4, 94,376; see the box), superseding 71,513, **+0.72%**, essentially all of it the shore effort process no longer being over-imputed at daily. C1 read 72,032 on the same fits; the 5-crab difference is the port assembly resampling permuted draws, which is why the gate is judged on per-fit summaries and never on the total.

| fit | `p_loo_frac` | Pareto k>0.7 | `cov50` dev | miscalibrated |
|---|---:|---:|---:|---|
| shore all-gear | 0.0949 | **0** | 0.0659 | FALSE |
| shore pot closure | 0.1129 | 1 | 0.0769 | FALSE |
| private boat all-gear | 0.0801 | 1 | 0.0954 | FALSE |
| private boat pot closure | 0.2120 | 1 | 0.1136 | FALSE |

**No fit carries the miscalibration flag**, the first run for which that is true.

### The cross-check widened, and the reason is the interesting part

The gear track reads 71,026, **-1.39%** apart, against -0.78% for the last pair. I predicted it would narrow. It did not, and the number that explains why is:

| comparison | pooled | gear | gap |
|---|---:|---:|---:|
| shore all-gear, as configured (weekly vs monthly) | 21,489 | 20,754 | -3.42% |
| shore all-gear, **both at monthly** | 20,771 | 20,754 | **-0.08%** |

**At the same resolution the two tracks agree to 0.08%**, 17 crab, between a pooled CPUE structure and one carrying per-gear CPUE at `G = 5`. That is the strongest agreement this cross-check has produced. The -1.39% is therefore a resolution difference, not a disagreement about the fishery: the gear track's `ar_max_resolution` was deliberately left untouched so the cross-check would measure the pooled change alone. `tau_bar` agrees to 0.05% (2.6128 vs 2.6114).

**Two controls held.** The gear track is **bit-identical to its own 2026-09-01 baseline** on all four components, confirming `estimate_catch_zi` is genuinely inert there. That inertness is itself a fact worth stating: `crab_bss_gear_resolved.stan` has no `theta_C`, so the two tracks now differ in the shore CATCH LIKELIHOOD as well as the resolution, worth about -0.3%.

**The question it raises**: should the gear track's shore cap move to weekly too? Probably, but not on this evidence. Its shore fits carry a thinner per-gear likelihood and may genuinely need a coarser AR; that needs its own ladder, and the machinery now records per-rung adequacy.

### One defect, the third of its shape

**`verdict_A2` crashed after 4.1 h of fitting and the run's verdicts were never written.** The pooled driver labels the port-total row `"Expected_Catch"`, the gear driver `"Catch"`; reading one label against the other file gives `numeric(0)`, and `is.finite(numeric(0)) && ...` is an **error** in R, not FALSE. Nothing was lost from the run, and re-scoring from disk takes seconds, but this is the third time verdict code has aborted a batch after all the expensive work succeeded. Fixed three ways: match the row by pattern; **wrap every verdict block so a defect records itself and the verdicts are written regardless**; assert both in the harness, including that the two tracks really do still use different labels so the pattern match cannot be simplified away.

---

## 1n. Season-portability audit (2026-09-09): the program goal restated, and what it changed

**The program goal, stated by Matt.** This branch supersedes `main`. **2024-25 was the development TEST season, not the target**: the model must run on any user-selected window (full season, part-season, multi-season span) and, on a naive new window, hand the user the information to tune the fits: the AR ladder finds the minimum resolution at which each fit passes the gate, and per-rung adequacy says what a finer or coarser choice gains and loses. Optimizing settings for 2024-25 was never the point; the settings are the worked example, the WORKFLOW is the product. Recorded in `CLAUDE.md`, at the top of `run_config.R`, and as the framing of the new `07_documentation/NEW_SEASON_GUIDE.md`.

**A fresh audit of the branch against that goal found the architecture largely ready** (the sub-season builder was generalized in July, the calendar indices are span-safe sequential factors, absent OSP/I-E degrade gracefully, the ladder + per-rung adequacy + `ar_force`/cap pair are exactly the tuning loop the goal describes) **and four concrete gaps, all fixed:**

1. **A window that missed the pot closure could not run at all.** Clamping the closure into the window inverted it and hit a stop(): a summer-only window and a winter window both ERRORED instead of yielding a single all-gear sub-season. Fixed in `build_subseasons.R` (a non-intersecting closure now yields one all-gear sub-season; a window inside the closure already worked) and covered by harness tests.
2. **Season selection was dual-keyed with no consistency check.** The window (`est_date_start/end`) and the data filter (`season_filter`) had to agree and nothing verified it: a stale `season_filter` produced empty readers and a wall of PE fallbacks with the cause nowhere on screen. New `validate_season_window.R` runs inside `fetch_crab_data()` for every driver and batch: it prints, per season, the date range and row counts in and out of the window, STOPS when the window captures zero effort and zero interviews, and WARNS on an unmatched season label or a census window wholly outside the estimation window.
3. **Multi-season spans could not be expressed.** `season_filter` now accepts a character vector (readers filter `%in%`; the holidays reader requires rows for EVERY listed season, loudly). The remaining structural limit is documented rather than hidden: ONE closure window per run (`pot_closure_start/end` are scalars), so a span containing two closures runs per season and combines outside the model, tracked as CHANGE_REGISTER **D8**. *(Superseded 2026-09-10: the closures table exists, `pot_closures`; see Section 1o / A14.)*
4. **The 2024-25-derived settings were not distinguishable from architecture.** Every such key in `run_config.R` is now tagged `SEASON-DERIVED` (turnover prior centers, `kappa_OSP` center, the ZI Beta(1,9), `shared_tau_min_obs`, `crab_fraction_set`, the whole `ar_max_resolution` map), the file opens with a NEW SEASON CHECKLIST, and `ar_force` is reframed from "experiment lever" to what the workflow makes it: the per-fit resolution PIN used while deciding (the only way to pin FINER than the selector, since caps only coarsen), migrated into `ar_max_resolution` once settled, exactly as the 2026-09-07 adoption did.

**The deliverable is `07_documentation/NEW_SEASON_GUIDE.md`**: inputs a season needs, the config checklist, the naive-run ladder settings and their honest cost (one full fit per rung; coarsening barely reduces per-fit cost), the decision files in reading order, how to choose among rungs that all pass (the 2024-25 ladder as the labeled worked example: gate silent, adequacy decisive, interval width only breaking ties between already-adequate rungs), pinning via `ar_force` then caps, part-season behavior, multi-season limits, and a failure-mode table where every stop and fallback is explained in the user's terms.

Harness at this section: build_subseasons window cases, vector season_filter, validator stop/pass, guide presence and links.

---

## 1o. Multi-season spans built; the 2023-25 run staged and blocked on data (2026-09-10)

**The request:** run 2023-24 + 2024-25 together, with (1) a pot-closure window fit per season, (2) the whole span in the monthly figures and fit diagnostics, (3) season-level summaries in the HTML.

**What was built (CHANGE_REGISTER A14; closes the D8 architecture gap):**

1. **`pot_closures`** in `run_config.R`: a list of `list(season =, start =, end =)`, one per closure, outranking the scalar pair. `build_subseasons()` emits, per closure, `ring_net_only_<season>` (biweekly, pots excluded) and the following all-gear gap as `all_gear_<season>`; a window starting before the first closure gets `all_gear_pre_<season1>`. Overlaps, missing season labels, and end-before-start stop loudly. Zero or one in-window closure falls through to the scalar path BYTE-FOR-BYTE, so every single-season name, fit label, and filename is unchanged.
2. **`census_windows`**: a NAMED list season -> `c(start, end)`. `estimate_comm_charter()` recurses per window with the scalar keys substituted (the expansion logic exists once), sums the totals, and keeps `by_season` for the report.
3. **Season totals in the report**: a `season-totals` chunk in the pooled Rmd sums catch/effort draws per sub-season `season` tag (PE substituted where the gate failed, census from `by_season`), writes `season_totals.csv` on EVERY run, and renders the table when the span has >1 season, captioned (corrected 2026-09-08, review item 7) that every sub-season is an INDEPENDENT fit, so a season total from a span equals a standalone run up to the pooled I/E day-length regression and the config priors; the first caption wrongly claimed a shared process.
4. **Validator sharpened**: per-season effort/interview counts split; a season with interviews but ZERO effort counts warns in the exact terms that matter ("its effort process would be pure imputation"); warnings fire even under `quiet = TRUE` (found by writing the test: the warning sat inside the print guard).

Requirement (2) needed nothing: month labels are year-qualified (`%Y-%m`) and calendar indices are span-safe sequential factors since July.

**The staged config** (`run_config.R` live, with a single-season rollback block): window 2023-09-16 to 2025-09-15, `season_filter = c("2023-24","2024-25")`, both closures (Sep 16 to Nov 30 each), both census windows, `run_tag = "two-season-2023-25"`.

**Why the run cannot start yet (data, not code; details in D8):** `effort_combined.xlsx` has NO 2023-24 rows (13,629 interviews vs zero effort counts: pure imputation), `wes_commercial_tally.xlsx` and `crabbing_holidays.xlsx` have NO 2023-24 rows (the holiday reader STOPS by design), WBL/OSP boat counts exist only for 2024-25 (survivable; weaker boat stream). Two flagged assumptions: the 2023-24 `census_windows` dates MIRROR 2024-25 (not from records; confirm), and `pot_open_date` remains a single scalar feeding `estimate_L_effective()`, exact only for the first season (A14 known approximation).

Harness section 46 (18 assertions, total 485): the exact 2023-25 four-sub-season shape, legacy fall-through including the single-closure list form, out-of-window drop, overlap and label stops, the census guard, the validator warning under quiet, and staged-config self-consistency (window spans closures; `season_filter` matches closure seasons; `census_windows` keys match).

---

## 1p. The 2026-09-06 review, answered: eight items, nine patches, a ladder to run (2026-09-08)

**The review** (`claude/branch-review-2026-09-06.md` in the project) raised eight items; Matt's decisions and the patches (`0001`-`0009`, applied in order on `bf4be01`):

| # | Item | Decision | Patch | Ships |
|---|---|---|---|---|
| 1 | Crabbing fraction f | use the sampler contacts (every boat approached; non-crabbing boats carry `crabbers = 0`) to update f dynamically, with OSP's crabbing-only counts when they arrive | 0006 (contacts, month strata), 0007 (dynamic f in both Stan models) | ON |
| 2 | Shore turnover | a `time` column was added to `ingress_egress.xlsx` | 0005 (diel profile, derived prior, per-population shared-tau keys) | OFF (`tau_shore_prior_mu = 1.7`; R4 is the case) |
| 3 | `tau_bar` prior | re-centre on the 61-day overlap calibration | 0001 | ON (`"calibration"`, 3.03 on 2024-25) |
| 4 | Census | "should not have any error" | 0008: exact on the 47 tally days, imputed on 23; SE 426 reported, `census_uncertainty = "none"` | reported, not carried |
| 5 | Boat PE turnover | update it | 0001 (the PE reads the resolved prior) | ON |
| 6 | Gear-split Dirichlet on crab counts | build on interviews | 0003 (trip-level bootstrap, per sub-season) | ON |
| 7 | Multi-season caption | correct it | 0004 | ON |
| 8 | Hours filter | remove it if irrelevant under the deployment unit | 0002 (unit-aware filters) | ON |

**Where the review pushed back, and what the data said.** (a) The census is exact on the days it was taken, but the roster left 23 of 70 days unsampled and the day-type means fill them: 3,869 of 11,753 crab, imputation SE 426 (3.6%). The number does not move; the split and the SE are now on every run (`census_daily.csv`, `census_variance.csv`) and the option to carry the SE exists. The better fix is a data request (D10). (b) The shore counts are taken at 9:00-14:00, when presence is about 0.7 of the daily peak, so the peak-based turnover 1.7 under-expands the shore by ~1.47x; the derived value is 2.48 on 40 I/E days (weekday 2.47, weekend 2.48). It is the largest single mover in the series and ships OFF until R4 shows it. (c) The contacts happen during sampler shifts; if finfish boats return later than crab boats the shift-time share overstates f in summer. Carried as a config caveat, with OSP's all-day column as the check (D11).

**The dynamic f, checked on real data before shipping** (boat all-gear 2024-25, 4 chains x 1200, `adapt_delta` 0.95): f tracks the monthly contact shares (Dec 0.947 [0.87, 0.99] on 46 of 47; Feb 0.90 [0.73, 0.98] on 5 contacts, borrowed from its neighbours), `sigma_f` 0.66 [0.34, 1.38], `cfi_kappa` prior-dominated as expected (most contact days hold 1-4 boats), every non-f parameter within Monte Carlo error of the legacy Phase A control (|z| < 1.3, the factorization proof), component +17.5% against Phase A (36.5k -> 42.9k), whose kappa = 20 prior held December at 0.78 and sent Feb/Mar/Jun back to 0.30. Divergences (62 vs 44 of 2,400 at `adapt_delta` 0.95) sit in the pooled boat model's known `sigma_mu_E` funnel in both fits, not in the f block.

**Two corrections to the plan as first written.** The plan's R3 criterion said "non-f parameters bit-identical to R2"; that is impossible once the parameter vector changes (HMC moves every coordinate jointly), and the correct test, agreement in distribution within Monte Carlo error, is now `fit_agreement()` (B13). The plan's census figure (~4,800 imputed of 11,821) was from the earlier run's means; the measured split is 3,869 of 11,753.

**Not yet run:** `06_diagnostics/run_improvements_2026-09-08.R` (B14, D12), about 22 h: R0 desk (ran clean here), R1 filters, R2 calibration prior, R3a monthly f (legacy), R3 dynamic f (the shipped configuration), R4 derived shore turnover, R5 gear cross-check. Harness at 628 assertions.

---

## 1q. The four follow-up decisions, answered in code (2026-09-09)

Matt's answers to the plan's four questions, and what each became (patches applied on `aa68bc6`):

**Item 4, the 23 unsampled census days.** There are no charters or commercial vessels on unsampled days: the samplers are scheduled on the days those operations are confirmed. So the day-type expansion was filling non-operating days with an operating-day mean. `census_expansion = "none"` (shipped) makes the census the exact sum over the tally days: **7,884** crab on 2024-25 (the expansion gave 11,753; the 2026-09-04 baseline carried 11,821). The per-vessel-mean SE that remains (131, 1.7%) is reported and stays out of the port interval. The expansion is kept as `"day_type"` for reproduction and for a season whose roster did not track operations (A22; D10 closed).

**Item 2, the derived shore turnover.** Adopted now: the count hours were not in the I/E workbook when 1.7 was derived, so the old value was the wrong quantity, not a competing estimate. `tau_shore_prior_mu` and `_sigma` ship `"derived"` (2.48 on 40 I/E days, weekday 2.47 / weekend 2.48, log-SE 0.044 floored at 0.10). Rung R4 of the ladder still isolates the change (R3 pins 1.7 / 0.3), and the gear cross-check follows R4 (A19 ADOPTED).

**Item 1, trip types.** The creel-database export (2022-23 to 2025-26, 37,460 rows) carries a `Trip Type` for every boat and jetty interview: crab only, `<fishery> & Crab`, finfish only, non fishing. `interview_combined.xlsx` is now BUILT from that export by `04_input_files/build_interview_combined.R` (the export is kept in `raw/`; the 17 legacy columns match the previous workbook row for row apart from source corrections made since July; `gear_tampered` is now filled, 13 Grays Harbor interviews of 2024-25 drop; 2025-26 extends to 2026-08-01). The contact classification reads the trip type (crabbing = crab-only + combo), is restricted to the launch sites the trailer and OSP counts measure, and the combo-trip share `c` is OBSERVED: 43 combos of 143 crabbing boats at the launches in 2024-25, none in Dec-Feb, half or more in the salmon and bottomfish months. In both Stan models `c` gets its own per-stratum logit walk (same order and innovations as `f`, its own step SD), observed per day by the combos among the typed crabbing boats; the OSP stream reads `f(1 - c[k])` and `f_lower_out = f(1 - c)` is the model's prediction of OSP's crabbing-only column, which is what the OSP data will check when they arrive (A21). Desk fit, boat all-gear 2024-25: `c` 0.02 in Dec-Feb rising to 0.4-0.8, `sigma_c` 1.05 [0.39, 2.40]; `f` Dec 0.92, Jan 0.89, Feb 0.87 (5 contacts), Mar 0.75, Apr 0.62, May 0.51, Jun 0.47, Jul 0.42, Aug 0.30, Sep 0.20; every non-f parameter within MC error of the Phase A control (max |z| 1.8); component 46,511 against 35,698 for Phase A, whose kappa = 20 prior and 20-contact floor sent five of ten months to 0.30.

**One choice this exposed, for Matt (D13).** Private boats interviewed at Float 20 / Float 17-21 (88 in 2024-25, finfish-heavy) and at the marina (34, all crab-only) are excluded from the classification by default because the counts measure boats launched at the ramp. If the dock-interviewed boats are trailered boats that stopped to unload, they belong in the classification and their exclusion biases the summer f up: July reads 0.47 at the launches against 0.31 with the docks, about +8% on the boat all-gear component. `contact_area_sensitivity.csv` (R0 desk) carries the monthly f under all three sets; `crab_fraction_contact_areas` selects.

**Item 1, sampler shift times.** `sampler_shifts.xlsx` (built from the survey export by `build_sampler_shifts.R`; 11 rows flagged for a 12-hour-clock check-out or a missing one, held out of the arithmetic) links to the interviews by survey id, and `interview_time` places every 2024-25 contact inside its shift (419 of 422). The shift-coverage diagnostic (A23) says: shifts run 09:44 to 15:52 (medians, 6.1 h) and cover about 75% of a day's boat returns (boat I/E profile, all launch sites pooled: WBL alone has 5 days), with 15% of returns after the last check-out; inside the shift the contact-hour distributions of combo, crab-only and finfish boats are indistinguishable (medians 12:59 to 13:18; non-fishing launches 11:41) and the crabbing share is flat by hour bin (0.42 to 0.58). That is the evidence available on the shift-time caveat; the returns outside the shift are unclassified, and only OSP's all-day count can say whether their mix differs. The return-time weighting is a documented hook (`crab_fraction_shift_weighting = "none"`), not implemented, because inside the shift it would move nothing and outside it there is nothing to weight.

**Two things the rebuild found (D14, D15).** Moored private boats (marina, docks) are in the CPUE interviews but in neither effort count, so the boat component omits their effort. And 2023-24 has no gear count in the export (`Number of Gear` blank, `Number Of Gear` a 0/1 flag), so under the deployment unit its 13,629 interviews have no CPUE denominator: with the missing 2023-24 effort counts (D8) that is a second blocker for the two-season run.

**Not yet run:** the ladder (`run_improvements_2026-09-08.R`, about 22 h; R0 desk ran clean with the new rows). Harness at 665 assertions; both Stan models compile under rstan 2.32.5.


## 1r. The inputs rebuilt from the per-season creel workbooks (2026-09-10)

Matt supplied the four per-season creel workbooks (`2223`, `2324`, `2425`, `2526 rec crab harvest data.xlsx`: the database exported one fishery season at a time, every sheet) and the I/E export, answered D13, and asked for the 2023-24 gear count to be found, every tab reviewed for modelling value, and every input carried through the most recent date.

**Every model workbook is now built from the season workbooks** (`04_input_files/raw/<YYYY><YY>_rec_crab_harvest_data.xlsx`; six builders, `build_helpers.R`, `build_all_inputs.R`; A24). The two pasted exports of 2026-09-09 are retired. Coverage is 2022-23 through **2026-09-08** for interviews (38,962 rows), effort counts (12,056; the 2026-07-16 workbook had 2024-25 only), shifts (3,077), the vessel tally (2024-25 and 2025-26) and a new charter trip roster; the holiday calendar runs 2022-23 through 2026-27 by one rule. The I/E workbook in the repository (through 2026-08-28) is newer than the export supplied and identical to it on every shared interval, so it stands. For 2024-25 the rebuilt rows reproduce the earlier workbooks (interviews: no difference on the legacy columns beyond database corrections; effort: none on 3,256 rows; tally: identical), so nothing fitted for 2024-25 changes; readers go through `read_input_workbook()` (B17).

**D15 was a misreading, corrected.** The 2023-24 gear count was in the export under the 2023-24 sheet's header spelling "Number Of Gear"; the builder had read only "Number of Gear" and I described the column as 0/1 flags, which it never was. Rebuilt from the season workbook every 2023-24 interview carries its count (0 to 23). With the 2023-24 effort counts (4,116 on 360 days) and holidays now present, the staged two-season run is unblocked on data except for its census: no vessel tally or charter roster was kept before 2024-25, so the 2023-24 component is 0 and the function warns that its frame is missing (D8). The pre-fit driver code and every Stan data set build for 2023-24, 2025-26 and the two-season config; a short real fit of the 2025-26 boat all-gear component ran to completion (below).

**D13 closed.** The dock-interviewed private boats are docked at the floats; their trailers sit in a lot that is not censused. The launch-site restriction stays the default, and D14 (moored boats outside the effort counts) is confirmed structural: the protocol's "Boats Entering Marina Count" was never recorded at Westport.

**What the tabs added (the review Matt asked for).** Modelling-relevant and now in the workbooks: the interview's fishing area (`bay_or_ocean` / `river_or_ocean`: Inside Bay against Ocean north or south of Point Chehalis, the basis for any bay-versus-ocean apportionment of the Westport boat catch), `boat_name` (per-vessel means for the census), `total_vehicles`, the released-crab fields (`crab_released`, `dungeness_returned` and reason, red rock), the gear label as recorded beside the harmonised one; every effort count the protocols recorded (`vehicle_count`, `boats_entering_marina`, `buoy_count`, `crabber_count`, `jetty_people_count`); the sampler's on-site conditions per survey (tide stage 2022-24, rain, cloud, wind, wind direction, the weather station site from 2024-25; `special_conditions` carries "Small Craft Advisory", "Bar Restrictions", "Halibut Opener", "Razor Clam Opener", the strongest candidates for a boat-effort covariate, unused for now); the charter roster; the survey-sheet vessel tallies of 2025-26. Reviewed and left out: the notes columns (free text with vehicle descriptions and personal details; kept in `raw/`), the per-season summary and pivot tabs (`boat summary`, `jetty summary`, the holiday recaps, `days sampled`, `Beach info`, `summary`: derived from the primary sheets, no input value), and the `location_key` / `gear_key` tabs, whose content is encoded in `build_helpers.R` (`area_map`) and `build_interview_combined.R` (`gear_map`).

**Three findings with a modelling consequence.**

1. *The gear vocabulary changed for 2024-25.* "Collapsible trap or ring" became "Ring Net", "Fishing rod with foldable trap" became "Trap (foldable, star)", "Fishing rod with snare" became "Snare". The gear-resolved classifier matches by regex, so on 2023-24 data every ring net (1,906 single-gear rows) was a TRAP. `gear_type` now carries the 2024-25 vocabulary in every season (`gear_type_raw` keeps the label as recorded).
2. *The Feb 2024 PFD stand-down.* From Feb 3 to 13, 2024 the samplers could not walk the floats and the Float 20 / 17-21 "gear counts" are the day's total gear from interviews ("Total from Interviews" in the notes; 123 on Feb 10 against a February mean of 14), a day's deployments rather than an instantaneous count. Twenty rows carry `qc_flag = "gear_count_from_interviews"` and the reader holds them out (`effort_qc_drop`, B16).
3. *The charter roster contradicts the census premise for charters* (A25, D16). The "charter trips" sheet lists every charter crab trip the operators reported. On 2024-25, 31 Westport trips sailed in the census window and **8 of them, all marked "missed", fell on days without a tally**; on three tally days the tally has one charter where the roster has two. Under the exact-sum census those trips were zero. `charter_frame = "roster"` (shipped) takes, per day, the larger of the roster's trips and the tally's charter column, on every day of the window: 2024-25 census **7,884 -> 8,592** (+708: 514 on the 8 roster-only days, 193 on the three disagreeing tally days), SE 131 -> 171; 2025-26 census 3,419 on Dec 1 to Jan 3. `charter_frame = "tally"` reproduces 7,884. The roster is Matt's own method (the sheet estimates the unobserved crab of the missed trips the same way), which is why it ships as the default rather than as an option to weigh; what it cannot do is test the premise for the commercial vessels, which have no roster (D16).

**The 2023-24 season on its own data.** Grays Harbor shifts were shorter (median 4.1 h against 6.0 h in 2024-25; 3.5 h in 2022-23), one count per site per day rather than three, and OSP covers it from 2024-03-09 only, so its boat stream leans on the trailer counts and the contact-based f covers less of the day. The 2025-26 season is complete through 2026-09-08 (the season ends Sep 15) except for OSP, which stops at 2025-10-18 (D18); in July and August 2026 the samplers contacted about 400 launched boats a month, so the summer f there (0.16) is measured to a precision the 2024-25 season never had.

**Not run:** the ladder (unchanged in structure; its R0 desk rung ran clean and now states both census numbers). Harness at 702 assertions; both drivers purl-parse; both Stan models unchanged since 2026-09-09.


## 1s. The charter component becomes an expansion; three holidays; a PE question (2026-09-11)

Matt's corrections after the 2026-09-10 rebuild, and one thing they exposed.

**The census statement applies to the commercial boats only.** "My statement of a complete census without error applies to the commercial boats. The charter vessels are not 100% sampled and need expansion." So `estimate_comm_charter()` now runs two estimators (A26). The commercial part is unchanged: the exact sum over the tally days, the unsampled days having no operation, with the per-vessel catch mean reported as a near-census sample mean and not carried. The charter part is an expansion over the charter TRIP frame: `N x` (mean catch per interviewed trip), stratified by vessel, variance `N^2 (1 - n/N) s^2 / n`. On 2024-25 that is 34 trips with 20 interviewed (59%), 1,286 crab observed and **2,133** expanded with SE 73; the census total goes 7,884 to **8,538** and, for the first time, part of the census SE (73) enters the port interval (`census_uncertainty = "charter"`). The estimator agrees with the spreadsheet's own per-vessel method to about 1%, which is the point: the sheet already expands the missed trips, sometimes at the observed mean and sometimes at crabbers x the daily limit, and the two rules disagree by up to 45% per vessel with no variance attached. `charter_frame = "tally"` reproduces 7,884 and the tally frame with `census_expansion = "day_type"` reproduces 11,753 exactly, so every earlier number is still recoverable.

**Three holidays added** (A27): Thanksgiving Day, Veterans Day (observed) and Juneteenth, on the same federal-observed-day rule as Independence Day, applied by rule to every season. The evidence was the samplers' own Holiday? flag against the Float 20 gear count: 2.4x, 2.3x and 1.7-2.3x the same month's weekday mean. Christmas Eve, Presidents Day and MLK Day stay out because their counts do not support it. Three days are re-typed per season.

**The ring-net question, answered.** There was no ring-net option before 2024-25. The 2022-23 and 2023-24 forms offered "Collapsible trap or ring" and, separately, "Fishing rod with foldable trap" and "Star trap", and the 2023-24 workbook's own `gear_key` sheet maps "Collapsible trap or ring" to "Ring net", which is what the builder implements. The caveat worth keeping: the old label is one option where the current form has two, and because the old trap option was rod-mounted, a hand-held foldable trap had nowhere else to go. The combined ring-plus-trap share of gear tokens is stable across the vocabulary change (43% then 45%) while its internal split moves from 66/34 to 47/53, which is a moved boundary rather than a change in gear use. It affects a gear-resolved fit of 2022-24 only (D20).

**What the holiday change exposed, and it is the biggest open number on this branch (D19).** The PE builds its effort strata as (week x day-type) from the SAMPLED days and left-joins the calendar, so a cell with calendar days but no sampled day expands at ZERO effort under the shipped `pe_empty_effort_stratum = "zero"`. That mechanism was flagged on 2026-08-25 as a thin-component problem; on the rebuilt 2024-25 inputs it is **45 of 289 shore all-gear days (15.6%)** and 48 of 289 boat all-gear days, three times the threshold that prints the run's own warning. The PE port total reads 72,224 under "zero", 85,243 under "day_type" (+18%) and 90,861 under the new "local_day_type" (+26%, B18): the empty cells sit in the HIGH months, so a month-local fill is larger than a sub-season one. Which is right cannot be settled against `20260904/pooled-CPUE-AD-A1-adopted`, because that run predates the derived shore turnover and the boat recentring; on that stale reference the shipped "zero" PE happens to sit closest to the BSS. The default is therefore unchanged and the decision belongs in the ladder. **This is a PE-only defect**: the BSS gives every calendar day the weekend (B1) and holiday (B2) effort effects and imputes unsampled days from its AR, so it never zeroes one. It matters because the PE is the fallback for a component whose gate fails and the reference the PE-vs-BSS table is read against.

**Verification:** harness 725 assertions, both drivers purl-parse, the R0 desk rung runs and reports the split census and the holiday re-typing, and the census reconciles (`census_daily.csv` sums to the component) in all four frame x expansion combinations.

---

## 1t. Set up for the full ladder run: the PE's thin cells, and proving the rungs comparable (2026-09-10)

Matt: *"Set up the files for a full ladder rung run. Ship with local_day_type option on. Review all files involved; make improvements. Set levers to support the full ladder rung run and promote comparability of the results. The full ladder rung run will only be for the 2024-25 season and will use the new method for estimating f."*

### The PE's unsampled cells: shipping the fill exposed two more defects behind it

The instruction was one lever. Shipping it honestly took three, because **`local_day_type` alone would have imputed 29% more shore effort at zero variance and then multiplied it by a season-average catch rate.**

**The CPUE half (B19).** The effort of an unsampled cell could be filled month-locally while its CPUE was always the sub-season-wide ratio-of-sums, so the seasonal gradient was counted once in the effort and again in the rate. The imputed BOAT cells are in the summer and the pooled boat CPUE is dominated by the high-CPUE winter: boat all-gear PE catch 42,841 under `"pooled"` against 37,018 under `"local"`, **-13.6%**. The shore is almost unmoved (29,705 to 29,737), so this is a boat correction. **This corrects the +26% figure in the 2026-09-11 D19 note to +17.8%**, and the decomposition is worth keeping because it is the shape of the error:

| effort fill | CPUE fill | PE port catch | PE port effort |
|---|---|---|---|
| `zero` | `pooled` | 72,224 | 45,175 |
| `local_day_type` | `pooled` | 90,861 | 56,793 |
| `day_type` | `local` | 81,160 | 53,846 |
| **`local_day_type`** | **`local`** (shipped) | **85,076** | **56,793** |

`"zero"` was omitting a fifth of the season's effort.

**The variance half (B20, D21), and this is the larger finding.** The stratum SE was `sqrt(N^2 sd^2 / max(n, 1))` with `sd` from the cell's own sampled days, and **`sd()` of one observation is NA, which the code replaced with 0.** With weekly strata and roughly 50% day coverage that is not a corner case:

| component | cells | unsampled | singleton | effort with NO variance | SE as reported |
|---|---|---|---|---|---|
| shore all-gear | 93 | 21 (45 days) | 37 (82 days) | 54% | 410 on 27,345 = **1.5%** |
| boat all-gear | 93 | 23 (48 days) | 34 (77 days) | 49% | 573 on 10,482 = 5.5% |
| shore pot-closure | 24 | 1 (1 day) | 10 (19 days) | 41% | 238 on 6,751 = 3.5% |
| boat pot-closure | 24 | 6 (11 days) | 7 (14 days) | 29% | 29 on 398 = 7.2% |

A 1.5% SE on a design whose within-cell CV runs 0.2 to 0.7 is not a measurement of anything. And imputed cells were treated the same way, so **the shore all-gear effort SE was bit-identical at 410 whether the fill added 0, 6,240 or 7,947 effort units.** Under `pe_variance = "impute_aware"` a singleton cell borrows its donor's spread with the divisor still 1 (the collapsed-stratum estimator standard for one-per-stratum designs: Hansen, Hurwitz & Madow 1953; Cochran, *Sampling Techniques* 3rd ed. 1977 sec. 5A.12; Wolter, *Introduction to Variance Estimation* 2nd ed. 2007 ch. 2) and an imputed cell carries the donor mean's sampling variance plus a between-cell surrogate, `s^2 (1/n_donor + 1)`. The point estimate does not move; the 2024-25 effort SEs go 410 to 1,379 (shore all-gear), 573 to 1,310 (boat all-gear), 238 to 550 and 29 to 79. **Note which half is bigger**: fixing the singletons alone takes the shore SE to 1,099 even under the retired `"zero"` fill, so most of the understatement was never about imputation. `"sampled_only"` restores the old arithmetic exactly, and every run reports `effort_se_sampled_only` beside `effort_se` whatever the setting.

**And it is not PE-only, contrary to what D19 said.** A component whose convergence gate fails reports its PE point in the port total, as a **constant with no interval** (section 7.7 of either driver). The thin boat pot closure is the component most likely to fall back, `bss_min_interviews` was lowered to 15 specifically to let it attempt a fit, and the honest expectation recorded at the time was that it may fail the gate anyway. So these levers can reach the headline, and calling them PE-only was wrong.

All three now live in **one** implementation (`03_R_functions/pe_effort_strata.R`, B21). Both PE runners had carried their own copy of the stratification, the fill and the SE formula, so every defect above existed in two places; on 2024-25 the pooled and gear tracks' PEs are now numerically identical component by component, which they were not before.

**What is still open (D22), both deliberately.** The PE applies **no finite-population correction** anywhere, although `estimate_comm_charter()` does, so the two estimators are inconsistent on that point; adding it would reduce the SE and would have mixed an opposite-signed change into B20. And the **PE catch has no SE at all** -- only effort does -- so the PE side of `pe_vs_bss_comparison.csv` is a point against a posterior interval, and a gate-failed component enters the port total with no uncertainty. The honest version is a ratio-estimator variance, a larger piece of work than everything in this patch.

### The ladder: comparability proven before the MCMC, not after

Four changes, none of them a new rung (B22).

**1. `F_METHOD`.** The ladder as designed fits R1 and R2 on the RETIRED f (one scalar at 0.30) precisely so R3a and R3 can attribute the f change. The instruction says the run uses the new f. Both readings are defensible and they cost differently, so it is a switch: `"new_throughout"` (shipped) gives R0, R1, R2, R4, R5 -- three pooled fits and one gear fit, about 12 to 14 h, with every rung's port total citable and the rung-to-rung deltas isolating the TURNOVER levers on a common f. What it costs is the f attribution and the factorization proof, and the optional `R2f` rung buys the proof back for one extra fit: R2 with the f block rolled back is the same R3-vs-R3a test, at half the cost and against the current turnover prior rather than the retired one. `"ladder"` restores the original seven-rung design. The mode is in every `run_tag` and every output filename, so the two can never land in the same verdicts file.

**2. The window pin was incomplete.** `WINDOW` pinned five keys and left four per-season keys to be inherited from the shipped two-season config, so **every 2024-25 rung would have run with `pot_open_date = 2023-12-01`.** That one turns out to be cosmetic, for a reason worth recording: `estimate_L_effective()` took `pot_open_date` as an argument and never referenced it, so the "I/E regression split" that `run_config` and CHANGE_REGISTER A14 both credited it with feeding **never existed**. The dead argument is removed, A14 is corrected, and the real multi-season approximation (the regression pools every season's I/E days into one `yday -> L` curve) is now documented at the function. `pot_closure_start`/`end` were right only by luck and are not cosmetic: `build_subseasons()` splits the season on them. All nine per-season keys are pinned, and `validate_season_window()` now warns when any of them falls outside the estimation window.

**3. The PE levers are pinned too**, for the reason above: they can reach the port total of any rung whose gate fails, and `run_config`'s defaults changed today.

**4. Comparability is now a preflight, not a post-mortem.** The runner resolves every rung's configuration before anything runs and fails if any non-delta key differs between rungs, if a rung differs from its comparison rung in a key it did not declare, or if a pinned key drifts; it writes `improvements_2026-09-08_manifest*.csv` with every rung's resolved levers and a config digest. Previously the only check was `config_delta()` reading `run_parameters.txt` **after** the fits, when the MCMC had already been spent. And `RESUME` now requires a matching digest (`IMP_STAGE.txt`) rather than reusing any folder with the right name and a `run_parameters.txt`; since this file's deltas changed on 2026-09-11 and again today, a partially-run ladder would otherwise have been silently completed with a mix of configurations.

**R0 now runs the PE.** The desk rung was called "the PE and reporting changes" and never ran the PE. It now runs all four unsampled-cell arms (a minute each) and reports the singleton shares and both SEs, so the D19 choice can be made against R4's own BSS instead of against the stale 2026-09-04 reference, at no MCMC cost.

**Verification:** harness **789** assertions (new section 59 recomputes the fill, the collapsed-stratum variance, the between-cell term and the CPUE fill by hand on a four-cell fixture, and asserts `"sampled_only"` reproduces the old arithmetic exactly); the legacy settings reproduce the pre-2026-09-12 PE to the digit (72,224 catch, 45,175 effort, SEs 410/573/238/29); both drivers purl-parse; the gear track's PE now matches the pooled track's exactly; the dry run is clean under `F_METHOD` = `"new_throughout"`, `"ladder"` and with `R2f` added.

---

## 1u. The two-pass run, and the two things that would have wasted it (2026-09-10)

Matt: *"set up the ladder run based on your recommendation: run the 4-rung version now and then follow up with the new `R2f` control rung."*

`LADDER_PASS` makes that one number: `1` fits R0/R1/R2/R4/R5 (~12-14 h, and every one of those rungs' port totals is citable), `2` adds R2f and re-runs, whereupon `RESUME` matches the four pass-1 fits by their config digest, skips them, fits R2f alone (~4 h) and recomputes every verdict from the folders on disk. Verified before shipping that adding R2f leaves the other four digests **unchanged**, so pass 2 cannot accidentally re-fit them.

Reviewing the plan before it runs turned up two defects that would each have cost real time.

**The f/c exclusion list the factorization verdict uses had a real gap, and finding out what it actually costs corrected my own first claim (D23).** `fit_agreement()` compares every shared row of two runs' summaries. `crab_bss_pooled.stan` declares `sigma_f_out`, `cfi_kappa_out`, `sigma_c_out` and `cfc_kappa_out` **unconditionally** and sets each to exactly `0.0` when its walk is off, and `sigma_c_out`, `cfc_kappa_out` and `z_c` were **never added** to the hand-typed exclusion list when the combo-share walk arrived on 2026-09-09 — in either of the two copies it existed in.

I first wrote that this would have made R2f's factorization verdict FAIL after four hours of fitting. **That was wrong, and the smoke fit below is what showed it.** Under the R2f configuration those quantities report `sd` 0 and therefore `se_mean` **NaN** — not 0; `n_eff` and `Rhat` are NaN too — and `fit_agreement()` skips any row whose combined `se` is not finite. So they are skipped with or without the list. The rows that *would* have been compared and would have produced an enormous z are `f_crab_out[*]` (a real posterior on **both** sides: 0.307 under R2f against the winter monthly values under R2), plus `E_sum` and `C_expected_sum` — and all three were already excluded before today.

The gap is therefore closed as **hardening, not a fix**: the masking rests on `rstan::summary()`'s NaN convention for a zero-variance parameter, an implementation detail rather than a property of the design. The durable part is not the three names but harness section 60, which asserts the regex covers every f/c quantity **both** Stan models report, asserts it does *not* cover the parameters the proof must compare (`B1`, `B2`, `tau_bar_out`, `R_G_boat_out`, `mu_mu_E`, `sigma_eps_C`, `sigma_IE_out`), and pins the NaN-skip behaviour the reasoning above depends on. The list maintains itself the next time an f output is added, and a future rstan that reports 0 instead of NaN cannot turn this into the failure it currently is not.

**A patch applied between the two passes would have invalidated the comparison, silently.** The config digest says a folder was built from this rung's *configuration*. It says nothing about the code. Demonstrated: appending a line to `crab_bss_pooled.stan` left every digest unchanged and `RESUME` still reused R2 — so R2f would have been fitted by one version of the pipeline and measured against an R2 fitted by another, with nothing on the page to say so. `IMP_STAGE.txt` now records a three-layer code fingerprint (Stan models / drivers / `03_R_functions`) and the rstan/StanHeaders versions. It is **reported, not enforced**: re-fitting 12 h because a comment changed in an R function would be worse than the problem. What it does instead is name the layer that moved in the RESUME table and at the rung, record a REVIEW row, and route every cross-rung claim through `V1cross()`, which downgrades PASS to REVIEW and says why — *the premise, not the result, is what failed*. Seven verdicts are wrapped: every `fit_exactness()` bit-identity claim and both `fit_agreement()` factorization proofs.

**Also checked, because under `new_throughout` nothing exercises it until R2f is reached twelve hours in: the legacy-f Stan path samples cleanly.** A 2-chain 300-iteration fit of the 2024-25 boat all-gear component under the R2f configuration (`crab_fraction_dynamic = FALSE`, `crab_fraction_strata = "none"`) built valid Stan data (`D` 289, `P_n` 10, `IntC` 130, `CFI_n` 0, `CFC_n` 0, `n_f_strata` 1, the retired Beta(6,14) at 0.30), passed `bss_assert_stan_data()` and ran in 0.7 min with every R-hat at 1.00: `f_crab_out` 0.307 (prior-driven, as it must be with `CFI_n` 0), `tau_bar_out` 2.99 (on the calibration, matching the 2.97 likelihood-only figure quoted for R2), `R_G_boat_out` 3.56. This is also the fit that produced the `se_mean` NaN result above, which is why the D23 claim could be corrected before it was acted on.

**Verification:** harness **817** assertions (section 60 covers the digest stability across the two passes, `F_EXCLUDE`'s coverage against both Stan models and its non-coverage of the parameters the proof compares, the pinned NaN-skip behaviour of `fit_agreement()`, the fingerprint's determinism and layer attribution, and the non-fatal rung loop); the code-drift warning fires and names `stan` when only a Stan model changes; a seeded set of pass-1 folders makes pass 2 report "RESUME will reuse 4 of 5 fitted rung(s); R2f will be fitted."

---

## 1v. THE LADDER RAN. The 2024-25 estimate, and what moved it (2026-09-11)

All five fitted rungs completed on 2026-09-11 (`LADDER_PASS <- 2`, R4 fitted first). Every code fingerprint matches across the five folders, so every cross-rung comparison is valid. **14 PASS, 7 INFO, 4 READ, 2 REVIEW**, and both REVIEWs are diagnosed below.

### The headline

| | port total | 95% CI | CI width / median | PE port |
|---|---|---|---|---|
| 2026-09-04 A1 baseline (superseded) | 72,027 | [53,018, 101,364] | 67% | 45,105 (-37% vs BSS) |
| **R4, the shipped configuration** | **94,376** | **[77,566, 118,602]** | **43%** | 85,076 (-10% vs BSS) |

**+31.0% on the median, and the interval tightened from 67% to 43% of it.** The PE-vs-BSS reconciliation went from a 37% gap to a 10% gap.

### What moved it, attributed

This is what the ladder was built to produce, and the parts sum to the whole within 24 crab:

| change | rung step | effect on the port |
|---|---|---|
| the dynamic monthly f (boat only) | R2f -> R2 | **+11,963** (boat all-gear +35.6% at a fixed turnover) |
| the derived shore turnover | R2 -> R4 | **+10,327** (shore all-gear x1.356, pot closure x1.413) |
| the boat turnover recentring | R1 -> R2 | +2,621 |
| the census split (commercial census + charter expansion) | in R4 | -3,283 against the baseline's 11,821, by design |
| the filters and the rebuilt inputs | A1 -> R1, shore | +0.1% shore all-gear, +0.4% pot closure |

### The design guarantees held exactly

- **Shore bit-identical across R1, R2 and R2f** (6,343 / 21,502 in all three, same divergence counts, same R-hat): the boat turnover prior and the f block have no path to the shore, and now that is measured rather than argued. 6,270 shared parameter rows identical at full precision.
- **Boat bit-identical between R2 and R4** (1,372 / 45,604): the shore turnover has no path to the boat. 4,170 rows identical.
- **THE FACTORIZATION PROOF IS ESTABLISHED.** R2f is R2 with the f block rolled back, so the pair differs only in f. Max |z| over 1,957 shared boat parameters is **3.20** (`L_raw[13]` of the pot-closure fit), with 2 rows above 3 (0.10%). f enters the boat generated quantities and nothing else; effort and CPUE posteriors are unchanged in distribution. This had never been run before.
- **The two tracks agree**: gear port 93,274 against pooled 94,376 (-1.17%, inside the 2% criterion), tau_bar to 0.02%, monthly f to 0.002.

### Sampler health: clean everywhere

Every component reports BSS in every rung -- **no PE fallbacks anywhere**, which also means the PE's unsampled-cell levers did not reach the reported total. Worst divergence fraction 2.17% (shore all-gear) against a 5% backstop; treedepth saturation 0%; every R-hat within 1.0007 of 1; n_eff 4,600 to 21,200 against a 400 floor; scale-aware divergence impact at most 0.003 posterior SD against a 0.10 threshold. `model_adequacy.csv` flags nothing overparameterised, nothing miscalibrated, no PIT bias; `p_loo` runs 7-18% of observations (the retired daily shore AR was 35%), and the only `flag_loo_unreliable` cases are one bad Pareto k each on shore all-gear and boat pot closure.

### D19 is settled, for the shore

**Every gap in this subsection is stated PE RELATIVE TO BSS**, so a negative number is a PE below the BSS. That is the opposite of `pe_vs_bss_comparison.csv`'s own `effort_diff_pct` / `catch_diff_pct`, which are BSS relative to PE and carry no direction in the column name (D27); the first version of this section quoted the two directions in one sentence.

**Shore all-gear: PE 35,292 against BSS 34,837, +1.3% on effort; PE 29,737 against BSS 29,210, +1.8% on catch.** Under the retired `zero` fill the effort comparison was **-21.5%** (PE 27,345); under `day_type`, **-3.6%** (PE 33,585). Two estimators with nothing in common -- a design-based stratified expansion with a month-local donor, and a Bayesian AR(1) state-space with day-type effects imputing every unsampled day -- landing within 1.3% is the strongest available validation of both, and it is what the `local_day_type` fill was adopted to achieve. Shore pot closure: **-2.9%** on effort, **-4.2%** on catch (PE 6,878 / 8,591 against BSS 7,086 / 8,963).

**For the boat it is not settled, and the evidence runs against the CPUE half of B19.** Boat all-gear PE 37,018 against BSS 45,604, -18.8%. The `pooled` CPUE fill would have given 42,841, **-6.1%** -- closer. The PE's own internal target points the same way: the boat interview ratio-of-sums is 3.276 crab per deployment, and the PE's implied CPUE is 2.650 under `local` (0.81x) against 3.067 under `pooled` (0.94x), with the BSS at 2.928 (0.89x). Against that, the theoretical argument for `local` (a month-local effort fill multiplied by a season-pooled rate counts the seasonal gradient twice) remains sound, and a correct month-local fill SHOULD pull the effort-weighted implied CPUE below an interview ratio-of-sums that is not weighted by the calendar. So this is genuinely unresolved: `local_day_type` for the effort fill is vindicated, `pe_empty_stratum = "local"` is not, and because no component fell back to PE the choice currently affects only the cross-check.

### The two REVIEWs

**1. R2, the boat rose 4.4% where the threshold expected 8-20%. The threshold was wrong (D25).** tau_bar moved 14.0% and the component moved 4.4%, because `mu_mu_E` fell 6.0% multiplicatively while `kappa_OSP` held still: the boat effort scale is jointly identified by the trailer counts, the OSP counts and the catch likelihood, so the prior move does not pass through. That is a better property than the threshold assumed. Band corrected, pass-through ratio now reported.

**2. R4, sigma_IE grew 55%, and this one is real (D24).** 0.373 -> 0.577 on the shore all-gear fit when the derived turnover was adopted; the pot-closure fit, which has its own I/E days, did not move (0.206 -> 0.202). Chasing it found that **`estimate_shore_turnover()` never filtered the I/E interval rows to the estimation window.** The shipped 2.477 rests on 40 days spanning 2023-08 to 2026-08, of which 6 are in the season and 4 in the sub-season it scales by 1.36 -- and it is the same number whichever season you run, while the guide described it as derived from the window. Three independent indications now point the same way:

| indication | value |
|---|---|
| in-window ratio of sums vs pooled | **2.225** vs 2.477 (**-10.2%**) |
| fitted posterior turnover | 2.394, between the two, with the log-SD 0.1 prior holding the rest |
| the 4 in-window I/E days | predicted arrivals 1.25x observed at tau 1.7, **1.64x** at 2.478 |

What this does **not** establish is that 2.225 is right: 6 days with a between-day log-SD of 0.20 gives a bootstrap log-SE of 0.068, so it is itself +/-15%, and both values sit far above the retired peak-based 1.669. **The direction of the 2026-09-09 adoption is well supported; the magnitude is worth about 3,500 crab (3.7%) on the port total.** `tau_shore_derive_window_only = TRUE` prices the alternative (it also needs `tau_shore_derive_min_days` lowered from 10, or the resolver rejects a 6-day derivation and falls back to 1.7). The field fix is a paired gear count on every I/E day.

### The monthly f, which is the largest mover, behaves as designed

| month | contacts | share | f | | month | contacts | share | f |
|---|---|---|---|---|---|---|---|---|
| 2024-12 | 31 | 0.97 | 0.92 | | 2025-05 | 17 | 0.47 | 0.51 |
| 2025-01 | 26 | 0.89 | 0.89 | | 2025-06 | 13 | 0.46 | 0.47 |
| 2025-02 | **5** | 1.00 | **0.87** | | 2025-07 | 17 | 0.47 | 0.42 |
| 2025-03 | 13 | 0.69 | 0.75 | | 2025-08 | 44 | 0.27 | 0.30 |
| 2025-04 | 26 | 0.62 | 0.62 | | 2025-09 | 29 | 0.14 | 0.20 |

No month sits at the retired 0.30 anchor; every informed month (n >= 20) tracks its own contact share within 0.06; the five thin months carry visibly wider intervals (mean width 0.360 against 0.249). February is the clearest demonstration: 5 contacts, all 5 crabbing, and the walk reports 0.87 rather than 1.00, shrinking toward January (0.89) and March (0.75) -- where the retired Beta(6,14) with kappa 20 would have pulled it toward 0.30. `sigma_f` 0.701, away from zero, so the walk is identified; f_Rhat within 1.0002, n_eff at least 10,108. **The flat f = 0.30 was wrong by roughly 3x in the winter months that carry most of the boat catch**, and that is the +35.6% R2f isolates.

### A caveat on the PE-vs-BSS agreement that the totals hide

The shore all-gear annual totals agree to 1.8%, but the MONTHLY distributions do not: the PE over-allocates January to March by about 2.2x and under-allocates June and July (0.36x, 0.59x), and the errors cancel. The annual agreement is real and worth having, but it is not month-by-month agreement and should not be cited as such.

**CORRECTED 2026-09-12 on the boat half.** This paragraph also said "the boat is worse (0.20x in June, 3.23x in September)". Most of that was an artefact of the reporting code, not a disagreement between the estimators: `pe_monthly_effort_share()` omitted the per-day crabbing fraction on the boat branch, which was exact for as long as `f` was the scalar 0.30 and stopped being exact the moment `f` became a month-varying walk. With `f` running 0.95 in December to 0.14 in September, the boat's monthly PE share was over-weighted **2.39x** in September and pulled to **0.36x** in December, which is the same shape and most of the size of the "disagreement" reported here. Fixed; the boat monthly PE-vs-BSS comparison needs re-reading on the next render before any statement is made about it. **No total moves**: the share is normalized and every component total comes from `run_pe_*()`, which always applied `f`.

### Two things the post-run patches themselves exposed (B24, B25)

The code fingerprint B23 added hashes raw file text, so the very patch that records this run flagged all five of its folders as "code changed" and would have downgraded four PASS verdicts to REVIEW on the next re-run -- for edits that cannot alter a fit. Comments and blank lines are now stripped before hashing, and `CODE_EQUIVALENT` carries an auditable fingerprint-to-reason map for changes whose default code path is provably identical. The 2026-09-11 run's fingerprint (`stan:6e4aca2a drivers:3de636d0 fns:44079dee`) is the list's first entry.

**And auditing that patch the next day found two defects in it, one silent and permanent (B25).** The equivalence entry was keyed on the RECORDED fingerprint alone, and `.code_delta()` returns early on that key, so the declaration excused its folder against **whatever the tree later became**: measured, `.code_delta("stan:6e4aca2a ...", a fingerprint with a changed Stan layer)` returned `""`. A Stan edit would have gone unflagged on all five committed rungs, forever, and every bit-identity and agreement claim above would have kept its PASS on a premise that had quietly stopped being true. The key now names both ends, so a declaration lapses when either moves. Separately, B24 changed the hash function itself, so the five rung stamps can never be re-emitted and differ in ALL THREE layers from anything the current function produces, including layers whose bytes never moved; `CODE_LEGACY` maps each to `code_fingerprint()` on the same tree (627a831: `stan:523f4e63 drivers:4c2ce454 fns:30ed14fb`, verified on a worktree) before any comparison. **Verification:** harness **836** assertions, including the lapse case that would have caught it.

---
