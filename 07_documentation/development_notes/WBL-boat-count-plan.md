# Plan of Attack: Using OSP Westport Boat-Launch Counts to Improve Private-Boat Effort/Catch

> **STATUS 2026-09-08: EXECUTED.** The OSP stream is production (`osp_scale_is_tau = TRUE`; the dense OSP series identifies the shared boat turnover `tau_bar`, adopted 2026-09-01), and the crab-only lower-bound machinery is built and inert pending the `WestportCrabOnlyEffort` column from OSP. See `CHANGE_REGISTER.md` items A2, A7, A8 and `PIPELINE_STATUS.md`. This document is the original plan, kept as a design record.

**Date:** 2026-07-28
**Author:** analysis session (for Matt George)
**Repo state:** `Coastal-Rec-Crab-BSS`, `main` @ `724eead` ("WBL count data")
**New input:** `04_input_files/WBL_boat_counts.xlsx` (OSP daily private-boat **total** counts, WBL, 2024–2025)
**Revision:** v2, 2026-07-28. `WestportPrivateEffort` confirmed to be a **daily boat total** (all private boats), which resolves the snapshot-vs-flow question; the duplicated 2025 rows are **fixed** (de-duplicated file, 459 → 306 rows, delivered as a drop-in replacement).
**Scope:** how to integrate the OSP boat counts and the crab-creel ingress/egress (I/E) crabbing-classification work into the private-boat effort/catch estimate; answers to the two design questions (static crabbing-% toggle vs. raw I/E in the model; automatic fallback to trailer counts outside the OSP window).

---

## 0. TL;DR

1. **Reframe.** There are two separate problems hiding in this task, and they should not be solved with one knob:
   - **(A) A denser, better effort *index*.** OSP counts are a high-signal **daily boat total** (median 24, max 329) that on their 61 overlapping days track the current trailer snapshot at **r = 0.98**. They attack the single most-cited weakness of the boat estimate (summer-extrapolation dominance, T1.1) and the weakest-identified boat parameter (`tau_boat`, GR-12).
   - **(B) A *directed-fishery fraction* the pipeline does not currently have.** OSP counts (and the current trailer counts) are **all private boats**, not crab boats. The crab creel's I/E classification is the only thing that can supply the crabbing fraction *f*. The current model has **no *f* at all**, it implicitly sets *f* = 1 (every trailer = a crab boat). Introducing *f* is therefore not just "how we use OSP data"; it is a **correction to a standing upward bias in the current boat catch**, and it will most likely move the boat number (and the port total) **down**, whichever effort series feeds the model.

2. **Q1, static toggle vs. raw I/E?** Neither in its pure form. Put *f* in the model as a **parameter with an informative Beta prior set from the I/E classification counts**, exposed through a `run_config` toggle exactly like `R_G_prior_mu` / `tau_boat_prior_mu`, with a pin-to-constant override for sensitivity. The static toggle is the right *first step and sensitivity tool*; raw I/E in the likelihood is the right *production* mechanism because **f is only identifiable from the I/E classification** (see §4.3). Allow *f* to vary by month/day-type as soon as the pilot supports it, a single annual *f* will over-state summer crab effort.

3. **Q2, automatic fallback to trailer counts when OSP is dark?** The instinct (trailer must carry the non-OSP window) is correct, but implement it as a **joint two-stream observation model**, not a conditional swap. The trailer series already spans the whole Sep→Sep window; OSP is a seasonal high-quality augmentation. Feed **both** as observations of the one latent boat-effort process, cross-calibrated on the 61-day overlap; the "fallback" then emerges continuously (tight intervals in the OSP window, appropriately wider in winter) with no boundary discontinuity. A hard swap either creates a scale jump at mid-October or, once you calibrate to avoid it, *is* the joint model minus the uncertainty propagation. **Caveat that bounds the ambition:** the OSP-dark window (mid-Oct→early-Mar) carries almost no boat crab effort (winter boat trailers: mean 1.3, max 5), so the fallback machinery buys *defensibility and correct winter intervals*, not a big number change. The headline mover is summer densification + the *f* correction.

4. **Data quality (fixed).** The 2025 rows in `WBL_boat_counts.xlsx` were **exact duplicates**, 153 dates each appearing twice with identical values (0 discrepant pairs); 2024 was single-row. Left as-is this double-weights every 2025 day in the effort likelihood and over-weights 2025 vs. 2024. **Fixed:** a de-duplicated file (459 → **306 rows**, 153 per year, every per-day value and the coverage identical, lossless) has been produced as a drop-in replacement. Still a standing design rule for the ingest: **NS = absent row = latent day the AR imputes; observed 0 = a real datum to keep**; OSP-dark days must be *missing*, never zero.

---

## 1. What the new file is, exactly

`WBL_boat_counts.xlsx` (sheet `Sheet1`), columns `Year, Month, Day, WestportPrivateEffort`. As delivered it held 459 rows; **de-duplicated to 306** (the 2025 dates were each duplicated). All stats below are on the de-duplicated data.

| Property | Value |
|---|---|
| Date span | 2024-03-09 → 2025-10-18 |
| Unique sampling days | **306** (153 in 2024, 153 in 2025) |
| Rows | **306** after de-dup (was 459: 2025 had two identical rows per day; fixed) |
| Month range | 3–10 only (mid-March → mid-October; OSP bottomfish-season dock coverage) |
| Value `WestportPrivateEffort` | integer **daily boat total** (all private boats); min 0, median 24, mean 56, max 329 |
| Observed zeros | 15 (real no-effort days, mostly March/October shoulders) |
| Missing values | 0 (missing days are simply absent rows = NS) |
| Weekend vs weekday | mean 72 vs 45; strong day-type signal (model already day-types) |
| Monthly means | Mar 7, Apr 7, May 26, Jun 40, **Jul 122, Aug 138**, Sep 50, Oct 9 |

**Semantics (resolved).** The field is a **daily boat total** (confirmed by Matt): the count of *all private boats* using the launch that day, across every fishery. It is **not** crab-specific and **not** an instantaneous snapshot. Two consequences follow. (1) It needs the crabbing fraction *f* to become crab effort (§4), exactly like the trailer series. (2) It sits at a **different scale from the trailer count**, which is an instantaneous snapshot: OSP is a full-day aggregate, so it attaches to the *expanded daily* boat effort (`E = lambda_E·tau`) while the trailer count attaches to the instantaneous latent rate `lambda_E`. The 61-day overlap (§3) bridges the two and, as a bonus, measures the within-day turnover that `tau_boat` encodes (§5.5).

---

## 2. How boat effort/catch is built today (and the two gaps this exposes)

**Chain (pooled and gear-resolved share it):**

```
boat_trailer_count (Westport Boat Launch, Ocean Shores)      fetch_crab_data.R:142-151
   → T_I[i] ~ NegBinomial2(lambda_E[day]/R_G_boat, r_E)      crab_bss_pooled.stan:403-407
        lambda_E = latent daily crab-gear-deployment rate, AR(1) over P_n periods
        R_G_boat ~ lognormal(log 4, 0.5)  (gear per boat group)   stan:370
   → E_boat = lambda_E · tau_boat                             bss_effort_spec.R:82-92 (L = tau_boat, ~1.2)
   → Catch  = E_boat · lambda_C                               stan:516  (lambda_C = crab per deployment, from crab interviews)
```

Authoritative boat numbers (2026-07-11 pooled, monthly AR): **catch 43,475; effort 14,716 deployments; CPUE 2.95**; port BSS total 83,914 (PIPELINE_STATUS §1).

**Gap A, no directed-crabbing fraction (the big one).** `fetch_crab_data.R` feeds `boat_trailer_count` straight in as `count_quantity` with **no crab filter** (L146), and the likelihood asserts `trailers = crab-boat groups = lambda_E/R_G_boat` (stan:400-407). Because `R_G_boat` (~4 crab pots/boat) is learned only from *crab* interviews, the model assigns ~4 crab pots to **every** trailer, including tuna/salmon/halibut/bottomfish boats that deploy zero crab gear. So latent crab effort is inflated by roughly (all boats / crab boats) = **1/f**, and boat catch scales with it. There is no `f`, `crab_fraction`, `directed`, or `target`-fishery term anywhere in the code or docs (grep of `07_documentation/`, `run_config.R`, both `.stan`). This is a genuine adaptation gap, not an oversight in the usual sense: the freshwater lineage it is built on (`docs_equations.Rmd`, the `BSS_creel_model_02*.stan` V/T/A/B streams) never needed a directed fraction because a river's boats are ~all the target fishery. A multi-fishery marine launch is different.

> **Consequence to confirm by run, not to take from my arithmetic:** applying an *f* < 1 will revise boat catch (and the port total) **downward**, in direct proportion to *f*. If the summer crabbing share is, say, 0.3–0.6, that is a large move on the dominant component. This is arguably **more consequential than the density improvement**, and it is orthogonal to the existing "boat BSS +97% over PE" debate (both PE and BSS use the same all-boat trailers, so both are biased together and their gap is unaffected, see §6).

**Gap B, the effort series is thin exactly where the harvest lives.** Boat trailer counts at the launch are sparse and low (whole-window median 0, mean 3.5). PIPELINE_STATUS's top credibility threat (T1.1) is that the boat harvest is summer-dominated and rests on extrapolation, and the boat catch is "proportional to `tau_boat`, whose prior rests on only 2 WBL I/E days" (GR-12; `next-steps-post-run5`). The documentation already names the fix: *"a more informative effort series (for example, access-point or camera exit counts)"* (pooled documentation L411). **The OSP series is that more-informative series**, and its dense months (Jul mean 122, Aug 138) sit exactly on the summer peak that currently drives the total on thin data.

---

## 3. The overlap: why a joint model is well-posed

Joining OSP (de-duplicated) to the current Westport Boat Launch trailer counts on shared dates:

| Quantity | Value |
|---|---|
| Overlapping days | 61 |
| Correlation (trailer snapshot vs. OSP) | **0.983** |
| OLS: trailer ≈ | **0.367 · OSP − 0.31** (through-origin slope 0.365) |
| OSP-zero days | 5, of which trailer also 0 on **5/5** |

Two series, two agencies, two protocols, **near-proportional through the origin with r = 0.98 and clean zero-agreement.** That is close to the ideal setup for treating them as two noisy observations of one latent effort process with a fixed scale ratio. Because OSP is a daily total and the trailer is an instantaneous snapshot, the ~0.37 factor is a direct read on **within-day boat turnover**: one trailer count catches ~37% of the day's boats, i.e. ~2.7 boats pass through the launch per boat present at the count instant. That is exactly the turnover `tau_boat` encodes (§5.5). This overlap is the empirical bridge that makes both the joint-effort model (§5) and the fraction calibration (§4) identifiable rather than assumed.

---

## 4. Q1, the crabbing fraction: static toggle, raw I/E, or hybrid?

### 4.1 The three options

| Option | What it is | Pros | Cons |
|---|---|---|---|
| **A. Static toggle** | `crab_fraction` in `run_config`; multiply OSP (and/or trailer) effort by it | Trivial; transparent; you own the number; ideal sensitivity lever | Throws away classification uncertainty; can't be validated/updated; a single value is biased if *f* varies seasonally; inconsistent with how the pipeline treats every other nuisance (all have priors + propagation) |
| **B. Raw I/E in the model** | Feed per-stratum classification counts (n_crab, n_total) as a Binomial likelihood; *f* is a parameter | Propagates classification uncertainty into the boat total; updates as the pilot grows; can be time-varying; statistically coherent with the BSS | More Stan/prep work; needs enough classified days to bind; identifiability must be handled (§4.3) |
| **C. Hybrid (recommended)** | *f* is a **parameter with a Beta prior set from the I/E counts**; `run_config` exposes the prior mean/precision and a pin-to-constant override | Best of both: degrades to "you set it" when data are thin (prior-dominated, transparent), sharpens to "data set it" as classification days accumulate; mirrors `R_G`/`tau_boat`/`L_effective` exactly | Same as B; requires deciding the stratification of *f* |

### 4.2 Recommendation

**Adopt C.** It is the same pattern the project already uses three times over: `R_G` (data-driven lognormal prior, `prep_bss_crab_pooled.R:179-189`), `tau_boat` (prior-centered parameter, `bss_effort_spec.R:90-91`), `L_effective` (regression prior with propagated `L_sigma`, `bss_day_length.R`). A point `crab_fraction` would be the *only* first-order multiplier on the headline boat number carried with **zero** uncertainty, indefensible next to the care taken on `tau_boat`'s ±0.3 log-SD. Concretely:

- Model: `f ~ Beta(a0 + n_crab, b0 + n_noncrab)` per stratum (or a logistic-with-prior if you want covariates), entering the boat effort observation as a multiplier on the OSP/trailer mean. `a0,b0` a weak prior (e.g. Beta(1,1) or a mildly informative center).
- Config, mirroring the existing toggles: `crab_fraction_prior_mu`, `crab_fraction_prior_kappa` (precision), and `crab_fraction_fixed = NA|value` to pin it for sensitivity, the same "unset = data-driven, set = override" contract as `R_G_prior_mu`.
- **Start simple, sharpen later:** one *f* (or two: summer vs. shoulder) until the pilot has enough classified days for a monthly/day-type *f*.

### 4.3 The identifiability point (why this is not merely stylistic)

With OSP total-boat counts as the effort index, a **constant** *f* enters as a single multiplicative scalar on the effort observation. That scalar is **confounded** with the OSP catchability / `R_G_boat`: the effort+catch data cannot separate "few boats, all crabbing" from "many boats, some crabbing." **The I/E classification is the only information that pins *f*.** Therefore:

- The static toggle is really *"assert f and accept the confound."* Fine for a sensitivity sweep; not fine as the production number, because the boat total then scales linearly with a number you asserted.
- Putting the I/E counts in the model (B/C) is *"estimate f from the one data source that identifies it."* That is the whole reason to bring I/E in rather than a slider.
- **Corollary:** a *time-varying* *f* is partially identifiable even without perfect I/E coverage, because its *shape* (summer-low, shoulder-high) is informed both by the classification and by the mismatch between the OSP total-boat seasonality and the crab-interview seasonality. Still, lead with the I/E data; don't lean on the effort/catch data to infer *f*.

### 4.4 Should *f* be time-varying? Almost certainly yes

Westport's summer boats are heavily tuna/salmon/halibut/bottomfish; crab's *share* of boats is very plausibly **lower in July–August than in the fall/late shoulder**. Because the harvest is summer-dominated, a single annual *f* (which would sit above the true summer *f*) would **over-state summer crab effort**, re-introducing a bias in the same place the model is already weakest. So the I/E pilot should **stratify its classification by month and day-type** (§7). Begin with a single *f*; move to monthly/day-type as data allow.

---

## 5. Q2, the OSP window is seasonal: fallback, or joint model?

### 5.1 The coverage arithmetic (one crab window, Sep 16 2024 → Sep 15 2025)

| | Days |
|---|---|
| Crab estimation window | 365 |
| OSP-covered (within window) | **148** (Sep 18–Oct 19 2024, then Mar 8–Sep 15 2025) |
| OSP-dark (mid-Oct → early-Mar) | **217** |
| Trailer-count coverage | **entire window** (Sep 17 2024 – Sep 13 2025) |

So it is **not** "primary series with holes + fallback." It is **always-on cheap index (trailer) + seasonal high-quality index (OSP)**. That reframing is the key to the answer.

### 5.2 Recommendation: joint two-stream, not a conditional swap

Feed **both** series as observations of the single latent boat-effort process, each with its own scaling coefficient and overdispersion, cross-calibrated on the 61-day overlap. This is exactly the **index + census/tie-in** pattern already present in the inherited freshwater Stan (`BSS_creel_model_02*.stan`: separate `V`/`T`/`A` index streams **and** a `B` census/tie-in stream). The crab model simply collapsed to a single boat stream; re-widening it is within the original design, not a new invention.

Under the joint model the "fallback" is **emergent and continuous**: in the OSP window the latent effort is pinned tightly by the dense high-signal OSP series; in the OSP-dark winter it is carried by the trailer stream with appropriately wider credible intervals; there is **no discontinuity** at the mid-October boundary because it is one AR process and the OSP↔trailer scale was learned on the overlap.

### 5.3 Why not a hard swap (challenging the simpler instinct)

| Hard swap ("use OSP if present, else trailer") | Problem |
|---|---|
| Two series on different scales (0.37 ratio) | Creates a **scale discontinuity** at mid-October unless you first cross-calibrate, and once you cross-calibrate, you have built the joint model minus its uncertainty propagation |
| Discards trailer during the OSP window | Loses the overlap that calibrates the two series **and** an independent cross-check of the effort model |
| Treats winter trailer as good as summer OSP | Does **not** widen intervals for the thin trailer-only winter; understates uncertainty exactly where it is largest |

### 5.4 The honest counter-weight: winter stakes are low

Boat crab effort in the OSP-dark window is nearly nil, winter (Dec–Feb, pots legal) boat-launch trailers total **76 across 60 sampling rows, mean 1.3, max 5**, vs. summer (Jun–Sep 15) **1,417, mean 27, max 122**. So the elegant winter-fallback behavior buys **defensibility and correct winter intervals, not a materially different total.** Do the joint model because it is the correct and only-slightly-harder thing and it is the right home for the *f* correction and the summer densification, but do not oversell the winter handling as a headline mover. If you want a strictly minimal first cut, §7 Phase 1 gets ~80% of the benefit with a much smaller change.

### 5.5 A bonus: OSP can help identify `tau_boat`

`tau_boat` (boat catch is proportional to it) currently rests on 2 winter WBL I/E days and its prior (GR-12). Now that OSP is confirmed a **daily total** and the trailer a **snapshot**, their overlap ratio (~2.7 = 1/0.37) is an empirical measure of within-day boat **turnover**, the very quantity `tau_boat` encodes. Feeding both streams in at their correct attachment points (OSP to the expanded daily effort, trailer to the instantaneous rate) lets the overlap inform that turnover directly. One caveat on separability: the overlap pins the *product* of turnover and 1/*f* (a day's total boats vs. a snapshot of crab boats mixes both), so `tau_boat` is fully identified only once the I/E classification pins *f* (§4.3). The two data sources are complementary; together they can move `tau_boat` from prior-dominated to data-identified, addressing GR-12.

---

## 6. How this lands against the existing backlog

- **T1.1 (summer-extrapolation dominance, the #1 credibility threat):** directly mitigated. Dense OSP summer observations replace extrapolation with data on the months that drive the total.
- **T1.1(b) (external validation, open, "needs a benchmark the repo does not contain"):** OSP is an **independent** effort measurement (different agency/protocol). Agreement between the trailer-driven latent effort and OSP on the 61 overlap days is the closest thing to an external effort check the project has. Report it as one.
- **GR-12 (`tau_boat` weakly identified):** addressable once OSP (turnover, now that it is a confirmed daily total) and the I/E fraction are both in the model (§5.5).
- **P0 / "boat BSS +97% over PE":** **orthogonal to *f*.** Both PE and BSS use the same all-boat trailers, so the *f* correction moves them **together** and does not explain the BSS-vs-PE gap. Keep the two debates separate; do not let a reviewer conflate them. (It does mean the *absolute* level of **both** is likely high today.)
- **Method-of-record note:** the *f* correction changes the estimand's *level*, so it must be versioned deliberately (it is not inference-neutral). Treat it like the effort-unit change: isolate, run, compare against the confirmed 43,475 baseline.

---

## 7. Phased plan of attack

Ordered to respect the project's "validate by run, one change at a time" rule (the "pin lesson," PIPELINE_STATUS §5). Each phase is a separate run compared against the confirmed baseline.

| Phase | Change | Effort | What it buys | Gate / check |
|---|---|---|---|---|
| **0. Data hygiene + EDA** | De-dup 2025 (**done**, 306 rows); formalize NS=missing / 0=kept; build the OSP↔trailer overlap calibration (the §3 regression) as a committed diagnostic. Semantics resolved: OSP is a daily boat total | Low | Correct inputs; the empirical scale ratio | 306 unique days; overlap r reproduced; zero/NS map audited |
| **1. OSP as a second effort index (minimal)** | Add an OSP boat-count stream to the boat effort likelihood over Mar–Oct, same NB2 form as `T_I`, with its own scale + overdispersion; **no swap logic** (trailer already covers winter) | Low–Med (one Stan stream + prep, mirrors `Gear_I`/`T_I`) | ~80% of the benefit: dense summer effort, tighter summer intervals, T1.1 mitigation | Boat effort intervals shrink in summer; convergence still passes; overlap-day posterior effort consistent with OSP |
| **2. Crabbing fraction *f* (hybrid prior)** | Add `f` parameter + Beta prior from I/E counts (`crab_fraction_prior_mu/_kappa`, `crab_fraction_fixed`); apply as the directed-effort multiplier | Med | Corrects Gap A; boat catch onto a directed-crab basis (expect **downward**) | Sensitivity sweep over fixed *f* ∈ {0.3,…,1.0}; posterior *f* vs. prior; boat total vs. baseline documented |
| **3. Time-varying *f* + `tau_boat` from OSP** | *f* by month/day-type once the pilot supports it; attach OSP as a daily total so its overlap turnover informs `tau_boat` (with *f*, addresses GR-12) | Med | Removes the summer-*f* bias; data-identified turnover | *f*(month) intervals; `tau_boat` posterior no longer prior-pinned; tau-sensitivity (`diagnose_tau_boat_sensitivity.R`) flattens |
| **4. Consolidate + document** | Fold into both pooled and gear-resolved (shared effort model); update method-of-record; PE parity (PE must use the same *f* and OSP index so it stays a fair cross-check) | Med | One estimand, both tracks; publishable | Pooled ≈ gear-resolved boat still reconciles; PE/BSS comparison stays unit- and *f*-consistent |

**Config surface (all mirroring existing patterns in `run_config.R`):**

- `boat_effort_file = "WBL_boat_counts.xlsx"`, `boat_effort_source = c("trailer","osp","both")` (default `both`); `osp_effort_areas`/location key like `boat_launch_areas`.
- `crab_fraction_prior_mu`, `crab_fraction_prior_kappa`, `crab_fraction_fixed = NA`, `crab_fraction_strata = c("none","season_block","month","day_type")`.
- Keep `use_boat_ie`, `tau_boat_prior_*` as-is; wire OSP→`tau_boat` behind a flag (OSP is a daily total, so it attaches to the expanded daily effort, not the snapshot rate).

---

## 8. Risks, open decisions, and data-quality items

1. **[Fixed] 2025 duplicate rows.** De-duplicated (459 → 306; the duplicate pairs were identical, so the collapse is lossless and the 2025 doubling was a copy artifact, not two sub-daily counts). A drop-in `WBL_boat_counts.xlsx` was produced. Only residual action: make sure the next OSP data export does not re-introduce the duplication.
2. **[Resolved] `WestportPrivateEffort` is a daily boat total** (confirmed). It attaches to the expanded daily boat effort (→ `E = lambda_E·tau`), not the instantaneous rate the trailer snapshot feeds; the joint model must carry the extra within-day turnover between the two, which the 61-day overlap calibrates (§5.5).
3. **[Confirm] Trailer counts are all-boat**, the 0.98 correlation, 0.37 through-origin slope, and the absence of any crab filter in `fetch_crab_data.R` all say yes, but the I/E pilot will settle it definitively. If, against expectation, samplers were already recording crab-directed trailers, then *f* applies to OSP only and Gap A shrinks. This single fact governs how large the *f* correction is; treat it as the pilot's first deliverable.
4. **[Design] I/E pilot must yield stratified classification counts**, for *f* to be identifiable and time-varying, the pilot needs (n_crab, n_total) **by month and day-type**, ideally at the same launch/among the same boats OSP counts, so *f* multiplies the right denominator. A handful of pooled days gives only a single prior-dominated *f*.
5. **[Watch] Level change is not inference-neutral**, *f* < 1 changes the published boat total. Version it explicitly against the 43,475 baseline; do not slip it in as a "fix."
6. **[Watch] Don't conflate with the BSS-vs-PE gap** (§6): *f* moves PE and BSS together.
7. **[Minor] A third private-boat signal exists**, `wes_commercial_tally.xlsx` `private_tally` (census-window, currently plot-only). Worth a glance as an extra cross-check, not a primary input.

---

## 9. One-paragraph answer to the two questions as posed

**"Would a static estimate of the crabbing percentage (a config toggle) be useful, or is raw I/E in the model better?"** Build it as a **model parameter with an I/E-derived prior**, the toggle and the raw-I/E approach are the two ends of one design, and the right production choice is the hybrid that starts as your toggle (prior-dominated, transparent) and becomes I/E-identified as the pilot grows, because **the crabbing fraction is only identifiable from the I/E classification** and it multiplies the headline boat number linearly, so it must carry propagated uncertainty rather than be asserted. Keep the static toggle, as the prior mean and the sensitivity lever, not as the estimate. **"Given OSP isn't year-round, is automatic fallback to trailer counts the right approach?"** The instinct is right but implement it as a **joint two-stream model** (trailer always on, OSP as a seasonal high-quality index, cross-calibrated on their r = 0.98 overlap) so the fallback is a continuous emergent property with correct uncertainty, not a conditional swap that risks a scale jump at the window edge, while remembering the OSP-dark window holds almost no boat crab effort, so the real prize is summer densification and the crabbing-fraction correction, not the winter handoff.
