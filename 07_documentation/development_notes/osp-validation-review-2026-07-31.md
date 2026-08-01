# OSP boat-count validation: review of all 14 runs (2026-07-31)

**Branch:** `OSP-boat-count-incorporation`. **Source of truth:** `05_output/osp_validation_summary.csv` plus the per-run CSVs under `05_output/{20260729,20260730,20260731}/`. Numbers below are medians, Dungeness Kept (matching the summary and the `port_total` `BSS_median` convention); credible intervals and coefficients of variation use the posterior mean and 2.5/97.5% from `bss_summary_private_boat_all_gear_Dungeness_Kept.csv` and say so. Convention: no em dashes.

## Verdict

All 14 runs completed and converged. The batch is mechanically clean: the six validation gates all pass, both models reconcile at every step, and nothing in the OSP or f machinery leaks into the shore mode. The two questions the batch surfaced have now been decided:

1. **The crabbing fraction f is the largest single lever on the private-boat harvest, and 0.3 is still a placeholder.** The boat catch scales exactly linearly with f (0.2 to 8.7k, 0.5 to 21.6k, 1.0 to 43.2k). Production ships f = 0.3 pending the WBL egress pilot, which is the only thing that turns 0.3 into data. This is the one remaining open input.
2. **The OSP-vs-trailer within-day turnover conflict is resolved.** With OSP on and its scale free, the model estimates that scale (`kappa_OSP`) at 3.15, 95% [2.50, 3.91], while the trailer-based `tau_boat` sat at about 1.2. The crab-creel trailer count was confirmed (2026-07-31) to be an instantaneous snapshot, so the 2.7 to 3.0 turnover is real and `tau_boat` about 1.2 was a roughly 2x under-count of boat effort. Production therefore adopts `osp_scale_is_tau = TRUE`, which raises the boat about 2x (Step 6 boat about 27.7k vs Step 3b's 13.8k).

**Production configuration (adopted 2026-07-31):** `use_osp_boat_counts = TRUE`, `use_crab_fraction = TRUE` with `crab_fraction_set = 0.3`, `osp_scale_is_tau = TRUE`. Production private-boat harvest about 27,700 Dungeness kept (pooled), about 27,500 (gear-resolved).

## The run matrix

Medians, Dungeness Kept. "boat" = `private_boat (All gear)`; "port" = whole-fishery total. All 14 rows have `boat_pass = TRUE`; boat divergence fraction is 0.13 to 0.82% throughout.

| Step | Model | OSP | f | tau | Boat catch | Boat effort | Port catch (BSS) | f_out | kappa_OSP |
|---|---|---|---|---|---|---|---|---|---|
| 1 baseline | pooled | off | off | off | 43,333 | 14,716 | 83,599 | 1 | (prior) |
| 1 baseline | gear | off | off | off | 43,120 | 14,774 | 82,825 | 1 | 3.145* |
| 2 osp | pooled | on | off | off | 47,260 | 16,243 | 87,560 | 1 | 3.147 |
| 2 osp | gear | on | off | off | 47,126 | 16,266 | 86,817 | 1 | 3.147 |
| 3 f030 | pooled | off | 0.3 | off | 12,597 | 4,316 | 52,131 | 0.300 | (prior) |
| 3 f030 | gear | off | 0.3 | off | 12,740 | 4,314 | 51,592 | 0.302 | 3.128* |
| 3b ospf | pooled | on | 0.3 | off | 13,787 | 4,763 | 53,295 | 0.300 | 3.144 |
| 3b ospf | gear | on | 0.3 | off | 13,800 | 4,792 | 52,634 | 0.300 | 3.141 |
| 4 fix020 | pooled | off | 0.2 pin | off | 8,682 | 2,952 | 48,045 | 0.200 | (prior) |
| 4 fix050 | pooled | off | 0.5 pin | off | 21,644 | 7,364 | 61,383 | 0.500 | (prior) |
| 4 fix100 | pooled | off | 1.0 pin | off | 43,193 | 14,742 | 83,590 | 1.000 | (prior) |
| 5 fmonth | pooled | off | 0.3 x10 | off | 12,932 | 4,397 | 52,399 | ~0.300 | (prior) |
| 5 fmonth | gear | off | 0.3 x10 | off | 12,886 | 4,392 | 51,751 | ~0.300 | 3.130* |
| 6 osptau | pooled | on | 0.3 | on | 27,668 | 9,745 | 67,223 | 0.301 | (prior) |
| 6 osptau | gear | on | 0.3 | on | 27,524 | 9,748 | 66,447 | 0.301 | 3.153* |

\* When OSP is off or `osp_scale_is_tau` is on, `kappa_OSP` is not in the likelihood and reverts to its prior (mean ~3.14, 95% [1.66, 5.43]). It is data-informed only in Steps 2 and 3b (95% narrows to about [2.50, 3.91]). This is the correct, expected behavior and is itself a useful check. **Step 6 is the production configuration.**

## Part 1: mechanical validation, every gate passed

**Behavior-neutrality holds.** Step 4 `fix100` (f pinned to 1.0) reproduces the f-off boat level to 0.1% (boat catch 43,193 vs baseline 43,333; effort 14,742 vs 14,716). Step 5 with `strata = "month"` and no pilot data returns the scalar-0.3 totals with `f_crab_out` a length-10 per-month vector all at about 0.300.

**f is exactly linear and CPUE is invariant.** f acts as a pure post-hoc multiplier on boat catch and effort, leaving CPUE untouched, exactly as the "generated quantities only, decoupled from sampling" design promised: f = 0.2/0.3/0.5/1.0 give boat catch/effort ratios to the f-off baseline of 0.200/0.300/0.500/0.999, with CPUE holding at 2.96 to 2.97 across the whole 5x range. No convergence cost (every f-run under 0.85% divergences).

**Shore is fully insulated.** Shore catch and effort are byte-identical (`C_sum` 20,821, `E_sum` 24,162) across every batch run, OSP on, f on, tau on, f swept. The fits are seeded, so any leak of the boat-only machinery into shore would show as a changed number; it does not.

**Convergence is uniform.** All 14 boat fits pass (divergence fraction 0.13 to 0.82%, Rhat ~1.00, n_eff in the thousands); both shore fits pass in every run. The one `pass_convergence = FALSE` row in every `convergence_report.csv`, including the baseline, is `private_boat_ring_net_only`, which lacks the data for a BSS fit and falls back to PE by design.

**Both models reconcile.** Pooled and gear-resolved agree within 1 to 2% at every step (Step 2 boat 47,260 vs 47,126; Step 3 12,597 vs 12,740; Step 6 27,668 vs 27,524).

## Part 2: OSP raises the boat about 10% and tightens it modestly

Turning OSP on (Step 2, f off) moves the boat up about 10% in both models (pooled effort 14.7k to 16.2k) and improves precision modestly (boat-effort CV 12.2% to 11.2% pooled, 11.9% to 11.0% gear; roughly 8% relative, 16% in variance). The season-total 95% interval width barely changes because the level rose at the same time. The level rise is prior-influenced: with `osp_scale_is_tau = FALSE` the free `kappa_OSP` (prior mean 3.0) absorbs most of the OSP/trailer discrepancy and `lambda_E` rises about 10% to reconcile. The more durable value of OSP is the 87 "OSP only" days (Part 4) that carry a direct effort observation the trailer series lacked.

## Part 3: the turnover conflict, resolved

When OSP is on and free, the model's posterior for `kappa_OSP` is 3.15, 95% [2.50, 3.91] (both models, Step 2), much tighter than the prior [1.66, 5.43], so the 61 overlap days genuinely inform it. `kappa_OSP` and `tau_boat` are structurally the same quantity (within-day boat turnover = OSP daily total / instantaneous boats), and they disagreed by about 2.6x. The raw overlap agrees with the model and shows the answer depends on what the trailer count measures (`osp_trailer_overlap_calibration.csv`, n = 61 paired days):

| Trailer metric treated as the snapshot | corr | implied turnover |
|---|---|---|
| mean per visit | 0.965 | 3.03 |
| max per day | 0.983 | 2.74 |
| sum per day | 0.899 | 2.01 |

**Resolution (2026-07-31): the crab-creel trailer count is an instantaneous snapshot.** Under that protocol the OSP/trailer ratio is the true within-day turnover, so 2.7 to 3.0 is real and `tau_boat` about 1.2 was a roughly 2x under-count of boat effort. `osp_scale_is_tau = TRUE` therefore corrects a real bias rather than importing an artifact, and is adopted for production. The consequence is the Step 6 doubling: boat catch 13,787 (Step 3b) to 27,668 (Step 6), 2.01x, in both models. A corroborating signal: the boat `distortion_E` diagnostic rises from about 0.008 (baseline) and 0.018 (Step 2) to about 0.10 (pooled) and 0.16 (gear) in Step 6, the fingerprint of OSP overriding the weak trailer-based `tau_boat` prior. (The metric name is read here, not its exact internal definition, so treat that as corroboration; the 2x level change is not in doubt.)

**One caveat retained honestly.** The exact multiplier assumes the trailer snapshot is timed representatively. If it is taken at a fixed near-peak time, occupancy is overstated and the ratio understates true turnover, so 2.7 to 3.0 is a floor rather than a ceiling (the correction would be at least this large, not smaller). Confirming the snapshot timing tightens the multiplier; it does not change its direction.

## Part 4: coverage and the trailer fallback, confirmed

From `osp_coverage_audit.csv` (365 days): OSP operates 245 days (mid-March to mid-October) and is present on 148 days, every one inside the operating window and zero outside. Within the window: 61 OSP+trailer days (these identify `kappa_OSP` / `tau_boat`), 87 OSP-only days (OSP is the sole effort observation), 49 trailer-only, 48 neither. The 120 non-operating days run on trailer-only or nothing, so the model degrades to the pre-OSP behavior with no special-casing: the fallback is an emergent property, not a switch.

## Part 5: scope of the batch, and the two interview fixes that ride along

- Every f-run uses the set value 0.3 because the WBL egress pilot has not delivered per-stratum counts; Step 5's per-month vector is therefore all 0.3. The time-varying machinery is validated as plumbing, not as an estimate, and cannot bind until the pilot classifies enough boats per stratum.
- The f sweep (Step 4) is pooled-only; defensible because f is a purely multiplicative post-processing step with no sampling interaction and both models already agree at the f = 0.3 point.
- OSP precision was assessed at the season-total level (about 1 point of CV). The daily/monthly in-window tightening on the 87 OSP-only days was not quantified; that is the natural next diagnostic if OSP is to be justified on precision rather than level.
- **The two interview data fixes postdate these runs.** Any interview with `number_of_gear == 0` (no crab gear deployed, regardless of trip-completion status) and any gear-tampered interview (`gear_tampered == 1`) are now dropped in both readers. On the 2024-25 Grays Harbor set the non-crabbing rule targets 233 rows but only 2 currently survive the existing `crabbers > 0` / `min_fishing_time` filters, and every targeted row carries zero Dungeness (including the 10 incomplete-trip `completed_trip == 0` rows added when this was widened from the original complete/blank-only spec); `gear_tampered` is all blank today. So the net effect on the numbers above is about 2 interviews (negligible), but a confirmatory production run with the filters in place is the clean certification, consistent with the "validate by run" discipline.

## Part 6: production recommendation and the open decision

**On OSP and tau.** Adopted for production: `use_osp_boat_counts = TRUE`, `osp_scale_is_tau = TRUE`. The instantaneous-snapshot confirmation makes the turnover correction real, so leaving tau off would knowingly ship a roughly 2x low bias on boat effort. Production boat harvest about 27,700 (pooled) / 27,500 (gear-resolved). The residual work is to confirm the snapshot is timed representatively (Part 3); it can only raise the correction, not reverse it.

**On f.** Ship f = 0.3 with the caveat surfaced prominently and the Step 4 sweep shown as the sensitivity band, because it is the least-wrong point estimate and the sweep already brackets it. This is a placeholder resting on zero observations and driving the entire boat harvest linearly; the WBL egress pilot's f_hat is the real unblock, and time-varying f by month can bind only once the pilot classifies enough boats per stratum.

**The one remaining question to route to the field/ops side:** when will the WBL egress pilot deliver per-stratum crab-vs-total classification counts, and at what stratum resolution? Everything in the code is ready for that answer the moment it arrives.
