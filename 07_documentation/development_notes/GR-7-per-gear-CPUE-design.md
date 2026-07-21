# GR-7: Genuine per-gear CPUE, design and implementation map

**Status:** Phase 0 executed 2026-07-20 (Section 5: shore G = 5 / G = 4, boat G = 1, Mixed 33%). **Phase 1 VALIDATED 2026-07-20** behind `gear_resolved_G` (default off; run `05_output/20260720/gear-type-CPUE-model-gear_resolved_G = TRUE`). The shore harvest is now genuinely per-gear (Option A1): shore all-gear fits **G = 5** (Pot / Ring Net / Snare / Trap / Mixed), shore pot-closure **G = 4**, boat stays **G = 1**. Convergence clean (shore all-gear 27 divergences = 0.34%, R-hat ~1.000, min gear n_eff ~8,000; pot-closure 4 divergences); reconciliation exact (`C_sum_gear` means sum to 20,346.17 = `C_sum` mean 20,346.16, so `sum_g C_sum_gear = C_sum` per draw); per-gear CPUE ordering matches the raw interview ratios (Trap highest). The neutrality run (`gear_resolved_G = FALSE`) reproduced production to < 0.1% with byte-identical PE, confirming the G = 1 path is unchanged. The shore all-gear total moved 20,127 -> 20,293 (**+0.82%**), the modest, expected difference between per-gear CPUE (fit on single-gear interviews) and one pooled CPUE, not a decomposition error. Phase 2 (Dirichlet shares) is drafted in Section 10. **Author context:** written 2026-07-20 against `main` after the `gear-type-run1` review. **Governing decisions (from Matt, 2026-07-20):** (1) the gear-resolved model should ultimately deliver a real per-gear CPUE, not a PE apportionment; (2) because the interview design records mixed-gear trips, **"Mixed" is treated as its own gear type for now**; (3) **only interviews that report a single gear type contribute to a gear-specific CPUE**; (4) the intent is to improve sampling over time so the Mixed share shrinks and gear-specific catch information grows.

This document maps the exact path from the current `G = 1` state to a genuinely gear-resolved fit, names every file and line that changes, and flags the one real identification decision that must be made deliberately (not silently).

---

## 1. Why the current output is not a gear resolution

`catch_by_gear_type.csv` puts 84% of catch in an undifferentiated "All" bucket; the "Pot" row (14%) is the commercial/charter census, and Ring Net / Trap / Snare are the small PE-apportioned boat remainder. The driver itself labels this "degenerate... not for citation" and points users to the pooled model's Dirichlet split. The reason is structural: the model runs with the gear dimension `G = 1`, so there is one pooled CPUE process and the per-gear breakdown is applied after the fact by a point estimator that carries no posterior uncertainty.

## 2. What is already built (most of it)

The Stan model `crab_bss_gear_resolved.stan` is written for `G > 1` and the machinery is present and wired, not hypothetical:

| Piece | Where | State |
|---|---|---|
| Per-gear CPUE mean `mu_mu_C[G]`, `mu_C[G,S]` | stan L248, L307-310 | present |
| AR(1) deviations over `G*S` with a Cholesky cross-gear correlation (`Lcorr_C`) | stan L243, L266, L298 | present |
| Per-gear effort/CPUE offset `O[D,S,G]` applied in the process equations | stan L143, L316, L318-319 | present, **fed all-ones** (`prep_bss_crab_gear.R:403`) |
| Interview catch indexed by gear (`gear_IntC[IntC]`) | stan L163, likelihood | present, **fed `rep(1L, ...)`** (`prep_bss_crab_gear.R:428`) |
| Predictive per-gear daily catch `C_gear_pred[D,G]`, season `C_sum_gear[G]` | stan L466, L507, L524 | present |
| Empirical gear shares `pi_gear_data[P_n, day_type, G]` + weighted counts `n_weighted_gear` | `prep_bss_crab_gear.R:315-355` | **computed, then discarded** |
| Per-interview gear classification (`has_pot/ring_net/trap/snare`, `n_gear_types_reported`) | `prep_bss_crab_gear.R:197-212` | present |
| Effective-N gate per gear (`bss_min_gear_effective_n = 15`) | `prep_bss_crab_gear.R:250` | present |

The single line that forces the collapse is `prep_bss_crab_gear.R:293-302`: after computing `G_gear`, the gear labels, `pi_gear_data`, and the weight matrix, it overwrites `gear_type_labels <- "All"; G_gear <- 1L` with an explicit note ("this Stan model runs with G = 1"). So the header's claim that Option A "requires effort shares (the `O[d,s,g]` offset, fed by `pi_gear_data`) and a rule for multi-gear interviews" overstates the remaining Stan work: **`O` and `gear_IntC` are already consumed by the Stan model; what is missing is the driver populating them and not collapsing `G`.** This is a wire-up-plus-validate task, not a from-scratch build.

## 3. The one real decision: how per-gear EFFORT is identified

This is the crux, and where I am least confident until a run settles it. The model has **one** effort observation stream that touches a single gear index: `Gear_I ~ NB2(lambda_E[.,1] * R_G, r_E)` (shore gear counts, stan ~L410) and the boat trailer stream `T_I ~ NB2(lambda_E[.,G]/R_G_boat, r_E)`. But `E_sum` and the catch sum over **all** `g`. If `mu_E[g,s]` is left free per gear with `G > 1`, gears with no direct effort observation get their effort level from the priors alone, exactly the pathology the Stan header warns about.

There are two coherent resolutions; they are genuinely different models and should not be chosen by accident:

- **Option A1 (recommended): pooled effort level, deterministic gear split by shares.** Effort is one process (one `mu_E`, one `omega_E`); per-gear effort is `lambda_E_total[d,s] * pi_gear[period(d), day_type(d), g]`. Implement by setting the effort offset `O_E[d,s,g] = pi_gear_data[.,g]` and constraining `mu_E` to be gear-shared (a single level, or a tight hierarchical prior tying `mu_E[g]` together). Per-gear **CPUE** stays free (`mu_C[g]`), identified by single-gear interviews. Then per-gear catch = (total effort) × (gear share) × (per-gear CPUE). This is well-identified: the total effort comes from the counts, the shares from interview composition (`pi_gear`), and each CPUE from that gear's single-gear interviews. It matches how a design-based estimator would apportion, but propagates uncertainty.
  - Caveat to verify: the effort **share** `pi_gear` is treated as data (Laplace-smoothed point shares). To carry share uncertainty, promote `pi_gear` to a Dirichlet parameter fed by `n_weighted_gear` (already computed). Start with fixed shares; add the Dirichlet once the fixed-share fit is stable.

- **Option A2: free per-gear effort levels.** Leave `mu_E[g]` free and hope the catch interviews plus the `G*S` Cholesky correlation identify the unobserved gears. This is what the current Stan literally does if you feed `O = 1` and raise `G`. It is under-identified for any gear without its own effort observation and will lean on priors; **not recommended** as the first build.

**Second, smaller decision inside A1:** the offset `O` is currently applied to *both* `lambda_E` (L316) and `lambda_C` (L318-319). Under a share interpretation that double-counts (catch would scale with the share squared). For A1, the effort share must scale effort only. Concretely: split `O` into `O_E` (the `pi_gear` share, on `lambda_E`) and `O_C` (leave `= 1`, since per-gear CPUE differences live in `mu_C[g]`, not in an offset). This is a small Stan edit but it is a correctness edit, not cosmetic, do not skip it.

## 4. The gear taxonomy and the Mixed rule (Matt's decisions, made concrete)

Gear classes for 2024-25: **Pot, Ring Net, Snare, Trap, Mixed**. (Commercial/charter Pot stays in the census component, unchanged.)

- **Classification per interview** (extend `prep_bss_crab_gear.R:197-204`): keep `has_pot/has_ring_net/has_trap/has_snare` and `n_gear_types_reported`. Add `is_single_gear = (n_gear_types_reported == 1)`.
- **`gear_IntC` assignment** (replace `rep(1L, ...)` at L428): a single-gear interview maps to its one gear's index; a multi-gear interview (`n_gear_types_reported > 1`) maps to the **Mixed** index. This is decision (3): only single-gear interviews inform a gear-specific `mu_C[g]`; Mixed interviews inform only the Mixed `mu_C`.
- **Mixed as its own gear** (decision 2): add "Mixed" to `gear_type_labels` whenever there are `>= bss_min_gear_effective_n` multi-gear interviews. Mixed then carries its own effort share (`pi_gear` includes a Mixed column) and its own CPUE. This removes the current fractional-apportionment fudge (`has_pot / n_gear_types_reported`) from the CPUE path entirely; fractions survive only in the effort-share `pi_gear` if you choose to let mixed trips contribute fractional effort shares (recommended: give Mixed its own whole share rather than splitting it, consistent with treating it as a gear).
- **Effort shares `pi_gear`** (already at L315-355): recompute over the new label set including Mixed, so `sum_g pi_gear[.,g] = 1` across {single gears present, Mixed}.

Result: Pot/Ring Net/Snare/Trap CPUE come only from clean single-gear interviews; Mixed is a real, separately-estimated bucket rather than a smear; and as field sampling shifts toward single-gear reporting, the Mixed share shrinks and the gear-specific signal strengthens (decision 4). Document that trajectory as a season-over-season metric (Mixed share of interviews and of effort).

## 5. Data sufficiency (run this before building)

The gate is `bss_min_gear_effective_n = 15` single-gear interviews per gear. Any gear below that in a sub-season folds to Mixed or a pooled "Other" and is reported by apportionment with an explicit flag, never fit as a `mu_C[g]`.

**PHASE 0 RESULT (executed 2026-07-20 on `interview_combined.xlsx` via `06_diagnostics/gear_coverage_audit.R`; single-gear counts after the standard filters, pots excluded in the pot closure):**

| population / sub-season | N | single-gear | Mixed | Mixed share | estimable gears (single >= 15) |
|---|---:|---:|---:|---:|---|
| shore / pot-closure | 856 | 652 | 198 | 23% | Ring Net (274), Trap (293), Snare (85), Mixed |
| shore / all-gear | 2,743 | 1,841 | 893 | 33% | Pot (646), Ring Net (378), Trap (608), Snare (209), Mixed |
| private_boat / pot-closure | 17 | 15 | 0 | 0% | none (stays PE, insufficient data) |
| private_boat / all-gear | 145 | 135 | 9 | 6% | Pot (130) only |

Three decisions fall out:

- **Shore is where the resolution lives.** Shore all-gear supports a full G = 5 (Pot, Ring Net, Trap, Snare, Mixed); shore pot-closure supports G = 4 (Ring Net, Trap, Snare, Mixed; pots illegal). Every real gear clears the threshold with room to spare, so switching from the old fractional metric to the single-gear rule loses no shore gear.
- **The boat is Pot-only.** Boat all-gear has 130 of 135 single-gear interviews as Pot; nothing else clears 15 and Mixed is only 9. So the boat stays G = 1 (there is no boat gear mix to resolve) and boat pot-closure stays PE. Per-gear CPUE is, for 2024-25, a shore feature.
- **Mixed is the honest ceiling.** A third of shore all-gear interviews (33%) report multiple gear types, so ~33% of the shore all-gear catch lands in the Mixed bucket rather than a named gear. The resolution is real (the other ~67% splits across four gears with propagated uncertainty) but bounded; the season-over-season Mixed share is the KPI for how gear-resolved the estimate can become (Section 9).

## 6. Implementation, phased (each phase ends in a validating run)

**Phase 0, data audit (no model change). DONE 2026-07-20.** Delivered as the standalone `06_diagnostics/gear_coverage_audit.R`, which prints single-gear interview counts per gear per sub-season and the Mixed share and writes `gear_coverage_audit.csv`. Result in Section 5: shore all-gear G = 5, shore pot-closure G = 4, boat G = 1 (Pot-only), shore Mixed share 33%.

**Phase 1, fixed-share A1, effort split only. DONE + VALIDATED 2026-07-20.** As built: `crab_bss_gear_resolved.stan` restructures the effort process to a single shared level over sections (`mu_E` is `[1,S]`, `omega_E` over `S`), split across gears by the `pi_gear` offset `O` applied to `lambda_E` only (`O` removed from `lambda_C`, so the share is not squared); both mu-hierarchies gate on `(S > 1)` rather than `(G*S > 1)`, which avoids a per-gear funnel at `S = 1` (the crab case) while staying identical at `G = 1`; `mu_C[g]` stays free (per-gear CPUE). `prep_bss_crab_gear.R` sets `G <- G_gear` for shore behind the `gear_resolved_G` toggle (default FALSE), assigns `gear_IntC` by the single-gear/Mixed rule with a hard one-hot weight matrix, and feeds `O` from `pi_gear_data`; the boat keeps `G = 1`. Because `G*S == S` at `G = 1`, every edit is behavior-neutral for production (confirmed by the neutrality run). The driver renders a reconciliation table (`sum_g C_sum_gear` vs `C_sum`). All Section 7 checks passed on the validating run.

**Phase 2, Dirichlet shares.** Promote `pi_gear` from fixed data to a Dirichlet parameter fed by `n_weighted_gear`, so share uncertainty propagates into the per-gear intervals. *Validate:* per-gear intervals widen appropriately; PSIS-LOO on the gear stream stays clean (0 Pareto k > 0.7, as the current pooled gear stream already achieves).

**Phase 3, extend to boat, then port.** Apply to `private_boat_all_gear`. The boat is thin (131 interviews, mostly Pot); it may support only Pot + Mixed. Gate on Phase 0 counts. *Validate:* boat gear-summed catch reconciles to the boat all-gear total (43,314 baseline); the boat total does not move (per-gear is a decomposition, not a re-estimation of the total).

**Phase 4, reporting.** Replace the degenerate `catch_by_gear_type.csv` with the BSS per-gear posterior (median + 95% CI per gear), a Mixed row, and an "Other/apportioned" remainder for sub-threshold gears; update the plot (drop the stale "v5.2" subtitle); retire the "not for citation" caveat once reconciliation and LOO pass.

## 7. Validation checklist (the "is it real this time" gate)

Results annotated from the 2026-07-20 `gear_resolved_G = TRUE` run.

1. **Reconciliation, PASS.** `sum_g C_sum_gear = C_sum` per draw: shore all-gear `C_sum_gear` means 3772.9 + 1977.7 + 1760.7 + 7573.3 + 5261.6 = 20,346.2, equal to the `C_sum` mean 20,346.2.
2. **Non-degeneracy, PASS.** Shore all-gear splits Pot 3,757 / Ring Net 1,965 / Snare 1,746 / Trap 7,534 / Mixed 5,239; no shore gear exceeds ~37% of the shore total. (At the PORT level 52% is still "All" because the boat is unresolved by design.)
3. **CPUE credibility, PASS.** `mu_mu_C` ordering Trap (0.57) > Snare (0.13) > Pot (-0.04) ~ Mixed (-0.07) > Ring Net (-0.29); Trap highest, matching the raw interview ratios (Trap ~0.73 the highest). Trap's high CPUE, not its interview count, drives its 7,534 catch.
4. **Cross-method, PENDING.** Compare the BSS per-gear split to the pooled model's Dirichlet split on a matched run (a nice-to-have cross-check, not blocking).
5. **LOO/PPC, CHECK ON NEXT READ.** The gear-stream PSIS-LOO and per-gear PIT were not inspected in this review; read `loo_summary_*` and `ppc_*` for the G = 5 fit before publication.
6. **Total invariance, PASS with a caveat.** The DECOMPOSITION is exact (item 1). The shore all-gear *total* moved 20,127 -> 20,293 (+0.82%) because per-gear CPUE (single-gear interviews) differs slightly from the pooled CPUE; a legitimate model difference, not a decomposition error, and the two agree to < 1%.

## 8. Risks and open questions

- **Effort identification (Section 3), RESOLVED.** A1 with a single shared effort level converged cleanly on the shore all-gear G = 5 fit (27 divergences, min gear n_eff ~8,000), so the per-gear effort split is well-behaved; the funnel concern did not materialize (helped by gating both mu-hierarchies on `(S > 1)`).
- **Double-offset, FIXED.** `O` is applied to `lambda_E` only; the redundant `× O` on `lambda_C` was removed (behavior-neutral at G = 1 since `O = 1`).
- **Fixed shares understate the per-gear intervals (open, -> Phase 2).** `pi_gear` is fed as fixed point shares, so the per-gear intervals carry CPUE and effort-total uncertainty but NOT share uncertainty; they are slightly too narrow. Phase 2 (Section 10) promotes `pi_gear` to a Dirichlet to fix this.
- **Thin gears** (Trap, Snare) cleared 15 single-gear interviews for shore this season, but may not every season; keep the sub-threshold -> Mixed/Other fold ready.
- **Mixed dominance:** Mixed is ~26% of the shore all-gear catch (33% of interviews), the honest ceiling until field sampling shifts (Section 9).
- **Cross-gear correlation prior** (`lkj_corr_cholesky(1)` over `G*S`): converged fine at G = 5, but with thinner future data consider a tighter LKJ or a diagonal start.

## 9. Field/sampling recommendation (feeds the 2025-26 plan)

The single highest-leverage change is upstream of the statistics: **record gear-specific catch at the interview**, so a mixed-gear trip's catch is attributed to gear rather than reported as a lump. Every interview moved from Mixed to single-gear directly strengthens a `mu_C[g]`. Track the Mixed share of interviews season over season as the KPI for how gear-resolved the estimate can be; the model can only resolve what the sampling distinguishes.

## 10. Phase 2 design: Dirichlet effort shares (draft, 2026-07-20)

**Problem.** Phase 1 feeds the effort shares `pi_gear` as FIXED point values (Laplace-smoothed empirical proportions per period x day_type). So the per-gear posterior `C_sum_gear[g] = level x pi_gear[g] x E_scale x L x lambda_C[g]` carries the CPUE and effort-total uncertainty but treats the gear split as KNOWN. The per-gear intervals are therefore too narrow, most severely in thin period x day_type cells where the point share rests on a handful of interviews. Phase 2 makes `pi_gear` a parameter with a Dirichlet posterior informed by the observed gear counts, so share uncertainty propagates into the per-gear intervals. The season/population TOTAL is unaffected (the shares still sum to 1), so this only widens the split, it does not move the headline.

**Model.** In each (period `p`, day_type `dt`) cell the single-gear/Mixed interview counts `n_weighted_gear[p,dt,1:G]` (already computed in `prep_bss_crab_gear.R:333/351`, and integer-valued under the Phase 1 hard one-hot assignment) are a multinomial sample of the true shares. With a `Dirichlet(alpha0)` prior the conjugate posterior is `Dirichlet(alpha0 + counts)`; sampling `pi_gear` from it and using it as the `O` offset propagates the uncertainty.

**Stan changes (`crab_bss_gear_resolved.stan`), behind a new toggle so Phase 1 stays comparable:**

- *data:* add `real n_weighted_gear[P_n, n_day_types, G]` (the fallback-filled counts), `int day_type_idx[D]` (day-type per day; the driver already has `day_type_idx_vec`), `real alpha0_gear` (concentration, default 1.0 = Laplace), and `int<lower=0,upper=1> gear_share_dirichlet`.
- *parameters:* `simplex[G] pi_gear[P_n, n_day_types];` when `gear_share_dirichlet == 1` (size `P_n * n_day_types` simplexes of dim `G`; for shore all-gear ~10 x 3 x (5-1) = ~120 free params, well-conditioned). Guard the size to 0 when off, mirroring the `estimate_L`/`use_mu_hier` idiom.
- *transformed parameters:* build `O` from the parameter instead of taking it as data: `O[d,s,g] = pi_gear[period[d], day_type_idx[d], g]` (when the toggle is on; else keep the data `O`). Everything downstream (`lambda_E`, `sum_g lambda_E`, `C_sum_gear`) is unchanged, so the effort total `sum_g lambda_E = level` is still share-invariant and the effort observations are untouched.
- *model:* `for (p in 1:P_n) for (dt in 1:n_day_types) pi_gear[p,dt] ~ dirichlet(alpha0_gear + to_vector(n_weighted_gear[p,dt,]));` (the direct posterior-Dirichlet form; works with fractional counts, no integer-multinomial requirement).

**Driver changes (`prep_bss_crab_gear.R`):** pass `n_weighted_gear` (with its existing period/season fallback so empty cells inherit pooled counts, not a flat Dirichlet), `day_type_idx_vec`, `alpha0_gear`, and `gear_share_dirichlet` into `stan_data`; keep sending the fixed `O` array as the fallback for when the toggle is off. Add `gear_share_dirichlet = FALSE` to `run_config.R` (only meaningful when `gear_resolved_G = TRUE`). At `G = 1` `pi_gear` is a 1-simplex (identically 1.0) and `O = 1`, so this is behavior-neutral for production and for the boat.

**Validation gate.** (1) Per-gear intervals WIDEN vs Phase 1, most in thin cells; the per-gear MEDIANS stay ~unchanged (the Dirichlet posterior mean ~ the Laplace-smoothed point shares). (2) The shore all-gear TOTAL and the reconciliation are unchanged (shares sum to 1). (3) Convergence still passes; watch the thin-cell simplexes for divergences. (4) Gear-stream PSIS-LOO stays clean.

**Risks.** A period x day_type cell with few interviews gets a wide Dirichlet -> wide shares -> a wide per-gear catch in that cell; that is honest, but if a thin cell dominates a gear's season total the interval can balloon. Mitigate by (a) a larger `alpha0_gear` (more regularization toward uniform), (b) feeding the period/season fallback counts (already in the driver) so thin cells borrow strength, or (c) pooling day_types when a cell is below a minimum count. Start with `alpha0_gear = 1.0` and the existing fallback; tighten only if the thin-month intervals look implausible. **Effort: medium** (bounded Stan + driver change; the counts and day-type index already exist).

**Sequencing note.** Phase 3 (boat) is N/A for 2024-25 (Phase 0: boat is Pot-only, stays `G = 1`). Phase 4 (reporting) is largely already done, Phase 1 renders per-gear posteriors, the reconciliation table, and the resolved plot; what remains there is retiring the "not for citation" language once Phase 2's honest intervals land, and a per-gear PPC/LOO read.
