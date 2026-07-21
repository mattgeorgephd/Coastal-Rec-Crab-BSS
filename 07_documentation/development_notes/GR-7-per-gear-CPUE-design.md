# GR-7: Genuine per-gear CPUE, design and implementation map

**Status:** Phase 0 executed 2026-07-20 (Section 5: shore G = 5 / G = 4, boat G = 1, Mixed 33%). **Phase 1 CODED 2026-07-20** behind `gear_resolved_G` (default off): the Stan effort process is now a single shared level split by the `pi_gear` share offset `O` with per-gear CPUE (Option A1), and `prep_bss_crab_gear.R` un-collapses `G` for shore, assigns `gear_IntC` by the single-gear/Mixed rule, and feeds `O` from `pi_gear_data`. NOT yet run: R/Stan could not be compiled in the authoring environment, so a `gear_resolved_G = TRUE` shore render is the validation gate (Section 7). At G = 1 every edit is behavior-neutral (because `G*S == S`), so production is unchanged. **Author context:** written 2026-07-20 against `main` after the `gear-type-run1` review. **Governing decisions (from Matt, 2026-07-20):** (1) the gear-resolved model should ultimately deliver a real per-gear CPUE, not a PE apportionment; (2) because the interview design records mixed-gear trips, **"Mixed" is treated as its own gear type for now**; (3) **only interviews that report a single gear type contribute to a gear-specific CPUE**; (4) the intent is to improve sampling over time so the Mixed share shrinks and gear-specific catch information grows.

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

**Phase 1, fixed-share A1, effort split only.** In `prep_bss_crab_gear.R`: stop the L293-302 collapse (guard it behind a `params$gear_resolved_G` toggle, default FALSE so production is unchanged); populate `O_E[d,s,g] = pi_gear_data[period(d), day_type(d), g]`; set `gear_IntC` per the Mixed rule; tie `mu_E[g]` (shared level). In `crab_bss_gear_resolved.stan`: split `O` into `O_E` (on `lambda_E`) and `O_C = 1` (on `lambda_C`); keep `mu_C[g]` free. Fit shore all-gear only (best-sampled). *Validate:* `sum_g C_sum_gear == C_sum` (gear catches reconcile to the all-gear total to within Monte-Carlo error); per-gear CPUE ordering is physically sensible; convergence gate passes; compare the gear split to the pooled Dirichlet split (they should agree in central tendency, with the BSS carrying wider, honest intervals).

**Phase 2, Dirichlet shares.** Promote `pi_gear` from fixed data to a Dirichlet parameter fed by `n_weighted_gear`, so share uncertainty propagates into the per-gear intervals. *Validate:* per-gear intervals widen appropriately; PSIS-LOO on the gear stream stays clean (0 Pareto k > 0.7, as the current pooled gear stream already achieves).

**Phase 3, extend to boat, then port.** Apply to `private_boat_all_gear`. The boat is thin (131 interviews, mostly Pot); it may support only Pot + Mixed. Gate on Phase 0 counts. *Validate:* boat gear-summed catch reconciles to the boat all-gear total (43,314 baseline); the boat total does not move (per-gear is a decomposition, not a re-estimation of the total).

**Phase 4, reporting.** Replace the degenerate `catch_by_gear_type.csv` with the BSS per-gear posterior (median + 95% CI per gear), a Mixed row, and an "Other/apportioned" remainder for sub-threshold gears; update the plot (drop the stale "v5.2" subtitle); retire the "not for citation" caveat once reconciliation and LOO pass.

## 7. Validation checklist (the "is it real this time" gate)

1. **Reconciliation:** `sum_g C_sum_gear` equals the all-gear `C_sum` (same posterior), so the decomposition conserves the total.
2. **Non-degeneracy:** no single gear holds > ~70% by construction; the split reflects `pi_gear` × per-gear CPUE, not an "All" bucket.
3. **CPUE credibility:** per-gear `mu_C[g]` orderings match the interview ratio-of-sums by gear (`sensitivity_incomplete_by_gear.csv` already shows Pot ~0.30, Ring Net ~0.38, Trap ~0.73, Snare ~0.44 crab/deployment on complete trips, use as the sniff test).
4. **Cross-method:** BSS per-gear split brackets the pooled Dirichlet split.
5. **LOO/PPC:** gear stream Pareto k clean; per-gear PIT ~0.5.
6. **Total invariance:** the port total and each population total are unchanged vs the `G = 1` run (a decomposition must not move the sum).

## 8. Risks and open questions

- **Effort identification (Section 3)** is the main risk; A1 with a shared effort level is the defensible default, but confirm on the shore all-gear fit before trusting per-gear effort.
- **Thin gears** (Trap, Snare) may never clear 15 single-gear interviews; be willing to publish only Pot, Ring Net, Mixed, and an apportioned remainder rather than force a fit.
- **Mixed dominance:** if Mixed is a large share, the "resolution" is limited until field sampling shifts; report the Mixed share as the honest ceiling on how gear-resolved the season actually is.
- **Cross-gear correlation prior** (`lkj_corr_cholesky(1)` over `G*S`): with few gears and thin data the correlation is weakly identified; consider a tighter LKJ or a diagonal start.
- **Double-offset bug** (`O` on both `lambda_E` and `lambda_C`): must be fixed for the share interpretation (Section 3), or per-gear catch is mis-scaled.

## 9. Field/sampling recommendation (feeds the 2025-26 plan)

The single highest-leverage change is upstream of the statistics: **record gear-specific catch at the interview**, so a mixed-gear trip's catch is attributed to gear rather than reported as a lump. Every interview moved from Mixed to single-gear directly strengthens a `mu_C[g]`. Track the Mixed share of interviews season over season as the KPI for how gear-resolved the estimate can be; the model can only resolve what the sampling distinguishes.
