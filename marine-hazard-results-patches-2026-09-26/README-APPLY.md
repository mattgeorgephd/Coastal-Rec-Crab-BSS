# Marine hazard ladder: results review and next steps (2026-09-26)

Two patches on top of `8b3f661` (the tip of `OSP-boat-count-incorporation` after your run commit). Nothing pushed.

```sh
git checkout OSP-boat-count-incorporation        # at 8b3f661
git am --3way 000*.patch
```

Verified here: both apply cleanly on a clean branch at `8b3f661` and the tree is identical to the one the patches were cut from.

## What the run said (full reading: `VALIDATION_CAMPAIGN.md` Section 1x)

| term | identified? | paired effort elpd (sampled days) | estimate | rule's verdict |
|---|---|---|---|---|
| shore SCA | no (+0.04 [-0.14, +0.21]); doubles shore divergences | -0.6 nats (-1.3 SE) | +0.2% shore | do not adopt |
| boat SCA | yes, 8 SE (rate ratio 0.31) | **+7.6 nats at 1.46 SE** (rule asks > 2) | boat +3.8%, port +1.6% | identified and harmless, no sampled-day gain |
| bar tick beyond the archive | by a hair (0.70 [0.49, 1.00]) | +0.3 nats (0.20 SE) | | redundant with the archive |

M1 (`off`) is bit-identical to R4 on every parameter row and component. `marine_hazard_mode` stays `"off"`. Nothing is adopted under the pre-committed rule; Section 1x.5 sets out the three defensible readings.

## What the patches contain

**0001, B38 + B39 (code).**
- `bss_stan_fit()` now rebuilds each fit's draw permutation from the Stan seed after the fit returns, and both drivers seed the census draw, so the port total is reproducible between identical renders. The M1-vs-R4 port difference (94,497 vs 94,376, 0.128%) was this jitter, not the model. Consequence: the next render of the R4 configuration will differ from 94,376 once, by about 0.1%, with every component identical; the register and the status box say so.
- `03_R_functions/bss_block_cv.R` + `06_diagnostics/run_marine_block_cv_2026-09-26.R`: the leave-one-week-out block cross-validation D31 asked for, by PSIS on the saved `ppc_draws_<fit>.rds`, no refit, with the rebuild checked against each rung's own pointwise LOO before anything is scored. Both drivers also write `loo_block_<stream>_<fit>.csv` per fit from now on.

**0002, the review (docs + runner fixes).** Section 1x; register A30/B35/D31 updates, two C rows (the runner's 0.05% port tolerance failed a correct baseline; the M5 display row), B38, B39; status document; design note Section 10; method document 13 and 21a; runner fixes (`DRY_RUN <- TRUE` restored, tolerance 0.3%, M5 row, `CODE_EQUIVALENT_MH`); harness section 77 (evidence pinned to the committed CSVs; block CV against exact conjugate answers).

## What to run next (minutes, on the machine that rendered the rungs)

```sh
Rscript 06_diagnostics/run_marine_block_cv_2026-09-26.R   # dry run: inventory + reconstruction check (R1)
# then set DRY_RUN <- FALSE in the file and source it again
```

It needs the git-ignored `ppc_draws_<fit>.rds` in the five rung folders (written by the run because `save_ppc_draws = TRUE`). It writes `05_output/marine_hazard_2026-09-26_blockcv{,_pairs,_verdicts}.csv` and `loo_block_*.csv` into each rung folder. Read it with its R2 rule in mind: weeks with Pareto k > 0.7 are excluded and counted; if fewer than 70% of weeks are reliable the comparison is REVIEW and the honest next step is refitting those weeks, not a verdict. Expect the shore fits (weekly AR) to lose many weeks; the boat (monthly AR) is where the answer should be readable.

If the boat trailer AND OSP streams clear +2 paired SE on held-out weeks, the adoption edit is `marine_hazard_mode = "manual"`, `marine_hazard_manual_boat = "nws_sca_any"` (the bar column adds nothing beyond the archive). D30 (a boat term describes all private boats) stands either way.

Optional, one minute: re-source `run_marine_hazard_batch_2026-09-25.R` with `DRY_RUN <- FALSE` after applying: every rung RESUMEs from its folder and the three ladder CSVs are rewritten with the corrected port row (PASS) and the M5 display fix. Restore `DRY_RUN <- TRUE` before committing.

## Harness

`Rscript 06_diagnostics/test_improvements_2026-08-25.R`: 1,100 assertions pass here; the 12 failures in this sandbox are the environmental ones (most `05_output` date folders absent from its sparse checkout) that fail identically on the untouched tip. In a full checkout expect 0.
