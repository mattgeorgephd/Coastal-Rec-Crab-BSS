# Block cross-validation: results read, next step built (2026-09-27)

One patch on top of `c7e8cd5` (the tip of `OSP-boat-count-incorporation` after your block-CV results commit). Nothing pushed.

```sh
git checkout OSP-boat-count-incorporation        # at c7e8cd5
git am --3way 0001-*.patch
Rscript 06_diagnostics/test_improvements_2026-08-25.R   # expect 0 failures in a full checkout (1,149 pass here; 12 sandbox-environmental)
```

Verified here: it applies cleanly on a clean worktree at `c7e8cd5` and the resulting tree is identical to the one it was cut from. Do not commit this tarball or this README into the repository (the `.gitignore` in the patch now excludes both; `c7e8cd5` had swept the previous ones in, and the patch removes them).

## What the run said (full reading: `VALIDATION_CAMPAIGN.md` Section 1y)

| | trailer (34 of 37 weeks) | OSP (25 of 29 weeks) | both streams (38 of 42 weeks) |
|---|---|---|---|
| boat SCA, M2 vs M1 | +7.9 nats, 1.20 SE (22 of 34 weeks positive) | **+18.9, 3.90 SE** (21 of 25 positive) | **+26.8, 2.86 SE** |
| bar beyond the archive, M4 vs M2 | -2.0 (-1.0 SE) | -0.4 | -2.4 |

The reconstruction check held on all 20 fits. The shore fits were not evaluable at weekly AR (8 to 16 of 38 weeks reliable), as the runner's header predicted; that does not matter for the decision (the shore term is unidentified). M5 = M4.

**The clause I pre-committed (each boat stream separately above +2 SE) was not met, and it was under-powered for the trailer stream:** its paired SE is 6.6 nats, so it needed +13, and the same effect scored 3.9 SE on the OSP counts and 1.6 SE on the trailer counts over the same weeks. That is recorded as a defect in the rule (C row), and the next run is judged on the joint row (R3').

**The winter.** 46% of the boat change under SCA falls in the weeks no held-out test could see (December to 2 March; the OSP counts begin in March; the trailer's 13 winter weeks read -0.6 +- 3.9). A desk screen of every sampled launch day over four seasons (794 days; `06_diagnostics/desk_sca_season_split_2026-09-27.R`, seconds) puts the winter SCA effect at **0.43 [0.29, 0.62]** against **0.22 [0.18, 0.27]** for the rest, the interaction at 3.0 SE (negative binomial) or 1.4 SE (quasi-Poisson), the same sign in every season. The season-constant term (0.31 in the fit) therefore probably over-corrects the winter by about a factor of two on the log scale. That is D32.

## What the patch contains

**B40 (fixes and the joint score).** The runner FAILED a comparison at "-2.31 SE" whose difference was 2e-14 nats (two copies of the same shore fit): identical tables now read INFO (R7). `bss_block_cv_fit()` adds the exact **joint** table (all held-out effort observations of a week together), which the runner judges the boat on (R3'). Per-week differences are written to `marine_hazard_2026-09-26_blockcv_weeks.csv`. `write_loo_diagnostics()` writes the OSP stream's pointwise LOO from now on, and the block-CV runner checks it exactly where it exists. `DRY_RUN <- TRUE` restored; the tarball and README removed from the tree and ignored.

**B41 (the season split).** `bss_marine_hazard_covariates.R` offers `nws_sca_any_winter` (December to February, `marine_hazard_winter_months` in `run_config.R`) and `nws_sca_any_rest`, one NWS definition in two columns, not in the shipped candidate lists (M2 to M5 resolve identically before and after, measured on the real inputs). The ladder runner gains rung **M6** (boat = the pair, no shore term) with **rule 9** stated before the render in three ordered branches; per-stage declared keys keep the M1 to M5 digests exactly as their folders recorded (checked), so RESUME reads them back and only M6 renders. The block-CV runner gains M6 and the pair M6 vs M2 on the boat fits. The desk screen and its table (`05_output/marine_hazard_2026-09-27_season_split_screen.csv`) are committed.

**Docs.** Campaign Section 1y; register A30 / B39 / D31 updated, B40, B41, three C rows, D32; status document; design note Section 11; both READMEs; `run_config.R` 4.4b comment; `CLAUDE.md`; `PULL_REQUEST.md`. Harness section 78.

## What to run next (on the machine that rendered the rungs)

```sh
# 1. the ladder: M1 to M5 RESUME from their folders (a minute), M6 renders (about 3.5 h)
#    set DRY_RUN <- FALSE in 06_diagnostics/run_marine_hazard_batch_2026-09-25.R, source it, restore TRUE
# 2. the block CV: joint and INFO rows and the weeks file for M1 to M5, every M6 row (a minute)
#    set DRY_RUN <- FALSE in 06_diagnostics/run_marine_block_cv_2026-09-26.R, source it, restore TRUE
# 3. commit the M6 rung folder and the CSVs (not the .rds); the harness asserts DRY_RUN TRUE on both runners
```

Rule 9 says in advance what M6 can settle. The likely outcome is branch (a): the winter coefficient not identified (36 sampled winter advisory days of five trailers or fewer). M6 then delivers the identified rest coefficient, the first exact OSP reconstruction check on a rendered fit, and the recorded fact that this season cannot place the winter effect that the four-season record can. Branch (b) (both identified and different; the split earns its place if the catch stream is unmoved against M2 and the block CV's joint boat all-gear row for M6 vs M2 is not below -2 SE) and branch (c) (no detectable difference; the constant term stands) are stated in the runner's header and printed by its recommendation.

## The decision (Section 1y.5)

1. Keep `off`: no winter correction at all, when the four-season record says an advisory halves the winter launch count.
2. Adopt the constant term (`manual`, `nws_sca_any`; not the bar column): boat +3.8%, port +1.6%, with a winter correction that is probably about twice too strong on the log scale (D32), D30 unresolved.
3. Adopt `auto`: M4's fits, a redundant bar column, no gain over 2.
4. Render M6 first, then adopt the shape it supports; if the winter coefficient stays unidentified, take 2 with D32 as the standing caveat.

My input is 4. Nothing is adopted by the patch; `marine_hazard_mode` still ships `"off"`.
