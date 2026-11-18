> **ARCHIVED: MERGED 2026-11-18.** This was the description of
> [pull request #8](https://github.com/mattgeorgephd/Coastal-Rec-Crab-BSS/pull/8), which merged
> `OSP-boat-count-incorporation` into `main` (merge commit `27d68c7`, from head `fa0186b`; 59
> commits, opened 2026-09-29, merged 2026-11-18 by Matt). It is kept as the record of what the
> merge carried. Its body is a stack of dated status paragraphs, newest first, each written when
> the branch reached that state (the last names `0075c8f`, B76); every figure in it is the figure
> of its date, and the only current one is the box at the top of
> `07_documentation/development_notes/PIPELINE_STATUS.md`. The head also carried `fa0186b`
> (2026-11-02, `04_input_files/razor-clam-dig-dates-2021-2026.xlsx`), which the description below
> predates (CHANGE_REGISTER B77).

# Merge `OSP-boat-count-incorporation` into `main`, the second span: A31 to A41, B50 to B76

## Summary

**Status (2026-10-20, `0075c8f`): A41 ADOPTED AND RENDERED, D49 AND D19 CLOSED, on Matt's decisions of 2026-10-20 ("Adopt D49 option (iii) and close D19 keeping local"). A41 (`e7dea7c`, the rule written first in `407379a`): a fit below the shared-turnover floor keeps ONE tau_bar at its prior that no day observes (`shared_tau_refusal = "shared_prior"`; the fit's `osp_scale_is_tau` set to 0 so its OSP totals, if any, observe presence through `kappa_OSP`, its boat I/E stream dropped), where it fell back to per-day draws whose season sum sat 13% above the calibration centre with the turnover's uncertainty averaged away; `"per_day"` restores the former fallback bit for bit. No 2024-25 fit is below the floor: the Stan data of every 2024-25 fit on both tracks are identical before and after (dump), so THE AUTHORITATIVE RUN STAYS B70 (86,497). B76, the 2025-26 readiness render re-rendered under A41 (`05_output/20261020/`): every clause of the rule met on both tracks, tau_bar at its prior (3.503 [1.312, 9.338]), the boat all-gear effort 11,865 to 10,652 (-10.2%, the arithmetic said -10.8%) with its interval's upper-to-lower ratio 2.0 to 7.7, the boat catch 23,370 to 21,045 [7,703, 59,870], the port 47,480 to 45,233 [30,855, 83,978], cross-check -0.47% PASS, the presence process unmoved; the D50 boat gap +27.3% to +14.3% on effort as predicted. The 2024-25 cache render could not load B70's fits (saved under the pre-B71 key); a full re-sample of 2024-25 on the committed tree (`0075c8f`'s predecessor `fa4977b`, 08:10Z to 10:30Z) then reproduced B70 draw for draw: all eight saved-draw files byte-identical, the port 86,497 [72,869, 104,595] and the cross-check 86,285 (-0.25%, PASS) exactly (`05_output/20261020/A41-B70-reproduction/`), so A41's clause (1) is met by run. D19 closed keeping `pe_empty_stratum = "local"` on the B70 re-measurement. Nothing from 2025-26 is citable (D18). Recorded in CHANGE_REGISTER A41, B76, D49, D19, D50, campaign 1z.19, the status document, the guide, CLAUDE.md, the method document (limitations 4 and 16) and both histories; 1,658 harness assertions passing.**

Earlier (2026-10-19, `64af6b9`): B75 READ, the 2025-26 readiness render. The novel-season workflow ran end to end on both tracks from the guide's section 7.1 block on the container, three attempts for one set of fits (the first failed at the gear driver's last chunk on a B73 column; the second rendered clean; the third is the record after the pooled report's turnover table was fixed to show a fit whose shared turnover is refused), every pre-written check R1 to R7 met, the cross-check PASS at -0.24% (gear 47,365 against pooled 47,480; readings, not estimates: nothing from 2025-26 is citable, D18). The day-effect-off refit of the boat all-gear fit sits 0.42 posterior SD below the render's, inside the rule, with the turnover at its prior in both arms. Two findings for Matt: D49, a fit whose shared turnover is REFUSED (no OSP day) runs independent per-day prior draws whose season sum sits at the lognormal MEAN, 13% above the calibration centre at the shipped log-SD 0.5, with the turnover's uncertainty averaged away; three options in the register, (iii) recommended (keep one tau_bar under the refusal with no observation); no 2024-25 fit is affected. D50, the season's PE/BSS gaps (-22.5% shore all-gear, +28.9% boat all-gear), part explained by D49 and by the shore fit's turnover following its 16 in-window I/E days (2.30 against the derived 2.52 the PE applies), the rest to examine on the complete data. Three report gaps fixed (`f352787`, `cee229b`, `b03743a`). Recorded in CHANGE_REGISTER B75, D49, D50, D18 and two C rows, campaign 1z.18, the status document, the guide (section 7.1), CLAUDE.md, the method document (limitation 16) and both histories; 1,644 harness assertions passing.

Earlier (2026-10-19, `143e888`): B74 READ, the D39/D2 pricing batch, every gate passing. D39 (the OSP days with boats and no crab-only value) is worth about 100 crab on the port whichever way OSP answers: read as zeros, the boat all-gear fit +0.1% (35,491 against 35,462) and the pot-closure fit -9.5% (1,225 against 1,354), the port -0.1%. D2 (the day's actual sampling rate) is large: at n = every returning boat the boat all-gear fit falls to 26,492 [19,131, 37,779] (-25.3%) and the port to 77,494 (-10.4%), with the turnover and the day effects unmoved, so the whole difference is the crabbing share. OSP's actual daily sampling rates (or the number of private boats sampled per day) are now the most consequential open data item, ahead of D39 and of the remaining 2025-26 OSP days. OSP's observed crab-only share sits below the model's f(1 - c) interval in four of nine months on B70 and never above (the contacts hold f up, D11). Recorded in CHANGE_REGISTER B74, D2, D39, D11, campaign 1z.17, the status document, CLAUDE.md and the method document. B75, the 2025-26 readiness render (not citable, D18), is running on the container from `8a2ff70` with the guide's block pasted into run_config.R; its checks are CHANGE_REGISTER B75.

Earlier (2026-10-19, `8a2ff70`): WAITING FOR OSP, the work Matt chose meanwhile. B72 (`585e41a`): the OSP absent-day reader key `osp_absent_day` (D39 pre-work; "unclassified" ships, the reading since B64; "zero" reads a day with boats and no crab-only value as 0 crab-only of a schedule n; 17 such days on 2024-25) and D41 (a) to (d) (the schedule-minimum and absent-day facts in run_warnings.csv and the report; OSP's observed share beside the model's f(1 - c) by month in section 12.3, inside the model's interval in 5 of 9 months on B70, OSP below it in the other four; an ingest error stops the run where it happens; `pe_f_source`), D41 CLOSED. B73 (`591897a`): the PE has a catch SE (`pe_catch_variance()`: the stratum effort variance with a day-level ratio variance of the CPUE and the donor covariances; on 2024-25 959, 3,102, 432 and 6,864 for the four components, every B70 BSS median inside +/- 1.96 SE), `PE_catch_se`, `PE_se`, `Dungeness_SE` written, and a PE-reported component enters the port total as seeded draws with its SE instead of a constant; D22 part 2 CLOSED; B70 bit-identical. D19 re-measured against B70: the evidence now favours the shipped `local` CPUE fill (the boat all-gear PE +4.2% on the BSS, against +20.6% at `pooled` and -35.8% at `zero`); closure is Matt's call. B74, the D39/D2 pricing batch (the two boat fits under the "zero" reading and under n = every boat), and B75, the 2025-26 readiness render (not citable, D18), are running on the container under rules written first (CHANGE_REGISTER B74, B75; campaign 1z.17, 1z.18). First results: the "zero" reading moves the boat all-gear fit +0.1% (35,491 against 35,462) and the boat pot-closure fit -9.5% (1,225 against 1,354), both gates passing: D39 is worth about -100 crab on the port. Harness 1,641 passed, 0 failed.**

Earlier (2026-10-19, `a608e9d`): B70 ACCEPTED (Matt: "Accept B70"). THE AUTHORITATIVE RUN IS `05_output/20261019/pooled-CPUE-2024-25`, port 86,497 [72,869, 104,595], cross-check `05_output/20261019/gear-type-CPUE-model-2024-25` at 86,285 (-0.25%, PASS); A39 (the per-day effort effect, D45) and A40 (its identification floor) are the method of record; B67 (86,113 [72,407, 103,557]) is the run it replaced, +0.4% on the port, the half-width +1.8%.** The box, the register, CLAUDE.md, the method document, campaign 1z.16, both histories, the README run table, the guide and the `run_config.R` header moved with it, and the shipped `run_config.R` is the authoritative run's configuration again. Recorded with the move: the shared boat turnover `tau_bar` on the boat all-gear fit is 2.97 in the posterior mean on B70 (median 2.95 [2.24, 3.82]) against 3.16 on B67; the per-day deviations on the trailer counts take up part of what the turnover absorbed, and the boat total, which OSP pins, does not move (35,403 to 35,462). The BSS sits 2.5% below the PE (88,758). Harness: 1,622 passed, 0 failed.

Earlier (2026-10-19, `7ce754d`): A40, the identification floor on the per-day effort effect (Matt's option (a) for B69), BUILT, RULE WRITTEN, RENDERED (B70): every clause of the rule met; the port 86,497 [72,869, 104,595], cross-check 86,285 (-0.25%, PASS). It is the candidate authoritative run pending Matt's acceptance, so B67 (86,113 [72,407, 103,557]) stays the authoritative run and every document says so.** `effort_day_effect = TRUE` again, with `effort_day_effect_min_repeat_days = 5`: each fit switches the term on only where at least five of its sampled days carry two or more counts (2024-25: shore all-gear 83 and boat all-gear 45, on; boat pot-closure 1 and shore pot-closure 0, off, with a run warning, and `fit_data_summary.csv` records the switch and the count for every fit). The render (`3fbecc5`, clean tree, the B66 fit cache, R 4.3.3 and rstan 2.32.5, 62 + 72 min): (1) 8 of 8 fits pass the gate; (2) no stuck chain; (3) the cross-check -0.25%; (4) `sigma_uE_out` clean on the four fits where the term is on (shore all-gear 0.42 [0.36, 0.50], n_eff 3,456 and 2,342; boat all-gear 0.84 [0.70, 1.01], 3,175 and 3,661; R-hat at most 1.0012); (5) the four pot-closure fits reproduce B67's draw for draw and the four all-gear fits reproduce B69's. Against B67: the port +0.4% and the half-width +1.8%; the shore all-gear 30,804 to 31,130 (+1.1%), the boat all-gear 35,403 to 35,462 (+0.2%), the pot-closure components unchanged. **B71:** every fit of the render was a cache miss because `bss_fit_cache_key()` digested the Stan data list's top-level attributes (the preps' metadata; A40 added one), so the unchanged all-gear fits were re-sampled rather than loaded and reproduced B69 draw for draw; the key now ignores those attributes (`bss_fit_cache_strip()`), and both live batch runners re-declare their fingerprints. Harness: 1,617 passed, 0 failed. **Matt's decision (CHANGE_REGISTER A40, B70):** whether B70 becomes the authoritative run; accepting it moves the box, the register, CLAUDE.md, the method document, the README and the guide to 86,497 [72,869, 104,595].

Earlier (2026-10-18, `81a83f3`): A39 (the D45 per-day effort effect on both tracks, Matt's decision) BUILT, RULE WRITTEN, RENDERED (B69): clauses (1) to (3) met, clause (4) FAILED on one fit; `effort_day_effect` is back to FALSE per the rule, so the shipped file is B67's configuration again and B67 (86,113 [72,407, 103,557]) stays the authoritative run.** The render (`b1c1cf5`, clean tree, the B66 fit cache, relaunched once after a container restart): 8 of 8 fits pass the gate, no stuck chain, cross-check 86,405 against 86,456 (-0.06%, PASS); the port 86,456 [72,750, 104,679] against B67's 86,113 (+0.4%, the half-width +2.5%), the shore all-gear +1.1%, the boat all-gear +0.2%. Clause (4), `sigma_uE_out` sampled cleanly on all eight fits, failed on the pooled shore pot-closure fit (n_eff 301; R-hat 1.0057): that fit has one gear count per sampled day, so the day deviation and the count noise are one quantity and `sigma_uE` trades off with `sigma_r_E`; the gear track's same fit reads 413 and every other fit 906 to 3,661; the gate quantities are untouched. **Matt's decision (CHANGE_REGISTER B69):** (a) an identification floor, the term only on fits with at least five sampled days carrying repeat counts (the two all-gear fits; both pot-closure fits return to B67's program bit for bit; the re-render reuses the all-gear fits from the cache), under a re-stated rule, recommended; (b) accept B69 with a documented exception to clause (4); (c) more iterations on the pot-closure fits; (d) leave the term off. D24 CLOSED by decision (keep 2.477; the window-only price, -4.5% on the port, quoted as the sensitivity); D45 CLOSED by adoption. Earlier the same day: D24, D29 and D45 priced (1z.14), D27 and D42 closed, B68 built.

Earlier (2026-10-18, `9867793`): D24, D29 and D45 PRICED on the container (campaign Section 1z.14); D27 and D42 CLOSED; the per-day effort effect built as an experiment lever (B68, `effort_day_effect`, ships FALSE and bit-identical).** The twelve-arm batch (04:50Z to 08:25Z) read: **D24**, the window-only shore turnover centre 2.225 prices at -9.5% and -10.4% on the two shore fits, -3,882 on the port (82,232 against 86,114), the fitted turnover following the centre; **D29**, the boat all-gear AR span re-measured under the current model, monthly 35,403, biweekly 37,947, weekly 38,767 (+9.5%), all gate-passing with no divergences, the port 86,113 to 89,387 (+3.8%, against +6.6% before A33 to A36), no longer stale, the decision (monthly kept, span quoted) standing; **D45**, the shore all-gear interval half-width +7.5% (SD +8.9%) with the point +1.1%, the boat +4.0% and +0.2%, the port SD +2.5% (86,483 [72,870, 104,654] against 86,114 [72,377, 103,516]). The lever's first cut lacked the mean-one correction and read the shore all-gear effort 5.8% low; the second cut (`4ad0dde`) is what shipped, and both cuts' OFF paths reproduce the A38 refits draw for draw. Two decisions for Matt: D24's centre, and whether to adopt the D45 lever on both tracks (a method change under an A36-type rule) or carry the understatement as a stated limitation. The authoritative run is unchanged (B67, 86,113 [72,407, 103,557]).

Earlier (2026-10-18, `d7ee759`): **THE AUTHORITATIVE RUN IS NOW B67.** Matt accepted the container render of the shipped file (A35 to A38; `05_output/20261017/pooled-CPUE-2024-25`, port **86,113 [72,407, 103,557]**, cross-check 85,878, -0.27%, PASS; code `91ab2a2`, clean tree, R 4.3.3 and rstan 2.32.5). The box, the register, the method document, CLAUDE.md, run_config.R, the README and the guide now name it; B65 (87,932) is the run it replaced. Open items and next steps are in the PIPELINE_STATUS backlog and the register's D rows (D39 with OSP, D43 to D47 sampling-design items, D24, D29's stale span, D41/D42 text items, the renv.lock audit). Earlier (`1d61fef`): the shipped file RENDERED on the cloud container (B67; `--model both` from `91ab2a2`, clean tree, per-fit checkpoints via the new `BSS_FIT_CACHE_DIR`, B66; R 4.3.3 and rstan 2.32.5 against the lock's R 4.2.2 and rstan 2.32.7). The A35/A36 rule is met on every clause: 8 of 8 fits pass the gate, no stuck chain, cross-check 85,878 against 86,113 (-0.27%), `B1_sun` clean on all eight fits. Port **86,113 [72,407, 103,557]** against B65's 87,932 (-2.1%): the shore all-gear fit 32,455 to 30,804 (a Sunday at 0.71 of a Saturday on four sampled Sundays), the boat all-gear 35,211 to 35,403 (the boat's Sunday 0.89 beside OSP's 0.87), A38 +4.8% on the port SD. Outputs committed as produced (`05_output/20261017/`). B65 stays the authoritative run until Matt says whether this container render stands or he re-renders on the lock versions. Earlier (`8ac877d`): the 17-arm container batch for A35 and A36 is read (VALIDATION_CAMPAIGN 1z.12): every arm clean; the holiday prior moves nothing (-0.2% on the port, seed noise); the Sunday split moves the shore all-gear fit about -5% (a Sunday at 0.71 of a Saturday on the four sampled Sundays; the boat's 0.89 on both tracks, beside OSP's raw 0.87) and the port about -2.2% (container base 87,973 to 86,059 under both); `B1_sun` clean on all eight fits; the gear boat all-gear arm -0.70% on the pooled. A38 measured on the D44 shore refit at +4.8% on the port SD (7,679 to 8,048), the median +0.06%. A35 to A38 stay PENDING Matt's render of the shipped file (`Rscript run_estimation.R --model both`) under the A35/A36 rule in CHANGE_REGISTER A36. Earlier the same day: A38 (`933da04`, Matt's option (i) for D44): the port total pairs the two shore fits' draws comonotonically in `tau_bar` (`port_tau_coupling = "comonotone"`), the interval only; D44 and D48 closed (the Pollock attribution dropped). Before that (`83d3972`): A35 (the holiday prior at normal(0, 0.5), `holiday_prior_sd`), A36 (the Saturday/Sunday split, `sunday_split`) and A37 (the pooled prep's boat I/E path; D35 closed) are adopted in `run_config.R` at `a4350dc` and PENDING ONE RENDER under the A35/A36 rule written in CHANGE_REGISTER before any fit; the container refits are running and their results will be written to the register and here. The 2026-10-16 branch status report's documentation items are applied: the branch text (every commit since 2026-09-29 is on this branch, not `main`), D43 to D48 registered from the 14 to 15 October sampling-design findings and measured against the repository's data, and this PR retitled. The authoritative run is unchanged: B65, 87,932 [75,193, 105,271], cross-check -0.03%.**

**Status (2026-10-04, `9129535`): A33 and A34 rendered and adopted by rule (B65); the authoritative run is 87,932 [75,193, 105,271], cross-check -0.03%. Open: D39 (with OSP), D29's span (stale).** Before that, 2026-10-02 (`0941849`): the batch is finished, every decision is made, and the documentation is audited (B63). Stage A is the authoritative run (99,822 [82,090, 124,718]). D3 and D6 are adopted together on the gear track (A32), which moves the cross-check to -0.53%. D29 keeps the boat all-gear fit monthly, and its rung span is reported. The runner now refuses to fit. Details are in "Update 2026-10-02 (B62)" below; the sections after it are the record of how the batch got here.

Everything needed for the overnight run. The runner renders the first authoritative run at the A31 code (stage A), plus the runs that close D29, D3 and D6, and applies decision rules written in its header before the run. This PR also carries a documentation audit and the `run_config.R` preparation.

**On `init_r = 0.5`.** Matt's render at 0.5 (B57, `05_output/20260928/pooled-CPUE-2024-25-222347`) failed the shore pot-closure gate and was never the authoritative run. The radius did not cause the failure. The pre-A31 model had unidentified single-section level funnels that trapped a chain at either radius: at 2, the shore all-gear fit ran at 4.07% divergences in the authoritative run, and a container refit got stuck at 23.3%. Moving to 0.5 moved the stuck chain from one fit to another. `main` now carries A31 (both funnels collapsed) with 0.5, and that combination has never had a production render. Stage A was that render and is now the authoritative run (99,822); R2 re-rendered it at 2 and found the radius is a precaution, not a crutch.

### Update 2026-10-04 (B65): the A33/A34 render is the authoritative run

- **Render.** Matt rendered from `4490fbe` on a clean tree with the `renv.lock` versions, and the outputs are committed as `c2cb6af`.
- **It meets the A33/A34 rule on every clause:**
  - all 8 fits pass the gate;
  - no chain is stuck;
  - the cross-check is **-0.03%** (87,903 against 87,932);
  - `beta_fvol` is -0.758 [-0.998, -0.522], excluding zero on both tracks.
- **What moved.** The port total goes from 99,822 to **87,932 [75,193, 105,271]** (-11.9%), and the interval is 29% narrower.
  - All of the change is in the boat: the boat all-gear fit goes from 47,105 to 35,211.
  - The shore fits are byte-identical to the previous run.
  - The BSS now sits 0.9% below the PE; before, it was 12.5% above.
- **Still open:**
  - D39 is with OSP.
  - D29's resolution span is stale: it was measured before A33 and A34.
  - In June and July, OSP's own crab-only share sits just below the model's interval (D11).
- **Documentation and harness:** box, register (B65), campaign 1z.11, CLAUDE.md, `run_config.R` header, method document, guide and README all updated (`9129535`). The harness passes all 1,520 checks.

### Update 2026-10-04 (A34): D40 fixed in code; D39 is a question to OSP

- **A34: the crabbing share's day-level volume term (`3b20967` to `4490fbe`).** The share is now `logit f_d = eta_f[month] + beta (log V_d - the month's centre)`.
  - The classification streams use observed volumes: OSP's daily total, or else trailers x turnover centre.
  - The totals use the model's own daily boat volume, so the share they apply is boat-weighted.
  - It is the same on both tracks and ships on (`crab_fraction_volume = TRUE`).
- **Exact when off:** the boat all-gear fit reproduces the pre-A34 fit draw for draw.
- **Container refits** (production sampler settings; these are not renders):

  | Fit | With A34 | Without A34 |
  |---|---|---|
  | Pooled boat all-gear | 35,519 [25,613, 50,023] | 46,538 (A33 alone) |
  | Gear-track boat all-gear | 35,133 (-1.1% vs pooled) | |
  | Pooled boat pot closure | 1,381 | 1,378 |

  - `beta` = -0.76 [-1.00, -0.52] on the pooled fit, and the same on the gear track.
  - No divergences on either all-gear fit.
  - The port total would fall by roughly 11%.
- **The A33 rule was revised before any render.**
  - The winter-`f` clause is retired: the volume term changes `f` in every month by design, and the container refit moves February +0.6 SD.
  - It is replaced by: `beta_fvol`'s 95% interval excludes zero on both tracks' boat all-gear fits.
- **Runners:** the two live batch runners now pin the new key off, and the improvement ladder also pins `use_osp_crab_lower` off, which A33 had leaked into it.
- **D39:** a question for OSP: is a day missing from the crab sheet a zero, or unclassified?
- The harness passes all 1,511 checks.

### Update 2026-10-04 (A33, D39 to D42): `use_osp_crab_lower` on, and two defects in what it fits

- **A33 (`31fd054`).** On Matt's instruction `use_osp_crab_lower = TRUE`: OSP's `WPTPrivateCrabOnly` (delivered 2026-10-02, B64) enters both tracks' private-boat fits as `BetaBinomial(n, f (1 - c), kappa_O)`, with n the boats sampled at OSP's minimum rate. Rule A33 was written before the render: every fit on both tracks passes the gate with no stuck chain; the cross-check is within 2%; December to February `f` moves less than 0.25 posterior SD. Container refits of the boat all-gear fit: 46,947 off and 46,538 on; 43,006 with every returning boat as n. The authoritative run (99,822) stands until the render passes, and the box says the shipped file differs from it in this one key.
- **D39 and D40 (`7aea940`), found by the review after the switch, both checked against OSP's data. Rule A33 tests neither.**
  - **D39.** Only 1 of OSP's 253 crab-sheet rows is a zero, and the absent days lean low-volume. The sheet looks built from crab trips, so the beta-binomial may be fitting zero-truncated data. This is a question for OSP.
  - **D40.** The per-day likelihood fits a day-averaged share, but the harvest needs the boat-weighted one. OSP's share falls on busy days, so in June to September the share the stream targets is 1.41x to 1.59x the boat-weighted share, and those months carry most of the boat catch.
  - **Recommendation:** settle or accept both before the render.
- **D41 and D42 (`cdbbf75`).** Reporting and robustness gaps on the OSP stream, plus stale comments and specification text. None of them affects a fit.
- **Corrected:** "the winter months carry most of the boat catch". December to February are about 4,700 of the boat's 46,000; May to August are about 36,000.
- The harness passes all 1,491 checks.

### Update 2026-10-02 (B63): housekeeping and documentation audit
- **Deleted:** the batch's six unstamped partial output folders (`58b830b`).
- **Documentation audit:** about 90 stale statements fixed across every document. The box now quotes every component on the full posterior: the batch tables use the 2,000 saved draws, so they gave 47,192 for a fit `component_totals.csv` gives as 47,105.
- **D18 corrected.** On 2025-26 as delivered, the boat all-gear fit has no OSP day, because all 23 OSP days fall in the pot closure. Its shared turnover is therefore refused, not fired as the docs had said since 2026-09-28. A refusal now reaches `run_warnings.csv`.
- **Harness:** 1,491 passed, 0 failed.

### Update 2026-10-02 (B62, A32): the batch finished; D3 and D6 adopted together; D29 kept monthly
- **The batch did not stall.** D3R, D6R, D29W and D29B rendered overnight in 35, 35, 244 and 225 min, and all five summary files were written.
- **D3: adopted only together with D6 (A32).**
  - D3R passed every clause, so the rule read ADOPT.
  - But the shipped `bss_seed` (20260619) is the seed at which D3 alone trapped chain 1 of the gear shore all-gear fit (1,546 of 2,000 draws divergent).
  - Renders are reproducible (B38), so adopting D3 alone would have put that fit on its PE in every cross-check.
  - D3 with D6 passed at both seeds.
  - B61's rule should have caught this; it is now a C row.
- **D6: adopted.** D6R against D3R, a converged pair: +10.1 nats at 2.14 paired SE, with every clause passing.
- **New cross-check: 99,294 [82,154, 123,380], -0.53%** (it was -1.44%).
  - It is stage D6, which is the shipped file's gear render: all 190 config keys are identical, and only `run_tag` differs.
  - The pooled render, and so the port total, does not change.
- **D29: KEEP monthly** (rule D29-4 b, Matt's option a).
  - The adequate rungs give 47,192 (monthly), 50,888 (biweekly) and 53,546 (weekly), on the batch's saved draws (full posterior: 47,105, 50,800, 53,651). Held-out weeks give no evidence either way.
  - The box quotes the port span, 99,822 to 106,380 (+6.6% at weekly), as resolution uncertainty outside the interval.
- **Closed out (B62).**
  - The batch runner and the 2026-09-13 D3/D6 ladder now refuse to fit.
  - Two output-identical fixes that had been held until the batch ended:
    - the dplyr 1.2.0 `case_when()` deprecation in `prep_days_crab()`;
    - a stage process that inherited `RENV_PROJECT` would have restored packages into its parent's library.
- **Harness:** 1,488 passed, 0 failed (section 95 is new).

### How to run it (historical: the runner now refuses to fit)
After `git pull`, restart R (the one-time `renv::restore()` was done 2026-09-30), do the dry run, then in the RStudio **Console** (the route that works on Matt's machine: `Rscript` is not on the Terminal's PATH there, and the stages inherit the Console's pandoc and package library):
```r
Sys.setenv(BSS_BATCH_GO = "1"); source("06_diagnostics/run_authoritative_batch_2026-09-29.R")
```
Or, where `Rscript` and pandoc are on the PATH: `Rscript 06_diagnostics/run_authoritative_batch_2026-09-29.R --go` (without `--go`, a dry run).
- The Console stays busy for the whole batch. Each stage's progress goes to its log in `05_output/authoritative_batch_2026-09-29_logs/`, not to the Console. RStudio's "session suspend timeout paused: a child process is running" is expected: the stage is that child process.
- Keep the machine from sleeping overnight. Do not pull or edit tracked files while it runs.
- Start it with the environment variable or `--go`. **Do not edit `DRY_RUN`**: stage A's manifest records any tracked file that differs from the commit.
- Stop it any time: press Stop, then end any leftover `Rscript.exe` in Task Manager. Running the same line again resumes.
- Outputs: `05_output/authoritative_batch_2026-09-29_{stages,components,verdicts,recommendation,d29_rungs}.csv`, plus a log per stage in `05_output/authoritative_batch_2026-09-29_logs/`.

### Update 2026-10-02: the batch's results so far, the authoritative run moves (B61), and two fixes (B60, B61)
- **Stage A is now the authoritative run.** `05_output/20260929/pooled-CPUE-2024-25-220449`, port **99,822 [82,090, 124,718]**, cross-check -1.44% PASS, rendered from `2523e8e` on a clean tree.
  - It met rule A: every fit on both tracks BSS, and shore all-gear 0.49% divergent (4.07% before A31).
  - Every component is within 0.02 posterior SD of the 2026-09-28 run it replaces (99,873).
  - The box, the register, campaign 1z.9, CLAUDE.md, `run_config.R`'s header and the method documents now name it. Only comments changed in `run_config.R`: its values are identical, so every completed stage still resumes.
- **R2: PASS.** At `init_r = 2` every gate passes and every component is within 0.08 SD of A, so the radius is a precaution, not a crutch.
- **B60: the pooled D29 ladder ran out of memory at its second rung**, in every attempt. The process held five posteriors where A, R2 and D29D hold four. It is replaced by forced weekly and biweekly renders (D29W, D29B), with stage A as the monthly rung. The rule's thresholds are unchanged.
- **B61: D3 and D6 are re-rendered once at a second seed** (D3R, D6R), with the rule written into the header before they run.
  - The first D3 failed D3-1 on one stuck chain in the gear shore all-gear fit.
  - D6 was judged against that unconverged fit, so it gains clause D6-0 (a converged baseline, else UNDECIDED); the first D6 reads UNDECIDED.
- **Left to render:** D3R, D6R, D29W and D29B (rendered overnight; see B62 above).

### The first overnight attempt (2026-09-29, `be906de`) failed before any fit; fixed (B59)
- **Cause.** Each stage is an `Rscript` started in the project folder, so the project's `.Rprofile` activated renv in it. The stages therefore loaded packages from `renv/library`, which on Matt's machine is incomplete and does not match the lock: no tidyverse, rstan or loo, and lubridate 1.9.5 against the lock's 1.9.1.
- The authoritative run's recorded versions (2026-09-28 manifest) match `renv.lock` in every package but four patch-level utilities. The project's renv library did not.
- **What failed.** Stage A restored about 90 packages (24.7 min) and then died in `bss_load_packages()`, whose check had loaded the very namespaces the restore replaced. The other six stages died in seconds at the `here` check. S0 had passed because it runs in the Console, not in a stage's process.
- **Fix.**
  - Every stage starts with renv's autoloader off and exactly the starting session's `.libPaths()`.
  - A new S0 check starts a process the same way and requires every package to load there at the session's version; a FAIL blocks every stage.
  - `bss_load_packages()` checks whether a package is installed without loading it.
  - `DRY_RUN` is back to `TRUE`.
- **Tested in the container:**
  - a mismatched library FAILs S0 and blocks all seven stages;
  - Matt's situation (a parent without renv, started from the project folder) PASSes with the fix, where an unfixed child tries to bootstrap renv;
  - a live D3 render through the new path completed in 8.0 min with no renv activity.
- **Second attempt (2026-09-30), in a fresh renv-active session: S0 stopped it before any fit, as designed.** `here` loaded nowhere: the renv library held here 1.0.2 beside rprojroot 2.0.3, left by the first attempt's partial restore.
  - S0 now also compares the session's library with `renv.lock`. A Stan-toolchain difference FAILs and blocks every stage; any other difference is listed as REVIEW.
  - The lock is the right reference: the authoritative run's sessionInfo matches it in every package but four patch-level utilities.
  - **Before the next attempt:** run `renv::restore()` in a fresh R session in the project, restart R, then do the dry run.

### The stages (each in its own R process, under a wall-clock limit)
| stage | what it renders | expected time on your machine |
|---|---|---|
| S0 | desk checks on the real inputs, no MCMC | seconds |
| **A** | `run_estimation.R --model both` on `run_config.R` as shipped; done, and the authoritative run since B61 | ~1.5-2 h |
| D3 | gear track with a per-population period; only the shore all-gear fit moves (monthly to weekly) | ~0.7 h |
| D6 | D3 + the gear-track zero-inflated shore catch, judged against D3 | ~0.7 h |
| D3R, D6R | D3 and D6 again at `bss_seed + 1`; they decide D3 and D6 (B61) | ~0.5 h each |
| D29W, D29B | pooled boat all-gear fit forced to weekly / biweekly; stage A is the monthly rung (B60; replaced the in-process ladder D29P) | ~4 h each |
| D29G | the same ladder on the gear track (corroboration only) | ~1.5 h |
| R2 | pooled at `init_r = 2` | ~1-1.5 h |
| D29D | pooled boat all-gear fit forced to daily | ~2-4 h |

Measured on Matt's machine: A 5.1 h, D3 0.5, D6 0.6, D29G 1.2, R2 3.7, D29D 5.0, D3R 0.6, D6R 0.6, D29W 4.1, D29B 3.7. All done.

### The decision rules (full text in the runner's header)
- **A:** every fit on both tracks passes the gate; shore all-gear under 1% divergent; no stuck chain; cross-check within 2% → ELIGIBLE to replace the box. A4 (movement against 2026-09-28, in posterior SD) is read and never blocks.
- **R2:** at `init_r = 2`, all gates pass, shore all-gear under 1%, no stuck chain, and every component within 0.2 SD of A → the radius is a precaution; otherwise keep 0.5.
- **D3:** every gate passes, the shore all-gear fit at weekly is adequate, and every fit ran at the period asked for → ADOPT the per-population `gear_period_bss`. Decided by D3R (B61).
- **D6:** paired catch elpd gain > 2 paired SE, both count bins improve, and the gain is not bought entirely at the zeros → ADOPT `catch_zi_tracks = c("pooled", "gear_resolved")`. Any clause at REVIEW → UNDECIDED. B61 adds D6-0: the D3 baseline's shore all-gear fit must have converged, else UNDECIDED. Decided by D6R against D3R.
- **D29:** the decision is made on the pooled track only, by the leave-one-week-out block CV of each finer rung against the monthly rung, with the gate and adequacy as preconditions.
  - (a) The finest eligible, adequate rung that interpolates BETTER, with nothing finer left undecided → move the cap to it.
  - (b) No finer rung does better → KEEP monthly, and report the span across the adequate rungs as this component's resolution uncertainty. Applying the shore precedent is Matt's decision.
  - (c) Otherwise → REVIEW.
- Every clause is PASS, FAIL or REVIEW. A file that is missing or a number that cannot be computed gives REVIEW, never a pass or a fail.

### Driver additions (inert outside a ladder)
- Under a ladder, both drivers write every rung's block CV (`ladder_block_<stream>_<fit>_<rung>.csv`) and log `cov50_trailer` per rung. Before, only the kept rung had a block CV.

### `run_config.R`
- Ships the method of record as stage A renders it: `model = "both"`, `run_tag = "2024-25"`, `mu_hier_collapse_single = "both"`, `bss_init_r = 0.5`.
- Since A32 (2026-10-02): `gear_period_bss` per population (shore weekly / biweekly, boat month / biweekly) and `catch_zi_tracks = c("pooled", "gear_resolved")`; the resolved file equals stage D6's configuration on all 190 keys.
- Header rewritten: B51 and A31 shipped after the 2026-09-28 run; B57 failed; stage A (99,822) is now the authoritative run (B61).
- Stale statements corrected: the four v2.0 changes, the `params_model` merge being intended-disjoint rather than guaranteed, and the -1.29% cross-check (B50).

### Documentation audit
- About 55 findings fixed across the method and gear documents, `PIPELINE_STATUS`, `CHANGE_REGISTER` (B57 and B58 rows; D3, D6 and D29 marked staged), `VALIDATION_CAMPAIGN` (1z.8), `CLAUDE.md`, `NEW_SEASON_GUIDE`, both histories and the READMEs.
- Most were statements that B51, A31, B55 or B57 had made stale. For example, several places said a render of the shipped file reproduces 99,873 to the crab; it does not, because A31 and 0.5 change every fit's draws.

## Validation
- **Dry run in the container:** every desk check but the new lock check PASSes on the real 2024-25 inputs, clean tree included. The lock check FAILs there, correctly: the container is not the lock's environment.
- **Harness:** 1,491 passed, 0 failed, 0 skipped at `0941849` (1,469 at B61). Sections 95 and 96 pin B62, A32 and B63; section 93 is new with the runner; it runs the D29 decision code itself on seven synthetic ladders and the manifest parser on a multi-line error. It also starts a process the way a stage is and compares its library paths with the session's (B59).
- **Adversarial review of the runner:** 13 defects found and fixed before handover. Among them:
  - the D29 rule could reach KEEP vacuously when a finer rung's adequacy was NA;
  - D6 could read ADOPT with a clause at REVIEW;
  - a desk FAIL for one stage stopped every stage;
  - a two-line error hid the second model from the manifest parser;
  - there was no timeout, no partial resume of A, and no pandoc check.
- **Timeout check:** on Linux, `system2(timeout =)` kills the stage's forked Stan workers too; I checked for orphans and found none. On Windows this was **not** tested: if a stage times out there, check Task Manager for stray `Rscript` processes before the next stage's fits compete with them.
- **Container smoke run of every stage**: 300 iterations, so these numbers test the machinery, not the estimates. The runner was killed mid-batch and restarted with `--go`: A was reused through RESUME, then every other stage rendered and every verdict was computed from real output.

| stage | smoke result (300 iterations; not evidence for the estimates) |
|---|---|
| A | A1-A6 PASS, ELIGIBLE; every component within 0.05 SD of 2026-09-28; cross-check -1.70% |
| D3 | D3-1..3 PASS, ADOPT; the three unchanged fits are byte-identical controls; gap -1.70% → -0.76% |
| D6 | D6-1..4 PASS, ADOPT; +10.5 nats at 2.25 SE (2026-09-14: +11.3 at 2.29); zero z 3.65 → 2.09, one z -6.10 → -3.49 |
| D29P | (replaced by D29W/D29B in B60: on Matt's machine it ran out of memory) complete ladder; weekly and biweekly NO EVIDENCE against monthly → KEEP (b); adequate rungs span 47,052 to 52,851; ladder monthly = A's monthly exactly (control) |
| D29G | gear ladder monthly = A's gear boat fit exactly (control); biweekly failed the gate at 300 iterations and is marked "out" |
| R2 | R2-1..4 PASS; components within 0.12 SD of A |
| D29D | ran at daily as forced; inadequate (p_loo_frac 0.241, 13 bad k), so "out"; D29 overall → KEEP (b) |

https://claude.ai/code/session_01QVJ6LiXnWdnXF5eY8VtCdL

---
_Generated by [Claude Code](https://claude.ai/code)_
