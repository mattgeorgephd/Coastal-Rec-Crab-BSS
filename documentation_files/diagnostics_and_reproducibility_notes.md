# Diagnostics, reproducibility, and documentation: what to add before the run

This sits on top of B1.6 (the Stan files and `B1.6_change_notes.md`). Everything
here is additive and `tryCatch`-wrapped: none of it can break a run or alter the
estimates. If you would rather run lean, the only thing you must add beyond the
Stan swap is the distortion-aware gate from `B1.6_change_notes.md`; the rest of
this file is evaluation and reproducibility polish.

## 0. Reminder: the Stan swap alone does not move shore

B1.6 changes `sigma_IE`, which only matters when `IE_n = 0` (the boat). Shore has
I/E data, so its `sigma_IE` prior was already active and shore is unchanged by the
Stan swap. Shore moves from PE to BSS only via the distortion-aware gate. Minimum
for the gate experiment = Stan swap + the gate edit.

## 1. The diagnostics module

Place `model_diagnostics.R` in `R_functions/`. It is auto-sourced by the RMD's
existing `purrr::walk(list.files(here("R_functions"), ...), source)`. It adds, per
real BSS fit:

- `structural_params_<label>.csv`: `sigma_eps`, `phi`, `sigma_r`/`r`, `sigma_mu`,
  `sigma_IE`, `R_G`, `R_T`, `B1`/`B2`/`B1_C` with 95% CI, n_eff, R-hat. The tuning
  knobs in one place; this is where you would have read `phi_E = 0.83` or
  `R_T = 1.0` off directly instead of running the ad-hoc diagnostic.
- `divergence_localization_<label>.csv`: where divergent draws sit relative to the
  bulk (the smd ranking) plus the divergent-vs-bulk shift in the totals. Makes the
  divergence diagnostic a standard output of every run.
- `ppc_calibration_<label>.csv` and `ppc_pit_<label>.png`: posterior predictive
  central-interval coverage and PIT for effort counts and interview catches.
  Coverage near 50/95 and uniform PITs mean the observation model fits; a U-shaped
  PIT means underdispersed, a central hump means overdispersed. This is the check
  that tells you whether the negative-binomial effort/catch model is right, which
  nothing currently tests.

## 2. The one chunk to add

After the fit loop and before (or after) the convergence report, add:

````r
```{r model_diagnostics, eval = run_bss}
# Additive per-fit diagnostics (safe to remove). Sourced from
# R_functions/model_diagnostics.R; each piece is tryCatch-wrapped.
cat("\n=== Model diagnostics ===\n")
for (label in names(bss_all)) {
  b <- bss_all[[label]]
  if (is.null(b$fit)) next                  # PE-only entry: nothing to diagnose
  cat(sprintf("  %s:\n", label))
  write_bss_diagnostics(b$fit,
                        if (!is.null(b$bss_data)) b$bss_data else NULL,
                        label, output_dir)
}
```
````

The PPC needs the Stan data. If `bss_all[[label]]` does not already carry it, add
one field where the fit is stored in the fit loop, e.g.:

```r
bss_all[[label]] <- list(fit = fit, bss_data = bss_data, <existing fields...>)
```

If you skip that one line, the structural-parameter and divergence-localization
CSVs still write; only the PPC is skipped (with a printed note). Nothing errors.

## 3. Fixed seed (reproducibility)

The `rstan::stan()` calls set no seed, so runs are not reproducible and part of the
292-vs-1669 divergence swing we saw across runs is just seed variation. For a model
others run, and to compare cleanly before/after the gate, fix it.

Add to `params`:

```r
  bss_seed = 20260619,   # fixed RNG seed for reproducible fits (B1.6)
```

Add to the `rstan::stan()` call (one argument):

```r
    seed = params$bss_seed,
```

rstan seeds each chain from `seed + chain_id`, so chains still differ and R-hat is
still meaningful; only the run-to-run randomness is removed. Honest tradeoff: a
fixed seed locks in one realization, so if you ever suspect a pathological seed,
change the number and re-run. For production and debugging, fixed is the right
default, and exposing it as a parameter keeps it changeable.

## 4. Session and version capture (reproducibility)

At the end of the RMD, add:

```r
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
cat(sprintf("\nrstan %s; StanHeaders %s; seed %s\n",
            utils::packageVersion("rstan"),
            utils::packageVersion("StanHeaders"),
            params$bss_seed),
    file = file.path(output_dir, "session_info.txt"), append = TRUE)
```

Pins package and Stan versions and the seed with each output set, so a result is
reproducible from the committed artifacts.

## 5. Documentation: proposed, not baked in (does not need a run)

The diagnostics above cover "evaluation." The "clear assumptions and choices for
others" goal is documentation, and it does not depend on this run, so it should
not gate it. The highest-value artifact is a single decisions log that states the
why behind each major choice in one place:

- three populations estimated independently and summed; comm/charter is census
- two gear sub-seasons split at the Dec 1 pot open, with the per-sub-season gear
  exclusion lists and why phantom pot interviews caused divergences
- gear-hours (not crabber-hours) for boats, L = 24h, and why the civil-twilight
  cap is wrong for crab
- `L_effective` from I/E vs the civil-twilight assumption, and how I/E enters
- the priors and their rationale
- AR resolution: daily for shore, weekly for the boat, and why trailer counts
  cannot identify a daily process
- the PE -> BSS fallback and, now, the distortion-aware gate, with the caveat that
  the distortion check is a safeguard not a proof
- the boat data limitation (the durable fix is 2025-26 exit counts, not parameter
  surgery)
- incomplete-trip filtering and the gear-dependent CPUE bias
- the B-number history (B1.3 omega_0, B1.5 NB marginalization, B1.6 sigma_IE)

Much of this already lives in the change-notes and the documentation.md version
histories; the work is consolidation into one legible file, not new content. I
have not audited the existing documentation.md yet, so before writing anything I
would read it to avoid duplication and instead extend it. Say the word and I will
do that as a parallel track while this run is going, since it needs no run output.

## 6. Runtime

The diagnostics operate on already-fitted objects. The PPC draws from the
predictive for each observation across at most 400 posterior draws, which is well
under a second even for the largest shore fit. No meaningful addition to run time.
