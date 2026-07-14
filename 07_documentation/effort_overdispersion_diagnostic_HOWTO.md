# Effort over-dispersion diagnostic: how to run and read it

- **Companion to:** `R_functions/diagnose_effort_overdispersion.R`
- **Purpose:** T1.5 step 2 (see `development_notes/PIPELINE_STATUS.md`, Section 4). Decompose the effort-count posterior predictive variance into its latent and observation parts so the lever behind the PPC effort over-dispersion is identified before any prior or model change.
- **Date:** 2026-06-22
- **No em dashes per project convention.**

---

## 1. What it answers

The v7.0 PPC showed the negative-binomial effort model is over-dispersed (gear and trailer PIT panels are central humps; `coverage_50` 0.63 to 0.75 against the nominal 0.50; `pit_sd` 0.222 to 0.260 against the 0.289 ideal). That tells us the effort predictive is too wide. It does not tell us which part of the model to tune. This diagnostic splits each effort observation's predictive variance into three additive parts and reports their shares, which is the decisive "which lever" read.

## 2. Why it needs the fits, not the CSVs

The decomposition needs the joint posterior draws of the latent effort intensity `lambda_E_S` at each effort-observation day, plus the `r_E`, `R_G`, `R_T` draws. The committed outputs carry only summaries. So it runs on the in-memory `stanfit` objects. The function placed in `R_functions/` is auto-sourced and is callable in two modes.

## 3. How to run

**Pipeline mode (recommended).** Add this chunk after section 7.12 (the per-fit model diagnostics) in `BSS-GH-pooled-CPUE-model.Rmd`. It writes a committed output every run, alongside the PPC.

````r
```{r effort_overdispersion, eval = run_bss}
cat("\n=== Effort over-dispersion decomposition (T1.5) ===\n")
for (label in names(bss_all)) {
  b <- bss_all[[label]]
  if (is.null(b$fit)) next                 # PE-only entry: nothing to fit
  write_effort_overdispersion_diag(
    b$fit,
    if (!is.null(b$bss_data)) b$bss_data else NULL,
    label, output_dir)
}
```
````

**Post-run / same-session mode (no re-run).** With `bss_all` still in memory after a run:

```r
write_effort_overdispersion_diag(
  bss_all[["shore_all_gear_Dungeness_Kept"]]$fit,
  bss_all[["shore_all_gear_Dungeness_Kept"]]$bss_data,
  "shore_all_gear_Dungeness_Kept", output_dir)
```

Both call the same function. It writes only CSVs and changes nothing in the model, so it is safe to add to any run. It is `tryCatch`-wrapped and cannot break a run.

## 4. The math (verified)

Each effort count is `NB2(mu_i, r_E)` with `mu_i = lambda_E[d_i] * R` (R is `R_G` for gear, `R_T` for trailer). `NB2(mu, r)` has variance `mu + mu^2/r`. Integrating the predictive over the posterior of `(mu_i, r_E)` and applying the law of total variance:

```text
Var(Y_i) = E[mu_i]            (V_poisson:  irreducible Poisson floor)
         + E[mu_i^2 / r_E]    (V_nb:       NB observation over-dispersion)
         + Var(mu_i)          (V_latent:   process + parameter uncertainty)
```

The three components are summed across observations and normalized to shares. The implementation was checked against a brute-force Monte Carlo predictive variance (draw `Y ~ NB2(mu, r_E)` across the posterior, take the variance); the analytic total matches within Monte Carlo noise (mean ratio 1.007 over a test set, shares sum to 1.000).

## 5. How to read the output and what to do next

`effort_overdispersion_decomp_<label>.csv`, one row per `data_type` (gear or trailer):

| Column | Meaning |
|---|---|
| `share_poisson` | Fraction of predictive variance that is the irreducible Poisson floor. Not a lever. |
| `share_nb_overdisp` | Fraction from the NB observation over-dispersion. Controlled by `r_E = 1/sigma_r_E^2`. |
| `share_latent` | Fraction from latent process and parameter uncertainty. Controlled by the AR innovation scale `sigma_eps_E` and the single-cell level redundancy (T2.5). |
| `pred_vmr_mean` | Model predictive variance-to-mean ratio, averaged over observations. A calibrated Poisson would be 1.0; values above 1 are the model's total over-dispersion. |
| `nb_vmr_mean` | The NB-only variance-to-mean ratio, `1 + E[mu^2/r]/E[mu]`. Isolates how much the observation NB inflates beyond Poisson. |
| `ppc_coverage_50`, `ppc_pit_sd` | Pulled from `ppc_calibration_<label>.csv` if present, so the cause (this table) and the symptom (the PPC) sit side by side. |
| `lever` | The verdict string, based on whether `share_nb_overdisp` or `share_latent` is larger. |

**Decision rule (this is the point of the exercise):**

- **If `share_nb_overdisp` > `share_latent`**, the over-wide effort predictive is mostly the observation NB. The lever is the `r_E` / `sigma_r_E` prior. Proceed to T1.5 option 5(a): tighten the `sigma_r_E` prior (a half-normal or exponential informed by the posterior `r_E`), test on shore all-gear first, and run it as the `sigma_r_E` arm of the T1.3 prior sensitivity. Then re-run the PPC (T1.5 step 6) and confirm `coverage_50` moves toward 0.50.

- **If `share_latent` > `share_nb_overdisp`**, the spread is mostly latent. The lever is the AR innovation scale `sigma_eps_E` or the single-cell level redundancy (T2.5), not `r_E`. Tightening `r_E` would not fix it; addressing the latent scale would. Note the T2.5 caveat: the obvious level collapse already broke shore once (B1.7), so a latent-side change needs the careful, guarded approach in T2.5.

- **If the shares are close**, both contribute; sequence 5(a) first (it is the cheaper, exact-prior change) and re-measure before touching the latent side.

**Two standing caveats, both consequential.** Any correction is a prior or inference change, so it follows the design-before-coding norm and needs a guarded test run; the B1.7 episode is the reminder that this model's geometry is delicate. And tightening the effort dispersion narrows the reported effort intervals, including the headline summer intervals (T1.1), which is a change to the reported uncertainty, so it needs the same explicit sign-off the gate change received. The target is calibration (`coverage_50` near 0.50), not zero over-dispersion; some over-dispersion is real.

## 6. Expectation (a hypothesis, not a conclusion)

Because the shore fits are well-identified (n_eff above 4,000, so the latent `lambda_E` is tightly pinned on observed days), I expect `share_latent` to be modest on the observed days and `share_nb_overdisp` to dominate, pointing to the `r_E` / `sigma_r_E` prior as the lever for shore. The boat, being weakly identified, may show a larger latent share. The diagnostic measures this rather than assuming it; read the verdict column and act per section 5.

## 7. Outputs

- `effort_overdispersion_decomp_<label>.csv`: the summary table above, one row per data type.
- `effort_overdispersion_byobs_<label>.csv`: per-observation components (`day`, `mean_mu`, `V_poisson`, `V_nb_overdisp`, `V_latent`, `V_total`, `pred_vmr`). Small (observation counts are in the low hundreds); useful for plotting the components against `mean_mu` or against day, to see whether the over-dispersion is uniform or concentrated in high-count days.
