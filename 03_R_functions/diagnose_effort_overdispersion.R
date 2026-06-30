# =============================================================================
# diagnose_effort_overdispersion.R
#
# T1.5 step 2: decompose the effort-count posterior predictive variance into its
# latent (process/parameter) and observation (negative-binomial) parts, so the
# lever behind the effort over-dispersion seen in the PPC (central-hump PITs,
# coverage_50 ~ 0.63-0.75 vs the nominal 0.50) is identified before any prior or
# model change is made.
#
# WHY THIS CANNOT RUN ON THE COMMITTED CSVs
#   The decomposition needs the JOINT posterior draws of the latent effort
#   intensity lambda_E_S at each effort-observation day together with the r_E,
#   R_G and R_T draws. The committed outputs carry only summaries (structural_*
#   has r_E as a mean; bss_daily_effort_* has E = lambda_E * L as daily
#   quantiles), so this operates on the in-memory stanfit objects.
#
# THE MATH (law of total variance)
#   Each gear count is   Gear_I[i] ~ NB2(mu_i, r_E),  mu_i = lambda_E[d_i] * R_G
#   Each trailer count is T_I[i]   ~ NB2(mu_i, r_E),  mu_i = lambda_E[d_i] * R_T
#   NB2(mu, r) has mean mu and variance mu + mu^2 / r. Integrating the predictive
#   over the posterior of (mu_i, r_E):
#       Var(Y_i) = E[ Var(Y_i | mu_i, r_E) ] + Var( E[Y_i | mu_i, r_E] )
#                = E[mu_i] + E[mu_i^2 / r_E] + Var(mu_i)
#   giving three additive components per observation:
#       V_pois   = E[mu_i]            irreducible Poisson sampling
#       V_nb     = E[mu_i^2 / r_E]    NB observation over-dispersion  (r_E lever)
#       V_latent = Var(mu_i)          process + parameter uncertainty (sigma_eps_E
#                                     / single-cell level lever, T2.5)
#   Summed across observations and normalized, the shares say which reducible
#   component dominates the over-wide predictive:
#       share_nb  > share_latent  ->  the r_E / sigma_r_E prior is the lever (T1.5
#                                     option 5a; pair with the T1.3 sweep)
#       share_latent > share_nb   ->  the latent scale sigma_eps_E or the T2.5
#                                     single-cell level redundancy is the lever
#   V_pois is irreducible (Poisson floor); it is reported for context, not as a
#   lever.
#
# WHAT THIS DOES NOT DO
#   It does not, by itself, prove the model is over-dispersed; the PPC
#   (ppc_calibration_*) establishes that. This attributes the over-dispersion to
#   a tunable part. It changes nothing in the model and writes only CSVs, so it
#   is safe to add to any run. A correction (tightening sigma_r_E, or addressing
#   the latent scale) is a separate, sign-off-gated step (T1.5 step 5).
#
# USAGE
#   Pipeline mode (recommended; writes a committed output every run). Add a chunk
#   after section 7.12 (the per-fit model diagnostics):
#
#     ```{r effort_overdispersion, eval = run_bss}
#     cat("\n=== Effort over-dispersion decomposition (T1.5) ===\n")
#     for (label in names(bss_all)) {
#       b <- bss_all[[label]]
#       if (is.null(b$fit)) next                 # PE-only entry: nothing to fit
#       write_effort_overdispersion_diag(
#         b$fit,
#         if (!is.null(b$bss_data)) b$bss_data else NULL,
#         label, output_dir)
#     }
#     ```
#
#   Post-run / same-session mode (no re-run; bss_all still in memory):
#
#     write_effort_overdispersion_diag(
#       bss_all[["shore_all_gear_Dungeness_Kept"]]$fit,
#       bss_all[["shore_all_gear_Dungeness_Kept"]]$bss_data,
#       "shore_all_gear_Dungeness_Kept", output_dir)
#
# OUTPUT
#   effort_overdispersion_decomp_<label>.csv  one row per data_type (gear/trailer)
#     with n_obs, n_draws_used, rE_median, R_median, mean_mu, the three variance
#     shares, the predictive and NB-only variance-to-mean ratios, the attached PPC
#     coverage_50 / pit_sd if ppc_calibration_<label>.csv is present, and the
#     lever verdict.
#   effort_overdispersion_byobs_<label>.csv   per-observation components (small;
#     n_obs is in the low hundreds) for plotting or closer inspection.
#
# Auto-sourced by the RMD's purrr::walk(list.files(here("R_functions"), ...),
# source). Every piece is tryCatch-wrapped so it cannot break a run.
# =============================================================================

# Decompose one set of effort observations (a single data_type) into the three
# variance components. obs_days is an integer vector into 1:D; Rdraws is the
# length-nd expansion-ratio draw vector (R_G for gear, R_T for trailer); rEdraws
# is the length-nd r_E draw vector; lamE is the [nd, D] latent-intensity matrix.
.eod_decompose <- function(lamE, obs_days, Rdraws, rEdraws, min_valid = 20) {
  nd <- nrow(lamE)
  lam_obs <- lamE[, obs_days, drop = FALSE]          # [nd, n_obs]
  # mu[i, j] = lam_obs[i, j] * Rdraws[i]. Column-wise recycling of a length-nd
  # vector down an nd-row matrix gives exactly this per-draw scaling (the same
  # recycling verified in the PPC fix). Rdraws and rEdraws are length nd.
  mu  <- lam_obs * Rdraws
  mu[!is.finite(mu)] <- NA                            # weakly-identified fits can
                                                      # overflow exp(); drop those
                                                      # draws per observation
  mu2 <- mu * mu
  mu2_over_r <- mu2 / rEdraws                         # mu2[i,j] / rEdraws[i]

  n_valid <- colSums(is.finite(mu))
  E_mu  <- colMeans(mu,         na.rm = TRUE)
  E_mu2 <- colMeans(mu2,        na.rm = TRUE)
  E_nb  <- colMeans(mu2_over_r, na.rm = TRUE)
  Var_mu <- pmax(E_mu2 - E_mu^2, 0)                   # clamp tiny negative (fp)

  keep <- is.finite(E_mu) & is.finite(E_nb) & is.finite(Var_mu) &
          (n_valid >= min_valid) & (E_mu > 0)
  V_pois <- E_mu[keep]; V_nb <- E_nb[keep]; V_latent <- Var_mu[keep]
  days_k <- obs_days[keep]
  V_total <- V_pois + V_nb + V_latent

  S_pois <- sum(V_pois); S_nb <- sum(V_nb); S_lat <- sum(V_latent)
  S_tot  <- S_pois + S_nb + S_lat

  list(
    n_obs        = length(V_pois),
    n_dropped    = sum(!keep),
    mean_mu      = mean(V_pois),
    share_poisson         = if (S_tot > 0) S_pois / S_tot else NA_real_,
    share_nb_overdisp     = if (S_tot > 0) S_nb  / S_tot else NA_real_,
    share_latent          = if (S_tot > 0) S_lat / S_tot else NA_real_,
    pred_vmr_mean = mean(V_total / V_pois),           # predictive variance / mean
    nb_vmr_mean   = mean((V_pois + V_nb) / V_pois),   # NB-only (1 + E[mu^2/r]/E[mu])
    byobs = data.frame(
      day = days_k, mean_mu = V_pois,
      V_poisson = V_pois, V_nb_overdisp = V_nb, V_latent = V_latent,
      V_total = V_total,
      pred_vmr = V_total / V_pois,
      stringsAsFactors = FALSE
    )
  )
}

write_effort_overdispersion_diag <- function(fit, stan_data, label, output_dir,
                                             n_draws_use = 2000) {
  res <- tryCatch({
    if (is.null(fit)) { cat(sprintf("  %s: no fit; skipped.\n", label)); return(invisible(NULL)) }
    if (is.null(stan_data)) {
      cat(sprintf("  %s: stan_data not carried on the fit; skipped (add bss_data to bss_all).\n", label))
      return(invisible(NULL))
    }

    ex <- rstan::extract(fit, pars = c("lambda_E_S", "r_E", "R_G", "R_T"))
    ndraw <- length(ex$r_E)
    use <- if (ndraw > n_draws_use) sort(sample.int(ndraw, n_draws_use)) else seq_len(ndraw)
    nd  <- length(use)

    # Same hardened reduction as the PPC: lambda_E_S is array[S] matrix[D,G];
    # rstan returns [iter, S, D, G]; reduce to [draws, D] for the S=1, G=1 pooled
    # model and guard dropped size-1 dims.
    get_lam <- function(arr) {
      d <- dim(arr)
      m <- if (length(d) == 4) arr[use, 1, , 1]
           else if (length(d) == 3) arr[use, , 1]
           else if (length(d) == 2) arr[use, ]
           else stop(sprintf("unexpected lambda dims: %s", paste(d, collapse = "x")))
      matrix(m, nrow = nd)
    }
    lamE <- get_lam(ex$lambda_E_S)                    # [draws, D]

    # Scalars come back as length-1-dim arrays; as.numeric strips the stray dim
    # (the Fix-3 lesson) so the column-wise recycling in .eod_decompose is correct.
    rE <- as.numeric(ex$r_E[use])
    RG <- as.numeric(ex$R_G[use])
    RT <- as.numeric(ex$R_T[use])

    # Optionally attach the PPC miscalibration for side-by-side reading.
    ppc_path <- file.path(output_dir, sprintf("ppc_calibration_%s.csv", label))
    ppc <- if (file.exists(ppc_path)) {
      tryCatch(utils::read.csv(ppc_path, stringsAsFactors = FALSE), error = function(e) NULL)
    } else NULL
    ppc_lookup <- function(dt) {
      if (is.null(ppc)) return(c(coverage_50 = NA_real_, pit_sd = NA_real_))
      row <- ppc[ppc$data_type == dt, , drop = FALSE]
      if (nrow(row) == 0) return(c(coverage_50 = NA_real_, pit_sd = NA_real_))
      c(coverage_50 = row$coverage_50[1], pit_sd = row$pit_sd[1])
    }

    verdict <- function(d) {
      if (!is.finite(d$share_nb_overdisp) || !is.finite(d$share_latent)) return("indeterminate")
      if (d$share_nb_overdisp >= d$share_latent)
        "observation NB over-dispersion dominates: lever is r_E / sigma_r_E prior (T1.5 5a; pair with T1.3)"
      else
        "latent process variance dominates: lever is sigma_eps_E or the single-cell level (T2.5), not r_E"
    }

    types <- list()
    if (!is.null(stan_data$Gear_n) && stan_data$Gear_n > 0)
      types[["gear"]]    <- list(days = as.integer(stan_data$day_Gear), R = RG)
    if (!is.null(stan_data$T_n) && stan_data$T_n > 0)
      types[["trailer"]] <- list(days = as.integer(stan_data$day_T),    R = RT)

    if (length(types) == 0) {
      cat(sprintf("  %s: no gear or trailer effort observations; skipped.\n", label))
      return(invisible(NULL))
    }

    summ_rows <- list(); byobs_all <- list()
    D <- ncol(lamE)
    for (dt in names(types)) {
      days <- types[[dt]]$days
      days <- days[is.finite(days) & days >= 1 & days <= D]
      if (length(days) == 0) next
      d <- .eod_decompose(lamE, days, types[[dt]]$R, rE)
      pc <- ppc_lookup(dt)
      summ_rows[[dt]] <- data.frame(
        fit = label, data_type = dt,
        n_obs = d$n_obs, n_obs_dropped = d$n_dropped, n_draws_used = nd,
        rE_median = stats::median(rE), R_median = stats::median(types[[dt]]$R),
        mean_mu = round(d$mean_mu, 3),
        share_poisson      = round(d$share_poisson, 4),
        share_nb_overdisp  = round(d$share_nb_overdisp, 4),
        share_latent       = round(d$share_latent, 4),
        pred_vmr_mean = round(d$pred_vmr_mean, 3),
        nb_vmr_mean   = round(d$nb_vmr_mean, 3),
        ppc_coverage_50 = round(unname(pc["coverage_50"]), 4),
        ppc_pit_sd      = round(unname(pc["pit_sd"]), 4),
        lever = verdict(d),
        stringsAsFactors = FALSE
      )
      bo <- d$byobs; bo$data_type <- dt
      byobs_all[[dt]] <- bo
      cat(sprintf(
        "  %s [%s]: n=%d | shares pois/NB/latent = %.2f/%.2f/%.2f | pred VMR %.2f (NB-only %.2f)%s\n      -> %s\n",
        label, dt, d$n_obs, d$share_poisson, d$share_nb_overdisp, d$share_latent,
        d$pred_vmr_mean, d$nb_vmr_mean,
        if (is.finite(pc["coverage_50"])) sprintf(" | PPC cov50 %.2f, pit_sd %.3f", pc["coverage_50"], pc["pit_sd"]) else "",
        verdict(d)))
    }

    if (length(summ_rows) > 0) {
      summ <- do.call(rbind, summ_rows)
      utils::write.csv(summ, file.path(output_dir,
        sprintf("effort_overdispersion_decomp_%s.csv", label)), row.names = FALSE)
    }
    if (length(byobs_all) > 0) {
      bo <- do.call(rbind, byobs_all)
      bo <- bo[, c("data_type", "day", "mean_mu", "V_poisson", "V_nb_overdisp",
                   "V_latent", "V_total", "pred_vmr")]
      utils::write.csv(bo, file.path(output_dir,
        sprintf("effort_overdispersion_byobs_%s.csv", label)), row.names = FALSE)
    }
    invisible(TRUE)
  }, error = function(e) {
    cat(sprintf("  [eod] %s skipped: %s\n", label, conditionMessage(e)))
    invisible(NULL)
  })
  res
}
