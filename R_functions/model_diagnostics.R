# =============================================================================
# R_functions/model_diagnostics.R
# Per-fit model-behavior diagnostics for the pooled BSS. Auto-sourced by the
# RMD (it sources everything in R_functions/). All functions are read-only on
# the fit and safe to call on any stanfit from this model lineage.
#
# Produced per fit (see write_bss_diagnostics):
#   structural_params_<label>.csv   key scale/structural parameters with CI,
#                                   n_eff, R-hat (the tuning knobs)
#   divergence_localization_<label>.csv  where divergent draws sit vs the bulk,
#                                   plus the divergent-vs-bulk shift in the totals
#   ppc_calibration_<label>.csv     posterior predictive coverage + PIT for
#                                   effort counts and interview catches
#   ppc_pit_<label>.png             PIT histograms (uniform => well calibrated)
# =============================================================================

# --- 1. Structural / scale parameter summary --------------------------------
bss_structural_summary <- function(fit) {
  pars <- c("mu_mu_E", "mu_mu_C",
            "sigma_eps_E", "sigma_eps_C",
            "phi_E", "phi_C",
            "sigma_r_E", "sigma_r_C", "r_E", "r_C",
            "sigma_mu_E", "sigma_mu_C",
            "sigma_IE", "R_G", "R_T",
            "B1", "B2", "B1_C")
  pars <- pars[pars %in% fit@model_pars]
  s <- summary(fit, pars = pars)$summary
  out <- data.frame(parameter = rownames(s),
                    mean   = s[, "mean"],
                    lo95   = s[, "2.5%"],
                    median = s[, "50%"],
                    hi95   = s[, "97.5%"],
                    n_eff  = round(s[, "n_eff"]),
                    Rhat   = round(s[, "Rhat"], 4),
                    row.names = NULL)
  out
}

# --- 2. Divergence localization (where do divergences sit?) ------------------
bss_divergence_localization <- function(fit, candidate_pars = NULL) {
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- unlist(lapply(sp, function(x) x[, "divergent__"])) == 1
  n_div <- sum(divergent); n_tot <- length(divergent)
  res <- list(n_divergent = n_div, n_total = n_tot,
              frac_divergent = n_div / n_tot,
              distortion_C = NA_real_, distortion_E = NA_real_, table = NULL)
  if (n_div == 0 || n_div == n_tot) return(res)

  if (is.null(candidate_pars))
    candidate_pars <- c("sigma_eps_E", "sigma_eps_C", "sigma_r_E", "sigma_r_C",
                        "r_E", "r_C", "phi_E", "phi_C", "sigma_mu_E",
                        "sigma_mu_C", "sigma_IE", "R_T", "R_G")
  candidate_pars <- candidate_pars[candidate_pars %in% fit@model_pars]
  flat <- function(p) as.vector(as.array(fit, pars = p)[, , 1])  # chain-major

  rows <- lapply(candidate_pars, function(p) {
    x <- flat(p); xd <- x[divergent]; xn <- x[!divergent]; s <- stats::sd(x)
    smd <- if (is.finite(s) && s > 0) (mean(xd) - mean(xn)) / s else NA_real_
    data.frame(param = p, median_bulk = stats::median(xn),
               median_div = stats::median(xd), smd = smd, abs_smd = abs(smd))
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(-tab$abs_smd), c("param", "median_bulk", "median_div", "smd")]
  res$table <- tab

  for (q in c("C_expected_sum", "E_sum")) {
    if (q %in% fit@model_pars) {
      v <- flat(q); mb <- stats::median(v[!divergent]); md <- stats::median(v[divergent])
      sh <- abs(md - mb) / mb
      if (q == "C_expected_sum") res$distortion_C <- sh else res$distortion_E <- sh
    }
  }
  res
}

# --- 3. Posterior predictive calibration ------------------------------------
# For effort counts (NB(lambda_E * R, r_E)) and interview catches
# (NB(lambda_C * h, r_C)): draw from the posterior predictive for each
# observation, then report central-interval coverage and the PIT. Coverage near
# nominal and uniform PITs indicate the observation model fits. All in R from
# extracted quantities (no RNG added to Stan), capped at n_draws_use for speed.
bss_ppc_calibration <- function(fit, stan_data, n_draws_use = 400, seed = 1) {
  set.seed(seed)
  ex <- rstan::extract(fit, pars = c("lambda_E_S", "lambda_C_S",
                                     "r_E", "r_C", "R_G", "R_T"))
  ndraw <- length(ex$r_E)
  use <- if (ndraw > n_draws_use) sort(sample.int(ndraw, n_draws_use)) else seq_len(ndraw)
  nd  <- length(use)
  lamE <- matrix(ex$lambda_E_S[use, 1, , 1], nrow = nd)  # [draws, D]
  lamC <- matrix(ex$lambda_C_S[use, 1, , 1], nrow = nd)
  rE <- ex$r_E[use]; rC <- ex$r_C[use]; RG <- ex$R_G[use]; RT <- ex$R_T[use]

  calib <- function(mu_mat, y, size_vec) {
    nobs <- length(y); cov50 <- cov95 <- pit <- numeric(nobs)
    for (i in seq_len(nobs)) {
      yp <- stats::rnbinom(nd, mu = pmax(mu_mat[, i], 1e-8), size = size_vec)
      qq <- stats::quantile(yp, c(.025, .25, .75, .975), names = FALSE)
      cov50[i] <- y[i] >= qq[2] && y[i] <= qq[3]
      cov95[i] <- y[i] >= qq[1] && y[i] <= qq[4]
      pit[i]   <- mean(yp < y[i]) + 0.5 * mean(yp == y[i])
    }
    list(summary = data.frame(coverage_50 = mean(cov50), coverage_95 = mean(cov95),
                              pit_mean = mean(pit), pit_sd = stats::sd(pit), n = nobs),
         pit = pit)
  }

  parts <- list()
  if (stan_data$Gear_n > 0)
    parts$gear    <- calib(lamE[, stan_data$day_Gear, drop = FALSE] * RG,
                           stan_data$Gear_I, rE)
  if (stan_data$T_n > 0)
    parts$trailer <- calib(lamE[, stan_data$day_T, drop = FALSE] * RT,
                           stan_data$T_I, rE)
  if (stan_data$IntC > 0)
    parts$catch   <- calib(sweep(lamC[, stan_data$day_IntC, drop = FALSE], 2,
                                 stan_data$h, "*"),
                           stan_data$c, rC)
  if (length(parts) == 0) return(NULL)

  summ <- do.call(rbind, lapply(names(parts), function(nm)
    cbind(data_type = nm, parts[[nm]]$summary)))
  pit_long <- do.call(rbind, lapply(names(parts), function(nm)
    data.frame(data_type = nm, pit = parts[[nm]]$pit)))
  list(summary = summ, pit = pit_long)
}

# --- 4. Write all diagnostics for one fit -----------------------------------
write_bss_diagnostics <- function(fit, stan_data, label, output_dir) {
  ok <- function(expr) tryCatch(expr, error = function(e) {
    cat(sprintf("    [diag] %s skipped: %s\n", label, conditionMessage(e))); NULL })

  ok({
    sp <- bss_structural_summary(fit)
    utils::write.csv(sp, file.path(output_dir, sprintf("structural_params_%s.csv", label)),
                     row.names = FALSE)
  })

  ok({
    dl <- bss_divergence_localization(fit)
    cat(sprintf("    [diag] %s: %d/%d divergent (%.1f%%); distortion C/E = %s/%s\n",
                label, dl$n_divergent, dl$n_total, 100 * dl$frac_divergent,
                ifelse(is.na(dl$distortion_C), "NA", sprintf("%.1f%%", 100 * dl$distortion_C)),
                ifelse(is.na(dl$distortion_E), "NA", sprintf("%.1f%%", 100 * dl$distortion_E))))
    if (!is.null(dl$table)) {
      dl$table$distortion_C <- dl$distortion_C
      dl$table$distortion_E <- dl$distortion_E
      utils::write.csv(dl$table,
                       file.path(output_dir, sprintf("divergence_localization_%s.csv", label)),
                       row.names = FALSE)
    }
  })

  ok({
    if (is.null(stan_data)) {
      cat("    [diag] PPC skipped: stan_data not stored with this fit\n")
    } else {
    ppc <- bss_ppc_calibration(fit, stan_data)
    if (!is.null(ppc)) {
      utils::write.csv(ppc$summary,
                       file.path(output_dir, sprintf("ppc_calibration_%s.csv", label)),
                       row.names = FALSE)
      cat("    [diag] PPC coverage (nominal 50/95):\n")
      for (k in seq_len(nrow(ppc$summary)))
        cat(sprintf("      %-8s 50%%=%.0f%%  95%%=%.0f%%  (n=%d)\n",
                    ppc$summary$data_type[k], 100 * ppc$summary$coverage_50[k],
                    100 * ppc$summary$coverage_95[k], ppc$summary$n[k]))
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        p <- ggplot2::ggplot(ppc$pit, ggplot2::aes(x = pit)) +
          ggplot2::geom_histogram(boundary = 0, bins = 10, fill = "steelblue",
                                  colour = "white") +
          ggplot2::geom_hline(yintercept = 0, colour = NA) +
          ggplot2::facet_wrap(~data_type, scales = "free_y") +
          ggplot2::labs(title = sprintf("PPC PIT: %s", label),
                        subtitle = "Uniform => calibrated; U-shape => underdispersed; hump => overdispersed",
                        x = "PIT", y = "count") +
          ggplot2::theme_bw()
        ggplot2::ggsave(file.path(output_dir, sprintf("ppc_pit_%s.png", label)),
                        p, width = 10, height = 4)
      }
    }
    }
  })
  invisible(NULL)
}
