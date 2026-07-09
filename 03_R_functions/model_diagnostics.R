# =============================================================================
# 03_R_functions/model_diagnostics.R
# Per-fit model-behavior diagnostics for the pooled BSS. Auto-sourced by the
# RMD (it sources everything in 03_R_functions/). All functions are read-only on
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
  # Model-agnostic trailer expansion: pooled carries R_T, gear-resolved carries
  # R_G_boat. bss_extract_pars() requests only the one this fit declares.
  trailer_par <- bss_trailer_par(fit)
  ex <- rstan::extract(fit, pars = bss_extract_pars(fit, c("lambda_E_S", "lambda_C_S",
                                                           "r_E", "r_C", "R_G")))
  ndraw <- length(ex$r_E)
  use <- if (ndraw > n_draws_use) sort(sample.int(ndraw, n_draws_use)) else seq_len(ndraw)
  nd  <- length(use)
  # lambda_*_S is array[S] matrix[D,G]; rstan returns [iter, S, D, G]. get_lam
  # reduces it to [draws, D] for the S=1, G=1 pooled model and guards against
  # dropped size-1 dims so a shape change cannot abort the whole PPC.
  get_lam <- function(arr) {
    d <- dim(arr)
    m <- if (length(d) == 4) arr[use, 1, , 1]
         else if (length(d) == 3) arr[use, , 1]
         else if (length(d) == 2) arr[use, ]
         else stop(sprintf("unexpected lambda dims: %s", paste(d, collapse = "x")))
    matrix(m, nrow = nd)
  }
  lamE <- get_lam(ex$lambda_E_S)   # [draws, D]
  lamC <- get_lam(ex$lambda_C_S)

  # ROOT-CAUSE FIX (supersedes the earlier B1.7 attempt, which hardened lambda
  # extraction and the rnbinom NA path but did not touch the scalars and so left
  # the failure in place: every fit still aborted with "non-conformable arrays").
  # rstan::extract(permuted = TRUE) returns a SCALAR parameter (r_E, r_C, R_G,
  # R_T) as a 1-D ARRAY -- it carries a length-1 `dim` attribute, it is not a
  # bare vector. Single-bracket indexing of a 1-D array PRESERVES that dim
  # (dim(ex$R_G[use]) is length(use), not NULL), unlike a plain vector (no dim)
  # or an [iter,1] matrix (single-index drops the dim). The predictive means
  # below then evaluate `lamE[, day, drop=FALSE] * RG`, multiplying a 2-D array
  # by a 1-D array: both operands have a `dim`, the dims differ, and R throws
  # "non-conformable arrays". This fired on the gear branch (shore fits, R_G)
  # and the trailer branch (boat fit, R_T) -- exactly the fits that failed --
  # while the catch branch never did, because sweep() rebuilds STATS to match
  # dim(x) and cannot raise this error. as.numeric() strips the stray dim; the
  # multiply then recycles column-wise with the correct per-draw scaling
  # (verified: element [i, j] = lamE[i, day[j]] * scalar[i]). as.numeric() is a
  # no-op on an already-bare vector, so the fix is safe across rstan versions.
  # A two-line check on any fit: dim(rstan::extract(fit, "R_G")$R_G) returns the
  # iteration count (length-1 dim), not NULL.
  rE <- as.numeric(ex$r_E[use]); rC <- as.numeric(ex$r_C[use])
  RG <- as.numeric(ex$R_G[use])
  # RT is the per-draw MULTIPLIER m with mu_trailer = lambda_E * m:
  #   pooled        m = R_T
  #   gear-resolved m = 1 / R_G_boat
  # NULL when the fit declares no trailer expansion parameter.
  RT <- bss_trailer_multiplier(ex, trailer_par, use)

  # B1.7 fix: score each observation on its finite predictive draws only. An
  # extreme lambda draw (exp() overflow in a weakly-identified fit) yields a
  # non-finite mu, and rnbinom(mu = Inf) returns NA, which previously aborted
  # quantile() and the entire PPC. Non-finite mu and non-finite draws are now
  # dropped; an observation with < 20 usable draws is recorded NA and excluded.
  calib <- function(mu_mat, y, size_vec) {
    nobs <- length(y); cov50 <- cov95 <- pit <- rep(NA_real_, nobs)
    for (i in seq_len(nobs)) {
      mu_i <- pmax(mu_mat[, i], 1e-8)
      keep <- is.finite(mu_i) & is.finite(size_vec)
      if (sum(keep) < 20) next
      yp <- stats::rnbinom(sum(keep), mu = mu_i[keep], size = size_vec[keep])
      yp <- yp[is.finite(yp)]
      if (length(yp) < 20) next
      qq <- stats::quantile(yp, c(.025, .25, .75, .975), names = FALSE, na.rm = TRUE)
      cov50[i] <- y[i] >= qq[2] && y[i] <= qq[3]
      cov95[i] <- y[i] >= qq[1] && y[i] <= qq[4]
      pit[i]   <- mean(yp < y[i]) + 0.5 * mean(yp == y[i])
    }
    usable <- is.finite(pit)
    list(summary = data.frame(coverage_50 = mean(cov50, na.rm = TRUE),
                              coverage_95 = mean(cov95, na.rm = TRUE),
                              pit_mean = mean(pit, na.rm = TRUE),
                              pit_sd = stats::sd(pit, na.rm = TRUE),
                              n = sum(usable), n_obs = nobs),
         pit = pit[usable])
  }

  parts <- list()
  if (stan_data$Gear_n > 0)
    parts$gear    <- calib(lamE[, stan_data$day_Gear, drop = FALSE] * RG,
                           stan_data$Gear_I, rE)
  if (stan_data$T_n > 0 && !is.null(RT))
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
