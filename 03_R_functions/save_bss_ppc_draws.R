###############################################################################
# PERSIST THE POSTERIOR DRAWS THE PPC NEEDS, SO A DIAGNOSTIC FIX IS NOT A RE-FIT
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   Nothing in this pipeline persisted a stanfit or any subset of its draws. Every
#   posterior-predictive number lived only as the CSV the run happened to write, so a
#   defect discovered in the DIAGNOSTIC CODE after a run could only be corrected by
#   refitting the model. That has now happened three times:
#
#     2026-08-31  ppc_calibration scored coverage against a quantile interval of
#                 simulated draws, which over-covers small counts by construction and
#                 produced a phantom "trailer over-coverage" open item.
#     2026-09-01  the tau_bar row key did not match rstan's "tau_bar[1]", which silently
#                 dropped whole prior_vs_posterior files.
#     2026-09-03  pit_block() and calib() scored the catch stream as plain NB2 while the
#                 fit was ZINB, reporting 419 expected zeros against 676 observed where
#                 the mixture gives 637.9. The ZINB feature was very nearly rejected on it.
#
#   In each case the FIT was fine and only the arithmetic downstream of it was wrong. A
#   shore all-gear daily fit costs about 3.7 h, so the standing tax for a one-line
#   diagnostic fix is several hours of compute and a day of turnaround.
#
# WHAT IS SAVED, AND WHY NOT THE WHOLE STANFIT
#   A full stanfit for a 289-day daily fit is several GB, mostly the per-day latent
#   series nothing downstream reads. What every R-side predictive statistic needs is the
#   handful of draw objects the likelihood is built from. Those come to roughly
#   8000 draws x 289 days x 2 matrices plus a few scalar vectors, i.e. tens of MB
#   compressed, which is committable alongside the CSVs.
#
#   The saved object is deliberately SELF-DESCRIBING: it carries stan_data, the draw
#   subset index and the parameter names, so a recomputation can be checked against the
#   run it came from rather than assumed to match.
#
# WHAT IT IS NOT
#   It is not a substitute for the fit. Anything that needs the full posterior (loo on a
#   parameter not in this list, a new generated quantity, a different likelihood) still
#   needs a re-fit. It covers exactly the PPC and zero-bin family, which is where the
#   defects have actually been.
###############################################################################

# Draw objects the PPC and zero-bin diagnostics read. Anything absent from a given fit
# (theta_C when zi_catch = 0, R_T when the fit uses R_G_boat) is skipped silently.
.sbd_pars <- c("lambda_E_S", "lambda_C_S", "r_E", "r_C", "R_G", "R_T", "R_G_boat",
               "L_out", "kappa_OSP", "r_OSP", "theta_C")

save_bss_ppc_draws <- function(fit, stan_data, label, output_dir,
                               pars = .sbd_pars, max_draws = NULL, quiet = FALSE) {
  if (is.null(fit) || is.null(stan_data)) return(invisible(NULL))
  out <- tryCatch({
    have <- intersect(pars, fit@sim$pars_oi %||% pars)
    dr <- list()
    for (p in have) {
      x <- try(rstan::extract(fit, pars = p)[[p]], silent = TRUE)
      if (!inherits(x, "try-error") && !is.null(x)) dr[[p]] <- x
    }
    if (!length(dr)) return(invisible(NULL))
    nd <- dim(dr[[1]])[1]
    # Optional thinning. The PPC averages over draws, so a subset is an unbiased but
    # noisier estimate; the default keeps everything because the point of the file is to
    # reproduce the run's numbers EXACTLY, not approximately.
    idx <- if (!is.null(max_draws) && nd > max_draws)
             sort(sample.int(nd, max_draws)) else seq_len(nd)
    if (length(idx) < nd) dr <- lapply(dr, function(x)
      if (is.null(dim(x))) x[idx] else if (length(dim(x)) == 2) x[idx, , drop = FALSE]
      else x[idx, , , drop = FALSE])
    obj <- list(label = label, saved_at = Sys.time(), n_draws = length(idx),
                thinned = length(idx) < nd, pars = names(dr), draws = dr,
                stan_data = stan_data,
                note = paste("PPC draw subset. Recompute ppc_byobs / ppc_calibration from",
                             "this without refitting; see 03_R_functions/save_bss_ppc_draws.R."))
    f <- file.path(output_dir, sprintf("ppc_draws_%s.rds", label))
    saveRDS(obj, f, compress = "xz")
    if (!isTRUE(quiet))
      cat(sprintf("    ppc draws saved: %s (%d draws, %d pars, %.1f MB)\n",
                  basename(f), length(idx), length(dr), file.size(f) / 1e6))
    f
  }, error = function(e) { if (!isTRUE(quiet))
      cat(sprintf("    ppc draws NOT saved for %s: %s\n", label, conditionMessage(e))); NULL })
  invisible(out)
}
