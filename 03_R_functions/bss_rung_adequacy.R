###############################################################################
# PER-RUNG ADEQUACY, COMPUTED INSIDE THE LADDER LOOP
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   The 2026-09-04 ladder fitted shore all-gear at weekly, biweekly and monthly and
#   all three passed the convergence gate. The gate therefore decided nothing: under
#   the production rule ("report the finest rung that PASSES") the answer is simply
#   the finest rung RUN, whatever its adequacy. The statistics that actually separate
#   the rungs are p_loo as a fraction of n_obs, the Pareto k > 0.7 count, and the
#   posterior-predictive coverage, and on that run they were available for exactly ONE
#   rung, because only the fit the loop KEEPS is passed to write_bss_diagnostics().
#
#   So the ladder answered "does daily overfit" (emphatically yes: p_loo 35.2% of
#   n_obs and 41 bad Pareto k at daily against 9.6% and 1 at weekly) and could not
#   answer "and where does it stop being too fine", which is the question a ladder is
#   for. Biweekly and monthly were fitted, cost 174 minutes, and left no evidence.
#
# WHAT THIS COMPUTES, AND WHAT IT COSTS
#   A compact adequacy row per rung, appended to ar_escalation_log.csv: p_loo and its
#   fraction of n_obs, the bad-k count, and coverage_50 per stream. It re-uses the
#   log_lik the model already wrote (loo is seconds to a minute on these sizes) and
#   the same randomized-PIT coverage the aggregate PPC uses, via model_diagnostics.R,
#   so the numbers are directly comparable with model_adequacy.csv for the kept rung.
#   Everything is wrapped: a rung whose adequacy cannot be computed logs NA and the
#   ladder continues, because losing a multi-hour fit to a diagnostic error is the
#   trade this project has already made twice and does not want again.
#
# WHAT IT DELIBERATELY DOES NOT DO
#   It does not gate. The convergence gate decides PE vs BSS and these numbers sit
#   beside it, exactly as model_adequacy.csv does. A rung that fails the gate is not
#   reportable whatever its p_loo; a rung that passes the gate with p_loo at a third
#   of n_obs is reportable and probably wrong, and the point of writing both down is
#   that a reader can see the difference.
###############################################################################

bss_rung_adequacy <- function(fit, stan_data, quiet = TRUE) {
  na_row <- list(p_loo = NA_real_, p_loo_frac = NA_real_, p_loo_worst_stream = NA_character_,
                 n_obs_loo = NA_integer_, n_pareto_bad = NA_integer_,
                 cov50_gear = NA_real_, cov50_catch = NA_real_, pit_mean_catch = NA_real_)
  if (is.null(fit) || is.null(stan_data)) return(na_row)
  # 2026-09-07 DEFECT FIX. bss_ppc_calibration() calls set.seed() unconditionally, and both
  # it and save_run_diagnostics() take a RANDOM SUBSET of draws. Calling this helper inside
  # the fitting loop therefore shifted the global RNG stream for every diagnostic computed
  # afterwards: on the 2026-09-04 run, C2's kept fit is BIT-IDENTICAL to L1's across 10,251
  # parameter rows and yet reported gear coverage_50 0.5659 against L1's 0.5595, purely
  # because a diagnostic-only toggle moved the seed. Bit-comparability is the tool that has
  # caught four defects in this project; a helper that quietly breaks it is worse than the
  # gap it fills. Save the RNG state and put it back.
  .seed_before <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (!is.null(.seed_before)) assign(".Random.seed", .seed_before, envir = globalenv())
    else suppressWarnings(rm(".Random.seed", envir = globalenv()))
  }, add = TRUE)
  out <- tryCatch({
    r <- na_row
    # --- loo, summed over streams EXACTLY as write_loo_diagnostics() does ------
    # Same stream list, same absent-stream guard, same loo::loo() call on the plain
    # draws-by-observation matrix (no r_eff). Matching it is the point: these numbers
    # have to be comparable with loo_summary_*.csv and model_adequacy.csv for the rung
    # the ladder keeps, or the ladder's own rows cannot be read against them.
    streams <- list(
      gear    = list(par = "log_lik_gear",    n = stan_data$Gear_n %||% 0),
      trailer = list(par = "log_lik_trailer", n = stan_data$T_n    %||% 0),
      catch   = list(par = "log_lik_catch",   n = stan_data$IntC   %||% 0))
    p_loo <- 0; n_obs <- 0L; bad <- 0L; any_ok <- FALSE
    frac <- setNames(rep(NA_real_, length(streams)), names(streams))
    if (requireNamespace("loo", quietly = TRUE)) {
      for (sn in names(streams)) {
        st <- streams[[sn]]
        if ((st$n %||% 0) == 0) next
        ll <- try(as.matrix(rstan::extract(fit, pars = st$par)[[1]]), silent = TRUE)
        if (inherits(ll, "try-error") || is.null(ll) || nrow(ll) == 0 || ncol(ll) == 0) next
        lo <- try(suppressWarnings(loo::loo(ll)), silent = TRUE)
        if (inherits(lo, "try-error") || is.null(lo)) next
        pl <- lo$estimates["p_loo", "Estimate"]
        frac[sn] <- pl / ncol(ll)
        p_loo <- p_loo + pl
        n_obs <- n_obs + ncol(ll)
        bad   <- bad + sum(lo$diagnostics$pareto_k > 0.7, na.rm = TRUE)
        any_ok <- TRUE
      }
    }
    if (any_ok && n_obs > 0) {
      # 2026-09-07 DEFECT FIX. p_loo_frac is the WORST STREAM's p_loo / n_obs, which is
      # what bss_model_adequacy.R reports and therefore the only definition that can be
      # read against model_adequacy.csv and against the daily baseline. The first version
      # summed every stream and divided by the summed n_obs, which on the 2026-09-04 C2
      # ladder printed 2.8% for a rung whose model_adequacy figure is 9.6%: two different
      # statistics under one name, in the same run, being compared with each other.
      w <- which.max(replace(frac, is.na(frac), -Inf))
      r$p_loo <- p_loo; r$n_obs_loo <- n_obs; r$n_pareto_bad <- bad
      if (is.finite(frac[w])) {
        r$p_loo_frac <- unname(frac[w]); r$p_loo_worst_stream <- names(frac)[w]
      }
    }
    # --- coverage ---------------------------------------------------------------
    # This is the ppc_calibration_*.csv statistic: the exact randomized PIT computed on a
    # 400-draw subset. NOTE it is NOT the number model_adequacy.csv reports, which prefers
    # ppc_byobs (a different draw subset) when that file exists and can differ by a few
    # thousandths. Compare these RUNG TO RUNG, where the statistic is identical throughout;
    # do not read a single rung's value against model_adequacy.csv.
    ppc <- try(bss_ppc_calibration(fit, stan_data), silent = TRUE)
    if (!inherits(ppc, "try-error") && !is.null(ppc$summary)) {
      sm <- as.data.frame(ppc$summary, stringsAsFactors = FALSE)
      pick <- function(st, col) {
        v <- suppressWarnings(as.numeric(sm[[col]][sm$data_type == st]))
        if (length(v)) v[1] else NA_real_
      }
      r$cov50_gear     <- pick("gear",  "coverage_50")
      r$cov50_catch    <- pick("catch", "coverage_50")
      r$pit_mean_catch <- pick("catch", "pit_mean")
    }
    r
  }, error = function(e) {
    if (!isTRUE(quiet)) cat("    rung adequacy unavailable:", conditionMessage(e), "\n")
    na_row
  })
  out
}
