###############################################################################
# diagnose_tau_boat_sensitivity.R
#
# Single-run prior-sensitivity diagnostic for the boat deployment-turnover
# parameter tau_boat (the gear-resolved boat's `L`). Answers: "how much does the
# boat catch, and the port total, move if the tau_boat prior is wrong?" WITHOUT
# refitting the model. Toggle with run_config$diagnose_tau_sensitivity.
#
# WHY A PROJECTION IS VALID HERE (and when it is not)
#   In crab_bss_gear_resolved.stan the boat season catch is
#       C_sum = sum_d ( lambda_E[d] * E_scale * L[d] * lambda_C[d] ),
#   with, for the boat, E_scale = 1 and L[d] = tau_boat_prior_mu * exp(sigma * L_raw[d]).
#   lambda_E is identified by the trailer counts (T_I ~ NB2(lambda_E / R_G_boat, r_E))
#   and lambda_C by the interview catch (c ~ NB2(lambda_C * number_of_gear, r_C));
#   NEITHER depends on tau. So scaling the tau_boat prior mean by a factor f scales
#   every L[d] by f and therefore scales C_sum by f, draw for draw. The catch is
#   proportional to tau_boat (elasticity 1) *as long as tau is prior-dominated*,
#   i.e. as long as the boat ingress/egress (I/E) stream carries no in-window days
#   to pin L against the effort. When boat I/E days ARE present, L becomes partly
#   data-identified and this proportional projection OVERSTATES the sensitivity;
#   in that case run the exact multi-refit check in 06_diagnostics/run_tau_sweep.R.
#
#   The diagnostic reports whether tau was prior-dominated this run (posterior vs
#   prior spread on log L) so the projection's assumption is auditable, not hidden.
#
# INPUTS
#   boat          : the bss_all element for private_boat_all_gear (needs $fit for
#                   the prior-vs-posterior check; the projection itself needs only
#                   boat_C_draws). Pass the list element, not the stanfit.
#   port_C_draws  : the combined port-total catch draws (bss_C_total) for the group.
#   boat_C_draws  : the boat's aligned contribution to port_C_draws (same length,
#                   same draw order as port_C_draws).
#   params        : the merged run params (uses tau_boat_prior_mu, tau_boat_prior_sigma,
#                   tau_sensitivity_grid, ie_min_obs_boat).
#   n_ie_boat     : number of in-window boat I/E days used by the fit (0 = inert).
#   output_dir    : if set, writes tau_boat_sensitivity.csv there.
#
# RETURNS (invisibly) a list: $table (data.frame), $prior_dominated (logical/NA),
#   $note (character), $csv (path or NA). Never stops the run; on any error it
#   returns a note and NULL table so the driver can render a graceful message.
###############################################################################

diagnose_tau_boat_sensitivity <- function(boat,
                                          port_C_draws,
                                          boat_C_draws,
                                          params,
                                          n_ie_boat = 0L,
                                          output_dir = NULL) {

  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

  fail <- function(msg) invisible(list(table = NULL, prior_dominated = NA,
                                       note = msg, csv = NA_character_))

  tau_mu  <- suppressWarnings(as.numeric(params$tau_boat_prior_mu %||% 1.2))
  tau_sig <- suppressWarnings(as.numeric(params$tau_boat_prior_sigma %||% 0.3))
  grid    <- params$tau_sensitivity_grid %||% c(0.9, 1.0, 1.2, 1.5, 1.8)
  grid    <- sort(unique(suppressWarnings(as.numeric(grid))))
  grid    <- grid[is.finite(grid) & grid > 0]

  if (!is.finite(tau_mu) || tau_mu <= 0)      return(fail("tau_boat_prior_mu not usable; tau sensitivity skipped."))
  if (length(grid) == 0)                      return(fail("tau_sensitivity_grid empty; tau sensitivity skipped."))
  if (is.null(boat_C_draws) || is.null(port_C_draws) ||
      length(boat_C_draws) == 0 || length(port_C_draws) == 0)
    return(fail("Boat contributed no BSS draws (PE fallback or insufficient data); tau sensitivity not applicable."))
  if (length(boat_C_draws) != length(port_C_draws))
    return(fail("boat_C_draws and port_C_draws lengths differ; tau sensitivity skipped (driver wiring)."))

  # --- Prior-vs-posterior check on L (is tau prior-dominated this run?) --------
  prior_dominated <- NA
  pd_note <- "prior-vs-posterior on L not evaluated (no fit passed)."
  fit <- tryCatch(boat$fit, error = function(e) NULL)
  if (!is.null(fit)) {
    Lsumm <- tryCatch({
      Ld <- rstan::extract(fit, "L")$L            # [draws, D]
      # log-scale spread of the day-mean turnover across draws vs the prior sigma
      day_mean <- rowMeans(Ld)                    # per-draw mean turnover
      post_sd_log <- stats::sd(log(day_mean[day_mean > 0]))
      list(post_mean = mean(day_mean),
           ratio_post_prior_sd = post_sd_log / tau_sig)
    }, error = function(e) NULL)
    if (!is.null(Lsumm)) {
      # If the posterior retains most of the prior's log-spread, the data barely
      # moved tau -> prior-dominated -> the proportional projection is exact.
      prior_dominated <- isTRUE(Lsumm$ratio_post_prior_sd > 0.80) &&
                         (n_ie_boat < (params$ie_min_obs_boat %||% 2L))
      pd_note <- sprintf(
        "L posterior day-mean %.2f vs prior mean %.2f; posterior/prior log-SD ratio %.2f; boat I/E days = %d.",
        Lsumm$post_mean, tau_mu, Lsumm$ratio_post_prior_sd, n_ie_boat)
    }
  } else if (n_ie_boat < (params$ie_min_obs_boat %||% 2L)) {
    prior_dominated <- TRUE
    pd_note <- sprintf("No fit passed for the L check, but boat I/E days = %d (< min), so tau is prior-driven.", n_ie_boat)
  }

  # --- Project boat & port totals across the grid -----------------------------
  base_boat <- stats::median(boat_C_draws)
  base_port <- stats::median(port_C_draws)
  q <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)

  rows <- lapply(grid, function(tp) {
    f        <- tp / tau_mu
    boat_p   <- boat_C_draws * f
    port_p   <- port_C_draws + boat_C_draws * (f - 1)   # swap only the boat term
    bq <- q(boat_p); pq <- q(port_p)
    data.frame(
      tau_boat_prior_mu = tp,
      scale_factor      = round(f, 4),
      boat_catch_median = round(bq[2]),
      boat_lo95         = round(bq[1]),
      boat_hi95         = round(bq[3]),
      boat_pct_change   = round(100 * (bq[2] - base_boat) / base_boat, 1),
      port_total_median = round(pq[2]),
      port_lo95         = round(pq[1]),
      port_hi95         = round(pq[3]),
      port_pct_change   = round(100 * (pq[2] - base_port) / base_port, 1),
      is_current_prior  = isTRUE(abs(tp - tau_mu) < 1e-8),
      stringsAsFactors  = FALSE)
  })
  tbl <- do.call(rbind, rows)

  csv_path <- NA_character_
  if (!is.null(output_dir)) {
    csv_path <- file.path(output_dir, "tau_boat_sensitivity.csv")
    tryCatch(utils::write.csv(tbl, csv_path, row.names = FALSE),
             error = function(e) { csv_path <<- NA_character_ })
  }

  # Elasticity of the port total to tau across the grid (should be < the boat's 1.0,
  # diluted by the fixed shore + census components).
  span <- range(grid)
  port_span_pct <- 100 * (tbl$port_total_median[which.max(grid == span[2])] -
                          tbl$port_total_median[which.max(grid == span[1])]) / base_port

  note <- paste0(
    sprintf("tau_boat prior mean = %.2f (sigma %.2f). Boat catch is proportional to tau_boat under prior dominance; ",
            tau_mu, tau_sig),
    if (isTRUE(prior_dominated)) "this run IS prior-dominated, so the projection is ~exact. "
    else if (isFALSE(prior_dominated)) "this run is NOT prior-dominated (boat I/E informs tau), so treat the projection as an UPPER bound and confirm with run_tau_sweep.R. "
    else "prior dominance could not be confirmed; treat the projection as approximate. ",
    sprintf("Across tau_boat in [%.2f, %.2f] the port total spans ~%.0f%% of its base. ", span[1], span[2], abs(port_span_pct)),
    pd_note)

  invisible(list(table = tbl, prior_dominated = prior_dominated,
                 note = note, csv = csv_path))
}
