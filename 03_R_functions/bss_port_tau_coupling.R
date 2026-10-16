###############################################################################
# bss_port_tau_coupling.R -- 2026-10-16 (CHANGE_REGISTER A38; closes D44 by Matt's option (i))
#
# THE PORT TOTAL PAIRS THE COMPONENT FITS' DRAWS ROW BY ROW. Two fits of one population
# (the pot-closure and all-gear sub-seasons, or the sub-seasons of a span) draw their
# turnover tau_bar around ONE prior centre: on the shore the derived I/E centre (2.477 on
# 2024-25, log-SD floored at 0.10), whose uncertainty is one number for both fits, and
# both shore fits are prior-dominated on it (B65: 2.38 [1.98, 2.87] and 2.48 [2.03, 3.02]).
# Pairing the fits' draws independently sums two draws of that one uncertainty as if they
# were independent, which understates the port interval (D44: on the authoritative run the
# bound is about +5% on the port SD; the mean of the sum is unchanged exactly, the median
# moves only by the components' asymmetry).
#
# bss_couple_component_draws() re-pairs the draws of the listed populations' components
# COMONOTONICALLY in tau_bar. Within a population and catch group, the component with the
# largest median expected catch is the REFERENCE and keeps its draws in their order; every
# other component is re-indexed so that row j carries the draw whose tau_bar quantile is the
# reference's at row j. Each component's own joint posterior is untouched (a row is still one
# of its draws, every draw vector moved together); only the cross-component pairing changes,
# and only within the population. Between populations the pairing is as before. It is the
# bound a fully shared centre implies, chosen over a shared factor or a joint fit on
# 2026-10-16 (CHANGE_REGISTER D44). port_tau_coupling = "none" returns the input unchanged,
# bit for bit, with no attribute.
#
# Draw counts may differ between fits (the drivers sample the pot-closure fits shorter): the
# re-indexed component gets the reference's length n_ref, row j taking the component's draw
# at tau-quantile rank ceiling(rank_ref[j] * n_k / n_ref), a quantile-matched resample of its
# own draws. The drivers' rep_len() then recycles the reference and the re-indexed component
# identically to n_draws_max. The index vector is kept on the component as draw_perm (into the
# ORIGINAL fit draws) so the monthly block, which re-reads per-day draws from the fit object,
# applies the same pairing; the original expected-catch draws are kept as C_exp_draws_indep so
# the driver can report the port interval under both pairings (port_tau_coupling.csv).
#
# A component is eligible when it reports BSS (bss_use_pe_for() is FALSE) and carries
# tau_draws of the same length as its C_exp_draws (the drivers extract tau_bar beside the
# totals). A population with fewer than two eligible components is left as it is and the
# fact is recorded in the attribute. Pure: reads params$port_tau_coupling ("comonotone" |
# "none") and params$port_tau_coupling_populations (character; "shore" ships), writes
# nothing; the console line is suppressed with quiet = TRUE.
###############################################################################
bss_couple_component_draws <- function(bss_all, params = list(), quiet = FALSE) {
  mode <- params$port_tau_coupling %||% "none"
  if (!is.character(mode) || length(mode) != 1L || !mode %in% c("comonotone", "none"))
    stop(sprintf("bss_couple_component_draws(): port_tau_coupling = '%s' is not one of comonotone | none.",
                 paste(mode, collapse = ",")), call. = FALSE)
  if (identical(mode, "none") || !length(bss_all)) return(bss_all)
  pops <- as.character(params$port_tau_coupling_populations %||% "shore")
  draw_names <- c("C_draws", "C_exp_draws", "C_pred_draws", "E_draws", "tau_draws")

  is_bss <- function(b) !is.null(b) && is.list(b) && !isTRUE(bss_use_pe_for(b))
  eligible <- function(b) is_bss(b) && is.numeric(b$tau_draws) && is.numeric(b$C_exp_draws) &&
    length(b$tau_draws) == length(b$C_exp_draws) && length(b$tau_draws) >= 2L && all(is.finite(b$tau_draws))

  out <- bss_all
  rows <- list()
  for (pop in pops) {
    keys <- names(out)[vapply(out, function(b) is.list(b) && identical(b$population %||% "", pop), logical(1))]
    ok <- keys[vapply(out[keys], eligible, logical(1))]
    if (length(ok) < 2L) {
      rows[[length(rows) + 1L]] <- data.frame(population = pop, catch_group = NA_character_, component = NA_character_,
        role = sprintf("no coupling: %d eligible BSS component(s) with tau_bar draws", length(ok)),
        reference = NA_character_, n_draws = NA_integer_, n_ref = NA_integer_, tau_median = NA_real_, stringsAsFactors = FALSE)
      next
    }
    cgs <- unique(vapply(out[ok], function(b) b$catch_group %||% "", character(1)))
    for (cg in cgs) {
      k_cg <- ok[vapply(out[ok], function(b) identical(b$catch_group %||% "", cg), logical(1))]
      if (length(k_cg) < 2L) next
      med <- vapply(out[k_cg], function(b) stats::median(b$C_exp_draws), numeric(1))
      ref <- k_cg[which.max(med)]
      tau_ref <- out[[ref]]$tau_draws; n_ref <- length(tau_ref)
      r_ref <- rank(tau_ref, ties.method = "first")
      rows[[length(rows) + 1L]] <- data.frame(population = pop, catch_group = cg, component = ref, role = "reference (keeps its draws)",
        reference = ref, n_draws = n_ref, n_ref = n_ref, tau_median = stats::median(tau_ref), stringsAsFactors = FALSE)
      for (k in setdiff(k_cg, ref)) {
        b <- out[[k]]; n_k <- length(b$tau_draws)
        ord_k <- order(b$tau_draws, method = "radix")          # the component's draws, ascending in tau_bar
        perm <- ord_k[pmax(1L, pmin(n_k, as.integer(ceiling(r_ref * n_k / n_ref))))]
        b$C_exp_draws_indep <- b$C_exp_draws
        b$draw_perm <- perm
        for (nm in draw_names) if (is.numeric(b[[nm]]) && length(b[[nm]]) == n_k) b[[nm]] <- b[[nm]][perm]
        if (is.matrix(b$C_gear_draws) && nrow(b$C_gear_draws) == n_k) b$C_gear_draws <- b$C_gear_draws[perm, , drop = FALSE]
        out[[k]] <- b
        rows[[length(rows) + 1L]] <- data.frame(population = pop, catch_group = cg, component = k, role = "re-indexed to the reference's tau_bar quantiles",
          reference = ref, n_draws = n_k, n_ref = n_ref, tau_median = stats::median(b$tau_draws), stringsAsFactors = FALSE)
      }
    }
  }
  tab <- do.call(rbind, rows)
  attr(out, "tau_coupling") <- tab
  if (!isTRUE(quiet)) {
    done <- tab[grepl("^re-indexed", tab$role), , drop = FALSE]
    if (nrow(done)) cat(sprintf("  Port total pairing (A38, port_tau_coupling = comonotone): %s re-indexed to %s's tau_bar quantiles.\n",
                                paste(done$component, collapse = ", "), paste(unique(done$reference), collapse = ", ")))
    else cat("  Port total pairing (A38, port_tau_coupling = comonotone): no population had two BSS components with tau_bar draws; nothing re-indexed.\n")
  }
  out
}
