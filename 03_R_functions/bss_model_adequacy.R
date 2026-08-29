# -----------------------------------------------------------------------------
# Part of Coastal-Rec-Crab-BSS: recreational Dungeness crab creel estimation
# for Grays Harbor / Westport (WDFW).
# Copyright (C) 2024-2026 Washington Department of Fish and Wildlife.
#
# Adapted from CreelEstimates, the WDFW freshwater creel estimation framework:
#   https://github.com/dfw-wa/CreelEstimates   (licensed GPL-3.0).
# Substantial portions of the methodology, structure, and R/Stan code originate
# in CreelEstimates and remain (C) their authors under GPL-3.0; changes for
# recreational crab are by WDFW.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License, version 3, as published by the Free
# Software Foundation. It is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program (see the LICENSE file);
# if not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
###############################################################################
# bss_model_adequacy.R  --  improvement-plan item 5.5.
#
# WHY THIS EXISTS
#
# The convergence gate answers ONE question: did this fit sample? R-hat, n_eff, the
# divergence fraction and the SD-normalized divergence impact are all properties of the
# SAMPLER. The project has been reading the answer as "is this model right?", and the
# 2026-08-29 improvement batch made the cost of that conflation measurable:
#
#   SIX configurations, all passing every gate criterion, spanning
#   25,883 to 37,392 on the private-boat all-gear component (44%)
#   and 66,237 to 78,250 on the port total (18%).
#
# For scale, the entire 2026-08-25 improvement batch moved the port total 1.6%. Nothing in
# the pipeline currently stops a future run adopting the daily boat AR on the strength of
# four green ticks, even though at that resolution the trailer stream's p_loo triples
# (4.9 -> 31.8 effective parameters on 195 observations), two Pareto k exceed 0.7, the CATCH
# stream's elpd gets WORSE, and the OSP observation-model overdispersion collapses (r_OSP
# 1.62 -> 10.61 with a 95% upper bound of 1,405, its parent sigma_r_OSP at n_eff 307).
#
# Every one of those quantities was already being computed and written to
# loo_summary_*.csv and ppc_calibration_*.csv. None of it reached the table a reader
# consults to decide whether to trust a fit. This assembles them into that table.
#
# WHAT IT DOES NOT DO
#
# It does not gate. `pass_convergence` keeps its exact meaning and keeps deciding PE vs BSS;
# these columns sit BESIDE it so the two questions are visibly different, and so a
# configuration comparison has something to be argued with. Turning any of this into a hard
# gate is a modelling decision that needs its own justification and its own run, not a
# default someone inherits.
#
# THE FLAGS
#   p_loo_frac        p_loo / n_obs, worst stream. A latent process with one state per day
#                     can drive this arbitrarily high; > 0.25 means the effective parameter
#                     count is a quarter of the data.
#   n_pareto_bad      count of Pareto k > 0.7 across streams. Above zero the LOO estimate
#                     itself is unreliable, so an elpd COMPARISON at that k is not evidence.
#   pit_worst_bias    max |pit_mean - 0.5| across streams. A systematic offset means the
#                     model predicts the wrong level for that stream, which is what a scale
#                     conflict between two effort streams looks like.
#   disp_neff_min     smallest n_eff among the observation-model dispersion parameters
#                     (r_E, r_C, r_OSP and their sigma_* parents). The gate checks n_eff on
#                     C_sum and E_sum only, so a dispersion parameter can sit below the
#                     gate's own floor without the gate noticing.
###############################################################################

# Read a per-fit diagnostic CSV written earlier in the same run. Returns NULL when absent,
# so this never depends on ordering and never aborts a run.
.bma_read <- function(output_dir, pattern, label) {
  f <- file.path(output_dir, sprintf(pattern, label))
  if (!file.exists(f)) return(NULL)
  tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
}

# The pure computation, shared by the live path (bss_model_adequacy, which has the stanfit)
# and the retro path (annotate_model_adequacy_run, which has only committed CSVs). Kept as
# one body on purpose: a reconstruction that computed its numbers slightly differently from
# the live writer would be worse than no reconstruction at all, because the two would be
# quoted side by side in the same table.
#
# `disp_tbl` is a data.frame of the NON-decoupled dispersion parameters with columns
# `parameter` and `n_eff`; each caller builds it from whatever it has.
.bma_core <- function(label, loo, ppc, disp_tbl,
                      neff_floor = 400, p_loo_frac_warn = 0.25, pit_bias_warn = 0.05,
                      source_tag = "live") {
  p_loo_frac <- NA_real_; n_bad_k <- NA_integer_; worst_stream <- NA_character_
  if (!is.null(loo) && all(c("p_loo", "n_obs") %in% names(loo))) {
    fr <- suppressWarnings(as.numeric(loo$p_loo) / pmax(as.numeric(loo$n_obs), 1))
    if (any(is.finite(fr))) {
      i <- which.max(replace(fr, !is.finite(fr), -Inf))
      p_loo_frac <- round(fr[i], 4)
      worst_stream <- as.character(loo$stream[i])
    }
    if ("n_pareto_k_gt_0.7" %in% names(loo))
      n_bad_k <- sum(suppressWarnings(as.integer(loo[["n_pareto_k_gt_0.7"]])), na.rm = TRUE)
  }

  pit_worst <- NA_real_; pit_stream <- NA_character_
  if (!is.null(ppc) && "pit_mean" %in% names(ppc)) {
    b <- abs(suppressWarnings(as.numeric(ppc$pit_mean)) - 0.5)
    if (any(is.finite(b))) {
      i <- which.max(replace(b, !is.finite(b), -Inf))
      pit_worst <- round(b[i], 4)
      pit_stream <- as.character(ppc$data_type[i])
    }
  }

  disp_neff <- NA_real_; disp_par <- NA_character_
  if (!is.null(disp_tbl) && nrow(disp_tbl)) {
    ne <- suppressWarnings(as.numeric(disp_tbl$n_eff))
    if (any(is.finite(ne))) {
      i <- which.min(replace(ne, !is.finite(ne), Inf))
      disp_neff <- round(ne[i]); disp_par <- as.character(disp_tbl$parameter[i])
    }
  }

  data.frame(
    fit = label,
    p_loo_frac = p_loo_frac,
    p_loo_worst_stream = worst_stream,
    n_pareto_bad = n_bad_k,
    pit_worst_bias = pit_worst,
    pit_worst_stream = pit_stream,
    disp_neff_min = disp_neff,
    disp_neff_min_par = disp_par,
    flag_overparameterised = isTRUE(p_loo_frac > p_loo_frac_warn),
    flag_loo_unreliable    = isTRUE(n_bad_k > 0),
    flag_pit_bias          = isTRUE(pit_worst > pit_bias_warn),
    flag_dispersion_neff   = isTRUE(disp_neff < neff_floor),
    source = source_tag,
    stringsAsFactors = FALSE)
}

# THE DISPERSION PARAMETERS THE GATE DOES NOT LOOK AT.
.BMA_DISP_PARS <- c("r_E", "r_C", "r_OSP", "sigma_r_E", "sigma_r_C", "sigma_r_OSP")

# Adequacy row for one fit, live, from the stanfit. Every field is NA when its source file
# is absent, so this never aborts a run that has already spent hours.
bss_model_adequacy <- function(fit, stan_data, label, output_dir,
                               neff_floor = 400,
                               p_loo_frac_warn = 0.25,
                               pit_bias_warn = 0.05) {
  # Decoupled dispersion parameters are excluded: an unused r_OSP in a shore fit has no
  # n_eff worth reporting and would fire the flag on every run.
  have <- .BMA_DISP_PARS[.BMA_DISP_PARS %in% fit@model_pars]
  if (length(have)) {
    dec  <- bss_decoupled_reasons(have, stan_data)
    have <- have[is.na(dec)]
  }
  disp_tbl <- NULL
  if (length(have)) {
    sm <- tryCatch(rstan::summary(fit, pars = have)$summary, error = function(e) NULL)
    if (!is.null(sm) && "n_eff" %in% colnames(sm))
      disp_tbl <- data.frame(parameter = rownames(sm), n_eff = sm[, "n_eff"],
                             stringsAsFactors = FALSE)
  }
  .bma_core(label,
            loo = .bma_read(output_dir, "loo_summary_%s.csv", label),
            ppc = .bma_read(output_dir, "ppc_calibration_%s.csv", label),
            disp_tbl = disp_tbl,
            neff_floor = neff_floor, p_loo_frac_warn = p_loo_frac_warn,
            pit_bias_warn = pit_bias_warn, source_tag = "live")
}

# ---------------------------------------------------------------------------
# RETRO PATH. Rebuild the adequacy table for a run folder that predates this file, from the
# CSVs it already committed. No stanfit, no refit, nothing rewritten.
#
# WHY IT EARNS ITS PLACE. Plan item 5.3 asks for the boat estimate to be attributed to one
# mechanism "with p_loo, r_OSP and the two PIT means quoted for each cell" of a 2x2. Two of
# those cells are ARCHIVED runs from the 2026-08-29 batch, fitted before this diagnostic
# existed. Without a reconstruction the comparison table would carry adequacy for the new
# cells and blanks for the old ones, which is the shape of an argument that quietly favours
# whichever cells happen to be new.
#
# The dispersion n_eff comes from structural_params_*.csv, which already carries both n_eff
# and the `decoupled` flag, so the same exclusion rule applies without needing the stanfit.
# Runs written before the decoupled column existed have no flag; those rows are kept and the
# file records source = "reconstructed (no decoupled column)" so the difference is visible
# rather than assumed away.
# ---------------------------------------------------------------------------
annotate_model_adequacy_run <- function(dir, overwrite = FALSE, quiet = FALSE) {
  if (!dir.exists(dir)) { if (!quiet) cat("  skip (no such folder):", dir, "\n"); return(invisible(NULL)) }
  out <- file.path(dir, "model_adequacy_reconstructed.csv")
  if (file.exists(out) && !isTRUE(overwrite)) {
    if (!quiet) cat("  skip (already annotated):", basename(dir), "\n")
    return(invisible(utils::read.csv(out, stringsAsFactors = FALSE)))
  }
  labs <- sub("^loo_summary_", "", sub("\\.csv$", "",
              list.files(dir, pattern = "^loo_summary_.*\\.csv$")))
  if (!length(labs)) { if (!quiet) cat("  skip (no loo_summary_*.csv):", basename(dir), "\n")
                       return(invisible(NULL)) }
  rows <- list()
  for (lab in labs) {
    sp <- .bma_read(dir, "structural_params_%s.csv", lab)
    tag <- "reconstructed"
    disp_tbl <- NULL
    if (!is.null(sp) && "parameter" %in% names(sp) && "n_eff" %in% names(sp)) {
      keep <- sp$parameter %in% .BMA_DISP_PARS
      if ("decoupled" %in% names(sp)) {
        # `%in% TRUE` and not `!sp$decoupled`: a row whose flag is NA is UNKNOWN, not
        # known-decoupled, and dropping it would silently shrink the set being minimised over.
        keep <- keep & !(as.logical(sp$decoupled) %in% TRUE)
      } else {
        tag <- "reconstructed (no decoupled column)"
      }
      if (any(keep))
        disp_tbl <- data.frame(parameter = as.character(sp$parameter[keep]),
                               n_eff = suppressWarnings(as.numeric(sp$n_eff[keep])),
                               stringsAsFactors = FALSE)
    }
    r <- tryCatch(.bma_core(lab,
                            loo = .bma_read(dir, "loo_summary_%s.csv", lab),
                            ppc = .bma_read(dir, "ppc_calibration_%s.csv", lab),
                            disp_tbl = disp_tbl, source_tag = tag),
                  error = function(e) NULL)
    if (!is.null(r)) rows[[length(rows) + 1]] <- r
  }
  if (!length(rows)) return(invisible(NULL))
  df <- do.call(rbind, rows)
  # Method selection is recoverable from the committed convergence report.
  cr <- tryCatch(utils::read.csv(file.path(dir, "convergence_report.csv"),
                                 stringsAsFactors = FALSE), error = function(e) NULL)
  df$method_selected <- if (!is.null(cr) && all(c("fit", "method_selected") %in% names(cr)))
    cr$method_selected[match(df$fit, cr$fit)] else NA_character_
  utils::write.csv(df, out, row.names = FALSE)
  if (!quiet)
    cat(sprintf("  %-46s %d fit(s) -> model_adequacy_reconstructed.csv\n", basename(dir), nrow(df)))
  invisible(df)
}

# Assemble the adequacy table for every fitted component and write it beside the
# convergence report. Returns the data frame invisibly.
#
# `bss_all` is the driver's per-fit list; entries with no fit are skipped.
write_model_adequacy <- function(bss_all, output_dir, params = list()) {
  if (!length(bss_all) || is.null(output_dir)) return(invisible(NULL))
  rows <- list()
  for (label in names(bss_all)) {
    b <- bss_all[[label]]
    if (is.null(b$fit)) next
    r <- tryCatch(bss_model_adequacy(b$fit, b$bss_data, label, output_dir,
                                     neff_floor = params$min_n_eff %||% 400),
                  error = function(e) NULL)
    if (is.null(r)) next
    # The gate is the single authority on method selection, so quote ITS wording rather
    # than inventing a second one: a fit that reads "PE (convergence fail)" in
    # convergence_report.csv must not read "PE (gate fail)" here.
    r$method_selected <- b$gate_info$method_selected %||%
      (if (isTRUE(b$use_bss)) "BSS" else "PE (convergence fail)")
    rows[[length(rows) + 1]] <- r
  }
  if (!length(rows)) return(invisible(NULL))
  df <- do.call(rbind, rows)
  utils::write.csv(df, file.path(output_dir, "model_adequacy.csv"), row.names = FALSE)

  flagged <- df[df$flag_overparameterised | df$flag_loo_unreliable |
                df$flag_pit_bias | df$flag_dispersion_neff, , drop = FALSE]
  if (nrow(flagged)) {
    cat("\n  MODEL ADEQUACY (separate from the convergence gate; nothing here changes",
        "PE-vs-BSS selection):\n")
    for (i in seq_len(nrow(flagged))) {
      r <- flagged[i, ]
      msgs <- c(
        if (isTRUE(r$flag_overparameterised))
          sprintf("p_loo is %.0f%% of n_obs on the %s stream", 100 * r$p_loo_frac, r$p_loo_worst_stream),
        if (isTRUE(r$flag_loo_unreliable))
          sprintf("%d Pareto k > 0.7, so the LOO estimate is unreliable", r$n_pareto_bad),
        if (isTRUE(r$flag_pit_bias))
          sprintf("PIT mean is %.3f off nominal on the %s stream", r$pit_worst_bias, r$pit_worst_stream),
        if (isTRUE(r$flag_dispersion_neff))
          sprintf("%s has n_eff %.0f, below the gate's own floor", r$disp_neff_min_par, r$disp_neff_min))
      cat(sprintf("    %s [%s]: %s\n", r$fit, r$method_selected, paste(msgs, collapse = "; ")))
    }
    cat("    A fit can pass the gate and still be flagged here. The gate asks whether the\n")
    cat("    sampler worked; these ask whether the model is carrying the data.\n")
  }
  invisible(df)
}
