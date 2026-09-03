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
# annotate_decoupled_run.R  --  retro-fit the decoupled flags onto an OLD run folder.
#
# WHY THIS EXISTS
#
# `structural_params_*.csv` has carried `decoupled` / `decoupled_reason` since 2026-08-27
# and an `estimate` column (NA when decoupled) since 2026-08-30. Every run folder written
# BEFORE those dates has neither, and those folders are not scratch: they are the
# authoritative baselines the whole project compares against, and they are committed to git
# on purpose so past estimates are preserved as produced.
#
# So every decoupling trap the 2026-08-29 review documented sits UNMARKED in exactly the
# folders a reader is most likely to open:
#   * `kappa_OSP` reading a tidy ~3.0 in every fit of every run, which is its
#     lognormal(log 3, 0.3) prior, because the production `osp_scale_is_tau = TRUE` makes
#     the OSP mean use L instead. Read as an estimate that is a claim the model measured the
#     boat turnover at 3.0. It did not; it was told 3.0.
#   * shore `r_OSP` with a posterior MEAN of 2.5 million and a median that is no better.
#   * shore `R_G_boat` printing "4.02 gear/group (1.51-10.54)" beside the boat fits' real 3.55.
#   * `sigma_IE` on any fit whose I/E stream the ie_min_obs guard dropped.
#
# The run outputs must NOT be rewritten - a preserved estimate that gets edited later is
# worse than one that needs a companion file. This writes a NEW file, `decoupled_audit.csv`,
# alongside them, and leaves everything else untouched.
#
# WHAT IT CAN AND CANNOT RECONSTRUCT
#
# The Stan data list is not saved with a run, so the flags are rebuilt from what is: the fit
# label (which gives the population), `fit_data_summary.csv` (`n_ie_obs`), and
# `run_parameters.txt` (the `str()` dump of the run config, which carries
# `osp_scale_is_tau`, `estimate_cpue_density`, `use_crab_fraction`, `use_osp_crab_lower`,
# `opener_covariate_mode` and `shared_tau`). That covers every case in the list above.
#
# It CANNOT reconstruct the two window-dependent rules that `bss_decoupled_reasons()` applies
# from the Stan data directly - whether a fit's window contained a weekend day (B1, B1_C) or
# a holiday (B2, B2_C). Those are reported as `unknown` rather than guessed. Every row also
# carries `source = "reconstructed"` so an audit file is never mistaken for the real thing.
#
# USAGE
#   annotate_decoupled_run("05_output/20260826/pooled-CPUE-PV4-minint")
#   invisible(lapply(list.dirs("05_output", recursive = TRUE), annotate_decoupled_run))
###############################################################################

# Was this dump written by a driver that truncated str()? Every run_parameters.txt written
# before 2026-09-04 was capped at str()'s list.len = 99 while run_config carried 120 keys, so
# roughly 21 keys are simply absent and every lookup for one of them returns NA. NA then
# flows into isTRUE() and reads as FALSE, which is indistinguishable from a genuine FALSE.
# On the 2026-09-01 and later folders `opener_covariate_mode` fell past the cut, so B_open
# would have been left UNFLAGGED in a run where it is decoupled: the same class of silent
# missing-flag as the unflagged R_G_boat raised on 2026-08-25. Callers get a warning rather
# than a quiet wrong answer; the driver-side annotation reads `params` directly and is not
# affected, so this only matters for post-hoc backfill of old folders.
.adr_truncated <- function(txt) any(grepl("list output truncated", txt, fixed = TRUE))

# Pull one scalar out of a run_parameters.txt str() dump. Returns NA when absent.
.adr_param <- function(txt, key) {
  hit <- grep(sprintf("^\\s*\\$ %s\\s*:", key), txt, value = TRUE)
  if (!length(hit)) return(NA)
  v <- trimws(sub("^.*:", "", hit[1]))
  if (grepl("^logi", v)) return(as.logical(trimws(sub("^logi", "", v))))
  if (grepl("^(num|int)", v)) return(suppressWarnings(as.numeric(trimws(sub("^(num|int)", "", v)))))
  if (grepl("^chr", v)) return(trimws(gsub('"', "", sub("^chr", "", v))))
  NA
}

annotate_decoupled_run <- function(dir, overwrite = FALSE, quiet = FALSE) {
  sp <- list.files(dir, pattern = "^structural_params_.*\\.csv$", full.names = TRUE)
  if (!length(sp)) return(invisible(NULL))
  out_path <- file.path(dir, "decoupled_audit.csv")
  if (file.exists(out_path) && !isTRUE(overwrite)) {
    if (!isTRUE(quiet)) cat("  already annotated:", basename(dir), "\n")
    return(invisible(out_path))
  }

  rp   <- if (file.exists(file.path(dir, "run_parameters.txt")))
            readLines(file.path(dir, "run_parameters.txt"), warn = FALSE) else character(0)
  fds  <- if (file.exists(file.path(dir, "fit_data_summary.csv")))
            utils::read.csv(file.path(dir, "fit_data_summary.csv"), stringsAsFactors = FALSE) else NULL

  if (.adr_truncated(rp) && !isTRUE(quiet))
    warning(sprintf(paste0("annotate_decoupled_run: run_parameters.txt in %s is TRUNCATED ",
                           "(str() list.len cap). Run flags absent from the dump read as ",
                           "FALSE, so decoupled flags derived from them may be missing. ",
                           "Re-render the run, or set the flags by hand."), basename(dir)),
            call. = FALSE)
  tau_sw  <- isTRUE(.adr_param(rp, "osp_scale_is_tau"))
  dens    <- isTRUE(.adr_param(rp, "estimate_cpue_density"))
  use_f   <- isTRUE(.adr_param(rp, "use_crab_fraction"))
  osp_lo  <- isTRUE(.adr_param(rp, "use_osp_crab_lower"))
  shared  <- isTRUE(.adr_param(rp, "shared_tau"))
  op_mode <- .adr_param(rp, "opener_covariate_mode")

  rows <- list()
  for (f in sp) {
    label <- sub("^structural_params_", "", tools::file_path_sans_ext(basename(f)))
    d <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || !"parameter" %in% names(d)) next
    is_shore <- grepl("^shore", label)

    n_ie <- NA_integer_
    meth <- NA_character_
    if (!is.null(fds) && "fit" %in% names(fds)) {
      i <- fds$fit == label
      if (any(i)) {
        if ("n_ie_obs" %in% names(fds)) n_ie <- suppressWarnings(as.integer(fds$n_ie_obs[i][1]))
        if ("method"   %in% names(fds)) meth <- as.character(fds$method[i][1])
      }
    }

    base <- sub("\\[.*$", "", d$parameter)
    reason <- rep(NA_character_, nrow(d))
    unknown <- rep(FALSE, nrow(d))
    set <- function(mask, why) reason[mask & is.na(reason)] <<- why

    if (isTRUE(n_ie == 0)) set(base == "sigma_IE",
      "no I/E observations in this fit (n_ie_obs = 0); sigma_IE is its prior")
    else if (is.na(n_ie))  unknown[base == "sigma_IE"] <- TRUE

    if (is_shore) {
      set(base %in% c("kappa_OSP", "sigma_r_OSP", "r_OSP"),
          "shore fit: no OSP observations (OSP_n = 0)")
      set(base %in% c("R_G_boat", "R_T"),
          "shore fit: no boat count stream; the boat gear ratio is its prior")
      set(base %in% c("f_crab", "f_lower"),
          "shore fit: apply_crab_fraction = 0, f is pinned at 1")
    } else {
      if (tau_sw) set(base == "kappa_OSP",
        "osp_scale_is_tau = TRUE: the OSP mean uses L, not kappa_OSP; this is the prior")
      if (!use_f) set(base == "f_crab", "use_crab_fraction = FALSE: f is pinned")
      if (!osp_lo) set(base == "f_lower",
        "use_osp_crab_lower = FALSE: the OSP lower bound is off and f_lower is pinned at 0")
    }
    if (!dens) set(base == "gamma_C", "estimate_cpue_density = FALSE: the density term is inert")
    if (identical(op_mode, "off")) set(base == "B_open",
      "opener_covariate_mode = 'off': no opener covariate is active")
    if (!shared) set(base %in% c("tau_bar", "tau_bar_out"),
      "shared_tau = FALSE: L is per-day independent draws; there is no shared turnover")
    # Window-dependent rules cannot be rebuilt from a run folder.
    unknown[base %in% c("B1", "B2", "B1_C", "B2_C")] <- is.na(reason[base %in% c("B1","B2","B1_C","B2_C")])

    rows[[length(rows) + 1]] <- data.frame(
      fit = label, parameter = d$parameter,
      median = if ("median" %in% names(d)) d$median else NA_real_,
      mean   = if ("mean"   %in% names(d)) d$mean   else NA_real_,
      decoupled = !is.na(reason),
      decoupled_reason = reason,
      status = ifelse(!is.na(reason), "prior only",
               ifelse(unknown, "unknown (window-dependent; not reconstructible)", "estimate")),
      estimate = ifelse(is.na(reason) & !unknown,
                        if ("median" %in% names(d)) d$median else NA_real_, NA_real_),
      fit_method = meth,
      source = "reconstructed",
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(invisible(NULL))
  df <- do.call(rbind, rows)
  utils::write.csv(df, out_path, row.names = FALSE)
  if (!isTRUE(quiet))
    cat(sprintf("  %-52s %d rows, %d prior-only, %d unknown\n", basename(dir), nrow(df),
                sum(df$decoupled), sum(df$status == "unknown (window-dependent; not reconstructible)")))
  invisible(out_path)
}
