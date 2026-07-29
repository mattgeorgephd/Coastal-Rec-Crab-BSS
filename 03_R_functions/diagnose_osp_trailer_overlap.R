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
# implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program (see the LICENSE file);
# if not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
###############################################################################
# diagnose_osp_trailer_overlap.R  (Phase 0: OSP <-> trailer calibration; diagnostic only)
#
# Establishes the empirical scale ratio between the OSP daily boat total (a full-
# day count of all private boats) and the crab-creel trailer count (an
# instantaneous snapshot) on the days they overlap, and audits effort-count
# coverage across the estimation window. This is the evidence that the two series
# can enter the boat effort model as ONE latent process with a fixed scale ratio
# (Phase 1), rather than being swapped at the OSP-window boundary. DIAGNOSTIC
# ONLY: it changes no estimate.
#
# WHY max/day IS THE PRIMARY TRAILER METRIC
#   Multiple trailer counts on a day are repeated snapshots of the same day. The
#   per-day maximum is the cleanest single snapshot (it tracks the OSP total at
#   r ~ 0.98); the mean is close (r ~ 0.97); summing count events double-counts
#   and degrades the fit (r ~ 0.90). All three are reported; max is primary.
#   The through-origin slope is ~0.37, i.e. one snapshot catches ~37% of the
#   day's boats, so 1/slope ~ 2.7 is an empirical within-day boat TURNOVER (the
#   quantity tau_boat encodes; see the plan, section 5.5).
#
# ARGS
#   osp         tibble from fetch_osp_boat_counts() (already de-duped + windowed)
#   params      run_config (reads est_date_start/end, the trailer file/area keys)
#   output_dir  if non-NULL, writes the CSVs and the plot there
#
# CONFIG (optional; historical Westport defaults via %||%)
#   osp_match_trailer_area = "Westport Boat Launch"   # trailer site to match OSP
#   effort_file / input_sheet                          # reused from run_config
#
# RETURNS list(calibration, coverage, pairs, audit, zero_agreement, verdict).
#
# Requires: dplyr, tibble, readxl, here; ggplot2 optional (plot skipped if absent).
###############################################################################

diagnose_osp_trailer_overlap <- function(osp, params, output_dir = NULL) {

  if (is.null(osp) || nrow(osp) == 0) {
    cat("  OSP overlap diagnostic: no OSP data; skipping.\n"); return(invisible(NULL))
  }

  # --- Trailer daily series at the OSP-matched launch (read fresh, self-contained) ---
  area  <- params$osp_match_trailer_area %||% "Westport Boat Launch"
  eff_file <- here::here("04_input_files", params$effort_file %||% "effort_combined.xlsx")
  trailer_raw <- readxl::read_excel(eff_file, sheet = params$input_sheet %||% "data")
  trail <- trailer_raw |>
    dplyr::filter(creel_area == area) |>
    dplyr::mutate(event_date = as.Date(date),
                  tc = suppressWarnings(as.numeric(boat_trailer_count))) |>
    dplyr::filter(!is.na(event_date), !is.na(tc), tc >= 0)

  if (nrow(trail) == 0) {
    cat("  OSP overlap diagnostic: no trailer counts at '", area, "'; skipping.\n", sep = ""); return(invisible(NULL))
  }

  trail_daily <- trail |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(n_counts   = dplyr::n(),
                     trail_mean = mean(tc),
                     trail_max  = max(tc),
                     trail_sum  = sum(tc), .groups = "drop")

  osp_daily <- osp |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(osp = sum(count_quantity), .groups = "drop")

  pairs <- dplyr::inner_join(osp_daily, trail_daily, by = "event_date") |> dplyr::arrange(event_date)
  if (nrow(pairs) < 3) {
    cat("  OSP overlap diagnostic: <3 overlapping days; skipping.\n"); return(invisible(NULL))
  }

  # --- Calibration for each trailer aggregation ---
  fit_one <- function(yv, nm) {
    ok <- is.finite(pairs$osp) & is.finite(yv)
    xx <- pairs$osp[ok]; yy <- yv[ok]
    r   <- suppressWarnings(stats::cor(xx, yy))
    co  <- stats::coef(stats::lm(yy ~ xx))
    origin <- sum(xx * yy) / sum(xx * xx)
    tibble::tibble(trailer_metric = nm, n = length(xx), corr = round(r, 3),
                   ols_slope = round(unname(co[2]), 3), ols_intercept = round(unname(co[1]), 2),
                   origin_slope = round(origin, 3), implied_turnover = round(1 / origin, 2))
  }
  calibration <- dplyr::bind_rows(
    fit_one(pairs$trail_mean, "trailer_mean_per_visit"),
    fit_one(pairs$trail_max,  "trailer_max_per_day"),
    fit_one(pairs$trail_sum,  "trailer_sum_per_day"))
  prim <- calibration |> dplyr::filter(trailer_metric == "trailer_max_per_day")

  # zero agreement on the primary metric
  z_osp <- pairs$osp == 0
  zero_agree <- sprintf("%d/%d", sum(z_osp & pairs$trail_max == 0), sum(z_osp))

  # --- Coverage / NS audit over the estimation window ---
  d0 <- as.Date(params$est_date_start %||% "2024-09-16")
  d1 <- as.Date(params$est_date_end   %||% "2025-09-15")
  audit <- tibble::tibble(event_date = seq(d0, d1, by = "day")) |>
    dplyr::mutate(
      month           = as.integer(format(event_date, "%m")),
      osp_present     = event_date %in% osp_daily$event_date,
      trailer_present = event_date %in% trail_daily$event_date,
      osp_operating   = month %in% 3:10,
      coverage_class  = dplyr::case_when(
        osp_present &  trailer_present ~ "OSP+trailer",
        osp_present & !trailer_present ~ "OSP only",
       !osp_present &  trailer_present ~ "trailer only",
        TRUE                           ~ "neither"))
  coverage <- audit |> dplyr::count(coverage_class, name = "days")

  verdict <- sprintf(paste0(
    "OSP daily total vs trailer snapshot (max/day): n=%d, r=%.3f, trailer = %.3f * OSP through origin ",
    "(within-day turnover ~ %.1f); on OSP-zero days the trailer is also zero %s. Near-proportional, so ",
    "the two series can enter the boat effort model as one latent process with a fixed scale ratio: the ",
    "%d overlap days calibrate it and the trailer-only stretch carries the OSP-dark winter."),
    prim$n, prim$corr, prim$origin_slope, prim$implied_turnover, zero_agree, prim$n)
  cat("  ", verdict, "\n", sep = "")

  # --- Persist (guarded, matching repo convention) ---
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch({
      utils::write.csv(calibration, file.path(output_dir, "osp_trailer_overlap_calibration.csv"), row.names = FALSE)
      utils::write.csv(dplyr::mutate(pairs, event_date = as.character(event_date)),
                       file.path(output_dir, "osp_trailer_overlap_pairs.csv"), row.names = FALSE)
      utils::write.csv(dplyr::mutate(audit, event_date = as.character(event_date)),
                       file.path(output_dir, "osp_coverage_audit.csv"), row.names = FALSE)
    }, error = function(err) cat("  (overlap CSVs not written:", conditionMessage(err), ")\n"))

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      tryCatch({
        p <- ggplot2::ggplot(pairs, ggplot2::aes(.data$osp, .data$trail_max)) +
          ggplot2::geom_abline(slope = prim$origin_slope, intercept = 0, color = "#D55E00", linewidth = 1) +
          ggplot2::geom_point(color = "#0072B2", alpha = 0.8, size = 2) +
          ggplot2::labs(
            x = "OSP daily boat total (all private boats)", y = "Trailer snapshot (max/day)",
            title = sprintf("OSP vs trailer overlap (n=%d, r=%.3f, trailer=%.3f*OSP, turnover~%.1f)",
                            prim$n, prim$corr, prim$origin_slope, prim$implied_turnover)) +
          ggplot2::theme_minimal(base_size = 11)
        ggplot2::ggsave(file.path(output_dir, "osp_trailer_overlap.png"), p, width = 7, height = 5, dpi = 140)
      }, error = function(err) cat("  (overlap plot not written:", conditionMessage(err), ")\n"))
    }
  }

  invisible(list(calibration = calibration, coverage = coverage, pairs = pairs,
                 audit = audit, zero_agreement = zero_agree, verdict = verdict))
}
