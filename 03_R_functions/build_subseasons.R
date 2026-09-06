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
# build_subseasons.R  (shared by the pooled and gear-resolved drivers)
#
# Build the list of within-season sub-seasons, split by the pot-closure window.
#
# The "pot closure" is the period when pots are not legal, so only non-pot gear
# (ring nets, snares, traps) may be used. It is NOT ring-net-only, so the internal
# label "ring_net_only" is kept only for output-filename continuity; the display_name
# is "Pot closure" and user-facing strings should use it. Outside the closure window,
# pots are allowed (all-gear).
#
# The closure window is given explicitly by params$pot_closure_start /
# params$pot_closure_end (added 2026-07-13). It is no longer assumed to start at the
# season start. This builder therefore handles the general case:
#   - an optional pre-closure all-gear period  [est_date_start, pot_closure_start - 1]
#   - the pot closure                          [pot_closure_start, pot_closure_end]
#   - an optional post-closure all-gear period [pot_closure_end + 1, est_date_end]
# A zero-length period is dropped.
#
# Backward compatibility: when pot_closure_start == est_date_start (the historical
# assumption, and the 2024-25 config), there is no pre-closure period and the result
# is the historical two-element list with the SAME internal names ("ring_net_only",
# "all_gear"), so keys, Stan fit labels, and output filenames are unchanged. Only
# when a mid-season closure is configured do the two all-gear periods take distinct
# names ("all_gear_pre" / "all_gear_post") to keep the population keys unique.
#
# Defaults: if params$pot_closure_start / _end are absent, they fall back to the old
# derivation (season start, and the day before params$pot_open_date), so an old config
# still runs. Each element also carries gear-model fields (period_bss, gear_exclude)
# that are inert for the pooled model.
#
# The model logic keys off gear_regime ("pot_closure" | "all_gear"), never off name.
###############################################################################

build_subseasons <- function(params) {
  est_start <- as.Date(params$est_date_start)
  est_end   <- as.Date(params$est_date_end)

  # -------------------------------------------------------------------------
  # 2026-09-10: MULTIPLE CLOSURE WINDOWS (multi-season spans). params$pot_closures, when
  # non-NULL, is a list of list(season =, start =, end =), one entry per pot closure in
  # the span, and OUTRANKS the scalar pot_closure_start/end keys. Each closure yields a
  # pot-closure sub-season named ring_net_only_<season> and the all-gear gap FOLLOWING it
  # is all_gear_<season> (the Dec-Sep block belongs to the season whose closure opened
  # it); a gap BEFORE the first closure is all_gear_pre_<season1>. Every sub-season
  # carries a `season` field so the report can roll fits up to season totals. Closures
  # outside the window are dropped; overlapping or unordered closures stop loudly.
  # With NULL (the default) or a single in-window closure, the SCALAR path below runs
  # unchanged, byte-for-byte, so every historical config keeps its exact sub-season
  # names, fit labels and output filenames.
  # -------------------------------------------------------------------------
  pcs <- params$pot_closures
  if (!is.null(pcs) && length(pcs)) {
    cl <- lapply(pcs, function(x) list(
      season = as.character(x$season %||% NA_character_),
      start  = as.Date(x$start), end = as.Date(x$end)))
    if (any(vapply(cl, function(x) is.na(x$season) || !nzchar(x$season), logical(1))))
      stop("build_subseasons(): every pot_closures entry needs a season label.", call. = FALSE)
    if (any(vapply(cl, function(x) x$end < x$start, logical(1))))
      stop("build_subseasons(): a pot_closures entry has end before start.", call. = FALSE)
    cl <- cl[order(vapply(cl, function(x) as.numeric(x$start), numeric(1)))]
    # drop closures that do not intersect the window; clamp the rest into it
    cl <- Filter(function(x) !(x$end < est_start || x$start > est_end), cl)
    cl <- lapply(cl, function(x) { x$start <- max(x$start, est_start)
                                   x$end   <- min(x$end,   est_end); x })
    if (length(cl) >= 2)
      for (k in 2:length(cl))
        if (cl[[k]]$start <= cl[[k - 1]]$end)
          stop("build_subseasons(): pot_closures windows overlap.", call. = FALSE)
    san <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)
    mk_ag <- function(nm, s, e, disp, season) list(
      name = nm, display_name = disp, gear_regime = "all_gear", season = season,
      start = s, end = e, period_bss = "month", gear_exclude = character(0))
    mk_cl <- function(x) list(
      name = paste0("ring_net_only_", san(x$season)), display_name = paste0("Pot closure ", x$season),
      gear_regime = "pot_closure", season = x$season,
      start = x$start, end = x$end, period_bss = "biweekly", gear_exclude = c("Pot"))
    if (length(cl) >= 2) {
      ss <- list()
      if (cl[[1]]$start > est_start)
        ss <- c(ss, list(mk_ag(paste0("all_gear_pre_", san(cl[[1]]$season)),
                               est_start, cl[[1]]$start - 1,
                               paste0("All gear (pre-closure ", cl[[1]]$season, ")"), cl[[1]]$season)))
      for (k in seq_along(cl)) {
        ss <- c(ss, list(mk_cl(cl[[k]])))
        gap_start <- cl[[k]]$end + 1
        gap_end   <- if (k < length(cl)) cl[[k + 1]]$start - 1 else est_end
        if (gap_end >= gap_start)
          ss <- c(ss, list(mk_ag(paste0("all_gear_", san(cl[[k]]$season)), gap_start, gap_end,
                                 paste0("All gear ", cl[[k]]$season), cl[[k]]$season)))
      }
      return(ss)
    }
    # 0 or 1 in-window closure: fall through to the scalar path so names stay legacy,
    # feeding it the one closure (or a non-intersecting sentinel for zero).
    if (length(cl) == 1) {
      params$pot_closure_start <- as.character(cl[[1]]$start)
      params$pot_closure_end   <- as.character(cl[[1]]$end)
    } else {
      params$pot_closure_start <- as.character(est_end + 365)
      params$pot_closure_end   <- as.character(est_end + 366)
    }
  }

  pc_start  <- as.Date(params$pot_closure_start %||% params$est_date_start)
  pc_end    <- as.Date(params$pot_closure_end   %||% (as.Date(params$pot_open_date) - 1))

  # A genuinely inverted CONFIG is an error worth stopping on.
  if (pc_end < pc_start)
    stop("build_subseasons(): pot_closure_end is before pot_closure_start.", call. = FALSE)

  # 2026-09-09: an estimation window that does not intersect the closure at all is a
  # single all-gear sub-season, not an error. Before this fix, clamping the closure into
  # the window inverted it and hit the stop() above, so a partial-season run that excluded
  # the closure (a summer-only window; a winter window starting after pots opened) could
  # not run at all. That contradicts the design goal that the model runs on ANY window the
  # user selects, part-season included; 2024-25 full-season configs are unaffected because
  # their closure always intersects the window.
  if (pc_end < est_start || pc_start > est_end) {
    return(list(list(
      name = "all_gear", display_name = "All gear", gear_regime = "all_gear",
      season = as.character(params$season_filter %||% NA_character_)[1],
      start = est_start, end = est_end, period_bss = "month",
      gear_exclude = character(0))))
  }

  if (pc_start < est_start) pc_start <- est_start
  if (pc_end   > est_end)   pc_end   <- est_end

  has_pre  <- pc_start > est_start
  has_post <- pc_end   < est_end

  season_tag <- as.character(params$season_filter %||% NA_character_)[1]
  allgear <- function(nm, s, e, disp) list(
    name = nm, display_name = disp, gear_regime = "all_gear", season = season_tag,
    start = s, end = e, period_bss = "month", gear_exclude = character(0))
  closure <- list(
    name = "ring_net_only", display_name = "Pot closure", gear_regime = "pot_closure",
    season = season_tag,
    start = pc_start, end = pc_end, period_bss = "biweekly", gear_exclude = c("Pot"))

  ss <- list()
  if (has_pre)
    ss <- c(ss, list(allgear(if (has_post) "all_gear_pre" else "all_gear",
                             est_start, pc_start - 1,
                             if (has_post) "All gear (pre-closure)" else "All gear")))
  ss <- c(ss, list(closure))
  if (has_post)
    ss <- c(ss, list(allgear(if (has_pre) "all_gear_post" else "all_gear",
                             pc_end + 1, est_end,
                             if (has_pre) "All gear (post-closure)" else "All gear")))
  ss
}
