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
# osp_sampling_rates.R  (2026-09-28, B48; read by fetch_osp_boat_counts())
#
# HOW MANY BOATS OSP CLASSIFIED ON A DAY: the binomial n of the crabbing-only share.
#
# OSP does not interview every boat. Per its sampling manual (Erica, OSP, 2026-09-28), a
# sampler picks a landmark every returning boat must pass and samples every k-th boat all
# day at a rate fixed in the morning from the ANTICIPATED effort and the staff on the docks;
# private boats and charters carry independent rates. The manual's schedule gives the
# MINIMUM rate by the day's exit/entrance count (04_input_files/osp_sampling_rates.xlsx):
#   < 30 boats 100%, 30-50 80%, 51-75 2/3, 76-100 50%, 101-150 40%, 151-200 1/3, > 200 25%.
# So the "crabbing_only" count OSP delivers is a count among the SAMPLED boats, and the
# share it measures is crabbing_only / n_sampled with n_sampled = rate x total, not
# crabbing_only / total. Dividing by the total would read a 50%-rate day's share at half.
#
# WHERE n COMES FROM, per day, first available wins:
#   1. osp_crab_checked_col   the number of private boats OSP sampled that day (exact)
#   2. osp_sample_rate_col    the rate OSP used that day, times the day's private total
#   3. the schedule           the manual's rate for the day's count: a MINIMUM, because a
#                             sampler may run above it ("these may be adjusted"), so n is
#                             then a lower bound and the share an UPPER bound on the truth;
#                             every such day is labelled "schedule (minimum rate)"
#   osp_sampling_rate_source = "none" restores the pre-B48 reading (n = the day's total),
#   for a delivery known to cover every boat.
# The rate column may be a fraction (0.5) or a percent (50); a column with any value above 1
# is read as percent, and one outside (0, 100] stops.
###############################################################################

# The schedule: one row per count band, rate = numerator / denominator (2/3 is exact, not
# 0.667). count_max NA means no upper limit.
read_osp_sampling_rates <- function(params = list()) {
  f <- here::here("04_input_files", params$osp_sampling_rates_file %||% "osp_sampling_rates.xlsx")
  if (!file.exists(f))
    stop("read_osp_sampling_rates(): the OSP sampling-rate schedule is missing (", f, ").", call. = FALSE)
  tab <- read_input_workbook(f, sheet = params$osp_sampling_rates_sheet %||% "data")
  need <- c("count_min", "count_max", "rate_numerator", "rate_denominator")
  if (!all(need %in% names(tab)))
    stop("read_osp_sampling_rates(): ", basename(f), " needs columns ", paste(need, collapse = ", "), call. = FALSE)
  tab <- data.frame(count_min = as.numeric(tab$count_min), count_max = suppressWarnings(as.numeric(tab$count_max)),
                    rate = as.numeric(tab$rate_numerator) / as.numeric(tab$rate_denominator))
  tab <- tab[order(tab$count_min), , drop = FALSE]
  if (any(!is.finite(tab$rate)) || any(tab$rate <= 0 | tab$rate > 1))
    stop("read_osp_sampling_rates(): every rate must lie in (0, 1].", call. = FALSE)
  hi <- ifelse(is.na(tab$count_max), Inf, tab$count_max)
  if (tab$count_min[1] > 0 || any(tab$count_min[-1] != hi[-nrow(tab)] + 1) || is.finite(hi[nrow(tab)]))
    stop("read_osp_sampling_rates(): the count bands must start at 0, be contiguous (each band starts one ",
         "above the last one's end) and leave the last band open (count_max blank).", call. = FALSE)
  tab
}

# The schedule's minimum rate for each day's count.
osp_rate_for_count <- function(count, tab) {
  k <- findInterval(as.numeric(count), tab$count_min)
  out <- rep(NA_real_, length(count))
  ok <- is.finite(count) & k >= 1
  out[ok] <- tab$rate[k[ok]]
  out
}

# A rate column as fractions: fractions kept, percents (any value above 1) divided by 100.
osp_rate_as_fraction <- function(x, col = "the rate column") {
  v <- suppressWarnings(as.numeric(x))
  f <- v[is.finite(v)]
  if (!length(f)) return(v)
  if (any(f <= 0) || any(f > 100))
    stop("OSP sampling rate: ", col, " holds a value outside (0, 100]; a rate is a fraction (0.5) or a percent (50).",
         call. = FALSE)
  if (any(f > 1)) {
    message("  OSP sampling rate: ", col, " read as PERCENT (a value above 1).")
    v <- v / 100
  }
  v
}

# Per day: n sampled, the rate behind it, and where it came from. `total` is the day's private
# boat count; `classified` and `rate` are the optional per-day columns (NA where absent).
osp_resolve_sample_n <- function(total, classified = NA_real_, rate = NA_real_, params = list(), tab = NULL) {
  src <- params$osp_sampling_rate_source %||% "auto"
  if (!src %in% c("auto", "column", "schedule", "none"))
    stop("osp_sampling_rate_source must be auto | column | schedule | none (got '", src, "').", call. = FALSE)
  m <- length(total)
  classified <- rep_len(as.numeric(classified), m); rate <- rep_len(as.numeric(rate), m)
  n <- rep(NA_real_, m); r <- rep(NA_real_, m); how <- rep(NA_character_, m)
  if (identical(src, "none"))
    return(data.frame(n = total, rate = ifelse(total > 0, 1, NA_real_), source = "none (every boat)"))
  if (identical(src, "schedule")) { classified[] <- NA_real_; rate[] <- NA_real_ }   # the schedule even when columns exist
  has_n <- is.finite(classified) & classified > 0
  n[has_n] <- classified[has_n]; r[has_n] <- pmin(classified[has_n] / pmax(total[has_n], 1), 1)
  how[has_n] <- "sampled count"
  has_r <- !has_n & is.finite(rate) & rate > 0
  n[has_r] <- round(total[has_r] * rate[has_r]); r[has_r] <- rate[has_r]; how[has_r] <- "rate column"
  rest <- is.na(how)
  if (any(rest) && identical(src, "column"))
    stop(sprintf("osp_sampling_rate_source = \"column\": %d OSP day(s) carry neither a sampled count nor a rate.", sum(rest)),
         call. = FALSE)
  if (any(rest)) {
    if (is.null(tab)) tab <- read_osp_sampling_rates(params)
    r[rest] <- osp_rate_for_count(total[rest], tab)
    n[rest] <- round(total[rest] * r[rest]); how[rest] <- "schedule (minimum rate)"
  }
  # a day with boats always has at least one sampled boat under every-k-th sampling
  n <- ifelse(total > 0, pmax(n, 1), 0)
  data.frame(n = n, rate = r, source = how)
}
