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
# gear_share_bootstrap.R  (review item 6, 2026-09-08)
#
# Catch shares by gear type with TRIP-LEVEL bootstrap uncertainty, replacing the
# Dirichlet(catch counts + 0.5) that the pooled report used until 2026-09-08.
#
# WHY THE DIRICHLET WAS TOO PRECISE. Dirichlet(n_1 + 0.5, ..., n_G + 0.5) with n_g the
# total CRAB caught on gear g is the posterior for a multinomial in which every crab is
# an independent draw. Crab are not independent: they arrive in trips (a pot boat lands
# 25 crab on one interview), so the effective sample size is closer to the number of
# TRIPS than to the number of crab, and the shares' intervals were several times too
# narrow. The nonparametric bootstrap resamples trips with replacement and recomputes
# the catch-weighted share vector each time, which respects the clustering exactly and
# needs no distributional assumption.
#
# PER SUB-SEASON. The old apportionment pooled every interview of a population across
# the season and applied one share vector to the population's whole-season total, which
# attributed "Pot" catch to the pot-closure period, where pots are illegal. Shares are
# now computed per population x sub-season and applied to that sub-season's own total
# draws (BSS draws where the gate passed, the PE point otherwise).
#
# gear_primary_class()  the same regex classification the report always used
# gear_share_bootstrap() n_boot x G matrix of share draws for one interview set
###############################################################################

gear_primary_levels <- c("Pot", "Ring Net", "Trap", "Snare", "Other")

gear_primary_class <- function(gear_type) {
  gt <- as.character(gear_type)
  dplyr::case_when(
    stringr::str_detect(gt, "(?i)^pot|, pot|pot,") & !stringr::str_detect(gt, "(?i)slip ring") ~ "Pot",
    stringr::str_detect(gt, "(?i)ring net")     ~ "Ring Net",
    stringr::str_detect(gt, "(?i)trap|star")    ~ "Trap",
    stringr::str_detect(gt, "(?i)snare")        ~ "Snare",
    TRUE                                        ~ "Other"
  )
}

# interviews: a frame with gear_primary (character) and catch (numeric >= 0), one row per
#             trip. Zero-catch trips are legitimate rows: resampling them changes how many
#             catch-bearing trips a replicate holds, which is part of the sampling variance.
# Returns an n_boot x length(levels) matrix of shares (rows sum to 1). If no trip carries
# catch the shares are undefined: every row is NA and the caller decides the fallback.
gear_share_bootstrap <- function(interviews, n_boot = 4000, levels = gear_primary_levels,
                                 seed = NULL) {
  stopifnot(is.data.frame(interviews), all(c("gear_primary", "catch") %in% names(interviews)))
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(interviews)
  out <- matrix(NA_real_, n_boot, length(levels), dimnames = list(NULL, levels))
  ct <- as.numeric(interviews$catch); ct[!is.finite(ct) | ct < 0] <- 0
  g  <- match(as.character(interviews$gear_primary), levels)
  g[is.na(g)] <- match("Other", levels)
  if (n == 0 || sum(ct) <= 0) return(out)
  # Each trip contributes its whole catch to exactly one gear class, so a replicate's
  # share vector is the gear-wise sum of the resampled trips' catch, normalized.
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    tot <- vapply(seq_along(levels), function(k) sum(ct[idx][g[idx] == k]), numeric(1))
    s <- sum(tot)
    out[b, ] <- if (s > 0) tot / s else NA_real_
  }
  out
}

# Point shares (catch-weighted, no resampling) for the same interview set, used for the
# table's median column cross-check and as the fallback when a replicate has no catch.
gear_share_point <- function(interviews, levels = gear_primary_levels) {
  ct <- as.numeric(interviews$catch); ct[!is.finite(ct) | ct < 0] <- 0
  g  <- match(as.character(interviews$gear_primary), levels); g[is.na(g)] <- match("Other", levels)
  tot <- vapply(seq_along(levels), function(k) sum(ct[g == k]), numeric(1))
  if (sum(tot) <= 0) return(setNames(rep(NA_real_, length(levels)), levels))
  setNames(tot / sum(tot), levels)
}
