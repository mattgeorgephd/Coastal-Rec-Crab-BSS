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
# desk_advisory_transfer_2026-11-19.R
#
# A DESK TEST, SECONDS, NO MCMC: does the boat advisory effect measured on ONE season predict
# the launch trailer counts of ANOTHER, on weeks neither it nor the target's own fit has seen?
# (CHANGE_REGISTER B80, D31 and D32; VALIDATION_CAMPAIGN Section 1z.21.)
#
# WHY. The boat SCA term (A30) was adopted on 2024-25. The 2025-26 fit gives it a coefficient
# of the same size (-1.32 against -1.15), but agreement of two in-sample coefficients is not a
# prediction, and the four-season desk screen of 2026-09-27 had already looked at the 2025-26
# trailer days. What nobody has run is the transfer: take the coefficient from the SOURCE
# season, do NOT re-estimate it, and ask whether it improves the prediction of held-out weeks
# of the TARGET season.
#
# WHAT. The daily sum of the launch trailer counts per sampled day (the series of
# desk_sca_season_split_2026-09-27.R), one season at a time as the target. Leave one ISO week
# out, fit the nuisance terms (calendar month, weekend, negative-binomial dispersion) on the
# target's other weeks, and score the held-out week's days under each model:
#   M0   no advisory term
#   TC   the CONSTANT coefficient of the source, as an offset (not re-estimated)
#   TS   the SPLIT coefficients of the source (winter, rest), as offsets
#   F    a constant coefficient estimated on the target's own training weeks (the ceiling a
#        carried coefficient is measured against)
# Two sources per target: the OTHER in-scope season alone, and every other season pooled (the
# out-of-scope seasons 2022-23 and 2023-24 lend their trailer days as the 2026-09-27 screen let
# them; the pooled source never contains the target).
#
# THE STATISTICS, per target and source, over the weeks every model could score:
#   TC - M0   does a carried constant coefficient predict the target?     (D31, out of sample)
#   TS - TC   does the carried split predict better than the carried constant?   (D32)
#   F  - TC   is the target's own coefficient detectably different from the carried one?
# each the sum over weeks of the held-out log predictive density difference, with a paired SE
# (sqrt(n) times the SD of the weekly differences), all weeks and winter / rest weeks apart.
#
# THE READING RULES are CHANGE_REGISTER B80's, written before this script was run on either
# in-scope season as a target. Its mechanics were dry-run on 2023-24 as the target, a season
# the rule does not judge.
#
# WRITES (merged by key) 05_output/advisory_replication_2026-11-19_transfer.csv (one row per
# target, source, comparison and week set) and _transfer_weeks.csv (the weekly scores).
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tibble); library(tidyr); library(lubridate) })
if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS is needed for the negative-binomial fits (install.packages('MASS')).")
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
source("run_config.R")
source("03_R_functions/batch_verdict_helpers.R")   # merge_csv_by

# TARGETS: the seasons scored. The default is the two in-scope seasons; --target=2023-24 (or
# ADVISORY_TRANSFER_TARGETS) dry-runs the mechanics on a season the rule does not judge.
.arg <- grep("^--target=", commandArgs(trailingOnly = TRUE), value = TRUE)
TARGETS <- if (length(.arg)) strsplit(sub("^--target=", "", .arg[1]), ",")[[1]] else {
  ev <- Sys.getenv("ADVISORY_TRANSFER_TARGETS", ""); if (nzchar(ev)) strsplit(ev, ",")[[1]] else c("2024-25", "2025-26") }
IN_SCOPE <- c("2024-25", "2025-26")

# ---- the series: the desk screen's, unchanged --------------------------------------------------
e <- read_excel("04_input_files/effort_combined.xlsx", sheet = "data") |> mutate(date = as.Date(date))
areas <- run_config$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch")
b <- e |> filter(creel_area %in% areas, !is.na(boat_trailer_count)) |>
  group_by(date, season) |> summarise(y = sum(boat_trailer_count), n_counts = n(), .groups = "drop")
rng <- range(b$date)
fl <- marine_hazard_flag_series(modifyList(run_config, list(est_date_start = rng[1], est_date_end = rng[2])), date_start = rng[1], date_end = rng[2])
wm <- run_config$marine_hazard_winter_months %||% c(12, 1, 2)
d <- b |> left_join(fl |> dplyr::select(event_date, nws_sca_any), by = c("date" = "event_date")) |> filter(!is.na(nws_sca_any)) |>
  mutate(month = factor(format(date, "%m")), wk = as.integer(wday(date) %in% c(1, 7)), season = as.character(season),
         winter = as.integer(as.integer(format(date, "%m")) %in% wm), sca = as.numeric(nws_sca_any),
         sca_w = sca * winter, sca_r = sca * (1 - winter), week = format(date, "%G-W%V"))
cat(sprintf("Launch trailer counts: %d sampled days, %s to %s; seasons %s\n", nrow(d), rng[1], rng[2], paste(sort(unique(d$season)), collapse = ", ")))

# ---- the source coefficients: a negative-binomial GLM on the source days, never on the target --
nb <- function(f, dd) suppressWarnings(MASS::glm.nb(f, data = dd))
src_coef <- function(dd, label) {
  multi <- length(unique(dd$season)) > 1
  fc <- if (multi) y ~ season + month + wk + sca else y ~ month + wk + sca
  fs <- if (multi) y ~ season + month + wk + sca_w + sca_r else y ~ month + wk + sca_w + sca_r
  mc <- nb(fc, dd); ms <- nb(fs, dd); cc <- summary(mc)$coefficients; cs <- summary(ms)$coefficients
  tibble(source = label, n_days = nrow(dd), n_advisory = sum(dd$sca), n_winter_advisory = sum(dd$sca_w),
         b_const = cc["sca", 1], se_const = cc["sca", 2], b_winter = cs["sca_w", 1], se_winter = cs["sca_w", 2],
         b_rest = cs["sca_r", 1], se_rest = cs["sca_r", 2])
}

# ---- leave one week out on the target ------------------------------------------------------------
lpd_nb <- function(m, new, off = 0) {
  mu <- exp(as.numeric(stats::predict(m, newdata = new, type = "link")))   # the offset is in the formula, so predict() carries it
  sum(stats::dnbinom(new$y, mu = mu, size = m$theta, log = TRUE))
}
score_target <- function(tg, src) {
  dt <- d[d$season == tg, , drop = FALSE]
  dt$off_c <- src$b_const * dt$sca; dt$off_s <- src$b_winter * dt$sca_w + src$b_rest * dt$sca_r
  out <- lapply(sort(unique(dt$week)), function(w) {
    te <- dt[dt$week == w, , drop = FALSE]; tr <- dt[dt$week != w, , drop = FALSE]
    if (!all(as.character(te$month) %in% as.character(tr$month)))
      return(tibble(week = w, n_days = nrow(te), scored = FALSE, why = "a held-out day's month has no training day"))
    tr$month <- droplevels(tr$month)
    fits <- tryCatch(list(M0 = nb(y ~ month + wk, tr), TC = nb(y ~ month + wk + offset(off_c), tr),
                          TS = nb(y ~ month + wk + offset(off_s), tr),
                          F  = if (sum(tr$sca) >= 3 && sum(1 - tr$sca) >= 3) nb(y ~ month + wk + sca, tr) else NULL),
                     error = function(e) NULL)
    if (is.null(fits) || is.null(fits$F)) return(tibble(week = w, n_days = nrow(te), scored = FALSE, why = "a training fit failed"))
    tibble(week = w, n_days = nrow(te), scored = TRUE, why = "", n_advisory = sum(te$sca), winter_week = as.integer(mean(te$winter) >= 0.5),
           y_sum = sum(te$y), M0 = lpd_nb(fits$M0, te), TC = lpd_nb(fits$TC, te), TS = lpd_nb(fits$TS, te), F = lpd_nb(fits$F, te))
  })
  bind_rows(out) |> mutate(target = tg, source = src$source, .before = 1)
}
paired <- function(x) { n <- length(x); if (n < 2) return(c(diff = sum(x), se = NA_real_, n = n, pos = sum(x > 0)))
  c(diff = sum(x), se = sqrt(n) * stats::sd(x), n = n, pos = sum(x > 0)) }
summarise_pairs <- function(wk) {
  s <- wk[wk$scored, , drop = FALSE]
  sets <- list(`all weeks` = rep(TRUE, nrow(s)), `winter weeks` = s$winter_week == 1, `rest weeks` = s$winter_week == 0)
  cmp <- list(`TC - M0` = s$TC - s$M0, `TS - TC` = s$TS - s$TC, `F - TC` = s$F - s$TC, `TS - M0` = s$TS - s$M0)
  bind_rows(lapply(names(cmp), function(cn) bind_rows(lapply(names(sets), function(sn) {
    p <- paired(cmp[[cn]][sets[[sn]]])
    tibble(comparison = cn, weeks = sn, n_weeks = unname(p["n"]), diff = unname(p["diff"]), se = unname(p["se"]),
           ratio = unname(p["diff"] / p["se"]), n_positive = unname(p["pos"]))
  }))))
}

res <- list(); wks <- list(); srcs <- list()
for (tg in TARGETS) {
  other <- setdiff(IN_SCOPE, tg)
  sources <- list()
  if (length(other) == 1 && tg %in% IN_SCOPE) sources[[other]] <- d[d$season == other, , drop = FALSE]
  if (!tg %in% IN_SCOPE) sources[["2024-25"]] <- d[d$season == "2024-25", , drop = FALSE]     # the dry run's single-season source
  sources[["all other seasons"]] <- d[d$season != tg, , drop = FALSE]
  for (sn in names(sources)) {
    sc <- src_coef(sources[[sn]], sn); wk <- score_target(tg, sc)
    srcs[[paste(tg, sn)]] <- sc |> mutate(target = tg, .before = 1)
    wks[[paste(tg, sn)]] <- wk
    res[[paste(tg, sn)]] <- summarise_pairs(wk) |> mutate(target = tg, source = sn, judged = tg %in% IN_SCOPE,   # 2023-24 is the dry run
                                                        n_weeks_target = length(unique(wk$week)), n_weeks_unscored = sum(!wk$scored), .before = 1)
  }
}
res <- bind_rows(res); wks <- bind_rows(wks); srcs <- bind_rows(srcs)

op <- options(width = 220); on.exit(options(op), add = TRUE)
cat("\n== THE SOURCE COEFFICIENTS (negative binomial; estimated without the target season)\n")
print(as.data.frame(srcs |> mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
cat("\n== HELD-OUT WEEKS OF THE TARGET: sum of the weekly log predictive density differences, paired SE, ratio\n")
print(as.data.frame(res |> mutate(across(c(diff, se, ratio), ~ round(.x, 2)))), row.names = FALSE)

merge_csv_by(as.data.frame(res |> left_join(srcs, by = c("target", "source"))),
             file.path("05_output", "advisory_replication_2026-11-19_transfer.csv"), key = c("target", "source", "comparison", "weeks"))
merge_csv_by(as.data.frame(wks), file.path("05_output", "advisory_replication_2026-11-19_transfer_weeks.csv"), key = c("target", "source", "week"))
cat("\nwrote 05_output/advisory_replication_2026-11-19_transfer.csv and _transfer_weeks.csv\n")
cat("READ BY THE RULES OF CHANGE_REGISTER B80. A screen on sampled days with a calendar-month adjustment, not the BSS.\n")
