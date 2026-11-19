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
# desk_advisory_brief_figures_2026-11-19.R
#
# THE FIGURES OF THE ADVISORY TERM BRIEF, SECONDS, NO MCMC.
# (07_documentation/development_notes/A30-advisory-term-brief-2026-11-19.docx and its .md
#  source; CHANGE_REGISTER B80, A30, D31 and D32; VALIDATION_CAMPAIGN Section 1z.21.)
#
# Draws PNGs into 07_documentation/development_notes/figures/A30/ from material that is
# committed: the launch trailer counts and the archived advisory flag, the tables of the
# transfer test and of the six arms (05_output/advisory_replication_2026-11-19_*.csv), the
# arms' kept evidence folders, and the B70 and B76 run folders. Fits nothing. A figure whose
# table is not there yet is skipped with a line saying so. With --docx it also builds the
# brief's Word file from its markdown (pandoc and figures/brief-reference.docx).
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tibble); library(tidyr); library(lubridate); library(ggplot2) })
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
source("run_config.R")

OUT <- file.path("07_documentation", "development_notes", "figures", "A30")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
tab <- function(kind) { f <- file.path("05_output", sprintf("advisory_replication_2026-11-19_%s.csv", kind))
  if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL }
C_ADV <- "#D55E00"; C_CLEAR <- "#0072B2"; C_GREY <- "#6E6E6E"; C_LIGHT <- "#B8B8B8"; C_BAND <- "#EFEFEF"; C_GREEN <- "#009E73"; C_PURPLE <- "#7A4EAB"
theme_brief <- function(base = 11.5) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = base + 1.5),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.subtitle = element_text(colour = "#444444"), plot.caption = element_text(colour = "#666666", hjust = 0, size = base - 3),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(10, 14, 8, 10))
save_fig <- function(p, name, w = 7, h = 4.3) { ggsave(file.path(OUT, name), p, width = w, height = h, dpi = 200, bg = "white"); cat("wrote", file.path(OUT, name), "\n") }
skip <- function(name, why) cat("SKIPPED", name, "(", why, ")\n")
IN_SCOPE <- c("2024-25", "2025-26")
wm <- as.integer(run_config$marine_hazard_winter_months %||% c(12, 1, 2))

# ---- the series the term is about: launch trailer counts per sampled day, and the flag ----------
e <- read_excel("04_input_files/effort_combined.xlsx", sheet = "data") |> mutate(date = as.Date(date))
areas <- run_config$boat_launch_areas %||% c("Westport Boat Launch", "Ocean Shores Boat Launch")
b <- e |> filter(creel_area %in% areas, !is.na(boat_trailer_count)) |> group_by(date, season) |> summarise(y = sum(boat_trailer_count), .groups = "drop")
rng <- range(b$date)
fl <- marine_hazard_flag_series(modifyList(run_config, list(est_date_start = rng[1], est_date_end = rng[2])), date_start = rng[1], date_end = rng[2])
d <- b |> left_join(fl |> dplyr::select(event_date, nws_sca_any), by = c("date" = "event_date")) |> filter(!is.na(nws_sca_any), season %in% IN_SCOPE) |>
  mutate(sca = nws_sca_any == 1, mon = as.integer(format(date, "%m")), winter = mon %in% wm,
         mlab = factor(format(date, "%b"), levels = c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")))

# ---- Figure 1: the raw contrast, month by month -----------------------------------------------
d1 <- d |> group_by(season, mlab, sca) |> summarise(n = n(), mean_y = mean(y), .groups = "drop") |>
  mutate(grp = ifelse(sca, "advisory in effect (04:00 to 16:00)", "no advisory"))
p1 <- ggplot(d1, aes(mlab, mean_y, colour = grp, group = grp)) +
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 0, ymax = Inf, fill = C_BAND) +
  geom_line(linewidth = 0.7) + geom_point(aes(size = n)) +
  facet_wrap(~season, ncol = 1) +
  scale_colour_manual(values = c("advisory in effect (04:00 to 16:00)" = C_ADV, "no advisory" = C_CLEAR)) +
  scale_size_area(max_size = 5, guide = "none") + scale_y_continuous(trans = "log1p", breaks = c(0, 2, 5, 10, 25, 50, 100)) +
  labs(title = "On an advisory day there are fewer trailers at the launch",
       subtitle = "Mean daily trailer count on sampled launch days, by month; point area is the number of days.\nShaded: December to February, where the two lines come closest.",
       x = NULL, y = "trailers per sampled day (log scale)",
       caption = "Source: effort_combined.xlsx (Westport and Ocean Shores launches) and nws_marine_hazards.xlsx,\nthrough marine_hazard_flag_series(), the production definition of an advisory day.") +
  theme_brief()
save_fig(p1, "fig1_trailers_by_month_and_advisory.png", 7, 5.4)

# ---- Figure 2: the coefficient, season by season ------------------------------------------------
fs_row <- function(dir, label, season, kind) {
  f <- file.path(dir, "bss_full_summary_private_boat_all_gear_Dungeness_Kept.csv"); if (!file.exists(f)) return(NULL)
  x <- utils::read.csv(f, row.names = 1, check.names = FALSE)["B_open_out[1]", ]
  tibble(label = label, season = season, kind = kind, est = x[["mean"]], lo = x[["2.5%"]], hi = x[["97.5%"]])
}
arms <- tab("arms")
d2 <- bind_rows(
  fs_row("05_output/20261019/pooled-CPUE-2024-25", "2024-25, pooled model (B70, the authoritative run)", "2024-25", "fitted model"),
  fs_row("05_output/20261019/gear-type-CPUE-model-2024-25", "2024-25, gear-resolved cross-check", "2024-25", "fitted model"),
  fs_row("05_output/20261020/pooled-CPUE-season-2025-26", "2025-26, pooled model (B76, readiness render)", "2025-26", "fitted model"),
  fs_row("05_output/20261020/gear-type-CPUE-model-season-2025-26", "2025-26, gear-resolved cross-check", "2025-26", "fitted model"))
if (!is.null(arms)) {
  a <- arms[arms$arm == "const" & is.finite(arms$b_const), , drop = FALSE]
  if (nrow(a)) d2 <- bind_rows(d2, tibble(label = sprintf("%s, pooled model re-rendered on this container (B80)", a$season), season = a$season, kind = "fitted model",
                                         est = a$b_const, lo = a$b_const_lo, hi = a$b_const_hi))
}
tr <- tab("transfer")
if (!is.null(tr)) {
  g <- tr |> filter(source %in% IN_SCOPE) |> distinct(source, b_const, se_const)
  d2 <- bind_rows(d2, tibble(label = sprintf("%s, desk regression on sampled days", g$source), season = g$source, kind = "desk regression",
                             est = g$b_const, lo = g$b_const - 1.96 * g$se_const, hi = g$b_const + 1.96 * g$se_const))
}
d2 <- d2 |> arrange(season, kind != "fitted model") |> mutate(label = factor(label, levels = rev(label)))
p2 <- ggplot(d2, aes(y = label, colour = season)) +
  geom_vline(xintercept = 1, colour = "#333333", linewidth = 0.4) +
  geom_errorbar(aes(xmin = exp(lo), xmax = exp(hi)), width = 0.25, linewidth = 0.6, orientation = "y") +
  geom_point(aes(x = exp(est), shape = kind), size = 2.8) +
  geom_text(aes(x = exp(hi) * 1.06, label = sprintf("%.2f [%.2f, %.2f]", exp(est), exp(lo), exp(hi))), hjust = 0, size = 3.1, colour = "#222222") +
  scale_x_log10(limits = c(0.13, 1.6), breaks = c(0.15, 0.2, 0.3, 0.5, 0.75, 1)) +
  scale_colour_manual(values = c("2024-25" = C_CLEAR, "2025-26" = C_ADV), guide = "none") + scale_shape_manual(values = c("fitted model" = 16, "desk regression" = 1)) +
  labs(title = "The advisory effect is the same size on the second season",
       subtitle = "Trailers on an advisory day as a share of an ordinary day's, with 95% intervals (log scale).",
       x = "rate ratio (1 = no effect)", y = NULL,
       caption = "Source: bss_full_summary_private_boat_all_gear_*.csv of the named runs (B_open_out[1]);\n05_output/advisory_replication_2026-11-19_arms.csv and _transfer.csv.") +
  theme_brief() + theme(panel.grid.major.y = element_blank())
save_fig(p2, "fig2_coefficient_by_season.png", 7.6, 0.9 + 0.42 * nrow(d2))

# ---- Figures 3 and 4: the transfer test -----------------------------------------------------------
if (is.null(tr)) { skip("fig3, fig4", "the transfer table is not written yet") } else {
  d3 <- tr |> filter(target %in% IN_SCOPE, comparison %in% c("TC - M0", "TS - TC")) |>
    mutate(what = ifelse(comparison == "TC - M0", "The carried coefficient against no term", "The carried split against the carried constant"),
           src = ifelse(source %in% IN_SCOPE, sprintf("from %s alone", source), "from every other season"),
           row = sprintf("target %s, %s", target, src), weeks = factor(weeks, levels = c("all weeks", "rest weeks", "winter weeks")))
  d3$row <- factor(d3$row, levels = rev(unique(d3$row[order(d3$target != "2025-26", d3$src != "from 2024-25 alone" & d3$src != "from 2025-26 alone")])))
  p3 <- ggplot(d3, aes(y = row, x = diff, colour = weeks)) +
    geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.4) +
    geom_errorbar(aes(xmin = diff - 2 * se, xmax = diff + 2 * se), width = 0.3, linewidth = 0.55, orientation = "y", position = position_dodge(width = 0.7)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.7)) +
    facet_wrap(~what, ncol = 1, scales = "free_x") +
    scale_colour_manual(values = c("all weeks" = "#222222", "rest weeks" = C_ADV, "winter weeks" = C_CLEAR), guide = guide_legend(reverse = TRUE)) +
    labs(title = "A coefficient measured on one season predicts the other's held-out weeks",
         subtitle = "Held-out log predictive density gained, summed over the target season's weeks.\nWhiskers are 2 paired standard errors, so a whisker clear of zero is a difference at more than 2 SE.",
         x = "gain in nats (positive: the first model predicts better)", y = NULL,
         caption = "Source: 05_output/advisory_replication_2026-11-19_transfer.csv. Desk regression on sampled launch days;\nthe coefficient is carried over from the source season, not re-estimated on the target.") +
    theme_brief() + theme(panel.grid.major.y = element_blank())
  save_fig(p3, "fig3_transfer_test.png", 7.4, 5.6)

  tw <- tab("transfer_weeks")
  iso_monday <- function(week) {                       # the Monday of an ISO week label, 2025-W49
    yr <- as.integer(substr(week, 1, 4)); wk <- as.integer(sub("^.*-W", "", week))
    jan4 <- as.Date(sprintf("%d-01-04", yr)); jan4 - (as.integer(format(jan4, "%u")) - 1L) + 7L * (wk - 1L) }
  d4 <- tw |> filter(target == "2025-26", source == "2024-25", scored) |>
    mutate(gain = TC - M0, start = iso_monday(week), grp = ifelse(winter_week == 1, "December to February", "rest of the season"))
  p4 <- ggplot(d4, aes(start, gain, fill = grp)) +
    geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.4) + geom_col(width = 6) +
    scale_fill_manual(values = c("December to February" = C_CLEAR, "rest of the season" = C_ADV)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(title = "Where the 2024-25 coefficient helps on 2025-26: mostly outside the winter",
         subtitle = sprintf("Held-out gain per week of 2025-26 from carrying over the 2024-25 coefficient.\nSum %+.1f nats over %d weeks, %d of them positive.",
                            sum(d4$gain), nrow(d4), sum(d4$gain > 0)),
         x = NULL, y = "gain in nats",
         caption = "Source: 05_output/advisory_replication_2026-11-19_transfer_weeks.csv (TC - M0, target 2025-26, source 2024-25).") +
    theme_brief()
  save_fig(p4, "fig4_transfer_by_week.png", 7, 3.9)
}

# ---- shared readers for the arms' evidence folders ------------------------------------------------
pairs <- tab("pairs"); monthly <- tab("monthly"); verdicts <- tab("verdicts"); bweeks <- tab("blockint_weeks"); bpairs <- tab("blockint")
EV <- file.path("05_output", "advisory_replication_2026-11-19_arms")
FITB <- "private_boat_all_gear_Dungeness_Kept"
rdf <- function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
iso_monday <- function(week) { yr <- as.integer(substr(week, 1, 4)); wk <- as.integer(sub("^.*-W", "", week))
  jan4 <- as.Date(sprintf("%d-01-04", yr)); jan4 - (as.integer(format(jan4, "%u")) - 1L) + 7L * (wk - 1L) }
top_stream <- function(stage) if (file.exists(file.path(EV, stage, sprintf("loo_block_joint_%s.csv", FITB)))) "joint" else "trailer"
SLAB <- c(S25 = "2024-25", S26 = "2025-26")

# ---- Figure 5: the held-out-week check stopped working with the per-day effect; the repair --------
S_OLD <- "before the per-day effect (B67 render, 2026-10-17)"
S_CND <- "current model, as the pipeline computes it"
S_INT <- "current model, day deviation integrated out"
k5 <- list()
for (sn in names(SLAB)) { st <- paste0(sn, "-const"); sm <- top_stream(st)
  cnd <- rdf(file.path(EV, st, sprintf("loo_block_%s_%s.csv", sm, FITB))); int <- rdf(file.path(EV, st, sprintf("loo_blockint_%s_%s.csv", sm, FITB)))
  if (is.null(cnd) || is.null(int)) next
  k5[[sn]] <- bind_rows(tibble(season = SLAB[[sn]], series = S_CND, week = cnd$block, k = cnd$pareto_k),
                        tibble(season = SLAB[[sn]], series = S_INT, week = int$block, k = int$pareto_k))
  if (sn == "S25") { b67 <- rdf(file.path("05_output", "20261017", "pooled-CPUE-2024-25", sprintf("loo_block_%s_%s.csv", sm, FITB)))
    if (!is.null(b67)) k5[[sn]] <- bind_rows(k5[[sn]], tibble(season = SLAB[[sn]], series = S_OLD, week = b67$block, k = b67$pareto_k)) }
}
k5 <- bind_rows(k5)
if (!nrow(k5)) { skip("fig5", "no const arm has both of its held-out-week tables yet") } else {
  n5 <- k5 |> group_by(season, series) |> summarise(ok = sum(is.finite(k) & k <= 0.7), n = n(), .groups = "drop")
  lab5 <- n5 |> group_by(season) |> summarise(lab = paste0(first(season), ": usable weeks ", paste(sprintf("%d of %d %s", ok, n,
      c("before the per-day effect", "as computed now", "with the repair")[match(series, c(S_OLD, S_CND, S_INT))])[order(match(series, c(S_OLD, S_CND, S_INT)))], collapse = "; ")), .groups = "drop")
  d5 <- k5 |> left_join(lab5, by = "season") |> mutate(start = iso_monday(week), kk = pmin(ifelse(is.finite(k), k, 1.6), 1.6),
                                                       series = factor(series, levels = c(S_OLD, S_CND, S_INT)))
  p5 <- ggplot(d5, aes(start, kk, colour = series, shape = series)) +
    annotate("rect", xmin = as.Date(-Inf), xmax = as.Date(Inf), ymin = 0.7, ymax = Inf, fill = C_BAND) +
    geom_hline(yintercept = 0.7, colour = "#333333", linewidth = 0.4, linetype = "22") +
    geom_point(size = 2, alpha = 0.9) +
    facet_wrap(~lab, ncol = 1, scales = "free_x") +
    scale_colour_manual(values = setNames(c(C_GREY, C_ADV, C_CLEAR), c(S_OLD, S_CND, S_INT))) +
    scale_shape_manual(values = setNames(c(1, 17, 16), c(S_OLD, S_CND, S_INT))) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") + scale_y_continuous(breaks = c(0, 0.35, 0.7, 1.05, 1.4)) +
    labs(title = "The held-out-week check, before and after the per-day effect",
         subtitle = "Each point is one held-out week of the boat all-gear fit. Pareto k says whether the\nweek's leave-out posterior can be approximated without a refit; above 0.7 (shaded)\nit cannot.",
         x = "week held out", y = "Pareto k",
         caption = "Source: loo_block_* and loo_blockint_* of the const arms (05_output/advisory_replication_2026-11-19_arms/)\nand loo_block_* of the B67 render (05_output/20261017/pooled-CPUE-2024-25). The 2024-25 const arm is the\nauthoritative configuration; the authoritative run's own table (B70) has the same count as computed now.") +
    theme_brief() + guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
  save_fig(p5, "fig5_held_out_week_check_reliability.png", 7.2, 6.2)
}

# ---- Figure 6: the fitted model on held-out weeks, with the term against without (R3a) ------------
d6 <- if (is.null(bweeks)) NULL else bweeks |> filter(pair == "const - off", grepl("integrated", table)) |>
  group_by(season) |> filter(stream == (if ("joint" %in% stream) "joint" else "trailer")) |> ungroup()
if (is.null(d6) || !nrow(d6)) { skip("fig6", "a const and an off arm of one season do not both have their integrated tables yet") } else {
  lab6 <- d6 |> group_by(season) |> summarise(n = n(), n_used = sum(used), s = sum(diff[used]), se = sqrt(sum(used)) * stats::sd(diff[used]), stream = first(stream), .groups = "drop") |>
    mutate(lab = sprintf("%s (%s): %+.1f nats over %d of %d weeks, paired SE %.1f, ratio %+.2f", season,
                         ifelse(stream == "joint", "trailer and OSP counts together", "trailer counts"), s, n_used, n, se, s / se))
  d6 <- d6 |> left_join(lab6 |> select(season, lab), by = "season") |>
    mutate(start = iso_monday(block), grp = ifelse(!used, "week not reliable (not counted)", ifelse(winter, "December to February", "rest of the season")))
  p6 <- ggplot(d6, aes(start, diff, fill = grp)) +
    geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.4) + geom_col(width = 6) +
    facet_wrap(~lab, ncol = 1, scales = "free") +
    scale_fill_manual(values = c("December to February" = C_CLEAR, "rest of the season" = C_ADV, "week not reliable (not counted)" = C_LIGHT)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(title = "The fitted model on held-out weeks: with the advisory term against without it",
         subtitle = "Each bar is one week held out of the boat all-gear fit: the gain in its predictive\ndensity from the term, with the day deviation integrated out (the repaired check of Figure 5).",
         x = "week held out", y = "gain in nats",
         caption = "Source: 05_output/advisory_replication_2026-11-19_blockint_weeks.csv (const minus off, the integrated table).") +
    theme_brief()
  save_fig(p6, "fig6_block_cv_by_week.png", 7.2, 5.6)
}

# ---- Figure 7: the winter question, coefficient by coefficient ------------------------------------
scr <- rdf("05_output/marine_hazard_2026-09-27_season_split_screen.csv")
sp <- if (is.null(arms)) NULL else arms[arms$arm == "split" & is.finite(arms$b_winter) & is.finite(arms$b_rest), , drop = FALSE]
if (is.null(sp) || !nrow(sp)) { skip("fig7", "no split arm is rendered yet") } else {
  long <- function(lab, grp, w, wsd, r, rsd) tibble(row = lab, group = grp, part = c("December to February", "rest of the season"),
                                                    est = c(w, r), lo = c(w - 1.96 * wsd, r - 1.96 * rsd), hi = c(w + 1.96 * wsd, r + 1.96 * rsd))
  G_FIT <- "Fitted model, split arm (B80)"; G_DESK <- "Desk regression on sampled days (the 2026-09-27 screen)"
  d7 <- bind_rows(lapply(seq_len(nrow(sp)), function(i) with(sp[i, ], tibble(row = season, group = G_FIT,
              part = c("December to February", "rest of the season"), est = c(b_winter, b_rest), lo = c(b_winter_lo, b_rest_lo), hi = c(b_winter_hi, b_rest_hi)))))
  if (nrow(sp) >= 2) { pw <- function(m, s) { w <- 1 / s^2; c(sum(w * m) / sum(w), sqrt(1 / sum(w))) }
    w <- pw(sp$b_winter, sp$b_winter_sd); r <- pw(sp$b_rest, sp$b_rest_sd)
    d7 <- bind_rows(d7, long("the two seasons combined", G_FIT, w[1], w[2], r[1], r[2])) }
  if (!is.null(scr)) { s <- scr[scr$model == "split" & scr$family == "negbin" & scr$seasons %in% c(IN_SCOPE, "all"), ]
    for (ss in unique(s$seasons)) { x <- s[s$seasons == ss, ]; w <- x[x$term == "nws_sca_any_winter", ]; r <- x[x$term == "nws_sca_any_rest", ]
      if (nrow(w) == 1 && nrow(r) == 1) d7 <- bind_rows(d7, long(if (ss == "all") "all four seasons" else ss, G_DESK, w$estimate, w$se, r$estimate, r$se)) } }
  d7 <- d7 |> mutate(row = factor(row, levels = rev(c("2024-25", "2025-26", "the two seasons combined", "all four seasons"))),
                     group = factor(group, levels = unique(group)))
  p7 <- ggplot(d7, aes(y = row, colour = part)) +
    geom_vline(xintercept = 1, colour = "#333333", linewidth = 0.4) +
    geom_errorbar(aes(xmin = exp(lo), xmax = exp(hi)), width = 0.28, linewidth = 0.6, orientation = "y", position = position_dodge(width = 0.6)) +
    geom_point(aes(x = exp(est)), size = 2.6, position = position_dodge(width = 0.6)) +
    facet_wrap(~group, ncol = 1, scales = "free_y") +
    scale_x_log10(breaks = c(0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5)) +
    scale_colour_manual(values = c("December to February" = C_CLEAR, "rest of the season" = C_ADV), guide = guide_legend(reverse = TRUE)) +
    labs(title = "The advisory effect in the winter months and in the rest of the season",
         subtitle = "Trailers on an advisory day as a share of an ordinary day's, estimated separately for December to\nFebruary and for the other months. 95% intervals, log scale.",
         x = "rate ratio (1 = no effect)", y = NULL,
         caption = "Source: 05_output/advisory_replication_2026-11-19_arms.csv (b_winter, b_rest) and marine_hazard_2026-09-27_season_split_screen.csv.\nThe combined row is the precision-weighted mean of the two seasons' coefficients.") +
    theme_brief() + theme(panel.grid.major.y = element_blank())
  save_fig(p7, "fig7_winter_against_rest.png", 7.2, 4.9)
}

# ---- Figure 8: what the shape does to the boat estimate, month by month ---------------------------
if (is.null(monthly) || !nrow(monthly)) { skip("fig8", "no arm is rendered yet") } else {
  m <- monthly |> mutate(mon = as.integer(substr(month, 6, 7)), mlab = factor(month.abb[mon], levels = c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")),
                         armlab = c(const = "constant term (shipped)", split = "winter and rest split", off = "no term")[arm]) |>
    filter(!is.na(mlab)) |> group_by(season, arm, armlab, mlab) |> summarise(v = sum(boat_catch_median), .groups = "drop")
  p8 <- ggplot(m, aes(mlab, v, colour = armlab, group = armlab)) +
    annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 0, ymax = Inf, fill = C_BAND) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
    facet_wrap(~season, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c("constant term (shipped)" = "#222222", "winter and rest split" = C_ADV, "no term" = C_CLEAR)) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Private-boat catch by month under each form of the term",
         subtitle = "Posterior medians. Shaded: December to February.\n2025-26 is a readiness render and none of its figures is an estimate.",
         x = NULL, y = "private-boat catch, crab",
         caption = "Source: 05_output/advisory_replication_2026-11-19_monthly.csv (monthly_estimates_by_mode.csv of each arm).") +
    theme_brief()
  save_fig(p8, "fig8_monthly_boat_catch_by_arm.png", 7, 5.2)
}

# ---- Figure 9: each arm against the shipped one ---------------------------------------------------
if (is.null(arms) || !any(arms$arm == "const")) { skip("fig9", "no const arm is rendered yet") } else {
  a9 <- arms |> group_by(season) |> mutate(ref_tot = boat_all_gear_catch[arm == "const"][1], ref_win = boat_winter_months_catch[arm == "const"][1],
                                           ref_port = port_total[arm == "const"][1]) |> ungroup() |> filter(arm != "const") |>
    transmute(season, armlab = c(split = "winter and rest split", off = "no term")[arm],
              `boat catch, December to February` = 100 * (boat_winter_months_catch / ref_win - 1),
              `boat all-gear catch, whole season` = 100 * (boat_all_gear_catch / ref_tot - 1),
              `port total` = 100 * (port_total / ref_port - 1)) |>
    pivot_longer(-c(season, armlab), names_to = "what", values_to = "pct") |>
    mutate(what = factor(what, levels = rev(c("boat catch, December to February", "boat all-gear catch, whole season", "port total"))))
  p9 <- ggplot(a9, aes(pct, what, fill = armlab)) +
    geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.4) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%+.1f%%", pct), hjust = ifelse(pct >= 0, -0.12, 1.12)), position = position_dodge(width = 0.75), size = 3.1, colour = "#222222") +
    facet_wrap(~season, ncol = 1) +
    scale_fill_manual(values = c("winter and rest split" = C_ADV, "no term" = C_CLEAR)) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    labs(title = "What changing the term would do, against the shipped constant term",
         subtitle = "Percent change in the posterior median, by season.",
         x = "change against the constant term, percent", y = NULL,
         caption = "Source: 05_output/advisory_replication_2026-11-19_arms.csv. 2025-26 is a readiness render (D18): its changes are shown, its levels are not estimates.") +
    theme_brief() + theme(panel.grid.major.y = element_blank())
  save_fig(p9, "fig9_arms_against_shipped.png", 7, 4.8)
}

# ---- The Word file ---------------------------------------------------------------------------------
# --docx builds the brief's .docx from its markdown with pandoc and the shared reference template. The
# markdown is the source; the Word file is never edited by hand.
if ("--docx" %in% commandArgs(trailingOnly = TRUE)) {
  dn <- file.path("07_documentation", "development_notes"); md <- "A30-advisory-term-brief-2026-11-19.md"
  if (!nzchar(Sys.which("pandoc"))) stop("pandoc is not on the PATH")
  if (!file.exists(file.path(dn, md))) stop("the brief's markdown is not there: ", file.path(dn, md))
  old <- setwd(dn)
  rc <- system2("pandoc", c(md, "-f", "gfm", "-t", "docx", "--reference-doc=figures/brief-reference.docx", "--shift-heading-level-by=-1",
                            "-o", sub("\\.md$", ".docx", md)))
  setwd(old)
  if (!identical(rc, 0L)) stop("pandoc failed with status ", rc)
  cat("wrote", file.path(dn, sub("\\.md$", ".docx", md)), "\n")
}
cat("done\n")
