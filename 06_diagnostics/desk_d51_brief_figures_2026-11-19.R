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
# desk_d51_brief_figures_2026-11-19.R
#
# THE FIGURES OF THE D51 DECISION BRIEF, SECONDS, NO MCMC.
# (07_documentation/development_notes/D51-shore-ie-likelihood-brief-2026-11-19.md;
#  CHANGE_REGISTER D51 and B78; VALIDATION_CAMPAIGN Section 1z.20.)
#
# Draws eight PNGs into 07_documentation/development_notes/figures/D51/ from material that is
# already committed: the four tables of the D50 shore turnover desk comparison
# (05_output/shore_turnover_d50_2026-11-18_*.csv, written by
# desk_shore_turnover_d50_2026-11-18.R), the shore gear counts through fetch_crab_data(), and
# the WDF20 tally in ingress_egress.xlsx. It computes nothing the tables do not already hold,
# except the count-hour shares and the observed two-float pairs it plots. Fits nothing.
#
# Re-run it after re-running the desk comparison; the figures carry no timestamp, so an
# unchanged input rewrites identical files.
###############################################################################
.root <- getwd()
if (!dir.exists(file.path(.root, "03_R_functions")) && dir.exists(file.path(.root, "..", "03_R_functions")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "03_R_functions"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
suppressPackageStartupMessages({ library(here); library(readxl); library(dplyr); library(tibble); library(tidyr)
                                 library(purrr); library(stringr); library(lubridate); library(ggplot2) })
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
invisible(lapply(list.files("03_R_functions", full.names = TRUE), function(f) source(f)))
source("run_config.R")

OUT <- file.path("07_documentation", "development_notes", "figures", "D51")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
tab <- function(kind) utils::read.csv(file.path("05_output", sprintf("shore_turnover_d50_2026-11-18_%s.csv", kind)), stringsAsFactors = FALSE)
by_day <- tab("by_day"); factors <- tab("factors"); by_count <- tab("by_count"); by_season <- tab("by_season")
quietly <- function(expr) { r <- NULL; utils::capture.output(r <- expr); r }

# one palette for the whole brief: the derivation blue, the fit vermillion, everything else grey
C_DER <- "#0072B2"; C_FIT <- "#D55E00"; C_GREY <- "#6E6E6E"; C_LIGHT <- "#B8B8B8"; C_BAND <- "#EFEFEF"; C_GREEN <- "#009E73"
theme_brief <- function(base = 12) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = base + 2),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.subtitle = element_text(colour = "#444444"), plot.caption = element_text(colour = "#666666", hjust = 0, size = base - 3),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(12, 16, 10, 12))
save_fig <- function(p, name, w = 8, h = 4.8) {
  ggsave(file.path(OUT, name), p, width = w, height = h, dpi = 200, bg = "white")
  cat("wrote", file.path(OUT, name), "\n")
}
hm <- function(h) sprintf("%02d:%02d", floor(h), round(60 * (h - floor(h))))

d26 <- by_day |> filter(season_window == "2025-26", fit == "shore_all_gear") |> mutate(event_date = as.Date(event_date))
f26 <- factors |> filter(season_window == "2025-26", fit == "shore_all_gear")
f25 <- factors |> filter(season_window == "2024-25", fit == "shore_all_gear")

# ---- Figure 1: what each reading compares (a schematic) ------------------------------------
box <- function(xmin, xmax, ymin, ymax, fill, colour) annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, colour = colour, linewidth = 0.7)
arr <- function(x, y, xend, yend, colour) annotate("segment", x = x, y = y, xend = xend, yend = yend, colour = colour, linewidth = 0.8,
                                                   arrow = arrow(length = unit(0.10, "inches"), type = "closed"))
p1 <- ggplot() + coord_cartesian(xlim = c(0, 10), ylim = c(0, 6.3), expand = FALSE) +
  box(0.4, 4.8, 3.7, 6.0, "#F7F7F7", C_GREY) + box(5.2, 9.6, 3.7, 6.0, "#F7F7F7", C_GREY) +
  annotate("text", x = 2.6, y = 5.65, label = "Float 20", fontface = "bold", size = 4.6) +
  annotate("text", x = 7.4, y = 5.65, label = "Floats 17-21", fontface = "bold", size = 4.6) +
  annotate("text", x = 2.6, y = 4.95, label = "I/E tally (WDF20): arrivals A, presence P", size = 3.7, colour = C_DER) +
  annotate("text", x = 2.6, y = 4.30, label = "gear count G20", size = 3.7) +
  annotate("text", x = 7.4, y = 4.95, label = "no tally", size = 3.7, colour = C_GREY, fontface = "italic") +
  annotate("text", x = 7.4, y = 4.30, label = "gear count G17, about 0.20 x G20", size = 3.7) +
  box(0.4, 4.8, 0.4, 2.5, "#E8F1F8", C_DER) + box(5.2, 9.6, 0.4, 2.5, "#FBEADF", C_FIT) +
  annotate("text", x = 2.6, y = 2.12, label = "The derivation", fontface = "bold", size = 4.3, colour = C_DER) +
  annotate("text", x = 2.6, y = 1.35, label = "turnover = A / P at the count hours\nFloat 20 over Float 20", size = 3.7, lineheight = 1.05) +
  annotate("text", x = 2.6, y = 0.68, label = "sets the PE's turnover and the fit's prior centre", size = 3.2, colour = "#444444") +
  annotate("text", x = 7.4, y = 2.12, label = "The fit's I/E likelihood", fontface = "bold", size = 4.3, colour = C_FIT) +
  annotate("text", x = 7.4, y = 1.35, label = "A compared with lambda x L,\nlambda calibrated to (G20 + G17) / R_G", size = 3.7, lineheight = 1.05) +
  annotate("text", x = 7.4, y = 0.68, label = "Float 20 over both floats: reads about 18% low", size = 3.2, colour = "#444444") +
  arr(2.0, 3.7, 2.0, 2.55, C_DER) + arr(3.6, 3.7, 6.2, 2.55, C_FIT) + arr(7.9, 3.7, 7.9, 2.55, C_FIT) +
  annotate("text", x = 1.85, y = 3.12, label = "A and P", size = 3.3, colour = C_DER, hjust = 1) +
  annotate("text", x = 5.1, y = 3.33, label = "A", size = 3.3, colour = C_FIT, hjust = 0) +
  annotate("text", x = 8.05, y = 3.12, label = "G20 + G17", size = 3.3, colour = C_FIT, hjust = 0) +
  labs(title = "The two readings of the shore turnover use the same tally, against different things") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13.5, margin = margin(b = 8)), plot.title.position = "plot",
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(12, 12, 10, 12))
save_fig(p1, "fig1_what_each_reading_compares.png", 8, 5.0)

# ---- Figure 2: the same 16 days read twice ---------------------------------------------------
d2 <- d26 |> arrange(event_date) |> mutate(day = factor(format(event_date, "%Y-%m-%d"), levels = rev(format(sort(event_date), "%Y-%m-%d"))))
gm <- function(x) exp(mean(log(x)))
p2 <- ggplot(d2) +
  geom_vline(xintercept = gm(d2$tau_derivation), colour = C_DER, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = gm(d2$tau_likelihood), colour = C_FIT, linetype = "dashed", linewidth = 0.5) +
  geom_segment(aes(x = tau_likelihood, xend = tau_derivation, y = day, yend = day), colour = C_LIGHT, linewidth = 0.9) +
  geom_point(aes(x = tau_derivation, y = day, colour = "The derivation: arrivals / the tally's presence at the count hours"), size = 2.6) +
  geom_point(aes(x = tau_likelihood, y = day, colour = "The fit's likelihood: arrivals / the fitted level"), size = 2.6) +
  annotate("text", x = gm(d2$tau_derivation) * 1.03, y = 16.75, label = sprintf("geometric mean %.3f", gm(d2$tau_derivation)), colour = C_DER, hjust = 0, size = 3.3) +
  annotate("text", x = gm(d2$tau_likelihood) / 1.03, y = 16.75, label = sprintf("geometric mean %.3f", gm(d2$tau_likelihood)), colour = C_FIT, hjust = 1, size = 3.3) +
  scale_x_log10(breaks = c(0.75, 1, 1.5, 2, 3, 4, 6), limits = c(0.7, 7.5)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.3))) +
  scale_colour_manual(values = c("The derivation: arrivals / the tally's presence at the count hours" = C_DER, "The fit's likelihood: arrivals / the fitted level" = C_FIT)) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(title = "The same 16 I/E days read twice: the fit reads the turnover 22% lower",
       subtitle = "2025-26 shore all-gear fit (readiness render B76). One row per I/E day; turnover on a log scale.",
       x = "turnover (trips per counted gear slot)", y = NULL,
       caption = "Source: 05_output/shore_turnover_d50_2026-11-18_by_day.csv (tau_derivation, tau_likelihood).") +
  theme_brief()
save_fig(p2, "fig2_two_readings_by_day.png", 8, 5.6)

# ---- Figure 3: where the gap comes from -------------------------------------------------------
d3 <- tibble(
  part  = c("The whole gap between the two readings", "Float coverage: both floats' gear over Float 20's", "Everything else together",
            "... the fitted level over the day's own counts", "... Float 20 gear count over the tally's presence", "... count timing on the I/E days"),
  value = c(f26$g, f26$f_cover, f26$f_all_but_cover, f26$f_smooth, f26$f_gear_B, f26$f_time_B),
  se    = c(f26$g_se, f26$f_cover_se, f26$f_all_but_cover_se, f26$f_smooth_se, f26$f_gear_B_se, f26$f_time_B_se),
  kind  = c("total", "cover", "rest", "sub", "sub", "sub")) |>
  mutate(part = factor(part, levels = rev(part)), label = sprintf("%+.3f (SE %.3f)", value, se))
p3 <- ggplot(d3, aes(y = part)) +
  geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.4) +
  geom_col(aes(x = value, fill = kind), width = 0.62) +
  geom_errorbar(aes(xmin = value - 2 * se, xmax = value + 2 * se), width = 0.18, colour = "#222222", linewidth = 0.45, orientation = "y") +
  geom_text(aes(x = pmax(value + 2 * se, 0) + 0.02, label = label), hjust = 0, size = 3.4, colour = "#222222") +
  scale_fill_manual(values = c(total = "#4D4D4D", cover = C_FIT, rest = C_GREY, sub = C_LIGHT), guide = "none") +
  scale_x_continuous(limits = c(-0.36, 0.72), breaks = seq(-0.3, 0.5, 0.1), labels = function(x) sprintf("%.1f", x)) +
  labs(title = "Coverage is the one part of the gap that is not noise",
       subtitle = "Mean log ratio over the 16 I/E days of the 2025-26 shore all-gear fit. Whiskers are 2 standard errors.\nThe last three rows sum to the third; the second and third sum to the first.",
       x = "mean log ratio (0.20 is a factor of 1.22)", y = NULL,
       caption = "Source: 05_output/shore_turnover_d50_2026-11-18_factors.csv (g, f_cover, f_all_but_cover, f_smooth, f_gear_B, f_time_B).") +
  theme_brief() + theme(panel.grid.major.y = element_blank())
save_fig(p3, "fig3_where_the_gap_comes_from.png", 8.6, 4.6)

# ---- the shore counts of the two seasons (for Figures 4 and 6) -------------------------------
season_window <- function(y) list(
  est_date_start = sprintf("%d-09-16", y), est_date_end = sprintf("%d-09-15", y + 1), season_filter = sprintf("%d-%02d", y, (y + 1) %% 100),
  pot_closures = NULL, census_windows = NULL,
  pot_closure_start = sprintf("%d-09-16", y), pot_closure_end = sprintf("%d-11-30", y), pot_open_date = sprintf("%d-12-01", y),
  census_start_date = sprintf("%d-12-01", y), census_end_date = sprintf("%d-01-03", y + 1), commercial_opener = sprintf("%d-01-04", y + 1),
  run_tag = sprintf("season-%d-%02d", y, (y + 1) %% 100))
with_window <- function(w) { p <- run_config; for (k in names(w)) p[k] <- list(w[[k]])
  p$boat_require_gear_time <- TRUE; p$crabbing_holiday_dates <- quietly(read_crabbing_holidays(p)); p }
counts <- bind_rows(lapply(2024:2025, function(y) quietly(fetch_crab_data(with_window(season_window(y))))$shore_effort |>
                             mutate(season = sprintf("%d-%02d", y, (y + 1) %% 100))))
ie_days <- as.Date(unique(by_count$event_date))

# ---- Figure 4: the two floats ------------------------------------------------------------------
d4 <- counts |> filter(f17_source == "observed")
r4 <- d4 |> group_by(season) |> summarise(n = n(), r = sum(f17_gear) / sum(f20_gear), xmax = max(f20_gear), .groups = "drop") |>
  mutate(lab = sprintf("%d observed pairs\nFloats 17-21 = %.3f x Float 20\nboth floats = %.3f x Float 20", n, r, 1 + r))
p4 <- ggplot(d4, aes(f20_gear, f17_gear)) +
  geom_point(colour = C_GREY, alpha = 0.35, size = 1.6, position = position_jitter(width = 0.25, height = 0.25, seed = 1)) +
  geom_abline(data = r4, aes(slope = r, intercept = 0), colour = C_FIT, linewidth = 0.9) +
  geom_text(data = r4, aes(x = 2, y = 74, label = lab), hjust = 0, vjust = 1, size = 3.4, colour = "#222222", lineheight = 1.05) +
  facet_wrap(~season) + coord_cartesian(ylim = c(0, 75)) +
  labs(title = "Floats 17-21 hold about a fifth as much gear as Float 20, in both seasons",
       subtitle = "Every Float 20 gear count with a Float 17-21 count taken beside it. The line is the ratio of sums.",
       x = "gear counted on Float 20", y = "gear counted on Floats 17-21",
       caption = "Source: effort_combined.xlsx through fetch_crab_data() and shore_dock_counts(); f17_source == \"observed\".") +
  theme_brief()
save_fig(p4, "fig4_the_two_floats.png", 8, 4.6)

# ---- Figure 5: is a Float 20 gear count the tally's presence? ---------------------------------
d5all <- by_count |> filter(season %in% c("2024-25", "2025-26"), is.finite(presence_balance))
d5 <- d5all |> filter(presence_balance > 0)                # a ratio needs someone present; the bars below keep every count
bins <- d5all |> mutate(cls = factor(hour_class, levels = c("before 10:00", "10:00 to 15:00", "15:00 and after"))) |>
  group_by(season, cls) |> summarise(n = n(), gpp = sum(f20_gear) / sum(presence_balance), x0 = min(count_hour), x1 = max(count_hour), .groups = "drop") |>
  mutate(x0 = c(`before 10:00` = 5.5, `10:00 to 15:00` = 10, `15:00 and after` = 15)[as.character(cls)],
         x1 = c(`before 10:00` = 10, `10:00 to 15:00` = 15, `15:00 and after` = 21)[as.character(cls)],
         lab = sprintf("%.2f (%d)", gpp, n))
rg <- by_season |> filter(season %in% c("2024-25", "2025-26")) |> transmute(season, R_G = R_G_all_gear_fit, both = gpp_if_tally_covered_both_floats)
p5 <- ggplot(d5, aes(count_hour, gear_per_present)) +
  annotate("rect", xmin = 10, xmax = 15, ymin = 0, ymax = Inf, fill = C_BAND) +
  geom_hline(data = rg, aes(yintercept = R_G), colour = C_DER, linetype = "dashed", linewidth = 0.6) +
  geom_hline(data = rg, aes(yintercept = both), colour = C_GREEN, linetype = "dotted", linewidth = 0.8) +
  geom_point(aes(size = presence_balance), colour = C_GREY, alpha = 0.45) +
  geom_segment(data = bins, aes(x = x0, xend = x1, y = gpp, yend = gpp), colour = C_FIT, linewidth = 1.3) +
  geom_text(data = bins, aes(x = (x0 + x1) / 2, y = gpp * 1.22, label = lab), colour = C_FIT, size = 3.3, fontface = "bold") +
  geom_text(data = rg, aes(x = 5.6, y = R_G * 1.10, label = sprintf("one float: R_G %.2f", R_G)), colour = C_DER, hjust = 0, size = 3.1) +
  geom_text(data = rg, aes(x = 5.6, y = both * 0.90, label = sprintf("both floats: %.2f", both)), colour = C_GREEN, hjust = 0, size = 3.1) +
  facet_wrap(~season) +
  scale_y_log10(breaks = c(0.5, 1, 2, 4), limits = c(0.4, 5)) + scale_x_continuous(breaks = seq(6, 20, 2), labels = function(h) sprintf("%02d:00", h), limits = c(5.5, 21.2)) +
  scale_size_area(max_size = 6, guide = "none") +
  labs(title = "At the regular count hours the tally reads as Float 20's alone; early in the day gear outruns it",
       subtitle = "Gear on Float 20 per crabber the tally has present, at each gear count taken on an I/E day (log scale).\nPoint area is presence. Red bars: the ratio of sums in each time band, with the number of counts.\nShaded: the regular count hours. Dashed: what a Float 20 tally implies. Dotted: what a two-float tally would imply.",
       x = "time of the gear count", y = "gear per present crabber",
       caption = "Source: 05_output/shore_turnover_d50_2026-11-18_by_count.csv and _by_season.csv. Presence is the tally's running balance at the count time.\nA count taken with nobody yet on the float has no point but enters its bar (one in each season, both before 10:00).") +
  theme_brief(11.5)
save_fig(p5, "fig5_gear_count_against_tally_presence.png", 9, 5.0)

# ---- Figure 6: the count schedule on I/E days against the other days -------------------------
d6 <- counts |> filter(season == "2025-26", count_sequence <= run_config$bss_max_count_seq) |>
  mutate(grp = ifelse(event_date %in% ie_days, "I/E days", "all other count days"), hbin = floor(count_hour)) |>
  count(grp, hbin) |> group_by(grp) |> mutate(share = n / sum(n), N = sum(n)) |> ungroup()
lab6 <- d6 |> group_by(grp) |> summarise(N = first(N), off = sum(n[hbin < 10 | hbin >= 15]), .groups = "drop") |>
  mutate(lab = sprintf("%s: %d counts, %d (%.0f%%) outside 10:00 to 15:00", grp, N, off, 100 * off / N))
d6 <- d6 |> left_join(lab6 |> select(grp, lab), by = "grp") |> mutate(lab = factor(lab, levels = lab6$lab[order(lab6$grp != "I/E days")]))
p6 <- ggplot(d6, aes(hbin + 0.5, share)) +
  annotate("rect", xmin = 10, xmax = 15, ymin = 0, ymax = Inf, fill = C_BAND) +
  geom_col(aes(fill = hbin < 10 | hbin >= 15), width = 0.9) +
  facet_wrap(~lab, ncol = 1) +
  scale_fill_manual(values = c(`TRUE` = C_FIT, `FALSE` = C_GREY), guide = "none") +
  scale_x_continuous(breaks = seq(6, 20, 2), labels = function(h) sprintf("%02d:00", h), limits = c(5, 21)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "On an I/E day the gear counts follow the tally, not the regular schedule",
       subtitle = "2025-26 shore gear counts by hour, as a share of each group's counts. Shaded: 10:00 to 15:00.",
       x = "hour of the gear count", y = "share of counts",
       caption = "Source: effort_combined.xlsx through fetch_crab_data(); count sequences 1 to 3, the ones the fit reads.") +
  theme_brief() + theme(strip.text = element_text(hjust = 0, face = "bold"))
save_fig(p6, "fig6_count_schedule.png", 8, 5.2)

# ---- Figure 7: one I/E day ---------------------------------------------------------------------
day7 <- as.Date("2026-04-03")
tally <- read_input_workbook(here("04_input_files", run_config$ie_data_file), sheet = run_config$ie_sheet) |>
  mutate(event_date = as.Date(date)) |> filter(location_name == run_config$ie_shore_location, event_date == day7)
tally$hour <- .ie_hour_of(tally$time)
tally <- tally |> arrange(hour) |> mutate(on = replace_na(as.numeric(crabbers_on), 0), off = replace_na(as.numeric(crabbers_off), 0), balance = cumsum(on - off))
c7 <- counts |> filter(event_date == day7, count_sequence <= run_config$bss_max_count_seq) |>
  mutate(lab = sprintf("gear count %d\n%s\n%d on Float 20", count_sequence, hm(count_hour), f20_gear))
p7 <- ggplot(tally, aes(hour, balance)) +
  annotate("rect", xmin = 10, xmax = 15, ymin = 0, ymax = Inf, fill = C_BAND) +
  geom_step(colour = C_DER, linewidth = 0.9) +
  geom_vline(data = c7, aes(xintercept = count_hour), colour = C_FIT, linewidth = 0.7, linetype = "dashed") +
  geom_label(data = c7, aes(x = count_hour, y = c(46, 56, 14), label = lab, hjust = c(0, 0.5, 1)), colour = C_FIT, size = 3.2, lineheight = 0.95,
             fill = "white", label.size = 0.3) +
  annotate("text", x = 12.5, y = 3.2, label = "the regular count hours", colour = "#555555", size = 3.3) +
  scale_x_continuous(breaks = seq(8, 18, 2), labels = function(h) sprintf("%02d:00", h)) +
  labs(title = "One I/E day: gear counts at the tally's start, at 11:11 and at its end",
       subtitle = sprintf("%s. Blue: crabbers the WDF20 tally has on Float 20 (arrivals less departures so far).\n%d arrivals in all. The afternoon, when most were present, has no gear count.", format(day7), sum(tally$on)),
       x = NULL, y = "crabbers present on Float 20",
       caption = "Source: ingress_egress.xlsx (WDF20) and effort_combined.xlsx.\nThe workbook's own flow column reads 0 on the first and last rows of this day; the line is the running sum of on less off.") +
  theme_brief()
save_fig(p7, "fig7_one_ie_day.png", 8, 4.6)

# ---- Figure 8: what each option does (desk approximation) --------------------------------------
opt <- function(f, season) tibble(
  season = season,
  option = c("(a) leave it", "(b) drop the shore I/E stream", "(c) Float 20 share in the prediction", "(d) observe the day's own turnover"),
  tau = c(f$tau_bar_fit, f$prior_mu, f$tau_bar_fit * f$tau_bar_desk_cover_over_fit, f$tau_bar_fit * f$tau_bar_desk_own_over_fit),
  centre = f$prior_mu, fitted = f$tau_bar_fit)
d8 <- bind_rows(opt(f25, "2024-25, the authoritative run (B70): 4 I/E days"), opt(f26, "2025-26, the readiness render (B76): 16 I/E days")) |>
  mutate(option = factor(option, levels = rev(unique(option))), change = tau / fitted - 1,
         lab = ifelse(abs(change) < 1e-9, sprintf("%.3f (as fitted)", tau), sprintf("%.3f (%+.1f%%)", tau, 100 * change)))
p8 <- ggplot(d8, aes(y = option)) +
  geom_vline(aes(xintercept = centre), colour = C_DER, linetype = "dashed", linewidth = 0.6) +
  geom_segment(aes(x = fitted, xend = tau, yend = option), colour = C_LIGHT, linewidth = 1.1) +
  geom_point(aes(x = tau, colour = option == "(a) leave it"), size = 3.2) +
  geom_text(aes(x = tau, label = lab), vjust = -1.1, size = 3.3, colour = "#222222") +
  geom_text(data = distinct(d8, season, centre), aes(x = centre, y = 0.42, label = sprintf("derived centre %.3f", centre)), colour = C_DER, size = 3.1, hjust = -0.06) +
  facet_wrap(~season, ncol = 1) +
  scale_colour_manual(values = c(`TRUE` = C_FIT, `FALSE` = "#222222"), guide = "none") +
  scale_x_continuous(limits = c(2.2, 2.72), breaks = seq(2.2, 2.7, 0.1)) + scale_y_discrete(expand = expansion(add = c(0.95, 0.95))) +
  labs(title = "What each option would do to the fitted shore turnover (desk approximation, not a fit)",
       subtitle = "The shore all-gear fit's tau_bar; that fit's effort and catch scale with it.\nPercent is the change from the fitted value. Dashed: the derived centre, which is also the PE's turnover.",
       x = "tau_bar, shore all-gear fit", y = NULL,
       caption = "Source: 05_output/shore_turnover_d50_2026-11-18_factors.csv (tau_bar_fit, prior_mu, tau_bar_desk_cover_over_fit, tau_bar_desk_own_over_fit).") +
  theme_brief() + theme(strip.text = element_text(hjust = 0, face = "bold"), panel.grid.major.y = element_blank())
save_fig(p8, "fig8_what_each_option_does.png", 8.4, 5.6)
cat("done\n")
