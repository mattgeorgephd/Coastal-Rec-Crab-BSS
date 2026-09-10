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
# build_helpers.R  (2026-09-10)
#
# Shared helpers for the input-workbook builders in this folder. Sourced (not run) by
# build_interview_combined.R, build_effort_combined.R, build_sampler_shifts.R,
# build_comm_charter_tally.R, build_charter_trips.R and build_crabbing_holidays.R.
#
# THE SOURCES. Since 2026-09-10 every model workbook is built from the PER-SEASON creel
# workbooks in raw/ (<YYYY><YY>_rec_crab_harvest_data.xlsx: the creel database exported
# one fishery season at a time, every sheet). They replace the two pasted exports of
# 2026-09-09 (interviewdata20222026.xlsx, surveydata20222026.xlsx), whose column
# pasting had left the 2023-24 gear count under a differently-capitalised header and
# whose 2025-26 rows stopped at 2026-08-01. The season workbooks differ from each other
# in header spelling and case ("Number of Gear" / "Number Of Gear", "creel Area",
# "Crab Released?" / "Crab released?"), in which columns exist (the 2024-25 protocol
# split the trailer and vehicle counts, added the released-crab fields and the
# vehicle count), and in how Excel stored dates and clock times (a date cell, ISO
# text, M/D/YYYY text; a time cell, "H:MM:SS" text). Everything is therefore read as
# TEXT and parsed here by one rule, so a builder never depends on readxl's type guess.
#
# CONVENTIONS every builder follows (the reader contract in ../README.md):
#   - one sheet "data"; dates as ISO yyyy-mm-dd text; clock times as text ("HH:MM" or
#     "HH:MM:SS", stated per column); a `season` column with the fishery-season label
#     (Sep 16 boundary), taken from the workbook's file name and checked against the
#     row's date.
#   - source values are never corrected silently: a value that cannot be parsed becomes
#     NA and is counted in the builder's report; a row that must be held out of the
#     arithmetic gets a qc_flag, not a deletion.
###############################################################################

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tibble); library(stringr); library(purrr); library(writexl)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- repository root and the season workbooks ----------------------------------------
build_root <- function() {
  r <- getwd()
  for (up in c(".", "..", "../..")) {
    cand <- normalizePath(file.path(r, up), mustWork = FALSE)
    if (dir.exists(file.path(cand, "04_input_files", "raw"))) return(cand)
  }
  stop("run the builders from the repository root (04_input_files/raw not found from ", r, ")")
}

# tibble(tag, season, path) for every raw/<YYYY><YY>_rec_crab_harvest_data.xlsx, oldest first
season_workbooks <- function(root = build_root(), raw_dir = file.path(root, "04_input_files", "raw")) {
  files <- list.files(raw_dir, pattern = "^[0-9]{4}_rec_crab_harvest_data\\.xlsx$", full.names = TRUE)
  if (!length(files)) stop("no season workbooks (raw/<YYYY><YY>_rec_crab_harvest_data.xlsx) under ", raw_dir)
  tag <- str_extract(basename(files), "^[0-9]{4}")
  y0 <- 2000L + as.integer(substr(tag, 1, 2)); y1 <- 2000L + as.integer(substr(tag, 3, 4))
  if (any(y1 != y0 + 1L)) stop("season workbook tag must be two consecutive years, e.g. 2425: ", paste(basename(files)[y1 != y0 + 1L], collapse = ", "))
  tibble(tag = tag, season = sprintf("%d-%02d", y0, y1 %% 100), path = files) |> arrange(tag)
}

# ---- sheets and columns by normalised name --------------------------------------------
norm_name <- function(x) {
  x <- tolower(as.character(x))
  x <- str_replace_all(x, "[?#().,/'-]", " ")
  str_squish(x)
}

# the workbook's sheet whose normalised name equals one of `wanted` (first hit), else NULL
sheet_named <- function(path, wanted) {
  sh <- readxl::excel_sheets(path)
  hit <- sh[norm_name(sh) %in% norm_name(wanted)]
  if (length(hit)) hit[1] else NULL
}

# read a sheet with every cell as text (dates come back as Excel serial-day text, times
# as day-fraction text) so the parsers below apply one rule; header row = `skip` + 1
read_sheet_text <- function(path, sheet, skip = 0) {
  df <- suppressWarnings(readxl::read_excel(path, sheet = sheet, col_types = "text", skip = skip,
                                            .name_repair = "minimal"))
  names(df) <- str_squish(names(df))
  df <- df[, nzchar(names(df)), drop = FALSE]            # unnamed spill-over columns
  df[rowSums(!is.na(df)) > 0, , drop = FALSE]             # fully blank rows
}

# the first column of `df` whose normalised name is one of `aliases`; NA vector when none
pick_col <- function(df, aliases, default = NA_character_) {
  nn <- norm_name(names(df)); want <- norm_name(aliases)
  hit <- which(nn %in% want)
  if (length(hit)) df[[hit[1]]] else rep(default, nrow(df))
}
has_col <- function(df, aliases) any(norm_name(names(df)) %in% norm_name(aliases))

# ---- parsers ------------------------------------------------------------------------
num <- function(x) { if (is.logical(x)) return(as.numeric(x)); suppressWarnings(as.numeric(as.character(x))) }

# Excel serial-day text ("45291"), ISO ("2024-01-01", with or without a time), or
# M/D/YYYY text -> Date. Anything else NA.
parse_date_any <- function(x) {
  x <- str_squish(as.character(x))
  out <- rep(as.Date(NA), length(x))
  serial <- grepl("^[0-9]{5}(\\.[0-9]+)?$", x)
  out[serial] <- as.Date(floor(as.numeric(x[serial])), origin = "1899-12-30")
  iso <- !serial & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}", x)
  out[iso] <- as.Date(substr(x[iso], 1, 10))
  us <- !serial & !iso & grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}", x)
  out[us] <- as.Date(sub("\\s.*$", "", x[us]), format = "%m/%d/%Y")
  out
}

# Day-fraction text ("0.4479"), "H:MM(:SS)" text (24-hour), or a serial datetime
# ("45291.4479") -> decimal hours; anything else NA.
parse_clock_hours <- function(x) {
  x <- str_squish(as.character(x))
  out <- rep(NA_real_, length(x))
  frac <- grepl("^[0-9]*\\.[0-9]+$", x) | grepl("^[01]$", x)
  v <- suppressWarnings(as.numeric(x[frac])); v <- (v - floor(v)) * 24
  out[frac] <- v
  hms <- !frac & grepl("^[0-9]{1,2}:[0-9]{1,2}(:[0-9]{1,2})?$", x)
  m <- str_match(x[hms], "^([0-9]{1,2}):([0-9]{1,2})(?::([0-9]{1,2}))?$")
  h <- as.numeric(m[, 2]); mi <- as.numeric(m[, 3]); s <- as.numeric(m[, 4]); s[is.na(s)] <- 0
  hh <- h + mi / 60 + s / 3600
  hh[h > 24 | mi >= 60 | s >= 60] <- NA_real_
  out[hms] <- hh
  out[!is.na(out) & (out < 0 | out >= 24)] <- NA_real_
  out
}
hhmm   <- function(hrs) ifelse(is.na(hrs), NA_character_, sprintf("%02d:%02d", floor(hrs + 1e-9), round((hrs - floor(hrs + 1e-9)) * 60)))
hhmmss <- function(hrs) {
  s <- round(hrs * 3600)
  ifelse(is.na(hrs), NA_character_, sprintf("%02d:%02d:%02d", s %/% 3600, (s %% 3600) %/% 60, s %% 60))
}

# fishery-season label by the Sep 16 boundary: a date on/after Sep 16 of year Y is "Y-(Y+1)"
season_of <- function(d) {
  y <- as.integer(format(d, "%Y")); after <- format(d, "%m%d") >= "0916"
  y0 <- ifelse(after, y, y - 1L)
  ifelse(is.na(d), NA_character_, sprintf("%d-%02d", y0, (y0 + 1L) %% 100))
}

# "S<n>" survey ids from the sheets' ID column ("S4949", "4949", "4949.0")
survey_id_of <- function(x) {
  x <- str_squish(as.character(x))
  n <- suppressWarnings(as.integer(round(num(sub("^[Ss]", "", x)))))
  ifelse(is.na(n), NA_character_, paste0("S", n))
}

# ---- creel-area names -----------------------------------------------------------------
# The 2022-23 workbook names some sites differently from 2023-24 onward (the location_key
# sheet of the 2023-24 workbook is the current vocabulary); the same mapping the 2026-09-09
# export carried. No Westport / Grays Harbor site is affected, so the run_config site
# filters match every season. Applied to interviews and effort counts alike.
area_map <- c(
  "Chinook Boat Launch" = "Chinook Boat Launch and Marina",
  "Ilwaco Boat Launch and Marina" = "Ilwaco Boat Launch",
  "Long Beach North" = "Long Beach North (North of Bay Ave.)",
  "Long Beach South" = "Long Beach South (South of Bay Ave.)",
  "Long Beach" = "Long Beach North (North of Bay Ave.)",
  "Tokeland boat Launch" = "Tokeland Boat Launch",
  "Twin Harbors Beach" = "Twin Harbors"
)
harmonise_area <- function(x) { x <- str_squish(as.character(x)); ifelse(x %in% names(area_map), unname(area_map[x]), x) }

# ---- reporting and writing -----------------------------------------------------------
say <- function(...) cat(sprintf(...), "\n", sep = "")

# count of source values that were non-blank but did not parse, per column, for the report
report_unparsed <- function(label, raw, parsed) {
  n <- sum(!is.na(raw) & nzchar(str_squish(as.character(raw))) & is.na(parsed))
  if (n) say("  NOTE %s: %d non-blank value(s) did not parse and are NA: %s", label, n,
             paste(head(unique(raw[!is.na(raw) & nzchar(str_squish(as.character(raw))) & is.na(parsed)]), 5), collapse = " | "))
  invisible(n)
}

write_data_sheet <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writexl::write_xlsx(list(data = as.data.frame(df)), path)
  say("Wrote %s: %d rows x %d columns.", path, nrow(df), ncol(df))
  invisible(path)
}

# column-by-column comparison of a rebuilt workbook with the previous one on a key
compare_with_previous <- function(old, new, keys, cols, tol = 1e-6, label = "previous workbook") {
  cols <- intersect(cols, intersect(names(old), names(new)))
  j <- inner_join(old |> select(all_of(c(keys, cols))), new |> select(all_of(c(keys, cols))),
                  by = keys, suffix = c(".old", ".new"))
  say("  %s: %d rows; rebuilt: %d rows; %d matched on %s. Differences:", label, nrow(old), nrow(new), nrow(j), paste(keys, collapse = "+"))
  for (cn in cols) {
    a <- j[[paste0(cn, ".old")]]; b <- j[[paste0(cn, ".new")]]
    if (is.numeric(b) || is.numeric(a)) { a <- num(a); b <- num(b); d <- !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < tol)) }
    else { a <- str_squish(as.character(a)); b <- str_squish(as.character(b)); d <- !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)) }
    say("    %-26s %d", cn, sum(d))
  }
  invisible(j)
}
