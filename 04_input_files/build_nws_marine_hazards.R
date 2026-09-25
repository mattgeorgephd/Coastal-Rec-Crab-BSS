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
# build_nws_marine_hazards.R  (2026-09-25)
#
# Builds 04_input_files/nws_marine_hazards.xlsx (sheet "data"): every National Weather
# Service marine watch / warning / advisory (VTEC event) issued for the marine zones that
# cover Westport, with the time each took effect and expired, from the Iowa Environmental
# Mesonet's VTEC archive (https://mesonet.agron.iastate.edu/info/datasets/vtec.html, parsed
# from the NWS text stream by pyIEM; NCEI keeps the raw text). Read by
# 03_R_functions/bss_marine_hazard_covariates.R, which turns it into the per-day
# Small-Craft-Advisory-or-higher flags the marine hazard effort covariates use
# (run_config.R section 4.4b; off in production).
#
# THE ZONES (NWS Seattle, https://www.weather.gov/marine/sewmz):
#   PZZ110  Grays Harbor Bar                                        (the bar itself)
#   PZZ156  Coastal Waters from Point Grenville to Cape Shoalwater out 10 NM
# PZZ176 (10 to 60 NM) is offshore of anything a recreational crabber crosses and is not
# pulled by default. Add zones with --zones.
#
# THE CODES. An SCA is VTEC phenomena "SC", significance "Y". BEFORE 2019-12-03 (Service
# Change Notice 19-83) the bar, seas and wind advisories had their own codes, RB.Y (Small
# Craft Advisory for Rough Bar), SW.Y (for Hazardous Seas) and SI.Y (for Winds); a multi-year
# series must read them as SCAs, and run_config's marine_hazard_codes does by default. The
# gale (GL.W), storm (SR.W), hazardous seas (SE.W) and hurricane-force (HF.W) warnings
# supersede an SCA and are read as "SCA or higher". Watches (*.A), marine dense fog (MF.Y)
# and special marine warnings (MA.W) are carried but not counted. The archive holds the bar
# zone from 2008-08-16 (its first RB.Y is 2008-08-20).
#
# WHAT THE JSON CARRIES, and the two times that matter. The service
#   https://mesonet.agron.iastate.edu/json/vtec_events_byugc.py?ugc=PZZ110&sdate=...&edate=...
# returns one row per VTEC segment: `issue` is when the hazard TAKES EFFECT, `expire` when it
# ends, `product_id` (yyyymmddHHMM-KSEW-WHUS76-MWWSEW) carries the time the product was
# ANNOUNCED. A row whose expire is at or before its issue was cancelled before it began; it is
# kept with in_effect = FALSE. The same eventid can appear in several rows (a re-issue after
# a gap); each row is its own period, so rows are never de-duplicated by eventid.
#
# COLUMNS WRITTEN (sheet "data"; times as ISO text, UTC)
#   ugc, zone_name, phenomena, significance, hazard, eventid, start_utc, end_utc,
#   product_id, product_issued_utc, in_effect, pull_start, pull_end, pulled_at, source
# pull_start / pull_end are the window that was REQUESTED for that zone, on every row: they
# are the archive's coverage, which the reader checks against the estimation window (the
# absence of an event is not evidence of coverage).
#
# PROVENANCE OF THE COMMITTED WORKBOOK (2026-09-25). The IEM host was not reachable from the
# session that built the covariate machinery, so the committed workbook was rebuilt with
# --from-csv from raw/nws_marine_hazards_iem_transcription_2026-09-25.csv: the JSON rows
# transcribed twice, in different window partitions, with the two copies agreeing row for
# row (147 PZZ110 rows, 470 PZZ156 rows; 2023-01-01 to 2026-09-26). Its `source` column says
# so. A direct pull with this script replaces it; compare_with_previous() prints the diff.
#
# USAGE (from the repository root; needs jsonlite, network access to mesonet.agron.iastate.edu)
#   Rscript 04_input_files/build_nws_marine_hazards.R                      # 2008-01-01 to today
#   Rscript 04_input_files/build_nws_marine_hazards.R --start 2023-01-01 --end 2026-09-26
#   Rscript 04_input_files/build_nws_marine_hazards.R --zones PZZ110,PZZ156,PZZ176
#   Rscript 04_input_files/build_nws_marine_hazards.R --from-csv raw/<file>.csv --start ... --end ...
# The pull is one request per zone per calendar year (the service is happiest with short
# windows), retried once. Nothing in a model run calls the network: the run reads the
# committed workbook only.
###############################################################################

source(file.path(if (dir.exists("04_input_files")) "04_input_files" else ".", "build_helpers.R"))
root <- build_root()
out_path <- file.path(root, "04_input_files", "nws_marine_hazards.xlsx")

# ---- arguments ---------------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.arg <- function(flag, default = NULL) { i <- which(.args == flag); if (length(i) && i[1] < length(.args)) .args[i[1] + 1] else default }
ZONES    <- strsplit(.arg("--zones", "PZZ110,PZZ156"), ",", fixed = TRUE)[[1]]
START    <- as.Date(.arg("--start", "2008-01-01"))
END      <- as.Date(.arg("--end", format(Sys.Date(), "%Y-%m-%d")))
FROM_CSV <- .arg("--from-csv", NULL)
if (is.na(START) || is.na(END) || END < START) stop("--start / --end must be ISO dates with start <= end")
# --from-csv rows carry no window of their own, and the workbook's pull_start / pull_end are
# what the reader trusts as COVERAGE: a re-import without an explicit window would stamp the
# defaults (2008-01-01 to today) on rows that may begin years later, and a run on a window the
# rows do not cover would then get silently-zero flags instead of the coverage stop.
if (!is.null(FROM_CSV) && (is.null(.arg("--start")) || is.null(.arg("--end"))))
  stop("--from-csv needs an explicit --start and --end: the window the rows were pulled for, ",
       "which becomes the coverage the reader enforces (the committed transcription is 2023-01-01 to 2026-09-26)")

ZONE_NAMES <- c(PZZ110 = "Grays Harbor Bar",
                PZZ153 = "Coastal Waters from James Island to Point Grenville out 10 NM",
                PZZ156 = "Coastal Waters from Point Grenville to Cape Shoalwater out 10 NM",
                PZZ173 = "Waters from James Island to Point Grenville from 10 to 60 NM",
                PZZ176 = "Waters from Point Grenville to Cape Shoalwater from 10 to 60 NM")
HAZARD_NAMES <- c("SC.Y" = "Small Craft Advisory", "RB.Y" = "Small Craft Advisory for Rough Bar",
                  "SW.Y" = "Small Craft Advisory for Hazardous Seas", "SI.Y" = "Small Craft Advisory for Winds",
                  "GL.W" = "Gale Warning", "GL.A" = "Gale Watch", "SR.W" = "Storm Warning", "SR.A" = "Storm Watch",
                  "SE.W" = "Hazardous Seas Warning", "SE.A" = "Hazardous Seas Watch",
                  "HF.W" = "Hurricane Force Wind Warning", "HF.A" = "Hurricane Force Wind Watch",
                  "MF.Y" = "Dense Fog Advisory (marine)", "MA.W" = "Special Marine Warning",
                  "MA.S" = "Marine Statement", "UP.Y" = "Freezing Spray Advisory", "UP.W" = "Heavy Freezing Spray Warning")

# ---- the pull ----------------------------------------------------------------------------
# One zone, one window -> the JSON's `events` array as a data frame (zero rows when empty).
# Pure: takes the URL text through `fetch` so the parser can be tested without a network.
iem_vtec_url <- function(ugc, sdate, edate)
  sprintf("https://mesonet.agron.iastate.edu/json/vtec_events_byugc.py?ugc=%s&sdate=%s&edate=%s",
          ugc, format(as.Date(sdate), "%Y-%m-%d"), format(as.Date(edate), "%Y-%m-%d"))

parse_iem_vtec <- function(json_text, ugc) {
  x <- jsonlite::fromJSON(json_text, simplifyVector = TRUE)
  ev <- x$events
  if (is.null(ev) || (is.data.frame(ev) && !nrow(ev)) || (is.list(ev) && !length(ev)))
    return(tibble(ugc = character(0), phenomena = character(0), significance = character(0), eventid = integer(0),
                  issue = character(0), expire = character(0), product_id = character(0), name = character(0)))
  ev <- as.data.frame(ev, stringsAsFactors = FALSE)
  pick <- function(col, default = NA_character_) if (col %in% names(ev)) as.character(ev[[col]]) else rep(default, nrow(ev))
  tibble(ugc = toupper(pick("ugc", ugc)), phenomena = pick("phenomena"), significance = pick("significance"),
         eventid = suppressWarnings(as.integer(pick("eventid"))), issue = pick("issue"), expire = pick("expire"),
         product_id = pick("product_id"), name = pick("name"))
}

fetch_iem_vtec <- function(ugc, sdate, edate, fetch = function(u) paste(readLines(url(u), warn = FALSE), collapse = "\n")) {
  u <- iem_vtec_url(ugc, sdate, edate)
  txt <- tryCatch(fetch(u), error = function(e) { Sys.sleep(3); fetch(u) })   # one retry
  parse_iem_vtec(txt, ugc)
}

pull_zone <- function(ugc, sdate, edate) {
  yrs <- seq(as.integer(format(sdate, "%Y")), as.integer(format(edate, "%Y")))
  parts <- lapply(yrs, function(y) {
    s <- max(sdate, as.Date(sprintf("%d-01-01", y))); e <- min(edate, as.Date(sprintf("%d-12-31", y)))
    say("  %s %s to %s ...", ugc, s, e)
    fetch_iem_vtec(ugc, s, e)
  })
  bind_rows(parts) |> distinct(ugc, phenomena, significance, eventid, issue, expire, product_id, .keep_all = TRUE)
}

# ---- the workbook ------------------------------------------------------------------------
# `raw` has the JSON columns (ugc, phenomena, significance, eventid, issue, expire, product_id,
# optional name). Times are ISO-8601 UTC text as the service emits them.
build_nws_workbook <- function(raw, zones, sdate, edate, source_text) {
  parse_iso <- function(x) as.POSIXct(sub("Z$", "", x), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  st <- parse_iso(raw$issue); en <- parse_iso(raw$expire)
  report_unparsed("issue time", raw$issue, st); report_unparsed("expire time", raw$expire, en)
  ps <- paste0(raw$phenomena, ".", raw$significance)
  prod_t <- as.POSIXct(substr(as.character(raw$product_id), 1, 12), format = "%Y%m%d%H%M", tz = "UTC")
  iso <- function(t) ifelse(is.na(t), NA_character_, format(t, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  tibble(
    ugc = raw$ugc,
    zone_name = ifelse(raw$ugc %in% names(ZONE_NAMES), unname(ZONE_NAMES[raw$ugc]), NA_character_),
    phenomena = raw$phenomena, significance = raw$significance,
    hazard = ifelse(ps %in% names(HAZARD_NAMES), unname(HAZARD_NAMES[ps]),
                    if ("name" %in% names(raw)) as.character(raw$name) else NA_character_),
    eventid = raw$eventid,
    start_utc = iso(st), end_utc = iso(en),
    product_id = as.character(raw$product_id), product_issued_utc = iso(prod_t),
    in_effect = as.integer(!is.na(st) & !is.na(en) & en > st),
    pull_start = format(sdate, "%Y-%m-%d"), pull_end = format(edate, "%Y-%m-%d"),
    pulled_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    source = source_text
  ) |>
    filter(ugc %in% zones) |>
    arrange(ugc, start_utc, product_id)
}

# ---- main --------------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  if (!is.null(FROM_CSV)) {
    csv <- if (file.exists(FROM_CSV)) FROM_CSV else file.path(root, "04_input_files", FROM_CSV)
    if (!file.exists(csv)) stop("--from-csv file not found: ", csv)
    raw <- suppressWarnings(readr::read_csv(csv, col_types = readr::cols(.default = readr::col_character()), progress = FALSE))
    need <- c("ugc", "phenomena", "significance", "eventid", "issue", "expire", "product_id")
    if (!all(need %in% names(raw))) stop("--from-csv needs columns: ", paste(need, collapse = ", "))
    raw <- raw |> mutate(ugc = toupper(ugc), eventid = suppressWarnings(as.integer(eventid)))
    src <- sprintf("IEM vtec_events_byugc.py rows re-imported from %s (not a live pull)", basename(csv))
    say("Rebuilding from %s: %d rows", basename(csv), nrow(raw))
  } else {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("install.packages('jsonlite') to pull the archive")
    say("Pulling the IEM VTEC archive for %s, %s to %s", paste(ZONES, collapse = ", "), START, END)
    raw <- bind_rows(lapply(ZONES, pull_zone, sdate = START, edate = END))
    src <- "IEM vtec_events_byugc.py (https://mesonet.agron.iastate.edu/info/datasets/vtec.html)"
  }
  wb <- build_nws_workbook(raw, ZONES, START, END, src)
  cat("\nEvents by zone x hazard:\n"); print(as.data.frame(wb |> count(ugc, phenomena, significance, hazard)), row.names = FALSE)
  cat(sprintf("\nCancelled-before-start rows (in_effect = 0): %d\n", sum(wb$in_effect == 0)))
  if (file.exists(out_path)) {
    old <- suppressWarnings(readxl::read_excel(out_path, sheet = "data", guess_max = 1e5, col_types = "text"))
    old$in_effect <- suppressWarnings(as.integer(old$in_effect))
    # the key carries the hazard too: one product can issue a watch AND a warning for the
    # same zone at the same instant (a Gale Watch upgraded to a Gale Warning), so
    # ugc + product_id + start_utc alone is not unique
    key <- c("ugc", "phenomena", "significance", "product_id", "start_utc")
    cat("\nComparison with the previous workbook on ugc + phenomena + significance + product_id + start_utc:\n")
    compare_with_previous(old, wb |> mutate(across(everything(), as.character)) |> mutate(in_effect = as.integer(in_effect)),
                          key, c("eventid", "end_utc", "in_effect"))
    only_old <- anti_join(old, wb |> mutate(across(everything(), as.character)), by = key)
    only_new <- anti_join(wb |> mutate(across(everything(), as.character)), old, by = key)
    say("  rows only in the previous workbook: %d; only in the rebuilt one: %d", nrow(only_old), nrow(only_new))
  }
  write_data_sheet(wb, out_path)
}
