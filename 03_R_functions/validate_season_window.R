###############################################################################
# LOUD SEASON / WINDOW CONSISTENCY CHECK
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   Selecting a season takes TWO independent settings: the estimation window
#   (est_date_start / est_date_end) and the data filter (season_filter, matched
#   against the `season` column of the effort and interview workbooks; the crabbing
#   holidays are matched the same way and STOP on a missing season). Nothing forced
#   them to agree. Point the window at 2025-26 while season_filter still says
#   "2024-25" and every reader silently keeps only 2024-25 rows, the window filter
#   then empties them, and the run limps to a wall of PE fallbacks and empty fits
#   whose real cause, one stale string, is nowhere on screen.
#
#   The design goal is that this pipeline runs on ANY window the user selects, full
#   season, part-season, or a multi-season span (season_filter is a vector as of
#   2026-09-09). The first thing a naive run owes the user is a plain statement of
#   what data the selection actually captured, and a hard stop when it captured
#   nothing.
#
# WHAT IT DOES
#   Prints, per season label found in the (already season-filtered) frames, the date
#   range and row counts inside and outside the estimation window; STOPS when the
#   window contains zero effort AND zero interview rows; and WARNS (does not stop)
#   when the commercial/charter census window falls entirely outside the estimation
#   window, when a requested season label matched nothing at all, or when a large
#   share of the season's data lies outside the window (expected for a deliberate
#   part-season run, so it is a note, not an error).
#
# WHERE IT RUNS
#   Called at the end of fetch_crab_data(), so both drivers and every batch runner
#   get it without any of them opting in. It reads only what is already in memory;
#   cost is milliseconds.
###############################################################################

validate_season_window <- function(effort, interview, params, quiet = FALSE) {
  ws <- as.Date(params$est_date_start); we <- as.Date(params$est_date_end)
  seasons_requested <- as.character(params$season_filter)

  eff_dates <- as.Date(effort$date %||% effort$event_date)
  int_dates <- as.Date(interview$event_date %||% interview$date)
  eff_season <- as.character(effort$season %||% rep(NA_character_, length(eff_dates)))
  int_season <- as.character(interview$season %||% rep(NA_character_, length(int_dates)))

  # 2026-09-10: quiet gates the PRINTS only. The per-season warnings fire either way;
  # a safety warning silenced by a verbosity flag is no safety warning at all.
  if (!isTRUE(quiet))
    cat(sprintf("\nSeason/window check: est window %s to %s; season_filter = %s\n",
                ws, we, paste(seasons_requested, collapse = " + ")))
  for (sn in seasons_requested) {
    ed <- eff_dates[eff_season == sn]; id <- int_dates[int_season == sn]
    if (!length(ed) && !length(id)) {
      warning(sprintf("season_filter '%s' matched NO effort and NO interview rows; check the season column spelling.", sn),
              call. = FALSE)
      if (!isTRUE(quiet)) cat(sprintf("  %-10s NO ROWS MATCHED\n", sn))
      next
    }
    rng <- range(c(ed, id), na.rm = TRUE)
    ne_in <- sum(ed >= ws & ed <= we, na.rm = TRUE)
    ni_in <- sum(id >= ws & id <= we, na.rm = TRUE)
    n_out <- length(ed) + length(id) - ne_in - ni_in
    if (!isTRUE(quiet))
      cat(sprintf("  %-10s data %s to %s | in window: %d effort counts, %d interviews | outside %d%s\n",
                  sn, rng[1], rng[2], ne_in, ni_in, n_out,
                  if ((ne_in + ni_in) > 0 && n_out > (ne_in + ni_in)) "  <- most of this season is OUTSIDE the window (fine for a deliberate part-season run)" else ""))
    # 2026-09-10: the asymmetric case that motivated splitting the counts. The first
    # two-season staging found 13,629 INTERVIEWS for 2023-24 and ZERO effort counts:
    # the interview workbook is the full export while effort_combined was trimmed to
    # 2024-25 during development. A season with interviews but no effort counts has a
    # CPUE and no effort scale, so every BSS effort likelihood for it is empty and the
    # PE has nothing to expand; the season would be imputed, not estimated.
    if (ni_in > 0 && ne_in == 0)
      warning(sprintf(paste0("Season '%s' has %d interviews in the window but ZERO effort counts. ",
                             "Its effort process would be pure imputation. Add the season's rows to ",
                             "effort_combined.xlsx (and wes_commercial_tally.xlsx) before fitting."),
                      sn, ni_in), call. = FALSE)
    if (ne_in > 0 && ni_in == 0)
      warning(sprintf("Season '%s' has effort counts but ZERO interviews in the window; CPUE for it has no data.", sn),
              call. = FALSE)
  }

  n_eff_in <- sum(eff_dates >= ws & eff_dates <= we, na.rm = TRUE)
  n_int_in <- sum(int_dates >= ws & int_dates <= we, na.rm = TRUE)
  if (n_eff_in == 0 && n_int_in == 0)
    stop(sprintf(paste0("validate_season_window(): the estimation window %s to %s contains NO effort ",
                        "counts and NO interviews for season_filter = %s. The usual cause is a stale ",
                        "season_filter after moving est_date_start/est_date_end to a new season; the two ",
                        "must be updated together (see 07_documentation/NEW_SEASON_GUIDE.md)."),
                 ws, we, paste(sprintf("'%s'", seasons_requested), collapse = " + ")), call. = FALSE)

  cs <- suppressWarnings(as.Date(params$census_start_date %||% NA))
  ce <- suppressWarnings(as.Date(params$census_end_date   %||% NA))
  if (!is.na(cs) && !is.na(ce) && (ce < ws || cs > we))
    warning(sprintf(paste0("The commercial/charter census window (%s to %s) lies entirely outside the ",
                           "estimation window; the census component will be empty. census_start_date / ",
                           "census_end_date are per-season settings."), cs, ce), call. = FALSE)
  invisible(TRUE)
}
