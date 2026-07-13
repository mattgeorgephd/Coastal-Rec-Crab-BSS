###############################################################################
# prep_fishery_events.R  (shared; read by the fishery-opener diagnostic and the
# razor_dig effort term)
#
# Returns a per-date table of other-fishery OPEN/CLOSED flags so the spillover
# diagnostic (diagnose_fishery_spillover.R) and the razor_dig shore-effort term
# (item 1) can use them.
#
# PRIMARY SOURCE (2026-07-13): one consolidated daily calendar CSV,
# 04_input_files/fishery_opener_dates.csv, with columns
#   date, ma2_bottomfish, ma2_halibut, ma2_salmon,
#   razor_long_beach, razor_twin_harbors, razor_copalis, razor_mocrocks, razor_kalaloch, razor_any
# (each OPEN/CLOSED, one row per day). This replaces the two source workbooks
# (MA2-fishing-dates*.xlsx, razor-clam-dig-dates*.xlsx), which are kept only as a
# fallback when the CSV is absent. Override the filename with params$fishery_opener_dates_file.
#
# razor_nearby_dig = any of params$razor_nearby_beaches open (default Twin Harbors +
# Copalis + Mocrocks, the beaches nearest Grays Harbor). In the 2024-25 data Twin Harbors
# is open on every dig day, so nearby == any that season; the code is general.
###############################################################################

prep_fishery_events <- function(params) {
  is_open   <- function(x) toupper(trimws(as.character(x))) == "OPEN"
  beach_col <- function(nm) paste0("razor_", gsub("[^a-z0-9]+", "_", tolower(trimws(nm))))
  nearby_names <- params$razor_nearby_beaches %||% c("Twin Harbors", "Copalis", "Mocrocks")

  csv_path <- here::here("04_input_files",
                         params$fishery_opener_dates_file %||% "fishery_opener_dates.csv")

  if (file.exists(csv_path)) {
    d <- utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
    names(d) <- tolower(trimws(names(d)))

    ma2_flags <- tibble::tibble(
      event_date          = as.Date(d$date),
      ma2_salmon_open     = is_open(d$ma2_salmon),
      ma2_halibut_open    = is_open(d$ma2_halibut),
      ma2_bottomfish_open = is_open(d$ma2_bottomfish)
    ) |>
      dplyr::filter(!is.na(event_date)) |>
      dplyr::distinct(event_date, .keep_all = TRUE)

    beach_all   <- setdiff(grep("^razor_", names(d), value = TRUE), "razor_any")
    nearby_cols <- intersect(beach_col(nearby_names), beach_all)
    open_mat    <- if (length(beach_all))
                     do.call(cbind, lapply(beach_all, function(cc) is_open(d[[cc]])))
                   else matrix(FALSE, nrow(d), 0)
    if (length(beach_all)) colnames(open_mat) <- beach_all

    razor_flags <- tibble::tibble(
      event_date       = as.Date(d$date),
      razor_any_dig    = if ("razor_any" %in% names(d)) is_open(d$razor_any)
                         else rowSums(open_mat, na.rm = TRUE) > 0,
      razor_nearby_dig = if (length(nearby_cols)) rowSums(open_mat[, nearby_cols, drop = FALSE], na.rm = TRUE) > 0
                         else rep(FALSE, nrow(d))
    ) |>
      dplyr::filter(!is.na(event_date)) |>
      dplyr::group_by(event_date) |>
      dplyr::summarise(razor_any_dig    = any(razor_any_dig),
                       razor_nearby_dig = any(razor_nearby_dig), .groups = "drop")

    return(list(ma2 = ma2_flags, razor = razor_flags,
                beach_cols = beach_all, nearby_cols = nearby_cols, source = "csv"))
  }

  # ---- Fallback: the two source workbooks (pre-2026-07-13 consolidation) ------------
  ma2_path   <- here::here("04_input_files",
                           params$ma2_dates_file   %||% "MA2-fishing-dates-2023-2026.xlsx")
  razor_path <- here::here("04_input_files",
                           params$razor_dates_file %||% "razor-clam-dig-dates-2021-2025.xlsx")

  ma2 <- readxl::read_excel(ma2_path, sheet = "data")
  names(ma2) <- toupper(trimws(names(ma2)))
  ma2_flags <- tibble::tibble(
    event_date          = as.Date(ma2$DATE),
    ma2_salmon_open     = is_open(ma2$SALMON),
    ma2_halibut_open    = is_open(ma2$HALIBUT),
    ma2_bottomfish_open = is_open(ma2$BOTTOMFISH)
  ) |>
    dplyr::filter(!is.na(event_date)) |>
    dplyr::distinct(event_date, .keep_all = TRUE)

  razor <- readxl::read_excel(razor_path, sheet = "data")
  names(razor) <- toupper(trimws(names(razor)))
  beach_all    <- intersect(c("LONG BEACH", "TWIN HARBORS", "COPALIS", "MOCROCKS", "KALALOCH"), names(razor))
  nearby_up    <- toupper(trimws(nearby_names))
  nearby_cols  <- intersect(nearby_up, beach_all)
  open_mat <- do.call(cbind, lapply(beach_all, function(cc) is_open(razor[[cc]])))
  colnames(open_mat) <- beach_all
  razor_flags <- tibble::tibble(
    event_date       = as.Date(razor$DATE),
    razor_any_dig    = rowSums(open_mat, na.rm = TRUE) > 0,
    razor_nearby_dig = rowSums(open_mat[, nearby_cols, drop = FALSE], na.rm = TRUE) > 0
  ) |>
    dplyr::filter(!is.na(event_date)) |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(razor_any_dig    = any(razor_any_dig),
                     razor_nearby_dig = any(razor_nearby_dig), .groups = "drop")

  list(ma2 = ma2_flags, razor = razor_flags,
       beach_cols = beach_all, nearby_cols = nearby_cols, source = "xlsx")
}
