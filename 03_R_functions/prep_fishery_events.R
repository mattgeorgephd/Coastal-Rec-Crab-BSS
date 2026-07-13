###############################################################################
# prep_fishery_events.R  (shared; read by the pooled fishery-opener diagnostic)
#
# Reads the two other-fishery opener calendars in 04_input_files/ and returns a
# per-date table of OPEN/CLOSED flags, so the spillover diagnostic
# (03_R_functions/diagnose_fishery_spillover.R) can test whether crab effort/CPUE
# differs on those dates. Diagnostic only; it does NOT change the estimate.
#
#   MA2-fishing-dates:  a FULL daily calendar with OPEN/CLOSED for Marine Area 2
#                       BOTTOMFISH, HALIBUT, SALMON (Westport finfish seasons).
#   razor-clam-dig:     a SPARSE list of dig days (rows exist only for open digs),
#                       with OPEN/CLOSED per beach (Long Beach, Twin Harbors, Copalis,
#                       Mocrocks, Kalaloch) plus an ANY column.
#
# Filenames come from params (ma2_dates_file / razor_dates_file); the "nearby" beach
# set is params$razor_nearby_beaches (default Twin Harbors + Copalis + Mocrocks, the
# beaches closest to Grays Harbor). NOTE on identifiability: in the 2024-25 data Twin
# Harbors and Long Beach are open on EVERY listed dig day, so razor_nearby_dig equals
# razor_any_dig that season; the code is general (a future season where Twin Harbors is
# not always open would separate them), and the diagnostic reports the overlap so the
# collinearity is visible rather than hidden.
###############################################################################

prep_fishery_events <- function(params) {
  is_open <- function(x) toupper(trimws(as.character(x))) == "OPEN"

  ma2_path   <- here::here("04_input_files",
                           params$ma2_dates_file   %||% "MA2-fishing-dates-2023-2026.xlsx")
  razor_path <- here::here("04_input_files",
                           params$razor_dates_file %||% "razor-clam-dig-dates-2021-2025.xlsx")

  # --- MA2 finfish (full daily calendar) -----------------------------------
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

  # --- Razor-clam digs (sparse; only dig days present) ---------------------
  razor <- readxl::read_excel(razor_path, sheet = "data")
  names(razor) <- toupper(trimws(names(razor)))
  beach_all    <- c("LONG BEACH", "TWIN HARBORS", "COPALIS", "MOCROCKS", "KALALOCH")
  beach_cols   <- intersect(beach_all, names(razor))
  nearby_names <- toupper(trimws(params$razor_nearby_beaches %||%
                                 c("TWIN HARBORS", "COPALIS", "MOCROCKS")))
  nearby_cols  <- intersect(nearby_names, beach_cols)

  # Per-beach OPEN logical matrix (rows = razor rows, cols = beaches present).
  open_mat <- do.call(cbind, lapply(beach_cols, function(cc) is_open(razor[[cc]])))
  colnames(open_mat) <- beach_cols

  razor_flags <- tibble::tibble(
    event_date       = as.Date(razor$DATE),
    razor_any_dig    = rowSums(open_mat, na.rm = TRUE) > 0,
    razor_nearby_dig = rowSums(open_mat[, nearby_cols, drop = FALSE], na.rm = TRUE) > 0
  ) |>
    dplyr::filter(!is.na(event_date)) |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(razor_any_dig    = any(razor_any_dig),
                     razor_nearby_dig = any(razor_nearby_dig), .groups = "drop")

  list(
    ma2         = ma2_flags,
    razor       = razor_flags,
    beach_cols  = beach_cols,
    nearby_cols = nearby_cols
  )
}
