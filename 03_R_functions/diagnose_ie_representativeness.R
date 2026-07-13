###############################################################################
# diagnose_ie_representativeness.R  (pooled report; diagnostic only)
#
# Item 4 / GR-9: is the shore all-gear sigma_IE (~1.07 in the 2026-07-12 run, vs ~0.60
# in the gear-resolved fit) a symptom of the I/E observation days being unrepresentative
# (measured on peak-effort days, so the I/E-derived effective day length disagrees with
# the effort counts), or just sparse-data scatter across a handful of I/E days?
#
# This compares the shore daily GEAR COUNTS on I/E observation days against all other
# sampled days. A mean percentile rank near 0.5 (with a non-significant rank test) means
# the I/E days sit in the middle of the effort distribution, i.e. they are representative,
# so a large sigma_IE is sparse-data noise, not peak-day sampling bias. A percentile well
# above 0.6 would mean the I/E days skew high and the I/E stream should be treated with
# caution (or reweighted). DIAGNOSTIC ONLY: it changes no estimate.
###############################################################################

diagnose_ie_representativeness <- function(dwg, ie_data, params) {
  if (is.null(ie_data) || nrow(ie_data) == 0) return(NULL)
  ie_days <- ie_data |>
    dplyr::filter(population == "shore") |>
    dplyr::distinct(event_date) |>
    dplyr::pull(event_date)
  if (length(ie_days) == 0) return(NULL)

  daily <- dwg$shore_effort |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(gear = sum(count_quantity), .groups = "drop") |>
    dplyr::mutate(is_ie = event_date %in% ie_days)

  g_ie  <- daily$gear[daily$is_ie]
  g_non <- daily$gear[!daily$is_ie]
  if (length(g_ie) == 0 || length(g_non) == 0) return(NULL)

  # Mean percentile rank of the I/E days within the full gear-count distribution
  # (0.5 = representative; > 0.6 = skewed toward higher-effort days).
  pctile <- mean(vapply(g_ie, function(v) mean(daily$gear < v, na.rm = TRUE), numeric(1)))
  pval   <- tryCatch(stats::wilcox.test(g_ie, g_non)$p.value, error = function(e) NA_real_)

  verdict <- if (!is.finite(pctile)) "n/a"
    else if (pctile > 0.6)
      sprintf("I/E days skew toward HIGHER-effort days (mean percentile %.2f, rank-test p = %.2f); the I/E-derived L may overstate effective day length and inflate sigma_IE. Treat the shore I/E stream with caution.", pctile, pval)
    else if (pctile < 0.4)
      sprintf("I/E days skew toward LOWER-effort days (mean percentile %.2f, rank-test p = %.2f).", pctile, pval)
    else
      sprintf("I/E days are representative in effort-count terms (mean percentile %.2f, rank-test p = %.2f); a large shore sigma_IE is sparse-data scatter across few I/E days, not peak-day sampling bias.", pctile, pval)

  list(
    table = tibble::tibble(
      group       = c("I/E days", "non-I/E days"),
      n           = c(length(g_ie), length(g_non)),
      mean_gear   = c(mean(g_ie), mean(g_non)),
      median_gear = c(stats::median(g_ie), stats::median(g_non))),
    pctile = pctile, pval = pval, n_ie = length(g_ie), verdict = verdict)
}
