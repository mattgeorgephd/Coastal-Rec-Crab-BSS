###############################################################################
# diagnose_fishery_spillover.R  (called by the pooled report; diagnostic only)
#
# Tests whether crab EFFORT (shore gear counts, boat trailer counts) and CPUE
# (interview Dungeness catch per crabber-hour) differ on other-fishery opener days:
# Marine Area 2 salmon / halibut / bottomfish, and coastal razor-clam digs on the
# nearby beaches (Twin Harbors / Copalis / Mocrocks). The question is whether these
# openers spill over onto crabbing enough to deserve their own day category (like the
# crabbing holidays). This is a DIAGNOSTIC: it changes no estimate and adds no term to
# the BSS or PE; it only reports associations for the analyst to judge on the next run.
#
# Two views per response, because the openers are confounded with day-type (they fall
# on weekends) and season (salmon is a summer fishery, digs are fall/winter), which
# already drive crab effort:
#   RAW       mean on opener days vs non-opener days, with a Welch t-test.
#   ADJUSTED  the opener's coefficient in lm(response ~ day_type + month + opener),
#             i.e. the marginal association AFTER removing the weekend and monthly
#             signal the model already captures. Each opener is fit in its own model
#             (adjusted for day-type + month, not for the other openers), since the
#             openers are largely disjoint in time. Salmon is nearly collinear with
#             summer, so read its adjusted estimate with care; the note column flags a
#             term that could not be identified.
#
# Runs on the in-memory dwg (same series as the Section 3 input plots), so it needs no
# extra data prep and can run before the multi-hour fits.
###############################################################################

diagnose_fishery_spillover <- function(dwg, params) {
  events <- prep_fishery_events(params)

  flag_cols   <- c("ma2_salmon_open", "ma2_halibut_open", "ma2_bottomfish_open", "razor_nearby_dig")
  flag_labels <- c(ma2_salmon_open     = "MA2 salmon open",
                   ma2_halibut_open    = "MA2 halibut open",
                   ma2_bottomfish_open = "MA2 bottomfish open",
                   razor_nearby_dig    = "Razor dig (nearby beaches)")
  all_flags   <- c("ma2_salmon_open", "ma2_halibut_open", "ma2_bottomfish_open",
                   "razor_any_dig", "razor_nearby_dig")

  # Attach opener flags + day_type + month to a per-date frame. The razor calendar is
  # sparse (dig days only), so a date not present means "no dig" -> FALSE, not NA.
  attach_ctx <- function(df) {
    df <- df |>
      dplyr::left_join(events$ma2,   by = "event_date") |>
      dplyr::left_join(events$razor, by = "event_date")
    for (fc in all_flags) df[[fc]][is.na(df[[fc]])] <- FALSE
    df |>
      dplyr::mutate(day_type = classify_day_type(event_date, params),
                    month    = factor(format(event_date, "%Y-%m")))
  }

  # --- Response frames -----------------------------------------------------
  shore_daily <- dwg$shore_effort |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(value = sum(count_quantity), .groups = "drop") |>
    attach_ctx() |> dplyr::mutate(series = "Shore gear (effort)")
  boat_daily <- dwg$boat_effort |>
    dplyr::group_by(event_date) |>
    dplyr::summarise(value = sum(count_quantity), .groups = "drop") |>
    attach_ctx() |> dplyr::mutate(series = "Boat trailers (effort)")

  cpue_all <- dwg$interview |>
    dplyr::filter(fishing_time_total > 0, dungeness_kept >= 0,
                  population %in% c("shore", "private_boat")) |>
    dplyr::mutate(value = dungeness_kept / fishing_time_total) |>
    dplyr::select(event_date, population, value) |>
    attach_ctx()
  cpue_shore <- cpue_all |> dplyr::filter(population == "shore")        |> dplyr::mutate(series = "Shore CPUE")
  cpue_boat  <- cpue_all |> dplyr::filter(population == "private_boat") |> dplyr::mutate(series = "Boat CPUE")

  responses <- list(shore_daily, boat_daily, cpue_shore, cpue_boat)

  # --- RAW: opener vs non-opener means + Welch t-test ----------------------
  raw_one <- function(df) {
    lab <- df$series[1]
    dplyr::bind_rows(lapply(flag_cols, function(fc) {
      x <- df$value; g <- df[[fc]] & !is.na(df$value)
      h <- !df[[fc]] & !is.na(df$value)
      me <- mean(x[g], na.rm = TRUE); mn <- mean(x[h], na.rm = TRUE)
      pv <- tryCatch(stats::t.test(x[g], x[h])$p.value, error = function(e) NA_real_)
      tibble::tibble(
        Series = lab, Opener = unname(flag_labels[fc]),
        n_opener = sum(g), mean_opener = me,
        n_other = sum(h), mean_other = mn,
        raw_diff = me - mn,
        raw_pct  = if (is.finite(mn) && mn != 0) 100 * (me - mn) / mn else NA_real_,
        t_p = pv)
    }))
  }

  # --- ADJUSTED: opener coefficient in lm(value ~ day_type + month + opener) -
  adj_one <- function(df) {
    lab <- df$series[1]
    dplyr::bind_rows(lapply(flag_cols, function(fc) {
      d <- df; d$resp <- d$value; d$opener_flag <- as.integer(d[[fc]])
      base <- tibble::tibble(Series = lab, Opener = unname(flag_labels[fc]),
                             adj_estimate = NA_real_, se = NA_real_,
                             ci_lo = NA_real_, ci_hi = NA_real_, adj_p = NA_real_, note = "")
      ok <- !is.na(d$resp)
      if (length(unique(d$opener_flag[ok])) < 2) { base$note <- "opener constant in sample"; return(base) }
      # Plain (non-dotted) column names: leading-dot names can interact badly with
      # formula/tidy-eval parsing (this project hit that once with a `.fit` column).
      fit <- tryCatch(stats::lm(resp ~ day_type + month + opener_flag, data = d), error = function(e) NULL)
      if (is.null(fit)) { base$note <- "model failed"; return(base) }
      co <- stats::coef(summary(fit))
      if (!"opener_flag" %in% rownames(co) || is.na(co["opener_flag", "Estimate"])) {
        base$note <- "not identified (collinear with day-type/month)"; return(base)
      }
      est <- co["opener_flag", "Estimate"]; se <- co["opener_flag", "Std. Error"]
      base$adj_estimate <- est; base$se <- se
      base$ci_lo <- est - 1.96 * se; base$ci_hi <- est + 1.96 * se
      base$adj_p <- co["opener_flag", "Pr(>|t|)"]
      base
    }))
  }

  raw      <- dplyr::bind_rows(lapply(responses, raw_one))
  adjusted <- dplyr::bind_rows(lapply(responses, adj_one))

  # --- Identifiability / overlap notes (surface the confounds) -------------
  overlap_nearby_any <- shore_daily |>
    dplyr::summarise(n_nearby = sum(razor_nearby_dig),
                     n_any    = sum(razor_any_dig),
                     n_diff   = sum(razor_nearby_dig != razor_any_dig)) |> as.list()
  n_open <- sapply(flag_cols, function(fc) sum(shore_daily[[fc]]))
  notes <- c(
    sprintf("Nearby-dig vs any-dig agree on %d of %d sampled shore days (differ on %d); when 0, the nearby flag is not separable from a generic dig day this season.",
            nrow(shore_daily) - overlap_nearby_any$n_diff, nrow(shore_daily), overlap_nearby_any$n_diff),
    sprintf("Opener OPEN-day counts in the sampled shore-effort series: salmon %d, halibut %d, bottomfish %d, nearby dig %d (of %d sampled days).",
            n_open["ma2_salmon_open"], n_open["ma2_halibut_open"], n_open["ma2_bottomfish_open"],
            n_open["razor_nearby_dig"], nrow(shore_daily)),
    "Bottomfish is open most of the year (near-constant), and salmon is open only in summer (near-collinear with month); read those adjusted estimates with care."
  )

  list(
    raw        = raw,
    adjusted   = adjusted,
    effort_long = dplyr::bind_rows(shore_daily, boat_daily),
    cpue_long   = dplyr::bind_rows(cpue_shore, cpue_boat),
    flag_cols  = flag_cols,
    flag_labels = flag_labels,
    notes      = notes
  )
}
