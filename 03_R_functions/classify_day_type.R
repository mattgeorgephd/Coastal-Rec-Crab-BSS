###############################################################################
# classify_day_type.R
#
# Standalone day-type classifier: assign holiday / weekend / weekday to ANY date,
# independent of the estimation calendar (used for diagnostic plots that may show
# data outside the estimation window). Extracted from the pooled driver. weekends
# and holidays are read from params (params$days_wkend, params$crabbing_holiday_dates).
###############################################################################

classify_day_type <- function(dates, params) {
  weekends <- params$days_wkend
  holidays <- params$crabbing_holiday_dates
  case_when(
    dates %in% holidays ~ "holiday",
    weekdays(dates) %in% weekends ~ "weekend",
    TRUE ~ "weekday"
  )
}
