# 04_input_files

The raw season data the drivers read. These are the only inputs a local run needs; everything in `05_output/` is regenerated from these input files plus the model and parameter settings. Reads are done with `here("04_input_files", <name>)`, so the files are found regardless of which driver knits or where it sits.

For what the estimate is and how these feed the PE and BSS paths, see the [root README](../README.md).

## Canonical input set (2026-07-16): six workbooks

All inputs are now `.xlsx` workbooks with a single `data` sheet (converted from CSV). Every file carries a fishery-season column, dates are stored as ISO `yyyy-mm-dd` text, and each file has been trimmed to the columns the pipeline actually uses.

| File | Role | Read by | Reader | Key columns |
|---|---|---|---|---|
| `effort_combined.xlsx` | model | both drivers | `fetch_crab_data` / `_v2` | `season, date, survey_id, creel_area, count_time, total_gear_count, boat_trailer_count` |
| `interview_combined.xlsx` | model | both drivers | `fetch_crab_data` / `_v2` | 16 cols incl. `survey_id, date, crabbing_mode, boat_type, crabbers, gear_type, number_of_gear, dungeness_kept, red_rock_kept, hours_fished, crabber_hours, gear_hours, completed_trip, season, creel_location` |
| `wes_commercial_tally.xlsx` | model | both drivers | `estimate_comm_charter` (+ a gear-report plot) | `date, season, private_tally, commercial_tally, charter_tally` |
| `ingress_egress.xlsx` | model | both drivers | `fetch_ie_data` | `location_name, date, season, day_type, crabbers_on, crabbers_off, crabber_flow, boats_in, boats_out, boat_flow` |
| `crabbing_holidays.xlsx` | day-typing | both drivers (+ weather module) | `read_crabbing_holidays` | `season, date, holiday_name` |
| `fishery_opener_dates.xlsx` | diagnostic | pooled only, via `prep_fishery_events` | `readxl::read_excel` | `date, season, ma2_bottomfish, ma2_halibut, ma2_salmon, razor_long_beach, razor_twin_harbors, razor_copalis, razor_mocrocks, razor_kalaloch, razor_any` |

The former per-fishery opener workbooks (`MA2-fishing-dates-2023-2026.xlsx`, `razor-clam-dig-dates-2021-2025.xlsx`) are retired; their content is represented in `fishery_opener_dates.xlsx`.

## Season columns (fishery season)

Every workbook carries a `season` column holding the **fishery season** label (`2024-25` = the window Sep 16, 2024 through Sep 15, 2025), not a calendar season. For the multi-year files (`ingress_egress`, `fishery_opener_dates`) the label is derived per row from the date, using a Sep 16 boundary (a date on or after Sep 16 of year Y is season `Y-(Y+1)`). This replaces the old `ingress_egress` `season` column, which held a calendar season (`Fall`/`Winter`/`Spring`/`Summer`) that no code used.

Which readers filter on it:

- `effort_combined.xlsx` and `interview_combined.xlsx` filter to `season == run_config$season_filter` (existing behavior).
- `wes_commercial_tally.xlsx` carries the column but is scoped by the census date window in `estimate_comm_charter` (unchanged); the column is provenance.
- `fishery_opener_dates.xlsx` carries the column; the spillover diagnostic joins by date, so all rows are kept (a lookup calendar).
- `ingress_egress.xlsx` carries the column but is **not** filtered by default: the `L_effective` day-length regression intentionally pools all seasons of I/E history. Set `run_config$ie_filter_by_season = TRUE` to restrict it to the current season.

## Columns: trimmed to modeling-essential

The files carry only the columns the pipeline uses (a whole-repo reference scan drove the cut). Notable drops: the effort weather/vehicle/buoy/jetty counts and `notes`; the interview locality/`*_returned`/`notes`/`total_vehicles`/`crabbing_holiday`/`interview_time` fields; the tally `*_interviewed` columns and `notes`; and the `ingress_egress` weather/metadata columns (the weather module fetches tide/weather from external APIs, not from this file). The `ingress_egress` file still carries `crabbers_on`, `crabbers_off`, and `boats_out` (the raw ingress/egress survey tallies), though only `crabber_flow` (shore) and `boats_in` (boat) reach a model.

## Input options in run_config

Input selection is surfaced in `run_config.R` so a different port, site set, or workbook can be run without editing the readers:

- Filenames: `effort_file`, `interview_file`, `tally_file`, `input_sheet`; `ie_data_file` / `ie_sheet`; `crabbing_holidays_file` / `crabbing_holidays_sheet`; `fishery_opener_dates_file` / `fishery_opener_sheet`.
- Row/site selection: `gh_creel_location` (interview), `gh_effort_areas` (effort site whitelist), `shore_dock_float20` / `shore_dock_float17` (the paired gear-count floats), `boat_launch_areas` (trailer-count sites), and `ie_shore_location` / `ie_boat_location` (the I/E `location_name` series).

Each reader uses the config value with the historical Grays Harbor / Westport value as a fallback, so a run works even if a key is not set.

## Notes on the workbook format

- **Dates are ISO text.** Each `date` column stores `yyyy-mm-dd` as text (not an Excel date cell) so `as.Date()` parses it identically in any timezone. Keep that format on re-export. This also resolves a prior defect where `effort_combined.csv` and `fishery_opener_dates.csv` had been re-exported as `M/D/YYYY`, which their `as.Date(date)` readers would have read as `NA`.
- **Sheet name.** Every workbook uses a single sheet named `data` (override per file via the sheet params above).
- **Interview gear column.** `number_of_gear` is the intended gear count (historically column N in the raw iForm export); keep the header name intact on regeneration.
- **Boat-type typo.** The source spells the commercial category "Commerical" (one 'm'); the prep matches it by regex, so do not "correct" it without updating the matcher.

## Regenerating an input

Rebuild a workbook keeping the exact column names, the `data` sheet, ISO-text dates, and a `season` column matching `run_config$season_filter` for the season being run, then re-run a driver from `01_BSS_models/`. Adding extra columns is harmless (readers reference columns by name). The full column-by-column rationale is in the input-file audit saved with the project docs.
