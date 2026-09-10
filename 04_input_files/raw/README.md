# 04_input_files/raw

The creel-database workbooks the model inputs are BUILT from. Nothing in the pipeline reads these directly; the builders in `04_input_files/` do (`build_all_inputs.R` runs the set), and the workbooks they write are what the readers use.

## The per-season creel workbooks (since 2026-09-10)

One workbook per fishery season, the creel database exported with every sheet, received 2026-09-10 from Matt George (file names had spaces; renamed with underscores, content untouched):

| Workbook | Season (Sep 16 to Sep 15) | Primary sheets used | Rows (interviews / effort counts / surveys) | Through |
|---|---|---|---|---|
| `2223_rec_crab_harvest_data.xlsx` | 2022-23 (data from 2023-01-03) | `Harvest Creel`, `Effort count`, `crab creel survey data` | 3,104 / 535 / 230 | 2023-09-15 |
| `2324_rec_crab_harvest_data.xlsx` | 2023-24 | the same, plus `location_key`, `gear_key` | 13,629 / 4,116 / 1,369 | 2024-09-15 |
| `2425_rec_crab_harvest_data.xlsx` | 2024-25 | the same, plus `wes commercial tally`, `charter trips` | 12,237 / 3,256 / 805 | 2025-09-13 |
| `2526_rec_crab_harvest_data.xlsx` | 2025-26 | the same, plus `wes commercial tally`, `charter trips` (the survey sheet carries the vessel tallies too) | 9,992 / 4,149 / 673 | 2026-09-08 |

What each sheet is, and which builder reads it:

- **`Harvest Creel`** (one row per interview) -> `build_interview_combined.R`. Headers differ by season ("Number of Gear" / "Number Of Gear"; "Crab Released?" / "Crab released?"), the 2022-23 sheet lacks the tampered flag, the released-crab fields and the vehicle count, and the 2024-25 sheet adds `Total Vehicles`. The builder reads by normalised header name.
- **`Effort count`** (one row per instantaneous count at a site) -> `build_effort_combined.R`. 2022-24: `Truck/Boat Trailer Count`, `Buoy Count`; 2024-25 on: `Boat Trailer Count` and `Vehicle Count` separately. Notes record the Feb 2024 interview-derived totals the builder flags.
- **`crab creel survey data`** (one row per sampler site-visit) -> `build_sampler_shifts.R` (check-in / check-out, the on-site conditions, the sampler's holiday flag) and, for 2025-26, `build_comm_charter_tally.R` (the `Commercial / Private / Charter Boat Tally` columns).
- **`wes commercial tally`** (2024-25 on; a two-row header, the column ORDER differs between seasons) -> `build_comm_charter_tally.R`, by header name.
- **`charter trips`** (2024-25 on; the operators' per-trip roster with interviewed / missed / canceled, plus a hand summary at the right that is not read) -> `build_charter_trips.R`.
- **`location_key`** (creel area -> creel location, I/E abbreviation, north of Point Chehalis) and **`gear_key`** (the 2022-24 gear labels against the 2024-25 vocabulary): encoded in `build_helpers.R` (`area_map`) and `build_interview_combined.R` (`gear_map`).
- Not read: the derived summary and pivot tabs (`boat summary`, `jetty summary`, `thanksgiving summary`, `Christmas New Year's Summary`, `last clam dig recap`, `Beach info`, `summary`, `days sampled`, `effort pivot table`, `interview pivot table`, `labor day weekend`). They are computed from the primary sheets and carry nothing the model needs.

Known blemishes in the sources, all handled by the builders without editing the files: one 2023-24 Dungeness-kept cell holds a formula fragment ("3+N6847:N13622"); one 2023-24 completed-trip value is 24; one 2023-24 red-rock-returned value is "unknown"; a stray fifth column of notes in the 2025-26 harvest sheet; dates and clock times stored as Excel cells in some rows and as text in others; a duplicated survey id (1677) in 2023-24; eight check-outs typed on a 12-hour clock.

## Retired (2026-09-10)

`interviewdata20222026.xlsx` and `surveydata20222026.xlsx` (received 2026-09-09) were column-pasted compilations of the season workbooks. The pasting had left the 2023-24 gear count under a second, differently-capitalised header, and the 2025-26 rows stopped at 2026-08-01. They are in the git history and are no longer read.

## Not here

`ingress_egress.xlsx` is maintained from the I/E database export (`ingressegressdata<YYYY><YYYY>.xlsx`); the copy in the repository (through 2026-08-28) is newer than the last export supplied and was verified identical to it on every shared interval on 2026-09-10. When a newer I/E export arrives, add it here and rebuild that workbook (see `../README.md`). `WBL_boat_counts.xlsx` (OSP) and `fishery_opener_dates.xlsx` have their own sources.

Keep the workbooks as received (the builders never modify them). A new season is: add its workbook here under the same naming pattern, run `Rscript 04_input_files/build_all_inputs.R` from the repository root, read the builders' reports, and commit the workbook and the rebuilt inputs together.
