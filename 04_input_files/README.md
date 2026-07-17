# 04_input_files

The raw season data the drivers read. These are the only inputs a local run needs; everything in `05_output/` is regenerated from these input files plus the model and parameter settings. Reads are done with `here("04_input_files", <name>)`, so the files are found regardless of which driver knits or where it sits.

For what the estimate is and how these feed the PE and BSS paths, see the [root README](../README.md).

## Canonical input set (2026-07-16)

Six files, in three roles. Two changes landed on 2026-07-16: the former per-fishery opener workbooks (`MA2-fishing-dates-2023-2026.xlsx`, `razor-clam-dig-dates-2021-2025.xlsx`) are **retired** (`prep_fishery_events` now reads only the consolidated `fishery_opener_dates.csv` and stops if it is absent, with no silent fallback), and the crabbing holidays moved out of `run_config.R` into the editable workbook `crabbing_holidays.xlsx`.

| File | Role | Read by | Via | Contents |
|---|---|---|---|---|
| `effort_combined.csv` | model | both `01_BSS_models/` drivers | `fetch_crab_data` / `fetch_crab_data_v2` | Effort counts (dock gear counts, boat-trailer counts) by date and location. The effort observation series for the BSS and the PE. |
| `interview_combined.csv` | model | both drivers | `fetch_crab_data` / `fetch_crab_data_v2` | Creel interviews: catch, gear, party size, fishing time. Source of CPUE and the gear proportions. |
| `wes_commercial_tally.csv` | model | both drivers | `estimate_comm_charter` (+ a gear-report vessel-tally plot) | Westport commercial/charter vessel tally. Folded into the port total by census expansion (no separate BSS fit for this component). |
| `ingress_egress.xlsx` | model | **both drivers** (`fetch_ie_data`) | `read_excel(here("04_input_files", params$ie_data_file), sheet = params$ie_sheet)` | Ingress/egress trip-timing observations. Feeds the `L_effective` regression that anchors shore effective day length. **Both** the pooled and gear-resolved drivers call `fetch_ie_data` (an earlier version of this note said pooled-only; that was incorrect). |
| `crabbing_holidays.xlsx` | day-typing | both drivers (and the weather module) | `read_crabbing_holidays` | High-effort non-weekend days treated as weekend for day-typing. Columns `season, date, holiday_name`, filtered to `season_filter`. Replaces the hardcoded `crabbing_holiday_dates` vector formerly in `run_config.R` (and a duplicated copy in the weather module). |
| `fishery_opener_dates.csv` | diagnostic | pooled driver only, via `prep_fishery_events` | `read.csv(here("04_input_files","fishery_opener_dates.csv"))` | Consolidated daily MA2 finfish (salmon / halibut / bottomfish) and coastal razor-clam-dig OPEN/CLOSED calendar, one row per day. The **only** source for the fishery-spillover diagnostic (pooled report Section 3.5) and the optional razor-dig effort term. Required whenever `run_fishery_spillover_diag = TRUE`. |

**Retired** (no longer read by any driver; retained for provenance unless deleted): `MA2-fishing-dates-2023-2026.xlsx`, `razor-clam-dig-dates-2021-2025.xlsx`. Their content is fully represented in `fishery_opener_dates.csv`.

The pooled driver names the I/E workbook and the opener/holiday files through parameters (`ie_data_file`, `ie_sheet`, `fishery_opener_dates_file`, `crabbing_holidays_file`, `crabbing_holidays_sheet`), so a different file or sheet can be swapped in without touching code. The three core CSVs (effort, interview, tally) are referenced by literal name in both drivers.

## Schema quirks (carried from the iForm exports)

These are real properties of the raw exports, not bugs in the pipeline. They are handled in the data-prep chunks, but they matter if you regenerate or edit an input by hand:

- **Interview gear column mapping.** `number_of_gear` is read from column **N**, not column W, in the raw iForm interview export. The export has a duplicate field name and the intended values live in column N. Re-exporting from a different form version can silently shift this.
- **Effort CSV quoting.** Re-export `effort_combined.csv` with `QUOTE_ALL`. The notes field contains commas; without full quoting the column alignment breaks on read.
- **Interview date format.** Dates are M/D/YYYY and are parsed with `col_date(format = "%m/%d/%Y")`. A locale or spreadsheet that rewrites these to ISO or two-digit years will not parse.
- **Boat-type typo.** The iForm export spells the commercial category "Commerical" (one 'm'). The prep code matches it by regex, so do not "correct" the spelling in the raw file without also updating the matcher.
- **Holiday workbook dates.** `crabbing_holidays.xlsx` stores `date` as ISO `yyyy-mm-dd` text (not Excel date cells) to avoid a timezone off-by-one on read; keep that format, and keep a `season` value that matches `season_filter` (for example, `2024-25`). The reader stops if no row matches the season, so a new season needs its own rows before that season is run.

## Regenerating an input

If you rebuild any of these from a fresh iForm pull, keep the column names and the quirks above intact, then re-run a driver from `01_BSS_models/`. No path changes are needed; the drivers always read from this folder by name. For which columns each file actually needs (several are unused and can be trimmed at the export step), see the input-file audit saved with the project docs.
