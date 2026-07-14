# 04_input_files

The raw season data the drivers read. These are the only inputs a local run needs; everything in `05_output/` is regenerated from these input files plus the model and parameter settings. Reads are done with `here("04_input_files", <name>)`, so the files are found regardless of which driver knits or where it sits.

For what the estimate is and how these feed the PE and BSS paths, see the [root README](../README.md).

## Files

| File | Read by | Via | Contents |
|---|---|---|---|
| `effort_combined.csv` | both `01_BSS_models/` drivers | `read_csv(here("04_input_files","effort_combined.csv"))` | Effort counts (census and index) by date, mode, and location. The effort observation series for the BSS and the PE. |
| `interview_combined.csv` | both `01_BSS_models/` drivers | `read_csv(here("04_input_files","interview_combined.csv"))` | Creel interviews: catch, gear, party size, fishing time. Source of CPUE and the gear proportions. |
| `wes_commercial_tally.csv` | both `01_BSS_models/` drivers | `read_csv(here("04_input_files","wes_commercial_tally.csv"))` | Westport commercial/charter vessel tally. Folded into the port total by census expansion (no separate BSS fit for this component). |
| `ingress_egress.xlsx` | pooled driver only | `read_excel(here("04_input_files", params$ie_data_file), sheet = params$ie_sheet)` | Ingress/egress trip-timing observations. Feeds the `L_effective` regression that anchors effective day length (shore) in the pooled model. The gear-resolved driver does not read it. |

The pooled driver names the I/E workbook through parameters (`ie_data_file = "ingress_egress.xlsx"`, `ie_sheet = "data"`), so a different I/E file or sheet can be swapped in without touching code. The three CSVs are referenced by literal name in both drivers.

## Schema quirks (carried from the iForm exports)

These are real properties of the raw exports, not bugs in the pipeline. They are handled in the data-prep chunks, but they matter if you regenerate or edit an input by hand:

- **Interview gear column mapping.** `number_of_gear` is read from column **N**, not column W, in the raw iForm interview export. The export has a duplicate field name and the intended values live in column N. Re-exporting from a different form version can silently shift this.
- **Effort CSV quoting.** Re-export `effort_combined.csv` with `QUOTE_ALL`. The notes field contains commas; without full quoting the column alignment breaks on read.
- **Interview date format.** Dates are M/D/YYYY and are parsed with `col_date(format = "%m/%d/%Y")`. A locale or spreadsheet that rewrites these to ISO or two-digit years will not parse.
- **Boat-type typo.** The iForm export spells the commercial category "Commerical" (one 'm'). The prep code matches it by regex, so do not "correct" the spelling in the raw file without also updating the matcher.

## Regenerating an input

If you rebuild any of these from a fresh iForm pull, keep the column names and the four quirks above intact, then re-run a driver from `01_BSS_models/`. No path changes are needed; the drivers always read from this folder by name.
