# 04_input_files/raw

The creel-database exports the model workbooks are BUILT from. Nothing in the pipeline reads these directly; the two builders do, and the workbooks they write are what `fetch_crab_data()` and `fetch_sampler_shifts()` read.

| Export | Received | Builder | Writes |
|---|---|---|---|
| `interviewdata20222026.xlsx` (sheet `data`, 37,460 rows, 44 columns, seasons 2022-23 to 2025-26 through 2026-08-01) | 2026-09-09 from Matt George | `../build_interview_combined.R` | `../interview_combined.xlsx` |
| `surveydata20222026.xlsx` (`Sheet1`, 3,077 surveys, check-in / check-out in 24-hour clock text) | 2026-09-09 from Matt George | `../build_sampler_shifts.R` | `../sampler_shifts.xlsx` |

Keep the exports as received (the builders never modify them). When a new export arrives, add it here under its own name, point the builder at it (first argument), and commit both the export and the rebuilt workbook together, so the workbook's provenance is always one script run from a file in this folder.
