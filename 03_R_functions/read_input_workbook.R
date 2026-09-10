# -----------------------------------------------------------------------------
# Part of Coastal-Rec-Crab-BSS: recreational Dungeness crab creel estimation
# for Grays Harbor / Westport (WDFW).
# Copyright (C) 2024-2026 Washington Department of Fish and Wildlife.
#
# Adapted from CreelEstimates, the WDFW freshwater creel estimation framework:
#   https://github.com/dfw-wa/CreelEstimates   (licensed GPL-3.0).
# Substantial portions of the methodology, structure, and R/Stan code originate
# in CreelEstimates and remain (C) their authors under GPL-3.0; changes for
# recreational crab are by WDFW.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License, version 3, as published by the Free
# Software Foundation. It is distributed WITHOUT ANY WARRANTY; without even the
# implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program (see the LICENSE file);
# if not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
###############################################################################
# read_input_workbook.R  (2026-09-10)
#
# The one way the pipeline reads a workbook in 04_input_files/.
#
# WHY. readxl guesses each column's type from its first `guess_max` rows (default
# 1,000). Every model workbook is now multi-season and sorted by date, and several
# columns are blank for a whole early season (gear_tampered, total_vehicles and the
# released-crab fields do not exist in the 2022-23 sheet; vehicle_count and
# weather_location start in 2024-25). A column whose first 1,000 rows are blank is
# guessed LOGICAL, and the numbers below the guess window are then coerced: 0/1 survive
# as FALSE/TRUE (the 2026-09-09 gear_tampered column was read that way and happened to
# work), but a count of 2, 3, 4 becomes TRUE and text becomes NA. This wrapper reads with
# guess_max = 100,000 (every workbook is far smaller), so the type is guessed from the
# whole column. Nothing else changes: the same readxl call, the same sheet.
#
# read_input_workbook(file, sheet = "data", ...): `file` is a name in 04_input_files/ or
# an absolute path; extra arguments go to readxl::read_excel.
###############################################################################

read_input_workbook <- function(file, sheet = "data", ..., guess_max = 100000L) {
  path <- if (file.exists(file)) file else here::here("04_input_files", file)
  if (!file.exists(path)) stop("input workbook not found: ", path, call. = FALSE)
  readxl::read_excel(path, sheet = sheet, guess_max = guess_max, ...)
}
