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
# desk_a42_render_check_2026-11-25.R
#
# DESK, SECONDS: reads a finished render against the A42 rule (CHANGE_REGISTER A42, committed
# with the switch in ebc4038 before the render). The rule is in the register and in run_config.R;
# this script only reads it, clause by clause, against B70:
#   (1) every fit on both tracks passes the gate and no chain is stuck (convergence_report.csv)
#   (2) the six fits the switch cannot reach reproduce B70's saved draws byte for byte
#   (3) the two boat all-gear fits' structural_params_*.csv equal B70's, row for row
#   (4) the pooled boat all-gear expected catch equals B86's exact replay, 34,361.4, within 0.5 crab
#   (5) the cross-check is within cross_check_tolerance (its own verdict)
# and reports the port total, the components against B70 and the boat catch by month.
#
# HOW TO RUN   RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_a42_render_check_2026-11-25.R <render date folder>
# WRITES       05_output/a42_render_check_2026-11-25.csv, _components.csv, _months.csv
###############################################################################
.args <- commandArgs(trailingOnly = TRUE)
if (!length(.args)) stop("give the render's dated folder, e.g. 05_output/20261125")
R <- .args[1]
B70 <- c(pooled = "05_output/20261019/pooled-CPUE-2024-25", gear = "05_output/20261019/gear-type-CPUE-model-2024-25")
NEW <- c(pooled = file.path(R, "pooled-CPUE-2024-25"), gear = file.path(R, "gear-type-CPUE-model-2024-25"))
stopifnot(all(dir.exists(NEW)))
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
UNREACHED <- c("shore_all_gear", "shore_ring_net_only", "private_boat_ring_net_only")
rows <- list(); add <- function(clause, what, observed, pass) rows[[length(rows) + 1L]] <<- data.frame(clause = clause, what = what, observed = observed, pass = pass, stringsAsFactors = FALSE)

for (tr in names(NEW)) {
  cr <- rd(file.path(NEW[[tr]], "convergence_report.csv"))
  fitted <- cr[!is.na(cr$ar_resolution), , drop = FALSE]
  gate <- as.logical(fitted$pass_convergence) & fitted$method_selected == "BSS"
  add("1", sprintf("%s: every fit passes the gate", tr), sprintf("%d of %d fits pass", sum(gate, na.rm = TRUE), nrow(fitted)), nrow(fitted) == 4L && all(gate %in% TRUE))
  for (f in UNREACHED) {
    a <- file.path(B70[[tr]], sprintf("bss_draws_summed_%s_Dungeness_Kept.csv", f)); b <- file.path(NEW[[tr]], basename(a))
    same <- file.exists(b) && identical(unname(tools::md5sum(a)), unname(tools::md5sum(b)))
    add("2", sprintf("%s %s: saved draws byte-identical to B70's", tr, f), if (same) "identical" else "DIFFER", same)
  }
  a <- rd(file.path(B70[[tr]], "structural_params_private_boat_all_gear_Dungeness_Kept.csv")); b <- rd(file.path(NEW[[tr]], "structural_params_private_boat_all_gear_Dungeness_Kept.csv"))
  m <- merge(a, b, by = "parameter"); num <- intersect(c("estimate", "median"), names(a))
  gap <- max(unlist(lapply(num, function(k) abs(m[[paste0(k, ".x")]] - m[[paste0(k, ".y")]]))), na.rm = TRUE)
  add("3", sprintf("%s boat all-gear: structural parameters equal B70's", tr), sprintf("%d of %d rows; max gap %.2e", nrow(m), nrow(a), gap), nrow(m) == nrow(a) && nrow(a) == nrow(b) && gap == 0)
}
pv <- rd(file.path(NEW[["pooled"]], "pe_vs_bss_comparison.csv"))
boat <- pv$BSS_catch[pv$component == "private_boat (All gear)"]
add("4", "pooled boat all-gear expected catch against B86's replay (34,361.4)", sprintf("%.1f", boat), isTRUE(abs(boat - 34361.4) <= 1))   # the file rounds to whole crab, so 34,361 is 0.4 from the replay
cc <- rd(Sys.glob(file.path(R, "cross_check_*.csv"))[1])
add("5", "the cross-check within tolerance", sprintf("pooled %s, gear %s, %+.2f%%, %s", format(cc$pooled_median, big.mark = ","), format(cc$gear_median, big.mark = ","), cc$gear_minus_pooled_pct, cc$verdict), identical(cc$verdict, "PASS"))
res <- do.call(rbind, rows); res$verdict <- if (all(res$pass)) "ADOPTED: every clause met" else "NOT MET: the key goes back to \"level\""
print(res[, c("clause", "what", "observed", "pass")], row.names = FALSE, right = FALSE); cat("\nA42:", res$verdict[1], "\n")
comp <- merge(rd(file.path(B70[["pooled"]], "pe_vs_bss_comparison.csv"))[, c("component", "BSS_catch")], pv[, c("component", "BSS_catch")], by = "component", suffixes = c("_B70", "_A42"))
comp$change_pct <- 100 * (comp$BSS_catch_A42 / comp$BSS_catch_B70 - 1)
pt <- function(d) { x <- rd(Sys.glob(file.path(d, "port_total_*.csv"))[1]); x[x$Estimate == "Expected_Catch", c("BSS_median", "BSS_lo95", "BSS_hi95")] }
port <- rbind(cbind(run = "B70", pt(B70[["pooled"]])), cbind(run = "A42", pt(NEW[["pooled"]])))
print(comp, row.names = FALSE); print(port, row.names = FALSE)
utils::write.csv(res, "05_output/a42_render_check_2026-11-25.csv", row.names = FALSE)
utils::write.csv(cbind(comp, port_B70 = port$BSS_median[1], port_A42 = port$BSS_median[2], port_A42_lo95 = port$BSS_lo95[2], port_A42_hi95 = port$BSS_hi95[2]),
                 "05_output/a42_render_check_2026-11-25_components.csv", row.names = FALSE)
mo <- merge(rd(file.path(B70[["pooled"]], "monthly_estimates_by_mode.csv")), rd(file.path(NEW[["pooled"]], "monthly_estimates_by_mode.csv")), by = c("Mode", "Month"), suffixes = c("_B70", "_A42"))
mo$change_pct <- 100 * (mo$median_A42 / mo$median_B70 - 1)
print(mo[mo$Mode != "Shore", c("Mode", "Month", "median_B70", "median_A42", "change_pct")], row.names = FALSE, digits = 4)
utils::write.csv(mo, "05_output/a42_render_check_2026-11-25_months.csv", row.names = FALSE)
cat("Wrote 05_output/a42_render_check_2026-11-25*.csv\n")
