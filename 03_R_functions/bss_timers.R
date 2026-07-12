###############################################################################
# bss_timers.R
#
# Section timers used by both drivers to log per-stage wall-clock time. Extracted
# from the pooled and gear-resolved drivers (identical inline copies) so the two
# share one implementation. Auto-sourced by both drivers via the 03_R_functions
# walk.
#
# Usage in a driver:
#   timer_start("data prep"); ...; timer_stop("data prep")
#   bss_timer_log()   # named list of {start, end, elapsed} per label, consumed by
#                     # the end-of-run TIMING SUMMARY table.
#
# The log lives in a module-local environment (not a driver global), so these
# functions carry no dependency on a `timer_log` object in the render environment.
# Sourcing this file resets the log, which is the correct behavior at the start of
# each render (each driver sources 03_R_functions in its setup chunk).
###############################################################################

.bss_timer <- new.env(parent = emptyenv())
.bss_timer$log <- list()

timer_start <- function(label) {
  .bss_timer$log[[label]] <- list(start = Sys.time())
  invisible()
}

timer_stop <- function(label) {
  .bss_timer$log[[label]]$end <- Sys.time()
  elapsed <- difftime(.bss_timer$log[[label]]$end,
                      .bss_timer$log[[label]]$start, units = "mins")
  .bss_timer$log[[label]]$elapsed <- elapsed
  cat(sprintf("  [TIMER] %s: %.1f min\n", label, as.numeric(elapsed)))
  invisible()
}

# Accessor for the collected timings (named list keyed by label).
bss_timer_log <- function() .bss_timer$log
