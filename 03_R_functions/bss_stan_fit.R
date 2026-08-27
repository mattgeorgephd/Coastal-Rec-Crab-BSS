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
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program (see the LICENSE file);
# if not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
###############################################################################
# bss_stan_fit.R  --  guarded wrapper around rstan::stan().
#
# WHY THIS FILE EXISTS (2026-08-26).
# The 2026-08-25 patch added five data variables to both .stan models
# (OSPF_n, osp_f_stratum, osp_f_total, osp_f_crab, osp_f_kappa_prior_mu) and
# forwarded only three of the eight new crab_fraction_stan_data() fields out of
# prep_bss_crab_pooled()/prep_bss_crab_gear(). Stan therefore threw
#
#   Exception: variable does not exist; processing stage=data initialization;
#   variable name=OSPF_n; base type=int
#
# on EVERY fit of EVERY run. rstan does not raise that as an R error: it prints it
# to the MESSAGE stream, returns an empty stanfit (mode 2), and lets the caller
# carry on. The pipeline then walked three layers downstream before dying:
#
#   summary(fit, pars=...)          -> rstan returns invisible(NULL) for mode 2, so
#                                      $summary is NULL and bss_summary_*.csv is 3 bytes
#   quantile(NULL, c(...))          -> NA NA NA, no error
#   apply(extract(fit,"E")$E[,1,,1], 2, median)
#                                   -> "dim(X) must have a positive length"
#
# ...which is what the operator sees, in a chunk marked results='hide', with the
# real cause 300 lines upstream and on a stream knitr had already consumed. Six
# validation runs and roughly an hour of wall clock were spent on that message.
#
# Two guards, both cheap, both run on every fit:
#   1. bss_assert_stan_data() parses the .stan data block and refuses to call the
#      sampler when a declared variable is absent from the data list, or when a
#      supplied value is non-finite. Fails BEFORE the compile/sample, naming the
#      variable.
#   2. bss_assert_fit_usable() rejects an empty stanfit immediately and reports the
#      Stan message text that rstan sent to the message stream, instead of letting
#      NULL propagate into the extraction code.
#
# Neither guard changes any fit that succeeds: on the happy path bss_stan_fit()
# calls rstan::stan() with exactly the arguments it was given and returns exactly
# what rstan returned. The withCallingHandlers() message handler RECORDS and does
# not muffle, so sampler progress still reaches knitr as before.
###############################################################################

# ---------------------------------------------------------------------------
# Parse the `data { ... }` block of a Stan program and return the declared
# variable names. Handles the repo's old-style array declarations
# (`int period[D];`, `real O[D,S,G];`), the newer sized types
# (`vector<lower=0>[D] w;`, `matrix[D,K] X;`) and line/block comments.
#
# Deliberately a text parse, not a stanc call: it must work on any machine that
# can run the pipeline, without a stanc binary, and it only needs the names.
# ---------------------------------------------------------------------------
bss_stan_data_names <- function(stan_file) {
  if (!file.exists(stan_file)) return(character(0))
  txt <- readLines(stan_file, warn = FALSE)
  txt <- sub("//.*$", "", txt)                              # line comments
  s   <- paste(txt, collapse = "\n")
  s   <- gsub("/\\*.*?\\*/", "", s, perl = TRUE)            # block comments

  i <- regexpr("(^|\\n)\\s*data\\s*\\{", s, perl = TRUE)
  if (i < 0) return(character(0))
  start <- i + attr(i, "match.length") - 1L                 # index of the '{'
  ch <- strsplit(substring(s, start), "")[[1]]
  depth <- 0L; end <- NA_integer_
  for (k in seq_along(ch)) {
    if (ch[k] == "{") depth <- depth + 1L
    else if (ch[k] == "}") { depth <- depth - 1L; if (depth == 0L) { end <- k; break } }
  }
  if (is.na(end)) return(character(0))
  body  <- substring(s, start + 1L, start + end - 2L)
  stmts <- strsplit(body, ";", fixed = TRUE)[[1]]

  nm <- vapply(stmts, function(z) {
    z <- trimws(gsub("\\s+", " ", z))
    if (!nzchar(z)) return(NA_character_)
    z <- sub("\\[[^]]*\\]\\s*$", "", z)                     # trailing array dims
    m <- regmatches(z, regexpr("[A-Za-z_][A-Za-z0-9_]*\\s*$", z))
    if (!length(m)) NA_character_ else trimws(m)
  }, character(1), USE.NAMES = FALSE)

  unique(nm[!is.na(nm)])
}

# ---------------------------------------------------------------------------
# Refuse to sample with an incomplete or non-finite data list.
#
# Extra names in `stan_data` are fine and expected: the prep functions attach
# metadata (.effort_unit, .L_unit, .gear_type_labels, ...) that Stan ignores.
# ---------------------------------------------------------------------------
bss_assert_stan_data <- function(stan_data, stan_file, label = NULL) {
  need <- bss_stan_data_names(stan_file)
  tag  <- if (is.null(label)) "" else sprintf(" [%s]", label)

  if (!length(need)) {
    warning(sprintf("bss_assert_stan_data%s: could not parse a data block from %s; ",
                    tag, basename(stan_file)),
            "skipping the completeness check.", call. = FALSE)
    return(invisible(FALSE))
  }

  miss <- setdiff(need, names(stan_data))
  if (length(miss))
    stop(sprintf(paste0("bss_assert_stan_data%s: %d variable(s) declared in %s are ",
                        "missing from the data list:\n    %s\n",
                        "  Stan would fail at data initialization, return an EMPTY stanfit, ",
                        "and the failure would only surface downstream. Add them in the prep ",
                        "function that builds this list."),
                 tag, length(miss), basename(stan_file), paste(miss, collapse = ", ")),
         call. = FALSE)

  supplied <- stan_data[intersect(need, names(stan_data))]
  bad <- names(supplied)[vapply(supplied, function(z)
    is.numeric(z) && length(z) > 0L && any(!is.finite(z)), logical(1))]
  if (length(bad))
    stop(sprintf(paste0("bss_assert_stan_data%s: non-finite value(s) (NA/NaN/Inf) in Stan ",
                        "data variable(s): %s. Stan would reject these at data ",
                        "initialization."),
                 tag, paste(bad, collapse = ", ")), call. = FALSE)

  not_num <- names(supplied)[!vapply(supplied, function(z)
    is.numeric(z) || is.integer(z) || is.logical(z), logical(1))]
  if (length(not_num))
    stop(sprintf(paste0("bss_assert_stan_data%s: non-numeric value(s) supplied for Stan data ",
                        "variable(s): %s."),
                 tag, paste(not_num, collapse = ", ")), call. = FALSE)

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# An empty stanfit is a failure, not a result. rstan signals it by returning an
# object with mode 2 (or mode 1 for test_grad) and an empty @sim; summary() on such
# an object returns invisible(NULL) and extract() returns invisible(NULL), so every
# downstream reader silently degrades instead of stopping.
# ---------------------------------------------------------------------------
bss_fit_has_draws <- function(fit) {
  inherits(fit, "stanfit") &&
    isTRUE(fit@mode == 0L) &&
    is.list(fit@sim) && length(fit@sim) > 0L &&
    !is.null(fit@sim$samples) && length(fit@sim$samples) > 0L
}

# ---------------------------------------------------------------------------
# Recover the real Stan error after an empty fit.
#
# WHY THIS IS NEEDED. At cores > 1 rstan runs each chain in a FORKED WORKER, and the
# lines that actually name the cause ('Error : Initialization failed.', the per-chain
# exception text) are printed by those workers straight to the terminal file
# descriptor. No R-level sink or condition handler in the parent can see them, so a
# parallel run's most useful diagnostic is unrecoverable by construction.
#
# rstan's own advice in that situation is "consider specifying chains = 1 to debug".
# This does exactly that, automatically: one chain, in THIS process, two iterations. A
# deterministic failure (a bad data value, a zero prior scale, an out-of-range index)
# reproduces instantly and its text is captured; the compiled model is already cached,
# so the cost is seconds. A failure that does NOT reproduce is itself informative and
# is reported as such rather than papered over.
# ---------------------------------------------------------------------------
bss_diagnose_empty_fit <- function(file, data) {
  msgs <- character(0)
  printed <- tryCatch(
    utils::capture.output(
      withCallingHandlers(
        suppressWarnings(rstan::stan(file = file, data = data, chains = 1, cores = 1,
                                     iter = 2, warmup = 1, seed = 1,
                                     open_progress = FALSE)),
        message = function(m) msgs <<- c(msgs, conditionMessage(m))),
      type = "output", split = FALSE),
    error = function(e) paste("diagnostic refit threw:", conditionMessage(e)))
  out <- sub("\n+$", "", c(printed, msgs))
  out <- out[nzchar(trimws(out))]
  # An init failure repeats its exception once per attempt (rstan retries 100 times), so
  # the raw text is the same two lines a hundred times over. Distinct lines only, in order.
  n_raw <- length(out)
  out   <- unique(out)
  if (length(out) > 25L)
    out <- c(out[1:25], sprintf("... (%d further distinct lines suppressed)", length(out) - 25L))
  if (n_raw > length(out))
    out <- c(out, sprintf("(%d repeated lines collapsed; rstan retries initialization 100 times)",
                          n_raw - length(out)))
  if (!length(out))
    out <- paste("(the single-chain refit produced no error text and may have succeeded;",
                 "the parallel failure is intermittent, not structural)")
  out
}

bss_assert_fit_usable <- function(fit, label = NULL, console = character(0)) {
  tag <- if (is.null(label)) "" else sprintf(" [%s]", label)

  if (bss_fit_has_draws(fit)) return(invisible(fit))

  said <- console[grepl("error|exception|fail|does not exist|Initialization",
                        console, ignore.case = TRUE)]
  if (!length(said)) said <- utils::tail(console, 12)
  said <- unique(said)
  if (length(said) > 20L)
    said <- c(said[1:20], sprintf("... (%d more lines; see the fit's stan_console log)",
                                  length(said) - 20L))
  detail <- if (length(said))
    paste0("\n  Stan said:\n    ", paste(sub("\n+$", "", said), collapse = "\n    "))
  else
    "\n  (rstan printed nothing to the message stream; run the fit outside knitr to see it.)"

  stop(sprintf(paste0("bss_assert_fit_usable%s: rstan returned an EMPTY stanfit (mode %s); ",
                      "no chain produced draws.%s"),
               tag,
               if (inherits(fit, "stanfit")) as.character(fit@mode) else "n/a",
               detail),
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# The call the drivers make. Signature mirrors rstan::stan() for the arguments the
# pipeline uses; `label` and `log_file` are ours.
#
#   label     fit label, used only in guard messages
#   log_file  optional path; whatever rstan reports in THIS process during the call is
#             written there, so a failure stays diagnosable from output_dir after the
#             fact even when the chunk is results='hide'
#
# WHAT THE LOG DOES AND DOES NOT CONTAIN. rstan splits its reporting across both R
# streams, and the useful half is not where you would expect:
#   * MESSAGE stream (parent): the data-initialization exception and the per-chain
#     "Initialization between (-2, 2) failed" text.
#   * OUTPUT stream (parent): 'Error : Initialization failed.', "error occurred during
#     calling the sampler; sampling not done", and the per-chain "does not contain
#     samples" list -- printed, not messaged, so a message handler alone misses them.
# Both are captured here. With cores > 1 the per-chain "Chain n: Iteration ..." progress
# is written by FORKED WORKERS straight to the terminal, bypasses both parent streams and
# is therefore NOT in the log; that is progress, not diagnosis, and is not worth
# serialising the run to collect. capture.output(split = TRUE) tees rather than swallows,
# so nothing that used to appear on the console stops appearing.
#
# CONSEQUENCE WORTH KNOWING: at cores > 1 a fit that SUCCEEDS usually reports nothing in
# this process, so no log file is written for it. A log file present next to a fit's other
# outputs means rstan had something to say; a failure always produces one.
# ---------------------------------------------------------------------------
bss_stan_fit <- function(file, data, label = NULL, log_file = NULL,
                         diagnose_empty = TRUE, ...) {

  bss_assert_stan_data(data, file, label = label)

  console <- character(0)
  fit <- NULL
  printed <- utils::capture.output(
    fit <- withCallingHandlers(
      rstan::stan(file = file, data = data, ...),
      # record only; no invokeRestart("muffleMessage"), so sampler messages and warnings
      # still reach the console / knitr exactly as before.
      message = function(m) console <<- c(console, conditionMessage(m))),
    type = "output", split = TRUE)

  console <- c(printed, console)

  # At cores > 1 the lines that name the cause were printed by forked workers and never
  # reached this process. Reproduce the failure in-process on one chain and keep THAT.
  if (isTRUE(diagnose_empty) && !bss_fit_has_draws(fit))
    console <- c(console, "",
                 "--- single-chain diagnostic refit (chains = 1, cores = 1, iter = 2) ---",
                 bss_diagnose_empty_fit(file, data))

  if (!is.null(log_file) && length(console))
    try(writeLines(sub("\n+$", "", console), log_file), silent = TRUE)

  bss_assert_fit_usable(fit, label = label, console = console)
  fit
}
