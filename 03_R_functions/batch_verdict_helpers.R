###############################################################################
# SHARED BATCH-VERDICT HELPERS
# -----------------------------------------------------------------------------
# Extracted 2026-09-04 from run_shore_ar_zi_2026-09-03.R, which was the fourth batch
# runner to carry its own copy. Each copy was a place a defect could be fixed in one
# runner and left in the others, and two of them already were:
#
#   merge_csv_by()   exists because an UNGATED overwrite of the verdicts file during a
#                    dry run truncated the 2026-09-01 V-stage verdicts to 11 desk rows,
#                    destroying results that had cost hours of fitting; and because
#                    appending on RESUME had left osp_validation_summary.csv with 28 rows
#                    for 14 jobs. Merge by key, never append and never blind-overwrite.
#   config_delta()   exists because an exactness comparison is only meaningful when the
#                    two runs differ in ONE thing, and that was got wrong by hand twice in
#                    three weeks, both times producing a FAIL that read as a code defect.
#                    `expect_delta` turns the assumption into a check.
#
#   READ THIS BEFORE TRUSTING A config_delta RESULT ON AN OLD FOLDER. It reads
#   run_parameters.txt, and every such file written before 2026-09-04 is TRUNCATED:
#   str() caps a list at list.len = 99 and run_config carries 120+ keys. That yields
#   false positives (a key that crossed the cut between two runs reads as differing;
#   the 2026-09-03 Z0 stage reported 9 such keys on a correct config) and false
#   NEGATIVES (a key past entry 99 in both dumps is never compared). .cfg_keys() now
#   reports the truncation so a stale comparison says so instead of looking clean.
#
# fit_exactness() never compares PORT TOTALS. rstan::extract(permuted = TRUE) permutes
# draws, so a port total assembled by resampling component draws moves by about 0.2%
# between bit-identical fits. Exactness is always tested on per-fit summaries.
###############################################################################

# GUARDED, and NULL-only on purpose. The drivers source this whole folder AFTER
# library(tidyverse), so an unconditional definition here would SHADOW rlang's `%||%`
# (and base R's since 4.4) for the entire run. Both are null-only; a length-0 variant
# would change behaviour anywhere a config value is legitimately zero-length, and
# `opener_manual_shore` / `opener_manual_boat` are exactly chr(0) in production. The
# guard is the same idiom save_run_diagnostics.R has used since it was written.
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

merge_csv_by <- function(new, path, key) {
  if (is.null(new) || !nrow(new)) return(invisible(NULL))
  old <- if (file.exists(path))
    tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e) NULL) else NULL
  if (!is.null(old) && nrow(old) && all(key %in% names(old)) && all(key %in% names(new))) {
    k_old <- do.call(paste, c(old[key], sep = "\r"))
    k_new <- do.call(paste, c(new[key], sep = "\r"))
    old <- old[!(k_old %in% k_new), , drop = FALSE]
    cols <- intersect(names(old), names(new))
    # A schema change between runs is possible (a new column added mid-project). Keep the
    # NEW schema and carry the old rows across on the shared columns rather than dropping
    # them, so an old row survives with blanks instead of vanishing.
    if (nrow(old)) {
      pad <- old[, cols, drop = FALSE]
      for (cn in setdiff(names(new), cols)) pad[[cn]] <- NA
      new <- rbind(pad[, names(new), drop = FALSE], new)
    }
  }
  utils::write.csv(new, path, row.names = FALSE)
  invisible(new)
}

.cfg_keys <- function(dir) {
  p <- file.path(dir, "run_parameters.txt")
  if (!file.exists(p)) return(NULL)
  kv <- list()
  txt <- readLines(p, warn = FALSE)
  # 2026-09-04: surface a truncated dump rather than silently comparing a partial config.
  if (any(grepl("list output truncated", txt, fixed = TRUE)))
    attr(kv, "truncated") <- TRUE
  for (l in txt) {
    m <- regmatches(l, regexec("^\\s*\\$ ([A-Za-z0-9_.]+)\\s*:(.*)$", l))[[1]]
    if (length(m) == 3) kv[[m[2]]] <- trimws(m[3])
  }
  if (any(grepl("list output truncated", txt, fixed = TRUE))) attr(kv, "truncated") <- TRUE
  kv
}

config_delta <- function(dir_a, dir_b) {
  a <- .cfg_keys(dir_a); b <- .cfg_keys(dir_b)
  if (is.null(a) || is.null(b)) return(NA_character_)
  trunc <- isTRUE(attr(a, "truncated")) || isTRUE(attr(b, "truncated"))
  ks <- setdiff(union(names(a), names(b)), c("run_tag", "model"))
  out <- sort(ks[vapply(ks, function(k) !identical(a[[k]] %||% "<absent>", b[[k]] %||% "<absent>"), logical(1))])
  # A truncated dump on either side makes an "absent" key uninterpretable, so say so in
  # the result rather than letting the caller read it as a real difference.
  if (trunc) attr(out, "truncated") <- TRUE
  out
}

fit_exactness <- function(new_dir, ref_dir, pat = NULL, what = "fits", expect_delta = NULL) {
  if (!dir.exists(new_dir %||% "") || !dir.exists(ref_dir %||% ""))
    return(list(observed = "run or reference folder missing", verdict = "REVIEW"))
  fs <- c(list.files(new_dir, pattern = "^bss_summary_.*\\.csv$"),
          list.files(new_dir, pattern = "^bss_full_summary_.*\\.csv$"))
  if (!is.null(pat)) fs <- fs[grepl(pat, fs)]
  fs <- fs[file.exists(file.path(ref_dir, fs))]
  n <- 0L; bad <- character(0)
  for (f in fs) {
    a <- tryCatch(utils::read.csv(file.path(ref_dir, f), row.names = 1, check.names = FALSE), error = function(e) NULL)
    b <- tryCatch(utils::read.csv(file.path(new_dir, f), row.names = 1, check.names = FALSE), error = function(e) NULL)
    if (is.null(a) || is.null(b)) { bad <- c(bad, f); next }
    cr <- intersect(rownames(a), rownames(b)); cc <- intersect(names(a), names(b))
    n <- n + length(cr)
    if (!isTRUE(all.equal(a[cr, cc], b[cr, cc], tolerance = 0))) bad <- c(bad, f)
  }
  cd <- config_delta(ref_dir, new_dir)
  extra <- if (length(cd) && !all(is.na(cd))) setdiff(cd, expect_delta %||% character(0)) else character(0)
  if (isTRUE(attr(cd, "truncated")))
    extra <- setdiff(extra, extra)   # a truncated dump cannot support an UNEXPECTED claim
  note <- if (!length(cd) || all(is.na(cd))) "" else
    sprintf("; configs differ in %d key(s): %s%s%s", length(cd), paste(cd, collapse = ", "),
            if (isTRUE(attr(cd, "truncated")))
              "  [one or both run_parameters.txt dumps are TRUNCATED; absent-key differences are not interpretable]" else "",
            if (length(extra)) sprintf("  <-- %d UNEXPECTED: %s", length(extra), paste(extra, collapse = ", ")) else "")
  list(observed = sprintf("%s: %d shared parameter rows across %d summaries %s%s", what, n, length(fs),
                          if (!length(bad)) "identical at full precision"
                          else sprintf("DIFFER in %s", paste(basename(bad), collapse = ", ")), note),
       verdict = if (!length(bad) && n > 0) "PASS" else "FAIL", unexpected_delta = extra)
}
