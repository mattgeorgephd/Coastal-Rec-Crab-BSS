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
# desk_stan_array_syntax_2026-11-25.R
#
# DESK, A MINUTE, NO MCMC (CHANGE_REGISTER D53 -> B84).
#
# THE QUESTION. Both Stan programs declared arrays the pre-2.26 way (`int<lower=1> period[D];`),
# a form Stan removed at 2.33. CRAN now serves StanHeaders 2.39.1.1, whose parser rejects it, so a
# machine set up without renv cannot compile either model (D53). The repair is to write every
# such declaration as `array[D] int<lower=1> period;`, which every Stan from 2.26 accepts. It is a
# change of syntax only, and this script is the proof that it is nothing more.
#
# THE CONVERSION. Mechanical and complete by construction: the lock's parser (stanc3 2.32.2, the
# one StanHeaders 2.32.10 ships) reports the line and column of every deprecated declaration, and
# each reported declaration, and nothing else, is rewritten. The programs before the conversion
# are read from git at OLD_REF below; the programs after it are the working tree.
#
# THE RULE, WRITTEN AND COMMITTED BEFORE THE PROGRAMS WERE CONVERTED.
#   (1) Baseline. The programs at OLD_REF parse under the 2.32 parser with exactly 50 (pooled)
#       and 53 (gear-resolved) deprecated-array warnings and no other warning, and fail under the
#       2.39 parser. (This is the defect as recorded; if it does not reproduce, stop.)
#   (2) The converted programs parse under the 2.32 parser with NO warning of any kind.
#   (3) The converted programs parse under the 2.39 parser (CRAN's StanHeaders 2.39.1.1, its
#       stanc.js read from the CRAN source tarball, SHA-256 recorded) with no error. Its warnings
#       are reported.
#   (4) THE PROOF. The C++ the 2.32 parser generates from the converted program is identical to
#       the C++ it generates from the old one, character for character, after one normalisation:
#       each source-location string "(in '<name>', line L, column a to column b)" is replaced by a
#       placeholder, because the conversion moves columns on the rewritten lines and those strings
#       feed only the text of a runtime error message. The comparison is made with rstan's own
#       default stanc flags and with the optimisation flag (O1) as well. Identical C++ compiles to
#       the same model, so every fit under either program is the same fit draw for draw.
#   ACCEPTED if (1) to (4) hold, for both programs. Anything else in (4), however small, is a
#   FAIL: the conversion is reverted and the difference investigated.
#   (5) CONFIRMATION, by run (the project's rule that a change is not done until a run confirms
#       it): `Rscript run_estimation.R --model both` on run_config.R as shipped, which is the
#       authoritative run's configuration, reproduces the committed bss_draws_summed_*.csv of B70
#       (05_output/20261019/pooled-CPUE-2024-25 and gear-type-CPUE-model-2024-25) byte for byte on
#       all eight fits, and its port totals. The render may carry D54's diagnostic-only changes;
#       they change no estimate, and this clause would catch it if they did. On this container the
#       B80 arm S25-const reproduced B70 byte for byte on 2026-11-19, so a difference would be
#       this change's. Read by compare_render() below once the render exists.
#   (6) REPORTED, NOT GATING: the converted pooled program compiled and sampled briefly under
#       StanHeaders 2.39.1.1 with rstan 2.32.7 rebuilt against it, in a separate library. It
#       shows the programs work on CRAN's current toolchain; it cannot gate, because a failure
#       there could be rstan's own and not the programs'.
#   EXPECTED: (1) to (5) hold.
#
# READ 2026-11-25: (1) to (4) hold on both programs (ACCEPTED). Then clause (6) ran.
#
# ADDENDUM (4b), WRITTEN AFTER CLAUSE (6) RAN AND BEFORE ANY RENAME. Under CRAN's pair
# (StanHeaders 2.39.1.1, rstan 2.32.7 rebuilt against it) the converted pooled program COMPILES
# (2.7 minutes) but cannot start: "variable does not exist; processing stage=data
# initialization; variable name=sun". Stan 2.39's parser writes the data variable `sun` as
# `_stan_sun` in its C++ (sun is a predefined macro on some compilers), and rstan 2.32.7 then does
# not find it. A three-line model reproduces it and runs when the variable is called `sunday`;
# on the lock's pair it runs either way. `sun` is the only identifier 2.39 renames in either
# program. So the second barrier is the NAME, and the repair is to call the A36 Sunday indicator
# `sunday` in both programs and both preps (and in the one diagnostic that reads it by name).
#   (4b) The C++ the 2.32 parser generates from the renamed program is identical, up to
#        source-location strings, to the C++ it generates from the program as converted (the
#        commit whose message starts "B84 conversion") once the identifier sun is written
#        sunday there (whole word, and the data name string "sun"). With rstan's default flags
#        and with O1. Anything else is a FAIL and the rename is reverted.
#   (6) is then re-run on the renamed program, and REPORTED: the pooled program samples 300
#        iterations of the 2024-25 boat pot-closure fit on CRAN's pair.
#   Clause (5), the confirmation render, is run on the renamed programs, so it covers both.
#
# WRITES 05_output/stan_array_syntax_2026-11-25.csv (one row per clause and program).
#
# HOW TO RUN (from the repository root; needs network once for the CRAN tarball)
#   RENV_ACTIVATE_PROJECT=FALSE Rscript 06_diagnostics/desk_stan_array_syntax_2026-11-25.R
#   ... --convert             rewrite the working-tree programs first (done once, 2026-11-25)
#   ... --render=<pooled dir>,<gear dir>   clause (5) against a finished render
###############################################################################
.args <- commandArgs(trailingOnly = TRUE)
.root <- getwd()
if (!dir.exists(file.path(.root, "02_stan_models")) && dir.exists(file.path(.root, "..", "02_stan_models")))
  .root <- normalizePath(file.path(.root, ".."))
if (!dir.exists(file.path(.root, "02_stan_models"))) stop("Run this from the repository root (or from 06_diagnostics/).")
setwd(.root)
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

OLD_REF <- "d963eb3"                     # the branch tip before the conversion
PROGRAMS <- c(pooled = "02_stan_models/crab_bss_pooled.stan", gear_resolved = "02_stan_models/crab_bss_gear_resolved.stan")
EXPECTED_DEPRECATED <- c(pooled = 50L, gear_resolved = 53L)
SH239_URL <- "https://cloud.r-project.org/src/contrib/StanHeaders_2.39.1.1.tar.gz"
B70 <- c(pooled = "05_output/20261019/pooled-CPUE-2024-25", gear_resolved = "05_output/20261019/gear-type-CPUE-model-2024-25")
OUT <- "05_output/stan_array_syntax_2026-11-25.csv"
DEP_MSG <- "brackets after a variable name"

# ---------------------------------------------------------------------------
# A stanc.js of any version, run the way rstan runs its own (QuickJSR, the same call).
# ---------------------------------------------------------------------------
stanc_run <- function(js, code, flags = character(0), name = "model") {
  ctx <- QuickJSR::JSContext$new(stack_size = 4 * 1024 * 1024)
  ctx$source(js)
  r <- ctx$call("stanc", name, code, as.array(flags))
  list(errors = unlist(r$errors), warnings = unlist(r$warnings), cpp = r$result %||% "")
}
stanc_version <- function(js) { x <- readLines(js, warn = FALSE); m <- regmatches(x, regexpr("stanc3 v[0-9.]+", x)); if (length(m)) m[1] else NA_character_ }
read_old <- function(path) paste(system2("git", c("show", sprintf("%s:%s", OLD_REF, path)), stdout = TRUE), collapse = "\n")
# The programs the rule certifies are those AT THE RENAME COMMIT (the message starting "B84 rename"), not the
# working tree: later deliberate changes (B86 added the D55 switch in 2442411) are not this rule's to judge, and
# a re-run on the working tree after them read as differences (read 2026-11-25). Before the rename exists, and
# with --convert, the working tree is read.
.new_ref <- local({ x <- suppressWarnings(system2("git", c("log", "-1", "--format=%h", shQuote("--grep=^B84 rename")), stdout = TRUE)); if (length(x) == 1L && nzchar(x)) x else NA_character_ })
read_new <- function(path) if (!is.na(.new_ref) && !("--convert" %in% .args)) paste(system2("git", c("show", sprintf("%s:%s", .new_ref, path)), stdout = TRUE), collapse = "\n") else paste(readLines(path, warn = FALSE), collapse = "\n")
# the identifier sun as the 2.32 parser writes it in C++: sun, and sun_data__ / sun_flat__ (a trailing
# underscore is still the same name), never B1_sun or sunday_split
SUN_MAP <- "(?<![A-Za-z0-9_])sun(?![A-Za-z0-9])"
# and because sunday is three characters longer, stanc's C++ printer indents the continuation lines
# of a call that names it three columns further: whitespace outside string literals means nothing to
# C++, so the renamed comparisons collapse every run of white space first (read 2026-11-25; the first
# two reads of 4b failed on the map and then on this indentation, both the comparison's, not the C++'s)
ws <- function(x) gsub("[[:space:]]+", " ", x)
norm_cpp <- function(cpp) gsub("\\(in '[^']*', line [0-9]+, column [0-9]+ to (line [0-9]+, )?column [0-9]+\\)", "(in <loc>)", cpp)

js232 <- system.file("stanc.js", package = "StanHeaders")
if (!nzchar(js232)) stop("StanHeaders is not installed")
if (as.character(utils::packageVersion("StanHeaders")) != "2.32.10") stop("this check needs the lock's StanHeaders, 2.32.10")
fetch_js239 <- function() {
  d <- file.path(tempdir(), "sh239"); dir.create(d, showWarnings = FALSE)
  tb <- file.path(d, "StanHeaders_2.39.1.1.tar.gz")
  if (!file.exists(tb)) utils::download.file(SH239_URL, tb, quiet = TRUE, mode = "wb")
  utils::untar(tb, files = "StanHeaders/inst/stanc.js", exdir = d)
  list(js = file.path(d, "StanHeaders/inst/stanc.js"), sha256 = if (requireNamespace("digest", quietly = TRUE)) digest::digest(file = tb, algo = "sha256") else NA_character_)
}

# ---------------------------------------------------------------------------
# --convert: rewrite exactly the declarations the 2.32 parser reports, nothing else.
# A declaration at (line, col) reads  TYPE NAME[DIMS]  where TYPE is an identifier, then an
# optional <...> and an optional [...] (both bracket-balanced); it becomes  array[DIMS] TYPE NAME.
# ---------------------------------------------------------------------------
balanced <- function(s, i, open, close) {      # s a character vector of single chars, s[i] == open; index of its match
  depth <- 0L
  for (k in i:length(s)) { if (s[k] == open) depth <- depth + 1L else if (s[k] == close) { depth <- depth - 1L; if (depth == 0L) return(k) } }
  stop("unbalanced ", open)
}
convert_line <- function(line, col0) {
  s <- strsplit(line, "")[[1]]; i <- col0 + 1L
  st <- i
  while (i <= length(s) && grepl("[A-Za-z0-9_]", s[i])) i <- i + 1L          # the base type
  if (i <= length(s) && s[i] == "<") i <- balanced(s, i, "<", ">") + 1L       # constraints
  if (i <= length(s) && s[i] == "[") i <- balanced(s, i, "[", "]") + 1L       # vector / matrix size
  type <- paste(s[st:(i - 1L)], collapse = "")
  while (i <= length(s) && s[i] == " ") i <- i + 1L
  ns <- i
  while (i <= length(s) && grepl("[A-Za-z0-9_]", s[i])) i <- i + 1L
  name <- paste(s[ns:(i - 1L)], collapse = "")
  while (i <= length(s) && s[i] == " ") i <- i + 1L
  if (i > length(s) || s[i] != "[") stop(sprintf("no array dimensions after '%s' in: %s", name, line))
  e <- balanced(s, i, "[", "]")
  dims <- paste(s[(i + 1L):(e - 1L)], collapse = "")
  paste0(if (st > 1L) paste(s[1:(st - 1L)], collapse = "") else "", "array[", dims, "] ", type, " ", name,
         if (e < length(s)) paste(s[(e + 1L):length(s)], collapse = "") else "")
}
if ("--convert" %in% .args) {
  for (m in names(PROGRAMS)) {
    lines <- readLines(PROGRAMS[[m]], warn = FALSE)
    w <- stanc_run(js232, paste(lines, collapse = "\n"))$warnings
    dep <- grep(DEP_MSG, w, value = TRUE)
    pos <- do.call(rbind, lapply(regmatches(dep, regexec("line ([0-9]+), column ([0-9]+)", dep)), function(x) as.integer(x[2:3])))
    if (is.null(pos)) { cat(m, ": nothing to convert\n"); next }
    if (anyDuplicated(pos[, 1])) stop("two deprecated declarations on one line: convert by hand")
    for (r in seq_len(nrow(pos))) lines[pos[r, 1]] <- convert_line(lines[pos[r, 1]], pos[r, 2])
    writeLines(lines, PROGRAMS[[m]])
    cat(sprintf("%s: %d declarations rewritten\n", m, nrow(pos)))
  }
}

# ---------------------------------------------------------------------------
# Clauses (1) to (4).
# ---------------------------------------------------------------------------
j239 <- fetch_js239()
cat(sprintf("parsers: %s (StanHeaders 2.32.10) and %s (StanHeaders 2.39.1.1, tarball sha256 %s)\n", stanc_version(js232), stanc_version(j239$js), j239$sha256))
rows <- list()
add <- function(program, clause, what, observed, pass) rows[[length(rows) + 1L]] <<- data.frame(program = program, clause = clause, what = what, observed = observed, pass = pass, stringsAsFactors = FALSE)
for (m in names(PROGRAMS)) {
  old <- read_old(PROGRAMS[[m]]); new <- read_new(PROGRAMS[[m]])
  o32 <- stanc_run(js232, old); o39 <- stanc_run(j239$js, old)
  ndep <- sum(grepl(DEP_MSG, o32$warnings)); nother <- sum(!grepl(DEP_MSG, o32$warnings))
  add(m, "1", "the old program: 2.32 parses with the recorded deprecated-array warnings and no other; 2.39 rejects it",
      sprintf("2.32: %d errors, %d deprecated-array warnings, %d other; 2.39: %d errors", length(o32$errors), ndep, nother, length(o39$errors)),
      length(o32$errors) == 0 && ndep == EXPECTED_DEPRECATED[[m]] && nother == 0 && length(o39$errors) > 0)
  n32 <- stanc_run(js232, new)
  add(m, "2", "the converted program: 2.32 parses it with no warning", sprintf("%d errors, %d warnings", length(n32$errors), length(n32$warnings)),
      length(n32$errors) == 0 && length(n32$warnings) == 0)
  n39 <- stanc_run(j239$js, new)
  add(m, "3", "the converted program: 2.39 parses it", sprintf("%d errors, %d warnings%s", length(n39$errors), length(n39$warnings),
      if (length(n39$warnings)) paste0(": ", paste(unique(sub("\\n.*", "", sub("^Warning in [^:]*: ", "", n39$warnings))), collapse = " | ")) else ""),
      length(n39$errors) == 0)
  for (fl in list(character(0), "O1")) {
    co <- stanc_run(js232, old, fl)$cpp; cn <- n32$cpp
    if (length(fl)) cn <- stanc_run(js232, new, fl)$cpp
    # once the rename of (4b) is in the tree, the old program's identifier sun is mapped to sunday
    # before the comparison (4b itself checks the rename against the converted program alone)
    renamed <- grepl("\\bsunday\\b", new, perl = TRUE)
    raw_same <- identical(co, cn)
    nsame <- if (renamed) identical(ws(gsub(SUN_MAP, "sunday", norm_cpp(co), perl = TRUE)), ws(norm_cpp(cn))) else identical(norm_cpp(co), norm_cpp(cn))
    nloc <- length(gregexpr("\\(in '[^']*', line [0-9]+", co)[[1]])
    first <- if (nsame) "" else { a <- strsplit(if (renamed) gsub(SUN_MAP, "sunday", norm_cpp(co), perl = TRUE) else norm_cpp(co), "\n")[[1]]; b <- strsplit(norm_cpp(cn), "\n")[[1]]; k <- which(a[seq_len(min(length(a), length(b)))] != b[seq_len(min(length(a), length(b)))])[1]
      sprintf("; first difference at C++ line %s: [%s] vs [%s]", k, substr(a[k], 1, 120), substr(b[k], 1, 120)) }
    add(m, "4", sprintf("generated C++ identical up to source-location strings%s (flags: %s)", if (renamed) ", the identifier sun -> sunday and white space" else "", if (length(fl)) fl else "rstan default"),
        sprintf("%s characters; raw text %s; %d location strings normalised; normalised text %s%s", format(nchar(co), big.mark = ","),
                if (raw_same) "identical" else "differs", nloc, if (nsame) "IDENTICAL" else "DIFFERS", first),
        nchar(co) > 0 && nsame)
  }
}
# (4b): the rename of sun to sunday, against the program as converted
conv_ref <- system2("git", c("log", "-1", "--format=%h", shQuote("--grep=^B84 conversion")), stdout = TRUE)   # system2 does not quote: the space needs it
if (length(conv_ref) == 1L && nzchar(conv_ref)) for (m in names(PROGRAMS)) {
  conv <- paste(system2("git", c("show", sprintf("%s:%s", conv_ref, PROGRAMS[[m]])), stdout = TRUE), collapse = "\n")
  new <- read_new(PROGRAMS[[m]])
  if (!grepl("\\bsunday\\b", new, perl = TRUE)) next                     # not renamed yet
  for (fl in list(character(0), "O1")) {
    cc <- stanc_run(js232, conv, fl)$cpp; cn <- stanc_run(js232, new, fl)$cpp
    cc_m <- gsub(SUN_MAP, "sunday", norm_cpp(cc), perl = TRUE)
    same <- identical(ws(cc_m), ws(norm_cpp(cn)))
    add(m, "4b", sprintf("the rename sun -> sunday: C++ identical to the converted program's with the identifier mapped (flags: %s)", if (length(fl)) fl else "rstan default"),
        sprintf("against %s; %s", conv_ref, if (same) "IDENTICAL" else "DIFFERS"), nchar(cc) > 0 && same)
  }
}
res <- do.call(rbind, rows)

# ---------------------------------------------------------------------------
# Clause (5), against a finished render: --render=<pooled dir>,<gear dir>
# ---------------------------------------------------------------------------
.rd <- sub("^--render=", "", grep("^--render=", .args, value = TRUE))
if (length(.rd)) {
  dirs <- setNames(strsplit(.rd[1], ",")[[1]], names(B70))
  for (m in names(B70)) {
    ref <- Sys.glob(file.path(B70[[m]], "bss_draws_summed_*.csv"))
    same <- vapply(ref, function(f) { g <- file.path(dirs[[m]], basename(f)); file.exists(g) && identical(unname(tools::md5sum(f)), unname(tools::md5sum(g))) }, logical(1))
    # the port totals: every column B70's file has, value for value. (Read 2026-11-25: the render's file also
    # carries PE_se, the column B73 added on 2026-10-19 after B70 was rendered; it is reported, not compared.)
    pt <- function(d) { f <- Sys.glob(file.path(d, "port_total_*.csv")); if (length(f)) utils::read.csv(f[1], stringsAsFactors = FALSE, check.names = FALSE) else NULL }
    a <- pt(B70[[m]]); b <- pt(dirs[[m]])
    pt_same <- !is.null(a) && !is.null(b) && all(names(a) %in% names(b)) && nrow(a) == nrow(b) && identical(a, b[, names(a)])
    add(m, "5", "the confirmation render reproduces B70's saved draws byte for byte, and its port totals value for value",
        sprintf("%d of %d bss_draws_summed files identical; port totals %s%s", sum(same), length(ref), if (pt_same) "identical in every column B70's file has" else "DIFFER",
                if (!is.null(b) && length(setdiff(names(b), names(a)))) sprintf(" (the render's file adds %s)", paste(setdiff(names(b), names(a)), collapse = ", ")) else ""),
        all(same) && length(ref) == 4L && pt_same)
    res <- do.call(rbind, rows)
  }
}
verdict <- if (all(res$pass[res$clause %in% c("1", "2", "3", "4", "4b")])) "ACCEPTED" else "NOT ACCEPTED: revert the conversion"
op <- options(width = 250); print(res[, c("program", "clause", "observed", "pass")], row.names = FALSE, right = FALSE); options(op)
cat("\nD53 rule, clauses (1) to (4):", verdict, "\n")
res$verdict <- verdict
res$parser_232 <- stanc_version(js232); res$parser_239 <- stanc_version(j239$js); res$sh239_tarball_sha256 <- j239$sha256
utils::write.csv(res, OUT, row.names = FALSE)
cat("Wrote", OUT, "\n")
