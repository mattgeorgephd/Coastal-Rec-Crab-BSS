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
# ---------------------------------------------------------------------------
# bss_sampler_override.R
#
# THE ONE SANCTIONED WAY FOR run_config TO REACH A SAMPLER SETTING.
#
# Both drivers merge their internal tuning ON TOP of run_config:
#
#     params <- modifyList(run_config, params_model)     # params_model WINS
#
# That direction is deliberate and documented (07_documentation/CLAUDE.md): per-fit
# sampler settings legitimately differ between the two tracks (pooled runs the boat
# all-gear fit at 5,000/2,500, gear-resolved runs everything at 2,000/1,000), and a flat
# run_config cannot hold two values for one key. The cost is that a key present in BOTH
# lists is decided silently by the driver. That is not hypothetical: `bss_min_interviews`
# sat at 20 while run_config appeared to own it, undetected until the 2026-08-25 audit.
#
# It bites again the moment an experiment needs a different iteration count. Plan item 5.2
# calls for re-running the GEAR track at 2,500 post-warmup draws, because its boat all-gear
# fit runs at 1,000 and that is where chain 3 stalled at the phi_E unit-root boundary
# (n_eff 40-49, R-hat 1.07-1.17, gate REJECTED, 2026-08-29 batch stage E). Setting
# `bss_iter_default` in a run_config delta does nothing at all: params_model overwrites it
# and the run looks like it honoured the request. A batch runner would then compare a
# 1,000-draw fit against a 1,000-draw fit and report that more draws did not help.
#
# `bss_sampler_override` is an explicit, narrow, LOUD escape hatch. It is a named list in
# run_config, applied AFTER the merge, and it is restricted to sampler keys: it cannot be
# used to smuggle a structural toggle past the merge, which would defeat the point of
# having the merge direction be a rule.
#
#   run_config$bss_sampler_override <- list(bss_iter_default = 5000,
#                                           bss_warmup_default = 2500)
#
# Deliberate design choices:
#   * UNRECOGNISED KEYS ERROR, they are not dropped. Silence is the failure mode this
#     function exists to remove; a typo that quietly does nothing is the same bug wearing
#     a different hat.
#   * It PRINTS what it changed, old value to new, into the run log, so a folder's console
#     records that the run was not at the driver's defaults.
#   * NULL / absent is a no-op, so production is untouched. It ships unset.
#
# NOT for production. Production sampler settings belong in the driver's params_model where
# they are version-controlled next to the fits they serve. This is for experiments.
# ---------------------------------------------------------------------------

# Sampler keys only. Matches bss_iter_*, bss_warmup_*, bss_treedepth_*, bss_delta_*,
# bss_max_treedepth_*, bss_adapt_delta*, bss_chains, bss_cores, bss_max_interviews.
.BSS_SAMPLER_OVERRIDE_PATTERN <-
  "^bss_(iter|warmup|treedepth|delta|max_treedepth|adapt_delta)(_|$)|^bss_(chains|cores|max_interviews)$"

bss_apply_sampler_override <- function(params, override = NULL, model_label = "",
                                       quiet = FALSE) {
  ov <- override %||% params$bss_sampler_override
  if (is.null(ov) || !length(ov)) return(params)
  if (!is.list(ov) || is.null(names(ov)) || any(!nzchar(names(ov))))
    stop("bss_sampler_override must be a NAMED list, e.g. list(bss_iter_default = 5000).",
         call. = FALSE)

  bad <- names(ov)[!grepl(.BSS_SAMPLER_OVERRIDE_PATTERN, names(ov))]
  if (length(bad))
    stop("bss_sampler_override accepts sampler keys only; refusing: ",
         paste(bad, collapse = ", "),
         ". Structural toggles must go through run_config and the documented merge, not ",
         "around it. Sampler keys match ", .BSS_SAMPLER_OVERRIDE_PATTERN, ".",
         call. = FALSE)

  if (!isTRUE(quiet))
    cat(sprintf("\n  SAMPLER OVERRIDE (%s): run_config is overriding %d driver default(s).\n",
                model_label, length(ov)))
  for (k in names(ov)) {
    old <- params[[k]]
    if (!isTRUE(quiet))
      cat(sprintf("    %-30s %s -> %s\n", k,
                  if (is.null(old)) "unset" else paste(format(old), collapse = ", "),
                  paste(format(ov[[k]]), collapse = ", ")))
    params[[k]] <- ov[[k]]
  }
  if (!isTRUE(quiet))
    cat("    These are EXPERIMENT settings. A run that used them is not comparable to one\n",
        "   that did not, on anything the sampler touches.\n", sep = "")
  params
}
