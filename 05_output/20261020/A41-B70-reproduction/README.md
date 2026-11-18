# A41's reproduction of B70 (B76, 2026-10-20)

What is kept of the 2024-25 `--model both` re-sample run under A41 on the committed tree `fa4977b`
(08:10Z to 10:30Z on the cloud container; the manifest `run_manifest_20261020_081004.txt` and the
cross-check `cross_check_20261020_081004.csv` are copies of the ones at `05_output/20261020/`).

Every one of its eight `bss_draws_summed_*.csv` files was byte-identical to B70's
(`05_output/20261019/pooled-CPUE-2024-25`, `gear-type-CPUE-model-2024-25`): the port 86,497
[72,869, 104,595] and the cross-check 86,285 (-0.25%, PASS) exactly. That is A41's clause (1),
that no 2024-25 fit changes, met by run. The full output folders were 160 MB of duplicates of
B70's and were deleted; the fits sit in the container's fit cache under the current (B71) key.

Files: `run_warnings_pooled.csv`, `run_warnings_gear.csv` (the pooled stage loaded its shore
pot-closure fit from the save of the stopped cache render and sampled the other three; the gear
stage sampled all four), `port_total_pooled.csv`, `port_total_gear.csv`. See CHANGE_REGISTER A41
and B76, and campaign Section 1z.19.
