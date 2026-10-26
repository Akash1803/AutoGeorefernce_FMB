# Evaluation rows: schema and how to run (PR 0)

Every run of the pipeline that places parcels can append one row per placed parcel to
`<project>\_logs\eval\rows.csv`. Rows are the input to every report and, if the gate in
PR 3 is met, to the confidence model. Nothing in a row comes from pixels except scalar
features and the imagery window's hash.

## Producing rows

```
python -m autogeoref.georef_village 35_04_077 --topology --rows      # production run, mode "full"
python -c "from autogeoref import evaluate; evaluate.leave_one_out('35_04_077')"     # mode "loo"
python -c "from autogeoref import evaluate; evaluate.seed_run('35_04_077', ['48A'], rings=2)"   # mode "seed"
python -m autogeoref.georef_village 35_04_077 --report                # report from all rows so far
python -m autogeoref.georef_village 35_04_077 --purge-cache           # delete the imagery cache
```

Configuration: `configs\default.json`, or `AUTOGEOREF_CONFIG=<path>`, or `--config <path>`.
Unknown keys are an error.

## Files under `_logs\eval\`

| File | Content |
|---|---|
| `rows.csv` | one row per parcel per run, appended, header once |
| `references.csv` | every hand-placed parcel with its two references and stretch measures |
| `excluded_stretched.csv` | parcels whose hand placement fails the stretch rule, with the reason |
| `excluded_disputed.csv` | parcels excluded for any other reason (Puvi identity conflict, railway conflict, the village's `disputed.csv`) |
| `report_<id>.md` | the comparison table and worst ten, boilerplate first |

## Row schema

Identity: `run_id`, `run_at`, `tool_git_sha`, `village`, `survey`, `mode` (loo, seed, full),
`seeds`, `pass_no`, `imagery_source`, `window_sha256`, `backfilled`.

Placement features: `method`, `n_candidates`, `support_m`, `support_next_m`, `runner_up_ratio`,
`boundary_rms_m`, `n_neighbours`, `anchor_partners`, `side_ok`, `side_bad`, `contradictions`,
`printed_far_n`, `share`, `margin`, `observable`, `ambiguous`, `in_window`, `sigma_pos_m`,
`sigma_head_deg`, `heading_deg`, `pose_tx`, `pose_ty`, `topo_max_move_m`, `topo_clipped`.

Rule-based verdict: `colour`, `confidence_rule`.

Reference: `ref_source` (points, geometry, none), `ref_rms_m`, `ref_area_scale`,
`ref_stretch_long_pct`, `ref_stretch_short_pct`, `ref_disagree_mean_m`, `ref_disagree_p95_m`,
`ref_disagree_corner_mean_m`, `stretch_band`, `stretched`, `stretched_reason`, `disputed`,
`disputed_reason`.

Errors: `err_fmb_m` (centroid, against the FMB-exact reference), `err_hand_m` (centroid, against
the saved hand geometry), `heading_err_deg`, `max_vertex_err_m` (Hausdorff, FMB-exact).

Label: `within_3m` is 1 or 0 when the parcel has a usable reference, empty when it is stretched,
disputed or has no hand placement. `acceptance_m` in the config sets the limit.

## Stretch rules (config `stretch.rule`)

- `pct`: area scale `sqrt(hand area / sheet area)` off 1 by more than `pct` percent. Default.
- `axis`: either side of the minimum rotated rectangle off by more than `pct` percent.
- `displacement`: the larger of the outline-sample mean and the corner mean of the distance
  between the FMB-exact and the hand placement above `displacement_mean_m`.

A rigid fit worse than `rms_m` excludes under every rule. Parcels between `band_pct` and `pct`
stay labelled and are reported separately (`stretch_band`). A survey listed in the village's
`FMB_Vector\<village>\disputed.csv` (`survey,reason`) is excluded whatever the rule.

## Split

`evalrows.folds(rows)` is leave-one-village-out over villages with at least one labelled row.
`evalrows.split_by_village(rows, village)` is the only split function; there is no parcel-level
split. Rows of one village never appear on both sides.

## Imagery

Windows are written to `<project>\_cache\imagery\<village>\` (gitignored) and logged in
`_cache\imagery\imagery_log.csv` with source, zoom, bounds, sha256, time and run id. The pin file
`FMB_Georef\<village>\raster.json` points into the cache. `imagery.source` is `google_xyz`
(inference only) or `geotiff:<path>`.

## Changes of 2026-09-21

- `append_rows` never writes a (run_id, village, survey) that is already in the file, and migrates
  the header when the schema grows (`migrate_rows`).
- New column `sheet_qc_reason`: a sheet the converter got wrong is logged (`mode = sheet_qc`,
  `qc_failure_row`) with the reason, no features and no label, and stays out of every figure
  (`metrics`, `worst10`, `distinct_summary`). Today: 613 (drawn 9.9 m short) and 47A (offset
  arrows cut plot 3).
- The report and `--report` print distinct labelled parcels with per-class counts (latest row per
  parcel) next to the row count; the PR 3 gate is 300 distinct labelled parcels across two
  villages (`report.pr3_gate`), never rows. On 2026-09-21: 14 distinct labelled parcels
  (9 within, 5 beyond) across two villages; gate not met.
