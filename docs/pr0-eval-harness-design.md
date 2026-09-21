# PR 0 design: evaluation harness and row logging

Status: for approval, 2026-09-21. No ML in this PR. The rule-based pipeline is unchanged and
its results must come out byte-identical (checked by diffing `georef_status.csv` before and
after). Depends on [ml-step3-data-plan.md](ml-step3-data-plan.md) for what a label is.

## Files

| File | Change | Responsibility |
|---|---|---|
| `autogeoref/references.py` | new | Reference A and B per hand-placed parcel; stretched and disputed flags with reasons; writes the two exclusion CSVs |
| `autogeoref/evalrows.py` | new | `FeatureRow` schema; rows from a run's status; append, load, split by village, folds |
| `autogeoref/report.py` | new | comparison table, worst-10 with causes, boilerplate; markdown to `_logs\eval\` |
| `autogeoref/config.py` | new | `Config` from `configs\default.json`; `imagery.source`, `imagery.cache_dir`, `imagery.zoom`, `ml.confidence = rule` |
| `autogeoref/raster.py` | modify | export into the cache dir; append `imagery_log.csv`; `purge_cache()` |
| `autogeoref/evaluate.py` | modify | `leave_one_out` and `seed_run` write rows after each run, errors against both references |
| `autogeoref/georef_village.py` | modify | `--config`, `--rows`, `--purge-cache` |
| `configs/default.json` | new | the defaults above |
| `requirements.txt` | new | pinned to what is installed here: geopandas 1.1.4, shapely 2.1.2, rasterio 1.5.0, scipy 1.18.0, numpy 2.5.1, scikit-learn 1.9.0, opencv-python 5.0.0, pyproj 3.7.2, openpyxl (version to read at commit time) |
| `.gitignore` | modify | `_cache/` |
| `tests/test_references.py`, `tests/test_evalrows.py`, `tests/test_report.py`, `tests/test_config.py` | new | |
| `docs/eval-rows.md` | new | schema and how to run |

## Row schema: `_logs\eval\rows.csv`, one row per parcel per run

Identity: `run_id`, `run_at`, `tool_git_sha`, `village`, `survey`, `mode` (loo, seed, full),
`seeds`, `pass`, `imagery_source`, `window_sha256`, `backfilled`.

Placement features (from the status row): `method`, `n_candidates`, `support_m`,
`support_next_m`, `runner_up_ratio`, `boundary_rms_m`, `n_neighbours`, `anchor_partners`,
`side_ok`, `side_bad`, `contradictions`, `printed_far_n`, `share`, `margin`, `observable`,
`ambiguous`, `in_window`, `sigma_pos_m`, `sigma_head_deg`, `heading_deg`, `pose_tx`, `pose_ty`,
`topo_max_move_m`, `topo_clipped`.

Rule-based verdict: `colour`, `confidence_rule`.

Reference: `ref_source` (points, geometry, none), `ref_rms_m`, `ref_scale`, `ref_aniso`,
`stretched`, `stretched_reason`, `disputed`, `disputed_reason`.

Errors: `err_fmb_m`, `err_hand_m`, `heading_err_deg`, `max_vertex_err_m`.

Label: `within_3m` (1, 0, or empty when excluded or unreferenced).

## Split logic

`evalrows.folds(rows)` yields, for every village with at least one labelled row, the pair
(train = all other villages, test = that village). `split_by_village(rows, test_village)` is the
only split function; a parcel-level split is deliberately not provided. Stretched and disputed
rows keep their features and errors and carry no label.

## Report: `_logs\eval\report_<run_id>.md`

1. Boilerplate sentence from Step 3 section 6.
2. Comparison table, one line per fold and per method (rule-based now; a candidate column
   stays empty until PR 3): rows, labelled rows, share of green that is within 3 m, share of
   within-3 m parcels ranked green or amber, share of red that is beyond 3 m, median and 90th
   percentile of `err_fmb_m`, same for `err_hand_m`.
3. Worst 10 by `err_fmb_m`, with survey, run, both errors, colour and the cause read from the
   status notes: contradiction, one neighbour, tie, outside window, stretched reference,
   no transcription, clipped by railway.

## Runs that feed it

- Backfill: the existing three leave-one-out CSVs and two seed runs on Kizhikaranai are read once
  by `evalrows.backfill()`; fields the old CSVs lack stay empty and the rows carry `backfilled = 1`.
- Fresh: one leave-one-out and the two seed rings on Kizhikaranai with the current code.
- Kolathur rows arrive after PR 1 gives it a transcription.

## Imagery handling in this PR

`raster.export` writes to `config.imagery.cache_dir` and appends a line to `imagery_log.csv`;
the pinned file for each village moves there; `georef_village --purge-cache` deletes the
directory. Only scalars and `window_sha256` reach the rows.

## Acceptance

- Fast test suite green, new tests included.
- `rows.csv`, both exclusion CSVs and one report produced from a leave-one-out run.
- `georef_status.csv` of a village run is identical before and after the PR.
- Runtime of the harness on this laptop reported in the PR description (unmeasured until run).

## Measured and unmeasured

Measured here: topology re-run 3 to 4 s; seed ring 2: 928 s; leave-one-out about 40 s per
parcel. Unmeasured: harness runtime, row counts, whether the PR 3 gate can be met.
