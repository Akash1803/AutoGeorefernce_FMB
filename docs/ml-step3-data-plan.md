# Step 3: data and labelling plan

Status: for approval, 2026-09-21. No code in this step. Follows the approved
[ml-hybrid-plan.md](ml-hybrid-plan.md) and the constraints of the 2026-09-21 brief.

## 0. Constraints this plan obeys

- Laptop CPU only; nothing is trained in this phase; no PyTorch. Labels are collected now so that
  PR 3 (scikit-learn confidence model) can use them if the 300-row gate is met.
- Google Satellite tiles are read for inference only. No training dataset is built from tiles;
  nothing derived from pixels is stored except scalar features and a window hash.
- There are no GNSS or independent check points. Every score is against hand placements and is
  capped by their mutual consistency of 1 to 3 m. This sentence goes into every report.
- Acceptance: the FMB-exact placement lands within 3 m of the hand placement.

## 1. What exists (checked on this laptop, 2026-09-21)

| Item | Count | Where |
|---|---|---|
| FMB PDFs in the 60 m corridor | 656, 19 villages | `FMB_Sketches\<village>\` |
| Converted sheets (sheet metres) | 656 | `FMB_Vector\<village>\<survey>_parcels.geojson` |
| Hand-placed parcels, Kizhikaranai 35_04_077 | 15 | `FMB_Vector\35_04_077\*_parcels_modified*.gpkg` |
| Hand-placed parcels, Kolathur 35_04_073 | 11 (count from the 2026-09-20 audit, not re-opened today: unverified) | `FMB_Vector\35_04_073\` |
| Team `.points` files | 26 | beside the sheets |
| Two-reader neighbour transcription | Kizhikaranai only | `neighbour_transcription\nb_override_35_04_077.json` |
| Run records | 3 leave-one-out CSVs, 2 seed-ring runs, all Kizhikaranai | `_logs\` |
| Independent check points | none | |

Puvi is reference only and never a label source.

## 2. What a label is and where it comes from

The unit is one survey sheet in one run: all its plots together, one row.

Two references per hand-placed parcel, both kept on every row:

- **Reference A, FMB-exact.** The rigid pose of the sheet fitted to the team's `.points`
  (`ref_source = points`), or to the hand geometry when that fit is worse than 3 m rms
  (`ref_source = geometry`). This is `evaluate.truth_pose` today. The sheet at this pose keeps
  every printed length.
- **Reference B, as saved.** The hand-placed geometry exactly as the team saved it (affine,
  vertex-edited). Newest file per survey; older versions are listed as superseded and unused.

Measures per placed parcel: centroid distance to A (`err_fmb_m`) and to B (`err_hand_m`),
heading difference to A (`heading_err_deg`), largest vertex distance to A (`max_vertex_err_m`),
plus the reference's own rms, scale and anisotropy.

**Label** (for PR 3): `within_3m = 1 if err_fmb_m <= 3.0 else 0`, per decision 1. Null when the
parcel is excluded (section 3) or has no reference (a parcel the team has not placed).

Why centroid plus heading: the tool moves a rigid body, so those two numbers describe the whole
error; `max_vertex_err_m` is kept so a long strip rotated by a degree is not hidden by a good
centroid.

## 3. Exclusions

**Stretched** (decision 2): the hand placement changed the sheet's dimensions. Proposed rule,
any of:

- area scale `sqrt(hand_area / sheet_area)` outside 0.90 to 1.10;
- similarity-fit anisotropy above 1.10;
- any outline edge of the hand geometry differing from the printed FMB length by more than 10 %;
- rigid fit of sheet to hand geometry worse than 3 m rms.

Effect on the 15 Kizhikaranai parcels with the area scales recovered on 2026-09-19: at plus or
minus 5 % five would be excluded (40B 0.934, 43B 0.917, 48B 0.902, 170 0.950, 47B 1.063); at
plus or minus 10 % none by area, and 47B by the edge rule (124 m drawn where the sheet prints
109 m). 169's `.points` affine is strongly anisotropic (1.21 by 0.59) and would be excluded by
the anisotropy rule. **Please confirm 10 %** or set another figure. Stretched parcels stay in
the rows with `stretched = 1`, `within_3m` null, and both errors reported; they are listed in
`_logs\eval\excluded_stretched.csv` with the triggering rule.

**Disputed**: the 22 surveys where Puvi and the portal describe different parcels (569B and its
kind), sheets marked No Sketch, parcels in `rail_conflicts.csv`, and any survey the team flags
in the tracker. Listed in `_logs\eval\excluded_disputed.csv` with the reason; never labelled.

## 4. Leakage-safe split

- **The unit of splitting is the village.** Rows of one village are never divided between train
  and test: support, neighbours and sides are computed from the same anchors, and rows from
  different runs of the same sheet share geometry. A sheet is atomic inside its village.
- With two labelled villages (Kizhikaranai now, Kolathur after PR 1) the split is
  leave-one-village-out, two folds, each reported on its own. No pooled number.
- Runs multiply rows (leave-one-out, seed rings, full runs); all rows of a village travel with it.
- Isotonic calibration is fitted on training folds only.
- PR 3 gate: 300 or more labelled rows across at least two villages. Today there is one village;
  the harness will count the rows (expected of the order of 75 from existing runs: unmeasured
  until PR 0 counts them). If the gate is not met, PR 3 stops and says so.

## 5. What the imagery constraint does to storage

- Tiles are read through the GDAL descriptor on the fly. Window GeoTIFFs go only to
  `<project>\_cache\imagery\<village>\`, gitignored and deletable with one command; each window
  is logged in `imagery_log.csv` (source, zoom, bounds, sha256, fetched at, run id).
- The current exports under `FMB_Georef\<village>\satellite_*.tif` move to that cache in PR 0
  and the pin file is updated. `*.tif` stays gitignored.
- Rows store scalars only (share, margin, observable, in_window) and the window's sha256.
  No crops, masks or tiles are kept as data. PR 2's SAM overlays for review are written to the
  cache and die with it.
- `imagery.source` in the config: `google_xyz` (default, inference only) or `geotiff:<path>` for
  licensed orthos. Same cache and log for both.

## 6. Boilerplate for every report

"Errors are measured against the team's hand placements, which disagree with each other by
1 to 3 m. Results are therefore capped by anchor consistency and are not absolute accuracy.
Stretched parcels are scored against both the FMB-exact and the saved placement and carry no
label."

## 7. Measured and unmeasured

Measured on this laptop: topology re-run 3 to 4 s per village; seed ring 2 run 928 s
(2026-09-20); leave-one-out about 40 s per parcel.

Unmeasured or unverified: the Kolathur parcel count and quality; the number of labelled rows
available for PR 3; everything about SAM and PaddleOCR timing and accuracy.
