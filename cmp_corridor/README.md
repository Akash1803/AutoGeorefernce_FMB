# CMP corridor parcel layer (Oct 2026)

Scripts used to turn the FMB sheets of the CMP road-proposal corridor (3 roads, 32 km,
Medavakkam to Echur; 17 villages of Chengalpattu district) into one clean parcel GeoJSON
for LPS submission: `CMP_Final_Deliverable.geojson` (EPSG:32644, 5,384 parcels on
2026-10-26) = Akash's part (`CMP_two_Parts`) + the team's part (`MONI_pART`).

Code only. Project data (sheets, vectors, registers, satellite renders, the deliverable)
lives under `D:\Projects\CMP_Roads_LPS\` and is never committed: it carries owner names
and patta numbers from Tamilnilam and Google imagery.

Rules the pipeline follows (Akash, 21-26 Oct 2026):

- one GeoJSON, edited only through QGIS edit sessions on his open layer (never reload or
  replace a layer he has open); every change is a small fix file applied by `qgis/*.py`
- his approved (green) and hand-made parcels never move; georeferenced sheets stay rigid;
  topology is closed by moving the non-manual side only
- no new polygons to hide gaps: a gap is closed by edge matching (one side extended onto
  the other, dividers carried straight across, growth <= 30 %), by sharing a seam strip at
  its midline, or by a PDF face that the vector conversion had lost - never by invention
- parcels are judged on Google Satellite; Puvi is at most a prior, never a target
  (for the Kolathur question Akash ruled: PDF data only)
- a register value is copied only when the key matches exactly; no rate on Reginet means
  every Reginet field stays empty
- run both validators (QGIS internal for duplicate vertices, GEOS for validity) and the
  union-hole / overlap report before showing anything

## Folders

`topology/` - exact coverage and cleaning

- `coverage_fix_gen.py` - node every boundary on a 2 mm grid, polygonize, give each face
  back to its parcel; hairline overlaps to the bigger parcel, hairline gaps split at the
  junctions along them. Env: `CF_FILE`, `CF_PFX` (default = CMP_two_Parts), `CF_MAXGAP`
  (m2, default 0.5) and `CF_WMIN` (m, default 0.05) raise what counts as a sliver.
  `coverage_fix2.py` is the EPSG:4326 version used on `CMP_my_Part`.
- `resolve_overlaps_gen.py` - real overlaps: green > georeferenced sheet > team parcel;
  same rank -> the split that leaves both nearer their Tamilnilam extent; no register ->
  the bigger parcel gives up the overlap. Loser cut only.
- `edge_match_gen.py` - seam gaps (<= 3 m wide) closed by extending one side; env
  `EM_MAXW`, `EM_MAXGROW`, `EM_RECV` (surveys allowed to receive), `EM_GAPS` (gap pieces
  prepared elsewhere), optional window `x0 y0 x1 y1`.
- `half_assign.py` - a seam strip split at its equidistant midline, each half to the side
  it touches, cut perpendicular to the strip at the receiving plots' dividers.
- `despike_two.py`, `tails_two.py`, `clean_nodes.py`, `collinear.py` - spikes (< 8 deg),
  thin hanging tails, 1 cm duplicate vertices on flagged parcels only, 5 mm collinear
  clutter.
- `topo_two.py`, `topo_strict.py`, `gap_report.py`, `gap_pics.py`, `validfix.py` - the
  checks: overlaps, union holes with width/corridor share/touching parcels, satellite
  pictures of selected gaps (`GP_ADDS` draws proposed fills in cyan).

`placement/` - rigid sheet placement

- `georef_survey.py` (sheet plots + rigid fit), `close_gap.py` (ICP onto a neighbour),
  `fit_on_team.py` (pose on the team footprint, north-up within +-20 deg),
  `rebuild_batch.py` / `replace_in_merged.py` (team shapes restored, sheets re-inserted
  in order), `chain_v2.py` / `prior.py` / `gap_close.py` (Ponmar chain: Puvi polygon as
  the sheet-shape prior, satellite as the judge), `icp_lib.py`, `georef_lib.py`.

`annotation_qc/` - dual check of every plot number

- `qc_tiles.py` - PDF tiles with clean copies and `plots.csv` for two independent readers
  (`reads/A_groupN.csv`, `B_groupN.csv`).
- `ann_check.py` - machine read of the printed numbers (glyph reader), rigid fit then a
  label-anchored thin-plate-spline warp of the sheet onto the layer; verdict per plot.
- `tile_poly_map.py`, `recon3.py` - three-way reconciliation (machine + 2 readers):
  OK / FIX / FIX (readers) / CHECK / NOT READ.
- `settle_crops.py`, `settle_crops2.py` - 6-panel PDF crops with the layer outline
  inverse-mapped, for the CHECK rows; `build_labels.py`, `patch_labels.py` turn the
  settled decisions into a fix file and refuse duplicate numbered plots.

`missing_plots/` - plots the vector conversion lost

- `fill_kind2.py` - a warped PDF plot clipped to the ONE host parcel holding >= 85 % of it.
- `fill_holes.py` - a union hole cut along the warped PDF plot boundaries, each piece a
  parcel of that survey, numbered only when the printed number is certain, in the
  register and not yet in the layer; `FH_CLIP` keeps to the corridor buffer. Cut existing
  parcels out of every fill first - a hole can hold island parcels.
- `add_roads.py` - unnumbered PDF faces on asphalt (satellite test + visual pick).
- `georef124.py`, `fill124.py` - local rigid poses for the lumps of a sheet that would not
  fit as a whole (Kattour 124), relabel only.
- `header_areas.py` - printed survey area of every sheet (title block) for the area check.

`registers_merge/`

- `register_update.py` - Tamilnilam + Reginet columns with the master rules
  (`REG_FILE`, `REG_OUT`).
- `merge_final.py` - the two parts into one schema (`dataset` says which), seam hairlines
  cut from the team side.
- `area_check_excel.py` - survey-wise Excel: PDF (title block) area vs layer area.
- `apply_offline.py` - the same fix-file format applied to a GeoJSON outside QGIS.

`qgis/` - helpers run inside QGIS through the qgis_mcp socket (`q.py`, localhost:9876):
apply a fix file (`change` / `add` / `delete`, matched by village + kide + a point),
attribute-only applies, validators, the final layer loader with its rule style.

## Fix file

```json
{"change": [{"vc": "35_11_047", "old_kide": "905/1B", "old_pt": [x, y], "wkt": "...", "attrs": {...}}],
 "add":    [{"wkt": "...", "attrs": {...}}],
 "delete": [{"vc": "...", "old_kide": "...", "old_pt": [x, y]}]}
```

Write the fix file immediately before applying it and delete it afterwards: on 2026-10-26
a crashed script left a stale file and the chain applied an old fix to 105 parcels (restored
from the backup, verified identical). Never mask a script's exit code with `| tail`.
