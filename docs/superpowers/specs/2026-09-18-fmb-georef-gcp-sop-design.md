# FMB georeferencing SOP v2: frozen anchors, satellite evidence, GCP files — design

Date: 2026-09-18, revision 2 (after the trials in §3 and a 40-finding adversarial review, 21 confirmed).
Status: for the user's review. Supersedes the placement engine of the 2026-09-17 design; the measurement
attributes (§11 there), the sheet-reading pipeline, the neighbour transcription and the topology fix are
reused. Project: Tambaram – Chengalpattu rail corridor, `D:\Projects\Tambaram_Chengalpattu`. Code home:
`D:\code\FMB_to_GeoJSON\autogeoref\` (the 2026-09-17 scripts are in `autogeoref\scratch_2026-09-17\`).

## 1. Why v2

The 2026-09-17 engine placed sheets by matching shared boundary lengths between neighbours. Measured today
against the six parcels the team hand-placed this morning (46A, 43A, 43B, 42A, 42B, 169), its placements
were 1–7 m off for ordinary parcels and 14–24 m off for the railway-side strips 42B and 169, and its block
adjustment moved the team's earlier parcels by 2.7–7.8 m. The user rejected that and set two rules:

1. Manually georeferenced parcels are never touched. Automation only places their neighbours.
2. Position evidence comes from a georeferenced satellite raster with the Puvi vector as the location
   reference, from which ground control points (GCPs) are derived automatically; the team corrects GCPs
   only where the tool is unsure. Approach approved on 2026-09-18: **B + C** (automatic GCPs, team GCPs
   as fallback).

Two facts found during the trials shape the design:

* **The team already works with GCPs.** Every manual parcel has a QGIS Georeferencer file
  `<survey>_parcels.geojson.points` beside its sheet (3–4 GCPs, sheet metres → EPSG:32644). The team's
  saved polygons are affine/polynomial transforms of those GCPs: implied scales run from 0.83 to 1.15 and,
  for the railway strip 169, 1.21 × 0.59. A scale-1 fit of the team's own GCPs differs from their saved
  polygon by 1–9 m (169: 40 m). So the `.points` file is the natural channel for auto GCPs and team
  corrections, and "frozen anchor" needs a definition (§2) because the frozen files are not rigid copies
  of the sheets.
* **Imagery fixes rotation and cross-track offset, not along-track position.** 92 % of the outline length
  that coincides with image edges lies on edges parallel to the corridor; the transverse boundaries
  between strips are visible for 16 % of their length and for 0 % on five of seven parcels. Along-track
  position must come from frozen neighbours, the printed neighbour numbers and shared boundaries.

## 2. Decisions

| Question | Decision |
|---|---|
| Rigidity | Rotation + shift only, scale 1, for everything the tool writes. Every FMB length and angle survives placement. |
| Frozen anchor | A survey with a manual file. The file on disk is never modified, moved, re-fitted, topology-edited or overwritten. For constraints the tool uses the anchor's **pose**: the scale-1 fit of the sheet to the team's `.points` GCPs when that fit's rms ≤ 3 m, else the scale-1 fit to the manual geometry (`manual_poses.py`). `anchors.csv` records file, pose, GCP rms, implied scale and anisotropy so a stretched anchor is visible. |
| Which manual file | Any `<survey>_parcels_modified*.gpkg` (and any hand GeoJSON the team names) is a manual file. One anchor per survey: the newest by modification time unless `anchors.csv` names another. Other variants are listed as superseded and never loaded. A file is manual iff its geometry fingerprint is not recorded as tool-written in `georef_status.csv`; files with no record, and every v1 `_AUTOv2` / `_AUTODEMO` output, count as manual. |
| Evidence order | 1. shared boundaries with frozen anchors (point-to-line, σ = max(0.30 m, anchor rms)); 2. shared boundaries with placed auto parcels and the printed neighbour numbers; 3. satellite image edges (refine ≤ 1.5 m / ≤ 2°, or decide alone only when unambiguous and observable in two directions); 4. Puvi shape and position (window and start pose only). |
| Puvi | Location reference only. Never a geometry source, never a gate. Exact survey key first (47B), numeric parent second. **Area gate:** when the Puvi polygon's area differs from the sheet's by more than a factor of 3 the two sources describe different parcels (22 cases corridor-wide, section 14) - Puvi is dropped for that survey, `puvi_trusted` is false, and the parcel takes the neighbour-window path. |
| Satellite raster | The QGIS "Google Satellite" XYZ layer exported per village buffer to a GeoTIFF in EPSG:32644 at zoom 20 (0.15 m/px) through GDAL's tile driver; pinned in `raster.json` (path, sha256, fetch date); `--raster` archives the previous file. No screenshots, no ECW. |
| GCPs | Auto GCPs are written in the team's own format, `<survey>_parcels.geojson.points`, only when no team file exists. The team opens the sheet in the QGIS Georeferencer with those GCPs pre-loaded, moves or adds points, saves; `--refit` reads any `.points` newer than the tool's record as team GCPs and fits rigidly. Team GCPs persist in `gcp\<survey>_gcp.csv` with `source = team`; no run ever rewrites them. |
| Output naming | Auto parcels: `FMB_Vector\<village>\<survey>_parcels_modified.gpkg`, with `georef_method` (anchors / image / team-gcp / puvi-only) and a colour, never written while a `_modified<N>` variant exists for that survey. The tool deletes nothing: a superseded auto file is moved to `_logs\georef_superseded_<run>\`. |
| Railway land | Overlaps between railway-land sheets (169, 170, 171 in Kizhikaranai) and parcel sheets are listed, never clipped automatically. |

## 3. Evidence (2026-09-18, `autogeoref\scratch_2026-09-17\spike_gcp\`)

**Raster export.** `rasterio` opened a GDAL_WMS/TMS description of the Google tile URL and wrote a
366 × 540 m window at zoom 20 to a GeoTIFF in EPSG:32644 in under a minute (2441 × 3601 px, 0.15 m/px).

**Edge coverage.** OpenCV's line segment detector found 2643 segments ≥ 4 m. 53 % of the team's seven
outlines (Kizhikaranai) lie within 1.5 m of a same-direction segment (±12°): 48A 81 %, 171 79 %, 47B 47 %,
170 46 %, 47A 43 %, 48B 43 %, 46B 12 %. Of the matched length, 92 % is corridor-parallel; transverse
boundaries are matched 16 % overall and 0 % on 170, 171, 46B, 47A, 48B. Along-track score profiles are flat
over ±15 m for those five parcels.

**Pose recovery, leave-one-out** (parcel hidden, pose searched, compared with the team's rigid pose):

| Parcel | Puvi centroid off | Puvi shape alone | Image alone, free search | Puvi shape, then image ±10 m/±10° |
|---|---|---|---|---|
| 48A | 1.6 m | 0.3° / 4.4 m | 0.7° / 0.9 m ✔ | 0.7° / 1.1 m ✔ |
| 47B | 3.6 m | 3.8° / 3.6 m | 0.3° / 2.6 m ✔ | 0.3° / 2.6 m ✔ |
| 47A | 5.6 m | 3.8° / 5.6 m | flipped 175° | 0.3° / 8.7 m slid |
| 48B | 4.5 m | 1.5° / 6.9 m | flipped 173° | 2.2° / 8.5 m slid |
| 46B | 8.1 m | 1.8° / 8.1 m | slid 28 m | 1.3° / 18.7 m slid |
| 171 (railway) | 11.5 m | flipped 178° | flipped | flipped |
| 170 (railway) | 36.8 m | flipped 180° | flipped | flipped |

**The team's GCP files.** Rigid (scale-1) fit rms of the team's own GCPs: 46A 0.19 m, 43A 0.77, 170 1.27,
46B 1.31, 171 2.11, 42B 2.21, 48B 3.31, 47A 3.42, 48A 3.52, 43B 5.33, 42A 5.65, 47B 6.41, 169 29.4 m.
Implied affine scales 0.83–1.15 (169: 1.21 × 0.59).

**What this fixes.** Imagery decides alone only for parcels with distinctive walls (48A, 47B) and, even
there, only cross-track and rotation are well observed. Near-symmetric strips are oriented only by the
printed neighbour numbers and frozen neighbours. Anchor constraints must be tolerant of the anchors' own
1–6 m non-rigidity. Thresholds must be calibrated on the production score, not on the spike's binary
share (which is 0.14 at the spike's own tolerance under the Gaussian score).

## 4. Scope

In scope: one village at a time; every buffer parcel with a sheet and no manual file; raster export;
start pose; boundary evidence; image alignment; observability and margin tests; rigid fit; block
adjustment of auto parcels with anchors fixed; `.points` output and refit; review group; topology between
auto parcels; status CSV; tracker columns; a leave-one-out test harness on a copy of the village folder.
Out of scope: changing any manual file or `.points` file the team wrote; ECW; other imagery; deciding
railway-land versus parcel conflicts; villages without sheets. 152 of the 655 corridor sheets have no
Puvi polygon: they take the neighbour-window path from the start (§5 step 4).

## 5. SOP (what runs, in order, per village)

1. **Inputs**: `FMB_Vector\<village>\<survey>_parcels.geojson` (sheet metres); the Puvi village
   shapefile; `Railway_Buffer_Vector_Plots.gpkg`; manual files and their `.points`; the QGIS satellite
   layer; `FMB_Vector\<village>\neighbour_transcription\nb_override_<village>.json` (two-reader
   transcription of the printed neighbour numbers, made once per village).
2. **Freeze anchors**: choose one manual file per survey (§2), fingerprint it through an OGR read,
   recover its pose from `.points` or geometry, write `anchors.csv` (file, pose, gcp_rms, scale,
   anisotropy, superseded variants). Report anchor pairs whose shared boundary disagrees by more than
   0.30 m + both rms (48A–171: 10.7 m) — both stay.
3. **Reference raster**: `FMB_Georef\<village>\satellite_z20_<date>.tif`, buffer extent + 100 m,
   overviews, nodata outside the buffer, pinned in `raster.json`.
4. **Search window and start poses** per auto parcel: Puvi polygon (exact key, then numeric parent)
   buffered 40 m; start pose by symmetric-difference alignment, keeping the 180° twin when it scores
   within 15 %; the printed neighbour sides choose between twins. No Puvi polygon → window = union of the
   placed neighbours' buffers, start pose from any shared boundary. Where the parcel touches a frozen
   anchor, the anchor-implied pose is a further hypothesis and the search is seeded from it.
5. **Boundary evidence**: the 2026-09-17 chain matcher against the anchors' rigid sheets at their
   recovered poses (never against hand vertices) and against placed auto parcels; only between sheets that
   print each other's number; label-side test; ≥ 3 vertex pairs or two anchors before any single edge.
   Observations are point-to-line (auto boundary sample to the anchor's boundary line), σ = max(0.30 m,
   anchor rms); vertex pairs only where both sheets have a real corner (turn ≥ 20°).
6. **Image alignment**: segments ≥ 3 m with direction inside the window; score = share of 1 m outline
   samples within 1.5 m of a same-direction segment (±12°), the spike's metric; search ±10 m / ±10°
   around every hypothesis (coarse 1 m / 1°, fine 0.25 m / 0.25°). **Margin test**: the best pose must
   beat the best pose more than 3 m or 3° away by ≥ 0.15 and every local maximum within 75 % of the best
   must lie within 1 m of it; otherwise the image result is "ambiguous" and only verifies.
   **Observability**: ≥ 40 m of outline matched on edges of one direction (rotation, cross-track) and
   ≥ 20 m matched on edges > 30° from it, or an anchor/auto boundary that fixes the along-edge position;
   without the second the parcel is at best amber "along-track from neighbours only".
7. **Decision**: with an anchor boundary the reference pose stands and imagery may move it ≤ 1.5 m / 2°
   when unambiguous; without an anchor an unambiguous, observable image pose decides; otherwise the
   parcel is written at the start pose as puvi-only (red). Corner GCPs are then derived from the decided
   pose (FMB corners with turn ≥ 20°, map position = posed corner, `matched` flag from the image) and
   written as `<survey>_parcels.geojson.points` when no team file exists, and to `gcp\<survey>_gcp.csv`.
8. **Block adjustment** of auto parcels only, anchors and team-gcp parcels fixed: auto–auto boundaries
   σ 0.30 m, anchor boundaries σ as in step 5, image evidence σ 1.0 m (unambiguous only), start pose
   σ 10 m; red parcels' boundaries enter at their prior σ (2.5 m) so they cannot drag greens. Trusted-first
   chain admission and demotion as in `block_adjust.py`, solved in a frame centred on the block.
9. **Confidence and colour**, computed after the adjustment: green = every anchor/auto boundary within
   0.30 m + anchor rms, observable in both directions, image share ≥ green threshold; amber = one
   condition short, or moved > 0.30 m by the adjustment since the team accepted it ("moved by
   adjustment"); red = puvi-only, ambiguous image, or Puvi distance > 40 m (reported, never used to move).
   Thresholds are calibrated before coding on the 20 manual parcels (7 + 13) with the production score
   and recorded in §3.
10. **Outputs**: `<survey>_parcels_modified.gpkg` (parcels + edges layers, §11 attributes; fingerprint
    taken after the last write of the run, topology included, stored as `fp_placed` and `fp_final`),
    `.points` files, `gcp\<survey>_gcp.csv`, `georef_status.csv` (colour, method, share, margin,
    observability, boundary residuals, Puvi distance, notes, alternative pose), read-only display layer
    `gcp_review.gpkg` (one point per corner, `survey, corner_id, sheet_x, sheet_y, auto_x, auto_y,
    matched, source`), QGIS group "Georef review – <village>" with the raster, Puvi outline (grey dashed),
    anchors (black), auto parcels by colour, GCP points, and dotted outlines for twin/alternative poses.
11. **Team fallback (C)**: for red or amber parcels the team opens the sheet in the Georeferencer with the
    tool's `.points`, moves or adds GCPs (two suffice for a rigid pose), saves. `--refit` reads `.points`
    files newer than the tool's record, fits rigidly (team GCPs σ 0.20 m, no ≥ 3-GCP rule, anchor
    boundaries reported but not enforced), writes `georef_method = team-gcp`, copies the rows to the CSV
    with `source = team`, and holds those parcels fixed in every later adjustment until `--reset-gcp`.
    A parcel whose geometry the team edits directly is treated as a team pose when the sheet fits it at
    scale 1 ± 0.2 % with rms < 0.05 m, else it becomes a frozen anchor by the §2 rule.
12. **Topology** between auto parcels (snap 0.30 m, clip, fill) with anchors and team-gcp parcels fixed;
    railway-land conflicts listed in `rail_conflicts.csv`; fingerprints re-recorded afterwards.
13. **Tracker and worklist**: the team's Status column is never written; tool-owned columns Auto colour,
    Auto note, Auto run are added to the Tracker sheet keyed by village code + survey number, Status moves
    only from Not Started to In Progress for red/amber. When the workbook is open the rows go to
    `tracker_pending.csv` and are applied at the next run. `worklist.csv` (village, survey, colour, note)
    is written for the team to open in QGIS.

## 6. Components (`autogeoref\`)

| Module | Job | Interface |
|---|---|---|
| `anchors.py` | choose, fingerprint and pose manual files; `.points` reader | `load_anchors(village) -> {survey: Anchor}`; `read_points(path)`; `is_tool_written(path)` |
| `raster.py` | export and pin the satellite GeoTIFF | `export(village, extent, zoom=20) -> path` (rasterio + GDAL_WMS) |
| `window.py` | Puvi window and start poses with twin | `start_poses(survey) -> [(theta, t)], window` |
| `neighbours.py` | adjacency and sides from the transcription | unchanged from 2026-09-17 |
| `match.py` | chains, label-side test, corroboration, point-to-line observations | from `autogeoref_demo.py` |
| `edges.py` | line segments with direction | `segments(raster, window) -> [Segment]` (OpenCV LSD, ≥ 3 m) |
| `align.py` | constrained search, share score, margin and observability tests | `align(outline, segments, hypotheses) -> AlignResult` |
| `gcp.py` | corner GCPs, `.points` writer/reader, CSV, refit | `derive`, `write_points`, `refit(village)` |
| `fit.py` | rigid fit; block adjustment with fixed anchors | from `block_adjust.py` |
| `files.py` | held-layer detection over the socket, temp-write-and-rename, WAL cleanup, superseded moves | `replace(path, writer)` |
| `review.py` | QGIS group (with `--project` guard, non-fatal when the socket is down), status CSV, tracker, worklist | `build_group`, `write_status`, `update_tracker` |
| `topology.py` | existing fix, anchors excluded, railway conflicts listed, fingerprints after | `fix(village, movable=auto_only)` |
| `georef_village.py` | command | `python georef_village.py <village> [--root DIR] [--raster] [--refit] [--reset-gcp SV] [--topology] [--review] [--project QGZ]` |

## 7. Algorithms

**Start pose from Puvi.** θ in 0…359° step 2°: rotate the outline about its centroid, shift to the Puvi
centroid, score = symmetric-difference area; refine ±3° / 0.5°, ±6 m / 2 m. Twin kept when within 15 %.
Sides from the transcription (N/E/S/W of each printed number) choose between twins.

**Share score.** Outline sampled every 1 m; a sample counts when a segment within 1.5 m has direction
within ±12°. Search grid 1 m / 1° over ±10 m / ±10°, then 0.25 m / 0.25°. Margin and observability as in
§5 step 6; the along-track σ of a one-direction match is reported.

**Point-to-line anchor observations.** For each auto boundary sample within 2 m of an anchor's boundary
line (anchor sheet at its recovered pose), the perpendicular distance to that line, σ = max(0.30 m,
anchor rms). Vertex pairs only at real corners on both sheets. Green requires the mean perpendicular offset
≤ 0.30 m + anchor rms.

**Rigid fit and block adjustment.** Unknowns θ, tx, ty per auto parcel; anchors and team-gcp parcels
fixed; Huber loss; centred frame; trusted-first admission (multi-vertex chains first, demotion above 1.0 m,
single edges by trial); per-parcel covariance from the normal matrix reported as position and heading σ.

**Refit from team GCPs.** Rigid (scale 1) least squares on the team's GCP pairs (sheet → map), σ 0.20 m;
residual per GCP written back to the `.points` row; two GCPs suffice; the anchor boundaries are reported
as offsets only.

## 8. Data flow

sheets + transcription → anchors chosen and posed → raster exported and pinned → per auto parcel: window
and start poses → boundary evidence → image alignment with margin and observability → decision → corner
GCPs → block adjustment (auto only) → files, `.points`, status, review group → team edits `.points` in the
Georeferencer for red/amber → `--refit` → topology (auto only) → fingerprints → tracker columns and
worklist.

## 9. Error handling

| Situation | Behaviour |
|---|---|
| Manual file edited since last run | Fingerprint differs → manual by the §2 rule (team pose if scale-1 fit is exact, else anchor); nothing deleted. |
| Several `_modified*` variants | Newest by mtime is the anchor unless `anchors.csv` says otherwise; the rest superseded; no auto file written to that name. |
| Two anchors disagree on a shared boundary | Reported; neither moves; auto parcels between them amber at best. |
| Anchor `.points` rms > 3 m (169: 29 m) | Pose from geometry instead; `anchors.csv` flags the stretch; the anchor still never moves. |
| Puvi polygon missing (152 sheets) | Neighbour-window path; red only when no neighbour is placed either. |
| Near-symmetric outline, sides unknown | Both twins written (solid + dotted); red; team decides in the Georeferencer. |
| Image ambiguous or observable in one direction only | Image verifies only; note with bearing and along-track σ; never moves the pose. |
| Team `.points` with one GCP | Refit refuses with "one GCP for <survey>"; parcel unchanged. |
| GeoPackage held by QGIS (any group, any project) | `files.py` finds the layers over the socket, aborts if any is in edit mode, else removes them, deletes stale `-wal/-shm`, writes temp + rename, re-adds with saved style; without the socket a `PermissionError` is listed at the end of the run. |
| Review layer or `.points` being edited | `--refit` and rewrites refuse until saved. |
| QGIS project open is not the `--project` one | Review group not built; message names the open project; files, status and worklist still written. |
| Tile download fails | Export retried with 4 connections, then aborts with the failed tile list; nothing else runs. |
| Sheet overlaps railway land | Listed in `rail_conflicts.csv`; not clipped. |
| Tracker open | Rows queued in `tracker_pending.csv`. |
| Any exception on one parcel | Recorded as failed in `georef_status.csv`; the run continues. |

## 10. Testing

* **Truth pose** = the scale-1 fit of the sheet to the manual file (or to the team's `.points` when its rms
  ≤ 3 m), with its rms and a jackknife heading spread recorded per parcel. Tolerances per parcel:
  max(1.5 m, rigid rms) and max(1.5°, heading spread); parcels with rigid rms > 3 m (47B, 48B, 42A, 43B,
  169; Thailavaram 55, 56) are reported but not counted.
* **Leave-one-out on a copy** of the village folder (`--root`), two modes per hidden parcel: full (its
  manual neighbours frozen) and neighbours-as-window-only (anchor weight 0), for the 13 Kizhikaranai and
  13 Thailavaram manual parcels. Every parcel reported green must pass its tolerance in the mode that
  produced it; no flipped or slid parcel may be green; the per-parcel table (colour, share, margin,
  observability, errors, anchor residuals) is a required artefact and the threshold calibration is derived
  from it before the engine is tuned.
* Unit tests (pytest): raster export writes a georeferenced tile of known size; `.points` round trip
  against a team file; start pose recovers a synthetic rotation and keeps the twin for a rectangle; share
  score is maximal at the true pose on a synthetic edge map and the margin test fires on parallel lines;
  observability flags a one-direction match; point-to-line observations against a stretched anchor
  converge to the sheet pose; refit: 0 auto GCPs + 2 team GCPs → pose within 0.3 m of the team points,
  5 wrong auto GCPs + 2 team GCPs → same, a neighbouring auto parcel wrong by 3 m leaves a team-gcp parcel
  unmoved, 1 team GCP → refused; fingerprints after topology; `files.py` on a GeoPackage held by a live
  QGIS layer.
* Integration: full Kizhikaranai run leaves every manual file and `.points` byte-identical, classifies all
  buffer parcels, orients 169/170 from the transcription, and every auto parcel's boundary residual to each
  frozen anchor is ≤ 0.30 m + that anchor's rms or the parcel is not green.

## 11. Measurement attributes

As §11 of the 2026-09-17 design: `area_sqm/are/hect/acre/cent`, `perimeter_m`, `fmb_area_*` preserved,
per-edge `length_m`, `bearing_grid_deg`, `bearing_true_deg`, invariants asserted after every placement.
`georef_status.csv` adds position σ and heading σ from the fit covariance.

## 12. Risks

| Risk | Handling |
|---|---|
| The team's manual files are affine (scales 0.83–1.15; 169 1.21 × 0.59) while the user requires FMB-exact dimensions | Anchors frozen as instructed; their stretch is reported in `anchors.csv`; the tool can write, on request, a scale-1 re-fit of the team's own GCPs beside the manual file (`<survey>_parcels_rigid.gpkg`) for the team to compare — never in place. **Needs the user's decision.** |
| Google imagery is itself 1–3 m off | Accepted as the visual base by the user; anchors override imagery wherever they touch. |
| Image refinement slides along parallel features (46B slid 18.7 m) | ≤ 1.5 m / 2° from the reference pose; margin and observability tests; never decides alone on a ridge. |
| Near-symmetric strips flip 180° | Twins; orientation from printed neighbour sides and frozen neighbours; red if undecided. |
| Vegetation hides boundaries (46B: 12 %) | Red → team GCPs in the Georeferencer. |
| Puvi far off for railway land (170: 37 m) | 40 m window; anchor-implied hypotheses; Puvi never moves a pose. |
| Windows file locks and WAL sidecars | `files.py` (§9). |
| Tile terms of use | Export for the team's internal review only, dated, not redistributed. |
| GDAL tile driver missing on a machine | Fallback: PyQGIS `QgsRasterFileWriter` on the XYZ layer via the socket. |

## 13. Estimate and order of work

Five to seven working days. Day 1: move the 2026-09-17 modules into `autogeoref\`, `anchors.py` with the
`.points` reader, threshold calibration on the 26 manual parcels. Days 2–3: raster, window, edges, align
with margin and observability. Day 4: fit and block adjustment with point-to-line anchors, `.points`
writer and refit. Day 5: `files.py`, review group, status, tracker, worklist, topology integration.
Days 6–7: leave-one-out harness and the two test modes, fixes.

## 14. Puvi disagrees with the portal for 22 surveys (added 2026-09-18, after the 569B report)

The user reported that Thirukatchur 569B is a large shape in Puvi and a small, different shape on
the FMB sheet. Both datasets are internally consistent, and they describe different parcels:

* **Portal:** giscode `S3504074569B` draws a 0.98 acre wedge, 34 x 218 m, printing the neighbours
  613, 614, 584, 101A, 565, 564, 563 - one tight cluster beside the GST Road. Re-fetching the same
  giscode at 1:2000 returns the same parcel at the same ground size, so the conversion and the
  scale handling are correct and the drawing is not cut off. The portal also serves a **569C**
  (0.96 acre, the same 34 x 218 m outline) that is absent from Puvi, and refuses plain `569`.
* **Puvi:** 569B is the 259 acre village tank (1064 x 1648 m, 96 vertices, two holes). Its
  neighbours are a completely different set and it is not adjacent to 569A. Puvi 569A + 569B
  equals the project's recorded area for base survey 569, so Puvi is self-consistent too.
* A 259 acre parcel cannot be drawn on an A0 page at 1:500, yet the portal serves 569B at 1:500
  without complaint - another sign that the portal's 569B is the small parcel.

Corridor-wide audit of every sheet against its Puvi polygon
(`FMB_Vector\puvi_vs_sheet_qc.csv`):

| Verdict | Sheets |
|---|---|
| agree within a factor of 1.5 | 559 |
| partial disagreement (1.5x to 3x) | 49 |
| identity mismatch (>= 3x or <= 1/3) | 22 |
| sheet vectorised to a fragment (< 150 m2) | 22 |

Eighteen of the 22 mismatches are letter subdivisions, and the portal often subdivides further
than Puvi: village 35_04_052 has separate sheets 53B and 53C against a single Puvi polygon, exactly
as survey 569 has A, B and C on the portal against A and B in Puvi.

**Consequences for this design.** The area gate in section 2 is required, not optional: without it
the search window and start pose for 569B would sit 700 m from where the sheet belongs.
`georef_status.csv` gains `puvi_trusted` and `puvi_ratio`. The 22 fragment sheets are a separate
defect in the PDF-to-vector conversion, tracked in the same QC file; they must not be
georeferenced until they are re-converted.

