# FMB sketch auto-georeferencing by mosaic chaining — design

Date: 2026-09-17 (revised the same evening after the critics' review and the Thailavaram test). Status: awaiting written review.
Project: Tambaram – Chengalpattu rail corridor, `D:\Projects\Tambaram_Chengalpattu`.
Code home: `D:\code\FMB_to_GeoJSON`.

## 1. Goal

Place the 655 FMB sketch polygon files (currently in local sheet metres, one file per
portal survey unit) onto the ground in EPSG:32644 automatically, so that the team only
seeds one sketch per block by hand and reviews the few placements the tool is not sure
about. Every sketch keeps its FMB dimensions and angles exactly: a placement is a
rotation plus a shift, never a stretch.

## 2. Decisions taken with the user

| Question | Decision |
|---|---|
| Rigidity | Sketches are only moved and rotated (scale fixed at 1). Residual gaps/overlaps between neighbours are resolved afterwards by the existing topology fix. |
| Position references | Neighbour sketches (adjacency) drive placement. The Puvi village vector is used only as adjacency fallback and for the low-confidence approximate fallback position. Google Satellite is the visual review basemap only. |
| Anchor | The team places one seed sketch per connected block by hand (as today). Existing hand-placed sheets are seeds, used only as a **pose** (rotation + shift recovered by a rigid fit of the sheet to the hand file), never as geometry: hand files carry scale 0.95–1.09 and edited vertices. |
| Approach | A + C: mosaic chaining outward from seeds reaches new sheets; a whole-block rigid least-squares adjustment (C) then settles every block before the topology fix. C was promoted from “later” to “required” after chain drift and loop-closure conflicts appeared on Thailavaram (§3). |
| Output naming | Same as the team's hand-placed files: `FMB_Vector\<village_code>\<survey>_parcels_modified.gpkg`. |
| Neighbour graph | Black neighbour numbers printed on each sheet, with Puvi adjacency as fallback. Sample verified on 5 sheets. |
| Run form | A command per village; QGIS review group. |

## 3. Evidence the design rests on

* Hand placements recovered from the seven Kizhikaranai sheets: rotations −22° to −26°,
  but the hand work also applied scales of 0.95–1.09 and 1.4–3.7 m residual vertex edits.
  The tool removes that drift by keeping scale at 1.
* **FMB sheets of one village do not share a drawing frame.** The team's 13-sheet
  Thailavaram file (35_04_052) gave recovered headings from −76° to +6°. Kizhikaranai only
  looked uniform by chance. No village-wide rotation prior or heading gate is used anywhere.
* Block adjustment on those 13 sheets (scratch `block_adjust.py`, 2026-09-17): 12 shared
  boundaries fitted to 0.02–0.61 m rms, shifts from the hand poses 0.5–2.5 m, headings
  changed ≤ 2.3°; pairwise overlap between surveys fell from 1323 m² (rigid sheets at hand
  poses) to 769 m² before the topology fix and 0 m² after it, with 0.2 m² of enclosed gaps.
  One 3-vertex chain (52–53A) could not be fitted better than 2.8 m rms: the two sheets
  disagree on a corner angle, so both are flagged. Surveys 55 and 56 have no equal-length
  shared edge with any sheet and stay at the hand pose; as rigid sheets they overlap by
  660 m², so at least one hand position is wrong by about 5 m.
* Adjacent sketches share identical boundary edge lengths (47A–171: 94.72/26.63 vs
  94.56/26.64 m; 47B–171: 29.03/80.46 vs 29.28/80.43 m; 47B–48A: 48.21/5.58 vs
  48.51/5.47 m). Two matched vertices fix rotation and shift completely.
* Neighbour survey numbers are printed in black just outside each survey outline
  (sheet 145: 223 N, 461 NE, 463 SE, 138 S; 47B: 48A N, 171 NE, 76 W, 46A S). First-pass
  recall ≈ two thirds; dense sheets and slashes ("103/5A") need library additions.

## 4. Scope

In scope: one village at a time; sketches under `FMB_Vector\<village_code>\`; seeds;
placement; confidence; review layer; status CSV; tracker update; topology fix trigger.
Out of scope: imagery matching, Puvi-driven positioning, editing hand-placed files,
corporation areas without sketches.

## 5. Inputs and outputs

Inputs (all existing):
* `FMB_Vector\<village_code>\<survey>_parcels.geojson` — sketch polygons in sheet metres
  (from `fmb_to_geojson.py`), plus `_lines.geojson`.
* `FMB_Sketches\<village_code>\<survey>.pdf` — for reading neighbour numbers.
* `Railway_Buffer_Vector_Plots.gpkg` / Puvi village vectors — adjacency fallback and
  approximate fallback position only.
* Seeds: any `FMB_Vector\<village_code>\<survey>_parcels_modified*.gpkg` not written by the
  tool, or written by the tool and since edited (fingerprint mismatch).

Outputs:
* `FMB_Vector\<village_code>\<survey>_parcels_modified.gpkg` (EPSG:32644), attributes as the
  sheet file plus `georef_method` (seed | chained | approximate | manual),
  `georef_confidence` (0–100), `georef_residual_m`, `anchored_to` (comma list of
  neighbours used), `georef_run` (ISO timestamp of the run).
* `FMB_Vector\<village_code>\georef_status.csv` — one row per sketch: status
  (seed / placed / approximate / waiting / failed), confidence, residual, anchors,
  neighbours read, notes, fingerprint of the file the tool wrote.
* `FMB_Vector\<village_code>\neighbours.csv` — sketch, neighbour, side, source
  (sheet | puvi | both).
* QGIS group "Georef review – <village>" (see §9).

Rule: the tool never overwrites a `_modified` file whose geometry fingerprint differs from
the one it recorded, or that it has no record of. Such files are anchors. Fingerprints live
only in `georef_status.csv`, never inside the GeoPackage. Any `_modified*` variant the team
saves (`_modified_1`, `_modified (2)`, …) counts as a hand file for that survey; the newest
by modification time is the pose. Files that QGIS holds open are never written in place: the
tool writes a temporary file beside it and renames, and skips (with a message) when the
rename fails.

## 6. Components

Each unit is a module in `D:\code\FMB_to_GeoJSON\autogeoref\` with one job and a testable
interface.

### 6.1 `outline.py` — sketch outline
`outline(features) -> Outline` : dissolves a sketch's polygons, takes the exterior ring and
returns every vertex with edge lengths and bearings in sheet coordinates. **No collinear
merge**: merging moved vertices differently on the two sheets of a pair and broke matches;
near-straight corners (turn < 5°) are instead passed through by the matcher. Works
identically on placed (ground) geometry.

### 6.2 `neighbours.py` — neighbour graph
`read_sheet_neighbours(pdf) -> [(number, side_bearing, distance_m)]` : black glyph runs
(same clustering as `build_glyphs`, black instead of blue), outside the outline, within
60 m; recognised with `glyphlib.json` extended with "/" and letters. Runs containing "?"
are dropped, not guessed. Labels with a slash ("103/5A") are reduced to the survey part
before the slash. A number that is not a survey of this village (neighbouring village,
road, tank) is kept in `neighbours.csv` with source `sheet-external` and never used for
placement.
`puvi_neighbours(village_code, survey) -> [number]` : Puvi polygons touching the survey's
Puvi polygon (survey number matched on the numeric part, so 47A → 47).
`build_graph(village_code) -> neighbours.csv` : union of both sources with the source
recorded; sheet labels win on side.

### 6.3 `match.py` — shared-boundary matching and rigid fit
`common_chain(outline_q, outline_p, side_hint) -> [(length, [(vertex_q, vertex_p)])]` :
walks both outlines by chainage and returns **every** run (≥ 8 m) of consecutive equal-length
edges (tolerance 0.30 m + 1 %), in both traversal directions, passing through near-straight
corners (< 5°). Returning only the longest run hid true matches. Candidates are then
filtered: the placed sheet must not overlap the anchor by more than 0.6 m mean width over the
matched length (relative, not an absolute area cap), the anchor's printed neighbour label
must lie on the placed side within 30 m, and the Puvi centroid of the same survey number
must be within 100 m (numeric part of the number only; Puvi is a gate, never a position).
`rigid_fit(pairs) -> (theta, tx, ty, rms, max_residual)` : least-squares rotation + shift
with scale fixed at 1 (Procrustes without scale). Pairs from several neighbours are pooled.
A chain of two vertices (one edge) can slide or mirror: such placements are capped at
confidence 79, and a sheet placed on one edge only is "tentative" until a second neighbour
or the block adjustment confirms it.

### 6.4 `engine.py` — placement
`place_village(village_code, seeds) -> status` :
1. Load outlines; recover each seed's pose from its hand file (rigid fit of the sheet to
   the hand geometry, poly_id pairing, nearest-vertex refinement within 4 m) and load it as
   placed with `georef_method` seed/manual.
2. Chain loop: for every unplaced sketch with ≥ 1 placed neighbour, match against all placed
   neighbours (§6.3), pool pairs, fit, score. Place the highest-scoring candidate above the
   threshold (50). Repeat until none qualifies. No rotation prior of any kind.
3. Block adjustment (approach C, `adjust.py`): unknowns θ, tx, ty per sheet; observations =
   all shared-boundary vertex pairs (σ 0.30 m) plus each sheet's current pose as a weak
   prior (σ 2.5 m position, 3° heading; seeds get the hand pose). Solve in a frame centred
   on the block (raw UTM unknowns stall the solver), robust Huber loss. Chains are admitted
   trusted-first: multi-vertex chains (≥ 3 pairs) go in first, any that still exceed 1.0 m
   rms is demoted (the two sheets disagree), then single-edge chains are tried one at a
   time and kept only if they fit within 1.0 m without worsening the others. Rejected
   chains are reported as sheet disagreements.
4. Sketches with a placed neighbour but no chain: approximate placement = Puvi centroid
   of the same survey number, heading from the best-matching neighbour's edge bearing,
   `georef_method` approximate, confidence ≤ 30. Sketches with no placed neighbour: waiting.
5. Write outputs (§5) and the status CSV.

Confidence (0–100): 40 × f(boundary rms ≤ 0.30 m) + 30 × f(matched length ≥ 20 m) + 30 ×
(number of neighbours matched, capped at 2) / 2; capped at 79 when only one edge was
matched, capped at 40 if the placed sketch overlaps any placed sketch by more than 0.6 m
mean width over the shared length, and 40 flat for a sheet kept at its hand pose without
boundary evidence. f() is linear between full and zero credit.

### 6.5 `review.py` — QGIS review group and tracker
`build_review_group(village_code)` via the QGIS socket: loads all `_modified` files of the
village into "Georef review – <village>", styles green (≥ 80), amber (50–79), red (< 50 or
approximate), labels `survey (confidence)`, lists waiting sketches in the layer group name.
`update_tracker(village_code)` : sets Status = Georeferenced and Remarks = confidence for
placed sketches in `Georeferencing_Tracker.xlsx` (only when the workbook is not locked).

### 6.6 Topology fix (existing)
`topo_georef_fix` logic from 2026-09-17 generalised to take a list of layers: snap shared
edges (0.30 m), clip overlaps (larger yields), fill enclosed gaps, recompute area fields.
Run per village once no reds remain, by explicit command.

### 6.7 `georef_village.py` — command
`python georef_village.py 35_04_077 [--review] [--topology]` runs §6.2–§6.5; `--topology`
runs §6.6 afterwards.

## 7. Data flow

sheet GeoJSONs + PDFs → outlines, neighbour graph → (seeds) → engine: match → fit →
score → place → repeat → `_modified.gpkg` + `georef_status.csv` → QGIS review group →
team corrects reds → rerun (fingerprints turn edited files into anchors) → topology fix →
tracker.

## 8. Error handling

| Situation | Behaviour |
|---|---|
| Neighbour named on sheet has no sketch (outside corridor, Natham A) | Edge skipped; noted in status. |
| Conflicting transforms from two placed neighbours (> 0.30 m apart) | Place by the better-matching one; confidence −20; note conflict; the block adjustment settles it. |
| Two sheets disagree on a shared corner (chain cannot be fitted within 1.0 m rms, e.g. 52–53A) | Chain demoted; both sheets keep their other evidence; note "sheet disagreement" on both; review. |
| Only one edge matches (chain of two vertices) | Placed tentatively, confidence ≤ 79; sliding/mirror ruled out by the label-side and overlap tests; confirmed or moved by the block adjustment. |
| Sheet read at 1:2000 (strip survey) | Same path; length tolerance and snap tolerance scaled ×4 (1.2 m) for that sheet only. |
| Puvi has 47 but sheets have 47A/47B | Puvi used only as the 100 m centroid gate for both; never as a position. |
| No placed neighbour | Status waiting; tool lists the block so the team seeds it. |
| No chain found but a placed neighbour exists | Approximate placement, confidence ≤ 30. |
| Strip sheet (1:2000) or unsubdivided sheet | Same path; only the outline is used. |
| Hand-edited `_modified` file | Detected by fingerprint; becomes an anchor; never overwritten. |
| Locked tracker / CSV | Skipped with a message; rerun later. |
| Any exception on one sketch | Recorded as failed with the message; village run continues. |

## 9. Testing

* Ground truth: the hand placements are affine and vertex-edited, so they are a loose
  reference, not truth. Seven hand-placed Kizhikaranai sheets (backup copies in
  `_logs\backup_georef_20260917`): seed with 47B only; the other six must land within
  3 m of the hand pose, with every shared-boundary residual ≤ 0.30 m rms. Repeat seeding
  with 171. Thirteen hand-placed Thailavaram sheets (`Downloads\thailavaram.geojson`):
  the block adjustment must reproduce the 2026-09-17 result (12 boundaries ≤ 0.61 m rms,
  0 m² overlap after topology) and flag 52–53A, 55 and 56.
* Unit tests (pytest, `tests\`): outline + collinear merge on sheets 145 and 47A;
  `common_chain` on 47A–171 and 47B–48A; `rigid_fit` on a synthetically rotated sketch
  (recover θ, t within 1e-6); neighbour reading on 145 and 47B (expected sets above);
  fingerprint rule (edited file is not overwritten); confidence formula edge cases.
* Integration: full village run on 35_04_077 with all seven seeds must produce zero
  changes (fingerprints match) and a status CSV with seven seeds.

## 10. Future work (not now)

Whole-block least-squares adjustment (approach C) using the same matched pairs; reading
the railway land polygon against the OSM track as an additional constraint; a QGIS
Processing wrapper for the command.

## 11. Measurement attributes (dimensions, bearings, area) — added 2026-09-17 at the user's request

These values feed later work and must be correct and traceable.

**Invariants.** A rigid placement (rotation + shift, scale 1) leaves every edge length and every
area exactly as on the sheet. The tool asserts this after each placement: per polygon
|area_after − area_before| < 1e-6 m² and per edge |length_after − length_before| < 1e-6 m;
a violation fails that sketch loudly rather than writing it.

**Per-polygon attributes** (in `<survey>_parcels_modified.gpkg`): `area_sqm`, `area_are`,
`area_hect` (existing), plus `area_acre` = area_sqm × 0.000247105381 and `area_cent`
= area_acre × 100 (Tamil Nadu practice: 1 acre = 100 cents), `perimeter_m`, `n_vertices`.

**Per-edge attributes** (new companion layer `<survey>_edges` in the same GeoPackage, one
LineString per outline and subdivision edge): `poly_id`, `plot_no`, `edge_no` (clockwise from the
northernmost vertex), `length_m` (sheet value, unchanged by placement), `bearing_grid_deg`
(clockwise from grid north of EPSG:32644, 0–360, measured after placement),
`bearing_true_deg` = grid bearing + grid convergence at the edge midpoint (computed from the
UTM 44N projection; ≈ +0.1° to −0.4° in this corridor), `sheet_dim_m` (the dimension printed on
the FMB for that edge where the glyph reader found one, else null) and `dim_diff_m` =
length_m − sheet_dim_m. A `dim_diff_m` above 0.3 m flags a drawing/extraction discrepancy for
review. Bearings are given to 0.01° with the placement's heading uncertainty recorded on the
polygon as `heading_uncertainty_deg` (from the rigid-fit residual and chain length).

**Units and datum stated in every file's metadata**: lengths in metres, areas in square
metres/ares/hectares/acres/cents, bearings in degrees clockwise from grid north (EPSG:32644)
and true north, CRS EPSG:32644, vertical datum not applicable.

## 12. Risk register (from the 37-agent critics' review, 2026-09-17)

| # | Risk | Where handled |
|---|---|---|
| 1 | Seeds are affine (scale 0.95–1.09) and vertex-edited; using them as geometry poisons every chain | §2 anchor row, §6.4 step 1: seeds give a pose only |
| 2 | Collinear merge moves vertices differently on the two sheets of a pair | §6.1: no merge; matcher passes near-straight corners |
| 3 | Coincidental single-edge matches (equal length by chance) | §6.3: label-side, overlap and Puvi gates; cap 79; tentative until confirmed |
| 4 | Fingerprints stored in the GeoPackage are lost on "save as" | §5: fingerprints in `georef_status.csv` only |
| 5 | Team saves `_modified_1` / `_modified (2)` variants | §5: any `_modified*` counts; newest is the pose |
| 6 | GeoPackage locked by QGIS while the tool writes | §5: write temp + rename; skip with message |
| 7 | Absolute overlap cap (2 m²) rejects long boundaries and accepts short bad ones | §6.3: 0.6 m mean width over matched length |
| 8 | No single village rotation (headings −76° to +6° on Thailavaram) | §3, §6.4: no rotation prior or heading gate anywhere |
| 9 | Chain drift and loop-closure conflicts | §6.4 step 3: block adjustment C, trusted-first admission |
| 10 | Slash labels ("103/5A") read as a different survey | §6.2: reduce to the survey part |
| 11 | Black glyph clustering picks up other-village numbers and road names | §6.2: `sheet-external`, never used for placement |
| 12 | Puvi numeric keys collide (47 vs 47A/47B) | §8: Puvi is a 100 m gate only |
| 13 | 1:2000 sheets need a wider length and snap tolerance | §8: ×4 for that sheet only |

## 13. Kizhikaranai full-village run (added 2026-09-17, evening)

All 15 sheets in the buffer were placed (7 team seeds as poses, 8 chained, whole block adjusted; 16 shared
boundaries fitted to 0.01-0.39 m rms). What the run changed in the design:

* **Neighbour numbers are transcribed, not glyph-matched.** Two independent readers per sheet (quadrant crops,
  page-up = N) agreed on every number; the glyph reader had returned fragments ("17" for 170, "4A" for 42A) and
  the village number of the adjoining village ("V.No. 74 Thirukachur") as a survey. `nb_override_<village>.json`
  is merged into the neighbour graph; with a complete transcription, two sheets that do not name each other are
  not neighbours and no chain between them is used (placement and block adjustment).
* **Placement order**: candidates with >= 3 vertex pairs or >= 2 anchors are always placed before any
  single-edge candidate. Strips beside the railway share 205 m straight edges that cross-match; ordering by
  confidence alone let a single-edge placement of the railway strip block the corroborated placements.
* **Railway land** (surveys whose Puvi plot contains the track: 169, 170, 171) overlaps the older parcel sheets
  by up to 2700 m². The topology fix must not clip these (it removed 30 % of 42A): railway-vs-parcel overlaps are
  left in place and listed for the team; parcel FMB areas equal their Puvi areas within a few percent, so the
  sheets are consistent and the conflict is a land-record question. The traced rail line is a score term only.
* **Puvi is a reference column** (`puvi_reference_m`), never a gate or a position (user rule).
* Overlap allowance 8 % of sheet area; a group with >= 3 pairs and both numbers printed at the boundary passes the
  overlap test with a 10-point penalty.

