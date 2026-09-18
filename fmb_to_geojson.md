# fmb_to_geojson — Tamil Nadu cadastral PDF → GeoJSON

Converts a Tamil Nadu **FMB / survey sketch** or a **CMDA / DTCP layout-approval
plan** (PDF) into polygonised GeoJSON with plot-number attributes.

These PDFs are *vector* drawings, not scans. The parcel boundaries are real line
segments, so the tool recovers **true geometry** — angles, lengths and shapes as
drawn — instead of tracing pixels. That is the whole reason the output is
trustworthy at the ~0.1 % level (see *Accuracy*).

---

## Install & run

Needs Python 3 with `pdfplumber` (and `matplotlib` for the preview PNG). No
network required.

```bash
# 1. read the title block to get scale / area / village
python fmb_to_geojson.py sketch.pdf --show-titleblock

# 2a. FMB survey sketch — plot numbers from glyph outlines
python fmb_to_geojson.py sketch.pdf --scale 500 --area 11853 \
    --village "Koladi [70]" --taluk Poonamallee --survey 8 \
    --glyphlib glyphlib.json -o out/

# 2b. CMDA / DTCP layout plan — geometry only, drawing is in a corner of the sheet
python fmb_to_geojson.py layout.pdf --scale 800 --kind layout \
    --clip 280 780 850 1220 -o out/
```

### Key options

| flag | meaning |
|---|---|
| `--scale N` | sheet scale denominator (`500` = 1:500). **Required.** Read it off the title block. |
| `--area M2` | stated area in m², enables the closure check in the run summary |
| `--kind` | `auto` (default), `sketch`, or `layout`. `auto` picks `layout` when the PDF has a real font layer (>2000 chars) |
| `--glyphlib F` | glyph-exemplar JSON, required to read plot numbers on **sketches** |
| `--clip X0 X1 Y0 Y1` | keep only geometry inside this box (PDF points, y-up). Use on layout sheets where the drawing occupies one corner |
| `--snap P` | node-snap tolerance in PDF points at A0 size; auto-scaled for A4. Default 2.7 |
| `--cls upper\|lower` | which of the two text sizes on a sheet holds the plot numbers (see *Plot numbers*) |
| `--survey`, `--village`, `--taluk` | written into the output metadata |

### Outputs (into `--out`)

* `<stem>_parcels.geojson` — the polygons, with `plot_no`, `area_sqm/are/hect`,
  `perimeter_m`, `n_vertices`, `n_holes`, `review_needed`, `label_source`
* `<stem>_parcels.csv` — the same attribute table, sorted by `plot_no`
* `<stem>_lines.geojson` — every classified line segment (`layer`, `length_m`, …)
* `<stem>_preview.png` — labelled parcel map for a visual QC glance

---

## How it works

### 1. Line extraction & classification (`extract_lines`)
Drops page furniture (border, title block) by **frame geometry**, with
thresholds that scale to the sheet so A0 and A4 prints behave identically. Each
surviving stroke is classified by the FMB drawing convention:

| layer | how it's identified | used for polygons? |
|---|---|---|
| `survey_boundary_main` | heaviest solid-black stroke on the sheet | ✅ |
| `subdivision_line` | lighter solid black | ✅ |
| `tie_line` | dash-dot-dot pattern (surveyor check measurements) | ❌ excluded |
| `blue_feature_line` | blue (channels, streams, flow arrows) | ❌ excluded |

Excluding tie lines and blue features is essential: treating tie diagonals as
boundaries shattered whole fields into phantom parcels in early runs.

### 2. Polygonisation (`polygonize`)
Split segments at crossings and T-junctions → snap coincident nodes within
tolerance → prune dangling (dead-end) edges → walk the planar half-edge graph
into faces → drop the unbounded outer face → re-attach any **nested island**
(a separate component enclosed by a face) as a **hole**, so it isn't
double-counted. The last step matters: a missed island inflated one sheet's
area by 6.5 % before it was handled.

### 3. Plot numbers
* **Layout plans** have a real text layer → numbers come from
  `page.extract_words()`.
* **FMB sketches** have *no* font layer → digits are vector glyph outlines. The
  tool normalises each glyph to a small bitmap and matches it against a library
  of verified exemplars (`--glyphlib`). A run is treated as a plot label only if
  it is horizontal (dimensions run tilted along their edge), bracket-free
  (bracketed values are tie distances), and every glyph matched.

Two text sizes coexist on a sheet — plot numbers and edge dimensions — and
**which size is the plot number is not constant**. It's the larger text on the
Koladi sheets but the smaller on some Ayanambakkam sheets. `--cls upper|lower`
selects the class; a 1-D k-means on glyph diameter finds the split.

### 4. Spatial join
Each label is assigned to the smallest polygon whose interior contains its
centroid. One label per polygon:
* sketch: the glyph recogniser already returns whole labels
* layout: a plot holds many dimension digits, so numeric tokens are **not**
  concatenated or guessed — only a single unambiguous named reserve
  (PARK / PP1 / PP2 / OSR / EWS) is taken. See *Limitations*.

### 5. Coordinates
Two conformal (shape-preserving) transforms only: a **Y-flip** (PDF y is
top-down) and a single **uniform scale** to ground metres
(`metres_per_pdf_point = 25.4/72 × scale / 1000`). No rotation, grid-snapping or
shape simplification. The one bounded exception is node-snapping in step 2,
whose maximum actual displacement is reported every run.

---

## Coordinate reference & georeferencing

Output is a **local planar CRS in ground metres**, origin at the page
bottom-left. It is **not georeferenced** — GeoJSON nominally expects WGS84, so
QGIS will load these as raw planar coordinates. To place them on the ground,
apply a 2-point **Helmert (similarity) fit** from known ground control points;
that keeps the geometry rigid.

---

## Accuracy

On sheets where every edge and tie is dimensioned, the recovered lengths match
the printed values to **~0.1 m** and reconstructed perimeters to **~0.06 %**.
So when a parcel's *area* differs from the register by a few percent, that's a
**sheet-vs-register** difference (FMB sketches are not precise cadastral
surveys), not an extraction error. Always sanity-check with `--area`.

Two register conventions seen, worth settling per village before reconciling
areas:
* **Koladi**: roads carry their own subdivision numbers and are counted in the
  registered extent.
* **Ayapakkam / layout plans**: the registered "plottable" area **excludes**
  road poromboke, so the raw polygon total runs high by the road fraction
  (e.g. +13 % to +26 %). Subtract the road polygon(s) before comparing.

---

## Limitations (read before trusting attributes)

* **Geometry is reliable across the board** — every scale tried (1:334 to
  1:1177), A0 and A4. Plot-number recovery is what varies.
* **Plot-number coverage tracks drawing density.** Sparse sheets come out fully
  labelled; dense ones come out partially labelled, because plot numbers and
  dimension text overlap and can't be separated by position or size alone.
* **Layout plans: numbers are not auto-extracted** when dense. The tool labels
  named reserves and leaves numbered plots blank rather than emit a wrong guess.
  Read those numbers off the sheet, or add a verification pass.
* **`--glyphlib` is required for sketch plot numbers** and must contain
  exemplars for every character that appears (digits plus the subdivision
  letters A–F, G–K, etc.). Build it once from a sheet you've verified by eye,
  then reuse it across that series.
* **DWG files can't be read here.** Export to **DXF (ASCII, R2013 or older)**
  from the CAD app, or convert with the free ODA File Converter, then feed the
  DXF to a DXF reader — a DXF carries real layers, text and possibly true
  coordinates, which removes most of the guesswork above.
* **`review_needed = true`** flags a polygon where more than one candidate label
  landed inside — usually a missing dividing line. Check these against the sheet.

---

## Building a glyph library (sketches only)

The library is a JSON list of `[char, bitmap, aspect_ratio]` exemplars. Bootstrap
it from one sheet whose numbers you've confirmed visually:

1. run `build_glyphs()` on the verified sheet to get glyph runs,
2. pair each run with the known label string,
3. store each glyph's `bmp` and `ar` under its character.

Add new characters (e.g. new subdivision letters G/H/J/K on a later sheet) by
verifying those runs and appending their exemplars. Re-check an already-verified
sheet after any change to confirm no regression.

---

## Changes made on 2026-09-16 for the Collabland portal sketches (Tambaram - Chengalpattu corridor)

The portal PDFs (fetched with fetch_fmb_gui.py, A0, 1:500) differ from the Koladi/Ayapakkam sheets the tool was written for:

* **No text layer at all** - even the title block is glyph outlines, so `--show-titleblock` prints nothing. Scale is what was requested from the portal (500); `batch_fmb_to_vector.py` cross-checks each sheet's polygon total against the village vector layer and re-runs at another standard scale if the area ratio demands it.
* **Three text sizes, only ~6 % apart** (tiny marks ~2.4-3.6 pt, edge dimensions ~5.9-6.9 pt, plot numbers ~6.7-8.1 pt). The old 2-means split with a 15 % ratio test merged dimensions into the plot-number class. New `size_groups()` finds size bands by gaps in the sorted glyph diameters; `--cls auto` (now the default) tries every band plus the merged text band and keeps the one that yields one label per polygon (score = polygons with exactly one label minus polygons with several).
* **Every dimension carries a decimal point** (a ~1 pt round blob just after the digit group, near the baseline); plot numbers never do. Runs with such a dot are dropped before the spatial join.
* **Inner counters of 0/6/8/9** are sometimes emitted as separate blue shapes; the library has them as the pseudo-character `~`, which `recognise()` skips.
* When several labels still land in one polygon, `plot_no` is the one nearest the polygon centroid (dimensions hug the edges) and the rest go to the new `alt_labels` field with `review_needed = true`.
* A sheet with no closed polygons no longer crashes the CSV writer; `run()` exposes `run.last_stats` and records `text_class` in the metadata.

`glyphlib.json` was built from 30 of these sheets: 6,900 blue glyph shapes clustered by bitmap similarity, clusters labelled by eye from an oriented render, thinned to 213 exemplars (digits, brackets, A/B/C, `~`). Letters beyond C were not present in the sample; add exemplars if a village uses them.

### Batch tooling (same folder)

* `batch_fmb_to_vector.py` - runs every sketch under `FMB_Sketches\<dd_tt_vvv>\`, one output set per sheet under `FMB_Vector\<dd_tt_vvv>\`, plus `fmb_vector_summary.csv`. File names are plain (`<survey>.pdf`, one per survey unit); a sheet fetched at another scale is listed in `FMB_Sketches\sheet_scales.csv` and converted at that scale. `--sheets a/b,c/d --merge` re-runs selected sheets and patches their rows into the summary.
* `scale_from_dimensions.py` - reads the edge dimensions (digits + decimal point, stitched across the point) and compares annotated metres with drawn length. Used as a **verification only**: `verified` (ratio within 6 % of 1), `unverified` (fewer than 3 dimensions read - e.g. unsubdivided surveys or dimensions written along a tilted edge, which the upright glyph library cannot read), or `MISMATCH - check manually`. The portal honours the requested scale, so the sheet is converted at the requested scale regardless.
* `offpage_check.py` - lists sheets whose survey boundary reaches the page margin (clipped strip surveys). `fetch_offpage_at_2000.py` re-fetches those at 1:2000 (saved as `<survey>_s2000.pdf`; after checking they fit, rename over the clipped file and record the scale in `sheet_scales.csv`).
* `enrich_attributes.py` - adds district/taluk/village codes, giscode and sketch id to every polygon, rewrites the per-sheet CSVs, and builds `all_parcels_attributes.csv/.xlsx` and `FMB_Vector_all_sketches.gpkg` (one layer per sketch; coordinates are still local sheet metres). Single-polygon sheets (unsubdivided surveys) get `plot_no` = the survey number with `label_source = whole_survey_no_subdivision`.

### Project folder layout (Tambaram - Chengalpattu corridor)

    D:\Projects\Tambaram_Chengalpattu      Tambaram_Chengalpattu_Railway.gpkg      rail_line, rail_buffer_30m
      Railway_Buffer_Vector_Plots.gpkg        vector_in_buffer_30m (677 plots), buffer_no_vector_data
      Tambaram_Chengalpattu.qgz               standalone project
      Georeferencing_Tracker.xlsx             one row per portal survey unit, links to PDF and GeoJSON
      FMB_Sketches\<dd_tt_vvv>\<survey>.pdf   627 sketches; fmb_sketch_summary.csv, fmb_missing_sketches.csv, sheet_scales.csv
      FMB_Vector\<dd_tt_vvv>\<survey>_parcels.geojson / _parcels.csv / _lines.geojson / _preview.png
                                              all_parcels_attributes.csv/.xlsx, FMB_Vector_all_sketches.gpkg, fmb_vector_summary.csv
      _logs\                                  download logs, batch inputs and run logs (safe to delete)
