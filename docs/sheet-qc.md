# Sheet QC: check the converted sheets before placing them

Two failures on 2026-09-21 in Thirukatchur came from the converted sheets, not from placement:

- **613** was drawn with one east-side segment at 60.9 m where the sheet (and 569C) print
  70.8 m. The converter keeps drawn geometry, so 613 came out 9.9 m short and 614, 584 and
  92A, hung from it, sat 9.9 m too far north. Fixed by correcting the outline to the printed
  length (recorded in `FMB_Vector/<village>/sheet_corrections.csv`, original under `_logs`).
- **47A** plot 3 came out as an 8 m wedge plus a 0.6 m² sliver where the sheet draws one
  triangle (J, 14, 15). The offset arrows ("2.0", "18.4") are thin black lines with a
  two-stroke arrowhead; polygonised, their feet snapped onto nearby lines and cut the plot.

## What changed in the converter (`fmb_to_geojson.py`)

- `drop_arrows(segs, tol)`: a line with two short strokes at one end, leaning back along it
  by 4 to 60 degrees on opposite sides, is an annotation arrow. It is kept in the lines file
  as layer `annotation_arrow` and left out of the polygonisation.
- `merge_slivers(polys, 2.5)`: faces under 2.5 m² are folded into the neighbour sharing the
  longest edge.
- Both counts go to the GeoJSON metadata (`arrows_dropped`, `slivers_merged`).
- `enrich_attributes.py` now keeps properties it does not know (such as `sheet_correction`).

Re-conversion of Thirukatchur (2026-09-21): 15 sheets adopted (outline unchanged or changed
only by a strip along the boundary that an arrow had cut off), 54 unchanged, 613 skipped
because its outline is hand-corrected. Backups under `_logs/sheet_regen_<stamp>/before`.

## The QC command

```
python -m autogeoref.georef_village <village> --sheet-qc
```

Writes `FMB_Vector/<village>/sheet_qc.csv` and lists, per sheet: slivers (< 3 m²), very thin
plots (compactness < 0.08), spikes (< 8 degrees with arms over 3 m), unnumbered plots,
invalid or multi-part outlines, overlapping plots, holes between plots, and thin lines that
cross the survey boundary away from their endpoints. Thin plots and spikes are often real
(a 183 m strip such as 564, the J-14-15 triangle of 47A); unnumbered plots are label
recognition gaps; crossing lines are the arrow signature or a drawing overshoot. The
printed-length check that would have caught 613 needs the dimension labels and is not part
of this command yet.

Run it before any placement run on a village; treat every flagged sheet as a question to
answer against the PDF, not as a defect to fix blindly.
