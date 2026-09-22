# SOP: georeferencing the FMB parcels onto the fine-tuned Puvi base

Date 2026-09-22. Status: for Akash's approval. No code written for this yet.
Every number below was measured today on the corridor data, not estimated.

## 1. What this produces

For each survey in the 30 m rail buffer, the parcel as its FMB sheet draws it, placed on the
ground. The sheet gives the geometry, which is the legal record and carries the subdivision plots;
the fine-tuned Puvi layer gives the position. Two layers are delivered side by side, as agreed:
the placed FMB parcels, and the fine-tuned Puvi beside them, so nothing is lost where the two
disagree.

## 2. What we have

| | count |
|---|---|
| parcels in the fine-tuned base | 678 rows, 664 distinct surveys, 20 villages |
| with a converted FMB sheet | 625 |
| with no sheet | 39, including all 9 of Natham A |
| already placed by your hand | 53, in Thirukatchur (38) and Kizhikaranai (15) |

The 39 without a sheet are listed and left for later, as agreed. They stay in the Puvi layer.

## 3. The method

### Stage 1: place each sheet on its base parcel

The sheet is fitted onto its fine-tuned Puvi parcel by rotation and shift only, with scale fixed at
1 so every printed length survives, choosing the turn that makes the two overlap most.

Measured on the three villages where your placements give the answer:

| Village | sheets fitted | shape match | area ratio, sheet to base | placed error against your placement | the base's own error |
|---|---|---|---|---|---|
| Thirukatchur | 66 | 0.86 | 1.00 | 2.8 m | 2.8 m |
| Kizhikaranai | 15 | 0.90 | 1.01 | 4.2 m | 4.2 m |
| Thailavaram | 31 | 0.91 | 0.96 | 2.7 m | 2.7 m |

Two things follow. The sheet and the Puvi parcel agree on shape closely enough to fit without
ambiguity. And the placed parcel inherits the base's position exactly: this stage cannot be more
accurate than the base, and it is not less.

### Stage 1b: where the shapes disagree, the base shape stands

Akash's rule, 2026-09-22: try the shape match first; where a parcel's shape would have to change to
fit, do not force it. The base parcel's shape is then taken as fixed and delivered, and the sheet is
recorded as disagreeing rather than placed on top of it.

The test is the same one that flags a parcel: a shape match under 0.6, or an area ratio outside 0.8
to 1.25. Those parcels carry `geometry_source = puvi base` in the output; everything else carries
`geometry_source = fmb sheet`. Nothing is silently reshaped either way.

### Stage 1c: ground control points from the base layer

Control points come from the base layer itself, as instructed: its corners, where its outline really
turns by 20 degrees or more, paired with the matching corner of the placed sheet within 5 m.

Measured on the placed sheets:

| Village | parcels | corners on the sheet | corners on the base | paired | residual | parcels with 3 or more pairs |
|---|---|---|---|---|---|---|
| Thirukatchur | 66 | 4 | 4 | 2 | 2.7 m | 31 of 66 |
| Oorapakkam | 64 | 4 | 4 | 3 | 2.6 m | 33 of 64 |
| Kizhikaranai | 15 | 3 | 4 | 2 | 2.9 m | 1 of 15 |

Read that honestly: the residual is the base's own accuracy, so these points do not make the
placement better than the base. What they give is a placement you can open in the QGIS Georeferencer
and drag, in the format the tool already writes (`<survey>_parcels.geojson.points`). About half the
parcels get three or more pairs; the rest get the pairs that exist and the count is recorded, so a
parcel resting on one point is never mistaken for a fitted one.

Where a parcel has no usable pair, the sheet keeps the pose from stage 1 and is marked
`gcp_points = 0`.

### Stage 2: make neighbours agree

Sheets placed one at a time do not meet. Measured after stage 1:

| Village | neighbouring pairs | overlap, median and 90th | gap, median and 90th |
|---|---|---|---|
| Oorapakkam | 145 | 0 m², 90 m² | 0.00 m, 6.6 m |
| Thirukatchur | 167 | 0 m², 86 m² | 0.09 m, 66 m |

Most pairs already meet. The tail does not, because the base's error varies from parcel to parcel
and because some sheets describe a different extent from their Puvi parcel.

So each parcel is adjusted against the parcels it touches: which parcels are neighbours comes from
the Puvi fabric, the shared boundary geometry comes from the sheets, and the adjustment is the
block adjustment already in the tool, rotation and shift per parcel, no scaling, your 53 parcels
held fixed as control. This is the existing `fit.block_adjust`, used with a different source of
neighbour pairs.

**Gate:** adopt stage 2 only if, measured on the 53 parcels you placed, it reduces the error rather
than merely closing the gaps. If it closes gaps but moves parcels away from your placements, it is
recorded as a failure and stage 1 is delivered alone.

### Stage 3: topology

The clip and fill already used on this corridor: overlaps clipped out of the yielding side, gaps
under 25 m² given to the parcel bordering them most, every edit refused if it would empty a parcel
or take more than a quarter of it. Your parcels are never edited.

### Stage 4: evidence, QC and scoring

Every parcel carries what it rests on, as the fine-tuned layer does now:

- `evidence`: `hand placed` for your 53, `measured` where the base is within 150 m of one of your
  placements, `not evidenced` otherwise.
- `shape_match` and `area_ratio` from stage 1, `geometry_source` (`fmb sheet` or `puvi base`), and
  `gcp_points`, how many control points the placement rests on. A parcel whose sheet and base
  describe different ground keeps the base shape and comes to you, rather than being resolved by
  the machine.
- `moved_m`, how far stage 2 moved it after stage 1.

Scoring is the existing `shift_score` against your 53 placements, reported before and after each
stage, with the standing caveat that your own placements disagree with each other by 1 to 3 m, so
nothing can be read finer than that.

## 4. What comes out

In `D:\Projects\Tambaram_Chengalpattu\FMB_on_Puvi_<date>\`:

| File | Contents |
|---|---|
| `FMB_parcels_georeferenced.geojson` | the deliverable: every placed FMB parcel, one row per plot |
| `FMB_parcels_by_survey.geojson` | the same dissolved to one row per survey, for comparison with Puvi |
| `Puvi_Vector_Finetunned.geojson` | the base, copied in unchanged, so the pair travels together |
| `placement_report.csv` | one row per survey: shape match, area ratio, geometry source, GCP count, moves, flags, evidence |
| `gcp\<survey>_parcels.geojson.points` | the control points, in QGIS Georeferencer format, for you to drag |
| `unplaced.csv` | the 39 with no sheet, and any parcel the run refused, with the reason |

## 5. What could go wrong, and what we do about it

| Risk | What we do |
|---|---|
| The base is wrong for a village with no control of yours, so every sheet on it is wrong together | Nothing in this pipeline can detect that. The `evidence` field says so, parcel by parcel. |
| The sheet and the Puvi parcel are different ground (Thirukatchur has known cases) | Flagged by shape match and area ratio, delivered as both layers, left for you |
| A sheet fits its base parcel two ways, turned 180 degrees | The fit reports both scores; where they are within 0.05 the parcel is flagged rather than guessed |
| Stage 2 closes gaps by moving parcels off the truth | The gate in stage 2. It is measured against your 53 before it is adopted |
| A rail parcel slides along the track | Known and unfixable from geometry; those parcels keep the base's along-track position |

## 6. Effort

| Step | Estimate |
|---|---|
| Stage 1 and the report, with tests | half a day |
| Stage 2, the neighbour adjustment and its gate | one day |
| Stage 3 and 4, reusing what exists | half a day |
| Full corridor run and your review | half a day of machine time, then your review |

## 7. Running it

```text
python -m autogeoref.fmb_on_base --corridor --base "D:\Projects\Tambaram_Chengalpattu\Puvi vector fine tunned\Puvi_Vector_Finetunned.geojson"
python -m autogeoref.shift_score <output>.geojson --rail
```

## 8. What this SOP does not do

It does not re-georeference from the sheet's own evidence: no satellite, no neighbour transcription,
no chain matching from a hand seed. Those were measured on this corridor and are recorded in the
fine-tuning README as rejected or as needing work the corridor does not have. This SOP moves the
accurate geometry onto the best positions we currently hold, and says plainly where those positions
are evidenced and where they are not.
