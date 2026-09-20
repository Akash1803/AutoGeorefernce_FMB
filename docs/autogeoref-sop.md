# Auto-georeferencing SOP for FMB sheets (Tambaram to Chengalpattu corridor)

Version 3, 2026-09-20. Code: `D:\code\FMB_to_GeoJSON\autogeoref`. Project:
`D:\Projects\Tambaram_Chengalpattu`.

## The idea in one paragraph

A parcel you placed by hand is the truth and is never moved. Every other sheet is placed by
laying it against a neighbour that is already on the ground, matching the shared boundary
edge for edge. Only rotation and shift are applied, so every FMB length stays exactly as
printed. Google Satellite confirms or doubts the result; it never decides it. Each parcel
comes back green, amber or red so you know what to check.

## Before you start a village

1. Convert the sheets. Every survey needs `<survey>_parcels.geojson` in
   `FMB_Vector\<village>\` (the existing fmb_to_geojson batch).
2. Place one or two seed parcels by hand in the QGIS Georeferencer, touching the 60 m
   buffer. Save the result as `<survey>_parcels_modified.gpkg` in the village folder and keep
   the `.points` file beside the GeoJSON. A village with no seed cannot start.
3. Make sure the two-reader neighbour transcription exists:
   `FMB_Vector\<village>\neighbour_transcription\nb_override_<village>.json`.

## The process (agreed 2026-09-20): grow from a seed, one ring at a time

1. **Seed.** Pick one parcel you have placed by hand (Kizhikaranai: 48A). It is fixed; its
   pose is read from your `.points` file. Nothing else of yours is used for the automation.
2. **Local satellite window.** The tool exports Google Satellite only around the seed and its
   ring (about 300 m across, 12 MB), never the whole village, so the tile reprojection error
   stays local.
3. **Ring 1.** The sheets print the survey numbers around each parcel and on which side. The
   parcels the seed's sheet names, and those whose sheets name the seed, are placed first:
   laid against the seed edge for edge, checked against the printed sides read from both
   sheets, tied to the seed's edge lines, adjusted, then topology-cleaned.
4. **Automatic ground control points.** Every corner of a placed parcel is written as a GCP
   in the QGIS Georeferencer format, matched to a satellite edge where one exists. Your points
   files are never read for the automated parcels.
5. **Review in your project.** A group is added to your open project: the seed in cyan, the
   automated parcels as red outlines at exactly your line width, the GCPs as yellow crosses,
   the satellite window switched off. Your layers are not touched, reloaded or restyled.
6. **Your verdict, then ring 2.** Parcels you accept become fixed for the next ring. The village
   grows outward until every sheet in the buffer is placed or flagged.

Command for one step (a work copy is made under `_logs\seed_<village>_<seed>_ring<n>_<date>`):

```
python -c "from autogeoref import evaluate; evaluate.seed_run('35_04_077', ['48A'], rings=1)"
```

## Run (whole village, once it has enough accepted parcels)

```
cd D:\code\FMB_to_GeoJSON
python -m autogeoref.georef_village 35_04_077 --raster --topology
```

| Flag | What it does |
|---|---|
| `--raster` | exports Google Satellite for the village once as a GeoTIFF and pins it. Needed the first time. |
| `--topology` | closes sub-metre gaps and overlaps between the new parcels. Your parcels are never edited. |
| `--no-review` | skips loading the review group into QGIS. |
| `--project <file.qgz>` | only load the review group if this project is the one open. |
| `--tracker` | adds three columns after Remarks in the tracker: Auto colour, Auto note, Auto run. Off unless asked. Status and formatting are never touched. |
| `--refit 46B,47A` | re-place those parcels from the `.points` files you edited. |
| `--reset-gcp 46B` | forget the tool's record of your GCPs for one parcel. |

A village of 15 parcels takes about 10 minutes. The first run with `--raster` adds one
download of roughly 65 MB.

## What comes out

In `FMB_Vector\<village>\`:

- `<survey>_parcels_modified.gpkg` for every parcel placed, with a `parcels` layer and an
  `edges` layer (length, grid bearing, true bearing per edge). FMB areas and lengths are
  carried as attributes next to the placed ones.
- `<survey>_parcels.geojson.points` in the QGIS Georeferencer format, so you can open the
  sheet in the Georeferencer and drag the points.
- `georef_status.csv` with one row per parcel: colour, method, residual, neighbours used,
  satellite share, notes.
- `anchors.csv` listing your hand-placed parcels and how the tool read their pose.
- `rail_conflicts.csv` where a parcel overlaps railway land; these are left for you.

In `FMB_Vector\worklist.csv`: the whole corridor, one row per parcel, newest verdicts.

In QGIS: a group `Georef review - <village>` with one layer per placed parcel coloured by
verdict.

## The colours

| Colour | Meaning | What to do |
|---|---|---|
| Green | fits three or more placed neighbours and no other pose comes close, or fits two neighbours with a clear satellite match | spot-check |
| Amber | fits one or two neighbours, or three with another pose nearly as good | look at it |
| Red | one neighbour and nothing else, contradicts its neighbours, or lands away from two parcels its own sheet names | fix by hand |
| Waiting | no placed neighbour shares a boundary with it yet | place a neighbour first or seed it |

On the 15 hand-placed parcels of Kizhikaranai, hiding each one in turn and letting the tool
place it from the other 14: 10 landed within 3 m of your placement, 13 within 5 m, both greens
were under 1.6 m, and the one red was the 16 m miss. Your own placements disagree with each
other by 1 to 3 m, so treat anything within 3 m as agreement.

Three rules override everything else:

- A parcel whose boundary residual against its neighbours exceeds 5 m is red.
- A parcel that lands more than 10 m from a placed parcel its own sheet names as a neighbour
  is red. This uses the two-reader transcription of the printed survey numbers.
- When the two best candidate poses explain nearly the same length of shared boundary (the
  runner-up within 80 % of the best), the neighbours have not really chosen, and the parcel
  is at most amber unless clear satellite edges break the tie. This is what happens when a
  parcel slides along a long straight boundary such as the railway strip.

## Fixing a parcel

1. Open the sheet GeoJSON in the QGIS Georeferencer. It loads the `.points` file.
2. Drag the wrong points onto the right ground positions. Keep at least two.
3. Save the points file, then run `--refit <survey>`. The tool re-places the sheet rigidly
   from your points and marks it as team-fixed.
4. Green parcels and refitted parcels become seeds for the next run, so the village grows
   outward from your work.

## Data cleaning before anything is shown

Every parcel the tool edits in the topology step is cleaned: gap strips are built with straight
mitred joins and limited to the space between the two parcels, hairline gaps are closed so a
filled strip merges into its plot (no double lines), needles thinner than 0.2 m are removed,
each plot is one ring, and collinear vertices left by clipping are dropped. Plots the topology
step did not touch keep the sheet's exact vertices. A railway strip keeps its shape; a parcel
the tool placed is clipped to it and the conflict is listed in `rail_conflicts.csv`.

## When your manual parcel and the automation disagree

The automation keeps every printed length. If your hand placement was stretched to fit (47B on
2026-09-20: 124 m along the strip where the sheet prints 109 m, area 7056 m² for a 6243 m²
parcel), the tool's version is the FMB parcel and the report's "distance from your placement"
is measuring your stretch, not an error of the tool. Decide per parcel which one stands.

## What the tool will not do

- Move, rewrite or delete any `_parcels_modified*.gpkg` you made. A newer hand file for a
  survey supersedes an older one; the older one is moved to `_logs`, never deleted.
- Stretch a sheet. Scale is always exactly 1.
- Use Puvi for position. Puvi appears only as a reference distance in the report.
- Write the tracker's Status column or change its styling.

## How the tool decides (for the curious)

1. Reads your hand-placed parcels and recovers each one's pose from its `.points` file, or
   from the geometry when the points fit worse than 3 m.
2. For each remaining sheet, walks its outline against every placed neighbour's outline,
   comparing edge lengths and turn angles. The walk stops at the corner where the two rings
   turn different ways, which is where the shared boundary ends. Every matching run of edges
   gives a candidate pose.
3. Drops candidates that sit on top of a placed parcel (after shaving 3 m off both, because
   hand placements disagree with each other by about that much), then ranks the rest by how
   much shared boundary they explain. The 180 degree twin of a rectangle loses here, and a
   pose that does not reach a parcel the sheet names as its neighbour ranks last.
4. Offers the best three to the satellite edge search. A pose is accepted from imagery only
   when it beats the runner-up clearly and the parcel has edges in two directions.
5. Adjusts all new parcels together, your parcels held fixed, so shared boundaries close.
6. Writes the files, the GCP points, the status CSV and the review group.

Calibration record: see the leave-one-out table in
`docs/superpowers/specs/2026-09-18-fmb-georef-gcp-sop-design.md`, section 3.
