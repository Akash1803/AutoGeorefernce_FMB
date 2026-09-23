# AutoGeorefernce_FMB

Tamil Nadu Field Measurement Book (FMB) sheets, from the Collabland portal PDF to a
georeferenced parcel on the ground, with every printed dimension kept exactly.

Built for the Tambaram to Chengalpattu railway corridor at FarmwiseAI (2026). Two parts:

1. **Sheet to vector** (`fmb_to_geojson.py`): reads the portal's FMB PDF, recovers the drawn
   polygons and the printed dimensions, and writes `<survey>_parcels.geojson` in the sheet's
   own metres.
2. **Vector to ground** (`autogeoref/`): places each sheet by laying it against a neighbour
   the team has already placed by hand, matching the shared boundary edge for edge. Rotation
   and shift only, scale exactly 1. Google Satellite confirms; Puvi is reference only. Every
   parcel comes back green, amber or red.

The team-facing procedure is in [docs/autogeoref-sop.md](docs/autogeoref-sop.md). The design,
its revision history and the calibration record are in
[docs/superpowers/specs/2026-09-18-fmb-georef-gcp-sop-design.md](docs/superpowers/specs/2026-09-18-fmb-georef-gcp-sop-design.md).

## Run

```
python -m autogeoref.georef_village 35_04_077 --raster --topology
```

Needs Python 3.12 with geopandas, shapely, rasterio, scipy, opencv, pyproj, openpyxl.
QGIS 3.40 is reached over the qgis_mcp socket (localhost:9876) when it is open; the tool
works without it. Project data lives outside the repo under `D:\Projects\<corridor>`
(`autogeoref/paths.py`). The portal fetch scripts read the session cookie from
`FMB_JSESSIONID`.

Sheet QC before placing (see [docs/sheet-qc.md](docs/sheet-qc.md)):

```
python -m autogeoref.georef_village 35_04_074 --sheet-qc
```

Neighbour transcription (PR 1, two readers) and the seed plan:

```
python -m autogeoref.georef_village 35_04_074 --render-sheets   # crops for the readers
python -m autogeoref.georef_village 35_04_074 --transcribe      # merge A.json + B.json -> table, review CSV, report, seeds
```

See [docs/pr1-transcription.md](docs/pr1-transcription.md). Every run appends evaluation rows
(`--no-rows` to skip).

## Tests

```
python -m pytest -m "not slow"      # 123 unit tests, about 90 s
python -m pytest -m slow            # leave-one-out and acceptance on the real village data
```

## Result on the calibration village (Kizhikaranai, 15 hand-placed parcels, 2026-09-19)

Each parcel hidden in turn and re-placed from the other 14: 10 within 3 m of the team's
placement, 13 within 5 m, all three railway strips oriented correctly. Greens both under
1.6 m, the one red the 16 m miss, no false green. The team's own placements disagree with each
other by 1 to 3 m, which sets the floor.

## Status

Research and development stage. Corridor rollout (19 villages, 656 sheets in the 60 m buffer)
waits on approval of the Kizhikaranai result and on one or two hand-placed seed parcels per
village.
