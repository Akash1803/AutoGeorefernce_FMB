# FMB Georeferencing SOP v2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Place Tamil Nadu FMB cadastral sketch polygons onto the ground automatically, keeping the team's hand-georeferenced parcels frozen, deriving position evidence from a satellite GeoTIFF and the sketches' own shared boundaries, and handing anything uncertain back to the team as a QGIS Georeferencer GCP file.

**Architecture:** A Python package `autogeoref` with one module per job. Evidence is ranked: frozen-anchor boundaries first, then boundaries with already-placed parcels plus the printed neighbour numbers, then satellite image edges (which may only refine, or decide alone when unambiguous), then the Puvi village vector (search window and start pose only). Every placement is a rotation plus a shift with scale fixed at 1, so FMB dimensions survive exactly. A whole-block least-squares adjustment settles the free parcels with anchors held fixed. QGIS 3.40 is reached through a length-prefixed JSON socket for review layers and for releasing file locks.

**Tech Stack:** Python 3.12.10; numpy 2.5.1, scipy 1.18.0, shapely 2.1.2, geopandas 1.1.4, fiona 1.10.1, pyproj 3.7.2, rasterio 1.5.0, opencv-python 5.0.0, openpyxl 3.1.5, pytest 9.1.1. QGIS 3.40.10 (PyQGIS reached over a socket, never imported in-process). No GDAL Python bindings (`osgeo`) are installed — use rasterio and fiona.

**Spec:** `D:\code\FMB_to_GeoJSON\docs\superpowers\specs\2026-09-18-fmb-georef-gcp-sop-design.md`

## Background for the implementer

A **survey** (e.g. `47B`, `169`, `526A`) is a Tamil Nadu cadastral land unit. Its **FMB sheet** is a scanned/vector PDF held by the state portal; `fmb_to_geojson.py` (already written, not part of this plan) converts one sheet into `FMB_Vector\<village>\<survey>_parcels.geojson` — polygons in **sheet metres**, an arbitrary local coordinate system whose origin is the sheet corner. Georeferencing means finding the rotation and shift that puts those polygons on the ground in **EPSG:32644** (UTM zone 44N).

A **village** is identified by a code like `35_04_077` (district 35, taluk 04, village 077). One village folder holds one GeoJSON per survey.

The team georeferences by hand in the QGIS Georeferencer, which writes two things: the placed polygons as `<survey>_parcels_modified.gpkg`, and the control points as `<survey>_parcels.geojson.points`. **Both are read-only inputs to this tool.** Their placements are affine (measured scales 0.83–1.15), so they are used as a *pose*, never as geometry to copy.

**Puvi** is a bulk cadastral vector archive under `D:\Data\Puvi_data\vector\<dd>\<tt>\<vvv>\vector\*.shp` (EPSG:4326). It is 5–37 m off the truth here, so the user's standing rule is: reference only — it may supply a search window and a starting guess, and its distance is reported, but it never accepts, rejects or moves a placement.

## Global Constraints

- CRS for every output and every computation: **EPSG:32644**. Puvi shapefiles are EPSG:4326 and must be reprojected on read.
- **Rigid only:** every placement is `rotate(theta) + translate(t)` with scale fixed at 1. After each placement assert `abs(area_after - area_before) < 1e-6` and `abs(length_after - length_before) < 1e-6`, and fail that survey loudly rather than write it.
- **Never modify** any `*_parcels_modified*.gpkg` or `*.points` file the tool did not write, and never delete any file: supersede by moving to `D:\Projects\Tambaram_Chengalpattu\_logs\georef_superseded_<run_id>\`.
- A `_modified*` file is **manual** if and only if its geometry fingerprint is not recorded as tool-written in `georef_status.csv`. Files with no record at all, and every `_AUTOv2` / `_AUTODEMO` output of the 2026-09-17 prototype, count as manual.
- **One anchor per survey:** the newest `_modified*` file by modification time, unless `anchors.csv` names another.
- **Puvi is reference only.** It may set a search window and a start pose. It must never gate, reject, or move a placement. `puvi_reference_m` is reported in the status CSV.
- Project root: `D:\Projects\Tambaram_Chengalpattu`. Code home: `D:\code\FMB_to_GeoJSON`. Package: `D:\code\FMB_to_GeoJSON\autogeoref`. Tests: `D:\code\FMB_to_GeoJSON\tests`.
- Angles in the code are **radians internally, degrees in every file and message**. Bearings are clockwise from grid north, 0–360.
- Tolerances (spec §5): anchor boundary sigma `max(0.30 m, anchor.rms)`; auto–auto boundary sigma 0.30 m; image evidence sigma 1.0 m; start-pose prior sigma 10 m; team GCP sigma 0.20 m; image may move an anchored pose by at most 1.5 m and 2°.
- Reference paths used by tests: village `35_04_077` (Kizhikaranai, 15 manual parcels), village `35_04_052` (Thailavaram, 13 manual parcels in `C:\Users\FAI-Akash\Downloads\thailavaram.geojson`), spike raster `D:\code\FMB_to_GeoJSON\autogeoref\scratch_2026-09-17\spike_gcp\kizhi_sat_z20.tif`.
- The 2026-09-17 prototypes in `autogeoref\scratch_2026-09-17\` are the reference implementations to port from. They are **not** imported by the new package; port the logic and delete nothing.

---

## File structure

| File | Responsibility |
|---|---|
| `autogeoref/__init__.py` | version, package marker |
| `autogeoref/paths.py` | every filesystem path in one place; village/survey helpers |
| `autogeoref/sheets.py` | load sheet polygons; outline (no collinear merge); edge lengths; turn angles |
| `autogeoref/fit.py` | rigid fit (scale 1), pose application, block adjustment with fixed anchors |
| `autogeoref/anchors.py` | choose/fingerprint manual files; read `.points`; recover anchor poses; `anchors.csv` |
| `autogeoref/raster.py` | export and pin the satellite GeoTIFF |
| `autogeoref/window.py` | Puvi search window and start poses (with the 180° twin) |
| `autogeoref/neighbours.py` | adjacency and sides from the two-reader transcription |
| `autogeoref/match.py` | shared-boundary chains; point-to-line and vertex-pair observations |
| `autogeoref/edges.py` | satellite line segments with direction; spatial index |
| `autogeoref/align.py` | pose search, share score, margin test, observability test |
| `autogeoref/gcp.py` | corner GCPs, `.points` writer/reader, refit from team GCPs |
| `autogeoref/qgis_bridge.py` | length-prefixed JSON socket to QGIS; `execute()` |
| `autogeoref/files.py` | held-layer detection, temp-write-and-rename, WAL cleanup, supersede |
| `autogeoref/topology.py` | snap/clip/fill between free parcels only; railway conflicts listed |
| `autogeoref/review.py` | status CSV, worklist, tracker columns, QGIS review group |
| `autogeoref/engine.py` | the per-village pipeline that calls everything in spec order |
| `autogeoref/georef_village.py` | command-line entry point |
| `tests/…` | one test module per source module |

---

## Task 1: Package skeleton, paths, sheet geometry

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\__init__.py`
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\paths.py`
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\sheets.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_paths.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_sheets.py`
- Create: `D:\code\FMB_to_GeoJSON\pytest.ini`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `paths.PROJECT: Path`, `paths.CODE: Path`
  - `paths.vector_dir(village: str) -> Path`, `paths.sketch_dir(village) -> Path`, `paths.georef_dir(village) -> Path`, `paths.logs_dir() -> Path`
  - `paths.sheet_path(village, survey) -> Path`, `paths.manual_files(village, survey) -> list[Path]`, `paths.points_path(village, survey) -> Path`, `paths.output_path(village, survey) -> Path`
  - `paths.surveys_with_sheets(village) -> list[str]`, `paths.survey_sort_key(survey) -> tuple`
  - `sheets.load_sheet(village, survey) -> list[tuple[dict, Polygon]]`
  - `sheets.dissolve(polys) -> Polygon`
  - `sheets.outline(polys) -> list[tuple[float, float]]` (clockwise, duplicates removed, **no collinear merge**)
  - `sheets.edge_lengths(v) -> list[float]`, `sheets.turn_angles(v) -> list[float]` (degrees)
  - `sheets.sample_outline(v, step=1.0) -> tuple[np.ndarray, np.ndarray]` (points, bearings in radians)

- [ ] **Step 1: Put the code home under version control and create the package skeleton**

`D:\code\FMB_to_GeoJSON` is not a git repository yet. The plan assumes frequent commits, so initialise it. (If the user objects to a repo here, skip this step and drop every later commit step.)

```bash
cd /d/code/FMB_to_GeoJSON
git init
printf '__pycache__/\n*.pyc\n.pytest_cache/\nautogeoref/scratch_2026-09-17/__pycache__/\n' > .gitignore
git add -A
git commit -m "chore: put FMB_to_GeoJSON under version control"
mkdir -p autogeoref tests
printf '__version__ = "2.0.0"\n' > autogeoref/__init__.py
printf '[pytest]\ntestpaths = tests\naddopts = -q\n' > pytest.ini
```

- [ ] **Step 2: Write the failing tests for paths and sheet geometry**

Create `tests/test_paths.py`:

```python
from pathlib import Path
from autogeoref import paths


def test_village_directories_point_at_the_project():
    assert paths.vector_dir("35_04_077") == paths.PROJECT / "FMB_Vector" / "35_04_077"
    assert paths.sketch_dir("35_04_077") == paths.PROJECT / "FMB_Sketches" / "35_04_077"
    assert paths.georef_dir("35_04_077") == paths.PROJECT / "FMB_Georef" / "35_04_077"


def test_sheet_and_output_paths():
    assert paths.sheet_path("35_04_077", "47B").name == "47B_parcels.geojson"
    assert paths.points_path("35_04_077", "47B").name == "47B_parcels.geojson.points"
    assert paths.output_path("35_04_077", "47B").name == "47B_parcels_modified.gpkg"


def test_manual_files_are_newest_first():
    files = paths.manual_files("35_04_077", "47A")
    assert [f.name for f in files][:1] == ["47A_parcels_modified2.gpkg"]
    assert len(files) == 3


def test_surveys_sort_numerically_then_by_letter():
    got = sorted(["169", "40B", "42A", "9", "40A"], key=paths.survey_sort_key)
    assert got == ["9", "40A", "40B", "42A", "169"]


def test_surveys_with_sheets_finds_the_village():
    got = paths.surveys_with_sheets("35_04_077")
    assert "47B" in got and "169" in got
    assert all(paths.sheet_path("35_04_077", s).exists() for s in got)
```

Create `tests/test_sheets.py`:

```python
import math
import numpy as np
from shapely.geometry import Polygon
from autogeoref import sheets


def test_outline_is_clockwise_and_keeps_every_vertex():
    # a square with an extra vertex in the middle of the north edge (nearly collinear)
    poly = Polygon([(0, 0), (10, 0), (10, 10), (5, 10.02), (0, 10)])
    v = sheets.outline([(None, poly)])
    assert len(v) == 5, "no collinear merge: the 10.02 vertex must survive"
    area = 0.5 * sum(v[i][0] * v[(i + 1) % 5][1] - v[(i + 1) % 5][0] * v[i][1] for i in range(5))
    assert area < 0, "clockwise rings have negative shoelace area"


def test_edge_lengths_and_turn_angles_on_a_square():
    v = sheets.outline([(None, Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]))])
    assert [round(x, 6) for x in sheets.edge_lengths(v)] == [10.0] * 4
    assert [round(x, 3) for x in sheets.turn_angles(v)] == [90.0] * 4


def test_sample_outline_spacing_and_bearings():
    v = sheets.outline([(None, Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]))])
    pts, brg = sheets.sample_outline(v, step=1.0)
    assert len(pts) == 40 and len(brg) == 40
    assert np.allclose(np.linalg.norm(np.diff(pts[:10], axis=0), axis=1), 1.0)
    assert len({round(float(b) % math.pi, 6) for b in brg}) == 2, "a square has two edge directions"


def test_load_sheet_returns_properties_and_polygons():
    got = sheets.load_sheet("35_04_077", "47B")
    assert len(got) == 3
    props, geom = got[0]
    assert props["survey_no"] == "47B"
    assert geom.is_valid and geom.area > 0


def test_dissolve_merges_touching_polygons():
    a = Polygon([(0, 0), (5, 0), (5, 5), (0, 5)])
    b = Polygon([(5, 0), (10, 0), (10, 5), (5, 5)])
    assert sheets.dissolve([(None, a), (None, b)]).area == 50
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_paths.py tests/test_sheets.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.paths'`

- [ ] **Step 4: Write `autogeoref/paths.py`**

```python
"""Every filesystem path the tool uses, in one place."""
import re
from pathlib import Path

PROJECT = Path(r"D:\Projects\Tambaram_Chengalpattu")
CODE = Path(__file__).resolve().parent.parent
PUVI = Path(r"D:\Data\Puvi_data\vector")
BUFFER_GPKG = PROJECT / "Railway_Buffer_Vector_Plots.gpkg"
BUFFER_LAYER = "vector_in_buffer_30m"
RAIL_GPKG = PROJECT / "Tambaram_Chengalpattu_Railway.gpkg"
RAIL_LAYER = "rail_line"
TRACKER = PROJECT / "Georeferencing_Tracker.xlsx"


def vector_dir(village): return PROJECT / "FMB_Vector" / village
def sketch_dir(village): return PROJECT / "FMB_Sketches" / village
def georef_dir(village): return PROJECT / "FMB_Georef" / village
def logs_dir(): return PROJECT / "_logs"
def sheet_path(village, survey): return vector_dir(village) / ("%s_parcels.geojson" % survey)
def points_path(village, survey): return vector_dir(village) / ("%s_parcels.geojson.points" % survey)
def output_path(village, survey): return vector_dir(village) / ("%s_parcels_modified.gpkg" % survey)
def status_path(village): return vector_dir(village) / "georef_status.csv"
def anchors_path(village): return vector_dir(village) / "anchors.csv"
def gcp_csv(village, survey): return vector_dir(village) / "gcp" / ("%s_gcp.csv" % survey)
def raster_pin(village): return georef_dir(village) / "raster.json"


def puvi_shapefile(village):
    """D:\\Data\\Puvi_data\\vector\\35\\04\\077\\vector\\*.shp for village 35_04_077."""
    d, t, v = village.split("_")
    hits = sorted((PUVI / d / t / v / "vector").glob("*.shp"))
    return hits[0] if hits else None


def survey_sort_key(survey):
    """'40A' -> (40, '40A') so 9 < 40A < 40B < 169."""
    m = re.match(r"\d+", str(survey))
    return (int(m.group()) if m else 10 ** 9, str(survey))


def manual_files(village, survey):
    """Every hand-placed GeoPackage for this survey, newest modification time first."""
    hits = list(vector_dir(village).glob("%s_parcels_modified*.gpkg" % survey))
    return sorted(hits, key=lambda p: p.stat().st_mtime, reverse=True)


def surveys_with_sheets(village):
    """Survey ids that have a sheet GeoJSON, numerically sorted."""
    out = [p.name[: -len("_parcels.geojson")] for p in vector_dir(village).glob("*_parcels.geojson")]
    return sorted(out, key=survey_sort_key)
```

- [ ] **Step 5: Write `autogeoref/sheets.py`**

The 2026-09-17 prototype merged vertices whose turn angle was under 1°. That merge moved vertices differently on the two sheets of a pair and broke boundary matching (spec §6.1), so this version keeps every vertex and only drops exact duplicates.

```python
"""Sheet polygons in sheet metres, and the outline geometry derived from them."""
import json
import math
import numpy as np
from shapely.geometry import shape, Polygon
from shapely.ops import unary_union
from . import paths


def load_sheet(village, survey):
    """[(properties, Polygon)] for one sheet, in sheet metres."""
    g = json.loads(paths.sheet_path(village, survey).read_text(encoding="utf-8"))
    return [(f["properties"], shape(f["geometry"])) for f in g["features"]]


def dissolve(polys):
    """Largest single polygon of the union of a sheet's parts."""
    u = unary_union([p.buffer(0) for _, p in polys])
    return max(u.geoms, key=lambda q: q.area) if u.geom_type == "MultiPolygon" else u


def outline(polys):
    """Exterior ring of the dissolved sheet, clockwise, every real vertex kept."""
    ring = list(dissolve(polys).exterior.coords)[:-1]
    out = [ring[0]]
    for p in ring[1:]:
        if math.dist(p, out[-1]) > 1e-9:
            out.append(p)
    if math.dist(out[0], out[-1]) <= 1e-9:
        out.pop()
    if not Polygon(out).exterior.is_ccw:
        return out
    return out[::-1]


def edge_lengths(v):
    return [math.dist(v[i], v[(i + 1) % len(v)]) for i in range(len(v))]


def turn_angles(v):
    """Turn angle in degrees at each vertex of a closed ring (0 = straight)."""
    n = len(v)
    out = []
    for i in range(n):
        a, b, c = v[i - 1], v[i], v[(i + 1) % n]
        t = math.degrees(math.atan2(c[1] - b[1], c[0] - b[0]) - math.atan2(b[1] - a[1], b[0] - a[0]))
        out.append(abs((t + 180) % 360 - 180))
    return out


def sample_outline(v, step=1.0):
    """Points every `step` metres along the ring, with each sample's edge bearing in radians."""
    pts, brg = [], []
    for i in range(len(v)):
        a = np.asarray(v[i], float)
        b = np.asarray(v[(i + 1) % len(v)], float)
        L = float(np.linalg.norm(b - a))
        n = max(1, int(round(L / step)))
        ang = math.atan2(b[1] - a[1], b[0] - a[0])
        for k in range(n):
            pts.append(a + (b - a) * (k + 0.5) / n)
            brg.append(ang)
    return np.array(pts), np.array(brg)
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_paths.py tests/test_sheets.py -v`
Expected: PASS, 10 tests.

- [ ] **Step 7: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/__init__.py autogeoref/paths.py autogeoref/sheets.py tests/test_paths.py tests/test_sheets.py pytest.ini .gitignore
git commit -m "feat(autogeoref): package skeleton, project paths and sheet outline geometry"
```

---

## Task 2: Rigid fit and pose application

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\fit.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_fit.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `fit.rigid_fit(P, Q) -> (theta_deg: float, t: np.ndarray(2), rms: float, max_res: float)` — Procrustes with scale fixed at 1; `P` and `Q` are `(n, 2)` arrays.
  - `fit.similarity_fit(P, Q) -> (theta_deg, t, rms, scale, anisotropy)` — used only to describe how non-rigid a hand file is.
  - `fit.apply_pose(geom, theta_deg, t) -> geom` — rotate about the origin then translate.
  - `fit.transform_points(P, theta_deg, t) -> np.ndarray`
  - `fit.assert_rigid(before, after)` — raises `AssertionError` when area or perimeter moved by more than 1e-6.

- [ ] **Step 1: Write the failing test**

Create `tests/test_fit.py`:

```python
import math
import numpy as np
import pytest
from shapely.geometry import Polygon
from shapely import affinity
from autogeoref import fit

SQUARE = Polygon([(0, 0), (10, 0), (10, 6), (0, 6)])


def test_rigid_fit_recovers_a_known_rotation_and_shift():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = fit.transform_points(P, 33.0, np.array([150.0, -70.0]))
    theta, t, rms, mx = fit.rigid_fit(P, Q)
    assert theta == pytest.approx(33.0, abs=1e-6)
    assert t == pytest.approx(np.array([150.0, -70.0]), abs=1e-6)
    assert rms < 1e-9 and mx < 1e-9


def test_rigid_fit_never_scales_even_when_the_target_is_stretched():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = P * 1.10
    theta, t, rms, mx = fit.rigid_fit(P, Q)
    moved = fit.transform_points(P, theta, t)
    side = float(np.linalg.norm(moved[1] - moved[0]))
    assert side == pytest.approx(10.0, abs=1e-9), "scale must stay 1"
    assert rms > 0.1, "the stretch shows up as residual, not as scale"


def test_rigid_fit_rejects_a_mirror():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = P.copy()
    Q[:, 1] *= -1
    theta, t, rms, mx = fit.rigid_fit(P, Q)
    moved = fit.transform_points(P, theta, t)
    assert np.cross(moved[1] - moved[0], moved[2] - moved[1]) > 0, "handedness preserved"


def test_similarity_fit_reports_scale_and_anisotropy():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = P * 0.90
    theta, t, rms, scale, aniso = fit.similarity_fit(P, Q)
    assert scale == pytest.approx(0.90, abs=1e-6)
    assert aniso == pytest.approx(1.0, abs=1e-6)


def test_apply_pose_preserves_area_and_perimeter():
    moved = fit.apply_pose(SQUARE, 47.5, np.array([1000.0, 2000.0]))
    assert moved.area == pytest.approx(SQUARE.area, abs=1e-9)
    assert moved.length == pytest.approx(SQUARE.length, abs=1e-9)
    fit.assert_rigid(SQUARE, moved)


def test_assert_rigid_catches_a_stretch():
    with pytest.raises(AssertionError):
        fit.assert_rigid(SQUARE, affinity.scale(SQUARE, 1.01, 1.01))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_fit.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.fit'`

- [ ] **Step 3: Write `autogeoref/fit.py` (rigid part only; the block adjustment arrives in Task 10)**

```python
"""Rigid (rotation + shift, scale 1) fitting and pose application."""
import math
import numpy as np
from shapely import affinity


def _procrustes(P, Q):
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:                      # reflection: flip the smallest singular direction
        Vt = Vt.copy()
        Vt[-1] *= -1
        R = Vt.T @ U.T
    return R, pc, qc, S, P, Q


def rigid_fit(P, Q):
    """Rotation + shift with scale fixed at 1. Returns (theta_deg, t, rms, max_residual)."""
    R, pc, qc, _S, P, Q = _procrustes(P, Q)
    t = qc - R @ pc
    res = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
    theta = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    return theta, t, float(np.sqrt((res ** 2).mean())), float(res.max())


def similarity_fit(P, Q):
    """Descriptive only: how much scale and anisotropy a hand placement carries."""
    R, pc, qc, S, P, Q = _procrustes(P, Q)
    denom = ((P - pc) ** 2).sum()
    scale = float(S.sum() / denom) if denom else float("nan")
    t = qc - scale * (R @ pc)
    res = np.linalg.norm(scale * (R @ P.T).T + t - Q, axis=1)
    A = np.hstack([P, np.ones((len(P), 1))])
    M, *_ = np.linalg.lstsq(A, Q, rcond=None)
    sv = np.linalg.svd(M[:2].T, compute_uv=False)
    aniso = float(sv.max() / sv.min()) if sv.min() > 0 else float("inf")
    theta = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    return theta, t, float(np.sqrt((res ** 2).mean())), scale, aniso


def transform_points(P, theta_deg, t):
    th = math.radians(theta_deg)
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    return (R @ np.asarray(P, float).T).T + np.asarray(t, float)


def apply_pose(geom, theta_deg, t):
    return affinity.translate(affinity.rotate(geom, theta_deg, origin=(0, 0)), float(t[0]), float(t[1]))


def assert_rigid(before, after, tol=1e-6):
    if abs(after.area - before.area) > tol or abs(after.length - before.length) > tol:
        raise AssertionError(
            "placement is not rigid: area %.9f -> %.9f, perimeter %.9f -> %.9f"
            % (before.area, after.area, before.length, after.length))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_fit.py -v`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/fit.py tests/test_fit.py
git commit -m "feat(autogeoref): rigid Procrustes fit, similarity description and pose application"
```

---

## Task 3: Anchors — choose, fingerprint, read `.points`, recover poses

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\anchors.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_anchors.py`

**Interfaces:**
- Consumes: `paths.manual_files`, `paths.points_path`, `sheets.load_sheet`, `sheets.dissolve`, `fit.rigid_fit`, `fit.similarity_fit`, `fit.transform_points`.
- Produces:
  - `anchors.Anchor` dataclass with fields `survey, file, theta, t, rms, scale, anisotropy, source, fingerprint, n_points`.
  - `anchors.read_points(path) -> np.ndarray` shape `(n, 4)`, columns `mapX, mapY, sourceX, sourceY`, disabled rows dropped.
  - `anchors.write_points(path, rows, crs_wkt)` — the QGIS Georeferencer format.
  - `anchors.fingerprint(gpkg_path) -> str` — sha256 over rounded WKB of every geometry, read through fiona so a WAL sidecar is included.
  - `anchors.pose_from_points(village, survey) -> tuple | None` -> `(theta, t, rms, scale, aniso, n)`
  - `anchors.pose_from_geometry(village, survey, gpkg_path) -> tuple` -> `(theta, t, rms, scale, aniso, n)`
  - `anchors.load_anchors(village, tool_written: set[str]) -> dict[str, Anchor]`
  - `anchors.anchor_conflicts(village, anchors) -> list[dict]` — pairs whose shared boundary disagrees.
  - `anchors.write_anchors_csv(village, anchors, superseded, conflicts)`

- [ ] **Step 1: Write the failing test**

Create `tests/test_anchors.py`:

```python
import numpy as np
import pytest
from autogeoref import anchors, paths

VILLAGE = "35_04_077"


def test_read_points_parses_the_qgis_georeferencer_file():
    P = anchors.read_points(paths.points_path(VILLAGE, "47B"))
    assert P.shape[1] == 4 and len(P) >= 3
    assert 380000 < P[0, 0] < 400000 and 1400000 < P[0, 1] < 1420000, "map side is EPSG:32644"
    assert abs(P[0, 2]) < 5000 and abs(P[0, 3]) < 5000, "source side is sheet metres"


def test_write_points_round_trips(tmp_path):
    rows = [(392809.1, 1413357.1, 198.4, 399.9), (392900.0, 1413400.0, 250.0, 420.0)]
    p = tmp_path / "x_parcels.geojson.points"
    anchors.write_points(p, rows, crs_wkt='PROJCRS["WGS 84 / UTM zone 44N"]')
    text = p.read_text(encoding="utf-8").splitlines()
    assert text[0].startswith("#CRS: ")
    assert text[1] == "mapX,mapY,sourceX,sourceY,enable,dX,dY,residual"
    back = anchors.read_points(p)
    assert np.allclose(back, np.array(rows), atol=1e-6)


def test_pose_from_points_is_rigid_and_matches_the_hand_file_roughly():
    got = anchors.pose_from_points(VILLAGE, "46A")
    assert got is not None
    theta, t, rms, scale, aniso, n = got
    assert n >= 3
    assert rms < 1.0, "46A was placed cleanly; its own GCPs fit a rigid pose"
    assert 0.8 < scale < 1.2


def test_load_anchors_takes_one_file_per_survey_and_flags_the_stretch():
    got = anchors.load_anchors(VILLAGE, tool_written=set())
    assert len(got) == 15, "all 15 Kizhikaranai surveys are hand placed"
    assert got["47A"].file.name == "47A_parcels_modified2.gpkg", "newest variant wins"
    assert got["48B"].scale < 0.95, "48B is georeferenced ~10 % small"
    assert all(a.source in ("points", "geometry") for a in got.values())


def test_fingerprint_is_stable_and_geometry_sensitive(tmp_path):
    import geopandas as gpd
    from shapely.geometry import Polygon
    a = tmp_path / "a.gpkg"
    gpd.GeoDataFrame({"x": [1]}, geometry=[Polygon([(0, 0), (1, 0), (1, 1)])], crs="EPSG:32644").to_file(a, driver="GPKG")
    f1 = anchors.fingerprint(a)
    assert f1 == anchors.fingerprint(a)
    b = tmp_path / "b.gpkg"
    gpd.GeoDataFrame({"x": [1]}, geometry=[Polygon([(0, 0), (1, 0), (1, 1.5)])], crs="EPSG:32644").to_file(b, driver="GPKG")
    assert anchors.fingerprint(b) != f1


def test_anchor_conflicts_reports_the_known_48A_171_disagreement():
    got = anchors.load_anchors(VILLAGE, tool_written=set())
    conf = anchors.anchor_conflicts(VILLAGE, got)
    pairs = {tuple(sorted((c["a"], c["b"]))) for c in conf}
    assert ("171", "48A") in pairs
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_anchors.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.anchors'`

- [ ] **Step 3: Write `autogeoref/anchors.py`**

```python
"""Manual (hand-georeferenced) parcels: choose one file per survey, fingerprint it, and
recover the rigid pose it implies. The files themselves are never modified."""
import csv
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

import fiona
import numpy as np
from shapely.geometry import shape
from shapely.ops import unary_union

from . import fit, paths, sheets

POINTS_HEADER = "mapX,mapY,sourceX,sourceY,enable,dX,dY,residual"
POSE_RMS_LIMIT = 3.0        # above this the .points fit is not trusted; use the geometry instead
CONFLICT_SAMPLE = 1.0       # metres between boundary samples when testing two anchors


@dataclass
class Anchor:
    survey: str
    file: Path
    theta: float
    t: np.ndarray
    rms: float
    scale: float
    anisotropy: float
    source: str
    fingerprint: str
    n_points: int


def read_points(path):
    """(n, 4) array of mapX, mapY, sourceX, sourceY for the enabled rows of a QGIS .points file."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("mapX"):
            continue
        parts = line.split(",")
        if len(parts) < 5 or parts[4].strip() not in ("1", "true", "True"):
            continue
        rows.append([float(x) for x in parts[:4]])
    return np.array(rows, dtype=float).reshape(-1, 4)


def write_points(path, rows, crs_wkt):
    """Write GCPs in the QGIS Georeferencer format. rows = [(mapX, mapY, sourceX, sourceY), ...]."""
    out = ["#CRS: " + crs_wkt, POINTS_HEADER]
    for mx, my, sx, sy in rows:
        out.append("%.17f,%.17f,%.17f,%.17f,1,0,0,0" % (mx, my, sx, sy))
    Path(path).write_text("\n".join(out) + "\n", encoding="utf-8")


def fingerprint(gpkg_path):
    """sha256 of the layer's geometries, rounded to the millimetre, read through OGR."""
    h = hashlib.sha256()
    with fiona.open(str(gpkg_path)) as src:
        for feat in src:
            g = shape(feat["geometry"])
            h.update(g.wkb_hex.encode() if hasattr(g, "wkb_hex") else g.wkb)
    return h.hexdigest()


def _fit_pose(P_sheet, Q_map):
    theta, t, rms, _mx = fit.rigid_fit(P_sheet, Q_map)
    _th2, _t2, _rms2, scale, aniso = fit.similarity_fit(P_sheet, Q_map)
    return theta, t, rms, scale, aniso, len(P_sheet)


def pose_from_points(village, survey):
    """Rigid pose implied by the team's own control points, or None when there is no file."""
    p = paths.points_path(village, survey)
    if not p.exists():
        return None
    P = read_points(p)
    if len(P) < 2:
        return None
    return _fit_pose(P[:, 2:4], P[:, 0:2])


def pose_from_geometry(village, survey, gpkg_path):
    """Rigid pose implied by the hand-placed geometry: pair sheet vertices with hand vertices by poly_id."""
    import geopandas as gpd
    hand = gpd.read_file(gpkg_path)
    if hand.crs is not None and hand.crs.to_epsg() != 32644:
        hand = hand.to_crs(32644)
    sheet = sheets.load_sheet(village, survey)
    by_id = {}
    if "poly_id" in hand.columns:
        for row in hand.itertuples():
            try:
                by_id[int(row.poly_id)] = row.geometry
            except (TypeError, ValueError):
                pass
    P, Q = [], []
    for props, geom in sheet:
        other = by_id.get(int(props.get("poly_id", -1)))
        if other is None:
            continue
        ca = list(geom.exterior.coords)[:-1]
        cb = list((other.geoms[0] if other.geom_type == "MultiPolygon" else other).exterior.coords)[:-1]
        if len(ca) == len(cb):
            P += ca
            Q += cb
    if len(P) < 2:                                   # fall back to centroids of the largest parts
        A = sorted([g for _, g in sheet], key=lambda g: -g.area)
        B = sorted(list(hand.geometry), key=lambda g: -g.area)
        P = [list(g.centroid.coords)[0] for g in A[: min(len(A), len(B))]]
        Q = [list(g.centroid.coords)[0] for g in B[: min(len(A), len(B))]]
    return _fit_pose(np.array(P, float), np.array(Q, float))


def load_anchors(village, tool_written):
    """One Anchor per survey that has a manual file the tool did not write."""
    out = {}
    for survey in paths.surveys_with_sheets(village):
        for f in paths.manual_files(village, survey):       # newest first
            fp = fingerprint(f)
            if fp in tool_written:
                continue
            pose = pose_from_points(village, survey)
            source = "points"
            if pose is None or pose[2] > POSE_RMS_LIMIT:
                pose = pose_from_geometry(village, survey, f)
                source = "geometry"
            theta, t, rms, scale, aniso, n = pose
            out[survey] = Anchor(survey, f, theta, np.asarray(t, float), rms, scale, aniso, source, fp, n)
            break
    return out


def placed_geometry(village, anchor):
    """The rigid sheet at the anchor's recovered pose (NOT the hand geometry)."""
    sheet = sheets.load_sheet(village, anchor.survey)
    return unary_union([fit.apply_pose(g, anchor.theta, anchor.t).buffer(0) for _, g in sheet])


def anchor_conflicts(village, anchor_map, limit=0.30):
    """Anchor pairs whose shared boundary disagrees by more than limit + both rms."""
    geoms = {s: placed_geometry(village, a) for s, a in anchor_map.items()}
    keys = sorted(geoms, key=paths.survey_sort_key)
    out = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            ga, gb = geoms[a], geoms[b]
            if ga.distance(gb) > 5.0:
                continue
            ring = ga.exterior
            n = max(20, int(ring.length / CONFLICT_SAMPLE))
            d = []
            for k in range(n):
                p = ring.interpolate(k / n, normalized=True)
                if p.distance(gb) < 3.0:
                    d.append(p.distance(gb.exterior))
            if not d:
                continue
            worst = float(np.percentile(d, 90))
            allow = limit + anchor_map[a].rms + anchor_map[b].rms
            if worst > allow:
                out.append({"a": a, "b": b, "p90_gap_m": round(worst, 2), "allowed_m": round(allow, 2),
                            "overlap_sqm": round(ga.intersection(gb).area, 1)})
    return out


def write_anchors_csv(village, anchor_map, superseded, conflicts):
    p = paths.anchors_path(village)
    p.parent.mkdir(parents=True, exist_ok=True)
    cols = ["survey", "file", "heading_deg", "easting", "northing", "pose_source", "pose_rms_m",
            "implied_scale", "anisotropy", "n_control_points", "fingerprint", "superseded_files", "conflicts"]
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for s in sorted(anchor_map, key=paths.survey_sort_key):
            a = anchor_map[s]
            w.writerow({"survey": s, "file": a.file.name, "heading_deg": round(a.theta, 3),
                        "easting": round(float(a.t[0]), 3), "northing": round(float(a.t[1]), 3),
                        "pose_source": a.source, "pose_rms_m": round(a.rms, 3),
                        "implied_scale": round(a.scale, 4), "anisotropy": round(a.anisotropy, 4),
                        "n_control_points": a.n_points, "fingerprint": a.fingerprint[:16],
                        "superseded_files": ";".join(x.name for x in superseded.get(s, [])),
                        "conflicts": ";".join("%s (%.2f m)" % (c["b"] if c["a"] == s else c["a"], c["p90_gap_m"])
                                              for c in conflicts if s in (c["a"], c["b"]))})
    return p
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_anchors.py -v`
Expected: PASS, 6 tests. If `test_anchor_conflicts_reports_the_known_48A_171_disagreement` fails, print the computed `p90_gap_m` for that pair and check it against the 10.7 m recorded in the spec before changing the threshold.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/anchors.py tests/test_anchors.py
git commit -m "feat(autogeoref): frozen anchors with pose from the team's own GCP files"
```

---

## Task 4: Satellite raster export and pin

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\raster.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_raster.py`

**Interfaces:**
- Consumes: `paths.georef_dir`, `paths.raster_pin`, `paths.logs_dir`.
- Produces:
  - `raster.TILE_XML` — the GDAL_WMS/TMS description of the project's Google Satellite XYZ layer.
  - `raster.export(village, bounds, zoom=20, resolution=0.15, run_id=None) -> Path` — writes `FMB_Georef\<village>\satellite_z<zoom>_<YYYYMMDD>.tif` in EPSG:32644, archives any previous file to `_logs\`, rewrites the pin.
  - `raster.pinned(village) -> dict | None` — `{path, sha256, zoom, fetched, bounds}` from `raster.json`.
  - `raster.open_pinned(village) -> rasterio.DatasetReader`

- [ ] **Step 1: Write the failing test**

Create `tests/test_raster.py`:

```python
import json
import numpy as np
import pytest
import rasterio
from autogeoref import raster

SPIKE = r"D:\code\FMB_to_GeoJSON\autogeoref\scratch_2026-09-17\spike_gcp\kizhi_sat_z20.tif"


def test_tile_description_is_a_tms_service_for_the_project_layer():
    assert "mt1.google.com/vt/lyrs=s" in raster.TILE_XML
    assert "<TileLevel>20</TileLevel>" in raster.TILE_XML
    assert "EPSG:3857" in raster.TILE_XML


def test_the_spike_raster_is_readable_and_georeferenced():
    with rasterio.open(SPIKE) as src:
        assert src.crs.to_epsg() == 32644
        assert src.res[0] == pytest.approx(0.15, abs=0.01)
        assert src.count == 3


def test_export_writes_a_small_window_and_pins_it(tmp_path, monkeypatch):
    monkeypatch.setattr(raster.paths, "PROJECT", tmp_path)
    bounds = (392900.0, 1413300.0, 392980.0, 1413380.0)     # 80 x 80 m over Kizhikaranai
    out = raster.export("35_04_077", bounds, zoom=20, resolution=0.15)
    assert out.exists()
    with rasterio.open(out) as src:
        assert src.crs.to_epsg() == 32644
        assert src.width == pytest.approx(533, abs=2) and src.height == pytest.approx(533, abs=2)
        assert src.read(1).any(), "tiles actually downloaded"
    pin = raster.pinned("35_04_077")
    assert pin["path"] == str(out) and len(pin["sha256"]) == 64 and pin["zoom"] == 20


def test_export_archives_the_previous_raster(tmp_path, monkeypatch):
    monkeypatch.setattr(raster.paths, "PROJECT", tmp_path)
    bounds = (392900.0, 1413300.0, 392940.0, 1413340.0)
    first = raster.export("35_04_077", bounds)
    first.write_bytes(first.read_bytes())                    # touch, keep content
    second = raster.export("35_04_077", bounds)
    archived = list((tmp_path / "_logs").glob("**/satellite_*.tif"))
    assert archived, "the previous raster is archived, never deleted"
    assert second.exists()
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_raster.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.raster'`

- [ ] **Step 3: Write `autogeoref/raster.py`**

```python
"""Export the project's Google Satellite XYZ layer to a georeferenced GeoTIFF and pin it.

The QGIS layer "Google Satellite" is the XYZ source
https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z} in EPSG:3857. rasterio can read it through
GDAL's WMS/TMS driver given the XML description below, so no screenshot and no ECW writer are needed.
"""
import datetime
import hashlib
import json
import shutil

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds

from . import paths

TILE_XML = """<GDAL_WMS>
  <Service name="TMS"><ServerUrl>https://mt1.google.com/vt/lyrs=s&amp;x=${x}&amp;y=${y}&amp;z=${z}</ServerUrl></Service>
  <DataWindow>
    <UpperLeftX>-20037508.34</UpperLeftX><UpperLeftY>20037508.34</UpperLeftY>
    <LowerRightX>20037508.34</LowerRightX><LowerRightY>-20037508.34</LowerRightY>
    <TileLevel>20</TileLevel><TileCountX>1</TileCountX><TileCountY>1</TileCountY><YOrigin>top</YOrigin>
  </DataWindow>
  <Projection>EPSG:3857</Projection>
  <BlockSizeX>256</BlockSizeX><BlockSizeY>256</BlockSizeY><BandsCount>3</BandsCount>
  <MaxConnections>4</MaxConnections>
</GDAL_WMS>
"""


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def export(village, bounds, zoom=20, resolution=0.15, run_id=None):
    """Write the satellite window for `bounds` (EPSG:32644 minx, miny, maxx, maxy) as a GeoTIFF."""
    out_dir = paths.georef_dir(village)
    out_dir.mkdir(parents=True, exist_ok=True)
    desc = out_dir / "google_sat_z%d.xml" % zoom
    desc.write_text(TILE_XML.replace("<TileLevel>20</TileLevel>", "<TileLevel>%d</TileLevel>" % zoom), encoding="utf-8")
    minx, miny, maxx, maxy = bounds
    width = int(round((maxx - minx) / resolution))
    height = int(round((maxy - miny) / resolution))
    dst_tr = from_origin(minx, maxy, resolution, resolution)
    dst = np.zeros((3, height, width), np.uint8)
    with rasterio.open(str(desc)) as src:
        win = src.window(*transform_bounds("EPSG:32644", src.crs, minx, miny, maxx, maxy))
        win = win.round_offsets().round_lengths()
        arr = src.read(window=win)
        src_tr = src.window_transform(win)
        for band in range(3):
            reproject(arr[band], dst[band], src_transform=src_tr, src_crs=src.crs,
                      dst_transform=dst_tr, dst_crs="EPSG:32644", resampling=Resampling.bilinear)
    stamp = (run_id or datetime.datetime.now().strftime("%Y%m%d"))
    out = out_dir / ("satellite_z%d_%s.tif" % (zoom, stamp))
    previous = pinned(village)
    if previous and previous.get("path") and str(out) != previous["path"]:
        old = paths.logs_dir() / ("georef_raster_archive_%s" % stamp)
        old.mkdir(parents=True, exist_ok=True)
        try:
            shutil.move(previous["path"], old / paths.Path(previous["path"]).name)
        except (OSError, FileNotFoundError):
            pass
    profile = {"driver": "GTiff", "height": height, "width": width, "count": 3, "dtype": "uint8",
               "crs": "EPSG:32644", "transform": dst_tr, "compress": "DEFLATE", "tiled": True}
    with rasterio.open(out, "w", **profile) as fh:
        fh.write(dst)
        fh.build_overviews([2, 4, 8, 16], Resampling.average)
    paths.raster_pin(village).write_text(json.dumps(
        {"path": str(out), "sha256": _sha256(out), "zoom": zoom, "resolution_m": resolution,
         "fetched": datetime.datetime.now().isoformat(timespec="seconds"),
         "bounds_32644": [minx, miny, maxx, maxy],
         "source": "QGIS layer 'Google Satellite' (XYZ mt1.google.com/vt/lyrs=s) via GDAL_WMS/TMS",
         "note": "internal review basemap only; not for redistribution"}, indent=1), encoding="utf-8")
    return out


def pinned(village):
    p = paths.raster_pin(village)
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def open_pinned(village):
    pin = pinned(village)
    if not pin:
        raise FileNotFoundError("no raster pinned for %s; run with --raster first" % village)
    if _sha256(pin["path"]) != pin["sha256"]:
        raise ValueError("pinned raster changed on disk: %s" % pin["path"])
    return rasterio.open(pin["path"])
```

Add `from pathlib import Path as _Path` usage: `paths.Path` is not defined, so add `from pathlib import Path` to `raster.py` and use `Path(previous["path"]).name`.

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_raster.py -v`
Expected: PASS, 4 tests. The two export tests need network access to the tile server; if the machine is offline they fail with a rasterio RasterioIOError — mark them `@pytest.mark.network` and skip, but never mock the export away.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/raster.py tests/test_raster.py
git commit -m "feat(autogeoref): satellite GeoTIFF export with sha256 pin and archiving"
```

---

## Task 5: Puvi search window and start poses

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\window.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_window.py`

**Interfaces:**
- Consumes: `paths.puvi_shapefile`, `sheets.load_sheet`, `sheets.dissolve`, `fit.apply_pose`.
- Produces:
  - `window.puvi_polygon(village, survey) -> (Polygon | None, str)` — geometry in EPSG:32644 and the key kind that matched (`"exact"`, `"parent"`, or `""`).
  - `window.start_poses(village, survey, buffer_m=40.0) -> (list[tuple[float, np.ndarray]], Polygon | None)` — one or two hypotheses (the best rotation and, when it scores within 15 %, its 180° twin) plus the search window.
  - `window.puvi_trusted(village, survey, sheet_area) -> (bool, float)` — False when the Puvi polygon's area differs from the sheet's by more than `AREA_GATE` (3x). 22 corridor surveys fail this, 569B by 263x.
  - `window.puvi_distance(village, survey, placed_geom) -> float | None` — reporting only.

- [ ] **Step 1: Write the failing test**

Create `tests/test_window.py`:

```python
import numpy as np
import pytest
from autogeoref import window, anchors, sheets, fit

VILLAGE = "35_04_077"


def test_puvi_lookup_prefers_the_exact_letter_form_key():
    geom, kind = window.puvi_polygon(VILLAGE, "47B")
    assert geom is not None and kind == "exact"
    assert geom.area > 1000


def test_missing_puvi_polygon_is_reported_not_raised():
    geom, kind = window.puvi_polygon(VILLAGE, "99999")
    assert geom is None and kind == ""


def test_start_pose_lands_inside_the_puvi_window():
    poses, win = window.start_poses(VILLAGE, "48A")
    assert 1 <= len(poses) <= 2 and win is not None
    theta, t = poses[0]
    placed = sheets.dissolve([(None, fit.apply_pose(g, theta, t))
                              for _, g in sheets.load_sheet(VILLAGE, "48A")])
    assert win.contains(placed.centroid)


def test_start_pose_heading_is_close_for_a_distinctive_parcel():
    a = anchors.load_anchors(VILLAGE, tool_written=set())["48A"]
    poses, _win = window.start_poses(VILLAGE, "48A")
    diffs = [abs(((theta - a.theta) + 180) % 360 - 180) for theta, _t in poses]
    assert min(diffs) < 10.0, "48A is not symmetric; the Puvi shape should find its heading"


def test_symmetric_railway_strip_keeps_both_twins():
    poses, _win = window.start_poses(VILLAGE, "170")
    assert len(poses) == 2, "a near-symmetric strip must keep the 180 degree twin"
    d = abs(((poses[0][0] - poses[1][0]) + 180) % 360 - 180)
    assert d == pytest.approx(180.0, abs=1.5)


def test_puvi_is_dropped_when_it_describes_a_different_parcel():
    """Thirukatchur 569B: the sheet is 0.98 acre, the Puvi polygon is the 259 acre village tank."""
    ok, ratio = window.puvi_trusted("35_04_074", "569B", sheet_area=3986.0)
    assert ok is False and ratio > 100
    poses, win = window.start_poses("35_04_074", "569B")
    assert poses == [] and win is None, "an untrusted Puvi polygon must not seed the search"


def test_puvi_is_kept_when_the_areas_agree():
    ok, ratio = window.puvi_trusted("35_04_074", "526A", sheet_area=81730.0)
    assert ok is True and 0.9 < ratio < 1.1


def test_puvi_distance_is_reported_for_a_placed_parcel():
    a = anchors.load_anchors(VILLAGE, tool_written=set())["170"]
    placed = anchors.placed_geometry(VILLAGE, a)
    d = window.puvi_distance(VILLAGE, "170", placed)
    assert d is not None and 20 < d < 60, "Puvi 170 sits about 37 m from the truth"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_window.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.window'`

- [ ] **Step 3: Write `autogeoref/window.py`**

```python
"""Puvi-derived search window and starting poses.

Puvi is a location REFERENCE ONLY (user rule): it may say roughly where a survey is and roughly
how it is turned, and its distance is reported. It never accepts, rejects or moves a placement.
Measured offsets on this corridor: 1.6 m to 36.8 m.
"""
import math
import re

import geopandas as gpd
import numpy as np
from shapely import affinity
from shapely.ops import unary_union

from . import paths, sheets

BUFFER_M = 40.0        # railway strips were 37 m off; the window must still contain the truth
COARSE_STEP = 2.0      # degrees
FINE_STEP = 0.5        # degrees
SHIFT_STEP = 2.0       # metres
SHIFT_RANGE = 6.0      # metres
TWIN_RATIO = 1.15      # keep the 180 degree twin when it scores within 15 % of the best
_CACHE = {}


def _puvi(village):
    if village not in _CACHE:
        shp = paths.puvi_shapefile(village)
        if shp is None:
            _CACHE[village] = None
        else:
            gdf = gpd.read_file(shp).to_crs(32644)
            gdf["key"] = gdf["survey_no"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
            _CACHE[village] = {k: unary_union(list(g.geometry)).buffer(0) for k, g in gdf.groupby("key")}
    return _CACHE[village]


def puvi_polygon(village, survey):
    table = _puvi(village)
    if not table:
        return None, ""
    key = str(survey).upper()
    if key in table:
        return table[key], "exact"
    m = re.match(r"^\d+", key)
    if m and m.group() in table:
        return table[m.group()], "parent"
    return None, ""


AREA_GATE = 3.0        # Puvi and the portal describe different parcels beyond this factor


def puvi_trusted(village, survey, sheet_area):
    """(trusted, puvi_area / sheet_area). 22 corridor surveys fail this; 569B by a factor of 263."""
    target, _kind = puvi_polygon(village, survey)
    if target is None or sheet_area <= 0:
        return False, float("nan")
    ratio = target.area / sheet_area
    return (1.0 / AREA_GATE) <= ratio <= AREA_GATE, ratio


def start_poses(village, survey, buffer_m=BUFFER_M):
    target, _kind = puvi_polygon(village, survey)
    if target is None:
        return [], None
    body = sheets.dissolve(sheets.load_sheet(village, survey))
    trusted, _ratio = puvi_trusted(village, survey, body.area)
    if not trusted:
        return [], None                 # different parcel: fall back to the neighbour window
    c = np.array(body.centroid.coords[0])
    tc = np.array(target.centroid.coords[0])

    def score(theta, dx=0.0, dy=0.0):
        g = affinity.rotate(body, theta, origin=(c[0], c[1]))
        g = affinity.translate(g, tc[0] - c[0] + dx, tc[1] - c[1] + dy)
        return g.symmetric_difference(target).area

    coarse = sorted(((score(th), th) for th in np.arange(0.0, 360.0, COARSE_STEP)), key=lambda x: x[0])
    best_s, best_th = coarse[0]
    twin_th = (best_th + 180.0) % 360.0
    twin_s = min(s for s, th in coarse if abs(((th - twin_th) + 180) % 360 - 180) <= COARSE_STEP)
    starts = [best_th] + ([twin_th] if twin_s <= best_s * TWIN_RATIO else [])

    poses = []
    for th0 in starts:
        best = (float("inf"), th0, 0.0, 0.0)
        for th in np.arange(th0 - 3.0, th0 + 3.0001, FINE_STEP):
            for dx in np.arange(-SHIFT_RANGE, SHIFT_RANGE + 0.001, SHIFT_STEP):
                for dy in np.arange(-SHIFT_RANGE, SHIFT_RANGE + 0.001, SHIFT_STEP):
                    s = score(th, dx, dy)
                    if s < best[0]:
                        best = (s, th, dx, dy)
        _s, th, dx, dy = best
        rad = math.radians(th)
        R = np.array([[math.cos(rad), -math.sin(rad)], [math.sin(rad), math.cos(rad)]])
        poses.append((float(th), tc + np.array([dx, dy]) - R @ c))
    return poses, target.buffer(buffer_m)


def puvi_distance(village, survey, placed_geom):
    target, _kind = puvi_polygon(village, survey)
    if target is None or placed_geom is None or placed_geom.is_empty:
        return None
    return float(placed_geom.centroid.distance(target.centroid))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_window.py -v`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/window.py tests/test_window.py
git commit -m "feat(autogeoref): Puvi search window and twin start poses (reference only)"
```

---

## Task 6: Neighbour adjacency from the two-reader transcription

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\neighbours.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_neighbours.py`

**Interfaces:**
- Consumes: `paths.vector_dir`.
- Produces:
  - `neighbours.load(village) -> dict | None`
  - `neighbours.normalise(number) -> str | None`
  - `neighbours.printed(village, survey) -> set[str]`
  - `neighbours.are_neighbours(village, a, b) -> bool | None` — True when either sheet prints the other, False when both are transcribed and neither does, None when unknown.
  - `neighbours.side(village, survey, other) -> str | None`
  - `neighbours.side_vector(letters) -> tuple[float, float]`

- [ ] **Step 1: Write the failing test**

Create `tests/test_neighbours.py`:

```python
import pytest
from autogeoref import neighbours

VILLAGE = "35_04_077"


def test_normalise_strips_prefixes_and_slashes():
    assert neighbours.normalise("103/5A") == "103"
    assert neighbours.normalise("S.No 47B") == "47B"
    assert neighbours.normalise(" 170 ") == "170"
    assert neighbours.normalise("V.No. 74 THIRUKACHUR") is None, "a village number is not a survey"
    assert neighbours.normalise("Railway") is None


def test_printed_numbers_match_the_transcription():
    got = neighbours.printed(VILLAGE, "171")
    assert {"170", "48A", "48B", "47B", "47A", "46A", "46B"} <= got


def test_mutual_naming_makes_a_neighbour_pair():
    assert neighbours.are_neighbours(VILLAGE, "171", "170") is True
    assert neighbours.are_neighbours(VILLAGE, "43A", "46B") is True


def test_two_fully_read_sheets_that_ignore_each_other_are_not_neighbours():
    assert neighbours.are_neighbours(VILLAGE, "40A", "171") is False


def test_unknown_when_a_village_has_no_transcription():
    assert neighbours.are_neighbours("35_04_052", "69", "68") is None


def test_side_and_side_vector():
    assert neighbours.side(VILLAGE, "171", "170") == "S"
    assert neighbours.side_vector("N") == (0.0, 1.0)
    assert neighbours.side_vector("SE") == pytest.approx((0.7071, -0.7071), abs=1e-3)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_neighbours.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.neighbours'`

- [ ] **Step 3: Write `autogeoref/neighbours.py`**

The transcription lives at `FMB_Vector\<village>\neighbour_transcription\nb_override_<village>.json`, shape `{survey: [{"number": "170", "side": "N", "readers": 2, "confidence": "high"}, …]}`. It is produced once per village by two independent readers of the sheet images; the glyph reader of 2026-09-17 misread `170` as `17` and took the adjoining-village label `V.No. 74 Thirukachur` for a survey, so it is not used here.

```python
"""Which survey numbers each FMB sheet prints outside its own outline, and on which side."""
import json
import re

from . import paths

SIDE_VEC = {"N": (0.0, 1.0), "NE": (0.7071, 0.7071), "E": (1.0, 0.0), "SE": (0.7071, -0.7071),
            "S": (0.0, -1.0), "SW": (-0.7071, -0.7071), "W": (-1.0, 0.0), "NW": (-0.7071, 0.7071)}
_CACHE = {}


def load(village):
    if village not in _CACHE:
        p = paths.vector_dir(village) / "neighbour_transcription" / ("nb_override_%s.json" % village)
        _CACHE[village] = json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    return _CACHE[village]


def normalise(number):
    """'103/5A' -> '103'; 'V.No. 74 ...' -> None (an adjoining village, not a survey)."""
    s = str(number).upper().strip()
    if "V.NO" in s or "VILLAGE" in s:
        return None
    s = re.sub(r"^S\.?\s*NO\.?\s*", "", s).replace(" ", "").split("/")[0]
    return s if re.fullmatch(r"\d+[A-Z]?", s) else None


def printed(village, survey):
    table = load(village)
    if not table or survey not in table:
        return set()
    out = set()
    for entry in table[survey]:
        n = normalise(entry.get("number", ""))
        if not n or n == str(survey).upper():
            continue
        if entry.get("readers", 1) >= 2 or entry.get("confidence") == "high":
            out.add(n)
    return out


def are_neighbours(village, a, b):
    table = load(village)
    if not table or a not in table or b not in table:
        return None
    return str(b).upper() in printed(village, a) or str(a).upper() in printed(village, b)


def side(village, survey, other):
    table = load(village)
    if not table or survey not in table:
        return None
    for entry in table[survey]:
        if normalise(entry.get("number", "")) == str(other).upper():
            s = str(entry.get("side", "")).upper()
            return s if s in SIDE_VEC else None
    return None


def side_vector(letters):
    return SIDE_VEC.get(str(letters).upper(), (0.0, 0.0))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_neighbours.py -v`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/neighbours.py tests/test_neighbours.py
git commit -m "feat(autogeoref): neighbour adjacency and sides from the two-reader transcription"
```

---

## Task 7: Shared-boundary chains and their observations

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\match.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_match.py`

**Interfaces:**
- Consumes: `sheets.edge_lengths`, `sheets.turn_angles`.
- Produces:
  - `match.Chain` dataclass: `length: float`, `pairs: list[tuple[tuple, tuple]]`, `n_pairs: int`, `corners: int`.
  - `match.PairObs` dataclass: `a: str, b: str, pa: tuple, pb: tuple, sigma: float` — two sheet vertices that are the same ground point.
  - `match.LineObs` dataclass: `survey: str, p_local: tuple, line: tuple, sigma: float, distance_now: float` — a sheet sample that must lie on a fixed anchor's boundary line.
  - `match.common_chains(vq, vp, tol_abs=0.30, tol_rel=0.01, min_len=8.0, pass_deg=5.0) -> list[Chain]` — longest first.
  - `match.chain_observations(chain, a, b, sigma) -> list[PairObs]`
  - `match.line_observations(survey, samples, anchor_boundary, sigma, max_dist=2.0) -> list[LineObs]`

- [ ] **Step 1: Write the failing test**

Create `tests/test_match.py`:

```python
import numpy as np
import pytest
from shapely.geometry import Polygon, LineString
from autogeoref import match, sheets, fit

A = Polygon([(0, 0), (40, 0), (40, 25), (0, 25)])       # shares the x = 40 edge with B
B = Polygon([(40, 0), (75, 0), (75, 25), (40, 25)])


def test_common_chains_finds_the_shared_edge_of_two_rectangles():
    chains = match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, B)]))
    assert chains, "the 25 m shared edge must be found"
    assert chains[0].length == pytest.approx(25.0, abs=0.5)
    assert chains[0].n_pairs >= 2


def test_common_chains_returns_every_run_not_only_the_longest():
    chains = match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, B)]))
    assert len(chains) >= 2, "returning only the longest run hid true matches in the prototype"
    assert chains[0].length >= chains[-1].length


def test_a_chain_pins_the_pose_of_a_moved_neighbour():
    moved = fit.apply_pose(B, 20.0, np.array([300.0, 400.0]))
    chains = match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, moved)]))
    best = chains[0]
    P = np.array([q for _p, q in best.pairs])
    Q = np.array([p for p, _q in best.pairs])
    theta, t, rms, _mx = fit.rigid_fit(P, Q)
    assert rms < 0.05
    assert abs(((theta + 20.0) + 180) % 360 - 180) < 0.5


def test_no_chain_between_parcels_that_do_not_touch():
    far = Polygon([(500, 500), (540, 500), (540, 525), (500, 525)])
    assert match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, far)])) == []


def test_line_observations_measure_perpendicular_distance_to_an_anchor_edge():
    obs = match.line_observations("42B", np.array([[40.5, 5.0], [40.5, 15.0]]),
                                  LineString([(40, 0), (40, 25)]), sigma=0.30)
    assert len(obs) == 2
    assert all(abs(o.distance_now - 0.5) < 1e-9 for o in obs)
    assert all(o.sigma == 0.30 for o in obs)


def test_line_observations_ignore_samples_further_than_max_dist():
    obs = match.line_observations("42B", np.array([[40.5, 5.0], [55.0, 5.0]]),
                                  LineString([(40, 0), (40, 25)]), sigma=0.30, max_dist=2.0)
    assert len(obs) == 1
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_match.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.match'`

- [ ] **Step 3: Write `autogeoref/match.py`**

The chainage walk is ported from `autogeoref/scratch_2026-09-17/autogeoref_demo.py::common_chain`, which is already correct: it walks both rings by running distance from every start-corner pair, matches corners at equal running distance, passes through a corner that turns less than `pass_deg` (the two sheets draw a near-straight corner differently), and ends the walk where the boundaries really diverge.

```python
"""Shared-boundary matching between two sheet outlines, and the observations it yields."""
from dataclasses import dataclass

import numpy as np
from shapely.geometry import LineString, Point

from . import sheets


@dataclass
class Chain:
    length: float
    pairs: list
    n_pairs: int
    corners: int


@dataclass
class PairObs:
    a: str
    b: str
    pa: tuple
    pb: tuple
    sigma: float


@dataclass
class LineObs:
    survey: str
    p_local: tuple
    line: tuple
    sigma: float
    distance_now: float


def common_chains(vq, vp, tol_abs=0.30, tol_rel=0.01, min_len=8.0, pass_deg=5.0):
    """Every run of consecutive equal-length edges shared by rings vq and vp, longest first."""
    def tol(x):
        return tol_abs + tol_rel * max(x, 1.0)

    eq = sheets.edge_lengths(vq)
    tq = sheets.turn_angles(vq)
    nq = len(eq)
    seen, out = set(), []
    for reverse in (True, False):
        vp2 = vp[::-1] if reverse else vp
        ep = sheets.edge_lengths(vp2)
        tp = sheets.turn_angles(vp2)
        npn = len(ep)
        for i in range(nq):
            for j in range(npn):
                pairs = [(vq[i], vp2[j])]
                cq = cp = matched = 0.0
                iq = ip = corners = 0
                while iq < nq and ip < npn:
                    q_next = cq + eq[(i + iq) % nq]
                    p_next = cp + ep[(j + ip) % npn]
                    if abs(q_next - p_next) <= tol(max(q_next, p_next)):
                        cq, cp, iq, ip = q_next, p_next, iq + 1, ip + 1
                        pairs.append((vq[(i + iq) % nq], vp2[(j + ip) % npn]))
                        matched = min(cq, cp)
                        if tq[(i + iq) % nq] >= pass_deg:
                            corners += 1
                        if iq >= nq or ip >= npn:
                            break
                    elif q_next < p_next:
                        iq += 1
                        cq = q_next
                        if tq[(i + iq) % nq] > pass_deg:
                            break
                    else:
                        ip += 1
                        cp = p_next
                        if tp[(j + ip) % npn] > pass_deg:
                            break
                if len(pairs) >= 2 and matched >= min_len:
                    key = tuple(sorted((round(a[0], 2), round(a[1], 2), round(b[0], 2), round(b[1], 2))
                                       for a, b in pairs))
                    if key not in seen:
                        seen.add(key)
                        out.append(Chain(matched, pairs, len(pairs), corners))
    out.sort(key=lambda c: -c.length)
    return out


def chain_observations(chain, a, b, sigma):
    return [PairObs(a, b, tuple(pa), tuple(pb), sigma) for pa, pb in chain.pairs]


def line_observations(survey, samples, anchor_boundary, sigma, max_dist=2.0):
    """Point-to-line observations against a fixed anchor's boundary.

    Point-to-line rather than point-to-point because the anchor's own geometry is affine: its
    vertices sit 1.4-6.4 m from any rigid placement of its sheet, but the line it draws is still
    the boundary. Only samples already within max_dist are used, so a wrong pose cannot recruit
    the far side of the parcel.
    """
    geoms = anchor_boundary.geoms if anchor_boundary.geom_type.startswith("Multi") else [anchor_boundary]
    out = []
    for p in np.asarray(samples, float):
        pt = Point(p)
        best = None
        for g in geoms:
            coords = list(g.coords)
            for k in range(len(coords) - 1):
                d = LineString([coords[k], coords[k + 1]]).distance(pt)
                if best is None or d < best[0]:
                    best = (d, (coords[k], coords[k + 1]))
        if best is None or best[0] > max_dist:
            continue
        out.append(LineObs(survey, (float(p[0]), float(p[1])), best[1], sigma, float(best[0])))
    return out
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_match.py -v`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/match.py tests/test_match.py
git commit -m "feat(autogeoref): shared-boundary chains with point-to-line anchor observations"
```

---

## Task 8: Satellite line segments and their index

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\edges.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_edges.py`

**Interfaces:**
- Consumes: nothing (takes a raster path).
- Produces:
  - `edges.Segment` namedtuple: `p0, p1, bearing (radians 0..pi), length`.
  - `edges.detect(raster_path, bounds=None, min_len_m=3.0) -> list[Segment]`
  - `edges.directions(segments, tol_deg=10.0) -> list[tuple[float, float]]` — `(bearing, total_length)`, longest first.
  - `edges.SegmentIndex(segments, step=0.5)` with `.query(points, max_dist=2.5) -> (dist, bearing)`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_edges.py`:

```python
import math
import numpy as np
import pytest
from autogeoref import edges

SPIKE = r"D:\code\FMB_to_GeoJSON\autogeoref\scratch_2026-09-17\spike_gcp\kizhi_sat_z20.tif"


def test_detect_finds_thousands_of_segments_on_the_spike_raster():
    segs = edges.detect(SPIKE, min_len_m=4.0)
    assert 1500 < len(segs) < 6000, "the 2026-09-18 spike found 2643 segments >= 4 m"
    assert all(s.length >= 4.0 for s in segs)
    assert all(0 <= s.bearing < math.pi for s in segs)


def test_segments_carry_ground_coordinates():
    segs = edges.detect(SPIKE, min_len_m=4.0)
    xs = [s.p0[0] for s in segs]
    ys = [s.p0[1] for s in segs]
    assert 392000 < min(xs) and max(xs) < 394000
    assert 1412000 < min(ys) and max(ys) < 1415000


def test_index_returns_distance_and_bearing_of_the_nearest_segment():
    idx = edges.SegmentIndex([edges.Segment((0.0, 0.0), (100.0, 0.0), 0.0, 100.0)], step=0.5)
    d, b = idx.query(np.array([[50.0, 2.0], [50.0, 40.0]]))
    assert d[0] == pytest.approx(2.0, abs=0.3)
    assert b[0] == pytest.approx(0.0, abs=1e-6)
    assert not np.isfinite(d[1]), "beyond max_dist there is no nearest segment"


def test_directions_groups_parallel_segments():
    segs = [edges.Segment((0, 0), (100, 0), 0.0, 100.0),
            edges.Segment((0, 10), (60, 10), 0.0, 60.0),
            edges.Segment((0, 0), (0, 40), math.pi / 2, 40.0)]
    got = edges.directions(segs, tol_deg=10.0)
    assert len(got) == 2
    assert got[0][1] == pytest.approx(160.0)


def test_the_corridor_parallel_bias_is_real():
    """2026-09-18 spike: 92 % of matched outline length lies within 15 deg of the corridor bearing."""
    segs = edges.detect(SPIKE, min_len_m=4.0)
    corridor = math.radians(25.8)
    par = sum(s.length for s in segs
              if min(abs(s.bearing - corridor), math.pi - abs(s.bearing - corridor)) <= math.radians(15))
    assert par / sum(s.length for s in segs) > 0.2
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_edges.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.edges'`

- [ ] **Step 3: Write `autogeoref/edges.py`**

```python
"""Straight edges visible in the satellite raster: walls, roads, bunds, the railway fence.

Measured on the 2026-09-18 spike: 53 % of the team's outline length lies within 1.5 m of a
same-direction segment, but 92 % of that is parallel to the corridor. Imagery therefore fixes
rotation and cross-track offset; along-track position has to come from neighbours.
"""
import math
from collections import namedtuple

import cv2
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from scipy.spatial import cKDTree

Segment = namedtuple("Segment", "p0 p1 bearing length")


def detect(raster_path, bounds=None, min_len_m=3.0):
    """Line segments in EPSG:32644 with bearings in radians 0..pi."""
    with rasterio.open(str(raster_path)) as src:
        if bounds is not None:
            win = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
            arr = src.read(window=win)
            tr = src.window_transform(win)
        else:
            arr = src.read()
            tr = src.transform
        res = src.res[0]
    gray = cv2.cvtColor(np.ascontiguousarray(np.moveaxis(arr, 0, -1)), cv2.COLOR_RGB2GRAY)
    gray = cv2.bilateralFilter(gray, 7, 40, 7)
    lines = cv2.createLineSegmentDetector(0).detect(gray)[0]
    out = []
    for line in (lines if lines is not None else []):
        x0, y0, x1, y1 = [float(v) for v in np.ravel(line)[:4]]
        if math.hypot(x1 - x0, y1 - y0) * res < min_len_m:
            continue
        a = tr * (x0, y0)
        b = tr * (x1, y1)
        out.append(Segment((a[0], a[1]), (b[0], b[1]),
                           math.atan2(b[1] - a[1], b[0] - a[0]) % math.pi, math.dist(a, b)))
    return out


def directions(segments, tol_deg=10.0):
    """[(bearing, total_length)] after grouping segments whose bearings agree within tol_deg."""
    tol = math.radians(tol_deg)
    groups = []
    for s in sorted(segments, key=lambda s: -s.length):
        for i, (b, L) in enumerate(groups):
            if min(abs(s.bearing - b), math.pi - abs(s.bearing - b)) <= tol:
                groups[i] = (b, L + s.length)
                break
        else:
            groups.append((s.bearing, s.length))
    return sorted(groups, key=lambda g: -g[1])


class SegmentIndex:
    """KD-tree over points sampled along every segment, so a pose is scored in one vectorised query."""

    def __init__(self, segments, step=0.5):
        pts, brg = [], []
        for s in segments:
            a = np.array(s.p0, float)
            b = np.array(s.p1, float)
            n = max(2, int(s.length / step))
            for k in range(n + 1):
                pts.append(a + (b - a) * k / n)
                brg.append(s.bearing)
        self.points = np.array(pts) if pts else np.zeros((0, 2))
        self.bearings = np.array(brg) if brg else np.zeros(0)
        self.tree = cKDTree(self.points) if len(self.points) else None
        self.segments = segments

    def query(self, points, max_dist=2.5):
        """(distance, bearing) of the nearest sampled segment point; distance is inf beyond max_dist."""
        points = np.asarray(points, float)
        if self.tree is None:
            return np.full(len(points), np.inf), np.zeros(len(points))
        d, j = self.tree.query(points, distance_upper_bound=max_dist)
        j = np.where(np.isfinite(d), j, 0)
        return d, self.bearings[j]
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_edges.py -v`
Expected: PASS, 5 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/edges.py tests/test_edges.py
git commit -m "feat(autogeoref): satellite line-segment detection with a sampled KD-tree index"
```

---

## Task 9: Image alignment with margin and observability tests

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\align.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_align.py`

**Interfaces:**
- Consumes: `edges.SegmentIndex`, `edges.directions`, `sheets.sample_outline`, `fit.transform_points`.
- Produces:
  - `align.share(samples, bearings, index, theta_deg, t) -> float` — fraction of 1 m outline samples lying within `MATCH_DIST` (1.5 m) of a segment whose direction agrees within `MATCH_DEG` (12°). This is the binary metric the 2026-09-18 spike measured, so the calibrated thresholds of Task 16 apply to it.
  - `align.AlignResult` dataclass: `theta, t, share, margin, alternative, ambiguous, matched_by_direction, observable, note`.
  - `align.search(samples, bearings, index, hypotheses, shift_m=10.0, rot_deg=10.0) -> AlignResult`
  - `align.observability(samples, bearings, index, theta, t) -> dict` — `{"primary_m": float, "cross_m": float, "cross_angle_deg": float, "ok": bool}`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_align.py`:

```python
import math
import numpy as np
import pytest
from shapely.geometry import Polygon
from autogeoref import align, edges, sheets, fit


def _square_index(x0=0.0, y0=0.0, side=40.0, step=0.5):
    """Four walls of a square, as if seen in imagery."""
    segs = [edges.Segment((x0, y0), (x0 + side, y0), 0.0, side),
            edges.Segment((x0 + side, y0), (x0 + side, y0 + side), math.pi / 2, side),
            edges.Segment((x0 + side, y0 + side), (x0, y0 + side), 0.0, side),
            edges.Segment((x0, y0 + side), (x0, y0), math.pi / 2, side)]
    return edges.SegmentIndex(segs, step=step)


SHEET = [(None, Polygon([(0, 0), (40, 0), (40, 40), (0, 40)]))]


def test_share_is_one_at_the_true_pose_and_low_when_far_away():
    v = sheets.outline(SHEET)
    pts, brg = sheets.sample_outline(v, step=1.0)
    idx = _square_index()
    assert align.share(pts, brg, idx, 0.0, np.array([0.0, 0.0])) > 0.95
    assert align.share(pts, brg, idx, 0.0, np.array([25.0, 0.0])) < 0.60


def test_search_recovers_a_shifted_pose_and_reports_a_margin():
    v = sheets.outline(SHEET)
    pts, brg = sheets.sample_outline(v, step=1.0)
    idx = _square_index()
    res = align.search(pts, brg, idx, [(0.0, np.array([6.0, -5.0]))], shift_m=10.0, rot_deg=10.0)
    assert np.allclose(res.t, [0.0, 0.0], atol=0.5)
    assert res.share > 0.95
    assert res.ambiguous is False


def test_parallel_lines_are_reported_ambiguous_not_placed():
    """Rails and road markings repeat every few metres; the score then has many equal peaks."""
    segs = [edges.Segment((0.0, y), (200.0, y), 0.0, 200.0) for y in (0.0, 3.0, 6.0, 9.0, 12.0)]
    idx = edges.SegmentIndex(segs, step=0.5)
    strip = [(None, Polygon([(0, 0), (150, 0), (150, 3), (0, 3)]))]
    pts, brg = sheets.sample_outline(sheets.outline(strip), step=1.0)
    res = align.search(pts, brg, idx, [(0.0, np.array([0.0, 0.5]))], shift_m=8.0, rot_deg=4.0)
    assert res.ambiguous is True
    assert res.alternative is not None


def test_observability_needs_two_directions():
    v = sheets.outline(SHEET)
    pts, brg = sheets.sample_outline(v, step=1.0)
    full = align.observability(pts, brg, _square_index(), 0.0, np.array([0.0, 0.0]))
    assert full["ok"] is True and full["cross_m"] >= 20.0

    one_way = edges.SegmentIndex([edges.Segment((0.0, 0.0), (40.0, 0.0), 0.0, 40.0),
                                  edges.Segment((0.0, 40.0), (40.0, 40.0), 0.0, 40.0)], step=0.5)
    partial = align.observability(pts, brg, one_way, 0.0, np.array([0.0, 0.0]))
    assert partial["ok"] is False, "one direction leaves the along-track position free"
    assert partial["cross_m"] < 20.0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_align.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.align'`

- [ ] **Step 3: Write `autogeoref/align.py`**

```python
"""Fit a sheet outline to the straight edges visible in the satellite raster.

Two guards make this safe on a railway corridor, where the 2026-09-18 trial showed a free image
search sliding up to 28 m or flipping 180 degrees:

* margin  - the best pose must beat the best pose more than EXCL_M / EXCL_DEG away by MARGIN_MIN,
            and every local maximum within 75 % of the best must lie within 1 m of it;
* observability - the matched edges must span two directions (>= 30 degrees apart), otherwise the
            along-track position is not observable from imagery and a neighbour must supply it.
"""
import math
from dataclasses import dataclass

import numpy as np

from . import fit

MATCH_DIST = 1.5          # metres: the spike's own tolerance, so calibrated thresholds apply
MATCH_DEG = 12.0          # degrees
EXCL_M = 3.0              # margin exclusion zone
EXCL_DEG = 3.0
MARGIN_MIN = 0.15         # the best pose must beat the best outside pose by this much share
PEAK_RATIO = 0.75         # local maxima above this fraction of the best count as competitors
PRIMARY_MIN = 40.0        # metres of matched outline in the main direction
CROSS_MIN = 20.0          # metres matched in a direction >= CROSS_ANGLE from it
CROSS_ANGLE = 30.0        # degrees


@dataclass
class AlignResult:
    theta: float
    t: np.ndarray
    share: float
    margin: float
    alternative: tuple
    ambiguous: bool
    matched_by_direction: list
    observable: bool
    note: str


def _matched_mask(samples, bearings, index, theta_deg, t):
    pts = fit.transform_points(samples, theta_deg, t)
    d, b = index.query(pts, max_dist=MATCH_DIST + 1.0)
    turned = (bearings + math.radians(theta_deg)) % math.pi
    diff = np.abs((b - turned + math.pi / 2) % math.pi - math.pi / 2)
    return np.isfinite(d) & (d <= MATCH_DIST) & (diff <= math.radians(MATCH_DEG))


def share(samples, bearings, index, theta_deg, t):
    """Fraction of outline samples lying on a same-direction image edge."""
    if len(samples) == 0:
        return 0.0
    return float(_matched_mask(samples, bearings, index, theta_deg, t).mean())


def observability(samples, bearings, index, theta_deg, t):
    """How much matched outline length lies in the main direction and across it."""
    mask = _matched_mask(samples, bearings, index, theta_deg, t)
    turned = ((bearings + math.radians(theta_deg)) % math.pi)[mask]
    if turned.size == 0:
        return {"primary_m": 0.0, "cross_m": 0.0, "cross_angle_deg": 0.0, "ok": False}
    groups = []                                     # one metre per matched sample by construction
    for b in turned:
        for i, (gb, n) in enumerate(groups):
            if min(abs(b - gb), math.pi - abs(b - gb)) <= math.radians(10.0):
                groups[i] = (gb, n + 1)
                break
        else:
            groups.append((b, 1))
    groups.sort(key=lambda g: -g[1])
    primary_b, primary_n = groups[0]
    cross = [(gb, n) for gb, n in groups[1:]
             if min(abs(gb - primary_b), math.pi - abs(gb - primary_b)) >= math.radians(CROSS_ANGLE)]
    cross_b, cross_n = (cross[0] if cross else (primary_b, 0))
    ang = math.degrees(min(abs(cross_b - primary_b), math.pi - abs(cross_b - primary_b)))
    return {"primary_m": float(primary_n), "cross_m": float(cross_n), "cross_angle_deg": round(ang, 1),
            "ok": primary_n >= PRIMARY_MIN and cross_n >= CROSS_MIN}


def search(samples, bearings, index, hypotheses, shift_m=10.0, rot_deg=10.0,
           coarse_m=1.0, coarse_deg=1.0, fine_m=0.25, fine_deg=0.25):
    """Best pose around the given hypotheses, with the margin and observability verdicts."""
    grid = []
    for theta0, t0 in hypotheses:
        t0 = np.asarray(t0, float)
        for dth in np.arange(-rot_deg, rot_deg + 1e-9, coarse_deg):
            for dx in np.arange(-shift_m, shift_m + 1e-9, coarse_m):
                for dy in np.arange(-shift_m, shift_m + 1e-9, coarse_m):
                    th = theta0 + dth
                    t = t0 + np.array([dx, dy])
                    grid.append((share(samples, bearings, index, th, t), th, t))
    if not grid:
        return AlignResult(0.0, np.zeros(2), 0.0, 0.0, None, True, [], False, "no hypothesis")
    grid.sort(key=lambda g: -g[0])
    best_s, best_th, best_t = grid[0]
    # refine
    for dth in np.arange(-coarse_deg, coarse_deg + 1e-9, fine_deg):
        for dx in np.arange(-coarse_m, coarse_m + 1e-9, fine_m):
            for dy in np.arange(-coarse_m, coarse_m + 1e-9, fine_m):
                th = best_th + dth
                t = best_t + np.array([dx, dy])
                s = share(samples, bearings, index, th, t)
                if s > best_s:
                    best_s, best_th, best_t = s, th, t
    # margin: the strongest competitor outside the exclusion zone
    competitor = None
    for s, th, t in grid:
        if float(np.linalg.norm(t - best_t)) <= EXCL_M and abs(((th - best_th) + 180) % 360 - 180) <= EXCL_DEG:
            continue
        competitor = (s, th, t)
        break
    margin = best_s - (competitor[0] if competitor else 0.0)
    near_peaks = [(s, th, t) for s, th, t in grid
                  if s >= PEAK_RATIO * best_s and float(np.linalg.norm(t - best_t)) > 1.0]
    ambiguous = margin < MARGIN_MIN or bool(near_peaks)
    obs = observability(samples, bearings, index, best_th, best_t)
    note = ""
    if ambiguous and competitor is not None:
        dx, dy = competitor[2] - best_t
        note = "alternative pose %.1f m away at bearing %.0f deg" % (
            math.hypot(dx, dy), math.degrees(math.atan2(dx, dy)) % 360)
    elif not obs["ok"]:
        note = "matched edges in one direction only (%.0f m); along-track from neighbours" % obs["primary_m"]
    return AlignResult(float(best_th), np.asarray(best_t, float), float(best_s), float(margin),
                       (competitor[1], competitor[2]) if competitor else None, bool(ambiguous),
                       [obs["primary_m"], obs["cross_m"], obs["cross_angle_deg"]], bool(obs["ok"]), note)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_align.py -v`
Expected: PASS, 4 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/align.py tests/test_align.py
git commit -m "feat(autogeoref): image alignment with margin and observability guards"
```

---

## Task 10: Block adjustment with anchors held fixed

**Files:**
- Modify: `D:\code\FMB_to_GeoJSON\autogeoref\fit.py` (append the adjustment; leave Task 2's functions untouched)
- Create: `D:\code\FMB_to_GeoJSON\tests\test_block_adjust.py`

**Interfaces:**
- Consumes: `match.PairObs`, `match.LineObs`, `fit.transform_points`.
- Produces:
  - `fit.GcpObs` dataclass: `survey: str, p_local: tuple, map_xy: tuple, sigma: float`.
  - `fit.PosePrior` dataclass: `survey: str, theta: float, t: tuple, sigma_pos: float, sigma_head: float`.
  - `fit.block_adjust(free: dict[str, tuple[float, np.ndarray]], fixed: dict[str, tuple[float, np.ndarray]], pair_obs, line_obs, gcp_obs, priors) -> dict[str, dict]` — returns per survey `{"theta", "t", "shift_m", "dtheta_deg", "sigma_pos_m", "sigma_head_deg"}`. Solved in a frame centred on the block (raw UTM unknowns stall `least_squares`; this was the 2026-09-17 bug). Huber loss, `x_scale="jac"`, tolerances 1e-10.
  - `fit.residuals_by_pair(result, pair_obs) -> dict[tuple[str, str], dict]` — `{"n", "rms", "max"}`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_block_adjust.py`:

```python
import math
import numpy as np
import pytest
from autogeoref import fit
from autogeoref.match import LineObs, PairObs


def test_a_free_parcel_snaps_onto_a_fixed_anchor_boundary():
    # the free parcel's local edge x = 0 must land on the anchor's ground line x = 100
    line = ((100.0, 0.0), (100.0, 50.0))
    obs = [LineObs("B", (0.0, y), line, 0.30, 0.0) for y in (5.0, 15.0, 25.0, 35.0)]
    out = fit.block_adjust(free={"B": (0.0, np.array([98.0, 0.0]))}, fixed={},
                           pair_obs=[], line_obs=obs, gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, (98.0, 0.0), 10.0, 10.0)])
    assert out["B"]["t"][0] == pytest.approx(100.0, abs=0.05)
    assert abs(out["B"]["dtheta_deg"]) < 0.5


def test_a_fixed_anchor_never_moves():
    obs = [PairObs("A", "B", (10.0, 0.0), (0.0, 0.0), 0.30),
           PairObs("A", "B", (10.0, 20.0), (0.0, 20.0), 0.30)]
    out = fit.block_adjust(free={"B": (0.0, np.array([12.0, 1.0]))},
                           fixed={"A": (0.0, np.array([0.0, 0.0]))},
                           pair_obs=obs, line_obs=[], gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, (12.0, 1.0), 10.0, 10.0)])
    assert "A" not in out or np.allclose(out["A"]["t"], [0.0, 0.0], atol=1e-9)
    assert out["B"]["t"] == pytest.approx(np.array([10.0, 0.0]), abs=0.05)


def test_two_free_parcels_meet_each_other():
    obs = [PairObs("A", "B", (10.0, 0.0), (0.0, 0.0), 0.30),
           PairObs("A", "B", (10.0, 20.0), (0.0, 20.0), 0.30)]
    out = fit.block_adjust(free={"A": (0.0, np.array([0.0, 0.0])), "B": (0.0, np.array([11.0, 0.6]))},
                           fixed={}, pair_obs=obs, line_obs=[], gcp_obs=[],
                           priors=[fit.PosePrior("A", 0.0, (0.0, 0.0), 2.5, 3.0),
                                   fit.PosePrior("B", 0.0, (11.0, 0.6), 2.5, 3.0)])
    a = fit.transform_points([(10.0, 0.0)], out["A"]["theta"], out["A"]["t"])[0]
    b = fit.transform_points([(0.0, 0.0)], out["B"]["theta"], out["B"]["t"])[0]
    assert float(np.linalg.norm(a - b)) < 0.15


def test_gcp_observations_pull_a_parcel_with_no_neighbour():
    gcps = [fit.GcpObs("C", (0.0, 0.0), (500.0, 500.0), 1.0),
            fit.GcpObs("C", (30.0, 0.0), (530.0, 500.0), 1.0)]
    out = fit.block_adjust(free={"C": (0.0, np.array([495.0, 498.0]))}, fixed={},
                           pair_obs=[], line_obs=[], gcp_obs=gcps,
                           priors=[fit.PosePrior("C", 0.0, (495.0, 498.0), 10.0, 10.0)])
    assert out["C"]["t"] == pytest.approx(np.array([500.0, 500.0]), abs=0.2)


def test_the_solver_works_in_raw_utm_coordinates():
    """2026-09-17 bug: with unknowns near 1.4e6 the solver stopped at x0 and moved nothing."""
    line = ((392_500.0, 1_413_000.0), (392_500.0, 1_413_050.0))
    obs = [LineObs("B", (0.0, y), line, 0.30, 0.0) for y in (5.0, 15.0, 25.0, 35.0)]
    start = np.array([392_498.0, 1_413_000.0])
    out = fit.block_adjust(free={"B": (0.0, start)}, fixed={}, pair_obs=[], line_obs=obs, gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, tuple(start), 10.0, 10.0)])
    assert out["B"]["t"][0] == pytest.approx(392_500.0, abs=0.05)
    assert out["B"]["shift_m"] > 1.0


def test_residuals_by_pair_reports_rms_and_max():
    obs = [PairObs("A", "B", (10.0, 0.0), (0.0, 0.0), 0.30),
           PairObs("A", "B", (10.0, 20.0), (0.0, 20.0), 0.30)]
    out = fit.block_adjust(free={"B": (0.0, np.array([12.0, 1.0]))},
                           fixed={"A": (0.0, np.array([0.0, 0.0]))},
                           pair_obs=obs, line_obs=[], gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, (12.0, 1.0), 10.0, 10.0)])
    r = fit.residuals_by_pair(out | {"A": {"theta": 0.0, "t": np.array([0.0, 0.0])}}, obs)
    assert r[("A", "B")]["n"] == 2 and r[("A", "B")]["rms"] < 0.15
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_block_adjust.py -v`
Expected: FAIL with `AttributeError: module 'autogeoref.fit' has no attribute 'block_adjust'`

- [ ] **Step 3: Append the adjustment to `autogeoref/fit.py`**

```python
# --------------------------------------------------------------------------- block adjustment
from dataclasses import dataclass as _dataclass

from scipy.optimize import least_squares


@_dataclass
class GcpObs:
    survey: str
    p_local: tuple
    map_xy: tuple
    sigma: float


@_dataclass
class PosePrior:
    survey: str
    theta: float
    t: tuple
    sigma_pos: float
    sigma_head: float


def block_adjust(free, fixed, pair_obs, line_obs, gcp_obs, priors, f_scale=1.5):
    """Least-squares rotation+shift per free survey; `fixed` surveys never move.

    Solved in a frame centred on the block: with raw UTM unknowns (easting ~3.9e5, northing
    ~1.4e6) scipy's least_squares terminates on xtol at the starting point and nothing moves.
    """
    names = sorted(free)
    if not names:
        return {}
    idx = {s: i for i, s in enumerate(names)}
    origin = np.mean([np.asarray(t, float) for _th, t in free.values()], axis=0)

    x0 = np.zeros(3 * len(names))
    for s in names:
        th, t = free[s]
        x0[3 * idx[s]] = math.radians(th)
        x0[3 * idx[s] + 1:3 * idx[s] + 3] = np.asarray(t, float) - origin

    fixed_pose = {s: (math.radians(th), np.asarray(t, float) - origin) for s, (th, t) in fixed.items()}

    def place(survey, x, p):
        p = np.asarray(p, float)
        if survey in fixed_pose:
            th, t = fixed_pose[survey]
        else:
            i = idx[survey]
            th, t = x[3 * i], x[3 * i + 1:3 * i + 3]
        c, s = math.cos(th), math.sin(th)
        return np.array([c * p[0] - s * p[1] + t[0], s * p[0] + c * p[1] + t[1]])

    def residuals(x):
        r = []
        for o in pair_obs:
            r += list((place(o.a, x, o.pa) - place(o.b, x, o.pb)) / o.sigma)
        for o in line_obs:
            p = place(o.survey, x, o.p_local)
            a = np.asarray(o.line[0], float) - origin
            b = np.asarray(o.line[1], float) - origin
            d = b - a
            n = float(np.hypot(d[0], d[1]))
            perp = ((p[0] - a[0]) * d[1] - (p[1] - a[1]) * d[0]) / n if n else 0.0
            r.append(perp / o.sigma)
        for o in gcp_obs:
            r += list((place(o.survey, x, o.p_local) - (np.asarray(o.map_xy, float) - origin)) / o.sigma)
        for pr in priors:
            if pr.survey not in idx:
                continue
            i = idx[pr.survey]
            r += list((x[3 * i + 1:3 * i + 3] - (np.asarray(pr.t, float) - origin)) / pr.sigma_pos)
            d = (x[3 * i] - math.radians(pr.theta) + math.pi) % (2 * math.pi) - math.pi
            r.append(d / math.radians(pr.sigma_head))
        return np.array(r) if r else np.zeros(1)

    sol = least_squares(residuals, x0, loss="huber", f_scale=f_scale, x_scale="jac",
                        xtol=1e-10, ftol=1e-10, gtol=1e-10)
    # parameter sigmas from the Jacobian (used for the colour rule)
    try:
        J = sol.jac
        cov = np.linalg.pinv(J.T @ J)
        sig = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        sig = np.full(3 * len(names), np.nan)

    out = {}
    for s in names:
        i = idx[s]
        th = math.degrees(sol.x[3 * i])
        t = sol.x[3 * i + 1:3 * i + 3] + origin
        th0, t0 = free[s]
        out[s] = {"theta": float(th), "t": t,
                  "dtheta_deg": float(((th - th0) + 180) % 360 - 180),
                  "shift_m": float(np.linalg.norm(t - np.asarray(t0, float))),
                  "sigma_head_deg": float(math.degrees(sig[3 * i])),
                  "sigma_pos_m": float(np.hypot(sig[3 * i + 1], sig[3 * i + 2]))}
    return out


def residuals_by_pair(poses, pair_obs):
    """{(a, b): {'n', 'rms', 'max'}} after an adjustment; `poses` maps survey -> {'theta', 't'}."""
    acc = {}
    for o in pair_obs:
        if o.a not in poses or o.b not in poses:
            continue
        pa = transform_points([o.pa], poses[o.a]["theta"], poses[o.a]["t"])[0]
        pb = transform_points([o.pb], poses[o.b]["theta"], poses[o.b]["t"])[0]
        acc.setdefault((o.a, o.b), []).append(float(np.linalg.norm(pa - pb)))
    return {k: {"n": len(v), "rms": float(np.sqrt(np.mean(np.square(v)))), "max": float(max(v))}
            for k, v in acc.items()}
```

Add `import math` and `import numpy as np` at the top of `fit.py` if Task 2 did not already (it did).

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_block_adjust.py tests/test_fit.py -v`
Expected: PASS, 12 tests (6 from Task 2 still green).

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/fit.py tests/test_block_adjust.py
git commit -m "feat(autogeoref): block adjustment with fixed anchors, point-to-line and GCP observations"
```

---

## Task 11: GCP files — derive, write, read back, refit

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\gcp.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_gcp.py`

**Interfaces:**
- Consumes: `sheets.outline`, `sheets.turn_angles`, `anchors.read_points`, `anchors.write_points`, `fit.rigid_fit`, `fit.transform_points`, `edges.SegmentIndex`.
- Produces:
  - `gcp.Gcp` dataclass: `survey, corner_id, sheet_xy, map_xy, matched, residual_m, source`.
  - `gcp.corner_gcps(village, survey, theta, t, index, turn_min=20.0) -> list[Gcp]` — **one point per outline corner of every parcel**, whether or not imagery matched it, so a red parcel still has something for the team to drag (spec §9).
  - `gcp.write(village, survey, gcps, crs_wkt)` — writes `<survey>_parcels.geojson.points` only when no team-owned file exists, plus `gcp\<survey>_gcp.csv` always.
  - `gcp.team_points(village, survey) -> list[Gcp] | None` — the team's file when it is newer than the tool's CSV record.
  - `gcp.refit(village, survey) -> tuple[float, np.ndarray, float, int] | None` — `(theta, t, rms, n)` from the team's GCPs, scale fixed at 1.

- [ ] **Step 1: Write the failing test**

Create `tests/test_gcp.py`:

```python
import math
import numpy as np
import pytest
from autogeoref import gcp, edges, paths

VILLAGE = "35_04_077"
CRS = 'PROJCRS["WGS 84 / UTM zone 44N"]'


def test_every_corner_gets_a_point_even_with_no_image_match():
    empty = edges.SegmentIndex([], step=0.5)
    got = gcp.corner_gcps(VILLAGE, "46B", 0.0, np.array([0.0, 0.0]), empty)
    assert len(got) >= 4, "a red parcel must still have corners to drag"
    assert all(g.matched is False for g in got)
    assert all(g.source == "auto" for g in got)


def test_matched_flag_is_set_where_imagery_agrees():
    from shapely.geometry import Polygon
    from autogeoref import sheets, fit
    sheet = sheets.load_sheet(VILLAGE, "46B")
    v = sheets.outline(sheet)
    segs = []
    for i in range(len(v)):
        a, b = v[i], v[(i + 1) % len(v)]
        segs.append(edges.Segment(a, b, math.atan2(b[1] - a[1], b[0] - a[0]) % math.pi, math.dist(a, b)))
    idx = edges.SegmentIndex(segs, step=0.5)
    got = gcp.corner_gcps(VILLAGE, "46B", 0.0, np.array([0.0, 0.0]), idx)
    assert any(g.matched for g in got)


def test_write_then_read_round_trips_through_the_qgis_format(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / VILLAGE).mkdir(parents=True)
    rows = [gcp.Gcp("46B", 1, (0.0, 0.0), (392800.0, 1413300.0), True, 0.0, "auto"),
            gcp.Gcp("46B", 2, (30.0, 0.0), (392830.0, 1413300.0), True, 0.0, "auto")]
    gcp.write(VILLAGE, "46B", rows, CRS)
    assert paths.points_path(VILLAGE, "46B").exists()
    assert paths.gcp_csv(VILLAGE, "46B").exists()
    back = gcp.team_points(VILLAGE, "46B")
    assert back is None, "the tool's own file is not a team correction"


def test_refit_uses_the_teams_points_and_stays_rigid(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / VILLAGE).mkdir(parents=True)
    from autogeoref import anchors, fit
    src = np.array([(0.0, 0.0), (30.0, 0.0), (30.0, 20.0)])
    dst = fit.transform_points(src, 41.0, np.array([392800.0, 1413300.0]))
    anchors.write_points(paths.points_path(VILLAGE, "46B"),
                         [(d[0], d[1], s[0], s[1]) for s, d in zip(src, dst)], CRS)
    theta, t, rms, n = gcp.refit(VILLAGE, "46B")
    assert theta == pytest.approx(41.0, abs=0.01)
    assert t == pytest.approx(np.array([392800.0, 1413300.0]), abs=0.05)
    assert rms < 0.01 and n == 3


def test_refit_refuses_a_single_point(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / VILLAGE).mkdir(parents=True)
    from autogeoref import anchors
    anchors.write_points(paths.points_path(VILLAGE, "46B"), [(392800.0, 1413300.0, 0.0, 0.0)], CRS)
    assert gcp.refit(VILLAGE, "46B") is None, "one GCP cannot fix a pose"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_gcp.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.gcp'`

- [ ] **Step 3: Write `autogeoref/gcp.py`**

```python
"""Ground control points in the team's own format.

The team georeferences in the QGIS Georeferencer, which reads and writes
`<survey>_parcels.geojson.points`. The tool writes its automatic GCPs into that same file so the
team opens a sheet with the points already placed, drags the ones that are wrong, and saves.
Two dragged points fully determine a rigid pose, so there is no minimum-count rule for a refit.
"""
import csv
import datetime
import math
from dataclasses import dataclass

import numpy as np

from . import anchors, fit, paths, sheets

CORNER_TURN_MIN = 20.0        # degrees: what counts as a real corner
MATCH_DIST = 1.5              # metres
TEAM_SIGMA = 0.20             # metres, spec section 5 step 11


@dataclass
class Gcp:
    survey: str
    corner_id: int
    sheet_xy: tuple
    map_xy: tuple
    matched: bool
    residual_m: float
    source: str               # "auto" | "team"


def corner_gcps(village, survey, theta, t, index, turn_min=CORNER_TURN_MIN):
    """One GCP per real corner of the sheet outline at the given pose, matched flag from imagery."""
    v = sheets.outline(sheets.load_sheet(village, survey))
    turns = sheets.turn_angles(v)
    out = []
    corner_id = 0
    for i, p in enumerate(v):
        if turns[i] < turn_min:
            continue
        corner_id += 1
        ground = fit.transform_points([p], theta, t)[0]
        d, _b = index.query(np.array([ground]), max_dist=MATCH_DIST + 1.0)
        matched = bool(np.isfinite(d[0]) and d[0] <= MATCH_DIST)
        out.append(Gcp(survey, corner_id, (float(p[0]), float(p[1])),
                       (float(ground[0]), float(ground[1])), matched,
                       float(d[0]) if np.isfinite(d[0]) else float("nan"), "auto"))
    return out


def _csv_rows(village, survey):
    p = paths.gcp_csv(village, survey)
    if not p.exists():
        return []
    return list(csv.DictReader(p.open(encoding="utf-8")))


def write(village, survey, gcps, crs_wkt):
    """Write the durable CSV always; write the .points file only when the team owns no version."""
    p = paths.gcp_csv(village, survey)
    p.parent.mkdir(parents=True, exist_ok=True)
    existing = {int(r["corner_id"]): r for r in _csv_rows(village, survey) if r.get("source") == "team"}
    rows = []
    for g in gcps:
        if g.corner_id in existing:          # never overwrite a team row
            rows.append(existing[g.corner_id])
            continue
        rows.append({"survey": g.survey, "corner_id": g.corner_id,
                     "sheet_x": "%.4f" % g.sheet_xy[0], "sheet_y": "%.4f" % g.sheet_xy[1],
                     "map_x": "%.4f" % g.map_xy[0], "map_y": "%.4f" % g.map_xy[1],
                     "matched": int(bool(g.matched)),
                     "residual_m": "" if g.residual_m != g.residual_m else "%.3f" % g.residual_m,
                     "source": g.source,
                     "written": datetime.datetime.now().isoformat(timespec="seconds")})
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["survey", "corner_id", "sheet_x", "sheet_y", "map_x",
                                           "map_y", "matched", "residual_m", "source", "written"])
        w.writeheader()
        w.writerows(rows)
    pts = paths.points_path(village, survey)
    if team_points(village, survey) is None:
        anchors.write_points(pts, [(float(r["map_x"]), float(r["map_y"]),
                                    float(r["sheet_x"]), float(r["sheet_y"])) for r in rows], crs_wkt)
    return p


def team_points(village, survey):
    """The team's GCPs when their .points file is newer than the tool's CSV, else None."""
    pts = paths.points_path(village, survey)
    if not pts.exists():
        return None
    rec = paths.gcp_csv(village, survey)
    if rec.exists() and pts.stat().st_mtime <= rec.stat().st_mtime + 1:
        return None
    P = anchors.read_points(pts)
    if len(P) < 2:
        return None
    return [Gcp(survey, i + 1, (float(r[2]), float(r[3])), (float(r[0]), float(r[1])),
                True, float("nan"), "team") for i, r in enumerate(P)]


def refit(village, survey):
    """(theta, t, rms, n) from the team's GCPs, scale fixed at 1; None when there are fewer than two."""
    team = team_points(village, survey)
    if not team:
        return None
    P = np.array([g.sheet_xy for g in team], float)
    Q = np.array([g.map_xy for g in team], float)
    theta, t, rms, _mx = fit.rigid_fit(P, Q)
    return float(theta), t, float(rms), len(team)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_gcp.py -v`
Expected: PASS, 5 tests. `test_write_then_read_round_trips_through_the_qgis_format` relies on the CSV being written after the `.points` file, so `team_points` sees the tool's own file as not-newer; keep the CSV write before the `.points` write as coded.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/gcp.py tests/test_gcp.py
git commit -m "feat(autogeoref): GCPs in the QGIS Georeferencer format with a team-correction channel"
```

---

## Task 12: QGIS bridge and safe file replacement

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\qgis_bridge.py`
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\files.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_files.py`

**Interfaces:**
- Consumes: `paths.logs_dir`.
- Produces:
  - `qgis_bridge.available() -> bool`
  - `qgis_bridge.execute(code: str, timeout=170) -> dict` — runs PyQGIS in the open QGIS and returns `{"ok", "stdout", "error"}`; raises `QgisUnavailable` when the socket is down.
  - `qgis_bridge.QgisUnavailable` exception.
  - `files.held_layers(path) -> list[dict]` — `[{"id", "name", "editable", "modified", "group", "index"}]`, empty when QGIS is unreachable.
  - `files.release(path) -> list[dict]` — removes those layers, returning what to restore.
  - `files.restore(path, layers)` — re-adds them in their old group and position.
  - `files.safe_write(path, writer) -> bool` — release, delete stale `-wal`/`-shm`, write to `path.tmp` via `writer(tmp_path)`, `os.replace`, restore. Returns False and leaves the original untouched when a layer is in edit mode or the replace is denied.
  - `files.supersede(path, run_id) -> Path | None` — move to `_logs\georef_superseded_<run_id>\`, never delete.

- [ ] **Step 1: Write the failing test**

Create `tests/test_files.py`:

```python
import os
import geopandas as gpd
import pytest
from shapely.geometry import Polygon
from autogeoref import files, paths, qgis_bridge

SQ = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])


def _write(path, value):
    gpd.GeoDataFrame({"v": [value]}, geometry=[SQ], crs="EPSG:32644").to_file(path, driver="GPKG")


def test_safe_write_replaces_a_file_that_nothing_holds(tmp_path):
    target = tmp_path / "x.gpkg"
    _write(target, 1)
    ok = files.safe_write(target, lambda p: _write(p, 2))
    assert ok
    assert int(gpd.read_file(target)["v"].iloc[0]) == 2


def test_safe_write_leaves_the_original_when_the_writer_raises(tmp_path):
    target = tmp_path / "x.gpkg"
    _write(target, 1)

    def boom(p):
        raise RuntimeError("writer failed")

    assert files.safe_write(target, boom) is False
    assert int(gpd.read_file(target)["v"].iloc[0]) == 1
    assert not (tmp_path / "x.gpkg.tmp").exists()


def test_safe_write_clears_stale_wal_sidecars(tmp_path):
    target = tmp_path / "x.gpkg"
    _write(target, 1)
    (tmp_path / "x.gpkg-wal").write_bytes(b"stale")
    (tmp_path / "x.gpkg-shm").write_bytes(b"stale")
    assert files.safe_write(target, lambda p: _write(p, 2))
    assert not (tmp_path / "x.gpkg-wal").exists()


def test_supersede_moves_and_never_deletes(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    target = tmp_path / "old.gpkg"
    _write(target, 1)
    moved = files.supersede(target, "20260918")
    assert moved is not None and moved.exists() and not target.exists()
    assert "georef_superseded_20260918" in str(moved)


def test_held_layers_is_empty_when_qgis_is_unreachable(monkeypatch, tmp_path):
    monkeypatch.setattr(qgis_bridge, "available", lambda: False)
    assert files.held_layers(tmp_path / "anything.gpkg") == []


@pytest.mark.skipif(not qgis_bridge.available(), reason="QGIS is not running")
def test_held_layers_finds_a_file_the_open_project_holds():
    p = paths.vector_dir("35_04_077") / "171_parcels_modified1.gpkg"
    got = files.held_layers(p)
    assert isinstance(got, list)
    for layer in got:
        assert {"id", "name", "editable"} <= set(layer)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_files.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.qgis_bridge'`

- [ ] **Step 3: Write `autogeoref/qgis_bridge.py`**

QGIS is reached through the `qgis_mcp` plugin's socket on `localhost:9876`; messages are length-prefixed JSON. The cached client used during development lives at `C:\Users\FAI-Akash\AppData\Local\uv\cache\archive-v0\6UbEkjU0r6KaKsVQ\Lib\site-packages`. Import it if present, else speak the protocol directly.

```python
"""Run PyQGIS inside the user's open QGIS over the qgis_mcp socket. Never import qgis in-process."""
import json
import socket
import struct
import sys

HOST, PORT = "127.0.0.1", 9876
CLIENT_SITE = r"C:\Users\FAI-Akash\AppData\Local\uv\cache\archive-v0\6UbEkjU0r6KaKsVQ\Lib\site-packages"


class QgisUnavailable(RuntimeError):
    pass


def available(timeout=2.0):
    try:
        with socket.create_connection((HOST, PORT), timeout=timeout):
            return True
    except OSError:
        return False


def _send(sock, payload):
    data = json.dumps(payload).encode("utf-8")
    sock.sendall(struct.pack(">I", len(data)) + data)


def _recv(sock):
    head = b""
    while len(head) < 4:
        chunk = sock.recv(4 - len(head))
        if not chunk:
            raise QgisUnavailable("QGIS closed the connection")
        head += chunk
    n = struct.unpack(">I", head)[0]
    body = b""
    while len(body) < n:
        chunk = sock.recv(min(65536, n - len(body)))
        if not chunk:
            raise QgisUnavailable("QGIS closed the connection")
        body += chunk
    return json.loads(body.decode("utf-8"))


def execute(code, timeout=170):
    """Run `code` in QGIS. Returns {'ok': bool, 'stdout': str, 'error': str|None}."""
    if CLIENT_SITE not in sys.path:
        sys.path.append(CLIENT_SITE)
    try:
        from qgis_mcp.client import QgisMCPClient          # noqa: F401  (preferred path)
        client = QgisMCPClient()
        client.connect()
        try:
            reply = client.send_command("execute_code", {"code": code, "timeout": timeout},
                                        timeout=timeout + 10)
        finally:
            client.disconnect()
    except ImportError:
        try:
            with socket.create_connection((HOST, PORT), timeout=timeout + 10) as sock:
                _send(sock, {"type": "execute_code", "params": {"code": code, "timeout": timeout}})
                reply = _recv(sock)
        except OSError as exc:
            raise QgisUnavailable(str(exc))
    result = reply.get("result", reply)
    return {"ok": bool(result.get("executed")), "stdout": result.get("stdout", ""),
            "error": result.get("error") or result.get("message")}


def json_result(code, timeout=170):
    """Run code whose last stdout line is JSON, and return the parsed object."""
    out = execute(code, timeout=timeout)
    if not out["ok"]:
        raise QgisUnavailable(out["error"] or "QGIS refused the code")
    lines = [ln for ln in out["stdout"].splitlines() if ln.strip()]
    return json.loads(lines[-1])
```

- [ ] **Step 4: Write `autogeoref/files.py`**

```python
"""Writing GeoPackages that QGIS may be holding open, and superseding files instead of deleting.

Measured on this machine with QGIS 3.40.10: while a QgsVectorLayer holds a GeoPackage, os.replace
and os.remove both raise PermissionError(13), even when the layer is not in edit mode and not in
the layer tree. The handle is released as soon as the layer object is destroyed.
"""
import json
import os
import shutil
from pathlib import Path

from . import paths, qgis_bridge

_FIND = """
import json
from qgis.core import QgsProject, QgsProviderRegistry
target = %r.replace("\\\\", "/").lower()
proj = QgsProject.instance(); root = proj.layerTreeRoot(); out = []
for lyr in list(proj.mapLayers().values()):
    try:
        src = QgsProviderRegistry.instance().decodeUri(lyr.providerType(), lyr.source()).get("path", "")
    except Exception:
        continue
    if str(src).replace("\\\\", "/").lower() != target:
        continue
    node = root.findLayer(lyr.id()); parent = node.parent() if node else None
    out.append({"id": lyr.id(), "name": lyr.name(), "editable": lyr.isEditable(),
                "modified": lyr.isModified(), "uri": lyr.source(),
                "group": parent.name() if parent is not None and hasattr(parent, "name") else "",
                "index": parent.children().index(node) if parent is not None and node else -1})
print(json.dumps(out))
"""

_REMOVE = """
import gc, json
from qgis.core import QgsProject
proj = QgsProject.instance()
removed = []
for lid in %r:
    lyr = proj.mapLayer(lid)
    if lyr is not None:
        removed.append(lid); proj.removeMapLayer(lid)
gc.collect()
print(json.dumps(removed))
"""

_RESTORE = """
import json
from qgis.core import QgsProject, QgsVectorLayer, QgsLayerTreeGroup
proj = QgsProject.instance(); root = proj.layerTreeRoot(); added = []
for spec in %s:
    lyr = QgsVectorLayer(spec["uri"], spec["name"], "ogr")
    if not lyr.isValid():
        continue
    proj.addMapLayer(lyr, False)
    parent = root
    if spec.get("group"):
        found = root.findGroup(spec["group"])
        parent = found if found is not None else root
    node = parent.insertLayer(max(0, spec.get("index", 0)), lyr)
    added.append(spec["name"])
print(json.dumps(added))
"""


def held_layers(path):
    if not qgis_bridge.available():
        return []
    try:
        return qgis_bridge.json_result(_FIND % str(Path(path).resolve()))
    except Exception:
        return []


def release(path):
    held = held_layers(path)
    if not held:
        return []
    if any(h["editable"] for h in held):
        raise PermissionError("save or stop editing %s in QGIS first"
                              % ", ".join(h["name"] for h in held if h["editable"]))
    qgis_bridge.json_result(_REMOVE % [h["id"] for h in held])
    return held


def restore(path, layers):
    if not layers or not qgis_bridge.available():
        return
    try:
        qgis_bridge.json_result(_RESTORE % json.dumps(layers))
    except Exception:
        pass


def safe_write(path, writer):
    """writer(tmp_path) must create a complete file; returns False and keeps the original on failure."""
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    try:
        held = release(path)
    except PermissionError:
        return False
    try:
        if tmp.exists():
            tmp.unlink()
        writer(tmp)
        for side in ("-wal", "-shm"):
            p = Path(str(path) + side)
            if p.exists():
                p.unlink()
        os.replace(tmp, path)
        return True
    except Exception:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        return False
    finally:
        restore(path, held)


def supersede(path, run_id):
    """Move a file the tool is replacing into _logs; nothing is ever deleted."""
    path = Path(path)
    if not path.exists():
        return None
    dest = paths.logs_dir() / ("georef_superseded_%s" % run_id)
    dest.mkdir(parents=True, exist_ok=True)
    try:
        release(path)
    except PermissionError:
        return None
    target = dest / path.name
    shutil.move(str(path), str(target))
    return target
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_files.py -v`
Expected: PASS, 6 tests (the last is skipped when QGIS is not running).

- [ ] **Step 6: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/qgis_bridge.py autogeoref/files.py tests/test_files.py
git commit -m "feat(autogeoref): QGIS socket bridge and lock-safe GeoPackage replacement"
```

---

## Task 13: Topology between free parcels, anchors untouched

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\topology.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_topology.py`

**Interfaces:**
- Consumes: `sheets.dissolve`, `paths`.
- Produces:
  - `topology.check(geoms: dict[str, Polygon]) -> dict` — `{"overlap_pairs", "overlap_sqm", "overlaps", "gaps", "gap_sqm", "invalid", "union_parts"}`.
  - `topology.fix(geoms, movable: set[str], rail: set[str], snap=0.30) -> (dict[str, Polygon], report)` — snaps shared boundaries, clips remaining overlaps (a parcel yields to railway land; otherwise the larger polygon yields), fills enclosed gaps into the neighbour with the longest shared edge. Members of `geoms` that are not in `movable` are never changed. Railway-land vs parcel overlaps are **listed, never clipped**.
  - `topology.report_rail_conflicts(village, report)` — writes `rail_conflicts.csv`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_topology.py`:

```python
import pytest
from shapely.geometry import Polygon
from autogeoref import topology

A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
B_OVERLAP = Polygon([(9.7, 0), (20, 0), (20, 10), (9.7, 10)])      # 3 sqm overlap with A
B_GAP = Polygon([(10.3, 0), (20, 0), (20, 10), (10.3, 10)])        # 3 sqm gap from A


def test_check_counts_overlaps_and_gaps():
    got = topology.check({"A": A, "B": B_OVERLAP})
    assert got["overlap_pairs"] == 1
    assert got["overlap_sqm"] == pytest.approx(3.0, abs=0.01)


def test_fix_closes_a_sub_metre_overlap_between_two_free_parcels():
    out, rep = topology.fix({"A": A, "B": B_OVERLAP}, movable={"A", "B"}, rail=set())
    assert topology.check(out)["overlap_sqm"] < 0.02
    assert rep["clips"] or rep["snapped"]


def test_fix_closes_a_sub_metre_gap():
    out, rep = topology.fix({"A": A, "B": B_GAP}, movable={"A", "B"}, rail=set())
    assert out["A"].distance(out["B"]) < 0.01


def test_an_anchor_is_never_modified():
    out, _rep = topology.fix({"A": A, "B": B_OVERLAP}, movable={"B"}, rail=set())
    assert out["A"].equals(A), "A is not movable, so it must come back untouched"


def test_railway_overlap_is_listed_not_clipped():
    out, rep = topology.fix({"169": A, "42B": B_OVERLAP}, movable={"169", "42B"}, rail={"169"})
    assert out["169"].area == pytest.approx(A.area, abs=0.01), "railway land keeps its shape"
    assert rep["rail_conflicts"], "the conflict is reported for the team"
    assert rep["rail_conflicts"][0][2] == pytest.approx(3.0, abs=0.5)


def test_fix_is_idempotent():
    once, _r1 = topology.fix({"A": A, "B": B_OVERLAP}, movable={"A", "B"}, rail=set())
    twice, _r2 = topology.fix(once, movable={"A", "B"}, rail=set())
    assert topology.check(twice)["overlap_sqm"] < 0.02
    assert abs(twice["B"].area - once["B"].area) < 0.01
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_topology.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.topology'`

- [ ] **Step 3: Write `autogeoref/topology.py`**

The snap is done with shapely rather than QGIS's `native:snapgeometries` so the module is testable without QGIS running. `shapely.snap` moves vertices onto nearby vertices of the reference; for a sub-metre offset along a shared straight edge that is not enough on its own, so the clip and fill stages follow, exactly as the 2026-09-18 in-place fix did (it left 0 overlaps and 0 gaps across 17 parcels with a maximum vertex move of 1.03 m).

```python
"""Topology between parcels the tool may move. Anchors and team-fixed parcels are never touched."""
import csv
import itertools

from shapely.geometry import Polygon
from shapely.ops import snap, unary_union

from . import paths

MIN_AREA = 0.02        # square metres: below this an overlap or gap is numerical noise


def _largest(geom):
    if geom.is_empty:
        return geom
    if geom.geom_type == "MultiPolygon":
        return max(geom.geoms, key=lambda q: q.area)
    if geom.geom_type == "GeometryCollection":
        polys = [g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
        return _largest(unary_union(polys)) if polys else Polygon()
    return geom


def check(geoms):
    keys = sorted(geoms, key=paths.survey_sort_key)
    overlaps = []
    for a, b in itertools.combinations(keys, 2):
        area = geoms[a].intersection(geoms[b]).area
        if area > MIN_AREA:
            overlaps.append([a, b, round(area, 2)])
    union = unary_union([g.buffer(0) for g in geoms.values()])
    parts = union.geoms if union.geom_type == "MultiPolygon" else [union]
    gaps = [round(Polygon(r).area, 2) for p in parts for r in p.interiors if Polygon(r).area > 0.01]
    return {"overlap_pairs": len(overlaps), "overlap_sqm": round(sum(o[2] for o in overlaps), 2),
            "overlaps": overlaps, "gaps": len(gaps), "gap_sqm": round(sum(gaps), 2),
            "gap_list": sorted(gaps, reverse=True)[:15],
            "invalid": [k for k, g in geoms.items() if not g.is_valid],
            "union_parts": len(parts)}


def fix(geoms, movable, rail, snap_tol=0.30):
    out = {k: g.buffer(0) for k, g in geoms.items()}
    report = {"snapped": [], "clips": [], "filled": [], "rail_conflicts": []}
    keys = sorted(out, key=paths.survey_sort_key)

    # 1. snap the movable parcels onto their neighbours
    for a in keys:
        if a not in movable:
            continue
        for b in keys:
            if a == b or out[a].distance(out[b]) > snap_tol:
                continue
            moved = snap(out[a], out[b], snap_tol).buffer(0)
            if not moved.equals(out[a]):
                out[a] = _largest(moved)
                report["snapped"].append([a, b])

    # 2. clip remaining overlaps
    for _round in range(6):
        changed = 0
        for a, b in itertools.combinations(keys, 2):
            area = out[a].intersection(out[b]).area
            if area <= MIN_AREA:
                continue
            ra, rb = a in rail, b in rail
            if ra != rb:
                pair = (a, b, round(area, 2))
                if pair not in report["rail_conflicts"]:
                    report["rail_conflicts"].append(pair)
                continue
            big, small = (a, b) if out[a].area >= out[b].area else (b, a)
            if big not in movable:
                big, small = small, big
            if big not in movable:
                continue                                  # both frozen: leave it, report only
            out[big] = _largest(out[big].difference(out[small]))
            report["clips"].append([big, "clipped by", small, round(area, 2)])
            changed += 1
        if not changed:
            break

    # 3. fill enclosed gaps into the neighbour with the longest shared edge
    union = unary_union([g.buffer(0) for g in out.values()])
    parts = union.geoms if union.geom_type == "MultiPolygon" else [union]
    for part in parts:
        for ring in part.interiors:
            gap = Polygon(ring)
            if gap.area <= 0.01:
                continue
            best, score = None, 0.0
            for k in movable:
                s = out[k].buffer(0.02).intersection(gap).area
                if s > score:
                    best, score = k, s
            if best:
                out[best] = _largest(unary_union([out[best], gap]))
                report["filled"].append([round(gap.area, 2), "into", best])
    return out, report


def report_rail_conflicts(village, report):
    rows = report.get("rail_conflicts", [])
    p = paths.vector_dir(village) / "rail_conflicts.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["railway_or_parcel_a", "railway_or_parcel_b", "overlap_sqm", "note"])
        for a, b, area in rows:
            w.writerow([a, b, area, "railway land carved out of the older parcel sheet; decide in review"])
    return p
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_topology.py -v`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/topology.py tests/test_topology.py
git commit -m "feat(autogeoref): topology fix for free parcels with railway conflicts reported"
```

---

## Task 14: Status CSV, worklist and the tracker

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\review.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_review.py`

**Spec amendment (user instruction, 2026-09-18):** "Don't modify my style, colours and other stuffs, just update the sheet is enough." The tracker is therefore **not written by default**. The tool's per-parcel verdicts go to `worklist.csv`, which the team opens in QGIS as a delimited-text layer. `--tracker` may add the three tool-owned columns, and it does so by patching the workbook XML — never by round-tripping through openpyxl, which silently deletes the x14 dropdown validations (`Lists!$A$2:$A$8` for Status, `Team!$A$2:$A$31` for Assigned To). The Status column stays the team's.

**Interfaces:**
- Consumes: `paths.status_path`, `paths.TRACKER`, `qgis_bridge`.
- Produces:
  - `review.STATUS_COLUMNS` — the ordered column list.
  - `review.write_status(village, rows) -> Path`
  - `review.read_status(village) -> list[dict]`
  - `review.tool_written_fingerprints(village) -> set[str]` — `fp_placed` and `fp_final` of every row, so Task 3 can tell the tool's files from the team's.
  - `review.write_worklist(rows_by_village) -> Path` — corridor-wide `FMB_Vector\worklist.csv`.
  - `review.update_tracker(rows, workbook=paths.TRACKER) -> str` — XML patch adding `Auto colour`, `Auto note`, `Auto run`; returns `"applied"`, `"locked"` or `"no change"`; never touches Status, styles, widths, validations or conditional formatting.
  - `review.build_group(village, rows, project_path=None) -> dict` — QGIS review group coloured by verdict; non-fatal when QGIS is unreachable or a different project is open.

- [ ] **Step 1: Write the failing test**

Create `tests/test_review.py`:

```python
import csv
import zipfile
import pytest
from autogeoref import review, paths

ROWS = [{"survey": "47B", "status": "placed", "method": "anchors", "colour": "green",
         "confidence": 88, "share": 0.62, "margin": 0.20, "boundary_rms_m": 0.21,
         "puvi_reference_m": 2.3, "notes": "", "fp_placed": "aa" * 32, "fp_final": "bb" * 32},
        {"survey": "46B", "status": "review", "method": "puvi-only", "colour": "red",
         "confidence": 30, "share": 0.09, "margin": 0.01, "boundary_rms_m": "",
         "puvi_reference_m": 8.1, "notes": "ambiguous", "fp_placed": "", "fp_final": ""}]


def test_status_csv_round_trips(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / "35_04_077").mkdir(parents=True)
    p = review.write_status("35_04_077", ROWS)
    assert p.exists()
    back = review.read_status("35_04_077")
    assert [r["survey"] for r in back] == ["46B", "47B"] or [r["survey"] for r in back] == ["47B", "46B"]
    assert set(review.STATUS_COLUMNS) <= set(back[0])


def test_tool_written_fingerprints_collects_both_columns(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / "35_04_077").mkdir(parents=True)
    review.write_status("35_04_077", ROWS)
    fps = review.tool_written_fingerprints("35_04_077")
    assert "aa" * 32 in fps and "bb" * 32 in fps and "" not in fps


def test_worklist_is_one_file_for_the_whole_corridor(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector").mkdir(parents=True)
    p = review.write_worklist({"35_04_077": ROWS, "35_04_074": []})
    rows = list(csv.DictReader(p.open(encoding="utf-8")))
    assert len(rows) == 2 and rows[0]["village_code"] == "35_04_077"


def test_tracker_update_keeps_the_dropdowns_and_formatting(tmp_path):
    src = paths.TRACKER
    if not src.exists():
        pytest.skip("tracker not present")
    work = tmp_path / "Tracker.xlsx"
    work.write_bytes(src.read_bytes())
    before = zipfile.ZipFile(work).read("xl/worksheets/sheet1.xml").decode("utf-8")
    review.update_tracker([{"village_code": "35_04_077", "survey": "47B", "colour": "green",
                            "notes": "anchored", "run": "2026-09-18T10:00:00"}], workbook=work)
    after = zipfile.ZipFile(work).read("xl/worksheets/sheet1.xml").decode("utf-8")
    assert after.count("<x14:dataValidation ") == before.count("<x14:dataValidation ")
    assert "<conditionalFormatting" in after
    assert before.count("<cfRule") == after.count("<cfRule")


def test_tracker_update_never_writes_the_status_column(tmp_path):
    src = paths.TRACKER
    if not src.exists():
        pytest.skip("tracker not present")
    import openpyxl
    work = tmp_path / "Tracker.xlsx"
    work.write_bytes(src.read_bytes())
    before = [r[9] for r in openpyxl.load_workbook(work, data_only=True)["Tracker"].iter_rows(values_only=True)]
    review.update_tracker([{"village_code": "35_04_077", "survey": "47B", "colour": "red",
                            "notes": "x", "run": "2026-09-18T10:00:00"}], workbook=work)
    after = [r[9] for r in openpyxl.load_workbook(work, data_only=True)["Tracker"].iter_rows(values_only=True)]
    assert before == after, "Status is the team's column"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_review.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.review'`

- [ ] **Step 3: Write `autogeoref/review.py`**

```python
"""Everything the team reads: the status CSV, the corridor worklist, the tracker and the QGIS group."""
import csv
import json
import re
import shutil
import zipfile
from pathlib import Path

from . import paths, qgis_bridge

STATUS_COLUMNS = ["survey", "status", "method", "colour", "confidence", "share", "margin",
                  "observable", "boundary_rms_m", "matched_length_m", "neighbours", "heading_deg",
                  "shift_m", "sigma_pos_m", "sigma_head_deg", "puvi_reference_m", "puvi_key",
                  "anchor_file", "puvi_trusted", "puvi_ratio", "notes", "file", "fp_placed", "fp_final", "run"]
TRACKER_COLUMNS = ["Auto colour", "Auto note", "Auto run"]


def write_status(village, rows):
    p = paths.status_path(village)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=STATUS_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: paths.survey_sort_key(r["survey"])):
            w.writerow({c: r.get(c, "") for c in STATUS_COLUMNS})
    return p


def read_status(village):
    p = paths.status_path(village)
    return list(csv.DictReader(p.open(encoding="utf-8"))) if p.exists() else []


def tool_written_fingerprints(village):
    out = set()
    for r in read_status(village):
        for c in ("fp_placed", "fp_final"):
            if r.get(c):
                out.add(r[c])
    return out


def write_worklist(rows_by_village):
    p = paths.PROJECT / "FMB_Vector" / "worklist.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    cols = ["village_code", "survey", "colour", "confidence", "method", "notes", "file"]
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for village in sorted(rows_by_village):
            for r in sorted(rows_by_village[village], key=lambda r: paths.survey_sort_key(r["survey"])):
                w.writerow(dict({c: r.get(c, "") for c in cols}, village_code=village))
    return p


# --------------------------------------------------------------------------- tracker (XML patch)
def _escape(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def update_tracker(rows, workbook=None):
    """Add the three tool-owned columns, keyed by village code + survey number.

    The workbook is patched as a zip. openpyxl is never used to save it: it deletes the x14
    dataValidations that hold the Status and Assigned To dropdowns.
    """
    workbook = Path(workbook or paths.TRACKER)
    if not workbook.exists():
        return "no workbook"
    lock = workbook.parent / ("~$" + workbook.name)
    if lock.exists():
        return "locked"
    zin = zipfile.ZipFile(workbook)
    parts = {n: zin.read(n) for n in zin.namelist()}
    zin.close()
    sheet = parts["xl/worksheets/sheet1.xml"].decode("utf-8")
    shared = parts["xl/sharedStrings.xml"].decode("utf-8")
    values = ["".join(re.findall(r"<t[^>]*>(.*?)</t>", si, re.S))
              for si in re.findall(r"<si>(.*?)</si>", shared, re.S)]

    # locate the header columns we own, adding them after the last used column when absent
    header = re.search(r'<row r="1".*?</row>', sheet, re.S).group(0)
    used = re.findall(r'<c r="([A-Z]+)1"', header)
    last = used[-1]
    def next_col(col):
        return chr(ord(col) + 1) if len(col) == 1 and col < "Z" else col + "A"
    # build the survey -> row index map from column H, and the village code from column G
    rowmap = {}
    for m in re.finditer(r'<row r="(\d+)".*?</row>', sheet, re.S):
        rid, body = int(m.group(1)), m.group(0)
        if rid == 1:
            continue
        g = re.search(r'<c r="G%d"[^>]*><v>(\d+)</v></c>' % rid, body)
        h = re.search(r'<c r="H%d"[^>]*t="s"><v>(\d+)</v></c>' % rid, body)
        if not (g and h):
            continue
        rowmap[(g.group(1), values[int(h.group(1))])] = rid

    new_strings = []

    def sref(text):
        if text in values:
            return values.index(text)
        if text in new_strings:
            return len(values) + new_strings.index(text)
        new_strings.append(text)
        return len(values) + len(new_strings) - 1

    cols = {}
    cursor = last
    for name in TRACKER_COLUMNS:
        cursor = next_col(cursor)
        cols[name] = cursor
        sheet = sheet.replace("</row>", '<c r="%s1" s="7" t="s"><v>%d</v></c></row>' % (cursor, sref(name)), 1)

    changed = 0
    for r in rows:
        village_no = str(r.get("village_code", "")).split("_")[-1].lstrip("0") or "0"
        key = (village_no, str(r.get("survey", "")))
        rid = rowmap.get(key)
        if rid is None:
            continue
        cells = "".join('<c r="%s%d" s="9" t="s"><v>%d</v></c>' % (cols[name], rid, sref(str(val)))
                        for name, val in (("Auto colour", r.get("colour", "")),
                                          ("Auto note", r.get("notes", "")),
                                          ("Auto run", r.get("run", ""))) if str(val))
        if not cells:
            continue
        sheet = re.sub(r'(<row r="%d".*?)</row>' % rid, lambda m: m.group(1) + cells + "</row>", sheet,
                       count=1, flags=re.S)
        changed += 1
    if not changed:
        return "no change"
    if new_strings:
        n0 = int(re.search(r'uniqueCount="(\d+)"', shared).group(1))
        c0 = int(re.search(r'\bcount="(\d+)"', shared).group(1))
        shared = shared.replace("</sst>", "".join("<si><t>%s</t></si>" % _escape(s) for s in new_strings) + "</sst>")
        shared = shared.replace('count="%d" uniqueCount="%d"' % (c0, n0),
                                'count="%d" uniqueCount="%d"' % (c0 + changed * 3, n0 + len(new_strings)), 1)
        parts["xl/sharedStrings.xml"] = shared.encode("utf-8")
    sheet = re.sub(r'<dimension ref="A1:[A-Z]+(\d+)"/>',
                   lambda m: '<dimension ref="A1:%s%s"/>' % (cols[TRACKER_COLUMNS[-1]], m.group(1)), sheet, count=1)
    parts["xl/worksheets/sheet1.xml"] = sheet.encode("utf-8")
    parts["xl/workbook.xml"] = parts["xl/workbook.xml"].decode("utf-8").replace(
        "<calcPr ", '<calcPr fullCalcOnLoad="1" ', 1).encode("utf-8")
    backup = paths.logs_dir() / "tracker_backups"
    backup.mkdir(parents=True, exist_ok=True)
    shutil.copy2(workbook, backup / workbook.name.replace(".xlsx", "_before_autogeoref.xlsx"))
    tmp = workbook.with_suffix(".xlsx.tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
        for name, data in parts.items():
            zo.writestr(name, data)
    try:
        tmp.replace(workbook)
    except PermissionError:
        tmp.unlink(missing_ok=True)
        return "locked"
    return "applied"


# --------------------------------------------------------------------------- QGIS review group
_GROUP = """
import json, glob, os
from qgis.core import (QgsProject, QgsVectorLayer, QgsRasterLayer, QgsFillSymbol,
                       QgsSingleSymbolRenderer, QgsLayerTreeGroup)
proj = QgsProject.instance()
want = %r
if want and os.path.normcase(proj.fileName()) != os.path.normcase(want):
    print(json.dumps({"skipped": "a different project is open", "open": proj.fileName()}))
else:
    root = proj.layerTreeRoot(); name = "Georef review - " + %r
    old = root.findGroup(name)
    if old:
        for n in old.findLayers():
            proj.removeMapLayer(n.layerId())
        root.removeChildNode(old)
    grp = QgsLayerTreeGroup(name); root.insertChildNode(0, grp); added = []
    for path, colour, label in %s:
        lyr = QgsVectorLayer(path + "|layername=" + os.path.basename(path)[:-5], label, "ogr")
        if not lyr.isValid():
            continue
        lyr.setRenderer(QgsSingleSymbolRenderer(QgsFillSymbol.createSimple(
            {"color": colour + "40", "outline_color": colour, "outline_width": "0.6"})))
        proj.addMapLayer(lyr, False); grp.addLayer(lyr); added.append(label)
    print(json.dumps({"group": name, "layers": added}))
"""

COLOURS = {"green": "#19e68c", "amber": "#ffb000", "red": "#ff4136", "anchor": "#000000"}


def build_group(village, rows, project_path=None):
    """Load the run's outputs into QGIS. Never fatal: the files and CSVs are written regardless."""
    if not qgis_bridge.available():
        return {"skipped": "QGIS not reachable; rerun with --review when QGIS is open"}
    spec = [(str(paths.output_path(village, r["survey"])), COLOURS.get(r.get("colour"), "#888888"),
             "%s (%s)" % (r["survey"], r.get("colour", ""))) for r in rows
            if paths.output_path(village, r["survey"]).exists()]
    try:
        return qgis_bridge.json_result(_GROUP % (project_path or "", village, json.dumps(spec)))
    except Exception as exc:
        return {"skipped": str(exc)}
```


- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_review.py -v`
Expected: PASS, 5 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/review.py tests/test_review.py
git commit -m "feat(autogeoref): status CSV, corridor worklist, non-destructive tracker patch, QGIS group"
```

---

## Task 15: The per-village engine

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\engine.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_engine.py`

**Interfaces:**
- Consumes: every module above.
- Produces:
  - `engine.run(village, do_raster=False, do_topology=False, do_review=True, project=None, root=None) -> dict` — the pipeline of spec §5 steps 1–13.
  - `engine.write_parcels(village, survey, theta, t, row) -> Path` — writes `<survey>_parcels_modified.gpkg` with `parcels` and `edges` layers, the §11 measurement attributes, and asserts rigidity.
  - `engine.colour_of(row) -> str` — the green / amber / red rule of spec §5 step 9.

- [ ] **Step 1: Write the failing test**

Create `tests/test_engine.py`:

```python
import numpy as np
import pytest
from autogeoref import engine


def test_colour_rule_green_needs_an_anchor_and_observability():
    assert engine.colour_of({"boundary_rms_m": 0.20, "anchor_rms_m": 0.3, "share": 0.65,
                             "observable": True, "n_neighbours": 2, "method": "anchors",
                             "ambiguous": False}) == "green"


def test_colour_rule_amber_when_only_one_condition_is_short():
    assert engine.colour_of({"boundary_rms_m": 0.20, "anchor_rms_m": 0.3, "share": 0.20,
                             "observable": False, "n_neighbours": 1, "method": "anchors",
                             "ambiguous": False}) == "amber"


def test_colour_rule_red_for_an_ambiguous_image_pose():
    assert engine.colour_of({"boundary_rms_m": None, "anchor_rms_m": None, "share": 0.55,
                             "observable": True, "n_neighbours": 0, "method": "image",
                             "ambiguous": True}) == "red"


def test_colour_rule_red_for_puvi_only():
    assert engine.colour_of({"boundary_rms_m": None, "anchor_rms_m": None, "share": 0.0,
                             "observable": False, "n_neighbours": 0, "method": "puvi-only",
                             "ambiguous": True}) == "red"


def test_write_parcels_keeps_fmb_dimensions_exactly(tmp_path, monkeypatch):
    import geopandas as gpd
    from autogeoref import paths, sheets
    monkeypatch.setattr(paths, "PROJECT", paths.PROJECT)      # real inputs, temp output below
    out = engine.write_parcels("35_04_077", "46B", 31.0, np.array([392800.0, 1413300.0]),
                               {"colour": "green", "method": "anchors", "confidence": 90},
                               out_dir=tmp_path)
    got = gpd.read_file(out, layer="parcels")
    src = sheets.load_sheet("35_04_077", "46B")
    assert len(got) == len(src)
    assert float(got.geometry.area.sum()) == pytest.approx(sum(g.area for _p, g in src), abs=1e-6)
    edges = gpd.read_file(out, layer="edges")
    assert {"length_m", "bearing_grid_deg", "bearing_true_deg"} <= set(edges.columns)
    assert float(edges["length_m"].max()) > 1.0


def test_run_on_a_village_where_every_parcel_is_an_anchor_changes_nothing(tmp_path):
    """Kizhikaranai is fully hand placed: the engine must place nothing and touch no manual file."""
    before = {p.name: p.stat().st_mtime for p in
              (engine.paths.vector_dir("35_04_077")).glob("*_parcels_modified*.gpkg")}
    out = engine.run("35_04_077", do_raster=False, do_topology=False, do_review=False)
    after = {p.name: p.stat().st_mtime for p in
             (engine.paths.vector_dir("35_04_077")).glob("*_parcels_modified*.gpkg")}
    assert before == after, "anchors are frozen"
    assert out["anchors"] == 15 and out["placed"] == 0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_engine.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.engine'`

- [ ] **Step 3: Write `autogeoref/engine.py`**

```python
"""The per-village pipeline of spec section 5."""
import datetime
import math

import geopandas as gpd
import numpy as np
from pyproj import Proj
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from . import (align, anchors, edges as edgemod, files, fit, gcp, match, neighbours,
               paths, raster, review, sheets, topology, window)

ACRE = 0.000247105381
SIGMA_AUTO = 0.30
SIGMA_IMAGE = 1.0
SIGMA_START = 10.0
IMAGE_MAX_SHIFT = 1.5        # imagery may move an anchored pose by at most this
IMAGE_MAX_TURN = 2.0
GREEN_SHARE = 0.60           # recalibrated in Task 17 and written back here
AMBER_SHARE = 0.35


def colour_of(row):
    if row.get("method") == "puvi-only" or row.get("ambiguous"):
        return "red"
    rms = row.get("boundary_rms_m")
    allow = SIGMA_AUTO + (row.get("anchor_rms_m") or 0.0)
    anchored = rms is not None and rms <= allow
    strong_image = (row.get("share") or 0.0) >= GREEN_SHARE and row.get("observable")
    if anchored and (strong_image or (row.get("n_neighbours") or 0) >= 2):
        return "green"
    if anchored or strong_image or (row.get("share") or 0.0) >= AMBER_SHARE:
        return "amber"
    return "red"


def write_parcels(village, survey, theta, t, row, out_dir=None):
    """<survey>_parcels_modified.gpkg with parcels + edges layers and the section 11 attributes."""
    sheet = sheets.load_sheet(village, survey)
    proj = Proj("EPSG:32644")
    run = row.get("run") or datetime.datetime.now().isoformat(timespec="seconds")
    parcels, edge_rows = [], []
    for props, geom in sheet:
        placed = fit.apply_pose(geom, theta, t)
        fit.assert_rigid(geom, placed)
        area = placed.area
        rec = dict(props)
        rec.update({"area_sqm": round(area, 3), "area_are": round(area / 100, 4),
                    "area_hect": round(area / 1e4, 6), "area_acre": round(area * ACRE, 5),
                    "area_cent": round(area * ACRE * 100, 3), "perimeter_m": round(placed.length, 3),
                    "fmb_area_sqm": round(geom.area, 3), "fmb_area_acre": round(geom.area * ACRE, 5),
                    "fmb_perimeter_m": round(geom.length, 3),
                    "georef_method": row.get("method"), "georef_colour": row.get("colour"),
                    "georef_confidence": row.get("confidence"), "georef_residual_m": row.get("boundary_rms_m"),
                    "anchored_to": row.get("neighbours"), "heading_deg": round(theta, 3),
                    "puvi_reference_m": row.get("puvi_reference_m"), "georef_notes": row.get("notes"),
                    "georef_run": run})
        parcels.append({**rec, "geometry": placed})
        ring = list(placed.exterior.coords)
        if Polygon(ring).exterior.is_ccw:
            ring = ring[::-1]
        start = max(range(len(ring) - 1), key=lambda i: ring[i][1])
        ring = ring[start:-1] + ring[:start]
        for i in range(len(ring)):
            (x0, y0), (x1, y1) = ring[i], ring[(i + 1) % len(ring)]
            grid = math.degrees(math.atan2(x1 - x0, y1 - y0)) % 360
            lon, lat = proj((x0 + x1) / 2, (y0 + y1) / 2, inverse=True)
            conv = proj.get_factors(lon, lat).meridian_convergence
            edge_rows.append({"survey_no": survey, "poly_id": props.get("poly_id"),
                              "plot_no": props.get("plot_no"), "edge_no": i + 1,
                              "length_m": round(math.dist((x0, y0), (x1, y1)), 3),
                              "bearing_grid_deg": round(grid, 2),
                              "bearing_true_deg": round((grid + conv) % 360, 2),
                              "grid_convergence_deg": round(conv, 3),
                              "geometry": LineString([(x0, y0), (x1, y1)])})
    target = (out_dir / ("%s_parcels_modified.gpkg" % survey)) if out_dir else paths.output_path(village, survey)

    def writer(path):
        gpd.GeoDataFrame(parcels, geometry="geometry", crs="EPSG:32644").to_file(path, layer="parcels", driver="GPKG")
        gpd.GeoDataFrame(edge_rows, geometry="geometry", crs="EPSG:32644").to_file(path, layer="edges", driver="GPKG")

    if target.exists():
        if not files.safe_write(target, writer):
            raise PermissionError("cannot replace %s (open in QGIS?)" % target)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        writer(target)
    return target


def run(village, do_raster=False, do_topology=False, do_review=True, project=None):
    """Steps 1-13 of spec section 5. Returns a summary dict; details go to georef_status.csv."""
    run_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    tool_fps = review.tool_written_fingerprints(village)
    anchor_map = anchors.load_anchors(village, tool_fps)
    conflicts = anchors.anchor_conflicts(village, anchor_map)
    superseded = {s: paths.manual_files(village, s)[1:] for s in anchor_map}
    anchors.write_anchors_csv(village, anchor_map, superseded, conflicts)

    todo = [s for s in paths.surveys_with_sheets(village) if s not in anchor_map]
    placed = {s: (a.theta, a.t) for s, a in anchor_map.items()}
    rows = []

    if do_raster or raster.pinned(village) is None:
        geoms = [anchors.placed_geometry(village, a) for a in anchor_map.values()]
        if geoms:
            b = unary_union(geoms).bounds
            raster.export(village, (b[0] - 100, b[1] - 100, b[2] + 100, b[3] + 100))

    index = None
    if raster.pinned(village) is not None:
        index = edgemod.SegmentIndex(edgemod.detect(raster.pinned(village)["path"], min_len_m=3.0))

    for survey in todo:
        row = {"survey": survey, "run": stamp, "method": "puvi-only", "ambiguous": True,
               "notes": "", "n_neighbours": 0, "file": paths.output_path(village, survey).name}
        hyps, win = window.start_poses(village, survey)
        pose = hyps[0] if hyps else None
        if index is not None and hyps:
            pts, brg = sheets.sample_outline(sheets.outline(sheets.load_sheet(village, survey)), 1.0)
            res = align.search(pts, brg, index, hyps)
            row.update({"share": round(res.share, 3), "margin": round(res.margin, 3),
                        "observable": res.observable, "ambiguous": res.ambiguous, "notes": res.note})
            if not res.ambiguous and res.observable:
                pose = (res.theta, res.t)
                row["method"] = "image"
        if pose is None:
            rows.append(dict(row, colour="red", confidence=0, status="waiting",
                             notes=(row["notes"] + " | no Puvi polygon and no placed neighbour").strip(" |")))
            continue
        row["heading_deg"] = round(pose[0], 3)
        placed[survey] = pose
        rows.append(row)

    # block adjustment: anchors fixed, everything else free
    free = {s: placed[s] for s in todo if s in placed}
    pair_obs, line_obs, gcp_obs, priors = [], [], [], []
    for s in free:
        priors.append(fit.PosePrior(s, free[s][0], tuple(free[s][1]), SIGMA_START, 10.0))
        v_s = sheets.outline(sheets.load_sheet(village, s))
        for other, (th_o, t_o) in placed.items():
            if other == s or neighbours.are_neighbours(village, s, other) is False:
                continue
            v_o = [tuple(p) for p in fit.transform_points(
                sheets.outline(sheets.load_sheet(village, other)), th_o, t_o)]
            v_s_ground = [tuple(p) for p in fit.transform_points(v_s, free[s][0], free[s][1])]
            for chain in match.common_chains(v_s_ground, v_o)[:1]:
                sigma = SIGMA_AUTO + (anchor_map[other].rms if other in anchor_map else 0.0)
                pair_obs += match.chain_observations(chain, s, other, sigma)
    adjusted = fit.block_adjust(free, {s: (a.theta, a.t) for s, a in anchor_map.items()},
                                pair_obs, line_obs, gcp_obs, priors)
    residuals = fit.residuals_by_pair(
        dict({s: {"theta": v["theta"], "t": v["t"]} for s, v in adjusted.items()},
             **{s: {"theta": a.theta, "t": a.t} for s, a in anchor_map.items()}), pair_obs)

    written = []
    for row in rows:
        s = row["survey"]
        if s not in adjusted:
            continue
        theta, t = adjusted[s]["theta"], adjusted[s]["t"]
        mine = [v for k, v in residuals.items() if s in k]
        row["boundary_rms_m"] = round(min(v["rms"] for v in mine), 3) if mine else None
        row["n_neighbours"] = len({k[0] if k[1] == s else k[1] for k in residuals if s in k})
        row["neighbours"] = ",".join(sorted({k[0] if k[1] == s else k[1] for k in residuals if s in k},
                                            key=paths.survey_sort_key))
        row["anchor_rms_m"] = max([anchor_map[n].rms for n in row["neighbours"].split(",")
                                   if n in anchor_map] or [0.0])
        row["shift_m"] = round(adjusted[s]["shift_m"], 2)
        row["sigma_pos_m"] = round(adjusted[s]["sigma_pos_m"], 3)
        row["sigma_head_deg"] = round(adjusted[s]["sigma_head_deg"], 3)
        if row["boundary_rms_m"] is not None and row["method"] != "image":
            row["method"] = "anchors"
        body = unary_union([fit.apply_pose(g, theta, t) for _p, g in sheets.load_sheet(village, s)])
        row["puvi_reference_m"] = round(window.puvi_distance(village, s, body) or 0.0, 1)
        trusted, ratio = window.puvi_trusted(village, s, sheets.dissolve(sheets.load_sheet(village, s)).area)
        row["puvi_trusted"], row["puvi_ratio"] = trusted, round(ratio, 2) if ratio == ratio else ""
        row["colour"] = colour_of(row)
        row["confidence"] = {"green": 90, "amber": 60, "red": 30}[row["colour"]]
        row["status"] = "placed"
        target = write_parcels(village, s, theta, t, row)
        row["fp_placed"] = anchors.fingerprint(target)
        if index is not None:
            gcp.write(village, s, gcp.corner_gcps(village, s, theta, t, index), _crs_wkt())
        written.append(s)

    if do_topology and written:
        geoms = {s: unary_union([fit.apply_pose(g, adjusted[s]["theta"], adjusted[s]["t"])
                                 for _p, g in sheets.load_sheet(village, s)]) for s in written}
        geoms.update({s: anchors.placed_geometry(village, a) for s, a in anchor_map.items()})
        rail = {s for s in geoms if s in _rail_parcels(village)}
        fixed, report = topology.fix(geoms, movable=set(written), rail=rail)
        topology.report_rail_conflicts(village, report)

    for row in rows:
        if row.get("fp_placed"):
            row["fp_final"] = row["fp_placed"]
    review.write_status(village, rows)
    if do_review:
        review.build_group(village, rows, project)
    return {"village": village, "anchors": len(anchor_map), "placed": len(written),
            "waiting": sum(1 for r in rows if r.get("status") == "waiting"),
            "conflicts": len(conflicts), "run": run_id}


def _crs_wkt():
    from pyproj import CRS
    return CRS.from_epsg(32644).to_wkt()


def _rail_parcels(village):
    """Surveys whose buffer polygon contains the traced rail centreline (reference knowledge)."""
    try:
        buf = gpd.read_file(paths.BUFFER_GPKG, layer=paths.BUFFER_LAYER).to_crs(32644)
        rail = unary_union(list(gpd.read_file(paths.RAIL_GPKG, layer=paths.RAIL_LAYER).to_crs(32644).geometry))
    except Exception:
        return set()
    sub = buf[buf["village_code"].astype(str).str.strip() == village]
    return {str(r.survey_no).strip() for r in sub.itertuples()
            if r.geometry.intersection(rail).length >= 50.0}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_engine.py -v`
Expected: PASS, 6 tests. The last test is the important one: on Kizhikaranai every survey is an anchor, so the engine must place nothing and leave every manual file byte-identical.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/engine.py tests/test_engine.py
git commit -m "feat(autogeoref): per-village engine wiring anchors, imagery, adjustment and outputs"
```

---

## Task 16: Command line

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\georef_village.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_cli.py`

**Interfaces:**
- Consumes: `engine.run`, `gcp.refit`, `review`.
- Produces: `georef_village.main(argv=None) -> int` and the console behaviour of `python -m autogeoref.georef_village <village> [--raster] [--refit] [--reset-gcp SURVEY] [--topology] [--review] [--project QGZ] [--tracker]`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_cli.py`:

```python
import pytest
from autogeoref import georef_village


def test_parser_defaults():
    args = georef_village.parse(["35_04_077"])
    assert args.village == "35_04_077"
    assert args.raster is False and args.topology is False and args.review is True
    assert args.tracker is False, "the tracker is not written unless asked (user rule)"


def test_parser_flags():
    args = georef_village.parse(["35_04_052", "--raster", "--topology", "--no-review", "--tracker"])
    assert args.raster and args.topology and args.review is False and args.tracker


def test_refit_requires_a_survey_list():
    args = georef_village.parse(["35_04_077", "--refit", "46B,47A"])
    assert args.refit == "46B,47A"


def test_main_returns_nonzero_for_an_unknown_village(capsys):
    assert georef_village.main(["99_99_999"]) != 0
    assert "no sheets" in capsys.readouterr().out.lower()
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_cli.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.georef_village'`

- [ ] **Step 3: Write `autogeoref/georef_village.py`**

```python
"""python -m autogeoref.georef_village 35_04_077 [--raster] [--topology] [--refit 46B,47A] ..."""
import argparse
import sys

from . import engine, gcp, paths, review


def parse(argv=None):
    ap = argparse.ArgumentParser(prog="georef_village", description="Auto-georeference one village.")
    ap.add_argument("village", help="village code, e.g. 35_04_077")
    ap.add_argument("--raster", action="store_true", help="re-export and pin the satellite GeoTIFF")
    ap.add_argument("--topology", action="store_true", help="run the topology fix on the placed parcels")
    ap.add_argument("--refit", metavar="SURVEYS", help="comma list: re-fit from the team's edited .points")
    ap.add_argument("--reset-gcp", metavar="SURVEY", dest="reset_gcp", help="drop the team GCP rows for one survey")
    ap.add_argument("--review", dest="review", action="store_true", default=True)
    ap.add_argument("--no-review", dest="review", action="store_false", help="skip the QGIS review group")
    ap.add_argument("--project", help="path to the .qgz the review group belongs in")
    ap.add_argument("--tracker", action="store_true",
                    help="also add the Auto colour/note/run columns to the tracker workbook")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse(argv)
    if not paths.vector_dir(args.village).exists() or not paths.surveys_with_sheets(args.village):
        print("no sheets found for village %s under %s" % (args.village, paths.vector_dir(args.village)))
        return 2
    if args.reset_gcp:
        p = paths.gcp_csv(args.village, args.reset_gcp)
        if p.exists():
            p.unlink()
        print("team GCP record cleared for", args.reset_gcp)
    if args.refit:
        for survey in [s.strip() for s in args.refit.split(",") if s.strip()]:
            got = gcp.refit(args.village, survey)
            if got is None:
                print("%-6s no team GCPs newer than the tool's record (or fewer than two)" % survey)
                continue
            theta, t, rms, n = got
            row = {"survey": survey, "method": "team-gcp", "colour": "green", "confidence": 95,
                   "notes": "refit from %d team GCPs, rms %.2f m" % (n, rms)}
            engine.write_parcels(args.village, survey, theta, t, row)
            print("%-6s refitted from %d team GCPs: heading %.2f, rms %.2f m" % (survey, n, theta, rms))
        return 0
    summary = engine.run(args.village, do_raster=args.raster, do_topology=args.topology,
                         do_review=args.review, project=args.project)
    rows = review.read_status(args.village)
    print("%s: %d anchors frozen, %d placed, %d waiting, %d anchor conflicts"
          % (summary["village"], summary["anchors"], summary["placed"], summary["waiting"],
             summary["conflicts"]))
    for colour in ("green", "amber", "red"):
        names = [r["survey"] for r in rows if r.get("colour") == colour]
        if names:
            print("  %-5s %s" % (colour, ", ".join(names)))
    if args.tracker:
        print("  tracker:", review.update_tracker(
            [dict(r, village_code=args.village, run=r.get("run", "")) for r in rows]))
    review.write_worklist({args.village: rows})
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_cli.py -v`
Expected: PASS, 4 tests.

- [ ] **Step 5: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/georef_village.py tests/test_cli.py
git commit -m "feat(autogeoref): command line for the village run, refit and reset"
```

---

## Task 17: Threshold calibration and the leave-one-out harness

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\autogeoref\evaluate.py`
- Create: `D:\code\FMB_to_GeoJSON\tests\test_evaluate.py`
- Modify: `D:\code\FMB_to_GeoJSON\autogeoref\engine.py` (write the calibrated `GREEN_SHARE` / `AMBER_SHARE` back)
- Modify: `D:\code\FMB_to_GeoJSON\docs\superpowers\specs\2026-09-18-fmb-georef-gcp-sop-design.md` (record the calibration in §3)

**Interfaces:**
- Consumes: `engine`, `align`, `anchors`, `edges`, `raster`.
- Produces:
  - `evaluate.truth_pose(village, survey) -> dict` — `{"theta", "t", "rms", "heading_spread", "source"}`; the rigid fit of the sheet to the manual file, plus a jackknife heading spread.
  - `evaluate.tolerances(truth) -> (float, float)` — `max(1.5, rms)` metres and `max(1.5, heading_spread)` degrees.
  - `evaluate.leave_one_out(village, surveys=None, mode="full") -> list[dict]` — hides one manual parcel at a time on a **copy** of the village folder under `_logs\loo_<village>_<date>\`, runs the engine, compares. `mode="gcp_only"` sets the anchor weight to zero so the image path is what is being measured.
  - `evaluate.calibrate(results) -> dict` — `{"green_share", "amber_share", "table"}` from the score at the truth pose against the best wrong pose.

- [ ] **Step 1: Write the failing test**

Create `tests/test_evaluate.py`:

```python
import pytest
from autogeoref import evaluate

VILLAGE = "35_04_077"


def test_truth_pose_reports_its_own_uncertainty():
    t = evaluate.truth_pose(VILLAGE, "47B")
    assert t["source"] in ("points", "geometry")
    assert t["rms"] > 3.0, "47B is the worst hand placement on this village"
    assert t["heading_spread"] >= 0.0


def test_tolerances_never_tighter_than_the_truth_itself():
    m, deg = evaluate.tolerances({"rms": 5.85, "heading_spread": 1.6})
    assert m == pytest.approx(5.85) and deg == pytest.approx(1.6)
    m2, deg2 = evaluate.tolerances({"rms": 0.19, "heading_spread": 0.3})
    assert m2 == pytest.approx(1.5) and deg2 == pytest.approx(1.5)


def test_calibrate_separates_truth_from_the_best_wrong_pose():
    rows = [{"survey": "48A", "share_truth": 0.72, "share_best_wrong": 0.44},
            {"survey": "47B", "share_truth": 0.48, "share_best_wrong": 0.30},
            {"survey": "46B", "share_truth": 0.12, "share_best_wrong": 0.46}]
    cal = evaluate.calibrate(rows)
    assert 0.0 < cal["amber_share"] < cal["green_share"] <= 1.0
    assert cal["green_share"] > 0.46, "a threshold below the best wrong pose would pass 46B wrongly"


@pytest.mark.slow
def test_leave_one_out_runs_on_a_copy_and_never_touches_the_project():
    from autogeoref import paths
    before = {p.name: p.stat().st_mtime for p in paths.vector_dir(VILLAGE).glob("*_parcels_modified*.gpkg")}
    out = evaluate.leave_one_out(VILLAGE, surveys=["43A"], mode="full")
    after = {p.name: p.stat().st_mtime for p in paths.vector_dir(VILLAGE).glob("*_parcels_modified*.gpkg")}
    assert before == after, "the harness must work on a copy"
    assert out and set(out[0]) >= {"survey", "colour", "centroid_error_m", "heading_error_deg", "passed"}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_evaluate.py -v -m "not slow"`
Expected: FAIL with `ModuleNotFoundError: No module named 'autogeoref.evaluate'`

- [ ] **Step 3: Write `autogeoref/evaluate.py`**

```python
"""Leave-one-out evaluation against the team's own placements, and threshold calibration.

The team's files are affine and vertex-edited, so the truth pose is the rigid fit of the sheet to
the manual file and its tolerance is that fit's own rms - never tighter than 1.5 m / 1.5 degrees.
Parcels whose rigid rms exceeds 3 m (47B, 48B, 42A, 43B, 169 here; 55 and 56 on Thailavaram) are
reported but not counted.
"""
import datetime
import shutil

import numpy as np

from . import anchors, engine, fit, paths, sheets

REPORT_ONLY_RMS = 3.0


def truth_pose(village, survey):
    pose = anchors.pose_from_points(village, survey)
    source = "points"
    if pose is None or pose[2] > REPORT_ONLY_RMS:
        hand = paths.manual_files(village, survey)
        if not hand:
            return None
        pose = anchors.pose_from_geometry(village, survey, hand[0])
        source = "geometry"
    theta, t, rms, scale, aniso, n = pose
    spread = _jackknife_heading(village, survey)
    return {"theta": theta, "t": np.asarray(t, float), "rms": rms, "scale": scale,
            "anisotropy": aniso, "n": n, "heading_spread": spread, "source": source}


def _jackknife_heading(village, survey):
    """How much the recovered heading moves when one control point is dropped."""
    p = paths.points_path(village, survey)
    if not p.exists():
        return 0.0
    P = anchors.read_points(p)
    if len(P) < 4:
        return 0.0
    out = []
    for k in range(len(P)):
        keep = np.delete(P, k, axis=0)
        theta, _t, _rms, _mx = fit.rigid_fit(keep[:, 2:4], keep[:, 0:2])
        out.append(theta)
    return float(max(out) - min(out))


def tolerances(truth):
    return max(1.5, float(truth["rms"])), max(1.5, float(truth.get("heading_spread", 0.0)))


def leave_one_out(village, surveys=None, mode="full"):
    """Hide one manual parcel at a time on a copy of the village folder and score the placement."""
    stamp = datetime.datetime.now().strftime("%Y%m%d")
    work_root = paths.logs_dir() / ("loo_%s_%s" % (village, stamp))
    src = paths.vector_dir(village)
    results = []
    targets = surveys or sorted({p.name.split("_")[0] for p in src.glob("*_parcels_modified*.gpkg")},
                                key=paths.survey_sort_key)
    real_project = paths.PROJECT
    for survey in targets:
        truth = truth_pose(village, survey)
        if truth is None:
            continue
        work = work_root / survey
        if work.exists():
            shutil.rmtree(work)
        (work / "FMB_Vector" / village).mkdir(parents=True)
        for f in src.iterdir():
            if f.is_file() and not f.name.startswith(survey + "_parcels_modified") \
                    and f.name != ("%s_parcels.geojson.points" % survey):
                shutil.copy2(f, work / "FMB_Vector" / village / f.name)
        for extra in ("FMB_Georef", "FMB_Sketches"):
            s = real_project / extra / village
            if s.exists():
                shutil.copytree(s, work / extra / village, dirs_exist_ok=True)
        try:
            paths.PROJECT = work
            if mode == "gcp_only":
                engine.SIGMA_AUTO = 1e6           # anchors contribute nothing
            summary = engine.run(village, do_raster=False, do_topology=False, do_review=False)
            from . import review
            row = next((r for r in review.read_status(village) if r["survey"] == survey), None)
        finally:
            paths.PROJECT = real_project
            engine.SIGMA_AUTO = 0.30
        if row is None:
            results.append({"survey": survey, "colour": "missing", "passed": False, "mode": mode})
            continue
        got_theta = float(row.get("heading_deg") or 0.0)
        placed = sheets.dissolve([(None, fit.apply_pose(g, got_theta, np.array([0.0, 0.0])))
                                  for _p, g in sheets.load_sheet(village, survey)])
        tol_m, tol_deg = tolerances(truth)
        err_deg = abs(((got_theta - truth["theta"]) + 180) % 360 - 180)
        err_m = float(row.get("shift_m") or 0.0)
        results.append({"survey": survey, "mode": mode, "colour": row.get("colour"),
                        "share": row.get("share"), "margin": row.get("margin"),
                        "centroid_error_m": round(err_m, 2), "heading_error_deg": round(err_deg, 2),
                        "tolerance_m": round(tol_m, 2), "tolerance_deg": round(tol_deg, 2),
                        "truth_rms_m": round(truth["rms"], 2),
                        "counted": truth["rms"] <= REPORT_ONLY_RMS,
                        "passed": err_m <= tol_m and err_deg <= tol_deg})
    return results


def calibrate(rows):
    """Green must sit above every best-wrong score; amber above the median wrong score."""
    wrong = [r["share_best_wrong"] for r in rows if r.get("share_best_wrong") is not None]
    truth = [r["share_truth"] for r in rows if r.get("share_truth") is not None]
    green = max(wrong) + 0.05 if wrong else 0.60
    amber = float(np.median(wrong)) if wrong else 0.35
    return {"green_share": round(min(green, 0.95), 3), "amber_share": round(min(amber, green - 0.05), 3),
            "table": sorted(rows, key=lambda r: -(r.get("share_truth") or 0))}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_evaluate.py -v -m "not slow"`
Expected: PASS, 3 tests. Then run the slow one: `python -m pytest tests/test_evaluate.py -v -m slow` (add `markers = slow` to `pytest.ini`).

- [ ] **Step 5: Run the real calibration and write the numbers back**

```bash
cd /d/code/FMB_to_GeoJSON
python - <<'PY'
from autogeoref import evaluate
import json
rows = evaluate.leave_one_out("35_04_077", mode="full") + evaluate.leave_one_out("35_04_052", mode="full")
print(json.dumps(rows, indent=1, default=str))
print("green passes:", sum(1 for r in rows if r["colour"] == "green" and r["passed"]),
      "of", sum(1 for r in rows if r["colour"] == "green"))
PY
```

Acceptance before the plan is done: every parcel reported green passes its own tolerance, and no flipped or slid parcel is green. Put the resulting `green_share` and `amber_share` into `engine.GREEN_SHARE` / `engine.AMBER_SHARE`, and append the per-parcel table to §3 of the spec.

- [ ] **Step 6: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add autogeoref/evaluate.py autogeoref/engine.py tests/test_evaluate.py pytest.ini docs/superpowers/specs/2026-09-18-fmb-georef-gcp-sop-design.md
git commit -m "feat(autogeoref): leave-one-out harness and calibrated colour thresholds"
```

---

## Task 18: Full-village acceptance run

**Files:**
- Create: `D:\code\FMB_to_GeoJSON\tests\test_acceptance.py`

**Interfaces:**
- Consumes: `engine.run`, `review.read_status`, `anchors.fingerprint`.
- Produces: nothing new; this is the gate that the whole pipeline behaves on real data.

- [ ] **Step 1: Write the acceptance test**

Create `tests/test_acceptance.py`:

```python
import pytest
from autogeoref import anchors, engine, paths, review

pytestmark = pytest.mark.slow
VILLAGE = "35_04_077"


def test_a_full_run_leaves_every_manual_file_and_points_file_byte_identical():
    manual = {p: p.read_bytes() for p in paths.vector_dir(VILLAGE).glob("*_parcels_modified*.gpkg")}
    points = {p: p.read_bytes() for p in paths.vector_dir(VILLAGE).glob("*.points")}
    engine.run(VILLAGE, do_raster=False, do_topology=True, do_review=False)
    for p, data in manual.items():
        assert p.read_bytes() == data, "%s was modified" % p.name
    for p, data in points.items():
        assert p.read_bytes() == data, "%s was modified" % p.name


def test_every_buffer_parcel_is_classified():
    engine.run(VILLAGE, do_raster=False, do_topology=False, do_review=False)
    rows = review.read_status(VILLAGE)
    got = {r["survey"] for r in rows} | set(anchors.load_anchors(VILLAGE, set()))
    assert set(paths.surveys_with_sheets(VILLAGE)) <= got
    assert all(r.get("colour") in ("green", "amber", "red") for r in rows)


def test_railway_strips_are_not_flipped_when_they_are_placed():
    """169, 170 and 171 are near-symmetric; only the transcription can orient them."""
    rows = {r["survey"]: r for r in review.read_status(VILLAGE)}
    for s in ("169", "170", "171"):
        if s in rows and rows[s].get("status") == "placed":
            assert rows[s]["colour"] != "green" or rows[s].get("observable") == "True"
```

- [ ] **Step 2: Run it**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest tests/test_acceptance.py -v -m slow`
Expected: PASS, 3 tests. Kizhikaranai is fully hand-placed, so the run should place nothing; the value of the test is that it proves the tool cannot damage the team's work.

- [ ] **Step 3: Run the whole suite**

Run: `cd /d/code/FMB_to_GeoJSON && python -m pytest -v`
Expected: PASS, about 70 tests.

- [ ] **Step 4: Commit**

```bash
cd /d/code/FMB_to_GeoJSON
git add tests/test_acceptance.py
git commit -m "test(autogeoref): full-village acceptance run on Kizhikaranai"
```

---

## Self-review notes for the implementer

Three places where the plan knowingly leaves a judgement to you, with the evidence to decide:

1. **`engine.run` places at most one chain per neighbour pair** (`common_chains(...)[:1]`). If the leave-one-out table in Task 17 shows parcels landing with a good boundary rms but a wrong along-track position, admit more chains and apply the trusted-first rule from `scratch_2026-09-17/block_adjust.py`: multi-vertex chains first, demote any that still exceed 1.0 m, then admit single-edge chains one at a time by trial solve.
2. **`GREEN_SHARE` starts at 0.60**, which is the spike's binary-metric number. Task 17 replaces it with a calibrated value; do not ship the placeholder.
3. **Image evidence is currently only used to choose a pose, not as `GcpObs` in the adjustment.** The hooks (`gcp_obs`, `SIGMA_IMAGE`) exist. Add them once Task 17 shows the image poses are trustworthy on this corridor, and re-run the acceptance test.
4. **Spec section 5 step 11 separates a team *pose* from a frozen *anchor***, and Task 3 does not implement it yet: a parcel the team merely dragged (its sheet still fits at scale 1 plus or minus 0.2 % with vertex rms below 0.05 m) can stay in the block adjustment as a strong prior, while one whose dimensions changed must be frozen. In `anchors.load_anchors`, after the pose is recovered, set `kind = "team-pose" if (abs(scale - 1.0) <= 0.002 and rms < 0.05) else "anchor"`, carry `kind` on the `Anchor` dataclass, and in `engine.run` put `team-pose` surveys into `free` with a `PosePrior` of sigma 0.30 m and 0.5 degrees instead of into `fixed`. Add a test that a hand file which is an exact rigid copy of its sheet comes back as `team-pose`, and one scaled by 1.05 comes back as `anchor`.
