"""Manual (hand-georeferenced) parcels: choose one file per survey, fingerprint it, and recover the
rigid pose it implies. The files themselves are never modified."""
import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

import fiona
import numpy as np
from shapely.geometry import Point, shape
from shapely.ops import unary_union

from . import fit, paths, sheets, visible

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
    """sha256 of the layer's geometries, read through OGR so a WAL sidecar is included."""
    h = hashlib.sha256()
    with fiona.open(str(gpkg_path)) as src:
        for feat in src:
            geom = feat["geometry"]
            if geom is None:
                h.update(b"None")
                continue
            h.update(shape(geom).wkb)
    return h.hexdigest()


def _fit_pose(P_sheet, Q_map):
    theta, t, rms, _mx = fit.rigid_fit(P_sheet, Q_map)
    _th2, _t2, _rms2, scale, aniso = fit.similarity_fit(P_sheet, Q_map)
    return theta, t, rms, scale, aniso, len(P_sheet)


def pose_from_points(village, survey):
    """Rigid pose implied by the team's own control points, or None when there is no usable file."""
    p = paths.points_path(village, survey)
    if not p.exists():
        return None
    P = read_points(p)
    if len(P) < 2:
        return None
    return _fit_pose(P[:, 2:4], P[:, 0:2])


VERTEX_REACH = 4.0          # metres: how far a sheet vertex may look for its hand-edited twin


def pose_from_geometry(village, survey, gpkg_path):
    """Rigid pose implied by the hand-placed geometry.

    The team edits vertices, so matching by vertex order fails on most sheets. This is the method
    proven on the 2026-09-17 Thailavaram file: align by poly_id centroids (or, for a single-polygon
    sheet, by a rotation search), then refine on nearest vertices within VERTEX_REACH.
    """
    import geopandas as gpd
    from shapely import affinity

    hand = visible.read_hand(gpkg_path)
    if hand.crs is not None and hand.crs.to_epsg() != 32644:
        hand = hand.to_crs(32644)
    sheet = {int(p.get("poly_id", i)): g for i, (p, g) in enumerate(sheets.load_sheet(village, survey))}
    by_id = {}
    for row in hand.itertuples():
        if row.geometry is None or row.geometry.is_empty:
            continue
        try:
            by_id[int(row.poly_id)] = row.geometry
        except (TypeError, ValueError, AttributeError):
            pass
    if not by_id:                                     # no poly_id: pair the largest parts by size
        A = sorted(sheet.items(), key=lambda kv: -kv[1].area)
        B = sorted([g for g in hand.geometry if g is not None and not g.is_empty], key=lambda g: -g.area)
        by_id = {k: b for (k, _a), b in zip(A, B)}
    common = [k for k in sheet if k in by_id]
    if not common:
        return float("nan"), np.zeros(2), float("nan"), float("nan"), float("nan"), 0

    if len(common) >= 2:                              # stage 1a: centroids give the pose directly
        theta, t, _rms, _mx = fit.rigid_fit(
            np.array([sheet[k].centroid.coords[0] for k in common], float),
            np.array([by_id[k].centroid.coords[0] for k in common], float))
    else:                                             # stage 1b: one polygon -> search the rotation
        k = common[0]
        src, dst = sheet[k], by_id[k]
        best = None
        for deg in np.arange(0.0, 360.0, 0.5):
            g = affinity.rotate(src, deg, origin="centroid")
            g = affinity.translate(g, dst.centroid.x - g.centroid.x, dst.centroid.y - g.centroid.y)
            d = g.symmetric_difference(dst).area
            if best is None or d < best[0]:
                best = (d, deg)
        theta = best[1]
        c = np.array(src.centroid.coords[0])
        t = np.array(dst.centroid.coords[0]) - fit.transform_points([c], theta, np.zeros(2))[0]

    # stage 2: refine on nearest vertices, starting wide because a stretched hand file (48B is
    # scaled 0.90) puts its vertices ten metres or more from any rigid placement of the sheet
    size = float(np.sqrt(max(sum(g.area for g in sheet.values()), 1.0)))
    best_P, best_Q = [], []
    for reach in (max(VERTEX_REACH, 0.15 * size), max(VERTEX_REACH, 0.07 * size), VERTEX_REACH):
        P, Q = [], []
        for k in common:
            other = by_id[k]
            if other.geom_type in ("MultiPolygon", "GeometryCollection"):
                parts = [q for q in other.geoms if q.geom_type == "Polygon" and q.area > 0]
                if not parts:
                    continue            # a degenerate sliver anchors nothing
                other = max(parts, key=lambda q: q.area)
            if other.geom_type != "Polygon":
                continue
            ring = other.exterior
            mv = np.array(list(ring.coords)[:-1], float)
            for v in list(sheet[k].exterior.coords)[:-1]:
                w = fit.transform_points([v], theta, t)[0]
                d = np.linalg.norm(mv - w, axis=1)
                j = int(d.argmin())
                if d[j] <= reach:
                    P.append(v)
                    Q.append(mv[j])
        if len(P) < 3:
            continue
        best_P, best_Q = P, Q
        theta, t, _rms, _mx = fit.rigid_fit(np.array(P, float), np.array(Q, float))
    if len(best_P) >= 3:
        return _fit_pose(np.array(best_P, float), np.array(best_Q, float))

    # no vertex correspondence survived: keep the stage 1 pose, describe it from areas and offsets
    sheet_area = sum(g.area for g in sheet.values())
    hand_area = sum(by_id[k].area for k in common)
    scale = float(np.sqrt(hand_area / sheet_area)) if sheet_area > 0 else float("nan")
    res = []
    for k in common:
        other = by_id[k]
        if other.geom_type in ("MultiPolygon", "GeometryCollection"):
            parts = [q for q in other.geoms if q.geom_type == "Polygon" and q.area > 0]
            if not parts:
                continue            # a degenerate sliver describes nothing
            other = max(parts, key=lambda q: q.area)
        if other.geom_type != "Polygon":
            continue
        boundary = other.exterior
        for v in list(sheet[k].exterior.coords)[:-1]:
            w = fit.transform_points([v], theta, t)[0]
            res.append(boundary.distance(Point(w)))
    rms = float(np.sqrt(np.mean(np.square(res)))) if res else float("nan")
    return theta, np.asarray(t, float), rms, scale, 1.0, len(common)


def load_anchors(village, tool_written):
    """One Anchor per survey that has a manual file the tool did not write."""
    from . import visible
    out = {}
    for survey in paths.surveys_with_sheets(village):
        for f in paths.manual_files(village, survey):       # the file his project shows first
            fp = fingerprint(f)
            if fp in tool_written and not visible.is_shown(f):
                continue                                    # an approved parcel is his, whoever drew it
            pose = pose_from_points(village, survey)
            source = "points"
            if pose is None or pose[2] > POSE_RMS_LIMIT:
                pose = pose_from_geometry(village, survey, f)
                source = "geometry"
            theta, t, rms, scale, aniso, n = pose
            # the scale that matters to the team is how much the hand placement stretched the
            # parcel as a whole, so measure it from the areas rather than from a handful of points
            area_scale = _area_scale(village, survey, f)
            if area_scale == area_scale:
                scale = area_scale
            out[survey] = Anchor(survey, f, theta, np.asarray(t, float), rms, scale, aniso, source, fp, n)
            break
    return out


def _area_scale(village, survey, gpkg_path):
    """sqrt(hand area / sheet area): 1.00 means the hand placement kept the FMB size."""
    import geopandas as gpd
    try:
        hand = visible.read_hand(gpkg_path)
    except Exception:
        return float("nan")
    if hand.crs is not None and hand.crs.to_epsg() != 32644:
        hand = hand.to_crs(32644)
    hand_area = sum(g.area for g in hand.geometry if g is not None and not g.is_empty)
    sheet_area = sum(g.area for _p, g in sheets.load_sheet(village, survey))
    return float(np.sqrt(hand_area / sheet_area)) if sheet_area > 0 else float("nan")


def placed_geometry(village, anchor):
    """The rigid sheet at the anchor's recovered pose (NOT the hand geometry)."""
    sheet = sheets.load_sheet(village, anchor.survey)
    return unary_union([fit.apply_pose(g, anchor.theta, anchor.t).buffer(0) for _p, g in sheet])


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
            ring = ga.exterior if ga.geom_type == "Polygon" else ga.boundary   # a hand file may dissolve to several parts
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
