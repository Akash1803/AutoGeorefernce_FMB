"""Write a shift-corrected copy of the Puvi survey polygons for one or more villages.

    python -m autogeoref.shift_puvi 35_04_074 35_04_052 --buffer-only
    python -m autogeoref.shift_puvi --stretch          # Thailavaram to Thirukatchur

Control is the team's own placements: every `<survey>_parcels_modified.gpkg` in the village folder,
plus any file named with `--control`, which is how the Thailavaram merged layer comes in. The
displacement between a control parcel and its Puvi twin is measured at the centroid, the field is
fitted by `shiftfit`, and every distinct vertex in the village is moved through it once, so two
parcels that share a boundary keep sharing it. Nothing existing is edited, and a village with no
control is written out unchanged and marked so.
"""
import argparse
import datetime
import itertools
import logging
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity, set_precision
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union
from shapely.validation import make_valid

from . import engine, paths, shiftfit, visible

log = logging.getLogger(__name__)

MAX_CONTROL_M = 200.0   # a control parcel further than this from its Puvi twin is a key clash
STRETCH = ["35_04_052", "35_04_054", "35_04_056", "35_04_076", "35_04_077", "35_04_074"]

# every revenue village the 30 m rail buffer passes through that has a Puvi vector, ordered along
# the line from the Tambaram end. Tambaram itself (35_05_010) has no Puvi vector in the archive.
CORRIDOR = ["03_14_194", "35_05_136", "35_05_134", "35_15_002", "35_15_003", "35_15_004",
            "35_15_006", "35_15_005", "35_04_052", "35_04_054", "35_04_056", "35_04_076",
            "35_04_077", "35_04_074", "35_04_073", "35_04_085", "35_04_086", "35_04_071",
            "35_04_221", "35_04_224", "35_04_225", "35_04_226"]
UTM = 32644


NODE_GRID_M = 0.01   # two vertices this close are the same drawn point, and move together


def _nodes(geoms, grid=NODE_GRID_M):
    """Every distinct vertex of the parcels, on a centimetre grid.

    Puvi's neighbours draw a shared boundary with coordinates that agree to a few millimetres but
    are not the same numbers. Warping each polygon on its own then leaves hairline slivers along
    every shared line (594 of them, averaging 0.15 m2, on the first warp). Snapping to a grid and
    moving each distinct node once makes a shared line one line again.
    """
    out = {}
    for g in geoms:
        for part in _polygons_of(g):
            for ring in [part.exterior, *part.interiors]:
                for c in ring.coords:                       # some Puvi polygons carry a Z value
                    x, y = float(c[0]), float(c[1])
                    out[(round(x / grid), round(y / grid))] = (x, y)
    return out


def _warp_table(nodes, points, disp, grid=NODE_GRID_M, **fit_kw):
    """Where each distinct node lands. One field evaluation per node, not per polygon."""
    keys = list(nodes)
    coords = np.array([nodes[k] for k in keys], float)
    if not len(coords):
        return {}
    shifts, _ = shiftfit.fit(points, disp, coords, **fit_kw)
    return {k: tuple(c + s) for k, c, s in zip(keys, coords, shifts)}


def _snap_or_keep(geom, grid=1e-9):
    """Snap to a fine grid so the parcel survives being written in degrees; keep it if that fails."""
    if geom is None or geom.is_empty:
        return geom
    try:
        snapped = _repair(set_precision(geom, grid))
        if snapped is not None and not snapped.is_empty:
            return snapped
    except Exception:
        pass
    return _repair(geom)


def _repair(geom):
    """Make a parcel valid again after a difference, without dropping any of it."""
    if geom is None or geom.is_empty:
        return geom
    g = geom if geom.is_valid else make_valid(geom)
    parts = _polygons_of(g)
    if parts:
        g = unary_union(parts)
    if g is not None and not g.is_empty and not g.is_valid:
        g = g.buffer(0)
    return g


def _polygons_of(geom):
    """Every polygon inside a geometry, whatever container it arrived in.

    Puvi survey 11 of Peramanur is a GeometryCollection, and the first whole-village run dropped it
    because the warp only looked at Polygon and MultiPolygon.
    """
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == "Polygon":
        return [geom]
    if geom.geom_type in ("MultiPolygon", "GeometryCollection"):
        return [x for g in geom.geoms for x in _polygons_of(g)]
    return []


def _warp_with(geom, table, grid=NODE_GRID_M, field=None):
    """Rebuild a parcel from the warped node table, so shared vertices stay shared.

    A parcel is never lost: if the rebuilt shape comes back empty or broken, the original is moved
    rigidly by the field at its centroid instead.
    """
    def ring(coords):
        out = []
        for c in coords:
            x, y = float(c[0]), float(c[1])
            out.append(table.get((round(x / grid), round(y / grid)), (x, y)))
        return out

    def poly(p):
        return Polygon(ring(p.exterior.coords), [ring(h.coords) for h in p.interiors])

    parts = _polygons_of(geom)
    out = _valid(MultiPolygon([poly(p) for p in parts])) if len(parts) > 1 else (
        _valid(poly(parts[0])) if parts else None)
    if out is not None and not out.is_empty:
        return out
    if field is not None and geom is not None and not geom.is_empty:
        d = field(np.array(geom.centroid.coords[0]))
        return _valid(affinity.translate(geom, float(d[0]), float(d[1])))
    return _valid(geom)


MEASURED_REACH_M = 150.0  # beyond this from a control parcel the correction is not evidenced

# Measured 2026-09-22 by predicting each control parcel from the others, grouped by how far the
# nearest remaining one is: under 150 m the fit works (Thirukatchur 24.0 -> 3.2 m); at 150 to 400 m
# it did nothing (40.2 -> 39.4 m). A whole village reaches far past that (Thirukatchur's median
# parcel is 759 m from the nearest placement, 35 % are over a kilometre), so those parcels get the
# village's average correction and are labelled as not evidenced rather than scored.
# The revenue village boundary is not a second source: its outline sits on Puvi's to the metre in
# three villages, so it was digitised from the same survey and carries the same error.

SEAM_REACH_M = 120.0   # how far across a village boundary to look for the neighbour's edge
SEAM_STEP_M = 15.0     # sampling step along the shared frontage


RAIL_TOL_M = 5.0          # leave a village alone if its railway land is already this close to the norm
RAIL_MAX_WIDTH_M = 45.0   # a plain strip this wide or less is railway land on its width alone
RAIL_WIDE_MAX_M = 85.0    # a wider one counts only if its own sheet draws the track through it
# Peramanur 77 is 57 m wide and 76 is 48 m: station land, wider than a running corridor. The width
# rule alone threw them out and left 77 sixteen metres off the track, which is what Akash saw on
# 2026-09-22. Their sheets draw the track (211 m and 199 m of it), so the sheet decides instead.
# Survey 43 stays out: 610 m across and no track drawn on its sheet at all.

# Why a tolerance at all: the norm is +3.0 m, but the three villages where Akash's placements give
# the truth sit at +1.7, +3.0 and +6.1 m. The norm is therefore only good to about +-2.5 m, and
# forcing every village onto exactly +3.0 m claims a precision the evidence does not have. Where
# Puvi already has the railway land within RAIL_TOL_M of the norm, Puvi is as good as anything we
# can say, and it is left alone (he spotted this on 2026-09-22: "Puvi vector is good compare to our
# results?").
RAIL_MIN_SAMPLES = 8


def _rail_line():
    import shapely.ops as _ops
    g = gpd.read_file(paths.RAIL_GPKG, layer="rail_line").to_crs(UTM).geometry.iloc[0]
    return _ops.linemerge(g) if g.geom_type == "MultiLineString" else g


def across_track(geom, line, step=8.0):
    """Signed distances from the track to a parcel's outline: negative left, positive right."""
    from shapely.geometry import Point
    out = []
    for g in _polygons_of(geom):
        c = np.array(g.exterior.coords)
        for i in range(len(c) - 1):
            a, b = np.array(c[i][:2], float), np.array(c[i + 1][:2], float)
            L = float(np.hypot(*(b - a)))
            for k in range(max(1, int(L / step))):
                pt = a + (b - a) * k / max(1, int(L / step))
                sdist = line.project(Point(pt))
                q = line.interpolate(sdist)
                t0 = line.interpolate(max(0.0, sdist - 15.0))
                t1 = line.interpolate(min(line.length, sdist + 15.0))
                d = np.array([t1.x - t0.x, t1.y - t0.y])
                n = np.linalg.norm(d)
                if n < 1e-9:
                    continue
                d /= n
                out.append(float((pt - np.array([q.x, q.y])) @ np.array([-d[1], d[0]])))
    return np.array(out)


def sheet_track(village, survey, min_span=0.5):
    """The railway line the FMB sheet itself draws inside a rail parcel, in sheet metres.

    A rail sheet draws the track as a ticked line running the length of the strip, between the two
    long boundary edges. Taking it from the sheet is better than assuming the strip is centred on
    the track: the sheet says where the track sits inside this particular parcel.

    Returns (LineString, parcel outline) in sheet coordinates, or (None, None).
    """
    from shapely.geometry import LineString
    pdir = paths.vector_dir(village)
    pf, lf = pdir / ("%s_parcels.geojson" % survey), pdir / ("%s_lines.geojson" % survey)
    if not pf.exists() or not lf.exists():
        return None, None
    try:
        poly = unary_union([_valid(g) for g in gpd.read_file(pf).geometry if g is not None])
        lines = gpd.read_file(lf)
    except Exception:
        return None, None
    if poly is None or poly.is_empty:
        return None, None
    pts = np.array(poly.exterior.coords)[:, :2] if poly.geom_type == "Polygon" else np.array(
        max(poly.geoms, key=lambda g: g.area).exterior.coords)[:, :2]
    centre = pts.mean(axis=0)
    u, sv, vt = np.linalg.svd(pts - centre)
    axis = vt[0]                       # along the strip
    normal = np.array([-axis[1], axis[0]])
    span = float(np.ptp((pts - centre) @ axis))
    best, best_score = None, None
    for g in lines.geometry:
        if g is None or g.is_empty:
            continue
        for part in (g.geoms if g.geom_type.startswith("Multi") else [g]):
            if part.geom_type != "LineString" or part.length < min_span * span:
                continue
            c = np.array(part.coords)[:, :2]
            off = (c - centre) @ normal
            # the track runs the length of the strip and stays near its middle, unlike the two
            # long boundary edges which sit at the extremes
            score = float(np.mean(np.abs(off)))
            if best_score is None or score < best_score:
                best, best_score = part, score
    if best is None:
        return None, None
    return LineString(np.array(best.coords)[:, :2]), poly


def rigid_fit(src, dst, coarse=2.0, fine=0.25):
    """Rotation and shift (no scale) putting `src` over `dst`, by best overlap."""
    from shapely import affinity
    sc = np.array(src.centroid.coords[0])
    dc = np.array(dst.centroid.coords[0])
    base = affinity.translate(src, dc[0] - sc[0], dc[1] - sc[1])
    best, best_ang = None, 0.0
    for ang in np.arange(0.0, 360.0, coarse):
        g = affinity.rotate(base, ang, origin="centroid")
        inter = g.intersection(dst).area
        score = inter / (g.area + dst.area - inter)
        if best is None or score > best:
            best, best_ang = score, ang
    for ang in np.arange(best_ang - coarse, best_ang + coarse + 1e-9, fine):
        g = affinity.rotate(base, ang, origin="centroid")
        inter = g.intersection(dst).area
        score = inter / (g.area + dst.area - inter)
        if score > best:
            best, best_ang = score, ang
    return best_ang, dc - sc, float(best)


def rail_strips(gdf, village, line):
    """The parcels that ARE the railway land here, with where each sits across the track.

    A survey counts only if it carries at least 50 m of the traced centreline (the rule the engine
    already uses) and is no wider than RAIL_MAX_WIDTH_M across the track. Without the width test a
    240 m block merely crossed by the line is counted and the measurement is meaningless.
    """
    known = engine._rail_parcels(village)
    out = []
    for i, row in gdf.iterrows():
        if str(row["survey_no"]) not in known:
            continue
        g = _valid(row.geometry)
        if g is None or g.distance(line) > 5.0:
            continue
        o = across_track(g, line)
        if len(o) < RAIL_MIN_SAMPLES:
            continue
        lo, hi = float(np.percentile(o, 5)), float(np.percentile(o, 95))
        width = hi - lo
        centre, source = (lo + hi) / 2.0, "strip centre"
        if width > RAIL_MAX_WIDTH_M:
            if width > RAIL_WIDE_MAX_M:
                continue
            # too wide for its centre to mean anything: ask the sheet where the track runs
            drawn = _drawn_track_offset(village, str(row["survey_no"]), g, line)
            if drawn is None:
                continue
            centre, source = drawn, "track drawn on the sheet"
        out.append({"index": i, "survey": str(row["survey_no"]), "centre": centre,
                    "width": width, "source": source, "at": np.array(g.centroid.coords[0])})
    return out


def _drawn_track_offset(village, survey, placed, line, samples=21):
    """How far the track this parcel's own sheet draws lies from the real track, in metres."""
    from shapely import affinity
    from shapely.geometry import Point
    t, poly = sheet_track(village, survey)
    if t is None or poly is None or placed is None or placed.is_empty:
        return None
    ang, shift, iou = rigid_fit(poly, placed)
    if iou < 0.2:
        return None
    tw = affinity.rotate(affinity.translate(t, shift[0], shift[1]), ang,
                         origin=Point(*placed.centroid.coords[0]))
    off = []
    for k in range(samples):
        q = tw.interpolate(k / (samples - 1.0), normalized=True)
        sdist = line.project(q)
        r = line.interpolate(sdist)
        t0 = line.interpolate(max(0.0, sdist - 15.0))
        t1 = line.interpolate(min(line.length, sdist + 15.0))
        d = np.array([t1.x - t0.x, t1.y - t0.y])
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        d /= n
        off.append(float((np.array([q.x, q.y]) - np.array([r.x, r.y])) @ np.array([-d[1], d[0]])))
    return float(np.median(off)) if off else None


def track_offsets(gdf, village, line, min_iou=0.45, samples=21):
    """For each rail parcel, how far the track its own sheet draws lies from the real track.

    The sheet is fitted rigidly onto the placed parcel and the drawn track carried across with it.
    This beats assuming the strip is centred on the track: in Kizhikaranai 171 the sheet puts the
    track 9.7 m off the strip's centre, and in Vandalur 273 7.2 m off.
    """
    from shapely import affinity
    from shapely.geometry import Point
    out = []
    for st in rail_strips(gdf, village, line):
        t, poly = sheet_track(village, st["survey"])
        if t is None or poly is None:
            continue
        placed = _valid(gdf.loc[st["index"], "geometry"])
        if placed is None or placed.is_empty:
            continue
        ang, shift, iou = rigid_fit(poly, placed)
        if iou < min_iou:
            continue
        tw = affinity.rotate(affinity.translate(t, shift[0], shift[1]), ang,
                             origin=Point(*placed.centroid.coords[0]))
        off = []
        for k in range(samples):
            q = tw.interpolate(k / (samples - 1.0), normalized=True)
            sdist = line.project(q)
            r = line.interpolate(sdist)
            t0 = line.interpolate(max(0.0, sdist - 15.0))
            t1 = line.interpolate(min(line.length, sdist + 15.0))
            d = np.array([t1.x - t0.x, t1.y - t0.y])
            n = np.linalg.norm(d)
            if n < 1e-9:
                continue
            d /= n
            off.append(float((np.array([q.x, q.y]) - np.array([r.x, r.y])) @ np.array([-d[1], d[0]])))
        if off:
            out.append({**st, "drawn_off": float(np.median(off)), "iou": iou})
    return out


def rail_control_from_sheets(gdf, village, line):
    """Control that puts the track each sheet draws onto the track that is really there.

    Across-track only, as ever: along the track a strip slides invisibly.
    """
    from shapely.geometry import Point
    rows = track_offsets(gdf, village, line)
    if not rows:
        return np.zeros((0, 2)), np.zeros((0, 2)), None
    med = float(np.median([r["drawn_off"] for r in rows]))
    P, D = [], []
    for r in rows:
        sdist = line.project(Point(r["at"]))
        t0 = line.interpolate(max(0.0, sdist - 15.0))
        t1 = line.interpolate(min(line.length, sdist + 15.0))
        d = np.array([t1.x - t0.x, t1.y - t0.y])
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        d /= n
        P.append(r["at"])
        D.append(np.array([-d[1], d[0]]) * (-med))
    return np.array(P), np.array(D), med


PARCEL_RAIL_ROUNDS = 6      # correcting one strip nudges the next, so repeat until settled
PARCEL_RAIL_TOL_M = 5.0     # a strip this far from the norm is corrected on its own; inside
                            # that the strips of a settled village scatter by about +-2.5 m anyway
PARCEL_RAIL_MAX_M = 18.0    # beyond this it is not a placement error but a different parcel
PARCEL_RAIL_WIDTH = (15.0, 36.0)


def place_rail_parcels_on_track(gdf, village, line, target_centre, tol=3.0, cap=30.0):
    """Move each railway parcel bodily onto the track. Returns (gdf, notes).

    The smooth field cannot do this where a village's rail parcels disagree with each other:
    Peramanur's sit 28 m apart across the track, so correcting one drags its neighbours off. Each
    is therefore translated on its own, across the track only, and the overlaps that creates
    between it and its neighbours are settled afterwards. Shapes are not touched.
    """
    from shapely import affinity
    from shapely.geometry import Point
    moved = []
    geoms = list(gdf.geometry)
    for st in rail_strips(gdf, village, line):
        off = st["centre"] - target_centre
        if abs(off) < tol or abs(off) > cap:
            continue
        sdist = line.project(Point(st["at"]))
        t0 = line.interpolate(max(0.0, sdist - 15.0))
        t1 = line.interpolate(min(line.length, sdist + 15.0))
        d = np.array([t1.x - t0.x, t1.y - t0.y])
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        d /= n
        shift = np.array([-d[1], d[0]]) * (-off)
        i = st["index"]
        pos = list(gdf.index).index(i)
        g = _valid(geoms[pos])
        if g is None:
            continue
        geoms[pos] = affinity.translate(g, float(shift[0]), float(shift[1]))
        moved.append((st["survey"], round(float(off), 1), st["source"]))
    if not moved:
        return gdf, []
    out = gdf.copy()
    out["geometry"] = geoms
    # a parcel put on the track on evidence outranks the neighbour it now lies on, so the
    # neighbour yields whatever it takes; the usual quarter-parcel guard would refuse 52 % and
    # leave the overlap standing
    placed = {s for s, _o, _src in moved}
    keep = [_valid(g) for g in out.geometry]
    for i, (idx_i, row) in enumerate(out.iterrows()):
        if str(row["survey_no"]) not in placed or keep[i] is None:
            continue
        for j, other in enumerate(keep):
            if j == i or other is None or str(out.iloc[j]["survey_no"]) in placed:
                continue
            inter = keep[i].intersection(other)
            if inter.is_empty or inter.area <= CLIP_MIN_M2:
                continue
            cut = _valid(other.difference(keep[i]))
            if cut is not None and not cut.is_empty:
                keep[j] = cut
    out["geometry"] = keep
    out, _fixed = clip_siblings(out)
    return out, moved


def rail_control_per_parcel(gdf, village, line, target_centre,
                            tol=PARCEL_RAIL_TOL_M, cap=PARCEL_RAIL_MAX_M):
    """Control for the individual strips still sitting off the track after the village is aligned.

    A village-wide correction moves the median; a strip that disagrees with its own neighbours
    stays off. 25 of the 102 rail parcels were 4 m or more off after the village stage. Each is
    given its own across-track control point, so it moves and the parcels beside it follow through
    the same field, which keeps the boundaries they share.

    Strips outside `PARCEL_RAIL_WIDTH` or more than `cap` off are left alone: at that point the
    parcel is not a mis-placed railway strip but a different piece of ground.
    """
    from shapely.geometry import Point
    P, D, fixed = [], [], []
    for st in rail_strips(gdf, village, line):
        off = st["centre"] - target_centre
        if abs(off) < tol or abs(off) > cap:
            continue
        if st.get("source") != "track drawn on the sheet" and not (
                PARCEL_RAIL_WIDTH[0] <= st["width"] <= PARCEL_RAIL_WIDTH[1]):
            continue
        sdist = line.project(Point(st["at"]))
        t0 = line.interpolate(max(0.0, sdist - 15.0))
        t1 = line.interpolate(min(line.length, sdist + 15.0))
        d = np.array([t1.x - t0.x, t1.y - t0.y])
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        d /= n
        P.append(st["at"])
        D.append(np.array([-d[1], d[0]]) * (-off))
        fixed.append((st["survey"], round(off, 1)))
    return np.array(P) if P else np.zeros((0, 2)), np.array(D) if D else np.zeros((0, 2)), fixed


def rail_control(gdf, village, line, target_centre):
    """Control that puts this village's railway land back across the track where it belongs.

    Only the across-track direction is used. Along the track a strip slides invisibly, so nothing
    can be said about it and nothing is claimed. Measured on 2026-09-22: villages anchored by the
    team's placements sit within a few metres of the corridor norm, while a village corrected only
    by its seams can be 10 to 21 m off (Peramanur).
    """
    from shapely.geometry import Point
    strips = rail_strips(gdf, village, line)
    if not strips:
        return np.zeros((0, 2)), np.zeros((0, 2)), None
    off = float(np.median([s["centre"] for s in strips])) - target_centre
    P, D = [], []
    for s in strips:
        sdist = line.project(Point(s["at"]))
        t0 = line.interpolate(max(0.0, sdist - 15.0))
        t1 = line.interpolate(min(line.length, sdist + 15.0))
        d = np.array([t1.x - t0.x, t1.y - t0.y])
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        d /= n
        P.append(s["at"])
        D.append(np.array([-d[1], d[0]]) * (-off))
    return np.array(P), np.array(D), off


def seam_control(target, fixed, reach=SEAM_REACH_M, step=SEAM_STEP_M, share=1.0):
    """Control points that pull one village's edge onto its neighbour's.

    Puvi digitises each village on its own, so along a shared boundary the two layers overlap by
    thousands of square metres or stand apart: 8363 m2 between Thailavaram and Potheri, a 19 m gap
    between Kattankulathur and Peramanur. Sampling the frontage and asking each sampled point to
    move onto the neighbour's edge gives the field the same kind of observation a hand placement
    gives it, so the seam closes by moving parcels rather than by stretching them.

    `share` is how much of the distance this side takes: 1.0 when the neighbour is already correct
    and must not move, 0.5 when neither side is an authority and they should meet in the middle.
    """
    if target is None or fixed is None or target.is_empty or fixed.is_empty:
        return np.zeros((0, 2)), np.zeros((0, 2))
    tb, fb = target.boundary, fixed.boundary
    n = max(2, int(tb.length / step))
    pts, disp = [], []
    for i in range(n + 1):
        t = i / n
        p = tb.interpolate(t, normalized=True)
        q = fb.interpolate(fb.project(p))
        d = p.distance(q)
        if d < 0.05 or d > reach:
            continue
        # the neighbour must lie across this piece of boundary, not along it: without this test a
        # point just past a shared corner is dragged sideways onto that corner
        ahead = tb.interpolate(min(1.0, t + step / max(tb.length, 1.0)), normalized=True)
        tangent = np.array([ahead.x - p.x, ahead.y - p.y])
        to_q = np.array([q.x - p.x, q.y - p.y])
        nt, nq = np.linalg.norm(tangent), np.linalg.norm(to_q)
        if nt < 1e-9 or nq < 1e-9:
            continue
        if abs(float(tangent @ to_q) / (nt * nq)) > 0.5:
            continue
        pts.append((p.x, p.y))
        disp.append((to_q[0] * share, to_q[1] * share))
    return np.array(pts) if pts else np.zeros((0, 2)), np.array(disp) if disp else np.zeros((0, 2))


CLIP_GUARD = 0.25      # never take more than this share of a parcel to settle a village seam
CLIP_MIN_M2 = 0.05


def clip_siblings(gdf):
    """Settle any overlap between two parcels of the same village; the larger one yields.

    A field strong enough to close a village seam can fold slightly where it turns, which left
    Peramanur with five overlapping pairs totalling 118 m2 on the whole-village run.
    """
    geoms = [_valid(g) if g is not None and not g.is_empty else g for g in gdf.geometry]
    idx = gpd.GeoSeries([g for g in geoms if g is not None]).sindex
    order = [i for i, g in enumerate(geoms) if g is not None]
    fixed = 0
    for pos, i in enumerate(order):
        for k in idx.query(geoms[i], predicate="intersects"):
            j = order[k]
            if j <= i:
                continue
            a, b = geoms[i], geoms[j]
            if a is None or b is None:
                continue
            inter = a.intersection(b)
            if inter.is_empty or inter.area <= CLIP_MIN_M2:
                continue
            big, small = (i, j) if a.area >= b.area else (j, i)
            if geoms[big].area > 0 and inter.area / geoms[big].area > CLIP_GUARD:
                continue
            cut = _valid(geoms[big].difference(geoms[small]))
            if cut is None or cut.is_empty:
                continue
            geoms[big] = cut
            fixed += 1
    gdf["geometry"] = [g if g is not None and not g.is_empty else o
                       for g, o in zip(geoms, gdf.geometry)]
    return gdf, fixed


def clip_village_overlaps(layers, authority):
    """Settle the land two villages both claim. Returns (layers, notes).

    Puvi draws each village separately, so along a shared boundary they overlap: 1071 m2 between
    Kizhikaranai and Thirukatchur before anything was moved. Where both villages are anchored to
    our own placements neither can be moved to fix it, so the disputed strip is clipped out of one
    side, parcel by parcel, and never by more than a quarter of a parcel.
    """
    notes = []
    bodies = {v: unary_union([g for g in gdf.geometry if g is not None and not g.is_empty])
              for v, gdf in layers.items()}
    for a, b in itertools.combinations(sorted(layers), 2):
        inter = bodies[a].intersection(bodies[b])
        if inter.is_empty or inter.area <= CLIP_MIN_M2:
            continue
        # the village with our placements keeps its ground; otherwise the larger one yields
        if authority.get(a) and not authority.get(b):
            yielder, keeper = b, a
        elif authority.get(b) and not authority.get(a):
            yielder, keeper = a, b
        else:
            yielder, keeper = (a, b) if bodies[a].area >= bodies[b].area else (b, a)
        gdf = layers[yielder]
        taken = 0.0
        refused = 0
        geoms = list(gdf.geometry)
        for i, g in enumerate(geoms):
            if g is None or g.is_empty:
                continue
            piece = g.intersection(bodies[keeper])
            if piece.is_empty or piece.area <= CLIP_MIN_M2:
                continue
            if g.area > 0 and piece.area / g.area > CLIP_GUARD:
                refused += 1
                continue
            cut = _valid(g.difference(bodies[keeper]))
            if cut is None or cut.is_empty:
                refused += 1
                continue
            geoms[i] = cut
            taken += piece.area
        gdf["geometry"] = geoms
        bodies[yielder] = unary_union([g for g in geoms if g is not None and not g.is_empty])
        notes.append("%s yielded %.0f m2 to %s%s"
                     % (yielder, taken, keeper, ", %d parcel(s) refused by the guard" % refused if refused else ""))
    return layers, notes


def _valid(geom):
    if geom is None or geom.is_empty:
        return None
    return geom if geom.is_valid else make_valid(geom)


def puvi_polygons(village):
    """Every Puvi survey polygon of a village, dissolved per survey, in EPSG:32644."""
    files = sorted((paths.PUVI / village[:2] / village[3:5] / village[6:] / "vector").glob("*_vector.shp"))
    if not files:
        return gpd.GeoDataFrame(columns=["key", "geometry"], geometry="geometry", crs=UTM)
    g = gpd.read_file(files[0]).to_crs(UTM)
    g["key"] = g["survey_no"].map(shiftfit.normalise)
    g["geometry"] = [_valid(x) for x in g.geometry]
    g = g[g.geometry.notna()]
    return g


def control_parcels(village, extra_files=()):
    """The team's placements for this village as {survey key: geometry}, in EPSG:32644."""
    out = {}
    for f in sorted(paths.vector_dir(village).glob("*_parcels_modified.gpkg")):
        key = shiftfit.normalise(f.name.split("_parcels_modified")[0])
        g = visible.read_hand(f)
        if g.crs is not None and g.crs.to_epsg() != UTM:
            g = g.to_crs(UTM)
        body = unary_union([x for x in (_valid(v) for v in g.geometry) if x is not None])
        if not body.is_empty:
            out[key] = body
    for path in extra_files:
        g = gpd.read_file(path)
        if g.crs is not None and g.crs.to_epsg() != UTM:
            g = g.to_crs(UTM)
        if "survey_no" not in g.columns:
            log.warning("%s has no survey_no column, skipped as control", path)
            continue
        if "village_code" in g.columns:
            # survey numbers repeat across villages: a control file only speaks for its own village
            g = g[g["village_code"].astype(str) == village]
            if not len(g):
                continue
        g["key"] = g["survey_no"].map(shiftfit.normalise)
        for key, sub in g.groupby("key"):
            if key in ("NAN", "NONE", ""):
                continue
            body = unary_union([x for x in (_valid(v) for v in sub.geometry) if x is not None])
            if not body.is_empty:
                out[key] = body
    return out


def in_buffer_keys(village):
    """Survey keys of the village that touch the 30 m rail buffer, from the team's own layer."""
    b = gpd.read_file(paths.BUFFER_GPKG, layer="vector_in_buffer_30m")
    b = b[b["village_code"].astype(str) == village]
    return {shiftfit.normalise(s) for s in b["survey_no"]}


def village_name(village, puvi):
    for col in ("village", "vill_name", "village_na"):
        if col in puvi.columns and len(puvi):
            return str(puvi[col].iloc[0])
    return village


def run_village(village, out_dir, buffer_only=True, extra_control=(), stamp=None, seams=(),
                rail_target=None):
    """Write the corrected polygons for one village. Returns the report row.

    `seams` is a list of (neighbour geometry, share) already placed, used only when the village has
    no hand control of its own: its edge is pulled onto theirs so the corridor is continuous.
    """
    puvi = puvi_polygons(village)
    if not len(puvi):
        log.warning("%s: no Puvi vector found", village)
        return None
    control = control_parcels(village, extra_control)
    scope = in_buffer_keys(village) if buffer_only else set(puvi["key"])
    targets = puvi[puvi["key"].isin(scope)].copy()
    if not len(targets):
        log.warning("%s: no Puvi polygon is in scope", village)
        return None

    pairs, rejected = [], []
    by_key = {k: sub.geometry.iloc[0] for k, sub in puvi.groupby("key")}
    for key, hand_geom in control.items():
        if key not in by_key:
            continue
        d = float(np.linalg.norm(np.array(hand_geom.centroid.coords[0])
                                 - np.array(by_key[key].centroid.coords[0])))
        if d > MAX_CONTROL_M:
            # the same survey number in two villages, or a mis-placed parcel: not control
            rejected.append((key, round(d)))
            continue
        pairs.append((key, hand_geom))
    if rejected:
        log.warning("%s: %d control parcel(s) ignored, more than %.0f m from their Puvi twin: %s",
                    village, len(rejected), MAX_CONTROL_M, sorted(rejected, key=lambda r: -r[1])[:5])
    cpoints, cdisp = [], []
    for key, hand_geom in pairs:
        p = by_key[key]
        cpoints.append(np.array(p.centroid.coords[0]))
        cdisp.append(np.array(hand_geom.centroid.coords[0]) - np.array(p.centroid.coords[0]))
    cpoints = np.array(cpoints) if cpoints else np.zeros((0, 2))
    cdisp = np.array(cdisp) if cdisp else np.zeros((0, 2))

    acc = shiftfit.accuracy(cpoints, cdisp)
    hand_points = cpoints.copy()          # before any seam or rail point joins the control
    seam_note = ""
    rail_note = ""
    # Stage one: railway land. A handful of rail strips cannot outvote hundreds of seam points in
    # one fit (the first run moved them 0.3 m of the 12 m needed), so the across-track correction is
    # applied on its own first and the seams then adjust what is left.
    if len(cpoints) == 0 and rail_target is not None:
        # the railway land of this village should straddle the track the way it does everywhere
        # else on the corridor. Only the across-track direction is claimed.
        try:
            rp, rd, off = rail_control(targets, village, _rail_line(), rail_target)
        except Exception as exc:
            rp, rd, off = np.zeros((0, 2)), np.zeros((0, 2)), None
            log.warning("%s: rail alignment skipped (%s)", village, exc)
        if len(rp) and abs(off) <= RAIL_TOL_M:
            log.info("%s: railway land is %+.1f m across the track, within %.0f m of the norm; left alone",
                     village, off + rail_target, RAIL_TOL_M)
            rp = np.zeros((0, 2))
        if len(rp):
            rail_table = _warp_table(_nodes(list(targets.geometry)), rp, rd,
                                     k=shiftfit.SEAM_K, smooth=shiftfit.SEAM_SMOOTH_M)
            before = [_valid(g) for g in targets.geometry]
            targets = targets.copy()
            targets["geometry"] = [_warp_with(g, rail_table) for g in before]
            moved = [float(np.linalg.norm(np.array(_valid(b).centroid.coords[0])
                                          - np.array(_valid(a).centroid.coords[0])))
                     for a, b in zip(before, targets.geometry) if _valid(a) is not None and _valid(b) is not None]
            rail_note = "%d rail strips, %+.1f m across the track, parcels moved %.1f m" % (
                len(rp), -off, float(np.median(moved)) if moved else 0.0)
            log.info("%s: railway land aligned, %s", village, rail_note)

    if len(hand_points) == 0 and seams:
        # no placement of ours anywhere in this village: the only thing known about it is that its
        # edge must meet the villages already corrected beside it
        body = unary_union(list(targets.geometry))
        sp, sd = [], []
        for neighbour, share in seams:
            a, b = seam_control(body, neighbour, share=share)
            if len(a):
                sp.append(a); sd.append(b)
        if sp:
            # add to whatever control this village already has: the railway-land points must not be
            # thrown away here, which is what happened on the first run and left the strips untouched
            add_p, add_d = np.vstack(sp), np.vstack(sd)
            cpoints = np.vstack([cpoints, add_p]) if len(cpoints) else add_p
            cdisp = np.vstack([cdisp, add_d]) if len(cdisp) else add_d
            seam_note = "%d seam points against %d corrected neighbour(s)" % (len(add_p), len(sp))
            log.info("%s: no control of our own, %s", village, seam_note)

    tpoints = np.array([g.centroid.coords[0] for g in targets.geometry])
    fit_kw = ({"k": shiftfit.SEAM_K, "smooth": shiftfit.SEAM_SMOOTH_M} if (seam_note or rail_note) else {})
    shifts, method = shiftfit.fit(cpoints, cdisp, tpoints, **fit_kw)
    if seam_note and rail_note:
        method = "railway land and seams"
    elif seam_note:
        method = "seam to corrected neighbour"
    elif rail_note:
        method = "railway land"

    # every parcel is warped through the same field, including the ones the team has placed by
    # hand: mixing two sources of geometry in one layer is what draws a double line along every
    # boundary between them. Their own files stay the authority and are untouched.
    control_keys = {k for k, _ in pairs}
    table = _warp_table(_nodes(targets.geometry), cpoints, cdisp, **fit_kw) if len(cpoints) else {}
    moved, source = [], []
    for (i, row), s in zip(targets.iterrows(), shifts):
        field = (lambda q: shiftfit.fit(cpoints, cdisp, q[None, :], **fit_kw)[0][0]) if len(cpoints) else None
        moved.append(_warp_with(row.geometry, table, field=field) if table else row.geometry)
        source.append("hand placed by the team" if row["key"] in control_keys else "puvi")
    out = targets.copy()
    out["geometry"] = moved
    out["village_code"] = village
    out["survey_no"] = out["key"]
    out["shift_x_m"] = [round(float(s[0]), 2) for s in shifts]
    out["shift_y_m"] = [round(float(s[1]), 2) for s in shifts]
    out["shift_m"] = [round(float(np.hypot(*s)), 2) for s in shifts]
    # evidence means nearness to a parcel the team placed, never to a seam point: a seam says two
    # villages touch, which is a statement about the fabric and not about where either one truly is
    if len(hand_points):
        near = np.array([float(np.min(np.linalg.norm(hand_points - q, axis=1))) for q in tpoints])
    else:
        near = np.full(len(tpoints), np.inf)
    out["source"] = source
    out["control_within_m"] = [round(float(d), 1) if np.isfinite(d) else "" for d in near]
    out["evidence"] = ["measured" if d <= MEASURED_REACH_M else "not evidenced" for d in near]
    out["hand_placed"] = ["yes" if k in control_keys else "no" for k in targets["key"]]
    out["fit_method"] = method
    out["control_parcels"] = len(cpoints)
    out["expected_error_m"] = acc["after_median_m"] if acc["after_median_m"] is not None else ""
    out["puvi_error_before_m"] = acc["before_median_m"] if acc["before_median_m"] is not None else ""
    out["run_at"] = stamp or shiftfit.run_stamp()

    keep = ["village_code", "survey_no", "source", "hand_placed", "evidence", "control_within_m",
            "fit_method", "control_parcels",
            "shift_x_m", "shift_y_m", "shift_m", "expected_error_m", "puvi_error_before_m",
            "run_at", "geometry"]
    # Stage three: the seams pull on the railway land too, so check it once more and correct what
    # is left. Without this a village can end further off the track than it started (Peramanur went
    # +14.7 to +15.7 m, Vinchiyambakkam +2.8 to +10.7 m on the 2026-09-22 run).
    # run whatever the first stage decided: a village left alone because Puvi already had its
    # railway land right still needs checking once the seams have pulled on it. Skipping this when
    # the first stage did nothing left Peramanur at +14.7 m across the track when Puvi had it at +0.2.
    if rail_target is not None and len(hand_points) == 0:
        try:
            line2 = _rail_line()
            rp2, rd2, off2 = rail_control(out, village, line2, rail_target)
            if len(rp2) and abs(off2) > RAIL_TOL_M:
                t2 = _warp_table(_nodes(list(out.geometry)), rp2, rd2,
                                 k=shiftfit.SEAM_K, smooth=shiftfit.SEAM_SMOOTH_M)
                out = out.copy()
                out["geometry"] = [_warp_with(_valid(g), t2) for g in out.geometry]
                rail_note += "; settled %+.1f m after the seams" % -off2
                log.info("%s: railway land settled %+.1f m after the seams", village, -off2)
            # and finally the strips that still disagree with their own neighbours, one by one.
            # Correcting one moves the parcels beside it, which nudges the next strip, so this is
            # repeated until nothing is left outside the tolerance.
            done = 0
            for _round in range(PARCEL_RAIL_ROUNDS):
                pp, pd, fixed = rail_control_per_parcel(out, village, line2, rail_target)
                if not len(pp):
                    break
                t3 = _warp_table(_nodes(list(out.geometry)), pp, pd)
                out = out.copy()
                out["geometry"] = [_warp_with(_valid(g), t3) for g in out.geometry]
                done += len(fixed)
                log.info("%s: round %d, %d rail strip(s) corrected on their own: %s",
                         village, _round + 1, len(fixed), fixed[:6])
            if done:
                rail_note += "; %d strip correction(s) over %d round(s)" % (done, _round + 1)
        except Exception as exc:
            log.warning("%s: second rail pass skipped (%s)", village, exc)
        out["fit_method"] = method

    out = out[[c for c in keep if c in out.columns]]
    out, fixed = clip_siblings(out)
    if fixed:
        log.info("%s: %d overlap(s) between parcels of this village settled", village, fixed)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / ("%s_puvi_shifted.geojson" % village)
    out.to_crs(4326).to_file(dest, driver="GeoJSON", COORDINATE_PRECISION=8)
    log.info("%s: %d polygons -> %s (%s, %d control)", village, len(out), dest, method, len(cpoints))

    return {"village_code": village, "village_name": village_name(village, puvi),
            "control": len(cpoints), "targets": len(out), "method": method,
            "seam_note": "; ".join(x for x in (rail_note, seam_note) if x),
            "mean_shift_x_m": round(float(cdisp.mean(0)[0]), 2) if len(cdisp) else "",
            "mean_shift_y_m": round(float(cdisp.mean(0)[1]), 2) if len(cdisp) else "",
            "run_at": out["run_at"].iloc[0], **{k: v for k, v in acc.items() if k != "control"}}


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="shift_puvi", description=__doc__.splitlines()[0])
    ap.add_argument("villages", nargs="*", help="village codes, e.g. 35_04_074")
    ap.add_argument("--stretch", action="store_true", help="the Thailavaram to Thirukatchur villages")
    ap.add_argument("--buffer-only", action="store_true", default=True,
                    help="only surveys touching the 30 m rail buffer (default)")
    ap.add_argument("--whole-village", dest="buffer_only", action="store_false",
                    help="every survey of the village, not just the buffer")
    ap.add_argument("--control", action="append", default=[],
                    help="extra control file with a survey_no column (repeatable)")
    ap.add_argument("--no-rail", action="store_true",
                    help="do not align a village's railway land across the track")
    ap.add_argument("--no-seams", action="store_true",
                    help="do not pull a village with no control onto its corrected neighbours")
    ap.add_argument("--corridor", action="store_true",
                    help="every village the rail buffer passes through that has a Puvi vector")
    ap.add_argument("--buffer-layer", action="store_true",
                    help="deliver only the parcels touching the rail buffer; whole villages are still "
                         "processed, into a _working folder, because the seams need them")
    ap.add_argument("--no-clip", action="store_true",
                    help="leave land that two villages both claim as it is")
    ap.add_argument("--out", default="", help="output folder (default: <project>\\Puvi_Shifted\\<date>)")
    args = ap.parse_args(argv)

    villages = list(args.villages) + (STRETCH if args.stretch else []) + (CORRIDOR if args.corridor else [])
    if not villages:
        ap.error("name at least one village, or pass --stretch or --corridor")
    stamp = shiftfit.run_stamp()
    out_dir = __import__("pathlib").Path(args.out) if args.out else shiftfit.output_dir()
    villages = list(dict.fromkeys(villages))
    has_control = {}
    for v in villages:
        puvi = puvi_polygons(v)
        keys = set(puvi["key"])
        has_control[v] = any(k in keys for k in control_parcels(v, args.control))
    order = [v for v in villages if has_control[v]] + [v for v in villages if not has_control[v]]
    if order != villages:
        log.info("villages with our own placements go first: %s", order)

    # whole villages are processed even when only the buffer parcels are wanted: the seams and the
    # topology between villages only come out right on the complete fabric. Those intermediate
    # files go to _working and the deliverable is written beside it.
    work_dir = (out_dir / "_working") if args.buffer_layer else out_dir
    work_dir.mkdir(parents=True, exist_ok=True)

    rail_target = None

    def corridor_norm():
        """Where railway land sits across the track in the villages our placements anchor.

        Measured on their corrected geometry, not on raw Puvi: the whole point is to bring the
        other villages to where the anchored ones ended up.
        """
        if args.no_rail:
            return None
        try:
            line = _rail_line()
        except Exception as exc:
            log.warning("rail alignment unavailable: %s", exc)
            return None
        seen = []
        for v in order:
            if not has_control[v]:
                continue
            f = work_dir / ("%s_puvi_shifted.geojson" % v)
            if not f.exists():
                continue
            g = gpd.read_file(f).to_crs(UTM)
            strips = rail_strips(g, v, line)
            if strips:
                seen.append(float(np.median([x["centre"] for x in strips])))
        if not seen:
            return None
        t = float(np.median(seen))
        log.info("railway land sits %+.2f m across the track where your placements anchor it "
                 "(%d village(s)); villages without control are brought to that", t, len(seen))
        return t

    rows, placed = [], {}
    for v in order:
        if not has_control[v] and rail_target is None:
            rail_target = corridor_norm()
        seams = []
        if not has_control[v] and not args.no_seams:
            body = unary_union(list(puvi_polygons(v).geometry))
            for other, geom in placed.items():
                if body.distance(geom) <= SEAM_REACH_M:
                    seams.append((geom, 1.0 if has_control[other] else 0.5))
        row = run_village(v, work_dir, buffer_only=args.buffer_only,
                          extra_control=args.control, stamp=stamp, seams=seams, rail_target=rail_target)
        if row:
            rows.append(row)
            written = gpd.read_file(work_dir / ("%s_puvi_shifted.geojson" % v)).to_crs(UTM)
            placed[v] = unary_union([x for x in (_valid(g) for g in written.geometry) if x is not None])
    if not rows:
        log.error("nothing written")
        return 1
    # a village placed before its neighbours only saw some of them: run the seam pass again now
    # that every village has a position
    if not args.no_seams:
        for v in [x for x in order if not has_control[x]]:
            body = unary_union(list(puvi_polygons(v).geometry))
            seams = [(geom, 1.0 if has_control[other] else 0.5)
                     for other, geom in placed.items()
                     if other != v and body.distance(geom) <= SEAM_REACH_M]
            if not seams:
                continue
            row = run_village(v, work_dir, buffer_only=args.buffer_only,
                              extra_control=args.control, stamp=stamp, seams=seams, rail_target=rail_target)
            if row:
                rows = [r for r in rows if r["village_code"] != v] + [row]
                written = gpd.read_file(work_dir / ("%s_puvi_shifted.geojson" % v)).to_crs(UTM)
                placed[v] = unary_union([x for x in (_valid(g) for g in written.geometry) if x is not None])

    # land two villages both claim, which no amount of moving can settle
    layers = {r["village_code"]: gpd.read_file(work_dir / ("%s_puvi_shifted.geojson" % r["village_code"])).to_crs(UTM)
              for r in rows}
    if not args.no_clip:
        layers, notes = clip_village_overlaps(layers, has_control)
        for n in notes:
            log.info("village seam: %s", n)
        # settling land two villages both claim moves their railway strips too, and nothing had
        # checked that afterwards: it left Perunkalathur and Peerkkararanai further off the track
        # than Puvi had them. So the track is checked once more, last of all.
        if rail_target is not None:
            try:
                line3 = _rail_line()
            except Exception:
                line3 = None
            for v, gdf in layers.items():
                if line3 is None or has_control.get(v):
                    continue
                for _round in range(PARCEL_RAIL_ROUNDS):
                    rp3, rd3, off3 = rail_control(gdf, v, line3, rail_target)
                    if not len(rp3) or abs(off3) <= RAIL_TOL_M:
                        break
                    t4 = _warp_table(_nodes(list(gdf.geometry)), rp3, rd3,
                                     k=shiftfit.SEAM_K, smooth=shiftfit.SEAM_SMOOTH_M)
                    gdf = gdf.copy()
                    gdf["geometry"] = [_warp_with(_valid(g), t4) for g in gdf.geometry]
                    layers[v] = gdf
                    log.info("%s: railway land settled %+.1f m after the village seams were clipped",
                             v, -off3)
                # and the individual strips the village median cannot reach. These are moved
                # bodily: a field smooth enough to keep the fabric would drag their neighbours off
                # the track with them.
                for _round in range(PARCEL_RAIL_ROUNDS):
                    gdf, moved = place_rail_parcels_on_track(gdf, v, line3, rail_target)
                    if not moved:
                        break
                    layers[v] = gdf
                    log.info("%s: %d rail parcel(s) moved onto the track: %s", v, len(moved), moved[:5])

        # moving a rail parcel bodily leaves it lying on its neighbours, so the claims inside each
        # village and between villages are settled again, and anything the differences broke is
        # repaired before it is written
        for v, gdf in layers.items():
            gdf, _n = clip_siblings(gdf)
            gdf = gdf.copy()
            gdf["geometry"] = [_repair(_valid(g)) if g is not None else g for g in gdf.geometry]
            layers[v] = gdf
        layers, notes2 = clip_village_overlaps(layers, has_control)
        for n in notes2:
            log.info("village seam, second pass: %s", n)

        for v, gdf in layers.items():
            out_g = gdf.to_crs(4326)
            # a parcel valid in metres can come back self-intersecting once written as degrees,
            # so it is snapped to a millimetre grid and repaired in the CRS it is written in
            out_g["geometry"] = [_snap_or_keep(g) for g in out_g.geometry]
            out_g["geometry"] = [g if g is None or g.is_valid else _repair(g.buffer(0))
                                 for g in out_g.geometry]
            out_g.to_file(work_dir / ("%s_puvi_shifted.geojson" % v),
                          driver="GeoJSON", COORDINATE_PRECISION=8)

    rows = sorted(rows, key=lambda r: order.index(r["village_code"]))
    if args.buffer_layer:
        parts = []
        for r in rows:
            v = r["village_code"]
            f = work_dir / ("%s_puvi_shifted.geojson" % v)
            if not f.exists():
                continue
            g = gpd.read_file(f)
            try:
                keys = in_buffer_keys(v)
            except Exception:
                keys = set()
            sub = g[g["survey_no"].astype(str).isin(keys)]
            if not len(sub):
                continue
            sub.to_file(out_dir / ("%s_buffer_parcels.geojson" % v), driver="GeoJSON",
                        COORDINATE_PRECISION=8)
            parts.append(sub)
            r["buffer_parcels"] = len(sub)
        if parts:
            buf = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=4326)
            buf.to_file(out_dir / "rail_buffer_parcels.geojson", driver="GeoJSON",
                        COORDINATE_PRECISION=8)
            log.info("deliverable: %d cadastral parcels in the rail buffer over %d villages -> %s",
                     len(buf), buf["village_code"].nunique(), out_dir / "rail_buffer_parcels.geojson")
    else:
        merged = pd.concat([gpd.read_file(out_dir / ("%s_puvi_shifted.geojson" % r["village_code"]))
                            for r in rows], ignore_index=True)
        gpd.GeoDataFrame(merged, crs=4326).to_file(out_dir / "puvi_shifted_all.geojson",
                                                   driver="GeoJSON", COORDINATE_PRECISION=8)
    shiftfit.write_report(rows, out_dir / "shift_report.csv")
    print("\n%-16s %8s %8s %-12s %10s %10s" % ("village", "control", "parcels", "method", "before", "after"))
    for r in rows:
        print("%-16s %8d %8d %-28s %9s m %9s m"
              % (r["village_name"][:16], r["control"], r["targets"], r["method"],
                 r["before_median_m"] if r["before_median_m"] is not None else "-",
                 r["after_median_m"] if r["after_median_m"] is not None else "-"))
    print("\nwritten to %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
