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
from shapely import affinity
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union
from shapely.validation import make_valid

from . import paths, shiftfit

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


def _warp_table(nodes, points, disp, grid=NODE_GRID_M):
    """Where each distinct node lands. One field evaluation per node, not per polygon."""
    keys = list(nodes)
    coords = np.array([nodes[k] for k in keys], float)
    if not len(coords):
        return {}
    shifts, _ = shiftfit.fit(points, disp, coords)
    return {k: tuple(c + s) for k, c, s in zip(keys, coords, shifts)}


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
        g = gpd.read_file(f)
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


def run_village(village, out_dir, buffer_only=True, extra_control=(), stamp=None, seams=()):
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
    hand_points = cpoints.copy()          # before any seam point joins the control
    seam_note = ""
    if len(cpoints) == 0 and seams:
        # no placement of ours anywhere in this village: the only thing known about it is that its
        # edge must meet the villages already corrected beside it
        body = unary_union(list(targets.geometry))
        sp, sd = [], []
        for neighbour, share in seams:
            a, b = seam_control(body, neighbour, share=share)
            if len(a):
                sp.append(a); sd.append(b)
        if sp:
            cpoints = np.vstack(sp)
            cdisp = np.vstack(sd)
            seam_note = "%d seam points against %d corrected neighbour(s)" % (len(cpoints), len(sp))
            log.info("%s: no control of our own, %s", village, seam_note)

    tpoints = np.array([g.centroid.coords[0] for g in targets.geometry])
    shifts, method = shiftfit.fit(cpoints, cdisp, tpoints)
    if seam_note:
        method = "seam to corrected neighbour"

    # every parcel is warped through the same field, including the ones the team has placed by
    # hand: mixing two sources of geometry in one layer is what draws a double line along every
    # boundary between them. Their own files stay the authority and are untouched.
    control_keys = {k for k, _ in pairs}
    table = _warp_table(_nodes(targets.geometry), cpoints, cdisp) if len(cpoints) else {}
    moved, source = [], []
    for (i, row), s in zip(targets.iterrows(), shifts):
        field = (lambda q: shiftfit.fit(cpoints, cdisp, q[None, :])[0][0]) if len(cpoints) else None
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
            "seam_note": seam_note,
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
    ap.add_argument("--no-seams", action="store_true",
                    help="do not pull a village with no control onto its corrected neighbours")
    ap.add_argument("--corridor", action="store_true",
                    help="every village the rail buffer passes through that has a Puvi vector")
    ap.add_argument("--buffer-layer", action="store_true",
                    help="also write just the parcels that touch the rail buffer, cut from the result")
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

    rows, placed = [], {}
    for v in order:
        seams = []
        if not has_control[v] and not args.no_seams:
            body = unary_union(list(puvi_polygons(v).geometry))
            for other, geom in placed.items():
                if body.distance(geom) <= SEAM_REACH_M:
                    seams.append((geom, 1.0 if has_control[other] else 0.5))
        row = run_village(v, out_dir, buffer_only=args.buffer_only,
                          extra_control=args.control, stamp=stamp, seams=seams)
        if row:
            rows.append(row)
            written = gpd.read_file(out_dir / ("%s_puvi_shifted.geojson" % v)).to_crs(UTM)
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
            row = run_village(v, out_dir, buffer_only=args.buffer_only,
                              extra_control=args.control, stamp=stamp, seams=seams)
            if row:
                rows = [r for r in rows if r["village_code"] != v] + [row]
                written = gpd.read_file(out_dir / ("%s_puvi_shifted.geojson" % v)).to_crs(UTM)
                placed[v] = unary_union([x for x in (_valid(g) for g in written.geometry) if x is not None])

    # land two villages both claim, which no amount of moving can settle
    layers = {r["village_code"]: gpd.read_file(out_dir / ("%s_puvi_shifted.geojson" % r["village_code"])).to_crs(UTM)
              for r in rows}
    if not args.no_clip:
        layers, notes = clip_village_overlaps(layers, has_control)
        for n in notes:
            log.info("village seam: %s", n)
        for v, gdf in layers.items():
            gdf.to_crs(4326).to_file(out_dir / ("%s_puvi_shifted.geojson" % v),
                                     driver="GeoJSON", COORDINATE_PRECISION=8)

    if args.buffer_layer:
        parts = []
        for v in order:
            f = out_dir / ("%s_puvi_shifted.geojson" % v)
            if not f.exists():
                continue
            g = gpd.read_file(f)
            try:
                keys = in_buffer_keys(v)
            except Exception:
                keys = set()
            sub = g[g["survey_no"].astype(str).isin(keys)]
            if len(sub):
                parts.append(sub)
        if parts:
            buf = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=4326)
            buf.to_file(out_dir / "puvi_shifted_rail_buffer.geojson", driver="GeoJSON",
                        COORDINATE_PRECISION=8)
            log.info("parcels in the rail buffer: %d -> %s", len(buf),
                     out_dir / "puvi_shifted_rail_buffer.geojson")

    rows = sorted(rows, key=lambda r: order.index(r["village_code"]))
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
