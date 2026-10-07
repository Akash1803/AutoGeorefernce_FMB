"""Topology for a placed ring, and its two independent checks.

clean_batch: the batch's plots close onto Akash's parcels with neat.clean - his parcels are fixed and
never change; a batch corner near his corner or edge moves onto it, overlaps are cut along his line,
wedges between PDF neighbours close. After the round trip to EPSG:4326 the batch's corners are put back
exactly on his vertices so the shared seams are identical coordinate for coordinate.

Dual verification (Akash, 2026-10-07: "always dual verification"):
  1. verify_files  - shapely, on the files as written;
  2. verify_qgis   - QGIS's own geometry engine, on the layers as QGIS shows them.
Both count overlaps between any two plots and the gaps (holes) the batch leaves against its neighbours.
"""
import json

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.ops import unary_union

from . import neat, qgis_bridge

UTM = 32644
M2_PER_DEG2 = 110540.0 * 108600.0   # near 12.7 N
SEAM_SNAP_DEG = 2e-7             # ~2 cm: batch vertices this close to his line take his vertices
SNAP_DEG = 1e-8          # ~1 mm: a batch corner this close to one of his corners becomes exactly that corner


def clean_batch(seed_plots, batch_plots, agreed):
    """(cleaned batch GeoDataFrame in EPSG:4326, neat report). Both inputs in UTM with survey_no."""
    seeds = seed_plots.assign(_seed=True)
    batch = batch_plots.assign(_seed=False)
    both = gpd.GeoDataFrame(pd.concat([seeds, batch], ignore_index=True), geometry="geometry", crs=UTM)
    both["survey_no"] = both["survey_no"].astype(str)
    printed = {k: set(v) for k, v in agreed.items()}
    out, rep = neat.clean(both, printed, fixed_surveys=sorted(set(seeds.survey_no.astype(str))),
                          first_surveys=())
    for i in np.where(both._seed.values)[0]:                 # his parcels: untouched, or stop
        assert out.geometry.iloc[i].equals(both.geometry.iloc[i]), "a seed parcel changed"
    cleaned = out[~out._seed.values].drop(columns="_seed")
    rep = rep[~both._seed.values]
    return cleaned, rep


def snap_to_seeds(batch_ll, seed_ll):
    """Batch corners within SNAP_DEG of a seed corner take that corner's exact coordinates (EPSG:4326)."""
    pts = shapely.get_coordinates(seed_ll.geometry.values)
    if not len(pts):
        return batch_ll
    from scipy.spatial import cKDTree
    tree = cKDTree(pts)

    def fix(g):
        def f(c):
            d, k = tree.query(c)
            c2 = c.copy()
            near = d <= SNAP_DEG
            c2[near] = pts[k[near]]
            return c2
        return shapely.transform(g, f)
    out = batch_ll.copy()
    out["geometry"] = [fix(g) for g in batch_ll.geometry]
    return out


SEAM_SLIVER_M = 0.6     # a hairline strip between two batch plots (any length) closes into the one it borders most
MICRO_GAP_M2 = 1.0      # a hole this small between a batch plot and its neighbours joins the batch plot


def final_pass(seed_ll, batch_ll):
    """Exact last pass in the layers' own CRS: what a batch plot still overlaps (his parcel or another batch
    plot) is cut from it along the other plot's line, and hairline holes against the batch are filled
    into the batch plot that borders them most. His parcels are never touched."""
    b = batch_ll.copy().reset_index(drop=True)
    g = list(b.geometry.values)
    sg = list(seed_ll.geometry.values)
    st = shapely.STRtree(sg)
    for i in range(len(g)):                 # seams take his exact vertices (a hair-width residue otherwise)
        near = [sg[j] for j in st.query(g[i].buffer(SEAM_SNAP_DEG))]
        if near:
            s2 = shapely.snap(g[i], unary_union([x.boundary for x in near]), SEAM_SNAP_DEG)
            if s2.is_valid and abs(s2.area - g[i].area) < 1e-3 / M2_PER_DEG2 + 1e-4 * g[i].area:
                g[i] = s2
    for i in range(len(g)):
        for j in st.query(g[i], predicate="intersects"):
            x = shapely.intersection(g[i], sg[j])
            if sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") > 0:
                g[i] = _largest(shapely.difference(g[i], sg[j]))
    for i in range(len(g)):
        for j in range(i + 1, len(g)):
            if g[i].intersects(g[j]):
                x = shapely.intersection(g[i], g[j])
                if sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") > 0:
                    small, big = (i, j) if g[i].area < g[j].area else (j, i)
                    g[small] = _largest(shapely.difference(g[small], g[big]))
    m2 = 110540.0 * 108600.0
    for _round in range(3):                  # filling one strip can leave a hairline beside it
        g = [despike(x) for x in g]
        U = unary_union(g + sg)
        filled = 0
        for p in getattr(U, "geoms", [U]):
            for r in p.interiors:
                h = shapely.Polygon(r).difference(U)
                if h.is_empty:
                    continue
                width_m = 2 * h.area * m2 / max(h.length * 110000.0, 1e-9)
                between_batch_only = not any(x.distance(h) < 1e-9 for x in sg)
                if h.area * m2 > MICRO_GAP_M2 and not (between_batch_only and width_m <= SEAM_SLIVER_M):
                    continue                    # real open ground, or a gap against his parcels: left for him
                touch = [(h.boundary.intersection(g[i].buffer(1e-9)).length, i) for i in range(len(g)) if g[i].distance(h) < 1e-9]
                if touch:
                    _l, i = max(touch)
                    g[i] = _largest(shapely.union(g[i], h))
                    filled += 1
        if not filled:
            break
    b["geometry"] = [despike(x) for x in g]
    return b


def despike(poly, tol=2e-8):
    """Drop zero-width spikes (a -> b -> a) and repeated corners left by snapping (264, 2026-10-07: a 4 mm
    spike made the plot invalid once projected to metres). tol is in degrees (~2 mm)."""
    if poly.geom_type != "Polygon":
        return poly
    pts = [tuple(p) for p in shapely.get_coordinates(poly.exterior)][:-1]
    changed = True
    while changed and len(pts) > 3:
        changed = False
        for i in range(len(pts)):
            a, b, c = pts[i - 1], pts[i], pts[(i + 1) % len(pts)]
            if np.hypot(a[0] - b[0], a[1] - b[1]) < tol:
                del pts[i]; changed = True; break
            if np.hypot(a[0] - c[0], a[1] - c[1]) < tol:
                del pts[i]                      # b is the tip of a spike; a and c are the same corner
                del pts[i % len(pts)]
                changed = True; break
    out = shapely.Polygon(pts)
    return out if out.is_valid and abs(out.area - poly.area) < 1e-12 + 1e-6 * poly.area else poly


def _largest(geom):
    parts = [p for p in getattr(geom, "geoms", [geom]) if p.geom_type == "Polygon" and not p.is_empty]
    return max(parts, key=lambda p: p.area) if parts else geom


def verify_files(seed_ll, batch_ll):
    """Check 1 (shapely on the written data): overlaps, gaps against neighbours, validity."""
    s = seed_ll.assign(_who="his")
    b = batch_ll.assign(_who="auto")
    allp = gpd.GeoDataFrame(pd.concat([s, b], ignore_index=True), geometry="geometry", crs=4326)
    g = allp.geometry.values          # the layers' own CRS, as QGIS's Topology Checker sees them
    k2 = M2_PER_DEG2
    is_b = (allp._who == "auto").values
    tree = shapely.STRtree(g)
    a_, b_ = tree.query(g, predicate="intersects")
    keep = (a_ < b_) & (is_b[a_] | is_b[b_])
    ov = []
    for i, j in zip(a_[keep], b_[keep]):
        x = shapely.intersection(g[i], g[j])
        ar = sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon")
        if ar > 0:
            ov.append((str(allp.survey_no.iat[i]), str(allp.survey_no.iat[j]), round(ar * k2, 4)))
    near = unary_union([g[k] for k in np.where(is_b)[0]]).buffer(50 / 110000.0)
    zone = [k for k in range(len(g)) if g[k].intersects(near)]
    U = unary_union([g[k] for k in zone])
    holes = []
    for p in getattr(U, "geoms", [U]):
        for r in p.interiors:
            h = shapely.Polygon(r).difference(U)
            if h.area > 0 and h.distance(unary_union([g[k] for k in np.where(is_b)[0]])) < 1e-7:
                holes.append(round(h.area * k2, 4))
    bad = [str(allp.survey_no.iat[k]) for k in np.where(is_b)[0] if not g[k].is_valid or g[k].geom_type != "Polygon"]
    return {"overlaps": len(ov), "overlap_m2": round(sum(x[2] for x in ov), 4), "overlap_pairs": ov[:20],
            "gaps": len(holes), "gap_m2": round(sum(holes), 4), "invalid_or_multipart": bad}


def verify_qgis(village_name, batch_surveys):
    """Check 2 (QGIS geometry engine, on the layers as QGIS shows them)."""
    code = """
import json
from qgis.core import (QgsProject, QgsCoordinateTransform, QgsCoordinateReferenceSystem, QgsGeometry,
                       QgsSpatialIndex, QgsFeature)
P = QgsProject.instance()
sub = P.layerTreeRoot().findGroup("Manual_Georeferenced").findGroup(%(sub)r)
batch = set(%(batch)r)
K2 = 110540.0 * 108600.0
feats = []
for n in sub.findLayers():
    L = n.layer()
    if L is None or "_parcels" not in L.name():
        continue
    s = L.name().split("_parcels")[0]
    is_auto = L.name().endswith("_parcels_auto") and s in batch     # other red layers are neighbours like his
    for f in L.getFeatures():
        g = QgsGeometry(f.geometry())
        if g.isEmpty():
            continue
        if L.crs().authid() != "EPSG:4326":
            g.transform(QgsCoordinateTransform(L.crs(), QgsCoordinateReferenceSystem("EPSG:4326"), P))
        feats.append((s, is_auto, g))
idx = QgsSpatialIndex()
for k, (s, a, g) in enumerate(feats):
    f = QgsFeature(k); f.setGeometry(g); idx.addFeature(f)
ov, bad = [], []
for i, (s, a, g) in enumerate(feats):
    if not a:
        continue
    if not g.isGeosValid() or g.isMultipart():
        bad.append(s)
    for j in idx.intersects(g.boundingBox()):
        if j == i or (feats[j][1] and j < i):
            continue
        x = g.intersection(feats[j][2])
        if x and not x.isEmpty() and x.area() > 0:
            ov.append([s, feats[j][0], round(x.area() * K2, 4)])
zone = QgsGeometry.unaryUnion([g for s, a, g in feats if a]).buffer(50 / 110000.0, 4)
U = QgsGeometry.unaryUnion([g for s, a, g in feats if g.intersects(zone)])
B = QgsGeometry.unaryUnion([g for s, a, g in feats if a])
holes = []
for part in U.asGeometryCollection():
    poly = part.asPolygon()
    for ring in poly[1:]:
        h = QgsGeometry.fromPolygonXY([ring]).difference(U)
        if h.area() > 0 and h.distance(B) < 1e-7:
            holes.append(round(h.area() * K2, 4))
print(json.dumps({"overlaps": len(ov), "overlap_m2": round(sum(x[2] for x in ov), 4), "overlap_pairs": ov[:20],
                  "gaps": len(holes), "gap_m2": round(sum(holes), 4), "invalid_or_multipart": bad,
                  "plots_checked": len(feats)}))
""" % {"sub": "%s_Village" % village_name, "batch": list(batch_surveys)}
    return qgis_bridge.json_result(code, timeout=300)
