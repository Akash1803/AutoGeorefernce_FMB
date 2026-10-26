"""One Stretch, step 3: the cleaned plots on the layer's coordinates (lon/lat, 10 decimals = 0.01 mm), seams noded so that
both sides of a seam carry the same vertices. Manual_Georeference plots keep the exact coordinates of their own file.
-> out/one_stretch_geom.gpkg (EPSG:4326) + the check QGIS's Topology Checker would make"""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
R = r"D:\Projects\Tambaram_Chengalpattu"
G = gpd.read_file("out/merged_clean.gpkg"); M = gpd.read_file(os.path.join(R, "Manual_Georeferenced", "Manual_Georeference.geojson")); n = len(G)
frozen = (G.src == "Manual_Georeference").values
g4 = list(G.to_crs(4326).geometry); g4 = [shapely.transform(g, lambda c: np.round(c, 10)) for g in g4]
bad0 = [i for i, g in enumerate(g4) if not g.is_valid]; print("plots not valid after rounding to 10 decimals:", len(bad0), [(G.src[i], G.sheet[i], G.first_read[i], shapely.is_valid_reason(g4[i])[:40]) for i in bad0][:6])
for i in bad0:
    ps = [q for q in shapely.get_parts(shapely.make_valid(g4[i])) if q.geom_type == "Polygon"]; q = max(ps, key=lambda q: q.area); assert not (G.src[i] == "Manual_Georeference") and abs(q.area - abs(g4[i].area)) < 1e-12, i; g4[i] = q
for i in np.where(frozen)[0]:
    m = M.geometry[int(G.src_row[i])]; m = m if m.geom_type == "Polygon" else shapely.get_parts(m)[0]
    assert shapely.hausdorff_distance(m, g4[i]) < 2e-10, (i, shapely.hausdorff_distance(m, g4[i])); g4[i] = m                      # the file's own coordinates
polys = lambda g: [q for q in shapely.get_parts(shapely.make_valid(g)) if q.geom_type == "Polygon" and q.area > 0]
EPS = 1.2e-9                      # 0.13 mm: a neighbour's vertex this close to an edge belongs on it


def node(ring, P, tree):
    c = np.array(ring.coords); out = []
    for a, b in zip(c[:-1], c[1:]):
        out.append(tuple(a)); L = np.hypot(*(b - a))
        if L < 1e-12: continue
        seg = shapely.LineString([a, b]); mids = []
        for k in tree.query(seg.buffer(EPS)):
            v = P[k]; da, db = np.hypot(*(v - a)), np.hypot(*(v - b))
            if da > EPS and db > EPS and seg.distance(shapely.Point(v)) <= EPS: mids.append((da, tuple(v)))
        out += [v for _d, v in sorted(set(mids))]
    out.append(tuple(c[-1])); return out


def micro():
    """crumbs under 0.01 m2 left at the ends of the closed hairlines: an overlap is taken off the plot that may move, a gap joins the plot it lies along"""
    k = 0; tr = shapely.STRtree(g4); LIM = 0.01 / 1.2e10
    for i in range(n):
        for j in tr.query(g4[i]):
            if j > i:
                x = g4[i].intersection(g4[j])
                if 0 < x.area <= LIM and not (frozen[i] and frozen[j]):
                    m, o = (j, i) if frozen[i] else (i, j) if frozen[j] else ((i, j) if g4[i].area >= g4[j].area else (j, i))
                    ps = polys(g4[m].difference(g4[o]))
                    if ps and abs(max(p.area for p in ps) - g4[m].area) <= 2 * LIM: g4[m] = max(ps, key=lambda p: p.area); k += 1
    tr = shapely.STRtree(g4); U = unary_union(g4)
    for q in polys(U):
        for r in q.interiors:
            h = shapely.Polygon(r)
            if h.area > LIM: continue
            c = sorted(((g4[t].boundary.intersection(h.boundary.buffer(1e-11)).length, t) for t in tr.query(h.buffer(1e-9)) if not frozen[t]), reverse=True)
            for _l, t in c:
                ps = polys(unary_union([g4[t], h.buffer(1e-12, join_style=2)]))
                if len(ps) == 1 and len(ps[0].interiors) == len(g4[t].interiors) and abs(ps[0].area - g4[t].area) <= 2 * LIM: g4[t] = shapely.transform(ps[0], lambda c_: np.round(c_, 10)); k += 1; break
    return k


for rnd in range(6):
    crumbs = micro(); g4[:] = [g if g.is_valid else max(polys(g), key=lambda p: p.area) for g in g4]; print("round", rnd, "| crumbs settled:", crumbs)
    P = np.unique(np.vstack([shapely.get_coordinates(g) for g in g4]), axis=0); tree = shapely.STRtree(shapely.points(P)); added = 0
    for i in range(n):
        ext = node(g4[i].exterior, P, tree); ints = [node(r, P, tree) for r in g4[i].interiors]; k = len(ext) - len(g4[i].exterior.coords) + sum(len(a) - len(r.coords) for a, r in zip(ints, g4[i].interiors))
        if k:
            g = shapely.Polygon(ext, ints)
            if g.is_valid and g.geom_type == "Polygon" and abs(g.area - g4[i].area) < 1e-13: g4[i] = g; added += k
    # vertices of two plots a hair apart (under 0.13 mm) become one: the mover is never a Manual plot
    P = np.unique(np.vstack([shapely.get_coordinates(g) for g in g4]), axis=0); tree = shapely.STRtree(shapely.points(P)); ref = {}
    fz = set(map(tuple, np.vstack([shapely.get_coordinates(g4[i]) for i in np.where(frozen)[0]]))) if frozen.any() else set()
    for k, v in enumerate(P):
        near = [j for j in tree.query(shapely.Point(v).buffer(EPS)) if j != k and np.hypot(*(P[j] - v)) <= EPS]
        if near and tuple(v) not in fz:
            c = [tuple(P[j]) for j in near]; f = [t for t in c if t in fz]; ref[tuple(v)] = f[0] if f else min(c + [tuple(v)])
    ref = {k: v for k, v in ref.items() if k != v}; snapped = 0
    if ref:
        for i in range(n):
            if frozen[i]: continue
            c = shapely.get_coordinates(g4[i]); hit = [tuple(x) in ref for x in c]
            if any(hit):
                g = shapely.transform(g4[i], lambda cc: np.array([ref.get(tuple(x), tuple(x)) for x in cc]))
                if g.is_valid and g.geom_type == "Polygon": g4[i] = g; snapped += 1
    print("round", rnd, "| vertices put on seams:", added, "| plots whose vertices were made common:", snapped)
    if not added and not snapped and not crumbs: break
for i in np.where(frozen)[0]:
    m = M.geometry[int(G.src_row[i])]; m = m if m.geom_type == "Polygon" else shapely.get_parts(m)[0]; assert shapely.hausdorff_distance(m, g4[i]) < 1e-12 and abs(m.area - g4[i].area) < 1e-18, i
assert all(g.is_valid and g.geom_type == "Polygon" for g in g4), [i for i, g in enumerate(g4) if not g.is_valid]
O = G.copy(); O["geometry"] = g4; O = O.set_crs(4326, allow_override=True); O.to_file("out/one_stretch_geom.gpkg", driver="GPKG")
# the checks of the Topology Checker: must not overlap, must not have gaps, no duplicates, valid, single part
tr = shapely.STRtree(g4); ov = []
for i in range(n):
    for j in tr.query(g4[i]):
        if j > i:
            a = g4[i].intersection(g4[j]).area
            if a > 0: ov.append(a * 1.2e10)
U = unary_union(g4); H = [shapely.Polygon(r).area * 1.2e10 for q in polys(U) for r in q.interiors]
b = [0, 1e-6, 1e-3, 0.01, 1, 1e9]; print("overlapping pairs:", len(ov), pd.cut(pd.Series(ov), b).value_counts().sort_index().to_dict()); print("gaps:", len(H), pd.cut(pd.Series(H), b).value_counts().sort_index().to_dict())
print("invalid:", sum(not g.is_valid for g in g4), "| multi-part:", sum(g.geom_type != "Polygon" for g in g4), "| duplicates:", n - len({shapely.normalize(g).wkb for g in g4}), "| fabric pieces:", len(polys(U)))
