"""After a corner merge: what now overlaps is cut from the lower-ranked plot along the other's line (red before his,
then the smaller), holes < FILL_M2 go to the bordering plot with the longest border (red first), seams get shared
vertices. Usage: python settle.py <in.gpkg> <out.gpkg>"""
import sys, os, math
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from autogeoref import ring_topo
src, OUT = sys.argv[1], sys.argv[2]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
FILL_M2, K2 = 25.0, 110540.0 * 108600.0
B = gpd.read_file(src, layer="plan")
red = [l.endswith("_auto") for l in B.layer]
g = list(B.geometry)
for _round in range(3):
    changed = 0
    t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects")
    for i, j in zip(a_, b_):
        if i >= j or B.survey.iat[i] == B.survey.iat[j]:
            continue
        x = shapely.intersection(g[i], g[j])
        if sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") <= 0:
            continue
        lose = i if (red[i] and not red[j]) else j if (red[j] and not red[i]) else (i if g[i].area < g[j].area else j)
        keep = j if lose == i else i
        g[lose] = ring_topo._largest(shapely.difference(g[lose], g[keep])); changed += 1
    U = unary_union(g)
    for p in getattr(U, "geoms", [U]):
        for r in p.interiors:
            h = Polygon(r).difference(U)
            if h.is_empty or h.area * K2 > FILL_M2:
                continue
            cand = [(red[k], h.boundary.intersection(g[k].buffer(1e-9)).length, k) for k in range(len(g)) if g[k].distance(h) < 1e-9]
            if cand:
                k = max(cand)[2]; u = unary_union([g[k], h])
                if u.geom_type == "Polygon":
                    g[k] = u; changed += 1
    if not changed:
        break
for i in range(len(g)):
    others = [g[j] for j in range(len(g)) if j != i and g[j].distance(g[i]) < 1e-8]
    g[i] = ring_topo.despike(ring_topo.insert_nodes(g[i], others))
B["geometry"] = g
t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects"); ov = []
for i, j in zip(a_, b_):
    if i < j and B.survey.iat[i] != B.survey.iat[j]:
        x = shapely.intersection(g[i], g[j]); ar = sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon")
        if ar > 0: ov.append((B.survey.iat[i], B.survey.iat[j], ar * K2))
U = unary_union(g); holes = [round(Polygon(r).difference(U).area * K2, 4) for p in getattr(U, "geoms", [U]) for r in p.interiors]
flags = 0
for gm in B.to_crs(32644).geometry:
    c = np.array(gm.exterior.coords)[:-1]; m = len(c)
    for i in range(m):
        a, b, d = c[i - 1], c[i], c[(i + 1) % m]
        e1, e2 = np.linalg.norm(b - a), np.linalg.norm(d - b)
        if e1 < 1e-9 or e2 < 1e-9: flags += 1; continue
        ang = math.degrees(math.acos(max(-1, min(1, np.dot(a - b, d - b) / (e1 * e2)))))
        if (ang < 25 and min(e1, e2) < 6) or (min(e1, e2) < 1.0 and ang < 120): flags += 1
print("overlaps", ov, "| holes", holes, "| invalid", int((~B.is_valid).sum()), "| flags", flags)
B.to_file(OUT, layer="plan", driver="GPKG")
