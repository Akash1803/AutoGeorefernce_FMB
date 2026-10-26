"""Corners of a plan closer than R_M (metres) become one corner, the same in every plot that has them; then the
seams get the same vertices on both sides again. Usage: python merge_close.py <in.gpkg> <out.gpkg> [R_M]"""
import sys, os, math
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from scipy.spatial import cKDTree
from autogeoref import ring_topo
src, OUT = sys.argv[1], sys.argv[2]
R_M = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
KX, KY = 108600.0, 110540.0
B = gpd.read_file(src, layer="plan")
geoms = list(B.geometry)
use = {}
for k, g in enumerate(geoms):
    for p in map(tuple, shapely.get_coordinates(g.exterior)[:-1]):
        use.setdefault(p, set()).add(k)
pts = list(use); arr = np.array(pts) * [KX, KY]
parent = list(range(len(pts)))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i, j in cKDTree(arr).query_pairs(R_M):
    parent[find(i)] = find(j)
groups = {}
for i in range(len(pts)):
    groups.setdefault(find(i), []).append(i)
canon = {}
fixed_rank = {s: r for r, s in enumerate(["267", "266", "68B", "69"])}
def rank(i):
    ks = use[pts[i]]
    return (len(ks), -min(fixed_rank.get(B.survey.iat[k], 9) for k in ks))
for m in groups.values():
    if len(m) > 1:
        best = max(m, key=rank)
        for i in m:
            canon[pts[i]] = pts[best]
out = []
for g in geoms:
    r = [canon.get(tuple(p), tuple(p)) for p in shapely.get_coordinates(g.exterior)[:-1]]
    r = [p for k, p in enumerate(r) if p != r[k - 1]]
    pg = Polygon(r)
    out.append(pg if pg.is_valid and len(r) >= 3 else g)
for i in range(len(out)):
    others = [out[j] for j in range(len(out)) if j != i and out[j].distance(out[i]) < 1e-8]
    out[i] = ring_topo.insert_nodes(out[i], others)
B["geometry"] = out
K2 = KX * KY
g = B.geometry.values; t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects")
ov = []
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
print("merged groups %d | overlaps %s | holes %s | invalid %d | spike/step flags %d" % (
    sum(1 for m in groups.values() if len(m) > 1), ov, holes, int((~B.is_valid).sum()), flags))
B.to_file(OUT, layer="plan", driver="GPKG")
