"""Clean spikes / steps in every Paranur parcel (his + red; his old 70A/70B skipped) - none of them exists in the PDFs.
1. corners of all parcels closer than MERGE_M become one corner (a step between two neighbours' corners closes);
2. per ring: repeated corners and needle spikes (angle < 25 deg on an edge < 6 m) go;
3. whatever then overlaps is cut from the automation parcel (or the smaller one); holes < 15 m2 join a neighbour;
4. untouched corners keep their exact original coordinates. Writes a plan to OUT only.
Usage: python despike3.py <out.gpkg>"""
import sys, os, math
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from scipy.spatial import cKDTree
from autogeoref import ring_io
OUT = sys.argv[1]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
SKIP_HIS = {"70A", "70B"}
MERGE_M, SPIKE_DEG, SPIKE_EDGE_M, FILL_M2 = 1.0, 25.0, 6.0, 15.0
K2 = 110540.0 * 108600.0
RAIL = {"264", "265", "266", "267"}

rows = []
for L in ring_io.village_layers("Paranur"):
    n = L["name"]
    if "_parcels" not in n:
        continue
    s = n.split("_parcels")[0]
    red = n.endswith("_parcels_auto")
    if s in SKIP_HIS and not red:
        continue
    g = ring_io._read_source(L["source"])
    for pos in range(len(g)):
        if g.geometry.iat[pos] is not None and not g.geometry.iat[pos].is_empty:
            rows.append({"layer": n, "survey": s, "red": red, "pos": pos,
                         "plot_no": str(g.plot_no.iat[pos]) if "plot_no" in g else "", "geometry": g.geometry.iat[pos]})
A = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
U = A.to_crs(32644)
rings = [[tuple(p) for p in np.array(g.exterior.coords)[:-1]] for g in U.geometry]

# 1. merge corners closer than MERGE_M (connected through short edges); the corner used by most parcels wins,
#    a rail strip's or one of his corners before an automation corner
use = {}
for k, r in enumerate(rings):
    for p in r:
        use.setdefault(p, set()).add(k)
pts = list(use)
kd = cKDTree(np.array(pts))
pairs = kd.query_pairs(MERGE_M)
parent = list(range(len(pts)))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i, j in pairs:
    parent[find(i)] = find(j)
groups = {}
for i in range(len(pts)):
    groups.setdefault(find(i), []).append(i)
canon = {}
def rank(i):
    ks = use[pts[i]]
    return (len(ks), any(A.survey.iat[k] in RAIL for k in ks), any(not A.red.iat[k] for k in ks))
for members in groups.values():
    if len(members) == 1:
        continue
    best = max(members, key=rank)
    for i in members:
        canon[pts[i]] = pts[best]

def clean_ring(r):
    r = [canon.get(p, p) for p in r]
    changed = True
    while changed and len(r) > 3:
        changed = False
        for i in range(len(r)):
            a, b, d = np.array(r[i - 1]), np.array(r[i]), np.array(r[(i + 1) % len(r)])
            e1, e2 = np.linalg.norm(b - a), np.linalg.norm(d - b)
            if e1 < 1e-6:
                del r[i]; changed = True; break
            if e2 < 1e-6:
                continue
            ang = math.degrees(math.acos(max(-1, min(1, np.dot(a - b, d - b) / (e1 * e2)))))
            if ang < SPIKE_DEG and min(e1, e2) < SPIKE_EDGE_M:
                del r[i]; changed = True; break
    return r
new = []
for k, r in enumerate(rings):
    c = clean_ring(r)
    p = Polygon(c)
    if not p.is_valid:
        p2 = p.buffer(0)
        p = max(getattr(p2, "geoms", [p2]), key=lambda q: q.area) if not p2.is_empty else U.geometry.iat[k]
    new.append(p)
# back to EPSG:4326, untouched corners exact
ll = list(gpd.GeoSeries(new, crs=32644).to_crs(4326))
opts = np.vstack([shapely.get_coordinates(g_) for g_ in A.geometry])
okd = cKDTree(opts)
def back(gm):
    def f(c):
        d, k = okd.query(c); c2 = c.copy(); near = d <= 1e-8; c2[near] = opts[k[near]]; return c2
    return shapely.transform(gm, f)
geoms = [back(g_) for g_ in ll]
# 3a. overlaps created by the merge: the automation parcel gives way (else the smaller)
for _round in range(2):
    t = shapely.STRtree(geoms); a_, b_ = t.query(geoms, predicate="intersects")
    for i, j in zip(a_, b_):
        if i >= j or A.survey.iat[i] == A.survey.iat[j]:
            continue
        x = shapely.intersection(geoms[i], geoms[j])
        if sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") <= 0:
            continue
        ri, rj = A.red.iat[i], A.red.iat[j]
        loser = i if (ri and not rj) else j if (rj and not ri) else (i if geoms[i].area < geoms[j].area else j)
        other = j if loser == i else i
        d = shapely.difference(geoms[loser], geoms[other])
        geoms[loser] = max(getattr(d, "geoms", [d]), key=lambda q: q.area if q.geom_type == "Polygon" else 0)
# 3b. holes < FILL_M2 join the bordering parcel with the longest border (automation first)
for _round in range(3):
    Ua = unary_union(geoms); filled = 0
    for p in getattr(Ua, "geoms", [Ua]):
        for r in p.interiors:
            h = Polygon(r).difference(Ua)
            if h.is_empty or h.area * K2 > FILL_M2:
                continue
            cand = [(A.red.iat[i], h.boundary.intersection(geoms[i].buffer(1e-9)).length, i) for i in range(len(geoms)) if geoms[i].distance(h) < 1e-9]
            if cand:
                _r, _l, i = max(cand)
                u = unary_union([geoms[i], h])
                if u.geom_type == "Polygon":
                    geoms[i] = u; filled += 1
    if not filled:
        break
# 3c. the cuts can leave twin corners a hair apart: corners within ~1 mm become one (same in every parcel)
allp = np.vstack([shapely.get_coordinates(g_) for g_ in geoms])
uniq = np.unique(allp, axis=0)
kd2 = cKDTree(uniq)
grp = {}
for i, j in kd2.query_pairs(5e-8):
    grp.setdefault(tuple(uniq[i]), tuple(uniq[j]))
def twin(gm):
    c = [tuple(p) for p in shapely.get_coordinates(gm.exterior)][:-1]
    c = [grp.get(p, p) for p in c]
    out = [p for k, p in enumerate(c) if p != c[k - 1]]
    pg = Polygon(out)
    return pg if pg.is_valid and abs(pg.area - gm.area) < 1e-12 + 1e-6 * gm.area else gm
geoms = [twin(g_) for g_ in geoms]
B = A.copy(); B["geometry"] = geoms
# 4. report: spikes / steps left, area and shape change, topology
def flags(g):
    c = np.array(g.exterior.coords)[:-1]; m = len(c); n_ = 0
    for i in range(m):
        a, b, d = c[i - 1], c[i], c[(i + 1) % m]
        e1, e2 = np.linalg.norm(b - a), np.linalg.norm(d - b)
        if e1 < 1e-9 or e2 < 1e-9: n_ += 1; continue
        ang = math.degrees(math.acos(max(-1, min(1, np.dot(a - b, d - b) / (e1 * e2)))))
        if (ang < SPIKE_DEG and min(e1, e2) < SPIKE_EDGE_M) or (min(e1, e2) < 1.0 and ang < 120): n_ += 1
    return n_
BU = B.to_crs(32644)
rep = pd.DataFrame({"layer": A.layer, "pos": A.pos, "plot_no": A.plot_no,
                    "flags_before": [flags(g) for g in U.geometry], "flags_after": [flags(g) for g in BU.geometry],
                    "area_change_%": [round(100 * (b.area / a.area - 1), 3) for a, b in zip(U.geometry, BU.geometry)],
                    "max_move_m": [round(shapely.hausdorff_distance(a, b), 2) for a, b in zip(U.geometry, BU.geometry)],
                    "valid": [g.is_valid and g.geom_type == "Polygon" for g in B.geometry],
                    "changed": [not a.equals_exact(b, 0) for a, b in zip(A.geometry, B.geometry)]})
pd.set_option("display.width", 200)
print(rep[rep.changed | (rep.flags_after > 0)].to_string(index=False))
print("flags before %d, after %d" % (rep.flags_before.sum(), rep.flags_after.sum()))
g = B.geometry.values; t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects"); k = a_ < b_
ov = [(B.survey.iat[i], B.survey.iat[j], round(sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") * K2, 4))
      for i, j in zip(a_[k], b_[k]) if B.survey.iat[i] != B.survey.iat[j] for x in [shapely.intersection(g[i], g[j])]]
ov = [o for o in ov if o[2] > 0]
Ua = unary_union(g); holes = [round(Polygon(r).difference(Ua).area * K2, 2) for p in getattr(Ua, "geoms", [Ua]) for r in p.interiors]
print("topology after: overlaps", len(ov), ov[:6], "| holes", holes)
B.drop(columns=["red"]).to_file(OUT, layer="plan", driver="GPKG"); rep.to_csv(OUT.replace(".gpkg", "_changes.csv"), index=False)
