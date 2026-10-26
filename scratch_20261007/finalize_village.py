"""One pass over the whole Paranur group, then one complete check. Plan only (nothing of his is written here).
His layers are fixed (except tiny corner merges <= MERGE_M); red parcels give way in the order 74, 73, 72, 70B, 70A
(70A last: it closes the junction). His old 70A / 70B are replaced by the red versions.
Checks: overlaps, gaps, validity / single part, sharp corners not in the PDFs, short steps, twin corners,
PDF neighbours that do not touch, parcels touching that the PDFs do not name, shape vs PDF.
Usage: python finalize_village.py <out.gpkg>"""
import sys, os, math, json
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union, split
from scipy.spatial import cKDTree
from autogeoref import ring_io, ring_topo, printed_numbers as pn, paths
OUT = sys.argv[1]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
V = "35_04_221"; KX, KY = 108600.0, 110540.0; K2 = KX * KY
RED_ORDER = ["74", "73", "72", "70B", "70A"]
REPLACE = {"70A", "70B"}            # his old versions are replaced by the red ones
MERGE_M, FILL_M2, SEAM_W = 0.7, 25.0, 3.5

rows = []
for L in ring_io.village_layers("Paranur"):
    n = L["name"]; s = n.split("_parcels")[0]
    if "_parcels" not in n:
        continue
    red = n.endswith("_parcels_auto")
    if s in REPLACE and not red:
        continue
    g = ring_io._read_source(L["source"]).to_crs(4326)
    for pos in range(len(g)):
        rows.append({"layer": n, "survey": s, "red": red, "pos": pos, "plot_no": str(g.plot_no.iat[pos]), "geometry": g.geometry.iat[pos]})
A = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
G = list(A.geometry)
red = list(A.red)
# 1. red parcels give way in order, from existing lines only
fixed = [G[k] for k in range(len(G)) if not red[k]]
for s in RED_ORDER:
    ks = [k for k in range(len(G)) if A.survey.iat[k] == s and red[k]]
    new = ring_topo.exact_seams([G[k] for k in ks], fixed, SEAM_W)
    for k, g_ in zip(ks, new):
        G[k] = g_
    fixed += new
# 2. corners closer than MERGE_M become one corner (his / rail first, then the corner used by most plots)
use = {}
for k, g_ in enumerate(G):
    for p in map(tuple, shapely.get_coordinates(g_.exterior)[:-1]):
        use.setdefault(p, set()).add(k)
pts = list(use); arr = np.array(pts) * [KX, KY]
parent = list(range(len(pts)))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i, j in cKDTree(arr).query_pairs(MERGE_M):
    parent[find(i)] = find(j)
groups = {}
for i in range(len(pts)):
    groups.setdefault(find(i), []).append(i)
canon = {}
for m in groups.values():
    if len(m) > 1:
        best = max(m, key=lambda i: (any(not red[k] for k in use[pts[i]]), len(use[pts[i]])))
        for i in m:
            canon[pts[i]] = pts[best]
for k in range(len(G)):
    r = [canon.get(tuple(p), tuple(p)) for p in shapely.get_coordinates(G[k].exterior)[:-1]]
    r = [p for i, p in enumerate(r) if p != r[i - 1]]
    pg = Polygon(r)
    if pg.is_valid and len(r) >= 3:
        G[k] = pg

def angles(poly):
    c = np.array(poly.exterior.coords)[:-1] * [KX, KY]; m = len(c); out = []
    for i in range(m):
        a, b, d = c[i - 1], c[i], c[(i + 1) % m]
        e1, e2 = np.linalg.norm(b - a), np.linalg.norm(d - b)
        out.append((180.0 if e1 < 1e-9 or e2 < 1e-9 else math.degrees(math.acos(max(-1, min(1, np.dot(a - b, d - b) / (e1 * e2))))), e1, e2))
    return out
_sharp = {}
def pdf_sharp(s):
    if s not in _sharp:
        res = []
        for ft in json.load(open(paths.sheet_path(V, s), encoding="utf-8"))["features"]:
            res += [a for a, _e1, _e2 in angles(Polygon(np.array(ft["geometry"]["coordinates"][0]) / [KX, KY])) if a < 25]
        _sharp[s] = res
    return _sharp[s]

def needles():
    n = 0
    for k in range(len(G)):
        ang = angles(G[k]); c = [tuple(p) for p in np.array(G[k].exterior.coords)[:-1]]
        drop = [i for i, (a, _e1, _e2) in enumerate(ang) if a < 3 or (a < 10 and not any(abs(a - x) < 1.5 for x in pdf_sharp(A.survey.iat[k])))]
        if drop:
            p = Polygon([q for i, q in enumerate(c) if i not in drop])
            if p.is_valid:
                G[k] = p; n += 1
    return n

def cut_overlaps():
    n = 0
    t = shapely.STRtree(G); a_, b_ = t.query(G, predicate="intersects")
    for i, j in zip(a_, b_):
        if i >= j or A.survey.iat[i] == A.survey.iat[j]:
            continue
        x = shapely.intersection(G[i], G[j])
        if sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") <= 0:
            continue
        ri = RED_ORDER.index(A.survey.iat[i]) if red[i] else -1
        rj = RED_ORDER.index(A.survey.iat[j]) if red[j] else -1
        lose = i if ri > rj else j if rj > ri else (i if G[i].area < G[j].area else j)
        G[lose] = ring_topo._largest(shapely.difference(G[lose], G[j if lose == i else i])); n += 1
    return n

def fill_holes():
    n = 0
    U = unary_union(G)
    for p in getattr(U, "geoms", [U]):
        for r in p.interiors:
            h = Polygon(r).difference(U)
            if h.is_empty or h.geom_type != "Polygon" or h.area * K2 > FILL_M2:
                continue
            share = sorted([(h.boundary.intersection(G[k].buffer(1e-9)).length, red[k], k) for k in range(len(G)) if G[k].distance(h) < 1e-9], reverse=True)
            if not share:
                continue
            if share[0][0] >= 0.6 * h.length / 2:
                u = unary_union([G[share[0][2]], h])
                if u.geom_type == "Polygon":
                    G[share[0][2]] = u; n += 1
                continue
            if len(share) >= 2 and share[0][0] + share[1][0] >= 0.8 * h.length:
                ka, kb = share[0][2], share[1][2]
                ang = [a for a, _e1, _e2 in angles(h)]; c = np.array(h.exterior.coords)[:-1]
                tip = c[int(np.argmin(ang))]
                rr = np.array(h.minimum_rotated_rectangle.exterior.coords)[:-1]
                e = [rr[1] - rr[0], rr[2] - rr[1]]; ax = max(e, key=lambda v: np.hypot(*(v * [KX, KY]))); ax = ax / np.hypot(*ax)
                L_ = 10 * max(np.ptp(c[:, 0]), np.ptp(c[:, 1]))
                pieces = [q for q in getattr(split(h, LineString([tuple(tip - ax * L_), tuple(tip + ax * L_)])), "geoms", []) if q.area > 0]
                if len(pieces) == 2:
                    for q in pieces:
                        k = ka if q.boundary.intersection(G[ka].buffer(1e-9)).length >= q.boundary.intersection(G[kb].buffer(1e-9)).length else kb
                        u = unary_union([G[k], q])
                        if u.geom_type == "Polygon":
                            G[k] = u
                    n += 1
    return n

for _round in range(5):
    if not (needles() + cut_overlaps() + fill_holes()):
        break
for i in range(len(G)):
    others = [G[j] for j in range(len(G)) if j != i and G[j].distance(G[i]) < 1e-8]
    G[i] = ring_topo.despike(ring_topo.insert_nodes(G[i], others))
B = A.copy(); B["geometry"] = G

# ---------------------------------------------------------------- the complete check
def full_check(B):
    G = list(B.geometry); rep = {}
    t = shapely.STRtree(G); a_, b_ = t.query(G, predicate="intersects"); ov = []
    for i, j in zip(a_, b_):
        if i < j and B.survey.iat[i] != B.survey.iat[j]:
            x = shapely.intersection(G[i], G[j]); ar = sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon")
            if ar > 0: ov.append((B.survey.iat[i], B.survey.iat[j], round(ar * K2, 4)))
    rep["overlaps"] = ov
    U = unary_union(G)
    rep["gaps"] = [round(Polygon(r).difference(U).area * K2, 3) for p in getattr(U, "geoms", [U]) for r in p.interiors]
    rep["invalid_or_multipart"] = [(B.layer.iat[k], B.plot_no.iat[k]) for k in range(len(G)) if not G[k].is_valid or G[k].geom_type != "Polygon"]
    sharp, steps, twins = [], [], []
    for k in range(len(G)):
        for a, e1, e2 in angles(G[k]):
            if min(e1, e2) < 0.01: twins.append((B.layer.iat[k], B.plot_no.iat[k]))
            elif a < 25 and not any(abs(a - x) < 1.5 for x in pdf_sharp(B.survey.iat[k])): sharp.append((B.layer.iat[k], B.plot_no.iat[k], round(a, 1)))
            elif min(e1, e2) < 0.5 and a < 150: steps.append((B.layer.iat[k], B.plot_no.iat[k], round(a), round(min(e1, e2), 2)))
    rep["sharp_corners_not_in_pdf"] = sharp; rep["short_steps"] = steps; rep["twin_corners"] = twins
    table = pn.village_table(V, str(paths.PROJECT / "Auto_Georef" / V / "printed_numbers.json"))
    agreed = pn.agreed_sets(table)
    S = {s: unary_union([G[k] for k in range(len(G)) if B.survey.iat[k] == s]) for s in set(B.survey)}
    rep["pdf_neighbours_not_touching"] = sorted({tuple(sorted((a, b))) for a in S for b in agreed.get(a, set()) if b in S and S[a].distance(S[b]) > 1e-7})
    rep["touching_not_named"] = sorted({tuple(sorted((a, b))) for a in S for b in S if a < b and b not in agreed.get(a, set()) and a not in agreed.get(b, set())
                                        and S[a].boundary.intersection(S[b].buffer(1e-7)).length * 110000 > 2})
    return rep
rep = full_check(B)
for k, v in rep.items():
    print("%-28s %s" % (k, v if v else "none"))
area = pd.DataFrame({"layer": A.layer, "plot": A.plot_no,
                     "before_m2": gpd.GeoSeries(A.geometry, crs=4326).to_crs(32644).area.round(1),
                     "after_m2": gpd.GeoSeries(G, crs=4326).to_crs(32644).area.round(1)})
area["change_m2"] = (area.after_m2 - area.before_m2).round(1)
print(area[area.change_m2.abs() > 0.05].to_string(index=False))
B.drop(columns=["red"]).to_file(OUT, layer="plan", driver="GPKG")
json.dump(rep, open(OUT.replace(".gpkg", "_check.json"), "w"), indent=1, default=str)
