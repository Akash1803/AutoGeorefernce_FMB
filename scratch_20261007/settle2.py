"""Final settle of a plan, without spikes:
* needle corners go (< 3 deg always; < 10 deg when the sheet has no such corner);
* overlaps are cut from the lower-ranked plot along the other's line (red before his, then the smaller);
* a small hole (<= FILL_M2) joins the plot it runs ALONG (>= 60 % of half its perimeter); a wedge between two
  plots is split between them along its long axis through its sharpest corner, so neither gets a spike;
* seams get shared vertices. Usage: python settle2.py <in.gpkg> <out.gpkg>"""
import sys, os, math, json
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.geometry import Polygon, LineString, Point
from shapely.ops import unary_union, split
from autogeoref import ring_topo
src, OUT = sys.argv[1], sys.argv[2]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
FILL_M2, K2, KX, KY = 25.0, 110540.0 * 108600.0, 108600.0, 110540.0
B = gpd.read_file(src, layer="plan")
red = [l.endswith("_auto") for l in B.layer]
g = list(B.geometry)

def angles_m(poly):
    c = np.array(poly.exterior.coords)[:-1] * [KX, KY]; m = len(c); out = []
    for i in range(m):
        a, b, d = c[i - 1], c[i], c[(i + 1) % m]
        e1, e2 = np.linalg.norm(b - a), np.linalg.norm(d - b)
        out.append(180.0 if e1 < 1e-9 or e2 < 1e-9 else math.degrees(math.acos(max(-1, min(1, np.dot(a - b, d - b) / (e1 * e2))))))
    return out
sheet_sharp = {}
def pdf_sharp(s):
    if s not in sheet_sharp:
        f = r"D:/Projects/Tambaram_Chengalpattu/FMB_Vector/35_04_221/%s_parcels.geojson" % s
        res = []
        for ft in json.load(open(f, encoding="utf-8"))["features"]:
            res += [a for a in angles_m(Polygon(np.array(ft["geometry"]["coordinates"][0]) / [KX, KY])) if a < 25]
        sheet_sharp[s] = res
    return sheet_sharp[s]

# 1. needles
for k in range(len(g)):
    ang = angles_m(g[k]); c = [tuple(p) for p in np.array(g[k].exterior.coords)[:-1]]
    drop = [i for i, a in enumerate(ang) if a < 3 or (a < 10 and not any(abs(a - x) < 1.5 for x in pdf_sharp(B.survey.iat[k])))]
    if drop:
        p = Polygon([q for i, q in enumerate(c) if i not in drop])
        if p.is_valid:
            print("needle corner(s) removed:", B.layer.iat[k], B.plot_no.iat[k], [round(ang[i], 1) for i in drop]); g[k] = p

def cut_overlaps():
    n = 0
    t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects")
    for i, j in zip(a_, b_):
        if i >= j or B.survey.iat[i] == B.survey.iat[j]:
            continue
        x = shapely.intersection(g[i], g[j])
        if sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") <= 0:
            continue
        lose = i if (red[i] and not red[j]) else j if (red[j] and not red[i]) else (i if g[i].area < g[j].area else j)
        g[lose] = ring_topo._largest(shapely.difference(g[lose], g[j if lose == i else i])); n += 1
    return n

def fill_holes():
    n = 0
    U = unary_union(g)
    for p in getattr(U, "geoms", [U]):
        for r in p.interiors:
            h = Polygon(r).difference(U)
            if h.is_empty or h.geom_type != "Polygon" or h.area * K2 > FILL_M2:
                continue
            half = h.length / 2
            share = sorted([(h.boundary.intersection(g[k].buffer(1e-9)).length, k) for k in range(len(g)) if g[k].distance(h) < 1e-9], reverse=True)
            if not share:
                continue
            if share[0][0] >= 0.6 * half:
                u = unary_union([g[share[0][1]], h])
                if u.geom_type == "Polygon":
                    g[share[0][1]] = u; n += 1
                continue
            if len(share) >= 2 and share[0][0] + share[1][0] >= 0.8 * h.length:
                ka, kb = share[0][1], share[1][1]
                ang = angles_m(h); c = np.array(h.exterior.coords)[:-1]
                tip = c[int(np.argmin(ang))]
                rr = np.array(h.minimum_rotated_rectangle.exterior.coords)[:-1]
                e = [rr[1] - rr[0], rr[2] - rr[1]]; axis = max(e, key=lambda v: np.hypot(*(v * [KX, KY])))
                axis = axis / np.hypot(*axis)
                L = 10 * max(np.ptp(c[:, 0]), np.ptp(c[:, 1]))
                cutter = LineString([tuple(tip - axis * L), tuple(tip + axis * L)])
                pieces = [q for q in getattr(split(h, cutter), "geoms", []) if q.area > 0]
                if len(pieces) == 2:
                    for q in pieces:
                        k = ka if q.boundary.intersection(g[ka].buffer(1e-9)).length >= q.boundary.intersection(g[kb].buffer(1e-9)).length else kb
                        u = unary_union([g[k], q])
                        if u.geom_type == "Polygon":
                            g[k] = u
                    n += 1
    return n

for _round in range(4):
    if not (cut_overlaps() + fill_holes()):
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
U = unary_union(g); holes = [round(Polygon(r).difference(U).area * K2, 3) for p in getattr(U, "geoms", [U]) for r in p.interiors]
flags = []
for k in range(len(g)):
    for a in angles_m(g[k]):
        if a < 25 and not any(abs(a - x) < 1.5 for x in pdf_sharp(B.survey.iat[k])):
            flags.append((B.layer.iat[k], B.plot_no.iat[k], round(a, 1)))
print("overlaps", ov, "| holes", holes, "| invalid", int((~B.is_valid).sum()), "| sharp corners not in the PDFs:", flags)
B.to_file(OUT, layer="plan", driver="GPKG")
