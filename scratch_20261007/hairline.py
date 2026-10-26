"""Close hairline holes (mean width < 0.6 m) the way this morning's block fix did: the corners of the plots of one side
that lie on the hole move onto the other side's line (its vertex if within 5 cm, else a new vertex inserted into
it), shared corners move in every plot that has them, then the moved side follows the other side's line vertex by
vertex. Also removes back-tracking needles (< 3 deg) first. Usage: python hairline.py <in.gpkg> <out.gpkg>"""
import sys, os, math
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union, nearest_points
from autogeoref import ring_topo
src, OUT = sys.argv[1], sys.argv[2]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
KX, KY = 108600.0, 110540.0; K2 = KX * KY; DEG = 1 / 110000.0
B = gpd.read_file(src, layer="plan")
G = list(B.geometry)
def ring(i): return [tuple(p) for p in np.array(G[i].exterior.coords)[:-1]]
def setring(i, pts):
    out = []
    for p in pts:
        if not out or p != out[-1]: out.append(p)
    while len(out) > 3 and out[0] == out[-1]: out.pop()
    P = Polygon(out); assert P.is_valid, (B.layer.iat[i], shapely.is_valid_reason(P)); G[i] = P
def d2(p, q): return ((p[0] - q[0]) * KX) ** 2 + ((p[1] - q[1]) * KY) ** 2
def angle(a, b, c):
    a, b, c = np.array(a) * [KX, KY], np.array(b) * [KX, KY], np.array(c) * [KX, KY]
    e1, e2 = np.linalg.norm(a - b), np.linalg.norm(c - b)
    return 180.0 if e1 < 1e-9 or e2 < 1e-9 else math.degrees(math.acos(max(-1, min(1, np.dot(a - b, c - b) / (e1 * e2)))))
# 1. back-tracking needles: drop the tip and the corners that walk back along the same line
for i in range(len(G)):
    r = ring(i); changed = True
    while changed and len(r) > 4:
        changed = False
        for k in range(len(r)):
            if angle(r[k - 1], r[k], r[(k + 1) % len(r)]) < 3:
                cand = r[:k] + r[k + 1:]
                if Polygon(cand).is_valid:
                    r = cand; changed = True; break
    if r != ring(i):
        print("needle removed:", B.layer.iat[i], B.plot_no.iat[i], "area change %.2f m2" % ((Polygon(r).area - G[i].area) * K2))
        setring(i, r)
def holes():
    U = unary_union(G)
    return [Polygon(rr).difference(U) for p in getattr(U, "geoms", [U]) for rr in p.interiors]
for it in range(4):
    done = 0
    for h in holes():
        if h.is_empty or h.geom_type != "Polygon":
            continue
        width = 2 * h.area * K2 / (h.length * 110000)
        if width >= 0.6 or h.area * K2 > 200:
            continue
        touch = [k for k in range(len(G)) if G[k].distance(h) < 1e-9 and h.boundary.intersection(G[k].buffer(1e-9)).length * 110000 > 0.5]
        by_s = {}
        for k in touch: by_s.setdefault(B.survey.iat[k], []).append(k)
        # the side with the most plots on the hole moves; the rest is the reference line
        mover_s = max(by_s, key=lambda s: (len(by_s[s]), sum(G[k].area for k in by_s[s])))
        movers = by_s[mover_s]; refs = [k for k in touch if k not in movers]
        if not refs:
            continue
        RB = unary_union([G[k].exterior for k in refs]); RV = [p for k in refs for p in ring(k)]
        moves = {}
        for k in movers:
            for v in ring(k):
                if v in RV or v in moves or Point(v).distance(h) > 1e-10 or RB.distance(Point(v)) < 1e-12:
                    continue
                q = nearest_points(RB, Point(v))[0]; t = (q.x, q.y)
                near = min(RV, key=lambda w: d2(w, v))
                if d2(near, v) < 0.05 ** 2: t = near
                if d2(t, v) <= 0.8 ** 2: moves[v] = t
        for v, t in moves.items():
            for k in range(len(G)):
                if k in refs: continue
                r = ring(k)
                if v in r: r[r.index(v)] = t; setring(k, r)
        # movers follow the reference outline between corners now on it; references get the new corners
        for k in refs:
            G[k] = ring_topo.insert_nodes(G[k], [G[m] for m in movers])
        for k in movers:
            G[k] = ring_topo.insert_nodes(G[k], [G[m] for m in refs])
        done += 1
        print("hairline %.1f m2 (%.2f m wide): %s moved onto %s, %d corners" % (h.area * K2, width, mover_s, sorted({B.survey.iat[k] for k in refs}), len(moves)))
    if not done:
        break
B["geometry"] = G
g = B.geometry.values; t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects"); ov = []
for i, j in zip(a_, b_):
    if i < j and B.survey.iat[i] != B.survey.iat[j]:
        x = shapely.intersection(g[i], g[j]); ar = sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon")
        if ar > 0: ov.append((B.survey.iat[i], B.survey.iat[j], ar * K2))
U = unary_union(g); hs = [round(Polygon(rr).difference(U).area * K2, 3) for p in getattr(U, "geoms", [U]) for rr in p.interiors]
print("overlaps", ov, "| holes", hs, "| invalid", int((~B.is_valid).sum()))
B.to_file(OUT, layer="plan", driver="GPKG")
