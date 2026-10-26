"""Micro-spike + hairline clean of CMP_my_Part (all parcels, his edits included; real sharp PDF corners kept):
1. a vertex is a spike when its angle < 8 deg and the triangle it cuts off is < MAXA m2 -> vertex removed (repeat).
2. hairline overlaps (< 1 m2) trimmed from the parcel that is not green / not the bigger one.
3. hairline gaps (< 2 m2) created by 1-2 go WHOLE to the touching parcel with the longest contact.
-> fix_despike.json. Usage: python despike_all.py [MAXA=2.0]"""
import sys, json, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace"); sys.path.insert(0, ".")
import numpy as np, geopandas as gpd, shapely
from shapely import make_valid
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
from validfix import valid4326
MAXA = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
F = r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"
P4 = gpd.read_file(F); P4["kide"] = P4.kide.astype(str); P = P4
G = {i: make_valid(g) for i, g in zip(P.index, P.geometry)}; orig = dict(G)
green = {i for i in P.index if str(P.status[i]).startswith("green")}


def despike_ring(coords):
    c = [tuple(p) for p in coords[:-1]]; n0 = len(c); removed = []
    again = True
    while again and len(c) > 3:
        again = False
        for k in range(len(c)):
            a, b, d = np.array(c[k - 1]), np.array(c[k]), np.array(c[(k + 1) % len(c)])
            u, v = a - b, d - b; nu, nv = np.linalg.norm(u), np.linalg.norm(v)
            if nu < 0.01 or nv < 0.01:
                c.pop(k); again = True; break
            ang = math.degrees(math.acos(max(-1, min(1, float(u @ v) / (nu * nv)))))
            tri = 0.5 * nu * nv * math.sin(math.radians(ang))
            if ang < 8 and tri < MAXA:
                removed.append(round(tri, 2)); c.pop(k); again = True; break
    return c + [c[0]], removed


log = []
for i in P.index:
    g = G[i]
    if g.geom_type != "Polygon": continue
    ext, r1 = despike_ring(list(g.exterior.coords))
    ints = []
    for r in g.interiors:
        rr, r2 = despike_ring(list(r.coords)); r1 += r2
        if len(rr) >= 4 and Polygon(rr).area > 0.5: ints.append(rr)
    if not r1: continue
    t = Polygon(ext, ints)
    if not t.is_valid:
        v = make_valid(t); t = max([q for q in getattr(v, "geoms", [v]) if q.geom_type == "Polygon"], key=lambda q: q.area)
    G[i] = t; log.append((P.kide[i], r1))
print("spikes removed:", log)
idx = list(P.index); T = STRtree([orig[i] for i in idx])
chg = {i for i in idx if not G[i].equals(orig[i])}
# 2. hairline overlaps around the changed parcels
for i in sorted(chg):
    for j in T.query(G[i].buffer(0.1)):
        j = idx[j]
        if j == i: continue
        x = G[i].intersection(G[j]).area
        if 0 < x < 1.0:
            cut = j if (j not in green and (i in green or G[j].area < G[i].area)) else i
            other = i if cut == j else j
            d = make_valid(G[cut].difference(G[other])); G[cut] = max([q for q in getattr(d, "geoms", [d]) if q.geom_type == "Polygon"], key=lambda q: q.area)
# 3. hairline gaps next to the changed parcels
before = make_valid(unary_union([orig[i] for i in idx])); after = make_valid(unary_union([G[i] for i in idx]))
newgap = make_valid(before.difference(after))
gaps = [p for p in getattr(newgap, "geoms", [newgap]) if p.geom_type == "Polygon" and p.area > 1e-4]
for p in gaps:
    if p.area > 2: print("  gap left (> 2 m2):", round(p.area, 2)); continue
    c = [(p.boundary.intersection(G[idx[j]].buffer(0.02)).length, idx[j]) for j in T.query(p.buffer(0.05))]
    c = [t for t in c if t[0] > 0.01]
    if not c: continue
    _, j = max(c, key=lambda t: (t[1] not in green, t[0]))
    u = make_valid(unary_union([G[j], p.buffer(0.0005, join_style=2)])); G[j] = max([q for q in getattr(u, "geoms", [u]) if q.geom_type == "Polygon"], key=lambda q: q.area)
chg = [i for i in idx if not G[i].equals(orig[i])]
mx = max((abs(G[i].area - orig[i].area) for i in chg), default=0)
T2 = STRtree([G[i] for i in idx]); ov = sum(G[idx[a]].intersection(G[idx[b]]).area for a in range(len(idx)) for b in T2.query(G[idx[a]]) if b > a)
print(f"parcels changed {len(chg)} | largest area change {mx:.2f} m2 | total overlap {ov:.3f} m2 | new gaps {len(gaps)} ({sum(p.area for p in gaps):.2f} m2)")
out = {"change": [], "add": [], "delete": []}
for i in chg:
    g = G[i]; g4 = shapely.remove_repeated_points(g, 0.001); a = g.area
    out["change"].append({"vc": str(P.village_code[i]), "old_kide": P.kide[i], "old_pt": list(make_valid(P4.geometry[i]).representative_point().coords[0]),
                          "wkt": g4.wkt, "attrs": {"area_sqm": round(a, 1), "area_acre": round(a / 4046.856, 4)}})
json.dump(out, open("fix_two_despike.json", "w")); print("prepared", len(out["change"]))
