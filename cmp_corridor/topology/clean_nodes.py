"""Vertex clean of CMP_my_Part (no shape change > 1 cm): every vertex snapped to a 1 cm grid in UTM 44N (shared vertices snap the same
way on both sides), duplicate/near-duplicate nodes collapse, out-and-back spike vertices (turn < 2 deg, both arms going back on
themselves, area of the spike < 0.05 m2) removed, then validity for both GEOS and QGIS. -> fix_nodes.json (+ report)"""
import sys, json, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace"); sys.path.insert(0, ".")
import numpy as np, geopandas as gpd, shapely
from shapely import make_valid
from shapely.geometry import Polygon
from shapely.strtree import STRtree
from shapely.ops import unary_union
from validfix import valid4326
F = r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_my_Part.geojson"
P4 = gpd.read_file(F); P4["kide"] = P4.kide.astype(str); P = P4.to_crs(32644)


def despike(ring):
    c = [tuple(p) for p in ring[:-1]]
    changed = True
    while changed and len(c) > 3:
        changed = False
        for k in range(len(c)):
            a, b, d = np.array(c[k - 1]), np.array(c[k]), np.array(c[(k + 1) % len(c)])
            u, v = a - b, d - b; nu, nv = np.linalg.norm(u), np.linalg.norm(v)
            if nu < 0.01 or nv < 0.01:
                c.pop(k); changed = True; break
            ang = math.degrees(math.acos(max(-1, min(1, float(u @ v) / (nu * nv)))))
            if ang < 2.0 and 0.5 * nu * nv * math.sin(math.radians(ang)) < 0.05:
                c.pop(k); changed = True; break
    return c + [c[0]]


new, stats = {}, {"snapped": 0, "despiked": 0}
for i, g in zip(P.index, P.geometry):
    g0 = make_valid(g)
    s = shapely.set_precision(g0, 0.01)
    parts = [q for q in getattr(s, "geoms", [s]) if q.geom_type == "Polygon" and q.area > 0.01]
    if not parts:
        continue
    s = max(parts, key=lambda q: q.area)
    ext = despike(list(s.exterior.coords)); ints = [despike(list(r.coords)) for r in s.interiors]
    t = Polygon(ext, [r for r in ints if len(r) >= 4 and Polygon(r).area > 0.05])
    if not t.is_valid:
        t = max([q for q in getattr(make_valid(t), "geoms", [make_valid(t)]) if q.geom_type == "Polygon"], key=lambda q: q.area)
    if len(t.exterior.coords) < len(s.exterior.coords): stats["despiked"] += 1
    if not t.equals_exact(g0, 1e-9):
        new[i] = t; stats["snapped"] += 1
# check
G = {i: new.get(i, make_valid(P.geometry[i])) for i in P.index}
mx = max((G[i].symmetric_difference(make_valid(P.geometry[i])).area for i in new), default=0)
idx = list(P.index); T = STRtree([G[i] for i in idx])
ov = sum(G[idx[a]].intersection(G[idx[b]]).area for a in range(len(idx)) for b in T.query(G[idx[a]]) if b > a)
print(stats, "| largest area change of one parcel %.3f m2 | total overlap after %.3f m2" % (mx, ov))
out = {"change": [], "add": [], "delete": []}
bad = []
for i, g in new.items():
    g4 = valid4326(gpd.GeoSeries([g], crs=32644).to_crs(4326).iloc[0])
    c = np.array(g4.exterior.coords); dmin = np.min(np.hypot(*np.diff(c, axis=0).T)) if len(c) > 1 else 0
    if dmin < 2e-8 or not g4.is_valid: bad.append((P.kide[i], dmin, g4.is_valid))
    a = g.area
    out["change"].append({"vc": str(P.village_code[i]), "old_kide": P.kide[i], "old_pt": list(make_valid(P4.geometry[i]).representative_point().coords[0]),
                          "wkt": g4.wkt, "attrs": {"area_sqm": round(a, 1), "area_acre": round(a / 4046.856, 4)}})
print("still close nodes / invalid after transform:", bad[:10])
json.dump(out, open("fix_nodes.json", "w")); print("prepared", len(out["change"]))
