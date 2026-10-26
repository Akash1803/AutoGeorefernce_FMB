"""New thin tails in CMP_two_Parts (not present in the file Akash gave, twoparts_v0_backup): parts narrower than 2R that hang on a
parcel (long, short base). Each new tail (+ crumbs) goes WHOLE to the parcel it lines; green never changes. -> fix_two_tails.json"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd, geopandas as gpd
from shapely import make_valid
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
R = 0.75
P = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); P["kide"] = P.kide.astype(str)
O = gpd.read_file("twoparts_v0_backup.geojson"); O["kide"] = O.kide.astype(str)
OT = STRtree([make_valid(g) for g in O.geometry])
G = {i: make_valid(g) for i, g in zip(P.index, P.geometry)}; orig = dict(G)
green = {i for i in P.index if str(P.status[i]).startswith("green")}
idx = list(P.index); T = STRtree([G[i] for i in idx])
def thin_parts(g):
    op = g.buffer(-R, join_style=2).buffer(R, join_style=2)
    d = make_valid(g.difference(op)); return [p for p in getattr(d, "geoms", [d]) if p.geom_type == "Polygon" and p.area > 0.5], op
rows = []
for i in idx:
    g = G[i]
    if g.area < 20 or 2 * g.area / g.length < 3.0: continue
    parts, op = thin_parts(g)
    for p in parts:
        xs = list(p.minimum_rotated_rectangle.exterior.coords); e = sorted([((xs[k][0] - xs[k + 1][0]) ** 2 + (xs[k][1] - xs[k + 1][1]) ** 2) ** 0.5 for k in range(2)])
        if e[1] < 4 or e[1] / max(e[0], 0.1) < 4: continue
        if p.boundary.intersection(op.buffer(0.05)).length > 0.5 * e[1]: continue
        # was it already in the original file? (same parcel there covers most of the tail)
        oc = [make_valid(O.geometry.iloc[int(j)]) for j in OT.query(p) if O.kide.iloc[int(j)] == P.kide[i]]
        if oc and max(o.intersection(p).area for o in oc) > 0.6 * p.area: continue
        nb = [(p.boundary.intersection(G[idx[j]].buffer(0.05)).length, idx[j]) for j in T.query(p.buffer(0.1)) if idx[j] != i and idx[j] not in green]
        nb = [t for t in nb if t[0] > 0.5]
        if not nb or i in green: rows.append((P.kide[i], round(p.area, 1), round(e[1], 1), "left (no non-approved parcel lines it)")); continue
        _, j = max(nb)
        a = make_valid(G[i].difference(p.buffer(0.001, join_style=2))); ps = sorted([z for z in getattr(a, "geoms", [a]) if z.geom_type == "Polygon"], key=lambda z: -z.area)
        crumbs = ps[1:]
        if sum(z.area for z in crumbs) > 35: rows.append((P.kide[i], round(p.area, 1), round(e[1], 1), "left (would split the parcel)")); continue
        G[i] = ps[0]
        u = make_valid(unary_union([G[j], p.buffer(0.001, join_style=2)] + [z.buffer(0.001, join_style=2) for z in crumbs]))
        u = max([z for z in getattr(u, "geoms", [u]) if z.geom_type == "Polygon"], key=lambda z: z.area); G[j] = Polygon(u.exterior, [h for h in u.interiors if Polygon(h).area > 0.5])
        rows.append((P.kide[i], round(p.area, 1), round(e[1], 1), f"moved to {P.kide[j]}"))
print(pd.DataFrame(rows, columns=["parcel", "tail_m2", "length_m", "action"]).to_string(index=False))
chg = [i for i in idx if not G[i].equals(orig[i])]
out = {"change": [{"vc": str(P.village_code[i]), "old_kide": P.kide[i], "old_pt": list(orig[i].representative_point().coords[0]), "wkt": G[i].wkt,
                   "attrs": {"area_sqm": round(G[i].area, 1), "area_acre": round(G[i].area / 4046.856, 4)}} for i in chg], "add": [], "delete": []}
json.dump(out, open("fix_two_tails.json", "w")); print("parcels changed", len(chg))
