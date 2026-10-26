"""Resolve the real overlaps of CMP_two_Parts (EPSG:32644). For each overlapping pair the overlap stays with ONE parcel:
 1. green (approved) > red georeferenced PDF sheet > team parcel
 2. two of the same rank: the assignment that leaves both closer to their Tamilnilam extent (tamilnilam_extent_sqm)
 3. no register area: the bigger parcel gives up the overlap (smallest relative change)
The loser is CUT only (nothing is added anywhere). -> fix_two_overlaps.json + two_overlaps_report.csv"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd, geopandas as gpd
from shapely import make_valid
from shapely.strtree import STRtree
import os
A = gpd.read_file(os.environ.get("CF_FILE", r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson")); PFX = os.environ.get("CF_PFX", "two")
A["kide"] = A.kide.astype(str)
G = {i: make_valid(g) for i, g in zip(A.index, A.geometry)}; orig = dict(G)
rank = {i: (2 if str(A.status[i]).startswith("green") else 1 if str(A.source[i]).startswith("georef") else 0) for i in A.index}
reg = {i: (float(A.tamilnilam_extent_sqm[i]) if str(A.tamilnilam_extent_sqm[i]) not in ("", "None", "nan") else None) for i in A.index}
poly = lambda g: g if g.geom_type == "Polygon" else max([z for z in getattr(g, "geoms", [g]) if z.geom_type == "Polygon"], key=lambda z: z.area)
idx = list(A.index); T = STRtree([G[i] for i in idx])
pairs = sorted({(min(idx[a], idx[b]), max(idx[a], idx[b])) for a in range(len(idx)) for b in T.query(G[idx[a]], predicate="overlaps") if idx[b] != idx[a]})
rep = []
for a, b in pairs:
    ov = G[a].intersection(G[b]).area
    if ov < 0.01: continue
    if rank[a] != rank[b]: win, why = (a, b)[rank[b] > rank[a]], "approved / rigid sheet keeps it"
    elif reg[a] and reg[b]:
        ea = abs(G[a].area - reg[a]) + abs(G[b].area - ov - reg[b]); eb = abs(G[a].area - ov - reg[a]) + abs(G[b].area - reg[b])
        win, why = (a if ea <= eb else b), "closer to the register extents"
    else: win, why = (a if G[a].area >= G[b].area else b), "no register: the bigger parcel gives it up"
    if why.startswith("no register"): win = b if win == a else a      # the bigger one LOSES
    lose = b if win == a else a
    cut = make_valid(G[lose].difference(G[win])); parts = sorted([z for z in getattr(cut, "geoms", [cut]) if z.geom_type == "Polygon"], key=lambda z: -z.area)
    if not parts: rep.append((A.kide[a], A.kide[b], round(ov, 1), "SKIPPED: loser would vanish")); continue
    G[lose] = parts[0]
    rep.append((A.kide[a], A.kide[b], round(ov, 1), f"{A.kide[win]} keeps ({why}); {A.kide[lose]} cut" + (f"; crumbs {sum(z.area for z in parts[1:]):.1f} m2" if len(parts) > 1 else "")))
R = pd.DataFrame(rep, columns=["parcel_a", "parcel_b", "overlap_m2", "decision"]); R.to_csv(f"{PFX}_overlaps_report.csv", index=False, encoding="utf-8-sig")
print(R.to_string(index=False))
# checks: remaining overlaps, thin parts created on the cut parcels
T2 = STRtree([G[i] for i in idx]); left = sum(G[idx[x]].intersection(G[idx[y]]).area for x in range(len(idx)) for y in T2.query(G[idx[x]], predicate="overlaps") if y > x)
chg = [i for i in idx if not G[i].equals(orig[i])]
thin = [(A.kide[i], round(G[i].difference(G[i].buffer(-0.25, join_style=2).buffer(0.25, join_style=2)).area - orig[i].difference(orig[i].buffer(-0.25, join_style=2).buffer(0.25, join_style=2)).area, 2)) for i in chg]
print("overlap left %.2f m2 | parcels cut %d | new thin parts (> 0.2 m2): %s" % (left, len(chg), [t for t in thin if t[1] > 0.2]))
out = {"change": [{"vc": str(A.village_code[i]), "old_kide": A.kide[i], "old_pt": list(orig[i].representative_point().coords[0]), "wkt": G[i].wkt,
                   "attrs": {"area_sqm": round(G[i].area, 1), "area_acre": round(G[i].area / 4046.856, 4)}} for i in chg], "add": [], "delete": []}
json.dump(out, open(f"fix_{PFX}_overlaps.json", "w"))
gpd.GeoDataFrame(A.drop(columns="geometry"), geometry=[G[i] for i in A.index], crs=A.crs).to_file(f"{PFX}_overlaps_preview.gpkg", driver="GPKG")
