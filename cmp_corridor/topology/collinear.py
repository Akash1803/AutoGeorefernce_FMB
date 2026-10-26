"""Drop redundant vertices (within 5 mm of the straight line through their neighbours) from every parcel of CMP_my_Part.
Shapes do not move (max 5 mm); shared edges stay where they are. -> fix_collinear.json"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace"); sys.path.insert(0, ".")
import numpy as np, geopandas as gpd, shapely
from shapely import make_valid
from validfix import valid4326
TOL = 0.005
P4 = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_my_Part.geojson"); P4["kide"] = P4.kide.astype(str); P = P4.to_crs(32644)
out = {"change": [], "add": [], "delete": []}; rep = []; tot0 = tot1 = 0
for i in P.index:
    g = make_valid(P.geometry[i]); s = g.simplify(TOL, preserve_topology=True)
    if s.geom_type != "Polygon" or not s.is_valid: continue
    n0, n1 = len(g.exterior.coords) + sum(len(r.coords) for r in g.interiors), len(s.exterior.coords) + sum(len(r.coords) for r in s.interiors)
    tot0 += n0; tot1 += n1
    if n1 >= n0: continue
    rep.append((P.kide[i], n0, n1, round(s.hausdorff_distance(g), 4), round(s.area - g.area, 3)))
    g4 = valid4326(gpd.GeoSeries([s], crs=32644).to_crs(4326).iloc[0])
    out["change"].append({"vc": str(P.village_code[i]), "old_kide": P.kide[i], "old_pt": list(make_valid(P4.geometry[i]).representative_point().coords[0]),
                          "wkt": g4.wkt, "attrs": {"area_sqm": round(s.area, 1), "area_acre": round(s.area / 4046.856, 4)}})
rep.sort(key=lambda t: -(t[1] - t[2]))
print(f"parcels simplified {len(rep)} | vertices {tot0} -> {tot1} | max shift {max(r[3] for r in rep):.4f} m | max area change {max(abs(r[4]) for r in rep):.3f} m2")
print("most cleaned:", rep[:15])
json.dump(out, open("fix_collinear.json", "w"))
