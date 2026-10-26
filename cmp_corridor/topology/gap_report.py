"""Holes in the union of a parcel layer (= QGIS Topology Checker 'gaps'), one row per hole with area, width, centroid,
touching parcels, and whether it lies on land the corridor has no vector data for. Usage: gap_report.py <geojson 32644> <out.csv> [<second geojson>]"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, numpy as np, pandas as pd, shapely
from shapely import make_valid
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
F, OUT = sys.argv[1], sys.argv[2]
P = gpd.read_file(F)
if len(sys.argv) > 3:
    Q = gpd.read_file(sys.argv[3]).to_crs(32644); Q["part"] = "B"; P["part"] = "A"
    P = pd.concat([P, Q], ignore_index=True)
else:
    P["part"] = "A"
P["kide"] = P.kide.astype(str); P["vc"] = P.village_code.astype(str)
G = [make_valid(g) for g in P.geometry]; T = STRtree(G)
U = unary_union(G)
holes = [Polygon(r) for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for r in pp.interiors]
NV = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\CMP_Roads_Buffer_Vector_Plots.gpkg", layer="buffer_no_vector_data").to_crs(32644)
BUF = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\CMP_Roads_Buffer.gpkg", layer="road_buffer_30m").to_crs(32644)
nv = unary_union(NV.geometry.values); buf = unary_union(BUF.geometry.values)
rows = []
for h in holes:
    r = shapely.maximum_inscribed_circle(h, 0.01); width = 2 * r.length
    nb = [i for i in T.query(h.buffer(0.02), predicate="intersects")]
    touch = sorted({f"{P.vc[i][-3:]}:{P.kide[i]}" + ("" if P.part[i] == "A" else "(B)") for i in nb})
    green = any(str(P.status[i]).startswith("green") for i in nb if "status" in P)
    c = h.representative_point()
    rows.append({"area_m2": round(h.area, 2), "width_m": round(width, 2), "x": round(c.x, 1), "y": round(c.y, 1), "n_touch": len(nb),
                 "parts": "".join(sorted({P.part[i] for i in nb})), "in_buffer_pct": round(100 * h.intersection(buf).area / h.area, 0) if h.area else 0,
                 "no_vector_pct": round(100 * h.intersection(nv).area / h.area, 0) if h.area else 0, "touches_green": green, "touch": "|".join(touch)[:200]})
D = pd.DataFrame(rows).sort_values("area_m2", ascending=False).reset_index(drop=True)
D.to_csv(OUT, index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300); pd.set_option("display.max_colwidth", 90)
print("holes:", len(D)); print(D.to_string())
