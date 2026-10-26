"""Strict topology count on the stored EPSG:4326 coordinates (as QGIS Topology Checker sees them, zero tolerance)."""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, numpy as np
from shapely import make_valid
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
F = sys.argv[1] if len(sys.argv) > 1 else r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_my_Part.geojson"
P = gpd.read_file(F); P["kide"] = P.kide.astype(str)
G = list(P.geometry); T = STRtree(G)
k = 1 / 8.1e-11   # deg^2 -> ~m2 near 12.8 N
ov = []
for i, g in enumerate(G):
    for j in T.query(g, predicate="overlaps"):
        if j > i:
            a = g.intersection(G[j]).area * k
            ov.append((P.kide[i], P.kide[j], a))
ov = np.array([o[2] for o in ov]) if ov else np.array([])
print("rows", len(P), "| invalid:", int((~P.is_valid).sum()), "| multipart:", int((P.geom_type != "Polygon").sum()))
print("overlapping pairs:", len(ov), "| of them < 0.01 m2:", int((ov < 0.01).sum()), "| 0.01-1 m2:", int(((ov >= 0.01) & (ov < 1)).sum()), "| >= 1 m2:", int((ov >= 1).sum()))
U = unary_union(G)
holes = [Polygon(r).area * k for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for r in pp.interiors]
h = np.array(holes)
print("gaps (holes in the union):", len(h), "| < 0.01 m2:", int((h < 0.01).sum()), "| 0.01-1:", int(((h >= 0.01) & (h < 1)).sum()), "| 1-100:", int(((h >= 1) & (h < 100)).sum()), "| >= 100:", int((h >= 100).sum()), [round(x) for x in sorted(h)[-6:]])
dup = P.geometry.to_wkb().duplicated().sum(); print("exact duplicate geometries:", int(dup))
