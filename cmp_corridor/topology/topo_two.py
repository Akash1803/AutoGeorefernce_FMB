import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, numpy as np, pandas as pd
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
F = sys.argv[1] if len(sys.argv) > 1 else r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"
P = gpd.read_file(F)
print("rows", len(P), "| crs", P.crs, "| columns", list(P.columns)[:12], "...")
for c in ("village_code", "source", "status"):
    if c in P: print(c, P[c].astype(str).value_counts().head(8).to_dict())
print("invalid:", int((~P.is_valid).sum()), "| empty:", int(P.geometry.is_empty.sum()), "| multipart:", int((P.geom_type != "Polygon").sum()), P.geom_type.value_counts().to_dict())
from shapely import make_valid
G = [make_valid(g) for g in P.geometry]; T = STRtree(G); ov = []
for i, g in enumerate(G):
    if g is None or g.is_empty: continue
    for j in T.query(g, predicate="overlaps"):
        if j > i: ov.append(g.intersection(G[j]).area)
ov = np.array(ov)
print("overlapping pairs:", len(ov), "| < 0.01 m2:", int((ov < 0.01).sum()), "| 0.01-1:", int(((ov >= .01) & (ov < 1)).sum()), "| 1-100:", int(((ov >= 1) & (ov < 100)).sum()), "| >= 100:", int((ov >= 100).sum()), "| total %.1f m2" % ov.sum())
U = unary_union([g for g in G if g is not None])
h = np.array([Polygon(r).area for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for r in pp.interiors])
print("gaps (holes):", len(h), "| < 0.01:", int((h < .01).sum()), "| 0.01-1:", int(((h >= .01) & (h < 1)).sum()), "| 1-100:", int(((h >= 1) & (h < 100)).sum()), "| >= 100:", int((h >= 100).sum()), sorted(h.round(1))[-8:])
w = P.geometry.to_wkb(); print("exact duplicate geometries:", int(w.duplicated().sum()))
if "kide" in P:
    key = P.village_code.astype(str) + "|" + P.kide.astype(str); d = key[P.kide.astype(str).str.contains("/")].value_counts(); print("same village+kide on several rows:", int((d > 1).sum()), d[d > 1].head(10).to_dict())
