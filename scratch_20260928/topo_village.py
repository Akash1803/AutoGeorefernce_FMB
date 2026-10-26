"""Whole-village topology on the layers QGIS shows (+ named extra files): overlaps between surveys and between plots,
enclosed holes, invalid / multipart plots, duplicate plot numbers. Read-only."""
import warnings, sys; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
from autogeoref import visible, paths, fabric_check as fc
V = sys.argv[1]; extra = sys.argv[2:]; D = paths.vector_dir(V); src = dict(visible.shown(V))
for s in extra: src.setdefault(s, (D / f"{s}_parcels_auto.gpkg", "parcels"))
parts = []
for s, (f, l) in src.items():
    g = gpd.read_file(f, layer=l); g = g.to_crs(32644) if g.crs and g.crs.to_epsg() != 32644 else g
    g["s"] = s; g["raw_valid"] = g.geometry.is_valid; g["parts"] = [len(fc._polys(x)) for x in g.geometry]; parts.append(g[["s", "plot_no", "raw_valid", "parts", "geometry"]])
P = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=32644)
print(f"{V}: {len(src)} layers, {len(P)} plots, {P.area.sum() / 1e4:.2f} ha")
print("invalid plots:", [(r.s, r.plot_no) for _, r in P.iterrows() if not r.raw_valid], "| multipart plots:", [(r.s, r.plot_no, r.parts) for _, r in P.iterrows() if r.parts != 1])
P["geometry"] = P.geometry.buffer(0); sidx = P.sindex; ov = []
for i, g in enumerate(P.geometry):
    for j in sidx.query(g):
        if j > i:
            a = g.intersection(P.geometry[j]).area
            if a > 0.01: ov.append((round(a, 3), f"{P.s[i]}/{P.plot_no[i]}", f"{P.s[j]}/{P.plot_no[j]}"))
print("plot overlaps > 0.01 m2:", sorted(ov, reverse=True))
L = {s: unary_union(list(g.geometry)) for s, g in P.groupby("s")}; geo = unary_union(list(L.values()))
holes = sorted([shapely.Polygon(r) for q in fc._polys(geo) for r in q.interiors if shapely.Polygon(r).area > 0.01], key=lambda h: -h.area)
print("fabric pieces:", len(fc._polys(geo)), "| enclosed holes > 0.01 m2:", len(holes))
for h in holes:
    nb = sorted(s for s, g in L.items() if g.distance(h) < 0.01); w = 2 * h.area / h.length
    print(f"   {h.area:8.3f} m2 at ({h.centroid.x:.1f},{h.centroid.y:.1f}) mean width {1000 * w:6.1f} mm, length {h.length / 2:5.1f} m, between {nb}")
dup = {s: [k for k, n in g.plot_no.dropna().astype(str).value_counts().items() if n > 1] for s, g in P.groupby("s")}
print("duplicate plot numbers:", {s: d for s, d in dup.items() if d}, "| plots without a number:", {s: int(g.plot_no.isna().sum()) for s, g in P.groupby("s") if g.plot_no.isna().sum()})
