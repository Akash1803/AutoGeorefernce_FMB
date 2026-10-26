"""3VILLAGES: what the topology errors are, by kind and size. Reads only. usage: python t1_audit.py [file]"""
import sys, os, warnings, json, collections; warnings.filterwarnings("ignore")
import geopandas as gpd, pandas as pd, numpy as np, shapely
from shapely.ops import unary_union
from shapely.validation import explain_validity
p = sys.argv[1] if len(sys.argv) > 1 else r"C:/Users/FAI-Akash/Downloads/3VILLAGES.geojson"
G = gpd.read_file(p); G["fid"] = range(len(G)); M = G.to_crs(32644)
emp = M.geometry.isna() | M.geometry.is_empty; print("features", len(G), "| without geometry:", int(emp.sum()), list(G[emp][["survey_no", "plot_no", "village"]].itertuples(index=False, name=None))[:10])
M = M[~emp].copy(); inv = M[~M.is_valid]; print("invalid:", len(inv), [(r.fid, r.survey_no, r.plot_no, explain_validity(r.geometry)[:60]) for r in inv.itertuples()])
M["parts"] = [len(shapely.get_parts(g)) for g in M.geometry]; mp = M[M.parts > 1]; print("features made of several parts:", len(mp), "| of which a part under 1 m2:", sum(1 for g in mp.geometry if min(x.area for x in shapely.get_parts(g)) < 1))
M["holes"] = [sum(len(x.interiors) for x in shapely.get_parts(g)) for g in M.geometry]; print("features with holes inside:", int((M.holes > 0).sum()))
w = collections.defaultdict(list)
for r in M.itertuples(): w[shapely.normalize(r.geometry).wkb].append(r.fid)
dup = [v for v in w.values() if len(v) > 1]; print("exact duplicate shapes:", len(dup), "groups,", sum(len(v) - 1 for v in dup), "extra features")
V = M.copy(); V["geometry"] = V.geometry.buffer(0); tr = shapely.STRtree(list(V.geometry)); ids = list(V.fid); ov = []
for i, g in enumerate(V.geometry):
    for j in tr.query(g):
        if j > i:
            x = g.intersection(V.geometry.iloc[j])
            if x.area > 0: ov.append((ids[i], ids[j], x.area, x.area / min(g.area, V.geometry.iloc[j].area), max((q.length for q in shapely.get_parts(x) if q.geom_type == "Polygon"), default=0)))
O = pd.DataFrame(ov, columns=["a", "b", "sqm", "share", "len"]); O["width"] = 2 * O.sqm / O.len.clip(lower=1e-9)
print("overlapping pairs:", len(O), "| total %.1f m2" % O.sqm.sum()); b = [0, 0.001, 0.01, 0.1, 1, 10, 100, 1e9]; print("   by size (m2):", pd.cut(O.sqm, b).value_counts().sort_index().to_dict())
print("   one plot inside the other almost fully (share > 0.95):", int((O.share > 0.95).sum()), "| hairline (mean width < 2 cm):", int((O.width < 0.02).sum()), "| width 2 cm - 0.5 m:", int(((O.width >= 0.02) & (O.width < 0.5)).sum()), "| wider:", int((O.width >= 0.5).sum()))
U = unary_union(list(V.geometry)); P = list(shapely.get_parts(U)); H = [shapely.Polygon(r) for q in P for r in q.interiors]
Hd = pd.DataFrame({"sqm": [h.area for h in H], "width": [2 * h.area / h.length for h in H]}); print("fabric pieces:", len(P), "| gaps (holes of the fabric):", len(H), "| total %.1f m2" % Hd.sqm.sum()); print("   by size (m2):", pd.cut(Hd.sqm, b).value_counts().sort_index().to_dict())
print("   hairline (mean width < 2 cm):", int((Hd.width < 0.02).sum()), "| 2 cm - 0.5 m:", int(((Hd.width >= 0.02) & (Hd.width < 0.5)).sum()), "| wider:", int((Hd.width >= 0.5).sum()))
na = G[G.village.isna() & ~emp.reindex(G.index, fill_value=False)]; print("features with no attributes (village blank):", len(na), "| their area %.1f ha" % (M[M.fid.isin(na.fid)].area.sum() / 1e4))
O.to_csv("out/overlaps.csv", index=False); gpd.GeoDataFrame(Hd, geometry=H, crs=32644).to_file("out/gaps.gpkg", driver="GPKG"); print("QGIS-style error count (overlaps + gaps + duplicates + invalid + multi-part):", len(O) + len(H) + sum(len(v) - 1 for v in dup) + len(inv) + len(mp))
