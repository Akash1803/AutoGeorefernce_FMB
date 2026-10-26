"""One Stretch, step 5: write Tambaram_Chengalpattu_One_Stretch.geojson (a NEW file, nothing of Akash's is overwritten),
read it back, and prove it: geometry against the three layers as received, attributes against the table, Tamilnilam
against the download. usage: python s5_write.py [OUT]"""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import mapping
from shapely.ops import unary_union
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(R, "Manual_Georeferenced", "Tambaram_Chengalpattu_One_Stretch.geojson")
G = gpd.read_file(os.path.join(HERE, "out", "one_stretch_geom.gpkg")); F = pd.read_pickle(os.path.join(HERE, "out", "one_stretch_table.pkl")); assert len(G) == len(F) and (F._i.values == np.arange(len(G))).all()
cols = [c for c in F.columns if not c.startswith("_")]
def val(c, v):
    if v is None or (isinstance(v, float) and np.isnan(v)) or v == "": return None if c in ("tamilnilam_extent_sqm",) else ("" if c not in ("area_acre", "area_sqm", "seam_moved_mm", "uid_tn_avail") else None)
    if c in ("area_acre", "area_sqm", "seam_moved_mm", "tamilnilam_extent_sqm"): return float(v)
    if c == "uid_tn_avail": return int(v)
    return str(v)
feats = [{"type": "Feature", "properties": {c: val(c, r[c]) for c in cols}, "geometry": mapping(g)} for (_, r), g in zip(F.iterrows(), G.geometry)]
fc = {"type": "FeatureCollection", "name": "Tambaram_Chengalpattu_One_Stretch", "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}}, "features": feats}
assert not os.path.exists(OUT) or "One_Stretch" in os.path.basename(OUT)
json.dump(fc, open(OUT, "w", encoding="utf-8"), ensure_ascii=False); print("written", OUT, os.path.getsize(OUT), "bytes")
# ---- read back
B = gpd.read_file(OUT); Bm = B.to_crs(32644); print("plots", len(B), "| fields", len(B.columns) - 1, "| types", B.geom_type.value_counts().to_dict(), "| invalid", int((~B.is_valid).sum()))
assert all(a.equals_exact(b, 0) for a, b in zip(B.geometry, G.geometry)), "geometry changed on the way to the file"
RAW = gpd.read_file(os.path.join(HERE, "out", "merged_raw.gpkg")); assert len(RAW) == len(B)
mv = np.array([a.hausdorff_distance(b) for a, b in zip(RAW.geometry, Bm.geometry)]); da = Bm.area.values - RAW.area.values
print("against the layers as received: plots whose shape is untouched (under 0.05 mm): %d | moved under 1 mm: %d | 1 mm to 5 cm: %d | more than 5 cm: %d | largest move %.1f mm" % ((mv < 5e-5).sum(), ((mv >= 5e-5) & (mv < 1e-3)).sum(), ((mv >= 1e-3) & (mv <= 0.05)).sum(), (mv > 0.05).sum(), 1000 * mv.max()))
print("area: as received %.1f m2 | now %.1f m2 | largest change on one plot %.2f m2" % (RAW.area.sum(), Bm.area.sum(), np.abs(da).max()))
M = gpd.read_file(os.path.join(R, "Manual_Georeferenced", "Manual_Georeference.geojson")); mb = B[B.source_layer == "Manual_Georeference"].reset_index(drop=True)
print("Manual_Georeference plots: %d | identical to their own layer, vertex for vertex: %d | kide, survey, subdivision identical: %d" % (len(mb), sum(a.equals_exact(b, 0) for a, b in zip(mb.geometry, M.geometry)), int(((mb.kide.values == M.kide.values) & (mb.survey_no.values == M.survey_no.values) & (mb.subdiv_no.fillna("").values == M.subdiv_no.fillna("").values)).sum())))
s = lambda x: x.fillna("").astype(str).str.replace(r"\.0$", "", regex=True).replace({"nan": "", "None": ""})
print("attributes that differ from the table:", {c: int((s(B[c]).values != s(F[c]).values).sum()) for c in cols if (s(B[c]).values != s(F[c]).values).any() and c not in ("area_acre", "area_sqm", "seam_moved_mm", "tamilnilam_extent_sqm")})
print("area fields off by more than 0.15 m2:", int((abs(Bm.area - B.area_sqm.astype(float)) > 0.15).sum()), "| blank kide:", int((B.kide.fillna("") == "").sum()), "| kide == survey/subdivision:", bool((B.kide == B.survey_no + np.where(B.subdiv_no.fillna("") != "", "/" + B.subdiv_no.fillna(""), "")).all()))
T = pd.read_excel(os.path.join(R, "Merged_Data", "02_Tamilnilam_output.xlsx"), dtype=str).fillna(""); T["u"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3) + "_" + T.fmb_survey_no + "/" + T.sub_division; T = T.set_index("u")
mm = B[B.uid_tn_avail == 1]; ok = sum(1 for r in mm.itertuples() if r.unique_id in T.index and T.loc[r.unique_id].patta_number == str(r.tamilnilam_patta_number) and T.loc[r.unique_id].owner_details == str(r.tamilnilam_owner_details) and T.loc[r.unique_id].village == r.village and T.loc[r.unique_id].tamilnilam_status == "Success" and T.loc[r.unique_id].clip_survey_no == r.survey_no)
print("Tamilnilam: plots with data %d | identical to the download under the same district_taluk_village_survey/subdivision: %d | data on an unmatched plot: %d" % (len(mm), ok, int(((B.uid_tn_avail != 1) & (B.tamilnilam_patta_number.fillna("") != "")).sum())))
# ---- topology, my side
g = list(B.geometry); tr = shapely.STRtree(g); ov = []
for i in range(len(g)):
    for j in tr.query(g[i]):
        if j > i:
            a = g[i].intersection(g[j]).area
            if a > 0: ov.append(a * 1.2e10)
U = unary_union(g); H = [shapely.Polygon(r).area * 1.2e10 for q in shapely.get_parts(U) for r in q.interiors]
print("TOPOLOGY (my check): overlaps %d (smallest %.3f m2) | gaps %d (smallest %.3f m2) | invalid %d | multi-part %d | duplicates %d" % (len(ov), min(ov or [0]), len(H), min(H or [0]), int((~B.is_valid).sum()), int((B.geom_type != "Polygon").sum()), len(g) - len({shapely.normalize(x).wkb for x in g})))
