"""Verification of the rewritten FMB_Fine_Tunned against the saved original and the Tamilnilam download. Reads only."""
import warnings, json, os; warnings.filterwarnings("ignore")
import geopandas as gpd, pandas as pd, numpy as np, shapely
from shapely.ops import unary_union
HERE = os.path.dirname(os.path.abspath(__file__))
A = gpd.read_file(os.path.join(HERE, "out", "layer_before.geojson")); B = gpd.read_file(r"D:/Projects/Tambaram_Chengalpattu/FMB_Georeferenced_Final/FMB_Fine_Tunned.geojson"); P = pd.read_csv(os.path.join(HERE, "out", "plan.csv"), dtype=str).fillna("")
Am, Bm = A.to_crs(32644), B.to_crs(32644)
print("plots", len(A), "->", len(B), "| fields", len(A.columns) - 1, "->", len(B.columns) - 1)
print("sum of plot areas: before %.1f m2 | now %.1f m2 | difference %.2f m2" % (Am.area.sum(), Bm.area.sum(), Bm.area.sum() - Am.area.sum()))
wb = {}
for i, g in enumerate(B.geometry.normalize()): wb.setdefault(g.wkb, []).append(i)
kept = []; gone = []
for i, g in enumerate(A.geometry.normalize()):
    if g.wkb in wb and wb[g.wkb]: kept.append((i, wb[g.wkb].pop(0)))
    else: gone.append(i)
newi = sorted(j for v in wb.values() for j in v); print("old plots kept with the identical shape:", len(kept), "| old rows cut:", len(gone), "| new features:", len(newi))
tr = shapely.STRtree(list(B.geometry.iloc[newi])); worst = 0; n_bad = 0
for i in gone:
    c = [newi[k] for k in tr.query(A.geometry[i], predicate="intersects")]; c = [j for j in c if B.geometry[j].representative_point().within(A.geometry[i].buffer(1e-9))]
    d = unary_union([B.geometry[j] for j in c]).symmetric_difference(A.geometry[i]).area if c else A.geometry[i].area; worst = max(worst, d); n_bad += d > 1e-13
print("every cut row is filled exactly by its pieces: rows off = %d, largest difference %.2e sq deg (= %.4f m2)" % (n_bad, worst, worst * 1.2e10))
print("invalid:", int((~B.is_valid).sum()), "| empty:", int(B.geometry.is_empty.sum()), "| types", B.geom_type.value_counts().to_dict())
s = lambda x: x.fillna("").astype(str).str.replace(r"\.0$", "", regex=True).replace({"nan": "", "None": "", "NaT": ""})
ka = [i for i, _j in kept]; kb = [j for _i, j in kept]; untouched = ["village_code", "village", "area_acre", "area_sqm", "evidence", "geometry_source", "shape_match", "area_ratio", "gcp_points", "conform_max_m", "topo_note"]
print("kept plots, old fields that changed (only kide, survey_no, subdiv_no, plot_no may):", {c: int((s(A[c].iloc[ka]).values != s(B[c].iloc[kb]).values).sum()) for c in [c for c in A.columns if c != "geometry"] if (s(A[c].iloc[ka]).values != s(B[c].iloc[kb]).values).any()})
pk = P.iloc[ka]; Bk = B.iloc[kb]; keepk = (pk.status == "no sheet (puvi base)").values
print("kept plots: survey_no equals the sheet the plot was drawn on:", int(((Bk.survey_no.values == pk.sheet.values) | keepk).sum()), "of", len(Bk))
print("area fields off by more than 0.15 m2:", int((abs(Bm.area - B.area_sqm.astype(float)) > 0.15).sum()))
k = B.survey_no + np.where(B.subdiv_no.fillna("") != "", "/" + B.subdiv_no.fillna(""), ""); print("kide == survey_no/subdiv_no on every plot:", bool((B.kide == k).all()), "| blank kide", int((B.kide.fillna("") == "").sum()), "| unique_id built from the same fields:", bool((B.unique_id == B.village_code + "_" + B.survey_no + "/" + B.subdiv_no.fillna("").replace("", "-")).all()))
T = pd.read_excel(r"D:/Projects/Tambaram_Chengalpattu/Merged_Data/02_Tamilnilam_output.xlsx", dtype=str).fillna(""); T["u"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3) + "_" + T.fmb_survey_no + "/" + T.sub_division; T = T.set_index("u")
mm = B[B.uid_tn_avail == 1]; ok = sum(1 for r in mm.itertuples() if r.unique_id in T.index and T.loc[r.unique_id].patta_number == str(r.tamilnilam_patta_number) and T.loc[r.unique_id].owner_details == str(r.tamilnilam_owner_details) and T.loc[r.unique_id].village == r.village and T.loc[r.unique_id].tamilnilam_status == "Success" and T.loc[r.unique_id].clip_survey_no == r.survey_no and T.loc[r.unique_id].land_type_english == str(r.tamilnilam_land_type))
print("Tamilnilam: plots with data %d | identical to the download under the same district_taluk_village_survey/subdivision: %d" % (len(mm), ok))
print("plots carrying Tamilnilam data without a match:", int(((B.uid_tn_avail != 1) & (B.tamilnilam_patta_number.fillna("") != "")).sum()), "| matched plots whose number is unverified:", int((mm.subdiv_status.str.startswith(("to be checked", "not checked", "sliver"))).sum()))
print(B.subdiv_status.str.replace(r" \(.*| - readers.*| - row merges.*| - number would.*| - one reader.*", "", regex=True).value_counts().to_string())
print(B.tamilnilam_match.str.replace(r" \(.*", "", regex=True).str.replace(r"(only the parent plot|under the number).*", r"\1", regex=True).value_counts().to_string())
