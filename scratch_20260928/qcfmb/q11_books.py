"""The two workbooks, written from the layers as they now stand on disk: the dual QC record beside FMB_Fine_Tunned, and
the master Excel (Akash's format: fmb + Tamilnilam sheets with Unique ID 1 / Unique ID 2) in Merged_Data."""
import os, sys, json, re, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
A = gpd.read_file(os.path.join(HERE, "out", "layer_before.geojson")); B = gpd.read_file(os.path.join(R, "FMB_Georeferenced_Final", "FMB_Fine_Tunned.geojson")); M = gpd.read_file(os.path.join(R, "Manual_Georeferenced", "Manual_Georeference.geojson"), ignore_geometry=True).fillna("")
wb = {}
for i, g in enumerate(B.geometry.normalize()): wb.setdefault(g.wkb, []).append(i)
kept = [(i, wb[g.wkb].pop(0)) for i, g in enumerate(A.geometry.normalize()) if g.wkb in wb and wb[g.wkb]]
A = pd.DataFrame(A.drop(columns="geometry")).fillna(""); F = pd.DataFrame(B.drop(columns="geometry")).fillna("")
grp = lambda s: s.str.replace(r" \(.*", "", regex=True).str.replace(r"(to be checked by hand).*", r"\1", regex=True)
def fit(w, name, d):
    ws = w.sheets[name]; ws.freeze_panes(1, 0)
    for i, c in enumerate(d.columns): ws.set_column(i, i, min(max(len(str(c)), int(d[c].astype(str).str.len().quantile(0.95)) if len(d) else 0) + 2, 70))
S1 = F.groupby(grp(F.subdiv_status)).size().rename("Plots").reset_index().rename(columns={"subdiv_status": "Sub Div No - how it was verified"}).sort_values("Plots", ascending=False); S1.loc[len(S1)] = ["Total", len(F)]
S2 = pd.crosstab(F.village, grp(F.subdiv_status)); S2["Total"] = S2.sum(axis=1); S2 = S2.reset_index().rename(columns={"village": "Village"})
ka = [i for i, _j in kept]; kb = [j for _i, j in kept]; o = A.iloc[ka].reset_index(drop=True); n = F.iloc[kb].reset_index(drop=True)
ch = (o.kide != n.kide) | (o.survey_no != n.survey_no) | (o.subdiv_no.astype(str) != n.subdiv_no.astype(str))
C = pd.DataFrame({"Village Code": n.village_code, "Village": n.village, "Kide (before)": o.kide, "Kide (now)": n.kide, "Survey No (before)": o.survey_no, "Survey No (now)": n.survey_no, "Sub Div No (before)": o.subdiv_no, "Sub Div No (now)": n.subdiv_no, "Area (acre)": n.area_acre, "Verified by": n.subdiv_status})[ch].sort_values(["Village Code", "Survey No (now)", "Kide (now)"])
cut = F[F.topo_note.astype(str).str.contains("cut from the merged row")]; SP = pd.DataFrame({"Village Code": cut.village_code, "Village": cut.village, "Row before (several plots under one number)": cut.topo_note.astype(str).str.extract(r"cut from the merged row (.*)$")[0], "Kide (now)": cut.kide, "Survey No": cut.survey_no, "Sub Div No": cut.subdiv_no,
                                                                                               "Area (acre)": cut.area_acre, "Area (sq m)": cut.area_sqm, "Verified by": cut.subdiv_status}).sort_values(["Village Code", "Row before (several plots under one number)", "Kide (now)"])
H = F[F.subdiv_status.str.startswith(("to be checked", "not checked", "sliver"))]; H = pd.DataFrame({"Village Code": H.village_code, "Village": H.village, "Kide": H.kide, "Survey No": H.survey_no, "Sub Div No (unverified)": H.subdiv_no, "Area (acre)": H.area_acre, "Area (sq m)": H.area_sqm, "Why": H.subdiv_status}).sort_values(["Village Code", "Survey No", "Kide"])
out = os.path.join(R, "FMB_Georeferenced_Final", "FMB_Fine_Tunned_Dual_QC.xlsx"); books = (("Summary", S1), ("By village", S2), ("Numbers corrected", C), ("Rows cut into their plots", SP), ("Still to check by hand", H))
with pd.ExcelWriter(out, engine="xlsxwriter") as w:
    for nm, d in books: d.to_excel(w, sheet_name=nm, index=False); fit(w, nm, d)
print("QC workbook:", {nm: len(d) for nm, d in books})
def fmb(D, layer):
    d = pd.DataFrame({"Layer": layer, "KIDE": D.kide, "Area": D.area_acre, "DC": D.village_code.str.split("_").str[0], "TC": D.village_code.str.split("_").str[1], "VC": D.village_code.str.split("_").str[2], "District": D.district, "Taluk": D.taluk, "Village": D.village,
                      "Survey Number": D.survey_no, "Subdivision Number": [x if x else "-" for x in D.subdiv_no.astype(str)]})
    d["Unique ID 1"] = d.DC + "_" + d.TC + "_" + d.VC + "_" + d["Survey Number"]; d["Unique ID 2"] = d["Unique ID 1"] + "/" + d["Subdivision Number"]; assert (d["Unique ID 2"].values == D.unique_id.values).all()
    d["Sub Div No verified by"] = D.subdiv_status.values; d["Tamilnilam match"] = D.tamilnilam_match.values; d["uid_tn_avail"] = D.uid_tn_avail.values
    for c in [c for c in D.columns if c.startswith("tamilnilam_") and c != "tamilnilam_match"]: d[c] = D[c].values
    return d
FM = pd.concat([fmb(F, "FMB_Fine_Tunned"), fmb(M, "Manual_Georeference")], ignore_index=True)
T = pd.read_excel(os.path.join(R, "Merged_Data", "02_Tamilnilam_output.xlsx"), dtype=str).fillna("")
TN = pd.DataFrame({"survey no": T.fmb_survey_no, "subdivision": T.sub_division, "district id": T.district_id.str.zfill(2), "taluk id": T.taluk_id.str.zfill(2), "village id": T.village_id.str.zfill(3), "response status": T.tamilnilam_status, "village": T.village, "asked in CLIP as": T.clip_survey_no})
for c in ["patta_number", "government_priority", "extent_hectares", "extent_ares", "land_type_english", "poramboke", "cultivable", "owner_details", "owner_update_dates", "owner_number", "remarks", "remarks1"]: TN[c] = T[c]
TN["Unique ID 1"] = TN["district id"] + "_" + TN["taluk id"] + "_" + TN["village id"] + "_" + TN["survey no"]; TN["Unique ID 2"] = TN["Unique ID 1"] + "/" + TN["subdivision"]
for L in ("FMB_Fine_Tunned", "Manual_Georeference"):
    used = set(FM[(FM.uid_tn_avail.astype(str).str.replace(".0", "", regex=False) == "1") & (FM.Layer == L)]["Unique ID 2"]); TN["in " + L] = [1 if u in used else 0 for u in TN["Unique ID 2"]]
g2 = FM["Tamilnilam match"].str.replace(r" \(.*", "", regex=True).str.replace(r"(only the parent plot|under the number).*", r"\1", regex=True)
S = pd.crosstab([FM.Layer, FM.Village], g2); S["Total"] = S.sum(axis=1); S = S.reset_index()
NA = FM[FM.uid_tn_avail.astype(str).str.replace(".0", "", regex=False) != "1"][["Layer", "Village", "KIDE", "Survey Number", "Subdivision Number", "Area", "Unique ID 2", "Sub Div No verified by", "Tamilnilam match"]]
mo = os.path.join(R, "Merged_Data", "05_Master_Excel.xlsx")
with pd.ExcelWriter(mo, engine="xlsxwriter") as w:
    for nm, d in (("fmb", FM), ("Tamilnilam", TN), ("Match summary", S), ("Not appended", NA)): d.to_excel(w, sheet_name=nm, index=False); fit(w, nm, d)
    pd.DataFrame({"Note": ["Reginet: not downloaded yet (2026-09-28)"]}).to_excel(w, sheet_name="Reginet", index=False)
print("master Excel:", {"fmb": len(FM), "Tamilnilam": len(TN), "not appended": len(NA), "register rows used by FMB_Fine_Tunned": int(TN["in FMB_Fine_Tunned"].sum()), "by Manual_Georeference": int(TN["in Manual_Georeference"].sum())})
print(S.groupby("Layer").sum(numeric_only=True).T.to_string())
