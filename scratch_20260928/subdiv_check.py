"""Subdivision_Check.xlsx beside the Manual_Georeference layer: what every Sub Div No rests on (FMB sheet / A-register).
Reads the layer and the cross-check workbook, writes nothing else."""
import warnings, os, shutil, datetime, json; warnings.filterwarnings("ignore")
import pandas as pd, geopandas as gpd
P = r"D:/Projects/Tambaram_Chengalpattu"; OUT = P + "/Manual_Georeferenced/Subdivision_Check.xlsx"
assert not os.path.exists(P + "/Manual_Georeferenced/~$Subdivision_Check.xlsx"), "the workbook is open in Excel"
G = gpd.read_file(P + "/Manual_Georeferenced/Manual_Georeference.geojson", ignore_geometry=True).fillna("")
X = pd.read_excel(P + "/Merged_Data/03_Crosscheck_Manual_Georeference.xlsx", sheet_name="Plots", dtype=str).fillna(""); X = X[X.kide != ""].reset_index(drop=True)
assert len(X) == len(G)
G["k"] = G.village_code + "|" + G.kide + "|" + G.area_sqm.astype(float).round(0).astype(int).astype(str)
T = pd.read_excel(P + "/Merged_Data/02_Tamilnilam_output.xlsx", dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3)
T = T.drop_duplicates(["vc", "fmb_survey_no", "sub_division"]).set_index(["vc", "fmb_survey_no", "sub_division"])
# the cross-check rows follow the layer grouped by village and survey: join on village, kide and its own measured area, one to one
X["k0"] = X.village_code + "|" + X.kide; G["k0"] = G.village_code + "|" + G.kide
res, reg, land, patta = [], [], [], []
pool = {k: g.to_dict("records") for k, g in X.groupby("k0")}
for r in G.itertuples():
    c = pool[r.k0]; j = min(range(len(c)), key=lambda i: abs(float(c[i]["our_area_sqm"]) - float(r.area_sqm))); m = c.pop(j)
    assert abs(float(m["our_area_sqm"]) - float(r.area_sqm)) < 1.0, (r.kide, m["our_area_sqm"], r.area_sqm)
    res.append(m["result"]); reg.append(m["register_area_sqm"])
    t = T.loc[(r.village_code, r.survey_no, r.subdiv_no or "-")] if m["result"] == "match" and (r.village_code, r.survey_no, r.subdiv_no or "-") in T.index else None
    land.append("" if t is None else t.land_type); patta.append("" if t is None else t.patta_number)
G["A-register"] = res; G["A-register area (sq m)"] = reg; G["Land type (A-register)"] = land; G["Patta No"] = patta
G["buf"] = G.in_rail_buffer_30m.astype(str).str.lower().isin(["true", "1", "yes"])
ren = {"village": "Village", "kide": "Kide", "survey_no": "Survey No", "subdiv_no": "Sub Div No", "area_acre": "Area (acre)", "area_sqm": "Area (sq m)", "subdiv_status": "Sub Div No - source"}
S = G.groupby(["subdiv_status", "buf"]).size().unstack(fill_value=0).reindex(columns=[True, False], fill_value=0); S.columns = ["In 30 m rail buffer", "Outside buffer"]; S["Total"] = S.sum(axis=1)
S.loc["Total"] = S.sum(); S = S.reset_index().rename(columns={"subdiv_status": "Sub Div No - source"})
cols = list(ren.values()) + ["A-register", "A-register area (sq m)", "Land type (A-register)", "Patta No"]
B = G[G.buf].rename(columns=ren).sort_values(["Village", "Survey No", "Sub Div No"])[cols]
N = G[(G.subdiv_no == "") & ~G.subdiv_status.str.startswith("whole")].rename(columns=ren).assign(**{"In 30 m rail buffer": lambda d: d.buf.map({True: "yes", False: ""})}).sort_values(["In 30 m rail buffer", "Village", "Survey No"], ascending=[False, True, True])[list(ren.values()) + ["A-register", "In 30 m rail buffer"]]
old = pd.read_excel(OUT, sheet_name="Numbers corrected", dtype=str).fillna("")
add = pd.DataFrame([
    ("Thirukatchur", "42", "6", "3A2", "printed 3A2 (the 6 read earlier is a digit of a dimension); the A-register has 42/3A2"),
    ("Thirukatchur", "89A", "6", "2", "printed 2 (the 6 read earlier is a digit of a dimension); the A-register has 89A/2"),
    ("Thirukatchur", "526A", "", "1J2", "1J2 is printed on the seam of three plots; the A-register gives 1J2 0.05 ares, which is this plot"),
    ("Thirukatchur", "192", "", "2", "no number printed; taken from the A-register: 192/2 is 5 ha 40.50 ares of Poramboke, which is this block"),
    ("Thirukatchur", "39B", "", "1", "no number printed; 1 is the only number missing on the sheet and the A-register has 39B/1 (1.00 are)"),
    ("Thirukatchur", "92A", "", "4C", "no number printed; set by Akash on 2026-09-28 (the A-register has 92A/4C, 2.50 ares)")], columns=old.columns)
key = lambda d: d["Village"] + "|" + d["Survey No"] + "|" + d["Sub Div No (now)"]
C = pd.concat([old[~key(old).isin(key(add))], add], ignore_index=True).sort_values(["Village", "Survey No", "Sub Div No (now)"])
# every corrected number must be in the layer today
have = set(G.village + "|" + G.survey_no + "|" + G.subdiv_no); gone = C[~(key(C)).isin(have) & (C["Sub Div No (now)"] != "")]; assert gone.empty, gone
if not os.environ.get("NO_BACKUP"): shutil.copy2(OUT, P + "/_logs/Subdivision_Check_before_%s.xlsx" % datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
with pd.ExcelWriter(OUT, engine="xlsxwriter") as w:
    S.to_excel(w, sheet_name="Summary", index=False); B.to_excel(w, sheet_name="Buffer zone plots", index=False); C.to_excel(w, sheet_name="Numbers corrected", index=False); N.to_excel(w, sheet_name="No number on sheet", index=False)
    for n, d in (("Summary", S), ("Buffer zone plots", B), ("Numbers corrected", C), ("No number on sheet", N)):
        for i, c in enumerate(d.columns): w.sheets[n].set_column(i, i, min(max(len(str(c)), int(d[c].astype(str).str.len().quantile(0.95))) + 2, 70))
        w.sheets[n].freeze_panes(1, 0)
print(S.to_string()); print("buffer plots", len(B), B["A-register"].str.split(":").str[0].str.replace(r"parent plot .*", "parent plot", regex=True).value_counts().to_dict()); print("corrected", len(C), "| no number", len(N), N["In 30 m rail buffer"].value_counts().to_dict())
print(B[B["A-register"] != "match"].to_string())
