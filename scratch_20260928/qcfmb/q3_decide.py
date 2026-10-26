"""Dual QC of FMB_Fine_Tunned, step 3: one verdict per layer row.
Sources: (1) the number the layer carries (first read, at conversion), (2) the second read of the PDF (read3) inside the
sheet polygon(s) the row was made from, (3) the A-register's subdivision list of that sheet number with its extent.
-> out/decide_all.csv"""
import os, sys, re, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
LAYER = os.path.join(HERE, "out", "layer_before.geojson") if len(sys.argv) < 2 else sys.argv[1]
clean = lambda x: "" if x is None or str(x).strip().lower() in ("nan", "none", "null", "0", "") else str(x).strip()
G = gpd.read_file(LAYER); area = G.to_crs(32644).area; G = pd.DataFrame(G.drop(columns="geometry")).fillna(""); G["area"] = area.round(1).values
G["unit"] = G.kide.astype(str).str.split("/").str[0]; G["sub"] = [clean(x) for x in G.subdiv_no]
K = pd.read_csv(os.path.join(HERE, "out", "link_all.csv"), dtype=str).fillna(""); K = K[K.state == "ok"]; K["share"] = K.share.astype(float); K["row"] = K.row.astype(int)
L = pd.read_csv(os.path.join(HERE, "out", "labels_all.csv"), dtype=str).fillna(""); num = L.text.str.isdigit()
L = L[(L.pid != "-1") & (((~num) & (L.isdim == "False")) | (num & (L.strict == "False") & (L.text.str.len() <= 3)))]
# a label printed as survey/subdivision, or with the sheet number in front, never reaches here (read3 joins characters, "/" is not a character)
r2 = L.groupby(["village_code", "sheet", "pid"]).text.apply(lambda s: sorted(set(s))).to_dict()
T = pd.read_excel(os.path.join(R, "Merged_Data", "02_Tamilnilam_output.xlsx"), dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3); T["sqm"] = T.extent_hectares.astype(float) * 10000 + T.extent_ares.astype(float) * 100
T = T[T.fmb_survey_no == T.clip_survey_no]                       # only what the register holds under the sheet's own number
reg = {k: dict(zip(g.sub_division, g.sqm)) for k, g in T.groupby(["vc", "fmb_survey_no"])}
by_row = {k: g for k, g in K[K.share >= 0.5].groupby("row")}; nrows = G.groupby(["village_code", "unit"]).size().to_dict()
out = []
for i, r in G.iterrows():
    k = by_row.get(i); pids = [] if k is None else list(k.pid); l1 = [] if k is None else sorted({x for x in k.l1 if x}); c = sorted({t for p in pids for t in r2.get((r.village_code, r.unit, p), [])})
    # the sheet number itself printed inside the plot is not a subdivision
    c = [t for t in c if t != r.unit and t != re.sub(r"[A-Z]+$", "", r.unit)]
    RG = reg.get((r.village_code, r.unit)); inreg = lambda t: RG is not None and t in RG
    fits = lambda t: inreg(t) and RG[t] > 0 and (0.5 <= r.area / RG[t] <= 2.0 or abs(r.area - RG[t]) <= 50)
    strong = [t for t in c if len(t) >= 2 or inreg(t)]
    if r.geometry_source != "fmb sheet": st, prop = "no sheet (puvi base)", r["sub"]
    elif k is None: st, prop = "not linked to a sheet polygon", r["sub"]
    elif nrows[(r.village_code, r.unit)] == 1 and not r["sub"]: st, prop = "whole survey", ""
    elif r["sub"] and r["sub"] in c: st, prop = "confirmed by both reads", r["sub"]
    elif r["sub"] and c: st, prop = "reads differ", "/".join(c)
    elif r["sub"]: st, prop = "second read found nothing", r["sub"]
    elif len(strong) == 1: st, prop = "number found by the second read", strong[0]
    elif len(strong) > 1: st, prop = "several numbers in the plot", "/".join(strong)
    else: st, prop = "no number printed", "/".join(c)
    out.append({"row": i, "village_code": r.village_code, "village": r.village, "sheet": r.unit, "kide": r.kide, "survey_no": r.survey_no, "subdiv_no": r["sub"], "area": r.area, "pids": "/".join(pids), "share": round(k.share.min(), 2) if k is not None else "",
                "first_read": "/".join(l1), "second_read": "/".join(c), "status": st, "proposed": prop, "register": "no register for this sheet" if RG is None else "/".join(f"{t}:{'in register' + (', area fits' if fits(t) else ', area differs') if inreg(t) else 'not in register'}" for t in (prop.split("/") if prop else [])),
                "reg_has_layer_no": "" if not r["sub"] else ("yes" if inreg(r["sub"]) else "no"), "reg_area_fits_layer_no": "" if not r["sub"] else ("yes" if fits(r["sub"]) else "no")})
D = pd.DataFrame(out); D.to_csv(os.path.join(HERE, "out", "decide_all.csv"), index=False)
print(len(D), D.status.value_counts().to_dict())
print(pd.crosstab(D.village, D.status).to_string())
x = D[D.status == "confirmed by both reads"]; print("confirmed by both reads: in the register %d, area fits %d, of %d" % ((x.reg_has_layer_no == "yes").sum(), (x.reg_area_fits_layer_no == "yes").sum(), len(x)))
x = D[D.status == "second read found nothing"]; print("second read found nothing: in the register %d, area fits %d, of %d" % ((x.reg_has_layer_no == "yes").sum(), (x.reg_area_fits_layer_no == "yes").sum(), len(x)))
x = D[D.status == "reads differ"]; print("reads differ: layer number in register %d of %d" % ((x.reg_has_layer_no == "yes").sum(), len(x)))
x = D[D.status == "number found by the second read"]; print("found by the second read: in register with fitting area %d, in register %d, of %d" % (x.register.str.contains("area fits").sum(), x.register.str.contains("in register").sum(), len(x)))
