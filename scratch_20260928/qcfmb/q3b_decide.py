"""Dual QC of FMB_Fine_Tunned, step 3: one verdict per layer row.
Sources: A = the number the layer carries (first read, made at conversion); B = the any-size second read (read5) inside
the sheet polygon(s) the row was made from; C = the A-register subdivision list of that sheet number, with extents.
-> out/decide5.csv"""
import os, sys, re, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
LAYER = os.path.join(HERE, "out", "layer_before.geojson") if len(sys.argv) < 2 else sys.argv[1]
clean = lambda x: "" if x is None or str(x).strip().lower() in ("nan", "none", "null", "0", "") else str(x).strip()
G = gpd.read_file(LAYER); area = G.to_crs(32644).area; G = pd.DataFrame(G.drop(columns="geometry")).fillna(""); G["area"] = area.round(1).values
_h = pd.read_csv(os.path.join(HERE, "out", "row_home.csv"), dtype=str, index_col=0).fillna(""); G["unit"] = _h.home.values; G["sub"] = [clean(x) for x in _h["sub"].values]
K = pd.read_csv(os.path.join(HERE, "out", "link_all.csv"), dtype=str).fillna(""); K = K[K.state == "ok"]; K["share"] = K.share.astype(float); K["row"] = K.row.astype(int)
L = pd.read_csv(os.path.join(HERE, "out", "labels5_all.csv"), dtype=str).fillna("")
SH = pd.read_csv(os.path.join(HERE, "out", "sheets5_all.csv"), dtype=str).fillna("").set_index(["village_code", "sheet"]).state.to_dict()
ins = L[L.pid != "-1"]; good = ins[ins.plotno == "True"]
b_good = good.groupby(["village_code", "sheet", "pid"]).text.apply(lambda s: sorted(set(s))).to_dict()
b_any = ins.groupby(["village_code", "sheet", "pid"]).text.apply(lambda s: sorted(set(s))).to_dict()
b_near = L[(L.pid == "-1") & (L.near_pid != "-1") & (L.plotno == "True")].groupby(["village_code", "sheet", "near_pid"]).text.apply(lambda s: sorted(set(s))).to_dict()
T = pd.read_excel(os.path.join(R, "Merged_Data", "02_Tamilnilam_output.xlsx"), dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3); T["sqm"] = T.extent_hectares.astype(float) * 10000 + T.extent_ares.astype(float) * 100
T = T[T.fmb_survey_no == T.clip_survey_no]                       # only what the register holds under the sheet's own number
reg = {k: dict(zip(g.sub_division, g.sqm)) for k, g in T.groupby(["vc", "fmb_survey_no"])}
by_row = {k: g for k, g in K[K.share >= 0.5].groupby("row")}; nrows = G.groupby(["village_code", "unit"]).size().to_dict()
out = []
for i, r in G.iterrows():
    k = by_row.get(i); pids = [] if k is None else list(k.pid); l1 = [] if k is None else sorted({x for x in k.l1 if x}); key = lambda p: (r.village_code, r.unit, p)
    notunit = lambda c: [t for t in c if t != r.unit and t != re.sub(r"[A-Z]+$", "", r.unit)]          # the sheet number printed in the plot is not a subdivision
    B = notunit(sorted({t for p in pids for t in b_good.get(key(p), [])})); Ball = sorted({t for p in pids for t in b_any.get(key(p), [])}); Bn = notunit(sorted({t for p in pids for t in b_near.get(key(p), [])}))
    RG = reg.get((r.village_code, r.unit)); inreg = lambda t: RG is not None and t in RG
    fits = lambda t: inreg(t) and RG[t] > 0 and (0.5 <= r.area / RG[t] <= 2.0 or abs(r.area - RG[t]) <= 50)
    A = r["sub"]; strong = [t for t in B if len(t) >= 2 or inreg(t)]
    if r.geometry_source != "fmb sheet": st, prop = "no sheet (puvi base)", A
    elif k is None: st, prop = "not linked to a sheet polygon", A
    elif A and A in B: st, prop = "confirmed by both reads", A
    elif A and A in Ball: st, prop = "confirmed by both reads (number overlapped by other print)", A
    elif A and B: st, prop = "reads differ", "/".join(B)
    elif A and Bn and A in Bn: st, prop = "confirmed by both reads (number printed beside the plot)", A
    elif A: st, prop = "second read found nothing", A
    elif nrows[(r.village_code, r.unit)] == 1: st, prop = "whole survey", ""
    elif len(strong) == 1: st, prop = "number found by the second read", strong[0]
    elif len(strong) > 1: st, prop = "several numbers in the plot", "/".join(strong)
    elif len(Bn) == 1: st, prop = "number printed beside the plot", Bn[0]
    else: st, prop = "no number printed", "/".join(B)
    cand = [t for t in prop.split("/") if t]
    out.append({"row": i, "village_code": r.village_code, "village": r.village, "sheet": r.unit, "kide": r.kide, "survey_no": r.survey_no, "subdiv_no": A, "area": r.area, "pids": "/".join(pids), "share": round(k.share.min(), 2) if k is not None else "",
                "first_read": "/".join(l1), "second_read": "/".join(B), "second_read_all": "/".join(Ball), "beside": "/".join(Bn), "status": st, "proposed": prop, "sheet_read": SH.get((r.village_code, r.unit), ""),
                "register": "none" if RG is None else "yes", "cand_in_register": "/".join(t for t in cand if inreg(t)), "cand_area_fits": "/".join(t for t in cand if fits(t)),
                "layer_no_in_register": "" if not A else ("yes" if inreg(A) else "no"), "layer_no_area_fits": "" if not A else ("yes" if fits(A) else "no")})
D = pd.DataFrame(out); D.to_csv(os.path.join(HERE, "out", "decide5.csv"), index=False)
print(len(D), D.status.value_counts().to_dict())
print(pd.crosstab(D.village, D.status.str.replace(r" \(.*", "", regex=True)).to_string())
for s in ["second read found nothing", "reads differ", "number found by the second read", "number printed beside the plot", "several numbers in the plot"]:
    x = D[D.status == s]; print(f"{s}: {len(x)} | layer number in register {(x.layer_no_in_register == 'yes').sum()}, area fits {(x.layer_no_area_fits == 'yes').sum()} | a candidate in register {(x.cand_in_register != '').sum()}, area fits {(x.cand_area_fits != '').sum()}")
