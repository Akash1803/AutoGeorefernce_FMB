"""Dual QC, step 8: the corrected attribute table of FMB_Fine_Tunned (and the pieces of the split rows), then the
strict Tamilnilam match. Nothing is written to the layer here. -> out/final_rows.pkl, out/final_summary.txt"""
import os, re, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
TALUK = {"04": "Chengalpattu", "05": "Tambaram", "15": "Vandalur"}                       # as printed in the sheet headers
G = gpd.read_file(os.path.join(HERE, "out", "layer_before.geojson")); P = pd.read_csv(os.path.join(HERE, "out", "plan.csv"), dtype=str).fillna("")
assert len(G) == len(P) and (G.kide.fillna("").values == P.kide.values).all()
S = pd.read_csv(os.path.join(HERE, "out", "split_pieces.csv"), dtype=str).fillna("")
V1 = pd.read_csv(os.path.join(HERE, "out", "visual_read.csv"), dtype=str).fillna("") if os.path.exists(os.path.join(HERE, "out", "visual_read.csv")) else pd.DataFrame(columns=["row", "final", "how"])
V2 = pd.read_csv(os.path.join(HERE, "out", "visual_read2.csv"), dtype=str).fillna("") if os.path.exists(os.path.join(HERE, "out", "visual_read2.csv")) else pd.DataFrame(columns=["row", "final", "how", "settled"])
v1 = {r.row: (r.final, r.how) for r in V1.itertuples()}; v2 = {r.row: (r.final, r.how) for r in V2.itertuples() if str(r.settled) == "True"}
def status(src, n):
    if src.startswith("both reads agree"): return "printed on the FMB sheet - both reads agree"
    if src.startswith("second read (printed beside"): return "printed beside the plot on the FMB sheet - second read, A-register agrees"
    if src.startswith("second read + A-register"): return "printed on the FMB sheet - second read, A-register agrees"
    if src.startswith("first read + A-register"): return "printed on the FMB sheet - first read, A-register agrees"
    if src.startswith("two readers: no number"): return "no number printed on the FMB sheet - checked on the PDF by two readers"
    if src.startswith("two readers agree") or src.startswith("one reader and the second read"): return "printed on the FMB sheet - checked on the PDF (%s)" % src.split(" (")[0]
    if src.startswith("whole survey"): return "whole survey - the FMB sheet has no subdivision"
    if src.startswith("no number printed"): return "no number printed on the FMB sheet"
    if src.startswith("no FMB sheet"): return "no FMB sheet - Puvi parcel"
    if src.startswith("not checked"): return "not checked - the plot could not be tied to a sheet polygon"
    return "to be checked by hand - " + src
rows = []; split_rows = set(S.row.astype(int))
for i, p in P.iterrows():
    base = {c: G.at[i, c] for c in G.columns if c != "geometry"}; base = {k: ("" if v is None or (isinstance(v, float) and np.isnan(v)) else v) for k, v in base.items()}
    if i in split_rows:
        for s in S[S.row.astype(int) == i].itertuples():
            n, src = s.sub, s.source
            if s.piece in v2: n, src = v2[s.piece]
            elif str(s.visual) == "True": n, src = "", "readers differ or not read: " + s.second_read
            g = shapely.from_wkb(bytes.fromhex(s.wkb)); r = dict(base); r.update({"kide": p.sheet + "/" + n if n else p.sheet, "survey_no": p.sheet, "subdiv_no": n, "area_sqm": float(s.area), "area_acre": round(float(s.area) / 4046.8564224, 4),
                         "first_read": base.get("plot_no", ""), "plot_no": n, "subdiv_status": status(src, n), "topo_note": (str(base.get("topo_note", "")) + "; cut from the merged row " + p.kide).strip("; "), "_row": i, "_piece": s.piece, "_geom": g, "_new": True})
            rows.append(r)
        continue
    n, src = p.new_sub, p.source
    if p.row in v1: n, src = v1[p.row]
    elif p.visual == "True": n, src = p.layer_subdiv if p.layer_subdiv != p.sheet else "", "readers differ or not read: second read %s" % (p.second_read or p.beside or "nothing")
    keepk = p.status == "no sheet (puvi base)"
    if src.startswith("row merges"): n = p.layer_subdiv if p.layer_subdiv != p.sheet else ""          # not cut: the row keeps its number, flagged; its survey number is the sheet's
    r = dict(base); r.update({"kide": p.kide if keepk else (p.sheet + "/" + n if n else p.sheet), "survey_no": p.survey_no if keepk else p.sheet, "subdiv_no": p.layer_subdiv if keepk else n, "first_read": base.get("plot_no", ""), "plot_no": base.get("plot_no", "") if keepk else n,
              "subdiv_status": status(src, n), "_row": i, "_piece": "", "_geom": None, "_new": False}); rows.append(r)
F = pd.DataFrame(rows); F["district_code"] = F.village_code.str.split("_").str[0]; F["taluk_code"] = F.village_code.str.split("_").str[1]; F["village_no"] = F.village_code.str.split("_").str[2]
F["district"] = "Chengalpattu"; F["taluk"] = F.taluk_code.map(TALUK); assert F.taluk.notna().all()
F["unique_id"] = [f"{v}_{s}/{n if n else '-'}" for v, s, n in zip(F.village_code, F.survey_no, F.subdiv_no)]
F.to_pickle(os.path.join(HERE, "out", "final_rows.pkl"))
print("rows", len(F), "| new pieces", int(F._new.sum()), "| from", len(split_rows), "merged rows"); print(F.subdiv_status.str.replace(r" \(.*| - readers.*| - row merges.*| - number would.*", "", regex=True).value_counts().to_string())
bad = F[(F.kide != F.survey_no + "/" + F.subdiv_no) & (F.subdiv_no != "")]; bad2 = F[(F.subdiv_no == "") & (F.kide != F.survey_no)]; print("kide not survey/subdiv:", len(bad), "| kide not survey where no subdivision:", len(bad2), "| blank kide:", int((F.kide == "").sum()))
d = F[(F.subdiv_no != "") & F.duplicated(["village_code", "kide"], keep=False)]; print("numbers that repeat in a village:", len(d), "rows,", d.drop_duplicates(["village_code", "kide"]).shape[0], "numbers")
