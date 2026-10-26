import sys, json, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd, geopandas as gpd
from shapely import make_valid
A = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); A["kide"] = A.kide.astype(str); A["vc"] = A.village_code.astype(str); A["sv"] = A.kide.str.split("/").str[0]
T = pd.read_excel("D:/Projects/CMP_Roads_LPS/Merged_Data/02_Tamilnilam_output.xlsx", dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3)
reg367 = set(T[(T.vc == "35_11_047") & (T.fmb_survey_no == "367")].sub_division); print("367 register has 4C3:", "4C3" in reg367, "| 4C1A1:", "4C1A1" in reg367)
f = json.load(open("fix_two_labels.json")); ch = f["change"]
rp = lambda i: list(make_valid(A.geometry[i]).representative_point().coords[0])
def entry(i, new, note):
    kide = f"{A.sv[i]}/{new}" if new else A.sv[i]
    return {"vc": A.vc[i], "old_kide": A.kide[i], "old_pt": rp(i), "fid_hint": int(i), "attrs": {"kide": kide, "subdiv_no": new, "unique_id": f"{A.vc[i]}_{kide}" + ("" if new else "/-"), "pdf_check": note}}
# 1. drop the two wrong sliver fixes (668 sliver -> 56, 260 strip -> 44B)
before = len(ch); ch = [c for c in ch if not (c["old_kide"] == "668" and c["attrs"]["kide"] == "668/56") and not (c["old_kide"] == "260" and c["attrs"]["kide"] == "260/44B")]
print("dropped sliver fixes:", before - len(ch))
# 2. 367/4C1A1 existing plot: both readers read 4C3
i = A.index[(A.vc == "35_11_047") & (A.kide == "367/4C1A1")][0]
ch.append(entry(i, "4C3" if "4C3" in reg367 else "", "both independent readers read 4C3 printed for this plot (4C1A1 is the plot to its west) - settled 2026-10-26" if "4C3" in reg367 else "label 4C1A1 belongs to the plot to its west; printed number unclear - cleared 2026-10-26"))
# 3. 904/2A2A piece that lies outside the sheet (the one not being relabelled to 1A) -> cleared
fixed = {c["fid_hint"] for c in ch}
for i in A.index[(A.vc == "35_11_047") & (A.kide == "904/2A2A")]:
    if i not in fixed: ch.append(entry(i, "", "second 904/2A2A piece lies outside the 904 sheet (east of the 979 line); the real 2A2A is the former 2A3 plot - cleared 2026-10-26"))
# 4. 43/2 thin sliver (team artifact) -> cleared; the real plot 2 is the unnumbered 43 piece being labelled
for i in A.index[(A.vc == "35_11_117") & (A.kide == "43/2")]:
    if i not in fixed: ch.append(entry(i, "", "thin sliver carrying label 2; the real plot 2 is the large piece beside it - cleared 2026-10-26"))
f["change"] = ch; json.dump(f, open("fix_two_labels.json", "w")); print("label changes now:", len(ch))
# duplicate check
L = A[["vc", "kide", "sv"]].copy()
for c in ch: L.loc[c["fid_hint"], "kide"] = c["attrs"]["kide"]
p = A[(A.vc == "35_11_047") & (A.kide == "344")]; L.loc[p.geometry.area.idxmax(), "kide"] = "344/2"; L = L.drop(index=[p.geometry.area.idxmin()])
for k in ("905/1A", "943/4B1A4", "36/3B2"): L = pd.concat([L, pd.DataFrame([{"vc": "35_11_047" if k != "36/3B2" else "35_11_117", "kide": k, "sv": k.split("/")[0]}])], ignore_index=True)
L["key"] = L.vc + "|" + L.kide; dup = L[L.kide.str.contains("/")].key.value_counts(); dup = dup[dup > 1]
print("duplicate numbered kide after all changes:", list(dup.index))
sys.exit(1 if len(dup) else 0)
