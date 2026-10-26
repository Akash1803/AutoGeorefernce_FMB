"""Survey-wise area check: printed sheet area (title block) vs layer area, one row per village+survey of the final layer.
Usage: area_check_excel.py <final.geojson> <out.xlsx>"""
import sys, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd, geopandas as gpd
F, OUT = sys.argv[1:3]
P = gpd.read_file(F); P["kide"] = P.kide.astype(str); P["vc"] = P.village_code.astype(str); P["sv"] = P.kide.str.split("/").str[0]
H = pd.read_csv("pdf_header_areas.csv", dtype={"village_code": str, "survey": str}); H["pdf_area_m2"] = pd.to_numeric(H.pdf_area_m2, errors="coerce")
hd = {(r.village_code, r.survey): r.pdf_area_m2 for r in H.itertuples()}; have_pdf = set(hd)
g = P.groupby(["vc", "sv"]).agg(taluk=("taluk", "first"), village=("village", "first"), plots=("kide", "size"), layer_m2=("geometry", lambda s: s.area.sum()),
                               numbered=("kide", lambda s: s.str.contains("/").sum())).reset_index()
rows = []
def unit_key(vc, sv):
    if (vc, sv) in have_pdf: return (vc, sv)
    m = re.match(r"(\d+)([A-Z]+)$", sv)                      # 109A -> sheet 109 when only the parent sheet exists
    if m and (vc, m.group(1)) in have_pdf: return (vc, m.group(1))
    return None
for r in g.itertuples():
    k = unit_key(r.vc, r.sv); pdf = hd.get(k) if k else None
    if k and k[1] != r.sv: pdf_note = f"sheet {k[1]} (letter units share one sheet)"
    elif k and pd.isna(pdf): pdf_note = "area not printed on the sheet"; pdf = None
    elif k is None: pdf_note = "no sketch downloaded"
    else: pdf_note = ""
    diff = None if pdf is None or not pdf else round((r.layer_m2 / pdf - 1) * 100, 1)
    if diff is None: rem = pdf_note
    elif abs(diff) <= 5: rem = "OK (within 5 %)"
    elif abs(diff) <= 15: rem = "check (5-15 %)"
    else: rem = "MISMATCH (> 15 %)" + (" - only part of the survey lies in the corridor layer" if r.layer_m2 < pdf else "")
    rows.append({"Taluk": r.taluk, "Village": r.village, "Village code": r.vc, "Survey No": r.sv, "PDF area (sq m)": None if pdf is None else round(pdf, 1),
                 "Layer area (sq m)": round(r.layer_m2, 1), "PDF area (acre)": None if pdf is None else round(pdf / 4046.856, 3), "Layer area (acre)": round(r.layer_m2 / 4046.856, 3),
                 "Difference (%)": diff, "Plots in layer": r.plots, "Numbered plots": r.numbered, "Remark": rem})
X = pd.DataFrame(rows).sort_values(["Taluk", "Village", "Survey No"], key=lambda s: s.map(lambda v: (int(re.match(r"\d+", str(v)).group()), str(v)) if re.match(r"\d+", str(v)) else (10**9, str(v))) if s.name == "Survey No" else s)
X.insert(0, "Sl.No", range(1, len(X) + 1))
with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    X.to_excel(xw, index=False, sheet_name="Survey area check")
    ws = xw.sheets["Survey area check"]
    for col, w in zip("ABCDEFGHIJKLM", [6, 16, 18, 12, 10, 16, 16, 14, 14, 13, 13, 14, 60]): ws.column_dimensions[col].width = w
    ws.freeze_panes = "A2"
print("rows", len(X), "| with PDF area", X["PDF area (sq m)"].notna().sum(), "| remarks:", X.Remark.str.split(" ").str[0].value_counts().to_dict())
