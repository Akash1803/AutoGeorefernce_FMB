"""Turn the settled annotation decisions into QGIS fix files for CMP_two_Parts:
 fix_two_labels.json   - attribute changes (kide / subdiv_no / unique_id / pdf_check), FIX and CLEAR
 fix_two_merge344.json - geometry: the two unnumbered 344 pieces -> one 344/2
 fix_two_fill_ok.json  - the approved fills that still stand (905/1A, 943/4B1A4, 36/3B2)
Each decision row is tied to ONE layer feature: via the warped PDF plot centroid (chk_warped.gpkg) when the row has a pdf_plot,
else by (village, kide, layer area). Prints the duplicate-label check for the surveys touched."""
import sys, json, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd, geopandas as gpd, numpy as np
from shapely import make_valid
from shapely.ops import unary_union


def norm(s):
    s = str(s).strip()
    return "" if s.lower() in ("", "-", "nan", "none", "?", "null") else re.sub(r"[^0-9A-Z]", "", s.upper())


A = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson")
A["kide"] = A.kide.astype(str); A["vc"] = A.village_code.astype(str); A["key"] = A.vc + "|" + A.kide; A["sv"] = A.kide.str.split("/").str[0]
W = gpd.read_file("chk_warped.gpkg"); W["pdf_plot"] = W.pdf_plot.astype(int)
R1 = pd.read_csv("qc/two_read1.csv", dtype=str).fillna("")
M = R1[(R1.kide != "") & ~R1.verdict.isin(["LABEL NOT IN REGISTER", "PRINTED NUMBER NOT IN LAYER"]) & ~R1.verdict.str.startswith("DUPLICATE")].reset_index(drop=True)
D = pd.read_csv("qc/two_decisions.csv", dtype=str).fillna("")
assert len(M) == len(D), (len(M), len(D))
D["layer_area"] = M.layer_area.values; D["pdf_plot"] = M.pdf_plot.values
green = set(A.key[A.status.astype(str).str.startswith("green")])
done1 = pd.read_csv("qc/settle/decided_keys.csv", dtype=str).fillna("")
D2 = D[~D.status.isin(["OK"]) & ~D.key.isin(green) & ~D.key.isin(set(done1.key))].reset_index(drop=True)
I2 = pd.read_csv("qc/settle2/index.csv", dtype=str).fillna("")
assert len(I2) == len(D2), (len(I2), len(D2))
# index.csv lists the no-pose rows first, then the panels in D2 order -> align panels with D2 minus the no-pose keys
nopose = set(I2[I2.image == ""].key); pan = I2[I2.image != ""].reset_index(drop=True)
D2 = D2[~D2.key.isin(nopose)].reset_index(drop=True)
assert len(pan) == len(D2) and (pan.key.values == D2.key.values).all(), "panel/decision order mismatch"
dec2 = pd.read_csv("qc/settle2/decisions.csv", dtype=str).fillna("").rename(columns={"new_label": "settled"}); dec2["image"] = "qc/settle2/" + dec2.image
pan = pan.merge(dec2, on=["image", "panel"], how="left")
assert len(pan) == len(D2)
D2["verdict"] = pan.verdict.values; D2["settled"] = pan.settled.values; D2["note"] = pan.note.values


def feature_for(row):
    cands = A[A.key == row.key]
    if len(cands) == 1:
        return cands.index[0]
    if row.pdf_plot not in ("", "nan"):
        w = W[(W.village == row.village) & (W.survey == row.survey) & (W.pdf_plot == int(float(row.pdf_plot)))]
        if len(w):
            c = w.geometry.iloc[0].representative_point(); inside = cands[cands.contains(c)]
            if len(inside) == 1:
                return inside.index[0]
            if len(inside) == 0:
                j = cands.distance(c).idxmin()
                if cands.loc[j].geometry.distance(c) < 3:
                    return j
    if row.layer_area:
        ar = float(row.layer_area); d = (cands.geometry.area - ar).abs()
        if d.min() < 0.6:
            return d.idxmin()
    return None


def rp(i):
    return list(make_valid(A.geometry[i]).representative_point().coords[0])


labels, errs = [], []
for r in D2.itertuples():
    if r.verdict not in ("FIX", "CLEAR"):
        continue
    i = feature_for(r)
    if i is None:
        errs.append((r.key, r.verdict, r.new_label, "feature not identified")); continue
    sv = A.sv[i]; new = norm(r.settled) if r.verdict == "FIX" else ""
    if r.verdict == "FIX" and not new:
        errs.append((r.key, "FIX without label")); continue
    kide = f"{sv}/{new}" if new else sv
    note = (f"plot number {new} printed on the PDF - settled by machine read + two independent readers + PDF check 2026-10-26" if new
            else "no plot number for this piece on the PDF (its earlier label belonged to the plot beside it) - settled 2026-10-26")
    labels.append({"vc": A.vc[i], "old_kide": A.kide[i], "old_pt": rp(i), "fid_hint": int(i),
                   "attrs": {"kide": kide, "subdiv_no": new, "unique_id": f"{A.vc[i]}_{kide}" + ("" if new else "/-"), "pdf_check": note}})
for r in done1[done1.verdict == "FIX"].itertuples():
    c = A[A.key == r.key]
    if len(c) != 1:
        errs.append((r.key, "round-1 FIX: feature not unique")); continue
    i = c.index[0]; new = norm(r.new_label); kide = f"{A.sv[i]}/{new}"
    labels.append({"vc": A.vc[i], "old_kide": A.kide[i], "old_pt": rp(i), "fid_hint": int(i),
                   "attrs": {"kide": kide, "subdiv_no": new, "unique_id": f"{A.vc[i]}_{kide}",
                             "pdf_check": f"plot number {new} printed on the PDF - settled by two independent readers + PDF check 2026-10-26"}})
print("label changes:", len(labels), "| problems:", errs)
json.dump({"change": labels}, open("fix_two_labels.json", "w"))

p = A[(A.vc == "35_11_047") & (A.kide == "344")]
assert len(p) == 2, len(p)
u = make_valid(unary_union([g.buffer(0.001, join_style=2) for g in p.geometry])).buffer(-0.001, join_style=2)
u = max([z for z in getattr(u, "geoms", [u]) if z.geom_type == "Polygon"], key=lambda z: z.area)
keep, drop = p.geometry.area.idxmax(), p.geometry.area.idxmin()
merge = {"change": [{"vc": "35_11_047", "old_kide": "344", "old_pt": rp(keep), "wkt": u.wkt,
                     "attrs": {"kide": "344/2", "subdiv_no": "2", "unique_id": "35_11_047_344/2", "area_sqm": round(u.area, 1), "area_acre": round(u.area / 4046.856, 4),
                               "pdf_check": "PDF plot 2 (register 1150 m2) was split in two unnumbered pieces by the vector conversion - merged 2026-10-26",
                               "merge_note": "2026-10-26: two unnumbered 344 pieces merged into plot 2"}}],
         "add": [], "delete": [{"vc": "35_11_047", "old_kide": "344", "old_pt": rp(drop)}]}
print("344 merge:", [round(g.area) for g in p.geometry], "->", round(u.area), "valid", u.is_valid)
json.dump(merge, open("fix_two_merge344.json", "w"))

F = json.load(open("fix_two_fill.json")); keepn = {"905/1A", "943/4B1A4", "36/3B2"}; hosts = {"905/1B", "943/4B1B", "36/3A"}
adds = [a for a in F["add"] if a["attrs"]["kide"] in keepn]; chg = [c for c in F["change"] if c["old_kide"] in hosts]
json.dump({"change": chg, "add": adds, "delete": []}, open("fix_two_fill_ok.json", "w"))
print("fills kept:", [a["attrs"]["kide"] for a in adds], "hosts cut:", [c["old_kide"] for c in chg])

L = A[["vc", "kide", "sv"]].copy()
for c in labels:
    L.loc[c["fid_hint"], "kide"] = c["attrs"]["kide"]
L.loc[keep, "kide"] = "344/2"; L = L.drop(index=[drop])
for a in adds:
    L = pd.concat([L, pd.DataFrame([{"vc": a["attrs"]["village_code"], "kide": a["attrs"]["kide"], "sv": a["attrs"]["survey_no"]}])], ignore_index=True)
L["key"] = L.vc + "|" + L.kide
touched = {A.vc[c["fid_hint"]] + "|" + A.sv[c["fid_hint"]] for c in labels}
dup = L[L.kide.str.contains("/")].key.value_counts(); dup = dup[dup > 1]
dup_t = [k for k in dup.index if k.split("|")[0] + "|" + k.split("|")[1].split("/")[0] in touched]
print("duplicate numbered kide after the change, in touched surveys:", dup_t)
print("duplicates elsewhere (pre-existing):", [k for k in dup.index if k not in dup_t])
