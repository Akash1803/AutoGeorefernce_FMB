"""Dual QC of FMB_Fine_Tunned, step 4: the change plan. A number is set or changed only when two independent sources
agree (second read of the PDF + A-register, or second read + a look at the PDF crop); everything else goes to the
visual list. -> out/plan.csv (one row per layer row), out/visual.csv (rows to look at on the PDF)"""
import os, re, json, warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
D = pd.read_csv(os.path.join(HERE, "out", "decide5.csv"), dtype=str).fillna(""); D["area"] = D.area.astype(float)
V = {}
if os.path.exists(os.path.join(HERE, "out", "visual_read.csv")):                        # what the look at the PDF said: row -> number ("" = nothing printed)
    V = pd.read_csv(os.path.join(HERE, "out", "visual_read.csv"), dtype=str).fillna("").set_index("row").final.to_dict()
MR = pd.read_csv(os.path.join(HERE, "out", "merged_rows.csv"), dtype=str).fillna(""); merged = set(MR[MR.distinct.astype(int) >= 2].row)
new, src, need = [], [], []
for r in D.itertuples():
    if str(r.row) in merged and r.status != "no sheet (puvi base)":          # one row, several different plots: no verdict for the row, its plots are dealt with one by one
        new.append(r.subdiv_no); src.append("row merges several plots of the sheet - to be split"); need.append(False); continue
    A = r.subdiv_no; B = [t for t in r.second_read.split("/") if t]; reg = set(t for t in r.cand_in_register.split("/") if t); fit = set(t for t in r.cand_area_fits.split("/") if t); st = r.status; n, s, v = A, "", False
    if str(r.row) in V: n, s = V[str(r.row)], "second read + checked on the PDF crop" if V[str(r.row)] else "checked on the PDF crop: no number printed"
    elif st.startswith("confirmed"): s = "both reads agree" + (" (" + st.split("(")[1] if "(" in st else "")
    elif st == "whole survey": n, s = "", "whole survey - the sheet has one plot"
    elif st == "no sheet (puvi base)": s = "no FMB sheet - Puvi parcel"
    elif st == "not linked to a sheet polygon": s = "not checked - the row could not be tied to a sheet polygon"
    elif st == "no number printed": n, s = "", "no number printed on the sheet (both reads)"
    elif st == "number found by the second read":
        if r.proposed in fit: n, s = r.proposed, "second read + A-register"
        else: v = True
    elif st == "number printed beside the plot":
        if r.proposed in fit: n, s = r.proposed, "second read (printed beside the plot) + A-register"
        else: v = True
    elif st == "reads differ":
        c = [t for t in B if t in fit]
        if len(c) == 1 and r.layer_no_in_register != "yes": n, s = c[0], "second read + A-register (first read was %s)" % A
        else: v = True
    elif st == "several numbers in the plot":
        c = [t for t in B if t in fit]
        if len(c) == 1: n, s = c[0], "second read + A-register (other print in the plot: %s)" % "/".join(t for t in B if t != c[0])
        else: v = True
    elif st == "second read found nothing":
        if r.layer_no_area_fits == "yes": s = "first read + A-register"
        else: v = True
    new.append(n); src.append(s); need.append(v)
D["new_sub"], D["source"], D["visual"] = new, src, need
# one number, one plot: a number set on two rows of a sheet is looked at, unless the rows are near pieces kept apart before
for (v, u), g in D[(D.new_sub != "") & ~D.visual].groupby(["village_code", "sheet"]):
    for n, h in g.groupby("new_sub"):
        if len(h) > 1 and not (h.subdiv_no == n).all():
            ch = h[h.subdiv_no != n]; D.loc[ch.index, "visual"] = True; D.loc[ch.index, "source"] = "number would repeat on the sheet (%s)" % n; D.loc[ch.index, "new_sub"] = ch.subdiv_no
D["merged"] = D.source.str.startswith("row merges"); D["new_survey"] = D.sheet; D["new_kide"] = [u + "/" + n if n else u for u, n in zip(D.sheet, D.new_sub)]
keep_kide = D.status.isin(["no sheet (puvi base)"]); D.loc[keep_kide, ["new_kide", "new_survey", "new_sub"]] = D.loc[keep_kide, ["kide", "survey_no", "subdiv_no"]].values
import geopandas as gpd
_L = gpd.read_file(os.path.join(HERE, "out", "layer_before.geojson"), ignore_geometry=True).fillna(""); assert len(_L) == len(D) and (_L.kide.values == D.kide.values).all()
D["layer_subdiv"] = _L.subdiv_no.astype(str).values
D["changes"] = [(a != b) or (c != d) or (e != f) for a, b, c, d, e, f in zip(D.kide, D.new_kide, D.survey_no, D.new_survey, D.layer_subdiv, D.new_sub)]
D.to_csv(os.path.join(HERE, "out", "plan.csv"), index=False); D[D.visual].to_csv(os.path.join(HERE, "out", "visual.csv"), index=False)
print("rows", len(D), "| to look at on the PDF:", int(D.visual.sum()), "| rows that change (apart from those):", int((D.changes & ~D.visual).sum()))
print(D[~D.visual].source.str.replace(r"\(first read was .*\)|\(other print.*\)|\(.*overlapped.*\)|\(number printed beside.*\)", "", regex=True).str.strip().value_counts().to_string())
print("to look at, by reason:", D[D.visual].status.where(~D[D.visual].source.str.startswith("number would repeat"), "number would repeat").value_counts().to_dict())
c = D[D.changes & ~D.visual]; print("changes: number set where there was none", ((c.layer_subdiv == "") & (c.new_sub != "")).sum(), "| number replaced", ((c.layer_subdiv != "") & (c.new_sub != "") & (c.layer_subdiv != c.new_sub)).sum(), "| the sheet number taken out of the subdivision field", ((c.layer_subdiv != "") & (c.new_sub == "")).sum(), "| survey_no corrected", (c.survey_no != c.new_survey).sum())
