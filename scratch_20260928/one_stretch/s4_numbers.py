"""One Stretch, step 4: the attribute table. Plot numbers of the 3Villages and Peramanur plots are set by the same
two-source rule as FMB_Fine_Tunned (both reads of the PDF agree, or second read + A-register with a fitting extent);
nothing is cut here - a plot that covers several plots of its sheet stays one plot and says so.
Manual_Georeference rows keep the attributes of their own layer. Then the strict Tamilnilam match.
-> out/one_stretch_table.pkl"""
import os, re, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"; QC = os.path.join(HERE, "..", "qcone", "out"); sys.path.insert(0, os.path.join(HERE, "..", "qcfmb"))
import q9_tamilnilam as q9
TALUK = {"04": "Chengalpattu", "05": "Tambaram", "15": "Vandalur"}
G = gpd.read_file(os.path.join(HERE, "out", "one_stretch_geom.gpkg")); Gm = G.to_crs(32644)
M = gpd.read_file(os.path.join(R, "Manual_Georeferenced", "Manual_Georeference.geojson"), ignore_geometry=True).fillna("")
D = pd.read_csv(os.path.join(QC, "decide5.csv"), dtype=str).fillna(""); K = pd.read_csv(os.path.join(QC, "link_all.csv"), dtype=str).fillna(""); K = K[K.state == "ok"]; K = K[K.share.astype(float) >= 0.5]
L = pd.read_csv(os.path.join(QC, "labels5_all.csv"), dtype=str).fillna(""); P5 = pd.read_csv(os.path.join(HERE, "..", "qcfmb", "out", "polys5_all.csv"), dtype=str).fillna("")
nsheet = P5.groupby(["village_code", "sheet"]).size().to_dict()
good = L[(L.plotno == "True") & (L.pid != "-1")].groupby(["village_code", "sheet", "pid"]).text.apply(lambda s: sorted(set(s))).to_dict()
BUF = unary_union(list(gpd.read_file(os.path.join(R, "Tambaram_Chengalpattu_Railway.gpkg"), layer="rail_buffer_30m").to_crs(32644).geometry))
Q = G[G.src != "Manual_Georeference"].reset_index(); assert len(Q) == len(D) and (Q.sheet.values == D.sheet.values).all()
rows = []
for k, (q, d) in enumerate(zip(Q.itertuples(), D.itertuples())):
    A = d.subdiv_no; B = [t for t in d.second_read.split("/") if t]; fit = set(t for t in d.cand_area_fits.split("/") if t); st = d.status; pids = [p for p in d.pids.split("/") if p]
    per = [tuple(t for t in good.get((d.village_code, d.sheet, p), []) if t != d.sheet) for p in pids]; distinct = len({x for x in per if x}); ns = nsheet.get((d.village_code, d.sheet), 0)
    n, s = A, ""
    if q.settled: s0 = q.settled
    else: s0 = ""
    if st == "no sheet (puvi base)": n, s = "", "no FMB sheet on the portal for this survey number"
    elif len(pids) >= 2 and distinct >= 2: n, s = "", "to be checked by hand - the plot covers %d plots of its sheet (%s ...), their lines are not drawn in this layer" % (len(pids), ", ".join(sorted({x[0] for x in per if x})[:4]))
    elif st == "whole survey" and ns > 1: n, s = "", "to be checked by hand - the plot is the outline of sheet %s, which has %d plots; their lines are not drawn in this layer" % (d.sheet, ns)
    elif st.startswith("confirmed"): s = "printed on the FMB sheet - both reads agree"
    elif st == "whole survey": n, s = "", "whole survey - the FMB sheet has no subdivision"
    elif st == "not linked to a sheet polygon": s = "not checked - the plot could not be tied to a sheet polygon"
    elif st == "no number printed": n, s = "", "no number printed on the FMB sheet"
    elif st in ("number found by the second read", "number printed beside the plot") and d.proposed in fit: n, s = d.proposed, "printed on the FMB sheet - second read, A-register agrees"
    elif st == "reads differ" and len([t for t in B if t in fit]) == 1 and d.layer_no_in_register != "yes": n, s = [t for t in B if t in fit][0], "printed on the FMB sheet - second read, A-register agrees (first read was %s)" % A
    elif st == "several numbers in the plot" and len([t for t in B if t in fit]) == 1: n, s = [t for t in B if t in fit][0], "printed on the FMB sheet - second read, A-register agrees"
    elif st == "second read found nothing" and d.layer_no_area_fits == "yes": s = "printed on the FMB sheet - first read, A-register agrees"
    else: s = "to be checked by hand - %s (layer %s, second read %s)" % (st, A or "none", d.second_read or d.beside or "nothing")
    if s0: s = ("to be checked by hand - " + s0) if not s.startswith("to be checked") else s + "; " + s0
    rows.append({"i": q.index, "survey_no": d.sheet, "subdiv_no": n, "first_read": A, "subdiv_status": s, "status": "Georeferenced by hand (team layer %s)" % q.src, "source_pdf": "FMB_Sketches/%s/%s.pdf" % (d.village_code, d.sheet) if os.path.exists(os.path.join(R, "FMB_Sketches", d.village_code, d.sheet + ".pdf")) else ""})
T = pd.DataFrame(rows).set_index("i")
# one number, one plot on a sheet: a number newly set on two plots is looked at instead
T["vc"] = G.village_code[T.index].values
for (v, s), g in T[T.subdiv_no != ""].groupby(["vc", "survey_no"]):
    for n_, h in g.groupby("subdiv_no"):
        if len(h) > 1 and (h.first_read != n_).any():
            for i in h.index[h.first_read != n_]: T.at[i, "subdiv_status"] = "to be checked by hand - number %s would repeat on the sheet" % n_; T.at[i, "subdiv_no"] = T.at[i, "first_read"]
out = []
for i, g in G.iterrows():
    a = float(Gm.geometry[i].area)
    if g.src == "Manual_Georeference":
        m = M.iloc[int(g.src_row)]; r = {"kide": m.kide, "survey_no": m.survey_no, "subdiv_no": m.subdiv_no, "first_read": "", "subdiv_status": m.subdiv_status, "status": m.status, "source_pdf": m.source_pdf, "village": m.village, "village_code": m.village_code}
    else:
        t = T.loc[i]; r = {"kide": t.survey_no + "/" + t.subdiv_no if t.subdiv_no else t.survey_no, "survey_no": t.survey_no, "subdiv_no": t.subdiv_no, "first_read": t.first_read, "subdiv_status": t.subdiv_status, "status": t.status, "source_pdf": t.source_pdf, "village": g.village, "village_code": g.village_code}
    r.update({"area_acre": round(a / 4046.8564224, 4), "area_sqm": round(a, 1), "district": "Chengalpattu", "taluk": TALUK[r["village_code"].split("_")[1]], "district_code": r["village_code"].split("_")[0], "taluk_code": r["village_code"].split("_")[1], "village_no": r["village_code"].split("_")[2],
              "in_rail_buffer_30m": "yes" if Gm.geometry[i].intersects(BUF) else "no", "source_layer": g.src, "unique_id": "%s_%s/%s" % (r["village_code"], r["survey_no"], r["subdiv_no"] or "-"), "seam_moved_mm": float(g.moved_mm), "_i": i})
    out.append(r)
F = pd.DataFrame(out); F = q9.match(F)
first = ["kide", "survey_no", "subdiv_no", "area_acre", "area_sqm", "village", "village_code", "taluk", "taluk_code", "district", "district_code", "village_no", "unique_id", "status", "subdiv_status", "first_read", "in_rail_buffer_30m", "source_layer", "source_pdf", "seam_moved_mm"]
F = F[first + [c for c in F.columns if c not in first]]; F.to_pickle(os.path.join(HERE, "out", "one_stretch_table.pkl"))
print("plots", len(F)); print(F.groupby(["source_layer", "village"]).size().to_string())
print(F.subdiv_status.str.replace(r" \(.*| - the plot.*| - number .*| - reads.*| - second.*| - several.*| - survey.*| - village.*| - not linked.*", "", regex=True).value_counts().to_string())
q9.summary(F, "\nTamilnilam"); b = F[F.in_rail_buffer_30m == "yes"]; print("\nin the 30 m rail buffer:", len(b), "| with Tamilnilam data:", int((b.uid_tn_avail == 1).sum()))
assert (F.kide != "").all() and (F.kide == F.survey_no + np.where(F.subdiv_no != "", "/" + F.subdiv_no, "")).all()
m = F[F.source_layer == "Manual_Georeference"].reset_index(drop=True); print("Manual rows: Tamilnilam verdict equal to the Manual layer's own:", int((m.tamilnilam_match.str.split(" \\(").str[0].values == M.tamilnilam_match.str.split(" \\(").str[0].values).sum()), "of", len(m))
