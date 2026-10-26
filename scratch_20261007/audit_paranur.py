"""Read-only audit of Paranur (35_04_221): every sheet PDF vs vector vs register vs placed geometry; topology of
everything placed (checked twice: shapely on the files, QGIS on the layers); PDF-neighbour seams; causes.
Writes Auto_Georef/35_04_221/Paranur_Audit_<date>.xlsx and Paranur_Audit_<date>.gpkg (overlaps, seam gaps)."""
import sys, os, json, math, datetime
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
from autogeoref import ring, ring_io, printed_numbers as pn, paths, match, fit, qgis_bridge
V, VN = "35_04_221", "Paranur"; PR = paths.PROJECT
K2 = 110540.0 * 108600.0
EYE = {"264": 7900, "71B": 5000, "71A": 300, "72": 3900, "74": 5330, "66A": 12150, "66B": 2050, "67A": 6750, "67B": 3300,
       "68A": 1800, "68B": 1800, "69": 750, "70A": 2550, "70B": 10250, "265": 2600, "266": 3950, "267": 5200, "62A": 156900,
       "62B": 14950, "73": 16500, "76A": 11450, "76B": 3550, "77": 1800, "78A": 16650, "78B": 5150, "79A": 7450, "79B": 1850,
       "80": 16100, "103": 17100, "268": 3950, "269": 6050}
table = pn.village_table(V, str(PR / "Auto_Georef" / V / "printed_numbers.json"))
agreed = pn.agreed_sets(table); printed = pn.printed_sets(table)
# placed parcels as QGIS shows them (his + red)
lay = ring_io.village_layers(VN)
placed, who = {}, {}
for L in lay:
    if "_parcels" not in L["name"]:
        continue
    s = L["name"].split("_parcels")[0]
    g = ring_io._read_source(L["source"]); g = g[g.geometry.notna()]
    placed[s] = g.to_crs(4326); who[s] = "red (auto)" if L["name"].endswith("_auto") else "yours"
F = gpd.read_file(PR / "FMB_Georeferenced_Final" / "FMB_Fine_Tunned.geojson")
reg = F[F.village_code == V].groupby(F.survey_no.astype(str)).tamilnilam_extent_sqm.sum(min_count=1)
regtype = F[F.village_code == V].groupby(F.survey_no.astype(str)).tamilnilam_remarks.agg(lambda x: "; ".join(sorted({str(v) for v in x if str(v) not in ("", "nan", "None")}))[:60])

rows = []
for s in sorted(table, key=paths.survey_sort_key):
    sh = ring.load_sheet(V, s, table)
    pdf_auto = pn.read_header_area(paths.sketch_dir(V) / ("%s.pdf" % s))
    vec = sum(p.area for _pr, p in sh.plots)
    r = {"survey": s, "placed": who.get(s, "not placed"), "pdf_title_m2 (read 1)": pdf_auto, "pdf_title_m2 (read 2, eye)": EYE.get(s),
         "register_m2": reg.get(s), "register_note": regtype.get(s, ""), "vector_m2": round(vec, 1),
         "plots": len(sh.plots), "plots_without_number": sum(1 for pr, _p in sh.plots if not pr.get("numbered") and len(sh.plots) > 1),
         "plots_to_confirm": sum(1 for pr, _p in sh.plots if pr.get("review_needed")),
         "printed_neighbours": ", ".join(sorted(printed.get(s, set()))), "two_way_neighbours": ", ".join(sorted(agreed.get(s, set())))}
    if s in placed:
        G = placed[s].to_crs(32644)
        U = unary_union(G.geometry.values).buffer(0)
        r["placed_m2"] = round(U.area, 1)
        r["placed_vs_vector_%"] = round(100 * (U.area / vec - 1), 1)
        r["placed_vs_register_%"] = round(100 * (U.area / reg[s] - 1), 1) if s in reg and reg[s] else None
        poly = ring.largest(U)
        ch = match.common_chains(sh.outline, ring.clockwise(poly.exterior.coords), tol_abs=4.0, tol_rel=0.12, min_len=5.0)
        if ch:
            th, t, rms, mx = fit.rigid_fit([a for a, _b in ch[0].pairs], [b for _a, b in ch[0].pairs])
            Rg = ring.ground_polygon(sh, th, t)
            # refine the rigid pose onto his outline (all vertices to nearest boundary) for a fair shape score
            for _ in range(3):
                P = np.array(Rg.exterior.coords)[:-1]
                Q = np.array([shapely.ops.nearest_points(poly.exterior, Point(p))[0].coords[0] for p in P])
                th2, t2, rms2, mx2 = fit.rigid_fit(sh.outline, Q)
                Rg = ring.ground_polygon(sh, th2, t2)
            r["shape_iou_vs_pdf"] = round(Rg.intersection(poly).area / Rg.union(poly).area, 3)
            r["max_outline_dev_m"] = round(shapely.hausdorff_distance(Rg.exterior, poly.exterior), 1)
            sim = fit.similarity_fit(sh.outline, Q)
            r["scale_vs_pdf_%"] = round(100 * (sim[3] - 1), 1)
    rows.append(r)
S = pd.DataFrame(rows)

# ---- topology of everything placed, twice
allp = gpd.GeoDataFrame(pd.concat([g.assign(survey=s, who=who[s]) for s, g in placed.items()], ignore_index=True), geometry="geometry", crs=4326)
g = allp.geometry.values
t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects"); k = a_ < b_
ov, ovgeo = [], []
for i, j in zip(a_[k], b_[k]):
    if allp.survey.iat[i] == allp.survey.iat[j]:
        continue
    x = shapely.intersection(g[i], g[j])
    parts = [p for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon"]
    ar = sum(p.area for p in parts) * K2
    if ar > 0:
        ov.append((allp.survey.iat[i], allp.survey.iat[j], ar)); ovgeo.append(unary_union(parts))
O = pd.DataFrame(ov, columns=["a", "b", "m2"])
O = O.groupby(["a", "b"], as_index=False).m2.sum() if len(O) else O
code = """
import json
from qgis.core import QgsProject, QgsGeometry
P = QgsProject.instance(); sub = P.layerTreeRoot().findGroup("Manual_Georeferenced").findGroup("%s_Village")
F = []
for n in sub.findLayers():
    L = n.layer()
    if L is None or "_parcels" not in L.name(): continue
    s = L.name().split("_parcels")[0]
    F += [(s, QgsGeometry(f.geometry())) for f in L.getFeatures() if not f.geometry().isEmpty()]
out = {}
for i in range(len(F)):
    for j in range(i + 1, len(F)):
        if F[i][0] == F[j][0] or not F[i][1].boundingBox().intersects(F[j][1].boundingBox()): continue
        x = F[i][1].intersection(F[j][1])
        if x and not x.isEmpty() and x.area() > 0:
            key = "|".join(sorted([F[i][0], F[j][0]])); out[key] = out.get(key, 0) + x.area() * %r
U = QgsGeometry.unaryUnion([g for s, g in F]); holes = []
for part in U.asGeometryCollection():
    for r in part.asPolygon()[1:]:
        h = QgsGeometry.fromPolygonXY([r]).difference(U)
        if h.area() > 0: holes.append(h.area() * %r)
print(json.dumps({"pairs": out, "holes": holes}))
""" % (VN, K2, K2)
Q2 = qgis_bridge.json_result(code, timeout=300)
O["m2_qgis"] = [Q2["pairs"].get("|".join(sorted([x, y])), 0.0) for x, y in zip(O.a, O.b)] if len(O) else []

def cause(a, b):
    out = []
    for s in (a, b):
        r = S[S.survey == s].iloc[0]
        if pd.notna(r.get("placed_vs_vector_%")) and abs(r["placed_vs_vector_%"]) > 5:
            out.append("%s is %+.0f%% vs its PDF drawing" % (s, r["placed_vs_vector_%"]))
        elif pd.notna(r.get("max_outline_dev_m")) and r["max_outline_dev_m"] > 3:
            out.append("%s shape differs from its PDF by up to %.0f m" % (s, r["max_outline_dev_m"]))
    nb = b in agreed.get(a, set())
    out.append("PDF neighbours" if nb else "PDFs do NOT name each other")
    return "; ".join(out)
if len(O):
    O["pdf_neighbours"] = [b in agreed.get(a, set()) for a, b in zip(O.a, O.b)]
    O["cause"] = [cause(a, b) for a, b in zip(O.a, O.b)]
    O = O.sort_values("m2", ascending=False)

# ---- seams between PDF neighbours that are both placed
U_s = {s: unary_union(g.to_crs(32644).geometry.values).buffer(0) for s, g in placed.items()}
seam, gapgeo = [], []
for a in U_s:
    for b in agreed.get(a, set()):
        if b not in U_s or a > b:
            continue
        A, B = U_s[a], U_s[b]
        d = A.distance(B)
        pts, _ = ring.sheets.sample_outline(list(ring.largest(A).exterior.coords)[:-1], 1.0)
        dd = shapely.distance(shapely.points(pts), B.boundary)
        near = dd[dd < 15]
        both = A.union(B)
        wedge = both.buffer(8, join_style=2).buffer(-8, join_style=2).difference(both)
        wedge = unary_union([p for p in getattr(wedge, "geoms", [wedge]) if p.area > 1 and p.distance(A) < 0.05 and p.distance(B) < 0.05])
        seam.append({"a": a, "b": b, "touch": d < 0.01, "min_gap_m": round(d, 2),
                     "shared_run_m": int((near < 0.05).sum()), "run_within_15m": len(near),
                     "max_separation_on_run_m": round(float(near.max()), 2) if len(near) else None,
                     "open_wedge_m2": round(wedge.area, 1) if not wedge.is_empty else 0.0})
        if not wedge.is_empty and wedge.area > 1:
            gapgeo.append((a, b, wedge))
SE = pd.DataFrame(seam).sort_values(["touch", "open_wedge_m2"], ascending=[True, False])
# contacts the PDFs do not name
nonnb = []
for a in U_s:
    for b in U_s:
        if a < b and b not in agreed.get(a, set()) and U_s[a].distance(U_s[b]) < 0.05:
            L = U_s[a].boundary.intersection(U_s[b].buffer(0.05)).length
            if L > 5:
                nonnb.append({"a": a, "b": b, "shared_m": round(L, 1), "a_prints_b": b.upper() in printed.get(a, set()), "b_prints_a": a.upper() in printed.get(b, set())})
NN = pd.DataFrame(nonnb)

stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
out_x = PR / "Auto_Georef" / V / ("Paranur_Audit_%s.xlsx" % stamp)
with pd.ExcelWriter(out_x) as xw:
    S.to_excel(xw, sheet_name="Sheets", index=False)
    O.to_excel(xw, sheet_name="Overlaps", index=False)
    SE.to_excel(xw, sheet_name="PDF neighbour seams", index=False)
    NN.to_excel(xw, sheet_name="Contacts not in PDFs", index=False)
out_g = PR / "Auto_Georef" / V / ("Paranur_Audit_%s.gpkg" % stamp)
if ovgeo:
    og = gpd.GeoDataFrame([{"a": x, "b": y, "m2": round(m, 2)} for x, y, m in ov], geometry=ovgeo, crs=4326)
    og.to_file(out_g, layer="overlaps", driver="GPKG")
if gapgeo:
    gg = gpd.GeoDataFrame([{"a": a, "b": b, "m2": round(w.area, 1)} for a, b, w in gapgeo], geometry=[w for _a, _b, w in gapgeo], crs=32644).to_crs(4326)
    gg.to_file(out_g, layer="seam_gaps", driver="GPKG")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 60)
cols = ["survey", "placed", "pdf_title_m2 (read 1)", "pdf_title_m2 (read 2, eye)", "register_m2", "vector_m2", "placed_m2",
        "placed_vs_vector_%", "placed_vs_register_%", "shape_iou_vs_pdf", "max_outline_dev_m", "scale_vs_pdf_%", "plots_to_confirm"]
print(S[[c for c in cols if c in S]].to_string(index=False))
print("\nOVERLAPS (check 1 shapely / check 2 QGIS):"); print(O.round(2).to_string(index=False) if len(O) else "none")
print("holes (QGIS):", [round(h, 2) for h in Q2["holes"]])
print("\nSEAMS between PDF neighbours:"); print(SE.to_string(index=False))
print("\nCONTACTS NOT IN PDFs:"); print(NN.to_string(index=False) if len(NN) else "none")
print(out_x); print(out_g)
