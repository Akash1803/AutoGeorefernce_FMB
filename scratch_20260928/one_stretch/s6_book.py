"""One Stretch, step 6: the record beside the layer - what was closed, what is left for hand work, which numbers wait."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
R = r"D:\Projects\Tambaram_Chengalpattu"; P = os.path.join(R, "Manual_Georeferenced", "Tambaram_Chengalpattu_One_Stretch.geojson")
B = gpd.read_file(P); M = B.to_crs(32644); g = list(M.geometry); n = len(g); tr = shapely.STRtree(g); nm = lambda i: "%s %s" % (B.village[i], B.kide[i])
wmax = lambda x: 2 * shapely.maximum_inscribed_circle(x).length if x.area > 0 else 0.0; rows = []
for i in range(n):
    for j in tr.query(g[i]):
        if j > i:
            x = g[i].intersection(g[j])
            for q in [q for q in shapely.get_parts(x) if q.geom_type == "Polygon" and q.area > 1e-6] if x.area > 0 else []:
                c = gpd.GeoSeries([q.representative_point()], crs=32644).to_crs(4326).iloc[0]; rows.append({"Kind": "overlap", "Between": nm(i) + "  |  " + nm(j), "Layers": " | ".join(sorted({B.source_layer[i], B.source_layer[j]})), "Area (sq m)": round(q.area, 2), "Widest (m)": round(wmax(q), 2), "Length (m)": round(q.length / 2, 1), "Longitude": round(c.x, 7), "Latitude": round(c.y, 7)})
U = unary_union(g)
for q in shapely.get_parts(U):
    for r in q.interiors:
        h = shapely.Polygon(r); near = [k for k in tr.query(h.buffer(0.01)) if g[k].distance(h) < 1e-6]; c = gpd.GeoSeries([h.representative_point()], crs=32644).to_crs(4326).iloc[0]
        rows.append({"Kind": "gap", "Between": "  |  ".join(sorted(nm(k) for k in near)[:6]) + (" ..." if len(near) > 6 else ""), "Layers": " | ".join(sorted({B.source_layer[k] for k in near})), "Area (sq m)": round(h.area, 2), "Widest (m)": round(wmax(h), 2), "Length (m)": round(h.length / 2, 1), "Longitude": round(c.x, 7), "Latitude": round(c.y, 7)})
L = pd.DataFrame(rows).sort_values(["Kind", "Area (sq m)"], ascending=[True, False])
log = pd.read_csv("out/topo_log.csv").rename(columns={"what": "What", "plot": "Plot that moved", "onto": "Onto the line of", "sqm": "Area (sq m)", "widest_mm": "Widest (mm)"})
notes = pd.read_csv("out/merge_notes.csv").rename(columns={"layer": "Layer", "row": "Row in that layer", "survey": "Survey No", "note": "Note"})
H = B[B.subdiv_status.str.startswith(("to be checked", "not checked"))][["source_layer", "village", "kide", "survey_no", "subdiv_no", "area_acre", "area_sqm", "in_rail_buffer_30m", "subdiv_status"]].rename(columns={"source_layer": "Layer", "village": "Village", "kide": "Kide", "survey_no": "Survey No", "subdiv_no": "Sub Div No (unverified)", "area_acre": "Area (acre)", "area_sqm": "Area (sq m)", "in_rail_buffer_30m": "In 30 m rail buffer", "subdiv_status": "Why"})
ch = B[B.first_read.fillna("") != B.subdiv_no.fillna("")]; ch = ch[ch.source_layer != "Manual_Georeference"][["village", "survey_no", "first_read", "subdiv_no", "kide", "area_acre", "subdiv_status"]].rename(columns={"village": "Village", "survey_no": "Survey No", "first_read": "Sub Div No (in the team layer)", "subdiv_no": "Sub Div No (now)", "kide": "Kide (now)", "area_acre": "Area (acre)", "subdiv_status": "Verified by"})
mv = B.seam_moved_mm.astype(float)
S = pd.DataFrame([("Plots in the layer", len(B)), ("   from 3VILLAGES", int((B.source_layer == "3VILLAGES").sum())), ("   from Peramanur_Merged_Fmb", int((B.source_layer == "Peramanur_Merged_Fmb").sum())), ("   from Manual_Georeference (not one vertex changed)", int((B.source_layer == "Manual_Georeference").sum())),
                  ("Hairline overlaps closed (under 5 cm at the widest)", int((log.What == "hairline overlap closed").sum())), ("Hairline gaps closed (under 5 cm at the widest)", int((log.What == "hairline gap closed").sum())), ("Plots whose line moved by 1 mm or more", int((mv >= 1).sum())),
                  ("Overlaps left for hand work", int((L.Kind == "overlap").sum())), ("Gaps left for hand work", int((L.Kind == "gap").sum())), ("Invalid / multi-part / duplicate shapes", 0),
                  ("Plots with Tamilnilam data", int((B.uid_tn_avail == 1).sum())), ("Plots in the 30 m rail buffer", int((B.in_rail_buffer_30m == "yes").sum())), ("   of these with Tamilnilam data", int(((B.in_rail_buffer_30m == "yes") & (B.uid_tn_avail == 1)).sum())),
                  ("Plot numbers still to check by hand", len(H))], columns=["Item", "Count"])
out = os.path.join(R, "Manual_Georeferenced", "Tambaram_Chengalpattu_One_Stretch_Record.xlsx")
def fit(w, name, d):
    ws = w.sheets[name]; ws.freeze_panes(1, 0)
    for i, c in enumerate(d.columns): ws.set_column(i, i, min(max(len(str(c)), int(d[c].astype(str).str.len().quantile(0.95)) if len(d) else 0) + 2, 80))
with pd.ExcelWriter(out, engine="xlsxwriter") as w:
    for nm_, d in (("Summary", S), ("Left for hand work", L), ("Seams closed", log), ("Numbers corrected", ch), ("Numbers to check by hand", H), ("Rows of the source layers", notes)): d.to_excel(w, sheet_name=nm_, index=False); fit(w, nm_, d)
print(S.to_string(index=False)); print(L.groupby("Kind")["Area (sq m)"].agg(["count", "sum", "max"]).round(1).to_string()); print(L.groupby(["Kind", "Layers"]).size().to_string()); print("numbers corrected", len(ch))
