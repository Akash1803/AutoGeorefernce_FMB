"""One Stretch, step 1: the three layers in one table, as received (no geometry touched), with the village of every plot.
-> out/merged_raw.gpkg (EPSG:32644, single polygons; src = layer, src_row = row in that layer)"""
import os, re, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
R = r"D:\Projects\Tambaram_Chengalpattu"
NAMES = {"35_04_052": "Thailavaram", "35_04_054": "Potheri", "35_04_056": "Kattankulathur", "35_04_076": "Peramanur", "35_04_074": "Thirukatchur", "35_04_077": "Kizhikaranai"}
clean = lambda x: "" if x is None or (isinstance(x, float) and np.isnan(x)) or str(x).strip().lower() in ("nan", "none", "null", "-", "0", "") else str(x).strip()
V = gpd.read_file("out/3VILLAGES_as_received.geojson"); P = gpd.read_file("out/Peramanur_as_received.geojson"); M = gpd.read_file(os.path.join(R, "Manual_Georeferenced", "Manual_Georeference.geojson"))
rows = []; notes = []
for name, G in (("3VILLAGES", V), ("Peramanur_Merged_Fmb", P)):
    for i, r in G.iterrows():
        if r.geometry is None or r.geometry.is_empty: notes.append((name, i, clean(r.get("survey_no")), "no geometry - left out")); continue
        s = clean(r.get("survey_no")); p = clean(r.get("plot_no")); lay = clean(r.get("layer"))
        if not s and lay: s = re.match(r"([0-9]+[A-Z]*)", lay).group(1); notes.append((name, i, s, "survey number taken from the source layer name " + lay))
        vc = clean(r.get("village_code")); parts = [q for q in shapely.get_parts(shapely.make_valid(r.geometry)) if q.geom_type == "Polygon" and q.area > 0]
        if not r.geometry.is_valid: notes.append((name, i, s, "invalid shape (%s) repaired without moving a vertex" % shapely.is_valid_reason(r.geometry)[:40]))
        for k, q in enumerate(sorted(parts, key=lambda q: -q.area)):
            rows.append({"src": name, "src_row": i, "part": k, "village_code": vc, "sheet": s, "first_read": "" if p == s else p, "geometry": q})
for i, r in M.iterrows():
    rows.append({"src": "Manual_Georeference", "src_row": i, "part": 0, "village_code": r.village_code, "sheet": r.survey_no, "first_read": clean(r.subdiv_no), "geometry": r.geometry})
G = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326).to_crs(32644)
# village of the plots that came without attributes: the village whose FMB sheet of that number exists and whose plots they touch
have = lambda v, s: os.path.exists(os.path.join(R, "FMB_Vector", v, s + "_parcels.geojson")) or os.path.exists(os.path.join(R, "FMB_Sketches", v, s + ".pdf"))
known = G[G.village_code != ""]; tr = shapely.STRtree(list(known.geometry)); kv = list(known.village_code)
for i in G.index[G.village_code == ""]:
    s = G.at[i, "sheet"]; cand = [v for v in ("35_04_052", "35_04_054", "35_04_056", "35_04_076") if have(v, s)]
    near = pd.Series([kv[j] for j in tr.query(G.geometry[i].buffer(30))]).value_counts()
    pick = [v for v in near.index if v in cand]
    if len(cand) == 1: G.at[i, "village_code"] = cand[0]
    elif pick: G.at[i, "village_code"] = pick[0]
    else: notes.append((G.src[i], G.src_row[i], s, "village not settled: sheets %s, neighbours %s" % (cand, near.to_dict())))
# four plots that came without a usable survey number, settled from their position, their area and the layer they were taken from
G["settled"] = ""
def settle(src, row, vc, sheet, why):
    k = G.index[(G.src == src) & (G.src_row == row)]; assert len(k) == 1, (src, row, len(k)); G.at[k[0], "village_code"] = vc; G.at[k[0], "sheet"] = sheet; G.at[k[0], "settled"] = why
settle("Peramanur_Merged_Fmb", 1, "35_04_076", "43B", "survey number not in the layer: second plot of the source layer 43A_modified, 12.01 ha against 12.39 ha on sheet 43B, lies on 43B")
settle("3VILLAGES", 509, "35_04_054", "226A", "survey number not in the layer: 3.43 ha against 3.17 ha on sheet 226A, lies on 226A")
settle("3VILLAGES", 434, "35_04_056", "332", "village not in the layer: taken from the plots around it; survey 332 has no sheet on the portal")
settle("3VILLAGES", 729, "35_04_052", "50", "village not in the layer: taken from the plots around it; the plot is the whole of survey 50 (sheets 50A + 50B)")
notes = [n for n in notes if not str(n[3]).startswith("village not settled")]
G["village"] = G.village_code.map(NAMES).fillna(""); print(G.groupby(["src", "village_code", "village"]).size().to_string()); print("plots", len(G), "| village unknown", int((G.village_code == "").sum()))
w = {}
for i, g in zip(G.index, G.geometry.normalize()): w.setdefault(g.wkb, []).append(i)
dup = [v for v in w.values() if len(v) > 1]; drop = [i for v in dup for i in v[1:]]
for v in dup: notes.append((G.src[v[0]], list(G.src_row[v]), G.sheet[v[0]], "the same shape %d times - one kept" % len(v)))
G = G.drop(index=drop).reset_index(drop=True); G["area"] = G.area.round(2)
G.to_file("out/merged_raw.gpkg", driver="GPKG"); N = pd.DataFrame(notes, columns=["layer", "row", "survey", "note"]); N.to_csv("out/merge_notes.csv", index=False); print(N.note.str.replace(r"\(.*", "", regex=True).str[:50].value_counts().to_string()); print(N[N.note.str.contains("village not|same shape|invalid|taken from")].to_string())
s = G[(G.src != "Manual_Georeference")]; print("sheets without a PDF:", sorted({(a, b) for a, b in zip(s.village_code, s.sheet) if not os.path.exists(os.path.join(R, "FMB_Sketches", a, b + ".pdf"))}))
